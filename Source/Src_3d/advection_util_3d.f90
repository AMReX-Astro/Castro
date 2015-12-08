module advection_util_module

  implicit none

  private

  public enforce_minimum_density, normalize_new_species, normalize_species_fluxes

contains

  subroutine normalize_species_fluxes(flux1,flux1_l1,flux1_l2,flux1_l3, &
                                      flux1_h1,flux1_h2,flux1_h3, &
                                      flux2,flux2_l1,flux2_l2,flux2_l3, &
                                      flux2_h1,flux2_h2,flux2_h3, &
                                      flux3,flux3_l1,flux3_l2,flux3_l3, &
                                      flux3_h1,flux3_h2,flux3_h3, &
                                      lo,hi)
    
    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(3),hi(3)
    integer          :: flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer          :: flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer          :: flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    double precision :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision :: flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
    
    ! Local variables
    integer          :: i,j,k,n
    double precision :: sum,fac
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux1(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux1(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux1(i,j,k,n) = flux1(i,j,k,n) * fac
             end do
          end do
       end do
    end do

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux2(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux2(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux2(i,j,k,n) = flux2(i,j,k,n) * fac
             end do
          end do
       end do
    end do

    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux3(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux3(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux3(i,j,k,n) = flux3(i,j,k,n) * fac
             end do
          end do
       end do
    end do

  end subroutine normalize_species_fluxes

! ::
! :: ----------------------------------------------------------
! ::

  subroutine enforce_minimum_density(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                                     uout,uout_l1,uout_l2,uout_l3, &
                                     uout_h1,uout_h2,uout_h3, &
                                     lo,hi,mass_added,eint_added,eden_added,verbose)
    
    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, &
                                   small_dens, small_temp, npassive, upass_map
    use bl_constants_module
    use eos_module

    implicit none

    integer          :: lo(3), hi(3), verbose
    integer          ::  uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
    integer          :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
    double precision ::  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
    double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision :: mass_added, eint_added, eden_added
    
    ! Local variables
    integer          :: i,ii,j,jj,k,kk,n,ipassive
    integer          :: i_set, j_set, k_set
    double precision :: max_dens
    
    double precision :: initial_mass, final_mass
    double precision :: initial_eint, final_eint
    double precision :: initial_eden, final_eden

    type (eos_t) :: eos_state
    
    initial_mass = ZERO
      final_mass = ZERO

    initial_eint = ZERO
      final_eint = ZERO

    initial_eden = ZERO
      final_eden = ZERO

    max_dens = ZERO

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             initial_mass = initial_mass + uout(i,j,k,URHO )
             initial_eint = initial_eint + uout(i,j,k,UEINT)
             initial_eden = initial_eden + uout(i,j,k,UEDEN)
             
             if (uout(i,j,k,URHO) .eq. ZERO) then
                
                print *,'DENSITY EXACTLY ZERO AT CELL ',i,j,k
                print *,'  in grid ',lo(1),lo(2),lo(3),hi(1),hi(2),hi(3)
                call bl_error("Error:: Castro_3d.f90 :: enforce_minimum_density")
                
             else if (uout(i,j,k,URHO) < small_dens) then
                
                max_dens = uout(i,j,k,URHO)
                i_set = i
                j_set = j
                k_set = k
                do kk = -1,1
                   do jj = -1,1
                      do ii = -1,1
                         if (i+ii.ge.lo(1) .and. j+jj.ge.lo(2) .and. k+kk.ge.lo(3) .and. &
                             i+ii.le.hi(1) .and. j+jj.le.hi(2) .and. k+kk.le.hi(3)) then
                              if (uout(i+ii,j+jj,k+kk,URHO) .gt. max_dens) then
                                  i_set = i+ii
                                  j_set = j+jj
                                  k_set = k+kk
                                  max_dens = uout(i_set,j_set,k_set,URHO)
                              endif
                         endif
                      end do
                   end do
                end do

                ! If no neighboring zones are above small_dens, our only recourse 
                ! is to set the density equal to small_dens, and the temperature 
                ! equal to small_temp. We set the velocities to zero, 
                ! though any choice here would be arbitrary.

                if (max_dens < small_dens) then

                   i_set = i
                   j_set = j
                   k_set = k
                   
                   do ipassive = 1, npassive
                      n = upass_map(ipassive)
                      uout(i,j,k,n) = uout(i,j,k,n) * (small_dens / uout(i,j,k,URHO))
                   end do

                   eos_state % rho = small_dens
                   eos_state % T   = small_temp
                   eos_state % xn  = uout(i,j,k,UFS:UFS+nspec-1) / small_dens

                   call eos(eos_input_rt, eos_state)

                   uout(i,j,k,URHO ) = eos_state % rho
                   uout(i,j,k,UTEMP) = eos_state % T

                   uout(i,j,k,UMX  ) = ZERO
                   uout(i,j,k,UMY  ) = ZERO
                   uout(i,j,k,UMZ  ) = ZERO

                   uout(i,j,k,UEINT) = eos_state % rho * eos_state % e
                   uout(i,j,k,UEDEN) = uout(i,j,k,UEINT)

                endif
                
                if (verbose .gt. 0) then
                   if (uout(i,j,k,URHO) < ZERO) then
                      print *,'   '
                      print *,'>>> RESETTING NEG.  DENSITY AT ',i,j,k
                      print *,'>>> FROM ',uout(i,j,k,URHO),' TO ',uout(i_set,j_set,k_set,URHO)
                      print *,'   '
                   else
                      print *,'   '
                      print *,'>>> RESETTING SMALL DENSITY AT ',i,j,k
                      print *,'>>> FROM ',uout(i,j,k,URHO),' TO ',uout(i_set,j_set,k_set,URHO)
                      print *,'   '
                   end if
                end if
                
                uout(i,j,k,URHO ) = uout(i_set,j_set,k_set,URHO )
                uout(i,j,k,UTEMP) = uout(i_set,j_set,k_set,UTEMP)
                uout(i,j,k,UEINT) = uout(i_set,j_set,k_set,UEINT)
                uout(i,j,k,UEDEN) = uout(i_set,j_set,k_set,UEDEN)
                uout(i,j,k,UMX  ) = uout(i_set,j_set,k_set,UMX  )
                uout(i,j,k,UMY  ) = uout(i_set,j_set,k_set,UMY  )
                uout(i,j,k,UMZ  ) = uout(i_set,j_set,k_set,UMZ  )
   
                do ipassive = 1, npassive
                   n = upass_map(ipassive)
                   uout(i,j,k,n) = uout(i_set,j_set,k_set,n)
                end do
                
             end if

             final_mass = final_mass + uout(i,j,k,URHO )
             final_eint = final_eint + uout(i,j,k,UEINT)
             final_eden = final_eden + uout(i,j,k,UEDEN)                
             
          enddo
       enddo
    enddo

    if ( max_dens /= ZERO ) then
       mass_added = mass_added + final_mass - initial_mass
       eint_added = eint_added + final_eint - initial_eint
       eden_added = eden_added + final_eden - initial_eden
    endif

  end subroutine enforce_minimum_density

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine normalize_new_species(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
    double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
    
    ! Local variables
    integer          :: i,j,k,n
    double precision :: fac,sum
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + u(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = u(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                u(i,j,k,n) = u(i,j,k,n) * fac
             end do
          end do
       end do
    end do
    
  end subroutine normalize_new_species

end module advection_util_module

