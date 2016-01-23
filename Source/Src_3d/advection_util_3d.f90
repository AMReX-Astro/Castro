module advection_util_module

  implicit none

  private

  public enforce_minimum_density, normalize_new_species, &
         normalize_species_fluxes, divu

contains

  subroutine normalize_species_fluxes(flux1,flux1_lo,flux1_hi, &
                                      flux2,flux2_lo,flux2_hi, &
                                      flux3,flux3_lo,flux3_hi, &
                                      lo, hi)
    
    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: flux1_lo(3), flux1_hi(3)
    integer          :: flux2_lo(3), flux2_hi(3)
    integer          :: flux3_lo(3), flux3_hi(3)
    double precision :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    double precision :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    double precision :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    
    ! Local variables
    integer          :: i, j, k, n
    double precision :: sum, fac
    
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

  subroutine enforce_minimum_density(uin,uin_lo,uin_hi, &
                                     uout,uout_lo,uout_hi, &
                                     lo,hi,mass_added,eint_added,eden_added,verbose) bind(C)
    
    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, &
                                   small_dens, small_temp, npassive, upass_map
    use bl_constants_module
    use eos_module

    implicit none

    integer          :: lo(3), hi(3), verbose
    integer          ::  uin_lo(3),  uin_hi(3)
    integer          :: uout_lo(3), uout_hi(3)
    double precision ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
    double precision :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
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

  subroutine normalize_new_species(u,u_lo,u_hi,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    
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


! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine divu(lo,hi,q,q_lo,q_hi,dx,div,div_lo,div_hi)
    
    use meth_params_module, only : QU, QV, QW, QVAR
    use bl_constants_module
    
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: q_lo(3), q_hi(3)
    integer          :: div_lo(3), div_hi(3)
    double precision :: dx(3)
    double precision :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    double precision :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    integer          :: i, j, k
    double precision :: ux, vy, wz, dxinv, dyinv, dzinv

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1
             
             ux = FOURTH*( &
                    + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                    + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                    + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                    + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) ) * dxinv

             vy = FOURTH*( &
                    + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                    + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                    + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                    + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) ) * dyinv

             wz = FOURTH*( &
                    + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                    + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                    + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                    + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) ) * dzinv

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo
    
  end subroutine divu

end module advection_util_module

