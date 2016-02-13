module advection_util_module

  implicit none

  private

  public enforce_minimum_density, &
         normalize_species_fluxes, divu

contains

  subroutine normalize_species_fluxes(  &
                    flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                    flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                    lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module
    
    implicit none

    integer          :: lo(2),hi(2)
    integer          :: flux1_l1,flux1_l2,flux1_h1,flux1_h2
    integer          :: flux2_l1,flux2_l2,flux2_h1,flux2_h2
    double precision :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
    double precision :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
    
    ! Local variables
    integer          :: i,j,n
    double precision :: sum,fac
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + flux1(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = flux1(i,j,URHO) / sum
          else
             fac = ONE
          end if
          do n = UFS, UFS+nspec-1
             flux1(i,j,n) = flux1(i,j,n) * fac
          end do
       end do
    end do
    do j = lo(2),hi(2)+1
       do i = lo(1),hi(1)
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + flux2(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = flux2(i,j,URHO) / sum
          else
             fac = ONE
          end if
          do n = UFS, UFS+nspec-1
             flux2(i,j,n) = flux2(i,j,n) * fac
          end do
       end do
    end do
    
  end subroutine normalize_species_fluxes

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine enforce_minimum_density(uin,uin_lo,uin_hi, &
                                     uout,uout_lo,uout_hi, &
                                     lo,hi,mass_added,eint_added, &
                                     eden_added,frac_change,verbose) &
                                     bind(C, name="enforce_minimum_density")
    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEINT, UEDEN, UTEMP, &
                                   UFS, small_dens, smalL_temp, npassive, upass_map
    use bl_constants_module
    use eos_type_module, only : eos_t
    use eos_module, only : eos_input_rt, eos

    implicit none

    integer          :: lo(2), hi(2), verbose
    integer          ::  uin_lo(2),  uin_hi(2)
    integer          :: uout_lo(2), uout_hi(2)
    double precision ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2),NVAR)
    double precision :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),NVAR)
    double precision :: mass_added, eint_added, eden_added, frac_change
    
    ! Local variables
    integer                       :: i,ii,j,jj,n,ipassive
    double precision              :: max_dens
    integer                       :: i_set, j_set
    double precision              :: initial_mass, final_mass
    double precision              :: initial_eint, final_eint
    double precision              :: initial_eden, final_eden

    type (eos_t) :: eos_state
    
    initial_mass = ZERO
    final_mass = ZERO
    initial_eint = ZERO
    final_eint = ZERO
    initial_eden = ZERO
    final_eden = ZERO

    max_dens = ZERO
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          initial_mass = initial_mass + uout(i,j,URHO)
          initial_eint = initial_eint + uout(i,j,UEINT)
          initial_eden = initial_eden + uout(i,j,UEDEN)

          if (uout(i,j,URHO) .eq. ZERO) then
             
             print *,'   '
             print *,'>>> Error: Castro_2d::enforce_minimum_density ',i,j
             print *,'>>> ... density exactly zero in grid ',lo(1),hi(1),lo(2),hi(2)
             print *,'    '
             call bl_error("Error:: Castro_2d.f90 :: enforce_minimum_density")
             
          else if (uout(i,j,URHO) < small_dens) then

             ! Store the maximum (negative) fractional change in the density
             
             if ( uout(i,j,URHO) < ZERO .and. &
                  (uout(i,j,URHO) - uin(i,j,URHO)) / uin(i,j,URHO) < frac_change) then

                frac_change = (uout(i,j,URHO) - uin(i,j,URHO)) / uin(i,j,URHO)

             endif

             max_dens = uout(i,j,URHO)
             do jj = -1,1
                do ii = -1,1
                   if (i+ii.ge.lo(1) .and. j+jj.ge.lo(2) .and. &
                       i+ii.le.hi(1) .and. j+jj.le.hi(2)) then
                        if (uout(i+ii,j+jj,URHO) > max_dens) then
                           i_set = i+ii
                           j_set = j+jj
                           max_dens = uout(i_set,j_set,URHO)
                        endif
                   endif
                end do
             end do

             ! If no neighboring zones are above small_dens, our only recourse 
             ! is to set the density equal to small_dens, and the temperature 
             ! equal to small_temp. We set the velocities to zero, 
             ! though any choice here would be arbitrary.

             if (max_dens < small_dens) then

                i_set = i
                j_set = j
                
                do ipassive = 1, npassive
                   n = upass_map(ipassive)
                   uout(i,j,n) = uout(i,j,n) * (small_dens / uout(i,j,URHO))
                end do

                eos_state % rho = small_dens
                eos_state % T   = small_temp
                eos_state % xn  = uout(i,j,UFS:UFS+nspec-1) / small_dens

                call eos(eos_input_rt, eos_state)

                uout(i,j,URHO ) = eos_state % rho
                uout(i,j,UTEMP) = eos_state % T

                uout(i,j,UMX  ) = ZERO
                uout(i,j,UMY  ) = ZERO

                uout(i,j,UEINT) = eos_state % rho * eos_state % e
                uout(i,j,UEDEN) = uout(i,j,UEINT)

             endif
             
             if (verbose .gt. 0) then
                if (uout(i,j,URHO) < ZERO) then
                   print *,'   '
                   print *,'>>> Warning: Castro_2d::enforce_minimum_density ',i,j
                   print *,'>>> ... resetting negative density '
                   print *,'>>> ... from ',uout(i,j,URHO),' to ',uout(i_set,j_set,URHO)
                   print *,'    '
                else
                   print *,'   '
                   print *,'>>> Warning: Castro_2d::enforce_minimum_density ',i,j
                   print *,'>>> ... resetting small density '
                   print *,'>>> ... from ',uout(i,j,URHO),' to ',uout(i_set,j_set,URHO)
                   print *,'    '
                end if
             end if

             uout(i,j,URHO ) = uout(i_set,j_set,URHO )
             uout(i,j,UTEMP) = uout(i_set,j_set,UTEMP)
             uout(i,j,UEINT) = uout(i_set,j_set,UEINT)
             uout(i,j,UEDEN) = uout(i_set,j_set,UEDEN)
             uout(i,j,UMX  ) = uout(i_set,j_set,UMX  )
             uout(i,j,UMY  ) = uout(i_set,j_set,UMY  )

             do ipassive = 1, npassive
                n = upass_map(ipassive)
                uout(i,j,n) = uout(i_set,j_set,n)
             end do
             
          end if

          final_mass = final_mass + uout(i,j,URHO )
          final_eint = final_eint + uout(i,j,UEINT)
          final_eden = final_eden + uout(i,j,UEDEN)                
          
       enddo
    enddo

    if (max_dens /= ZERO) then
    
       mass_added = mass_added + (final_mass - initial_mass)
       eint_added = eint_added + (final_eint - initial_eint)
       eden_added = eden_added + (final_eden - initial_eden)

    endif
    
  end subroutine enforce_minimum_density

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine normalize_new_species(u,u_l1,u_l2,u_h1,u_h2,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module
    
    implicit none

    integer          :: lo(2), hi(2)
    integer          :: u_l1,u_l2,u_h1,u_h2
    double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
    
    ! Local variables
    integer          :: i,j,n
    double precision :: fac,sum
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + u(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = u(i,j,URHO) / sum
          else
             fac = ONE
          end if
          do n = UFS, UFS+nspec-1
             u(i,j,n) = u(i,j,n) * fac
          end do
       end do
    end do
    
  end subroutine normalize_new_species

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2,dx, &
                  div,div_l1,div_l2,div_h1,div_h2)

    use prob_params_module, only : coord_type
    use meth_params_module, only : QU, QV
    use bl_constants_module
    
    implicit none
    
    integer          :: lo(2),hi(2)
    integer          :: q_l1,q_l2,q_h1,q_h2
    integer          :: div_l1,div_l2,div_h1,div_h2
    double precision :: q(q_l1:q_h1,q_l2:q_h2,*)
    double precision :: div(div_l1:div_h1,div_l2:div_h2)
    double precision :: dx(2)
    
    integer          :: i, j
    double precision :: rl, rr, rc, ul, ur
    double precision :: vb, vt
    double precision :: ux,vy
    
    if (coord_type .eq. 0) then
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1
             ux = HALF*(q(i,j,QU)-q(i-1,j,QU)+q(i,j-1,QU)-q(i-1,j-1,QU))/dx(1)
             vy = HALF*(q(i,j,QV)-q(i,j-1,QV)+q(i-1,j,QV)-q(i-1,j-1,QV))/dx(2)
             div(i,j) = ux + vy
          enddo
       enddo
    else
       do i=lo(1),hi(1)+1
          
          if (i.eq.0) then
             
             div(i,lo(2):hi(2)+1) = ZERO
             
          else 

             rl = (dble(i)-HALF) * dx(1)
             rr = (dble(i)+HALF) * dx(1)
             rc = (dble(i)     ) * dx(1)
             
             do j=lo(2),hi(2)+1
                ! These are transverse averages in the y-direction
                ul = HALF * (q(i-1,j,QU)+q(i-1,j-1,QU))
                ur = HALF * (q(i  ,j,QU)+q(i  ,j-1,QU))
                
                ! Take 1/r d/dr(r*u)
                div(i,j) = (rr*ur - rl*ul) / dx(1) / rc
                
                ! These are transverse averages in the x-direction
                vb = HALF * (q(i,j-1,QV)+q(i-1,j-1,QV))
                vt = HALF * (q(i,j  ,QV)+q(i-1,j  ,QV))
                
                div(i,j) = div(i,j) + (vt - vb) / dx(2)
             enddo
             
          end if
       enddo
    end if

  end subroutine divu

end module advection_util_module
