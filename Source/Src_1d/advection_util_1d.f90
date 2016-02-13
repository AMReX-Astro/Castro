module advection_util_module

  implicit none

  private

  public enforce_minimum_density, normalize_species_fluxes

contains

  subroutine normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none
    
    integer          :: lo(1),hi(1)
    integer          :: flux_l1,flux_h1
    double precision :: flux(flux_l1:flux_h1,NVAR)
    
    ! Local variables
    integer          :: i,n
    double precision :: sum,fac
    
    do i = lo(1),hi(1)+1
       sum = ZERO
       do n = UFS, UFS+nspec-1
          sum = sum + flux(i,n)
       end do
       if (sum .ne. ZERO) then
          fac = flux(i,URHO) / sum
       else
          fac = ONE
       end if
       do n = UFS, UFS+nspec-1
          flux(i,n) = flux(i,n) * fac
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
    use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, &
                                   UFS, UFX, small_dens, small_temp, npassive, upass_map
    use bl_constants_module
    use eos_type_module, only : eos_t
    use eos_module, only : eos_input_rt, eos

    implicit none
    
    integer          :: lo(1), hi(1), verbose
    integer          ::  uin_lo(1),  uin_hi(1)
    integer          :: uout_lo(1), uout_hi(1)
    double precision ::  uin( uin_lo(1): uin_hi(1),NVAR)
    double precision :: uout(uout_lo(1):uout_hi(1),NVAR)
    double precision :: mass_added, eint_added, eden_added, frac_change
    
    ! Local variables
    integer          :: i,ii,n,ipassive
    double precision :: max_dens
    integer          :: i_set
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
    
    do i = lo(1),hi(1)
       
       initial_mass = initial_mass + uout(i,URHO)
       initial_eint = initial_eint + uout(i,UEINT)
       initial_eden = initial_eden + uout(i,UEDEN)
       
       if (uout(i,URHO) .eq. ZERO) then
          
          print *,'   '
          print *,'>>> Error: Castro_1d::enforce_minimum_density ',i
          print *,'>>> ... density exactly zero in grid ',lo(1),hi(1)
          print *,'    '
          call bl_error("Error:: Castro_1d.f90 :: enforce_minimum_density")
          
       else if (uout(i,URHO) < small_dens) then

          ! Store the maximum (negative) fractional change in the density
          
          if ( uout(i,URHO) < ZERO .and. &
               (uout(i,URHO) - uin(i,URHO)) / uin(i,URHO) < frac_change) then

             frac_change = (uout(i,URHO) - uin(i,URHO)) / uin(i,URHO)

          endif

          max_dens = uout(i,URHO)
          do ii = -1,1
             if (i+ii.ge.lo(1) .and. i+ii.le.hi(1)) then
                if (uout(i+ii,URHO) .gt. max_dens) then
                   i_set = i+ii
                   max_dens = uout(i_set,URHO)
                endif
             endif
          end do

          ! If no neighboring zones are above small_dens, our only recourse 
          ! is to set the density equal to small_dens, and the temperature 
          ! equal to small_temp. We set the velocities to zero, 
          ! though any choice here would be arbitrary.

          if (max_dens < small_dens) then

             i_set = i
             
             do ipassive = 1, npassive
                n = upass_map(ipassive)
                uout(i,n) = uout(i,n) * (small_dens / uout(i,URHO))
             end do

             eos_state % rho = small_dens
             eos_state % T   = small_temp
             eos_state % xn  = uout(i,UFS:UFS+nspec-1) / small_dens

             call eos(eos_input_rt, eos_state)

             uout(i,URHO ) = eos_state % rho
             uout(i,UTEMP) = eos_state % T

             uout(i,UMX  ) = ZERO

             uout(i,UEINT) = eos_state % rho * eos_state % e
             uout(i,UEDEN) = uout(i,UEINT)

          endif
          
          if (verbose .gt. 0) then
             if (uout(i,URHO) < ZERO) then
                print *,'   '
                print *,'>>> Warning: Castro_1d::enforce_minimum_density ',i
                print *,'>>> ... resetting negative density '
                print *,'>>> ... from ',uout(i,URHO),' to ',uout(i_set,URHO)
                print *,'    '
             else
                print *,'   '
                print *,'>>> Warning: Castro_1d::enforce_minimum_density ',i
                print *,'>>> ... resetting small density '
                print *,'>>> ... from ',uout(i,URHO),' to ',uout(i_set,URHO)
                print *,'    '
             end if
          end if

          uout(i,URHO ) = uout(i_set,URHO )
          uout(i,UTEMP) = uout(i_set,UTEMP)
          uout(i,UEINT) = uout(i_set,UEINT)
          uout(i,UEDEN) = uout(i_set,UEDEN)
          uout(i,UMX  ) = uout(i_set,UMX  )

          do ipassive = 1, npassive
             n = upass_map(ipassive)
             uout(i,n) = uout(i_set,n)
          end do
          
       end if

       final_mass = final_mass + uout(i,URHO)
       final_eint = final_eint + uout(i,UEINT)
       final_eden = final_eden + uout(i,UEDEN)

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

  subroutine normalize_new_species(u,u_l1,u_h1,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module    

    implicit none

    integer          :: lo(1), hi(1)
    integer          :: u_l1,u_h1
    double precision :: u(u_l1:u_h1,NVAR)
    
    ! Local variables
    integer          :: i,n
    double precision :: fac,sum
    
    do i = lo(1),hi(1)
       sum = ZERO
       do n = UFS, UFS+nspec-1
          sum = sum + u(i,n)
       end do
       if (sum .ne. ZERO) then
          fac = u(i,URHO) / sum
       else
          fac = ONE
       end if
       do n = UFS, UFS+nspec-1
          u(i,n) = u(i,n) * fac
       end do
    end do
    
  end subroutine normalize_new_species
  
end module advection_util_module
