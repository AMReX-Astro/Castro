module hybrid_advection_module

  implicit none

contains

  ! Takes the initial linear momentum data in a state and converts it
  ! to the hybrid momenta.

  subroutine init_hybrid_momentum(lo, hi, state, s_lo, s_hi) bind(C,name='init_hybrid_momentum')

    use meth_params_module, only: NVAR, UMR, UMP, UMX, UMZ
    use castro_util_module, only: position
    use prob_params_module, only: center

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer          :: i, j, k
    double precision :: loc(3)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             state(i,j,k,UMR:UMP) = linear_to_hybrid_momentum(loc, state(i,j,k,UMX:UMZ))

          enddo
       enddo
    enddo

  end subroutine init_hybrid_momentum



  ! Fill a sources array with the source terms in the hybrid momentum equations.

  subroutine ca_hybrid_hydro_source(lo, hi, state, s_lo, s_hi, ext_src, e_lo, e_hi) bind(C,name='ca_hybrid_hydro_source')

    use meth_params_module, only: NVAR, URHO, UEINT, UTEMP, UFS, UFX, UMR, UML
    use prob_params_module, only: center
    use castro_util_module, only: position
    use network, only: nspec, naux
    use eos_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    integer          :: e_lo(3), e_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: ext_src(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),NVAR)

    integer          :: i, j, k
    double precision :: loc(3), R, rhoInv, rho, pres

    type (eos_t)     :: eos_state

    eos_state % reset = .true.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             R = sqrt( loc(1)**2 + loc(2)**2 )

             rho = state(i,j,k,URHO)
             rhoInv = ONE / rho

             ! Calculate the pressure

             eos_state % rho = rho
             eos_state % T   = state(i,j,k,UTEMP)
             eos_state % e   = state(i,j,k,UEINT) * rhoInv
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_re, eos_state)

             pres = eos_state % p

             ! Now we can calculate the source term

             ext_src(i,j,k,UMR) = ext_src(i,j,k,UMR) &
                                + (rhoInv / R**3) * (state(i,j,k,UML)**2 + pres * rho * R**2)

          enddo
       enddo
    enddo

  end subroutine ca_hybrid_hydro_source



  ! Convert a linear momentum into a "hybrid" momentum that has
  ! an angular momentum component.

  function linear_to_hybrid_momentum(loc, mom_in) result(mom_out)

    implicit none

    double precision :: loc(3), mom_in(3), mom_out(3)

    double precision :: R, mom(3)

    R = sqrt( loc(1)**2 + loc(2)**2 )

    mom = mom_in
    
    mom_out(1) = mom(1) * (loc(1) / R) + mom(2) * (loc(2) / R)
    mom_out(2) = mom(2) * loc(1)       - mom(1) * loc(2)
    mom_out(3) = mom(3)

  end function linear_to_hybrid_momentum



  ! Convert a "hybrid" momentum into a linear momentum.

  function hybrid_to_linear_momentum(loc, mom_in) result(mom_out)

    implicit none

    double precision :: loc(3), mom_in(3), mom_out(3)

    double precision :: R, mom(3)

    mom = mom_in

    R = sqrt( loc(1)**2 + loc(2)**2 )
    
    mom_out(1) = mom(1) * (loc(1) / R)    - mom(2) * (loc(2) / R**2)
    mom_out(2) = mom(2) * (loc(1) / R**2) + mom(1) * (loc(2) / R)
    mom_out(3) = mom(3)
    
  end function hybrid_to_linear_momentum



  ! Update momentum to account for source term

  subroutine add_hybrid_momentum_source(loc, mom, source)

    implicit none

    double precision :: loc(3), mom(3), source(3)

    double precision :: R

    R = sqrt( loc(1)**2 + loc(2)**2 )

    mom(1) = mom(1) - source(1) * (loc(1) / R) - source(2) * (loc(2) / R)
    mom(2) = mom(2) + source(1) * loc(2) - source(2) * loc(1)
    mom(3) = mom(3) + source(3)

  end subroutine add_hybrid_momentum_source



  subroutine compute_hybrid_flux(state, flux, idir, idx)

    use meth_params_module, only: NVAR, NGDNV, GDRHO, GDU, GDV, GDW, GDPRES, UMR, UML, UMP
    use bl_error_module, only: bl_error
    use prob_params_module, only: center
    use castro_util_module, only: position

    implicit none

    double precision :: state(NGDNV)
    double precision :: flux(NVAR)
    integer          :: idir, idx(3)

    double precision :: linear_mom(3), hybrid_mom(3)
    double precision :: loc(3), R

    if (idir .eq. 1) then
       loc = position(idx(1),idx(2),idx(3),ccx=.false.) - center
    else if (idir .eq. 2) then
       loc = position(idx(1),idx(2),idx(3),ccy=.false.) - center
    else if (idir .eq. 3) then
       loc = position(idx(1),idx(2),idx(3),ccz=.false.) - center
    else
       call bl_error("Error: unknown direction in compute_hybrid_flux.")
    endif

    R = sqrt(loc(1)**2 + loc(2)**2)

    linear_mom = state(GDRHO) * state(GDU:GDW)

    hybrid_mom = linear_to_hybrid_momentum(loc, linear_mom)

    if (idir .eq. 1) then

       flux(UMR) = hybrid_mom(1) * state(GDU) + (loc(1) / R) * state(GDPRES)
       flux(UML) = hybrid_mom(2) * state(GDU) + loc(2) * state(GDPRES)
       flux(UMP) = hybrid_mom(3) * state(GDU)

    else if (idir .eq. 2) then

       flux(UMR) = hybrid_mom(1) * state(GDV) + (loc(2) / R) * state(GDPRES)
       flux(UML) = hybrid_mom(2) * state(GDV) - loc(1) * state(GDPRES)
       flux(UMP) = hybrid_mom(3) * state(GDV)

    else if (idir .eq. 3) then

       flux(UMR) = hybrid_mom(1) * state(GDW)
       flux(UML) = hybrid_mom(2) * state(GDW)
       flux(UMP) = hybrid_mom(3) * state(GDW)

    else

       call bl_error("Error: unknown direction in compute_hybrid_flux.")

    endif
       
  end subroutine compute_hybrid_flux



  ! Update state to account for hybrid advection.

  subroutine hybrid_update(lo, hi, state, state_lo, state_hi) bind(C,name='hybrid_update')

    use meth_params_module, only: UMR, UML, UMP, UMX, UMZ, NVAR, hybrid_hydro
    use castro_util_module, only: position

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: state_lo(3), state_hi(3)
    double precision :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

    integer          :: i, j, k
    double precision :: loc(3)

    ! If we're doing the hybrid advection scheme, update the momenta accordingly.

    if (hybrid_hydro .eq. 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k)

                state(i,j,k,UMX:UMZ) = hybrid_to_linear_momentum(loc, state(i,j,k,UMR:UMP))

             enddo
          enddo
       enddo

    endif

  end subroutine hybrid_update

end module hybrid_advection_module

