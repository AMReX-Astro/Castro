module hybrid_advection_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains


  subroutine ca_linear_to_hybrid_momentum(lo, hi, state, s_lo, s_hi) bind(C, name='ca_linear_to_hybrid_momentum')
    ! Convert linear momentum to hybrid momentum.

    use meth_params_module, only: NVAR, UMR, UMP, UMX, UMZ
    use castro_util_module, only: position ! function
    use prob_params_module, only: center

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer  :: i, j, k
    real(rt) :: loc(3), mom(3)

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center
             mom = state(i,j,k,UMX:UMZ)

             state(i,j,k,UMR:UMP) = linear_to_hybrid(loc, mom)

          enddo
       enddo
    enddo

  end subroutine ca_linear_to_hybrid_momentum



  subroutine ca_hybrid_to_linear_momentum(lo, hi, state, s_lo, s_hi) bind(C, name='ca_hybrid_to_linear_momentum')
    ! Convert hybrid momentum to linear momentum.

    use meth_params_module, only: NVAR, URHO, UMR, UMP, UMX, UMZ
    use castro_util_module, only: position ! function
    use prob_params_module, only: center

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer  :: i, j, k
    real(rt) :: loc(3), mom(3)

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center
             mom = state(i,j,k,UMR:UMP)

             state(i,j,k,UMX:UMZ) = hybrid_to_linear(loc, mom)

          enddo
       enddo
    enddo

  end subroutine ca_hybrid_to_linear_momentum



  subroutine ca_hybrid_hydro_source(lo, hi, state, s_lo, s_hi, ext_src, e_lo, e_hi, mult_factor) bind(C,name='ca_hybrid_hydro_source')
    ! Fill a sources array with the source terms in the hybrid momentum
    ! equations.

    use amrex_constants_module, only: ONE
    use meth_params_module, only: NVAR, URHO, UMR, UML
    use prob_params_module, only: center
    use castro_util_module, only: position ! function
    use network, only: nspec, naux

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: e_lo(3), e_hi(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: ext_src(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),NVAR)
    real(rt), intent(in   ), value :: mult_factor

    integer  :: i, j, k
    real(rt) :: loc(3), R, rhoInv

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             R = sqrt( loc(1)**2 + loc(2)**2 )

             rhoInv = ONE / state(i,j,k,URHO)

             ext_src(i,j,k,UMR) = ext_src(i,j,k,UMR) + mult_factor * (rhoInv / R**3) * state(i,j,k,UML)**2

          enddo
       enddo
    enddo

  end subroutine ca_hybrid_hydro_source



  function linear_to_hybrid(loc, mom_in) result(mom_out)
    ! Convert a linear momentum into the "hybrid" scheme
    ! that has radial and angular components.

    use amrex_constants_module, only: ZERO

    implicit none

    real(rt), intent(in) :: loc(3), mom_in(3)
    real(rt) :: mom_out(3)

    real(rt) :: R

    !$gpu

    R = sqrt( loc(1)**2 + loc(2)**2 )

    ! This conversion is Eqs. 25 and 26 in Byerly et al. 2014.
    ! Note that we expect the linear momentum to be consistent
    ! with which frame we're measuring the fluid quantities in.
    ! So we're effectively always using the first form of those
    ! equalities, not the second. If state_in_rotating_frame = 1,
    ! then we're not including the centrifugal term in the angular
    ! momentum anyway, and if state_in_rotating_frame = 0, then
    ! the linear momenta are already expressed in the inertial frame,
    ! so we don't need to explicitly take rotation into account.

    mom_out(1) = mom_in(1) * (loc(1) / R) + mom_in(2) * (loc(2) / R)
    mom_out(2) = mom_in(2) * loc(1) - mom_in(1) * loc(2)
    mom_out(3) = mom_in(3)

  end function linear_to_hybrid


  function hybrid_to_linear(loc, mom_in) result(mom_out)
    ! Convert a "hybrid" momentum into a linear one.

    use amrex_constants_module, only: ZERO

    implicit none

    real(rt), intent(in) :: loc(3), mom_in(3)
    real(rt) :: mom_out(3)

    real(rt) :: R

    !$gpu

    R = sqrt( loc(1)**2 + loc(2)**2 )

    ! This is the inverse of Byerly et al., Equations 25 and 26.

    mom_out(1) = mom_in(1) * (loc(1) / R)    - mom_in(2) * (loc(2) / R**2)
    mom_out(2) = mom_in(2) * (loc(1) / R**2) + mom_in(1) * (loc(2) / R)
    mom_out(3) = mom_in(3)

  end function hybrid_to_linear



  subroutine add_hybrid_momentum_source(loc, mom, source)
    ! Update hybrid momenta to account for source term to linear momenta.

    implicit none

    real(rt), intent(in   ) :: loc(3), source(3)
    real(rt), intent(inout) :: mom(3)

    real(rt) :: R

    !$gpu

    R = sqrt( loc(1)**2 + loc(2)**2 )

    ! This is analogous to the conversion of linear momentum to hybrid momentum.

    mom(1) = mom(1) + source(1) * (loc(1) / R) + source(2) * (loc(2) / R)
    mom(2) = mom(2) + source(2) * loc(1) - source(1) * loc(2)
    mom(3) = mom(3) + source(3)

  end subroutine add_hybrid_momentum_source


  subroutine set_hybrid_momentum_source(loc, mom, source)

    implicit none

    real(rt), intent(in   ) :: loc(3), source(3)
    real(rt), intent(inout) :: mom(3)

    real(rt) :: R

    !$gpu

    R = sqrt( loc(1)**2 + loc(2)**2 )

    ! This is analogous to the conversion of linear momentum to hybrid momentum.

    mom(1) =  source(1) * (loc(1) / R) + source(2) * (loc(2) / R)
    mom(2) =  source(2) * loc(1) - source(1) * loc(2)
    mom(3) =  source(3)

  end subroutine set_hybrid_momentum_source



  subroutine compute_hybrid_flux(state, flux, idir, idx, cell_centered)

    use meth_params_module, only: NVAR, NGDNV, GDRHO, GDU, GDV, GDW, GDPRES, UMR, UML, UMP
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use prob_params_module, only: center
    use castro_util_module, only: position ! function

    implicit none

    real(rt)          :: state(NGDNV)
    real(rt)          :: flux(NVAR)
    integer           :: idir, idx(3)
    logical, optional :: cell_centered

    real(rt) :: linear_mom(3), hybrid_mom(3)
    real(rt) :: loc(3), R

    real(rt) :: u_adv
    logical  :: cc

    !$gpu

    cc = .false.

    if (present(cell_centered)) then
       if (cell_centered) then
          cc = .true.
       endif
    endif

    if (idir .eq. 1) then
       loc = position(idx(1),idx(2),idx(3),ccx=cc) - center
       u_adv = state(GDU)
    else if (idir .eq. 2) then
       loc = position(idx(1),idx(2),idx(3),ccy=cc) - center
       u_adv = state(GDV)
    else if (idir .eq. 3) then
       loc = position(idx(1),idx(2),idx(3),ccz=cc) - center
       u_adv = state(GDW)
#ifndef AMREX_USE_CUDA
    else
       call castro_error("Error: unknown direction in compute_hybrid_flux.")
#endif
    endif

    R = sqrt(loc(1)**2 + loc(2)**2)

    linear_mom = state(GDRHO) * state(GDU:GDW)

    hybrid_mom = linear_to_hybrid(loc, linear_mom)

    if (idir .eq. 1) then

       flux(UMR) = hybrid_mom(1) * u_adv
       flux(UML) = hybrid_mom(2) * u_adv - loc(2) * state(GDPRES)
       flux(UMP) = hybrid_mom(3) * u_adv

    else if (idir .eq. 2) then

       flux(UMR) = hybrid_mom(1) * u_adv
       flux(UML) = hybrid_mom(2) * u_adv + loc(1) * state(GDPRES)
       flux(UMP) = hybrid_mom(3) * u_adv

    else if (idir .eq. 3) then

       flux(UMR) = hybrid_mom(1) * u_adv
       flux(UML) = hybrid_mom(2) * u_adv
       flux(UMP) = hybrid_mom(3) * u_adv + state(GDPRES)

#ifndef AMREX_USE_CUDA
    else

       call castro_error("Error: unknown direction in compute_hybrid_flux.")
#endif
    endif

  end subroutine compute_hybrid_flux



  subroutine add_hybrid_advection_source(lo, hi, dt, &
                                         update, u_lo, u_hi, &
                                         qx, qx_lo, qx_hi, &
                                         qy, qy_lo, qy_hi, &
                                         qz, qz_lo, qz_hi)

    use meth_params_module, only: NVAR, NGDNV, GDPRES, UMR
    use prob_params_module, only: center, dx_level
    use castro_util_module, only: position ! function
    use amrinfo_module, only: amr_level

    implicit none

    integer  :: lo(3), hi(3)
    integer  :: u_lo(3), u_hi(3)
    integer  :: qx_lo(3), qx_hi(3)
    integer  :: qy_lo(3), qy_hi(3)
    integer  :: qz_lo(3), qz_hi(3)
    real(rt) :: update(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt) :: dt

    integer  :: i, j, k
    real(rt) :: loc(3), R, dx(3)

    !$gpu

    dx = dx_level(:,amr_level)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             R = sqrt( loc(1)**2 + loc(2)**2 )

             update(i,j,k,UMR) = update(i,j,k,UMR) - ( (loc(1) / R) * (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx(1) + &
                  (loc(2) / R) * (qy(i,j+1,k,GDPRES) - qy(i,j,k,GDPRES)) / dx(2) )

          enddo
       enddo
    enddo

  end subroutine add_hybrid_advection_source

end module hybrid_advection_module
