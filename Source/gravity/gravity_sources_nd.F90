module gravity_sources_module

  implicit none

  public

contains

  subroutine ca_gsrc(lo, hi, &
                     domlo, domhi, &
                     uold, uold_lo, uold_hi, &
                     phi, phi_lo, phi_hi, &
                     grav, grav_lo, grav_hi, &
                     source, src_lo, src_hi, &
                     dx, dt, time) bind(C, name="ca_gsrc")

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, grav_source_type, NSRC
    use castro_util_module, only: position ! function
    use prob_params_module, only: center
#ifdef HYBRID_MOMENTUM
    use meth_params_module, only: UMR, UMP
    use hybrid_advection_module, only: set_hybrid_momentum_source
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uold_lo(3), uold_hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: src_lo(3), src_hi(3)

    real(rt), intent(in   ) :: uold(uold_lo(1):uold_hi(1),uold_lo(2):uold_hi(2),uold_lo(3):uold_hi(3),NVAR)
    real(rt), intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3),3)
    real(rt), intent(inout) :: source(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NSRC)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, time

    real(rt) :: rho, rhoInv
    real(rt) :: Sr(3), SrE
    real(rt) :: old_ke, new_ke
    real(rt) :: loc(3)

    integer  :: i, j, k

    ! Temporary array for holding the update to the state.
    
    real(rt) :: src(NSRC)

    ! Temporary array for seeing what the new state would be if the update were applied here.

    real(rt) :: snew(NVAR)

    !$gpu

    ! Initialize the update and temporary state to zero. We only need to do this once outside
    ! the loop, since the array access pattern is consistent across loop iterations.

    Sr(:) = ZERO
    src(:) = ZERO
    snew(:) = ZERO

    ! Gravitational source options for how to add the work to (rho E):
    ! grav_source_type =
    ! 1: Original version ("does work")
    ! 2: Modification of type 1 that updates the momentum before constructing the energy corrector
    ! 3: Puts all gravitational work into KE, not (rho e)
    ! 4: Conservative energy formulation

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rho    = uold(i,j,k,URHO)
             rhoInv = ONE / rho

             loc = position(i,j,k) - center

             src = ZERO
             snew = uold(i,j,k,:)

             old_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoInv

             Sr = rho * grav(i,j,k,:)

             src(UMX:UMZ) = Sr

             snew(UMX:UMZ) = snew(UMX:UMZ) + dt * src(UMX:UMZ)

#ifdef HYBRID_MOMENTUM
             call set_hybrid_momentum_source(loc, src(UMR:UMP), Sr)

             snew(UMR:UMP) = snew(UMR:UMP) + dt * src(UMR:UMP)
#endif

             if (grav_source_type == 1 .or. grav_source_type == 2) then

                ! Src = rho u dot g, evaluated with all quantities at t^n

                SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)

             else if (grav_source_type .eq. 3) then

                new_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoInv
                SrE = new_ke - old_ke

             else if (grav_source_type .eq. 4) then

                ! The conservative energy formulation does not strictly require
                ! any energy source-term here, because it depends only on the
                ! fluid motions from the hydrodynamical fluxes which we will only
                ! have when we get to the 'corrector' step. Nevertheless we add a
                ! predictor energy source term in the way that the other methods
                ! do, for consistency. We will fully subtract this predictor value
                ! during the corrector step, so that the final result is correct.
                ! Here we use the same approach as grav_source_type == 2.

                SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)

#ifndef AMREX_USE_CUDA
             else

                call castro_error("Error:: gravity_sources_nd.F90 :: invalid grav_source_type")
#endif

             end if

             src(UEDEN) = SrE

             snew(UEDEN) = snew(UEDEN) + dt * SrE

             ! Add to the outgoing source array.

             source(i,j,k,:) = source(i,j,k,:) + src

          enddo
       enddo
    enddo

  end subroutine ca_gsrc

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine ca_corrgsrc(lo, hi, &
                         domlo, domhi, &
                         uold, uo_lo, uo_hi, &
                         unew, un_lo, un_hi, &
                         pold, po_lo, po_hi, &
                         pnew, pn_lo, pn_hi, &
                         gold, go_lo, go_hi, &
                         gnew, gn_lo, gn_hi, &
                         vol, vol_lo, vol_hi, &
                         flux1, f1_lo, f1_hi, &
                         flux2, f2_lo, f2_hi, &
                         flux3, f3_lo, f3_hi, &
                         source, sr_lo, sr_hi, &
                         dx, dt, time) bind(C, name="ca_corrgsrc")

    use amrex_fort_module, only: rt => amrex_real
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use amrex_constants_module, only: ZERO, HALF, ONE, TWO
    use amrex_mempool_module, only: bl_allocate, bl_deallocate
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, NSRC, &
                                 grav_source_type, gravity_type_int, PoissonGrav, &
                                 MonopoleGrav, get_g_from_phi, ambient_safety_factor
    use prob_params_module, only: dg, center, physbc_lo, physbc_hi, Symmetry
    use castro_util_module, only: position ! function
    use ambient_module, only: ambient_state
#ifdef HYBRID_MOMENTUM
    use meth_params_module, only: UMR, UMP
    use hybrid_advection_module, only: set_hybrid_momentum_source
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uo_lo(3), uo_hi(3)
    integer,  intent(in   ) :: un_lo(3), un_hi(3)
    integer,  intent(in   ) :: po_lo(3), po_hi(3)
    integer,  intent(in   ) :: pn_lo(3), pn_hi(3)
    integer,  intent(in   ) :: go_lo(3), go_hi(3)
    integer,  intent(in   ) :: gn_lo(3), gn_hi(3)
    integer,  intent(in   ) :: vol_lo(3), vol_hi(3)
    integer,  intent(in   ) :: f1_lo(3), f1_hi(3)
    integer,  intent(in   ) :: f2_lo(3), f2_hi(3)
    integer,  intent(in   ) :: f3_lo(3), f3_hi(3)

    integer,  intent(in   ) :: sr_lo(3), sr_hi(3)

    ! Old and new time state data

    real(rt), intent(in   ) :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),NVAR)
    real(rt), intent(in   ) :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),NVAR)

    ! Old and new time gravitational potential

    real(rt), intent(in   ) :: pold(po_lo(1):po_hi(1),po_lo(2):po_hi(2),po_lo(3):po_hi(3))
    real(rt), intent(in   ) :: pnew(pn_lo(1):pn_hi(1),pn_lo(2):pn_hi(2),pn_lo(3):pn_hi(3))

    ! Old and new time gravitational acceleration

    real(rt), intent(in   ) :: gold(go_lo(1):go_hi(1),go_lo(2):go_hi(2),go_lo(3):go_hi(3),3)
    real(rt), intent(in   ) :: gnew(gn_lo(1):gn_hi(1),gn_lo(2):gn_hi(2),gn_lo(3):gn_hi(3),3)

    ! Cell volume

    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))

    ! Hydrodynamical mass fluxes

    real(rt), intent(in   ) :: flux1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3))
    real(rt), intent(in   ) :: flux2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3))
    real(rt), intent(in   ) :: flux3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3))

    ! The source term to send back

    real(rt), intent(inout) :: source(sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3),NSRC)

    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, time

    integer  :: i, j, k

    real(rt) :: Sr_old(3), Sr_new(3), Srcorr(3)
    real(rt) :: vold(3), vnew(3)
    real(rt) :: SrE_old, SrE_new, SrEcorr, dSrE_cons, dSrE_non_cons, rho_c, delta_rho
    real(rt) :: rhoo, rhooinv, rhon, rhoninv

    real(rt) :: old_ke, new_ke
    real(rt) :: loc(3)

    real(rt) :: hdtInv

    real(rt) :: phi, phixl, phixr, phiyl, phiyr, phizl, phizr
    real(rt) :: g(3), gxl, gxr, gyl, gyr, gzl, gzr

    real(rt) :: src(NSRC)

    ! Temporary array for seeing what the new state would be if the update were applied here.

    real(rt) :: snew(NVAR)

    !$gpu

    Sr_old(:) = ZERO
    Sr_new(:) = ZERO
    Srcorr(:) = ZERO
    src(:) = ZERO
    snew(:) = ZERO

    hdtInv = HALF / dt

    ! Gravitational source options for how to add the work to (rho E):
    ! grav_source_type =
    ! 1: Original version ("does work")
    ! 2: Modification of type 1 that updates the U before constructing SrEcorr
    ! 3: Puts all gravitational work into KE, not (rho e)
    ! 4: Conservative gravity approach (discussed in first white dwarf merger paper).

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             rhoo    = uold(i,j,k,URHO)
             rhooinv = ONE / uold(i,j,k,URHO)

             rhon    = unew(i,j,k,URHO)
             rhoninv = ONE / unew(i,j,k,URHO)

             snew = unew(i,j,k,:)

             old_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoninv

             ! Define old source terms

             vold = uold(i,j,k,UMX:UMZ) * rhooinv

             Sr_old = rhoo * gold(i,j,k,:)
             SrE_old = dot_product(vold, Sr_old)

             ! Define new source terms

             vnew = snew(UMX:UMZ) * rhoninv

             Sr_new = rhon * gnew(i,j,k,:)
             SrE_new = dot_product(vnew, Sr_new)

             ! Define corrections to source terms

             Srcorr = HALF * (Sr_new - Sr_old)

             ! Correct momenta

             src(UMX:UMZ) = Srcorr

             snew(UMX:UMZ) = snew(UMX:UMZ) + dt * src(UMX:UMZ)

#ifdef HYBRID_MOMENTUM
             call set_hybrid_momentum_source(loc, src(UMR:UMP), Srcorr)

             snew(UMR:UMP) = snew(UMR:UMP) + dt * src(UMR:UMP)
#endif

             ! Correct energy

             if (grav_source_type .eq. 1) then

                ! If grav_source_type == 1, then we calculated SrEcorr before updating the velocities.

                SrEcorr = HALF * (SrE_new - SrE_old)

             else if (grav_source_type .eq. 2) then

                ! For this source type, we first update the momenta
                ! before we calculate the energy source term.

                vnew = snew(UMX:UMZ) * rhoninv
                SrE_new = dot_product(vnew, Sr_new)

                SrEcorr = HALF * (SrE_new - SrE_old)

             else if (grav_source_type .eq. 3) then

                ! Instead of calculating the energy source term explicitly,
                ! we simply update the kinetic energy.

                new_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoninv
                SrEcorr = new_ke - old_ke

             else if (grav_source_type .eq. 4) then

                ! First, subtract the predictor step we applied earlier.

                SrEcorr = - SrE_old

                ! For an explanation of this approach, see wdmerger paper I.
                ! The main idea is that we are evaluating the change of the
                ! potential energy at zone edges and applying that in an equal
                ! and opposite sense to the gas energy. The physics is described
                ! in Section 2.4; we are using a version of the formula similar to
                ! Equation 94 in Springel (2010) based on the gradient rather than
                ! the potential because the gradient-version works for all forms
                ! of gravity we use, some of which do not explicitly calculate phi.

                ! Construct the time-averaged edge-centered gravity.

                g(:) = HALF * (gnew(i,j,k,:) + gold(i,j,k,:))

                gxl = HALF * (g(1) + HALF * (gnew(i-1*dg(1),j,k,1) + gold(i-1*dg(1),j,k,1)))
                gxr = HALF * (g(1) + HALF * (gnew(i+1*dg(1),j,k,1) + gold(i+1*dg(1),j,k,1)))

                gyl = HALF * (g(2) + HALF * (gnew(i,j-1*dg(2),k,2) + gold(i,j-1*dg(2),k,2)))
                gyr = HALF * (g(2) + HALF * (gnew(i,j+1*dg(2),k,2) + gold(i,j+1*dg(2),k,2)))

                gzl = HALF * (g(3) + HALF * (gnew(i,j,k-1*dg(3),3) + gold(i,j,k-1*dg(3),3)))
                gzr = HALF * (g(3) + HALF * (gnew(i,j,k+1*dg(3),3) + gold(i,j,k+1*dg(3),3)))

                dSrE_cons = hdtInv * (flux1(i        ,j,k) * gxl * dx(1) + &
                                      flux1(i+1*dg(1),j,k) * gxr * dx(1) + &
                                      flux2(i,j        ,k) * gyl * dx(2) + &
                                      flux2(i,j+1*dg(2),k) * gyr * dx(2) + &
                                      flux3(i,j,k        ) * gzl * dx(3) + &
                                      flux3(i,j,k+1*dg(3)) * gzr * dx(3)) / vol(i,j,k)

                ! The conservative scheme above is suboptimal for low-density, ambient material.
                ! Because it ties energy updates to mass flow into/out of the zone at the edges,
                ! it has the possibility of generating an energy update that is quite inconsistent
                ! with the zone's current mass and energy. This is most likely to occur at the
                ! interface between high and low density material, such as the edge of a star. So
                ! we also construct here the standard non-conservative update and apply it for
                ! the low density material. This sacrifices a small amount of energy conservation
                ! in return for avoiding significant changes in the internal energy of the
                ! ambient material, which could result in, for example, very high sound speeds.

                vnew = snew(UMX:UMZ) * rhoninv
                SrE_new = dot_product(vnew, Sr_new)

                dSrE_non_cons = HALF * (SrE_new - SrE_old)

                ! Smooth between the two schemes, using ambient_safety_factor as the cutoff point
                ! between the normal material and the ambient material.

                rho_c = (ONE + HALF * (ambient_safety_factor - ONE)) * ambient_state(URHO)
                delta_rho = (HALF * (ambient_safety_factor - ONE)) * ambient_state(URHO)

                dSrE_non_cons = dSrE_non_cons * HALF * (ONE + tanh((rhon - rho_c) / delta_rho))
                dSrE_cons = dSrE_cons * HALF * (ONE - tanh((rhon - rho_c) / delta_rho))

                SrEcorr = SrEcorr + dSrE_non_cons + dSrE_cons

#ifndef AMREX_USE_CUDA
             else

                call castro_error("Error:: gravity_sources_nd.F90 :: invalid grav_source_type")
#endif
             end if

             src(UEDEN) = SrEcorr

             snew(UEDEN) = snew(UEDEN) + dt * SrEcorr

             ! Add to the outgoing source array.

             source(i,j,k,:) = source(i,j,k,:) + src

          enddo
       enddo
    enddo

  end subroutine ca_corrgsrc

end module gravity_sources_module
