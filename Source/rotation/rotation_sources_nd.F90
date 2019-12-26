module rotation_sources_module

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_rsrc(lo,hi,domlo,domhi,phi,phi_lo,phi_hi,rot,rot_lo,rot_hi, &
                     uold,uold_lo,uold_hi, &
                     source,src_lo,src_hi,vol,vol_lo,vol_hi, &
                     dx,dt,time) bind(C, name="ca_rsrc")

    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, rot_source_type, NSRC
    use prob_params_module, only: center
    use amrex_constants_module
    use castro_util_module, only: position ! function
#ifdef HYBRID_MOMENTUM
    use meth_params_module, only: UMR, UMP, state_in_rotating_frame
    use hybrid_advection_module, only: set_hybrid_momentum_source
#endif
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: domlo(3), domhi(3)
    integer         , intent(in   ) :: phi_lo(3), phi_hi(3)
    integer         , intent(in   ) :: rot_lo(3), rot_hi(3)
    integer         , intent(in   ) :: uold_lo(3), uold_hi(3)
    integer         , intent(in   ) :: src_lo(3), src_hi(3)
    integer         , intent(in   ) :: vol_lo(3), vol_hi(3)

    real(rt)        , intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt)        , intent(in   ) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3),3)
    real(rt)        , intent(in   ) :: uold(uold_lo(1):uold_hi(1),uold_lo(2):uold_hi(2),uold_lo(3):uold_hi(3),NVAR)
    real(rt)        , intent(inout) :: source(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NSRC)
    real(rt)        , intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt)        , intent(in   ) :: dx(3)
    real(rt), value , intent(in   ) :: dt, time

    integer          :: i, j ,k
    real(rt)         :: Sr(3),SrE
    real(rt)         :: rho, rhoInv

    real(rt)         :: old_ke, new_ke
    real(rt)         :: loc(3)

    real(rt)         :: src(NSRC)

    ! Temporary array for seeing what the new state would be if the update were applied here.

    real(rt)         :: snew(NVAR)

    !$gpu

    Sr(:) = ZERO
    src(:) = ZERO
    snew(:) = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             rho = uold(i,j,k,URHO)
             rhoInv = ONE / rho

             src(:) = ZERO
             snew = uold(i,j,k,:)

             old_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoInv

             Sr = rho * rot(i,j,k,:)

             src(UMX:UMZ) = Sr

             snew(UMX:UMZ) = snew(UMX:UMZ) + dt * src(UMX:UMZ)

#ifdef HYBRID_MOMENTUM
             if (state_in_rotating_frame == 1) then
                call set_hybrid_momentum_source(loc, src(UMR:UMP), Sr)

                snew(UMR:UMP) = snew(UMR:UMP) + dt * src(UMR:UMP)
             endif
#endif

             ! The conservative energy formulation does not strictly require
             ! any energy source-term here, because it depends only on the
             ! fluid motions from the hydrodynamical fluxes which we will only
             ! have when we get to the 'corrector' step. Nevertheless we add a
             ! predictor energy source term in the way that the other methods
             ! do, for consistency. We will fully subtract this predictor value
             ! during the corrector step, so that the final result is correct.

             if (rot_source_type == 1 .or. rot_source_type == 2) then

                SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)

             else if (rot_source_type .eq. 3) then

                new_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoInv
                SrE = new_ke - old_ke

             else if (rot_source_type .eq. 4) then

                ! The conservative energy formulation does not strictly require
                ! any energy source-term here, because it depends only on the
                ! fluid motions from the hydrodynamical fluxes which we will only
                ! have when we get to the 'corrector' step. Nevertheless we add a
                ! predictor energy source term in the way that the other methods
                ! do, for consistency. We will fully subtract this predictor value
                ! during the corrector step, so that the final result is correct.
                ! Here we use the same approach as rot_source_type == 2.

                SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)

             else
#ifndef AMREX_USE_GPU
                call castro_error("Error:: rotation_sources_nd.F90 :: invalid rot_source_type")
#endif
             end if

             src(UEDEN) = src(UEDEN) + SrE

             snew(UEDEN) = snew(UEDEN) + dt * src(UEDEN)

             ! Add to the outgoing source array.

             source(i,j,k,:) = source(i,j,k,:) + src

          enddo
       enddo
    enddo

  end subroutine ca_rsrc



  subroutine ca_corrrsrc(lo,hi,domlo,domhi, &
                         phi_old,po_lo,po_hi, &
                         phi_new,pn_lo,pn_hi, &
                         rold,ro_lo,ro_hi, &
                         rnew,rn_lo,rn_hi, &
                         uold,uo_lo,uo_hi, &
                         unew,un_lo,un_hi, &
                         source,sr_lo,sr_hi, &
                         flux1,f1_lo,f1_hi, &
                         flux2,f2_lo,f2_hi, &
                         flux3,f3_lo,f3_hi, &
                         dx,dt,time, &
                         vol,vol_lo,vol_hi) bind(C, name="ca_corrrsrc")
    ! Corrector step for the rotation source terms. This is applied
    ! after the hydrodynamics update to fix the time-level n
    ! prediction and add the time-level n+1 data.  This subroutine
    ! exists outside of the Fortran module above because it needs to
    ! be called directly from C++.

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, rot_source_type, NSRC, &
                                  implicit_rotation_update, rotation_include_coriolis, state_in_rotating_frame
    use prob_params_module, only: center, dg
    use amrex_constants_module
    use math_module, only: cross_product ! function
    use rotation_module, only: rotational_acceleration ! function
    use rotation_frequency_module, only: get_omega ! function
    use rotation_frequency_module, only: get_domegadt ! function
    use castro_util_module, only: position ! function
#ifdef HYBRID_MOMENTUM
    use meth_params_module, only : UMR, UMP
    use hybrid_advection_module, only: set_hybrid_momentum_source
#endif
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: domlo(3), domhi(3)

    integer          :: po_lo(3),po_hi(3)
    integer          :: pn_lo(3),pn_hi(3)
    integer          :: ro_lo(3),ro_hi(3)
    integer          :: rn_lo(3),rn_hi(3)
    integer          :: uo_lo(3),uo_hi(3)
    integer          :: un_lo(3),un_hi(3)
    integer          :: sr_lo(3),sr_hi(3)
    integer          :: f1_lo(3),f1_hi(3)
    integer          :: f2_lo(3),f2_hi(3)
    integer          :: f3_lo(3),f3_hi(3)
    integer          :: vol_lo(3),vol_hi(3)

    ! Time centered rotational potential

    real(rt)         :: phi_old(po_lo(1):po_hi(1),po_lo(2):po_hi(2),po_lo(3):po_hi(3))
    real(rt)         :: phi_new(pn_lo(1):pn_hi(1),pn_lo(2):pn_hi(2),pn_lo(3):pn_hi(3))

    ! Old and new time rotational acceleration

    real(rt)         :: rold(ro_lo(1):ro_hi(1),ro_lo(2):ro_hi(2),ro_lo(3):ro_hi(3),3)
    real(rt)         :: rnew(rn_lo(1):rn_hi(1),rn_lo(2):rn_hi(2),rn_lo(3):rn_hi(3),3)

    ! Old and new time state data

    real(rt)         :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),NVAR)
    real(rt)         :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),NVAR)

    ! The source term to send back

    real(rt)         :: source(sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3),NSRC)

    ! Hydrodynamical mass fluxes

    real(rt)         :: flux1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3))
    real(rt)         :: flux2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3))
    real(rt)         :: flux3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3))

    real(rt)         :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))

    real(rt)         :: dx(3)
    real(rt), value  :: dt, time

    integer          :: i,j,k
    real(rt)         :: loc(3)
    real(rt)         :: vnew(3),vold(3),omega_old(3),omega_new(3),domegadt_old(3),domegadt_new(3)
    real(rt)         :: Sr_old(3), Sr_new(3), Srcorr(3), SrEcorr, SrE_old, SrE_new
    real(rt)         :: rhoo, rhon, rhooinv, rhoninv

    real(rt)         :: old_ke, new_ke
    real(rt)         :: dt_omega_matrix(3,3), dt_omega(3), new_mom(3)

    real(rt)         :: src(NSRC)

    real(rt)         :: phi, phixl, phixr, phiyl, phiyr, phizl, phizr

    ! Temporary array for seeing what the new state would be if the update were applied here.

    real(rt)         :: snew(NVAR)

    ! Note that the time passed to this function
    ! is the new time at time-level n+1.

    !$gpu

    Sr_old(:) = ZERO
    Sr_new(:) = ZERO
    Srcorr(:) = ZERO
    src(:) = ZERO
    snew(:) = ZERO

    omega_old = get_omega(time-dt)
    omega_new = get_omega(time   )

    domegadt_old = get_domegadt(time-dt)
    domegadt_new = get_domegadt(time   )

    if (implicit_rotation_update == 1) then

       ! Don't do anything here if we've got the Coriolis force disabled.

       if (rotation_include_coriolis == 1) then

          ! If the state variables are in the inertial frame, then we are doing
          ! an implicit solve using (dt / 2) multiplied by the standard Coriolis term.
          ! If not, then the rotation source term to the linear momenta (Equations 16
          ! and 17 in Byerly et al., 2014) still retains a Coriolis-like form, with
          ! the only difference being that the magnitude is half as large. Consequently
          ! we can still do an implicit solve in that case.

          if (state_in_rotating_frame == 1) then

             dt_omega = dt * omega_new

          else

             dt_omega = HALF * dt * omega_new

          endif

       else

          dt_omega = ZERO

       endif

       dt_omega_matrix(1,1) = ONE + dt_omega(1)**2
       dt_omega_matrix(1,2) = dt_omega(1) * dt_omega(2) + dt_omega(3)
       dt_omega_matrix(1,3) = dt_omega(1) * dt_omega(3) - dt_omega(2)

       dt_omega_matrix(2,1) = dt_omega(2) * dt_omega(1) - dt_omega(3)
       dt_omega_matrix(2,2) = ONE + dt_omega(2)**2
       dt_omega_matrix(2,3) = dt_omega(2) * dt_omega(3) + dt_omega(1)

       dt_omega_matrix(3,1) = dt_omega(3) * dt_omega(1) + dt_omega(2)
       dt_omega_matrix(3,2) = dt_omega(3) * dt_omega(2) - dt_omega(1)
       dt_omega_matrix(3,3) = ONE + dt_omega(3)**2

       dt_omega_matrix = dt_omega_matrix / (ONE + dt_omega(1)**2 + dt_omega(2)**2 + dt_omega(3)**2)

    endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             rhoo = uold(i,j,k,URHO)
             rhooinv = ONE / uold(i,j,k,URHO)

             rhon = unew(i,j,k,URHO)
             rhoninv = ONE / unew(i,j,k,URHO)

             src = ZERO
             snew = unew(i,j,k,:)

             old_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoninv

             ! Define old source terms

             vold = uold(i,j,k,UMX:UMZ) * rhooinv

             Sr_old = rhoo * rold(i,j,k,:)
             SrE_old = dot_product(vold, Sr_old)

             ! Define new source terms

             vnew = unew(i,j,k,UMX:UMZ) * rhoninv

             Sr_new = rhon * rnew(i,j,k,:)
             SrE_new = dot_product(vnew, Sr_new)

             ! Define correction terms

             Srcorr = HALF * (Sr_new - Sr_old)

             if (implicit_rotation_update == 1) then

                ! Coupled/implicit momentum update (wdmerger paper I; Section 2.4)
                ! http://adsabs.harvard.edu/abs/2016ApJ...819...94K

                ! Do the full corrector step with the old contribution (subtract 1/2 times the old term) and do
                ! the non-Coriolis parts of the new contribution (add 1/2 of the new term).

                new_mom = unew(i,j,k,UMX:UMZ) - HALF * Sr_old * dt + &
                          HALF * rhon * rotational_acceleration(loc, vnew, time, coriolis = .false.) * dt

                ! The following is the general solution to the 3D coupled system,
                ! assuming that the rotation vector has components along all three
                ! axes, obtained using Cramer's rule (the coefficient matrix is
                ! defined above). In practice the user will probably only be using
                ! one axis for rotation; if it's the z-axis, then this reduces to
                ! Equations 25 and 26 in the wdmerger paper. Note that this will
                ! have the correct form regardless of whether the state variables are
                ! measured in the rotating frame or not; we handled that in the construction
                ! of the dt_omega_matrix. It also has the correct form if we have disabled
                ! the Coriolis force entirely; at that point it reduces to the identity matrix.

                new_mom = matmul(dt_omega_matrix, new_mom)

                ! Obtain the effective source term; remember that we're ultimately going
                ! to multiply the source term by dt to get the update to the state.

                Srcorr = (new_mom - unew(i,j,k,UMX:UMZ)) / dt

             endif

             ! Correct momenta

             src(UMX:UMZ) = Srcorr

             snew(UMX:UMZ) = snew(UMX:UMZ) + dt * src(UMX:UMZ)

#ifdef HYBRID_MOMENTUM
             ! The source terms vanish if the state variables are measured in the
             ! inertial frame; see wdmerger paper III.

             if (state_in_rotating_frame == 1) then
                call set_hybrid_momentum_source(loc, src(UMR:UMP), Srcorr)

                snew(UMR:UMP) = snew(UMR:UMP) + dt * src(UMR:UMP)
             endif
#endif

             ! Correct energy

             if (rot_source_type == 1) then

                ! If rot_source_type == 1, then we calculated SrEcorr before updating the velocities.

                SrEcorr = HALF * (SrE_new - SrE_old)

             else if (rot_source_type == 2) then

                ! For this source type, we first update the momenta
                ! before we calculate the energy source term.

                vnew = snew(UMX:UMZ) * rhoninv
                Sr_new = rhon * rotational_acceleration(loc, vnew, time)
                SrE_new = dot_product(vnew, Sr_new)

                SrEcorr = HALF * (SrE_new - SrE_old)

             else if (rot_source_type == 3) then

                ! Instead of calculating the energy source term explicitly,
                ! we simply update the kinetic energy.

                new_ke = HALF * sum(snew(UMX:UMZ)**2) * rhoninv
                SrEcorr = new_ke - old_ke

             else if (rot_source_type == 4) then

                ! Conservative energy update

                ! First, subtract the predictor step we applied earlier.

                SrEcorr = - SrE_old

                ! The change in the gas energy is equal in magnitude to, and opposite in sign to,
                ! the change in the rotational potential energy, rho * phi.
                ! This must be true for the total energy, rho * E_gas + rho * phi, to be conserved.
                ! Consider as an example the zone interface i+1/2 in between zones i and i + 1.
                ! There is an amount of mass drho_{i+1/2} leaving the zone. From this zone's perspective
                ! it starts with a potential phi_i and leaves the zone with potential phi_{i+1/2} =
                ! (1/2) * (phi_{i-1}+phi_{i}). Therefore the new rotational energy is equal to the mass
                ! change multiplied by the difference between these two potentials.
                ! This is a generalization of the cell-centered approach implemented in
                ! the other source options, which effectively are equal to
                ! SrEcorr = - drho(i,j,k) * phi(i,j,k),
                ! where drho(i,j,k) = HALF * (unew(i,j,k,URHO) - uold(i,j,k,URHO)).

                ! Note that in the hydrodynamics step, the fluxes used here were already
                ! multiplied by dA and dt, so dividing by the cell volume is enough to
                ! get the density change (flux * dt * dA / dV). We then divide by dt
                ! so that we get the source term and not the actual update, which will
                ! be applied later by multiplying by dt.

                phi = HALF * (phi_new(i,j,k) + phi_old(i,j,k))
                phixl = HALF * (phi_new(i-1*dg(1),j,k) + phi_old(i-1*dg(1),j,k))
                phixr = HALF * (phi_new(i+1*dg(1),j,k) + phi_old(i+1*dg(1),j,k))
                phiyl = HALF * (phi_new(i,j-1*dg(2),k) + phi_old(i,j-1*dg(2),k))
                phiyr = HALF * (phi_new(i,j+1*dg(2),k) + phi_old(i,j+1*dg(2),k))
                phizl = HALF * (phi_new(i,j,k-1*dg(3)) + phi_old(i,j,k-1*dg(3)))
                phizr = HALF * (phi_new(i,j,k+1*dg(3)) + phi_old(i,j,k+1*dg(3)))

                SrEcorr = SrEcorr - (HALF / dt) * ( flux1(i        ,j,k) * (phi - phixl) - &
                                                    flux1(i+1*dg(1),j,k) * (phi - phixr) + &
                                                    flux2(i,j        ,k) * (phi - phiyl) - &
                                                    flux2(i,j+1*dg(2),k) * (phi - phiyr) + &
                                                    flux3(i,j,k        ) * (phi - phizl) - &
                                                    flux3(i,j,k+1*dg(3)) * (phi - phizr) ) / vol(i,j,k)

                ! Correct for the time rate of change of the potential, which acts
                ! purely as a source term. This is only necessary for this source type;
                ! it is captured automatically for the others since the time rate of change
                ! of omega also appears in the velocity source term.

                Sr_old = - rhoo * cross_product(domegadt_old, loc)
                Sr_new = - rhon * cross_product(domegadt_new, loc)

                vnew = snew(UMX:UMZ) * rhoninv

                SrEcorr = SrEcorr + HALF * (dot_product(vold, Sr_old) + dot_product(vnew, Sr_new))

             else
#ifndef AMREX_USE_GPU
                call castro_error("Error:: rotation_sources_nd.F90 :: invalid rot_source_type")
#endif
             end if

             src(UEDEN) = SrEcorr

             ! Add to the outgoing source array.

             source(i,j,k,:) = source(i,j,k,:) + src

          enddo
       enddo
    enddo

  end subroutine ca_corrrsrc

end module rotation_sources_module
