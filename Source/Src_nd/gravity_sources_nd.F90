module gravity_sources_module

  implicit none

  public

contains

  subroutine ca_gsrc(lo,hi,domlo,domhi, &
                     uold,uold_lo,uold_hi, &
#ifdef SELF_GRAVITY
                     grav,grav_lo,grav_hi, &
#endif
                     source,src_lo,src_hi, &
                     dx,dt,time) bind(C, name="ca_gsrc")

    use bl_fort_module, only: rt => c_real
    use bl_constants_module, only: ZERO, HALF, ONE
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, grav_source_type
    use math_module, only: cross_product
    use castro_util_module, only: position
    use prob_params_module, only: center
#ifdef HYBRID_MOMENTUM
    use meth_params_module, only: UMR, UMP
    use hybrid_advection_module, only: add_hybrid_momentum_source
#endif
#ifndef SELF_GRAVITY
    use meth_params_module, only: const_grav
    use prob_params_module, only: dim
#endif

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: uold_lo(3), uold_hi(3)
#ifdef SELF_GRAVITY
    integer, intent(in)     :: grav_lo(3), grav_hi(3)
#endif
    integer, intent(in)     :: src_lo(3), src_hi(3)

    real(rt), intent(in)    :: uold(uold_lo(1):uold_hi(1),uold_lo(2):uold_hi(2),uold_lo(3):uold_hi(3),NVAR)
#ifdef SELF_GRAVITY
    real(rt), intent(in)    :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3),3)
#endif
    real(rt), intent(inout) :: source(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    real(rt), intent(in)    :: dx(3), dt, time

    real(rt) :: rho, rhoInv
    real(rt) :: Sr(3), SrE
    real(rt) :: old_ke, new_ke
    real(rt) :: loc(3)

    integer  :: i, j, k

    ! Temporary array for holding the update to the state.
    
    real(rt) :: src(NVAR)

    ! Temporary array for seeing what the new state would be if the update were applied here.

    real(rt) :: snew(NVAR)

    ! Initialize the update and temporary state to zero. We only need to do this once outside
    ! the loop, since the array access pattern is consistent across loop iterations.

    Sr(:) = ZERO
    src(:) = ZERO
    snew(:) = ZERO

    ! For constant gravity, we can just initialize the gravitational acceleration here.

#ifndef SELF_GRAVITY
    Sr(dim) = const_grav
#endif

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

#ifdef SELF_GRAVITY
             Sr = rho * grav(i,j,k,:)
#else
             Sr(dim) = rho * const_grav
#endif

             src(UMX:UMZ) = Sr

             snew(UMX:UMZ) = snew(UMX:UMZ) + dt * src(UMX:UMZ)

#ifdef HYBRID_MOMENTUM
             call add_hybrid_momentum_source(loc, src(UMR:UMP), Sr)

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

             else

                call bl_error("Error:: gravity_sources_nd.F90 :: invalid grav_source_type")

             end if

             src(UEDEN) = SrE

             snew(UEDEN) = snew(UEDEN) + dt * SrE

             ! Add to the outgoing source array.

             source(i,j,k,:) = src

          enddo
       enddo
    enddo

  end subroutine ca_gsrc

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine ca_corrgsrc(lo,hi,domlo,domhi, &
                         uold,uo_lo,uo_hi, &
                         unew,un_lo,un_hi, &
#ifdef SELF_GRAVITY
                         gold,go_lo,go_hi, &
                         gnew,gn_lo,gn_hi, &
                         gpold1,gpo1_lo,gpo1_hi, &
                         gpold2,gpo2_lo,gpo2_hi, &
                         gpold3,gpo3_lo,gpo3_hi, &
                         gpnew1,gpn1_lo,gpn1_hi, &
                         gpnew2,gpn2_lo,gpn2_hi, &
                         gpnew3,gpn3_lo,gpn3_hi, &
#endif
                         vol,vol_lo,vol_hi, &
                         flux1,f1_lo,f1_hi, &
                         flux2,f2_lo,f2_hi, &
                         flux3,f3_lo,f3_hi, &
                         source,sr_lo,sr_hi, &
                         dx,dt,time) bind(C, name="ca_corrgsrc")

    use bl_fort_module, only: rt => c_real
    use bl_constants_module, only: ZERO, FOURTH, HALF, ONE
    use mempool_module, only: bl_allocate, bl_deallocate
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, &
                                  grav_source_type, gravity_type, get_g_from_phi
    use prob_params_module, only: dg, center
    use fundamental_constants_module, only: Gconst
    use castro_util_module, only: position
#ifdef HYBRID_MOMENTUM
    use meth_params_module, only: UMR, UMP
    use hybrid_advection_module, only: add_hybrid_momentum_source
#endif
#ifndef SELF_GRAVITY
    use meth_params_module, only: const_grav
    use prob_params_module, only: dim
#endif

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: uo_lo(3), uo_hi(3)
    integer, intent(in)     :: un_lo(3), un_hi(3)
#ifdef SELF_GRAVITY
    integer, intent(in)     :: go_lo(3), go_hi(3)
    integer, intent(in)     :: gn_lo(3), gn_hi(3)
    integer, intent(in)     :: gpo1_lo(3), gpo1_hi(3)
    integer, intent(in)     :: gpo2_lo(3), gpo2_hi(3)
    integer, intent(in)     :: gpo3_lo(3), gpo3_hi(3)
    integer, intent(in)     :: gpn1_lo(3), gpn1_hi(3)
    integer, intent(in)     :: gpn2_lo(3), gpn2_hi(3)
    integer, intent(in)     :: gpn3_lo(3), gpn3_hi(3)
#endif
    integer, intent(in)     :: vol_lo(3), vol_hi(3)
    integer, intent(in)     :: f1_lo(3), f1_hi(3)
    integer, intent(in)     :: f2_lo(3), f2_hi(3)
    integer, intent(in)     :: f3_lo(3), f3_hi(3)

    integer, intent(in)     :: sr_lo(3), sr_hi(3)

    ! Old and new time state data

    real(rt), intent(in)    :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),NVAR)
    real(rt), intent(in)    :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),NVAR)

#ifdef SELF_GRAVITY
    ! Old and new time gravitational acceleration

    real(rt), intent(in)    :: gold(go_lo(1):go_hi(1),go_lo(2):go_hi(2),go_lo(3):go_hi(3),3)
    real(rt), intent(in)    :: gnew(gn_lo(1):gn_hi(1),gn_lo(2):gn_hi(2),gn_lo(3):gn_hi(3),3)

    ! Edge centered gravitational acceleration

    real(rt), intent(in)    :: gpold1(gpo1_lo(1):gpo1_hi(1),gpo1_lo(2):gpo1_hi(2),gpo1_lo(3):gpo1_hi(3))
    real(rt), intent(in)    :: gpold2(gpo2_lo(1):gpo2_hi(1),gpo2_lo(2):gpo2_hi(2),gpo2_lo(3):gpo2_hi(3))
    real(rt), intent(in)    :: gpold3(gpo3_lo(1):gpo3_hi(1),gpo3_lo(2):gpo3_hi(2),gpo3_lo(3):gpo3_hi(3))

    real(rt), intent(in)    :: gpnew1(gpn1_lo(1):gpn1_hi(1),gpn1_lo(2):gpn1_hi(2),gpn1_lo(3):gpn1_hi(3))
    real(rt), intent(in)    :: gpnew2(gpn2_lo(1):gpn2_hi(1),gpn2_lo(2):gpn2_hi(2),gpn2_lo(3):gpn2_hi(3))
    real(rt), intent(in)    :: gpnew3(gpn3_lo(1):gpn3_hi(1),gpn3_lo(2):gpn3_hi(2),gpn3_lo(3):gpn3_hi(3))
#endif

    ! Cell volume

    real(rt), intent(in)    :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))

    ! Hydrodynamics fluxes

    real(rt), intent(in)    :: flux1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NVAR)
    real(rt), intent(in)    :: flux2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3),NVAR)
    real(rt), intent(in)    :: flux3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3),NVAR)

    ! The source term to send back

    real(rt), intent(inout) :: source(sr_lo(1):sr_hi(1),sr_lo(2):sr_hi(2),sr_lo(3):sr_hi(3),NVAR)

    real(rt), intent(in)    :: dx(3), dt, time

    integer  :: i, j, k

    real(rt) :: Sr_old(3), Sr_new(3), Srcorr(3)
    real(rt) :: vold(3), vnew(3)
    real(rt) :: SrE_old, SrE_new, SrEcorr
    real(rt) :: rhoo, rhooinv, rhon, rhoninv

    real(rt) :: old_ke, new_ke
    real(rt) :: loc(3)

    real(rt) :: hdtInv

    real(rt) :: src(NVAR)

    ! Temporary array for seeing what the new state would be if the update were applied here.

    real(rt) :: snew(NVAR)

    real(rt), pointer :: gravx(:,:,:)
    real(rt), pointer :: gravy(:,:,:)
    real(rt), pointer :: gravz(:,:,:)

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
    ! 4: Conservative gravity approach (discussed in Katz et al. (2016), ApJ, 819, 94).

    if (grav_source_type .eq. 4) then

       ! Allocate space for time-centered, edge-centered gravity.

       call bl_allocate(gravx, lo(1), hi(1)+1*dg(1), lo(2), hi(2), lo(3), hi(3))
       gravx(:,:,:) = ZERO

       call bl_allocate(gravy, lo(1), hi(1), lo(2), hi(2)+1*dg(2), lo(3), hi(3))
       gravy(:,:,:) = ZERO

       call bl_allocate(gravz, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3)+1*dg(3))
       gravz(:,:,:) = ZERO

       ! Construct the gravity array. Note that for Poisson gravity,
       ! grad(phi) is equal to -g. For non-Poisson gravity, the edge-centered
       ! gravity array is not filled with valid data, so we must rely on the
       ! cell-centered gravity and construct a second-order approximation.

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1*dg(1)

#ifdef SELF_GRAVITY
                if (gravity_type == "PoissonGrav") then

                   gravx(i,j,k) = -HALF * (gpold1(i,j,k) + gpnew1(i,j,k))

                else

                   gravx(i,j,k) = FOURTH * (gold(i,j,k,1) + gnew(i,j,k,1) + gold(i-1,j,k,1) + gnew(i-1,j,k,1))

                endif
#else
                ! For constant gravity, the only contribution is from the direction the gravity is
                ! pointing in (which, by convention, is the last spatial dimension).

                if (dim .eq. 1) then

                   gravx(i,j,k) = const_grav

                endif
#endif

             enddo
          enddo
       enddo

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1*dg(2)
             do i = lo(1), hi(1)

#ifdef SELF_GRAVITY
                if (gravity_type == "PoissonGrav") then

                   gravy(i,j,k) = -HALF * (gpold2(i,j,k) + gpnew2(i,j,k))

                else

                   gravy(i,j,k) = FOURTH * (gold(i,j,k,2) + gnew(i,j,k,2) + gold(i,j-1,k,2) + gnew(i,j-1,k,2))

                endif
#else
                if (dim .eq. 2) then

                   gravy(i,j,k) = const_grav

                endif
#endif

             enddo
          enddo
       enddo

       do k = lo(3), hi(3)+1*dg(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

#ifdef SELF_GRAVITY
                if (gravity_type == "PoissonGrav") then

                   gravz(i,j,k) = -HALF * (gpold3(i,j,k) + gpnew3(i,j,k))

                else

                   gravz(i,j,k) = FOURTH * (gold(i,j,k,3) + gnew(i,j,k,3) + gold(i,j,k-1,3) + gnew(i,j,k-1,3))

                endif
#else
                if (dim .eq. 3) then

                   gravz(i,j,k) = const_grav

                endif
#endif

             enddo
          enddo
       enddo

    endif

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

#ifdef SELF_GRAVITY
             Sr_old = rhoo * gold(i,j,k,:)
#else
             Sr_old(dim) = rhoo * const_grav
#endif
             SrE_old = dot_product(vold, Sr_old)

             ! Define new source terms

             vnew = snew(UMX:UMZ) * rhoninv

#ifdef SELF_GRAVITY
             Sr_new = rhon * gnew(i,j,k,:)
#else
             Sr_new(dim) = rhon * const_grav
#endif
             SrE_new = dot_product(vnew, Sr_new)

             ! Define corrections to source terms

             Srcorr = HALF * (Sr_new - Sr_old)

             ! Correct momenta

             src(UMX:UMZ) = Srcorr

             snew(UMX:UMZ) = snew(UMX:UMZ) + dt * src(UMX:UMZ)

#ifdef HYBRID_MOMENTUM
             call add_hybrid_momentum_source(loc, src(UMR:UMP), Srcorr)

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

                ! The change in the gas energy is equal in magnitude to, and opposite in sign to,
                ! the change in the gravitational potential energy, rho * phi.
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

                SrEcorr = SrEcorr + hdtInv * ( flux1(i        ,j,k,URHO) * gravx(i        ,j,k) * dx(1) + &
                                               flux1(i+1*dg(1),j,k,URHO) * gravx(i+1*dg(1),j,k) * dx(1) + &
                                               flux2(i,j        ,k,URHO) * gravy(i,j        ,k) * dx(2) + &
                                               flux2(i,j+1*dg(2),k,URHO) * gravy(i,j+1*dg(2),k) * dx(2) + &
                                               flux3(i,j,k        ,URHO) * gravz(i,j,k        ) * dx(3) + &
                                               flux3(i,j,k+1*dg(3),URHO) * gravz(i,j,k+1*dg(3)) * dx(3) ) / vol(i,j,k)

             else

                call bl_error("Error:: gravity_sources_nd.F90 :: invalid grav_source_type")

             end if

             src(UEDEN) = SrEcorr

             snew(UEDEN) = snew(UEDEN) + dt * SrEcorr

             ! Add to the outgoing source array.

             source(i,j,k,:) = src

          enddo
       enddo
    enddo

    if (grav_source_type .eq. 4) then

       call bl_deallocate(gravx)
       call bl_deallocate(gravy)
       call bl_deallocate(gravz)

    endif

  end subroutine ca_corrgsrc

end module gravity_sources_module
