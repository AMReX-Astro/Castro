module advection_util_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_enforce_minimum_density(lo, hi, &
                                        state, s_lo, s_hi, &
                                        frac_change, verbose) bind(c,name='ca_enforce_minimum_density')

    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, small_dens, density_reset_method
    use amrex_constants_module, only: ZERO
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif
    use amrex_fort_module, only: rt => amrex_real
    use reduction_module, only: reduce_min

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: frac_change
    integer,  intent(in   ), value :: verbose

    ! Local variables
    integer  :: i, j, k
    integer  :: ii, jj, kk
    integer  :: i_set, j_set, k_set
    real(rt) :: max_dens, old_rho
    real(rt) :: uold(NVAR), unew(NVAR)
    integer  :: num_positive_zones

    !$gpu

    max_dens = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,URHO) .eq. ZERO) then

#ifndef AMREX_USE_GPU
                print *,'DENSITY EXACTLY ZERO AT CELL ', i, j, k
                print *,'  in grid ',lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
                call castro_error("Error :: ca_enforce_minimum_density")
#endif

             else if (state(i,j,k,URHO) < small_dens) then

                old_rho = state(i,j,k,URHO)

                if (density_reset_method == 1) then

                   ! Reset to the characteristics of the adjacent state with the highest density.

                   max_dens = state(i,j,k,URHO)
                   i_set = i
                   j_set = j
                   k_set = k
                   do kk = -1, 1
                      do jj = -1, 1
                         do ii = -1, 1

                            if (i+ii >= s_lo(1) .and. j+jj >= s_lo(2) .and. k+kk >= s_lo(3) .and. &
                                i+ii <= s_hi(1) .and. j+jj <= s_hi(2) .and. k+kk <= s_hi(3)) then

                               if (state(i+ii,j+jj,k+kk,URHO) .gt. max_dens) then

                                  i_set = i+ii
                                  j_set = j+jj
                                  k_set = k+kk
                                  max_dens = state(i_set,j_set,k_set,URHO)

                               end if

                            end if

                         end do
                      end do
                   end do

                   if (max_dens < small_dens) then

                      ! We could not find any nearby zones with sufficient density.

                      uold = state(i,j,k,:)
                      call reset_to_small_state(uold, [i, j, k], s_lo, s_hi, verbose)
                      state(i,j,k,:) = uold

                   else

                      uold = state(i,j,k,:)
                      unew = state(i_set,j_set,k_set,:)

                      call reset_to_zone_state(uold, unew, [i, j, k], s_lo, s_hi, verbose)

                      state(i,j,k,:) = uold

                   endif

                else if (density_reset_method == 2) then

                   ! Reset to the average of adjacent zones. The median is independently calculated for each variable.

                   num_positive_zones = 0
                   unew(:) = ZERO

                   do kk = -1, 1
                      do jj = -1, 1
                         do ii = -1, 1

                            if (i+ii >= s_lo(1) .and. j+jj >= s_lo(2) .and. k+kk >= s_lo(3) .and. &
                                i+ii <= s_hi(1) .and. j+jj <= s_hi(2) .and. k+kk <= s_hi(3)) then

                               if (state(i+ii,j+jj,k+kk,URHO) .ge. small_dens) then

                                  unew(:) = unew(:) + state(i+ii,j+jj,k+kk,:)
                                  num_positive_zones = num_positive_zones + 1

                               end if

                            end if

                         end do
                      end do
                   end do

                   if (num_positive_zones == 0) then

                      ! We could not find any nearby zones with sufficient density.

                      uold = state(i,j,k,:)
                      call reset_to_small_state(uold, [i, j, k], s_lo, s_hi, verbose)
                      state(i,j,k,:) = uold

                   else

                      uold = state(i,j,k,:)
                      unew(:) = unew(:) / num_positive_zones

                      call reset_to_zone_state(uold, unew, [i, j, k], s_lo, s_hi, verbose)

                      state(i,j,k,:) = uold

                   endif

#ifndef AMREX_USE_CUDA
                else

                   call castro_error("Unknown density_reset_method in subroutine ca_enforce_minimum_density.")
#endif
                endif

                ! Store the maximum (negative) fractional change in the density from this reset.

                if (old_rho < ZERO) then
                   call reduce_min(frac_change, (state(i,j,k,URHO) - old_rho) / old_rho)
                end if

             end if

          end do
       end do
    end do

  end subroutine ca_enforce_minimum_density


  subroutine reset_to_small_state(state, idx, lo, hi, verbose)
    ! If no neighboring zones are above small_dens, our only recourse
    ! is to set the density equal to small_dens, and the temperature
    ! equal to small_temp. We set the velocities to zero,
    ! though any choice here would be arbitrary.
    !

    use amrex_constants_module, only: ZERO
    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEINT, UEDEN, UFS, small_temp, small_dens, npassive, upass_map
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_module, only: eos
    use castro_util_module, only: position ! function
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only: linear_to_hybrid ! function
    use meth_params_module, only: UMR, UMP
#endif

    use amrex_fort_module, only: rt => amrex_real
    implicit none

    real(rt)         :: state(NVAR)
    integer          :: idx(3), lo(3), hi(3), verbose

    integer          :: n, ipassive
    type (eos_t)     :: eos_state

#ifdef HYBRID_MOMENTUM
    real(rt)         :: loc(3)
#endif

    !$gpu

#ifndef AMREX_USE_CUDA
    if (verbose .gt. 0) then
       print *,'   '
       if (state(URHO) < ZERO) then
          print *,'>>> RESETTING NEG.  DENSITY AT ', idx(1), idx(2), idx(3)
       else
          print *,'>>> RESETTING SMALL DENSITY AT ', idx(1), idx(2), idx(3)
       endif
       print *,'>>> FROM ', state(URHO), ' TO ', small_dens
       print *,'>>> IN GRID ', lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
       print *,'   '
    end if
#endif

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       state(n) = state(n) * (small_dens / state(URHO))
    end do

    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = state(UFS:UFS+nspec-1) / small_dens
    eos_state % aux = state(UFS:UFS+naux-1) / small_dens

    call eos(eos_input_rt, eos_state)

    state(URHO ) = eos_state % rho
    state(UTEMP) = eos_state % T

    state(UMX  ) = ZERO
    state(UMY  ) = ZERO
    state(UMZ  ) = ZERO

    state(UEINT) = eos_state % rho * eos_state % e
    state(UEDEN) = state(UEINT)

#ifdef HYBRID_MOMENTUM
    loc = position(idx(1),idx(2),idx(3))
    state(UMR:UMP) = linear_to_hybrid(loc, state(UMX:UMZ))
#endif

  end subroutine reset_to_small_state



  subroutine reset_to_zone_state(state, input_state, idx, lo, hi, verbose)

    use amrex_constants_module, only: ZERO
    use meth_params_module, only: NVAR, URHO
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt) :: state(NVAR), input_state(NVAR)
    integer  :: idx(3), lo(3), hi(3), verbose

    !$gpu

#ifndef AMREX_USE_CUDA
    if (verbose .gt. 0) then
       if (state(URHO) < ZERO) then
          print *,'   '
          print *,'>>> RESETTING NEG.  DENSITY AT ',idx(1), idx(2), idx(3)
          print *,'>>> FROM ', state(URHO) ,' TO ', input_state(URHO)
          print *,'>>> IN GRID ', lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
          print *,'   '
       else
          print *,'   '
          print *,'>>> RESETTING SMALL DENSITY AT ', idx(1), idx(2), idx(3)
          print *,'>>> FROM ', state(URHO), ' TO ', input_state(URHO)
          print *,'>>> IN GRID ', lo(1), lo(2), lo(3), hi(1), hi(2), hi(3)
          print *,'   '
       end if
    end if
#endif

    state(:) = input_state(:)

  end subroutine reset_to_zone_state


  subroutine ca_compute_cfl(lo, hi, &
                            q, q_lo, q_hi, &
                            qaux, qa_lo, qa_hi, &
                            dt, dx, courno, verbose) &
                            bind(C, name = "ca_compute_cfl")
    ! Compute running max of Courant number over grids
    !

    use amrex_constants_module, only: ZERO, ONE
    use meth_params_module, only: NQ, QRHO, QU, QV, QW, QC, NQAUX, time_integration_method
    use prob_params_module, only: dim
    use amrex_fort_module, only: rt => amrex_real
    use reduction_module, only: reduce_max

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)

    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in   ) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in   ), value :: dt
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: courno
    integer,  intent(in   ), value :: verbose

    real(rt) :: courx, coury, courz, courmx, courmy, courmz, courtmp
    real(rt) :: dtdx, dtdy, dtdz
    integer  :: i, j, k

    !$gpu

    courmx = courno
    courmy = courno
    courmz = courno

    dtdx = dt / dx(1)

    if (dim .ge. 2) then
       dtdy = dt / dx(2)
    else
       dtdy = ZERO
    endif

    if (dim .eq. 3) then
       dtdz = dt / dx(3)
    else
       dtdz = ZERO
    endif

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             courx = ( qaux(i,j,k,QC) + abs(q(i,j,k,QU)) ) * dtdx
             coury = ( qaux(i,j,k,QC) + abs(q(i,j,k,QV)) ) * dtdy
             courz = ( qaux(i,j,k,QC) + abs(q(i,j,k,QW)) ) * dtdz

             courmx = max( courmx, courx )
             courmy = max( courmy, coury )
             courmz = max( courmz, courz )

             if (time_integration_method == 0) then

                ! CTU integration constraint

#ifndef AMREX_USE_CUDA
                if (verbose == 1) then

                   if (courx .gt. ONE) then
                      print *,'   '
                      print *, "Warning:: advection_util_nd.F90 :: CFL violation in compute_cfl"
                      print *,'>>> ... (u+c) * dt / dx > 1 ', courx
                      print *,'>>> ... at cell (i,j,k)   : ', i, j, k
                      print *,'>>> ... u, c                ', q(i,j,k,QU), qaux(i,j,k,QC)
                      print *,'>>> ... density             ', q(i,j,k,QRHO)
                   end if

                   if (coury .gt. ONE) then
                      print *,'   '
                      print *, "Warning:: advection_util_nd.F90 :: CFL violation in compute_cfl"
                      print *,'>>> ... (v+c) * dt / dx > 1 ', coury
                      print *,'>>> ... at cell (i,j,k)   : ', i,j,k
                      print *,'>>> ... v, c                ', q(i,j,k,QV), qaux(i,j,k,QC)
                      print *,'>>> ... density             ', q(i,j,k,QRHO)
                   end if

                   if (courz .gt. ONE) then
                      print *,'   '
                      print *, "Warning:: advection_util_nd.F90 :: CFL violation in compute_cfl"
                      print *,'>>> ... (w+c) * dt / dx > 1 ', courz
                      print *,'>>> ... at cell (i,j,k)   : ', i, j, k
                      print *,'>>> ... w, c                ', q(i,j,k,QW), qaux(i,j,k,QC)
                      print *,'>>> ... density             ', q(i,j,k,QRHO)
                   end if

                end if
#endif

             else

                ! method-of-lines constraint
                courtmp = courx
                if (dim >= 2) then
                   courtmp = courtmp + coury
                endif
                if (dim == 3) then
                   courtmp = courtmp + courz
                endif

#ifndef AMREX_USE_CUDA
                if (verbose == 1) then

                   ! note: it might not be 1 for all RK integrators
                   if (courtmp > ONE) then
                      print *,'   '
                      print *, "Warning:: advection_util_nd.F90 :: CFL violation in compute_cfl"
                      print *,'>>> ... at cell (i,j,k)   : ', i, j, k
                      print *,'>>> ... u,v,w, c            ', q(i,j,k,QU), q(i,j,k,QV), q(i,j,k,QW), qaux(i,j,k,QC)
                      print *,'>>> ... density             ', q(i,j,k,QRHO)
                   endif

                end if
#endif

                call reduce_max(courno, courtmp)
             endif
          enddo
       enddo
    enddo

    if (time_integration_method == 0) then
       call reduce_max(courno, max(courmx, courmy, courmz))
    endif

  end subroutine ca_compute_cfl



  subroutine ca_ctoprim(lo, hi, &
                        uin, uin_lo, uin_hi, &
#ifdef RADIATION
                        Erin, Erin_lo, Erin_hi, &
                        lam, lam_lo, lam_hi, &
#endif
                        q,     q_lo,   q_hi, &
                        qaux, qa_lo,  qa_hi) bind(c,name='ca_ctoprim')

    use actual_network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, &
         UEDEN, UEINT, UTEMP, &
         QRHO, QU, QV, QW, &
         QREINT, QPRES, QTEMP, QGAME, QFS, QFX, &
         NQ, QC, QGAMC, QGC, QDPDR, QDPDE, NQAUX, &
#ifdef RADIATION
         QCG, QGAMCG, QLAMS, &
         QPTOT, QRAD, QRADHI, QREITOT, &
#endif
         npassive, upass_map, qpass_map, dual_energy_eta1, &
         small_dens

    use amrex_constants_module, only: ZERO, HALF, ONE
    use castro_error_module
#ifdef ROTATION
    use meth_params_module, only: do_rotation, state_in_rotating_frame
    use rotation_module, only: inertial_to_rotational_velocity
    use amrinfo_module, only: amr_time
#endif
#ifdef RADIATION
    use rad_params_module, only: ngroups
    use rad_util_module, only: compute_ptot_ctot
#endif

    use amrex_fort_module, only: rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
#ifdef RADIATION
    integer, intent(in) :: Erin_lo(3), Erin_hi(3)
    integer, intent(in) :: lam_lo(3), lam_hi(3)
#endif
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
#ifdef RADIATION
    real(rt), intent(in   ) :: Erin(Erin_lo(1):Erin_hi(1),Erin_lo(2):Erin_hi(2),Erin_lo(3):Erin_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
#endif

    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), parameter :: small = 1.e-8_rt

    integer  :: i, j, k, g
    integer  :: n, iq, ipassive
    real(rt) :: kineng, rhoinv
    real(rt) :: vel(3)

    type (eos_t) :: eos_state

#ifdef RADIATION
    real(rt) :: ptot, ctot, gamc_tot
    real(rt) :: lams(0:ngroups-1), qs(NQ)
#endif

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)

#ifndef AMREX_USE_CUDA
          do i = lo(1), hi(1)
             if (uin(i,j,k,URHO) .le. ZERO) then
                print *,'   '
                print *,'>>> Error: advection_util_nd.F90::ctoprim ',i, j, k
                print *,'>>> ... negative density ', uin(i,j,k,URHO)
                call castro_error("Error:: advection_util_nd.f90 :: ctoprim")
             else if (uin(i,j,k,URHO) .lt. small_dens) then
                print *,'   '
                print *,'>>> Error: advection_util_nd.F90::ctoprim ',i, j, k
                print *,'>>> ... small density ', uin(i,j,k,URHO)
                call castro_error("Error:: advection_util_nd.f90 :: ctoprim")
             endif
          end do
#endif
          do i = lo(1), hi(1)

             q(i,j,k,QRHO) = uin(i,j,k,URHO)
             rhoinv = ONE/q(i,j,k,QRHO)

             vel = uin(i,j,k,UMX:UMZ) * rhoinv

             q(i,j,k,QU:QW) = vel

             ! Get the internal energy, which we'll use for
             ! determining the pressure.  We use a dual energy
             ! formalism. If (E - K) < eta1 and eta1 is suitably
             ! small, then we risk serious numerical truncation error
             ! in the internal energy.  Therefore we'll use the result
             ! of the separately updated internal energy equation.
             ! Otherwise, we'll set e = E - K.

             kineng = HALF * q(i,j,k,QRHO) * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)

             if ( (uin(i,j,k,UEDEN) - kineng) / uin(i,j,k,UEDEN) .gt. dual_energy_eta1) then
                q(i,j,k,QREINT) = (uin(i,j,k,UEDEN) - kineng) * rhoinv
             else
                q(i,j,k,QREINT) = uin(i,j,k,UEINT) * rhoinv
             endif

             ! If we're advecting in the rotating reference frame,
             ! then subtract off the rotation component here.

#ifdef ROTATION
             if (do_rotation == 1 .and. state_in_rotating_frame /= 1) then
                vel = q(i,j,k,QU:QW)
                call inertial_to_rotational_velocity([i, j, k], amr_time, vel)
                q(i,j,k,QU:QW) = vel
             endif
#endif

             q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)
#ifdef RADIATION
             q(i,j,k,qrad:qradhi) = Erin(i,j,k,:)
#endif

             ! Load passively advected quatities into q
             do ipassive = 1, npassive
                n  = upass_map(ipassive)
                iq = qpass_map(ipassive)
                q(i,j,k,iq) = uin(i,j,k,n) * rhoinv
             enddo

             ! get gamc, p, T, c, csml using q state
             eos_state % T   = q(i,j,k,QTEMP )
             eos_state % rho = q(i,j,k,QRHO  )
             eos_state % e   = q(i,j,k,QREINT)
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)
             eos_state % aux = q(i,j,k,QFX:QFX+naux-1)

             call eos(eos_input_re, eos_state)

             q(i,j,k,QTEMP)  = eos_state % T
             q(i,j,k,QREINT) = eos_state % e * q(i,j,k,QRHO)
             q(i,j,k,QPRES)  = eos_state % p
             q(i,j,k,QGAME)  = q(i,j,k,QPRES) / q(i,j,k,QREINT) + ONE
             q(i,j,k,QGC) = eos_state % gam1

             qaux(i,j,k,QDPDR)  = eos_state % dpdr_e
             qaux(i,j,k,QDPDE)  = eos_state % dpde

#ifdef RADIATION
             qaux(i,j,k,QGAMCG)   = eos_state % gam1
             qaux(i,j,k,QCG)      = eos_state % cs

             lams(:) = lam(i,j,k,:)
             qs(:) = q(i,j,k,:)
             call compute_ptot_ctot(lams, qs, qaux(i,j,k,QCG), &
                                    ptot, ctot, gamc_tot)

             q(i,j,k,QPTOT) = ptot

             qaux(i,j,k,QC)    = ctot
             qaux(i,j,k,QGAMC) = gamc_tot

             do g = 0, ngroups-1
                qaux(i,j,k,QLAMS+g) = lam(i,j,k,g)
             enddo

             q(i,j,k,qreitot) = q(i,j,k,QREINT) + sum(q(i,j,k,qrad:qradhi))
#else
             qaux(i,j,k,QGAMC)  = eos_state % gam1
             qaux(i,j,k,QC   )  = eos_state % cs
#endif

          enddo
       enddo
    enddo

  end subroutine ca_ctoprim


  function dflux(u, q, dir, idx) result(flux)
    ! Given a conservative state and its corresponding primitive state, calculate the
    ! corresponding flux in a given direction.
    !

    use amrex_constants_module, only: ZERO
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, &
         NQ, QU, QPRES, &
         npassive, upass_map
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only: compute_hybrid_flux
    use meth_params_module, only: NGDNV, GDRHO, GDU, GDW, GDPRES, QRHO, QW
#endif
    use prob_params_module, only: mom_flux_has_p
    use amrex_fort_module, only: rt => amrex_real
    implicit none

    integer :: dir, idx(3)
    real(rt)         :: u(NVAR), q(NQ), flux(NVAR)

    real(rt)         :: v_adv
    integer :: ipassive, n
#ifdef HYBRID_MOMENTUM
    real(rt)         :: qgdnv(NGDNV)
    logical :: cell_centered
#endif

    !$gpu

    ! Set everything to zero; this default matters because some
    ! quantities like temperature are not updated through fluxes.

    flux = ZERO

    ! Determine the advection speed based on the flux direction.

    v_adv = q(QU + dir - 1)

    ! Core quantities (density, momentum, energy).

    flux(URHO) = u(URHO) * v_adv
    flux(UMX:UMZ) = u(UMX:UMZ) * v_adv
    flux(UEDEN) = (u(UEDEN) + q(QPRES)) * v_adv
    flux(UEINT) = u(UEINT) * v_adv

    ! Optionally include the pressure term in the momentum flux.
    ! It is optional because for some geometries we cannot write
    ! the pressure term in a conservative form.

    if (mom_flux_has_p(dir)%comp(UMX+dir-1)) then
       flux(UMX + dir - 1) = flux(UMX + dir - 1) + q(QPRES)
    endif

    ! Hybrid flux.

#ifdef HYBRID_MOMENTUM
    ! Create a temporary edge-based q for this routine.
    qgdnv(:) = ZERO
    qgdnv(GDRHO) = q(QRHO)
    qgdnv(GDU:GDW) = q(QU:QW)
    qgdnv(GDPRES) = q(QPRES)
    cell_centered = .true.
    call compute_hybrid_flux(qgdnv, flux, dir, idx, cell_centered)
#endif

    ! Passively advected quantities.

    do ipassive = 1, npassive

       n = upass_map(ipassive)
       flux(n) = u(n) * v_adv

    enddo

  end function dflux


  subroutine limit_hydro_fluxes_on_small_dens(lo, hi, &
                                              idir, &
                                              u, u_lo, u_hi, &
                                              q, q_lo, q_hi, &
                                              vol, vol_lo, vol_hi, &
                                              flux, flux_lo, flux_hi, &
                                              area, area_lo, area_hi, &
                                              dt, dx) bind(c, name="limit_hydro_fluxes_on_small_dens")
    ! The following algorithm comes from Hu, Adams, and Shu (2013), JCP, 242, 169,
    ! "Positivity-preserving method for high-order conservative schemes solving
    ! compressible Euler equations." It has been modified to enforce not only positivity
    ! but also the stronger requirement that rho > small_dens. We do not limit on pressure
    ! (or, similarly, internal energy) because those cases are easily fixed by calls to
    ! reset_internal_energy that enforce a thermodynamic floor. The density limiter, by
    ! contrast, is very important because calls to enforce_minimum_density can yield
    ! hydrodynamic states that are inconsistent (there is no clear strategy for what to do
    ! when a density is negative).

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE, TWO
    use meth_params_module, only: NVAR, NQ, URHO, UTEMP, USHK, small_dens, cfl
    use prob_params_module, only: dim
    use amrex_mempool_module, only: bl_allocate, bl_deallocate

    implicit none

    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in), value :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: flux_lo(3), flux_hi(3)
    integer, intent(in) :: area_lo(3), area_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: flux(flux_lo(1):flux_hi(1),flux_lo(2):flux_hi(2),flux_lo(3):flux_hi(3),NVAR)
    real(rt), intent(in   ) :: area(area_lo(1):area_hi(1),area_lo(2):area_hi(2),area_lo(3):area_hi(3))

    integer  :: i, j, k

    real(rt) :: rhoL, rhoR, drhoL, drhoR, fluxLF(NVAR), fluxL(NVAR), fluxR(NVAR), rhoLF, drhoLF, dtdx, theta, alpha
    real(rt) :: uL(NVAR), uR(NVAR), qL(NQ), qR(NQ), volL, volR, flux_coefL, flux_coefR
    integer  :: idxL(3), idxR(3)

    real(rt), parameter :: density_floor_tolerance = 1.1_rt
    real(rt) :: density_floor

    !$gpu

    ! The density floor is the small density, modified by a small factor.
    ! In practice numerical error can cause the density that is created
    ! by this flux limiter to be slightly lower than the target density,
    ! so we set the target to be slightly larger than the real density floor
    ! to avoid density resets.

    density_floor = small_dens * density_floor_tolerance

    dtdx = dt / dx(idir)
    alpha = ONE / DIM

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! Grab the states on either side of the interface we are working with,
             ! depending on which dimension we're currently calling this with.

             uR = u(i,j,k,:)
             qR = q(i,j,k,:)
             volR = vol(i,j,k)
             idxR = [i,j,k]

             if (idir == 1) then
                uL = u(i-1,j,k,:)
                qL = q(i-1,j,k,:)
                volL = vol(i-1,j,k)
                idxL = [i-1,j,k]
             else if (idir == 2) then
                uL = u(i,j-1,k,:)
                qL = q(i,j-1,k,:)
                volL = vol(i,j-1,k)
                idxL = [i,j-1,k]
             else
                uL = u(i,j,k-1,:)
                qL = q(i,j,k-1,:)
                volL = vol(i,j,k-1)
                idxL = [i,j,k-1]
             end if

             ! If an adjacent zone has a floor-violating density, set the flux to zero and move on.
             ! At that point, the only thing to do is wait for a reset at a later point.

             if (uR(URHO) < density_floor .or. uL(URHO) < density_floor) then

                flux(i,j,k,:) = ZERO
                cycle

             endif

             ! Construct cell-centered fluxes.

             fluxL = dflux(uL, qL, idir, idxL)
             fluxR = dflux(uR, qR, idir, idxR)

             ! Construct the Lax-Friedrichs flux on the interface (Equation 12).
             ! Note that we are using the information from Equation 9 to obtain the
             ! effective maximum wave speed, (|u| + c)_max = CFL / lambda where
             ! lambda = dt/(dx * alpha); alpha = 1 in 1D and may be chosen somewhat
             ! freely in multi-D as long as alpha_x + alpha_y + alpha_z = 1.

             fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha) * (uL(:) - uR(:)))

             ! Coefficients of fluxes on either side of the interface.

             flux_coefR = TWO * (dt / alpha) * area(i,j,k) / volR
             flux_coefL = TWO * (dt / alpha) * area(i,j,k) / volL

             ! Obtain the one-sided update to the density, based on Hu et al., Eq. 11.
             ! If we would violate the floor, then we need to limit the flux. Since the
             ! flux adds to the density on one side and subtracts from the other, the floor
             ! can only be violated in at most one direction, so we'll do an if-else test
             ! below. This means that we can simplify the approach of Hu et al. -- whereas
             ! they constructed two thetas for each interface (corresponding to either side)
             ! we can complete the operation in one step with a single theta.

             drhoL = flux_coefL * flux(i,j,k,URHO)
             rhoL = uL(URHO) - drhoL

             drhoR = flux_coefR * flux(i,j,k,URHO)
             rhoR = uR(URHO) + drhoR

             theta = ONE

             if (rhoL < density_floor) then

                ! Obtain the final density corresponding to the LF flux.

                drhoLF = flux_coefL * fluxLF(URHO)
                rhoLF = uL(URHO) - drhoLF

                ! Solve for theta from (1 - theta) * rhoLF + theta * rho = density_floor.

                theta = (density_floor - rhoLF) / (rhoL - rhoLF)

                ! Limit theta to the valid range (this will deal with roundoff issues).

                theta = min(ONE, max(theta, ZERO))

             else if (rhoR < density_floor) then

                drhoLF = flux_coefR * fluxLF(URHO)
                rhoLF = uR(URHO) + drhoLF

                theta = (density_floor - rhoLF) / (rhoR - rhoLF)

                theta = min(ONE, max(theta, ZERO))

             endif

             ! Assemble the limited flux (Equation 16).

             flux(i,j,k,:) = (ONE - theta) * fluxLF(:) + theta * flux(i,j,k,:)

             ! Zero out fluxes for quantities that don't advect.

             flux(i,j,k,UTEMP) = ZERO
#ifdef SHOCK_VAR
             flux(i,j,k,USHK) = ZERO
#endif

             ! Now, apply our requirement that the final flux cannot violate the density floor.

             drhoR = flux_coefR * flux(i,j,k,URHO)
             drhoL = flux_coefL * flux(i,j,k,URHO)

             if (uR(URHO) + drhoR < density_floor) then
                flux(i,j,k,:) = flux(i,j,k,:) * abs((density_floor - uR(URHO)) / drhoR)
             else if (uL(URHO) - drhoL < density_floor) then
                flux(i,j,k,:) = flux(i,j,k,:) * abs((density_floor - uL(URHO)) / drhoL)
             endif

          enddo
       enddo
    enddo

  end subroutine limit_hydro_fluxes_on_small_dens



  subroutine limit_hydro_fluxes_on_large_vel(lo, hi, &
                                             idir, &
                                             u, u_lo, u_hi, &
                                             q, q_lo, q_hi, &
                                             vol, vol_lo, vol_hi, &
                                             flux, flux_lo, flux_hi, &
                                             area, area_lo, area_hi, &
                                             dt, dx) bind(c, name="limit_hydro_fluxes_on_large_vel")
    ! This limiter is similar to the density-based limiter above, but limits
    ! on velocities that are too large instead. The comments are minimal since
    ! the algorithm is effectively the same.

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE, TWO
    use meth_params_module, only: NVAR, NQ, URHO, UMX, UTEMP, USHK, cfl, speed_limit
    use prob_params_module, only: dim

    implicit none

    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in), value :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: flux_lo(3), flux_hi(3)
    integer, intent(in) :: area_lo(3), area_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: flux(flux_lo(1):flux_hi(1),flux_lo(2):flux_hi(2),flux_lo(3):flux_hi(3),NVAR)
    real(rt), intent(in   ) :: area(area_lo(1):area_hi(1),area_lo(2):area_hi(2),area_lo(3):area_hi(3))

    integer  :: i, j, k, n

    real(rt) :: fluxLF(NVAR), fluxL(NVAR), fluxR(NVAR), dtdx, theta, alpha
    real(rt) :: rhouL, rhouR, drhouL, drhouR, rhouLF, drhouLF
    real(rt) :: rhoL, rhoR, drhoL, drhoR
    real(rt) :: uL(NVAR), uR(NVAR), qL(NQ), qR(NQ), volL, volR, flux_coefL, flux_coefR
    integer  :: idxL(3), idxR(3), UMOM

    real(rt) :: momentum_ceiling

    !$gpu

    dtdx = dt / dx(idir)
    alpha = ONE / DIM

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             uR = u(i,j,k,:)
             qR = q(i,j,k,:)
             volR = vol(i,j,k)
             idxR = [i,j,k]

             if (idir == 1) then
                uL = u(i-1,j,k,:)
                qL = q(i-1,j,k,:)
                volL = vol(i-1,j,k)
                idxL = [i-1,j,k]
             else if (idir == 2) then
                uL = u(i,j-1,k,:)
                qL = q(i,j-1,k,:)
                volL = vol(i,j-1,k)
                idxL = [i,j-1,k]
             else
                uL = u(i,j,k-1,:)
                qL = q(i,j,k-1,:)
                volL = vol(i,j,k-1)
                idxL = [i,j,k-1]
             end if

             ! Construct cell-centered fluxes.

             fluxL = dflux(uL, qL, idir, idxL)
             fluxR = dflux(uR, qR, idir, idxR)

             ! Construct the Lax-Friedrichs flux on the interface.

             fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx / alpha) * (uL(:) - uR(:)))

             ! Coefficients of fluxes on either side of the interface.

             flux_coefR = TWO * (dt / alpha) * area(i,j,k) / volR
             flux_coefL = TWO * (dt / alpha) * area(i,j,k) / volL

             theta = ONE

             ! Loop over all three momenta, and choose the strictest
             ! limiter among them.

             do n = 1, 3

                UMOM = UMX + n - 1

                ! Obtain the one-sided update to the momentum.

                drhouL = flux_coefL * flux(i,j,k,UMOM)
                rhouL = abs(uL(UMOM) - drhouL)

                drhoL = flux_coefL * flux(i,j,k,URHO)
                rhoL = uL(URHO) - drhoL

                drhouR = flux_coefR * flux(i,j,k,UMOM)
                rhouR = abs(uR(UMOM) + drhouR)

                drhoR = flux_coefR * flux(i,j,k,URHO)
                rhoR = uR(uRHO) + drhoR

                if (abs(rhouL) > rhoL * speed_limit) then

                   ! Obtain the final density corresponding to the LF flux.

                   drhouLF = flux_coefL * fluxLF(UMOM)
                   rhouLF = abs(uL(UMOM) - drhouLF)

                   ! Solve for theta from (1 - theta) * rhouLF + theta * rhou = momentum_ceiling.

                   theta = (momentum_ceiling - rhouLF) / (rhouL - rhouLF)

                   ! Limit theta to the valid range (this will deal with roundoff issues).

                   theta = min(ONE, max(theta, ZERO))

                else if (abs(rhouR) > rhoR * speed_limit) then

                   drhouLF = flux_coefR * fluxLF(UMOM)
                   rhouLF = abs(uR(UMOM) + drhouLF)

                   theta = (momentum_ceiling - abs(rhouLF)) / (abs(rhouR - rhouLF))

                   theta = min(ONE, max(theta, ZERO))

                endif

             end do

             ! Assemble the limited flux (Equation 16).

             flux(i,j,k,:) = (ONE - theta) * fluxLF(:) + theta * flux(i,j,k,:)

             ! Zero out fluxes for quantities that don't advect.

             flux(i,j,k,UTEMP) = ZERO
#ifdef SHOCK_VAR
             flux(i,j,k,USHK) = ZERO
#endif

             ! Now, apply our requirement that the final flux cannot violate the momentum ceiling.

             do n = 1, 3

                UMOM = UMX + n - 1

                drhouR = flux_coefR * flux(i,j,k,UMOM)
                drhouL = flux_coefL * flux(i,j,k,UMOM)

                drhoR = flux_coefR * flux(i,j,k,URHO)
                drhoL = flux_coefL * flux(i,j,k,URHO)

                if (abs(uR(UMOM) + drhouR) > (uR(URHO) + drhoR) * speed_limit) then
                   flux(i,j,k,:) = flux(i,j,k,:) * abs((momentum_ceiling - abs(uR(UMOM))) / drhouR)
                else if (abs(uL(UMOM) - drhouL) > (uL(URHO) - drhoL) * speed_limit) then
                   flux(i,j,k,:) = flux(i,j,k,:) * abs((momentum_ceiling - abs(uL(UMOM))) / drhouL)
                endif

             end do

          enddo
       enddo
    enddo

  end subroutine limit_hydro_fluxes_on_large_vel



  subroutine ca_shock(lo, hi, &
                      q, qd_lo, qd_hi, &
                      shk, s_lo, s_hi, &
                      dx) bind(C, name="ca_shock")
    ! This is a basic multi-dimensional shock detection algorithm.
    ! This implementation follows Flash, which in turn follows
    ! AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
    !
    ! The spirit of this follows the shock detection in Colella &
    ! Woodward (1984)
    !

    use meth_params_module, only: QPRES, QU, QV, QW, NQ
    use prob_params_module, only: coord_type
    use amrex_constants_module, only: ZERO, HALF, ONE
    use castro_error_module
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(inout) :: shk(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    integer :: i, j, k

    real(rt) :: dxinv, dyinv, dzinv
    real(rt) :: div_u
    real(rt) :: px_pre, px_post, py_pre, py_post, pz_pre, pz_post
    real(rt) :: e_x, e_y, e_z, d
    real(rt) :: p_pre, p_post, pjump

    real(rt), parameter :: small = 1.e-10_rt
    real(rt), parameter :: eps = 0.33e0_rt

    real(rt) :: rm, rc, rp

    !$gpu

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

#ifndef AMREX_USE_CUDA
    if (coord_type /= 0) then
       call castro_error("ERROR: invalid geometry in shock()")
    endif
#endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! construct div{U}
             if (coord_type == 0) then
                ! Cartesian
                div_u = HALF*(q(i+1,j,k,QU) - q(i-1,j,k,QU))*dxinv
#if (AMREX_SPACEDIM >= 2)
                div_u = div_u + HALF*(q(i,j+1,k,QV) - q(i,j-1,k,QV))*dyinv
#endif
#if (AMREX_SPACEDIM == 3)
                div_u = div_u + HALF*(q(i,j,k+1,QW) - q(i,j,k-1,QW))*dzinv
#endif
             elseif (coord_type == 1) then
                ! r-z
                rc = dble(i + HALF)*dx(1)
                rm = dble(i - 1 + HALF)*dx(1)
                rp = dble(i + 1 + HALF)*dx(1)

#if (AMREX_SPACEDIM == 1)
                div_u = HALF*(rp*q(i+1,j,k,QU) - rm*q(i-1,j,k,QU))/(rc*dx(1))
#endif
#if (AMREX_SPACEDIM == 2)
                div_u = HALF*(rp*q(i+1,j,k,QU) - rm*q(i-1,j,k,QU))/(rc*dx(1)) + &
                     HALF*(q(i,j+1,k,QV) - q(i,j-1,k,QV)) * dyinv
#endif

             elseif (coord_type == 2) then
                ! 1-d spherical
                rc = dble(i + HALF)*dx(1)
                rm = dble(i - 1 + HALF)*dx(1)
                rp = dble(i + 1 + HALF)*dx(1)

                div_u = HALF*(rp**2*q(i+1,j,k,QU) - rm**2*q(i-1,j,k,QU))/(rc**2*dx(1))

#ifndef AMREX_USE_CUDA
             else
                call castro_error("ERROR: invalid coord_type in shock")
#endif
             endif


             ! find the pre- and post-shock pressures in each direction
             if (q(i+1,j,k,QPRES) - q(i-1,j,k,QPRES) < ZERO) then
                px_pre  = q(i+1,j,k,QPRES)
                px_post = q(i-1,j,k,QPRES)
             else
                px_pre  = q(i-1,j,k,QPRES)
                px_post = q(i+1,j,k,QPRES)
             endif

#if (AMREX_SPACEDIM >= 2)
             if (q(i,j+1,k,QPRES) - q(i,j-1,k,QPRES) < ZERO) then
                py_pre  = q(i,j+1,k,QPRES)
                py_post = q(i,j-1,k,QPRES)
             else
                py_pre  = q(i,j-1,k,QPRES)
                py_post = q(i,j+1,k,QPRES)
             endif
#else
             py_pre = ZERO
             py_post = ZERO
#endif

#if (AMREX_SPACEDIM == 3)
             if (q(i,j,k+1,QPRES) - q(i,j,k-1,QPRES) < ZERO) then
                pz_pre  = q(i,j,k+1,QPRES)
                pz_post = q(i,j,k-1,QPRES)
             else
                pz_pre  = q(i,j,k-1,QPRES)
                pz_post = q(i,j,k+1,QPRES)
             endif
#else
             pz_pre = ZERO
             pz_post = ZERO
#endif

             ! use compression to create unit vectors for the shock direction
             e_x = (q(i+1,j,k,QU) - q(i-1,j,k,QU))**2
#if (AMREX_SPACEDIM >= 2)
             e_y = (q(i,j+1,k,QV) - q(i,j-1,k,QV))**2
#else
             e_y = ZERO
#endif
#if (AMREX_SPACEDIM == 3)
             e_z = (q(i,j,k+1,QW) - q(i,j,k-1,QW))**2
#else
             e_z = ZERO
#endif
             d = ONE/(e_x + e_y + e_z + small)

             e_x = e_x*d
             e_y = e_y*d
             e_z = e_z*d

             ! project the pressures onto the shock direction
             p_pre  = e_x*px_pre + e_y*py_pre + e_z*pz_pre
             p_post = e_x*px_post + e_y*py_post + e_z*pz_post

             ! test for compression + pressure jump to flag a shock
             if (p_pre == ZERO) then
                ! this can happen if U = 0, so e_x, ... = 0
                pjump = ZERO
             else
                pjump = eps - (p_post - p_pre)/p_pre
             endif

             if (pjump < ZERO .and. div_u < ZERO) then
                shk(i,j,k) = ONE
             else
                shk(i,j,k) = ZERO
             endif

          enddo
       enddo
    enddo

  end subroutine ca_shock

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::


  subroutine divu(lo, hi, &
                  q, q_lo, q_hi, &
                  dx, div, div_lo, div_hi) bind(C, name='divu')
    ! this computes the *node-centered* divergence
    !

    use meth_params_module, only: QU, QV, QW, NQ
    use amrex_constants_module, only: HALF, FOURTH, ONE, ZERO
    use prob_params_module, only: dg, coord_type, problo
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: div_lo(3), div_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

    integer  :: i, j, k
    real(rt) :: ux, vy, wz, dxinv, dyinv, dzinv
    real(rt) :: rl, rr, rc, ul, ur, vt, vb

    !$gpu

    dxinv = ONE/dx(1)
#if AMREX_SPACEDIM >= 2
    dyinv = ONE/dx(2)
#else
    dyinv = ZERO
#endif
#if AMREX_SPACEDIM == 3
    dzinv = ONE/dx(3)
#else
    dzinv = ZERO
#endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

#if AMREX_SPACEDIM == 1
             if (coord_type == 0) then
                div(i,j,k) = (q(i,j,k,QU) - q(i-1*dg(1),j,k,QU)) * dxinv

             else if (coord_type == 1) then
                ! axisymmetric
                if (i == 0) then
                   div(i,j,k) = ZERO
                else
                   rl = (dble(i)-HALF) * dx(1) + problo(1)
                   rr = (dble(i)+HALF) * dx(1) + problo(1)
                   rc = (dble(i)     ) * dx(1) + problo(1)

                   div(i,j,k) = (rr*q(i,j,k,QU) - rl*q(i-1*dg(1),j,k,QU)) * dxinv / rc
                endif
             else
                ! spherical
                if (i == 0) then
                   div(i,j,k) = ZERO
                else
                   rl = (dble(i)-HALF) * dx(1) + problo(1)
                   rr = (dble(i)+HALF) * dx(1) + problo(1)
                   rc = (dble(i)     ) * dx(1) + problo(1)

                   div(i,j,k) = (rr**2*q(i,j,k,QU) - rl**2*q(i-1*dg(1),j,k,QU)) * dxinv / rc**2
                endif

             endif

#endif

#if AMREX_SPACEDIM == 2
             if (coord_type == 0) then
                ux = HALF*(q(i,j,k,QU) - q(i-1*dg(1),j        ,k,QU) + q(i        ,j-1*dg(2),k,QU) - q(i-1*dg(1),j-1*dg(2),k,QU)) * dxinv
                vy = HALF*(q(i,j,k,QV) - q(i        ,j-1*dg(2),k,QV) + q(i-1*dg(1),j        ,k,QV) - q(i-1*dg(1),j-1*dg(2),k,QV)) * dyinv

             else

                if (i == 0) then
                   ux = ZERO
                   vy = ZERO   ! is this part correct?

                else
                   rl = (dble(i)-HALF) * dx(1) + problo(1)
                   rr = (dble(i)+HALF) * dx(1) + problo(1)
                   rc = (dble(i)     ) * dx(1) + problo(1)

                   ! These are transverse averages in the y-direction
                   ul = HALF * (q(i-1*dg(1),j,k,QU) + q(i-1*dg(1),j-1*dg(2),k,QU))
                   ur = HALF * (q(i        ,j,k,QU) + q(i        ,j-1*dg(2),k,QU))

                   ! Take 1/r d/dr(r*u)
                   ux = (rr*ur - rl*ul) * dxinv / rc

                   ! These are transverse averages in the x-direction
                   vb = HALF * (q(i,j-1*dg(2),k,QV) + q(i-1*dg(1),j-1*dg(2),k,QV))
                   vt = HALF * (q(i,j        ,k,QV) + q(i-1*dg(1),j        ,k,QV))

                   vy = (vt - vb) * dyinv

                end if

             endif

             div(i,j,k) = ux + vy
#endif

#if AMREX_SPACEDIM == 3
             ux = FOURTH*( &
                  + q(i        ,j        ,k        ,QU) - q(i-1*dg(1),j        ,k        ,QU) &
                  + q(i        ,j        ,k-1*dg(3),QU) - q(i-1*dg(1),j        ,k-1*dg(3),QU) &
                  + q(i        ,j-1*dg(2),k        ,QU) - q(i-1*dg(1),j-1*dg(2),k        ,QU) &
                  + q(i        ,j-1*dg(2),k-1*dg(3),QU) - q(i-1*dg(1),j-1*dg(2),k-1*dg(3),QU) ) * dxinv

             vy = FOURTH*( &
                  + q(i        ,j        ,k        ,QV) - q(i        ,j-1*dg(2),k        ,QV) &
                  + q(i        ,j        ,k-1*dg(3),QV) - q(i        ,j-1*dg(2),k-1*dg(3),QV) &
                  + q(i-1*dg(1),j        ,k        ,QV) - q(i-1*dg(1),j-1*dg(2),k        ,QV) &
                  + q(i-1*dg(1),j        ,k-1*dg(3),QV) - q(i-1*dg(1),j-1*dg(2),k-1*dg(3),QV) ) * dyinv

             wz = FOURTH*( &
                  + q(i        ,j        ,k        ,QW) - q(i        ,j        ,k-1*dg(3),QW) &
                  + q(i        ,j-1*dg(2),k        ,QW) - q(i        ,j-1*dg(2),k-1*dg(3),QW) &
                  + q(i-1*dg(1),j        ,k        ,QW) - q(i-1*dg(1),j        ,k-1*dg(3),QW) &
                  + q(i-1*dg(1),j-1*dg(2),k        ,QW) - q(i-1*dg(1),j-1*dg(2),k-1*dg(3),QW) ) * dzinv

             div(i,j,k) = ux + vy + wz
#endif

          enddo
       enddo
    enddo

  end subroutine divu


  subroutine avisc(lo, hi, &
                   q, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   dx, avis, a_lo, a_hi, idir)
    ! this computes the *face-centered* artifical viscosity using the
    ! 4th order expression from McCorquodale & Colella (Eq. 35)

    use meth_params_module, only: QU, QV, QW, QC, NQ, NQAUX
    use amrex_constants_module, only: HALF, FOURTH, ONE, ZERO
    use prob_params_module, only: dg
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: a_lo(3), a_hi(3)
    integer, intent(in) :: idir
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(inout) :: avis(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    real(rt) :: coeff, cmin

    real(rt), parameter :: beta = 0.3_rt

    integer :: i, j, k

    real(rt) :: dxinv, dyinv, dzinv

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    do k = lo(3), hi(3)+dg(3)
       do j = lo(2), hi(2)+dg(2)
          do i = lo(1), hi(1)+1

             if (idir == 1) then

                ! normal direction
                avis(i,j,k) = (q(i,j,k,QU) - q(i-1,j,k,QU)) * dxinv
#if BL_SPACEDIM >= 2
                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i,j+1,k,QV) - q(i,j-1,k,QV) + &
                     q(i-1,j+1,k,QV) - q(i-1,j-1,k,QV)) * dyinv
#endif
#if BL_SPACEDIM >= 3
                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i,j,k+1,QW) - q(i,j,k-1,QW) + &
                     q(i-1,j,k+1,QW) - q(i-1,j,k-1,QW)) * dzinv
#endif

                cmin = min(qaux(i,j,k,QC), qaux(i-1,j,k,QC))

             else if (idir == 2) then

                ! normal direction
                avis(i,j,k) = (q(i,j,k,QV) - q(i,j-1,k,QV)) * dyinv

                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i+1,j,k,QU) - q(i-1,j,k,QU) + &
                     q(i+1,j-1,k,QU) - q(i-1,j-1,k,QU)) * dxinv

#if BL_SPACEDIM >= 3
                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i,j,k+1,QW) - q(i,j,k-1,QW) + &
                     q(i,j-1,k+1,QW) - q(i,j-1,k-1,QW)) * dzinv
#endif

                cmin = min(qaux(i,j,k,QC), qaux(i,j-1,k,QC))

             else

                ! normal direction
                avis(i,j,k) = (q(i,j,k,QW) - q(i,j,k-1,QW)) * dzinv

                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i,j+1,k,QV) - q(i,j-1,k,QV) + &
                     q(i,j+1,k-1,QV) - q(i,j-1,k-1,QV)) * dyinv

                avis(i,j,k) = avis(i,j,k) + 0.25_rt*( &
                     q(i+1,j,k,QU) - q(i-1,j,k,QU) + &
                     q(i+1,j,k-1,QU) - q(i-1,j,k-1,QU)) * dxinv

                cmin = min(qaux(i,j,k,QC), qaux(i,j,k-1,QC))

             endif

             ! MC Eq. 36
             coeff = min(ONE, (dx(idir)*avis(i,j,k))**2/(beta * cmin**2))

             if (avis(i,j,k) < ZERO) then
                avis(i,j,k) = dx(idir)*avis(i,j,k)*coeff
             else
                avis(i,j,k) = ZERO
             endif

          enddo
       enddo
    enddo

  end subroutine avisc


  subroutine normalize_species_fluxes(lo, hi, flux, f_lo, f_hi) bind(c, name="normalize_species_fluxes")
    ! Normalize the fluxes of the mass fractions so that
    ! they sum to 0.  This is essentially the CMA procedure that is
    ! defined in Plewa & Muller, 1999, A&A, 342, 179.
    !

    use network, only: nspec
    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO, UFS

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    ! Local variables
    integer  :: i, j, k, n
    real(rt) :: sum, fac

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             sum = ZERO

             do n = UFS, UFS+nspec-1
                sum = sum + flux(i,j,k,n)
             end do

             if (sum .ne. ZERO) then
                fac = flux(i,j,k,URHO) / sum
             else
                fac = ONE
             end if

             do n = UFS, UFS+nspec-1
                flux(i,j,k,n) = flux(i,j,k,n) * fac
             end do

          end do
       end do
    end do

  end subroutine normalize_species_fluxes



  subroutine apply_av(lo, hi, idir, dx, &
                      div, div_lo, div_hi, &
                      uin, uin_lo, uin_hi, &
                      flux, f_lo, f_hi) bind(c, name="apply_av")

    use amrex_constants_module, only: ZERO, FOURTH
    use meth_params_module, only: NVAR, UTEMP, USHK, difmag
    use prob_params_module, only: dg

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idir

    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    integer :: i, j, k, n

    real(rt) :: div1

    !$gpu

    do n = 1, NVAR

       if ( n == UTEMP ) cycle
#ifdef SHOCK_VAR
       if ( n == USHK  ) cycle
#endif

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (idir .eq. 1) then

                   div1 = FOURTH * (div(i,j,k        ) + div(i,j+1*dg(2),k        ) + &
                        div(i,j,k+1*dg(3)) + div(i,j+1*dg(2),k+1*dg(3)))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i-1*dg(1),j,k,n))

                else if (idir .eq. 2) then

                   div1 = FOURTH * (div(i,j,k        ) + div(i+1*dg(1),j,k        ) + &
                        div(i,j,k+1*dg(3)) + div(i+1*dg(1),j,k+1*dg(3)))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i,j-1*dg(2),k,n))

                else

                   div1 = FOURTH * (div(i,j        ,k) + div(i+1*dg(1),j        ,k) + &
                        div(i,j+1*dg(2),k) + div(i+1*dg(1),j+1*dg(2),k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i,j,k-1*dg(3),n))

                end if

                flux(i,j,k,n) = flux(i,j,k,n) + dx(idir) * div1

             end do
          end do
       end do

    end do

  end subroutine apply_av

#ifdef RADIATION
  subroutine apply_av_rad(lo, hi, idir, dx, &
                          div, div_lo, div_hi, &
                          Erin, Ein_lo, Ein_hi, &
                          radflux, rf_lo, rf_hi) bind(c, name="apply_av_rad")

    use amrex_constants_module, only: ZERO, FOURTH
    use meth_params_module, only: NVAR, UTEMP, USHK, difmag
    use prob_params_module, only: dg
    use rad_params_module, only: ngroups
    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: Ein_lo(3), Ein_hi(3)
    integer,  intent(in   ) :: rf_lo(3), rf_hi(3)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idir

    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: Erin(Ein_lo(1):Ein_hi(1),Ein_lo(2):Ein_hi(2),Ein_lo(3):Ein_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: radflux(rf_lo(1):rf_hi(1),rf_lo(2):rf_hi(2),rf_lo(3):rf_hi(3),0:ngroups-1)

    integer :: i, j, k, n

    real(rt) :: div1

    !$gpu

    do n = 0, ngroups-1

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (idir .eq. 1) then

                   div1 = FOURTH * (div(i,j,k        ) + div(i,j+1*dg(2),k        ) + &
                        div(i,j,k+1*dg(3)) + div(i,j+1*dg(2),k+1*dg(3)))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (Erin(i,j,k,n) - Erin(i-1*dg(1),j,k,n))

                else if (idir .eq. 2) then

                   div1 = FOURTH * (div(i,j,k        ) + div(i+1*dg(1),j,k        ) + &
                        div(i,j,k+1*dg(3)) + div(i+1*dg(1),j,k+1*dg(3)))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (Erin(i,j,k,n) - Erin(i,j-1*dg(2),k,n))

                else

                   div1 = FOURTH * (div(i,j        ,k) + div(i+1*dg(1),j        ,k) + &
                        div(i,j+1*dg(2),k) + div(i+1*dg(1),j+1*dg(2),k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (Erin(i,j,k,n) - Erin(i,j,k-1*dg(3),n))

                end if

                radflux(i,j,k,n) = radflux(i,j,k,n) + dx(idir) * div1

             end do
          end do
       end do

    end do

  end subroutine apply_av_rad
#endif


  subroutine scale_flux(lo, hi, &
#if AMREX_SPACEDIM == 1
                        qint, qi_lo, qi_hi, &
#endif
                        flux, f_lo, f_hi, &
                        area, a_lo, a_hi, dt) bind(c, name="scale_flux")

    use meth_params_module, only: NVAR, GDPRES, UMX, NGDNV
#if AMREX_SPACEDIM == 1
    use prob_params_module, only: coord_type
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
#if AMREX_SPACEDIM == 1
    integer,  intent(in   ) :: qi_lo(3), qi_hi(3)
#endif
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    real(rt), intent(in   ), value :: dt

    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
#if AMREX_SPACEDIM == 1
    real(rt), intent(in   ) :: qint(qi_lo(1):qi_hi(1), qi_lo(2):qi_hi(2), qi_lo(3):qi_hi(3), NGDNV)
#endif
    integer :: i, j, k, n

    !$gpu

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux(i,j,k,n) = dt * flux(i,j,k,n) * area(i,j,k)
#if AMREX_SPACEDIM == 1
                ! Correct the momentum flux with the grad p part.
                if (coord_type == 0 .and. n == UMX) then
                   flux(i,j,k,n) = flux(i,j,k,n) + dt * area(i,j,k) * qint(i,j,k,GDPRES)
                endif
#endif
             enddo
          enddo
       enddo
    enddo

  end subroutine scale_flux

#ifdef RADIATION
  subroutine scale_rad_flux(lo, hi, &
                            rflux, rf_lo, rf_hi, &
                            area, a_lo, a_hi, dt) bind(c, name="scale_rad_flux")

    use rad_params_module, only: ngroups

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: rf_lo(3), rf_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    real(rt), intent(in   ), value :: dt

    real(rt), intent(inout) :: rflux(rf_lo(1):rf_hi(1),rf_lo(2):rf_hi(2),rf_lo(3):rf_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    integer :: i, j, k, g

    !$gpu

    do g = 0, ngroups-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rflux(i,j,k,g) = dt * rflux(i,j,k,g) * area(i,j,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine scale_rad_flux
#endif

end module advection_util_module
