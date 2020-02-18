#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

// add the transverse flux difference in direction idir_t to the
// interface states in direction idir_n

void
Castro::trans_single(const Box& bx,
                     const int idir_t, const int idir_n,
                     Array4<Real const> const qm,
                     Array4<Real> const qmo,
                     Array4<Real const> const qp,
                     Array4<Real> const qpo,
                     Array4<Real const> const qaux,
                     Array4<Real const> const flux_t,
#ifdef RADIATION
                     Array4<Real const> const rflux_t,
#endif
                     Array4<Real const> const q_t,
#if AMREX_SPACEDIM == 2
                     Array4<Real const> const area_t,
                     Array4<Real const> const vol,
#endif
                     const Real hdt, const Real cdtdx)
{


  //       qm|qp
  //         |
  // --------+--------
  //   i-1       i
  //        i-1/2
  //
  // the qp state will see the transverse flux in zone i
  // the qm state will see the transverse flux in zone i-1

  // we account for this with the 'd' variable loop below
  // d = 0 will do qp and d = -1 will do qm

  // idir_t is the transverse direction and we set il,jl,kl
  // and ir,jr,kr to be the face-centered indices needed for
  // the transverse flux difference

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    for (int d = -1; d <= 0; d++) {

      // We are handling the states at the interface of
      // (i, i+1) in the x-direction, and similarly for
      // the y- and z- directions.

      int il = i;
      int jl = j;
      int kl = k;

      // set the face indices in the transverse direction

      if (idir_t == 1) {
        int ir = i+1;
        int jr = j;
        int kr = k;

      } else if (idir_t == 2) {
        int ir = i;
        int jr = j+1;
        int kr = k;

      } else {
        int ir = i;
        int jr = j;
        int kr = k+1;
      }

      // We're handling both the plus and minus states;
      // for the minus state we're shifting one zone to
      // the left in our chosen direction.

      if (idir_n == 1) {
        il += d;
        ir += d;

      } else if (idir_n == 2) {
        jl += d;
        jr += d;

      } else {
        kl += d;
        kr += d;
      }

      // store a local copy of the current interface state
      // this is qp or qm depending on the "d" loop

      GpuArray<Real, NQ> lqn;
      GpuArray<Real, NQ> lqno;

      if (d == -1) {
        for (int n = 0; n < NQ; n++) {
          lqn(n) = qm(i,j,k,n);
        }
      } else {
        for (int n = 0; n < NQ; n++) {
          lqn(n) = qp(i,j,k,n);
        }
      }

      // update all of the passively-advected quantities with the
      // transverse term and convert back to the primitive quantity

#if AMREX_SPACEDIM == 2
      const Real volinv = 1.0_rt/vol(il,jl,kl);
#endif

      for (int ipassive = 0; ipassive < npassive; ipassive++) {
        int n  = upass_map[ipassive];
        int nqp = qpass_map[ipassive];

#if AMREX_SPACEDIM == 2
        Real rrnew = lqn(QRHO) - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,URHO) -
                                      area_t(il,jl,kl)*flux_t(il,jl,kl,URHO)) * volinv;
        Real compu = lqn(QRHO)*lqn(nqp) -
          hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,n) -
               area_t(il,jl,kl)*flux_t(il,jl,kl,n)) * volinv;
        lqno(nqp) = compu/rrnew;
#else
        Real rrnew = lqn(QRHO) - cdtdx*(flux_t(ir,jr,kr,URHO) - flux_t(il,jl,kl,URHO));
        Real compu = lqn(QRHO)*lqn(nqp) - cdtdx*(flux_t(ir,jr,kr,n) - flux_t(il,jl,kl,n));
        Real lqno(nqp) = compu/rrnew;
#endif
      }

      Real pgp  = q_t(ir,jr,kr,GDPRES);
      Real pgm  = q_t(il,jl,kl,GDPRES);
      Real ugp  = q_t(ir,jr,kr,GDU+idir_t);
      Real ugm  = q_t(il,jl,kl,GDU+idir_t);
      Real gegp = q_t(ir,jr,kr,GDGAME);
      Real gegm = q_t(il,jl,kl,GDGAME);

#ifdef RADIATION
      GpuArray<Real, ngroups> lambda;
      GpuArray<Real, ngroups> ergp;
      GpuArray<Real, ngroups> ergm;

      for (int g = 0; g < ngroups; g++) {
        lambda(g) = qaux(il,jl,kl,QLAMS+g);
        ergp(g) = q_t(ir,jr,kr,GDERADS+g);
        ergm(g) = q_t(il,jl,kl,GDERADS+g);
      }
#endif

      // we need to augment our conserved system with either a p
      // equation or gammae (if we have ppm_predict_gammae = 1) to
      // be able to deal with the general EOS

#if AMREX_SPACEDIM == 2
      Real dup = area_t(ir,jr,kr)*pgp*ugp - area_t(il,jl,kl)*pgm*ugm;
      Real du = area_t(ir,jr,kr)*ugp-area_t(il,jl,kl)*ugm;
#else
      Real dup = pgp*ugp - pgm*ugm;
      Real du = ugp-ugm;
#endif
      Real pav = HALF*(pgp+pgm);
      Real uav = HALF*(ugp+ugm);
      Real geav = HALF*(gegp+gegm);
      Real dge = gegp-gegm;

      // this is the gas gamma_1
#ifdef RADIATION
      Real gamc = qaux(il,jl,kl,QGAMCG);
#else
      Real gamc = qaux(il,jl,kl,QGAMC);
#endif

#ifdef RADIATION
      GpuArray<Real, ngroups> lamge;
      GpuArray<Real, ngroups> luge;
      GpuArray<Real, ngroups> der;
      Real dmom = 0.0_rt;
      Real dre = 0.0_rt;
      for (int g = 0; g < ngroups; g++) {
        lamge(g) = lambda(g) * (ergp(g)-ergm(g));
        dmom += - cdtdx*lamge(g);
        luge(g) = uav * lamge(g);
        dre += -cdtdx*luge(g);
      }

      if (fspace_type == 1 && comoving) {
        for (int g=0; g < ngroups; g++) {
          Real eddf = Edd_factor(lambda(g));
          Real f1 = 0.5_rt*(1.0_rt - eddf);
          der(g) = cdtdx * uav * f1 * (ergp(g) - ergm(g));
        }

      } else if (fspace_type == 2) {
#if AMREX_SPACEDIM == 2
        Real divu = (area_t(ir,jr,kr)*ugp - area_t(il,jl,kl)*ugm) * volinv;
        for (int g=0; g < ngroups; g++) {
          Real eddf = Edd_factor(lambda(g));
          Real f1 = 0.5_rt*(1.0_rt - eddf);
          der(g) = -hdt * f1 * 0.5_rt*(ergp(g)+ergm(g)) * divu;
        }
#else
        for (int g=0, g < ngroups; g++) {
          Real eddf = Edd_factor(lambda(g));
          Real f1 = 0.5_rt*(1.0_rt-eddf);
          der(g) = cdtdx * f1 * 0.5_rt*(ergp(g)+ergm(g)) * (ugm-ugp);
        }
#endif
      } else { // mixed frame
        for (int g = 0; g < ngroups; g++) {
          der(g) = cdtdx * luge(g);
        }
      }
#endif

      // Convert to conservation form
      Real rrn = lqn(QRHO);
      Real run = rrn*lqn(QU);
      Real rvn = rrn*lqn(QV);
      Real rwn = rrn*lqn(QW);
      Real ekenn = 0.5_rt * rrn * (lqn(QU)*lqn(QU) + lqn(QV)*lqn(QV) + lqn(QW)*lqn(QW));
      Real ren = lqn(QREINT) + ekenn;
#ifdef RADIATION
      for (int g = 0; g < ngroups; g++) {
        ern(g) = lqn(QRAD+g);
      }
#endif

#if AMREX_SPACEDIM == 2
      // Add transverse predictor
      Real rrnewn = rrn - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,URHO) -
                               area_t(il,jl,kl)*flux_t(il,jl,kl,URHO)) * volinv;

      // Note that pressure may be treated specially here, depending on
      // the geometry.  Our y-interface equation for (rho u) is:
      //
      //  d(rho u)/dt + d(rho u v)/dy = - 1/r d(r rho u u)/dr - dp/dr
      //
      // in cylindrical coords -- note that the p term is not
      // in a divergence for UMX in the x-direction, so there
      // are no area factors.  For this geometry, we do not
      // include p in our definition of the flux in the
      // x-direction, for we need to fix this now.
      Real runewn = run - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,UMX) -
                               area_t(il,jl,kl)*flux_t(il,jl,kl,UMX)) * volinv;
      if (idir_t == 0 && ! momx_flux_has_p(idir_t)) {
        runewn = runewn - cdtdx * (pgp-pgm);
      }
      Real rvnewn = rvn - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,UMY) -
                               area_t(il,jl,kl)*flux_t(il,jl,kl,UMY)) * volinv;
      Real rwnewn = rwn - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,UMZ) -
                               area_t(il,jl,kl)*flux_t(il,jl,kl,UMZ)) * volinv;
      Real renewn = ren - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,UEDEN) -
                               area_t(il,jl,kl)*flux_t(il,jl,kl,UEDEN)) * volinv;

#ifdef RADIATION
      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      runewn = runewn - HALF*hdt*(area_t(ir,jr,kr)+area_t(il,jl,kl))*sum(lamge) * volinv;
      renewn = renewn + dre;
      for (int g = 0; g < ngroups; g++) {
        ernewn(g) = ern(g) - hdt*(area_t(ir,jr,kr)*rflux_t(ir,jr,kr,g) -
                                  area_t(il,jl,kl)*rflux_t(il,jl,kl,g)) * volinv + der(g);
      }
#endif

#else
      // Add transverse predictor
      rrnewn = rrn - cdtdx*(flux_t(ir,jr,kr,URHO) - flux_t(il,jl,kl,URHO));
      runewn = run - cdtdx*(flux_t(ir,jr,kr,UMX) - flux_t(il,jl,kl,UMX));
      rvnewn = rvn - cdtdx*(flux_t(ir,jr,kr,UMY) - flux_t(il,jl,kl,UMY));
      rwnewn = rwn - cdtdx*(flux_t(ir,jr,kr,UMZ) - flux_t(il,jl,kl,UMZ));
      renewn = ren - cdtdx*(flux_t(ir,jr,kr,UEDEN) - flux_t(il,jl,kl,UEDEN));
#ifdef RADIATION
      runewn = runewn + dmom;
      renewn = renewn + dre;
      for (int g = 0; g < ngroups; g++) {
        ernewn(g)  = ern(g) - cdtdx*(rflux_t(ir,jr,kr,g) - rflux_t(il,jl,kl,g)) + der(g);
      }
#endif
#endif

      // Reset to original value if adding transverse terms made density negative
      bool reset_state = false;
      if (transverse_reset_density == 1 && rrnewn < 0.0_rt) {
        rrnewn = rrn;
        runewn = run;
        rvnewn = rvn;
        rwnewn = rwn;
        renewn = ren;
#ifdef RADIATION
        for (int g = 0; g < ngroups; g++) {
          ernewn(g) = ern(g);
        }
#endif
        reset_state = true;
      }

      // Convert back to primitive form
      lqno(QRHO) = rrnewn;
      Real rhoinv = 1.0_rt/rrnewn;
      lqno(QU) = runewn*rhoinv;
      lqno(QV) = rvnewn*rhoinv;
      lqno(QW) = rwnewn*rhoinv;

      // note: we run the risk of (rho e) being negative here
      Real rhoekenn = HALF*(runewn**2 + rvnewn**2 + rwnewn**2)*rhoinv;
                lqno(QREINT) = renewn - rhoekenn

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. lqno(QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating (rho e).
#if AMREX_SPACEDIM == 2
                      lqno(QREINT) = lqn(QREINT) - &
                           hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,UEINT) - &
                                area_t(il,jl,kl)*flux_t(il,jl,kl,UEINT) + pav*du) * volinv
#else
                      lqno(QREINT) = lqn(QREINT) - &
                           cdtdx*(flux_t(ir,jr,kr,UEINT) - flux_t(il,jl,kl,UEINT) + pav*du)
#endif
                   end if

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
#if AMREX_SPACEDIM == 2
                      ! the divergences here, dup and du, already have area factors
                      pnewn = lqn(QPRES) - hdt*(dup + pav*du*(gamc - ONE)) * volinv
#else
                      pnewn = lqn(QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
#endif
                      lqno(QPRES) = max(pnewn, small_pres)
                   else
                      ! Update gammae with its transverse terms
#if AMREX_SPACEDIM == 2
                      lqno(QGAME) = lqn(QGAME) + &
                           hdt*( (geav-ONE)*(geav - gamc)*du) * volinv - cdtdx*uav*dge
#else
                      lqno(QGAME) = lqn(QGAME) + &
                           cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )
#endif
                      ! and compute the p edge state from this and (rho e)
                      lqno(QPRES) = lqno(QREINT)*(lqno(QGAME) - ONE)
                      lqno(QPRES) = max(lqno(QPRES), small_pres)
                   end if
                else
                   lqno(QPRES) = lqn(QPRES)
                   lqno(QGAME) = lqn(QGAME)
                endif

#ifdef RADIATION
                lqno(qrad:qrad-1+ngroups) = ernewn(:)
                lqno(qptot  ) = sum(lambda(:)*ernewn(:)) + lqno(QPRES)
                lqno(qreitot) = sum(lqno(qrad:qrad-1+ngroups)) + lqno(QREINT)
#endif

                ! store the state
                if (d == -1) then
                   qmo(i,j,k,:) = lqno(:)
                else
                   qpo(i,j,k,:) = lqno(:)
                end if

             end do   ! d loop

          end do
       end do
    end do

  end subroutine trans_single

  subroutine trans_final(lo, hi, &
                         idir_n, idir_t1, idir_t2, &
                         qm, qm_lo, qm_hi, &
                         qmo, qmo_lo, qmo_hi, &
                         qp, qp_lo, qp_hi, &
                         qpo, qpo_lo, qpo_hi, &
                         qaux, qa_lo, qa_hi, &
                         flux_t1, ft1_lo, ft1_hi, &
#ifdef RADIATION
                         rflux_t1, rft1_lo, rft1_hi, &
#endif
                         flux_t2, ft2_lo, ft2_hi, &
#ifdef RADIATION
                         rflux_t2, rft2_lo, rft2_hi, &
#endif
                         q_t1, qt1_lo, qt1_hi, &
                         q_t2, qt2_lo, qt2_hi, &
                         hdt, cdtdx_n, cdtdx_t1, cdtdx_t2) bind(C, name="trans_final")

    ! here, lo and hi are the bounds of the x interfaces we are looping over

    use amrex_constants_module, only : ZERO, ONE, HALF

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qptot, qreitot, &
                                   fspace_type, comoving, &
                                   GDERADS, GDLAMS, &
                                   QCG, QGAMCG, QLAMS, &
#endif
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                   small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_predict_gammae

#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor ! function
#endif

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qmo_lo(3), qmo_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qpo_lo(3), qpo_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: ft1_lo(3), ft1_hi(3)
    integer, intent(in) :: ft2_lo(3), ft2_hi(3)
    integer, intent(in) :: qt1_lo(3), qt1_hi(3)
    integer, intent(in) :: qt2_lo(3), qt2_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    integer, intent(in), value :: idir_n, idir_t1, idir_t2

    real(rt), intent(in), value :: hdt, cdtdx_n, cdtdx_t1, cdtdx_t2

#ifdef RADIATION
    integer, intent(in) :: rft1_lo(3), rft1_hi(3)
    integer, intent(in) :: rft2_lo(3), rft2_hi(3)
    real(rt), intent(in) :: rflux_t1(rft1_lo(1):rft1_hi(1),rft1_lo(2):rft1_hi(2),rft1_lo(3):rft1_hi(3),0:ngroups-1)
    real(rt), intent(in) :: rflux_t2(rft2_lo(1):rft2_hi(1),rft2_lo(2):rft2_hi(2),rft2_lo(3):rft2_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(in) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)
    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),NQ)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),NQ)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: flux_t1(ft1_lo(1):ft1_hi(1),ft1_lo(2):ft1_hi(2),ft1_lo(3):ft1_hi(3),NVAR)
    real(rt), intent(in) :: flux_t2(ft2_lo(1):ft2_hi(1),ft2_lo(2):ft2_hi(2),ft2_lo(3):ft2_hi(3),NVAR)
    real(rt), intent(in) :: q_t1(qt1_lo(1):qt1_hi(1),qt1_lo(2):qt1_hi(2),qt1_lo(3):qt1_hi(3),NGDNV)
    real(rt), intent(in) :: q_t2(qt2_lo(1):qt2_hi(1),qt2_lo(2):qt2_hi(2),qt2_lo(3):qt2_hi(3),NGDNV)

    integer i, j, k, n, nqp, ipassive

    real(rt)         rrn, run, rvn, rwn, ren, ekenn, rhoekenn
    real(rt)         rrnewn, runewn, rvnewn, rwnewn, renewn
    real(rt)         pnewn
    real(rt)         pgt1p, ugt1p, gegt1p, pgt1m, ugt1m, gegt1m, dut1, dupt1, pt1av, pt1new, get1new
    real(rt)         pgt2p, ugt2p, gegt2p, pgt2m, ugt2m, gegt2m, dut2, dupt2, pt2av, pt2new, get2new
    real(rt)         ut1av, get1av, dget1, ut2av, get2av, dget2
    real(rt)         compn, compnn

#ifdef RADIATION
    real(rt) :: dmt1, dmt2, dre
    real(rt), dimension(0:ngroups-1) :: der, lambda, luget1, luget2, lget1, lget2, &
         ern, ernewn, ergt1m, ergt1p, ergt2m, ergt2p
    real(rt) :: eddf, f1
    integer :: g
#endif

    real(rt) :: lqn(NQ), lqno(NQ)

    logical :: reset_state

    integer :: d
    integer :: iln, jln, kln
    integer :: il_t1, jl_t1, kl_t1
    integer :: ir_t1, jr_t1, kr_t1
    integer :: il_t2, jl_t2, kl_t2
    integer :: ir_t2, jr_t2, kr_t2

    !$gpu

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             do d = -1, 0

                ! the normal state
                iln = i
                jln = j
                kln = k

                ! the first transverse state
                il_t1 = i
                jl_t1 = j
                kl_t1 = k

                ! the second transverse state
                il_t2 = i
                jl_t2 = j
                kl_t2 = k

                if (idir_n == 1) then
                   ! x is the normal direction

                   ! y is the first transverse state
                   ir_t1 = i+d
                   jr_t1 = j+1
                   kr_t1 = k

                   ! z is the second transverse state
                   ir_t2 = i+d
                   jr_t2 = j
                   kr_t2 = k+1

                   ! offset for the plus/minus state
                   iln = i+d
                   il_t1 = i+d
                   il_t2 = i+d

                else if (idir_n == 2) then
                   ! y is the normal direction

                   ! x is the first transverse state
                   ir_t1 = i+1
                   jr_t1 = j+d
                   kr_t1 = k

                   ! z is the second transverse state
                   ir_t2 = i
                   jr_t2 = j+d
                   kr_t2 = k+1

                   ! offset for the plus/minus state
                   jln = j+d
                   jl_t1 = j+d
                   jl_t2 = j+d

                else
                   ! z is the normal direction

                   ! x is the first transverse state
                   ir_t1 = i+1
                   jr_t1 = j
                   kr_t1 = k+d

                   ! y is the second transverse state
                   ir_t2 = i
                   jr_t2 = j+1
                   kr_t2 = k+d

                   ! offset for the plus/minus state
                   kln = k+d
                   kl_t1 = k+d
                   kl_t2 = k+d
                end if

                if (d == -1) then
                   lqn(:) = qm(i,j,k,:)
                else
                   lqn(:) = qp(i,j,k,:)
                end if

                ! update all of the passively-advected quantities with the
                ! transerse term and convert back to the primitive quantity

                do ipassive = 1,npassive
                   n  = upass_map(ipassive)
                   nqp = qpass_map(ipassive)

                   rrn = lqn(QRHO)
                   compn = rrn*lqn(nqp)
                   rrnewn = rrn - cdtdx_t1*(flux_t1(ir_t1,jr_t1,kr_t1,URHO) - &
                                            flux_t1(il_t1,jl_t1,kl_t1,URHO)) &
                                - cdtdx_t2*(flux_t2(ir_t2,jr_t2,kr_t2,URHO) - &
                                            flux_t2(il_t2,jl_t2,kl_t2,URHO))
                   compnn = compn - cdtdx_t1*(flux_t1(ir_t1,jr_t1,kr_t1,n) - &
                                              flux_t1(il_t1,jl_t1,kl_t1,n)) &
                                  - cdtdx_t2*(flux_t2(ir_t2,jr_t2,kr_t2,n) - &
                                              flux_t2(il_t2,jl_t2,kl_t2,n))

                   lqno(nqp) = compnn/rrnewn
                end do

                ! add the transverse differences to the normal states for the
                ! fluid variables

                pgt1p  = q_t1(ir_t1,jr_t1,kr_t1,GDPRES)
                pgt1m  = q_t1(il_t1,jl_t1,kl_t1,GDPRES)
                ugt1p  = q_t1(ir_t1,jr_t1,kr_t1,GDU+idir_t1-1)
                ugt1m  = q_t1(il_t1,jl_t1,kl_t1,GDU+idir_t1-1)
                gegt1p = q_t1(ir_t1,jr_t1,kr_t1,GDGAME)
                gegt1m = q_t1(il_t1,jl_t1,kl_t1,GDGAME)
#ifdef RADIATION
                ergt1p = q_t1(ir_t1,jr_t1,kr_t1,GDERADS:GDERADS-1+ngroups)
                ergt1m = q_t1(il_t1,jl_t1,kl_t1,GDERADS:GDERADS-1+ngroups)
#endif

                pgt2p  = q_t2(ir_t2,jr_t2,kr_t2,GDPRES)
                pgt2m  = q_t2(il_t2,jl_t2,kl_t2,GDPRES)
                ugt2p  = q_t2(ir_t2,jr_t2,kr_t2,GDU+idir_t2-1)
                ugt2m  = q_t2(il_t2,jl_t2,kl_t2,GDU+idir_t2-1)
                gegt2p = q_t2(ir_t2,jr_t2,kr_t2,GDGAME)
                gegt2m = q_t2(il_t2,jl_t2,kl_t2,GDGAME)
#ifdef RADIATION
                ergt2p = q_t2(ir_t2,jr_t2,kr_t2,GDERADS:GDERADS-1+ngroups)
                ergt2m = q_t2(il_t2,jl_t2,kl_t2,GDERADS:GDERADS-1+ngroups)
#endif

                dupt1 = pgt1p*ugt1p - pgt1m*ugt1m
                pt1av = HALF*(pgt1p + pgt1m)
                ut1av = HALF*(ugt1p + ugt1m)
                get1av = HALF*(gegt1p + gegt1m)
                dut1 = ugt1p - ugt1m
                dget1 = gegt1p - gegt1m
#ifdef RADIATION
                pt1new = cdtdx_t1*(dupt1 + pt1av*dut1*(qaux(iln,jln,kln,QGAMCG) - ONE))
                get1new = cdtdx_t1*( (get1av-ONE)*(get1av - qaux(iln,jln,kln,QGAMCG))*dut1 - ut1av*dget1 )
#else
                pt1new = cdtdx_t1*(dupt1 + pt1av*dut1*(qaux(iln,jln,kln,QGAMC) - ONE))
                get1new = cdtdx_t1*( (get1av-ONE)*(get1av - qaux(iln,jln,kln,QGAMC))*dut1 - ut1av*dget1 )
#endif

                dupt2 = pgt2p*ugt2p - pgt2m*ugt2m
                pt2av = HALF*(pgt2p + pgt2m)
                ut2av = HALF*(ugt2p + ugt2m)
                get2av = HALF*(gegt2p + gegt2m)
                dut2 = ugt2p - ugt2m
                dget2 = gegt2p - gegt2m
#ifdef RADIATION
                pt2new = cdtdx_t2*(dupt2 + pt2av*dut2*(qaux(iln,jln,kln,QGAMCG) - ONE))
                get2new = cdtdx_t2*( (get2av-ONE)*(get2av - qaux(iln,jln,kln,QGAMCG))*dut2 - ut2av*dget2 )
#else
                pt2new = cdtdx_t2*(dupt2 + pt2av*dut2*(qaux(iln,jln,kln,QGAMC) - ONE))
                get2new = cdtdx_t2*( (get2av-ONE)*(get2av - qaux(iln,jln,kln,QGAMC))*dut2 - ut2av*dget2 )
#endif

#ifdef RADIATION
                lambda(:) = qaux(iln,jln,kln,QLAMS:QLAMS+ngroups-1)

                lget1 = lambda(:) * (ergt1p(:) - ergt1m(:))
                lget2 = lambda(:) * (ergt2p(:) - ergt2m(:))
                dmt1 = - cdtdx_t1*sum(lget1)
                dmt2 = - cdtdx_t2*sum(lget2)
                luget1 = HALF*(ugt1p + ugt1m) * lget1(:)
                luget2 = HALF*(ugt2p + ugt2m) * lget2(:)
                dre = -cdtdx_t1*sum(luget1) - cdtdx_t2*sum(luget2)

                if (fspace_type .eq. 1 .and. comoving) then
                   do g=0, ngroups-1
                      eddf = Edd_factor(lambda(g))
                      f1 = HALF*(ONE-eddf)
                      der(g) = f1*(cdtdx_t1*HALF*(ugt1p + ugt1m)*(ergt1p(g) - ergt1m(g)) + &
                                   cdtdx_t2*HALF*(ugt2p + ugt2m)*(ergt2p(g) - ergt2m(g)) )
                   end do
                else if (fspace_type .eq. 2) then
                   do g=0, ngroups-1
                      eddf = Edd_factor(lambda(g))
                      f1 = HALF*(ONE-eddf)
                      der(g) = f1*(cdtdx_t1*HALF*(ergt1p(g) + ergt1m(g))*(ugt1m - ugt1p) + &
                                   cdtdx_t2*HALF*(ergt2p(g) + ergt2m(g))*(ugt2m - ugt2p) )
                   end do
                else ! mixed frame
                   der(:) = cdtdx_t1*luget1 + cdtdx_t2*luget2
                end if
#endif

                ! Convert to conservation form
                rrn = lqn(QRHO)
                run = rrn*lqn(QU)
                rvn = rrn*lqn(QV)
                rwn = rrn*lqn(QW)
                ekenn = HALF*rrn*sum(lqn(QU:QW)**2)
                ren = lqn(QREINT) + ekenn
#ifdef RADIATION
                ern(:) = lqn(qrad:qrad-1+ngroups)
#endif

                ! Add transverse predictor
                rrnewn = rrn - cdtdx_t1*(flux_t1(ir_t1,jr_t1,kr_t1,URHO) - &
                                         flux_t1(il_t1,jl_t1,kl_t1,URHO)) &
                             - cdtdx_t2*(flux_t2(ir_t2,jr_t2,kr_t2,URHO) - &
                                         flux_t2(il_t2,jl_t2,kl_t2,URHO))
                runewn = run - cdtdx_t1*(flux_t1(ir_t1,jr_t1,kr_t1,UMX) - &
                                         flux_t1(il_t1,jl_t1,kl_t1,UMX)) &
                             - cdtdx_t2*(flux_t2(ir_t2,jr_t2,kr_t2,UMX) - &
                                         flux_t2(il_t2,jl_t2,kl_t2,UMX))
                rvnewn = rvn - cdtdx_t1*(flux_t1(ir_t1,jr_t1,kr_t1,UMY) - &
                                         flux_t1(il_t1,jl_t1,kl_t1,UMY)) &
                             - cdtdx_t2*(flux_t2(ir_t2,jr_t2,kr_t2,UMY) - &
                                         flux_t2(il_t2,jl_t2,kl_t2,UMY))
                rwnewn = rwn - cdtdx_t1*(flux_t1(ir_t1,jr_t1,kr_t1,UMZ) - &
                                         flux_t1(il_t1,jl_t1,kl_t1,UMZ)) &
                             - cdtdx_t2*(flux_t2(ir_t2,jr_t2,kr_t2,UMZ) - &
                                         flux_t2(il_t2,jl_t2,kl_t2,UMZ))
                renewn = ren - cdtdx_t1*(flux_t1(ir_t1,jr_t1,kr_t1,UEDEN) - &
                                         flux_t1(il_t1,jl_t1,kl_t1,UEDEN)) &
                             - cdtdx_t2*(flux_t2(ir_t2,jr_t2,kr_t2,UEDEN) - &
                                         flux_t2(il_t2,jl_t2,kl_t2,UEDEN))
#ifdef RADIATION
                if (idir_n == 1) then
                   rvnewn = rvnewn + dmt1
                   rwnewn = rwnewn + dmt2
                else if (idir_n == 2) then
                   runewn = runewn + dmt1
                   rwnewn = rwnewn + dmt2
                else
                   runewn = runewn + dmt1
                   rvnewn = rvnewn + dmt2
                end if
                renewn = renewn + dre
                ernewn = ern(:) - cdtdx_t1*(rflux_t1(ir_t1,jr_t1,kr_t1,:) - &
                                            rflux_t1(il_t1,jl_t1,kl_t1,:)) &
                                - cdtdx_t2*(rflux_t2(ir_t2,jr_t2,kr_t2,:) - &
                                            rflux_t2(il_t2,jl_t2,kl_t2,:)) &
                                + der(:)
#endif

                ! Reset to original value if adding transverse terms
                ! made density negative
                reset_state = .false.
                if (transverse_reset_density == 1 .and. rrnewn < ZERO) then
                   rrnewn = rrn
                   runewn = run
                   rvnewn = rvn
                   rwnewn = rwn
                   renewn = ren
#ifdef RADIATION
                   ernewn = ern(:)
#endif
                   reset_state = .true.
                end if

                lqno(QRHO) = rrnewn
                lqno(QU) = runewn/rrnewn
                lqno(QV) = rvnewn/rrnewn
                lqno(QW) = rwnewn/rrnewn

                ! note: we run the risk of (rho e) being negative here
                rhoekenn = HALF*(runewn**2 + rvnewn**2 + rwnewn**2)/rrnewn
                lqno(QREINT) = renewn - rhoekenn

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. lqno(QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating
                      ! (rho e).
                      lqno(QREINT) = lqn(QREINT) &
                           - cdtdx_t1*(flux_t1(ir_t1,jr_t1,kr_t1,UEINT) - &
                                       flux_t1(il_t1,jl_t1,kl_t1,UEINT) + pt1av*dut1) &
                           - cdtdx_t2*(flux_t2(ir_t2,jr_t2,kr_t2,UEINT) - &
                                       flux_t2(il_t2,jl_t2,kl_t2,UEINT) + pt2av*dut2)
                   end if

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewn = lqn(QPRES) - pt1new - pt2new
                      lqno(QPRES) = pnewn
                   else
                      ! Update gammae with its transverse terms
                      lqno(QGAME) = lqn(QGAME) + get1new + get2new

                      ! and compute the p edge state from this and (rho e)
                      lqno(QPRES) = lqno(QREINT)*(lqno(QGAME)-ONE)
                   end if
                else
                   lqno(QPRES) = lqn(QPRES)
                   lqno(QGAME) = lqn(QGAME)
                endif

                lqno(QPRES) = max(lqno(QPRES), small_pres)

#ifdef RADIATION
                lqno(qrad:qrad-1+ngroups) = ernewn(:)
                lqno(qptot  ) = sum(lambda(:)*ernewn(:)) + lqno(QPRES)
                lqno(qreitot) = sum(lqno(qrad:qrad-1+ngroups)) + lqno(QREINT)
#endif

                if (d == -1) then
                   qmo(i,j,k,:) = lqno(:)
                else
                   qpo(i,j,k,:) = lqno(:)
                end if

             end do

          end do
       end do
    end do

  end subroutine trans_final


end module transverse_module
