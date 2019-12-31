
module rad_nd_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), parameter, private :: tiny = 1.e-50_rt
  real(rt), parameter, private :: BIGKR = 1.e25_rt

contains

  subroutine flxlim(lo, hi, &
                    lambda, l_lo, l_hi, &
                    limiter) &
                    bind(C, name="flxlim")

    use amrex_fort_module, only: rt => amrex_real
    use rad_util_module, only: FLDlambda ! function

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: l_lo(3), l_hi(3)
    real(rt), intent(inout) :: lambda(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3))
    integer,  intent(in   ), value :: limiter

    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             lambda(i,j,k) = FLDlambda(lambda(i,j,k), limiter)

          end do
       end do
    end do

  end subroutine flxlim



  subroutine scgrd(lo, hi, &
                   r, r_lo, r_hi, &
                   idir, dx, &
                   kappar, k_lo, k_hi, &
                   er, e_lo, e_hi, &
                   include_cross_terms) &
                   bind(C, name="scgrd")

    use amrex_constants_module, only: FOURTH, HALF
    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    integer,  intent(in   ) :: e_lo(3), e_hi(3)
    real(rt), intent(inout) :: r(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: kappar(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3))
    real(rt), intent(in   ) :: er(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idir, include_cross_terms

    integer  :: i, j, k, d
    real(rt) :: kap
    real(rt) :: rg, dal, dar, dbl, dbr
    real(rt) :: dxInv(3)

    !$gpu

    do d = 1, dim
       dxInv(d) = 1.e0_rt / dx(d)
    end do
    do d = dim+1, 3
       dxInv(d) = 0.e0_rt
    end do

    if (idir == 0) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (include_cross_terms == 1) then

#if (BL_SPACEDIM >= 2)
                   dal = er(i-1,j+1,k) - er(i-1,j-1,k)
                   dar = er(i  ,j+1,k) - er(i  ,j-1,k)

                   if      (er(i-1,j-1,k) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i-1,j+1,k) - er(i-1,j  ,k))
                   else if (er(i-1,j+1,k) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i-1,j  ,k) - er(i-1,j-1,k))
                   end if

                   if      (er(i  ,j-1,k) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i  ,j+1,k) - er(i  ,j  ,k))
                   else if (er(i  ,j+1,k) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i  ,j  ,k) - er(i  ,j-1,k))
                   end if
#else
                   dal = 0.e0_rt
                   dar = 0.e0_rt
#endif

#if (BL_SPACEDIM == 3)
                   dbl = er(i-1,j,k+1) - er(i-1,j,k-1)
                   dbr = er(i  ,j,k+1) - er(i  ,j,k-1)

                   if      (er(i-1,j,k-1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i-1,j,k+1) - er(i-1,j,k  ))
                   else if (er(i-1,j,k+1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i-1,j,k  ) - er(i-1,j,k-1))
                   end if

                   if      (er(i  ,j,k-1) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j,k+1) - er(i  ,j,k  ))
                   else if (er(i  ,j,k+1) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j,k  ) - er(i  ,j,k-1))
                   end if
#else
                   dbl = 0.e0_rt
                   dbr = 0.e0_rt
#endif

                else

                   dal = 0.e0_rt
                   dar = 0.e0_rt
                   dbl = 0.e0_rt
                   dbr = 0.e0_rt

                end if
                   

                if (er(i-1,j,k) == -1.e0_rt) then

                   rg = ((er(i+1,j,k) - er(i,j,k)) * dxInv(1))**2 + &
                         (HALF * dar * dxInv(2))**2 + (HALF * dbr * dxInv(3))**2

                else if (er(i,j,k) == -1.e0_rt) then

                   rg = ((er(i-1,j,k) - er(i-2,j,k)) * dxInv(1))**2 + &
                        (HALF * dal * dxInv(2))**2 + (HALF * dbl * dxInv(3))**2

                else

                   rg = ((er(i,j,k) - er(i-1,j,k)) * dxInv(1))**2 + &
                        (FOURTH * (dal + dar) * dxInv(2))**2 + (FOURTH * (dbl + dbr) * dxInv(3))**2

                end if

                kap = kavg(kappar(i-1,j,k), kappar(i,j,k), dx(1), -1)
                r(i,j,k) = sqrt(rg) / (kap * max(er(i-1,j,k), er(i,j,k), tiny))

             end do
          end do
       end do

    else if (idir == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (include_cross_terms == 1) then

                   dal = er(i+1,j-1,k  ) - er(i-1,j-1,k  )
                   dar = er(i+1,j  ,k  ) - er(i-1,j  ,k  )

                   if      (er(i-1,j-1,k  ) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i+1,j-1,k  ) - er(i  ,j-1,k  ))
                   else if (er(i+1,j-1,k  ) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i  ,j-1,k  ) - er(i-1,j-1,k  ))
                   end if

                   if      (er(i-1,j  ,k  ) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i+1,j  ,k  ) - er(i  ,j  ,k  ))
                   else if (er(i+1,j  ,k  ) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i  ,j  ,k  ) - er(i-1,j  ,k  ))
                   end if

#if (BL_SPACEDIM == 3)
                   dbl = er(i  ,j-1,k+1) - er(i  ,j-1,k-1)
                   dbr = er(i  ,j  ,k+1) - er(i  ,j  ,k-1)

                   if      (er(i  ,j-1,k-1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i  ,j-1,k+1) - er(i  ,j-1,k  ))
                   else if (er(i  ,j-1,k+1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i  ,j-1,k  ) - er(i  ,j-1,k-1))
                   end if

                   if      (er(i  ,j  ,k-1) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j  ,k+1) - er(i  ,j,  k  ))
                   else if (er(i  ,j  ,k+1) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j  ,k  ) - er(i  ,j,  k-1))
                   end if
#else
                   dbl = 0.e0_rt
                   dbr = 0.e0_rt
#endif

                else

                   dal = 0.e0_rt
                   dar = 0.e0_rt
                   dbl = 0.e0_rt
                   dbr = 0.e0_rt

                end if
                   

                if (er(i,j-1,k) == -1.e0_rt) then

                   rg = ((er(i,j+1,k) - er(i,j,k)) * dxInv(2))**2 + &
                        (HALF * dar * dxInv(1))**2 + (HALF * dbr * dxInv(3))**2

                else if (er(i,j,k) == -1.e0_rt) then

                   rg = ((er(i,j-1,k) - er(i,j-2,k)) * dxInv(2))**2 + &
                        (HALF * dal * dxInv(1))**2 + (HALF * dbl * dxInv(3))**2

                else

                   rg = ((er(i,j,k) - er(i,j-1,k)) * dxInv(2))**2 + &
                        (FOURTH * (dal + dar) * dxInv(1))**2 + (FOURTH * (dbl + dbr) * dxInv(3))**2

                end if

                kap = kavg(kappar(i,j-1,k), kappar(i,j,k), dx(2), -1)
                r(i,j,k) = sqrt(rg) / (kap * max(er(i,j-1,k), er(i,j,k), tiny))

             end do
          end do
       end do

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (include_cross_terms == 1) then

                   dal = er(i+1,j  ,k-1) - er(i-1,j  ,k-1)
                   dar = er(i+1,j  ,k  ) - er(i-1,j  ,k  )

                   if      (er(i-1,j  ,k-1) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i+1,j  ,k-1) - er(i  ,j  ,k-1))
                   else if (er(i+1,j  ,k-1) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i  ,j  ,k-1) - er(i-1,j  ,k-1))
                   end if

                   if      (er(i-1,j  ,k  ) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i+1,j  ,k  ) - er(i  ,j  ,k  ))
                   else if (er(i+1,j  ,k  ) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i  ,j  ,k  ) - er(i-1,j  ,k  ))
                   end if

                   dbl = er(i  ,j+1,k-1) - er(i  ,j-1,k-1)
                   dbr = er(i  ,j+1,k  ) - er(i  ,j-1,k  )

                   if      (er(i  ,j-1,k-1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i  ,j+1,k-1) - er(i  ,j  ,k-1))
                   else if (er(i  ,j+1,k-1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i  ,j  ,k-1) - er(i  ,j-1,k-1))
                   end if

                   if      (er(i  ,j-1,k  ) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j+1,k  ) - er(i  ,j,  k  ))
                   else if (er(i  ,j+1,k  ) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j  ,k  ) - er(i  ,j-1,k  ))
                   end if

                else

                   dal = 0.e0_rt
                   dar = 0.e0_rt
                   dbl = 0.e0_rt
                   dbr = 0.e0_rt

                end if

                if (er(i,j,k-1) == -1.e0_rt) then

                   rg = ((er(i,j,k+1) - er(i,j,k)) * dxInv(3))**2 + &
                        (HALF * dar * dxInv(1))**2 + (HALF * dbr * dxInv(2))**2

                else if (er(i,j,k) == -1.e0_rt) then

                   rg = ((er(i,j,k-1) - er(i,j,k-2)) * dxInv(3))**2 + &
                        (HALF * dal * dxInv(1))**2 + (HALF * dbl * dxInv(2))**2

                else

                   rg = ((er(i,j,k) - er(i,j,k-1)) * dxInv(3))**2 + &
                        (FOURTH * (dal + dar) * dxInv(1))**2 + (FOURTH * (dbl + dbr) * dxInv(2))**2

                end if

                kap = kavg(kappar(i,j,k-1), kappar(i,j,k), dx(3), -1)
                r(i,j,k) = sqrt(rg) / (kap * max(er(i,j,k-1), er(i,j,k), tiny))

             end do
          end do
       end do

    end if

  end subroutine scgrd



  function kavg(a, b, d, opt) result(k)

    use amrex_fort_module, only: rt => amrex_real
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    implicit none

    real(rt), intent(in   ) :: a, b, d
    integer,  intent(in   ) :: opt

    real(rt) :: k

    !$gpu

#ifndef AMREX_USE_GPU
    if (opt > 2) then
       call castro_error("kavg: invalid averaging option")
    end if
#endif

    if (opt == 0) then

       ! arithmetic average, geometrically correct(?) but underestimates surface flux
       k = 0.5e0_rt * (a + b + tiny)

    else if (opt == 1) then

       ! harmonic average, overestimates surface flux
       k = (2.e0_rt * a * b) / (a + b + tiny) + tiny

    else

       ! chooses arithmetic where optically thin, harmonic where optically thick,
       ! surface flux approximation at a thick/thin boundary
       k = min(0.5e0_rt * (a + b + tiny), &
               max((2.e0_rt * a * b) / (a + b + tiny) + tiny, &
                   4.e0_rt / (3.e0_rt * d)))

    end if

  end function kavg



  subroutine lrhs(lo, hi, &
                  rhs, r_lo, r_hi, &
                  temp, t_lo, t_hi, &
                  fkp, fk_lo, fk_hi, &
                  eta, et_lo, et_hi, &
                  etainv, ei_lo, ei_hi, &
                  frhoem, fm_lo, fm_hi, &
                  frhoes, fs_lo, fs_hi, &
                  dfo, d_lo, d_hi, &
                  ero, er_lo, er_hi, &
                  edot, ed_lo, ed_hi, &
                  dt, dx, sigma, c, theta) &
                  bind(C, name="lrhs")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
    use habec_nd_module, only: cell_center_metric

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: fk_lo(3), fk_hi(3)
    integer,  intent(in   ) :: et_lo(3), et_hi(3)
    integer,  intent(in   ) :: ei_lo(3), ei_hi(3)
    integer,  intent(in   ) :: fm_lo(3), fm_hi(3)
    integer,  intent(in   ) :: fs_lo(3), fs_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: er_lo(3), er_hi(3)
    integer,  intent(in   ) :: ed_lo(3), ed_hi(3)
    real(rt), intent(inout) :: rhs(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: fkp(fk_lo(1):fk_hi(1),fk_lo(2):fk_hi(2),fk_lo(3):fk_hi(3))
    real(rt), intent(in   ) :: eta(et_lo(1):et_hi(1),et_lo(2):et_hi(2),et_lo(3):et_hi(3))
    real(rt), intent(in   ) :: etainv(ei_lo(1):ei_hi(1),ei_lo(2):ei_hi(2),ei_lo(3):ei_hi(3))
    real(rt), intent(in   ) :: frhoem(fm_lo(1):fm_hi(1),fm_lo(2):fm_hi(2),fm_lo(3):fm_hi(3))
    real(rt), intent(in   ) :: frhoes(fs_lo(1):fs_hi(1),fs_lo(2):fs_hi(2),fs_lo(3):fs_hi(3))
    real(rt), intent(in   ) :: dfo(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    real(rt), intent(in   ) :: ero(er_lo(1):er_hi(1),er_lo(2):er_hi(2),er_lo(3):er_hi(3))
    real(rt), intent(in   ) :: edot(ed_lo(1):ed_hi(1),ed_lo(2):ed_hi(2),ed_lo(3):ed_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, sigma, c, theta

    integer  :: i, j, k
    real(rt) :: dtm, ek, bs, es, ekt, r, s

    !$gpu

    dtm = 1.e0_rt / dt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ek = fkp(i,j,k) * eta(i,j,k)
             bs = etainv(i,j,k) * 4.e0_rt * sigma * fkp(i,j,k) * temp(i,j,k)**4
             es = eta(i,j,k) * (frhoem(i,j,k) - frhoes(i,j,k))
             ekt = (1.e0_rt - theta) * eta(i,j,k)

             call cell_center_metric(i, j, k, dx, r, s)

             if (dim == 1) then
                s = 1.e0_rt
             end if

             rhs(i,j,k) = (rhs(i,j,k) + r * s * &
                           (bs + dtm * (ero(i,j,k) + es) + &
                            ek * c * edot(i,j,k) - &
                            ekt * dfo(i,j,k))) / (1.e0_rt - ekt)

          end do
       end do
    end do

  end subroutine lrhs



  subroutine lacoef(lo, hi, &
                    a, a_lo, a_hi, &
                    fkp, f_lo, f_hi, &
                    eta, et_lo, et_hi, &
                    etainv, ei_lo, ei_hi, &
                    dx, c, dt, theta) &
                    bind(C, name="lacoef")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
    use habec_nd_module, only: cell_center_metric

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: et_lo(3), et_hi(3)
    integer,  intent(in   ) :: ei_lo(3), ei_hi(3)
    real(rt), intent(inout) :: a(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    real(rt), intent(in   ) :: fkp(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in   ) :: eta(et_lo(1):et_hi(1),et_lo(2):et_hi(2),et_lo(3):et_hi(3))
    real(rt), intent(in   ) :: etainv(ei_lo(1):ei_hi(1),ei_lo(2):ei_hi(2),ei_lo(3):ei_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: c, dt, theta

    integer  :: i, j, k
    real(rt) :: dtm, r, s

    !$gpu

    dtm = 1.e0_rt / dt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             call cell_center_metric(i, j, k, dx, r, s)

             if (dim == 1) then
                s = 1.e0_rt
             end if

             a(i,j,k) = r * s * &
                  (fkp(i,j,k) * etainv(i,j,k) * c + dtm) / &
                  (1.e0_rt - (1.e0_rt - theta) * eta(i,j,k))

          end do
       end do
    end do

  end subroutine lacoef



  subroutine bclim(lo, hi, &
                   b, b_lo, b_hi, &
                   lambda, l_lo, l_hi, &
                   n, &
                   kappar, k_lo, k_hi, &
                   c, dx) bind(C, name="bclim")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
    use habec_nd_module, only: edge_center_metric

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    integer,  intent(in   ) :: l_lo(3), l_hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    real(rt), intent(inout) :: b(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    real(rt), intent(in   ) :: lambda(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3))
    real(rt), intent(in   ) :: kappar(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: n
    real(rt), intent(in   ), value :: c

    integer  :: i, j, k
    real(rt) :: kap, r, s

    !$gpu

    if (n == 0) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                call edge_center_metric(i, j, k, n+1, dx, r, s)

                if (dim == 1) then
                   s = 1.e0_rt
                end if

                kap = kavg(kappar(i-1,j,k), kappar(i,j,k), dx(1), -1)
                b(i,j,k) = r * s * c * lambda(i,j,k) / kap

             end do
          end do
       end do

    else if (n == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                call edge_center_metric(i, j, k, n+1, dx, r, s)

                if (dim == 1) then
                   s = 1.e0_rt
                end if

                kap = kavg(kappar(i,j-1,k), kappar(i,j,k), dx(2), -1)
                b(i,j,k) = r * s * c * lambda(i,j,k) / kap

             end do
          end do
       end do

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                call edge_center_metric(i, j, k, n+1, dx, r, s)

                if (dim == 1) then
                   s = 1.e0_rt
                end if
                
                kap = kavg(kappar(i,j,k-1), kappar(i,j,k), dx(3), -1)
                b(i,j,k) = r * s * c * lambda(i,j,k) / kap

             end do
          end do
       end do

    end if

  end subroutine bclim



  subroutine ceup(lo, hi, &
                  relres, absres, &
                  frhoes, fs_lo, fs_hi, &
                  frhoem, fm_lo, fm_hi, &
                  eta, et_lo, et_hi, &
                  etainv, ei_lo, ei_hi, &
                  dfo, do_lo, do_hi, &
                  dfn, dn_lo, dn_hi, &
                  exch, ex_lo, ex_hi, &
                  dt, theta) &
                  bind(C, name="ceup")

    use amrex_fort_module, only: rt => amrex_real
    use reduction_module, only: reduce_max

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: fs_lo(3), fs_hi(3)
    integer,  intent(in   ) :: fm_lo(3), fm_hi(3)
    integer,  intent(in   ) :: et_lo(3), et_hi(3)
    integer,  intent(in   ) :: ei_lo(3), ei_hi(3)
    integer,  intent(in   ) :: do_lo(3), do_hi(3)
    integer,  intent(in   ) :: dn_lo(3), dn_hi(3)
    integer,  intent(in   ) :: ex_lo(3), ex_hi(3)
    real(rt), intent(inout) :: frhoes(fs_lo(1):fs_hi(1),fs_lo(2):fs_hi(2),fs_lo(3):fs_hi(3))
    real(rt), intent(in   ) :: frhoem(fm_lo(1):fm_hi(1),fm_lo(2):fm_hi(2),fm_lo(3):fm_hi(3))
    real(rt), intent(in   ) :: eta(et_lo(1):et_hi(1),et_lo(2):et_hi(2),et_lo(3):et_hi(3))
    real(rt), intent(in   ) :: etainv(ei_lo(1):ei_hi(1),ei_lo(2):ei_hi(2),ei_lo(3):ei_hi(3))
    real(rt), intent(in   ) :: dfo(do_lo(1):do_hi(1),do_lo(2):do_hi(2),do_lo(3):do_hi(3))
    real(rt), intent(in   ) :: dfn(dn_lo(1):dn_hi(1),dn_lo(2):dn_hi(2),dn_lo(3):dn_hi(3))
    real(rt), intent(in   ) :: exch(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2),ex_lo(3):ex_hi(3))
    real(rt), intent(inout) :: relres, absres
    real(rt), intent(in   ), value :: dt, theta

    integer  :: i, j, k
    real(rt) :: tmp, chg, tot

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             chg = 0.e0_rt
             tot = 0.e0_rt

             tmp = eta(i,j,k) * frhoes(i,j,k) + &
                   etainv(i,j,k) * &
                   (frhoem(i,j,k) - &
                    dt * ((1.e0_rt - theta) * &
                    (dfo(i,j,k) - dfn(i,j,k)) + &
                    exch(i,j,k)))

             chg = abs(tmp - frhoes(i,j,k))
             tot = abs(frhoes(i,j,k))

             frhoes(i,j,k) = tmp

             call reduce_max(absres, chg)
             call reduce_max(relres, chg / (tot + tiny))

          end do
       end do
    end do

  end subroutine ceup



 
  subroutine cetot(lo, hi, &
                   state, s_lo, s_hi, &
                   frhoe, f_lo, f_hi) &
                   bind(C, name="cetot")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, UEINT, UEDEN

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: frhoe(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))

    integer  :: i, j, k
    real(rt) :: kin

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             kin = state(i,j,k,UEDEN) - state(i,j,k,UEINT)
             state(i,j,k,UEINT) = frhoe(i,j,k)
             state(i,j,k,UEDEN) = frhoe(i,j,k) + kin
          end do
       end do
    end do

  end subroutine cetot



  subroutine ca_compute_rosseland(lo, hi, &
                                  kpr, k_lo, k_hi, &
                                  state, s_lo, s_hi) &
                                  bind(C, name="ca_compute_rosseland")

    use rad_params_module, only: ngroups, nugroup
    use opacity_table_module, only: get_opacities
    use network, only: naux
    use meth_params_module, only: NVAR, URHO, UTEMP, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: kpr(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer  :: i, j, k, g
    real(rt) :: kp, kr, nu, rho, temp, Ye
    logical, parameter :: comp_kp = .false. 
    logical, parameter :: comp_kr = .true.

    !$gpu

    do g = 0, ngroups-1

       nu = nugroup(g)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                rho = state(i,j,k,URHO)
                temp = state(i,j,k,UTEMP)
                if (naux > 0) then
                   Ye = state(i,j,k,UFX)
                else
                   Ye = 0.e0_rt
                end if

                call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)

                kpr(i,j,k,g) = kr

             end do
          end do
       end do
    end do

  end subroutine ca_compute_rosseland



  subroutine ca_compute_planck(lo, hi, &
                               kpp, k_lo, k_hi, &
                               state, s_lo, s_hi, &
                               temp_offset) &
                               bind(C, name="ca_compute_planck")

    use rad_params_module, only: ngroups, nugroup
    use opacity_table_module, only: get_opacities
    use network, only: naux
    use meth_params_module, only: NVAR, URHO, UTEMP, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: kpp(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ), value :: temp_offset

    integer  :: i, j, k, g
    real(rt) :: kp, kr, nu, rho, temp, Ye
    logical, parameter :: comp_kp = .true. 
    logical, parameter :: comp_kr = .false.

    !$gpu

    do g = 0, ngroups-1

       nu = nugroup(g)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                rho = state(i,j,k,URHO)
                temp = state(i,j,k,UTEMP) + temp_offset
                if (naux > 0) then
                   Ye = state(i,j,k,UFX)
                else
                   Ye = 0.e0_rt
                end if

                call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)

                kpp(i,j,k,g) = kp

             end do
          end do
       end do
    end do

  end subroutine ca_compute_planck



  subroutine fkpn(lo, hi, &
                  fkp, f_lo, f_hi, &
                  const, em, en, &
                  ep, nu, tf, &
                  temp, t_lo, t_hi, &
                  state, s_lo, s_hi, &
                  temp_offset) &
                  bind(C, name="fkpn")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: fkp(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ), value :: const, em, en, tf, ep, nu, temp_offset

    real(rt) :: teff
    integer  :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             teff = max(temp(i,j,k) + temp_offset, tiny)
             teff = teff + tf * exp(-teff / (tf + tiny))

             fkp(i,j,k) = const * &
                          (state(i,j,k,URHO)**em) * &
                          (teff**(-en)) * &
                          (nu**(ep))

          end do
       end do
    end do

  end subroutine fkpn



  subroutine nfloor(lo, hi, &
                    dest, d_lo, d_hi, &
                    flr, nvar) &
                    bind(C, name="nfloor")

    use amrex_fort_module, only: rt => amrex_real

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    real(rt), intent(inout) :: dest(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nvar-1)
    integer,  intent(in   ), value :: nvar
    real(rt), intent(in   ), value :: flr

    integer :: i, j, k, n

    !$gpu

    do n = 0, nvar-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (dest(i,j,k,n) < flr) then
                   dest(i,j,k,n) = flr
                end if
             end do
          end do
       end do
    end do

  end subroutine nfloor



  subroutine ceta2(lo, hi, &
                   eta, et_lo, et_hi, &
                   etainv, ei_lo, ei_hi, &
                   frho, fr_lo, fr_hi, &
                   temp, t_lo, t_hi, &
                   cv, c_lo, c_hi, &
                   fkp, fk_lo, fk_hi, &
                   er, er_lo, er_hi, &
                   dtemp, dtime, sigma, &
                   c, underr, lagpla) &
                   bind(C, name="ceta2")

    use amrex_fort_module, only: rt => amrex_real

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: et_lo(3), et_hi(3)
    integer,  intent(in   ) :: ei_lo(3), ei_hi(3)
    integer,  intent(in   ) :: fr_lo(3), fr_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: fk_lo(3), fk_hi(3)
    integer,  intent(in   ) :: er_lo(3), er_hi(3)
    real(rt), intent(inout) :: eta(et_lo(1):et_hi(1),et_lo(2):et_hi(2),et_lo(3):et_hi(3))
    real(rt), intent(inout) :: etainv(ei_lo(1):ei_hi(1),ei_lo(2):ei_hi(2),ei_lo(3):ei_hi(3))
    real(rt), intent(in   ) :: frho(fr_lo(1):fr_hi(1),fr_lo(2):fr_hi(2),fr_lo(3):fr_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: cv(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(in   ) :: fkp(fk_lo(1):fk_hi(1),fk_lo(2):fk_hi(2),fk_lo(3):fk_hi(3))
    real(rt), intent(in   ) :: er(er_lo(1):er_hi(1),er_lo(2):er_hi(2),er_lo(3):er_hi(3))
    real(rt), intent(in   ), value :: dtemp, dtime, sigma, c, underr
    integer,  intent(in   ), value :: lagpla

    real(rt) :: d, frc, fac0, fac1, fac2
    integer  :: i, j, k

    !$gpu

    fac1 = 16.e0_rt * sigma * dtime
    if (lagpla == 0) then
       fac0 = 0.25e0_rt * fac1 / dtemp
       fac2 = dtime * c / dtemp
    endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (lagpla /= 0) then

                ! assume eta and fkp are the same
                d = fac1 * fkp(i,j,k) * temp(i,j,k) ** 3

             else

                d = fac0 * (eta(i,j,k) * (temp(i,j,k) + dtemp) ** 4 - &
                            fkp(i,j,k) * (temp(i,j,k)        ) ** 4) - &
                            fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
                ! alternate form, sometimes worse, sometimes better:
                !                  d = fac1 * fkp(i,j,k) * temp(i,j,k) ** 3 +
                !     @                fac0 * (eta(i,j,k) - fkp(i,j,k)) * temp(i,j,k) ** 4 -
                !     @                fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
                ! another alternate form (much worse):
                !                  d = fac1 * fkp(i,j,k) * (temp(i,j,k) + dtemp) ** 3 +
                !     @                fac0 * (eta(i,j,k) - fkp(i,j,k))
                !     @                     * (temp(i,j,k) + dtemp) ** 4 -
                !     @                fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
             end if

             frc = frho(i,j,k) * cv(i,j,k) + tiny
             eta(i,j,k) = d / (d + frc)
             etainv(i,j,k) = underr * frc / (d + frc)
             eta(i,j,k) = 1.e0_rt - etainv(i,j,k)
             !               eta(i,j,k) = 1.e0_rt - underr * (1.e0_rt - eta(i,j,k))

          end do
       end do
    end do

  end subroutine ceta2



  subroutine gcv(lo, hi, &
                 cv, c_lo, c_hi, &
                 temp, t_lo, t_hi, &
                 const, em, en, tf, &
                 state, s_lo, s_hi) bind(C, name="gcv")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: cv(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! temp contains temp on input
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ), value :: const, em, en, tf

    real(rt) :: alpha, teff, frhoal
    integer  :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (em == 0.e0_rt) then
                alpha = const
             else
                alpha = const * state(i,j,k,URHO)**em
             end if

             frhoal = state(i,j,k,URHO) * alpha + tiny

             if (en == 0.e0_rt) then
                cv(i,j,k) = alpha
             else
                teff = max(temp(i,j,k), tiny)
                teff = teff + tf * exp(-teff / (tf + tiny))
                cv(i,j,k) = alpha * teff**(-en)
             end if

          end do
       end do
    end do

  end subroutine gcv



  subroutine ca_compute_c_v(lo, hi, &
                            cv, c_lo, c_hi, &
                            temp, t_lo, t_hi, &
                            state, s_lo, s_hi) &
                            bind(C, name="ca_compute_c_v")

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UFS, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: cv(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer     :: i, j, k
    real(rt)    :: rhoInv
    type(eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = 1.e0_rt / state(i,j,k,URHO)
             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   = temp(i,j,k)
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_rt, eos_state)

             cv(i,j,k) = eos_state % cv

          end do
       end do
    end do

  end subroutine ca_compute_c_v



  subroutine ca_get_rhoe(lo, hi, &
                         rhoe, r_lo, r_hi, &
                         temp, t_lo, t_hi, &
                         state, s_lo, s_hi) &
                         bind(C, name="ca_get_rhoe")

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UFS, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: rhoe(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer     :: i, j, k
    real(rt)    :: rhoInv
    type(eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = 1.e0_rt / state(i,j,k,URHO)
             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   =  temp(i,j,k)
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = state(i,j,k,UFX:UFX+naux -1) * rhoInv

             call eos(eos_input_rt, eos_state)

             rhoe(i,j,k) = eos_state % rho * eos_state % e

          end do
       end do
    end do

  end subroutine ca_get_rhoe



  subroutine gtemp(lo, hi, &
                   temp, t_lo, t_hi, &
                   const, em, en, &
                   state, s_lo, s_hi, &
                   update_state) bind(C, name="gtemp")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO, UTEMP
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))  ! temp contains frhoe on input
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ), value :: const, em, en
    integer,  intent(in   ), value :: update_state

    real(rt) :: alpha, teff, ex, frhoal
    integer  :: i, j, k

    !$gpu

#ifndef AMREX_USE_GPU
    if (en >= 1.e0_rt) then
       call castro_error("Bad exponent for cv calculation")
    end if
#endif

    ex = 1.e0_rt / (1.e0_rt - en)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (em == 0.e0_rt) then
                alpha = const
             else
                alpha = const * state(i,j,k,URHO)**em
             end if

             frhoal = state(i,j,k,URHO) * alpha + tiny

             if (en == 0.e0_rt) then
                temp(i,j,k) = temp(i,j,k) / frhoal
             else
                teff = max(temp(i,j,k), tiny)
                temp(i,j,k) = ((1.e0_rt - en) * teff / frhoal)**ex
             end if

             if (update_state == 1) then
                state(i,j,k,UTEMP) = temp(i,j,k)
             end if

          end do
       end do
    end do

  end subroutine gtemp



  subroutine ca_compute_temp_given_rhoe(lo, hi, &
                                        temp, t_lo, t_hi, &
                                        state, s_lo, s_hi, &
                                        update_state) &
                                        bind(C, name="ca_compute_temp_given_rhoe")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: NVAR, URHO, UTEMP, UFS, UFX, small_temp
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! temp contains rhoe as input
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer,  intent(in   ), value :: update_state

    integer      :: i, j, k
    real(rt)     :: rhoInv
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (temp(i,j,k) .le. 0.e0_rt) then

                temp(i,j,k) = small_temp

             else

                rhoInv = 1.e0_rt / state(i,j,k,URHO)
                eos_state % rho = state(i,j,k,URHO)
                eos_state % T   = state(i,j,k,UTEMP)
                eos_state % e   =  temp(i,j,k)*rhoInv 
                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
                eos_state % aux = state(i,j,k,UFX:UFX+naux -1) * rhoInv

                call eos(eos_input_re, eos_state)

                temp(i,j,k) = eos_state % T

             end if

             if (update_state == 1) then
                state(i,j,k,UTEMP) = temp(i,j,k)
             end if

          end do
       end do
    end do

  end subroutine ca_compute_temp_given_rhoe



  subroutine ca_compute_temp_given_cv(lo, hi, &
                                      temp, t_lo, t_hi, &
                                      state, s_lo, s_hi, &
                                      const_c_v, c_v_exp_m, c_v_exp_n, &
                                      update_state) &
                                      bind(C, name="ca_compute_temp_given_cv")

    use meth_params_module, only: NVAR, URHO, UTEMP
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! temp contains rhoe as input
    real(rt), intent(in   ), value :: const_c_v, c_v_exp_m, c_v_exp_n
    integer,  intent(in   ), value :: update_state

    integer  :: i, j, k
    real(rt) :: ex, alpha, rhoal, teff

    !$gpu

    ex = 1.e0_rt / (1.e0_rt - c_v_exp_n)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (c_v_exp_m .eq. 0.e0_rt) then
                alpha = const_c_v
             else
                alpha = const_c_v * state(i,j,k,URHO)**c_v_exp_m
             endif

             rhoal = state(i,j,k,URHO) * alpha + 1.e-50_rt

             if (c_v_exp_n .eq. 0.e0_rt) then
                temp(i,j,k) = temp(i,j,k) / rhoal
             else
                teff = max(temp(i,j,k), 1.e-50_rt)
                temp(i,j,k) = ((1.e0_rt - c_v_exp_n) * teff / rhoal)**ex
             end if

             if (update_state == 1) then
                state(i,j,k,UTEMP) = temp(i,j,k)
             end if

          end do
       end do
    end do

  end subroutine ca_compute_temp_given_cv



  subroutine cfrhoe(lo, hi, &
                    frhoe, f_lo, f_hi, &
                    state, s_lo, s_hi) &
                    bind(C, name='cfrhoe')

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, UEINT

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: frhoe(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ! kin = 0.5e0_rt * (state(i,j,k,XMOM)   ** 2 +
             !                   state(i,j,k,XMOM+1) ** 2 +
             !                   state(i,j,k,XMOM+2) ** 2) /
             !                   state(i,j,k,DEN)
             ! frhoe(i,j,k) = state(i,j,k,EDEN) - kin
             frhoe(i,j,k) = state(i,j,k,UEINT)
          end do
       end do
    end do

  end subroutine cfrhoe



  subroutine rosse1(lo, hi, &
                    const, em, en, &
                    ep, nu, &
                    tf, kfloor, &
                    state, s_lo, s_hi, &
                    kappar, k_lo, k_hi) bind(C, name="rosse1")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO, UTEMP

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: kappar(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3))
    real(rt), intent(in   ), value :: const, em, en, ep, nu, tf, kfloor

    real(rt) :: kf, teff
    integer  :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             teff = max(state(i,j,k,UTEMP), tiny)
             teff = teff + tf * exp(-teff / (tf + tiny))

             kf = const * &
                  (state(i,j,k,URHO) ** em) * &
                  (teff ** (-en)) * &
                  (nu ** (ep))

             kappar(i,j,k) = max(kf, kfloor)

          end do
       end do
    end do

  end subroutine rosse1



  subroutine rosse1s(lo, hi, &
                     const, em, en, &
                     ep, sconst, &
                     sem, sen, &
                     sep, nu, &
                     tf, kfloor, &
                     state, s_lo, s_hi, &
                     kappar, k_lo, k_hi) bind(C, name="rosse1s")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO, UTEMP

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    real(rt), intent(inout) :: kappar(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ), value :: const, em, en, ep, sconst, sem, sen, sep, nu, tf, kfloor

    real(rt) :: kf, teff, sct
    integer  :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             teff = max(state(i,j,k,UTEMP), tiny)
             teff = teff + tf * exp(-teff / (tf + tiny))

             kf = const * &
                  (state(i,j,k,URHO) ** em) * &
                  (teff ** (-en)) * &
                  (nu ** (ep))

             sct = sconst * &
                  (state(i,j,k,URHO) ** sem) * &
                  (teff ** (-sen)) * &
                  (nu ** (sep))

             kappar(i,j,k) = max(kf + sct, kfloor)

          end do
       end do
    end do

  end subroutine rosse1s



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The following routined are used by neutrinos only.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ca_compute_temp_given_reye(lo, hi, &
                                        temp, t_lo, t_hi, &
                                        rhoe, r_lo, r_hi, &
                                        ye, y_lo, y_hi, &
                                        state, s_lo, s_hi) &
                                        bind(C, name='ca_compute_temp_given_reye')

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UFS, UFX, small_temp
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: y_lo(3), y_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: rhoe(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: ye(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer      :: i, j, k
    real(rt)     :: rhoInv
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (rhoe(i,j,k) .le. 0.e0_rt) then

                temp(i,j,k) = small_temp

             else

                rhoInv = 1.e0_rt / state(i,j,k,URHO)
                eos_state % rho = state(i,j,k,URHO)
                ! set initial guess of temperature
                eos_state % T = temp(i,j,k)
                eos_state % e = rhoe(i,j,k)*rhoInv 
                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
                if (naux > 0) then
                   eos_state % aux = ye(i,j,k)
                end if

                call eos(eos_input_re, eos_state)

                temp(i,j,k) = eos_state % T

#ifndef AMREX_USE_GPU
                if (temp(i,j,k) .lt. 0.e0_rt) then
                   print *, 'negative temp in compute_temp_given_reye ', temp(i,j,k)
                   call castro_error("Error:: ca_compute_temp_given_reye")
                endif
#endif

             end if

          end do
       end do
    end do

  end subroutine ca_compute_temp_given_reye



  subroutine ca_compute_reye_given_ty(lo, hi, &
                                      rhoe, re_lo, re_hi, &
                                      rhoY, rY_lo, rY_hi, &
                                      temp, t_lo, t_hi, &
                                      ye, y_lo, y_hi, &
                                      state, s_lo, s_hi) &
                                      bind(C, name='ca_compute_reye_given_ty')

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use meth_params_module, only: NVAR, URHO, UFS, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: re_lo(3), re_hi(3)
    integer,  intent(in   ) :: rY_lo(3), ry_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: y_lo(3), y_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: rhoe(re_lo(1):re_hi(1),re_lo(2):re_hi(2),re_lo(3):re_hi(3))
    real(rt), intent(inout) :: rhoY(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: ye(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer      :: i, j, k
    real(rt)     :: rhoInv
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = 1.e0_rt / state(i,j,k,URHO)
             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   = temp(i,j,k)
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv

             if (naux > 0) then
                eos_state % aux = ye(i,j,k)
                rhoY(i,j,k) = state(i,j,k,URHO)*ye(i,j,k)        
             end if

             call eos(eos_input_rt, eos_state)

             rhoe(i,j,k) = eos_state % rho * eos_state % e

          end do
       end do
    end do

  end subroutine ca_compute_reye_given_ty

end module rad_nd_module



!! -----------------------------------------------------------
!> @brief This routine is called at problem setup time and is used
!! to initialize values of physical constants used by the
!! radiation package.
!! -----------------------------------------------------------
subroutine ca_initradconstants(p, c, h, k, s, a, m, J_is_used) bind(C, name="ca_initradconstants")

  use fundamental_constants_module, only: c_fcm=>c_light, h_fcm=>hplanck, &
                                          k_fcm=>k_B, s_fcm=>sigma_SB, a_fcm=>n_A, ev2erg_fcm=>ev2erg
  use rad_params_module, only: pi, clight, hplanck
  use rad_params_module, only: kboltz, stefbol, arad, avogadro
  use rad_params_module, only: Hz2MeV, mev2erg, tiny
  use rad_params_module, only: radtoE  !, radtoJ, Etorad, radfluxtoF
  use rad_params_module, only: etafactor
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt) :: p, c, h, k, s, a, m
  integer :: J_is_used

  c = c_fcm
  h = h_fcm
  k = k_fcm
  s = s_fcm
  a = a_fcm
  m = 1.e6_rt * ev2erg_fcm

  pi       = p
  clight   = c
  hplanck  = h
  kboltz   = k
  stefbol  = s
  arad     = 4.*stefbol/clight
  avogadro = a
  mev2erg  = m
  Hz2MeV   = h / m
  tiny     = 1.e-50_rt

  if (J_is_used > 0) then
     radtoE = 4.e0_rt*pi/clight
     !           radtoJ = 1.0e0_rt
     !           Etorad = 1.e0_rt/radtoE
     !           radfluxtoF = 4.e0_rt*pi
     etafactor = 1.e0_rt
  else
     radtoE = 1.0e0_rt
     !           radtoJ = clight/(4.e0_rt*pi)
     !           Etorad = 1.0e0_rt
     !           radfluxtoF = 1.e0_rt
     etafactor = 4.e0_rt*pi/clight
  end if

end subroutine ca_initradconstants

! For single group, let set ngroups to 1.
subroutine ca_initsinglegroup(ngr) bind(C, name="ca_initsinglegroup")

  use rad_params_module, only: ngroups, nugroup, dnugroup, ng0, ng1
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer :: ngr

  ! Local variables
  integer :: i

  ng0 = 0
  ng1 = 0

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))

  do i = 0, ngroups-1
     nugroup(i)  = 1.e0_rt  ! dummy
     dnugroup(i) = 1.e0_rt
  enddo

end subroutine ca_initsinglegroup

!! -----------------------------------------------------------
!> @brief This routine is called at problem setup time and is used
!! to initialize the arrays nugroup and dnugroup in
!! probdata with the neutrino group energies and widths.
!!
!! The widths are used to derive neutrino spectrum for plot files
!! -----------------------------------------------------------
subroutine ca_initgroups(nugr, dnugr, ngr, ngr0, ngr1)

  use rad_params_module, only: ngroups, ng0, ng1, nugroup, dnugroup, &
                               current_group, nnuspec
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt) :: nugr(0:ngr-1), dnugr(0:ngr-1)
  integer :: ngr, ngr0, ngr1

  ! Local variables
  integer :: i

  allocate(current_group, ng0, ng1, nnuspec)

  ng0     = ngr0
  ng1     = ngr1

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))

  do i = 0, ngroups-1
     nugroup(i)  = nugr(i)
     dnugroup(i) = dnugr(i)
  enddo

end subroutine ca_initgroups

subroutine ca_initgroups2(nugr, dnugr, xnugr, ngr)

  use rad_params_module, only: ngroups, nugroup, dnugroup, xnu, dlognu, lognugroup, &
                               current_group, ng0, ng1, nnuspec
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), intent(in) :: nugr(0:ngr-1), dnugr(0:ngr-1), xnugr(0:ngr)
  integer :: ngr

  ! Local variables
  integer   :: i

  allocate(current_group, ng0, ng1, nnuspec)

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))
  allocate(xnu(0:ngroups))
  allocate(dlognu(0:ngroups-1))
  allocate(lognugroup(0:ngroups-1))

  nugroup(:) = nugr(:)
  dnugroup(:) = dnugr(:)
  xnu(:) = xnugr(:)
  lognugroup(:) = log(nugroup)

  dlognu(0:ngroups-1) = log(xnu(1:ngroups)) - log(xnu(0:ngroups-1))

end subroutine ca_initgroups2

subroutine ca_initgroups3(nugr, dnugr, dlognugr, xnugr, ngr, ngr0, ngr1)
  ! used by MGFLDSolver

  use rad_params_module, only: ngroups, ng0, ng1, nnuspec, nradspec, nugroup, dnugroup, &
                               xnu, dlognu, lognugroup, erg2rhoYe, avogadro, hplanck, &
                               current_group
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), intent(in) :: nugr(0:ngr-1), dnugr(0:ngr-1), dlognugr(0:ngr-1), xnugr(0:ngr+2)
  integer :: ngr, ngr0, ngr1

  ! Local variables
  integer :: i

  allocate(current_group, ng0, ng1, nnuspec)

  ng0     = ngr0
  ng1     = ngr1

  if (ng0 > 0) then
     if (ng1 .eq. 0) then
        nnuspec = 1  ! one neutrino species
     else if (ngroups .eq. ng0+ng1) then
        nnuspec = 2  ! two neutrino species
     else
        nnuspec = 3  ! three neutrino species
     end if
  else
     nnuspec = 0
  end if

  nradspec = max(nnuspec, 1)

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))
  allocate(xnu(0:ngroups+2))
  allocate(dlognu(0:ngroups-1))
  allocate(erg2rhoYe(0:ngroups-1))
  allocate(lognugroup( 0:ngroups-1))

  nugroup(:) = nugr(:)
  dnugroup(:) = dnugr(:)
  xnu(:) = xnugr(:)
  dlognu(:) = dlognugr(:)
  lognugroup(:) = log(nugroup)

  erg2rhoYe = 0.e0_rt
  if (ng0 > 0) then
     erg2rhoYe(0:ng0-1) = 1.e0_rt / (avogadro*hplanck*nugroup(0:ng0-1))
     if (ng1 > 0) then
        erg2rhoYe(ng0:ng0+ng1-1) = -1.e0_rt / (avogadro*hplanck*nugroup(ng0:ng0+ng1-1))
     end if
  end if

end subroutine ca_initgroups3

!! -----------------------------------------------------------

subroutine ca_setgroup(igroup)

  use rad_params_module, only: current_group
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer :: igroup

  current_group = igroup

end subroutine ca_setgroup

!! -----------------------------------------------------------

subroutine ca_inelastic_sct(lo, hi, &
                            uu,uu_lo,uu_hi, &
                            Er,Er_lo,Er_hi, &
                            ks,ks_lo,ks_hi, &
                            dt) bind(C)

  use meth_params_module, only: NVAR, UEDEN, UEINT, UTEMP
  use rad_params_module, only: ngroups, nugroup, dlognu
  use radhydro_nd_module, only: inelastic_scatter
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: uu_lo(3),uu_hi(3)
  integer, intent(in) :: Er_lo(3),Er_hi(3)
  integer, intent(in) :: ks_lo(3),ks_hi(3)
  real(rt), intent(inout) :: uu(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3),NVAR)
  real(rt), intent(inout) :: Er(Er_lo(1):Er_hi(1),Er_lo(2):Er_hi(2),Er_lo(3):Er_hi(3),0:ngroups-1)
  real(rt), intent(in   ) :: ks(ks_lo(1):ks_hi(1),ks_lo(2):ks_hi(2),ks_lo(3):ks_hi(3))
  real(rt), intent(in) :: dt

  integer :: i, j, k
  real(rt) :: Ertotold, Ertmp(0:ngroups-1), dEr
  real(rt) :: Erscale(0:ngroups-1)

  Erscale = nugroup*dlognu

  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           Ertmp = Er(i,j,k,:)
           Ertotold = sum(Ertmp)
           Ertmp = Ertmp / Erscale

           call inelastic_scatter(uu(i,j,k,UTEMP), Ertmp, ks(i,j,k), dt, (/i,j,k/))

           Ertmp = Ertmp * Erscale
           dEr = sum(Ertmp) - Ertotold
           Er(i,j,k,:) = Ertmp
           uu(i,j,k,UEINT) = uu(i,j,k,UEINT) - dEr
           uu(i,j,k,UEDEN) = uu(i,j,k,UEDEN) - dEr
        end do
     end do
  end do

end subroutine ca_inelastic_sct

!! -----------------------------------------------------------

subroutine ca_compute_scattering(lo, hi, &
                                 kps,kps_lo,kps_hi, &
                                 sta,sta_lo,sta_hi)

  use rad_params_module, only : ngroups, nugroup
  use opacity_table_module, only : get_opacities
  use network, only : naux
  use meth_params_module, only : NVAR, URHO, UTEMP, UFX
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: kps_lo(3),kps_hi(3)
  integer, intent(in) :: sta_lo(3),sta_hi(3)
  real(rt), intent(inout) :: kps(kps_lo(1):kps_hi(1),kps_lo(2):kps_hi(2),kps_lo(3):kps_hi(3))
  real(rt), intent(in   ) :: sta(sta_lo(1):sta_hi(2),sta_lo(2):sta_hi(2),sta_lo(3):sta_hi(3),NVAR)

  integer :: i, j, k
  real(rt) :: kp, kr, nu, rho, temp, Ye
  logical, parameter :: comp_kp = .true.
  logical, parameter :: comp_kr = .true.

  ! scattering is assumed to be independent of nu.

  nu = nugroup(0)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           rho = sta(i,j,k,URHO)
           temp = sta(i,j,k,UTEMP)
           if (naux > 0) then
              Ye = sta(i,j,k,UFX)
           else
              Ye = 0.e0_rt
           end if

           call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)

           kps(i,j,k) = max(kr - kp, 0.e0_rt)
        end do
     end do
  end do

end subroutine ca_compute_scattering

subroutine ca_compute_scattering_2(lo, hi, &
                                   kps,kps_lo,kps_hi, &
                                   sta,sta_lo,sta_hi, &
                                   k0_p, m_p, n_p, &
                                   k0_r, m_r, n_r, &
                                   Tfloor, kfloor)

  use rad_params_module, only: ngroups, nugroup
  use meth_params_module, only: NVAR, URHO, UTEMP
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: kps_lo(3),kps_hi(3)
  integer, intent(in) :: sta_lo(3),sta_hi(3)
  real(rt), intent(inout) :: kps(kps_lo(1):kps_hi(1),kps_lo(2):kps_hi(2),kps_lo(3):kps_hi(3))
  real(rt), intent(in   ) :: sta(sta_lo(1):sta_hi(2),sta_lo(2):sta_hi(2),sta_lo(3):sta_hi(3),NVAR)
  real(rt), intent(in) :: k0_p, m_p, n_p
  real(rt), intent(in) :: k0_r, m_r, n_r
  real(rt), intent(in) :: Tfloor, kfloor

  integer :: i, j, k
  real(rt), parameter :: tiny = 1.0e-50_rt
  real(rt) :: Teff, k_p, k_r

  ! scattering is assumed to be independent of nu.

  if ( m_p.eq.0.e0_rt .and. n_p.eq.0.e0_rt .and. &
       m_r.eq.0.e0_rt .and. n_r.eq.0.e0_rt ) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              kps(i,j,k) = k0_r - k0_p
           end do
        end do
     end do
  else
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              Teff = max(sta(i,j,k,UTEMP), tiny)
              Teff = Teff + Tfloor * exp(-Teff / (Tfloor + tiny))
              k_p = k0_p * (sta(i,j,k,URHO) ** m_p) * (Teff ** (-n_p))
              k_r = k0_r * (sta(i,j,k,URHO) ** m_r) * (Teff ** (-n_r))
              kps(i,j,k) = max(k_r-k_p, kfloor)
           end do
        end do
     end do
  end if

end subroutine ca_compute_scattering_2
