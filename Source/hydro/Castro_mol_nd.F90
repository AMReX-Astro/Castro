! advection routines in support of method of lines integration
!

module castro_mol_module

  implicit none

contains


  subroutine ca_mol_consup(lo, hi, &
                           shk, shk_lo, shk_hi, &
                           uin, uin_lo, uin_hi, &
                           srcU, srU_lo, srU_hi, &
                           update, updt_lo, updt_hi, &
                           dx, dt, &
                           flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
                           flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                           flux3, flux3_lo, flux3_hi, &
#endif
                           area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                           area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                           area3, area3_lo, area3_hi, &
#endif
                           q1, q1_lo, q1_hi, &
#if AMREX_SPACEDIM >= 2
                           q2, q2_lo, q2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                           q3, q3_lo, q3_hi, &
#endif
                           vol, vol_lo, vol_hi) bind(C, name="ca_mol_consup")

    use castro_error_module
    use meth_params_module, only : NQ, NVAR, NGDNV, NSRC, GDPRES, &
                                   UTEMP, UMX, &
                                   QPRES, &
                                   QTEMP, QFS, QFX, QREINT, QRHO, &
                                   first_order_hydro, difmag, hybrid_riemann, &
                                   limit_fluxes_on_small_dens, ppm_type, ppm_temp_fix, do_hydro
    use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
    use amrex_fort_module, only : rt => amrex_real
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : add_hybrid_advection_source
#endif
    use network, only : nspec, naux
    use prob_params_module, only : dg, coord_type


    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: shk_lo(3), shk_hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: srU_lo(3), srU_hi(3)
    integer, intent(in) :: updt_lo(3), updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: q1_lo(3), q1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) :: q2_lo(3), q2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)
    integer, intent(in) :: q3_lo(3), q3_hi(3)
#endif

    integer, intent(in) :: vol_lo(3), vol_hi(3)

    real(rt), intent(in) :: shk(shk_lo(1):shk_hi(1), shk_lo(2):shk_hi(2), shk_lo(3):shk_hi(3))
    real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(in) :: srcU(srU_lo(1):srU_hi(1), srU_lo(2):srU_hi(2), srU_lo(3):srU_hi(3), NSRC)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)

    real(rt), intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
    real(rt), intent(in) :: q1(q1_lo(1):q1_hi(1), q1_lo(2):q1_hi(2), q1_lo(3):q1_hi(3), NGDNV)

#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
    real(rt), intent(in) :: q2(q2_lo(1):q2_hi(1), q2_lo(2):q2_hi(2), q2_lo(3):q2_hi(3), NGDNV)
#endif

#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)
    real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
    real(rt), intent(in) :: q3(q3_lo(1):q3_hi(1), q3_lo(2):q3_hi(2), q3_lo(3):q3_hi(3), NGDNV)
#endif

    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

    integer :: i, j, k, n

    !$gpu

    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

#if AMREX_SPACEDIM == 1
                update(i,j,k,n) = update(i,j,k,n) + &
                     (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) ) / vol(i,j,k)

#elif AMREX_SPACEDIM == 2
                update(i,j,k,n) = update(i,j,k,n) + &
                     (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                     flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) ) / vol(i,j,k)

#else
                update(i,j,k,n) = update(i,j,k,n) + &
                     (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                     flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) + &
                     flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) ) / vol(i,j,k)
#endif

#if AMREX_SPACEDIM == 1
                if (do_hydro == 1) then
                   if (n == UMX) then
                      update(i,j,k,UMX) = update(i,j,k,UMX) - &
                           ( q1(i+1,j,k,GDPRES) - q1(i,j,k,GDPRES) ) / dx(1)
                   end if
                end if
#endif

#if AMREX_SPACEDIM == 2
                if (do_hydro == 1) then
                   if (n == UMX) then
                      ! add the pressure source term for axisymmetry
                      if (coord_type > 0) then
                         update(i,j,k,n) = update(i,j,k,n) - (q1(i+1,j,k,GDPRES) - q1(i,j,k,GDPRES))/ dx(1)
                      end if
                   end if
                end if
#endif

                if (n <= NSRC) then
                   update(i,j,k,n) = update(i,j,k,n) + srcU(i,j,k,n)
                end if
             enddo
          enddo
       enddo
    enddo

#ifdef SHOCK_VAR
    ! We'll update the shock data for future use in the burning step.
    ! For the update, we are starting from USHK == 0 (set at the
    ! beginning of the timestep) and we need to divide by dt since
    ! we'll be multiplying that for the update calculation.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             update(i,j,k,USHK) = shk(i,j,k) / dt
          enddo
       enddo
    enddo
#endif

#if AMREX_SPACEDIM == 3
#ifdef HYBRID_MOMENTUM
    call add_hybrid_advection_source(lo, hi, dt, &
                                     update, updt_lo, updt_hi, &
                                     q1, q1_lo, q1_hi, &
                                     q2, q2_lo, q2_hi, &
                                     q3, q3_lo, q3_hi)
#endif
#endif

  end subroutine ca_mol_consup

  subroutine ca_mol_diffusive_flux(lo, hi, &
                                   idir, &
                                   U, U_lo, U_hi, &
                                   cond, c_lo, c_hi, &
                                   flux, f_lo, f_hi, &
                                   dx) bind(C, name="ca_mol_diffusive_flux")

    use castro_error_module
    use meth_params_module, only : NVAR, UTEMP, UEINT, UEDEN
    use amrex_constants_module, only : HALF
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3)
    integer, intent(in) :: c_lo(3), c_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in), value :: idir

    real(rt), intent(in) :: U(U_lo(1):U_hi(1), U_lo(2):U_hi(2), U_lo(3):U_hi(3), NVAR)
    real(rt), intent(in) :: cond(c_lo(1):c_hi(1), c_lo(2):c_hi(2), c_lo(3):c_hi(3))
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(in) :: dx(3)

    integer :: i, j, k

    real(rt) :: diff_term, cond_int

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (idir == 1) then
                cond_int = HALF * (cond(i,j,k) + cond(i-1,j,k))
                diff_term = -cond_int * (U(i,j,k,UTEMP) - U(i-1,j,k,UTEMP))/dx(1)

             else if (idir == 2) then
                cond_int = HALF * (cond(i,j,k) + cond(i,j-1,k))
                diff_term = -cond_int * (U(i,j,k,UTEMP) - U(i,j-1,k,UTEMP))/dx(2)

             else
                cond_int = HALF * (cond(i,j,k) + cond(i,j,k-1))
                diff_term = -cond_int * (U(i,j,k,UTEMP) - U(i,j,k-1,UTEMP))/dx(3)
             end if

             flux(i,j,k,UEINT) = flux(i,j,k,UEINT) + diff_term
             flux(i,j,k,UEDEN) = flux(i,j,k,UEDEN) + diff_term

          end do
       end do
    end do

  end subroutine ca_mol_diffusive_flux

end module castro_mol_module
