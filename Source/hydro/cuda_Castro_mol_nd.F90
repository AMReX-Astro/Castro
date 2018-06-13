! advection routines in support of method of lines integration

module mol_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  AMREX_DEVICE subroutine ca_construct_flux(lo, hi, domlo, domhi, dx, dt, idir, &
                                            uin, uin_lo, uin_hi, &
                                            div, div_lo, div_hi, &
                                            qaux, qa_lo, qa_hi, &
                                            qm, qm_lo, qm_hi, &
                                            qp, qp_lo, qp_hi, &
                                            qint, qe_lo, qe_hi, &
                                            flux, f_lo, f_hi, &
                                            area, a_lo, a_hi) &
                                            bind(c,name='ca_construct_flux')

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, NGDNV, NQAUX, NQ
    use advection_util_module, only: apply_av, normalize_species_fluxes, scale_flux
    use riemann_module, only: cmpflx

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ), value :: idir
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)
    integer,  intent(in   ) :: qe_lo(3), qe_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1), div_lo(2):div_hi(2), div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ,3)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ,3)
    real(rt), intent(inout) :: qint(qe_lo(1):qe_hi(1), qe_lo(2):qe_hi(2), qe_lo(3):qe_hi(3), NGDNV)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt

    call cmpflx(lo, hi, domlo, domhi, idir, qm, qm_lo, qm_hi, qp, qp_lo, qp_hi, &
                qint, qe_lo, qe_hi, flux, f_lo, f_hi, qaux, qa_lo, qa_hi)
    call apply_av(lo, hi, idir, dx, div, div_lo, div_hi, uin, uin_lo, uin_hi, flux, f_lo, f_hi)
    call normalize_species_fluxes(lo, hi, flux, f_lo, f_hi)
    call scale_flux(lo, hi, flux, f_lo, f_hi, area, a_lo, a_hi, dt)

  end subroutine ca_construct_flux

end module mol_module
