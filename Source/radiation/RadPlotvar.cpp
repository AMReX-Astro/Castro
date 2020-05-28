
#include <Radiation.H>

#include <RAD_F.H>

using namespace amrex;

void Radiation::save_lambda_in_plotvar(
    int level, const Array<MultiFab, BL_SPACEDIM>& lambda) {
    int nlambda = lambda[0].nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*plotvar[level], true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        int scomp = 0;
        ca_face2center(
            bx.loVect(), bx.hiVect(), scomp, icomp_lambda, nlambda, nlambda,
            nplotvar,
            D_DECL(BL_TO_FORTRAN(lambda[0][mfi]), BL_TO_FORTRAN(lambda[1][mfi]),
                   BL_TO_FORTRAN(lambda[2][mfi])),
            BL_TO_FORTRAN((*plotvar[level])[mfi]));
    }
}

void Radiation::save_lab_Er_in_plotvar(int level, const MultiFab& Snew,
                                       const MultiFab& Ecom, const MultiFab& F,
                                       int iflx) {
    int nflx = F.nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*plotvar[level], true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        BL_FORT_PROC_CALL(CA_ER_COM2LAB, ca_er_com2lab)
        (bx.loVect(), bx.hiVect(), BL_TO_FORTRAN(Snew[mfi]),
         BL_TO_FORTRAN(Ecom[mfi]), BL_TO_FORTRAN(F[mfi]), iflx, nflx,
         BL_TO_FORTRAN((*plotvar[level])[mfi]), icomp_lab_Er, nplotvar);
    }
}

void Radiation::save_lab_flux_in_plotvar(
    int level, const MultiFab& Snew, const Array<MultiFab, BL_SPACEDIM>& lambda,
    const MultiFab& Er, const MultiFab& Fr, int iflx) {
    const Real flag = 1.0;  // comovinng --> lab

    int nflx = Fr.nComp();
    int nlambda = lambda[0].nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox f;
        for (MFIter mfi(*plotvar[level], true); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            f.resize(bx, nGroups);

            BL_FORT_PROC_CALL(CA_COMPUTE_FCC, ca_compute_fcc)
            (bx.loVect(), bx.hiVect(),
             D_DECL(BL_TO_FORTRAN(lambda[0][mfi]),
                    BL_TO_FORTRAN(lambda[1][mfi]),
                    BL_TO_FORTRAN(lambda[2][mfi])),
             nlambda, BL_TO_FORTRAN(f));

            BL_FORT_PROC_CALL(CA_TRANSFORM_FLUX, ca_transform_flux)
            (bx.loVect(), bx.hiVect(), flag, BL_TO_FORTRAN(Snew[mfi]),
             BL_TO_FORTRAN(f), BL_TO_FORTRAN(Er[mfi]), BL_TO_FORTRAN(Fr[mfi]),
             iflx, nflx, BL_TO_FORTRAN((*plotvar[level])[mfi]), icomp_lab_Fr,
             nplotvar);
        }
    }
}

void Radiation::save_com_flux_in_plotvar(
    int level, const MultiFab& Snew, const Array<MultiFab, BL_SPACEDIM>& lambda,
    const MultiFab& Er, const MultiFab& Fr, int iflx) {
    const Real flag = -1.0;  // lab --> comoving

    int nflx = Fr.nComp();
    int nlambda = lambda[0].nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox f;
        for (MFIter mfi(*plotvar[level], true); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            f.resize(bx, nGroups);

            BL_FORT_PROC_CALL(CA_COMPUTE_FCC, ca_compute_fcc)
            (bx.loVect(), bx.hiVect(),
             D_DECL(BL_TO_FORTRAN(lambda[0][mfi]),
                    BL_TO_FORTRAN(lambda[1][mfi]),
                    BL_TO_FORTRAN(lambda[2][mfi])),
             nlambda, BL_TO_FORTRAN(f));

            BL_FORT_PROC_CALL(CA_TRANSFORM_FLUX, ca_transform_flux)
            (bx.loVect(), bx.hiVect(), flag, BL_TO_FORTRAN(Snew[mfi]),
             BL_TO_FORTRAN(f), BL_TO_FORTRAN(Er[mfi]), BL_TO_FORTRAN(Fr[mfi]),
             iflx, nflx, BL_TO_FORTRAN((*plotvar[level])[mfi]), icomp_com_Fr,
             nplotvar);
        }
    }
}
