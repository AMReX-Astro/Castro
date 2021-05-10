
#include <Radiation.H>

#include <RAD_F.H>

using namespace amrex;

void Radiation::save_lambda_in_plotvar(int level, const Array<MultiFab,BL_SPACEDIM>& lambda)
{
    int nlambda = lambda[0].nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*plotvar[level],true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto lamx = lambda[0][mfi].array();
#if AMREX_SPACEDIM >= 2
        auto lamy = lambda[1][mfi].array();
#endif
#if AMREX_SPACEDIM == 3
        auto lamz = lambda[2][mfi].array();
#endif

        auto lamc = (*plotvar[level])[mfi].array();

        int dcomp = icomp_lambda;

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            for (int n = 0; n < nlambda; ++n) {
#if AMREX_SPACEDIM == 1
                lamc(i,j,k,dcomp+n) = (lamx(i,j,k,n) + lamx(i+1,j,k,n)) * 0.5_rt;
#elif AMREX_SPACEDIM == 2
                lamc(i,j,k,dcomp+n) = (lamx(i,j,k,n) + lamx(i+1,j,k,n) +
                                       lamy(i,j,k,n) + lamy(i,j+1,k,n)) * 0.25_rt;
#else
                lamc(i,j,k,dcomp+n) = (lamx(i,j,k,n) + lamx(i+1,j,k,n) +
                                       lamy(i,j,k,n) + lamy(i,j+1,k,n) +
                                       lamz(i,j,k,n) + lamz(i,j,k+1,n)) * (1.0_rt / 6.0_rt);
#endif
            }
        });
    }
}

void Radiation::save_lab_Er_in_plotvar(int level, const MultiFab& Snew, 
                                       const MultiFab& Ecom, const MultiFab& F, int iflx)
{
    int nflx = F.nComp();

    int ifx = iflx;

#if AMREX_SPACEDIM >= 2
    int ify = iflx + NGROUPS;
#else
    int ify = 0;
#endif

#if AMREX_SPACEDIM == 3
    int ifz = iflx + NGROUPS * NGROUPS;
#else
    int ifz = 0;
#endif

    int ier = icomp_lab_Er;

    Real c2 = 1.e0_rt / (Radiation::clight * Radiation::clight);

    GpuArray<Real, NGROUPS> dlognu = {0.0};

    if (NGROUPS > 1) {
        ca_get_dlognu(dlognu.begin());
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*plotvar[level],true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto Snew_arr = Snew[mfi].array();
        auto Ecom_arr = Ecom[mfi].array();
        auto F_arr = F[mfi].array();
        auto Elab_arr = (*plotvar[level])[mfi].array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            Real rhoInv = 1.0_rt / Snew_arr(i,j,k,URHO);

            Real vxc2 = Snew_arr(i,j,k,UMX) * rhoInv * c2;

#if AMREX_SPACEDIM >= 2
            Real vyc2 = Snew_arr(i,j,k,UMY) * rhoInv * c2;
#else
            Real vyc2 = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
            Real vzc2 = Snew_arr(i,j,k,UMZ) * rhoInv * c2;
#else
            Real vzc2 = 0.0_rt;
#endif

            for (int g = 0; g < NGROUPS; ++g) {
                Elab_arr(i,j,k,g+ier) = Ecom_arr(i,j,k,g) + 2.0_rt * (vxc2 * F_arr(i,j,k,ifx+g) +
                                                                      vyc2 * F_arr(i,j,k,ify+g) +
                                                                      vzc2 * F_arr(i,j,k,ifz+g));
            }

            if (NGROUPS > 1) {
                Array1D<Real, -1, NGROUPS> nufnux, nufnuy, nufnuz;

                for (int g = 0; g < NGROUPS; ++g) {
                    nufnux(g) = F_arr(i,j,k,ifx+g) / dlognu[g];
                    nufnuy(g) = F_arr(i,j,k,ify+g) / dlognu[g];
                    nufnuz(g) = F_arr(i,j,k,ifz+g) / dlognu[g];
                }

                nufnux(-1) = -nufnux(0);
                nufnuy(-1) = -nufnuy(0);
                nufnuz(-1) = -nufnuz(0);
                nufnux(NGROUPS) = -nufnux(NGROUPS-1);
                nufnuy(NGROUPS) = -nufnuy(NGROUPS-1);
                nufnuz(NGROUPS) = -nufnuz(NGROUPS-1);

                for (int g = 0; g < NGROUPS; ++g) {
                    Elab_arr(i,j,k,g+ier) -= vxc2 * 0.5_rt * (nufnux(g+1) - nufnux(g-1)) +
                                             vyc2 * 0.5_rt * (nufnuy(g+1) - nufnuy(g-1)) +
                                             vzc2 * 0.5_rt * (nufnuz(g+1) - nufnuz(g-1));
                }
            }
        });
    }
}

void Radiation::save_lab_flux_in_plotvar(int level, const MultiFab& Snew, 
                                         const Array<MultiFab,BL_SPACEDIM>& lambda,
                                         const MultiFab& Er, const MultiFab& Fr, int iflx)
{
    const Real flag = 1.0;  // comovinng --> lab

    int nflx = Fr.nComp();
    int nlambda = lambda[0].nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox f;
        for (MFIter mfi(*plotvar[level],true); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            f.resize(bx, nGroups);

            BL_FORT_PROC_CALL(CA_COMPUTE_FCC, ca_compute_fcc)
                (bx.loVect(), bx.hiVect(),
                 D_DECL(BL_TO_FORTRAN(lambda[0][mfi]),
                        BL_TO_FORTRAN(lambda[1][mfi]),
                        BL_TO_FORTRAN(lambda[2][mfi])), nlambda,
                 BL_TO_FORTRAN(f));

            BL_FORT_PROC_CALL(CA_TRANSFORM_FLUX, ca_transform_flux)
                (bx.loVect(), bx.hiVect(), flag,
                 BL_TO_FORTRAN(Snew[mfi]),
                 BL_TO_FORTRAN(f),
                 BL_TO_FORTRAN(Er[mfi]),
                 BL_TO_FORTRAN(Fr[mfi]), iflx, nflx,
                 BL_TO_FORTRAN((*plotvar[level])[mfi]), icomp_lab_Fr, nplotvar);
        }
    }
}

void Radiation::save_com_flux_in_plotvar(int level, const MultiFab& Snew, 
                                         const Array<MultiFab,BL_SPACEDIM>& lambda,
                                         const MultiFab& Er, const MultiFab& Fr, int iflx)
{
    const Real flag = -1.0;  // lab --> comoving

    int nflx = Fr.nComp();
    int nlambda = lambda[0].nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox f;
        for (MFIter mfi(*plotvar[level],true); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            f.resize(bx, nGroups);

            BL_FORT_PROC_CALL(CA_COMPUTE_FCC, ca_compute_fcc)
                (bx.loVect(), bx.hiVect(),
                 D_DECL(BL_TO_FORTRAN(lambda[0][mfi]),
                        BL_TO_FORTRAN(lambda[1][mfi]),
                        BL_TO_FORTRAN(lambda[2][mfi])), nlambda,
                 BL_TO_FORTRAN(f));

            BL_FORT_PROC_CALL(CA_TRANSFORM_FLUX, ca_transform_flux)
                (bx.loVect(), bx.hiVect(), flag,
                 BL_TO_FORTRAN(Snew[mfi]),
                 BL_TO_FORTRAN(f),
                 BL_TO_FORTRAN(Er[mfi]),
                 BL_TO_FORTRAN(Fr[mfi]), iflx, nflx,
                 BL_TO_FORTRAN((*plotvar[level])[mfi]), icomp_com_Fr, nplotvar);
        }
    }
}

