#ifndef _HABEC_H_
#define _HABEC_H_

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

using namespace amrex;

namespace HABEC
{
    // habec is Hypre abec, where abec is the form of the linear equation
    // we are solving:
    // 
    // alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

    AMREX_INLINE
    void hbflx (Array4<Real> const flux,
                Array4<Real> const er,
                const Box& reg,
                int cdir, int bct, int bho, Real bcl,
                Array4<Real const> const bcval,
                Array4<int const> const mask,
                Array4<Real const> const b,
                Real beta, const Real* dx, int inhom)
    {
        Real h;

        bool x_left = false;
        bool x_right = false;
        bool y_left = false;
        bool y_right = false;
        bool z_left = false;
        bool z_right = false;

#if AMREX_SPACEDIM == 1
        if (cdir == 0) {
            h = dx[0];
            x_left = true;
        }
        else if (cdir == 1) {
            h = dx[0];
            x_right = true;
        }
#elif AMREX_SPACEDIM == 2
        if (cdir == 0) {
            h = dx[0];
            x_left = true;
        }
        else if (cdir == 2) {
            h = dx[0];
            x_right = true;
        }
        else if (cdir == 1) {
            h = dx[1];
            y_left = true;
        }
        else if (cdir == 3) {
            h = dx[1];
            y_right = true;
        }
#else
        if (cdir == 0) {
            h = dx[0];
            x_left = true;
        }
        else if (cdir == 3) {
            h = dx[0];
            x_right = true;
        }
        else if (cdir == 1) {
            h = dx[1];
            y_left = true;
        }
        else if (cdir == 4) {
            h = dx[1];
            y_right = true;
        }
        else if (cdir == 2) {
            h = dx[2];
            z_left = true;
        }
        else if (cdir == 5) {
            h = dx[2];
            z_right = true;
        }
#endif

        Real bfv, bfm, bfm2;

        if (bct == LO_DIRICHLET) {
            if (bho >= 1) {
                Real h2 = 0.5e0_rt * h;
                Real th2 = 3.e0_rt * h2;
                bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2));
                bfm = (beta / h) * (th2 - bcl) / (bcl + h2);
                bfm2 = (beta / h) * (bcl - h2) / (bcl + th2);
            }
            else {
                bfv = beta / (0.5e0_rt * h + bcl);
                bfm = bfv;
            }
        }
        else {
            amrex::Error("hbflx: unsupported boundary type");
        }

        if (inhom == 0) {
            bfv = 0.e0_rt;
        }

        int reg_ilo = reg.loVect3d()[0];
        int reg_ihi = reg.hiVect3d()[0];
        int reg_jlo = reg.loVect3d()[1];
        int reg_jhi = reg.hiVect3d()[1];
        int reg_klo = reg.loVect3d()[2];
        int reg_khi = reg.hiVect3d()[2];

        if (x_left) {
            int i = reg_ilo;
            for (int k = reg_klo; k <= reg_khi; ++k) {
                for (int j = reg_jlo; j <= reg_jhi; ++j) {
                    if (mask(i-1,j,k) > 0) {
                        flux(i,j,k) = b(i,j,k) * (bfv * bcval(i-1,j,k) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i,j,k) = flux(i,j,k) - b(i,j,k) * bfm2 * er(i+1,j,k);
                        }
                    }
                }
            }
        }
        else if (x_right) {
            int i = reg_ihi;
            for (int k = reg_klo; k <= reg_khi; ++k) {
                for (int j = reg_jlo; j <= reg_jhi; ++j) {
                    if (mask(i+1,j,k) > 0) {
                        flux(i+1,j,k) = -b(i+1,j,k) * (bfv * bcval(i+1,j,k) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i+1,j,k) = flux(i+1,j,k) + b(i+1,j,k) * bfm2 * er(i-1,j,k);
                        }
                    }
                }
            }
        }
        else if (y_left) {
            int j = reg_jlo;
            for (int k = reg_klo; k <= reg_khi; ++k) {
                for (int i = reg_ilo; i <= reg_ihi; ++i) {
                    if (mask(i,j-1,k) > 0) {
                        flux(i,j,k) = b(i,j,k) * (bfv * bcval(i,j-1,k) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i,j,k) = flux(i,j,k) - b(i,j,k) * bfm2 * er(i,j+1,k);
                        }
                    }
                }
            }
        }
        else if (y_right) {
            int j = reg_jhi;
            for (int k = reg_klo; k <= reg_khi; ++k) {
                for (int i = reg_ilo; i <= reg_ihi; ++i) {
                    if (mask(i,j+1,k) > 0) {
                        flux(i,j+1,k) = -b(i,j+1,k) * (bfv * bcval(i,j+1,k) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i,j+1,k) = flux(i,j+1,k) + b(i,j+1,k) * bfm2 * er(i,j-1,k);
                        }
                    }
                }
            }
        }
        else if (z_left) {
            int k = reg_klo;
            for (int j = reg_jlo; j <= reg_jhi; ++j) {
                for (int i = reg_ilo; i <= reg_ihi; ++i) {
                    if (mask(i,j,k-1) > 0) {
                        flux(i,j,k) = b(i,j,k) * (bfv * bcval(i,j,k-1) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i,j,k) = flux(i,j,k) - b(i,j,k) * bfm2 * er(i,j,k+1);
                        }
                    }
                }
            }
        }
        else if (z_right) {
            int k = reg_khi;
            for (int j = reg_jlo; j <= reg_jhi; ++j) {
                for (int i = reg_ilo; i <= reg_ihi; ++i) {
                    if (mask(i,j,k+1) > 0) {
                        flux(i,j,k+1) = -b(i,j,k+1) * (bfv * bcval(i,j,k+1) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i,j,k+1) = flux(i,j,k+1) + b(i,j,k+1) * bfm2 * er(i,j,k-1);
                        }
                    }
                }
            }
        }
        else {
            std::cout << "hbflx: impossible face orientation" << std::endl;
        }
    }

    AMREX_INLINE
    void hbflx3 (Array4<Real> const flux,
                 Array4<Real const> const er,
                 const Box& reg,
                 int cdir, int bctype,
                 Array4<int const> const tf,
                 int bho, Real bcl,
                 Array4<Real const> const bcval,
                 Array4<int const> const mask,
                 Array4<Real const> const b,
                 Real beta, const Real* const dx, Real c,
                 const Real* const r, int inhom,
                 Array4<Real const> const spa)
    {
        int reg_ilo = reg.loVect3d()[0];
        int reg_ihi = reg.hiVect3d()[0];
        int reg_jlo = reg.loVect3d()[1];
        int reg_jhi = reg.hiVect3d()[1];
        int reg_klo = reg.loVect3d()[2];
        int reg_khi = reg.hiVect3d()[2];

        Real h;

        bool x_left = false;
        bool x_right = false;
        bool y_left = false;
        bool y_right = false;
        bool z_left = false;
        bool z_right = false;

#if AMREX_SPACEDIM == 1
        if (cdir == 0) {
            h = dx[0];
            x_left = true;
        }
        else if (cdir == 1) {
            h = dx[0];
            x_right = true;
        }

        // r is passed as an array but actually has only one element, which is
        // appropriate for the face we are doing here.
        Real r0 = r[0];
#elif AMREX_SPACEDIM == 2
        if (cdir == 0) {
            h = dx[0];
            x_left = true;
        }
        else if (cdir == 2) {
            h = dx[0];
            x_right = true;
        }
        else if (cdir == 1) {
            h = dx[1];
            y_left = true;
        }
        else if (cdir == 3) {
            h = dx[1];
            y_right = true;
        }

        // For the left and right faces, r is constant and the array actually has
        // only one element.  This following "r[0]" looks wrong, but what it
        // really means is that the array dimensions are meaningless for this case.
        Real r0 = r[0];
#else
        if (cdir == 0) {
            h = dx[0];
            x_left = true;
        }
        else if (cdir == 3) {
            h = dx[0];
            x_right = true;
        }
        else if (cdir == 1) {
            h = dx[1];
            y_left = true;
        }
        else if (cdir == 4) {
            h = dx[1];
            y_right = true;
        }
        else if (cdir == 2) {
            h = dx[2];
            z_left = true;
        }
        else if (cdir == 5) {
            h = dx[2];
            z_right = true;
        }

        Real r0 = 1.0_rt;
#endif

        Real bfv, bfm, bfm2, h2, th2;
        int bct;

        if (x_left) {
            int i = reg_ilo;
            for (int k = reg_klo; k <= reg_khi; ++k) {
                for (int j = reg_jlo; j <= reg_jhi; ++j) {
                    if (mask(i-1,j,k) > 0) {
                        if (bctype == -1) {
                            bct = tf(i-1,j,k);
                        }
                        else {
                            bct = bctype;
                        }
                        if (bct == LO_DIRICHLET) {
                            if (bho >= 1) {
                                h2 = 0.5e0_rt * h;
                                th2 = 3.e0_rt * h2;
                                bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k);
                                bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k);
                                bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k);
                            }
                            else {
                                bfv = beta / (0.5e0_rt * h + bcl) * b(i,j,k);
                                bfm = bfv;
                            }
                        }
                        else if (bct == LO_NEUMANN) {
                            bfv  = beta * r0;
                            bfm  = 0.e0_rt;
                            bfm2 = 0.e0_rt;
                        }
                        else if (bct == LO_MARSHAK) {
                            bfv = 2.e0_rt * beta * r0;
                            if (bho >= 1) {
                                bfm  =  0.75e0_rt * beta * c;
                                bfm2 = -0.25e0_rt * beta * c;
                            }
                            else {
                                bfm = 0.5e0_rt * beta * c;
                            }
                        }
                        else if (bct == LO_SANCHEZ_POMRANING) {
                            bfv = 2.e0_rt * beta * r0;
                            if (bho >= 1) {
                                bfm  =  3.0e0_rt * spa(i,j,k) * beta * c;
                                bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c;
                            }
                            else {
                                bfm = 2.0e0_rt * spa(i,j,k) * beta * c;
                            }
                        }
                        else {
                            amrex::Error("hbflx3: unsupported boundary type");
                        }
                        if (inhom == 0) {
                            bfv = 0.e0_rt;
                        }
                        flux(i,j,k) = (bfv * bcval(i-1,j,k) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i,j,k) = flux(i,j,k) - bfm2 * er(i+1,j,k);
                        }
                    }
                }
            }
        }
        else if (x_right) {
            int i = reg_ihi;
            for (int k = reg_klo; k <= reg_khi; ++k) {
                for (int j = reg_jlo; j <= reg_jhi; ++j) {
                    if (mask(i+1,j,k) > 0) {
                        if (bctype == -1) {
                            bct = tf(i+1,j,k);
                        }
                        else {
                            bct = bctype;
                        }
                        if (bct == LO_DIRICHLET) {
                            if (bho >= 1) {
                                h2 = 0.5e0_rt * h;
                                th2 = 3.e0_rt * h2;
                                bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i+1,j,k);
                                bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i+1,j,k);
                                bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i+1,j,k);
                            }
                            else {
                                bfv = beta / (0.5e0_rt * h + bcl) * b(i+1,j,k);
                                bfm = bfv;
                            }
                        }
                        else if (bct == LO_NEUMANN) {
                            bfv  = beta * r0;
                            bfm  = 0.e0_rt;
                            bfm2 = 0.e0_rt;
                        }
                        else if (bct == LO_MARSHAK) {
                            bfv = 2.e0_rt * beta * r0;
                            if (bho >= 1) {
                                bfm  =  0.75e0_rt * beta * c;
                                bfm2 = -0.25e0_rt * beta * c;
                            }
                            else {
                                bfm = 0.5e0_rt * beta * c;
                            }
                        }
                        else if (bct == LO_SANCHEZ_POMRANING) {
                            bfv = 2.e0_rt * beta * r0;
                            if (bho >= 1) {
                                bfm  =  3.0e0_rt * spa(i,j,k) * beta * c;
                                bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c;
                            }
                            else {
                                bfm = 2.0e0_rt * spa(i,j,k) * beta * c;
                            }
                        }
                        else {
                            amrex::Error("hbflx3: unsupported boundary type");
                        }
                        if (inhom == 0) {
                            bfv = 0.e0_rt;
                        }
                        flux(i+1,j,k) = -(bfv * bcval(i+1,j,k) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i+1,j,k) = flux(i+1,j,k) + bfm2 * er(i-1,j,k);
                        }
                    }
                }
            }
        }
        else if (y_left) {
            int j = reg_jlo;
            for (int k = reg_klo; k <= reg_khi; ++k) {
                for (int i = reg_ilo; i <= reg_ihi; ++i) {
                    if (mask(i,j-1,k) > 0) {
                        if (bctype == -1) {
                            bct = tf(i,j-1,k);
                        }
                        else {
                            bct = bctype;
                        }
                        if (bct == LO_DIRICHLET) {
                            if (bho >= 1) {
                                h2 = 0.5e0_rt * h;
                                th2 = 3.e0_rt * h2;
                                bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k);
                                bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k);
                                bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k);
                            }
                            else {
                                bfv = beta / (0.5e0_rt * h + bcl) * b(i,j,k);
                                bfm = bfv;
                            }
                        }
                        else if (bct == LO_NEUMANN) {
                            bfv  = beta * r[i - reg_ilo];
                            bfm  = 0.e0_rt;
                            bfm2 = 0.e0_rt;
                        }
                        else if (bct == LO_MARSHAK) {
                            bfv = 2.e0_rt * beta * r[i - reg_ilo];
                            if (bho >= 1) {
                                bfm  =  0.75e0_rt * beta * c;
                                bfm2 = -0.25e0_rt * beta * c;
                            }
                            else {
                                bfm = 0.5e0_rt * beta * c;
                            }
                        }
                        else if (bct == LO_SANCHEZ_POMRANING) {
                            bfv = 2.e0_rt * beta * r[i - reg_ilo];
                            if (bho >= 1) {
                                bfm  =  3.0e0_rt * spa(i,j,k) * beta * c;
                                bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c;
                            }
                            else {
                                bfm = 2.0e0_rt * spa(i,j,k) * beta * c;
                            }
                        }
                        else {
                            amrex::Error("hbflx3: unsupported boundary type");
                        }
                        if (inhom == 0) {
                            bfv = 0.e0_rt;
                        }
                        flux(i,j,k) = (bfv * bcval(i,j-1,k) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i,j,k) = flux(i,j,k) - bfm2 * er(i,j+1,k);
                        }
                    }
                }
            }
        }
        else if (y_right) {
            int j = reg_jhi;
            for (int k = reg_klo; k <= reg_khi; ++k) {
                for (int i = reg_ilo; i <= reg_ihi; ++i) {
                    if (mask(i,j+1,k) > 0) {
                        if (bctype == -1) {
                            bct = tf(i,j+1,k);
                        }
                        else {
                            bct = bctype;
                        }
                        if (bct == LO_DIRICHLET) {
                            if (bho >= 1) {
                                h2 = 0.5e0_rt * h;
                                th2 = 3.e0_rt * h2;
                                bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j+1,k);
                                bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j+1,k);
                                bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j+1,k);
                            }
                            else {
                                bfv = beta / (0.5e0_rt * h + bcl) * b(i,j+1,k);
                                bfm = bfv;
                            }
                        }
                        else if (bct == LO_NEUMANN) {
                            bfv  = beta * r[i - reg_ilo];
                            bfm  = 0.e0_rt;
                            bfm2 = 0.e0_rt;
                        }
                        else if (bct == LO_MARSHAK) {
                            bfv = 2.e0_rt * beta * r[i - reg_ilo];
                            if (bho >= 1) {
                                bfm  =  0.75e0_rt * beta * c;
                                bfm2 = -0.25e0_rt * beta * c;
                            }
                            else {
                                bfm = 0.5e0_rt * beta * c;
                            }
                        }
                        else if (bct == LO_SANCHEZ_POMRANING) {
                            bfv = 2.e0_rt * beta * r[i - reg_ilo];
                            if (bho >= 1) {
                                bfm  =  3.0e0_rt * spa(i,j,k) * beta * c;
                                bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c;
                            }
                            else {
                                bfm = 2.0e0_rt * spa(i,j,k) * beta * c;
                            }
                        }
                        else {
                            amrex::Error("hbflx3: unsupported boundary type");
                        }
                        if (inhom == 0) {
                            bfv = 0.e0_rt;
                        }
                        flux(i,j+1,k) = -(bfv * bcval(i,j+1,k) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i,j+1,k) = flux(i,j+1,k) + bfm2 * er(i,j-1,k);
                        }
                    }
                }
            }
        }
        else if (z_left) {
            int k = reg_klo;
            for (int j = reg_jlo; j <= reg_jhi; ++j) {
                for (int i = reg_ilo; i <= reg_ihi; ++i) {
                    if (mask(i,j,k-1) > 0) {
                        if (bctype == -1) {
                            bct = tf(i,j,k-1);
                        }
                        else {
                            bct = bctype;
                        }
                        if (bct == LO_DIRICHLET) {
                            if (bho >= 1) {
                                h2 = 0.5e0_rt * h;
                                th2 = 3.e0_rt * h2;
                                bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k);
                                bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k);
                                bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k);
                            }
                            else {
                                bfv = beta / (0.5e0_rt * h + bcl) * b(i,j,k);
                                bfm = bfv;
                            }
                        }
                        else if (bct == LO_NEUMANN) {
                            bfv  = beta;
                            bfm  = 0.e0_rt;
                            bfm2 = 0.e0_rt;
                        }
                        else if (bct == LO_MARSHAK) {
                            bfv = 2.e0_rt * beta;
                            if (bho >= 1) {
                                bfm  =  0.75e0_rt * beta * c;
                                bfm2 = -0.25e0_rt * beta * c;
                            }
                            else {
                                bfm = 0.5e0_rt * beta * c;
                            }
                        }
                        else if (bct == LO_SANCHEZ_POMRANING) {
                            bfv = 2.e0_rt * beta;
                            if (bho >= 1) {
                                bfm  =  3.0e0_rt * spa(i,j,k) * beta * c;
                                bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c;
                            }
                            else {
                                bfm = 2.0e0_rt * spa(i,j,k) * beta * c;
                            }
                        }
                        else {
                            amrex::Error("hbflx3: unsupported boundary type");
                        }
                        if (inhom == 0) {
                            bfv = 0.e0_rt;
                        }
                        flux(i,j,k) = (bfv * bcval(i,j,k-1) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i,j,k) = flux(i,j,k) - bfm2 * er(i,j,k+1);
                        }
                    }
                }
            }
        }
        else if (z_right) {
            int k = reg_khi;
            for (int j = reg_jlo; j <= reg_jhi; ++j) {
                for (int i = reg_ilo; i <= reg_ihi; ++i) {
                    if (mask(i,j,k+1) > 0) {
                        if (bctype == -1) {
                            bct = tf(i,j,k+1);
                        }
                        else {
                            bct = bctype;
                        }
                        if (bct == LO_DIRICHLET) {
                            if (bho >= 1) {
                                h2 = 0.5e0_rt * h;
                                th2 = 3.e0_rt * h2;
                                bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j,k+1);
                                bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j,k+1);
                                bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j,k+1);
                            }
                            else {
                                bfv = beta / (0.5e0_rt * h + bcl) * b(i,j,k+1);
                                bfm = bfv;
                            }
                        }
                        else if (bct == LO_NEUMANN) {
                            bfv  = beta;
                            bfm  = 0.e0_rt;
                            bfm2 = 0.e0_rt;
                        }
                        else if (bct == LO_MARSHAK) {
                            bfv = 2.e0_rt * beta;
                            if (bho >= 1) {
                                bfm  =  0.75e0_rt * beta * c;
                                bfm2 = -0.25e0_rt * beta * c;
                            }
                            else {
                                bfm = 0.5e0_rt * beta * c;
                            }
                        }
                        else if (bct == LO_SANCHEZ_POMRANING) {
                            bfv = 2.e0_rt * beta;
                            if (bho >= 1) {
                                bfm  =  3.0e0_rt * spa(i,j,k) * beta * c;
                                bfm2 = -1.0e0_rt * spa(i,j,k) * beta * c;
                            }
                            else {
                                bfm = 2.0e0_rt * spa(i,j,k) * beta * c;
                            }
                        }
                        else {
                            amrex::Error("hbflx3: unsupported boundary type");
                        }
                        if (inhom == 0) {
                            bfv = 0.e0_rt;
                        }
                        flux(i,j,k+1) = -(bfv * bcval(i,j,k+1) - bfm * er(i,j,k));
                        if (bho >= 1) {
                            flux(i,j,k+1) = flux(i,j,k+1) + bfm2 * er(i,j,k-1);
                        }
                    }
                }
            }
        }
        else {
            std::cout << "hbflx3: impossible face orientation" << std::endl;
        }
    }

} // namespace HABEC

#endif
