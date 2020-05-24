#ifndef CASTRO_MHD_F_H
#define CASTRO_MHD_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void mhd_flatten
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn));

  void plm
    (const int lo[], const int hi[],
     const int idir,
     const BL_FORT_FAB_ARG_3D(s),
     const BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(flatn),
     const BL_FORT_FAB_ARG_3D(bx),
     const BL_FORT_FAB_ARG_3D(by),
     const BL_FORT_FAB_ARG_3D(bz),
     BL_FORT_FAB_ARG_3D(Ip),
     BL_FORT_FAB_ARG_3D(Im),
     const BL_FORT_FAB_ARG_3D(srcq),
     const amrex::Real dx[], const amrex::Real dt);

  void electric_edge_x(const int lo[], const int hi[],
                       const BL_FORT_FAB_ARG_3D(q),
                       BL_FORT_FAB_ARG_3D(E),
                       BL_FORT_FAB_ARG_3D(flxy),
                       BL_FORT_FAB_ARG_3D(flxz));

  void electric_edge_y(const int lo[], const int hi[],
                       const BL_FORT_FAB_ARG_3D(q),
                       BL_FORT_FAB_ARG_3D(E),
                       BL_FORT_FAB_ARG_3D(flxx),
                       BL_FORT_FAB_ARG_3D(flxz));

  void electric_edge_z(const int lo[], const int hi[],
                       const BL_FORT_FAB_ARG_3D(q),
                       BL_FORT_FAB_ARG_3D(E),
                       BL_FORT_FAB_ARG_3D(flxx),
                       BL_FORT_FAB_ARG_3D(flxy));

  void corner_couple(const int lo[], const int hi[],
                     BL_FORT_FAB_ARG_3D(qr_out),
                     BL_FORT_FAB_ARG_3D(ql_out),
                     const BL_FORT_FAB_ARG_3D(ur),
                     const BL_FORT_FAB_ARG_3D(ul),
                     const BL_FORT_FAB_ARG_3D(flxd2),
                     const BL_FORT_FAB_ARG_3D(Ed1),
                     const BL_FORT_FAB_ARG_3D(Ed3),
                     const int d1, const int d2, const int d3,
                     const Real dx, const Real dt);

  void half_step(const int lo[], const int hi[],
                 BL_FORT_FAB_ARG_3D(qr_out),
                 BL_FORT_FAB_ARG_3D(ql_out),
                 const BL_FORT_FAB_ARG_3D(ur),
                 const BL_FORT_FAB_ARG_3D(ul),
                 const BL_FORT_FAB_ARG_3D(flxd1),
                 const BL_FORT_FAB_ARG_3D(flxd2),
                 const BL_FORT_FAB_ARG_3D(Ed),
                 const BL_FORT_FAB_ARG_3D(Ed1),
                 const BL_FORT_FAB_ARG_3D(Ed2),
                 const int d, const int d1, const int d2,
                 const Real dx, const Real dt);

  void hlld(const int lo[], const int hi[],
            const BL_FORT_FAB_ARG_3D(qleft),
            const BL_FORT_FAB_ARG_3D(qright),
            BL_FORT_FAB_ARG_3D(flx),
            const int dir);


#ifdef __cplusplus
}
#endif

#endif

