void
Castro::advance_mhd(Real time, Real dt)
{




  FArrayBox flux[AMREX_SPACEDIM], E[AMREX_SPACEDIM];
  FArrayBox cs[AMREX_SPACEDIM];

  FArrayBox bcc;
  FArrayBox q;
  FArrayBox srcQ;

  FArrayBox qm;
  FArrayBox qp;

  FArrayBox flxx;
  FArrayBox flxy;
  FArrayBox flxz;

  FArrayBox Extmp;
  FArrayBox Eytmp;
  FArrayBox Eztmp;

  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
  {

    const Box& bx = mfi.tilebox();

    const Box& obx = amrex::grow(bx, 1);

    // box with NUM_GROW ghost cells for PPM stuff
    const Box& bx_gc = amrex::grow(bx, NUM_GROW);

    const int* lo = bx.loVect();
    const int* hi = bx.hiVect();

    FArrayBox &statein  = Sborder[mfi];
    FArrayBox &stateout = S_new[mfi];

    FArrayBox &source_in  = sources_for_hydro[mfi];
    FArrayBox &source_out = hydro_source[mfi]; 

    FArrayBox& Bx  = Bx_old_tmp[mfi];
    FArrayBox& By  = By_old_tmp[mfi]; 
    FArrayBox& Bz  = Bz_old_tmp[mfi];

    FArrayBox& Bxout = Bx_new[mfi];
    FArrayBox& Byout = By_new[mfi];
    FArrayBox& Bzout = Bz_new[mfi];


    // allocate fabs for fluxes
    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
      const Box& bxtmp = amrex::surroundingNodes(bx, i);
      flux[i].resize(bxtmp, NUM_STATE);
      E[i].resize(bxtmp, NUM_STATE);
    }

    // Calculate primitives based on conservatives
    bcc.resize(bc_gc, 3);
    q.resize(bc_gc, NQ);
    srcQ.resize(obx, NQSRC);

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
      cs[idir].resize(bc_gc, 1);
    }

    const int* lo_gc = bx_gc.loVect();
    const int* hi_gc = bx_gc.hiVect();

    ctoprim_mhd(lo_gc, hi_gc,
                BL_TO_FORTRAN_ANYD(statein),
                BL_TO_FORTRAN_ANYD(bcc),
                BL_TO_FORTRAN_ANYD(Bx),
                BL_TO_FORTRAN_ANYD(By),
                BL_TO_FORTRAN_ANYD(Bz),
                BL_TO_FORTRAN_ANYD(q),
                BL_TO_FORTRAN_ANYD(cs[0]),
                BL_TO_FORTRAN_ANYD(cs[1]),
                BL_TO_FORTRAN_ANYD(cs[2]),
                BL_TO_FORTRAN_ANYD(source_in),
                BL_TO_FORTRAN_ANYD(srcQ));

    const int* lo1 = obx.loVect();
    const int* hi1 = obx.hiVect();

    srctoprim_mhd(lo1, hi1,
                  BL_TO_FORTRAN_ANYD(q),
                  BL_TO_FORTRAN_ANYD(source_in),
                  BL_TO_FORTRAN_ANYD(srcQ));

    check_for_mhd_cfl_violation(lo, hi,
                                BL_TO_FORTRAN_ANYD(q),
                                BL_TO_FORTRAN_ANYD(cs[0]),
                                BL_TO_FORTRAN_ANYD(cs[1]),
                                BL_TO_FORTRAN_ANYD(cs[2]));

    flatn.resize(obx, 1);
    auto flatn_arr = flatn.array();
    amrex::ParallelFor(obx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      flatn_arr(i,j,k) = 0.0;
    });

    mhd_flatten(lo1, hi1,
                BL_TO_FORTRAN_ANYD(q),
                BL_TO_FORTRAN_ANYD(flatn));
 

    // Interpolate Cell centered values to faces
    qp.resize(bx_gc, NQ * AMREX_SPACEDIM);
    qm.resize(bx_gc, NQ * AMREX_SPACEDIM);

    plm(lo, hi,
        BL_TO_FORTRAN_ANYD(q)
        BL_TO_FORTRAN_ANYD(flatn),
        BL_TO_FORTRAN_ANYD(Bx),
        BL_TO_FORTRAN_ANYD(By),
        BL_TO_FORTRAN_ANYD(Bz),
        BL_TO_FORTRAN_ANYD(qp),
        BL_TO_FORTRAN_ANYD(qm),
        BL_TO_FORTRAN_ANYD(srcQ),
        dx, dt);


    // Corner Couple and find the correct fluxes + electric fields

    // need to revisit these box sizes 
    nbx = amrex::surroundingNodes(bx, 0);
    nbxf = amrex::grow(nbx, IntVect(2, 3, 3));

    nby = amrex::surroundingNodes(bx, 1);
    nbyf = amrex::grow(nby, IntVect(3, 2, 3));

    nbz = amrex::surroundingNodes(bx, 2);
    nbzf = amrex::grow(nbz, IntVect(3, 3, 2));

    flxx.resize(nbxf, NVAR+3);
    Extmp.resize(nbxf);

    flxy.resize(nbyf, NVAR+3);
    Eytmp.resize(nbyf);

    flxz.resize(nbzf, NVAR+3);
    Eztmp.resize(nbzf);

    amrex::ParallelFor(nbxf, NVAR+3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      flxx(i,j,k,n) = 0.0_rt;
    });

    amrex::ParallelFor(nbyf, NVAR+3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      flxy(i,j,k,n) = 0.0_rt;
    });

    amrex::ParallelFor(nbzf, NVAR+3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      flxz(i,j,k,n) = 0.0_rt;
    });

    corner_transport(BL_TO_FORTRAN_ANYD(q),
                     BL_TO_FORTRAN_ANYD(qm),
                     BL_TO_FORTRAN_ANYD(qp),
                     BL_TO_FORTRAN_ANYD(flxx),
                     BL_TO_FORTRAN_ANYD(flxy),
                     BL_TO_FORTRAN_ANYD(flxz),
                     BL_TO_FORTRAN_ANYD(Extmp),
                     BL_TO_FORTRAN_ANYD(Eytmp),
                     BL_TO_FORTRAN_ANYD(Eztmp),
                     dx, dt);

    // Conservative update
    consup(uin,  uin_lo, uin_hi, &
              uout, uout_lo, uout_hi, &
              bcc, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
              flxx,flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
              flxy,flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
              flxz,flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
              lo ,hi ,dx ,dy ,dz ,dt)
 
  bxout = 0.0
  byout = 0.0
  bzout = 0.0

  !Step Five Magnetic Update
  call magup(bxin, bxin_lo, bxin_hi, &
             byin, byin_lo, byin_hi, &
             bzin, bzin_lo, bzin_hi, &
             bxout, bxout_lo, bxout_hi, &
             byout, byout_lo, byout_hi, &
             bzout, bzout_lo, bzout_hi, &
             Extemp, extemp_l1,extemp_l2,extemp_l3,extemp_h1,extemp_h2,extemp_h3, &
             Eytemp, eytemp_l1,eytemp_l2,eytemp_l3,eytemp_h1,eytemp_h2,eytemp_h3, &
             Eztemp, eztemp_l1,eztemp_l2,eztemp_l3,eztemp_h1,eztemp_h2,eztemp_h3, &
             lo, hi, dx, dy, dz, dt)

  
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),URHO) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),URHO)
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),UMX) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),UMX)
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),UMY) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),UMY)
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),UMZ) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),UMZ)
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),UEDEN) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),UEDEN)
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),UFS:UFS+nspec-1) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),UFS:UFS+nspec-1)

  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),URHO) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),URHO)
  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMX) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMX)
  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMY) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMY)
  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMZ) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMZ)
  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UEDEN) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UEDEN)
  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UFS:UFS+nspec-1) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UFS:UFS+nspec-1)


  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),URHO) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),URHO)
  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMX) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMX)
  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMY) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMY)
  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMZ) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMZ)
  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UEDEN) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UEDEN)
  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UFS:UFS+nspec-1) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UFS:UFS+nspec-1)


  Ex(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2), ex_lo(3):ex_hi(3)) = Extemp(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2),ex_lo(3):ex_hi(3))
  Ey(ey_lo(1):ey_hi(1),ey_lo(2):ey_hi(2), ey_lo(3):ey_hi(3)) = Eytemp(ey_lo(1):ey_hi(1),ey_lo(2):ey_hi(2),ey_lo(3):ey_hi(3))
  Ez(ez_lo(1):ez_hi(1),ez_lo(2):ez_hi(2), ez_lo(3):ez_hi(3)) = Eztemp(ez_lo(1):ez_hi(1),ez_lo(2):ez_hi(2),ez_lo(3):ez_hi(3))

  ! We are done with these here so can go ahead and free up the space
  call bl_deallocate(q)
  call bl_deallocate(bcc)
  call bl_deallocate(flatn)
  call bl_deallocate(cx)
  call bl_deallocate(cy)
  call bl_deallocate(cz)
  call bl_deallocate(csml)
  call bl_deallocate(srcQ)

  deallocate(qm)
  deallocate(qp)

  deallocate(flxx,flxy,flxz)
  deallocate(Extemp,Eytemp,Eztemp)



end subroutine ca_advance_mhd

