module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_filcc_module, only: filccn

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  AMREX_LAUNCH subroutine hypfill(adv, adv_l1, adv_h1, domlo, domhi, delta, xlo, time, bc)

    use meth_params_module, only: NVAR
    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_h1
    integer,  intent(in   ) :: bc(1,2,NVAR)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,NVAR)

    integer  :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = adv_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, adv, lo, hi, NVAR, domlo, domhi, delta, xlo, bc)

  end subroutine hypfill


  subroutine ca_hypfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_h1
    integer,  intent(in   ) :: bc(1,2,NVAR)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,NVAR)

#ifdef CUDA
    attributes(device) :: adv, adv_l1, adv_h1, bc, delta, xlo, time, domlo, domhi
#endif

    call hypfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_l1, adv_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_hypfill
  

  AMREX_LAUNCH subroutine denfill(adv, adv_l1, adv_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_h1
    integer,  intent(in   ) :: bc(1,2,1)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = adv_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, adv, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine denfill


  subroutine ca_denfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_h1
    integer,  intent(in   ) :: bc(1,2,1)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1)

#ifdef CUDA
    attributes(device) :: adv, adv_l1, adv_h1, bc, delta, xlo, time, domlo, domhi
#endif

    call denfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_l1, adv_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_denfill



#ifdef GRAVITY
  AMREX_LAUNCH subroutine phigravfill(phi, phi_l1, phi_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_l1, phi_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = phi_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = phi_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, phi, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phigravfill


  subroutine ca_phigravfill(phi,phi_l1,phi_h1,domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    
    implicit none

    integer,  intent(in   ) :: phi_l1, phi_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1)

#ifdef CUDA
    attributes(device) :: phi, phi_l1, phi_h1, bc, domlo, domhi, delta, xlo, time
#endif

    call phigravfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (phi, phi_l1, phi_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phigravfill



  AMREX_LAUNCH subroutine gravxfill(grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = grav_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravxfill


  subroutine ca_gravxfill(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravxfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

#ifdef CUDA
    attributes(device) :: grav, grav_l1, grav_h1
    attributes(device) :: domlo, domhi, delta, xlo, time, bc
#endif

    call gravxfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravxfill



  AMREX_LAUNCH subroutine gravyfill(grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = grav_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravyfill


  subroutine ca_gravyfill(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravyfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

#ifdef CUDA
    attributes(device) :: grav, grav_l1, grav_h1
    attributes(device) :: domlo, domhi, delta, xlo, time, bc
#endif

    call gravyfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravyfill


  AMREX_LAUNCH subroutine gravzfill(grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = grav_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravzfill


  subroutine ca_gravzfill(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravzfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

#ifdef CUDA
    attributes(device) :: grav, grav_l1, grav_h1
    attributes(device) :: domlo, domhi, delta, xlo, time, bc
#endif

    call gravzfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravzfill
#endif



#ifdef ROTATION
  AMREX_LAUNCH subroutine phirotfill(phi, phi_l1, phi_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_l1, phi_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = phi_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = phi_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, phi, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phirotfill


  subroutine ca_phirotfill(phi,phi_l1,phi_h1,domlo,domhi,delta,xlo,time,bc) &
                           bind(C, name="ca_phirotfill")

    implicit none

    integer,  intent(in   ) :: phi_l1, phi_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1)

#ifdef CUDA
    attributes(device) :: phi, phi_l1, phi_h1, bc, domlo, domhi, delta, xlo, time
#endif

    call phirotfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (phi, phi_l1, phi_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phirotfill



  AMREX_LAUNCH subroutine rotxfill(rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = rot_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = rot_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotxfill


  subroutine ca_rotxfill(rot,rot_l1,rot_h1,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotxfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

#ifdef CUDA
    attributes(device) :: rot, rot_l1, rot_h1, bc, domlo, domhi, delta, xlo, time
#endif

    call rotxfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotxfill



  AMREX_LAUNCH subroutine rotyfill(rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = rot_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = rot_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotyfill


  subroutine ca_rotyfill(rot,rot_l1,rot_h1,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotyfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

#ifdef CUDA
    attributes(device) :: rot, rot_l1, rot_h1, bc, domlo, domhi, delta, xlo, time
#endif

    call rotyfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotyfill



  AMREX_LAUNCH subroutine rotzfill(rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = rot_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = rot_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotzfill


  subroutine ca_rotzfill(rot,rot_l1,rot_h1,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotzfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

#ifdef CUDA
    attributes(device) :: rot, rot_l1, rot_h1, bc, domlo, domhi, delta, xlo, time
#endif

    call rotzfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotzfill
#endif



#ifdef REACTIONS
  AMREX_LAUNCH subroutine reactfill(react, react_l1, react_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: react_l1, react_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: react(react_l1:react_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = react_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = react_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, react, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine reactfill


  subroutine ca_reactfill(react,react_l1,react_h1,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    implicit none

    integer,  intent(in   ) :: react_l1, react_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: react(react_l1:react_h1)

#ifdef CUDA
    attributes(device) :: react, react_l1, react_h1, bc, domlo, domhi, delta, xlo, time
#endif

    call reactfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (react, react_l1, react_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_reactfill
#endif



#ifdef RADIATION
  AMREX_LAUNCH subroutine radfill(rad, rad_l1, rad_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rad_l1, rad_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rad(rad_l1:rad_h1)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = rad_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = rad_h1
    hi(2) = 0
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, rad, lo, hi, 1, domlo, domhi, delta, xlo, time, bc)

  end subroutine radfill


  subroutine ca_radfill(rad,rad_l1,rad_h1,domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_radfill")

    implicit none

    integer,  intent(in   ) :: rad_l1, rad_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rad(rad_l1:rad_h1)

#ifdef CUDA
    attributes(device) :: rad, rad_l1, rad_h1, bc, domlo, domhi, delta, xlo, time
#endif

    call radfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rad, rad_l1, rad_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_radfill
#endif

end module bc_fill_module
