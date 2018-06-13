module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_filcc_module, only: amrex_filccn
  use prob_params_module, only: dim
#ifdef CUDA
  use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  ! All subroutines in this file must be threadsafe because they are called
  ! inside OpenMP parallel regions.

  AMREX_LAUNCH subroutine hypfill(adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc)

    use meth_params_module, only: NVAR
    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    integer  :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, adv_lo, adv_hi)

    call amrex_filccn(blo, bhi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, delta, xlo, bc)

  end subroutine hypfill


  subroutine ca_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

#ifdef CUDA
    attributes(device) :: adv, adv_lo, adv_hi, bc, delta, xlo, domlo, domhi, time
#endif

    call hypfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_hypfill


  AMREX_LAUNCH subroutine denfill(adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(3,2,1)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, adv_lo, adv_hi)

    call amrex_filccn(blo, bhi, adv, adv_lo, adv_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine denfill


  subroutine ca_denfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_denfill")

    implicit none

    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

#ifdef CUDA
    attributes(device) :: adv, adv_lo, adv_hi, bc, delta, xlo, time, domlo, domhi
#endif

    call denfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_denfill


  
#ifdef GRAVITY
  AMREX_LAUNCH subroutine phigravfill(phi, phi_lo, phi_hi, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

    integer :: blo(3), bhi(3), lo(3), hi(3)

    call get_loop_bounds(blo, bhi, phi_lo, phi_hi)

    call amrex_filccn(blo, bhi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phigravfill


  subroutine ca_phigravfill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    implicit none

    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

#ifdef CUDA
    attributes(device) :: phi, phi_lo, phi_hi, domlo, domhi, delta, xlo, time, bc
#endif

    call phigravfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (phi, phi_lo, phi_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phigravfill
  

  AMREX_LAUNCH subroutine gravxfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, grav_lo, grav_hi)

    call amrex_filccn(blo, bhi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravxfill


  subroutine ca_gravxfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravxfill")

    implicit none

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

#ifdef CUDA
    attributes(device) :: grav, grav_lo, grav_hi
    attributes(device) :: domlo, domhi, delta, xlo, time, bc
#endif

    call gravxfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravxfill


  AMREX_LAUNCH subroutine gravyfill(grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, grav_lo, grav_hi)

    call amrex_filccn(blo, bhi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravyfill


  subroutine ca_gravyfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravyfill")

    implicit none

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

#ifdef CUDA
    attributes(device) :: grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc
#endif

    call gravyfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravyfill


  AMREX_LAUNCH subroutine gravzfill(grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, grav_lo, grav_hi)

    call amrex_filccn(blo, bhi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravzfill


  subroutine ca_gravzfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravzfill")

    implicit none

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

#ifdef CUDA
    attributes(device) :: grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc
#endif

    call gravzfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravzfill
#endif

  

#ifdef ROTATION
  AMREX_LAUNCH subroutine phirotfill(phi, phi_lo, phi_hi, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, phi_lo, phi_hi)

    call amrex_filccn(blo, bhi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phirotfill


  subroutine ca_phirotfill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc) &
                           bind(C, name="ca_phirotfill")

    implicit none

    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

#ifdef CUDA
    attributes(device) :: phi, phi_lo, phi_hi, bc, domlo, domhi, delta, xlo, time
#endif

    call phirotfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (phi, phi_lo, phi_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phirotfill
  

  AMREX_LAUNCH subroutine rotxfill(rot, rot_lo, rot_hi, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_lo(3), rot_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, rot_lo, rot_hi)

    call amrex_filccn(blo, bhi, rot, rot_lo, rot_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotxfill


  subroutine ca_rotxfill(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotxfill")

    implicit none

    integer,  intent(in   ) :: rot_lo(3), rot_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

#ifdef CUDA
    attributes(device) :: rot, rot_lo, rot_hi, bc, domlo, domhi, delta, xlo, time
#endif

    call rotxfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rot, rot_lo, rot_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotxfill


  AMREX_LAUNCH subroutine rotyfill(rot, rot_lo, rot_hi, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_lo(3), rot_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, rot_lo, rot_hi)

    call amrex_filccn(blo, bhi, rot, rot_lo, rot_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotyfill


  subroutine ca_rotyfill(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotyfill")

    implicit none

    integer,  intent(in   ) :: rot_lo(3), rot_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

#ifdef CUDA
    attributes(device) :: rot, rot_lo, rot_hi, bc, domlo, domhi, delta, xlo, time
#endif

    call rotyfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rot, rot_lo, rot_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotyfill


  AMREX_LAUNCH subroutine rotzfill(rot, rot_lo, rot_hi, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_lo(3), rot_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, rot_lo, rot_hi)

    call amrex_filccn(blo, bhi, rot, rot_lo, rot_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotzfill


  subroutine ca_rotzfill(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotzfill")

    implicit none

    integer,  intent(in   ) :: rot_lo(3), rot_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

#ifdef CUDA
    attributes(device) :: rot, rot_lo, rot_hi, bc, domlo, domhi, delta, xlo, time
#endif

    call rotzfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rot, rot_lo, rot_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotzfill
#endif  


#ifdef REACTIONS
  AMREX_LAUNCH subroutine reactfill(react, react_lo, react_hi, domlo, domhi, delta, xlo, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: react_lo(3), react_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: react(react_lo(1):react_hi(1),react_lo(2):react_hi(2),react_lo(3):react_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, react_lo, react_hi)

    call amrex_filccn(blo, bhi, react, react_lo, react_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine reactfill


  subroutine ca_reactfill(react,react_lo,react_hi,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    implicit none

    integer,  intent(in   ) :: react_lo(3), react_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: react(react_lo(1):react_hi(1),react_lo(2):react_hi(2),react_lo(3):react_hi(3))

#ifdef CUDA
    attributes(device) :: react, react_lo, react_hi, bc, domlo, domhi, delta, xlo
#endif

    call reactfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (react, react_lo, react_hi, domlo, domhi, delta, xlo, bc)

  end subroutine ca_reactfill
#endif



#ifdef RADIATION
  AMREX_LAUNCH subroutine radfill(rad, rad_lo, rad_hi, domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rad_lo(3), rad_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rad(rad_lo(1):rad_hi(1),rad_lo(2):rad_hi(2),rad_lo(3):rad_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, rad_lo, rad_hi)

    call amrex_filccn(blo, bhi, rad, rad_lo, rad_hi, 1, domlo, domhi, delta, xlo, time, bc)

  end subroutine radfill


  subroutine ca_radfill(rad,rad_lo,rad_hi,domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_radfill")

    implicit none

    integer,  intent(in   ) :: rad_lo(3), rad_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rad(rad_lo(1):rad_hi(1),rad_lo(2):rad_hi(2),rad_lo(3):rad_hi(3))

#ifdef CUDA
    attributes(device) :: rad, rad_lo, rad_hi, bc, domlo, domhi, delta, xlo, time
#endif

    call radfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rad, rad_lo, rad_hi, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_radfill
#endif
  
end module bc_fill_module
