module bc_fill_module

#ifndef AMREX_USE_CUDA
  use bc_ext_fill_module
#endif
  use amrex_fort_module, only: rt => amrex_real
#ifdef CUDA
  use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  AMREX_LAUNCH subroutine hypfill(adv, adv_l1, adv_l2, adv_h1, adv_h2, &
                                  domlo, domhi, delta, xlo, time, bc)

    use meth_params_module, only: NVAR
    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_h1, adv_h2
    integer,  intent(in   ) :: bc(2,2,NVAR)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer  :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = adv_l2
    lo(3) = 0
    hi(1) = adv_h1
    hi(2) = adv_h2
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, adv, lo, hi, NVAR, domlo, domhi, delta, xlo, bc)

  end subroutine hypfill


  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_h1, adv_h2
    integer,  intent(in   ) :: bc(2,2,NVAR)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

#ifdef CUDA
    attributes(device) :: adv, adv_l1, adv_l2, adv_h1, adv_h2, bc, delta, xlo, time, domlo, domhi
#endif

    call hypfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_l1, adv_l2, adv_h1, adv_h2, domlo, domhi, delta, xlo, time, bc)

#ifndef AMREX_USE_CUdA
    ! process the external BCs here
    call ext_fill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc)
#endif
    
  end subroutine ca_hypfill


  AMREX_LAUNCH subroutine denfill(adv, adv_l1, adv_l2, adv_h1, adv_h2, &
                                  domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_h1, adv_h2
    integer,  intent(in   ) :: bc(2,2,1)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = adv_l2
    lo(3) = 0
    hi(1) = adv_h1
    hi(2) = adv_h2
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, adv, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine denfill


  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_h1, adv_h2
    integer,  intent(in   ) :: bc(2,2,1)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2)

#ifdef CUDA
    attributes(device) :: adv, adv_l1, adv_l2, adv_h1, adv_h2, bc, delta, xlo, time, domlo, domhi
#endif

    call denfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_l1, adv_l2, adv_h1, adv_h2, domlo, domhi, delta, xlo, time, bc)

#ifndef AMREX_USE_CUDA
    ! process the external BCs here
    call ext_denfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc)
#endif
    
  end subroutine ca_denfill


#ifdef GRAVITY
  AMREX_LAUNCH subroutine phigravfill(phi, phi_l1, phi_l2, phi_h1, phi_h2, &
                                      domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_h1, phi_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)

    lo(1) = phi_l1
    lo(2) = phi_l2
    lo(3) = 0
    hi(1) = phi_h1
    hi(2) = phi_h2
    hi(3) = 0

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, phi, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phigravfill


  subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                            domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    
    implicit none

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_h1, phi_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

#ifdef CUDA
    attributes(device) :: phi, phi_l1, phi_l2, phi_h1, phi_h2, bc, domlo, domhi, delta, xlo, time
#endif

    call phigravfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (phi, phi_l1, phi_l2, phi_h1, phi_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phigravfill

  
  AMREX_LAUNCH subroutine gravxfill(grav, grav_l1, grav_l2, grav_h1, grav_h2, &
                                      domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)
    integer :: bc_temp(2,2)

    lo(1) = grav_l1
    lo(2) = grav_l2
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = grav_h2
    hi(3) = 0
    
    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    if ( bc(2,1) == EXT_DIR .and. grav_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine gravxfill


  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravxfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)
    
#ifdef CUDA
    attributes(device) :: grav, grav_l1, grav_l2, grav_h1, grav_h2, bc
    attributes(device) :: domlo, domhi, delta, xlo, time
#endif

    call gravxfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (grav, grav_l1, grav_l2, grav_h1, grav_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravxfill


  AMREX_LAUNCH subroutine gravyfill(grav, grav_l1, grav_l2, grav_h1, grav_h2, &
                                    domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)
    integer :: bc_temp(2,2)

    lo(1) = grav_l1
    lo(2) = grav_l2
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = grav_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    if ( bc(2,1) == EXT_DIR .and. grav_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    
    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine gravyfill


  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravyfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

#ifdef CUDA
    attributes(device) :: grav, grav_l1, grav_l2, grav_h1, grav_h2
    attributes(device) :: domlo, domhi, delta, xlo, time, bc
#endif

    call gravyfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (grav, grav_l1, grav_l2, grav_h1, grav_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravyfill


  AMREX_LAUNCH subroutine gravzfill(grav, grav_l1, grav_l2, grav_h1, grav_h2, &
                                    domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)
    integer :: bc_temp(2,2)
    
    lo(1) = grav_l1
    lo(2) = grav_l2
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = grav_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    if ( bc(2,1) == EXT_DIR .and. grav_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif
    
    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine gravzfill


  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravzfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

#ifdef CUDA
    attributes(device) :: grav, grav_l1, grav_l2, grav_h1, grav_h2
    attributes(device) :: domlo, domhi, delta, xlo, time, bc
#endif

    call gravzfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (grav, grav_l1, grav_l2, grav_h1, grav_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravzfill
#endif



#ifdef ROTATION
  AMREX_LAUNCH subroutine phirotfill(phi, phi_l1, phi_l2, phi_h1, phi_h2, &
                                    domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_h1, phi_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)
    integer :: bc_temp(2,2)

    lo(1) = phi_l1
    lo(2) = phi_l2
    lo(3) = 0
    hi(1) = phi_h1
    hi(2) = phi_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    if ( bc(2,1) == EXT_DIR .and. phi_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif
    
    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, phi, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine phirotfill


  subroutine ca_phirotfill(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                           domlo,domhi,delta,xlo,time,bc) &
                           bind(C, name="ca_phirotfill")

    implicit none

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_h1, phi_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

#ifdef CUDA
    attributes(device) :: phi, phi_l1, phi_l2, phi_h1, phi_h2, bc, domlo, domhi, delta, xlo, time
#endif

    call phirotfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (phi, phi_l1, phi_l2, phi_h1, phi_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phirotfill


  AMREX_LAUNCH subroutine rotxfill(rot, rot_l1, rot_l2, rot_h1, rot_h2, &
                                    domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_h1, rot_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)
    integer :: bc_temp(2,2)

    lo(1) = rot_l1
    lo(2) = rot_l2
    lo(3) = 0
    hi(1) = rot_h1
    hi(2) = rot_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    if ( bc(2,1) == EXT_DIR .and. rot_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif
    
    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine rotxfill


  subroutine ca_rotxfill(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotxfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_h1, rot_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2)

#ifdef CUDA
    attributes(device) :: rot, rot_l1, rot_l2, rot_h1, rot_h2, bc, domlo, domhi, delta, xlo, time
#endif

    call rotxfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rot, rot_l1, rot_l2, rot_h1, rot_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotxfill


  AMREX_LAUNCH subroutine rotyfill(rot, rot_l1, rot_l2, rot_h1, rot_h2, &
                                    domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_h1, rot_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)
    integer :: bc_temp(2,2)

    lo(1) = rot_l1
    lo(2) = rot_l2
    lo(3) = 0
    hi(1) = rot_h1
    hi(2) = rot_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    if ( bc(2,1) == EXT_DIR .and. rot_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif
    
    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine rotyfill


  subroutine ca_rotyfill(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotyfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_h1, rot_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2)

#ifdef CUDA
    attributes(device) :: rot, rot_l1, rot_l2, rot_h1, rot_h2, bc, domlo, domhi, delta, xlo, time
#endif

    call rotyfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rot, rot_l1, rot_l2, rot_h1, rot_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotyfill

  AMREX_LAUNCH subroutine rotzfill(rot, rot_l1, rot_l2, rot_h1, rot_h2, &
                                   domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_h1, rot_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)
    integer :: bc_temp(2,2)

    lo(1) = rot_l1
    lo(2) = rot_l2
    lo(3) = 0
    hi(1) = rot_h1
    hi(2) = rot_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    if ( bc(2,1) == EXT_DIR .and. rot_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif
    
    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine rotzfill


  subroutine ca_rotzfill(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotzfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_h1, rot_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2)

#ifdef CUDA
    attributes(device) :: rot, rot_l1, rot_l2, rot_h1, rot_h2, bc, domlo, domhi, delta, xlo, time
#endif

    call rotzfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rot, rot_l1, rot_l2, rot_h1, rot_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotzfill
#endif



#ifdef REACTIONS
  AMREX_LAUNCH subroutine reactfill(react, react_l1, react_l2, react_h1, react_h2, &
                                    domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: react_l1, react_l2, react_h1, react_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: react(react_l1:react_h1,react_l2:react_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)
    integer :: bc_temp(2,2)

    lo(1) = react_l1
    lo(2) = react_l2
    lo(3) = 0
    hi(1) = react_h1
    hi(2) = react_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    if ( bc(2,1) == EXT_DIR .and. react_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif
    
    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, react, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine reactfill


  subroutine ca_reactfill(react,react_l1,react_l2,react_h1,react_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    implicit none

    integer,  intent(in   ) :: react_l1, react_l2, react_h1, react_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: react(react_l1:react_h1,react_l2:react_h2)

#ifdef CUDA
    attributes(device) :: react, react_l1, react_l2, react_h1, react_h2, bc, domlo, domhi, delta, xlo, time
#endif

    call reactfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (react, react_l1, react_l2, react_h1, react_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_reactfill
#endif



#ifdef RADIATION
  AMREX_LAUNCH subroutine radfill(rad, rad_l1, rad_l2, rad_h1, rad_h2, &
                                  domlo, domhi, delta, xlo, time, bc)

    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rad_l1, rad_l2, rad_h1, rad_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: rad(rad_l1:rad_h1,rad_l2:rad_h2)

    integer :: blo(3), bhi(3), lo(3), hi(3)
    integer :: bc_temp(2,2)

    lo(1) = rad_l1
    lo(2) = rad_l2
    lo(3) = 0
    hi(1) = rad_h1
    hi(2) = rad_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    if ( bc(2,1) == EXT_DIR .and. rad_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif
    
    call get_loop_bounds(blo, bhi, lo, hi)

    call amrex_filccn(blo, bhi, rad, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine radfill


  subroutine ca_radfill(rad,rad_l1,rad_l2,rad_h1,rad_h2, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_radfill")

    implicit none

    integer,  intent(in   ) :: rad_l1, rad_l2, rad_h1, rad_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: rad(rad_l1:rad_h1,rad_l2:rad_h2)

#ifdef CUDA
    attributes(device) :: rad, rad_l1, rad_l2, rad_h1, rad_h2, bc, domlo, domhi, delta, xlo, time
#endif

    call radfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (rad, rad_l1, rad_l2, rad_h1, rad_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_radfill
#endif

end module bc_fill_module
