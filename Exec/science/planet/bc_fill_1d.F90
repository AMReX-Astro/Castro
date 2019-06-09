module bc_fill_module
  use bc_ext_fill_module
  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR
  use amrex_constants_module
  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine hypfill(adv, adv_l1, adv_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn


    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_h1
    integer,  intent(in   ) :: bc(1,2,NVAR)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,NVAR)

    integer  :: lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = adv_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, adv, lo, hi, NVAR, domlo, domhi, delta, xlo, bc)

  end subroutine hypfill


  subroutine ca_hypfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, &
                                   UFS, UTEMP, const_grav, &
                                   hse_zero_vels, hse_interp_temp, hse_reflect_vels, &
                                   xl_ext,xr_ext, EXT_HSE, EXT_INTERP
    use amrex_fort_module, only : rt => amrex_real
    use interpolate_module
    use eos_module
    use network, only: nspec
    use model_parser_module
    use amrex_error_module
    use eos_type_module

    include 'AMReX_bc_types.fi'




    integer,  intent(in   ) :: adv_l1, adv_h1
    integer,  intent(in   ) :: bc(1,2,NVAR)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,NVAR)


    call hypfill(adv, adv_l1, adv_h1, domlo, domhi, delta, xlo, time, bc)
    if ( bc(1,1,1).eq. FOEXTRAP .and. adv_l1.lt.domlo(1)) then
      if (xl_ext == EXT_HSE) then
         call ext_fill(adv,adv_l1,adv_h1, &
                    domlo,domhi,delta,xlo,time,bc)
      end if
    end if


  end subroutine ca_hypfill
  

  subroutine denfill(adv, adv_l1, adv_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1)

    integer :: lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = adv_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, adv, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine denfill


  subroutine ca_denfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1)

    call denfill(adv, adv_l1, adv_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_denfill



#ifdef GRAVITY
  subroutine phigravfill(phi, phi_l1, phi_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_l1, phi_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1)

    integer :: lo(3), hi(3)

    lo(1) = phi_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = phi_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, phi, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phigravfill


  subroutine ca_phigravfill(phi,phi_l1,phi_h1,domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    
    implicit none

    integer,  intent(in   ) :: phi_l1, phi_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1)

    call phigravfill(phi, phi_l1, phi_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phigravfill



  subroutine gravxfill(grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

    integer :: lo(3), hi(3)

    lo(1) = grav_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravxfill


  subroutine ca_gravxfill(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravxfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

    call gravxfill(grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravxfill



  subroutine gravyfill(grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

    integer :: lo(3), hi(3)

    lo(1) = grav_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravyfill


  subroutine ca_gravyfill(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravyfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

    call gravyfill(grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravyfill


  subroutine gravzfill(grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

    integer :: lo(3), hi(3)

    lo(1) = grav_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravzfill


  subroutine ca_gravzfill(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravzfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1)

    call gravzfill(grav, grav_l1, grav_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravzfill
#endif



#ifdef ROTATION
  subroutine phirotfill(phi, phi_l1, phi_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_l1, phi_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1)

    integer :: lo(3), hi(3)

    lo(1) = phi_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = phi_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, phi, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phirotfill


  subroutine ca_phirotfill(phi,phi_l1,phi_h1,domlo,domhi,delta,xlo,time,bc) &
                           bind(C, name="ca_phirotfill")

    implicit none

    integer,  intent(in   ) :: phi_l1, phi_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1)

    call phirotfill(phi, phi_l1, phi_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phirotfill



  subroutine rotxfill(rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

    integer :: lo(3), hi(3)

    lo(1) = rot_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = rot_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotxfill


  subroutine ca_rotxfill(rot,rot_l1,rot_h1,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotxfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

    call rotxfill(rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotxfill



  subroutine rotyfill(rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

    integer :: lo(3), hi(3)

    lo(1) = rot_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = rot_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotyfill


  subroutine ca_rotyfill(rot,rot_l1,rot_h1,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotyfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

    call rotyfill(rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotyfill



  subroutine rotzfill(rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

    integer :: lo(3), hi(3)

    lo(1) = rot_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = rot_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotzfill


  subroutine ca_rotzfill(rot,rot_l1,rot_h1,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotzfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1)

    call rotzfill(rot, rot_l1, rot_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotzfill
#endif



#ifdef REACTIONS
  subroutine reactfill(react, react_l1, react_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: react_l1, react_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: react(react_l1:react_h1)

    integer :: lo(3), hi(3)

    lo(1) = react_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = react_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, react, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine reactfill


  subroutine ca_reactfill(react,react_l1,react_h1,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    implicit none

    integer,  intent(in   ) :: react_l1, react_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: react(react_l1:react_h1)

    call reactfill(react, react_l1, react_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_reactfill
#endif



#ifdef RADIATION
  subroutine radfill(rad, rad_l1, rad_h1, domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rad_l1, rad_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rad(rad_l1:rad_h1)

    integer :: lo(3), hi(3)

    lo(1) = rad_l1
    lo(2) = 0
    lo(3) = 0
    hi(1) = rad_h1
    hi(2) = 0
    hi(3) = 0

    call amrex_filccn(lo, hi, rad, lo, hi, 1, domlo, domhi, delta, xlo, bc)
  end subroutine radfill


  subroutine ca_radfill(rad,rad_l1,rad_h1,domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_radfill")

    implicit none

    integer,  intent(in   ) :: rad_l1, rad_h1
    integer,  intent(in   ) :: bc(1,2)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    real(rt), intent(in   ) :: delta(1), xlo(1), time
    real(rt), intent(inout) :: rad(rad_l1:rad_h1)

    call radfill(rad, rad_l1, rad_h1, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_radfill
#endif

end module bc_fill_module
