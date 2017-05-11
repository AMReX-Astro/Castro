module bc_fill_module

  use bc_ext_fill_module

  use bl_fort_module, only : rt => c_real
  implicit none

  include 'bc_types.fi'

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR

    use bl_fort_module, only : rt => c_real
    integer          :: adv_l1,adv_l2,adv_h1,adv_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer          :: n

    do n = 1,NVAR
       call filcc(adv(:,:,n),adv_l1,adv_l2,adv_h1,adv_h2, &
            domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

    ! process the external BCs here
    call ext_fill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                  domlo,domhi,delta,xlo,time,bc)

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_denfill")

    use bl_fort_module, only : rt => c_real
    integer          :: adv_l1,adv_l2,adv_h1,adv_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2)

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    ! process the external BCs here
    call ext_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                     domlo,domhi,delta,xlo,time,bc)

  end subroutine ca_denfill



  subroutine ca_phigravfill(phi,phi_l1,phi_l2, &
                            phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    use bl_fort_module, only : rt => c_real
    integer          :: phi_l1,phi_l2,phi_h1,phi_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill



  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravxfill")

    use bl_fort_module, only : rt => c_real
    integer          :: grav_l1,grav_l2,grav_h1,grav_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: bc_temp(2,2)
    
    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:,1)
    
    if (bc(2,1,1) == EXT_DIR .and. grav_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif
    
   call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)
   
  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravyfill")

    use bl_fort_module, only : rt => c_real
    integer          :: grav_l1,grav_l2,grav_h1,grav_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: bc_temp(2,2)

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1) == EXT_DIR .and. grav_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravzfill")

    use bl_fort_module, only : rt => c_real
    integer          :: grav_l1,grav_l2,grav_h1,grav_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: bc_temp(2,2)

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1) == EXT_DIR .and. grav_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravzfill



  subroutine ca_reactfill(react,react_l1,react_l2, &
                          react_h1,react_h2,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    use bl_fort_module, only : rt => c_real
    integer          :: react_l1,react_l2,react_h1,react_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: react(react_l1:react_h1,react_l2:react_h2)

    integer :: bc_temp(2,2)

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1) == EXT_DIR .and. react_l2 < domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(react,react_l1,react_l2,react_h1,react_h2, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_reactfill



  subroutine ca_radfill(rad,rad_l1,rad_l2, &
       rad_h1,rad_h2,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_radfill")

    use bl_fort_module, only : rt => c_real
    use prob_rad_params_module, only : bcval_lo, bcval_hi
    use fluxlimiter_module, only : FLDlambda, limiter
    use probdata_module
    use parmparse_module, only : parmparse_build, parmparse_destroy, ParmParse
    
    
    type (ParmParse) :: pp 

    integer :: rad_l1,rad_l2,rad_h1,rad_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         rad(rad_l1:rad_h1,rad_l2:rad_h2)
    
    real(rt):: r(rad_l1:rad_h1,rad_l2:rad_h2)  

    integer :: bc_temp(2,2)

    integer :: i, j


     real:: c_light = 2.99e10
     real:: F_in
     real(rt) :: const_k, exp_m, exp_n 
     real:: D 
     real:: T, rho, k

     call parmparse_build(pp,"radiation")
     
     call pp%query("const_kappa_r", const_k)
     call pp%query("kappa_r_exp_m", exp_m)
     call pp%query("kappa_r_exp_n", exp_n)

     call parmparse_destroy(pp)
     F_in = bcval_lo(2)


     ! handle an external BC via extrapolation here 
     bc_temp(:,:) = bc(:,:,1)

     if ( bc(2,1,1) == EXT_DIR .and. rad_l2 < domlo(2)) then
        bc_temp(2,1) = FOEXTRAP
     endif

     call filcc(rad,rad_l1,rad_l2,rad_h1,rad_h2,domlo,domhi,delta,xlo,bc)

     !compute the gradient of Er to then get R ??  .... 


     !for now using the last one? cause it is the same at initial time at end of domain
     !maybe do an interpolation and find opacity in each cell?
     T = T_0(npts)
     rho = rho_0(npts)
     k = const_k * rho**exp_m * T**(-exp_n) !rosseland mean opacity
      

     D = FLDlambda(0.0e0_rt, limiter) / k  !R is 0? at initialization  

      !fill ghost cell 
     rad(rad_l1:rad_h1,rad_l2) = 2.0*D/(2.0*D + delta(2)) * (2.0*delta(2)/(c_light*D)*F_in + rad(rad_l1:rad_h1,rad_l2+1))


  end subroutine ca_radfill

end module bc_fill_module
