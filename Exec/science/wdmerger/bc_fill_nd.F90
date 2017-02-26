module bc_fill_module

  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_hypfill")

     use bl_constants_module, only: HALF
     use meth_params_module, only: NVAR
     use prob_params_module, only: center, dim
     use probdata_module, only: fill_ambient_bc
     use wdmerger_util_module, only: fill_ambient

     implicit none

     include 'bc_types.fi'

     integer          :: adv_lo(3), adv_hi(3)
     integer          :: bc(dim,2,NVAR)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

     integer          :: i, j, k, n
     double precision :: loc(3)

     ! First, use the generic filling routines to make sure we have valid data everywhere
     ! on physical domain ghost cells.
     
     do n = 1,NVAR
        call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
     enddo
     
     ! Override the generic routine at the physical boundaries by
     ! setting the material to the ambient state. Note that we
     ! don't want to do this for interior/periodic boundaries, which have bc == 0.

     if (fill_ambient_bc) then

        do k = adv_lo(3), adv_hi(3)
           loc(3) = xlo(3) + dble(k - adv_lo(3) + HALF)*delta(3) - center(3)
           do j = adv_lo(2), adv_hi(2)
              loc(2) = xlo(2) + dble(j - adv_lo(2) + HALF)*delta(2) - center(2)
              do i = adv_lo(1), adv_hi(1)
                 loc(1) = xlo(1) + dble(i - adv_lo(1) + HALF)*delta(1) - center(1)

                 if (i .lt. domlo(1) .and. bc(1,1,1) .ne. 0) then
                    call fill_ambient(adv(i,j,k,:), loc, time)
                 else if (i .gt. domhi(1) .and. bc(1,2,1) .ne. 0) then
                    call fill_ambient(adv(i,j,k,:), loc, time)
                 else if (j .lt. domlo(2) .and. bc(2,1,1) .ne. 0) then
                    call fill_ambient(adv(i,j,k,:), loc, time)
                 else if (j .gt. domhi(2) .and. bc(2,2,1) .ne. 0) then
                    call fill_ambient(adv(i,j,k,:), loc, time)
                 else if (dim .eq. 3) then
                    if (k .lt. domlo(3) .and. bc(3,1,1) .ne. 0) then
                       call fill_ambient(adv(i,j,k,:), loc, time)
                    else if (k .gt. domhi(3) .and. bc(3,2,1) .ne. 0) then
                       call fill_ambient(adv(i,j,k,:), loc, time)
                    endif
                 endif                 

              enddo
           enddo
        enddo

     endif
        
   end subroutine ca_hypfill



   subroutine ca_denfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
        bind(C, name="ca_denfill")

     use prob_params_module, only: dim
     
     implicit none
     include 'bc_types.fi'

     integer          :: adv_lo(3), adv_hi(3)
     integer          :: bc(dim,2)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

     call filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

   end subroutine ca_denfill



#ifdef GRAVITY
   subroutine ca_gravxfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
        bind(C, name="ca_gravxfill")

     use prob_params_module, only: dim     
     
     implicit none
     include 'bc_types.fi'

     integer          :: grav_lo(3), grav_hi(3)
     integer          :: bc(dim,2)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

     call filcc_nd(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,bc)

   end subroutine ca_gravxfill



   subroutine ca_gravyfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
        bind(C, name="ca_gravyfill")

     use prob_params_module, only: dim     
     
     implicit none
     include 'bc_types.fi'

     integer          :: grav_lo(3), grav_hi(3)
     integer          :: bc(dim,2)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

     call filcc_nd(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,bc)

   end subroutine ca_gravyfill



   subroutine ca_gravzfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
        bind(C, name="ca_gravzfill")

     use prob_params_module, only: dim     
     
     implicit none
     include 'bc_types.fi'

     integer          :: grav_lo(3), grav_hi(3)
     integer          :: bc(dim,2)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

     call filcc_nd(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,bc)

   end subroutine ca_gravzfill



   subroutine ca_phigravfill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc) &
        bind(C, name="ca_phigravfill")

     use prob_params_module, only: dim     
     
     implicit none
     include 'bc_types.fi'

     integer          :: phi_lo(3), phi_hi(3)
     integer          :: bc(dim,2)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

     call filcc_nd(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,bc)

   end subroutine ca_phigravfill
#endif



#ifdef ROTATION
   subroutine ca_rotxfill(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,time,bc) &
        bind(C, name="ca_rotxfill")

     use prob_params_module, only: dim     
     
     implicit none
     include 'bc_types.fi'

     integer          :: rot_lo(3), rot_hi(3)
     integer          :: bc(dim,2)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

     call filcc_nd(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,bc)

   end subroutine ca_rotxfill



   subroutine ca_rotyfill(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,time,bc) &
        bind(C, name="ca_rotyfill")

     use prob_params_module, only: dim     
     
     implicit none
     include 'bc_types.fi'

     integer          :: rot_lo(3), rot_hi(3)
     integer          :: bc(dim,2)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

     call filcc_nd(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,bc)

   end subroutine ca_rotyfill



   subroutine ca_rotzfill(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,time,bc) &
        bind(C, name="ca_rotzfill")

     use prob_params_module, only: dim     
     
     implicit none
     include 'bc_types.fi'

     integer          :: rot_lo(3), rot_hi(3)
     integer          :: bc(dim,2)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

     call filcc_nd(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,bc)

   end subroutine ca_rotzfill



   subroutine ca_phirotfill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc) &
        bind(C, name="ca_phirotfill")

     use prob_params_module, only: dim     
     
     implicit none
     include 'bc_types.fi'

     integer          :: phi_lo(3), phi_hi(3)
     integer          :: bc(dim,2)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

     call filcc_nd(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,bc)

   end subroutine ca_phirotfill
#endif



#ifdef REACTIONS
   subroutine ca_reactfill(react,react_lo,react_hi,domlo,domhi,delta,xlo,time,bc) &
        bind(C, name="ca_reactfill")

     use prob_params_module, only: dim     
     
     implicit none
     include 'bc_types.fi'

     integer          :: react_lo(3), react_hi(3)
     integer          :: bc(dim,2)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: react(react_lo(1):react_hi(1),react_lo(2):react_hi(2),react_lo(3):react_hi(3))

     call filcc_nd(react,react_lo,react_hi,domlo,domhi,delta,xlo,bc)

   end subroutine ca_reactfill
#endif

end module bc_fill_module
