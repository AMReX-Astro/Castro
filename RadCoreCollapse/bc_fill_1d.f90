
      subroutine ca_hypfill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)

      use meth_params_module, only : NVAR

      implicit none

      include 'bc_types.fi'
      integer adv_l1,adv_h1
      integer bc(1,2,*)
      integer domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision adv(adv_l1:adv_h1,NVAR)

      double precision state(NVAR)
      double precision staten(NVAR)

      integer n

      do n = 1,NVAR
         call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
             domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      end

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_denfill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_h1
      integer bc(1,2,*)
      integer domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision adv(adv_l1:adv_h1)

      call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

      end

!-----------------------------------------------------------------------

      subroutine ca_gravxfill(grav,grav_l1,grav_h1,&
                              domlo,domhi,delta,xlo,time,bc)
 
      use probdata_module
      implicit none
      include 'bc_types.fi'
 
      integer :: grav_l1,grav_h1
      integer :: bc(1,2,*)
      integer :: domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision grav(grav_l1:grav_h1)
 
      double precision :: ri,rim1
      integer          :: i
 
      call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)
 
!     Outflow boundary condition for gravity at hi r
      if ( bc(1,2,1).eq. FOEXTRAP .and. grav_h1.gt.domhi(1)) then
         do i = domhi(1)+1, grav_h1
            ri   = (dble(i  )+0.5d0) * delta(1)
            rim1 = (dble(i-1)+0.5d0) * delta(1)
            grav(i) = grav(i-1) * (ri/rim1)**2
         end do
      end if
 
      end subroutine ca_gravxfill

! ::: -----------------------------------------------------------

      subroutine ca_radfill(rad,rad_l1,rad_h1,&
                              domlo,domhi,delta,xlo,time,bc)
 
      use probdata_module
      implicit none
      include 'bc_types.fi'
 
      integer :: rad_l1,rad_h1
      integer :: bc(1,2,*)
      integer :: domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision rad(rad_l1:rad_h1)
 
      call filcc(rad,rad_l1,rad_h1,domlo,domhi,delta,xlo,bc)
 
      end subroutine ca_radfill

! ::: -----------------------------------------------------------
