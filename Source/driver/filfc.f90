module fc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for face centered data. 
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: q          <=  array to fill
! ::: q_l1,q_l2,q_l3   => index lower bounds of q
! ::: q_h1,q_h2,q_h3   => index upper bounds of q
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of q array
! ::: bc        => array of boundary flags bc(SPACEDIM,lo:hi)
! ::: dir       => direction of face centered data (i.e. x-direction)
! ::: 
! ::: -----------------------------------------------------------

      subroutine filfc(q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,domlo,domhi,dx,xlo,bc,dir)

      integer    domlo(3), domhi(3)
      integer    q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
      real(rt)   xlo(3), dx(3)
      real(rt)   q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      integer    bc(3,2)
      integer    dir

      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    is,  ie,  js,  je,  ks,  ke
      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    i, j, k, n

      is = max(q_l1,domlo(1))
        if(dir.eq.1) then
      ie = min(q_h1,domhi(1)+1)
        else
      ie = min(q_h1,domhi(1))
        endif
      js = max(q_l2,domlo(2))
        if(dir.eq.2) then
      je = min(q_h2,domhi(2)+1)
        else
      je = min(q_h2,domhi(2))
        endif
      ks = max(q_l3,domlo(3))
        if(dir.eq.3) then
      ke = min(q_h3,domhi(3)+1)
        else
      ke = min(q_h3,domhi(3))
        endif

      nlft = max(0,domlo(1)-q_l1)
      nrgt = max(0,q_h1-(domhi(1)))
      nbot = max(0,domlo(2)-q_l2)
      ntop = max(0,q_h2-(domhi(2)))
      ndwn = max(0,domlo(3)-q_l3)
      nup  = max(0,q_h3-(domhi(3)))

      if (dir.eq.1) then
!
!     ::::: X-face extending in lo-x direction
!
      if (q_l1 .lt. domlo(1)) then
         ilo = domlo(1)

         if (bc(1,1) .eq. FOEXTRAP .or. &
             bc(1,1) .eq. HOEXTRAP) then

            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-1,j,k) = 2.d0*q(ilo,j,k) - q(ilo+1,j,k)
            end do
            end do

            do n = 2, nlft
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-n,j,k) = q(ilo-1,j,k)
               end do
               end do
            end do

         else if (bc(1,1) .eq. REFLECT_EVEN) then

            do n = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-n,j,k) = q(ilo+n,j,k)
            end do
            end do
            end do
         else !periodic
                if(q_h1.gt.domhi(1)) then
                 ihi = domhi(1)
                   do n = 1,nlft
                        do k = q_l3, q_h3
                        do j = q_l2, q_h2
                        q(ilo-n,j,k) = q(ihi-n+1,j,k)
                        enddo
                        enddo
                   enddo
                endif
         end if
      end if
!
!     ::::: X-face extending in hi-x direction
!
      if (q_h1 .gt. domhi(1)+1) then
         ihi = domhi(1)!+1

         if (bc(1,2) .eq. FOEXTRAP .or. &
             bc(1,2) .eq. HOEXTRAP) then

            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+1,j,k) = 2.d0*q(ihi,j,k) - q(ihi-1,j,k)
            end do
            end do

            do n = 2, nrgt
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+n,j,k) = q(ihi+1,j,k)
               end do
               end do
            end do

         else if (bc(1,2) .eq. REFLECT_EVEN) then

            do n = 1, nrgt
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+n,j,k) = q(ihi-n,j,k)
            end do
            end do
            end do
         else 
            if(q_l1.lt.domlo(1)) then
            ilo = domlo(1)
            do n = 1, nrgt
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+n,j,k) = q(ilo+(n-1),j,k)
            end do
            end do
            end do
            endif
         end if
      end if
!
!     ::::: X-face extending in lo-y direction
!
      if (q_l2 .lt. domlo(2)) then

         jlo = domlo(2)

         if (bc(2,1) .eq. FOEXTRAP .or. &
             bc(2,1) .eq. HOEXTRAP) then

            do n = 1, nbot
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-n,k) = q(i,jlo,k)
               end do
               end do
            end do

         else if (bc(2,1) .eq. REFLECT_EVEN) then

            do n = 1, nbot
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-n,k) = q(i,jlo+n-1,k)
            end do
            end do
            end do
         else 
            if(q_h2.gt.domhi(2)) then
            jhi = domhi(2)
            do n = 1, nbot
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-n,k) = q(i,jhi-(n-1),k)
            end do
            end do
            end do
            endif
         end if
      end if

!
!     ::::: X-face extending in hi-y direction
!
      if (q_h2 .gt. domhi(2)) then

         jhi = domhi(2)

         if (bc(2,2) .eq. FOEXTRAP .or. &
             bc(2,2) .eq. HOEXTRAP) then

            do n = 1, ntop
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+n,k) = q(i,jhi,k)
               end do
               end do
            end do

         else if (bc(2,2) .eq. REFLECT_EVEN) then

            do n = 1, ntop
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+n,k) = q(i,jhi-n+1,k)
            end do
            end do
            end do
         else 
            if(q_l2.lt.domlo(2)) then
            jlo = domlo(2)
            do n = 1, ntop
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+n,k) = q(i,jlo+(n-1),k)
            end do
            end do
            end do
            endif
         end if
      end if

!
!     ::::: X-face extending in lo-z direction
!
      if (q_l3 .lt. domlo(3)) then

         klo = domlo(3)

         if (bc(3,1) .eq. FOEXTRAP .or. &
             bc(3,1) .eq. HOEXTRAP) then

            do n = 1, ndwn
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-n) = q(i,j,klo)
               end do
               end do
            end do

         else if (bc(3,1) .eq. REFLECT_EVEN) then

            do n = 1, ndwn
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-n) = q(i,j,klo+n-1)
            end do
            end do
            end do
         else 
            if(q_h3.gt.domhi(3)) then
            khi = domhi(3)
            do n = 1, ndwn
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-n) = q(i,j,khi-(n-1))
            end do
            end do
            end do
            endif
         end if
      end if

!
!     ::::: X-face extending in hi-z direction
!
      if (q_h3 .gt. domhi(3)) then

         khi = domhi(3)

         if (bc(3,2) .eq. FOEXTRAP .or. &
             bc(3,2) .eq. HOEXTRAP) then

            do n = 1, nup
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+n) = q(i,j,khi)
               end do
               end do
            end do

         else if (bc(3,2) .eq. REFLECT_EVEN) then

            do n = 1, nup
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+n) = q(i,j,khi-n+1)
            end do
            end do
            end do
         else 
            if(q_l3.lt.domlo(3)) then
            klo = domlo(3)
            do n = 1, nup
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+n) = q(i,j,klo+(n-1))
            end do
            end do
            end do
            endif
         end if
      end if

      else if (dir.eq.2) then
!
!     ::::: Y-face extending in lo-x direction
!
      if (q_l1 .lt. domlo(1)) then
         ilo = domlo(1)

         if (bc(1,1) .eq. FOEXTRAP .or. &
             bc(1,1) .eq. HOEXTRAP) then

            do n = 1, nlft
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-n,j,k) = q(ilo,j,k)
               end do
               end do
            end do

         else if (bc(1,1) .eq. REFLECT_EVEN) then

            do n = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-n,j,k) = q(ilo+n-1,j,k)
            end do
            end do
            end do
         else 
        if(q_h1.gt.domhi(1)) then
                 ihi = domhi(1)
                   do n = 1,nlft
                        do k = q_l3, q_h3
                        do j = q_l2, q_h2
                        q(ilo-n,j,k) = q(ihi-n+1,j,k)
                        enddo
                        enddo
                   enddo
                endif
         end if
      end if
!
!     ::::: Y-face extending in hi-x direction
!
      if (q_h1 .gt. domhi(1)) then
         ihi = domhi(1)

         if (bc(1,2) .eq. FOEXTRAP .or. &
             bc(1,2) .eq. HOEXTRAP) then

            do n = 1, nrgt
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+n,j,k) = q(ihi,j,k)
               end do
               end do
            end do

         else if (bc(1,2) .eq. REFLECT_EVEN) then

            do n = 1, nrgt
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+n,j,k) = q(ihi-n+1,j,k)
            end do
            end do
            end do
         else 
            if(q_l1.lt.domlo(1)) then
            ilo = domlo(1)
            do n = 1, nrgt
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+n,j,k) = q(ilo+(n-1),j,k)
            end do
            end do
            end do
            endif
         end if
      end if
!
!     ::::: Y-face extending in lo-y direction
!
      if (q_l2 .lt. domlo(2)) then

         jlo = domlo(2)

         if (bc(2,1) .eq. FOEXTRAP .or. &
             bc(2,1) .eq. HOEXTRAP) then

            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-1,k) = 2.d0*q(i,jlo,k) - q(i,jlo+1,k)
            end do
            end do

            do n = 2, nbot
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-n,k) = q(i,jlo-1,k)
               end do
               end do
            end do

         else if (bc(2,1) .eq. REFLECT_EVEN) then

            do n = 1, nbot
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-n,k) = q(i,jlo+n,k)
            end do
            end do
            end do
         else 
            if(q_h2.gt.domhi(2)) then
            jhi = domhi(2)
            do n = 1, nbot
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-n,k) = q(i,jhi-(n-1),k)
            end do
            end do
            end do
            endif
         end if
      end if

!
!     ::::: Y-face extending in hi-y direction
!
      if (q_h2 .gt. domhi(2)+1) then

         jhi = domhi(2)

         if (bc(2,2) .eq. FOEXTRAP .or. &
             bc(2,2) .eq. HOEXTRAP) then

            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+1,k) = 2.d0*q(i,jhi,k) - q(i,jhi-1,k)
            end do
            end do

            do n = 2, ntop
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+n,k) = q(i,jhi+1,k)
               end do
               end do
            end do

         else if (bc(2,2) .eq. REFLECT_EVEN) then

            do n = 1, ntop
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+n,k) = q(i,jhi-n,k)
            end do
            end do
            end do
         else 
            if(q_l2.lt.domlo(2)) then
            jlo = domlo(2)
            do n = 1, ntop
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+n,k) = q(i,jlo+(n-1),k)
            end do
            end do
            end do
            endif
         end if
      end if

!
!     ::::: Y-face extending in lo-z direction
!
      if (q_l3 .lt. domlo(3)) then

         klo = domlo(3)

         if (bc(3,1) .eq. FOEXTRAP .or. &
             bc(3,1) .eq. HOEXTRAP) then

            do n = 1, ndwn
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-n) = q(i,j,klo)
               end do
               end do
            end do

         else if (bc(3,1) .eq. REFLECT_EVEN) then

            do n = 1, ndwn
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-n) = q(i,j,klo+n-1)
            end do
            end do
            end do
         else 
            if(q_h3.gt.domhi(3)) then
            khi = domhi(3)
            do n = 1, ndwn
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-n) = q(i,j,khi-(n-1))
            end do
            end do
            end do
            endif
         end if
      end if

!
!     ::::: Y-face extending in hi-z direction
!
      if (q_h3 .gt. domhi(3)) then

         khi = domhi(3)

         if (bc(3,2) .eq. FOEXTRAP .or. &
             bc(3,2) .eq. HOEXTRAP) then

            do n = 1, nup
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+n) = q(i,j,khi)
               end do
               end do
            end do

         else if (bc(3,2) .eq. REFLECT_EVEN) then

            do n = 1, nup
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+n) = q(i,j,khi-n+1)
            end do
            end do
            end do
         else 
            if(q_l3.lt.domlo(3)) then
            klo = domlo(3)
            do n = 1, nup
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+n) = q(i,j,klo+(n-1))
            end do
            end do
            end do
            endif
         end if
      end if

      else if (dir.eq.3) then
!
!     ::::: Z-face extending in lo-x direction
!
      if (q_l1 .lt. domlo(1)) then
         ilo = domlo(1)

         if (bc(1,1) .eq. FOEXTRAP .or. &
             bc(1,1) .eq. HOEXTRAP) then

            do n = 1, nlft
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-n,j,k) = q(ilo,j,k)
               end do
               end do
            end do

         else if (bc(1,1) .eq. REFLECT_EVEN) then

            do n = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-n,j,k) = q(ilo+n-1,j,k)
            end do
            end do
            end do
         else 
        if(q_h1.gt.domhi(1)) then
                 ihi = domhi(1)
                   do n = 1,nlft
                        do k = q_l3, q_h3
                        do j = q_l2, q_h2
                        q(ilo-n,j,k) = q(ihi-n+1,j,k)
                        enddo
                        enddo
                   enddo
                endif
         end if
      end if
!
!     ::::: Z-face extending in hi-x direction
!
      if (q_h1 .gt. domhi(1)) then
         ihi = domhi(1)

         if (bc(1,2) .eq. FOEXTRAP .or. &
             bc(1,2) .eq. HOEXTRAP) then

            do n = 1, nrgt
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+n,j,k) = q(ihi,j,k)
               end do
               end do
            end do

         else if (bc(1,2) .eq. REFLECT_EVEN) then

            do n = 1, nrgt
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+n,j,k) = q(ihi-n+1,j,k)
            end do
            end do
            end do
         else 
            if(q_l1.lt.domlo(1)) then
            ilo = domlo(1)
            do n = 1, nrgt
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+n,j,k) = q(ilo+(n-1),j,k)
            end do
            end do
            end do
            endif
         end if
      end if
!
!     ::::: Z-face extending in lo-y direction
!
      if (q_l2 .lt. domlo(2)) then

         jlo = domlo(2)

         if (bc(2,1) .eq. FOEXTRAP .or. &
             bc(2,1) .eq. HOEXTRAP) then

            do n = 1, nbot
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-n,k) = q(i,jlo,k)
               end do
               end do
            end do

         else if (bc(2,1) .eq. REFLECT_EVEN) then

            do n = 1, nbot
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-n,k) = q(i,jlo+n-1,k)
            end do
            end do
            end do
         else 
            if(q_h2.gt.domhi(2)) then
            jhi = domhi(2)
            do n = 1, nbot
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-n,k) = q(i,jhi-(n-1),k)
            end do
            end do
            end do
            endif
         end if
      end if

!
!     ::::: Z-face extending in hi-y direction
!
      if (q_h2 .gt. domhi(2)+1) then

         jhi = domhi(2)

         if (bc(2,2) .eq. FOEXTRAP .or. &
             bc(2,2) .eq. HOEXTRAP) then

            do n = 1, ntop
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+n,k) = q(i,jhi,k)
               end do
               end do
            end do

         else if (bc(2,2) .eq. REFLECT_EVEN) then

            do n = 1, ntop
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+n,k) = q(i,jhi-n+1,k)
            end do
            end do
            end do
         else 
            if(q_l2.lt.domlo(2)) then
            jlo = domlo(2)
            do n = 1, ntop
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+n,k) = q(i,jlo+(n-1),k)
            end do
            end do
            end do
            endif
         end if
      end if

!
!     ::::: Z-face extending in lo-z direction
!
      if (q_l3 .lt. domlo(3)) then

         klo = domlo(3)

         if (bc(3,1) .eq. FOEXTRAP .or. &
             bc(3,1) .eq. HOEXTRAP) then

            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-1) = 2.d0*q(i,j,klo) - q(i,j,klo+1)
            end do
            end do

            do n = 2, ndwn
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-n) = q(i,j,klo-1)
               end do
               end do
            end do

         else if (bc(3,1) .eq. REFLECT_EVEN) then

            do n = 1, ndwn
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-n) = q(i,j,klo+n)
            end do
            end do
            end do
         else 
            if(q_h3.gt.domhi(3)) then
            khi = domhi(3)
            do n = 1, ndwn
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-n) = q(i,j,khi-(n-1))
            end do
            end do
            end do
            endif
         end if
      end if

!
!     ::::: Z-face extending in hi-z direction
!
      if (q_h3 .gt. domhi(3)+1) then

         khi = domhi(3)!+1

         if (bc(3,2) .eq. FOEXTRAP .or. &
             bc(3,2) .eq. HOEXTRAP) then

            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+1) = 2.d0*q(i,j,khi) - q(i,j,khi-1)
            end do
            end do

            do n = 2, nup
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+n) = q(i,j,khi+1)
               end do
               end do
            end do

         else if (bc(3,2) .eq. REFLECT_EVEN) then

            do n = 1, nup
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+n) = q(i,j,khi-n)
            end do
            end do
            end do
         else 
            if(q_l3.lt.domlo(3)) then
            klo = domlo(3)
            do n = 1, nup
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+n) = q(i,j,klo+(n-1))
            end do
            end do
            end do
            endif
         end if
      end if
 end if

!
!    First correct the i-j edges and all corners for x-dir
!
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP .or. bc(1,1) .eq. FOEXTRAP) .and.&
          (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP .or. bc(2,1) .eq. FOEXTRAP) ) then
         if (jlo+2 .le. je) then
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = 0.5d0 * 0.125d0 * &
                   (15*q(ilo-1,jlo,k) - 10*q(ilo-1,jlo+1,k) + &
                   3*q(ilo-1,jlo+2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = 0.25d0* &
                   (3*q(ilo-1,jlo,k) - q(ilo-1,jlo+1,k))
            end do
         end if

         if (ilo+2 .le. ie) then
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + 0.5d0 * 0.125d0 * &
                   (15*q(ilo,jlo-1,k) - 10*q(ilo+1,jlo-1,k) + &
                   3*q(ilo+2,jlo-1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + 0.25d0 * &
                   (3*q(ilo,jlo-1,k) - q(ilo+1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP  .or. bc(3,1) .eq. FOEXTRAP) then
            if (klo+2 .le. ke) then
               q(ilo-1,jlo-1,klo-1) = 0.125d0 * ( &
                   (15*q(ilo-1,jlo-1,klo) - 10*q(ilo-1,jlo-1,klo+1) +&
                   3*q(ilo-1,jlo-1,klo+2)) )
            else
               q(ilo-1,jlo-1,klo-1) = 0.5d0 * &
                   (3*q(ilo-1,jlo-1,klo) - q(ilo-1,jlo-1,klo+1))
             endif
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP  .or. bc(3,2) .eq. FOEXTRAP) then
            if (khi-2 .ge. ks) then
               q(ilo-1,jlo-1,khi+1) = 0.125d0 * ( &
                   (15*q(ilo-1,jlo-1,khi) - 10*q(ilo-1,jlo-1,khi-1) + &
                   3*q(ilo-1,jlo-1,khi-2)) )
            else
               q(ilo-1,jlo-1,khi+1) = 0.5d0 * &
                   (3*q(ilo-1,jlo-1,khi) - q(ilo-1,jlo-1,khi-1))
            end if
         end if    
      end if  
!
! ****************************************************************************
!
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP  .or. bc(1,1) .eq. FOEXTRAP) .and. &
          (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP .or. bc(2,2) .eq. FOEXTRAP) ) then
         if (jhi-2 .ge. js) then 
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = 0.5d0 * 0.125d0 * &
                   (15*q(ilo-1,jhi,k) - 10*q(ilo-1,jhi-1,k) + &
                   3*q(ilo-1,jhi-2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = 0.5d0 * 0.5d0 * &
                   (3*q(ilo-1,jhi,k) - q(ilo-1,jhi-1,k)) 
            end do
         end if

         if (ilo+2 .le. ie) then 
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + 0.5d0 * 0.125d0 * &
                   (15*q(ilo,jhi+1,k) - 10*q(ilo+1,jhi+1,k) + &
                   3*q(ilo+2,jhi+1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + 0.5d0 * 0.5d0 * &
                   (3*q(ilo,jhi+1,k) - q(ilo+1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP  .or. bc(3,1) .eq. FOEXTRAP) then
            if (klo+2 .le. ke) then
               q(ilo-1,jhi+1,klo-1) = 0.125d0 * ( &
                   (15*q(ilo-1,jhi+1,klo) - 10*q(ilo-1,jhi+1,klo+1) + &
                   3*q(ilo-1,jhi+1,klo+2)) )
            else
               q(ilo-1,jhi+1,klo-1) = 0.5d0 * &
                   (3*q(ilo-1,jhi+1,klo) - q(ilo-1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP  .or. bc(3,2) .eq. FOEXTRAP) then
            if (khi-2 .ge. ks) then
               q(ilo-1,jhi+1,khi+1) = 0.125d0 * ( &
                   (15*q(ilo-1,jhi+1,khi) - 10*q(ilo-1,jhi+1,khi-1) + &
                   3*q(ilo-1,jhi+1,khi-2)) )
            else
               q(ilo-1,jhi+1,khi+1) = 0.5d0 * &
                   (3*q(ilo-1,jhi+1,khi) - q(ilo-1,jhi+1,khi-1))
            end if
         end if
       end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP  .or. bc(1,2) .eq. FOEXTRAP) .and. &
          (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP  .or. bc(2,1) .eq. FOEXTRAP) ) then
         if (jlo+2 .le. je) then 
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = 0.5d0 * 0.125d0 * &
                   (15*q(ihi+1,jlo,k) - 10*q(ihi+1,jlo+1,k) + &
                   3*q(ihi+1,jlo+2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = 0.5d0 * 0.5d0 * &
                   (3*q(ihi+1,jlo,k) - q(ihi+1,jlo+1,k))
            end do
         end if

         if (ihi-2 .ge. is) then 
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + 0.5d0 * 0.125d0 * &
                   (15*q(ihi,jlo-1,k) - 10*q(ihi-1,jlo-1,k) + &
                   3*q(ihi-2,jlo-1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + 0.5d0 * 0.5d0 * &
                   (3*q(ihi,jlo-1,k) - q(ihi-1,jlo-1,k)) 
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP  .or. bc(3,1) .eq. FOEXTRAP) then
             if (klo+2 .le. ke) then
               q(ihi+1,jlo-1,klo-1) = 0.125d0 * &
                   (15*q(ihi+1,jlo-1,klo) - 10*q(ihi+1,jlo-1,klo+1) + &
                   3*q(ihi+1,jlo-1,klo+2))
            else
               q(ihi+1,jlo-1,klo-1) = 0.5d0 * &
                   (3*q(ihi+1,jlo-1,klo) - q(ihi+1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP  .or. bc(3,2) .eq. FOEXTRAP) then
            if (khi-2 .ge. ks) then
               q(ihi+1,jlo-1,khi+1) = 0.125d0 * & 
                   (15*q(ihi+1,jlo-1,khi) - 10*q(ihi+1,jlo-1,khi-1) + &
                   3*q(ihi+1,jlo-1,khi-2))
            else
               q(ihi+1,jlo-1,khi+1) = 0.5d0 * &
                   (3*q(ihi+1,jlo-1,khi) - q(ihi+1,jlo-1,khi-1))
            end if
         end if
         end if

!
! ****************************************************************************
!
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP  .or. bc(1,2) .eq. FOEXTRAP) .and. &
          (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP  .or. bc(2,2) .eq. FOEXTRAP) ) then
         if (jhi-2 .ge. js) then
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = 0.5d0 * 0.125d0 * &
                   (15*q(ihi+1,jhi,k) - 10*q(ihi+1,jhi-1,k) + &
                   3*q(ihi+1,jhi-2,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = 0.5d0 * 0.5d0 * &
                   (3*q(ihi+1,jhi,k) - q(ihi+1,jhi-1,k))
            end do
         end if

         if (ihi-2 .ge. is) then
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + 0.5d0 * 0.125d0 * &
                   (15*q(ihi,jhi+1,k) - 10*q(ihi-1,jhi+1,k) + &
                   3*q(ihi-2,jhi+1,k))
            end do
         else
            do k = q_l3,q_h3
               q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + 0.5d0 * 0.5d0 * &
                   (3*q(ihi,jhi+1,k) - q(ihi-1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP  .or. bc(3,1) .eq. FOEXTRAP) then
            if (klo+2 .le. ke) then
               q(ihi+1,jhi+1,klo-1) = 0.125d0 * &
                   (15*q(ihi+1,jhi+1,klo) - 10*q(ihi+1,jhi+1,klo+1) + &
                   3*q(ihi+1,jhi+1,klo+2))
            else
               q(ihi+1,jhi+1,klo-1) = 0.5d0 * &
                   (3*q(ihi+1,jhi+1,klo) - q(ihi+1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP  .or. bc(3,2) .eq. FOEXTRAP) then
            if (khi-2 .ge. ks) then
               q(ihi+1,jhi+1,khi+1) = 0.125d0 * &
                   (15*q(ihi+1,jhi+1,khi) - 10*q(ihi+1,jhi+1,khi-1) + &
                   3*q(ihi+1,jhi+1,khi-2))
            else
               q(ihi+1,jhi+1,khi+1) = 0.5d0 * &
                   (3*q(ihi+1,jhi+1,khi) - q(ihi+1,jhi+1,khi-1))
            end if
         end if
      endif
      
!
!    Next correct the i-k edges
!
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP  .or. bc(1,1) .eq. FOEXTRAP) .and. &
          (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP  .or. bc(3,1) .eq. FOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = 0.5d0 * 0.125d0 * &
                   (15*q(ilo-1,j,klo) - 10*q(ilo-1,j,klo+1) + &
                   3*q(ilo-1,j,klo+2))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = 0.5d0 * 0.5d0 * &
                   (3*q(ilo-1,j,klo) - q(ilo-1,j,klo+1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + 0.5d0 * 0.125d0 * &
                   (15*q(ilo,j,klo-1) - 10*q(ilo+1,j,klo-1) + &
                   3*q(ilo+2,j,klo-1))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + 0.5d0 * 0.5d0 * &
                   (3*q(ilo,j,klo-1) - q(ilo+1,j,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP  .or. bc(1,1) .eq. FOEXTRAP) .and. &
          (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP  .or. bc(3,2) .eq. FOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = 0.5d0 * 0.125d0 * &
                   (15*q(ilo-1,j,khi) - 10*q(ilo-1,j,khi-1) + &
                   3*q(ilo-1,j,khi-2))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = 0.5d0 * 0.5d0 * &
                   (3*q(ilo-1,j,khi) - q(ilo-1,j,khi-1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + 0.5d0 * 0.125d0 * &
                   (15*q(ilo,j,khi+1) - 10*q(ilo+1,j,khi+1) + &
                   3*q(ilo+2,j,khi+1))
            end do
         else
            do j = q_l2,q_h2
               q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + 0.5d0 * 0.5d0 * &
                   (3*q(ilo,j,khi+1) - q(ilo+1,j,khi+1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP  .or. bc(1,2) .eq. FOEXTRAP) .and. &
          (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP .or. bc(3,1) .eq. FOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = 0.5d0 * 0.125d0 * &
                   (15*q(ihi+1,j,klo) - 10*q(ihi+1,j,klo+1) + &
                   3*q(ihi+1,j,klo+2))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = 0.5d0 * 0.5d0 * &
                   (3*q(ihi+1,j,klo) - q(ihi+1,j,klo+1))
            end do
         end if

         if (ihi-2 .ge. is) then
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + 0.5d0 * 0.125d0 * &
                   (15*q(ihi,j,klo-1) - 10*q(ihi-1,j,klo-1) + &
                   3*q(ihi-2,j,klo-1))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + 0.5d0 * 0.5d0 * &
                   (3*q(ihi,j,klo-1) - q(ihi-1,j,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP  .or. bc(1,2) .eq. FOEXTRAP) .and. &
          (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP  .or. bc(3,2) .eq. FOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = 0.5d0 * 0.125d0 * &
                   (15*q(ihi+1,j,khi) - 10*q(ihi+1,j,khi-1) + &
                   3*q(ihi+1,j,khi-2))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = 0.5d0 * 0.5d0 * &
                   (3*q(ihi+1,j,khi) - q(ihi+1,j,khi-1))
            end do
         end if
         if (ihi-2 .ge. is) then
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + 0.5d0 * 0.125d0 * &
                   (15*q(ihi,j,khi+1) - 10*q(ihi-1,j,khi+1) + &
                   3*q(ihi-2,j,khi+1))
            end do
         else
            do j = q_l2,q_h2
               q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + 0.5d0 * 0.5d0 * &
                   (3*q(ihi,j,khi+1) - q(ihi-1,j,khi+1))
            end do
         end if
      end if

!
!    Next correct the j-k edges
!
      if ((nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP  .or. bc(2,1) .eq. FOEXTRAP) .and. &
          (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP .or. bc(3,1) .eq. FOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = 0.5d0 * 0.125d0 * &
                   (15*q(i,jlo-1,klo) - 10*q(i,jlo-1,klo+1) + &
                   3*q(i,jlo-1,klo+2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = 0.5d0 * 0.5d0 * &
                   (3*q(i,jlo-1,klo) - q(i,jlo-1,klo+1))
            end do
         end if
         if (jlo+2 .le. je) then
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + 0.5d0 * 0.125d0 * &
                   (15*q(i,jlo,klo-1) - 10*q(i,jlo+1,klo-1) + &
                   3*q(i,jlo+2,klo-1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + 0.5d0 * 0.5d0 * &
                   (3*q(i,jlo,klo-1) - q(i,jlo+1,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP  .or. bc(2,1) .eq. FOEXTRAP) .and. &
          (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP  .or. bc(3,2) .eq. FOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = 0.5d0 * 0.125d0 * &
                   (15*q(i,jlo-1,khi) - 10*q(i,jlo-1,khi-1) + &
                   3*q(i,jlo-1,khi-2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = 0.5d0 * 0.5d0 * &
                   (3*q(i,jlo-1,khi) - q(i,jlo-1,khi-1))
            end do
         end if

         if (jlo+2 .le. je) then
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + 0.5d0 * 0.125d0 * &
                   (15*q(i,jlo,khi+1) - 10*q(i,jlo+1,khi+1) + &
                   3*q(i,jlo+2,khi+1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + 0.5d0 * 0.5d0 * &
                   (3*q(i,jlo,khi+1) - q(i,jlo+1,khi+1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP  .or. bc(2,2) .eq. FOEXTRAP) .and. &
          (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP  .or. bc(3,1) .eq. FOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = 0.5d0 * 0.125d0 * &
                   (15*q(i,jhi+1,klo) - 10*q(i,jhi+1,klo+1) + &
                   3*q(i,jhi+1,klo+2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = 0.5d0 * 0.5d0 * &
                   (3*q(i,jhi+1,klo) - q(i,jhi+1,klo+1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + 0.5d0 * 0.125d0 * &
                   (15*q(i,jhi,klo-1) - 10*q(i,jhi-1,klo-1) + &
                   3*q(i,jhi-2,klo-1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + 0.5d0 * 0.5d0 * &
                   (3*q(i,jhi,klo-1) - q(i,jhi-1,klo-1))
            end do
         end if
      end if
!
! ****************************************************************************
!
      if ((ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP  .or. bc(2,2) .eq. FOEXTRAP) .and. &
          (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP  .or. bc(3,2) .eq. FOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = 0.5d0 * 0.125d0 * &
                   (15*q(i,jhi+1,khi) - 10*q(i,jhi+1,khi-1) + &
                   3*q(i,jhi+1,khi-2))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = 0.5d0 * 0.5d0 * &
                   (3*q(i,jhi+1,khi) - q(i,jhi+1,khi-1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + 0.5d0 * 0.125d0 * &
                   (15*q(i,jhi,khi+1) - 10*q(i,jhi-1,khi+1) + &
                   3*q(i,jhi-2,khi+1))
            end do
         else
            do i = q_l1,q_h1
               q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + 0.5d0 * 0.5d0 * &
                   (3*q(i,jhi,khi+1) - q(i,jhi-1,khi+1))
            end do
         end if
      end if
      
     end subroutine filfc

end module fc_fill_module
