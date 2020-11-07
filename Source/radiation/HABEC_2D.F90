#include <AMReX_LO_BCTYPES.H>
#include <AMReX_ArrayLim.H>


module habec_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

subroutine hbflx(flux, &
                 DIMS(fbox), &
                 er, DIMS(ebox), &
                 DIMS(reg), &
                 cdir, bct, bho, bcl, &
                 bcval, DIMS(bcv), &
                 mask, DIMS(msk), &
                 b, DIMS(bbox), &
                 beta, dx, inhom) bind(C, name="hbflx")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(fbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: cdir, bct, bho, inhom
  real(rt)         :: bcl, beta, dx(2)
  real(rt)         :: flux(DIMV(fbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: h, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i, j
  if (cdir == 0 .OR. cdir == 2) then
     h = dx(1)
  else
     h = dx(2)
  endif
  if (bct == LO_DIRICHLET) then
     if (bho >= 1) then
        h2 = 0.5e0_rt * h
        th2 = 3.e0_rt * h2
        bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2))
        bfm = (beta / h) * (th2 - bcl) / (bcl + h2)
        bfm2 = (beta / h) * (bcl - h2) / (bcl + th2)
     else
        bfv = beta / (0.5e0_rt * h + bcl)
        bfm = bfv
     endif
  else
     print *, "hbflx: unsupported boundary type"
     stop
  endif
  if (inhom == 0) then
     bfv = 0.e0_rt
  endif
  if (cdir == 0) then
     ! Left face of grid
     i = reg_l1
     do j = reg_l2, reg_h2
        if (mask(i-1,j) > 0) then
           flux(i,j) = b(i,j) * (bfv * bcval(i-1,j) - bfm * er(i,j))
           if (bho >= 1) then
              flux(i,j) = flux(i,j) - b(i,j) * bfm2 * er(i+1,j)
           endif
        endif
     enddo
  else if (cdir == 2) then
     ! Right face of grid
     i = reg_h1
     do j = reg_l2, reg_h2
        if (mask(i+1,j) > 0) then
           flux(i+1,j) = -b(i+1,j) * (bfv * bcval(i+1,j) - bfm * er(i,j))
           if (bho >= 1) then
              flux(i+1,j) = flux(i+1,j) + b(i+1,j) * bfm2 * er(i-1,j)
           endif
        endif
     enddo
  else if (cdir == 1) then
     ! Bottom face of grid
     j = reg_l2
     do i = reg_l1, reg_h1
        if (mask(i,j-1) > 0) then
           flux(i,j) = b(i,j) * (bfv * bcval(i,j-1) - bfm * er(i,j))
           if (bho >= 1) then
              flux(i,j) = flux(i,j) - b(i,j) * bfm2 * er(i,j+1)
           endif
        endif
     enddo
  else if (cdir == 3) then
     ! Top face of grid
     j = reg_h2
     do i = reg_l1, reg_h1
        if (mask(i,j+1) > 0) then
           flux(i,j+1) = -b(i,j+1) * (bfv * bcval(i,j+1) - bfm * er(i,j))
           if (bho >= 1) then
              flux(i,j+1) = flux(i,j+1) + b(i,j+1) * bfm2 * er(i,j-1)
           endif
        endif
     enddo
  else
     print *, "hbflx: impossible face orientation"
  endif
end subroutine hbflx

subroutine hbflx3(flux, &
                  DIMS(fbox), &
                  er, DIMS(ebox), &
                  DIMS(reg), &
                  cdir, bctype, tf, bho, bcl, &
                  bcval, DIMS(bcv), &
                  mask, DIMS(msk), &
                  b, DIMS(bbox), &
                  beta, dx, c, r, inhom, &
                  spa, DIMS(spabox)) bind(C, name="hbflx3")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(fbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(bbox)
  integer :: DIMDEC(spabox)
  integer :: cdir, bctype, tf(DIMV(bcv)), bho, inhom
  real(rt)         :: bcl, beta, dx(2), c
  real(rt)         :: flux(DIMV(fbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: spa(DIMV(spabox))
  real(rt)         :: r(reg_l1:reg_h1)
  real(rt)         :: h, bfm, bfv, r0
  real(rt)         :: bfm2, h2, th2
  integer :: i, j, bct
  if (cdir == 0 .OR. cdir == 2) then
     h = dx(1)
     ! For the left and right faces, r is constant and the array actually has
     ! only one element.  This following "r(reg_l1)" looks wrong, but what it
     ! really means is that the array dimensions are meaningless for this case.
     r0 = r(reg_l1)
  else
     h = dx(2)
  endif
  if (cdir == 0) then
     ! Left face of grid
     i = reg_l1
     do j = reg_l2, reg_h2
        if (mask(i-1,j) > 0) then
           if (bctype == -1) then
              bct = tf(i-1,j)
           else
              bct = bctype
           endif
           if (bct == LO_DIRICHLET) then
              if (bho >= 1) then
                 h2 = 0.5e0_rt * h
                 th2 = 3.e0_rt * h2
                 bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j)
                 bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j)
                 bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j)
              else
                 bfv = beta / (0.5e0_rt * h + bcl) * b(i,j)
                 bfm = bfv
              endif
           else if (bct == LO_NEUMANN) then
              bfv  = beta * r0
              bfm  = 0.e0_rt
              bfm2 = 0.e0_rt
           else if (bct == LO_MARSHAK) then
              bfv = 2.e0_rt * beta * r0
              if (bho >= 1) then
                 bfm  =  0.375e0_rt * c * bfv
                 bfm2 = -0.125e0_rt * c * bfv
              else
                 bfm = 0.25e0_rt * c * bfv
              endif
           else if (bct == LO_SANCHEZ_POMRANING) then
              bfv = 2.e0_rt * beta * r0
              if (bho >= 1) then
                 bfm  =  1.5e0_rt * spa(i,j) * c * bfv
                 bfm2 = -0.5e0_rt * spa(i,j) * c * bfv
              else
                 bfm = spa(i,j) * c * bfv
              endif
           else
              print *, "hbflx3: unsupported boundary type"
              stop
           endif
           if (inhom == 0) then
              bfv = 0.e0_rt
           endif
           flux(i,j) = (bfv * bcval(i-1,j) - bfm * er(i,j))
           if (bho >= 1) then
              flux(i,j) = flux(i,j) - bfm2 * er(i+1,j)
           endif
        endif
     enddo
  else if (cdir == 2) then
     ! Right face of grid
     i = reg_h1
     do j = reg_l2, reg_h2
        if (mask(i+1,j) > 0) then
           if (bctype == -1) then
              bct = tf(i+1,j)
           else
              bct = bctype
           endif
           if (bct == LO_DIRICHLET) then
              if (bho >= 1) then
                 h2 = 0.5e0_rt * h
                 th2 = 3.e0_rt * h2
                 bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i+1,j)
                 bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i+1,j)
                 bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i+1,j)
              else
                 bfv = beta / (0.5e0_rt * h + bcl) * b(i+1,j)
                 bfm = bfv
              endif
           else if (bct == LO_NEUMANN) then
              bfv  = beta * r0
              bfm  = 0.e0_rt
              bfm2 = 0.e0_rt
           else if (bct == LO_MARSHAK) then
              bfv = 2.e0_rt * beta * r0
              if (bho >= 1) then
                 bfm  =  0.375e0_rt * c * bfv
                 bfm2 = -0.125e0_rt * c * bfv
              else
                 bfm = 0.25e0_rt * c * bfv
              endif
           else if (bct == LO_SANCHEZ_POMRANING) then
              bfv = 2.e0_rt * beta * r0
              if (bho >= 1) then
                 bfm  =  1.5e0_rt * spa(i,j) * c * bfv
                 bfm2 = -0.5e0_rt * spa(i,j) * c * bfv
              else
                 bfm = spa(i,j) * c * bfv
              endif
           else
              print *, "hbflx3: unsupported boundary type"
              stop
           endif
           if (inhom == 0) then
              bfv = 0.e0_rt
           endif
           flux(i+1,j) = -(bfv * bcval(i+1,j) - bfm * er(i,j))
           if (bho >= 1) then
              flux(i+1,j) = flux(i+1,j) + bfm2 * er(i-1,j)
           endif
        endif
     enddo
  else if (cdir == 1) then
     ! Bottom face of grid
     j = reg_l2
     do i = reg_l1, reg_h1
        if (mask(i,j-1) > 0) then
           if (bctype == -1) then
              bct = tf(i,j-1)
           else
              bct = bctype
           endif
           if (bct == LO_DIRICHLET) then
              if (bho >= 1) then
                 h2 = 0.5e0_rt * h
                 th2 = 3.e0_rt * h2
                 bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j)
                 bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j)
                 bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j)
              else
                 bfv = beta / (0.5e0_rt * h + bcl) * b(i,j)
                 bfm = bfv
              endif
           else if (bct == LO_NEUMANN) then
              bfv  = beta * r(i)
              bfm  = 0.e0_rt
              bfm2 = 0.e0_rt
           else if (bct == LO_MARSHAK) then
              bfv = 2.e0_rt * beta * r(i)
              if (bho >= 1) then
                 bfm  =  0.375e0_rt * c * bfv
                 bfm2 = -0.125e0_rt * c * bfv
              else
                 bfm = 0.25e0_rt * c * bfv
              endif
           else if (bct == LO_SANCHEZ_POMRANING) then
              bfv = 2.e0_rt * beta * r(i)
              if (bho >= 1) then
                 bfm  =  1.5e0_rt * spa(i,j) * c * bfv
                 bfm2 = -0.5e0_rt * spa(i,j) * c * bfv
              else
                 bfm = spa(i,j) * c * bfv
              endif
           else
              print *, "hbflx3: unsupported boundary type"
              stop
           endif
           if (inhom == 0) then
              bfv = 0.e0_rt
           endif
           flux(i,j) = (bfv * bcval(i,j-1) - bfm * er(i,j))
           if (bho >= 1) then
              flux(i,j) = flux(i,j) - bfm2 * er(i,j+1)
           endif
        endif
     enddo
  else if (cdir == 3) then
     ! Top face of grid
     j = reg_h2
     do i = reg_l1, reg_h1
        if (mask(i,j+1) > 0) then
           if (bctype == -1) then
              bct = tf(i,j+1)
           else
              bct = bctype
           endif
           if (bct == LO_DIRICHLET) then
              if (bho >= 1) then
                 h2 = 0.5e0_rt * h
                 th2 = 3.e0_rt * h2
                 bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i,j+1)
                 bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i,j+1)
                 bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i,j+1)
              else
                 bfv = beta / (0.5e0_rt * h + bcl) * b(i,j+1)
                 bfm = bfv
              endif
           else if (bct == LO_NEUMANN) then
              bfv  = beta * r(i)
              bfm  = 0.e0_rt
              bfm2 = 0.e0_rt
           else if (bct == LO_MARSHAK) then
              bfv = 2.e0_rt * beta * r(i)
              if (bho >= 1) then
                 bfm  =  0.375e0_rt * c * bfv
                 bfm2 = -0.125e0_rt * c * bfv
              else
                 bfm = 0.25e0_rt * c * bfv
              endif
           else if (bct == LO_SANCHEZ_POMRANING) then
              bfv = 2.e0_rt * beta * r(i)
              if (bho >= 1) then
                 bfm  =  1.5e0_rt * spa(i,j) * c * bfv
                 bfm2 = -0.5e0_rt * spa(i,j) * c * bfv
              else
                 bfm = spa(i,j) * c * bfv
              endif
           else
              print *, "hbflx3: unsupported boundary type"
              stop
           endif
           if (inhom == 0) then
              bfv = 0.e0_rt
           endif
           flux(i,j+1) = -(bfv * bcval(i,j+1) - bfm * er(i,j))
           if (bho >= 1) then
              flux(i,j+1) = flux(i,j+1) + bfm2 * er(i,j-1)
           endif
        endif
     enddo
  else
     print *, "hbflx3: impossible face orientation"
  endif
end subroutine hbflx3

subroutine hdterm(dterm, &
                  DIMS(dtbox), &
                  er, DIMS(ebox), &
                  DIMS(reg), &
                  cdir, bct, bcl, &
                  bcval, DIMS(bcv), &
                  mask, DIMS(msk), &
                  d, DIMS(dbox), &
                  dx) bind(C, name="hdterm")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(dtbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(dbox)
  integer :: cdir, bct
  real(rt)         :: bcl, dx(2)
  real(rt)         :: dterm(DIMV(dtbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: h
  integer :: i, j
  if (cdir == 0 .OR. cdir == 2) then
     h = dx(1)
  else
     h = dx(2)
  endif
  if (bct == LO_DIRICHLET) then
     if (cdir == 0) then
        !     Left face of grid
        i = reg_l1
        do j = reg_l2, reg_h2
           if (mask(i-1,j) > 0) then
              dterm(i,j) = d(i,j)*(er(i,j)-bcval(i-1,j))/(0.5e0_rt*h+bcl)
           endif
        enddo
     else if (cdir == 2) then
        !     Right face of grid
        i = reg_h1
        do j = reg_l2, reg_h2
           if (mask(i+1,j) > 0) then
              dterm(i+1,j) = d(i+1,j)*(bcval(i+1,j)-er(i,j))/(0.5e0_rt*h+bcl)
           endif
        enddo
     else if (cdir == 1) then
        !     Bottom face of grid
        j = reg_l2
        do i = reg_l1, reg_h1
           if (mask(i,j-1) > 0) then
              dterm(i,j) = d(i,j)*(er(i,j)-bcval(i,j-1))/(0.5e0_rt*h+bcl)
           endif
        enddo
     else if (cdir == 3) then
        !     Top face of grid
        j = reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j+1) > 0) then
              dterm(i,j+1) = d(i,j+1)*(bcval(i,j+1)-er(i,j))/(0.5e0_rt*h+bcl)
           endif
        enddo
     else
        print *, "hdterm: impossible face orientation"
     endif
  else
     print *, "hdterm: unsupported boundary type"
     stop
  endif
end subroutine hdterm

subroutine hdterm3(dterm, &
                   DIMS(dtbox), &
                   er, DIMS(ebox), &
                   DIMS(reg), &
                   cdir, bctype, tf, bcl, &
                   bcval, DIMS(bcv), &
                   mask, DIMS(msk), &
                   d, DIMS(dbox), &
                   dx) bind(C, name="hdterm3")

  use amrex_fort_module, only : rt => amrex_real
  integer :: DIMDEC(dtbox)
  integer :: DIMDEC(ebox)
  integer :: DIMDEC(reg)
  integer :: DIMDEC(bcv)
  integer :: DIMDEC(msk)
  integer :: DIMDEC(dbox)
  integer :: cdir, bctype, tf(DIMV(bcv))
  real(rt)         :: bcl, dx(2)
  real(rt)         :: dterm(DIMV(dtbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: h
  integer :: i, j, bct
  if (cdir == 0 .OR. cdir == 2) then
     h = dx(1)
  else
     h = dx(2)
  endif
  if (cdir == 0) then
     ! Left face of grid
     i = reg_l1
     do j = reg_l2, reg_h2
        if (mask(i-1,j) > 0) then
           if (bctype == -1) then
              bct = tf(i-1,j)
           else
              bct = bctype
           endif
           if (bct == LO_DIRICHLET) then
              dterm(i,j) = d(i,j)*(er(i,j)-bcval(i-1,j))/(0.5e0_rt*h+bcl)
           else if (bct == LO_NEUMANN .AND. bcval(i-1,j) == 0.e0_rt) then
              dterm(i,j) = 0.e0_rt
           else
              print *, "hdterm3: unsupported boundary type"
              stop
           endif
        endif
     enddo
  else if (cdir == 2) then
     ! Right face of grid
     i = reg_h1
     do j = reg_l2, reg_h2
        if (mask(i+1,j) > 0) then
           if (bctype == -1) then
              bct = tf(i+1,j)
           else
              bct = bctype
           endif
           if (bct == LO_DIRICHLET) then
              dterm(i+1,j) = d(i+1,j)*(bcval(i+1,j)-er(i,j))/(0.5e0_rt*h+bcl)
           else if (bct == LO_NEUMANN .AND. bcval(i+1,j) == 0.e0_rt) then
              dterm(i+1,j) = 0.e0_rt
           else
              print *, "hdterm3: unsupported boundary type"
              stop
           endif
        endif
     enddo
  else if (cdir == 1) then
     ! Bottom face of grid
     j = reg_l2
     do i = reg_l1, reg_h1
        if (mask(i,j-1) > 0) then
           if (bctype == -1) then
              bct = tf(i,j-1)
           else
              bct = bctype
           endif
           if (bct == LO_DIRICHLET) then
              dterm(i,j) = d(i,j)*(er(i,j)-bcval(i,j-1))/(0.5e0_rt*h+bcl)
           else if (bct == LO_NEUMANN .AND. bcval(i,j-1) == 0.e0_rt) then
              dterm(i,j) = 0.e0_rt
           else
              print *, "hdterm3: unsupported boundary type"
              stop
           endif
        endif
     enddo
  else if (cdir == 3) then
     ! Top face of grid
     j = reg_h2
     do i = reg_l1, reg_h1
        if (mask(i,j+1) > 0) then
           if (bctype == -1) then
              bct = tf(i,j+1)
           else
              bct = bctype
           endif
           if (bct == LO_DIRICHLET) then
              dterm(i,j+1) = d(i,j+1)*(bcval(i,j+1)-er(i,j))/(0.5e0_rt*h+bcl)
           else if (bct == LO_NEUMANN .AND. bcval(i,j+1) == 0.e0_rt) then
              dterm(i,j+1) = 0.e0_rt
           else
              print *, "hdterm3: unsupported boundary type"
              stop
           endif
        endif
     enddo
  else
     print *, "hdterm3: impossible face orientation"
  endif
end subroutine hdterm3

end module habec_module
