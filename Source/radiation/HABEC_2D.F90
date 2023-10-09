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
