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
  real(rt)         :: bcl, dx(1)
  real(rt)         :: dterm(DIMV(dtbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: h
  integer :: i
  h = dx(1)
  if (bct == LO_DIRICHLET) then
     if (cdir == 0) then
        !     Left face of grid
        i = reg_l1
        if (mask(i-1) > 0) then
           dterm(i) = d(i)*(er(i) - bcval(i-1))/(0.5e0_rt*h+bcl)
        endif
     else if (cdir == 1) then
        !     Right face of grid
        i = reg_h1
        if (mask(i+1) > 0) then
           dterm(i+1) = d(i+1)*(bcval(i+1)-er(i))/(0.5e0_rt*h+bcl)
        endif
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
  real(rt)         :: bcl, dx(1)
  real(rt)         :: dterm(DIMV(dtbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: d(DIMV(dbox))
  real(rt)         :: h
  integer :: i, bct
  h = dx(1)
  if (cdir == 0) then
     ! Left face of grid
     i = reg_l1
     if (mask(i-1) > 0) then
        if (bctype == -1) then
           bct = tf(i-1)
        else
           bct = bctype
        endif
        if (bct == LO_DIRICHLET) then
           dterm(i) = d(i)*(er(i) - bcval(i-1))/(0.5e0_rt*h+bcl)
        else if (bct == LO_NEUMANN .AND. bcval(i-1) == 0.e0_rt) then
           dterm(i) = 0.e0_rt
        else
           print *, "hdterm3: unsupported boundary type"
           stop
        endif
     endif
  else if (cdir == 1) then
     ! Right face of grid
     i = reg_h1
     if (mask(i+1) > 0) then
        if (bctype == -1) then
           bct = tf(i+1)
        else
           bct = bctype
        endif
        if (bct == LO_DIRICHLET) then
           dterm(i+1) = d(i+1)*(bcval(i+1)-er(i))/(0.5e0_rt*h+bcl)
        else if (bct == LO_NEUMANN .AND. bcval(i+1) == 0.e0_rt) then
           dterm(i+1) = 0.e0_rt
        else
           print *, "hbterm3: unsupported boundary type"
           stop
        endif
     endif
  else
     print *, "hdterm3: impossible face orientation"
  endif
end subroutine hdterm3

end module habec_module
