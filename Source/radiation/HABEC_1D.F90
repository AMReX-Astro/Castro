#include "AMReX_LO_BCTYPES.H"
#include "AMReX_ArrayLim.H"

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
  real(rt)         :: bcl, beta, dx(1)
  real(rt)         :: flux(DIMV(fbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: h, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i
  h = dx(1)
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
     if (mask(i-1) > 0) then
        flux(i) = b(i) * (bfv * bcval(i-1) - bfm * er(i))
        if (bho >= 1) then
           flux(i) = flux(i) - b(i) * bfm2 * er(i+1)
        endif
     endif
  else if (cdir == 1) then
     ! Right face of grid
     i = reg_h1
     if (mask(i+1) > 0) then
        flux(i+1) = -b(i+1) * (bfv * bcval(i+1) - bfm * er(i))
        if (bho >= 1) then
           flux(i+1) = flux(i+1) + b(i+1) * bfm2 * er(i-1)
        endif
     endif
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
  real(rt)         :: bcl, beta, dx(1), c
  real(rt)         :: flux(DIMV(fbox))
  real(rt)         :: er(DIMV(ebox))
  real(rt)         :: bcval(DIMV(bcv))
  integer :: mask(DIMV(msk))
  real(rt)         :: b(DIMV(bbox))
  real(rt)         :: spa(DIMV(spabox))
  real(rt)         :: r(1)
  real(rt)         :: h, bfm, bfv, r0
  real(rt)         :: bfm2, h2, th2
  integer :: i, bct
  h = dx(1)
  ! r is passed as an array but actually has only one element, which is
  ! appropriate for the face we are doing here.
  r0 = r(1)
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
           if (bho >= 1) then
              h2 = 0.5e0_rt * h
              th2 = 3.e0_rt * h2
              bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i)
              bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i)
              bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i)
           else
              bfv = beta / (0.5e0_rt * h + bcl) * b(i)
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
              bfm  =  1.5e0_rt * spa(i) * c * bfv
              bfm2 = -0.5e0_rt * spa(i) * c * bfv
           else
              bfm = spa(i) * c * bfv
           endif
        else
           print *, "hbflx3: unsupported boundary type"
           stop
        endif
        if (inhom == 0) then
           bfv = 0.e0_rt
        endif
        flux(i) = (bfv * bcval(i-1) - bfm * er(i))
        if (bho >= 1) then
           flux(i) = flux(i) - bfm2 * er(i+1)
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
           if (bho >= 1) then
              h2 = 0.5e0_rt * h
              th2 = 3.e0_rt * h2
              bfv = 2.e0_rt * beta * h / ((bcl + h2) * (bcl + th2)) * b(i+1)
              bfm = (beta / h) * (th2 - bcl) / (bcl + h2)  * b(i+1)
              bfm2 = (beta / h) * (bcl - h2) / (bcl + th2) * b(i+1)
           else
              bfv = beta / (0.5e0_rt * h + bcl) * b(i+1)
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
              bfm  =  1.5e0_rt * spa(i) * c * bfv
              bfm2 = -0.5e0_rt * spa(i) * c * bfv
           else
              bfm = spa(i) * c * bfv
           endif
        else
           print *, "hbflx3: unsupported boundary type"
           stop
        endif
        if (inhom == 0) then
           bfv = 0.e0_rt
        endif
        flux(i+1) = -(bfv * bcval(i+1) - bfm * er(i))
        if (bho >= 1) then
           flux(i+1) = flux(i+1) + bfm2 * er(i-1)
        endif
     endif
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
