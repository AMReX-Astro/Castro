! arithmetic average, geometrically correct(?) but underestimates surface flux

#define tiny 1.d-50

#define KAVG0(a,b) (0.5d0 * (a + b + tiny))

! define KAVG(a,b,d) KAVG0(a,b)

! harmonic average, overestimates surface flux

#define KAVG1(a,b) ((2.d0 * a * b) / (a + b + tiny) + tiny)

! define KAVG(a,b,d) KAVG1(a,b)

! chooses arithmetic where optically thin, harmonic where optically thick,
! surface flux approximation at a thick/thin boundary

#define KAVG2(a,b,d) min(KAVG0(a,b), max(KAVG1(a,b), 4.d0 / (3.d0 * d)))

! define KAVG(a,b,d) KAVG2(a,b,d)

#define KAVG(a,b,d) kavg(a,b,d,-1)

real*8 function kavg(a, b, d, iopt)
  implicit none
  real*8 :: a, b, d
  integer :: iopt
  integer, save :: opt=100
  if (iopt >= 0) then
     opt = iopt
     if (opt > 2) then
        print *, "Fortran KAVG: invalid averaging option"
     endif
     return
  endif
  if (opt == 0) then
     kavg = KAVG0(a,b)
  else if (opt == 1) then
     kavg = KAVG1(a,b)
  else
     kavg = KAVG2(a,b,d)
  endif
end function kavg

