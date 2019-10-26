subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use probdata_module
  use castro_error_module
  use amrex_constants_module
  use riemann_module
  use meth_params_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  real(rt), pointer :: ql(:,:,:,:), qr(:,:,:,:)
  real(rt), pointer :: qint(:,:,:,:), qaux(:,:,:,:)

  integer :: lo(3), hi(3), loa(3), hia(3)

  call probdata_init(name, namlen)

  ! call the Riemann solver
  lo(:) = [1, 1, 0]
  hi(:) = [1, 1, 0]

  loa(:) = [0, 1, 0]
  hia(:) = [1, 1, 0]

  call bl_allocate(ql, lo, hi, NQ)
  call bl_allocate(qr, lo, hi, NQ)

  ! Riemann indexes i-1 in qaux
  call bl_allocate(qaux, loa, hia, NQAUX)
  call bl_allocate(qint, lo, hi, NQ)

  ! set the Riemann arrays
  ql(:,:,:,:) = ZERO
  ql(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), QRHO) = rho_l
  ql(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), QU) = u_l
  ql(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), QPRES) = p_l
  ql(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), QREINT) = re_l

  qaux(:,:,:,:) = ZERO

  qaux(lo(1)-1,:,:,QGAMC) = gc_l

  qr(:,:,:,:) = ZERO
  qr(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), QRHO) = rho_r
  qr(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), QU) = u_r
  qr(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), QPRES) = p_r
  qr(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), QREINT) = re_r

  qaux(lo(1),:,:,QGAMC) = gc_r

  qaux(:,:,:,QC) = cav_s

  ! call the Riemann solver
  idir = 1

  call riemanncg(ql, lo, hi, &
                 qr, lo, hi, 1, 1, &
                 qaux, loa, hia, &
                 qint, lo, hi, &
                 idir, lo, hi, &
                 [lo(1)-1, lo(2)-1, 0], [hi(1)+1, hi(2)+1, 0])

  ! we're done -- abort the code
  call castro_error("done with Riemann")

end subroutine amrex_probinit


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       delta, xlo, xhi)

  use meth_params_module , only: NVAR
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

end subroutine ca_initdata
