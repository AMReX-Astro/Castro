subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use probdata_module
  use amrex_error_module
  use amrex_constants_module
  use riemann_module
  use meth_params_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt) :: problo(2), probhi(2)

  real(rt), pointer :: ql(:,:,:,:), qr(:,:,:,:)
  real(rt), pointer :: qint(:,:,:,:), qaux(:,:,:,:)

  integer :: lo(3), hi(3), loa(3), hia(3)
  integer :: idir

  integer :: untin, i

  namelist /fortin/ rho_l, u_l, p_l, re_l, gc_l, &
                    rho_r, u_r, p_r, re_r, gc_r, &
                    cav_s, smallc_s, idir


  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) then
     call amrex_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Set namelist defaults

  rho_l = 1.0
  u_l = 1.0
  p_l = 1.0
  re_l = 1.0
  gc_l = 4.0/3.0

  rho_r = 1.0
  u_r = 1.0
  p_r = 1.0
  re_r = 1.0
  gc_r = 4.0e0_rt/3.0e0_rt

  cav_s = 1.0
  smallc_s = 1.e-10_rt

  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(unit=untin)

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
  call amrex_error("done with Riemann")

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
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use probdata_module
  use amrex_constants_module, only: M_PI, FOUR3RD
  use meth_params_module , only: NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS
  use prob_params_module, only : center
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  real(rt)         :: xlo(2), xhi(2), time, delta(2)
  real(rt)         :: state(state_l1:state_h1,state_l2:state_h2,NVAR)


end subroutine ca_initdata
