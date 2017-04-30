subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: init, namlen
  integer          :: name(namlen)
  real(rt)         :: problo(2), probhi(2)

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
! ::: nvar      => number of state components.
! ::: state     <= scalar array
! ::: dx        => cell size
! ::: xlo, xhi  => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------

subroutine ca_initdata(level,time,lo,hi,nvar, &
                       state,state_l1,state_l2,state_h1,state_h2, &
                       dx,xlo,xhi)

  use bl_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: level, nvar
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  real(rt)         :: xlo(2), xhi(2), time, dx(2)
  real(rt)         :: state(state_l1:state_h1,state_l2:state_h2,nvar)

  ! Remove this call if you're defining your own problem; it is here to 
  ! ensure that you cannot run Castro if you haven't got your own copy of this function.

  call bl_error("Prob_2d.f90 has not been defined for this problem!")

end subroutine ca_initdata

