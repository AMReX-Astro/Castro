
!> @brief Compute the external sources for all the conservative equations.
!!
!! This is called twice in the evolution:
!!
!! First, for the predictor, it is called with (old, old) states.
!!
!! This is also used in the first pass of the conservative update
!! (adding dt * S there).
!!
!! Next we correct the source terms in the conservative update to
!! time-center them.  Here we call ext_src(old, new), and then
!! in time_center_source_terms we subtract off 1/2 of the first S
!! and add 1/2 of the new S.
!!
!! Therefore, to get a properly time-centered source, generally
!! speaking, you always want to use the "new" state here.  That
!! will be the time n state in the first call and the n+1 in the
!! second call.
!!
subroutine ca_ext_src(lo, hi, &
                      old_state, os_lo, os_hi, &
                      new_state, ns_lo, ns_hi, &
                      src, src_lo, src_hi, &
                      problo, dx, time, dt) bind(C, name='ca_ext_src')

  use amrex_constants_module, only: HALF
  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR, NSRC, URHO, UMX, UEDEN
  use prob_params_module, only: center, probhi
  use probdata_module, only: grav_acceleration, center_T

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: os_lo(3), os_hi(3)
  integer,  intent(in   ) :: ns_lo(3), ns_hi(3)
  integer,  intent(in   ) :: src_lo(3), src_hi(3)
  real(rt), intent(in   ) :: old_state(os_lo(1):os_hi(1),os_lo(2):os_hi(2),os_lo(3):os_hi(3),NVAR)
  real(rt), intent(in   ) :: new_state(ns_lo(1):ns_hi(1),ns_lo(2):ns_hi(2),ns_lo(3):ns_hi(3),NVAR)
  real(rt), intent(inout) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NSRC)
  real(rt), intent(in   ) :: problo(3), dx(3)
  real(rt), intent(in   ), value :: time, dt

  integer  :: i, j, k
  real(rt) :: x, g, c_T

  !$gpu

  ! Add a mock gravitational acceleration which points to the center
  ! with uniform magnitude on either side of the center.

  c_T = problo(1) + center_T * (probhi(1) - problo(1))

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           x = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           if (x > c_T) then
              g = -grav_acceleration
           else
              g = grav_acceleration
           end if

           src(i,j,k,UMX) = new_state(i,j,k,URHO) * g

           ! Energy source is v . momentum source

           src(i,j,k,UEDEN) = new_state(i,j,k,UMX) * g

        end do
     end do
  end do

end subroutine ca_ext_src
