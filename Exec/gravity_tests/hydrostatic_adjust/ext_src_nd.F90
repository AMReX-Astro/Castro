
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

  use amrex_constants_module
  use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, UFS
  use probdata_module   , only : xmin, ymin, zmin, &
                                 heating_time, heating_rad, &
                                 heating_peak, heating_sigma, prob_type
  use prob_params_module, only: center
  use network
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: os_lo(3), os_hi(3)
  integer,  intent(in   ) :: ns_lo(3), ns_hi(3)
  integer,  intent(in   ) :: src_lo(3), src_hi(3)
  real(rt), intent(in   ) :: old_state(os_lo(1):os_hi(1),os_lo(2):os_hi(2),os_lo(3):os_hi(3),NVAR)
  real(rt), intent(in   ) :: new_state(ns_lo(1):ns_hi(1),ns_lo(2):ns_hi(2),ns_lo(3):ns_hi(3),NVAR)
  real(rt), intent(inout) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
  real(rt), intent(in   ) :: problo(3), dx(3)
  real(rt), intent(in   ), value :: time, dt

  integer  :: i, j, k
  real(rt) :: x, y, z, dist, r_0, H_0, W_0, Hext, t_stop

  integer :: ihe4_p

  r_0 = heating_rad

  H_0 = heating_peak
  W_0 = heating_sigma

  t_stop = heating_time

  if (prob_type == 1) then

     ! For heating at the center
     if (time <= t_stop) then
        if (lo(1) == 0 .and. lo(2) == 0 .and. lo(3) == 0) print *,'TIME vs TSTOP ',time, t_stop
        do k = lo(3), hi(3)
           z = zmin + (dble(k) + HALF)*dx(3) - center(3)

           do j = lo(2), hi(2)
              y = ymin + (dble(j) + HALF)*dx(2) - center(2)

              do i = lo(1), hi(1)
                 x = xmin + (dble(i) + HALF)*dx(1) - center(1)

                 dist = sqrt(x*x + y*y + z*z)

                 Hext = H_0*exp(-((dist - r_0)**2)/W_0**2)

                 src(i,j,k,UEINT) = old_state(i,j,k,URHO)*Hext
                 src(i,j,k,UEDEN) = old_state(i,j,k,URHO)*Hext

              end do
           end do
        end do
     endif

  else if (prob_type == 3) then

       ihe4_p = network_species_index("helium-4")

       ! sub-chandra heating -- modulate by He
       if (time <= t_stop) then
          if (lo(1) == 0 .and. lo(2) == 0 .and. lo(3) == 0) print *,'TIME vs TSTOP ',time, t_stop

        do k = lo(3), hi(3)
           z = zmin + (dble(k) + HALF)*dx(3) - center(3)

           do j = lo(2), hi(2)
              y = ymin + (dble(j) + HALF)*dx(2) - center(2)

              do i = lo(1), hi(1)
                 x = xmin + (dble(i) + HALF)*dx(1) - center(1)

                 dist = sqrt(x*x + y*y + z*z)

                 Hext = H_0*exp(-((dist - r_0)**2)/W_0**2) * &
                      old_state(i,j,k,UFS-1+ihe4_p)/old_state(i,j,k,URHO)

                 src(i,j,k,UEINT) = old_state(i,j,k,URHO)*Hext
                 src(i,j,k,UEDEN) = old_state(i,j,k,URHO)*Hext

              end do
           end do
        end do

     endif

  end if

end subroutine ca_ext_src
