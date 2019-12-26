subroutine ca_ext_src(lo, hi, &
                      old_state, os_lo, os_hi, &
                      new_state, ns_lo, ns_hi, &
                      src, src_lo, src_hi, &
                      problo, dx, time, dt) bind(C, name='ca_ext_src')

  use amrex_constants_module, only: THIRD, M_PI
  use meth_params_module, only : NVAR, NSRC, UEDEN, UEINT, URHO, UTEMP, UFS
  use network, only: network_species_index
  use amrex_fort_module, only : rt => amrex_real

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

  integer          :: i,j, k
  real(rt)         :: x, y

  real(rt)        :: H, ey, y_layer, H_0, L_x, pi
  real(rt)        :: H_max

  integer, save :: ic12, im24, io16

  logical, save :: firstCall = .true.

  if (firstCall) then
     ic12 = network_species_index("carbon-12")
     io16 = network_species_index("oxygen-16")
     im24 = network_species_index("magnesium-24")
     firstCall = .false.
  endif

  src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.e0_rt
  H_max = 0.e0_rt

  y_layer = 1.25e8_rt
  L_x = 2.5e8_rt
  H_0 = 2.5e16_rt

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        y = (dble(j)+0.5e0_rt)*dx(2) + problo(2)
        ey = exp(-(y-y_layer)*(y-y_layer)/1.e14_rt)

        do i = lo(1), hi(1)
           x = (dble(i)+0.5e0_rt)*dx(1) + problo(1)

           H = ey*(1.e0_rt + &
                .00625_rt * sin( 2*M_PI*x/L_x) &
              + .01875_rt * sin((6*M_PI*x/L_x) + M_PI/3.e0_rt) &
              + .01250_rt * sin((8*M_PI*x/L_x) + M_PI/5.e0_rt))

           ! Source terms
           src(i,j,k,UEDEN) = new_state(i,j,k,URHO) * H * 2.5e16_rt
           src(i,j,k,UEINT) = src(i,j,k,UEDEN)

        end do
     end do
  end do

end subroutine ca_ext_src
