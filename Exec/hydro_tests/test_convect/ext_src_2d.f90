subroutine ca_ext_src(lo,hi, &
                      old_state,old_state_l1,old_state_l2,old_state_h1,old_state_h2,&
                      new_state,new_state_l1,new_state_l2,new_state_h1,new_state_h2,&
                      src,src_l1,src_l2,src_h1,src_h2,problo,dx,time,dt)

  use amrex_constants_module, only: THIRD
  use meth_params_module, only : NVAR, UEDEN, UEINT, URHO, UTEMP, UFS
  use network, only: network_species_index

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer         , intent(in   ) :: lo(2),hi(2)
  integer         , intent(in   ) :: old_state_l1,old_state_l2,old_state_h1,old_state_h2
  integer         , intent(in   ) :: new_state_l1,new_state_l2,new_state_h1,new_state_h2
  integer         , intent(in   ) :: src_l1,src_l2,src_h1,src_h2
  real(rt)        , intent(in   ) :: old_state(old_state_l1:old_state_h1,old_state_l2:old_state_h2,NVAR)
  real(rt)        , intent(in   ) :: new_state(new_state_l1:new_state_h1,new_state_l2:new_state_h2,NVAR)
  real(rt)        , intent(  out) :: src(    src_l1:  src_h1,  src_l2:src_h2  ,NVAR)
  real(rt)        , intent(in   ) :: problo(2),dx(2),time,dt

  integer          :: i,j
  real(rt)         :: x,y

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

  src(lo(1):hi(1),lo(2):hi(2),:) = 0.e0_rt
  H_max = 0.e0_rt

  y_layer = 1.25e8_rt
  L_x = 2.5e8_rt
  pi = 3.1415926535897932384626433e0_rt
  H_0 = 2.5e16_rt

  do j = lo(2),hi(2)
      y = (dble(j)+0.5e0_rt)*dx(2)
      ey = exp(-(y-y_layer)*(y-y_layer)/1.e14_rt)
      do i = lo(1),hi(1)
         x =  (dble(i)+0.5e0_rt)*dx(1) 

         H = ey*(1.e0_rt + &
                .00625_rt * sin( 2*pi*x/L_x) &
              + .01875_rt * sin((6*pi*x/L_x) + pi/3.e0_rt) &
              + .01250_rt * sin((8*pi*x/L_x) + pi/5.e0_rt))

        ! Source terms
        src(i,j,UEDEN) = new_state(i,j,URHO) * H * 2.5e16_rt
        src(i,j,UEINT) = src(i,j,UEDEN)

        ! H_max = max(H_max, src(i,j,UEDEN))

      end do
  end do
  
end subroutine ca_ext_src
