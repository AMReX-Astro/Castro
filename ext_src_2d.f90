subroutine ca_ext_src(lo,hi, &
                      old_state,old_state_l1,old_state_l2,old_state_h1,old_state_h2,&
                      new_state,new_state_l1,new_state_l2,new_state_h1,new_state_h2,&
                      src,src_l1,src_l2,src_h1,src_h2,problo,dx,time,dt)

  use bl_types, only: dp_t
  use bl_constants_module, only: THIRD
  use meth_params_module, only : NVAR, UEDEN, UEINT, URHO, UTEMP, UFS
  use network, only: network_species_index

  implicit none
  integer         , intent(in   ) :: lo(2),hi(2)
  integer         , intent(in   ) :: old_state_l1,old_state_l2,old_state_h1,old_state_h2
  integer         , intent(in   ) :: new_state_l1,new_state_l2,new_state_h1,new_state_h2
  integer         , intent(in   ) :: src_l1,src_l2,src_h1,src_h2
  real(kind=dp_t) , intent(in   ) :: old_state(old_state_l1:old_state_h1,old_state_l2:old_state_h2,NVAR)
  real(kind=dp_t) , intent(in   ) :: new_state(new_state_l1:new_state_h1,new_state_l2:new_state_h2,NVAR)
  real(kind=dp_t) , intent(  out) :: src(    src_l1:  src_h1,  src_l2:src_h2  ,NVAR)
  real(kind=dp_t) , intent(in   ) :: problo(2),dx(2),time,dt

  integer          :: i,j
  double precision :: x,y

  real(kind=dp_t) :: H, ey, y_layer, H_0, L_x, pi
  real(kind=dp_t) :: H_max

  integer, save :: ic12, im24, io16

  logical, save :: firstCall = .true.

  if (firstCall) then
     ic12 = network_species_index("carbon-12")
     io16 = network_species_index("oxygen-16")
     im24 = network_species_index("magnesium-24")
     firstCall = .false.
  endif

  src(lo(1):hi(1),lo(2):hi(2),:) = 0.d0
  H_max = 0.d0

  y_layer = 1.25d8
  L_x = 2.5d8
  pi = 3.1415926535897932384626433d0
  H_0 = 2.5d16

  do j = lo(2),hi(2)
      y = (dble(j)+0.5d0)*dx(2)
      ey = exp(-(y-y_layer)*(y-y_layer)/1.d14)
      do i = lo(1),hi(1)
         x =  (dble(i)+0.5d0)*dx(1) 

         H = ey*(1.d0 + &
                .00625_dp_t * sin( 2*pi*x/L_x) &
              + .01875_dp_t * sin((6*pi*x/L_x) + pi/3.d0) &
              + .01250_dp_t * sin((8*pi*x/L_x) + pi/5.d0))

        ! Source terms
        src(i,j,UEDEN) = new_state(i,j,URHO) * H * 2.5d16
        src(i,j,UEINT) = src(i,j,UEDEN)

        ! H_max = max(H_max, src(i,j,UEDEN))

      end do
  end do
  
end subroutine ca_ext_src
