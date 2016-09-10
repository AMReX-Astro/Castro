! reference: R. J. Purser (J. of Clim. and Apld. Meteorology, 1987)

module filter_module

  ! 3-point filter: T=1, R=S=0
  double precision, dimension(0:1), parameter :: ff1 = &
       (/ 0.5d0, 0.25d0 /)
  ! For boundary cell
  double precision, dimension(0:1), parameter :: ff1b = &
       (/ 0.75d0, 0.25d0 /)
  

  ! 5-point filter: T=2, R+S+1=T
  ! S = 0, 1
  double precision, dimension(0:2,0:1), parameter :: ff2 = reshape( &
       [ 0.625d0, 0.25d0, -0.0625d0, &
       & 0.375d0, 0.25d0,  0.0625d0 ], &
       [3,2] )
  ! For first bondary cell
  double precision, dimension(0:2), parameter :: ff2b0 = & 
       (/ 17.d0/16.d0, -2.d0/16.d0, 1.d0/16.d0 /)
  ! for second boundary cell
  double precision, dimension(-1:2), parameter :: ff2b1 = & 
       (/ -2.d0/16.d0, 21.d0/16.d0, -4.d0/16.d0, 1.d0/16.d0 /)


  ! 7-point filter: T=3, R+S+1=T
  ! S = 0, 1, 2
  double precision, dimension(0:3,0:2), parameter :: ff3 = reshape( &
       [ 44.d0/64.d0, 15.d0/64.d0, -6.d0/64.d0,  1.d0/64.d0, &
       & 32.d0/64.d0, 18.d0/64.d0,  0.d0      , -2.d0/64.d0, &
       & 20.d0/64.d0, 15.d0/64.d0,  6.d0/64.d0,  1.d0/64.d0 ], &
       [4,3] ) 
  ! for boundary cells
  double precision, dimension(0:3), parameter :: ff3b0 = & 
       (/ 63.d0/64.d0, 3.d0/64.d0, -3.d0/64.d0, 1.d0/64.d0 /)
  double precision, dimension(-1:3), parameter :: ff3b1 = & 
       (/ 3.d0/64.d0, &
       54.d0/64.d0, 12.d0/64.d0, -6.d0/64.d0, 1.d0/64.d0 /)
  double precision, dimension(-2:3), parameter :: ff3b2 = & 
       (/ -3.d0/64.d0, 12.d0/64.d0, &
       45.d0/64.d0, 15.d0/64.d0, -6.d0/64.d0, 1.d0/64.d0 /)
  

  ! 9-point filter: T=4, R+S+1=T
  ! S = 0, 1, 2, 3
  double precision, dimension(0:4,0:3), parameter :: ff4 = reshape( &
       [ 186.d0/256.d0, 56.d0/256.d0, -28.d0/256.d0,  8.d0/256.d0, -1.d0/256.d0, &
       & 146.d0/256.d0, 72.d0/256.d0, -12.d0/256.d0, -8.d0/256.d0,  3.d0/256.d0, &
       & 110.d0/256.d0, 72.d0/256.d0,  12.d0/256.d0, -8.d0/256.d0, -3.d0/256.d0, &
       &  70.d0/256.d0, 56.d0/256.d0,  28.d0/256.d0,  8.d0/256.d0,  1.d0/256.d0 ], &
       [5,4] )
  ! for boundary cells
  double precision, dimension(0:4), parameter :: ff4b0 = & 
       (/ 257.d0/256.d0, -4.d0/256.d0, 6.d0/256.d0, -4.d0/256.d0, 1.d0/256.d0 /)
  double precision, dimension(-1:4), parameter :: ff4b1 = & 
       (/ -4.d0/256.d0,  &
       273.d0/256.d0, -28.d0/256.d0, 22.d0/256.d0, -8.d0/256.d0, 1.d0/256.d0 /)
  double precision, dimension(-2:4), parameter :: ff4b2 = & 
       (/ 6.d0/256.d0, -28.d0/256.d0, &
       309.d0/256.d0, -52.d0/256.d0, 28.d0/256.d0, -8.d0/256.d0, 1.d0/256.d0 /)
  double precision, dimension(-3:4), parameter :: ff4b3 = & 
       (/ -4.d0/256.d0, 22.d0/256.d0, -52.d0/256.d0, &
       325.d0/256.d0, -56.d0/256.d0, 28.d0/256.d0, -8.d0/256.d0, 1.d0/256.d0 /)


end module filter_module
