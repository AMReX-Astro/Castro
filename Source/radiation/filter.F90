! reference: R. J. Purser (J. of Clim. and Apld. Meteorology, 1987)

module filter_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  ! 3-point filter: T=1, R=S=0
  real(rt), dimension(0:1), parameter :: ff1 = [0.5e0_rt, 0.25e0_rt]

  ! For boundary cell
  real(rt), dimension(0:1), parameter :: ff1b = [0.75e0_rt, 0.25e0_rt]

  ! 5-point filter: T=2, R+S+1=T
  ! S = 0, 1
  real(rt), dimension(0:2,0:1), parameter :: ff2 = reshape( &
                                                   [ 0.625e0_rt, 0.25e0_rt, &
                                                    -0.0625e0_rt, 0.375e0_rt, &
                                                     0.25e0_rt,  0.0625e0_rt ], &
                                                     [3,2] )

  ! For first boundary cell
  real(rt), dimension(0:2), parameter :: ff2b0 = [17.e0_rt/16.e0_rt, &
                                                  -2.e0_rt/16.e0_rt, &
                                                  1.e0_rt/16.e0_rt]

  ! For second boundary cell
  real(rt), dimension(-1:2), parameter :: ff2b1 = [-2.e0_rt/16.e0_rt, &
                                                   21.e0_rt/16.e0_rt, &
                                                   -4.e0_rt/16.e0_rt, &
                                                   1.e0_rt/16.e0_rt]

  ! 7-point filter: T=3, R+S+1=T
  ! S = 0, 1, 2
  real(rt), dimension(0:3,0:2), parameter :: ff3 = reshape( &
                                                   [44.e0_rt/64.e0_rt, &
                                                    15.e0_rt/64.e0_rt, &
                                                    -6.e0_rt/64.e0_rt, &
                                                    1.e0_rt/64.e0_rt, &
                                                    32.e0_rt/64.e0_rt, &
                                                    18.e0_rt/64.e0_rt, &
                                                    0.e0_rt, -2.e0_rt/64.e0_rt, &
                                                    & 20.e0_rt/64.e0_rt, &
                                                    15.e0_rt/64.e0_rt, &
                                                    6.e0_rt/64.e0_rt, &
                                                    1.e0_rt/64.e0_rt ], &
                                                    [4,3] )
  ! For boundary cells
  real(rt), dimension(0:3), parameter :: ff3b0 = [63.e0_rt/64.e0_rt, &
                                                  3.e0_rt/64.e0_rt, &
                                                  -3.e0_rt/64.e0_rt, &
                                                  1.e0_rt/64.e0_rt]

  real(rt), dimension(-1:3), parameter :: ff3b1 = [3.e0_rt/64.e0_rt, &
                                                   54.e0_rt/64.e0_rt, &
                                                   12.e0_rt/64.e0_rt, &
                                                   -6.e0_rt/64.e0_rt, &
                                                   1.e0_rt/64.e0_rt]

  real(rt), dimension(-2:3), parameter :: ff3b2 = [-3.e0_rt/64.e0_rt, &
                                                   12.e0_rt/64.e0_rt, &
                                                   45.e0_rt/64.e0_rt, &
                                                   15.e0_rt/64.e0_rt, &
                                                   -6.e0_rt/64.e0_rt, &
                                                   1.e0_rt/64.e0_rt]

  ! 9-point filter: T=4, R+S+1=T
  ! S = 0, 1, 2, 3
  real(rt), dimension(0:4,0:3), parameter :: ff4 = reshape( &
                                                   [ 186.e0_rt/256.e0_rt, &
                                                     56.e0_rt/256.e0_rt, &
                                                    -28.e0_rt/256.e0_rt, &
                                                     8.e0_rt/256.e0_rt, &
                                                    -1.e0_rt/256.e0_rt, &
                                                     146.e0_rt/256.e0_rt, &
                                                     72.e0_rt/256.e0_rt, &
                                                    -12.e0_rt/256.e0_rt, &
                                                    -8.e0_rt/256.e0_rt, &
                                                     3.e0_rt/256.e0_rt, &
                                                     110.e0_rt/256.e0_rt, &
                                                     72.e0_rt/256.e0_rt, &
                                                     12.e0_rt/256.e0_rt, &
                                                    -8.e0_rt/256.e0_rt, &
                                                    -3.e0_rt/256.e0_rt, &
                                                     70.e0_rt/256.e0_rt, &
                                                     56.e0_rt/256.e0_rt, &
                                                     28.e0_rt/256.e0_rt, &
                                                     8.e0_rt/256.e0_rt, &
                                                     1.e0_rt/256.e0_rt ], &
                                                     [5,4] )

  ! For boundary cells
  real(rt), dimension(0:4), parameter :: ff4b0 = [257.e0_rt/256.e0_rt, &
                                                  -4.e0_rt/256.e0_rt, &
                                                  6.e0_rt/256.e0_rt, &
                                                  -4.e0_rt/256.e0_rt, &
                                                  1.e0_rt/256.e0_rt]

  real(rt), dimension(-1:4), parameter :: ff4b1 = [-4.e0_rt/256.e0_rt, &
                                                   273.e0_rt/256.e0_rt, &
                                                   -28.e0_rt/256.e0_rt, &
                                                   22.e0_rt/256.e0_rt, &
                                                   -8.e0_rt/256.e0_rt, &
                                                   1.e0_rt/256.e0_rt]

  real(rt), dimension(-2:4), parameter :: ff4b2 = [6.e0_rt/256.e0_rt, &
                                                   -28.e0_rt/256.e0_rt, &
                                                   309.e0_rt/256.e0_rt, &
                                                   -52.e0_rt/256.e0_rt, &
                                                   28.e0_rt/256.e0_rt, &
                                                   -8.e0_rt/256.e0_rt, &
                                                   1.e0_rt/256.e0_rt]

  real(rt), dimension(-3:4), parameter :: ff4b3 = [-4.e0_rt/256.e0_rt, &
                                                   22.e0_rt/256.e0_rt, &
                                                   -52.e0_rt/256.e0_rt, &
                                                   325.e0_rt/256.e0_rt, &
                                                   -56.e0_rt/256.e0_rt, &
                                                   28.e0_rt/256.e0_rt, &
                                                   -8.e0_rt/256.e0_rt, &
                                                   1.e0_rt/256.e0_rt]

end module filter_module
