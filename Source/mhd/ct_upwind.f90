module ct_upwind

 use amrex_fort_module, only : rt => amrex_real
 use hlld_solver, only : hlld
 use meth_params_module

 implicit none

 private primtocons
 public corner_transport

 contains

subroutine corner_transport(lo, hi, &
                            q, q_lo, q_hi, &
                            qm, qm_lo, qm_hi, &
                            qp, qp_lo, qp_hi, &
                            flxx, flxx_lo , flxx_hi, &
                            flxy, flxy_lo , flxy_hi, &
                            flxz, flxz_lo , flxz_hi, &
                            Ex, ex_lo, ex_hi, &
                            Ey, ey_lo, ey_hi, &
                            Ez, ez_lo, ez_hi, &
                            dx, dt) bind(C, name="corner_transport")

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module, only : NQ
  use electric_field
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qm_lo(3), qm_hi(3)
  integer, intent(in) :: qp_lo(3), qp_hi(3)

  integer, intent(in) :: ex_lo(3), ex_hi(3)
  integer, intent(in) :: ey_lo(3), ey_hi(3)
  integer, intent(in) :: ez_lo(3), ez_hi(3)

  integer, intent(in) :: flxx_lo(3), flxx_hi(3)
  integer, intent(in) :: flxy_lo(3), flxy_hi(3)
  integer, intent(in) :: flxz_lo(3), flxz_hi(3)

  real(rt), intent(in) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
  real(rt), intent(in) :: qm(qm_lo(1):qm_hi(1), qm_lo(2):qm_hi(2), qm_lo(3):qm_hi(3), NQ, 3)
  real(rt), intent(in) :: qp(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), NQ, 3)

  ! fluxes should be NVAR+3
  real(rt), intent(out) :: flxx(flxx_lo(1):flxx_hi(1),flxx_lo(2):flxx_hi(2),flxx_lo(3):flxx_hi(3), NVAR+3)   !Half Step Fluxes
  real(rt), intent(out) :: flxy(flxy_lo(1):flxy_hi(1),flxy_lo(2):flxy_hi(2),flxy_lo(3):flxy_hi(3), NVAR+3)   !Half Step Fluxes
  real(rt), intent(out) :: flxz(flxz_lo(1):flxz_hi(1),flxz_lo(2):flxz_hi(2),flxz_lo(3):flxz_hi(3), NVAR+3)   !Half Step Fluxes

  real(rt), intent(out)  :: Ex(ex_lo(1):ex_hi(1), ex_lo(2):ex_hi(2), ex_lo(3):ex_hi(3))
  real(rt), intent(out)  :: Ey(ey_lo(1):ey_hi(1), ey_lo(2):ey_hi(2), ey_lo(3):ey_hi(3))
  real(rt), intent(out)  :: Ez(ez_lo(1):ez_hi(1), ez_lo(2):ez_hi(2), ez_lo(3):ez_hi(3))

  ! these are conserved + magnetic field (cell centered)

        real(rt)  :: um(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3) ,NVAR+3,3) !PtoC Vars
        real(rt)  :: up(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3) ,NVAR+3,3)
        real(rt)  :: cons_temp_M(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3) ,NVAR+3,3,2) !2D Temporary Conservative Vars
        real(rt)  :: cons_temp_P(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3) ,NVAR+3,3,2) !2D Temporary Conservative Vars
        real(rt)  :: cons_half_M(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3) ,NVAR+3,3) !Flux Corrected Conservative Vars
        real(rt)  :: cons_half_P(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3) ,NVAR+3,3)

        real(rt)  :: q_temp_M(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3) ,NQ,3,2) !2D Temporary Primitive Vars
        real(rt)  :: q_temp_P(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3) ,NQ,3,2) !2D Temporary Primitive Vars
        real(rt)  :: q_half_M(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3) ,NQ,3) !Flux Corrected Primitive Vars
        real(rt)  :: q_half_P(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3) ,NQ,3) !Flux Corrected Primitive Vars


        real(rt)  :: flxx1D(flxx_lo(1):flxx_hi(1),flxx_lo(2):flxx_hi(2),flxx_lo(3):flxx_hi(3) ,NVAR+3)
        real(rt)  :: flxy1D(flxy_lo(1):flxy_hi(1),flxy_lo(2):flxy_hi(2),flxy_lo(3):flxy_hi(3), NVAR+3)
        real(rt)  :: flxz1D(flxz_lo(1):flxz_hi(1),flxz_lo(2):flxz_hi(2),flxz_lo(3):flxz_hi(3), NVAR+3) !Flux1d for all directions

        real(rt)  :: flxx2D(flxx_lo(1):flxx_hi(1),flxx_lo(2):flxx_hi(2),flxx_lo(3):flxx_hi(3) ,NVAR+3, 2) !Flux2d for all directions 2 perpendicular directions
        real(rt)  :: flxy2D(flxy_lo(1):flxy_hi(1),flxy_lo(2):flxy_hi(2),flxy_lo(3):flxy_hi(3), NVAR+3, 2) !Flux2d for all directions 2 perpendicular directions
        real(rt)  :: flxz2D(flxz_lo(1):flxz_hi(1),flxz_lo(2):flxz_hi(2),flxz_lo(3):flxz_hi(3), NVAR+3, 2) !Flux2d for all directions 2 perpendicular directions

        real(rt)  :: q2D(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
        real(rt)  :: dx(3)
        real(rt), intent(in), value :: dt

        integer  :: i, work_lo(3), work_hi(3)
        integer  :: UMAGX, UMAGY, UMAGZ



        UMAGX = NVAR+1
        UMAGY = NVAR+2
        UMAGZ = NVAR+3

        um = 0.d0
        up = 0.d0
        cons_temp_M = 0.d0
        cons_temp_P = 0.d0
        q_temp_M = 0.d0
        q_temp_P = 0.d0
        cons_half_M = 0.d0
        cons_half_P = 0.d0
        q_half_M = 0.d0
        q_half_P = 0.d0
        q2D = 0.d0

        flxx1D = 0.d0
        flxy1D = 0.d0
        flxz1D = 0.d0

        flxx2D = 0.d0
        flxy2D = 0.d0
        flxz2D = 0.d0

        Ex = 0.d0
        Ey = 0.d0
        Ez = 0.d0

        !Calculate Flux 1D, eq.35
        !x-dir
        ![lo(1)-2, lo(2)-3, lo(3)-3] [hi(1)+3, hi(2)+3, hi(3)+3]
        work_lo = (/ lo(1)-2, lo(2)-3, lo(3)-3 /)   
        work_hi = (/ hi(1)+3, hi(2)+3, hi(3)+3 /) 
        call hlld(work_lo, work_hi, qm, qm_lo, qm_hi, qp, qp_lo, qp_hi, &
                  flxx1D(:,:,:,:), flxx_lo, flxx_hi, 1)

        !y-dir
        ![lo(1)-3, lo(2)-2, lo(3)-3] [hi(1)+3, hi(2)+3, hi(3)+3]
        work_lo = (/ lo(1)-3, lo(2)-2, lo(3)-3 /) 
        work_hi = (/ hi(1)+3, hi(2)+3, hi(3)+3 /)
        call hlld(work_lo, work_hi, qm, qm_lo, qm_hi, qp, qp_lo, qp_hi, &
                  flxy1D(:,:,:,:), flxy_lo, flxy_hi, 2)

        !z-dir
        ![lo(1)-3, lo(2)-3, lo(3)-2] [hi(1)+3, hi(2)+3, hi(3)+3]
        work_lo = (/ lo(1)-3, lo(2)-3, lo(3)-2 /)
        work_hi = (/ hi(1)+3, hi(2)+3, hi(3)+3 /)
        call hlld(work_lo, work_hi, qm, qm_lo, qm_hi, qp, qp_lo, qp_hi, &
                  flxz1D(:,:,:,:), flxz_lo, flxz_hi, 3)

        !Prim to Cons
        do i = 1,3
           call PrimToCons(qm(:,:,:,:,i), qm_lo, qm_hi, &
                           um(:,:,:,:,i), q_lo, q_hi)
           call PrimToCons(qp(:,:,:,:,i), qp_lo, qp_hi, &
                           up(:,:,:,:,i), q_lo, q_hi)
        enddo


        !Use "1D" fluxes To interpolate Temporary Edge Centered Electric Fields, eq.36

        ![lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+2, hi(2)+3, hi(3)+3]
        work_lo = (/ lo(1)-2, lo(2)-2, lo(3)-2 /)
        work_hi = (/ hi(1)+2, hi(2)+3, hi(3)+3 /)
        call electric_edge_x(work_lo, work_hi, &
                             q, q_lo, q_hi, &
                             Ex, ex_lo, ex_hi, &
                             flxy1D, flxy_lo, flxy_hi, &
                             flxz1D, flxz_lo, flxz_hi)

        ![lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+3, hi(2)+2, hi(3)+3]
        work_lo = (/ lo(1)-2, lo(2)-2, lo(3)-2 /) 
        work_hi = (/ hi(1)+3, hi(2)+2, hi(3)+3 /)
        call electric_edge_y(work_lo, work_hi, &
                             q, q_lo, q_hi, &
                             Ey, ey_lo, ey_hi, &
                             flxx1D, flxx_lo, flxx_hi, &
                             flxz1D, flxz_lo, flxz_hi)

        ![lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+3, hi(2)+3, hi(3)+2]
        work_lo = (/ lo(1)-2, lo(2)-2, lo(3)-2 /)
        work_hi = (/ hi(1)+3, hi(2)+3, hi(3)+2 /)
        call electric_edge_z(work_lo, work_hi, &
                             q, q_lo, q_hi, &
                             Ez, ez_lo, ez_hi, &
                             flxx1D, flxx_lo, flxx_hi, &
                             flxy1D, flxy_lo, flxy_hi)



       !Corner Couple, eq. 37, 38 and 39 Correct Conservative vars using Transverse Fluxes

        !X direction
        ! affected by Y Flux  
        ![lo(1)-2, lo(2)-2, lo(3)-2] [hi(1)+2, hi(2)+2, hi(2)+2] 
        work_lo = (/ lo(1)-2 , lo(2)-2, lo(3)-2 /)
        work_hi = (/ hi(1)+2 , hi(2)+2, hi(3)+2 /) 
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, q_lo, q_hi, &
                           cons_temp_P, q_lo, q_hi, &
                           um, q_lo, q_hi, &
                           up, q_lo, q_hi, &
                           flxy1D, flxy_lo, flxy_hi, &
                           !d1 = x, d2 = y, dir2 = 1 last component in the cons_temp arrays
                           1 , 2, 1, &
                           dx(1), dt) !qmpxy

        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, q_lo, q_hi, &
                               cons_temp_P, q_lo, q_hi, &
                               um, q_lo, q_hi, &
                               up, q_lo, q_hi, &
                               Ex, ex_lo, ex_hi, &
                               Ez, ez_lo, ez_hi, &
                               !x,y,z,dir2=1, sgn=+, UMAGD1, UMAGD2, UMAGD3
                               1, 2, 3, 1, 1, UMAGX, UMAGY, UMAGZ, dx(1), dt)

        call ConsToPrim(work_lo, work_hi, &
                        q_temp_M(:,:,:,:,1,1), q_lo, q_hi, cons_temp_M(:,:,:,:,1,1), q_lo, q_hi)
        call ConsToPrim(work_lo, work_hi, &
                        q_temp_P(:,:,:,:,1,1), q_lo, q_hi, cons_temp_P(:,:,:,:,1,1), q_lo, q_hi)

        ! affected by Z Flux
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, q_lo, q_hi, &
                           cons_temp_P, q_lo, q_hi, &
                           um, q_lo, q_hi, &
                           up, q_lo, q_hi, &
                           flxz1D, flxz_lo, flxz_hi, &
                           1, 3, 2, &
                           dx(1), dt) !qmpxz

        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, q_lo, q_hi, &
                               cons_temp_P, q_lo, q_hi, &
                               um, q_lo, q_hi, &
                               up, q_lo, q_hi, &
                               Ex, ex_lo, ex_hi, &
                               Ey, ey_lo, ey_hi, &
                               1, 3, 2, 2, -1, UMAGX, UMAGZ, UMAGY, dx(1), dt)

        call ConsToPrim(work_lo, work_hi, &
                        q_temp_M(:,:,:,:,1,2), q_lo, q_hi, cons_temp_M(:,:,:,:,1,2), q_lo, q_hi)
        call ConsToPrim(work_lo, work_hi, &
                        q_temp_P(:,:,:,:,1,2), q_lo, q_hi, cons_temp_P(:,:,:,:,1,2), q_lo, q_hi)

        !Y direction
        ! affected by X Flux
        ![lo(1)-2, lo(2)-2, lo(3)-2] [hi(1)+2, hi(2)+2, hi(3)+2] 
        work_lo = (/ lo(1)-2, lo(2)-2, lo(3)-2 /)
        work_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)                
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, q_lo, q_hi, &
                           cons_temp_P, q_lo, q_hi, &
                           um, q_lo, q_hi, &
                           up, q_lo, q_hi, &
                           flxx1D, flxx_lo, flxx_hi, &
                           2, 1, 1, &
                           dx(2), dt) !qmpyx

        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, q_lo, q_hi, &
                               cons_temp_P, q_lo, q_hi, &
                               um, q_lo, q_hi, &
                               up, q_lo, q_hi, &
                               Ey, ey_lo, ey_hi, &
                               Ez, ez_lo, ez_hi, &
                               2, 1, 3, 1, -1, UMAGY, UMAGX, UMAGZ, dx(2), dt)

        call ConsToPrim(work_lo, work_hi, &
                        q_temp_M(:,:,:,:,2,1), q_lo, q_hi, cons_temp_M(:,:,:,:,2,1), q_lo, q_hi)
        call ConsToPrim(work_lo, work_hi, &
                        q_temp_P(:,:,:,:,2,1), q_lo, q_hi, cons_temp_P(:,:,:,:,2,1), q_lo, q_hi)

        ! affected by Z Flux
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, q_lo, q_hi, &
                           cons_temp_P, q_lo, q_hi, &
                           um, q_lo, q_hi, &
                           up, q_lo, q_hi, &
                           flxz1D, flxz_lo, flxz_hi, &
                           2, 3, 2, &
                           dx(2), dt) !qmpyz

        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, q_lo, q_hi, &
                               cons_temp_P, q_lo, q_hi, &
                               um, q_lo, q_hi, &
                               up, q_lo, q_hi, &
                               Ey, ey_lo, ey_hi, &
                               Ex, ex_lo, ex_hi, &
                               2, 3, 1, 2, 1, UMAGY, UMAGZ, UMAGX, dx(2), dt)

        call ConsToPrim(work_lo, work_hi, &
                        q_temp_M(:,:,:,:,2,2), q_lo, q_hi, cons_temp_M(:,:,:,:,2,2), q_lo, q_hi)
        call ConsToPrim(work_lo, work_hi, &
                        q_temp_P(:,:,:,:,2,2), q_lo, q_hi, cons_temp_P(:,:,:,:,2,2), q_lo, q_hi)

        !Z direction
        ! affected by X Flux 
        ![lo(1)-2, lo(2)-2, lo(3)-2] [hi(1)+2, hi(2)+2, hi(3)+2] 
        work_lo = (/ lo(1)-2, lo(2)-2, lo(3)-2 /)
        work_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2/)         
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, q_lo, q_hi, &
                           cons_temp_P, q_lo, q_hi, &
                           um, q_lo, q_hi, &
                           up, q_lo, q_hi, &
                           flxx1D, flxx_lo, flxx_hi, &
                           3, 1, 1, &
                           dx(3), dt) !qmpzx

        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, q_lo, q_hi, &
                               cons_temp_P, q_lo, q_hi, &
                               um, q_lo, q_hi, &
                               up, q_lo, q_hi, &
                               Ez, ez_lo, ez_hi, &
                               Ey, ey_lo, ey_hi, &
                               3, 1, 2, 1, 1, UMAGZ, UMAGX, UMAGY, dx(3), dt)

        call ConsToPrim(work_lo, work_hi, &
                        q_temp_M(:,:,:,:,3,1), q_lo, q_hi, cons_temp_M(:,:,:,:,3,1), q_lo, q_hi)
        call ConsToPrim(work_lo, work_hi, &
                        q_temp_P(:,:,:,:,3,1), q_lo, q_hi, cons_temp_P(:,:,:,:,3,1), q_lo, q_hi)

        ! affected by Y Flux
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, q_lo, q_hi, &
                           cons_temp_P, q_lo, q_hi, &
                           um, q_lo, q_hi, &
                           up, q_lo, q_hi, &
                           flxy1D, flxy_lo, flxy_hi, &
                           3, 2, 2, &
                           dx(3), dt) !qmpzy

        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, q_lo, q_hi, &
                               cons_temp_P, q_lo, q_hi, &
                               um, q_lo, q_hi, &
                               up, q_lo, q_hi, &
                               Ez, ez_lo, ez_hi, &
                               Ex, ex_lo, ex_hi, &
                               3, 2, 1, 2, -1, UMAGZ, UMAGY, UMAGX, dx(3), dt)

        call ConsToPrim(work_lo, work_hi, &
                        q_temp_M(:,:,:,:,3,2), q_lo, q_hi, cons_temp_M(:,:,:,:,3,2), q_lo, q_hi)
        call ConsToPrim(work_lo, work_hi, &
                        q_temp_P(:,:,:,:,3,2), q_lo, q_hi, cons_temp_P(:,:,:,:,3,2), q_lo, q_hi)

        !Calculate Flux 2D eq. 40
        ![lo(1)-1, lo(2)-2, lo(3)-2][hi(1)+2,hi(2)+2,hi(3)+2]
        work_lo = (/ lo(1)-1, lo(2)-2, lo(3)-2 /)
        work_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)
        !x-dir
        call hlld(work_lo, work_hi, &
                  q_temp_M(:,:,:,:,:,1), q_lo, q_hi, q_temp_P(:,:,:,:,:,1), q_lo, q_hi, &
                  flxx2D(:,:,:,:,1), flxx_lo, flxx_hi, 1) !F^{x|y}
        call hlld(work_lo, work_hi, &
                  q_temp_M(:,:,:,:,:,2), q_lo, q_hi, q_temp_P(:,:,:,:,:,2), q_lo, q_hi, &
                  flxx2D(:,:,:,:,2), flxx_lo, flxx_hi, 1) !F^{x|z}
        !y-dir
        ![lo(1)-2, lo(2)-1, lo(3)-2][hi(1)+2,hi(2)+2,hi(3)+2]
        work_lo = (/ lo(1)-2, lo(2)-1, lo(3)-2 /)
        work_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)
        call hlld(work_lo, work_hi, &
                  q_temp_M(:,:,:,:,:,1), q_lo, q_hi, q_temp_P(:,:,:,:,:,1), q_lo, q_hi, &
                  flxy2D(:,:,:,:,1), flxy_lo, flxy_hi, 2) !F^{y|x}
        call hlld(work_lo, work_hi, &
                  q_temp_M(:,:,:,:,:,2), q_lo, q_hi, q_temp_P(:,:,:,:,:,2), q_lo, q_hi, &
                  flxy2D(:,:,:,:,2), flxy_lo, flxy_hi, 2) !F^{y|z}

        !z-dir
        ![lo(1)-2,lo(2)-2,lo(3)-1][h1(1)+2, h1(2)+2, h1(3)+2]
        work_lo = (/ lo(1)-2, lo(2)-2, lo(3)-1 /)
        work_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)
        call hlld(work_lo, work_hi, &
                  q_temp_M(:,:,:,:,:,1), q_lo, q_hi, q_temp_P(:,:,:,:,:,1), q_lo, q_hi, &
                  flxz2D(:,:,:,:,1), flxz_lo, flxz_hi, 3) !F^{z|x}
        call hlld(work_lo, work_hi, &
                  q_temp_M(:,:,:,:,:,2), q_lo, q_hi, q_temp_P(:,:,:,:,:,2), q_lo, q_hi, &
                  flxz2D(:,:,:,:,2), flxz_lo, flxz_hi, 3) !F^{z|y}


        !Use Averaged 2D fluxes to interpolate temporary Edge Centered Electric Fields, reuse "flx1D"
        ! eq.  42 and 43 ?
        flxx1D(:,:,:,:) = 0.5d0*(flxx2D(:,:,:,:,1) + flxx2D(:,:,:,:,2))
        flxy1D(:,:,:,:) = 0.5d0*(flxy2D(:,:,:,:,1) + flxy2D(:,:,:,:,2))
        flxz1D(:,:,:,:) = 0.5d0*(flxz2D(:,:,:,:,1) + flxz2D(:,:,:,:,2))

        !eq. 41
        ![lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+2, hi(3)+2]  
        work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
        work_hi = (/ hi(1)+1, hi(2)+2, hi(3)+2 /)
        call electric_edge_x(work_lo, work_hi, &
                             q, q_lo, q_hi, &
                             Ex, ex_lo, ex_hi, &
                             flxy1D, flxy_lo, flxy_hi, &
                             flxz1D, flxz_lo, flxz_hi)

        ![lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+2, hi(2)+1, hi(3)+2]
        work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
        work_hi = (/ hi(1)+2, hi(2)+1, hi(3)+2 /)
        call electric_edge_y(work_lo, work_hi, &
                             q, q_lo, q_hi, &
                             Ey, ey_lo, ey_hi, &
                             flxx1D, flxx_lo, flxx_hi, &
                             flxz1D, flxz_lo, flxz_hi)

        ![lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+2, hi(2)+2, hi(3)+1]
        work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
        work_hi = (/ hi(1)+2, hi(2)+2, hi(3)+1 /)
        call electric_edge_z(work_lo, work_hi, &
                             q, q_lo, q_hi, &
                             Ez, ez_lo, ez_hi, &
                             flxx1D, flxx_lo, flxx_hi, &
                             flxy1D, flxy_lo, flxy_hi)

        !Half Step conservative vars eq.44, eq.45, eq.46

        !for x direction
        ![lo(1)-1,lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
        work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
        work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
        call half_step(work_lo, work_hi, &
                       cons_half_M, q_lo, q_hi, &
                       cons_half_P, q_lo, q_hi, &
                       um, q_lo, q_hi, &
                       up, q_lo, q_hi, &
                       flxy2D, flxy_lo, flxy_hi, &
                       flxz2D, flxz_lo, flxz_hi, &
                       !dir = x, d1 =y, d2 =z
                       1, 2, 3, dx(1), dt)

        call half_step_mag(work_lo, work_hi, &
                           cons_half_M, q_lo, q_hi, &
                           cons_half_P, q_lo, q_hi, &
                           um, q_lo, q_hi, &
                           up, q_lo, q_hi, &
                           Ex, ex_lo, ex_hi, &
                           Ey, ey_lo, ey_hi, &
                           Ez, ez_lo, ez_hi, &
                           !d=x, d1=y, d2=z, UMAGD UMAGD1, UMAGD2, sgn,
                           1, 2, 3, UMAGX, UMAGY, UMAGZ, -1, &
                           dx(1), dt)

        call ConsToPrim(work_lo, work_hi, &
                        q_half_M(:,:,:,:,1), q_lo, q_hi, cons_half_M(:,:,:,:,1), q_lo, q_hi)
        call ConsToPrim(work_lo, work_hi, &
                        q_half_P(:,:,:,:,1), q_lo, q_hi, cons_half_P(:,:,:,:,1), q_lo, q_hi)


        !for y direction
        ![lo(1)-1,lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
        work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
        work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
        call half_step(work_lo, work_hi, &
                       cons_half_M, q_lo, q_hi, &
                       cons_half_P, q_lo, q_hi, &
                       um, q_lo, q_hi, &
                       up, q_lo, q_hi, &
                       flxx2D, flxx_lo, flxx_hi, &
                       flxz2D, flxz_lo, flxz_hi, &
                       2, 1, 3, dx(2), dt)

        call half_step_mag(work_lo, work_hi, &
                           cons_half_M, q_lo, q_hi, &
                           cons_half_P, q_lo, q_hi, &
                           um, q_lo, q_hi, &
                           up, q_lo, q_hi, &
                           Ey, ey_lo, ey_hi, &
                           Ex, ex_lo, ex_hi, &
                           Ez, ez_lo, ez_hi, &
                           !d, d1, d2, UMAGD UMAGD1, UMAGD2, sgn,
                           2, 1, 3, UMAGY, UMAGX, UMAGZ, 1, &
                           dx(2), dt)

        call ConsToPrim(work_lo, work_hi, &
                        q_half_M(:,:,:,:,2), q_lo, q_hi, cons_half_M(:,:,:,:,2), q_lo, q_hi)
        call ConsToPrim(work_lo, work_hi, &
                        q_half_P(:,:,:,:,2), q_lo, q_hi, cons_half_P(:,:,:,:,2), q_lo, q_hi)

        !for z direction
        ![lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
        work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
        work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
        call half_step(work_lo, work_hi, &
                       cons_half_M, q_lo, q_hi, &
                       cons_half_P, q_lo, q_hi, &
                       um, q_lo, q_hi, &
                       up, q_lo, q_hi, &
                       flxx2D, flxx_lo, flxx_hi, &
                       flxy2D, flxy_lo, flxy_hi, &
                       3, 1, 2, dx(3), dt)

        call half_step_mag(work_lo, work_hi, &
                           cons_half_M, q_lo, q_hi, &
                           cons_half_P, q_lo, q_hi, &
                           um, q_lo, q_hi, &
                           up, q_lo, q_hi, &
                           Ez, ez_lo, ez_hi, &
                           Ex, ex_lo, ex_hi, &
                           Ey, ey_lo, ey_hi, &
                           !d, d1, d2, UMAGD UMAGD1, UMAGD2, sgn,
                           3, 1, 2, UMAGZ, UMAGX, UMAGY, -1, &
                           dx(3), dt)

        call ConsToPrim(work_lo, work_hi, &
                        q_half_M(:,:,:,:,3), q_lo, q_hi, cons_half_M(:,:,:,:,3), q_lo, q_hi)
        call ConsToPrim(work_lo, work_hi, &
                        q_half_P(:,:,:,:,3), q_lo, q_hi, cons_half_P(:,:,:,:,3), q_lo, q_hi)


        !Final Fluxes eq.47

        !x-dir
        ![lo(1), lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
        work_lo = (/ lo(1), lo(2)-1, lo(3)-1 /)
        work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
        call hlld(work_lo, work_hi, &
                  q_half_M, q_lo, q_hi, q_half_P, q_lo, q_hi, &
                  flxx, flxx_lo, flxx_hi, 1)

        !y-dir
        ![lo(1)-1, lo(2), lo(3)-1][hi(1)+1,hi(2)+1, hi(3)+1]
        work_lo = (/ lo(1)-1, lo(2), lo(3)-1 /)
        work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
        call hlld(work_lo, work_hi, &
                  q_half_M, q_lo, q_hi, q_half_P, q_lo, q_hi, &
                  flxy, flxy_lo, flxy_hi, 2)

        !z-dir
        ![lo(1)-1,lo(2)-1,lo(3)][hi(1)+1, hi(2)+1, hi(3)+1]
        work_lo = (/ lo(1)-1, lo(2)-1, lo(3) /)
        work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
        call hlld(work_lo, work_hi, &
                  q_half_M, q_lo, q_hi, q_half_P, q_lo, q_hi, &
                  flxz, flxz_lo, flxz_hi, 3)

        !Primitive update eq. 48
        ![lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
        work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
        work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
        call prim_half(work_lo, work_hi, &
                       q2D, q_lo, q_hi, &
                       q, q_lo, q_hi, &
                       flxx1D, flxx_lo, flxx_hi, &
                       flxy1D, flxy_lo, flxy_hi, &
                       flxz1D, flxz_lo, flxz_hi, &
                       dx(1), dx(2), dx(3), dt)

        !Final Electric Field Update eq.48
        ![lo(1), lo(2), lo(3)][hi(1), hi(2)+1, hi(3)+1]
        work_lo = (/ lo(1), lo(2), lo(3) /)
        work_hi = (/ hi(1), hi(2)+1, hi(3)+1 /)
        call electric_edge_x(work_lo, work_hi, &
                             q2D, q_lo, q_hi, &
                             Ex, ex_lo, ex_hi, &
                             flxy, flxy_lo, flxy_hi, &
                             flxz, flxz_lo, flxz_hi)

        ![lo(1), lo(2), lo(3)][hi(1)+1, hi(2), hi(3)+1]
        work_lo = (/ lo(1), lo(2), lo(3) /)
        work_hi = (/ hi(1)+1, hi(2), hi(3)+1 /)
        call electric_edge_y(work_lo, work_hi, &
                             q2D, q_lo, q_hi, &
                             Ey, ey_lo, ey_hi, &
                             flxx, flxx_lo, flxx_hi, &
                             flxz, flxz_lo, flxz_hi)

        ![lo(1), lo(2), lo(3)][hi(1)+1, hi(2)+1 ,hi(3)]
        work_lo = (/ lo(1), lo(2), lo(3) /)
        work_hi = (/ hi(1)+1, hi(2)+1, hi(3) /)
        call electric_edge_z(work_lo, work_hi, &
                             q2D, q_lo, q_hi, &
                             Ez, ez_lo, ez_hi, &
                             flxx, flxx_lo, flxx_hi, &
                             flxy, flxy_lo, flxy_hi)


end subroutine corner_transport

!================================================= Calculate the Conservative Variables ===============================================

subroutine PrimToCons(q, q_lo, q_hi, u, u_lo, u_hi)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module
 use eos_type_module, only : eos_t, eos_input_rp
 use eos_module, only: eos
 use network, only: nspec

 implicit none

 integer, intent(in) :: q_lo(3), q_hi(3)
 integer, intent(in) :: u_lo(3), u_hi(3)

 real(rt), intent(in) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
 real(rt), intent(out) :: u(u_lo(1):u_hi(1), u_lo(2):u_hi(2), u_lo(3):u_hi(3), NVAR+3)

 integer :: i ,j ,k

 type(eos_t) :: eos_state

 do k = q_lo(3)+1, q_hi(3)-1
    do j = q_lo(2)+1, q_hi(2)-1
       do i = q_lo(1)+1, q_hi(1)-1

          u(i,j,k,URHO)  = q(i,j,k,QRHO)
          u(i,j,k,UMX)    = q(i,j,k,QRHO)*q(i,j,k,QU)
          u(i,j,k,UMY)    = q(i,j,k,QRHO)*q(i,j,k,QV)
          u(i,j,k,UMZ)    = q(i,j,k,QRHO)*q(i,j,k,QW)

          eos_state % rho = q(i, j, k, QRHO)
          eos_state % p   = q(i, j, k, QPRES)
          eos_state % T   = 100.d0 !dummy initial g.
          eos_state % xn  = q(i, j, k, QFS:QFS+nspec-1)

          call eos(eos_input_rp, eos_state)

          u(i,j,k,UEDEN) = eos_state % rho * eos_state % e &
               + 0.5d0*q(i,j,k,QRHO)*dot_product(q(i,j,k,QU:QW),q(i,j,k,QU:QW)) &
               + 0.5d0*(dot_product(q(i,j,k,QMAGX:QMAGZ),q(i,j,k,QMAGX:QMAGZ)))

          u(i,j,k,UEINT) = eos_state % rho * eos_state % e

          u(i,j,k,NVAR+1:NVAR+3) = q(i,j,k,QMAGX:QMAGZ)

          ! species
          u(i,j,k,UFS:UFS-1+nspec) = q(i,j,k,QRHO) * q(i,j,k,QFS:QFS-1+nspec)

       enddo
    enddo
 enddo
end subroutine PrimToCons

!================================================= Calculate the Primitve Variables ===============================================

subroutine ConsToPrim(w_lo, w_hi, q, q_lo, q_hi, u, u_lo, u_hi)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module
 use eos_module, only : eos
 use eos_type_module, only : eos_t, eos_input_re
 use network, only : nspec

 implicit none

 integer, intent(in) :: w_lo(3), w_hi(3)
 integer, intent(in) :: q_lo(3), q_hi(3)
 integer, intent(in) :: u_lo(3), u_hi(3)
 real(rt), intent(in) :: u(q_lo(1):q_hi(1), q_lo(2):q_hi(2),q_lo(3):q_hi(3),NVAR+3)
 real(rt), intent(out) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3),NQ)

 integer :: i ,j ,k
 integer :: UMAGX, UMAGZ

 type (eos_t) :: eos_state

 UMAGX = NVAR+1
 UMAGZ = NVAR+3

 !q = u
 do k = w_lo(3), w_hi(3)
    do j = w_lo(2), w_hi(2)
       do i = w_lo(1), w_hi(1)
          q(i,j,k,QRHO)  = u(i,j,k,URHO)
          q(i,j,k,QU)    = u(i,j,k,UMX)/q(i,j,k,QRHO)
          q(i,j,k,QV)    = u(i,j,k,UMY)/q(i,j,k,QRHO)
          q(i,j,k,QW)    = u(i,j,k,UMZ)/q(i,j,k,QRHO)
          q(i,j,k,QREINT) = u(i,j,k,UEDEN) - 0.5d0*q(i,j,k,QRHO)*dot_product(q(i,j,k,QU:QW),q(i,j,k,QU:QW)) &
               - 0.5d0*dot_product(u(i,j,k,UMAGX:UMAGZ), u(i,j,k,UMAGX:UMAGZ)) !NVAR+1->UMAGX

          ! species
          q(i,j,k,QFS:QFS-1+nspec) = u(i,j,k,UFS:UFS-1+nspec)/u(i,j,k,URHO)

          eos_state % rho = q(i, j, k, QRHO)
          eos_state % e   = q(i, j, k, QREINT)/eos_state % rho
          eos_state % xn  = q(i, j, k, QFS:QFS+nspec-1)
          eos_state % T   = 100.d0  ! initial guess

          call eos(eos_input_re, eos_state)

          q(i,j,k,QTEMP) = eos_state % T
          q(i,j,k,QPRES) = eos_state % p
          q(i,j,k,QMAGX:QMAGZ) = u(i,j,k,UMAGX:UMAGZ)
       enddo
    enddo
 enddo
end subroutine ConsToPrim

!======================================= Update the Temporary Conservative Variables with Transverse 1D Fluxes ========================
subroutine corner_couple(w_lo, w_hi, &
                         uL, ul_lo, ul_hi, &
                         uR, ur_lo, ur_hi, &
                         um, um_lo, um_hi, &
                         up, up_lo, up_hi, &
                         flxd2, flxd2_lo, flxd2_hi, &
                         d1, d2, dir2, &
                         dx, dt)

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module
  use network, only : nspec

  implicit none

  integer, intent(in) :: w_lo(3), w_hi(3)
  integer, intent(in) :: ul_lo(3), ul_hi(3)
  integer, intent(in) :: ur_lo(3), ur_hi(3)
  integer, intent(in) :: um_lo(3), um_hi(3)
  integer, intent(in) :: up_lo(3), up_hi(3)
  integer, intent(in) :: flxd2_lo(3), flxd2_hi(3)
  integer, intent(in) :: d1, d2, dir2
  real(rt), intent(in) :: dx, dt

  real(rt), intent(in) :: um(um_lo(1):um_hi(1), um_lo(2):um_hi(2), um_lo(3):um_hi(3), NVAR+3,3)
  real(rt), intent(in) :: up(up_lo(1):up_hi(1), up_lo(2):up_hi(2), up_lo(3):up_hi(3), NVAR+3,3)

  real(rt), intent(out) :: flxd2(flxd2_lo(1):flxd2_hi(1),flxd2_lo(2):flxd2_hi(2),flxd2_lo(3):flxd2_hi(3),NVAR+3)

  real(rt), intent(out) :: uL(ul_lo(1):ul_hi(1),ul_lo(2):ul_hi(2),ul_lo(3):ul_hi(3),NVAR+3,3,2)
  real(rt), intent(out) :: uR(ur_lo(1):ur_hi(1),ur_lo(2):ur_hi(2),ur_lo(3):ur_hi(3),NVAR+3,3,2)

  real(rt) :: u, v, w
  integer  :: i ,j ,k
  integer  :: d(3) !for the addition of +1 to either i,j,k depending on d2

  uL(:,:,:,:,d1,dir2) = um(:,:,:,:,d1)
  uR(:,:,:,:,d1,dir2) = up(:,:,:,:,d1)


  d = 0

  !the first term of the flxd2 substraction is shifted by 1 on the direction d2
  d(d2) = 1

  do k = w_lo(3), w_hi(3)
     do j = w_lo(2), w_hi(2)
        do i = w_lo(1), w_hi(1)

           ! eq. 37 from Miniati paper, for both + and -

           !Left Corrected States
           uL(i,j,k,URHO,d1,dir2) = um(i,j,k,URHO,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),URHO) - flxd2(i,j,k,URHO))
           uL(i,j,k,UMX,d1,dir2) = um(i,j,k,UMX,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),UMX) - flxd2(i,j,k,UMX))
           uL(i,j,k,UMY,d1,dir2) = um(i,j,k,UMY,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),UMY) - flxd2(i,j,k,UMY))
           uL(i,j,k,UMZ,d1,dir2) = um(i,j,k,UMZ,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),UMZ) - flxd2(i,j,k,UMZ))
           uL(i,j,k,UEDEN,d1,dir2) = um(i,j,k,UEDEN,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),UEDEN) - flxd2(i,j,k,UEDEN))
           uL(i,j,k,UFS:UFS+nspec-1,d1,dir2) = um(i,j,k,UFS:UFS+nspec-1,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),UFS:UFS+nspec-1) &
                                             - flxd2(i,j,k,UFS:UFS+nspec-1))


           u = uL(i,j,k,UMX,d1,dir2)/uL(i,j,k,URHO,d1,dir2)
           v = uL(i,j,k,UMY,d1,dir2)/uL(i,j,k,URHO,d1,dir2)
           w = uL(i,j,k,UMZ,d1,dir2)/uL(i,j,k,URHO,d1,dir2)
           uL(i,j,k,UEINT,d1,dir2) = uL(i,j,k,UEDEN,d1,dir2) - 0.5d0*uL(i,j,k,URHO,d1,dir2)*(u**2 + v**2 + w**2)


           !Right Corrected States
           uR(i,j,k,URHO,d1,dir2) = up(i,j,k,URHO,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),URHO) - flxd2(i,j,k,URHO))
           uR(i,j,k,UMX,d1,dir2) = up(i,j,k,UMX,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),UMX) - flxd2(i,j,k,UMX))
           uR(i,j,k,UMY,d1,dir2) = up(i,j,k,UMY,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),UMY) - flxd2(i,j,k,UMY))
           uR(i,j,k,UMZ,d1,dir2) = up(i,j,k,UMZ,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),UMZ) - flxd2(i,j,k,UMZ))
           uR(i,j,k,UEDEN,d1,dir2) = up(i,j,k,UEDEN,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),UEDEN) - flxd2(i,j,k,UEDEN))
           uR(i,j,k,UFS:UFS+nspec-1,d1,dir2) = up(i,j,k,UFS:UFS+nspec-1,d1) - dt/(3.d0*dx)*(flxd2(i+d(1),j+d(2),k+d(3),UFS:UFS+nspec-1) &
                                             - flxd2(i,j,k,UFS:UFS+nspec-1))


           u = uR(i,j,k,UMX,d1,dir2)/uR(i,j,k,URHO,d1,dir2)
           v = uR(i,j,k,UMY,d1,dir2)/uR(i,j,k,URHO,d1,dir2)
           w = uR(i,j,k,UMZ,d1,dir2)/uR(i,j,k,URHO,d1,dir2)
           uR(i,j,k,UEINT,d1,dir2) = uR(i,j,k,UEDEN,d1,dir2) - 0.5d0*uR(i,j,k,URHO,d1,dir2)*(u**2 + v**2 + w**2)

        enddo
     enddo
  enddo
end subroutine corner_couple

!================================== Use 1D Electric Fields to Transverse correct the Temporary Magnetic Fields ===========================
subroutine corner_couple_mag(w_lo, w_hi, &
                             uL, ul_lo, ul_hi, &
                             uR, ur_lo, ur_hi, &
                             um, um_lo, um_hi, &
                             up, up_lo, up_hi, &
                             Ed1, ed1_lo, ed1_hi, &
                             Ed3, ed3_lo, ed3_hi, &
                             d1, d2, d3, dir2, sgn, UMAGD1, UMAGD2, UMAGD3, &
                             dx, dt)
  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module, only : NVAR, UEINT

  !Correction using Faraday's Law
  implicit none

  integer, intent(in) :: w_lo(3), w_hi(3)
  integer, intent(in) :: ul_lo(3), ul_hi(3)
  integer, intent(in) :: ur_lo(3), ur_hi(3)
  integer, intent(in) :: um_lo(3), um_hi(3)
  integer, intent(in) :: up_lo(3), up_hi(3)
  integer, intent(in) :: ed1_lo(3), ed1_hi(3)
  integer, intent(in) :: ed3_lo(3), ed3_hi(3)
  integer, intent(in) :: d1, d2, d3, dir2, sgn, UMAGD1, UMAGD2, UMAGD3   !UMAGD1 corresponds to d1, and UMAGD2 to dir2

  real(rt), intent(inout) :: uL(ul_lo(1):ul_hi(1), ul_lo(2):ul_hi(2), ul_lo(3):ul_hi(3), NVAR+3,3,2)
  real(rt), intent(inout) :: uR(ur_lo(1):ur_hi(1), ur_lo(2):ur_hi(2), ur_lo(3):ur_hi(3), NVAR+3,3,2)
  real(rt), intent(in)    :: um(um_lo(1):um_hi(1), um_lo(2):um_hi(2), um_lo(3):um_hi(3), NVAR+3,3)
  real(rt), intent(in)    :: up(up_lo(1):up_hi(1), up_lo(2):up_hi(2), up_lo(3):up_hi(3), NVAR+3,3)

  real(rt), intent(in) :: Ed1(ed1_lo(1):ed1_hi(1),ed1_lo(2):ed1_hi(2),ed1_lo(3):ed1_hi(3))
  real(rt), intent(in) :: Ed3(ed3_lo(1):ed3_hi(1),ed3_lo(2):ed3_hi(2),ed3_lo(3):ed3_hi(3))

  real(rt) :: dx, dt
  integer :: i ,j ,k
  integer :: d(3), a1(3), a2(3), a3(3), d_2(3) !for the additions of +1 to i,j,k
  integer :: UMAGX, UMAGY, UMAGZ

  UMAGX = NVAR+1
  UMAGY = NVAR+2
  UMAGZ = NVAR+3

  d   = 0
  a1  = 0
  a2  = 0
  a3  = 0
  d_2 = 0


  d(d2)  = 1   !for example if d2 = y, j+1 on Ed3
  a1(d2) = 1   !j+1 on first and third term of addition Ed1
  a3(d2) = 1

  a1(d3) = 1
  a2(d3) = 1   !the second term of addition Ed1 increments by 1 the i,j,k


  do k = w_lo(3), w_hi(3)
     do j = w_lo(2), w_hi(2)
        do i = w_lo(1), w_hi(1)
           !Left State

           ! d1 -direction
           ! affected by dir2 flux

           ! eq. 38 and 39 of Miniati for -

           uL(i,j,k,UMAGD1,d1,dir2) = um(i,j,k,UMAGD1,d1) + sgn*dt/(3.d0*dx)*(Ed3(i+d(1),j+d(2),k+d(3)) - Ed3(i,j,k))
           uL(i,j,k,UMAGD3,d1,dir2) = um(i,j,k,UMAGD3,d1) + sgn*dt/(6.d0*dx)* &
                ((Ed1(i+a1(1),j+a1(2),k+a1(3)) - Ed1(i+a2(1),j+a2(2),k+a2(3))) + &
                 (Ed1(i+a3(1),j+a3(2),k+a3(3)) - Ed1(i,j,k)))

           uL(i,j,k,UMAGD2,d1,dir2) = um(i,j,k,UMAGD2,d1)


           uL(i,j,k,UEINT,d1,dir2) = uL(i,j,k,UEINT,d1,dir2) -0.5d0*dot_product(uL(i,j,k,UMAGX:UMAGZ,d1,dir2),uL(i,j,k,UMAGX:UMAGZ,d1,dir2))

        enddo
     enddo
  enddo


  d(d1)  = 1
  d_2(d1) = 1


  do k = w_lo(3), w_hi(3)
     do j = w_lo(2), w_hi(2)
        do i = w_lo(1), w_hi(1)
           !Right State

           !d1 -direction
           !-> Affected by dir2 flux

           ! eq. 38 and 39 of Miniati for +

           uR(i,j,k,UMAGD1,d1,dir2) = up(i,j,k,UMAGD1,d1) + sgn*dt/(3.d0*dx)*(Ed3(i+d(1),j+d(2),k+d(3)) - Ed3(i+d_2(1),j+d_2(2),k+d_2(3)))
           uR(i,j,k,UMAGD3,d1,dir2) = up(i,j,k,UMAGD3,d1) + sgn*dt/(6.d0*dx)*&
                ((Ed1(i+a1(1),j+a1(2),k+a1(3)) - Ed1(i+a2(1),j+a2(2),k+a2(3))) + &
                 (Ed1(i+a3(1),j+a3(2),k+a3(3)) - Ed1(i,j,k)))


           uR(i,j,k,UMAGD2,d1,dir2) = up(i,j,k,UMAGD2,d1)


           uR(i,j,k,UEINT,d1,dir2) = uR(i,j,k,UEINT,d1,dir2) -0.5d0*dot_product(uR(i,j,k,UMAGX:UMAGZ,d1,dir2),uR(i,j,k,UMAGX:UMAGZ,d1,dir2))

        enddo
     enddo
  enddo

end subroutine corner_couple_mag

!====================================================== Final Conservative Corrections================================================================
subroutine half_step(w_lo, w_hi, &
                     uL, ul_lo, ul_hi, &
                     uR, ur_lo, ur_hi, &
                     um, um_lo, um_hi, &
                     up, up_lo, up_hi, &
                     flxd1, flxd1_lo, flxd1_hi, &
                     flxd2, flxd2_lo, flxd2_hi, &
                     dir, d1, d2, dx, dt)

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UMZ, URHO, UEINT, UFS
  use network, only: nspec

  implicit none

  integer, intent(in) :: w_lo(3), w_hi(3)
  integer, intent(in) :: ul_lo(3), ul_hi(3)
  integer, intent(in) :: ur_lo(3), ur_hi(3)
  integer, intent(in) :: um_lo(3), um_hi(3)
  integer, intent(in) :: up_lo(3), up_hi(3)

  integer, intent(in) :: flxd1_lo(3), flxd1_hi(3)
  integer, intent(in) :: flxd2_lo(3), flxd2_hi(3)

  real(rt), intent(inout) :: uL(ul_lo(1):ul_hi(1), ul_lo(2):ul_hi(2), ul_lo(3):ul_hi(3), NVAR+3,3)
  real(rt), intent(inout) :: uR(ur_lo(1):ur_hi(1), ur_lo(2):ur_hi(2), ur_lo(3):ur_hi(3), NVAR+3,3)
  real(rt), intent(in)    :: um(um_lo(1):um_hi(1), um_lo(2):um_hi(2), um_lo(3):um_hi(3), NVAR+3,3)
  real(rt), intent(in)    :: up(up_lo(1):up_hi(1), up_lo(2):up_hi(2), up_lo(3):up_hi(3), NVAR+3,3)

  real(rt), intent(in)  :: dx, dt !dx will be dx, dy or dz
  integer, intent(in)   :: dir, d1, d2 ! following notation of eq. 44

  real(rt), intent(in)  :: flxd1(flxd1_lo(1):flxd1_hi(1),flxd1_lo(2):flxd1_hi(2),flxd1_lo(3):flxd1_hi(3),NVAR+3,2)
  real(rt), intent(in)  :: flxd2(flxd2_lo(1):flxd2_hi(1),flxd2_lo(2):flxd2_hi(2),flxd2_lo(3):flxd2_hi(3),NVAR+3,2)

  real(rt) :: u, v, w
  integer  :: i ,j ,k, n
  integer  :: d(3), d_2(3) !for the shift in i,j,k
  integer  :: flxd1c, flxd2c !last component of flxd1 and flxd2

  uL(:,:,:,:,dir) = um(:,:,:,:,dir)
  uR(:,:,:,:,dir) = up(:,:,:,:,dir)

  d = 0
  d_2 = 0

  d(d1) = 1   ! add +1 to the d1 direction in the first flxd1 term of the substraction
  d_2(d2) = 1 ! add +1 to the d2 direction in the first flxd2 term of the substraction

  if (dir .eq. 1) then  !d=x, d1=y and d2=z
     flxd1c = 2  !Fy, z
     flxd2c = 2  !Fz, y
  endif

  if (dir .eq. 2) then !d=y, d1=x  and d2=z
     flxd1c = 2 !Fx, z
     flxd2c = 1 !Fz, x
  endif

  if (dir .eq. 3)  then !d=z, d1=x and d2=y
     flxd1c = 1 !Fy, x
     flxd2c = 1 !Fx, y
  endif

  do k = w_lo(3), w_hi(3)
     do j = w_lo(2), w_hi(2)
        do i = w_lo(1), w_hi(1)

           ! eq. 44 of Miniati paper for + and -

!left state
           uL(i,j,k,URHO,dir) = um(i,j,k,URHO,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),URHO,flxd1c) - flxd1(i,j,k,URHO,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),URHO,flxd2c) - flxd2(i,j,k,URHO,flxd2c))
           uL(i,j,k,UMX,dir) = um(i,j,k,UMX,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),UMX,flxd1c) - flxd1(i,j,k,UMX,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),UMX,flxd2c) - flxd2(i,j,k,UMX,flxd2c))
           uL(i,j,k,UMY,dir) = um(i,j,k,UMY,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),UMY,flxd1c) - flxd1(i,j,k,UMY,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),UMY,flxd2c) - flxd2(i,j,k,UMY,flxd2c))
           uL(i,j,k,UMZ,dir) = um(i,j,k,UMZ,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),UMZ,flxd1c) - flxd1(i,j,k,UMZ,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),UMZ,flxd2c) - flxd2(i,j,k,UMZ,flxd2c))
           uL(i,j,k,UEDEN,dir) = um(i,j,k,UEDEN,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),UEDEN,flxd1c) - flxd1(i,j,k,UEDEN,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),UEDEN,flxd2c) - flxd2(i,j,k,UEDEN,flxd2c))
           uL(i,j,k,UFS:UFS+nspec-1,dir) = um(i,j,k,UFS:UFS+nspec-1,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),UFS:UFS+nspec-1,flxd1c) &
                - flxd1(i,j,k,UFS:UFS+nspec-1,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),UFS:UFS+nspec-1,flxd2c) &
                - flxd2(i,j,k,UFS:UFS+nspec-1,flxd2c))


           u = uL(i,j,k,UMX,dir)/uL(i,j,k,URHO,dir)
           v = uL(i,j,k,UMY,dir)/uL(i,j,k,URHO,dir)
           w = uL(i,j,k,UMZ,dir)/uL(i,j,k,URHO,dir)

           uL(i,j,k,UEINT,dir) = uL(i,j,k,UEDEN,dir) - 0.5d0*uL(i,j,k,URHO,dir)*(u**2 + v**2 + w**2)

           !right state
           uR(i,j,k,URHO,dir) = up(i,j,k,URHO,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),URHO,flxd1c) - flxd1(i,j,k,URHO,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),URHO,flxd2c) - flxd2(i,j,k,URHO,flxd2c))
           uR(i,j,k,UMX,dir) = up(i,j,k,UMX,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),UMX,flxd1c) - flxd1(i,j,k,UMX,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),UMX,flxd2c) - flxd2(i,j,k,UMX,flxd2c))
           uR(i,j,k,UMY,dir) = up(i,j,k,UMY,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),UMY,flxd1c) - flxd1(i,j,k,UMY,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),UMY,flxd2c) - flxd2(i,j,k,UMY,flxd2c))
           uR(i,j,k,UMZ,dir) = up(i,j,k,UMZ,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),UMZ,flxd1c) - flxd1(i,j,k,UMZ,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),UMZ,flxd2c) - flxd2(i,j,k,UMZ,flxd2c))
           uR(i,j,k,UEDEN,dir) = up(i,j,k,UEDEN,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),UEDEN,flxd1c) - flxd1(i,j,k,UEDEN,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),UEDEN,flxd2c) - flxd2(i,j,k,UEDEN,flxd2c))
           uR(i,j,k,UFS:UFS+nspec-1,dir) = up(i,j,k,UFS:UFS+nspec-1,dir) - 0.5d0*dt/dx*(flxd1(i+d(1),j+d(2),k+d(3),UFS:UFS+nspec-1,flxd1c) &
                - flxd1(i,j,k,UFS:UFS+nspec-1,flxd1c)) &
                - 0.5d0*dt/dx*(flxd2(i+d_2(1),j+d_2(2),k+d_2(3),UFS:UFS+nspec-1,flxd2c) &
                                                               - flxd2(i,j,k,UFS:UFS+nspec-1,flxd2c))


           u = uR(i,j,k,UMX,dir)/uR(i,j,k,URHO,dir)
           v = uR(i,j,k,UMY,dir)/uR(i,j,k,URHO,dir)
           w = uR(i,j,k,UMZ,dir)/uR(i,j,k,URHO,dir)
           uR(i,j,k,UEINT,dir) = uR(i,j,k,UEDEN,dir) - 0.5d0*uR(i,j,k,URHO,dir)*(u**2 + v**2 + w**2)

        enddo
     enddo
  enddo

end subroutine

!================================================= Final Magnetic Corrections ========================================================================
subroutine half_step_mag(w_lo, w_hi, &
                         uL, ul_lo, ul_hi, &
                         uR, ur_lo, ur_hi, &
                         um, um_lo, um_hi, &
                         up, up_lo, up_hi, &
                         Ed, ed_lo, ed_hi, &
                         Ed1, ed1_lo, ed1_hi, &
                         Ed2, ed2_lo, ed2_hi, &
                         d, d1, d2, UMAGD, UMAGD1, UMAGD2, sgn, &
                         dx, dt)

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module, only : NVAR, UEINT

  !Correction using Faraday's Law
  implicit none

  integer, intent(in) :: w_lo(3), w_hi(3)
  integer, intent(in) :: ul_lo(3), ul_hi(3)
  integer, intent(in) :: ur_lo(3), ur_hi(3)
  integer, intent(in) :: um_lo(3), um_hi(3)
  integer, intent(in) :: up_lo(3), up_hi(3)
  integer, intent(in) :: ed_lo(3), ed_hi(3)
  integer, intent(in) :: ed1_lo(3), ed1_hi(3)
  integer, intent(in) :: ed2_lo(3), ed2_hi(3)

  integer, intent(in)   :: d, d1, d2, UMAGD, UMAGD1, UMAGD2, sgn

  real(rt), intent(inout) :: uL(ul_lo(1):ul_hi(1), ul_lo(2):ul_hi(2), ul_lo(3):ul_hi(3), NVAR+3,3)
  real(rt), intent(inout) :: uR(ur_lo(1):ur_hi(1), ur_lo(2):ur_hi(2), ur_lo(3):ur_hi(3), NVAR+3,3)
  real(rt), intent(in)    :: um(um_lo(1):um_hi(1), um_lo(2):um_hi(2), um_lo(3):um_hi(3), NVAR+3,3)
  real(rt), intent(in)    :: up(up_lo(1):up_hi(1), up_lo(2):up_hi(2), up_lo(3):up_hi(3), NVAR+3,3)

  real(rt), intent(in) :: Ed(ed_lo(1):ed_hi(1),ed_lo(2):ed_hi(2),ed_lo(3):ed_hi(3))
  real(rt), intent(in) :: Ed1(ed1_lo(1):ed1_hi(1),ed1_lo(2):ed1_hi(2),ed1_lo(3):ed1_hi(3))
  real(rt), intent(in) :: Ed2(ed2_lo(1):ed2_hi(1),ed2_lo(2):ed2_hi(2),ed2_lo(3):ed2_hi(3))

  real(rt), intent(in)  :: dx, dt

  integer :: i ,j ,k
  integer :: a1(3), a2(3), b1(3), b2(3), b3(3), b4(3), b5(3), b6(3) !to manage the +1 shifts on  i,j,k
  integer      :: UMAGX, UMAGY, UMAGZ

  UMAGX = NVAR+1
  UMAGY = NVAR+2
  UMAGZ = NVAR+3

  a1 = 0
  a2 = 0
  b1 = 0
  b2 = 0
  b3 = 0
  b4 = 0
  b5 = 0
  b6 = 0

  !for Bd
  a1(d2) = 1 !shift on first term of Ed1 substraction, in d2 component
  a2(d1) = 1 !shift on first term of Ed2 substraction, in d1 component

  !for Bd1 and Bd2
  b1(d1) = 1 !shift on 1st term for Bd1 and Bd2
  b1(d2) = 1 ! in d1 and d2 components
  b2(d1) = 1 !shift on 2nd and 6th term for Bd1, and 3rd term for Bd2
  b3(d2) = 1 !shift on 3rd term for Bd1, and 2nd and 6th term for Bd2
  b4(d)  = 1 !shift on 5th term for Bd1, on d and d1 components
  b4(d1) = 1
  b5(d)  = 1 !shift on 7th term for Bd1 and Bd2
  b6(d)  = 1 !shift on 5th term for Bd2, on d and d2 components
  b6(d2) = 1


  do k = w_lo(3), w_hi(3)
     do j = w_lo(2), w_hi(2)
        do i = w_lo(1), w_hi(1)
           !---------------------------------------left state-----------------------------------------------------

           !d-Direction
           !Bd eq.45 in Miniati for -
           uL(i,j,k,UMAGD,d) = um(i,j,k,UMAGD,d) + sgn*0.5d0*dt/dx*((Ed1(i+a1(1),j+a1(2),k+a1(3)) - Ed1(i,j,k)) &
                - (Ed2(i+a2(1),j+a2(2),k+a2(3)) - Ed2(i,j,k)))
           !Bd1 eq.46 in Miniati
           uL(i,j,k,UMAGD1,d) = um(i,j,k,UMAGD1,d) + sgn*0.5d0*dt/dx*((Ed(i+b1(1),j+b1(2),k+b1(3)) - Ed(i+b2(1),j+b2(2),k+b2(3))) &
                + (Ed(i+b3(1),j+b3(2),k+b3(3)) - Ed(i,j,k)) &
                - (Ed2(i+b4(1),j+b4(2),k+b4(3)) - Ed2(i+b2(1),j+b2(2),k+b2(3))) &
                - (Ed2(i+b5(1),j+b5(2),k+b5(3)) - Ed2(i,j,k)))
           !Bd2 eq. 46 in Miniati
           uL(i,j,k,UMAGD2,d) = um(i,j,k,UMAGD2,d) - sgn* 0.5d0*dt/dx*((Ed(i+b1(1),j+b1(2),k+b1(3)) - Ed(i+b3(1),j+b3(2),k+b3(3))) &
                + (Ed(i+b2(1),j+b2(2),k+b2(3)) - Ed(i,j,k)) &
                - (Ed1(i+b6(1),j+b6(2),k+b6(3)) - Ed1(i+b3(1),j+b3(2),k+b3(3))) &
                - (Ed1(i+b5(1),j+b5(2),k+b5(3)) - Ed1(i,j,k)))

           uL(i,j,k,UEINT,d) = uL(i,j,k,UEINT,d) - 0.5d0*dot_product(uL(i,j,k,UMAGX:UMAGZ,d),uL(i,j,k,UMAGX:UMAGZ,d))


           !---------------------------------------right state-----------------------------------------------------

           !d-Direction
           !Bd eq. 45 in Miniati for +
           ! for the + case, the shifts mentioned above in b6, b5, and b4
           ! also correspond to the 1st, 2nd and 4th, and 3rd term respectevely
           uR(i,j,k,UMAGD,d) = up(i,j,k,UMAGD,d) + sgn*0.5d0*dt/dx*((Ed1(i+b6(1),j+b6(2),k+b6(3)) - Ed1(i+b5(1),j+b5(2),k+b5(3))) &
                - (Ed2(i+b4(1),j+b4(2),k+b4(3)) - Ed2(i+b5(1),j+b5(2),k+b5(3))))
           !Bd1 eq. 46 in Miniati
           uR(i,j,k,UMAGD1,d) = up(i,j,k,UMAGD1,d) + sgn*0.5d0*dt/dx*((Ed(i+b1(1),j+b1(2),k+b1(3)) - Ed(i+b2(1),j+b2(2),k+b2(3))) &
                + (Ed(i+b3(1),j+b3(2),k+b3(3)) - Ed(i,j,k)) &
                - (Ed2(i+b4(1),j+b4(2),k+b4(3)) - Ed2(i+b2(1),j+b2(2),k+b2(3))) &
                - (Ed2(i+b5(1),j+b5(2),k+b5(3)) - Ed2(i,j,k)))
           !Bd2 eq. 46 in Miniati
           uR(i,j,k,UMAGD2,d) = up(i,j,k,UMAGD2,d) - sgn*0.5d0*dt/dx*((Ed(i+b1(1),j+b1(2),k+b1(3)) - Ed(i+b3(1),j+b3(2),k+b3(3))) &
                + (Ed(i+b2(1),j+b2(2),k+b2(3)) - Ed(i,j,k)) &
                - (Ed1(i+b6(1),j+b6(2),k+b6(3)) - Ed1(i+b3(1),j+b3(2),k+b3(3))) &
                - (Ed1(i+b5(1),j+b5(2),k+b5(3)) - Ed1(i,j,k)))
           uR(i,j,k,UEINT,d) = uR(i,j,k,UEINT,d) -0.5d0*dot_product(uR(i,j,k,UMAGX:UMAGZ,d),uR(i,j,k,UMAGX:UMAGZ,d))

        enddo
     enddo
  enddo

end subroutine half_step_mag

!================================== Find the 2D corrected primitive variables =======================================
subroutine prim_half(w_lo, w_hi, &
                     q2D, q2_lo, q2_hi, &
                     q, q_lo, q_hi, &
                     flxx, flxx_lo, flxx_hi, &
                     flxy, flxy_lo, flxy_hi, &
                     flxz, flxz_lo, flxz_hi, &
                     dx, dy, dz, dt)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module, only : NVAR

 implicit none

 integer, intent(in) :: w_lo(3), w_hi(3)
 integer, intent(in) :: q_lo(3), q_hi(3)
 integer, intent(in) :: q2_lo(3), q2_hi(3)
 integer, intent(in) :: flxx_lo(3), flxx_hi(3)
 integer, intent(in) :: flxy_lo(3), flxy_hi(3)
 integer, intent(in) :: flxz_lo(3), flxz_hi(3)

 real(rt), intent(in)  :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
 real(rt), intent(in)  :: flxx(flxx_lo(1):flxx_hi(1),flxx_lo(2):flxx_hi(2),flxx_lo(3):flxx_hi(3),NVAR+3)
 real(rt), intent(in)  :: flxy(flxy_lo(1):flxy_hi(1),flxy_lo(2):flxy_hi(2),flxy_lo(3):flxy_hi(3),NVAR+3)
 real(rt), intent(in)  :: flxz(flxz_lo(1):flxz_hi(1),flxz_lo(2):flxz_hi(2),flxz_lo(3):flxz_hi(3),NVAR+3)

 real(rt), intent(out) :: q2D(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NQ)

 real(rt)  :: flx_sum(NVAR+3)
 real(rt)  :: qflx(NQ)
 real(rt)  :: dx, dy, dz, dt
 integer   :: i, j, k

 do k = w_lo(3),w_hi(3)
    do j = w_lo(2),w_hi(2)
       do i = w_lo(1),w_hi(1)
          flx_sum = (flxx(i+1,j,k,:) - flxx(i,j,k,:)) + (flxy(i,j+1,k,:) - flxy(i,j,k,:)) + (flxz(i,j,k+1,:) - flxz(i,j,k,:))
          call qflux(qflx,flx_sum,q(i,j,k,:))
          !Right below eq. 48
          q2D(i,j,k,:) = q(i,j,k,:) - 0.5d0*dt/dx*qflx
       enddo
    enddo
 enddo
end subroutine prim_half


!================================= Calculate the C to P Jacobian applied to the fluxes ===================================

subroutine qflux(qflx,flx,q)
 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module !,only : QRHO, QU, QV, QW, QPRES, QMAGX, QMAGY, QMAGZ, QVAR, NVAR
 use eos_module, only : eos
 use eos_type_module, only: eos_t, eos_input_rp
 use network, only : nspec

 implicit none

 ! this is step 10 in the paper, just after Eq. 48

 ! this seems to implement dW/dU . qflux, where dW/dU is the Jacobian of
 ! the primitive quantities (W) with respect to conserved quantities (U)

 real(rt), intent(in)           ::flx(NVAR+3), q(NQ)
 real(rt), intent(out)          ::qflx(NQ)
 real(rt) :: dedp, dedrho, totalE

 integer :: UMAGX, UMAGY, UMAGZ

 type (eos_t) :: eos_state

        UMAGX = NVAR+1
        UMAGY = NVAR+2
        UMAGZ = NVAR+3

        qflx = 0.d0
        qflx(QRHO)  = flx(URHO)
        qflx(QU)    = ( flx(UMX) - flx(URHO) * q(QU) )/q(QRHO)
        qflx(QV)    = ( flx(UMY) - flx(URHO) * q(QV) )/q(QRHO)
        qflx(QW)    = ( flx(UMZ) - flx(URHO) * q(QW) )/q(QRHO)

        qflx(QFS:QFS+nspec-1)  = ( flx(UFS:UFS+nspec-1) - flx(URHO) * q(QFS:QFS+nspec-1) )/q(QRHO)

        eos_state % rho = q(QRHO)
        eos_state % p   = q(QPRES)
        eos_state % T   = 100.d0 !dummy initial guess
        eos_state % xn  = q(QFS:QFS+nspec-1)

        call eos(eos_input_rp, eos_state)

        dedrho  = eos_state % dedr - eos_state % dedT * eos_state % dPdr * 1.0d0/eos_state % dPdT
        dedp    = eos_state % dedT * 1.0d0/eos_state % dPdT

        qflx(QPRES) = ( -q(QMAGX)*flx(UMAGX) - q(QMAGY)*flx(UMAGY) - q(QMAGZ)*flx(UMAGZ) + &
                         flx(UEDEN) - flx(UMX)*q(QU) - flx(UMY)*q(QV) - &
                         flx(UMZ)*q(QW) + flx(URHO)*(0.5*(q(QU)**2+q(QV)**2+q(QW)**2) - &
                         eos_state % e -q(QRHO)*dedrho) ) / ( dedp * q(QRHO))
        qflx(QMAGX) = flx(UMAGX)
        qflx(QMAGY) = flx(UMAGY)
        qflx(QMAGZ) = flx(UMAGZ)

end subroutine qflux


end module ct_upwind
