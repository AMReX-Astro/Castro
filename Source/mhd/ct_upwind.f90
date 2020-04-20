module ct_upwind

 use amrex_fort_module, only : rt => amrex_real
 use hlld_solver, only : hlld
 use meth_params_module

 implicit none 

 private primtocons
 public corner_transport
 
 
interface checkisnan
       module procedure checkisnanmult
       module procedure checkisnans
end interface checkisnan

 contains

subroutine corner_transport( q, qm, qp, q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3, &      
                                flxx, flxx_l1 , flxx_l2 , flxx_l3 , flxx_h1 , flxx_h2 , flxx_h3, &
                                flxy, flxy_l1 , flxy_l2 , flxy_l3 , flxy_h1 , flxy_h2 , flxy_h3, &
                                flxz, flxz_l1 , flxz_l2 , flxz_l3 , flxz_h1 , flxz_h2 , flxz_h3, &
                                Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                                Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                                Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                                dx, dy, dz, dt)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module, only : NQ
 use electric_field
implicit none

        integer, intent(in)   :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3

        integer, intent(in)   :: ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3
        integer, intent(in)   :: ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3
        integer, intent(in)   :: ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3

        integer, intent(in)   :: flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
        integer, intent(in)   :: flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3
        integer, intent(in)   :: flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3

        real(rt), intent(in)  :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ) !prim vars at time t^n
        real(rt), intent(in)  :: qm(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ,3)
        real(rt), intent(in)  :: qp(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ,3)

        ! fluxes should be NVAR+3
        real(rt), intent(out) :: flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,NVAR+3)   !Half Step Fluxes
        real(rt), intent(out) :: flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,NVAR+3)   !Half Step Fluxes
        real(rt), intent(out) :: flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,NVAR+3)   !Half Step Fluxes

        real(rt), intent(out)  :: Ex(ex_l1:ex_h1,ex_l2:ex_h2,ex_l3:ex_h3)
        real(rt), intent(out)  :: Ey(ey_l1:ey_h1,ey_l2:ey_h2,ey_l3:ey_h3)
        real(rt), intent(out)  :: Ez(ez_l1:ez_h1,ez_l2:ez_h2,ez_l3:ez_h3)
 
        ! these are conserved + magnetic field (cell centered)
 
        real(rt)  :: um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3) !PtoC Vars
        real(rt)  :: up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)
        real(rt)  :: cons_temp_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3,2) !2D Temporary Conservative Vars
        real(rt)  :: cons_temp_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3,2) !2D Temporary Conservative Vars
        real(rt)  :: cons_half_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3) !Flux Corrected Conservative Vars
        real(rt)  :: cons_half_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)

        real(rt)  :: q_temp_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ,3,2) !2D Temporary Primitive Vars
        real(rt)  :: q_temp_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ,3,2) !2D Temporary Primitive Vars
        real(rt)  :: q_half_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ,3) !Flux Corrected Primitive Vars
        real(rt)  :: q_half_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ,3) !Flux Corrected Primitive Vars


        real(rt)  :: flxx1D(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,NVAR+3)
        real(rt)  :: flxy1D(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,NVAR+3) 
        real(rt)  :: flxz1D(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,NVAR+3) !Flux1d for all directions

        real(rt)  :: flxx2D(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,NVAR+3, 2) !Flux2d for all directions 2 perpendicular directions
        real(rt)  :: flxy2D(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,NVAR+3, 2) !Flux2d for all directions 2 perpendicular directions
        real(rt)  :: flxz2D(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,NVAR+3, 2) !Flux2d for all directions 2 perpendicular directions

        real(rt)  :: q2D(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ)
        real(rt)  :: dx, dy, dz, dt

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

!Calculate Flux 1D
        !x-dir
        work_lo(1) = flxx_l1
        work_lo(2) = flxx_l2
        work_lo(3) = flxx_l3
        work_hi(1) = flxx_h1
        work_hi(2) = flxx_h2
        work_hi(3) = flxx_h3
        
        call hlld(work_lo, work_hi, qm,qp,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                  flxx1D(:,:,:,:),flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, 1)



        work_lo(1) = flxy_l1
        work_lo(2) = flxy_l2
        work_lo(3) = flxy_l3
        work_hi(1) = flxy_h1
        work_hi(2) = flxy_h2
        work_hi(3) = flxy_h3

        !y-dir  
        call hlld(work_lo, work_hi, qm,qp,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                  flxy1D(:,:,:,:),flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, 2)



        work_lo(1) = flxz_l1
        work_lo(2) = flxz_l2
        work_lo(3) = flxz_l3
        work_hi(1) = flxz_h1
        work_hi(2) = flxz_h2
        work_hi(3) = flxz_h3
                  
        !z-dir
        call hlld(work_lo, work_hi, qm,qp,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                  flxz1D(:,:,:,:),flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, 3)


!Prim to Cons
do i = 1,3
        call PrimToCons(qm(:,:,:,:,i), um(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
        call PrimToCons(qp(:,:,:,:,i), up(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
enddo


!Use "1D" fluxes To interpolate Temporary Edge Centered Electric Fields


        work_lo(1) = ex_l1+1
        work_lo(2) = ex_l2+1
        work_lo(3) = ex_l3+1
        work_hi(1) = ex_h1-1
        work_hi(2) = ex_h2-1
        work_hi(3) = ex_h3-1
        call electric_edge_x(work_lo, work_hi, &
                             q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                             Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                             flxy1D, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                             flxz1D, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

        work_lo(1) = ey_l1+1
        work_lo(2) = ey_l2+1
        work_lo(3) = ey_l3+1
        work_hi(1) = ey_h1-1
        work_hi(2) = ey_h2-1
        work_hi(3) = ey_h3-1
        call electric_edge_y(work_lo, work_hi, &
                             q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                             Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                             flxx1D, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                             flxz1D, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

        work_lo(1) = ez_l1+1
        work_lo(2) = ez_l2+1
        work_lo(3) = ez_l3+1
        work_hi(1) = ez_h1-1
        work_hi(2) = ez_h2-1
        work_hi(3) = ez_h3-1
        call electric_edge_z(work_lo, work_hi, &
                             q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                             Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                             flxx1D, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                             flxy1D, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3)

!Corner Couple
        
        work_lo(1) = q_l1+1
        work_lo(2) = q_l2+1
        work_lo(3) = q_l3+1
        work_hi(1) = q_h1-1
        work_hi(2) = q_h2-1
        work_hi(3) = q_h3-1
        
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                           flxy1D, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                           !d1 = x, d2 = y, dir2 = 1 last component in the cons_temp arrays 
                           1 , 2, 1, &    !y corrected x  
                           dx, dt) !Correct Conservative vars using Transverse Fluxes
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                           flxz1D, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                           1, 3, 2, &    !z corrected x  
                           dx, dt)
 
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                           flxx1D, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                           2, 1, 1, &    !x corrected y  
                           dy, dt) 
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                           flxz1D, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                           2, 3, 2, &    !z corrected y 
                           dy, dt) 

        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                           flxx1D, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                           3, 1, 1, &    !x corrected z  
                           dz, dt) 
        call corner_couple(work_lo, work_hi, &
                           cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                           flxy1D, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                           3, 2, 2, &    !y corrected z 
                           dz, dt) 


        
        !X direction
        ! affected by Y Flux                                                    
        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                               Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                               Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                              !x,y,z,dir2=1, sgn=+, UMAGD1, UMAGD2, UMAGD3  
                               1, 2, 3, 1, 1, UMAGX, UMAGY, UMAGZ, dx, dt) !Correct Conservative vars using Transverse Fluxes

      
       ! affected by Z Flux                                                      
        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                               Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                               Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &  
                               1, 3, 2, 2, -1, UMAGX, UMAGZ, UMAGY, dx, dt)   
       
        !Y direction
        ! affected by X Flux                                                    
        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                               Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                               Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &  
                               2, 1, 3, 1, -1, UMAGY, UMAGX, UMAGZ, dy, dt) 

      
       ! affected by Z Flux                                                      
        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                               Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                               Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                               2, 3, 1, 2, 1, UMAGY, UMAGZ, UMAGX, dy, dt) 

        !Z direction
        ! affected by X Flux                                                    
        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                               Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                               Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &  
                               3, 1, 2, 1, 1, UMAGZ, UMAGX, UMAGY, dz, dt) 

      
       ! affected by Y Flux                                                      
        call corner_couple_mag(work_lo, work_hi, &
                               cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                               Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                               Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                               3, 2, 1, 2, -1, UMAGZ, UMAGY, UMAGX, dz, dt) 


!Cons To Prim
do i = 1,3
        call ConsToPrim(q_temp_M(:,:,:,:,i,1), cons_temp_M(:,:,:,:,i,1), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
        call ConsToPrim(q_temp_P(:,:,:,:,i,1), cons_temp_P(:,:,:,:,i,1), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
        call ConsToPrim(q_temp_M(:,:,:,:,i,2), cons_temp_M(:,:,:,:,i,2), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
        call ConsToPrim(q_temp_P(:,:,:,:,i,2), cons_temp_P(:,:,:,:,i,2), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
enddo

!Calculate Flux 2D
do i = 1,2
        work_lo(1) = flxx_l1+1
        work_lo(2) = flxx_l2+1
        work_lo(3) = flxx_l3+1
        work_hi(1) = flxx_h1-1
        work_hi(2) = flxx_h2-1
        work_hi(3) = flxx_h3-1
        !x-dir
        call hlld(work_lo, work_hi, q_temp_M(:,:,:,:,:,i),q_temp_P(:,:,:,:,:,i),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                  flxx2D(:,:,:,:,i),flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, 1)

        !y-dir  
        work_lo(1) = flxy_l1+1
        work_lo(2) = flxy_l2+1
        work_lo(3) = flxy_l3+1
        work_hi(1) = flxy_h1-1
        work_hi(2) = flxy_h2-1
        work_hi(3) = flxy_h3-1

        call hlld(work_lo, work_hi, q_temp_M(:,:,:,:,:,i),q_temp_P(:,:,:,:,:,i),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                  flxy2D(:,:,:,:,i),flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, 2)

        !z-dir

        work_lo(1) = flxz_l1+1
        work_lo(2) = flxz_l2+1
        work_lo(3) = flxz_l3+1
        work_hi(1) = flxz_h1-1
        work_hi(2) = flxz_h2-1
        work_hi(3) = flxz_h3-1
        
        call hlld(work_lo, work_hi, q_temp_M(:,:,:,:,:,i),q_temp_P(:,:,:,:,:,i),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                  flxz2D(:,:,:,:,i),flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, 3)

enddo

!Use Averaged 2D fluxes to interpolate temporary Edge Centered Electric Fields, reuse "flx1D"
        flxx1D(:,:,:,:) = 0.5d0*(flxx2D(:,:,:,:,1) + flxx2D(:,:,:,:,2))
        flxy1D(:,:,:,:) = 0.5d0*(flxy2D(:,:,:,:,1) + flxy2D(:,:,:,:,2))
        flxz1D(:,:,:,:) = 0.5d0*(flxz2D(:,:,:,:,1) + flxz2D(:,:,:,:,2))


        work_lo(1) = ex_l1+1!+2
        work_lo(2) = ex_l2+1!+2
        work_lo(3) = ex_l3+1!+2
        work_hi(1) = ex_h1-1!-2
        work_hi(2) = ex_h2-1!-2
        work_hi(3) = ex_h3-1!-2
        call electric_edge_x(work_lo, work_hi, &
                         q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                             Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                             flxy1D, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                             flxz1D, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

        work_lo(1) = ey_l1+1!+2
        work_lo(2) = ey_l2+1!+2
        work_lo(3) = ey_l3+1!+2
        work_hi(1) = ey_h1-1!-2
        work_hi(2) = ey_h2-1!-2
        work_hi(3) = ey_h3-1!-2
        call electric_edge_y(work_lo, work_hi, &
                         q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                             Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                             flxx1D, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                             flxz1D, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

        work_lo(1) = ez_l1+1!+2
        work_lo(2) = ez_l2+1!+2
        work_lo(3) = ez_l3+1!+2
        work_hi(1) = ez_h1-1!-2
        work_hi(2) = ez_h2-1!-2
        work_hi(3) = ez_h3-1!-2
        call electric_edge_z(work_lo, work_hi, &
                         q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                             Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                             flxx1D, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                             flxy1D, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3)

!Half Step conservative vars
        work_lo(1) = q_l1+1
        work_lo(2) = q_l2+1
        work_lo(3) = q_l3+1
        work_hi(1) = q_h1-2
        work_hi(2) = q_h2-2
        work_hi(3) = q_h3-2
                                     
        
        !for x direction   
        call half_step(work_lo, work_hi, &
                       cons_half_M, cons_half_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                       flxy2D, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                       flxz2D, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                       !dir = x, d1 =y, d2 =z
                       1, 2, 3, dx, dt)
        
        !for y direction
        call half_step(work_lo, work_hi, &
                       cons_half_M, cons_half_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                       flxx2D, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                       flxz2D, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                       2, 1, 3, dy, dt)

        !for z direction
        call half_step(work_lo, work_hi, &
                       cons_half_M, cons_half_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                       flxx2D, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                       flxy2D, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                       3, 1, 2, dz, dt)
 

        !x direction 
        call half_step_mag(work_lo, work_hi, &
                           cons_half_M, cons_half_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                           Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                           Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                           Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                           !d=x, d1=y, d2=z, UMAGD UMAGD1, UMAGD2, sgn,
                           1, 2, 3, UMAGX, UMAGY, UMAGZ, -1, &
                           dx, dt)
        !y direction
        call half_step_mag(work_lo, work_hi, &
                           cons_half_M, cons_half_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                           Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                           Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                           Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                           !d, d1, d2, UMAGD UMAGD1, UMAGD2, sgn,
                           2, 1, 3, UMAGY, UMAGX, UMAGZ, 1, &
                           dy, dt)

        !z direction
        call half_step_mag(work_lo, work_hi, &
                           cons_half_M, cons_half_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                           Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                           Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                           Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                           !d, d1, d2, UMAGD UMAGD1, UMAGD2, sgn,
                           3, 1, 2, UMAGZ, UMAGX, UMAGY, -1, &
                           dz, dt)


do i = 1,3
        call ConsToPrim(q_half_M(:,:,:,:,i), cons_half_M(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
        call ConsToPrim(q_half_P(:,:,:,:,i), cons_half_P(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
enddo

!Final Fluxes


        !x-dir
        work_lo(1) = flxx_l1+2
        work_lo(2) = flxx_l2+2
        work_lo(3) = flxx_l3+2
        work_hi(1) = flxx_h1-2
        work_hi(2) = flxx_h2-2
        work_hi(3) = flxx_h3-2

        call hlld(work_lo, work_hi, q_half_M,q_half_P,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                   flxx,flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, 1)

        !y-dir  
        work_lo(1) = flxy_l1+2
        work_lo(2) = flxy_l2+2
        work_lo(3) = flxy_l3+2
        work_hi(1) = flxy_h1-2
        work_hi(2) = flxy_h2-2
        work_hi(3) = flxy_h3-2
        
        call hlld(work_lo, work_hi, q_half_M,q_half_P,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                   flxy,flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, 2)

        !z-dir
        work_lo(1) = flxz_l1+2
        work_lo(2) = flxz_l2+2
        work_lo(3) = flxz_l3+2
        work_hi(1) = flxz_h1-2
        work_hi(2) = flxz_h2-2
        work_hi(3) = flxz_h3-2
        
        call hlld(work_lo, work_hi, q_half_M,q_half_P,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                  flxz,flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, 3)
        
!Primitive update
       call prim_half(q2D,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                      flxx1D, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                      flxy1D, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                      flxz1D, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                      dx, dy, dz, dt)

!Final Electric Field Update


        work_lo(1) = ex_l1+2!+3
        work_lo(2) = ex_l2+2!+3
        work_lo(3) = ex_l3+2!+3
        work_hi(1) = ex_h1-2!-3
        work_hi(2) = ex_h2-2!-3
        work_hi(3) = ex_h3-2!-3
        call electric_edge_x(work_lo, work_hi, &
                         q2D, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                                             Ex, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                                             flxy, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                                             flxz, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

        work_lo(1) = ey_l1+2!+3
        work_lo(2) = ey_l2+2!+3
        work_lo(3) = ey_l3+2!+3
        work_hi(1) = ey_h1-2!-3
        work_hi(2) = ey_h2-2!-3
        work_hi(3) = ey_h3-2!-3
        call electric_edge_y(work_lo, work_hi, &
                         q2D, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                                             Ey, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                                             flxx, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                                             flxz, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3)

        work_lo(1) = ez_l1+2!+3
        work_lo(2) = ez_l2+2!+3
        work_lo(3) = ez_l3+2!+3
        work_hi(1) = ez_h1-2!-3
        work_hi(2) = ez_h2-2!-3
        work_hi(3) = ez_h3-2!-3
        call electric_edge_z(work_lo, work_hi, &
                         q2D, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                                             Ez, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                                             flxx, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                                             flxy, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3)

end subroutine corner_transport

!================================================= Calculate the Conservative Variables ===============================================

subroutine PrimToCons(q, u, q_l1 ,q_l2 ,q_l3 ,q_h1 ,q_h2 ,q_h3)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module
 use eos_type_module, only : eos_t, eos_input_rp
 use eos_module, only: eos
 use network, only: nspec

implicit none

        integer, intent(in)             ::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
        real(rt), intent(in)    ::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ)
        real(rt), intent(out)   ::u(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3)
        integer                                 :: i ,j ,k

  type(eos_t) :: eos_state        

 do k = q_l3,q_h3
    do j = q_l2,q_h2
       do i = q_l1, q_h1
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

subroutine ConsToPrim(q, u, q_l1 ,q_l2 ,q_l3 ,q_h1 ,q_h2 ,q_h3)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module
 use eos_module, only : eos
 use eos_type_module, only : eos_t, eos_input_re
 use network, only : nspec

 implicit none

 integer, intent(in)            ::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
 real(rt), intent(in)   ::u(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3)
 real(rt), intent(out)  ::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ)
 integer                                        :: i ,j ,k
 integer                :: UMAGX, UMAGZ

 type (eos_t) :: eos_state

 UMAGX = NVAR+1
 UMAGZ = NVAR+3
 !q = u
 do k = q_l3,q_h3
    do j = q_l2,q_h2
       do i = q_l1, q_h1
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
                         uL, uR, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                         flxd2, flxd2_l1,flxd2_l2,flxd2_l3,flxd2_h1,flxd2_h2,flxd2_h3, &
                         d1, d2, dir2, &
                         dx, dt)

use amrex_fort_module, only : rt => amrex_real
use meth_params_module
use network, only : nspec

implicit none
        
        integer, intent(in)      :: w_lo(3), w_hi(3)
        integer, intent(in)      :: q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
        integer, intent(in)      :: flxd2_l1,flxd2_l2,flxd2_l3,flxd2_h1,flxd2_h2,flxd2_h3
        integer, intent(in)      :: d1, d2, dir2
        real(rt), intent(in)     :: dx, dt
        
        real(rt), intent(in)    ::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)
        real(rt), intent(in)    ::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)

        real(rt), intent(out) :: flxd2(flxd2_l1:flxd2_h1,flxd2_l2:flxd2_h2,flxd2_l3:flxd2_h3,NVAR+3)

        real(rt), intent(out)   ::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3,2)
        real(rt), intent(out)   ::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3,2)

        real(rt)                :: u, v, w
        integer                 :: i ,j ,k
        integer                 :: d(3) !for the addition of +1 to either i,j,k depending on d2    

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
                             uL, uR, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                             Ed1, ed1_l1,ed1_l2,ed1_l3,ed1_h1,ed1_h2,ed1_h3, &
                             Ed3, ed3_l1,ed3_l2,ed3_l3,ed3_h1,ed3_h2,ed3_h3, &
                             d1, d2, d3, dir2, sgn, UMAGD1, UMAGD2, UMAGD3, &           
                             dx, dt)
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : NVAR, UEINT

!Correction using Faraday's Law
implicit none

        integer, intent(in)     :: w_lo(3), w_hi(3)
        integer, intent(in)     :: q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
        integer, intent(in)     :: ed1_l1,ed1_l2,ed1_l3,ed1_h1,ed1_h2, ed1_h3
        integer, intent(in)     :: ed3_l1,ed3_l2,ed3_l3,ed3_h1,ed3_h2, ed3_h3
        integer, intent(in)     :: d1, d2, d3, dir2, sgn, UMAGD1, UMAGD2, UMAGD3   !UMAGD1 corresponds to d1, and UMAGD2 to dir2 
        real(rt), intent(inout) :: uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3,2)
        real(rt), intent(inout) :: uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3,2)
        real(rt), intent(in)    :: um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)
        real(rt), intent(in)    :: up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)   

        real(rt), intent(in)    :: Ed1(ed1_l1:ed1_h1,ed1_l2:ed1_h2,ed1_l3:ed1_h3)
        real(rt), intent(in)    :: Ed3(ed3_l1:ed3_h1,ed3_l2:ed3_h2,ed3_l3:ed3_h3)

        real(rt)                :: dx, dt
        integer                 :: i ,j ,k
        integer                 :: d(3), a1(3), a2(3), a3(3), d_2(3) !for the additions of +1 to i,j,k 
        integer                 ::UMAGX, UMAGY, UMAGZ

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
                         uL, uR, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
                     flxd1, flxd1_l1,flxd1_l2,flxd1_l3,flxd1_h1,flxd1_h2,flxd1_h3, &
                     flxd2, flxd2_l1,flxd2_l2,flxd2_l3,flxd2_h1,flxd2_h2,flxd2_h3, &
                     dir, d1, d2, dx, dt)
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UMZ, URHO, UEINT, UFS
use network, only: nspec

implicit none
        
        integer, intent(in)   :: w_lo(3), w_hi(3),q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
        integer, intent(in)   :: flxd1_l1,flxd1_l2,flxd1_l3,flxd1_h1,flxd1_h2,flxd1_h3
        integer, intent(in)   :: flxd2_l1,flxd2_l2,flxd2_l3,flxd2_h1,flxd2_h2,flxd2_h3
        
        real(rt), intent(in)  ::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)
        real(rt), intent(in)  ::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)

        real(rt), intent(in)  :: dx, dt !dx will be dx, dy or dz
        integer, intent(in)   :: dir, d1, d2 ! following notation of eq. 44 

        real(rt), intent(in)  :: flxd1(flxd1_l1:flxd1_h1,flxd1_l2:flxd1_h2,flxd1_l3:flxd1_h3,NVAR+3,2)
        real(rt), intent(in)  :: flxd2(flxd2_l1:flxd2_h1,flxd2_l2:flxd2_h2,flxd2_l3:flxd2_h3,NVAR+3,2)

        real(rt), intent(out)   ::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)
        real(rt), intent(out)   ::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)

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
                          uL, uR, um, up, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                          Ed, ed_l1,ed_l2,ed_l3,ed_h1,ed_h2,ed_h3, &
                          Ed1, ed1_l1,ed1_l2,ed1_l3,ed1_h1,ed1_h2,ed1_h3, &
                          Ed2, ed2_l1,ed2_l2,ed2_l3,ed2_h1,ed2_h2,ed2_h3, &
                          d, d1, d2, UMAGD, UMAGD1, UMAGD2, sgn, &
                          dx, dt)
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : NVAR, UEINT

!Correction using Faraday's Law
implicit none

        integer, intent(in)   :: w_lo(3), w_hi(3), q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
        integer, intent(in)   :: ed_l1,ed_l2,ed_l3,ed_h1,ed_h2,ed_h3
        integer, intent(in)   :: ed1_l1,ed1_l2,ed1_l3,ed1_h1,ed1_h2,ed1_h3
        integer, intent(in)   :: ed2_l1,ed2_l2,ed2_l3,ed2_h1,ed2_h2,ed2_h3
        integer, intent(in)   :: d, d1, d2, UMAGD, UMAGD1, UMAGD2, sgn

        real(rt), intent(inout)   ::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)
        real(rt), intent(inout)   ::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)
        real(rt), intent(in)      ::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)
        real(rt), intent(in)      ::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NVAR+3,3)    

        real(rt), intent(in)  :: Ed(ed_l1:ed_h1,ed_l2:ed_h2,ed_l3:ed_h3)
        real(rt), intent(in)  :: Ed1(ed1_l1:ed1_h1,ed1_l2:ed1_h2,ed1_l3:ed1_h3)
        real(rt), intent(in)  :: Ed2(ed2_l1:ed2_h1,ed2_l2:ed2_h2,ed2_l3:ed2_h3)

        real(rt), intent(in)  :: dx, dt

        integer               :: i ,j ,k
        integer               :: a1(3), a2(3), b1(3), b2(3), b3(3), b4(3), b5(3), b6(3) !to manage the +1 shifts on  i,j,k  
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
subroutine prim_half(q2D,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                     flxx, flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                     flxy, flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                     flxz, flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                     dx, dy, dz, dt)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module, only : NVAR

implicit none

        integer, intent(in)             ::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3

        integer, intent(in)   :: flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
        integer, intent(in)   :: flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3
        integer, intent(in)   :: flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3

        real(rt), intent(in)    :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ)
        real(rt), intent(in)    :: flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,NVAR+3)
        real(rt), intent(in)    :: flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,NVAR+3)
        real(rt), intent(in)    :: flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,NVAR+3)

        real(rt), intent(out)   ::q2D(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ)

        real(rt)                                ::flx_sum(NVAR+3)
        real(rt)                                ::qflx(NQ)
        real(rt)                                :: dx, dy, dz, dt       
        integer                                 ::i, j, k

        do k = q_l3+1,q_h3-1
                do j = q_l2+1,q_h2-1
                        do i = q_l1+1,q_h1-1
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


!============================================ Debug code =====================================================
        subroutine checkisnanmult(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, num)
           use amrex_fort_module, only : rt => amrex_real

        implicit none
        integer, intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, num
        real(rt), intent(in) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,num)

        integer :: i,j,k,n


        do n = 1,num
                do k = uout_l3,uout_h3
                        do j = uout_l2, uout_h2
                                do i = uout_l1,uout_h1
                                        if(isnan(uout(i,j,k,n)).or.(abs(uout(i,j,k,n)).ge. 1d14)) then
                                                write(*,*) "Bad values ",  uout(i,j,k,:)
                                                write(*,*) "Failure to converge ", "i, j, k, n = ", i, j, k, n
                                                stop
                                        endif
                                enddo
                        enddo
                enddo
        enddo
        end subroutine checkisnanmult
!============ single =====================      

        subroutine checkisnans(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
           use amrex_fort_module, only : rt => amrex_real

        implicit none
        integer, intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
        real(rt), intent(in) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3)

        integer :: i,j,k

                do k = uout_l3,uout_h3
                        do j = uout_l2, uout_h2
                                do i = uout_l1,uout_h1
                                        if(isnan(uout(i,j,k)).or.(abs(uout(i,j,k)).ge. 1d14)) then
                                                write(*,*) "Bad values ",  uout(i,j,k)
                                                write(*,*) "Failure to converge ", "i, j, k = ", i, j, k
                                                stop
                                        endif
                                enddo
                        enddo
                enddo
        end subroutine checkisnans

!====================================== Density Check ========================================
        subroutine checknegdens(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
           use amrex_fort_module, only : rt => amrex_real

        implicit none
        integer, intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
        real(rt), intent(in) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3)

        integer :: i,j,k

                do k = uout_l3,uout_h3
                        do j = uout_l2, uout_h2
                                do i = uout_l1,uout_h1
                                        if(uout(i,j,k).le. 0.d0) then
                                                write(*,*) "Non-Positive Density ",  uout(i,j,k)
                                                write(*,*) "i, j, k = ", i, j, k
                                                stop
                                        endif
                                enddo
                        enddo
                enddo
        end subroutine checknegdens
end module ct_upwind
