
       subroutine ca_ext_src(lo,hi, &
                             old_state,old_state_l1,old_state_l2,old_state_h1,old_state_h2,&
                             new_state,new_state_l1,new_state_l2,new_state_h1,new_state_h2,&
                             src,src_l1,src_l2,src_h1,src_h2,problo,dx,time,dt)

       use meth_params_module, only : NVAR, UEDEN, UEINT, URHO
       use probdata_module, only: pert_rad_factor

       implicit none
       integer         , intent(in   ) :: lo(2),hi(2)
       integer         , intent(in   ) :: old_state_l1,old_state_l2,old_state_h1,old_state_h2
       integer         , intent(in   ) :: new_state_l1,new_state_l2,new_state_h1,new_state_h2
       integer         , intent(in   ) :: src_l1,src_l2,src_h1,src_h2
       double precision, intent(in   ) :: old_state(old_state_l1:old_state_h1,old_state_l2:old_state_h2,NVAR)
       double precision, intent(in   ) :: new_state(new_state_l1:new_state_h1,new_state_l2:new_state_h2,NVAR)
       double precision, intent(  out) :: src(    src_l1:  src_h1,  src_l2:src_h2  ,NVAR)
       double precision, intent(in   ) :: problo(2),dx(2),time,dt

       integer          :: i,j
       double precision :: x,y,x1,y1,x2,y2,x3,y3,r1,r2,r3,den,fac

       src(lo(1):hi(1),lo(2):hi(2),:) = 0.d0

       if (time .lt. 0.1d0) then

          if (time + dt .gt. 0.1d0) then
             fac = (0.1d0 - time) / dt
          else
             fac = 1.d0
          end if

          x1 = 5.0d7
          y1 = 6.5d7

          x2 = 1.2d8
          y2 = 8.5d7

          x3 = 2.0d8
          y3 = 7.5d7

          do j=lo(2),hi(2)
             y = problo(2) + dx(2)*(float(j) + 0.5d0)
             do i=lo(1),hi(1)
                x = problo(1) + dx(1)*(float(i) + 0.5d0)
           
                r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / (2.5d6*pert_rad_factor)
                r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / (2.5d6*pert_rad_factor)
                r3 = sqrt( (x-x3)**2 +(y-y3)**2 ) / (2.5d6*pert_rad_factor)

                den = 0.5d0 * (old_state(i,j,URHO) + new_state(i,j,URHO))

                src(i,j,UEDEN) = fac * den * 1.d18 * ( 0.d0*(1.d0 + tanh(2.d0-r1)) &
                                                      +(1.d0 + tanh(2.d0-r2)) &
                                                      +0.d0*(1.d0 + tanh(2.d0-r3)) )
                src(i,j,UEINT) = src(i,j,UEDEN)
                
             end do
          end do

       end if

       end subroutine ca_ext_src
