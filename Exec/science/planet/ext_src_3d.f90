  ! compute the external sources for all the conservative equations.
  ! 
  ! This is called twice in the evolution:
  ! 
  ! first, for the predictor, it is called with (old, old) states
  !
  ! This is also used in the first pass of the conservative update
  ! (adding dt * S there)
  !
  ! Next we correct the source terms in the conservative update to
  ! time-center them.  Here we call ext_src(old, new), and then
  ! in time_center_source_terms we subtract off 1/2 of the first S
  ! and add 1/2 of the new S.
  !
  ! Therefore, to get a properly time-centered source, generally
  ! speaking, you always want to use the "new" state here.  That
  ! will be the time n state in the first call and the n+1 in the
  ! second call.
  !
  
  subroutine ca_ext_src(lo,hi,&
                        old_state,old_state_l1,old_state_l2,old_state_l3, &
                        old_state_h1,old_state_h2,old_state_h3,&
                        new_state,new_state_l1,new_state_l2,new_state_l3, &
                        new_state_h1,new_state_h2,new_state_h3,&
                        src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,problo,dx,time,dt)
    use amrex_constants_module, only : ZERO, TWO, HALF
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, NVAR, UTEMP, UFS, UFX, UEINT
    use prob_params_module, only: center, probhi
    use actual_network, only: nspec, naux
    use eos_type_module, only: eos_t, eos_input_rt, eos_input_re
    use eos_module, only: eos

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer         ,intent(in   ) :: lo(3),hi(3)
    integer         ,intent(in   ) :: old_state_l1,old_state_l2,old_state_l3, &
         old_state_h1,old_state_h2,old_state_h3
    integer         ,intent(in   ) :: new_state_l1,new_state_l2,new_state_l3, &
         new_state_h1,new_state_h2,new_state_h3
    integer         ,intent(in   ) :: src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    real(rt)        ,intent(in   ) :: old_state(old_state_l1:old_state_h1, &
         old_state_l2:old_state_h2, &
         old_state_l3:old_state_h3,NVAR)
    real(rt)        ,intent(in   ) :: new_state(new_state_l1:new_state_h1, &
         new_state_l2:new_state_h2, &
         new_state_l3:new_state_h3,NVAR)
    real(rt)        ,intent(  out) :: src(src_l1:src_h1, &
         src_l2:src_h2, &
         src_l3:src_h3,NVAR)
    real(rt)        ,intent(in   ) :: problo(3),dx(3),time,dt
    real(rt) :: source(3)
    integer  :: i, j, k
    real(rt) :: ang_vel, y, beta, p_orb_vel, p_radius,reset_center
    p_orb_vel = 1D-3
    p_radius = 1D10
    source = ZERO
    ang_vel = ZERO
    beta = p_orb_vel * TWO / p_radius
    reset_center = (problo(2) + probhi(2)) * HALF
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          y = problo(2) + dx(2)*(dble(j) + HALF) - reset_center
          do i = lo(1), hi(1)
             ang_vel = 0.5 * beta * y
       
             source(1) = + TWO * new_state(i,j,k,UMY) * ang_vel
             source(2) = - TWO * new_state(i,j,k,UMX) * ang_vel
             source(3) = ZERO
             src(i,j,k,UMX) = src(i,j,k,UMX) + source(1)
             src(i,j,k,UMY) = src(i,j,k,UMY) + source(2)
             src(i,j,k,UMZ) = src(i,j,k,UMZ) + source(3)

          end do
        end do
     end do


  end subroutine ca_ext_src
