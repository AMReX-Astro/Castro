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
  
  subroutine ca_ext_src(lo,hi, &
                        old_state,old_state_l1,old_state_h1, &
                        new_state,new_state_l1,new_state_h1, &
                        src,src_l1,src_h1,problo,dx,time,dt)

    use meth_params_module, only : NVAR, NSRC
    use amrex_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer         ,intent(in   ) :: lo(1),hi(1)
    integer         ,intent(in   ) :: old_state_l1,old_state_h1
    integer         ,intent(in   ) :: new_state_l1,new_state_h1
    integer         ,intent(in   ) :: src_l1,src_h1
    real(rt)        ,intent(in   ) :: old_state(old_state_l1:old_state_h1,NVAR)
    real(rt)        ,intent(in   ) :: new_state(new_state_l1:new_state_h1,NVAR)
    real(rt)        ,intent(  out) :: src(src_l1:src_h1,NSRC)
    real(rt)        ,intent(in   ) :: problo(1),dx(1),time,dt

    ! lo and hi specify work region
    src(lo(1):hi(1),:) = ZERO ! Fill work region only

  end subroutine ca_ext_src

