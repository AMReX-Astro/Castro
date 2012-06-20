
     subroutine ca_ext_src(lo,hi,&
                           old_state,old_state_l1,old_state_l2,old_state_l3,old_state_h1,old_state_h2,old_state_h3,&
                           new_state,new_state_l1,new_state_l2,new_state_l3,new_state_h1,new_state_h2,new_state_h3,&
                           src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,problo,dx,time,dt)

     use meth_params_module, only : NVAR
!    use probdata_module   , only : center
 
     implicit none

     integer         ,intent(in   ) :: lo(3),hi(3)
     integer         ,intent(in   ) :: old_state_l1,old_state_l2,old_state_l3,old_state_h1,old_state_h2,old_state_h3
     integer         ,intent(in   ) :: new_state_l1,new_state_l2,new_state_l3,new_state_h1,new_state_h2,new_state_h3
     integer         ,intent(in   ) :: src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
     double precision,intent(in   ) :: old_state(old_state_l1:old_state_h1,old_state_l2:old_state_h2, &
                                                 old_state_l3:old_state_h3,NVAR)
     double precision,intent(in   ) :: new_state(new_state_l1:new_state_h1,new_state_l2:new_state_h2, &
                                                 new_state_l3:new_state_h3,NVAR)
     double precision,intent(  out) :: src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
     double precision,intent(in   ) :: problo(3),dx(3),time,dt

     src = 0.d0
 
     end subroutine ca_ext_src

