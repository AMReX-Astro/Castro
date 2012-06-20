
      subroutine ca_ext_src(lo,hi, &
                            old_state,old_state_l1,old_state_h1, &
                            new_state,new_state_l1,new_state_h1, &
                            src,src_l1,src_h1,problo,dx,time,dt)

      use meth_params_module, only : NVAR, UMX
      use probdata_module   , only : center
 
      implicit none
      integer         ,intent(in   ) :: lo(1),hi(1)
      integer         ,intent(in   ) :: old_state_l1,old_state_h1
      integer         ,intent(in   ) :: new_state_l1,new_state_h1
      integer         ,intent(in   ) :: src_l1,src_h1
      double precision,intent(in   ) :: old_state(old_state_l1:old_state_h1,NVAR)
      double precision,intent(in   ) :: new_state(new_state_l1:new_state_h1,NVAR)
      double precision,intent(  out) :: src(src_l1:src_h1,NVAR)
      double precision,intent(in   ) :: problo(1),dx(1),time,dt

      integer          :: i

      src = 0.d0
 
      end subroutine ca_ext_src

