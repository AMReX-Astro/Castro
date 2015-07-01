
subroutine ca_tempdiffextrap(lo, hi, tdif, t_l1, t_h1)
  integer :: lo(1), hi(1)
  integer :: t_l1, t_h1
  double precision :: tdif(t_l1:t_h1)

  tdif(lo(1)-1) = tdif(lo(1))
  tdif(hi(1)+1) = tdif(hi(1))

end subroutine ca_tempdiffextrap
