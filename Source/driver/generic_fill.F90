module generic_fill_module

  implicit none

contains

  ! Used for a generic fill of any StateData.

  subroutine ca_generic_single_fill(state, s_lo, s_hi, &
                                    domlo, domhi, delta, xlo, time, bc) &
                                    bind(C, name="ca_generic_single_fill")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(dim), xlo(dim), time
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    call filcc_nd(state,s_lo,s_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_generic_single_fill




  subroutine ca_generic_multi_fill(state, s_lo, s_hi, &
                                   domlo, domhi, delta, xlo, time, bc) &
                                   bind(C, name="ca_generic_multi_fill")

    use meth_params_module, only: NVAR
    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(dim), xlo(dim), time
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer :: n

    do n = 1, NVAR
       call filcc_nd(state(:,:,:,n),s_lo,s_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    end do

  end subroutine ca_generic_multi_fill

end module generic_fill_module
