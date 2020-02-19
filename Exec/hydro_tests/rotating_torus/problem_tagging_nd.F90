
module problem_tagging_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains

  subroutine set_problem_tags(lo, hi, &
                              tag, tag_lo, tag_hi, &
                              state, state_lo, state_hi, &
                              dx, problo, &
                              set, clear, time, level) &
                              bind(C, name='set_problem_tags')

    use meth_params_module, only: NVAR
    use prob_params_module, only: center
    use castro_util_module, only: position ! function
    use probdata_module, only: torus_center, torus_width
    use amrex_fort_module, only : rt => amrex_real
    use iso_c_binding, only : c_int8_t

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: tag_lo(3),tag_hi(3)
    integer,    intent(in   ) :: state_lo(3), state_hi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    real(rt),   intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    real(rt),   intent(in   ) :: dx(3), problo(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: level
    real(rt),   intent(in   ), value :: time

    integer  :: i, j, k
    real(rt) :: loc(3), R, Z

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             R = sqrt(loc(1)**2 + loc(2)**2)

             Z = loc(3)

             if ( (torus_center - R)**2 + Z**2 < torus_width**2 ) then
                tag(i,j,k) = set
             else
                tag(i,j,k) = clear
             endif

          enddo
       enddo
    enddo

  end subroutine set_problem_tags

end module problem_tagging_module
