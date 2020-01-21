module rad_source_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_rad_source(lo, hi, &
                           rhs, rhs_lo, rhs_hi, &
                           dx, dt, time, igroup) &
                           bind(C, name="ca_rad_source")

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3) 
    integer,  intent(in   ) :: rhs_lo(3), rhs_hi(3)
    real(rt), intent(inout) :: rhs(rhs_lo(1):rhs_hi(1),rhs_lo(2):rhs_hi(2),rhs_lo(3):rhs_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, time
    integer,  intent(in   ), value :: igroup

    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! This is a dummy routine that does nothing, but could be implemented
             ! by a problem setup to inject a radiation source to group igroup.
             ! This routine must add to rhs.

          end do
       end do
    end do

  end subroutine ca_rad_source

end module rad_source_module
