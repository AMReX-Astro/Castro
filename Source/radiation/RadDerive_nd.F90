
      subroutine ca_derertot(Et, Et_lo, Et_hi, ncomp_Et, &
                             Er, Er_lo, Er_hi, ncomp_Er, &
                             lo,hi, domlo, domhi, &
                             dx, xlo, time, dt, bc, level, grid_no) bind(C)

        use rad_params_module, only: radtoE
        use amrex_fort_module, only: rt => amrex_real

        implicit none

        integer,  intent(in   ) :: Et_lo(3), Et_hi(3), ncomp_Et
        integer,  intent(in   ) :: Er_lo(3), Er_hi(3), ncomp_Er
        integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
        real(rt), intent(inout) :: Et(Et_lo(1):Et_hi(1),Et_lo(2):Et_hi(2),Et_lo(3):Et_hi(3),ncomp_Et)
        real(rt), intent(in   ) :: Er(Er_lo(1):Er_hi(1),Er_lo(2):Er_hi(2),Er_lo(3):Er_hi(3),ncomp_Er)
        real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
        integer,  intent(in   ) :: bc(3,2,ncomp_Er), level, grid_no
      
        integer :: i, j, k, g

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Et(i,j,k,1) = 0.e0_rt
              end do
           end do
        end do

        do g = 1, ncomp_Er
           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)
                    Et(i,j,k,1) = Et(i,j,k,1) + Er(i,j,k,g)*radtoE
                 end do
              end do
           end do
        end do

      end subroutine ca_derertot

