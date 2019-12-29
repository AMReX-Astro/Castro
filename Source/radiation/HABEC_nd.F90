module habec_nd_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine set_abec_flux(lo, hi, &
                           dir, &
                           density, d_lo, d_hi, &
                           dcoef, c_lo, c_hi, &
                           beta, &
                           dx, &
                           flux, f_lo, f_hi) &
                           bind(C, name="set_abec_flux")

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: density(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    real(rt), intent(in   ) :: dcoef(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    integer,  intent(in   ), value :: dir
    real(rt), intent(in   ), value :: beta
    real(rt), intent(in   ) :: dx(3)

    integer  :: i, j, k
    real(rt) :: fac

    !$gpu

    if (dir == 0) then

       fac = -beta / dx(1)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux(i,j,k) = dcoef(i,j,k) * (density(i,j,k) - density(i-1,j,k)) * fac
             end do
          end do
       end do

    else if (dir == 1) then

       fac = -beta / dx(2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux(i,j,k) = dcoef(i,j,k) * (density(i,j,k) - density(i,j-1,k)) * fac
             end do
          end do
       end do

    else

       fac = -beta / dx(3)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux(i,j,k) = dcoef(i,j,k) * (density(i,j,k) - density(i,j,k-1)) * fac
             end do
          end do
       end do

    end if

  end subroutine set_abec_flux

end module habec_nd_module
