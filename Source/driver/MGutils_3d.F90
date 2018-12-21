! This is a dummy module - the routines here do nothing and should never be called,
! but we need to create 3d versions of the routines that have been offloaded to
! gpu in the 1d and 2d modules for 3d codes to compile properly.

module MGutils_3D_module

  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

    subroutine ca_weight_cc(lo, hi, &
         cc, clo, chi,  &
         dx, coord_type) bind(C, name="ca_weight_cc")

      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer, intent(in) :: lo(3),hi(3)
      integer, intent(in) :: clo(3), chi(3)
      integer, value, intent(in) :: coord_type
      real(rt), intent(inout) ::cc(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
      real(rt), intent(in) :: dx(3)

      !$gpu

  #ifndef AMREX_USE_CUDA
         print *,'Should not be here '
         call amrex_error("Error:: MGutils_3d.F90 :: ca_weight_cc")
  #endif

    end subroutine ca_weight_cc


  subroutine ca_unweight_cc(lo, hi, &
       cc, clo, chi,  &
       dx, coord_type) bind(C, name="ca_unweight_cc")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3),hi(3)
    integer, intent(in) :: clo(3), chi(3)
    integer, value, intent(in) :: coord_type
    real(rt), intent(inout) ::cc(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in) :: dx(3)

    !$gpu

#ifndef AMREX_USE_CUDA
       print *,'Should not be here '
       call amrex_error("Error:: MGutils_3d.F90 :: ca_unweight_cc")
#endif

  end subroutine ca_unweight_cc


  subroutine ca_unweight_edges(lo, hi, &
       ec, eclo, echi, dx, coord_type, idir) &
       bind(C, name="ca_unweight_edges")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) ::  lo(3),hi(3)
    integer, intent(in) :: eclo(3), echi(3)
    integer, value, intent(in) :: coord_type, idir
    real(rt), intent(inout) :: ec(eclo(1):echi(1),eclo(2):echi(2),eclo(3):echi(3))
    real(rt), intent(in) :: dx(3)

    !$gpu

#ifndef AMREX_USE_CUDA
       print *,'Should not be here '
       call amrex_error("Error:: MGutils_3d.f90 :: ca_unweight_edges")
#endif

  end subroutine ca_unweight_edges


end module MGutils_3D_module
