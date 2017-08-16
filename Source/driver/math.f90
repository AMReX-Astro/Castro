module math_module

  public

contains

  ! Compute the standard cross-product of two three-vectors.

  function cross_product(A,B) result(C)

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: A(3), B(3)
    real(rt)         :: C(3)

    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1) - A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)

  end function cross_product

end module math_module
