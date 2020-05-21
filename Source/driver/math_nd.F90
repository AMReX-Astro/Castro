module math_module

  public

contains

  function cross_product(A,B) result(C)
    ! Compute the standard cross-product of two three-vectors.

    use amrex_fort_module, only : rt => amrex_real
    use prob_params_module, only : coord_type

    implicit none

    real(rt)         :: A(3), B(3)
    real(rt)         :: C(3)

    !$gpu

    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1) - A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)

    ! for axisymmetry, we store the coordinates as (r, z, theta),
    ! since the simulation plane is r-z.  But this is a left-handed
    ! system.  However because of the ordering of our state, we get a
    ! Coriolis force that has the correct signs for a right-handed
    ! system, so no change is needed.

  end function cross_product

end module math_module
