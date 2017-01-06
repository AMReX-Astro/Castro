! Problem-specific Fortran routines that are designed to interact with C++

subroutine problem_checkpoint(int_dir_name, len) bind(c)

  ! called by the IO processor during checkpoint

  use bl_fort_module, only : rt => c_real
  implicit none

  integer :: len
  integer :: int_dir_name(len)

end subroutine problem_checkpoint



subroutine problem_restart(int_dir_name, len) bind(c)

  ! called by ALL processors during restart 

  use bl_fort_module, only : rt => c_real
  implicit none

  integer :: len
  integer :: int_dir_name(len)

end subroutine problem_restart



! Return the problem type.

subroutine get_problem_number(problem_out) bind(C,name='get_problem_number')

  use probdata_module, only: problem

  use bl_fort_module, only : rt => c_real
  implicit none

  integer :: problem_out

  problem_out = problem

end subroutine get_problem_number



! Return the diameter.

subroutine get_diameter(diameter_out) bind(C,name='get_diameter')

  use probdata_module, only: diameter

  use bl_fort_module, only : rt => c_real
  implicit none

  real(rt)         :: diameter_out

  diameter_out = diameter

end subroutine get_diameter



! Return the density.

subroutine get_density(density_out) bind(C,name='get_density')

  use probdata_module, only: density

  use bl_fort_module, only : rt => c_real
  implicit none

  real(rt)         :: density_out

  density_out = density

end subroutine get_density



! Update the density field. This ensures
! that the sum of the mass on the domain
! is what we intend it to be.

subroutine update_density(lo, hi, dx, &
                          state, s_lo, s_hi, &
                          update_factor) bind(C, name='update_density')

  use bl_constants_module, only: HALF
  use network, only: nspec
  use meth_params_module, only: NVAR, URHO, UFS
  use prob_params_module, only: problo, center
  use probdata_module, only: problem, diameter

  use bl_fort_module, only : rt => c_real
  implicit none

  integer         , intent(in   ) :: lo(3), hi(3)
  integer         , intent(in   ) :: s_lo(3), s_hi(3)

  real(rt)        , intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

  real(rt)        , intent(in   ) :: dx(3)
  real(rt)        , intent(in   ) :: update_factor

  integer          :: i, j, k
  real(rt)         :: xx, yy, zz
  
  if (problem .eq. 2) then

     do k = lo(3), hi(3)
        zz = problo(3) + dx(3) * (dble(k)+HALF) - center(3)

        do j = lo(2), hi(2)
           yy = problo(2) + dx(2) * (dble(j)+HALF) - center(2)

           do i = lo(1), hi(1)
              xx = problo(1) + dx(1) * (dble(i)+HALF) - center(1)

              if ((xx**2 + yy**2 + zz**2)**0.5 < diameter / 2) then

                 state(i,j,k,URHO) = state(i,j,k,URHO) * update_factor
                 state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,UFS:UFS-1+nspec) * update_factor

              endif

           enddo
        enddo
     enddo

  endif

end subroutine update_density
