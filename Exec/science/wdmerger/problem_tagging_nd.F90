module problem_tagging_module

  implicit none

  public

contains

  subroutine set_problem_tags(lo, hi, &
                              tag, tag_lo, tag_hi, &
                              state, state_lo, state_hi, &
                              dx, problo, &
                              set, clear, time, level) &
                              bind(C,name='set_problem_tags')

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, TWO
    use meth_params_module, only: NVAR, URHO, UTEMP
    use prob_params_module, only: center, probhi, dim, Symmetry, physbc_lo, physbc_hi, &
                                  n_error_buf, ref_ratio, blocking_factor, domlo_level, domhi_level
    use probdata_module, only: max_tagging_radius, &
                               max_stellar_tagging_level, &
                               max_temperature_tagging_level, &
                               max_center_tagging_level, &
                               roche_tagging_factor, &
                               stellar_density_threshold, &
                               temperature_tagging_threshold, &
                               center_tagging_radius, &
                               com_P, com_S, roche_rad_P, roche_rad_S, &
                               problem

    implicit none

    integer,    intent(in   ) :: lo(3), hi(3)
    integer,    intent(in   ) :: tag_lo(3), tag_hi(3)
    integer,    intent(in   ) :: state_lo(3), state_hi(3)
    integer(1), intent(inout) :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
    real(rt),   intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    real(rt),   intent(in   ) :: dx(3), problo(3)
    integer(1), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: level
    real(rt),   intent(in   ), value :: time

    integer  :: i, j, k, n
    real(rt) :: loc(3), r, r_P, r_S

    logical  :: outer_boundary_test(3)
    integer  :: boundary_buf(3), idx(3)

    !$gpu

    do k = lo(3), hi(3)
       loc(3) = problo(3) + (dble(k) + HALF) * dx(3)

       do j = lo(2), hi(2)
          loc(2) = problo(2) + (dble(j) + HALF) * dx(2)

          do i = lo(1), hi(1)
             loc(1) = problo(1) + (dble(i) + HALF) * dx(1)

             if (level < max_stellar_tagging_level) then

                if (problem .eq. 0 .or. problem .eq. 2) then

                   ! For the collision, free-fall, and TDE problems, we just want to tag every
                   ! zone that meets the density criterion; we don't want to bother with
                   ! the Roche lobe radius as that doesn't mean much in these cases.

                   if (state(i,j,k,URHO) > stellar_density_threshold) then

                      tag(i,j,k) = set

                   endif

                else

                   if (level == 0) then

                      ! On the coarse grid, tag all regions within the Roche radii of each star.
                      ! We'll add a buffer around each star to double the Roche
                      ! radius to ensure there aren't any sharp gradients in regions of
                      ! greater than ambient density.

                      r_P = ( sum((loc-com_P)**2) )**HALF
                      r_S = ( sum((loc-com_S)**2) )**HALF

                      if (r_P <= roche_tagging_factor * roche_rad_P) then
                         tag(i,j,k) = set
                      endif

                      if (r_S <= roche_tagging_factor * roche_rad_S) then
                         tag(i,j,k) = set
                      endif

                   else if (level >= 1) then

                      ! On more refined levels, tag all regions within the stars themselves (defined as 
                      ! areas where the density is greater than some threshold).

                      if (state(i,j,k,URHO) > stellar_density_threshold) then

                         tag(i,j,k) = set

                      endif

                   endif

                endif

             endif

             ! Tag all zones at all levels that are hotter than a specified temperature threshold.

             if (level < max_temperature_tagging_level) then

                if (state(i,j,k,UTEMP) > temperature_tagging_threshold) then

                   tag(i,j,k) = set

                endif

             endif

             ! Tag all zones within a specified distance from the center.

             r = ( sum((loc-center)**2) )**HALF

             if (level < max_center_tagging_level) then

                if (r + n_error_buf(level) * minval(dx(1:dim)) <= center_tagging_radius) then

                   tag(i,j,k) = set

                endif

             endif

             ! Clear all tagging that occurs outside the radius set by max_tagging_radius.

             if (r .gt. max_tagging_radius * max(maxval(abs(problo-center)), maxval(abs(probhi-center)))) then

                tag(i,j,k) = clear

             endif

             ! We must ensure that the outermost zones are untagged
             ! due to the Poisson equation boundary conditions.
             ! We currently do not know how to fill the boundary conditions
             ! for fine levels that touch the physical boundary.
             ! (Note that this does not apply for interior/symmetry boundaries.)
             ! To do this properly we need to be aware of AMReX's strategy
             ! for tagging, which is not cell-based, but rather chunk-based.
             ! The size of the chunk on the coarse grid is given by
             ! blocking_factor / ref_ratio -- the idea here being that
             ! blocking_factor is the smallest possible group of cells on a
             ! given level, so the smallest chunk of cells possible on the
             ! coarse grid is given by that number divided by the refinement ratio.
             ! So we cannot tag anything within that distance from the boundary.
             ! Additionally we need to stay a further amount n_error_buf away,
             ! since n_error_buf zones are always added as padding around
             ! tagged zones.

             outer_boundary_test = .false.
             boundary_buf(1:dim) = n_error_buf(level) + blocking_factor(level+1) / ref_ratio(1:dim, level)
             idx = [i, j, k]

             do n = 1, dim

                if ((physbc_lo(n) .ne. Symmetry) .and. (idx(n) .le. domlo_level(n, level) + boundary_buf(n))) then
                   outer_boundary_test(n) = .true.
                endif

                if ((physbc_hi(n) .ne. Symmetry) .and. (idx(n) .ge. domhi_level(n, level) - boundary_buf(n))) then
                   outer_boundary_test(n) = .true.
                endif

             enddo

             if ( any(outer_boundary_test) ) then

                tag(i,j,k) = clear

             endif

          enddo
       enddo
    enddo

  end subroutine set_problem_tags

end module problem_tagging_module
