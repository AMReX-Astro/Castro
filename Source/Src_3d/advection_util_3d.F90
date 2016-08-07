module advection_util_3d_module

  implicit none

  private

  public normalize_species_fluxes, divu, &
         limit_hydro_fluxes_on_small_dens

contains

  subroutine normalize_species_fluxes(flux1,flux1_lo,flux1_hi, &
                                      flux2,flux2_lo,flux2_hi, &
                                      flux3,flux3_lo,flux3_hi, &
                                      lo, hi)

    ! here we normalize the fluxes of the mass fractions so that
    ! they sum to 0.  This is essentially the CMA procedure that is
    ! defined in Plewa & Muller, 1999, A&A, 342, 179

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    double precision, intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    double precision, intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    double precision, intent(inout) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)

    ! Local variables
    integer          :: i, j, k, n
    double precision :: sum, fac

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux1(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux1(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux1(i,j,k,n) = flux1(i,j,k,n) * fac
             end do
          end do
       end do
    end do

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux2(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux2(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux2(i,j,k,n) = flux2(i,j,k,n) * fac
             end do
          end do
       end do
    end do

    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux3(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux3(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux3(i,j,k,n) = flux3(i,j,k,n) * fac
             end do
          end do
       end do
    end do

  end subroutine normalize_species_fluxes

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine limit_hydro_fluxes_on_small_dens(state,state_lo,state_hi, &
                                              flux1,flux1_lo,flux1_hi, &
                                              flux2,flux2_lo,flux2_hi, &
                                              flux3,flux3_lo,flux3_hi, &
                                              vol,  vol_lo,  vol_hi,   &
                                              lo,hi)

    use bl_constants_module, only: ZERO, ONE
    use meth_params_module, only: NVAR, URHO, small_dens
    use mempool_module, only: bl_allocate, bl_deallocate

    implicit none

    integer, intent(in) :: state_lo(3), state_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    integer, intent(in) :: vol_lo(3),   vol_hi(3)
    integer, intent(in) :: lo(3),       hi(3)

    double precision, intent(in   ) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    double precision, intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    double precision, intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    double precision, intent(inout) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    double precision, intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))

    double precision, pointer :: limit1(:,:,:), limit2(:,:,:), limit3(:,:,:)

    integer          :: i, j, k, n
    double precision :: drho, drho_new, limit_factor, volInv
    double precision, parameter :: neg_rho_eps_factor = 1.d-12

    call bl_allocate(limit1,flux1_lo(1),flux1_hi(1),flux1_lo(2),flux1_hi(2),flux1_lo(3),flux1_hi(3))
    call bl_allocate(limit2,flux2_lo(1),flux2_hi(1),flux2_lo(2),flux2_hi(2),flux2_lo(3),flux2_hi(3))
    call bl_allocate(limit3,flux3_lo(1),flux3_hi(1),flux3_lo(2),flux3_hi(2),flux3_lo(3),flux3_hi(3))

    limit1 = ONE
    limit2 = ONE
    limit3 = ONE

    ! The approach we are using is the following: we check every zone
    ! and ask if the sum of the fluxes going into that zone that would
    ! cause the density to decrease would cause the density to become
    ! less than small_dens. If so, we multiply every negative flux going
    ! into that zone so that the total sum of the negative fluxes yields
    ! at least small_dens when applied (and we'll add a small factor on
    ! top of this to ensure that we don't hit roundoff level issues in later
    ! comparisons to small_dens). We cannot do this approach using all of
    ! the fluxes (that is, including the positive ones) and try to hit
    ! small_dens exactly because an adjacent zone might also have limiting and
    ! limit the positive fluxes such that there isn't enough material coming in
    ! to stay above small_dens, so limiting only the negative fluxes seems
    ! like a reasonable compromise. To make this independent of the order in
    ! which we check the fluxes, we store these changes in a multiplicative
    ! array with the same indexing as the fluxes, and only multiply the
    ! fluxes in a second step.

    ! This is NOT thread-safe, so if we ever do away with tiling and decide
    ! to do this routine using explicit loop-level OpenMP, do not thread this
    ! loop. Also, if the limiting occurs on the edge of a box, it may create
    ! a flux which is inconsistent with the corresponding flux from the box
    ! which is adjacent to it. At present I do not know how to deal with this.
    ! This routine is still worthwhile because it is less non-conservative than
    ! the approach of negative density resets (this subroutine is conservative
    ! everywhere away from box edges) and because negative density resets can
    ! introduce some ugliness depending on how severe they are.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             volInv = ONE / vol(i,j,k)

             drho = ZERO

             if (flux1(i,j,k,URHO) < ZERO) then
                drho = drho + flux1(i,j,k,URHO) * volInv
             endif

             if (flux1(i+1,j,k,URHO) > ZERO) then
                drho = drho - flux1(i+1,j,k,URHO) * volInv
             endif

             if (flux2(i,j,k,URHO) < ZERO) then
                drho = drho + flux2(i,j,k,URHO) * volInv
             endif

             if (flux2(i,j+1,k,URHO) > ZERO) then
                drho = drho - flux2(i,j+1,k,URHO) * volInv
             endif

             if (flux3(i,j,k,URHO) < ZERO) then
                drho = drho + flux3(i,j,k,URHO) * volInv
             endif

             if (flux3(i,j,k+1,URHO) > ZERO) then
                drho = drho - flux3(i,j,k+1,URHO) * volInv
             endif

             if (state(i,j,k,URHO) < small_dens) then

                limit_factor = ZERO

             elseif (state(i,j,k,URHO) + drho < small_dens) then

                drho_new = state(i,j,k,URHO) - small_dens
                limit_factor = abs(drho_new / drho) * (ONE - neg_rho_eps_factor)

             else

                limit_factor = ONE

             endif

             if (flux1(i,j,k,URHO) < ZERO) then
                limit1(i,j,k) = limit1(i,j,k) * limit_factor
             endif

             if (flux1(i+1,j,k,URHO) > ZERO) then
                limit1(i+1,j,k) = limit1(i+1,j,k) * limit_factor
             endif

             if (flux2(i,j,k,URHO) < ZERO) then
                limit2(i,j,k) = limit2(i,j,k) * limit_factor
             endif

             if (flux2(i,j+1,k,URHO) > ZERO) then
                limit2(i,j+1,k) = limit2(i,j+1,k) * limit_factor
             endif

             if (flux3(i,j,k,URHO) < ZERO) then
                limit3(i,j,k) = limit3(i,j,k) * limit_factor
             endif

             if (flux3(i,j,k+1,URHO) > ZERO) then
                limit3(i,j,k+1) = limit3(i,j,k+1) * limit_factor
             endif

          enddo
       enddo
    enddo

    ! Now limit every flux by this factor we have computed.
    ! This relies on every flux being proportional to the
    ! density. The only state variable this is not true for
    ! is the temperature, and the flux is zero for that one
    ! anyway. However if another state variable is added that
    ! is not linear in density then this may need adjustment.

    do n = 1, NVAR

       flux1(:,:,:,n) = flux1(:,:,:,n) * limit1(:,:,:)
       flux2(:,:,:,n) = flux2(:,:,:,n) * limit2(:,:,:)
       flux3(:,:,:,n) = flux3(:,:,:,n) * limit3(:,:,:)

    enddo

    call bl_deallocate(limit1)
    call bl_deallocate(limit2)
    call bl_deallocate(limit3)

  end subroutine limit_hydro_fluxes_on_small_dens

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine divu(lo,hi,q,q_lo,q_hi,dx,div,div_lo,div_hi)

    use meth_params_module, only : QU, QV, QW, QVAR
    use bl_constants_module

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: div_lo(3), div_hi(3)
    double precision, intent(in) :: dx(3)
    double precision, intent(inout) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    double precision, intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    integer          :: i, j, k
    double precision :: ux, vy, wz, dxinv, dyinv, dzinv

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1

             ux = FOURTH*( &
                    + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                    + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                    + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                    + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) ) * dxinv

             vy = FOURTH*( &
                    + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                    + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                    + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                    + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) ) * dyinv

             wz = FOURTH*( &
                    + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                    + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                    + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                    + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) ) * dzinv

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo

  end subroutine divu

end module advection_util_3d_module
