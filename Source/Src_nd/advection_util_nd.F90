module advection_util_module

  implicit none

  private

  public enforce_minimum_density, compute_cfl, ctoprim, srctoprim

contains

  subroutine enforce_minimum_density(uin,uin_lo,uin_hi, &
                                     uout,uout_lo,uout_hi, &
                                     vol,vol_lo,vol_hi, &
                                     lo,hi,mass_added,eint_added, &
                                     eden_added,frac_change,verbose) &
                                     bind(C, name="enforce_minimum_density")
    
    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UEINT, UEDEN, small_dens, density_reset_method
    use bl_constants_module, only : ZERO

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) ::  uin_lo(3),  uin_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) ::  vol_lo(3),  vol_hi(3)

    double precision, intent(in) ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
    double precision, intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    double precision, intent(in) ::  vol( vol_lo(1): vol_hi(1), vol_lo(2): vol_hi(2), vol_lo(3): vol_hi(3))
    double precision, intent(inout) :: mass_added, eint_added, eden_added, frac_change
    
    ! Local variables
    integer          :: i,ii,j,jj,k,kk
    integer          :: i_set, j_set, k_set
    double precision :: max_dens
    double precision :: unew(NVAR)
    integer          :: num_positive_zones
    
    double precision :: initial_mass, final_mass
    double precision :: initial_eint, final_eint
    double precision :: initial_eden, final_eden

    initial_mass = ZERO
      final_mass = ZERO

    initial_eint = ZERO
      final_eint = ZERO

    initial_eden = ZERO
      final_eden = ZERO

    max_dens = ZERO

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             initial_mass = initial_mass + uout(i,j,k,URHO ) * vol(i,j,k)
             initial_eint = initial_eint + uout(i,j,k,UEINT) * vol(i,j,k)
             initial_eden = initial_eden + uout(i,j,k,UEDEN) * vol(i,j,k)

             if (uout(i,j,k,URHO) .eq. ZERO) then

                print *,'DENSITY EXACTLY ZERO AT CELL ',i,j,k
                print *,'  in grid ',lo(1),lo(2),lo(3),hi(1),hi(2),hi(3)
                call bl_error("Error:: advection_util_nd.f90 :: enforce_minimum_density")

             else if (uout(i,j,k,URHO) < small_dens) then

                ! Store the maximum (negative) fractional change in the density

                if ( uout(i,j,k,URHO) < ZERO .and. &
                     (uout(i,j,k,URHO) - uin(i,j,k,URHO)) / uin(i,j,k,URHO) < frac_change) then

                   frac_change = (uout(i,j,k,URHO) - uin(i,j,k,URHO)) / uin(i,j,k,URHO)

                endif

                if (density_reset_method == 1) then

                   ! Reset to the characteristics of the adjacent state with the highest density.
                   
                   max_dens = uout(i,j,k,URHO)
                   i_set = i
                   j_set = j
                   k_set = k
                   do kk = -1,1
                      do jj = -1,1
                         do ii = -1,1
                            if (i+ii.ge.lo(1) .and. j+jj.ge.lo(2) .and. k+kk.ge.lo(3) .and. &
                                 i+ii.le.hi(1) .and. j+jj.le.hi(2) .and. k+kk.le.hi(3)) then
                               if (uout(i+ii,j+jj,k+kk,URHO) .gt. max_dens) then
                                  i_set = i+ii
                                  j_set = j+jj
                                  k_set = k+kk
                                  max_dens = uout(i_set,j_set,k_set,URHO)
                               endif
                            endif
                         end do
                      end do
                   end do

                   if (max_dens < small_dens) then

                      ! We could not find any nearby zones with sufficient density.

                      call reset_to_small_state(uin(i,j,k,:), uout(i,j,k,:), [i, j, k], verbose)

                   else

                      unew = uout(i_set,j_set,k_set,:)

                      call reset_to_zone_state(uin(i,j,k,:), uout(i,j,k,:), unew(:), [i, j, k], verbose)

                   endif

                else if (density_reset_method == 2) then

                   ! Reset to the average of adjacent zones. The median is independently calculated for each variable.

                   num_positive_zones = 0
                   unew(:) = ZERO

                   do kk = -1, 1
                      do jj = -1, 1
                         do ii = -1, 1
                            if (i+ii.ge.lo(1) .and. j+jj.ge.lo(2) .and. k+kk.ge.lo(3) .and. &
                                i+ii.le.hi(1) .and. j+jj.le.hi(2) .and. k+kk.le.hi(3)) then
                               if (uout(i+ii,j+jj,k+kk,URHO) .ge. small_dens) then
                                  unew(:) = unew(:) + uout(i+ii,j+jj,k+kk,:)
                                  num_positive_zones = num_positive_zones + 1
                               endif
                            endif
                         enddo
                      enddo
                   enddo

                   if (num_positive_zones == 0) then

                      ! We could not find any nearby zones with sufficient density.

                      call reset_to_small_state(uin(i,j,k,:), uout(i,j,k,:), [i, j, k], verbose)

                   else

                      unew(:) = unew(:) / num_positive_zones

                      call reset_to_zone_state(uin(i,j,k,:), uout(i,j,k,:), unew(:), [i, j, k], verbose)

                   endif

                elseif (density_reset_method == 3) then

                   ! Reset to the original zone state.

                   if (uin(i,j,k,URHO) < small_dens) then

                      call reset_to_small_state(uin(i,j,k,:), uout(i,j,k,:), [i, j, k], verbose)

                   else

                      unew(:) = uin(i,j,k,:)

                      call reset_to_zone_state(uin(i,j,k,:), uout(i,j,k,:), unew(:), [i, j, k], verbose)

                   endif

                else

                   call bl_error("Unknown density_reset_method in subroutine enforce_minimum_density.")

                endif

             end if

             final_mass = final_mass + uout(i,j,k,URHO ) * vol(i,j,k)
             final_eint = final_eint + uout(i,j,k,UEINT) * vol(i,j,k)
             final_eden = final_eden + uout(i,j,k,UEDEN) * vol(i,j,k)

          enddo
       enddo
    enddo

    if ( max_dens /= ZERO ) then
       mass_added = mass_added + final_mass - initial_mass
       eint_added = eint_added + final_eint - initial_eint
       eden_added = eden_added + final_eden - initial_eden
    endif

  end subroutine enforce_minimum_density



  subroutine reset_to_small_state(old_state, new_state, idx, verbose)

    use bl_constants_module, only: ZERO
    use network, only: nspec
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEINT, UEDEN, UFS, small_temp, small_dens, npassive, upass_map
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_module, only: eos
    use castro_util_module, only: position
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only: linear_to_hybrid
    use meth_params_module, only: UMR, UMP
#endif

    implicit none

    double precision :: old_state(NVAR), new_state(NVAR)
    integer          :: idx(3), verbose

    integer          :: n, ipassive
    type (eos_t)     :: eos_state

#ifdef HYBRID_MOMENTUM
    double precision :: loc(3)
#endif

    ! If no neighboring zones are above small_dens, our only recourse 
    ! is to set the density equal to small_dens, and the temperature 
    ! equal to small_temp. We set the velocities to zero, 
    ! though any choice here would be arbitrary.
    
    if (verbose .gt. 0) then
       print *,'   '
       if (new_state(URHO) < ZERO) then
          print *,'>>> RESETTING NEG.  DENSITY AT ',idx(1),idx(2),idx(3)
       else
          print *,'>>> RESETTING SMALL DENSITY AT ',idx(1),idx(2),idx(3)
       endif
       print *,'>>> FROM ',new_state(URHO),' TO ',small_dens
       print *,'>>> ORIGINAL DENSITY FOR OLD STATE WAS ',old_state(URHO)
       print *,'   '
    end if

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       new_state(n) = new_state(n) * (small_dens / new_state(URHO))
    end do

    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = new_state(UFS:UFS+nspec-1) / small_dens

    call eos(eos_input_rt, eos_state)

    new_state(URHO ) = eos_state % rho
    new_state(UTEMP) = eos_state % T

    new_state(UMX  ) = ZERO
    new_state(UMY  ) = ZERO
    new_state(UMZ  ) = ZERO

    new_state(UEINT) = eos_state % rho * eos_state % e
    new_state(UEDEN) = new_state(UEINT)

#ifdef HYBRID_MOMENTUM
    loc = position(idx(1),idx(2),idx(3))
    new_state(UMR:UMP) = linear_to_hybrid(loc, new_state(UMX:UMZ))
#endif

  end subroutine reset_to_small_state


 
  subroutine reset_to_zone_state(old_state, new_state, input_state, idx, verbose)

    use bl_constants_module, only: ZERO
    use meth_params_module, only: NVAR, URHO

    implicit none

    double precision :: old_state(NVAR), new_state(NVAR), input_state(NVAR)
    integer          :: idx(3), verbose

    if (verbose .gt. 0) then
       if (new_state(URHO) < ZERO) then
          print *,'   '
          print *,'>>> RESETTING NEG.  DENSITY AT ',idx(1),idx(2),idx(3)
          print *,'>>> FROM ',new_state(URHO),' TO ',input_state(URHO)
          print *,'>>> ORIGINAL DENSITY FOR OLD STATE WAS ',old_state(URHO)
          print *,'   '
       else
          print *,'   '
          print *,'>>> RESETTING SMALL DENSITY AT ',idx(1),idx(2),idx(3)
          print *,'>>> FROM ',new_state(URHO),' TO ',input_state(URHO)
          print *,'>>> ORIGINAL DENSITY FOR OLD STATE WAS ',old_state(URHO)
          print *,'   '
       end if
    end if

    new_state(:) = input_state(:)

  end subroutine reset_to_zone_state



  subroutine compute_cfl(q, q_lo, q_hi, lo, hi, dt, dx, courno) &
                         bind(C, name = 'ca_compute_cfl')

    use bl_constants_module, only: ZERO, ONE
    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QC
    use prob_params_module, only: dim

    implicit none

    integer :: lo(3), hi(3)
    integer :: q_lo(3), q_hi(3)

    double precision :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    double precision :: dt, dx(3), courno

    double precision :: courx, coury, courz, courmx, courmy, courmz
    double precision :: dtdx, dtdy, dtdz
    integer          :: i, j, k

    ! Compute running max of Courant number over grids

    courmx = courno
    courmy = courno
    courmz = courno

    dtdx = dt / dx(1)

    if (dim .ge. 2) then
       dtdy = dt / dx(2)
    else
       dtdy = ZERO
    endif

    if (dim .eq. 3) then
       dtdz = dt / dx(3)
    else
       dtdz = ZERO
    endif

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             courx = ( q(i,j,k,QC) + abs(q(i,j,k,QU)) ) * dtdx
             coury = ( q(i,j,k,QC) + abs(q(i,j,k,QV)) ) * dtdy
             courz = ( q(i,j,k,QC) + abs(q(i,j,k,QW)) ) * dtdz

             courmx = max( courmx, courx )
             courmy = max( courmy, coury )
             courmz = max( courmz, courz )

             if (courx .gt. ONE) then
                print *,'   '
                call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ca_compute_cfl")
                print *,'>>> ... (u+c) * dt / dx > 1 ', courx
                print *,'>>> ... at cell (i,j,k)   : ', i, j, k
                print *,'>>> ... u, c                ', q(i,j,k,QU), q(i,j,k,QC)
                print *,'>>> ... density             ', q(i,j,k,QRHO)
             end if

             if (coury .gt. ONE) then
                print *,'   '
                call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ca_compute_cfl")
                print *,'>>> ... (v+c) * dt / dx > 1 ', coury
                print *,'>>> ... at cell (i,j,k)   : ', i,j,k
                print *,'>>> ... v, c                ', q(i,j,k,QV), q(i,j,k,QC)
                print *,'>>> ... density             ', q(i,j,k,QRHO)
             end if

             if (courz .gt. ONE) then
                print *,'   '
                call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ca_compute_cfl")
                print *,'>>> ... (w+c) * dt / dx > 1 ', courz
                print *,'>>> ... at cell (i,j,k)   : ', i, j, k
                print *,'>>> ... w, c                ', q(i,j,k,QW), q(i,j,k,QC)
                print *,'>>> ... density             ', q(i,j,k,QRHO)
             end if

          enddo
       enddo
    enddo

    courno = max( courmx, courmy, courmz )

  end subroutine compute_cfl



  subroutine ctoprim(lo, hi, &
                     uin, uin_lo, uin_hi, &
                     q,     q_lo,   q_hi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use actual_network, only : nspec, naux
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_re
    use meth_params_module, only : NVAR, URHO, UMX, UMZ, &
                                   UEDEN, UEINT, UTEMP, &
                                   QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QTEMP, QGAME, QFS, QFX, &
                                   QC, QCSML, QGAMC, QDPDR, QDPDE, &
                                   npassive, upass_map, qpass_map, dual_energy_eta1, &
                                   allow_negative_energy
    use bl_constants_module, only: ZERO, HALF, ONE
    use castro_util_module, only: position

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)

    double precision, intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    double precision, intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)

    double precision, parameter :: small = 1.d-8

    integer          :: i, j, k
    integer          :: n, nq, ipassive
    double precision :: kineng, rhoinv
    double precision :: loc(3), vel(3)

    type (eos_t) :: eos_state

    !
    ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
    ! The temperature is used as an initial guess for the eos call and will be overwritten.
    !
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (uin(i,j,k,URHO) .le. ZERO) then
                print *,'   '
                print *,'>>> Error: Castro_advection_3d::ctoprim ',i,j,k
                print *,'>>> ... negative density ',uin(i,j,k,URHO)
                call bl_error("Error:: Castro_advection_3d.f90 :: ctoprim")
             end if
          end do

          do i = lo(1),hi(1)

             loc = position(i,j,k)

             q(i,j,k,QRHO) = uin(i,j,k,URHO)
             rhoinv = ONE/q(i,j,k,QRHO)

             vel = uin(i,j,k,UMX:UMZ) * rhoinv

             q(i,j,k,QU:QW) = vel

             ! Get the internal energy, which we'll use for determining the pressure.
             ! We use a dual energy formalism. If (E - K) < eta1 and eta1 is suitably small,
             ! then we risk serious numerical truncation error in the internal energy.
             ! Therefore we'll use the result of the separately updated internal energy equation.
             ! Otherwise, we'll set e = E - K.

             kineng = HALF * q(i,j,k,QRHO) * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)

             if ( (uin(i,j,k,UEDEN) - kineng) / uin(i,j,k,UEDEN) .gt. dual_energy_eta1) then
                q(i,j,k,QREINT) = (uin(i,j,k,UEDEN) - kineng) * rhoinv
             else
                q(i,j,k,QREINT) = uin(i,j,k,UEINT) * rhoinv
             endif

             q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)

          enddo
       enddo
    enddo

    ! Load passively advected quatities into q
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
             enddo
          enddo
       enddo
    enddo

    if (allow_negative_energy .eq. 0) eos_state % reset = .true.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state % T   = q(i,j,k,QTEMP )
             eos_state % rho = q(i,j,k,QRHO  )
             eos_state % e   = q(i,j,k,QREINT)
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)
             eos_state % aux = q(i,j,k,QFX:QFX+naux-1)

             call eos(eos_input_re, eos_state)

             q(i,j,k,QTEMP)  = eos_state % T
             q(i,j,k,QREINT) = eos_state % e
             q(i,j,k,QPRES)  = eos_state % p
             q(i,j,k,QDPDR)  = eos_state % dpdr_e
             q(i,j,k,QDPDE)  = eos_state % dpde
             q(i,j,k,QC)     = eos_state % cs
             q(i,j,k,QGAMC)  = eos_state % gam1
             q(i,j,k,QCSML)  = max(small, small * q(i,j,k,QC))
             q(i,j,k,QREINT) = q(i,j,k,QREINT) * q(i,j,k,QRHO)
             q(i,j,k,QGAME)  = q(i,j,k,QPRES) / q(i,j,k,QREINT) + ONE

          enddo
       enddo
    enddo

  end subroutine ctoprim



  subroutine srctoprim(lo, hi, &
                       q,     q_lo,   q_hi, &
                       src, src_lo, src_hi, &
                       srcQ,srQ_lo, srQ_hi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use actual_network, only : nspec, naux
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_re
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, &
                                   QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QDPDR, QDPDE, &
                                   npassive, upass_map, qpass_map
    use bl_constants_module, only: ZERO, HALF, ONE
    use castro_util_module, only: position

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: srQ_lo(3), srQ_hi(3)

    double precision, intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    double precision, intent(in   ) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    double precision, intent(inout) :: srcQ(srQ_lo(1):srQ_hi(1),srQ_lo(2):srQ_hi(2),srQ_lo(3):srQ_hi(3),QVAR)

    double precision, parameter :: small = 1.d-8

    integer          :: i, j, k
    integer          :: n, nq, ipassive
    double precision :: rhoinv

    srcQ = ZERO

    ! compute srcQ terms
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoinv = ONE / q(i,j,k,QRHO)

             srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
             srcQ(i,j,k,QU    ) = (src(i,j,k,UMX) - q(i,j,k,QU) * srcQ(i,j,k,QRHO)) * rhoinv
             srcQ(i,j,k,QV    ) = (src(i,j,k,UMY) - q(i,j,k,QV) * srcQ(i,j,k,QRHO)) * rhoinv
             srcQ(i,j,k,QW    ) = (src(i,j,k,UMZ) - q(i,j,k,QW) * srcQ(i,j,k,QRHO)) * rhoinv
             srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) &
                                                   - q(i,j,k,QV)*src(i,j,k,UMY) &
                                                   - q(i,j,k,QW)*src(i,j,k,UMZ) &
                                    + HALF * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2) * srcQ(i,j,k,QRHO)

             srcQ(i,j,k,QPRES ) = q(i,j,k,QDPDE)*(srcQ(i,j,k,QREINT) - &
                                  q(i,j,k,QREINT)*srcQ(i,j,k,QRHO)*rhoinv) * rhoinv + &
                                  q(i,j,k,QDPDR)*srcQ(i,j,k,QRHO)

          enddo
       enddo
    enddo

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                srcQ(i,j,k,nq) = ( src(i,j,k,n) - q(i,j,k,nq) * srcQ(i,j,k,QRHO) ) / &
                                 q(i,j,k,QRHO)
             enddo
          enddo
       enddo

    enddo

  end subroutine srctoprim

end module advection_util_module
