
  subroutine interpolate_fluid (lo, hi, &
                                S , s_lo, s_hi, ns , &
                                n_fluid_dof, ng) 

    use amrex_constants_module, only : fourth, half, zero, one, two
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEDEN,UFX,UTEMP
    use FluidFieldsModule, only : uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use EquationOfStateModule_TABLE, only : ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
                                            ComputeThermodynamicStates_Primitive_TABLE
    use UnitsModule, only : Gram, Centimeter, Second, AtomicMassUnit, Erg, Kelvin

    use ReferenceElementModuleX, only: NodesX_q

    implicit none
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3),  s_hi(3)
    integer, intent(in) :: ns
    integer, intent(in) ::  n_fluid_dof
    integer, intent(in) :: ng

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ns) 

    ! Temporary variables
    integer  :: i,j,k,n
    integer  :: ic,jc,kc
    integer  :: ind
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne

    ! For interpolation
    real(rt) :: xslope(ns), yslope(ns), zslope(ns)
    real(rt) :: xyslope(ns), xzslope(ns), yzslope(ns)
    real(rt) :: Sval(ns)
    real(rt) :: dlft(ns), drgt(ns), dcen(ns), dlim, dsgn

    real(rt) :: rho_in(1)
    real(rt) :: Ye_in(1)
    real(rt) :: E_in(1)
    real(rt) :: T_in(1)

    real(rt) :: T_out(1)
    real(rt) :: Ne_out(1)
    real(rt) :: Epervol_out(1)
    real(rt) :: Epermass_out(1)

    real(rt) :: fac(ns),fac_all
    real(rt) :: delta_S_slope(ns)
    real(rt) :: delta_S_val(ns)
    real(rt) :: S_lll(ns),S_llh(ns),S_lhl(ns),S_lhh(ns)
    real(rt) :: S_hll(ns),S_hlh(ns),S_hhl(ns),S_hhh(ns)

    real(rt) :: x(n_fluid_dof)
    real(rt) :: y(n_fluid_dof)
    real(rt) :: z(n_fluid_dof)

    integer  :: interp_type

    ! No interpolation
    ! interp_type = 0 

    ! Use conserved variables for interpolation
    ! interp_type = 1 

    ! Use primitive variables (rho, E, Ye) for interpolation
    ! interp_type = 2 

    ! Use primitive variables (rho, T, Ye) for interpolation
    ! interp_type = 3

    ! Use primitive variables (rho, T, Ye) for interpolation -- first fill the corners
    !     of the cell then use bilinear interpolation
    interp_type = 4

    ! Conversion to thornado units
    conv_dens = Gram / Centimeter**3
    conv_mom  = Gram / Centimeter**2 / Second
    conv_enr  = Erg / Centimeter**3
    conv_ne   = 1.d0 / Centimeter**3

    ! These are the locations of the DG nodes in the space [-.5:.5]
    ! Note that in 3D with 2 nodes in each direction, they are ordered 
    ! ind = 1:  (Lo,Lo,Lo)
    ! ind = 2:  (Hi,Lo,Lo)
    ! ind = 3:  (Lo,Hi,Lo)
    ! ind = 4:  (Hi,Hi,Lo)
    ! ind = 5:  (Lo,Lo,Hi)
    ! ind = 6:  (Hi,Lo,Hi)
    ! ind = 7:  (Lo,Hi,Hi)
    ! ind = 8:  (Hi,Hi,Hi)
     do ind = 1, n_fluid_dof
        x(ind) = NodesX_q(1,ind)
        y(ind) = NodesX_q(2,ind)
        z(ind) = NodesX_q(3,ind)
     end do

     ! The uCF array was allocated in CreateFluidFieldsConserved with 
     !     ALLOCATE( uCF &
     !      (1:nDOFX, &
     !       1-swX(1):nX(1)+swX(1), &
     !       1-swX(2):nX(2)+swX(2), &
     !       1-swX(3):nX(3)+swX(3), &
     !       1:nCF) )

    ! ************************************************************************************
    ! No interpolation
    ! ************************************************************************************
    if (interp_type .eq. 0) then
       do kc = lo(3)-ng,hi(3)+ng
       do jc = lo(2)-ng,hi(2)+ng
       do ic = lo(1)-ng,hi(1)+ng

           !   S spatial indices start at lo - (number of ghost zones)
           ! uCF spatial indices start at 1 - (number of ghost zones)
           i = ic - lo(1) + 1
           j = jc - lo(2) + 1
           k = kc - lo(3) + 1

          ! Thornado uses units where c = G = k = 1, Meter = 1
          uCF(1:n_fluid_dof,i,j,k,iCF_D)  = S(ic,jc,kc,URHO)  * conv_dens
          uCF(1:n_fluid_dof,i,j,k,iCF_S1) = S(ic,jc,kc,UMX)   * conv_mom
          uCF(1:n_fluid_dof,i,j,k,iCF_S2) = S(ic,jc,kc,UMY)   * conv_mom
          uCF(1:n_fluid_dof,i,j,k,iCF_S3) = S(ic,jc,kc,UMZ)   * conv_mom
          uCF(1:n_fluid_dof,i,j,k,iCF_E)  = S(ic,jc,kc,UEDEN) * conv_enr
          uCF(1:n_fluid_dof,i,j,k,iCF_Ne) = S(ic,jc,kc,UFX)   * conv_ne

       end do
       end do
       end do

    ! ************************************************************************************
    ! Interpolate from the Castro "S" arrays into Thornado "uCF" arrays 
    !    using the conserved variables
    ! ************************************************************************************
    else if (interp_type .eq. 1) then

       do kc = lo(3)-ng,hi(3)+ng
       do jc = lo(2)-ng,hi(2)+ng
       do ic = lo(1)-ng,hi(1)+ng

         ! Define limited second-order slopes in x-direction
         dlft(:) = (S(ic  ,jc,kc,:) - S(ic-1,jc,kc,:)) * two
         drgt(:) = (S(ic+1,jc,kc,:) - S(ic  ,jc,kc,:)) * two
         dcen(:) = (S(ic+1,jc,kc,:) - S(ic-1,jc,kc,:)) * half

         do n = 1, ns
            dsgn = sign(one, dcen(n))
            xslope(n) = min( abs(dlft(n)), abs(drgt(n)) )
            if (dlft(n) * drgt(n) .ge. zero) then
               dlim = xslope(n)
            else
               dlim = zero
            endif
            xslope(n) = dsgn * min( dlim, abs(dcen(n)) )
         end do

         ! Define limited second-order slopes in y-direction
         dlft(:) = (S(ic,jc  ,kc,:) - S(ic,jc-1,kc,:)) * two
         drgt(:) = (S(ic,jc+1,kc,:) - S(ic,jc  ,kc,:)) * two
         dcen(:) = (S(ic,jc+1,kc,:) - S(ic,jc-1,kc,:)) * half

         do n = 1, ns
            dsgn = sign(one, dcen(n))
            yslope(n) = min( abs(dlft(n)), abs(drgt(n)) )
            if (dlft(n) * drgt(n) .ge. zero) then
               dlim = yslope(n)
            else
               dlim = zero
            endif
            yslope(n) = dsgn * min( dlim, abs(dcen(n)) )
         end do

         ! Define limited second-order slopes in z-direction
         dlft(:) = (S(ic,jc,kc  ,:) - S(ic,jc,kc-1,:)) * two
         drgt(:) = (S(ic,jc,kc+1,:) - S(ic,jc,kc  ,:)) * two
         dcen(:) = (S(ic,jc,kc+1,:) - S(ic,jc,kc-1,:)) * half

         do n = 1, ns
            dsgn = sign(one, dcen(n))
            zslope(n) = min( abs(dlft(n)), abs(drgt(n)) )
            if (dlft(n) * drgt(n) .ge. zero) then
               dlim = zslope(n)
            else
               dlim = zero
            endif
            zslope(n) = dsgn * min( dlim, abs(dcen(n)) )
         end do

         fac(:) = 1.d0
  
         ! Limit the slopes so that we can't overshoot corner values
         do n = 1, ns
            if (n.eq.URHO .or. n.eq.UEDEN .or. n.eq.UFX) then

               ! *********************************************
               ! Limit relative to (hi,hi,hi) 
               ! *********************************************
               delta_S_slope(n) = half*xslope(n) + half*yslope(n) + half*zslope(n)
               delta_S_val(n)   = S(ic+1,jc+1,kc+1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (lo,hi,hi) 
               ! *********************************************
               delta_S_slope(n) = -half*xslope(n) + half*yslope(n) + half*zslope(n)
               delta_S_val(n)   = S(ic-1,jc+1,kc+1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (hi,lo,hi) 
               ! *********************************************
               delta_S_slope(n) = half*xslope(n) - half*yslope(n) + half*zslope(n)
               delta_S_val(n)   = S(ic+1,jc-1,kc+1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (lo,lo,hi) 
               ! *********************************************
               delta_S_slope(n) = -half*xslope(n) - half*yslope(n) + half*zslope(n)
               delta_S_val(n)   = S(ic-1,jc-1,kc+1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (hi,hi,lo) 
               ! *********************************************
               delta_S_slope(n) = half*xslope(n) + half*yslope(n) - half*zslope(n)
               delta_S_val(n)   = S(ic+1,jc+1,kc-1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (lo,hi,lo) 
               ! *********************************************
               delta_S_slope(n) = -half*xslope(n) + half*yslope(n) - half*zslope(n)
               delta_S_val(n)   = S(ic-1,jc+1,kc-1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (hi,lo,lo) 
               ! *********************************************
               delta_S_slope(n) = half*xslope(n) - half*yslope(n) - half*zslope(n)
               delta_S_val(n)   = S(ic+1,jc-1,kc-1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (lo,lo,lo) 
               ! *********************************************
               delta_S_slope(n) = -half*xslope(n) - half*yslope(n) - half*zslope(n)
               delta_S_val(n)   = S(ic-1,jc-1,kc-1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if
 
               if (fac(n).lt.0.d0 .or. fac(n) .gt. 1.d0) then
                  print *,'WRONG ',ic,jc,kc,n,fac(n)
                  stop
               end if

            end if ! if n = URHO or UEDEN or UFX
         end do ! n = 1,ns

         fac_all = fac(1) 
         do n = 2,ns
            fac_all = min(fac_all,fac(n))
         end do

         xslope(:) = fac_all * xslope(:)
         yslope(:) = fac_all * yslope(:)
         zslope(:) = fac_all * zslope(:)

         !   S spatial indices start at lo - (number of ghost zones)
         ! uCF spatial indices start at 1 - (number of ghost zones)
         i = ic - lo(1) + 1
         j = jc - lo(2) + 1
         k = kc - lo(3) + 1

         do ind = 1, n_fluid_dof

            Sval(:) = S(ic,jc,kc,:) + x(ind)*xslope(:) + y(ind)*yslope(:) + z(ind)*zslope(:)

            ! Sanity check to make sure interpolation keeps us in range EOS is happy with
            rho_in(1) = Sval(URHO )  * conv_dens
            E_in  (1) = Sval(UEDEN) / Sval(URHO) * (conv_enr / conv_dens)
            Ye_in (1) = Sval(UFX)  / Sval(URHO)  * (AtomicMassUnit/Gram)

            call ComputeTemperatureFromSpecificInternalEnergy_TABLE( rho_in, E_in, Ye_in, T_out)

            ! Thornado uses units where c = G = k = 1, Meter = 1
            uCF(ind,i,j,k,iCF_D)  = Sval(URHO)  * conv_dens
            uCF(ind,i,j,k,iCF_S1) = Sval(UMX)   * conv_mom
            uCF(ind,i,j,k,iCF_S2) = Sval(UMY)   * conv_mom
            uCF(ind,i,j,k,iCF_S3) = Sval(UMZ)   * conv_mom
            uCF(ind,i,j,k,iCF_E)  = Sval(UEDEN) * conv_enr
            uCF(ind,i,j,k,iCF_Ne) = Sval(UFX)   * conv_ne

         end do ! do ind

    end do
    end do
    end do

    ! ************************************************************************************
    ! Interpolate from the Castro "S" arrays into Thornado "uCF" arrays 
    !    using the primitive variables (rho, E, Ye) if (interp_type .eq. 2)
    !    using the primitive variables (rho, T, Ye) if (interp_type .eq. 3 or 4)
    ! ************************************************************************************
    else if (interp_type .eq. 2 .or. &
             interp_type .eq. 3 .or. &
             interp_type .eq. 4) then

       ! Convert (rho E) to E and (rho Y) to Y for interpolation
       do kc = lo(3)-ng-1,hi(3)+ng+1
       do jc = lo(2)-ng-1,hi(2)+ng+1
       do ic = lo(1)-ng-1,hi(1)+ng+1
          S(ic,jc,kc,UEDEN) = S(ic,jc,kc,UEDEN) / S(ic,jc,kc,URHO)
          S(ic,jc,kc,UFX  ) = S(ic,jc,kc,UFX  ) / S(ic,jc,kc,URHO)
       end do
       end do
       end do

       ! Create temperature of the original (un-interpolated) data
       if (interp_type .eq. 3) then
          do kc = lo(3)-ng-1,hi(3)+ng+1
          do jc = lo(2)-ng-1,hi(2)+ng+1
          do ic = lo(1)-ng-1,hi(1)+ng+1

             rho_in(1) = S(ic,jc,kc,URHO ) * conv_dens
             E_in  (1) = S(ic,jc,kc,UEDEN) * (conv_enr / conv_dens)
             Ye_in (1) = S(ic,jc,kc,UFX)   * (AtomicMassUnit/Gram)
 
             call ComputeTemperatureFromSpecificInternalEnergy_TABLE( rho_in, E_in, Ye_in, T_out)

             S(ic,jc,kc,UTEMP) = T_out(1) / Kelvin

          end do
          end do
          end do
       end if

       if (interp_type.eq.2 .or. interp_type.eq.3) then

          do kc = lo(3)-ng,hi(3)+ng
          do jc = lo(2)-ng,hi(2)+ng
          do ic = lo(1)-ng,hi(1)+ng

            ! Define limited second-order slopes in x-direction
            dlft(:) = (S(ic  ,jc,kc,:) - S(ic-1,jc,kc,:)) * two
            drgt(:) = (S(ic+1,jc,kc,:) - S(ic  ,jc,kc,:)) * two
            dcen(:) = (S(ic+1,jc,kc,:) - S(ic-1,jc,kc,:)) * half

            do n = 1, ns
               dsgn = sign(one, dcen(n))
               xslope(n) = min( abs(dlft(n)), abs(drgt(n)) )
               if (dlft(n) * drgt(n) .ge. zero) then
                  dlim = xslope(n)
               else
                  dlim = zero
               endif
               xslope(n) = dsgn * min( dlim, abs(dcen(n)) )
            end do

            ! Define limited second-order slopes in y-direction
            dlft(:) = (S(ic,jc  ,kc,:) - S(ic,jc-1,kc,:)) * two
            drgt(:) = (S(ic,jc+1,kc,:) - S(ic,jc  ,kc,:)) * two
            dcen(:) = (S(ic,jc+1,kc,:) - S(ic,jc-1,kc,:)) * half
   
            do n = 1, ns
               dsgn = sign(one, dcen(n))
               yslope(n) = min( abs(dlft(n)), abs(drgt(n)) )
               if (dlft(n) * drgt(n) .ge. zero) then
                  dlim = yslope(n)
               else
                  dlim = zero
               endif
               yslope(n) = dsgn * min( dlim, abs(dcen(n)) )
            end do

            ! Define limited second-order slopes in z-direction
            dlft(:) = (S(ic,jc,kc  ,:) - S(ic,jc,kc-1,:)) * two
            drgt(:) = (S(ic,jc,kc+1,:) - S(ic,jc,kc  ,:)) * two
            dcen(:) = (S(ic,jc,kc+1,:) - S(ic,jc,kc-1,:)) * half

            do n = 1, ns
               dsgn = sign(one, dcen(n))
               zslope(n) = min( abs(dlft(n)), abs(drgt(n)) )
               if (dlft(n) * drgt(n) .ge. zero) then
                  dlim = zslope(n)
               else
                  dlim = zero
               endif
               zslope(n) = dsgn * min( dlim, abs(dcen(n)) )
            end do

            fac(:) = 1.d0
   
            ! Limit the slopes so that we can't overshoot corner values
            do n = 1, ns

            if (n.eq.URHO .or. n.eq.UFX .or. &
                  (interp_type.eq.2 .and. n.eq.UEDEN) .or. &
                  (interp_type.eq.3 .and. n.eq.UTEMP) ) then

               ! *********************************************
               ! Limit relative to (hi,hi,hi) 
               ! *********************************************
               delta_S_slope(n) = half*xslope(n) + half*yslope(n) + half*zslope(n)
               delta_S_val(n)   = S(ic+1,jc+1,kc+1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (lo,hi,hi) 
               ! *********************************************
               delta_S_slope(n) = -half*xslope(n) + half*yslope(n) + half*zslope(n)
               delta_S_val(n)   = S(ic-1,jc+1,kc+1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (hi,lo,hi) 
               ! *********************************************
               delta_S_slope(n) = half*xslope(n) - half*yslope(n) + half*zslope(n)
               delta_S_val(n)   = S(ic+1,jc-1,kc+1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (lo,lo,hi) 
               ! *********************************************
               delta_S_slope(n) = -half*xslope(n) - half*yslope(n) + half*zslope(n)
               delta_S_val(n)   = S(ic-1,jc-1,kc+1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (hi,hi,lo) 
               ! *********************************************
               delta_S_slope(n) = half*xslope(n) + half*yslope(n) - half*zslope(n)
               delta_S_val(n)   = S(ic+1,jc+1,kc-1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (lo,hi,lo) 
               ! *********************************************
               delta_S_slope(n) = -half*xslope(n) + half*yslope(n) - half*zslope(n)
               delta_S_val(n)   = S(ic-1,jc+1,kc-1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (hi,lo,lo) 
               ! *********************************************
               delta_S_slope(n) = half*xslope(n) - half*yslope(n) - half*zslope(n)
               delta_S_val(n)   = S(ic+1,jc-1,kc-1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if

               ! *********************************************
               ! Limit relative to (lo,lo,lo) 
               ! *********************************************
               delta_S_slope(n) = -half*xslope(n) - half*yslope(n) - half*zslope(n)
               delta_S_val(n)   = S(ic-1,jc-1,kc-1,n) - S(ic,jc,kc,n)

               ! If opposite signs then zero the slopes
               if ( (delta_S_slope(n) *  delta_S_val(n)) .le. 0.d0 ) then
                  fac(n) = 0.d0

               ! Else make sure that the corner value extrapolated  using slopes
               !      is no further from the value at (ic,jc) than the actual corner value
               else if (delta_S_slope(n) .gt. delta_S_val(n)) then
                  fac(n) = min(fac(n), delta_S_val(n)/delta_S_slope(n))
               end if
 
               if (fac(n).lt.0.d0 .or. fac(n) .gt. 1.d0) then
                  print *,'WRONG ',ic,jc,n,fac(n)
                  stop
               end if

            end if ! if n = URHO or UEDEN/UTEMP or UFX
            end do ! n = 1,ns

            fac_all = fac(1) 
            do n = 2,ns
               fac_all = min(fac_all,fac(n))
            end do

            xslope(:) = fac_all * xslope(:)
            yslope(:) = fac_all * yslope(:)
            zslope(:) = fac_all * zslope(:)

            !   S spatial indices start at lo - (number of ghost zones)
            ! uCF spatial indices start at 1 - (number of ghost zones)
            i = ic - lo(1) + 1
            j = jc - lo(2) + 1
            k = kc - lo(3) + 1

            do ind = 1, n_fluid_dof

               Sval(:) = S(ic,jc,kc,:) + x(ind)*xslope(:) + y(ind)*yslope(:) + z(ind)*zslope(:)

               rho_in(1) = Sval(URHO ) * conv_dens
               Ye_in (1) = Sval(UFX)   * (AtomicMassUnit/Gram)

               ! Sanity check to make sure interpolation keeps us in range EOS is happy with
               if (interp_type .eq. 2) then
   
                  E_in  (1) = Sval(UEDEN) * (conv_enr / conv_dens)
                  call ComputeTemperatureFromSpecificInternalEnergy_TABLE( rho_in, E_in, Ye_in, T_out)
                  Sval(UEDEN) = Sval(UEDEN) * Sval(URHO)

               ! Create E from interpolated T
               else if (interp_type .eq. 3) then
   
                  T_in  (1) = Sval(UTEMP) * Kelvin
                  call ComputeThermodynamicStates_Primitive_TABLE(rho_in, T_in, Ye_in, &
                                                                  Epervol_out, Epermass_out, Ne_out )
                  Sval(UEDEN) = Epervol_out(1) / conv_enr
               end if

               ! Make sure to pass the conserved variables to thornado
               Sval(UFX  ) = Sval(UFX  ) * Sval(URHO)

               ! Thornado uses units where c = G = k = 1, Meter = 1
               uCF(ind,i,j,k,iCF_D)  = Sval(URHO)  * conv_dens
               uCF(ind,i,j,k,iCF_S1) = Sval(UMX)   * conv_mom
               uCF(ind,i,j,k,iCF_S2) = Sval(UMY)   * conv_mom
               uCF(ind,i,j,k,iCF_S3) = Sval(UMZ)   * conv_mom
               uCF(ind,i,j,k,iCF_E)  = Sval(UEDEN) * conv_enr
               uCF(ind,i,j,k,iCF_Ne) = Sval(UFX)   * conv_ne

            end do ! do ind

          end do
          end do
          end do

       else

          if (interp_type.ne.4) then
            print *,"IF I GOT HERE THEN -- OOPS!"
            stop
          end if

          do kc = lo(3)-ng,hi(3)+ng
          do jc = lo(2)-ng,hi(2)+ng
          do ic = lo(1)-ng,hi(1)+ng

            ! Define values at all four corners of the cell

            S_hhh(:) = 0.125d0 * ( S(ic,jc  ,kc  ,:) + S(ic+1,jc  ,kc  ,:) &
                                  +S(ic,jc+1,kc  ,:) + S(ic+1,jc+1,kc  ,:) &
                                  +S(ic,jc  ,kc+1,:) + S(ic+1,jc  ,kc+1,:) &
                                  +S(ic,jc+1,kc+1,:) + S(ic+1,jc+1,kc+1,:) )
            S_lhh(:) = 0.125d0 * ( S(ic,jc  ,kc  ,:) + S(ic-1,jc  ,kc  ,:) &
                                  +S(ic,jc+1,kc  ,:) + S(ic-1,jc+1,kc  ,:) &
                                  +S(ic,jc  ,kc+1,:) + S(ic-1,jc  ,kc+1,:) &
                                  +S(ic,jc+1,kc+1,:) + S(ic-1,jc+1,kc+1,:) )
            S_hlh(:) = 0.125d0 * ( S(ic,jc  ,kc  ,:) + S(ic+1,jc  ,kc  ,:) &
                                  +S(ic,jc-1,kc  ,:) + S(ic+1,jc-1,kc  ,:) &
                                  +S(ic,jc  ,kc+1,:) + S(ic+1,jc  ,kc+1,:) &
                                  +S(ic,jc-1,kc+1,:) + S(ic+1,jc-1,kc+1,:) )
            S_llh(:) = 0.125d0 * ( S(ic,jc  ,kc  ,:) + S(ic-1,jc  ,kc  ,:) &
                                  +S(ic,jc-1,kc  ,:) + S(ic-1,jc-1,kc  ,:) &
                                  +S(ic,jc  ,kc+1,:) + S(ic-1,jc  ,kc+1,:) &
                                  +S(ic,jc-1,kc+1,:) + S(ic-1,jc-1,kc+1,:) )
            S_hhl(:) = 0.125d0 * ( S(ic,jc  ,kc  ,:) + S(ic+1,jc  ,kc  ,:) &
                                  +S(ic,jc+1,kc  ,:) + S(ic+1,jc+1,kc  ,:) &
                                  +S(ic,jc  ,kc+1,:) + S(ic+1,jc  ,kc-1,:) &
                                  +S(ic,jc+1,kc+1,:) + S(ic+1,jc+1,kc-1,:) )
            S_lhl(:) = 0.125d0 * ( S(ic,jc  ,kc  ,:) + S(ic-1,jc  ,kc  ,:) &
                                  +S(ic,jc+1,kc  ,:) + S(ic-1,jc+1,kc  ,:) &
                                  +S(ic,jc  ,kc+1,:) + S(ic-1,jc  ,kc-1,:) &
                                  +S(ic,jc+1,kc+1,:) + S(ic-1,jc+1,kc-1,:) )
            S_hll(:) = 0.125d0 * ( S(ic,jc  ,kc  ,:) + S(ic+1,jc  ,kc  ,:) &
                                  +S(ic,jc-1,kc  ,:) + S(ic+1,jc-1,kc  ,:) &
                                  +S(ic,jc  ,kc+1,:) + S(ic+1,jc  ,kc-1,:) &
                                  +S(ic,jc-1,kc+1,:) + S(ic+1,jc-1,kc-1,:) )
            S_lll(:) = 0.125d0 * ( S(ic,jc  ,kc  ,:) + S(ic-1,jc  ,kc  ,:) &
                                  +S(ic,jc-1,kc  ,:) + S(ic-1,jc-1,kc  ,:) &
                                  +S(ic,jc  ,kc+1,:) + S(ic-1,jc  ,kc-1,:) &
                                  +S(ic,jc-1,kc+1,:) + S(ic-1,jc-1,kc-1,:) )

             xslope(:) = 0.25d0 * (S_hhh(:) + S_hlh(:) - S_lhh(:) - S_llh(:) + &
                                   S_hhl(:) + S_hll(:) - S_lhl(:) - S_lll(:))
             yslope(:) = 0.25d0 * (S_hhh(:) + S_lhh(:) - S_hlh(:) - S_llh(:) + &
                                   S_hhl(:) + S_lhl(:) - S_hll(:) - S_lll(:))
             zslope(:) = 0.25d0 * (S_hhh(:) + S_lhh(:) + S_hlh(:) + S_llh(:) - &
                                   S_hhl(:) - S_lhl(:) - S_hll(:) - S_llh(:))

            xyslope(:) = 0.5d0  * (S_hhh(:) + S_llh(:) - S_hlh(:) - S_lhh(:) + &
                                   S_hhl(:) + S_lll(:) - S_hll(:) - S_lhl(:) )
            xzslope(:) = 0.5d0  * (S_hhh(:) + S_lhl(:) - S_hhl(:) - S_lhh(:) + &
                                   S_hlh(:) + S_lll(:) - S_hll(:) - S_llh(:) )
            yzslope(:) = 0.5d0  * (S_hhh(:) + S_hll(:) - S_hhl(:) - S_hlh(:) + &
                                   S_lhh(:) + S_lll(:) - S_lhl(:) - S_llh(:) )

            !   S spatial indices start at lo - (number of ghost zones)
            ! uCF spatial indices start at 1 - (number of ghost zones)
            i = ic - lo(1) + 1
            j = jc - lo(2) + 1
            k = kc - lo(3) + 1
   
            do ind = 1, n_fluid_dof

               Sval(:) = S(ic,jc,kc,:) + x(ind)*xslope(:) + y(ind)*yslope(:) + z(ind)*zslope(:) + &
                                         x(ind)*y(ind)*xyslope(:) + &
                                         x(ind)*z(ind)*xzslope(:) + &
                                         y(ind)*z(ind)*yzslope(:) 
   
               rho_in(1) = Sval(URHO ) * conv_dens
               Ye_in (1) = Sval(UFX)   * (AtomicMassUnit/Gram)
   
               ! Create E from interpolated T
   
               T_in  (1) = Sval(UTEMP) * Kelvin
               call ComputeThermodynamicStates_Primitive_TABLE(rho_in, T_in, Ye_in, &
                                                               Epervol_out, Epermass_out, Ne_out )
               Sval(UEDEN) = Epervol_out(1) / conv_enr

               ! Make sure to pass the conserved variables to thornado
               Sval(UFX  ) = Sval(UFX  ) * Sval(URHO)
   
               ! Thornado uses units where c = G = k = 1, Meter = 1
               uCF(ind,i,j,k,iCF_D)  = Sval(URHO)  * conv_dens
               uCF(ind,i,j,k,iCF_S1) = Sval(UMX)   * conv_mom
               uCF(ind,i,j,k,iCF_S2) = Sval(UMY)   * conv_mom
               uCF(ind,i,j,k,iCF_S3) = Sval(UMZ)   * conv_mom
               uCF(ind,i,j,k,iCF_E)  = Sval(UEDEN) * conv_enr
               uCF(ind,i,j,k,iCF_Ne) = Sval(UFX)   * conv_ne

            end do ! do ind

          end do
          end do
          end do

       end if ! interp_type = 4

       ! Convert E back to (rho E) and Y back to (rho Y) 
       do kc = lo(3)-ng-1,hi(3)+ng+1
       do jc = lo(2)-ng-1,hi(2)+ng+1
       do ic = lo(1)-ng-1,hi(1)+ng+1
          S(ic,jc,kc,UEDEN) = S(ic,jc,kc,UEDEN) * S(ic,jc,kc,URHO)
          S(ic,jc,kc,UFX  ) = S(ic,jc,kc,UFX  ) * S(ic,jc,kc,URHO)
       end do
       end do
       end do

    else
      call  amrex_abort("Dont know this interp_type in interpolate fluid!")
    end if

  end subroutine interpolate_fluid
