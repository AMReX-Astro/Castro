module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use prob_params_module, only: center
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ , UEDEN, UEINT, UFS, UTEMP, const_grav
    use interpolate_module
    use eos_module
    use eos_type_module
    use actual_eos_module, only : gamma_const
    use network, only: nspec

    use amrex_fort_module     , only : rt => amrex_real
    use amrex_constants_module, only : M_PI, sixth
    implicit none
    
    include 'AMReX_bc_types.fi'

    real(rt), parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
    real(rt), parameter :: pi_over_3 = M_PI / 3.d0
    real(rt), parameter :: ff = 0.25d0
    
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer i,j,ii,jj,n
    real(rt)         x,y,xcen,ycen
    real(rt)         X_in(nspec)
    real(rt)         shockfront 

    integer npts_1d
    real(rt)         const

    type (eos_t) :: eos_state

    do n=1,NVAR
       call filcc(adv(adv_l1,adv_l2,n),adv_l1,adv_l2,adv_h1,adv_h2, &
            domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    !        XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then

       do j=adv_l2,adv_h2
          do i=domlo(1)-1,adv_l1,-1
                ! Need to stay as initial post shock conditions.
                adv(i,j,URHO) = rho_l
                adv(i,j,UMX) = rho_l*u_l
                adv(i,j,UMY) = rho_l*v_l
                adv(i,j,UMZ) = 0.0e0_rt 
                adv(i,j,UEDEN) = rhoe_l + 0.5e0_rt*rho_l*(u_l*u_l + v_l*v_l)
                adv(i,j,UEINT) = rhoe_l
                adv(i,j,UTEMP) = T_l
                adv(i,j,UFS:UFS-1+nspec) = 0.0e0_rt 
                adv(i,j,UFS) = adv(i,j,URHO)
          end do
       end do

    end if

    !        XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then

       do j=adv_l2,adv_h2
          do i=domhi(1)+1,adv_h1
                ! extrapolate all vars for outflow on left boundary
                adv(i,j,:) = adv(domhi(1),j,:)
          end do
       end do

    end if

    !        YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       ! this do loop counts backwards since we want to work downward
       do j=domlo(2)-1,adv_l2,-1
          do i=adv_l1,adv_h1
             x = xlo(1) + delta(1)*(float(i-adv_l1) + 0.5e0_rt)
                if(x.lt.1.e0_rt/6.e0_rt) then 
                        ! ICs
                        adv(i,j,URHO) = rho_l
                        adv(i,j,UMX) = rho_l*u_l
                        adv(i,j,UMY) = rho_l*v_l 
                        adv(i,j,UEDEN) = rhoe_l + 0.5e0_rt*rho_l*(u_l*u_l + v_l*v_l)
                        adv(i,j,UEINT) = rhoe_l
                        adv(i,j,UTEMP) = T_l
                else 
                        ! reflect!
                        adv(i,j,URHO) =  adv(i,domlo(2),URHO)
                        adv(i,j,UMX)  =  adv(i,domlo(2),UMX)
                        adv(i,j,UMY)  =  -1.e0_rt*adv(i,domlo(2),UMY)
                        adv(i,j,UEDEN) = adv(i,domlo(2),UEDEN)
                        adv(i,j,UEINT) = adv(i,domlo(2),UEINT)
                        adv(i,j,UTEMP) = adv(i,domlo(2),UTEMP)
               endif
                adv(i,j,UMZ) = 0.0e0_rt 
                adv(i,j,UFS:UFS-1+nspec) = 0.0e0_rt 
                adv(i,j,UFS) = adv(i,j,URHO)
        end do
       end do
    end if

    !        YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then

       do j=domhi(2)+1,adv_h2
          ycen = delta(2)*(dble(j)+0.5d0)

       do i=adv_l1,adv_h1
          xcen = delta(1)*(dble(i) + 0.5d0)

          adv(i,j,URHO ) = 0.d0
          adv(i,j,UMX  ) = 0.d0
          adv(i,j,UMY  ) = 0.d0
          adv(i,j,UEDEN) = 0.d0
          adv(i,j,UEINT) = 0.d0
          adv(i,j,UTEMP) = 0.d0

        do jj = 1,2
          y = ycen + 0.5d0*delta(2)*gp(jj)

          shockfront = sixth + y/tan(pi_over_3) + (10.d0/sin(pi_over_3))*time
!         shockfront = sixth + (10.d0)*time

          do ii = 1,2

             x = xcen + 0.5d0*delta(1)*gp(ii)

             if (x.lt.shockfront) then 
             ! Post shock ICs
                     adv(i,j,URHO ) = adv(i,j,URHO ) + ff*rho_l
                     adv(i,j,UMX  ) = adv(i,j,UMX  ) + ff*rho_l*u_l
                     adv(i,j,UMY  ) = adv(i,j,UMY  ) + ff*rho_l*v_l
                     adv(i,j,UEDEN) = adv(i,j,UEDEN) + ff*(rhoe_l + 0.5e0_rt*rho_l*(u_l*u_l + v_l*v_l))
                     adv(i,j,UEINT) = adv(i,j,UEINT) + ff*rhoe_l
                     adv(i,j,UTEMP) = adv(i,j,UTEMP) + ff*T_l
              else
             !Pre Shock ICs
                     adv(i,j,URHO ) = adv(i,j,URHO ) + ff*rho_r
                     adv(i,j,UMX  ) = adv(i,j,UMX  ) + ff*rho_r*u_r
                     adv(i,j,UMY  ) = adv(i,j,UMY  ) + ff*rho_r*v_r 
                     adv(i,j,UEDEN) = adv(i,j,UEDEN) + ff*(rhoe_r + 0.5e0_rt*rho_r*(u_r*u_r + v_r*v_r))
                     adv(i,j,UEINT) = adv(i,j,UEINT) + ff*rhoe_r
                     adv(i,j,UTEMP) = adv(i,j,UTEMP) + ff*T_r
             endif

           end do
           end do

           adv(i,j,UFS:UFS-1+nspec) = 0.0e0_rt 
           adv(i,j,UFS) = adv(i,j,URHO)
           adv(i,j,UMZ  ) = 0.0e0_rt 

         end do
       end do

    end if

  end subroutine ca_hypfill

  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use interpolate_module
    use actual_eos_module, only: gamma_const
    use meth_params_module, only : const_grav

    use amrex_constants_module, only : M_PI, sixth
    use amrex_fort_module     , only : rt => amrex_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    real(rt), parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
    real(rt), parameter :: pi_over_3 = M_PI / 3.d0
    real(rt), parameter :: ff = 0.25d0
    
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2)

    integer i,j,ii,jj
    real(rt)       x,y,xcen,ycen,shockfront

    !     Note: this function should not be needed, technically, but is provided
    !     to filpatch because there are many times in the algorithm when just
    !     the density is needed.  We try to rig up the filling so that the same
    !     function is called here and in hypfill where all the states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       do j=adv_l2,adv_h2
          do i=adv_l1,domlo(1)-1
                ! Need to stay as initial post shock conditions.
                adv(i,j) = rho_l
          end do
       end do
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       do j=adv_l2,adv_h2
          do i=domhi(1)+1,adv_h1
               adv(i,j) = adv(domhi(1), j)
          end do
       end do
    end if

    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       do j=domlo(2)-1,adv_l2,-1
          do i=adv_l1,adv_h1
             x = xlo(1) + delta(1)*(float(i-adv_l1)+0.5e0_rt)
                if(x.le.1.0e0_rt/6.0e0_rt) then 
                        ! IC
                        adv(i,j) = rho_l
                else 
                        ! reflect!
                        adv(i,j) =  adv(i,domlo(2)-j)
                endif
          end do
       end do
    end if

    !     YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then

       do j=domhi(2)+1,adv_h2
          ycen = delta(2)*(dble(j)+0.5d0)

         do i=adv_l1,adv_h1
         xcen = delta(1)*(dble(i) + 0.5d0)

            adv(i,j) = 0.d0

            do jj = 1, 2
              y = ycen + 0.5d0*delta(2)*gp(jj)

              shockfront = sixth + y/tan(pi_over_3) + (10.d0/sin(pi_over_3))*time
!             shockfront = sixth + (10.d0)*time

              do ii = 1, 2
              x = xcen + 0.5d0*delta(1)*gp(ii)

                if (x.lt.shockfront) then 
                 ! Post shock ICs
                     adv(i,j) = adv(i,j) + ff*rho_l
                  else
                 !Pre Shock ICs
                     adv(i,j) = adv(i,j) + ff*rho_r
                  endif
              end do
            end do
         end do
       end do
    end if

  end subroutine ca_denfill

end module bc_fill_module
