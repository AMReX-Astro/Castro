module bc_fill_module

  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C)
 
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, &
                                   UFS, UTEMP, const_grav
    use interpolate_module
    use eos_module
    use network, only: nspec
    use model_parser_module
    use bl_error_module
    use eos_type_module

    include 'AMReX_bc_types.fi'

    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    double precision delta(2), xlo(2), time
    double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)
    
    integer i,j,q,n
    double precision y
    double precision pres_above,p_want,pres_zone
    double precision temp_zone,X_zone(nspec),dens_zone
    double precision :: y_base, dens_base, slope
  
    type (eos_t) :: eos_state
!print *, 'const_grav=',const_grav!   const_grav=-1.0d3
!print *,'bc:AMReX_FILCC_2D.f, ghost cell upper right(x,y)',adv_h1,adv_h2
!print *,'bc:AMReX_FILCC_2D.f, ghost cell bottom left(y)',adv_l2
!print *,'bc:AMReX_FILCC_2D.f, ghost cell bottom left(x,y)',adv_l1,adv_l2
!print *,'bc:AMReX_FILCC_2D.f, ghost cell upper right(y)',adv_h2
!print *,'bc:AMReX_FILCC_2D.f, normalcell upper right(x,y)',domhi(1),domhi(2)

!print *,'bc:AMReX_FILCC_2D.f, normalcell bottom left(x,y)',domlo(1),domlo(2)
!print *, 'bc:AMRex_FICC_2D.f',xlo(1),xlo(2),NVAR,EXT_DIR
!print *,'BCs(lower):',bc(1,1,1),bc(2,1,1)
!print *,'BCs(upper):',bc(1,2,1),bc(2,2,1)
!          do j=adv_h2,adv_l2,-1
!          if(adv_l1<=adv_h1) then
!          if(j==adv_h2)then
!         print*, 'BEFORE----------------------------------------------------------------'
! write(*,"(A10,f6.3,A10,I3,A10,I3,A10,I3)") 'time1:',time,'ghost',j, 'ghost',adv_l1, 'normal',adv_h1
!         write(*,"(A10,f6.3,A10,I3,10000f15.4)"), 'time2:',time,'ghost',j,(adv(i,j,UTEMP),i=adv_l1,adv_h1)
!          else if(j==adv_l2)then
!         write(*,"(A10,f6.3,A10,I3,10000f15.4)"), 'time3:',time,'ghost',j,(adv(i,j,UTEMP),i=adv_l1,adv_h1)
!         print*, 'BEFORE----------------------------------------------------------------'
!          else
!         write(*,"(A10,f6.3,A10,I3,10000f15.4)"), 'time4:',time,'normal',j,(adv(i,j,UTEMP),i=adv_l1,adv_h1)
!       end if   
!       else
!          if(j==adv_h2)then
!         print*, '----------------------------------------------------------------'
! print *, 'time1:',time,'y',xlo(2),'ghost',j, 'ghost',adv_l1, '         normal',domlo(1)
!  print *, 'time2:',time,'y',xlo(2),'ghost',j,(adv(i,j,UTEMP),i=domlo(1),adv_l1)
!    else if(j==adv_l2)then
!         print *, 'time3:',time,'y',xlo(2),'ghost',j,(adv(i,j,UTEMP),i=domlo(1),adv_l1)
!         print*, '----------------------------------------------------------------'
!          else
!          print *, 'time4:',time,'y',xlo(2),'normal',j,(adv(i,j,UTEMP),i=domlo(1),adv_l1)
!       end if


!          end if
!       end do

  


    do n = 1,NVAR
       call filcc(adv(adv_l1,adv_l2,n),adv_l1,adv_l2,adv_h1,adv_h2, &
                  domlo,domhi,delta,xlo,bc(1,1,n))
!print *, 'bc3:AMRex_FICC_2D.f',n,bc(1,1,n),bc(1,2,n),bc(2,1,n),bc(2,2,n)
    enddo

!            do j=domhi(2)+1,adv_h2  
!            adv(adv_l1:adv_h1,j,UTEMP)=0.0d0
!          end do

    do n = 1, NVAR
         
       ! XLO
       if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then

          ! we are periodic in x -- we should never get here
          call bl_error("ERROR: invalid BC in Prob_2d.f90")
        
       end if
     
       ! XHI
       if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
          
          ! we are periodic in x -- we should never get here
          call bl_error("ERROR: invalid BC in Prob_2d.f90")
        
       end if
    enddo

    ! YLO -- HSE with linear density profile, T found via iteration
    ! we do all variables at once here
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then

       y_base = xlo(2) + delta(2)*(float(domlo(2)-adv_l2) + 0.5d0)
!print *, '1:',adv_l2,URHO        
       do i=adv_l1,adv_h1

          dens_base = adv(i,domlo(2),URHO)

          ! density slope
          slope = (adv(i,domlo(2)+1,URHO) - adv(i,domlo(2),URHO))/delta(2)
!print*,'2:',i,adv(i,domlo(2)+1,URHO),adv(i,domlo(2),URHO),slope           
          ! this do loop counts backwards since we want to work downward
          do j=domlo(2)-1,adv_l2,-1
             y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5d0)
!print*,'3:',j,y,xlo(2),j-adv_l2           
             ! zero-gradient catch-all -- this will get the radiation
             ! energy
             adv(i,j,:) = adv(i,j+1,:)

             ! HSE integration to get temperature, pressure
                    
             ! density is linear from the last two zones
             dens_zone = dens_base + slope*(y - y_base)

             ! temperature guess and species held constant in BCs
             temp_zone = adv(i,j+1,UTEMP)
!if(i==adv_l1 .and. j==adv_l2)then
!print *, 'temp_zone=',temp_zone
!print *, 'dens_zone=',dens_base,dens_zone,slope,y/1D5,y_base/1D5,adv(i,domlo(2)+1,URHO),adv(i,domlo(2),URHO)
!end if

             X_zone(:) = adv(i,j+1,UFS:UFS-1+nspec)/adv(i,j+1,URHO)
             
             ! get pressure in zone above
             eos_state%rho = adv(i,j+1,URHO)
             eos_state%T = adv(i,j+1,UTEMP)
             eos_state%xn(:) = adv(i,j+1,UFS:UFS-1+nspec)/adv(i,j+1,URHO)
             
             call eos(eos_input_rt, eos_state)
             
             pres_above = eos_state%p

             ! pressure needed from HSE
             p_want = pres_above - &
                  delta(2)*0.5d0*(dens_zone + adv(i,j+1,URHO))*const_grav

             ! EOS with HSE pressure + linear density profile yields T, e, ...
             eos_state%rho = dens_zone
             eos_state%T = temp_zone   ! guess
             eos_state%xn(:) = X_zone(:)
             eos_state%p = p_want
 
!if(i==adv_l1 .and. j==adv_l2)then
!print *, 'temp_zone=',temp_zone
!print *, 'dens_zone=',dens_base,dens_zone,slope,y/1D5,y_base/1D5,adv(i,domlo(2)+1,URHO),adv(i,domlo(2),URHO\
!)
!print *, 'pres_zone=',p_want,pres_above,delta(2),(dens_zone + adv(i,j+1,URHO)),const_grav
!end if


            
             call eos(eos_input_rp, eos_state)

             ! velocity
             if (zero_vels) then
                
                ! zero normal momentum causes pi waves to pass through
                adv(i,j,UMY) = 0.d0
                
                ! zero transverse momentum
                adv(i,j,UMX) = 0.d0
             else
                
                ! zero gradient velocity
                adv(i,j,UMX) = dens_zone*(adv(i,domlo(2),UMX)/adv(i,domlo(2),URHO))
                adv(i,j,UMY) = dens_zone*(adv(i,domlo(2),UMY)/adv(i,domlo(2),URHO))
             endif
           
             adv(i,j,URHO) = dens_zone
             adv(i,j,UEINT) = dens_zone*eos_state%e
             adv(i,j,UEDEN) = dens_zone*eos_state%e + & 
                  0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/dens_zone
             adv(i,j,UTEMP) = eos_state%T
             adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
             
          end do
       end do
    end if

  
  
    ! YHI
    do n = 1, nvar
       if ( bc(2,2,n).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
          
          do j=domhi(2)+1,adv_h2
             y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5d0)
!print *,'1_1:',j,y,xlo(2),j-adv_l2
             ! zero-gradient catch-all -- this will get the radiation
             ! energy
!print *,'2_1:',adv(adv_l1,j,UTEMP),adv(adv_l1,j-1,UTEMP)

             adv(adv_l1:adv_h1,j,:) = adv(adv_l1:adv_h1,j-1,:)
             adv(adv_l1:adv_h1,j,UTEMP)=min(adv(adv_l1:adv_h1,j-1,UTEMP),temp_zone)
             do i=adv_l1,adv_h1
                
                ! set all the variables even though we're testing on URHO
                if (n .eq. URHO) then
                     
                   dens_zone = interpolate(y,npts_model,model_r, &
                                           model_state(:,idens_model)) 

                   temp_zone = interpolate(y,npts_model,model_r, &
                                           model_state(:,itemp_model))
!if(j==domhi(2)+1 .and. i==domhi(1))then
!print *,n,'Upper,T  =', temp_zone,adv(domhi(1),domhi(2),UTEMP)
!print *,n,'Upper,rho=', dens_zone,adv(domhi(1),domhi(2),n)
!end if 
                   do q = 1, nspec
                      X_zone(q) = interpolate(y,npts_model,model_r, &
                                              model_state(:,ispec_model-1+q))
                   enddo

                   ! extrap normal momentum
                   adv(i,j,UMY) = max(0.d0,adv(i,domhi(2),UMY))
                   
                   ! zero transverse momentum
                   adv(i,j,UMX) = 0.d0
                   
                   eos_state%rho = dens_zone
                   eos_state%T = temp_zone
                   eos_state%xn(:) = X_zone

                   
                   call eos(eos_input_rt, eos_state)
                   
                   adv(i,j,URHO) = dens_zone
                   adv(i,j,UEINT) = dens_zone*eos_state%e
                   adv(i,j,UEDEN) = dens_zone*eos_state%e + &
                        0.5d0*(adv(i,j,UMX)**2+adv(i,j,UMY)**2)/dens_zone
                   adv(i,j,UTEMP) = temp_zone
                   adv(i,j,UFS:UFS-1+nspec) = dens_zone*X_zone(:)
               
                end if

!if(j==domhi(2)+1 .and. i==domhi(1))then
!print *,n,'upper,T  =', temp_zone,adv(domhi(1),domhi(2),UTEMP),adv(domhi(1),domhi(2)+1,UTEMP)
!print *,n,'upper,rho=', dens_zone,adv(domhi(1),domhi(2),URHO),adv(domhi(1),domhi(2)+1,URHO)
!end if
              
             end do
          end do
       end if
     
    end do

!         do j=adv_h2,adv_l2,-1
     !     if(adv_l1<=adv_h1) then                                                                          
!          if(j==adv_h2)then
!         print*, 'AFTER----------------------------------------------------------------'
! write(*,"(A10,f6.3,A10,I3,A10,I3,A10,I3)") 'time1:',time,'ghost',j, 'ghost',adv_l1, 'normal',adv_h1
!         write(*,"(A10,f6.3,A10,I3,10000f15.4)"), 'time2:',time,'ghost',j,(adv(i,j,UTEMP),i=adv_l1,adv_h1)
!          else if(j==adv_l2)then
!         write(*,"(A10,f6.3,A10,I3,10000f15.4)"), 'time3:',time,'ghost',j,(adv(i,j,UTEMP),i=adv_l1,adv_h1)
!         print*, 'AFTER----------------------------------------------------------------'
!          else
!         write(*,"(A10,f6.3,A10,I3,10000f15.4)"), 'time4:',time,'normal',j,(adv(i,j,UTEMP),i=adv_l1,adv_h1)
!       end if
!       end do

   ! do n = 1, nvar
   !    if ( bc(2,2,n).eq.FOEXTRAP .and. adv_h2.gt.domhi(2)) then
          
!          do j=domhi(2)+1,adv_h2
!           do j=adv_l2,adv_h2 
            !adv(adv_l1:adv_h1,j,UTEMP)=0.0d0
!            adv(adv_l1:adv_h1,j,URHO)=1D7
!          end do
   !    end if
!print *,(adv(adv_l1:adv_h1,j,UTEMP),j=domhi(2)+1,adv_h2)
   ! end do

  end subroutine ca_hypfill

  
  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C)
    
    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, &
         UFS, UTEMP, const_grav
    use bl_error_module
    use interpolate_module
    use model_parser_module

    implicit none
    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    double precision delta(2), xlo(2), time
    double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

    integer i,j,q,n
    double precision y
    double precision :: y_base, dens_base, slope
    double precision TOL

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       call bl_error("We shoundn't be here (xlo denfill)")
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       call bl_error("We shoundn't be here (xlo denfill)")
    endif


    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then

       y_base = xlo(2) + delta(2)*(float(domlo(2)-adv_l2) + 0.5d0)

       do i=adv_l1,adv_h1

          dens_base = adv(i,domlo(2))

          ! density slope
          slope = (adv(i,domlo(2)+1) - adv(i,domlo(2)))/delta(2)

          ! this do loop counts backwards since we want to work downward
          do j=domlo(2)-1,adv_l2,-1
             y = xlo(2) + delta(2)*(float(j-adv_l2) + 0.5d0)

             ! density is linear from the last two zones
             adv(i,j) = dens_base + slope*(y - y_base)

          end do
       end do
    end if

    !     YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
       do j=domhi(2)+1,adv_h2
          y = xlo(2) + delta(2)*(float(j-adv_l2)+ 0.5d0)
          do i=adv_l1,adv_h1
             adv(i,j) = interpolate(y,npts_model,model_r,model_state(:,idens_model))
          end do
       end do
    end if

  end subroutine ca_denfill

  
  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    double precision delta(2), xlo(2), time
    double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)
    integer :: i, j

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

    ! our lower boundary is inflow, so we need to make sure the
    ! gravitational acceleration is set correctly there
    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. grav_l2.lt.domlo(2)) then
       do j=grav_l2,domlo(2)-1
          do i=grav_l1,grav_h1
             grav(i,j) = 0.0
          end do
       end do
    end if

  end subroutine ca_gravxfill


  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use meth_params_module, only: const_grav

    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    double precision delta(2), xlo(2), time
    double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)
    integer :: i, j

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

    ! our lower boundary is inflow, so we need to make sure the
    ! gravitational acceleration is set correctly there
    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. grav_l2.lt.domlo(2)) then
       do j=grav_l2,domlo(2)-1
          do i=grav_l1,grav_h1
             grav(i,j) = const_grav
          end do
       end do
    end if

  end subroutine ca_gravyfill


  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use meth_params_module, only: const_grav

    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    double precision delta(2), xlo(2), time
    double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)
    integer :: i, j

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

    ! our lower boundary is inflow, so we need to make sure the
    ! gravitational acceleration is set correctly there
    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. grav_l2.lt.domlo(2)) then
       do j=grav_l2,domlo(2)-1
          do i=grav_l1,grav_h1
             grav(i,j) = 0.0
          end do
       end do
    end if

  end subroutine ca_gravzfill


  subroutine ca_reactfill(react,react_l1,react_l2, &
                          react_h1,react_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: react_l1,react_l2,react_h1,react_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    double precision delta(2), xlo(2), time
    double precision react(react_l1:react_h1,react_l2:react_h2)

    call filcc(react,react_l1,react_l2,react_h1,react_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_reactfill


  subroutine ca_radfill(rad,rad_l1,rad_l2, &
                        rad_h1,rad_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: rad_l1,rad_l2,rad_h1,rad_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    double precision delta(2), xlo(2), time
    double precision rad(rad_l1:rad_h1,rad_l2:rad_h2)

    integer :: j

    call filcc(rad,rad_l1,rad_l2,rad_h1,rad_h2,domlo,domhi,delta,xlo,bc)

    ! we are inflow at the lower boundary, so we need to take the appropriate
    ! action for the radiation here (during the hydro step)

    ! this do loop counts backwards since we want to work downward
    if ( bc(2,1,1).eq.EXT_DIR .and. rad_l2.lt.domlo(2)) then
       do j=domlo(2)-1,rad_l2,-1

          ! zero-gradient catch-all -- this will get the radiation
          ! energy
          rad(rad_l1:rad_h1,j) = rad(rad_l1:rad_h1,j+1)
       enddo
    endif

  end subroutine ca_radfill

  
  subroutine ca_phigravfill(phi,phi_l1,phi_l2, &
                            phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_h1,phi_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
         domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill

end module bc_fill_module
