module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")

    use meth_params_module, only : NVAR

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    include 'AMReX_bc_types.fi'
    
    integer :: adv_l1,adv_l2,adv_h1,adv_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    real(rt)         :: state(NVAR)
    real(rt)         :: staten(NVAR)

    integer :: i, j, n
    logical rho_only

    do n = 1,NVAR
       call filcc(adv(adv_l1,adv_l2,n), &
                  adv_l1,adv_l2,adv_h1,adv_h2, &
                  domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    ! The strategy here is to set Dirichlet condition for inflow and
    ! outflow boundaries, and let the Riemann solver sort out the proper
    ! upwinding.  However, this decision makes this routine look
    ! somewhat non-orthodox, in that we need to set external values in
    ! either case....how do we know it's Outflow?  We have to assume
    ! that the setup routines converted Outflow to FOEXTRAP.

    ! Set flag for bc function
    rho_only = .FALSE.

    !     XLO
    if ( (bc(1,1,1).eq.EXT_DIR.or.bc(1,1,1).eq.FOEXTRAP).and. &
         adv_l1.lt.domlo(1)) then
       do i = adv_l1, domlo(1)-1
          do j = adv_l2, adv_h2
             do n=1,NVAR
                state(n) = adv(domlo(1),j,n)
             enddo
             call bcnormal(state,staten,1,+1,rho_only)
             do n=1,NVAR
                adv(i,j,n) = staten(n)
             enddo
          end do
       end do
    end if

    !     XHI
    if ( (bc(1,2,1).eq.EXT_DIR.or.bc(1,2,1).eq.FOEXTRAP).and. &
         adv_h1.gt.domhi(1)) then
       do i = domhi(1)+1, adv_h1
          do j = adv_l2, adv_h2
             do n=1,NVAR
                state(n) = adv(domhi(1),j,n)
             enddo
             call bcnormal(state,staten,1,-1,rho_only)
             do n=1,NVAR
                adv(i,j,n) = staten(n)
             enddo
          end do
       end do
    end if

    !     YLO
    if ( (bc(2,1,1).eq.EXT_DIR.or.bc(2,1,1).eq.FOEXTRAP).and. &
         adv_l2.lt.domlo(2)) then
       do i = adv_l1, adv_h1
          do j = adv_l2, domlo(2)-1
             do n=1,NVAR
                state(n) = adv(i,domlo(2),n)
             enddo
             call bcnormal(state,staten,2,+1,rho_only)
             do n=1,NVAR
                adv(i,j,n) = staten(n)
             enddo
          end do
       end do
    end if

    !     YHI
    if ( (bc(2,2,1).eq.EXT_DIR.or.bc(2,2,1).eq.FOEXTRAP).and. &
         adv_h2.gt.domhi(2)) then
       do i = adv_l1, adv_h1
          do j = domhi(2)+1, adv_h2
             do n=1,NVAR
                state(n) = adv(i,domhi(2),n)
             enddo
             call bcnormal(state,staten,2,-1,rho_only)
             do n=1,NVAR
                adv(i,j,n) = staten(n)
             enddo
          end do
       end do
    end if

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    include 'AMReX_bc_types.fi'
    
    integer :: adv_l1,adv_l2,adv_h1,adv_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2)
    logical rho_only
    integer :: i,j

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    rho_only = .TRUE.

    !     XLO
    if ( (bc(1,1,1).eq.EXT_DIR.or.bc(1,1,1).eq.FOEXTRAP).and.adv_l1.lt.domlo(1)) then
       do i = adv_l1, domlo(1)-1
          do j = adv_l2, adv_h2
             call bcnormal(adv(domlo(1),j),adv(i,j),1,+1,rho_only)
          end do
       end do
    end if

    !     XHI
    if ( (bc(1,2,1).eq.EXT_DIR.or.bc(1,2,1).eq.FOEXTRAP).and.adv_h1.gt.domhi(1)) then
       do i = domhi(1)+1, adv_h1
          do j = adv_l2, adv_h2
             call bcnormal(adv(domhi(1),j),adv(i,j),1,-1,rho_only)
          end do
       end do
    end if

    !     YLO
    if ( (bc(2,1,1).eq.EXT_DIR.or.bc(2,1,1).eq.FOEXTRAP).and.adv_l2.lt.domlo(2)) then
       do i = adv_l1, adv_h1
          do j = adv_l2, domlo(2)-1
             call bcnormal(adv(i,domlo(2)),adv(i,j),2,+1,rho_only)
          end do
       end do
    end if

    !     YHI
    if ( (bc(2,2,1).eq.EXT_DIR.or.bc(2,2,1).eq.FOEXTRAP).and.adv_h2.gt.domhi(2)) then
       do i = adv_l1, adv_h1
          do j = domhi(2)+1, adv_h2
             call bcnormal(adv(i,domhi(2)),adv(i,j),2,-1,rho_only)
          end do
       end do
    end if
  end subroutine ca_denfill



  subroutine bcnormal(u_int,u_ext,dir,sgn,rho_only)

    use probdata_module
    use eos_module, only : gamma_const
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT
    
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: u_int(*),u_ext(*)
    logical rho_only
    integer :: dir,sgn
    real(rt)         :: rho, rhou(2), eden, T
    integer :: n,t1,i

    ! for the Sedov problem, we will always set the state to the ambient conditions
    if (rho_only .EQV. .TRUE. ) then

       u_ext(1) = dens_ambient

    else

       ! First set everything from internal data (this is probably a bad
       ! thing to do...)  That is, we should have explicit boundary data
       ! for advected fields and species

       do i=1,NVAR
          u_ext(i) = u_int(i)
       enddo

       u_ext(URHO)   = dens_ambient
       u_ext(UMX)    = 0.e0_rt
       u_ext(UMY)    = 0.e0_rt
       u_ext(UEDEN)  = p_ambient/(gamma_const-1.e0_rt)
       u_ext(UEINT)  = u_ext(UEDEN)

    endif

  end subroutine bcnormal

end module bc_fill_module
