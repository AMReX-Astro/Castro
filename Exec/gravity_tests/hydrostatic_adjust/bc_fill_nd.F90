module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only: dim
   use amrex_filcc_module, only: amrex_filccn

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine ca_hypfill(adv, adv_lo, adv_hi, &
                        domlo, domhi, delta, xlo, time, bc) bind(C, name="ca_hypfill")

    use amrex_constants_module
    use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UMZ, UTEMP, UEINT, UFS
    use probdata_module, only: hse_rho_top, hse_t_top, hse_X_top, &
         hse_eint_top, hse_p_top
    use network, only: nspec

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: adv_lo(3), adv_hi(3)
    integer, intent(in) :: bc(dim,2,NVAR)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    integer :: n, i, j, k
    real(rt) :: vel

    ! call the generic ghostcell filling routine
    call amrex_filccn(adv_lo, adv_hi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, delta, xlo, bc)

    ! override the generic routine at the top physical boundary
    ! by resetting the velocity to zero there.

    if (dim == 1) then

       if (adv_hi(1) > domhi(1)) then
          if (bc(1,2,UMX) == FOEXTRAP) then
             do k = adv_lo(3), adv_hi(3)
                do j = adv_lo(2), adv_hi(2)
                   do i = domhi(1)+1, adv_hi(1)

                      !adv(i,UMX) = adv(domhi(1),UMX)
                      vel = max(adv(i,j,k,UMX)/adv(i,j,k,URHO), ZERO)
                      adv(i,j,k,URHO)  = hse_rho_top
                      adv(i,j,k,UMX)   = adv(i,j,k,URHO)*vel
                      adv(i,j,k,UMY)   = ZERO
                      adv(i,j,k,UMZ)   = ZERO
                      adv(i,j,k,UTEMP) = hse_T_top
                      adv(i,j,k,UEINT) = hse_rho_top*hse_eint_top
                      adv(i,j,k,UEDEN) = hse_rho_top*hse_eint_top + &
                           HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                      adv(i,j,k,UFS:UFS+nspec-1) = hse_rho_top*hse_X_top(:)
                   end do
                end do
             end do
          end if
       end if

    else if (dim == 2) then

       if (adv_hi(2) > domhi(2)) then
          if (bc(2,2,UMY) == FOEXTRAP) then
             do k = adv_lo(3), adv_hi(3)
                do i = adv_lo(1), adv_hi(1)
                   do j = domhi(2)+1, adv_hi(2)

                      !adv(i,UMX) = adv(domhi(1),UMX)
                      vel = max(adv(i,j,k,UMY)/adv(i,j,k,URHO), ZERO)
                      adv(i,j,k,URHO)  = hse_rho_top
                      adv(i,j,k,UMX)   = ZERO
                      adv(i,j,k,UMY)   = adv(i,j,k,URHO)*vel
                      adv(i,j,k,UMZ)   = ZERO
                      adv(i,j,k,UTEMP) = hse_T_top
                      adv(i,j,k,UEINT) = hse_rho_top*hse_eint_top
                      adv(i,j,k,UEDEN) = hse_rho_top*hse_eint_top + &
                           HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                      adv(i,j,k,UFS:UFS+nspec-1) = hse_rho_top*hse_X_top(:)
                   end do
                end do
             end do
          end if
       end if

    else

       if (adv_hi(3) > domhi(3)) then
          if (bc(3,2,UMZ) == FOEXTRAP) then
             do j = adv_lo(2), adv_hi(2)
                do i = adv_lo(1), adv_hi(1)
                   do k = domhi(3)+1, adv_hi(3)

                      !adv(i,UMX) = adv(domhi(1),UMX)
                      vel = max(adv(i,j,k,UMZ)/adv(i,j,k,URHO), ZERO)
                      adv(i,j,k,URHO)  = hse_rho_top
                      adv(i,j,k,UMX)   = ZERO
                      adv(i,j,k,UMY)   = ZERO
                      adv(i,j,k,UMZ)   = adv(i,j,k,URHO)*vel
                      adv(i,j,k,UTEMP) = hse_T_top
                      adv(i,j,k,UEINT) = hse_rho_top*hse_eint_top
                      adv(i,j,k,UEDEN) = hse_rho_top*hse_eint_top + &
                           HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                      adv(i,j,k,UFS:UFS+nspec-1) = hse_rho_top*hse_X_top(:)
                   end do
                end do
             end do
          end if
       end if

    end if

  end subroutine ca_hypfill

  subroutine ca_denfill(adv, adv_lo, adv_hi, &
                        domlo, domhi, delta, xlo, time, bc) bind(C, name="ca_denfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: adv_lo(3), adv_hi(3)
    integer, intent(in) :: bc(dim,2,*)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

    call amrex_filccn(adv_lo, adv_hi, adv, adv_lo, adv_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_denfill



#ifdef GRAVITY
  subroutine ca_gravxfill(grav, grav_lo, grav_hi, &
                          domlo, domhi, delta, xlo, time, bc) bind(C, name="ca_gravxfill")

    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call amrex_filccn(grav_lo, grav_hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_gravxfill

  subroutine ca_gravyfill(grav, grav_lo, grav_hi, &
                          domlo, domhi, delta, xlo, time, bc) bind(C, name="ca_gravyfill")

    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call amrex_filccn(grav_lo, grav_hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav, grav_lo, grav_hi, &
                          domlo, domhi, delta, xlo, time, bc) bind(C, name="ca_gravzfill")

    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call amrex_filccn(grav_lo, grav_hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_gravzfill

  subroutine ca_phigravfill(phi, phi_lo, phi_hi, &
                            domlo, domhi, delta, xlo, time, bc) bind(C, name="ca_phigravfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

    call amrex_filccn(phi_lo, phi_hi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_phigravfill
#endif


#ifdef REACTIONS
  subroutine ca_reactfill(adv, adv_lo, adv_hi, &
                          domlo, domhi, delta, xlo, time, bc) bind(C, name="ca_reactfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: react_lo(3), react_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: react(react_lo(1):react_hi(1),react_lo(2):react_hi(2),react_lo(3):react_hi(3))

    call amrex_filccn(react_lo, react_hi, react, react_lo, react_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_reactfill
#endif

end module bc_fill_module
