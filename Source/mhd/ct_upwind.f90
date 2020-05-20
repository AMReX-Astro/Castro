module ct_upwind
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_fort_module, only : rt => amrex_real
  use hlld_solver, only : hlld
  use meth_params_module
  use mhd_state_module

  implicit none

  private primtocons

  ! note: in this module, we use left and right to mean with respect
  ! to the interface.  So qleft, uleft, ul, ... are the left state on
  ! an interface and qright, uright, ur, ... are the right state on an
  ! interface.

  ! Miniati and Martin use + and - with respect to a zone center, so
  ! u+ would be on the right edge of a zone, which is the left state
  ! to an interface.

contains

  subroutine PrimToCons(lo, hi, q, q_lo, q_hi, u, u_lo, u_hi) bind(C, name="PrimToCons")

    ! calculate the conserved variables from the primitive

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module
    use eos_type_module, only : eos_t, eos_input_rp
    use eos_module, only: eos
    use network, only: nspec

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)

    real(rt), intent(in) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(out) :: u(u_lo(1):u_hi(1), u_lo(2):u_hi(2), u_lo(3):u_hi(3), NVAR+3)

    integer :: i ,j ,k

    type(eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             u(i,j,k,URHO)  = q(i,j,k,QRHO)
             u(i,j,k,UMX)    = q(i,j,k,QRHO)*q(i,j,k,QU)
             u(i,j,k,UMY)    = q(i,j,k,QRHO)*q(i,j,k,QV)
             u(i,j,k,UMZ)    = q(i,j,k,QRHO)*q(i,j,k,QW)

             eos_state % rho = q(i, j, k, QRHO)
             eos_state % p   = q(i, j, k, QPRES)
             eos_state % T   = 100.d0 !dummy initial g.
             eos_state % xn  = q(i, j, k, QFS:QFS+nspec-1)

             call eos(eos_input_rp, eos_state)

             u(i,j,k,UEDEN) = eos_state % rho * eos_state % e &
                  + 0.5d0*q(i,j,k,QRHO)*dot_product(q(i,j,k,QU:QW),q(i,j,k,QU:QW)) &
                  + 0.5d0*(dot_product(q(i,j,k,QMAGX:QMAGZ),q(i,j,k,QMAGX:QMAGZ)))

             u(i,j,k,UEINT) = eos_state % rho * eos_state % e
             u(i,j,k,UTEMP) = eos_state % T

             u(i,j,k,UMAGX:UMAGZ) = q(i,j,k,QMAGX:QMAGZ)

             ! species
             u(i,j,k,UFS:UFS-1+nspec) = q(i,j,k,QRHO) * q(i,j,k,QFS:QFS-1+nspec)

          enddo
       enddo
    enddo
  end subroutine PrimToCons


  subroutine ConsToPrim(q, u)

    ! calculate the primitive variables from the conserved

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_re
    use network, only : nspec

    implicit none

    real(rt), intent(in) :: u(NVAR+3)
    real(rt), intent(out) :: q(NQ)


    type (eos_t) :: eos_state

    q(QRHO)  = u(URHO)
    q(QU)    = u(UMX)/q(QRHO)
    q(QV)    = u(UMY)/q(QRHO)
    q(QW)    = u(UMZ)/q(QRHO)
    q(QREINT) = u(UEDEN) - &
         0.5d0*q(QRHO)*dot_product(q(QU:QW), q(QU:QW)) - &
         0.5d0*dot_product(u(UMAGX:UMAGZ), u(UMAGX:UMAGZ)) 

    ! species
    q(QFS:QFS-1+nspec) = u(UFS:UFS-1+nspec)/u(URHO)

    eos_state % rho = q(QRHO)
    eos_state % e   = q(QREINT)/eos_state % rho
    eos_state % xn  = q(QFS:QFS+nspec-1)
    eos_state % T   = 100.d0  ! initial guess

    call eos(eos_input_re, eos_state)

    q(QTEMP) = eos_state % T
    q(QPRES) = eos_state % p
    q(QMAGX:QMAGZ) = u(UMAGX:UMAGZ)

  end subroutine ConsToPrim

  subroutine corner_couple(w_lo, w_hi, &
                           qr_out, qro_lo, qro_hi, &
                           ql_out, qlo_lo, qlo_hi, &
                           ur, ur_lo, ur_hi, &
                           ul, ul_lo, ul_hi, &
                           flxd2, flxd2_lo, flxd2_hi, &
                           Ed1, ed1_lo, ed1_hi, &
                           Ed3, ed3_lo, ed3_hi, &
                           d1, d2, d3, &
                           dx, dt) bind(C, name="corner_couple")

    ! take conservative interface states ul and ur and update them
    ! with with the transverse flux difference (corner coupling) to
    ! produce ul_out and ur_out
    !
    ! This implements MM step 3 of the CTU algorithm.
    ! the normal direction (for the interface states) is d1
    ! the transverse direction (for the flux difference) is d2

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module
    use network, only : nspec

    implicit none

    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: qlo_lo(3), qlo_hi(3)
    integer, intent(in) :: qro_lo(3), qro_hi(3)
    integer, intent(in) :: ur_lo(3), ur_hi(3)
    integer, intent(in) :: ul_lo(3), ul_hi(3)
    integer, intent(in) :: flxd2_lo(3), flxd2_hi(3)
    integer, intent(in) :: ed1_lo(3), ed1_hi(3)
    integer, intent(in) :: ed3_lo(3), ed3_hi(3)
    integer, intent(in), value :: d1, d2, d3
    real(rt), intent(in), value :: dx
    real(rt), intent(in), value :: dt

    real(rt), intent(in) :: ur(ur_lo(1):ur_hi(1), ur_lo(2):ur_hi(2), ur_lo(3):ur_hi(3), NVAR+3)
    real(rt), intent(in) :: ul(ul_lo(1):ul_hi(1), ul_lo(2):ul_hi(2), ul_lo(3):ul_hi(3), NVAR+3)

    real(rt), intent(out) :: flxd2(flxd2_lo(1):flxd2_hi(1),flxd2_lo(2):flxd2_hi(2),flxd2_lo(3):flxd2_hi(3),NVAR+3)

    real(rt), intent(in) :: Ed1(ed1_lo(1):ed1_hi(1),ed1_lo(2):ed1_hi(2),ed1_lo(3):ed1_hi(3))
    real(rt), intent(in) :: Ed3(ed3_lo(1):ed3_hi(1),ed3_lo(2):ed3_hi(2),ed3_lo(3):ed3_hi(3))

    real(rt), intent(out) :: qr_out(qro_lo(1):qro_hi(1),qro_lo(2):qro_hi(2),qro_lo(3):qro_hi(3),NQ)
    real(rt), intent(out) :: ql_out(qlo_lo(1):qlo_hi(1),qlo_lo(2):qlo_hi(2),qlo_lo(3):qlo_hi(3),NQ)

    real(rt) :: cdtdx
    integer  :: i, j, k, n

    ! cl and cr are the offsets to the indices for the conserved state fluxes
    ! they will be offset in d2 to capture the flux difference
    integer :: cl(3), cr(3)

    ! b is for indexing into the magnetic field that is in the d1
    ! (normal) direction
    integer :: b(3)

    ! for indexing the electric field
    integer :: err(3), elr(3), erl(3), ell(3)

    integer :: UMAGD1, UMAGD2, UMAGD3   !UMAGD1 corresponds to d1, and UMAGD2 to d2, UMAGD3 to d3
    integer :: sgn

    real(rt) :: utmp(NVAR+3)

    ! update the state on interface direction d1 with the input flux in direction d2

    ! for the flux difference, F_r - F_l, we need to shift the indices in the first flux (F_r)
    ! in d2 to get a difference across the interface.  We also need to shift by a zone in d1
    ! for the left interface.  cr(:) and cl(:) will hold these shifts.
    cr(:) = 0
    cl(:) = 0

    ! the first term of the flxd2 substraction is shifted by 1 on the direction d2
    cr(d2) = 1

    ! for the normal B component
    b(:) = 0
    b(d2) = 1

    ! err will capture the right state in both transverse directions
    ! (e.g. Ez_{i+1/2,j+1/2,k})
    err(:) = 0
    err(d2) = 1
    err(d3) = 1

    ! elr will capture the right state in the second transverse direction
    ! (e.g. Ez_{i-1/2,j+1/2,k})
    elr(:) = 0
    elr(d3) = 1

    ! erl will capture the right state in the first transverse direction
    ! (e.g. Ez_{i+1/2,j-1/2,k})
    erl(:) = 0
    erl(d2) = 1

    ! ell is the lower-left E in both directions
    ! (e.g. Ez_{i-1/2,j-1/2,k})
    ell(:) = 0

    sgn = epsilon_ijk(d1, d2, d3)
    cdtdx = dt/(3.d0*dx)

    UMAGD1 = UMAGX - 1 + d1
    UMAGD2 = UMAGX - 1 + d2
    UMAGD3 = UMAGX - 1 + d3

    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             ! first the conserved state

             ! right interface (e.g. U_{i-1/2,j,k,R} or the "-" state in MM notation)
             ! MM Eq. 37

             do n = 1, NVAR
                if (n == UTEMP) cycle

                utmp(n) = ur(i,j,k,n) - &
                     cdtdx*(flxd2(i+cr(1),j+cr(2),k+cr(3),n) - flxd2(i,j,k,n))
             end do

             ! now magnetic fields

             ! right state on the interface (e.g. B_{i-1/2,j,k,R} or `-` in MM notation)

             ! d1 -- this is perpendicular to the face, MM Eq. 38
             ! (note MM Eq. 38 has a sign error) e.g., for d1 = x and
             ! d2 = y, this gets updated as
             !
             ! Bx|y_{i-1/2,j,k,R} = Bx_{i-1/2,j,k) -
             !     dt/3dx (Ez_{i-1/2,j+1/2,k) - Ez_{i-1/2,j-1/2,k})
             !
             ! we use b(:) to captured the j+1/2 indexing into Ez

             utmp(UMAGD1) = ur(i,j,k,UMAGD1) - sgn * cdtdx * &
                  (Ed3(i+b(1),j+b(2),k+b(3)) - Ed3(i,j,k))

             ! d3 -- this is in the plane of the face, MM Eq. 39
             ! e.g.,g for d1 = x, and d2 = y, this gets updated as
             !
             ! Bz|y_{i-1/2,j,k,R} = Bz_{i-1/2,j,k} + 1/2 dt/3dx
             !     (Ex_{i,j+1/2,k+1/2} - Ex_{i,j-1/2,k+1/2} +
             !      Ex_{i,j+1/2,k-1/2} - Ex_{i,j-1/2,k-1/2})
             !
             ! we use err(:) for the first E term, elr(:) for the
             ! second, and erl(:) for the third

             utmp(UMAGD3) = ur(i,j,k,UMAGD3) + sgn * 0.5_rt * cdtdx * &
                  ((Ed1(i+err(1),j+err(2),k+err(3)) - Ed1(i+elr(1),j+elr(2),k+elr(3))) + &
                   (Ed1(i+erl(1),j+erl(2),k+erl(3)) - Ed1(i,j,k)))

             ! the component pointing in the transverse update direction, d2, is unchanged
             utmp(UMAGD2) = ur(i,j,k,UMAGD2)

             call ConsToPrim(qr_out(i,j,k,:), utmp)

          end do
       end do
    end do


    ! left interface (e.g., U_{i-1/2,j,k,L} or the "+" state in MM notation)
    ! note: this uses information one zone to the left in d1
    cl(d1) = -1
    cr(d1) = -1

    ! The in-plane B component at B_{i-1/2,j,k,L} uses the information one zone to the left
    ! in direction d1
    err(d1) = -1
    elr(d1) = -1
    erl(d1) = -1
    ell(d1) = -1

    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             ! conservative state

             do n = 1, NVAR
                if (n == UTEMP) cycle

                utmp(n) = ul(i,j,k,n) - &
                     cdtdx*(flxd2(i+cr(1),j+cr(2),k+cr(3),n) - &
                            flxd2(i+cl(1),j+cl(2),k+cl(3),n))
             end do

             ! left state on the interface (e.g. B_{i-1/2,j,k,L} or `+` in MM notation)

             utmp(UMAGD1) = ul(i,j,k,UMAGD1) - sgn * cdtdx * &
                  (Ed3(i+b(1),j+b(2),k+b(3)) - Ed3(i,j,k))

             utmp(UMAGD3) = ul(i,j,k,UMAGD3) + sgn * 0.5_rt * cdtdx * &
                  ((Ed1(i+err(1),j+err(2),k+err(3)) - Ed1(i+elr(1),j+elr(2),k+elr(3))) + &
                   (Ed1(i+erl(1),j+erl(2),k+erl(3)) - Ed1(i+ell(1),j+ell(2),k+ell(3))))

             utmp(UMAGD2) = ul(i,j,k,UMAGD2)

             call ConsToPrim(ql_out(i,j,k,:), utmp)

          enddo
       enddo
    enddo
  end subroutine corner_couple


  subroutine half_step(w_lo, w_hi, &
                       qr_out, qro_lo, qro_hi, &
                       ql_out, qlo_lo, qlo_hi, &
                       ur, ur_lo, ur_hi, &
                       ul, ul_lo, ul_hi, &
                       flxd1, flxd1_lo, flxd1_hi, &
                       flxd2, flxd2_lo, flxd2_hi, &
                       Ed, ed_lo, ed_hi, &
                       Ed1, ed1_lo, ed1_hi, &
                       Ed2, ed2_lo, ed2_hi, &
                       d, d1, d2, dx, dt) bind(C, name="half_step")

    ! Final transverse flux corrections to the conservative state

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UMZ, URHO, UEINT, UFS
    use network, only: nspec

    implicit none

    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: qro_lo(3), qro_hi(3)
    integer, intent(in) :: qlo_lo(3), qlo_hi(3)
    integer, intent(in) :: ur_lo(3), ur_hi(3)
    integer, intent(in) :: ul_lo(3), ul_hi(3)
    integer, intent(in) :: ed_lo(3), ed_hi(3)
    integer, intent(in) :: ed1_lo(3), ed1_hi(3)
    integer, intent(in) :: ed2_lo(3), ed2_hi(3)

    integer, intent(in) :: flxd1_lo(3), flxd1_hi(3)
    integer, intent(in) :: flxd2_lo(3), flxd2_hi(3)

    real(rt), intent(inout) :: qr_out(qro_lo(1):qro_hi(1), qro_lo(2):qro_hi(2), qro_lo(3):qro_hi(3), NQ)
    real(rt), intent(inout) :: ql_out(qlo_lo(1):qlo_hi(1), qlo_lo(2):qlo_hi(2), qlo_lo(3):qlo_hi(3), NQ)
    real(rt), intent(in)    :: ur(ur_lo(1):ur_hi(1), ur_lo(2):ur_hi(2), ur_lo(3):ur_hi(3), NVAR+3)
    real(rt), intent(in)    :: ul(ul_lo(1):ul_hi(1), ul_lo(2):ul_hi(2), ul_lo(3):ul_hi(3), NVAR+3)

    real(rt), intent(in) :: Ed(ed_lo(1):ed_hi(1),ed_lo(2):ed_hi(2),ed_lo(3):ed_hi(3))
    real(rt), intent(in) :: Ed1(ed1_lo(1):ed1_hi(1),ed1_lo(2):ed1_hi(2),ed1_lo(3):ed1_hi(3))
    real(rt), intent(in) :: Ed2(ed2_lo(1):ed2_hi(1),ed2_lo(2):ed2_hi(2),ed2_lo(3):ed2_hi(3))

    real(rt), intent(in), value :: dx, dt !dx will be dx, dy or dz
    integer, intent(in), value :: d, d1, d2 ! following notation of eq. 44

    real(rt), intent(in)  :: flxd1(flxd1_lo(1):flxd1_hi(1),flxd1_lo(2):flxd1_hi(2),flxd1_lo(3):flxd1_hi(3),NVAR+3)
    real(rt), intent(in)  :: flxd2(flxd2_lo(1):flxd2_hi(1),flxd2_lo(2):flxd2_hi(2),flxd2_lo(3):flxd2_hi(3),NVAR+3)

    integer  :: i ,j ,k, n

    ! for the shift in i,j,k

    ! c1l, c1r are for indexing flxd1 offsets, c2l, c2r are for flxd2
    integer  :: c1l(3), c1r(3), c2l(3), c2r(3)
    integer :: a1(3), a2(3)
    integer :: err(3), erl(3), elr(3), ell(3)
    integer :: e1rr(3), e1rl(3), e1lr(3)
    integer :: e2rr(3), e2rl(3), e2lr(3)

    real(rt) :: hdtdx
    integer :: UMAGD, UMAGD1, UMAGD2
    integer :: sgn

    real(rt) :: utmp(NVAR+3)

    hdtdx = 0.5_rt * dt/dx
    sgn = -1 * epsilon_ijk(d, d1, d2)

    UMAGD = UMAGX - 1 + d
    UMAGD1 = UMAGX - 1 + d1
    UMAGD2 = UMAGX - 1 + d2

    c1l(:) = 0
    c1r(:) = 0
    c2l(:) = 0
    c2r(:) = 0

    c1r(d1) = 1  ! add +1 to the d1 direction in the first flxd1 term of the subtraction
    c2r(d2) = 1  ! add +1 to the d2 direction in the first flxd2 term of the subtraction

    ! for Ed
    err(:) = 0
    err(d1) = 1
    err(d2) = 1

    erl(:) = 0
    erl(d1) = 1

    elr(:) = 0
    elr(d2) = 1

    ell(:) = 0

    ! for Ed1
    e1rr(:) = 0
    e1rr(d) = 1
    e1rr(d2) = 1

    e1lr(:) = 0
    e1lr(d2) = 1

    e1rl(:) = 0
    e1rl(d) = 1

    ! for Ed2
    e2rr(:) = 0
    e2rr(d) = 1
    e2rr(d1) = 1

    e2lr(:) = 0
    e2lr(d1) = 1

    e2rl(:) = 0
    e2rl(d) = 1


    ! for the normal component of B
    a1(:) = 0
    a2(:) = 0

    a1(d2) = 1 ! shift on first term of Ed1 substraction, in d2 direction
    a2(d1) = 1 ! shift on first term of Ed2 substraction, in d1 direction


    ! right interface (e.g. U_{i-1/2,j,k,R} or the "-" state in MM notation)
    ! MM Eq. 44

    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             ! first the conservative state

             do n = 1, NVAR
                if (n == UTEMP) cycle

                utmp(n) = ur(i,j,k,n) - &
                     hdtdx * (flxd1(i+c1r(1),j+c1r(2),k+c1r(3),n) - flxd1(i,j,k,n)) - &
                     hdtdx * (flxd2(i+c2r(1),j+c2r(2),k+c2r(3),n) - flxd2(i,j,k,n))

             end do

             ! now the magnetic field

             ! right state on the interface (e.g. B_{i-1/2,j,k,R} or `-` in MM notation)

             ! Bd -- this is perpendicular to the face. Note MM eq.45
             ! in Miniati has a sign error in the epsilon term

             utmp(UMAGD) = ur(i,j,k,UMAGD) - sgn * hdtdx * &
                  ((Ed1(i+a1(1),j+a1(2),k+a1(3)) - Ed1(i,j,k)) - &
                   (Ed2(i+a2(1),j+a2(2),k+a2(3)) - Ed2(i,j,k)))

             ! Bd1 -- this is one of the components of B in the plane of the face d
             ! Eq.46 in Miniati

             utmp(UMAGD1) = ur(i,j,k,UMAGD1) + sgn * hdtdx * &
                  ((Ed(i+err(1),j+err(2),k+err(3)) - Ed(i+erl(1),j+erl(2),k+erl(3))) + &
                   (Ed(i+elr(1),j+elr(2),k+elr(3)) - Ed(i+ell(1),j+ell(2),k+ell(3))) - &
                   (Ed2(i+e2rr(1),j+e2rr(2),k+e2rr(3)) - Ed2(i+e2lr(1),j+e2lr(2),k+e2lr(3))) - &
                   (Ed2(i+e2rl(1),j+e2rl(2),k+e2rl(3)) - Ed2(i+ell(1),j+ell(2),k+ell(3))))

             ! Bd2 -- this is the other component of B in the plane of the face d
             ! Eq. 46 in Miniati

             utmp(UMAGD2) = ur(i,j,k,UMAGD2) - sgn * hdtdx * &
                  ((Ed(i+err(1),j+err(2),k+err(3)) - Ed(i+elr(1),j+elr(2),k+elr(3))) + &
                   (Ed(i+erl(1),j+erl(2),k+erl(3)) - Ed(i+ell(1),j+ell(2),k+ell(3))) - &
                   (Ed1(i+e1rr(1),j+e1rr(2),k+e1rr(3)) - Ed1(i+e1lr(1),j+e1lr(2),k+e1lr(3))) - &
                   (Ed1(i+e1rl(1),j+e1rl(2),k+e1rl(3)) - Ed1(i+ell(1),j+ell(2),k+ell(3))))

             ! convert to primitive
             call ConsToPrim(qr_out(i,j,k,:), utmp)

          end do
       end do
    end do

    ! for the left state U components on the face dir, the flux
    ! difference is in the zone to the left
    c1r(d) = -1
    c1l(d) = -1
    c2r(d) = -1
    c2l(d) = -1

    ! left state on the interface (e.g., B_{i-1/2,j,k,L} or `+` in MM notation)

    ! The in-plane B component at B_{i-1/2,j,k,L} uses the information one zone to the left
    ! in direction d1

    err(d) = err(d) - 1
    erl(d) = erl(d) - 1
    elr(d) = elr(d) - 1
    ell(d) = ell(d) - 1

    e1rr(d) = e1rr(d) - 1
    e1rl(d) = e1rl(d) - 1
    e1lr(d) = e1lr(d) - 1

    e2rr(d) = e2rr(d) - 1
    e2rl(d) = e2rl(d) - 1
    e2lr(d) = e2lr(d) - 1


    ! left interface (e.g., U_{i+1/2,j,k,L} or the "+" state in MM notation)
    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             ! first the conserved state

             do n = 1, NVAR
                if (n == UTEMP) cycle

                utmp(n) = ul(i,j,k,n) - &
                     hdtdx * (flxd1(i+c1r(1),j+c1r(2),k+c1r(3),n) - &
                              flxd1(i+c1l(1),j+c1l(2),k+c1l(3),n)) - &
                     hdtdx * (flxd2(i+c2r(1),j+c2r(2),k+c2r(3),n) - &
                              flxd2(i+c2l(1),j+c2l(2),k+c2l(3),n))

             end do

             ! now the B fields

             ! for the + case, the shifts mentioned above in b6, b5, and b4
             ! also correspond to the 1st, 2nd and 4th, and 3rd term respectevely

             ! Bd -- this is perpendicular to the face (MM Eq. 45 with sign fix)

             ! this is the same face as the right state, so the update the identical
             utmp(UMAGD) = ul(i,j,k,UMAGD) - sgn * hdtdx * &
                  ((Ed1(i+a1(1),j+a1(2),k+a1(3)) - Ed1(i,j,k)) - &
                   (Ed2(i+a2(1),j+a2(2),k+a2(3)) - Ed2(i,j,k)))

             ! Bd1 -- first component on face d, eq. 46 in Miniati

             utmp(UMAGD1) = ul(i,j,k,UMAGD1) + sgn * hdtdx * &
                  ((Ed(i+err(1),j+err(2),k+err(3)) - Ed(i+erl(1),j+erl(2),k+erl(3))) + &
                   (Ed(i+elr(1),j+elr(2),k+elr(3)) - Ed(i+ell(1),j+ell(2),k+ell(3))) - &
                   (Ed2(i+e2rr(1),j+e2rr(2),k+e2rr(3)) - Ed2(i+e2lr(1),j+e2lr(2),k+e2lr(3))) - &
                   (Ed2(i+e2rl(1),j+e2rl(2),k+e2rl(3)) - Ed2(i+ell(1),j+ell(2),k+ell(3))))

             ! Bd2 -- second component on face d, eq. 46 in Miniati

             utmp(UMAGD2) = ul(i,j,k,UMAGD2) - sgn * hdtdx * &
                  ((Ed(i+err(1),j+err(2),k+err(3)) - Ed(i+elr(1),j+elr(2),k+elr(3))) + &
                   (Ed(i+erl(1),j+erl(2),k+erl(3)) - Ed(i+ell(1),j+ell(2),k+ell(3))) - &
                   (Ed1(i+e1rr(1),j+e1rr(2),k+e1rr(3)) - Ed1(i+e1lr(1),j+e1lr(2),k+e1lr(3))) - &
                   (Ed1(i+e1rl(1),j+e1rl(2),k+e1rl(3)) - Ed1(i+ell(1),j+ell(2),k+ell(3))))

             ! convert to primitive
             call ConsToPrim(ql_out(i,j,k,:), utmp)

          end do
       end do
    end do

  end subroutine


  subroutine prim_half(w_lo, w_hi, &
                       q2D, q2_lo, q2_hi, &
                       q, q_lo, q_hi, &
                       flxx, flxx_lo, flxx_hi, &
                       flxy, flxy_lo, flxy_hi, &
                       flxz, flxz_lo, flxz_hi, &
                       dx, dy, dz, dt) bind(C, name="prim_half")

    ! Find the 2D corrected primitive variables

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR

    implicit none

    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: flxx_lo(3), flxx_hi(3)
    integer, intent(in) :: flxy_lo(3), flxy_hi(3)
    integer, intent(in) :: flxz_lo(3), flxz_hi(3)

    real(rt), intent(in)  :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in)  :: flxx(flxx_lo(1):flxx_hi(1),flxx_lo(2):flxx_hi(2),flxx_lo(3):flxx_hi(3),NVAR+3)
    real(rt), intent(in)  :: flxy(flxy_lo(1):flxy_hi(1),flxy_lo(2):flxy_hi(2),flxy_lo(3):flxy_hi(3),NVAR+3)
    real(rt), intent(in)  :: flxz(flxz_lo(1):flxz_hi(1),flxz_lo(2):flxz_hi(2),flxz_lo(3):flxz_hi(3),NVAR+3)

    real(rt), intent(out) :: q2D(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NQ)

    real(rt)  :: divF(NVAR+3)
    real(rt)  :: divF_q(NQ)
    real(rt), value  :: dx, dy, dz, dt

    integer   :: i, j, k

    do k = w_lo(3),w_hi(3)
       do j = w_lo(2),w_hi(2)
          do i = w_lo(1),w_hi(1)

             divF(:) = (flxx(i+1,j,k,:) - flxx(i,j,k,:)) / dx + &
                       (flxy(i,j+1,k,:) - flxy(i,j,k,:)) / dy + &
                       (flxz(i,j,k+1,:) - flxz(i,j,k,:)) / dz

             ! that is a flux of conserved variables -- transform it to primitive
             call qflux(divF_q, divF, q(i,j,k,:))

             ! Right below eq. 48
             q2D(i,j,k,:) = q(i,j,k,:) - 0.5d0*dt * divF_q

          enddo
       enddo
    enddo
  end subroutine prim_half


  subroutine qflux(qflx,flx,q)

    ! Calculate the C to P Jacobian applied to the fluxes

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module !,only : QRHO, QU, QV, QW, QPRES, QMAGX, QMAGY, QMAGZ, QVAR, NVAR
    use eos_module, only : eos
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only : nspec

    implicit none

    ! this is step 10 in the paper, just after Eq. 48

    ! this implements dW/dU . qflux, where dW/dU is the Jacobian of
    ! the primitive quantities (W) with respect to conserved quantities (U)

    real(rt), intent(in) :: flx(NVAR+3), q(NQ)
    real(rt), intent(out) :: qflx(NQ)
    real(rt) :: dedp, dedrho, totalE

    type (eos_t) :: eos_state

    qflx = 0.d0
    qflx(QRHO) = flx(URHO)
    qflx(QU) = ( flx(UMX) - flx(URHO) * q(QU) )/q(QRHO)
    qflx(QV) = ( flx(UMY) - flx(URHO) * q(QV) )/q(QRHO)
    qflx(QW) = ( flx(UMZ) - flx(URHO) * q(QW) )/q(QRHO)

    qflx(QFS:QFS+nspec-1) = ( flx(UFS:UFS+nspec-1) - flx(URHO) * q(QFS:QFS+nspec-1) )/q(QRHO)

    eos_state % rho = q(QRHO)
    eos_state % p   = q(QPRES)
    eos_state % T   = 100.d0 !dummy initial guess
    eos_state % xn  = q(QFS:QFS+nspec-1)

    call eos(eos_input_rp, eos_state)

    dedrho = eos_state % dedr - eos_state % dedT * eos_state % dPdr * 1.0d0/eos_state % dPdT
    dedp = eos_state % dedT * 1.0d0/eos_state % dPdT

    qflx(QPRES) = ( -q(QMAGX)*flx(UMAGX) - q(QMAGY)*flx(UMAGY) - q(QMAGZ)*flx(UMAGZ) + &
         flx(UEDEN) - flx(UMX)*q(QU) - flx(UMY)*q(QV) - &
         flx(UMZ)*q(QW) + flx(URHO)*(0.5*(q(QU)**2+q(QV)**2+q(QW)**2) - &
         eos_state % e -q(QRHO)*dedrho) ) / ( dedp * q(QRHO))

    qflx(QMAGX) = flx(UMAGX)
    qflx(QMAGY) = flx(UMAGY)
    qflx(QMAGZ) = flx(UMAGZ)

    qflx(QTEMP) = 0.0_rt

  end subroutine qflux

end module ct_upwind
