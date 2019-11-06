! problem-specific Fortran derive routines go here

subroutine derextheating(h, h_lo, h_hi, ncomp_h, &
                         s, s_lo, s_hi, ncomp_s, &
                         lo, hi, domlo, domhi, &
                         dx, time) &
                         bind(C, name='derextheating')
    ! Calculate the external heating

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE, M_PI
    use prob_params_module, only : problo
    use probdata_module, only : heating_factor

    implicit none

    integer,  intent(in   ) :: h_lo(3), h_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt), intent(inout) :: h(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3),ncomp_h)
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ncomp_s)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: ncomp_h ! == 1
    integer,  intent(in   ), value :: ncomp_s
    real(rt), intent(in   ), value :: time

    integer :: i, j, k
    real(rt) :: fheat, y

    !$gpu

    h(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:ncomp_h) = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j) + HALF) * dx(2) 
          if (y < 1.125e0_rt * 4.e8_rt) then 

            fheat = sin(8.e0_rt * M_PI * (y/ 4.e8_rt - ONE))

            do i = lo(1), hi(1)
    
               ! Source terms
               h(i,j,k,1) = heating_factor * fheat
    
            end do
          endif

       end do
    end do

end subroutine derextheating


subroutine deradinvariant(A, A_lo, A_hi, ncomp_A, &
                         s, s_lo, s_hi, ncomp_s, &
                         lo, hi, domlo, domhi, &
                         dx, time) &
                         bind(C, name='deradinvariant')
    ! calculate the adiabatic invariant

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE, M_PI
    use network, only : nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX

    implicit none

    integer,  intent(in   ) :: A_lo(3), A_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt), intent(inout) :: A(A_lo(1):A_hi(1),A_lo(2):A_hi(2),A_lo(3):A_hi(3),ncomp_A)
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ncomp_s)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: ncomp_A ! == 1
    integer,  intent(in   ), value :: ncomp_s
    real(rt), intent(in   ), value :: time

    integer :: i, j, k
    real(rt) :: rhoInv
    type(eos_t) :: eos_state

    !$gpu

    A(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:ncomp_A) = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                rhoInv = ONE / s(i,j,k,URHO)
   
                eos_state % e     = s(i,j,k,UEINT) * rhoInv
                eos_state % T     = s(i,j,k,UTEMP)
                eos_state % rho   = s(i,j,k,URHO)
                eos_state % xn  = s(i,j,k,UFS:UFS+nspec-1) * rhoInv
                eos_state % aux = s(i,j,k,UFX:UFX+naux-1) * rhoInv
   
                call eos(eos_input_re, eos_state)

                A(i,j,k,1) = eos_state % p / s(i,j,k,URHO)**(5.0e0_rt / 3.0e0_rt)
    
            end do
       end do
    end do

end subroutine deradinvariant


subroutine derenthalpyfluxy(f, f_lo, f_hi, ncomp_f, &
                         s, s_lo, s_hi, ncomp_s, &
                         lo, hi, domlo, domhi, &
                         dx, time) &
                         bind(C, name='derenthalpyfluxy')
    ! vertical enthalpy flux

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE, M_PI
    use network, only : nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX, UMY

    implicit none

    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt), intent(inout) :: f(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),ncomp_f)
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ncomp_s)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: ncomp_f ! == 1
    integer,  intent(in   ), value :: ncomp_s
    real(rt), intent(in   ), value :: time

    integer :: i, j, k
    real(rt) :: rhoInv
    type(eos_t) :: eos_state

    !$gpu

    f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:ncomp_f) = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                rhoInv = ONE / s(i,j,k,URHO)
   
                eos_state % e     = s(i,j,k,UEINT) * rhoInv
                eos_state % T     = s(i,j,k,UTEMP)
                eos_state % rho   = s(i,j,k,URHO)
                eos_state % xn  = s(i,j,k,UFS:UFS+nspec-1) * rhoInv
                eos_state % aux = s(i,j,k,UFX:UFX+naux-1) * rhoInv
   
                call eos(eos_input_re, eos_state)

                f(i,j,k,1) = eos_state % h * s(i,j,k,UMY) / s(i,j,k,URHO)
    
            end do
       end do
    end do

end subroutine derenthalpyfluxy


subroutine derkinengfluxy(f, f_lo, f_hi, ncomp_f, &
                         s, s_lo, s_hi, ncomp_s, &
                         lo, hi, domlo, domhi, &
                         dx, time) &
                         bind(C, name='derkinengfluxy')
    ! vertical kinetic energy flux

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF

    implicit none

    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt), intent(inout) :: f(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),ncomp_f)
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ncomp_s)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: ncomp_f ! == 1
    integer,  intent(in   ), value :: ncomp_s
    real(rt), intent(in   ), value :: time

    integer :: i, j, k

    !$gpu

    f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:ncomp_f) = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
            do i = lo(1), hi(1)

                f(i,j,k,1) = HALF / s(i,j,k,1) * sum(s(i,j,k,2:4)**2) *s(i,j,k,3) 
    
            end do
       end do
    end do

end subroutine derkinengfluxy
