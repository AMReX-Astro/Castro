subroutine ca_ext_src(lo, hi, &
                      old_state, os_lo, os_hi, &
                      new_state, ns_lo, ns_hi, &
                      src, src_lo, src_hi, &
                      problo, dx, time, dt) bind(C, name='ca_ext_src')

  use amrex_constants_module, only: THIRD
  use meth_params_module, only : NVAR, NSRC, UEDEN, UEINT, URHO, UTEMP, UFS
  use network, only: network_species_index
  use amrex_fort_module, only : rt => amrex_real
  use probdata_module, only : ih1, ic12, in14, io16

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: os_lo(3), os_hi(3)
  integer,  intent(in   ) :: ns_lo(3), ns_hi(3)
  integer,  intent(in   ) :: src_lo(3), src_hi(3)
  real(rt), intent(in   ) :: old_state(os_lo(1):os_hi(1),os_lo(2):os_hi(2),os_lo(3):os_hi(3),NVAR)
  real(rt), intent(in   ) :: new_state(ns_lo(1):ns_hi(1),ns_lo(2):ns_hi(2),ns_lo(3):ns_hi(3),NVAR)
  real(rt), intent(inout) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NSRC)
  real(rt), intent(in   ) :: problo(3), dx(3)
  real(rt), intent(in   ), value :: time, dt

  integer :: i, j, k
  real(rt) :: rho, temp, T6, T613, X_CNO, X_1, g14, eps_CNO

  src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.e0_rt

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           rho = new_state(i,j,k,URHO)

           T6 = new_state(i,j,k,UTEMP)/1.0e6_rt
           T613 = T6**THIRD

           ! CNO abundance
           X_CNO = (new_state(i,j,k,UFS-1+ic12) + &
                    new_state(i,j,k,UFS-1+in14) + &
                    new_state(i,j,k,UFS-1+io16))/rho

           ! H abundance
           X_1 = new_state(i,j,k,UFS-1+ih1)/rho

           ! CNO heating from Kippenhahn & Weigert, Eq. 18.65
           g14 = 1.0_rt + 2.7e-3_rt*T613 - 7.78e-3_rt*T613**2 - 1.49e-4_rt*T6
           eps_CNO = 8.67e27_rt * g14 * X_CNO * X_1 * rho * exp(-152.28_rt/T613) / T613**2

           ! source terms
           src(i,j,k,UEDEN) = rho*eps_CNO
           src(i,j,k,UEINT) = rho*eps_CNO

        end do
     end do
  end do

end subroutine ca_ext_src
