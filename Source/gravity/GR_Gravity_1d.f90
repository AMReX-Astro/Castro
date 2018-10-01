! :::
! ::: ----------------------------------------------------------------
! :::

subroutine ca_compute_avgpres (lo,hi,dx,dr,&
     var,r_lo,r_hi, &
     radial_pres,problo, &
     n1d,drdxfac,level) bind(C,name='ca_compute_avgpres')
  use prob_params_module, only : center, coord_type, Symmetry, physbc_lo
  use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, UFS, UFX
  use eos_module, only : eos
  use eos_type_module, only : eos_input_re, eos_t
  use network, only : nspec, naux
  use amrex_constants_module, only: HALF, FOUR3RD, M_PI
  use amrex_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer , intent(in   ) :: lo(3),hi(3)
  real(rt), intent(in   ) :: dx(3)
  real(rt), value, intent(in   ) :: dr
  real(rt), intent(in   ) :: problo(3)

  integer , value, intent(in   ) :: n1d,drdxfac,level
  real(rt), intent(inout) :: radial_pres(0:n1d-1)

  integer , intent(in   ) :: r_lo(3),r_hi(3)
  real(rt), intent(in   ) :: var(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),NVAR)

  integer          :: i,j,k,n,index
  integer          :: ii
  real(rt)         :: r
  real(rt)         :: fac,dx_frac,vol
  real(rt)         :: lo_i,rlo,rhi

  type (eos_t) :: eos_state

#ifndef AMREX_USE_CUDA
  if (physbc_lo(1) .ne. Symmetry) then
     call amrex_error("Error: GR_Gravity_1d.f90 :: 1D gravity assumes symmetric lower boundary.")
  endif

  if (coord_type .ne. 2) then
     call amrex_error("Error: GR_Gravity_1d.f90 :: 1D gravity assumes spherical coordinates.")
  endif
#endif

  fac  = dble(drdxfac)
  dx_frac = dx(1) / fac

  j = lo(2)
  k = lo(3)

  do i = lo(1), hi(1)

     r = abs(problo(1) + (dble(i) + HALF) * dx(1) - center(1))

     index = int(r / dr)

     if (index .gt. n1d-1) then
#ifndef AMREX_USE_CUDA
        if (level .eq. 0) then
           print *,'   '
           print *,'>>> Error: GR_Gravity_1d::ca_compute_avgpres ',i
           print *,'>>> ... index too big: ', index,' > ',n1d-1
           print *,'>>> ... at i     : ',i
           print *,'    '
           call amrex_error("Error:: GR_Gravity_1d.f90 :: ca_compute_avgpres")
        end if
#endif

     else

        eos_state % rho = var(i,j,k,URHO)
        eos_state % e   = var(i,j,k,UEINT) / eos_state % rho
        eos_state % T   = var(i,j,k,UTEMP)
        eos_state % xn  = var(i,j,k,UFS:UFS+nspec-1) / eos_state % rho
        eos_state % aux = var(i,j,k,UFX:UFX+naux-1) / eos_state % rho

        ! Compute pressure from the EOS

        call eos(eos_input_re, eos_state)

        ! Note that we assume we are in spherical coordinates in 1d or we wouldn't be
        ! doing monopole gravity.

        lo_i = problo(1) + dble(i) * dx(1) - center(1)

        do ii = 0, drdxfac-1

           r   = abs(lo_i + (dble(ii  ) + HALF) * dx_frac)
           rlo = abs(lo_i +  dble(ii  )         * dx_frac)
           rhi = abs(lo_i +  dble(ii+1)         * dx_frac)

           vol = FOUR3RD * M_PI * (rhi**3 - rlo**3)

           index = int(r / dr)

           if (index .le. n1d-1) then
              radial_pres(index) = radial_pres(index) + vol * eos_state % P
           end if

        end do

     end if
  enddo

end subroutine ca_compute_avgpres
