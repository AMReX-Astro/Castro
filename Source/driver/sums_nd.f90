module castro_sums_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains
   
  subroutine ca_summass(lo,hi,rho,r_lo,r_hi,dx,&
                        vol,v_lo,v_hi,mass) bind(C, name="ca_summass")

    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: r_lo(3), r_hi(3)
    integer          :: v_lo(3), v_hi(3)
    real(rt)         :: mass, dx(3)
    real(rt)         :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt)         :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer          :: i, j, k

    mass = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mass = mass + rho(i,j,k) * vol(i,j,k)
          enddo
       enddo
    enddo

  end subroutine ca_summass



  subroutine ca_sumsquared(lo,hi,rho,r_lo,r_hi,dx,&
                           vol,v_lo,v_hi,mass) bind(C, name="ca_sumsquared")

    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: r_lo(3), r_hi(3)
    integer          :: v_lo(3), v_hi(3)
    real(rt)         :: mass, dx(3)
    real(rt)         :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt)         :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer          :: i, j, k

    mass = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mass = mass + rho(i,j,k)*rho(i,j,k)*vol(i,j,k)
          enddo
       enddo
    enddo

  end subroutine ca_sumsquared

  

  subroutine ca_sumlocmass(lo,hi,rho,r_lo,r_hi,dx,&
                           vol,v_lo,v_hi,mass,idir) bind(C, name="ca_sumlocmass")

    use prob_params_module, only: problo, center, probhi, dim, physbc_lo, physbc_hi, Symmetry
    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: idir
    integer          :: lo(3), hi(3)
    integer          :: r_lo(3), r_hi(3)
    integer          :: v_lo(3), v_hi(3)
    real(rt)         :: mass, dx(3)
    real(rt)         :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt)         :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer          :: i, j, k
    real(rt)         :: x, y, z
    real(rt)         :: symlo, symhi

    mass = ZERO

    symlo = ZERO
    symhi = ZERO

    if (idir .eq. 0) then
       do i = lo(1), hi(1)
          x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
          if (physbc_lo(1) .eq. Symmetry) then
             symlo = problo(1) - x
          endif
          if (physbc_hi(1) .eq. Symmetry) then
             symhi = x - probhi(1)
          endif
          x = x + symlo + symhi
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                mass = mass + rho(i,j,k) * vol(i,j,k) * x
             enddo
          enddo
       enddo
    else if (idir .eq. 1 .and. dim .ge. 2) then
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
          if (physbc_lo(2) .eq. Symmetry) then
             symlo = problo(2) - y
          endif
          if (physbc_hi(2) .eq. Symmetry) then
             symhi = y - probhi(2)
          endif
          y = y + symlo + symhi
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)
                mass = mass + rho(i,j,k) * vol(i,j,k) * y
             enddo
          enddo
       enddo
    else if (dim .eq. 3) then
       do k = lo(3), hi(3)
          z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
          if (physbc_lo(3) .eq. Symmetry) then
             symlo = problo(3) - z
          endif
          if (physbc_hi(3) .eq. Symmetry) then
             symhi = z - probhi(3)
          endif
          z = z + symlo + symhi
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mass = mass + rho(i,j,k) * vol(i,j,k) * z
             enddo
          enddo
       enddo
    end if

  end subroutine ca_sumlocmass



  subroutine ca_sumlocmass2d(lo,hi,rho,r_lo,r_hi,dx,&
                             vol,v_lo,v_hi,mass,idir1,idir2) bind(C, name="ca_sumlocmass2d")

    use prob_params_module, only: problo, center, probhi, dim, physbc_lo, physbc_hi, Symmetry
    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: idir1, idir2
    integer          :: lo(3), hi(3)
    integer          :: r_lo(3), r_hi(3)
    integer          :: v_lo(3), v_hi(3)
    real(rt)         :: mass, dx(3)
    real(rt)         :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt)         :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer          :: i, j, k
    real(rt)         :: x, y, z
    real(rt)         :: symlo1, symhi1, symlo2, symhi2

    mass = ZERO

    symlo1 = ZERO
    symlo2 = ZERO
    symhi1 = ZERO
    symhi2 = ZERO

    if (idir1 .eq. 0) then
       do i = lo(1), hi(1)
          x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
          if (physbc_lo(1) .eq. Symmetry) then
             symlo1 = problo(1) - x
          endif
          if (physbc_hi(1) .eq. Symmetry) then
             symhi1 = x - probhi(1)
          endif
          x = x + symlo1 + symhi1
          if (idir2 .eq. 0) then
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   mass = mass + rho(i,j,k) * vol(i,j,k) * x * x
                enddo
             enddo
          elseif (idir2 .eq. 1 .and. dim .ge. 2) then
             do j = lo(2), hi(2)
                y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
                if (physbc_lo(2) .eq. Symmetry) then
                   symlo2 = problo(2) - y
                endif
                if (physbc_hi(2) .eq. Symmetry) then
                   symhi2 = y - probhi(2)
                endif
                y = y + symlo2 + symhi2
                do k = lo(3), hi(3)
                   mass = mass + rho(i,j,k) * vol(i,j,k) * x * y
                enddo
             enddo
          elseif (dim .eq. 3) then
             do k = lo(3), hi(3)
                z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
                if (physbc_lo(3) .eq. Symmetry) then
                   symlo2 = problo(3) - z
                endif
                if (physbc_hi(3) .eq. Symmetry) then
                   symhi2 = z - probhi(3)
                endif
                z = z + symlo2 + symhi2
                do j = lo(2), hi(2)
                   mass = mass + rho(i,j,k) * vol(i,j,k) * x * z
                enddo
             enddo
          endif
       enddo
    else if (idir1 .eq. 1 .and. dim .ge. 2) then
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
          if (physbc_lo(2) .eq. Symmetry) then
             symlo1 = problo(2) - y
          endif
          if (physbc_hi(2) .eq. Symmetry) then
             symhi1 = y - probhi(2)
          endif
          y = y + symlo1 + symhi1
          if (idir2 .eq. 0) then
             do i = lo(1), hi(1)
                x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
                if (physbc_lo(1) .eq. Symmetry) then
                   symlo2 = problo(1) - x
                endif
                if (physbc_hi(1) .eq. Symmetry) then
                   symhi2 = x - probhi(1)
                endif
                x = x + symlo2 + symhi2
                do k = lo(3), hi(3)
                   mass = mass + rho(i,j,k) * vol(i,j,k) * y * x
                enddo
             enddo
          elseif (idir2 .eq. 1) then
             do i = lo(1), hi(1)
                do k = lo(3), hi(3)
                   mass = mass + rho(i,j,k) * vol(i,j,k) * y * y
                enddo
             enddo
          elseif (dim .eq. 3) then
             do k = lo(3), hi(3)
                z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
                if (physbc_lo(3) .eq. Symmetry) then
                   symlo2 = problo(3) - z
                endif
                if (physbc_hi(3) .eq. Symmetry) then
                   symhi2 = z - probhi(3)
                endif
                z = z + symlo2 + symhi2                   
                do i = lo(1), hi(1)
                   mass = mass + rho(i,j,k) * vol(i,j,k) * y * z
                enddo
             enddo
          endif
       enddo
    elseif (dim .eq. 3) then
       do k = lo(3), hi(3)
          z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
          if (physbc_lo(3) .eq. Symmetry) then
             symlo1 = problo(3) - z
          endif
          if (physbc_hi(3) .eq. Symmetry) then
             symhi1 = z - probhi(3)
          endif
          z = z + symlo1 + symhi1
          if (idir2 .eq. 0) then
             do i = lo(1), hi(1)
                x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
                if (physbc_lo(1) .eq. Symmetry) then
                   symlo2 = problo(1) - x
                endif
                if (physbc_hi(1) .eq. Symmetry) then
                   symhi2 = x - probhi(1)
                endif
                x = x + symlo2 + symhi2
                do j = lo(2), hi(2)
                   mass = mass + rho(i,j,k) * vol(i,j,k) * z * x
                enddo
             enddo
          elseif (idir2 .eq. 1) then
             do j = lo(2), hi(2)
                y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
                if (physbc_lo(2) .eq. Symmetry) then
                   symlo2 = problo(2) - y
                endif
                if (physbc_hi(2) .eq. Symmetry) then
                   symhi2 = y - probhi(2)
                endif
                y = y + symlo2 + symhi2
                do i = lo(1), hi(1)
                   mass = mass + rho(i,j,k) * vol(i,j,k) * z * y
                enddo
             enddo
          else
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   mass = mass + rho(i,j,k) * vol(i,j,k) * z * z
                enddo
             enddo
          endif
       enddo
    end if

  end subroutine ca_sumlocmass2d



  subroutine ca_sumlocsquaredmass(lo,hi,rho,r_lo,r_hi,dx,&
                                  vol,v_lo,v_hi,mass,idir) bind(C, name="ca_sumlocsquaredmass")

    use prob_params_module, only: problo, center, dim
    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: idir
    integer          :: lo(3), hi(3)
    integer          :: r_lo(3), r_hi(3)
    integer          :: v_lo(3), v_hi(3)
    real(rt)         :: mass, dx(3)
    real(rt)         :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt)         :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer          :: i, j, k
    real(rt)         :: x, y, z

    mass = ZERO

    if (idir .eq. 0) then
       do i = lo(1), hi(1)
          x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                mass = mass + rho(i,j,k) * vol(i,j,k) * (x**2)
             enddo
          enddo
       enddo
    else if (idir .eq. 1 .and. dim .ge. 2) then
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)
                mass = mass + rho(i,j,k) * vol(i,j,k) * (y**2)
             enddo
          enddo
       enddo
    else if (dim .eq. 3) then
       do k = lo(3), hi(3)
          z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mass = mass + rho(i,j,k) * vol(i,j,k) * (z**2)
             enddo
          enddo
       enddo
    end if

  end subroutine ca_sumlocsquaredmass


  
  subroutine ca_sumproduct(lo,hi,f1,f1_lo,f1_hi,f2,f2_lo,f2_hi,dx,&
                           vol,v_lo,v_hi,product) bind(C, name="ca_sumproduct")

    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: f1_lo(3), f1_hi(3)
    integer          :: f2_lo(3), f2_hi(3)
    integer          :: v_lo(3), v_hi(3)
    real(rt)         :: product
    real(rt)         :: dx(3)
    real(rt)         :: f1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3))
    real(rt)         :: f2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3))
    real(rt)         :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer          :: i, j, k

    product = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             product = product + f1(i,j,k) * f2(i,j,k) * vol(i,j,k)
          enddo
       enddo
    enddo

  end subroutine ca_sumproduct

end module castro_sums_module
