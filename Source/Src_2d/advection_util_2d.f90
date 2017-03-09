module advection_util_2d_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public normalize_species_fluxes, divu

contains

  subroutine normalize_species_fluxes(  &
                    flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                    flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                    lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module
    
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(2),hi(2)
    integer          :: flux1_l1,flux1_l2,flux1_h1,flux1_h2
    integer          :: flux2_l1,flux2_l2,flux2_h1,flux2_h2
    real(rt)         :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
    real(rt)         :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
    
    ! Local variables
    integer          :: i,j,n
    real(rt)         :: sum,fac
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + flux1(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = flux1(i,j,URHO) / sum
          else
             fac = ONE
          end if
          do n = UFS, UFS+nspec-1
             flux1(i,j,n) = flux1(i,j,n) * fac
          end do
       end do
    end do
    do j = lo(2),hi(2)+1
       do i = lo(1),hi(1)
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + flux2(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = flux2(i,j,URHO) / sum
          else
             fac = ONE
          end if
          do n = UFS, UFS+nspec-1
             flux2(i,j,n) = flux2(i,j,n) * fac
          end do
       end do
    end do
    
  end subroutine normalize_species_fluxes

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine normalize_new_species(u,u_l1,u_l2,u_h1,u_h2,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module
    
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(2), hi(2)
    integer          :: u_l1,u_l2,u_h1,u_h2
    real(rt)         :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
    
    ! Local variables
    integer          :: i,j,n
    real(rt)         :: fac,sum
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + u(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = u(i,j,URHO) / sum
          else
             fac = ONE
          end if
          do n = UFS, UFS+nspec-1
             u(i,j,n) = u(i,j,n) * fac
          end do
       end do
    end do
    
  end subroutine normalize_new_species

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2,dx, &
                  div,div_l1,div_l2,div_h1,div_h2)

    use prob_params_module, only : coord_type
    use meth_params_module, only : QU, QV
    use bl_constants_module
    
    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    integer          :: lo(2),hi(2)
    integer          :: q_l1,q_l2,q_h1,q_h2
    integer          :: div_l1,div_l2,div_h1,div_h2
    real(rt)         :: q(q_l1:q_h1,q_l2:q_h2,*)
    real(rt)         :: div(div_l1:div_h1,div_l2:div_h2)
    real(rt)         :: dx(2)
    
    integer          :: i, j
    real(rt)         :: rl, rr, rc, ul, ur
    real(rt)         :: vb, vt
    real(rt)         :: ux,vy
    
    if (coord_type .eq. 0) then
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1
             ux = HALF*(q(i,j,QU)-q(i-1,j,QU)+q(i,j-1,QU)-q(i-1,j-1,QU))/dx(1)
             vy = HALF*(q(i,j,QV)-q(i,j-1,QV)+q(i-1,j,QV)-q(i-1,j-1,QV))/dx(2)
             div(i,j) = ux + vy
          enddo
       enddo
    else
       do i=lo(1),hi(1)+1
          
          if (i.eq.0) then
             
             div(i,lo(2):hi(2)+1) = ZERO
             
          else 

             rl = (dble(i)-HALF) * dx(1)
             rr = (dble(i)+HALF) * dx(1)
             rc = (dble(i)     ) * dx(1)
             
             do j=lo(2),hi(2)+1
                ! These are transverse averages in the y-direction
                ul = HALF * (q(i-1,j,QU)+q(i-1,j-1,QU))
                ur = HALF * (q(i  ,j,QU)+q(i  ,j-1,QU))
                
                ! Take 1/r d/dr(r*u)
                div(i,j) = (rr*ur - rl*ul) / dx(1) / rc
                
                ! These are transverse averages in the x-direction
                vb = HALF * (q(i,j-1,QV)+q(i-1,j-1,QV))
                vt = HALF * (q(i,j  ,QV)+q(i-1,j  ,QV))
                
                div(i,j) = div(i,j) + (vt - vb) / dx(2)
             enddo
             
          end if
       enddo
    end if

  end subroutine divu

end module advection_util_2d_module
