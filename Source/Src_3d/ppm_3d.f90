module ppm_module

  implicit none

  private

  public ppm

contains
  !
  ! characteristics based on u
  !
  subroutine ppm(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                 u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                 flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                 Ip,Im, &
                 ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)

    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals

    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2

    double precision     s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    double precision     u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    double precision  cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn( f_l1: f_h1, f_l2: f_h2, f_l3:f_h3)
    double precision Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    double precision dx,dy,dz,dt
    integer          k3d,kc
   
    if (ppm_type .eq. 1) then

        call ppm_type1(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)

    else if (ppm_type .eq. 2) then

        call ppm_type2(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)

    end if

  end subroutine ppm

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::
  
  subroutine ppm_type1(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)

    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals

    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2

    double precision     s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    double precision     u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    double precision  cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn( f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)

    double precision Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    double precision dx,dy,dz,dt
    integer          k3d,kc

    ! local
    integer i,j,k

    double precision dsl, dsr, dsc
    double precision sigma, s6

    ! s_{\ib,+}, s_{\ib,-}
    double precision, allocatable :: sp(:,:)
    double precision, allocatable :: sm(:,:)

    ! \delta s_{\ib}^{vL}
    double precision, allocatable :: dsvl(:,:)
    double precision, allocatable :: dsvlm(:,:)
    double precision, allocatable :: dsvlp(:,:)

    ! s_{i+\half}^{H.O.}
    double precision, allocatable :: sedge(:,:)
    double precision, allocatable :: sedgez(:,:,:)

    ! cell-centered indexing
    allocate(sp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(sm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    if (ppm_type .ne. 1) &
         call bl_error("Should have ppm_type = 1 in ppm_type1")

    if (s_l1 .gt. ilo1-3 .or. s_l2 .gt. ilo2-3) then
         print *,'Low bounds of array: ',s_l1, s_l2
         print *,'Low bounds of  loop: ',ilo1 , ilo2
         call bl_error("Need more ghost cells on array in ppm_type1")
    end if

    if (s_h1 .lt. ihi1+3 .or. s_h2 .lt. ihi2+3) then
         print *,'Hi  bounds of array: ',s_h1, s_h2
         print *,'Hi  bounds of  loop: ',ihi1 , ihi2
         call bl_error("Need more ghost cells on array in ppm_type1")
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(ilo1-2:ihi1+2,ilo2-1:ihi2+1))

    ! edge-centered indexing for x-faces -- ppm_type = 1 only
    allocate(sedge(ilo1-1:ihi1+2,ilo2-1:ihi2+1))

    ! compute s at x-edges

    ! compute van Leer slopes in x-direction
    dsvl = 0.d0
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+2
          dsc = 0.5d0 * (s(i+1,j,k3d) - s(i-1,j,k3d))
          dsl = 2.d0  * (s(i  ,j,k3d) - s(i-1,j,k3d))
          dsr = 2.d0  * (s(i+1,j,k3d) - s(i  ,j,k3d))
          if (dsl*dsr .gt. 0.d0) &
               dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do
    end do

    ! interpolate s to x-edges
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+2
          sedge(i,j) = 0.5d0*(s(i,j,k3d)+s(i-1,j,k3d)) &
               - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i-1,j))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i-1,j,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i-1,j,k3d)))
       end do
    end do

    !$OMP PARALLEL DO PRIVATE(i,j,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! copy sedge into sp and sm
          sp(i,j) = sedge(i+1,j)
          sm(i,j) = sedge(i  ,j)

          if (ppm_flatten_before_integrals == 1) then
             ! flatten the parabola BEFORE doing the other                     
             ! monotonization -- this is the method that Flash does       
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          ! modify using quadratic limiters -- note this version of the limiting comes
          ! from Colella and Sekora (2008), not the original PPM paper.
          if ((sp(i,j)-s(i,j,k3d))*(s(i,j,k3d)-sm(i,j)) .le. 0.d0) then
             sp(i,j) = s(i,j,k3d)
             sm(i,j) = s(i,j,k3d)

          else if (abs(sp(i,j)-s(i,j,k3d)) .ge. 2.d0*abs(sm(i,j)-s(i,j,k3d))) then
          !else if (-(sp(i,j)-sm(i,j))**2/6.0d0 > &
          !     (sp(i,j) - sm(i,j))*(s(i,j,k3d) - 0.5d0*(sm(i,j) + sp(i,j)))) then
             sp(i,j) = 3.d0*s(i,j,k3d) - 2.d0*sm(i,j)

          else if (abs(sm(i,j)-s(i,j,k3d)) .ge. 2.d0*abs(sp(i,j)-s(i,j,k3d))) then
          !else if ((sp(i,j)-sm(i,j))*(s(i,j,k3d) - 0.5d0*(sm(i,j) + sp(i,j))) > &
          !     (sp(i,j) - sm(i,j))**2/6.0d0) then
             sm(i,j) = 3.d0*s(i,j,k3d) - 2.d0*sp(i,j)
          end if

          if (ppm_flatten_before_integrals == 2) then
             ! flatten the parabola AFTER doing the monotonization --
             ! this is the method that Miller & Colella do
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          ! compute x-component of Ip and Im
          s6 = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))

          ! Ip/m is the integral under the parabola for the extent
          ! that a wave can travel over a timestep
          !
          ! Ip integrates to the right edge of a cell
          ! Im integrates to the left edge of a cell

          ! u-c wave
          sigma = abs(u(i,j,k3d,1)-cspd(i,j,k3d))*dt/dx

          if (u(i,j,k3d,1)-cspd(i,j,k3d) <= 0.0d0) then
             Ip(i,j,kc,1,1) = sp(i,j)
          else
             Ip(i,j,kc,1,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          if (u(i,j,k3d,1)-cspd(i,j,k3d) >= 0.0d0) then
             Im(i,j,kc,1,1) = sm(i,j) 
          else
             Im(i,j,kc,1,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          ! u wave
          sigma = abs(u(i,j,k3d,1))*dt/dx

          if (u(i,j,k3d,1) <= 0.0d0) then
             Ip(i,j,kc,1,2) = sp(i,j) 
          else
             Ip(i,j,kc,1,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif
             
          if (u(i,j,k3d,1) >= 0.0d0) then
             Im(i,j,kc,1,2) = sm(i,j) 
          else
             Im(i,j,kc,1,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          ! u+c wave
          sigma = abs(u(i,j,k3d,1)+cspd(i,j,k3d))*dt/dx

          if (u(i,j,k3d,1)+cspd(i,j,k3d) <= 0.0d0) then
             Ip(i,j,kc,1,3) = sp(i,j) 
          else
             Ip(i,j,kc,1,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          if (u(i,j,k3d,1)+cspd(i,j,k3d) >= 0.0d0) then
             Im(i,j,kc,1,3) = sm(i,j) 
          else
             Im(i,j,kc,1,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-2:ihi2+2))

    ! edge-centered indexing for y-faces
    allocate(sedge(ilo1-1:ihi1+1,ilo2-1:ihi2+2))

    ! compute s at y-edges

    ! compute van Leer slopes in y-direction
    dsvl = 0.d0
    do j=ilo2-2,ihi2+2
       do i=ilo1-1,ihi1+1
          dsc = 0.5d0 * (s(i,j+1,k3d) - s(i,j-1,k3d))
          dsl = 2.d0  * (s(i,j  ,k3d) - s(i,j-1,k3d))
          dsr = 2.d0  * (s(i,j+1,k3d) - s(i,j  ,k3d))
          if (dsl*dsr .gt. 0.d0) &
               dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do
    end do

    ! interpolate s to y-edges
    do j=ilo2-1,ihi2+2
       do i=ilo1-1,ihi1+1
          sedge(i,j) = 0.5d0*(s(i,j,k3d)+s(i,j-1,k3d)) &
               - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i,j-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i,j-1,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i,j-1,k3d)))
       end do
    end do

    !$OMP PARALLEL DO PRIVATE(i,j,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! copy sedge into sp and sm
          sp(i,j) = sedge(i,j+1)
          sm(i,j) = sedge(i,j  )

          if (ppm_flatten_before_integrals == 1) then
             ! flatten the parabola BEFORE doing the other                     
             ! monotonization -- this is the method that Flash does       
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          ! modify using quadratic limiters
          if ((sp(i,j)-s(i,j,k3d))*(s(i,j,k3d)-sm(i,j)) .le. 0.d0) then
             sp(i,j) = s(i,j,k3d)
             sm(i,j) = s(i,j,k3d)

          else if (abs(sp(i,j)-s(i,j,k3d)) .ge. 2.d0*abs(sm(i,j)-s(i,j,k3d))) then
          !else if (-(sp(i,j)-sm(i,j))**2/6.0d0 > &
          !     (sp(i,j) - sm(i,j))*(s(i,j,k3d) - 0.5d0*(sm(i,j) + sp(i,j)))) then
             sp(i,j) = 3.d0*s(i,j,k3d) - 2.d0*sm(i,j)

          else if (abs(sm(i,j)-s(i,j,k3d)) .ge. 2.d0*abs(sp(i,j)-s(i,j,k3d))) then
          !else if ((sp(i,j)-sm(i,j))*(s(i,j,k3d) - 0.5d0*(sm(i,j) + sp(i,j))) > &
          !     (sp(i,j) - sm(i,j))**2/6.0d0) then
             sm(i,j) = 3.d0*s(i,j,k3d) - 2.d0*sp(i,j)
          end if

          if (ppm_flatten_before_integrals == 2) then
             ! flatten the parabola AFTER doing the monotonization --
             ! this is the method that Miller & Colella do
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          ! compute y-component of Ip and Im
          s6 = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))

          ! v-c wave
          sigma = abs(u(i,j,k3d,2)-cspd(i,j,k3d))*dt/dy

          if (u(i,j,k3d,2)-cspd(i,j,k3d) <= 0.0d0) then
             Ip(i,j,kc,2,1) = sp(i,j)
          else
             Ip(i,j,kc,2,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          if (u(i,j,k3d,2)-cspd(i,j,k3d) >= 0.0d0) then
             Im(i,j,kc,2,1) = sm(i,j) 
          else
             Im(i,j,kc,2,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          ! v wave
          sigma = abs(u(i,j,k3d,2))*dt/dy

          if (u(i,j,k3d,2) <= 0.0d0) then
             Ip(i,j,kc,2,2) = sp(i,j) 
          else
             Ip(i,j,kc,2,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          if (u(i,j,k3d,2) >= 0.0d0) then
             Im(i,j,kc,2,2) = sm(i,j) 
          else
             Im(i,j,kc,2,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          ! v+c wave
          sigma = abs(u(i,j,k3d,2)+cspd(i,j,k3d))*dt/dy

          if (u(i,j,k3d,2)+cspd(i,j,k3d) <= 0.0d0) then
             Ip(i,j,kc,2,3) = sp(i,j) 
          else
             Ip(i,j,kc,2,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          if (u(i,j,k3d,2)+cspd(i,j,k3d) >= 0.0d0) then
             Im(i,j,kc,2,3) = sm(i,j) 
          else
             Im(i,j,kc,2,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(dsvl,sedge)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    allocate(sedgez(ilo1-1:ihi1+1,ilo2-2:ihi2+3,k3d-1:k3d+2))

    ! compute s at z-edges

    ! compute van Leer slopes in z-direction
    dsvl  = 0.d0
    dsvlm = 0.d0
    dsvlp = 0.d0

    !$OMP PARALLEL DO PRIVATE(i,j,k,dsc,dsl,dsr,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! compute on slab below
          k = k3d-1
          dsc = 0.5d0 * (s(i,j,k+1) - s(i,j,k-1))
          dsl = 2.0d0 * (s(i,j,k  ) - s(i,j,k-1))
          dsr = 2.0d0 * (s(i,j,k+1) - s(i,j,k  ))
          if (dsl*dsr .gt. 0.d0) &
               dsvlm(i,j) = sign(1.0d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

          ! compute on slab above
          k = k3d+1
          dsc = 0.5d0 * (s(i,j,k+1) - s(i,j,k-1))
          dsl = 2.0d0 * (s(i,j,k  ) - s(i,j,k-1))
          dsr = 2.0d0 * (s(i,j,k+1) - s(i,j,k  ))
          if (dsl*dsr .gt. 0.d0) &
               dsvlp(i,j) = sign(1.0d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

          ! compute on current slab
          k = k3d
          dsc = 0.5d0 * (s(i,j,k+1) - s(i,j,k-1))
          dsl = 2.0d0 * (s(i,j,k  ) - s(i,j,k-1))
          dsr = 2.0d0 * (s(i,j,k+1) - s(i,j,k  ))
          if (dsl*dsr .gt. 0.d0) &
               dsvl(i,j) = sign(1.0d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

          ! interpolate to lo face
          k = k3d
          sm(i,j) = 0.5d0*(s(i,j,k)+s(i,j,k-1)) - (1.0d0/6.0d0)*(dsvl(i,j)-dsvlm(i,j))
          ! make sure sedge lies in between adjacent cell-centered values
          sm(i,j) = max(sm(i,j),min(s(i,j,k),s(i,j,k-1)))
          sm(i,j) = min(sm(i,j),max(s(i,j,k),s(i,j,k-1)))

          ! interpolate to hi face
          k = k3d+1
          sp(i,j) = 0.5d0*(s(i,j,k)+s(i,j,k-1)) - (1.0d0/6.0d0)*(dsvlp(i,j)-dsvl(i,j))

          ! make sure sedge lies in between adjacent cell-centered values
          sp(i,j) = max(sp(i,j),min(s(i,j,k),s(i,j,k-1)))
          sp(i,j) = min(sp(i,j),max(s(i,j,k),s(i,j,k-1)))

          if (ppm_flatten_before_integrals == 1) then
             ! flatten the parabola BEFORE doing the other                     
             ! monotonization -- this is the method that Flash does       
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
          endif


          ! modify using quadratic limiters
          if ((sp(i,j)-s(i,j,k3d))*(s(i,j,k3d)-sm(i,j)) .le. 0.0d0) then
             sp(i,j) = s(i,j,k3d)
             sm(i,j) = s(i,j,k3d)

          else if (abs(sp(i,j)-s(i,j,k3d)) .ge. 2.0d0*abs(sm(i,j)-s(i,j,k3d))) then
          !else if (-(sp(i,j)-sm(i,j))**2/6.0d0 > &
          !     (sp(i,j) - sm(i,j))*(s(i,j,k3d) - 0.5d0*(sm(i,j) + sp(i,j)))) then
             sp(i,j) = 3.0d0*s(i,j,k3d) - 2.0d0*sm(i,j)

          else if (abs(sm(i,j)-s(i,j,k3d)) .ge. 2.0d0*abs(sp(i,j)-s(i,j,k3d))) then
          !else if ((sp(i,j)-sm(i,j))*(s(i,j,k3d) - 0.5d0*(sm(i,j) + sp(i,j))) > &
          !     (sp(i,j) - sm(i,j))**2/6.0d0) then
             sm(i,j) = 3.0d0*s(i,j,k3d) - 2.0d0*sp(i,j)
          end if

          if (ppm_flatten_before_integrals == 2) then
             ! flatten the parabola AFTER doing the monotonization --
             ! this is the method that Miller & Colella do
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (1.d0-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          ! compute z-component of Ip and Im
          s6 = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))

          ! w-c wave
          sigma = abs(u(i,j,k3d,3)-cspd(i,j,k3d))*dt/dz

          if (u(i,j,k3d,3)-cspd(i,j,k3d) <= 0.0d0) then
             Ip(i,j,kc,3,1) = sp(i,j) 
          else
             Ip(i,j,kc,3,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          if (u(i,j,k3d,3)-cspd(i,j,k3d) >= 0.0d0) then
             Im(i,j,kc,3,1) = sm(i,j) 
          else
             Im(i,j,kc,3,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          ! w wave
          sigma = abs(u(i,j,k3d,3))*dt/dz

          if (u(i,j,k3d,3) <= 0.0d0) then
             Ip(i,j,kc,3,2) = sp(i,j) 
          else
             Ip(i,j,kc,3,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          if (u(i,j,k3d,3) >= 0.0d0) then
             Im(i,j,kc,3,2) = sm(i,j) 
          else
             Im(i,j,kc,3,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          ! w+c wave
          sigma = abs(u(i,j,k3d,3)+cspd(i,j,k3d))*dt/dz

          if (u(i,j,k3d,3)+cspd(i,j,k3d) <= 0.0d0) then
             Ip(i,j,kc,3,3) = sp(i,j) 
          else
             Ip(i,j,kc,3,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

          if (u(i,j,k3d,3)+cspd(i,j,k3d) >= 0.0d0) then
             Im(i,j,kc,3,3) = sm(i,j) 
          else
             Im(i,j,kc,3,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          endif

       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(dsvl,dsvlm,dsvlp,sp,sm,sedgez)

  end subroutine ppm_type1

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine ppm_type2(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)

    use meth_params_module, only : ppm_type

    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2
    double precision s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    double precision u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,1:3)
    double precision cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)
    double precision Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision dx,dy,dz,dt
    integer          k3d,kc

    ! local
    integer i,j,k
    logical extremum, bigp, bigm

    double precision D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    double precision sgn, sigma, s6
    double precision dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
    double precision dachkm, dachkp
    double precision amax, delam, delap

    ! s_{\ib,+}, s_{\ib,-}
    double precision, allocatable :: sp(:,:)
    double precision, allocatable :: sm(:,:)

    ! \delta s_{\ib}^{vL}
    double precision, allocatable :: dsvl(:,:)
    double precision, allocatable :: dsvlm(:,:)
    double precision, allocatable :: dsvlp(:,:)

    ! s_{i+\half}^{H.O.}
    double precision, allocatable :: sedge(:,:)
    double precision, allocatable :: sedgez(:,:,:)

    ! constant used in Colella 2008
    double precision, parameter :: C = 1.25d0

    ! cell-centered indexing
    allocate(sp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(sm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    if (ppm_type .ne. 2) &
         call bl_error("Should have ppm_type = 2 in ppm_type2")

    if (s_l1 .gt. ilo1-3 .or. s_l2 .gt. ilo2-3) then
         print *,'Low bounds of array: ',s_l1, s_l2
         print *,'Low bounds of  loop: ',ilo1 , ilo2
         call bl_error("Need more ghost cells on array in ppm_type2")
    end if

    if (s_h1 .lt. ihi1+3 .or. s_h2 .lt. ihi2+3) then
         print *,'Hi  bounds of array: ',s_h1, s_h2
         print *,'Hi  bounds of  loop: ',ihi1 , ihi2
         call bl_error("Need more ghost cells on array in ppm_type2")
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(ilo1-2:ihi1+2,ilo2-1:ihi2+1))

    ! edge-centered indexing for x-faces
    allocate(sedge(ilo1-2:ihi1+3,ilo2-1:ihi2+1))

    ! compute s at x-edges

    ! interpolate s to x-edges
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+3
          sedge(i,j) = (7.d0/12.d0)*(s(i-1,j,k3d)+s(i  ,j,k3d)) &
               - (1.d0/12.d0)*(s(i-2,j,k3d)+s(i+1,j,k3d))
          !
          ! limit sedge
          !
          if ((sedge(i,j)-s(i-1,j,k3d))*(s(i,j,k3d)-sedge(i,j)) .lt. 0.d0) then
             D2  = 3.d0*(s(i-1,j,k3d)-2.d0*sedge(i,j)+s(i,j,k3d))
             D2L = s(i-2,j,k3d)-2.d0*s(i-1,j,k3d)+s(i,j,k3d)
             D2R = s(i-1,j,k3d)-2.d0*s(i,j,k3d)+s(i+1,j,k3d)
             sgn = sign(1.d0,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
             sedge(i,j) = 0.5d0*(s(i-1,j,k3d)+s(i,j,k3d)) - (1.d0/6.d0)*D2LIM
          end if
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,alphap,alpham,bigp,bigm,extremum,dafacem) &
    !$OMP PRIVATE(dafacep,dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L) &
    !$OMP PRIVATE(D2R,D2C,sgn,D2LIM,amax,delam,delap,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedge(i+1,j)-s(i,j,k3d)
          alpham   = sedge(i  ,j)-s(i,j,k3d)
          bigp     = abs(alphap).gt.2.d0*abs(alpham)
          bigm     = abs(alpham).gt.2.d0*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. 0.d0) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedge(i,j) - sedge(i-1,j)
             dafacep   = sedge(i+2,j) - sedge(i+1,j)
             dabarm    = s(i,j,k3d) - s(i-1,j,k3d)
             dabarp    = s(i+1,j,k3d) - s(i,j,k3d)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2     = 6.d0*(alpham + alphap)
             D2L    = s(i-2,j,k3d)-2.d0*s(i-1,j,k3d)+s(i,j,k3d)
             D2R    = s(i,j,k3d)-2.d0*s(i+1,j,k3d)+s(i+2,j,k3d)
             D2C    = s(i-1,j,k3d)-2.d0*s(i,j,k3d)+s(i+1,j,k3d)
             sgn    = sign(1.d0,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(1.d0,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1,j,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -2.d0*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(1.d0,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1,j,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -2.d0*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k3d) + alpham
          sp(i,j) = s(i,j,k3d) + alphap
          !
          ! Compute x-component of Ip and Im.
          !
          s6    = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))

          sigma = abs(u(i,j,k3d,1)-cspd(i,j,k3d))*dt/dx
          Ip(i,j,kc,1,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,1,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)

          sigma = abs(u(i,j,k3d,1))*dt/dx
          Ip(i,j,kc,1,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,1,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)

          sigma = abs(u(i,j,k3d,1)+cspd(i,j,k3d))*dt/dx
          Ip(i,j,kc,1,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,1,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)

       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-2:ihi2+2))

    ! edge-centered indexing for y-faces
    allocate(sedge(ilo1-1:ihi1+1,ilo2-2:ihi2+3))

    ! compute s at y-edges

    ! interpolate s to y-edges
    do j=ilo2-2,ihi2+3
       do i=ilo1-1,ihi1+1
          sedge(i,j) = (7.d0/12.d0)*(s(i,j-1,k3d)+s(i,j,k3d)) &
               - (1.d0/12.d0)*(s(i,j-2,k3d)+s(i,j+1,k3d))
          !
          ! limit sedge
          !
          if ((sedge(i,j)-s(i,j-1,k3d))*(s(i,j,k3d)-sedge(i,j)) .lt. 0.d0) then
             D2  = 3.d0*(s(i,j-1,k3d)-2.d0*sedge(i,j)+s(i,j,k3d))
             D2L = s(i,j-2,k3d)-2.d0*s(i,j-1,k3d)+s(i,j,k3d)
             D2R = s(i,j-1,k3d)-2.d0*s(i,j,k3d)+s(i,j+1,k3d)
             sgn = sign(1.d0,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
             sedge(i,j) = 0.5d0*(s(i,j-1,k3d)+s(i,j,k3d)) - (1.d0/6.d0)*D2LIM
          end if
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep,dabarm,dabarp,dafacemin) &
    !$OMP PRIVATE(dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax,delam,delap,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedge(i,j+1)-s(i,j,k3d)
          alpham   = sedge(i,j  )-s(i,j,k3d)
          bigp     = abs(alphap).gt.2.d0*abs(alpham)
          bigm     = abs(alpham).gt.2.d0*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. 0.d0) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedge(i,j) - sedge(i,j-1)
             dafacep   = sedge(i,j+2) - sedge(i,j+1)
             dabarm    = s(i,j,k3d) - s(i,j-1,k3d)
             dabarp    = s(i,j+1,k3d) - s(i,j,k3d)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2     = 6.d0*(alpham + alphap)
             D2L    = s(i,j-2,k3d)-2.d0*s(i,j-1,k3d)+s(i,j,k3d)
             D2R    = s(i,j,k3d)-2.d0*s(i,j+1,k3d)+s(i,j+2,k3d)
             D2C    = s(i,j-1,k3d)-2.d0*s(i,j,k3d)+s(i,j+1,k3d)
             sgn    = sign(1.d0,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(1.d0,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i,j-1,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -2.d0*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(1.d0,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i,j+1,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -2.d0*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k3d) + alpham
          sp(i,j) = s(i,j,k3d) + alphap
          !
          ! Compute y-component of Ip and Im.
          !
          s6    = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))

          sigma = abs(u(i,j,k3d,2)-cspd(i,j,k3d))*dt/dy
          Ip(i,j,kc,2,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,2,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)

          sigma = abs(u(i,j,k3d,2))*dt/dy
          Ip(i,j,kc,2,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,2,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)

          sigma = abs(u(i,j,k3d,2)+cspd(i,j,k3d))*dt/dy
          Ip(i,j,kc,2,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,2,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(dsvl,sedge)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    allocate(sedgez(ilo1-1:ihi1+1,ilo2-2:ihi2+3,k3d-1:k3d+2))

    ! compute s at z-edges

    ! interpolate s to z-edges
    !$OMP PARALLEL DO PRIVATE(i,j,k,D2,D2L,D2R,sgn,D2LIM)
    do k=k3d-1,k3d+2
       do j=ilo2-1,ihi2+1
          do i=ilo1-1,ihi1+1
             sedgez(i,j,k) = (7.d0/12.d0)*(s(i,j,k-1)+s(i,j,k)) &
                  - (1.d0/12.d0)*(s(i,j,k-2)+s(i,j,k+1))
             !
             ! limit sedgez
             !
             if ((sedgez(i,j,k)-s(i,j,k-1))*(s(i,j,k)-sedgez(i,j,k)) .lt. 0.d0) then
                D2  = 3.d0*(s(i,j,k-1)-2.d0*sedgez(i,j,k)+s(i,j,k))
                D2L = s(i,j,k-2)-2.d0*s(i,j,k-1)+s(i,j,k)
                D2R = s(i,j,k-1)-2.d0*s(i,j,k)+s(i,j,k+1)
                sgn = sign(1.d0,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                sedgez(i,j,k) = 0.5d0*(s(i,j,k-1)+s(i,j,k)) - (1.d0/6.d0)*D2LIM
             end if
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    k = k3d
    !$OMP PARALLEL DO PRIVATE(i,j,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep,dabarm,dabarp,dafacemin) &
    !$OMP PRIVATE(dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax,delam,delap,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedgez(i,j,k+1)-s(i,j,k)
          alpham   = sedgez(i,j,k  )-s(i,j,k)
          bigp     = abs(alphap).gt.2.d0*abs(alpham)
          bigm     = abs(alpham).gt.2.d0*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. 0.d0) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedgez(i,j,k) - sedgez(i,j,k-1)
             dafacep   = sedgez(i,j,k+2) - sedgez(i,j,k+1)
             dabarm    = s(i,j,k) - s(i,j,k-1)
             dabarp    = s(i,j,k+1) - s(i,j,k)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2     = 6.d0*(alpham + alphap)
             D2L    = s(i,j,k-2)-2.d0*s(i,j,k-1)+s(i,j,k)
             D2R    = s(i,j,k)-2.d0*s(i,j,k+1)+s(i,j,k+2)
             D2C    = s(i,j,k-1)-2.d0*s(i,j,k)+s(i,j,k+1)
             sgn    = sign(1.d0,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(1.d0,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i,j,k-1) - s(i,j,k)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -2.d0*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(1.d0,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i,j,k+1) - s(i,j,k)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -2.d0*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k) + alpham
          sp(i,j) = s(i,j,k) + alphap
          !
          ! Compute z-component of Ip and Im.
          !
          s6    = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))

          sigma = abs(u(i,j,k3d,3)-cspd(i,j,k3d))*dt/dz
          Ip(i,j,kc,3,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,3,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)

          sigma = abs(u(i,j,k3d,3))*dt/dz
          Ip(i,j,kc,3,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,3,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)

          sigma = abs(u(i,j,k3d,3)+cspd(i,j,k3d))*dt/dz
          Ip(i,j,kc,3,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,3,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)

       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(dsvl,dsvlm,dsvlp,sp,sm,sedgez)

  end subroutine ppm_type2

end module ppm_module

