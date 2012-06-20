      subroutine ls_phiupd(flag,phi,phi_l1,phi_l2,phi_h1,phi_h2, & 
           phin,phin_l1,phin_l2,phin_h1,phin_h2, & 
           uadv,uadv_l1,uadv_l2,uadv_h1,uadv_h2, & 
           vadv,vadv_l1,vadv_l2,vadv_h1,vadv_h2, & 
           nband, nbandsize, mine, minesize, &
           lo, hi, dt, dx, &
           type,type_l1,type_l2,type_h1,type_h2)
        
        use LS_probdata_module, only: LSorder, kapa, kapb, nbandwidth, mineloc, LARGEINT
        
        implicit none
        
        integer  flag
        integer  phi_l1,phi_l2,phi_h1,phi_h2
        integer  phin_l1,phin_l2,phin_h1,phin_h2
        integer  uadv_l1,uadv_l2,uadv_h1,uadv_h2
        integer  vadv_l1,vadv_l2,vadv_h1,vadv_h2
        integer  type_l1,type_l2,type_h1,type_h2
        
        double precision ::  phi(phi_l1:phi_h1,phi_l2:phi_h2)
        double precision :: phin(phin_l1:phin_h1,phin_l2:phin_h2)
        double precision :: uadv(uadv_l1:uadv_h1,uadv_l2:uadv_h2)
        double precision :: vadv(vadv_l1:vadv_h1,vadv_l2:vadv_h2)
        integer          :: type(type_l1:type_h1,type_l2:type_h2)
        
        integer  nbandsize, minesize
        integer  nband(nbandsize,2)
        integer  mine(minesize,2)
        integer  lo(2), hi(2)
        double precision   dt
        double precision   dx(2)
        
        integer  i, j, p,s,t
        double precision   Dxm, Dxp
        double precision   Dym, Dyp
        double precision   phix, phiy, phixx, phiyy, phixy, kappa
        double precision   uavg, vavg, Fo, Fkappa, Fadv
        double precision   SWITCH
        double precision   Dxpp, Dxmm, Dxpm, Dypp, Dymm, Dypm
        double precision   Ad, Bd, Cd, Dd, Delp, Delm
        
        p=1

        do while (nband(p,1) .gt. -LARGEINT)

           i = nband(p,1)
           j = nband(p,2)
           p = p + 1

           if (max(type(i+1,j), type(i-1,j), type(i,j+1), type(i,j-1), &
                type(i+1,j+1), type(i-1,j-1), type(i+1,j-1), type(i-1,j+1)) .le. 1 ) then

              Dxm   = ( phi(i,j) - phi(i-1,j) ) / dx(1)
              Dxp   = ( phi(i+1,j) - phi(i,j) ) / dx(1)
              Dym   = ( phi(i,j) - phi(i,j-1) ) / dx(2)
              Dyp   = ( phi(i,j+1) - phi(i,j) ) / dx(2)

              phix  = ( phi(i+1,j) - phi(i-1,j) ) / (2*dx(1))
              phiy  = ( phi(i,j+1) - phi(i,j-1) ) / (2*dx(2))

              phixx = ( phi(i+1,j) - 2*phi(i,j) + phi(i-1,j) ) &
                   / (dx(1)**2)
              phiyy = ( phi(i,j+1) - 2*phi(i,j) + phi(i,j-1) )&
                   / (dx(2)**2)
              phixy = ( phi(i+1,j+1) + phi(i-1,j-1) &
                   - phi(i+1,j-1) - phi(i-1,j+1) ) &
                   / (4*dx(1)*dx(2))

              if (phix**2 + phiy**2 .gt. 0 ) then
                 Fkappa = -kapb*( phixx*phiy**2 - 2*phiy*phix*phixy &
                      + phiyy*phix**2 ) / (phix**2 + phiy**2) 
              else
                 Fkappa = 0
              endif


              if (LSorder .eq. 2) then

                 Dxpp = 0
                 if (i+2 .le. hi(1) +1) then
                    if(type(i+2,j) .lt. 2) then
                       Dxpp = (phi(i+2,j) - 2*phi(i+1,j) + phi(i,j))/dx(1)
                    endif
                 endif

                 Dxmm = 0
                 if (i-2 .ge. lo(1)) then
                    if (type(i-2,j) .lt. 2) then
                       Dxmm = (phi(i,j) - 2* phi(i-1,j) + phi(i-2,j))/dx(1)
                    endif
                 endif

                 Dxpm = (phi(i+1,j) - 2* phi(i,j) + phi(i-1,j))/dx(1)

                 Dypp = 0
                 if(j+2 .le. hi(2)) then
                    if (type(i,j+2) .lt. 2) then
                       Dypp = (phi(i,j+2) - 2* phi(i,j+1) + phi(i,j))/dx(2)
                    endif
                 endif

                 Dymm = 0
                 if(j-2 .ge. lo(2)) then
                    if(type(i,j-2) .lt. 2) then
                       Dymm = (phi(i,j) - 2* phi(i,j-1) + phi(i,j-2))/dx(2)
                    endif
                 endif

                 Dypm = (phi(i,j+1) - 2* phi(i,j) + phi(i,j-1))/dx(2)


                 Ad = Dxm + .5*SWITCH(Dxmm,Dxpm)
                 Bd = Dxp + .5*SWITCH(Dxpp,Dxpm)
                 Cd = Dym + .5*SWITCH(Dymm,Dypm)
                 Dd = Dyp + .5*SWITCH(Dypp,Dypm)

                 Delp = (max(Ad,0.d0)**2 + min(Bd,0.d0)**2 + max(Cd,0.d0)**2 + min(Dd,0.d0)**2)**(.5)

                 Delm = (max(Bd,0.d0)**2 + min(Ad,0.d0)**2 + max(Dd,0.d0)**2 + min(Cd,0.d0)**2)**(.5)

              endif

              Fo = kapa

              if (LSorder .eq. 1) then

                 if (Fo .gt. 0) then
                    Fo = Fo*( ( max(Dxm,0.d0) + min(Dxp,0.d0) )**2  &
                         + ( max(Dym,0.d0) + min(Dyp,0.d0) )**2  )**(1./2.) 
                 else
                    Fo = Fo*( ( min(Dxm,0.d0) + max(Dxp,0.d0) )**2  &
                         + ( min(Dym,0.d0) + max(Dyp,0.d0) )**2  )**(1./2.) 
                 endif

              else if (LSorder .eq. 2) then

                 if (Fo .gt. 0) then
                    Fo = Fo*Delp
                 else
                    Fo = Fo*Delm
                 endif
              endif

              Fadv = 0
              uavg = ( uadv(i,j) + uadv(i+1,j) ) / 2.
              vavg = ( vadv(i,j) + vadv(i,j+1) ) / 2.

              if (LSorder .eq. 1) then

                 if (uavg .gt. 0) then
                    Fadv = Fadv + uavg*Dxm
                 else 
                    Fadv = Fadv + uavg*Dxp
                 endif

                 if (vavg .gt. 0) then
                    Fadv = Fadv + vavg*Dym
                 else 
                    Fadv = Fadv + vavg*Dyp
                 endif

              else if (LSorder .eq. 2) then
                 Fadv = uavg*phix + vavg*phiy
              endif

              phin(i,j) = phi(i,j) - dt*( Fo + Fkappa + Fadv )

           endif
        enddo

        flag = 0

        p = 1

        do while (mine(p,1) .gt. -LARGEINT)
           i = mine(p,1)
           j = mine(p,2)
           p = p + 1

           if (sign(1.d0,phi(i,j))*sign(1.d0,phin(i,j)) .le. 0) then
              flag = 1
              exit
           endif
        enddo

        return 

      end subroutine ls_phiupd

      double precision function SWITCH(x,y)
        implicit none      
        double precision x,y
        
        if (x*y .ge. 0) then
           if (abs(x) .le. abs(y)) then
              SWITCH = x
           else 
              SWITCH = y
           endif
        else
           SWITCH = 0
        endif
      end function SWITCH
      

      subroutine ls_cfl(lscfl,phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                        uadv,uadv_l1,uadv_l2,uadv_h1,uadv_h2, & 
                        vadv,vadv_l1,vadv_l2,vadv_h1,vadv_h2, & 
                        nband, nbandsize, mine, minesize, &
                        lo, hi, phit, dx, &
                        type,type_l1,type_l2,type_h1,type_h2)

      use LS_probdata_module, only: kapa, kapb, nbandwidth, LARGEINT

      implicit none

      integer  phi_l1,phi_l2,phi_h1,phi_h2
      integer  uadv_l1,uadv_l2,uadv_h1,uadv_h2
      integer  vadv_l1,vadv_l2,vadv_h1,vadv_h2
      integer  type_l1,type_l2,type_h1,type_h2

      double precision ::  phi(phi_l1:phi_h1,phi_l2:phi_h2)
      double precision :: uadv(uadv_l1:uadv_h1,uadv_l2:uadv_h2)
      double precision :: vadv(vadv_l1:vadv_h1,vadv_l2:vadv_h2)
      integer          :: type(type_l1:type_h1,type_l2:type_h2)

      integer  nbandsize, minesize
      integer  nband(nbandsize,2)
      integer  mine(minesize,2)
      integer  lo(2), hi(2)
      double precision :: lscfl
      double precision ::  dx(2), phit

      integer i, j, p
      double precision phix, phiy, phixx, phiyy, phixy, kappa, speed
      double precision phidt

      phidt = phit
      
      p = 1
      do while (nband(p,1) .gt. -LARGEINT)
         
         i = nband(p,1)
         j = nband(p,2)
         p = p + 1
         
         if(max(type(i+1,j),type(i-1,j),type(i,j+1),type(i,j-1), &
              type(i+1,j+1),type(i-1,j-1), type(i+1,j-1), type(i-1,j+1)) .le. 1 ) then    
          
            phix = (phi(i+1,j)-phi(i-1,j))/(2*dx(1))
            phiy = (phi(i,j+1)-phi(i,j-1))/(2*dx(2))
            phixx = (phi(i+1,j)-2*phi(i,j)+phi(i-1,j))/(dx(1)**2)
            phiyy = ( phi(i,j+1) - 2*phi(i,j) + phi(i,j-1) )/(dx(2)**2)
            phixy = (phi(i+1,j+1) + phi(i-1,j-1) - phi(i+1,j-1) - phi(i-1,j+1))/(4*dx(1)*dx(2))
            if (phix**2 + phiy**2 .gt. 0 ) then
               kappa = (phixx*phiy**2 - 2*phiy*phix*phixy + phiyy*phix**2)/((phix**2+phiy**2)**(3./2.))
            else
               kappa = 0.d0
            endif
            
            speed = ( (.5*( uadv(i,j) + uadv(i+1,j) ) )**2 + ( .5*( vadv(i,j) +vadv(i,j+1) ) )**2)**.5 &
                 + abs(kapa-kapb*kappa)
             
            phidt = min( phit, &
                         1/(max(1./phidt,speed/(.8*min(dx(1),dx(2))), 4*abs(kapb)/(.8*min(dx(1),dx(2))**2) ))  )     

         endif
      enddo
        
      lscfl = phidt

      end subroutine ls_cfl
      

      subroutine ls_findinterface(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                                  phin,phin_l1,phin_l2,phin_h1,phin_h2, & 
                                  type,type_l1,type_l2,type_h1,type_h2, &
                                  lo, hi, dx, intfacenump, intfacenumn, intfacep,intfacen, &
                                  nband, nbandsize, intfacesize)

      use LS_probdata_module, only: LARGEINT, BOGUS

      implicit none

      integer  phi_l1,phi_l2,phi_h1,phi_h2
      integer  phin_l1,phin_l2,phin_h1,phin_h2
      integer  type_l1,type_l2,type_h1,type_h2

      double precision ::  phi(phi_l1:phi_h1,phi_l2:phi_h2)
      double precision :: phin(phin_l1:phin_h1,phin_l2:phin_h2)
      integer          :: type(type_l1:type_h1,type_l2:type_h2)
      
      integer    lo(2), hi(2)    
      double precision     dx(2)
      integer    intfacenump, intfacenumn, intfacesize
      integer    intfacep(intfacesize,2), intfacen(intfacesize,2)
      integer    nbandsize
      integer    nband(nbandsize,2)
      
!     Local variables
      integer    i, j, r
      
      intfacenump=0
      intfacenumn=0

      r = 1
      do while(nband(r,1) .gt. -LARGEINT)
        i = nband(r,1)
        j = nband(r,2)
        r = r + 1      
        
        phin(i,j) = sign(BOGUS,phi(i,j))
      enddo
      
      r = 1
      do while(nband(r,1) .gt. -LARGEINT)
        i = nband(r,1)
        j = nband(r,2)
        r = r + 1
        
        call UPDATEF(i, j, phi, phi_l1, phi_l2, phi_h1, phi_h2, &
                     phin, phin_l1, phin_l2, phin_h1, phin_h2, &
                     type, type_l1, type_l2, type_h1, type_h2, &
                     intfacep, intfacen, intfacenump, intfacenumn, &
                     dx, intfacesize, lo, hi) 
           
        if( (i .eq. lo(1)) .OR. (j .eq. lo(2)) ) then
           
           if(i .eq. lo(1)) then
                call UPDATEF(i-1, j, phi, phi_l1, phi_l2, phi_h1, phi_h2, &
                             phin, phin_l1, phin_l2, phin_h1, phin_h2, &
                             type, type_l1, type_l2, type_h1, type_h2, &
                             intfacep, intfacen, intfacenump, intfacenumn,  &
                             dx, intfacesize, lo, hi )
             endif
             
             if(j .eq. lo(2)) then 
                call UPDATEF(i, j-1, phi, phi_l1, phi_l2, phi_h1, phi_h2, &
                             phin, phin_l1, phin_l2, phin_h1, phin_h2, &
                             type, type_l1, type_l2, type_h1, type_h2, &
                             intfacep, intfacen, intfacenump, intfacenumn,  &
                             dx, intfacesize, lo, hi )
             endif
             
             if (i .eq. lo(1) .and. j .eq. lo(2)) then
                call UPDATEF(i-1, j-1, phi, phi_l1, phi_l2, phi_h1, phi_h2, &
                             phin, phin_l1, phin_l2, phin_h1, phin_h2, &
                             type, type_l1, type_l2, type_h1, type_h2, &
                             intfacep, intfacen, intfacenump, intfacenumn,  &
                             dx, intfacesize, lo, hi )
             endif
          endif
      enddo
      	
      nband(1,1) = -LARGEINT
      nband(1,2) = -LARGEINT
      end


      subroutine UPDATEF(i,j,phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                         phin,phin_l1,phin_l2,phin_h1,phin_h2, & 
                         type,type_l1,type_l2,type_h1,type_h2, &
                         intfacep,intfacen,intfacenump,intfacenumn,  &
                         dx,intfacesize,lo,hi)
      
      use LS_probdata_module, only: BOGUS

      implicit none      

      integer  phi_l1,phi_l2,phi_h1,phi_h2
      integer  phin_l1,phin_l2,phin_h1,phin_h2
      integer  type_l1,type_l2,type_h1,type_h2

      double precision :: phi (phi_l1:phi_h1,phi_l2:phi_h2)
      double precision :: phin(phin_l1:phin_h1,phin_l2:phin_h2)
      integer          :: type(type_l1:type_h1,type_l2:type_h2)

      integer    intfacenump, intfacenumn, intfacesize
      integer    intfacep(intfacesize,2), intfacen(intfacesize,2)
      double precision     dx(2)
      integer    lo(2), hi(2)
      
!     Local variables
      integer    c,d,i, j, ii, jj, iii, jjj, k,l, m,n,p,q,r,s,t
      double precision     A(16,16)
      double precision     B(16)
      double precision     x,y
      double precision     FINDDIST, distance
      double precision     work(16)
      double precision     grad(2)
      integer              iwork(16)
      integer              max0
      integer              nrow(16)
      double precision     geppX(16)
       
      if(max( abs(phi(i+1,j)),abs(phi(i,j+1)),abs(phi(i+1,j+1))) .lt. BOGUS) then
         
         
         if ( phi(i,j)*phi(i+1,j+1)   .lt. 0 &
             .OR. phi(i,j)*phi(i+1,j) .lt. 0  &
             .OR. phi(i,j)*phi(i,j+1) .lt. 0 ) then


            m=0
            do ii=0,1
               do jj=0,1
                  x = ii*dx(1)
                  y = jj*dx(2)
                  do n = 0,15
		     c = n/4
	       	     d = n-4*(n/4)
	 	     A(m+1,n+1) = x**(c) * y**(d)
	 	     A(m+2,n+1) = c * x**(max(c-1,0)) * y**(d)
	 	     A(m+3,n+1) = d * x**(c) * y**(max(d-1,0))
	 	     A(m+4,n+1) = c * d * x**(max(c-1,0)) * y**(max(d-1,0))                    
                  enddo
	     	  B(m+1) = phi(i+ii,j+jj)
	     	  B(m+2) = (phi(i+ii+1,j+jj) - phi(i+ii-1,j+jj))/(2*dx(1))
	     	  B(m+3) = (phi(i+ii,j+jj+1) - phi(i+ii,j+jj-1))/(2*dx(2))
	     	  B(m+4) = (phi(i+ii+1,j+jj+1) - phi(i+ii-1,j+jj+1) - phi(i+ii+1,j+jj-1) + phi(i+ii-1,j+jj-1))/(4*dx(1)*dx(2))
	     	  m = m + 4
               enddo
            enddo
            
            call gepp(A,B,16,geppX,nrow)
            
            do ii=0,1
               do jj=0,1
                  
                  iii = i + ii
                  jjj = j + jj
                  
                  distance = FINDDIST(grad,B,sign(1.d0,phi(iii,jjj)),ii*dx(1),jj*dx(2),dx)              
                  
                  if (type(iii,jjj) .NE. 0 .and. phi(iii,jjj) .ge. 0 .and. distance .ge. 0 &
                       .and. iii .ge. lo(1) .and. iii .le. hi(1) &
                       .and. jjj .ge. lo(2) .and. jjj .le. hi(2)) then
                     
                     intfacenump = intfacenump + 1
                     intfacep(intfacenump,1)=iii
                     intfacep(intfacenump,2)=jjj
                     
                     
                  else if (type(iii,jjj) .NE. 0 .and. distance .ge. 0 &
                          .and. iii .ge. lo(1) .and. iii .le. hi(1) &
                          .and. jjj .ge. lo(2) .and. jjj .le. hi(2))  then
                     
                     intfacenumn = intfacenumn + 1
                     intfacen(intfacenumn,1)=iii
                     intfacen(intfacenumn,2)=jjj
                     
                  endif
                  
                  if (distance .ge. 0) then
                     type(iii,jjj) = 0
                     phin(iii,jjj) = min(abs(phin(iii,jjj)),distance)*sign(1.d0,phi(iii,jjj))
                  endif
                  
               enddo
            enddo
            
         endif
      endif
      return
      end

      double precision function FINDDIST(grad,B,sgn,x0,y0,dx)

      use LS_probdata_module, only: BOGUS

      implicit none
      double precision :: grad(2)
      double precision :: B(16)
      double precision :: sgn
      double precision :: x0,y0
      double precision :: dx(2)

!     Local variables
      double precision t,tp
      double precision POLYVAL
      double precision DPOLYVAL
      double precision FA,FP
      double precision a,d,p
      double precision x,y
      double precision delta1(2), delta2(2)
      integer i    
      integer ITERMAX
      parameter (ITERMAX=30)
      
      x = x0
      y = y0 
      
      i = 0

      delta1(1) = BOGUS
      delta1(2) = BOGUS
      delta2(1) = BOGUS
      delta2(2) = BOGUS
      
      do while  ( (delta1(1)**2 + delta1(2)**2 + delta2(1)**2 + delta2(2)**2)**(.5) &
           .GT. 10.0**(-6.0)*dx(1)*dx(2) .AND. i .LT. ITERMAX)
        
        CALL GRADPVAL(B,grad,x,y)
        
      	delta1(1) = -polyval(B,x,y) * grad(1)/(grad(1)**2 + grad(2)**2)
      	delta1(2) = -polyval(B,x,y) * grad(2)/(grad(1)**2 + grad(2)**2)
      	
      	delta2(1) = (x0 - x) - grad(1)*( (x0 - x)*grad(1) + (y0 - y)*grad(2) ) /(grad(1)**2 + grad(2)**2)
      	delta2(2) = (y0 - y) - grad(2)*( (x0 - x)*grad(1) + (y0 - y)*grad(2) ) /(grad(1)**2 + grad(2)**2) 
      	
      	x = x + delta1(1) + delta2(1)
      	y = y + delta1(2) + delta2(2)
      	
      	i = i + 1
      
      enddo
      
      if (i .GE. 30 .OR. x .LT. 0 .OR. x .GT. dx(1) .OR. y .LT. 0 .OR. y .GT. dx(2)) then
      
      	FINDDIST = -1
      	
      else
      
      	FINDDIST = ( (x - x0)**2 + (y-y0)**2 )**(.5)
      
      endif
      

      end function      
      
      double precision function POLYVAL(B,x,y)
      double precision B(16)
      double precision x,y
      integer c,d,n  
      POLYVAL=0.d0
      do n=0,15
        c = n/4
        d = n-4*(n/4)
        POLYVAL = POLYVAL + B(n+1)*x**(c)*y**(d)
      enddo
      
      end function
      
      subroutine GRADPVAL(B,grad,x,y)

      implicit none

      double precision :: B(16)
      double precision :: grad(2)
      double precision :: x,y      

!     Local variables
      integer c,d,n
      
      grad(1) = 0
      grad(2) = 0
      
      do n = 0,15
      
	c = n/4
	d = n-4*(n/4)
      
      	grad(1) = grad(1) + B(n+1) * c * x**(max(c-1,0)) * y**(d)
      	grad(2) = grad(2) + B(n+1) * d * x**(c) * y**(max(d-1,0))
      
      enddo

      return
      end       
      
      double precision function DPOLYVAL(B,grad,x,y)
      implicit none
      double precision :: B(16)
      double precision :: grad(2)
      double precision :: x,y      
      integer n
      DPOLYVAL= (3*x**2*(B(1)*y**3 + B(2)*y**2 + B(3)*y + B(4)) + 2*x*(B(5)*y**3 + B(6)*y**2 &
                 + B(7)*y + B(8))  + (B(9)*y**3 + B(10)*y**2 + B(11)*y + B(12)) )*grad(1)  &
                 + (3*y**2*(B(1)*x**3 + B(5)*x**2 + B(9)*x + B(13)) + 2*y*(B(2)*x**3 &
                 + B(6)*x**2 + B(10)*x + B(14))  + (B(3)*x**3 + B(7)*x**2 + B(11)*x + B(15)) )*grad(2)
      end function      
      
      subroutine gepp(A,B,N,X,nrow)
      implicit none
      integer          :: N
      double precision :: A(N,N)
      double precision :: B(N)
      double precision :: X(N)
      integer          :: nrow(N)

      integer          :: i,j,k,p,ncopy
      double precision :: m,mi,mj
      
      do i = 1,N
        nrow(i) = i
      enddo
      
      do i = 1,N-1
         mi = abs(A(nrow(i),i))
         p=i
         do j = i+1, N
            mj = abs(A(nrow(j),i))
            if (mi < mj) then
               mi = mj
               p = j
            endif
         enddo
         
         if (nrow(i) .NE. nrow(p)) then
            ncopy = nrow(i)
            nrow(i) = nrow(p)
            nrow(p) = ncopy
         endif

         do j = i+1, N
            m = A(nrow(j),i)/A(nrow(i),i)
            do k = 1,N
               A(nrow(j),k) = A(nrow(j),k)- m*A(nrow(i),k)        
            enddo
            B(nrow(j)) = B(nrow(j))- m*B(nrow(i))
        enddo
      enddo
        
      X(N) = B(nrow(N))/A(nrow(N),N)
      do i= N, 1, -1
        X(i) = B(nrow(i))
        do j = i+1,N
          X(i) = X(i) - A(nrow(i),j)*X(j)
        enddo
        X(i) = X(i)/A(nrow(i),i)
      enddo
      
      do i = 1, N
        B(i) = X(i)
      enddo
      end 

      subroutine ls_narrowband(type,type_l1,type_l2,type_h1,type_h2,& 
      		               nband,nbandsize,mine,minesize,lo,hi)

      use LS_probdata_module, only: LARGEINT
        
      implicit none     
     
      integer type_l1,type_l2,type_h1,type_h2
      integer type(type_l1:type_h1,type_l2:type_h2)

      integer nbandsize, minesize
      integer lo(2), hi(2)
      integer nband(nbandsize,2)
      integer  mine(minesize,2)
     
      integer numband
      integer nummine
      integer i,j
     
      numband = 0
      nummine = 0
     
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            
            if(type(i,j) .eq. 0  ) then
               numband = numband + 1
               nband(numband,1) = i
               nband(numband,2) = j
            else if (type(i,j) .eq. 1) then
               numband = numband + 1
               nband(numband,1) = i
               nband(numband,2) = j          
               
               nummine = nummine + 1
               mine(nummine,1) = i
               mine(nummine,2) = j
            endif
         enddo
      enddo
     
      nband(numband + 1,1)  = -LARGEINT
      nband(numband + 1,2)  = -LARGEINT
    
      mine(nummine + 1,1) = -LARGEINT
      mine(nummine + 1,2) = -LARGEINT
      
      end
     
      subroutine ls_retypify(type,type_l1,type_l2,type_h1,type_h2,nband,nbandsize)     

      use LS_probdata_module, only: LARGEINT

      implicit none  

      integer type_l1,type_l2,type_h1,type_h2
      integer type(type_l1:type_h1,type_l2:type_h2)

      integer nbandsize
      integer nband(nbandsize,2)
      integer i, j, p

      p = 1
      do while (nband(p,1) .gt. -LARGEINT)
        i = nband(p,1)
        j = nband(p,2)
        p = p + 1
        type(i,j) = 3
      enddo
      end

      subroutine ls_fastmarch(phi,phi_l1,phi_l2,phi_h1,phi_h2, & 
                              type,type_l1,type_l2,type_h1,type_h2, &
                              lo, hi, dx, intfacenum, intface,  &
                              nband, nbandsize, nbandnum, mine, &
                              sgn, intfacesize,heap, heaploc)

      use LS_probdata_module, only: nbandwidth, mineloc, BOGUS, LARGEINT

      implicit none

      integer phi_l1,phi_l2,phi_h1,phi_h2
      integer type_l1,type_l2,type_h1,type_h2

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
      integer          :: type(type_l1:type_h1,type_l2:type_h2)

      integer heaploc(type_l1:type_h1,type_l2:type_h2)

      integer     lo(2), hi(2)
      double precision      dx(2)
      integer     intfacenum, intfacesize
      integer     intface(intfacesize,2)
      integer     nbandsize
      integer     nband(nbandsize,2)
      integer     mine(nbandsize,2)
      integer     nbandnum
      integer     sgn
      integer heap(nbandsize,2)

      integer i, j, n, p
      integer numtent
      
      numtent = 0
      do n = 1, intfacenum
        i = intface(n,1)
        j = intface(n,2)
      
	CALL UPDATE(phi,i,j,sgn, type, heap,numtent, &
                    phi_l1,phi_l2,phi_h1,phi_h2, &
                    type_l1,type_l2,type_h1,type_h2, &
     &    	    nbandsize, lo, hi, dx, heaploc)
	
	nbandnum = nbandnum + 1
	nband(nbandnum,1) = i
	nband(nbandnum,2) = j
      enddo
      
      do while (numtent .gt. 0  )

         CALL RMVNODE(heap, i,j,numtent, &
                      phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                      nbandsize,heaploc, &
                      type_l1,type_l2,type_h1,type_h2)
     
         if (abs(phi(i,j)) .lt. nbandwidth) then
            nbandnum = nbandnum + 1
            nband(nbandnum,1) = i
            nband(nbandnum,2) = j
            type(i,j) = 0
         else
            type(i,j) = 3
            phi(i,j) = sign(BOGUS,phi(i,j))
            exit
         endif
         
         if (abs(phi(i,j)) .gt. mineloc &
              .and. abs(phi(i,j)) .lt. nbandwidth ) then
          
            type(i,j) = 1
            
         endif
            
         CALL UPDATE(phi,i,j,sgn, type, heap, numtent,&
                     phi_l1,phi_l2,phi_h1,phi_h2, &
                     type_l1,type_l2,type_h1,type_h2, &
                     nbandsize, lo, hi, dx,heaploc)
         
      enddo
      
      nband(nbandnum+1,1) = -LARGEINT
      nband(nbandnum+1,2) = -LARGEINT      
      do while (numtent .gt. 0)

         CALL RMVNODE(heap, i, j, numtent, &
                      phi,phi_l1,phi_l2,phi_h1,phi_h2,nbandsize,heaploc, &
                      type_l1,type_l2,type_h1,type_h2)
         type(i,j) =3
         phi(i,j) = sign(BOGUS,phi(i,j))

      enddo
      end

      subroutine UPDATE(phi, i, j, sgn, type, heap, numtent, &
                        phi_l1,phi_l2,phi_h1,phi_h2, & 
                        type_l1,type_l2,type_h1,type_h2, &
                        nbandsize, lo, hi, dx, heaploc)
      
      implicit none

      integer phi_l1,phi_l2,phi_h1,phi_h2
      integer type_l1,type_l2,type_h1,type_h2

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
      integer          :: type(type_l1:type_h1,type_l2:type_h2)
      integer          :: heaploc(type_l1:type_h1,type_l2:type_h2)

      integer i,j
      integer nbandsize
      integer heap(nbandsize,2)
      integer lo(2), hi(2)
      integer numtent
      integer sgn
      integer n,ii,jj
      double precision dx(2)
      
      do n = 1,4
         
         ii = i - 1 + n/2
         jj = j + 2*(n/3) - n/2
         
         if ( ii .ge. lo(1) .and. ii .le. hi(1) &
              .and. jj .ge. lo(2) .and. jj .le. hi(2)) then
            
            if (type(ii,jj) .gt. 1   .and.   sgn*phi(ii,jj) .ge. 0) then
               
               CALL EVAL(phi, ii, jj, phi_l1, phi_l2, phi_h1, phi_h2, &
                         type_l1, type_l2, type_h1, type_h2, lo, hi, type, sgn, dx)
               
               if (type(ii,jj) .gt. 2) then
                  
                  type(ii,jj) = 2
                  CALL ADDNODE(heap, ii, jj, numtent, &
                               phi, phi_l1, phi_l2, phi_h1, phi_h2, &
                               lo, hi, nbandsize,heaploc, &
                               type_l1, type_l2, type_h1, type_h2)
                  
               else

                  CALL UPDATENODE(heap, ii, jj, numtent, &
                                  phi, phi_l1, phi_l2, phi_h1, phi_h2, &
                                  lo, hi, nbandsize,heaploc, &
                                  type_l1, type_l2, type_h1, type_h2)
                  
               endif
            endif    
         endif
      enddo
      end
      
      
      
      subroutine EVAL(phi,i,j,phi_l1,phi_l2,phi_h1,phi_h2, &
                      type_l1,type_l2,type_h1,type_h2,lo,hi,type,sgn,dx)
      
      implicit none
      integer i,j     

      integer  phi_l1,phi_l2,phi_h1,phi_h2
      integer  type_l1,type_l2,type_h1,type_h2

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
      integer          :: type(type_l1:type_h1,type_l2:type_h2)

      integer lo(2), hi(2)
      double precision a,b,c
      double precision dx(2)
      integer sgn
      integer  left,right,up,down
      LOGICAL  lok, rok, uok, dok
      
      a = 0.d0
      b = 0.d0
      c = -dx(1)*dx(2)
      
      left  = i - 1
      right = i + 1
      up    = j + 1
      down  = j - 1

      lok  = left  .ge. lo(1) .and. sgn*phi(left,j)  .ge. 0 .and. type(left,j ) .le. 1
      rok  = right .le. hi(1) .and. sgn*phi(right,j) .ge. 0 .and. type(right,j) .le. 1 
      uok  = up    .le. hi(2) .and. sgn*phi(i,up)    .ge. 0 .and. type(i,up   ) .le. 1
      dok  = down  .ge. lo(2) .and. sgn*phi(i,down)  .ge. 0 .and. type(i,down ) .le. 1
!     llok = FMMorder .eq. 2 . and. ll .ge. lo(1) .and. sgn*phi(ll,j) .ge. 0 .and. type(ll,j) .le. 1
!     rrok = FMMorder .eq. 2 . and. rr. le. hi(1) .and. sgn*phi(rr,j) .ge. 0 . and. type(rr,j) .le. 1 
!     uuok = FMMorder .eq. 2 . and. uu .le. hi(2) .and. sgn*phi(i,uu) .ge. 0  .and. type(i,uu) .le. 1
!     ddok = FMMorder .eq. 2 . and. dd . ge. lo(2) .and. sgn*phi(i,dd) .ge. 0  .and. type(i,dd) .le. 1      
            
!     FIXME: The following had a sign(right,...), mistake?
      if (lok .and. rok) then
                    
        a = a + dx(2)/dx(1)
        b = b + min(sgn*phi(left,j),sgn*phi(right,j))*(dx(2)/dx(1))
        c = c + min(phi(left,j)**2,phi(right,j)**2)*(dx(2)/dx(1))
        
      else if (lok) then
      
        a = a + dx(2)/dx(1)
	b = b + sgn*phi(left,j)*(dx(2)/dx(1))
        c = c + phi(left,j)**2*(dx(2)/dx(1))
      
      else if (rok ) then

        a = a + dx(2)/dx(1)
	b = b + sgn*phi(right,j)*(dx(2)/dx(1))
        c = c + phi(right,j)**2*(dx(2)/dx(1))
        
      endif
      
      
      if (dok .and. uok) then
      
        a = a + dx(1)/dx(2)
        b = b + min(sgn*phi(i,down),sgn*phi(i,up))*(dx(1)/dx(2))
        c = c + min(phi(i,down)**2,phi(i,up)**2)*(dx(1)/dx(2))
        
      else if (dok) then
      
        a = a + dx(1)/dx(2)
	b = b + sgn*phi(i,down)*(dx(1)/dx(2))
        c = c + phi(i,down)**2*(dx(1)/dx(2))
      
      else if (uok ) then

        a = a + dx(1)/dx(2)
	b = b + sgn*phi(i,up)*(dx(1)/dx(2))
        c = c + phi(i,up)**2*(dx(1)/dx(2))
        
      endif      
        
      b = -2*b
      
      if (b**2 - 4*a*c .lt. 0) then
      
         phi(i,j) = sgn*(-b)/(2*a)
         return
      
      endif
      
      phi(i,j) = sgn*(-b + SQRT(b**2-4*a*c))/(2*a)

      end
 
 
 
      subroutine ADDNODE(heap,i,j,n,phi,phi_l1,phi_l2,phi_h1,phi_h2,lo,hi, &
                         nbandsize, heaploc,type_l1,type_l2,type_h1,type_h2)
      implicit none

      integer phi_l1,phi_l2,phi_h1,phi_h2
      integer type_l1,type_l2,type_h1,type_h2

      double precision phi(phi_l1:phi_h1,phi_l2:phi_h2)
      integer heaploc(type_l1:type_h1,type_l2:type_h2)

      integer nbandsize
      integer lo(2), hi(2)
      integer heap(nbandsize,2)
      integer i,j,n

      integer index
      integer parent
      
      index = n + 1
      parent = index/2

      if (n .eq. 0) then
         
         heap(index,1) = i
         heap(index,2) = j
         heaploc(i,j) = index        
         
         n=n+1
         
         return
      endif
      
      do while ( ABS(phi(heap(parent,1), heap(parent,2))) .gt. ABS(phi(i,j)))
            
        heap(index,1) = heap(parent,1)
        heap(index,2) = heap(parent,2)
        heaploc(heap(index,1),heap(index,2)) = index        
        index = parent
        parent = index/2
         
        if (parent .eq. 0) exit
      enddo
      
      heap(index,1) = i
      heap(index,2) = j
      heaploc(i,j) = index
      n = n + 1
      end


      subroutine UPDATENODE(heap,i,j,n,phi,phi_l1,phi_l2,phi_h1,phi_h2,lo,hi, &
                           nbandsize, heaploc,type_l1,type_l2,type_h1,type_h2)
      implicit none

      integer  phi_l1,phi_l2,phi_h1,phi_h2
      integer  type_l1,type_l2,type_h1,type_h2

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
      integer          :: heaploc(type_l1:type_h1,type_l2:type_h2)

      integer    nbandsize
      integer    lo(2), hi(2)
      integer    heap(nbandsize,2)
      integer    i,j, n

      integer index
      integer parent
      
      index = heaploc(i,j)
      parent = index/2

      if (index .eq. 1) then
        return
      endif
      
      do while ( ABS(phi(heap(parent,1 ), heap(parent,2))) .gt. ABS(phi(i,j)))
            
        heap(index,1) = heap(parent,1)
        heap(index,2) = heap(parent,2)
        heaploc(heap(index,1),heap(index,2)) = index
        index = parent
        parent = index/2
         
        if (parent .eq. 0) exit         
      enddo
      
      heap(index,1) = i
      heap(index,2) = j
      heaploc(i,j) = index
      end
            

      subroutine RMVNODE(heap,i,j,n,phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                         nbandsize, heaploc,type_l1,type_l2,type_h1,type_h2)
      implicit none

      integer  phi_l1,phi_l2,phi_h1,phi_h2
      integer  type_l1,type_l2,type_h1,type_h2

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
      integer          :: heaploc(type_l1:type_h1,type_l2:type_h2)

      integer    nbandsize
      integer    heap(nbandsize,2)
      integer    i,j, n
      integer    index, left, right

      i=heap(1,1)
      j=heap(1,2)
      heaploc(i,j) = -1
      
      index = 1
      left  = 2*index
      right = 2*index + 1
   
      do while (.TRUE.)
         
         if (left .le. n-1) then
            
            if( ABS(phi(heap(left,1),heap(left,2))) .lt. ABS(phi(heap(n,1),heap(n,2)))) then
               
               if (right .le. n-1) then
                  if (ABS(phi(heap(left,1),heap(left,2))) .lt. ABS(phi(heap(right,1),heap(right,2))) ) then
                     heap(index,1) = heap(left,1)
                     heap(index,2) = heap(left,2)
                     heaploc(heap(index,1),heap(index,2)) = index
                     index = left
                     left=2*index
                     right=2*index+1
                  else
                     heap(index,1) = heap(right,1)
                     heap(index,2) = heap(right,2)
                     heaploc(heap(index,1),heap(index,2)) = index	      
                     index = right
                     left=2*index
                     right=2*index+1  
                  endif
               else                  
                  heap(index,1) = heap(left,1)
                  heap(index,2) = heap(left,2)
                  heaploc(heap(index,1),heap(index,2)) = index	    
                  
                  index = left
                  left=2*index
                  right=2*index+1  
               endif     
               
            else if (right .le. n-1) then
         
               if (abs(phi(heap(right,1),heap(right,2))) .lt. abs( phi(heap(n,1),heap(n,2)))) then
                  heap(index,1) = heap(right,1)
                  heap(index,2) = heap(right,2)
                  heaploc(heap(index,1),heap(index,2)) = index	      
                  index = right
                  left=2*index
                  right=2*index+1
               else       
                  exit       
               endif
            else    
               exit      
            endif 
         else       
            exit      
         endif
      enddo
      
      heap(index,1) = heap(n,1)
      heap(index,2) = heap(n,2)
      if(n .gt. 1) then
        heaploc(heap(index,1),heap(index,2)) = index
      endif
      n = n - 1
      end
      
      subroutine ls_fastmarch2(flag,phi,phi_l1,phi_l2,phi_h1,phi_h2,&
                               type,type_l1,type_l2,type_h1,type_h2,&
                               lo, hi, dx, nband, nbandsize, nbandnum, &
                               sgn, heaploc)

      use LS_probdata_module, only: nbandwidth, mineloc, LARGEINT, BOGUS

      implicit none

      integer  flag
      integer  phi_l1,phi_l2,phi_h1,phi_h2
      integer  type_l1,type_l2,type_h1,type_h2

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
      integer          :: type(type_l1:type_h1,type_l2:type_h2)
      integer          :: heaploc(type_l1:type_h1,type_l2:type_h2)

      integer     lo(2), hi(2)
      double precision      dx(2)
      integer     nbandsize
      integer     nband(nbandsize,2)
      integer     nbandnum
      integer     sgn
      
      integer i, j, n
      integer heap(nbandsize,2)
      integer numtent
      
      numtent = 0
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            heaploc(i,j) = -1
         enddo
      enddo
      
      i = lo(1)-1 
      do j = lo(2), hi(2)
         if((type(i,j) .eq. 0 .OR. type(i,j) .eq. 1) .and. sgn*phi(i,j) .ge. 0) then
            CALL UPDATE2(phi,i,j,sgn, type,heap,numtent,phi_l1,phi_l2,phi_h1,phi_h2, &
                         type_l1,type_l2,type_h1,type_h2, &
                         nbandsize, lo, hi, dx,heaploc)
         endif
      enddo
      
      i = hi(1)+1 
      do j = lo(2), hi(2) 
         if((type(i,j) .eq. 0 .OR. type(i,j) .eq. 1) .and. sgn*phi(i,j) .ge. 0) then
            CALL UPDATE2(phi,i,j,sgn, type, heap,numtent,phi_l1,phi_l2,phi_h1,phi_h2, &
                         type_l1,type_l2,type_h1,type_h2, &
                         nbandsize, lo, hi, dx,heaploc)
         endif
      enddo
      
      j = lo(2)-1 
      do i = lo(1) , hi(1)
         if((type(i,j) .eq. 0 .OR. type(i,j) .eq. 1) .and. sgn*phi(i,j) .ge. 0) then
            CALL UPDATE2(phi,i,j,sgn, type, heap,numtent,phi_l1,phi_l2,phi_h1,phi_h2, &
                         type_l1,type_l2,type_h1,type_h2, &
                         nbandsize, lo, hi, dx,heaploc)
         endif
      enddo
      
      j = hi(2)+1 
      do i = lo(1) , hi(1) 
         if((type(i,j) .eq. 0 .OR. type(i,j) .eq. 1) .and. sgn*phi(i,j) .ge. 0) then
            CALL UPDATE2(phi,i,j,sgn, type, heap,numtent,phi_l1,phi_l2,phi_h1,phi_h2, &
                         type_l1,type_l2,type_h1,type_h2, &
                         nbandsize, lo, hi, dx,heaploc)
         endif
      enddo
      
      flag = 0
      do while (numtent .gt. 0 )
         CALL RMVNODE(heap, i,j,numtent, phi, phi_l1,phi_l2,phi_h1,phi_h2, &
                      nbandsize,heaploc, type_l1,type_l2,type_h1,type_h2)

         if (abs(phi(i,j)) .lt. nbandwidth ) then
         
            flag = 1

            if (type(i,j) .ge. 2) then
               nbandnum = nbandnum + 1         
               nband(nbandnum,1) = i
               nband(nbandnum,2) = j            
            endif
          
            type(i,j) = 0
            CALL UPDATE2(phi,i,j,sgn, type, heap,numtent,phi_l1,phi_l2,phi_h1,phi_h2, &
                         type_l1,type_l2,type_h1,type_h2, &
                         nbandsize, lo, hi, dx,heaploc)

         else

            type(i,j) = 3
            phi(i,j) = sign(BOGUS,phi(i,j))
            exit

         endif
         
         if (abs(phi(i,j)) .gt. mineloc &
              .and. abs(phi(i,j)) .lt. nbandwidth) then
          
            type(i,j) = 1
        endif
      enddo
      
      nband(nbandnum + 1,1) = -LARGEINT
      nband(nbandnum + 1,2) = -LARGEINT
      
      do while (numtent .gt. 0)
            
         CALL RMVNODE(heap, i,j,numtent, phi, phi_l1,phi_l2,phi_h1,phi_h2, &
                      nbandsize,heaploc,type_l1,type_l2,type_h1,type_h2)
         type(i,j) = 3
         phi(i,j) = sign(BOGUS,phi(i,j))    
      enddo

      end



      subroutine UPDATE2(phi,i,j,sgn, type, heap,numtent,phi_l1,phi_l2,phi_h1,phi_h2,&
                         type_l1,type_l2,type_h1,type_h2,nbandsize, lo, hi, dx, heaploc)      

      use LS_probdata_module, only: BOGUS

      implicit none

      integer  phi_l1,phi_l2,phi_h1,phi_h2
      integer  type_l1,type_l2,type_h1,type_h2

      double precision :: phi (phi_l1:phi_h1,phi_l2:phi_h2)
      integer          :: type(type_l1:type_h1,type_l2:type_h2)
      integer          :: heaploc(type_l1:type_h1,type_l2:type_h2)

      integer i,j
      integer nbandsize
      integer heap(nbandsize,2)
      integer lo(2), hi(2)
      integer numtent
      integer sgn
      integer n,ii,jj
      double precision dx(2)

      logical EVAL2
      logical isnew
      integer min
      
      do n = 1,4
      
        ii = i -1 +n/2
        jj = j + 2*(n/3) - n/2
              
        if (ii .ge. lo(1) .and. ii .le. hi(1) &
             .and. jj .ge. lo(2) .and. jj .le. hi(2)) then
        
           if (sgn*sign(1.d0,phi(ii,jj)) .ge. 0 .and. abs(phi(ii,jj)) .gt. abs(phi(i,j))  &
                .and. (abs(phi(ii,jj)) .ge. BOGUS .OR. .NOT.(sgn*phi(ii+1,jj) .le. 0   &
                .OR. sgn*phi(ii-1,jj) .le. 0 .OR. sgn*phi(ii,jj+1) .le. 0 .OR. sgn*phi(ii,jj-1) .le. 0) )) then
!                .OR. sgn*phi(ii+1,jj+1) .le. 0 .OR. sgn*phi(ii-1,jj-1) .le. 0 .OR. sgn*phi(ii+1,jj-1) .le. 0  &
!                .OR. sgn*phi(ii-1,jj+1) .le. 0) )) then
              
              isnew = EVAL2(phi,ii,jj,phi_l1,phi_l2,phi_h1,phi_h2, &
                            type_l1,type_l2,type_h1,type_h2, &
                            lo,hi, type,sgn,dx, phi(i,j))
           
              isnew = isnew .OR. EVAL2(phi,ii,jj,phi_l1,phi_l2,phi_h1,phi_h2, &
                                       type_l1,type_l2,type_h1,type_h2, &
                                       lo,hi, type,sgn,dx, phi(ii,jj))
              
              if (isnew ) then
                 type(ii,jj) = min(2,type(ii,jj))
                 if(heaploc(ii,jj) .eq. -1) then
                    CALL ADDNODE(heap,ii,jj, numtent,phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                                 lo, hi, nbandsize, heaploc, &
                                 type_l1,type_l2,type_h1,type_h2)
                 else
                    CALL UPDATENODE(heap,ii,jj, numtent,phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                                    lo, hi, nbandsize, heaploc, &
                                    type_l1,type_l2,type_h1,type_h2)
                 endif
              endif
           endif    
        endif        
      enddo
      end
      
      
      LOGICAL function EVAL2(phi,i,j,phi_l1,phi_l2,phi_h1,phi_h2, &
                             type_l1,type_l2,type_h1,type_h2, &
                             lo,hi, type,sgn,dx,phisrc)      
      implicit none

      integer  phi_l1,phi_l2,phi_h1,phi_h2
      integer  type_l1,type_l2,type_h1,type_h2

      double precision :: phi (phi_l1:phi_h1,phi_l2:phi_h2)
      integer          :: type(type_l1:type_h1,type_l2:type_h2)

      integer i,j     
      integer lo(2), hi(2)
      double precision a,b,c
      double precision dx(2)
      integer sgn
      integer  left,right,up,down
      double precision phisrc
      
      a = 0.0d0
      b = 0.0d0
      c = -dx(1)*dx(2)
      
      left  = i - 1
      right = i + 1
      up    = j + 1
      down  = j - 1
      
      if ( (sgn*phi(left,j) .ge. 0) .and. (type(left,j) .le. 1) .and. &
           (abs(phi(left,j)) .le. abs(phisrc)) .and.   &
           (sgn*phi(right,j) .ge. 0) .and. &
           (type(right,j) .le. 1) .and. (abs(phi(right,j)) .le. abs(phisrc)) ) then
         
         a = a + dx(2)/dx(1)        
         b = b + min(sgn*phi(left,j),sgn*phi(right,j))*(dx(2)/dx(1))
         c = c + min(phi(left,j)**2,phi(right,j)**2)*(dx(2)/dx(1))
         
      else if ( (sgn*phi(left,j) .ge. 0)  .and. (type(left,j) .le. 1) .and. &
           (abs(phi(left,j)) .le. abs(phisrc)) ) then
         
         a = a + dx(2)/dx(1)
         b = b + sgn*phi(left,j)*(dx(2)/dx(1))
         c = c + phi(left,j)**2*(dx(2)/dx(1))
         
      else if ( (sgn*phi(right,j) .ge. 0)  .and. (type(right,j) .le. 1) .and. &
           (abs(phi(right,j)) .le. abs(phisrc)) ) then
         
         a = a + dx(2)/dx(1)
         b = b + sgn*phi(right,j)*(dx(2)/dx(1))
         c = c + phi(right,j)**2*(dx(2)/dx(1))
         
      endif

      if ( sgn*phi(i,down) .ge. 0  .and. type(i,down) .le. 1 .and. abs(phi(i,down)) .le. abs(phisrc)  &
           .and. sgn*phi(i,up) .ge. 0  .and. type(i,up) .le. 1 .and. abs(phi(i,up)) .le. abs(phisrc) ) then
         
         a = a + dx(1)/dx(2)
         b = b + min(sgn*phi(i,down),sgn*phi(i,up))*(dx(1)/dx(2))
         c = c + min(phi(i,down)**2,phi(i,up)**2)*(dx(1)/dx(2))
         
      else if (sgn*phi(i,down) .ge. 0  .and. type(i,down) .le. 1 .and. abs(phi(i,down)) .le. abs(phisrc) ) then
         
         a = a + dx(1)/dx(2)
         b = b + sgn*phi(i,down)*(dx(1)/dx(2))
         c = c + phi(i,down)**2*(dx(1)/dx(2))
         
      else if (sgn*phi(i,up) .ge. 0  .and. type(i,up) .le. 1 .and. abs(phi(i,up)) .le. abs(phisrc) ) then
         
         a = a + dx(1)/dx(2)
         b = b + sgn*phi(i,up)*(dx(1)/dx(2))
         c = c + phi(i,up)**2*(dx(1)/dx(2))
         
      endif
      
      b = -2*b
      if (a .eq. 0.d0) then

        EVAL2 = .FALSE.
        return
       
      endif
      
      EVAL2 = .FALSE.
      
      if (b**2 - 4*a*c .lt. 0.d0) then
         if (ABS(phi(i,j)) .gt. (-b)/(2*a) + 1.d-10) then
            phi(i,j) = sgn*(-b)/(2*a)
            EVAL2 = .TRUE.
            return
         endif
      endif
      
      if (abs(phi(i,j)) .gt. (-b+SQRT(b**2-4*a*c))/(2*a)+1.d-10) then
        phi(i,j) = sgn*(-b+SQRT(b**2-4*a*c))/(2*a)
        EVAL2 = .TRUE.
      endif
      
      end


     
      subroutine ls_mine(type,type_l1,type_l2,type_h1,type_h2,nband, nbandsize,mine, minesize,lo, hi)    

      use LS_probdata_module, only: LARGEINT

      implicit none     

      integer type_l1,type_l2,type_h1,type_h2
      integer type(type_l1:type_h1,type_l2:type_h2)

      integer nbandsize, minesize
      integer lo(2), hi(2)
      integer  nband(nbandsize,2)
      integer  mine(minesize,2)
      
      integer i, j, p
      integer nummine
      
      nummine = 0
      p = 1
      do while(nband(p,1) .gt. -LARGEINT)
      
        i = nband(p,1)
        j = nband(p,2)        
        p = p + 1
        
        if (type(i,j) .eq. 1) then
          nummine = nummine +1
          mine(nummine,1) = i
          mine(nummine,2) = j
        endif
        
      enddo
      
      mine(nummine+1,1) = -LARGEINT
      mine(nummine+1,2) = -LARGEINT      
      end
     

     
      subroutine ls_nbandnumify(nband, nbandsize,nbandnum)

      use LS_probdata_module, only: LARGEINT

      implicit none      
      integer nbandsize, nbandnum
      integer nband(nbandsize,2)
      integer p
      
      p = 0
      do while(nband(p+1,1) .gt. -LARGEINT)
        p = p + 1
      enddo
      nbandnum = p
      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the density
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: ls        => level set array
! ::: nd        => number of components in den array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_lserror(tag,tagl1,tagl2,tagh1,tagh2, &
                            set,clear, &
                            ls,lsl1,lsl2,lsh1,lsh2, &
                            lo,hi,nd,domlo,domhi, &
                            delta,xlo,problo,time,level)
      use LS_probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer lsl1,lsl2,lsh1,lsh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision ls(lsl1:lsh1,lsl2:lsh2,nd)
      double precision delta(2), xlo(2), problo(2), time

      integer i, j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            if ( abs(ls(i,j,1)) .lt. lvlerr ) then
               tag(i,j) = set
            else
               tag(i,j) = tag(i,j)
            end if
         end do
      end do
      end
