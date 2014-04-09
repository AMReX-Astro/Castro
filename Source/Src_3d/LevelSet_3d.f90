      subroutine ls_phiupd(flag,phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, & 
                           phin,phin_l1,phin_l2,phin_l3,phin_h1,phin_h2,phin_h3, & 
                           uadv,uadv_l1,uadv_l2,uadv_l3,uadv_h1,uadv_h2,uadv_h3, & 
                           vadv,vadv_l1,vadv_l2,vadv_l3,vadv_h1,vadv_h2,vadv_h3, & 
                           wadv,wadv_l1,wadv_l2,wadv_l3,wadv_h1,wadv_h2,wadv_h3, & 
                           nband,nbandsize,mine,minesize,lo,hi,dt,dx, &
                           type,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
        
        use LS_probdata_module, only: LSorder, kapa, kapb, nbandwidth, mineloc, LARGEINT
        use bl_constants_module
        
        implicit none
        
        integer  flag
        integer  phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
        integer  phin_l1,phin_l2,phin_l3,phin_h1,phin_h2,phin_h3
        integer  uadv_l1,uadv_l2,uadv_l3,uadv_h1,uadv_h2,uadv_h3
        integer  vadv_l1,vadv_l2,vadv_l3,vadv_h1,vadv_h2,vadv_h3
        integer  wadv_l1,wadv_l2,wadv_l3,wadv_h1,wadv_h2,wadv_h3
        integer  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3
        
        double precision ::  phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
        double precision :: phin(phin_l1:phin_h1,phin_l2:phin_h2,phin_l3:phin_h3)
        double precision :: uadv(uadv_l1:uadv_h1,uadv_l2:uadv_h2,uadv_l3:uadv_h3)
        double precision :: vadv(vadv_l1:vadv_h1,vadv_l2:vadv_h2,vadv_l3:vadv_h3)
        double precision :: wadv(wadv_l1:wadv_h1,wadv_l2:wadv_h2,wadv_l3:wadv_h3)
        integer          :: type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)
        
        integer  nbandsize, minesize
        integer  nband(nbandsize,3)
        integer  mine(minesize,3)
        integer  lo(3), hi(3)
        double precision   dt
        double precision   dx(3)

        integer  i, j,k, p,s,t
        double precision   Dxm, Dxp
        double precision   Dym, Dyp
        double precision   Dzm, Dzp
        double precision   phix, phiy,phiz, phixx, phiyy, phizz, phixy, phixz, phiyz, kappa
        double precision   uavg, vavg, wavg, Fo, Fkappa, Fadv
        double precision   SWITCH
        double precision   Dxpp, Dxmm, Dxpm, Dypp, Dymm, Dypm, Dzpp, Dzmm, Dzpm
        double precision   Ad, Bd, Cd, Dd, Ed, Fd, Delp, Delm  
    
        p=1

        do while (nband(p,1) .GT. -LARGEINT)

           i = nband(p,1)
           j = nband(p,2)
           k = nband(p,3)
           p = p + 1

           if (max(type(i+1,j,k), type(i-1,j,k), type(i,j+1,k), type(i,j-1,k), &
                type(i,j,k+1), type(i,j,k-1), type(i+1,j+1,k), type(i-1,j-1,k), &
                type(i+1,j-1,k), type(i-1,j+1,k), type(i+1,j,k+1), type(i-1,j,k+1), &
                type(i+1,j,k-1), type(i-1,j,k-1), type(i,j+1,k+1), type(i,j-1,k-1), &
                type(i,j-1,k+1), type(i,j+1,k-1), type(i+1,j+1,k+1), type(i-1,j-1,k-1), &
                type(i+1,j-1,k-1), type(i-1,j+1,k+1), type(i+1,j+1,k-1), type(i-1,j-1,k+1), &
                type(i+1,j-1,k+1), type(i-1,j+1,k-1)) .LE. 1 ) then

              Dxm   = ( phi(i,j,k) - phi(i-1,j,k) ) / dx(1)
              Dxp   = ( phi(i+1,j,k) - phi(i,j,k) ) / dx(1)
              Dym   = ( phi(i,j,k) - phi(i,j-1,k) ) / dx(2)
              Dyp   = ( phi(i,j+1,k) - phi(i,j,k) ) / dx(2)
              Dzm   = ( phi(i,j,k) - phi(i,j,k-1) ) / dx(3)
              Dzp   = ( phi(i,j,k+1) - phi(i,j,k) ) / dx(3)

              phix  = ( phi(i+1,j,k) - phi(i-1,j,k) ) / (2*dx(1))
              phiy  = ( phi(i,j+1,k) - phi(i,j-1,k) ) / (2*dx(2))
              phiz  = ( phi(i,j,k+1) - phi(i,j,k-1) ) / (2*dx(3))

              phixx = ( phi(i+1,j,k) - 2*phi(i,j,k) + phi(i-1,j,k) ) / (dx(1)**2)
              phiyy = ( phi(i,j+1,k) - 2*phi(i,j,k) + phi(i,j-1,k) ) / (dx(2)**2)
              phizz = ( phi(i,j,k+1) - 2*phi(i,j,k) + phi(i,j,k-1) ) / (dx(3)**2)     
              phixy = ( phi(i+1,j+1,k) + phi(i-1,j-1,k) &
                   - phi(i+1,j-1,k) - phi(i-1,j+1,k) ) / (4*dx(1)*dx(2))
              phixz = ( phi(i+1,j,k+1) + phi(i-1,j,k-1) &
                   - phi(i+1,j,k-1) - phi(i-1,j,k+1) ) / (4*dx(1)*dx(3))     
              phiyz = ( phi(i,j+1,k+1) + phi(i,j-1,k-1) &
                   - phi(i,j-1,k+1) - phi(i,j+1,k-1) ) / (4*dx(1)*dx(3))     


              if (phix**2 + phiy**2 + phiz**2 .GT. 0 ) then
                 Fkappa = -kapb*( (phiyy + phizz)*phix**2 + (phixx + phizz)*phiy**2 &
                      + (phixx + phiyy)*phiz**2 - 2*phix*phiy*phixy - 2*phix*phiz*phixz &
                      - 2*phiy*phiz*phiyz) / (phix**2 + phiy**2 + phiz**2)      
              else
                 Fkappa = 0
              endif


              if (LSorder .EQ. 2) then

                 Dxpp = 0
                 if (i+2 .LE. hi(1) +1) then
                    if(type(i+2,j,k) .LT. 2) then
                       Dxpp = (phi(i+2,j,k) - 2*phi(i+1,j,k) + phi(i,j,k))/dx(1)
                    endif
                 endif

                 Dxmm = 0
                 if (i-2 .GE. lo(1)) then
                    if (type(i-2,j,k) .LT. 2) then
                       Dxmm = (phi(i,j,k) - 2* phi(i-1,j,k) + phi(i-2,j,k))/dx(1)
                    endif
                 endif

                 Dxpm = (phi(i+1,j,k) - 2* phi(i,j,k) + phi(i-1,j,k))/dx(1)

                 Dypp = 0
                 if(j+2 .LE. hi(2)) then
                    if (type(i,j+2,k) .LT. 2) then
                       Dypp = (phi(i,j+2,k) - 2* phi(i,j+1,k) + phi(i,j,k))/dx(2)
                    endif
                 endif

                 Dymm = 0
                 if(j-2 .GE. lo(2)) then
                    if(type(i,j-2,k) .LT. 2) then
                       Dymm = (phi(i,j,k) - 2* phi(i,j-1,k) + phi(i,j-2,k))/dx(2)
                    endif
                 endif

                 Dypm = (phi(i,j+1,k) - 2* phi(i,j,k) + phi(i,j-1,k))/dx(2)

                 Dzpp = 0
                 if (k+2 .LE. hi(3) +1) then
                    if(type(i,j,k+2) .LT. 2) then
                       Dzpp = (phi(i,j,k+2) - 2*phi(i,j,k+1) + phi(i,j,k))/dx(3)
                    endif
                 endif

                 Dzmm = 0
                 if (k-2 .GE. lo(3)) then
                    if (type(i,j,k-2) .LT. 2) then
                       Dzmm = (phi(i,j,k) - 2* phi(i,j,k-1) + phi(i,j,k-2))/dx(3)
                    endif
                 endif

                 Dzpm = (phi(i,j,k+1) - 2* phi(i,j,k) + phi(i,j,k-1))/dx(3)



                 Ad = Dxm + HALF*SWITCH(Dxmm,Dxpm)
                 Bd = Dxp + HALF*SWITCH(Dxpp,Dxpm)
                 Cd = Dym + HALF*SWITCH(Dymm,Dypm)
                 Dd = Dyp + HALF*SWITCH(Dypp,Dypm)
                 Ed = Dzm + HALF*SWITCH(Dzmm,Dzpm)
                 Fd = Dzp + HALF*SWITCH(Dzpp,Dzpm)

                 Delp = (max(Ad,ZERO)**2 + min(Bd,ZERO)**2 + max(Cd,ZERO)**2 &
                      + min(Dd,ZERO)**2 + max(Ed,ZERO)**2 + min(Fd,ZERO)**2)**(.5)

                 Delm = (max(Bd,ZERO)**2 + min(Ad,ZERO)**2 + max(Dd,ZERO)**2 &
                      + min(Cd,ZERO)**2 + max(Fd,ZERO)**2 + min(Ed,ZERO)**2)**(.5)

              endif

              Fo = kapa
              if (LSorder .EQ. 1) then
                 if (Fo .GT. 0) then
                    Fo = Fo*( ( max(Dxm,ZERO) + min(Dxp,ZERO) )**2 &
                         + ( max(Dym,ZERO) + min(Dyp,ZERO) )**2 &
                         + ( max(Dzm,ZERO) + min(Dzp,ZERO) )**2)**(1./2.) 
                 else
                    Fo = Fo*( ( min(Dxm,ZERO) + max(Dxp,ZERO) )**2 &
                         + ( min(Dym,ZERO) + max(Dyp,ZERO) )**2 &
                         + ( min(Dzm,ZERO) + max(Dzp,ZERO) )**2)**(1./2.) 
                 endif
              else if (LSorder .EQ. 2) then
                 if (Fo .GT. 0) then
                    Fo = Fo*Delp
                 else
                    Fo = Fo*Delm
                 endif
              endif

              Fadv = 0
              uavg = ( uadv(i,j,k) + uadv(i+1,j,k) ) / 2.
              vavg = ( vadv(i,j,k) + vadv(i,j+1,k) ) / 2.
              wavg = ( wadv(i,j,k) + wadv(i,j,k+1) ) / 2.

              if (LSorder .EQ. 1) then

                 if (uavg .GT. 0) then
                    Fadv = Fadv + uavg*Dxm
                 else 
                    Fadv = Fadv + uavg*Dxp
                 endif

                 if (vavg .GT. 0) then
                    Fadv = Fadv + vavg*Dym
                 else 
                    Fadv = Fadv + vavg*Dyp
                 endif

                 if (wavg .GT. 0) then
                    Fadv = Fadv + wavg*Dzm
                 else 
                    Fadv = Fadv + wavg*Dzp
                 endif

              else if (LSorder .EQ. 2) then
                 Fadv = uavg*phix + vavg*phiy + wavg*phiz
              endif

              phin(i,j,k) = phi(i,j,k) - dt*( Fo + Fkappa + Fadv )

           endif
        enddo

        flag = 0

        p = 1

        do while (mine(p,1) .GT. -LARGEINT)
           i = mine(p,1)
           j = mine(p,2)
           k = mine(p,3)
           p = p + 1

           if (sign(ONE,phi(i,j,k))*sign(ONE,phin(i,j,k)) .LE. 0) then
              flag = 1
              exit
           endif
        enddo

        return 

      end subroutine ls_phiupd


      double precision function SWITCH(x,y)
        implicit none      
        double precision x,y
        
        if (x*y .GE. 0) then
           if (abs(x) .LE. abs(y)) then
              SWITCH = x
           else 
              SWITCH = y
           endif
        else
           SWITCH = 0
        endif
      end function SWITCH

      subroutine ls_cfl(lscfl,phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                        uadv,uadv_l1,uadv_l2,uadv_l3,uadv_h1,uadv_h2,uadv_h3, & 
                        vadv,vadv_l1,vadv_l2,vadv_l3,vadv_h1,vadv_h2,vadv_h3, & 
                        wadv,wadv_l1,wadv_l2,wadv_l3,wadv_h1,wadv_h2,wadv_h3, & 
                        nband,nbandsize,mine,minesize,lo,hi,phit,dx, &
                        type,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
        
        use LS_probdata_module, only: kapa, kapb, nbandwidth, LARGEINT
        use bl_constants_module

        implicit none

        integer  phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
        integer  uadv_l1,uadv_l2,uadv_l3,uadv_h1,uadv_h2,uadv_h3
        integer  vadv_l1,vadv_l2,vadv_l3,vadv_h1,vadv_h2,vadv_h3
        integer  wadv_l1,wadv_l2,wadv_l3,wadv_h1,wadv_h2,wadv_h3
        integer  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

        double precision ::  phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
        double precision :: uadv(uadv_l1:uadv_h1,uadv_l2:uadv_h2,uadv_l3:uadv_h3)
        double precision :: vadv(vadv_l1:vadv_h1,vadv_l2:vadv_h2,vadv_l3:vadv_h3)
        double precision :: wadv(wadv_l1:wadv_h1,wadv_l2:wadv_h2,wadv_l3:wadv_h3)
        integer          :: type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

        integer  nbandsize, minesize
        integer  nband(nbandsize,3)
        integer  mine(minesize,3)
        integer  lo(3), hi(3)
        double precision :: lscfl
        double precision ::  dx(3), phit

        integer i, j, k, p
        double precision phix, phiy, phiz, phixx, phiyy, phizz, phixy, phixz, phiyz 
        double precision kappa, speed    
        double precision phidt

        phidt = phit

        p = 1
        do while (nband(p,1) .GT. -LARGEINT)

           i = nband(p,1)
           j = nband(p,2)
           k = nband(p,3)
           p = p + 1

           if (max(type(i+1,j,k), type(i-1,j,k), type(i,j+1,k), type(i,j-1,k), &
                type(i,j,k+1), type(i,j,k-1), type(i+1,j+1,k), type(i-1,j-1,k), &
                type(i+1,j-1,k), type(i-1,j+1,k), type(i+1,j,k+1), type(i-1,j,k+1), &
                type(i+1,j,k-1), type(i-1,j,k-1), type(i,j+1,k+1), type(i,j-1,k-1), &
                type(i,j-1,k+1), type(i,j+1,k-1), type(i+1,j+1,k+1), type(i-1,j-1,k-1), &
                type(i+1,j-1,k-1), type(i-1,j+1,k+1), type(i+1,j+1,k-1), type(i-1,j-1,k+1), &
                type(i+1,j-1,k+1), type(i-1,j+1,k-1)) .LE. 1 ) then

              phix  = ( phi(i+1,j,k) - phi(i-1,j,k) ) / (2*dx(1))
              phiy  = ( phi(i,j+1,k) - phi(i,j-1,k) ) / (2*dx(2))
              phiz  = ( phi(i,j,k+1) - phi(i,j,k-1) ) / (2*dx(3))

              phixx = ( phi(i+1,j,k) - 2*phi(i,j,k) + phi(i-1,j,k) ) / (dx(1)**2)
              phiyy = ( phi(i,j+1,k) - 2*phi(i,j,k) + phi(i,j-1,k) ) / (dx(2)**2)
              phizz = ( phi(i,j,k+1) - 2*phi(i,j,k) + phi(i,j,k-1) ) / (dx(3)**2)     
              phixy = ( phi(i+1,j+1,k) + phi(i-1,j-1,k) - phi(i+1,j-1,k) - phi(i-1,j+1,k) ) &
                   / (4*dx(1)*dx(2))
              phixz = ( phi(i+1,j,k+1) + phi(i-1,j,k-1) - phi(i+1,j,k-1) - phi(i-1,j,k+1) ) &
                   / (4*dx(1)*dx(3))     
              phiyz = ( phi(i,j+1,k+1) + phi(i,j-1,k-1) - phi(i,j-1,k+1) - phi(i,j+1,k-1) ) &
                   / (4*dx(1)*dx(3))     


              if (phix**2 + phiy**2 + phiz**2 .GT. 0 ) then
                 kappa = ( (phiyy + phizz)*phix**2 + (phixx + phizz)*phiy**2 &
                      + (phixx + phiyy)*phiz**2 - 2*phix*phiy*phixy - 2*phix*phiz*phixz &
                      - 2*phiy*phiz*phiyz) / ((phix**2 + phiy**2 + phiz**2)**(3./2.))

              else
                 kappa = ZERO
              endif


              speed = ( (HALF*( uadv(i,j,k) + uadv(i+1,j,k) ) )**2 &
                   + ( HALF*( vadv(i,j,k) +vadv(i,j+1,k) ) )**2 &
                   + ( HALF*( wadv(i,j,k) +wadv(i,j,k+1) ) )**2)**.5 + abs(kapa-kapb*kappa)

              phidt = min( phit,1/(max(1./phidt,speed/(.8*min(dx(1),dx(2),dx(3))), &
                   6*abs(kapb)/(.8*min(dx(1),dx(2),dx(3))**2) ))  )     

           endif
        enddo

        lscfl = phidt

      end subroutine ls_cfl

      


      subroutine ls_findinterface(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                                  phin,phin_l1,phin_l2,phin_l3,phin_h1,phin_h2,phin_h3, & 
                                  type,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                                  lo,hi,dx,intfacenump,intfacenumn,intfacep,intfacen, &
                                  nband,nbandsize,intfacesize)

        use LS_probdata_module, only: LARGEINT, BOGUS

        implicit none

        integer  phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
        integer  phin_l1,phin_l2,phin_l3,phin_h1,phin_h2,phin_h3
        integer  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

        double precision ::  phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
        double precision :: phin(phin_l1:phin_h1,phin_l2:phin_h2,phin_l3:phin_h3)
        integer          :: type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

        integer    lo(3), hi(3)    
        double precision     dx(3)
        integer    intfacenump, intfacenumn, intfacesize
        integer    intfacep(intfacesize,3), intfacen(intfacesize,3)
        integer    nbandsize
        integer    nband(nbandsize,3)


        !     Local variables
        INTEGER    i, j, k, r

        intfacenump=0
        intfacenumn=0

        r = 1
        do while(nband(r,1) .GT. -LARGEINT)
           i = nband(r,1)
           j = nband(r,2)
           k = nband(r,3)
           r = r + 1      

           phin(i,j,k) = sign(BOGUS,phi(i,j,k))
        enddo


        r = 1
        do while(nband(r,1) .GT. -LARGEINT)
           i = nband(r,1)
           j = nband(r,2)
           k = nband(r,3)
           r = r + 1

           call UPDATEF(i, j, k, &
                        phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                        phin, phin_l1, phin_l2, phin_l3, phin_h1, phin_h2, phin_h3, &
                        type, type_l1, type_l2, type_l3, type_h1, type_h2, type_h3, &
                        intfacep, intfacen, intfacenump, intfacenumn, &
                        dx, intfacesize, lo, hi) 

           if( (i .EQ. lo(1)) .OR. (j .EQ. lo(2)) .OR. k .EQ. lo(3)) then

              if(i .EQ. lo(1)) then
                 call UPDATEF(i-1, j, k, &
                              phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                              phin, phin_l1, phin_l2, phin_l3, phin_h1, phin_h2, phin_h3, &
                              type, type_l1, type_l2, type_l3, type_h1, type_h2, type_h3, &
                              intfacep, intfacen, intfacenump, intfacenumn, &
                              dx, intfacesize, lo, hi )
              endif

              if(j .EQ. lo(2)) then 
                 call UPDATEF(i, j-1, k, &
                              phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                              phin, phin_l1, phin_l2, phin_l3, phin_h1, phin_h2, phin_h3, &
                              type, type_l1, type_l2, type_l3, type_h1, type_h2, type_h3, &
                              intfacep, intfacen, intfacenump, intfacenumn, &
                              dx, intfacesize, lo, hi )
              endif

              if(k .EQ. lo(3)) then 
                 call UPDATEF(i, j, k-1, &
                              phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                              phin, phin_l1, phin_l2, phin_l3, phin_h1, phin_h2, phin_h3, &
                              type, type_l1, type_l2, type_l3, type_h1, type_h2, type_h3, &
                              intfacep, intfacen, intfacenump, intfacenumn, &
                              dx, intfacesize, lo, hi )
              endif

              if (i .EQ. lo(1) .AND. j .EQ. lo(2)) then
                 call UPDATEF(i-1, j-1, k, &
                              phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                              phin, phin_l1, phin_l2, phin_l3, phin_h1, phin_h2, phin_h3, &
                              type, type_l1, type_l2, type_l3, type_h1, type_h2, type_h3, &
                              intfacep, intfacen, intfacenump, intfacenumn, &
                              dx, intfacesize, lo, hi )
              endif

              if (i .EQ. lo(1) .AND. k .EQ. lo(3)) then
                 call UPDATEF(i-1, j, k-1, &
                              phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                              phin, phin_l1, phin_l2, phin_l3, phin_h1, phin_h2, phin_h3, &
                              type, type_l1, type_l2, type_l3, type_h1, type_h2, type_h3, &
                              intfacep, intfacen, intfacenump, intfacenumn, &
                              dx, intfacesize, lo, hi )
              endif

              if (j .EQ. lo(2) .AND. k .EQ. lo(3)) then
                 call UPDATEF(i, j-1, k-1, &
                              phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                              phin, phin_l1, phin_l2, phin_l3, phin_h1, phin_h2, phin_h3, &
                              type, type_l1, type_l2, type_l3, type_h1, type_h2, type_h3, &
                              intfacep, intfacen, intfacenump, intfacenumn, &
                              dx, intfacesize, lo, hi )
              endif


              if (i .EQ. lo(1) .AND. j .EQ. lo(2) .AND. k .EQ. lo(3)) then
                 call UPDATEF(i-1, j-1, k-1, &
                              phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                              phin, phin_l1, phin_l2, phin_l3, phin_h1, phin_h2, phin_h3, &
                              type, type_l1, type_l2, type_l3, type_h1, type_h2, type_h3, &
                              intfacep, intfacen, intfacenump, intfacenumn, &
                              dx, intfacesize, lo, hi )
              endif

           endif
        enddo

        nband(1,1) = -LARGEINT
        nband(1,2) = -LARGEINT
        nband(1,3) = -LARGEINT

      end subroutine ls_findinterface
      

      subroutine UPDATEF(i,j,k, &
                         phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                         phin, phin_l1, phin_l2, phin_l3, phin_h1, phin_h2, phin_h3, &
                         type, type_l1, type_l2, type_l3, type_h1, type_h2, type_h3, &
                         intfacep, intfacen, intfacenump, intfacenumn, &
                         dx, intfacesize, lo, hi  )
      
      use LS_probdata_module, only: BOGUS
      use bl_constants_module

      implicit none      

      integer  phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer  phin_l1,phin_l2,phin_l3,phin_h1,phin_h2,phin_h3
      integer  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

      double precision :: phi (phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
      double precision :: phin(phin_l1:phin_h1,phin_l2:phin_h2,phin_l3:phin_h3)
      integer          :: type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer          intfacenump, intfacenumn, intfacesize
      integer          intfacep(intfacesize,3), intfacen(intfacesize,3)
      double precision dx(3)
      integer          lo(3), hi(3)
      
!     Local variables
      INTEGER    c,d,e,i, j, k, ii, jj, kk, iii, jjj, kkk,l, m,n,p,q,r,s,t
      double precision     A(64,64)
      double precision     B(64)
      double precision     x,y,z
      double precision     FINDDIST, distance
      double precision     grad(3)
      integer              max0
      integer              nrow(64)
      double precision     geppX(64)
      
      if(max( abs(phi(i+1,j,k)),abs(phi(i,j+1,k)),abs(phi(i+1,j+1,k)),abs(phi(i,j,k+1)), &
           abs(phi(i+1,j+1,k+1)), abs(phi(i+1,j,k+1)), abs(phi(i,j+1,k+1)) ) .LT. BOGUS) then

     
         if ( phi(i,j,k)*phi(i+1,j+1,k) .LT. 0 .OR. phi(i,j,k)*phi(i+1,j,k) .LT. 0 &
              .OR. phi(i,j,k)*phi(i,j+1,k) .LT. 0 .OR. phi(i,j,k)*phi(i,j,k+1) .LT. 0 &
              .OR. phi(i,j,k)*phi(i,j+1,k+1) .LT. 0 .OR. phi(i,j,k)*phi(i+1,j,k+1) .LT. 0 &
              .OR. phi(i,j,k)*phi(i+1,j+1,k+1) .LT. 0 ) then

            m=0
            do ii=0,1
               do jj=0,1
                  do kk =0,1
                     x = ii*dx(1)
                     y = jj*dx(2)
                     z = kk*dx(3)
                     do n = 0,63
		        c = n/16
	         	d = n/4 - 4*(n/16)
	         	e = n   - 4*(n/4)
	 	        A(m+1,n+1) = x**c * y**d * z**e
	 	        A(m+2,n+1) = c * x**(max(c-1,0)) * y**d * z**e
	 	        A(m+3,n+1) = d * x**c * y**(max(d-1,0)) * z**e
	 	        A(m+4,n+1) = e * x**c * y**d * z**(max(e-1,0))
	 	        A(m+5,n+1) = c * d * x**(max(c-1,0)) * y**(max(d-1,0)) * z**e
                        A(m+6,n+1) = c * e * x**(max(c-1,0)) * y**d * z**(max(e-1,0))  
	 	        A(m+7,n+1) = d * e * x**c * y**(max(d-1,0)) * z**(max(e-1,0))  
	 	        A(m+8,n+1) = c * d * e &
                             * x**(max(c-1,0)) * y**(max(d-1,0)) * z**(max(e-1,0)) 
                     enddo
     
	     	     B(m+1) = phi(i+ii,j+jj,k+kk)
	     	     B(m+2) = (phi(i+ii+1,j+jj,k+kk) - phi(i+ii-1,j+jj,k+kk))/(2*dx(1))
	     	     B(m+3) = (phi(i+ii,j+jj+1,k+kk) - phi(i+ii,j+jj-1,k+kk))/(2*dx(2))
	     	     B(m+4) = (phi(i+ii,j+jj,k+kk+1) - phi(i+ii,j+jj,k+kk-1))/(2*dx(3))
	     	     B(m+5) = (phi(i+ii+1,j+jj+1,k+kk) - phi(i+ii-1,j+jj+1,k+kk) &
                          - phi(i+ii+1,j+jj-1,k+kk) &
                          + phi(i+ii-1,j+jj-1,k+kk))/(4*dx(1)*dx(2))
	     	     B(m+6) = (phi(i+ii+1,j+jj,k+kk+1) - phi(i+ii-1,j+jj,k+kk+1) &
                          - phi(i+ii+1,j+jj,k+kk-1) &
                          + phi(i+ii-1,j+jj,k+kk-1))/(4*dx(1)*dx(3))
	     	     B(m+7) = (phi(i+ii,j+jj+1,k+kk+1) - phi(i+ii,j+jj-1,k+kk+1) &
                          - phi(i+ii,j+jj+1,k+kk-1) &
                          + phi(i+ii,j+jj-1,k+kk-1))/(4*dx(2)*dx(3))	     	     
	     	     B(m+8) = (phi(i+ii+1,j+jj+1,k+kk+1) - phi(i+ii+1,j+jj+1,k+kk-1) &
                          + phi(i+ii+1,j+jj-1,k+kk-1) - phi(i+ii+1,j+jj-1,k+kk+1) &
                          - phi(i+ii-1,j+jj+1,k+kk+1) + phi(i+ii-1,j+jj+1,k+kk-1) &
                          - phi(i+ii-1,j+jj-1,k+kk-1) + phi(i+ii-1,j+jj-1,k+kk+1)) &
                          /(8*dx(1)*dx(2)*dx(3))
	     	     m = m + 8
	     	  enddo
               enddo
            enddo
            

            CALL GEPP(A,B,64,geppX,nrow)

            
            do ii=0,1
               do jj=0,1
                  do kk=0,1
                     
                     iii = i + ii
                     jjj = j + jj
                     kkk = k + kk
                  
                     distance = FINDDIST(grad,B, &
                          sign(ONE,phi(iii,jjj,kkk)),ii*dx(1),jj*dx(2),kk*dx(3),dx)     


                     if (type(iii,jjj,kkk) .NE. 0 .AND. phi(iii,jjj,kkk) .GE. 0 &
                          .AND. distance .GE. 0 &
                          .AND. iii .GE. lo(1) .AND. iii .LE. hi(1) &
                          .AND. jjj .GE. lo(2) .AND. jjj .LE. hi(2) &
     		          .AND. kkk .GE. lo(3) .AND. kkk .LE. hi(3)) then
                       
                        intfacenump = intfacenump + 1
                        intfacep(intfacenump,1)=iii
                        intfacep(intfacenump,2)=jjj
                        intfacep(intfacenump,3)=kkk
                     
                     
                     else if (type(iii,jjj,kkk) .NE. 0 .AND. distance .GE. 0 &
                          .AND. iii .GE. lo(1) .AND. iii .LE. hi(1) &
                          .AND. jjj .GE. lo(2) .AND. jjj .LE. hi(2) &
                          .AND. kkk .GE. lo(3) .AND. kkk .LE. hi(3))  then
                     
                        intfacenumn = intfacenumn + 1
                        intfacen(intfacenumn,1)=iii
                        intfacen(intfacenumn,2)=jjj
                        intfacen(intfacenumn,3)=kkk
                     
                     endif
                  
                     if (distance .GE. 0) then
                        type(iii,jjj,kkk) = 0
                        phin(iii,jjj,kkk) = min(abs(phin(iii,jjj,kkk)),distance) &
                             *sign(ONE,phi(iii,jjj,kkk))
                     endif
                  enddo   
               enddo
            enddo
            
         endif
      endif
      return

    end subroutine UPDATEF

      double precision function FINDDIST(grad,B,sgn,x0,y0,z0,dx)

      use LS_probdata_module, only: BOGUS

      implicit none
      double precision :: grad(3)
      double precision :: B(64)
      double precision :: sgn
      double precision :: x0,y0,z0
      double precision :: dx(3)

!     Local variables
      double precision t,tp
      double precision POLYVAL
      double precision DPOLYVAL
      double precision FA,FP
      double precision a,d,p
      double precision x,y,z
      double precision delta1(3), delta2(3)

      INTEGER i    
      integer ITERMAX
      parameter (ITERMAX=30)
      
      x = x0
      y = y0 
      z = z0
      i = 0
      
      delta1(1) = BOGUS
      delta1(2) = BOGUS
      delta1(3) = BOGUS
      delta2(1) = BOGUS
      delta2(2) = BOGUS
      delta2(3) = BOGUS
      
      do while  ( (delta1(1)**2 + delta1(2)**2 + delta1(3)**2 &
           + delta2(1)**2 + delta2(2)**2 + delta2(3)**2)**(.5) &
           .GT. 10.0**(-5.0)*dx(1)*dx(2)*dx(3) .AND. i .LT. ITERMAX)
        
        CALL GRADPVAL(B,grad,x,y,z)
        
      	delta1(1) = -polyval(B,x,y,z) * grad(1)/(grad(1)**2 + grad(2)**2 + grad(3)**2)
      	delta1(2) = -polyval(B,x,y,z) * grad(2)/(grad(1)**2 + grad(2)**2 + grad(3)**2)
      	delta1(3) = -polyval(B,x,y,z) * grad(3)/(grad(1)**2 + grad(2)**2 + grad(3)**2)
      	
      	delta2(1) = (x0 - x) - grad(1)*( (x0 - x)*grad(1) + (y0 - y)*grad(2) &
             + (z0-z)*grad(3) ) /(grad(1)**2 + grad(2)**2 + grad(3)**2)
      	delta2(2) = (y0 - y) - grad(2)*( (x0 - x)*grad(1) + (y0 - y)*grad(2) &
             + (z0-z)*grad(3) ) /(grad(1)**2 + grad(2)**2 + grad(3)**2) 
      	delta2(3) = (z0 - z) - grad(3)*( (x0 - x)*grad(1) + (y0 - y)*grad(2) &
             + (z0-z)*grad(3) ) /(grad(1)**2 + grad(2)**2 + grad(3)**2)       	
      	
      	x = x + delta1(1) + delta2(1)
      	y = y + delta1(2) + delta2(2)
      	z = z + delta1(3) + delta2(3)
      	
      	i = i + 1
      
      enddo
      
      
      if (i .GE. 30 .OR. x .LT. 0 .OR. x .GT. dx(1) .OR. y .LT. 0 .OR. y .GT. dx(2) &
           .OR. z .LT. 0 .OR. z .GT. dx(3) ) then
      
      	FINDDIST = -1
      	
      else
      
      	FINDDIST = ( (x - x0)**2 + (y-y0)**2 +(z-z0)**2 )**(.5)

      endif

    end function FINDDIST


    double precision function POLYVAL(B,x,y,z)

      use bl_constants_module

      implicit none

      double precision B(64)
      double precision x,y,z
      integer c,d,e,n  
      POLYVAL=ZERO
      do n=0,63
         c = n/16
         d = n/4 - 4*(n/16)
         e = n   - 4*(n/4)
         POLYVAL = POLYVAL + B(n+1) * x**c * y**d * z**e
      enddo
    end function POLYVAL
    
      
      
    subroutine GRADPVAL(B,grad,x,y,z)
      implicit none
      double precision B(64)
      double precision grad(3)
      double precision x,y,z      
      
!     Local variables
      INTEGER c,d,e,n
      
      grad(1) = 0
      grad(2) = 0
      grad(3) = 0
      
      do n = 0,63
      
        c = n/16
	d = n/4 - 4*(n/16)
	e = n   - 4*(n/4)
      
      	grad(1) = grad(1) + B(n+1) * c * x**(max(c-1,0)) * y**d * z**e
      	grad(2) = grad(2) + B(n+1) * d * x**c * y**(max(d-1,0)) * z**e
      	grad(3) = grad(3) + B(n+1) * e * x**c * y**d * z**(max(e-1,0))
      
      enddo
      return
    end subroutine GRADPVAL
      
      subroutine GEPP(A,B,N,X,nrow)
      implicit none
      INTEGER N
      double precision  A(N,N)
      double precision  B(N)
      INTEGER i,j,k,p
      INTEGER ncopy
      INTEGER nrow(N)
      double precision  maximum
      double precision  X(N)
      double precision  m
      
      do i = 1,N
        nrow(i) = i
      enddo
      
      do i = 1,N-1
         maximum = abs(A(nrow(i),i))
         p=i
         do j = i+1, N
            if (maximum < abs(A(nrow(j),i))) then
               maximum = abs(A(nrow(j),i))
               p=j
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
      
      
      subroutine ls_narrowband(type,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, & 
      		               nband,nbandsize,mine,minesize,lo,hi)

      use LS_probdata_module, only: LARGEINT
        
      implicit none     
     
      integer type_l1,type_l2,type_l3,type_h1,type_h2,type_h3
      integer type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer nbandsize, minesize
      integer lo(3), hi(3)
      integer nband(nbandsize,3)
      integer  mine(minesize,3)


      integer numband
      integer nummine
      integer i,j,k
     
      numband = 0
      nummine = 0
     
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)            
             if(type(i,j,k) .eq. 0  ) then
                numband = numband + 1
                nband(numband,1) = i
                nband(numband,2) = j
                nband(numband,3) = k
             else if (type(i,j,k) .eq. 1) then
                numband = numband + 1
                nband(numband,1) = i
                nband(numband,2) = j          
                nband(numband,3) = k
                
                nummine = nummine + 1
                mine(nummine,1) = i
                mine(nummine,2) = j
                mine(nummine,3) = k
             endif
          enddo
        enddo
      enddo
     
      nband(numband + 1,1)  = -LARGEINT
      nband(numband + 1,2)  = -LARGEINT
      nband(numband + 1,3)  = -LARGEINT      
      
    
      mine(nummine + 1,1) = -LARGEINT
      mine(nummine + 1,2) = -LARGEINT
      mine(nummine + 1,3) = -LARGEINT      
      
    end subroutine ls_narrowband
     
 
    subroutine ls_retypify(type,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                           nband,nbandsize)     

      use LS_probdata_module, only: LARGEINT

      implicit none  

      integer type_l1,type_l2,type_l3,type_h1,type_h2,type_h3
      integer type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer nbandsize
      integer nband(nbandsize,3)

      integer i, j, k, p

      p = 1
      do while (nband(p,1) .GT. -LARGEINT)
        i = nband(p,1)
        j = nband(p,2)
        k = nband(p,3)
        p = p + 1
        type(i,j,k) = 3
      enddo

    end subroutine ls_retypify
      
     
    subroutine ls_fastmarch(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, & 
                            type,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                            lo, hi, dx, intfacenum, intface,  &
                            nband, nbandsize, nbandnum, mine, &
                            sgn, intfacesize,heap, heaploc)

      use LS_probdata_module, only: nbandwidth, mineloc, BOGUS, LARGEINT

      implicit none

      integer phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
      integer          :: type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer heaploc(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer     lo(3), hi(3)
      double precision      dx(3)
      integer     intfacenum, intfacesize
      integer     intface(intfacesize,3)
      integer     nbandsize
      integer     nband(nbandsize,3)
      integer     mine(nbandsize,3)
      integer     nbandnum
      integer     sgn
      integer heap(nbandsize,3)

      integer i, j, k, n, p
      integer numtent
      
      numtent = 0
      do n = 1, intfacenum
        i = intface(n,1)
        j = intface(n,2)
        k = intface(n,3)
      
	CALL UPDATE(phi,i,j,k,sgn,type,heap,numtent, &
                    phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                    type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                    nbandsize, lo, hi, dx, heaploc)
	
	nbandnum = nbandnum + 1
	nband(nbandnum,1) = i
	nband(nbandnum,2) = j
	nband(nbandnum,3) = k
      enddo
      
      do while (numtent .GT. 0  )

         CALL RMVNODE(heap,i,j,k,numtent,phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                      nbandsize,heaploc,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
     
         if (abs(phi(i,j,k)) .LT. nbandwidth) then
            nbandnum = nbandnum + 1
            nband(nbandnum,1) = i
            nband(nbandnum,2) = j
            nband(nbandnum,3) = k
            type(i,j,k) = 0
         else
            type(i,j,k) = 3
            phi(i,j,k) = sign(BOGUS,phi(i,j,k))
            exit
         endif
         
         if (abs(phi(i,j,k)) .GT. mineloc .AND. abs(phi(i,j,k)) .LT. nbandwidth ) then
          
            type(i,j,k) = 1
            
         endif
         
         CALL UPDATE(phi,i,j,k,sgn, type, heap, numtent, &
                     phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                     type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                     nbandsize, lo, hi, dx,heaploc)
         
      enddo

      nband(nbandnum+1,1) = -LARGEINT
      nband(nbandnum+1,2) = -LARGEINT      
      nband(nbandnum+1,3) = -LARGEINT
      
      do while (numtent .GT. 0)

         CALL RMVNODE(heap, i, j, k,numtent, &
                      phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                      nbandsize,heaploc,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
         type(i,j,k) = 3
         phi(i,j,k) = sign(BOGUS,phi(i,j,k))
         
      enddo

    end subroutine ls_fastmarch
      
      subroutine UPDATE(phi, i, j, k, sgn, type, heap, numtent, &
                        phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                        type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                        nbandsize, lo, hi, dx, heaploc)
      
      implicit none

      integer phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
      integer          :: type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)
      integer          :: heaploc(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer i,j,k
      integer nbandsize
      integer heap(nbandsize,3)
      integer lo(3), hi(3)
      integer numtent
      integer sgn
      integer n,ii,jj,kk
      double precision dx(3)
      
      do n = 1,6
         
         ii = i + (2*n-3)*(1-(n/3)+(n/6))
         jj = j + (2*n-7)*(n/3-(n/5)-(n/6))
         kk = k + (2*n-11)*(n/5)
         
         if ( ii .GE. lo(1) .AND. ii .LE. hi(1) .AND. jj .GE. lo(2) .AND. jj .LE. hi(2) &
              .AND. kk .GE. lo(3) .AND. kk .LE. hi(3)) then
            
            if (type(ii,jj,kk) .GT. 1   .AND.   sgn*phi(ii,jj,kk) .GE. 0) then
               
               CALL EVAL(phi, ii, jj, kk, phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                         type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                         lo, hi, type, sgn, dx)
               
               if (type(ii,jj,kk) .GT. 2) then
                  
                  type(ii,jj,kk) = 2
                  CALL ADDNODE(heap, ii, jj, kk, numtent, &
                       phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, lo, hi, &
                       nbandsize,heaploc,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
                  
               else
                  
                  CALL UPDATENODE(heap, ii, jj, kk, numtent, &
                                  phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,lo,hi, &
                                  nbandsize,heaploc, &
                                  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
                  
               endif
            endif    
         endif
      enddo
      end
            

      subroutine EVAL(phi,i,j,k,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                      type_l1,type_l2,type_l3,type_h1,type_h2,type_h3,lo,hi,type,sgn,dx)

      use bl_constants_module      

      implicit none

      integer i,j,k

      integer  phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
      integer          :: type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer lo(3), hi(3)
      double precision a,b,c
      double precision dx(3)

      integer sgn
      integer  left,right,up,down,front,back
      LOGICAL  lok, rok, uok, dok, fok, bok
      
      a = ZERO
      b = ZERO
      c = - 1
      
      left  = i - 1
      right = i + 1
      up    = j + 1
      down  = j - 1
      front = k + 1
      back  = k - 1
      
      lok  = left .GE. lo(1) .AND. sgn*phi(left,j,k) .GE. 0 .AND. type(left,j,k) .LE. 1
      rok  = right .LE. hi(1) .AND. sgn*phi(right,j,k) .GE. 0  .AND. type(right,j,k) .LE. 1
      uok  = up .LE. hi(2) .AND. sgn*phi(i,up,k) .GE. 0  .AND. type(i,up,k) .LE. 1
      dok  = down  .GE. lo(2) .AND. sgn*phi(i,down,k) .GE. 0  .AND. type(i,down,k) .LE. 1
      fok  = front .LE. hi(3) .AND. sgn*phi(i,j,front) .GE. 0  .AND. type(i,j,front) .LE. 1
      bok  = back  .GE. lo(3) .AND. sgn*phi(i,j,back) .GE. 0  .AND. type(i,j,back) .LE. 1

!     FIXME: The following had a sign(right,...), mistake?
      if (lok .AND. rok) then
                    
        a = a + 1/dx(1)**2
        b = b + min(sgn*phi(left,j,k),sgn*phi(right,j,k))/(dx(1)**2)
        c = c + min(phi(left,j,k)**2,phi(right,j,k)**2)/(dx(1)**2)
        
      else if (lok) then
      
        a = a + 1/dx(1)**2
	b = b + sgn*phi(left,j,k)/(dx(1)**2)
        c = c + phi(left,j,k)**2/(dx(1)**2)
      
      else if (rok ) then

        a = a + 1/dx(1)**2
	b = b + sgn*phi(right,j,k)/(dx(1)**2)
        c = c + phi(right,j,k)**2/(dx(1)**2)
        
      endif
      
      
      if (dok .AND. uok) then
      
        a = a + 1/dx(2)**2
        b = b + min(sgn*phi(i,down,k),sgn*phi(i,up,k))/(dx(2)**2)
        c = c + min(phi(i,down,k)**2,phi(i,up,k)**2)/(dx(2)**2)
        
      else if (dok) then
      
        a = a + 1/dx(2)**2
	b = b + sgn*phi(i,down,k)/(dx(2)**2)
        c = c + phi(i,down,k)**2/(dx(2)**2)
      
      else if (uok ) then

        a = a + 1/dx(2)**2
	b = b + sgn*phi(i,up,k)/(dx(2)**2)
        c = c + phi(i,up,k)**2/(dx(2)**2)
        
      endif   
      
      if (fok .AND. bok) then
      
        a = a + 1/dx(3)**2
        b = b + min(sgn*phi(i,j,front),sgn*phi(i,j,back))/(dx(3)**2)
        c = c + min(phi(i,j,front)**2,phi(i,j,back)**2)/(dx(3)**2)
        
      else if (fok) then
      
        a = a + 1/dx(3)**2
	b = b + sgn*phi(i,j,front)/(dx(3)**2)
        c = c + phi(i,j,front)**2/(dx(3)**2)
      
      else if (bok ) then

        a = a + 1/dx(3)**2
	b = b + sgn*phi(i,j,back)/(dx(3)**2)
        c = c + phi(i,j,back)**2/(dx(3)**2)
        
      endif       
        
      b = -2*b
      
      if (b**2 - 4*a*c .LT. 0) then
      
         phi(i,j,k) = sgn*(-b)/(2*a)
         return
      
      endif
      phi(i,j,k) = sgn*(-b + SQRT(b**2-4*a*c))/(2*a)

    end subroutine EVAL
 

      subroutine ADDNODE(heap,i,j,k,n,phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,lo,hi, &
                         nbandsize, heaploc,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
      implicit none

      integer phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

      double precision phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
      integer heaploc(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer nbandsize
      integer lo(3), hi(3)
      integer heap(nbandsize,3)
      integer i,j,k,n

      integer index
      integer parent
      
      index = n + 1
      parent = index/2

      if (n .EQ. 0) then
         
         heap(index,1) = i
         heap(index,2) = j
         heap(index,3) = k
         heaploc(i,j,k) = index        
         
         n=n+1
         
         return
      endif
      
      do while ( ABS(phi(heap(parent,1),heap(parent,2),heap(parent,3))) .GT. ABS(phi(i,j,k)))
            
        heap(index,1) = heap(parent,1)
        heap(index,2) = heap(parent,2)
        heap(index,3) = heap(parent,3)
        heaploc(heap(index,1),heap(index,2),heap(index,3)) = index        
        index = parent
        parent = index/2
         
        if (parent .EQ. 0) exit
      enddo
      
      heap(index,1) = i
      heap(index,2) = j
      heap(index,3) = k
      heaploc(i,j,k) = index
      n = n + 1

    end subroutine ADDNODE
      

      subroutine UPDATENODE(heap,i,j,k,n, &
                            phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,lo,hi,nbandsize, &
                            heaploc,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
      implicit none

      integer  phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
      integer          :: heaploc(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer    nbandsize
      integer    lo(3), hi(3)
      integer    heap(nbandsize,3)
      integer    i,j,k, n

      integer index
      integer parent
      
      index = heaploc(i,j,k)
      parent = index/2

      if (index .EQ. 1) then
        return
      endif
      
      do while ( ABS(phi(heap(parent,1 ), heap(parent,2), heap(parent,3) )) &
           .GT. ABS(phi(i,j,k)))
            
        heap(index,1) = heap(parent,1)
        heap(index,2) = heap(parent,2)
        heap(index,3) = heap(parent,3)
        heaploc(heap(index,1),heap(index,2),heap(index,3)) = index
        index = parent
        parent = index/2
         
        if (parent .EQ. 0) exit         
      enddo
      
      heap(index,1) = i
      heap(index,2) = j
      heap(index,3) = k
      heaploc(i,j,k) = index
    end subroutine UPDATENODE
      

      subroutine RMVNODE(heap,i,j,k,n,phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                         nbandsize, heaploc,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
      implicit none

      integer  phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
      integer          :: heaploc(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer    nbandsize
      integer    heap(nbandsize,3)
      integer    i,j,k, n
      integer    index, left, right

      i=heap(1,1)
      j=heap(1,2)
      k=heap(1,3)
      
      heaploc(i,j,k) = -1
      
      index = 1
      left  = 2*index
      right = 2*index + 1
   
      do while (.TRUE.)
         
         if (left .LE. n-1) then
            
            if( ABS(phi(heap(left,1),heap(left,2),heap(left,3))) .LT. &
                 ABS(phi(heap(n,1),heap(n,2),heap(n,3)))) then
               
               if (right .LE. n-1) then
                  if (ABS(phi(heap(left,1),heap(left,2),heap(left,3))) .LT. &
                       ABS(phi(heap(right,1),heap(right,2),heap(right,3))) ) then
                     heap(index,1) = heap(left,1)
                     heap(index,2) = heap(left,2)
                     heap(index,3) = heap(left,3)
                     heaploc(heap(index,1),heap(index,2),heap(index,3)) = index
                     index = left
                     left=2*index
                     right=2*index+1
                  else
                     heap(index,1) = heap(right,1)
                     heap(index,2) = heap(right,2)
                     heap(index,3) = heap(right,3)
                     heaploc(heap(index,1),heap(index,2),heap(index,3)) = index	      
                     index = right
                     left=2*index
                     right=2*index+1  
                  endif
               else                  
                  heap(index,1) = heap(left,1)
                  heap(index,2) = heap(left,2)
                  heap(index,3) = heap(left,3)
                  heaploc(heap(index,1),heap(index,2),heap(index,3)) = index	    
                  
                  index = left
                  left=2*index
                  right=2*index+1  
               endif     
               
            else if (right .LE. n-1) then
         
               if (abs(phi(heap(right,1),heap(right,2),heap(right,3))) &
                    .LT. abs( phi(heap(n,1),heap(n,2),heap(n,3)))) then
                  heap(index,1) = heap(right,1)
                  heap(index,2) = heap(right,2)
                  heap(index,3) = heap(right,3)
                  heaploc(heap(index,1),heap(index,2),heap(index,3)) = index	      
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
      heap(index,3) = heap(n,3)
      if(n .GT. 1) then
        heaploc(heap(index,1),heap(index,2),heap(index,3)) = index
      endif
      n = n - 1
    end subroutine RMVNODE

      subroutine ls_fastmarch2(flag,phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                               type,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                               lo, hi, dx, nband, nbandsize, nbandnum, &
                               sgn, heaploc)

      use LS_probdata_module, only: nbandwidth, mineloc, LARGEINT, BOGUS

      implicit none

      integer  flag
      integer  phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

      double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
      integer          :: type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)
      integer          :: heaploc(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer     lo(3), hi(3)
      double precision      dx(3)
      integer     nbandsize
      integer     nband(nbandsize,3)
      integer     nbandnum
      integer     sgn
      
      integer i, j, k, n
      integer heap(nbandsize,3)
      integer numtent
      
      numtent = 0
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             heaploc(i,j,k) = -1
          enddo
        enddo
      enddo
     
      i = lo(1)-1 
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent, &
                           phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                           type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                           nbandsize, lo, hi, dx,heaploc)             
           endif
        enddo 
      enddo
      
      i = hi(1)+1 
      do k = lo(3), hi(3)
        do j = lo(2), hi(2)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent, &
                           phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                           type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                           nbandsize, lo, hi, dx,heaploc)    
           endif
        enddo 
      enddo   
      
      
      j = lo(2)-1 
      do k = lo(3), hi(3)
        do i = lo(1), hi(1)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0)  then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent, &
                           phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                           type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                           nbandsize, lo, hi, dx,heaploc)   
           endif
        enddo 
      enddo
      
      
      j = hi(2)+1       
      do k = lo(3), hi(3)
        do i = lo(1), hi(1)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0 ) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent, &
                           phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                           type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                           nbandsize, lo, hi, dx,heaploc)   
           endif
        enddo 
      enddo
      
      k = lo(3)-1       
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent, &
                           phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                           type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                           nbandsize, lo, hi, dx,heaploc)  
           endif
        enddo 
      enddo      

      k = hi(3) + 1       
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)      
           if((type(i,j,k) .EQ. 0 .OR. type(i,j,k) .EQ. 1) .AND. sgn*phi(i,j,k) .GE. 0) then
              CALL UPDATE2(phi,i,j,k,sgn, type,heap,numtent, &
                           phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                           type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                           nbandsize, lo, hi, dx,heaploc)  
           endif
        enddo 
      enddo      
      
      flag = 0
      do while (numtent .GT. 0 )
         CALL RMVNODE(heap, i,j,k,numtent, phi, phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                      nbandsize,heaploc,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
         if (abs(phi(i,j,k)) .LT. nbandwidth ) then
         
            flag = 1

            if (type(i,j,k) .GE. 2) then
               nbandnum = nbandnum + 1         
               nband(nbandnum,1) = i
               nband(nbandnum,2) = j            
	       nband(nbandnum,3) = k                           
            endif
          
            type(i,j,k) = 0
            CALL UPDATE2(phi,i,j,k,sgn, type, heap,numtent, &
                         phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                         type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                         nbandsize, lo, hi, dx,heaploc)

         else

            type(i,j,k) = 3
            phi(i,j,k) = sign(BOGUS,phi(i,j,k))
            exit

         endif
         
         if (abs(phi(i,j,k)) .GT. mineloc .AND. abs(phi(i,j,k)) .LT. nbandwidth) then
          
            type(i,j,k) = 1
        endif
      enddo
      
      nband(nbandnum + 1,1) = -LARGEINT
      nband(nbandnum + 1,2) = -LARGEINT
      nband(nbandnum + 1,3) = -LARGEINT
      
      do while (numtent .GT. 0)
            
         CALL RMVNODE(heap, i,j,k,numtent, phi, phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                      nbandsize,heaploc,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
         type(i,j,k) = 3
         phi(i,j,k) = sign(BOGUS,phi(i,j,k))    
      enddo
    end subroutine ls_fastmarch2
      
      
      
      subroutine UPDATE2(phi,i,j,k,sgn, &
                         type,heap,numtent,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,&
                         type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                         nbandsize, lo, hi, dx, heaploc)

      use LS_probdata_module, only: BOGUS
      use bl_constants_module

      implicit none

      integer  phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

      double precision :: phi (phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
      integer          :: type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)
      integer          :: heaploc(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer i,j,k
      integer nbandsize
      integer heap(nbandsize,3)
      integer lo(3), hi(3)
      integer numtent
      integer sgn
      integer n,ii,jj,kk
      double precision dx(3)

      logical EVAL2
      logical isnew
      integer min

      
      do n = 1,6
         
         ii = i + (2*n-3)*(1-(n/3)+(n/6))
         jj = j + (2*n-7)*(n/3-(n/5)-(n/6))
         kk = k + (2*n-11)*(n/5)
         
              
        if (ii .GE. lo(1) .AND. ii .LE. hi(1) .AND. jj .GE. lo(2) .AND. jj .LE. hi(2) &
             .AND. kk .GE. lo(3) .AND. kk .LE. hi(3)) then
        
           if (sgn*sign(ONE,phi(ii,jj,kk)) .GE. 0 .AND. &
                abs(phi(ii,jj,kk)) .GT. abs(phi(i,j,k)) &
                .AND. (abs(phi(ii,jj,kk)) .GE. BOGUS .OR. .NOT.(sgn*phi(ii+1,jj,kk) .LE. 0  &
                .OR. sgn*phi(ii-1,jj,kk) .LE. 0 &
                .OR. sgn*phi(ii,jj+1,kk) .LE. 0 .OR. sgn*phi(ii,jj-1,kk) .LE. 0 &
!                .OR. sgn*phi(ii+1,jj+1,kk) .LE. 0 &
!                .OR. sgn*phi(ii-1,jj-1,kk) .LE. 0 .OR. sgn*phi(ii+1,jj-1,kk) .LE. 0 &
!                .OR. sgn*phi(ii-1,jj+1,kk) .LE. 0 &
                .OR. sgn*phi(ii,jj,kk+1) .LE. 0 .OR. sgn*phi(ii,jj,kk-1) .LE. 0))) then
!                .OR. sgn*phi(ii+1,jj,kk+1) .LE. 0 &
!                .OR. sgn*phi(ii+1,jj,kk-1) .LE. 0 .OR. sgn*phi(ii-1,jj,kk+1) .LE. 0 &
!                .OR. sgn*phi(ii-1,jj,kk-1) .LE. 0 .OR. sgn*phi(ii,jj+1,kk+1) .LE. 0 &
!                .OR. sgn*phi(ii,jj+1,kk-1) .LE. 0 .OR. sgn*phi(ii,jj-1,kk+1) .LE. 0 &
!                .OR. sgn*phi(ii,jj-1,kk-1) .LE. 0 .OR. sgn*phi(ii+1,jj+1,kk+1) .LE. 0 &
!                .OR. sgn*phi(ii+1,jj+1,kk-1) .LE. 0 .OR. sgn*phi(ii+1,jj-1,kk+1) .LE. 0 &
!                .OR. sgn*phi(ii+1,jj-1,kk-1) .LE. 0 .OR. sgn*phi(ii-1,jj+1,kk+1) .LE. 0 &
!                .OR. sgn*phi(ii-1,jj+1,kk-1) .LE. 0 .OR. sgn*phi(ii-1,jj-1,kk+1) .LE. 0 &
!                .OR. sgn*phi(ii-1,jj-1,kk-1) .LE. 0))) then             
              isnew = EVAL2(phi,ii,jj,kk,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                   type_l1,type_l2,type_l3,type_h1,type_h2,type_h3,lo,hi, &
                   type,sgn,dx, phi(i,j,k))
           
              isnew = isnew .OR. EVAL2(phi,ii,jj,kk, &
                   phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                   type_l1,type_l2,type_l3,type_h1,type_h2,type_h3,lo,hi, &
                   type,sgn,dx, phi(ii,jj,kk))
              
              if (isnew ) then
                 type(ii,jj,kk) = min(2,type(ii,jj,kk))
                 if(heaploc(ii,jj,kk) .EQ. -1) then
                    CALL ADDNODE(heap,ii,jj, kk, numtent, &
                         phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,lo,hi,nbandsize, &
                         heaploc,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)
                 else
                    CALL UPDATENODE(heap,ii,jj, kk, numtent, &
                         phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,lo,hi,nbandsize, &
                         heaploc,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3)   
                 endif
              endif
           endif    
        endif        
      enddo
    end subroutine UPDATE2
      
    
      LOGICAL function EVAL2(phi,i,j,k,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                             type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                             lo,hi, type,sgn,dx,phisrc)      

      use bl_constants_module

      implicit none

      integer  phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer  type_l1,type_l2,type_l3,type_h1,type_h2,type_h3

      double precision :: phi (phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
      integer          :: type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer i,j,k
      integer lo(3), hi(3)
      double precision a,b,c
      double precision dx(3)

      integer sgn
      integer  left,right,up,down,front, back
      double precision phisrc
      LOGICAL  lok, rok, uok, dok, fok, bok
      
      a = ZERO
      b = ZERO
      c = - 1
      
      left  = i - 1
      right = i + 1
      up    = j + 1
      down  = j - 1
      front = k + 1
      back  = k - 1
      
      lok  = sgn*phi(left,j,k) .GE. 0 .AND. type(left,j,k) .LE. 1 &
           .AND. abs(phi(left,j,k)) .LE. abs(phisrc)
      rok  = sgn*phi(right,j,k) .GE. 0  .AND. type(right,j,k) .LE. 1 &
           .AND. abs(phi(right,j,k)) .LE. abs(phisrc)
      uok  = sgn*phi(i,up,k) .GE. 0  .AND. type(i,up,k) .LE. 1 &
           .AND. abs(phi(i,up,k)) .LE. abs(phisrc)
      dok  = sgn*phi(i,down,k) .GE. 0  .AND. type(i,down,k) .LE. 1 &
           .AND. abs(phi(i,down,k)) .LE. abs(phisrc)
      fok  = sgn*phi(i,j,front) .GE. 0  .AND. type(i,j,front) .LE. 1 &
           .AND. abs(phi(i,j,front)) .LE. abs(phisrc)
      bok  = sgn*phi(i,j,back) .GE. 0  .AND. type(i,j,back) .LE. 1 &
           .AND. abs(phi(i,j,back)) .LE. abs(phisrc)         

      if (lok .AND. rok ) then
         
         a = a + 1/dx(1)**2        
         b = b + min(sgn*phi(left,j,k),sgn*phi(right,j,k))/(dx(1)**2)
         c = c + min(phi(left,j,k)**2,phi(right,j,k)**2)/(dx(1)**2)
         
      else if (lok) then
         
         a = a + 1/dx(1)**2
         b = b + sgn*phi(left,j,k)/(dx(1)**2)
         c = c + phi(left,j,k)**2/(dx(1)**2)
         
      else if (rok ) then
         
         a = a + 1/dx(1)**2
         b = b + sgn*phi(right,j,k)/(dx(1)**2)
         c = c + phi(right,j,k)**2/(dx(1)**2)
         
      endif

      if ( uok .AND. dok ) then
         
         a = a + 1/dx(2)**2
         b = b + min(sgn*phi(i,down,k),sgn*phi(i,up,k))/(dx(2)**2)
         c = c + min(phi(i,down,k)**2,phi(i,up,k)**2)/(dx(2)**2)
         
      else if (dok) then
         
         a = a + 1/dx(2)**2
         b = b + sgn*phi(i,down,k)/(dx(2)**2)
         c = c + phi(i,down,k)**2/(dx(2)**2)
         
      else if (uok ) then
         
         a = a + 1/dx(2)**2
         b = b + sgn*phi(i,up,k)/(dx(2)**2)
         c = c + phi(i,up,k)**2/(dx(2)**2)
         
      endif
         
      if ( fok .AND. bok ) then
         
         a = a + 1/dx(3)**2
         b = b + min(sgn*phi(i,j,back),sgn*phi(i,j,front))/(dx(3)**2)
         c = c + min(phi(i,j,back)**2,phi(i,j,front)**2)/(dx(3)**2)
         
      else if (bok) then
         
         a = a + 1/dx(3)**2
         b = b + sgn*phi(i,j,back)/(dx(3)**2)
         c = c + phi(i,j,back)**2/(dx(3)**2)
         
      else if (fok ) then
         
         a = a + 1/dx(3)**2
         b = b + sgn*phi(i,j,front)/(dx(3)**2)
         c = c + phi(i,j,front)**2/(dx(3)**2)
         
      endif
      
      b = -2*b
      if (a .EQ. ZERO) then

        EVAL2 = .FALSE.
        return
       
      endif
      
      EVAL2 = .FALSE.
      
      if (b**2 - 4*a*c .LT. ZERO) then
         if (ABS(phi(i,j,k)) .GT. (-b)/(2*a) + 1.d-10) then
            phi(i,j,k) = sgn*(-b)/(2*a)
            EVAL2 = .TRUE.
            return
         endif
      endif
      
      if (abs(phi(i,j,k)) .GT. (-b+SQRT(b**2-4*a*c))/(2*a)+1.d-10) then
        phi(i,j,k) = sgn*(-b+SQRT(b**2-4*a*c))/(2*a)
        EVAL2 = .TRUE.
      endif
      
    end function EVAL2
      
      



     
      subroutine ls_mine(type,type_l1,type_l2,type_l3,type_h1,type_h2,type_h3, &
                         nband, nbandsize,mine, minesize,lo, hi)    

      use LS_probdata_module, only: LARGEINT

      implicit none     

      integer type_l1,type_l2,type_l3,type_h1,type_h2,type_h3
      integer type(type_l1:type_h1,type_l2:type_h2,type_l3:type_h3)

      integer nbandsize, minesize
      integer lo(3), hi(3)
      integer  nband(nbandsize,3)
      integer  mine(minesize,3)

      integer i, j,k, p
      integer nummine
      
      nummine = 0
      p = 1
      do while(nband(p,1) .GT. -LARGEINT)
      
        i = nband(p,1)
        j = nband(p,2)
        k = nband(p,3)
        p = p + 1
        
        if (type(i,j,k) .EQ. 1) then
          nummine = nummine +1
          mine(nummine,1) = i
          mine(nummine,2) = j
          mine(nummine,3) = k
        endif        
      enddo
      
      mine(nummine+1,1) = -LARGEINT
      mine(nummine+1,2) = -LARGEINT      
      mine(nummine+1,3) = -LARGEINT
      
    end subroutine ls_mine
           
     

     
    subroutine ls_nbandnumify(nband, nbandsize,nbandnum)

      use LS_probdata_module, only: LARGEINT

      implicit none      
      integer nbandsize, nbandnum
      integer nband(nbandsize,3)
      integer p
      
      p = 0
      do while(nband(p+1,1) .gt. -LARGEINT)
        p = p + 1
      enddo
      nbandnum = p
    end subroutine ls_nbandnumify
      

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
      subroutine ca_lserror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                            set,clear, &
                            ls,lsl1,lsl2,lsl3,lsh1,lsh2,lsh3, &
                            lo,hi,nd,domlo,domhi, &
                            delta,xlo,problo,time,level)
      use LS_probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer lsl1,lsl2,lsl3,lsh1,lsh2,lsh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision ls(lsl1:lsh1,lsl2:lsh2,lsl3:lsh3,nd)
      double precision delta(3), xlo(3), problo(3), time

      integer i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( abs(ls(i,j,k,1)) .lt. lvlerr ) then
                  tag(i,j,k) = set
               else
                  tag(i,j,k) = tag(i,j,k)
               end if
            end do
         end do
      end do
      end
