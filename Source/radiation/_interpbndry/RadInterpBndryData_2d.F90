#include <AMReX_BC_TYPES.H>
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <RADINTERPBNDRYDATA_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2
#define NUMDERIV 2
#define XDER   1
#define X2DER  2
      
      
! ---------------------------------------------------------------
! ::  FORT_BDINTERPXLO : Interpolation on Xlo Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  DIMS(bdry)  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  DIMS(cb)    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(2)   => refinement ratios
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  DIMS(mask)  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  DIMS(crse)  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

subroutine FORT_BDINTERPXLO (bdry,DIMS(bdry), &
                             lo,hi,DIMS(cb),nvar,ratios,not_covered, &
                             mask,DIMS(mask),crse,DIMS(crse),derives)

  integer  nvar, ratios(2), not_covered
  integer  lo(SDIM), hi(SDIM)
  integer  DIMDEC(bdry)
  integer  DIMDEC(mask)
  integer  DIMDEC(crse)
  integer  DIMDEC(cb)
  REAL_T   bdry(DIMV(bdry),nvar)
  REAL_T   derives(DIM2(cb),NUMDERIV)      
  integer  mask(DIMV(mask))
  REAL_T   crse(DIMV(crse),nvar)
  
  REAL_T   xx, ainterp
  integer  i, j, ic, jc, off, n
  integer  jclo, jchi, ratioy

  ratioy = ratios(2)

  jclo = ARG_L2(cb)
  jchi = ARG_H2(cb)
  ic   = ARG_L1(cb)-1
  i    = lo(1)-1

  do n = 1, nvar
     ! define interp coefs
     do jc = jclo, jchi
        j = ratioy*jc
        derives(jc,XDER)  = half * (crse(ic,jc+1,n) - crse(ic,jc-1,n))
        derives(jc,X2DER) = half * (crse(ic,jc+1,n) &
                                  - crse(ic,jc  ,n) * two &
                                  + crse(ic,jc-1,n))

        if (mask(i,j-1) .ne. not_covered) then
           derives(jc,XDER)  = crse(ic,jc+1,n) - crse(ic,jc  ,n)
           derives(jc,X2DER) = zero
        endif
        if (mask(i,j+ratioy) .ne. not_covered) then
           derives(jc,XDER)  = crse(ic,jc  ,n) - crse(ic,jc-1,n)
           derives(jc,X2DER) = zero
        endif
        if (mask(i,j-1)     .ne. not_covered .and. &
            mask(i,j+ratioy) .ne. not_covered) then
           derives(jc,XDER)  = zero
        endif
     enddo

     ! interpolate to fine grid
     do off = 0, ratioy - 1
        xx = (off - half * ratioy + half)/ratioy
        do jc = jclo, jchi
           j = ratioy*jc + off
           bdry(i,j,n) = crse(ic,jc,n) &
                + derives(jc,XDER) *xx &
                + derives(jc,X2DER)*xx**2
        enddo
     enddo
  enddo

  return
end subroutine FORT_BDINTERPXLO


! ---------------------------------------------------------------
! ::  FORT_BDINTERPXHI : Interpolation on Xhi Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  DIMS(bdry)  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  DIMS(cb)    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(2)   => refinement ratios
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  DIMS(mask)  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  DIMS(crse)  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

subroutine FORT_BDINTERPXHI (bdry,DIMS(bdry), &
                             lo,hi,DIMS(cb),nvar,ratios,not_covered, &
                             mask,DIMS(mask),crse,DIMS(crse),derives)

  integer  nvar, ratios(2), not_covered
  integer  lo(SDIM), hi(SDIM)
  integer  DIMDEC(bdry)
  integer  DIMDEC(mask)
  integer  DIMDEC(cb)
  integer  DIMDEC(crse)
  REAL_T   bdry(DIMV(bdry),nvar)
  REAL_T   derives(DIM2(cb),NUMDERIV)      
  integer  mask(DIMV(mask))
  REAL_T   crse(DIMV(crse),nvar)
  
  REAL_T   xx, ainterp
  integer  i, j, ic, jc, off, n
  integer  jclo, jchi, ratioy
  
  ratioy = ratios(2)

  jclo = ARG_L2(cb)
  jchi = ARG_H2(cb)
  ic   = ARG_H1(cb)+1
  i    = hi(1)+1
      
  do n = 1, nvar
     ! define interp coefs
     do jc = jclo, jchi
        j = ratioy*jc
        derives(jc,XDER)  = half * (crse(ic,jc+1,n) - crse(ic,jc-1,n))
        derives(jc,X2DER) = half * (crse(ic,jc+1,n) &
                                  - crse(ic,jc  ,n) * two &
                                  + crse(ic,jc-1,n))

        if (mask(i,j-1) .ne. not_covered) then
           derives(jc,XDER)  = crse(ic,jc+1,n) - crse(ic,jc  ,n)
           derives(jc,X2DER) = zero
        endif
        if (mask(i,j+ratioy) .ne. not_covered) then
           derives(jc,XDER)  = crse(ic,jc  ,n) - crse(ic,jc-1,n)
           derives(jc,X2DER) = zero
        endif
        if (mask(i,j-1) .ne. not_covered .and. &
            mask(i,j+ratioy) .ne. not_covered) then
           derives(jc,XDER)  = zero
        endif
     enddo

     ! interpolate to fine grid
     do off = 0, ratioy - 1
        xx = (off - half * ratioy + half)/ratioy
        do jc = jclo, jchi
           j = ratioy*jc + off
           bdry(i,j,n) = crse(ic,jc,n) &
                + derives(jc,XDER) *xx &
                + derives(jc,X2DER)*xx**2
        enddo
     enddo
  enddo
  
  return
end subroutine FORT_BDINTERPXHI


! ---------------------------------------------------------------
! ::  FORT_BDINTERPYLO : Interpolation on Ylo Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  DIMS(bdry)  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  DIMS(cb)    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(2)   => refinement ratios
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  DIMS(mask)  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  DIMS(crse)  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

subroutine FORT_BDINTERPYLO (bdry,DIMS(bdry), &
                             lo,hi,DIMS(cb),nvar,ratios,not_covered, &
                             mask,DIMS(mask),crse,DIMS(crse),derives)

  integer  nvar, ratios(2), not_covered
  integer  lo(SDIM), hi(SDIM)
  integer  DIMDEC(bdry)
  integer  DIMDEC(mask)
  integer  DIMDEC(cb)
  integer  DIMDEC(crse)
  REAL_T   bdry(DIMV(bdry),nvar)
  REAL_T   derives(DIM1(cb),NUMDERIV)
  integer  mask(DIMV(mask))
  REAL_T   crse(DIMV(crse),nvar)
  REAL_T   xx, ainterp
  integer  i, j, ic, jc, off, n
  integer  iclo, ichi, ratiox

  ratiox = ratios(1)

  iclo = ARG_L1(cb)
  ichi = ARG_H1(cb)
  jc   = ARG_L2(cb)-1
  j    = lo(2)-1
  
  do n = 1, nvar
     ! define interp coefs
     do ic = iclo, ichi
        i = ratiox*ic
        derives(ic,XDER)  = half * (crse(ic+1,jc,n) - crse(ic-1,jc,n))
        derives(ic,X2DER) = half * (crse(ic+1,jc,n) &
                                  - crse(ic  ,jc,n) * two &
                                  + crse(ic-1,jc,n))
        if (mask(i-1,j) .ne. not_covered) then
           derives(ic,XDER)  = crse(ic+1,jc,n) - crse(ic  ,jc,n)
           derives(ic,X2DER) = zero
        endif
        if (mask(i+ratiox,j) .ne. not_covered) then
           derives(ic,XDER)  = crse(ic  ,jc,n) - crse(ic-1,jc,n)
           derives(ic,X2DER) = zero
        endif
        if (mask(i-1,j) .ne. not_covered .and. &
            mask(i+ratiox,j) .ne. not_covered) then
           derives(ic,XDER)  = zero
        endif
     enddo

     ! interpolate to fine grid
     do off = 0, ratiox - 1
        xx = (off - half * ratiox + half)/ratiox
        do ic = iclo, ichi
           i = ratiox*ic + off
           bdry(i,j,n) = crse(ic,jc,n) &
                + derives(ic,XDER) *xx &
                + derives(ic,X2DER)*xx**2
        enddo
     enddo
  enddo
  
  return
end subroutine FORT_BDINTERPYLO


! ---------------------------------------------------------------
! ::  FORT_BDINTERPYHI : Interpolation on Yhi Face
! ::       Quadratic Interpolation from crse data
! ::       in directions transverse to face of grid
! ::
! ::  Inputs/Outputs:
! ::  bdry       <=  fine grid bndry data strip
! ::  DIMS(bdry)  => index limits of bdry
! ::  lo,hi       => index limits of grd interior
! ::  DIMS(cb)    => index limits of coarsened grid interior
! ::  nvar        => number of variables to interpolate
! ::  ratios(2)   => refinement ratios
! ::  not_covered => mask is set to this value if cell is not
! ::                 covered by another fine grid and not outside the domain.
! ::  mask        => fine grid mask bndry strip
! ::  DIMS(mask)  => index limits of mask array
! ::  crse        => crse grid bndry data strip
! ::  DIMS(crse)  => index limits of crse array
! ::  derives     => crse grid tmp array
! ---------------------------------------------------------------

subroutine FORT_BDINTERPYHI (bdry,DIMS(bdry), &
                             lo,hi,DIMS(cb),nvar,ratios,not_covered, &
                             mask,DIMS(mask),crse,DIMS(crse),derives)

  integer  nvar, ratios(2), not_covered
  integer  lo(SDIM), hi(SDIM)
  integer  DIMDEC(bdry)
  integer  DIMDEC(mask)
  integer  DIMDEC(cb)
  integer  DIMDEC(crse)
  REAL_T   bdry(DIMV(bdry),nvar)
  REAL_T   derives(DIM1(cb),NUMDERIV)
  integer  mask(DIMV(mask))
  REAL_T   crse(DIMV(crse),nvar)
  REAL_T   xx, ainterp
  integer  i, j, ic, jc, off, n
  integer  iclo, ichi, ratiox

  ratiox = ratios(1)

  iclo = ARG_L1(cb)
  ichi = ARG_H1(cb)
  jc   = ARG_H2(cb)+1
  j    = hi(2)+1

  do n = 1, nvar
     ! define interp coefs
     do ic = iclo, ichi
        i = ratiox*ic
        derives(ic,XDER)  = half * (crse(ic+1,jc,n) - crse(ic-1,jc,n))
        derives(ic,X2DER) = half * (crse(ic+1,jc,n) &
                                  - crse(ic  ,jc,n) * two &
                                  + crse(ic-1,jc,n))
        if (mask(i-1,j) .ne. not_covered) then
           derives(ic,XDER)  = crse(ic+1,jc,n) - crse(ic  ,jc,n)
           derives(ic,X2DER) = zero
        endif
        if (mask(i+ratiox,j) .ne. not_covered) then
           derives(ic,XDER)  = crse(ic  ,jc,n) - crse(ic-1,jc,n)
           derives(ic,X2DER) = zero
        endif
        if (mask(i-1,j) .ne. not_covered .and. &
            mask(i+ratiox,j) .ne. not_covered) then
           derives(ic,XDER)  = zero
        endif
     enddo

     ! interpolate to fine grid
     do off = 0, ratiox - 1
        xx = (off - half * ratiox + half)/ratiox
        do ic = iclo, ichi
           i = ratiox*ic + off
           bdry(i,j,n) = crse(ic,jc,n) &
                + derives(ic,XDER) *xx &
                + derives(ic,X2DER)*xx**2
        enddo
     enddo
  enddo
  
  return
end subroutine FORT_BDINTERPYHI

