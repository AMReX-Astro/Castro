module tagging_module

  implicit none

  double precision, save ::    denerr,   dengrad
  double precision, save ::    enterr,   entgrad
  double precision, save ::    velerr,   velgrad
  double precision, save ::   temperr,  tempgrad
  double precision, save ::  presserr, pressgrad
  double precision, save ::    raderr,   radgrad
  integer         , save ::  max_denerr_lev,   max_dengrad_lev
  integer         , save ::  max_enterr_lev,   max_entgrad_lev
  integer         , save ::  max_velerr_lev,   max_velgrad_lev
  integer         , save ::  max_temperr_lev,  max_tempgrad_lev
  integer         , save ::  max_presserr_lev, max_pressgrad_lev
  integer         , save ::  max_raderr_lev,   max_radgrad_lev

  public

contains

  ! For a given error estimator, check whether we fail to meet any of
  ! the criteria for tagging on this level. If so, set flag to zero.

  subroutine check_tagging_criteria(int_name, int_name_length, level, time, flag) bind(C)

    use meth_params_module, only: dxnuc

    implicit none

    integer, intent(in)          :: int_name(int_name_length)
    integer, intent(in)          :: int_name_length
    integer, intent(in)          :: level
    double precision, intent(in) :: time
    integer, intent(inout)       :: flag

    character(len=:), allocatable :: name
    integer                       :: i

    ! First, convert the integer name to a Fortran character type.

    allocate(character(len=int_name_length) :: name)

    do i = 1, int_name_length
       name(i:i) = char(int_name(i))
    enddo

    ! Now check the criteria for this type.

    if (name == "density") then

       if (level .ge. max_denerr_lev .and. level .ge. max_dengrad_lev) then

          flag = 0

       endif

    else if (name == "x_velocity" .or. name == "y_velocity" .or. name == "z_velocity") then

       if (level .ge. max_velerr_lev .and. level .ge. max_velgrad_lev) then

          flag = 0

       endif

    else if (name == "entropy") then

       if (level .ge. max_enterr_lev .and. level .ge. max_entgrad_lev) then

          flag = 0

       endif

    else if (name == "Temp") then

       if (level .ge. max_temperr_lev .and. level .ge. max_tempgrad_lev) then

          flag = 0

       endif

    else if (name == "pressure") then

       if (level .ge. max_presserr_lev .and. level .ge. max_pressgrad_lev) then

          flag = 0

       endif

    else if (name == "t_sound_t_enuc") then

       if (dxnuc > 1.d199) then

          flag = 0

       endif

    else if (name == "rad") then

       if (level .ge. max_raderr_lev .and. level .ge. max_radgrad_lev) then

          flag = 0

       endif

    endif

  end subroutine check_tagging_criteria



  ! All tagging subroutines in this file must be threadsafe because
  ! they are called inside OpenMP parallel regions.

! ::: -----------------------------------------------------------
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => variable to test on
! ::: np        => number of components in temp array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------


  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the Laplacian.
  ! ::: -----------------------------------------------------------

  subroutine ca_laplac_error(tag,taglo,taghi, &
                             set,clear, &
                             var,varlo,varhi, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level) bind(C)

    use prob_params_module, only: dg, dim

    implicit none

    integer          :: set, clear, nd, level
    integer          :: taglo(3), taghi(3)
    integer          :: varlo(3), varhi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: var(varlo(1):varhi(1),varlo(2):varhi(2),varlo(3):varhi(3))
    double precision :: delta(3), xlo(3), problo(3), time

    integer          :: i, j, k
    double precision ::  delu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
    double precision :: delua(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
    double precision :: delu2(9), delu3(9), delu4(9)
    double precision :: num, denom, error

    ! This value is  taken from FLASH
    double precision, parameter :: ctore=0.8
    double precision, parameter :: epsil=0.02

    ! adapted from ref_marking.f90 in FLASH2.5

    delu = 0.0
    delua = 0.0

    ! d/dx
    do k=lo(3)-1*dg(3),hi(3)+1*dg(3)
       do j=lo(2)-1*dg(2),hi(2)+1*dg(2)
          do i=lo(1)-1*dg(1),hi(1)+1*dg(1)
             delu(i,j,k,1)  =     var(i+1*dg(1),j,k)  -     var(i-1*dg(1),j,k)
             delua(i,j,k,1) = abs(var(i+1*dg(1),j,k)) + abs(var(i-1*dg(1),j,k))
          end do
       end do
    end do

    ! d/dy
    if (dim .ge. 2) then
       do k=lo(3)-1*dg(3),hi(3)+1*dg(3)
          do j=lo(2)-1*dg(2),hi(2)+1*dg(2)
             do i=lo(1)-1*dg(1),hi(1)+1*dg(1)
                delu(i,j,k,2)  =     var(i,j+1*dg(2),k)  -     var(i,j-1*dg(2),k) 
                delua(i,j,k,2) = abs(var(i,j+1*dg(2),k)) + abs(var(i,j-1*dg(2),k))
             end do
          end do
       end do
    endif

    ! d/dz
    if (dim .eq. 3) then
       do k=lo(3)-1*dg(3),hi(3)+1*dg(3)
          do j=lo(2)-1*dg(2),hi(2)+1*dg(2)
             do i=lo(1)-1*dg(1),hi(1)+1*dg(1)
                delu(i,j,k,3)  =     var(i,j,k+1*dg(3))  -     var(i,j,k-1*dg(3))
                delua(i,j,k,3) = abs(var(i,j,k+1*dg(3))) + abs(var(i,j,k-1*dg(3)))
             end do
          end do
       end do
    endif

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             ! d/dxdx
             delu2(1) =     delu(i+1,j,k,1)  -     delu(i-1,j,k,1)
             delu3(1) = abs(delu(i+1,j,k,1)) + abs(delu(i-1,j,k,1))
             delu4(1) =    delua(i+1,j,k,1)  +    delua(i-1,j,k,1)

             ! d/dydx
             delu2(2) =     delu(i,j+1,k,1)  -     delu(i,j-1,k,1)
             delu3(2) = abs(delu(i,j+1,k,1)) + abs(delu(i,j-1,k,1))
             delu4(2) =    delua(i,j+1,k,1)  +    delua(i,j-1,k,1)

             ! d/dxdy
             delu2(3) =     delu(i+1,j,k,2)  -     delu(i-1,j,k,2)
             delu3(3) = abs(delu(i+1,j,k,2)) + abs(delu(i-1,j,k,2))
             delu4(3) =    delua(i+1,j,k,2)  +    delua(i-1,j,k,2)

             ! d/dydy
             delu2(4) =     delu(i,j+1,k,2)  -     delu(i,j-1,k,2)
             delu3(4) = abs(delu(i,j+1,k,2)) + abs(delu(i,j-1,k,2))
             delu4(4) =    delua(i,j+1,k,2)  +    delua(i,j-1,k,2)

             ! d/dzdx
             delu2(5) =     delu(i,j,k+1,1)  -     delu(i,j,k-1,1)
             delu3(5) = abs(delu(i,j,k+1,1)) + abs(delu(i,j,k-1,1))
             delu4(5) =    delua(i,j,k+1,1)  +    delua(i,j,k-1,1)

             ! d/dzdy
             delu2(6) =     delu(i,j,k+1,2)  -     delu(i,j,k-1,2)
             delu3(6) = abs(delu(i,j,k+1,2)) + abs(delu(i,j,k-1,2))
             delu4(6) =    delua(i,j,k+1,2)  +    delua(i,j,k-1,2)

             ! d/dxdz
             delu2(7) =     delu(i+1,j,k,3)  -     delu(i-1,j,k,3)
             delu3(7) = abs(delu(i+1,j,k,3)) + abs(delu(i-1,j,k,3))
             delu4(7) =    delua(i+1,j,k,3)  +    delua(i-1,j,k,3)

             ! d/dydz
             delu2(8) =     delu(i,j+1,k,3)  -     delu(i,j-1,k,3)
             delu3(8) = abs(delu(i,j+1,k,3)) + abs(delu(i,j-1,k,3))
             delu4(8) =    delua(i,j+1,k,3)  +    delua(i,j-1,k,3)

             ! d/dzdz
             delu2(9) =     delu(i,j,k+1,3)  -     delu(i,j,k-1,3)
             delu3(9) = abs(delu(i,j,k+1,3)) + abs(delu(i,j,k-1,3))
             delu4(9) =    delua(i,j,k+1,3)  +    delua(i,j,k-1,3)

             ! compute the error
             num   = sum(delu2**2)

             denom = sum((delu3 + (epsil*delu4+1.d-99))**2)

             error = sqrt(num/denom)

             if (error .gt. ctore) tag(i,j,k) = set

          end do
       end do
    end do

  end subroutine ca_laplac_error

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the density
  ! ::: -----------------------------------------------------------

  subroutine ca_denerror(tag,taglo,taghi, &
                         set,clear, &
                         den,denlo,denhi, &
                         lo,hi,nd,domlo,domhi, &
                         delta,xlo,problo,time,level) bind(C)

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, nd, level
    integer          :: taglo(3), taghi(3)
    integer          :: denlo(3), denhi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: den(denlo(1):denhi(1),denlo(2):denhi(2),denlo(3):denhi(3),nd)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high density
    if (level .lt. max_denerr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (den(i,j,k,1) .ge. denerr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high density gradient
    if (level .lt. max_dengrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(den(i+1*dg(1),j,k,1) - den(i,j,k,1))
                ay = ABS(den(i,j+1*dg(2),k,1) - den(i,j,k,1))
                az = ABS(den(i,j,k+1*dg(3),1) - den(i,j,k,1))
                ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. dengrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_denerror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the temperature
  ! ::: -----------------------------------------------------------

  subroutine ca_temperror(tag,taglo,taghi, &
                          set,clear, &
                          temp,templo,temphi, &
                          lo,hi,np,domlo,domhi, &
                          delta,xlo,problo,time,level) bind(C)

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, np, level
    integer          :: taglo(3), taghi(3)
    integer          :: templo(3), temphi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: temp(templo(1):temphi(1),templo(2):temphi(2),templo(3):temphi(3),np)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high temperature
    if (level .lt. max_temperr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (temp(i,j,k,1) .ge. temperr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high temperature gradient
    if (level .lt. max_tempgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(temp(i+1*dg(1),j,k,1) - temp(i,j,k,1))
                ay = ABS(temp(i,j+1*dg(2),k,1) - temp(i,j,k,1))
                az = ABS(temp(i,j,k+1*dg(3),1) - temp(i,j,k,1))
                ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(temp(i,j,k,1) - temp(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(temp(i,j,k,1) - temp(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. tempgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_temperror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the pressure
  ! ::: -----------------------------------------------------------

  subroutine ca_presserror(tag,taglo,taghi, &
                           set,clear, &
                           press,presslo,presshi, &
                           lo,hi,np,domlo,domhi, &
                           delta,xlo,problo,time,level) bind(C)

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, np, level
    integer          :: taglo(3), taghi(3)
    integer          :: presslo(3), presshi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: press(presslo(1):presshi(1),presslo(2):presshi(2),presslo(3):presshi(3),np)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high pressure
    if (level .lt. max_presserr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (press(i,j,k,1) .ge. presserr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high pressure gradient
    if (level .lt. max_pressgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(press(i+1*dg(1),j,k,1) - press(i,j,k,1))
                ay = ABS(press(i,j+1*dg(2),k,1) - press(i,j,k,1))
                az = ABS(press(i,j,k+1*dg(3),1) - press(i,j,k,1))
                ax = MAX(ax,ABS(press(i,j,k,1) - press(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(press(i,j,k,1) - press(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(press(i,j,k,1) - press(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. pressgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_presserror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the velocity
  ! ::: -----------------------------------------------------------

  subroutine ca_velerror(tag,taglo,taghi, &
                         set,clear, &
                         vel,vello,velhi, &
                         lo,hi,nv,domlo,domhi, &
                         delta,xlo,problo,time,level) bind(C)

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, nv, level
    integer          :: taglo(3), taghi(3)
    integer          :: vello(3), velhi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3),nv)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high velocity
    if (level .lt. max_velerr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (vel(i,j,k,1) .ge. velerr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high velocity gradient
    if (level .lt. max_velgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(vel(i+1*dg(1),j,k,1) - vel(i,j,k,1))
                ay = ABS(vel(i,j+1*dg(2),k,1) - vel(i,j,k,1))
                az = ABS(vel(i,j,k+1*dg(3),1) - vel(i,j,k,1))
                ax = MAX(ax,ABS(vel(i,j,k,1) - vel(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(vel(i,j,k,1) - vel(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(vel(i,j,k,1) - vel(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. velgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_velerror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the radiation
  ! ::: -----------------------------------------------------------

  subroutine ca_raderror(tag,taglo,taghi, &
                         set,clear, &
                         rad,radlo,radhi, &
                         lo,hi,nr,domlo,domhi, &
                         delta,xlo,problo,time,level) bind(C)

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, nr, level
    integer          :: taglo(3), taghi(3)
    integer          :: radlo(3), radhi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: rad(radlo(1):radhi(1),radlo(2):radhi(2),radlo(3):radhi(3),nr)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high radiation
    if (level .lt. max_raderr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (rad(i,j,k,1) .ge. raderr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high radiation gradient
    if (level .lt. max_radgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(rad(i+1*dg(1),j,k,1) - rad(i,j,k,1))
                ay = ABS(rad(i,j+1*dg(2),k,1) - rad(i,j,k,1))
                az = ABS(rad(i,j,k+1*dg(3),1) - rad(i,j,k,1))
                ax = MAX(ax,ABS(rad(i,j,k,1) - rad(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(rad(i,j,k,1) - rad(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(rad(i,j,k,1) - rad(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. radgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_raderror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag high error cells based on the entropy
  ! ::: -----------------------------------------------------------

  subroutine ca_enterror(tag,taglo,taghi, &
                         set,clear, &
                         ent,entlo,enthi, &
                         lo,hi,nr,domlo,domhi, &
                         delta,xlo,problo,time,level) bind(C)

    use prob_params_module, only: dg

    implicit none

    integer          :: set, clear, nr, level
    integer          :: taglo(3), taghi(3)
    integer          :: entlo(3), enthi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: ent(entlo(1):enthi(1),entlo(2):enthi(2),entlo(3):enthi(3),nr)
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: ax, ay, az
    integer          :: i, j, k

    !     Tag on regions of high radiation
    if (level .lt. max_enterr_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (ent(i,j,k,1) .ge. enterr) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

    !     Tag on regions of high radiation gradient
    if (level .lt. max_entgrad_lev) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ax = ABS(ent(i+1*dg(1),j,k,1) - ent(i,j,k,1))
                ay = ABS(ent(i,j+1*dg(2),k,1) - ent(i,j,k,1))
                az = ABS(ent(i,j,k+1*dg(3),1) - ent(i,j,k,1))
                ax = MAX(ax,ABS(ent(i,j,k,1) - ent(i-1*dg(1),j,k,1)))
                ay = MAX(ay,ABS(ent(i,j,k,1) - ent(i,j-1*dg(2),k,1)))
                az = MAX(az,ABS(ent(i,j,k,1) - ent(i,j,k-1*dg(3),1)))
                if ( MAX(ax,ay,az) .ge. entgrad) then
                   tag(i,j,k) = set
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine ca_enterror

  ! ::: -----------------------------------------------------------
  ! ::: This routine will tag cells based on the sound crossing time
  ! ::: relative to the nuclear energy injection timescale.
  ! ::: At present we tag for maximal refinement since this
  ! ::: criterion is necessary for numerical burning stability.
  ! ::: -----------------------------------------------------------

  subroutine ca_nucerror(tag,taglo,taghi, &
                         set,clear, &
                         t,tlo,thi, &
                         lo,hi,nr,domlo,domhi, &
                         delta,xlo,problo,time,level) bind(C)

    use meth_params_module, only: dxnuc

    implicit none

    integer          :: set, clear, nr, level
    integer          :: taglo(3), taghi(3)
    integer          :: tlo(3), thi(3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    double precision :: t(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),nr) ! t_sound / t_e
    double precision :: delta(3), xlo(3), problo(3), time

    double precision :: t_sound, t_enuc
    integer          :: i, j, k

    ! Disable if we're not utilizing this tagging

    if (dxnuc > 1.d199) return

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (t(i,j,k,1) > dxnuc) then

                tag(i,j,k) = set

             endif

          enddo
       enddo
    enddo

  end subroutine ca_nucerror

end module tagging_module
