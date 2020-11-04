module tagging_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding, only : c_int8_t


  implicit none

  real(rt), allocatable ::    denerr,   dengrad, dengrad_rel
  real(rt), allocatable ::    enterr,   entgrad, entgrad_rel
  real(rt), allocatable ::    velerr,   velgrad, velgrad_rel
  real(rt), allocatable ::   temperr,  tempgrad, tempgrad_rel
  real(rt), allocatable ::  presserr, pressgrad, pressgrad_rel
  real(rt), allocatable ::    raderr,   radgrad, radgrad_rel
  real(rt), allocatable ::   enucerr

  integer, allocatable ::  max_denerr_lev,   max_dengrad_lev, max_dengrad_rel_lev
  integer, allocatable ::  max_velerr_lev,   max_velgrad_lev, max_velgrad_rel_lev
  integer, allocatable ::  max_temperr_lev,  max_tempgrad_lev, max_tempgrad_rel_lev
  integer, allocatable ::  max_presserr_lev, max_pressgrad_lev, max_pressgrad_rel_lev
  integer, allocatable ::  max_raderr_lev,   max_radgrad_lev, max_radgrad_rel_lev
  integer, allocatable ::  max_enucerr_lev

  ! limit the zone size based on how much the burning can change the
  ! internal energy of a zone. The zone size on the finest level must
  ! be smaller than dxnuc * c_s * (e/ \dot{e}) where c_s is the sound
  ! speed.  This ensures that the sound-crossing time is smaller than
  ! the nuclear energy injection timescale.
  real(rt), allocatable :: dxnuc_min

  ! Disable limiting based on dxnuc above this threshold. This allows
  !  zones that have already ignited or are about to ignite to be
  !  de-refined.
  real(rt), allocatable :: dxnuc_max

  ! Disable limiting based on dxnuc above this AMR level.
  integer, allocatable :: max_dxnuc_lev

  public

#if (defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA))
  attributes(managed) ::    denerr,   dengrad, dengrad_rel
  attributes(managed) ::    velerr,   velgrad, velgrad_rel
  attributes(managed) ::   temperr,  tempgrad, tempgrad_rel
  attributes(managed) ::  presserr, pressgrad, pressgrad_rel
  attributes(managed) ::    raderr,   radgrad, radgrad_rel
  attributes(managed) ::   enucerr

  attributes(managed) ::  max_denerr_lev,   max_dengrad_lev, max_dengrad_rel_lev
  attributes(managed) ::  max_velerr_lev,   max_velgrad_lev, max_velgrad_rel_lev
  attributes(managed) ::  max_temperr_lev,  max_tempgrad_lev, max_tempgrad_rel_lev
  attributes(managed) ::  max_presserr_lev, max_pressgrad_lev, max_pressgrad_rel_lev
  attributes(managed) ::  max_raderr_lev,   max_radgrad_lev, max_radgrad_rel_lev
  attributes(managed) ::  max_enucerr_lev

  attributes(managed) :: dxnuc_min
  attributes(managed) :: dxnuc_max
  attributes(managed) :: max_dxnuc_lev
#endif

contains

  AMREX_CUDA_FORT_DEVICE subroutine ca_denerror(i, j, k, &
                                                tag, taglo, taghi, &
                                                den, denlo, denhi, nd, &
                                                delta, problo, &
                                                set, clear, time, level) &
                                                bind(C, name="ca_denerror")
     !
     ! This routine will tag high error cells based on the density
     !

    use prob_params_module, only: dg

    implicit none

    integer,    intent(in   ), value :: i, j, k
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: denlo(3), denhi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: den(denlo(1):denhi(1),denlo(2):denhi(2),denlo(3):denhi(3),nd)
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: nd, level
    real(rt),   intent(in   ), value :: time

    real(rt) :: ax, ay, az

    !     Tag on regions of high density
    if (level .lt. max_denerr_lev) then
       if (den(i,j,k,1) .ge. denerr) then
          tag(i,j,k) = set
       endif
    endif

    !     Tag on regions of high density gradient
    if (level .lt. max_dengrad_lev .or. level .lt. max_dengrad_rel_lev) then
       ax = ABS(den(i+1*dg(1),j,k,1) - den(i,j,k,1))
       ay = ABS(den(i,j+1*dg(2),k,1) - den(i,j,k,1))
       az = ABS(den(i,j,k+1*dg(3),1) - den(i,j,k,1))
       ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1*dg(1),j,k,1)))
       ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1*dg(2),k,1)))
       az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1*dg(3),1)))
       if (MAX(ax,ay,az) .ge. dengrad .or. MAX(ax,ay,az) .ge. ABS(dengrad_rel * den(i,j,k,1))) then
          tag(i,j,k) = set
       endif
    endif

  end subroutine ca_denerror



  AMREX_CUDA_FORT_DEVICE subroutine ca_temperror(i, j, k, &
                                                 tag, taglo, taghi, &
                                                 temp, templo, temphi, np, &
                                                 delta, problo, &
                                                 set, clear, time, level) &
                                                 bind(C, name="ca_temperror")
  !
  ! This routine will tag high error cells based on the temperature
  !

    use prob_params_module, only: dg

    implicit none

    integer,    intent(in   ), value :: i, j, k
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: templo(3), temphi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: temp(templo(1):temphi(1),templo(2):temphi(2),templo(3):temphi(3),np)
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: np, level
    real(rt),   intent(in   ), value :: time

    real(rt) :: ax, ay, az

    !     Tag on regions of high temperature
    if (level .lt. max_temperr_lev) then
       if (temp(i,j,k,1) .ge. temperr) then
          tag(i,j,k) = set
       endif
    endif

    !     Tag on regions of high temperature gradient
    if (level .lt. max_tempgrad_lev .or. level .lt. max_tempgrad_rel_lev) then
       ax = ABS(temp(i+1*dg(1),j,k,1) - temp(i,j,k,1))
       ay = ABS(temp(i,j+1*dg(2),k,1) - temp(i,j,k,1))
       az = ABS(temp(i,j,k+1*dg(3),1) - temp(i,j,k,1))
       ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1*dg(1),j,k,1)))
       ay = MAX(ay,ABS(temp(i,j,k,1) - temp(i,j-1*dg(2),k,1)))
       az = MAX(az,ABS(temp(i,j,k,1) - temp(i,j,k-1*dg(3),1)))
       if (MAX(ax,ay,az) .ge. tempgrad .or. MAX(ax,ay,az) .ge. ABS(tempgrad_rel * temp(i,j,k,1))) then
          tag(i,j,k) = set
       endif
    endif

  end subroutine ca_temperror



  AMREX_CUDA_FORT_DEVICE subroutine ca_presserror(i, j, k, &
                                                  tag, taglo, taghi, &
                                                  press, presslo, presshi, np, &
                                                  delta, problo, &
                                                  set, clear, time, level) &
                                                  bind(C, name="ca_presserror")
   !
   ! This routine will tag high error cells based on the pressure
   !

    use prob_params_module, only: dg
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,    intent(in   ), value :: i, j, k
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: presslo(3), presshi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: press(presslo(1):presshi(1),presslo(2):presshi(2),presslo(3):presshi(3),np)
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: np, level
    real(rt),   intent(in   ), value :: time

    real(rt) :: ax, ay, az

    !     Tag on regions of high pressure
    if (level .lt. max_presserr_lev) then
       if (press(i,j,k,1) .ge. presserr) then
          tag(i,j,k) = set
       endif
    endif

    !     Tag on regions of high pressure gradient
    if (level .lt. max_pressgrad_lev .or. level .lt. max_pressgrad_rel_lev) then
       ax = ABS(press(i+1*dg(1),j,k,1) - press(i,j,k,1))
       ay = ABS(press(i,j+1*dg(2),k,1) - press(i,j,k,1))
       az = ABS(press(i,j,k+1*dg(3),1) - press(i,j,k,1))
       ax = MAX(ax,ABS(press(i,j,k,1) - press(i-1*dg(1),j,k,1)))
       ay = MAX(ay,ABS(press(i,j,k,1) - press(i,j-1*dg(2),k,1)))
       az = MAX(az,ABS(press(i,j,k,1) - press(i,j,k-1*dg(3),1)))
       if (MAX(ax,ay,az) .ge. pressgrad .or. MAX(ax,ay,az) .ge. ABS(pressgrad_rel * press(i,j,k,1))) then
          tag(i,j,k) = set
       endif
    endif

  end subroutine ca_presserror



  AMREX_CUDA_FORT_DEVICE subroutine ca_velerror(i, j, k, &
                                                tag, taglo, taghi, &
                                                vel, vello, velhi, nv, &
                                                delta, problo, &
                                                set, clear, time, level) &
                                                bind(C, name="ca_velerror")
    !
    ! This routine will tag high error cells based on the velocity
    !

    use prob_params_module, only: dg

    implicit none

    integer,    intent(in   ), value :: i, j, k
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: vello(3), velhi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3),nv)
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: nv, level
    real(rt),   intent(in   ), value :: time

    real(rt) :: ax, ay, az

    !     Tag on regions of high velocity
    if (level .lt. max_velerr_lev) then
       if (abs(vel(i,j,k,1)) .ge. velerr) then
          tag(i,j,k) = set
       endif
    endif

    !     Tag on regions of high velocity gradient
    if (level .lt. max_velgrad_lev .or. level .lt. max_velgrad_rel_lev) then
       ax = ABS(vel(i+1*dg(1),j,k,1) - vel(i,j,k,1))
       ay = ABS(vel(i,j+1*dg(2),k,1) - vel(i,j,k,1))
       az = ABS(vel(i,j,k+1*dg(3),1) - vel(i,j,k,1))
       ax = MAX(ax,ABS(vel(i,j,k,1) - vel(i-1*dg(1),j,k,1)))
       ay = MAX(ay,ABS(vel(i,j,k,1) - vel(i,j-1*dg(2),k,1)))
       az = MAX(az,ABS(vel(i,j,k,1) - vel(i,j,k-1*dg(3),1)))
       if (MAX(ax,ay,az) .ge. velgrad .or. MAX(ax,ay,az) .ge. ABS(velgrad_rel * vel(i,j,k,1))) then
          tag(i,j,k) = set
       endif
    endif

  end subroutine ca_velerror



  AMREX_CUDA_FORT_DEVICE subroutine ca_raderror(i, j, k, &
                                                tag, taglo, taghi, &
                                                rad, radlo, radhi, nr, &
                                                delta, problo, &
                                                set, clear, time, level) &
                                                bind(C, name="ca_raderror")
    !
    ! This routine will tag high error cells based on the radiation
    !

    use prob_params_module, only: dg

    implicit none

    integer,    intent(in   ), value :: i, j, k
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: radlo(3), radhi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: rad(radlo(1):radhi(1),radlo(2):radhi(2),radlo(3):radhi(3),nr)
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: nr, level
    real(rt),   intent(in   ), value :: time

    real(rt) :: ax, ay, az

    !     Tag on regions of high radiation
    if (level .lt. max_raderr_lev) then
       if (rad(i,j,k,1) .ge. raderr) then
          tag(i,j,k) = set
       endif
    endif

    !     Tag on regions of high radiation gradient
    if (level .lt. max_radgrad_lev .or. level .lt. max_radgrad_rel_lev) then
       ax = ABS(rad(i+1*dg(1),j,k,1) - rad(i,j,k,1))
       ay = ABS(rad(i,j+1*dg(2),k,1) - rad(i,j,k,1))
       az = ABS(rad(i,j,k+1*dg(3),1) - rad(i,j,k,1))
       ax = MAX(ax,ABS(rad(i,j,k,1) - rad(i-1*dg(1),j,k,1)))
       ay = MAX(ay,ABS(rad(i,j,k,1) - rad(i,j-1*dg(2),k,1)))
       az = MAX(az,ABS(rad(i,j,k,1) - rad(i,j,k-1*dg(3),1)))
       if (MAX(ax,ay,az) .ge. radgrad .or. MAX(ax,ay,az) .ge. ABS(radgrad_rel * rad(i,j,k,1))) then
          tag(i,j,k) = set
       endif
    endif

  end subroutine ca_raderror



#ifdef REACTIONS
  AMREX_CUDA_FORT_DEVICE subroutine ca_nucerror(i, j, k, &
                                                tag, taglo, taghi, &
                                                t, tlo, thi, nr, &
                                                delta, problo, &
                                                set, clear, time, level) &
                                                bind(C, name="ca_nucerror")
    !
    ! This routine will tag cells based on the sound crossing time
    ! relative to the nuclear energy injection timescale.
    !

    implicit none

    integer,    intent(in   ), value :: i, j, k
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: tlo(3), thi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: t(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),nr) ! t_sound / t_e
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: nr, level
    real(rt),   intent(in   ), value :: time

    ! Disable if we're not utilizing this tagging

    if (dxnuc_min > 1.e199_rt) return

    if (level .lt. max_dxnuc_lev) then
       if (t(i,j,k,1) > dxnuc_min .and. t(i,j,k,1) < dxnuc_max) then
          tag(i,j,k) = set
       endif
    end if

  end subroutine ca_nucerror



  AMREX_CUDA_FORT_DEVICE subroutine ca_enucerror(i, j, k, &
                                                 tag, taglo, taghi, &
                                                 enuc,enuclo,enuchi, nd, &
                                                 delta, problo, &
                                                 set, clear, time, level) &
                                                 bind(C, name="ca_enucerror")
    !
    ! This routine will tag high error cells based on the nuclear
    ! energy generation rate
    !

    implicit none

    integer,    intent(in   ), value :: i, j, k
    integer,    intent(in   ) :: taglo(3), taghi(3)
    integer,    intent(in   ) :: enuclo(3), enuchi(3)
    integer(kind=c_int8_t), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(rt),   intent(in   ) :: enuc(enuclo(1):enuchi(1),enuclo(2):enuchi(2),enuclo(3):enuchi(3),nd)
    real(rt),   intent(in   ) :: delta(3), problo(3)
    integer(kind=c_int8_t), intent(in   ), value :: set, clear
    integer,    intent(in   ), value :: nd, level
    real(rt),   intent(in   ), value :: time

    ! Tag on regions of high nuclear energy generation rate
    if (level .lt. max_enucerr_lev) then
       if (enuc(i,j,k,1) .ge. enucerr) then
          tag(i,j,k) = set
       endif
    endif

  end subroutine ca_enucerror
#endif



  ! Routines for retrieving the maximum tagging level.

  subroutine get_max_denerr_lev(lev) bind(c, name='get_max_denerr_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_denerr_lev

  end subroutine get_max_denerr_lev



  subroutine get_max_dengrad_lev(lev) bind(c, name='get_max_dengrad_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_dengrad_lev

  end subroutine get_max_dengrad_lev



  subroutine get_max_dengrad_rel_lev(lev) bind(c, name='get_max_dengrad_rel_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_dengrad_rel_lev

  end subroutine get_max_dengrad_rel_lev



  subroutine get_max_velerr_lev(lev) bind(c, name='get_max_velerr_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_velerr_lev

  end subroutine get_max_velerr_lev



  subroutine get_max_velgrad_lev(lev) bind(c, name='get_max_velgrad_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_velgrad_lev

  end subroutine get_max_velgrad_lev



  subroutine get_max_velgrad_rel_lev(lev) bind(c, name='get_max_velgrad_rel_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_velgrad_rel_lev

  end subroutine get_max_velgrad_rel_lev



  subroutine get_max_temperr_lev(lev) bind(c, name='get_max_temperr_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_temperr_lev

  end subroutine get_max_temperr_lev



  subroutine get_max_tempgrad_lev(lev) bind(c, name='get_max_tempgrad_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_tempgrad_lev

  end subroutine get_max_tempgrad_lev



  subroutine get_max_tempgrad_rel_lev(lev) bind(c, name='get_max_tempgrad_rel_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_tempgrad_rel_lev

  end subroutine get_max_tempgrad_rel_lev



  subroutine get_max_presserr_lev(lev) bind(c, name='get_max_presserr_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_presserr_lev

  end subroutine get_max_presserr_lev



  subroutine get_max_pressgrad_lev(lev) bind(c, name='get_max_pressgrad_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_pressgrad_lev

  end subroutine get_max_pressgrad_lev



  subroutine get_max_pressgrad_rel_lev(lev) bind(c, name='get_max_pressgrad_rel_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_pressgrad_rel_lev

  end subroutine get_max_pressgrad_rel_lev



  subroutine get_max_raderr_lev(lev) bind(c, name='get_max_raderr_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_raderr_lev

  end subroutine get_max_raderr_lev



  subroutine get_max_radgrad_lev(lev) bind(c, name='get_max_radgrad_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_radgrad_lev

  end subroutine get_max_radgrad_lev



  subroutine get_max_radgrad_rel_lev(lev) bind(c, name='get_max_radgrad_rel_lev')

    implicit none

    integer, intent(out) :: lev

    lev = max_radgrad_rel_lev

  end subroutine get_max_radgrad_rel_lev

end module tagging_module
