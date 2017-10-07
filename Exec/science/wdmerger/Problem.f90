! Problem-specific Fortran routines that are designed to interact with C++

subroutine problem_checkpoint(int_dir_name, len) bind(C, name="problem_checkpoint")

  ! called by the IO processor during checkpoint

  use bl_IO_module
  use probdata_module, only: com_P, com_S, vel_P, vel_S, mass_P, mass_S, t_ff_P, t_ff_S, &
                             T_global_max, rho_global_max, ts_te_global_max
  use prob_params_module, only: center
  use meth_params_module, only: rot_period
  use probdata_module, only: jobIsDone, signalJobIsNotDone, num_previous_ener_timesteps, total_ener_array, &
                             relaxation_is_done, num_zones_ignited, ignition_level

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i, un

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

  un = unit_new()
  open (unit=un, file=trim(dir)//"/COM", status="unknown")

100 format(1x, g30.20)
200 format(1x, g30.20, 1x, g30.20)
300 format(1x, g30.20, 1x, g30.20, 1x, g30.20)

  write (un,300) center(1), center(2), center(3)
  write (un,200) mass_P,   mass_S
  write (un,200) com_P(1), com_S(1)
  write (un,200) com_P(2), com_S(2)
  write (un,200) com_P(3), com_S(3)
  write (un,200) vel_P(1), vel_S(1)
  write (un,200) vel_P(2), vel_S(2)
  write (un,200) vel_P(3), vel_S(3)
  write (un,200) t_ff_P,   t_ff_S

  close (un)



  open (unit=un, file=trim(dir)//"/Rotation", status="unknown")

  write (un,100) rot_period

  close (un)



  open (unit=un, file=trim(dir)//"/Extrema", status="unknown")

  write (un,100) T_global_max
  write (un,100) rho_global_max
  write (un,100) ts_te_global_max

  close (un)



  open (unit=un, file=trim(dir)//"/Energy", status="unknown")

  do i = 1, num_previous_ener_timesteps

     write (un,100) total_ener_array(i)

  enddo

  close (un)



  open (unit=un, file=trim(dir)//"/Relaxation", status="unknown")

  write (un,100) relaxation_is_done

  close (un)



  open (unit=un, file=trim(dir)//"/Ignition", status="unknown")

  write (un,*) num_zones_ignited, ignition_level

  close (un)


  ! If the job is done, write a file in the checkpoint indicating this.

  if (jobIsDone) then

     open  (unit=un, file=trim(dir)//"/jobIsDone", status="unknown")
     close (un)

  else

     if (signalJobIsNotDone .and. .not. jobIsDone) then

        open  (unit=un, file=trim(dir)//"/jobIsNotDone", status="unknown")
        close (un)

     endif

  endif

end subroutine problem_checkpoint



subroutine problem_restart(int_dir_name, len) bind(C, name="problem_restart")

  ! called by ALL processors during restart 

  use bl_IO_module
  use probdata_module, only: com_P, com_S, vel_P, vel_S, mass_P, mass_S, t_ff_P, t_ff_S, problem, &
                             T_global_max, rho_global_max, ts_te_global_max, &
                             jobIsDone, num_previous_ener_timesteps, total_ener_array, &
                             problem, relaxation_is_done, num_zones_ignited, ignition_level
  use wdmerger_util_module, only: turn_off_relaxation
  use problem_io_module, only: ioproc
  use prob_params_module, only: center
  use meth_params_module, only: rot_period

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i, un, stat
  logical :: fileExists

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

  un = unit_new()
  open (unit=un, file=trim(dir)//"/COM", status="old")

100 format(1x, g30.20)
200 format(1x, g30.20, 1x, g30.20)
300 format(1x, g30.20, 1x, g30.20, 1x, g30.20)

  read (un,300) center(1), center(2), center(3)
  read (un,200) mass_P,   mass_S
  read (un,200) com_P(1), com_S(1)
  read (un,200) com_P(2), com_S(2)
  read (un,200) com_P(3), com_S(3)
  read (un,200) vel_P(1), vel_S(1)
  read (un,200) vel_P(2), vel_S(2)
  read (un,200) vel_P(3), vel_S(3)
  read (un,200) t_ff_P,   t_ff_S

  close (un)



  if (problem .eq. 3) then

     open (unit=un, file=trim(dir)//"/Rotation", status="old", IOSTAT = stat)

     if (stat .eq. 0) then

        read (un,100) rot_period

        if (ioproc) then
           print *, ""
           write(*,1001) rot_period
           print *, ""

1001       format ("  Based on the checkpoint, updating the rotational period to ", F7.3 " s.")
        endif

        close (un)

     else

        rot_period = -1.0d0

     endif

  endif



  open (unit=un, file=trim(dir)//"/Extrema", status="old", IOSTAT = stat)

  if (stat .eq. 0) then

     read (un,100) T_global_max
     read (un,100) rho_global_max
     read (un,100) ts_te_global_max

     close(un)

  else

     T_global_max = 0.0d0
     rho_global_max = 0.0d0
     ts_te_global_max = 0.0d0

  endif



  open (unit=un, file=trim(dir)//"/Energy", status="old", IOSTAT = stat)

  if (stat .eq. 0) then

     do i = 1, num_previous_ener_timesteps

        read (un,100) total_ener_array(i)

     enddo

     close (un)

  else

     total_ener_array(:) = -1.d200

  endif



  open (unit=un, file=trim(dir)//"/Relaxation", status="old", IOSTAT = stat)

  if (stat .eq. 0) then

     read (un,100) relaxation_is_done

     close (un)

  else

     if (problem == 3) then
        call bl_error("Error: no Relaxation file found in the checkpoint.")
     endif

  endif

  if (relaxation_is_done == 1) then

     call turn_off_relaxation(-1.0d0)

  endif



  open (unit=un, file=trim(dir)//"/Ignition", status="old", IOSTAT = stat)

  if (stat .eq. 0) then

     read (un,*) num_zones_ignited, ignition_level

     close (un)

  end if



  ! Pick up whether the job has been completed.

  inquire(file=trim(dir)//"/jobIsDone", EXIST = fileExists)

  if (fileExists) then

     jobIsDone = .true.

  else

     inquire(file=trim(dir)//"/jobIsNotDone", EXIST = fileExists)

     if (fileExists) then

        jobIsDone = .false.

     endif

  endif

end subroutine problem_restart



! Return the mass-weighted center of mass and velocity
! for the primary and secondary, for a given FAB.
! Since this will rely on a sum over processors,
! we should only add to the relevant variables
! in anticipation of a MPI reduction, and not
! overwrite them. Note that ultimately what we
! are doing here is to use an old guess at the
! effective potential of the primary and secondary
! to generate a new estimate.

subroutine wdcom(rho,  r_lo, r_hi, &
                 xmom, px_lo, px_hi, &
                 ymom, py_lo, py_hi, &
                 zmom, pz_lo, pz_hi, &
                 pmask, pm_lo, pm_hi, &
                 smask, sm_lo, sm_hi, &
                 vol,  vo_lo, vo_hi, &
                 lo, hi, dx, time, &
                 com_p_x, com_p_y, com_p_z, &
                 com_s_x, com_s_y, com_s_z, &
                 vel_p_x, vel_p_y, vel_p_z, &
                 vel_s_x, vel_s_y, vel_s_z, &
                 m_p, m_s) bind(C,name='wdcom')

  use bl_constants_module, only: HALF, ZERO, ONE
  use prob_params_module, only: problo, probhi, physbc_lo, physbc_hi, Symmetry
  use castro_util_module, only: position

  implicit none

  integer         , intent(in   ) :: r_lo(3), r_hi(3)
  integer         , intent(in   ) :: px_lo(3), px_hi(3)
  integer         , intent(in   ) :: py_lo(3), py_hi(3)
  integer         , intent(in   ) :: pz_lo(3), pz_hi(3)
  integer         , intent(in   ) :: pm_lo(3), pm_hi(3)
  integer         , intent(in   ) :: sm_lo(3), sm_hi(3)
  integer         , intent(in   ) :: vo_lo(3), vo_hi(3)

  double precision, intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
  double precision, intent(in   ) :: xmom(px_lo(1):px_hi(1),px_lo(2):px_hi(2),px_lo(3):px_hi(3))
  double precision, intent(in   ) :: ymom(py_lo(1):py_hi(1),py_lo(2):py_hi(2),py_lo(3):py_hi(3))
  double precision, intent(in   ) :: zmom(pz_lo(1):pz_hi(1),pz_lo(2):pz_hi(2),pz_lo(3):pz_hi(3))
  double precision, intent(in   ) :: pmask(pm_lo(1):pm_hi(1),pm_lo(2):pm_hi(2),pm_lo(3):pm_hi(3))
  double precision, intent(in   ) :: smask(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
  double precision, intent(in   ) :: vol(vo_lo(1):vo_hi(1),vo_lo(2):vo_hi(2),vo_lo(3):vo_hi(3))

  integer         , intent(in   ) :: lo(3), hi(3)
  double precision, intent(in   ) :: dx(3), time
  double precision, intent(inout) :: com_p_x, com_p_y, com_p_z
  double precision, intent(inout) :: com_s_x, com_s_y, com_s_z
  double precision, intent(inout) :: vel_p_x, vel_p_y, vel_p_z
  double precision, intent(inout) :: vel_s_x, vel_s_y, vel_s_z
  double precision, intent(inout) :: m_p, m_s

  integer          :: i, j, k
  double precision :: r(3), dm

  ! Add to the COM locations and velocities of the primary and secondary
  ! depending on which potential dominates, ignoring unbound material.
  ! Note that in this routine we actually are
  ! summing mass-weighted quantities for the COM and the velocity; 
  ! we will account for this at the end of the calculation in 
  ! post_timestep() by dividing by the mass.

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           ! Our convention is that the COM locations for the WDs are 
           ! absolute positions on the grid, not relative to the center.

           r = position(i,j,k)

           ! We account for symmetric boundaries in this sum as usual,
           ! by adding to the position the locations that would exist
           ! on the opposite side of the symmetric boundary.

           r = merge(r + (problo - r), r, physbc_lo(:) .eq. Symmetry)
           r = merge(r + (r - probhi), r, physbc_hi(:) .eq. Symmetry)

           dm = rho(i,j,k) * vol(i,j,k)

           if (pmask(i,j,k) > ZERO) then

              m_p = m_p + dm

              com_p_x = com_p_x + dm * r(1)
              com_p_y = com_p_y + dm * r(2)
              com_p_z = com_p_z + dm * r(3)

              vel_p_x = vel_p_x + xmom(i,j,k) * vol(i,j,k)
              vel_p_y = vel_p_y + ymom(i,j,k) * vol(i,j,k)
              vel_p_z = vel_p_z + zmom(i,j,k) * vol(i,j,k)

           else if (smask(i,j,k) > ZERO) then

              m_s = m_s + dm

              com_s_x = com_s_x + dm * r(1)
              com_s_y = com_s_y + dm * r(2)
              com_s_z = com_s_z + dm * r(3)

              vel_s_x = vel_s_x + xmom(i,j,k) * vol(i,j,k)
              vel_s_y = vel_s_y + ymom(i,j,k) * vol(i,j,k)
              vel_s_z = vel_s_z + zmom(i,j,k) * vol(i,j,k)

           endif

        enddo
     enddo
  enddo

end subroutine wdcom



! This function uses the known center of mass of the two white dwarfs,
! and given a density cutoff, computes the total volume of all zones
! whose density is greater or equal to that density cutoff.
! We also impose a distance requirement so that we only look
! at zones within the Roche lobe of the white dwarf.

subroutine ca_volumeindensityboundary(rho,r_lo,r_hi, &
                                      pmask,pm_lo,pm_hi, &
                                      smask,sm_lo,sm_hi, &
                                      vol,v_lo,v_hi, &
                                      lo,hi,dx, &
                                      volp,vols,rho_cutoff) &
                                      bind(C,name='ca_volumeindensityboundary')

  use bl_constants_module

  implicit none

  integer          :: r_lo(3), r_hi(3)
  integer          :: pm_lo(3), pm_hi(3)
  integer          :: sm_lo(3), sm_hi(3)
  integer          :: v_lo(3), v_hi(3)
  integer          :: lo(3), hi(3)
  double precision :: volp, vols, rho_cutoff, dx(3)
  double precision :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
  double precision :: pmask(pm_lo(1):pm_hi(1),pm_lo(2):pm_hi(2),pm_lo(3):pm_hi(3))
  double precision :: smask(sm_lo(1):sm_hi(1),sm_lo(2):sm_hi(2),sm_lo(3):sm_hi(3))
  double precision :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

  integer          :: i, j, k

  volp = ZERO
  vols = ZERO

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           if (rho(i,j,k) > rho_cutoff) then

              if (pmask(i,j,k) > ZERO) then

                 volp = volp + vol(i,j,k)

              else if (smask(i,j,k) > ZERO) then

                 vols = vols + vol(i,j,k)

              endif

           endif

        enddo
     enddo
  enddo

end subroutine ca_volumeindensityboundary



! Determine the critical Roche potential at the Lagrange point L1.
! We will use a tri-linear interpolation that gets a contribution
! from all the zone centers that bracket the Lagrange point.

subroutine get_critical_roche_potential(phiEff,p_lo,p_hi,lo,hi,L1,potential) &
                                        bind(C,name='get_critical_roche_potential')

  use bl_constants_module, only: ZERO, HALF, ONE
  use castro_util_module, only: position
  use prob_params_module, only: dim, dx_level
  use amrinfo_module, only: amr_level

  implicit none

  integer          :: lo(3), hi(3)
  integer          :: p_lo(3), p_hi(3)
  double precision :: phiEff(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
  double precision :: L1(3), potential

  double precision :: r(3), dx(3)
  integer          :: i, j, k

  dx = dx_level(:,amr_level)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           r = position(i,j,k) - L1

           ! Scale r by dx (in dimensions we're actually simulating).

           r(1:dim) = r(1:dim) / dx(1:dim)
           r(dim+1:3) = ZERO

           ! We want a contribution from this zone if it is
           ! less than one zone width away from the Lagrange point.

           if (sum(r**2) < ONE) then

              potential = potential + product(ONE - abs(r)) * phiEff(i,j,k)

           endif

        enddo
     enddo
  enddo

end subroutine get_critical_roche_potential



! Given state data in the rotating frame, transform it to the inertial frame.

subroutine transform_to_inertial_frame(state, s_lo, s_hi, lo, hi, time) &
                                       bind(C,name='transform_to_inertial_frame')

  use meth_params_module, only: NVAR, URHO, UMX, UMZ
  use wdmerger_util_module, only: inertial_velocity
  use castro_util_module, only: position

  implicit none

  integer          :: lo(3), hi(3)
  integer          :: s_lo(3), s_hi(3)
  double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
  double precision :: time

  double precision :: loc(3), vel(3)
  integer          :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           loc = position(i,j,k)
           vel = state(i,j,k,UMX:UMZ) / state(i,j,k,URHO)

           state(i,j,k,UMX:UMZ) = state(i,j,k,URHO) * inertial_velocity(loc, vel, time)

        enddo
     enddo
  enddo

end subroutine transform_to_inertial_frame



! Calculate the second time derivative of the quadrupole moment tensor,
! according to the formula in Equation 6.5 of Blanchet, Damour and Schafer 1990.
! It involves integrating the mass distribution and then taking the symmetric 
! trace-free part of the tensor. We can do the latter operation here since the 
! integral is a linear operator and each part of the domain contributes independently.

subroutine quadrupole_tensor_double_dot(rho, r_lo, r_hi, &
                                        xmom, px_lo, px_hi, ymom, py_lo, py_hi, zmom, pz_lo, pz_hi, &
                                        gravx, gx_lo, gx_hi, gravy, gy_lo, gy_hi, gravz, gz_lo, gz_hi, &
                                        vol, vo_lo, vo_hi, &
                                        lo, hi, dx, time, Qtt) bind(C,name='quadrupole_tensor_double_dot')

  use bl_constants_module, only: ZERO, THIRD, HALF, ONE, TWO, M_PI
  use prob_params_module, only: center, dim
  use wdmerger_util_module, only: inertial_rotation, inertial_velocity
  use castro_util_module, only: position

  implicit none

  integer         , intent(in   ) :: r_lo(3), r_hi(3)
  integer         , intent(in   ) :: px_lo(3), px_hi(3)
  integer         , intent(in   ) :: py_lo(3), py_hi(3)
  integer         , intent(in   ) :: pz_lo(3), pz_hi(3)
  integer         , intent(in   ) :: gx_lo(3), gx_hi(3)
  integer         , intent(in   ) :: gy_lo(3), gy_hi(3)
  integer         , intent(in   ) :: gz_lo(3), gz_hi(3)
  integer         , intent(in   ) :: vo_lo(3), vo_hi(3)

  double precision, intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
  double precision, intent(in   ) :: xmom(px_lo(1):px_hi(1),px_lo(2):px_hi(2),px_lo(3):px_hi(3))
  double precision, intent(in   ) :: ymom(py_lo(1):py_hi(1),py_lo(2):py_hi(2),py_lo(3):py_hi(3))
  double precision, intent(in   ) :: zmom(pz_lo(1):pz_hi(1),pz_lo(2):pz_hi(2),pz_lo(3):pz_hi(3))
  double precision, intent(in   ) :: gravx(gx_lo(1):gx_hi(1),gx_lo(2):gx_hi(2),gx_lo(3):gx_hi(3))
  double precision, intent(in   ) :: gravy(gy_lo(1):gy_hi(1),gy_lo(2):gy_hi(2),gy_lo(3):gy_hi(3))
  double precision, intent(in   ) :: gravz(gz_lo(1):gz_hi(1),gz_lo(2):gz_hi(2),gz_lo(3):gz_hi(3))
  double precision, intent(in   ) :: vol(vo_lo(1):vo_hi(1),vo_lo(2):vo_hi(2),vo_lo(3):vo_hi(3))
  integer         , intent(in   ) :: lo(3), hi(3)
  double precision, intent(in   ) :: dx(3), time
  double precision, intent(inout) :: Qtt(3,3)

  integer          :: i, j, k, l, m
  double precision :: r(3), pos(3), vel(3), g(3), rhoInv, dm
  double precision :: dQtt(3,3)

  dQtt(:,:) = ZERO

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           r = position(i,j,k) - center

           if (rho(i,j,k) > ZERO) then
              rhoInv = ONE / rho(i,j,k)
           else
              rhoInv = ZERO
           endif

           ! Account for rotation, if there is any. These will leave 
           ! r and vel and changed, if not.

           pos = inertial_rotation(r, time)

           ! For constructing the velocity in the inertial frame, we need to 
           ! account for the fact that we have rotated the system already, so that 
           ! the r in omega x r is actually the position in the inertial frame, and 
           ! not the usual position in the rotating frame. It has to be on physical 
           ! grounds, because for binary orbits where the stars aren't moving, that 
           ! r never changes, and so the contribution from rotation would never change.
           ! But it must, since the motion vector of the stars changes in the inertial 
           ! frame depending on where we are in the orbit.

           vel(1) = xmom(i,j,k) * rhoInv
           vel(2) = ymom(i,j,k) * rhoInv
           vel(3) = zmom(i,j,k) * rhoInv

           vel = inertial_velocity(pos, vel, time)

           g(1) = gravx(i,j,k)
           g(2) = gravy(i,j,k)
           g(3) = gravz(i,j,k)

           ! We need to rotate the gravitational field to be consistent with the rotated position.

           g = inertial_rotation(g, time)

           ! Absorb the factor of 2 outside the integral into the zone mass, for efficiency.

           dm = TWO * rho(i,j,k) * vol(i,j,k)

           if (dim .eq. 3) then

              do m = 1, 3
                 do l = 1, 3
                    dQtt(l,m) = dQtt(l,m) + dM * (vel(l) * vel(m) + pos(l) * g(m))
                 enddo
              enddo

           else

              ! For axisymmetric coordinates we need to be careful here.
              ! We want to calculate the quadrupole tensor in terms of
              ! Cartesian coordinates but our coordinates are cylindrical (R, z).
              ! What we can do is to first express the Cartesian coordinates
              ! as (x, y, z) = (R cos(phi), R sin(phi), z). Then we can integrate
              ! out the phi coordinate for each component. The off-diagonal components
              ! all then vanish automatically. The on-diagonal components xx and yy
              ! pick up a factor of cos**2(phi) which when integrated from (0, 2*pi)
              ! yields pi. Note that we're going to choose that the cylindrical z axis
              ! coincides with the Cartesian x-axis, which is our default choice.

              ! We also need to then divide by the volume by 2*pi since
              ! it has already been integrated out.

              dm = dm / (TWO * M_PI)

              dQtt(1,1) = dQtt(1,1) + dm * (TWO * M_PI) * (vel(2)**2 + pos(2) * g(2))
              dQtt(2,2) = dQtt(2,2) + dm * M_PI * (vel(1)**2 + pos(1) * g(1))
              dQtt(3,3) = dQtt(3,3) + dm * M_PI * (vel(1)**2 + pos(1) * g(1))

           endif

        enddo
     enddo
  enddo

  ! Now take the symmetric trace-free part of the quadrupole moment.
  ! The operator is defined in Equation 6.7 of Blanchet et al. (1990):
  ! STF(A^{ij}) = 1/2 A^{ij} + 1/2 A^{ji} - 1/3 delta^{ij} sum_{k} A^{kk}.

  do l = 1, 3
     do m = 1, 3

        Qtt(l,m) = Qtt(l,m) + HALF * dQtt(l,m) + HALF * dQtt(m,l)
        Qtt(l,l) = Qtt(l,l) - THIRD * dQtt(m,m)

     enddo
  enddo

end subroutine quadrupole_tensor_double_dot



! Given the above quadrupole tensor, calculate the strain tensor.

subroutine gw_strain_tensor(h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3, Qtt, time) &
                            bind(C,name='gw_strain_tensor')

  use bl_constants_module, only: ZERO, HALF, ONE, TWO
  use fundamental_constants_module, only: Gconst, c_light, parsec
  use probdata_module, only: gw_dist, axis_1, axis_2, axis_3
  use prob_params_module, only: dim

  implicit none

  double precision, intent(inout) :: h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3
  double precision, intent(in   ) :: Qtt(3,3)
  double precision, intent(in   ) :: time

  integer          :: i, j, k, l, dir
  double precision :: h(3,3), proj(3,3,3,3), delta(3,3), n(3), r
  double precision :: dist(3)

  ! Standard Kronecker delta.

  delta(:,:) = ZERO

  do i = 1, 3
     delta(i,i) = ONE
  enddo

  ! Unit vector for the wave is simply the distance 
  ! vector to the observer normalized by the total distance.
  ! We are going to repeat this process by looking along 
  ! all three coordinate axes.

  do dir = 1, 3

     dist(:) = ZERO
     dist(dir) = gw_dist

     r = sqrt(sum(dist**2))

     n(:) = dist(:) / r

     h = ZERO

     ! Projection operator onto the unit vector n.

     do l = 1, 3
        do k = 1, 3
           do j = 1, 3
              do i = 1, 3
                 proj(i,j,k,l) = (delta(i,k) - n(i) * n(k)) * (delta(j,l) - n(j) * n(l)) &
                               - HALF * (delta(i,j) - n(i) * n(j)) * (delta(k,l) - n(k) * n(l))
              enddo
           enddo
        enddo
     enddo

     ! Now we can calculate the strain tensor.

     do l = 1, 3
        do k = 1, 3
           do j = 1, 3
              do i = 1, 3
                 h(i,j) = h(i,j) + proj(i,j,k,l) * Qtt(k,l)
              enddo
           enddo
        enddo
     enddo

     ! Finally multiply by the coefficients.

     r = r * parsec * 1d3 ! Convert from kpc to cm.

     h(:,:) = h(:,:) * TWO * Gconst / (c_light**4 * r)

     if (dim .eq. 3) then

        ! If rot_axis == 3, then h_+ = h_{11} = -h_{22} and h_x = h_{12} = h_{21}.
        ! Analogous statements hold along the other axes.

        ! We are adding here so that this calculation makes sense on multiple levels.

        if (dir .eq. axis_1) then

           h_plus_1  = h_plus_1  + h(axis_2,axis_2)
           h_cross_1 = h_cross_1 + h(axis_2,axis_3)

        else if (dir .eq. axis_2) then

           h_plus_2  = h_plus_2  + h(axis_3,axis_3)
           h_cross_2 = h_cross_2 + h(axis_3,axis_1)

        else if (dir .eq. axis_3) then

           h_plus_3  = h_plus_3  + h(axis_1,axis_1)
           h_cross_3 = h_cross_3 + h(axis_1,axis_2)

        endif

     else

        ! In 2D axisymmetric coordinates, enforce that axis_1 is the x-axis,
        ! axis_2 is the y-axis, and axis_3 is the z-axis.

        if (dir .eq. 1) then

           h_plus_1  = h_plus_1  + h(2,2)
           h_cross_1 = h_cross_1 + h(2,3)

        else if (dir .eq. 2) then

           h_plus_2  = h_plus_2  + h(3,3)
           h_cross_2 = h_cross_2 + h(3,1)

        else if (dir .eq. 3) then

           h_plus_3  = h_plus_3  + h(1,1)
           h_cross_3 = h_cross_3 + h(1,2)

        endif

     endif

  enddo

end subroutine gw_strain_tensor



subroutine update_center(time) bind(C,name='update_center')

  use bl_constants_module, only: ZERO
  use probdata_module, only: bulk_velx, bulk_vely, bulk_velz, &
                             center_fracx, center_fracy, center_fracz
  use prob_params_module, only: center, problo, probhi, dim

  implicit none

  double precision, intent(in) :: time

  ! Determine the original location of the center.

  if (dim .eq. 3) then

     center(1) = problo(1) + center_fracx * (probhi(1) - problo(1))
     center(2) = problo(2) + center_fracy * (probhi(2) - problo(2))
     center(3) = problo(3) + center_fracz * (probhi(3) - problo(3))

  else if (dim .eq. 2) then

     center(1) = problo(1)
     center(2) = problo(2) + center_fracz * (probhi(2) - problo(2))
     center(3) = ZERO

  else if (dim .eq. 1) then

     center(1) = problo(1) + center_fracx * (probhi(1) - problo(1))
     center(2) = ZERO
     center(3) = ZERO

  else

     call bl_error("Error: unknown dim in subroutine update_center.")

  endif

  ! Now update using the time passed since the beginning of the simulation.

  center(1) = center(1) + bulk_velx * time
  center(2) = center(2) + bulk_vely * time
  center(3) = center(3) + bulk_velz * time

end subroutine update_center



! Updates the CASTRO rotational period.

subroutine set_period(period) bind(C,name='set_period')

  use meth_params_module, only: rot_period

  implicit none

  double precision :: period

  rot_period = period

end subroutine set_period



! Returns the CASTRO rotational period.

subroutine get_period(period) bind(C,name='get_period')

  use meth_params_module, only: rot_period

  implicit none

  double precision :: period

  period = rot_period

end subroutine get_period



! Returns the CASTRO rotation frequency vector.

subroutine get_omega_vec(omega_in, time) bind(C,name='get_omega_vec')

  use rotation_frequency_module, only: get_omega

  implicit none

  double precision :: omega_in(3), time

  omega_in = get_omega(time)

end subroutine get_omega_vec



! Updates the global extrema.

subroutine set_extrema(T_max, rho_max, ts_te_max) bind(C,name='set_extrema')

  use probdata_module, only: T_global_max, rho_global_max, ts_te_global_max

  implicit none

  double precision, intent(in) :: T_max, rho_max, ts_te_max

  T_global_max     = T_max
  rho_global_max   = rho_max
  ts_te_global_max = ts_te_max

end subroutine set_extrema



! Retrieves the global extrema.

subroutine get_extrema(T_max, rho_max, ts_te_max) bind(C,name='get_extrema')

  use probdata_module, only: T_global_max, rho_global_max, ts_te_global_max

  implicit none

  double precision, intent(inout) :: T_max, rho_max, ts_te_max

  T_max     = T_global_max
  rho_max   = rho_global_max
  ts_te_max = ts_te_global_max

end subroutine get_extrema



! Returns whether the simulation is done.

subroutine get_job_status(jobDoneStatus) bind(C,name='get_job_status')

  use probdata_module, only: jobIsDone

  implicit none

  integer, intent(inout) :: jobDoneStatus

  if (jobIsDone) then
     jobDoneStatus = 1
  else
     jobDoneStatus = 0
  endif

end subroutine get_job_status



! Sets whether the simulation is done.

subroutine set_job_status(jobDoneStatus) bind(C,name='set_job_status')

  use probdata_module, only: jobIsDone

  implicit none

  integer, intent(in) :: jobDoneStatus

  if (jobDoneStatus == 1) then
     jobIsDone = .true.
  else
     jobIsDone = .false.
  endif

end subroutine set_job_status



! Gets whether the relaxation is done.

subroutine get_relaxation_status(relaxation_status) bind(C,name='get_relaxation_status')

  use probdata_module, only: relaxation_is_done

  implicit none

  integer, intent(inout) :: relaxation_status

  relaxation_status = relaxation_is_done

end subroutine get_relaxation_status



! Sets whether the relaxation is done.

subroutine set_relaxation_status(relaxation_status) bind(C,name='set_relaxation_status')

  use probdata_module, only: relaxation_is_done

  implicit none

  integer, intent(in) :: relaxation_status

  relaxation_is_done = relaxation_status

end subroutine set_relaxation_status



! Retrieve the total energy array.

subroutine get_total_ener_array(ener_array_in) bind(C,name='get_total_ener_array')

  use probdata_module, only: num_previous_ener_timesteps, total_ener_array

  implicit none

  double precision, intent(inout) :: ener_array_in(num_previous_ener_timesteps)

  ener_array_in(:) = total_ener_array(:)

end subroutine get_total_ener_array



! Set the total energy array.

subroutine set_total_ener_array(ener_array_in) bind(C,name='set_total_ener_array')

  use probdata_module, only: num_previous_ener_timesteps, total_ener_array

  implicit none

  double precision, intent(in) :: ener_array_in(num_previous_ener_timesteps)

  total_ener_array(:) = ener_array_in(:)

end subroutine set_total_ener_array



! Set the number of ignited zones.

subroutine set_num_zones_ignited(num_zones, lev) bind(C,name='set_num_zones_ignited')

  use probdata_module, only: num_zones_ignited, ignition_level

  integer, intent(in) :: num_zones, lev

  num_zones_ignited = num_zones
  ignition_level = lev

end subroutine set_num_zones_ignited



! Get the number of ignited zones.

subroutine get_num_zones_ignited(num_zones, lev) bind(C,name='get_num_zones_ignited')

  use probdata_module, only: num_zones_ignited, ignition_level

  integer, intent(inout) :: num_zones, lev

  num_zones = num_zones_ignited
  lev = ignition_level

end subroutine get_num_zones_ignited
