module radhydro_params_module

  integer, save :: QRADVAR, qrad, qradhi, qptot, qreitot 
  integer, save :: fspace_type
  logical, save :: comoving

  logical, save :: first_order_hydro = .false.
  double precision, save :: flatten_pp_threshold = -1.d0

contains

  subroutine init_radhydro_params(fsp_type_in, com_in)
    use meth_params_module, only : QVAR
    use rad_params_module, only : ngroups

    implicit none

    integer, intent(in) :: fsp_type_in, com_in

    qptot  = QVAR+1
    qreitot = QVAR+2
    qrad   = QVAR+3
    qradhi = qrad+ngroups-1

    QRADVAR = QVAR + 2 + ngroups

    if (ngroups .eq. 1) then
       fspace_type = 1
    else
       fspace_type = fsp_type_in
    end if

    if (fsp_type_in .ne. 1 .and. fsp_type_in .ne. 2) then
       print *, "Unknown fspace_type", fspace_type
       stop 
    end if

    if (com_in .eq. 1) then
       comoving = .true.
    else if (com_in .eq. 0) then
       comoving = .false.
    else
       print *, "Wrong value for comoving", fspace_type
       stop        
    end if

  end subroutine init_radhydro_params

end module radhydro_params_module

subroutine ca_init_more_radhydro_pars(foh,fppt)
  use radhydro_params_module
  implicit none
  integer, intent(in) :: foh
  double precision, intent(in) :: fppt
  if (foh .eq. 1) then
     first_order_hydro = .true.
  else
     first_order_hydro = .false.
  end if
  flatten_pp_threshold = fppt
end subroutine ca_init_more_radhydro_pars
