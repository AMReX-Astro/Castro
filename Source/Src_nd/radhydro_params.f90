module radhydro_params_module

  implicit none

  integer, save :: QRADVAR, qrad, qradhi, qptot, qreitot 
  integer, save :: fspace_type
  logical, save :: do_inelastic_scattering
  logical, save :: comoving

  double precision, save :: flatten_pp_threshold = -1.d0

contains
  
  subroutine get_qradvar(qradvar_in) bind(C, name="get_qradvar")

    implicit none

    integer, intent(inout) :: qradvar_in

    qradvar_in = QRADVAR

  end subroutine get_qradvar

end module radhydro_params_module


subroutine ca_init_radhydro_pars(fsp_type_in, do_is_in, com_in,fppt)
  use meth_params_module, only : QVAR
  use rad_params_module, only : ngroups
  use radhydro_params_module
  implicit none
  integer, intent(in) :: fsp_type_in, do_is_in, com_in
  double precision, intent(in) :: fppt

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
  
  do_inelastic_scattering = (do_is_in .ne. 0)

  if (com_in .eq. 1) then
     comoving = .true.
  else if (com_in .eq. 0) then
     comoving = .false.
  else
     print *, "Wrong value for comoving", fspace_type
     stop        
  end if
  
  flatten_pp_threshold = fppt

end subroutine ca_init_radhydro_pars
