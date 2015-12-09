
! This module stores the runtime parameters and integer names for 
! tagging zones for refinement.
!
! These parameters are initialized in get_tagging_params().

module tagging_params_module

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

end module tagging_params_module
