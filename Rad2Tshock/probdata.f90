module probdata_module
  
  ! Ensman test variables -- we set the defaults here
  double precision, save :: rho0 = 5.4588672836d-13
  double precision, save :: T0 = 100.d0
  double precision, save :: v0 = 235435.371882d0
  double precision, save :: rho1 = 1.24793794736d-12
  double precision, save :: T1 = 207.756999533d0
  double precision, save :: v1 = 102986.727159d0
  
  namelist /fortin/ rho0, T0, v0, rho1, T1, v1

  ! for convenience
  double precision, save :: xmin, xmax, ymin, ymax, zmin, zmax
  
end module probdata_module
