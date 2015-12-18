module probdata_module

  ! 1D model file name
  character (len=80), save :: model_name = "wd_hse.dat"

  ! Ambient material
  double precision, save :: dens_fluff  = 1.0d-3
  double precision, save :: temp_fluff  = 3.0d7
  double precision, save :: xc12_fluff  = 0.5d0
  double precision, save :: xne22_fluff = 0.02d0

  double precision, save :: rep_ne_frac = 0.03d0
  
end module probdata_module
