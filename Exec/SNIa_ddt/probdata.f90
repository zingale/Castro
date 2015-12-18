module probdata_module

  ! 1D model file name
  character (len=80), save :: model_name = "wd_hse.dat"

  ! Ambient material
  double precision, save :: dens_fluff  = 1.0d-3
  double precision, save :: temp_fluff  = 3.0d7
  double precision, save :: xc12_fluff  = 0.5d0
  double precision, save :: xne22_fluff = 0.02d0

  double precision, save :: rep_ne_frac = 0.03d0

  ! Ignition data

  logical,          save :: sim_ignite = .false.
  double precision, save :: sim_ignX = 0.0d0
  double precision, save :: sim_ignY = 0.0d0
  double precision, save :: sim_ignZ = 0.0d0
  double precision, save :: sim_ignR = 0.0d0

  logical, save :: sim_ignitionFile = .false.    
  character (len=80), save :: ignition_file_name = "ignition_points.dat"

  ! Multipole ignition data  
  
  logical,          save :: sim_ignMPole = .false.
  logical,          save :: sim_ignMPoleSym = .false.
  double precision, save :: sim_ignMPoleA = 0.3d7
  integer,          save :: sim_ignMPoleMinL = 4
  integer,          save :: sim_ignMPoleMaxL = 16
  integer,          save :: sim_ignMPoleSeed = 169332

  ! Sinusoidal ignition surface perturbation

  logical,          save :: sim_ignSin = .false.
  double precision, save :: sim_ignSinN = 6.0d0
  double precision, save :: sim_ignSinA = 30.0d5
  
  ! Problem tagging criteria for custom refinement

  ! refinement threshold for derefinement in Fluff
  ! below this density a maximum refinement level of refFluffLevel is enforced

  double precision, save :: sim_refFluffDensThresh = 1.0d5
  
  ! Padding margin for thresholding (in fractional units)
  ! in margin above threshold ...

  double precision, save :: sim_refFluffMargin = 0.2d0  
  
  ! Maximum refinement level of fluff

  integer,          save :: sim_refFluffLevel = 1

  double precision, save :: sim_refNogenEnucThresh = 1.0d16
  double precision, save :: sim_refNogenFldtThresh = 0.02d0
  double precision, save :: sim_refNogenMargin = 0.5d0
  integer,          save :: sim_refNogenLevel = 2

  double precision, save :: sim_laminarWidth  
  
  ! Other parameters
  
  logical,          save :: useBurn = .true.
  logical,          save :: bn_thermalReact = .false.


  ! Local data to save

  double precision, save :: sim_wd_dr_inv
  integer,          save :: sim_wd_npnts

  double precision, allocatable, dimension(:,:), save  :: mp_A, mp_delta

end module probdata_module
