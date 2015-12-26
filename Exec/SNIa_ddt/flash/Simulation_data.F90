! Dean M. Townley 2009
!
! This is the static module data for the Simulation Unit for the SNIa_ddt setup
! Note that some of these are allocatable, and are allocated by Simulation_init()

module Simulation_data
#include "Flash.h"
#include "Eos.h"

  real,allocatable,dimension(:),save :: sim_wd_dens_tab, sim_wd_temp_tab, sim_wd_c12_tab, sim_wd_ne22_tab
  real, save :: sim_wd_dr_inv
  integer, save :: sim_wd_npnts

  real, save :: sim_densFluff, sim_tempFluff, sim_xc12Fluff, sim_xne22Fluff

  logical, save :: sim_ignite
  real, save :: sim_ignX, sim_ignY, sim_ignZ, sim_ignR

  logical, save :: sim_ignitionFile

  logical, save :: sim_ignMPole, sim_ignMPoleSym
  real, save :: sim_ignMPoleA
  integer, save :: sim_ignMpoleMinL, sim_ignMpoleMaxL, sim_ignMpoleSeed
  integer, save :: sim_geom
  real, allocatable, dimension(:,:), save  :: mp_A, mp_delta

  logical, save :: sim_ignSin
  real, save    :: sim_ignSinN, sim_ignSinA

  real,    save           :: sim_laminarWidth

  real, save :: sim_refFluffDensThresh, sim_refFluffMargin
  real, save :: sim_refNogenEnucThresh, sim_refNogenFldtThresh, sim_refNogenMargin
  integer, save :: sim_refFluffLevel, sim_refNogenLevel

end module Simulation_data
