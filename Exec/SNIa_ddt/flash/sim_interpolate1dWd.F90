!!
!! Dean M. Townsley 2009
!!
!! This subroutine interpolates within a 1-dimensional WD profile
!! by averaging over a cell width.  The WD profile, stored in the
!! Simulation Unit's module data area, is averaged over an interval of
!! size "dr" at radius "radius".  Returned are the density, temperature,
!! and c12 and ne22 mass fractions.  Temperature averaging is just
!! multiplicatively mass-weighted for lack of a better thing to do.
!!
!! Averaging over an interval is used to allow a graceful mapping onto
!! the unrefined grid during initial refinement.

subroutine sim_interpolate1dWd( radius, dr, dens, temp, xc12, xne22)

  use Simulation_data, ONLY : sim_wd_dens_tab, sim_wd_temp_tab, sim_wd_c12_tab, &
                              sim_wd_ne22_tab, sim_wd_dr_inv, sim_wd_npnts, &
                              sim_densFluff, sim_tempFluff, sim_xc12Fluff, &
                              sim_xne22Fluff
  implicit none

  real, intent(in) :: radius, dr
  real, intent(out):: dens, temp, xc12, xne22

  real :: left_edge, right_edge
  real :: shell_left, shell_right, shellvol
  real :: vol, mass, masstemp, massc12, massne22
  integer :: imin, imax, i

!  print *, radius, dr

  ! want to construct a mapping that is direct if the grids are the same, and
  ! otherwise averages the cells in a reasonable way
  ! to do this we assume an evenly spaced grid for which the value
  !  of each variable is constant on the intervals

  ! calculate boundaries of averaging region in grid index coordinates
  left_edge = (radius-0.5*dr) * sim_wd_dr_inv
  right_edge = (radius+0.5*dr) * sim_wd_dr_inv

  ! if below bounds we redefine our averaging region to cover
  ! the last point on  that end
  if (floor(right_edge) <= 0) then
     left_edge = 0.0
     right_edge = 1.0
  endif
  ! above use fluff
  if (floor(left_edge) >= sim_wd_npnts) then
     dens = sim_densFluff
     temp = sim_tempFluff
     xc12 = sim_xc12Fluff
     xne22= sim_xne22Fluff
!  print *, 'fluff', dens, temp, xc12, xne22
     return
  endif

  vol = 0.0
  mass = 0.0
  masstemp = 0.0
  massc12  = 0.0
  massne22 = 0.0
  if ( right_edge > real(sim_wd_npnts) ) then
     shellvol = right_edge**3 - sim_wd_npnts**3
     
     vol = vol + shellvol
     mass = mass + sim_densFluff*shellvol
     masstemp = masstemp + sim_densFluff*sim_tempFluff*shellvol
     massc12  = massc12  + sim_densFluff*sim_xc12Fluff*shellvol
     massne22 = massne22 + sim_densFluff*sim_xne22Fluff*shellvol
  endif

  ! sum through cells in 1-d grid that this region overlaps
  imin = max(0, floor(left_edge))
  imax = min(sim_wd_npnts-1, floor(right_edge))
!  print *, left_edge, right_edge
!  print *, imin, imax
  do i = imin, imax
     ! average over just the portion of this cell which overlaps the averaging region
     shell_left = max(real(i), left_edge)
     shell_right = min(real(i+1), right_edge)
!     print *, shell_left, shell_right
     ! assume 1d profile has spherecial geometry to do mass averages
!     shellvol = (shell_right-shell_left)*0.25*(shell_right+shell_left)**2 ! * 4*pi
     shellvol = shell_right**3-shell_left**3

     vol  = vol  + shellvol
     mass = mass + sim_wd_dens_tab(i+1)*shellvol
     masstemp = masstemp + sim_wd_dens_tab(i+1)*sim_wd_temp_tab(i+1)*shellvol
     massc12  = massc12  + sim_wd_dens_tab(i+1)*sim_wd_c12_tab(i+1)*shellvol
     massne22 = massne22 + sim_wd_dens_tab(i+1)*sim_wd_ne22_tab(i+1)*shellvol
  enddo
  dens = mass/vol
  temp = masstemp/mass
  xc12 = massc12/mass
  xne22 = massne22/mass

!  print *, dens, temp, xc12, xne22
end subroutine
