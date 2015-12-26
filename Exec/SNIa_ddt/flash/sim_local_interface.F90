module sim_local_interface

  interface sim_interpolate1dWd
     subroutine sim_interpolate1dWd( radius, dr, dens, temp, xc12, xne22)
     implicit none
     real, intent(in) :: radius, dr
     real, intent(out):: dens, temp, xc12, xne22
     end subroutine
  end interface

  interface sim_LCGRandomIterate
    subroutine sim_LCGRandomIterate(state)
    implicit none
    integer, intent(inout) :: state
    end subroutine
  end interface

end module

