!  Dean Townsley 2008

module flame_module

  implicit none

  ! runtime parameters
  logical, save :: fl_useFlame
  real, save :: fl_epsilon_0, fl_epsilon_1, fl_kpp_fact, fl_b
  real, save :: fl_initProfileAdjustWidth

  ! some constants for diffusion-reaction equation
  real, save :: fl_R_over_s, fl_kappa_over_s, fl_width

  real, save :: fl_effsmlrho

contains

  subroutine Flame_getProfile(x, f)

    implicit none

    double precision, intent(in)  :: x
    double precision, intent(out) :: f

    ! This is an approximate profile form based on widths determined in
    ! Vladimirova et al.  Here the "width" is approximately the
    ! distance between where phi=0.12 and 0.88.  Over twice width, phi
    ! goes from 0.02 to 0.98.
    ! The fl_initProfileAdjustmentWidth is to allow compatibility with
    ! slight variations if necessary

    ! tanh is slow, but no need to be super efficient here since this
    ! should only be called at init
    f = 0.5 * (1.0 - tanh(x/fl_width/0.5/fl_initProfileAdjustWidth))

  end subroutine Flame_getProfile



  subroutine Flame_getWidth(laminarWidth)

    implicit none

    double precision, intent(out) :: laminarWidth

    laminarWidth=fl_width

  end subroutine Flame_getWidth

end module flame_module
