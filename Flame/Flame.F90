module flame_module

  implicit none

  logical,          save :: initialized = .false.

  ! runtime parameters

  double precision, save :: epsilon_0, epsilon_1, kpp_fact, b
  double precision, save :: initProfileAdjustWidth

  ! some constants for diffusion-reaction equation
  double precision, save :: R_over_s, kappa_over_s, width

  double precision, save :: effsmlrho

  ! Laminar flame speed data
  logical,          save :: fsUseConstFlameSpeed, fsUseTFI
  double precision, save :: fsConstFlameSpeed, fsConstFlameWidth

contains

  subroutine flame_init() bind(C)

    use amrinfo_module, only: amr_level
    use prob_params_module, only: dx_level

    use extern_probin_module, only: fl_epsilon_0, fl_epsilon_1, &
                                    fl_kpp_fact, fl_b, fl_initProfileAdjustWidth, &
                                    fl_fsUseConstFlameSpeed, fl_fsUseTFI, &
                                    fl_fsConstFlameSpeed, fl_fsConstFlameWidth, &
                                    fl_effsmlrho

    implicit none

    double precision :: dx

    if (initialized) return

    epsilon_0 = fl_epsilon_0
    epsilon_1 = fl_epsilon_1
    kpp_fact = fl_kpp_fact
    b = fl_b
    initProfileAdjustWidth = fl_initProfileAdjustWidth

    effsmlrho = fl_effsmlrho

    fsUseConstFlameSpeed = fl_fsUseConstFlameSpeed
    fsUseTFI             = fl_fsUseTFI
    fsConstFlameSpeed    = fl_fsConstFlameSpeed
    fsConstFlameWidth    = fl_fsConstFlameWidth

    dx = minval(dx_level(:,amr_level))

    R_over_s = kpp_fact * 4.0 / fl_b / dx
    kappa_over_s = fl_b * dx / 16.0

    width = b * dx

    initialized = .true.

  end subroutine flame_init
  


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
    f = 0.5 * (1.0 - tanh(x/width/0.5/initProfileAdjustWidth))

  end subroutine Flame_getProfile



  subroutine Flame_getWidth(laminarWidth)

    implicit none

    double precision, intent(out) :: laminarWidth

    laminarWidth=width

  end subroutine Flame_getWidth



  subroutine ca_flame_step(lo, hi, state, s_lo, s_hi, grav, g_lo, g_hi, time, dt) bind(C)

    use bl_constants_module
    use mempool_module
    use meth_params_module, only: NVAR, UFX
    use network
    use prob_params_module, only: dg

    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer :: g_lo(3), g_hi(3)

    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision ::  grav(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),3)
    double precision :: time, dt

    integer :: istat
    integer :: n, bid

    double precision, pointer, dimension(:,:,:) :: flam, flamdot, flamespeed, lapl

    double precision :: f, inv_dt
    integer :: i,j,k
    integer :: sizeI, sizeJ, sizeK

    ! These only need two ghost cells but we'll give them the same size as the state data

    call bl_allocate(flam      , s_lo(1), s_hi(1), s_lo(2), s_hi(2), s_lo(3), s_hi(3))
    call bl_allocate(flamdot   , s_lo(1), s_hi(1), s_lo(2), s_hi(2), s_lo(3), s_hi(3))
    call bl_allocate(flamespeed, s_lo(1), s_hi(1), s_lo(2), s_hi(2), s_lo(3), s_hi(3))
    call bl_allocate(lapl      , s_lo(1), s_hi(1), s_lo(2), s_hi(2), s_lo(3), s_hi(3))

    inv_dt = ONE / dt

    ! Copy the flame scalar into a separate array

    flam = state(:,:,:,UFX+UFLAM-1)

    call fl_flameSpeed(lo, hi, state, flamespeed, s_lo, s_hi, grav, g_lo, g_hi, 2)

    do k = lo(3) - 2 * dg(3), hi(3) + 2 * dg(3)
       do j = lo(2) - 2 * dg(2), hi(2) + 2 * dg(2)
          do i = lo(1) - 2 * dg(1), hi(1) + 2 * dg(1)
             f = flam(i,j,k)
             flam(i,j,k) = f + dt*R_over_s*flamespeed(i,j,k)*(f-epsilon_0)*(ONE+epsilon_1-f)
          enddo
       enddo
    enddo

    ! 1 specifies the step size should be 1 grid cell
    call fl_laplacian(lo, hi, lapl, flam, s_lo, s_hi, 1)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             flam(i,j,k) = max(ZERO, min(ONE, flam(i,j,k) + dt*kappa_over_s*flamespeed(i,j,k)*lapl(i,j,k) ) )
             flamdot(i,j,k) = (flam(i,j,k) - state(i,j,k,UFX+UFLAM-1))*inv_dt
             state(i,j,k,UFX+UFLAM-1) = flam(i,j,k)
             state(i,j,k,UFX+UFLDT-1) = flamdot(i,j,k)
          enddo
       enddo
    enddo

    call bl_deallocate(lapl)
    call bl_deallocate(flamespeed)
    call bl_deallocate(flam)
    call bl_deallocate(flamdot)

  end subroutine ca_flame_step



  subroutine fl_laplacian(lo, hi, lapl, flam, s_lo, s_hi, h)

    use amrinfo_module, only: amr_level
    use prob_params_module, only: coord_type, dx_level, dim
    use castro_util_module, only: position
    use bl_constants_module
    use mempool_module

    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    
    double precision :: lapl(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision :: flam(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    integer :: h

    double precision :: inv_12_dx2, inv_12_dy2, inv_12_dz2
    double precision :: inv_12_dr2, inv_12_dr, inv_12_dtheta, inv_12_dtheta2, inv_12_dphi2
    double precision :: d2, d1
    double precision :: inv_r, inv_r2, two_over_r, ctan, inv_sin2
    integer :: i, j, k
    integer :: h2
    double precision :: dx(3), loc(3)

    dx = dx_level(:,amr_level)

    h2 = 2*h

    select case (coord_type)

    case (0) ! Cartesian

       ! seems like any self-respecting optimizing compiler would be able to
       ! pull these out of the loop but we won't trust that
       inv_12_dx2 = TWELFTH / dx(1)**2
       if (dim >= 2) inv_12_dy2 = TWELFTH / dx(2)**2
       if (dim == 3) inv_12_dz2 = TWELFTH / dx(3)**2

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                lapl(i,j,k) = ( -flam(i-h2,j,k) + 16*flam(i-h,j,k) -30*flam(i,j,k) &
                     + 16*flam(i+h,j,k) - flam(i+h2,j,k) ) * inv_12_dx2
                if (dim >= 2) then
                   lapl(i,j,k) = lapl(i,j,k) + &
                        ( -flam(i,j-h2,k) + 16*flam(i,j-h,k) -30*flam(i,j,k) &
                        + 16*flam(i,j+h,k) - flam(i,j+h2,k) ) * inv_12_dy2
                endif
                if (dim == 3) then
                   lapl(i,j,k) = lapl(i,j,k) + &
                        ( -flam(i,j,k-h2) + 16*flam(i,j,k-h) -30*flam(i,j,k) &
                        + 16*flam(i,j,k+h) - flam(i,j,k+h2) ) * inv_12_dz2
                endif

             enddo
          enddo
       enddo

    case (1) ! Cylindrical

       inv_12_dr  = TWELFTH / dx(1)
       inv_12_dr2 = inv_12_dr / dx(1)

       if (dim >= 2) inv_12_dz2 = TWELFTH / dx(2)**2
       if (dim == 3) inv_12_dtheta2 = TWELFTH / dx(3)**2

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k)
                
                inv_r = ONE / loc(1)
                inv_r2 = inv_r**2

                d2 = ( -flam(i-h2,j,k) + 16*flam(i-h,j,k) -30*flam(i,j,k) &
                     + 16*flam(i+h,j,k) - flam(i+h2,j,k) ) * inv_12_dr2
                d1 = ( flam(i-h2,j,k) -8*flam(i-h,j,k)  +8*flam(i+h,j,k) -flam(i+h2,j,k) )*inv_12_dr

                lapl(i,j,k) = d2 + inv_r*d1

                if (dim >= 2) then
                   d2 = ( -flam(i,j-h2,k) + 16*flam(i,j-h,k) -30*flam(i,j,k) &
                        + 16*flam(i,j+h,k) - flam(i,j+h2,k) ) * inv_12_dz2
                   lapl(i,j,k) = lapl(i,j,k) +  d2
                endif
                if (dim == 3) then
                   d2 = ( -flam(i,j,k-h2) + 16*flam(i,j,k-h) -30*flam(i,j,k) &
                        + 16*flam(i,j,k+h) - flam(i,j,k+h2) ) * inv_12_dtheta2
                   lapl(i,j,k) = lapl(i,j,k) +  inv_r2*d2
                endif

             enddo
          enddo
       enddo

    case (2) ! Spherical

       inv_12_dr  = TWELFTH / dx(1)
       inv_12_dr2 = inv_12_dr / dx(1)

       if (dim >= 2) then
          inv_12_dtheta = TWELFTH / dx(2)
          inv_12_dtheta2 = inv_12_dtheta / dx(2)
       endif

       if (dim == 3) then
          inv_12_dphi2 = TWELFTH / dx(3)**2
       endif

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k)
                
                two_over_r = TWO / loc(1)
                inv_r2 = ONE / loc(1)**2

                if (dim >= 2) ctan = ONE / tan(loc(2))
                if (dim == 3) inv_sin2 = ONE + ctan**2                   

                d2 = ( -flam(i-h2,j,k) + 16*flam(i-h,j,k) -30*flam(i,j,k) &
                     + 16*flam(i+h,j,k) - flam(i+h2,j,k) ) * inv_12_dr2
                d1 = ( flam(i-h2,j,k)-8*flam(i-h,j,k)+8*flam(i+h,j,k) -flam(i+h2,j,k) )*inv_12_dr
                lapl(i,j,k) = d2 + two_over_r*d1

                if (dim>=2) then
                   d2 = ( -flam(i,j-h2,k) + 16*flam(i,j-h,k) -30*flam(i,j,k) &
                        + 16*flam(i,j+h,k) - flam(i,j+h2,k) ) * inv_12_dtheta2
                   d1 = ( flam(i,j-h2,k)-8*flam(i,j-h,k)+8*flam(i,j+h,k)-flam(i,j+h2,k) )*inv_12_dtheta
                   lapl(i,j,k) = lapl(i,j,k) + inv_r2*( d2 + ctan*d1 )
                endif
                if (dim==3) then
                   d2 = ( -flam(i,j,k-h2) + 16*flam(i,j,k-h) -30*flam(i,j,k) &
                        + 16*flam(i,j,k+h) - flam(i,j,k+h2) ) * inv_12_dphi2
                   lapl(i,j,k) = lapl(i,j,k) + inv_r2*inv_sin2*d2
                endif
             enddo
          enddo
       enddo

    end select

  end subroutine fl_laplacian



  subroutine Flame_rhJumpReactive(eosData_u, qbar_u, eosData_b, qbar_b, eos_mode)

    use eos_module
    !  use Eos_helmData, ONLY : eos_tol
    use NSE_data, ONLY: NSE_finalAtDens 
    use bl_constants_module, only: ONE, TWO

    implicit none

    type (eos_t)    , intent(inout) :: eosData_u
    double precision, intent(in)    :: qbar_u
    type (eos_t)    , intent(out)   :: eosData_b
    double precision, intent(out)   :: qbar_b
    integer,          intent(in)    :: eos_mode

    integer, parameter :: max_newton = 50
    double precision, parameter :: eos_tol = 1.0d-8

    type (eos_t)     :: eosData
    double precision :: ye, pres_u, hmq_u
    double precision :: dens_n, emq, pres_n, qbar, sumyi, tempguess, edot, yedot
    double precision :: error, dd, dpdd, f, dfdd, dens_n_old

    double precision, parameter :: cgsMeVperAmu = 9.6485e17  

    integer :: niters

    ! calculate thermodynamic info about initial state according to mode argument
    call eos(eos_mode, eosData_u)
    pres_u = eosData_u % p
    ye = eosData_u % zbar / eosData_u % abar
    hmq_u = eosData_u % e + eosData_u % p / eosData_u % rho - cgsMeVperAmu * qbar_u

    ! A bit unusual here
    !  we will use the (dens,emq) table to construct the endpoint of
    ! a constant pressure burn.  Although this is what is stored in
    ! the (pres,hmq) table, due to interpolation accuracy and the fact
    ! that fully burned cells use the (dens,emq) table, we will instead use it
    ! to construct the fully burned state so that there is no "jolt" in
    ! fully burned cells at the beginning of the simulation.

    ! we solve for the density such that
    !  pres = pres_u
    ! with
    !  emq = hmq_u - pres_u / dens

    ! our initial guess
    dens_n = eosData_u % rho

    error = TWO * eos_tol
    niters = 0

    do while ( (error > eos_tol) .and. (niters < max_newton) )

       ! get the nse state information
       emq = hmq_u - pres_u/dens_n
       call NSE_finalAtDens(qbar, sumyi, tempguess, edot, yedot, ye, dens_n, emq)
       ! and pressure
       eosData % rho  = dens_n
       eosData % T    = tempguess
       eosData % xn   = eosData_u % xn
       eosData % e    = emq + cgsMeVperAmu * qbar
       eosData % abar = ONE / sumyi
       eosData % zbar = ye * eosData % abar
       eosData % y_e  = ye
       call eos(eos_input_re, eosData)

       pres_n = eosData % p

       ! derivative
       dd = 1.0d-7 * dens_n
       ! get the nse state information
       emq = hmq_u - pres_u / (dens_n+dd)
       call NSE_finalAtDens(qbar, sumyi, tempguess, edot, yedot, ye, dens_n+dd, emq)

       ! and pressure
       eosData % rho  = dens_n + dd
       eosData % T    = tempguess
       eosData % xn   = eosData_u % xn
       eosData % e    = emq + cgsMeVperAmu * qbar
       eosData % abar = ONE / sumyi
       eosData % zbar = ye * eosData % abar
       eosData % y_e  = ye
       call eos(eos_input_re, eosData)

       dpdd = (eosData % p - pres_n) / dd

       ! function we are zeroing and deriv
       f = pres_n - pres_u
       dfdd = dpdd

       dens_n_old = dens_n
       dens_n = dens_n - f/dfdd

       if (dens_n .lt. effsmlrho) then
          write (6,*) 'small density in nseJump'
          dens_n = 0.5*dens_n_old
       endif

       error = abs( (dens_n-dens_n_old) / dens_n )

       niters = niters + 1

    enddo

    if (niters >= max_newton) then
       write(6,*) 'exceeded number of newton steps in nsejump'
       write(6,*) 'ended with dens ', dens_n
    endif

    ! now fill output information
    emq = hmq_u - pres_u / dens_n
    call NSE_finalAtDens(qbar_b, sumyi, tempguess, edot, yedot, ye, dens_n, emq)

    ! and pressure
    eosData_b % rho  = dens_n
    eosData_b % T    = tempguess
    eosData_b % xn   = eosData_u % xn
    eosData_b % e    = emq + cgsMeVperAmu * qbar_b
    eosData_b % abar = ONE / sumyi
    eosData_b % zbar = ye * eosData_b % abar
    eosData_b % y_e  = ye
    call eos(eos_input_re, eosData_b)

  end subroutine Flame_rhJumpReactive



  subroutine Flame_rhjump(eosData_u, eosData_b, q, s, mode)

    use eos_module
    use network

    implicit none

    type (eos_t),     intent(inout) :: eosData_u
    type (eos_t),     intent(inout) :: eosData_b
    double precision, intent(in   ) :: q, s
    integer,          intent(in   ) :: mode  !! This is the EOS mode

    double precision :: dfd, dft, dgd, dgt, f, g
    double precision :: sq_s_dens, d_inv_dens, determinant_inv
    double precision :: dens_b_old, temp_b_old
    double precision :: dens_u, temp_u, ener_u, pres_u, dens_b, temp_b, ener_b, pres_b
    double precision ::  error
    integer          :: niters

!1   format(5(2X, E10.4))

    ! first calculate the rest of the thermodynamic state of the
    ! unburned material according to the "mode"
    ! we only need dens, pres and eint

    call eos(mode, eosData_u)

    !  write (*,*) '------------------------------------------------------'
    !  write (*,*) ''
    !  write (*,*) 'rhjump: flame speed, heat release:'
    !  write (*,1) s, q
    !  write (*,*) 'rhjump: unburned temp, dens, ener, pres:'
    !  write (*,1) temp_u, dens_u, ener_u, pres_u
    !  write (*,*) ''

    !.. initial (not very good) guess for burned state:

    eosData_b % rho = eosData_u % rho
    eosData_b % e   = eosData_u % e + q
    ! need a guess temperature
    eosData_b % T = eosData_u % T

    call eos(eos_input_rt, eosData_b)

    ! some renames for readability
    dens_u = eosData_u % rho
    temp_u = eosData_u % T
    ener_u = eosData_u % e
    pres_u = eosData_u % p
    sq_s_dens = (s*dens_u)**2
    dens_b = eosData_b % rho
    temp_b = eosData_b % T
    ener_b = eosData_b % e
    pres_b = eosData_b % p

    ! 2d newton loop to find dens_b and temp_b which solve
    !    pres_b = pres_u - (dens_u*s)**2*(1/dens_b - 1/dens_u)
    !    ener_b = ener_u + q - 0.5 ( pres_b+ pres_u)*(1/dens_b-1/dens_u)

    !  the ordering of evaluations and tests here is a little weird -- should be fixed (DMT 2007/4/14)
    error = 1
    niters = 0
    do while (error > 1.e-8 .and. niters < 100)

       d_inv_dens = 1./dens_b - 1./dens_u

       f = pres_b - pres_u + sq_s_dens * d_inv_dens 
       g = ener_b - ener_u - q + 0.5*(pres_b + pres_u) * d_inv_dens

       dfd = eosData_b % dpdr - sq_s_dens/dens_b**2
       dft = eosData_b % dpdt
       dgd = eosData_b % dedr + 0.5*d_inv_dens*eosData_b % dpdr - 0.5*(pres_b + pres_u) / dens_b**2
       dgt = eosData_b % dedT + 0.5*d_inv_dens*eosData_b % dpdT

       determinant_inv = 1./(dfd*dgt - dft*dgd)

       dens_b_old = dens_b
       temp_b_old = temp_b

       dens_b = dens_b - (f*dgt - g*dft) * determinant_inv
       temp_b = temp_b + (f*dgd - g*dfd) * determinant_inv

       if (dens_b .lt. effsmlrho) then

          dens_b = 0.5*dens_b_old
          temp_b = temp_b_old

       elseif (temp_b .lt. temp_u) then

          temp_b = temp_u
          dens_b = dens_b_old

       endif

       ! un-rename new values
       eosData_b % rho = dens_b
       eosData_b % T   = temp_b

       !     write (6,*) 'in loop call', dens_b, temp_b, mfrac_b

       call eos(eos_input_rt, eosData_b)

       ! rename
       pres_b=eosData_b % p
       ener_b=eosData_b % e

       error = abs(f/pres_u) + abs(g/ener_u)

       !     write(*,1) temp_b, dens_b, ener_b, pres_b, error

       niters = niters + 1

    enddo

    !  write (*,*) ''
    !  write (*,*) 'rhjump: burned temp, dens, ener, pres:'
    !  write (*,1) temp_b, dens_b, ener_b, pres_b
    !  write (*,*) 'rhjump:   niters, error:', niters, error
    !  write (*,*) ''
    !  write (*,*) '-------------------------------------------------'
    if (niters >= 100) write (6,*) 'rhjump did not converge'

    !------------------------------------------------------------------

  end subroutine Flame_rhjump



  subroutine fl_flameSpeed(lo, hi, state, flamespeed, s_lo, s_hi, grav, g_lo, g_hi, nlayers)

    use meth_params_module, ONLY : NVAR, UFX
    use network

    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer :: g_lo(3), g_hi(3)

    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3),NVAR)
    double precision :: flamespeed(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    double precision :: grav(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),3)
    integer          :: nlayers

    double precision, pointer, dimension(:,:,:) :: dens, flamewidth
    integer :: sizeX, sizeY, sizeZ

    integer :: comp_lo(3), comp_hi(3)

    flamespeed(:,:,:) = fsConstFlameSpeed
    state(:,:,:,UFX+UFLSP-1) = fsConstFlameSpeed

  end subroutine fl_flameSpeed

end module flame_module
