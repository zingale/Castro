! calculate post-flame NSE state given unburned material with properties
! eosData_u and qbar_u
! eosData_u is updated once with Eos() called with eos_mode, allowing other
! combinations of parameters other than density and temperature of the unburned
! material to be specified
!
! Dean Townsley 2007,2008

! note that this routine assumes that eosData_u(EOS_ABAR) and eosData_u(EOS_ZBAR)
! are already set appropriately.  i.e. this is a non-species-based routine

subroutine Flame_rhJumpReactive(eosData_u, qbar_u, eosData_b, qbar_b, eos_mode)

  use eos_module
  use flame_module, only: fl_effsmlrho
!  use Eos_helmData, ONLY : eos_tol

  use NSE_data, ONLY: NSE_finalAtDens 
  
  implicit none

  type (eos_t)    , intent(inout) :: eosData_u
  double precision, intent(in)    :: qbar_u
  type (eos_t)    , intent(out)   :: eosData_b
  double precision, intent(out)   :: qbar_b
  integer,          intent(in)    :: eos_mode

  integer, parameter :: max_newton = 50
  double precision, parameter :: eos_tol = 1.0d-8
  
  type (eos_t)     :: eosData
  double precision ::  ye, pres_u, hmq_u
  double precision ::  dens_n, emq, pres_n, qbar, sumyi, tempguess, edot, yedot
  double precision ::  error, dd, dpdd, f, dfdd, dens_n_old

  integer :: niters

  ! calculate thermodynamic info about initial state according to mode argument
  call eos(eos_mode, eosData_u)
  pres_u = eosData_u % p
  ye = eosData_u % zbar / eosData_u % abar
  hmq_u = eosData_u % e + eosData_u % p / eosData_u % rho - 9.6485e17*qbar_u

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

  error = 2.0* eos_tol
  niters = 0
  
  do while ( (error > eos_tol) .and. (niters < max_newton) )

     ! get the nse state information
     emq = hmq_u - pres_u/dens_n
     call NSE_finalAtDens(qbar, sumyi, tempguess, edot, yedot, ye, dens_n, emq)
     ! and pressure
     eosData % rho = dens_n
     eosData % T = tempguess
     eosData % e = emq + 9.6485e17*qbar
     eosData % abar = 1.e0/sumyi
     eosData % zbar = ye * eosData % abar
     call eos(eos_input_re, eosData)

     pres_n = eosData % p

     ! derivative
     dd = 1.0e-7*dens_n
     ! get the nse state information
     emq = hmq_u - pres_u/(dens_n+dd)
     call NSE_finalAtDens(qbar, sumyi, tempguess, edot, yedot, ye, dens_n+dd, emq)

     ! and pressure
     eosData % rho = dens_n + dd
     eosData % T = tempguess
     eosData % e = emq + 9.6485e17*qbar
     eosData % abar = 1.e0/sumyi
     eosData % zbar = ye * eosData % abar
     call eos(eos_input_re, eosData)

     dpdd = (eosData % p - pres_n)/dd

     ! function we are zeroing and deriv
     f = pres_n - pres_u
     dfdd = dpdd

     dens_n_old = dens_n
     dens_n = dens_n - f/dfdd

     if (dens_n .lt. fl_effsmlrho) then
        write (6,*) 'small density in nseJump'
        dens_n = 0.5*dens_n_old
     endif

     error = abs( (dens_n-dens_n_old)/dens_n )

     niters = niters + 1

  enddo
 
  if (niters >= max_newton) then
     write(6,*) 'exceeded number of newton steps in nsejump'
     write(6,*) 'ended with dens ', dens_n
  endif

  ! now fill output information
  emq = hmq_u - pres_u/dens_n
  call NSE_finalAtDens(qbar_b, sumyi, tempguess, edot, yedot, ye, dens_n, emq)

  ! and pressure
  eosData_b % rho = dens_n
  eosData_b % T = tempguess
  eosData_b % e = emq + 9.6485e17*qbar_b
  eosData_b % abar = 1.e0/sumyi
  eosData_b % zbar = ye * eosData_b % abar
  call eos(eos_input_re, eosData_b)

end subroutine Flame_rhJumpReactive
