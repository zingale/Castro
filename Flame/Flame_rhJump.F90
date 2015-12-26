
! Computes burned state from a Rankine-Hugoniot jump condition across the flame front
! 
! computes the thermodynamic properties of the ash (dens,temp) given
! those of the fuel, the energy release q, the front speed s (w.r.t. the fuel), and the
! composition of the ash material
!
! notes:  the unburned state is detemined by calling the eos in mode "mode"
!           using the supplied eosData_u and mfrac_u
!         This routine works with or without USE_EOS_LITE.  But the "compositon"
!           of the burned state is set in different ways:
!           with USE_EOS_LITE: set eosData_b(EOS_ABAR) and eosData_b(EOS_ZBAR)
!           without:           set mfrac_b
subroutine Flame_rhjump(eosData_u, eosData_b, q, s, mode)

  use flame_module, ONLY : fl_effsmlrho
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

1 format(5(2X, E10.4))

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

     if (dens_b .lt. fl_effsmlrho) then

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
