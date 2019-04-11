subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_error_module
  use network
  use prob_params_module, only: center
  use probdata_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ inlet_mach, do_stratified, do_isentropic

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !
  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call amrex_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  inlet_mach = 1.e-1_rt
  do_stratified = .false.
  do_isentropic = .false.

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! set center variable in prob_params_module
  center(1) = 0.5_rt*(problo(1)+probhi(1))
  center(2) = 0.5_rt*(problo(2)+probhi(2))

end subroutine amrex_probinit


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use network, only: nspec
  use probdata_module
  use prob_params_module, only: center
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, &
                                 UFS, const_grav

  use amrex_constants_module, only : ZERO, HALF, ONE
  use eos_module
  use eos_type_module
  use inlet_bc_module
  use actual_eos_module, only : gamma_const
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), time, delta(2)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)

  real(rt) :: x,y,z
  real(rt) :: xn(nspec)
  real(rt) :: H,pres_base,dens_base
  integer  :: i,j,npts_1d
  real(rt), allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state

  ! state(state_l1:state_h1,state_l2:state_h2,1:NVAR) = ZERO

  if (do_stratified) then

      ! first make a 1D initial model for the entire domain
      npts_1d = (2.e0_rt*center(2)+1.e-8_rt) / delta(2)

      allocate(pressure(0:npts_1d-1))
      allocate(density (0:npts_1d-1))
      allocate(temp    (0:npts_1d-1))
      allocate(eint    (0:npts_1d-1))

      pres_base = 1.e6_rt
      dens_base = 1.e-3_rt

      ! only initialize the first species
      xn(:) = ZERO
      xn(1) = ONE

      ! compute the pressure scale height (for an isothermal, ideal-gas
      ! atmosphere)
      H = pres_base / dens_base / abs(const_grav)

      ! const = pres_base/dens_base**gamma_const

      pressure(0) = pres_base
      density(0)  = dens_base

      do j=0,npts_1d-1

         ! initial guess
         temp(j) = 1000.e0_rt

         if (do_isentropic) then
            z = dble(j) * delta(2)
            density(j) = dens_base*(const_grav*dens_base*(gamma_const - ONE)*z/ &
                 (gamma_const*pres_base) + ONE)**(ONE/(gamma_const - ONE))
         else
            z = (dble(j)+HALF) * delta(2)
            density(j) = dens_base * exp(-z/H)
         end if

         if (j .gt. 0) then
            pressure(j) = pressure(j-1) - &
                 delta(2) * HALF * (density(j)+density(j-1)) * abs(const_grav)
         end if

         eos_state%p = pressure(j)
         eos_state%T = temp(j)
         eos_state%rho = density(j)
         eos_state%xn(:) = xn(:)

         call eos(eos_input_rp, eos_state)

         eint(j) = eos_state%e
         temp(j) = eos_state%T

      end do

  endif

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        if (do_stratified) then

            state(i,j,URHO   ) = density(j)
            state(i,j,UMX:UMY) = ZERO
            state(i,j,UEDEN  ) = density(j) * eint(j)
            state(i,j,UEINT  ) = density(j) * eint(j)
            state(i,j,UTEMP  ) = temp(j)

        else

            ! use the EOS to make the state consistent
            eos_state%T    = 10.d0
            eos_state%rho   = 1.d-3
            eos_state%p     = 1.d6
            eos_state%xn(:) = 1.d0

            ! (rho,p) --> T, h
            call eos(eos_input_rp, eos_state)

            state(i,j,URHO   ) = eos_state%rho
            state(i,j,UMX:UMY) = ZERO
            state(i,j,UEDEN  ) = eos_state%rho * eos_state%e
            state(i,j,UEINT  ) = eos_state%rho * eos_state%e
            state(i,j,UTEMP  ) = eos_state%T

        endif

        state(i,j,UFS:UFS-1+nspec) = ZERO
        state(i,j,UFS  ) = state(i,j,URHO)

     enddo
  enddo

  call set_inlet_bcs()

end subroutine ca_initdata
