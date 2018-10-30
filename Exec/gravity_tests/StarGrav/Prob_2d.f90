subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use model_parser_module
  use amrex_error_module
  use prob_params_module, only : center

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer init, namlen
  integer name(namlen)
  real(rt)         problo(2), probhi(2)

  integer untin,i,j,k,dir

  namelist /fortin/ &
       model_name

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !
  integer, parameter :: maxlen = 127
  character probin*(maxlen)
  character model*(maxlen)
  integer ipp, ierr, ipp1

  if (namlen .gt. maxlen) call amrex_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)


  ! read initial model
  call read_model_file(model_name)

  ! assume axisymmetric
  center(1) = 0.e0_rt
  center(2) = 0.5e0_rt*(problo(2)+probhi(2))

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

  use probdata_module
  use interpolate_module
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UTEMP,&
       UEDEN, UEINT, UFS, small_temp
  use network, only : nspec
  use model_parser_module
  use prob_params_module, only : center
  use eos_type_module
  use eos_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), time, delta(2)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)

  real(rt)         xcen,ycen,dist,pres,total
  integer i,j,n
  real(rt)      :: model_eden(npts_model)
  real(rt) :: smallx = 1.d-10
  logical :: status

  type(eos_t) :: eos_state

  ! first calculate the eden across the model
  do i = 1, npts_model

      eos_state%rho = model_state(i,idens_model)
      eos_state%T = model_state(i,itemp_model)
      eos_state%xn(:) = model_state(i,ispec_model:ispec_model-1+nspec)

      call eos(eos_input_rt, eos_state)

      model_eden(i) = eos_state%e * model_state(i,idens_model)
  enddo

  do j = lo(2), hi(2)
     ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5e0_rt) - center(2)

     do i = lo(1), hi(1)
        xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5e0_rt) - center(1)

        dist = sqrt(xcen**2 + ycen**2)

        call conservative_interpolate(state(i,j,UEDEN), dist,npts_model,model_r,model_eden, delta, status)

        if (state(i,j,UEDEN) < 0.0d0 .or. state(i,j,UEDEN) /= state(i,j,UEDEN) .or. (.not. status)) then
            print *, "conservative interpolate of eden_model failed :("
            state(i,j,UEDEN) = interpolate(dist,npts_model,model_r,model_eden)
        endif

        do n = 1, nspec
           call conservative_interpolate(state(i,j,UFS-1+n),dist,npts_model,model_r,model_state(:,ispec_model-1+n)*model_state(:,idens_model), delta, status)

           if (state(i,j,UFS-1+n) < 0.0d0 .or. state(i,j,UFS-1+n) /= state(i,j,UFS-1+n) .or. (.not. status)) then
               print *, "conservative interpolate of X_i*dens_model failed :("
               state(i,j,UFS-1+n) = interpolate(dist,npts_model,model_r,model_state(:,ispec_model-1+n)*model_state(:,idens_model))
           endif
        enddo

        state(i,j,URHO) = sum(state(i,j,UFS:UFS+nspec-1))

        if (state(i,j,URHO) < 0.0d0 .or. state(i,j,URHO) /= state(i,j,URHO)) then
            print *, "summing of species' partial densities failed"
            print *, "species : ", state(i,j,UFS:UFS+nspec-1)
            state(i,j,URHO) = interpolate(dist,npts_model,model_r,model_state(:,idens_model))
        endif

     enddo
  enddo

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        eos_state%rho = state(i,j,URHO)
        eos_state%e = state(i,j,UEDEN) / state(i,j,URHO)
        eos_state%xn(:) = state(i,j,UFS:UFS-1+nspec) / state(i,j,URHO)
        eos_state % T   = 10000.0e0_rt

        call eos(eos_input_re, eos_state)

        state(i,j,UTEMP) = eos_state%T

        state(i,j,UEDEN) = eos_state%e

        state(i,j,UEINT) = state(i,j,URHO) * state(i,j,UEDEN)
        state(i,j,UEDEN) = state(i,j,URHO) * state(i,j,UEDEN)

        ! make sure that the species (mass fractions) sum to 1
        total = 0.0d0
        do n = 1, nspec
           total = total + state(i,j,UFS+n-1)
        enddo

        do n = 1,nspec
           state(i,j,UFS+n-1) = max(smallx, state(i,j,UFS+n-1)/total * state(i,j,URHO))
        enddo

     enddo
  enddo

  ! Initial velocities = 0
  state(state_l1:state_h1,state_l2:state_h2,UMX:UMY) = 0.e0_rt

end subroutine ca_initdata
