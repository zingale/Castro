subroutine PROBINIT (init,name,namlen,problo,probhi)

  use bl_constants_module
  use probdata_module
  use model_parser_module
  use bl_error_module
  use prob_params_module, only: dim, center, coord_type

  implicit none
  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  integer untin,i,j,k,dir

  namelist /fortin/ model_name, dens_fluff, temp_fluff, &
                    xc12_fluff, xne22_fluff, rep_ne_frac

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer, parameter :: maxlen = 127
  character probin*(maxlen)
  character model*(maxlen)
  integer ipp, ierr, ipp1

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

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

  if (dim .eq. 1) then

     call bl_error("This problem does not makes sense in 1D.")
     
  else if (dim .eq. 2) then

     if (coord_type .ne. 1) then
        call bl_error("In 2D this problem only makes sense in cylindrical coordinates.")
     endif
     
     ! Axisymmetric
     
     center(1) = ZERO
     center(2) = HALF * (problo(2) + probhi(2))
     center(3) = ZERO

  else

     call bl_error("This problem is not yet configured in 3D.")

  endif

end subroutine PROBINIT

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
                       state,s_lo,s_hi, &
                       dx,xlo,xhi)

  use bl_constants_module
  use probdata_module
  use interpolate_module
  use eos_module
  use meth_params_module, only: NVAR, URHO, UMX, UMZ, UTEMP,&
       UEDEN, UEINT, UFS
  use network, only : nspec
  use model_parser_module
  use prob_params_module, only: center
  use castro_util_module, only: position
  use eos_type_module
  use eos_module

  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: s_lo(3), s_hi(3)
  double precision :: xlo(3), xhi(3), time, dx(3)
  double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

  double precision :: loc(3), vel(3), r
  double precision :: rho, T, xn(nspec)
  integer :: i, j, k, n

  type (eos_t) :: eos_state

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)   

           loc = position(i, j, k) - center

           r = sqrt(sum(loc**2))

           ! Interpolate from 1D model

           rho = interpolate(r,npts_model,model_r,model_state(:,idens_model))
           T   = interpolate(r,npts_model,model_r,model_state(:,itemp_model))
           do n = 1, nspec
              xn(n) = interpolate(r,npts_model,model_r,model_state(:,ispec_model-1+n))
           enddo

           ! Thermodynamics

           eos_state % rho = rho
           eos_state % T   = T
           eos_state % xn  = xn

           call eos(eos_input_rt, eos_state)

           ! Model is initially stationary

           vel = ZERO

           ! Save data to the state array

           state(i,j,k,URHO)    = rho
           state(i,j,k,UTEMP)   = T
           state(i,j,k,UMX:UMZ) = rho * vel        
           state(i,j,k,UEINT)   = rho * eos_state % e        
           state(i,j,k,UEDEN)   = state(i,j,k,UEINT) + HALF * rho * sum(vel**2)
           state(i,j,k,UFS:UFS+nspec-1) = rho * xn

        enddo
     enddo
  enddo

end subroutine ca_initdata
