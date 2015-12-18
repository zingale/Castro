subroutine PROBINIT (init,name,namlen,problo,probhi)

  use bl_constants_module
  use probdata_module
  use model_parser_module
  use bl_error_module
  use prob_params_module, only: dim, center, coord_type
  use random_iterate_module, only: sim_LCGRandomIterate
  
  implicit none
  
  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(1:dim), probhi(1:dim)

  integer :: untin

  namelist /fortin/ model_name, dens_fluff, temp_fluff, &
                    xc12_fluff, xne22_fluff, rep_ne_frac, &
                    ignite, x_match, y_match, z_match, r_match, &
                    sim_ignMpole, sim_ignMpoleA, sim_ignMpoleMinL, &
                    sim_ignMpoleMaxL, sim_ignMPoleSym, sim_ignMPoleSeed, &
                    sim_ignSin, sim_ignSinN, sim_ignSinA, &
                    sim_refFluffDensThresh, sim_refFluffMargin, sim_refFluffLevel, &
                    sim_refNogenEnucThresh, sim_refNogenFldtThresh, &
                    sim_refNogenMargin, sim_refNogenLevel, &
                    useBurn, bn_thermalReact

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer, parameter :: maxlen = 127
  character probin*(maxlen)
  character model*(maxlen)

  integer          :: i, j, jstart, jend, k
  double precision :: lastradius, radius

  integer          :: try, l, m, sgn
  logical          :: accept
  double precision :: u, v, r, ir
  double precision :: deltaplus, deltaminus, costheta, sintheta, phi
  double precision :: P_lm1_m, P_l_m, hold
  double precision :: factlmm, factlpm, fact  
  
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



  ! Set up multipole ignition
  
  if (sim_ignMPole) then
     
     ! set up coefficients of multipole contributions
     ! check some things first to catch mistakes
     
     if (dim == 1) call bl_error("multipole ignition not sensible in 1-dimension")

     if (dim == 2) then
        
        if (coord_type .ne. 1) &
             call bl_error("multipole ignition assumes cylindrical geometry in 2D")
        
        ! in 2D the ignition condition must be symmetric
        sim_ignMPoleSym = .true.
        
     endif

     if ( (dim == 3) .and. (coord_type .ne. 0) ) &
        call bl_error("multipole ignition assumes cartesian geometry in 3D")

     ! the j-index is for m values ordered 0,-1,1,-2,2...
     jstart = 0
     if (sim_ignMPoleSym) then
        jend = 0 ! only m=0 for symmetric
     else
        jend = 2*sim_ignMpoleMaxL ! only support max_m = max_l
     endif

     allocate(mp_A(sim_ignMpoleMinL:sim_ignMpoleMaxL,jstart:jend))
     allocate(mp_delta(sim_ignMpoleMinL:sim_ignMpoleMaxL,jstart:jend))

     try = 0
     accept = .false.
     do while ( .not. accept )

        try = try + 1

        ! generate set of random coefficients
        do j = jstart, jend, 2
           m = j / 2
           do i = sim_ignMpoleMinL, sim_ignMpoleMaxL

              ! generate a normally distributed number, r
              ! via the Box-Muller method
              
              call sim_LCGRandomIterate(sim_ignMPoleSeed)
              
              u = dble(sim_ignMPoleSeed)/2147483646.0
              
              ! avoid taking log of zero
              if (u==0) u = 1
              
              call sim_LCGRandomIterate(sim_ignMPoleSeed)
              
              v = dble(sim_ignMPoleSeed)/2147483646.0
              
              r = sqrt(-2*log(u))*cos(2*M_PI*v)
              
              if (.not. sim_ignMPoleSym) then

                 ! need a complex coefficient, uniformly distributed

                 call sim_LCGRandomIterate(sim_ignMPoleSeed)

                 u = dble(sim_ignMPoleSeed)/2147483646.0

                 if (u==0) u = 1

                 ! map into range -pi to pi
                 mp_delta(i,j) = 2*M_PI*u - M_PI
                 
              else
                 
                 mp_delta(i,j) = 0.0
                 
              endif

              
              ! now set amplitude of this component, this will be multiplied by
              ! the legendre polynomial of order i to construct the initial 
              ! flame surface
              ! normalization is chosen to match spherical harmonics as
              ! given in Jackson.
              factlmm = 1.0   !!  (l - m)!
              if ( i - m > 0 ) then
                 do k = 1, i - m
                    factlmm = factlmm*k
                 enddo
              endif

              factlpm = 1.0   !!  (l + m)!
              if ( i + m > 0 ) then
                 do k = 1, i + m
                    factlpm = factlpm*k
                 enddo
              endif

              fact = sqrt((2*i+1)/4.0/M_PI * factlmm/factlpm)

              mp_A(i,j) = r*sim_ignMpoleA * fact

              ! Need to account for 3D, there are 2l+1 terms as opposed
              ! to just one.
              if (.not. sim_ignMPoleSym) mp_A(i,j) = mp_A(i,j) / sqrt(dble(2*i+1))


              ! this will only happen for sim_ignMPoleSym = .false.
              if (m > 0) then
                 call sim_LCGRandomIterate(sim_ignMPoleSeed)
                 u = dble(sim_ignMPoleSeed)/2147483646.0
                 if (u==0) u = 1
                 call sim_LCGRandomIterate(sim_ignMPoleSeed)
                 v = dble(sim_ignMPoleSeed)/2147483646.0
                 r = sqrt(-2*log(u))*cos(2*M_PI*v)
                 ! need a complex coefficient
                 call sim_LCGRandomIterate(sim_ignMPoleSeed)
                 u = dble(sim_ignMPoleSeed)/2147483646.0
                 if (u==0) u = 1
                 ! this is for negative m which is equivalent to the
                 ! complex conjugate (c.c.) of positive m up to a sign.
                 ! the coefficient itself is not c.c.'ed but it will
                 ! be multiplied by the c.c.'ed part of Ylm.
                 ! Account for sign change here.
                 mp_delta(i,j-1) = M_PI - 2*M_PI*u

                 if (mod(m,2) == 0) then 
                    sgn = 1
                 else 
                    sgn = -1
                 endif

                 mp_A(i,j-1) = sgn * r*sim_ignMpoleA * fact

                 ! Need to account for 3D, there are 2l+1 terms as opposed
                 ! to just one.
                 mp_A(i,j-1) = mp_A(i,j-1) / sqrt(dble(2*i+1))
              endif

           ! another initialization which turned out to be non-stardard
           ! this was used for some early testing and is just kept for reference
           !call sim_LCGRandomIterate(sim_ignMPoleSeed)
           ! random phase angle, between 0 and pi since we will not use imaginary part
           !alpha = dble(sim_ignMPoleSeed)/2147483646.0 * M_PI
           ! coefficients have amplitudes of equal complex amplitude, with random
           ! phase.  We then use real part.  Normalization is for spherical
           ! harmonics in Jackson
           !mp_A(i) = sqrt((2*i+1)/4.0/M_PI)*cos(alpha)*sim_ignMpoleA*2/(sim_ignMpoleMaxL-sim_ignMpoleMinL+1)

           !print *, i,alpha/M_PI, mp_A(i)
           enddo
        enddo

        ! now check that this is negative on + and - axis to help prevent pathologies
        if (sim_ignMPoleSym) then

           deltaplus = 0.0
           costheta = 1.0
           sintheta = 0.0

           ! Legendre polynomial
           m = 0 ! only possible value for symmetric
           do l = m, sim_ignMPoleMaxL
              if ( l == m ) then
                 P_l_m = 1.0
                 P_lm1_m = 0.0
              else
                 hold = P_l_m
                 P_l_m = (costheta*(2*l-1)*P_l_m - (l+m-1)*P_lm1_m)/dble(l-m)
                 P_lm1_m = hold
              endif

              if ( l >= sim_ignMPoleMinL ) then
                    deltaplus = deltaplus + mp_A(l,m)*P_l_m
              endif
           enddo

           deltaminus = 0.0
           costheta = -1.0
           sintheta = 0.0
           ! Legendre polynomial
           do l = m, sim_ignMPoleMaxL
              if ( l == m ) then
                 P_l_m = 1.0
                 P_lm1_m = 0.0
              else
                 hold = P_l_m
                 P_l_m = (costheta*(2*l-1)*P_l_m - (l+m-1)*P_lm1_m)/dble(l-m)
                 P_lm1_m = hold
              endif

              if ( l >= sim_ignMPoleMinL ) then
                 deltaminus = deltaminus + mp_A(l,m)*P_l_m
              endif
           enddo

           if (deltaplus < 0.0 .and. deltaminus < 0.0) accept = .true.

        else ! not symmetric so accept as is

           accept = .true.

        endif

     enddo

  endif



  ! Read initial model
  
  call read_model_file(model_name)



  ! Make sure we're dimensionally safe.
  
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

     center = HALF * (problo + probhi)

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
