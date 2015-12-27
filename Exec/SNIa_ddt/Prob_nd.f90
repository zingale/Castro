subroutine PROBINIT (init,name,namlen,problo,probhi)

  use bl_constants_module
  use probdata_module
  use model_parser_module
  use bl_error_module
  use flame_module
  use prob_params_module, only: dim, center, coord_type
  use random_iterate_module, only: sim_LCGRandomIterate
  
  implicit none
  
  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(1:dim), probhi(1:dim)

  integer :: untin

  namelist /fortin/ model_name, dens_fluff, temp_fluff, &
                    xc12_fluff, xne22_fluff, rep_ne_frac, &
                    sim_ignite, sim_ignX, sim_ignY, sim_ignZ, sim_ignR, &
                    sim_ignMpole, sim_ignMpoleA, sim_ignMpoleMinL, &
                    sim_ignMpoleMaxL, sim_ignMPoleSym, sim_ignMPoleSeed, &
                    sim_ignSin, sim_ignSinN, sim_ignSinA, &
                    sim_refFluffDensThresh, sim_refFluffMargin, sim_refFluffLevel, &
                    sim_refNogenEnucThresh, sim_refNogenFldtThresh, &
                    sim_refNogenMargin, sim_refNogenLevel

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer, parameter :: maxlen = 127
  character probin*(maxlen)

  integer          :: i, j, jstart, jend, k

  integer          :: try, l, m, sgn
  logical          :: accept
  double precision :: u, v, r
  double precision :: deltaplus, deltaminus, costheta, sintheta
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

  ! only need to get width of artificial flame once
  call Flame_getWidth(sim_laminarWidth)

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



subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,s_lo,s_hi, &
                       dx,xlo,xhi)

  use bl_constants_module
  use probdata_module
  use interpolate_module
  use eos_module
  use meth_params_module
  use network
  use model_parser_module
  use flame_module
  use prob_params_module, only: center, dim
  use castro_util_module, only: position

  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: s_lo(3), s_hi(3)
  double precision :: xlo(3), xhi(3), time, dx(3)
  double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

  double precision :: loc(3), vel(3), radius
  integer :: i, j, k, l, m, n, p, nmax

  type (eos_t) :: zone_state, unburned_state, nse_state

  double precision :: costheta, sintheta, theta, phi
  double precision :: yi, yi_f, yi_a, ye, ye_f, ye_a
  double precision :: flam, qbar_nse, dqbar_qn, dyi_qn, enuc, fact
  double precision :: P_l_m, P_lm1_m, ign_dist, flame_radius, fsurf_distance, hold
  double precision :: qbar_f, qbar_a

  double precision :: cgsMeVperAmu = 9.6485e17  

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)   

           loc = position(i, j, k) - center

           radius = sqrt(sum(loc**2))

           ! Interpolate from 1D model

           unburned_state % rho = interpolate(radius,npts_model,model_r,model_state(:,idens_model))
           unburned_state % T   = interpolate(radius,npts_model,model_r,model_state(:,itemp_model))
           do n = 1, nspec
              unburned_state % xn(n) = interpolate(radius,npts_model,model_r,model_state(:,ispec_model-1+n))
           enddo

           call paraFuelAshProperties(unburned_state % xn(iC12), unburned_state % xn(iNe22), ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)

           unburned_state % abar = ONE / yi_f
           unburned_state % zbar = ye_f / yi_f

           if (.not. sim_ignite) then
              
              ! no burned material, only unburned
              
              zone_state = unburned_state
              call eos(eos_input_rt, zone_state)
              
              flam     = 0.0              
              ye       = zone_state % zbar / zone_state % abar
              dyi_qn   = 0.0
              dqbar_qn = 0.0
              enuc     = 0.0
              
           else
              
              !-----------------------------------------------
              ! initialize, including a burned region
              !-----------------------------------------------

              ! find distance from flame surface (positive is in front of flame)

              ! default to a spherical region centered at specified coordinates
              ! distance from center of ignition region
              
              ign_dist = sqrt( (loc(1) - sim_ignX)**2 + (loc(2) - sim_ignY)**2 + (loc(3) - sim_ignZ)**2 )

              flame_radius = sim_ignR

              if (sim_ignMpole .and. sim_ignSin)  &
                 call bl_error("ca_initdata: multipole and sinusoidal ignition are exclusive")

              if (sim_ignMpole) then
                 if (dim == 2) then
                    ! assume 2d is cylindrical
                    costheta = (loc(2)-sim_ignY)/ign_dist
                    phi = 0.0
                 else if (dim == 3) then
                    ! assume 3d is cartesian
                    costheta = (loc(3)-sim_ignZ)/ign_dist
                    phi = atan2( loc(2)-sim_ignY, loc(1)-sim_ignX )
                 endif
                 if (sim_ignMPoleSym) then
                    nmax = 0
                 else
                    nmax = 2*sim_ignMPoleMaxL
                 endif
                 sintheta = sqrt(1.0-costheta*costheta)
                 ! Legendre polynomial
                 do n = 0, nmax, 2
                    m = n / 2
                    do l = m, sim_ignMPoleMaxL
                       if ( l == m ) then
                          P_l_m = 1.0
                          if ( m > 0 ) then
                             fact = 1.0
                             do p = 1, m
                                P_l_m = -P_l_m*fact*sintheta
                                fact = fact + 2.0
                             enddo
                          endif
                          P_lm1_m = 0.0
                       else
                          hold = P_l_m
                          P_l_m = (costheta*(2*l-1)*P_l_m - (l+m-1)*P_lm1_m) / &
                                  dble(l-m)
                          P_lm1_m = hold
                       endif
                          
                       if ( l >= sim_ignMPoleMinL ) then
                          flame_radius = flame_radius +  &
                                mp_A(l,n)*P_l_m*cos(m*phi+mp_delta(l,n))
                          ! coefficients for -m are stored in n-1 index.
                          if (m /= 0) flame_radius = flame_radius +  &
                                mp_A(l,n-1)*P_l_m*cos(-m*phi+mp_delta(l,n-1))
                       endif
                    enddo
                 enddo
              endif

              if (sim_ignSin) then
                 ! sinusoidal (in theta) border between burned and unburned region
                 ! sim_ignitionSinN = number of periods between theta = 0 and pi
                 if (dim == 2) then
                    theta = acos( (loc(2)-sim_ignY)/ign_dist)
                    flame_radius = flame_radius - sim_ignSinA*cos(2*theta*sim_ignSinN)
                 else if (dim == 3) then
                    theta = acos( (loc(3)-sim_ignZ)/ign_dist)
                    flame_radius = flame_radius - sim_ignSinA*cos(2*theta*sim_ignSinN)
                 endif
              endif


              fsurf_distance = ign_dist - flame_radius

              ! determine local state in this zone
              ! assume deltas are equal
              if ( (fsurf_distance-0.5*dx(1)) > 1.5*sim_laminarWidth ) then
                 
                 ! whole cell unburned material
                 zone_state = unburned_state
                 call eos(eos_input_rt, zone_state)
                 flam = 0.0
                 ye = zone_state % zbar / zone_state % abar
                 dyi_qn = 0.0
                 dqbar_qn = 0.0
                 enuc = 0.0
                 
              else if ( (fsurf_distance+0.5*dx(1)) < -1.5*sim_laminarWidth ) then
                 
                 ! fully burned to NSE
                 call Flame_rhJumpReactive(unburned_state, qbar_f, zone_state, dqbar_qn, eos_input_rt)
                 flam   = ONE
                 dyi_qn = ONE / zone_state % abar
                 ye     = dyi_qn * zone_state % zbar
                 enuc   = 0.0
                 
              else
                 
                 ! partially burned
                 ! at least one cell will fall here (necessary to get initial refinement right)
                 call Flame_getProfile(fsurf_distance, flam)
                 
                 ! calculate properties of NSE final state
                 call Flame_rhJumpReactive(unburned_state, qbar_f, nse_state, qbar_nse, eos_input_rt)

                 ! calculate properties for partially burned material
                 ! note, in fact ye_f and ye_a should be equal
                 
                 yi = yi_f * (ONE-flam) + (ONE / nse_state % abar) * flam
                 
                 ye = ye_f * (ONE-flam) + (nse_state % zbar / nse_state % abar) * flam
                 
                 zone_state = unburned_state
                 zone_state % abar = ONE / yi
                 zone_state % zbar = ye / yi

                 ! put this in pressure equilibrium with unburned material
                 call Flame_rhJump(unburned_state, zone_state, flam*(qbar_nse-qbar_f)*cgsMeVperAmu, ZERO, eos_input_rt)

                 dyi_qn   = flam * ONE / nse_state % abar
                 dqbar_qn = flam * qbar_nse

                 ! to trigger refinement
                 enuc = 1.1*sim_refNogenEnucThresh
                 
              endif

           endif ! sim_ignite
           
           ! Model is initially stationary

           vel = ZERO
           
           ! Save data to the state array

           state(i,j,k,URHO)    = zone_state % rho
           state(i,j,k,UTEMP)   = zone_state % T
           state(i,j,k,UMX:UMZ) = zone_state % rho * vel        
           state(i,j,k,UEINT)   = zone_state % rho * zone_state % e
           state(i,j,k,UEDEN)   = state(i,j,k,UEINT) + HALF * zone_state % rho * sum(vel**2)
           state(i,j,k,UFS:UFS+nspec-1) = zone_state % rho * zone_state % xn

           state(i,j,k,UFX+UFLAM-1) = flam
           state(i,j,k,UFX+UFLDT-1) = ZERO
           state(i,j,k,UFX+UFLSP-1) = ZERO
           state(i,j,k,UFX+UCI  -1) = zone_state % xn(iC12)
           state(i,j,k,UFX+UNEI -1) = zone_state % xn(iNe22)
           state(i,j,k,UFX+UPHFA-1) = flam
           state(i,j,k,UFX+UPHAQ-1) = flam
           state(i,j,k,UFX+UPHQN-1) = flam
           state(i,j,k,UFX+UYE  -1) = ye
           state(i,j,k,UFX+UDYQN-1) = dyi_qn
           state(i,j,k,UFX+UDQQN-1) = dqbar_qn

        enddo
     enddo
  enddo

end subroutine ca_initdata
