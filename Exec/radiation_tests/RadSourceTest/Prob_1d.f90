
      subroutine PROBINIT (init,name,namlen,problo,probhi)

      use probdata_module
      use network, only : network_init
      use eos_module, only: eos_given_ReX
      use network, only : nspec
      implicit none

      integer :: init, namlen
      integer :: name(namlen)
      double precision :: problo(1), probhi(1)

      integer :: untin,i

      double precision :: gamma_eos, p_eos, c_eos, dpdr, dpde
      double precision :: X_in(nspec), e_0

      namelist /fortin/ rho_0, rhoe_0, E_rad


      !-----------------------------------------------------------------------
      ! Read in any runtime parameters
      !-----------------------------------------------------------------------

      ! Build "probin" filename -- the name of file containing fortin 
      ! namelist.

      integer maxlen
      parameter (maxlen=256)
      character probin*(maxlen)

      call network_init()

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
      end if

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do
         
      ! set namelist defaults

      ! initialize the refinement criteria
      raderr = 1.d20
      radgrad = 1.d20
      max_raderr_lev = -1
      max_radgrad_lev = -1


      ! read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)


      !-----------------------------------------------------------------------
      ! Set local variable defaults
      !-----------------------------------------------------------------------
      
      ! get T_0 corresponding to rhoe_0 and rho_0 through the EOS
      e_0 = rhoe_0/rho_0
      X_in(:) = 0.0
      X_in(1) = 1.d0

      call eos_given_ReX(gamma_eos, p_eos, c_eos, T_0, dpdr, dpde, rho_0, e_0, X_in)


      ! the 'center' variables are the location of the middle of the 
      ! domain -- this is where we put the interface
      center(1) = 0.5d0*(problo(1)+probhi(1))

      ! domain extrema
      xmin = problo(1)
      xmax = probhi(1)

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
                            state,state_l1,state_h1,delta,xlo,xhi)

     use probdata_module
     use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UFS, UTEMP
     use network, only : nspec
     use eos_module

     implicit none

     integer level, nscal
     integer lo(1), hi(1)
     integer state_l1,state_h1
     double precision state(state_l1:state_h1,NVAR)
     double precision time, delta(1)
     double precision xlo(1), xhi(1)

     double precision xcen
     integer i

      do i = lo(1), hi(1)
         xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)
            
         state(i,URHO ) = rho_0
         state(i,UMX  ) = 0.d0
         state(i,UEDEN) = rhoe_0 
         state(i,UEINT) = rhoe_0

         ! set the composition to be all in the first species
         state(i,UFS:UFS-1+nspec) = 0.d0
         state(i,UFS) = state(i,URHO)

         state(i,UTEMP) = T_0

      enddo

      end subroutine ca_initdata

! ::: 
! ::: -----------------------------------------------------------
! :::
      subroutine ca_initrad(level,time,lo,hi,nrad, &
           rad_state,rad_state_l1,rad_state_h1, &
           delta,xlo,xhi)

        use probdata_module,   only: E_rad
        use rad_params_module, only: nugroup, dnugroup
        use rad_params_module, only: pi, clight, hplanck, kboltz, stefbol

        implicit none
        integer level, nrad
        integer lo(1), hi(1)
        integer rad_state_l1,rad_state_h1

        double precision xlo(1), xhi(1), time, delta(1)
        double precision rad_state(rad_state_l1:rad_state_h1,0:nrad-1)

        ! local variables
        double precision trad, mgfac, nu, radtmp
        integer i, n

        if (nrad == 1) then 

           do i = lo(1), hi(1)  
              rad_state(i,0) = E_rad
           enddo

        else

           trad = (E_rad * clight * 0.25d0 / stefbol) ** (0.25d0)

           mgfac = 8.d0 * pi * hplanck / clight**3

           do n = 0, nrad-1
              nu     = nugroup(n)
              radtmp = exp(hplanck * nu / (kboltz * trad)) - 1.d0
              radtmp = (mgfac * nu**3 / radtmp) * dnugroup(n)

              do i = lo(1), hi(1)
                 rad_state(i,n) = radtmp
              enddo

           enddo

        endif

      end subroutine ca_initrad

! ::: 
! ::: -----------------------------------------------------------
! :::

     subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                           domlo,domhi,delta,xlo,time,bc)

     use meth_params_module, only : NVAR
     implicit none
     include 'bc_types.fi'

     integer          :: adv_l1,adv_h1
     integer          :: bc(1,2,*)
     integer          :: domlo(1), domhi(1)
     double precision :: delta(1), xlo(1), time
     double precision :: adv(adv_l1:adv_h1,NVAR)

     integer n

     do n = 1,NVAR
        call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
                   domlo,domhi,delta,xlo,bc(1,1,n))
     enddo

     do n = 1, NVAR

!        XLO
         if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
            print *,'SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) '
            stop
         end if

!        XHI
         if ( bc(1,2,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
            print *,'SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) '
            stop
         end if

     end do

     end subroutine ca_hypfill

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_denfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc)

      implicit none
      include 'bc_types.fi'

      integer          :: adv_l1,adv_h1
      integer          :: bc(1,2,*)
      integer          :: domlo(1), domhi(1)
      double precision :: delta(1), xlo(1), time
      double precision :: adv(adv_l1:adv_h1)

!     Note: this function should not be needed, technically, but is provided
!     to filpatch because there are many times in the algorithm when just
!     the density is needed.  We try to rig up the filling so that the same
!     function is called here and in hypfill where all the states are filled.

      call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

!     XLO
      if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
         print *,'SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
         stop
      end if

!     XHI
      if ( bc(1,2,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
         print *,'SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
         stop
      end if

      end subroutine ca_denfill
