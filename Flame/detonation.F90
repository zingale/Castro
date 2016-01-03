module detonation_module

  implicit none

  logical, save           :: initialized = .false.

  integer, parameter      :: infile_unit = 2

  ! Logical file unit for saving detonation ignition points
  integer, save           :: burn_lun = 31
  character(len=80), save :: detIgnFileName
  
  logical,          save :: thermalReact, autoDDT
  double precision, save :: thermalReactInFlameThreshold

  double precision, allocatable, save :: IgnTime(:), IgnLoc(:,:), IgnR(:)

  ! table for detonation ignition points
  double precision, save :: IgnRho, IgnRhoFact, IgnPhfa, IgnDist, IgnRad, IgnSep
  integer,          save :: IgnNum, IgnNumMax, IgnRhoCalib

  ! variable for saving processor-local and global neutrino loss energy integrals
  double precision, save :: neutLossThisProcStep, neutLoss

contains

  subroutine detonation_init() bind(C)

    use extern_probin_module, only: bn_detIgnFile, bn_detIgnFileName, bn_autoDDT, &
                                    bn_thermalReact, bn_thermalReactInFlameThreshold, bn_neutLoss, &
                                    bn_detIgnFile, bn_detIgnFileName, bn_autoDDT, &
                                    pbIgnDist, pbIgnNumMax, pbIgnPhfa, pbIgnRad, pbIgnRho, &
                                    pbIgnRhoCalib, pbIgnRhoFact, pbIgnSep

    implicit none

    integer :: i, istat

    logical :: detIgnFile

    if (initialized) return

!    if (.not. restart) then
       neutLoss = 0.0
!    endif
    neutLossThisProcStep = 0.0

    IgnDist = pbIgnDist
    IgnNumMax = pbIgnNumMax
    IgnPhfa = pbIgnPhfa
    IgnRad = pbIgnRad
    IgnRho = pbIgnRho
    IgnRhoCalib = pbIgnRhoCalib
    IgnRhoFact = pbIgnRhoFact
    IgnSep = pbIgnSep

    thermalReact = bn_thermalReact
    thermalReactInFlameThreshold = bn_thermalReactInFlameThreshold
    autoDDT = bn_autoDDT

    !  determine whether we are detonating manually or automatically
    detIgnFile     = bn_detIgnFile
    detIgnFileName = bn_detIgnFileName
    autoDDT        = bn_autoDDT

    if (detIgnFile .and. autoDDT) then
       call bl_error("Cannot detonate manually and automatically!")
    endif

    !  read detonation ignition points from file
    if (detIgnFile) then

       open(unit=infile_unit,file=detIgnFileName,status='OLD',iostat=istat)
       if (istat /= 0) call bl_error("Unable to open detonation ignition points file")

       ! eat header
       read(infile_unit,*)
       read(infile_unit,*) IgnNum
       allocate(IgnTime(IgnNum))
       allocate(IgnLoc(IgnNum,3))
       allocate(IgnR(IgnNum))
       do i = 1, IgnNum
          read(infile_unit,*) IgnTime(i), IgnLoc(i,1), IgnLoc(i,2), IgnLoc(i,3), IgnR(i)
       enddo
       close(unit=infile_unit)

    else if (autoDDT) then

       IgnRho      = pbIgnRho
       IgnRhoCalib = pbIgnRhoCalib
       IgnRhoFact  = pbIgnRhoFact
       IgnPhfa     = pbIgnPhfa
       IgnDist     = pbIgnDist
       IgnRad      = pbIgnRad
       IgnSep      = pbIgnSep
       IgnNumMax   = pbIgnNumMax

       IgnRho = exp( log(IgnRho) - IgnRhoCalib )
       if (IgnRhoFact >= 1.0)  &
            call bl_error("pbIgnRhoFact must be less than 1")
       IgnRhoFact = 10.e0**( IgnRhoFact * log10(IgnRho) ) 

       ! allocate global arrays for detonation points
       IgnNum = 0
       allocate(IgnTime(IgnNumMax),stat=istat)
       if (istat/=0) call bl_error("Cannot allocate IgnTime in detonation init")
       allocate(IgnLoc(IgnNumMax,3),stat=istat)
       if (istat/=0) call bl_error("Cannot allocate IgnLoc in detonation init")
       allocate(IgnR(IgnNumMax),stat=istat)
       if (istat/=0) call bl_error("Cannot allocate IgnR in detonation init")

       ! if (restart) then

       !    open(unit=infile_unit,file=detIgnFileName,status='unknown',IOSTAT=istat)
       !    if (istat /= 0) call bl_error("Unable to open detonation ignition points file")

       !    istat = 0
       !    do i = 1, IgnNumMax

       !       read(infile_unit,*,IOSTAT=istat) IgnTime(i), IgnLoc(i,1),  &
       !            IgnLoc(i,2), IgnLoc(i,3)
       !       if (istat/=0) exit

       !    enddo

       !    if (istat > 0) then
       !       call bl_error("Unable to read detonation ignition points file")
       !    else if (istat < 0) then !EOF reached
       !       IgnNum = i - 1
       !    else
       !       call bl_error("IgnNumMax is too small to read in previous detonation ignition points")
       !    endif

       !    close(unit=infile_unit)
       ! endif

    else

       IgnNum = 0

    endif

    initialized = .true.

  end subroutine detonation_init



  subroutine ca_get_ign_num_max(IgnNumMax_in) bind(C)

    implicit none

    integer, intent(inout) :: IgnNumMax_in

    IgnNumMax_in = IgnNumMax

  end subroutine ca_get_ign_num_max



  subroutine ca_get_ignsep(IgnSep_in) bind(C)

    implicit none

    integer, intent(inout) :: IgnSep_in

    IgnSep_in = IgnSep

  end subroutine ca_get_ignsep



  subroutine ca_check_ignition(lo, hi, state, s_lo, s_hi, ignition_coords) bind(C)

    use meth_params_module, only: NVAR, URHO, UEINT, UFS, UFX
    use bl_constants_module, only: ZERO, ONE
    use network, only: nspec
    use eos_module
    use castro_util_module, only: position
    
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: ignition_coords(3 * ignNumMax)

    integer :: ignition_num

    logical :: spark_box
    double precision :: spark_box_radius
    double precision :: rhoInv, radius

    logical :: ignition_test
    double precision :: det_loc(3), loc(3)

    type (eos_t) :: eos_state

    integer :: i, j, k

    if (autoDDT .and. thermalReact) then

       ignition_num = 0

       spark_box = .false.
       spark_box_radius = ZERO

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k)

                rhoInv = ONE / state(i,j,k,URHO)

                ! Get the pressure from the EOS

                eos_state % rho = state(i,j,k,URHO)
                eos_state % e   = state(i,j,k,UEINT) * rhoInv
                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
                eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

                call eos(eos_input_re, eos_state)

                call paraSpark(loc, state(i,j,k,URHO), eos_state % p,  &
                     state(i,j,k,UPHFA), state(i,j,k,UFLAM), &
                     state(i,j,k,UCI), state(i,j,k,UNEI),    &
                     det_loc, ignition_test )

                if ( ignition_test ) then

                   radius = sqrt( sum(det_loc**2) )

                   ! check to see if we have already ignited in this box
                   if (spark_box) then

                      ! check which is further out from the center
                      if (radius > spark_box_radius) then

                         ! if current radius is larger, then overwrite previous
                         ! ignition point
                         spark_box_radius = radius
                         ignition_coords((ignition_num-1) * 3 + 1:(ignition_num-1) * 3 + 3) = det_loc

                      endif !! (radius > spark_box_radius)

                      ! otherwise, our previous ignition point is further out
                      ! and we will keep it

                   else !! (spark_box)

                      ! we have not already found an ignition condition on this
                      ! box so, lets set a new one

                      ! make sure we have enough space
                      if ( ignition_num < IgnNumMax ) then

                         spark_box = .true.
                         ignition_coords( ignition_num * 3 + 1:ignition_num * 3 + 3 ) = det_loc

                         ignition_num = ignition_num + 1

                         ! we have too many ignition points, 
                         ! need to create more space
                      else !! ( ignition_num < pbIgnNumMax )

                         call bl_error("Not enough space to store all ignition points")

                      endif !! ( ignition_num < pbIgnNumMax )

                   endif !! (spark_box)

                endif !! (ignition_test)

             enddo
          enddo
       enddo

    end if

  end subroutine ca_check_ignition



  subroutine ca_check_valid_ignition(lo, hi, state, s_lo, s_hi, det_num, ignition_conditions, &
                                     det_xCoord, det_yCoord, det_zCoord) bind(C)

    use meth_params_module, only: NVAR, UFX
    use network, only: UPHFA
    use castro_util_module, only: position
    
    implicit none

    integer :: det_num
    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: det_xCoord(det_num), det_yCoord(det_num), det_zCoord(det_num)
    integer :: ignition_conditions(det_num)

    integer :: i, j, k, l
    double precision :: dist, loc(3)

    ! First, assume all points are valid

    ignition_conditions(:) = 1

    det_search: do l = 1, det_num
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k)

                ! Check if we are near a detonation point
                
                dist = sqrt( (det_xCoord(l) - loc(1))**2 + (det_yCoord(l) - loc(2))**2 + (det_zCoord(l) - loc(3))**2 )
                
                if (dist <= IgnRad) then
                   
                   ! Check if we are in ash, and we are at ignition density
                   if ( state(i,j,k,UFX+UPHFA-1) > IgnPhfa ) then
                      
                      ! we are in ash 
                      ! so invalidate detonation point
                      ignition_conditions(l) = 0

                      ! move to the next detonation point
                      cycle det_search
                   endif
                endif
             
             enddo
          enddo
       enddo
    enddo det_search

  end subroutine ca_check_valid_ignition



  subroutine ca_set_ignition_points(time, det_num, ignition_conditions, det_xCoord, det_yCoord, det_zCoord) bind(C)

    use bl_error_module

    implicit none

    integer :: det_num
    double precision :: det_xCoord(det_num), det_yCoord(det_num), det_zCoord(det_num)
    integer :: ignition_conditions(det_num)
    double precision :: time

    integer :: l
    logical :: ignition_conditions_logical(det_num)
    double precision :: detX, detY, detZ
    integer :: istat

    if (det_num > 0) then

       do l = 1, det_num
          if (ignition_conditions(l) == 1) then
             ignition_conditions_logical(l) = .true.
          else
             ignition_conditions_logical(l) = .false.
          endif
       enddo

       ! Announce detonation points and save to global list
       if (IgnNum + count(ignition_conditions_logical) < IgnNumMax) then
          do l = 1, det_num
             if ( ignition_conditions_logical(l) ) then

                ! if (ioproc) then
                !    detX = det_xCoord(l)
                !    detY = det_yCoord(l)
                !    detZ = det_zCoord(l)

                !    open (burn_lun, file=bn_detIgnFileName, &
                !         position='append', iostat=istat)
                !    if (istat/=0) call Driver_abortFlash("Unable to open detonation ignition points file")

                !    write (burn_lun,*) time, detX, detY, detZ

                !    close(unit=burn_lun)

                ! endif

                IgnNum = IgnNum + 1
                IgnTime( IgnNum ) = time
                IgnLoc( IgnNum, 1 ) = det_xCoord(l)
                IgnLoc( IgnNum, 2 ) = det_yCoord(l)
                IgnLoc( IgnNum, 3 ) = det_zCoord(l)

             endif
          enddo
       else
          call bl_error("Not enough space to save detonation points")
       endif

    endif

  end subroutine ca_set_ignition_points



  ! Aaron Jackson 2009
  !
  ! This subroutine checks for local ignition conditions.
  ! If found, we calculate the desired detonation point from these
  ! ignition conditions.  We compare the new detonation point with previous
  ! detonation points and make sure we are not detonating too close to a
  ! previous one.
  !
  ! If ignition_test is .true. then, det_loc contains the detonation
  ! coords.
  ! If ignition_test is .false. then they do not matter.

  subroutine paraSpark(loc, dens, pres, phfa, flame, c12, ne22, &
                       det_loc, ignition_test)

    use bl_constants_module, only : ZERO

    implicit none

    double precision, intent(in)  :: loc(3), dens, pres, phfa, flame, c12, ne22
    double precision, intent(out) :: det_loc(3)
    logical, intent(out)          :: ignition_test

    double precision, parameter :: ye12 = 0.5
    double precision, parameter :: ye16 = 0.5
    double precision, parameter :: ye22 = 10.0/22.0

    integer :: i
    double precision :: r, costheta, sintheta, detD, test_dens, yei

    ignition_test = .false.

    det_loc = ZERO

    ! only if we are at the flame front edge 
    if ( flame > 0.001 .and. flame < 0.01 ) then

       ! **********************************************************************
       !    This density is an estimate of the unburned density from
       ! Hansen & Kawaler.  This estimate uses the material composition and the
       ! pressure. This correspondence is good for dens >~ 1.5e5 (degenerate)
       ! and 0.92 is adjusted for fit.
       !    Since we are in the flame, and our ddt-density is in the
       ! fuel, this estimate produces more consistant results with placing
       ! the detonation points by hand.
       ! **********************************************************************
       !     yei = c12*ye12 + ne22*ye22 + (1.0-c12-ne22)*ye16
       !     yei = 1.0e0/yei

       !     test_dens = 0.92*yei*sqrt( (pres/1.243e15)**(6.0/4.0) + &
       !                                (pres/1.004e13)**(6.0/5.0) )

       ! **********************************************************************
       ! This density we test for is the local density in the flame.
       ! However, as material is burned, it becomes less dense and we could meet
       ! these conditions sooner than we want.
       ! **********************************************************************
       test_dens = dens

       ! check that rhofact > dens >= rho
       ! rhofact needs to be less than 1
       if ( test_dens > IgnRhoFact .and.  &
            test_dens <= IgnRho ) then

          ignition_test = .true.

          ! calculate the detonation coordinates
          r = sqrt( sum(loc**2) )

          ! Note that this is written for 2D cylindrical coordinates
          ! and will need to be updated for Cartesian.

          sintheta = loc(1) / r
          costheta = loc(2) / r

          r = r + IgnDist

          det_loc(1) = r * sintheta
          det_loc(2) = r * costheta
          det_loc(3) = ZERO

          ! now check that we are not igniting near something we've
          ! already detonated
          if (IgnNum > 0) then

             do i = 1, IgnNum

                ! distance to previous det point i
                detD = sqrt( sum( (det_loc - IgnLoc(i,:))**2 ) )

                ! if we're too close, then we fail spark test
                if ( detD <= IgnSep ) then

                   ignition_test = .false.
                   det_loc = ZERO

                   ! this detonation is too close to a previous one
                   ! fail ignition_test and exit loop.
                   exit

                endif ! end separation check

             enddo

          endif ! end previous detonation number check

       endif ! end density check

    endif ! end in-flame check

  end subroutine paraSpark

end module detonation_module
