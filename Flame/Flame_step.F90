! see Flame_interface.F90 at top level for function description
!
! Dean Townsley 2008
!

! Implementation details
!
! This implementation is operator split betwee reaction and diffusion
! but dimensionally unsplit.  Diffusion is treated by directly computing
! the Laplacian instead of differencing fluxes, so that this is not a
! conservative diffusion operator.
!
! Computation of the Laplacian itself is done in a subroutine, as its
! form depends heavily on the mesh geometry

subroutine Flame_step( num_blocks, blockList, dt  )    

  use Flame_data

  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getBlkPtr, Grid_releaseBlkPtr, &
                             Grid_getBlkIndexLimits, Grid_fillGuardCells
  use fl_fsInterface, only : fl_flameSpeed
  use fl_effInterface, only: fl_effects
  use fl_interface, only : fl_laplacian
  use Driver_interface, only : Driver_abortFlash
  use Timers_interface, only : Timers_start, Timers_stop
  use Logfile_interface, only : Logfile_stamp
       
  implicit none
  integer, INTENT(in)                        :: num_blocks
  integer, INTENT(in), DIMENSION(num_blocks) :: blockList
  real,    INTENT(in)                        :: dt

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC, fspeedLimits
  integer :: istat
  integer :: n, bid

  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:,:,:) :: flam, flamdot, flamespeed, lapl

  real :: f, inv_dt
  integer :: i,j,k
  integer :: sizeI, sizeJ, sizeK

  inv_dt = 1.0/dt

     ! extract flam variable, should make cache work better
     ! need two layers in GCs becausee of RD splitting
     ! flam( blkLimits(LOW,IAXIS)-2 : blkLimits(HIGH,IAXIS)+2 , &
     !       blkLimits(LOW,JAXIS)-2*K2D : blkLimits(HIGH,JAXIS)+2*K2D , &
     !       blkLimits(LOW,KAXIS)-2*K3D : blkLimits(HIGH,KAXIS)+2*K3D) &
     !     = solnData(FLAM_MSCALAR, blkLimits(LOW,IAXIS)-2 : blkLimits(HIGH,IAXIS)+2 , &
     !             blkLimits(LOW,JAXIS)-2*K2D : blkLimits(HIGH,JAXIS)+2*K2D , &
     !             blkLimits(LOW,KAXIS)-2*K3D : blkLimits(HIGH,KAXIS)+2*K3D)


     call fl_flameSpeed(solnData, flamespeed, bid, 2)

     do k = blkLimits(LOW,KAXIS)-2*K3D, blkLimits(HIGH,KAXIS)+2*K3D
        do j = blkLimits(LOW,JAXIS)-2*K2D, blkLimits(HIGH,JAXIS)+2*K2D
           do i = blkLimits(LOW,IAXIS)-2, blkLimits(HIGH,IAXIS)+2
              f = flam(i,j,k)
              flam(i,j,k) = f + dt*fl_R_over_s*flamespeed(i,j,k)*(f-fl_epsilon_0)*(1.0+fl_epsilon_1-f)
           enddo
        enddo
     enddo

     ! 1 specifies the step size should be 1 grid cell
     ! cannot be any larger because flam is filled with only 2 guard cell layers
     call fl_laplacian(lapl, flam, 1, bid)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              flam(i,j,k) = max(0.0, min(1.0, flam(i,j,k) + dt*fl_kappa_over_s*flamespeed(i,j,k)*lapl(i,j,k) ) )
              flamdot(i,j,k) = (flam(i,j,k) - solnData(FLAM_MSCALAR,i,j,k))*inv_dt
              solnData(FLAM_MSCALAR, i,j,k) = flam(i,j,k)
           enddo
        enddo
     enddo

     deallocate(lapl)
     deallocate(flamespeed)
     deallocate(flam)

     call fl_effects( solnData, flamdot, dt, bid)

     deallocate(flamdot)

     call Grid_releaseBlkPtr(bid, solnData)

  enddo

  call Timers_stop("flame")

  return
end subroutine Flame_step
