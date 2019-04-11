module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use prob_params_module, only: center
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ , UEDEN, UEINT, UFS, UTEMP, const_grav
    use interpolate_module
    use eos_module
    use eos_type_module
    use actual_eos_module, only : gamma_const
    use network, only: nspec
    use inlet_bc_module

    use amrex_fort_module     , only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer i,j,n
    real(rt)         x,y
    real(rt)         A, B, vx, vy
    real(rt)         X_in(nspec)

    type (eos_t) :: eos_state

    A = 4.5e-2_rt
    B = 1.e2_rt

    do n=1,NVAR
       call filcc(adv(adv_l1,adv_l2,n),adv_l1,adv_l2,adv_h1,adv_h2, &
            domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    !        XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then

       do j=adv_l2,adv_h2
          do i=domlo(1)-1,adv_l1,-1
             adv(i,j,:) = adv(domlo(1),j,:)
          end do
       end do

    end if

    !        XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then

       do j=adv_l2,adv_h2
          do i=domhi(1)+1,adv_h1
             adv(i,j,:) = adv(domhi(1),j,:)
          end do
       end do

    end if

    ! write(*,*) "prefactor", (inlet_mach/1.e-1_rt)* INLET_CS*1.e-2_rt, INLET_CS

    !        YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       ! this do loop counts backwards since we want to work downward
       do j=domlo(2)-1,adv_l2,-1
          do i=adv_l1,adv_h1

             x = xlo(1) + delta(1)*(dble(i-adv_l1) + HALF)
             adv(i,j,URHO) = INLET_RHO
             vx = ZERO
             vy = (inlet_mach/1.e-1_rt)* &
                  INLET_CS*(1.e-2_rt + A*(tanh(B*(x-0.40e0_rt)) + tanh(B*(0.6e0_rt-x))))
             adv(i,j,UMX) = vx
             adv(i,j,UMY) = INLET_RHO * vy
             adv(i,j,UTEMP) = INLET_TEMP

             adv(i,j,UEINT) = INLET_RHO * INLET_E
             adv(i,j,UEDEN) = INLET_RHO * INLET_E + HALF * INLET_RHO * (vx*vx + vy*vy)

             adv(i,j,UMZ) = ZERO
             adv(i,j,UFS:UFS-1+nspec) = ZERO
             adv(i,j,UFS) = adv(i,j,URHO)
          end do
       end do

       ! write(*,*) adv(:,domlo(2)-1,UMY)
    end if

    !        YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then

       do j=domhi(2)+1,adv_h2
          do i=adv_l1,adv_h1
             adv(i,j,: ) = adv(i,domhi(2),:)
             adv(i,j,UMX) = ZERO
             adv(i,j,UMY) = ZERO
             adv(i,j,UMZ) = ZERO
          end do
       end do

    end if

  end subroutine ca_hypfill

  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
       domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    use interpolate_module
    use actual_eos_module, only: gamma_const
    use meth_params_module, only : const_grav
    use inlet_bc_module

    use amrex_constants_module, only : M_PI, sixth
    use amrex_fort_module     , only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2)

    integer i,j

    !     Note: this function should not be needed, technically, but is provided
    !     to filpatch because there are many times in the algorithm when just
    !     the density is needed.  We try to rig up the filling so that the same
    !     function is called here and in hypfill where all the states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       do j=adv_l2,adv_h2
          do i=adv_l1,domlo(1)-1
             adv(i,j) = adv(domlo(1),j)
          end do
       end do
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       do j=adv_l2,adv_h2
          do i=domhi(1)+1,adv_h1
             adv(i,j) = adv(domhi(1),j)
          end do
       end do
    end if

    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       do j=domlo(2)-1,adv_l2,-1
          do i=adv_l1,adv_h1
             adv(i,j) = INLET_RHO
          end do
       end do
    end if

    !     YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
       do j=domhi(2)+1,adv_h2
          do i=adv_l1,adv_h1
             adv(i,j) = adv(i,domhi(2))
          end do
       end do
    end if

  end subroutine ca_denfill

end module bc_fill_module
