module bc_mask_module

  private
  public set_bc_mask

contains

  subroutine set_bc_mask(lo, hi, domlo, domhi, &
                         x_bcMask, x_bcMask_lo, x_bcMask_hi, &
                         y_bcMask, y_bcMask_lo, y_bcMask_hi, &
                         z_bcMask, z_bcMask_lo, z_bcMask_hi) &
                         bind(C, name="set_bc_mask")

    use prob_params_module, only : physbc_lo, physbc_hi, problo, probhi

    implicit none


    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    integer, intent(in) :: x_bcMask_lo(3), x_bcMask_hi(3)
    integer, intent(in) :: y_bcMask_lo(3), y_bcMask_hi(3)
    integer, intent(in) :: z_bcMask_lo(3), z_bcMask_hi(3)

    integer, intent(inout) :: x_bcMask(x_bcMask_lo(1):x_bcMask_hi(1),x_bcMask_lo(2):x_bcMask_hi(2),x_bcMask_lo(3):x_bcMask_hi(3))
    integer, intent(inout) :: y_bcMask(y_bcMask_lo(1):y_bcMask_hi(1),y_bcMask_lo(2):y_bcMask_hi(2),y_bcMask_lo(3):y_bcMask_hi(3))
    integer, intent(inout) :: z_bcMask(z_bcMask_lo(1):z_bcMask_hi(1),z_bcMask_lo(2):z_bcMask_hi(2),z_bcMask_lo(3):z_bcMask_hi(3))

    if (x_bcMask_lo(1) == domlo(1)) then
    ! Left x face
      x_bcMask(domlo(1),x_bcMask_lo(2):x_bcMask_hi(2),x_bcMask_lo(3):x_bcMask_hi(3)) = physbc_lo(1)
    endif

    if (x_bcMask_hi(1) == domhi(1)+1) then
    ! Right x face
      x_bcMask(domhi(1)+1,x_bcMask_lo(2):x_bcMask_hi(2),x_bcMask_lo(3):x_bcMask_hi(3)) = physbc_hi(1)
    end if

#if AMREX_SPACEDIM >= 2
    if (y_bcMask_lo(2) == domlo(2)) then
    ! Left y face
      y_bcMask(y_bcMask_lo(1):y_bcMask_hi(1),domlo(2),y_bcMask_lo(3):y_bcMask_hi(3)) = physbc_lo(2)
    end if

    if (y_bcMask_hi(2) == domhi(2)+1) then
      y_bcMask(y_bcMask_lo(1):y_bcMask_hi(1),domhi(2)+1,y_bcMask_lo(3):y_bcMask_hi(3)) = physbc_hi(2)
    end if
#endif

#if AMREX_SPACEDIM == 3
    if (z_bcMask_lo(3) == domlo(3)) then
      z_bcMask(z_bcMask_lo(1):z_bcMask_hi(1),z_bcMask_lo(2):z_bcMask_hi(2), domlo(3)) = physbc_lo(3)
    end if

    if (z_bcMask_hi(3) == domhi(3)+1) then
      z_bcMask(z_bcMask_lo(1):z_bcMask_hi(1),z_bcMask_lo(2):z_bcMask_hi(2), domhi(3)+1) = physbc_hi(3)
    end if
#endif

  end subroutine set_bc_mask

end module bc_mask_module
