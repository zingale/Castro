! These routines do the characteristic tracing under the parabolic
! profiles in each zone to the edge / half-time.

module trace_ppm_module

  use prob_params_module, only : dg
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public tracexy_ppm, tracez_ppm

contains

  subroutine tracexy_ppm(q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                         qxm, qxp, qym, qyp, qs_lo, qs_hi, &
#if (AMREX_SPACEDIM < 3)
                         dloga, dloga_lo, dloga_hi, &
#endif
                         ilo1, ilo2, ihi1, ihi2, domlo, domhi, &
                         dx, dt, kc, k3d)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QTEMP, QGAME, QC, QGAMC, QFX, QFS, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   npassive, qpass_map, ppm_temp_fix, &
                                   fix_mass_flux
    use amrex_constants_module, only : ZERO, HALF, ONE
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qs_lo(3),qs_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: I_lo(3), I_hi(3)
#if (AMREX_SPACEDIM < 3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: kc, k3d
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,NQ)
    real(rt), intent(in) :: Im(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,NQ)

    real(rt), intent(in) :: Ip_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,QVAR)
    real(rt), intent(in) :: Im_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,QVAR)

    real(rt), intent(in) :: Ip_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,1)
    real(rt), intent(in) :: Im_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3,1)

    real(rt), intent(inout) :: qxm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qxp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qym(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qyp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
#if (AMREX_SPACEDIM < 3)
    real(rt), intent(in) ::  dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: dt, dx(3)

    ! Local variables
    integer :: i, j
    integer :: n, ipassive

    type(eos_t) :: eos_state

    real(rt) :: hdt

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! gam_g : the gas Gamma_1
    ! game : gas gamma_e
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    real(rt) :: cc, csq, cgassq, Clag
    real(rt) :: rho, u, v, w, p, rhoe_g, h_g, temp
    real(rt) :: gam_g

    real(rt) :: drho, dptot, drhoe_g
    real(rt) :: de, dge
    real(rt) :: dup, dvp, dptotp
    real(rt) :: dum, dvm, dptotm

    real(rt) :: rho_ref, u_ref, v_ref, p_ref, rhoe_g_ref, h_g_ref

    real(rt) :: cc_ref, csq_ref, Clag_ref, gam_g_ref
    real(rt) :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, h_g_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g
    real(rt) :: sourcr, sourcp, source, courn, eta, dlogatmp
    real(rt) :: tau_s, e_s

    logical :: fix_mass_flux_lo, fix_mass_flux_hi

#ifndef AMREX_USE_CUDA    
    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_nd.f90 :: tracexy_ppm")
    end if
#endif

    hdt = HALF * dt


    fix_mass_flux_lo = (fix_mass_flux == 1) .and. (physbc_lo(1) == Outflow) &
         .and. (ilo1 == domlo(1))
    fix_mass_flux_hi = (fix_mass_flux == 1) .and. (physbc_hi(1) == Outflow) &
         .and. (ihi1 == domhi(1))


    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    !
    ! We come in with the Im and Ip arrays -- these are the averages
    ! of the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    !
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right ("p") state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left ("m") state at that interface).
    !
    ! The indices are: Ip(i, j, k, dim, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.


    !-------------------------------------------------------------------------
    ! x-direction
    !-------------------------------------------------------------------------

    ! Trace to left and right edges using upwind PPM

    do j = ilo2-dg(2), ihi2+dg(2)
       do i = ilo1-1, ihi1+1

          rho = q(i,j,k3d,QRHO)

          cc = qaux(i,j,k3d,QC)
          csq = cc**2
          Clag = rho*cc

          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)

          p = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          h_g = ( (p + rhoe_g)/rho)/csq
          temp = q(i,j,k3d,QTEMP)

          gam_g = qaux(i,j,k3d,QGAMC)


          !-------------------------------------------------------------------
          ! plus state on face i
          !-------------------------------------------------------------------

          if (i >= ilo1) then

             ! Set the reference state
             ! This will be the fastest moving state to the left --
             ! this is the method that Miller & Colella and Colella &
             ! Woodward use
             rho_ref  = Im(i,j,kc,1,1,QRHO)
             u_ref    = Im(i,j,kc,1,1,QU)

             p_ref    = Im(i,j,kc,1,1,QPRES)
             rhoe_g_ref = Im(i,j,kc,1,1,QREINT)

             gam_g_ref  = Im_gc(i,j,kc,1,1,1)

             rho_ref = max(rho_ref, small_dens)
             p_ref = max(p_ref, small_pres)

             ! For tracing (optionally)
             csq_ref = gam_g_ref*p_ref/rho_ref
             cc_ref = sqrt(csq_ref)
             Clag_ref = rho_ref*cc_ref
             h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c

             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)


             ! we also add the sources here so they participate in the tracing
             dum = u_ref - Im(i,j,kc,1,1,QU) - hdt*Im_src(i,j,kc,1,1,QU)
             dptotm = p_ref - Im(i,j,kc,1,1,QPRES) - hdt*Im_src(i,j,kc,1,1,QPRES)

             drho = rho_ref - Im(i,j,kc,1,2,QRHO) - hdt*Im_src(i,j,kc,1,2,QRHO)
             dptot = p_ref - Im(i,j,kc,1,2,QPRES) - hdt*Im_src(i,j,kc,1,2,QPRES)
             drhoe_g = rhoe_g_ref - Im(i,j,kc,1,2,QREINT) - hdt*Im_src(i,j,kc,1,2,QREINT)

             dup = u_ref - Im(i,j,kc,1,3,QU) - hdt*Im_src(i,j,kc,1,3,QU)
             dptotp = p_ref - Im(i,j,kc,1,3,QPRES) - hdt*Im_src(i,j,kc,1,3,QPRES)

             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             h_g_ev = h_g_ref
             p_ev    = p_ref

             ! we are generally working with p in our eigensystem

             ! (rho, u, p, (rho e) eigensystem

             ! These are analogous to the beta's from the original PPM
             ! paper (except we work with rho instead of tau).  This is
             ! simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dum)*(rho_ev/cc_ev)
             alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dup)*(rho_ev/cc_ev)
             alpha0r = drho - dptot/csq_ev
             alpha0e_g = drhoe_g - dptot*h_g_ev  ! note h_g has a 1/c**2 in it

             alpham = merge(ZERO, -alpham, u-cc > ZERO)
             alphap = merge(ZERO, -alphap, u+cc > ZERO)
             alpha0r = merge(ZERO, -alpha0r, u > ZERO)
             alpha0e_g = merge(ZERO, -alpha0e_g, u > ZERO)

             ! The final interface states are just
             ! q_s = q_ref - sum(l . dq) r
             ! note that the a{mpz}right as defined above have the minus already
             qxp(i,j,kc,QRHO  ) =  rho_ref +  alphap + alpham + alpha0r
             qxp(i,j,kc,QU    ) =    u_ref + (alphap - alpham)*cc_ev/rho_ev
             qxp(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
             qxp(i,j,kc,QPRES ) =    p_ref + (alphap + alpham)*csq_ev

             ! Enforce small_*
             qxp(i,j,kc,QRHO ) = max(qxp(i,j,kc,QRHO ), small_dens)
             qxp(i,j,kc,QPRES) = max(qxp(i,j,kc,QPRES), small_pres)

             ! Transverse velocities -- there's no projection here, so
             ! we don't need a reference state.  We only care about
             ! the state traced under the middle wave

             ! Recall that I already takes the limit of the parabola
             ! in the event that the wave is not moving toward the
             ! interface
             qxp(i,j,kc,QV) = Im(i,j,kc,1,2,QV) + hdt*Im_src(i,j,kc,1,2,QV)
             qxp(i,j,kc,QW) = Im(i,j,kc,1,2,QW) + hdt*Im_src(i,j,kc,1,2,QW)

          end if


          !-------------------------------------------------------------------
          ! minus state on face i + 1
          !-------------------------------------------------------------------
          if (i <= ihi1) then

             ! Set the reference state
             ! This will be the fastest moving state to the right
             rho_ref  = Ip(i,j,kc,1,3,QRHO)
             u_ref    = Ip(i,j,kc,1,3,QU)

             p_ref    = Ip(i,j,kc,1,3,QPRES)
             rhoe_g_ref = Ip(i,j,kc,1,3,QREINT)

             gam_g_ref  = Ip_gc(i,j,kc,1,3,1)

             rho_ref = max(rho_ref,small_dens)
             p_ref = max(p_ref, small_pres)

             ! For tracing (optionally)
             csq_ref = gam_g_ref*p_ref/rho_ref
             cc_ref = sqrt(csq_ref)
             Clag_ref = rho_ref*cc_ref
             h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c

             dum = u_ref - Ip(i,j,kc,1,1,QU) - hdt*Ip_src(i,j,kc,1,1,QU)
             dptotm  = p_ref - Ip(i,j,kc,1,1,QPRES) - hdt*Ip_src(i,j,kc,1,1,QPRES)

             drho = rho_ref - Ip(i,j,kc,1,2,QRHO) - hdt*Ip_src(i,j,kc,1,2,QRHO)
             dptot = p_ref - Ip(i,j,kc,1,2,QPRES) - hdt*Ip_src(i,j,kc,1,2,QPRES)
             drhoe_g = rhoe_g_ref - Ip(i,j,kc,1,2,QREINT) - hdt*Ip_src(i,j,kc,1,2,QREINT)

             dup = u_ref - Ip(i,j,kc,1,3,QU) - hdt*Ip_src(i,j,kc,1,3,QU)
             dptotp = p_ref - Ip(i,j,kc,1,3,QPRES) - hdt*Ip_src(i,j,kc,1,3,QPRES)

             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             h_g_ev = h_g_ref
             p_ev    = p_ref

             ! (rho, u, p, (rho e)) eigensystem

             ! These are analogous to the beta's from the original PPM
             ! paper (except we work with rho instead of tau).  This is
             ! simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dum)*(rho_ev/cc_ev)
             alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dup)*(rho_ev/cc_ev)
             alpha0r = drho - dptot/csq_ev
             alpha0e_g = drhoe_g - dptot*h_g_ev  ! h_g has a 1/c**2 in it

             alpham = merge(-alpham, ZERO, u-cc > ZERO)
             alphap = merge(-alphap, ZERO, u+cc > ZERO)
             alpha0r = merge(-alpha0r, ZERO, u > ZERO)
             alpha0e_g = merge(-alpha0e_g, ZERO, u > ZERO)

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}left as defined above have the minus already

             qxm(i+1,j,kc,QRHO  ) =  rho_ref +  alphap + alpham + alpha0r
             qxm(i+1,j,kc,QU    ) =    u_ref + (alphap - alpham)*cc_ev/rho_ev
             qxm(i+1,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
             qxm(i+1,j,kc,QPRES ) =    p_ref + (alphap + alpham)*csq_ev

             ! Enforce small_*
             qxm(i+1,j,kc,QRHO ) = max(qxm(i+1,j,kc,QRHO ),small_dens)
             qxm(i+1,j,kc,QPRES) = max(qxm(i+1,j,kc,QPRES),small_pres)

             ! transverse velocities
             qxm(i+1,j,kc,QV    ) = Ip(i,j,kc,1,2,QV) + hdt*Ip_src(i,j,kc,1,2,QV)
             qxm(i+1,j,kc,QW    ) = Ip(i,j,kc,1,2,QW) + hdt*Ip_src(i,j,kc,1,2,QW)

          end if

          !-------------------------------------------------------------------
          ! geometry source terms
          !-------------------------------------------------------------------

#if (AMREX_SPACEDIM < 3)
          if (dloga(i,j,k3d) /= 0) then
             courn = dt/dx(1)*(cc+abs(u))
             eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k3d)))
             dlogatmp = min(eta, ONE)*dloga(i,j,k3d)
             sourcr = -HALF*dt*rho*dlogatmp*u
             sourcp = sourcr*csq
             source = sourcp*h_g

             if (i <= ihi1) then
                qxm(i+1,j,kc,QRHO) = qxm(i+1,j,kc,QRHO) + sourcr
                qxm(i+1,j,kc,QRHO) = max(qxm(i+1,j,kc,QRHO), small_dens)
                qxm(i+1,j,kc,QPRES) = qxm(i+1,j,kc,QPRES) + sourcp
                qxm(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QREINT) + source
             end if

             if (i >= ilo1) then
                qxp(i,j,kc,QRHO) = qxp(i,j,kc,QRHO) + sourcr
                qxp(i,j,kc,QRHO) = max(qxp(i,j,kc,QRHO), small_dens)
                qxp(i,j,kc,QPRES) = qxp(i,j,kc,QPRES) + sourcp
                qxp(i,j,kc,QREINT) = qxp(i,j,kc,QREINT) + source
             end if

          endif
#endif

#if (AMREX_SPACEDIM == 1)
          ! Enforce constant mass flux rate if specified
          if (fix_mass_flux_lo) then
             qxm(ilo1,j,kc,QRHO  ) = q(domlo(1)-1,j,k3d,QRHO)
             qxm(ilo1,j,kc,QU    ) = q(domlo(1)-1,j,k3d,QU  )
             qxm(ilo1,j,kc,QPRES ) = q(domlo(1)-1,j,k3d,QPRES)
             qxm(ilo1,j,kc,QREINT) = q(domlo(1)-1,j,k3d,QREINT)
          end if

          ! Enforce constant mass flux rate if specified
          if (fix_mass_flux_hi) then
             qxp(ihi1+1,j,kc,QRHO  ) = q(domhi(1)+1,j,k3d,QRHO)
             qxp(ihi1+1,j,kc,QU    ) = q(domhi(1)+1,j,k3d,QU  )
             qxp(ihi1+1,j,kc,QPRES ) = q(domhi(1)+1,j,k3d,QPRES)
             qxp(ihi1+1,j,kc,QREINT) = q(domhi(1)+1,j,k3d,QREINT)
          end if

#endif
       end do
    end do


    !-------------------------------------------------------------------------
    ! passively advected quantities
    !-------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  But we already dealt with all 3 velocity
       ! components above, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       do j = ilo2-dg(2), ihi2+dg(2)

          ! Plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,QU)

             ! We have
             !
             ! q_l = q_ref - Proj{(q_ref - I)}
             !
             ! and Proj{} represents the characteristic projection.
             ! But for these, there is only 1-wave that matters, the u
             ! wave, so no projection is needed.  Since we are not
             ! projecting, the reference state doesn't matter

             qxp(i,j,kc,n) = merge(q(i,j,k3d,n), Im(i,j,kc,1,2,n), u > ZERO)
          enddo

          ! Minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)

             qxm(i+1,j,kc,n) = merge(Ip(i,j,kc,1,2,n), q(i,j,k3d,n), u > ZERO)
          enddo

#if (AMREX_SPACEDIM == 1)
          if (fix_mass_flux_hi) qxp(ihi1+1,j,kc,n) = q(ihi1+1,j,k3d,n)
          if (fix_mass_flux_lo) qxm(ilo1,j,kc,n) = q(ilo1-1,j,k3d,n)
#endif

       enddo
    enddo


#if (AMREX_SPACEDIM >= 2)
    !-------------------------------------------------------------------------
    ! y-direction
    !-------------------------------------------------------------------------

    ! Trace to bottom and top edges using upwind PPM

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho = q(i,j,k3d,QRHO)

          cc = qaux(i,j,k3d,QC)
          csq = cc**2
          Clag = rho*cc

          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)

          p = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          h_g = ( (p + rhoe_g)/rho )/csq
          temp = q(i,j,k3d,QTEMP)

          gam_g = qaux(i,j,k3d,QGAMC)

          !-------------------------------------------------------------------
          ! plus state on face j
          !-------------------------------------------------------------------

          if (j >= ilo2) then

             ! Set the reference state
             ! This will be the fastest moving state to the left
             rho_ref  = Im(i,j,kc,2,1,QRHO)
             v_ref    = Im(i,j,kc,2,1,QV)

             p_ref    = Im(i,j,kc,2,1,QPRES)
             rhoe_g_ref = Im(i,j,kc,2,1,QREINT)

             gam_g_ref  = Im_gc(i,j,kc,2,1,1)

             rho_ref = max(rho_ref, small_dens)
             p_ref = max(p_ref, small_pres)

             ! For tracing (optionally)
             csq_ref = gam_g_ref*p_ref/rho_ref
             cc_ref = sqrt(csq_ref)
             Clag_ref = rho_ref*cc_ref
             h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c

             dvm = v_ref - Im(i,j,kc,2,1,QV) - hdt*Im_src(i,j,kc,2,1,QV)
             dptotm = p_ref - Im(i,j,kc,2,1,QPRES) - hdt*Im_src(i,j,kc,2,1,QPRES)

             drho = rho_ref - Im(i,j,kc,2,2,QRHO) - hdt*Im_src(i,j,kc,2,2,QRHO)
             dptot = p_ref - Im(i,j,kc,2,2,QPRES) - hdt*Im_src(i,j,kc,2,2,QPRES)
             drhoe_g = rhoe_g_ref - Im(i,j,kc,2,2,QREINT) - hdt*Im_src(i,j,kc,2,2,QREINT)

             dvp = v_ref - Im(i,j,kc,2,3,QV) - hdt*Im_src(i,j,kc,2,3,QV)
             dptotp = p_ref - Im(i,j,kc,2,3,QPRES) - hdt*Im_src(i,j,kc,2,3,QPRES)

             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             h_g_ev = h_g_ref
             p_ev    = p_ref

             ! (rho, u, p, (rho e) eigensystem

             ! These are analogous to the beta's from the original PPM
             ! paper (except we work with rho instead of tau).  This
             ! is simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dvm)*(rho_ev/cc_ev)
             alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dvp)*(rho_ev/cc_ev)
             alpha0r = drho - dptot/csq_ev
             alpha0e_g = drhoe_g - dptot*h_g_ev

             alpham = merge(ZERO, -alpham, v-cc > ZERO)
             alphap = merge(ZERO, -alphap, v+cc > ZERO)
             alpha0r = merge(ZERO, -alpha0r, v > ZERO)
             alpha0e_g = merge(ZERO, -alpha0e_g, v > ZERO)

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}right as defined above have the minus already
             qyp(i,j,kc,QRHO  ) = rho_ref + alphap + alpham + alpha0r
             qyp(i,j,kc,QV    ) = v_ref + (alphap - alpham)*cc_ev/rho_ev
             qyp(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
             qyp(i,j,kc,QPRES ) = p_ref + (alphap + alpham)*csq_ev

             ! Enforce small_*
             qyp(i,j,kc,QRHO ) = max(qyp(i,j,kc,QRHO ), small_dens)
             qyp(i,j,kc,QPRES) = max(qyp(i,j,kc,QPRES), small_pres)

             ! transverse velocities
             qyp(i,j,kc,QU    ) = Im(i,j,kc,2,2,QU) + hdt*Im_src(i,j,kc,2,2,QU)
             qyp(i,j,kc,QW    ) = Im(i,j,kc,2,2,QW) + hdt*Im_src(i,j,kc,2,2,QW)

          end if

          !-------------------------------------------------------------------
          ! minus state on face j+1
          !-------------------------------------------------------------------

          if (j <= ihi2) then

             ! Set the reference state
             ! This will be the fastest moving state to the right
             rho_ref  = Ip(i,j,kc,2,3,QRHO)
             v_ref    = Ip(i,j,kc,2,3,QV)

             p_ref    = Ip(i,j,kc,2,3,QPRES)
             rhoe_g_ref = Ip(i,j,kc,2,3,QREINT)

             gam_g_ref  = Ip_gc(i,j,kc,2,3,1)

             rho_ref = max(rho_ref, small_dens)
             p_ref = max(p_ref, small_pres)

             ! For tracing (optionally)
             csq_ref = gam_g_ref*p_ref/rho_ref
             cc_ref = sqrt(csq_ref)
             Clag_ref = rho_ref*cc_ref
             h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c

             dvm = v_ref - Ip(i,j,kc,2,1,QV) - hdt*Ip_src(i,j,kc,2,1,QV)
             dptotm = p_ref - Ip(i,j,kc,2,1,QPRES) - hdt*Ip_src(i,j,kc,2,1,QPRES)

             drho = rho_ref - Ip(i,j,kc,2,2,QRHO) - hdt*Ip_src(i,j,kc,2,2,QRHO)
             dptot = p_ref - Ip(i,j,kc,2,2,QPRES) - hdt*Ip_src(i,j,kc,2,2,QPRES)
             drhoe_g = rhoe_g_ref - Ip(i,j,kc,2,2,QREINT) - hdt*Ip_src(i,j,kc,2,2,QREINT)

             dvp = v_ref - Ip(i,j,kc,2,3,QV) - hdt*Ip_src(i,j,kc,2,3,QV)
             dptotp = p_ref - Ip(i,j,kc,2,3,QPRES) - hdt*Ip_src(i,j,kc,2,3,QPRES)

             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             h_g_ev = h_g_ref
             p_ev    = p_ref

             ! (rho, u, p, (rho e) eigensystem

             ! These are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dvm)*(rho_ev/cc_ev)
             alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dvp)*(rho_ev/cc_ev)
             alpha0r = drho - dptot/csq_ev
             alpha0e_g = drhoe_g - dptot*h_g_ev

             alpham = merge(-alpham , ZERO, v-cc > ZERO)
             alphap = merge(-alphap, ZERO, v+cc > ZERO)
             alpha0r = merge(-alpha0r, ZERO, v > ZERO)
             alpha0e_g = merge(-alpha0e_g, ZERO, v > ZERO)

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}left as defined above has the minus already

             qym(i,j+1,kc,QRHO  ) = rho_ref + alphap + alpham + alpha0r
             qym(i,j+1,kc,QV    ) = v_ref + (alphap - alpham)*cc_ev/rho_ev
             qym(i,j+1,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
             qym(i,j+1,kc,QPRES ) = p_ref + (alphap + alpham)*csq_ev


             ! Enforce small_*
             qym(i,j+1,kc,QRHO ) = max(qym(i,j+1,kc,QRHO ), small_dens)
             qym(i,j+1,kc,QPRES) = max(qym(i,j+1,kc,QPRES), small_pres)

             ! transverse velocities
             qym(i,j+1,kc,QU    ) = Ip(i,j,kc,2,2,QU) + hdt*Ip_src(i,j,kc,2,2,QU)
             qym(i,j+1,kc,QW    ) = Ip(i,j,kc,2,2,QW) + hdt*Ip_src(i,j,kc,2,2,QW)
          end if

       end do
    end do

    !-------------------------------------------------------------------------
    ! Passively advected quantities
    !-------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  But we already dealt with all 3 velocity
       ! components above, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       ! Plus state on face j
       do j = ilo2, ihi2+1
          do i = ilo1-1, ihi1+1
             v = q(i,j,k3d,QV)
             qyp(i,j,kc,n) = merge(q(i,j,k3d,n), Im(i,j,kc,2,2,n), v > ZERO)
          enddo
       end do

       ! Minus state on face j+1
       do j = ilo2-1, ihi2
          do i = ilo1-1, ihi1+1
             v = q(i,j,k3d,QV)
             qym(i,j+1,kc,n) = merge(Ip(i,j,kc,2,2,n), q(i,j,k3d,n), v > ZERO)
          enddo

       enddo
    enddo

#endif

  end subroutine tracexy_ppm



  subroutine tracez_ppm(q, qd_lo, qd_hi, &
                        qaux, qa_lo, qa_hi, &
                        Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                        qzm, qzp, qs_lo, qs_hi, &
                        ilo1, ilo2, ihi1, ihi2, domlo, domhi, &
                        dt, km, kc, k3d)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQAUX, QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QPRES, QTEMP, QGAME, QC, QGAMC, QFS, QFX, &
                                   small_dens, small_pres, &
                                   ppm_type, &
                                   npassive, qpass_map
    use amrex_constants_module, only : ZERO, HALF, ONE
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qs_lo(3),qs_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: I_lo(3), I_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: km, kc, k3d
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: Ip(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,NQ)
    real(rt), intent(in) :: Im(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,NQ)

    real(rt), intent(in) :: Ip_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)
    real(rt), intent(in) :: Im_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)

    real(rt), intent(in) :: Ip_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,1)
    real(rt), intent(in) :: Im_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,1)

    real(rt), intent(inout) :: qzm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)
    real(rt), intent(inout) :: qzp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),NQ)

    real(rt), intent(in) :: dt

    !     Local variables
    integer :: i, j
    integer :: n, ipassive

    real(rt) :: hdt

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! gam_g : the gas Gamma_1
    ! game : gas gamma_e
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    real(rt) :: cc, csq, cgassq, Clag
    real(rt) :: rho, u, v, w, p, rhoe_g, h_g, temp
    real(rt) :: gam_g

    real(rt) :: drho, dptot, drhoe_g
    real(rt) :: de, dge
    real(rt) :: dwp, dptotp
    real(rt) :: dwm, dptotm

    real(rt) :: rho_ref, w_ref, p_ref, rhoe_g_ref, h_g_ref

    real(rt) :: cc_ref, csq_ref, Clag_ref, gam_g_ref
    real(rt) :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, h_g_ev

    real(rt) :: alpham, alphap, alpha0r, alpha0e_g

    real(rt) :: tau_s, e_s

    type(eos_t) :: eos_state

    hdt = HALF * dt

#ifndef AMREX_USE_CUDA    
    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call amrex_error("Error:: trace_ppm_3d.f90 :: tracez_ppm")
    end if
#endif

    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! Trace to left and right edges using upwind PPM
    !
    ! Note: in contrast to the above code for x and y, here the loop
    ! is over interfaces, not over cell-centers.


    !-------------------------------------------------------------------------
    ! construct qzp  -- plus state on face kc
    !-------------------------------------------------------------------------
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          rho  = q(i,j,k3d,QRHO)

          cc   = qaux(i,j,k3d,QC)
          csq  = cc**2
          Clag = rho*cc

          u    = q(i,j,k3d,QU)
          v    = q(i,j,k3d,QV)
          w    = q(i,j,k3d,QW)

          p    = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          h_g = ( (p+rhoe_g)/rho )/csq
          temp = q(i,j,k3d,QTEMP)

          gam_g = qaux(i,j,k3d,QGAMC)

          ! Set the reference state
          ! This will be the fastest moving state to the left
          rho_ref  = Im(i,j,kc,3,1,QRHO)
          w_ref    = Im(i,j,kc,3,1,QW)

          p_ref    = Im(i,j,kc,3,1,QPRES)
          rhoe_g_ref = Im(i,j,kc,3,1,QREINT)

          gam_g_ref  = Im_gc(i,j,kc,3,1,1)

          rho_ref = max(rho_ref, small_dens)
          p_ref = max(p_ref, small_pres)

          ! For tracing (optionally)
          csq_ref = gam_g_ref*p_ref/rho_ref
          cc_ref = sqrt(csq_ref)
          Clag_ref = rho_ref*cc_ref
          h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm = w_ref - Im(i,j,kc,3,1,QW) - hdt*Im_src(i,j,kc,3,1,QW)
          dptotm = p_ref - Im(i,j,kc,3,1,QPRES) - hdt*Im_src(i,j,kc,3,1,QPRES)

          drho = rho_ref - Im(i,j,kc,3,2,QRHO) - hdt*Im_src(i,j,kc,3,2,QRHO)
          dptot = p_ref - Im(i,j,kc,3,2,QPRES) - hdt*Im_src(i,j,kc,3,2,QPRES)
          drhoe_g = rhoe_g_ref - Im(i,j,kc,3,2,QREINT) - hdt*Im_src(i,j,kc,3,2,QREINT)

          dwp = w_ref - Im(i,j,kc,3,3,QW) - hdt*Im_src(i,j,kc,3,3,QW)
          dptotp = p_ref - Im(i,j,kc,3,3,QPRES) - hdt*Im_src(i,j,kc,3,3,QPRES)

          rho_ev  = rho_ref
          cc_ev   = cc_ref
          csq_ev  = csq_ref
          Clag_ev = Clag_ref
          h_g_ev = h_g_ref
          p_ev    = p_ref

          ! (rho, u, p, (rho e) eigensystem

          ! These are analogous to the beta's from the original PPM
          ! paper.  This is simply (l . dq), where dq = qref - I(q)
          alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dwm)*(rho_ev/cc_ev)
          alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dwp)*(rho_ev/cc_ev)
          alpha0r = drho - dptot/csq_ev
          alpha0e_g = drhoe_g - dptot*h_g_ev

          alpham = merge(ZERO, -alpham, w-cc > ZERO)
          alphap = merge(ZERO, -alphap, w+cc > ZERO)
          alpha0r = merge(ZERO, -alpha0r, w > ZERO)
          alpha0e_g = merge(ZERO, -alpha0e_g, w > ZERO)

          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          ! note that the a{mpz}right as defined above have the minus already

          qzp(i,j,kc,QRHO  ) = rho_ref + alphap + alpham + alpha0r
          qzp(i,j,kc,QW    ) = w_ref + (alphap - alpham)*cc_ev/rho_ev
          qzp(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
          qzp(i,j,kc,QPRES ) = p_ref + (alphap + alpham)*csq_ev

          ! Enforce small_*
          qzp(i,j,kc,QRHO ) = max(qzp(i,j,kc,QRHO ), small_dens)
          qzp(i,j,kc,QPRES) = max(qzp(i,j,kc,QPRES), small_pres)

          ! transverse velocities
          qzp(i,j,kc,QU    ) = Im(i,j,kc,3,2,QU) + hdt*Im_src(i,j,kc,3,2,QU)
          qzp(i,j,kc,QV    ) = Im(i,j,kc,3,2,QV) + hdt*Im_src(i,j,kc,3,2,QV)


          !-------------------------------------------------------------------
          ! This is all for qzm -- minus state on face kc
          !-------------------------------------------------------------------

          ! Note this is different from how we do 1D, 2D, and the
          ! x and y-faces in 3D, where the analogous thing would have
          ! been to find the minus state on face kc+1

          rho  = q(i,j,k3d-1,QRHO)

          cc   = qaux(i,j,k3d-1,QC)
          csq  = cc**2
          Clag = rho*cc

          u    = q(i,j,k3d-1,QU)
          v    = q(i,j,k3d-1,QV)
          w    = q(i,j,k3d-1,QW)

          p    = q(i,j,k3d-1,QPRES)
          rhoe_g = q(i,j,k3d-1,QREINT)
          h_g = ( (p + rhoe_g)/rho)/csq
          temp = q(i,j,k3d-1,QTEMP)

          gam_g = qaux(i,j,k3d-1,QGAMC)

          ! Set the reference state
          ! This will be the fastest moving state to the right
          rho_ref  = Ip(i,j,km,3,3,QRHO)
          w_ref    = Ip(i,j,km,3,3,QW)

          p_ref    = Ip(i,j,km,3,3,QPRES)
          rhoe_g_ref = Ip(i,j,km,3,3,QREINT)

          gam_g_ref  = Ip_gc(i,j,km,3,3,1)

          rho_ref = max(rho_ref,small_dens)
          p_ref = max(p_ref,small_pres)

          ! For tracing (optionally)
          csq_ref = gam_g_ref*p_ref/rho_ref
          cc_ref = sqrt(csq_ref)
          Clag_ref = rho_ref*cc_ref
          h_g_ref = ( (p_ref + rhoe_g_ref)/rho_ref)/csq_ref

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm = w_ref - Ip(i,j,km,3,1,QW) - hdt*Ip_src(i,j,km,3,1,QW)
          dptotm = p_ref - Ip(i,j,km,3,1,QPRES) - hdt*Ip_src(i,j,km,3,1,QPRES)

          drho = rho_ref - Ip(i,j,km,3,2,QRHO) - hdt*Ip_src(i,j,km,3,2,QRHO)
          dptot = p_ref - Ip(i,j,km,3,2,QPRES) - hdt*Ip_src(i,j,km,3,2,QPRES)
          drhoe_g = rhoe_g_ref - Ip(i,j,km,3,2,QREINT) - hdt*Ip_src(i,j,km,3,2,QREINT)

          dwp = w_ref - Ip(i,j,km,3,3,QW) - hdt*Ip_src(i,j,km,3,3,QW)
          dptotp = p_ref - Ip(i,j,km,3,3,QPRES) - hdt*Ip_src(i,j,km,3,3,QPRES)

          rho_ev  = rho_ref
          cc_ev   = cc_ref
          csq_ev  = csq_ref
          Clag_ev = Clag_ref
          h_g_ev = h_g_ref
          p_ev    = p_ref

          ! (rho, u, p, (rho e) eigensystem

          ! These are analogous to the beta's from the original PPM
          ! paper.  This is simply (l . dq), where dq = qref - I(q)
          alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dwm)*(rho_ev/cc_ev)
          alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dwp)*(rho_ev/cc_ev)
          alpha0r = drho - dptot/csq_ev
          alpha0e_g = drhoe_g - dptot*h_g_ev

          alpham = merge(-alpham, ZERO, w-cc > ZERO)
          alphap = merge(-alphap, ZERO, w+cc > ZERO)
          alpha0r = merge(-alpha0r, ZERO, w > ZERO)
          alpha0e_g = merge(-alpha0e_g, ZERO, w > ZERO)

          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          ! note that the a{mpz}left as defined above have the minus already

          qzm(i,j,kc,QRHO  ) = rho_ref + alphap + alpham + alpha0r
          qzm(i,j,kc,QW    ) = w_ref + (alphap - alpham)*cc_ev/rho_ev
          qzm(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ev*csq_ev + alpha0e_g
          qzm(i,j,kc,QPRES ) = p_ref + (alphap + alpham)*csq_ev

          ! Enforce small_*
          qzm(i,j,kc,QRHO ) = max(qzm(i,j,kc,QRHO ),small_dens)
          qzm(i,j,kc,QPRES) = max(qzm(i,j,kc,QPRES),small_pres)

          ! Transverse velocity
          qzm(i,j,kc,QU    ) = Ip(i,j,km,3,2,QU) + hdt*Ip_src(i,j,km,3,2,QU)
          qzm(i,j,kc,QV    ) = Ip(i,j,km,3,2,QV) + hdt*Ip_src(i,j,km,3,2,QV)

       end do
    end do

    !-------------------------------------------------------------------------
    ! passively advected quantities
    !-------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! For DIM < 3, the velocities are included in the passive
       ! quantities.  But we already dealt with all 3 velocity
       ! components above, so don't process them here.
       if (n == QU .or. n == QV .or. n == QW) cycle

       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! Plus state on face kc
             w = q(i,j,k3d,QW)
             qzp(i,j,kc,n) = merge(q(i,j,k3d,n), Im(i,j,kc,3,2,n), w > ZERO)

             ! Minus state on face k
             w = q(i,j,k3d-1,QW)
             qzm(i,j,kc,n) = merge(Ip(i,j,km,3,2,n), q(i,j,k3d-1,n), w > ZERO)

          enddo
       enddo
    enddo

  end subroutine tracez_ppm

end module trace_ppm_module
