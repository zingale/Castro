module trace_ppm_module

  implicit none

  private

  public tracexy_ppm, tracez_ppm

contains

  subroutine tracexy_ppm(q,c,flatn,qd_lo,qd_hi, &
                         Ip,Im,Ip_src,Im_src,Ip_gc,Im_gc,I_lo,I_hi, &
                         qxm,qxp,qym,qyp,qs_lo,qs_hi, &
                         gamc,gc_lo,gc_hi, &
                         ilo1,ilo2,ihi1,ihi2,dt,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QPRES, QGAME, &
         small_dens, small_pres, &
         ppm_type, ppm_reference, ppm_trace_sources, &
         ppm_tau_in_tracing, ppm_reference_eigenvectors, &
         ppm_reference_edge_limit, &
         ppm_predict_gammae, &
         npassive, qpass_map
    use bl_constants_module

    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qs_lo(3),qs_hi(3)
    integer, intent(in) :: gc_lo(3), gc_hi(3)
    integer, intent(in) :: I_lo(3), I_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: kc, k3d

    double precision, intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    double precision, intent(in) ::     c(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision, intent(in) :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))

    double precision, intent(in) :: Ip(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)
    double precision, intent(in) :: Im(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)

    double precision, intent(in) :: Ip_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)
    double precision, intent(in) :: Im_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)

    double precision, intent(in) :: Ip_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,1)
    double precision, intent(in) :: Im_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,1)

    double precision, intent(inout) :: qxm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QVAR)
    double precision, intent(inout) :: qxp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QVAR)
    double precision, intent(inout) :: qym(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QVAR)
    double precision, intent(inout) :: qyp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QVAR)

    double precision, intent(in) :: gamc(gc_lo(1):gc_hi(1),gc_lo(2):gc_hi(2),gc_lo(3):gc_hi(3))

    double precision, intent(in) :: dt

    ! Local variables
    integer i, j
    integer n
    integer ipassive

    double precision cc, csq, Clag, rho, u, v, w, p, rhoe_g, enth

    double precision drho, du, dv, dw, dptot, drhoe_g, de, dge, dtau
    double precision dup, dvp, dptotp
    double precision dum, dvm, dptotm

    double precision :: rho_ref, u_ref, v_ref, p_ref, rhoe_g_ref, tau_ref
    double precision :: tau_s, e_s

    double precision :: cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref, game_ref, gfactor
    double precision :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev, tau_ev
    double precision :: gam, game

    double precision alpham, alphap, alpha0r, alpha0e

    double precision halfdt

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call bl_error("Error:: trace_ppm_3d.f90 :: tracexy_ppm")
    end if

    halfdt = HALF * dt


    !==========================================================================
    ! PPM CODE
    !==========================================================================

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

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          gfactor = ONE ! to help compiler resolve ANTI dependence

          rho = q(i,j,k3d,QRHO)

          cc = c(i,j,k3d)
          csq = cc**2
          Clag = rho*cc

          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)

          p = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          enth = (rhoe_g+p)/(rho*csq)

          game = q(i,j,k3d,QGAME)

          gam = gamc(i,j,k3d)


          !--------------------------------------------------------------------
          ! plus state on face i
          !--------------------------------------------------------------------

          if (i .ge. ilo1) then

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. u - cc >= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                u_ref    = u

                p_ref    = p
                rhoe_g_ref = rhoe_g

                tau_ref  = ONE/rho

                gam_ref  = gam

                game_ref = game

             else
                ! This will be the fastest moving state to the left --
                ! this is the method that Miller & Colella and Colella &
                ! Woodward use
                rho_ref  = Im(i,j,kc,1,1,QRHO)
                u_ref    = Im(i,j,kc,1,1,QU)

                p_ref    = Im(i,j,kc,1,1,QPRES)
                rhoe_g_ref = Im(i,j,kc,1,1,QREINT)

                tau_ref  = ONE/Im(i,j,kc,1,1,QRHO)

                gam_ref  = Im_gc(i,j,kc,1,1,1)

                game_ref = Im(i,j,kc,1,1,QGAME)
             endif

             rho_ref = max(rho_ref,small_dens)
             p_ref = max(p_ref,small_pres)

             ! For tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = (rhoe_g_ref+p_ref)/(rho_ref*csq_ref)

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c

             ! Note: for the transverse velocities, the jump is carried
             !       only by the u wave (the contact)

             dum   = u_ref    - Im(i,j,kc,1,1,QU)
             dptotm   = p_ref    - Im(i,j,kc,1,1,QPRES)

             drho  = rho_ref  - Im(i,j,kc,1,2,QRHO)
             dptot    = p_ref    - Im(i,j,kc,1,2,QPRES)
             drhoe_g = rhoe_g_ref - Im(i,j,kc,1,2,QREINT)
             dtau  = tau_ref  - ONE/Im(i,j,kc,1,2,QRHO)

             dup   = u_ref    - Im(i,j,kc,1,3,QU)
             dptotp   = p_ref    - Im(i,j,kc,1,3,QPRES)


             ! If we are doing source term tracing, then add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_sources .eq. 1) then
                dum = dum - halfdt*Im_src(i,j,kc,1,1,QU)
                dup = dup - halfdt*Im_src(i,j,kc,1,3,QU)
             endif

             ! Optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
                tau_ev  = ONE/rho
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
                tau_ev  = tau_ref
             endif

             if (ppm_tau_in_tracing == 0) then

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is
                ! simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dum)*(rho_ev/cc_ev)
                alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dup)*(rho_ev/cc_ev)
                alpha0r = drho - dptot/csq_ev
                alpha0e = drhoe_g - dptot*enth_ev  ! note enth has a 1/c**2 in it

             else

                ! (tau, u, p, e) eigensystem
                ! or
                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_g_ref/rho_ref - Im(i,j,kc,1,2,QREINT)/Im(i,j,kc,1,2,QRHO))
                dge   = game_ref - Im(i,j,kc,1,2,QGAME)

                alpham = HALF*( dum - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                alphap = HALF*(-dup - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                if (ppm_predict_gammae == 0) then
                   alpha0e = de - dptot*p_ev*(ONE/Clag_ev)**2
                else
                   gfactor = (game - ONE)*(game - gam)
                   alpha0e = gfactor*dptot/(tau_ev*Clag_ev**2) + dge
                endif

             endif    ! which tracing method

             if (u-cc .gt. ZERO) then
                alpham = ZERO
             else if (u-cc .lt. ZERO) then
                alpham = -alpham
             else
                alpham = -HALF*alpham
             endif

             if (u+cc .gt. ZERO) then
                alphap = ZERO
             else if (u+cc .lt. ZERO) then
                alphap = -alphap
             else
                alphap = -HALF*alphap
             endif

             if (u .gt. ZERO) then
                alpha0r = ZERO
                alpha0e = ZERO
             else if (u .lt. ZERO) then
                alpha0r = -alpha0r
                alpha0e = -alpha0e
             else
                alpha0r = -HALF*alpha0r
                alpha0e = -HALF*alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum(l . dq) r
             ! note that the a{mpz}right as defined above have the minus already
             if (ppm_tau_in_tracing == 0) then
                qxp(i,j,kc,QRHO  ) =  rho_ref +  alphap + alpham + alpha0r
                qxp(i,j,kc,QU    ) =    u_ref + (alphap - alpham)*cc_ev/rho_ev
                qxp(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*enth_ev*csq_ev + alpha0e
                qxp(i,j,kc,QPRES ) =    p_ref + (alphap + alpham)*csq_ev
             else
                tau_s = tau_ref + alphap + alpham + alpha0r
                qxp(i,j,kc,QRHO  ) = ONE/tau_s
                qxp(i,j,kc,QU    ) = u_ref + (alpham - alphap)*Clag_ev

                qxp(i,j,kc,QPRES ) = p_ref + (-alphap - alpham)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_g_ref/rho_ref + (alpha0e - p_ev*alpham - p_ev*alphap)
                   qxp(i,j,kc,QREINT) = e_s*(ONE/tau_s)
                else
                   qxp(i,j,kc,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e
                   qxp(i,j,kc,QREINT) = qxp(i,j,kc,QPRES )/(qxp(i,j,kc,QGAME) - ONE)
                endif
             endif

             ! Enforce small_*
             qxp(i,j,kc,QRHO ) = max(qxp(i,j,kc,QRHO ),small_dens)
             qxp(i,j,kc,QPRES) = max(qxp(i,j,kc,QPRES),small_pres)

             ! Transverse velocities -- there's no projection here, so
             ! we don't need a reference state.  We only care about
             ! the state traced under the middle wave

             ! Recall that I already takes the limit of the parabola
             ! in the event that the wave is not moving toward the
             ! interface
             if (u > ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qxp(i,j,kc,QV    ) = Im(i,j,kc,1,2,QV)
                   qxp(i,j,kc,QW    ) = Im(i,j,kc,1,2,QW)
                else
                   qxp(i,j,kc,QV    ) = v
                   qxp(i,j,kc,QW    ) = w
                endif
             else ! wave moving toward the interface
                qxp(i,j,kc,QV) = Im(i,j,kc,1,2,QV)
                qxp(i,j,kc,QW) = Im(i,j,kc,1,2,QW)
             endif

             if (ppm_trace_sources .eq. 1) then
                qxp(i,j,kc,QV) = qxp(i,j,kc,QV) + halfdt*Im_src(i,j,kc,1,2,QV)
                qxp(i,j,kc,QW) = qxp(i,j,kc,QW) + halfdt*Im_src(i,j,kc,1,2,QW)
             endif

          end if


          !--------------------------------------------------------------------
          ! minus state on face i + 1
          !--------------------------------------------------------------------
          if (i .le. ihi1) then

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. u + cc <= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc values
                rho_ref  = rho
                u_ref    = u

                p_ref    = p
                rhoe_g_ref = rhoe_g

                tau_ref  = ONE/rho

                gam_ref  = gam

                game_ref = game

             else
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,kc,1,3,QRHO)
                u_ref    = Ip(i,j,kc,1,3,QU)

                p_ref    = Ip(i,j,kc,1,3,QPRES)
                rhoe_g_ref = Ip(i,j,kc,1,3,QREINT)

                tau_ref  = ONE/Ip(i,j,kc,1,3,QRHO)

                gam_ref  = Ip_gc(i,j,kc,1,3,1)

                game_ref = Ip(i,j,kc,1,3,QGAME)
             endif

             rho_ref = max(rho_ref,small_dens)
             p_ref = max(p_ref,small_pres)

             ! For tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = (rhoe_g_ref+p_ref)/(rho_ref*csq_ref)

             ! *m are the jumps carried by u-c
             ! *p are the jumps carried by u+c

             dum   = u_ref    - Ip(i,j,kc,1,1,QU)
             dptotm   = p_ref    - Ip(i,j,kc,1,1,QPRES)

             drho  = rho_ref  - Ip(i,j,kc,1,2,QRHO)
             dptot    = p_ref    - Ip(i,j,kc,1,2,QPRES)
             drhoe_g = rhoe_g_ref - Ip(i,j,kc,1,2,QREINT)
             dtau  = tau_ref  - ONE/Ip(i,j,kc,1,2,QRHO)

             dup   = u_ref    - Ip(i,j,kc,1,3,QU)
             dptotp   = p_ref    - Ip(i,j,kc,1,3,QPRES)

             ! If we are doing source term tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_sources .eq. 1) then
                dum = dum - halfdt*Ip_src(i,j,kc,1,1,QU)
                dup = dup - halfdt*Ip_src(i,j,kc,1,3,QU)
             endif


             ! Optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
                tau_ev  = ONE/rho
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
                tau_ev  = tau_ref
             endif

             if (ppm_tau_in_tracing == 0) then

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is
                ! simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dum)*(rho_ev/cc_ev)
                alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dup)*(rho_ev/cc_ev)
                alpha0r = drho - dptot/csq_ev
                alpha0e = drhoe_g - dptot*enth_ev  ! enth has a 1/c**2 in it

             else
                ! (tau, u, p, e) eigensystem
                ! or
                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_g_ref/rho_ref - Ip(i,j,kc,1,2,QREINT)/Ip(i,j,kc,1,2,QRHO))
                dge = game_ref - Ip(i,j,kc,1,2,QGAME)

                alpham = HALF*( dum - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                alphap = HALF*(-dup - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                if (ppm_predict_gammae == 0) then
                   alpha0e = de - dptot*p_ev*(ONE/Clag_ev)**2
                else
                   gfactor = (game - ONE)*(game - gam)
                   alpha0e = gfactor*dptot/(tau_ev*Clag_ev**2) + dge
                endif

             end if

             if (u-cc .gt. ZERO) then
                alpham = -alpham
             else if (u-cc .lt. ZERO) then
                alpham = ZERO
             else
                alpham = -HALF*alpham
             endif

             if (u+cc .gt. ZERO) then
                alphap = -alphap
             else if (u+cc .lt. ZERO) then
                alphap = ZERO
             else
                alphap = -HALF*alphap
             endif

             if (u .gt. ZERO) then
                alpha0r = -alpha0r
                alpha0e = -alpha0e
             else if (u .lt. ZERO) then
                alpha0r = ZERO
                alpha0e = ZERO
             else
                alpha0r = -HALF*alpha0r
                alpha0e = -HALF*alpha0e
             endif


             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}left as defined above have the minus already
             if (ppm_tau_in_tracing == 0) then
                qxm(i+1,j,kc,QRHO  ) =  rho_ref +  alphap + alpham + alpha0r
                qxm(i+1,j,kc,QU    ) =    u_ref + (alphap - alpham)*cc_ev/rho_ev
                qxm(i+1,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*enth_ev*csq_ev + alpha0e
                qxm(i+1,j,kc,QPRES ) =    p_ref + (alphap + alpham)*csq_ev
             else
                tau_s = tau_ref + (alphap + alpham + alpha0r)
                qxm(i+1,j,kc,QRHO  ) = ONE/tau_s
                qxm(i+1,j,kc,QU    ) = u_ref + (alpham - alphap)*Clag_ev

                qxm(i+1,j,kc,QPRES ) = p_ref + (-alphap - alpham)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_g_ref/rho_ref + (alpha0e - p_ev*alpham - p_ev*alphap)
                   qxm(i+1,j,kc,QREINT) = e_s*(ONE/tau_s)
                else
                   qxm(i+1,j,kc,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e
                   qxm(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QPRES )/(qxm(i+1,j,kc,QGAME) - ONE)
                endif
             endif

             ! Enforce small_*
             qxm(i+1,j,kc,QRHO ) = max(qxm(i+1,j,kc,QRHO ),small_dens)
             qxm(i+1,j,kc,QPRES) = max(qxm(i+1,j,kc,QPRES),small_pres)

             ! transverse velocities
             if (u < ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qxm(i+1,j,kc,QV    ) = Ip(i,j,kc,1,2,QV)
                   qxm(i+1,j,kc,QW    ) = Ip(i,j,kc,1,2,QW)
                else
                   qxm(i+1,j,kc,QV    ) = v
                   qxm(i+1,j,kc,QW    ) = w
                endif
             else ! wave moving toward interface
                qxm(i+1,j,kc,QV    ) = Ip(i,j,kc,1,2,QV)
                qxm(i+1,j,kc,QW    ) = Ip(i,j,kc,1,2,QW)
             endif

             if (ppm_trace_sources .eq. 1) then
                qxm(i+1,j,kc,QV) = qxm(i+1,j,kc,QV) + halfdt*Ip_src(i,j,kc,1,2,QV)
                qxm(i+1,j,kc,QW) = qxm(i+1,j,kc,QW) + halfdt*Ip_src(i,j,kc,1,2,QW)
             endif

          end if

       end do
    end do


    !--------------------------------------------------------------------------
    ! passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1

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

             if (u .gt. ZERO) then
                qxp(i,j,kc,n) = q(i,j,k3d,n)
             else if (u .lt. ZERO) then
                qxp(i,j,kc,n) = Im(i,j,kc,1,2,n)
             else
                qxp(i,j,kc,n) = q(i,j,k3d,n) + HALF*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo

          ! Minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)

             if (u .gt. ZERO) then
                qxm(i+1,j,kc,n) = Ip(i,j,kc,1,2,n)
             else if (u .lt. ZERO) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n)
             else
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) + HALF*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo
       enddo
    enddo


    !--------------------------------------------------------------------------
    ! y-direction
    !--------------------------------------------------------------------------

    ! Trace to bottom and top edges using upwind PPM

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          gfactor = ONE ! to help compiler resolve ANTI dependence

          rho = q(i,j,k3d,QRHO)

          cc = c(i,j,k3d)
          csq = cc**2
          Clag = rho*cc

          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)

          p = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          enth = ( (rhoe_g+p)/rho )/csq

          gam = gamc(i,j,k3d)

          game = q(i,j,k3d,QGAME)

          !--------------------------------------------------------------------
          ! plus state on face j
          !--------------------------------------------------------------------

          if (j .ge. ilo2) then

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. v - cc >= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                v_ref    = v

                p_ref    = p
                rhoe_g_ref = rhoe_g

                tau_ref  = ONE/rho

                gam_ref  = gam

                game_ref = game

             else
                ! This will be the fastest moving state to the left
                rho_ref  = Im(i,j,kc,2,1,QRHO)
                v_ref    = Im(i,j,kc,2,1,QV)

                p_ref    = Im(i,j,kc,2,1,QPRES)
                rhoe_g_ref = Im(i,j,kc,2,1,QREINT)

                tau_ref  = ONE/Im(i,j,kc,2,1,QRHO)
                gam_ref  = Im_gc(i,j,kc,2,1,1)

                game_ref = Im(i,j,kc,2,1,QGAME)
             endif

             rho_ref = max(rho_ref,small_dens)
             p_ref = max(p_ref,small_pres)

             ! For tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = (rhoe_g_ref+p_ref)/(rho_ref*csq_ref)

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c

             dvm   = v_ref    - Im(i,j,kc,2,1,QV)
             dptotm   = p_ref    - Im(i,j,kc,2,1,QPRES)

             drho  = rho_ref  - Im(i,j,kc,2,2,QRHO)
             dptot    = p_ref    - Im(i,j,kc,2,2,QPRES)
             drhoe_g = rhoe_g_ref - Im(i,j,kc,2,2,QREINT)
             dtau  = tau_ref  - ONE/Im(i,j,kc,2,2,QRHO)

             dvp   = v_ref    - Im(i,j,kc,2,3,QV)
             dptotp   = p_ref    - Im(i,j,kc,2,3,QPRES)

             ! If we are doing source term tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_sources .eq. 1) then
                dvm = dvm - halfdt*Im_src(i,j,kc,2,1,QV)
                dvp = dvp - halfdt*Im_src(i,j,kc,2,3,QV)
             endif

             ! Optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
                tau_ev  = ONE/rho
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
                tau_ev  = tau_ref
             endif

             if (ppm_tau_in_tracing == 0) then

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This
                ! is simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dvm)*(rho_ev/cc_ev)
                alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dvp)*(rho_ev/cc_ev)
                alpha0r = drho - dptot/csq_ev
                alpha0e = drhoe_g - dptot*enth_ev

             else

                ! (tau, u, p, e) eigensystem
                ! or
                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_g_ref/rho_ref - Im(i,j,kc,2,2,QREINT)/Im(i,j,kc,2,2,QRHO))
                dge = game_ref - Im(i,j,kc,2,2,QGAME)

                alpham = HALF*( dvm - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                alphap = HALF*(-dvp - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                if (ppm_predict_gammae == 0) then
                   alpha0e = de - dptot*p_ev*(ONE/Clag_ev)**2
                else
                   gfactor = (game - ONE)*(game - gam)
                   alpha0e = gfactor*dptot/(tau_ev*Clag_ev**2) + dge
                endif

             end if

             if (v-cc .gt. ZERO) then
                alpham = ZERO
             else if (v-cc .lt. ZERO) then
                alpham = -alpham
             else
                alpham = -HALF*alpham
             endif

             if (v+cc .gt. ZERO) then
                alphap = ZERO
             else if (v+cc .lt. ZERO) then
                alphap = -alphap
             else
                alphap = -HALF*alphap
             endif

             if (v .gt. ZERO) then
                alpha0r = ZERO
                alpha0e = ZERO
             else if (v .lt. ZERO) then
                alpha0r = -alpha0r
                alpha0e = -alpha0e
             else
                alpha0r = -HALF*alpha0r
                alpha0e = -HALF*alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}right as defined above have the minus already
             if (ppm_tau_in_tracing == 0) then
                qyp(i,j,kc,QRHO  ) = rho_ref + alphap + alpham + alpha0r
                qyp(i,j,kc,QV    ) = v_ref + (alphap - alpham)*cc_ev/rho_ev
                qyp(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*enth_ev*csq_ev + alpha0e
                qyp(i,j,kc,QPRES ) = p_ref + (alphap + alpham)*csq_ev
             else
                tau_s = tau_ref + alphap + alpham + alpha0r
                qyp(i,j,kc,QRHO  ) = ONE/tau_s
                qyp(i,j,kc,QV    ) = v_ref + (alpham - alphap)*Clag_ev

                qyp(i,j,kc,QPRES ) = p_ref + (-alphap - alpham)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_g_ref/rho_ref + (alpha0e - p_ev*alpham - p_ev*alphap)
                   qyp(i,j,kc,QREINT) = e_s*(ONE/tau_s)
                else
                   qyp(i,j,kc,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e
                   qyp(i,j,kc,QREINT) = qyp(i,j,kc,QPRES )/(qyp(i,j,kc,QGAME) - ONE)
                endif
             endif

             ! Enforce small_*
             qyp(i,j,kc,QRHO ) = max(qyp(i,j,kc,QRHO ),small_dens)
             qyp(i,j,kc,QPRES) = max(qyp(i,j,kc,QPRES),small_pres)

             ! transverse velocities
             if (v > ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qyp(i,j,kc,QU    ) = Im(i,j,kc,2,2,QU)
                   qyp(i,j,kc,QW    ) = Im(i,j,kc,2,2,QW)
                else
                   qyp(i,j,kc,QU    ) = u
                   qyp(i,j,kc,QW    ) = w
                endif
             else ! wave moving toward the interface
                qyp(i,j,kc,QU    ) = Im(i,j,kc,2,2,QU)
                qyp(i,j,kc,QW    ) = Im(i,j,kc,2,2,QW)
             endif

             if (ppm_trace_sources .eq. 1) then
                qyp(i,j,kc,QU) = qyp(i,j,kc,QU) + halfdt*Im_src(i,j,kc,2,2,QU)
                qyp(i,j,kc,QW) = qyp(i,j,kc,QW) + halfdt*Im_src(i,j,kc,2,2,QW)
             endif

          end if

          !--------------------------------------------------------------------
          ! minus state on face j+1
          !--------------------------------------------------------------------

          if (j .le. ihi2) then

             ! Set the reference state
             if (ppm_reference == 0 .or. &
                  (ppm_reference == 1 .and. v + cc <= ZERO .and. &
                   ppm_reference_edge_limit == 0) ) then
                ! original Castro way -- cc value
                rho_ref  = rho
                v_ref    = v

                p_ref    = p
                rhoe_g_ref = rhoe_g

                tau_ref  = ONE/rho

                gam_ref  = gam

                game_ref = game

             else
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(i,j,kc,2,3,QRHO)
                v_ref    = Ip(i,j,kc,2,3,QV)

                p_ref    = Ip(i,j,kc,2,3,QPRES)
                rhoe_g_ref = Ip(i,j,kc,2,3,QREINT)

                tau_ref  = ONE/Ip(i,j,kc,2,3,QRHO)

                gam_ref  = Ip_gc(i,j,kc,2,3,1)

                game_ref = Ip(i,j,kc,2,3,QGAME)
             endif

             rho_ref = max(rho_ref,small_dens)
             p_ref = max(p_ref,small_pres)

             ! For tracing (optionally)
             cc_ref = sqrt(gam_ref*p_ref/rho_ref)
             csq_ref = cc_ref**2
             Clag_ref = rho_ref*cc_ref
             enth_ref = (rhoe_g_ref+p_ref)/(rho_ref*csq_ref)

             ! *m are the jumps carried by v-c
             ! *p are the jumps carried by v+c

             dvm   = v_ref    - Ip(i,j,kc,2,1,QV)
             dptotm   = p_ref    - Ip(i,j,kc,2,1,QPRES)

             drho  = rho_ref  - Ip(i,j,kc,2,2,QRHO)
             dptot    = p_ref    - Ip(i,j,kc,2,2,QPRES)
             drhoe_g = rhoe_g_ref - Ip(i,j,kc,2,2,QREINT)
             dtau  = tau_ref  - ONE/Ip(i,j,kc,2,2,QRHO)

             dvp   = v_ref    - Ip(i,j,kc,2,3,QV)
             dptotp   = p_ref    - Ip(i,j,kc,2,3,QPRES)

             ! If we are doing source term tracing, then we add the force
             ! to the velocity here, otherwise we will deal with this
             ! in the trans_X routines
             if (ppm_trace_sources .eq. 1) then
                dvm = dvm - halfdt*Ip_src(i,j,kc,2,1,QV)
                dvp = dvp - halfdt*Ip_src(i,j,kc,2,3,QV)
             endif


             ! Optionally use the reference state in evaluating the
             ! eigenvectors
             if (ppm_reference_eigenvectors == 0) then
                rho_ev  = rho
                cc_ev   = cc
                csq_ev  = csq
                Clag_ev = Clag
                enth_ev = enth
                p_ev    = p
                tau_ev  = ONE/rho
             else
                rho_ev  = rho_ref
                cc_ev   = cc_ref
                csq_ev  = csq_ref
                Clag_ev = Clag_ref
                enth_ev = enth_ref
                p_ev    = p_ref
                tau_ev  = tau_ref
             endif

             if (ppm_tau_in_tracing == 0) then

                ! These are analogous to the beta's from the original PPM
                ! paper.  This is simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dvm)*(rho_ev/cc_ev)
                alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dvp)*(rho_ev/cc_ev)
                alpha0r = drho - dptot/csq_ev
                alpha0e = drhoe_g - dptot*enth_ev

             else
                ! (tau, u, p, e) eigensystem
                ! or
                ! (tau, u, p, game) eigensystem

                ! This is the way things were done in the original PPM
                ! paper -- here we work with tau in the characteristic
                ! system
                de = (rhoe_g_ref/rho_ref - Ip(i,j,kc,2,2,QREINT)/Ip(i,j,kc,2,2,QRHO))
                dge = game_ref - Ip(i,j,kc,2,2,QGAME)

                alpham = HALF*( dvm - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
                alphap = HALF*(-dvp - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
                alpha0r = dtau + dptot*(ONE/Clag_ev)**2

                if (ppm_predict_gammae == 0) then
                   alpha0e = de - dptot*p_ev*(ONE/Clag_ev)**2
                else
                   gfactor = (game - ONE)*(game - gam)
                   alpha0e = gfactor*dptot/(tau_ev*Clag_ev**2) + dge
                endif

             end if

             if (v-cc .gt. ZERO) then
                alpham = -alpham
             else if (v-cc .lt. ZERO) then
                alpham = ZERO
             else
                alpham = -HALF*alpham
             endif

             if (v+cc .gt. ZERO) then
                alphap = -alphap
             else if (v+cc .lt. ZERO) then
                alphap = ZERO
             else
                alphap = -HALF*alphap
             endif

             if (v .gt. ZERO) then
                alpha0r = -alpha0r
                alpha0e = -alpha0e
             else if (v .lt. ZERO) then
                alpha0r = ZERO
                alpha0e = ZERO
             else
                alpha0r = -HALF*alpha0r
                alpha0e = -HALF*alpha0e
             endif

             ! The final interface states are just
             ! q_s = q_ref - sum (l . dq) r
             ! note that the a{mpz}left as defined above has the minus already
             if (ppm_tau_in_tracing == 0) then
                qym(i,j+1,kc,QRHO  ) = rho_ref + alphap + alpham + alpha0r
                qym(i,j+1,kc,QV    ) = v_ref + (alphap - alpham)*cc_ev/rho_ev
                qym(i,j+1,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*enth_ev*csq_ev + alpha0e
                qym(i,j+1,kc,QPRES ) = p_ref + (alphap + alpham)*csq_ev
             else
                tau_s = tau_ref + alphap + alpham + alpha0r
                qym(i,j+1,kc,QRHO  ) = ONE/tau_s
                qym(i,j+1,kc,QV    ) = v_ref + (alpham - alphap)*Clag_ev

                qym(i,j+1,kc,QPRES ) = p_ref + (-alphap - alpham)*Clag_ev**2

                if (ppm_predict_gammae == 0) then
                   e_s = rhoe_g_ref/rho_ref + (alpha0e - p_ev*alpham - p_ev*alphap)
                   qym(i,j+1,kc,QREINT) = e_s*(ONE/tau_s)
                else
                   qym(i,j+1,kc,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e
                   qym(i,j+1,kc,QREINT) = qym(i,j+1,kc,QPRES )/(qym(i,j+1,kc,QGAME) - ONE)
                endif

             endif

             ! Enforce small_*
             qym(i,j+1,kc,QRHO ) = max(qym(i,j+1,kc,QRHO ),small_dens)
             qym(i,j+1,kc,QPRES) = max(qym(i,j+1,kc,QPRES),small_pres)

             ! transverse velocities
             if (v < ZERO) then
                if (ppm_reference_edge_limit == 1) then
                   qym(i,j+1,kc,QU    ) = Ip(i,j,kc,2,2,QU)
                   qym(i,j+1,kc,QW    ) = Ip(i,j,kc,2,2,QW)
                else
                   qym(i,j+1,kc,QU    ) = u
                   qym(i,j+1,kc,QW    ) = w
                endif
             else ! wave is moving toward the interface
                qym(i,j+1,kc,QU    ) = Ip(i,j,kc,2,2,QU)
                qym(i,j+1,kc,QW    ) = Ip(i,j,kc,2,2,QW)
             endif

             if (ppm_trace_sources .eq. 1) then
                qym(i,j+1,kc,QU) = qym(i,j+1,kc,QU) + halfdt*Ip_src(i,j,kc,2,2,QU)
                qym(i,j+1,kc,QW) = qym(i,j+1,kc,QW) + halfdt*Ip_src(i,j,kc,2,2,QW)
             endif

          end if

       end do
    end do

    !--------------------------------------------------------------------------
    ! Passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! Plus state on face j
       do j = ilo2, ihi2+1
          do i = ilo1-1, ihi1+1
             v = q(i,j,k3d,QV)

             if (v .gt. ZERO) then
                qyp(i,j,kc,n) = q(i,j,k3d,n)
             else if (v .lt. ZERO) then
                qyp(i,j,kc,n) = Im(i,j,kc,2,2,n)
             else
                qyp(i,j,kc,n) = q(i,j,k3d,n) + HALF*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo
       end do

       ! Minus state on face j+1
       do j = ilo2-1, ihi2
          do i = ilo1-1, ihi1+1
             v = q(i,j,k3d,QV)

             if (v .gt. ZERO) then
                qym(i,j+1,kc,n) = Ip(i,j,kc,2,2,n)
             else if (v .lt. ZERO) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n)
             else
                qym(i,j+1,kc,n) = q(i,j,k3d,n) + HALF*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo

       enddo
    enddo

  end subroutine tracexy_ppm



  subroutine tracez_ppm(q,c,flatn,qd_lo,qd_hi, &
                        Ip,Im,Ip_src,Im_src,Ip_gc,Im_gc,I_lo,I_hi, &
                        qzm,qzp,qs_lo,qs_hi, &
                        gamc,gc_lo,gc_hi, &
                        ilo1,ilo2,ihi1,ihi2,dt,km,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
         QREINT, QPRES, QGAME, &
         small_dens, small_pres, &
         ppm_type, ppm_reference, ppm_trace_sources, &
         ppm_tau_in_tracing, ppm_reference_eigenvectors, &
         ppm_reference_edge_limit, &
         ppm_predict_gammae, &
         npassive, qpass_map
    use bl_constants_module

    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qs_lo(3),qs_hi(3)
    integer, intent(in) :: gc_lo(3), gc_hi(3)
    integer, intent(in) :: I_lo(3), I_hi(3)
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    integer, intent(in) :: km, kc, k3d

    double precision, intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    double precision, intent(in) ::     c(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision, intent(in) :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))

    double precision, intent(in) :: Ip(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)
    double precision, intent(in) :: Im(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)

    double precision, intent(in) :: Ip_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)
    double precision, intent(in) :: Im_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,QVAR)

    double precision, intent(in) :: Ip_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,1)
    double precision, intent(in) :: Im_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:3,1:3,1)

    double precision, intent(inout) :: qzm(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QVAR)
    double precision, intent(inout) :: qzp(qs_lo(1):qs_hi(1),qs_lo(2):qs_hi(2),qs_lo(3):qs_hi(3),QVAR)

    double precision, intent(in) :: gamc(gc_lo(1):gc_hi(1),gc_lo(2):gc_hi(2),gc_lo(3):gc_hi(3))

    double precision, intent(in) :: dt

    !     Local variables
    integer i, j
    integer n
    integer ipassive

    double precision cc, csq, Clag, rho, u, v, w, p, rhoe_g, enth

    double precision drho, du, dv, dptot, drhoe_g, de, dge, dtau
    double precision dwp, dptotp
    double precision dwm, dptotm

    double precision :: rho_ref, w_ref, p_ref, rhoe_g_ref, tau_ref
    double precision :: tau_s, e_s

    double precision :: cc_ref, csq_ref, Clag_ref, enth_ref, gam_ref, game_ref, gfactor
    double precision :: cc_ev, csq_ev, Clag_ev, rho_ev, p_ev, enth_ev, tau_ev
    double precision :: gam, game

    double precision alpham, alphap, alpha0r, alpha0e

    double precision halfdt

    halfdt = HALF * dt

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call bl_error("Error:: trace_ppm_3d.f90 :: tracez_ppm")
    end if

    !==========================================================================
    ! PPM CODE
    !==========================================================================

    ! Trace to left and right edges using upwind PPM
    !
    ! Note: in contrast to the above code for x and y, here the loop
    ! is over interfaces, not over cell-centers.


    !--------------------------------------------------------------------------
    ! construct qzp  -- plus state on face kc
    !--------------------------------------------------------------------------
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          gfactor = ONE ! to help compiler resolve ANTI dependence

          rho  = q(i,j,k3d,QRHO)

          cc   = c(i,j,k3d)
          csq  = cc**2
          Clag = rho*cc

          u    = q(i,j,k3d,QU)
          v    = q(i,j,k3d,QV)
          w    = q(i,j,k3d,QW)

          p    = q(i,j,k3d,QPRES)
          rhoe_g = q(i,j,k3d,QREINT)
          enth = ( (p+rhoe_g)/rho )/csq

          gam = gamc(i,j,k3d)

          game = q(i,j,k3d,QGAME)

          ! Set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. w - cc >= ZERO .and. &
                ppm_reference_edge_limit == 0) ) then
             ! original Castro way -- cc value
             rho_ref  = rho
             w_ref    = w

             p_ref    = p
             rhoe_g_ref = rhoe_g

             tau_ref  = ONE/rho

             gam_ref  = gam

             game_ref = game

          else
             ! This will be the fastest moving state to the left
             rho_ref  = Im(i,j,kc,3,1,QRHO)
             w_ref    = Im(i,j,kc,3,1,QW)

             p_ref    = Im(i,j,kc,3,1,QPRES)
             rhoe_g_ref = Im(i,j,kc,3,1,QREINT)

             tau_ref  = ONE/Im(i,j,kc,3,1,QRHO)
             gam_ref  = Im_gc(i,j,kc,3,1,1)

             game_ref = Im(i,j,kc,3,1,QGAME)
          endif

          rho_ref = max(rho_ref,small_dens)
          p_ref = max(p_ref,small_pres)

          ! For tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = (rhoe_g_ref+p_ref)/(rho_ref*csq_ref)

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   = w_ref    - Im(i,j,kc,3,1,QW)
          dptotm   = p_ref    - Im(i,j,kc,3,1,QPRES)

          drho  = rho_ref  - Im(i,j,kc,3,2,QRHO)
          dptot    = p_ref    - Im(i,j,kc,3,2,QPRES)
          drhoe_g = rhoe_g_ref - Im(i,j,kc,3,2,QREINT)
          dtau  = tau_ref  - ONE/Im(i,j,kc,3,2,QRHO)

          dwp   = w_ref    - Im(i,j,kc,3,3,QW)
          dptotp   = p_ref    - Im(i,j,kc,3,3,QPRES)

          ! If we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources .eq. 1) then
             dwm = dwm - halfdt*Im_src(i,j,kc,3,1,QW)
             dwp = dwp - halfdt*Im_src(i,j,kc,3,3,QW)
          endif

          ! Optionally use the reference state in evaluating the
          ! eigenvectors
          if (ppm_reference_eigenvectors == 0) then
             rho_ev  = rho
             cc_ev   = cc
             csq_ev  = csq
             Clag_ev = Clag
             enth_ev = enth
             p_ev    = p
             tau_ev  = ONE/rho
          else
             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             enth_ev = enth_ref
             p_ev    = p_ref
             tau_ev  = tau_ref
          endif

          if (ppm_tau_in_tracing == 0) then

             ! These are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)
             alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dwm)*(rho_ev/cc_ev)
             alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dwp)*(rho_ev/cc_ev)
             alpha0r = drho - dptot/csq_ev
             alpha0e = drhoe_g - dptot*enth_ev

          else
             ! (tau, u, p, e) eigensystem
             ! or
             ! (tau, u, p, game) eigensystem

             ! This is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system
             de = (rhoe_g_ref/rho_ref - Im(i,j,kc,3,2,QREINT)/Im(i,j,kc,3,2,QRHO))
             dge = game_ref - Im(i,j,kc,3,2,QGAME)

             alpham = HALF*( dwm - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
             alphap = HALF*(-dwp - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
             alpha0r = dtau + dptot*(ONE/Clag_ev)**2

             if (ppm_predict_gammae == 0) then
                alpha0e = de - dptot*p_ev*(ONE/Clag_ev)**2
             else
                gfactor = (game - ONE)*(game - gam)
                alpha0e = gfactor*dptot/(tau_ev*Clag_ev**2) + dge
             endif

          endif

          if (w-cc .gt. ZERO) then
             alpham = ZERO
          else if (w-cc .lt. ZERO) then
             alpham = -alpham
          else
             alpham = -HALF*alpham
          endif
          if (w+cc .gt. ZERO) then
             alphap = ZERO
          else if (w+cc .lt. ZERO) then
             alphap = -alphap
          else
             alphap = -HALF*alphap
          endif
          if (w .gt. ZERO) then
             alpha0r = ZERO
             alpha0e = ZERO
          else if (w .lt. ZERO) then
             alpha0r = -alpha0r
             alpha0e = -alpha0e
          else
             alpha0r = -HALF*alpha0r
             alpha0e = -HALF*alpha0e
          endif

          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          ! note that the a{mpz}right as defined above have the minus already
          if (ppm_tau_in_tracing == 0) then
             qzp(i,j,kc,QRHO  ) = rho_ref + alphap + alpham + alpha0r
             qzp(i,j,kc,QW    ) = w_ref + (alphap - alpham)*cc_ev/rho_ev
             qzp(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*enth_ev*csq_ev + alpha0e
             qzp(i,j,kc,QPRES ) = p_ref + (alphap + alpham)*csq_ev
          else
             tau_s = tau_ref + alphap + alpham + alpha0r
             qzp(i,j,kc,QRHO  ) = ONE/tau_s
             qzp(i,j,kc,QW    ) = w_ref + (alpham - alphap)*Clag_ev

             qzp(i,j,kc,QPRES ) = p_ref + (-alphap - alpham)*Clag_ev**2

             if (ppm_predict_gammae == 0) then
                e_s = rhoe_g_ref/rho_ref + (alpha0e - p_ev*alpham - p_ev*alphap)
                qzp(i,j,kc,QREINT) = e_s*(ONE/tau_s)
             else
                qzp(i,j,kc,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e
                qzp(i,j,kc,QREINT) = qzp(i,j,kc,QPRES )/(qzp(i,j,kc,QGAME) - ONE)
             endif

          endif

          ! Enforce small_*
          qzp(i,j,kc,QRHO ) = max(qzp(i,j,kc,QRHO ),small_dens)
          qzp(i,j,kc,QPRES) = max(qzp(i,j,kc,QPRES),small_pres)

          ! transverse velocities
          if (w > ZERO) then
             if (ppm_reference_edge_limit == 1) then
                qzp(i,j,kc,QU    ) = Im(i,j,kc,3,2,QU)
                qzp(i,j,kc,QV    ) = Im(i,j,kc,3,2,QV)
             else
                qzp(i,j,kc,QU    ) = u
                qzp(i,j,kc,QV    ) = v
             endif
          else ! wave moving toward the interface
             qzp(i,j,kc,QU    ) = Im(i,j,kc,3,2,QU)
             qzp(i,j,kc,QV    ) = Im(i,j,kc,3,2,QV)
          endif

          if (ppm_trace_sources .eq. 1) then
             qzp(i,j,kc,QU) = qzp(i,j,kc,QU) + halfdt*Im_src(i,j,kc,3,2,QU)
             qzp(i,j,kc,QV) = qzp(i,j,kc,QV) + halfdt*Im_src(i,j,kc,3,2,QV)
          endif


          !--------------------------------------------------------------------
          ! This is all for qzm -- minus state on face kc
          !--------------------------------------------------------------------

          ! Note this is different from how we do 1D, 2D, and the
          ! x and y-faces in 3D, where the analogous thing would have
          ! been to find the minus state on face kc+1

          rho  = q(i,j,k3d-1,QRHO)

          cc   = c(i,j,k3d-1)
          csq  = cc**2
          Clag = rho*cc

          u    = q(i,j,k3d-1,QU)
          v    = q(i,j,k3d-1,QV)
          w    = q(i,j,k3d-1,QW)

          p    = q(i,j,k3d-1,QPRES)
          rhoe_g = q(i,j,k3d-1,QREINT)
          enth = (rhoe_g+p)/(rho*csq)

          gam = gamc(i,j,k3d-1)

          game = q(i,j,k3d-1,QGAME)

          ! Set the reference state
          if (ppm_reference == 0 .or. &
               (ppm_reference == 1 .and. w + cc <= ZERO .and. &
                ppm_reference_edge_limit == 0) ) then
             rho_ref  = rho
             w_ref    = w

             p_ref    = p
             rhoe_g_ref = rhoe_g

             tau_ref  = ONE/rho

             gam_ref  = gam

             game_ref = game

          else
             ! This will be the fastest moving state to the right
             rho_ref  = Ip(i,j,km,3,3,QRHO)
             w_ref    = Ip(i,j,km,3,3,QW)

             p_ref    = Ip(i,j,km,3,3,QPRES)
             rhoe_g_ref = Ip(i,j,km,3,3,QREINT)

             tau_ref  = ONE/Ip(i,j,km,3,3,QRHO)

             gam_ref  = Ip_gc(i,j,km,3,3,1)

             game_ref = Ip(i,j,km,3,3,QGAME)
          endif

          rho_ref = max(rho_ref,small_dens)
          p_ref = max(p_ref,small_pres)

          ! For tracing (optionally)
          cc_ref = sqrt(gam_ref*p_ref/rho_ref)
          csq_ref = cc_ref**2
          Clag_ref = rho_ref*cc_ref
          enth_ref = (rhoe_g_ref+p_ref)/(rho_ref*csq_ref)

          ! *m are the jumps carried by w-c
          ! *p are the jumps carried by w+c

          ! Note: for the transverse velocities, the jump is carried
          !       only by the w wave (the contact)

          dwm   = w_ref    - Ip(i,j,km,3,1,QW)
          dptotm   = p_ref    - Ip(i,j,km,3,1,QPRES)

          drho  = rho_ref  - Ip(i,j,km,3,2,QRHO)
          dptot    = p_ref    - Ip(i,j,km,3,2,QPRES)
          drhoe_g = rhoe_g_ref - Ip(i,j,km,3,2,QREINT)
          dtau  = tau_ref  - ONE/Ip(i,j,km,3,2,QRHO)

          dwp   = w_ref    - Ip(i,j,km,3,3,QW)
          dptotp   = p_ref    - Ip(i,j,km,3,3,QPRES)

          ! If we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources .eq. 1) then
             dwm = dwm - halfdt*Ip_src(i,j,km,3,1,QW)
             dwp = dwp - halfdt*Ip_src(i,j,km,3,3,QW)
          endif

          ! Optionally use the reference state in evaluating the
          ! eigenvectors
          if (ppm_reference_eigenvectors == 0) then
             rho_ev  = rho
             cc_ev   = cc
             csq_ev  = csq
             Clag_ev = Clag
             enth_ev = enth
             p_ev    = p
             tau_ev  = ONE/rho
          else
             rho_ev  = rho_ref
             cc_ev   = cc_ref
             csq_ev  = csq_ref
             Clag_ev = Clag_ref
             enth_ev = enth_ref
             p_ev    = p_ref
             tau_ev  = tau_ref
          endif

          if (ppm_tau_in_tracing == 0) then

             ! These are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)
             alpham = HALF*(dptotm*(ONE/(rho_ev*cc_ev)) - dwm)*(rho_ev/cc_ev)
             alphap = HALF*(dptotp*(ONE/(rho_ev*cc_ev)) + dwp)*(rho_ev/cc_ev)
             alpha0r = drho - dptot/csq_ev
             alpha0e = drhoe_g - dptot*enth_ev

          else

             ! (tau, u, p, e) eigensystem
             ! or
             ! (tau, u, p, game) eigensystem

             ! This is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system
             de = (rhoe_g_ref/rho_ref - Ip(i,j,km,3,2,QREINT)/Ip(i,j,km,3,2,QRHO))
             dge = game_ref - Ip(i,j,km,3,2,QGAME)

             alpham = HALF*( dwm - dptotm*(ONE/Clag_ev))*(ONE/Clag_ev)
             alphap = HALF*(-dwp - dptotp*(ONE/Clag_ev))*(ONE/Clag_ev)
             alpha0r = dtau + dptot*(ONE/Clag_ev)**2

             if (ppm_predict_gammae == 0) then
                alpha0e = de - dptot*p_ev*(ONE/Clag_ev)**2
             else
                gfactor = (game - ONE)*(game - gam)
                alpha0e = gfactor*dptot/(tau_ev*Clag_ev**2) + dge
             endif

          endif

          if (w-cc .gt. ZERO) then
             alpham = -alpham
          else if (w-cc .lt. ZERO) then
             alpham = ZERO
          else
             alpham = -HALF*alpham
          endif
          if (w+cc .gt. ZERO) then
             alphap = -alphap
          else if (w+cc .lt. ZERO) then
             alphap = ZERO
          else
             alphap = -HALF*alphap
          endif
          if (w .gt. ZERO) then
             alpha0r = -alpha0r
             alpha0e = -alpha0e
          else if (w .lt. ZERO) then
             alpha0r = ZERO
             alpha0e = ZERO
          else
             alpha0r = -HALF*alpha0r
             alpha0e = -HALF*alpha0e
          endif

          ! The final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          ! note that the a{mpz}left as defined above have the minus already
          if (ppm_tau_in_tracing == 0) then
             qzm(i,j,kc,QRHO  ) = rho_ref + alphap + alpham + alpha0r
             qzm(i,j,kc,QW    ) = w_ref + (alphap - alpham)*cc_ev/rho_ev
             qzm(i,j,kc,QREINT) = rhoe_g_ref + (alphap + alpham)*enth_ev*csq_ev + alpha0e
             qzm(i,j,kc,QPRES ) = p_ref + (alphap + alpham)*csq_ev
          else
             tau_s = tau_ref + alphap + alpham + alpha0r
             qzm(i,j,kc,QRHO  ) = ONE/tau_s
             qzm(i,j,kc,QW    ) = w_ref + (alpham - alphap)*Clag_ev

             qzm(i,j,kc,QPRES ) = p_ref + (-alphap - alpham)*Clag_ev**2

             if (ppm_predict_gammae == 0) then
                e_s = rhoe_g_ref/rho_ref + (alpha0e - p_ev*alpham - p_ev*alphap)
                qzm(i,j,kc,QREINT) = e_s*(ONE/tau_s)
             else
                qzm(i,j,kc,QGAME) = game_ref + gfactor*(alpham + alphap)/tau_ev + alpha0e
                qzm(i,j,kc,QREINT) = qzm(i,j,kc,QPRES )/(qzm(i,j,kc,QGAME) - ONE)
             endif

          endif

          ! Enforce small_*
          qzm(i,j,kc,QRHO ) = max(qzm(i,j,kc,QRHO ),small_dens)
          qzm(i,j,kc,QPRES) = max(qzm(i,j,kc,QPRES),small_pres)

          ! Transverse velocity
          if (w < ZERO) then
             if (ppm_reference_edge_limit == 1) then
                qzm(i,j,kc,QU    ) = Ip(i,j,km,3,2,QU)
                qzm(i,j,kc,QV    ) = Ip(i,j,km,3,2,QV)
             else
                qzm(i,j,kc,QU    ) = u
                qzm(i,j,kc,QV    ) = v
             endif
          else ! wave moving toward the interface
             qzm(i,j,kc,QU    ) = Ip(i,j,km,3,2,QU)
             qzm(i,j,kc,QV    ) = Ip(i,j,km,3,2,QV)
          endif

          if (ppm_trace_sources .eq. 1) then
             qzm(i,j,kc,QU) = qzm(i,j,kc,QU) + halfdt*Ip_src(i,j,km,3,2,QU)
             qzm(i,j,kc,QV) = qzm(i,j,kc,QV) + halfdt*Ip_src(i,j,km,3,2,QV)
          endif

       end do
    end do

    !--------------------------------------------------------------------------
    ! passively advected quantities
    !--------------------------------------------------------------------------

    ! Do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! Plus state on face kc
             w = q(i,j,k3d,QW)

             if (w .gt. ZERO) then
                qzp(i,j,kc,n) = q(i,j,k3d,n)
             else if (w .lt. ZERO) then
                qzp(i,j,kc,n) = Im(i,j,kc,3,2,n)
             else
                qzp(i,j,kc,n) = q(i,j,k3d,n) + HALF*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
             endif

             ! Minus state on face k
             w = q(i,j,k3d-1,QW)

             if (w .gt. ZERO) then
                qzm(i,j,kc,n) = Ip(i,j,km,3,2,n)
             else if (w .lt. ZERO) then
                qzm(i,j,kc,n) = q(i,j,k3d-1,n)
             else
                qzm(i,j,kc,n) = q(i,j,k3d-1,n) + HALF*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
             endif

          enddo
       enddo
    enddo

  end subroutine tracez_ppm

end module trace_ppm_module
