module transverse_module

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg

  implicit none

contains

  !===========================================================================
  ! transx routines
  !===========================================================================
  subroutine trans1_on_2states(lo, hi, &
                               idir1, idir2, &
                               q2m, q2m_lo, q2m_hi, &
                               q2mo, q2mo_lo, q2mo_hi, &
                               q2p, q2p_lo, q2p_hi, &
                               q2po, q2po_lo, q2po_hi, &
                               qaux, qa_lo, qa_hi, &
                               f1, f1_lo, f1_hi, &
#ifdef RADIATION
                               rf1, rf1_lo, rf1_hi, &
#endif
                               q1, q1_lo, q1_hi, &
#if AMREX_SPACEDIM == 2
                               area1, area1_lo, area1_hi, &
                               vol, vol_lo, vol_hi, &
#endif
                               hdt, cdtdx) bind(C, name="trans1_on_2states")

    ! here, lo and hi are the bounds of the interfaces we are looping over

    use amrex_constants_module, only : ZERO, ONE, HALF
    use network, only : nspec, naux
    use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
                                   fspace_type, comoving, &
                                   GDERADS, GDLAMS, &
                                   QCG, QGAMCG, QLAMS, &
#endif
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                   small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_predict_gammae

#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t
#if AMREX_SPACEDIM == 2
    use prob_params_module, only : mom_flux_has_p
#endif

    integer, intent(in) :: q2m_lo(3), q2m_hi(3)
    integer, intent(in) :: q2p_lo(3), q2p_hi(3)
    integer, intent(in) :: q2mo_lo(3), q2mo_hi(3)
    integer, intent(in) :: q2po_lo(3), q2po_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: f1_lo(3), f1_hi(3)
#ifdef RADIATION
    integer, intent(in) :: rf1_lo(3), rf1_hi(3)
#endif
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: lo(3), hi(3)
#if AMREX_SPACEDIM == 2
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
#endif

#ifdef RADIATION
    real(rt) :: rf1(rf1_lo(1):rf1_hi(1),rf1_lo(2):rf1_hi(2),rf1_lo(3):rf1_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in), value :: hdt, cdtdx
    integer,  intent(in), value :: idir1, idir2

    real(rt), intent(in) :: q2m(q2m_lo(1):q2m_hi(1),q2m_lo(2):q2m_hi(2),q2m_lo(3):q2m_hi(3),NQ)
    real(rt), intent(in) :: q2p(q2p_lo(1):q2p_hi(1),q2p_lo(2):q2p_hi(2),q2p_lo(3):q2p_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: f1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NVAR)
    real(rt), intent(in) :: q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)

    real(rt), intent(out) :: q2mo(q2mo_lo(1):q2mo_hi(1),q2mo_lo(2):q2mo_hi(2),q2mo_lo(3):q2mo_hi(3),NQ)
    real(rt), intent(out) :: q2po(q2po_lo(1):q2po_hi(1),q2po_lo(2):q2po_hi(2),q2po_lo(3):q2po_hi(3),NQ)
#if AMREX_SPACEDIM == 2
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
#endif

    integer :: d, il, jl, kl, ir, jr, kr

    integer :: i, j, k, n, nqp, ipassive

    real(rt) :: lq2(NQ), lq2o(NQ)

    real(rt) :: rhoinv
    real(rt) :: rrnew
    real(rt) :: rrr2, rrl2
    real(rt) :: rur2, rul2
    real(rt) :: rvr2, rvl2
    real(rt) :: rwr2, rwl2
    real(rt) :: ekenr2, ekenl2
    real(rt) :: rer2, rel2
    real(rt) :: rrnewr2, rrnewl2
    real(rt) :: runewr2, runewl2
    real(rt) :: rvnewr2, rvnewl2
    real(rt) :: rwnewr2, rwnewl2
    real(rt) :: renewr2, renewl2
    real(rt) :: pnewr2, pnewl2
    real(rt) :: rhoekenr2, rhoekenl2
    real(rt) :: pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt) :: compu
    real(rt) :: gamc

#ifdef RADIATION
    real(rt) :: dre, dmom, divu
    real(rt), dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
                                        lamge, luge, der
    real(rt) :: eddf, f1, ugc
    integer  :: g
#endif

    logical :: reset_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------------
             ! update all of the passively-advected quantities with the
             ! transverse term and convert back to the primitive quantity
             !-------------------------------------------------------------------------

             !       qm|qp
             !         |
             ! --------+--------
             !   i-1       i
             !        i-1/2
             !
             ! the qm state will see the transverse flux in zone i-1

             ! Loop over plus and minus states
             
             do d = -1, 0

                ! We are handling the states at the interface of
                ! (i, i+1) in the x-direction, and similarly for
                ! the y- and z- directions.
                
                il = i
                jl = j
                kl = k
 
                if (idir1 == 1) then
                   ir = i+1
                   jr = j
                   kr = k
                else if (idir1 == 2) then
                   ir = i
                   jr = j+1
                   kr = k
                else
                   ir = i
                   jr = j
                   kr = k+1
                end if

                ! We're handling both the plus and minus states;
                ! for the minus state we're shifting one zone to
                ! the left in our chosen direction.

                if (idir2 == 1) then
                   il = il+d
                   ir = ir+d
                else if (idir2 == 2) then
                   jl = jl+d
                   jr = jr+d
                else
                   kl = kl+d
                   kr = kr+d
                end if

                if (d == -1) then
                   lq2(:) = q2m(i,j,k,:)
                else
                   lq2(:) = q2p(i,j,k,:)
                end if

                do ipassive = 1, npassive
                   n  = upass_map(ipassive)
                   nqp = qpass_map(ipassive)

#if AMREX_SPACEDIM == 2
                   rrnew = lq2(QRHO) - hdt*(area1(ir,jr,kr)*f1(ir,jr,kr,URHO) - &
                                            area1(il,jl,kl)*f1(il,jl,kl,URHO)) / vol(il,jl,kl)
                   compu = lq2(QRHO)*lq2(nqp) - &
                           hdt*(area1(ir,jr,kr)*f1(ir,jr,kr,n) - &
                                area1(il,jl,kl)*f1(il,jl,kl,n)) / vol(il,jl,kl)
#else
                   rrnew = lq2(QRHO) - cdtdx*(f1(ir,jr,kr,URHO) - f1(il,jl,kl,URHO))
                   compu = lq2(QRHO)*lq2(nqp) - cdtdx*(f1(ir,jr,kr,n) - f1(il,jl,kl,n))
#endif
                   lq2o(nqp) = compu/rrnew

                end do

                !-------------------------------------------------------------------
                ! add the transverse flux difference in the 1-direction to 2-states
                ! for the fluid variables
                !-------------------------------------------------------------------

                pgp  = q1(ir,jr,kr,GDPRES)
                pgm  = q1(il,jl,kl,GDPRES)
                ugp  = q1(ir,jr,kr,GDU+idir1-1)
                ugm  = q1(il,jl,kl,GDU+idir1-1)
                gegp = q1(ir,jr,kr,GDGAME)
                gegm = q1(il,jl,kl,GDGAME)

#ifdef RADIATION
                lambda = qaux(il,jl,kl,QLAMS:QLAMS+ngroups-1)
                ugc = HALF*(ugp+ugm)
                ergp = q1(ir,jr,kr,GDERADS:GDERADS-1+ngroups)
                ergm = q1(il,jl,kl,GDERADS:GDERADS-1+ngroups)
#endif

                ! we need to augment our conserved system with either a p
                ! equation or gammae (if we have ppm_predict_gammae = 1) to
                ! be able to deal with the general EOS

#if AMREX_SPACEDIM == 2
                dup = area1(ir,jr,kr)*pgp*ugp - area1(il,jl,kl)*pgm*ugm
                du = area1(ir,jr,kr)*ugp-area1(il,jl,kl)*ugm
#else
                dup = pgp*ugp - pgm*ugm
                du = ugp-ugm
#endif
                pav = HALF*(pgp+pgm)
                uav = HALF*(ugp+ugm)
                geav = HALF*(gegp+gegm)
                dge = gegp-gegm

                ! this is the gas gamma_1
#ifdef RADIATION
                gamc = qaux(il,jl,kl,QGAMCG)
#else
                gamc = qaux(il,jl,kl,QGAMC)
#endif

#ifdef RADIATION
                lamge = lambda(:) * (ergp(:)-ergm(:))
                dmom = - cdtdx*sum(lamge(:))
                luge = ugc * lamge(:)
                dre = -cdtdx*sum(luge)

                if (fspace_type .eq. 1 .and. comoving) then
                   do g=0, ngroups-1
                      eddf = Edd_factor(lambda(g))
                      f1 = HALF*(ONE-eddf)
                      der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
                   end do
                else if (fspace_type .eq. 2) then
#if AMREX_SPACEDIM == 2
                   divu = (area1(ir,jr,kr)*ugp-area1(il,jl,kl)*ugm)/vol(il,jl,kl)
                   do g=0, ngroups-1
                      eddf = Edd_factor(lambda(g))
                      f1 = 0.5e0_rt*(1.e0_rt-eddf)
                      der(g) = -hdt * f1 * 0.5e0_rt*(ergp(g)+ergm(g)) * divu
                   end do
#else
                   do g=0, ngroups-1
                      eddf = Edd_factor(lambda(g))
                      f1 = HALF*(ONE-eddf)
                      der(g) = cdtdx * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
                   end do
#endif
                else ! mixed frame
                   der(:) = cdtdx * luge
                end if
#endif

                ! Convert to conservation form
                rrl2 = lq2(QRHO)
                rul2 = rrl2*lq2(QU)
                rvl2 = rrl2*lq2(QV)
                rwl2 = rrl2*lq2(QW)
                ekenl2 = HALF*rrl2*sum(lq2(QU:QW)**2)
                rel2 = lq2(QREINT) + ekenl2
#ifdef RADIATION
                erl  = lq2(qrad:qradhi)
#endif

#if AMREX_SPACEDIM == 2
                rrnewl2 = rrl2 - hdt*(area1(ir,jr,kr)*f1(ir,jr,kr,URHO) -  &
                                      area1(il,jl,kl)*f1(il,jl,kl,URHO))/vol(il,jl,kl)
                runewl2 = rul2 - hdt*(area1(ir,jr,kr)*f1(ir,jr,kr,UMX) -  &
                                      area1(il,jl,kl)*f1(il,jl,kl,UMX))/vol(il,jl,kl)
                if (.not. mom_flux_has_p(1)%comp(UMX)) then
                   runewl2 = runewl2 - cdtdx *(pgp-pgm)
                endif
                rvnewl2 = rvl2 - hdt*(area1(ir,jr,kr)*f1(ir,jr,kr,UMY) -  &
                                      area1(il,jl,kl)*f1(il,jl,kl,UMY))/vol(il,jl,kl)
                rwnewl2 = rwl2 - hdt*(area1(ir,jr,kr)*f1(ir,jr,kr,UMZ) -  &
                                      area1(il,jl,kl)*f1(il,jl,kl,UMZ))/vol(il,jl,kl)
                renewl2 = rel2 - hdt*(area1(ir,jr,kr)*f1(ir,jr,kr,UEDEN) - &
                                      area1(il,jl,kl)*f1(il,jl,kl,UEDEN))/vol(il,jl,kl)

#ifdef RADIATION
                runewl2 = runewl2 - HALF*hdt*(area1(ir,jr,kr)+area1(il,jl,kl))*sum(lamge)/vol(il,jl,kl)
                renewl2 = renewl2 + dre
                ernewl(:) = erl(:) - hdt*(area1(ir,jr,kr)*rf1(ir,jr,kr,:) - &
                                          area1(il,jl,kl)*rf1(il,jl,kl,:))/vol(il,jl,kl) + der(:)
#endif

#else
                ! Add transverse predictor
                rrnewl2 = rrl2 - cdtdx*(f1(ir,jr,kr,URHO) - f1(il,jl,kl,URHO))
                runewl2 = rul2 - cdtdx*(f1(ir,jr,kr,UMX) - f1(il,jl,kl,UMX))
                rvnewl2 = rvl2 - cdtdx*(f1(ir,jr,kr,UMY) - f1(il,jl,kl,UMY))
                rwnewl2 = rwl2 - cdtdx*(f1(ir,jr,kr,UMZ) - f1(il,jl,kl,UMZ))
                renewl2 = rel2 - cdtdx*(f1(ir,jr,kr,UEDEN) - f1(il,jl,kl,UEDEN))
#ifdef RADIATION
                runewl2 = runewl2 + dmom
                renewl2 = renewl2 + dre
                ernewl  = erl(:) - cdtdx*(rf1(ir,jr,kr,:) - rf1(il,jl,kl,:)) + der(:)
#endif
#endif
                ! Reset to original value if adding transverse terms made density negative
                reset_state = .false.
                if (transverse_reset_density == 1 .and. rrnewl2 < ZERO) then
                   rrnewl2 = rrl2
                   runewl2 = rul2
                   rvnewl2 = rvl2
                   rwnewl2 = rwl2
                   renewl2 = rel2
#ifdef RADIATION
                   ernewl  = erl(:)
#endif
                   reset_state = .true.
                endif

                ! Convert back to primitive form
                lq2o(QRHO) = rrnewl2
                rhoinv = ONE/rrnewl2
                lq2o(QU) = runewl2*rhoinv
                lq2o(QV) = rvnewl2*rhoinv
                lq2o(QW) = rwnewl2*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenl2 = HALF*(runewl2**2 + rvnewl2**2 + rwnewl2**2)*rhoinv
                lq2o(QREINT) = renewl2 - rhoekenl2

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. lq2o(QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
#if AMREX_SPACEDIM == 2
                      lq2o(QREINT) = lq2(QREINT) - &
                                     hdt*(area1(ir,jr,kr)*f1(ir,jr,kr,UEINT) - &
                                          area1(il,jl,kl)*f1(il,jl,kl,UEINT) + pav*du)/vol(il,jl,kl)
#else
                      lq2o(QREINT) = lq2(QREINT) - &
                                     cdtdx*(f1(ir,jr,kr,UEINT) - f1(il,jl,kl,UEINT) + pav*du)
#endif
                   end if

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
#if AMREX_SPACEDIM == 2
                      pnewl2 = lq2(QPRES) - hdt*(dup + pav*du*(gamc - ONE))/vol(il,jl,kl)
#else
                      pnewl2 = lq2(QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
#endif
                      lq2o(QPRES) = max(pnewl2,small_pres)
                   else
                      ! Update gammae with its transverse terms
#if AMREX_SPACEDIM == 2
                      lq2o(QGAME) = lq2(QGAME) + &
                                    hdt*( (geav-ONE)*(geav - gamc)*du)/vol(il,jl,kl) - cdtdx*uav*dge
#else
                      lq2o(QGAME) = lq2(QGAME) + &
                                    cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )
#endif

                      ! and compute the p edge state from this and (rho e)
                      lq2o(QPRES) = lq2o(QREINT)*(lq2o(QGAME)-ONE)
                      lq2o(QPRES) = max(lq2o(QPRES), small_pres)
                   end if
                else
                   lq2o(QPRES) = lq2(QPRES)
                   lq2o(QGAME) = lq2(QGAME)
                endif

#ifdef RADIATION
                lq2o(qrad:qradhi) = ernewl(:)
                lq2o(qptot  ) = sum(lambda(:)*ernewl(:)) + lq2o(QPRES)
                lq2o(qreitot) = sum(lq2o(qrad:qradhi)) + lq2o(QREINT)
#endif

                if (d == -1) then
                   q2mo(i,j,k,:) = lq2o(:)
                   call reset_edge_state_thermo(q2mo, q2mo_lo, q2mo_hi, i, j, k)
                else
                   q2po(i,j,k,:) = lq2o(:)
                   call reset_edge_state_thermo(q2po, q2po_lo, q2po_hi, i, j, k)
                end if

             end do

          end do
       end do
    end do

  end subroutine trans1_on_2states



#if AMREX_SPACEDIM == 3
  !===========================================================================
  ! transyz
  !===========================================================================
  subroutine transyz(lo, hi, &
                     idir, &
                     qm, qm_lo, qm_hi, &
                     qmo, qmo_lo, qmo_hi, &
                     qp, qp_lo, qp_hi, &
                     qpo, qpo_lo, qpo_hi, &
                     qaux, qa_lo, qa_hi, &
                     fyz, fyz_lo, fyz_hi, &
#ifdef RADIATION
                     rfyz, rfyz_lo, rfyz_hi, &
#endif
                     fzy, fzy_lo, fzy_hi, &
#ifdef RADIATION
                     rfzy, rfzy_lo, rfzy_hi, &
#endif
                     qy, qy_lo, qy_hi, &
                     qz, qz_lo, qz_hi, &
                     hdt, cdtdy, cdtdz) bind(C, name="transyz")

    ! here, lo and hi are the bounds of the x interfaces we are looping over

    use amrex_constants_module, only : ZERO, ONE, HALF
    use network, only : nspec, naux
    use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
                                   fspace_type, comoving, &
                                   GDERADS, GDLAMS, &
                                   QCG, QGAMCG, QLAMS, &
#endif
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                   small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_predict_gammae

#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t


    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qmo_lo(3), qmo_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qpo_lo(3), qpo_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fyz_lo(3), fyz_hi(3)
    integer, intent(in) :: fzy_lo(3), fzy_hi(3)
    integer, intent(in) :: qy_lo(3), qy_hi(3)
    integer, intent(in) :: qz_lo(3), qz_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in), value :: idir

    real(rt), intent(in), value :: hdt, cdtdy, cdtdz

#ifdef RADIATION
    integer, intent(in) :: rfyz_lo(3), rfyz_hi(3)
    integer, intent(in) :: rfzy_lo(3), rfzy_hi(3)
    real(rt), intent(in) :: rfyz(fyz_lo(1):fyz_hi(1),fyz_lo(2):fyz_hi(2),fyz_lo(3):fyz_hi(3),0:ngroups-1)
    real(rt), intent(in) :: rfzy(fzy_lo(1):fzy_hi(1),fzy_lo(2):fzy_hi(2),fzy_lo(3):fzy_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(in) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)
    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),NQ)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),NQ)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fyz(fyz_lo(1):fyz_hi(1),fyz_lo(2):fyz_hi(2),fyz_lo(3):fyz_hi(3),NVAR)
    real(rt), intent(in) :: fzy(fzy_lo(1):fzy_hi(1),fzy_lo(2):fzy_hi(2),fzy_lo(3):fzy_hi(3),NVAR)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

    integer :: i, j, k, n, nqp, ipassive
    integer :: d, idir1, idir2, il, jl, kl, il1, jl1, kl1, ir1, jr1, kr1, il2, jl2, kl2, ir2, jr2, kr2

    real(rt) :: lqo(NQ), lq(NQ)

    real(rt) :: rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt) :: rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt) :: rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt) :: rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt) :: pnewr, pnewl
    real(rt) :: pgyp, pgym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    real(rt) :: pgzp, pgzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    real(rt) :: uyav, geyav, dgey, uzav, gezav, dgez
    real(rt) :: compr, compl, compnr, compnl

#ifdef RADIATION
    real(rt) :: dmy, dmz, dre
    real(rt), dimension(0:ngroups-1) :: der, lambda, lugey, lugez, lgey, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergyp, ergzm, ergym
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !$gpu

    !-------------------------------------------------------------------
    ! add the transverse yz and zy differences to the x-states for the
    ! fluid variables
    !-------------------------------------------------------------------

    if (idir == 1) then
       idir1 = 2
       idir2 = 3
    else if (idir == 2) then
       idir1 = 1
       idir2 = 3
    else
       idir1 = 1
       idir2 = 2
    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------------
             ! update all of the passively-advected quantities with the
             ! transerse term and convert back to the primitive quantity
             !-------------------------------------------------------------------------

             do d = -1, 0

                il = i
                jl = j
                kl = k
             
                il1 = i
                jl1 = j
                kl1 = k

                il2 = i
                jl2 = j
                kl2 = k

                if (idir == 1) then
                   ir1 = i+d
                   jr1 = j+1
                   kr1 = k

                   ir2 = i+d
                   jr2 = j
                   kr2 = k+1

                   il  = i+d
                   il1 = i+d
                   il2 = i+d
                else if (idir == 2) then
                   ir1 = i+1
                   jr1 = j+d
                   kr1 = k

                   ir2 = i
                   jr2 = j+d
                   kr2 = k+1

                   jl  = j+d
                   jl1 = j+d
                   jl2 = j+d
                else
                   ir1 = i+1
                   jr1 = j
                   kr1 = k+d

                   ir2 = i
                   jr2 = j+1
                   kr2 = k+d

                   kl  = k+d
                   kl1 = k+d
                   kl2 = k+d
                end if

                if (d == -1) then
                   lq(:) = qm(i,j,k,:)
                else
                   lq(:) = qp(i,j,k,:)
                end if

                do ipassive = 1,npassive
                   n  = upass_map(ipassive)
                   nqp = qpass_map(ipassive)

                   rrr = lq(QRHO)
                   compr = rrr*lq(nqp)
                   rrnewr = rrr - cdtdy*(fyz(ir1,jr1,kr1,URHO) - fyz(il1,jl1,kl1,URHO)) &
                                - cdtdz*(fzy(ir2,jr2,kr2,URHO) - fzy(il2,jl2,kl2,URHO))
                   compnr = compr - cdtdy*(fyz(ir1,jr1,kr1,n) - fyz(il1,jl1,kl1,n)) &
                                  - cdtdz*(fzy(ir2,jr2,kr2,n) - fzy(il2,jl2,kl2,n))

                   lqo(nqp) = compnr/rrnewr
                end do

                pgyp  = qy(ir1,jr1,kr1,GDPRES)
                pgym  = qy(il1,jl1,kl1,GDPRES)
                ugyp  = qy(ir1,jr1,kr1,GDU+idir1-1)
                ugym  = qy(il1,jl1,kl1,GDU+idir1-1)
                gegyp = qy(ir1,jr1,kr1,GDGAME)
                gegym = qy(il1,jl1,kl1,GDGAME)
#ifdef RADIATION
                ergyp = qy(ir1,jr1,kr1,GDERADS:GDERADS-1+ngroups)
                ergym = qy(il1,jl1,kl1,GDERADS:GDERADS-1+ngroups)
#endif

                pgzp  = qz(ir2,jr2,kr2,GDPRES)
                pgzm  = qz(il2,jl2,kl2,GDPRES)
                ugzp  = qz(ir2,jr2,kr2,GDU+idir2-1)
                ugzm  = qz(il2,jl2,kl2,GDU+idir2-1)
                gegzp = qz(ir2,jr2,kr2,GDGAME)
                gegzm = qz(il2,jl2,kl2,GDGAME)
#ifdef RADIATION
                ergzp = qz(ir2,jr2,kr2,GDERADS:GDERADS-1+ngroups)
                ergzm = qz(il2,jl2,kl2,GDERADS:GDERADS-1+ngroups)
#endif

                duyp = pgyp*ugyp - pgym*ugym
                pyav = HALF*(pgyp+pgym)
                uyav = HALF*(ugyp+ugym)
                geyav = HALF*(gegyp+gegym)
                duy = ugyp-ugym
                dgey = gegyp-gegym
#ifdef RADIATION
                pynew = cdtdy*(duyp + pyav*duy*(qaux(il,jl,kl,QGAMCG) - ONE))
                geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(il,jl,kl,QGAMCG))*duy - uyav*dgey )
#else
                pynew = cdtdy*(duyp + pyav*duy*(qaux(il,jl,kl,QGAMC) - ONE))
                geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(il,jl,kl,QGAMC))*duy - uyav*dgey )
#endif

                duzp = pgzp*ugzp - pgzm*ugzm
                pzav = HALF*(pgzp+pgzm)
                uzav = HALF*(ugzp+ugzm)
                gezav = HALF*(gegzp+gegzm)
                duz = ugzp-ugzm
                dgez = gegzp-gegzm
#ifdef RADIATION
                pznew = cdtdz*(duzp + pzav*duz*(qaux(il,jl,kl,QGAMCG) - ONE))
                geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(il,jl,kl,QGAMCG))*duz - uzav*dgez )
#else
                pznew = cdtdz*(duzp + pzav*duz*(qaux(il,jl,kl,QGAMC) - ONE))
                geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(il,jl,kl,QGAMC))*duz - uzav*dgez )
#endif

#ifdef RADIATION
                lambda(:) = qaux(il,jl,kl,QLAMS:QLAMS+ngroups-1)

                lgey = lambda(:) * (ergyp(:)-ergym(:))
                lgez = lambda(:) * (ergzp(:)-ergzm(:))
                dmy = - cdtdy*sum(lgey)
                dmz = - cdtdz*sum(lgez)
                lugey = HALF*(ugyp+ugym) * lgey(:)
                lugez = HALF*(ugzp+ugzm) * lgez(:)
                dre = -cdtdy*sum(lugey) - cdtdz*sum(lugez)

                if (fspace_type .eq. 1 .and. comoving) then
                   do g=0, ngroups-1
                      eddf = Edd_factor(lambda(g))
                      f1 = HALF*(ONE-eddf)
                      der(g) = f1*(cdtdy*HALF*(ugyp+ugym)*(ergyp(g)-ergym(g)) &
                           +       cdtdz*HALF*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
                   end do
                else if (fspace_type .eq. 2) then
                   do g=0, ngroups-1
                      eddf = Edd_factor(lambda(g))
                      f1 = HALF*(ONE-eddf)
                      der(g) = f1*(cdtdy*HALF*(ergyp(g)+ergym(g))*(ugym-ugyp) &
                           +       cdtdz*HALF*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
                   end do
                else ! mixed frame
                   der(:) = cdtdy*lugey + cdtdz*lugez
                end if
#endif

                ! Convert to conservation form
                rrr = lq(QRHO)
                rur = rrr*lq(QU)
                rvr = rrr*lq(QV)
                rwr = rrr*lq(QW)
                ekenr = HALF*rrr*sum(lq(QU:QW)**2)
                rer = lq(QREINT) + ekenr
#ifdef RADIATION
                err = lq(qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewr = rrr - cdtdy*(fyz(ir1,jr1,kr1,URHO) - fyz(il1,jl1,kl1,URHO)) &
                             - cdtdz*(fzy(ir2,jr2,kr2,URHO) - fzy(il2,jl2,kl2,URHO))
                runewr = rur - cdtdy*(fyz(ir1,jr1,kr1,UMX) - fyz(il1,jl1,kl1,UMX)) &
                             - cdtdz*(fzy(ir2,jr2,kr2,UMX) - fzy(il2,jl2,kl2,UMX))
                rvnewr = rvr - cdtdy*(fyz(ir1,jr1,kr1,UMY) - fyz(il1,jl1,kl1,UMY)) &
                             - cdtdz*(fzy(ir2,jr2,kr2,UMY) - fzy(il2,jl2,kl2,UMY))
                rwnewr = rwr - cdtdy*(fyz(ir1,jr1,kr1,UMZ) - fyz(il1,jl1,kl1,UMZ)) &
                             - cdtdz*(fzy(ir2,jr2,kr2,UMZ) - fzy(il2,jl2,kl2,UMZ))
                renewr = rer - cdtdy*(fyz(ir1,jr1,kr1,UEDEN) - fyz(il1,jl1,kl1,UEDEN)) &
                             - cdtdz*(fzy(ir2,jr2,kr2,UEDEN) - fzy(il2,jl2,kl2,UEDEN))
#ifdef RADIATION
                rvnewr = rvnewr + dmy
                rwnewr = rwnewr + dmz
                renewr = renewr + dre
                ernewr = err(:) - cdtdy*(rfyz(ir1,jr1,kr1,:) - rfyz(il1,jl1,kl1,:)) &
                                - cdtdz*(rfzy(ir2,jr2,kr2,:) - rfzy(il2,jl2,kl2,:)) &
                                + der(:)
#endif

                ! Reset to original value if adding transverse terms
                ! made density negative
                reset_state = .false.
                if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
                   rrnewr = rrr
                   runewr = rur
                   rvnewr = rvr
                   rwnewr = rwr
                   renewr = rer
#ifdef RADIATION
                   ernewr = err(:)
#endif
                   reset_state = .true.
                end if

                lqo(QRHO  ) = rrnewr
                lqo(QU    ) = runewr/rrnewr
                lqo(QV    ) = rvnewr/rrnewr
                lqo(QW    ) = rwnewr/rrnewr

                ! note: we run the risk of (rho e) being negative here
                rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
                lqo(QREINT) = renewr - rhoekenr

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. lqo(QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating
                      ! (rho e).
                      lqo(QREINT) = lq(QREINT) - &
                                    cdtdy*(fyz(ir1,jr1,kr1,UEINT) - fyz(il1,jl1,kl1,UEINT) + pyav*duy) - &
                                    cdtdz*(fzy(ir2,jr2,kr2,UEINT) - fzy(il2,jl2,kl2,UEINT) + pzav*duz)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewr = lq(QPRES) - pynew - pznew
                      lqo(QPRES) = pnewr
                   else
                      ! Update gammae with its transverse terms
                      lqo(QGAME) = lq(QGAME) + geynew + geznew

                      ! and compute the p edge state from this and (rho e)
                      lqo(QPRES) = lqo(QREINT)*(lqo(QGAME)-ONE)
                   end if
                else
                   lqo(QPRES) = lq(QPRES)
                   lqo(QGAME) = lq(QGAME)
                endif

                lqo(QPRES) = max(lqo(QPRES), small_pres)

#ifdef RADIATION
                lqo(qrad:qradhi) = ernewr(:)
                lqo(qptot  ) = sum(lambda(:)*ernewr(:)) + lqo(QPRES)
                lqo(qreitot) = sum(lqo(qrad:qradhi)) + lqo(QREINT)
#endif

                if (d == -1) then
                   qmo(i,j,k,:) = lqo(:)
                   call reset_edge_state_thermo(qmo, qmo_lo, qmo_hi, i, j, k)
                else
                   qpo(i,j,k,:) = lqo(:)
                   call reset_edge_state_thermo(qpo, qpo_lo, qpo_hi, i, j, k)
                end if

             end do

          end do
       end do
    end do

  end subroutine transyz

  !===========================================================================
  ! transxz
  !===========================================================================
  subroutine transxz(lo, hi, &
                     qm, qm_lo, qm_hi, &
                     qmo, qmo_lo, qmo_hi, &
                     qp, qp_lo, qp_hi, &
                     qpo, qpo_lo, qpo_hi, &
                     qaux, qa_lo, qa_hi, &
                     fxz, fxz_lo, fxz_hi, &
#ifdef RADIATION
                     rfxz, rfxz_lo, rfxz_hi, &
#endif
                     fzx, fzx_lo, fzx_hi, &
#ifdef RADIATION
                     rfzx, rfzx_lo, rfzx_hi, &
#endif
                     qx, qx_lo, qx_hi, &
                     qz, qz_lo, qz_hi, &
                     hdt, cdtdx, cdtdz) bind(C, name="transxz")

    ! here, lo and hi are the bounds of the y interfaces we are looping over

    use amrex_constants_module, only : ZERO, ONE, HALF
    use network, only : nspec, naux
    use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
                                   fspace_type, comoving, &
                                   GDERADS, GDLAMS, &
                                   QCG, QGAMCG, QLAMS, &
#endif
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                   small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_predict_gammae

#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t


    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qmo_lo(3), qmo_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qpo_lo(3), qpo_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fxz_lo(3), fxz_hi(3)
    integer, intent(in) :: fzx_lo(3), fzx_hi(3)
    integer, intent(in) :: qx_lo(3),qx_hi(3)
    integer, intent(in) :: qz_lo(3),qz_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in), value :: hdt, cdtdx, cdtdz

#ifdef RADIATION
    integer, intent(in) :: rfxz_lo(3), rfxz_hi(3)
    integer, intent(in) :: rfzx_lo(3), rfzx_hi(3)
    real(rt), intent(in) :: rfxz(rfxz_lo(1):rfxz_hi(1),rfxz_lo(2):rfxz_hi(2),rfxz_lo(3):rfxz_hi(3),0:ngroups-1)
    real(rt), intent(in) :: rfzx(rfzx_lo(1):rfzx_hi(1),rfzx_lo(2):rfzx_hi(2),rfzx_lo(3):rfzx_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(in) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)
    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),NQ)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),NQ)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fxz(fxz_lo(1):fxz_hi(1),fxz_lo(2):fxz_hi(2),fxz_lo(3):fxz_hi(3),NVAR)
    real(rt), intent(in) :: fzx(fzx_lo(1):fzx_hi(1),fzx_lo(2):fzx_hi(2),fzx_lo(3):fzx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

    integer :: i, j, k, n, nqp, ipassive

    real(rt) :: rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt) :: rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt) :: rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt) :: rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt) :: pnewr, pnewl
    real(rt) :: pgxp, pgxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    real(rt) :: pgzp, pgzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    real(rt) :: uxav, gexav, dgex, uzav, gezav, dgez
    real(rt) :: compr, compl, compnr, compnl

#ifdef RADIATION
    real(rt) :: dmx, dmz, dre
    real(rt), dimension(0:ngroups-1) :: der, lambda, lugex, lugez, lgex, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergxp, ergzm,  ergxm
    real(rt) :: eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !$gpu

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrr = qp(i,j,k,QRHO)
                compr = rrr*qp(i,j,k,nqp)
                rrnewr = rrr - cdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                             - cdtdz*(fzx(i  ,j,k+1,URHO) - fzx(i,j,k,URHO))
                compnr = compr - cdtdx*(fxz(i+1,j,k,n) - fxz(i,j,k,n)) &
                               - cdtdz*(fzx(i  ,j,k+1,n) - fzx(i,j,k,n))

                qpo(i,j  ,k,nqp) = compnr/rrnewr

                rrl = qm(i,j,k,QRHO)
                compl = rrl*qm(i,j,k,nqp)
                rrnewl = rrl - cdtdx*(fxz(i+1,j-1,k,URHO) - fxz(i,j-1,k,URHO)) &
                             - cdtdz*(fzx(i  ,j-1,k+1,URHO) - fzx(i,j-1,k,URHO))
                compnl = compl - cdtdx*(fxz(i+1,j-1,k,n) - fxz(i,j-1,k,n)) &
                               - cdtdz*(fzx(i  ,j-1,k+1,n) - fzx(i,j-1,k,n))

                qmo(i,j,k,nqp) = compnl/rrnewl
             end do
          end do
       end do
    end do

    !-------------------------------------------------------------------
    ! add the transverse xz and zx differences to the y-states for the
    ! fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qypo state
             !-------------------------------------------------------------------

             pgxp  = qx(i+1,j,k,GDPRES)
             pgxm  = qx(i,j,k,GDPRES)
             ugxp  = qx(i+1,j,k,GDU)
             ugxm  = qx(i,j,k,GDU)
             gegxp = qx(i+1,j,k,GDGAME)
             gegxm = qx(i,j,k,GDGAME)
#ifdef RADIATION
             ergxp = qx(i+1,j,k,GDERADS:GDERADS-1+ngroups)
             ergxm = qx(i  ,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             pgzp  = qz(i,j,k+1,GDPRES)
             pgzm  = qz(i,j,k,GDPRES)
             ugzp  = qz(i,j,k+1,GDW)
             ugzm  = qz(i,j,k,GDW)
             gegzp = qz(i,j,k+1,GDGAME)
             gegzm = qz(i,j,k,GDGAME)
#ifdef RADIATION
             ergzp = qz(i,j,k+1,GDERADS:GDERADS-1+ngroups)
             ergzm = qz(i,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm
#ifdef RADIATION
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMCG) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMCG))*dux - uxav*dgex )
#else
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMC))*dux - uxav*dgex )
#endif

             duzp = pgzp*ugzp - pgzm*ugzm
             pzav = HALF*(pgzp+pgzm)
             uzav = HALF*(ugzp+ugzm)
             gezav = HALF*(gegzp+gegzm)
             duz = ugzp-ugzm
             dgez = gegzp-gegzm
#ifdef RADIATION
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMCG) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMCG))*duz - uzav*dgez )
#else
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMC) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMC))*duz - uzav*dgez )
#endif

#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)

             lgex = lambda(:) * (ergxp(:)-ergxm(:))
             lgez = lambda(:) * (ergzp(:)-ergzm(:))
             dmx = - cdtdx*sum(lgex)
             dmz = - cdtdz*sum(lgez)
             lugex = HALF*(ugxp+ugxm) * lgex(:)
             lugez = HALF*(ugzp+ugzm) * lgez(:)
             dre = -cdtdx * sum(lugex) - cdtdz * sum(lugez)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = f1*(cdtdx*HALF*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                        +       cdtdz*HALF*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
                end do
             else if (fspace_type .eq. 2) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = f1*(cdtdx*HALF*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                        +       cdtdz*HALF*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
                end do
             else ! mixed frame
                der(:) = cdtdx*lugex + cdtdz*lugez
             end if
#endif

             ! Convert to conservation form
             rrr = qp(i,j,k,QRHO)
             rur = rrr*qp(i,j,k,QU)
             rvr = rrr*qp(i,j,k,QV)
             rwr = rrr*qp(i,j,k,QW)
             ekenr = HALF*rrr*sum(qp(i,j,k,QU:QW)**2)
             rer = qp(i,j,k,QREINT) + ekenr
#ifdef RADIATION
             err = qp(i,j,k,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewr = rrr - cdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                          - cdtdz*(fzx(i,j,k+1,URHO) - fzx(i,j,k,URHO))
             runewr = rur - cdtdx*(fxz(i+1,j,k,UMX) - fxz(i,j,k,UMX)) &
                          - cdtdz*(fzx(i,j,k+1,UMX) - fzx(i,j,k,UMX))
             rvnewr = rvr - cdtdx*(fxz(i+1,j,k,UMY) - fxz(i,j,k,UMY)) &
                          - cdtdz*(fzx(i,j,k+1,UMY) - fzx(i,j,k,UMY))
             rwnewr = rwr - cdtdx*(fxz(i+1,j,k,UMZ) - fxz(i,j,k,UMZ)) &
                          - cdtdz*(fzx(i,j,k+1,UMZ) - fzx(i,j,k,UMZ))
             renewr = rer - cdtdx*(fxz(i+1,j,k,UEDEN) - fxz(i,j,k,UEDEN)) &
                          - cdtdz*(fzx(i,j,k+1,UEDEN) - fzx(i,j,k,UEDEN))
#ifdef RADIATION
             runewr = runewr + dmx
             rwnewr = rwnewr + dmz
             renewr = renewr + dre
             ernewr = err(:) - cdtdx*(rfxz(i+1,j,k,:) - rfxz(i,j,k,:)) &
                             - cdtdz*(rfzx(i  ,j,k+1,:) - rfzx(i,j,k,:)) &
                             + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                rwnewr = rwr
                renewr = rer
#ifdef RADIATION
                ernewr = err(:)
#endif
                reset_state = .true.
             endif

             ! Convert back to primitive form
             qpo(i,j,k,QRHO  ) = rrnewr
             qpo(i,j,k,QU    ) = runewr/rrnewr
             qpo(i,j,k,QV    ) = rvnewr/rrnewr
             qpo(i,j,k,QW    ) = rwnewr/rrnewr

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,k,QREINT) = renewr - rhoekenr

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating
                   ! (rho e).
                   qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                        - cdtdx*(fxz(i+1,j,k,UEINT) - fxz(i,j,k,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i  ,j,k+1,UEINT) - fzx(i,j,k,UEINT) + pzav*duz)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,k,QPRES) - pxnew - pznew
                   qpo(i,j,k,QPRES) = pnewr
                else
                   ! Update gammae with its transverse terms
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + gexnew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                endif
             else
                qpo(i,j,k,QPRES) = qp(i,j,k,QPRES)
                qpo(i,j,k,QGAME) = qp(i,j,k,QGAME)
             endif

             qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qpo, qpo_lo, qpo_hi, i, j, k)

#ifdef RADIATION
             qpo(i,j,k,qrad:qradhi) = ernewr(:)
             qpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k,QPRES)
             qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qradhi)) + qpo(i,j,k,QREINT)
#endif

             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------

             pgxp  = qx(i+1,j-1,k,GDPRES)
             pgxm  = qx(i,j-1,k,GDPRES)
             ugxp  = qx(i+1,j-1,k,GDU)
             ugxm  = qx(i,j-1,k,GDU)
             gegxp = qx(i+1,j-1,k,GDGAME)
             gegxm = qx(i,j-1,k,GDGAME)
#ifdef RADIATION
             ergxp = qx(i+1,j-1,k,GDERADS:GDERADS-1+ngroups)
             ergxm = qx(i  ,j-1,k,GDERADS:GDERADS-1+ngroups)
#endif

             pgzp  = qz(i,j-1,k+1,GDPRES)
             pgzm  = qz(i,j-1,k,GDPRES)
             ugzp  = qz(i,j-1,k+1,GDW)
             ugzm  = qz(i,j-1,k,GDW)
             gegzp = qz(i,j-1,k+1,GDGAME)
             gegzm = qz(i,j-1,k,GDGAME)
#ifdef RADIATION
             ergzp = qz(i,j-1,k+1,GDERADS:GDERADS-1+ngroups)
             ergzm = qz(i,j-1,k,GDERADS:GDERADS-1+ngroups)
#endif

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm
#ifdef RADIATION
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j-1,k,QGAMCG) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j-1,k,QGAMCG))*dux - uxav*dgex )
#else
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j-1,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j-1,k,QGAMC))*dux - uxav*dgex )
#endif

             duzp = pgzp*ugzp - pgzm*ugzm
             pzav = HALF*(pgzp+pgzm)
             uzav = HALF*(ugzp+ugzm)
             gezav = HALF*(gegzp+gegzm)
             duz = ugzp-ugzm
             dgez = gegzp-gegzm
#ifdef RADIATION
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j-1,k,QGAMCG) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j-1,k,QGAMCG))*duz - uzav*dgez )
#else
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j-1,k,QGAMC) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j-1,k,QGAMC))*duz - uzav*dgez )
#endif

#ifdef RADIATION
             lambda(:) = qaux(i,j-1,k,QLAMS:QLAMS+ngroups-1)

             lgex = lambda(:) * (ergxp(:)-ergxm(:))
             lgez = lambda(:) * (ergzp(:)-ergzm(:))
             dmx = - cdtdx*sum(lgex)
             dmz = - cdtdz*sum(lgez)
             lugex = HALF*(ugxp+ugxm) * lgex(:)
             lugez = HALF*(ugzp+ugzm) * lgez(:)
             dre = -cdtdx * sum(lugex) - cdtdz * sum(lugez)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = f1*(cdtdx*HALF*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                        +       cdtdz*HALF*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
                end do
             else if (fspace_type .eq. 2) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = f1*(cdtdx*HALF*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                        +       cdtdz*HALF*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
                end do
             else ! mixed frame
                der(:) = cdtdx*lugex + cdtdz*lugez
             end if
#endif

             ! Convert to conservation form
             rrl = qm(i,j,k,QRHO)
             rul = rrl*qm(i,j,k,QU)
             rvl = rrl*qm(i,j,k,QV)
             rwl = rrl*qm(i,j,k,QW)
             ekenl = HALF*rrl*sum(qm(i,j,k,QU:QW)**2)
             rel = qm(i,j,k,QREINT) + ekenl
#ifdef RADIATION
             erl = qm(i,j,k,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewl = rrl - cdtdx*(fxz(i+1,j-1,k,URHO) - fxz(i,j-1,k,URHO)) &
                          - cdtdz*(fzx(i,j-1,k+1,URHO) - fzx(i,j-1,k,URHO))
             runewl = rul - cdtdx*(fxz(i+1,j-1,k,UMX) - fxz(i,j-1,k,UMX)) &
                          - cdtdz*(fzx(i,j-1,k+1,UMX) - fzx(i,j-1,k,UMX))
             rvnewl = rvl - cdtdx*(fxz(i+1,j-1,k,UMY) - fxz(i,j-1,k,UMY)) &
                          - cdtdz*(fzx(i,j-1,k+1,UMY) - fzx(i,j-1,k,UMY))
             rwnewl = rwl - cdtdx*(fxz(i+1,j-1,k,UMZ) - fxz(i,j-1,k,UMZ)) &
                          - cdtdz*(fzx(i,j-1,k+1,UMZ) - fzx(i,j-1,k,UMZ))
             renewl = rel - cdtdx*(fxz(i+1,j-1,k,UEDEN) - fxz(i,j-1,k,UEDEN)) &
                          - cdtdz*(fzx(i,j-1,k+1,UEDEN) - fzx(i,j-1,k,UEDEN))
#ifdef RADIATION
             runewl = runewl + dmx
             rwnewl = rwnewl + dmz
             renewl = renewl + dre
             ernewl = erl(:) - cdtdx*(rfxz(i+1,j-1,k,:) - rfxz(i,j-1,k,:)) &
                             - cdtdz*(rfzx(i  ,j-1,k+1,:) - rfzx(i,j-1,k,:)) &
                             + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewl < ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                rwnewl = rwl
                renewl = rel
#ifdef RADIATION
                ernewl = erl(:)
#endif
                reset_state = .true.
             endif

             ! Convert back to primitive form
             qmo(i,j,k,QRHO  ) = rrnewl
             qmo(i,j,k,QU    ) = runewl/rrnewl
             qmo(i,j,k,QV    ) = rvnewl/rrnewl
             qmo(i,j,k,QW    ) = rwnewl/rrnewl

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i,j,k,QREINT) = renewl - rhoekenl

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo(i,j,k,QREINT) = qm(i,j,k,QREINT) &
                        - cdtdx*(fxz(i+1,j-1,k,UEINT) - fxz(i,j-1,k,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i,j-1,k+1,UEINT) - fzx(i,j-1,k,UEINT) + pzav*duz)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i,j,k,QPRES) - pxnew - pznew
                   qmo(i,j,k,QPRES) = pnewl
                else
                   ! Update gammae with its transverse terms
                   qmo(i,j,k,QGAME) = qm(i,j,k,QGAME) + gexnew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j,k,QPRES) = qmo(i,j,k,QREINT)*(qmo(i,j,k,QGAME)-ONE)
                endif
             else
                qmo(i,j,k,QPRES) = qm(i,j,k,QPRES)
                qmo(i,j,k,QGAME) = qm(i,j,k,QGAME)
             endif

             qmo(i,j,k,QPRES) = max(qmo(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qmo, qmo_lo, qmo_hi, i, j, k)

#ifdef RADIATION
             qmo(i,j,k,qrad:qradhi) = ernewl(:)
             qmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j,k,QPRES)
             qmo(i,j,k,qreitot) = sum(qmo(i,j,k,qrad:qradhi)) + qmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transxz

  !===========================================================================
  ! transxy
  !===========================================================================
  subroutine transxy(lo, hi, &
                     qm, qm_lo, qm_hi, &
                     qmo, qmo_lo, qmo_hi, &
                     qp, qp_lo, qp_hi, &
                     qpo, qpo_lo, qpo_hi, &
                     qaux, qa_lo, qa_hi, &
                     fxy, fxy_lo, fxy_hi, &
#ifdef RADIATION
                     rfxy, rfxy_lo, rfxy_hi, &
#endif
                     fyx, fyx_lo, fyx_hi, &
#ifdef RADIATION
                     rfyx, rfyx_lo, rfyx_hi, &
#endif
                     qx, qx_lo, qx_hi, &
                     qy, qy_lo, qy_hi, &
                     hdt, cdtdx, cdtdy) bind(C, name="transxy")

    ! here, lo and hi are the bounds of edges we are looping over


    use amrex_constants_module, only : ZERO, ONE, HALF
    use network, only : nspec, naux
    use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
                                   fspace_type, comoving, &
                                   GDERADS, GDLAMS, &
                                   QCG, QGAMCG, QLAMS, &
#endif
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                   small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_predict_gammae

#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t


    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qmo_lo(3), qmo_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qpo_lo(3), qpo_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fxy_lo(3), fxy_hi(3)
    integer, intent(in) :: fyx_lo(3), fyx_hi(3)
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: qy_lo(3), qy_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in), value :: hdt, cdtdx, cdtdy

#ifdef RADIATION
    integer, intent(in) :: rfxy_lo(3), rfxy_hi(3)
    integer, intent(in) :: rfyx_lo(3), rfyx_hi(3)
    real(rt), intent(in) :: rfxy(rfxy_lo(1):rfxy_hi(1),rfxy_lo(2):rfxy_hi(2),rfxy_lo(3):rfxy_hi(3),0:ngroups-1)
    real(rt), intent(in) :: rfyx(rfyx_lo(1):rfyx_hi(1),rfyx_lo(2):rfyx_hi(2),rfyx_lo(3):rfyx_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) ::   qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(in) ::   qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)
    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),NQ)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),NQ)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fxy(fxy_lo(1):fxy_hi(1),fxy_lo(2):fxy_hi(2),fxy_lo(3):fxy_hi(3),NVAR)
    real(rt), intent(in) :: fyx(fyx_lo(1):fyx_hi(1),fyx_lo(2):fyx_hi(2),fyx_lo(3):fyx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)

    integer :: i, j, k, n, nqp, ipassive

    real(rt) :: rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt) :: rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt) :: rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt) :: rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt) :: pnewr, pnewl
    real(rt) :: pgxp, pgxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    real(rt) :: pgyp, pgym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    real(rt) :: uxav, gexav, dgex, uyav, geyav, dgey
    real(rt) :: compr, compl, compnr, compnl

#ifdef RADIATION
    real(rt) :: dmx, dmy, dre
    real(rt), dimension(0:ngroups-1) :: der, lambda, lugex, lugey, lgex, lgey, &
         err, ernewr, erl, ernewl, ergxp, ergyp, ergxm, ergym
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !$gpu

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrr = qp(i,j,k,QRHO)
                compr = rrr*qp(i,j,k,nqp)
                rrnewr = rrr - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                             - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
                compnr = compr - cdtdx*(fxy(i+1,j,k,n) - fxy(i,j,k,n)) &
                               - cdtdy*(fyx(i,j+1,k,n) - fyx(i,j,k,n))

                qpo(i,j,k,nqp) = compnr/rrnewr

                rrl = qm(i,j,k,QRHO)
                compl = rrl*qm(i,j,k,nqp)
                rrnewl = rrl - cdtdx*(fxy(i+1,j,k-1,URHO) - fxy(i,j,k-1,URHO)) &
                             - cdtdy*(fyx(i,j+1,k-1,URHO) - fyx(i,j,k-1,URHO))
                compnl = compl - cdtdx*(fxy(i+1,j,k-1,n) - fxy(i,j,k-1,n)) &
                               - cdtdy*(fyx(i,j+1,k-1,n) - fyx(i,j,k-1,n))

                qmo(i,j,k,nqp) = compnl/rrnewl
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse xy and yx differences to the z-states for the
    ! fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qzpo state
             !-------------------------------------------------------------------

             pgxp = qx(i+1,j,k,GDPRES)
             pgxm = qx(i,j,k,GDPRES)
             ugxp = qx(i+1,j,k,GDU)
             ugxm = qx(i,j,k,GDU)
             gegxp = qx(i+1,j,k,GDGAME)
             gegxm = qx(i,j,k,GDGAME)
#ifdef RADIATION
             ergxp = qx(i+1,j,k,GDERADS:GDERADS-1+ngroups)
             ergxm = qx(i  ,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             pgyp = qy(i,j+1,k,GDPRES)
             pgym = qy(i,j,k,GDPRES)
             ugyp = qy(i,j+1,k,GDV)
             ugym = qy(i,j,k,GDV)
             gegyp = qy(i,j+1,k,GDGAME)
             gegym = qy(i,j,k,GDGAME)
#ifdef RADIATION
             ergyp = qy(i,j+1,k,GDERADS:GDERADS-1+ngroups)
             ergym = qy(i,j  ,k,GDERADS:GDERADS-1+ngroups)
#endif

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm
#ifdef RADIATION
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMCG) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMCG))*dux - uxav*dgex )
#else
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMC))*dux - uxav*dgex )
#endif

             duyp = pgyp*ugyp - pgym*ugym
             pyav = HALF*(pgyp+pgym)
             uyav = HALF*(ugyp+ugym)
             geyav = HALF*(gegyp+gegym)
             duy = ugyp-ugym
             dgey = gegyp-gegym
#ifdef RADIATION
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMCG) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMCG))*duy - uyav*dgey )
#else
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMC) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMC))*duy - uyav*dgey )
#endif

#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)

             lgex = lambda(:) * (ergxp(:)-ergxm(:))
             lgey = lambda(:) * (ergyp(:)-ergym(:))
             dmx = - cdtdx*sum(lgex)
             dmy = - cdtdy*sum(lgey)
             lugex = HALF*(ugxp+ugxm) * lgex(:)
             lugey = HALF*(ugyp+ugym) * lgey(:)
             dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = f1*(cdtdx*HALF*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                        +       cdtdy*HALF*(ugyp+ugym)*(ergyp(g)-ergym(g)) )
                end do
             else if (fspace_type .eq. 2) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = f1*(cdtdx*HALF*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                        +       cdtdy*HALF*(ergyp(g)+ergym(g))*(ugym-ugyp) )
                end do
             else ! mixed frame
                der(:) = cdtdx * lugex + cdtdy * lugey
             end if
#endif

             ! Convert to conservation form
             rrr = qp(i,j,k,QRHO)
             rur = rrr*qp(i,j,k,QU)
             rvr = rrr*qp(i,j,k,QV)
             rwr = rrr*qp(i,j,k,QW)
             ekenr = HALF*rrr*sum(qp(i,j,k,QU:QW)**2)
             rer = qp(i,j,k,QREINT) + ekenr
#ifdef RADIATION
             err = qp(i,j,k,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewr = rrr - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                          - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
             runewr = rur - cdtdx*(fxy(i+1,j,k,UMX) - fxy(i,j,k,UMX)) &
                          - cdtdy*(fyx(i,j+1,k,UMX) - fyx(i,j,k,UMX))
             rvnewr = rvr - cdtdx*(fxy(i+1,j,k,UMY) - fxy(i,j,k,UMY)) &
                          - cdtdy*(fyx(i,j+1,k,UMY) - fyx(i,j,k,UMY))
             rwnewr = rwr - cdtdx*(fxy(i+1,j,k,UMZ) - fxy(i,j,k,UMZ)) &
                          - cdtdy*(fyx(i,j+1,k,UMZ) - fyx(i,j,k,UMZ))
             renewr = rer - cdtdx*(fxy(i+1,j,k,UEDEN) - fxy(i,j,k,UEDEN)) &
                          - cdtdy*(fyx(i,j+1,k,UEDEN) - fyx(i,j,k,UEDEN))
#ifdef RADIATION
             runewr = runewr + dmx
             rvnewr = rvnewr + dmy
             renewr = renewr + dre
             ernewr = err(:) - cdtdx*(rfxy(i+1,j,k,:) - rfxy(i,j,k,:)) &
                             - cdtdy*(rfyx(i,j+1,k,:) - rfyx(i,j,k,:))  &
                             + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                rwnewr = rwr
                renewr = rer
#ifdef RADIATION
                ernewr = err(:)
#endif
                reset_state = .true.
             endif

             ! Convert back to primitive form
             qpo(i,j,k,QRHO  ) = rrnewr
             qpo(i,j,k,QU    ) = runewr/rrnewr
             qpo(i,j,k,QV    ) = rvnewr/rrnewr
             qpo(i,j,k,QW    ) = rwnewr/rrnewr

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,k,QREINT) = renewr - rhoekenr

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating
                   ! (rho e).
                   qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                        - cdtdx*(fxy(i+1,j,k,UEINT) - fxy(i,j,k,UEINT) + pxav*dux) &
                        - cdtdy*(fyx(i,j+1,k,UEINT) - fyx(i,j,k,UEINT) + pyav*duy)
                endif


                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,k,QPRES) - pxnew - pynew
                   qpo(i,j,k,QPRES) = pnewr
                else
                   ! Update gammae with its transverse terms
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + gexnew + geynew

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                endif
             else
                qpo(i,j,k,QPRES) = qp(i,j,k,QPRES)
                qpo(i,j,k,QGAME) = qp(i,j,k,QGAME)
             endif

             qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qpo, qpo_lo, qpo_hi, i, j, k)

#ifdef RADIATION
             qpo(i,j,k,qrad:qradhi) = ernewr(:)
             qpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k,QPRES)
             qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qradhi)) + qpo(i,j,k,QREINT)
#endif


             !-------------------------------------------------------------------
             ! qzmo state
             !-------------------------------------------------------------------

             pgxp = qx(i+1,j,k-1,GDPRES)
             pgxm = qx(i,j,k-1,GDPRES)
             ugxp = qx(i+1,j,k-1,GDU)
             ugxm = qx(i,j,k-1,GDU)
             gegxp = qx(i+1,j,k-1,GDGAME)
             gegxm = qx(i,j,k-1,GDGAME)
#ifdef RADIATION
             ergxp = qx(i+1,j,k-1,GDERADS:GDERADS-1+ngroups)
             ergxm = qx(i  ,j,k-1,GDERADS:GDERADS-1+ngroups)
#endif

             pgyp = qy(i,j+1,k-1,GDPRES)
             pgym = qy(i,j,k-1,GDPRES)
             ugyp = qy(i,j+1,k-1,GDV)
             ugym = qy(i,j,k-1,GDV)
             gegyp = qy(i,j+1,k-1,GDGAME)
             gegym = qy(i,j,k-1,GDGAME)
#ifdef RADIATION
             ergyp = qy(i,j+1,k-1,GDERADS:GDERADS-1+ngroups)
             ergym = qy(i,j  ,k-1,GDERADS:GDERADS-1+ngroups)
#endif

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm
#ifdef RADIATION
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k-1,QGAMCG) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k-1,QGAMCG))*dux - uxav*dgex )
#else
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k-1,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k-1,QGAMC))*dux - uxav*dgex )
#endif

             duyp = pgyp*ugyp - pgym*ugym
             pyav = HALF*(pgyp+pgym)
             uyav = HALF*(ugyp+ugym)
             geyav = HALF*(gegyp+gegym)
             duy = ugyp-ugym
             dgey = gegyp-gegym
#ifdef RADIATION
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k-1,QGAMCG) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k-1,QGAMCG))*duy - uyav*dgey )
#else
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k-1,QGAMC) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k-1,QGAMC))*duy - uyav*dgey )
#endif

#ifdef RADIATION
             lambda(:) = qaux(i,j,k-1,QLAMS:QLAMS+ngroups-1)

             lgex = lambda(:) * (ergxp(:)-ergxm(:))
             lgey = lambda(:) * (ergyp(:)-ergym(:))
             dmx = - cdtdx*sum(lgex)
             dmy = - cdtdy*sum(lgey)
             lugex = HALF*(ugxp+ugxm) * lgex(:)
             lugey = HALF*(ugyp+ugym) * lgey(:)
             dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = f1*(cdtdx*HALF*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                        +       cdtdy*HALF*(ugyp+ugym)*(ergyp(g)-ergym(g)) )
                end do
             else if (fspace_type .eq. 2) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = f1*(cdtdx*HALF*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                        +       cdtdy*HALF*(ergyp(g)+ergym(g))*(ugym-ugyp) )
                end do
             else ! mixed frame
                der(:) = cdtdx * lugex + cdtdy * lugey
             end if
#endif

             ! Convert to conservation form
             rrl = qm(i,j,k,QRHO)
             rul = rrl*qm(i,j,k,QU)
             rvl = rrl*qm(i,j,k,QV)
             rwl = rrl*qm(i,j,k,QW)
             ekenl = HALF*rrl*sum(qm(i,j,k,QU:QW)**2)
             rel = qm(i,j,k,QREINT) + ekenl
#ifdef RADIATION
             erl = qm(i,j,k,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewl = rrl - cdtdx*(fxy(i+1,j,k-1,URHO) - fxy(i,j,k-1,URHO)) &
                          - cdtdy*(fyx(i,j+1,k-1,URHO) - fyx(i,j,k-1,URHO))
             runewl = rul - cdtdx*(fxy(i+1,j,k-1,UMX) - fxy(i,j,k-1,UMX)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMX) - fyx(i,j,k-1,UMX))
             rvnewl = rvl - cdtdx*(fxy(i+1,j,k-1,UMY) - fxy(i,j,k-1,UMY)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMY) - fyx(i,j,k-1,UMY))
             rwnewl = rwl - cdtdx*(fxy(i+1,j,k-1,UMZ) - fxy(i,j,k-1,UMZ)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMZ) - fyx(i,j,k-1,UMZ))
             renewl = rel - cdtdx*(fxy(i+1,j,k-1,UEDEN) - fxy(i,j,k-1,UEDEN)) &
                          - cdtdy*(fyx(i,j+1,k-1,UEDEN) - fyx(i,j,k-1,UEDEN))
#ifdef RADIATION
             runewl = runewl + dmx
             rvnewl = rvnewl + dmy
             renewl = renewl + dre
             ernewl = erl(:) - cdtdx*(rfxy(i+1,j  ,k-1,:) - rfxy(i,j,k-1,:)) &
                             - cdtdy*(rfyx(i  ,j+1,k-1,:) - rfyx(i,j,k-1,:)) &
                             + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewl < ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                rwnewl = rwl
                renewl = rel
#ifdef RADIATION
                ernewl = erl(:)
#endif
                reset_state = .true.
             endif

             ! Convert back to primitive form
             qmo(i,j,k,QRHO  ) = rrnewl
             qmo(i,j,k,QU    ) = runewl/rrnewl
             qmo(i,j,k,QV    ) = rvnewl/rrnewl
             qmo(i,j,k,QW    ) = rwnewl/rrnewl

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i,j,k,QREINT) = renewl - rhoekenl

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo(i,j,k,QREINT) = qm(i,j,k,QREINT) &
                        - cdtdx*(fxy(i+1,j,k-1,UEINT) - fxy(i,j,k-1,UEINT) + pxav*dux) &
                        - cdtdy*(fyx(i,j+1,k-1,UEINT) - fyx(i,j,k-1,UEINT) + pyav*duy)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i,j,k,QPRES) - pxnew - pynew
                   qmo(i,j,k,QPRES) = pnewl
                else
                   ! Update gammae with its transverse terms
                   qmo(i,j,k,QGAME) = qm(i,j,k,QGAME) + gexnew + geynew

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j,k,QPRES) = qmo(i,j,k,QREINT)*(qmo(i,j,k,QGAME)-ONE)
                endif
             else
                qmo(i,j,k,QPRES) = qm(i,j,k,QPRES)
                qmo(i,j,k,QGAME) = qm(i,j,k,QGAME)
             endif

             qmo(i,j,k,QPRES) = max(qmo(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qmo, qmo_lo, qmo_hi, i, j, k)

#ifdef RADIATION
             qmo(i,j,k,qrad:qradhi) = ernewl(:)
             qmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j,k,QPRES)
             qmo(i,j,k,qreitot) = sum(qmo(i,j,k,qrad:qradhi)) + qmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transxy
#endif

  subroutine reset_edge_state_thermo(qedge, qd_lo, qd_hi, ii, jj, kk)

    use amrex_constants_module, only : ZERO, ONE, HALF
    use network, only : nspec, naux
    use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, &
         QPRES, QREINT, QGAME, QFS, QFX, &
         small_pres, small_temp, &
         ppm_predict_gammae, &
         transverse_use_eos, transverse_reset_rhoe

    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t

    implicit none

    integer, intent(in) :: ii, jj, kk
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    real(rt)        , intent(inout) :: qedge(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    logical :: reset
    type (eos_t) :: eos_state

    !$gpu

    reset = .false.

    if (transverse_reset_rhoe == 1) then
       ! if we are still negative, then we need to reset
       if (qedge(ii,jj,kk,QREINT) < ZERO) then
          reset = .true.

          eos_state % rho = qedge(ii,jj,kk,QRHO)
          eos_state % T = small_temp
          eos_state % xn(:) = qedge(ii,jj,kk,QFS:QFS-1+nspec)
          eos_state % aux(:) = qedge(ii,jj,kk,QFX:QFX-1+naux)

          call eos(eos_input_rt, eos_state)

          qedge(ii,jj,kk,QREINT) = qedge(ii,jj,kk,QRHO)*eos_state % e
          qedge(ii,jj,kk,QPRES) = eos_state % p
       endif

    end if

    if (ppm_predict_gammae == 0 ) then

       if (transverse_use_eos == 1) then
          eos_state % rho = qedge(ii,jj,kk,QRHO)
          eos_state % e   = qedge(ii,jj,kk,QREINT) / qedge(ii,jj,kk,QRHO)
          eos_state % T   = small_temp
          eos_state % xn  = qedge(ii,jj,kk,QFS:QFS+nspec-1)
          eos_state % aux = qedge(ii,jj,kk,QFX:QFX+naux-1)

          call eos(eos_input_re, eos_state)

          qedge(ii,jj,kk,QREINT) = eos_state % e * eos_state % rho
          qedge(ii,jj,kk,QPRES) = max(eos_state % p, small_pres)
       end if

    else
       if (reset) then
          ! recompute the p edge state from this and (rho e), since we reset
          ! qreint  (actually, is this code even necessary?)
          qedge(ii,jj,kk,QPRES) = qedge(ii,jj,kk,QREINT)*(qedge(ii,jj,kk,QGAME)-ONE)
          qedge(ii,jj,kk,QPRES) = max(qedge(ii,jj,kk,QPRES), small_pres)
       end if
    end if

  end subroutine reset_edge_state_thermo

end module transverse_module
