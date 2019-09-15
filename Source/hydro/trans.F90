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

                call reset_edge_state_thermo(lq2o)
                
                if (d == -1) then
                   q2mo(i,j,k,:) = lq2o(:)
                else
                   q2po(i,j,k,:) = lq2o(:)
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

                call reset_edge_state_thermo(lqo)
                
                if (d == -1) then
                   qmo(i,j,k,:) = lqo(:)
                else
                   qpo(i,j,k,:) = lqo(:)
                end if

             end do

          end do
       end do
    end do

  end subroutine transyz
#endif



  subroutine reset_edge_state_thermo(qedge)

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

    real(rt), intent(inout) :: qedge(NQ)

    logical :: reset
    type (eos_t) :: eos_state

    !$gpu

    reset = .false.

    if (transverse_reset_rhoe == 1) then
       ! if we are still negative, then we need to reset
       if (qedge(QREINT) < ZERO) then
          reset = .true.

          eos_state % rho = qedge(QRHO)
          eos_state % T = small_temp
          eos_state % xn(:) = qedge(QFS:QFS-1+nspec)
          eos_state % aux(:) = qedge(QFX:QFX-1+naux)

          call eos(eos_input_rt, eos_state)

          qedge(QREINT) = qedge(QRHO)*eos_state % e
          qedge(QPRES) = eos_state % p
       endif

    end if

    if (ppm_predict_gammae == 0 ) then

       if (transverse_use_eos == 1) then
          eos_state % rho = qedge(QRHO)
          eos_state % e   = qedge(QREINT) / qedge(QRHO)
          eos_state % T   = small_temp
          eos_state % xn  = qedge(QFS:QFS+nspec-1)
          eos_state % aux = qedge(QFX:QFX+naux-1)

          call eos(eos_input_re, eos_state)

          qedge(QREINT) = eos_state % e * eos_state % rho
          qedge(QPRES) = max(eos_state % p, small_pres)
       end if

    else
       if (reset) then
          ! recompute the p edge state from this and (rho e), since we reset
          ! qreint  (actually, is this code even necessary?)
          qedge(QPRES) = qedge(QREINT)*(qedge(QGAME)-ONE)
          qedge(QPRES) = max(qedge(QPRES), small_pres)
       end if
    end if

  end subroutine reset_edge_state_thermo

end module transverse_module
