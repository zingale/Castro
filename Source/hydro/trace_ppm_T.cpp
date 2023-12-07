#include <Castro.H>
#include <Castro_util.H>

#include <cmath>

#include <ppm.H>
#include <reconstruction.H>
#include <flatten.H>

#include <network.H>

using namespace amrex;
using namespace reconstruction;

void
Castro::trace_ppm_T(const Box& bx,
                    const int idir,
                    Array4<Real const> const& U_arr,
                    Array4<Real const> const& rho_inv_arr,
                    Array4<Real const> const& q_arr,
                    Array4<Real const> const& qaux_arr,
                    Array4<Real const> const& srcQ,
                    Array4<Real> const& qm,
                    Array4<Real> const& qp,
#if (AMREX_SPACEDIM < 3)
                    Array4<Real const> const& dloga,
#endif
                    const Box& vbx,
                    const Real dt) {

  // here, lo and hi are the range we loop over -- this can include ghost cells
  // vlo and vhi are the bounds of the valid box (no ghost cells)



  // To allow for easy integration of radiation, we adopt the
  // following conventions:
  //
  // rho : mass density
  // u, v, w : velocities
  // p : gas (hydro) pressure
  // ptot : total pressure (note for pure hydro, this is
  //        just the gas pressure)
  // rhoe_g : gas specific internal energy
  // cgas : sound speed for just the gas contribution
  // cc : total sound speed (including radiation)
  // h_g : gas specific enthalpy / cc**2
  // gam_g : the gas Gamma_1
  // game : gas gamma_e
  //
  // for pure hydro, we will only consider:
  //    rho, u, v, w, ptot, rhoe_g, cc, h_g


  const auto dx = geom.CellSizeArray();

  Real hdt = 0.5_rt * dt;
  Real dtdx = dt / dx[idir];

#ifndef AMREX_USE_GPU
  auto lo = bx.loVect3d();
  auto hi = bx.hiVect3d();
#endif

  auto vlo = vbx.loVect3d();
  auto vhi = vbx.hiVect3d();


#ifndef AMREX_USE_GPU

  // if we're on the CPU, we preprocess the sources over the whole
  // tile up front -- we don't want to trace under a source that is
  // empty. This check only needs to be done over the tile we're
  // working on, since the PPM reconstruction and integration done
  // here is only local to this tile.

  GpuArray<int, NQSRC> do_source_trace;

  for (int n = 0; n < NQSRC; n++) {
    do_source_trace[n] = 0;

    for (int k = lo[2]-2*dg2; k <= hi[2]+2*dg2; k++) {
      for (int j = lo[1]-2*dg1; j <= hi[1]+2*dg1; j++) {
        for (int i = lo[0]-2; i <= hi[0]+2; i++) {
          if (std::abs(srcQ(i,j,k,n)) > 0.0_rt) {
            do_source_trace[n] = 1;
            break;
          }
        }
        if (do_source_trace[n] == 1) {
            break;
        }
      }
      if (do_source_trace[n] == 1) {
          break;
      }
    }
  }
#endif

  // This does the characteristic tracing to build the interface
  // states using the normal predictor only (no transverse terms).
  //
  // For each zone, we construct Im and Ip arrays -- these are the averages
  // of the various primitive state variables under the parabolic
  // interpolant over the region swept out by one of the 3 different
  // characteristic waves.
  //
  // Im is integrating to the left interface of the current zone
  // (which will be used to build the right ("p") state at that interface)
  // and Ip is integrating to the right interface of the current zone
  // (which will be used to build the left ("m") state at that interface).
  //
  //
  // The choice of reference state is designed to minimize the
  // effects of the characteristic projection.  We subtract the I's
  // off of the reference state, project the quantity such that it is
  // in terms of the characteristic variables, and then add all the
  // jumps that are moving toward the interface to the reference
  // state to get the full state on that interface.

  int QUN = -1;
  int QUT = -1;
  int QUTT = -1;

  if (idir == 0) {
    QUN = QU;
    QUT = QV;
    QUTT = QW;
  } else if (idir == 1) {
    QUN = QV;
    QUT = QW;
    QUTT = QU;
  } else if (idir == 2) {
    QUN = QW;
    QUT = QU;
    QUTT = QV;
  }

  Real lsmall_dens = small_dens;
  Real lsmall_pres = small_pres;
  Real lsmall_temp = small_temp;

  // Trace to left and right edges using upwind PPM
  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    Real cc = qaux_arr(i,j,k,QC);
    Real csq_inv = 1.0 / (cc * cc);

    Real un = q_arr(i,j,k,QUN);

    Real rho = q_arr(i,j,k,QRHO);
    Real T = q_arr(i,j,k,QTEMP);

    eos_t eos_state;
    eos_state.rho = rho;
    eos_state.T = T;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = U_arr(i,j,k,UFS+n) * rho_inv_arr(i,j,k);
    }

    eos(eos_input_rt, eos_state);

    // do the parabolic reconstruction and compute the
    // integrals under the characteristic waves
    Real s[nslp];

    Real flat = 1.0;

    if (castro::first_order_hydro) {
        flat = 0.0;
    }
    else if (castro::use_flattening) {
        flat = hydro::flatten(i, j, k, q_arr, QPRES);
    }

    Real sm;
    Real sp;

    // reconstruct density

    Real Ip_rho[3];
    Real Im_rho[3];

    load_stencil(q_arr, idir, i, j, k, QRHO, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile(sm, sp, s[i0], un, cc, dtdx, Ip_rho, Im_rho);

    // reconstruct normal velocity

    Real Ip_un_0;
    Real Im_un_0;
    Real Ip_un_2;
    Real Im_un_2;

    load_stencil(q_arr, idir, i, j, k, QUN, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile_single(sm, sp, s[i0], un-cc, dtdx, Ip_un_0, Im_un_0);
    ppm_int_profile_single(sm, sp, s[i0], un+cc, dtdx, Ip_un_2, Im_un_2);

    // reconstruct T

    Real Ip_T[3];
    Real Im_T[3];

    load_stencil(q_arr, idir, i, j, k, QTEMP, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile(sm, sp, s[i0], un, cc, dtdx, Ip_T, Im_T);

    // reconstruct transverse velocities

    Real Ip_ut_1;
    Real Im_ut_1;
    Real Ip_utt_1;
    Real Im_utt_1;

    load_stencil(q_arr, idir, i, j, k, QUT, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile_single(sm, sp, s[i0], un, dtdx, Ip_ut_1, Im_ut_1);

    load_stencil(q_arr, idir, i, j, k, QUTT, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile_single(sm, sp, s[i0], un, dtdx, Ip_utt_1, Im_utt_1);

    // source terms -- we only trace if the terms in the stencil are non-zero
    int do_trace;

    // density

    Real Ip_src_rho[3] = {0.0_rt};
    Real Im_src_rho[3] = {0.0_rt};

#ifndef AMREX_USE_GPU
    do_trace = do_source_trace[QRHO];
#else
    do_trace = check_trace_source(srcQ, idir, i, j, k, QRHO);
#endif

    if (do_trace) {
        load_stencil(srcQ, idir, i, j, k, QRHO, s);
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile(sm, sp, s[i0], un, cc, dtdx, Ip_src_rho, Im_src_rho);
    }

    // normal velocity

    Real Ip_src_un_0 = 0.0_rt;
    Real Im_src_un_0 = 0.0_rt;
    Real Ip_src_un_2 = 0.0_rt;
    Real Im_src_un_2 = 0.0_rt;

#ifndef AMREX_USE_GPU
    do_trace = do_source_trace[QUN];
#else
    do_trace = check_trace_source(srcQ, idir, i, j, k, QUN);
#endif

    if (do_trace) {
        load_stencil(srcQ, idir, i, j, k, QUN, s);
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile_single(sm, sp, s[i0], un-cc, dtdx, Ip_src_un_0, Im_src_un_0);
        ppm_int_profile_single(sm, sp, s[i0], un+cc, dtdx, Ip_src_un_2, Im_src_un_2);
    }

    // T would go here, but we don't currently fill temperature sources

    // transverse velocities

    Real Ip_src_ut_1 = 0.0_rt;
    Real Im_src_ut_1 = 0.0_rt;

#ifndef AMREX_USE_GPU
    do_trace = do_source_trace[QUT];
#else
    do_trace = check_trace_source(srcQ, idir, i, j, k, QUT);
#endif

    if (do_trace) {
        load_stencil(srcQ, idir, i, j, k, QUT, s);
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile_single(sm, sp, s[i0], un, dtdx, Ip_src_ut_1, Im_src_ut_1);
    }

    Real Ip_src_utt_1 = 0.0_rt;
    Real Im_src_utt_1 = 0.0_rt;

#ifndef AMREX_USE_GPU
    do_trace = do_source_trace[QUTT];
#else
    do_trace = check_trace_source(srcQ, idir, i, j, k, QUTT);
#endif

    if (do_trace) {
        load_stencil(srcQ, idir, i, j, k, QUTT, s);
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile_single(sm, sp, s[i0], un, dtdx, Ip_src_utt_1, Im_src_utt_1);
    }


    // do the passives separately

    // the passive stuff is the same regardless of the tracing

    Real Ip_passive;
    Real Im_passive;

    for (int ipassive = 0; ipassive < npassive; ipassive++) {

        int nc = upassmap(ipassive);
        int n = qpassmap(ipassive);

        load_passive_stencil(U_arr, rho_inv_arr, idir, i, j, k, nc, s);
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile_single(sm, sp, s[i0], un, dtdx, Ip_passive, Im_passive);

        // Plus state on face i

        if ((idir == 0 && i >= vlo[0]) ||
            (idir == 1 && j >= vlo[1]) ||
            (idir == 2 && k >= vlo[2])) {

            // We have
            //
            // q_l = q_ref - Proj{(q_ref - I)}
            //
            // and Proj{} represents the characteristic projection.
            // But for these, there is only 1-wave that matters, the u
            // wave, so no projection is needed.  Since we are not
            // projecting, the reference state doesn't matter

            qp(i,j,k,n) = Im_passive;
        }

        // Minus state on face i+1
        if (idir == 0 && i <= vhi[0]) {
            qm(i+1,j,k,n) = Ip_passive;
        } else if (idir == 1 && j <= vhi[1]) {
            qm(i,j+1,k,n) = Ip_passive;
        } else if (idir == 2 && k <= vhi[2]) {
            qm(i,j,k+1,n) = Ip_passive;
        }
    }

    // now the rest of the state

    // plus state on face i

    if ((idir == 0 && i >= vlo[0]) ||
        (idir == 1 && j >= vlo[1]) ||
        (idir == 2 && k >= vlo[2])) {

      // Set the reference state
      // This will be the fastest moving state to the left --
      // this is the method that Miller & Colella and Colella &
      // Woodward use. These papers don't include the effect of
      // the source term in constructing the reference state, but
      // we do because the source term could be large relative to
      // the current state (for example, consider the gravity term
      // acting on a fluid that is initially uniformly at rest --
      // the dt * g term will be the only non-zero contribution).
      // We ignore the effect of the source term for gamma.
      Real rho_ref = Im_rho[0] + hdt * Im_src_rho[0];
      Real un_ref = Im_un_0 + hdt * Im_src_un_0;

      Real T_ref = Im_T[0];

      rho_ref = amrex::max(rho_ref, lsmall_dens);
      T_ref = amrex::max(T_ref, lsmall_temp);


      // *m are the jumps carried by un-c
      // *p are the jumps carried by un+c

      // Note: for the transverse velocities, the jump is carried
      //       only by the u wave (the contact)


      // we also add the sources here so they participate in the tracing
      Real drhom = rho_ref - Im_rho[0] - hdt * Im_src_rho[0];
      Real dum = un_ref - Im_un_0 - hdt * Im_src_un_0;
      Real dTm = T_ref - Im_T[0];

      Real drho = rho_ref - Im_rho[1] - hdt * Im_src_rho[1];
      Real dT = T_ref - Im_T[1];

      Real drhop = rho_ref - Im_rho[2] - hdt * Im_src_rho[2];
      Real dup = un_ref - Im_un_2 - hdt*Im_src_un_2;
      Real dTp = T_ref - Im_T[2];

      // {rho, u, T} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper.  This is simply (l . dq), where dq = qref - I(q)

      Real alpham = 0.5_rt * (drhom * eos_state.dpdr + dTm * eos_state.dpdT - dum * cc * rho) * csq_inv;
      Real alphap = 0.5_rt * (drhop * eos_state.dpdr + dTp * eos_state.dpdT + dup * cc * rho) * csq_inv;
      Real alpha0r = (drho * (cc * cc - eos_state.dpdr) - dT * eos_state.dpdT) * csq_inv;

      alpham = un-cc > 0.0_rt ? 0.0_rt : -alpham;
      alphap = un+cc > 0.0_rt ? 0.0_rt : -alphap;
      alpha0r = un > 0.0_rt ? 0.0_rt : -alpha0r;

      // The final interface states are just
      // q_s = q_ref - sum(l . dq) r
      // note that the a{mpz}right as defined above have the minus already
      qp(i,j,k,QRHO) = amrex::max(lsmall_dens, rho_ref + alphap + alpham + alpha0r);
      qp(i,j,k,QUN) = un_ref + (alphap - alpham) * cc / rho;
      qp(i,j,k,QTEMP) = amrex::max(lsmall_temp,
                                   T_ref + (alphap + alpham) * (cc * cc - eos_state.dpdr) / eos_state.dpdT + alpha0r * eos_state.dpdr / eos_state.dpdT);

      // Transverse velocities -- there's no projection here, so we
      // don't need a reference state.  We only care about the state
      // traced under the middle wave

      // Recall that I already takes the limit of the parabola
      // in the event that the wave is not moving toward the
      // interface
      qp(i,j,k,QUT) = Im_ut_1 + hdt*Im_src_ut_1;
      qp(i,j,k,QUTT) = Im_utt_1 + hdt*Im_src_utt_1;

    }

    // minus state on face i + 1

    if ((idir == 0 && i <= vhi[0]) ||
        (idir == 1 && j <= vhi[1]) ||
        (idir == 2 && k <= vhi[2])) {

      // Set the reference state
      // This will be the fastest moving state to the right
      Real rho_ref = Ip_rho[2] + hdt * Ip_src_rho[2];
      Real un_ref = Ip_un_2 + hdt * Ip_src_un_2;

      Real T_ref = Ip_T[2];

      rho_ref = amrex::max(rho_ref, lsmall_dens);
      T_ref = amrex::max(T_ref, lsmall_temp);


      // *m are the jumps carried by u-c
      // *p are the jumps carried by u+c

      Real drhom = rho_ref - Ip_rho[0] - hdt * Ip_src_rho[0];
      Real dum = un_ref - Ip_un_0 - hdt * Ip_src_un_0;
      Real dTm  = T_ref - Ip_T[0];

      Real drho = rho_ref - Ip_rho[1] - hdt * Ip_src_rho[1];
      Real dT = T_ref - Ip_T[1];

      Real drhop = rho_ref - Ip_rho[2] - hdt * Ip_src_rho[2];
      Real dup = un_ref - Ip_un_2 - hdt * Ip_src_un_2;
      Real dTp = T_ref - Ip_T[2];

      // {rho, u, T} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper.  This is simply (l . dq), where dq = qref - I(q)

      Real alpham = 0.5_rt * (drhom * eos_state.dpdr + dTm * eos_state.dpdT - dum * cc * rho) * csq_inv;
      Real alphap = 0.5_rt * (drhop * eos_state.dpdr + dTp * eos_state.dpdT + dup * cc * rho) * csq_inv;
      Real alpha0r = (drho * (cc * cc - eos_state.dpdr) - dT * eos_state.dpdT) * csq_inv;

      alpham = un-cc > 0.0_rt ? -alpham : 0.0_rt;
      alphap = un+cc > 0.0_rt ? -alphap : 0.0_rt;
      alpha0r = un > 0.0_rt ? -alpha0r : 0.0_rt;

      // The final interface states are just
      // q_s = q_ref - sum (l . dq) r
      // note that the a{mpz}left as defined above have the minus already
      if (idir == 0) {
        qm(i+1,j,k,QRHO) = amrex::max(lsmall_dens, rho_ref +  alphap + alpham + alpha0r);
        qm(i+1,j,k,QUN) = un_ref + (alphap - alpham) * cc / rho;
        qm(i+1,j,k,QTEMP) = amrex::max(lsmall_temp,
                                       T_ref + (alphap + alpham) * (cc * cc - eos_state.dpdr) / eos_state.dpdT + alpha0r * eos_state.dpdr / eos_state.dpdT);

        // transverse velocities
        qm(i+1,j,k,QUT) = Ip_ut_1 + hdt*Ip_src_ut_1;
        qm(i+1,j,k,QUTT) = Ip_utt_1 + hdt*Ip_src_utt_1;

      } else if (idir == 1) {
        qm(i,j+1,k,QRHO) = amrex::max(lsmall_dens, rho_ref +  alphap + alpham + alpha0r);
        qm(i,j+1,k,QUN) = un_ref + (alphap - alpham) * cc / rho;
        qm(i,j+1,k,QTEMP) = amrex::max(lsmall_temp,
                                       T_ref + (alphap + alpham) * (cc * cc - eos_state.dpdr) / eos_state.dpdT + alpha0r * eos_state.dpdr / eos_state.dpdT);

        // transverse velocities
        qm(i,j+1,k,QUT) = Ip_ut_1 + hdt*Ip_src_ut_1;
        qm(i,j+1,k,QUTT) = Ip_utt_1 + hdt*Ip_src_utt_1;

      } else if (idir == 2) {
        qm(i,j,k+1,QRHO) = amrex::max(lsmall_dens, rho_ref +  alphap + alpham + alpha0r);
        qm(i,j,k+1,QUN) = un_ref + (alphap - alpham) * cc / rho;
        qm(i,j,k+1,QTEMP) = amrex::max(lsmall_temp,
                                       T_ref + (alphap + alpham) * (cc * cc - eos_state.dpdr) / eos_state.dpdT + alpha0r * eos_state.dpdr / eos_state.dpdT);

        // transverse velocities
        qm(i,j,k+1,QUT) = Ip_ut_1 + hdt*Ip_src_ut_1;
        qm(i,j,k+1,QUTT) = Ip_utt_1 + hdt*Ip_src_utt_1;
      }

    }

    // geometry source terms
    // TODO: what is the geometry term for temperature?

#if (AMREX_SPACEDIM < 3)
    // these only apply for x states (idir = 0)
    if (idir == 0 && dloga(i,j,k) != 0.0_rt) {

      Real courn = dt/dx[0]*(cc+std::abs(un));
      Real eta = (1.0_rt - courn)/(cc*dt*std::abs(dloga(i,j,k)));
      Real dlogatmp = amrex::min(eta, 1.0_rt)*dloga(i,j,k);
      Real sourcr = -0.5_rt*dt*rho*dlogatmp*un;
      Real sourcp = sourcr * cc * cc;
      Real source = sourcp*((q_arr(i,j,k,QPRES) + q_arr(i,j,k,QREINT))/rho) * csq_inv;

      if (i <= vhi[0]) {
        qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr;
        qm(i+1,j,k,QRHO) = amrex::max(qm(i+1,j,k,QRHO), lsmall_dens);
        qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp;
        qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source;
      }

      if (i >= vlo[0]) {
        qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr;
        qp(i,j,k,QRHO) = amrex::max(qp(i,j,k,QRHO), lsmall_dens);
        qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp;
        qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source;
      }
    }
#endif

    // now compute the interface states for p and (rho e) via the EOS

    if ((idir == 0 && i >= vlo[0]) ||
        (idir == 1 && j >= vlo[1]) ||
        (idir == 2 && k >= vlo[2])) {

        eos_state.rho = qp(i,j,k,QRHO);
        eos_state.T = qp(i,j,k,QTEMP);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = qp(i,j,k,QFS+n);
        }

        eos(eos_input_rt, eos_state);

        qp(i,j,k,QREINT) = qp(i,j,k,QRHO) * eos_state.e;
        qp(i,j,k,QPRES) = amrex::max(lsmall_pres, eos_state.p);

    }

    if ((idir == 0 && i <= vhi[0]) ||
        (idir == 1 && j <= vhi[1]) ||
        (idir == 2 && k <= vhi[2])) {

        if (idir == 0) {
            eos_state.rho = qm(i+1,j,k,QRHO);
            eos_state.T = qm(i+1,j,k,QTEMP);
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = qm(i+1,j,k,QFS+n);
            }

            eos(eos_input_rt, eos_state);

            qm(i+1,j,k,QREINT) = qm(i+1,j,k,QRHO) * eos_state.e;
            qm(i+1,j,k,QPRES) = amrex::max(lsmall_pres, eos_state.p);

        } else if (idir == 1) {
            eos_state.rho = qm(i,j+1,k,QRHO);
            eos_state.T = qm(i,j+1,k,QTEMP);
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = qm(i,j+1,k,QFS+n);
            }

            eos(eos_input_rt, eos_state);

            qm(i,j+1,k,QREINT) = qm(i,j+1,k,QRHO) * eos_state.e;
            qm(i,j+1,k,QPRES) = amrex::max(lsmall_pres, eos_state.p);

        } else if (idir == 2) {
            eos_state.rho = qm(i,j,k+1,QRHO);
            eos_state.T = qm(i,j,k+1,QTEMP);
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = qm(i,j,k+1,QFS+n);
            }

            eos(eos_input_rt, eos_state);

            qm(i,j,k+1,QREINT) = qm(i,j,k+1,QRHO) * eos_state.e;
            qm(i,j,k+1,QPRES) = amrex::max(lsmall_pres, eos_state.p);

        }

    }


  });
}


