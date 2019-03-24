#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

using namespace amrex;

void
Castro::do_sdc_update(int m_start, int m_end, Real dt_m) {

  BL_PROFILE("Castro::do_sdc_update()");
    
  // this routine needs to do the update from time node m to m+1
  //
  // We come in with:
  //   A_new[m_start] : this is the advective update at node m_start
  //   A_old[:] : this is the advective source for all nodes at the old iterate
  //
  //   R_old[:] : this is the reaction source for all nodes at the old iterate

  // If we do advection only, then the update is explicit.  If we do
  // reactions, then the update is implicit within a zone.

  // for 4th order reactive SDC, we need to first compute the source, C
  // and do a ghost cell fill on it

  MultiFab& k_start = get_new_data(SDC_k_state_start+m_start);
  MultiFab& k_end = get_new_data(SDC_k_state_start+m_end);

  MultiFab& A_new = get_new_data(SDC_A_state_start+m_start);

  MultiFab tmp;

  MultiFab& A_old_0 = get_old_data(SDC_A_state_start);
  MultiFab& A_old_1 = get_old_data(SDC_A_state_start+1);
  MultiFab& A_old_2 = (fourth_order == 1) ? get_old_data(SDC_A_state_start+2) : tmp;

#ifdef REACTIONS

  MultiFab& R_old_0 = get_old_data(SDC_R_state_start);
  MultiFab& R_old_1 = get_old_data(SDC_R_state_start+1);
  MultiFab& R_old_2 = (fourth_order == 1) ? get_old_data(SDC_R_state_start+2) : tmp;

  MultiFab& R_new = get_new_data(SDC_R_state_start+m_end);

  // C_source will only really be used by 4th order
  MultiFab& C_source = get_new_data(SDC_R_state_start+m_end);

  if (sdc_order == 4) {

    // for 4th order reacting flow, we need to create the "source" C
    // as averages and then convert it to cell centers.  The
    // cell-center version needs to have 2 ghost cells.  We're going
    // to store this in the "new" slot of the SDC "R" state data

    // Reactions are implicit, so they are always at the time node we
    // are integrating too

    for (MFIter mfi(C_source); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      ca_sdc_compute_C4(BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_3D((A_new[mfi]),
                        BL_TO_FORTRAN_3D(A_old_0[mfi]),
                        BL_TO_FORTRAN_3D(A_old_1[mfi]),
                        BL_TO_FORTRAN_3D(A_old_2[mfi]),
                        BL_TO_FORTRAN_3D(R_old_0[mfi]),
                        BL_TO_FORTRAN_3D(R_old_1[mfi]),
                        BL_TO_FORTRAN_3D(R_old_2[mfi]),
                        BL_TO_FORTRAN_3D(C_source[mfi]),
                        &m_start);
    }

    // need to construct the time for this stage -- but it is not really
    // at a single instance in time.  For single level this does not matter.
    Real time = state[SDC_Source_Type].curTime();
    AmrLevel::FillPatch(*this, C_source, C_source.nGrow(), time,
                        SDC_R_state_start+m_end, 0, NUM_STATE);

  }
#endif

  // main update loop -- we are updating k_start to k_end

  FArrayBox U_center;
  FArrayBox C_center;
  FArrayBox U_new_center;
  FArrayBox R_new_tmp;

  for (MFIter mfi(k_start); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();
    const Box& bx1 = mfi.growntilebox(1);

#ifdef REACTIONS
    // advection + reactions
    if (sdc_order == 2) {
      ca_sdc_update_o2(BL_TO_FORTRAN_BOX(bx), &dt_m,
                       BL_TO_FORTRAN_3D(k_start[mfi]),
                       BL_TO_FORTRAN_3D(k_end[mfi]),
                       BL_TO_FORTRAN_3D(A_new[mfi]),
                       BL_TO_FORTRAN_3D(A_old_0[mfi]),
                       BL_TO_FORTRAN_3D(A_old_1[mfi]),
                       BL_TO_FORTRAN_3D(R_old_0[mfi]),
                       BL_TO_FORTRAN_3D(R_old_1[mfi]),
                       &sdc_iteration,
                       &m_start);
    } else {

      // convert the starting U to cell-centered on a fab-by-fab basis
      // -- including one ghost cell
      U_center.resize(bx1, NUM_STATE);

      // TODO: we need k_start to have ghost cells -- is Sborder still valid?
      ca_make_cell_center(BL_TO_FORTRAN_BOX(bx1),
                          BL_TO_FORTRAN_FAB(Sborder[mfi]),
                          BL_TO_FORTRAN_FAB(U_center));

      // convert the C source to cell-centers
      C_center.resize(bx1, NUM_STATE);
      ca_make_cell_center(BL_TO_FORTRAN_BOX(bx1),
                          BL_TO_FORTRAN_FAB(C_source[mfi]),
                          BL_TO_FORTRAN_FAB(C_center));

      // solve for the updated cell-center U using our cell-centered C -- we
      // need to do this with one ghost cell
      U_new_center.resize(bx1, NUM_STATE);
      ca_sdc_update_centers_o4(BL_TO_FORTRAN_BOX(bx1), &dt_m,
                               BL_TO_FORTRAN_3D(U_center),
                               BL_TO_FORTRAN_3D(U_new_center),
                               BL_TO_FORTRAN_3D(C_center),
                               &sdc_iteration);

      // compute R_i and in 1 ghost cell and then convert to <R> in
      // place (only for the interior)
      R_new_tmp.resize(bx1, NUM_STATE);
      ca_instantaneous_react(BL_TO_FORTRAN_BOX(bx1),
                             BL_TO_FORTRAN_3D(U_new_center),
                             BL_TO_FORTRAN_3D(R_new_tmp));

      ca_make_fourth_in_place(BL_TO_FORTRAN_BOX(bx),
                              BL_TO_FORTRAN_FAB(R_new_tmp));

      // now do the conservative update using this <R> to get <U>
      // We'll also need to pass in <C>.  This will only work in
      // the interior
      ca_sdc_conservative_update(BL_TO_FORTRAN_BOX(bx), &dt_m,
                                 BL_TO_FORTRAN_3D(k_start[mfi]),
                                 BL_TO_FORTRAN_3D(k_end[mfi]),
                                 BL_TO_FORTRAN_3D(C_source[mfi]),
                                 BL_TO_FORTRAN_3D(R_new_tmp));

      // we are done with this FAB's C_source, so we can overwrite the
      // data with R_new_tmp -- the new cell-average reactive source

      // TODO

    }
#else
    // pure advection
    if (sdc_order == 2) {
      ca_sdc_update_advection_o2(BL_TO_FORTRAN_BOX(bx), &dt_m,
                                 BL_TO_FORTRAN_3D(k_start[mfi]),
                                 BL_TO_FORTRAN_3D(k_end[mfi]),
                                 BL_TO_FORTRAN_3D(A_new[mfi]),
                                 BL_TO_FORTRAN_3D(A_old_0[mfi]),
                                 BL_TO_FORTRAN_3D(A_old_1[mfi]),
                                 &m_start);
    } else {
      ca_sdc_update_advection_o4(BL_TO_FORTRAN_BOX(bx), &dt_m,
                                 BL_TO_FORTRAN_3D(k_start[mfi]),
                                 BL_TO_FORTRAN_3D(k_end[mfi]),
                                 BL_TO_FORTRAN_3D(A_new[mfi]),
                                 BL_TO_FORTRAN_3D(A_old_0[mfi]),
                                 BL_TO_FORTRAN_3D(A_old_1[mfi]),
                                 BL_TO_FORTRAN_3D(A_old_2[mfi]),
                                 &m_start);
    }
#endif

  }
}


#ifdef REACTIONS
void
Castro::construct_old_react_source(amrex::MultiFab& U_state,
                                   amrex::MultiFab& R_source) {

  BL_PROFILE("Castro::construct_old_react_source()");
    
  // this routine simply fills R_source with the reactive source from
  // state U_state.  Note: it is required that U_state have atleast 2
  // valid ghost cells for 4th order.

  // at this point, k_new has not yet been updated, so it represents
  // the state at the SDC nodes from the previous iteration
  for (MFIter mfi(U_state); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    // construct the reactive source term
    ca_instantaneous_react(BL_TO_FORTRAN_BOX(bx),
                           BL_TO_FORTRAN_3D(U_state[mfi]),
                           BL_TO_FORTRAN_3D(R_source[mfi]));

  }
}
#endif
