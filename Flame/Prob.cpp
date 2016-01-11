#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

#include "ParmParse.H"

#ifdef problem_source_term
void Castro::problem_source_half_dt(MultiFab& s, Real time, Real dt, int is_new)
{
    BL_PROFILE("Castro::problem_source_half_dt()");

    const Real strt_time = ParallelDescriptor::second();

    // Get the do_flame parameter from the inputs file, and return if we're not using the flame model.

    ParmParse pp("castro");

    int do_flame = 0;

    pp.query("do_flame", do_flame);

    if (do_flame == 0) return;

    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "\n" << "... Entering flame model and doing half-timestep of flame evolution." << "\n";

    int ngrow = s.nGrow();

    if (ngrow < 2 && ParallelDescriptor::IOProcessor())
      BoxLib::Abort("State MultiFab doesn't have enough ghost cells in problem_source_half_dt.");

    // Note that we are intentionally not tiling this due to the way the flame model is written.
   
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(s, false); mfi.isValid(); ++mfi)
    {

	// Note that box is *not* necessarily just the valid region!
	const Box& bx = mfi.growntilebox(ngrow);

	if (is_new) {
	  MultiFab& g = get_new_data(Gravity_Type);
	  ca_flame_step(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
			BL_TO_FORTRAN_3D(s[mfi]),
			BL_TO_FORTRAN_3D(g[mfi]),
			time, 0.5 * dt);
	} else {
	  MultiFab& g = get_old_data(Gravity_Type);
	  ca_flame_step(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
			BL_TO_FORTRAN_3D(s[mfi]),
			BL_TO_FORTRAN_3D(g[mfi]),
			time, 0.5 * dt);
	}

    }

    // Initialize the detonation data, if we have not already done so.

    detonation_init();

    // Determine the maximum number of ignition points to allow from the Fortran module.

    int ignNumMax;

    ca_get_ign_num_max(&ignNumMax);

    Real ignition_coords[3 * ignNumMax];

    for (int i = 0; i < 3 * ignNumMax; i++)
      ignition_coords[i] = 0.0;

    // Now determine the initial ignition points. Due to the logic used, we should not 
    // use OpenMP in this section, but this check is cheap so that's not so bad.

    for (MFIter mfi(s, true); mfi.isValid(); ++mfi)
    {

	// Note that box is *not* necessarily just the valid region!
	const Box& bx = mfi.validbox();

	ca_check_ignition(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
			  BL_TO_FORTRAN_3D(s[mfi]),
			  ignition_coords);

    }

    // Now communicate the ignition conditions to all processors

    int nprocs = ParallelDescriptor::NProcs();

    int bcast_order[nprocs];
    
    int det_num = nprocs;

    // MPI::Comm::Allgather( det_num, 1, MPI_INTEGER,    
    // 			  bcast_order, 1, MPI_INTEGER,
    // 			  ParallelDescriptor::Communicator() );

    // det_num = sum(bcast_order);

    // If there are no ignition conditions, then we are done.
    // Otherwise we need to allocate space for the detonation coordinates.

    if (det_num != 0) {

      Real recv_buf[3 * det_num];
      int displ[3 * nprocs];

      // displ[0] = 0;
      // bcast_order[0] = 3 * bcast_order[0];

      // for (int i = 1; i <= nprocs-1; i++) {
      // 	displ[i] = displ[i-1] + bcast_order[i-1];
      // 	bcast_order[i] = 3 * bcast_order[i];
      // }

      for (int i = 0; i < nprocs; i++) {
	displ[i] = ParallelDescriptor::MyProc();
	bcast_order[i] = 3;
      }

      MPI_Datatype typ = ParallelDescriptor::Mpi_typemap<Real>::type();

      MPI_Allgatherv( &(ignition_coords[0]),     
		      3,
		      typ,
		      &(recv_buf[0]),         
		      &(bcast_order[0]),      
		      &(displ[0]),
		      typ,
		      ParallelDescriptor::Communicator() );

      Real radius[det_num];

      Real det_xCoord[det_num];
      Real det_yCoord[det_num];
      Real det_zCoord[det_num];

      Real IgnSep;

      ca_get_ignsep(&IgnSep);

      for (int i = 0; i < det_num; i++) {

	det_xCoord[i] = recv_buf[ (i-1) * 3    ];
	radius[i] = det_xCoord[i] * det_xCoord[i];

	det_yCoord[i] = recv_buf[ (i-1) * 3 + 1];
	radius[i] = radius[i] + det_yCoord[i] * det_yCoord[i];
	  
	det_zCoord[i] = recv_buf[ (i-1) * 3 + 2];
	radius[i] = radius[i] + det_zCoord[i] * det_zCoord[i];

	radius[i] = std::sqrt( radius[i] );

      }

      bool detMask[det_num];
      Real dist;

      for (int i = 0; i < det_num; i++) detMask[i] = false;

      // Now determine if any detonation points overlap
      for (int l = 0; l < det_num; l++) {

	// Check with other detonation points found this time step
	for (int ll = l+1; ll < det_num; l++) {

	  // Make sure we are not comparing the same point and that it has
	  // not been previously removed
	  if ( ! (detMask[l] || detMask[ll]) ) {

	    dist = std::sqrt( ( det_xCoord[l] - det_xCoord[ll] ) * ( det_xCoord[l] - det_xCoord[ll] ) + 
			      ( det_yCoord[l] - det_yCoord[ll] ) * ( det_yCoord[l] - det_yCoord[ll] ) + 
			      ( det_zCoord[l] - det_zCoord[ll] ) * ( det_zCoord[l] - det_zCoord[ll] ) );

	    if ( dist <= IgnSep ) {

              if (radius[l] >= radius[ll])
		// get rid of 2nd radius point
		detMask[ll] = true;
	      else  
		// get rid of 1st radius point
		detMask[l] = true;

	    }

	  }

	}

      }

      // Now adjust and move arrays
      bool have_det = false;

      for (int l = 0; l < det_num; l++)
	if (detMask[l])
	  have_det = true;

      if (have_det) {

	int l = 0;

	while (l < det_num) {

	    if (detMask[l]) {

	      // Remove this point

	      det_num = det_num - 1;

	      for (int ll = l; ll < det_num; ll++) {

		det_xCoord[ll] = det_xCoord[ll+1];
		det_yCoord[ll] = det_yCoord[ll+1];
		det_zCoord[ll] = det_zCoord[ll+1];
		detMask[ll] = detMask[ll+1];

	      }

	      // Move l back 1
	      l = l - 1;

	    }

	    l = l + 1;

       }

      }

      // Now we need to check whether any of the detonation points are invalid
      // because they are in ash. Again, we are not parallelizing this region.

      int ignition_conditions[det_num];

      for (MFIter mfi(s, true); mfi.isValid(); ++mfi)
      {

	// Note that box is *not* necessarily just the valid region!
	const Box& bx = mfi.validbox();

	ca_check_valid_ignition(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
				BL_TO_FORTRAN_3D(s[mfi]),
				det_num,
				ignition_conditions,
				det_xCoord, det_yCoord, det_zCoord);

      }

      // Now do a global reduce so we can discard all invalidated ignition points.

      ParallelDescriptor::ReduceIntSum(ignition_conditions, det_num);
      
      for (int i = 0; i < det_num; i++)
	if (ignition_conditions[i] < det_num)
	  ignition_conditions[i] = 0;
	else
	  ignition_conditions[i] = 1;

      // Finally, save it to the Fortran module so that the burn unit can access it.

      ca_set_ignition_points(time, det_num, ignition_conditions, det_xCoord, det_yCoord, det_zCoord);

    }

    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "... Leaving flame model after completing half-timestep of flame evolution." << "\n\n";

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);
	if (ParallelDescriptor::IOProcessor()) 
	    std::cout << "problem_source_half_dt time = " << run_time << '\n';
#ifdef BL_LAZY
	});
#endif
    }
}
#endif
