#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

#include "Problem.H"

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
