#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

#ifdef FLAME
void Castro::flame_half_dt(MultiFab& s, Real time, Real dt, int ngrow)
{
    BL_PROFILE("Castro::flame_half_dt()");

    const Real strt_time = ParallelDescriptor::second();

    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "\n" << "... Entering flame model and doing half-timestep of burning." << "\n";

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(s, true); mfi.isValid(); ++mfi)
      {

	const Box& bx = mfi.growntilebox(ngrow);

	// Note that box is *not* necessarily just the valid region!
	ca_flame_step(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
		      BL_TO_FORTRAN_3D(s[mfi]),
		      time, 0.5 * dt);

      }

    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "... Leaving flame model after completing half-timestep of burning." << "\n\n";

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);
	if (ParallelDescriptor::IOProcessor()) 
	    std::cout << "flame_half_dt time = " << run_time << '\n';
#ifdef BL_LAZY
	});
#endif
    }
}
#endif
