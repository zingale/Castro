
#include <AMReX_BLFort.H>
#include <Castro.H>
#include <Castro_bc_fill_nd.H>
#include <Castro_bc_fill_nd_F.H>
#include <Castro_bc_ext_fill_nd.H>
#include <Castro_generic_fill.H>
#include <Castro_generic_fill_F.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

  void ca_statefill(Box const& bx, FArrayBox& data,
                    const int dcomp, const int numcomp,
                    Geometry const& geom, const Real time,
                    const Vector<BCRec>& bcr, const int bcomp,
                    const int scomp)
  {
      // Make a copy of the raw BCRec data in the format
      // our BC routines can handle (a contiguous array
      // of integers).

      Vector<int> bcrs(2 * AMREX_SPACEDIM * numcomp);

      for (int n = 0; n < numcomp; ++n)
          for (int k = 0; k < 2 * AMREX_SPACEDIM; ++k)
              bcrs[2 * AMREX_SPACEDIM * n + k] = bcr[n].vect()[k];

#ifdef AMREX_USE_CUDA
      int* bc_f = prepare_bc(bcrs.data(), NUM_STATE);
      set_bc_launch_config();
#else
      const int* bc_f = bcrs.data();
#endif

      // This routine either comes in with one component or all NUM_STATE.

      if (numcomp == 1) {

#pragma gpu box(bx)
          denfill(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                  BL_TO_FORTRAN_N_ANYD(data, dcomp),
                  AMREX_INT_ANYD(geom.Domain().loVect()), AMREX_INT_ANYD(geom.Domain().hiVect()),
                  AMREX_REAL_ANYD(geom.CellSize()), AMREX_REAL_ANYD(geom.ProbLo()), time, bc_f);

      }
      else {

          AMREX_ALWAYS_ASSERT(numcomp == NUM_STATE);

#pragma gpu box(bx)
          hypfill(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                  BL_TO_FORTRAN_ANYD(data),
                  AMREX_INT_ANYD(geom.Domain().loVect()), AMREX_INT_ANYD(geom.Domain().hiVect()),
                  AMREX_REAL_ANYD(geom.CellSize()), AMREX_REAL_ANYD(geom.ProbLo()), time, bc_f);

      }

#ifdef AMREX_USE_CUDA
      clean_bc_launch_config();
#endif

      if (numcomp == 1) {

          ca_ext_denfill(bx.loVect(), bx.hiVect(), data.dataPtr(dcomp), data.loVect(), data.hiVect(),
                         geom.Domain().loVect(), geom.Domain().hiVect(),
                         geom.CellSize(), geom.ProbLo(), &time, bc_f);

      }
      else {

          AMREX_ALWAYS_ASSERT(numcomp == NUM_STATE);

          ca_ext_fill(bx.loVect(), bx.hiVect(), data.dataPtr(), data.loVect(), data.hiVect(),
                      geom.Domain().loVect(), geom.Domain().hiVect(),
                      geom.CellSize(), geom.ProbLo(), &time, bc_f);

      }

#ifdef AMREX_USE_CUDA
      clean_bc(bc_f);
#endif
  }

#ifdef GRAVITY
  void ca_phigravfill(Real* phi, const int* phi_lo, const int* phi_hi,
                      const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                      const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = phi_lo[i];
      hi[i] = phi_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    phigravfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                phi, AMREX_INT_ANYD(phi_lo), AMREX_INT_ANYD(phi_hi),
                AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }

  void ca_gravxfill(Real* grav, const int* grav_lo, const int* grav_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = grav_lo[i];
      hi[i] = grav_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    gravxfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
              grav, AMREX_INT_ANYD(grav_lo), AMREX_INT_ANYD(grav_hi),
              AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
              AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
#endif

    ca_ext_gravxfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc(bc_f);
#endif

  }

  void ca_gravyfill(Real* grav, const int* grav_lo, const int* grav_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = grav_lo[i];
      hi[i] = grav_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    gravyfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
              grav, AMREX_INT_ANYD(grav_lo), AMREX_INT_ANYD(grav_hi),
              AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
              AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
#endif

    ca_ext_gravyfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc(bc_f);
#endif

  }

  void ca_gravzfill(Real* grav, const int* grav_lo, const int* grav_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = grav_lo[i];
      hi[i] = grav_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    gravzfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
              grav, AMREX_INT_ANYD(grav_lo), AMREX_INT_ANYD(grav_hi),
              AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
              AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
#endif

    ca_ext_gravzfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc(bc_f);
#endif

  }
#endif

#ifdef ROTATION
  void ca_phirotfill(Real* phi, const int* phi_lo, const int* phi_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = phi_lo[i];
      hi[i] = phi_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    phirotfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               phi, AMREX_INT_ANYD(phi_lo), AMREX_INT_ANYD(phi_hi),
               AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
               AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }

  void ca_rotxfill(Real* rot, const int* rot_lo, const int* rot_hi,
                   const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                   const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = rot_lo[i];
      hi[i] = rot_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    rotxfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
             rot, AMREX_INT_ANYD(rot_lo), AMREX_INT_ANYD(rot_hi),
             AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
             AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }

  void ca_rotyfill(Real* rot, const int* rot_lo, const int* rot_hi,
                   const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                   const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = rot_lo[i];
      hi[i] = rot_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    rotyfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
             rot, AMREX_INT_ANYD(rot_lo), AMREX_INT_ANYD(rot_hi),
             AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
             AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }

  void ca_rotzfill(Real* rot, const int* rot_lo, const int* rot_hi,
                   const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                   const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = rot_lo[i];
      hi[i] = rot_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    rotzfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
             rot, AMREX_INT_ANYD(rot_lo), AMREX_INT_ANYD(rot_hi),
             AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
             AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }
#endif

#ifdef REACTIONS
  void ca_reactfill(Real* react, const int* react_lo, const int* react_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = react_lo[i];
      hi[i] = react_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    reactfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
              react, AMREX_INT_ANYD(react_lo), AMREX_INT_ANYD(react_hi),
              AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
              AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }
#endif

#ifdef RADIATION
  void ca_radfill(Real* rad, const int* rad_lo, const int* rad_hi,
                  const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                  const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = rad_lo[i];
      hi[i] = rad_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    radfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
            rad, AMREX_INT_ANYD(rad_lo), AMREX_INT_ANYD(rad_hi),
            AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
            AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }
#endif

#ifdef __cplusplus
}
#endif
