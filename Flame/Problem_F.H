
BL_FORT_PROC_DECL(PROBLEM_CHECKPOINT,problem_checkpoint)(int * int_dir_name, int * len);

BL_FORT_PROC_DECL(PROBLEM_RESTART,problem_restart)(int * int_dir_name, int * len);

extern "C" {
  void ca_flame_step(const int* lo, const int* hi,
		     const BL_FORT_FAB_ARG_3D(state),
		     const BL_FORT_FAB_ARG_3D(grav),
		     const Real& time, const Real& dt_flame);
}

