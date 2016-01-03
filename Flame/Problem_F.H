
BL_FORT_PROC_DECL(PROBLEM_CHECKPOINT,problem_checkpoint)(int * int_dir_name, int * len);

BL_FORT_PROC_DECL(PROBLEM_RESTART,problem_restart)(int * int_dir_name, int * len);

extern "C" {
  void ca_flame_step(const int* lo, const int* hi,
		     const BL_FORT_FAB_ARG_3D(state),
		     const BL_FORT_FAB_ARG_3D(grav),
		     const Real& time, const Real& dt_flame);

  void detonation_init();

  void ca_check_ignition(const int* lo, const int* hi,
			 const BL_FORT_FAB_ARG_3D(state),
			 Real ignition_coords[]);

  void ca_check_valid_ignition(const int* lo, const int* hi,
			       const BL_FORT_FAB_ARG_3D(state),
			       const int& det_num,
			       int ignition_conditions[],
			       Real det_xCoord[],
			       Real det_yCoord[],
			       Real det_zCoord[]);

  void ca_set_ignition_points(const Real& time,
			      const int& det_num, int ignition_conditions[],
			      Real det_xCoord[],
			      Real det_yCoord[],
			      Real det_zCoord[]);

  void ca_get_ign_num_max(int* ignnum);

  void ca_get_ignsep(Real* ignsep);
}
