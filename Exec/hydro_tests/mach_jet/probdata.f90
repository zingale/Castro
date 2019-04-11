module probdata_module

    use amrex_fort_module, only : rt => amrex_real

    real(rt), save :: inlet_mach
    logical, save :: do_stratified, do_isentropic

end module probdata_module
