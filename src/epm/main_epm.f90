!
!  Entry point for theory='epm': local empirical pseudopotential method
!  ground-state calculation (Cohen-Bergstresser form factors), producing
!  SYSNAME_k.data / SYSNAME_eigen.data / SYSNAME_tm.data in the format
!  consumed by gs_info_ssbe (EPM -> SBE chain).
!
!  *** DEPRECATED ***
!  The maintained EPM in this repository is now the Python implementation
!  (epm_gaas_reference.py), which is more flexible (continuous q-space /
!  DeePseudopot Zunger & NN local potentials, DFT-fitted form factors via
!  tools/dft_to_epm). This Fortran path is kept only for backward
!  compatibility and is no longer the primary development target.
!
subroutine main_epm(icomm)
    use communication
    use salmon_global, only: sysname, base_directory
    use epm_solver
    implicit none
    integer, intent(in) :: icomm
    type(s_epm_info) :: epm
    integer :: irank, nproc

    call comm_get_groupinfo(icomm, irank, nproc)

    if (comm_is_root(irank)) then
        write(*,'(a)') '# EPM: local empirical pseudopotential calculation (Cohen-Bergstresser)'
        write(*,'(a)') '# WARNING: the Fortran theory=''epm'' path is DEPRECATED.'
        write(*,'(a)') '#          The maintained EPM is the Python epm_gaas_reference.py'
        write(*,'(a)') '#          (DFT-fitted form factors via tools/dft_to_epm).'
    end if

    call init_epm_info(epm, icomm)
    call run_epm_calculation(epm, icomm)

    if (comm_is_root(irank)) then
        call write_epm_files(epm, sysname, base_directory)
    end if

    call comm_sync_all(icomm)

    call finalize_epm_info(epm)

    return
end subroutine main_epm
