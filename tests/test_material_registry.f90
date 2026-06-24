!
!  test_material_registry.f90  -  central per-material constant registry.
!
!  Checks get_material_params(): every SBE dissipation channel auto-selects its
!  constants (dielectric, impact-ionization fit, electron-phonon table) through
!  this one struct, so adding a material is a single cited `case` block. The test
!  pins the GaAs and Si entries (the values existing runs depend on), the
!  Si_cb alias, the diamond/zincblende flag, the polar-LO convention, and that
!  an unknown material is reported as not-found (the callers then stop).
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_material_registry
    use sbe_superres_ssbe
    implicit none
    integer :: nfail
    type(s_material_params) :: ga, si, sicb, cds, unk
    real(8), parameter :: TOL = 1d-12

    nfail = 0
    ga   = get_material_params('GaAs')
    si   = get_material_params('Si')
    sicb = get_material_params('Si_cb')
    cds  = get_material_params('CdS')
    unk  = get_material_params('Ge')

    ! --- found / unknown ---------------------------------------------------
    if (.not. ga%found)   call bad('GaAs not found')
    if (.not. si%found)   call bad('Si not found')
    if (.not. sicb%found) call bad('Si_cb not found')
    if (.not. cds%found)  call bad('CdS not found')
    if (unk%found)        call bad('unknown material reported as found')

    ! cubic materials: al(1:3) box is isotropic = a
    call chk('GaAs cell x', ga%cell_au(1), 10.68d0)
    call chk('GaAs cell isotropic', ga%cell_au(3), ga%cell_au(1))
    call chk('Si cell isotropic', si%cell_au(2), si%cell_au(1))

    ! --- GaAs (zincblende, polar) -- the legacy defaults -------------------
    if (ga%is_diamond)             call bad('GaAs flagged diamond')
    if (.not. ga%eph_polar)        call bad('GaAs not flagged polar')
    call chk('GaAs eps0', ga%eps0, 12.9d0)
    call chk('GaAs eps_inf', ga%eps_inf, 10.89d0)
    if (trim(ga%ii_form) /= 'stobbe_quartic') call bad('GaAs ii_form')
    call chk('GaAs ii_exponent', ga%ii_exponent, 4d0)
    call chk('GaAs ii_threshold_ev', ga%ii_threshold_ev, 2.1d0)
    call ichk('GaAs eph_nph', ga%eph_nph, 6)
    call chk('GaAs LO hw', ga%eph_hw_mev(1), 36.0d0)
    ! polar-LO weight convention: raw LO weight = sum of the intervalley weights
    call chk('GaAs LO dominant weight', ga%eph_wraw(1), sum(ga%eph_wraw(2:6)))
    call chk('GaAs nu_sat', ga%eph_nu_sat_si, 1.0d14)

    ! --- Si (diamond, non-polar) -------------------------------------------
    if (.not. si%is_diamond)  call bad('Si not flagged diamond')
    if (si%eph_polar)         call bad('Si flagged polar')
    call chk('Si eps0', si%eps0, 11.7d0)
    call chk('Si eps_inf=eps0 (non-polar)', si%eps_inf, 11.7d0)
    if (trim(si%ii_form) /= 'keldysh_quadratic') call bad('Si ii_form')
    call chk('Si ii_exponent', si%ii_exponent, 2d0)
    call chk('Si ii_threshold_ev', si%ii_threshold_ev, 1.1d0)
    call ichk('Si eph_nph', si%eph_nph, 6)
    call chk('Si first IV hw', si%eph_hw_mev(1), 10.0d0)
    call chk('Si nu_sat', si%eph_nu_sat_si, 1.3d14)
    call chk('Si lattice', si%a_lattice_au, 10.26d0)

    ! Si_cb is an alias of Si (same diamond constants)
    call chk('Si_cb eps0 == Si', sicb%eps0, si%eps0)
    if (.not. sicb%is_diamond) call bad('Si_cb not diamond')

    ! --- CdS (wurtzite, polar, non-centrosymmetric, orthorhombic cell) ------
    if (cds%is_diamond)     call bad('CdS flagged diamond')
    if (.not. cds%eph_polar) call bad('CdS not flagged polar')
    call chk('CdS eps0', cds%eps0, 9.0d0)
    call chk('CdS eps_inf', cds%eps_inf, 5.3d0)
    if (trim(cds%ii_form) /= 'keldysh_quadratic') call bad('CdS ii_form')
    call chk('CdS ii_exponent', cds%ii_exponent, 2d0)
    call chk('CdS ii_threshold_ev (1.5 Eg)', cds%ii_threshold_ev, 3.6d0)
    call chk('CdS nu_sat', cds%eph_nu_sat_si, 2.0d14)
    call ichk('CdS eph_nph (single Frohlich LO)', cds%eph_nph, 1)
    call chk('CdS LO hw', cds%eph_hw_mev(1), 38.0d0)
    ! orthorhombic cell is anisotropic: b = a*sqrt(3), c independent
    if (abs(cds%cell_au(2) - cds%cell_au(1)) < 1d-6) call bad('CdS cell not anisotropic (b==a)')
    call chk('CdS cell b = a*sqrt3', cds%cell_au(2), cds%cell_au(1)*sqrt(3d0), 1d-3)
    if (abs(cds%cell_au(3) - cds%cell_au(1)) < 1d-6) call bad('CdS cell c == a')

    ! every entry's phonon weights must be positive (normalizable downstream)
    if (sum(ga%eph_wraw(1:ga%eph_nph)) <= 0d0) call bad('GaAs eph weights non-positive')
    if (sum(si%eph_wraw(1:si%eph_nph)) <= 0d0) call bad('Si eph weights non-positive')
    if (sum(cds%eph_wraw(1:cds%eph_nph)) <= 0d0) call bad('CdS eph weights non-positive')

    if (nfail == 0) then
        write(*,'(a)') 'PASS'; call exit(0)
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'; call exit(1)
    end if

contains
    subroutine chk(name, got, want, tol_in)
        character(*), intent(in) :: name
        real(8), intent(in) :: got, want
        real(8), intent(in), optional :: tol_in   ! NB: 'tol' would alias TOL (case-insensitive)
        real(8) :: t
        t = TOL
        if (present(tol_in)) t = tol_in
        if (abs(got-want) > t*max(1d0,abs(want))) then
            write(*,'(a,a,a,es16.8,a,es16.8)') '  FAIL: ',name,' got=',got,' want=',want
            nfail = nfail + 1
        end if
    end subroutine chk
    subroutine ichk(name, got, want)
        character(*), intent(in) :: name
        integer, intent(in) :: got, want
        if (got /= want) then
            write(*,'(a,a,a,i0,a,i0)') '  FAIL: ',name,' got=',got,' want=',want
            nfail = nfail + 1
        end if
    end subroutine ichk
    subroutine bad(name)
        character(*), intent(in) :: name
        write(*,'(a,a)') '  FAIL: ', name
        nfail = nfail + 1
    end subroutine bad
end program test_material_registry
