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
    type(s_material_params) :: ga, si, sicb, cds, gr, unk
    real(8), parameter :: TOL = 1d-12

    nfail = 0
    ga   = get_material_params('GaAs')
    si   = get_material_params('Si')
    sicb = get_material_params('Si_cb')
    cds  = get_material_params('CdS')
    gr   = get_material_params('graphene')
    unk  = get_material_params('Ge')

    ! --- found / unknown ---------------------------------------------------
    if (.not. ga%found)   call bad('GaAs not found')
    if (.not. si%found)   call bad('Si not found')
    if (.not. sicb%found) call bad('Si_cb not found')
    if (.not. cds%found)  call bad('CdS not found')
    if (.not. gr%found)   call bad('graphene not found')
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

    ! --- provenance gates: GaAs/Si channels are all cited (.true.) -----------
    if (.not. (ga%ii_ok .and. ga%eph_ok .and. ga%eeh_ok .and. ga%coulomb_ok)) &
        call bad('GaAs channels should all be supported')
    if (.not. (si%ii_ok .and. si%eph_ok .and. si%eeh_ok .and. si%coulomb_ok)) &
        call bad('Si channels should all be supported')

    ! --- CdS (wurtzite, orthorhombic cell): channels enabled per the md spec --
    if (cds%is_diamond)  call bad('CdS flagged diamond')
    ! orthorhombic cell is anisotropic: b = a*sqrt(3), c independent
    if (abs(cds%cell_au(2) - cds%cell_au(1)) < 1d-6) call bad('CdS cell not anisotropic (b==a)')
    call chk('CdS cell b = a*sqrt3', cds%cell_au(2), cds%cell_au(1)*sqrt(3d0), 1d-3)
    if (abs(cds%cell_au(3) - cds%cell_au(1)) < 1d-6) call bad('CdS cell c == a')
    ! cited dissipation channels from the md: e-ph (Frohlich), Coulomb, II
    if (.not. cds%eph_ok)     call bad('CdS e-ph (Frohlich) should be enabled (md cited)')
    if (.not. cds%coulomb_ok) call bad('CdS Coulomb should be enabled (md cited eps)')
    if (.not. cds%ii_ok)      call bad('CdS impact ionization should be enabled (md cited E_th)')
    ! carrier-carrier (e-e) is FORBIDDEN for CdS: there is no cited CdS rate (the
    ! generic 1e14 scale is GaAs/Si only) -> eeh_ok must stay .false. (enabling
    ! yn_sbe_eeh aborts). The CdS literature only fixes an e-e *timescale*.
    if (cds%eeh_ok)           call bad('CdS carrier-carrier must be forbidden (no cited CdS rate)')
    ! e-ph is the cited Frohlich polar-LO at 38 meV
    if (.not. cds%eph_polar)  call bad('CdS should be flagged polar (Frohlich)')
    call ichk('CdS eph_nph (single Frohlich LO)', cds%eph_nph, 1)
    call chk('CdS LO hw = 38 meV', cds%eph_hw_mev(1), 38.0d0)
    call chk('CdS eps0', cds%eps0, 8.9d0)
    call chk('CdS II threshold = 1.5 Eg', cds%ii_threshold_ev, 3.6d0)
    ! II prefactor is a fit parameter (md): registry holds a sentinel, not a value
    if (cds%ii_prefactor > 0d0) call bad('CdS II prefactor must be a sentinel (fit parameter)')
    ! Auger recombination is FORBIDDEN for CdS: there is NO verified CdS Auger
    ! coefficient. The previously-shipped "C = 2.0e-30 cm^6/s [Haury PRB 57,
    ! 11513 (1998)]" was a fabricated reference (the real Haury 1997 PRL is on
    ! CdMnTe ferromagnetism, unrelated) and was removed -> auger_ok must be
    ! .false. and the registry default C must be 0 (a user may still opt in via
    ! sbe_auger_c_cm6s).
    if (cds%auger_ok)  call bad('CdS Auger must be forbidden (no verified C; Haury ref was fabricated)')
    call chk('CdS Auger C default = 0 (no verified value)', cds%auger_c_cm6s, 0d0, 1d-40)
    ! No material ships a verified Auger C default yet (GaAs/Si/graphene cited
    ! coefficients are the nonlocal-Auger task -- wiki/07).
    if (ga%auger_ok) call bad('GaAs Auger must be forbidden (not yet cited)')
    if (si%auger_ok) call bad('Si Auger must be forbidden (not yet cited)')

    ! --- graphene (GAPLESS Dirac): ONLY e-ph cited; all gap-based channels off ---
    ! V^A = 0 (centrosymmetric) -> flagged diamond like Si.
    if (.not. gr%is_diamond)  call bad('graphene should be flagged diamond (V^A=0, centrosymmetric)')
    ! e-ph IS cited (Piscanec 2004): two optical modes E2g(196) + A1'(160 meV).
    if (.not. gr%eph_ok)      call bad('graphene e-ph should be enabled (Piscanec 2004)')
    if (gr%eph_polar)         call bad('graphene must NOT be flagged polar (no Frohlich LO)')
    call ichk('graphene eph_nph (E2g + A1'')', gr%eph_nph, 2)
    call chk('graphene E2g hw = 196 meV', gr%eph_hw_mev(1), 196.0d0)
    call chk('graphene A1'' hw = 160 meV', gr%eph_hw_mev(2), 160.0d0)
    ! A1' (K, Kohn anomaly, x2 GW) must dominate the E2g (Gamma) weight.
    if (gr%eph_wraw(2) <= gr%eph_wraw(1)) call bad('graphene A1'' weight should dominate E2g (Kohn anomaly)')
    if (sum(gr%eph_wraw(1:gr%eph_nph)) <= 0d0) call bad('graphene eph weights non-positive')
    ! Gap-based / many-body channels FORBIDDEN on the gapless cone -- EXCEPT
    ! the Auger, which is now the CITED 2D Rana branch (R07; wiki/00 TODO-1):
    if (gr%ii_ok)      call bad('graphene impact ionization must be forbidden (gapless: no E_th)')
    if (.not. gr%auger_ok) call bad('graphene Auger should be ENABLED (2D Rana branch, R07 cited)')
    if (.not. gr%auger_2d_rana) call bad('graphene Auger must be flagged 2D-Rana (not C n^3)')
    call chk('graphene Rana v_F = 1e8 cm/s (a.u.)', gr%rana_vf_au, 1.0d8/2.18769126364d8, 1d-12)
    call chk('graphene Rana eps_r = 10 (R07 Fig.4)', gr%rana_eps_r, 10.0d0)
    if (gr%eeh_ok)     call bad('graphene carrier-carrier must be forbidden (no cited FD rate)')
    if (gr%coulomb_ok) call bad('graphene Coulomb HF must be forbidden (needs 2D V(q), not yet wired)')
    ! the bulk materials must NOT carry the 2D flag; Si must NOT take the
    ! dynamic free-carrier lambda (Burt [L90]: lambda = 0 is the Si physics),
    ! GaAs DOES (registry dyn_lambda_ok gates the ring-kernel lambda^2(n(t))).
    if (ga%auger_2d_rana .or. si%auger_2d_rana .or. cds%auger_2d_rana) &
        call bad('bulk materials must not be flagged 2D-Rana')
    if (si%dyn_lambda_ok) call bad('Si must keep lambda=0 in the ring kernels (Burt/L90)')
    if (.not. ga%dyn_lambda_ok) call bad('GaAs should take the dynamic lambda^2(n(t)) (Part-G)')
    ! A2 Cp/Cn (cited): Si 0.99/2.8 [L90/Dziewior-Schmid], GaAs ~4.8 [S14 hhe~5x eeh]
    call chk('Si Cp/Cn = 0.354', si%ii_cpcn, 0.99d0/2.8d0, 1d-12)
    call chk('GaAs Cp/Cn = 4.8 (S14)', ga%ii_cpcn, 4.8d0, 1d-12)
    if (cds%ii_cpcn > 0d0) call bad('CdS must have no cited Cp/Cn')
    ! A4 acoustic constants: Si [Jacoboni-Reggiani], GaAs [Fischetti-Laux 7.0 eV],
    ! graphene [Hwang-Das Sarma 16 eV, input-overridable], CdS [Rode PRB 2,
    ! 1012 (1970): E1 = 14.5 eV -- maintainer-supplied; MUST run TF-screened]
    call chk('Si acoustic Xi_d = 9.0 eV', si%eph_ac_xi_ev, 9.0d0)
    call chk('GaAs acoustic Xi_d = 7.0 eV', ga%eph_ac_xi_ev, 7.0d0)
    call chk('graphene acoustic D = 16 eV', gr%eph_ac_xi_ev, 16.0d0)
    call chk('CdS acoustic E1 = 14.5 eV (Rode 1970)', cds%eph_ac_xi_ev, 14.5d0)
    if (si%eph_ac_cs_cmps <= 0d0 .or. ga%eph_ac_cs_cmps <= 0d0 .or. &
        gr%eph_ac_cs_cmps <= 0d0 .or. cds%eph_ac_cs_cmps <= 0d0) &
        call bad('acoustic sound velocities missing')

    ! every cited phonon table's weights must be positive (normalizable)
    if (sum(ga%eph_wraw(1:ga%eph_nph)) <= 0d0) call bad('GaAs eph weights non-positive')
    if (sum(si%eph_wraw(1:si%eph_nph)) <= 0d0) call bad('Si eph weights non-positive')

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
