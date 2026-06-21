!
!  test_superres_rates.f90  -  Part C primitives (sbe_superres_ssbe).
!
!  Standalone Fortran unit test of the pure rate / energy-bin-search functions
!  and the cited material data tables. Compiled directly with gfortran against
!  src/ssbe/sbe_superres_ssbe.f90 (self-contained, no SALMON deps). Prints
!  PASS/FAIL and exits 0 (pass) / 1 (fail), so tests/run_all.py can aggregate.
!
program test_superres_rates
    use sbe_superres_ssbe
    implicit none
    integer :: nfail
    real(8) :: v, w, area, dE, s
    integer :: i
    real(8), parameter :: TOL = 1d-9

    nfail = 0

    ! --- nu_saturation: limits and shape (n=1: nu(eps0)=nu_sat(1-1/e)) --------
    call check("nu(eps<=0)=0", nu_saturation(-1d0, 1d14, 0.8d0, 2d0), 0d0, TOL)
    call check("nu(eps0,n=1)=nu_sat(1-1/e)", &
               nu_saturation(0.8d0, 1d14, 0.8d0, 1d0), 1d14*(1d0-exp(-1d0)), 1d0)
    if (nu_saturation(1d3, 1d14, 0.8d0, 2d0) < 0.999d0*1d14) then
        call fail("nu(large eps) should saturate to nu_sat")
    end if

    ! --- bose_factor: N_B(hw=kT) = 1/(e-1); classical limit ------------------
    call check("N_B(hw=kT)=1/(e-1)", bose_factor(1d0, 1d0), 1d0/(exp(1d0)-1d0), TOL)
    call check("N_B(T=0)=0", bose_factor(0.036d0, 0d0), 0d0, TOL)
    ! classical limit hw<<kT: N_B ~ kT/hw - 1/2
    if (abs(bose_factor(1d-3, 1d0) - (1d0/(exp(1d-3)-1d0))) > 1d-6) then
        call fail("N_B classical/series mismatch")
    end if

    ! --- gaussian_bin: peak value and unit area ------------------------------
    s = 0.05d0
    call check("gaussian peak", gaussian_bin(0d0, s), &
               1d0/(sqrt(2d0*3.141592653589793d0)*s), 1d-6)
    area = 0d0; dE = -1d0
    do i = 0, 200000                     ! trapezoid over +-1 (>>sigma)
        area = area + gaussian_bin(-1d0 + 2d0*dble(i)/200000d0, s) * (2d0/200000d0)
    end do
    call check("gaussian unit area", area, 1d0, 1d-4)

    ! --- rect_bin: inside/outside and unit area ------------------------------
    call check("rect inside", rect_bin(0d0, 0.2d0), 1d0/0.2d0, TOL)
    call check("rect outside", rect_bin(0.2d0, 0.2d0), 0d0, TOL)

    ! --- frohlich_hi_factor: zero below hw0, asinh form above ----------------
    call check("frohlich below hw0", frohlich_hi_factor(0.02d0, 0.036d0), 0d0, TOL)
    v = frohlich_hi_factor(1d0, 0.036d0)       ! E=1, hw0=0.036
    w = log(sqrt(1d0/0.036d0-1d0) + sqrt((1d0/0.036d0-1d0)+1d0))/sqrt(1d0)
    call check("frohlich asinh value", v, w, TOL)

    ! --- ii_rate_general: threshold gate and 2^a scaling ---------------------
    call check("II below threshold=0", ii_rate_general(1.0d0, 1.1d0, 1d13, 2d0), 0d0, TOL)
    call check("II scaling 2^a (a=4)", &
               ii_rate_general(3d0,1d0,1d13,4d0)/ii_rate_general(2d0,1d0,1d13,4d0), &
               16d0, 1d-9)

    ! --- bgr_gap_shift_ev: -19 meV at 1e18, -41 meV at 1e19 (K=1.9e-8) -------
    call check("BGR @1e18 = -19 meV", bgr_gap_shift_ev(1d18, 1.9d-8), -0.019d0, 1d-4)
    call check("BGR @1e19 = -41 meV", bgr_gap_shift_ev(1d19, 1.9d-8), -0.0409d0, 1d-3)

    ! --- data tables: sizes and a couple of cited values ---------------------
    if (SI_N_PHONON /= 6)  call fail("Si phonon table size != 6")
    if (GAAS_N_IV /= 5)    call fail("GaAs intervalley table size != 5")
    call check("Si g-LO energy = 63 meV", SI_PHONON_E_MEV(3), 63d0, TOL)
    call check("Si g-LO D = 6.0e8 eV/cm", SI_PHONON_D_1E8EVCM(3), 6.0d0, TOL)
    call check("GaAs Gamma-L D = 10 eV/A", GAAS_IV_D_EVANG(1), 10d0, TOL)
    call check("GaAs hw_LO = 36 meV", GAAS_HW_LO_MEV, 36d0, TOL)

    ! --- unit conversions to a.u. (golden-rule e-ph) -------------------------
    call check("mev_to_ha(1000)=1eV in Ha", mev_to_ha(1000d0), 1d0/27.211386245988d0, 1d-9)
    call check("D[eV/A]/D[eV/cm] ratio = 1e8", &
               d_evang_to_au(1d0)/d_evcm_to_au(1d0), 1d8, 1d-6)
    call check("rho(2.33 g/cm3) in m_e/Bohr^3", rho_gcm3_to_au(2.33d0), 379.0d0, 2d0)

    ! --- golden_rule_prefactor pi D^2/(rho omega) ----------------------------
    call check("golden prefactor pi*4/(3*4)=pi/3", &
               golden_rule_prefactor(2d0,3d0,4d0), 3.141592653589793d0/3d0, 1d-9)
    call check("golden prefactor omega=0 guard", golden_rule_prefactor(2d0,3d0,0d0), 0d0, TOL)

    ! --- eph_thermal_split: normalization + detailed balance -----------------
    block
        real(8) :: fe, fa
        call eph_thermal_split(0.5d0, fe, fa)
        call check("split fe+fa=1", fe+fa, 1d0, TOL)
        call check("split fe/fa=(N+1)/N", fe/fa, 1.5d0/0.5d0, 1d-9)
        call eph_thermal_split(0d0, fe, fa)
        call check("split N=0: fe=1", fe, 1d0, TOL)
        call check("split N=0: fa=0", fa, 0d0, TOL)
    end block

    if (nfail == 0) then
        write(*,'(a)') 'PASS'
        call exit(0)
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'
        call exit(1)
    end if

contains
    subroutine check(name, got, want, tol)
        character(*), intent(in) :: name
        real(8), intent(in) :: got, want, tol
        if (abs(got - want) > tol * max(1d0, abs(want))) then
            write(*,'(a,a,a,es16.8,a,es16.8)') '  FAIL: ', name, '  got=', got, ' want=', want
            nfail = nfail + 1
        end if
    end subroutine check
    subroutine fail(name)
        character(*), intent(in) :: name
        write(*,'(a,a)') '  FAIL: ', name
        nfail = nfail + 1
    end subroutine fail
end program test_superres_rates
