!
!  validate_auger_c.f90 -- MANUAL validator (not part of run_all.py):
!  effective Auger coefficient C(n) of the ring Auger kernel vs the
!  SOURCE-VERIFIED literature tables (wiki/07 sec.7).
!
!  Prepares Fermi-Dirac electron/hole populations (n = p, quasi-Fermi levels
!  found by bisection) on a REAL Si primitive EPM spectrum (SYSNAME_eigen/
!  _k.data in the current directory), calls auger_interk_dpop in the linear
!  regime, and extracts
!      C_eff = R_vol / (n^2 p)   [cm^6/s]
!  for a ladder of carrier densities. Reference: Si eeh (n-type), experiment
!  Dziewior-Schmid C_n = 2.8e-31 cm^6/s; [L90] pure direct AR reproduces it.
!  Our kernel's magnitude is tied to the cited Keldysh II fit by detailed
!  balance (no separate C), so this is an ORDER-OF-MAGNITUDE consistency
!  check, sensitive to the grid (4^3 here) and the energy broadening sigma
!  (printed for two sigmas).
!
!  Usage:  gfortran -O2 -ffree-line-length-none src/ssbe/sbe_superres_ssbe.f90 \
!              tests/validate_auger_c.f90 -o validate_auger_c
!          cd <dir with Si_prim_eigen.data / Si_prim_k.data (4^3, nstate=20)>
!          ./validate_auger_c
!
program validate_auger_c
    use sbe_superres_ssbe, only: auger_interk_dpop, mp_grid_triple, &
                                 get_material_params, s_material_params
    implicit none
    integer, parameter :: NK = 64, NB = 20, NVB = 4
    integer, parameter :: KN(3) = (/ 4, 4, 4 /)
    real(8), parameter :: PI = 3.14159265358979324d0
    real(8), parameter :: A_LAT = 10.26d0                 ! Si a [Bohr]
    real(8), parameter :: AU_EV = 27.211386245988d0
    real(8), parameter :: AU_FS = 0.02418884326505d0      ! 1 a.u.t = this many fs... (inverse below)
    real(8), parameter :: A0_CM = 0.52917721067d-8
    real(8), parameter :: KB_HA = 3.166811563d-6          ! [Ha/K]
    real(8), parameter :: KT = 300d0 * KB_HA

    type(s_material_params) :: mp
    real(8) :: eval(NB, NK), kred(3, NK), f(NB, NK), dpop(NB, NK)
    integer :: kidx(3, NK), klut(0:NK-1)
    real(8) :: bmat(3,3), vol, n_val, kf, qtf2, wp2, q2reg, blen2
    real(8) :: ecbm, eth, pref_au, tau, resid, maxresid
    real(8) :: ncm3(3), n_au, ne_cell, mu_c, mu_v, rdot_cell, r_vol, c_au, c_cm6s
    real(8) :: sig(2)
    integer :: ik, ib, i, isig, id, fh
    character(64) :: dum

    ! --- lattice: FCC primitive, reciprocal rows b_i = (2pi/a)(-1,1,1)cyc ----
    bmat(1,:) = (2d0*PI/A_LAT) * (/ -1d0,  1d0,  1d0 /)
    bmat(2,:) = (2d0*PI/A_LAT) * (/  1d0, -1d0,  1d0 /)
    bmat(3,:) = (2d0*PI/A_LAT) * (/  1d0,  1d0, -1d0 /)
    vol = A_LAT**3 / 4d0

    ! --- read the EPM ground state (formats of gs_info_ssbe) ------------------
    open(newunit=fh, file='Si_prim_eigen.data', status='old', action='read')
    read(fh,'(a)') dum; read(fh,'(a)') dum; read(fh,'(a)') dum
    do ik = 1, NK
        read(fh,'(a)') dum
        do ib = 1, NB
            read(fh,*) i, eval(ib, ik)
        end do
    end do
    close(fh)
    open(newunit=fh, file='Si_prim_k.data', status='old', action='read')
    do
        read(fh,'(a)') dum
        if (dum(1:1) /= '#') then
            backspace(fh); exit
        end if
    end do
    do ik = 1, NK
        read(fh,*) i, kred(1,ik), kred(2,ik), kred(3,ik)
    end do
    close(fh)

    ! --- MP momentum map ------------------------------------------------------
    klut = 0;  maxresid = 0d0
    do ik = 1, NK
        call mp_grid_triple(kred(:,ik), KN, kidx(:,ik), resid)
        maxresid = max(maxresid, resid)
        klut(kidx(1,ik) + KN(1)*(kidx(2,ik) + KN(2)*kidx(3,ik))) = ik
    end do
    if (maxresid > 1d-6 .or. minval(klut) < 1) stop 'not an MP grid'

    ! --- screening + II constants (same recipe as apply_ii_interk_ring) ------
    mp = get_material_params('Si')
    n_val = 8d0 / vol
    kf    = (3d0*PI*PI*n_val)**(1d0/3d0)
    qtf2  = 4d0*kf/PI
    wp2   = 4d0*PI*n_val
    q2reg = huge(1d0)
    do id = 1, 3
        blen2 = dot_product(bmat(id,:), bmat(id,:))
        q2reg = min(q2reg, blen2/dble(KN(id))**2)
    end do
    q2reg = 0.25d0*q2reg
    ! P [s^-1 eV^-a] -> [1/(Ha^a a.u.t)]: x (a.u.t in s) x (eV/Ha)^a
    ! (same conversion as bloch_solver: * (au_fs*1e-15) * au_ev**a)
    pref_au = mp%ii_prefactor * (2.418884326505d-17) * AU_EV**mp%ii_exponent
    eth  = mp%ii_threshold_ev / AU_EV
    ecbm = minval(eval(NVB+1, :))
    tau  = 1d-2                                  ! linear regime [a.u.t]
    sig  = (/ 0.02d0, 0.01d0 /)                  ! broadening sensitivity [Ha]

    write(*,'(a)') '# effective ring-Auger C(n) on the Si primitive 4^3 EPM spectrum, T=300K, n=p'
    write(*,'(a)') '#   reference (source-verified, wiki/07 sec.7): Si eeh C_n = 2.8e-31 cm^6/s'
    write(*,'(a,f6.2,a,es10.3,a,f5.2,a)') '#   eps_inf=', mp%eps_inf, '  P=', mp%ii_prefactor, &
        ' s^-1 eV^-a  a=', mp%ii_exponent, '  (Keldysh II fit; Auger = detailed-balance partner)'
    ncm3 = (/ 1d18, 1d19, 1d20 /)
    do isig = 1, 2
        write(*,'(a,f6.3,a)') '# sigma = ', sig(isig), ' Ha'
        do i = 1, 3
            n_au = ncm3(i) * A0_CM**3
            ne_cell = n_au * vol
            call fill_fd(ne_cell, mu_c, mu_v)
            call auger_interk_dpop(NK, NB, eval, f, 2d0, 0d0, ecbm, eth, &
                                   pref_au, mp%ii_exponent, NVB, NVB+1, kidx, KN, klut, &
                                   bmat, mp%eps_inf, qtf2, wp2, 0d0, q2reg, sig(isig), tau, dpop)
            ! net CB electron loss per cell per time -> volume rate
            rdot_cell = -sum(dpop(NVB+1:NB, :)) / dble(NK) / tau
            r_vol = rdot_cell / vol
            c_au  = r_vol / max(n_au**3, 1d-300)
            c_cm6s = c_au * A0_CM**6 / 2.418884326505d-17
            write(*,'(a,es9.2,a,es11.3,a,f7.1)') '  n=p = ', ncm3(i), ' cm^-3   C_eff = ', &
                c_cm6s, ' cm^6/s   ratio to 2.8e-31: ', c_cm6s/2.8d-31
        end do
    end do

contains

    ! FD populations: electrons in CB with mu_c, holes in VB with mu_v (n = p),
    ! quasi-Fermi levels by bisection on the per-cell counts (occ_max = 2).
    subroutine fill_fd(target_cell, muc, muv)
        implicit none
        real(8), intent(in)  :: target_cell
        real(8), intent(out) :: muc, muv
        real(8) :: lo, hi, mid, cnt
        integer :: it2, ikk, ibb
        ! electrons in the conduction bands
        lo = ecbm - 1d0; hi = ecbm + 1d0
        do it2 = 1, 80
            mid = 0.5d0*(lo+hi); cnt = 0d0
            do ikk = 1, NK
                do ibb = NVB+1, NB
                    cnt = cnt + 2d0/(exp(min(max((eval(ibb,ikk)-mid)/KT,-60d0),60d0))+1d0)
                end do
            end do
            cnt = cnt/dble(NK)
            if (cnt < target_cell) then
                lo = mid
            else
                hi = mid
            end if
        end do
        muc = 0.5d0*(lo+hi)
        ! holes in the valence bands
        lo = maxval(eval(NVB,:)) - 1d0; hi = maxval(eval(NVB,:)) + 1d0
        do it2 = 1, 80
            mid = 0.5d0*(lo+hi); cnt = 0d0
            do ikk = 1, NK
                do ibb = 1, NVB
                    cnt = cnt + 2d0 - 2d0/(exp(min(max((eval(ibb,ikk)-mid)/KT,-60d0),60d0))+1d0)
                end do
            end do
            cnt = cnt/dble(NK)
            if (cnt > target_cell) then
                lo = mid
            else
                hi = mid
            end if
        end do
        muv = 0.5d0*(lo+hi)
        ! assemble populations
        do ikk = 1, NK
            do ibb = 1, NVB
                f(ibb,ikk) = 2d0/(exp(min(max((eval(ibb,ikk)-muv)/KT,-60d0),60d0))+1d0)
            end do
            do ibb = NVB+1, NB
                f(ibb,ikk) = 2d0/(exp(min(max((eval(ibb,ikk)-muc)/KT,-60d0),60d0))+1d0)
            end do
        end do
    end subroutine fill_fd

end program validate_auger_c
