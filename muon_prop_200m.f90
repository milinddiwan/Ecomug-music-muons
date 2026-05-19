program muon_prop_200m
    ! Generate cosmic muons at sea level on a horizontal surface using the
    ! ECoMUG flux parametrization, propagate each through 200 m of water
    ! with MUSIC, and write surviving muons to an ASCII kinematics file.
    !
    ! Generation: Metropolis-Hastings sampling of the sea-level differential
    !   flux dalt(pmu,costh) projected onto a horizontal surface (PDF ∝ flux×cosθ).
    !   Azimuthal angle drawn uniformly.
    !
    ! Propagation: MUSIC full 3-D Monte Carlo (bremsstrahlung, pair production,
    !   nuclear inelastic, ionisation, multiple scattering).
    !
    ! Output columns (muons_200m.dat):
    !   pmu0   – initial momentum [GeV/c]
    !   cth0   – initial cos(zenith angle)
    !   phi0   – initial azimuthal angle [rad]
    !   emu_f  – final total energy [GeV]
    !   x_f    – final x position [cm]
    !   y_f    – final y position [cm]
    !   z_f    – final z (depth) [cm]
    !   cx_f   – final direction cosine x
    !   cy_f   – final direction cosine y
    !   cz_f   – final direction cosine z (cz>0 = still downward)
    !
    ! Compile:
    !   gfortran -O2 -o muon_prop_200m \
    !       randn.f90 corset.f90 corgen.f90 music-sr.f90 ecomug.f90 muon_prop_200m.f90
    !
    ! Data files required in the working directory:
    !   music-eloss-sr.dat  music-double-diff-rock.dat  music-cross-sections-sr.dat
    use ecomug_mod
    implicit none

    integer,  parameter :: N_MU    = 100000          ! muons to generate
    integer,  parameter :: NWARMUP = 5000            ! Metropolis burn-in
    real(8),  parameter :: DEPTH_M  = 200.0d0        ! vertical depth [m]
    real(8),  parameter :: PI       = 3.141592653589793d0
    real(8),  parameter :: MU_MASS  = 0.105655d0     ! muon rest mass [GeV]
    ! Average ionisation energy loss (same value as in ecomug_mod)
    real(8),  parameter :: AAVE     = 2.305d0        ! MeV / (g/cm^2)
    ! Minimum cos(theta) accepted: avoids unreachably large slant depths
    real(8),  parameter :: CTH_MIN_PROP = 0.02d0
    character(len=*), parameter :: MEDIUM  = 'water'
    character(len=*), parameter :: OUTFILE = 'muons_200m.dat'

    ! MUSIC medium parameters (filled by initialize_music)
    real(8) :: n_a, n_z_med, n_rho, n_lambda
    common /rock2/ n_a, n_z_med, n_rho, n_lambda

    real(8), allocatable :: pmu0_arr(:), cth0_arr(:), phi0_arr(:)

    ! Per-muon working variables
    real(8) :: pmu0, cth0, phi0, sth0
    real(8) :: x, y, z, cx, cy, cz, emu, ttime, depth_max

    integer :: i, n_survive, n_above_thr, iout

    integer :: iseed(8) = [72889, 64746, 45, 567, 72828, 27282, 83838, 28293]

    ! ---- Initialisation -----------------------------------------------
    call initialize_music(MEDIUM)
    call random_seed(put=iseed)

    write(*, '(/, a)') '  Sea-level muon generation + 200 m water propagation'
    write(*, '(a, a)')   '  Medium      : ', MEDIUM
    write(*, '(a, f6.2, a)') '  Depth       : ', DEPTH_M, ' m'
    write(*, '(a, f4.2, a)') '  Density     : ', n_rho, ' g/cm^3'
    write(*, '(a, i0)')  '  Muons       : ', N_MU

    ! ---- Generate sea-level horizontal-surface kinematics -------------
    write(*, '(/, a)', advance='no') '  Generating muons ... '
    allocate(pmu0_arr(N_MU), cth0_arr(N_MU), phi0_arr(N_MU))
    call gen_horiz(0.0d0, NWARMUP, N_MU, pmu0_arr, cth0_arr, phi0_arr)
    write(*, '(a)') 'done.'

    ! ---- Open output file ----------------------------------------------
    open(newunit=iout, file=OUTFILE, status='replace', action='write')
    write(iout, '(a)') &
        '# Sea-level muons propagated 200 m in water (MUSIC Monte Carlo)'
    write(iout, '(a)') &
        '# pmu0[GeV/c]  cth0  phi0[rad]  emu_f[GeV]  x_f[cm]  y_f[cm]  z_f[cm]  cx_f  cy_f  cz_f'

    ! ---- Propagation loop ----------------------------------------------
    write(*, '(a)', advance='no') '  Propagating    ... '
    n_survive   = 0
    n_above_thr = 0

    do i = 1, N_MU

        pmu0 = pmu0_arr(i)
        cth0 = cth0_arr(i)
        phi0 = phi0_arr(i)

        ! Skip near-horizontal muons: they can never reach the vertical depth
        if (cth0 < CTH_MIN_PROP) cycle

        ! Slant path to the 200 m vertical plane: L = depth / cos(theta)
        depth_max = (DEPTH_M / cth0) * 100.0d0   ! [cm]
        emu       = sqrt(pmu0*pmu0 + MU_MASS*MU_MASS)   ! total energy [GeV]

        ! Pre-cut: discard muons whose total energy is below the ionisation-only
        ! energy loss along the full slant path.  Such muons cannot reach the
        ! detector plane regardless of fluctuations, so MUSIC would always stop them.
        !   E_loss [GeV] = AAVE [MeV/(g/cm^2)] * rho [g/cm^3] * L [cm] * 1e-3
        if (emu < AAVE * n_rho * depth_max * 1.0d-3) cycle
        n_above_thr = n_above_thr + 1

        sth0 = sqrt(max(0.0d0, 1.0d0 - cth0*cth0))

        ! Start at the surface (origin), muon going downward
        x  = 0.0d0;  y  = 0.0d0;  z  = 0.0d0
        cx = sth0 * cos(phi0)
        cy = sth0 * sin(phi0)
        cz = cth0
        ttime = 0.0d0

        call muon_transport(x, y, z, cx, cy, cz, emu, depth_max, ttime)

        ! emu = 0 means the muon stopped before reaching the 200 m plane
        if (emu > 0.0d0) then
            n_survive = n_survive + 1
            write(iout, '(10(es14.6,1x))') &
                pmu0, cth0, phi0, emu, x, y, z, cx, cy, cz
        end if

    end do

    write(*, '(a)') 'done.'
    close(iout)
    deallocate(pmu0_arr, cth0_arr, phi0_arr)

    ! ---- Summary -------------------------------------------------------
    write(*, '(/, a, i0, a, i0, a, f5.1, a)') &
        '  Above E_min : ', n_above_thr, ' / ', N_MU, &
        '  (', 100.0d0 * n_above_thr / N_MU, ' %)  [passed pre-cut]'
    write(*, '(a, i0, a, i0, a, f5.1, a)') &
        '  Survived    : ', n_survive, ' / ', N_MU, &
        '  (', 100.0d0 * n_survive / N_MU, ' %)  [reached 200 m]'
    write(*, '(a, a, /)') '  Output      : ', OUTFILE

end program muon_prop_200m
