program muon_prop_200m
    ! Generate cosmic muons at sea level on a horizontal surface using the
    ! ECoMUG flux parametrization, propagate each through 200 m of water
    ! with MUSIC, and write surviving muons to an ASCII kinematics file.
    !
    ! Muons are generated one at a time (streaming Metropolis chain) until
    ! N_TARGET survivors have been collected, avoiding large pre-allocation.
    !
    ! Generation: Metropolis-Hastings sampling of the sea-level differential
    !   flux dalt(pmu,costh) projected onto a horizontal surface (PDF ∝ flux×cosθ).
    !   Azimuthal angle drawn uniformly.
    !
    ! Propagation: MUSIC full 3-D Monte Carlo (bremsstrahlung, pair production,
    !   nuclear inelastic, ionisation, multiple scattering).
    !
    ! Pre-cut: muons whose total energy is below the ionisation-only energy
    !   loss along the full slant path are skipped before calling MUSIC.
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
    !   cz_f   – final direction cosine z  (>0 = still downward)
    !
    ! Compile:
    !   gfortran -O2 -o muon_prop_200m \
    !       randn.f90 corset.f90 corgen.f90 music-sr.f90 ecomug.f90 muon_prop_200m.f90
    !
    ! Data files required in the working directory:
    !   music-eloss-sr.dat  music-double-diff-rock.dat  music-cross-sections-sr.dat
    use ecomug_mod
    implicit none

    integer,  parameter :: N_TARGET     = 10000     ! stop after this many survivors
    integer,  parameter :: NWARMUP      = 5000      ! Metropolis burn-in steps
    real(8),  parameter :: DEPTH_M      = 200.0d0   ! vertical depth [m]
    real(8),  parameter :: MU_MASS      = 0.105655d0! muon rest mass [GeV]
    real(8),  parameter :: AAVE         = 2.305d0   ! ionisation loss [MeV/(g/cm^2)]
    real(8),  parameter :: CTH_MIN_PROP = 0.02d0    ! drop near-horizontal muons
    character(len=*), parameter :: MEDIUM  = 'water'
    character(len=*), parameter :: OUTFILE = 'muons_200m.dat'

    ! MUSIC medium parameters (filled by initialize_music)
    real(8) :: n_a, n_z_med, n_rho, n_lambda
    common /rock2/ n_a, n_z_med, n_rho, n_lambda

    ! Persistent Metropolis chain state
    real(8) :: s_p, s_x, s_f

    ! Per-muon working variables
    real(8) :: pmu0, cth0, phi0, sth0
    real(8) :: x, y, z, cx, cy, cz, emu, ttime, depth_max

    integer :: n_survive, iout
    integer :: n_total, n_above_thr

    integer :: iseed(8) = [72889, 64746, 45, 567, 72828, 27282, 83838, 28293]

    ! ---- Initialisation -----------------------------------------------
    call initialize_music(MEDIUM)
    call random_seed(put=iseed)

    write(*, '(/, a)') '  Sea-level muon generation + 200 m water propagation'
    write(*, '(a, a)')       '  Medium   : ', MEDIUM
    write(*, '(a, f6.2, a)') '  Depth    : ', DEPTH_M, ' m'
    write(*, '(a, f4.2, a)') '  Density  : ', n_rho, ' g/cm^3'
    write(*, '(a, i0, a)')   '  Target   : ', N_TARGET, ' survivors'

    ! ---- Warm up the Metropolis chain ---------------------------------
    write(*, '(/, a)', advance='no') '  Warming up generator ... '
    call gen_horiz_init(0.0d0, NWARMUP, s_p, s_x, s_f)
    write(*, '(a)') 'done.'

    ! ---- Open output file ----------------------------------------------
    open(newunit=iout, file=OUTFILE, status='replace', action='write')
    write(iout, '(a)') &
        '# Sea-level muons propagated 200 m in water (MUSIC Monte Carlo)'
    write(iout, '(a)') &
        '# pmu0[GeV/c]  cth0  phi0[rad]  emu_f[GeV]  x_f[cm]  y_f[cm]  z_f[cm]  cx_f  cy_f  cz_f'

    ! ---- Main loop: generate one at a time, stop when N_TARGET reached --
    write(*, '(a)', advance='no') '  Propagating    ... '
    n_survive   = 0
    n_above_thr = 0
    n_total     = 0

    do while (n_survive < N_TARGET)

        call gen_horiz_step(0.0d0, s_p, s_x, s_f, pmu0, cth0, phi0)
        n_total = n_total + 1

        ! Drop near-horizontal tracks (can never reach 200 m vertical depth)
        if (cth0 < CTH_MIN_PROP) cycle

        ! Slant path to the 200 m vertical plane
        depth_max = (DEPTH_M / cth0) * 100.0d0   ! [cm]
        emu       = sqrt(pmu0*pmu0 + MU_MASS*MU_MASS)

        ! Pre-cut: ionisation floor — skip guaranteed stoppers
        if (emu < AAVE * n_rho * depth_max * 1.0d-3) cycle
        n_above_thr = n_above_thr + 1

        sth0 = sqrt(max(0.0d0, 1.0d0 - cth0*cth0))
        x  = 0.0d0;  y  = 0.0d0;  z  = 0.0d0
        cx = sth0 * cos(phi0)
        cy = sth0 * sin(phi0)
        cz = cth0
        ttime = 0.0d0

        call muon_transport(x, y, z, cx, cy, cz, emu, depth_max, ttime)

        if (emu > 0.0d0) then
            n_survive = n_survive + 1
            write(iout, '(10(es14.6,1x))') &
                pmu0, cth0, phi0, emu, x, y, z, cx, cy, cz
        end if

    end do

    write(*, '(a)') 'done.'
    close(iout)

    ! ---- Summary -------------------------------------------------------
    write(*, '(/, a, i0)')              '  Generated    : ', n_total
    write(*, '(a, i0, a, f5.2, a)')    &
        '  Above E_min  : ', n_above_thr, &
        '  (', 100.0d0 * n_above_thr / n_total, ' %)'
    write(*, '(a, i0, a, f5.2, a)')    &
        '  Survived     : ', n_survive,   &
        '  (', 100.0d0 * n_survive / n_total, ' %)'
    write(*, '(a, a, /)')               '  Output       : ', OUTFILE

end program muon_prop_200m
