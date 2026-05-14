! ================================================================
! muon_range.f90
!
! Monte Carlo computation of muon range statistics using MUSIC.
!
! For each initial energy in E_BINS, propagates N_MU muons through
! the medium configured in initialize_music (default: pure water).
! Reports mean range, standard deviation, lateral spread, and
! mean final angle for each energy bin.
!
! A single call to muon_transport with DEPTH_MAX >> expected range
! propagates each muon all the way until it stops (or exhausts
! path_max inside music). emu = 0 on return indicates the muon stopped.
!
! Compile:
!   gfortran -O2 -o muon_range \
!       randn.f90 corset.f90 corgen.f90 music-sr.f90 muon_range.f90
!
! Data files required in the working directory:
!   music-eloss-sr.dat
!   music-double-diff-rock.dat
!   music-cross-sections-sr.dat
! ================================================================
program muon_range
    implicit none

    ! Number of muons per energy bin
    integer,  parameter :: N_MU     = 1000

    ! Initial energies [GeV]
    integer,  parameter :: N_BINS   = 8
    real(8),  parameter :: E_BINS(N_BINS) = &
        [5.d0, 10.d0, 50.d0, 100.d0, 500.d0, 1000.d0, 5000.d0, 10000.d0]

    ! Maximum depth passed to muon_transport [cm].
    ! Must exceed the range of the highest-energy muon in the medium.
    ! 5×10^6 cm = 50,000 mwe in water, far exceeding ~8,000 mwe for 10 TeV.
    real(8),  parameter :: DEPTH_MAX = 5.0d6   ! [cm]

    real(8),  parameter :: pi        = 3.141592653589793d0

    ! Medium selection: 'water', 'seawater', or 'rock'
    character(len=16), parameter :: MEDIUM = 'rock'

    ! Medium parameters (set by initialize_music, needed for mwe conversion)
    real(8) :: n_a, n_z_med, n_rho, n_lambda
    common /rock2/ n_a, n_z_med, n_rho, n_lambda

    ! Fixed random seed for reproducibility
    integer :: iseed(8) = [72889, 64746, 45, 567, 72828, 27282, 83838, 28293]

    ! Per-muon state
    real(8) :: x, y, z, cx, cy, cz, emu, ttime

    ! Per-energy-bin accumulators
    real(8) :: sum_z, sum_z2, sum_r, sum_theta
    integer :: n_stopped

    ! Derived statistics
    real(8) :: mean_z, std_z, mean_r_lat, mean_theta_f
    real(8) :: range_mwe, std_mwe

    integer :: i_e, i_mu

    ! ---- Initialisation ----
    call initialize_music(MEDIUM)
    call random_seed(put=iseed)

    write(*, '(/, a)')     '  MUSIC muon range simulation'
    write(*, '(a, f4.2, a)')   '  Medium density n_rho  = ', n_rho,    ' g/cm^3'
    write(*, '(a, f5.2, a)')   '  Radiation length X0   = ', n_lambda, ' g/cm^2'
    write(*, '(a, i0, a)') '  Muons per energy bin  : ', N_MU
    write(*, '(a, es8.2, a)') '  Max propagation depth : ', DEPTH_MAX, ' cm'
    write(*, *)
    write(*, '(a)') &
      '  E_mu     <R>        <R>      sigma_R  f_stop  <r_lat>  <theta_f>'
    write(*, '(a)') &
      '  [GeV]    [cm]      [mwe]     [mwe]            [cm]     [deg]'
    write(*, '(a)') &
      '  ----------------------------------------------------------------'

    ! ---- Main loop over energy bins ----
    do i_e = 1, N_BINS

        sum_z     = 0.0d0
        sum_z2    = 0.0d0
        sum_r     = 0.0d0
        sum_theta = 0.0d0
        n_stopped = 0

        do i_mu = 1, N_MU

            ! Initial state: at origin, going straight down, zero travel time
            x = 0.0d0;  y = 0.0d0;  z = 0.0d0
            cx = 0.0d0; cy = 0.0d0; cz = 1.0d0
            emu   = E_BINS(i_e)
            ttime = 0.0d0

            ! Propagate until stopped or DEPTH_MAX is reached
            call muon_transport(x, y, z, cx, cy, cz, emu, DEPTH_MAX, ttime)

            ! emu == 0 means the muon was stopped in the medium
            if (emu == 0.0d0) n_stopped = n_stopped + 1

            sum_z     = sum_z  + z
            sum_z2    = sum_z2 + z*z
            sum_r     = sum_r  + sqrt(x*x + y*y)
            ! cz is cos(theta_final); clamp for numerical safety before acos
            sum_theta = sum_theta + acos(max(-1.0d0, min(1.0d0, cz))) * 180.0d0/pi

        end do

        mean_z      = sum_z  / N_MU
        std_z       = sqrt(max(0.0d0, sum_z2/N_MU - mean_z**2))
        mean_r_lat  = sum_r  / N_MU
        mean_theta_f = sum_theta / N_MU

        ! Convert cm → mwe: 1 mwe = 100 g/cm² = 100/n_rho cm
        range_mwe = mean_z * n_rho / 100.0d0
        std_mwe   = std_z  * n_rho / 100.0d0

        write(*, '(2x, f7.0, 2x, es10.3, 2x, f9.2, 2x, f8.2, 2x, f5.3, 2x, f7.1, 2x, f7.3)') &
            E_BINS(i_e), mean_z, range_mwe, std_mwe, &
            real(n_stopped, 8)/N_MU, mean_r_lat, mean_theta_f

        ! Warn if any muons escaped without stopping (DEPTH_MAX too small)
        if (n_stopped < N_MU) then
            write(*, '(a, i0, a)') &
                '  ** WARNING: ', N_MU - n_stopped, &
                ' muons reached DEPTH_MAX without stopping — increase DEPTH_MAX'
        end if

    end do

    write(*, '(a)') &
      '  ----------------------------------------------------------------'
    write(*, *)
    write(*, '(a)') '  R     = z-displacement of stopped muon [range]'
    write(*, '(a)') '  mwe   = meters water equivalent = g/cm^2 / 100'
    write(*, '(a)') '  r_lat = sqrt(x^2+y^2), lateral displacement from axis'
    write(*, '(a)') '  f_stop = fraction of muons that stopped (should be 1.0)'
    write(*, *)

end program muon_range
