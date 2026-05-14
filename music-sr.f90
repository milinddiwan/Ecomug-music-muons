! ================================================================
! MUSIC — MUon SImulation Code
!
! Authors : P. Antonioli, C. Ghetti, E.V. Korolkova,
!           V.A. Kudryavtsev, G. Sartorelli
! Contact : v.kudryavtsev@sheffield.ac.uk
!
! Modified by V. Kudryavtsev, April 1999:
!   - NORMCO replaced by CORSET/CORGEN for correlated Gaussian sampling
!   - RNDM (obsolete) replaced by RANMAR; both RANLUX and RANMAR used
!   - All angles in radians (previously degrees) for Linux compatibility
!   - Sea water density and cross-sections used by default
!
! MUSIC is a 3D Monte Carlo muon propagation code.
! Simulated processes:
!   jp=1  Bremsstrahlung         (Van Ginneken model)
!   jp=2  Pair production
!   jp=3  Nuclear inelastic scattering
!   jp=4  Ionisation (continuous below threshold v = 10^-3)
!   Multiple Coulomb scattering  (Highland 1975 / Lipari & Stanev 1991)
!
! Units: path lengths in g/cm^2, angles in radians, energies in GeV.
!
! Converted from fixed-form Fortran 77 to free-form Fortran 90.
! dated May 14, 2026  
! ================================================================


! ----------------------------------------------------------------
! initialize_music
!
! Read lookup tables from data files and populate COMMON blocks.
! Must be called once before any call to muon_transport / music.
!
! Argument:
!   medium  case-insensitive string selecting the propagation medium:
!             'water'     pure water  (default physics)
!             'seawater'  sea water   (same Z/A, density 1.035 g/cm^3)
!             'rock'      standard rock
!           Any other value is rejected with an error stop.
!
! Data files expected in the working directory:
!   music-eloss-sr.dat           energy-loss coefficients (81×6)
!   music-double-diff-rock.dat   inelastic double-differential cross-section
!   music-cross-sections-sr.dat  integral cross-sections (100×71×4 × 2 arrays)
!
! COMMON blocks filled:
!   /rock2/   medium parameters (Z, A, density, radiation length)
!   /coenlo/  continuous energy-loss coefficient arrays
!   /sig/     cross-section arrays and energy-transfer index j0
!   /ang/     angular scattering distribution table
! ----------------------------------------------------------------
subroutine initialize_music(medium)
    implicit none
    character(len=*), intent(in) :: medium

    real(8) :: n_a, n_z, n_rho, n_lambda

    integer  :: j0
    real(4)  :: v2(100), cs1(100,71,4), cs2(100,71,4)
    real(4)  :: emulo(81,6), emc(81)
    real(4)  :: ema(61), va(30), anga(50), dang(50,30,61)

    common /sig/    v2, cs1, cs2, j0
    common /coenlo/ emulo, emc
    common /ang/    ema, va, anga, dang
    common /rock2/  n_a, n_z, n_rho, n_lambda

    integer :: i, il, j, k
    character(len=32) :: med

    ! Normalise to lowercase for case-insensitive comparison
    med = medium
    call str_lower(med)

    select case (trim(med))
    case ('water')
        n_z      = 6.60d0
        n_a      = 11.89d0
        n_rho    = 1.00d0
        n_lambda = 36.1d0
        print *, 'music: medium = pure water'
    case ('seawater', 'sea_water', 'sea water')
        n_z      = 6.60d0
        n_a      = 11.89d0
        n_rho    = 1.035d0
        n_lambda = 36.1d0
        print *, 'music: medium = sea water'
    case ('rock', 'standard_rock', 'standard rock')
        n_z      = 11.00d0
        n_a      = 22.00d0
        n_rho    = 2.65d0
        n_lambda = 26.48d0
        print *, 'music: medium = standard rock'
    case default
        print *, 'initialize_music: unknown medium "', trim(medium), '"'
        print *, '  valid choices: water, seawater, rock'
        error stop
    end select

    print '(a,f6.3,a,f6.3,a,f6.2,a,f6.2)', &
        '  n_z=', n_z, '  n_a=', n_a, '  rho=', n_rho, '  X0=', n_lambda

    ! Energy-loss coefficients: columns 1-4 are stochastic process coefficients,
    ! column 5 is the relative radiative loss (bs), column 6 is ionisation (b0).
    print *, 'music: loading energy losses table...'
    open(unit=4, file='music-eloss-sr.dat', form='formatted', status='old')
    read(4, 99) emulo
 99 format(16(5E12.5/),E12.5/)
    close(4)
c    write(6, 99) emulo

    ! Inelastic double-differential cross-section table (bilinear in E, v, angle):
    !   ema(61)       log10(E_mu/GeV),  0.0 to  6.0, step 0.1
    !   va(30)        log10(v),          0.0 to -3.0, step 0.1
    !   anga(50)      log10(theta/rad),  0.4 to -4.5, step 0.1
    !   dang(50,30,61) cumulative angular probability at each (E,v) node
    print *, 'music: loading inel. scatt. cross section...'
    open(unit=3, file='music-double-diff-rock.dat', form='formatted', status='old')
    read(3, '(5e15.7)') ema, va, anga, dang
    close(3)

    ! j0 — index into v2 at threshold v = 10^(-3*0.05*(j0-1)).
    ! j0=61 → threshold at log10(v) = -3.0: losses above this are sampled
    ! stochastically; below it they enter the continuous-loss term.
    j0 = 61

    ! v2(i) = log10 of relative energy transfer, from 0.0 to -4.95, step -0.05
    do i = 1, 100
        v2(i) = 0.0 - 0.05*(i-1)
    end do

    ! emc(il) = log10 of muon energy [GeV], from -1.0 to 7.0, step 0.1
    do il = 1, 81
        emc(il) = -1.0 + 0.1*(il-1)
    end do

    ! Integral cross-section table
    print *, 'music: loading integral cross sections...'
    open(unit=1, file='music-cross-sections-sr.dat', form='formatted', status='old')
    read(1, '(4(71(20(5E14.6/)/)/)/)')  cs2
    close(unit=1)

    ! Normalise: cs1(j,k,i) = cs2(j,k,i) / cs2(j0,k,i)
    ! cs1 is used in vcs to sample the fractional energy transfer uniformly.
    do i = 1, 4
        do k = 1, 71
            do j = 1, j0
                if (cs2(j0,k,i) > 0.0) &
                    cs1(j,k,i) = cs2(j,k,i) / cs2(j0,k,i)
            end do
        end do
    end do

end subroutine initialize_music


! ----------------------------------------------------------------
! muon_transport
!
! Propagate a single muon through material by one step (to depth zf),
! updating position, direction, and timing in-place.
!
! Arguments (all intent(inout)):
!   x0,y0,z0     position [cm / n_rho]
!   cx0,cy0,cz0  direction cosines
!   Emuin0       muon energy [GeV]; set to 0 if muon stops (< 0.106 GeV)
!   depth0       depth along z-axis [cm]
!   tmu0         accumulated travel time [ns]; updated by this call
! ----------------------------------------------------------------
subroutine muon_transport(x0, y0, z0, cx0, cy0, cz0, Emuin0, depth0, tmu0)
    implicit none
    real(8), intent(inout) :: x0, y0, z0, cx0, cy0, cz0, Emuin0, depth0, tmu0

    real(8) :: n_a, n_z, n_rho, n_lambda
    common /rock2/ n_a, n_z, n_rho, n_lambda

    real(8), parameter :: pi = 3.141592653589793d0

    real(8) :: emuin, emu_f, depth, path_max, emu
    real(8) :: theta, phi, dr, x, y, z, t
    real(8) :: emu_0, theta0, phi0
    real(8) :: gamma, beta, vmu, tmu
    real(8) :: thetac, phic
    real(8) :: x1, y1, z1

    emuin = Emuin0

    emu_f    = 0.0d0
    depth    = depth0 * n_rho       ! convert depth [cm] to [g/cm^2]
    path_max = 1.0d6 * 100.0d0      ! maximum muon path [g/cm^2]
    emu      = emuin

    theta = 0.0d0;  phi = 0.0d0
    dr    = 0.0d0
    x     = 0.0d0;  y   = 0.0d0;  z = 0.0d0
    t     = 0.0d0

    emu_0  = emu
    theta0 = theta
    phi0   = phi

    call music(emu, depth, emu_f, x, y, z, t, theta, phi, dr, path_max)

    ! Compute travel time from relativistic kinematics (averaged over step)
    gamma = (emu_0 + emu_f) / 2.0d0 / 0.105655d0
    beta  = sqrt(gamma*gamma - 1.0d0) / gamma
    vmu   = beta * 29.9792458d0         ! [cm/ns]
    tmu   = dr / vmu / n_rho
    tmu0  = tmu0 + tmu

    ! Polar angle of the incoming direction
    thetac = acos(cz0)
    if (cx0 /= 0.0d0) then
        phic = atan(cy0 / cx0)
        if (cx0 < 0.0d0) phic = pi + phic
    else
        phic = pi / 2.0d0
        if (cy0 < 0.0d0) phic = pi * 1.5d0
    end if
    if (phic >  2.0d0*pi) phic = phic - 2.0d0*pi
    if (phic <  0.0d0)    phic = phic + 2.0d0*pi

    ! Transform displacement from muon frame to lab frame
    call coord_transform(x1, y1, z1, x, y, z, thetac, phic)
    x0 = x1/n_rho + x0
    y0 = y1/n_rho + y0
    z0 = z1/n_rho + z0

    theta0 = thetac
    phi0   = phic
    call angle_transform(theta0, phi0, theta, phi)

    cz0 = cos(theta0)
    cx0 = sin(theta0) * cos(phi0)
    cy0 = sin(theta0) * sin(phi0)

    Emuin0 = emu_f
    if (Emuin0 <= 0.106d0) Emuin0 = 0.0d0  ! muon has stopped

end subroutine muon_transport


! ----------------------------------------------------------------
! music
!
! Core Monte Carlo propagation engine. Transports a muon from its
! current position to depth zf (or until it stops), sampling
! stochastic interactions above the continuous-loss threshold v=10^-3.
!
! Input:
!   Emu    initial muon energy [GeV]
!   zf     observation depth in z [g/cm^2]
!   t0     initial accumulated path length [g/cm^2]
!   theta  polar angle with respect to z-axis [rad]
!   phi    azimuthal angle [rad]
!   t10    maximum allowed muon path [g/cm^2]
! Input/output:
!   x,y,z  muon position [g/cm^2]
! Output:
!   Emu_f  final muon energy [GeV]; = mmu if stopped
!   dr     path length travelled in this step [g/cm^2]
!   theta  updated polar angle [rad]
!   phi    updated azimuthal angle [rad]
!
! Common blocks used: /sig/, /coenlo/, /ang/, /rock2/
!
! Control flags (compile-time parameters):
!   ms_flag=1  enable multiple scattering (0 to disable)
!   d_flag=1   enable deflection at stochastic interactions (0 to disable)
! ----------------------------------------------------------------
subroutine music(Emu, zf, Emu_f, x, y, z, t0, theta, phi, dr, t10)
    implicit none
    real(8), intent(in)    :: zf, t0, t10
    real(8), intent(inout) :: Emu, x, y, z, theta, phi
    real(8), intent(out)   :: Emu_f, dr

    real(8) :: n_a, n_z, n_rho, n_lambda
    common /rock2/ n_a, n_z, n_rho, n_lambda

    real(8),  parameter :: mmu    = 0.105655d0
    integer,  parameter :: ms_flag = 1   ! set to 0 to disable multiple scattering
    integer,  parameter :: d_flag  = 1   ! set to 0 to disable stochastic deflections

    integer  :: j0
    real(4)  :: v2(100), cs1(100,71,4), cs2(100,71,4)
    real(4)  :: emulo(81,6), emc(81)
    real(4)  :: ema(61), va(30), anga(50), dang(50,30,61)
    common /sig/    v2, cs1, cs2, j0
    common /coenlo/ emulo, emc
    common /ang/    ema, va, anga, dang

    ! Local working variables
    real(4)  :: yfl         ! uniform random number (single precision for common-block interop)
    integer  :: err         ! error flag from multiple scattering
    integer  :: J, J1, J2   ! energy index into cross-section tables
    integer  :: ip          ! interaction process (1=bremss, 2=pair, 3=inel, 4=ion)
    integer  :: NM          ! step-type flag (1=to depth, 2=sub-1GeV; informational only)
    real(8)  :: em          ! persistent muon energy across iterations
    real(8)  :: em1         ! muon energy at start of current iteration
    real(8)  :: em2         ! energy after continuous losses to interaction point
    real(8)  :: em3         ! energy after continuous losses to observation depth
    real(8)  :: emf         ! rough estimate of final energy (for midpoint eval)
    real(8)  :: emi         ! midpoint energy for vdem re-evaluation
    real(8)  :: Demf        ! energy deposited in this stochastic interaction
    real(8)  :: t1          ! running accumulated path [g/cm^2]
    real(8)  :: tpath       ! path to next stochastic interaction [g/cm^2]
    real(8)  :: t2          ! path from current position to observation depth
    real(8)  :: t3          ! t1 + tpath (path total after this step)
    real(8)  :: t4          ! path to depth for sub-1-GeV case
    real(8)  :: z1          ! z-projection of tpath
    real(8)  :: z4          ! remaining z distance to depth (sub-1 GeV)
    real(8)  :: dz          ! z-projection of t2
    real(8)  :: b0          ! ionisation loss coefficient [GeV cm^2/g]
    real(8)  :: bs          ! relative radiative loss coefficient [cm^2/g]
    real(8)  :: FP, FP1     ! mean free paths at J and J-1 energy bins
    real(8)  :: css, css1   ! total cross sections for interaction-type sampling
    real(8)  :: al1, al2, al3  ! fractional cross sections for bremss, pair, inel
    real(8)  :: v1          ! fractional energy transfer to secondary
    real(8)  :: theta1, phi1   ! deflection angles at stochastic interaction
    real(8)  :: theta2, phi2   ! deflection angles from multiple scattering
    real(8)  :: deltax, deltay ! lateral displacement [g/cm^2]
    real(8)  :: deltar         ! path length correction from multiple scattering
    real(8)  :: dx1, dy1, dz1  ! coordinate increments in lab frame
    logical  :: normal_exit    ! .true. when depth is reached normally (dr = t1-t0)

    em          = emu
    t1          = t0
    err         = 0
    normal_exit = .true.

    ! ---- Main propagation loop ----
    ! Each iteration advances the muon by one mean free path (above 1 GeV)
    ! or to the observation depth (below 1 GeV).
    muon_loop: do

        em1 = em
        J   = int(dlog10(em1)*10.0d0) + 2
        if (J >= 71) J = 71
        J1  = J - 1

        if (em1 <= 1.0d0) then
            ! ---- Sub-1 GeV: skip stochastic sampling, only ionisation ----
            call vdem(em1, b0, bs)
            z4  = zf - z
            t4  = z4 / dabs(dcos(theta))
            emf = em1 - (em1*bs + b0)*t4
            if (emf <= mmu) emf = mmu
            emi = 10.0d0**((dlog10(emf) + dlog10(em1))/2.0d0)
            call vdem(emi, b0, bs)
            call vem(em1, b0, bs, t4, em3)
            if (em3 <= mmu) then
                dr          = t1 - t0 + em1/b0  ! range formula for stopping
                em1         = mmu
                normal_exit = .false.
                exit muon_loop
            end if
            theta2 = 0.0d0;  phi2   = 0.0d0
            deltax = 0.0d0;  deltay = 0.0d0;  deltar = 0.0d0
            dx1    = 0.0d0;  dy1    = 0.0d0;  dz1    = z4
            if (ms_flag == 1) &
                call multiple(t4, em1, em3, b0, bs, &
                              theta2, phi2, deltax, deltay, deltar, err)
            call coord_transform(dx1, dy1, dz1, deltax, deltay, z4, theta, phi)
            em1 = em3
            if (em1 <= mmu) then
                normal_exit = .false.
                exit muon_loop
            end if
            NM = 2
            z  = z  + dz1;  x = x + dx1;  y = y + dy1
            t1 = t1 + t4  + deltar
            if (theta2 > 0.0d0) call angle_transform(theta, phi, theta2, phi2)
            exit muon_loop  ! depth reached normally

        end if

        ! ---- Sample next stochastic interaction point ----
        ! Mean free path FP interpolated linearly in log10(E) between bins J1 and J.
        call random_number(yfl)
        FP  = 1.0d0 / (cs2(j0,J,1)  + cs2(j0,J,2)  + cs2(j0,J,3)  + cs2(j0,J,4))
        FP1 = 1.0d0 / (cs2(j0,J1,1) + cs2(j0,J1,2) + cs2(j0,J1,3) + cs2(j0,J1,4))
        FP  = (FP - FP1)/(emc(J+10) - emc(J1+10)) * (dlog10(em1) - emc(J1+10)) + FP1
        tpath = -FP * log(real(yfl, 8))   ! exponential waiting-time
        t3    = t1 + tpath
        z1    = tpath * dcos(theta)

        if (z + z1 > zf) then
            ! ---- Interaction would be beyond observation depth ----
            ! Propagate continuously from here to the depth surface.

            if (t3 > t10) then
                ! Muon path exceeds maximum allowed
                em1 = mmu;  dr = t10
                normal_exit = .false.
                exit muon_loop
            end if

            t2  = (zf - z) / dabs(dcos(theta))
            dz  = t2 * dcos(theta)

            ! Continuous energy loss to depth (two-point midpoint correction)
            call vdem(em1, b0, bs)
            emf = em1 - (em1*bs + b0)*t2
            if (emf <= mmu) emf = mmu
            emi = 10.0d0**((dlog10(emf) + dlog10(em1))/2.0d0)
            call vdem(emi, b0, bs)
            call vem(em1, b0, bs, t2, em3)
            if (em3 <= mmu) then
                dr          = t1 - t0 + em1/b0
                em1         = mmu
                normal_exit = .false.
                exit muon_loop
            end if

            theta2 = 0.0d0;  phi2   = 0.0d0
            deltax = 0.0d0;  deltay = 0.0d0;  deltar = 0.0d0
            dx1    = 0.0d0;  dy1    = 0.0d0;  dz1    = dz
            if (ms_flag == 1) &
                call multiple(t2, em1, em3, b0, bs, &
                              theta2, phi2, deltax, deltay, deltar, err)
            call coord_transform(dx1, dy1, dz1, deltax, deltay, dz, theta, phi)
            if (theta2 > 0.0d0) call angle_transform(theta, phi, theta2, phi2)

            em1 = em3
            if (em3 <= mmu) then
                normal_exit = .false.
                exit muon_loop
            end if
            NM    = 1
            tpath = tpath - t2 + deltar   ! remaining step after depth surface
            z     = z  + dz1;  x = x + dx1;  y = y + dy1
            t1    = t1 + t2   + deltar
            exit muon_loop  ! depth reached normally

        else
            ! ---- Interaction within remaining path ----
            ! Continuous losses and multiple scattering up to interaction point.

            call vdem(em1, b0, bs)
            emf = em1 - (em1*bs + b0)*tpath
            if (emf <= mmu) emf = mmu
            emi = 10.0d0**((dlog10(emf) + dlog10(em1))/2.0d0)
            call vdem(emi, b0, bs)
            call vem(em1, b0, bs, tpath, em2)
            if (em2 <= mmu) then
                dr          = t1 - t0 + em1/b0
                em1         = mmu
                normal_exit = .false.
                exit muon_loop
            end if

            theta2 = 0.0d0;  phi2   = 0.0d0
            deltax = 0.0d0;  deltay = 0.0d0;  deltar = 0.0d0
            dx1    = 0.0d0;  dy1    = 0.0d0;  dz1    = z1
            if (ms_flag == 1) &
                call multiple(tpath, em1, em2, b0, bs, &
                              theta2, phi2, deltax, deltay, deltar, err)
            call coord_transform(dx1, dy1, dz1, deltax, deltay, z1, theta, phi)

            em1 = em2
            if (em1 <= mmu) then
                normal_exit = .false.
                exit muon_loop
            end if
            t1 = t3 + deltar
            z  = z  + dz1;  x = x + dx1;  y = y + dy1
            if (theta2 > 0.0d0) call angle_transform(theta, phi, theta2, phi2)

            if (z > zf) exit muon_loop  ! depth reached normally

            ! Update persistent energy before potential cycle.
            em = em1
            if (em1 <= 1.0d0) cycle muon_loop  ! next iter handles sub-GeV

            ! ---- Sample interaction type ----
            call random_number(yfl)
            J  = int(dlog10(em1)*10.0d0) + 2
            J2 = int(dlog10(em1)*10.0d0 + 0.5d0) + 1
            if (J  >= 71) J  = 71
            if (J2 >= 71) J2 = 71
            J1 = J - 1

            css  = cs2(j0,J,1)  + cs2(j0,J,2)  + cs2(j0,J,3)  + cs2(j0,J,4)
            css1 = cs2(j0,J1,1) + cs2(j0,J1,2) + cs2(j0,J1,3) + cs2(j0,J1,4)
            css  = (css - css1) / (emc(J+10) - emc(J1+10)) * &
                   (dlog10(em1) - emc(J1+10)) + css1

            al1 = ((cs2(j0,J,1) - cs2(j0,J1,1)) / (emc(J+10) - emc(J1+10)) * &
                   (dlog10(em1)  - emc(J1+10))   + cs2(j0,J1,1)) / css
            al2 = ((cs2(j0,J,2) - cs2(j0,J1,2)) / (emc(J+10) - emc(J1+10)) * &
                   (dlog10(em1)  - emc(J1+10))   + cs2(j0,J1,2)) / css
            al3 = ((cs2(j0,J,3) - cs2(j0,J1,3)) / (emc(J+10) - emc(J1+10)) * &
                   (dlog10(em1)  - emc(J1+10))   + cs2(j0,J1,3)) / css

            if (yfl <= al1) then
                ip = 1                       ! bremsstrahlung
            else if (yfl <= al1 + al2) then
                ip = 2                       ! pair production
            else if (yfl <= al1 + al2 + al3) then
                ip = 3                       ! nuclear inelastic
            else
                ip = 4                       ! ionisation above threshold
            end if

            ! ---- Sample fractional energy transfer to secondary ----
            call random_number(yfl)
            call vcs(J2, yfl, ip, v1)

            em     = em1
            theta1 = 0.0d0;  phi1 = 0.0d0
            if (ip >= 1 .and. ip <= 3 .and. d_flag == 1) &
                call mu_scatt(ip, em, v1, theta1, phi1)

            Demf = v1 * em1
            em   = em  - Demf
            em1  = em
            if (theta1 > 0.0d0) call angle_transform(theta, phi, theta1, phi1)

            if (em1 <= mmu) then
                normal_exit = .false.
                exit muon_loop
            end if
            cycle muon_loop   ! next interaction

        end if

    end do muon_loop

    ! Normal exit: observation depth was reached
    if (normal_exit) then
        dr    = t1 - t0
        emu_f = em1
    end if
    ! Ensure emu_f reflects a stopped muon regardless of exit path
    if (em1 <= mmu) emu_f = mmu

end subroutine music


! ----------------------------------------------------------------
! vcs
!
! Sample the fractional energy transfer v = dE/E to a secondary
! particle, given the muon energy index j, a uniform random number
! yfl, and the interaction process ip.
!
! Uses the normalised cumulative cross-section cs1(100,71,4) from
! COMMON /sig/ together with the v2 grid, via linear interpolation.
!
! Input:  j   energy bin index (1..71)
!         yfl uniform random variate in [0,1]
!         ip  process index (1=bremss, 2=pair, 3=inel, 4=ion)
! Output: v1  sampled fractional energy transfer
! ----------------------------------------------------------------
subroutine vcs(j, yfl, ip, v1)
    implicit none
    integer,  intent(in)  :: j, ip
    real(4),  intent(in)  :: yfl
    real(8),  intent(out) :: v1

    integer  :: j0
    real(4)  :: v2(100), cs1(100,71,4), cs2(100,71,4)
    common /sig/ v2, cs1, cs2, j0

    integer :: i
    real(8) :: v3

    ! Walk up the normalised CDF until cs1(i) >= yfl
    i = 1
    do
        i = i + 1
        if (i > 100) exit   ! safety clamp: yfl = 1.0 exactly
        if (yfl <= cs1(i,j,ip)) exit
    end do
    i = min(i, 100)

    ! Linear interpolation in log-v space
    v3 = (v2(i) - v2(i-1)) / (cs1(i,j,ip) - cs1(i-1,j,ip)) * &
         (yfl    - cs1(i-1,j,ip)) + v2(i-1)
    v1 = 10.0d0**v3

end subroutine vcs


! ----------------------------------------------------------------
! vdem
!
! Compute continuous energy-loss coefficients at muon energy EM:
!   b0 [GeV cm^2/g]  — ionisation loss rate (Bethe-Bloch)
!   bs [cm^2/g]      — fractional radiative loss rate for v < 10^-3
!
! Both are interpolated from the tabulated array in COMMON /coenlo/.
! Below 1 GeV (log10(E) <= 0), bs is forced to zero and b0 is
! clamped to tabulated boundary values.
! ----------------------------------------------------------------
subroutine vdem(EM, B0, BS)
    implicit none
    real(8), intent(in)  :: EM
    real(8), intent(out) :: B0, BS

    real(4) :: f(81,6), em0(81)
    common /coenlo/ f, em0

    real(8) :: e1, f1, f2
    integer :: i

    e1 = dlog10(em)

    ! Find the bracketing index: first i where em0(i) >= e1
    i = 2
    do while (i <= 81 .and. e1 > em0(i))
        i = i + 1
    end do
    if (i > 81) i = 81

    ! Radiative loss coefficient (column 5); zero below 1 GeV
    f1 = (f(i,5) - f(i-1,5)) / (em0(i) - em0(i-1)) * (e1 - em0(i-1)) + f(i-1,5)
    bs = f1
    if (bs <= 0.0d0 .or. e1 <= 0.0d0) bs = 0.0d0

    ! Ionisation loss coefficient (column 6) with boundary clamping
    f2 = (f(i,6) - f(i-1,6)) / (em0(i) - em0(i-1)) * (e1 - em0(i-1)) + f(i-1,6)
    b0 = f2
    if (e1 >= -0.1d0 .and. e1 < 0.0d0) b0 = f(10,6)  ! near-threshold clamp
    if (e1 <= -0.9d0)                    b0 = f(2,6)   ! low-energy clamp
    if (b0 <= 0.0d0) b0 = 0.0d0

end subroutine vdem


! ----------------------------------------------------------------
! vem
!
! Integrate the continuous energy-loss equation over path length T:
!   dE/dt = -(b0 + bs*E)
!
! Exact solution:  E(T) = (E0 + b0/bs)*exp(-bs*T) - b0/bs
! Linear approx:   E(T) = E0 - b0*T - bs*T*E0      (when bs*T << 1)
!
! Input:  EM0  initial energy [GeV]
!         B0   ionisation coefficient
!         BS   radiative coefficient
!         T    path length [g/cm^2]
! Output: EM   final energy [GeV]
! ----------------------------------------------------------------
subroutine vem(EM0, B0, BS, T, EM)
    implicit none
    real(8), intent(in)  :: EM0, B0, BS, T
    real(8), intent(out) :: EM

    if (bs*T > 1.0d-6) then
        em = em0 * dexp(-bs*T) - b0/bs * (1.0d0 - dexp(-bs*T))
    else
        em = em0 - b0*T - bs*T*em0
    end if

end subroutine vem


! ----------------------------------------------------------------
! defl
!
! Sample the muon deflection angle from nuclear inelastic scattering
! (ip=3 only; returns immediately for other process types).
!
! Uses 3D lookup table dang(50,30,61) — cumulative angular probability
! at each (log10(theta), log10(v), log10(E)) node — via bilinear
! interpolation in (log10(v), log10(E)).
!
! The same random variate yfl is used across all four interpolation
! nodes to sample a consistent quantile from the marginal distributions.
!
! Output: theta  polar deflection angle [rad] (0 if below -4.5 threshold)
!         phi    azimuthal angle, uniform on [0, 2pi]
! ----------------------------------------------------------------
subroutine defl(em, ip, v, theta, phi)
    implicit none
    real(8), intent(in)    :: em, v
    integer, intent(in)    :: ip
    real(8), intent(inout) :: theta, phi

    real(4) :: ema(61), va(30), anga(50), dang(50,30,61)
    common /ang/ ema, va, anga, dang

    real(8), parameter :: pi = 3.141592653589793d0

    real(4)  :: yfl
    real(8)  :: ang(2,2)           ! sampled log10(theta) at 4 (v,E) nodes
    real(8)  :: eml, vl            ! log10(E), log10(v)
    real(8)  :: ang1, ang2, angt   ! interpolated log10(theta) values
    integer  :: k1, k2, m1, m2    ! bracketing indices in ema, va grids
    integer  :: kk, mm, k, m, j

    if (ip /= 3) return

    call random_number(yfl)
    if (em <= 1.0d0) return

    eml = dlog10(em)
    k1  = int(10.0d0*eml) + 1
    k1  = max(1, min(60, k1))
    k2  = max(2, min(61, k1+1))

    vl  = dlog10(v)
    m1  = int(-vl*10.0d0)
    m1  = max(1, min(29, m1))
    m2  = max(2, min(30, m1+1))

    ! Sample angle at each of the four bracketing (m,k) nodes using the
    ! same quantile yfl so that bilinear interpolation stays consistent.
    do kk = 1, 2
        k = merge(k1, k2, kk == 1)
        do mm = 1, 2
            m = merge(m1, m2, mm == 1)
            j = 1
            ang_find: do
                j = j + 1
                if (j > 49) then
                    ang(mm,kk) = -5.0d0   ! below minimum angle in table
                    exit ang_find
                end if
                if (yfl <= dang(j,m,k)) then
                    ang(mm,kk) = (anga(j) - anga(j-1)) / &
                                 (dang(j,m,k) - dang(j-1,m,k)) * &
                                 (yfl - dang(j-1,m,k)) + anga(j-1)
                    exit ang_find
                end if
            end do ang_find
        end do
    end do

    ! Bilinear interpolation: first in v-axis, then in E-axis
    ang1 = (ang(2,1) - ang(1,1)) / (va(m2) - va(m1)) * (vl - va(m1)) + ang(1,1)
    ang2 = (ang(2,2) - ang(1,2)) / (va(m2) - va(m1)) * (vl - va(m1)) + ang(1,2)
    angt = (ang2 - ang1) / (ema(k2) - ema(k1)) * (eml - ema(k1)) + ang1

    ! -4.5 is the minimum log10(theta) in the table; treat as zero deflection
    if (angt > -4.5d0) then
        theta = 10.0d0**angt
    else
        theta = 0.0d0
    end if

    call random_number(yfl)
    phi = 2.0d0 * pi * yfl
    if (theta == 0.0d0) phi = 0.0d0

end subroutine defl


! ----------------------------------------------------------------
! mu_scatt
!
! Sample the muon deflection angle from a single stochastic interaction.
!
! Input:
!   jp    process: 1=bremsstrahlung, 2=pair production, 3=nuclear inel.
!   E     muon energy before interaction [GeV]
!   v     fractional energy transfer dE/E
! Output:
!   theta  polar deflection angle [rad], always >= 0
!   phi    azimuthal angle [rad], uniform on [0, 2pi]
!
! Bremsstrahlung: Van Ginneken parametrisation
! Pair production: parametric formula
! Nuclear inelastic: delegated to defl()
! ----------------------------------------------------------------
subroutine mu_scatt(jp, E, v, theta, phi)
    implicit none
    integer, intent(in)  :: jp
    real(8), intent(in)  :: E, v
    real(8), intent(out) :: theta, phi

    real(8) :: n_a, n_z, n_rho, n_lambda
    common /rock2/ n_a, n_z, n_rho, n_lambda

    real(8), parameter :: pi   = 3.141592653589793d0
    real(8), parameter :: mmu  = 0.105655d0          ! muon mass [GeV]
    real(8), parameter :: mel  = 0.511d-3             ! electron mass [GeV]
    ! conv = 180/pi [deg/rad]: kept for reference; angles are in radians
    real(8), parameter :: conv = 57.29577951d0

    real(4)  :: yfl
    real(8)  :: the_mean
    real(8)  :: ak1, ak2, ak3, ak4, ak5, an, the1, the2, the3, at
    real(8)  :: a_min

    if (jp == 1) then
        ! ---- Bremsstrahlung: Van Ginneken ----
        ak1 = 0.092d0 * E**(-1.0d0/3.0d0)
        ak3 = 0.22d0  * E**(-0.92d0)
        ak4 = 0.26d0  * E**(-0.91d0)
        an  = 0.81d0  * dsqrt(E) / (dsqrt(E) + 1.8d0)
        ak2 = 0.052d0 / E * n_z**(-0.25d0)
        at  = min(ak1*dsqrt(v), ak2)
        the1 = max(at, ak3*v)
        the2 = ak4 * v**(1.0d0 + an) * (1.0d0 - v)**(-an)
        ak5  = ak4 * 0.5d0**(1.0d0 + an) * 0.5d0**(-an) / 0.5d0**(-0.5d0)
        the3 = ak5 / dsqrt(1.0d0 - v)
        if (v <= 0.5d0) then
            the_mean = the1
        else if (the2 < 0.2d0) then
            the_mean = the2
        else
            the_mean = the3
        end if
        the_mean = the_mean * the_mean  ! now the_mean = <theta^2>

        ! Exponential distribution in theta^2; rejection for theta > pi
        do
            call random_number(yfl)
            theta = dsqrt(-the_mean * log(real(yfl, 8)))
            if (theta <= pi) exit
        end do
        call random_number(yfl)
        phi = 2.0d0 * pi * yfl

    else if (jp == 2) then
        ! ---- Pair production ----
        a_min = min(8.9d-4*sqrt(sqrt(real(v,8)))*(1.0d0 + 1.5d-5*E) + &
                    0.032d0*v/(v + 1.0d0), 0.1d0)
        the_mean = (2.3d0 + dlog(E)) / E / (1.0d0-v) * &
                   (v - 2.0d0*mel/E)**2 / v / v * a_min
        the_mean = the_mean * the_mean

        do
            call random_number(yfl)
            theta = dsqrt(-the_mean * log(real(yfl, 8)))
            if (theta <= pi) exit
        end do
        call random_number(yfl)
        phi = 2.0d0 * pi * yfl

    else if (jp == 3) then
        ! ---- Nuclear inelastic: use angular distribution table ----
        call defl(E, jp, v, theta, phi)

    else
        print *, 'Invalid jp argument in mu_scatt: ', jp
        stop
    end if

end subroutine mu_scatt


! ----------------------------------------------------------------
! multiple
!
! Simulate multiple Coulomb scattering over path length x0.
!
! Based on:
!   L. Highland, Nucl. Instr. Methods 129, 497 (1975)
!   P. Lipari & T. Stanev, Phys. Rev. D 44, 3543 (1991)
!
! The scattering widths sigma_theta and sigma_x are computed from an
! integral along the muon path. Correlated pairs (theta_x, x) and
! (theta_y, y) are drawn using corset/corgen.
!
! Input:
!   x0      path length [g/cm^2]
!   einmu   initial muon energy [GeV]
!   eoutmu  final muon energy [GeV]
!   alpha   ionisation loss coefficient (b0) [GeV cm^2/g]
!   beta    radiative loss coefficient  (bs) [cm^2/g]
! Output:
!   theta   polar deflection [rad]
!   phi     azimuthal deflection [rad]
!   deltax  lateral displacement x [g/cm^2]
!   deltay  lateral displacement y [g/cm^2]
!   deltar  approximate path-length increase [g/cm^2]
!   err     0=ok, 1=energy below mmu, 2=integral <= 0
! ----------------------------------------------------------------
subroutine multiple(x0, einmu, eoutmu, alpha, beta, &
                    theta, phi, deltax, deltay, deltar, err)
    implicit none
    real(8), intent(in)  :: x0, einmu, eoutmu, alpha, beta
    real(8), intent(out) :: theta, phi, deltax, deltay, deltar
    integer, intent(out) :: err

    real(8) :: n_a, n_z, n_rho, n_lambda
    common /rock2/ n_a, n_z, n_rho, n_lambda

    real(8), parameter :: pi   = 3.141592653589793d0
    real(8), parameter :: mmu  = 0.105655d0
    real(8), parameter :: ems  = 0.015d0   ! Highland constant [GeV] (tuned 1997)

    real(8) :: p1, g1, bet1, p2, g2, bet2
    real(8) :: sigt, sigma, sigmat, sigmax
    real(8) :: sz, sz1, step, zz, ro
    real(8) :: thetax, thetay
    real(8) :: vv(2,2), cc(2,2), vcx(2)   ! real(8) required by corset/corgen

    err = 0
    if (einmu <= mmu .or. eoutmu <= mmu) then
        err = 1
        return
    end if

    ! Relativistic momentum × velocity at entry and exit
    p1   = dsqrt(einmu*einmu   - mmu*mmu)
    g1   = einmu  / mmu
    bet1 = dsqrt(g1*g1 - 1.0d0) / g1
    p2   = dsqrt(eoutmu*eoutmu - mmu*mmu)
    g2   = eoutmu / mmu
    bet2 = dsqrt(g2*g2 - 1.0d0) / g2

    ! Integral of 1/(p*beta)^2 along path — Highland formula
    if (beta > 0.0d0) then
        sigt = 1.0d0/alpha * (1.0d0/(p2*bet2) - 1.0d0/(p1*bet1)) + &
               beta/alpha/alpha * &
               (dlog(dabs(1.0d0 + alpha/beta/p1/bet1)) - &
                dlog(dabs(1.0d0 + alpha/beta/p2/bet2)))
    else
        sigt = 1.0d0/alpha * (1.0d0/(p2*bet2) - 1.0d0/(p1*bet1))
    end if
    if (sigt <= 0.0d0) sigt = 1.0d0/(p1*bet1)**2  ! fallback

    sigma  = sigt * ems*ems / n_lambda
    sigmat = dsqrt(sigma)                ! sigma(theta_x) = sigma(theta_y)

    ! Override for test programs where einmu == eoutmu
    if (einmu == eoutmu) &
        sigmat = dsqrt(ems*ems * x0 / einmu**2 / n_lambda)

    ! Numerical integral of (z/p(z)/beta(z))^2 dz for sigma_x
    sz   = 0.0d0
    step = x0 / 100.0d0
    zz   = 0.0d0
    path_integral: do while (zz < x0)
        zz = zz + step/2.0d0
        if (beta > 0.0d0) then
            sz1 = p1*bet1 * dexp(-beta*zz) - alpha/beta*(1.0d0 - dexp(-beta*zz))
        else
            sz1 = p1*bet1 - alpha*zz
        end if
        if (sz1 <= mmu) then
            sz = 0.0d0
            exit path_integral
        end if
        sz = sz + zz*zz * step / (sz1*sz1)
        zz = zz + step/2.0d0
    end do path_integral

    if (sz <= 0.0d0) sz = x0**3 / 3.0d0 / (p1*bet1)**2  ! fallback
    if (sz <= 0.0d0) then
        err = 2
        return
    end if

    sigmax = dsqrt(sz * ems*ems / n_lambda)  ! sigma(x) = sigma(y)
    ro     = sqrt(3.0d0) / 2.0d0            ! correlation between theta and x

    ! Draw two correlated (theta, x) pairs — one for each transverse plane
    vv(1,1) = sigmat * sigmat
    vv(2,2) = sigmax * sigmax
    vv(1,2) = ro * sigmat * sigmax
    vv(2,1) = vv(1,2)
    call corset(vv, cc, 2)

    call corgen(cc, vcx, 2)
    thetax = vcx(1)
    deltax = vcx(2)

    call corgen(cc, vcx, 2)
    thetay = vcx(1)
    deltay = vcx(2)

    ! Convert component deflections to polar (theta, phi)
    theta = datan(dsqrt(dtan(thetax)**2 + dtan(thetay)**2))

    if (dsin(thetax) /= 0.0d0) then
        phi = datan(dtan(thetay) / dtan(thetax))
        if (dtan(thetax) < 0.0d0) phi = pi + phi
    else
        phi = pi / 2.0d0
        if (dtan(thetay) < 0.0d0) phi = 1.5d0 * pi
    end if
    if (phi >  2.0d0*pi) phi = phi - 2.0d0*pi
    if (phi <  0.0d0)    phi = phi + 2.0d0*pi

    deltar = sigmat * sigmat * x0   ! approximate extra path from deflection

end subroutine multiple


! ----------------------------------------------------------------
! coord_transform
!
! Rotate displacement vector (x1,y1,z1) from the muon reference frame
! (z'-axis along muon direction theta1,phi1) to the lab frame.
!
! The rotation matrix R(theta,phi) is constructed so that its third
! row points along the muon direction.
!
! Input:  x1,y1,z1   vector in muon frame
!         theta1,phi1 muon polar/azimuthal angles in lab frame [rad]
! Output: x,y,z       rotated vector in lab frame
! ----------------------------------------------------------------
subroutine coord_transform(x, y, z, x1, y1, z1, theta1, phi1)
    implicit none
    real(8), intent(out) :: x, y, z
    real(8), intent(in)  :: x1, y1, z1, theta1, phi1

    real(8) :: theta, phi
    real(8) :: t11, t12, t13, t21, t22, t23, t31, t32, t33

    theta = theta1
    phi   = phi1

    t11 = dcos(theta)*dcos(phi)*dcos(phi) + dsin(phi)*dsin(phi)
    t21 = dcos(theta)*dcos(phi)*dsin(phi) - dcos(phi)*dsin(phi)
    t31 = -dsin(theta)*dcos(phi)
    t12 = dcos(theta)*dcos(phi)*dsin(phi) - dcos(phi)*dsin(phi)
    t22 = dcos(theta)*dsin(phi)*dsin(phi) + dcos(phi)*dcos(phi)
    t32 = -dsin(theta)*dsin(phi)
    t13 = dsin(theta)*dcos(phi)
    t23 = dsin(theta)*dsin(phi)
    t33 = dcos(theta)

    x = t11*x1 + t12*y1 + t13*z1
    y = t21*x1 + t22*y1 + t23*z1
    z = t31*x1 + t32*y1 + t33*z1

end subroutine coord_transform


! ----------------------------------------------------------------
! angle_transform
!
! Compose a deflection (theta1,phi1) in the muon local frame with the
! muon's current direction (theta,phi) in the lab frame to produce
! the new lab-frame direction.
!
! The deflection vector is built in spherical coordinates, rotated by
! R(theta,phi) to the lab frame, then converted back to (theta,phi).
!
! Input/output: theta,phi  current muon direction [rad]; updated in-place
! Input:        theta1,phi1 deflection in muon frame [rad]
! ----------------------------------------------------------------
subroutine angle_transform(theta, phi, theta1, phi1)
    implicit none
    real(8), intent(inout) :: theta, phi
    real(8), intent(in)    :: theta1, phi1

    real(8), parameter :: pi = 3.141592653589793d0

    real(8) :: t11, t12, t13, t21, t22, t23, t31, t32, t33
    real(8) :: x2, y2, z2, x21, y21, z21

    ! Zero deflection: trivial case
    if (theta == 0.0d0) then
        theta = theta1
        phi   = phi1
        return
    end if

    t11 = dcos(theta)*dcos(phi)*dcos(phi) + dsin(phi)*dsin(phi)
    t21 = dcos(theta)*dcos(phi)*dsin(phi) - dcos(phi)*dsin(phi)
    t31 = -dsin(theta)*dcos(phi)
    t12 = dcos(theta)*dcos(phi)*dsin(phi) - dcos(phi)*dsin(phi)
    t22 = dcos(theta)*dsin(phi)*dsin(phi) + dcos(phi)*dcos(phi)
    t32 = -dsin(theta)*dsin(phi)
    t13 = dsin(theta)*dcos(phi)
    t23 = dsin(theta)*dsin(phi)
    t33 = dcos(theta)

    ! Direction cosines of deflection vector in muon frame
    x2 = dsin(theta1)*dcos(phi1)
    y2 = dsin(theta1)*dsin(phi1)
    z2 = dcos(theta1)

    ! Rotate to lab frame
    x21 = t11*x2 + t12*y2 + t13*z2
    y21 = t21*x2 + t22*y2 + t23*z2
    z21 = t31*x2 + t32*y2 + t33*z2

    ! Guard against floating-point drift outside [-1, 1] before acos
    z21 = max(-1.0d0, min(1.0d0, z21))
    theta = dacos(z21)

    if (x21 /= 0.0d0) then
        phi = datan(y21 / x21)
        if (x21 < 0.0d0) phi = pi + phi
    else
        phi = pi / 2.0d0
        if (y21 < 0.0d0) phi = 1.5d0 * pi
    end if
    if (phi >  2.0d0*pi) phi = phi - 2.0d0*pi
    if (phi <  0.0d0)    phi = phi + 2.0d0*pi

end subroutine angle_transform


! ----------------------------------------------------------------
! str_lower  —  convert a character string to lowercase in-place
! ----------------------------------------------------------------
subroutine str_lower(s)
    implicit none
    character(len=*), intent(inout) :: s
    integer :: i, c
    do i = 1, len(s)
        c = iachar(s(i:i))
        if (c >= iachar('A') .and. c <= iachar('Z')) &
            s(i:i) = achar(c + 32)
    end do
end subroutine str_lower
