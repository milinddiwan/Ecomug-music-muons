module ecomug_mod
    ! ECoMUG-based cosmic muon flux and Monte Carlo generator.
    !
    ! Provides sea-level and shallow-depth flux functions, plus
    ! Metropolis-Hastings generators for horizontal and vertical surfaces.
    !
    ! Reference: D. Pagano et al., NIM A 1014 (2021) 165732
    !
    ! Units: momentum p [GeV/c], cos(zenith) x in [0,1], depth in metres,
    !        flux in muons m^-2 sr^-1 s^-1 (GeV/c)^-1.
    implicit none
    private
    public :: dalt_flux, dmu_shallow, gen_horiz, gen_vert, &
              gen_horiz_init, gen_horiz_step

    real(8), parameter :: PI        = 3.141592653589793d0
    ! Average ionisation loss averaged over 10 and 100 GeV (apar values from notebook)
    real(8), parameter :: AAVE      = 2.305d0   ! MeV / (g/cm^2)
    real(8), parameter :: RHO_WATER = 1.0d0     ! g/cm^3
    ! Metropolis proposal widths (symmetric uniform box: Delta drawn from [-w/2, +w/2])
    real(8), parameter :: DP_STEP   = 10.0d0    ! GeV/c
    real(8), parameter :: DCT_STEP  = 0.10d0    ! dimensionless
    ! Hard generation boundaries
    real(8), parameter :: PMU_LO = 1.0d0,  PMU_HI = 100.0d0
    real(8), parameter :: CTH_LO = 0.0d0,  CTH_HI = 1.0d0

contains

    ! ------------------------------------------------------------------
    pure function dalt_flux(pmu, x) result(f)
        ! ECoMUG sea-level differential muon flux (Pagano et al. 2021).
        ! pmu : momentum [GeV/c],  x = cos(zenith), 0 <= x <= 1
        real(8), intent(in) :: pmu, x
        real(8) :: f, npo
        npo = max(0.1d0, 2.856d0 - 0.655d0 * log(pmu))
        f   = 1600.0d0 * (pmu + 2.68d0)**(-3.175d0) * pmu**0.279d0 * x**npo
    end function dalt_flux

    ! ------------------------------------------------------------------
    pure function dmu_shallow(pmu, x, depth_m) result(f)
        ! Attenuated flux at vertical depth depth_m [m] (water equivalent).
        !
        ! Propagates straight-line energy loss backwards: the muon observed at
        ! depth with momentum pmu needed energy E_surf = pmu + a*rho*L at the
        ! surface, where L is the slant path length.  Then evaluates the
        ! sea-level flux at E_surf.  Valid for depths of a few hundred mwe.
        !
        ! pmu     : muon momentum at depth [GeV/c]  (p ~ E for relativistic muons)
        ! x       : cos(zenith angle)
        ! depth_m : vertical depth [m]
        real(8), intent(in) :: pmu, x, depth_m
        real(8) :: f, slant_cm, eloss_gev, esurf
        if (x <= 0.0d0) then
            f = 0.0d0;  return
        end if
        slant_cm  = 100.0d0 * depth_m / x           ! slant path length [cm]
        eloss_gev = AAVE * RHO_WATER * slant_cm * 1.0d-3  ! MeV -> GeV
        esurf     = pmu + eloss_gev
        f         = dalt_flux(esurf, x)
    end function dmu_shallow

    ! ------------------------------------------------------------------
    ! Internal unnormalised PDF kernels
    ! ------------------------------------------------------------------

    pure function horpdf_val(pmu, x, depth_m) result(f)
        ! PDF kernel proportional to dmu_shallow * cos(theta) for a
        ! horizontal surface; the cos(theta) = x factor is the flux
        ! projection onto the surface normal.
        real(8), intent(in) :: pmu, x, depth_m
        real(8) :: f
        if (pmu < PMU_LO .or. pmu > PMU_HI .or. &
            x   < CTH_LO .or. x   > CTH_HI) then
            f = 0.0d0
        else
            f = dmu_shallow(pmu, x, depth_m) * x
        end if
    end function horpdf_val

    pure function vertpdf_val(pmu, x, depth_m) result(f)
        ! PDF kernel proportional to dmu_shallow * sin(theta) for a
        ! vertical surface; the sin(theta) = sqrt(1-x^2) factor is the
        ! flux projection onto the surface normal (azimuthal integral over
        ! the front hemisphere gives factor 2, absorbed into normalisation).
        real(8), intent(in) :: pmu, x, depth_m
        real(8) :: f
        if (pmu < PMU_LO .or. pmu > PMU_HI .or. &
            x   < CTH_LO .or. x   > CTH_HI) then
            f = 0.0d0
        else
            f = dmu_shallow(pmu, x, depth_m) * sqrt(1.0d0 - x*x)
        end if
    end function vertpdf_val

    ! ------------------------------------------------------------------
    subroutine metro_step(cur_p, cur_x, cur_f, pdf_id, depth_m)
        ! One Metropolis-Hastings step with a symmetric uniform-box proposal.
        ! pdf_id = 1 -> horpdf_val,  2 -> vertpdf_val
        real(8), intent(inout) :: cur_p, cur_x, cur_f
        integer, intent(in)    :: pdf_id
        real(8), intent(in)    :: depth_m
        real(8) :: dp, dx, new_p, new_x, new_f, u

        call random_number(dp);  dp = (dp - 0.5d0) * DP_STEP
        call random_number(dx);  dx = (dx - 0.5d0) * DCT_STEP
        new_p = cur_p + dp
        new_x = min(cur_x + dx, CTH_HI)   ! clamp: avoid sqrt(1-x^2) for x > 1

        if (pdf_id == 1) then
            new_f = horpdf_val(new_p, new_x, depth_m)
        else
            new_f = vertpdf_val(new_p, new_x, depth_m)
        end if

        if (new_f > 0.0d0) then
            if (cur_f == 0.0d0) then
                cur_p = new_p;  cur_x = new_x;  cur_f = new_f
            else
                call random_number(u)
                if (u < min(1.0d0, new_f / cur_f)) then
                    cur_p = new_p;  cur_x = new_x;  cur_f = new_f
                end if
            end if
        end if
    end subroutine metro_step

    ! ------------------------------------------------------------------
    subroutine gen_horiz(depth_m, nwarmup, n, pmu_arr, cth_arr, phi_arr)
        ! Metropolis generator for muons incident on a horizontal surface.
        !
        ! depth_m  - vertical depth [m]
        ! nwarmup  - MCMC burn-in steps before recording (suggest >= 5000)
        ! n        - number of muons to generate
        ! pmu_arr  - muon momenta [GeV/c]
        ! cth_arr  - cos(zenith angle), distributed as dmu_shallow * cos(theta)
        ! phi_arr  - azimuthal angle [rad], uniform in [0, 2*pi)
        real(8), intent(in)  :: depth_m
        integer, intent(in)  :: nwarmup, n
        real(8), intent(out) :: pmu_arr(n), cth_arr(n), phi_arr(n)
        real(8) :: cur_p, cur_x, cur_f, u
        integer :: i

        cur_p = 1.0d0;  cur_x = 0.999d0
        cur_f = horpdf_val(cur_p, cur_x, depth_m)
        do i = 1, nwarmup
            call metro_step(cur_p, cur_x, cur_f, 1, depth_m)
        end do
        do i = 1, n
            call metro_step(cur_p, cur_x, cur_f, 1, depth_m)
            pmu_arr(i) = cur_p
            cth_arr(i) = cur_x
            call random_number(u)
            phi_arr(i) = u * 2.0d0 * PI
        end do
    end subroutine gen_horiz

    ! ------------------------------------------------------------------
    subroutine gen_vert(depth_m, nwarmup, n, pmu_arr, cth_arr, phi_rel_arr)
        ! Metropolis generator for muons incident on a vertical surface.
        !
        ! depth_m     - vertical depth of the surface midpoint [m]
        ! nwarmup     - MCMC burn-in steps
        ! n           - number of muons
        ! pmu_arr     - muon momenta [GeV/c]
        ! cth_arr     - cos(zenith angle), distributed as dmu_shallow * sin(theta)
        ! phi_rel_arr - azimuthal angle relative to the inward surface normal [rad],
        !               sampled proportional to cos(phi_rel) over (-pi/2, +pi/2)
        !               via CDF inversion: phi = arcsin(2*u - 1), u uniform
        real(8), intent(in)  :: depth_m
        integer, intent(in)  :: nwarmup, n
        real(8), intent(out) :: pmu_arr(n), cth_arr(n), phi_rel_arr(n)
        real(8) :: cur_p, cur_x, cur_f, u
        integer :: i

        cur_p = 1.0d0;  cur_x = 0.999d0
        cur_f = vertpdf_val(cur_p, cur_x, depth_m)
        do i = 1, nwarmup
            call metro_step(cur_p, cur_x, cur_f, 2, depth_m)
        end do
        do i = 1, n
            call metro_step(cur_p, cur_x, cur_f, 2, depth_m)
            pmu_arr(i) = cur_p
            cth_arr(i) = cur_x
            call random_number(u)
            phi_rel_arr(i) = asin(2.0d0 * u - 1.0d0)
        end do
    end subroutine gen_vert

    ! ------------------------------------------------------------------
    subroutine gen_horiz_init(depth_m, nwarmup, s_p, s_x, s_f)
        ! Warm up the horizontal-surface Metropolis chain and return its state.
        ! Call once; then drive event-by-event with gen_horiz_step.
        real(8), intent(in)  :: depth_m
        integer, intent(in)  :: nwarmup
        real(8), intent(out) :: s_p, s_x, s_f
        integer :: i
        s_p = 1.0d0;  s_x = 0.999d0
        s_f = horpdf_val(s_p, s_x, depth_m)
        do i = 1, nwarmup
            call metro_step(s_p, s_x, s_f, 1, depth_m)
        end do
    end subroutine gen_horiz_init

    subroutine gen_horiz_step(depth_m, s_p, s_x, s_f, pmu, cth, phi)
        ! Advance the chain by one step and return one muon's kinematics.
        ! s_p, s_x, s_f are updated in place (persistent chain state).
        real(8), intent(in)    :: depth_m
        real(8), intent(inout) :: s_p, s_x, s_f
        real(8), intent(out)   :: pmu, cth, phi
        real(8) :: u
        call metro_step(s_p, s_x, s_f, 1, depth_m)
        pmu = s_p;  cth = s_x
        call random_number(u)
        phi = u * 2.0d0 * PI
    end subroutine gen_horiz_step

end module ecomug_mod
