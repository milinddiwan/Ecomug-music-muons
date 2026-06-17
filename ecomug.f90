module ecomug2_mod
    ! ECoMUG-based cosmic muon flux and Monte Carlo generator.
    !
    ! Provides sea-level and shallow-depth flux functions, plus
    ! Metropolis-Hastings generators for horizontal and vertical surfaces.
    !
    ! Reference: D. Pagano et al., NIM A 1014 (2021) 165732
    !
    ! Units: momentum p [GeV/c], cos(zenith) x in [0,1], depth in metres,
    !        flux in muons m^-2 sr^-1 s^-1 (GeV/c)^-1.
    !
    ! Changes vs ecomug_mod (v1):
    !   - metro_step uses a multiplicative log10(p) proposal instead of
    !     additive DP_STEP; the Metropolis acceptance ratio is corrected
    !     for the resulting Jacobian factor (new_p / cur_p).
    !   - NTHIN MH steps are discarded between each recorded sample to
    !     reduce chain autocorrelation.
    implicit none
    private
    public :: dalt_flux, dmu_shallow, gen_horiz, gen_vert, &
              gen_horiz_init, gen_horiz_step

    real(8), parameter :: PI        = 3.141592653589793d0
    real(8), parameter :: AAVE      = 2.305d0   ! MeV / (g/cm^2)
    real(8), parameter :: RHO_WATER = 1.0d0     ! g/cm^3

    ! Metropolis proposal widths
    ! DLGP_STEP : sigma of normal proposal in log10(p); 0.15 -> typical step x/÷ 1.41
    ! DCT_STEP  : sigma of normal proposal in cos(theta)
    real(8), parameter :: DLGP_STEP = 0.15d0
    real(8), parameter :: DCT_STEP  = 0.10d0

    ! Thinning: number of MH steps taken per recorded sample
    integer,  parameter :: NTHIN    = 10

    ! Hard generation boundaries
    real(8), parameter :: PMU_LO = 1.0d0,  PMU_HI = 2000.0d0
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
        real(8), intent(in) :: pmu, x, depth_m
        real(8) :: f, slant_cm, eloss_gev, esurf
        if (x <= 0.0d0) then
            f = 0.0d0;  return
        end if
        slant_cm  = 100.0d0 * depth_m / x
        eloss_gev = AAVE * RHO_WATER * slant_cm * 1.0d-3
        esurf     = pmu + eloss_gev
        f         = dalt_flux(esurf, x)
    end function dmu_shallow

    ! ------------------------------------------------------------------
    ! Internal unnormalised PDF kernels
    ! ------------------------------------------------------------------
    pure function horpdf_val(pmu, x, depth_m) result(f)
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
        ! One Metropolis-Hastings step with a log-space proposal for p.
        !
        ! Proposal: new_p = cur_p * 10^(rn1 * DLGP_STEP)  [multiplicative]
        !           new_x = cur_x + rn2 * DCT_STEP         [additive]
        !
        ! Because the proposal is symmetric in log10(p) but we evaluate the
        ! target density in p-space, detailed balance requires multiplying
        ! the acceptance ratio by the Jacobian factor new_p / cur_p.
        !
        ! pdf_id = 1 -> horpdf_val,  2 -> vertpdf_val
        real(8), intent(inout) :: cur_p, cur_x, cur_f
        integer, intent(in)    :: pdf_id
        real(8), intent(in)    :: depth_m
        real(8) :: dlgp, dx, new_p, new_x, new_f, ratio, u
        real(8) :: rn(2)

        call randn(2, rn)
        dlgp  = rn(1) * DLGP_STEP          ! step in log10(p)
        dx    = rn(2) * DCT_STEP
        new_p = cur_p * 10.0d0**dlgp       ! multiplicative; always positive
        new_x = cur_x + dx

        if (pdf_id == 1) then
            new_f = horpdf_val(new_p, new_x, depth_m)
        else
            new_f = vertpdf_val(new_p, new_x, depth_m)
        end if

        if (new_f > 0.0d0) then
            if (cur_f == 0.0d0) then
                cur_p = new_p;  cur_x = new_x;  cur_f = new_f
            else
                ! Jacobian correction: ratio *= new_p / cur_p
                ratio = (new_f * new_p) / (cur_f * cur_p)
                call random_number(u)
                if (u < min(1.0d0, ratio)) then
                    cur_p = new_p;  cur_x = new_x;  cur_f = new_f
                end if
            end if
        end if
    end subroutine metro_step

    ! ------------------------------------------------------------------
    subroutine gen_horiz(depth_m, nwarmup, n, pmu_arr, cth_arr, phi_arr)
        ! Metropolis generator for muons incident on a horizontal surface.
        real(8), intent(in)  :: depth_m
        integer, intent(in)  :: nwarmup, n
        real(8), intent(out) :: pmu_arr(n), cth_arr(n), phi_arr(n)
        real(8) :: cur_p, cur_x, cur_f, u
        integer :: i, j

        cur_p = 1.0d0;  cur_x = 0.999d0
        cur_f = horpdf_val(cur_p, cur_x, depth_m)

        do i = 1, nwarmup
            call metro_step(cur_p, cur_x, cur_f, 1, depth_m)
        end do

        do i = 1, n
            do j = 1, NTHIN
                call metro_step(cur_p, cur_x, cur_f, 1, depth_m)
            end do
            pmu_arr(i) = cur_p
            cth_arr(i) = cur_x
            call random_number(u)
            phi_arr(i) = u * 2.0d0 * PI
        end do
    end subroutine gen_horiz

    ! ------------------------------------------------------------------
    subroutine gen_vert(depth_m, nwarmup, n, pmu_arr, cth_arr, phi_rel_arr)
        ! Metropolis generator for muons incident on a vertical surface.
        real(8), intent(in)  :: depth_m
        integer, intent(in)  :: nwarmup, n
        real(8), intent(out) :: pmu_arr(n), cth_arr(n), phi_rel_arr(n)
        real(8) :: cur_p, cur_x, cur_f, u
        integer :: i, j

        cur_p = 1.0d0;  cur_x = 0.999d0
        cur_f = vertpdf_val(cur_p, cur_x, depth_m)

        do i = 1, nwarmup
            call metro_step(cur_p, cur_x, cur_f, 2, depth_m)
        end do

        do i = 1, n
            do j = 1, NTHIN
                call metro_step(cur_p, cur_x, cur_f, 2, depth_m)
            end do
            pmu_arr(i) = cur_p
            cth_arr(i) = cur_x
            call random_number(u)
            phi_rel_arr(i) = asin(2.0d0 * u - 1.0d0)
        end do
    end subroutine gen_vert

    ! ------------------------------------------------------------------
    subroutine gen_horiz_init(depth_m, nwarmup, s_p, s_x, s_f)
        ! Warm up the horizontal-surface Metropolis chain and return its state.
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
        ! Advance the chain by NTHIN steps and return one muon's kinematics.
        ! s_p, s_x, s_f are updated in place (persistent chain state).
        real(8), intent(in)    :: depth_m
        real(8), intent(inout) :: s_p, s_x, s_f
        real(8), intent(out)   :: pmu, cth, phi
        real(8) :: u
        integer :: i

        do i = 1, NTHIN
            call metro_step(s_p, s_x, s_f, 1, depth_m)
        end do
        pmu = s_p;  cth = s_x
        call random_number(u)
        phi = u * 2.0d0 * PI
    end subroutine gen_horiz_step

end module ecomug2_mod
