program test_ecomug
    ! Quick verification of the ecomug_mod routines.
    !
    ! Checks:
    !   1. dalt_flux reproduces known surface values
    !   2. dmu_shallow(p,x,0) == dalt_flux(p,x)
    !   3. gen_horiz / gen_vert produce physically sensible mean (p, cos_theta)
    !
    ! Compile:
    !   gfortran -O2 -o test_ecomug ecomug.f90 test_ecomug.f90
    use ecomug_mod
    implicit none

    integer, parameter :: NWARM  = 5000
    integer, parameter :: NEVT   = 100000
    real(8), parameter :: DEPTH_H = 188.0d0   ! top surface depth [m]
    real(8), parameter :: DEPTH_V = 198.0d0   ! side surface midpoint [m]

    real(8), allocatable :: pmu(:), cth(:), ang(:)
    real(8) :: sum_p, sum_c, sum_a, mean_p, mean_c, mean_a
    integer :: i

    integer :: iseed(8) = [72889, 64746, 45, 567, 72828, 27282, 83838, 28293]
    call random_seed(put=iseed)

    ! ---- Flux checks ------------------------------------------------
    write(*, '(a)') '=== Flux function checks ==='
    write(*, '(a,es12.5)') '  dalt_flux(1,   1)           = ', dalt_flux(1.0d0,   1.0d0)
    write(*, '(a,es12.5)') '  dalt_flux(10,  1)           = ', dalt_flux(10.0d0,  1.0d0)
    write(*, '(a,es12.5)') '  dalt_flux(100, 1)           = ', dalt_flux(100.0d0, 1.0d0)
    write(*, '(a,es12.5)') '  dalt_flux(10,  0.5)         = ', dalt_flux(10.0d0,  0.5d0)
    write(*, *)
    write(*, '(a,es12.5)') '  dmu_shallow(10, 1, 0) [=above?] = ', dmu_shallow(10.0d0, 1.0d0, 0.0d0)
    write(*, '(a,es12.5)') '  dmu_shallow(50, 1,188)          = ', dmu_shallow(50.0d0, 1.0d0, 188.0d0)
    write(*, '(a,es12.5)') '  dmu_shallow(50, 0.5,188)        = ', dmu_shallow(50.0d0, 0.5d0, 188.0d0)
    write(*, *)

    allocate(pmu(NEVT), cth(NEVT), ang(NEVT))

    ! ---- Horizontal surface -----------------------------------------
    write(*, '(a,f6.1,a,i0,a)') '=== Horizontal surface  depth=', DEPTH_H, &
        ' m   N=', NEVT, ' ==='
    call gen_horiz(DEPTH_H, NWARM, NEVT, pmu, cth, ang)
    sum_p = 0.0d0;  sum_c = 0.0d0;  sum_a = 0.0d0
    do i = 1, NEVT
        sum_p = sum_p + pmu(i)
        sum_c = sum_c + cth(i)
        sum_a = sum_a + ang(i)
    end do
    mean_p = sum_p / NEVT
    mean_c = sum_c / NEVT
    mean_a = sum_a / NEVT
    write(*, '(a,f8.3,a)') '  <pmu>       = ', mean_p, ' GeV/c'
    write(*, '(a,f8.4)')    '  <cos theta> = ', mean_c
    write(*, '(a,f8.4,a)')  '  <phi>       = ', mean_a, ' rad  (expect ~pi)'
    write(*, *)

    ! ---- Vertical surface -------------------------------------------
    write(*, '(a,f6.1,a,i0,a)') '=== Vertical surface    depth=', DEPTH_V, &
        ' m   N=', NEVT, ' ==='
    call gen_vert(DEPTH_V, NWARM, NEVT, pmu, cth, ang)
    sum_p = 0.0d0;  sum_c = 0.0d0;  sum_a = 0.0d0
    do i = 1, NEVT
        sum_p = sum_p + pmu(i)
        sum_c = sum_c + cth(i)
        sum_a = sum_a + ang(i)
    end do
    mean_p = sum_p / NEVT
    mean_c = sum_c / NEVT
    mean_a = sum_a / NEVT
    write(*, '(a,f8.3,a)') '  <pmu>       = ', mean_p, ' GeV/c'
    write(*, '(a,f8.4)')    '  <cos theta> = ', mean_c
    write(*, '(a,f8.5,a)')  '  <phi_rel>   = ', mean_a, ' rad  (expect ~0 by symmetry)'
    write(*, *)

    deallocate(pmu, cth, ang)

end program test_ecomug
