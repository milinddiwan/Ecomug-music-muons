subroutine randn(n, arr)
    implicit none
    integer, intent(in)  :: n
    real(8), intent(out) :: arr(n)

    real(8), parameter :: TWO_PI = 6.2831853071795864769d0
    real(8) :: u1, u2
    integer :: i

    call random_seed()

    do i = 1, n - 1, 2
        call random_number(u1)
        call random_number(u2)
        ! Avoid log(0)
        u1 = max(u1, 1.0d-300)
        arr(i)   = sqrt(-2.0d0 * log(u1)) * cos(TWO_PI * u2)
        arr(i+1) = sqrt(-2.0d0 * log(u1)) * sin(TWO_PI * u2)
    end do

    ! If n is odd, generate one extra pair and keep one value
    if (mod(n, 2) /= 0) then
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, 1.0d-300)
        arr(n) = sqrt(-2.0d0 * log(u1)) * cos(TWO_PI * u2)
    end if

end subroutine randn
