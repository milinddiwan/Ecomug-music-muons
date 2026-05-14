subroutine corgen(C, X, NP)
    ! Generates NP correlated Gaussian random numbers with covariance V = C*C'.
    ! Call corset first to compute C from V.
    implicit none
    integer, intent(in)  :: NP
    real(8), intent(in)  :: C(NP, NP)
    real(8), intent(out) :: X(NP)

    integer, parameter :: NMAX = 100
    real(8) :: Z(NMAX)
    integer :: i

    if (NP > NMAX) then
        write(*, '(a,i5,a,i5)') &
            'ERROR IN CORGEN. VECTOR LENGTH NP=', NP, &
            ', BUT MAXIMUM ALLOWED IS', NMAX
        return
    end if

    call randn(NP, Z(1:NP))

    do i = 1, NP
        X(i) = dot_product(C(i, 1:i), Z(1:i))
    end do

end subroutine corgen
