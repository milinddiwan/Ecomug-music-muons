subroutine corset(V, C, NP)
    ! Correlated Gaussian random number setup.
    ! Computes C = Cholesky lower-triangular square root of covariance matrix V.
    ! Call corgen after this to generate correlated samples.
    implicit none
    integer,  intent(in)    :: NP
    real(8),  intent(in)    :: V(NP, NP)
    real(8),  intent(out)   :: C(NP, NP)

    integer :: i, j
    real(8) :: ck

    C = 0.0d0

    do j = 1, NP
        ! Diagonal term
        ck = sum(C(j, 1:j-1)**2)
        C(j, j) = sqrt(abs(V(j, j) - ck))

        ! Off-diagonal terms
        do i = j + 1, NP
            ck = sum(C(i, 1:j-1) * C(j, 1:j-1))
            C(i, j) = (V(i, j) - ck) / C(j, j)
        end do
    end do

end subroutine corset
