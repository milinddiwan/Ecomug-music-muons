program test_corgen
    implicit none
    integer, parameter :: N = 1000
    real(8) :: samples(N)
    real(8) :: V(2,2), C(2,2), X(2)

    call randn(N, samples)
    write(*, *) 'First 5:', samples(1:5)

    V(1,1) = 100.0d0
    V(2,1) = -50.0d0
    V(1,2) = -50.0d0
    V(2,2) = 100.0d0

    call corset(V, C, 2)
    call corgen(C, X, 2)

    write(*, *) 'Correlated randn:', X(1), X(2)

end program test_corgen
