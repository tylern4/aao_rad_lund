!------------------------------------------------------------------------

SUBROUTINE splin2(x1a, x2a, ya, y2a, m, n, x1, x2, y)
    !
    ! Given an m by n tabulated function ya(1:m,1:n), and tabulated independent
    ! variables x1a(1:m), x2a(1:n), and tabulated 2nd-derivatives y2a(1:m,1:n) as
    ! produced by SPLIE2; and given a desired interpolating point x1, x2; this
    ! routine returns an interpolated function value y by 2-dimensional natural
    ! cubic spline interpolation.  From "Numerical Recipes" (FORTRAN, 1992 Ed.),
    ! Ch. 3.6, p. 121.  Call this routine sequentially after one call to SPLIE2.
    ! The process is of order  m x logn + m + logm + 1 ; i.e., roughly, of order
    ! m x logn.  Therefore, the optimum choice is n > m.
    IMPLICIT NONE
    ! Passed variables:
    INTEGER m, n
    REAL  x1a(m), x2a(n), ya(m, n), y2a(m, n), x1, x2, y
    ! Local variables:
    INTEGER i, j, NN
    REAL bound
    PARAMETER (NN = 100)           ! Maximum expected value of m, n
    REAL ytmp(NN), yytmp(NN), y2tmp(NN)
    ! Boundary condition threshold for natural cubic spline:
    SAVE bound
    DATA bound/1.e30/
    !
    do i = 1, m
        do j = 1, n
            ytmp(j) = ya(i, j)
            y2tmp(j) = y2a(i, j)
        enddo
        ! Perform n evaluations of the splines of the (1:n)-dimension constructed by
        ! SPLIE2, using the 1-dimensional spline evaluator SPLINT.
        call splint(x2a, ytmp, y2tmp, n, x2, yytmp(i))
    enddo
    ! Construct the 1-dimensional natural spline of the (1:m)-dimension and
    ! evaluate it.
    call spline(x1a, yytmp, m, bound, bound, y2tmp)
    call splint(x1a, yytmp, y2tmp, m, x1, y)
    !
    return
END