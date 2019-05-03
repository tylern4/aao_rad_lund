!------------------------------------------------------------------------
!EV The following two routines concern natural cubic spline interpolation
!   of a function of two independent variables.
!------------------------------------------------------------------------

SUBROUTINE splie2(x2a, ya, m, n, y2a)
    !
    ! Given an m by n tabulated function ya(1:m,1:n), and tabulated independent
    ! variable x2a(1:n), this routine constructs 1-dimensional natural cubic
    ! splines of the 2nd dimension (1:n) of ya and returns the 2nd-derivatives in
    ! the array y2a(1:m,1:n).  From "Numerical Recipes" (FORTRAN, 1992 Ed.),
    ! Ch. 3.6, p. 121.  Call this routine only once, before a sequence of calls to
    ! SPLIN2.  The process is of order  m x n.
    IMPLICIT NONE
    ! Passed variables:
    INTEGER m, n
    REAL x2a(n), ya(m, n), y2a(m, n)
    ! Local variables:
    INTEGER i, j, NN
    REAL bound
    PARAMETER (NN = 100)     ! Maximum expected value of n
    REAL ytmp(NN), y2tmp(NN)
    ! Boundary condition threshold for natural cubic spline:
    SAVE bound
    DATA bound/1.e30/
    !
    do i = 1, m
        do j = 1, n
            ytmp(j) = ya(i, j)
        enddo
        ! Natural spline:
        call spline(x2a, ytmp, n, bound, bound, y2tmp)
        do j = 1, n
            y2a(i, j) = y2tmp(j)
        enddo
    enddo
    !
    return
END