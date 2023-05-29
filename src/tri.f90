!--------------------------------------------------------------------------

! Tridiagonal solver (uses Thomas Algorithm)

    SUBROUTINE TRI(n,a,d,c,b,x)

    USE precisions

    INTEGER(int_p), INTENT(in) :: n
    INTEGER(int_p) :: i
    REAL(real_p) :: xmult
    REAL(real_p), DIMENSION(n) :: a,b,c,d,x


    DO i = 2,n
        xmult = a(i)/d(i-1)
        d(i) = d(i) - xmult*c(i-1)
        b(i) = b(i) - xmult*b(i-1)
    ENDDO

    x(n) = b(n)/d(n)
    DO i = n-1,1,-1
        x(i) = (b(i)-c(i)*x(i+1))/d(i)
    ENDDO

    END SUBROUTINE TRI