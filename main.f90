PROGRAM DRIVER

  USE, NON_INTRINSIC :: PRECISIONS
  USE, NON_INTRINSIC :: CONVOLUTION, ONLY: ACORR

  REAL(FPP), DIMENSION(10) :: X, A

  !---------------------------------------------------------------------------------------------------------------------------------

  X = [1.8, -4.1, 3.2, 0.8, 1.4, -2.3, -0.1, 2.8, 1.1, -0.9]

  CALL ACORR(X, A)

  OPEN(1, FILE = 'debug.txt')

  DO I = 1, SIZE(X)
    WRITE(1, *) X(I), A(I)
  ENDDO

  CLOSE(1)


END PROGRAM DRIVER
