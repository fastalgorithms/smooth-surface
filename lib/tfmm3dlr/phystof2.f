C***********************************************************************
      SUBROUTINE PHYSTOF(MEXPF,NLAMBS,NUMFOUR,NUMPHYS,
     1                      MEXPPHYS,FEXPBACK)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      COMPLEX *16 MEXPF(1)
      COMPLEX *16 MEXPPHYS(1),IMA
      COMPLEX *16 FEXPBACK(1)
      REAL *8     ALPHAS(0:100)
cccc      INTEGER *4  NLAMBS,NUMFOUR(NLAMBS),NUMPHYS(NLAMBS),NTHMAX
      INTEGER   NLAMBS,NUMFOUR(NLAMBS),NUMPHYS(NLAMBS),NTHMAX
      DATA IMA/(0.0D0,1.0D0)/
C***********************************************************************
C
C     This subroutine converts the discretized exponential moment function
C     into its Fourier expansion.
C
C     On INPUT:
C
C     MEXPPHYS(*):  Discrete values of the moment function 
C                   M(\lambda,\alpha), ordered as follows.
C
C         MEXPPHYS(1),...,MEXPPHYS(NUMPHYS(1)) = M(\lambda_1,0),..., 
C              M(\lambda_1, 2*PI*(NUMPHYS(1)-1)/NUMPHYS(1)).
C         MEXPPHYS(NUMPHYS(1)+1),...,MEXPPHYS(NUMPHYS(2)) = 
C              M(\lambda_2,0),...,
C                  M(\lambda_2, 2*PI*(NUMPHYS(2)-1)/NUMPHYS(2)).
C         etc.
C
C     NLAMBS:        number of discretization pts. in lambda integral
C     NUMFOUR(j):   number of Fourier modes in the expansion
C                      of the function M(\lambda_j,\alpha)
C     NTHMAX =      max_j NUMFOUR(j)
C
C     On OUTPUT:
C
C     MEXPF(*):     Fourier coefficients of the function 
C                   MEXP(lambda,alpha) for discrete lambda values. 
C                   They are ordered as follows:
C
C               MEXPF(1,...,NUMFOUR(1)) = Fourier modes for lambda_1
C               MEXPF(NUMFOUR(1)+1,...,NUMFOUR(2)) = Fourier modes
C                                              for lambda_2
C               etc.
C
C------------------------------------------------------------
      DONE=1.0d0
C
C
      PI=DATAN(DONE)*4
      NFTOT = 0
      NPTOT  = 0
      NEXT  = 1
      DO 2000 I=1,NLAMBS
        NALPHA = NUMPHYS(I)
        HALPHA=2*PI/NALPHA
        DO 200 J=1,NALPHA
          ALPHAS(J)=(J-1)*HALPHA
200     CONTINUE
        MEXPF(NFTOT+1) = 0.0D0
        DO 400 IVAL=1,NALPHA
           MEXPF(NFTOT+1) = MEXPF(NFTOT+1) + MEXPPHYS(NPTOT+IVAL) 
400     CONTINUE
        MEXPF(NFTOT+1) = MEXPF(NFTOT+1)/NALPHA
        DO 800 MM = 2,NUMFOUR(I)
           MEXPF(NFTOT+MM) = 0.0D0
           DO 600 IVAL=1,NALPHA
              MEXPF(NFTOT+MM) = MEXPF(NFTOT+MM) +
     1          FEXPBACK(NEXT)*MEXPPHYS(NPTOT+IVAL)
                NEXT = NEXT+1
600        CONTINUE
           MEXPF(NFTOT+MM) = MEXPF(NFTOT+MM)/NALPHA
800     CONTINUE
        NFTOT = NFTOT+NUMFOUR(I)
        NPTOT = NPTOT+NUMPHYS(I)
 2000 CONTINUE
      RETURN
      END
C
