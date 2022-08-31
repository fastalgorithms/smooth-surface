CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      SUBROUTINE FTOPHYS(MEXPF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1                      NTHMAX,MEXPPHYS,FEXPE,FEXPO)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      COMPLEX *16 MEXPF(1)
      COMPLEX *16 MEXPPHYS(1),IMA,CTMP
ccc      COMPLEX *16 FEXPS(40,100,40)
      COMPLEX *16 FEXPE(1)
      COMPLEX *16 FEXPO(1)
      REAL *8     RLAMS(NLAMBS)
      REAL *8     ALPHAS(0:300)
cccc      INTEGER *4  NLAMBS,NUMFOUR(NLAMBS),NUMPHYS(NLAMBS),NTHMAX
      INTEGER  NLAMBS,NUMFOUR(NLAMBS),NUMPHYS(NLAMBS),NTHMAX
      DATA IMA/(0.0D0,1.0D0)/
C***********************************************************************
C
C     This subroutine evaluates the Fourier expansion of the
C     exponential moment function M(\lambda,\alpha) at equispaced
C     nodes.
C
C     On INPUT:
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
C     NLAMBS:        number of discretization pts. in lambda integral
C     RLAMS(NLAMBS): discretization points in lambda integral.
C     NUMFOUR(j):   number of Fourier modes in the expansion
C                      of the function M(\lambda_j,\alpha)
C     NTHMAX =      max_j NUMFOUR(j)
C     FEXPE =      precomputed array of exponentials needed for
C                  Fourier series evaluation
C     FEXPO =      precomputed array of exponentials needed for
C                  Fourier series evaluation
C
C     On OUTPUT:
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
C------------------------------------------------------------
      DONE=1.0d0
C
C
      PI=DATAN(DONE)*4
      NFTOT = 0
      NPTOT  = 0
      NEXTE = 1
      NEXTO = 1
      DO 2000 I=1,NLAMBS
        DO 1200 IVAL=1,NUMPHYS(I)
           MEXPPHYS(NPTOT+IVAL) = MEXPF(NFTOT+1)
           DO 200 MM = 2,NUMFOUR(I),2
              RT1 = DIMAG(FEXPE(NEXTE))*DREAL(MEXPF(NFTOT+MM))
              RT2 = DREAL(FEXPE(NEXTE))*DIMAG(MEXPF(NFTOT+MM))
	      RTMP = 2*(RT1+RT2)
              NEXTE = NEXTE + 1
              MEXPPHYS(NPTOT+IVAL) = MEXPPHYS(NPTOT+IVAL) +
     1                DCMPLX(0.0D0,RTMP)
200        CONTINUE
           DO 400 MM = 3,NUMFOUR(I),2
              RT1 = DREAL(FEXPO(NEXTO))*DREAL(MEXPF(NFTOT+MM))
              RT2 = DIMAG(FEXPO(NEXTO))*DIMAG(MEXPF(NFTOT+MM))
	      RTMP = 2*(RT1-RT2)
              NEXTO = NEXTO + 1
              MEXPPHYS(NPTOT+IVAL) = MEXPPHYS(NPTOT+IVAL) +
     1                DCMPLX(RTMP,0.0D0)
400        CONTINUE
 1200   CONTINUE
        NFTOT = NFTOT+NUMFOUR(I)
        NPTOT = NPTOT+NUMPHYS(I)
 2000 CONTINUE
      RETURN
      END
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      SUBROUTINE EXPEVAL(XP,YP,ZP,CENTER,RLAMS,WHTS,NLAMBS,NUMTETS,
     1                   MEXPPHYS,CINT)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
cccc      INTEGER *4 NUMTETS(1),NLAMBS
      INTEGER  NUMTETS(1),NLAMBS
      COMPLEX *16 ALPHINTV,CINT
      COMPLEX *16 MEXPPHYS(1)
      DIMENSION CENTER(3)
      DIMENSION RLAMS(1),WHTS(1)
      DIMENSION ALPHAS(200)
C***********************************************************************
C
C     This subroutine evaluates the integral 
C       
C       1/2pi \int_0^{\infty} e^{-z * \lambda}                (1) 
C                                                                
C       \int_0^{2 \pi} e^{i rlam  (x cos(\alpha)+ y sin(\alpha) ) } * 
C        
C       M(lambda,alpha)  d\alpha  d\lambda
C
C       via Norman's quadratures.
C
C     ON INPUT:
C
C     XP,YP,ZP:      Cartesian coordinates of evaluation point
C     CENTER:        Cartesian coordinates of center (implicitely
C                    the origin in the above formula).
C     NLAMBS:        number of discretization pts. in lambda integral
C     RLAMS(NLAMBS): Norman's nodes for the lambda integral.
C     WHTS(NLAMBS):  Norman's weights for the lambda integral.
C     NUMTETS(NLAMBS): number of nodes in alpha integral for 
C                      corresponding lambda value.
C     MEXPPHYS(*):   Discrete values of the moment function 
C                    M(\lambda,\alpha), ordered as follows.
C
C         MEXPPHYS(1),...,MEXPPHYS(NUMPHYS(1)) = M(\lambda_1,0),..., 
C              M(\lambda_1, 2*PI*(NUMPHYS(1)-1)/NUMPHYS(1)).
C         MEXPPHYS(NUMPHYS(1)+1),...,MEXPPHYS(NUMPHYS(2)) = 
C              M(\lambda_2,0),...,
C                  M(\lambda_2, 2*PI*(NUMPHYS(2)-1)/NUMPHYS(2)).
C         etc.
C
C     ON OUTPUT:
C
C     CINT:         Value of integral.
C
C
C------------------------------------------------------------
c
      DONE=1.0d0
      PI=DATAN(DONE)*4
      CINT=0
      NTOT = 1
C
C----- (outer) lambda integral
C
      X = XP - CENTER(1)
      Y = YP - CENTER(2)
      Z = ZP - CENTER(3)
      DO 2000 I=1,NLAMBS
	NALPHA = NUMTETS(I)
        HALPHA=2*PI/NALPHA
        DO 1400 J=1,NALPHA
          ALPHAS(J)=(J-1)*HALPHA
 1400   CONTINUE
C
C----- (inner) alpha integral
C
        CALL ALPHINT(RLAMS(I),X,Y,NALPHA,ALPHAS,
     1               MEXPPHYS(NTOT),ALPHINTV)
        CINT=CINT+WHTS(I)*ALPHINTV*DEXP(-RLAMS(I)*Z)
        NTOT = NTOT+NUMTETS(I)
 2000 CONTINUE
        CINT=CINT/(2*PI)
        RETURN
        END
c
c
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
        SUBROUTINE ALPHINT(RLAM,X,Y,NALPHA,ALPHAS,MEXPPHYS,ALPHINTV)
C***********************************************************************
        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION ALPHAS(1)
        COMPLEX *16 ALPHINTV,CD,IMA,CDD
        COMPLEX *16 MEXPPHYS(1)
        DATA IMA/(0.0D0,1.0D0)/
C***********************************************************************
C
C       This subroutine evaluates the integral 
C
C       \int_0^{2 \pi} e^{i*rlam*(x*cos(\alpha)+y*sin(\alpha))} *
C           M(rlam,alpha)  d\alpha                                (1)
C
C       ON INPUT:
C
C       RLAM:    given value of lambda variable.
C       X,Y:     Cartesian coordinates of target point. 
C       NALPHA:  the number of nodes in the trapezoidal discretization
C                of the integral (1)
C       ALPHAS:  NALPHA equispaced nodes on the interval [0, 2 \pi]
C       MEXPPHYS(NALPHA): discrete values of function M(rlam,alpha)
C                         at equispaced nodes.
C
C       ON OUTPUT:
C
C       ALPHINTV: NALPHA-point trapezoidal approximation to integral (1)
C  
C----------------------------------------------------------------------
      CD=0
      DO 1200 I=1,NALPHA
         CDD=CDEXP(IMA*RLAM*(X*DCOS(ALPHAS(I))+Y*DSIN(ALPHAS(I)) )  )
         CDD=CDD*MEXPPHYS(I)
         CD=CD+CDD
 1200 CONTINUE
      ALPHINTV=CD*(ALPHAS(2)-ALPHAS(1))
      RETURN
      END
C
