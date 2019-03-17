C***********************************************************************
      SUBROUTINE ROTGEN(NTERMS,CARRAY,RDPI2,RDMPI2,RDSQ3,RDMSQ3,DC)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 PI
      REAL *8 CARRAY(0:4*NTERMS,0:4*NTERMS)
      REAL *8 DC(0:4*NTERMS,0:4*NTERMS)
      REAL *8 RDPI2(0:NTERMS,0:NTERMS,-NTERMS:NTERMS) 
      REAL *8 RDMPI2(0:NTERMS,0:NTERMS,-NTERMS:NTERMS) 
      REAL *8 RDSQ3(0:NTERMS,0:NTERMS,-NTERMS:NTERMS) 
      REAL *8 RDMSQ3(0:NTERMS,0:NTERMS,-NTERMS:NTERMS) 
cccc      INTEGER *4 NTERMS
      INTEGER NTERMS
c
c----- call initialization routines
c
      CALL PRINI(6,13)
      PI = 4*DATAN(1.0d0)
c
      CALL BNLCFT(CARRAY,DC,4*NTERMS)
      THETA = PI/2.0d0
      CALL FSTRTN(NTERMS,RDPI2,DC,THETA)
      THETA = -PI/2.0d0
      CALL FSTRTN(NTERMS,RDMPI2,DC,THETA)
      THETA = DACOS(DSQRT(3.0D0)/3.0D0)
      call prin2(' theta = *',theta,1)
      CALL FSTRTN(NTERMS,RDSQ3,DC,THETA)
      THETA = DACOS(-DSQRT(3.0D0)/3.0D0)
      call prin2(' theta = *',theta,1)
      CALL FSTRTN(NTERMS,RDMSQ3,DC,THETA)
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE ZFLIP(NTERMS,MPOLE,MROTATE)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
      COMPLEX *16 MROTATE(0:NTERMS,0:NTERMS)
C
      DO 200 L=0,NTERMS,2
         DO 100 M=0,L
            MROTATE(L,M)=DCONJG(MPOLE(L,M))
100      CONTINUE
200   CONTINUE
      DO 400 L=1,NTERMS,2
         DO 300 M=0,L
            MROTATE(L,M)=-DCONJG(MPOLE(L,M))
300      CONTINUE
400   CONTINUE
ccc      CALL PRINF(' MROTATE is *',NTERMS,0)
ccc      CALL PRINM(MROTATE,NTERMS)
      RETURN
      END
C***********************************************************************
      SUBROUTINE PRINOUT(MPOLE,LL,NTERMS)
C***********************************************************************
c
c     printing routine displays multipole moments one degree (l)
c     at a time , from m = 0 to m = l  .
c
C***********************************************************************
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
cccc      INTEGER *4 NTERMS
      INTEGER NTERMS
C
      DO 100 L = 0,LL
	 WRITE(6,1000)(MPOLE(L,M),M=0,LL)
	 WRITE(13,1000)(MPOLE(L,M),M=0,LL)
ccc	 WRITE(6,1001)
ccc	 WRITE(13,1001)
100   CONTINUE
1000  FORMAT(6D12.5)
1001  FORMAT(/)
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE BNLCFT(C, SQC, NTERMS)
C***********************************************************************
C
C     Usage:
C           Computes the binomial coefficients C_NTERMS^n, where
C           n=0,1,2,...,NTERMS.
C     Input:
C           NTERMS: an integer indicates the number we are going
C                 to choose from.
C     Output:
C           C:    an array consists of the binomial coefficients.
C           SQC:   an array consists of the square root of the
C                 binormial coefficients.
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NTERMS,N,M
      DOUBLE PRECISION C(0:NTERMS,0:NTERMS)
      DOUBLE PRECISION SQC(0:NTERMS,0:NTERMS)
C
      DO 100 N=0,NTERMS
            C(N,0)=1.0D0
            SQC(N,0)=1.0D0
100   CONTINUE
      DO 300 M=1,NTERMS
         C(M,M)=1.0D0
         SQC(M,M)=1.0D0
         DO 200 N=M+1,NTERMS
            C(N,M)=C(N-1,M)+C(N-1,M-1)
            SQC(N,M)=DSQRT(C(N,M))
200      CONTINUE
300   CONTINUE
C
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE FSTRTN(NTERMS,D,SQC,THETA)
C***********************************************************************
C
C     Usage:
C           Implement the fast version of rotation matrices from
C           the recurrences formulas.
C     Input:
C           NTERMS: an integer indicates the dimension of D.
C           SQC:    an array contains the square root of the
C                   binormial coefficients.
C           THETA:  the rotate angle about the Y-axis.
C     Output:
C           D:      an array which contains the rotation matrix.
C                   Note: only half of D are evaluated, the other
C                   half can be obtained by using the symmetricity.
C
C***********************************************************************
C
C---------------Declares of FSTRTN--------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NTERMS
C
ccc     PARAMETER (NTERMS=8)
ccc     DOUBLE PRECISION D(0:NTERMS,-NTERMS:NTERMS,-NTERMS:NTERMS)
ccc     DOUBLE PRECISION C(0:4*NTERMS, 0:4*NTERMS)
C
      DOUBLE PRECISION D(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      DOUBLE PRECISION SQC(0:4*NTERMS, 0:4*NTERMS)
      DOUBLE PRECISION CTHETA, STHETA, HSTHTA
      DOUBLE PRECISION CTHTAP, CTHTAN, THETA, PRECIS
C
      DATA PRECIS/1.0D-19/
      DATA WW/0.7071067811865476D+00/
C
C---------------Executions of FSTRTN------------------------------------
C
      CTHETA=DCOS(THETA)
      IF (DABS(CTHETA).LE.PRECIS) CTHETA=0.0D0
      STHETA=DSIN(-THETA)
      IF (DABS(STHETA).LE.PRECIS) STHETA=0.0D0
C     STHETA=DSQRT(1.0D0-CTHETA*CTHETA)
      HSTHTA=WW*STHETA
      CTHTAP=WW*(1.0D0+CTHETA)
      CTHTAN=-WW*(1.0D0-CTHETA)
C
C     Initial setup for some coefficient matrix.
C
C     CALL BNLCFT(C, SQC, 4*NTERMS)
C
      D(0,0,0)=1.0D0
C
      DO 1000 IJ=1,NTERMS
c        if (ij.ge.3) stop
C
C     Compute the result for m=0 case, use formula (1).
C
         DO 200 IM=-IJ,-1
            D(IJ,0,IM)=-SQC(IJ-IM,2)*D(IJ-1,0,IM+1)
            IF (IM.GT.(1-IJ)) THEN
               D(IJ,0,IM)=D(IJ,0,IM)+SQC(IJ+IM,2)*D(IJ-1,0,IM-1)
            ENDIF
            D(IJ,0,IM)=D(IJ,0,IM)*HSTHTA
            IF (IM.GT.-IJ) THEN
               D(IJ,0,IM)=D(IJ,0,IM)+
     1           D(IJ-1,0,IM)*CTHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
            ENDIF
            D(IJ,0,IM)=D(IJ,0,IM)/IJ
200      CONTINUE
C
         D(IJ,0,0)=D(IJ-1,0,0)*CTHETA
         IF (IJ.GT.1) THEN
            D(IJ,0,0)=D(IJ,0,0)+HSTHTA*SQC(IJ,2)*(D(IJ-1,0,-1)+
     1                 D(IJ-1,0,1))/IJ
         ENDIF
C
         DO 400 IM=1,IJ
            D(IJ,0,IM)=-SQC(IJ+IM,2)*D(IJ-1,0,IM-1)
            IF (IM.LT.(IJ-1)) THEN
               D(IJ,0,IM)=D(IJ,0,IM)+SQC(IJ-IM,2)*D(IJ-1,0,IM+1)
            ENDIF 
            D(IJ,0,IM)=D(IJ,0,IM)*HSTHTA
            IF (IM.LT.IJ) THEN
               D(IJ,0,IM)=D(IJ,0,IM)+
     1           D(IJ-1,0,IM)*CTHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
            ENDIF
            D(IJ,0,IM)=D(IJ,0,IM)/IJ
400      CONTINUE
C
C     Compute the result for 0<m<=j case, use formula (2).
C
         DO 800 IMP=1,IJ
            DO 500 IM=-IJ,-1
               D(IJ,IMP,IM)=D(IJ-1,IMP-1,IM+1)*CTHTAN*SQC(IJ-IM,2)
               IF (IM.GT.(1-IJ)) THEN
                  D(IJ,IMP,IM)=D(IJ,IMP,IM)-
     1             D(IJ-1,IMP-1,IM-1)*CTHTAP*SQC(IJ+IM,2)
               ENDIF
               IF (IM.GT.-IJ) THEN
                  D(IJ,IMP,IM)=D(IJ,IMP,IM)+
     1             D(IJ-1,IMP-1,IM)*STHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
               ENDIF 
               D(IJ,IMP,IM)=D(IJ,IMP,IM)*WW/SQC(IJ+IMP,2)
500         CONTINUE
C
            D(IJ,IMP,0)=IJ*STHETA*D(IJ-1,IMP-1,0)
            IF (IJ.GT.1) THEN
               D(IJ,IMP,0)=D(IJ,IMP,0)-SQC(IJ,2)*(
     1            D(IJ-1,IMP-1,-1)*CTHTAP+D(IJ-1,IMP-1,1)*CTHTAN)
            ENDIF
            D(IJ,IMP,0)=D(IJ,IMP,0)*WW/SQC(IJ+IMP,2)
C
            DO 600 IM=1,IJ
               D(IJ,IMP,IM)=D(IJ-1,IMP-1,IM-1)*CTHTAP*SQC(IJ+IM,2)
               IF (IM.LT.(IJ-1)) THEN
                  D(IJ,IMP,IM)=D(IJ,IMP,IM)-
     1             D(IJ-1,IMP-1,IM+1)*CTHTAN*SQC(IJ-IM,2)
               ENDIF
               IF (IM.LT.IJ) THEN
                  D(IJ,IMP,IM)=D(IJ,IMP,IM)+
     1             D(IJ-1,IMP-1,IM)*STHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
               ENDIF
               D(IJ,IMP,IM)=D(IJ,IMP,IM)*WW/SQC(IJ+IMP,2)
600         CONTINUE
C
C     Use symmetricity, i.e. formula (3.80) in Biedenharn & Loucks
C     book, to compute the lower part of the matrix
C
C           DO IM=-IJ,IJ
C              D(IJ,-IMP,IM)=D(IJ,IMP,-IM)
C           ENDDO
C
800      CONTINUE
1000  CONTINUE
C
C----------End of FSTRTN------------------------------------------------
C
      RETURN
      END
C
