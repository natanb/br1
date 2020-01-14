        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        REAL*4 GASDEV,RAN1
        CHARACTER*80 NFILE,SDUMMY
        PARAMETER (NMAX=10000,NBMAX=2000)
        DIMENSION VAL(NMAX),X(NMAX),Y(NMAX)
C
C PARAMETRI DELLA FUNZIONE DI COVARIANZA
C
c NCASE select MAF
c DA DB DC coefficients of MAF
c SN noise in the origin of the MAF
c XCUT Abscissa where the MAF is definitively to zero 
c NX NY LK  grid parameters and matrix width band  
c IDUM negative number for starting pseudo-numbers generator
c
        READ(*,*) N
        READ(*,'(a80)') NFILE
        
        OPEN(UNIT=20,FILE=NFILE,STATUS='UNKNOWN')

        READ(20,*) NP,PMED,SQM
        READ(20,*) NCASE,DA,DB,DC,SN,TGL
        CLOSE(20)

        READ(*,'(a80)') NFILE
        OPEN(UNIT=20,FILE=NFILE,STATUS='UNKNOWN')
        
        DO I=1,N
        READ(20,*) NP,X(I),Y(I),PHA
        ENDDO
        CLOSE(20)

        READ(*,'(a80)') NFILE
        OPEN(UNIT=20,FILE=NFILE,STATUS='UNKNOWN')
        
C AZZERA MATRICE NORMALE
C
C COSTRUZIONE MATRICE DEL SEGNALE Css
C
        DO I=1,N
        DO J=1,I
        XCOV=SQRT((X(I)-X(J))**2+(Y(I)-Y(J))**2)
        VAL(J)=CMAF(NCASE,DA,DB,DC,XCOV)
        ENDDO
c        VAL(I)=VAL(I)+SN             !con o senza noise???
        write(20,'(f20.12)') (VAL(J),J=1,I)
        ENDDO

        END

          DOUBLE PRECISION FUNCTION CMAF(NCASE,DA,DB,DC,XCOV)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          GOTO (5,15,25,35,45,55,65,75,85,95,105,115), NCASE
5         CMAF=DA*DEXP(-DB*XCOV)
          RETURN
15        CMAF=DA*DEXP(-DB*XCOV**2)
          RETURN
25        CMAF=DA*(1.D0-DC*XCOV**2)*DEXP(-DB*XCOV)
          RETURN
35        CMAF=DA*(1.D0-DC*XCOV**2)*DEXP(-DB*XCOV**2)
          RETURN
45        CMAF=DA*DEXP(-DB*XCOV)*SIN(DC*XCOV)/(DC*XCOV)
          if(XCOV.EQ.0.D0) CMAF=DA
          RETURN
55        CMAF=DA*DEXP(-DB*XCOV**2)*DSIN(DC*XCOV)/(DC*XCOV)
          if(XCOV.EQ.0.D0) CMAF=DA
          RETURN
65        CMAF=DA*DEXP(-DB*XCOV)*DCOS(DC*XCOV)
          RETURN
75        CMAF=DA*DEXP(-DB*XCOV**2)*DCOS(DC*XCOV)
          RETURN
85        CMAF=DA*DEXP(-DB*XCOV)*DFJ1(DC*XCOV)/(DC*XCOV)
          if(XCOV.EQ.0.D0) CMAF=DA
          RETURN
95        CMAF=DA*DEXP(-DB*XCOV**2)*DFJ1(DC*XCOV)/(DC*XCOV)
          if(XCOV.EQ.0.D0) CMAF=DA
          RETURN
105       CMAF=DA*DEXP(-DB*XCOV)*DFJ0(DC*XCOV)
          RETURN
115       CMAF=DA*DEXP(-DB*XCOV**2)*DFJ0(DC*XCOV)
          RETURN
          END


      SUBROUTINE CLASCOOR(NC,NX,NY,NCX)
      IMPLICIT INTEGER(A-Z)
      NX=MOD(NC - 1,NCX) + 1
      NY=(NC - NX)/NCX + 1
      RETURN
      END


          FUNCTION DFN01(Z)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          REAL*4 SRW,SRV,ERFREC
          SRW=DABS(Z)/DSQRT(2.D0)
          SRV=DBLE(ERFREC(SRW))
          DFN01=.5D0+.5D0*SRV*Z/DABS(Z)
          RETURN
          END
C*************************    ERF
      FUNCTION ERFREC(X)
      IF(X.LT.0.)THEN
        ERFREC=-GAMMP(.5,X**2)
      ELSE
        ERFREC=GAMMP(.5,X**2)
      ENDIF
      RETURN
      END
C*************************    GAMMP
      FUNCTION GAMMP(A,X)
      IF(X.LT.0..OR.A.LE.0.)PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END
C*************************    GSER
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END

C*************************    GCF
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
C*************************    GAMMLN
      FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C*************************    BESSJ0
          FUNCTION DFJ0(X)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          REAL*4 SRW,SRV,BESSJ0
          SRW=X
          SRV=BESSJ0(SRW)
          DFJ0=DBLE(SRV)
          RETURN
          END
C
      FUNCTION BESSJ0(X)
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     *    S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     *    -.2073370639D-5,.2093887211D-6/, Q1,Q2,Q3,Q4,Q5/-.1562499995D-
     *1,
     *    .1430488765D-3,-.6911147651D-5,.7621095161D-6,-.934945152D-7/
      DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,651619640.7D
     *0,
     *    -11214424.18D0,77392.33017D0,-184.9052456D0/,
     *    S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,
     *    9494680.718D0,59272.64853D0,267.8532712D0,1.D0/
      IF(ABS(X).LT.8.)THEN
        Y=X**2
        BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=ABS(X)
        Z=8./AX
        Y=Z**2
        XX=AX-.785398164
        BESSJ0=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
      ENDIF
      RETURN
      END
C*************************    BESSJ1
          FUNCTION DFJ1(X)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          REAL*4 SRV,SRW,BESSJ1
          SRV=X
          SRW=BESSJ1(SRV)
          DFJ1=2.D0*DBLE(SRW)
          RETURN
          END
C
      FUNCTION BESSJ1(X)
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     *    S1,S2,S3,S4,S5,S6
      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,242396853.1D0
     *,
     *    -2972611.439D0,15704.48260D0,-30.16036606D0/,
     *    S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,
     *    18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5
     *,
     *    -.240337019D-6/, Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3
     *,
     *    .8449199096D-5,-.88228987D-6,.105787412D-6/
      IF(ABS(X).LT.8.)THEN
        Y=X**2
        BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=ABS(X)
        Z=8./AX
        Y=Z**2
        XX=AX-2.356194491
        BESSJ1=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
     *      *SIGN(1.,X)
      ENDIF
      RETURN
      END

C*****************************     SCAL2
          FUNCTION SCAL2(Z)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          IF(Z.LT.2.D0) THEN
          PIG=4.D0*DATAN(1.D0)
          R=1.D0
          Y1=5.D0*R**6*PIG
          A=5*R**6*PIG
          B=(15.D0/2.D0)*R**4*Z**2*PIG
          SQ=DSQRT(DABS(R**2-(Z/2.D0)**2))
          AS=DASIN(Z/(2.D0*R))
          C=5.D0*R**4*Z+(20.D0/3.D0)*R**2*Z**3-(5.D0/12.D0)*Z**5
          D=15.D0*R**4*Z**2-10.D0*R**6
          Y=A-B+C*SQ+D*AS
          SCAL2=Y/Y1
          ELSE
          SCAL2=0.D0
          ENDIF
          RETURN
          END


C*************************    GASDEV
      FUNCTION GASDEV(IDUM)
      DATA ISET/0/
      DATA GSET /0./
      IF (ISET.EQ.0) THEN
1     V1=2.*RAN1(IDUM)-1.
      V2=2.*RAN1(IDUM)-1.
      R=V1**2+V2**2
      IF(R.GE.1.)GO TO 1
      FAC=SQRT(-2.*LOG(R)/R)
      GSET=V1*FAC
      GASDEV=V2*FAC
      ISET=1
      ELSE
      GASDEV=GSET
      ISET=0
      ENDIF
      RETURN
      END


C*************************    RAN1
      FUNCTION RAN1(IDUM)
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/, IX1 /0/, IX2 /0/, IX3 /0/ 
      DATA R/97*0./
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1)PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END


