          PARAMETER (NMAX=1000)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          CHARACTER*1 DATO(120)
          DIMENSION XCOV(NMAX),COV(NMAX),TCOV(NMAX),SRV(NMAX),
     .              TSRV(NMAX)
          READ(1,1030) NCOV
          READ(1,1040) (DATO(I),I=1,120)
          NCOV=NCOV+1
          DO 5 I=1,NCOV
          READ(1,1000) XCOV(I),COV(I)
5         CONTINUE
          WRITE(*,1050)
          READ(*,*) IAUCR
          IF(IAUCR.EQ.1) THEN
          WRITE(*,1060)
          READ(*,*) S2X,S2Y
          ENDIF
          WRITE(*,1061)
          READ(*,*) IFINI,STPF
          STPF=STPF/2.
10        WRITE(*,1070)
          READ(*,*) NCASE,DB,DC
          CALL SELCAS(NCASE,DA,DB,DC,NCOV,XCOV,COV,TCOV)
          IF(IAUCR.EQ.1) THEN
          VAR=DSQRT(S2X*S2Y)
          ELSE
          VAR=COV(1)
          ENDIF
          WRITE(*,1240) NCASE,DA,DB,DC
          WRITE(*,1250)
          DO 20 I=1,NCOV
          WRITE(*,1200) XCOV(I),COV(I),TCOV(I)
20        CONTINUE
          DO 30 I=1,NCOV
          SRV(I)=COV(I)/VAR
          TSRV(I)=TCOV(I)/VAR
30        CONTINUE
          STP=XCOV(3)-XCOV(2)
          CALL GRACO(SRV,TSRV,STP,NCOV)
          
          IF(IFINI.EQ.1) THEN
          DO 40 I=1,NCOV
c          TSRV(I)=TCOV(I)*SCAL1(XCOV(I)/STPF)
          TSRV(I)=TCOV(I)*SCAL3(XCOV(I)/STPF)
          WRITE(*,1200) XCOV(I),COV(I),TSRV(I)
40        CONTINUE
          DO 50 I=1,NCOV
          SRV(I)=COV(I)/VAR
          TSRV(I)=TCOV(I)*SCAL1(XCOV(I)/STPF)/VAR
50        CONTINUE
          STP=XCOV(3)-XCOV(2)
          CALL GRACO(SRV,TSRV,STP,NCOV)
          ENDIF

          WRITE(*,1080)
          READ(*,*) ICONT
          IF(ICONT.EQ.1) GOTO 10
          STOP
1000      FORMAT(F10.3,F15.6)
1030      FORMAT(3I5)
1040      FORMAT(120A1)
1050      FORMAT(//,5X,'FUNZIONE EMPIRICA DA INTERPOLARE',
     ./5X,'AUTOCOVARIANZA = 0',/5X,'CROSSCOVARIANZA = 1')
1060      FORMAT(//5X,'PROCESSO CROSSCOVARIANZA (X,Y)',
     ./5X,'VARIANZA X, VARIANZA Y')
1061      FORMAT(//5X,'FUNZIONE FINITA (NO=0,SI=1) E'
     ./5X,'VALORE ASCISSA (F MODELLO = 0)')
1070      FORMAT(//5X,'TIPO E PARAMETRI FUNZIONE TEORICA',
     ./5X,'NCASE  DB  DC')
1080      FORMAT(//5X,'NUOVA INTERPOLAZIONE (NO=0,SI=1)')
1240      FORMAT(//5X,'NCASE =',I5,3X,'DA =',D15.8,3X,'DB =',D15.8,
     .3X,'DC =',D15.8)
1250      FORMAT(//7X,'TIME INTERVAL',11X,'EMPYRICAL',8X,'INTERPOLATED',
     ./9X,'OR DISTANCE',10X,'COVARIANCE',10X,'COVARIANCE',/32X,
     .'FUNCTION',12X,'FUNCTION',/)
1200      FORMAT(5X,F15.3,5X,F15.6,5X,F15.6)
          END
          SUBROUTINE SELCAS(NCASE,DA,DB,DC,NCOV,XCOV,COV,TCOV)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION XCOV(*),COV(*),TCOV(*)
          GOTO (5,15,25,35,45,55,65,75,85,95,105,115), NCASE
5         DA=COV(2)*DEXP(DB*XCOV(2))
          print *,'A*EXP(-B*X)'
          DO 10 I=1,NCOV
          TCOV(I)=DA*DEXP(-DB*XCOV(I))
10        CONTINUE
          RETURN
15        DA=COV(2)*DEXP(DB*XCOV(2)**2)
          print *,'A*EXP(-B*X**2)'
          DO 30 I=1,NCOV
          TCOV(I)=DA*DEXP(-DB*XCOV(I)**2)
30        CONTINUE
          RETURN
25        SRV=DEXP(-DB*XCOV(2))
          print *,'A*(1-C*X**2)*EXP(-B*X)'
          DA=COV(2)/((1.D0-DC*XCOV(2)**2)*SRV)
          DO 20 I=1,NCOV
          SRV=DEXP(-DB*XCOV(I))
          TCOV(I)=DA*(1.D0-DC*XCOV(I)**2)*SRV
20        CONTINUE
          RETURN
35        SRV=DEXP(-DB*XCOV(2)**2)
          print *,'A*(1-C*X**2)*EXP(-B*X**2)'
          DA=COV(2)/((1.D0-DC*XCOV(2)**2)*SRV)
          DO 40 I=1,NCOV
          SRV=DEXP(-DB*XCOV(I)**2)
          TCOV(I)=DA*(1.D0-DC*XCOV(I)**2)*SRV
40        CONTINUE
          RETURN
45        SRV=DEXP(-DB*XCOV(2))
          print *,'A*EXP(-B*X)*SIN(C*X)/(C*X)'
          DA=COV(2)/(SRV*DSIN(DC*XCOV(2))/(DC*XCOV(2)))
          TCOV(1)=DA
          DO 80 I=2,NCOV
          TCOV(I)=DA*DEXP(-DB*XCOV(I))*DSIN(DC*XCOV(I))/(DC*XCOV(I))
80        CONTINUE
          RETURN
55        SRV=DEXP(-DB*XCOV(2)**2)
          print *,'A*EXP(-B*X**2)*SIN(C*X)/(C*X)'
          DA=COV(2)/(SRV*DSIN(DC*XCOV(2))/(DC*XCOV(2)))
          TCOV(1)=DA
          DO 70 I=2,NCOV
          TCOV(I)=DA*DEXP(-DB*XCOV(I)**2)*DSIN(DC*XCOV(I))/(DC*XCOV(I))
70        CONTINUE
          RETURN
65        DA=COV(2)/(DEXP(-DB*XCOV(2))*DCOS(DC*XCOV(2)))
          print *,'A*EXP(-B*X)*COS(C*X)'
          DO 50 I=1,NCOV
          TCOV(I)=DA*DEXP(-DB*XCOV(I))*DCOS(DC*XCOV(I))
50        CONTINUE
          RETURN
75        DA=COV(2)/(DEXP(-DB*XCOV(2)**2)*DCOS(DC*XCOV(2)))
          print *,'A*EXP(-B*X**2)*COS(C*X)'
          DO 60 I=1,NCOV
          TCOV(I)=DA*DEXP(-DB*XCOV(I)**2)*DCOS(DC*XCOV(I))
60        CONTINUE
          RETURN
85        SRV=DEXP(-DB*XCOV(2))
          print *,'A*EXP(-B*X)*DFJ1(C*X)/(C*X)'
          DA=COV(2)/(SRV*DFJ1(DC*XCOV(2))/(DC*XCOV(2)))
          TCOV(1)=DA
          DO 120 I=2,NCOV
          TCOV(I)=DA*DEXP(-DB*XCOV(I))*DFJ1(DC*XCOV(I))/(DC*XCOV(I))
120       CONTINUE
          RETURN
95        SRV=DEXP(-DB*XCOV(2)**2)
          print *,'A*EXP(-B*X**2)*DFJ1(C*X)/(C*X)'
          DA=COV(2)/(SRV*DFJ1(DC*XCOV(2))/(DC*XCOV(2)))
          TCOV(1)=DA
          DO 110 I=2,NCOV
          TCOV(I)=DA*DEXP(-DB*XCOV(I)**2)*DFJ1(DC*XCOV(I))/(DC*XCOV(I))
110       CONTINUE
          RETURN
105       DA=COV(2)/(DEXP(-DB*XCOV(2))*DFJ0(DC*XCOV(2)))
          print *,'A*EXP(-B*X)*DFJ0(C*X)'
          DO 90 I=1,NCOV
          TCOV(I)=DA*DEXP(-DB*XCOV(I))*DFJ0(DC*XCOV(I))
90        CONTINUE
          RETURN
115       DA=COV(2)/(DEXP(-DB*XCOV(2)**2)*DFJ0(DC*XCOV(2)))
          print *,'A*EXP(-B*X**2)*DFJ0(C*X)'
          DO 100 I=1,NCOV
          TCOV(I)=DA*DEXP(-DB*XCOV(I)**2)*DFJ0(DC*XCOV(I))
100       CONTINUE
          RETURN
          END
C*****************************     GRACO
          SUBROUTINE GRACO(COV,TCOV,STP,NCOV)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 LB,PUN,AST,PIU,BL
          DIMENSION LB(101),A(11),COV(1),TCOV(1)
          DATA PUN,AST,PIU,BL/'.','*','+',' '/ 
          DO 10 I=1,11 
          A(I)=-1.D0+(I-1)*.2D0
10        CONTINUE 
          WRITE(6,100) (A(I),I=1,11)
          DO 20 I=1,101
          LB(I)=PUN
20        CONTINUE 
          DO 30 I=1,101,10 
          LB(I)='I'
30        CONTINUE 
          WRITE(6,110) (LB(K),K=1,101)  
          X0=STP/2.D0  
          X=0.D0
          DO 40 I=1,NCOV
          DO 35 K=1,101
          LB(K)=BL 
35        CONTINUE 
          LB(51)=PUN
          II=IDNINT(COV(I)*50.D0)+51
          JJ=IDNINT(TCOV(I)*50.D0)+51  
          LB(JJ)=PIU
          LB(II)=AST
          WRITE(6,120) X,(LB(J),J=1,101)
          X=X0+(I-1)*STP
40        CONTINUE 
          RETURN
100       FORMAT(9X,11(F4.1,6X))
110       FORMAT(11X,101A1)
120       FORMAT(1X,F10.3,101A1)
          END  
C*****************************     SCAL1
          FUNCTION SCAL1(Z)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          IF(Z.LT.2.D0) THEN
          R=1.D0
          Y1=1.0667D0  
          A=(16.D0/15.D0)*R**5 
          B=(4.D0/3.D0)*R**3
          C=(2.D0/3.D0)*R**2
          D=(1.D0/30.D0)
          Y=A-B*Z**2+C*Z**3-D*Z**5 
          SCAL1=Y/Y1
          ELSE 
          SCAL1=0.D0
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
C*****************************     SCAL3
          FUNCTION SCAL3(Z)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          IF(Z.LT.2.D0) THEN
          IZ=IDINT(Z/2.D0)+3
          SRV=SPLT(Z,IZ)
          SCAL3=(3.D0*SRV/2.D0)**3
          ELSE 
          SCAL3=0.D0
          ENDIF
          RETURN
          END  
C*****************************     SPLT
          FUNCTION SPLT(X,IX)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          GOTO (10,20,30,40),IX
10        SRV=(X+2.D0)**3  
          SPLT=SRV/6.D0
          RETURN
20        SRV=(X+2.D0)**3-4.D0*(X+1.D0)**3 
          SPLT=SRV/6.D0
          RETURN
30        SRV=(2.D0-X)**3-4.D0*(1.D0-X)**3 
          SPLT=SRV/6.D0
          RETURN
40        SRV=(2.D0-X)**3  
          SPLT=SRV/6.D0
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
  
