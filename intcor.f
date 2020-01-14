          PARAMETER (MAX=500)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO 
          CHARACTER*4 TITOLO
          DIMENSION COV(MAX),XCOV(MAX),TCOV(MAX),X0(4),IX0(3),IXMM(2), 
     .SIGMA2(6),SGNL(6),XCOVT(6),IND(12), 
     .DATO(120),TITOLO(20) 
          DATA IND/0,0,0,0,1,2,0,0,1,4,2,4/
C
          CALL LINGUA(ILI)
          IF(ILI.EQ.1) THEN
           CALL PAINTC(TITOLO,ICASO,ITIPO,IFREE,ISCAL,IFILE,XCOVT)
           WRITE(6,1020) TITOLO 
          ELSE
           CALL PIINTC(TITOLO,ICASO,ITIPO,IFREE,ISCAL,IFILE,XCOVT)
           WRITE(6,1020) TITOLO 
          ENDIF
          CALL LEGEND(ICASO,ILI)
          NTER=0
          IF(ITIPO.EQ.1) NTER=1
          IF(ITIPO.EQ.2) NTER=3
          IF(ITIPO.EQ.3) NTER=6
          DO 200 ITER=1,NTER
          NCASE=0  
          DA=0.D0  
          DB=0.D0  
          DC=0.D0  
          SNCOV=0.D0
          SNCRO=0.D0
          DD=0.D0  
          IF(ILI.EQ.1) THEN
           IF(ITER.EQ.1) WRITE(6,2600)  
           IF(ITER.EQ.2) WRITE(6,2700)  
           IF(ITER.EQ.3) WRITE(6,2800)  
           IF(ITER.EQ.4) WRITE(6,2900)  
           IF(ITER.EQ.5) WRITE(6,3000)  
           IF(ITER.EQ.6) WRITE(6,3100)
          ELSE  
           IF(ITER.EQ.1) WRITE(6,2601)  
           IF(ITER.EQ.2) WRITE(6,2701)  
           IF(ITER.EQ.3) WRITE(6,2801)  
           IF(ITER.EQ.4) WRITE(6,2901)  
           IF(ITER.EQ.5) WRITE(6,3001)  
           IF(ITER.EQ.6) WRITE(6,3101)
          ENDIF
          READ(1,1030) NCOV
          READ(1,1040) (DATO(I),I=1,120)
          NCOV=NCOV+1  
          IF(NCOV.GT.MAX) THEN 
           IF(ILI.EQ.1) THEN
            WRITE(6,1005)
           ELSE
            WRITE(6,1006)
           ENDIF
           STOP 
          ENDIF
          DO 5 I=1,NCOV
          READ(1,1000) XCOV(I),COV(I)  
5         CONTINUE 
C
          IF(ISCAL.EQ.0) XCOVT(ITER)=1.D30
           SIGMA2(ITER)=COV(1)  
           PROD=COV(1)*COV(2)
           IF(((ITER.EQ.1.OR.ITER.EQ.2.OR.ITER.EQ.4).AND.
     .        (COV(1).GT.0.D0.AND.COV(2).GT.0.D0)).OR.
     .        ((ITER.EQ.3.OR.ITER.EQ.5.OR.ITER.EQ.6).AND.
     .        PROD.GT.0.D0)) THEN
           IDELT2=IDIF2(COV)
           CALL ZERO(COV,XCOV,IX0,X0,N0,NCOV)
           IF(N0.NE.0) THEN 
            CALL MIMA(COV,IXMM,IX0,NCOV)
            RAP0=(X0(2)-X0(1))/X0(1) 
           ENDIF 
           IF(IDELT2.EQ.1.AND.N0.EQ.0) THEN 
                    CALL CALCO1(COV,XCOV,DA,DB)  
                    NCASE=1
                    IF(ILI.EQ.1) THEN  
                    WRITE(6,1300)
                    ELSE
                    WRITE(6,1301)
                    ENDIF
                    WRITE(6,1090) DA,DB  
                    DO 10 I=1,NCOV
                    TCOV(I)=DA*DEXP(-DB*XCOV(I)) 
10                  CONTINUE 
          ENDIF
          IF(IDELT2.EQ.2.AND.N0.EQ.0) THEN 
                    CALL CALCO2(COV,XCOV,DA,DB)  
                    NCASE=2 
                    IF(ILI.EQ.1) THEN 
                    WRITE(6,1400)
                    ELSE
                    WRITE(6,1401)
                    ENDIF
                    WRITE(6,1090) DA,DB  
                    DO 20 I=1,NCOV
                    TCOV(I)=DA*DEXP(-DB*XCOV(I)**2)  
20                  CONTINUE 
          ENDIF
          IF(IDELT2.EQ.1.AND.N0.EQ.1) THEN 
                    CALL CALCO3(COV,XCOV,X0,IXMM,DA,DB,DD,ICASO) 
                    NCASE=3 
                    IF(ILI.EQ.1) THEN 
                    IF(ICASO.EQ.1) WRITE(6,1500)
                    IF(ICASO.EQ.2) WRITE(6,1502)
                    IF(ICASO.EQ.3) WRITE(6,1505)
                    ELSE
                    IF(ICASO.EQ.1) WRITE(6,1501)
                    IF(ICASO.EQ.2) WRITE(6,1503)
                    IF(ICASO.EQ.3) WRITE(6,1506)
                    ENDIF
                    WRITE(6,1095) DA,DB,DD
                    DO 30 I=1,NCOV
                    SRV=DEXP(-DB*XCOV(I))
                    TCOV(I)=DA*(1.D0-DD*XCOV(I)**2)*SRV  
30                  CONTINUE 
          ENDIF
          IF(IDELT2.EQ.2.AND.N0.EQ.1) THEN 
                    CALL CALCO4(COV,XCOV,X0,IXMM,DA,DB,DD,ICASO) 
                    NCASE=4
                    IF(ILI.EQ.1) THEN  
                    IF(ICASO.EQ.1) WRITE(6,1600)
                    IF(ICASO.EQ.2) WRITE(6,1602)
                    IF(ICASO.EQ.3) WRITE(6,1605)
                    ELSE
                    IF(ICASO.EQ.1) WRITE(6,1601)
                    IF(ICASO.EQ.2) WRITE(6,1603)
                    IF(ICASO.EQ.3) WRITE(6,1606)
                    ENDIF
                    WRITE(6,1095) DA,DB,DD
                    DO 40 I=1,NCOV
                    SRV=DEXP(-DB*XCOV(I)**2) 
                    TCOV(I)=DA*(1.D0-DD*XCOV(I)**2)*SRV  
40                  CONTINUE 
          ENDIF
          IF(N0.GT.1.AND.(ICASO.EQ.1.OR.ICASO.EQ.3)) THEN  
          IF((IDELT2.EQ.1.AND.RAP0.LT.1.5D0.AND.ICASO.EQ.1).OR.
     .    (IDELT2.EQ.1.AND.IFREE.EQ.1.AND.ICASO.EQ.3)) THEN
                    CALL CALCO5(COV,XCOV,X0,IXMM,IX0,DA,DB,DC)
                    NCASE=5
                    IF(ILI.EQ.1) THEN
                    WRITE(6,1700)
                    ELSE
                    WRITE(6,1701)
                    ENDIF
                    WRITE(6,1100) DA,DB,DC
                    TCOV(1)=DA
                    DO 50 I=2,NCOV
                    TCOV(I)=DA*DEXP(-DB*XCOV(I))*DSIN(DC*XCOV(I))
     .              /(DC*XCOV(I))
50                  CONTINUE 
          ENDIF
          IF((IDELT2.EQ.2.AND.RAP0.LT.1.5D0.AND.ICASO.EQ.1).OR.
     .    (IDELT2.EQ.2.AND.IFREE.EQ.1.AND.ICASO.EQ.3)) THEN
                    CALL CALCO6(COV,XCOV,X0,IXMM,IX0,DA,DB,DC)
                    NCASE=6
                    IF(ILI.EQ.1) THEN  
                    WRITE(6,1800)
                    ELSE
                    WRITE(6,1801)
                    ENDIF
                    WRITE(6,1100) DA,DB,DC
                    TCOV(1)=DA
                    DO 60 I=2,NCOV
                    TCOV(I)=DA*DEXP(-DB*XCOV(I)**2)*DSIN(DC*XCOV(I))
     .              /(DC*XCOV(I))
60                  CONTINUE 
          ENDIF
          ENDIF
          IF(N0.GT.1.AND.ICASO.EQ.1) THEN
          IF(IDELT2.EQ.1.AND.RAP0.GE.1.5D0) THEN
                    CALL CALCO7(COV,XCOV,X0,IXMM,IX0,DA,DB,DC)
                    NCASE=7
                    IF(ILI.EQ.1) THEN  
                    WRITE(6,1900)
                    ELSE
                    WRITE(6,1901)
                    ENDIF
                    WRITE(6,1100) DA,DB,DC
                    DO 70 I=1,NCOV
                    TCOV(I)=DA*DEXP(-DB*XCOV(I))*DCOS(DC*XCOV(I))
70                  CONTINUE 
          ENDIF
          IF(IDELT2.EQ.2.AND.RAP0.GE.1.5D0) THEN
                    CALL CALCO8(COV,XCOV,X0,IXMM,IX0,DA,DB,DC)
                    NCASE=8
                    IF(ILI.EQ.1) THEN  
                    WRITE(6,2000)
                    ELSE
                    WRITE(6,2001)
                    ENDIF
                    WRITE(6,1100) DA,DB,DC
                    DO 80 I=1,NCOV
                    TCOV(I)=DA*DEXP(-DB*XCOV(I)**2)*DCOS(DC*XCOV(I)) 
80                  CONTINUE 
          ENDIF
          ENDIF
          IF(N0.GT.1.AND.(ICASO.EQ.2.OR.ICASO.EQ.3)) THEN  
          IF((IDELT2.EQ.1.AND.RAP0.LT.1.063D0.AND.ICASO.EQ.2).OR.
     .    (IDELT2.EQ.1.AND.IFREE.EQ.2.AND.ICASO.EQ.3)) THEN 
                    CALL CALCO9(COV,XCOV,X0,IXMM,IX0,DA,DB,DC)
                    NCASE=9 
                    IF(ILI.EQ.1) THEN
                    WRITE(6,2100)
                    ELSE
                    WRITE(6,2101) 
                    ENDIF
                    WRITE(6,1100) DA,DB,DC
                    TCOV(1)=DA
                    DO 90 I=2,NCOV  
                    TCOV(I)=DA*DEXP(-DB*XCOV(I))*DFJ1(DC*XCOV(I))
     .              /(DC*XCOV(I))
90                  CONTINUE 
          ENDIF
          IF((IDELT2.EQ.2.AND.RAP0.LT.1.063D0.AND.ICASO.EQ.2).OR.
     .    (IDELT2.EQ.2.AND.IFREE.EQ.2.AND.ICASO.EQ.3)) THEN 
                    CALL CALC10(COV,XCOV,X0,IXMM,IX0,DA,DB,DC)
                    NCASE=10
                    IF(ILI.EQ.1) THEN
                    WRITE(6,2200)
                    ELSE
                    WRITE(6,2201)
                    ENDIF
                    WRITE(6,1100) DA,DB,DC
                    TCOV(1)=DA
                    DO 100 I=2,NCOV  
                    TCOV(I)=DA*DEXP(-DB*XCOV(I)**2)*DFJ1(DC*XCOV(I))
     .              /(DC*XCOV(I))
100                 CONTINUE 
          ENDIF
          ENDIF
          IF(N0.GT.1.AND.ICASO.EQ.2) THEN
          IF(IDELT2.EQ.1.AND.RAP0.GE.1.063D0) THEN 
                    CALL CALC11(COV,XCOV,X0,IXMM,IX0,DA,DB,DC)
                    NCASE=11
                    IF(ILI.EQ.1) THEN
                    WRITE(6,2300)
                    ELSE
                    WRITE(6,2301)
                    ENDIF
                    WRITE(6,1100) DA,DB,DC
                    DO 110 I=1,NCOV
                    TCOV(I)=DA*DEXP(-DB*XCOV(I))*DFJ0(DC*XCOV(I))
110                 CONTINUE 
          ENDIF
          IF(IDELT2.EQ.2.AND.RAP0.GE.1.063D0) THEN 
                    CALL CALC12(COV,XCOV,X0,IXMM,IX0,DA,DB,DC)
                    NCASE=12 
                    IF(ILI.EQ.1) THEN
                    WRITE(6,2400)
                    ELSE
                    WRITE(6,2401)
                    ENDIF
                    WRITE(6,1100) DA,DB,DC
                    DO 120 I=1,NCOV  
                    TCOV(I)=DA*DEXP(-DB*XCOV(I)**2)*DFJ0(DC*XCOV(I)) 
120                 CONTINUE 
          ENDIF
          ENDIF
C
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1040) (DATO(I),I=1,120)  
          WRITE(2,1030) ICASO,ITIPO,NCASE  
          IF(NCASE.EQ.1.OR.NCASE.EQ.2) WRITE(2,3300) DA,DB 
          IF(NCASE.EQ.3.OR.NCASE.EQ.4) WRITE(2,3400) DA,DB,DD  
          IF(NCASE.GE.5.AND.NCASE.LE.12) WRITE(2,3400) DA,DB,DC
          DO 125 I=1,NCOV  
          WRITE(2,1060) XCOV(I),COV(I),TCOV(I) 
125       CONTINUE 
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1250)
          ELSE
          WRITE(6,1251)
          ENDIF
          DO 130 I=1,NCOV  
          WRITE(6,1200) XCOV(I),COV(I),TCOV(I) 
130       CONTINUE 
          SGNL(ITER)=DA
          VAR=0.D0 
          IF(ITER.EQ.1.OR.ITER.EQ.2.OR.ITER.EQ.4) THEN 
          SQM=DSQRT(SIGMA2(ITER))  
          SIGNAL=DSQRT(SGNL(ITER)) 
          SNCOV=SIGMA2(ITER)-SGNL(ITER)
          SNOISE=DSQRT(SIGMA2(ITER)-SGNL(ITER))
          IF(ILI.EQ.1) THEN
          WRITE(6,2500) SQM,SIGNAL,SNOISE
          ELSE  
          WRITE(6,2501) SQM,SIGNAL,SNOISE
          ENDIF  
          VAR=SIGMA2(ITER) 
          ENDIF
          IF(ITER.EQ.3.OR.ITER.EQ.5.OR.ITER.EQ.6) THEN 
	  LOOP=0
          I1=IND(2*ITER-1) 
          I2=IND(2*ITER)
          SNCRO=SIGMA2(ITER)-SGNL(ITER)
132	  LOOP=LOOP+1
	  IF(LOOP.LE.100) THEN
          IF(SGNL(ITER)**2.GT.SGNL(I1)*SGNL(I2)) THEN  
          SGNL(ITER)=DSQRT(.98D0*SGNL(I1)*SGNL(I2))  
          SNCRO=SIGMA2(ITER)-SGNL(ITER)
          ENDIF
          IF(SNCRO**2.GT.(SIGMA2(I1)-SGNL(I1))*(SIGMA2(I2)-SGNL(I2)))  
     .THEN 
          SNCRO=DSQRT(.98D0*(SIGMA2(I1)-SGNL(I1))*(SIGMA2(I2)-SGNL(I2))  
     .)
          SGNL(ITER)=SIGMA2(ITER)-SNCRO
          GOTO 132 
	  ENDIF
          ENDIF
          ERRE=SIGMA2(ITER)/DSQRT(SIGMA2(I1)*SIGMA2(I2))
          ERRES=SGNL(ITER)/DSQRT(SGNL(I1)*SGNL(I2))
          ERREN=(SIGMA2(ITER)-SGNL(ITER))/DSQRT((SIGMA2(I1)-SGNL(I1))* 
     .(SIGMA2(I2)-SGNL(I2)))
          WRITE(6,2550) ERRE,ERRES,ERREN
          VAR=DSQRT(SIGMA2(I1)*SIGMA2(I2)) 
          ENDIF
          DO 135 I=1,NCOV  
          COV(I)=COV(I)/VAR
          TCOV(I)=TCOV(I)/VAR  
135       CONTINUE 
          STP=XCOV(3)-XCOV(2) 
          IF(COV(1) .GE. TCOV(1)) CALL GRACO(COV,TCOV,STP,NCOV)
C
          IF(ISCAL.EQ.1) THEN  
          XCVT=XCOVT(ITER) 
          IF(ILI.EQ.1) THEN
          WRITE(6,1280)
          ELSE
          WRITE(6,1281)
          ENDIF
          IF(ICASO.EQ.1) THEN  
                    DO 140 I=1,NCOV  
                    IF(XCOV(I).LE.2.D0*XCVT) THEN
                    ZCOV=XCOV(I)/XCVT
                    TCOV(I)=TCOV(I)*SCAL1(ZCOV)  
                    ELSE 
                    TCOV(I)=0.D0 
                    ENDIF
                    WRITE(6,1200) XCOV(I),COV(I)*VAR,TCOV(I)*VAR 
140                 CONTINUE 
          ENDIF
          IF(ICASO.EQ.2) THEN  
                   DO 150 I=1,NCOV  
                   IF(XCOV(I).LE.2.D0*XCVT) THEN
                   ZCOV=XCOV(I)/XCVT
                   TCOV(I)=TCOV(I)*SCAL2(ZCOV)  
                   ELSE 
                   TCOV(I)=0.D0 
                   ENDIF
                   WRITE(6,1200) XCOV(I),COV(I)*VAR,TCOV(I)*VAR 
150                CONTINUE 
          ENDIF
          IF(ICASO.EQ.3) THEN  
                   DO 160 I=1,NCOV  
                   IF(XCOV(I).LE.2.D0*XCVT) THEN
                   ZCOV=XCOV(I)/XCVT
                   TCOV(I)=TCOV(I)*SCAL3(ZCOV)  
                   ELSE 
                   TCOV(I)=0.D0 
                   ENDIF
                   WRITE(6,1200) XCOV(I),COV(I)*VAR,TCOV(I)*VAR 
160                CONTINUE
          ENDIF
           
          WRITE(6,3200)
          CALL GRACO(COV,TCOV,STP,NCOV)
          ENDIF
          ELSE
C
          IF(ILI.EQ.1) THEN
          WRITE(6,1050)
          ELSE
          WRITE(6,1051)
          ENDIF
          SNCRO=COV(1)
          NCASE=1  
          ENDIF
C
          IF(IFILE.EQ.1.AND.(ITER.EQ.1.OR.ITER.EQ.2.OR.ITER.EQ.4)) THEN
          WRITE(3,1040) (DATO(I),I=1,120)  
          IF(NCASE.EQ.1.OR.NCASE.EQ.2) 
     .WRITE(3,3500) NCASE,SGNL(ITER),DB,DC,SNCOV,XCOVT(ITER)
          IF(NCASE.EQ.3.OR.NCASE.EQ.4) 
     .WRITE(3,3500) NCASE,SGNL(ITER),DB,DD,SNCOV,XCOVT(ITER)
          IF(NCASE.GE.5.AND.NCASE.LE.12)
     .WRITE(3,3500) NCASE,SGNL(ITER),DB,DC,SNCOV,XCOVT(ITER)
          ENDIF
          IF(IFILE.EQ.1.AND.(ITER.EQ.3.OR.ITER.EQ.5.OR.ITER.EQ.6)) THEN
          WRITE(3,1040) (DATO(I),I=1,120)  
          IF(NCASE.EQ.1.OR.NCASE.EQ.2) 
     .WRITE(3,3500) NCASE,SGNL(ITER),DB,DC,SNCRO,XCOVT(ITER)
          IF(NCASE.EQ.3.OR.NCASE.EQ.4) 
     .WRITE(3,3500) NCASE,SGNL(ITER),DB,DD,SNCRO,XCOVT(ITER)
          IF(NCASE.GE.5.AND.NCASE.LE.12)
     .WRITE(3,3500) NCASE,SGNL(ITER),DB,DC,SNCRO,XCOVT(ITER)
          ENDIF
200       CONTINUE 
          STOP 
1000      FORMAT(F10.3,F15.6)  
1005      FORMAT(5X,'NUMERO DEI DATI ',/5X,
     .'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'///)  
1006      FORMAT(5X,'NUMBER OF DATA ',/5X,
     .'EXCEEDING THE PREVIOUS DIMENSION'///)  
1010      FORMAT(20A4) 
1020      FORMAT(/////5X,20A4,//)  
1030      FORMAT(3I5)  
1040      FORMAT(120A1)
1050      FORMAT(/////5X,'GAMMA(0) E/O GAMMA(1) MINORE DI ZERO: '  
     .'INTERPOLAZIONE IMPOSSIBILE',/////)
1051      FORMAT(/////5X,'GAMMA(0) AND/OR GAMMA(1) LESS THAN ZERO: ' 
     .'IMPOSSIBLE INTERPOLATION',/////)
1060      FORMAT(F10.3,2F15.6) 
1090      FORMAT(//5X,'A =',D15.8,/5X,'B =',D15.8,/5X) 
1095      FORMAT(//5X,'A =',D15.8,/5X,'B =',D15.8,/5X,'D =',D15.8,/)
1100      FORMAT(//5X,'A =',D15.8,/5X,'B =',D15.8,/5X,'C =',D15.8,/)
1200      FORMAT(5X,F15.3,5X,F15.6,5X,F15.6)
1250      FORMAT(//10X,'INTERVALLO',10X,'FUNZIONE',12X,'FUNZIONE', 
     ./10X,'O DISTANZA',10X,'COVARIANZA',10X,'COVARIANZA',/30X,
     .'EMPIRICA',12X,'TEORICA',/)  
1251      FORMAT(//9X,'INTERVAL',13X,'EMPYRICAL',10X,'THEORETICAL', 
     ./9X,'OR DISTANCE',10X,'COVARIANCE',9X,'COVARIANCE',/30X,
     .'FUNCTION',11X,'FUNCTION',/)  
1280      FORMAT(//10X,'INTERVALLO',10X,'FUNZIONE',12X,'FUNZIONE', 
     ./10X,'O DISTANZA',10X,'COVARIANZA',10X,'COVARIANZA',/30X,
     .'EMPIRICA',12X,'TEORICA',/50X,'FINITA',/)  
1281      FORMAT(//9X,'INTERVAL',13X,'EMPYRICAL',10X,'FINITE', 
     ./9X,'OR DISTANCE',10X,'COVARIANCE',9X,'THEORETICAL', 
     ./30X,'FUNCTION',11X,'COVARIANCE',/49X,'FUNCTION',/)
1300      FORMAT(5X,'CASO ESPONENZIALE',// 
     .5X,'A*EXP(-B*ABS(X))')
1301      FORMAT(5X,'EXPONENTIAL CASE',// 
     .5X,'A*EXP(-B*ABS(X))')
1400      FORMAT(5X,'CASO NORMALE',//  
     .5X,'A*EXP(-B*X**2)') 
1401      FORMAT(5X,'NORMAL CASE',//  
     .5X,'A*EXP(-B*X**2)') 
1500      FORMAT(5X,'CASO ESPONENZIALE PARABOLA'
     .,//5X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(2*D)') 
1501     FORMAT(5X,'EXPONENTIAL PARABOLA CASE'
     .,//5X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(2*D)') 
1502      FORMAT(5X,'CASO ESPONENZIALE PARABOLA'
     .,//5X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(6*D)') 
1503      FORMAT(5X,'EXPONENTIAL PARABOLA CASE'
     .,//5X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(6*D)') 
1505     FORMAT(5X,'CASO ESPONENZIALE PARABOLA'
     .,//5X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(24*D)') 
1506     FORMAT(5X,'EXPONENTIAL PARABOLA CASE'
     .,//5X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(24*D)') 
1600      FORMAT(5X,'CASO NORMALE PARABOLA',// 
     .5X,'A*(1-D*X**2)*EXP(-B*X**2)   B.GT.D/2')  
1601      FORMAT(5X,'NORMAL PARABOLA CASE',// 
     .5X,'A*(1-D*X**2)*EXP(-B*X**2)   B.GT.D/2')  
1602      FORMAT(5X,'CASO NORMALE PARABOLA',// 
     .5X,'A*(1-D*X**2)*EXP(-B*X**2)   B.GT.D')  
1603      FORMAT(5X,'NORMAL PARABOLA CASE',// 
     .5X,'A*(1-D*X**2)*EXP(-B*X**2)   B.GT.D')  
1605     FORMAT(5X,'CASO NORMALE PARABOLA',// 
     .5X,'A*(1-D*X**2)*EXP(-B*X**2)   B.GT.2*D')  
1606     FORMAT(5X,'NORMAL PARABOLA CASE',// 
     .5X,'A*(1-D*X**2)*EXP(-B*X**2)   B.GT.2*D')  
1700      FORMAT(5X,'CASO ESPONENZIALE SENO SU ASCISSA',// 
     .5X,'A*EXP(-B*ABS(X))*SIN(C*ABS(X))/(C*ABS(X))') 
1701      FORMAT(5X,'EXPONENTIAL SINE OVER X COORDINATE CASE',// 
     .5X,'A*EXP(-B*ABS(X))*SIN(C*ABS(X))/(C*ABS(X))')  
1800      FORMAT(5X,'CASO NORMALE SENO SU ASCISSA',//5X,
     .'A*EXP(-B*X**2)*SIN(C*ABS(X))/(C*ABS(X))')
1801      FORMAT(5X,'NORMAL SINE OVER X COORDINATE CASE',//5X,
     .'A*EXP(-B*X**2)*SIN(C*ABS(X))/(C*ABS(X))')
1900      FORMAT(5X,'CASO ESPONENZIALE COSENO',//  
     .5X,'A*EXP(-B*ABS(X))*COS(C*ABS(X))') 
1901      FORMAT(5X,'EXPONENTIAL COSINE CASE',//  
     .5X,'A*EXP(-B*ABS(X))*COS(C*ABS(X))') 
2000      FORMAT(5X,'CASO NORMALE COSENO',//5X,
     .'A*EXP(-B*X**2)*COS(C*ABS(X))')  
2001      FORMAT(5X,'NORMAL COSINE CASE',//5X,
     .'A*EXP(-B*X**2)*COS(C*ABS(X))')  
2100      FORMAT(5X,'CASO ESPONENZIALE J1 SU ASCISSA',//
     .5X,'2*A*EXP(-B*ABS(X))*J1(C*ABS(X))/(C*ABS(X))') 
2101      FORMAT(5X,'EXPONENTIAL J1 OVER X COORDINATE CASE',//
     .5X,'2*A*EXP(-B*ABS(X))*J1(C*ABS(X))/(C*ABS(X))') 
2200      FORMAT(5X,'CASO NORMALE J1 SU ASCISSA',//5X, 
     .'2*A*EXP(-B*X**2)*J1(C*ABS(X))/(C*ABS(X))')  
2201      FORMAT(5X,'NORMAL J1 OVER X COORDINATE',//5X, 
     .'2*A*EXP(-B*X**2)*J1(C*ABS(X))/(C*ABS(X))')  
2300      FORMAT(5X,'CASO ESPONENZIALE J0',
     .//5X,'A*EXP(-B*ABS(X))*J0(C*ABS(X))')
2301      FORMAT(5X,'EXPONENTIAL J0 CASE',
     .//5X,'A*EXP(-B*ABS(X))*J0(C*ABS(X))')
2400      FORMAT(5X,'CASO NORMALE J0',//5X,
     .'A*EXP(-B*X**2)*J0(C*ABS(X))')
2401      FORMAT(5X,'NORMAL J0 CASE',//5X,
     .'A*EXP(-B*X**2)*J0(C*ABS(X))')
2500      FORMAT(//7X,'SQM',10X,'SQM-S',10X,'SQM-N',//3(F10.3,5X),////)
2501      FORMAT(//8X,'SD',11X,'SD-S',11X,'SD-N',//3(F10.3,5X),////)
2550      FORMAT(//8X,'RO',11X,'RO-S',11X,'RO-N',//3(F10.3,5X),////)
2600      FORMAT(///5X,'FUNZIONE AUTOCOVARIANZA X',///)
2601      FORMAT(///5X,'X AUTOCOVARIANCE FUNCTION',///)
2700      FORMAT(///5X,'FUNZIONE AUTOCOVARIANZA Y',///)
2701      FORMAT(///5X,'Y AUTOCOVARIANCE FUNCTION',///)
2800      FORMAT(///5X,'FUNZIONE CROSSCOVARIANZA X Y',///) 
2801      FORMAT(///5X,'X Y CROSSCOVARIANCE FUNCTION',///) 
2900      FORMAT(///5X,'FUNZIONE AUTOCOVARIANZA Z',///)
2901      FORMAT(///5X,'Z AUTOCOVARIANCE FUNCTION',///)
3000      FORMAT(///5X,'FUNZIONE CROSSCOVARIANZA X Z',///) 
3001      FORMAT(///5X,'X Z CROSSCOVARIANCE FUNCTION',///) 
3100      FORMAT(///5X,'FUNZIONE CROSSCOVARIANZA Y Z',///) 
3101      FORMAT(///5X,'Y Z CROSSCOVARIANCE FUNCTION',///) 
3200      FORMAT(/////)
3300      FORMAT(2D15.8)
3400      FORMAT(3D15.8)
3500      FORMAT(I5,5D15.8)
3600      FORMAT(2I5,F10.3)
          END  
