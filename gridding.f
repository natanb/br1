          PARAMETER (NPMAX=500000,NCMAX=9000,NCMAX1=NCMAX+1)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION IA(NPMAX),EZ(NCMAX),NNG(NCMAX1),LC(NCMAX1)
          DIMENSION Z(NPMAX),SZ(NCMAX)
          CHARACTER*80 FORMA
          READ(5,*) X0,X1,Y0,Y1,STP
          READ(5,*) NMPC
          READ(5,'(A80)') FORMA
          NCX=(X1-X0)/STP
          NCY=(Y1-Y0)/STP
          NCT=NCX*NCY
          NMAG=NPMAX/NCT
          WRITE(6,'(5X,''NCX NCY NCT NMAG'',4I5)') NCX,NCY,NCT,NMAG
          IF(NMAG.LT.NMPC) THEN
          WRITE(6,'(5X,''NUMERO PUNTI PER CLASSE'',/5X,
     .''ECCEDENTE IL DIMENSIONAMENTO PREVISTO'')')
          ENDIF
          NP=0
          NNB=0
10        READ(1,*,END=20) ZZ,XX,YY
          NNB=NNB+1
          IF(XX.LT.X0.OR.XX.GT.X1) GOTO 10
          IF(YY.LT.Y0.OR.YY.GT.Y1) GOTO 10
          NP=NP+1
          NX=(XX-X0)/STP
          NY=(YY-Y0)/STP
          NC=NY*NCX+NX+1
          EZ(NC)=EZ(NC)+1.D0
          NNG(NC)=NNG(NC)+1
          IA((NC-1)*NMPC+NNG(NC))=NP
          Z(NP)=ZZ
          GOTO 10
20        XC=0.D0
          YC=0.D0
          DO 30 II=1,NCT
          IF(NNG(II).GT.0) THEN
          WRITE(6,'(5X,''CLASSE N. '',I5,3X,I10)') II,NNG(II)
          ENDIF
30        CONTINUE
C         CALCOLA LE MEDIE PER CLASSE
C          DO 40 I=1,NCT
C          IF(NNG(I).GT.0) THEN
C          EZ(I)=EZ(I)/NNG(I)
C          ENDIF
C40        CONTINUE
          LC(1)=1
          JJ=0
          DO 48 I=1,NCT
          LC(I+1)=LC(I)+NNG(I)
          L=NNG(I)
          IF(L.GT.0) THEN
          II=(I-1)*NMPC
          DO 45 J=1,L
          JJ=JJ+1
          IA(JJ)=IA(II+J)
45        CONTINUE
          ENDIF
48        CONTINUE
          WRITE(6,'(///5X,''NUMERO PUNTI TOTALI'',I10//)') LC(NCT+1)-1
          WRITE(6,'(//5X,''JJ, NNB, NP'',3I10//)') JJ,NNB,NP
C         CALCOLA LE VARIANZE PER CLASSE
C          DO 58 I=1,NCT
C          I1=LC(I)
C          I2=LC(I+1)-1
C          IF(I1.LE.I2) THEN
C          DO 50 J=I1,I2
C          JS=IA(J)
C          SZ(I)=SZ(I)+(Z(JS)-EZ(I))**2
C50        CONTINUE
C          SZ(I)=DSQRT(SZ(I)/(NNG(I)-1))
C          ENDIF
C58        CONTINUE
          DO 70 I=1,NCY
          DO 70 J=1,NCX
          II=(I-1)*NCX+J
          XC=J*STP-STP/2.D0+X0
          YC=I*STP-STP/2.D0+Y0
          WRITE(2,1000) II,XC,YC,EZ(II),SZ(II),NNG(II)
70        CONTINUE
          STOP
1000      FORMAT(I5,2F15.3,2F15.3,I10)
          END
