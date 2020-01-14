          PARAMETER (NCMAX=9000,NCMAX1=NCMAX+1)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION EZ(NCMAX),NNG(NCMAX1)
          READ(5,*) X0,X1,Y0,Y1,STP
          NCX=(X1-X0)/STP
          NCY=(Y1-Y0)/STP
          NCT=NCX*NCY
          NP=0
10        READ(1,*,END=20) ZZ,XX,YY
          IF(XX.LT.X0.OR.XX.GT.X1) GOTO 10
          IF(YY.LT.Y0.OR.YY.GT.Y1) GOTO 10
          NP=NP+1
          NX=(XX-X0)/STP
          NY=(YY-Y0)/STP
          NC=NY*NCX+NX+1
          EZ(NC)=EZ(NC)+1.D0
          NNG(NC)=NNG(NC)+1
          GOTO 10
20        XC=0.D0
          YC=0.D0
          DO 70 I=1,NCY
          DO 70 J=1,NCX
          II=(I-1)*NCX+J
          XC=J*STP-STP/2.D0+X0
          YC=I*STP-STP/2.D0+Y0
c          IF(EZ(II).GT.0.) WRITE(2,1000) II,XC,YC,EZ(II)
          WRITE(2,1000) II,XC,YC,EZ(II)
70        CONTINUE
          STOP
1000      FORMAT(I5,2F15.3,F15.3)
          END
