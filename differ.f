          PARAMETER (LMAX=100,NMAX=LMAX**2)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*4 TITOLO,FORMA 
          CHARACTER*10 VOLUME,BLANK
          DIMENSION Z(NMAX),TITOLO(20),FORMA(20)
          DATA BLANK/'          '/ 
          PIG=4.D0*DATAN(1.D0)
          CONV=200.D0/PIG
C          
          CALL LINGUA(ILI)
          IF(ILI.EQ.1) THEN
          WRITE(*,1000)
          READ(*,1010) TITOLO  
1         WRITE(*,1050)
          READ(*,*,ERR=1) IFILE  
2         WRITE(*,1100)
          READ(*,*,ERR=2) ALATO  
3         WRITE(*,1110)
          READ(*,*,ERR=3) SCX,SCY,SCZ
          WRITE(*,1140)
          READ(*,1010) FORMA
          WRITE(*,1150)
          ELSE
          WRITE(*,1001)
          READ(*,1010) TITOLO  
5         WRITE(*,1051)
          READ(*,*,ERR=5) IFILE  
7         WRITE(*,1101)
          READ(*,*,ERR=7) ALATO  
8         WRITE(*,1111)
          READ(*,*,ERR=8) SCX,SCY,SCZ
          WRITE(*,1141)
          READ(*,1010) FORMA
          WRITE(*,1151)
          ENDIF
C
          XMIN=1.D30  
          YMIN=1.D30  
          XMAX=-1.D30 
          YMAX=-1.D30 
          I=0  
10        READ(1,FORMA) NPT,XX,YY,ZZ
          IF(NPT.EQ.0) GO TO 20
          I=I+1     
          XMIN=DMIN1(XMIN,XX*SCX)
          XMAX=DMAX1(XMAX,XX*SCX)
          YMIN=DMIN1(YMIN,YY*SCY)          
          YMAX=DMAX1(YMAX,YY*SCY)
          Z(I)=SCZ*ZZ  
          GO TO 10 
20        N=I           
          K=IDNINT((XMAX-XMIN)/ALATO)+1
          L=IDNINT((YMAX-YMIN)/ALATO)+1
C
          WRITE(6,1020) TITOLO 
          IF(ILI.EQ.1) THEN
          WRITE(6,1200)
          WRITE(6,1250)
          ELSE
          WRITE(6,1201)
          WRITE(6,1251)
          ENDIF
          K1=K-1
          L1=L-1
          NPIU=0
          NMENO=0  
          VPIU=0.D0
          VMENO=0.D0
          XP=0.D0  
          YP=0.D0  
          XM=0.D0  
          YM=0.D0  
          DO 40 J=1,L1 
          IJ=J 
          IJ2=IJ+1 
          I2J=K+J  
          I2J2=I2J+1
          Z1=Z(IJ) 
          Z2=Z(I2J)
          Z3=Z(IJ2)
          Z4=Z(I2J2)
          XG=ALATO/2.D0  
          YG=((J-1)+0.5D0)*ALATO  
          IF((Z1.LT.0.D0.AND.Z2.LT.0.D0.AND.Z3.LT.0.D0.AND.Z4.LT.0.D0).
     .OR.(Z1.GT.0.D0.AND.Z2.GT.0.D0.AND.Z3.GT.0.D0.AND.Z4.GT.0.D0)) THEN
          A=DABS(Z1)
          B=(DABS(Z2)-A)/ALATO 
          C=(DABS(Z3)-A)/ALATO 
          D=(DABS(Z4)-A-B*ALATO-C*ALATO)/ALATO**2  
          VOL=A*ALATO**2+(B+C)*ALATO**3/2.D0+D*ALATO**4/4.D0
          IF(Z1.GT.0.D0) THEN  
          NPIU=NPIU+1  
          VPIU=VPIU+DABS(VOL)  
          XP=XP+DABS(VOL)*XG
          YP=YP+DABS(VOL)*YG
          ELSE 
          NMENO=NMENO+1
          VMENO=VMENO+DABS(VOL)
          XM=XM+DABS(VOL)*XG
          YM=YM+DABS(VOL)*YG
          ENDIF
          ENDIF
40        CONTINUE 
          DO 50 I=2,K1 
          IJ=(I-1)*K+1 
          IJ2=IJ+1 
          I2J=I*K+1
          I2J2=I2J+1
          Z1=Z(IJ) 
          Z2=Z(I2J)
          Z3=Z(IJ2)
          Z4=Z(I2J2)
          XG=((I-1)+0.5D0)*ALATO  
          YG=ALATO/2.D0  
          IF((Z1.LT.0.D0.AND.Z2.LT.0.D0.AND.Z3.LT.0.D0.AND.Z4.LT.0.D0).
     .OR.(Z1.GT.0.D0.AND.Z2.GT.0.D0.AND.Z3.GT.0.D0.AND.Z4.GT.0.D0)) THEN
          A=DABS(Z1)
          B=(DABS(Z2)-A)/ALATO 
          C=(DABS(Z3)-A)/ALATO 
          D=(DABS(Z4)-A-B*ALATO-C*ALATO)/ALATO**2  
          VOL=A*ALATO**2+(B+C)*ALATO**3/2.D0+D*ALATO**4/4.D0
          IF(Z1.GT.0.D0) THEN  
          NPIU=NPIU+1  
          VPIU=VPIU+DABS(VOL)  
          XP=XP+DABS(VOL)*XG
          YP=YP+DABS(VOL)*YG
          ELSE 
          NMENO=NMENO+1
          VMENO=VMENO+DABS(VOL)
          XM=XM+DABS(VOL)*XG
          YM=YM+DABS(VOL)*YG
          ENDIF
          ENDIF
50        CONTINUE 
          DO 60 I=2,L1 
          DO 60 J=2,K1 
          IJ=(I-1)*K+J 
          IJ1=IJ-1 
          IJ2=IJ+1 
          I1J=(I-2)*K+J
          I1J1=I1J-1
          I1J2=I1J+1
          I2J=I*K+J
          I2J1=I2J-1
          I2J2=I2J+1
          PI=(Z(I2J)-Z(I1J))/2.D0  
          PJ=(Z(IJ2)-Z(IJ1))/2.D0  
          PM=(PI*PJ)/DABS(PI*PJ)*DSQRT(PI**2+PJ**2)
          ANG=DATAN2(PJ,PI)
          CII=Z(I2J)-2.D0*Z(IJ)+Z(I1J) 
          CJJ=Z(IJ2)-2.D0*Z(IJ)+Z(IJ1) 
          CIJ=Z(I2J2)-Z(I2J1)-Z(I1J2)+Z(I1J1)  
          S=(CII+CJJ)/2.D0 
          R=DSQRT((CII-CJJ)**2+4.D0*CIJ**2)/2.D0
          C1=S-R
          C2=S+R
          CM=(C1*C2)/DABS(C1*C2)*DSQRT(DABS(C1)*DABS(C2))  
          VOL=0.D0 
          Z1=Z(IJ) 
          Z2=Z(I2J)
          Z3=Z(IJ2)
          Z4=Z(I2J2)
          XG=((I-1)+0.5D0)*ALATO  
          YG=((J-1)+0.5D0)*ALATO  
          IF((Z1.LT.0.D0.AND.Z2.LT.0.D0.AND.Z3.LT.0.D0.AND.Z4.LT.0.D0).
     .OR.(Z1.GT.0.D0.AND.Z2.GT.0.D0.AND.Z3.GT.0.D0.AND.Z4.GT.0.D0)) THEN
          A=DABS(Z1)
          B=(DABS(Z2)-A)/ALATO 
          C=(DABS(Z3)-A)/ALATO 
          D=(DABS(Z4)-A-B*ALATO-C*ALATO)/ALATO**2  
          VOL=A*ALATO**2+(B+C)*ALATO**3/2.D0+D*ALATO**4/4.D0
          IF(Z1.GT.0.D0) THEN  
          NPIU=NPIU+1  
          VPIU=VPIU+DABS(VOL)  
          XP=XP+DABS(VOL)*XG
          YP=YP+DABS(VOL)*YG
          ELSE 
          NMENO=NMENO+1
          VMENO=VMENO+DABS(VOL)
          XM=XM+DABS(VOL)*XG
          YM=YM+DABS(VOL)*YG
          ENDIF
          ENDIF
          IF(VOL.NE.0.D0) THEN 
          WRITE(VOLUME,1230) VOL  
          ELSE 
          VOLUME=BLANK 
          ENDIF
          WRITE(6,1220) I,J,Z(IJ),PI,PJ,PM,ANG*CONV,
     .CII,CJJ,CIJ,C1,C2,CM,VOLUME  
          IF(IFILE.EQ.1) WRITE(2,1210) I,J,Z(IJ),PI,PJ,PM, 
     .ANG*CONV,CII,CJJ,CIJ,C1,C2,CM,VOLUME  
60        CONTINUE 
          IF(IFILE.EQ.1) WRITE(2,1010) 
          XP=XP/VPIU
          YP=YP/VPIU
          XM=XM/VMENO  
          YP=YM/VMENO  
          IF(ILI.EQ.1) THEN
          WRITE(6,1280) VPIU,NPIU,XP,YP,VMENO,NMENO,XM,YM
          ELSE  
          WRITE(6,1281) VPIU,NPIU,XP,YP,VMENO,NMENO,XM,YM
          ENDIF  
          STOP 
1000      FORMAT(5X,'TITOLO')  
1001     FORMAT(5X,'TITLE')  
1010      FORMAT(20A4) 
1020      FORMAT(///1X,20A4)
1050      FORMAT(5X,'IFILE:',/15X,'SALVA RISULTATI = 1', 
     .5X,'(DATI OUTPUT FILE 2.)',/15X,'NON SALVA RISULTATI = 0')
1051      FORMAT(5X,'IFILE:',/15X,'SAVE THE RESULTS = 1', 
     .5X,'(OUTPUT DATA ON FILE 2.)',/15X,'NO SAVE THE RESULTS = 0')
1100      FORMAT(5X,'ALATO:',  
     ./15X,'LATO DELLA MAGLIA QUADRATA')
1101      FORMAT(5X,'ALATO:',  
     ./15X,'SIDE OF THE SQUARE GRID')
1110      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y',  
     ./15X,'SCALA COMPONENTE W')
1111      FORMAT(5X,'SCALE:',/15X,'SCALES OF X Y COORDINATES',  
     ./15X,'SCALE OF W COMPONENT')
1140      FORMAT(5X,'FORMATO DATI INPUT')  
1141      FORMAT(5X,'FORMAT OF INPUT DATA')  
1150      FORMAT(5X,'DATI INPUT FILE 1.')  
1151      FORMAT(5X,'INPUT DATA ON FILE 1.')  
1200      FORMAT(///5X,'DEM (QUOTA, PENDENZE, CURVATURE E VOLUME)'///) 
1201      FORMAT(///5X,'DEM (HEIGHT, SLOPES, CURVATURES AND VOLUME)'///) 
1210      FORMAT(2I5,3F10.3,F10.4,6F10.3,A10)  
1220      FORMAT(5X,2I5,3F10.3,F10.4,6F10.3,A10)
1230      FORMAT(F10.3)
1250      FORMAT(5X,'POS. PUNTI',5X,'QUOTA',6X,'P(X)',6X,'P(Y)',4X,
     .'P(MED)',4X,'ANGOLO',5X,'C(XX)',5X,'C(YY)',5X,'C(XY)',
     .4X,'C(MAX)',4X,'C(MIN)',4X,'C(MED)',4X,'VOLUME'/)
1251      FORMAT(5X,'POINT POS.',4X,'HEIGHT',6X,'S(X)',6X,'S(Y)',4X,
     .'S(MED)',5X,'ANGLE',5X,'C(XX)',5X,'C(YY)',5X,'C(XY)',
     .4X,'C(MAX)',4X,'C(MIN)',4X,'C(MED)',4X,'VOLUME'/)
1280      FORMAT(///5X,'VOLUMI POSITIVI',  
     .5X,'NUMEROSITA',14X,'COORDINATE BARICENTRO',  
     ./8X,'NEGATIVI'//1X,F15.3,5X,I10,5X,2F15.3,
     ./1X,F15.3,5X,I10,5X,2F15.3)  
1281     FORMAT(///5X,'POSITIVE VOLUMES',  
     .7X,'QUANTITY',6X,'COORDINATES OF THE BARICENTER',  
     ./5X,'NEGATIVE VOLUMES'//6X,F15.3,5X,I10,5X,2F15.3,
     ./6X,F15.3,5X,I10,5X,2F15.3)  
          END  
