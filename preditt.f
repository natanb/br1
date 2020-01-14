
          PARAMETER (NMAX=10000,MMAX=3*NMAX,NMAX1=NMAX+1,NSMAX=500,  
     .NPMAX=100) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 TITOLO,DATO,PCHAR,FORMA,XYZ,BLANK,PTO
          CHARACTER*120 CDATO,CFORM
          DIMENSION TEMPO(NMAX),NTIPO(NMAX),KA(NMAX1),COEF0(3,NPMAX),  
     .COEF1(3,NSMAX),CSPS(3,MMAX),Z(MMAX),U(3),U0(3),U1(3),U2(3),UC(3),
     .UD(3),UMAX(3),SMUD2(3),STP(3),T0(3),T1(3),NGRAD(3),NI(3),IND(7), 
     .XYZ(3),IPUNT(NMAX),NCASO(6),DA(6),DB(6),DC(6),SN(6),TGL(6),EC(3),
     .ES(3),NPC(3),ISEG(3,7),TITOLO(80),FORMA(80),DATO(121),PCHAR(120),
     .IPTO(6),ISTI(3,7),ISCO(7)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,BLANK/' '/
          DATA IND/1,1,2,1,2,2,3/  
          DATA XYZ/'X','Y','Z'/
          DATA ((ISEG(I,J),J=1,7),I=1,3) /1,0,1,0,1,0,1,
     .0,1,-1,0,0,1,-1,0,0,0,1,-1,-1,1/
          DATA ((ISTI(I,J),J=1,7),I=1,3) /1,0,-1,0,1,0,-1,
     .0,1,1,0,0,1,1, 0,0,0,-1,-1,-1,-1/
          DATA ISCO/-1,-1,1,1,-1,-1,1/
C
          CALL LINGUA(ILI)
          IF(ILI.EQ.1) THEN
          CALL PAPRET(TITOLO,ITIPO,IFILE,IPOLI,ISPLI,ICOLL,FORMA)
          ELSE
          CALL PIPRET(TITOLO,ITIPO,IFILE,IPOLI,ISPLI,ICOLL,FORMA)
          ENDIF
          WRITE(6,1020) TITOLO 
          IF(IPOLI.EQ.1) THEN  
          DO 12 I=1,ITIPO  
          READ(8,1095) T0(I),NGRAD(I)  
          N=NGRAD(I)+1 
          IF(N.GT.NPMAX) THEN
          IF(ILI.EQ.1) THEN  
          WRITE(6,1210)
          ELSE
          WRITE(6,1211)
          ENDIF
          STOP 
          ENDIF
          READ(8,1090) (COEF0(I,J),J=1,N)  
12        CONTINUE 
          ENDIF
          IF(ISPLI.EQ.1) THEN  
          DO 15 I=1,ITIPO  
          READ(9,1015) T1(I),STP(I),ES(I),NI(I)  
          IF(NI(I).GT.NSMAX) THEN
          IF(ILI.EQ.1) THEN  
          WRITE(6,1220)
          ELSE
          WRITE(6,1221)
          ENDIF
          STOP 
          ENDIF
          READ(9,1090) (COEF1(I,J),J=1,NI(I))  
15        CONTINUE 
          ENDIF
          IF(ICOLL.EQ.1) THEN  
          READ(10,*) M 
          IF(M.GT.MMAX) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,1230)
          ELSE 
          WRITE(6,1231)
          ENDIF
          STOP 
          ENDIF
          READ(10,1090) (Z(I),I=1,M)
          READ(10,*) N 
          IF(N.GT.NMAX) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,1240)
          ELSE
          WRITE(6,1241)
          ENDIF
          STOP 
          ENDIF
          READ(10,1095) (TEMPO(I),NTIPO(I),I=1,N)
          NTER=0
          IF(ITIPO.EQ.1) NTER=1
          IF(ITIPO.EQ.2) NTER=3
          IF(ITIPO.EQ.3) NTER=6
          DO 20 I=1,NTER
          READ(2,1000) (DATO(L),L=1,120)
          IF(I.EQ.1) READ(CDATO,1060) EU  
          IF(I.EQ.2) READ(CDATO,1060) EV  
          IF(I.EQ.4) READ(CDATO,1060) EW  
          READ(2,1050) NCASO(I),DA(I),DB(I),DC(I),SN(I),TGL(I) 
          SN(I)=0.D0
20        CONTINUE 
          EC(1)=EU 
          EC(2)=EV 
          EC(3)=EW 
          XCOVT=0.D0
          DO 22 I=1,NTER
          XCOVT=DMAX1(XCOVT,TGL(I))
22        CONTINUE
          KA(1)=0  
          DO 25 I=1,N  
          KA(I+1)=KA(I)+IND(NTIPO(I))  
          IPUNT(I)=I
25        CONTINUE 
          ENDIF

C
          DO 27 I=1,120
          PCHAR(I)=BLANK
27        CONTINUE
          J=0
28        J=J+1
          READ(1,1000) (DATO(I),I=1,120)
          K=0
          DO 30 I=1,120
          IF(DATO(I).NE.BLANK) K=K+1
30        CONTINUE
          IF(K.NE.0) THEN
          DO 32 I=1,120
          IF(DATO(I).EQ.PTO) PCHAR(I)=PTO
32        CONTINUE
          GOTO 28
          ENDIF
          NPT=J-1
          J=0
          DATO(121)=BLANK
          DO 33 I=1,6
          IPTO(I)=121
33        CONTINUE
          DO 35 I=1,120
          IF(PCHAR(I).EQ.PTO) THEN
          J=J+1
          IPTO(J)=I
          ENDIF
35        CONTINUE
C                
          REWIND 1
          DO 38 I=1,ITIPO
          NPC(I)=0
          UMAX(I)=0.D0
          SMUD2(I)=0.D0
38        CONTINUE
          IF(ILI.EQ.1) THEN
          IF(ITIPO.EQ.1) WRITE(6,1115)
          IF(ITIPO.EQ.2) WRITE(6,1112)
          IF(ITIPO.EQ.3) WRITE(6,1110)
          ELSE
          IF(ITIPO.EQ.1) WRITE(6,1116)
          IF(ITIPO.EQ.2) WRITE(6,1113)
          IF(ITIPO.EQ.3) WRITE(6,1111)
          ENDIF
          DO 100 J=1,NPT
          READ(1,1000) (DATO(I),I=1,120)
          IF(DATO(IPTO(2)).EQ.BLANK.AND.DATO(IPTO(3)).EQ.BLANK.AND.
     .DATO(IPTO(4)).EQ.BLANK) THEN
          READ(CDATO,FORMA) NP,TP
          ELSE 
          IF(ITIPO.EQ.1) READ(CDATO,FORMA) NP,TP,UC(1) 
          IF(ITIPO.EQ.2) READ(CDATO,FORMA) NP,TP,UC(1),UC(2)
          IF(ITIPO.EQ.3) READ(CDATO,FORMA) NP,TP,UC(1),UC(2)
     .,UC(3) 
          ENDIF
          DO 48 I=1,ITIPO  
          U0(I)=0.D0
          U1(I)=0.D0
          U2(I)=0.D0
48        CONTINUE 
          IF(IPOLI.EQ.1) THEN  
          CALL PREPLT(ITIPO,TP,T0,COEF0,NGRAD,U0)  
          ENDIF
          IF(ISPLI.EQ.1) THEN  
          CALL PRESPT(ITIPO,TP,T1,STP,COEF1,U1)
          ENDIF
          IF(ICOLL.EQ.1) THEN  
          CALL PRECLT(ITIPO,TP,XCOVT,Z,TEMPO,CSPS,IPUNT,NTIPO,KA,  
     .NMAX,N,DA,DB,DC,SN,NCASO,U2,ILI) 
          ENDIF
          IRIC=IPOLI+ISPLI*2+ICOLL*4
          DO 50 I=1,ITIPO
          U(I)=ISEG(1,IRIC)*U0(I)+ISEG(2,IRIC)*(U1(I)+ES(I))+
     .ISEG(3,IRIC)*(U2(I)+EC(I))
50        CONTINUE 
          IF(DATO(IPTO(2)).EQ.BLANK.AND.DATO(IPTO(3)).EQ.BLANK.AND.
     .DATO(IPTO(4)).EQ.BLANK) THEN
          WRITE(6,1120) NP,TP,(U(I),I=1,ITIPO)
          ELSE 
          DO 60 I=1,ITIPO
          IF(DATO(IPTO(1+I)).EQ.PTO) THEN  
          UD(I)=ISTI(1,IRIC)*U0(I)+ISTI(2,IRIC)*(U1(I)+ES(I))+
     .ISTI(3,IRIC)*(U2(I)+EC(I)) + ISCO(IRIC)*UC(I)
C          UD(I)=U(I)-UC(I)
          NPC(I)=NPC(I)+1 
          UMAX(I)=DMAX1(UMAX(I),DABS(UD(I)))
          SMUD2(I)=SMUD2(I)+UD(I)**2
          ENDIF
60        CONTINUE 
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IPTO(3)).EQ.BLANK.AND.
     .DATO(IPTO(4)).EQ.BLANK) THEN
          WRITE(6,1120) NP,TP,(U(I),I=1,ITIPO),UD(1)
          ENDIF
          IF(DATO(IPTO(2)).EQ.BLANK.AND.DATO(IPTO(3)).EQ.PTO.AND.
     .DATO(IPTO(4)).EQ.BLANK) THEN
          CFORM='(5X,I5,F15.3,$F12.3,12X,F12.3)'
          CFORM(15:15)=CHAR(ITIPO+48)
          WRITE(6,CFORM) NP,TP,(U(I),I=1,ITIPO),UD(2)
          ENDIF
          IF(DATO(IPTO(2)).EQ.BLANK.AND.DATO(IPTO(3)).EQ.BLANK.AND.
     .DATO(IPTO(4)).EQ.PTO) THEN
          CFORM='(5X,I5,F15.3,$F12.3,24X,F12.3)'
          CFORM(15:15)=CHAR(ITIPO+48)
          WRITE(6,CFORM) NP,TP,(U(I),I=1,ITIPO),UD(3)
          ENDIF
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IPTO(3)).EQ.PTO.AND.
     .DATO(IPTO(4)).EQ.BLANK) THEN
          WRITE(6,1120) NP,TP,(U(I),I=1,ITIPO),UD(1),UD(2)
          ENDIF
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IPTO(3)).EQ.BLANK.AND.
     .DATO(IPTO(4)).EQ.PTO) THEN
          CFORM='(5X,I5,F15.3,$F12.3,F12.3,12X,F12.3)'
          CFORM(15:15)=CHAR(ITIPO+48)
          WRITE(6,CFORM) NP,TP,(U(I),I=1,ITIPO),UD(1),UD(3)
          ENDIF
          IF(DATO(IPTO(2)).EQ.BLANK.AND.DATO(IPTO(3)).EQ.PTO.AND.
     .DATO(IPTO(4)).EQ.PTO) THEN
          CFORM='(5X,I5,F15.3,$F12.3,F12.3,12X,2F12.3)'
          CFORM(15:15)=CHAR(ITIPO+48)
          WRITE(6,CFORM) NP,TP,(U(I),I=1,ITIPO),UD(2),UD(3)
          ENDIF
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IPTO(3)).EQ.PTO.AND.
     .DATO(IPTO(4)).EQ.PTO) THEN
          WRITE(6,1120) NP,TP,(U(I),I=1,ITIPO),(UD(I),I=1,3)   
          ENDIF
          ENDIF    
          IF(IFILE.EQ.1) WRITE(11,1145) NP,TP,(U(I),I=1,ITIPO)
100       CONTINUE
          IF(ILI.EQ.1) THEN
          WRITE(6,1150)
          ELSE
          WRITE(6,1151)
          ENDIF
          DO 110 I=1,ITIPO 
          IF(NPC(I).GT.0) THEN
          WRITE(6,1155) NPC(I),DSQRT(SMUD2(I)/NPC(I)),UMAX(I),XYZ(I) 
          ENDIF
110       CONTINUE
          IF(IFILE.EQ.1) WRITE(11,1060)
          STOP 
1000      FORMAT(120A1)
1015      FORMAT(3D20.12,I5)
1020      FORMAT(///5X,120A1///)
1050      FORMAT(I5,5D15.8)
1060      FORMAT(5X,F15.6) 
1090      FORMAT(6D20.12)  
1095      FORMAT(D20.12,I5)
1110      FORMAT(///5X,'NOME ',10X,'TEMPO',3(2X,'COMPONENTE'),
     .4X,'SCARTO X',4X,'SCARTO Y',4X,'SCARTO Z' 
     ./5X,'PUNTO',17X,'STIMATA X ',2X,'STIMATA Y ',2X,'STIMATA Z'//)  
1111      FORMAT(///5X,'POINT',
     .11X,'TIME',2X,'X EXPECTED',2X,'Y EXPECTED', 
     .2X,'Z EXPECTED',2X,'X RESIDUAL',2X,'Y RESIDUAL',2X,'Z RESIDUAL' 
     ./5X'NAME ',15X,3(2X,'COMPONENTS')//)
1112      FORMAT(///5X,'NOME ',10X,'TEMPO',2(2X,'COMPONENTE'),
     .4X,'SCARTO X',4X,'SCARTO Y' 
     ./5X,'PUNTO',17X,'STIMATA X ',2X,'STIMATA Y '//) 
1113      FORMAT(///5X,'POINT',
     .11X,'TIME',2X,'X EXPECTED',2X,'Y EXPECTED', 
     .2X,'X RESIDUAL',2X,'Y RESIDUAL' 
     ./5X'NAME ',15X,2(2X,'COMPONENTS')//)
1115      FORMAT(///5X,'NOME ',10X,'TEMPO',2X,'COMPONENTE',
     .4X,'SCARTO X'/5X,'PUNTO',17X,'STIMATA X '//)
1116      FORMAT(///5X,'POINT',
     .11X,'TIME',2X,'X EXPECTED',2X,'X RESIDUAL' 
     ./5X'NAME ',15X,2X,'COMPONENTS'//)
1120      FORMAT(5X,I5,F15.3,6F12.3)
c1125      FORMAT(5X,I5,F15.3,<ITIPO>F12.3,12X,F12.3) 
c1130      FORMAT(5X,I5,F15.3,<ITIPO>F12.3,24X,F12.3) 
c1135      FORMAT(5X,I5,F15.3,<ITIPO>F12.3,F12.3,12X,F12.3)  
c1140      FORMAT(5X,I5,F15.3,<ITIPO>F12.3,12X,2F12.3)  
1145      FORMAT(I5,F15.3,3F15.5)
1150      FORMAT(///5X,'ANALISI STATISTICA',//8X,'N.',7X,'RMS',7X,'MAX',
     .//)  
1151      FORMAT(///5X,'STATISTICAL ANALYSIS',//8X,'N.',7X,'RMS',
     .7X,'MAX',//) 
1155      FORMAT(5X,I5,2F10.3,5X,A1,//)  
1210      FORMAT(//5X,'NUMERO COEFFICIENTI POLINOMIO', 
     ./5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO') 
1211      FORMAT(//5X,'NUMBER OF POLYNOMIAL COEFFICENTS', 
     ./5X,'EXCEEDING THE PREVIOUS DIMENSION') 
1220      FORMAT(//5X,'NUMERO COEFFICIENTI SPLINES',
     ./5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO') 
1221      FORMAT(//5X,'NUMBER OF SPLINE COEFFICENTS',
     ./5X,'EXCEEDING THE PREVIOUS DIMENSION') 
1230      FORMAT(//5X,'NUMERO ELEMENTI DELLA SOLUZIONE',
     ./5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO') 
1231      FORMAT(//5X,'NUMBER OF FILTERING COEFFICENTS',
     ./5X,'EXCEEDING THE PREVIOUS DIMENSION') 
1240      FORMAT(//5X,'NUMERO PUNTI FILTRATI', 
     ./5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO') 
1241      FORMAT(//5X,'NUMBER OF FILTERED POINTS', 
     ./5X,'EXCEEDING THE PREVIOUS DIMENSION') 
          END  
