          PARAMETER (NMAX=1000,KMAX=3,MMAX=KMAX*NMAX,NBMAX=300,
     .MBMAX=KMAX*NBMAX)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PCHAR,PTO,AST,BLANK 
          CHARACTER*4 FORMA,TITOLO 
          CHARACTER*120 CDATO
          DIMENSION V(MMAX),Z(MMAX),SRV(MMAX),C(MMAX,MBMAX),TEMPO(NMAX),
     .NOME(NMAX),NTIPO(NMAX),LASTP(NMAX),IPUNT(NMAX),DA(6),DB(6),DC(6),
     .SNOISE(6),NCASO(6),DATO(120),PCHAR(120),IPTO(4),IAST(3), 
     .IND(7),ICA(3,3,7),AA(3,3),BB(3,3),CC(3,3),D(3),E(3),EC(3), 
     .KP(3),LP(3),S1(3),S2(3),S3(3),S4(3),S5(3),S6(3),S7(3),S8(3), 
     .TGL(6),FORMA(20),TITOLO(20)  
          EQUIVALENCE (CDATO,DATO(1))
          DATA IND/1,1,2,1,2,2,3/  
          DATA PTO/'.'/,AST/'*'/,BLANK/' '/
          DATA IPTO/4*120/,IAST/116,118,120/
          DATA ((ICA(I,J,1),J=1,3),I=1,3)/1,0,0,0,0,0,0,0,0/
          DATA ((ICA(I,J,2),J=1,3),I=1,3)/2,0,0,0,0,0,0,0,0/
          DATA ((ICA(I,J,3),J=1,3),I=1,3)/1,3,0,0,2,0,0,0,0/
          DATA ((ICA(I,J,4),J=1,3),I=1,3)/4,0,0,0,0,0,0,0,0/
          DATA ((ICA(I,J,5),J=1,3),I=1,3)/1,5,0,0,4,0,0,0,0/
          DATA ((ICA(I,J,6),J=1,3),I=1,3)/2,6,0,0,4,0,0,0,0/
          DATA ((ICA(I,J,7),J=1,3),I=1,3)/1,3,5,0,2,6,0,0,4/
C    
          READ(5,*) ILI
          READ(5,1010) TITOLO  
          WRITE(6,1020) TITOLO 
          READ(5,*) ICASO,ITIPO
          READ(5,*) ISCRI,IFILE,IMATRI,TOLL
          READ(5,*) SCT,SCU,SCV,SCW
          READ(5,1010) FORMA
          CALL LEGEND(ICASO,ILI)
          IF(ILI.EQ.1) THEN
          WRITE(6,1070)
          ELSE
          WRITE(6,1071)
          ENDIF
          NTER=0
          IF(ITIPO.EQ.1) NTER=1
          IF(ITIPO.EQ.2) NTER=3
          IF(ITIPO.EQ.3) NTER=6
          DO 1 I=1,NTER
          READ(2,1000) (DATO(L),L=1,120)
          IF(I.EQ.1) READ(CDATO,1060) EU  
          IF(I.EQ.2) READ(CDATO,1060) EV  
          IF(I.EQ.4) READ(CDATO,1060) EW  
          READ(2,1050) NCASO(I),DA(I),DB(I),DC(I),SNOISE(I),TGL(I) 
          WRITE(6,1075) NCASO(I),DA(I),DB(I),DC(I),SNOISE(I),TGL(I)
1         CONTINUE 
 	  EC(1)=EU
	  EC(2)=EV
	  EC(3)=EW 
          XCOVT=0.D0
          DO 3 I=1,NTER
          XCOVT=DMAX1(XCOVT,TGL(I))
3         CONTINUE 
          DO 5 I=1,120 
          PCHAR(I)=BLANK
5         CONTINUE 
          J=0  
10        J=J+1   
          READ(1,1000) (DATO(I),I=1,120) 
          K=0
          DO 15 I=1,120
          IF(DATO(I).NE.BLANK) K=K+1
15        CONTINUE
          IF(K.NE.0) THEN 
          DO 20 I=1,120
          IF(DATO(I).EQ.PTO) PCHAR(I)=PTO  
20        CONTINUE 
          GOTO 10
          ENDIF 
30        N=J-1  
          J=0  
          DO 40 I=1,120
          IF(PCHAR(I).EQ.PTO) THEN 
          J=J+1
          IPTO(J)=I
          ENDIF
40        CONTINUE 
          NPTO=J
          IF(NPTO.NE.(ICASO+ITIPO)) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,1150)
          ELSE
          WRITE(6,1151)
          ENDIF
          STOP
          ENDIF
C           
          IF(ISCRI.EQ.1) THEN  
          REWIND 1 
          CALL TESTA(ICASO,ITIPO,0,IPTO,NPTO,ILI)
          DO 45 I=1,N  
          READ(1,1000) (DATO(J),J=1,120)
          WRITE(6,1005) (DATO(J),J=1,120)  
45        CONTINUE 
          ENDIF
C         
          REWIND 1
          READ(1,1000) (DATO(J),J=1,120)
          LAST=IPTO(1)
          DO 47 I=LAST+1,LAST+10
          IF(DATO(I).EQ.BLANK) GOTO 48
47        CONTINUE
48        LAST=I-1
C
          REWIND 1 
          K=0  
          J=0  
          ET=0.D0  
50        READ(1,1000)(DATO(I),I=1,120)
          IF(ITIPO.EQ.1) READ(CDATO,FORMA) NPT,TS,US  
          IF(ITIPO.EQ.2) READ(CDATO,FORMA) NPT,TS,US,VS
          IF(ITIPO.EQ.3) READ(CDATO,FORMA) NPT,TS,US,VS,WS
          IF(NPT.EQ.0) GOTO 60
          IP1=0
          IP2=0
          IP3=0
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          K=K+1
          V(K)=US*SCU-EU
          IP1=1
          ENDIF
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          K=K+1
          V(K)=VS*SCV-EV
          IP2=2
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(3)).NE.AST) THEN
          K=K+1
          V(K)=WS*SCW-EW
          IP3=4
          ENDIF
          IF((IP1+IP2+IP3).EQ.0) GOTO 50  
          J=J+1
          NOME(J)=NPT  
          IPUNT(J)=J
          ET=ET+TS*SCT 
          TEMPO(J)=TS*SCT  
          NTIPO(J)=IP1+IP2+IP3 
          GOTO 50 
60        N=J  
          M=K
          IF(ILI.EQ.1) THEN  
          WRITE(6,1008)
          ELSE
          WRITE(6,1009)
          ENDIF
          ET=ET/N  
          IF(ITIPO.EQ.1) WRITE(CDATO,FORMA) N,ET/SCT,EU/SCU
          IF(ITIPO.EQ.2) WRITE(CDATO,FORMA) N,ET/SCT,EU/SCU,EV/SCV 
          IF(ITIPO.EQ.3) WRITE(CDATO,FORMA) N,ET/SCT,EU/SCU,EV/SCV,
     .    EW/SCW
          WRITE(6,1005) (DATO(I),I=1,120)  
          IF(N.GT.NMAX.OR.M.GT.MMAX) THEN
          IF(ILI.EQ.1) THEN  
          WRITE(6,1130)
          ELSE
          WRITE(6,1131)
          ENDIF
          STOP 
          ENDIF
C   
          CALL BMAX(XCOVT,TEMPO,IPUNT,NTIPO,LASTP,NMAX,N,MB)
          IF(MB.GT.MBMAX) THEN
          IF(ILI.EQ.1) THEN 
          WRITE(6,1140)
          ELSE
          WRITE(6,1141)
          ENDIF
          STOP 
          ENDIF
          INR=0
          DO 100 I=1,N-1
          INC=0
          II=NTIPO(I)  
          IR=IND(II)
          DO 70 L=1,IR 
          DO 70 K=L,IR 
          IS=ICA(L,K,II)
          C(INR+L,INC+K-L+1)=SNOISE(IS)+DA(IS) 
70        CONTINUE 
          INC=INC+IR
          LNI=LASTP(I) 
          DO 90 J=I+1,LNI  
          JJ=NTIPO(J)  
          IC=IND(JJ)
          DELTA=DABS(TEMPO(I)-TEMPO(J))
          DO 80 L=1,IR 
          DO 80 K=1,IC 
          C(INR+L,INC+K-L+1)=GASS1(DELTA,XCOVT,II,JJ,L,K,  
     .DA,DB,DC,SNOISE,NCASO,ICASO)
80        CONTINUE 
          INC=INC+IC
90        CONTINUE 
          INR=INR+IR
100       CONTINUE 
          INC=0
          II=NTIPO(N)  
          IR=IND(II)
          DO 110 L=1,IR
          DO 110 K=L,IR
          IS=ICA(L,K,II)
          C(INR+L,INC+K-L+1)=SNOISE(IS)+DA(IS) 
110       CONTINUE 
          INR=INR+IR
          DO 120 I=1,M 
          Z(I)=V(I)
120       CONTINUE
          IF(ILI.EQ.1) THEN 
          WRITE(6,1110) N,M,MB,M**2,M*MB,M**2/(M*MB)
          ELSE
          WRITE(6,1111) N,M,MB,M**2,M*MB,M**2/(M*MB)
          ENDIF
          IF(IMATRI.EQ.1) THEN 
          WRITE(12,1090) ((C(I,J),J=1,MB),Z(I),I=1,M)
          ENDIF
C   
          F1=FMATB(IMAX,C,SRV,MMAX,M,MB)
          C(M+1,1)=DBLE(ILI)
          CALL TCHOLB(C,MMAX,M,MB) 
          CALL SOLSB(C,Z,MMAX,M,MB)
          CALL INVERB(C,SRV,MMAX,M,MB) 
          F2=FMATB(IMAX,C,SRV,MMAX,M,MB)
          RAP=F1/F2
          IF(ILI.EQ.1) THEN
          WRITE(6,1120) RAP,IMAX
          ELSE
          WRITE(6,1121) RAP,IMAX
          ENDIF
          IF(IFILE.EQ.1) THEN  
          WRITE(10,*) M
          WRITE(10,1090) (-Z(I),I=1,M) 
          WRITE(10,*) N
          WRITE(10,1095) (TEMPO(I),NTIPO(I),I=1,N) 
          ENDIF
C   
          REWIND 1
          IF(ILI.EQ.1) THEN
          WRITE(6,1080)
          ELSE
          WRITE(6,1081)
          ENDIF
          DO 125 I=1,3 
          KP(I)=0  
          LP(I)=0  
          S1(I)=0.D0
          S2(I)=0.D0
          S3(I)=0.D0
          S4(I)=0.D0
          S5(I)=0.D0
          S6(I)=0.D0
          S7(I)=0.D0
          S8(I)=0.D0
125       CONTINUE 
          INR=0
          DO 220 I=1,N
          READ(1,1000) (DATO(J),J=1,120)           
          NPT=NOME(I)
          DO 130 LL=1,3
          DO 130 JJ=1,3
          AA(LL,JJ)=0.D0
          BB(LL,JJ)=0.D0
          CC(LL,JJ)=0.D0
130       CONTINUE 
          II=NTIPO(I)  
          IR=IND(II)
          DO 150 L=1,IR
          DO 150 K=L,IR
          IS=ICA(L,K,II)
          AA(L,K)=SNOISE(IS)
          AA(K,L)=AA(L,K)  
          BB(L,K)=C(INR+L,K-L+1)
          BB(K,L)=BB(L,K)  
150       CONTINUE 
          DO 170 L=1,IR
          SUM=0.D0 
          DO 160 K=1,IR
          SUM=SUM-AA(L,K)*Z(INR+K) 
160       CONTINUE 
          E(L)=SUM
170       CONTINUE 
          DO 180 L=1,IR
          D(L)=V(INR+L)-E(L)+EC(L)
180       CONTINUE 
          INR=INR+IR
          CALL MAT3(AA,BB,CC)  
          DO 190 L=1,IR
          DO 190 K=L,IR
          CC(L,K)=AA(L,K)-BB(L,K)  
          IS=ICA(L,K,II)
          AA(L,K)=DA(IS)-CC(L,K)
190       CONTINUE 
          DO 210 L=1,IR
          AA(L,L)=DSQRT(AA(L,L))
          BB(L,L)=DSQRT(BB(L,L))
          KP(L)=KP(L)+1
          S1(L)=S1(L)+D(L) 
          S2(L)=S2(L)+D(L)**2  
          S3(L)=S3(L)+E(L) 
          S4(L)=S4(L)+E(L)**2  
          IF(DABS(E(L))/BB(L,L).LE.TOLL) THEN  
          LP(L)=LP(L)+1
          S5(L)=S5(L)+E(L) 
          S6(L)=S6(L)+E(L)**2  
          ENDIF
          S7(L)=S7(L)+CC(L,L)  
          S8(L)=DMAX1(S8(L),CC(L,L))
          CC(L,L)=DSQRT(CC(L,L))
          IF(L.GT.1) THEN  
          DO 200 K=1,L-1
	  SRW=0.D0
	  PRD=AA(L,L)*AA(K,K)
	  IF(DABS(PRD).GT.1.D-10) SRW=AA(K,L)/PRD
	  AA(K,L)=SRW
	  SRW=0.D0
	  PRD=BB(L,L)*BB(K,K)
	  IF(DABS(PRD).GT.1.D-10) SRW=BB(K,L)/PRD
	  BB(K,L)=SRW
	  SRW=0.D0
	  PRD=CC(L,L)*CC(K,K)
	  IF(DABS(PRD).GT.1.D-10) SRW=CC(K,L)/PRD
	  CC(K,L)=SRVW
200       CONTINUE 
          ENDIF
210       CONTINUE 
          CALL CODA(I,NPT,II,D,AA,E,BB,CC,IFILE,DATO,LAST,TOLL) 
          write(13,'(I5,f10.3,f15.8)') NPT,E(1),D(1) 
220       CONTINUE
          IF(ILI.EQ.1) THEN 
          WRITE(6,1100)
          ELSE
          WRITE(6,1101)
          ENDIF
          DO 230 I=1,3 
          K=KP(I)  
          IF(K.EQ.0) GOTO 230  
          S1(I)=S1(I)/K
          S2(I)=DSQRT((S2(I)/K-S1(I)**2)*K/(K-1))  
          S3(I)=S3(I)/K
          S4(I)=DSQRT((S4(I)/K-S3(I)**2)*K/(K-1))  
          S5(I)=S5(I)/LP(I)
          S6(I)=DSQRT((S6(I)/LP(I)-S5(I)**2)*LP(I)/(LP(I)-1))  
          S7(I)=DSQRT(S7(I)/K) 
          S8(I)=DSQRT(S8(I))
          WRITE(6,1105) I,K,K-LP(I),S1(I),S2(I),S3(I),S4(I),S5(I),S6(I),
     .S7(I),S8(I)
230       CONTINUE 
          ALFA=2.D0*(1.D0-DFN01(TOLL))
          IF(ILI.EQ.1) THEN 
          WRITE(6,1108) ALFA
          ELSE
          WRITE(6,1109) ALFA
          ENDIF
          IF(IFILE.EQ.1) WRITE(11,1005)
          STOP 
1000      FORMAT(120A1)
1005      FORMAT(5X,120A1) 
1008      FORMAT(///5X,'NUMERO DEI DATI - MEDIE DEL TEMPO ',
     .'E DELLE COMPONENTI'/)
1009     FORMAT(///5X,'NUMBER OF DATA - AVERAGE OF THE TIME ',
     .'AND THE COMPONENTS'/)
1010      FORMAT(20A4) 
1020      FORMAT(///5X,20A4)
1040      FORMAT(5X,3F10.3)
1050      FORMAT(I5,5D15.8)
1060      FORMAT(5X,F15.6) 
1070      FORMAT(//5X,'CASO:',4X,'VARIANZA O ',4X,'FATTORE DI ',
     .3X,'FATTORE DI  ',5X,'VARIANZA O',4X,'ASCISSA DI '/
     .14X,'COVARIANZA ',4X,'SMORZAMENTO',3X,'OSCILLAZIONE',
     .5X,'COVARIANZA',4X,'AZZERAMENTO'/  
     .14X,'DEL SEGNALE',35X,'DEL RUMORE'/)  
1071      FORMAT(//5X,'CASE:',4X,'VARIANCE OR  ',7X,'SMOOTING',
     .4X,'OSCILLATING',4X,'VARIANCE OR ',5X,'ZERO POINT'/
     .14X,'COVARIANCE   ',7X,'FACTOR  ',4X,'FACTOR     ',
     .4X,'COVARIANCE  ',5X,'COORDINATE'/  
     .14X,'OF THE SIGNAL',34X,'OF THE NOISE'/)  
1075      FORMAT(5X,I5,2X,F15.6,2D15.8,1X,F15.6,D15.8)
1080      FORMAT(///5X,'COLLOCAZIONE MIN. QUADRATI (FILTRAGGIO)'// 
     .8X,'N.',1X,'PNT.',7X,'SIGNAL X',5X,'SQM-S',1X,
     .'R XX',3X,'XY',3X,'XZ',8X,'NOISE X',10X,'SQM-N',1X,  
     .'R XX',3X,'XY',3X,'XZ',10X,'SQM-E',1X,
     .'R XX',3X,'XY',3X,'XZ',/29X,'Y',18X,'YY',3X,'YZ',
     .14X,'Y',23X,'YY',3X,'YZ',23X,'YY',3X,'YZ'/
     .29X,'Z',23X,'ZZ',14X,'Z',28X,'ZZ',28X,'ZZ'/) 
1081      FORMAT(///5X,'LEAST SQUARES COLLOCATION (FILTERING)'// 
     .8X,'N.',1X,'PNT.',7X,'SIGNAL X',6X,'SD-S',1X,
     .'R XX',3X,'XY',3X,'XZ',8X,'NOISE X',11X,'SD-N',1X,  
     .'R XX',3X,'XY',3X,'XZ',11X,'SD-E',1X,
     .'R XX',3X,'XY',3X,'XZ',/29X,'Y',18X,'YY',3X,'YZ',
     .14X,'Y',23X,'YY',3X,'YZ',23X,'YY',3X,'YZ'/
     .29X,'Z',23X,'ZZ',14X,'Z',28X,'ZZ',28X,'ZZ'/) 
1090      FORMAT(6D20.12)  
1095      FORMAT(D20.12,I5)
1100      FORMAT(///5X,'ANALISI STATISTICA'//  
     .5X,'COMP.',3X,'N.',4X,'-',6X,
     .'M(S)',4X,'SQM(S)',6X,'M(N)',4X,'SQM(N)',
     .5X,'M(N-)',3X,'SQM(N-)',4X,'RMS(E)',4X,'MAX(E)'/)
1101      FORMAT(///5X,'STATISTICAL ANALYSIS'//  
     .5X,'COMP.',3X,'N.',4X,'-',6X,
     .'M(S)',5X,'SD(S)',6X,'M(N)',5X,'SD(N)',
     .5X,'M(N-)',4X,'SD(N-)',4X,'RMS(E)',4X,'MAX(E)'/)
1105      FORMAT(5X,3I5,8F10.3)
1108      FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEL TEST '
     .'DI REIEZIONE DEGLI ERRORI GROSSOLANI: ALFA = ',F5.3, 
     .' (SU DUE CODE)')
1109      FORMAT(/5X,'SIGNIFICANCE LEVEL OF TEST '
     .'FOR OUTLIER REJECTION: ALFA = ',F5.3, 
     .' (ON TWO SIDES)')
1110      FORMAT(///5X,'OCCUPAZIONE DI MEMORIA'//  
     .5X,'NUMERO DATI VALIDI',5X,I5/
     .5X,'NUMERO COMPONENTI ',5X,I5//  
     .5X,'AMPIEZZA DI BANDA ',5X,I5//  
     .40X,'PIENA',12X,'COMPATTA',9X,'FATTORE  DI'/74X,'RIEMPIMENTO'/
     .5X,'MATRICE DI COVARIANZA',9X,I10,10X,I10,10X,'1:',1X,I5) 
1111     FORMAT(///5X,'MEMORY STORAGE'//  
     .5X,'NUMBER OF DATA       ',5X,I5/
     .5X,'NUMBER OF COMPONENTS ',5X,I5//  
     .5X,'BANDWIDTH            ',5X,I5//  
     .41X,'FULL',14X,'SPARSE',11X,'FILL-IN'/
     .5X,'COVARIANCE MATRIX',13X,I10,10X,I10,10X,'1:',1X,I5) 
1120      FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
1121      FORMAT(///5X,'CONDITION NUMBER', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
1130      FORMAT(///5X,'NUMERO DI ELEMENTI DELLA MATRICE DI COVARIANZA '
     .5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'/) 
1131      FORMAT(///5X,'NUMBER OF ELEMENTS OF THE COVARIANCE MATRIX '
     .5X,'EXCEEDING THE PREVIOUS DIMENSION'/) 
1140      FORMAT(///5X,'NUMERO DI PUNTI E/O COMPONENTI ECCEDENTE ',
     .'IL DIMENSIONAMENTO PREVISTO'/)  
1141      FORMAT(///5X,'NUMBER OF POINTS AND/OR COMPONENTS ',
     .'EXCEEDING THE PREVIOUS DIMENSION'/)  
1150      FORMAT(5X,'CASO E/O TIPO DEI DATI ERRATO'///)
1151      FORMAT(5X,'WRONG CASE AND/OR TYPE OF DATA'///)        
          END  

