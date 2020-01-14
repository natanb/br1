C*****************************     SPLINP2
          SUBROUTINE SPLBIL(X,Y,Z,A,IA,P,B,KA,JA,LA,IB,JB,
     .C,T,JC,LC,D,W,R,U,ISRV,NNG,RES,FORMA,STP,
     .MMAX,NMAX,SCX,SCY,SCU,KK,LL,TOLL,
     .IAUTO,ITABLE,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,AST,PTO,ASTI
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),Z(1),A(MMAX,1),IA(MMAX,1),KA(1),B(1),
     .P(1),JA(1),LA(1),IB(1),JB(1),
     .C(1),T(1),JC(1),LC(1),D(1),
     .W(1),R(1),U(1),ISRV(1),NNG(1),RES(1),FORMA(20),
     .DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA AST/'*'/,PTO/'.'/,ASTI/' '/
C
          REWIND 1
          DO 5 I=1,NMAX
          NNG(I)=0
5         CONTINUE
          X0= 1.D30
          X9=-1.D30
          Y0= 1.D30
          Y9=-1.D30
          EU=0.D0
          I=0
10        READ(1,1010) (DATO(L),L=1,120)
          READ(CDATO,FORMA) NPT,XX,YY,US
          IF(NPT.EQ.0) GOTO 15
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN
          EU=EU+US
          X0=DMIN1(X0,XX)
          X9=DMAX1(X9,XX)
          Y0=DMIN1(Y0,YY)
          Y9=DMAX1(Y9,YY)
          I=I+1
	  X(I)=XX
	  Y(I)=YY
	  Z(I)=US
          ENDIF
          GOTO 10
15        M=I
          EU=EU/M
          IF(M.GT.MMAX) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,1190)
          ELSE
          WRITE(6,1191)
          ENDIF
          STOP
          ENDIF
C
          IF(IAUTO.NE.0) THEN
          VOL=(X9-X0)*(Y9-Y0)
          UNI=DSQRT(VOL*SCX*SCY/M)
          STP=IAUTO*UNI
          ENDIF
          X0=X0-STP/(1000.D0*SCX)
          Y0=Y0-STP/(1000.D0*SCY)
          NCX=IDINT((X9-X0)*SCX/STP)+1
          NCY=IDINT((Y9-Y0)*SCY/STP)+1
c          X0=X0-STP/SCX
c          Y0=Y0-STP/SCY
          NTC=NCX*NCY
          NIX=NCX+1
          NIY=NCY+1
          NTI=NIX*NIY
          IF(ILI.EQ.1) THEN
          WRITE(6,1170) NCX,NCY,NTC,NTI,STP,X0,Y0
          ELSE
          WRITE(6,1170) NCX,NCY,NTC,NTI,STP,X0,Y0
          ENDIF

	  DO I=1,M
c          US=Z(I)-EU
	  us=z(i)
          P(I)=1.D0
          B(I)=-US*SCU
          KA(I)=4
          XX=(X(I)-X0)*SCX
          YY=(Y(I)-Y0)*SCY
          NX=IDINT(XX/STP)+1
          NY=IDINT(YY/STP)+1
          NC=(NY-1)*NCX+NX
          NNG(NC)=NNG(NC)+1
	  ISRV(I)=NC
          EPS=(XX-(NX-1)*STP)/STP
          ETA=(YY-(NY-1)*STP)/STP
          CALL COORCLAS(NI,NX,NY,NIX)
	  IA(I,1)= NI
	  IA(I,2)= NI+1
	  IA(I,3)= NI+NIX
	  IA(I,4)= NI+NIX+1
	  A(I,1) = 1.D0 - EPS - ETA + EPS*ETA
	  A(I,2) = EPS - EPS*ETA
	  A(I,3) = ETA - EPS*ETA
	  A(I,4) = EPS*ETA

          ENDDO 
          N=NTI

          IF(ILI.EQ.1) THEN
          WRITE(6,1110)
          ELSE
          WRITE(6,1111)
          ENDIF
          DO 65 K=1,NTC
          J=(K-1)/NCX  
          WRITE(6,1115) K,NNG(K),K+J,K+J+1,K+NIX+J,K+NIX+J+1
65        CONTINUE 

	  CALL CLASSI(ISRV,LA,JA,M,NTC)

          NV=0
          DO 80 I=1,NIX
          NV=NV+1
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=10000.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
	b(m)=-z(ipnear)*scu
80        CONTINUE 

	  DO 90 I=NIX+1,N-2*NIX+1,NIX
          NV=NV+1
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=10000.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
	b(m)=-z(ipnear)*scu
90        CONTINUE 

          DO 100 I=2*NIX,N-NIX,NIX 
          NV=NV+1  
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=10000.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
	b(m)=-z(ipnear)*scu
100       CONTINUE 
	
          DO 110 I=(NIY-1)*NIX+1,N 
          NV=NV+1  
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=10000.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
	b(m)=-z(ipnear)*scu
110       CONTINUE 


          DO 120 I=2,NIY-1 
          DO 120 J=2,NIX-1
          N1=(I-2)*NCX+(J-2)+1
          N2=N1+1  
          N3=N1+NCX
          N4=N3+1  
	  NNNN=NNG(N1)+NNG(N2)+NNG(N3)+NNG(N4)
          IF(NNNN.LT.4.AND.NNNN.GT.0) THEN
	  CALL COORCLAS(NI,J,I,NIX)
	  CALL VINMEDIA(N1,N2,N3,N4,LA,JA,Z,ELZ)  
          NV=NV+1  
	  M=M+1
          A(M,1)=1.D0
	  IA(M,1)=NI
	  KA(M)=1
	  P(M)=100.D0
c         B(M)=-(ELZ-EU)*SCU  
	b(m)=-elz*scu
C 
C CALCOLA LA MEDIA DEI PUNTI VICINI E VINCOLA IL NODO
C
          ENDIF

  	  IF(NNNN.EQ.0) THEN
	  CALL VININTER(N1,N2,N3,N4,LA,JA,IB,NCX,NCY,NP)
	  CALL COORCLAS(NI,J,I,NIX)
	  IPNEAR=NEARPOIN(NI,NCX,X,Y,X0,Y0,STP,IB,NP)
          NV=NV+1  
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=NI
	  KA(M)=1
	  P(M)=10000.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
	b(m)=-z(ipnear)*scu
	  ENDIF

120       CONTINUE 
       
          CALL NORMAI(A,IA,KA,JA,LA,IB,JB,B,P,C,JC,LC,D,W,
     .ISRV,MMAX,M,N)
          DO 70 K=1,N  
          ISRV(K)=K
70        CONTINUE
          ISRV(N+1)=ILI
          IF(ITABLE.EQ.1) CALL WGRAFO(ISRV,JA,LA,IB,JB,N,0)
          IF(ITABLE.EQ.1) CALL WGRAFO(ISRV,JC,LC,IB,JB,N,1)
          NN=LC(N+1)-1
          IF(ILI.EQ.1) THEN 
          WRITE(6,1020) M,N,NN/N,N*(N+1)/2,NN,N*(N+1)/(2*NN)
          ELSE
          WRITE(6,1021) M,N,NN/N,N*(N+1)/2,NN,N*(N+1)/(2*NN)
          ENDIF

          F1=FMATI1(IMAX,C,JC,LC,W,N)  
          CALL ALGITE(C,T,JC,LC,D, 
     .JA,LA,IB,JB, 
     .1.D-5,100,ITE,1,1,R,U,W,N)
          F2=FMATI2(IMAX,T,JA,JC,LA,W,N)
          RAP=F1/F2
          IF(ILI.EQ.1) THEN
          WRITE(6,1125) ITE
          WRITE(6,1130) RAP,IMAX
          ELSE
          WRITE(6,1126) ITE
          WRITE(6,1131) RAP,IMAX
          ENDIF
          CALL SZEROR(A,IA,KA,B,P,D,SZ,SZ2,MMAX,M,N,NV)
          CALL DQVVI(A,IA,KA,T,JA,JC,LA,P,W,MMAX,M)
C  
          REWIND 1 
          IF(ILI.EQ.1) THEN
          WRITE(6,1030)
          ELSE
          WRITE(6,1031)
          ENDIF
          DO 130 I=1,N 
          IS=LC(I+1)-1 
          IF(SZ2*C(IS).GT.1.D-10) C(IS)=DSQRT(SZ2*C(IS))
          WRITE(6,1040) I,D(I),C(IS)
130       CONTINUE 
          IF(IFILE.EQ.1) THEN
	eu=0.d0
          WRITE(3,1015) X0,Y0,STP,EU,N,NIX
          WRITE(3,1018) (D(I),I=1,N)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1050)
          ELSE
          WRITE(6,1051)
          ENDIF
          REWIND 1 
          I=0  
140       READ(1,1010) (DATO(L),L=1,120)
          READ(CDATO,FORMA) NPT,X,Y,US
          IF(NPT.EQ.0) GOTO 150
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          I=I+1
          UU=US
          UU=UU*SCU
          SGNL=UU+B(I) 
          SQMV=0.D0
          ASTI=' ' 
          IF(SZ2*P(I).GT.1.D-10) THEN  
          SQMV=DSQRT(SZ2*P(I)) 
          IF(DABS(B(I))/SQMV.GE.TOLL) ASTI='*' 
          ENDIF
          SRV=1.D0-P(I)
          SQMS=0.D0
          IF(SZ2*SRV.GT.1.D-10) SQMS=DSQRT(SZ2*SRV)
          WRITE(6,1060) NPT,US,SGNL,SQMS,B(I),SQMV,ASTI
          ENDIF
          RES(I)=B(I)
          GOTO 140 
150       NGL=M-N+NV
          IF(ILI.EQ.1) THEN
          WRITE(6,1070) SZ,M,N,NV,NGL  
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,1108) ALFA
          ELSE
          WRITE(6,1071) SZ,M,N,NV,NGL  
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,1109) ALFA
          ENDIF
          WRITE(6,1005)
          RETURN
1000      FORMAT(20A4) 
1005      FORMAT(///5X,20A4)
1010      FORMAT(120A1)
1015      FORMAT(4D20.12,2I5)  
1018      FORMAT(6D20.12)  
1020      FORMAT(///5X,'OCCUPAZIONE DI MEMORIA'//  
     .5X,'NUMERO DATI VALIDI',5X,I5/
     .5X,'NUMERO NODI       ',5X,I5/
     .5X,'DENSITA MEDIA     ',5X,I5//  
     .39X,'PIENA',11X,'COMPATTA',7X,'FATTORE  DI'/70X,'RIEMPIMENTO'/
     .5X,'MATRICE NORMALE      ',8X,I10,9X,I10,10X,'1:',1X,I5) 
1021     FORMAT(///5X,'MEMORY STORAGE'//  
     .5X,'NUMBER OF DATA ',5X,I5/
     .5X,'NUMBER OF KNOTS',5X,I5/
     .5X,'AVERAGE DENSITY',5X,I5//  
     .36X,'FULL',14X,'SPARSE',11X,'FILL-IN'/
     .5X,'NORMAL MATRIX  ',10X,I10,10X,I10,10X,'1:',1X,I5) 
1030      FORMAT(///8X,'NUMERO DEL NODO',9X,'COEFFICENTE',7X,'SQM',/)  
1031      FORMAT(///5X,'NUMBER OF THE KNOT',10X,'COEFFICENT',8X,'SD',/)  
1040      FORMAT(13X,I10,5X,F15.6,F10.3)
1050      FORMAT(///12X,'NOME PUNTO',  
     .5X,'COMPONENTE',5X,'COMPONENTE',7X,'SQM',8X,'RESIDUO',7X,'SQM'
     ./28X,'OSSERVATA',8X,'STIMATA',/)
1051      FORMAT(///5X,'NAME OF THE POINT', 
     .7X,'OBSERVED',7X,'EXPECTED',8X,'SD',7X,'RESIDUAL',8X,'SD'
     ./28X,'COMPONENT',6X,'COMPONENT',/)
1060      FORMAT(12X,I10,F15.6,F15.6,F10.3,F15.6,F10.3,4X,A1)
1070      FORMAT(////5X,'SIGMA ZERO   ',F10.3,//
     ./5X,'NUMERO EQUAZIONI  ',I5,/5X,'NUMERO INCOGNITE  ',I5, 
     ./5X,'NUMERO VINCOLI    ',I5,/5X,'GRADI DI LIBERTA  ',I5) 
1071      FORMAT(////5X,'SIGMA ZERO   ',F10.3,//
     ./5X,'NUMBER OF EQUATIONS  ',I5,/5X,'NUMBER OF UNKNOWNS   ',I5, 
     ./5X,'NUMBER OF CONSTRAINTS',I5,/5X,'DEGREES OF FREEDOM   ',I5) 
1108      FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEL TEST '
     .'DI REIEZIONE DEGLI ERRORI GROSSOLANI: ALFA = ',F5.3, 
     .' (SU DUE CODE)')
1109     FORMAT(/5X,'SIGNIFICANCE LEVEL OF THE TEST '
     .'FOR OUTLIER REJECTION: ALFA = ',F5.3, 
     .' (ON TWO SIDES)')
1110      FORMAT(///5X,'CLASSE -',10X,'NODI COINVOLTI DALLA CLASSE'/
     .5X,'NUMEROSITA',8X,'(VERTICI DEL QUADRATO)'/)
1111      FORMAT(///7X,'CLASS -',11X,'KNOTS CONNECTED THE CLASS'/
     .7X,'QUANTITY',10X,'(VERTICES OF THE SQUARE)'/)
1115      FORMAT(5X,2I5,15X,4I5)
1125      FORMAT(//5X,'NUMERO DI ITERAZIONI (METODO DEL GRADIENTE ',
     .'CONIUGATO)',I5//)  
1126     FORMAT(//5X,'ITERATION NUMBER (CONJUGATE GRADIENT METHOD)',
     .I5//)  
1130      FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
1131      FORMAT(///5X,'CONDITION NUMBER', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
1170      FORMAT(///5X,'NUMERO CLASSI IN X E Y    ',5X,2I5,
     .//5X,'NUMERO TOTALE CLASSI      ',10X,I5,
     .//5X,'NUMERO NODI               ',10X,I5,
     .//5X,'PASSO DELLA INTERPOLAZIONE',5X,F10.3, 
     .//5X,'COORDINATE ORIGINE        ',F15.3,5X,F15.3//) 
1171      FORMAT(///5X,'NUMBER OF CLASSES IN X AND Y',5X,2I5,
     .//5X,'NUMBER OF CLASSES      ',10X,I5,
     .//5X,'NUMBER OF KNOTS        ',10X,I5,
     .//5X,'STEP    ',20X,F10.3, 
     .//5X,'INITIAL COORDINATES    ',F15.3,5X,F15.3//) 
1190      FORMAT(//5X,'NUMERO DEI DATI',
     ./5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'//)
1191      FORMAT(//5X,'NUMBER OF DATA',
     ./5X,'EXCEEDING THE PREVIOUS DIMENSION'//)
          END  
