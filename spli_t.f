C*****************************     SPLINT
          SUBROUTINE SPLINT(A,IA,P,B,KA,C,D,W,NNG,RES,FORMA,STP,
     .MMAX,NMAX,NBMAX,SCT,SCX,KK,LL,TOLL,IAUTO,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          CHARACTER*1 DATO,PTO,AST,ASTI
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION A(MMAX,1),IA(MMAX,1),P(1),B(1),KA(1),C(NMAX,1),
     .D(1),W(1),NNG(1),RES(1),FORMA(1),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/,ASTI/' '/ 
C
          REWIND 1 
          DO 5 I=1,NMAX
          NNG(I)=0 
5         CONTINUE 
          T0= 1.D30
          T9=-1.D30
          EX=0.D0  
          I=0
10        READ(1,1010) (DATO(L),L=1,120)
          READ(CDATO,FORMA) NPT,TS,XS
          IF(NPT.EQ.0) GOTO 15 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          T0=DMIN1(T0,TS)
          T9=DMAX1(T9,TS)
          I=I+1
          EX=EX+XS
          ENDIF
          GOTO 10  
15        M=I
          EX=EX/M
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
          VOL=(T9-T0)
          UNI=(VOL*SCT/M)
          STP=IAUTO*UNI
          ENDIF
          T0=T0-STP/(1000.D0*SCT)
          NTC=IDINT((T9-T0)*SCT/STP)+3
          T0=T0-STP/SCT
          NTI=NTC+1
          IF(ILI.EQ.1) THEN
          WRITE(6,1130) NTC,NTI,STP,T0
          ELSE
          WRITE(6,1131) NTC,NTI,STP,T0
          ENDIF
          REWIND 1 
          I=0  
20        READ(1,1010) (DATO(L),L=1,120)
          READ(CDATO,FORMA) NPT,TS,XS
          IF(NPT.EQ.0) GOTO 30 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          I=I+1
          XS=XS-EX
          P(I)=1.D0
          B(I)=-XS*SCX 
          KA(I)=4  
          TS=(TS-T0)*SCT  
          NT=IDINT(TS/STP)+1  
          NNG(NT)=NNG(NT)+1  
          SRV=(TS-(NT-1)*STP)/STP
          NI=NT-2  
          DO 25 J=1,4  
          AJ=DBLE(2-J) 
          JS=5-J
          A(I,J)=SPLT(SRV+AJ,JS)
          IA(I,J)=NI+J
25        CONTINUE 
          ENDIF
          GOTO 20  
30        N=NTI
          WRITE(6,1140)
          DO 32 K=1,N-1
          WRITE(6,1150) K,NNG(K),K,K+1 
32        CONTINUE  
          CALL NORMAB(A,IA,KA,B,P,C,D,MMAX,NMAX,NBMAX,M,N,NB)  
          C(1,1)=C(1,1)+1.D10  
          C(N,1)=C(N,1)+1.D10  
          NV=2 
          DO 35 I=2,N-1
          IF(NNG(I-1).EQ.0.AND.NNG(I).EQ.0) THEN
          NV=NV+1  
          C(I,1)=C(I,1)+1.D10  
          ENDIF
35        CONTINUE 
          IF(ILI.EQ.1) THEN
          WRITE(6,1020) M,N,NB,N**2,N*NB,N**2/(N*NB)
          ELSE
          WRITE(6,1021) M,N,NB,N**2,N*NB,N**2/(N*NB)
          ENDIF
C
          F1=FMATB(IMAX,C,RES,NMAX,N,NB)
          C(N+1,1)=DBLE(ILI)
          CALL TCHOLB(C,NMAX,N,NB) 
          CALL SOLSB(C,D,NMAX,N,NB)
          CALL SZEROR(A,IA,KA,B,P,D,SZ,SZ2,MMAX,M,N,NV)
          CALL INVERB(C,W,NMAX,N,NB)
          F2=FMATB(IMAX,C,RES,NMAX,N,NB)
          RAP=F1/F2
          IF(ILI.EQ.1) THEN
          WRITE(6,1120) RAP,IMAX
          ELSE
          WRITE(6,1121) RAP,IMAX
          ENDIF
          CALL DQVVB(A,IA,KA,C,P,MMAX,NMAX,M)
C            
          DO 40 I=1,N  
          W(I)=0.D0
          IF(SZ2*C(I,1).GT.1.D-10) W(I)=DSQRT(SZ2*C(I,1))  
40        CONTINUE
          IF(ILI.EQ.1) THEN 
          WRITE(6,1030)
          ELSE
          WRITE(6,1031)
          ENDIF
          DO 60 I=1,N  
          WRITE(6,1040) I,D(I),W(I)
60        CONTINUE 
          IF(IFILE.EQ.1) THEN  
          WRITE(3,1015) T0,STP,EX,N  
          WRITE(3,1018) (D(I),I=1,N)
          ENDIF
          REWIND 1
          IF(ILI.EQ.1) THEN 
          WRITE(6,1050)
          ELSE
          WRITE(6,1051)
          ENDIF
          I=0  
70        READ(1,1010) (DATO(L),L=1,120)
          READ(CDATO,FORMA) NPT,TS,XS
          IF(NPT.EQ.0) GOTO 80 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          I=I+1
          XX=XS
          TS=TS*SCT
          XX=XX*SCX
          SGNL=XX+B(I) 
          SRV=0.D0 
          IF(SZ2*(1.D0-P(I)).GT.1.D-10) SRV=DSQRT(SZ2*(1.D0-P(I))) 
          SQMS=SRV 
          SRV=0.D0 
          IF(SZ2*P(I).GT.1.D-10) SRV=DSQRT(SZ2*P(I))
          P(I)=SRV 
          ASTI=' ' 
          IF(DABS(B(I))/P(I).GE.TOLL) ASTI='*' 
          WRITE(6,1060) NPT,XS,SGNL,SQMS,B(I),P(I),ASTI 
          ENDIF
          RES(I)=B(I)
          GOTO 70  
80        NGL=M-N+NV
          IF(ILI.EQ.1) THEN
          WRITE(6,1071) SZ,M,N,NV,NGL  
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,1100) ALFA
          ELSE
          WRITE(6,1071) SZ,M,N,NV,NGL  
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,1101) ALFA
          ENDIF
          WRITE(6,1005)
          RETURN
1000      FORMAT(20A4) 
1005      FORMAT(///5X,20A4)
1010      FORMAT(120A1)
1012      FORMAT(5X,120A1) 
1015      FORMAT(3D20.12,I5)
1018      FORMAT(6D20.12)  
1020      FORMAT(///5X,'OCCUPAZIONE DI MEMORIA'//  
     .5X,'NUMERO DATI VALIDI',5X,I5/
     .5X,'NUMERO NODI       ',5X,I5/
     .5X,'AMPIEZZA DI BANDA ',5X,I5//  
     .39X,'PIENA',11X,'COMPATTA',7X,'FATTORE  DI'/70X,'RIEMPIMENTO'/
     .5X,'MATRICE NORMALE      ',8X,I10,9X,I10,10X,'1:',1X,I5) 
1021     FORMAT(///5X,'MEMORY STORAGE'//  
     .5X,'NUMBER OF DATA ',5X,I5/
     .5X,'NUMBER OF KNOTS',5X,I5/
     .5X,'BANDWIDTH      ',5X,I5//  
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
1100      FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEL TEST '
     .'DI REIEZIONE DEGLI ERRORI GROSSOLANI: ALFA = ',F5.3, 
     .' (SU DUE CODE)')
1101     FORMAT(/5X,'SIGNIFICANCE LEVEL OF THE TEST '
     .'FOR OUTLIER REJECTION: ALFA = ',F5.3, 
     .' (ON TWO SIDES)')
1120      FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
1121      FORMAT(///5X,'CONDITION NUMBER', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
1130      FORMAT(///5X,'NUMERO CLASSI             ',10X,I5,
     .//5X,'NUMERO NODI               ',10X,I5,
     .//5X,'PASSO DELLA INTERPOLAZIONE',5X,F10.3, 
     .//5X,'TEMPO ORIGINE             ',F15.3//) 
[B1131      FORMAT(///5X,'NUMBER OF CLASSES',10X,I5,
     .//5X,'NUMBER OF KNOTS  ',10X,I5,
     .//5X,'STEP   ',15X,F10.3, 
     .//5X,'INITIAL TIME     ',F15.3//) 
1140      FORMAT(///5X,'CLASSE -',10X,'NODI COINVOLTI DALLA CLASSE'/
     .5X,'NUMEROSITA',8X,'(ESTREMI DELLO INTERVALLO)'/)
1141      FORMAT(///7X,'CLASS -',11X,'KNOTS CONNECTED THE CLASS'/
     .7X,'QUANTITY',10X,'(BOUNDS OF THE INTERVAL)'/)
1150      FORMAT(5X,2I5,25X,2I5)
1190      FORMAT(//5X,'NUMERO DEI DATI',
     ./5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'//)
1191      FORMAT(//5X,'NUMBER OF DATA',
     ./5X,'EXCEEDING THE PREVIOUS DIMENSION'//)
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
