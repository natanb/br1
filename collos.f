          PARAMETER (NMAX=4600,KMAX=1,MMAX=KMAX*NMAX,
     .NWCMAX=4600,MWCMAX=KMAX*NWCMAX,NWDMAX=4600,
     .NMAX1=NMAX+1,MMAX1=MMAX+1,NNWMAX=NMAX*NWCMAX,
     .MMWMAX=MMAX*MWCMAX,MMW2MX=MMWMAX/2,LLWMAX=NNWMAX)
c    ,LLWMAX=MAX0(NNWMAX,MMW2MX)) error in compilation
          IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
          CHARACTER*1 DATO,PCHAR,PTO,AST,BLANK 
          CHARACTER*4 FORMA,TITOLO 
          CHARACTER*120 CDATO
          DIMENSION V(MMAX),H(MMAX),LC(MMAX1),JC(LLWMAX),C(MMW2MX),
     .T(MMW2MX),ICOVA(MWCMAX,3),COVA(MWCMAX,3),KONT(NMAX1),
     .R(MMAX),U(MMAX),W(MMAX), 
     .X(NMAX),Y(NMAX),Z(NMAX),IPX(NMAX),IPY(NMAX),IPZ(NMAX),
     .LPX(NMAX),LPY(NMAX),LPZ(NMAX), 
     .JCRONA(MMAX1),JCOVAS(LLWMAX),IBANDA(NMAX),JBANDA(NMAX),  
     .NUMERO(NMAX),INDIR(NMAX),
     .NOME(NMAX),NTIPO(NMAX),IWD(NWDMAX),IPUNT(NWDMAX),
     .DA(6),DB(6),DC(6),SNOISE(6),NCASO(6),DATO(120),PCHAR(120),IPTO(6),
     .IAST(3),IND(7),ICA(3,3,7),ISGN(2,6),AA(3,3),BB(3,3),CC(3,3),D(3), 
     .E(3),KP(3),LP(3),S1(3),S2(3),S3(3),S4(3),S5(3),S6(3),S7(3),S8(3), 
     .TGL(6),FORMA(20),TITOLO(20),ND(7,3)  
          EQUIVALENCE (CDATO,DATO(1)) 
          DATA IND/1,1,2,1,2,2,3/  
          DATA PTO/'.'/,AST/'*'/,BLANK/' '/
          DATA IPTO/6*120/,IAST/116,118,120/
          DATA ((ND(I,J),I=1,7),J=1,3)/1,2,1,3,1,2,1,0,0,2,0,3,3,2,
     .0,0,0,0,0,0,3/
          DATA ((ICA(I,J,1),J=1,3),I=1,3)/1,0,0,0,0,0,0,0,0/
          DATA ((ICA(I,J,2),J=1,3),I=1,3)/2,0,0,0,0,0,0,0,0/
          DATA ((ICA(I,J,3),J=1,3),I=1,3)/1,3,0,0,2,0,0,0,0/
          DATA ((ICA(I,J,4),J=1,3),I=1,3)/4,0,0,0,0,0,0,0,0/
          DATA ((ICA(I,J,5),J=1,3),I=1,3)/1,5,0,0,4,0,0,0,0/
          DATA ((ICA(I,J,6),J=1,3),I=1,3)/2,6,0,0,4,0,0,0,0/
          DATA ((ICA(I,J,7),J=1,3),I=1,3)/1,3,5,0,2,6,0,0,4/
          DATA ((ISGN(I,J),J=1,6),I=1,2)/1,1,-1,1,-1,-1,1,1,1,1,1,1/
C  
          READ(5,*) ILI
          READ(5,1010) TITOLO  
          WRITE(6,1020) TITOLO 
          READ(5,*) ICASO,ITIPO
          READ(5,*) ISCRI,ITABLE,IFILE,IMATRI,IORDI,TOLL  
          READ(5,*) SCX,SCY,SCZ,SCU,SCV,SCW
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
          READ(2,1050) NCASO(I),DA(I),DB(I),DC(I),SNOISE(I),TGL(I)
          WRITE(6,1075) NCASO(I),DA(I),DB(I),DC(I),SNOISE(I),TGL(I)
1         CONTINUE 
          XCOVT=0.D0
          DO 3 I=1,NTER
          XCOVT=DMAX1(XCOVT,TGL(I))
3         CONTINUE 
          IDIAG=1
          IF(DA(3)+DA(5)+DA(6).NE.0.D0) IDIAG=2
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
          WRITE(6,1190)
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
          READ(1,1000) (DATO(I),I=1,120)
          LAST=IPTO(3)
          DO 47 I=LAST+1,LAST+10
          IF(DATO(I).EQ.BLANK) GOTO 48
47        CONTINUE
48        LAST=I-1        
C
          REWIND 1 
          K=0  
          J=0  
          EX=0.D0  
          EY=0.D0 
          EZ=0.D0 
50        READ(1,1000)(DATO(I),I=1,120)
          IF(ITIPO.EQ.1) READ(CDATO,FORMA) NPT,XS,YS,ZS,US
          IF(ITIPO.EQ.2) READ(CDATO,FORMA) NPT,XS,YS,ZS,US,VS
          IF(ITIPO.EQ.3) READ(CDATO,FORMA) NPT,XS,YS,ZS,US,VS,WS 
          IF(NPT.EQ.0) GOTO 55
          IP1=0
          IP2=0
          IP3=0
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          K=K+1
          V(K)=US*SCU
          IP1=1
          ENDIF
          IF(DATO(IPTO(5)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          K=K+1
          V(K)=VS*SCV
          IP2=2
          ENDIF
          IF(DATO(IPTO(6)).EQ.PTO.AND.DATO(IAST(3)).NE.AST) THEN
          K=K+1
          V(K)=WS*SCW
          IP3=4
          ENDIF
          IF((IP1+IP2+IP3).EQ.0) GOTO 50  
          J=J+1
          NOME(J)=NPT  
          EX=EX+XS*SCX 
          EY=EY+YS*SCY 
          EZ=EZ+ZS*SCZ
          X(J)=XS*SCX  
          Y(J)=YS*SCY 
          Z(J)=ZS*SCZ 
          NTIPO(J)=IP1+IP2+IP3 
          GOTO 50 
55        N=J  
          M=K
          IF(N.GT.NMAX.OR.M.GT.MMAX) THEN 
          IF(ILI.EQ.1) THEN 
          WRITE(6,1140)
          ELSE
          WRITE(6,1141)
          ENDIF
          STOP 
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1008)
          ELSE
          WRITE(6,1009)
          ENDIF
          EX=EX/N  
          EY=EY/N
          EZ=EZ/N
          DO 56 J=1,3
          E(J)=0.D0
          LP(J)=0
56        CONTINUE
          JJ=0
          DO 58 I=1,N
          KK=NTIPO(I)
          L=IND(KK)
          DO 57 J=1,L
          JJ=JJ+1
          IQ=ND(KK,J)
          LP(IQ)=LP(IQ)+1
          E(IQ)=E(IQ)+V(JJ)
57        CONTINUE
58        CONTINUE
          DO 59 I=1,ITIPO
          E(I)=E(I)/LP(I)
59        CONTINUE
          IF(ITIPO.EQ.1) WRITE(CDATO,FORMA) N,EX/SCX,EY/SCY,EZ/SCZ,
     .E(1)/SCU 
          IF(ITIPO.EQ.2) WRITE(CDATO,FORMA) N,EX/SCX,EY/SCY,EZ/SCZ,
     .E(1)/SCU,E(2)/SCV
          IF(ITIPO.EQ.3) WRITE(CDATO,FORMA) N,EX/SCX,EY/SCY,EZ/SCZ,
     .E(1)/SCU,E(2)/SCV,E(3)/SCW
          WRITE(6,1005) (DATO(I),I=1,120)
          JJ=0
          DO 60 I=1,N
          KK=NTIPO(I)
          L=IND(KK)
          DO 60 J=1,L
          JJ=JJ+1
          IQ=ND(KK,J)
          V(JJ)=V(JJ)-E(IQ)
60        CONTINUE  
C
          CALL SPARS3(XCOVT,X,Y,Z,IPX,IPY,IPZ,LPX,LPY,LPZ,
     .IWD,IPUNT,JCRONA,JCOVAS,
     .IBANDA,JBANDA,NUMERO,INDIR,LC,JC,NMAX,NWDMAX,N,NWD,IORDI,ILI)
          IF(ITABLE.EQ.1) THEN 
          IF(ILI.EQ.1) THEN
          WRITE(6,1150)
          ELSE
          WRITE(6,1151)
          ENDIF
          DO 61 I=1,N,20
          I1=I 
          I2=MIN0(N,I+19)  
          WRITE(6,1160) (J,J=I1,I2)
          WRITE(6,1160) (NUMERO(J),J=I1,I2)
          WRITE(6,1160) (NOME(J),J=I1,I2)  
          WRITE(6,1000)
61        CONTINUE 
          NOME(N+1)=ILI
          CALL WGRAFO(NOME,JCOVAS,JCRONA,IBANDA,JBANDA,N,1)
          ENDIF
          KONT(1)=0
          DO 62 I=1,N  
          J=INDIR(I)
          JJ=NTIPO(J)  
          II=IND(JJ)
          KONT(I+1)=KONT(I)+II 
62        CONTINUE 
          IF(ITABLE.EQ.1) THEN 
          IF(ILI.EQ.1) THEN
          WRITE(6,1170)
          ELSE
          WRITE(6,1171)
          ENDIF
          DO 64 I=1,N,20
          I1=I 
          I2=MIN0(N,I+19)  
          IS=I2-I+1
          DO 63 J=1,IS 
          IWD(J)=KONT(I+J)-KONT(I+J-1) 
63        CONTINUE 
          WRITE(6,1160) (J,J=I1,I2)
          WRITE(6,1160) (IWD(J),J=1,IS)
          WRITE(6,1160) (KONT(J+1),J=I1,I2)
          WRITE(6,1000)
64        CONTINUE 
          ENDIF
          NU=JCRONA(N+1)-1 
          IF(ITIPO**2*NU.GT.MMWMAX) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,1130)
          ELSE
          WRITE(6,1131)
          ENDIF
          STOP 
          ENDIF
          LC(1)=1  
          KK=0 
          INC=0
          INR=0
          IP=INDIR(1)  
          II=NTIPO(IP) 
          IR=IND(II)
          DO 70 L=1,IR 
          DO 65 K=1,L  
          IS=ICA(K,L,II)
          KK=KK+1  
          C(KK)=SNOISE(IS)+DA(IS)  
          JC(KK)=INR+K 
65        CONTINUE 
          LC(INC+L+1)=KK+1 
70        CONTINUE
          KK=KK+1 
          INC=INC+IR
          DO 110 I=2,N 
          INR=0
          IP=INDIR(I)  
          II=NTIPO(IP) 
          IC=IND(II)
          J1=JCRONA(I) 
          J2=JCRONA(I+1)-1 
          DO 75 JS=J1,J2
          IWD(JS+1-J1)=JCOVAS(JS)  
          IPUNT(JS+1-J1)=JS+1-J1
75        CONTINUE 
          NWD=J2+1-J1  
          CALL BINARY(I,IWD,IPUNT,NWDMAX,JM,JM,JM,JM,NWD)  
          IF(JM.GT.1) THEN 
          DO 90 JS=J1,J1+JM-2  
          J=JCOVAS(JS) 
          JP=INDIR(J)  
          JJ=NTIPO(JP) 
          IR=IND(JJ)
          DELTA=DSQRT((X(IP)-X(JP))**2+(Y(IP)-Y(JP))**2+
     .(Z(IP)-Z(JP))**2)     
          DO 80 L=1,IC 
          DO 80 K=1,IR 
          COVA(INR+K,L)=GASS3(DELTA,XCOVT,JJ,II,K,L,
     .DA,DB,DC,SNOISE,NCASO,ICASO)
          ICOVA(INR+K,L)=(KONT(J)+K)*ISGN(IDIAG,ICASO) 
80        CONTINUE 
          INR=INR+IR
90        CONTINUE 
          ENDIF
          DO 92 L=1,IC 
          DO 92 K=1,L  
          IS=ICA(K,L,II)
          COVA(INR+K,L)=SNOISE(IS)+DA(IS)  
          ICOVA(INR+K,L)=KONT(I)+K 
92        CONTINUE 
          INR=INR+IC         
          DO 100 L=1,IC
          DO 95 K=1,INR+L-IC 
          C(KK)=COVA(K,L)  
          JC(KK)=ICOVA(K,L)
          KK=KK+(ICOVA(K,L)/IABS(ICOVA(K,L))+1)/2
95        CONTINUE 
          LC(INC+L+1)=KK 
100       CONTINUE 
          INC=INC+IC
110       CONTINUE 
          IF(IORDI.EQ.1) THEN  
          K=0  
          DO 112 I=1,N 
          IPT=NUMERO(I)
          INR=KONT(IPT)
          II=NTIPO(I)  
          IR=IND(II)
          DO 112 J=1,IR
          K=K+1
          H(INR+J)=V(K)
112       CONTINUE 
          DO 115 I=1,M 
          V(I)=H(I)
115       CONTINUE 
          ELSE 
          DO 120 I=1,M 
          H(I)=V(I)
120       CONTINUE 
          ENDIF
          MU=LC(M+1)-1 
          IF(ILI.EQ.1) THEN
          WRITE(6,1110) N,M,MU/M,M*(M+1)/2,MU,M*(M+1)/(2*MU)
          ELSE
          WRITE(6,1111) N,M,MU/M,M*(M+1)/2,MU,M*(M+1)/(2*MU)
          ENDIF
          IF(IMATRI.EQ.1) THEN 
          WRITE(12,1095) LC(1) 
          DO 122 I=1,M 
          WRITE(12,1095) LC(I+1)
          L1=LC(I) 
          L2=LC(I+1)-1 
          WRITE(12,1095) (JC(J),J=L1,L2)
          WRITE(12,1090) (C(J),J=L1,L2),Z(I)
122       CONTINUE 
          ENDIF
C 
          F1=FMATI1(IMAX,C,JC,LC,W,M)  
          CALL ALGITE(C,T,JC,LC,H, 
     .JCOVAS,JCRONA,IBANDA,JBANDA, 
     .1.D-3,5,ITE,1,1,R,U,W,M)
          F2=FMATI2(IMAX,C,JCOVAS,JC,JCRONA,W,M)
          RAP=F1/F2
          IF(ILI.EQ.1) THEN
          WRITE(6,1180) ITE
          WRITE(6,1120) RAP,IMAX
          ELSE
          WRITE(6,1181) ITE
          WRITE(6,1121) RAP,IMAX
          ENDIF
          IF(IFILE.EQ.1) THEN  
          WRITE(10,*) M
          WRITE(10,1090) (-H(I),I=1,M) 
          WRITE(10,*) N
          WRITE(10,1093) (X(I),Y(I),Z(I),I=1,N) 
          WRITE(10,1095) (IPX(I),IPY(I),IPZ(I),
     .LPX(I),LPY(I),LPZ(I),NTIPO(I),INDIR(I),I=1,N)  
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
          DO 220 I=1,N 
          READ(1,1000) (DATO(J),J=1,120)
          NPT=NOME(I) 
          IPT=NUMERO(I)
          INR=KONT(IPT)
          II=NTIPO(I)  
          IR=IND(II)
          DO 130 L=1,3 
          DO 130 K=1,3 
          AA(L,K)=0.D0 
          BB(L,K)=0.D0 
          CC(L,K)=0.D0 
130       CONTINUE 
          DO 150 L=1,IR
          J=LC(INR+L+1)-1  
          DO 150 K=1,L 
          JJ=J+K-L 
          IS=ICA(K,L,II)
          AA(K,L)=SNOISE(IS)
          AA(L,K)=AA(K,L)  
          BB(K,L)=C(JJ)
          BB(L,K)=BB(K,L)  
150       CONTINUE
          DO 170 L=1,IR
          SUM=0.D0 
          DO 160 K=1,IR
          SUM=SUM-AA(L,K)*H(INR+K) 
160       CONTINUE 
          E(L)=SUM 
170       CONTINUE 
          DO 180 L=1,IR
          D(L)=V(INR+L)-E(L)
180       CONTINUE 
          CALL MAT3(AA,BB,CC)  
          DO 190 L=1,IR
          DO 190 K=L,IR
          CC(L,K)=AA(L,K)-BB(L,K)  
          IS=ICA(L,K,II)
          AA(L,K)=DA(IS)-CC(L,K)
190       CONTINUE 
          DO 210 L=1,IR
          SRV=AA(L,L)  
          AA(L,L)=0.D0 
          IF(SRV.GT.1.D-10) AA(L,L)=DSQRT(SRV) 
          SRV=BB(L,L)  
          BB(L,L)=0.D0 
          IF(SRV.GT.1.D-10) BB(L,L)=DSQRT(SRV) 
          U(INR+L)=BB(L,L) 
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
          SRV=CC(L,L)  
          CC(L,L)=0.D0 
          IF(SRV.GT.1.D-10) CC(L,L)=DSQRT(SRV) 
          IF(L.GT.1) THEN  
          DO 200 K=1,L-1
          SRV=0.D0
          PRD=AA(L,L)*AA(K,K)
          IF(DABS(PRD).GT.1.D-10) SRV=AA(K,L)/PRD
          AA(K,L)=SRV
          SRV=0.D0
          PRD=BB(L,L)*BB(K,K)
          IF(DABS(PRD).GT.1.D-10) SRV=BB(K,L)/PRD
          BB(K,L)=SRV
          SRV=0.D0
          PRD=CC(L,L)*CC(K,K)
          IF(DABS(PRD).GT.1.D-10) SRV=CC(K,L)/PRD
          CC(K,L)=SRV
200       CONTINUE 
          ENDIF
210       CONTINUE 
          CALL CODA(I,NPT,II,D,AA,E,BB,CC,IFILE,DATO,LAST,TOLL)  
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
          SRV=(S2(I)/K-S1(I)**2)*K/(K-1)
          S2(I)=0.D0
          IF(SRV.GT.1.D-10) S2(I)=DSQRT(SRV)
          S3(I)=S3(I)/K
          SRV=(S4(I)/K-S3(I)**2)*K/(K-1)
          S4(I)=0.D0
          IF(SRV.GT.1.D-10) S4(I)=DSQRT(SRV)
          S5(I)=S5(I)/LP(I)
          SRV=(S6(I)/LP(I)-S5(I)**2)*LP(I)/(LP(I)-1)
          S6(I)=0.D0
          IF(SRV.GT.1.D-10) S6(I)=DSQRT(SRV)
          SRV=S7(I)
          S7(I)=0.D0
          IF(SRV.GT.1.D-10) S7(I)=DSQRT(SRV/K) 
          SRV=S8(I)
          S8(I)=0.D0
          IF(SRV.GT.1.D-10) S8(I)=DSQRT(SRV)
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
1008      FORMAT(///5X,'NUMERO DEI DATI - MEDIE DELLE COORDINATE ',
     .'E DELLE COMPONENTI'/)
1009      FORMAT(///5X,'NUMBER OF DATA - AVERAGE OF THE COORDINATES ',
     .'AND THE COMPONENTS'/)
1010      FORMAT(20A4) 
1020      FORMAT(///5X,20A4)
1030      FORMAT(2I5,F10.3)
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
1093      FORMAT(3D20.12)  
1095      FORMAT(8I10) 
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
     .5X,'DENSITA MEDIA     ',5X,I5//  
     .40X,'PIENA',12X,'COMPATTA',9X,'FATTORE  DI'/74X,'RIEMPIMENTO'/
     .5X,'MATRICE DI COVARIANZA',9X,I10,10X,I10,10X,'1:',1X,I5) 
1111     FORMAT(///5X,'MEMORY STORAGE'//  
     .5X,'NUMBER OF DATA       ',5X,I5/
     .5X,'NUMBER OF COMPONENTS ',5X,I5//  
     .5X,'AVERAGE DENSITY      ',5X,I5//  
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
1150      FORMAT(///5X,'RIORDINO DELLA NUMERAZIONE'//  
     .5X,'(CONTATORE, NUMERO DI ORDINE E NOME DEI PUNTI) '/)
1151      FORMAT(///5X,'RIORDERING OF THE NUMERATION'//  
     .5X,'(COUNTER, ORDERING NUMBER AND NAME OF THE POINTS) '/)
1160      FORMAT(5X,20I5)  
1170      FORMAT(///5X,'CONTATORE DELLE INCOGNITE'//
     .5X,'(NUMERO DI ORDINE DEI PUNTI, N. DI INCOGNITE PER ',
     .'PUNTO, SOMMATORIA PROGRESSIVA DEL N. DI INCOGNITE)'/)  
1171      FORMAT(///5X,'COUNTER OF THE UNKNOWNS'//
     .5X,'(ORDERING NUMBER OF THE POINTS, NUMBER OF UNKNOWNS FOR ',
     .'EACH POINT, POINTER ARRAY)'/)  
1180      FORMAT(///5X,'NUMERO DI ITERAZIONI (METODO DEL GRADIENTE ',  
     .'CONIUGATO)',I5) 
1181      FORMAT(///5X,'ITERATION NUMBER (CONJUGATE GRADIENT METHOD)',
     .I5) 
1190      FORMAT(5X,'CASO E/O TIPO DEI DATI ERRATO'///)
1191      FORMAT(5X,'WRONG CASE AND/OR TYPE OF DATA'///)        
          END  

