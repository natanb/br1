	  PROGRAM SPLINEP4
          PARAMETER (MMAX=5000,NMAX=3000,NBMAX=16,MNBMAX=MMAX*NBMAX, 
     .NNBMAX=25*NMAX)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,RESI,RES2,RES3,PCHAR,BLANK,PTO 
          CHARACTER*4 FORMA,FORMA2,FORMA3,TITOLO
          CHARACTER*15 CRESI,CRES2,CRES3
          DIMENSION A(MMAX,NBMAX),IA(MMAX,NBMAX),KA(MMAX),B(MMAX), 
     .P(MMAX),JA(MNBMAX),LA(NMAX),IB(NMAX),JB(NMAX),
     .C(NNBMAX),T(NNBMAX),JC(NNBMAX),LC(NMAX),D(MMAX), 
     .H(NMAX),R(NMAX),S(NMAX),ISRV(MMAX),NNG(NMAX),
     .U(MMAX),V(MMAX),W(MMAX),
     .IPX(MMAX),IPY(MMAX),LPX(MMAX+1),LPY(MMAX+1),IWD(NMAX),
     .DATO(120),FORMA(20),FORMA2(20),FORMA3(20),TITOLO(20),
     .PCHAR(120),IPTO(6),IAST(3),RESI(15),RES2(15),RES3(15)
          EQUIVALENCE (CRESI,RESI(1)),(CRES2,RES2(1)),(CRES3,RES3(1))
          DATA BLANK/' '/,PTO/'.'/ 
          DATA IPTO/6*120/,IAST/116,118,120/
C
          CALL LINGUA(ILI)
          ICASO=2 
          IF(ILI.EQ.1) THEN 
          CALL PASPLI(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IQUO,IAUTO,
     .ITABLE,TOLL,STP,STP2,STP3,SCX,SCY,SCZ,SCU,SCV,SCW,  
     .FORMA,FORMA2,FORMA3)
          ELSE  
          CALL PISPLI(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IQUO,IAUTO,
     .ITABLE,TOLL,STP,STP2,STP3,SCX,SCY,SCZ,SCU,SCV,SCW,  
     .FORMA,FORMA2,FORMA3)
          ENDIF 
          WRITE(6,1020) TITOLO 
          DO 5 I=1,120 
          PCHAR(I)=BLANK
5         CONTINUE 
          J=0  
10        READ(1,1000,END=30) (DATO(I),I=1,120)
          J=J+1
          DO 20 I=1,120
          IF(DATO(I).EQ.PTO) PCHAR(I)=PTO  
20        CONTINUE 
          GOTO 10 
30        N=J-1
          IF(N.GT.MMAX) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,1050)
          ELSE
          WRITE(6,1051) 
          ENDIF
          STOP 
          ENDIF
          J=0  
          K=0  
          DO 40 I=1,120
          IF(PCHAR(I).EQ.PTO) THEN 
          J=J+1
          IPTO(J)=I
          ENDIF
40        CONTINUE 
          NPTO=J
          IF(NPTO.NE.(ICASO+ITIPO+IQUO)) THEN
          IF(ILI.EQ.1) THEN  
          WRITE(6,1040)
          ELSE
          WRITE(6,1041)
          ENDIF
          STOP 
          ENDIF
C
          IF(ISCRI.EQ.1) THEN
          REWIND 1
          CALL TESTA(ICASO,ITIPO,IQUO,IPTO,NPTO,ILI)
          DO 50 J=1,N
          READ(1,1000) (DATO(I),I=1,120)
          WRITE(6,1030) (DATO(I),I=1,120)
50        CONTINUE
          ENDIF
          REWIND 1
          READ(1,1000) (DATO(I),I=1,120)
          LAST=IPTO(2+IQUO)
          DO 60 I=LAST+1,LAST+10
          IF(DATO(I).EQ.blank) GOTO 70
60        CONTINUE
70        LAST=I-1
C
          IF(IQUO.EQ.1) THEN
          DO 80 I=4,6  
          IPTO(I-1)=IPTO(I)
80        CONTINUE 
          ENDIF
C
                    K=IPTO(3)
                    L=IAST(1)
                    CALL SPLINP4(U,V,W,A,IA,P,B,KA,JA,LA,IB,JB, 
     .              C,T,JC,LC,D,H,R,S,ISRV,NNG,U,FORMA,STP,
     .              MMAX,NMAX,SCX,SCY,SCU,K,L,TOLL,
     .              IPX,IPY,LPX,LPY,IWD,
     .              IAUTO,ITABLE,IFILE,ILI) 
          IF(ITIPO.GE.2) THEN  
                    K=IPTO(4)
                    L=IAST(2)
                    CALL SPLINP4(U,V,W,A,IA,P,B,KA,JA,LA,IB,JB, 
     .              C,T,JC,LC,D,H,R,S,ISRV,NNG,V,FORMA2,STP2,  
     .              MMAX,NMAX,SCX,SCY,SCV,K,L,TOLL,
     .              IPX,IPY,LPX,LPY,IWD,
     .              IAUTO,ITABLE,IFILE,ILI) 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
                    K=IPTO(5)
                    L=IAST(3)
                    CALL SPLINP4(U,V,W,A,IA,P,B,KA,JA,LA,IB,JB, 
     .              C,T,JC,LC,D,H,R,S,ISRV,NNG,W,FORMA3,STP3,  
     .              MMAX,NMAX,SCX,SCY,SCW,K,L,TOLL,
     .              IPX,IPY,LPX,LPY,IWD,
     .              IAUTO,ITABLE,IFILE,ILI) 
          ENDIF
C
          IF(IFILE.EQ.1) THEN
          KONTU=0
          KONTV=0
          KONTW=0
	  REWIND 1
          DO 100 J=1,N
          READ(1,1000) (DATO(I),I=1,120)
          DO 90 I=1,15
          RESI(I)=BLANK
          RES2(I)=BLANK
          RES3(I)=BLANK
90        CONTINUE
          IF(DATO(IPTO(3)).EQ.PTO) THEN
          KONTU=KONTU+1
          WRITE(CRESI,1010) U(KONTU)
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO) THEN
          KONTV=KONTV+1
          WRITE(CRES2,1010) V(KONTV)
          ENDIF
          IF(DATO(IPTO(5)).EQ.PTO) THEN
          KONTW=KONTW+1
          WRITE(CRES3,1010) W(KONTW)
          ENDIF
          WRITE(2,1000) (DATO(I),I=1,LAST),RESI,RES2,RES3
100       CONTINUE
          WRITE(2,1000)
          ENDIF
          STOP 
1000      FORMAT(120A1)
1010      FORMAT(F15.5)
1020      FORMAT(/////5X,20A4//)
1030      FORMAT(5X,120A1) 
1040      FORMAT(5X,'CASO E/O TIPO DEI DATI ERRATO'///)
1041      FORMAT(5X,'WRONG CASE AND/OR TYPE OF DATA'///)
1050      FORMAT(5X,'NUMERO DEI DATI ',/5X,
     .'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'///)  
1051      FORMAT(5X,'NUMBER OF DATA ',/5X,
     .'EXCEEDING THE PREVIOUS DIMENSION'///)  
          END  
C*****************************     SPLINP4
          SUBROUTINE SPLINP4(X,Y,Z,A,IA,P,B,KA,JA,LA,IB,JB,
     .C,T,JC,LC,D,W,R,U,ISRV,NNG,RES,FORMA,STP,
     .MMAX,NMAX,SCX,SCY,SCU,KK,LL,TOLL,
     .IPX,IPY,LPX,LPY,IWD,
     .IAUTO,ITABLE,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,AST,PTO,ASTI
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),Z(1),A(MMAX,1),IA(MMAX,1),KA(1),B(1),
     .P(1),JA(1),LA(1),IB(1),JB(1),
     .C(1),T(1),JC(1),LC(1),D(1),
     .IPX(1),IPY(1),LPX(1),LPY(1),IWD(1),
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
          NCX=IDINT((X9-X0)*SCX/STP)+3
          NCY=IDINT((Y9-Y0)*SCY/STP)+3
          X0=X0-STP/SCX
          Y0=Y0-STP/SCY
          NTC=NCX*NCY
          NIX=NCX+1
          NIY=NCY+1
          NTI=NIX*NIY
          IF(ILI.EQ.1) THEN
          WRITE(6,1171) M,NCX,NCY,NTC,NTI,STP,X0,Y0
          ELSE
          WRITE(6,1170) M,NCX,NCY,NTC,NTI,STP,X0,Y0
          ENDIF

	  MM=M
	  DO I=1,M
	  US=Z(I)
          P(I)=1.D0
          B(I)=-US*SCU
          KA(I)=16
          XX=(X(I)-X0)*SCX
          YY=(Y(I)-Y0)*SCY
          NX=IDINT(XX/STP)+1
          NY=IDINT(YY/STP)+1
          NC=(NY-1)*NCX+NX
          NNG(NC)=NNG(NC)+1
	  ISRV(I)=NC
          SRVX=(XX-(NX-1)*STP)/STP
          SRVY=(YY-(NY-1)*STP)/STP
          NI=(NY-2)*NIX+NX-2
          DO 50 J=1,4
          AY=DBLE(2-J)
          JY=5-J
          VALY=SPLT(SRVY+AY,JY)
          DO 50 K=1,4  
          JS=(J-1)*4+K 
          AX=DBLE(2-K) 
          JX=5-K
          A(I,JS)=VALY*SPLT(SRVX+AX,JX)
          IA(I,JS)=NI+(J-1)*NIX+K
50        CONTINUE 
          ENDDO 
          N=NTI

	  DO I =2,NCY-1
	   DO J =2,NCX-1
	    XC=X0+STP/2+(J-1)*STP
	    YC=Y0+STP/2+(I-1)*STP
  	    NC=(I-1)*NCX+J
	    IF(NNG(NC).EQ.0) WRITE(11,*) NC,NNG(NC),XC,YC
	   ENDDO
          ENDDO



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
	  P(M)=100.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
          b(m)=-z(ipnear)*scu
c	  write(*,'(2i5,3f15.6)')i,ipnear,z(ipnear),x(ipnear),y(ipnear)
80        CONTINUE 
          DO 85 I=NIX+2,2*NIX-1
          NV=NV+1
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=100.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
          b(m)=-z(ipnear)*scu
c	  write(*,'(2i5,3f15.6)')i,ipnear,z(ipnear),x(ipnear),y(ipnear)
85        CONTINUE 

	  DO 90 I=NIX+1,N-2*NIX+1,NIX
          NV=NV+1
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=100.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
	b(m)=-z(ipnear)*scu
c	  write(*,'(2i5,3f15.6)')i,ipnear,z(ipnear),x(ipnear),y(ipnear)
90        CONTINUE 

	  DO 95 I=NIX+2,N-2*NIX+2,NIX
          NV=NV+1
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=100.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
	b(m)=-z(ipnear)*scu
c	  write(*,'(2i5,3f15.6)')i,ipnear,z(ipnear),x(ipnear),y(ipnear)
95        CONTINUE 

          DO 100 I=2*NIX,N-NIX,NIX 
          NV=NV+1  
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=100.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
	b(m)=-z(ipnear)*scu
c	  write(*,'(2i5,3f15.6)')i,ipnear,z(ipnear),x(ipnear),y(ipnear)
100       CONTINUE 
          DO 105 I=2*NIX-1,N-NIX-1,NIX 
          NV=NV+1  
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=100.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
	b(m)=-z(ipnear)*scu
c	  write(*,'(2i5,3f15.6)')i,ipnear,z(ipnear),x(ipnear),y(ipnear)
105       CONTINUE 
	
          DO 110 I=(NIY-1)*NIX+1,N 
          NV=NV+1  
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=100.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
c	  write(*,'(2i5,3f15.6)')i,ipnear,z(ipnear),x(ipnear),y(ipnear)
	b(m)=-z(ipnear)*scu
110       CONTINUE 
          DO 115 I=(NIY-2)*NIX+2,N-NIX-1 
          NV=NV+1  
	  CALL VINBORDI(LA,JA,IB,NP,I,NCX,NCY)
	  IPNEAR=NEARPOIN(I,NCX,X,Y,X0,Y0,STP,IB,NP)
	  M=M+1
	  A(M,1)=1.D0
	  IA(M,1)=I
	  KA(M)=1
	  P(M)=100.D0
c          B(M)=-(Z(IPNEAR)-EU)*SCU  
c	  write(*,'(2i5,3f15.6)')i,ipnear,z(ipnear),x(ipnear),y(ipnear)
	b(m)=-z(ipnear)*scu
115       CONTINUE 


          NMIN=3  ! n. min di pti circostanti
          CALL BIORDA(X,Y,IPX,IPY,LPX,LPY,MM)
          
          DO 720 II=2,NCY-1
          DO 720 JJ=2,NCX-1
  	  NC=(II-1)*NCX+JJ
          IF(NNG(NC).EQ.0) THEN 
	  XC=X0+STP/2+(JJ-1)*STP
	  YC=Y0+STP/2+(II-1)*STP
          AMPX=1.5*STP
          AMPY=1.5*STP
          NWD=0
          DO 500 KL=1,MAX(NCX,NCY)
          CALL WINDOW(XC,YC,AMPX,AMPY,X,Y,IPX,IPY,LPX,LPY,IWD,
     .MMAX,NMAX,MM,NWD,ILI)
          IF(NWD.GT.NMIN) GOTO 505
          AMPX=AMPX+STP
          AMPY=AMPY+STP
500	  CONTINUE
       
505	  DISTMIN=1.D30
          PRINT *,'NWD = ',NWD,' NC = ',NC,' NNG(NC) = ',NNG(NC)
          IPT=0
 	  DISTTOT=0.
          ZMED=0.
          DO KL=1,NWD
          IS=IWD(KL)
          DIST=(XC-X(IS))**2+(YC-Y(IS))**2
          ZMED=ZMED+Z(IS)/DIST
          DISTTOT=DISTTOT+1./DIST
          IF(DIST.LT.DISTMIN) THEN
            DISTMIN=DIST
            IPT=IS
          ENDIF
          ENDDO
          ZMED=ZMED/DISTTOT
C         ZMED=Z(IPT)
          print *,' nx ny nc ',nx,ny,NC,' XC  YC ',XC,YC
          DO KKL=-1,1,2

C          PRINT *,'IPT =',IPT,X(IPT),Y(IPT),ZMED
          XX=(XC+KKL*0.25*STP-X0)*SCX
          NX=IDINT(XX/STP)+1
          SRVX=(XX-(NX-1)*STP)/STP

          DO KLL=-1,1,2

          NC=(NY-1)*NCX+NX

          NY=IDINT(YY/STP)+1
          YY=(YC+KLL*0.25*STP-Y0)*SCY
          SRVY=(YY-(NY-1)*STP)/STP
          NI=(NY-2)*NIX+NX-2       
          write(10,'(i5,5x,2(f15.4,10x),f15.5)') IPT,XX,YY,ZMED

          M=M+1
          P(M)=1.D0
          B(M)=-ZMED*SCU
          KA(M)=16

          DO 610 J=1,4
          AY=DBLE(2-J)
          JY=5-J
          VALY=SPLT(SRVY+AY,JY)
          DO 610 K=1,4  
          JS=(J-1)*4+K 
          AX=DBLE(2-K) 
          JX=5-K
          A(M,JS)=VALY*SPLT(SRVX+AX,JX)
          IA(M,JS)=NI+(J-1)*NIX+K
610       CONTINUE 
          ENDDO  ! KLL
          ENDDO  ! KKL
          ENDIF
720	  CONTINUE	  


          PRINT *,'MM = ', MM, ' M= ',M, ' N= ',N
       
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
c          WRITE(6,1020) M,N,NN/N,N*(N+1)/2,NN,N*(N+1)/(2*NN)
          WRITE(6,1020) M,N,NN/N,N*((N+1)/2),NN,N*(N+1)/(2*NN)
          ELSE
c          WRITE(6,1021) M,N,NN/N,N*(N+1)/2,NN,N*(N+1)/(2*NN)
          WRITE(6,1021) M,N,NN/N,N*((N+1)/2),NN,N*(N+1)/(2*NN)
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
          NV=0
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
          READ(CDATO,FORMA) NPT,XX,YY,US
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
          IF(ASTI.EQ.'*') THEN
          WRITE(2,'(I5,5x,2(f15.4,10x),2f15.5,a)')NPT,XX,YY,US,B(I),ASTI
          ELSE
          WRITE(2,'(I5,5x,2(f15.4,10x),2f15.5,a)')NPT,XX,YY,US
          ENDIF
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
1015      FORMAT(4D20.12,I10,I5)  
1018      FORMAT(6D20.12)  
1020      FORMAT(///5X,'OCCUPAZIONE DI MEMORIA'//  
     .5X,'NUMERO DATI VALIDI',5X,I7/
     .5X,'NUMERO NODI       ',5X,I7/
     .5X,'DENSITA MEDIA     ',5X,I7//  
     .39X,'PIENA',11X,'COMPATTA',7X,'FATTORE  DI'/70X,'RIEMPIMENTO'/
     .5X,'MATRICE NORMALE      ',8X,I10,9X,I10,10X,'1:',1X,I7) 
1021     FORMAT(///5X,'MEMORY STORAGE'//  
     .5X,'NUMBER OF DATA ',5X,I7/
     .5X,'NUMBER OF KNOTS',5X,I7/
     .5X,'AVERAGE DENSITY',5X,I7//  
     .36X,'FULL',14X,'SPARSE',11X,'FILL-IN'/
     .5X,'NORMAL MATRIX  ',10X,I10,10X,I10,10X,'1:',1X,I7) 
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
     ./5X,'NUMERO EQUAZIONI  ',I10,/5X,'NUMERO INCOGNITE  ',I10, 
     ./5X,'NUMERO VINCOLI    ',I10,/5X,'GRADI DI LIBERTA  ',I10) 
1071      FORMAT(////5X,'SIGMA ZERO   ',F10.3,//
     ./5X,'NUMBER OF EQUATIONS  ',I10,/5X,'NUMBER OF UNKNOWNS   ',I10, 
     ./5X,'NUMBER OF CONSTRAINTS',I10,/5X,'DEGREES OF FREEDOM   ',I10) 
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
1170      FORMAT(/
     .//5x,'NUMERO DI DATI VALIDI     ',10X,I10,
     .//5X,'NUMERO CLASSI IN X E Y    ',2I10,
     .//5X,'NUMERO TOTALE CLASSI      ',10X,I10,
     .//5X,'NUMERO NODI               ',10X,I10,
     .//5X,'PASSO DELLA INTERPOLAZIONE',10X,F10.3, 
     .//5X,'COORDINATE ORIGINE ',F15.3,5X,F15.3//) 
1171      FORMAT(/
     .//5X,'NUMBER OF VALID DATA        ',10X,I10,
     .//5X,'NUMBER OF CLASSES IN X AND Y',2I10,
     .//5X,'NUMBER OF CLASSES           ',10X,I10,
     .//5X,'NUMBER OF KNOTS             ',10X,I10,
     .//5X,'STEP                        ',10X,F10.3, 
     .//5X,'INITIAL COORDINATES ',F15.3,5X,F15.3//) 
1190      FORMAT(//5X,'NUMERO DEI DATI',
     ./5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'//)
1191      FORMAT(//5X,'NUMBER OF DATA',
     ./5X,'EXCEEDING THE PREVIOUS DIMENSION'//)
          END  
