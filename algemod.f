C*************************         ABINRY          *************************
          SUBROUTINE ABINRY(VALORE,HNFORM,IPUNT,NMAX,LNINF,LNSUP,
     .    LPINF,LPSUP,N)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION HNFORM(1),IPUNT(1)
          LPINF=0
          LPSUP=0
          LNINF=0
          LNSUP=0
          NI=1
          NF=N
          IPUNT(N+1)=NMAX
          HNFORM(NMAX)=1.D10
          LI=IPUNT(NI)
          LF=IPUNT(NF)
          IF(VALORE.LT.HNFORM(LI)) THEN
          LPINF=1
          LPSUP=1
          LNINF=IPUNT(1)
          LNSUP=IPUNT(1)
          IPUNT(N+1)=0
          HNFORM(NMAX)=0.D0
          RETURN
          ENDIF
          IF(VALORE.GT.HNFORM(LF)) THEN
          LPINF=N
          LPSUP=N
          LNINF=IPUNT(N)
          LNSUP=IPUNT(N)
          IPUNT(N+1)=0
          HNFORM(NMAX)=0.D0
          RETURN
          ENDIF
          NF=NF+1
5         NC=(NI+NF)/2
          LC=IPUNT(NC)
          IF(VALORE.EQ.HNFORM(LC)) THEN
          DO 10 I=NC,NF
          L=IPUNT(I)
          IF(VALORE.NE.HNFORM(L)) GO TO 12
          LPSUP=I
          LNSUP=IPUNT(I)
10        CONTINUE
12        DO 15 I=NC,NI,-1
          L=IPUNT(I)
          IF(VALORE.NE.HNFORM(L)) GO TO 18
          LPINF=I
          LNINF=IPUNT(I)
15        CONTINUE
18        IPUNT(N+1)=0
          HNFORM(NMAX)=0.D0
          RETURN
          ENDIF
          IF((NF-NI).LE.1) THEN
          LPINF=NI
          LPSUP=NF
          LNINF=IPUNT(NI)
          LNSUP=IPUNT(NF)
          IPUNT(N+1)=0
          HNFORM(NMAX)=0.D0
          RETURN
          ENDIF
          IF(VALORE.LT.HNFORM(LC)) THEN
          NF=NC
          ENDIF
          IF(VALORE.GT.HNFORM(LC)) THEN
          NI=NC
          ENDIF
          GOTO 5
          END
C************************************************************************** 
C Subroutine abinry: Contiene bug quando lavora su grigliati,
C restituisce valori di puntatori nulli!!!!!!!!
C 3 - aprile 1992
C*************************          ABINRY          *************************
c          SUBROUTINE ABINRY(VALORE,HNFORM,IPUNT,NMAX,LNINF,LNSUP,
c     .    LPINF,LPSUP,N)
c          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c          DIMENSION HNFORM(1),IPUNT(1)
c          LPINF=0
c          LPSUP=0
c          LNINF=0
c          LNSUP=0
c          NI=1
c          NF=N
c          IPUNT(N+1)=NMAX
c          HNFORM(NMAX)=1.D10
c          LI=IPUNT(NI)
c          LF=IPUNT(NF)
c          IF(VALORE.LT.HNFORM(LI)) THEN
c          LPINF=0
c          LPSUP=1
c          LNINF=0
c          LNSUP=IPUNT(1)
c          IPUNT(N+1)=0
c          HNFORM(NMAX)=0.D0
c          RETURN
c          ENDIF
c          IF(VALORE.GT.HNFORM(LF)) THEN
c          LPINF=N
c          LPSUP=N+1
c          LNINF=IPUNT(N)
c          LNSUP=NMAX+1
c          IPUNT(N+1)=0
c          HNFORM(NMAX)=0.D0
c          RETURN
c          ENDIF
c          NF=NF+1
c5         NC=(NI+NF)/2
c          LC=IPUNT(NC)
c          IF(VALORE.EQ.HNFORM(LC)) THEN
c          DO 10 I=NC,NF
c          L=IPUNT(I)
c          IF(VALORE.NE.HNFORM(L)) GO TO 12
c          LPSUP=I
c          LNSUP=IPUNT(I)
c10        CONTINUE
c12        DO 15 I=NC,NI,-1
c          L=IPUNT(I)
c          IF(VALORE.NE.HNFORM(L)) GO TO 18
c          LPINF=I
c          LNINF=IPUNT(I)
c15        CONTINUE
c18        IPUNT(N+1)=0
c          HNFORM(NMAX)=0.D0
c          RETURN
c          ENDIF
c          IF((NF-NI).LE.1) THEN
c          LPINF=NI
c          LPSUP=NF
c          LNINF=IPUNT(NI)
c          LNSUP=IPUNT(NF)
c          IPUNT(N+1)=0
c          HNFORM(NMAX)=0.D0
c          RETURN
c          ENDIF
c          IF(VALORE.LT.HNFORM(LC)) THEN
c          NF=NC
c          ENDIF
c          IF(VALORE.GT.HNFORM(LC)) THEN
c          NI=NC
c          ENDIF
c          GOTO 5
c          END
C*************************          ALGITE          *************************
          SUBROUTINE ALGITE(C,T,JC,LC,D,
     .IC,KC,IW,IR,
     .RERMAX,ITEMAX,ITE,ITER,IFINE,R,P,W,N)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION C(1),T(1),JC(1),LC(1),D(1),
     .IC(1),KC(1),IW(1),IR(1),
     .R(1),P(1),W(1)
          IF(ITER.EQ.1) THEN
          NN=LC(N+1)-1
          DO 10 I=1,NN
          T(I)=C(I)
10        CONTINUE
          CALL TCHOLI(T,JC,LC,W,N)
          ENDIF
          CALL CNJGRD(C,T,JC,LC,D,RERMAX,ITEMAX,ITE,R,P,W,N)
          IF(IFINE.NE.1) RETURN
          CALL INVERI(T,JC,LC,IC,KC,W,IW,R,IR,P,N)
          NN=LC(N+1)-1
          DO 20 I=1,NN
          C(I)=T(I)
20        CONTINUE
          RETURN
          END
C*************************          WINDOW          *************************
          SUBROUTINE WINDOW(XC,YC,AMPX,AMPY,X,Y,IPX,IPY,LPX,LPY,
     .IWD,NMAX,NWDMAX,N,NWD,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),IPX(1),IPY(1),LPX(1),LPY(1),IWD(1)
          XMI=XC-AMPX/2.D0
          XMA=XC+AMPX/2.D0
          YMI=YC-AMPY/2.D0
          YMA=YC+AMPY/2.D0
          CALL ABINRY(XMI,X,IPX,NMAX,LNI,LNS,LPI,LPS,N)
          LXI=LPS
          IF(XMI.EQ.X(LNI)) LXI=LPI
          CALL ABINRY(XMA,X,IPX,NMAX,LNI,LNS,LPI,LPS,N)
          LXS=LPI
          IF(XMA.EQ.X(LNS)) LXS=LPS
          CALL ABINRY(YMI,Y,IPY,NMAX,LNI,LNS,LPI,LPS,N)
          LYI=LPS
          IF(YMI.EQ.Y(LNI)) LYI=LPI
          CALL ABINRY(YMA,Y,IPY,NMAX,LNI,LNS,LPI,LPS,N)
          LYS=LPI
          IF(YMA.EQ.Y(LNS)) LYS=LPS
          DO 1 I=1,NWDMAX
          IWD(I)=0
1         CONTINUE
          NXS=LXS-LXI+1
          NYS=LYS-LYI+1
          J=0
          IF(NXS.GE.NYS) THEN
          DO 2 I=LYI,LYS
          L=IPY(I)
          LL=LPX(L)
          IF(LL.GE.LXI.AND.LL.LE.LXS) THEN
          J=J+1
          IWD(J)=IPY(I)
          ENDIF
2         CONTINUE
          NWD=J
          ELSE
          DO 3 I=LXI,LXS
          L=IPX(I)
          LL=LPY(L)
          IF(LL.GE.LYI.AND.LL.LE.LYS) THEN
          J=J+1
          IWD(J)=IPX(I)
          ENDIF
3         CONTINUE
          NWD=J
          ENDIF
          IF(NWD.GT.NWDMAX) THEN
          IF(ILI.EQ.1) THEN          
          WRITE(*,100) XC,YC,AMPX,AMPY,NWD
          ELSE
          WRITE(*,101) XC,YC,AMPX,AMPY,NWD
          ENDIF
          ENDIF         
          RETURN
100       FORMAT(////5X,'NUMERO DEI PUNTI NELLA FINESTRA',/
     .    5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'//
     .    14X,'COORDINATE DEL CENTRO',8X,'AMPIEZZA FINESTRA',5X,
     .    'NUMEROSITA'//5X,2F15.3,5X,2F10.3,5X,I10)
101       FORMAT(////5X,'NUMBER OF POINTS IN THE WINDOW',/
     .    5X,'EXCEEDING THE PREVIOUS DIMENSION'//
     .    10X,'COORDINATES OF THE CENTER',6X,'WIDTH OF THE WINDOW',7X,
     .    'QUANTITY'//5X,2F15.3,5X,2F10.3,5X,I10)
          END
C*************************          POCKET          *************************
          SUBROUTINE POCKET(XC,YC,ZC,AMPX,AMPY,AMPZ,X,Y,Z,
     .    IPX,IPY,IPZ,LPX,LPY,LPZ,IWD,NMAX,NWDMAX,N,NWD,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),Z(1),IPX(1),IPY(1),IPZ(1),
     .    LPX(1),LPY(1),LPZ(1),IWD(1)
          XMI=XC-AMPX/2.D0
          XMA=XC+AMPX/2.D0
          YMI=YC-AMPY/2.D0
          YMA=YC+AMPY/2.D0
          ZMI=ZC-AMPZ/2.D0
          ZMA=ZC+AMPZ/2.D0
          CALL ABINRY(XMI,X,IPX,NMAX,LNI,LNS,LPI,LPS,N)
          LXI=LPS
          IF(XMI.EQ.X(LNI)) LXI=LPI
          CALL ABINRY(XMA,X,IPX,NMAX,LNI,LNS,LPI,LPS,N)
          LXS=LPI
          IF(XMA.EQ.X(LNS)) LXS=LPS
          CALL ABINRY(YMI,Y,IPY,NMAX,LNI,LNS,LPI,LPS,N)
          LYI=LPS
          IF(YMI.EQ.Y(LNI)) LYI=LPI
          CALL ABINRY(YMA,Y,IPY,NMAX,LNI,LNS,LPI,LPS,N)
          LYS=LPI
          IF(YMA.EQ.Y(LNS)) LYS=LPS
          CALL ABINRY(ZMI,Z,IPZ,NMAX,LNI,LNS,LPI,LPS,N)
          LZI=LPS
          IF(ZMI.EQ.Z(LNI)) LZI=LPI
          CALL ABINRY(ZMA,Z,IPZ,NMAX,LNI,LNS,LPI,LPS,N)
          LZS=LPI
          IF(ZMA.EQ.Z(LNS)) LZS=LPS
          DO 1 I=1,NWDMAX
          IWD(I)=0
1         CONTINUE
          NXS=LXS-LXI+1
          NYS=LYS-LYI+1
          NZS=LZS-LZI+1
          J=0
          IF(NXS.GE.NYS.AND.NXS.GE.NZS) THEN
          IF(NZS.GE.NYS) THEN
          DO 10 I=LYI,LYS
          L=IPY(I)
          LL=LPZ(L)
          IF(LL.GE.LZI.AND.LL.LE.LZS) THEN
          KK=LPX(L)
          IF(KK.GE.LXI.AND.KK.LE.LXS) THEN
          J=J+1
          IWD(J)=L
          ENDIF
          ENDIF
10        CONTINUE
          NWD=J
          ELSE
          DO 20 I=LZI,LZS
          L=IPZ(I)
          LL=LPY(L)
          IF(LL.GE.LYI.AND.LL.LE.LYS) THEN
          KK=LPX(L)
          IF(KK.GE.LXI.AND.KK.LE.LXS) THEN          
          J=J+1
          IWD(J)=L
          ENDIF
          ENDIF          
20        CONTINUE
          NWD=J
          ENDIF
          ENDIF
          IF(NYS.GT.NXS.AND.NYS.GE.NZS) THEN
          IF(NZS.GE.NXS) THEN
          DO 30 I=LXI,LXS
          L=IPX(I)
          LL=LPZ(L)
          IF(LL.GE.LZI.AND.LL.LE.LZS) THEN
          KK=LPY(L)
          IF(KK.GE.LYI.AND.KK.LE.LYS) THEN
          J=J+1
          IWD(J)=L
          ENDIF
          ENDIF
30        CONTINUE
          NWD=J
          ELSE
          DO 40 I=LZI,LZS
          L=IPZ(I)
          LL=LPX(L)
          IF(LL.GE.LXI.AND.LL.LE.LXS) THEN
          KK=LPY(L)
          IF(KK.GE.LYI.AND.KK.LE.LYS) THEN          
          J=J+1
          IWD(J)=L
          ENDIF
          ENDIF          
40        CONTINUE
          NWD=J
          ENDIF
          ENDIF
          IF(NZS.GT.NXS.AND.NZS.GT.NYS) THEN
          IF(NYS.GE.NXS) THEN
          DO 50 I=LXI,LXS
          L=IPX(I)
          LL=LPY(L)
          IF(LL.GE.LYI.AND.LL.LE.LYS) THEN
          KK=LPZ(L)
          IF(KK.GE.LZI.AND.KK.LE.LZS) THEN
          J=J+1
          IWD(J)=L
          ENDIF
          ENDIF
50        CONTINUE
          NWD=J
          ELSE
          DO 60 I=LYI,LYS
          L=IPY(I)
          LL=LPX(L)
          IF(LL.GE.LXI.AND.LL.LE.LXS) THEN
          KK=LPZ(L)
          IF(KK.GE.LZI.AND.KK.LE.LZS) THEN          
          J=J+1
          IWD(J)=L
          ENDIF
          ENDIF          
60        CONTINUE
          NWD=J
          ENDIF
          ENDIF
          IF(NWD.GT.NWDMAX) THEN
          IF(ILI.EQ.1) THEN          
          WRITE(*,100) XC,YC,ZC,AMPX,AMPY,AMPZ,NWD
          ELSE
          WRITE(*,101) XC,YC,ZC,AMPX,AMPY,AMPZ,NWD
          ENDIF
          ENDIF         
          RETURN
100       FORMAT(////5X,'NUMERO DEI PUNTI NELLA FINESTRA',/
     .    5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'//
     .    29X,'COORDINATE DEL CENTRO',18X,'AMPIEZZA FINESTRA',5X,
     .    'NUMEROSITA'//5X,3F15.3,5X,3F10.3,5X,I10)
101       FORMAT(////5X,'NUMBER OF POINTS IN THE POCKET',/
     .    5X,'EXCEEDING THE PREVIOUS DIMENSION'//
     .    25X,'COORDINATES OF THE CENTER',16X,'WIDTH OF THE POCKET',7X,
     .    'QUANTITY'//5X,3F15.3,5X,3F10.3,5X,I10)
          END
C*************************          BIORDA          *************************
          SUBROUTINE BIORDA(X,Y,IPX,IPY,LPX,LPY,N)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),IPX(1),IPY(1),LPX(1),LPY(1)
          CALL ORDINA(X,IPX,N,1)
          CALL ORDINA(Y,IPY,N,1)
          CALL PUNTA(IPX,LPX,0,1,N,N)
          CALL PUNTA(IPY,LPY,0,1,N,N)
          RETURN
          ENTRY RBORDA(X,Y,IPX,IPY,LPX,LPY,N,NN)
          CALL RIORDA(X,IPX,N,NN,1)
          CALL RIORDA(Y,IPY,N,NN,1)
          CALL PUNTA(IPX,LPX,0,1,N,NN)
          CALL PUNTA(IPY,LPY,0,1,N,NN)
          RETURN
          END
C*************************          TIORDA          *************************
          SUBROUTINE TIORDA(X,Y,Z,IPX,IPY,IPZ,LPX,LPY,LPZ,N)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),Z(1),IPX(1),IPY(1),IPZ(1),
     .    LPX(1),LPY(1),LPZ(1)
          CALL ORDINA(X,IPX,N,1)
          CALL ORDINA(Y,IPY,N,1)
          CALL ORDINA(Z,IPZ,N,1)
          CALL PUNTA(IPX,LPX,0,1,N,N)
          CALL PUNTA(IPY,LPY,0,1,N,N)
          CALL PUNTA(IPZ,LPZ,0,1,N,N)
          RETURN
          ENTRY RTORDA(X,Y,Z,IPX,IPY,IPZ,LPX,LPY,LPZ,N,NN)
          CALL RIORDA(X,IPX,N,NN,1)
          CALL RIORDA(Y,IPY,N,NN,1)
          CALL RIORDA(Z,IPZ,N,NN,1)
          CALL PUNTA(IPX,LPX,0,1,N,NN)
          CALL PUNTA(IPY,LPY,0,1,N,NN)
          CALL PUNTA(IPZ,LPZ,0,1,N,NN)
          RETURN
          END
C*************************          CNJGRD          *************************
          SUBROUTINE CNJGRD(C,T,JC,LC,D,RERMAX,ITEMAX,ITE,R,P,W,N)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION C(1),T(1),JC(1),LC(1),D(1),R(1),P(1),W(1)
          DO 10 I=1,N
          R(I)=-D(I)
          P(I)=R(I)
10        CONTINUE
          CALL SOLSI(T,JC,LC,P,N)
          S=0.D0
          DO 20 I=1,N
          S=S+R(I)*P(I)
          D(I)=0.D0
20        CONTINUE
          RTR=S
          ITE=0
30        RMAX=0.D0
          DMAX=0.D0
          DO 40 I=1,N
          RMAX=DMAX1(RMAX,DABS(R(I)))
          DMAX=DMAX1(DMAX,DABS(D(I)))
40        CONTINUE
          IF(RMAX.LE.DMAX*RERMAX.OR.ITE.GE.ITEMAX) RETURN
          ITE=ITE+1
          CALL PMATI(C,JC,LC,P,W,N)
          S=0.D0
          DO 50 I=1,N
          S=S+P(I)*W(I)
50        CONTINUE
          PTCP=S
          ALFA=RTR/PTCP
          DO 60 I=1,N
          D(I)=D(I)+ALFA*P(I)
          R(I)=R(I)-ALFA*W(I)
          W(I)=R(I)
60        CONTINUE
          CALL SOLSI(T,JC,LC,W,N)
          S=0.D0
          DO 70 I=1,N
          S=S+R(I)*W(I)
70        CONTINUE
          BETA=S/RTR
          RTR=S
          DO 80 I=1,N
          P(I)=W(I)+BETA*P(I)
80        CONTINUE
          GO TO 30
          END
          SUBROUTINE SOLSI(T,JC,LC,D,N)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION T(1),JC(1),LC(1),D(1)
          DO 20 I=1,N
          I1=LC(I)
          I2=LC(I+1)-2
          S=D(I)
          IF(I1.LE.I2) THEN
          DO 10 JS=I1,I2
          J=JC(JS)
          S=S-T(JS)*D(J)
10        CONTINUE
          ENDIF
          I2=I2+1
          D(I)=S/T(I2)
20        CONTINUE
          DO 30 I=1,N
          IS=LC(I+1)-1
          D(I)=D(I)*T(IS)
30        CONTINUE
          DO 50 I=N,1,-1
          I1=LC(I)
          I2=LC(I+1)-1
          D(I)=D(I)/T(I2)
          I2=I2-1
          IF(I1.LE.I2) THEN
          DS=D(I)
          DO 40 JS=I2,I1,-1
          J=JC(JS)
          D(J)=D(J)-T(JS)*DS
40        CONTINUE
          ENDIF
50        CONTINUE
          RETURN
          END
          SUBROUTINE PMATI(C,JC,LC,P,W,N)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION C(1),JC(1),LC(1),P(1),W(1)
          DO 20 I=1,N
          I1=LC(I)
          I2=LC(I+1)-1
          S=0.D0
          DO 10 JS=I1,I2
          J=JC(JS)
          S=S+C(JS)*P(J)
10        CONTINUE
          W(I)=S
20        CONTINUE
          DO 40 I=N,1,-1
          I1=LC(I)
          I2=LC(I+1)-2
          IF(I1.LE.I2) THEN
          PS=P(I)
          DO 30 JS=I2,I1,-1
          J=JC(JS)
          W(J)=W(J)+C(JS)*PS
30        CONTINUE
          ENDIF
40        CONTINUE
          RETURN
          END
C*************************          DQVV          *************************
          SUBROUTINE DQVV(A,C,P,MMAX,NMAX,M,N)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION A(MMAX,1),C(NMAX,1),P(1)
          DO 10 K=1,M
          T=0.D0
          DO 8 I=1,N
          S=0.D0
          DO 5 J=I,N
          S=S+C(I,J)*A(K,J)
5         CONTINUE
          T=T+2.D0*A(K,I)*S-C(I,I)*A(K,I)**2
8         CONTINUE
          P(K)=1.D0/P(K)-T
10        CONTINUE
          RETURN
          END
C*************************          DQVVB          *************************
          SUBROUTINE DQVVB(A,IA,KA,C,P,MMAX,NMAX,M)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION A(MMAX,1),IA(MMAX,1),KA(1),C(NMAX,1),P(1)
          DO 10 K=1,M
          T=0.D0
          L=KA(K)
          DO 8 IS=1,L
          S=0.D0
          I=IA(K,IS)
          DO 5 JS=IS,L
          J=IA(K,JS)
          II=MIN0(I,J)
          JJ=I+J-II
          S=S+C(II,JJ-II+1)*A(K,JS)
5         CONTINUE
          T=T+2.D0*A(K,IS)*S-C(I,1)*A(K,IS)**2
8         CONTINUE
          P(K)=1.D0/P(K)-T
10        CONTINUE
          RETURN
          END
C*************************          DQVVI          ************************* 
      SUBROUTINE DQVVI(A,IA,KA,C,IC,JC,KC,P,W,MMAX,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MMAX,1),IA(MMAX,1),KA(1),C(1),JC(1),KC(1),IC(1),
     .W(1),P(1)
      DO 10 K=1,M
      T=0.D0
      L=KA(K)
      DO 8 IS=1,L
      S=0.D0
      I=IA(K,IS)
      I1=KC(I)+1
      I2=KC(I+1)
      DO 2 JSS=I1,I2
      JS=IC(JSS)
      J=JC(JS)
      W(J)=C(JS)
    2 CONTINUE
      DO 5 JS=1,L
      J=IA(K,JS)
      IF(J.LT.I) GO TO 5
      S=S+W(J)*A(K,JS)
    5 CONTINUE
      T=T+2.D0*A(K,IS)*S-W(I)*A(K,IS)**2
    8 CONTINUE
      P(K)=1.D0/P(K)-T
   10 CONTINUE
      RETURN
      END
C*************************          FMAT           *************************
          FUNCTION FMAT(IMAX,C,W,NMAX,N)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION C(NMAX,1),W(1)
          DO 10 I=1,N
          W(I)=0.D0
10        CONTINUE
          DO 20 I=1,N-1
          W(I)=W(I)+C(I,I)
          DO 20 J=I+1,N
          W(I)=W(I)+DABS(C(I,J))
          W(J)=W(J)+DABS(C(I,J))
20        CONTINUE
          W(N)=W(N)+C(N,N)
          IMAX=0
          WMAX=0.D0
          DO 30 I=1,N
          IF(WMAX.LT.W(I))THEN
          IMAX=I
          WMAX=W(I)
          ENDIF
30        CONTINUE
          FMAT=WMAX
          RETURN
          END
C*************************          FMATB          *************************
          FUNCTION FMATB(IMAX,C,W,NMAX,N,NB)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION C(NMAX,1),W(1)
          DO 10 I=1,N
          W(I)=0.D0
10        CONTINUE
          DO 20 I=1,N-1
          W(I)=W(I)+C(I,1)
          L=I-1+NB
          L=MIN0(L,N)
          DO 20 J=I+1,L
          W(I)=W(I)+DABS(C(I,J-I+1))
          W(J)=W(J)+DABS(C(I,J-I+1))
20        CONTINUE
          W(N)=W(N)+C(N,1)
          IMAX=0
          WMAX=0.D0
          DO 30 I=1,N
          IF(WMAX.LT.W(I))THEN
          IMAX=I
          WMAX=W(I)
          ENDIF
30        CONTINUE
          FMATB=WMAX
          RETURN
          END
C*************************          FMATI          *************************
          FUNCTION FMATI1(IMAX,C,JC,LC,W,N)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION C(1),JC(1),LC(1),W(1)
          DO 10 I=1,N
          W(I)=0.D0
10        CONTINUE
          W(1)=C(1)
          DO 20 I=2,N
          IS=LC(I+1)-1
          W(I)=W(I)+C(IS)
          INI=LC(I)
          IFI=IS-1
          IF(IFI.GE.INI) THEN
          DO 15 JS=INI,IFI
          J=JC(JS)
          W(I)=W(I)+DABS(C(JS))
          W(J)=W(J)+DABS(C(JS))
15        CONTINUE
          ENDIF
20        CONTINUE
          IMAX=0
          WMAX=0.D0
          DO 30 I=1,N
          IF(WMAX.LT.W(I))THEN
          IMAX=I
          WMAX=W(I)
          ENDIF
30        CONTINUE
          FMATI1=WMAX
          RETURN
          END
          FUNCTION FMATI2(IMAX,C,IC,JC,KC,W,N)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION C(1),IC(1),JC(1),KC(1),W(1)
          DO 40 I=1,N
          W(I)=0.D0
40        CONTINUE
          DO 50 I=1,N-1
          ISS=KC(I)+1
          IS=IC(ISS)
          W(I)=W(I)+C(IS)
          INI=ISS+1
          IFI=KC(I+1)
          IF(IFI.GE.INI) THEN
          DO 45 JSS=INI,IFI
          JS=IC(JSS)
          J=JC(JS)
          W(I)=W(I)+DABS(C(JS))
          W(J)=W(J)+DABS(C(JS))
45        CONTINUE
          ENDIF
50        CONTINUE
          NNS=KC(N)+1
          NS=IC(NNS)
          W(N)=W(N)+C(NS)
          IMAX=0
          WMAX=0.D0
          DO 60 I=1,N
          IF(WMAX.LT.W(I))THEN
          IMAX=I
          WMAX=W(I)
          ENDIF
60        CONTINUE
          FMATI2=WMAX
          RETURN
          END
C*************************          INVER          *************************
      SUBROUTINE INVER(C,W,NMAX,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(NMAX,1),W(1)
      DO 70 I=N,1,-1
      L=I+1
      IF(I.EQ.N) GO TO 40
      DO 30 J=N,L,-1
      S=0.D0
      DO 20 K=N,L,-1
      IF(J.LT.K) GO TO 10
      S=S+C(I,K)*C(K,J)
      GO TO 20
   10 S=S+C(I,K)*C(J,K)
   20 CONTINUE
      W(J)=-S/C(I,I)
   30 CONTINUE
   40 S=0.D0
      IF(I.EQ.N) GO TO 65
      DO 50 K=N,L,-1
      S=S+C(I,K)*W(K)
   50 CONTINUE
      DO 60 J=N,L,-1
      C(I,J)=W(J)
   60 CONTINUE
   65 C(I,I)=(1.D0/C(I,I)-S)/C(I,I)
   70 CONTINUE
      RETURN
      END
C*************************          INVERB          *************************
      SUBROUTINE INVERB(C,W,NMAX,N,NB)
c      VIRTUAL (ALL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION C(NMAX,1),W(1)
      DO 70 I=N,1,-1
      L=I+1
      MINO=MIN0(N,I+NB-1)
      IF(I.EQ.N) GO TO 40
      IF(MINO.LT.L) GO TO 40
      DO 30 J=MINO,L,-1
      S=0.D0
      DO 20 K=MINO,L,-1
      IF(J.LT.K) GO TO 10
      S=S+C(I,K-I+1)*C(K,J-K+1)
      GO TO 20
   10 S=S+C(I,K-I+1)*C(J,K-J+1)
   20 CONTINUE
      W(J)=-S/C(I,1)
   30 CONTINUE
   40 S=0.D0
      IF(I.EQ.N) GO TO 65
      IF(MINO.LT.L) GO TO 65
      DO 50 K=MINO,L,-1
      S=S+C(I,K-I+1)*W(K)
   50 CONTINUE
      DO 60 J=MINO,L,-1
      C(I,J-I+1)=W(J)
   60 CONTINUE
   65 C(I,1)=(1.D0/C(I,1)-S)/C(I,1)
   70 CONTINUE
      RETURN
      END
C*************************          INVERI          *************************
          SUBROUTINE INVERI(T,JC,LC,IC,KC,W,IW,R,IR,Q,N)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION T(1),JC(1),LC(1),IC(1),KC(1),
     .W(1),IW(1),R(1),IR(1),Q(1)
          NN=LC(N+1)-1
          T(NN)=1.D0/T(NN)
          CALL CLASSI(JC,KC,IC,NN,N)
          DO 10 I=1,N
          W(I)=0.D0
          R(I)=0.D0
          Q(I)=0.D0
          I1=LC(I)
          I2=LC(I+1)-1
          DO 10 JS=I1,I2
          JC(JS)=I
10        CONTINUE
          DO 80 I=N-1,1,-1
          I1=KC(I)+1
          I2=KC(I+1)
          IS=IC(I1)
          TT=T(IS)
          QQ=0.D0
          I1=I1+1
          IF(I1.LE.I2) THEN
          KK=0
          DO 20 KS=I1,I2
          JS=IC(KS)
          J=JC(JS)
          KK=KK+1
          IR(KK)=J
          R(J)=T(JS)
20        CONTINUE
          K=0
          JJ=0
          DO 50 KS=I1,I2
          K=K+1
          J=IR(K)
          J1=KC(J)+1
          J2=KC(J+1)
          LL=0
          DO 30 LSS=J1,J2
          LS=IC(LSS)
          L=JC(LS)
          LL=LL+1
          IW(LL)=L
          W(L)=T(LS)
30        CONTINUE
          JJ=MAX0(JJ,L)
          RR=R(J)
          S=0.D0
          DO 40 LS=1,LL
          L=IW(LS)
          S=S+R(L)*W(L)
          Q(L)=Q(L)+RR*W(L)
40        CONTINUE
          JS=IC(KS)
          T(JS)=-(Q(J)-RR*W(J)+S)/TT
          QQ=QQ+RR*T(JS)
          DO 50 LS=1,LL
          L=IW(LS)
          W(L)=0.D0
50        CONTINUE
          DO 60 KS=1,KK
          J=IR(KS)
          R(J)=0.D0
60        CONTINUE
          DO 70 J=I,JJ
          Q(J)=0.D0
70        CONTINUE
          ENDIF
          T(IS)=(1.D0-QQ)/TT
80        CONTINUE
          RETURN
          END
C*************************          NORMAB          *************************
      SUBROUTINE NORMAB(A,IA,KA,B,P,C,D,MMAX,NMAX,NBMAX,M,N,NB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(MMAX,1),IA(MMAX,1),KA(1),B(1),P(1),C(NMAX,1),D(1)
      NB=0
      DO 5 I=1,N
      DO 5 J=1,NBMAX
      C(I,J)=0.D0
    5 CONTINUE
      DO 10 K=1,M
      L=KA(K)
      DO 10 I=1,L
      L1=IA(K,I)
      DO 10 J=I,L
      L2=IA(K,J)
      LC=MIN0(L1,L2)
      LD=IABS(L2-L1)
      NB=MAX0(NB,LD+1)
      IF(NB.GT.NBMAX) GO TO 30
      C(LC,LD+1)=C(LC,LD+1)+A(K,I)*P(K)*A(K,J)
   10 CONTINUE
      GO TO 12
      ENTRY NORTNB(A,IA,KA,B,P,D,MMAX,M,N)
   12 DO 15 I=1,N
      D(I)=0.D0
   15 CONTINUE
      DO 20 K=1,M
      L=KA(K)
      DO 20 I=1,L
      L1=IA(K,I)
      D(L1)=D(L1)+A(K,I)*P(K)*B(K)
   20 CONTINUE
      RETURN
   30 ILI=KA(M+1)
      IF(ILI.EQ.1) THEN
      WRITE(6,100) K
      ELSE
      WRITE(6,101) K
      ENDIF
      STOP
  100 FORMAT(///,5X,'MATRICE NON RIDUCIBILE IN BANDA - EQUAZIONE',I5,///
     .)
  101 FORMAT(///,5X,'MATRIX NOT STORED IN BANDWIDTH TYPE - EQUATION',I5,
     .///)
      END
C*************************          NORMAI          *************************
          SUBROUTINE NORMAI(A,IA,KA,JA,LA,IB,JB,B,P,C,JC,LC,D,W,
     .INFORM,MMAX,M,N)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION A(MMAX,1),IA(MMAX,1),KA(1),JA(1),LA(1),IB(1),JB(1),
     .B(1),P(1),C(1),JC(1),LC(1),D(1),W(1),INFORM(1)
	  CALL BANDA(IA,KA,IB,JB,MMAX,M,N,IT,JT,LB)
          CALL TRASPS (IA,KA,LA,JA,MMAX,M,N)
	  JMIN=0
	  JMAX=0
          DO 1 I=1,N
          INFORM(I)=0
          W(I)=0.D0
1         CONTINUE
          KONT=0
          LC(1)=1
          DO 40 I=1,N
          L1=LA(I)+1
          L2=LA(I+1)
          IF(L1.LE.L2) THEN
	  JMIN=I-IB(I)
	  JMAX=I+JB(I)
          DO 20 KS=L1,L2
          K=JA(KS)
          LL=KA(K)
          DO 10 L=1,LL
          J=IA(K,L)
          INFORM(J)=L
10        CONTINUE
          II=INFORM(I)
          DO 20 L=1,LL
          J=IA(K,L)
          W(J)=W(J)+A(K,II)*P(K)*A(K,L)
20        CONTINUE
          DO 30 J=JMIN,I
          IF(INFORM(J).EQ.0) GOTO 30
          KONT=KONT+1
          JC(KONT)=J
          C(KONT)=W(J)
30        CONTINUE
	  LC(I+1)=KONT+1
	  DO 35 J=JMIN,JMAX
	  INFORM(J)=0
	  W(J)=0.D0
35	  CONTINUE
          ELSE
          KONT=KONT+1
          JC(KONT)=I
          LC(I+1)=KONT+1
	  ENDIF
40	  CONTINUE
          GOTO 50
          ENTRY NORTNI(A,IA,KA,B,P,D,MMAX,M,N)
50        DO 60 I=1,N
          D(I)=0.D0
60        CONTINUE
          DO 70 K=1,M
          LL=KA(K)
          DO 70 L=1,LL
          I=IA(K,L)
          D(I)=D(I)+A(K,L)*P(K)*B(K)
70        CONTINUE
          RETURN
          END
C*************************          NORMAL          *************************
      SUBROUTINE NORMAL(A,B,P,C,D,MMAX,NMAX,M,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(MMAX,1),B(1),P(1),C(NMAX,1),D(1)
      DO 5 I=1,N
      DO 5 J=I,N
      C(I,J)=0.D0
      DO 5 K=1,M
      C(I,J)=C(I,J)+A(K,I)*P(K)*A(K,J)
    5 CONTINUE
      GO TO 8
      ENTRY NORTN(A,B,P,D,M,N)
    8 DO 10 I=1,N
      D(I)=0.D0
      DO 10 K=1,M
      D(I)=D(I)+A(K,I)*P(K)*B(K)
   10 CONTINUE
      RETURN
      END
C*************************          ORDINA          *************************
          SUBROUTINE ORDINA(HNFORM,IPUNT,N,ISGN)
          DOUBLE PRECISION HNFORM
          DIMENSION HNFORM(1),IPUNT(1)
          NS=N
          DO 1 I=1,N
          IPUNT(I)=I
1         CONTINUE
          GO TO 2
          ENTRY RIORDA(HNFORM,IPUNT,N,NN,ISGN)
          NS=NN
2         IF(N.EQ.1) RETURN
          IF(ISGN.NE.1) THEN
          DO 5 I=1,N
          J=IPUNT(I)
          HNFORM(J)=ISGN*HNFORM(J)
5         CONTINUE
          ENDIF
          NM1=N-1
          K=0
          DO 8 I=1,NM1
          J1=IPUNT(I)
          J2=IPUNT(I+1)
          IF(HNFORM(J1).GT.HNFORM(J2)) K=K+1
8         CONTINUE
          IF(K.GT.NM1/2) THEN
          N2=N/2
          DO 10 I=1,N2
          ISRV=IPUNT(I)
          IPUNT(I)=IPUNT(N+1-I)
          IPUNT(N+1-I)=ISRV
10        CONTINUE
          ENDIF
          IF(K.NE.0.AND.K.NE.NM1) THEN
          DO 15 I=N,1,-1
          IPUNT(I+1)=IPUNT(I)
15        CONTINUE
          IPUNT(1)=NS+1
          HNFORM(NS+1)=-1000000000.D0
          DO 18 I=2,N
          J=I
17        L1=IPUNT(J)
          L2=IPUNT(J+1)
          IF(HNFORM(L1).LE.HNFORM(L2)) GO TO 18
          IPUNT(J)=L2
          IPUNT(J+1)=L1
          J=J-1
          GO TO 17
18        CONTINUE
          DO 20 I=1,N
          IPUNT(I)=IPUNT(I+1)
20        CONTINUE
          ENDIF
          IF(ISGN.EQ.1) RETURN
          DO 25 I=1,N
          J=IPUNT(I)
          HNFORM(J)=ISGN*HNFORM(J)
25        CONTINUE
          RETURN
          END
C*************************          SOLS          *************************
      SUBROUTINE SOLS(C,D,NMAX,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(NMAX,1),D(1)
      DO 10 I=1,N
      L=I-1
      S=0.D0
      IF(I.EQ.1) GO TO 8
      DO 5 K=1,L
      S=S-C(K,I)*D(K)
    5 CONTINUE
    8 D(I)=(S+D(I))/C(I,I)
   10 CONTINUE
      DO 20 I=N,1,-1
      II=I+1
      S=0.D0
      IF(I.EQ.N) GO TO 18
      DO 15 K=N,II,-1
      S=S+C(I,K)*D(K)
   15 CONTINUE
   18 D(I)=-(S+D(I))/C(I,I)
   20 CONTINUE
      RETURN
      END
C*************************          SOLSB          *************************
      SUBROUTINE SOLSB(C,D,NMAX,N,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(NMAX,1),D(1)
      DO 10 I=1,N
      L=I-1
      S=0.D0
      IF(I.EQ.1) GO TO 8
      MAXO=MAX0(1,I-NB+1)
      IF(MAXO.GT.L) GO TO 8
      DO 5 K=MAXO,L
      S=S-C(K,I-K+1)*D(K)
    5 CONTINUE
    8 D(I)=(S+D(I))/C(I,1)
   10 CONTINUE
      DO 20 I=N,1,-1
      II=I+1
      S=0.D0
      IF(I.EQ.N) GO TO 18
      MINO=MIN0(N,I+NB-1)
      IF(MINO.LT.II) GO TO 18
      DO 15 K=MINO,II,-1
      S=S+C(I,K-I+1)*D(K)
   15 CONTINUE
   18 D(I)=-(S+D(I))/C(I,1)
   20 CONTINUE
      RETURN
      END
C*************************          SZERO          *************************
          SUBROUTINE SZERO(A,B,P,D,SZ,SZ2,MMAX,M,N,L)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION A(MMAX,1),B(1),P(1),D(1)
          SRV=0.D0
          DO 10 I=1,M
          S=B(I)
          DO 5 J=1,N
          S=S+A(I,J)*D(J)
5         CONTINUE
          B(I)=S
          SRV=SRV+B(I)**2*P(I)
10        CONTINUE
          SZ=0.D0
          SZ2=0.D0
          IF((M-N+L).GT.0) THEN
          SZ2=SRV/(M-N+L)
          SZ=DSQRT(SZ2)
          ENDIF
          RETURN
          END
C*************************          SZEROR          *************************
          SUBROUTINE SZEROR(A,IA,KA,B,P,D,SZ,SZ2,MMAX,M,N,L)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION A(MMAX,1),IA(MMAX,1),KA(1),B(1),P(1),D(1)
          SRV=0.D0
          DO 10 I=1,M
          S=B(I)
          K=KA(I)
          DO 5 J=1,K
          KS=IA(I,J)
          S=S+A(I,J)*D(KS)
5         CONTINUE
          B(I)=S
          SRV=SRV+B(I)**2*P(I)
10        CONTINUE
          SZ=0.D0
          SZ2=0.D0
          IF((M-N+L).GT.0) THEN
          SZ2=SRV/(M-N+L)
          SZ=DSQRT(SZ2)
          ENDIF
          RETURN
          END
C*************************          TCHOL          *************************
      SUBROUTINE TCHOL(C,NMAX,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(NMAX,1)
      DO 40 I=1,N
      L=I-1
      S=0.D0
      IF(L.EQ.0) GO TO 15
      DO 10 K=1,L
      S=S+C(K,I)**2
   10 CONTINUE
   15 C(I,I)=C(I,I)-S
      IF(C(I,I).LE.0.D0) GO TO 50
      C(I,I)=DSQRT(C(I,I))
      IF(I.EQ.N) GO TO 40
      II=I+1
      DO 30 J=II,N
      S=0.D0
      IF(L.EQ.0) GO TO 25
      DO 20 K=1,L
      S=S+C(K,I)*C(K,J)
   20 CONTINUE
   25 C(I,J)=(C(I,J)-S)/C(I,I)
   30 CONTINUE
   40 CONTINUE
      RETURN
   50 ILI=IDNINT(C(N+1,N+1))
      IF(ILI.EQ.1) THEN
      WRITE(*,100) I
      ELSE
      WRITE(6,101) I
      ENDIF
      STOP
  100 FORMAT(///5X,'MATRICE NON DEFINITA POSITIVA - EL. DIAGONALE',I5/)
  101 FORMAT(///5X,'MATRIX NOT POSITIVE DEFINED - DIAGONAL ELEMENT',I5/)
      END
C*************************          TCHOLB          *************************
      SUBROUTINE TCHOLB(C,NMAX,N,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(NMAX,1)
      DO 40 I=1,N
      L=I-1
      S=0.D0
      IF(L.EQ.0) GO TO 15
      MAXO=MAX0(1,I-NB+1)
      IF(MAXO.GT.L) GO TO 15
      DO 10 K=MAXO,L
      S=S+C(K,I-K+1)**2
   10 CONTINUE
   15 C(I,1)=C(I,1)-S
      IF(C(I,1).LE.0.D0) GO TO 50
      C(I,1)=DSQRT(C(I,1))
      IF(NB.EQ.1) GO TO 40
      IF(I.EQ.N) GO TO 40
      II=I+1
      MINO=MIN0(N,I+NB-1)
      IF(MINO.LT.II) GO TO 40
      DO 30 J=II,MINO
      S=0.D0
      IF(L.EQ.0) GO TO 25
      MAXO=MAX0(1,J-NB+1)
      IF(MAXO.GT.L) GO TO 25
      DO 20 K=MAXO,L
      S=S+C(K,I-K+1)*C(K,J-K+1)
   20 CONTINUE
   25 C(I,J-I+1)=(C(I,J-I+1)-S)/C(I,1)
   30 CONTINUE
   40 CONTINUE
      RETURN
   50 ILI=IDNINT(C(N+1,1))
      IF(ILI.EQ.1) THEN
      WRITE(*,100) I
      ELSE
      WRITE(6,101) I
      ENDIF
      STOP
  100 FORMAT(///5X,'MATRICE NON DEFINITA POSITIVA - EL. DIAGONALE',I5/)
  101 FORMAT(///5X,'MATRIX NOT POSITIVE DEFINED - DIAGONAL ELEMENT',I5/)
      END
C*************************          TCHOLI          *************************
          SUBROUTINE TCHOLI(T,JC,LC,W,N)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          DIMENSION T(1),JC(1),LC(1),W(1)
          DO 10 I=1,N
          W(I)=0.D0
10        CONTINUE
          DO 50 I=1,N
          I1=LC(I)
          I2=LC(I+1)-1
          DO 40 JS=I1,I2
          J=JC(JS)
          J1=LC(J)
          J2=LC(J+1)-2
          S=T(JS)
          IF(J1.LE.J2) THEN
          DO 20 KS=J1,J2
          K=JC(KS)
          S=S-T(KS)*W(K)
20        CONTINUE
          ENDIF
          T(JS)=S
          J2=J2+1
       	  W(J)=S/T(J2)
40        CONTINUE
          DO 50 JS=I1,I2
          J=JC(JS)
          W(J)=0.D0
50        CONTINUE
          RETURN
          END
