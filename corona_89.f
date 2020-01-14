          IMPLICIT DOUBLE PRECISION(A-H,O-Z)   
          PARAMETER (NMAX=10000,NWDMAX=10000)  
          CHARACTER*1 DATO,PCHAR,PTO,BLANK 
          CHARACTER*4 FORMA(20),TITOLO(20) 
          DIMENSION X(NMAX),Y(NMAX),U(NMAX),V(NMAX),W(NMAX),
     .    IPXU(NMAX),IPYU(NMAX),LPXU(NMAX),LPYU(NMAX), 
     .    IPXV(NMAX),IPYV(NMAX),LPXV(NMAX),LPYV(NMAX), 
     .    IPXW(NMAX),IPYW(NMAX),LPXW(NMAX),LPYW(NMAX), 
     .    IWD(NWDMAX),IPTO(6),IAST(3),DATO(120),PCHAR(120) 
C         COMMON /COM/ IPXU,IPYU,LPXU,LPYU,IPXV,IPYV,LPXV,LPYV,
C    .IPXW,IPYW,LPXW,LPYW,IWD,X,Y,U,V,W
          DATA PTO/'.'/,BLANK/' '/ 
          DATA IPTO/6*120/,IAST/116,118,120/
          CALL PACORO(TITOLO,ICASO,ITIPO,ISCRI,ITRIM,IQUO,TOLL,
     .AINF,ASUP,STP,SCX,SCY,SCZ,SCU,SCV,SCW,FORMA) 
          WRITE(6,1020) TITOLO 
          DO 5 I=1,120 
          PCHAR(I)=BLANK
5         CONTINUE 
          J=0  
10        READ(1,1000,END=30) (DATO(I),I=1,120)
          WRITE(2,1000) (DATO(I),I=1,120)  
          J=J+1
          DO 20 I=1,120
          IF(DATO(I).EQ.PTO) PCHAR(I)=PTO  
20        CONTINUE 
          GOTO 10  
30        N=J  
          IF(N.GT.NMAX) THEN
          WRITE(6,1060)
          STOP 
          ENDIF
          J=0  
          DO 40 I=1,120
          IF(PCHAR(I).EQ.PTO) THEN 
          J=J+1
          IPTO(J)=I
          ENDIF
40        CONTINUE 
          NPTO=J
          IF(NPTO.NE.(ICASO+ITIPO+IQUO)) THEN  
          WRITE(6,1050)
          print *, npto,icaso,itipo,iquo
          STOP 
          ENDIF
C  
          IF(ISCRI.EQ.1) THEN  
          REWIND 1 
          CALL TESTA(ICASO,ITIPO,IQUO,IPTO,NPTO)
          DO 50 J=1,N  
          READ(1,1000) (DATO(I),I=1,120)
          WRITE(6,1030) (DATO(I),I=1,120)  
50        CONTINUE 
          ENDIF
C  
          IF(IQUO.EQ.1) THEN
          DO 60 I=4,6  
          IPTO(I-1)=IPTO(I)
60        CONTINUE 
          ENDIF
          IF(ICASO.EQ.1) THEN  
                    IF(ITIPO.EQ.1) CALL OTTIT1 
     .              (X,W,IPXW,LPXW,AINF,ASUP,STP,  
     .              NMAX,SCZ,SCW,IPTO,IAST,FORMA,ITRIM,TOLL)
                    IF(ITIPO.EQ.2) CALL OTTIT2 
     .              (X,U,V,IPXU,IPXV,LPXU,LPXV,AINF,ASUP,STP,  
     .              NMAX,SCZ,SCU,SCV,IPTO,IAST,FORMA,ITRIM,TOLL)
                    IF(ITIPO.EQ.3) CALL OTTIT3 
     .              (X,U,V,W,IPXU,IPXV,IPXW,LPXU,LPXV,LPXW,AINF,ASUP,  
     .              STP,NMAX,SCZ,SCU,SCV,SCW,IPTO,IAST,FORMA,ITRIM,TOLL)
          ENDIF
          IF(ICASO.EQ.2) THEN  
                    IF(ITIPO.EQ.1) CALL OTTIP1 
     .              (X,Y,W,IPXW,IPYW,LPXW,LPYW,AINF,ASUP,STP,IWD,  
     .              NMAX,NWDMAX,SCX,SCY,SCW,IPTO,IAST,FORMA,ITRIM,TOLL)
                    IF(ITIPO.EQ.2) CALL OTTIP2 
     .              (X,Y,U,V,IPXU,IPYU,IPXV,IPYV,  
     .              LPXU,LPYU,LPXV,LPYV,AINF,ASUP,STP, 
     .              IWD,NMAX,NWDMAX,SCX,SCY,SCU,SCV,
     .              IPTO,IAST,FORMA,ITRIM,TOLL)
                    IF(ITIPO.EQ.3) CALL OTTIP3 
     .              (X,Y,U,V,W,IPXU,IPYU,IPXV,IPYV,IPXW,IPYW,  
     .              LPXU,LPYU,LPXV,LPYV,LPXW,LPYW,AINF,ASUP,STP,
     .              IWD,NMAX,NWDMAX,SCX,SCY,SCU,SCV,SCW,
     .              IPTO,IAST,FORMA,ITRIM,TOLL)
          ENDIF
          IF(ITRIM.EQ.1) THEN  
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,1040) ALFA
          ENDIF
c          STOP 
1000      FORMAT(120A1)
1010      FORMAT(20A4) 
1020      FORMAT(//5X,20A4,//) 
1030      FORMAT(5X,120A1) 
1040      FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEL TEST '
     .'DI REIEZIONE DEGLI ERRORI GROSSOLANI: ALFA =',F5.3, 
     .' (SU DUE CODE)')
1050      FORMAT(5X,'CASO E/O TIPO DEI DATI ERRATO'///)
1060      FORMAT(5X,'NUMERO DEI DATI ',/5X,
     .'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'///)  
          END  
          SUBROUTINE BMAX(XCOVT,TEMPO,IPUNT,NTIPO,LASTP,NMAX,N,MB) 
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          DIMENSION TEMPO(1),IPUNT(1),NTIPO(1),LASTP(1),IND(7) 
          DATA IND/1,1,2,1,2,2,3/  
          MB=0 
          DO 20 I=1,N  
          TMAX=TEMPO(I)+2.D0*XCOVT 
          CALL ABINRY(TMAX,TEMPO,IPUNT,NMAX,LNI,LNS,LPI,LPS,N) 
          ISRV=0
          DO 10 J=I,LNI
          JJ=NTIPO(J)  
          ISRV=ISRV+IND(JJ)
10        CONTINUE 
          MB=MAX0(MB,ISRV) 
          LASTP(I)=LNI 
20        CONTINUE 
          RETURN
          END  
          SUBROUTINE CALCO1(DCOV,DXCOV,DA,DB)  
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*)
          DHELP = DXCOV(3)-DXCOV(2)
          DB = (DLOG(DCOV(2))-DLOG(DCOV(3))) / DHELP
          IF(DB.LT.0.D0) DB=0.D0
          DA = DCOV(2) * DEXP(DB*DXCOV(2)) 
          RETURN
          END  


          SUBROUTINE CALCO2(DCOV,DXCOV,DA,DB)  
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*)
          DHELP = DXCOV(3)**2-DXCOV(2)**2  
          DB = (DLOG(DCOV(2))-DLOG(DCOV(3))) / DHELP
          IF(DB.LT.0.D0) DB=0.D0
          DA = DCOV(2) * DEXP(DB*DXCOV(2)**2)  
          RETURN
          END  


          SUBROUTINE CALCO3(DCOV,DXCOV,DX0,IXMM,DA,DB,DD,ICASO)
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*),DX0(*)
          DIMENSION IXMM(*)
          DD = 1.0D0 / DX0(1)**2
          DHELP1 = (1.0D0 - DD*DXCOV(2)**2)
          DHELP2 = (1.0D0 - DD*DXCOV(IXMM(1))**2)  
          DHELP = DXCOV(IXMM(1)) - DXCOV(2)
          DB = DLOG(DCOV(2)/DHELP1) - DLOG(DCOV(IXMM(1))/DHELP2)
          DB = DB / DHELP  
          DHELP3 = DSQRT(2.D0*DD)  
          IF(DB.LT.DHELP3.AND.ICASO.EQ.1) DB = DHELP3  
          DHELP4 = DSQRT(6.D0*DD)  
          IF(DB.LT.DHELP4.AND.ICASO.EQ.2) DB = DHELP4  
          DA = DCOV(2) * DEXP(DB*DXCOV(2)) / DHELP1
          RETURN
          END  


          SUBROUTINE CALCO4(DCOV,DXCOV,DX0,IXMM,DA,DB,DD,ICASO)
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*),DX0(*)
          DIMENSION IXMM  (*)  
          DD = 1.D0 / DX0(1)**2
          DHELP1 = (1.D0 - DD*DXCOV(2)**2) 
          DHELP2 = (1.D0 - DD * DXCOV(IXMM(1))**2) 
          DHELP = DXCOV(IXMM(1))**2-DXCOV(2)**2
          DB = DLOG(DCOV(2)/DHELP1) - DLOG(DCOV(IXMM(1))/DHELP2)
          DB = DB / DHELP  
          DHELP3 = DD/2.D0 
          IF(DB.LT.DHELP3.AND.ICASO.EQ.1) DB = DHELP3  
          DHELP4 = DD  
          IF(DB.LT.DHELP4.AND.ICASO.EQ.2) DB = DHELP4  
          DA = DCOV(2) * DEXP(DB*DXCOV(2)**2) / DHELP1 
          RETURN
          END  


          SUBROUTINE CALCO5(DCOV,DXCOV,DX0,IXMM,IX0,DA,DB,DC)  
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*),DX0(*)
          DIMENSION IXMM(*),IX0(*) 
          DPI = 4.D0*DATAN(1.D0)
          DC = (DPI/DX0(1) + 2.D0*DPI/DX0(2)) / 2.D0
          IF(IX0(3).NE.0) DC = (2.D0*DC + 3.D0*DPI/DX0(3)) / 3.D0  
          DHELP = DPI / DC 
          DB = DLOG(-DCOV(IXMM(1)))-DLOG(DCOV(IXMM(2)))
          DB = DB / DHELP  
          IF(DB.LT.0.D0) DB=0.D0
          DA = DCOV(2) * DEXP(DB*DXCOV(2)) 
          DA = DA / (DSIN(DC*DXCOV(2)) / (DC*DXCOV(2)))
          RETURN
          END  


          SUBROUTINE CALCO6(DCOV,DXCOV,DX0,IXMM,IX0,DA,DB,DC)  
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*),DX0(*)
          DIMENSION IXMM(*),IX0(*) 
          DPI = 4.D0*DATAN(1.D0)
          DC = (DPI/DX0(1) + 2.D0*DPI/DX0(2)) / 2.D0
          IF(IX0(3).NE.0) DC = (2.D0*DC + 3.D0*DPI/DX0(3)) / 3.D0  
          DHELP = 4.D0*DPI**2  
          DHELP = DHELP / DC**2
          DB = DLOG(-DCOV(IXMM(1))) - DLOG(DCOV(IXMM(2)))  
          DB = DB / DHELP  
          IF(DB.LT.0.D0) DB=0.D0
          DA = DCOV(2) * DEXP(DB*DXCOV(2)**2)  
          DA = DA / (DSIN(DC*DXCOV(2)) / (DC*DXCOV(2)))
          RETURN
          END  


          SUBROUTINE CALCO7(DCOV,DXCOV,DX0,IXMM,IX0,DA,DB,DC)  
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*),DX0(*)
          DIMENSION IXMM(*),IX0(*) 
          DPI = 4.D0*DATAN(1.D0)
          DC = (DPI/(2.D0*DX0(1)) + 3.D0*DPI/(2.D0*DX0(2))) / 2.D0 
          IF(IX0(3).NE.0) DC= (2.D0*DC + 5.D0*DPI/(2.D0*DX0(3)))/ 3.D0 
          DHELP = DPI / DC 
          DB = DLOG(-DCOV(IXMM(1)))- DLOG(DCOV(IXMM(2)))
          DB = DB / DHELP  
          IF(DB.LT.0.D0) DB=0.D0
          DA = DCOV(2) * DEXP(DB*DXCOV(2))/(DCOS(DC*DXCOV(2))) 
          RETURN
          END  


          SUBROUTINE CALCO8(DCOV,DXCOV,DX0,IXMM,IX0,DA,DB,DC)  
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*),DX0(*)
          DIMENSION IXMM(*),IX0(*) 
          DPI = 4.D0*DATAN(1.D0)
          DC = (DPI/(2.D0*DX0(1)) + 3.D0*DPI/(2.D0*DX0(2))) / 2.D0 
          IF(IX0(3).NE.0) DC= (2.D0*DC + 5.D0*DPI/(2.D0*DX0(3)))/ 3.D0 
          DHELP = 3.D0*DPI**2  
          DHELP = DHELP / DC**2
          DB = DLOG(-DCOV(IXMM(1)))-DLOG(DCOV(IXMM(2)))
          DB = DB / DHELP  
          IF(DB.LT.0.D0) DB=0.D0
          DA = DCOV(2) * DEXP(DB*DXCOV(2)**2)/(DCOS(DC*DXCOV(2)))  
          RETURN
          END  


          SUBROUTINE CALCO9(DCOV,DXCOV,DX0,IXMM,IX0,DA,DB,DC)  
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*),DX0(*),DX0J1(3),DXMMJ1(2) 
          DIMENSION DMMJ1(2)
          DIMENSION IXMM(*),IX0(*) 
          DATA DX0J1/3.8317060D0,7.0155867D0,10.1734681D0/ 
          DATA DXMMJ1/5.331D0,8.536D0/ 
          DMMJ1(1)=DFJ1(DXMMJ1(1)) 
          DMMJ1(2)=DFJ1(DXMMJ1(2)) 
          DC = (DX0J1(1)/DX0(1) + DX0J1(2)/DX0(2)) / 2.D0  
          IF(IX0(3).NE.0) DC = (2.D0*DC + DX0J1(3)/DX0(3)) / 3.D0  
          DHELP = (DXMMJ1(2)-DXMMJ1(1)) / DC
          DB = DLOG(DCOV(IXMM(1))/DMMJ1(1))-DLOG(DCOV(IXMM(2))/DMMJ1(2))
          DB = DB / DHELP  
          IF(DB.LT.0.D0) DB=0.D0
          DA = DCOV(2) * DEXP(DB*DXCOV(2)) 
          DA = DA / (DFJ1(DC*DXCOV(2))/(DC*DXCOV(2)))  
          RETURN
          END  


          SUBROUTINE CALC10(DCOV,DXCOV,DX0,IXMM,IX0,DA,DB,DC)  
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*),DX0(*),DX0J1(3),DXMMJ1(2) 
          DIMENSION DMMJ1(2)
          DIMENSION IXMM(*),IX0(*) 
          DATA DX0J1/3.8317060D0,7.0155867D0,10.1734681D0/ 
          DATA DXMMJ1/5.331D0,8.536D0/ 
          DMMJ1(1)=DFJ1(DXMMJ1(1)) 
          DMMJ1(2)=DFJ1(DXMMJ1(2)) 
          DC = (DX0J1(1)/DX0(1) + DX0J1(2)/DX0(2)) / 2.D0  
          IF(IX0(3).NE.0) DC = (2.D0*DC + DX0J1(3)/DX0(3)) / 3.D0  
          DHELP = (DXMMJ1(2)**2-DXMMJ1(1)**2) / DC**2  
          DB = DLOG(DCOV(IXMM(1))/DMMJ1(1))-DLOG(DCOV(IXMM(2))/DMMJ1(2))
          DB = DB / DHELP  
          IF(DB.LT.0.D0) DB=0.D0
          DA = DCOV(2) * DEXP(DB*DXCOV(2)**2)  
          DA = DA / (DFJ1(DC*DXCOV(2)) / (DC*DXCOV(2)))
          RETURN
          END  


          SUBROUTINE CALC11(DCOV,DXCOV,DX0,IXMM,IX0,DA,DB,DC)  
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*),DX0(*),DX0J0(3),DXMMJ0(2) 
          DIMENSION DMMJ0(2)
          DIMENSION IXMM(*),IX0(*) 
          DATA DX0J0/2.4048256D0,5.5200781D0,8.6537279D0/  
          DATA DXMMJ0/3.8317060D0,7.0155867D0/ 
          DMMJ0(1)=DFJ0(DXMMJ0(1)) 
          DMMJ0(2)=DFJ0(DXMMJ0(2)) 
          DC = (DX0J0(1)/DX0(1) + DX0J0(2)/DX0(2)) / 2.D0  
          IF(IX0(3).NE.0) DC = (2.D0*DC + DX0J0(3)/DX0(3)) / 3.D0  
          DHELP = (DXMMJ0(2)-DXMMJ0(1)) / DC
          DB = DLOG(DCOV(IXMM(1))/DMMJ0(1))-DLOG(DCOV(IXMM(2))/DMMJ0(2))
          DB = DB / DHELP  
          IF(DB.LT.0.D0) DB=0.D0
          DA = DCOV(2) * DEXP(DB*DXCOV(2))/(DFJ0(DC*DXCOV(2))) 
          RETURN
          END  


          SUBROUTINE CALC12(DCOV,DXCOV,DX0,IXMM,IX0,DA,DB,DC)  
          IMPLICIT DOUBLE PRECISION (D)
          IMPLICIT INTEGER          (I)
          DIMENSION DCOV(*),DXCOV(*),DX0(*),DX0J0(3),DXMMJ0(2) 
          DIMENSION DMMJ0(2)
          DIMENSION IXMM(*),IX0(*) 
          DATA DX0J0/2.4048256D0,5.5200781D0,8.6537279D0/  
          DATA DXMMJ0/3.8317060D0,7.0155867D0/ 
          DMMJ0(1)=DFJ0(DXMMJ0(1)) 
          DMMJ0(2)=DFJ0(DXMMJ0(2)) 
          DC = (DX0J0(1)/DX0(1) + DX0J0(2)/DX0(2)) / 2.D0  
          IF(IX0(3).NE.0) DC = (2.D0*DC + DX0J0(3)/DX0(3)) / 3.D0  
          DHELP = (DXMMJ0(2)**2-DXMMJ0(1)**2) / DC**2  
          DB = DLOG(DCOV(IXMM(1))/DMMJ0(1))-DLOG(DCOV(IXMM(2))/DMMJ0(2))
          DB = DB / DHELP  
          IF(DB.LT.0.D0) DB=0.D0
          DA = DCOV(2) * DEXP(DB*DXCOV(2)**2)/(DFJ0(DC*DXCOV(2)))  
          RETURN
          END  
          SUBROUTINE CODA(I,NPT,II,D,AA,E,BB,CC,IFILE,TOLL)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          CHARACTER*1 AST,BLK,TX,TY,TZ 
          DIMENSION AA(3,1),BB(3,1),CC(3,1),D(1),E(1)  
          DATA AST/'*'/,BLK/' '/
          SRV=1.D0 
          GO TO (10,20,30,40,50,60,70) II  
10        TX=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TX=AST
          WRITE(6,100) I,NPT,D(1),AA(1,1),SRV,E(1),TX,BB(1,1),SRV, 
     .CC(1,1),SRV  
          IF(IFILE.EQ.1) WRITE(11,200) NPT,D(1),AA(1,1),
     .E(1),BB(1,1),TX  
          RETURN
20        TY=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TY=AST
          WRITE(6,110) I,NPT,D(1),AA(1,1),SRV,E(1),TY,BB(1,1),SRV, 
     .CC(1,1),SRV  
          IF(IFILE.EQ.1) WRITE(11,210) NPT,D(1),AA(1,1),
     .E(1),BB(1,1),TY  
          RETURN
30        TX=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TX=AST
          TY=BLK
          IF(DABS(E(2))/BB(2,2).GT.TOLL) TY=AST
          WRITE(6,120) I,NPT,D(1),AA(1,1),SRV,AA(1,2), 
     .E(1),TX,BB(1,1),SRV,BB(1,2),CC(1,1),SRV,CC(1,2)  
          WRITE(6,130) D(2),AA(2,2),SRV,E(2),TY,BB(2,2),SRV,CC(2,2),SRV
          IF(IFILE.EQ.1) WRITE(11,220) NPT,D(1),AA(1,1),D(2),AA(2,2),  
     .E(1),BB(1,1),E(2),BB(2,2),TX,TY  
          RETURN
40        TZ=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TZ=AST
          WRITE(6,140) I,NPT,D(1),AA(1,1),SRV,E(1),TZ,BB(1,1),SRV, 
     .CC(1,1),SRV  
          IF(IFILE.EQ.1) WRITE(11,230) NPT,D(1),AA(1,1),
     .E(1),BB(1,1),TZ  
          RETURN
50        TX=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TX=AST
          TZ=BLK
          IF(DABS(E(2))/BB(2,2).GT.TOLL) TZ=AST
          WRITE(6,150) I,NPT,D(1),AA(1,1),SRV,AA(1,2), 
     .E(1),TX,BB(1,1),SRV,BB(1,2),CC(1,1),SRV,CC(1,2)  
          WRITE(6,160) D(2),AA(2,2),SRV,E(2),TZ,BB(2,2),SRV,CC(2,2),SRV
          IF(IFILE.EQ.1) WRITE(11,240) NPT,D(1),AA(1,1),D(2),AA(2,2),  
     .E(1),BB(1,1),E(2),BB(2,2),TX,TZ  
          RETURN
60        TY=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TY=AST
          TZ=BLK
          IF(DABS(E(2))/BB(2,2).GT.TOLL) TZ=AST
          WRITE(6,170) I,NPT,D(1),AA(1,1),SRV,AA(1,2), 
     .E(1),TY,BB(1,1),SRV,BB(1,2),CC(1,1),SRV,CC(1,2)  
          WRITE(6,160) D(2),AA(2,2),SRV,E(2),TZ,BB(2,2),SRV,CC(2,2),SRV
          IF(IFILE.EQ.1) WRITE(11,250) NPT,D(1),AA(1,1),D(2),AA(2,2),  
     .E(1),BB(1,1),E(2),BB(2,2),TY,TZ  
          RETURN
70        TX=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TX=AST
          TY=BLK
          IF(DABS(E(2))/BB(2,2).GT.TOLL) TY=AST
          TZ=BLK
          IF(DABS(E(3))/BB(3,3).GT.TOLL) TZ=AST
          WRITE(6,190) I,NPT,D(1),AA(1,1),SRV,AA(1,2),AA(1,3), 
     .E(1),TX,BB(1,1),SRV,BB(1,2),BB(1,3),CC(1,1),SRV,CC(1,2),CC(1,3)  
          WRITE(6,180) D(2),AA(2,2),SRV,AA(2,3),
     .E(2),TY,BB(2,2),SRV,BB(2,3),CC(2,2),SRV,CC(2,3)  
          WRITE(6,160) D(3),AA(3,3),SRV,E(3),TZ,BB(3,3),SRV,CC(3,3),SRV
          IF(IFILE.EQ.1) WRITE(11,260) NPT,D(1),AA(1,1),D(2),AA(2,2),  
     .D(3),AA(3,3),E(1),BB(1,1),E(2),BB(2,2),E(3),BB(3,3),TX,TY,TZ 
          RETURN
100       FORMAT(5X,2I5,5X,2F10.3,F5.2,15X,F10.3,2X,A1,2X, 
     .F10.3,F5.2,15X,F10.3,F5.2/)  
110       FORMAT(5X,2I5,5X,2F10.3,5X,F5.2,10X,F10.3,2X,A1,2X,  
     .F10.3,5X,F5.2,10X,F10.3,5X,F5.2/)
120       FORMAT(5X,2I5,5X,2F10.3,2F5.2,10X,F10.3,2X,A1,2X,
     .F10.3,2F5.2,10X,F10.3,2F5.2) 
130       FORMAT(20X,2F10.3,5X,F5.2,10X,F10.3,2X,A1,2X,
     .F10.3,5X,F5.2,10X,F10.3,5X,F5.2/)
140       FORMAT(5X,2I5,5X,2F10.3,10X,F5.2,5X,F10.3,2X,A1,2X,  
     .F10.3,10X,F5.2,5X,F10.3,10X,F5.2/)
150       FORMAT(5X,2I5,5X,2F10.3,F5.2,5X,F5.2,5X,F10.3,2X,A1,2X,  
     .F10.3,F5.2,5X,F5.2,5X,F10.3,F5.2,5X,F5.2)
160       FORMAT(20X,2F10.3,10X,F5.2,5X,F10.3,2X,A1,2X,
     .F10.3,10X,F5.2,5X,F10.3,10X,F5.2/)
170       FORMAT(5X,2I5,5X,2F10.3,5X,2F5.2,5X,F10.3,2X,A1,2X,  
     .F10.3,5X,2F5.2,5X,F10.3,5X,2F5.2)
180       FORMAT(20X,2F10.3,5X,2F5.2,5X,F10.3,2X,A1,2X,
     .F10.3,5X,2F5.2,5X,F10.3,5X,2F5.2)
190       FORMAT(5X,2I5,5X,2F10.3,3F5.2,5X,F10.3,2X,A1,2X, 
     .F10.3,3F5.2,5X,F10.3,3F5.2)  
200       FORMAT(I5,2F10.3,40X,2F10.3,42X,A1)  
210       FORMAT(I5,20X,2F10.3,40X,2F10.3,23X,A1)  
220       FORMAT(I5,4F10.3,20X,4F10.3,22X,2A1) 
230       FORMAT(I5,40X,2F10.3,40X,2F10.3,4X,A1)
240       FORMAT(I5,2F10.3,20X,4F10.3,20X,2F10.3,2X,A1,1X,A1)  
250       FORMAT(I5,20X,4F10.3,20X,4F10.3,3X,2A1)  
260       FORMAT(I5,12F10.3,2X,3A1)
          END  
          SUBROUTINE CORRP1(X,Y,W,IPXW,IPYW,LPXW,LPYW,COV,SRV,NSRV,NCOR,
     .    RAGGIO,ALATO,IWD,NMAX,NWDMAX,SCX,SCY,SCW,
     .    IPTO,IAST,FORMA,IFILE)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20)
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),W(1),IPXW(1),IPYW(1),
     .    LPXW(1),LPYW(1),IWD(1),COV(1),SRV(1),NSRV(1),IPTO(1),IAST(1),
     .    DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
C         INIZIO CASO=2 TIPO=1 VARIABILI X,Y,W 
          SRV1=0.D0
          SRV2=1.D0
          ND=0 
          NW=0 
          REWIND 1 
10        READ(1,1000) (DATO(K),K=1,120)
C          DECODE(150,FORMA,DATO) NPT,XS,YS,WS  
          READ(CDATO,FORMA) NPT,XS,YS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY 
          W(ND)=WS*SCW 
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NW=NW+1  
          IPXW(NW)=ND  
          IPYW(NW)=ND  
          ENDIF
          GOTO 10  
20        CALL MEDVA2(W,IPXW,EW,S2W,SW,NW) 
          DO 30 I=1,ND 
          W(I)=W(I)-EW 
30        CONTINUE 
          CALL RBORDA(X,Y,IPXW,IPYW,LPXW,LPYW,NW,ND)
          CALL COVAP(X,Y,W,IPXW,IPYW,LPXW,LPYW,NW,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)  
          WRITE(6,1010) NW,EW,SW
          WRITE(6,1020) SRV1,S2W,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2W,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1050) NCOR
          WRITE(2,1050) NW,EW,SW
          WRITE(2,1040) SRV1,S2W,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2W,J=1,NCOR)
          ENDIF
          WRITE(6,1060)
          SRV(1)=1.D0  
          DO 40 I=1,NCOR
          SRV(I+1)=COV(I)/S2W  
40        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1060)
C         FINE CASO=2 TIPO=1 VARIABILI X,Y,W
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(I5,F15.6,F10.3)
1060      FORMAT(/////)
          END  
          SUBROUTINE CORRP2(X,Y,U,V,IPXU,IPYU,LPXU,LPYU,
     .    IPXV,IPYV,LPXV,LPYV,COV,SRV,NSRV,
     .    NCOR,RAGGIO,ALATO,IWD,NMAX,NWDMAX,SCX,SCY,SCU,SCV,
     .    IPTO,IAST,FORMA,IFILE)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20)
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),U(1),V(1),IPXU(1),IPYU(1),
     .    LPXU(1),LPYU(1),IPXV(1),IPYV(1),LPXV(1),LPYV(1), 
     .    IWD(1),COV(1),SRV(1),NSRV(1),IPTO(1),IAST(1),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          REWIND 1 
          SRV1=0.D0
          SRV2=1.D0
          ND=0 
          NU=0 
          NV=0 
10        READ(1,1000) (DATO(K),K=1,120)
C          DECODE(150,FORMA,DATO) NPT,XS,YS,US,VS
          READ(CDATO,FORMA) NPT,XS,YS,US,VS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPXU(NU)=ND  
          IPYU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPXV(NV)=ND  
          IPYV(NV)=ND  
          ENDIF
          GOTO 10  
20        CALL MEDVA2(U,IPXU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPXV,EV,S2V,SV,NV) 
          CALL RBORDA(X,Y,IPXU,IPYU,LPXU,LPYU,NU,ND)
          CALL RBORDA(X,Y,IPXV,IPYV,LPXV,LPYV,NV,ND)
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
30        CONTINUE 
C         INIZIO CASO=2 TIPO=1 VARIABILI X,Y,U 
          CALL COVAP(X,Y,U,IPXU,IPYU,LPXU,LPYU,NU,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)  
          WRITE(6,1010) NU,EU,SU
          WRITE(6,1020) SRV1,S2U,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2U,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1090) NCOR
          WRITE(2,1090) NU,EU,SU
          WRITE(2,1040) SRV1,S2U,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2U,J=1,NCOR)
          ENDIF
          WRITE(6,1100)
          SRV(1)=1.D0  
          DO 35 I=1,NCOR
          SRV(I+1)=COV(I)/S2U  
35        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1100)
C         FINE CASO=2 TIPO=1 VARIABILI X,Y,U
C         INIZIO CASO=2 TIPO=1 VARIABILI X,Y,V 
          CALL COVAP(X,Y,V,IPXV,IPYV,LPXV,LPYV,NV,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)  
          WRITE(6,1050) NV,EV,SV
          WRITE(6,1060) SRV1,S2V,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2V,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1090) NCOR
          WRITE(2,1090) NV,EV,SV
          WRITE(2,1040) SRV1,S2V,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2V,J=1,NCOR)
          ENDIF
          WRITE(6,1100)
          SRV(1)=1.D0  
          DO 37 I=1,NCOR
          SRV(I+1)=COV(I)/S2V  
37        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1100)
C         FINE CASO=2 TIPO=1 VARIABILI X,Y,V
          DO 40 I=1,ND 
          U(I)=U(I)+EU 
          V(I)=V(I)+EV 
40        CONTINUE 
C         INIZIO CASO=2 TIPO=2 VARIABILI X,Y,U,V
          CALL RCORR2(U,V,LPXU,LPXV,EU,EV,SU,SV,S2U,S2V,SUV,RUV,ND)
          DO 50 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
50        CONTINUE 
          CALL CROSSP(X,Y,U,V,IPXU,IPYU,LPXU,LPYU,IPXV,IPYV,LPXV,LPYV, 
     .    NU,NV,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)
          SPUV=SU*SV
          WRITE(6,1070) NU,EU,SU,NV,EV,SV,RUV  
          WRITE(6,1080) SRV1,SUV,RUV
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUV,J=1,NCOR)  
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1090) NCOR
          WRITE(2,1090) NU,EU,SU,NV,EV,SV,RUV  
          WRITE(2,1040) SRV1,SUV,RUV
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUV,J=1,NCOR)  
          ENDIF
          WRITE(6,1100)
          SRV(1)=RUV
          DO 65 I=1,NCOR
          SRV(I+1)=COV(I)/SPUV 
65        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1100)
C         FINE CASO=2 TIPO=2 VARIABILI X,Y,U,V 
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y)=',F10.3,5X,'S(Y)=',F10.3,//)
1060      FORMAT(//10X,'DELTA',7X,'GAMMA(Y)',5X,'RO(Y)',///,
     .    5X,F10.3,F15.6,F10.3)
1070      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//5X
     .    ,'N(Y) =',I5,5X,'M(Y)=',F10.3,5X,'S(Y)=',F10.3,5X,'R(X,Y)=', 
     .    F10.2,//)
1080      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Y)',3X,'RO(X,Y)',///5X, 
     .    F10.3,F15.6,F10.3)
1090      FORMAT(I5,F15.6,F10.3,5X,I5,F15.6,F10.3,F10.2)
1100      FORMAT(/////)
          END  
          SUBROUTINE CORRP3(X,Y,U,V,W,IPXU,IPYU,LPXU,LPYU, 
     .    IPXV,IPYV,LPXV,LPYV,IPXW,IPYW,LPXW,LPYW,COV,SRV,NSRV,
     .    NCOR,RAGGIO,ALATO,IWD,NMAX,NWDMAX,SCX,SCY,SCU,SCV,SCW,
     .    IPTO,IAST,FORMA,IFILE)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20)
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),U(1),V(1),W(1),IPXU(1),IPYU(1),  
     .    LPXU(1),LPYU(1),IPXV(1),IPYV(1),LPXV(1),LPYV(1), 
     .    IPXW(1),IPYW(1),LPXW(1),LPYW(1),IWD(1),COV(1),SRV(1),NSRV(1),
     .    IPTO(1),IAST(1),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          SRV1=0.D0
          SRV2=1.D0
          ND=0 
          NU=0 
          NV=0 
          NW=0 
          REWIND 1 
10        READ(1,1000) (DATO(K),K=1,120)
C          DECODE(150,FORMA,DATO) NPT,XS,YS,US,VS,WS
          READ(CDATO,FORMA) NPT,XS,YS,US,VS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          W(ND)=WS*SCW 
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPXU(NU)=ND  
          IPYU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPXV(NV)=ND  
          IPYV(NV)=ND  
          ENDIF
          IF(DATO(IPTO(5)).EQ.PTO.AND.DATO(IAST(3)).NE.AST) THEN
          NW=NW+1  
          IPXW(NW)=ND  
          IPYW(NW)=ND  
          ENDIF
          GOTO 10  
20        CALL MEDVA2(U,IPXU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPXV,EV,S2V,SV,NV) 
          CALL MEDVA2(W,IPXW,EW,S2W,SW,NW) 
          CALL RBORDA(X,Y,IPXU,IPYU,LPXU,LPYU,NU,ND)
          CALL RBORDA(X,Y,IPXV,IPYV,LPXV,LPYV,NV,ND)
          CALL RBORDA(X,Y,IPXW,IPYW,LPXW,LPYW,NW,ND)
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
30        CONTINUE 
C         INIZIO CASO=2 TIPO=1 VARIABILI X,Y,U 
          CALL COVAP(X,Y,U,IPXU,IPYU,LPXU,LPYU,NU,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)  
          WRITE(6,1010) NU,EU,SU
          WRITE(6,1020) SRV1,S2U,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2U,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NU,EU,SU
          WRITE(2,1040) SRV1,S2U,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2U,J=1,NCOR)
          ENDIF
          WRITE(6,1160)
          SRV(1)=1.D0  
          DO 35 I=1,NCOR
          SRV(I+1)=COV(I)/S2U  
35        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=2 TIPO=1 VARIABILI X,Y,U
C         INIZIO CASO=2 TIPO=1 VARIABILI X,Y,V 
          CALL COVAP(X,Y,V,IPXV,IPYV,LPXV,LPYV,NV,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)  
          WRITE(6,1050) NV,EV,SV
          WRITE(6,1060) SRV1,S2V,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2V,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NV,EV,SV
          WRITE(2,1040) SRV1,S2V,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2V,J=1,NCOR)
          ENDIF
          WRITE(6,1160)
          SRV(1)=1.D0  
          DO 37 I=1,NCOR
          SRV(I+1)=COV(I)/S2V  
37        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=2 TIPO=1 VARIABILI X,Y,V
          DO 40 I=1,ND 
          U(I)=U(I)+EU 
          V(I)=V(I)+EV 
40        CONTINUE 
C         INIZIO CASO=2 TIPO=2 VARIABILI X,Y,U,V
          CALL RCORR2(U,V,LPXU,LPXV,EU,EV,SU,SV,S2U,S2V,SUV,RUV,ND)
          DO 50 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
50        CONTINUE 
          CALL CROSSP(X,Y,U,V,IPXU,IPYU,LPXU,LPYU,IPXV,IPYV,LPXV,LPYV, 
     .    NU,NV,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)
          SPUV=SU*SV
          WRITE(6,1090) NU,EU,SU,NV,EV,SV,RUV  
          WRITE(6,1100) SRV1,SUV,RUV
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUV,J=1,NCOR)  
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NU,EU,SU,NV,EV,SV,RUV  
          WRITE(2,1040) SRV1,SUV,RUV
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUV,J=1,NCOR)  
          ENDIF
          WRITE(6,1160)
          SRV(1)=RUV
          DO 53 I=1,NCOR
          SRV(I+1)=COV(I)/SPUV 
53        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=2 TIPO=2 VARIABILI X,Y,U,V 
C         INIZIO CASO=2 TIPO=1 VARIABILI X,Y,W 
          CALL COVAP(X,Y,W,IPXW,IPYW,LPXW,LPYW,NW,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)  
          WRITE(6,1070) NW,EW,SW
          WRITE(6,1080) SRV1,S2W,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2W,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NW,EW,SW
          WRITE(2,1040) SRV1,S2W,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2W,J=1,NCOR)
          ENDIF
          WRITE(6,1160)
          SRV(1)=1.D0  
          DO 55 I=1,NCOR
          SRV(I+1)=COV(I)/S2W  
55        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=2 TIPO=1 VARIABILI X,Y,W
          DO 60 I=1,ND 
          U(I)=U(I)+EU 
          W(I)=W(I)+EW 
60        CONTINUE 
C         INIZIO CASO=2 TIPO=2 VARIABILI X,Y,U,W
          CALL RCORR2(U,W,LPXU,LPXW,EU,EW,SU,SW,S2U,S2W,SUW,RUW,ND)
          DO 70 I=1,ND 
          U(I)=U(I)-EU 
          W(I)=W(I)-EW 
70        CONTINUE 
          CALL CROSSP(X,Y,U,W,IPXU,IPYU,LPXU,LPYU,IPXW,IPYW,LPXW,LPYW, 
     .    NU,NW,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)
          SPUW=SU*SW
          WRITE(6,1110) NU,EU,SU,NW,EW,SW,RUW  
          WRITE(6,1120) SRV1,SUW,RUW
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUW,J=1,NCOR)  
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NU,EU,SU,NW,EW,SW,RUW  
          WRITE(2,1040) SRV1,SUW,RUW
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUW,J=1,NCOR)  
          ENDIF
          WRITE(6,1160)
          SRV(1)=RUW
          DO 75 I=1,NCOR
          SRV(I+1)=COV(I)/SPUW 
75        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=2 TIPO=2 VARIABILI X,Y,U,W 
          DO 80 I=1,ND 
          V(I)=V(I)+EV 
          W(I)=W(I)+EW 
80        CONTINUE 
C         INIZIO CASO=2 TIPO=2 VARIABILI X,Y,V,W
          CALL RCORR2(V,W,LPXV,LPXW,EV,EW,SV,SW,S2V,S2W,SVW,RVW,ND)
          DO 90 I=1,ND 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
90        CONTINUE 
          CALL CROSSP(X,Y,V,W,IPXV,IPYV,LPXV,LPYV,IPXW,IPYW,LPXW,LPYW, 
     .    NV,NW,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)
          SPVW=SV*SW
          WRITE(6,1130) NV,EV,SV,NW,EW,SW,RVW  
          WRITE(6,1140) SRV1,SVW,RVW
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPVW,J=1,NCOR)  
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NV,EV,SV,NW,EW,SW,RVW  
          WRITE(2,1040) SRV1,SVW,RVW
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPVW,J=1,NCOR)  
          ENDIF
          WRITE(6,1160)
          SRV(1)=RVW
          DO 95 I=1,NCOR
          SRV(I+1)=COV(I)/SPVW 
95        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=2 TIPO=2 VARIABILI X,Y,V,W 
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y)=',F10.3,5X,'S(Y)=',F10.3,//)
1060      FORMAT(//10X,'DELTA',7X,'GAMMA(Y)',5X,'RO(Y)',///,
     .    5X,F10.3,F15.6,F10.3)
1070      FORMAT(//5X,'N(Z) =',I5,5X,'M(Z)=',F10.3,5X,'S(Z)=',F10.3,//)
1080      FORMAT(//10X,'DELTA',7X,'GAMMA(Z)',5X,'RO(Z)',///,
     .    5X,F10.3,F15.6,F10.3)
1090      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//5X
     .    ,'N(Y) =',I5,5X,'M(Y)=',F10.3,5X,'S(Y)=',F10.3,5X,'R(X,Y)=', 
     .    F10.2,//)
1100      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Y)',3X,'RO(X,Y)',///5X, 
     .    F10.3,F15.6,F10.3)
1110      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//5X
     .    ,'N(Z) =',I5,5X,'M(Z)=',F10.3,5X,'S(Z)=',F10.3,5X,'R(X,Z)=', 
     .    F10.2,//)
1120      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Z)',3X,'RO(X,Z)',///5X, 
     .    F10.3,F15.6,F10.3)
1130      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y)=',F10.3,5X,'S(Y)=',F10.3,//5X
     .    ,'N(Z) =',I5,5X,'M(Z)=',F10.3,5X,'S(Z)=',F10.3,5X,'R(Y,Z)=', 
     .    F10.2,//)
1140      FORMAT(//10X,'DELTA',5X,'GAMMA(Y,Z)',3X,'RO(Y,Z)',///5X, 
     .    F10.3,F15.6,F10.3)
1150      FORMAT(I5,F15.6,F10.3,5X,I5,F15.6,F10.3,F10.2)
1160      FORMAT(/////)
          END  
          SUBROUTINE CORRT1(Z,W,IPZW,COV,SRV,NSRV,NCOR,RAGGIO,ALATO,
     .    NMAX,SCZ,SCW,IPTO,IAST,FORMA,IFILE)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20)
          CHARACTER*120 CDATO
          DIMENSION Z(1),W(1),IPZW(1),COV(1),SRV(1),NSRV(1),
     .    IPTO(1),IAST(1),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
C         INIZIO CASO=1 TIPO=1 VARIABILI Z,W
          SRV1=0.D0
          SRV2=1.D0
          REWIND 1 
          ND=0 
          NW=0 
10        READ(1,1000) (DATO(K),K=1,120)
C          DECODE(150,FORMA,DATO) NPT,TS,WS 
          READ(CDATO,FORMA) NPT,TS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          Z(ND)=TS*SCZ 
          W(ND)=WS*SCW 
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NW=NW+1  
          IPZW(NW)=ND  
          ENDIF
          GOTO 10  
20        CALL MEDVA2(W,IPZW,EW,S2W,SW,NW) 
          DO 30 I=1,ND 
          W(I)=W(I)-EW 
30        CONTINUE 
          CALL RIORDA(Z,IPZW,NW,ND,1)  
          CALL COVAT(Z,W,IPZW,NW,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX)  
          WRITE(6,1010) NW,EW,SW
          WRITE(6,1020) SRV1,S2W,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2W,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1050) NCOR
          WRITE(2,1050) NW,EW,SW
          WRITE(2,1040) SRV1,S2W,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2W,J=1,NCOR)
          ENDIF
          WRITE(6,1060)
          SRV(1)=1.D0  
          DO 40 I=1,NCOR
          SRV(I+1)=COV(I)/S2W  
40        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1060)
C         FINE CASO=1 TIPO=1 VARIABILI Z,W 
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(I5,F15.6,F10.3)
1060      FORMAT(/////)
          END  

          SUBROUTINE CORRT2(Z,U,V,IPZU,LPZU,IPZV,LPZV,COV,SRV,NSRV,NCOR,
     .    RAGGIO,ALATO,NMAX,SCZ,SCU,SCV,IPTO,IAST,FORMA,IFILE) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20)
          CHARACTER*120 CDATO
          DIMENSION Z(1),U(1),V(1),IPZU(1),LPZU(1),IPZV(1),LPZV(1),
     .    COV(1),SRV(1),NSRV(1),IPTO(1),IAST(1),DATO(120)  
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          SRV1=0.D0
          SRV2=1.D0
          ND=0 
          NU=0 
          NV=0 
          REWIND 1 
10        READ(1,1000) (DATO(K),K=1,120)
C          DECODE(150,FORMA,DATO) NPT,TS,US,VS  
          READ(CDATO,FORMA) NPT,TS,US,VS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          Z(ND)=TS*SCZ 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPZU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPZV(NV)=ND  
          ENDIF
          GOTO 10  
20        CALL RIORDA(Z,IPZU,NU,ND,1)  
          CALL RIORDA(Z,IPZV,NV,ND,1)  
          CALL PUNTA(IPZU,LPZU,0,1,NU,ND)  
          CALL PUNTA(IPZV,LPZV,0,1,NV,ND)  
          CALL MEDVA2(U,IPZU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPZV,EV,S2V,SV,NV) 
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
30        CONTINUE 
C         INIZIO CASO=1 TIPO=1 VARIABILI Z,U
          CALL COVAT(Z,U,IPZU,NU,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX)  
          WRITE(6,1010) NU,EU,SU
          WRITE(6,1020) SRV1,S2U,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2U,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1090) NCOR
          WRITE(2,1090) NU,EU,SU
          WRITE(2,1040) SRV1,S2U,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2U,J=1,NCOR)
          ENDIF
          WRITE(6,1100)
          SRV(1)=1.D0  
          DO 35 I=1,NCOR
          SRV(I+1)=COV(I)/S2U  
35        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1100)
C         FINE CASO=1 TIPO=1 VARIABILI Z,U 
C         INIZIO CASO=1 TIPO=1 VARIABILI Z,V
          CALL COVAT(Z,V,IPZV,NV,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX)  
          WRITE(6,1050) NV,EV,SV
          WRITE(6,1060) SRV1,S2V,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2V,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1090) NCOR
          WRITE(2,1090) NV,EV,SV
          WRITE(2,1040) SRV1,S2V,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2V,J=1,NCOR)
          ENDIF
          WRITE(6,1100)
          SRV(1)=1.D0  
          DO 37 I=1,NCOR
          SRV(I+1)=COV(I)/S2V  
37        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1100)
C         FINE CASO=1 TIPO=1 VARIABILI Z,V 
          DO 40 I=1,ND 
          U(I)=U(I)+EU 
          V(I)=V(I)+EV 
40        CONTINUE 
C         INIZIO CASO=1 TIPO=2 VARIABILI Z,U,V 
          CALL RCORR2(U,V,LPZU,LPZV,EU,EV,SU,SV,S2U,S2V,SUV,RUV,ND)
          DO 50 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
50        CONTINUE 
          CALL CROSST(Z,U,V,IPZU,LPZU,IPZV,LPZV,NU,NV, 
     .    NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX) 
          SPUV=SU*SV
          WRITE(6,1070) NU,EU,SU,NV,EV,SV,RUV  
          WRITE(6,1080) SRV1,SUV,RUV
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUV,J=1,NCOR)  
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1090) NCOR
          WRITE(2,1090) NU,EU,SU,NV,EV,SV,RUV  
          WRITE(2,1040) SRV1,SUV,RUV
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUV,J=1,NCOR)  
          ENDIF
          WRITE(6,1100)
          SRV(1)=RUV
          DO 65 I=1,NCOR
          SRV(I+1)=COV(I)/SPUV 
65        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1100)
C         FINE CASO=1 TIPO=2 VARIABILI Z,U,V
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y)=',F10.3,5X,'S(Y)=',F10.3,//)
1060      FORMAT(//10X,'DELTA',7X,'GAMMA(Y)',5X,'RO(Y)',///,
     .    5X,F10.3,F15.6,F10.3)
1070      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//5X
     .    ,'N(Y) =',I5,5X,'M(Y)=',F10.3,5X,'S(Y)=',F10.3,5X,'R(X,Y)=', 
     .    F10.2,//)
1080      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Y)',3X,'RO(X,Y)',///5X, 
     .    F10.3,F15.6,F10.3)
1090      FORMAT(I5,F15.6,F10.3,5X,I5,F15.6,F10.3,F10.2)
1100      FORMAT(/////)
          END  
          SUBROUTINE CORRT3(Z,U,V,W,IPZU,LPZU,IPZV,LPZV,IPZW,LPZW, 
     .    COV,SRV,NSRV,NCOR,RAGGIO,ALATO,NMAX, 
     .    SCZ,SCU,SCV,SCW,IPTO,IAST,FORMA,IFILE)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20)
          CHARACTER*120 CDATO
          DIMENSION Z(1),U(1),V(1),W(1),IPZU(1),LPZU(1),
     .    IPZV(1),LPZV(1),IPZW(1),LPZW(1), 
     .    COV(1),SRV(1),NSRV(1),IPTO(1),IAST(1),DATO(120)  
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          SRV1=0.D0
          SRV2=1.D0
          ND=0 
          NU=0 
          NV=0 
          NW=0 
          REWIND 1 
10        READ(1,1000) (DATO(K),K=1,120)
C          DECODE(150,FORMA,DATO) NPT,TS,US,VS,WS
          READ(CDATO,FORMA) NPT,TS,US,VS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          Z(ND)=TS*SCZ 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          W(ND)=WS*SCW 
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPZU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPZV(NV)=ND  
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(3)).NE.AST) THEN
          NW=NW+1  
          IPZW(NW)=ND  
          ENDIF
          GOTO 10  
20        CALL RIORDA(Z,IPZU,NU,ND,1)  
          CALL RIORDA(Z,IPZV,NV,ND,1)  
          CALL RIORDA(Z,IPZW,NW,ND,1)  
          CALL PUNTA(IPZU,LPZU,0,1,NU,ND)  
          CALL PUNTA(IPZV,LPZV,0,1,NV,ND)  
          CALL PUNTA(IPZW,LPZW,0,1,NW,ND)  
          CALL MEDVA2(U,IPZU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPZV,EV,S2V,SV,NU) 
          CALL MEDVA2(W,IPZW,EW,S2W,SW,NU) 
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
30        CONTINUE 
C         INIZIO CASO=1 TIPO=1 VARIABILI Z,U
          CALL COVAT(Z,U,IPZU,NU,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX)  
          WRITE(6,1010) NU,EU,SU
          WRITE(6,1020) SRV1,S2U,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2U,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NU,EU,SU
          WRITE(2,1040) SRV1,S2U,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2U,J=1,NCOR)
          ENDIF
          WRITE(6,1160)
          SRV(1)=1.D0  
          DO 35 I=1,NCOR
          SRV(I+1)=COV(I)/S2U  
35        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=1 TIPO=1 VARIABILI Z,U 
C         INIZIO CASO=1 TIPO=1 VARIABILI Z,V
          CALL COVAT(Z,V,IPZV,NV,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX)  
          WRITE(6,1050) NV,EV,SV
          WRITE(6,1060) SRV1,S2V,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2V,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NV,EV,SV
          WRITE(2,1040) SRV1,S2V,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2V,J=1,NCOR)
          ENDIF
          WRITE(6,1160)
          SRV(1)=1.D0  
          DO 37 I=1,NCOR
          SRV(I+1)=COV(I)/S2V  
37        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=1 TIPO=1 VARIABILI Z,V 
          DO 40 I=1,ND 
          U(I)=U(I)+EU 
          V(I)=V(I)+EV 
40        CONTINUE 
C         INIZIO CASO=1 TIPO=2 VARIABILI Z,U,V 
          CALL RCORR2(U,V,LPZU,LPZV,EU,EV,SU,SV,S2U,S2V,SUV,RUV,ND)
          DO 50 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
50        CONTINUE 
          CALL CROSST(Z,U,V,IPZU,LPZU,IPZV,LPZV,NU,NV, 
     .    NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX) 
          SPUV=SU*SV
          WRITE(6,1090) NU,EU,SU,NV,EV,SV,RUV  
          WRITE(6,1100) SRV1,SUV,RUV
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUV,J=1,NCOR)  
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NU,EU,SU,NV,EV,SV,RUV  
          WRITE(2,1040) SRV1,SUV,RUV
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUV,J=1,NCOR)  
          ENDIF
          WRITE(6,1160)
          SRV(1)=RUV
          DO 53 I=1,NCOR
          SRV(I+1)=COV(I)/SPUV 
53        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=1 TIPO=2 VARIABILI Z,U,V
C         INIZIO CASO=1 TIPO=1 VARIABILI Z,W
          CALL COVAT(Z,W,IPZW,NW,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX)  
          WRITE(6,1070) NW,EW,SW
          WRITE(6,1080) SRV1,S2W,SRV2  
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2W,J=1,NCOR)
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NW,EW,SW
          WRITE(2,1040) SRV1,S2W,SRV2  
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/S2W,J=1,NCOR)
          ENDIF
          WRITE(6,1160)
          SRV(1)=1.D0  
          DO 55 I=1,NCOR
          SRV(I+1)=COV(I)/S2W  
55        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=1 TIPO=1 VARIABILI Z,W 
          DO 60 I=1,ND 
          U(I)=U(I)+EU 
          W(I)=W(I)+EW 
60        CONTINUE 
C         INIZIO CASO=1 TIPO=2 VARIABILI Z,U,W 
          CALL RCORR2(U,W,LPZU,LPZW,EU,EW,SU,SW,S2U,S2W,SUW,RUW,ND)
          DO 70 I=1,ND 
          U(I)=U(I)-EU 
          W(I)=W(I)-EW 
70        CONTINUE 
          CALL CROSST(Z,U,W,IPZU,LPZU,IPZW,LPZW,NU,NW, 
     .    NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX) 
          SPUW=SU*SW
          WRITE(6,1110) NU,EU,SU,NW,EW,SW,RUW  
          WRITE(6,1120) SRV1,SUW,RUW
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUW,J=1,NCOR)  
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NU,EU,SU,NW,EW,SW,RUW  
          WRITE(2,1040) SRV1,SUW,RUW
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPUW,J=1,NCOR)  
          ENDIF
          WRITE(6,1160)
          SRV(1)=RUW
          DO 75 I=1,NCOR
          SRV(I+1)=COV(I)/SPUW 
75        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=1 TIPO=2 VARIABILI Z,U,W
          DO 80 I=1,ND 
          V(I)=V(I)+EV 
          W(I)=W(I)+EW 
80        CONTINUE 
C         INIZIO CASO=1 TIPO=2 VARIABILI Z,V,W 
          CALL RCORR2(V,W,LPZV,LPZW,EV,EW,SV,SW,S2V,S2W,SVW,RVW,ND)
          DO 90 I=1,ND 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
90        CONTINUE 
          CALL CROSST(Z,V,W,IPZV,LPZV,IPZW,LPZW,NV,NW, 
     .    NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX) 
          SPVW=SV*SW
          WRITE(6,1130) NV,EV,SV,NW,EW,SW,RVW  
          WRITE(6,1140) SRV1,SVW,RVW
          WRITE(6,1030) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPVW,J=1,NCOR)  
          IF(IFILE.EQ.1) THEN  
          WRITE(2,1150) NCOR
          WRITE(2,1150) NV,EV,SV,NW,EW,SW,RVW  
          WRITE(2,1040) SRV1,SVW,RVW
          WRITE(2,1040) (RAGGIO*(J-.5D0),COV(J),COV(J)/SPVW,J=1,NCOR)  
          ENDIF
          WRITE(6,1160)
          SRV(1)=RVW
          DO 95 I=1,NCOR
          SRV(I+1)=COV(I)/SPVW 
95        CONTINUE 
          NN=NCOR+1
          CALL GRACO(SRV,SRV,RAGGIO,NN)
          WRITE(6,1160)
C         FINE CASO=1 TIPO=2 VARIABILI Z,V,W
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y)=',F10.3,5X,'S(Y)=',F10.3,//)
1060      FORMAT(//10X,'DELTA',7X,'GAMMA(Y)',5X,'RO(Y)',///,
     .    5X,F10.3,F15.6,F10.3)
1070      FORMAT(//5X,'N(Z) =',I5,5X,'M(Z)=',F10.3,5X,'S(Z)=',F10.3,//)
1080      FORMAT(//10X,'DELTA',7X,'GAMMA(Z)',5X,'RO(Z)',///,
     .    5X,F10.3,F15.6,F10.3)
1090      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,// 
     .    5X,'N(Y) =',I5,5X,'M(Y)=',F10.3,5X,'S(Y)=',F10.3,5X,'R(X,Y)=',
     .    F10.2,//)
1100      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Y)',3X,'RO(X,Y)',///5X, 
     .    F10.3,F15.6,F10.3)
1110      FORMAT(//5X,'N(X) =',I5,5X,'M(X)=',F10.3,5X,'S(X)=',F10.3,//5X
     .    ,'N(Z) =',I5,5X,'M(Z)=',F10.3,5X,'S(Z)=',F10.3,5X,'R(X,Z)=', 
     .    F10.2,//)
1120      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Z)',3X,'RO(X,Z)',///5X, 
     .    F10.3,F15.6,F10.3)
1130      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y)=',F10.3,5X,'S(Y)=',F10.3,//,
     .    5X,'N(Z) =',I5,5X,'M(Z)=',F10.3,5X,'S(Z)=',F10.3,5X,'R(Y,Z)=',
     .    F10.2,//)
1140      FORMAT(//10X,'DELTA',5X,'GAMMA(Y,Z)',3X,'RO(Y,Z)',///5X, 
     .    F10.3,F15.6,F10.3)
1150      FORMAT(I5,F15.6,F10.3,5X,I5,F15.6,F10.3,F10.2)
1160      FORMAT(/////)
          END  
          SUBROUTINE COVAP(X,Y,U,IPX,IPY,LPX,LPY,N,NCOR,RGGO,ALATO,COV,
     .    SRV,NSRV,IWD,NMAX,NWDMAX)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),U(1),IPX(1),IPY(1),LPX(1),LPY(1),COV(1), 
     .    IWD(1),SRV(1),NSRV(1)
          DO 10 K=1,NCOR
          COV(K)=0.D0  
10        CONTINUE 
          DO 50 I=1,N  
          LL=IPX(I)
          DO 20 K=1,NCOR
          SRV(K)=0.D0  
          NSRV(K)=0
20        CONTINUE 
          XC=X(LL) 
          YC=Y(LL) 
          CALL AWINDW(XC,YC,ALATO,ALATO,X,Y,IPX,IPY,LPX,LPY,IWD,NMAX,  
     .    NWDMAX,N,NWD)
          DO 30 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+U(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          SRV(1)=SRV(1)-U(LL)  
          NSRV(1)=NSRV(1)-1
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+U(LL)*SRV(K)/NK
40        CONTINUE 
50        CONTINUE 
          DO 60 K=1,NCOR
          COV(K)=COV(K)/N  
60        CONTINUE 
          RETURN
          END  
          SUBROUTINE COVAT(Z,U,IPZ,N,NCOR,RGGO,ALATO,COV,SRV,NSRV, 
     .    NMAX)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION Z(1),U(1),IPZ(1),COV(1),SRV(1),NSRV(1) 
          ALATO2=ALATO/2.D0
          DO 10 K=1,NCOR
          COV(K)=0.D0  
10        CONTINUE 
          DO 50 I=1,N  
          LL=IPZ(I)
          DO 20 K=1,NCOR
          SRV(K)=0.D0  
          NSRV(K)=0
20        CONTINUE 
          ZC=Z(LL) 
          ZMI=ZC-ALATO2
          ZMA=ZC+ALATO2
          CALL ABINRY(ZMI,Z,IPZ,NMAX,LNI,LNP,LPI,LPS,N)
          LZI=LPS  
          CALL ABINRY(ZMA,Z,IPZ,NMAX,LNI,LNP,LPI,LPS,N)
          LZS=LPI  
          DO 30 J=LZI,LZS  
          L=IPZ(J) 
          DIJ=DABS(ZC-Z(L))
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+U(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          SRV(1)=SRV(1)-U(LL)  
          NSRV(1)=NSRV(1)-1
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+U(LL)*SRV(K)/NK
40        CONTINUE 
50        CONTINUE 
          DO 60 K=1,NCOR
          COV(K)=COV(K)/N  
60        CONTINUE 
          RETURN
          END  
          SUBROUTINE CROSSP(X,Y,V,W,IPX,IPY,LPX,LPY,IPX2,IPY2,LPX2,LPY2,
     .    N,M,NCOR,RGGO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),V(1),W(1),IPX(1),IPY(1),LPX(1),LPY(1),
     .    IPX2(1),IPY2(1),LPX2(1),LPY2(1),COV(1),SRV(1),NSRV(1),IWD(1) 
          DO 10 K=1,NCOR
          COV(K)=0.D0  
10        CONTINUE 
          DO 50 I=1,N  
          DO 20 K=1,NCOR
          SRV(K)=0.D0  
          NSRV(K)=0
20        CONTINUE 
          LL=IPX(I)
          XC=X(LL) 
          YC=Y(LL) 
          CALL AWINDW(XC,YC,ALATO,ALATO,X,Y,IPX2,IPY2,LPX2,LPY2,IWD,
     .    NMAX,NWDMAX,M,NWD)
          DO 30 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+W(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          IF(LPX2(LL).NE.0) THEN
          SRV(1)=SRV(1)-W(LL)  
          NSRV(1)=NSRV(1)-1
          ENDIF
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+V(LL)*SRV(K)/NK
40        CONTINUE 
50        CONTINUE 
          DO 90 I=1,M  
          DO 60 K=1,NCOR
          SRV(K)=0.D0  
          NSRV(K)=0
60        CONTINUE 
          LL=IPX2(I)
          XC=X(LL) 
          YC=Y(LL) 
          CALL AWINDW(XC,YC,ALATO,ALATO,X,Y,IPX,IPY,LPX,LPY,IWD,
     .    NMAX,NWDMAX,N,NWD)
          DO 70 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+V(L)
          NSRV(K)=NSRV(K)+1
70        CONTINUE 
          IF(LPX(LL).NE.0) THEN
          SRV(1)=SRV(1)-V(LL)  
          NSRV(1)=NSRV(1)-1
          ENDIF
          DO 80 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+W(LL)*SRV(K)/NK
80        CONTINUE 
90        CONTINUE 
          DO 100 K=1,NCOR  
          COV(K)=COV(K)/(N+M)  
100       CONTINUE 
          RETURN
          END  
          SUBROUTINE CROSST(Z,V,W,IPZ,LPZ,IPZ2,LPZ2,N,M,NCOR,  
     .    RGGO,ALATO,COV,SRV,NSRV,NMAX)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION Z(1),V(1),W(1),IPZ(1),IPZ2(1),LPZ(1),LPZ2(1),COV(1),
     .    SRV(1),NSRV(1)
          ALATO2=ALATO/2.D0
          DO 10 K=1,NCOR
          COV(K)=0.D0  
10        CONTINUE 
          DO 50 I=1,N  
          DO 20 K=1,NCOR
          SRV(K)=0.D0  
          NSRV(K)=0
20        CONTINUE 
          LL=IPZ(I)
          ZC=Z(LL) 
          ZMI=ZC-ALATO2
          ZMA=ZC+ALATO2
          CALL ABINRY(ZMI,Z,IPZ2,NMAX,LNI,LNP,LPI,LPS,M)
          LZI=LPS  
          CALL ABINRY(ZMA,Z,IPZ2,NMAX,LNI,LNP,LPI,LPS,M)
          LZS=LPI  
          DO 30 J=LZI,LZS  
          L=IPZ2(J)
          DIJ=DABS(ZC-Z(L))
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+W(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          IF(LPZ2(LL).NE.0) THEN
          SRV(1)=SRV(1)-W(LL)  
          NSRV(1)=NSRV(1)-1
          ENDIF
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+V(LL)*SRV(K)/NK
40        CONTINUE 
50        CONTINUE 
          DO 90 I=1,M  
          DO 60 K=1,NCOR
          SRV(K)=0.D0  
          NSRV(K)=0
60        CONTINUE 
          LL=IPZ2(I)
          ZC=Z(LL) 
          ZMI=ZC-ALATO2
          ZMA=ZC+ALATO2
          CALL ABINRY(ZMI,Z,IPZ,NMAX,LNI,LNP,LPI,LPS,N)
          LZI=LPS  
          CALL ABINRY(ZMA,Z,IPZ,NMAX,LNI,LNP,LPI,LPS,N)
          LZS=LPI  
          DO 70 J=LZI,LZS  
          L=IPZ(J) 
          DIJ=DABS(ZC-Z(L))
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+V(L)
          NSRV(K)=NSRV(K)+1
70        CONTINUE 
          IF(LPZ(LL).NE.0) THEN
          SRV(1)=SRV(1)-V(LL)  
          NSRV(1)=NSRV(1)-1
          ENDIF
          DO 80 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+W(LL)*SRV(K)/NK
80        CONTINUE 
90        CONTINUE 
          DO 100 K=1,NCOR  
          COV(K)=COV(K)/(N+M)  
100       CONTINUE 
          RETURN
          END  
          FUNCTION GASS1(DELTA,XCOVT,II,JJ,L,K,DA,DB,DC,SN,NCASO)  
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          DIMENSION DA(1),DB(1),DC(1),SN(1),NCASO(1),MATI(3,3,49)  
          DATA ((MATI(I,J,1),J=1,3),I=1,3)/1,0,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,2),J=1,3),I=1,3)/3,0,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,3),J=1,3),I=1,3)/1,3,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,4),J=1,3),I=1,3)/5,0,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,5),J=1,3),I=1,3)/1,5,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,6),J=1,3),I=1,3)/3,5,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,7),J=1,3),I=1,3)/1,3,5,0,0,0,0,0,0/  
          DATA ((MATI(I,J,8),J=1,3),I=1,3)/3,0,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,9),J=1,3),I=1,3)/2,0,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,10),J=1,3),I=1,3)/3,2,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,11),J=1,3),I=1,3)/6,0,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,12),J=1,3),I=1,3)/3,6,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,13),J=1,3),I=1,3)/2,6,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,14),J=1,3),I=1,3)/3,2,6,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,15),J=1,3),I=1,3)/1,0,0,3,0,0,0,0,0/ 
          DATA ((MATI(I,J,16),J=1,3),I=1,3)/3,0,0,2,0,0,0,0,0/ 
          DATA ((MATI(I,J,17),J=1,3),I=1,3)/1,3,0,3,2,0,0,0,0/ 
          DATA ((MATI(I,J,18),J=1,3),I=1,3)/5,0,0,6,0,0,0,0,0/ 
          DATA ((MATI(I,J,19),J=1,3),I=1,3)/1,5,0,3,6,0,0,0,0/ 
          DATA ((MATI(I,J,20),J=1,3),I=1,3)/3,5,0,2,6,0,0,0,0/ 
          DATA ((MATI(I,J,21),J=1,3),I=1,3)/1,3,5,3,2,6,0,0,0/ 
          DATA ((MATI(I,J,22),J=1,3),I=1,3)/5,0,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,23),J=1,3),I=1,3)/6,0,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,24),J=1,3),I=1,3)/5,6,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,25),J=1,3),I=1,3)/4,0,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,26),J=1,3),I=1,3)/5,4,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,27),J=1,3),I=1,3)/6,4,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,28),J=1,3),I=1,3)/5,6,4,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,29),J=1,3),I=1,3)/1,0,0,5,0,0,0,0,0/ 
          DATA ((MATI(I,J,30),J=1,3),I=1,3)/3,0,0,6,0,0,0,0,0/ 
          DATA ((MATI(I,J,31),J=1,3),I=1,3)/1,3,0,5,6,0,0,0,0/ 
          DATA ((MATI(I,J,32),J=1,3),I=1,3)/5,0,0,4,0,0,0,0,0/ 
          DATA ((MATI(I,J,33),J=1,3),I=1,3)/1,5,0,5,4,0,0,0,0/ 
          DATA ((MATI(I,J,34),J=1,3),I=1,3)/3,5,0,6,4,0,0,0,0/ 
          DATA ((MATI(I,J,35),J=1,3),I=1,3)/1,3,5,5,6,4,0,0,0/ 
          DATA ((MATI(I,J,36),J=1,3),I=1,3)/3,0,0,5,0,0,0,0,0/ 
          DATA ((MATI(I,J,37),J=1,3),I=1,3)/2,0,0,6,0,0,0,0,0/ 
          DATA ((MATI(I,J,38),J=1,3),I=1,3)/3,2,0,5,6,0,0,0,0/ 
          DATA ((MATI(I,J,39),J=1,3),I=1,3)/6,0,0,4,0,0,0,0,0/ 
          DATA ((MATI(I,J,40),J=1,3),I=1,3)/3,6,0,5,4,0,0,0,0/ 
          DATA ((MATI(I,J,41),J=1,3),I=1,3)/2,6,0,6,4,0,0,0,0/ 
          DATA ((MATI(I,J,42),J=1,3),I=1,3)/3,2,6,5,6,4,0,0,0/ 
          DATA ((MATI(I,J,43),J=1,3),I=1,3)/1,0,0,3,0,0,5,0,0/ 
          DATA ((MATI(I,J,44),J=1,3),I=1,3)/3,0,0,2,0,0,6,0,0/ 
          DATA ((MATI(I,J,45),J=1,3),I=1,3)/1,3,0,3,2,0,5,6,0/ 
          DATA ((MATI(I,J,46),J=1,3),I=1,3)/5,0,0,6,0,0,4,0,0/ 
          DATA ((MATI(I,J,47),J=1,3),I=1,3)/1,5,0,3,6,0,5,4,0/ 
          DATA ((MATI(I,J,48),J=1,3),I=1,3)/3,5,0,2,6,0,6,4,0/ 
          DATA ((MATI(I,J,49),J=1,3),I=1,3)/1,3,5,3,2,6,5,6,4/ 
          IND=(II-1)*7+JJ  
          IC=MATI(L,K,IND) 
          DDA=DA(IC)
          DDB=DB(IC)
          DDC=DC(IC)
          SSN=SN(IC)
          ICC=NCASO(IC)
          IF(DELTA.GT.0.D0) THEN
          GO TO(10,20,30,40,50,60,70,80,90),ICC
10        SRV=DDA*DEXP(-DDB*DELTA) 
          GASS1=SRV*SCAL1(DELTA/XCOVT) 
          RETURN
20        SRV=DDA*DEXP(-DDB*DELTA**2)  
          GASS1=SRV*SCAL1(DELTA/XCOVT) 
          RETURN
30        SRV=DDA*(1.D0-DDC*DELTA**2)*DEXP(-DDB*DELTA) 
          GASS1=SRV*SCAL1(DELTA/XCOVT) 
          RETURN
40        SRV=DDA*(1.D0-DDC*DELTA**2)*DEXP(-DDB*DELTA**2)  
          GASS1=SRV*SCAL1(DELTA/XCOVT) 
          RETURN
50        SRV=DDA*DEXP(-DDB*DELTA)*DSIN(DDC*DELTA)/(DDC*DELTA) 
          GASS1=SRV*SCAL1(DELTA/XCOVT) 
          RETURN
60        SRV=DDA*DEXP(-DDB*DELTA**2)*DSIN(DDC*DELTA)/(DDC*DELTA)  
          GASS1=SRV*SCAL1(DELTA/XCOVT) 
          RETURN
70        SRV=DDA*DEXP(-DDB*DELTA)*DCOS(DDC*DELTA) 
          GASS1=SRV*SCAL1(DELTA/XCOVT) 
          RETURN
80        SRV=DDA*DEXP(-DDB*DELTA**2)*DCOS(DDC*DELTA)  
          GASS1=SRV*SCAL1(DELTA/XCOVT) 
          RETURN
90        SRV=DDA*DEXP(-DDB*DELTA)*(DCOS(DDC*DELTA)+
     .    (DDB/DDC)*DSIN(DDC*DELTA))
          GASS1=SRV*SCAL1(DELTA/XCOVT) 
          RETURN
          ELSE 
          GASS1=DDA+SSN
          ENDIF
          RETURN
          END  
          FUNCTION GASS2(DELTA,XCOVT,II,JJ,L,K,DA,DB,DC,SN,NCASO)  
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          DIMENSION DA(1),DB(1),DC(1),SN(1),NCASO(1),MATI(3,3,49)  
          DATA ((MATI(I,J,1),J=1,3),I=1,3)/1,0,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,2),J=1,3),I=1,3)/3,0,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,3),J=1,3),I=1,3)/1,3,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,4),J=1,3),I=1,3)/5,0,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,5),J=1,3),I=1,3)/1,5,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,6),J=1,3),I=1,3)/3,5,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,7),J=1,3),I=1,3)/1,3,5,0,0,0,0,0,0/  
          DATA ((MATI(I,J,8),J=1,3),I=1,3)/3,0,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,9),J=1,3),I=1,3)/2,0,0,0,0,0,0,0,0/  
          DATA ((MATI(I,J,10),J=1,3),I=1,3)/3,2,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,11),J=1,3),I=1,3)/6,0,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,12),J=1,3),I=1,3)/3,6,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,13),J=1,3),I=1,3)/2,6,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,14),J=1,3),I=1,3)/3,2,6,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,15),J=1,3),I=1,3)/1,0,0,3,0,0,0,0,0/ 
          DATA ((MATI(I,J,16),J=1,3),I=1,3)/3,0,0,2,0,0,0,0,0/ 
          DATA ((MATI(I,J,17),J=1,3),I=1,3)/1,3,0,3,2,0,0,0,0/ 
          DATA ((MATI(I,J,18),J=1,3),I=1,3)/5,0,0,6,0,0,0,0,0/ 
          DATA ((MATI(I,J,19),J=1,3),I=1,3)/1,5,0,3,6,0,0,0,0/ 
          DATA ((MATI(I,J,20),J=1,3),I=1,3)/3,5,0,2,6,0,0,0,0/ 
          DATA ((MATI(I,J,21),J=1,3),I=1,3)/1,3,5,3,2,6,0,0,0/ 
          DATA ((MATI(I,J,22),J=1,3),I=1,3)/5,0,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,23),J=1,3),I=1,3)/6,0,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,24),J=1,3),I=1,3)/5,6,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,25),J=1,3),I=1,3)/4,0,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,26),J=1,3),I=1,3)/5,4,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,27),J=1,3),I=1,3)/6,4,0,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,28),J=1,3),I=1,3)/5,6,4,0,0,0,0,0,0/ 
          DATA ((MATI(I,J,29),J=1,3),I=1,3)/1,0,0,5,0,0,0,0,0/ 
          DATA ((MATI(I,J,30),J=1,3),I=1,3)/3,0,0,6,0,0,0,0,0/ 
          DATA ((MATI(I,J,31),J=1,3),I=1,3)/1,3,0,5,6,0,0,0,0/ 
          DATA ((MATI(I,J,32),J=1,3),I=1,3)/5,0,0,4,0,0,0,0,0/ 
          DATA ((MATI(I,J,33),J=1,3),I=1,3)/1,5,0,5,4,0,0,0,0/ 
          DATA ((MATI(I,J,34),J=1,3),I=1,3)/3,5,0,6,4,0,0,0,0/ 
          DATA ((MATI(I,J,35),J=1,3),I=1,3)/1,3,5,5,6,4,0,0,0/ 
          DATA ((MATI(I,J,36),J=1,3),I=1,3)/3,0,0,5,0,0,0,0,0/ 
          DATA ((MATI(I,J,37),J=1,3),I=1,3)/2,0,0,6,0,0,0,0,0/ 
          DATA ((MATI(I,J,38),J=1,3),I=1,3)/3,2,0,5,6,0,0,0,0/ 
          DATA ((MATI(I,J,39),J=1,3),I=1,3)/6,0,0,4,0,0,0,0,0/ 
          DATA ((MATI(I,J,40),J=1,3),I=1,3)/3,6,0,5,4,0,0,0,0/ 
          DATA ((MATI(I,J,41),J=1,3),I=1,3)/2,6,0,6,4,0,0,0,0/ 
          DATA ((MATI(I,J,42),J=1,3),I=1,3)/3,2,6,5,6,4,0,0,0/ 
          DATA ((MATI(I,J,43),J=1,3),I=1,3)/1,0,0,3,0,0,5,0,0/ 
          DATA ((MATI(I,J,44),J=1,3),I=1,3)/3,0,0,2,0,0,6,0,0/ 
          DATA ((MATI(I,J,45),J=1,3),I=1,3)/1,3,0,3,2,0,5,6,0/ 
          DATA ((MATI(I,J,46),J=1,3),I=1,3)/5,0,0,6,0,0,4,0,0/ 
          DATA ((MATI(I,J,47),J=1,3),I=1,3)/1,5,0,3,6,0,5,4,0/ 
          DATA ((MATI(I,J,48),J=1,3),I=1,3)/3,5,0,2,6,0,6,4,0/ 
          DATA ((MATI(I,J,49),J=1,3),I=1,3)/1,3,5,3,2,6,5,6,4/ 
          IND=(II-1)*7+JJ  
          IC=MATI(L,K,IND) 
          DDA=DA(IC)
          DDB=DB(IC)
          DDC=DC(IC)
          SSN=SN(IC)
          ICC=NCASO(IC)
          IF(DELTA.GT.0.D0) THEN
          GO TO(10,20,30,40,45,45,45,45,50,60,70,80),ICC
10        SRV=DDA*DEXP(-DDB*DELTA) 
          GASS2=SRV*SCAL2(DELTA/XCOVT) 
          RETURN
20        SRV=DDA*DEXP(-DDB*DELTA**2)  
          GASS2=SRV*SCAL2(DELTA/XCOVT) 
          RETURN
30        SRV=DDA*(1.D0-DDC*DELTA**2)*DEXP(-DDB*DELTA) 
          GASS2=SRV*SCAL2(DELTA/XCOVT) 
          RETURN
40        SRV=DDA*(1.D0-DDC*DELTA**2)*DEXP(-DDB*DELTA**2)  
          GASS2=SRV*SCAL2(DELTA/XCOVT) 
          RETURN
45        GASS2=0.D0
          RETURN
50        SRV=DDA*DEXP(-DDB*DELTA)*DFJ1(DDC*DELTA)/(DDC*DELTA) 
          GASS2=SRV*SCAL2(DELTA/XCOVT) 
          RETURN
60        SRV=DDA*DEXP(-DDB*DELTA**2)*DFJ1(DDC*DELTA)/(DDC*DELTA)  
          GASS2=SRV*SCAL2(DELTA/XCOVT) 
          RETURN
70        SRV=DDA*DEXP(-DDB*DELTA)*DFJ0(DDC*DELTA) 
          GASS2=SRV*SCAL2(DELTA/XCOVT) 
          RETURN
80        SRV=DDA*DEXP(-DDB*DELTA**2)*DFJ0(DDC*DELTA)  
          GASS2=SRV*SCAL2(DELTA/XCOVT) 
          RETURN
          ELSE 
          GASS2=DDA+SSN
          ENDIF
          RETURN
          END  
          SUBROUTINE GRACO(COV,TCOV,STP,NCOV)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 LB(101),PUN,AST,PIU,BL
          DIMENSION A(11),COV(1),TCOV(1)
          DATA PUN,AST,PIU,BL/'.','*','+',' '/ 
          DO 10 I=1,11 
          A(I)=-1.D0+(I-1)*.2D0
10        CONTINUE 
          WRITE(6,100) (A(I),I=1,11)
          DO 20 I=1,101
          LB(I)=PUN
20        CONTINUE 
          DO 30 I=1,101,10 
          LB(I)='I'
30        CONTINUE 
          WRITE(6,110)(LB(K),K=1,101)  
          X0=STP/2.D0  
          X=0.D0
          DO 40 I=1,NCOV
          DO 35 K=1,101
          LB(K)=BL 
35        CONTINUE 
          LB(51)=PUN
          II=IDNINT(COV(I)*50.D0)+51
          JJ=IDNINT(TCOV(I)*50.D0)+51  
          LB(JJ)=PIU
          LB(II)=AST
          WRITE(6,120) X,(LB(J),J=1,101)
          X=X0+(I-1)*STP
40        CONTINUE 
          RETURN
100       FORMAT(9X,11(F4.1,6X))
110       FORMAT(11X,101A1)
120       FORMAT(1X,F10.3,101A1)
          END  
      FUNCTION IDIF2(COV)  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      DIMENSION COV(1) 
      D=COV(2)-2.D0*COV(3)+COV(4)  
      IF(D.GE.0.D0) THEN
      IDIF2=1  
      ELSE 
      IDIF2=2  
      ENDIF
      RETURN
      END  
          SUBROUTINE LEGEND(ICASO) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          WRITE(6,1000)
          IF(ICASO.EQ.1) THEN  
          WRITE(6,1010)
          WRITE(6,1020)
          WRITE(6,1030)
          WRITE(6,1040)
          WRITE(6,1050)
          WRITE(6,1060)
          WRITE(6,1070)
          WRITE(6,1080)
          ENDIF
          IF(ICASO.EQ.2) THEN  
          WRITE(6,1010)
          WRITE(6,1020)
          WRITE(6,1090)
          WRITE(6,1040)
          WRITE(6,1100)
          WRITE(6,1110)
          WRITE(6,1120)
          WRITE(6,1130)
          ENDIF
          RETURN
1000      FORMAT(///5X,'CASO:',15X,'FUNZIONE TEORICA:'/)
1010      FORMAT(9X,'1',15X,'CASO ESPONENZIALE'
     ./25X,'A*EXP(-B*ABS(X))'/)
1020      FORMAT(9X,'2',15X,'CASO NORMALE' 
     ./25X,'A*EXP(-B*X**2)'/)  
1030      FORMAT(9X,'3',15X,'CASO ESPONENZIALE PARABOLA'
     ./25X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(2*D)'/)
1040      FORMAT(9X,'4',15X,'CASO NORMALE PARABOLA'
     ./25X,'A*(1-D*X**2)*EXP(-B*X**2)'/)
1050      FORMAT(9X,'5',15X,'CASO ESPONENZIALE SENO SU ASCISSA'
     ./25X,'A*EXP(-B*ABS(X))*SIN(C*ABS(X))/(C*ABS(X))'/)
1060      FORMAT(9X,'6',15X,'CASO NORMALE SENO SU ASCISSA' 
     ./25X,'A*EXP(-B*X**2)*SIN(C*ABS(X))/(C*ABS(X))'/) 
1070      FORMAT(9X,'7',15X,'CASO ESPONENZIALE COSENO' 
     ./25X,'A*EXP(-B*ABS(X))*COS(C*ABS(X))'/)  
1080      FORMAT(9X,'8',15X,'CASO NORMALE COSENO'  
     ./25X,'A*EXP(-B*X**2)*COS(C*ABS(X))'/)
1090      FORMAT(9X,'3',15X,'CASO ESPONENZIALE PARABOLA'
     ./25X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(6*D)'/)
1100      FORMAT(9X,'9',15X,'CASO ESPONENZIALE J1 SU ASCISSA'  
     ./25X,'2*A*EXP(-B*ABS(X))*J1(C*ABS(X))/(C*ABS(X))'/)  
1110      FORMAT(8X,'10',15X,'CASO NORMALE J1 SU ASCISSA'  
     ./25X,'2*A*EXP(-B*X**2)*J1(C*ABS(X))/(C*ABS(X))'/)
1120      FORMAT(8X,'11',15X,'CASO ESPONENZIALE J0'
     ./25X,'A*EXP(-B*ABS(X))*J0(C*ABS(X))'/)
1130      FORMAT(8X,'12',15X,'CASO NORMALE J0' 
     ./25X,'A*EXP(-B*X**2)*J0(C*ABS(X))'/) 
          END  
          SUBROUTINE MAT3(A,B,C)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION A(3,3),B(3,3),C(3,3)
          DO 10 I=1,3  
          DO 10 J=1,3  
          SUM=0.D0 
          DO 5 K=1,3
          SUM=SUM+A(I,K)*B(K,J)
5         CONTINUE 
          C(I,J)=SUM
10        CONTINUE 
          DO 20 I=1,3  
          DO 20 J=1,3  
          SUM=0.D0 
          DO 15 K=1,3  
          SUM=SUM+C(I,K)*A(K,J)
15        CONTINUE 
          B(I,J)=SUM
20        CONTINUE 
          RETURN
          END  
          SUBROUTINE MEDVA2(X,IPX,EX,S2X,SX,N) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),IPX(1)
          EX=0.D0  
          S2X=0.D0 
          DO 1 I=1,N
          L=IPX(I) 
          EX=EX+X(L)
          S2X=S2X+X(L)**2  
1         CONTINUE 
          EX=EX/N  
          S2X=(S2X/N-EX**2)*N/(N-1)
          SX=DSQRT(S2X)
          RETURN
          END  
          SUBROUTINE MIMA(COV,IMM,IX0,NCOV)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION COV(1),IX0(1),IMM(1)
          I1=IX0(1)
          I2=IX0(2)
          I3=IX0(3)
          GMIN=1.D10
          IF(I2.EQ.0) I2=NCOV  
          DO 10 I=I1,I2
          IF(COV(I).LT.GMIN) THEN  
          IMM(1)=I 
          GMIN=COV(I)  
          ENDIF
10        CONTINUE 
          IF(I2.EQ.NCOV) RETURN
          GMAX=-1.D10  
          IF(I3.EQ.0) I3=NCOV  
          DO 20 I=I2,I3
          IF(COV(I).GT.GMAX) THEN  
          IMM(2)=I 
          GMAX=COV(I)  
          ENDIF
20        CONTINUE 
          RETURN
          END  
          SUBROUTINE OTTIP1(X,Y,W,IPXW,IPYW,
     .    LPXW,LPYW,AINF,ASUP,STP,IWD, 
     .    NMAX,NWDMAX,SCX,SCY,SCW,IPTO,IAST,FORMA,ITRIM,TOLL)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20),WFORMA(30) 
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),W(1),IPXW(1),IPYW(1),
     .    LPXW(1),LPYW(1), 
     .    IWD(1),IPTO(1),IAST(1)
          DIMENSION SRV(2),NSRV(2),GMMW(2),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          DATA WFORMA/'(//6','X,''M','(X)''',',9X,','''SQM','(X)''',
     .    ',8X,','''RO-','1(X)',''',7X',',''SQ','M-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NW=0 
          ND=0 
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
C          DECODE(150,FORMA,DATO) NPT,XS,YS,WS  
          READ(CDATO,FORMA) NPT,XS,YS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY 
          W(ND)=WS*SCW 
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NW=NW+1  
          IPXW(NW)=ND  
          IPYW(NW)=ND  
          ENDIF
          GOTO 10  
20        CALL MEDVA2(W,IPXW,EW,S2W,SW,NW) 
          DO 30 I=1,ND 
          W(I)=W(I)-EW 
30        CONTINUE 
          CALL RBORDA(X,Y,IPXW,IPYW,LPXW,LPYW,NW,ND)
          S2NSE=1.D30 
          NSTP=IDINT((ASUP-AINF)/STP)+1
          WRITE(6,1200)
          DO 60 ISTP=1,NSTP
          RAGGIO=AINF+(ISTP-1)*STP 
          ALATO=RAGGIO*2.D0
          CALL COVAP(X,Y,W,IPXW,IPYW,LPXW,LPYW,NW,2,RAGGIO,ALATO,  
     .    GMMW,SRV,NSRV,IWD,NMAX,NWDMAX)
          SIGMA2=(S2W-GMMW(1)) 
          WRITE(6,1300) RAGGIO,DSQRT(SIGMA2)
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=ALATO  
          S2NSE=SIGMA2 
          SNW=DSQRT(S2W-GMMW(1))
          GMW=GMMW(1) 
          ENDIF
60        CONTINUE 
          NWM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMP(X,Y,W,IPXW,IPYW,LPXW,LPYW,IWD,SNW,OLATO,  
     .    1,NW,ND,NWM,NWDMAX,NMAX,IAST,1,2,TOLL)
          WRITE(2,1000)
          ENDIF
          WRITE(6,1100) OLATO/2.D0 
          IF(NWM.LT.9) THEN
          WFORMA(20)='I1,/'
          WFORMA(26)='I8) '
          ENDIF
          IF(NWM.GE.10.AND.NWM.LT.100) THEN
          WFORMA(20)='I2,/'
          WFORMA(26)='I9) '
          ENDIF
          IF(NWM.GE.100.AND.NWM.LT.1000) THEN  
          WFORMA(20)='I3,/'
          WFORMA(26)='I10)'
          ENDIF
          WRITE(6,WFORMA) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM  
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'RAGGIO OTTIMALE = ',F10.3,//)  
1200      FORMAT(///14X,'RAGGIO',7X,'SQM-N(X)',//) 
1300      FORMAT(5X,2F15.3)
          END  
          SUBROUTINE OTTIP2(X,Y,U,V,IPXU,IPYU,IPXV,IPYV,
     .    LPXU,LPYU,LPXV,LPYV,AINF,ASUP,STP,IWD,
     .    NMAX,NWDMAX,SCX,SCY,SCU,SCV,IPTO,IAST,FORMA,ITRIM,TOLL)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20),UFORMA(30),VFORMA(30)  
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),U(1),V(1),IPXU(1),IPYU(1),
     .    IPXV(1),IPYV(1),LPXU(1),LPYU(1), 
     .    LPXV(1),LPYV(1),IWD(1),IPTO(1),IAST(1)
          DIMENSION SRV(2),NSRV(2),GMMU(2),GMMV(2),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          DATA UFORMA/'(//6','X,''M','(X)''',',9X,','''SQM','(X)''',
     .    ',8X,','''RO-','1(X)',''',7X',',''SQ','M-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA VFORMA/'(//6','X,''M','(Y)''',',9X,','''SQM','(Y)''',
     .    ',8X,','''RO-','1(Y)',''',7X',',''SQ','M-N(','Y)'',',
     .    '11X,','''N(Y',')''  ',',7X,', 
     .    '''N(Y',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NU=0 
          NV=0 
          ND=0 
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
C          DECODE(150,FORMA,DATO) NPT,XS,YS,US,VS
          READ(CDATO,FORMA) NPT,XS,YS,US,VS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPXU(NU)=ND  
          IPYU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPXV(NV)=ND  
          IPYV(NV)=ND  
          ENDIF
          GOTO 10  
20        CALL MEDVA2(U,IPXU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPXV,EV,S2V,SV,NV) 
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
30        CONTINUE 
          CALL RBORDA(X,Y,IPXU,IPYU,LPXU,LPYU,NU,ND)
          CALL RBORDA(X,Y,IPXV,IPYV,LPXV,LPYV,NV,ND)
          S2NSE=1.D30 
          NSTP=IDINT((ASUP-AINF)/STP)+1
          WRITE(6,1200)
          DO 60 ISTP=1,NSTP
          RAGGIO=AINF+(ISTP-1)*STP 
          ALATO=RAGGIO*2.D0
          CALL COVAP(X,Y,U,IPXU,IPYU,LPXU,LPYU,NU,2,RAGGIO,ALATO,  
     .    GMMU,SRV,NSRV,IWD,NMAX,NWDMAX)
          CALL COVAP(X,Y,V,IPXV,IPYV,LPXV,LPYV,NV,2,RAGGIO,ALATO,  
     .    GMMV,SRV,NSRV,IWD,NMAX,NWDMAX)
          SIGMA2=(S2U-GMMU(1))+(S2V-GMMV(1))
          SIGMA2=SIGMA2/2.D0
          WRITE(6,1300) RAGGIO,DSQRT(S2U-GMMU(1)),DSQRT(S2V-GMMV(1)),  
     .    DSQRT(SIGMA2)
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=ALATO  
          S2NSE=SIGMA2 
          SNU=DSQRT(S2U-GMMU(1))
          SNV=DSQRT(S2V-GMMV(1))
          GMU=GMMU(1)
          GMV=GMMV(1)
          ENDIF
60        CONTINUE 
          NUM=0
          NVM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMP(X,Y,U,IPXU,IPYU,LPXU,LPYU,IWD,SNU,OLATO,  
     .    1,NU,ND,NUM,NWDMAX,NMAX,IAST,2,1,TOLL)
          CALL TRIMP(X,Y,V,IPXV,IPYV,LPXV,LPYV,IWD,SNV,OLATO,  
     .    2,NV,ND,NVM,NWDMAX,NMAX,IAST,1,2,TOLL)
          WRITE(2,1000)
          ENDIF
          WRITE(6,1100) OLATO/2.D0 
          IF(NUM.LT.9) THEN
          UFORMA(20)='I1,/'
          UFORMA(26)='I8) '
          ENDIF
          IF(NUM.GE.10.AND.NUM.LT.100) THEN
          UFORMA(20)='I2,/'
          UFORMA(26)='I9) '
          ENDIF
          IF(NUM.GE.100.AND.NUM.LT.1000) THEN  
          UFORMA(20)='I3,/'
          UFORMA(26)='I10)'
          ENDIF
          WRITE(6,UFORMA) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          IF(NVM.LT.9) THEN
          VFORMA(20)='I1,/'
          VFORMA(26)='I8) '
          ENDIF
          IF(NVM.GE.10.AND.NVM.LT.100) THEN
          VFORMA(20)='I2,/'
          VFORMA(26)='I9) '
          ENDIF
          IF(NVM.GE.100.AND.NVM.LT.1000) THEN  
          VFORMA(20)='I3,/'
          VFORMA(26)='I10)'
          ENDIF
          WRITE(6,VFORMA) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM  
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'RAGGIO OTTIMALE = ',F10.3,//)  
1200      FORMAT(///14X,'RAGGIO',7X,'SQM-N(X)',7X,'SQM-N(Y)',5X,
     .    'SQM-N(XUY)',//) 
1300      FORMAT(5X,4F15.3)
          END  
          SUBROUTINE OTTIP3(X,Y,U,V,W,IPXU,IPYU,IPXV,IPYV,IPXW,IPYW,
     .    LPXU,LPYU,LPXV,LPYV,LPXW,LPYW,AINF,ASUP,STP,IWD, 
     .    NMAX,NWDMAX,SCX,SCY,SCU,SCV,SCW,IPTO,IAST,FORMA,ITRIM,TOLL)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20),UFORMA(30),VFORMA(30),WFORMA(30)
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),U(1),V(1),W(1),IPXU(1),IPYU(1),  
     .    IPXV(1),IPYV(1),IPXW(1),IPYW(1),LPXU(1),LPYU(1), 
     .    LPXV(1),LPYV(1),LPXW(1),LPYW(1),IWD(1),IPTO(1),IAST(1)
          DIMENSION SRV(2),NSRV(2),GMMU(2),GMMV(2),GMMW(2),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          DATA UFORMA/'(//6','X,''M','(X)''',',9X,','''SQM','(X)''',
     .    ',8X,','''RO-','1(X)',''',7X',',''SQ','M-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA VFORMA/'(//6','X,''M','(Y)''',',9X,','''SQM','(Y)''',
     .    ',8X,','''RO-','1(Y)',''',7X',',''SQ','M-N(','Y)'',',
     .    '11X,','''N(Y',')''  ',',7X,', 
     .    '''N(Y',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA WFORMA/'(//6','X,''M','(Z)''',',9X,','''SQM','(Z)''',
     .    ',8X,','''RO-','1(Z)',''',7X',',''SQ','M-N(','Z)'',',
     .    '11X,','''N(Z',')''  ',',7X,', 
     .    '''N(Z',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NU=0 
          NV=0 
          NW=0 
          ND=0 
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
C          DECODE(150,FORMA,DATO) NPT,XS,YS,US,VS,WS
          READ(CDATO,FORMA) NPT,XS,YS,US,VS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          W(ND)=WS*SCW 
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPXU(NU)=ND  
          IPYU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPXV(NV)=ND  
          IPYV(NV)=ND  
          ENDIF
          IF(DATO(IPTO(5)).EQ.PTO.AND.DATO(IAST(3)).NE.AST) THEN
          NW=NW+1  
          IPXW(NW)=ND  
          IPYW(NW)=ND  
          ENDIF
          GOTO 10  
20        CALL MEDVA2(U,IPXU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPXV,EV,S2V,SV,NV) 
          CALL MEDVA2(W,IPXW,EW,S2W,SW,NW) 
          DO 30 I=1,NU 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
30        CONTINUE 
          CALL RBORDA(X,Y,IPXU,IPYU,LPXU,LPYU,NU,ND)
          CALL RBORDA(X,Y,IPXV,IPYV,LPXV,LPYV,NV,ND)
          CALL RBORDA(X,Y,IPXW,IPYW,LPXW,LPYW,NW,ND)
          OLATO=0.D0
          S2NSE=1.D30 
          NSTP=IDINT((ASUP-AINF)/STP)+1
          WRITE(6,1200)
          DO 60 ISTP=1,NSTP
          RAGGIO=AINF+(ISTP-1)*STP 
          ALATO=RAGGIO*2.D0
          CALL COVAP(X,Y,U,IPXU,IPYU,LPXU,LPYU,NU,2,RAGGIO,ALATO,  
     .    GMMU,SRV,NSRV,IWD,NMAX,NWDMAX)
          CALL COVAP(X,Y,V,IPXV,IPYV,LPXV,LPYV,NV,2,RAGGIO,ALATO,  
     .    GMMV,SRV,NSRV,IWD,NMAX,NWDMAX)
          CALL COVAP(X,Y,W,IPXW,IPYW,LPXW,LPYW,NW,2,RAGGIO,ALATO,  
     .    GMMW,SRV,NSRV,IWD,NMAX,NWDMAX)
          SIGMA2=(S2U-GMMU(1))+(S2V-GMMV(1))+(S2W-GMMW(1)) 
          SIGMA2=SIGMA2/3.D0
          WRITE(6,1300) RAGGIO,DSQRT(S2U-GMMU(1)),DSQRT(S2V-GMMV(1)),  
     .    DSQRT(S2W-GMMW(1)),DSQRT(SIGMA2) 
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=ALATO  
          S2NSE=SIGMA2 
          SNU=DSQRT(S2U-GMMU(1))
          SNV=DSQRT(S2V-GMMV(1))
          SNW=DSQRT(S2W-GMMW(1))
          GMU=GMMU(1)
          GMV=GMMV(1)
          GMW=GMMW(1)
          ENDIF
60        CONTINUE 
          NUM=0
          NVM=0
          NWM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMP(X,Y,U,IPXU,IPYU,LPXU,LPYU,IWD,SNU,OLATO,  
     .    1,NU,ND,NUM,NWDMAX,NMAX,IAST,1,2,TOLL)
          CALL TRIMP(X,Y,V,IPXV,IPYV,LPXV,LPYV,IWD,SNV,OLATO,  
     .    2,NV,ND,NVM,NWDMAX,NMAX,IAST,2,1,TOLL)
          CALL TRIMP(X,Y,W,IPXW,IPYW,LPXW,LPYW,IWD,SNW,OLATO,  
     .    3,NW,ND,NWM,NWDMAX,NMAX,IAST,1,2,TOLL)
          WRITE(2,1000)
          ENDIF
          WRITE(6,1100) OLATO/2.D0 
          IF(NUM.LT.9) THEN
          UFORMA(20)='I1,/'
          UFORMA(26)='I8) '
          ENDIF
          IF(NUM.GE.10.AND.NUM.LT.100) THEN
          UFORMA(20)='I2,/'
          UFORMA(26)='I9) '
          ENDIF
          IF(NUM.GE.100.AND.NUM.LT.1000) THEN  
          UFORMA(20)='I3,/'
          UFORMA(26)='I10)'
          ENDIF
          WRITE(6,UFORMA) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          IF(NVM.LT.9) THEN
          VFORMA(20)='I1,/'
          VFORMA(26)='I8) '
          ENDIF
          IF(NVM.GE.10.AND.NVM.LT.100) THEN
          VFORMA(20)='I2,/'
          VFORMA(26)='I9) '
          ENDIF
          IF(NVM.GE.100.AND.NVM.LT.1000) THEN  
          VFORMA(20)='I3,/'
          VFORMA(26)='I10)'
          ENDIF
          WRITE(6,VFORMA) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM  
          IF(NWM.LT.9) THEN
          WFORMA(20)='I1,/'
          WFORMA(26)='I8) '
          ENDIF
          IF(NWM.GE.10.AND.NWM.LT.100) THEN
          WFORMA(20)='I2,/'
          WFORMA(26)='I9) '
          ENDIF
          IF(NWM.GE.100.AND.NWM.LT.1000) THEN  
          WFORMA(20)='I3,/'
          WFORMA(26)='I10)'
          ENDIF
          WRITE(6,WFORMA) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM  
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'RAGGIO OTTIMALE = ',F10.3,//)  
1200      FORMAT(///14X,'RAGGIO',7X,'SQM-N(X)',7X,'SQM-N(Y)',7X,
     .    'SQM-N(Z)',3X,'SQM-N(XUYUZ)'//)  
1300      FORMAT(5X,5F15.3)
          END  
          SUBROUTINE OTTIT1(Z,W,IPZW,LPZW, 
     .    AINF,ASUP,STP,NMAX,SCZ,SCW,IPTO,IAST,FORMA,ITRIM,TOLL)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20),WFORMA(30) 
          CHARACTER*120 CDATO
          DIMENSION Z(1),W(1),IPZW(1),LPZW(1),IPTO(1),IAST(1)  
          DIMENSION SRV(2),NSRV(2),GMMW(2),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          DATA WFORMA/'(//6','X,''M','(X)''',',9X,','''SQM','(X)''',
     .    ',8X,','''RO-','1(X)',''',7X',',''SQ','M-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NW=0 
          ND=0 
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
C          DECODE(150,FORMA,DATO) NPT,TS,WS 
          READ(CDATO,FORMA) NPT,TS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          Z(ND)=TS*SCZ 
          W(ND)=WS*SCW 
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NW=NW+1  
          IPZW(NW)=ND  
          ENDIF
          GOTO 10  
20        CALL MEDVA2(W,IPZW,EW,S2W,SW,NW) 
          DO 30 I=1,ND 
          W(I)=W(I)-EW 
30        CONTINUE 
          CALL RIORDA(Z,IPZW,NW,ND,1)  
          CALL PUNTA(IPZW,LPZW,0,1,NW,ND)  
          S2NSE=1.D30 
          NSTP=IDINT((ASUP-AINF)/STP)+1
          OLATO=0.D0
          WRITE(6,1200)
          DO 60 ISTP=1,NSTP
          RAGGIO=AINF+(ISTP-1)*STP 
          ALATO=RAGGIO*2.D0
          CALL COVAT(Z,W,IPZW,NW,2,RAGGIO,ALATO,GMMW,SRV,NSRV,NMAX)
          SIGMA2=(S2W-GMMW(1)) 
          WRITE(6,1300) RAGGIO,DSQRT(SIGMA2)
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=ALATO  
          S2NSE=SIGMA2 
          SNW=DSQRT(S2W-GMMW(1))
          GMW=GMMW(1) 
          ENDIF
60        CONTINUE 
          NWM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMT(Z,W,IPZW,LPZW,SNW,OLATO,NMAX,1,NW,ND,NWM,IAST,1,2,
     .TOLL)
          WRITE(2,1000)
          ENDIF
          WRITE(6,1100) OLATO/2.D0 
          IF(NWM.LT.9) THEN
          WFORMA(20)='I1,/'
          WFORMA(26)='I8) '
          ENDIF
          IF(NWM.GE.10.AND.NWM.LT.100) THEN
          WFORMA(20)='I2,/'
          WFORMA(26)='I9) '
          ENDIF
          IF(NWM.GE.100.AND.NWM.LT.1000) THEN  
          WFORMA(20)='I3,/'
          WFORMA(26)='I10)'
          ENDIF
          WRITE(6,WFORMA) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM  
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'SEMI-INTERVALLO OTTIMALE = ',F10.3,//) 
1200      FORMAT(///5X,'SEMI-INTERVALLO',7X,'SQM-N(X)'//)  
1300      FORMAT(5X,2F15.3)
          END 
          SUBROUTINE OTTIT2(Z,U,V,IPZU,IPZV,LPZU,LPZV, 
     .    AINF,ASUP,STP,NMAX,SCZ,SCU,SCV,IPTO,IAST,FORMA,ITRIM,TOLL)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20),UFORMA(30),VFORMA(30)  
          CHARACTER*120 CDATO
          DIMENSION Z(1),U(1),V(1),IPZU(1),IPZV(1),
     .    LPZU(1),LPZV(1),IPTO(1),IAST(1)  
          DIMENSION SRV(2),NSRV(2),GMMU(2),GMMV(2),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          DATA UFORMA/'(//6','X,''M','(X)''',',9X,','''SQM','(X)''',
     .    ',8X,','''RO-','1(X)',''',7X',',''SQ','M-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA VFORMA/'(//6','X,''M','(Y)''',',9X,','''SQM','(Y)''',
     .    ',8X,','''RO-','1(Y)',''',7X',',''SQ','M-N(','Y)'',',
     .    '11X,','''N(Y',')''  ',',7X,', 
     .    '''N(Y',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NU=0 
          NV=0 
          ND=0 
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
C          DECODE(150,FORMA,DATO) NPT,TS,US,VS  
          READ(CDATO,FORMA) NPT,TS,US,VS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          Z(ND)=TS*SCZ 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPZU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPZV(NV)=ND  
          ENDIF
          GOTO 10  
20        CALL MEDVA2(U,IPZU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPZV,EV,S2V,SV,NV) 
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
30        CONTINUE 
          CALL RIORDA(Z,IPZU,NU,ND,1)  
          CALL RIORDA(Z,IPZV,NV,ND,1)  
          CALL PUNTA(IPZU,LPZU,0,1,NU,ND)  
          CALL PUNTA(IPZV,LPZV,0,1,NV,ND)  
          S2NSE=1.D30 
          NSTP=IDINT((ASUP-AINF)/STP)+1
          OLATO=0.D0
          WRITE(6,1200)
          DO 60 ISTP=1,NSTP
          RAGGIO=AINF+(ISTP-1)*STP 
          ALATO=RAGGIO*2.D0
          CALL COVAT(Z,U,IPZU,NU,2,RAGGIO,ALATO,GMMU,SRV,NSRV,NMAX)
          CALL COVAT(Z,V,IPZV,NV,2,RAGGIO,ALATO,GMMV,SRV,NSRV,NMAX)
          SIGMA2=(S2U-GMMU(1))+(S2V-GMMV(1))
          SIGMA2=SIGMA2/2.D0
          WRITE(6,1300) RAGGIO,DSQRT(S2U-GMMU(1)),DSQRT(S2V-GMMV(1)),  
     .    DSQRT(SIGMA2)
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=ALATO  
          S2NSE=SIGMA2 
          SNU=DSQRT(S2U-GMMU(1))
          SNV=DSQRT(S2V-GMMV(1))
          GMU=GMMU(1)
          GMV=GMMV(1)
          ENDIF
60        CONTINUE 
          NUM=0
          NVM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMT(Z,U,IPZU,LPZU,SNU,OLATO,NMAX,1,NU,ND,NUM,IAST,2,1,
     .TOLL)
          CALL TRIMT(Z,V,IPZV,LPZV,SNV,OLATO,NMAX,2,NV,ND,NVM,IAST,1,2,
     .TOLL)
          WRITE(2,1000)
          ENDIF
          WRITE(6,1100) OLATO/2.D0 
          IF(NUM.LT.9) THEN
          UFORMA(20)='I1,/'
          UFORMA(26)='I8) '
          ENDIF
          IF(NUM.GE.10.AND.NUM.LT.100) THEN
          UFORMA(20)='I2,/'
          UFORMA(26)='I9) '
          ENDIF
          IF(NUM.GE.100.AND.NUM.LT.1000) THEN  
          UFORMA(20)='I3,/'
          UFORMA(26)='I10)'
          ENDIF
          WRITE(6,UFORMA) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          IF(NVM.LT.9) THEN
          VFORMA(20)='I1,/'
          VFORMA(26)='I8) '
          ENDIF
          IF(NVM.GE.10.AND.NVM.LT.100) THEN
          VFORMA(20)='I2,/'
          VFORMA(26)='I9) '
          ENDIF
          IF(NVM.GE.100.AND.NVM.LT.1000) THEN  
          VFORMA(20)='I3,/'
          VFORMA(26)='I10)'
          ENDIF
          WRITE(6,VFORMA) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM  
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'SEMI-INTERVALLO OTTIMALE = ',F10.3,//) 
1200      FORMAT(///5X,'SEMI-INTERVALLO',7X,'SQM-N(X)',7X,'SQM-N(Y)',  
     .    5X,'SQM-N(XUY)',//)  
1300      FORMAT(5X,F15.3,3F15.3)  
          END 
          SUBROUTINE OTTIT3(Z,U,V,W,IPZU,IPZV,IPZW,LPZU,LPZV,LPZW, 
     .    AINF,ASUP,STP,NMAX,SCZ,SCU,SCV,SCW,IPTO,IAST,FORMA,ITRIM,TOLL)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA(20),UFORMA(30),VFORMA(30),WFORMA(30)
          CHARACTER*120 CDATO
          DIMENSION Z(1),U(1),V(1),W(1),IPZU(1),IPZV(1),IPZW(1),
     .    LPZU(1),LPZV(1),LPZW(1),IAST(1),IPTO(1)  
          DIMENSION SRV(2),NSRV(2),GMMU(2),GMMV(2),GMMW(2),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          DATA UFORMA/'(//6','X,''M','(X)''',',9X,','''SQM','(X)''',
     .    ',8X,','''RO-','1(X)',''',7X',',''SQ','M-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA VFORMA/'(//6','X,''M','(Y)''',',9X,','''SQM','(Y)''',
     .    ',8X,','''RO-','1(Y)',''',7X',',''SQ','M-N(','Y)'',',
     .    '11X,','''N(Y',')''  ',',7X,', 
     .    '''N(Y',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA WFORMA/'(//6','X,''M','(Z)''',',9X,','''SQM','(Z)''',
     .    ',8X,','''RO-','1(Z)',''',7X',',''SQ','M-N(','Z)'',',
     .    '11X,','''N(Z',')''  ',',7X,', 
     .    '''N(Z',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NU=0 
          NV=0 
          NW=0 
          ND=0 
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
C          DECODE(150,FORMA,DATO) NPT,TS,US,VS,WS
          READ(CDATO,FORMA) NPT,TS,US,VS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          Z(ND)=TS*SCZ 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          W(ND)=WS*SCW 
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPZU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPZV(NV)=ND  
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(3)).NE.AST) THEN
          NW=NW+1  
          IPZW(NW)=ND  
          ENDIF
          GOTO 10  
20        CALL MEDVA2(U,IPZU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPZV,EV,S2V,SV,NV) 
          CALL MEDVA2(W,IPZW,EW,S2W,SW,NW) 
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
30        CONTINUE 
          CALL RIORDA(Z,IPZU,NU,ND,1)  
          CALL RIORDA(Z,IPZV,NV,ND,1)  
          CALL RIORDA(Z,IPZW,NW,ND,1)  
          CALL PUNTA(IPZU,LPZU,0,1,NU,ND)  
          CALL PUNTA(IPZV,LPZV,0,1,NV,ND)  
          CALL PUNTA(IPZW,LPZW,0,1,NW,ND)  
          S2NSE=1.D30 
          NSTP=IDINT((ASUP-AINF)/STP)+1
          OLATO=0.D0
          WRITE(6,1200)
          DO 60 ISTP=1,NSTP
          RAGGIO=AINF+(ISTP-1)*STP 
          ALATO=RAGGIO*2.D0
          CALL COVAT(Z,U,IPZU,NU,2,RAGGIO,ALATO,GMMU,SRV,NSRV,NMAX)
          CALL COVAT(Z,V,IPZV,NV,2,RAGGIO,ALATO,GMMV,SRV,NSRV,NMAX)
          CALL COVAT(Z,W,IPZW,NW,2,RAGGIO,ALATO,GMMW,SRV,NSRV,NMAX)
          SIGMA2=(S2U-GMMU(1))+(S2V-GMMV(1))+(S2W-GMMW(1)) 
          SIGMA2=SIGMA2/3.D0
          WRITE(6,1300) RAGGIO,DSQRT(S2U-GMMU(1)),DSQRT(S2V-GMMV(1)),  
     .    DSQRT(S2W-GMMW(1)),DSQRT(SIGMA2) 
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=ALATO  
          S2NSE=SIGMA2 
          SNU=DSQRT(S2U-GMMU(1))
          SNV=DSQRT(S2V-GMMV(1))
          SNW=DSQRT(S2W-GMMW(1))
          GMU=GMMU(1)
          GMV=GMMV(1)
          GMW=GMMW(1)
          ENDIF
60        CONTINUE 
          NUM=0
          NVM=0
          NWM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMT(Z,U,IPZU,LPZU,SNU,OLATO,NMAX,1,NU,ND,NUM,IAST,1,2,
     .TOLL)
          CALL TRIMT(Z,V,IPZV,LPZV,SNV,OLATO,NMAX,2,NV,ND,NVM,IAST,2,1,
     .TOLL)
          CALL TRIMT(Z,W,IPZW,LPZW,SNW,OLATO,NMAX,3,NW,ND,NWM,IAST,1,2,
     .TOLL)
          WRITE(2,1000)
          ENDIF
          WRITE(6,1100) OLATO/2.D0 
          IF(NUM.LT.9) THEN
          UFORMA(20)='I1,/'
          UFORMA(26)='I8) '
          ENDIF
          IF(NUM.GE.10.AND.NUM.LT.100) THEN
          UFORMA(20)='I2,/'
          UFORMA(26)='I9) '
          ENDIF
          IF(NUM.GE.100.AND.NUM.LT.1000) THEN  
          UFORMA(20)='I3,/'
          UFORMA(26)='I10)'
          ENDIF
          WRITE(6,UFORMA) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          IF(NVM.LT.9) THEN
          VFORMA(20)='I1,/'
          VFORMA(26)='I8) '
          ENDIF
          IF(NVM.GE.10.AND.NVM.LT.100) THEN
          VFORMA(20)='I2,/'
          VFORMA(26)='I9) '
          ENDIF
          IF(NVM.GE.100.AND.NVM.LT.1000) THEN  
          VFORMA(20)='I3,/'
          VFORMA(26)='I10)'
          ENDIF
          WRITE(6,VFORMA) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM  
          IF(NWM.LT.9) THEN
          WFORMA(20)='I1,/'
          WFORMA(26)='I8) '
          ENDIF
          IF(NWM.GE.10.AND.NWM.LT.100) THEN
          WFORMA(20)='I2,/'
          WFORMA(26)='I9) '
          ENDIF
          IF(NWM.GE.100.AND.NWM.LT.1000) THEN  
          WFORMA(20)='I3,/'
          WFORMA(26)='I10)'
          ENDIF
          WRITE(6,WFORMA) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM  
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'SEMI-INTERVALLO OTTIMALE = ',F10.3,//) 
1200      FORMAT(///5X,'SEMI-INTERVALLO',7X,'SQM-N(X)',7X,'SQM-N(Y)',7X,
     .    'SQM-N(Z)',3X,'SQM-N(XUYUZ)',//) 
1300      FORMAT(5X,5F15.3)
          END 
          SUBROUTINE PACORO(TITOLO,ICASO,ITIPO,ISCRI,ITRIM,IQUO,TOLL,  
     .AINF,ASUP,STP,SCX,SCY,SCZ,SCU,SCV,SCW,FORMA) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*4 TITOLO,FORMA 
          DIMENSION TITOLO(20),FORMA(20)
          IQUO=0
          WRITE(*,1000)
          READ(*,1010) TITOLO  
10        WRITE(*,1020)
          READ(*,*,ERR=10) ICASO
15        WRITE(*,1030)
          READ(*,*,ERR=15) ITIPO
20        WRITE(*,1040)
          READ(*,*,ERR=20) ISCRI
25        WRITE(*,1050)
          READ(*,*,ERR=25) ITRIM
          IF(ITRIM.EQ.1) THEN  
28        WRITE(*,1055)
          READ(*,*,ERR=28) TOLL
          ENDIF
          IF(ICASO.EQ.2) THEN  
30        WRITE(*,1060)
          READ(*,*,ERR=30) IQUO
          ENDIF
35        WRITE(*,1070)
          READ(*,*,ERR=35) AINF,ASUP,STP
          IF(ICASO.EQ.1) THEN  
          IF(ITIPO.EQ.1) THEN  
40        WRITE(*,1080)
          READ(*,*,ERR=40) SCZ,SCW 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
45        WRITE(*,1090)
          READ(*,*,ERR=45) SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
50        WRITE(*,1100)
          READ(*,*,ERR=50) SCZ,SCU,SCV,SCW 
          ENDIF
          ENDIF
          IF(ICASO.EQ.2) THEN  
          IF(ITIPO.EQ.1) THEN  
55        WRITE(*,1110)
          READ(*,*,ERR=55) SCX,SCY,SCW 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
60        WRITE(*,1120)
          READ(*,*,ERR=60) SCX,SCY,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
65        WRITE(*,1130)
          READ(*,*,ERR=65) SCX,SCY,SCU,SCV,SCW 
          ENDIF
          ENDIF
          WRITE(*,1140)
          READ(*,1010) FORMA
          RETURN
1000      FORMAT(5X,'TITOLO')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASO:',/15X,'UNIDIMENSIONALE = 1',
     ./15X,'BIDIMENSIONALE  = 2')  
1030      FORMAT(5X,'TIPO:',/15X,'UNA COMPONENTE = 1', 
     ./15X,'DUE COMPONENTI = 2',
     ./15X,'TRE COMPONENTI = 3')
1040      FORMAT(5X,'ISCRI:',/15X,'SCRIVE DATI INPUT = 1', 
     ./15X,'NON SCRIVE DATI INPUT = 0')
1050      FORMAT(5X,'ITRIM:',/15X,'ELIMINA ERRORI GROSSOLANI = 1', 
     .'     (DATI OUTPUT FILE 2.)',
     ./15X,'NON ELIMINA ERRORI GROSSOLANI = 0')
1055      FORMAT(5X,'TOLL:',/15X,'VALORE CRITICO PER LA REIEZIONE ',
     .'DEGLI ERRORI GROSSOLANI')
1060      FORMAT(5X,'IQUO:',/15X,'PRESENZA QUOTE = 1', 
     ./15X,'ASSENZA QUOTE  = 0')
1070      FORMAT(5X,'AINF ASUP STP:',/15X,'ESTREMI INFERIORE E '
     .'SUPERIORE, PASSO',/15X,'(RICERCA MIGLIOR DELTA : '  
     .'GAMMA(1) = MAX)')
1080      FORMAT(5X,'SCALE:',/15X,'SCALA TEMPO T', 
     ./15X,'SCALA COMPONENTE W')
1090      FORMAT(5X,'SCALE:',/15X,'SCALA TEMPO T', 
     ./15X,'SCALE COMPONENTI U V') 
1100      FORMAT(5X,'SCALE:',/15X,'SCALA TEMPO T', 
     ./15X,'SCALE COMPONENTI U V W')
1110      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y',  
     ./15X,'SCALA COMPONENTE W')
1120      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y',  
     ./15X,'SCALE COMPONENTI U V') 
1130      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y',  
     ./15X,'SCALE COMPONENTI U V W')
1140      FORMAT(5X,'FORMATO DATI INPUT')  
1150      FORMAT(5X,'DATI INPUT FILE 1.')  
          END
          SUBROUTINE PACORR(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IQUO,
     .NCOR,RAGGIO,ALATO,SCX,SCY,SCZ,SCU,SCV,SCW,FORMA) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*4 TITOLO,FORMA 
          DIMENSION TITOLO(20),FORMA(20)
          WRITE(*,1000)
          READ(*,1010) TITOLO  
10        WRITE(*,1020)
          READ(*,*,ERR=10) ICASO
15        WRITE(*,1030)
          READ(*,*,ERR=15) ITIPO
20        WRITE(*,1040)
          READ(*,*,ERR=20) ISCRI
25        WRITE(*,1050)
          READ(*,*,ERR=25) IFILE
          IF(ICASO.EQ.2) THEN  
30        WRITE(*,1060)
          READ(*,*,ERR=30) IQUO
          ENDIF
35        WRITE(*,1070)
          READ(*,*,ERR=35) NCOR,RAGGIO 
          ALATO=2.D0*NCOR*RAGGIO
          IF(ICASO.EQ.1) THEN  
          IF(ITIPO.EQ.1) THEN  
40        WRITE(*,1080)
          READ(*,*,ERR=40) SCZ,SCW 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
45        WRITE(*,1090)
          READ(*,*,ERR=45) SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
50        WRITE(*,1100)
          READ(*,*,ERR=50) SCZ,SCU,SCV,SCW 
          ENDIF
          ENDIF
          IF(ICASO.EQ.2) THEN  
          IF(ITIPO.EQ.1) THEN  
55        WRITE(*,1110)
          READ(*,*,ERR=55) SCX,SCY,SCW 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
60        WRITE(*,1120)
          READ(*,*,ERR=60) SCX,SCY,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
65        WRITE(*,1130)
          READ(*,*,ERR=65) SCX,SCY,SCU,SCV,SCW 
          ENDIF
          ENDIF
          WRITE(*,1140)
          READ(*,1010) FORMA
          RETURN
1000      FORMAT(5X,'TITOLO')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASO:',/15X,'UNIDIMENSIONALE = 1',
     ./15X,'BIDIMENSIONALE  = 2')  
1030      FORMAT(5X,'TIPO:',/15X,'UNA COMPONENTE = 1', 
     ./15X,'DUE COMPONENTI = 2',
     ./15X,'TRE COMPONENTI = 3')
1040      FORMAT(5X,'ISCRI:',/15X,'SCRIVE DATI INPUT = 1', 
     ./15X,'NON SCRIVE DATI INPUT = 0')
1050      FORMAT(5X,'IFILE:',/15X,'SALVA CORRELOGRAMMI = 1',
     .5X,'(DATI OUTPUT FILE 2.)'
     ./15X,'NON SALVA CORRELOGRAMMI = 0')  
1060      FORMAT(5X,'IQUO:',/15X,'PRESENZA QUOTE = 1', 
     ./15X,'ASSENZA QUOTE  = 0')
1070      FORMAT(5X,'NCOR DELTA:',/15X,'NUMERO DI PASSI',  
     ./15X,'PASSO DEL CORRELOGRAMMA')  
1080      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO T'
     ./15X,'SCALA COMPONENTE W')
1090      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO T'
     ./15X,'SCALE COMPONENTI U V') 
1100      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO T'
     ./15X,'SCALE COMPONENTI U V W')
1110      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y',  
     ./15X,'SCALA COMPONENTE W')
1120      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y',  
     ./15X,'SCALE COMPONENTI U V') 
1130      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y',  
     ./15X,'SCALE COMPONENTI U V W')
1140      FORMAT(5X,'FORMATO DATI INPUT')  
1150      FORMAT(5X,'DATI INPUT FILE 1.')  
          END  
          SUBROUTINE PAINTC(TITOLO,ICASE,ITIPO,ISCAL,IFILE,XCOVT,DOWN) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*4 TITOLO
          DIMENSION XCOVT(1),DOWN(1),TITOLO(20)
          WRITE(*,1000)
          READ(*,1010) TITOLO  
          WRITE(*,1020)
          READ(*,*) ICASE  
          WRITE(*,1030)
          READ(*,*) ITIPO  
          WRITE(*,1040)
          READ(*,*) IFILE  
          WRITE(*,1050)
          READ(*,*) ISCAL  
          NTER=0
          IF(ITIPO.EQ.1) NTER=1
          IF(ITIPO.EQ.2) NTER=3
          IF(ITIPO.EQ.3) NTER=6
          IF(ISCAL.EQ.1) THEN  
          WRITE(*,1060)
          READ(*,*) (XCOVT(I),I=1,NTER)
          DO 1 I=1,NTER
          XCOVT(I)=XCOVT(I)/2.D0
1         CONTINUE 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
          WRITE(*,1070)
          READ(*,*) DOWN(1)
          ENDIF
          IF(ITIPO.EQ.3) THEN  
          WRITE(*,1070)
          READ(*,*) (DOWN(I),I=1,3)
          ENDIF
          RETURN
1000      FORMAT(5X,'TITOLO')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASO:'/15X,'UNIDIMENSIONALE = 1'  
     ./15X,'BIDIMENSIONALE = 2')
1030      FORMAT(5X,'TIPO:'/15X,'UNA AUTOCOVARIANZA = 1'
     ./15X,'DUE AUTO E UNA CROSS COVARIANZA = 2'
     ./15X,'TRE AUTO E TRE CROSS COVARIANZE = 3')  
1040      FORMAT(5X,'FILE:'/15X,'SALVA VALORI DELLE FUNZIONI EMPIRICHE '
     .'E TEORICHE'/15X,'SALVA I PARAMETRI DELLE FUNZIONI TEORICHE = 1' 
     ./15X,'ALTRIMENTI = 0'/15X,'(FILE DI OUTPUT 2. E 3.)')
1050      FORMAT(5X,'SCAL:'/15X,'PRESENZA FUNZIONE DI COVARIANZA FINITA'
     .' = 1'/15X,'ASSENZA FUNZIONE DI COVARIANZA FINITA = 0')  
1060      FORMAT(5X,'XCOVT:'/15X,'PUNTI DI AZZERAMENTO FUNZIONE EMPIRICA
     .'/15X,' = 1 VALORE PER TIPO UGUALE A UNO'
     ./15X,' = 3 VALORI PER TIPO UGUALE A DUE' 
     ./15X,' = 6 VALORI PER TIPO UGUALE A TRE')
1070      FORMAT(5X,'DOWN:'/15X,'EVENTUALE FATTORE DI RIDUZIONE'
     ./15X,'PER LE FUNZIONI DI CROSS COVARIANZA'
     ./15X,' = 1 VALORE PER TIPO UGUALE DUE'
     ./15X,' = 3 VALORI PER TIPO UGUALE TRE')  
          END  
          SUBROUTINE PAPOLI(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IQUO,
     .IGRA,NGRA,TOLL,SCX,SCY,SCZ,SCU,SCV,SCW,FORMA,FORMA2,FORMA3)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*4 TITOLO,FORMA,FORMA2,FORMA3
          DIMENSION TITOLO(20),FORMA(20),FORMA2(20),FORMA3(20) 
          WRITE(*,1000)
          READ(*,1010) TITOLO  
10        WRITE(*,1020)
          READ(*,*,ERR=10) ICASO
15        WRITE(*,1030)
          READ(*,*,ERR=15) ITIPO
20        WRITE(*,1040)
          READ(*,*,ERR=20) ISCRI
25        WRITE(*,1050)
          READ(*,*,ERR=25) IFILE
          IQUO=0
          IF(ICASO.EQ.2) THEN  
30        WRITE(*,1060)
          READ(*,*,ERR=30) IQUO
          ENDIF
35        WRITE(*,1070)
          READ(*,*,ERR=35) IGRA,NGRA,TOLL  
          NGRA=NGRA+1  
          IF(ICASO.EQ.1) THEN  
          IF(ITIPO.EQ.1) THEN  
40        WRITE(*,1080)
          READ(*,*,ERR=40) SCZ,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
45        WRITE(*,1090)
          READ(*,*,ERR=45) SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
50        WRITE(*,1100)
          READ(*,*,ERR=50) SCZ,SCU,SCV,SCW 
          ENDIF
          ENDIF
          IF(ICASO.EQ.2) THEN  
          IF(ITIPO.EQ.1) THEN  
55        WRITE(*,1110)
          READ(*,*,ERR=55) SCX,SCY,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
60        WRITE(*,1120)
          READ(*,*,ERR=60) SCX,SCY,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
65        WRITE(*,1130)
          READ(*,*,ERR=65) SCX,SCY,SCU,SCV,SCW 
          ENDIF
          ENDIF
          WRITE(*,1140)
          READ(*,1010) FORMA
          IF(ITIPO.GE.2) THEN  
          WRITE(*,1142)
          READ(*,1010) FORMA2  
          ENDIF
          IF(ITIPO.EQ.3) THEN  
          WRITE(*,1143)
          READ(*,1010) FORMA3  
          ENDIF
          RETURN
1000      FORMAT(5X,'TITOLO')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASO:',/15X,'UNIDIMENSIONALE = 1',
     ./15X,'BIDIMENSIONALE  = 2')  
1030      FORMAT(5X,'TIPO:',/15X,'UNA COMPONENTE = 1', 
     ./15X,'DUE COMPONENTI = 2',
     ./15X,'TRE COMPONENTI = 3')
1040      FORMAT(5X,'ISCRI:',/15X,'SCRIVE DATI INPUT = 1', 
     ./15X,'NON SCRIVE DATI INPUT = 0')
1050      FORMAT(5X,'IFILE:',/15X,'SALVA SCARTI-RESIDUI E ',
     .'COEFFICIENTI POLINOMIO = 1',/15X,'(DATI OUTPUT FILES 2. E 3.)', 
     ./15X,'NON SALVA SCARTI-RESIDUI E COEFFICIENTI POLINOMIO = 0')
1060      FORMAT(5X,'IQUO:',/15X,'PRESENZA QUOTE = 1', 
     ./15X,'ASSENZA QUOTE  = 0')
1070      FORMAT(5X,'IGRA NGRA TOLL:',/15X,
     .'IMPONE GRADO DEL POLINOMIO = 1',/15X,
     .'NON IMPONE GRADO DEL POLINOMIO = 0',/15X,
     .'GRADO MAX DEL POLINOMIO',/15X,  
     .'(GRADO IMPOSTO SE IGRA = 1)'
     ./15X,'VALORE CRITICO SULLA SIGNIFICATIVITA DEI PARAMETRI',
     ./15X,'E PER LA REIEZIONE DEGLI ERRORI GROSSOLANI')
1080      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALA COMPONENTE X')
1090      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALE COMPONENTI X Y') 
1100      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALE COMPONENTI X Y Z')
1110      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE PIANE',
     ./15X,'SCALA COMPONENTE X')
1120      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE PIANE',
     ./15X,'SCALE COMPONENTI X Y') 
1130      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE PIANE',
     ./15X,'SCALE COMPONENTI X Y Z')
1140      FORMAT(5X,'FORMATO DATI INPUT (I COMPONENTE)')
1142      FORMAT(5X,'FORMATO DATI INPUT (II COMPONENTE)')  
1143      FORMAT(5X,'FORMATO DATI INPUT (III COMPONENTE)') 
1150      FORMAT(5X,'DATI INPUT FILE 1.')  
          END  
          SUBROUTINE PASPLI(TITOLO,ICASO,ITIPO,IFILE,IQUO,ITABLE,  
     .TOLL,STP1,STP2,STP3,SCX,SCY,SCZ,SCU,SCV,SCW, 
     .FORMA,FORMA2,FORMA3) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*4 TITOLO,FORMA,FORMA2,FORMA3
          DIMENSION TITOLO(20),FORMA(20),FORMA2(20),FORMA3(20) 
          WRITE(*,1000)
          READ(*,1005) TITOLO  
5         WRITE(*,1025)
          READ(*,*,ERR=5) ITIPO
10        WRITE(*,1035)
          READ(*,*,ERR=10) IFILE
15        WRITE(*,1030)
          READ(*,*,ERR=15) TOLL
          IF(ICASO.EQ.1) THEN  
          IF(ITIPO.EQ.1) THEN  
20        WRITE(*,1010)
          READ(*,*,ERR=20) STP1
25        WRITE(*,1050)
          READ(*,*,ERR=25) SCZ,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
30        WRITE(*,1015)
          READ(*,*,ERR=30) STP1,STP2
35        WRITE(*,1055)
          READ(*,*,ERR=35) SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
40        WRITE(*,1020)
          READ(*,*,ERR=40) STP1,STP2,STP3  
45        WRITE(*,1060)
          READ(*,*,ERR=45) SCZ,SCU,SCV,SCW 
          ENDIF
          ENDIF
          IF(ICASO.EQ.2) THEN  
50        WRITE(*,1045)
          READ(*,*,ERR=50) IQUO
55        WRITE(*,1040)
          READ(*,*,ERR=55) ITABLE  
          IF(ITIPO.EQ.1) THEN  
60        WRITE(*,1010)
          READ(*,*,ERR=60) STP1
65        WRITE(*,1065)
          READ(*,*,ERR=65) SCX,SCY,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
70        WRITE(*,1015)
          READ(*,*,ERR=70) STP1,STP2
75        WRITE(*,1070)
          READ(*,*,ERR=75) SCX,SCY,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
80        WRITE(*,1020)
          READ(*,*,ERR=80) STP1,STP2,STP3  
85        WRITE(*,1075)
          READ(*,*,ERR=85) SCX,SCY,SCU,SCV,SCW 
          ENDIF
          ENDIF
          WRITE(*,1080)
          READ(*,1005) FORMA
          IF(ITIPO.GE.2) THEN  
          WRITE(*,1085)
          READ(*,1005) FORMA2  
          ENDIF
          IF(ITIPO.EQ.3) THEN  
          WRITE(*,1090)
          READ(*,1005) FORMA3  
          ENDIF
          WRITE(*,1095)
          RETURN
1000      FORMAT(5X,'TITOLO')  
1005      FORMAT(20A4) 
1010      FORMAT(5X,'PASSO:',/15X,'(UN VALORE)')
1015      FORMAT(5X,'PASSO:',/15X,'(DUE VALORI)')  
1020      FORMAT(5X,'PASSO:',/15X,'(TRE VALORI)')  
1025      FORMAT(5X,'TIPO:',/15X,'UNA COMPONENTE = 1', 
     ./15X,'DUE COMPONENTI = 2',
     ./15X,'TRE COMPONENTI = 3')
1030      FORMAT(5X,'TOLL:',/15X,'VALORE CRITICO PER LA REIEZIONE ',
     .'DEGLI ERRORI GROSSOLANI')
1035      FORMAT(5X,'IFILE:',/15X,'SALVA SCARTI-RESIDUI E ',
     .'COEFFICIENTI F. SPLINES = 1',/15X,'(DATI OUTPUT FILES 2. E 3.)',
     ./15X,'NON SALVA SCARTI-RESIDUI E COEFFICIENTI F. SPLINES = 0')
1040      FORMAT(5X,'TABLE:'/15X,'SCRIVE TABELLE'  
     .' INCIDENZA E CONNESSIONE = 1'/15X,'NON SCRIVE ' 
     .'TABELLE INCIDENZA E CONNESSIONE = 0')
1045      FORMAT(5X,'IQUO:',/15X,'PRESENZA QUOTE = 1', 
     ./15X,'ASSENZA QUOTE  = 0')
1050      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALA COMPONENTE X')
1055      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALE COMPONENTI X Y') 
1060      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALE COMPONENTI X Y Z')
1065      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE PIANE',
     ./15X,'SCALA COMPONENTE X')
1070      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE PIANE',
     ./15X,'SCALE COMPONENTI X Y') 
1075      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE PIANE',
     ./15X,'SCALE COMPONENTI X Y Z')
1080      FORMAT(5X,'FORMATO DATI INPUT (I COMPONENTE)')
1085      FORMAT(5X,'FORMATO DATI INPUT (II COMPONENTE)')  
1090      FORMAT(5X,'FORMATO DATI INPUT (III COMPONENTE)') 
1095      FORMAT(5X,'DATI INPUT FILE 1.')  
          END  
          SUBROUTINE POLINP(X,Y,W,A,B,P,C,D,S,IGRA,NGRA,TOLL,  
     .    MMAX,NMAX,SCX,SCY,SCW,K,L,FORMA,IFILE)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          CHARACTER*1 DATO,BLK,PTO,AST,SGN 
          CHARACTER*4 FORMA(20)
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),W(1),A(MMAX,1),B(1),P(1),C(NMAX,1),D(1), 
     .    S(1),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA BLK/' '/,PTO/'.'/,AST/'*'/  
          REWIND 1 
          I=0  
10        READ(1,500) (DATO(J),J=1,120)
C          DECODE(150,FORMA,DATO) NPT,XS,YS,WS  
          READ(CDATO,FORMA) NPT,XS,YS,WS
          IF(NPT.EQ.0) GOTO 20 
          IF(DATO(K).EQ.PTO.AND.DATO(L).NE.AST) THEN
          I=I+1
          X(I)=XS*SCX  
          Y(I)=YS*SCY  
          W(I)=WS*SCW  
          ENDIF
          GOTO 10  
20        M=I  
          CALL MEDVAR(X,EX,S2X,SX,M)
          CALL MEDVAR(Y,EY,S2Y,SY,M)
          DO 30 I=1,M  
          X(I)=X(I)-EX 
          Y(I)=Y(I)-EY 
30        CONTINUE 
          IF(IGRA.EQ.0) THEN
          NOTT=0
          SQMW=1.D30  
          DO 60 L=1,NGRA
          DO 50 K=1,M  
          N=1  
          A(K,1)=1.D0  
          IF(L.NE.1) THEN  
          DO 40 I=1,L-1
          N=N+1
          A(K,N)=X(K)**I
          IF(I.NE.1) THEN  
          DO 35 J=1,I-1
          N=N+1
          A(K,N)=X(K)**J*Y(K)**(I-J)
35        CONTINUE 
          ENDIF
          N=N+1
          A(K,N)=Y(K)**I
40        CONTINUE 
          ENDIF
          B(K)=-W(K)
          P(K)=1.D0
50        CONTINUE 
          CALL NORMAL(A,B,P,C,D,MMAX,NMAX,M,N) 
          CALL TCHOL(C,NMAX,N) 
          CALL SOLS(C,D,NMAX,N)
          CALL SZERO(A,B,P,D,SQM,SQM2,MMAX,M,N,0)  
          IF(SQM.LT.SQMW) THEN 
          NOTT=L
          SQMW=SQM 
          ENDIF
60        CONTINUE 
          L=NOTT
          ELSE 
          L=NGRA
          ENDIF
          DO 80 K=1,M  
          N=1  
          A(K,1)=1.D0  
          IF(L.NE.1) THEN  
          DO 70 I=1,L-1
          N=N+1
          A(K,N)=X(K)**I
          IF(I.NE.1) THEN  
          DO 65 J=1,I-1
          N=N+1
          A(K,N)=X(K)**J*Y(K)**(I-J)
65        CONTINUE 
          ENDIF
          N=N+1
          A(K,N)=Y(K)**I
70        CONTINUE 
          ENDIF
          B(K)=-W(K)
          P(K)=1.D0
80        CONTINUE 
          CALL NORMAL(A,B,P,C,D,MMAX,NMAX,M,N) 
          F1=FMAT(IMAX,C,S,NMAX,N) 
          CALL TCHOL(C,NMAX,N) 
          CALL SOLS(C,D,NMAX,N)
          CALL SZERO(A,B,P,D,SQM,SQM2,MMAX,M,N,0)  
          CALL INVER(C,S,NMAX,N)
          F2=FMAT(IMAX,C,S,NMAX,N) 
          CALL DQVV(A,C,P,MMAX,NMAX,M,N)
          RAP=F1/F2
          WRITE(6,690) RAP,IMAX
          DO 90 I=1,N  
          S(I)=0.D0
          IF(SQM2*C(I,I).GT.1.D-10) S(I)=DSQRT(SQM2*C(I,I))
90        CONTINUE 
          WRITE(6,600) 
          WRITE(6,610) 
          N=0  
          DO 110 I=1,L 
          DO 110 J=1,I 
          N=N+1
          IF(S(N).GT.0.D0) THEN 
          ZS=DABS(D(N)/S(N))
          ELSE
          ZS=0.D0
          ENDIF
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,620) D(N),S(N),I-1,J-1,I-J,ZS,SGN
110       CONTINUE 
          IF(IFILE.EQ.1) THEN  
          WRITE(3,630) EX,EY,L-1
          WRITE(3,640) (D(I),I=1,N)
          ENDIF
          WRITE(6,650) 
          DO 120 I=1,M 
          QQ=0.D0  
          IF(SQM2*P(I).GT.1.D-10) QQ=DSQRT(SQM2*P(I))  
          VAL=W(I)+B(I)
          SRV=1.D0 
          IF(SQM2*(1.D0-P(I)).GT.1.D-10) SRV=DSQRT(SQM2*(1.D0-P(I)))
          VINF=VAL-TOLL*SRV
          VSUP=VAL+TOLL*SRV
          ZS=DABS(B(I))/QQ 
          SGN=BLK  
          IF(ZS.GT.TOLL) SGN=AST
          WRITE(6,660) I,W(I),VAL,SRV,B(I),QQ,SGN,VINF,VSUP
          W(I)=B(I)/SCW
120       CONTINUE 
          WRITE(6,670) M,N,L-1,SQM 
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,680) ALFA
          RETURN
500       FORMAT(120A1)
600       FORMAT(///5X,'INTERPOLAZIONE POLINOMIALE')
610       FORMAT(///6X,'PARAMETRI',7X,'SQM',5X,'GRADI',7X,'Z N(0,1)'/) 
620       FORMAT(5X,2D10.3,4X,I2,2X,I2,2X,I2,3X,F6.3,1X,A1)
630       FORMAT(2D20.12,I5)
640       FORMAT(6D20.12)  
650       FORMAT(///,8X,'N.',1X,'VALORE',6X,'VALORE',8X,'SQM', 
     .    4X,'SCARTO',7X,'SQM',5X,'INTERVALLI FIDUCIARI'/  
     .    11X,'OSSERVATO',3X,'STIMATO'/)
660       FORMAT(5X,I5,5F10.3,2X,A1,2X,2F10.3) 
670       FORMAT(///5X,'N. DI OSSERVAZIONI =',I5/  
     .    5X,'N. DI PARAMETRI    =',I5//
     .    5X,'GRADO =',I2,10X,'SIGMA ZERO =',F9.3) 
680       FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEI TESTS: ', 
     .    'ALFA =',F5.3,' (SU DUE CODE)')  
690       FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
          END  
          SUBROUTINE POLINT(Z,W,A,B,P,C,D,S,IGRA,NGRA,TOLL,
     .    MMAX,NMAX,SCZ,SCW,K,L,FORMA,IFILE)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          CHARACTER*1 DATO,BLK,PTO,AST,SGN 
          CHARACTER*4 FORMA(20)
          CHARACTER*120 CDATO
          DIMENSION Z(1),W(1),A(MMAX,1),B(1),P(1),C(NMAX,1),D(1),  
     .    S(1),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA BLK/' '/,PTO/'.'/,AST/'*'/  
          REWIND 1 
          I=0  
10        READ(1,500) (DATO(J),J=1,120)
C          DECODE(150,FORMA,DATO) NPT,TS,WS 
          READ(CDATO,FORMA) NPT,TS,WS
          IF(NPT.EQ.0) GOTO 20 
          IF(DATO(K).EQ.PTO.AND.DATO(L).NE.AST) THEN
          I=I+1
          Z(I)=TS*SCZ  
          W(I)=WS*SCW  
          ENDIF
          GOTO 10  
20        M=I  
          CALL MEDVAR(Z,EZ,S2Z,SZ,M)
          DO 30 I=1,M  
          Z(I)=Z(I)-EZ 
30        CONTINUE 
          IF(IGRA.EQ.0) THEN
          NOTT=0
          SQMW=1.D30  
          DO 60 N=1,NGRA
          DO 50 K=1,M  
          A(K,1)=1.D0  
          IF(N.NE.1) THEN  
          DO 40 I=1,N-1
          A(K,I+1)=Z(K)**I 
40        CONTINUE 
          ENDIF
          B(K)=-W(K)
          P(K)=1.D0
50        CONTINUE 
          CALL NORMAL(A,B,P,C,D,MMAX,NMAX,M,N) 
          CALL TCHOL(C,NMAX,N) 
          CALL SOLS(C,D,NMAX,N)
          CALL SZERO(A,B,P,D,SQM,SQM2,MMAX,M,N,0)  
          IF(SQM.LT.SQMW) THEN 
          NOTT=N
          SQMW=SQM 
          ENDIF
60        CONTINUE 
          N=NOTT
          ELSE 
          N=NGRA
          ENDIF
          DO 80 K=1,M  
          A(K,1)=1.D0  
          IF(N.NE.1) THEN  
          DO 70 I=1,N-1
          A(K,I+1)=Z(K)**I 
70        CONTINUE 
          ENDIF
          B(K)=-W(K)
          P(K)=1.D0
80        CONTINUE 
          CALL NORMAL(A,B,P,C,D,MMAX,NMAX,M,N) 
          F1=FMAT(IMAX,C,S,NMAX,N) 
          CALL TCHOL(C,NMAX,N) 
          CALL SOLS(C,D,NMAX,N)
          CALL SZERO(A,B,P,D,SQM,SQM2,MMAX,M,N,0)  
          CALL INVER(C,S,NMAX,N)
          F2=FMAT(IMAX,C,S,NMAX,N) 
          CALL DQVV(A,C,P,MMAX,NMAX,M,N)
          RAP=F1/F2
          WRITE(6,690) RAP,IMAX
          DO 90 I=1,N  
          S(I)=0.D0
          IF(SQM2*C(I,I).GT.1.D-10) S(I)=DSQRT(SQM2*C(I,I))
90        CONTINUE 
          WRITE(6,600) 
          WRITE(6,610) 
          DO 110 I=1,N 
          IF(S(I).GT.0.D0) ZS=DABS(D(I)/S(I))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,620) D(I),S(I),I-1,ZS,SGN
110       CONTINUE 
          IF(IFILE.EQ.1) THEN  
          WRITE(3,630) EZ,N-1  
          WRITE(3,640) (D(I),I=1,N)
          ENDIF
          WRITE(6,650) 
          DO 120 I=1,M 
          QQ=0.D0  
          IF(SQM2*P(I).GT.1.D-10) QQ=DSQRT(SQM2*P(I))  
          VAL=W(I)+B(I)
          SRV=1.D0 
          IF(SQM2*(1.D0-P(I)).GT.1.D-10) SRV=DSQRT(SQM2*(1.D0-P(I)))
          VINF=VAL-TOLL*SRV
          VSUP=VAL+TOLL*SRV
          IF(QQ.NE.0.D0) ZS=DABS(B(I))/QQ 
          SGN=BLK  
          IF(ZS.GT.TOLL) SGN=AST
          WRITE(6,660) I,W(I),VAL,SRV,B(I),QQ,SGN,VINF,VSUP
          W(I)=B(I)/SCW
120       CONTINUE 
          WRITE(6,670) M,N,N-1,SQM 
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,680) ALFA
          RETURN
500       FORMAT(120A1)
600       FORMAT(///5X,'INTERPOLAZIONE POLINOMIALE')
610       FORMAT(///6X,'PARAMETRI',7X,'SQM',5X,'GRADO',6X,'Z N(0,1)'/) 
620       FORMAT(5X,2D10.3,4X,I2,6X,F10.3,1X,A1)
630       FORMAT(D20.12,I5)
640       FORMAT(6D20.12)  
650       FORMAT(///,8X,'N.',1X,'VALORE',6X,'VALORE',8X,'SQM', 
     .    4X,'SCARTO',7X,'SQM',5X,'INTERVALLI FIDUCIARI'/  
     .    11X,'OSSERVATO',3X,'STIMATO'/)
660       FORMAT(5X,I5,5F10.3,2X,A1,2X,2F10.3) 
670       FORMAT(///5X,'N. DI OSSERVAZIONI =',I5/  
     .    5X,'N. DI PARAMETRI    =',I5//
     .    5X,'GRADO =',I2,10X,'SIGMA ZERO =',F9.3) 
680       FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEI TESTS: ', 
     .    'ALFA =',F5.3,' (SU DUE CODE)')  
690       FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
          END
          SUBROUTINE PRECLP(NS,XS,YS,X,Y,Z,C,XCOV,INDIR,
     .IPX,IPY,LPX,LPY,IWD,NTIPO,KA,NMAX,NWDMAX,NWD,
     .N,DA,DB,DC,SN,NCASO,U)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION Z(1),C(3,1),X(1),Y(1),IPX(1),IPY(1),
     .LPX(1),LPY(1),IWD(1),NTIPO(1),KA(1),INDIR(1),
     .DA(1),DB(1),DC(1),SN(1),NCASO(1),U(1),INV(3),IND(7)  
          DATA INV/1,3,7/,IND/1,1,2,1,2,2,3/
          DO 1 I=1,NS  
          U(I)=0.D0
1         CONTINUE 
          AMPX=4.D0*XCOV
          AMPY=4.D0*XCOV
          CALL AWINDW(XS,YS,AMPX,AMPY,X,Y,IPX,IPY,LPX, 
     .LPY,IWD,NMAX,NWDMAX,N,NWD)
          II=INV(NS)
          IR=NS
          DO 10 I=1,NWD
          IS=IWD(I)
          JS=INDIR(IS) 
          DELTA=DSQRT((XS-X(JS))**2+(YS-Y(JS))**2) 
          JJ=NTIPO(JS) 
          IC=IND(JJ)
          INI=KA(JS)+1 
          DO 10 L=1,IR 
          DO 10 K=1,IC 
          C(L,INI+K-1)=GASS2(DELTA,XCOV,II,JJ,L,K,DA,DB,
     .DC,SN,NCASO) 
10        CONTINUE 
          DO 20 J=1,NWD
          IS=IWD(J)
          JS=INDIR(IS) 
          JJ=NTIPO(JS) 
          IC=IND(JJ)
          INI=KA(JS)+1 
          DO 20 I=1,IR 
          DO 20 K=1,IC 
          U(I)=U(I)+Z(INI+K-1)*C(I,INI+K-1)
          C(I,INI+K-1)=0.D0
20        CONTINUE 
          RETURN
          END  
          SUBROUTINE PRECLT(NS,TP,XCOVT,Z,TEMPO,C,IPUNT,
     .NTIPO,KA,NMAX,N,DA,DB,DC,SN,NCASO,U) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION Z(1),C(3,1),TEMPO(1),NTIPO(1),KA(1),U(1),IPUNT(1), 
     .DA(1),DB(1),DC(1),SN(1),NCASO(1),INV(3),IND(7)
          DATA INV/1,3,7/,IND/1,1,2,1,2,2,3/
          DO 1 I=1,NS  
          U(I)=0.D0
1         CONTINUE 
          TMI=TP-2.D0*XCOVT
          TMA=TP+2.D0*XCOVT
          CALL ABINRY(TMI,TEMPO,IPUNT,NMAX,LNI,LNS,LPI,LPS,N)  
          LI=LNS
          CALL ABINRY(TMA,TEMPO,IPUNT,NMAX,LNI,LNS,LPI,LPS,N)  
          LM=LNI
          INI=KA(LI)+1 
          IFI=KA(LM+1) 
          II=INV(NS)
          IR=NS
          INC=INI  
          DO 10 I=LI,LM
          JJ=NTIPO(I)  
          IC=IND(JJ)
          DELTA=DABS(TP-TEMPO(I))  
          DO 5 L=1,IR  
          DO 5 K=1,IC  
          C(L,INC+K-1)=GASS1(DELTA,XCOVT,II,JJ,L,K,DA,DB,DC,SN,NCASO)  
5         CONTINUE 
          INC=INC+IC
10        CONTINUE 
          DO 20 I=1,IR 
          DO 20 J=INI,IFI  
          U(I)=U(I)+Z(J)*C(I,J)
          C(I,J)=0.D0  
20        CONTINUE 
          RETURN
          END  
          SUBROUTINE PREPLP(NS,XS,YS,X0,Y0,COEF,NGR,U) 
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          DIMENSION COEF(3,1),X0(1),Y0(1),NGR(1),U(1)  
          DO 1 I=1,NS  
          U(I)=0.D0
1         CONTINUE 
          DO 10 I=1,NS 
          X=XS-X0(I)
          Y=YS-Y0(I)
          N=NGR(I)+1
          K=1  
          U(I)=COEF(I,1)
          IF(N.EQ.1) GO TO 10  
          DO 8 J=1,N-1 
          K=K+1
          U(I)=U(I)+COEF(I,K)*X**J 
          IF(J.EQ.1) GO TO 7
          DO 5 L=1,J-1 
          K=K+1
          U(I)=U(I)+COEF(I,K)*X**L*Y**(J-L)
5         CONTINUE 
7         K=K+1
          U(I)=U(I)+COEF(I,K)*Y**J 
8         CONTINUE 
10        CONTINUE 
          RETURN
          END  
          SUBROUTINE PREPLT(NS,TP,TT,COEF,NGR,U)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION COEF(3,1),NGR(1),U(1),TT(1)
          DO 5 I=1,NS  
          U(I)=0.D0
5         CONTINUE 
          DO 10 I=1,NS 
          T=TP-TT(I)
          N=NGR(I)+1
          U(I)=COEF(I,1)
          IF(N.EQ.1) GO TO 10  
          DO 8 J=1,N-1 
          U(I)=U(I)+COEF(I,J+1)*T**J
8         CONTINUE 
10        CONTINUE 
          RETURN
          END  
          SUBROUTINE PRESPP(NS,XS,YS,X0,Y0,STP,COEF,NI,U)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION COEF(3,1),X0(1),Y0(1),STP(1),NI(1),U(1)
          DO 1 I=1,NS  
          U(I)=0.D0
1         CONTINUE 
          DO 20 I=1,NS 
          X=XS-X0(I)
          Y=YS-Y0(I)
          ST=STP(I)
          NIX=NI(I)
          NX=IDINT(X/ST)+1 
          NY=IDINT(Y/ST)+1 
          NC=(NY-1)*(NIX-1)+NX 
          SRX=(X-(NX-1)*ST)/ST 
          SRY=(Y-(NY-1)*ST)/ST 
          NIC=(NY-2)*NIX+NX-1  
          DO 10 J=1,4  
          NII=NIC+(J-1)*NIX
          AY=DBLE(2-J) 
          JY=5-J
          VALY=SPLT(SRY+AY,JY) 
          DO 10 K=1,4  
          NJ=NII+K-1
          AX=DBLE(2-K) 
          JX=5-K
          U(I)=U(I)+COEF(I,NJ)*VALY*SPLT(SRX+AX,JX)
10        CONTINUE 
20        CONTINUE 
          RETURN
          END  
          SUBROUTINE PRESPT(NS,TP,TT,STP,COEF,U)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION COEF(3,1),STP(1),U(1),TT(1)
          DO 5 I=1,NS  
          U(I)=0.D0
5         CONTINUE 
          DO 10 I=1,NS 
          T=TP-TT(I)
          NC=IDINT(T/STP(I))
          TS=NC*STP(I) 
          SRV=(T-TS)/STP(I)
          DO 10 J=1,4  
          AJ=DBLE(2-J) 
          JS=5-J
          NI=NC+J-1
          U(I)=U(I)+COEF(I,NI)*SPLT(SRV+AJ,JS) 
10        CONTINUE 
          RETURN
          END  
          SUBROUTINE RCORR2(X,Y,LPX,LPY,EX,EY,SX,SY,S2X,S2Y,SXY,RXY,NN)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),LPX(1),LPY(1)
          NXY=0
          SRVX=0.D0
          SRVY=0.D0
          SRV2X=0.D0
          SRV2Y=0.D0
          SRVXY=0.D0
          DO 1 I=1,NN  
          IF(LPX(I).EQ.0.OR.LPY(I).EQ.0) GOTO 1
          NXY=NXY+1
          SRVX=SRVX+X(I)
          SRVY=SRVY+Y(I)
          SRV2X=SRV2X+X(I)**2  
          SRV2Y=SRV2Y+Y(I)**2  
          SRVXY=SRVXY+X(I)*Y(I)
1         CONTINUE 
          EX=SRVX/NXY  
          EY=SRVY/NXY  
          S2X=(SRV2X/NXY-EX**2)*NXY/(NXY-1)
          S2Y=(SRV2Y/NXY-EY**2)*NXY/(NXY-1)
          SX=DSQRT(S2X)
          SY=DSQRT(S2Y)
          SXY=SRVXY/NXY-EX*EY  
          RXY=SXY/(SX*SY)  
          RETURN
          END  

          FUNCTION SCAL1(Z)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          IF(Z.LT.2.D0) THEN
          R=1.D0
          Y1=1.0667D0  
          A=(16.D0/15.D0)*R**5 
          B=(4.D0/3.D0)*R**3
          C=(2.D0/3.D0)*R**2
          D=(1.D0/30.D0)
          Y=A-B*Z**2+C*Z**3-D*Z**5 
          SCAL1=Y/Y1
          ELSE 
          SCAL1=0.D0
          ENDIF
          RETURN
          END  
          FUNCTION SCAL2(Z)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          IF(Z.LT.2.D0) THEN
          PIG=4.D0*DATAN(1.D0) 
          R=1.D0
          Y1=5.D0*R**6*PIG 
          A=5*R**6*PIG 
          B=(15.D0/2.D0)*R**4*Z**2*PIG 
          SQ=DSQRT(DABS(R**2-(Z/2.D0)**2)) 
          AS=DASIN(Z/(2.D0*R)) 
          C=5.D0*R**4*Z+(20.D0/3.D0)*R**2*Z**3-(5.D0/12.D0)*Z**5
          D=15.D0*R**4*Z**2-10.D0*R**6 
          Y=A-B+C*SQ+D*AS  
          SCAL2=Y/Y1
          ELSE 
          SCAL2=0.D0
          ENDIF
          RETURN
          END  
          SUBROUTINE SPARSE(XCOVT,X,Y,IPX,IPY,LPX,LPY,IWD,IPUNT,
     .JCRONA,JCOVAS,IBANDA,JBANDA,NUMERO,INDIR,LC,JC,  
     .NMAX,NWDMAX,N,NWD,IORDI) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),IPX(1),IPY(1),LPX(1),LPY(1),IWD(1),  
     .IPUNT(1),JCRONA(1),JCOVAS(1),IBANDA(1),JBANDA(1),
     .NUMERO(1),INDIR(1),LC(1),JC(1)
          CALL BIORDA(X,Y,IPX,IPY,LPX,LPY,N)
          XCOVT2=2.D0*XCOVT
          AMPX=4.D0*XCOVT  
          AMPY=4.D0*XCOVT  
          JCRONA(1)=1  
          K=0  
          DO 20 I=1,N  
          XC=X(I)  
          YC=Y(I)  
          CALL AWINDW(XC,YC,AMPX,AMPY,X,Y,IPX,IPY,LPX,LPY,IWD, 
     .NMAX,NWDMAX,N,NWD)
          CALL ORDINI(IWD,IPUNT,NWD,1)
          DO 10 J=1,NWD
          JJ=IPUNT(J)  
          II=IWD(JJ)
          DIST=DSQRT((XC-X(II))**2+(YC-Y(II))**2)  
          IF(DIST.LE.XCOVT2) THEN  
          K=K+1
          JCOVAS(K)=II 
          ENDIF
10        CONTINUE 
          JCRONA(I+1)=K+1  
20        CONTINUE 
          IF(IORDI.EQ.1) THEN  
          N1=1+NMAX
          N2=N1+NMAX
          N3=N2+NMAX
          N4=N3+NMAX
          N5=N4+NMAX
          N6=N5+NMAX
          CALL GIBBS(NUMERO,INDIR,JCRONA,JCOVAS,JC,JC(N1),JC(N2),  
     .JC(N3),JC(N4),JC(N5),JC(N6),N)
          LC(1)=1  
          K=0  
          DO 50 I=1,N  
          IS=INDIR(I)  
          J1=JCRONA(IS)
          J2=JCRONA(IS+1)-1
          NWD=J2-J1+1  
          DO 30 J=1,NWD
          JS=JCOVAS(J1-1+J)
          IWD(J)=NUMERO(JS)
30        CONTINUE 
          CALL ORDINI(IWD,IPUNT,NWD,1) 
          DO 40 J=1,NWD
          K=K+1
          JS=IPUNT(J)  
          JC(K)=IWD(JS)
40        CONTINUE 
          LC(I+1)=LC(I)+NWD
50        CONTINUE 
          DO 60 I=1,N+1
          JCRONA(I)=LC(I)  
60        CONTINUE 
          NNWD=JCRONA(N+1)-1
          DO 70 I=1,NNWD
          JCOVAS(I)=JC(I)  
70        CONTINUE 
          ELSE 
          DO 80 I=1,N  
          NUMERO(I)=I  
          INDIR(I)=I
80        CONTINUE 
          ENDIF
          DO 90 I=1,N  
          JBANDA(I)=0  
90        CONTINUE 
          DO 100 I=1,N 
          J1=JCRONA(I) 
          J2=JCRONA(I+1)-1 
          IBANDA(I)=I-JCOVAS(J1)
          DO 100 JS=J1,J2  
          J=JCOVAS(JS) 
          JBANDA(J)=MAX0(JBANDA(J),I-J)
100       CONTINUE 
          RETURN
          END  
          SUBROUTINE SPLINP(A,IA,P,B,KA,JA,LA,IB,JB,
     .C,T,JC,LC,D,W,R,U,ISRV,NNG,FORMA,STP,
     .MMAX,NMAX,SCX,SCY,SCU,KK,LL,TOLL,ITABLE,IFILE)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,AST,PTO,ASTI
          CHARACTER*4 FORMA(20),FF1(8),FF2(7),FF3(5)
          CHARACTER*120 CDATO
          DIMENSION A(MMAX,1),IA(MMAX,1),KA(1),B(1),
     .P(1),JA(1),LA(1),IB(1),JB(1),
     .C(1),T(1),JC(1),LC(1),D(1),  
     .W(1),R(1),U(1),ISRV(1),NNG(1),NI(4), 
     .DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA AST/'*'/,PTO/'.'/,ASTI/' '/
          DATA FF1/'(I10',',2F1','5.3,','F15.','6,F1','0.3,',
     .'    ','A1) '/
          DATA FF2/'(I10',',2F1','5.3,','F15.','    ','X,A1','   )'/
          DATA FF3/'(I10',',2F1','5.3,','    ','A1) '/
          REWIND 1 
          IF(LL.EQ.116) THEN
          FF1(7)='50X,'
          FF2(5)='6,60'
          FF3(4)='75X,'
          ENDIF
          IF(LL.EQ.118) THEN
          FF1(7)='52X,'
          FF2(5)='6,62'
          FF3(4)='77X,'
          ENDIF
          IF(LL.EQ.120) THEN
          FF1(7)='54X,'
          FF2(5)='6,64'
          FF3(4)='79X,'
          ENDIF
          DO 5 I=1,NMAX
          NNG(I)=0 
5         CONTINUE 
C         CALL CRONO('INPUT  .')
          X0=1.D30
          Y0=1.D30
          EU=0.D0
          I=0  
10        READ(1,1010) (DATO(L),L=1,120)
C          DECODE(150,FORMA,DATO) NPT,X,Y,US
          READ(CDATO,FORMA) NPT,X,Y,US
          IF(NPT.EQ.0) GOTO 15 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          EU=EU+US
          X0=DMIN1(X,X0)
          Y0=DMIN1(Y,Y0)
          I=I+1
          ENDIF
          GOTO 10  
15        M=I  
          EU=EU/M
          IF(M.GT.MMAX) THEN
          WRITE(6,1190)
          STOP 
          ENDIF
          REWIND 1 
          NX=IDINT(X0) 
          NY=IDINT(Y0) 
          X0=DBLE(NX)-STP-STP/1000.D0  
          Y0=DBLE(NY)-STP-STP/1000.D0  
          NCX=0
          NCY=0
          I=0  
20        READ(1,1010) (DATO(L),L=1,120)
C          DECODE(150,FORMA,DATO) NPT,X,Y,US
          READ(CDATO,FORMA) NPT,X,Y,US
          IF(NPT.EQ.0) GOTO 30 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          I=I+1
          US=US-EU
          P(I)=1.D0
          B(I)=-US*SCU 
          KA(I)=16 
          NX=IDINT((X-X0)*SCX/STP)+1
          NY=IDINT((Y-Y0)*SCY/STP)+1
          IF(NX.GT.NCX) NCX=NX 
          IF(NY.GT.NCY) NCY=NY 
          ENDIF
          GOTO 20  
30        NCX=NCX+1
          NCY=NCY+1
          NTC=NCX*NCY  
          NIX=NCX+1
          WRITE(6,1170) NCX,NCY,NTC,(NCX+1)*(NCY+1),STP,X0,Y0  
          REWIND 1 
          I=0  
40        READ(1,1010) (DATO(L),L=1,120)
C          DECODE(150,FORMA,DATO) NPT,X,Y,US
          READ(CDATO,FORMA) NPT,X,Y,US
          IF(NPT.EQ.0) GOTO 60 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          I=I+1
          X=(X-X0)*SCX 
          Y=(Y-Y0)*SCY 
          NX=IDINT(X/STP)+1
          NY=IDINT(Y/STP)+1
          NC=(NY-1)*NCX+NX 
          NNG(NC)=NNG(NC)+1
          SRVX=X-(NX-1)*STP
          SRVY=Y-(NY-1)*STP
          SRVX=SRVX/STP
          SRVY=SRVY/STP
          NI(1)=(NY-2)*NIX+NX-1
          NI(2)=NI(1)+NIX  
          NI(3)=NI(2)+NIX  
          NI(4)=NI(3)+NIX  
          DO 50 J=1,4  
          AY=DBLE(2-J) 
          JY=5-J
          VALY=SPLT(SRVY+AY,JY)
          DO 50 K=1,4  
          JS=(J-1)*4+K 
          AX=DBLE(2-K) 
          JX=5-K
          A(I,JS)=VALY*SPLT(SRVX+AX,JX)
          IA(I,JS)=NI(J)+K-1
50        CONTINUE 
          ENDIF
          GOTO 40  
60        N=(NCX+1)*(NCY+1)
          WRITE(6,1110)
          DO 65 K=1,NTC
          J=(K-1)/NCX  
          WRITE(6,1115) K,NNG(K),K+J,K+J+1,K+NIX+J,K+NIX+J+1
65        CONTINUE 
          CALL NORMAI(A,IA,KA,JA,LA,IB,JB,B,P,C,JC,LC,D,W,
     .ISRV,MMAX,M,N)
          DO 70 K=1,N  
          ISRV(K)=K
70        CONTINUE 
          IF(ITABLE.EQ.1) CALL WGRAFO(ISRV,JA,LA,IB,JB,N,0)
          IF(ITABLE.EQ.1) CALL WGRAFO(ISRV,JC,LC,IB,JB,N,1)
          NN=LC(N+1)-1 
          WRITE(6,1020) M,N,NN/N,N*(N+1)/2,NN,N*(N+1)/(2*NN)
          NV=0 
          DO 80 I=1,NIX
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
80        CONTINUE 
          DO 90 I=NIX+1,N-NIX+1,NIX
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
90        CONTINUE 
          DO 100 I=2*NIX,N,NIX 
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
100       CONTINUE 
          DO 110 I= NCY*NIX+2,N-1  
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
110       CONTINUE 
          DO 120 I=1,NCY-1 
          NII=(NIX+2)+(I-1)*NIX
          NIF=NII+NIX-3
          DO 120 J=NII,NIF 
          N1=J-NIX-I
          N2=N1+1  
          N3=N1+NCX
          N4=N3+1  
c	  NPN=NNG(N1)+NNG(N2)+NNG(N3)+NNG(N4)
c	  IF(NPN.LT.9) THEN
          IF(NNG(N1).EQ.0.AND.NNG(N2).EQ.0.AND.
     .NNG(N3).EQ.0.AND.NNG(N4).EQ.0) THEN  
          NV=NV+1  
          IS=LC(J+1)-1 
          C(IS)=C(IS)+1.D10
          ENDIF
120       CONTINUE 
C         CALL CRONO('SISTEMA.')
          F1=FMATI1(IMAX,C,JC,LC,W,N)  
          CALL ALGITE(C,T,JC,LC,D, 
     .JA,LA,IB,JB, 
     .1.D-5,100,ITE,1,1,R,U,W,N)
          F2=FMATI2(IMAX,T,JA,JC,LA,W,N)
          RAP=F1/F2
          WRITE(6,1125) ITE
          WRITE(6,1130) RAP,IMAX
          CALL SZEROR(A,IA,KA,B,P,D,SZ,SZ2,MMAX,M,N,NV)
          CALL DQVVI(A,IA,KA,T,JA,JC,LA,P,W,MMAX,M)
C         CALL CRONO('ALGEBRA.')
          REWIND 1 
          WRITE(6,1030)
          DO 130 I=1,N 
          IS=LC(I+1)-1 
          IF(SZ2*C(IS).GT.1.D-10) C(IS)=DSQRT(SZ2*C(IS))
          WRITE(6,1040) I,D(I),C(IS)
130       CONTINUE 
          IF(IFILE.EQ.1) THEN  
          WRITE(3,1015) X0,Y0,STP,EU,N,NIX
          WRITE(3,1018) (D(I),I=1,N)
          ENDIF
          WRITE(6,1050)
          REWIND 1 
          I=0  
140       READ(1,1010) (DATO(L),L=1,120)
C          DECODE(150,FORMA,DATO) NPT,X,Y,US
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
          WRITE(6,1060) NPT,X,Y,US,SGNL,SQMS,B(I),SQMV,ASTI
          IF(IFILE.EQ.1) WRITE(2,FF1) NPT,X,Y,B(I),SQMV,ASTI  
          ELSE
          IF(IFILE.EQ.1) THEN
          IF(DATO(KK).EQ.PTO) WRITE(2,FF2) NPT,X,Y,US,DATO(LL)
          IF(DATO(KK).NE.PTO) WRITE(2,FF3) NPT,X,Y,DATO(LL)
          ENDIF
          ENDIF
          GOTO 140 
150       WRITE(2,1000)
          NGL=M-N+NV
          WRITE(6,1070) SZ,M,N,NV,NGL  
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,1108) ALFA
          WRITE(6,1005)
C         CALL CRONO('OUTPUT.')
C         CALL PCRONO  
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
1030      FORMAT(///4X,'NUMERO NODO',9X,'COEFFICENTE',7X,'SQM',/)  
1040      FORMAT(5X,I10,5X,F15.6,F10.3)
1050      FORMAT(///5X,'NOME PUNTO',3X,'COORDINATA X',3X,  
     .'COORDINATA Y',  
     .5X,'COMPONENTE',5X,'COMPONENTE',7X,'SQM',8X,'RESIDUO',7X,'SQM'
     ./65X,'STIMATA',/)
1060      FORMAT(5X,I10,2F15.3,F15.6,F15.6,F10.3,F15.6,F10.3,5X,A1)
1070      FORMAT(////5X,'SIGMA ZERO  ',F10.3,//
     ./5X,'NUMERO EQUAZIONI  ',I5,/5X,'NUMERO INCOGNITE  ',I5, 
     ./5X,'NUMERO VINCOLI    ',I5,/5X,'GRADI DI LIBERTA  ',I5) 
1108      FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEL TEST '
     .'DI REIEZIONE DEGLI ERRORI GROSSOLANI: ALFA =',F5.3, 
     .' (SU DUE CODE)')
1110      FORMAT(///5X,'CLASSE -',5X,'NODI COINVOLTI DALLA CLASSE'/
     .5X,'NUMEROSITA',3X,'(VERTICI DEL QUADRATO)'/)
1115      FORMAT(5X,2I5,10X,4I5)
1125      FORMAT(//5X,'ITERAZIONI=',I5//)  
1130      FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
1170      FORMAT(///5X,'NUMERO CLASSI IN X E Y    ',2I5,
     .//5X,'NUMERO TOTALE CLASSI      ',I5,
     .//5X,'NUMERO NODI               ',I5,
     .//5X,'PASSO DELLA INTERPOLAZIONE',F10.3, 
     .//5X,'COORDINATE ORIGINE        ',2F12.3,//) 
1190      FORMAT(//5X,'NUMERO DEI DATI',
     ./5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'//)
          END  
          SUBROUTINE SPLINT(A,IA,P,B,KA,C,D,W,NNG,SRW,FORMA,STP,
     .MMAX,NMAX,NBMAX,SCT,SCX,KK,LL,TOLL,IFILE)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          CHARACTER*1 DATO,PTO,AST,ASTI
          CHARACTER*4 FORMA(20),FF1(8),FF2(7),FF3(5)
          CHARACTER*120 CDATO
          DIMENSION A(MMAX,1),IA(MMAX,1),P(1),B(1),KA(1),C(NMAX,1),
     .D(1),W(1),SRW(1),NNG(1),DATO(120)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/,ASTI/' '/ 
          DATA FF1/'(I10',',1F1','5.3,','F15.','6,F1','0.3,',
     .'    ','A1) '/
          DATA FF2/'(I10',',1F1','5.3,','F15.','    ','X,A1','   )'/
          DATA FF3/'(I10',',1F1','5.3,','    ','A1) '/
          REWIND 1 
          IF(LL.EQ.116) THEN
          FF1(7)='65X,'
          FF2(5)='6,75'
          FF3(4)='80X,'
          ENDIF
          IF(LL.EQ.118) THEN
          FF1(7)='67X,'
          FF2(5)='6,77'
          FF3(4)='82X,'
          ENDIF
          IF(LL.EQ.120) THEN
          FF1(7)='69X,'
          FF2(5)='6,79'
          FF3(4)='84X,'
          ENDIF
          REWIND 1 
          DO 5 I=1,NMAX
          NNG(I)=0 
5         CONTINUE 
C         CALL CRONO('INPUT  .')
          TOR=1.D30
          I=0
          EX=0.D0  
10        READ(1,1010) (DATO(L),L=1,120)
C          DECODE(150,FORMA,DATO) NPT,TS,XS 
          READ(CDATO,FORMA) NPT,TS,XS
          IF(NPT.EQ.0) GOTO 15 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          TOR=DMIN1(TS,TOR)
          I=I+1
          EX=EX+XS
          ENDIF
          GOTO 10  
15        M=I
          EX=EX/M
          NO=IDINT(TOR)
          TOR=DBLE(NO)-STP-STP/1000.D0 
          IF(M.GT.MMAX) THEN
          WRITE(6,1190)
          STOP 
          ENDIF
          REWIND 1 
          I=0  
20        READ(1,1010) (DATO(L),L=1,120)
C          DECODE(150,FORMA,DATO) NPT,TS,XS 
          READ(CDATO,FORMA) NPT,TS,XS
          IF(NPT.EQ.0) GOTO 30 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          I=I+1
          XS=XS-EX
          P(I)=1.D0
          B(I)=-XS*SCX 
          KA(I)=4  
          TS=(TS-TOR)*SCT  
          N=IDINT(TS/STP)  
          NNG(N+1)=NNG(N+1)+1  
          TT=N*STP 
          SRV=(TS-TT)/STP  
          DO 25 J=1,4  
          AJ=DBLE(2-J) 
          JS=5-J
          A(I,J)=SPLT(SRV+AJ,JS)
          IA(I,J)=N+J-1
25        CONTINUE 
          ENDIF
          GOTO 20  
30        N=IA(M,4)
          WRITE(6,1130) N-1,N,STP,TOR  
          WRITE(6,1140)
          DO 33 K=1,N-1
          WRITE(6,1150) K,NNG(K),K,K+1 
33        CONTINUE 
C         CALL CRONO('SISTEMA.')
          CALL NORMAB(A,IA,KA,B,P,C,D,MMAX,NMAX,NBMAX,M,N,NB)  
          C(1,1)=C(1,1)+1.D10  
          C(N,1)=C(N,1)+1.D10  
          ISRV=0
          NV=2 
          DO 35 I=2,N-1
          IF(ISRV.EQ.0.AND.NNG(I).EQ.0) THEN
          NV=NV+1  
          C(I,1)=C(I,1)+1.D10  
          ENDIF
          ISRV=NNG(I)  
35        CONTINUE 
          WRITE(6,1020) M,N,NB,N**2,N*NB,N**2/(N*NB)
          F1=FMATB(IMAX,C,SRW,MMAX,M,MB)
          CALL TCHOLB(C,NMAX,N,NB) 
          CALL SOLSB(C,D,NMAX,N,NB)
          CALL SZEROR(A,IA,KA,B,P,D,SZ,SZ2,MMAX,M,N,NV)
          CALL INVERB(C,W,NMAX,N,NB)
          F2=FMATB(IMAX,C,SRW,MMAX,M,MB)
          RAP=F1/F2
          WRITE(6,1120) RAP,IMAX
          CALL DQVVB(A,IA,KA,C,P,MMAX,NMAX,M)  
          DO 40 I=1,N  
          W(I)=0.D0
          IF(SZ2*C(I,1).GT.1.D-10) W(I)=DSQRT(SZ2*C(I,1))  
40        CONTINUE 
C         CALL CRONO('ALGEBRA.')
          WRITE(6,1030)
          DO 60 I=1,N  
          WRITE(6,1040) I,D(I),W(I)
60        CONTINUE 
          IF(IFILE.EQ.1) THEN  
          WRITE(3,1015) TOR,STP,EX,N  
          WRITE(3,1018) (D(I),I=1,N)
          ENDIF
          REWIND 1 
          WRITE(6,1050)
          I=0  
70        READ(1,1010) (DATO(L),L=1,120)
C          DECODE(150,FORMA,DATO) NPT,TS,XS 
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
          WRITE(6,1060) NPT,TS,XS,SGNL,SQMS,B(I),P(I),ASTI 
          IF(IFILE.EQ.1) WRITE(2,FF1) NPT,TS,B(I),P(I),ASTI
          ELSE
          IF(IFILE.EQ.1) THEN
          IF(DATO(KK).EQ.PTO) WRITE(2,FF2) NPT,TS,XS,DATO(LL)
          IF(DATO(KK).NE.PTO) WRITE(2,FF3) NPT,TS,DATO(LL)
          ENDIF
          ENDIF
          GOTO 70  
80        WRITE(2,1000)
          NGL=M-N+NV
          WRITE(6,1070) SZ,M,N,NV,NGL  
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,1100) ALFA
          WRITE(6,1005)
C         CALL CRONO('OUTPUT.')
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
1030      FORMAT(///4X,'NUMERO NODO',9X,'COEFFICENTE',7X,'SQM',/)  
1040      FORMAT(5X,I10,5X,F15.6,F10.3)
1050      FORMAT(///5X,'NOME PUNTO',10X,'TEMPO',5X,'COMPONENTE',5X,
     .'COMPONENTE',7X,'SQM',8X,'RESIDUO',7X,'SQM'  
     ./50X,'STIMATA',/)
1060      FORMAT(5X,I10,F15.3,F15.6,F15.6,F10.3,F15.6,F10.3,20X,A1)
1070      FORMAT(////5X,'SIGMA ZERO  ',F10.3,//
     ./5X,'NUMERO EQUAZIONI  ',I5,/5X,'NUMERO INCOGNITE  ',I5, 
     ./5X,'NUMERO VINCOLI    ',I5,/5X,'GRADI DI LIBERTA  ',I5) 
1080      FORMAT(5X,D15.6) 
1100      FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEL TEST '
     .'DI REIEZIONE DEGLI ERRORI GROSSOLANI: ALFA =',F5.3, 
     .' (SU DUE CODE)')
1120      FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
1130      FORMAT(///5X,'NUMERO CLASSI             ',I5,
     .//5X,'NUMERO NODI               ',I5,
     .//5X,'PASSO DELLA INTERPOLAZIONE',F10.3, 
     .//5X,'TEMPO ORIGINE             ',F10.3,//)  
1140      FORMAT(///5X,'CLASSE -',5X,'NODI COINVOLTI DALLA CLASSE'/
     .5X,'NUMEROSITA',3X,'(ESTREMI DELLO INTERVALLO)'/)
1150      FORMAT(5X,2I5,20X,2I5)
1190      FORMAT(//5X,'NUMERO DEI DATI',
     ./5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'//)
          END  

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
          SUBROUTINE TESTA(ICASO,ITIPO,IQUO,IPTO,NPTO) 
          CHARACTER*1 DATO,PARA,PARC,BLANK 
          DIMENSION DATO(120),IPTO(1)  
          DATA PARA/'('/,PARC/')'/,BLANK/' '/  
          DO 5 I=1,120 
          DATO(I)=BLANK
5         CONTINUE 
          WRITE(6,1000)
          IF(ICASO.EQ.1) THEN  
          IF(ITIPO.EQ.1) WRITE(6,1010) 
          IF(ITIPO.EQ.2) WRITE(6,1020) 
          IF(ITIPO.EQ.3) WRITE(6,1030) 
          ELSE 
          IF(IQUO.EQ.0) THEN
          IF(ITIPO.EQ.1) WRITE(6,1040) 
          IF(ITIPO.EQ.2) WRITE(6,1050) 
          IF(ITIPO.EQ.3) WRITE(6,1060) 
          ELSE 
          IF(ITIPO.EQ.1) WRITE(6,1070) 
          IF(ITIPO.EQ.2) WRITE(6,1080) 
          IF(ITIPO.EQ.3) WRITE(6,1090) 
          ENDIF
          ENDIF
          DATO(3)=PARA 
          DATO(4)='1'  
          DATO(5)=PARC 
          DO 10 I=1,NPTO
          J=IPTO(I)
          DATO(J-1)=PARA
          DATO(J)=CHAR(I+49)
          DATO(J+1)=PARC
10        CONTINUE 
          WRITE(6,1100) (DATO(I),I=1,120)  
          RETURN
1000      FORMAT(///5X,'DATI DI INPUT',///5X,'LEGENDA:',/) 
1010      FORMAT(5X,'(1) NOME PUNTO',/ 
     .5X,'(2) TEMPO',/ 
     .5X,'(3) COMPONENTE X')
1020      FORMAT(5X,'(1) NOME PUNTO',/ 
     .5X,'(2) TEMPO',/ 
     .5X,'(3) COMPONENTE X',/  
     .5X,'(4) COMPONENTE Y')
1030      FORMAT(5X,'(1) NOME PUNTO',/ 
     .5X,'(2) TEMPO',/ 
     .5X,'(3) COMPONENTE X',/  
     .5X,'(4) COMPONENTE Y',/  
     .5X,'(5) COMPONENTE Z')
1040      FORMAT(5X,'(1) NOME PUNTO',/ 
     .5X,'(2) COORDINATA PIANA ASCISSA',/  
     .5X,'(3) COORDINATA PIANA ORDINATA',/ 
     .5X,'(4) COMPONENTE X')
1050      FORMAT(5X,'(1) NOME PUNTO',/ 
     .5X,'(2) COORDINATA PIANA ASCISSA',/  
     .5X,'(3) COORDINATA PIANA ORDINATA',/ 
     .5X,'(4) COMPONENTE X',/  
     .5X,'(5) COMPONENTE Y')
1060      FORMAT(5X,'(1) NOME PUNTO',/ 
     .5X,'(2) COORDINATA PIANA ASCISSA',/  
     .5X,'(3) COORDINATA PIANA ORDINATA',/ 
     .5X,'(4) COMPONENTE X',/  
     .5X,'(5) COMPONENTE Y',/  
     .5X,'(6) COMPONENTE Z')
1070      FORMAT(5X,'(1) NOME PUNTO',/ 
     .5X,'(2) COORDINATA PIANA ASCISSA',/  
     .5X,'(3) COORDINATA PIANA ORDINATA',/ 
     .5X,'(4) QUOTA',/ 
     .5X,'(5) COMPONENTE X')
1080      FORMAT(5X,'(1) NOME PUNTO',/ 
     .5X,'(2) COORDINATA PIANA ASCISSA',/  
     .5X,'(3) COORDINATA PIANA ORDINATA',/ 
     .5X,'(4) QUOTA',/ 
     .5X,'(5) COMPONENTE X',/  
     .5X,'(6) COMPONENTE Y')
1090      FORMAT(5X,'(1) NOME PUNTO',/ 
     .5X,'(2) COORDINATA PIANA ASCISSA',/  
     .5X,'(3) COORDINATA PIANA ORDINATA',/ 
     .5X,'(4) QUOTA',/ 
     .5X,'(5) COMPONENTE X',/  
     .5X,'(6) COMPONENTE Y',/  
     .5X,'(7) COMPONENTE Z')
1100      FORMAT(//5X,120A1/)  
          END 
          SUBROUTINE TRIMP(X,Y,U,IPX,IPY,LPX,LPY,IWD,SNOISE,OLATO, 
     .    IND,N,NN,NM,NWDMAX,NMAX,IAST,IU1,IU2,TOLL)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO 
          DIMENSION X(1),Y(1),U(1),IPX(1),IPY(1),LPX(1),LPY(1),
     .    IWD(1),IAST(1),DATO(120) 
          REWIND IU1
          REWIND IU2
          NM=0 
          DO 20 I=1,NN 
          READ(IU1,1000) (DATO(K),K=1,120) 
          IF(LPX(I).EQ.0) GOTO 15  
          XC=X(I)  
          YC=Y(I)  
          CALL AWINDW(XC,YC,OLATO,OLATO,X,Y,IPX,IPY,LPX,LPY,IWD,NMAX,  
     .    NWDMAX,N,NWD)
          SUM=0.D0 
          DO 10 J=1,NWD
          LL=IWD(J)
          SUM=SUM+U(LL)
10        CONTINUE 
          SUM=SUM-U(I) 
          EMED=SUM/NWD 
          ANOISE=U(I)-EMED 
          ZNOISE=ANOISE/SNOISE 
          IF(DABS(ZNOISE).GT.TOLL) THEN
          NM=NM+1  
          DATO(IAST(IND))='*'  
          ENDIF
15        WRITE(IU2,1000) (DATO(K),K=1,120)
20        CONTINUE 
          RETURN
1000      FORMAT(120A1)
          END  
           SUBROUTINE TRIMT(Z,U,IPZ,LPZ,SNOISE,OLATO,NMAX,IND,N,NN,NM,  
     .    IAST,IU1,IU2,TOLL)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO 
          DIMENSION Z(1),U(1),IPZ(1),LPZ(1),IAST(1),DATO(120)  
          OLATO2=OLATO/2.D0
          REWIND IU1
          REWIND IU2
          NM=0 
          DO 20 I=1,NN 
          READ(IU1,1000) (DATO(K),K=1,120) 
          IF(LPZ(I).EQ.0) GOTO 15  
          ZC=Z(I)  
          ZMI=ZC-OLATO2
          ZMA=ZC+OLATO2
          CALL ABINRY(ZMI,Z,IPZ,NMAX,LNI,LNS,LPI,LPS,N)
          LZI=LPS  
          CALL ABINRY(ZMA,Z,IPZ,NMAX,LNI,LNS,LPI,LPS,N)
          LZS=LPI  
          SUM=0.D0 
          DO 10 J=LZI,LZS  
          LL=IPZ(J)
          SUM=SUM+U(LL)
10        CONTINUE 
          SUM=SUM-U(I) 
          EMED=SUM/(LZS-LZI)
          NOISE=U(I)-EMED  
          ZNOISE=NOISE/SNOISE  
          IF(DABS(ZNOISE).GT.TOLL) THEN
          NM=NM+1  
          DATO(IAST(IND))='*'  
          ENDIF
15        WRITE(IU2,1000) (DATO(K),K=1,120)
20        CONTINUE 
          RETURN
1000      FORMAT(120A1)
          END  
          SUBROUTINE ZERO(COV,XCOV,IX0,X0,N0,NCOV) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION COV(1),XCOV(1),X0(1),IX0(1)
          N0=0 
          DO 1 I=1,3
          IX0(I)=0 
          X0(I)=0.D0
1         CONTINUE 
          DO 5 I=NCOV,1,-1 
          IF(COV(I).NE.0.D0) GOTO 8
5         CONTINUE 
8         NCOV=I
          DO 10 I=1,NCOV-1 
          IF(COV(I+1).EQ.0.D0) GOTO 10 
          IF(COV(I).EQ.0.D0) THEN  
          N0=N0+1  
          IX0(N0)=I
          X0(N0)=XCOV(I)
          ELSE 
          ISGN1=COV(I)/DABS(COV(I))
          ISGN2=COV(I+1)/DABS(COV(I+1))
          IF(ISGN1.NE.ISGN2) THEN  
          N0=N0+1  
          IX0(N0)=I
          X0(N0)=(XCOV(I+1)*DABS(COV(I))+XCOV(I)*DABS(COV(I+1)))/  
     .(DABS(COV(I))+DABS(COV(I+1)))
          ENDIF
          ENDIF
          IF(N0.GE.3) RETURN
10        CONTINUE 
          RETURN
          END  
C******************************     ADIACE
      SUBROUTINE ADIACE (NUMERO,INDIR,JCRONA,JGRAFO,LIVELW,LIV,
     .IPUNT,INFORM,IP,IC,IN,N,IS)
      DIMENSION NUMERO(1),INDIR(1),JCRONA(1),JGRAFO(1),
     .LIVELW(1),IPUNT(1),INFORM(1)
    1 IP=IP+1
      L=INDIR(IP)
      L1=JCRONA(L)
      L2=JCRONA(L+1)-1
      J=0
      DO 10 I=L1,L2
      L=JGRAFO(I)
      IF(LIVELW(L).NE.LIV.OR.NUMERO(L).NE.0) GO TO 10
      J=J+1
      IPUNT(J)=L
      INFORM(L)=JCRONA(L+1)-JCRONA(L)
   10 CONTINUE
      NL=J
      IF(NL.NE.0) THEN
      CALL RIORDI(INFORM,IPUNT,NL,N,1)
      DO 20 I=1,NL
      IN=IN+1
      L=IPUNT(I)
      NUMERO(L)=IN
      INDIR(IN)=L
   20 CONTINUE
      ENDIF
      IF(IP.LT.IS*IC+(1-IS)*IN) GO TO 1
      RETURN
      END
C******************************     BANDA
      SUBROUTINE BANDA(INDICI,KRIGA,IBANDA,JBANDA,
     .MMAX,M,N,IT,JT,LB)
      DIMENSION INDICI(MMAX,1),KRIGA(1),IBANDA(1),JBANDA(1)
      DO 1 I=1,N
      IBANDA(I)=0
      JBANDA(I)=0
    1 CONTINUE
      DO 5 I=1,M
      L=KRIGA(I)
      LL=1000000000
      DO 2 J=1,L
      LL=MIN0(LL,INDICI(I,J))
    2 CONTINUE
      IF(LL.LT.0) GO TO 5
      DO 3 J1=1,L
      K=INDICI(I,J1)
      DO 3 J2=1,L
      IBANDA(K)=MAX0(IBANDA(K),K-INDICI(I,J2))
      JBANDA(K)=MAX0(JBANDA(K),INDICI(I,J2)-K)
    3 CONTINUE
    5 CONTINUE
      IT=N
      JT=N
      LB=0
      DO 10 I=1,N
      IT=IT+IBANDA(I)
      JT=JT+JBANDA(I)
      LB=MAX0(LB,IBANDA(I))
   10 CONTINUE
      RETURN
      END
C******************************     BINARY
          SUBROUTINE BINARY(IVALRE,INFORM,IPUNT,NMAX,LNINF,LNSUP,
     .    LPINF,LPSUP,N)
          DIMENSION INFORM(1),IPUNT(1)
          LPINF=0
          LPSUP=0
          LNINF=0
          LNSUP=0
          NI=1
          NF=N
          IPUNT(N+1)=NMAX
          INFORM(NMAX)=100000000
          LI=IPUNT(NI)
          LF=IPUNT(NF)
          IF(IVALRE.LT.INFORM(LI)) THEN
          LPINF=0
          LPSUP=1
          LNINF=0
          LNSUP=IPUNT(1)
          IPUNT(N+1)=0
          INFORM(NMAX)=0
          RETURN
          ENDIF
          IF(IVALRE.GT.INFORM(LF)) THEN
          LPINF=N
          LPSUP=N+1
          LNINF=IPUNT(N)
          LNSUP=NMAX+1
          IPUNT(N+1)=0
          INFORM(NMAX)=0
          RETURN
          ENDIF
          NF=NF+1
5         NC=(NI+NF)/2
          LC=IPUNT(NC)
          IF(IVALRE.EQ.INFORM(LC)) THEN
          DO 10 I=NC,NF
          L=IPUNT(I)
          IF(IVALRE.NE.INFORM(L)) GO TO 12
          LPSUP=I
          LNSUP=IPUNT(I)
10        CONTINUE
12        DO 15 I=NC,NI,-1
          L=IPUNT(I)
          IF(IVALRE.NE.INFORM(L)) GO TO 18
          LPINF=I
          LNINF=IPUNT(I)
15        CONTINUE
18        IPUNT(N+1)=0
          INFORM(NMAX)=0
          RETURN
          ENDIF
          IF((NF-NI).LE.1) THEN
          LPINF=NI
          LPSUP=NF
          LNINF=IPUNT(NI)
          LNSUP=IPUNT(NF)
          IPUNT(N+1)=0
          INFORM(NMAX)=0
          RETURN
          ENDIF
          IF(IVALRE.LT.INFORM(LC)) THEN
          NF=NC
          ENDIF
          IF(IVALRE.GT.INFORM(LC)) THEN
          NI=NC
          ENDIF
          GOTO 5
          END
C******************************     CLASSI
      SUBROUTINE CLASSI(INFORM,KONT,IPUNT,M,N)
      DIMENSION INFORM(1),KONT(1),IPUNT(1)
      DO 1 I=1,N
      KONT(I+1)=0
    1 CONTINUE
      KONT(1)=0
      DO 5 I=1,M
      K=INFORM(I)
      IF(K.NE.0) THEN
      KONT(K+1)=KONT(K+1)+1
      ENDIF
    5 CONTINUE
      DO 10 I=1,N
      KONT(I+1)=KONT(I)+KONT(I+1)
   10 CONTINUE
      DO 15 I=1,M
      K1=INFORM(I)
      IF(K1.NE.0) THEN
      KONT(K1)=KONT(K1)+1
      K2=KONT(K1)
      IPUNT(K2)=I
      ENDIF
   15 CONTINUE
      DO 20 I=N,1,-1
      KONT(I+1)=KONT(I)
   20 CONTINUE
      KONT(1)=0
      RETURN
      END
C******************************     CONNEX
      SUBROUTINE CONNEX(JCRONA,JGRAFO,LIVEL,LIPRIM,LIULT,NOME,N)
      DIMENSION JCRONA(1),JGRAFO(1),LIVEL(1),LIPRIM(1),LIULT(1),NOME(1)
      DO 5 I=1,N
      LIVEL(I)=0
    5 CONTINUE
      IL=1
      LIV=1
      LIVEL(1)=1
      LIPRIM(1)=1
      NU=1
   10 LIV=LIV+1
      IN=0
      DO 20 I=1,NU
      K=LIPRIM(I)
      L1=JCRONA(K)
      L2=JCRONA(K+1)-1
      DO 20 J=L1,L2
      L=JGRAFO(J)
      IF(LIVEL(L).NE.0) GO TO 20
      IN=IN+1
      LIVEL(L)=LIV
      LIULT(IN)=L
   20 CONTINUE
      IF(IN.EQ.0) GO TO 40
      IL=IL+IN
      NU=IN
      DO 30 I=1,NU
      LIPRIM(I)=LIULT(I)
   30 CONTINUE
      GO TO 10
   40 IF(IL.LT.N) GO TO 50
      RETURN
   50 WRITE(6,100)
      DO 60 I=1,N,20
      I1=I
      I2=MIN0(N,I+19)
      WRITE(6,120) (J,J=I1,I2)
      WRITE(6,120) (NOME(J),J=I1,I2)
      WRITE(6,120) (LIVEL(J),J=I1,I2)
      WRITE(6,110)
   60 CONTINUE
      STOP
  100 FORMAT(///,5X,'GRAFO NON CONNESSO',///)
  110 FORMAT(/)
  120 FORMAT(5X,20I5)
      END
C******************************     CONTA
          SUBROUTINE CONTA(NPUNTO,NPUNT2,KPUNTO,NMS,NPS,NPS2,
     .NP,IFOTO,ITIPO,ISCRI)
          CHARACTER*1 DATO,PTO
          DIMENSION NPUNTO(1),NPUNT2(1),KPUNTO(1),DATO(132),KSERV(20)
          DATA PTO /'.'/
          KPUNTO(1)=0
          IF(ITIPO.EQ.2) GO TO 2
          REWIND 2
          DO 1  I=1,NPS
          READ(2,219) (DATO(J),J=1,132)
          IP=NPUNTO(I)
          IF(IP.EQ.0) GO TO 1
          KPUNTO(IP+1)=0
          IF(DATO(21).EQ.PTO) KPUNTO(IP+1)=KPUNTO(IP+1)+1
          IF(DATO(46).EQ.PTO) KPUNTO(IP+1)=KPUNTO(IP+1)+1
          IF(DATO(70).EQ.PTO) KPUNTO(IP+1)=KPUNTO(IP+1)+1
1         CONTINUE
2         IF(ITIPO.EQ.1) GO TO 20
          DO 5 I=1,NPS2
          IP=NPUNT2(I)
          IF(IP.EQ.0) GO TO 5
          KPUNTO(IP+1)=3
5         CONTINUE
          IF(IFOTO.EQ.1) GO TO 12
          DO 10 I=1,NMS
          IP=NPUNT2(NPS2+I)
          IF(IP.EQ.0) GO TO 10
          KPUNTO(IP+1)=7
10        CONTINUE
          GO TO 20
12        DO 15 I=1,NMS
          IP=NPUNT2(NPS2+I)
          IF(IP.EQ.0) GO TO 15
          KPUNTO(IP+1)=6
15        CONTINUE
20        IF(ITIPO.EQ.2) GO TO 30
          REWIND 2
          DO 25 I=1,NPS
          READ(2,219) (DATO(J),J=1,132)
          IP=NPUNTO(I)
          IF(IP.EQ.0) GO TO 25
          IF(DATO(90).EQ.PTO) KPUNTO(IP+1)=KPUNTO(IP+1)+1
25        CONTINUE
30        DO 40 I=1,NP
          KPUNTO(I+1)=KPUNTO(I+1)+KPUNTO(I)
40        CONTINUE
          IF(ISCRI.EQ.0) RETURN
          WRITE(6,105)
          DO 50 I=1,NP,20
          I2=MIN0(I+19,NP)
          I3=I2-I+1
          DO 45 J=1,I3
          KSERV(J)=KPUNTO(I+J)-KPUNTO(I+J-1)
45        CONTINUE
          WRITE(6,100) (J,J=I,I2)
          WRITE(6,100) (KSERV(J),J=1,I3)
          WRITE(6,103) (KPUNTO(J+1),J=I,I2)
50        CONTINUE
          RETURN
100       FORMAT(5X,20I5)
103       FORMAT(5X,20I5//)
105       FORMAT(///5X,'CONTATORE DELLE INCOGNITE',//5X,
     .'(NUMERO ORDINE DEL VERTICE, N. DI INCOGNITE DEL VERTICE',
     .', INTEGRATORE DEL N. DI INCOGNITE)'/)
219       FORMAT(132A1)
          END
C******************************     DIAME
      SUBROUTINE DIAME(JCRONA,JGRAFO,LIVELV,LIVELU,LIPRIM,LIULT,LIVLT,
     .ISRV,LIMAX,KAVMAX,KAUMAX,INI,IFI,NN,N)
      DIMENSION JCRONA(1),JGRAFO(1),LIVELV(1),LIVELU(1),LIPRIM(1),LIULT(
     .1),LIVLT(1),ISRV(1)
      L=1000000000
      INI=0
      DO 10 I=1,N
      KRONA=JCRONA(I+1)-JCRONA(I)
      IF(KRONA.GE.L.OR.KRONA.EQ.1) GO TO 10
      L=KRONA
      INI=I
   10 CONTINUE
      IF(INI.EQ.0) RETURN
      CALL LIVEL(INI,JCRONA,JGRAFO,LIVELV,LIPRIM,LIULT,0,
     .LIVMAX,KAVMAX,NV,NN,N)
      IF(NV.NE.1) THEN
      DO 15 I=1,NV
      J=LIULT(I)
      ISRV(J)=JCRONA(J+1)-JCRONA(J)
   15 CONTINUE
      CALL RIORDI(ISRV,LIULT,NV,N,1)
      ENDIF
   20 DO 25 I=1,NV
      LIVLT(I)=LIULT(I)
   25 CONTINUE
      L=N
      DO 30 I=1,NV
      K=LIVLT(I)
      CALL LIVEL(K,JCRONA,JGRAFO,LIVELU,LIPRIM,LIULT,0,
     .LIUMAX,KAUMAX,NU,NN,N)
      IF(NU.NE.1) THEN
      DO 28 IS=1,NU
      J=LIULT(IS)
      ISRV(J)=JCRONA(J+1)-JCRONA(J)
   28 CONTINUE
      CALL RIORDI(ISRV,LIULT,NU,N,1)
      ENDIF
      IF(LIUMAX.GT.LIVMAX) GO TO 40
      IF(L.LE.KAUMAX) GO TO 30
      L=KAUMAX
      IFI=K
   30 CONTINUE
      CALL LIVEL(IFI,JCRONA,JGRAFO,LIVELU,LIPRIM,LIULT,0,
     .LIUMAX,KAUMAX,NU,NN,N)
      LIMAX=LIUMAX
      IF(NU.NE.1) THEN
      DO 35 I=1,NU
      J=LIULT(I)
      ISRV(J)=JCRONA(J+1)-JCRONA(J)
   35 CONTINUE
      CALL RIORDI(ISRV,LIULT,NU,N,1)
      ENDIF
      RETURN
   40 INI=K
      LIVMAX=LIUMAX
      KAVMAX=KAUMAX
      NV=NU
      DO 50 I=1,N
      LIVELV(I)=LIVELU(I)
   50 CONTINUE
      GO TO 20
      END
C******************************     GIBBS
      SUBROUTINE GIBBS(NUMERO,INDIR,JCRONA,JGRAFO,IA,IB,IC,ID,IE,IF,IG,
     .N)
      DIMENSION NUMERO(1),INDIR(1),JCRONA(1),JGRAFO(1),IA(1),IB(1),
     .IC(1),ID(1),IE(1),IF(1),IG(1)
      CALL DIAME(JCRONA,JGRAFO,IA,IB,IC,ID,IE,IF,LIMAX,KAVMAX,
     .KAUMAX,INI,IFI,NN,N)
      IF(INI.NE.0) THEN
      CALL RIDUCE(JCRONA,JGRAFO,IA,IB,IC,ID,IE,IF,IG,NUMERO,INDIR,
     .LIMAX,KAVMAX,KAUMAX,ISC,NN,N)
      CALL MEGLIO(NUMERO,INDIR,JCRONA,JGRAFO,IC,IA,IB,ID,IE,INI,IFI,
     .LIMAX,ISC,NN,N)
      ELSE
      DO 10 I=1,N
      NUMERO(I)=I
      INDIR(I)=I
   10 CONTINUE
      ENDIF
      RETURN
      END
C******************************     LISTA
          SUBROUTINE LISTA(NPUNTO,NUMERO,KO,LS,N)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION NPUNTO(1),NUMERO(1)
          DO 1 I=1,N
          K=NPUNTO(I)
          IF(K.EQ.0) GO TO 1
          IND=(K-KO)/LS
          NPUNTO(I)=NUMERO(IND)
1         CONTINUE
          RETURN
          END
C******************************     LIVEL
      SUBROUTINE LIVEL(INI,JCRONA,JGRAFO,LIVELW,LIPRIM,LIULT,LI0,
     .LIV,KAMAX,NU,NN,N)
      DIMENSION JCRONA(1),JGRAFO(1),LIVELW(1),LIPRIM(1),LIULT(1)
      IF(LI0.EQ.0)THEN
      DO 1 I=1,N
      LIVELW(I)=0
    1 CONTINUE
      LIV=1
      ELSE
      LIV=LI0+1
      ENDIF
      LIVELW(INI)=LIV
      LIPRIM(1)=INI
      NU=1
      NN=1
      KAMAX=1
   10 LIV=LIV+1
      IN=0
      DO 20 I=1,NU
      K=LIPRIM(I)
      L1=JCRONA(K)
      L2=JCRONA(K+1)-1
      DO 20 J=L1,L2
      L=JGRAFO(J)
      IF(LIVELW(L).EQ.0) THEN
      IN=IN+1
      IF(KAMAX.LT.NU) KAMAX=NU
      LIVELW(L)=LIV
      LIULT(IN)=L
      ENDIF
   20 CONTINUE
      IF(IN.EQ.0) THEN
      LIV=LIV-1
      RETURN
      ENDIF
      NU=IN
      NN=NN+IN
      DO 30 I=1,NU
      LIPRIM(I)=LIULT(I)
   30 CONTINUE
      GO TO 10
      END
C******************************     MEGLIO
      SUBROUTINE MEGLIO(NUMERO,INDIR,JCRONA,JGRAFO,LIVELW,KLIVEL,JLIVEL,
     .IPUNT,INFORM,INI,IFI,LIMAX,ISC,NN,N)
      DIMENSION NUMERO(1),INDIR(1),JCRONA(1),JGRAFO(1),LIVELW(1),
     .KLIVEL(1),JLIVEL(1),IPUNT(1),INFORM(1)
      JSC=0
      KRONA=JCRONA(INI+1)-JCRONA(INI)
      LRONA=JCRONA(IFI+1)-JCRONA(IFI)
      IF(KRONA.GT.LRONA) THEN
      DO 10 I=1,N
      LIV=LIVELW(I)
      IF(LIV.LE.LIMAX) LIVELW(I)=LIMAX+1-LIV
   10 CONTINUE
      INI=IFI
      JSC=1
      ENDIF
      CALL CLASSI (LIVELW,KLIVEL,JLIVEL,N,LIMAX)
      DO 20 I=1,N
      NUMERO(I)=0
      INDIR(I)=0
   20 CONTINUE
      LIV=1
      NUMERO(INI)=1
      INDIR(1)=INI
      IN=1
      IL=0
      IK=0
   30 CALL ADIACE (NUMERO,INDIR,JCRONA,JGRAFO,LIVELW,LIV,IPUNT,INFORM,
     .IL,IN,IN,N,0)
      K1=KLIVEL(LIV)+1
      K2=KLIVEL(LIV+1)
      J=0
      DO 40 I=K1,K2
      L=JLIVEL(I)
      IF(NUMERO(L).NE.0) GO TO 40
      J=J+1
      IPUNT(J)=L
      INFORM(L)=JCRONA(L+1)-JCRONA(L)
   40 CONTINUE
      NL=J
      IF(NL.NE.0) THEN
      CALL RIORDI(INFORM,IPUNT,NL,N,1)
      IN=IN+1
      L=IPUNT(1)
      NUMERO(L)=IN
      INDIR(IN)=L
      GO TO 30
      ELSE
      IF(LIV.EQ.LIMAX) GO TO 50
      LIV=LIV+1
      CALL ADIACE (NUMERO,INDIR,JCRONA,JGRAFO,LIVELW,LIV,IPUNT,INFORM,
     .IK,IL,IN,N,1)
      GO TO 30
      ENDIF
   50 IN=NN
      DO 60 I=1,N
      IF(LIVELW(I).LE.LIMAX) GO TO 60
      IN=IN+1
      NUMERO(I)=IN
      INDIR(IN)=I
   60 CONTINUE
      IF(ISC.NE.JSC) RETURN
      NN2=NN/2
      DO 70 I=1,NN2
      ISRV=INDIR(I)
      INDIR(I)=INDIR(NN+1-I)
      INDIR(NN+1-I)=ISRV
   70 CONTINUE
      DO 80 I=1,N
      IF(LIVELW(I).GT.LIMAX) GO TO 80
      L=NUMERO(I)
      NUMERO(I)=NN+1-L
   80 CONTINUE
      RETURN
      END
C******************************     NORMAS
      SUBROUTINE NORMAS(IGRAFO,ICRONA,JGRAFO,JCRONA,INDICI,KRIGA,
     .IBANDA,JBANDA,KGRAFO,MMAX,N)
      DIMENSION IGRAFO(1),ICRONA(1),JGRAFO(1),JCRONA(1),INDICI(MMAX,1),
     .KRIGA(1),IBANDA(1),JBANDA(1),KGRAFO(1)
      DO 1 I=1,N
      KGRAFO(I)=0
    1 CONTINUE
      KCOLL=0
      JCRONA(1)=1
      DO 40 I=1,N
      I1=ICRONA(I)+1
      I2=ICRONA(I+1)
      IF(I1.LE.I2) THEN
      DO 10 J=I1,I2
      MISURA=IGRAFO(J)
      LL=KRIGA(MISURA)
      DO 10 L=1,LL
      NCOLL=INDICI(MISURA,L)
      KGRAFO(NCOLL)=NCOLL
   10 CONTINUE
      JMIN=I-IBANDA(I)
      JMAX=I+JBANDA(I)
      DO 20 J=JMIN,JMAX
      IF(KGRAFO(J).EQ.0) GOTO 20
      KCOLL=KCOLL+1
      JGRAFO(KCOLL)=J
   20 CONTINUE
      JCRONA(I+1)=KCOLL+1
      DO 30 J=JMIN,JMAX
      KGRAFO(J)=0
   30 CONTINUE
      ELSE
      KCOLL=KCOLL+1
      JGRAFO(KCOLL)=I
      JCRONA(I+1)=KCOLL+1
      ENDIF
   40 CONTINUE
      RETURN
      END
C******************************     NUMERA
      SUBROUTINE NUMERA(INDICI,KRIGA,KO,LS,NOME,NUMERO,
     .ISRV,MMAX,NPTMAX,M,N)
      DIMENSION INDICI(MMAX,1),KRIGA(1),NOME(1),NUMERO(1),
     .ISRV(1)
      DO 1 I=1,NPTMAX
      NOME(I)=0
      NUMERO(I)=0
    1 CONTINUE
      NPT=0
      DO 10 I=1,M
      L=KRIGA(I)
      DO 5 J=1,L
      K=INDICI(I,J)
      IF(K.LT.0) GO TO 5
      IND=(K-KO)/LS
      NOME(IND)=K
      NPT=MAX0(NPT,IND)
    5 CONTINUE
   10 CONTINUE
      J=0
      DO 20 I=1,NPT
      IF(NOME(I).EQ.0) GO TO 20
      J=J+1
      NUMERO(I)=J
      NOME(J)=NOME(I)
   20 CONTINUE
      N=J
      GO TO 30
      ENTRY RINUME(INDICI,KRIGA,KO,LS,NOME,NUMERO,
     .ISRV,MMAX,M,N)
      IF(N.EQ.0) GO TO 30
      DO 22 I=1,N
      L=NUMERO(I)
      ISRV(L)=NOME(I)
   22 CONTINUE
      DO 25 I=1,N
      NOME(I)=ISRV(I)
   25 CONTINUE
   30 DO 50 I=1,M
      L=KRIGA(I)
      DO 40 J=1,L
      K=INDICI(I,J)
      IF(K.LT.0) GO TO 40
      IND=(K-KO)/LS
      INDICI(I,J)=NUMERO(IND)
   40 CONTINUE
   50 CONTINUE
      RETURN
      END
C******************************     ORDINI
          SUBROUTINE ORDINI(INFORM,IPUNT,N,ISGN)
          DIMENSION INFORM(1),IPUNT(1)
          NS=N
          DO 1 I=1,N
          IPUNT(I)=I
1         CONTINUE
          GO TO 2
          ENTRY RIORDI(INFORM,IPUNT,N,NN,ISGN)
          NS=NN
2         IF(N.EQ.1) RETURN
          IF(ISGN.NE.1) THEN
          DO 5 I=1,N
          J=IPUNT(I)
          INFORM(J)=ISGN*INFORM(J)
5         CONTINUE
          ENDIF
          NM1=N-1
          K=0
          DO 8 I=1,NM1
          J1=IPUNT(I)
          J2=IPUNT(I+1)
          IF(INFORM(J1).GT.INFORM(J2)) K=K+1
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
          INFORM(NS+1)=-10000000
          DO 18 I=2,N
          J=I
17        L1=IPUNT(J)
          L2=IPUNT(J+1)
          IF(INFORM(L1).LE.INFORM(L2)) GO TO 18
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
          INFORM(J)=ISGN*INFORM(J)
25        CONTINUE
          RETURN
          END
          SUBROUTINE PUNTA(NPUNTO,LPUNTO,KO,LS,N,NN)
          DIMENSION LPUNTO(1),NPUNTO(1)
          DO 5 I=1,NN
          LPUNTO(I)=0
5         CONTINUE
          DO 10 I=1,N
          J=NPUNTO(I)
          IF(J.EQ.0) GO TO 10
          IND=(J-KO)/LS
          LPUNTO(IND)=I
10        CONTINUE
          RETURN
          END
C******************************     RIDUCE
      SUBROUTINE RIDUCE(JCRONA,JGRAFO,LIVELV,LIVELU,LIVELW,KOMPW,KARD,
     .KOMP,KONTV,KONTU,KONTW,LIMAX,KAVMAX,KAUMAX,ISC,NN,N)
      DIMENSION JCRONA(1),JGRAFO(1),LIVELV(1),LIVELU(1),LIVELW(1),
     .KOMPW(1),KARD(1),KOMP(1),KONTV(1),KONTU(1),KONTW(1)
      DO 10 I=1,N
      LIVELU(I)=LIMAX+1-LIVELU(I)
      LIVELW(I)=0
   10 CONTINUE
      IF(NN.LT.N) THEN
      DO 15 I=1,N
      IF(LIVELV(I).NE.0) GO TO 15
      LIVELV(I)=LIMAX+1
      LIVELU(I)=LIMAX+1
   15 CONTINUE
      ENDIF
      DO 20 I=1,LIMAX+1
      KONTW(I)=0
   20 CONTINUE
      IR=0
      DO 40 I=1,N
      IF(LIVELV(I).EQ.LIVELU(I)) THEN
      IR=IR+1
      LIV=LIVELV(I)
      LIVELW(I)=LIV
      KONTW(LIV)=KONTW(LIV)+1
      ENDIF
   40 CONTINUE
      IF(IR.EQ.N) THEN
      ISC=0
      ELSE
      ICO=0
      DO 60 I=1,N
      IF(LIVELW(I).EQ.0) THEN
      ICO=ICO+1
      LI0=100*ICO
      CALL LIVEL(I,JCRONA,JGRAFO,LIVELW,KONTV,KONTU,LI0,LIW,KAMAX,
     .NU,IR,N)
      KARD(ICO)=IR
      ENDIF
   60 CONTINUE
      CALL ORDINI(KARD,KONTV,ICO,-1)
      CALL PUNTA(KONTV,KONTU,0,1,ICO,ICO)
      DO 70 I=1,N
      KO=LIVELW(I)/100
      IF(KO.EQ.0) THEN
      KOMPW(I)=0
      ELSE
      KOMPW(I)=KONTU(KO)
      ENDIF
   70 CONTINUE
      CALL CLASSI(KOMPW,KARD,KOMP,N,ICO)
      DO 140 I=1,ICO
      DO 80 L=1,LIMAX
      KONTV(L)=0
      KONTU(L)=0
   80 CONTINUE
      K1=KARD(I)+1
      K2=KARD(I+1)
      DO 90 K=K1,K2
      J=KOMP(K)
      LIV=LIVELV(J)
      LIU=LIVELU(J)
      KONTV(LIV)=KONTV(LIV)+1
      KONTU(LIU)=KONTU(LIU)+1
   90 CONTINUE
      KOVMAX=KONTW(1)+KONTV(1)
      KOUMAX=KONTW(1)+KONTU(1)
      DO 100 L=2,LIMAX
      IF(KONTV(L).GT.0) KOVMAX=MAX0(KOVMAX,KONTW(L)+KONTV(L))
      IF(KONTU(L).GT.0) KOUMAX=MAX0(KOUMAX,KONTW(L)+KONTU(L))
  100 CONTINUE
      IF(KOVMAX.LE.KOUMAX.OR.(KOVMAX.EQ.KOUMAX.AND.KAVMAX.LE.KAUMAX))
     .THEN
      DO 110 K=K1,K2
      J=KOMP(K)
      LIV=LIVELV(J)
      KONTW(LIV)=KONTW(LIV)+1
      LIVELW(J)=LIV
  110 CONTINUE
      IF(I.EQ.1) ISC=0
      ELSE
      DO 130 K=K1,K2
      J=KOMP(K)
      LIU=LIVELU(J)
      KONTW(LIU)=KONTW(LIU)+1
      LIVELW(J)=LIU
  130 CONTINUE
      IF(I.EQ.1) ISC=1
      ENDIF
  140 CONTINUE
      ENDIF
      RETURN
      END
C******************************     REVERS
      SUBROUTINE REVERS(INDICI,KRIGA,NOME,NUMERO,IBANDA,JBANDA,
     .MMAX,M,N,IT,JT)
      DIMENSION INDICI(MMAX,1),KRIGA(1),NOME(1),NUMERO(1),
     .IBANDA(1),JBANDA(1)
      IF(IT.LE.JT) RETURN
      N2=N/2
      DO 5 I=1,N2
      ISRV=NOME(I)
      NOME(I)=NOME(N+1-I)
      NOME(N+1-I)=ISRV
    5 CONTINUE
      DO 10 I=1,N
      L=NUMERO(I)
      NUMERO(I)=N+1-L
   10 CONTINUE
      DO 15 I=1,M
      L=KRIGA(I)
      DO 12 J=1,L
      K=INDICI(I,J)
      IF(K.LT.0) GO TO 12
      INDICI(I,J)=N+1-K
   12 CONTINUE
   15 CONTINUE
      DO 20 I=1,N
      ISRV=IBANDA(I)
      JSRV=JBANDA(N-I+1)
      IBANDA(I)=JSRV
      JBANDA(N-I+1)=ISRV
   20 CONTINUE
      KSRV=IT
      IT=JT
      JT=KSRV
      RETURN
      END
C******************************     TRASPS
      SUBROUTINE TRASPS(INDICI,KRIGA,ICRONA,IGRAFO,MMAX,M,N)
      DIMENSION INDICI(MMAX,1),KRIGA(1),ICRONA(1),IGRAFO(1)
      DO 1 I=1,N
      ICRONA(I+1)=0
    1 CONTINUE
      ICRONA(1)=0
      DO 5 I=1,M
      L=KRIGA(I)
      LL=1000000000
      DO 2 J=1,L
      LL=MIN0(LL,INDICI(I,J))
    2 CONTINUE
      IF(LL.LT.0) GO TO 5
      DO 3 J=1,L
      K=INDICI(I,J)
      ICRONA(K+1)=ICRONA(K+1)+1
    3 CONTINUE
    5 CONTINUE
      DO 10 I=1,N
      ICRONA(I+1)=ICRONA(I)+ICRONA(I+1)
   10 CONTINUE
      DO 15 I=1,M
      L=KRIGA(I)
      LL=1000000000
      DO 12 J=1,L
      LL=MIN0(LL,INDICI(I,J))
   12 CONTINUE
      IF(LL.LT.0) GO TO 15
      DO 13 J=1,L
      K1=INDICI(I,J)
      ICRONA(K1)=ICRONA(K1)+1
      K2=ICRONA(K1)
      IGRAFO(K2)=I
   13 CONTINUE
   15 CONTINUE
      DO 20 I=N,1,-1
      ICRONA(I+1)=ICRONA(I)
   20 CONTINUE
      ICRONA(1)=0
      RETURN
      END
C******************************     WGRAFO
      SUBROUTINE WGRAFO(NOME,GRAFO,CRONA,BANDI,BANDA,N,INIZIO)
      INTEGER GRAFO,CRONA,BANDI,BANDA
      DIMENSION NOME(1),GRAFO(1),CRONA(1),BANDI(1),BANDA(1),IFORM(23),
     .ISRV(18)
c      DATA IFORM/'(5X,','2I5,','4X, ',18*'    ','4X, ','3I5)'/
c      INTERO='I5, '
c      ICS='5X, '
c      IF(INIZIO.EQ.0) WRITE(6,100)
c      IF(INIZIO.EQ.1) WRITE(6,110)
c      DO 5 L=1,18
c      ISRV(L)=0
c    5 CONTINUE
c      DO 30 I=1,N
c      DO 10 J=1,18
c      IFORM(J+3)=ICS
c   10 CONTINUE
c      KRONA=CRONA(I+1)-CRONA(I)
c      INI=CRONA(I)+1-INIZIO
c      IFI=CRONA(I+1)-INIZIO
c      K=MIN0(KRONA,18)
c      DO 20 J=1,K
c      IFORM(J+3)=INTERO
c   20 CONTINUE
c      IFJ=MIN0(IFI,INI+17)
c      LL=IFJ-INI+1
c      DO 21 L=1,LL
c      ISRV(L)=GRAFO(INI+L-1)
c   21 CONTINUE
c      WRITE(6,IFORM) I,NOME(I),(ISRV(J),J=1,LL),KRONA,BANDI(I),BANDA(I)
c      J=0
c   25 J=J+1
c      DO 23 L=1,18
c      ISRV(L)=0
c   23 CONTINUE
c      INJ=INI+18*J
c      IF(INJ.GT.IFI) GO TO 30
c      IFJ=MIN0(IFI,INJ+17)
c      LL=IFJ-INJ+1
c      DO 28 L=1,LL
c      ISRV(L)=GRAFO(INJ+L-1)
c   28 CONTINUE
c      WRITE(6,120) (ISRV(K),K=1,LL)
c      GO TO 25
c   30 CONTINUE
      RETURN
c  100 FORMAT(/////5X,'MATRICE DISEGNO SIMBOLICA',
c     .//,5X,'(NUMERO DI ORDINE E NOME DEL PUNTO - MISURE CORONA,',
c     .' CON TOTALE - BANDA INDIETRO E AVANTI)',/)
c  110 FORMAT(/////5X,'MATRICE NORMALE SIMBOLICA',
c     .//,5X,'(NUMERO DI ORDINE E NOME DEL PUNTO - PUNTI CORONA,',
c     .' CON TOTALE - BANDA INDIETRO E AVANTI)',/)
c  120 FORMAT(19X,18I5)
      END
C*************************          ABINRY          *************************
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
          LPINF=0
          LPSUP=1
          LNINF=0
          LNSUP=IPUNT(1)
          IPUNT(N+1)=0
          HNFORM(NMAX)=0.D0
          RETURN
          ENDIF
          IF(VALORE.GT.HNFORM(LF)) THEN
          LPINF=N
          LPSUP=N+1
          LNINF=IPUNT(N)
          LNSUP=NMAX+1
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
C*************************          AWINDW          *************************
          SUBROUTINE AWINDW(XC,YC,AMPX,AMPY,X,Y,IPX,IPY,LPX,LPY,
     .IWD,NMAX,NWDMAX,N,NWD)
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
          IF(NWD.GT.NWDMAX) WRITE(*,10) XC,YC,AMPX,AMPY,NWD
          RETURN
10        FORMAT(////5X,'NUMERO DEI PUNTI NELLA FINESTRA',/
     .    5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'//
     .    5X,'COORDINATE CENTRO',18X,'AMPIEZZA FINESTRA',8X,
     .    'NUMEROSITA'//5X,2F15.3,5X,2F10.3,5X,I10)
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
C*************************          BORDO          *************************
          SUBROUTINE BORDO(XC,YC,AMPX,AMPY,PASX,PASY,
     .X,Y,IPX,IPY,LPX,LPY,IWD,NMAX,NWDMAX,N,NWD)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),IPX(1),IPY(1),LPX(1),LPY(1),IWD(1)
          CALL AWINDW(XC,YC,AMPX+PASX/2.D0,AMPY+PASY/2.D0,
     .X,Y,IPX,IPY,LPX,LPY,IWD,NMAX,NWDMAX,N,NWD)
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
C*************************          FMAT
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
   30 WRITE(6,100) K
      STOP
  100 FORMAT(///,5X,'MATRICE NON RIDUCIBILE IN BANDA - EQUAZIONE',I5,///
     .)
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
   50 WRITE(*,100) I
      STOP
  100 FORMAT(///5X,'MATRICE NON DEFINITA POSITIVA - EL. DIAGONALE',I5/)
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
   50 WRITE(6,100) I
      STOP
  100 FORMAT(///,5X,'MATRICE NON DEFINITA POSITIVA - EL. DIAGONALE',I5/)
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
C*************************          TRASPC          *************************
      SUBROUTINE TRASPC(AR,IAR,H,AC,IAC,IC,IB,JB,MMAX,M,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER H
      DIMENSION AR(MMAX,1),IAR(MMAX,1),H(1),AC(1),IAC(1),IC(1),IB(1),JB
     .(1)
      IC(1)=0
      DO 5 I=1,N
      IC(I+1)=0
      IB(I)=0
      JB(I)=0
    5 CONTINUE
      DO 10 I=1,M
      L=H(I)
      DO 10 J1=1,L
      K=IAR(I,J1)
      IC(K+1)=IC(K+1)+1
      DO 10 J2=1,L
      IB(K)=MAX0(IB(K),K-IAR(I,J2))
      JB(K)=MAX0(JB(K),IAR(I,J2)-K)
   10 CONTINUE
      DO 20 I=1,N
      IC(I+1)=IC(I)+IC(I+1)
   20 CONTINUE
      DO 30 I=1,M
      L=H(I)
      DO 30 J=1,L
      K1=IAR(I,J)
      IC(K1)=IC(K1)+1
      K2=IC(K1)
      IAC(K2)=I
      AC(K2)=AR(I,J)
   30 CONTINUE
      DO 40 I=N,1,-1
      IC(I+1)=IC(I)
   40 CONTINUE
      IC(1)=0
      RETURN
      END
          FUNCTION DFN01(Z)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          REAL*4 SRW,SRV,ERF1
          SRW=DABS(Z)/DSQRT(2.D0)
          SRV=DBLE(ERF1(SRW))
          DFN01=.5D0+.5D0*SRV*Z/DABS(Z)
          RETURN
          END
C*************************    ERF1
      FUNCTION ERF1(X)
      IF(X.LT.0.)THEN
        ERF1=-GAMMP(.5,X**2)
      ELSE
        ERF1=GAMMP(.5,X**2)
      ENDIF
      RETURN
      END
C*************************    GAMMP
      FUNCTION GAMMP(A,X)
      IF(X.LT.0..OR.A.LE.0.)PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END

C*************************    GSER
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END

C*************************    GCF
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END

C*************************    GAMMLN
      FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END


          FUNCTION DFJ0(X)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          REAL*4 SRW,SRV,BESSJ0
          SRW=X
          SRV=BESSJ0(SRW)
          DFJ0=DBLE(SRV)
          RETURN
          END
C*************************    BESSJ0
      FUNCTION BESSJ0(X)
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     *    S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     *    -.2073370639D-5,.2093887211D-6/, Q1,Q2,Q3,Q4,Q5/-.1562499995D-
     *1,
     *    .1430488765D-3,-.6911147651D-5,.7621095161D-6,-.934945152D-7/
      DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,651619640.7D
     *0,
     *    -11214424.18D0,77392.33017D0,-184.9052456D0/,
     *    S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,
     *    9494680.718D0,59272.64853D0,267.8532712D0,1.D0/
      IF(ABS(X).LT.8.)THEN
        Y=X**2
        BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=ABS(X)
        Z=8./AX
        Y=Z**2
        XX=AX-.785398164
        BESSJ0=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
      ENDIF
      RETURN
      END

          FUNCTION DFJ1(X)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          REAL*4 SRV,SRW,BESSJ1
          SRV=X
          SRW=BESSJ1(SRV)
          DFJ1=2.D0*DBLE(SRW)
          RETURN
          END
C*************************    BESSJ1
      FUNCTION BESSJ1(X)
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     *    S1,S2,S3,S4,S5,S6
      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,242396853.1D0
     *,
     *    -2972611.439D0,15704.48260D0,-30.16036606D0/,
     *    S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,
     *    18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5
     *,
     *    -.240337019D-6/, Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3
     *,
     *    .8449199096D-5,-.88228987D-6,.105787412D-6/
      IF(ABS(X).LT.8.)THEN
        Y=X**2
        BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=ABS(X)
        Z=8./AX
        Y=Z**2
        XX=AX-2.356194491
        BESSJ1=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
     *      *SIGN(1.,X)
      ENDIF
      RETURN
      END
C************************** MEDVAR
          SUBROUTINE MEDVAR(X,EX,S2X,SX,N)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1)
          EX=0.D0
          S2X=0.D0
          DO 1 I=1,N
          EX=EX+X(I)
          S2X=S2X+X(I)**2
1         CONTINUE
          EX=EX/N
          S2X=(S2X/N-EX**2)*N/(N-1)
          SX=DSQRT(S2X)
          RETURN
          END
