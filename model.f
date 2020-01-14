C*****************************     BMAX
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
C*****************************     CALCO
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
C
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
C
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
          DHELP5 = DSQRT(24.D0*DD) 
          IF(DB.LT.DHELP5.AND.ICASO.EQ.3) DB = DHELP5 
          DA = DCOV(2) * DEXP(DB*DXCOV(2)) / DHELP1
          RETURN
          END  
C
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
          DHELP5 = 2.D0*DD  
          IF(DB.LT.DHELP5.AND.ICASO.EQ.3) DB = DHELP5  
          DA = DCOV(2) * DEXP(DB*DXCOV(2)**2) / DHELP1 
          RETURN
          END  
C
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
C
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
C
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
C
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
C
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
C
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
C
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
C
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
C*****************************     CENTRO
          FUNCTION CENTRO(RAGGIO,DEN,COV,NCOR)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          DIMENSION COV(1)
          CENTRO=0.D0
          DO 10 I=1,NCOR
          IF(COV(I)/DEN.LT.0.5D0) GOTO 20
          CENTRO=(I-0.5D0)*RAGGIO
10        CONTINUE
20        RETURN
          END  
C*****************************     CODA
          SUBROUTINE CODA(I,NPT,II,D,AA,E,BB,CC,IFILE,DATO,LAST,TOLL)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          CHARACTER*1 DATO,RESI,RES2,RES3,AST,BLK,TX,TY,TZ
          CHARACTER*15 CRESI,CRES2,CRES3          
          DIMENSION AA(3,*),BB(3,*),CC(3,*),D(*),E(*),DATO(*) 
          DIMENSION RESI(15),RES2(15),RES3(15)
          EQUIVALENCE (CRESI,RESI(1)),(CRES2,RES2(1)),(CRES3,RES3(1))
          DATA AST/'*'/,BLK/' '/
          SRV=1.D0
          DO 5 J=1,15
          RESI(J)=BLK
          RES2(J)=BLK
          RES3(J)=BLK
5         CONTINUE            
          GO TO (10,20,30,40,50,60,70) II  
10        TX=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TX=AST
          WRITE(6,100) I,NPT,D(1),AA(1,1),SRV,E(1),TX,BB(1,1),SRV, 
     .CC(1,1),SRV    
          IF(IFILE.EQ.1) THEN
          WRITE(CRESI,300) E(1)
          WRITE(11,200) (DATO(J),J=1,LAST),RESI,RES2,RES3
          ENDIF        
          RETURN
20        TY=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TY=AST
          WRITE(6,110) I,NPT,D(1),AA(1,1),SRV,E(1),TY,BB(1,1),SRV, 
     .CC(1,1),SRV  
          IF(IFILE.EQ.1) THEN
          WRITE(CRES2,300) E(2)
          WRITE(11,200) (DATO(J),J=1,LAST),RESI,RES2,RES3
          ENDIF
          RETURN
30        TX=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TX=AST
          TY=BLK
          IF(DABS(E(2))/BB(2,2).GT.TOLL) TY=AST
          WRITE(6,120) I,NPT,D(1),AA(1,1),SRV,AA(1,2), 
     .E(1),TX,BB(1,1),SRV,BB(1,2),CC(1,1),SRV,CC(1,2)  
          WRITE(6,130) D(2),AA(2,2),SRV,E(2),TY,BB(2,2),SRV,CC(2,2),SRV
          IF(IFILE.EQ.1) THEN 
          WRITE(CRESI,300) E(1)
          WRITE(CRES2,300) E(2)
          WRITE(11,200) (DATO(J),J=1,LAST),RESI,RES2,RES3
          ENDIF
          RETURN
40        TZ=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TZ=AST
          WRITE(6,140) I,NPT,D(1),AA(1,1),SRV,E(1),TZ,BB(1,1),SRV, 
     .CC(1,1),SRV  
          IF(IFILE.EQ.1) THEN
          WRITE(CRES3,300) E(3)
          WRITE(11,200) (DATO(J),J=1,LAST),RESI,RES2,RES3
          ENDIF
          RETURN
50        TX=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TX=AST
          TZ=BLK
          IF(DABS(E(2))/BB(2,2).GT.TOLL) TZ=AST
          WRITE(6,150) I,NPT,D(1),AA(1,1),SRV,AA(1,2), 
     .E(1),TX,BB(1,1),SRV,BB(1,2),CC(1,1),SRV,CC(1,2)  
          WRITE(6,160) D(2),AA(2,2),SRV,E(2),TZ,BB(2,2),SRV,CC(2,2),SRV
          IF(IFILE.EQ.1) THEN 
          WRITE(CRESI,300) E(1)
          WRITE(CRES3,300) E(3)
          WRITE(11,200) (DATO(J),J=1,LAST),RESI,RES2,RES3
          ENDIF
          RETURN
60        TY=BLK
          IF(DABS(E(1))/BB(1,1).GT.TOLL) TY=AST
          TZ=BLK
          IF(DABS(E(2))/BB(2,2).GT.TOLL) TZ=AST
          WRITE(6,170) I,NPT,D(1),AA(1,1),SRV,AA(1,2), 
     .E(1),TY,BB(1,1),SRV,BB(1,2),CC(1,1),SRV,CC(1,2)  
          WRITE(6,160) D(2),AA(2,2),SRV,E(2),TZ,BB(2,2),SRV,CC(2,2),SRV
          IF(IFILE.EQ.1) THEN 
          WRITE(CRES2,300) E(2)
          WRITE(CRES3,300) E(3)
          WRITE(11,200) (DATO(J),J=1,LAST),RESI,RES2,RES3
          ENDIF
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
          IF(IFILE.EQ.1) THEN 
          WRITE(CRESI,300) E(1)
          WRITE(CRES2,300) E(2)
          WRITE(CRES3,300) E(3)
          WRITE(11,200) (DATO(J),J=1,LAST),RESI,RES2,RES3
          ENDIF       
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
200       FORMAT(120A1)
300       FORMAT(F15.5)
          END  
C*****************************     CORRP
          SUBROUTINE CORRP1(X,Y,W,IPXW,IPYW,LPXW,LPYW,COV,SRV,NSRV,NCOR,
     .    RAGGIO,ALATO,IWD,NMAX,NWDMAX,SCX,SCY,SCW,
     .    IPTO,IAST,FORMA,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),W(*),IPXW(*),IPYW(*),
     .    LPXW(*),LPYW(*),IWD(*),COV(*),SRV(*),NSRV(*),IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          SRV1=0.D0
          SRV2=1.D0
          ND=0 
          NW=0 
          REWIND 1 
10        READ(1,1000) (DATO(K),K=1,120) 
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
          CALL RBORDA(X,Y,IPXW,IPYW,LPXW,LPYW,NW,ND)
          DO 30 I=1,ND
          W(I)=W(I)-EW
30        CONTINUE
C         VARIABILI: X,Y,W 
          CALL COVAP(X,Y,W,IPXW,IPYW,LPXW,LPYW,NW,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          VALE=CENTRO(RAGGIO,S2W,COV,NCOR) 
          ZERO=FINE(RAGGIO,S2W,COV,NCOR) 
          WRITE(6,1010) NW,EW,SW,COV(1)/S2W,ZERO,VALE 
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
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    5X,'RO-1(X) =',F10.3//
     .    63X,'ZERO(X) =',F10.3,5X,'L-1/2(X) =',F10.3//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(I5,F15.6,F10.3)
1060      FORMAT(/////)
          END  
C
          SUBROUTINE CORRP2(X,Y,U,V,IPXU,IPYU,LPXU,LPYU,
     .    IPXV,IPYV,LPXV,LPYV,COV,SRV,NSRV,
     .    NCOR,RAGGIO,ALATO,IWD,NMAX,NWDMAX,SCX,SCY,SCU,SCV,
     .    IPTO,IAST,FORMA,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),U(*),V(*),IPXU(*),IPYU(*),
     .    LPXU(*),LPYU(*),IPXV(*),IPYV(*),LPXV(*),LPYV(*), 
     .    IWD(*),COV(*),SRV(*),NSRV(*),IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          REWIND 1 
          SRV1=0.D0
          SRV2=1.D0
          ND=0 
          NU=0 
          NV=0 
10        READ(1,1000) (DATO(K),K=1,120)
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
C
20        CALL MEDVA2(U,IPXU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPXV,EV,S2V,SV,NV) 
          CALL RBORDA(X,Y,IPXU,IPYU,LPXU,LPYU,NU,ND)
          CALL RBORDA(X,Y,IPXV,IPYV,LPXV,LPYV,NV,ND)
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
30        CONTINUE 
C         VARIABILI: X,Y,U 
          CALL COVAP(X,Y,U,IPXU,IPYU,LPXU,LPYU,NU,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          VALE=CENTRO(RAGGIO,S2U,COV,NCOR)  
          ZERO=FINE(RAGGIO,S2U,COV,NCOR)  
          WRITE(6,1010) NU,EU,SU,COV(1)/S2U,ZERO,VALE
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
C         VARIABILI: X,Y,V 
          CALL COVAP(X,Y,V,IPXV,IPYV,LPXV,LPYV,NV,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          VALE=CENTRO(RAGGIO,S2V,COV,NCOR)   
          ZERO=FINE(RAGGIO,S2V,COV,NCOR)  
          WRITE(6,1050) NV,EV,SV,COV(1)/S2V,ZERO,VALE
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
C        
          DO 40 I=1,ND 
          U(I)=U(I)+EU 
          V(I)=V(I)+EV 
40        CONTINUE 
          CALL RCORR2(U,V,LPXU,LPXV,EU,EV,SU,SV,S2U,S2V,SUV,RUV,ND)
          DO 50 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
50        CONTINUE 
C         VARIABILI: X,Y,U,V
          CALL CROSSP(X,Y,U,V,IPXU,IPYU,LPXU,LPYU,IPXV,IPYV,LPXV,LPYV, 
     .    NU,NV,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          SPUV=SU*SV
          VALE=CENTRO(RAGGIO,SPUV,COV,NCOR)
          ZERO=FINE(RAGGIO,SPUV,COV,NCOR)
          WRITE(6,1070) NU,EU,SU,NV,EV,SV,RUV,ZERO,VALE  
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
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    5X,'RO-1(X) =',F10.3//
     .    63X,'ZERO(X) =',F10.3,5X,'L-1/2(X) =',F10.3//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y)= ',F10.3,
     .    5X,'RO-1(Y) =',F10.3//
     .    63X,'ZERO(Y) =',F10.3,5X,'L-1/2(Y) =',F10.3//)
1060      FORMAT(//10X,'DELTA',7X,'GAMMA(Y)',5X,'RO(Y)',///,
     .    5X,F10.3,F15.6,F10.3)
1070      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    //5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3, 
     .    5X,'R(X,Y)    =',F10.3//
     .    63X,'ZERO(X,Y) =',F10.3,5X,',L-1/2(X,Y) =',F10.3//)
1080      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Y)',3X,'RO(X,Y)',///5X, 
     .    F10.3,F15.6,F10.3)
1090      FORMAT(I5,F15.6,F10.3,5X,I5,F15.6,F10.3,F10.2)
1100      FORMAT(/////)
          END  
C
          SUBROUTINE CORRP3(X,Y,U,V,W,IPXU,IPYU,LPXU,LPYU, 
     .    IPXV,IPYV,LPXV,LPYV,IPXW,IPYW,LPXW,LPYW,COV,SRV,NSRV,
     .    NCOR,RAGGIO,ALATO,IWD,NMAX,NWDMAX,SCX,SCY,SCU,SCV,SCW,
     .    IPTO,IAST,FORMA,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),U(*),V(*),W(*),IPXU(*),IPYU(*),  
     .    LPXU(*),LPYU(*),IPXV(*),IPYV(*),LPXV(*),LPYV(*), 
     .    IPXW(*),IPYW(*),LPXW(*),LPYW(*),IWD(*),COV(*),SRV(*),NSRV(*),
     .    IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20)
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
C
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
C         VARIABILI: X,Y,U 
          CALL COVAP(X,Y,U,IPXU,IPYU,LPXU,LPYU,NU,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          VALE=CENTRO(RAGGIO,S2U,COV,NCOR)  
          ZERO=FINE(RAGGIO,S2U,COV,NCOR)  
          WRITE(6,1010) NU,EU,SU,COV(1)/S2U,ZERO,VALE
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
C         VARIABILI: X,Y,V 
          CALL COVAP(X,Y,V,IPXV,IPYV,LPXV,LPYV,NV,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          VALE=CENTRO(RAGGIO,S2V,COV,NCOR)  
          ZERO=FINE(RAGGIO,S2V,COV,NCOR)  
          WRITE(6,1050) NV,EV,SV,COV(1)/S2V,ZERO,VALE
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
C        
          DO 40 I=1,ND 
          U(I)=U(I)+EU 
          V(I)=V(I)+EV 
40        CONTINUE 
          CALL RCORR2(U,V,LPXU,LPXV,EU,EV,SU,SV,S2U,S2V,SUV,RUV,ND)
          DO 50 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
50        CONTINUE 
C         VARIABILI: X,Y,U,V
          CALL CROSSP(X,Y,U,V,IPXU,IPYU,LPXU,LPYU,IPXV,IPYV,LPXV,LPYV, 
     .    NU,NV,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          SPUV=SU*SV
          VALE=CENTRO(RAGGIO,SPUV,COV,NCOR)
          ZERO=FINE(RAGGIO,SPUV,COV,NCOR)
          WRITE(6,1090) NU,EU,SU,NV,EV,SV,RUV,ZERO,VALE 
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
C         VARIABILI: X,Y,W 
          CALL COVAP(X,Y,W,IPXW,IPYW,LPXW,LPYW,NW,NCOR,RAGGIO, 
     .    ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI) 
          VALE=CENTRO(RAGGIO,S2W,COV,NCOR) 
          ZERO=FINE(RAGGIO,S2W,COV,NCOR) 
          WRITE(6,1070) NW,EW,SW,COV(1)/S2W,ZERO,VALE
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
C       
          DO 60 I=1,ND 
          U(I)=U(I)+EU 
          W(I)=W(I)+EW 
60        CONTINUE 
          CALL RCORR2(U,W,LPXU,LPXW,EU,EW,SU,SW,S2U,S2W,SUW,RUW,ND)
          DO 70 I=1,ND 
          U(I)=U(I)-EU 
          W(I)=W(I)-EW 
70        CONTINUE 
C         VARIABILI: X,Y,U,W
          CALL CROSSP(X,Y,U,W,IPXU,IPYU,LPXU,LPYU,IPXW,IPYW,LPXW,LPYW, 
     .    NU,NW,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          SPUW=SU*SW
          VALE=CENTRO(RAGGIO,SPUW,COV,NCOR)
          ZERO=FINE(RAGGIO,SPUW,COV,NCOR)
          WRITE(6,1110) NU,EU,SU,NW,EW,SW,RUW,ZERO,VALE  
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
C      
          DO 80 I=1,ND 
          V(I)=V(I)+EV 
          W(I)=W(I)+EW 
80        CONTINUE 
          CALL RCORR2(V,W,LPXV,LPXW,EV,EW,SV,SW,S2V,S2W,SVW,RVW,ND)
          DO 90 I=1,ND 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
90        CONTINUE 
C         VARIABILI: X,Y,V,W
          CALL CROSSP(X,Y,V,W,IPXV,IPYV,LPXV,LPYV,IPXW,IPYW,LPXW,LPYW, 
     .    NV,NW,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          SPVW=SV*SW
          VALE=CENTRO(RAGGIO,SPVW,COV,NCOR)
          ZERO=FINE(RAGGIO,SPVW,COV,NCOR)
          WRITE(6,1130) NV,EV,SV,NW,EW,SW,RVW,ZERO,VALE  
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
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    5X,'RO-1(X) =',F10.3//
     .    63X,'ZERO(X) =',F10.3,5X,'L-1/2(X) =',F10.3//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,
     .    5X,'RO-1(Y) =',F10.3,5X//
     .    63X,'ZERO(Y) =',F10.3,'L-1/2(Y) =',F10.3//)
1060      FORMAT(//10X,'DELTA',7X,'GAMMA(Y)',5X,'RO(Y)',///,
     .    5X,F10.3,F15.6,F10.3)
1070      FORMAT(//5X,'N(Z) =',I5,5X,'M(Z) =',F10.3,5X,'S(Z) =',F10.3,
     .    5X,'RO-1(Z) =',F10.3//
     .    63X,'ZERO(Z) =',F10.3,5X,'L-1/2(Z) =',F10.3//)
1080      FORMAT(//10X,'DELTA',7X,'GAMMA(Z)',5X,'RO(Z)',///,
     .    5X,F10.3,F15.6,F10.3)
1090      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    //5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,5X, 
     .    'R(X,Y)    =',F10.3//
     .    63X,'ZERO(X,Y) =',F10.3,5X,'L-1/2(X,Y) =',F10.3//)
1100      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Y)',3X,'RO(X,Y)',///5X, 
     .    F10.3,F15.6,F10.3)
1110      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    //5X,'N(Z) =',I5,5X,'M(Z) =',F10.3,5X,'S(Z) =',F10.3,5X, 
     .    'R(X,Z)    =',F10.3//
     .    63X,'ZERO(X,Z) =',F10.3,5X,'L-1/2(X,Z) =',F10.3//)
1120      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Z)',3X,'RO(X,Z)',///5X, 
     .    F10.3,F15.6,F10.3)
1130      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y)  =',F10.3,5X,'S(Y)=',F10.3,
     .    //5X,'N(Z) =',I5,5X,'M(Z) =',F10.3,5X,'S(Z) =',F10.3,5X, 
     .    'R(Y,Z)    =',F10.3//
     .    63X,'ZERO(Y,Z) =',F10.3,5X,'L-1/2(Y,Z) =',F10.3//)
1140      FORMAT(//10X,'DELTA',5X,'GAMMA(Y,Z)',3X,'RO(Y,Z)',///5X, 
     .    F10.3,F15.6,F10.3)
1150      FORMAT(I5,F15.6,F10.3,5X,I5,F15.6,F10.3,F10.2)
1160      FORMAT(/////)
          END    
C*****************************     CORRS
          SUBROUTINE CORRS1(X,Y,Z,W,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,
     .    COV,SRV,NSRV,NCOR,
     .    RAGGIO,ALATO,IWD,NMAX,NWDMAX,SCX,SCY,SCZ,SCW,
     .    IPTO,IAST,FORMA,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),Z(*),W(*),IPXW(*),IPYW(*),IPZW(*),
     .    LPXW(*),LPYW(*),LPZW(*),IWD(*),COV(*),SRV(*),NSRV(*),
     .    IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          SRV1=0.D0
          SRV2=1.D0
          ND=0 
          NW=0 
          REWIND 1 
10        READ(1,1000) (DATO(K),K=1,120) 
          READ(CDATO,FORMA) NPT,XS,YS,ZS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY
          Z(ND)=ZS*SCZ 
          W(ND)=WS*SCW 
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NW=NW+1  
          IPXW(NW)=ND  
          IPYW(NW)=ND
          IPZW(NW)=ND  
          ENDIF
          GOTO 10  
C
20        CALL MEDVA2(W,IPXW,EW,S2W,SW,NW) 
          CALL RTORDA(X,Y,Z,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,NW,ND)
          DO 30 I=1,ND
          W(I)=W(I)-EW
30        CONTINUE
C         VARIABILI: X,Y,Z,W 
          CALL COVAS(X,Y,Z,W,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,
     .    NW,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          VALE=CENTRO(RAGGIO,S2W,COV,NCOR)  
          ZERO=FINE(RAGGIO,S2W,COV,NCOR)  
          WRITE(6,1010) NW,EW,SW,COV(1)/S2W,ZERO,VALE
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
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    5X,'RO-1(X) =',F10.3//
     .    63X,'ZERO(X) =',F10.3,5X,'L-1/2(X) =',F10.3//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(I5,F15.6,F10.3)
1060      FORMAT(/////)
          END  
C
          SUBROUTINE CORRS2(X,Y,Z,U,V,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,
     .    IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,COV,SRV,NSRV,
     .    NCOR,RAGGIO,ALATO,IWD,NMAX,NWDMAX,SCX,SCY,SCZ,SCU,SCV,
     .    IPTO,IAST,FORMA,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),Z(*),U(*),V(*),IPXU(*),IPYU(*),
     .    IPZU(*),LPXU(*),LPYU(*),LPZU(*),
     .    IPXV(*),IPYV(*),IPZV(*),LPXV(*),LPYV(*),LPZV(*), 
     .    IWD(*),COV(*),SRV(*),NSRV(*),IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          REWIND 1 
          SRV1=0.D0
          SRV2=1.D0
          ND=0 
          NU=0 
          NV=0 
10        READ(1,1000) (DATO(K),K=1,120)
          READ(CDATO,FORMA) NPT,XS,YS,ZS,US,VS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY
          Z(ND)=ZS*SCZ 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPXU(NU)=ND  
          IPYU(NU)=ND
          IPZU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPXV(NV)=ND  
          IPYV(NV)=ND
          IPZV(NV)=ND  
          ENDIF
          GOTO 10  
C
20        CALL MEDVA2(U,IPXU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPXV,EV,S2V,SV,NV) 
          CALL RTORDA(X,Y,Z,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,NU,ND)
          CALL RTORDA(X,Y,Z,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,NV,ND)
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
30        CONTINUE 
C         VARIABILI: X,Y,Z,U 
          CALL COVAS(X,Y,Z,U,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,
     .    NU,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          VALE=CENTRO(RAGGIO,S2U,COV,NCOR)  
          ZERO=FINE(RAGGIO,S2U,COV,NCOR)  
          WRITE(6,1010) NU,EU,SU,COV(1)/S2U,ZERO,VALE
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
C         VARIABILI: X,Y,Z,V 
          CALL COVAS(X,Y,Z,V,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,
     .    NV,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI) 
          VALE=CENTRO(RAGGIO,S2V,COV,NCOR) 
          ZERO=FINE(RAGGIO,S2V,COV,NCOR) 
          WRITE(6,1050) NV,EV,SV,COV(1)/S2V,ZERO,VALE
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
C     
          DO 40 I=1,ND 
          U(I)=U(I)+EU 
          V(I)=V(I)+EV 
40        CONTINUE 
          CALL RCORR2(U,V,LPXU,LPXV,EU,EV,SU,SV,S2U,S2V,SUV,RUV,ND)
          DO 50 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
50        CONTINUE 
C         VARIABILI: X,Y,Z,U,V
          CALL CROSSS(X,Y,Z,U,V,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,
     .    IPXV,IPYV,IPZV,LPXV,LPYV,LPZV, 
     .    NU,NV,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          SPUV=SU*SV
          VALE=CENTRO(RAGGIO,SPUV,COV,NCOR)
          ZERO=FINE(RAGGIO,SPUV,COV,NCOR)
          WRITE(6,1070) NU,EU,SU,NV,EV,SV,RUV,ZERO,VALE  
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
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    5X,'RO-1(X) =',F10.3//
     .    63X,'ZERO(X) =',F10.3,5X,'L-1/2(X) =',F10.3//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,
     .    5X,'RO-1(Y) =',F10.3//
     .    63X,'ZERO(Y) =',F10.3,5X,'L-1/2(Y) =',F10.3//)
1060      FORMAT(//10X,'DELTA',7X,'GAMMA(Y)',5X,'RO(Y)',///,
     .    5X,F10.3,F15.6,F10.3)
1070      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    //5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,5X, 
     .    'R(X,Y)    =',F10.3//
     .    63X,'ZERO(X,Y) =',F10.3,5X,'L-1/2(X,Y) =',F10.3//)
1080      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Y)',3X,'RO(X,Y)',///5X, 
     .    F10.3,F15.6,F10.3)
1090      FORMAT(I5,F15.6,F10.3,5X,I5,F15.6,F10.3,F10.2)
1100      FORMAT(/////)
          END  
C
          SUBROUTINE CORRS3(X,Y,Z,U,V,W,IPXU,IPYU,IPZU,
     .    LPXU,LPYU,LPZU,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,
     .    IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,COV,SRV,NSRV,
     .    NCOR,RAGGIO,ALATO,IWD,NMAX,NWDMAX,SCX,SCY,SCZ,SCU,SCV,SCW,
     .    IPTO,IAST,FORMA,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),Z(*),U(*),V(*),W(*),IPXU(*),IPYU(*),  
     .    IPZU(*),LPXU(*),LPYU(*),LPZU(*),IPXV(*),IPYV(*),IPZV(*),
     .    LPXV(*),LPYV(*),LPZV(*),IPXW(*),IPYW(*),IPZW(*),
     .    LPXW(*),LPYW(*),LPZW(*),IWD(*),COV(*),SRV(*),NSRV(*),
     .    IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20)
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
          READ(CDATO,FORMA) NPT,XS,YS,ZS,US,VS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY 
          Z(ND)=ZS*SCZ
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          W(ND)=WS*SCW 
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPXU(NU)=ND  
          IPYU(NU)=ND
          IPZU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPXV(NV)=ND  
          IPYV(NV)=ND
          IPZV(NV)=ND  
          ENDIF
          IF(DATO(IPTO(5)).EQ.PTO.AND.DATO(IAST(3)).NE.AST) THEN
          NW=NW+1  
          IPXW(NW)=ND  
          IPYW(NW)=ND
          IPZW(NW)=ND  
          ENDIF
          GOTO 10  
C
20        CALL MEDVA2(U,IPXU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPXV,EV,S2V,SV,NV) 
          CALL MEDVA2(W,IPXW,EW,S2W,SW,NW) 
          CALL RTORDA(X,Y,Z,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,NU,ND)
          CALL RTORDA(X,Y,Z,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,NV,ND)
          CALL RTORDA(X,Y,Z,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,NW,ND)
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
30        CONTINUE 
C         VARIABILI: X,Y,Z,U 
          CALL COVAS(X,Y,Z,U,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,
     .    NU,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI) 
          VALE=CENTRO(RAGGIO,S2U,COV,NCOR) 
          ZERO=FINE(RAGGIO,S2U,COV,NCOR) 
          WRITE(6,1010) NU,EU,SU,COV(1)/S2U,ZERO,VALE
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
C         VARIABILI: X,Y,Z,V 
          CALL COVAS(X,Y,Z,V,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,
     .    NV,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          VALE=CENTRO(RAGGIO,S2V,COV,NCOR)  
          ZERO=FINE(RAGGIO,S2V,COV,NCOR)  
          WRITE(6,1050) NV,EV,SV,COV(1)/S2V,ZERO,VALE
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
C       
          DO 40 I=1,ND 
          U(I)=U(I)+EU 
          V(I)=V(I)+EV 
40        CONTINUE 
          CALL RCORR2(U,V,LPXU,LPXV,EU,EV,SU,SV,S2U,S2V,SUV,RUV,ND)
          DO 50 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
50        CONTINUE 
C         VARIABILI: X,Y,Z,U,V
          CALL CROSSS(X,Y,Z,U,V,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,
     .    IPXV,IPYV,IPZV,LPXV,LPYV,LPZV, 
     .    NU,NV,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          SPUV=SU*SV
          VALE=CENTRO(RAGGIO,SPUV,COV,NCOR)
          ZERO=FINE(RAGGIO,SPUV,COV,NCOR)
          WRITE(6,1090) NU,EU,SU,NV,EV,SV,RUV,ZERO,VALE  
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
C         VARIABILI: X,Y,Z,W 
          CALL COVAS(X,Y,Z,W,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,
     .    NW,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          VALE=CENTRO(RAGGIO,S2W,COV,NCOR)         
          ZERO=FINE(RAGGIO,S2W,COV,NCOR)  
          WRITE(6,1070) NW,EW,SW,COV(1)/S2W,ZERO,VALE
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
C      
          DO 60 I=1,ND 
          U(I)=U(I)+EU 
          W(I)=W(I)+EW 
60        CONTINUE 
          CALL RCORR2(U,W,LPXU,LPXW,EU,EW,SU,SW,S2U,S2W,SUW,RUW,ND)
          DO 70 I=1,ND 
          U(I)=U(I)-EU 
          W(I)=W(I)-EW 
70        CONTINUE 
C         VARIABILI: X,Y,Z,U,W
          CALL CROSSS(X,Y,Z,U,W,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,
     .    IPXW,IPYW,IPZW,LPXW,LPYW,LPZW, 
     .    NU,NW,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          SPUW=SU*SW
          VALE=CENTRO(RAGGIO,SPUW,COV,NCOR)
          ZERO=FINE(RAGGIO,SPUW,COV,NCOR)
          WRITE(6,1110) NU,EU,SU,NW,EW,SW,RUW,ZERO,VALE 
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
C        
          DO 80 I=1,ND 
          V(I)=V(I)+EV 
          W(I)=W(I)+EW 
80        CONTINUE 
          CALL RCORR2(V,W,LPXV,LPXW,EV,EW,SV,SW,S2V,S2W,SVW,RVW,ND)
          DO 90 I=1,ND 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
90        CONTINUE 
C         VARIABIL: X,Y,Z,V,W
          CALL CROSSS(X,Y,Z,V,W,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,
     .    IPXW,IPYW,IPZW,LPXW,LPYW,LPZW, 
     .    NV,NW,NCOR,RAGGIO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          SPVW=SV*SW
          VALE=CENTRO(RAGGIO,SPVW,COV,NCOR)
          ZERO=FINE(RAGGIO,SPVW,COV,NCOR)
          WRITE(6,1130) NV,EV,SV,NW,EW,SW,RVW,ZERO,VALE 
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
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    5X,'RO-1(X) =',F10.3//
     .    63X,'ZERO(X) =',F10.3,5X,'L-1/2(X) =',F10.3//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,
     .    5X,'RO-1(Y) =',F10.3//
     .    63X,'ZERO(Y) =',F10.3,5X,'L-1/2(Y) =',F10.3//)
1060      FORMAT(//10X,'DELTA',7X,'GAMMA(Y)',5X,'RO(Y)',///,
     .    5X,F10.3,F15.6,F10.3)
1070      FORMAT(//5X,'N(Z) =',I5,5X,'M(Z) =',F10.3,5X,'S(Z) =',F10.3,
     .    5X,'RO-1(Z) =',F10.3//
     .    63X,'ZERO(Z) =',F10.3,5X,'L-1/2(Z) =',F10.3//)
1080      FORMAT(//10X,'DELTA',7X,'GAMMA(Z)',5X,'RO(Z)',///,
     .    5X,F10.3,F15.6,F10.3)
1090      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    //5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,5X, 
     .    'R(X,Y)    =',F10.3//
     .    63X,'ZERO(X,Y) =',F10.3,5X,'L-1/2(X,Y) =',F10.3//)
1100      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Y)',3X,'RO(X,Y)',///5X, 
     .    F10.3,F15.6,F10.3)
1110      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    //5X,'N(Z) =',I5,5X,'M(Z) =',F10.3,5X,'S(Z) =',F10.3,5X, 
     .    'R(X,Z)    =',F10.3//
     .    63X,'ZERO(X,Z) =',F10.3,5X,'L-1/2(X,Z) =',F10.3//)
1120      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Z)',3X,'RO(X,Z)',///5X, 
     .    F10.3,F15.6,F10.3)
1130      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,
     .    //5X,'N(Z) =',I5,5X,'M(Z) =',F10.3,5X,'S(Z) =',F10.3,5X, 
     .    'R(Y,Z)    =',F10.3//
     .    63X,'ZERO(Y,Z) =',F10.3,5X,'L-1/2(Y,Z) =',F10.3//)
1140      FORMAT(//10X,'DELTA',5X,'GAMMA(Y,Z)',3X,'RO(Y,Z)',///5X, 
     .    F10.3,F15.6,F10.3)
1150      FORMAT(I5,F15.6,F10.3,5X,I5,F15.6,F10.3,F10.2)
1160      FORMAT(/////)
          END  
C*****************************     CORRT
          SUBROUTINE CORRT1(Z,W,IPZW,COV,SRV,NSRV,NCOR,RAGGIO,ALATO,
     .    NMAX,SCZ,SCW,IPTO,IAST,FORMA,IFILE,ILI)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION Z(*),W(*),IPZW(*),COV(*),SRV(*),NSRV(*),
     .    IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20)
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          SRV1=0.D0
          SRV2=1.D0
          REWIND 1 
          ND=0 
          NW=0 
10        READ(1,1000) (DATO(K),K=1,120)
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
C
20        CALL MEDVA2(W,IPZW,EW,S2W,SW,NW) 
          CALL RIORDA(Z,IPZW,NW,ND,1)  
          DO 30 I=1,ND 
          W(I)=W(I)-EW 
30        CONTINUE 
C         VARIABILI: Z,W
          CALL COVAT(Z,W,IPZW,NW,NCOR,RAGGIO,ALATO,COV,SRV,
     .    NSRV,NMAX,ILI)  
          VALE=CENTRO(RAGGIO,S2W,COV,NCOR)
          ZERO=FINE(RAGGIO,S2W,COV,NCOR)
          WRITE(6,1010) NW,EW,SW,COV(1)/S2W,ZERO,VALE
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
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    5X,'RO-1(X) =',F10.3//
     .    63X,'ZERO(X) =',F10.3,5X,'L-1/2(X) =',F10.3//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(I5,F15.6,F10.3)
1060      FORMAT(/////)
          END  
C
          SUBROUTINE CORRT2(Z,U,V,IPZU,LPZU,IPZV,LPZV,COV,SRV,NSRV,NCOR,
     .    RAGGIO,ALATO,NMAX,SCZ,SCU,SCV,IPTO,IAST,FORMA,IFILE,ILI) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION Z(*),U(*),V(*),IPZU(*),LPZU(*),IPZV(*),LPZV(*),
     .    COV(*),SRV(*),NSRV(*),IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20)  
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          SRV1=0.D0
          SRV2=1.D0
          ND=0 
          NU=0 
          NV=0 
          REWIND 1 
10        READ(1,1000) (DATO(K),K=1,120)
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
C
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
C         VARIABILI: Z,U
          CALL COVAT(Z,U,IPZU,NU,NCOR,RAGGIO,ALATO,COV,SRV,
     .    NSRV,NMAX,ILI)  
          VALE=CENTRO(RAGGIO,S2U,COV,NCOR)
          ZERO=FINE(RAGGIO,S2U,COV,NCOR)
          WRITE(6,1010) NU,EU,SU,COV(1)/S2U,ZERO,VALE
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
C         VARIABILI: Z,V
          CALL COVAT(Z,V,IPZV,NV,NCOR,RAGGIO,ALATO,COV,SRV,
     .    NSRV,NMAX,ILI)  
          VALE=CENTRO(RAGGIO,S2V,COV,NCOR)
          ZERO=FINE(RAGGIO,S2V,COV,NCOR)
          WRITE(6,1050) NV,EV,SV,COV(1)/S2V,ZERO,VALE
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
C      
          DO 40 I=1,ND 
          U(I)=U(I)+EU 
          V(I)=V(I)+EV 
40        CONTINUE 
          CALL RCORR2(U,V,LPZU,LPZV,EU,EV,SU,SV,S2U,S2V,SUV,RUV,ND)
          DO 50 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
50        CONTINUE 
C         VARIABILI: Z,U,V
          CALL CROSST(Z,U,V,IPZU,LPZU,IPZV,LPZV,NU,NV, 
     .    NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX,ILI)
          SPUV=SU*SV
          VALE=CENTRO(RAGGIO,SPUV,COV,NCOR) 
          ZERO=FINE(RAGGIO,SPUV,COV,NCOR) 
          WRITE(6,1070) NU,EU,SU,NV,EV,SV,RUV,ZERO,VALE  
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
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    5X,'RO-1(X) =',F10.3//
     .    63X,'ZERO(X) =',F10.3,5X,'L-1/2(X) =',F10.3//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,
     .    5X,'RO-1(Y) =',F10.3//
     .    63X,'ZERO(Y) =',F10.3,5X,'L-1/2(Y) =',F10.3//)
1060      FORMAT(//10X,'DELTA',7X,'GAMMA(Y)',5X,'RO(Y)',///,
     .    5X,F10.3,F15.6,F10.3)
1070      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    //5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,5X, 
     .    'R(X,Y)    =',F10.3//
     .    63X,'ZERO(X,Y) =',F10.3,5X,'L-1/2(X,Y) =',F10.3//)
1080      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Y)',3X,'RO(X,Y)',///5X, 
     .    F10.3,F15.6,F10.3)
1090      FORMAT(I5,F15.6,F10.3,5X,I5,F15.6,F10.3,F10.2)
1100      FORMAT(/////)
          END  
C
          SUBROUTINE CORRT3(Z,U,V,W,IPZU,LPZU,IPZV,LPZV,IPZW,LPZW, 
     .    COV,SRV,NSRV,NCOR,RAGGIO,ALATO,NMAX, 
     .    SCZ,SCU,SCV,SCW,IPTO,IAST,FORMA,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION Z(*),U(*),V(*),W(*),IPZU(*),LPZU(*),
     .    IPZV(*),LPZV(*),IPZW(*),LPZW(*), 
     .    COV(*),SRV(*),NSRV(*),IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20)  
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
C
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
C         VARIABILI: Z,U
          CALL COVAT(Z,U,IPZU,NU,NCOR,RAGGIO,ALATO,COV,SRV,
     .    NSRV,NMAX,ILI)  
          VALE=CENTRO(RAGGIO,S2U,COV,NCOR)
          ZERO=FINE(RAGGIO,S2U,COV,NCOR)
          WRITE(6,1010) NU,EU,SU,COV(1)/S2U,ZERO,VALE
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
C         VARIABILI: Z,V
          CALL COVAT(Z,V,IPZV,NV,NCOR,RAGGIO,ALATO,COV,SRV,          
     .    NSRV,NMAX,ILI)  
          VALE=CENTRO(RAGGIO,S2V,COV,NCOR)
          ZERO=FINE(RAGGIO,S2V,COV,NCOR)
          WRITE(6,1050) NV,EV,SV,COV(1)/S2V,ZERO,VALE
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
C        
          DO 40 I=1,ND 
          U(I)=U(I)+EU 
          V(I)=V(I)+EV 
40        CONTINUE 

          CALL RCORR2(U,V,LPZU,LPZV,EU,EV,SU,SV,S2U,S2V,SUV,RUV,ND)
          DO 50 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
50        CONTINUE 
C         VARIABILI: Z,U,V 
          CALL CROSST(Z,U,V,IPZU,LPZU,IPZV,LPZV,NU,NV, 
     .    NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX,ILI)
          SPUV=SU*SV
          VALE=CENTRO(RAGGIO,SPUV,COV,NCOR) 
          ZERO=FINE(RAGGIO,SPUV,COV,NCOR) 
          WRITE(6,1090) NU,EU,SU,NV,EV,SV,RUV,ZERO,VALE  
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
C         VARIABILI: Z,W
          CALL COVAT(Z,W,IPZW,NW,NCOR,RAGGIO,ALATO,COV,SRV,          
     .    NSRV,NMAX,ILI)  
          VALE=CENTRO(RAGGIO,S2W,COV,NCOR)
          ZERO=FINE(RAGGIO,S2W,COV,NCOR)
          WRITE(6,1070) NW,EW,SW,COV(1)/S2W,ZERO,VALE
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
C       
          DO 60 I=1,ND 
          U(I)=U(I)+EU 
          W(I)=W(I)+EW 
60        CONTINUE 
          CALL RCORR2(U,W,LPZU,LPZW,EU,EW,SU,SW,S2U,S2W,SUW,RUW,ND)
          DO 70 I=1,ND 
          U(I)=U(I)-EU 
          W(I)=W(I)-EW 
70        CONTINUE 
C         VARIABILI: Z,U,W 
          CALL CROSST(Z,U,W,IPZU,LPZU,IPZW,LPZW,NU,NW, 
     .    NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX,ILI)
          SPUW=SU*SW
          VALE=CENTRO(RAGGIO,SPUW,COV,NCOR) 
          ZERO=FINE(RAGGIO,SPUW,COV,NCOR) 
          WRITE(6,1110) NU,EU,SU,NW,EW,SW,RUW,ZERO,VALE 
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
C       
          DO 80 I=1,ND 
          V(I)=V(I)+EV 
          W(I)=W(I)+EW 
80        CONTINUE 
          CALL RCORR2(V,W,LPZV,LPZW,EV,EW,SV,SW,S2V,S2W,SVW,RVW,ND)
          DO 90 I=1,ND 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
90        CONTINUE 
C         VARIABILI: Z,V,W 
          CALL CROSST(Z,V,W,IPZV,LPZV,IPZW,LPZW,NV,NW, 
     .    NCOR,RAGGIO,ALATO,COV,SRV,NSRV,NMAX,ILI) 
          SPVW=SV*SW
          VALE=CENTRO(RAGGIO,SPVW,COV,NCOR)
          ZERO=FINE(RAGGIO,SPVW,COV,NCOR)
          WRITE(6,1130) NV,EV,SV,NW,EW,SW,RVW,ZERO,VALE  
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
          RETURN
1000      FORMAT(120A1)
1010      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    5X,'RO-1(X) =',F10.3//
     .    63X,'ZERO(X) =',F10.3,5X,'L-1/2(X) =',F10.3//)
1020      FORMAT(//10X,'DELTA',7X,'GAMMA(X)',5X,'RO(X)',///,
     .    5X,F10.3,F15.6,F10.3)
1030      FORMAT(5X,F10.3,F15.6,F10.3) 
1040      FORMAT(F10.3,F15.6,F10.3)
1050      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,
     .    5X,'RO-1(Y) =',F10.3//
     .    63X,'ZERO(Y) =',F10.3,5X,'L-1/2(Y) =',F10.3//)
1060      FORMAT(//10X,'DELTA',7X,'GAMMA(Y)',5X,'RO(Y)',///,
     .    5X,F10.3,F15.6,F10.3)
1070      FORMAT(//5X,'N(Z) =',I5,5X,'M(Z) =',F10.3,5X,'S(Z) =',F10.3,
     .    5X,'RO-1(Z) ='F10.3//
     .    63X,'ZERO(Z) =',F10.3,5X,'L-/2(Z) =',F10.3//)
1080      FORMAT(//10X,'DELTA',7X,'GAMMA(Z)',5X,'RO(Z)',///,
     .    5X,F10.3,F15.6,F10.3)
1090      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3, 
     .    //5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,5X,
     .    'R(X,Y)    =',F10.3//
     .    63X,'ZERO(X,Y) =',F10.3,5X,'L-1/2(X,Y) =',F10.3//)
1100      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Y)',3X,'RO(X,Y)',///5X, 
     .    F10.3,F15.6,F10.3)
1110      FORMAT(//5X,'N(X) =',I5,5X,'M(X) =',F10.3,5X,'S(X) =',F10.3,
     .    //5X,'N(Z) =',I5,5X,'M(Z) =',F10.3,5X,'S(Z) =',F10.3,5X, 
     .    'R(X,Z)    =',F10.3//
     .    63X,'ZERO(X,Z) =',F10.3,5X,'L-1/2(X,Z) =',F10.3//)
1120      FORMAT(//10X,'DELTA',5X,'GAMMA(X,Z)',3X,'RO(X,Z)',///5X, 
     .    F10.3,F15.6,F10.3)
1130      FORMAT(//5X,'N(Y) =',I5,5X,'M(Y) =',F10.3,5X,'S(Y) =',F10.3,
     .    //5X,'N(Z) =',I5,5X,'M(Z) =',F10.3,5X,'S(Z) =',F10.3,5X,
     .    'R(Y,Z) =',F10.3//
     .    63X,'ZERO(Y,Z) =',F10.3,5X,'L-1/2(Y,Z) =',F10.3//)
1140      FORMAT(//10X,'DELTA',5X,'GAMMA(Y,Z)',3X,'RO(Y,Z)',///5X, 
     .    F10.3,F15.6,F10.3)
1150      FORMAT(I5,F15.6,F10.3,5X,I5,F15.6,F10.3,F10.2)
1160      FORMAT(/////)
          END  
C*****************************     COVAP
          SUBROUTINE COVAP(X,Y,W,IPX,IPY,LPX,LPY,N,NCOR,RGGO,ALATO,COV,
     .    SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(*),Y(*),W(*),IPX(*),IPY(*),LPX(*),LPY(*),COV(*), 
     .    IWD(*),SRV(*),NSRV(*)
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
          CALL WINDOW(XC,YC,ALATO,ALATO,X,Y,IPX,IPY,LPX,LPY,IWD,NMAX,  
     .    NWDMAX,N,NWD,ILI)
          DO 30 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+W(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          SRV(1)=SRV(1)-W(LL)  
          NSRV(1)=NSRV(1)-1
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+W(LL)*SRV(K)/NK
40        CONTINUE 
50        CONTINUE
          DO 60 K=1,NCOR
          COV(K)=COV(K)/N 
60        CONTINUE 
          RETURN
          END      
C*****************************     COROP
          SUBROUTINE COROP(X,Y,W,IPX,IPY,LPX,LPY,N,NCOR,RGGO,ALATO,COV,
     .    SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(*),Y(*),W(*),IPX(*),IPY(*),LPX(*),LPY(*),COV(*), 
     .    IWD(*),SRV(*),NSRV(*)
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
          CALL WINDOW(XC,YC,ALATO,ALATO,X,Y,IPX,IPY,LPX,LPY,IWD,NMAX,  
     .    NWDMAX,N,NWD,ILI)
          DO 30 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+W(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          SRV(1)=SRV(1)-W(LL)  
          NSRV(1)=NSRV(1)-1
          DO 35 K=2,NCOR
          SRV(K)=SRV(K-1)+SRV(K)
          NSRV(K)=NSRV(K-1)+NSRV(K)
35        CONTINUE
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+W(LL)*SRV(K)/NK
40        CONTINUE 
50        CONTINUE 
          DO 60 K=1,NCOR
          COV(K)=COV(K)/N  
60        CONTINUE 
          RETURN
          END      
C*****************************     COVAS
          SUBROUTINE COVAS(X,Y,Z,W,IPX,IPY,IPZ,LPX,LPY,LPZ,N,NCOR,
     .    RGGO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(*),Y(*),Z(*),W(*),IPX(*),IPY(*),IPZ(*),
     .    LPX(*),LPY(*),LPZ(*),COV(*),IWD(*),SRV(*),NSRV(*)
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
          ZC=Z(LL)
          CALL POCKET(XC,YC,ZC,ALATO,ALATO,ALATO,X,Y,Z,IPX,IPY,IPZ,
     .    LPX,LPY,LPZ,IWD,NMAX,NWDMAX,N,NWD,ILI)
          DO 30 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2+
     .    (Z(LL)-Z(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+W(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          SRV(1)=SRV(1)-W(LL)  
          NSRV(1)=NSRV(1)-1
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+W(LL)*SRV(K)/NK
40        CONTINUE 
50        CONTINUE 
          DO 60 K=1,NCOR
          COV(K)=COV(K)/N  
60        CONTINUE 
          RETURN
          END 
C*****************************     COROS
          SUBROUTINE COROS(X,Y,Z,W,IPX,IPY,IPZ,LPX,LPY,LPZ,N,NCOR,
     .    RGGO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(*),Y(*),Z(*),W(*),IPX(*),IPY(*),IPZ(*),
     .    LPX(*),LPY(*),LPZ(*),COV(*),IWD(*),SRV(*),NSRV(*)
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
          ZC=Z(LL)
          CALL POCKET(XC,YC,ZC,ALATO,ALATO,ALATO,X,Y,Z,IPX,IPY,IPZ,
     .    LPX,LPY,LPZ,IWD,NMAX,NWDMAX,N,NWD,ILI)
          DO 30 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2+
     .    (Z(LL)-Z(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+W(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          SRV(1)=SRV(1)-W(LL)  
          NSRV(1)=NSRV(1)-1
          DO 35 K=2,NCOR
          SRV(K)=SRV(K-1)+SRV(K)
          NSRV(K)=NSRV(K-1)+NSRV(K)
35        CONTINUE          
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+W(LL)*SRV(K)/NK
40        CONTINUE 
50        CONTINUE 
          DO 60 K=1,NCOR
          COV(K)=COV(K)/N  
60        CONTINUE
          RETURN
          END  
C*****************************     COVAT
          SUBROUTINE COVAT(Z,W,IPZ,N,NCOR,RGGO,ALATO,COV,SRV,NSRV, 
     .    NMAX,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION Z(1),W(1),IPZ(1),COV(1),SRV(1),NSRV(1) 
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
          SRV(K)=SRV(K)+W(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          SRV(1)=SRV(1)-W(LL)  
          NSRV(1)=NSRV(1)-1
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+W(LL)*SRV(K)/NK
40        CONTINUE 
50        CONTINUE 
          DO 60 K=1,NCOR
          COV(K)=COV(K)/N  
60        CONTINUE 
          RETURN
          END  
C*****************************     COROT
          SUBROUTINE COROT(Z,W,IPZ,N,NCOR,RGGO,ALATO,COV,SRV,NSRV, 
     .    NMAX,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION Z(1),W(1),IPZ(1),COV(1),SRV(1),NSRV(1) 
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
          SRV(K)=SRV(K)+W(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          SRV(1)=SRV(1)-W(LL)  
          NSRV(1)=NSRV(1)-1
          DO 35 K=2,NCOR
          SRV(K)=SRV(K-1)+SRV(K)
          NSRV(K)=NSRV(K-1)+NSRV(K)
35        CONTINUE 
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+W(LL)*SRV(K)/NK
40        CONTINUE 
50        CONTINUE 
          DO 60 K=1,NCOR
          COV(K)=COV(K)/N  
60        CONTINUE 
          RETURN
          END  
C*****************************     CROSSP
          SUBROUTINE CROSSP(X,Y,U,V,IPX,IPY,LPX,LPY,IPX2,IPY2,LPX2,LPY2,
     .    N,M,NCOR,RGGO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),U(1),V(1),IPX(1),IPY(1),LPX(1),LPY(1),
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
          CALL WINDOW(XC,YC,ALATO,ALATO,X,Y,IPX2,IPY2,LPX2,LPY2,IWD,
     .    NMAX,NWDMAX,M,NWD,ILI)
          DO 30 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+V(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          IF(LPX2(LL).NE.0) THEN
          SRV(1)=SRV(1)-V(LL)  
          NSRV(1)=NSRV(1)-1
          ENDIF
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+U(LL)*SRV(K)/NK
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
          CALL WINDOW(XC,YC,ALATO,ALATO,X,Y,IPX,IPY,LPX,LPY,IWD,
     .    NMAX,NWDMAX,N,NWD,ILI)
          DO 70 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+U(L)
          NSRV(K)=NSRV(K)+1
70        CONTINUE 
          IF(LPX(LL).NE.0) THEN
          SRV(1)=SRV(1)-U(LL)  
          NSRV(1)=NSRV(1)-1
          ENDIF
          DO 80 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+V(LL)*SRV(K)/NK
80        CONTINUE 
90        CONTINUE 
          DO 100 K=1,NCOR  
          COV(K)=COV(K)/(N+M)  
100       CONTINUE 
          RETURN
          END  
C*****************************     CROSSS
          SUBROUTINE CROSSS(X,Y,Z,U,V,IPX,IPY,IPZ,LPX,LPY,LPZ,
     .    IPX2,IPY2,IPZ2,LPX2,LPY2,LPZ2,
     .    N,M,NCOR,RGGO,ALATO,COV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),Z(1),U(1),V(1),IPX(1),IPY(1),IPZ(1),
     .    LPX(1),LPY(1),LPZ(1),IPX2(1),IPY2(1),IPZ2(1),
     .    LPX2(1),LPY2(1),LPZ2(1),COV(1),SRV(1),NSRV(1),IWD(1) 
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
          ZC=Z(LL) 
          CALL POCKET(XC,YC,ZC,ALATO,ALATO,ALATO,X,Y,Z,
     .    IPX2,IPY2,IPZ2,LPX2,LPY2,LPZ2,IWD,NMAX,NWDMAX,M,NWD,ILI)
          DO 30 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2+
     .    (Z(LL)-Z(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+V(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          IF(LPX2(LL).NE.0) THEN
          SRV(1)=SRV(1)-V(LL)  
          NSRV(1)=NSRV(1)-1
          ENDIF
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+U(LL)*SRV(K)/NK
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
          ZC=Z(LL)
          CALL POCKET(XC,YC,ZC,ALATO,ALATO,ALATO,X,Y,Z,
     .    IPX,IPY,IPZ,LPX,LPY,LPZ,IWD,NMAX,NWDMAX,N,NWD,ILI)
          DO 70 J=1,NWD
          L=IWD(J) 
          DIJ=DSQRT((X(LL)-X(L))**2+(Y(LL)-Y(L))**2+
     .    (Z(LL)-Z(L))**2)
          K=INT(DIJ/RGGO)+1
          SRV(K)=SRV(K)+U(L)
          NSRV(K)=NSRV(K)+1
70        CONTINUE 
          IF(LPX(LL).NE.0) THEN
          SRV(1)=SRV(1)-U(LL)  
          NSRV(1)=NSRV(1)-1
          ENDIF
          DO 80 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+V(LL)*SRV(K)/NK
80        CONTINUE 
90        CONTINUE 
          DO 100 K=1,NCOR  
          COV(K)=COV(K)/(N+M)  
100       CONTINUE 
          RETURN
          END  
C*****************************     CROSST
          SUBROUTINE CROSST(Z,U,V,IPZ,LPZ,IPZ2,LPZ2,N,M,NCOR,  
     .    RGGO,ALATO,COV,SRV,NSRV,NMAX,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION Z(1),U(1),V(1),IPZ(1),IPZ2(1),LPZ(1),LPZ2(1),COV(1),
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
          SRV(K)=SRV(K)+V(L)
          NSRV(K)=NSRV(K)+1
30        CONTINUE 
          IF(LPZ2(LL).NE.0) THEN
          SRV(1)=SRV(1)-V(LL)  
          NSRV(1)=NSRV(1)-1
          ENDIF
          DO 40 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+U(LL)*SRV(K)/NK
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
          SRV(K)=SRV(K)+U(L)
          NSRV(K)=NSRV(K)+1
70        CONTINUE 
          IF(LPZ(LL).NE.0) THEN
          SRV(1)=SRV(1)-U(LL)  
          NSRV(1)=NSRV(1)-1
          ENDIF
          DO 80 K=1,NCOR
          NK=NSRV(K)
          IF(NK.GT.0) COV(K)=COV(K)+V(LL)*SRV(K)/NK
80        CONTINUE 
90        CONTINUE 
          DO 100 K=1,NCOR  
          COV(K)=COV(K)/(N+M)  
100       CONTINUE 
          RETURN
          END  
C*****************************     FINE
          FUNCTION FINE(RAGGIO,DEN,COV,NCOR)  
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          DIMENSION COV(1)
          FINE=9999.999D0
          DO 10 I=NCOR,1,-1
          IF(DABS(COV(I))/DEN.GT.0.1D0) GOTO 20
          FINE=(I-0.5D0)*RAGGIO 
10        CONTINUE 
20        RETURN
          END  
C*****************************     GASS1
          FUNCTION GASS1(DELTA,XCOVT,II,JJ,L,K,DA,DB,DC,SN,NCASO,IC)  
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
C*****************************     GASS2
          FUNCTION GASS2(DELTA,XCOVT,II,JJ,L,K,DA,DB,DC,SN,NCASO,IC)  
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
C*****************************     GASS3
          FUNCTION GASS3(DELTA,XCOVT,II,JJ,L,K,DA,DB,DC,SN,NCASO,IC)  
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
          GO TO(10,20,30,40,42,43,45,45,50,60),ICC
10        SRV=DDA*DEXP(-DDB*DELTA) 
          GASS3=SRV*SCAL3(DELTA/XCOVT) 
          RETURN
20        SRV=DDA*DEXP(-DDB*DELTA**2)  
          GASS3=SRV*SCAL3(DELTA/XCOVT) 
          RETURN
30        SRV=DDA*(1.D0-DDC*DELTA**2)*DEXP(-DDB*DELTA) 
          GASS3=SRV*SCAL3(DELTA/XCOVT) 
          RETURN
40        SRV=DDA*(1.D0-DDC*DELTA**2)*DEXP(-DDB*DELTA**2)  
          GASS3=SRV*SCAL3(DELTA/XCOVT) 
          RETURN
42        SRV=DDA*DEXP(-DDB*DELTA)*DSIN(DDC*DELTA)/(DDC*DELTA) 
          GASS3=SRV*SCAL3(DELTA/XCOVT) 
          RETURN
43        SRV=DDA*DEXP(-DDB*DELTA**2)*DSIN(DDC*DELTA)/(DDC*DELTA)  
          GASS3=SRV*SCAL3(DELTA/XCOVT) 
          RETURN
45        SRV=DDA*DEXP(-DDB*DELTA)*DCOS(DDC*DELTA) 
          GASS3=SRV*SCAL1(DELTA/XCOVT)
          RETURN
50        SRV=DDA*DEXP(-DDB*DELTA)*DFJ1(DDC*DELTA)/(DDC*DELTA) 
          GASS3=SRV*SCAL3(DELTA/XCOVT) 
          RETURN
60        SRV=DDA*DEXP(-DDB*DELTA**2)*DFJ1(DDC*DELTA)/(DDC*DELTA)  
          GASS3=SRV*SCAL3(DELTA/XCOVT) 
          RETURN
          ELSE 
          GASS3=DDA+SSN
          ENDIF
          RETURN
          END  
C*****************************     GRACO
          SUBROUTINE GRACO(COV,TCOV,STP,NCOV)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 LB,PUN,AST,PIU,BL
          DIMENSION LB(101),A(11),COV(1),TCOV(1)
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
          WRITE(6,110) (LB(K),K=1,101)  
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
C*****************************     IDIF2
      FUNCTION IDIF2(COV)  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      DIMENSION COV(*) 
      D=COV(2)-2.D0*COV(3)+COV(4)  
      IF(D.GE.0.D0) THEN
      IDIF2=1  
      ELSE 
      IDIF2=2  
      ENDIF
      RETURN
      END  
C*****************************     LEGEND
          SUBROUTINE LEGEND(ICASO,ILI) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          IF(ILI.EQ.1) THEN
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
          WRITE(6,1035)
          WRITE(6,1045)
          WRITE(6,1100)
          WRITE(6,1110)
          WRITE(6,1120)
          WRITE(6,1130)
          ENDIF
          IF(ICASO.EQ.3) THEN  
          WRITE(6,1010)
          WRITE(6,1020)
          WRITE(6,1038)
          WRITE(6,1048)
          WRITE(6,1050)
          WRITE(6,1060)
          WRITE(6,1100)
          WRITE(6,1110)
          ENDIF
          ELSE
          WRITE(6,2000)
          IF(ICASO.EQ.1) THEN  
          WRITE(6,2010)
          WRITE(6,2020)
          WRITE(6,2030)
          WRITE(6,2040)
          WRITE(6,2050)
          WRITE(6,2060)
          WRITE(6,2070)
          WRITE(6,2080)
          ENDIF
          IF(ICASO.EQ.2) THEN  
          WRITE(6,2010)
          WRITE(6,2020)
          WRITE(6,2035)
          WRITE(6,2045)
          WRITE(6,2100)
          WRITE(6,2110)
          WRITE(6,2120)
          WRITE(6,2130)
          ENDIF
          IF(ICASO.EQ.3) THEN  
          WRITE(6,2010)
          WRITE(6,2020)
          WRITE(6,2038)
          WRITE(6,2048)
          WRITE(6,2050)
          WRITE(6,2060)
          WRITE(6,2100)
          WRITE(6,2110)
          ENDIF
          ENDIF
          RETURN
1000      FORMAT(///5X,'CASO:',15X,'FUNZIONE TEORICA:'/)
1010      FORMAT(9X,'1',15X,'CASO ESPONENZIALE'
     ./25X,'A*EXP(-B*ABS(X))'/)
1020      FORMAT(9X,'2',15X,'CASO NORMALE' 
     ./25X,'A*EXP(-B*X**2)'/)  
1030      FORMAT(9X,'3',15X,'CASO ESPONENZIALE PARABOLA'
     ./25X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(2*D)'/)
1035      FORMAT(9X,'3',15X,'CASO ESPONENZIALE PARABOLA'
     ./25X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(6*D)'/)
1038     FORMAT(9X,'3',15X,'CASO ESPONENZIALE PARABOLA'
     ./25X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(24*D)'/)
1040      FORMAT(9X,'4',15X,'CASO NORMALE PARABOLA'
     ./25X,'A*(1-D*X**2)*EXP(-B*X**2)     B.GT.D/2'/)
1045      FORMAT(9X,'4',15X,'CASO NORMALE PARABOLA'
     ./25X,'A*(1-D*X**2)*EXP(-B*X**2)     B.GT.D'/)
1048      FORMAT(9X,'4',15X,'CASO NORMALE PARABOLA'
     ./25X,'A*(1-D*X**2)*EXP(-B*X**2)     B.GT.2*D'/)
1050      FORMAT(9X,'5',15X,'CASO ESPONENZIALE SENO SU ASCISSA'
     ./25X,'A*EXP(-B*ABS(X))*SIN(C*ABS(X))/(C*ABS(X))'/)
1060      FORMAT(9X,'6',15X,'CASO NORMALE SENO SU ASCISSA' 
     ./25X,'A*EXP(-B*X**2)*SIN(C*ABS(X))/(C*ABS(X))'/) 
1070      FORMAT(9X,'7',15X,'CASO ESPONENZIALE COSENO' 
     ./25X,'A*EXP(-B*ABS(X))*COS(C*ABS(X))'/)  
1080      FORMAT(9X,'8',15X,'CASO NORMALE COSENO'  
     ./25X,'A*EXP(-B*X**2)*COS(C*ABS(X))'/)
1100      FORMAT(9X,'9',15X,'CASO ESPONENZIALE J1 SU ASCISSA'  
     ./25X,'2*A*EXP(-B*ABS(X))*J1(C*ABS(X))/(C*ABS(X))'/)  
1110      FORMAT(8X,'10',15X,'CASO NORMALE J1 SU ASCISSA'  
     ./25X,'2*A*EXP(-B*X**2)*J1(C*ABS(X))/(C*ABS(X))'/)
1120      FORMAT(8X,'11',15X,'CASO ESPONENZIALE J0'
     ./25X,'A*EXP(-B*ABS(X))*J0(C*ABS(X))'/)
1130      FORMAT(8X,'12',15X,'CASO NORMALE J0' 
     ./25X,'A*EXP(-B*X**2)*J0(C*ABS(X))'/) 
2000      FORMAT(///5X,'CASE:',15X,'THEORETICAL FUNCTION:'/)
2010      FORMAT(9X,'1',15X,'EXPONENTIAL CASE'
     ./25X,'A*EXP(-B*ABS(X))'/)
2020      FORMAT(9X,'2',15X,'NORMAL CASE' 
     ./25X,'A*EXP(-B*X**2)'/)  
2030      FORMAT(9X,'3',15X,'EXPONENTIAL PARABOLA CASE'
     ./25X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(2*D)'/)
2035      FORMAT(9X,'3',15X,'EXPONENTIAL PARABOLA CASE'
     ./25X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(6*D)'/)
2038     FORMAT(9X,'3',15X,'EXPONENTIAL PARABOLA CASE'
     ./25X,'A*(1-D*X**2)*EXP(-B*ABS(X))   B.GT.SQRT(24*D)'/)
2040      FORMAT(9X,'4',15X,'NORMAL PARABOLA CASE'
     ./25X,'A*(1-D*X**2)*EXP(-B*X**2)     B.GT.D/2'/)
2045      FORMAT(9X,'4',15X,'NORMAL PARABOLA CASE'
     ./25X,'A*(1-D*X**2)*EXP(-B*X**2)     B.GT.D'/)
2048      FORMAT(9X,'4',15X,'NORMAL PARABOLA CASE'
     ./25X,'A*(1-D*X**2)*EXP(-B*X**2)     B.GT.2*D'/)
2050      FORMAT(9X,'5',15X,'EXPONENTIAL SINE OVER X-COORDINATE CASE'
     ./25X,'A*EXP(-B*ABS(X))*SIN(C*ABS(X))/(C*ABS(X))'/)
2060      FORMAT(9X,'6',15X,'NORMAL SINE OVER X-COORDINATE CASE' 
     ./25X,'A*EXP(-B*X**2)*SIN(C*ABS(X))/(C*ABS(X))'/) 
2070      FORMAT(9X,'7',15X,'EXPONENTIAL COSINE CASE' 
     ./25X,'A*EXP(-B*ABS(X))*COS(C*ABS(X))'/)  
2080      FORMAT(9X,'8',15X,'NORMAL COSINE CASE'  
     ./25X,'A*EXP(-B*X**2)*COS(C*ABS(X))'/)
2100      FORMAT(9X,'9',15X,'EXPONENTIAL J1 OVER X-COORDINATE CASE'  
     ./25X,'2*A*EXP(-B*ABS(X))*J1(C*ABS(X))/(C*ABS(X))'/)  
2110      FORMAT(8X,'10',15X,'NORMAL J1 OVER X-COORDINATE CASE'  
     ./25X,'2*A*EXP(-B*X**2)*J1(C*ABS(X))/(C*ABS(X))'/)
2120      FORMAT(8X,'11',15X,'EXPONENTIAL J0 CASE'
     ./25X,'A*EXP(-B*ABS(X))*J0(C*ABS(X))'/)
2130      FORMAT(8X,'12',15X,'NORMAL J0 CASE' 
     ./25X,'A*EXP(-B*X**2)*J0(C*ABS(X))'/) 
          END  
C*****************************     LINGUA
          SUBROUTINE LINGUA(ILI)
10        WRITE(*,100)
          READ(*,*,ERR=10) ILI
          WRITE(*,'(///)')
          RETURN
100       FORMAT(///5X,'ITALIANO (ITALIAN): 1'/
     .    5X,'(INGLESE) ENGLISH:  0')
          END  
C*****************************     MAT3
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
C*****************************     MEDVAR
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
C*****************************     MEDVA2
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
C*****************************     MIMA
          SUBROUTINE MIMA(COV,IMM,IX0,NCOV)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION COV(*),IX0(*),IMM(*)
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
C*****************************     OTTIP
          SUBROUTINE OTTIP1(X,Y,W,IPXW,IPYW,
     .    LPXW,LPYW,GMMW,SRV,NSRV,AINF,ASUP,STP,IWD, 
     .    NMAX,NWDMAX,SCX,SCY,SCW,IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,
     .    ILI)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA,WFORMA,WFORM2 
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),W(*),IPXW(*),IPYW(*),
     .    LPXW(*),LPYW(*),GMMW(*),SRV(*),NSRV(*), 
     .    IWD(*),IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20),WFORMA(30),WFORM2(30) 
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          DATA WFORMA/'(//6','X,''M','(X)''',',9X,','''SQM','(X)''',
     .    ',8X,','''RO-','1(X)',''',7X',',''SQ','M-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA WFORM2/'(//6','X,''M','(X)''',',10X,','''SD','(X)''',
     .    ',8X,','''RO-','1(X)',''',8X',',''SD','-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NW=0 
          ND=0
          XMIN= 1.D30
          XMAX=-1.D30
          YMIN= 1.D30
          YMAX=-1.D30 
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
          READ(CDATO,FORMA) NPT,XS,YS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY 
          W(ND)=WS*SCW 
          XMIN=DMIN1(XMIN,X(ND))
          XMAX=DMAX1(XMAX,X(ND))
          YMIN=DMIN1(YMIN,Y(ND))
          YMAX=DMAX1(YMAX,Y(ND))
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NW=NW+1  
          IPXW(NW)=ND  
          IPYW(NW)=ND  
          ENDIF
          GOTO 10  
C
20        CALL MEDVA2(W,IPXW,EW,S2W,SW,NW) 
          DO 30 I=1,ND 
          W(I)=W(I)-EW 
30        CONTINUE 
          CALL RBORDA(X,Y,IPXW,IPYW,LPXW,LPYW,NW,ND)
          VOL=(XMAX-XMIN)*(YMAX-YMIN)
          UNI=DSQRT(VOL/ND)
          IF(IAUTO.EQ.1) THEN
          AINF=UNI/10.D0
          ASUP=5.D0*UNI
          STP=UNI/10.D0
          ENDIF
          NSTP=IDINT((ASUP-AINF)/STP)+1
          RAGGIO=AINF
          ALATO=2.D0*ASUP
          CALL COROP(X,Y,W,IPXW,IPYW,LPXW,LPYW,NW,NSTP,RAGGIO,ALATO,  
     .    GMMW,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          OLATO=0.D0
          S2NSE=1.D30 
          IF(ILI.EQ.1) THEN
          WRITE(6,1200)
          ELSE
          WRITE(6,1201)
          ENDIF
          DO 60 ISTP=1,NSTP
          SIGMA2=(S2W-GMMW(ISTP)) 
          WRITE(6,1300) ISTP*RAGGIO,DSQRT(SIGMA2)
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=2.D0*ISTP*RAGGIO  
          S2NSE=SIGMA2 
          SNW=DSQRT(S2W-GMMW(ISTP))
          GMW=GMMW(ISTP) 
          ENDIF
60        CONTINUE 
C
          NWM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMP(X,Y,W,IPXW,IPYW,LPXW,LPYW,IWD,SNW,OLATO,  
     .    1,NW,ND,NWM,NWDMAX,NMAX,IAST,3,4,TOLL,ILI)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1100) OLATO/2.D0
          ELSE 
          WRITE(6,1101) OLATO/2.D0
          ENDIF 
          IF(NWM.LT.9) THEN
          WFORMA(20)='I1,/'
          WFORMA(26)='I8) '
          WFORM2(20)='I1,/'
          WFORM2(26)='I8) '
          ENDIF
          IF(NWM.GE.10.AND.NWM.LT.100) THEN
          WFORMA(20)='I2,/'
          WFORMA(26)='I9) '
          WFORM2(20)='I2,/'
          WFORM2(26)='I9) '
          ENDIF
          IF(NWM.GE.100) THEN  
          WFORMA(20)='I3,/'
          WFORMA(26)='I10)'
          WFORM2(20)='I3,/'
          WFORM2(26)='I10)'
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,WFORMA) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM 
          ELSE 
          WRITE(6,WFORM2) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM 
          ENDIF 
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'RAGGIO OTTIMALE = ',F10.3,//)  
1101      FORMAT(///5X,'OPTIMAL RADIUS = ',F10.3,//)  
1200      FORMAT(///14X,'RAGGIO',7X,'SQM-N(X)',//) 
1201      FORMAT(///14X,'RADIUS',8X,'SD-N(X)',//) 
1300      FORMAT(5X,2F15.3)
          END  
C   
          SUBROUTINE OTTIP2(X,Y,U,V,IPXU,IPYU,IPXV,IPYV,
     .    LPXU,LPYU,LPXV,LPYV,GMMU,GMMV,SRV,NSRV,AINF,ASUP,STP,IWD,
     .    NMAX,NWDMAX,SCX,SCY,SCU,SCV,IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,
     .    ILI)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA,UFORMA,VFORMA,UFORM2,VFORM2  
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),U(*),V(*),IPXU(*),IPYU(*),
     .    IPXV(*),IPYV(*),LPXU(*),LPYU(*),LPXV(*),LPYV(*), 
     .    GMMU(*),GMMV(*),SRV(*),NSRV(*),IWD(*),IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20),UFORMA(30),VFORMA(30),
     .    UFORM2(30),VFORM2(30)
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
          DATA UFORM2/'(//6','X,''M','(X)''',',10X,','''SD','(X)''',
     .    ',8X,','''RO-','1(X)',''',8X',',''SD','-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA VFORM2/'(//6','X,''M','(Y)''',',10X,','''SD','(Y)''',
     .    ',8X,','''RO-','1(Y)',''',8X',',''SD','-N(','Y)'',',
     .    '11X,','''N(Y',')''  ',',7X,', 
     .    '''N(Y',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NU=0 
          NV=0 
          ND=0 
          XMIN= 1.D30
          XMAX=-1.D30
          YMIN= 1.D30
          YMAX=-1.D30 
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
          READ(CDATO,FORMA) NPT,XS,YS,US,VS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY 
          U(ND)=US*SCU 
          V(ND)=VS*SCV
          XMIN=DMIN1(XMIN,X(ND))
          XMAX=DMAX1(XMAX,X(ND))
          YMIN=DMIN1(YMIN,Y(ND))
          YMAX=DMAX1(YMAX,Y(ND))           
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
C
20        CALL MEDVA2(U,IPXU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPXV,EV,S2V,SV,NV) 
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
30        CONTINUE 
          CALL RBORDA(X,Y,IPXU,IPYU,LPXU,LPYU,NU,ND)
          CALL RBORDA(X,Y,IPXV,IPYV,LPXV,LPYV,NV,ND)
          VOL=(XMAX-XMIN)*(YMAX-YMIN)
          UNI=DSQRT(VOL/ND)
          IF(IAUTO.EQ.1) THEN
          AINF=UNI/10.D0
          ASUP=5.D0*UNI
          STP=UNI/10.D0
          ENDIF          
          NSTP=IDINT((ASUP-AINF)/STP)+1
          RAGGIO=AINF 
          ALATO=2.D0*ASUP
          CALL COROP(X,Y,U,IPXU,IPYU,LPXU,LPYU,NU,NSTP,RAGGIO,ALATO,  
     .    GMMU,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          CALL COROP(X,Y,V,IPXV,IPYV,LPXV,LPYV,NV,NSTP,RAGGIO,ALATO,  
     .    GMMV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          OLATO=0.D0
          S2NSE=1.D30 
          IF(ILI.EQ.1) THEN
          WRITE(6,1200)
          ELSE
          WRITE(6,1201)
          ENDIF
          DO 60 ISTP=1,NSTP
          SIGMA2=(S2U-GMMU(ISTP))+(S2V-GMMV(ISTP))
          SIGMA2=SIGMA2/2.D0
          WRITE(6,1300) ISTP*RAGGIO,DSQRT(S2U-GMMU(ISTP)),  
     .    DSQRT(S2V-GMMV(ISTP)),DSQRT(SIGMA2)
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=2.D0*ISTP*RAGGIO  
          S2NSE=SIGMA2 
          SNU=DSQRT(S2U-GMMU(ISTP))
          SNV=DSQRT(S2V-GMMV(ISTP))
          GMU=GMMU(ISTP)
          GMV=GMMV(ISTP)
          ENDIF
60        CONTINUE 
C
          NUM=0
          NVM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMP(X,Y,U,IPXU,IPYU,LPXU,LPYU,IWD,SNU,OLATO,  
     .    1,NU,ND,NUM,NWDMAX,NMAX,IAST,4,3,TOLL,ILI)
          CALL TRIMP(X,Y,V,IPXV,IPYV,LPXV,LPYV,IWD,SNV,OLATO,  
     .    2,NV,ND,NVM,NWDMAX,NMAX,IAST,3,4,TOLL,ILI)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1100) OLATO/2.D0
          ELSE 
          WRITE(6,1101) OLATO/2.D0
          ENDIF 
          IF(NUM.LT.9) THEN
          UFORMA(20)='I1,/'
          UFORMA(26)='I8) '
          UFORM2(20)='I1,/'
          UFORM2(26)='I8) '
          ENDIF
          IF(NUM.GE.10.AND.NUM.LT.100) THEN
          UFORMA(20)='I2,/'
          UFORMA(26)='I9) '
          UFORM2(20)='I2,/'
          UFORM2(26)='I9) '
          ENDIF
          IF(NUM.GE.100) THEN  
          UFORMA(20)='I3,/'
          UFORMA(26)='I10)'
          UFORM2(20)='I3,/'
          UFORM2(26)='I10)'
          ENDIF
          IF(NVM.LT.9) THEN
          VFORMA(20)='I1,/'
          VFORMA(26)='I8) '
          VFORM2(20)='I1,/'
          VFORM2(26)='I8) '
          ENDIF
          IF(NVM.GE.10.AND.NVM.LT.100) THEN
          VFORMA(20)='I2,/'
          VFORMA(26)='I9) '
          VFORM2(20)='I2,/'
          VFORM2(26)='I9) '
          ENDIF
          IF(NVM.GE.100) THEN  
          VFORMA(20)='I3,/'
          VFORMA(26)='I10)'
          VFORM2(20)='I3,/'
          VFORM2(26)='I10)'
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,UFORMA) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORMA) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM
          ELSE  
          WRITE(6,UFORM2) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORM2) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM 
          ENDIF 
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'RAGGIO OTTIMALE = ',F10.3,//)  
1101      FORMAT(///5X,'OPTIMAL RADIUS = ',F10.3,//)  
1200      FORMAT(///14X,'RAGGIO',7X,'SQM-N(X)',7X,'SQM-N(Y)',5X,
     .    'SQM-N(XUY)',//) 
1201      FORMAT(///14X,'RADIUS',8X,'SD-N(X)',8X,'SD-N(Y)',6X,
     .    'SD-N(XUY)',//) 
1300      FORMAT(5X,4F15.3)
          END  
C
          SUBROUTINE OTTIP3(X,Y,U,V,W,IPXU,IPYU,IPXV,IPYV,IPXW,IPYW,
     .    LPXU,LPYU,LPXV,LPYV,LPXW,LPYW,
     .    GMMU,GMMV,GMMW,SRV,NSRV,AINF,ASUP,STP,IWD, 
     .    NMAX,NWDMAX,SCX,SCY,SCU,SCV,SCW,IPTO,IAST,FORMA,IAUTO,ITRIM,
     .    TOLL,ILI)       
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA,UFORMA,VFORMA,WFORMA,UFORM2,VFORM2,WFORM2
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),U(*),V(*),W(*),IPXU(*),IPYU(*),  
     .    IPXV(*),IPYV(*),IPXW(*),IPYW(*),LPXU(*),LPYU(*), 
     .    LPXV(*),LPYV(*),LPXW(*),LPYW(*),GMMU(*),GMMV(*),GMMW(*),
     .    SRV(*),NSRV(*),IWD(*),IPTO(*),IAST(*)
          DIMENSION DATO(20),FORMA(20),UFORMA(30),VFORMA(30),
     .    WFORMA(30),UFORM2(30),VFORM2(30),WFORM2(30)
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
          DATA UFORM2/'(//6','X,''M','(X)''',',10X,','''SD','(X)''',
     .    ',8X,','''RO-','1(X)',''',8X',',''SD','-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA VFORM2/'(//6','X,''M','(Y)''',',10X,','''SD','(Y)''',
     .    ',8X,','''RO-','1(Y)',''',8X',',''SD','-N(','Y)'',',
     .    '11X,','''N(Y',')''  ',',7X,', 
     .    '''N(Y',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA WFORM2/'(//6','X,''M','(Z)''',',10X,','''SD','(Z)''',
     .    ',8X,','''RO-','1(Z)',''',8X',',''SD','-N(','Z)'',',
     .    '11X,','''N(Z',')''  ',',7X,', 
     .    '''N(Z',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NU=0 
          NV=0 
          NW=0 
          ND=0
          XMIN= 1.D30
          XMAX=-1.D30
          YMIN= 1.D30
          YMAX=-1.D30            
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
          READ(CDATO,FORMA) NPT,XS,YS,US,VS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          W(ND)=WS*SCW
          XMIN=DMIN1(XMIN,X(ND))
          XMAX=DMAX1(XMAX,X(ND))
          YMIN=DMIN1(YMIN,Y(ND))
          YMAX=DMAX1(YMAX,Y(ND))           
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
C
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
          VOL=(XMAX-XMIN)*(YMAX-YMIN)
          UNI=DSQRT(VOL/ND)
          IF(IAUTO.EQ.1) THEN
          AINF=UNI/10.D0
          ASUP=5.D0*UNI
          STP=UNI/10.D0
          ENDIF          
          NSTP=IDINT((ASUP-AINF)/STP)+1
          RAGGIO=AINF 
          ALATO=2.D0*ASUP
          CALL COROP(X,Y,U,IPXU,IPYU,LPXU,LPYU,NU,NSTP,RAGGIO,ALATO,  
     .    GMMU,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          CALL COROP(X,Y,V,IPXV,IPYV,LPXV,LPYV,NV,NSTP,RAGGIO,ALATO,  
     .    GMMV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          CALL COROP(X,Y,W,IPXW,IPYW,LPXW,LPYW,NW,NSTP,RAGGIO,ALATO,  
     .    GMMW,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          OLATO=0.D0
          S2NSE=1.D30 
          IF(ILI.EQ.1) THEN
          WRITE(6,1200)
          ELSE
          WRITE(6,1201)
          ENDIF
          DO 60 ISTP=1,NSTP
          SIGMA2=(S2U-GMMU(ISTP))+(S2V-GMMV(ISTP))+(S2W-GMMW(ISTP)) 
          SIGMA2=SIGMA2/3.D0
          WRITE(6,1300) ISTP*RAGGIO,DSQRT(S2U-GMMU(ISTP)),  
     .    DSQRT(S2V-GMMV(ISTP)),DSQRT(S2W-GMMW(ISTP)),DSQRT(SIGMA2) 
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=2.D0*ISTP*RAGGIO  
          S2NSE=SIGMA2 
          SNU=DSQRT(S2U-GMMU(ISTP))
          SNV=DSQRT(S2V-GMMV(ISTP))
          SNW=DSQRT(S2W-GMMW(ISTP))
          GMU=GMMU(ISTP)
          GMV=GMMV(ISTP)
          GMW=GMMW(ISTP)
          ENDIF
60        CONTINUE 
C
          NUM=0
          NVM=0
          NWM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMP(X,Y,U,IPXU,IPYU,LPXU,LPYU,IWD,SNU,OLATO,  
     .    1,NU,ND,NUM,NWDMAX,NMAX,IAST,3,4,TOLL,ILI)
          CALL TRIMP(X,Y,V,IPXV,IPYV,LPXV,LPYV,IWD,SNV,OLATO,  
     .    2,NV,ND,NVM,NWDMAX,NMAX,IAST,4,3,TOLL,ILI)
          CALL TRIMP(X,Y,W,IPXW,IPYW,LPXW,LPYW,IWD,SNW,OLATO,  
     .    3,NW,ND,NWM,NWDMAX,NMAX,IAST,3,4,TOLL,ILI)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1100) OLATO/2.D0
          ELSE 
          WRITE(6,1101) OLATO/2.D0
          ENDIF 
          IF(NUM.LT.9) THEN
          UFORMA(20)='I1,/'
          UFORMA(26)='I8) '
          UFORM2(20)='I1,/'
          UFORM2(26)='I8) '
          ENDIF
          IF(NUM.GE.10.AND.NUM.LT.100) THEN
          UFORMA(20)='I2,/'
          UFORMA(26)='I9) '
          UFORM2(20)='I2,/'
          UFORM2(26)='I9) '
          ENDIF
          IF(NUM.GE.100) THEN  
          UFORMA(20)='I3,/'
          UFORMA(26)='I10)'
          UFORM2(20)='I3,/'
          UFORM2(26)='I10)'
          ENDIF
          IF(NVM.LT.9) THEN
          VFORMA(20)='I1,/'
          VFORMA(26)='I8) '
          VFORM2(20)='I1,/'
          VFORM2(26)='I8) '
          ENDIF
          IF(NVM.GE.10.AND.NVM.LT.100) THEN
          VFORMA(20)='I2,/'
          VFORMA(26)='I9) '
          VFORM2(20)='I2,/'
          VFORM2(26)='I9) '
          ENDIF
          IF(NVM.GE.100) THEN  
          VFORMA(20)='I3,/'
          VFORMA(26)='I10)'
          VFORM2(20)='I3,/'
          VFORM2(26)='I10)'
          ENDIF
          IF(NWM.LT.9) THEN
          WFORMA(20)='I1,/'
          WFORMA(26)='I8) '
          WFORM2(20)='I1,/'
          WFORM2(26)='I8) '
          ENDIF
          IF(NWM.GE.10.AND.NWM.LT.100) THEN
          WFORMA(20)='I2,/'
          WFORMA(26)='I9) '
          WFORM2(20)='I2,/'
          WFORM2(26)='I9) '
          ENDIF
          IF(NWM.GE.100) THEN  
          WFORMA(20)='I3,/'
          WFORMA(26)='I10)'
          WFORM2(20)='I3,/'
          WFORM2(26)='I10)'
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,UFORMA) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORMA) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM  
          WRITE(6,WFORMA) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM
          ELSE  
          WRITE(6,UFORM2) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORM2) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM  
          WRITE(6,WFORM2) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM
          ENDIF  
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'RAGGIO OTTIMALE = ',F10.3,//)  
1101      FORMAT(///5X,'OPTIMAL RADIUS = ',F10.3,//)  
1200      FORMAT(///14X,'RAGGIO',7X,'SQM-N(X)',7X,'SQM-N(Y)',7X,
     .    'SQM-N(Z)',3X,'SQM-N(XUYUZ)'//) 
1201      FORMAT(///14X,'RADIUS',8X,'SD-N(X)',8X,'SD-N(Y)',8X,
     .    'SD-N(Z)',4X,'SD-N(XUYUZ)'//)  
1300      FORMAT(5X,5F15.3)
          END 
C*****************************     OTTIS
          SUBROUTINE OTTIS1(X,Y,Z,W,IPXW,IPYW,IPZW,
     .    LPXW,LPYW,LPZW,GMMW,SRV,NSRV,AINF,ASUP,STP,IWD, 
     .    NMAX,NWDMAX,SCX,SCY,SCZ,SCW,IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,
     .    ILI)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA,WFORMA,WFORM2 
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),Z(*),W(*),IPXW(*),IPYW(*),IPZW(*),
     .    LPXW(*),LPYW(*),LPZW(*),GMMW(*),SRV(*),NSRV(*),
     .    IWD(*),IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20),WFORMA(30),WFORM2(30) 
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          DATA WFORMA/'(//6','X,''M','(X)''',',9X,','''SQM','(X)''',
     .    ',8X,','''RO-','1(X)',''',7X',',''SQ','M-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA WFORM2/'(//6','X,''M','(X)''',',10X,','''SD','(X)''',
     .    ',8X,','''RO-','1(X)',''',8X',',''SD','-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NW=0 
          ND=0 
          XMIN= 1.D30
          XMAX=-1.D30
          YMIN= 1.D30
          YMAX=-1.D30   
          ZMIN= 1.D30
          ZMAX=-1.D30          
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120) 
          READ(CDATO,FORMA) NPT,XS,YS,ZS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY
          Z(ND)=ZS*SCZ 
          W(ND)=WS*SCW
          XMIN=DMIN1(XMIN,X(ND))
          XMAX=DMAX1(XMAX,X(ND))
          YMIN=DMIN1(YMIN,Y(ND))
          YMAX=DMAX1(YMAX,Y(ND))
          ZMIN=DMIN1(ZMIN,Z(ND))
          ZMAX=DMAX1(ZMAX,Z(ND))           
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NW=NW+1  
          IPXW(NW)=ND  
          IPYW(NW)=ND
          IPZW(NW)=ND  
          ENDIF
          GOTO 10  
C
20        CALL MEDVA2(W,IPXW,EW,S2W,SW,NW) 
          DO 30 I=1,ND 
          W(I)=W(I)-EW 
30        CONTINUE 
          CALL RTORDA(X,Y,Z,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,NW,ND)
          VOL=(XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN)
          UNI=DEXP(DLOG(VOL/ND)/3.D0)
          IF(IAUTO.EQ.1) THEN
          AINF=UNI/10.D0
          ASUP=5.D0*UNI
          STP=UNI/10.D0
          ENDIF         
          NSTP=IDINT((ASUP-AINF)/STP)+1
          RAGGIO=AINF 
          ALATO=2.D0*ASUP
          CALL COROS(X,Y,Z,W,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,NW,NSTP,
     .    RAGGIO,ALATO,GMMW,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          OLATO=0.D0
          S2NSE=1.D30
          IF(ILI.EQ.1) THEN 
          WRITE(6,1200)
          ELSE
          WRITE(6,1201)
          ENDIF
          DO 60 ISTP=1,NSTP
          SIGMA2=(S2W-GMMW(ISTP)) 
          WRITE(6,1300) ISTP*RAGGIO,DSQRT(SIGMA2)
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=2.D0*ISTP*RAGGIO  
          S2NSE=SIGMA2 
          SNW=DSQRT(S2W-GMMW(ISTP))
          GMW=GMMW(ISTP) 
          ENDIF
60        CONTINUE 
C
          NWM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMS(X,Y,Z,W,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,
     .    IWD,SNW,OLATO,1,NW,ND,NWM,NWDMAX,NMAX,IAST,3,4,TOLL,ILI)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1100) OLATO/2.D0 
          ELSE
          WRITE(6,1101) OLATO/2.D0
          ENDIF 
          IF(NWM.LT.9) THEN
          WFORMA(20)='I1,/'
          WFORMA(26)='I8) '
          WFORM2(20)='I1,/'
          WFORM2(26)='I8) '
          ENDIF
          IF(NWM.GE.10.AND.NWM.LT.100) THEN
          WFORMA(20)='I2,/'
          WFORMA(26)='I9) '
          WFORM2(20)='I2,/'
          WFORM2(26)='I9) '
          ENDIF
          IF(NWM.GE.100) THEN  
          WFORM2(20)='I3,/'
          WFORM2(26)='I10)'
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,WFORMA) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM 
          ELSE 
          WRITE(6,WFORM2) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM
          ENDIF  
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'RAGGIO OTTIMALE = ',F10.3,//)  
1101      FORMAT(///5X,'OPTIMAL RADIUS = ',F10.3,//)  
1200      FORMAT(///14X,'RAGGIO',7X,'SQM-N(X)',//) 
1201      FORMAT(///14X,'RADIUS',8X,'SD-N(X)',//) 
1300      FORMAT(5X,2F15.3)
          END  
C
          SUBROUTINE OTTIS2(X,Y,Z,U,V,IPXU,IPYU,IPZU,
     .    IPXV,IPYV,IPZV,LPXU,LPYU,LPZU,
     .    LPXV,LPYV,LPZV,GMMU,GMMV,SRV,NSRV,AINF,ASUP,STP,IWD,
     .    NMAX,NWDMAX,SCX,SCY,SCZ,SCU,SCV,
     .    IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,ILI)       
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA,UFORMA,VFORMA,UFORM2,VFORM2  
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),Z(*),U(*),V(*),IPXU(*),IPYU(*),
     .    IPZU(*),IPXV(*),IPYV(*),IPZV(*),LPXU(*),LPYU(*), 
     .    LPZU(*),LPXV(*),LPYV(*),LPZV(*),
     .    GMMU(*),GMMV(*),SRV(*),NSRV(*),IWD(*),IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20),UFORMA(30),VFORMA(30),
     .    UFORM2(30),VFORM2(30) 
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
          DATA UFORM2/'(//6','X,''M','(X)''',',10X,','''SD','(X)''',
     .    ',8X,','''RO-','1(X)',''',8X',',''SD','-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA VFORM2/'(//6','X,''M','(Y)''',',10X,','''SD','(Y)''',
     .    ',8X,','''RO-','1(Y)',''',8X',',''SD','-N(','Y)'',',
     .    '11X,','''N(Y',')''  ',',7X,', 
     .    '''N(Y',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NU=0 
          NV=0 
          ND=0 
          XMIN= 1.D30
          XMAX=-1.D30
          YMIN= 1.D30
          YMAX=-1.D30   
          ZMIN= 1.D30
          ZMAX=-1.D30  
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
          READ(CDATO,FORMA) NPT,XS,YS,ZS,US,VS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY
          Z(ND)=ZS*SCZ 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          XMIN=DMIN1(XMIN,X(ND))
          XMAX=DMAX1(XMAX,X(ND))
          YMIN=DMIN1(YMIN,Y(ND))
          YMAX=DMAX1(YMAX,Y(ND))
          ZMIN=DMIN1(ZMIN,Z(ND))
          ZMAX=DMAX1(ZMAX,Z(ND))            
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPXU(NU)=ND  
          IPYU(NU)=ND
          IPZU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPXV(NV)=ND  
          IPYV(NV)=ND
          IPZV(NV)=ND  
          ENDIF
          GOTO 10  
C
20        CALL MEDVA2(U,IPXU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPXV,EV,S2V,SV,NV) 
          DO 30 I=1,ND 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
30        CONTINUE 
          CALL RTORDA(X,Y,Z,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,NU,ND)
          CALL RTORDA(X,Y,Z,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,NV,ND)
          VOL=(XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN)
          UNI=DEXP(DLOG(VOL/ND)/3.D0)
          IF(IAUTO.EQ.1) THEN
          AINF=UNI/10.D0
          ASUP=5.D0*UNI
          STP=UNI/10.D0
          ENDIF            
          NSTP=IDINT((ASUP-AINF)/STP)+1
          RAGGIO=AINF 
          ALATO=2.D0*ASUP
          CALL COROS(X,Y,Z,U,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,NU,NSTP,
     .    RAGGIO,ALATO,GMMU,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          CALL COROS(X,Y,Z,V,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,NV,NSTP,
     .    RAGGIO,ALATO,GMMV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          OLATO=0.D0
          S2NSE=1.D30
          IF(ILI.EQ.1) THEN 
          WRITE(6,1200)
          ELSE
          WRITE(6,1201)
          ENDIF
          DO 60 ISTP=1,NSTP
          SIGMA2=(S2U-GMMU(ISTP))+(S2V-GMMV(ISTP))
          SIGMA2=SIGMA2/2.D0
          WRITE(6,1300) ISTP*RAGGIO,DSQRT(S2U-GMMU(ISTP)),  
     .    DSQRT(S2V-GMMV(ISTP)),DSQRT(SIGMA2)
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=2.D0*ISTP*RAGGIO  
          S2NSE=SIGMA2 
          SNU=DSQRT(S2U-GMMU(ISTP))
          SNV=DSQRT(S2V-GMMV(ISTP))
          GMU=GMMU(ISTP)
          GMV=GMMV(ISTP)
          ENDIF
60        CONTINUE 
C
          NUM=0
          NVM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMS(X,Y,Z,U,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,
     .    IWD,SNU,OLATO,1,NU,ND,NUM,NWDMAX,NMAX,IAST,4,3,TOLL,ILI)
          CALL TRIMS(X,Y,Z,V,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,
     .    IWD,SNV,OLATO,2,NV,ND,NVM,NWDMAX,NMAX,IAST,3,4,TOLL,ILI)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1100) OLATO/2.D0 
          ELSE
          WRITE(6,1101) OLATO/2.D0 
          ENDIF
          IF(NUM.LT.9) THEN
          UFORMA(20)='I1,/'
          UFORMA(26)='I8) '
          UFORM2(20)='I1,/'
          UFORM2(26)='I8) '
          ENDIF
          IF(NUM.GE.10.AND.NUM.LT.100) THEN
          UFORMA(20)='I2,/'
          UFORMA(26)='I9) '
          UFORM2(20)='I2,/'
          UFORM2(26)='I9) '
          ENDIF
          IF(NUM.GE.100) THEN  
          UFORMA(20)='I3,/'
          UFORMA(26)='I10)'
          UFORM2(20)='I3,/'
          UFORM2(26)='I10)'
          ENDIF
          IF(NVM.LT.9) THEN
          VFORMA(20)='I1,/'
          VFORMA(26)='I8) '
          VFORM2(20)='I1,/'
          VFORM2(26)='I8) '
          ENDIF
          IF(NVM.GE.10.AND.NVM.LT.100) THEN
          VFORMA(20)='I2,/'
          VFORMA(26)='I9) '
          VFORM2(20)='I2,/'
          VFORM2(26)='I9) '
          ENDIF
          IF(NVM.GE.100) THEN  
          VFORMA(20)='I3,/'
          VFORMA(26)='I10)'
          VFORM2(20)='I3,/'
          VFORM2(26)='I10)'
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,UFORMA) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORMA) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM
          ELSE  
          WRITE(6,UFORM2) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORM2) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM
          ENDIF  
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'RAGGIO OTTIMALE = ',F10.3,//)  
1101      FORMAT(///5X,'OPTIMAL RADIUS = ',F10.3,//)  
1200      FORMAT(///14X,'RAGGIO',7X,'SQM-N(X)',7X,'SQM-N(Y)',5X,
     .    'SQM-N(XUY)',//) 
1201      FORMAT(///14X,'RADIUS',8X,'SD-N(X)'8X,'SD-N(Y)',6X,
     .    'SD-N(XUY)',//) 
1300      FORMAT(5X,4F15.3)
          END  
C
          SUBROUTINE OTTIS3(X,Y,Z,U,V,W,IPXU,IPYU,IPZU,IPXV,IPZV,IPYV,
     .    IPXW,IPYW,IPZW,LPXU,LPYU,LPZU,LPXV,LPYV,LPZV,LPXW,LPYW,LPZW,
     .    GMMU,GMMV,GMMW,SRV,NSRV,AINF,ASUP,STP,IWD,NMAX,NWDMAX,
     .    SCX,SCY,SCZ,SCU,SCV,SCW,IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,
     .    ILI)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA,UFORMA,VFORMA,WFORMA,UFORM2,VFORM2,WFORM2
          CHARACTER*120 CDATO
          DIMENSION X(*),Y(*),Z(*),U(*),V(*),W(*),IPXU(*),IPYU(*),  
     .    IPZU(*),IPXV(*),IPYV(*),IPZV(*),IPXW(*),IPYW(*),IPZW(*),
     .    LPXU(*),LPYU(*),LPXV(*),LPYV(*),LPZV(*),LPXW(*),LPYW(*),
     .    LPZW(*),GMMU(*),GMMV(*),GMMW(*),SRV(*),NSRV(*),
     .    IWD(*),IPTO(*),IAST(*)
          DIMENSION DATO(120),FORMA(20),UFORMA(30),VFORMA(30),
     .    WFORMA(30),UFORM2(30),VFORM2(30),WFORM2(30)  
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
          DATA UFORM2/'(//6','X,''M','(X)''',',10X,','''SD','(X)''',
     .    ',8X,','''RO-','1(X)',''',8X',',''SD','-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA VFORM2/'(//6','X,''M','(Y)''',',10X,','''SD','(Y)''',
     .    ',8X,','''RO-','1(Y)',''',8X',',''SD','-N(','Y)'',',
     .    '11X,','''N(Y',')''  ',',7X,', 
     .    '''N(Y',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA WFORM2/'(//6','X,''M','(Z)''',',10X,','''SD','(Z)''',
     .    ',8X,','''RO-','1(Z)',''',8X',',''SD','-N(','Z)'',',
     .    '11X,','''N(Z',')''  ',',7X,', 
     .    '''N(Z',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NU=0 
          NV=0 
          NW=0 
          ND=0           
          XMIN= 1.D30
          XMAX=-1.D30
          YMIN= 1.D30
          YMAX=-1.D30   
          ZMIN= 1.D30
          ZMAX=-1.D30  
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
          READ(CDATO,FORMA) NPT,XS,YS,ZS,US,VS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          X(ND)=XS*SCX 
          Y(ND)=YS*SCY
          Z(ND)=ZS*SCZ 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          W(ND)=WS*SCW 
          XMIN=DMIN1(XMIN,X(ND))
          XMAX=DMAX1(XMAX,X(ND))
          YMIN=DMIN1(YMIN,Y(ND))
          YMAX=DMAX1(YMAX,Y(ND))
          ZMIN=DMIN1(ZMIN,Z(ND))
          ZMAX=DMAX1(ZMAX,Z(ND))            
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPXU(NU)=ND  
          IPYU(NU)=ND
          IPZU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(4)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPXV(NV)=ND  
          IPYV(NV)=ND  
          IPZV(NV)=ND
          ENDIF
          IF(DATO(IPTO(5)).EQ.PTO.AND.DATO(IAST(3)).NE.AST) THEN
          NW=NW+1  
          IPXW(NW)=ND  
          IPYW(NW)=ND
          IPZW(NW)=ND  
          ENDIF
          GOTO 10  
C
20        CALL MEDVA2(U,IPXU,EU,S2U,SU,NU) 
          CALL MEDVA2(V,IPXV,EV,S2V,SV,NV) 
          CALL MEDVA2(W,IPXW,EW,S2W,SW,NW) 
          DO 30 I=1,NU 
          U(I)=U(I)-EU 
          V(I)=V(I)-EV 
          W(I)=W(I)-EW 
30        CONTINUE 
          CALL RTORDA(X,Y,Z,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,NU,ND)
          CALL RTORDA(X,Y,Z,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,NV,ND)
          CALL RTORDA(X,Y,Z,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,NW,ND)
          VOL=(XMAX-XMIN)*(YMAX-YMIN)*(ZMAX-ZMIN)
          UNI=DEXP(DLOG(VOL/ND)/3.D0)
          IF(IAUTO.EQ.1) THEN
          AINF=UNI/10.D0
          ASUP=5.D0*UNI
          STP=UNI/10.D0
          ENDIF            
          NSTP=IDINT((ASUP-AINF)/STP)+1
          RAGGIO=AINF 
          ALATO=2.D0*ASUP
          CALL COROS(X,Y,Z,U,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,NU,NSTP,
     .    RAGGIO,ALATO,GMMU,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          CALL COROS(X,Y,Z,V,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,NV,NSTP,
     .    RAGGIO,ALATO,GMMV,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          CALL COROS(X,Y,Z,W,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,NW,NSTP,
     .    RAGGIO,ALATO,GMMW,SRV,NSRV,IWD,NMAX,NWDMAX,ILI)
          OLATO=0.D0
          S2NSE=1.D30
          IF(ILI.EQ.1) THEN 
          WRITE(6,1200)
          ELSE
          WRITE(6,1201)
          ENDIF
          DO 60 ISTP=1,NSTP
          SIGMA2=(S2U-GMMU(ISTP))+(S2V-GMMV(ISTP))+(S2W-GMMW(ISTP)) 
          SIGMA2=SIGMA2/3.D0
          WRITE(6,1300) ISTP*RAGGIO,DSQRT(S2U-GMMU(ISTP)),  
     .    DSQRT(S2V-GMMV(ISTP)),DSQRT(S2W-GMMW(ISTP)),DSQRT(SIGMA2) 
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=2.D0*ISTP*RAGGIO  
          S2NSE=SIGMA2 
          SNU=DSQRT(S2U-GMMU(ISTP))
          SNV=DSQRT(S2V-GMMV(ISTP))
          SNW=DSQRT(S2W-GMMW(ISTP))
          GMU=GMMU(ISTP)
          GMV=GMMV(ISTP)
          GMW=GMMW(ISTP)
          ENDIF
60        CONTINUE 
C
          NUM=0
          NVM=0
          NWM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMS(X,Y,Z,U,IPXU,IPYU,IPZU,LPXU,LPYU,LPZU,IWD,
     .    SNU,OLATO,1,NU,ND,NUM,NWDMAX,NMAX,IAST,3,4,TOLL,ILI)
          CALL TRIMS(X,Y,Z,V,IPXV,IPYV,IPZV,LPXV,LPYV,LPZV,IWD,
     .    SNV,OLATO,2,NV,ND,NVM,NWDMAX,NMAX,IAST,4,3,TOLL,ILI)
          CALL TRIMS(X,Y,Z,W,IPXW,IPYW,IPZW,LPXW,LPYW,LPZW,IWD,
     .    SNW,OLATO,3,NW,ND,NWM,NWDMAX,NMAX,IAST,3,4,TOLL,ILI)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1100) OLATO/2.D0
          ELSE 
          WRITE(6,1101) OLATO/2.D0
          ENDIF 
          IF(NUM.LT.9) THEN
          UFORMA(20)='I1,/'
          UFORMA(26)='I8) '
          UFORM2(20)='I1,/'
          UFORM2(26)='I8) '
          ENDIF
          IF(NUM.GE.10.AND.NUM.LT.100) THEN
          UFORMA(20)='I2,/'
          UFORMA(26)='I9) '
          UFORM2(20)='I2,/'
          UFORM2(26)='I9) '
          ENDIF
          IF(NUM.GE.100) THEN  
          UFORMA(20)='I3,/'
          UFORMA(26)='I10)'
          UFORM2(20)='I3,/'
          UFORM2(26)='I10)'
          ENDIF
          IF(NVM.LT.9) THEN
          VFORMA(20)='I1,/'
          VFORMA(26)='I8) '
          VFORM2(20)='I1,/'
          VFORM2(26)='I8) '
          ENDIF
          IF(NVM.GE.10.AND.NVM.LT.100) THEN
          VFORMA(20)='I2,/'
          VFORMA(26)='I9) '
          VFORM2(20)='I2,/'
          VFORM2(26)='I9) '
          ENDIF
          IF(NVM.GE.100) THEN  
          VFORMA(20)='I3,/'
          VFORMA(26)='I10)'
          VFORM2(20)='I3,/'
          VFORM2(26)='I10)'
          ENDIF
          IF(NWM.LT.9) THEN
          WFORMA(20)='I1,/'
          WFORMA(26)='I8) '
          WFORM2(20)='I1,/'
          WFORM2(26)='I8) '
          ENDIF
          IF(NWM.GE.10.AND.NWM.LT.100) THEN
          WFORMA(20)='I2,/'
          WFORMA(26)='I9) '
          WFORM2(20)='I2,/'
          WFORM2(26)='I9) '
          ENDIF
          IF(NWM.GE.100) THEN  
          WFORMA(20)='I3,/'
          WFORMA(26)='I10)'
          WFORM2(20)='I3,/'
          WFORM2(26)='I10)'
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,UFORMA) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORMA) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM  
          WRITE(6,WFORMA) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM 
          ELSE 
          WRITE(6,UFORM2) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORM2) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM  
          WRITE(6,WFORM2) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM 
          ENDIF 
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'RAGGIO OTTIMALE = ',F10.3,//) 
1101      FORMAT(///5X,'OPTIMAL RADIUS = ',F10.3,//)  
1200      FORMAT(///14X,'RAGGIO',7X,'SQM-N(X)',7X,'SQM-N(Y)',7X,
     .    'SQM-N(Z)',3X,'SQM-N(XUYUZ)'//)  
1201      FORMAT(///14X,'RADIUS',8X,'SD-N(X)',8X,'SD-N(Y)',8X,
     .    'SD-N(Z)',4X,'SD-N(XUYUZ)'//)  
1300      FORMAT(5X,5F15.3)
          END           
C*****************************     OTTIT       
          SUBROUTINE OTTIT1(Z,W,IPZW,LPZW,GMMW,SRV,NSRV, 
     .    AINF,ASUP,STP,NMAX,SCZ,SCW,IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL
     .    ILI)     
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA,WFORMA,WFORM2 
          CHARACTER*120 CDATO
          DIMENSION Z(*),W(*),IPZW(*),LPZW(*),
     .    GMMW(*),SRV(*),NSRV(*),IPTO(*),IAST(*)  
          DIMENSION DATO(120),FORMA(20),WFORMA(30),WFORM2(30) 
          EQUIVALENCE (CDATO,DATO(1))
          DATA PTO/'.'/,AST/'*'/
          DATA WFORMA/'(//6','X,''M','(X)''',',9X,','''SQM','(X)''',
     .    ',8X,','''RO-','1(X)',''',7X',',''SQ','M-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA WFORM2/'(//6','X,''M','(X)''',',10X,','''SD','(X)''',
     .    ',8X,','''RO-','1(X)',''',8X',',''SD','-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NW=0 
          ND=0 
          ZMIN= 1.D30
          ZMAX=-1.D30             
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
          READ(CDATO,FORMA) NPT,TS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          Z(ND)=TS*SCZ 
          W(ND)=WS*SCW 
          ZMIN=DMIN1(ZMIN,Z(ND))
          ZMAX=DMAX1(ZMAX,Z(ND))            
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NW=NW+1  
          IPZW(NW)=ND  
          ENDIF
          GOTO 10  
C
20        CALL MEDVA2(W,IPZW,EW,S2W,SW,NW) 
          DO 30 I=1,ND 
          W(I)=W(I)-EW 
30        CONTINUE 
          CALL RIORDA(Z,IPZW,NW,ND,1)  
          CALL PUNTA(IPZW,LPZW,0,1,NW,ND)
          VOL=(ZMAX-ZMIN)
          UNI=(VOL/ND)
          IF(IAUTO.EQ.1) THEN
          AINF=UNI/10.D0
          ASUP=5.D0*UNI
          STP=UNI/10.D0
          ENDIF              
          NSTP=IDINT((ASUP-AINF)/STP)+1
          RAGGIO=AINF 
          ALATO=2.D0*ASUP
          CALL COROT(Z,W,IPZW,NW,NSTP,RAGGIO,ALATO,GMMW,SRV,NSRV,NMAX,
     .    ILI)
          OLATO=0.D0
          S2NSE=1.D30
          IF(ILI.EQ.1) THEN 
          WRITE(6,1200)
          ELSE
          WRITE(6,1201)
          ENDIF
          DO 60 ISTP=1,NSTP
          SIGMA2=(S2W-GMMW(ISTP)) 
          WRITE(6,1300) ISTP*RAGGIO,DSQRT(SIGMA2)
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=2.D0*ISTP*RAGGIO  
          S2NSE=SIGMA2 
          SNW=DSQRT(S2W-GMMW(ISTP))
          GMW=GMMW(ISTP) 
          ENDIF
60        CONTINUE 
C
          NWM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMT(Z,W,IPZW,LPZW,SNW,OLATO,NMAX,1,NW,ND,NWM,IAST,3,4,
     .    TOLL,ILI)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1100) OLATO/2.D0 
          ELSE
          WRITE(6,1101) OLATO/2.D0
          ENDIF 
          IF(NWM.LT.9) THEN
          WFORMA(20)='I1,/'
          WFORMA(26)='I8) '
          WFORM2(20)='I1,/'
          WFORM2(26)='I8) '
          ENDIF
          IF(NWM.GE.10.AND.NWM.LT.100) THEN
          WFORMA(20)='I2,/'
          WFORMA(26)='I9) '
          WFORM2(20)='I2,/'
          WFORM2(26)='I9) '
          ENDIF
          IF(NWM.GE.100) THEN  
          WFORMA(20)='I3,/'
          WFORMA(26)='I10)'
          WFORM2(20)='I3,/'
          WFORM2(26)='I10)'
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,WFORMA) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM 
          ELSE 
          WRITE(6,WFORM2) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM
          ENDIF  
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'SEMI-INTERVALLO OTTIMALE = ',F10.3,//) 
1101      FORMAT(///5X,'OPTIMAL SEMI-INTERVAL = ',F10.3,//) 
1200      FORMAT(///5X,'SEMI-INTERVALLO',7X,'SQM-N(X)'//)  
1201      FORMAT(///5X,'SEMI-INTERVAL',10X,'SD-N(X)'//)  
1300      FORMAT(5X,2F15.3)
          END 
C
          SUBROUTINE OTTIT2(Z,U,V,IPZU,IPZV,LPZU,LPZV,
     .    GMMU,GMMV,SRV,NSRV,AINF,ASUP,STP,  
     .    NMAX,SCZ,SCU,SCV,IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA,UFORMA,VFORMA,UFORM2,VFORM2  
          CHARACTER*120 CDATO
          DIMENSION Z(*),U(*),V(*),IPZU(*),IPZV(*),LPZU(*),LPZV(*),
     .    GMMU(*),GMMV(*),SRV(*),NSRV(*),IPTO(*),IAST(*)  
          DIMENSION DATO(120),FORMA(20),UFORMA(30),VFORMA(30),
     .    UFORM2(30),VFORM2(30)
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
          DATA UFORM2/'(//6','X,''M','(X)''',',10X,','''SD','(X)''',
     .    ',8X,','''RO-','1(X)',''',8X',',''SD','-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA VFORM2/'(//6','X,''M','(Y)''',',10X,','''SD','(Y)''',
     .    ',8X,','''RO-','1(Y)',''',8X',',''SD','-N(','Y)'',',
     .    '11X,','''N(Y',')''  ',',7X,', 
     .    '''N(Y',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NU=0 
          NV=0 
          ND=0
          ZMIN= 1.D30
          ZMAX=-1.D30   
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
          READ(CDATO,FORMA) NPT,TS,US,VS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          Z(ND)=TS*SCZ 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          ZMIN=DMIN1(ZMIN,Z(ND))
          ZMAX=DMAX1(ZMAX,Z(ND))            
          IF(DATO(IPTO(2)).EQ.PTO.AND.DATO(IAST(1)).NE.AST) THEN
          NU=NU+1  
          IPZU(NU)=ND  
          ENDIF
          IF(DATO(IPTO(3)).EQ.PTO.AND.DATO(IAST(2)).NE.AST) THEN
          NV=NV+1  
          IPZV(NV)=ND  
          ENDIF
          GOTO 10  
C
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
          VOL=(ZMAX-ZMIN)
          UNI=(VOL/ND)
          IF(IAUTO.EQ.1) THEN
          AINF=UNI/10.D0
          ASUP=5.D0*UNI
          STP=UNI/10.D0
          ENDIF             
          NSTP=IDINT((ASUP-AINF)/STP)+1
          RAGGIO=AINF 
          ALATO=2.D0*ASUP
          CALL COROT(Z,U,IPZU,NU,NSTP,RAGGIO,ALATO,GMMU,SRV,NSRV,NMAX,
     .    ILI)   
          CALL COROT(Z,V,IPZV,NV,NSTP,RAGGIO,ALATO,GMMV,SRV,NSRV,NMAX,
     .    ILI)
          OLATO=0.D0
          S2NSE=1.D30 
          IF(ILI.EQ.1) THEN
          WRITE(6,1200)
          ELSE
          WRITE(6,1201)
          ENDIF
          DO 60 ISTP=1,NSTP
          SIGMA2=(S2U-GMMU(ISTP))+(S2V-GMMV(ISTP))
          SIGMA2=SIGMA2/2.D0
          WRITE(6,1300) ISTP*RAGGIO,DSQRT(S2U-GMMU(ISTP)),  
     .    DSQRT(S2V-GMMV(ISTP)),DSQRT(SIGMA2)
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=2.D0*ISTP*RAGGIO  
          S2NSE=SIGMA2 
          SNU=DSQRT(S2U-GMMU(ISTP))
          SNV=DSQRT(S2V-GMMV(ISTP))
          GMU=GMMU(ISTP)
          GMV=GMMV(ISTP)
          ENDIF
60        CONTINUE 
C 
          NUM=0
          NVM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMT(Z,U,IPZU,LPZU,SNU,OLATO,NMAX,1,NU,ND,NUM,IAST,4,3,
     .    TOLL,ILI)
          CALL TRIMT(Z,V,IPZV,LPZV,SNV,OLATO,NMAX,2,NV,ND,NVM,IAST,3,4,
     .    TOLL,ILI)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1100) OLATO/2.D0
          ELSE 
          WRITE(6,1101) OLATO/2.D0
          ENDIF 
          IF(NUM.LT.9) THEN
          UFORMA(20)='I1,/'
          UFORMA(26)='I8) '
          UFORM2(20)='I1,/'
          UFORM2(26)='I8) '
          ENDIF
          IF(NUM.GE.10.AND.NUM.LT.100) THEN
          UFORMA(20)='I2,/'
          UFORMA(26)='I9) '
          UFORM2(20)='I2,/'
          UFORM2(26)='I9) '
          ENDIF
          IF(NUM.GE.100) THEN  
          UFORMA(20)='I3,/'
          UFORMA(26)='I10)'
          UFORM2(20)='I3,/'
          UFORM2(26)='I10)'
          ENDIF
          IF(NVM.LT.9) THEN
          VFORMA(20)='I1,/'
          VFORMA(26)='I8) '
          VFORM2(20)='I1,/'
          VFORM2(26)='I8) '
          ENDIF
          IF(NVM.GE.10.AND.NVM.LT.100) THEN
          VFORMA(20)='I2,/'
          VFORMA(26)='I9) '
          VFORM2(20)='I2,/'
          VFORM2(26)='I9) '
          ENDIF
          IF(NVM.GE.100) THEN  
          VFORMA(20)='I3,/'
          VFORMA(26)='I10)'
          VFORM2(20)='I3,/'
          VFORM2(26)='I10)'
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,UFORMA) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORMA) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM
          ELSE  
          WRITE(6,UFORM2) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORM2) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM
          ENDIF
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'SEMI-INTERVALLO OTTIMALE = ',F10.3,//) 
1101      FORMAT(///5X,'OPTIMAL SEMI-INTERVAL = ',F10.3,//) 
1200      FORMAT(///5X,'SEMI-INTERVALLO',7X,'SQM-N(X)',7X,'SQM-N(Y)',  
     .    5X,'SQM-N(XUY)',//)  
1201      FORMAT(///5X,'SEMI-INTERVAL',10X,'SD-N(X)',8X,'SD-N(Y)',  
     .    6X,'SD-N(XUY)',//)  
1300      FORMAT(5X,F15.3,3F15.3)  
          END 
C
          SUBROUTINE OTTIT3(Z,U,V,W,IPZU,IPZV,IPZW,LPZU,LPZV,LPZW, 
     .    GMMU,GMMV,GMMW,SRV,NSRV,AINF,ASUP,STP,
     .    NMAX,SCZ,SCU,SCV,SCW,IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,PTO,AST 
          CHARACTER*4 FORMA,UFORMA,VFORMA,WFORMA,UFORM2,VFORM2,WFORM2
          CHARACTER*120 CDATO
          DIMENSION Z(*),U(*),V(*),W(*),IPZU(*),IPZV(*),IPZW(*),
     .    LPZU(*),LPZV(*),LPZW(*),GMMU(*),GMMV(*),GMMW(*),  
     .    SRV(*),NSRV(*),IAST(*),IPTO(*)
          DIMENSION DATO(120),FORMA(20),UFORMA(30),VFORMA(30),
     .    WFORMA(30),UFORM2(30),VFORM2(30),WFORM2(30)
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
          DATA UFORM2/'(//6','X,''M','(X)''',',10X,','''SD','(X)''',
     .    ',8X,','''RO-','1(X)',''',8X',',''SD','-N(','X)'',',
     .    '11X,','''N(X',')''  ',',7X,', 
     .    '''N(X',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA VFORM2/'(//6','X,''M','(Y)''',',10X,','''SD','(Y)''',
     .    ',8X,','''RO-','1(Y)',''',8X',',''SD','-N(','Y)'',',
     .    '11X,','''N(Y',')''  ',',7X,', 
     .    '''N(Y',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          DATA WFORM2/'(//6','X,''M','(Z)''',',10X,','''SD','(Z)''',
     .    ',8X,','''RO-','1(Z)',''',8X',',''SD','-N(','Z)'',',
     .    '11X,','''N(Z',')''  ',',7X,', 
     .    '''N(Z',')-'',','I$,/','/4(F','10.3',',5X)',',I10',
     .    ',5X,','I$$)','    ','    ','    ','    '/
          NU=0 
          NV=0 
          NW=0 
          ND=0
          ZMIN= 1.D30
          ZMAX=-1.D30             
          REWIND 1 
10        READ(1,1000) (DATO(I),I=1,120)
          READ(CDATO,FORMA) NPT,TS,US,VS,WS
          IF(NPT.EQ.0) GOTO 20 
          ND=ND+1  
          Z(ND)=TS*SCZ 
          U(ND)=US*SCU 
          V(ND)=VS*SCV 
          W(ND)=WS*SCW 
          ZMIN=DMIN1(ZMIN,Z(ND))
          ZMAX=DMAX1(ZMAX,Z(ND))            
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
C
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
          VOL=(ZMAX-ZMIN)
          UNI=(VOL/ND)
          IF(IAUTO.EQ.1) THEN
          AINF=UNI/10.D0
          ASUP=5.D0*UNI
          STP=UNI/10.D0
          ENDIF              
          NSTP=IDINT((ASUP-AINF)/STP)+1
          RAGGIO=AINF
          ALATO=2.D0*ASUP
          CALL COROT(Z,U,IPZU,NU,NSTP,RAGGIO,ALATO,GMMU,SRV,NSRV,NMAX,
     .    ILI)
          CALL COROT(Z,V,IPZV,NV,NSTP,RAGGIO,ALATO,GMMV,SRV,NSRV,NMAX,
     .    ILI)
          CALL COROT(Z,W,IPZW,NW,NSTP,RAGGIO,ALATO,GMMW,SRV,NSRV,NMAX,
     .    ILI)
          OLATO=0.D0
          S2NSE=1.D30
          IF(ILI.EQ.1) THEN
          WRITE(6,1200)
          ELSE
          WRITE(6,1201)
          ENDIF
          DO 60 ISTP=1,NSTP
          SIGMA2=(S2U-GMMU(ISTP))+(S2V-GMMV(ISTP))+(S2W-GMMW(ISTP)) 
          SIGMA2=SIGMA2/3.D0
          WRITE(6,1300) ISTP*RAGGIO,DSQRT(S2U-GMMU(ISTP)),  
     .    DSQRT(S2V-GMMV(ISTP)),DSQRT(S2W-GMMW(ISTP)),DSQRT(SIGMA2) 
          IF(SIGMA2.LT.S2NSE) THEN 
          OLATO=2.D0*ISTP*RAGGIO  
          S2NSE=SIGMA2 
          SNU=DSQRT(S2U-GMMU(ISTP))
          SNV=DSQRT(S2V-GMMV(ISTP))
          SNW=DSQRT(S2W-GMMW(ISTP))
          GMU=GMMU(ISTP)
          GMV=GMMV(ISTP)
          GMW=GMMW(ISTP)
          ENDIF
60        CONTINUE 
C 
          NUM=0
          NVM=0
          NWM=0
          IF(ITRIM.EQ.1) THEN  
          CALL TRIMT(Z,U,IPZU,LPZU,SNU,OLATO,NMAX,1,NU,ND,NUM,IAST,3,4,
     .    TOLL,ILI)
          CALL TRIMT(Z,V,IPZV,LPZV,SNV,OLATO,NMAX,2,NV,ND,NVM,IAST,4,3,
     .    TOLL,ILI)
          CALL TRIMT(Z,W,IPZW,LPZW,SNW,OLATO,NMAX,3,NW,ND,NWM,IAST,3,4,
     .    TOLL,ILI)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,1100) OLATO/2.D0 
          ELSE
          WRITE(6,1101) OLATO/2.D0
          ENDIF 
          IF(NUM.LT.9) THEN
          UFORMA(20)='I1,/'
          UFORMA(26)='I8) '
          UFORM2(20)='I1,/'
          UFORM2(26)='I8) '
          ENDIF
          IF(NUM.GE.10.AND.NUM.LT.100) THEN
          UFORMA(20)='I2,/'
          UFORMA(26)='I9) '
          UFORM2(20)='I2,/'
          UFORM2(26)='I9) '
          ENDIF
          IF(NUM.GE.100) THEN  
          UFORMA(20)='I3,/'
          UFORMA(26)='I10)'
          UFORM2(20)='I3,/'
          UFORM2(26)='I10)'
          ENDIF
          IF(NVM.LT.9) THEN
          VFORMA(20)='I1,/'
          VFORMA(26)='I8) '
          VFORM2(20)='I1,/'
          VFORM2(26)='I8) '
          ENDIF
          IF(NVM.GE.10.AND.NVM.LT.100) THEN
          VFORMA(20)='I2,/'
          VFORMA(26)='I9) '
          VFORM2(20)='I2,/'
          VFORM2(26)='I9) '
          ENDIF
          IF(NVM.GE.100) THEN  
          VFORMA(20)='I3,/'
          VFORMA(26)='I10)'
          VFORM2(20)='I3,/'
          VFORM2(26)='I10)'
          ENDIF
          IF(NWM.LT.9) THEN
          WFORMA(20)='I1,/'
          WFORMA(26)='I8) '
          WFORM2(20)='I1,/'
          WFORM2(26)='I8) '
          ENDIF
          IF(NWM.GE.10.AND.NWM.LT.100) THEN
          WFORMA(20)='I2,/'
          WFORMA(26)='I9) '
          WFORM2(20)='I2,/'
          WFORM2(26)='I9) '
          ENDIF
          IF(NWM.GE.100) THEN  
          WFORMA(20)='I3,/'
          WFORMA(26)='I10)'
          WFORM2(20)='I3,/'
          WFORM2(26)='I10)'
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,UFORMA) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORMA) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM  
          WRITE(6,WFORMA) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM 
          ELSE 
          WRITE(6,UFORM2) NUM,EU,SU,GMU/S2U,SNU,NU,NU-NUM  
          WRITE(6,VFORM2) NVM,EV,SV,GMV/S2V,SNV,NV,NV-NVM  
          WRITE(6,WFORM2) NWM,EW,SW,GMW/S2W,SNW,NW,NW-NWM
          ENDIF 
          RETURN
1000      FORMAT(120A1)
1100      FORMAT(///5X,'SEMI-INTERVALLO OTTIMALE = ',F10.3,//) 
1101      FORMAT(///5X,'OPTIMAL SEMI-INTERVAL = ',F10.3,//) 
1200      FORMAT(///5X,'SEMI-INTERVALLO',7X,'SQM-N(X)',7X,'SQM-N(Y)',7X,
     .    'SQM-N(Z)',3X,'SQM-N(XUYUZ)',//) 
1201      FORMAT(///5X,'SEMI-INTERVAL',10X,'SD-N(X)',8X,'SD-N(Y)',8X,
     .    'SD-N(Z)',4X,'SD-N(XUYUZ)',//) 
1300      FORMAT(5X,5F15.3)
          END 
C*****************************     PACORO
          SUBROUTINE PACORO(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IAUTO,ITRIM, 
     .IQUO,TOLL,AINF,ASUP,STP,SCX,SCY,SCZ,SCU,SCV,SCW,FORMA) 
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
22        WRITE(*,1160) 
          READ(*,*,ERR=22) IFILE        
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
32        WRITE(*,1045)
          READ(*,*,ERR=32) IAUTO          
          AINF=0.D0
          ASUP=0.D0
          STP=0.D0
          IF(IAUTO.EQ.0) THEN
35        WRITE(*,1070)
          READ(*,*,ERR=35) NCOR,DELTA
          AINF=DELTA
          ASUP=NCOR*DELTA
          STP=DELTA
          ENDIF
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
          IF(ICASO.EQ.3) THEN  
          IF(ITIPO.EQ.1) THEN  
70        WRITE(*,1115)
          READ(*,*,ERR=70) SCX,SCY,SCZ,SCW 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
75        WRITE(*,1125)
          READ(*,*,ERR=75) SCX,SCY,SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
80        WRITE(*,1135)
          READ(*,*,ERR=80) SCX,SCY,SCZ,SCU,SCV,SCW 
          ENDIF
          ENDIF
          WRITE(*,1140)
          READ(*,1010) FORMA
          WRITE(6,1150)
          RETURN
1000      FORMAT(5X,'TITOLO')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASO:',/15X,'UNIDIMENSIONALE = 1',
     ./15X,'BIDIMENSIONALE  = 2',/15X,'TRIDIMENSIONALE = 3')  
1030      FORMAT(5X,'TIPO:',/15X,'UNA COMPONENTE = 1', 
     ./15X,'DUE COMPONENTI = 2',
     ./15X,'TRE COMPONENTI = 3')
1040      FORMAT(5X,'ISCRI:',/15X,'SCRIVE DATI INPUT = 1', 
     ./15X,'NON SCRIVE DATI INPUT = 0')
1045      FORMAT(5X,'IAUTO:',/15X,'OTTIMIZZAZIONE AUTOMATICA = 1', 
     ./15X,'OTTIMIZZAZIONE MANUALE = 0')
1050      FORMAT(5X,'ITRIM:',/15X,'ELIMINA ERRORI GROSSOLANI = 1', 
     ./15X,'NON ELIMINA ERRORI GROSSOLANI = 0')
1055      FORMAT(5X,'TOLL:',/15X,'VALORE CRITICO PER LA REIEZIONE ',
     .'DEGLI ERRORI GROSSOLANI')
1060      FORMAT(5X,'IQUO:',/15X,'PRESENZA QUOTE = 1', 
     ./15X,'ASSENZA QUOTE  = 0')
1070      FORMAT(5X,'NCOR DELTA:',/15X,'NUMERO DI PASSI' 
     ./15X,'PASSO DEL CORRELOGRAMMA')  
1080      FORMAT(5X,'SCALE:',/15X,'SCALA TEMPO T', 
     ./15X,'SCALA COMPONENTE W')
1090      FORMAT(5X,'SCALE:',/15X,'SCALA TEMPO T', 
     ./15X,'SCALE COMPONENTI U V') 
1100      FORMAT(5X,'SCALE:',/15X,'SCALA TEMPO T', 
     ./15X,'SCALE COMPONENTI U V W')
1110      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y',  
     ./15X,'SCALA COMPONENTE W')
1115      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y Z',  
     ./15X,'SCALA COMPONENTE W')
1120      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y',  
     ./15X,'SCALE COMPONENTI U V') 
1125      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y Z',  
     ./15X,'SCALE COMPONENTI U V')      
1130      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y',  
     ./15X,'SCALE COMPONENTI U V W')
1135      FORMAT(5X,'SCALE:',/15X,'SCALE COORDINATE X Y Z',  
     ./15X,'SCALE COMPONENTI U V W')
1140      FORMAT(5X,'FORMATO DATI INPUT')  
1150      FORMAT(5X,'DATI INPUT FILE 1.') 
1160      FORMAT(5X,'IFILE:',/15X,'SALVA DATI POTATI = 1', 
     .'     (DATI OUTPUT FILE 2.)',
     ./15X,'NON SALVA DATI POTATI = 0')
          END
C
          SUBROUTINE PICORO(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IAUTO,ITRIM, 
     .IQUO,TOLL,AINF,ASUP,STP,SCX,SCY,SCZ,SCU,SCV,SCW,FORMA) 
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
22        WRITE(*,1160) 
          READ(*,*,ERR=22) IFILE        
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
32        WRITE(*,1045)
          READ(*,*,ERR=32) IAUTO          
          AINF=0.D0
          ASUP=0.D0
          STP=0.D0
          IF(IAUTO.EQ.0) THEN
35        WRITE(*,1070)
          READ(*,*,ERR=35) NCOR,DELTA
          AINF=DELTA
          ASUP=NCOR*DELTA
          STP=DELTA
          ENDIF
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
          IF(ICASO.EQ.3) THEN  
          IF(ITIPO.EQ.1) THEN  
70        WRITE(*,1115)
          READ(*,*,ERR=70) SCX,SCY,SCZ,SCW 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
75        WRITE(*,1125)
          READ(*,*,ERR=75) SCX,SCY,SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
80        WRITE(*,1135)
          READ(*,*,ERR=80) SCX,SCY,SCZ,SCU,SCV,SCW 
          ENDIF
          ENDIF
          WRITE(*,1140)
          READ(*,1010) FORMA
          WRITE(6,1150)
          RETURN
1000      FORMAT(5X,'TITLE')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASE:',/15X,'ONEDIMENSIONAL   = 1',
     ./15X,'TWODIMENSIONAL   = 2',/15X,'THREEDIMENSIONAL = 3')  
1030      FORMAT(5X,'TYPE:',/15X,'ONE COMPONENT    = 1', 
     ./15X,'TWO COMPONENTS   = 2',
     ./15X,'THREE COMPONENTS = 3')
1040      FORMAT(5X,'ISCRI:',/15X,'PRINT INPUT DATA = 1', 
     ./15X,'NO PRINT INPUT DATA = 0')
1045      FORMAT(5X,'IAUTO:',/15X,'AUTOMATIC OPTIMIZATION = 1', 
     ./15X,'MANUAL OPTIMIZATION = 0')
1050      FORMAT(5X,'ITRIM:',/15X,'BLUNDER REJECTION = 1', 
     ./15X,'NO BLUNDER REJECTION = 0')
1055      FORMAT(5X,'TOLL:',/15X,'CRITICAL VALUE FOR BLUNDER ',
     .'REJECTION')
1060      FORMAT(5X,'IQUO:',/15X,'HEIGHTS = 1', 
     ./15X,'NO HEIGHTS  = 0')
1070      FORMAT(5X,'NCOR DELTA:',/15X,'NUMBER OF STEPS' 
     ./15X,'STEP OF THE COVARIANCE FUNCTIONS')  
1080      FORMAT(5X,'SCALE:',/15X,'SCALE OF THE TIME', 
     ./15X,'SCALE OF W COMPONENT')
1090      FORMAT(5X,'SCALE:',/15X,'SCALE OF THE TIME', 
     ./15X,'SCALES OF U V COMPONENTS') 
1100      FORMAT(5X,'SCALE:',/15X,'SCALE OF THE TIME', 
     ./15X,'SCALES OF U V W COMPONENTS')
1110      FORMAT(5X,'SCALE:',/15X,'SCALES OF X Y COORDINATES',  
     ./15X,'SCALE OF W COMPONENT')
1115      FORMAT(5X,'SCALE:',/15X,'SCALES OF X Y Z  COORDINATES',  
     ./15X,'SCALE OF W COMPONENT')
1120      FORMAT(5X,'SCALE:',/15X,'SCALES  OF X Y COORDINATES',  
     ./15X,'SCALES OF U V W COMPONENTS') 
1125      FORMAT(5X,'SCALE:',/15X,'SCALES OF X Y Z  COORDINATES',  
     ./15X,'SCALES OF U V COMPONENTS')      
1130      FORMAT(5X,'SCALE:',/15X,'SCALES OF X Y COORDINATES',  
     ./15X,'SCALES OF U V W COMPONENTS')
1135      FORMAT(5X,'SCALE:',/15X,'SCALES  OF X Y Z COORDINATES',  
     ./15X,'SCALES OF U V W COMPONENTS')
1140      FORMAT(5X,'FORMAT OF INPUT DATA')  
1150      FORMAT(5X,'INPUT DATA ON FILE 1.') 
1160      FORMAT(5X,'IFILE:',/15X,'SAVE TRIMMED DATA = 1', 
     .'     (OUTPUT DATA ON FILE 2.)',
     ./15X,'NO SAVE TRIMMED DATA = 0')
          END
C*****************************     PACORR
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
          IF(ICASO.EQ.3) THEN  
          IF(ITIPO.EQ.1) THEN  
70        WRITE(*,1115)
          READ(*,*,ERR=70) SCX,SCY,SCZ,SCW 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
75        WRITE(*,1125)
          READ(*,*,ERR=75) SCX,SCY,SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
80        WRITE(*,1135)
          READ(*,*,ERR=80) SCX,SCY,SCZ,SCU,SCV,SCW 
          ENDIF
          ENDIF
          WRITE(*,1140)
          READ(*,1010) FORMA
          WRITE(6,1150)
          RETURN
1000      FORMAT(5X,'TITOLO')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASO:'/15X,'UNIDIMENSIONALE = 1'
     ./15X,'BIDIMENSIONALE  = 2'/15x,'TRIDIMENSIONALE = 3')  
1030      FORMAT(5X,'TIPO:'/15X,'UNA COMPONENTE = 1' 
     ./15X,'DUE COMPONENTI = 2'
     ./15X,'TRE COMPONENTI = 3')
1040      FORMAT(5X,'ISCRI:'/15X,'SCRIVE DATI INPUT = 1' 
     ./15X,'NON SCRIVE DATI INPUT = 0')
1050      FORMAT(5X,'IFILE:'/15X,'SALVA CORRELOGRAMMI = 1',
     .5X,'(DATI OUTPUT FILE 2.)'
     ./15X,'NON SALVA CORRELOGRAMMI = 0')  
1060      FORMAT(5X,'IQUO:'/15X,'PRESENZA QUOTE = 1' 
     ./15X,'ASSENZA QUOTE  = 0')
1070      FORMAT(5X,'NCOR DELTA:',/15X,'NUMERO DI PASSI' 
     ./15X,'PASSO DEL CORRELOGRAMMA')  
1080      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO T'
     ./15X,'SCALA COMPONENTE W')
1090      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO T'
     ./15X,'SCALE COMPONENTI U V') 
1100      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO T'
     ./15X,'SCALE COMPONENTI U V W')
1110      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE X Y'  
     ./15X,'SCALA COMPONENTE W')
1115      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE X Y Z'  
     ./15X,'SCALA COMPONENTE W')
1120      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE X Y'  
     ./15X,'SCALE COMPONENTI U V')
1125      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE X Y Z'  
     ./15X,'SCALE COMPONENTI U V') 
1130      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE X Y'  
     ./15X,'SCALE COMPONENTI U V W')
1135      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE X Y Z' 
     ./15X,'SCALE COMPONENTI U V W')
1140      FORMAT(5X,'FORMATO DATI INPUT')  
1150      FORMAT(5X,'DATI INPUT FILE 1.')  
          END  
C
          SUBROUTINE PICORR(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IQUO,
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
          IF(ICASO.EQ.3) THEN  
          IF(ITIPO.EQ.1) THEN  
70        WRITE(*,1115)
          READ(*,*,ERR=70) SCX,SCY,SCZ,SCW 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
75        WRITE(*,1125)
          READ(*,*,ERR=75) SCX,SCY,SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
80        WRITE(*,1135)
          READ(*,*,ERR=80) SCX,SCY,SCZ,SCU,SCV,SCW 
          ENDIF
          ENDIF
          WRITE(*,1140)
          READ(*,1010) FORMA
          WRITE(6,1150)
          RETURN
1000      FORMAT(5X,'TITLE')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASE:'/15X,'ONEDIMENSIONAL   = 1'
     ./15X,'TWODIMENSIONAL   = 2'/15x,'THREEDIMENSIONAL = 3')  
1030      FORMAT(5X,'TYPE:'/15X,'ONE COMPONENT    = 1' 
     ./15X,'TWO COMPONENTS   = 2'
     ./15X,'THREE COMPONENTS = 3')
1040      FORMAT(5X,'ISCRI:'/15X,'PRINT INPUT DATA = 1' 
     ./15X,'NO PRINT INPUT DATA = 0')
1050      FORMAT(5X,'IFILE:'/15X,'SAVE COVARIANCE FUNCTIONS = 1',
     .5X,'(OUTPUT DATA ON FILE 2.)'
     ./15X,'NO SAVE COVARIANCE FUNCTIONS = 0')  
1060      FORMAT(5X,'IQUO:'/15X,'HEIGHTS = 1' 
     ./15X,'NO HEIGHTS  = 0')
1070      FORMAT(5X,'NCOR DELTA:',/15X,'NUMBER OF STEPS' 
     ./15X,'STEP OF THE COVARIANCE FUNCTIONS')  
1080      FORMAT(5X,'SCALE:'/15X,'SCALE OF THE TIME'
     ./15X,'SCALE OF W COMPONENT')
1090      FORMAT(5X,'SCALE:'/15X,'SCALE OF THE TIME'
     ./15X,'SCALES OF U V COMPONENTS') 
1100      FORMAT(5X,'SCALE:'/15X,'SCALE OF THE TIME'
     ./15X,'SCALES OF U V W COMPONENTS')
1110      FORMAT(5X,'SCALE:'/15X,'SCALES OF X Y COORDINATES'  
     ./15X,'SCALE OF W COMPONENT')
1115      FORMAT(5X,'SCALE:'/15X,'SCALES OF X Y Z COORDINATES'  
     ./15X,'SCALE OF W COMPONENT')
1120      FORMAT(5X,'SCALE:'/15X,'SCALES OF X Y COORDINATES'  
     ./15X,'SCALES OF U V COMPONENTS')
1125      FORMAT(5X,'SCALE:'/15X,'SCALES OF X Y Z COORDINATES'  
     ./15X,'SCALES OF U V COMPONENTS') 
1130      FORMAT(5X,'SCALE:'/15X,'SCALES OF X Y COORDINATES'  
     ./15X,'SCALES OF U V W COMPONENTS')
1135      FORMAT(5X,'SCALE:'/15X,'SCALES OF X Y Z COORDINATES' 
     ./15X,'SCALES OF U V W COMPONENTS')
1140      FORMAT(5X,'FORMAT OF INPUT DATA')  
1150      FORMAT(5X,'INPUT DATA ON FILE 1.')  
          END  
C*****************************     PAINTC
          SUBROUTINE PAINTC(TITOLO,ICASE,ITIPO,IFREE,ISCAL,IFILE, 
     .    XCOVT)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*4 TITOLO
          DIMENSION XCOVT(*),TITOLO(20)
          IFREE=0
          WRITE(*,1000)
          READ(*,1010) TITOLO  
1         WRITE(*,1020)
          READ(*,*,ERR=1) ICASE  
2         WRITE(*,1030)
          READ(*,*,ERR=2) ITIPO  
3         WRITE(*,1040)
          READ(*,*,ERR=3) IFILE
          IF(ICASE.EQ.3) THEN
4         WRITE(6,1080)
          READ(*,*,ERR=4) IFREE
          ENDIF  
5         WRITE(*,1050)
          READ(*,*,ERR=5) ISCAL  
          NTER=0
          IF(ITIPO.EQ.1) NTER=1
          IF(ITIPO.EQ.2) NTER=3
          IF(ITIPO.EQ.3) NTER=6
          IF(ISCAL.EQ.1) THEN  
6         WRITE(*,1060)
          READ(*,*,ERR=6) (XCOVT(I),I=1,NTER)
          DO 7 I=1,NTER
          XCOVT(I)=XCOVT(I)/2.D0
7         CONTINUE 
          ENDIF
          WRITE(6,1090)
          RETURN
1000      FORMAT(5X,'TITOLO')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASO:'/15X,'UNIDIMENSIONALE = 1'  
     ./15X,'BIDIMENSIONALE  = 2'/15X,'TRIDIMENSIONALE = 3')
1030      FORMAT(5X,'TIPO:'/15X,'UNA AUTOCOVARIANZA = 1'
     ./15X,'DUE AUTO E UNA CROSS COVARIANZA = 2'
     ./15X,'TRE AUTO E TRE CROSS COVARIANZE = 3')  
1040      FORMAT(5X,'IFILE:'
     ./15X,'SALVA VALORI DELLE FUNZIONI EMPIRICHE '
     .'E TEORICHE'/15X,'SALVA I PARAMETRI DELLE FUNZIONI TEORICHE = 1' 
     ./15X,'ALTRIMENTI = 0'/15X,'(DATI OUTPUT FILES 2. E 3.)')
1050      FORMAT(5X,'SCAL:'/15X,'PRESENZA FUNZIONE DI COVARIANZA FINITA'
     .' = 1'/15X,'ASSENZA FUNZIONE DI COVARIANZA FINITA = 0')  
1060      FORMAT(5X,'XCOVT:'/15X,'PUNTI DI AZZERAMENTO DELLA FUNZIONI'
     ./15X,' = 1 VALORE PER TIPO UGUALE A UNO'
     ./15X,' = 3 VALORI PER TIPO UGUALE A DUE' 
     ./15X,' = 6 VALORI PER TIPO UGUALE A TRE')
1080      FORMAT(5X,'IFREE:'
     ./15X,'USO, A SCELTA, F. TRIGONOMETRICHE = 1'
     ./30X,'F. DI BESSEL',7X,'= 2')
1090      FORMAT(5X,'DATI INPUT FILE 1.')       
          END  
C
          SUBROUTINE PIINTC(TITOLO,ICASE,ITIPO,IFREE,ISCAL,IFILE, 
     .    XCOVT)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*4 TITOLO
          DIMENSION XCOVT(*),TITOLO(20)
          IFREE=0
          WRITE(*,1000)
          READ(*,1010) TITOLO  
1         WRITE(*,1020)
          READ(*,*,ERR=1) ICASE  
2         WRITE(*,1030)
          READ(*,*,ERR=2) ITIPO  
3         WRITE(*,1040)
          READ(*,*,ERR=3) IFILE
          IF(ICASE.EQ.3) THEN
4         WRITE(6,1080)
          READ(*,*,ERR=4) IFREE
          ENDIF  
5         WRITE(*,1050)
          READ(*,*,ERR=5) ISCAL  
          NTER=0
          IF(ITIPO.EQ.1) NTER=1
          IF(ITIPO.EQ.2) NTER=3
          IF(ITIPO.EQ.3) NTER=6
          IF(ISCAL.EQ.1) THEN  
6         WRITE(*,1060)
          READ(*,*,ERR=6) (XCOVT(I),I=1,NTER)
          DO 7 I=1,NTER
          XCOVT(I)=XCOVT(I)/2.D0
7         CONTINUE 
          ENDIF
          WRITE(6,1090)
          RETURN
1000      FORMAT(5X,'TITLE')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASE:'/15X'ONEDIMENSIONAL   = 1'  
     ./15X,'TWODIMENSIONAL   = 2'/15X,'THREEDIMENSIONAL = 3')
1030      FORMAT(5X,'TYPE:'/15X,'ONE AUTOCOVARIANCE FUNCTION = 1'
     ./15X,'TWO AUTO AND THREE CROSS COVARIANCE FUNCTIONS   = 2'
     ./15X,'THREE AUTO AND THREE CROSS COVARIANCE FUNCTIONS = 3')  
1040      FORMAT(5X,'IFILE:'
     ./15X,'SAVE EMPYRICAL AND THEORICAL FUNCTIONS'/
     .15X,'SAVE THE PARAMETERS OF THE THOERICAL FUNCTIONS = 1' 
     ./15X,'OTHERWISE = 0'/15X,'(OUTPUT DATA ON FILES 2. AND 3.)')
1050      FORMAT(5X,'SCAL:'/15X,'FINITE COVARIANCE FUNCTION'
     .' = 1'/15X,'NO FINITE COVARIANCE FUNCTION = 0')  
1060      FORMAT(5X,'XCOVT:'/15X,'ZERO OF THE FUNCTIONS'
     ./15X,' = 1 VALUE  IF TYPE EQUAL TO ONE'
     ./15X,' = 3 VALUES IF TYPE EQUAL TO TWO' 
     ./15X,' = 6 VALUES IF TYPE EQUAL TO THREE')
1080      FORMAT(5X,'IFREE:'
     ./15X,'FREE USE, TRIGONOMETRIC FUNCTIONS = 1'
     ./25X,'BESSEL FUNCTIONS = 2')
1090      FORMAT(5X,'INPUT DATA ON FILE 1.')       
          END  
C*****************************     PAPOLI
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
32        WRITE(*,1070)
          READ(*,*,ERR=32) IGRA
          IF(IGRA.EQ.1) THEN
35        WRITE(*,1075)
          READ(*,*,ERR=35) NGRA
          ELSE
          IF(ICASO.EQ.1) NGRA=20
          IF(ICASO.EQ.2) NGRA=10
          IF(ICASO.EQ.3) NGRA=5
          ENDIF 
          NGRA=NGRA+1  
38        WRITE(*,1078)
          READ(*,*,ERR=38) TOLL 
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
          IF(ICASO.EQ.3) THEN  
          IF(ITIPO.EQ.1) THEN  
70        WRITE(*,1115)
          READ(*,*,ERR=70) SCX,SCY,SCZ,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
75        WRITE(*,1125)
          READ(*,*,ERR=75) SCX,SCY,SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
80        WRITE(*,1135)
          READ(*,*,ERR=80) SCX,SCY,SCZ,SCU,SCV,SCW 
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
          WRITE(6,1150)
          RETURN
1000      FORMAT(5X,'TITOLO')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASO:'/15X,'UNIDIMENSIONALE = 1'
     ./15X,'BIDIMENSIONALE  = 2',/15X,'TRIDIMENSIONALE = 3')  
1030      FORMAT(5X,'TIPO:'/15X,'UNA COMPONENTE = 1' 
     ./15X,'DUE COMPONENTI = 2'
     ./15X,'TRE COMPONENTI = 3')
1040      FORMAT(5X,'ISCRI:'/15X,'SCRIVE DATI INPUT = 1' 
     ./15X,'NON SCRIVE DATI INPUT = 0')
1050      FORMAT(5X,'IFILE:'/15X,'SALVA SCARTI-RESIDUI E ',
     .'COEFFICIENTI POLINOMIO = 1'/15X,'(DATI OUTPUT FILES 2. E 3.)' 
     ./15X,'NON SALVA SCARTI-RESIDUI E COEFFICIENTI POLINOMIO = 0')
1060      FORMAT(5X,'IQUO:'/15X,'PRESENZA QUOTE = 1' 
     ./15X,'ASSENZA QUOTE  = 0')
1070      FORMAT(5X,'IGRA:'/15X,
     .'IMPONE GRADO DEL POLINOMIO = 1'/15X,
     .'NON IMPONE GRADO DEL POLINOMIO = 0')
1075      FORMAT(5X,'NGRA:'/15X,
     .'GRADO DEL POLINOMIO: MAX = 20 SE CASO = 1'/36X,
     .'MAX = 10 SE CASO = 2'/36X,'MAX =  5 SE CASO = 3 ') 
1078      FORMAT(5X,'TOLL:'
     ./15X,'VALORE CRITICO SULLA SIGNIFICATIVITA DEI PARAMETRI'
     ./15X,'E PER LA REIEZIONE DEGLI ERRORI GROSSOLANI')
1080      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALA COMPONENTE X')
1090      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALE COMPONENTI X Y') 
1100      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALE COMPONENTI X Y Z')
1110      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE'
     ./15X,'SCALA COMPONENTE X')
1115      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE E QUOTA'
     ./15X,'SCALA COMPONENTE X')
1120      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE'
     ./15X,'SCALE COMPONENTI X Y') 
1125      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE E QUOTA'
     ./15X,'SCALE COMPONENTI X Y') 
1130      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE'
     ./15X,'SCALE COMPONENTI X Y Z')
1135      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE E QUOTA'
     ./15X,'SCALE COMPONENTI X Y Z')
1140      FORMAT(5X,'FORMATO DATI INPUT (I COMPONENTE)')
1142      FORMAT(5X,'FORMATO DATI INPUT (II COMPONENTE)')  
1143      FORMAT(5X,'FORMATO DATI INPUT (III COMPONENTE)') 
1150      FORMAT(5X,'DATI INPUT FILE 1.')  
          END  
C
          SUBROUTINE PIPOLI(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IQUO,
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
32        WRITE(*,1070)
          READ(*,*,ERR=32) IGRA
          IF(IGRA.EQ.1) THEN
35        WRITE(*,1075)
          READ(*,*,ERR=35) NGRA
          ELSE
          IF(ICASO.EQ.1) NGRA=20
          IF(ICASO.EQ.2) NGRA=10
          IF(ICASO.EQ.3) NGRA=5
          ENDIF 
          NGRA=NGRA+1  
38        WRITE(*,1078)
          READ(*,*,ERR=38) TOLL 
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
          IF(ICASO.EQ.3) THEN  
          IF(ITIPO.EQ.1) THEN  
70        WRITE(*,1115)
          READ(*,*,ERR=70) SCX,SCY,SCZ,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN  
75        WRITE(*,1125)
          READ(*,*,ERR=75) SCX,SCY,SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
80        WRITE(*,1135)
          READ(*,*,ERR=80) SCX,SCY,SCZ,SCU,SCV,SCW 
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
          WRITE(6,1150)
          RETURN
1000      FORMAT(5X,'TITLE')  
1010      FORMAT(20A4) 
1020      FORMAT(5X,'CASE:'/15X,'ONEDIMENSIONAL   = 1'
     ./15X,'TWODIMENSIONAL   = 2',/15X,'THREEDIMENSIONAL = 3')  
1030      FORMAT(5X,'TYPE:'/15X,'ONE COMPONENT    = 1' 
     ./15X,'TWO COMPONENTS   = 2'
     ./15X,'THREE COMPONENTS = 3')
1040      FORMAT(5X,'ISCRI:'/15X,'PRINT INPUT DATA = 1' 
     ./15X,'NO PRINT INPUT DATA = 0')
1050      FORMAT(5X,'IFILE:'/15X,'SAVE RESIDUALS AND COEFFICENTS = 1'
     ./15X,'(OUTPUT DATA ON FILES 2. AND 3.)' 
     ./15X,'NO SAVE RESIDUALS AND COEFFICENTS = 0')
1060      FORMAT(5X,'IQUO:'/15X,'HEIGHTS = 1' 
     ./15X,'NO HEIGHTS  = 0')
1070      FORMAT(5X,'IGRA:'/15X,
     .'FIXED DEGREE FOR THE POLYNOMIAL = 1'/15X,
     .'FREE DEGREE FOR THE POLYNOMIAL = 0')
1075      FORMAT(5X,'NGRA:'/15X,
     .'DEGREE OF THE POLYNOMIAL: MAX = 20 IF CASE = 1'/41X,
     .'MAX = 10 IF CASE = 2'/41X,'MAX =  5 IF CASE = 3')  
1078      FORMAT(5X,'TOLL:'
     ./15X,'CRITICAL VALUE FOR PARAMETER TESTING ',
     .'AND OUTLIER REJECTION')
1080      FORMAT(5X,'SCALE:'/15X,'SCALE OF THE TIME' 
     ./15X,'SCALE OF X COMPONENT')
1090      FORMAT(5X,'SCALE:'/15X,'SCALE OF THE TIME' 
     ./15X,'SCALES OF X Y COMPONENTS') 
1100      FORMAT(5X,'SCALE:'/15X,'SCALE OF THE TIME' 
     ./15X,'SCALES OF X Y Z COMPONENTS')
1110      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES'
     ./15X,'SCALE OF X COMPONENT')
1115      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES '
     .'AND OF THE HEIGHT'
     ./15X,'SCALE  OF X COMPONENT')
1120      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES'
     ./15X,'SCALES OF X Y COMPONENTS') 
1125      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES '
     .'AND OF THE HEIGHT'
     ./15X,'SCALES OF X Y COMPONENTS') 
1130      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES'
     ./15X,'SCALES OF X Y Z COMPONENTS')
1135      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES ' 
     .'AND OF THE HEIGHT'
     ./15X,'SCALES OF X Y Z COMPONENTS')
1140      FORMAT(5X,'FORMAT OF INPUT DATA (1ST COMPONENT)')
1142      FORMAT(5X,'FORMAT OF INPUT DATA (2ND COMPONENT)')  
1143      FORMAT(5X,'FORMAT OF INPUT DATA (3RD COMPONENT)') 
1150      FORMAT(5X,'INPUT DATA ON FILE 1.')  
          END  
C*****************************     PASPLI
          SUBROUTINE PASPLI(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IQUO,  
     .IAUTO,ITABLE,TOLL,STP1,STP2,STP3,SCX,SCY,SCZ,SCU,SCV,SCW, 
     .FORMA,FORMA2,FORMA3) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*4 TITOLO,FORMA,FORMA2,FORMA3
          DIMENSION TITOLO(20),FORMA(20),FORMA2(20),FORMA3(20) 
          WRITE(*,1000)
          READ(*,1005) TITOLO  
5         WRITE(*,1025)
          READ(*,*,ERR=5) ITIPO
8         WRITE(*,1100)
          READ(*,*,ERR=8) ISCRI
10        WRITE(*,1035)
          READ(*,*,ERR=10) IFILE
15        WRITE(*,1030)
          READ(*,*,ERR=15) TOLL
          IF(ICASO.EQ.1) THEN 
          IQUO=0
          ITABLE=0
18        WRITE(*,1048)
          READ(*,*,ERR=18) IAUTO 
          IF(ITIPO.EQ.1) THEN 
          IF(IAUTO.EQ.0) THEN 
20        WRITE(*,1010)
          READ(*,*,ERR=20) STP1
          ELSE
          STP1=0.D0
          ENDIF
25        WRITE(*,1050)
          READ(*,*,ERR=25) SCZ,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN
          IF(IAUTO.EQ.0) THEN  
30        WRITE(*,1015)
          READ(*,*,ERR=30) STP1,STP2
          ELSE
          STP1=0.D0
          STP2=0.D0
          ENDIF
35        WRITE(*,1055)
          READ(*,*,ERR=35) SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN
          IF(IAUTO.EQ.0) THEN  
40        WRITE(*,1020)
          READ(*,*,ERR=40) STP1,STP2,STP3
          ELSE
          STP1=0.D0
          STP2=0.D0
          STP3=0.D0
          ENDIF  
45        WRITE(*,1060)
          READ(*,*,ERR=45) SCZ,SCU,SCV,SCW 
          ENDIF
          ENDIF
          IF(ICASO.EQ.2) THEN  
50        WRITE(*,1045)
          READ(*,*,ERR=50) IQUO
55        WRITE(*,1040)
          READ(*,*,ERR=55) ITABLE
58        WRITE(*,1048)
          READ(*,*,ERR=58) IAUTO   
          IF(ITIPO.EQ.1) THEN
          IF(IAUTO.EQ.0) THEN  
60        WRITE(*,1010)
          READ(*,*,ERR=60) STP1
          ELSE
          STP1=0.D0
          ENDIF
65        WRITE(*,1065)
          READ(*,*,ERR=65) SCX,SCY,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN 
          IF(IAUTO.EQ.0) THEN 
70        WRITE(*,1015)
          READ(*,*,ERR=70) STP1,STP2
          ELSE
          STP1=0.D0
          STP2=0.D0
          ENDIF
75        WRITE(*,1070)
          READ(*,*,ERR=75) SCX,SCY,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN
          IF(IAUTO.EQ.0) THEN  
80        WRITE(*,1020)
          READ(*,*,ERR=80) STP1,STP2,STP3
          ELSE
          STP1=0.D0
          STP2=0.D0
          STP3=0.D0
          ENDIF  
85        WRITE(*,1075)
          READ(*,*,ERR=85) SCX,SCY,SCU,SCV,SCW 
          ENDIF
          ENDIF
          IF(ICASO.EQ.3) THEN
          IQUO=0  
90        WRITE(*,1040)
          READ(*,*,ERR=90) ITABLE 
98        WRITE(*,1048)
          READ(*,*,ERR=98) IAUTO    
          IF(ITIPO.EQ.1) THEN
          IF(IAUTO.EQ.0) THEN  
95        WRITE(*,1015)
          READ(*,*,ERR=95) STP1
          ELSE
          STP1=0.D0
          ENDIF
100       WRITE(*,1068)
          READ(*,*,ERR=100) SCX,SCY,SCZ,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN
          IF(IAUTO.EQ.0) THEN  
105       WRITE(*,1010)
          READ(*,*,ERR=105) STP1,STP2
          ELSE
          STP1=0.D0
          STP2=0.D0
          ENDIF
110       WRITE(*,1072)
          READ(*,*,ERR=110) SCX,SCY,SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN
          IF(IAUTO.EQ.0) THEN  
115       WRITE(*,1020)
          READ(*,*,ERR=115) STP1,STP2,STP3
          ELSE
          STP1=0.D0
          STP2=0.D0
          STP3=0.D0
          ENDIF  
120       WRITE(*,1078)
          READ(*,*,ERR=120) SCX,SCY,SCZ,SCU,SCV,SCW 
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
1010      FORMAT(5X,'PASSO DELLE F. SPLINES (UN VALORE):')
1015      FORMAT(5X,'PASSO DELLE F. SPLINES (DUE VALORI):')  
1020      FORMAT(5X,'PASSO DELLE F. SPLINES (TRE VALORI):')  
1025      FORMAT(5X,'TIPO:'/15X,'UNA COMPONENTE = 1' 
     ./15X,'DUE COMPONENTI = 2'
     ./15X,'TRE COMPONENTI = 3')
1030      FORMAT(5X,'TOLL:'/15X,'VALORE CRITICO PER LA REIEZIONE ',
     .'DEGLI ERRORI GROSSOLANI')
1035      FORMAT(5X,'IFILE:'/15X,'SALVA SCARTI-RESIDUI E ',
     .'COEFFICIENTI F. SPLINES = 1'/15X,'(DATI OUTPUT FILES 2. E 3.)'
     ./15X,'NON SALVA SCARTI-RESIDUI E COEFFICIENTI F. SPLINES = 0')
1040      FORMAT(5X,'TABLE:'/15X,'SCRIVE TABELLE'  
     .' INCIDENZA E CONNESSIONE = 1'/15X,'NON SCRIVE ' 
     .'TABELLE INCIDENZA E CONNESSIONE = 0')
1045      FORMAT(5X,'IQUO:'/15X,'PRESENZA QUOTE = 1' 
     ./15X,'ASSENZA QUOTE  = 0')
1048      FORMAT(5X,'IAUTO:'
     ./15X,'INTERPOLAZIONE PREDEFINITA = (1),2,...,10,...' 
     ./15X,'(1: PASSO =  1 * DENSITA MEDIA: (?),'
     ./15X,' 2: PASSO =  2 * DENSITA MEDIA, ...,'
     ./15X,'10: PASSO = 10 * DENSITA MEDIA, ...)'
     ./15X,'INTERPOLAZIONE LIBERA      = 0')
1050      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALA COMPONENTE X')
1055      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALE COMPONENTI X Y') 
1060      FORMAT(5X,'SCALE:'/15X,'SCALA TEMPO' 
     ./15X,'SCALE COMPONENTI X Y Z')
1065      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE'
     ./15X,'SCALA COMPONENTE X')
1068      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE E QUOTA'
     ./15X,'SCALA COMPONENTE X')
1070      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE'
     ./15X,'SCALE COMPONENTI X Y') 
1072      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE E QUOTA'
     ./15X,'SCALE COMPONENTI X Y') 
1075      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE'
     ./15X,'SCALE COMPONENTI X Y Z')
1078      FORMAT(5X,'SCALE:'/15X,'SCALE COORDINATE PIANE E QUOTA'
     ./15X,'SCALE COMPONENTI X Y Z')
1080      FORMAT(5X,'FORMATO DATI INPUT (I COMPONENTE)')
1085      FORMAT(5X,'FORMATO DATI INPUT (II COMPONENTE)')  
1090      FORMAT(5X,'FORMATO DATI INPUT (III COMPONENTE)') 
1095      FORMAT(5X,'DATI INPUT FILE 1.')  
1100      FORMAT(5X,'ISCRI:'/15X,'SCRIVE DATI INPUT =1'
     ./15X,'NON SCRIVE DATI INPUT = 0')
          END  
C
          SUBROUTINE PISPLI(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IQUO,  
     .IAUTO,ITABLE,TOLL,STP1,STP2,STP3,SCX,SCY,SCZ,SCU,SCV,SCW, 
     .FORMA,FORMA2,FORMA3) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*4 TITOLO,FORMA,FORMA2,FORMA3
          DIMENSION TITOLO(20),FORMA(20),FORMA2(20),FORMA3(20) 
          WRITE(*,1000)
          READ(*,1005) TITOLO  
5         WRITE(*,1025)
          READ(*,*,ERR=5) ITIPO
8         WRITE(*,1100)
          READ(*,*,ERR=8) ISCRI
10        WRITE(*,1035)
          READ(*,*,ERR=10) IFILE
15        WRITE(*,1030)
          READ(*,*,ERR=15) TOLL
          IF(ICASO.EQ.1) THEN 
          IQUO=0
          ITABLE=0
18        WRITE(*,1048)
          READ(*,*,ERR=18) IAUTO 
          IF(ITIPO.EQ.1) THEN 
          IF(IAUTO.EQ.0) THEN 
20        WRITE(*,1010)
          READ(*,*,ERR=20) STP1
          ELSE
          STP1=0.D0
          ENDIF
25        WRITE(*,1050)
          READ(*,*,ERR=25) SCZ,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN
          IF(IAUTO.EQ.0) THEN  
30        WRITE(*,1015)
          READ(*,*,ERR=30) STP1,STP2
          ELSE
          STP1=0.D0
          STP2=0.D0
          ENDIF
35        WRITE(*,1055)
          READ(*,*,ERR=35) SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN
          IF(IAUTO.EQ.0) THEN  
40        WRITE(*,1020)
          READ(*,*,ERR=40) STP1,STP2,STP3
          ELSE
          STP1=0.D0
          STP2=0.D0
          STP3=0.D0
          ENDIF  
45        WRITE(*,1060)
          READ(*,*,ERR=45) SCZ,SCU,SCV,SCW 
          ENDIF
          ENDIF
          IF(ICASO.EQ.2) THEN  
50        WRITE(*,1045)
          READ(*,*,ERR=50) IQUO
55        WRITE(*,1040)
          READ(*,*,ERR=55) ITABLE
58        WRITE(*,1048)
          READ(*,*,ERR=58) IAUTO   
          IF(ITIPO.EQ.1) THEN
          IF(IAUTO.EQ.0) THEN  
60        WRITE(*,1010)
          READ(*,*,ERR=60) STP1
          ELSE
          STP1=0.D0
          ENDIF
65        WRITE(*,1065)
          READ(*,*,ERR=65) SCX,SCY,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN 
          IF(IAUTO.EQ.0) THEN 
70        WRITE(*,1015)
          READ(*,*,ERR=70) STP1,STP2
          ELSE
          STP1=0.D0
          STP2=0.D0
          ENDIF
75        WRITE(*,1070)
          READ(*,*,ERR=75) SCX,SCY,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN
          IF(IAUTO.EQ.0) THEN  
80        WRITE(*,1020)
          READ(*,*,ERR=80) STP1,STP2,STP3
          ELSE
          STP1=0.D0
          STP2=0.D0
          STP3=0.D0
          ENDIF  
85        WRITE(*,1075)
          READ(*,*,ERR=85) SCX,SCY,SCU,SCV,SCW 
          ENDIF
          ENDIF
          IF(ICASO.EQ.3) THEN
          IQUO=0  
90        WRITE(*,1040)
          READ(*,*,ERR=90) ITABLE 
98        WRITE(*,1048)
          READ(*,*,ERR=98) IAUTO    
          IF(ITIPO.EQ.1) THEN
          IF(IAUTO.EQ.0) THEN  
95        WRITE(*,1015)
          READ(*,*,ERR=95) STP1
          ELSE
          STP1=0.D0
          ENDIF
100       WRITE(*,1068)
          READ(*,*,ERR=100) SCX,SCY,SCZ,SCU 
          ENDIF
          IF(ITIPO.EQ.2) THEN
          IF(IAUTO.EQ.0) THEN  
105       WRITE(*,1010)
          READ(*,*,ERR=105) STP1,STP2
          ELSE
          STP1=0.D0
          STP2=0.D0
          ENDIF
110       WRITE(*,1072)
          READ(*,*,ERR=110) SCX,SCY,SCZ,SCU,SCV 
          ENDIF
          IF(ITIPO.EQ.3) THEN
          IF(IAUTO.EQ.0) THEN  
115       WRITE(*,1020)
          READ(*,*,ERR=115) STP1,STP2,STP3
          ELSE
          STP1=0.D0
          STP2=0.D0
          STP3=0.D0
          ENDIF  
120       WRITE(*,1078)
          READ(*,*,ERR=120) SCX,SCY,SCZ,SCU,SCV,SCW 
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
1000      FORMAT(5X,'TITLE')  
1005      FORMAT(20A4) 
1010      FORMAT(5X,'STEP OF THE SPLINE FUNCTIONS (ONE VALUE):')
1015      FORMAT(5X,'STEP OF THE SPLINE FUNCTIONS (TWO VALUES):')
1020      FORMAT(5X,'STEP OF THE SPLINE FUNCTIONS (THREE VALUES):') 
1025      FORMAT(5X,'TYPE:'/15X,'ONE COMPONENT    = 1' 
     ./15X,'TWO COMPONENTS   = 2'
     ./15X,'THREE COMPONENTS = 3')
1030      FORMAT(5X,'TOLL:'/15X,'CRITICAL VALUE FOR OUTLIER ',
     .'REJECTION')
1035      FORMAT(5X,'IFILE:'/15X,'SAVE RESIDUALS AND COEFFICENTS ',
     .'= 1'/15X,'(OUTPUT DATA ON FILES 2. E 3.)'
     ./15X,'NO SAVE RESIDUALS AND COEFFICENTS = 0')
1040      FORMAT(5X,'TABLE:'/15X,'PRINT POINT AND KNOT '  
     .'TABLES = 1'/15X,'NO PRINT POINT AND KNOT TABLES = 0') 
1045      FORMAT(5X,'IQUO:'/15X,'HEIGHTS = 1' 
     ./15X,'NO HEIGHTS  = 0')
1048      FORMAT(5X,'IAUTO:'
     ./15X,'FIXED INTERPOLATION = (1),2,...,10,...' 
     ./15X,'(1: STEP =  1 * AVERAGE DENSITY: (?),'
     ./15X,' 2: STEP =  2 * AVERAGE DENSITY, ...,'
     ./15X,'10: STEP = 10 * AVERAGE DENSITY, ...)'
     ./15X,'FREE INTERPOLATION  = 0')
1050      FORMAT(5X,'SCALE:'/15X,'SCALE OF THE TIME' 
     ./15X,'SCALE OF X COMPONENT')
1055      FORMAT(5X,'SCALE:'/15X,'SCALE OF THE TIME' 
     ./15X,'SCALES OF X Y COMPONENTS') 
1060      FORMAT(5X,'SCALE:'/15X,'SCALE OF THE TIME' 
     ./15X,'SCALES OF X Y Z COMPONENTS')
1065      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES'
     ./15X,'SCALE OF X COMPONENT')
1068      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES '
     .'AND OF THE HEIGTH'
     ./15X,'SCALE OF X COMPONENT')
1070      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES'
     ./15X,'SCALES OF X Y COMPONENTS') 
1072      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES '
     .'AND OF THE HEIGHT'
     ./15X,'SCALES OF X Y COMPONENTS') 
1075      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES'
     ./15X,'SCALES OF X Y Z COMPONENTS')
1078      FORMAT(5X,'SCALE:'/15X,'SCALES OF THE PLANE COORDINATES '
     .'AND OF THE HEIGHT'
     ./15X,'SCALES OF X Y Z COMPONENTS')
1080      FORMAT(5X,'FORMAT OF INPUT DATA (1ST COMPONENT)')
1085      FORMAT(5X,'FORMAT OF INPUT DATA (2ND COMPONENT)')  
1090      FORMAT(5X,'FORMAT OF INPUT DATA (3RD COMPONENT)') 
1095      FORMAT(5X,'INPUT DATA ON FILE 1.')  
1100      FORMAT(5X,'ISCRI:'/15X,'PRINT INPUT DATA =1'
     ./15X,'NO PRINT INPUT DATA = 0')
          END  
C*****************************     POLINP
          SUBROUTINE POLINP(X,Y,W,A,B,P,C,D,S,IGRA,NGRA,TOLL,  
     .    MMAX,NMAX,SCX,SCY,SCW,K,L,FORMA,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          CHARACTER*1 DATO,BLK,PTO,AST,SGN 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),W(1),A(MMAX,1),B(1),P(1),C(NMAX,1),D(1), 
     .    S(1)
          DIMENSION DATO(120),FORMA(20)
          EQUIVALENCE (CDATO,DATO(1))
          DATA BLK/' '/,PTO/'.'/,AST/'*'/  
          REWIND 1 
          I=0  
10        READ(1,500) (DATO(J),J=1,120)
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
          CALL MEDVAR(W,EW,S2W,SW,M)
          DO 30 I=1,M  
          X(I)=X(I)-EX 
          Y(I)=Y(I)-EY 
30        CONTINUE
          IF((NGRA-1).GT.10) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,700)
          ELSE
          WRITE(6,701) 
          ENDIF
          STOP
          ENDIF 
          IF(IGRA.EQ.0) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,720)
          ELSE
          WRITE(6,721)
          ENDIF
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
          IF(N.GT.M) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,700)
          ELSE
          WRITE(6,721)
          ENDIF
          STOP
          ENDIF 
          CALL NORMAL(A,B,P,C,D,MMAX,NMAX,M,N)
          C(N+1,N+1)=DBLE(ILI) 
          CALL TCHOL(C,NMAX,N) 
          CALL SOLS(C,D,NMAX,N)
          CALL SZERO(A,B,P,D,SQM,SQM2,MMAX,M,N,0)
          IF(ILI.EQ.1) THEN
          WRITE(6,710) L-1,SQM
          ELSE
          WRITE(6,711) L-1,SQM
          ENDIF  
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
          IF(N.GT.M) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,700)
          ELSE
          WRITE(6,721)
          ENDIF
          STOP
          ENDIF 
          CALL NORMAL(A,B,P,C,D,MMAX,NMAX,M,N) 
          F1=FMAT(IMAX,C,S,NMAX,N) 
          C(N+1,N+1)=DBLE(ILI)
          CALL TCHOL(C,NMAX,N) 
          CALL SOLS(C,D,NMAX,N)
          CALL SZERO(A,B,P,D,SQM,SQM2,MMAX,M,N,0)  
          CALL INVER(C,S,NMAX,N)
          F2=FMAT(IMAX,C,S,NMAX,N) 
          CALL DQVV(A,C,P,MMAX,NMAX,M,N)
          RAP=F1/F2
          IF(ILI.EQ.1) THEN
          WRITE(6,690) RAP,IMAX
          ELSE
          WRITE(6,691) RAP,IMAX
          ENDIF
          DO 90 I=1,N  
          S(I)=0.D0
          IF(SQM2*C(I,I).GT.1.D-10) S(I)=DSQRT(SQM2*C(I,I))
90        CONTINUE
          IF(ILI.EQ.1) THEN 
          WRITE(6,600) 
          WRITE(6,610)
          ELSE
          WRITE(6,601)
          WRITE(6,611)
          ENDIF 
          N=0  
          DO 110 I=1,L 
          DO 110 J=1,I 
          N=N+1
          ZS=0.D0
          IF(S(N).NE.0.D0) ZS=DABS(D(N)/S(N))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,620) D(N),S(N),I-1,J-1,I-J,ZS,SGN
110       CONTINUE 
          IF(IFILE.EQ.1) THEN  
          WRITE(3,630) EX,EY,L-1
          WRITE(3,640) (D(I),I=1,N)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,650) 
          ELSE
          WRITE(6,651)
          ENDIF
          DO 120 I=1,M 
          QQ=0.D0  
          IF(SQM2*P(I).GT.1.D-10) QQ=DSQRT(SQM2*P(I))  
          VAL=W(I)+B(I)
          SRV=1.D0 
          IF(SQM2*(1.D0-P(I)).GT.1.D-10) SRV=DSQRT(SQM2*(1.D0-P(I)))
          VINF=VAL-TOLL*SRV
          VSUP=VAL+TOLL*SRV
          ZS=0.D0
          IF(QQ.NE.0.D0) ZS=DABS(B(I))/QQ 
          SGN=BLK  
          IF(ZS.GT.TOLL) SGN=AST
          WRITE(6,660) I,W(I),VAL,SRV,B(I),QQ,SGN,VINF,VSUP
          W(I)=B(I)
120       CONTINUE 
          IF(ILI.EQ.1) THEN
          WRITE(6,670) M,N,L-1,SQM,SW 
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,680) ALFA
          ELSE
          WRITE(6,671) M,N,L-1,SQM,SW 
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,681) ALFA
          ENDIF
          RETURN
500       FORMAT(120A1)
600       FORMAT(///5X,'INTERPOLAZIONE POLINOMIALE')
601       FORMAT(///5X,'POLYNOMIAL INTERPOLATION')
610       FORMAT(///6X,'PARAMETRI',7X,'SQM',10X,'GRADI',7X,'Z N(0,1)'/) 
611       FORMAT(///5X,'PARAMETERS',8X,'SD',8X,'DEGREES',7X,'Z N(0,1)'/) 
620       FORMAT(5X,2D10.3,3X,I2,3X,I2,3X,I2,3X,F10.3,1X,A1)
630       FORMAT(2D20.12,I5)
640       FORMAT(6D20.12)  
650       FORMAT(///,8X,'N.',1X,'VALORE',6X,'VALORE',8X,'SQM', 
     .    4X,'SCARTO',7X,'SQM',5X,'INTERVALLI FIDUCIARI'/  
     .    11X,'OSSERVATO',3X,'STIMATO'/)
651       FORMAT(///,8X,'N.',2X,'OBSERVED',2X,'EXPECTED',8X,'SD', 
     .    2X,'RESIDUAL',9X,'SD',6X,'CONFIDENCE INTERVAL'/  
     .    12X,'VALUE',5X,'VALUE'/)
660       FORMAT(5X,I5,5F10.3,2X,A1,2X,2F10.3) 
670       FORMAT(///5X,'N. DI OSSERVAZIONI =',I5/  
     .    5X,'N. DI PARAMETRI    =',I5//
     .    5X,'GRADO =',I2,10X,'SIGMA ZERO   =',F9.3/
     .    24X,'SQM A PRIORI =',F9.3)
671       FORMAT(///5X,'N. OF OBSERVATIONS =',I5/  
     .    5X,'N. OF PARAMETERS   =',I5//
     .    5X,'DEGREE =',I2,10X,'SIGMA ZERO  =',F9.3/
     .    25X,'A PRIORI SD =',F9.3)
680       FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEI TESTS: ', 
     .    'ALFA = ',F5.3,' (SU DUE CODE)')  
681       FORMAT(/5X,'SIGNIFICANCE LEVEL OF TESTS: ', 
     .    'ALFA = ',F5.3,' (ON TWO SIDES)')  
690       FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
691       FORMAT(///5X,'CONDITION NUMBER', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
700       FORMAT(/5X,'GRADO DEL POLINOMIO ECCESSIVO')
701       FORMAT(/5X,'EXCESSIVE DEGREE OF THE POLYNOMIAL')
710       FORMAT(5X,I5,5X,F10.3)
711       FORMAT(6X,I5,5X,F10.3)
720      FORMAT(///5X,'GRADO',5X,'SIGMA ZERO'/)
721      FORMAT(///5X,'DEGREE',5X,'SIGMA ZERO'/)
          END  
C*****************************     POLINS
          SUBROUTINE POLINS(X,Y,Z,W,A,B,P,C,D,S,IGRA,NGRA,TOLL,  
     .    MMAX,NMAX,SCX,SCY,SCZ,SCW,K,L,FORMA,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          CHARACTER*1 DATO,BLK,PTO,AST,SGN 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION X(1),Y(1),Z(1),W(1),A(MMAX,1),B(1),P(1),C(NMAX,1),
     .    D(1),S(1)
          DIMENSION DATO(120),FORMA(20)
          EQUIVALENCE (CDATO,DATO(1))
          DATA BLK/' '/,PTO/'.'/,AST/'*'/  
          REWIND 1 
          I=0  
10        READ(1,500) (DATO(J),J=1,120)
          READ(CDATO,FORMA) NPT,XS,YS,ZS,WS
          IF(NPT.EQ.0) GOTO 20 
          IF(DATO(K).EQ.PTO.AND.DATO(L).NE.AST) THEN
          I=I+1
          X(I)=XS*SCX  
          Y(I)=YS*SCY  
          Z(I)=ZS*SCZ  
          W(I)=WS*SCW  
          ENDIF
          GOTO 10  
20        M=I  
          CALL MEDVAR(X,EX,S2X,SX,M)
          CALL MEDVAR(Y,EY,S2Y,SY,M)
          CALL MEDVAR(Z,EZ,S2Z,SZ,M)
          CALL MEDVAR(W,EW,S2W,SW,M)
          DO 30 I=1,M  
          X(I)=X(I)-EX 
          Y(I)=Y(I)-EY 
          Z(I)=Z(I)-EZ 
30        CONTINUE
          IF((NGRA-1).GT.5) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,700)
          ELSE
          WRITE(6,701)
          ENDIF
          STOP
          ENDIF 
          IF(IGRA.EQ.0) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,720)
          ELSE
          WRITE(6,721)
          ENDIF
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
          N=N+1
          A(K,N)=Y(K)**I
          N=N+1
          A(K,N)=Z(K)**I
          IF(I.NE.1) THEN  
          DO 32 J=1,I-1
          N=N+1
          A(K,N)=X(K)**J*Y(K)**(I-J)
32        CONTINUE 
          DO 33 J=1,I-1
          N=N+1
          A(K,N)=Y(K)**J*Z(K)**(I-J)
33        CONTINUE 
          DO 35 J=1,I-1
          N=N+1
          A(K,N)=Z(K)**J*X(K)**(I-J)
35        CONTINUE 
          ENDIF
          IF(I.GE.3) THEN
          DO 38 J=1,I-2
          DO 38 JS=1,J
          N=N+1
          A(K,N)=X(K)**(I-J-1)*Y(K)**(J-JS+1)*Z(K)**JS
38        CONTINUE
          ENDIF
40        CONTINUE 
          ENDIF
          B(K)=-W(K)
          P(K)=1.D0
50        CONTINUE 
          IF(N.GT.M) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,700)
          ELSE
          WRITE(6,701)
          ENDIF
          STOP
          ENDIF
          CALL NORMAL(A,B,P,C,D,MMAX,NMAX,M,N)
          C(N+1,N+1)=DBLE(ILI) 
          CALL TCHOL(C,NMAX,N) 
          CALL SOLS(C,D,NMAX,N)
          CALL SZERO(A,B,P,D,SQM,SQM2,MMAX,M,N,0)
          IF(ILI.EQ.1) THEN
          WRITE(6,710) L-1,SQM
          ELSE
          WRITE(6,711) L-1,SQM
          ENDIF  
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
          N=N+1
          A(K,N)=Y(K)**I
          N=N+1
          A(K,N)=Z(K)**I
          IF(I.NE.1) THEN  
          DO 62 J=1,I-1
          N=N+1
          A(K,N)=X(K)**J*Y(K)**(I-J)
62        CONTINUE 
          DO 63 J=1,I-1
          N=N+1
          A(K,N)=Y(K)**J*Z(K)**(I-J)
63        CONTINUE 
          DO 65 J=1,I-1
          N=N+1
          A(K,N)=Z(K)**J*X(K)**(I-J)
65       CONTINUE 
          ENDIF
          IF(I.GE.3) THEN
          DO 68 J=1,I-2
          DO 68 JS=1,J
          N=N+1
          A(K,N)=X(K)**(I-J-1)*Y(K)**(J-JS+1)*Z(K)**JS
68        CONTINUE
          ENDIF
70        CONTINUE 
          ENDIF
          B(K)=-W(K)
          P(K)=1.D0
80        CONTINUE         
          IF(N.GT.M) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,700)
          ELSE
          WRITE(6,701)
          ENDIF
          STOP
          ENDIF 
          CALL NORMAL(A,B,P,C,D,MMAX,NMAX,M,N) 
          F1=FMAT(IMAX,C,S,NMAX,N)
          C(N+1,N+1)=DBLE(ILI) 
          CALL TCHOL(C,NMAX,N) 
          CALL SOLS(C,D,NMAX,N)
          CALL SZERO(A,B,P,D,SQM,SQM2,MMAX,M,N,0)  
          CALL INVER(C,S,NMAX,N)
          F2=FMAT(IMAX,C,S,NMAX,N) 
          CALL DQVV(A,C,P,MMAX,NMAX,M,N)
          RAP=F1/F2
          IF(ILI.EQ.1) THEN
          WRITE(6,690) RAP,IMAX
          ELSE
          WRITE(6,691) RAP,IMAX
          ENDIF
          DO 90 I=1,N  
          S(I)=0.D0
          IF(SQM2*C(I,I).GT.1.D-10) S(I)=DSQRT(SQM2*C(I,I))
90        CONTINUE
          IF(ILI.EQ.1) THEN 
          WRITE(6,600) 
          WRITE(6,610)
          ELSE 
          WRITE(6,601) 
          WRITE(6,611) 
          ENDIF
          ZS=0.D0
          IF(S(1).NE.0.D0) ZS=DABS(D(1)/S(1))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,612) D(1),S(1),ZS,SGN
          IF(L.NE.1) THEN
          N=1  
          DO 110 I=1,L-1 
          N=N+1
          ZS=0.D0
          IF(S(N).NE.0.D0) ZS=DABS(D(N)/S(N))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,615) D(N),S(N),I,I,ZS,SGN
          N=N+1
          ZS=0.D0
          IF(S(N).NE.0.D0) ZS=DABS(D(N)/S(N))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,618) D(N),S(N),I,I,ZS,SGN
          N=N+1
          ZS=0.D0
          IF(S(N).NE.0.D0) ZS=DABS(D(N)/S(N))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,620) D(N),S(N),I,I,ZS,SGN
          IF(I.NE.1) THEN
          DO 102 J=1,I-1 
          N=N+1
          ZS=0.D0
          IF(S(N).NE.0.D0) ZS=DABS(D(N)/S(N))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,622) D(N),S(N),I,J,I-J,ZS,SGN
102       CONTINUE 
          DO 103 J=1,I-1 
          N=N+1
          ZS=0.D0
          IF(S(N).NE.0.D0) ZS=DABS(D(N)/S(N))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,623) D(N),S(N),I,J,I-J,ZS,SGN
103       CONTINUE 
          DO 105 J=1,I-1 
          N=N+1
          ZS=0.D0
          IF(S(N).NE.0.D0) ZS=DABS(D(N)/S(N))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,625) D(N),S(N),I,J,I-J,ZS,SGN
105       CONTINUE
          ENDIF
          IF(I.GE.3) THEN
          DO 108 J=1,I-1
          DO 108 JS=1,J 
          N=N+1
          ZS=0.D0
          IF(S(N).NE.0.D0) ZS=DABS(D(N)/S(N))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,628) D(N),S(N),I,I-J-1,J-JS+1,JS,ZS,SGN
108       CONTINUE
          ENDIF 
110       CONTINUE
          ENDIF 
          IF(IFILE.EQ.1) THEN  
          WRITE(3,630) EX,EY,L-1
          WRITE(3,640) (D(I),I=1,N)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,650) 
          ELSE
          WRITE(6,651)
          ENDIF
          DO 120 I=1,M 
          QQ=0.D0  
          IF(SQM2*P(I).GT.1.D-10) QQ=DSQRT(SQM2*P(I))  
          VAL=W(I)+B(I)
          SRV=1.D0 
          IF(SQM2*(1.D0-P(I)).GT.1.D-10) SRV=DSQRT(SQM2*(1.D0-P(I)))
          VINF=VAL-TOLL*SRV
          VSUP=VAL+TOLL*SRV
          ZS=0.D0
          IF(QQ.NE.0.D0) ZS=DABS(B(I))/QQ 
          SGN=BLK  
          IF(ZS.GT.TOLL) SGN=AST
          WRITE(6,660) I,W(I),VAL,SRV,B(I),QQ,SGN,VINF,VSUP
          W(I)=B(I)
120       CONTINUE
          IF(ILI.EQ.1) THEN 
          WRITE(6,670) M,N,L-1,SQM,SW 
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,680) ALFA
          ELSE
          WRITE(6,671) M,N,L-1,SQM,SW 
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,681) ALFA
          ENDIF
          RETURN
500       FORMAT(120A1)
600       FORMAT(///5X,'INTERPOLAZIONE POLINOMIALE')
601       FORMAT(///5X,'POLYNOMIAL INTERPOLATION')
610       FORMAT(///6X,'PARAMETRI',7X,'SQM',15X,'GRADI',7X,'Z N(0,1)'/) 
611       FORMAT(///5X,'PARAMETERS',8X,'SD',13X,'DEGREES',7X,'Z N(0,1)'
     .    /) 
612       FORMAT(5X,2D10.3,4X,'0',4X,'0',4X,'0',4X,'0',3X,F10.3,1X,A1)
615       FORMAT(5X,2D10.3,3X,I2,3X,I2,4X,'0',4X,'0',3X,F10.3,1X,A1)
618       FORMAT(5X,2D10.3,3X,I2,4X,'0',3X,I2,4X,'0',3X,F10.3,1X,A1)
620       FORMAT(5X,2D10.3,3X,I2,4X,'0',4X,'0',3X,I2,3X,F10.3,1X,A1)
622       FORMAT(5X,2D10.3,3X,I2,3X,I2,3X,I2,4X,'0',3X,F10.3,1X,A1)
623       FORMAT(5X,2D10.3,3X,I2,4X,'0',3X,I2,3X,I2,3X,F10.3,1X,A1)
625       FORMAT(5X,2D10.3,3X,I2,3X,I2,4X,'0',3X,I2,3X,F10.3,1X,A1)
628       FORMAT(5X,2D10.3,3X,I2,3X,I2,3X,I2,3X,I2,3X,F10.3,1X,A1)
630       FORMAT(2D20.12,I5)
640       FORMAT(6D20.12)  
650       FORMAT(///,8X,'N.',1X,'VALORE',6X,'VALORE',8X,'SQM', 
     .    4X,'SCARTO',7X,'SQM',5X,'INTERVALLI FIDUCIARI'/  
     .    11X,'OSSERVATO',3X,'STIMATO'/)
651       FORMAT(///,8X,'N.',2X,'OBSERVED',2X,'EXPECTED',8X,'SD', 
     .    2X,'RESIDUAL',8X,'SD',6X,'CONFIDENCE INTERVAL'/  
     .    12X,'VALUE',5X,'VALUE'/)
660       FORMAT(5X,I5,5F10.3,2X,A1,2X,2F10.3) 
670       FORMAT(///5X,'N. DI OSSERVAZIONI =',I5/  
     .    5X,'N. DI PARAMETRI    =',I5//
     .    5X,'GRADO =',I2,10X,'SIGMA ZERO =',F9.3/
     .    24X,'SQM A PRIORI =',F9.3)
671       FORMAT(///5X,'N. OF OBSERVATIONS =',I5/  
     .    5X,'N. OF PARAMETERS   =',I5//
     .    5X,'DEGREE =',I2,10X,'SIGMA ZERO  =',F9.3/
     .    25X,'A PRIORI SD =',F9.3)
680       FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEI TESTS: ', 
     .    'ALFA = ',F5.3,' (SU DUE CODE)')  
681       FORMAT(/5X,'SIGNIFICANCE LEVEL OF TESTS: ', 
     .    'ALFA = ',F5.3,' (ON TWO SIDES)')  
690       FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
691       FORMAT(///5X,'CONDITION NUMBER', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
700       FORMAT(/5X,'GRADO DEL POLINOMIO ECCESSIVO')
701       FORMAT(/5X,'EXCESSIVE DEGREE OF THE POLYNOMIAL')
710       FORMAT(5X,I5,5X,F10.3)
711       FORMAT(6X,I5,5X,F10.3)
720       FORMAT(///5X,'GRADO',5X,'SIGMA ZERO'/)
721       FORMAT(///5X,'DEGREE',5X,'SIGMA ZERO'/)
          END  
C*****************************     POLINT
          SUBROUTINE POLINT(Z,W,A,B,P,C,D,S,IGRA,NGRA,TOLL,
     .    MMAX,NMAX,SCZ,SCW,K,L,FORMA,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          CHARACTER*1 DATO,BLK,PTO,AST,SGN 
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION Z(1),W(1),A(MMAX,1),B(1),P(1),C(NMAX,1),D(1),  
     .    S(1)
          DIMENSION DATO(120),FORMA(20)
          EQUIVALENCE (CDATO,DATO(1))
          DATA BLK/' '/,PTO/'.'/,AST/'*'/  
          REWIND 1 
          I=0  
10        READ(1,500) (DATO(J),J=1,120)
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
          CALL MEDVAR(W,EW,S2W,SW,M)
          DO 30 I=1,M  
          Z(I)=Z(I)-EZ 
30        CONTINUE
          IF((NGRA-1).GT.20) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,700)
          ELSE
          WRITE(6,701) 
          ENDIF
          STOP
          ENDIF 
          IF(IGRA.EQ.0) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,720)
          ELSE
          WRITE(6,721)
          ENDIF
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
          IF(N.GT.M) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,700)
          ELSE
          WRITE(6,701)
          ENDIF
          STOP
          ENDIF 
          CALL NORMAL(A,B,P,C,D,MMAX,NMAX,M,N) 
          C(N+1,N+1)=DBLE(ILI)
          CALL TCHOL(C,NMAX,N) 
          CALL SOLS(C,D,NMAX,N)
          CALL SZERO(A,B,P,D,SQM,SQM2,MMAX,M,N,0)
          IF(ILI.EQ.1) THEN 
          WRITE(6,710) N-1,SQM
          ELSE
          WRITE(6,711) L-1,SQM
          ENDIF 
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
          IF(N.GT.M) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,700)
          ELSE
          WRITE(6,701)
          ENDIF
          STOP
          ENDIF  
          CALL NORMAL(A,B,P,C,D,MMAX,NMAX,M,N) 
          F1=FMAT(IMAX,C,S,NMAX,N)
          C(N+1,N+1)=DBLE(ILI) 
          CALL TCHOL(C,NMAX,N) 
          CALL SOLS(C,D,NMAX,N)
          CALL SZERO(A,B,P,D,SQM,SQM2,MMAX,M,N,0)  
          CALL INVER(C,S,NMAX,N)
          F2=FMAT(IMAX,C,S,NMAX,N) 
          CALL DQVV(A,C,P,MMAX,NMAX,M,N)
          RAP=F1/F2
          IF(ILI.EQ.1) THEN
          WRITE(6,690) RAP,IMAX
          ELSE
          WRITE(6,691) RAP,IMAX
          ENDIF
          DO 90 I=1,N  
          S(I)=0.D0
          IF(SQM2*C(I,I).GT.1.D-10) S(I)=DSQRT(SQM2*C(I,I))
90        CONTINUE
          IF(ILI.EQ.1) THEN 
          WRITE(6,600) 
          WRITE(6,610)
          ELSE   
          WRITE(6,601) 
          WRITE(6,611)
          ENDIF 
          DO 110 I=1,N
          ZS=0.D0 
          IF(S(I).NE.0.D0) ZS=DABS(D(I)/S(I))
          SGN=BLK  
          IF(ZS.LT.TOLL) SGN=AST
          WRITE(6,620) D(I),S(I),I-1,ZS,SGN
110       CONTINUE 
          IF(IFILE.EQ.1) THEN  
          WRITE(3,630) EZ,N-1  
          WRITE(3,640) (D(I),I=1,N)
          ENDIF
          IF(ILI.EQ.1) THEN
          WRITE(6,650)
          ELSE
          WRITE(6,651)
          ENDIF 
          DO 120 I=1,M 
          QQ=0.D0  
          IF(SQM2*P(I).GT.1.D-10) QQ=DSQRT(SQM2*P(I))  
          VAL=W(I)+B(I)
          SRV=1.D0 
          IF(SQM2*(1.D0-P(I)).GT.1.D-10) SRV=DSQRT(SQM2*(1.D0-P(I)))
          VINF=VAL-TOLL*SRV
          VSUP=VAL+TOLL*SRV
          ZS=0.D0
          IF(QQ.NE.0.D0) ZS=DABS(B(I))/QQ 
          SGN=BLK  
          IF(ZS.GT.TOLL) SGN=AST
          WRITE(6,660) I,W(I),VAL,SRV,B(I),QQ,SGN,VINF,VSUP
          W(I)=B(I)
120       CONTINUE
          IF(ILI.EQ.1) THEN 
          WRITE(6,670) M,N,N-1,SQM,SW 
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,680) ALFA
          ELSE
          WRITE(6,671) M,N,N-1,SQM,SW 
          ALFA=2.D0*(1.D0-DFN01(TOLL)) 
          WRITE(6,681) ALFA
          ENDIF
          RETURN
500       FORMAT(120A1)
600       FORMAT(///5X,'INTERPOLAZIONE POLINOMIALE')
601       FORMAT(///5X,'POLYNOMIAL INTERPOLATION')
610       FORMAT(///6X,'PARAMETRI',7X,'SQM',5X,'GRADO',7X,'Z N(0,1)'/)
611       FORMAT(///5X,'PARAMETERS',8X,'SD',3X,'DEGREES',7X,'Z N(0,1)'/) 
620       FORMAT(5X,2D10.3,8X,I2,3X,F10.3,1X,A1)
630       FORMAT(D20.12,I5)
640       FORMAT(6D20.12)  
650       FORMAT(///,8X,'N.',1X,'VALORE',6X,'VALORE',8X,'SQM', 
     .    4X,'SCARTO',7X,'SQM',5X,'INTERVALLI FIDUCIARI'/  
     .    11X,'OSSERVATO',3X,'STIMATO'/)
651       FORMAT(///,8X,'N.',2X,'OBSERVED',2X,'EXPECTED',8X,'SD', 
     .    2X,'RESIDUAL',8X,'SD',6X,'CONFIDENCE INTERVAL'/  
     .    12X,'VALUE',5X,'VALUE'/)
660       FORMAT(5X,I5,5F10.3,2X,A1,2X,2F10.3) 
670       FORMAT(///5X,'N. DI OSSERVAZIONI =',I5/  
     .    5X,'N. DI PARAMETRI    =',I5//
     .    5X,'GRADO =',I2,10X,'SIGMA ZERO   =',F9.3/
     .    24X,'SQM A PRIORI =',F9.3) 
671       FORMAT(///5X,'N. OF OBSERVATIONS =',I5/  
     .    5X,'N. OF PARAMETRES   =',I5//
     .    5X,'DEGREE =',I2,10X,'SIGMA ZERO   =',F9.3/
     .    25X,'A PRIORI SD =',F9.3) 
680       FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEI TESTS: ', 
     .    'ALFA = ',F5.3,' (SU DUE CODE)') 
681       FORMAT(/5X,'SIGNIFICANCE LEVEL OF TESTS: ', 
     .    'ALFA = ',F5.3,' (ON TWO SIDES)')  
690       FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
691       FORMAT(///5X,'CONDITION NUMBER', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
700       FORMAT(/5X,'GRADO DEL POLINOMIO ECCESSIVO')
701      FORMAT(/5X,'EXCESSIVE DEGREE OF THE POLYNOMIAL')
710       FORMAT(5X,I5,5X,F10.3)
711       FORMAT(6X,I5,5X,F10.3)
720       FORMAT(///5X,'GRADO',5X,'SIGMA ZERO'/)
721       FORMAT(///5X,'DEGREE',5X,'SIGMA ZERO'/)
          END

C
	  SUBROUTINE PAPREP(TITOLO,ITIPO,IFILE,IQUO,IPOLI,ISPLI,
     .ICOLL,FORMA)
	  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  CHARACTER*1 FORMA(80),TITOLO(80)

          WRITE(*,1005)
          READ(*,1000) TITOLO  
1         WRITE(*,1025)
          READ(*,*,ERR=1) ITIPO
2         WRITE(6,1027)
          READ(*,*,ERR=2) IFILE  
3         WRITE(*,1030)
          READ(*,*,ERR=3) IQUO
4         WRITE(*,1032)
          READ(*,*,ERR=4) IPOLI  
5         WRITE(*,1035)
          READ(*,*,ERR=5) ISPLI  
6         WRITE(*,1040)
          READ(*,*,ERR=6) ICOLL  
          WRITE(*,1045)
          READ(*,1000) FORMA
	  RETURN
1000      FORMAT(120A1)
1005      FORMAT(5X,'TITOLO ') 
1025      FORMAT(5X,'TIPO:',/15X,'UNA COMPONENTE = 1', 
     ./15X,'DUE COMPONENTI = 2',
     ./15X,'TRE COMPONENTI = 3')
1027      FORMAT(5X,'IFILE:',/15X,'SALVA RISULTATI = 1',
     .5X,'(DATI OUTPUT FILE 11.)'/
     .15X,'NON SALVA RISULTATI = 0')
1030      FORMAT(5X,'IQUO:',/15X,'PRESENZA QUOTE = 1', 
     ./15X,'ASSENZA QUOTE  = 0')
1032      FORMAT(5X,'PREDIZIONE DA INTERPOLAZIONE POLINOMIALE',
     .' (SI=1,NO=0)')  
1035      FORMAT(5X,'PREDIZIONE DA INTERPOLAZIONE CON F. SPLINES', 
     .' (SI=1,NO=0)')  
1040      FORMAT(5X,'PREDIZIONE DA COLLOCAZIONE (FILTRAGGIO)', 
     .' (SI=1,NO=0)') 
1045      FORMAT(5X,'FORMATO PUNTI DI PREDIZIONE O DI CONTROLLO') 
          END

	  SUBROUTINE PIPREP(TITOLO,ITIPO,IFILE,IQUO,IPOLI,ISPLI,
     .ICOLL,FORMA)
	  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  CHARACTER*1 FORMA(80),TITOLO(80)
          WRITE(*,1006)
          READ(*,1000) TITOLO  
7         WRITE(*,1026)
          READ(*,*,ERR=7) ITIPO
8         WRITE(6,1028)
          READ(*,*,ERR=8) IFILE  
9         WRITE(*,1031)
          READ(*,*,ERR=9) IQUO
10        WRITE(*,1033)
          READ(*,*,ERR=10) IPOLI  
11        WRITE(*,1036)
          READ(*,*,ERR=11) ISPLI  
12        WRITE(*,1041)
          READ(*,*,ERR=12) ICOLL  
          WRITE(*,1046)
          READ(*,1000) FORMA
	  RETURN
1000      FORMAT(120A1)
1006      FORMAT(5X,'TITLE ') 
1026      FORMAT(5X,'TYPE:',/15X,'ONE COMPONENT    = 1', 
     ./15X,'TWO COMPONENTS   = 2',
     ./15X,'THREE COMPONENTS = 3')
1028      FORMAT(5X,'IFILE:',/15X,'SAVE RESULTS = 1',
     .5X,'(OUTPUT DATA ON FILE 11.)'/
     .15X,'NO SAVE RESULTS = 0')
1031      FORMAT(5X,'IQUO:',/15X,'HEIGHTS = 1', 
     ./15X,'NO HEIGHTS = 0')
1033      FORMAT(5X,'PREDICTION BY POLYNOMIAL INTERPOLATION',
     .' (YES=1,NO=0)')  
1036      FORMAT(5X,'PREDICTION BY SPLINE INTERPOLATION', 
     .' (YES=1,NO=0)')  
1041      FORMAT(5X,'PREDICTION BY COLLOCATION FILTERING', 
     .' (YES=1,NO=0)')  
1046      FORMAT(5X,'FORMAT OF PREDICTED AND CONTROL POINTS')  
          END
c
c
	  SUBROUTINE PAPRES(TITOLO,IFILE,IPOLI,ISPLI,ICOLL,FORMA)
	  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  CHARACTER*1 FORMA(80),TITOLO(80)
          WRITE(*,1005)
          READ(*,1000) TITOLO
1         WRITE(*,1025)
          READ(*,*,ERR=1) ITIPO
2         WRITE(*,1027)
          READ(*,*,ERR=2) IFILE  
3         WRITE(*,1032)
          READ(*,*,ERR=3) IPOLI  
4         WRITE(*,1035)
          READ(*,*,ERR=4) ISPLI
5         WRITE(*,1040)
          READ(*,*,ERR=5) ICOLL  
          WRITE(*,1045)
          READ(*,1000) FORMA
	  RETURN
1000      FORMAT(120A1)
1005      FORMAT(5X,'TITOLO ')
1025      FORMAT(5X,'TIPO:',/15X,'UNA COMPONENTE = 1', 
     ./15X,'DUE COMPONENTI = 2',
     ./15X,'TRE COMPONENTI = 3')
1027      FORMAT(5X,'IFILE:',/15X,'SALVA RISULTATI = 1',
     .5X,'(DATI OUTPUT FILE 11.)'/
     .15X,'NON SALVA RISULTATI = 0')
1032      FORMAT(5X,'PREDIZIONE DA INTERPOLAZIONE POLINOMIALE',
     .' (SI=1,NO=0)')  
1035      FORMAT(5X,'PREDIZIONE DA INTERPOLAZIONE CON F. SPLINES', 
     .' (SI=1,NO=0)')  
1040      FORMAT(5X,'PREDIZIONE DA COLLOCAZIONE (FILTRAGGIO)', 
     .' (SI=1,NO=0)') 
1045      FORMAT(5X,'FORMATO PUNTI DI PREDIZIONE O DI CONTROLLO') 
          END
c
c
	  SUBROUTINE PIPRES(TITOLO,IFILE,IPOLI,ISPLI,ICOLL,FORMA)
	  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  CHARACTER*1 FORMA(80),TITOLO(80)
          WRITE(*,1006)
          READ(*,1000) TITOLO
6         WRITE(*,1026)
          READ(*,*,ERR=6) ITIPO
7         WRITE(*,1028)
          READ(*,*,ERR=7) IFILE  
8         WRITE(*,1033)
          READ(*,*,ERR=8) IPOLI  
9         WRITE(*,1036)
          READ(*,*,ERR=9) ISPLI
10        WRITE(*,1041)
          READ(*,*,ERR=10) ICOLL  
          WRITE(*,1046)
          READ(*,1000) FORMA
	  RETURN
1000      FORMAT(120A1)
1006      FORMAT(5X,'TITLE ') 
1026      FORMAT(5X,'TYPE:',/15X,'ONE COMPONENT    = 1', 
     ./15X,'TWO COMPONENTS   = 2',
     ./15X,'THREE COMPONENTS = 3')
1028      FORMAT(5X,'IFILE:',/15X,'SAVE RESULTS = 1',
     .5X,'(OUTPUT DATA ON FILE 11.)'/
     .15X,'NO SAVE RESULTS = 0')
1033      FORMAT(5X,'PREDICTION BY POLYNOMIAL INTERPOLATION',
     .' (YES=1,NO=0)')  
1036      FORMAT(5X,'PREDICTION BY SPLINE INTERPOLATION', 
     .' (YES=1,NO=0)')  
1041      FORMAT(5X,'PREDICTION BY COLLOCATION FILTERING', 
     .' (YES=1,NO=0)')  
1046      FORMAT(5X,'FORMAT OF PREDICTED AND CONTROL POINTS')  
          END
c
c
	  SUBROUTINE PAPRET(TITOLO,ITIPO,IFILE,IPOLI,ISPLI,
     .ICOLL,FORMA)
	  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  CHARACTER*1 FORMA(80),TITOLO(80)

          WRITE(*,1005)
          READ(*,1000) TITOLO
1         WRITE(*,1025)
          READ(*,*,ERR=1) ITIPO
2         WRITE(6,1027) 
          READ(*,*,ERR=2) IFILE  
3         WRITE(*,1032)
          READ(*,*,ERR=3) IPOLI  
4         WRITE(*,1035)
          READ(*,*,ERR=4) ISPLI  
5         WRITE(*,1040)
          READ(*,*,ERR=5) ICOLL  
          WRITE(*,1045)
          READ(*,1000) FORMA
	  RETURN
1000      FORMAT(120A1)
1005      FORMAT(5X,'TITOLO ') 
1025      FORMAT(5X,'TIPO:',/15X,'UNA COMPONENTE = 1', 
     ./15X,'DUE COMPONENTI = 2',
     ./15X,'TRE COMPONENTI = 3')
1027      FORMAT(5X,'IFILE:',/15X,'SALVA RISULTATI = 1',
     .5X,'(DATI OUTPUT FILE 11.)'/
     .15X,'NON SALVA RISULTATI = 0')
1032      FORMAT(5X,'PREDIZIONE DA INTERPOLAZIONE POLINOMIALE',
     .' (SI=1,NO=0)')  
1035      FORMAT(5X,'PREDIZIONE DA INTERPOLAZIONE CON F. SPLINES', 
     .' (SI=1,NO=0)')  
1040      FORMAT(5X,'PREDIZIONE DA COLLOCAZIONE (FILTRAGGIO)', 
     .' (SI=1,NO=0)') 
1045      FORMAT(5X,'FORMATO PUNTI DI PREDIZIONE O DI CONTROLLO') 
          END  

c
c
	  SUBROUTINE PIPRET(TITOLO,ITIPO,IFILE,IPOLI,ISPLI,
     .ICOLL,FORMA)
	  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  CHARACTER*1 FORMA(80),TITOLO(80)

          WRITE(*,1006)
          READ(*,1000) TITOLO
6         WRITE(*,1026)
          READ(*,*,ERR=6) ITIPO
7         WRITE(6,1028) 
          READ(*,*,ERR=7) IFILE  
8         WRITE(*,1033)
          READ(*,*,ERR=8) IPOLI  
9         WRITE(*,1036)
          READ(*,*,ERR=9) ISPLI  
10        WRITE(*,1041)
          READ(*,*,ERR=10) ICOLL  
          WRITE(*,1046)
          READ(*,1000) FORMA
	  RETURN
1000      FORMAT(120A1)
1006      FORMAT(5X,'TITLE ') 
1026      FORMAT(5X,'TYPE:',/15X,'ONE COMPONENT    = 1', 
     ./15X,'TWO COMPONENTS   = 2',
     ./15X,'THREE COMPONENTS = 3')
1028      FORMAT(5X,'IFILE:',/15X,'SAVE RESULTS = 1',
     .5X,'(OUTPUT DATA ON FILE 11.)'/
     .15X,'NO SAVE RESULTS = 0')
1033      FORMAT(5X,'PREDICTION BY POLYNOMIAL INTERPOLATION',
     .' (YES=1,NO=0)')  
1036      FORMAT(5X,'PREDICTION BY SPLINE INTERPOLATION', 
     .' (YES=1,NO=0)')  
1041      FORMAT(5X,'PREDICTION BY COLLOCATION FILTERING', 
     .' (YES=1,NO=0)')  
1046      FORMAT(5X,'FORMAT OF PREDICTED AND CONTROL POINTS')  
          END  
C
C*****************************     PRECLP
          SUBROUTINE PRECLP(NS,XS,YS,X,Y,Z,C,XCOV,INDIR,
     .IPX,IPY,LPX,LPY,IWD,NTIPO,KA,NMAX,NWDMAX,NWD,
     .N,DA,DB,DC,SN,NCASO,U,ILI)
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
          CALL WINDOW(XS,YS,AMPX,AMPY,X,Y,IPX,IPY,LPX, 
     .LPY,IWD,NMAX,NWDMAX,N,NWD,ILI)
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
     .DC,SN,NCASO,ICASO)
10        CONTINUE 
          DO 20 J=1,NWD
          IS=IWD(J)
          JS=INDIR(IS) 
          JJ=NTIPO(JS) 
          IC=IND(JJ)
          INI=KA(JS)+1 
          DO 20 I=1,IR 
          DO 20 K=1,IC 
          U(I)=U(I)+C(I,INI+K-1)*Z(INI+K-1)
          C(I,INI+K-1)=0.D0
20        CONTINUE 
          RETURN
          END  
C*****************************     PRECLS
          SUBROUTINE PRECLS(NS,XS,YS,ZS,X,Y,Z,H,C,XCOV,INDIR,
     .IPX,IPY,IPZ,LPX,LPY,LPZ,IWD,NTIPO,KA,NMAX,NWDMAX,NWD,
     .N,DA,DB,DC,SN,NCASO,U,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION H(1),C(3,1),X(1),Y(1),Z(1),
     .IPX(1),IPY(1),IPZ(1),LPX(1),LPY(1),LPZ(1),
     .IWD(1),NTIPO(1),KA(1),INDIR(1),
     .DA(1),DB(1),DC(1),SN(1),NCASO(1),U(1),INV(3),IND(7)  
          DATA INV/1,3,7/,IND/1,1,2,1,2,2,3/
          DO 1 I=1,NS  
          U(I)=0.D0
1         CONTINUE 
          AMPX=4.D0*XCOV
          AMPY=4.D0*XCOV
          AMPZ=4.D0*XCOV
          CALL POCKET(XS,YS,ZS,AMPX,AMPY,AMPZ,X,Y,Z, 
     .IPX,IPY,IPZ,LPX,LPY,LPZ,IWD,NMAX,NWDMAX,N,NWD,ILI)
          II=INV(NS)
          IR=NS
          DO 10 I=1,NWD
          IS=IWD(I)
          JS=INDIR(IS) 
          DELTA=DSQRT((XS-X(JS))**2+(YS-Y(JS))**2+
     .(ZS-Z(JS))**2) 
          JJ=NTIPO(JS) 
          IC=IND(JJ)
          INI=KA(JS)+1 
          DO 10 L=1,IR 
          DO 10 K=1,IC 
          C(L,INI+K-1)=GASS3(DELTA,XCOV,II,JJ,L,K,DA,DB,
     .DC,SN,NCASO,ICASO) 
10        CONTINUE 
          DO 20 J=1,NWD
          IS=IWD(J)
          JS=INDIR(IS) 
          JJ=NTIPO(JS) 
          IC=IND(JJ)
          INI=KA(JS)+1 
          DO 20 I=1,IR 
          DO 20 K=1,IC 
          U(I)=U(I)+C(I,INI+K-1)*H(INI+K-1)
          C(I,INI+K-1)=0.D0
20        CONTINUE 
          RETURN
          END  
C*****************************     PRECLT
          SUBROUTINE PRECLT(NS,TP,XCOVT,Z,TEMPO,C,IPUNT,
     .NTIPO,KA,NMAX,N,DA,DB,DC,SN,NCASO,U,ILI) 
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
          C(L,INC+K-1)=GASS1(DELTA,XCOVT,II,JJ,L,K,DA,DB,DC,SN,NCASO, 
     .ICASO) 
5         CONTINUE 
          INC=INC+IC
10        CONTINUE 
          DO 20 I=1,IR 
          DO 20 J=INI,IFI  
          U(I)=U(I)+C(I,J)*Z(J)
          C(I,J)=0.D0  
20        CONTINUE 
          RETURN
          END  
C*****************************     PREPLP
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
          IF(N.EQ.1) GOTO 10  
          DO 8 J=1,N-1 
          K=K+1
          U(I)=U(I)+COEF(I,K)*X**J 
          IF(J.EQ.1) GOTO 7
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
C*****************************     PREPLS
          SUBROUTINE PREPLS(NS,XS,YS,ZS,X0,Y0,Z0,COEF,NGR,U) 
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          DIMENSION COEF(3,1),X0(1),Y0(1),Z0(1),NGR(1),U(1)  
          DO 1 I=1,NS  
          U(I)=0.D0
1         CONTINUE 
          DO 12 I=1,NS 
          X=XS-X0(I)
          Y=YS-Y0(I)
          Z=ZS-Z0(I)
          N=NGR(I)+1
          K=1  
          U(I)=COEF(I,1)
          IF(N.EQ.1) GOTO 12  
          DO 10 J=1,N-1 
          K=K+1
          U(I)=U(I)+COEF(I,K)*X**J 
          K=K+1
          U(I)=U(I)+COEF(I,K)*Y**J 
          K=K+1
          U(I)=U(I)+COEF(I,K)*Z**J 
          IF(J.EQ.1) GOTO 7
          DO 2 L=1,J-1 
          K=K+1
          U(I)=U(I)+COEF(I,K)*X**L*Y**(J-L)
2         CONTINUE   
          DO 3 L=1,J-1 
          K=K+1
          U(I)=U(I)+COEF(I,K)*Y**L*Z**(J-L)
3         CONTINUE      
          DO 5 L=1,J-1 
          K=K+1
          U(I)=U(I)+COEF(I,K)*Z**L*X**(J-L)
5         CONTINUE 
7         IF(J.LT.3) GOTO 10
          DO 8 L=1,J-2
          DO 8 LS=1,L
          K=K+1
          U(I)=U(I)+COEF(I,K)*X**(J-L+1)*Y**(L-LS+1)*Z**LS
8         CONTINUE 
10        CONTINUE 
12        CONTINUE 
          RETURN
          END  
C*****************************     PREPLT
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
          IF(N.EQ.1) GOTO 10  
          DO 8 J=1,N-1 
          U(I)=U(I)+COEF(I,J+1)*T**J
8         CONTINUE 
10        CONTINUE 
          RETURN
          END  
C*****************************     PRESPP
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
          NJJ=NIC+(J-1)*NIX
          AY=DBLE(2-J) 
          JY=5-J
          VALY=SPLT(SRY+AY,JY) 
          DO 10 K=1,4  
          NKK=NJJ+K-1
          AX=DBLE(2-K) 
          JX=5-K
          U(I)=U(I)+COEF(I,NKK)*VALY*SPLT(SRX+AX,JX)
10        CONTINUE 
20        CONTINUE 
          RETURN
          END  
C*****************************     PRESPS
          SUBROUTINE PRESPS(NS,XS,YS,ZS,X0,Y0,Z0,STP,COEF,NI,N2,U)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION COEF(3,1),X0(1),Y0(1),Z0(1),STP(1),
     .NI(1),N2(1),U(1)
          DO 1 I=1,NS  
          U(I)=0.D0
1         CONTINUE 
          DO 20 I=1,NS 
          X=XS-X0(I)
          Y=YS-Y0(I)
          Z=ZS-Z0(I)
          ST=STP(I)
          NIX=NI(I)
          NIY=N2(I)
          NX=IDINT(X/ST)+1 
          NY=IDINT(Y/ST)+1 
          NZ=IDINT(Z/ST)+1 
          NC=(NZ-1)*(NIY-1)*(NIX-1)+(NY-1)*(NIX-1)+NX 
          SRX=(X-(NX-1)*ST)/ST 
          SRY=(Y-(NY-1)*ST)/ST 
          SRZ=(Z-(NZ-1)*ST)/ST 
          NIC=(NZ-2)*NIY*NIX+(NY-2)*NIX+NX-1  
          DO 10 II=1,4
          NII=NIC+(II-1)*NIY*NIX
          AZ=DBLE(2-II)
          JZ=5-II
          VALZ=SPLT(SRZ+AZ,JZ)
          DO 10 J=1,4  
          NJJ=NII+(J-1)*NIX
          AY=DBLE(2-J) 
          JY=5-J
          VALY=SPLT(SRY+AY,JY) 
          DO 10 K=1,4  
          NKK=NJJ+K-1
          AX=DBLE(2-K) 
          JX=5-K
          U(I)=U(I)+COEF(I,NKK)*VALZ*VALY*SPLT(SRX+AX,JX)
10        CONTINUE 
20        CONTINUE 
          RETURN
          END  
C*****************************     PRESPT
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
C*****************************     RCORR
          SUBROUTINE RCORR(X,Y,EX,EY,SX,SY,S2X,S2Y,SXY,RXY,N)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1)
          SRVX=0.D0
          SRVY=0.D0
          SRV2X=0.D0
          SRV2Y=0.D0
          SRVXY=0.D0
          DO 1 I=1,N  
          SRVX=SRVX+X(I)
          SRVY=SRVY+Y(I)
          SRV2X=SRV2X+X(I)**2  
          SRV2Y=SRV2Y+Y(I)**2  
          SRVXY=SRVXY+X(I)*Y(I)
1         CONTINUE 
          EX=SRVX/N  
          EY=SRVY/N  
          S2X=(SRV2X/N-EX**2)*N/(N-1)
          S2Y=(SRV2Y/N-EY**2)*N/(N-1)
          SX=DSQRT(S2X)
          SY=DSQRT(S2Y)
          SXY=SRVXY/N-EX*EY  
          RXY=SXY/(SX*SY)
          RETURN
          END
C*****************************     RCORR2
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
C*****************************     SCAL1
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
C*****************************     SCAL2
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
C*****************************     SCAL3
          FUNCTION SCAL3(Z)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          IF(Z.LT.2.D0) THEN
          IZ=IDINT(Z/2.D0)+3
          SRV=SPLT(Z,IZ)
          SCAL3=(3.D0*SRV/2.D0)**3
          ELSE 
          SCAL3=0.D0
          ENDIF
          RETURN
          END  
C*****************************     SPARSE
          SUBROUTINE SPARSE(XCOVT,X,Y,IPX,IPY,LPX,LPY,IWD,IPUNT,
     .JCRONA,JCOVAS,IBANDA,JBANDA,NUMERO,INDIR,LC,JC,  
     .NMAX,NWDMAX,N,NWD,IORDI,ILI) 
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
          CALL WINDOW(XC,YC,AMPX,AMPY,X,Y,IPX,IPY,LPX,LPY,IWD, 
     .NMAX,NWDMAX,N,NWD,ILI)
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
C*****************************     SPARS3
          SUBROUTINE SPARS3(XCOVT,X,Y,Z,IPX,IPY,IPZ,LPX,LPY,LPZ,
     .IWD,IPUNT,JCRONA,JCOVAS,IBANDA,JBANDA,NUMERO,INDIR,LC,JC,  
     .NMAX,NWDMAX,N,NWD,IORDI,ILI) 
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION X(1),Y(1),Z(1),IPX(1),IPY(1),IPZ(1), 
     .LPX(1),LPY(1),LPZ(1),IWD(1), 
     .IPUNT(1),JCRONA(1),JCOVAS(1),IBANDA(1),JBANDA(1),
     .NUMERO(1),INDIR(1),LC(1),JC(1)
          CALL TIORDA(X,Y,Z,IPX,IPY,IPZ,LPX,LPY,LPZ,N)
          XCOVT2=2.D0*XCOVT
          AMPX=4.D0*XCOVT  
          AMPY=4.D0*XCOVT  
          AMPZ=4.D0*XCOVT  
          JCRONA(1)=1  
          K=0  
          DO 20 I=1,N  
          XC=X(I)  
          YC=Y(I)  
          ZC=Z(I)  
          CALL POCKET(XC,YC,ZC,AMPX,AMPY,AMPZ,X,Y,Z,
     .IPX,IPY,IPZ,LPX,LPY,LPZ,IWD, 
     .NMAX,NWDMAX,N,NWD,ILI)
          CALL ORDINI(IWD,IPUNT,NWD,1)
          DO 10 J=1,NWD
          JJ=IPUNT(J)  
          II=IWD(JJ)
          DIST=DSQRT((XC-X(II))**2+(YC-Y(II))**2+(ZC-Z(II))**2)
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
C*****************************     SPLINP
          SUBROUTINE SPLINP(A,IA,P,B,KA,JA,LA,IB,JB,
     .C,T,JC,LC,D,W,R,U,ISRV,NNG,RES,FORMA,STP,
     .MMAX,NMAX,SCX,SCY,SCU,KK,LL,TOLL,
     .IAUTO,ITABLE,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,AST,PTO,ASTI
          CHARACTER*4 FORMA     
          CHARACTER*120 CDATO
          DIMENSION A(MMAX,1),IA(MMAX,1),KA(1),B(1),
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
          READ(CDATO,FORMA) NPT,X,Y,US
          IF(NPT.EQ.0) GOTO 15 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          EU=EU+US
          X0=DMIN1(X0,X)
          X9=DMAX1(X9,X)
          Y0=DMIN1(Y0,Y)
          Y9=DMAX1(Y9,Y)
          I=I+1
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
          WRITE(6,1170) NCX,NCY,NTC,NTI,STP,X0,Y0 
          ELSE 
          WRITE(6,1171) NCX,NCY,NTC,NTI,STP,X0,Y0 
          ENDIF 
          REWIND 1 
          I=0  
40        READ(1,1010) (DATO(L),L=1,120)
          READ(CDATO,FORMA) NPT,X,Y,US
          IF(NPT.EQ.0) GOTO 60 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          I=I+1        
          US=US-EU
          P(I)=1.D0
          B(I)=-US*SCU 
          KA(I)=16 
          X=(X-X0)*SCX 
          Y=(Y-Y0)*SCY 
          NX=IDINT(X/STP)+1
          NY=IDINT(Y/STP)+1
          NC=(NY-1)*NCX+NX 
          NNG(NC)=NNG(NC)+1
          SRVX=(X-(NX-1)*STP)/STP
          SRVY=(Y-(NY-1)*STP)/STP
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
          ENDIF
          GOTO 40  
60        N=NTI
          IF(ILI.EQ.1) THEN
          WRITE(6,1110)
          ELSE
          WRITE(6,1111)
          ENDIF
          DO 65 K=1,NTC
          J=(K-1)/NCX  
          WRITE(6,1115) K,NNG(K),K+J,K+J+1,K+NIX+J,K+NIX+J+1
65        CONTINUE 
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
          NV=0 
          DO 80 I=1,NIX
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
80        CONTINUE 
          DO 90 I=NIX+1,N-2*NIX+1,NIX
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
90        CONTINUE 
          DO 100 I=2*NIX,N-NIX,NIX 
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
100       CONTINUE 
          DO 110 I=(NIY-1)*NIX+1,N 
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
110       CONTINUE 
          DO 120 I=2,NIY-1 
          DO 120 J=2,NIX-1
          N1=(I-2)*NCX+(J-2)+1
          N2=N1+1  
          N3=N1+NCX
          N4=N3+1  
          IF(NNG(N1).EQ.0.AND.NNG(N2).EQ.0.AND.
     .NNG(N3).EQ.0.AND.NNG(N4).EQ.0) THEN 
          NV=NV+1 
          JJ=(I-1)*NIX+(J-1)+1  
          IS=LC(JJ+1)-1 
          C(IS)=C(IS)+1.D10
          ENDIF
120       CONTINUE 
C
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
C*****************************     SPLINS
          SUBROUTINE SPLINS(A,IA,P,B,KA,JA,LA,IB,JB,
     .C,T,JC,LC,D,W,R,U,ISRV,NNG,RES,FORMA,STP,
     .MMAX,NMAX,SCX,SCY,SCZ,SCU,KK,LL,TOLL,
     .IAUTO,ITABLE,IFILE,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,AST,PTO,ASTI
          CHARACTER*4 FORMA
          CHARACTER*120 CDATO
          DIMENSION A(MMAX,1),IA(MMAX,1),KA(1),B(1),
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
          Z0= 1.D30
          Z9=-1.D30
          EU=0.D0
          I=0  
10        READ(1,1010) (DATO(L),L=1,120)
          READ(CDATO,FORMA) NPT,X,Y,Z,US
          IF(NPT.EQ.0) GOTO 15 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          EU=EU+US
          X0=DMIN1(X0,X)
          X9=DMAX1(X9,X)
          Y0=DMIN1(Y0,Y)
          Y9=DMAX1(Y9,Y)
          Z0=DMIN1(Z0,Z)
          Z9=DMAX1(Z9,Z)
          I=I+1
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
          VOL=(X9-X0)*(Y9-Y0)*(Z9-Z0)
          UNI=DEXP(DLOG(VOL*SCX*SCY*SCZ/M)/3.D0)
          STP=IAUTO*UNI
          ENDIF
          X0=X0-STP/(1000.D0*SCX)  
          Y0=Y0-STP/(1000.D0*SCY)         
          Z0=Z0-STP/(1000.D0*SCZ) 
          NCX=IDINT((X9-X0)*SCX/STP)+3
          NCY=IDINT((Y9-Y0)*SCY/STP)+3
          NCZ=IDINT((Z9-Z0)*SCZ/STP)+3
          X0=X0-STP/SCX
          Y0=Y0-STP/SCY
          Z0=Z0-STP/SCZ
          NTC=NCX*NCY*NCZ
          NIX=NCX+1
          NIY=NCY+1
          NIZ=NCZ+1
          NTI=NIX*NIY*NIZ
          IF(ILI.EQ.1) THEN
          WRITE(6,1170) NCX,NCY,NCZ,NTC,NTI,STP,X0,Y0,Z0
          ELSE 
          WRITE(6,1171) NCX,NCY,NCZ,NTC,NTI,STP,X0,Y0,Z0
          ENDIF 
          REWIND 1 
          I=0  
40        READ(1,1010) (DATO(L),L=1,120)
          READ(CDATO,FORMA) NPT,X,Y,Z,US
          IF(NPT.EQ.0) GOTO 60 
          IF(DATO(KK).EQ.PTO.AND.DATO(LL).NE.AST) THEN 
          I=I+1
          US=US-EU
          P(I)=1.D0
          B(I)=-US*SCU 
          KA(I)=64
          X=(X-X0)*SCX 
          Y=(Y-Y0)*SCY 
          Z=(Z-Z0)*SCZ 
          NX=IDINT(X/STP)+1
          NY=IDINT(Y/STP)+1
          NZ=IDINT(Z/STP)+1
          NC=(NZ-1)*NCY*NCX+(NY-1)*NCX+NX 
          NNG(NC)=NNG(NC)+1
          SRVX=(X-(NX-1)*STP)/STP
          SRVY=(Y-(NY-1)*STP)/STP
          SRVZ=(Z-(NZ-1)*STP)/STP
          NI=(NZ-2)*NIY*NIX+(NY-2)*NIX+NX-2
          DO 50 II=1,4
          AZ=DBLE(2-II)
          JZ=5-II
          VALZ=SPLT(SRVZ+AZ,JZ)
          DO 50 J=1,4  
          AY=DBLE(2-J) 
          JY=5-J
          VALY=SPLT(SRVY+AY,JY)
          DO 50 K=1,4  
          JS=(II-1)*16+(J-1)*4+K 
          AX=DBLE(2-K) 
          JX=5-K
          A(I,JS)=VALZ*VALY*SPLT(SRVX+AX,JX)
          IA(I,JS)=NI+(II-1)*NIY*NIX+(J-1)*NIX+K
50        CONTINUE 
          ENDIF
          GOTO 40  
60        N=NTI
          IF(ILI.EQ.1) THEN
          WRITE(6,1110)
          ELSE
          WRITE(6,1111)
          ENDIF
          DO 65 K=1,NTC
          J=(K-1)/(NCY*NCX) 
          WRITE(6,1115) K,NNG(K),K+J,K+J+1,K+NIX+J,K+NIX+J+1,
     .K+NIY*NIX+J,K+NIY*NIX+J+1,K+NIY*NIX+NIX+J,K+NIY*NIX+NIX+J+1          
65        CONTINUE 
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
          NV=0 
          DO 80 I=1,NIY*NIX
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
80        CONTINUE
          DO 105 II=2,NIZ-1
          INI=(II-1)*NIY*NIX
          IFI=II*NIY*NIX
          DO 85 I=INI+1,INI+NIX 
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
85        CONTINUE
          DO 90 I=INI+NIX+1,IFI-2*NIX+1,NIX
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
90        CONTINUE 
          DO 100 I=INI+2*NIX,IFI-NIX,NIX 
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
100       CONTINUE
          DO 105 I=IFI-NIX+1,IFI       
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
105       CONTINUE 
          DO 110 I=(NIZ-1)*NIY*NIX+1,N  
          NV=NV+1  
          IS=LC(I+1)-1 
          C(IS)=C(IS)+1.D10
110       CONTINUE
          DO 120 II=2,NIZ-1 
          DO 120 I=2,NIY-1 
          DO 120 J=2,NIX-1 
          N1=(II-2)*NCY*NCX+(I-2)*NCX+(J-2)+1
          N2=N1+1  
          N3=N1+NCX
          N4=N3+1
          N5=N1+NCY*NCX
          N6=N5+1
          N7=N5+NCX
          N8=N7+1  
          IF(NNG(N1).EQ.0.AND.NNG(N2).EQ.0.AND.
     .NNG(N3).EQ.0.AND.NNG(N4).EQ.0.AND.NNG(N5).EQ.0.AND.   
     .NNG(N6).EQ.0.AND.NNG(N7).EQ.0.AND.NNG(N8).EQ.0) THEN
          NV=NV+1
          JJ=(II-1)*NIY*NIX+(I-1)*NIX+(J-1)+1
          IS=LC(JJ+1)-1 
          C(IS)=C(IS)+1.D10
          ENDIF
120       CONTINUE
C    
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
          WRITE(3,1015) X0,Y0,Z0,STP,EU,N,NIX,NIY,NIZ
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
          READ(CDATO,FORMA) NPT,X,Y,Z,US
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
1015      FORMAT(5D20.12,4I5)  
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
     .5X,'NUMEROSITA',8X,'(VERTICI DEL CUBO)'/)
1111      FORMAT(///7X,'CLASS -',11X,'KNOTS CONNECTED THE CLASS'/
     .7X,'QUANTITY',10X,'(VERTICES OF THE CUBE)'/)
1115      FORMAT(5X,2I5,15X,4I5/30X,4I5)
1125      FORMAT(//5X,'NUMERO DI ITERAZIONI (METODO DEL GRADIENTE ',
     .'CONIUGATO)',I5//)  
1126     FORMAT(//5X,'ITERATION NUMBER (CONJUGATE GRADIENT METHOD)',
     .I5//)  
1130      FORMAT(///5X,'NUMERO DI CONDIZIONE', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
1131      FORMAT(///5X,'CONDITION NUMBER', 
     .1X,'( CHI = MAX(ABS(C) / MAX(ABS(C**-1)) )',2X,1PD10.1,5X,I5)
1170      FORMAT(///5X,'NUMERO CLASSI IN X, Y E Z ',3I5,
     .//5X,'NUMERO TOTALE CLASSI      ',10X,I5,
     .//5X,'NUMERO NODI               ',10X,I5,
     .//5X,'PASSO DELLA INTERPOLAZIONE',5X,F10.3, 
     .//5X,'COORDINATE ORIGINE        ',F15.3,5X,F15.3,5X,F15.6//) 
1171      FORMAT(///5X,'NUMBER OF CLASSES IN X, Y AND Z',3I5,
     .//5X,'NUMBER OF CLASSES              ',10X,I5,
     .//5X,'NUMBER OF KNOTS                ',10X,I5,
     .//5X,'STEP  ',30X,F10.3, 
     .//5X,'INITIAL COORDINATES            ',F15.3,5X,F15.3,5X,F15.6//) 
1190      FORMAT(//5X,'NUMERO DEI DATI',
     ./5X,'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'//)
1191      FORMAT(//5X,'NUMBER OF DATA',
     ./5X,'EXCEEDING THE PREVIOUS DIMENSION'//)
          END  
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
1131      FORMAT(///5X,'NUMBER OF CLASSES',10X,I5,
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
C*****************************     TESTA
          SUBROUTINE TESTA(ICASO,ITIPO,IQUO,IPTO,NPTO,ILI) 
          CHARACTER*1 DATO,PARA,PARC,BLANK 
          DIMENSION DATO(120),IPTO(1)  
          DATA PARA/'('/,PARC/')'/,BLANK/' '/  
          DO 5 I=1,120 
          DATO(I)=BLANK
5         CONTINUE 
          IF(ILI.EQ.1) THEN
          WRITE(6,1000)
          IF(ICASO.EQ.1) THEN  
          IF(ITIPO.EQ.1) WRITE(6,1010) 
          IF(ITIPO.EQ.2) WRITE(6,1020) 
          IF(ITIPO.EQ.3) WRITE(6,1030) 
          ENDIF
          IF(ICASO.EQ.2) THEN 
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
          IF(ICASO.EQ.3) THEN  
          IF(ITIPO.EQ.1) WRITE(6,1070) 
          IF(ITIPO.EQ.2) WRITE(6,1080) 
          IF(ITIPO.EQ.3) WRITE(6,1090) 
          ENDIF
          ELSE
          WRITE(6,2000)
          IF(ICASO.EQ.1) THEN  
          IF(ITIPO.EQ.1) WRITE(6,2010) 
          IF(ITIPO.EQ.2) WRITE(6,2020) 
          IF(ITIPO.EQ.3) WRITE(6,2030) 
          ENDIF
          IF(ICASO.EQ.2) THEN 
          IF(IQUO.EQ.0) THEN
          IF(ITIPO.EQ.1) WRITE(6,2040) 
          IF(ITIPO.EQ.2) WRITE(6,2050) 
          IF(ITIPO.EQ.3) WRITE(6,2060) 
          ELSE 
          IF(ITIPO.EQ.1) WRITE(6,2070) 
          IF(ITIPO.EQ.2) WRITE(6,2080) 
          IF(ITIPO.EQ.3) WRITE(6,2090) 
          ENDIF
          ENDIF
          IF(ICASO.EQ.3) THEN  
          IF(ITIPO.EQ.1) WRITE(6,2070) 
          IF(ITIPO.EQ.2) WRITE(6,2080) 
          IF(ITIPO.EQ.3) WRITE(6,2090) 
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
1000      FORMAT(///5X,'DATI DI INPUT'///5X,'LEGEND:'/) 
1010      FORMAT(5X,'(1) NOME PUNTO'/ 
     .5X,'(2) TEMPO'/ 
     .5X,'(3) COMPONENTE X')
1020      FORMAT(5X,'(1) NOME PUNTO'/ 
     .5X,'(2) TEMPO'/ 
     .5X,'(3) COMPONENTE X'/  
     .5X,'(4) COMPONENTE Y')
1030      FORMAT(5X,'(1) NOME PUNTO'/ 
     .5X,'(2) TEMPO'/ 
     .5X,'(3) COMPONENTE X'/  
     .5X,'(4) COMPONENTE Y'/  
     .5X,'(5) COMPONENTE Z')
1040      FORMAT(5X,'(1) NOME PUNTO'/ 
     .5X,'(2) COORDINATA PIANA ASCISSA'/  
     .5X,'(3) COORDINATA PIANA ORDINATA'/ 
     .5X,'(4) COMPONENTE X')
1050      FORMAT(5X,'(1) NOME PUNTO'/ 
     .5X,'(2) COORDINATA PIANA ASCISSA'/  
     .5X,'(3) COORDINATA PIANA ORDINATA'/ 
     .5X,'(4) COMPONENTE X'/  
     .5X,'(5) COMPONENTE Y')
1060      FORMAT(5X,'(1) NOME PUNTO'/ 
     .5X,'(2) COORDINATA PIANA ASCISSA'/  
     .5X,'(3) COORDINATA PIANA ORDINATA'/ 
     .5X,'(4) COMPONENTE X'/  
     .5X,'(5) COMPONENTE Y'/  
     .5X,'(6) COMPONENTE Z')
1070      FORMAT(5X,'(1) NOME PUNTO'/ 
     .5X,'(2) COORDINATA PIANA ASCISSA'/  
     .5X,'(3) COORDINATA PIANA ORDINATA'/ 
     .5X,'(4) QUOTA'/ 
     .5X,'(5) COMPONENTE X')
1080      FORMAT(5X,'(1) NOME PUNTO'/ 
     .5X,'(2) COORDINATA PIANA ASCISSA'/  
     .5X,'(3) COORDINATA PIANA ORDINATA'/ 
     .5X,'(4) QUOTA'/ 
     .5X,'(5) COMPONENTE X'/  
     .5X,'(6) COMPONENTE Y')
1090      FORMAT(5X,'(1) NOME PUNTO'/ 
     .5X,'(2) COORDINATA PIANA ASCISSA'/  
     .5X,'(3) COORDINATA PIANA ORDINATA'/ 
     .5X,'(4) QUOTA'/ 
     .5X,'(5) COMPONENTE X'/  
     .5X,'(6) COMPONENTE Y'/  
     .5X,'(7) COMPONENTE Z')
1100      FORMAT(//5X,120A1/)  
2000      FORMAT(///5X,'INPUT DATA'///5X,'LEGENDA:'/) 
2010      FORMAT(5X,'(1) NAME OF THE POINT'/ 
     .5X,'(2) TIME'/ 
     .5X,'(3) X COMPONENT')
2020      FORMAT(5X,'(1) NAME OF THE POINT'/ 
     .5X,'(2) TIME'/ 
     .5X,'(3) X COMPONENT',/  
     .5X,'(4) Y COMPONENT')
2030      FORMAT(5X,'(1) NAME OF THE POINT'/ 
     .5X,'(2) TIME'/ 
     .5X,'(3) X COMPONENT'/  
     .5X,'(4) Y COMPONENT'/  
     .5X,'(5) Z COMPONENT')
2040      FORMAT(5X,'(1) NAME OF THE POINT'/ 
     .5X,'(2) X COORDINATE'/  
     .5X,'(3) Y COORDINATE'/ 
     .5X,'(4) X COMPONENT')
2050      FORMAT(5X,'(1) NAME OF THE POINT'/ 
     .5X,'(2) X COORDINATE'/  
     .5X,'(3) Y COORDINATE'/ 
     .5X,'(4) X COMPONENT'/  
     .5X,'(5) Y COMPONENT')
2060      FORMAT(5X,'(1) NAME OF THE POINT'/ 
     .5X,'(2) X COORDINATE'/  
     .5X,'(3) Y COORDINATE'/ 
     .5X,'(4) X COMPONENT'/  
     .5X,'(5) Y COMPONENT'/  
     .5X,'(6) Z COMPONENT')
2070      FORMAT(5X,'(1) NAME OF THE POINT'/ 
     .5X,'(2) X COORDINATE'/  
     .5X,'(3) Y COORDINATE'/ 
     .5X,'(4) Z COORDINATE'/ 
     .5X,'(5) X COMPONENT')
2080      FORMAT(5X,'(1) NAME OF THE POINT'/ 
     .5X,'(2) X COORDINATE'/  
     .5X,'(3) Y COORDINATE'/ 
     .5X,'(4) Z COORDINATE'/ 
     .5X,'(5) X COMPONENT'/  
     .5X,'(6) Y COMPONENT')
2090      FORMAT(5X,'(1) NAME OF THE POINT'/ 
     .5X,'(2) X COORDINATE'/  
     .5X,'(3) Y COORDINATE'/ 
     .5X,'(4) Z COORDINATE'/ 
     .5X,'(5) X COMPONENT'/  
     .5X,'(6) Y COMPONENT'/  
     .5X,'(7) Z COMPONENT')
2100      FORMAT(//5X,120A1/)  
          END 
C*****************************     TRIMP
          SUBROUTINE TRIMP(X,Y,U,IPX,IPY,LPX,LPY,IWD,SNOISE,OLATO, 
     .    IND,N,NN,NM,NWDMAX,NMAX,IAST,IU1,IU2,TOLL,ILI)
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
          CALL WINDOW(XC,YC,OLATO,OLATO,X,Y,IPX,IPY,LPX,LPY,IWD,NMAX,  
     .    NWDMAX,N,NWD,ILI)
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
C*****************************     TRIMS
          SUBROUTINE TRIMS(X,Y,Z,U,IPX,IPY,IPZ,LPX,LPY,LPZ,IWD, 
     .    SNOISE,OLATO,IND,N,NN,NM,NWDMAX,NMAX,IAST,IU1,IU2,TOLL,ILI)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO 
          DIMENSION X(1),Y(1),Z(1),U(1),IPX(1),IPY(1),IPZ(1),
     .    LPX(1),LPY(1),LPZ(1),IWD(1),IAST(1),DATO(120) 
          REWIND IU1
          REWIND IU2
          NM=0 
          DO 20 I=1,NN 
          READ(IU1,1000) (DATO(K),K=1,120) 
          IF(LPX(I).EQ.0) GOTO 15  
          XC=X(I)  
          YC=Y(I)
          ZC=Z(I)  
          CALL POCKET(XC,YC,ZC,OLATO,OLATO,OLATO,X,Y,Z,IPX,IPY,IPZ, 
     .    LPX,LPY,LPZ,IWD,NMAX, NWDMAX,N,NWD,ILI)
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
C*****************************     TRIMT
           SUBROUTINE TRIMT(Z,U,IPZ,LPZ,SNOISE,OLATO,NMAX,IND,N,NN,NM,  
     .    IAST,IU1,IU2,TOLL,ILI)
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
          IF(LZI.LT.LZS)THEN
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
          ENDIF
15        WRITE(IU2,1000) (DATO(K),K=1,120)
20        CONTINUE 
          RETURN
1000      FORMAT(120A1)
          END  
C*****************************     ZERO
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
