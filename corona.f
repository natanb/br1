          PARAMETER (NMAX=700000,MMAX=1000,NWDMAX=100000)  
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)   
          CHARACTER*1 DATO,PCHAR,PTO,BLANK 
          CHARACTER*4 FORMA,TITOLO 
          DIMENSION X(NMAX),Y(NMAX),Z(NMAX),
     .    U(NMAX),V(NMAX),W(NMAX),
     .    IPXU(NMAX),IPYU(NMAX),IPZU(NMAX),
     .    LPXU(NMAX),LPYU(NMAX),LPZU(NMAX),
     .    IPXV(NMAX),IPYV(NMAX),IPZV(NMAX),
     .    LPXV(NMAX),LPYV(NMAX),LPZV(NMAX),
     .    IPXW(NMAX),IPYW(NMAX),IPZW(NMAX),
     .    LPXW(NMAX),LPYW(NMAX),LPZW(NMAX),
     .    COVU(MMAX),COVV(MMAX),COVW(MMAX),
     .    SRV(MMAX),NSRV(MMAX),
     .    IWD(NWDMAX),IPTO(6),IAST(3),DATO(120),PCHAR(120),
     .    TITOLO(20),FORMA(20) 
          DATA PTO/'.'/,BLANK/' '/ 
          DATA IPTO/6*120/,IAST/116,118,120/
C
          OPEN(3)
          OPEN(4)
          CALL LINGUA(ILI)
          IF(ILI.EQ.1) THEN
          CALL PACORO(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IAUTO,ITRIM,
     .IQUO,TOLL,AINF,ASUP,STP,SCX,SCY,SCZ,SCU,SCV,SCW,FORMA)
          ELSE 
          CALL PICORO(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IAUTO,ITRIM,
     .IQUO,TOLL,AINF,ASUP,STP,SCX,SCY,SCZ,SCU,SCV,SCW,FORMA)
          ENDIF 
          WRITE(6,1020) TITOLO 
          DO 5 I=1,120 
          PCHAR(I)=BLANK
5         CONTINUE 
          J=0  
10        J=J+1
          READ(1,1000) (DATO(I),I=1,120)
          WRITE(3,1000) (DATO(I),I=1,120)
          WRITE(4,1000) (DATO(I),I=1,120)
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
          IF(N.GT.NMAX) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,1060)
          ELSE
          WRITE(6,1061)
          ENDIF
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
          IF(ILI.EQ.1) THEN  
          WRITE(6,1050)
          ELSE
          WRITE(6,1051)
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
C  
          IF(IQUO.EQ.1) THEN
          DO 60 I=4,6  
          IPTO(I-1)=IPTO(I)
60        CONTINUE 
          ENDIF
          WRITE(2,1030)
C          
          IF(ICASO.EQ.1) THEN  
                    IF(ITIPO.EQ.1) CALL OTTIT1 
     .              (X,W,IPXW,LPXW,  
     .              COVW,SRV,NSRV,AINF,ASUP,STP,
     .              NMAX,SCZ,SCW,IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,
     .              ILI)
                    IF(ITIPO.EQ.2) CALL OTTIT2 
     .              (X,U,V,IPXU,IPXV,LPXU,LPXV,
     .              COVU,COVV,SRV,NRSV,AINF,ASUP,STP,  
     .              NMAX,SCZ,SCU,SCV,IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,
     .              ILI)
                    IF(ITIPO.EQ.3) CALL OTTIT3 
     .              (X,U,V,W,IPXU,IPXV,IPXW,LPXU,LPXV,LPXW,
     .              COVU,COVV,COVW,SRV,NSRV,AINF,ASUP,  
     .              STP,NMAX,SCZ,SCU,SCV,SCW,IPTO,IAST,FORMA,
     .              IAUTO,ITRIM,TOLL,ILI)
          ENDIF
          IF(ICASO.EQ.2) THEN  
                    IF(ITIPO.EQ.1) CALL OTTIP1 
     .              (X,Y,W,IPXW,IPYW,LPXW,LPYW,
     .              COVW,SRV,NSRV,AINF,ASUP,STP,IWD,  
     .              NMAX,NWDMAX,SCX,SCY,SCW,IPTO,IAST,FORMA,
     .              IAUTO,ITRIM,TOLL,ILI)
                    IF(ITIPO.EQ.2) CALL OTTIP2 
     .              (X,Y,U,V,IPXU,IPYU,IPXV,IPYV,  
     .              LPXU,LPYU,LPXV,LPYV, 
     .              COVU,COVV,SRV,NSRV,AINF,ASUP,STP,
     .              IWD,NMAX,NWDMAX,SCX,SCY,SCU,SCV,
     .              IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,ILI)
                    IF(ITIPO.EQ.3) CALL OTTIP3 
     .              (X,Y,U,V,W,IPXU,IPYU,IPXV,IPYV,IPXW,IPYW,  
     .              LPXU,LPYU,LPXV,LPYV,LPXW,LPYW,
     .              COVU,COVV,COVWW,SRV,NSRV,AINF,ASUP,STP,
     .              IWD,NMAX,NWDMAX,SCX,SCY,SCU,SCV,SCW,
     .              IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,ILI)
          ENDIF
          IF(ICASO.EQ.3) THEN  
                    IF(ITIPO.EQ.1) CALL OTTIS1 
     .              (X,Y,Z,W,IPXW,IPYW,IPZW,
     .              LPXW,LPYW,LPZW,
     .              COVW,SRV,NSRV,AINF,ASUP,STP,  
     .              IWD,NMAX,NWDMAX,SCX,SCY,SCZ,SCW,
     .              IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,ILI)
                    IF(ITIPO.EQ.2) CALL OTTIS2 
     .              (X,Y,Z,U,V,IPXU,IPYU,IPZU,IPXV,IPYV,IPZV,  
     .              LPXU,LPYU,LPZU,LPXV,LPYV,LPZV,
     .              COVU,COVV,SRV,NSRV,AINF,ASUP,STP, 
     .              IWD,NMAX,NWDMAX,SCX,SCY,SCZ,SCU,SCV,
     .              IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,ILI)
                    IF(ITIPO.EQ.3) CALL OTTIS3 
     .              (X,Y,Z,U,V,W,IPXU,IPYU,IPZU,IPXV,IPYV,IPZV,
     .              IPXW,IPYW,IPZW,LPXU,LPYU,LPZU,  
     .              LPXV,LPYV,LPZV,LPXW,LPYW,LPZW,
     .              COVU,COVV,COVW,SRV,NSRV,AINF,ASUP,STP,
     .              IWD,NMAX,NWDMAX,SCX,SCY,SCZ,SCU,SCV,SCW,
     .              IPTO,IAST,FORMA,IAUTO,ITRIM,TOLL,ILI)
          ENDIF
          IF(ITRIM.EQ.1) THEN 
          IF(ILI.EQ.1) THEN 
          ALFA=2.D0*(1.D0-DFN01(TOLL))
          WRITE(6,1040) ALFA
          ELSE
          ALFA=2.D0*(1.D0-DFN01(TOLL))
          WRITE(6,1041) ALFA
          ENDIF
          ENDIF
          IF(IFILE.EQ.1) THEN
          REWIND 4
          DO 80 J=1,N
          READ(4,1000) (DATO(I),I=1,120)
          WRITE(2,1000) (DATO(I),I=1,120)
80        CONTINUE
          WRITE(2,1030)
          ENDIF
          CLOSE(3)
          CLOSE(4)
          STOP 
1000      FORMAT(120A1)
1010      FORMAT(20A4) 
1020      FORMAT(//5X,20A4,//) 
1030      FORMAT(5X,120A1) 
1040      FORMAT(/5X,'LIVELLO DI SIGNIFICATIVITA DEL TEST '
     .'DI REIEZIONE DEGLI ERRORI GROSSOLANI: ALFA =',F5.3, 
     .' (SU DUE CODE)')
1041      FORMAT(/5X,'SIGNIFICANCE LEVEL OF THE TEST '
     .'FOR BLUNDER REJECTION: ALFA =',F5.3, 
     .' (ON TWO SIDES)')
1050      FORMAT(5X,'WRONG CASE AND/OR TYPE OF DATA'///)
1051     FORMAT(5X,'CASO E/O TIPO DEI DATI ERRATO'///)
1060      FORMAT(5X,'NUMERO DEI DATI ',/5X,
     .'ECCEDENTE IL DIMENSIONAMENTO PREVISTO'///)  
1061      FORMAT(5X,'NUMBER OF DATA ',/5X,
     .'EXCEEDING THE PREVIOUS DIMENSION'///) 
          END  
