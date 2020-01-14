          PARAMETER (MMAX=30000,NMAX=40000,NBMAX=16,MNBMAX=MMAX*NBMAX, 
     .NNBMAX=25*NMAX)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          CHARACTER*1 DATO,RESI,RES2,RES3,PCHAR,BLANK,PTO 
          CHARACTER*4 FORMA,FORMA2,FORMA3,TITOLO
          CHARACTER*15 CRESI,CRES2,CRES3
          DIMENSION A(MMAX,NBMAX),IA(MMAX,NBMAX),KA(MMAX),B(MMAX), 
     .P(MMAX),JA(MNBMAX),LA(NMAX),IB(NMAX),JB(NMAX),
     .C(NNBMAX),T(NNBMAX),JC(NNBMAX),LC(NMAX),D(MMAX), 
     .H(NMAX),R(NMAX),S(NMAX),ISRV(NMAX),NNG(NMAX),
     .U(MMAX),V(MMAX),W(MMAX),
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
                    CALL SPLINP(A,IA,P,B,KA,JA,LA,IB,JB, 
     .              C,T,JC,LC,D,H,R,S,ISRV,NNG,U,FORMA,STP,
     .              MMAX,NMAX,SCX,SCY,SCU,K,L,TOLL,IAUTO,ITABLE,IFILE,
     .              ILI) 
          IF(ITIPO.GE.2) THEN  
                    K=IPTO(4)
                    L=IAST(2)
                    CALL SPLINP(A,IA,P,B,KA,JA,LA,IB,JB, 
     .              C,T,JC,LC,D,H,R,S,ISRV,NNG,V,FORMA2,STP2,  
     .              MMAX,NMAX,SCX,SCY,SCV,K,L,TOLL,IAUTO,ITABLE,IFILE,
     .              ILI) 
          ENDIF
          IF(ITIPO.EQ.3) THEN  
                    K=IPTO(5)
                    L=IAST(3)
                    CALL SPLINP(A,IA,P,B,KA,JA,LA,IB,JB, 
     .              C,T,JC,LC,D,H,R,S,ISRV,NNG,W,FORMA3,STP3,  
     .              MMAX,NMAX,SCX,SCY,SCW,K,L,TOLL,IAUTO,ITABLE,IFILE,
     .              ILI) 
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
