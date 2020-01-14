          PARAMETER (MMAX=200000,NMAX=100) 
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          CHARACTER*1 DATO,RESI,RES2,RES3,PCHAR,PTO,BLANK  
          CHARACTER*4 TITOLO,FORMA,FORMA2,FORMA3
          CHARACTER*15 CRESI,CRES2,CRES3
          DIMENSION X(MMAX),Y(MMAX),Z(MMAX),U(MMAX),V(MMAX),W(MMAX),
     .    A(MMAX,NMAX),B(MMAX),P(MMAX),C(NMAX,NMAX),D(NMAX),S(NMAX),
     .    IPTO(6),IAST(3),DATO(120),PCHAR(120),
     .    TITOLO(20),FORMA(20),FORMA2(20),FORMA3(20),
     .    RESI(15),RES2(15),RES3(15)
          EQUIVALENCE (CRESI,RESI(1)),(CRES2,RES2(1)),(CRES3,RES3(1))
          DATA PTO/'.'/,BLANK/' '/ 
          DATA IPTO/6*120/,IAST/116,118,120/
C
          CALL LINGUA(ILI)
          IF(ILI.EQ.1) THEN
          CALL PAPOLI(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IQUO, 
     .    IGRA,NGRA,TOLL,SCX,SCY,SCZ,SCU,SCV,SCW,FORMA,FORMA2,FORMA3)
          ELSE        
          CALL PIPOLI(TITOLO,ICASO,ITIPO,ISCRI,IFILE,IQUO, 
     .    IGRA,NGRA,TOLL,SCX,SCY,SCZ,SCU,SCV,SCW,FORMA,FORMA2,FORMA3)
          ENDIF  
          WRITE(6,1020) TITOLO 
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
          IF(N.GT.MMAX) THEN
          IF(ILI.EQ.1) THEN
          WRITE(6,1050)
          ELSE
          WRITE(6,1051)
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
C  
          REWIND 1 
          READ(1,1000) (DATO(I),I=1,120)
          LAST=IPTO(ICASO+IQUO)
          DO 60 I=LAST+1,LAST+10
          IF(DATO(I).EQ.BLANK) GOTO 70
60        CONTINUE 
70        LAST=I-1 
          IF(IQUO.EQ.1) THEN
          DO 80 I=4,6  
          IPTO(I-1)=IPTO(I)
80        CONTINUE 
          ENDIF
C
          IF(ICASO.EQ.1) THEN  
                    K=IPTO(2)  
                    L=IAST(1)  
                    CALL POLINT(X,U,A,B,P,C,D,S,IGRA,NGRA,TOLL,
     .              MMAX,NMAX,SCZ,SCU,K,L,FORMA,IFILE,ILI) 
                    IF(ITIPO.GE.2) THEN
                    K=IPTO(3)  
                    L=IAST(2)  
                    CALL POLINT(X,V,A,B,P,C,D,S,IGRA,NGRA,TOLL,
     .              MMAX,NMAX,SCZ,SCV,K,L,FORMA2,IFILE,ILI)
                    ENDIF  
                    IF(ITIPO.EQ.3) THEN
                    K=IPTO(4)  
                    L=IAST(3)  
                    CALL POLINT(X,W,A,B,P,C,D,S,IGRA,NGRA,TOLL,
     .              MMAX,NMAX,SCZ,SCW,K,L,FORMA3,IFILE,ILI)
                    ENDIF  
          ENDIF
          IF(ICASO.EQ.2) THEN  
                    K=IPTO(3)  
                    L=IAST(1)  
                    CALL POLINP(X,Y,U,A,B,P,C,D,S,IGRA,NGRA,TOLL,  
     .              MMAX,NMAX,SCX,SCY,SCU,K,L,FORMA,IFILE,ILI) 
                    IF(ITIPO.GE.2) THEN
                    K=IPTO(4)  
                    L=IAST(2)  
                    CALL POLINP(X,Y,V,A,B,P,C,D,S,IGRA,NGRA,TOLL,  
     .              MMAX,NMAX,SCX,SCY,SCV,K,L,FORMA2,IFILE,ILI)
                    ENDIF  
                    IF(ITIPO.EQ.3) THEN
                    K=IPTO(5)  
                    L=IAST(3)  
                    CALL POLINP(X,Y,W,A,B,P,C,D,S,IGRA,NGRA,TOLL,  
     .              MMAX,NMAX,SCX,SCY,SCW,K,L,FORMA3,IFILE,ILI)
                    ENDIF  
          ENDIF
          IF(ICASO.EQ.3) THEN  
                    K=IPTO(4)  
                    L=IAST(1)  
                    CALL POLINS(X,Y,Z,U,A,B,P,C,D,S,IGRA,NGRA,TOLL,  
     .              MMAX,NMAX,SCX,SCY,SCZ,SCU,K,L,FORMA,IFILE,ILI) 
                    IF(ITIPO.GE.2) THEN
                    K=IPTO(5)  
                    L=IAST(2)  
                    CALL POLINS(X,Y,Z,V,A,B,P,C,D,S,IGRA,NGRA,TOLL,  
     .              MMAX,NMAX,SCX,SCY,SCZ,SCV,K,L,FORMA2,IFILE,ILI)
                    ENDIF  
                    IF(ITIPO.EQ.3) THEN
                    K=IPTO(6)  
                    L=IAST(3)  
                    CALL POLINS(X,Y,Z,W,A,B,P,C,D,S,IGRA,NGRA,TOLL,  
     .              MMAX,NMAX,SCX,SCY,SCZ,SCW,K,L,FORMA3,IFILE,ILI)
                    ENDIF  
          ENDIF
C
          IF(IFILE.EQ.1) THEN
          REWIND 1 
          KONTU=0  
          KONTV=0  
          KONTW=0  
          DO 100 J=1,N 
          READ(1,1000) (DATO(I),I=1,120)
          DO 90 I=1,15 
          RESI(I)=BLANK
          RES2(I)=BLANK
          RES3(I)=BLANK
90        CONTINUE 
          IF(DATO(IPTO(ICASO+1)).EQ.PTO) THEN  
          KONTU=KONTU+1
          WRITE(CRESI,1010) U(KONTU)
          ENDIF
          IF(DATO(IPTO(ICASO+2)).EQ.PTO) THEN  
          KONTV=KONTV+1
          WRITE(CRES2,1010) V(KONTV)
          ENDIF
          IF(DATO(IPTO(ICASO+3)).EQ.PTO) THEN  
          KONTW=KONTW+1
          WRITE(CRES3,1010) W(KONTW)
          ENDIF
          WRITE(2,1000) (DATO(I),I=1,LAST),RESI,RES2,RES3  
100       CONTINUE 
          WRITE(2,1030)
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
