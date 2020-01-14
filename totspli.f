          PARAMETER (MMAX=10000,NMAX=5000,NBMAX=16,MNBMAX=MMAX*NBMAX, 
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
C*****************************     IDIF2
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
C*****************************     LINGUA
          SUBROUTINE LINGUA(ILI)
10        WRITE(*,100)
          READ(*,*,ERR=10) ILI
          WRITE(*,'(///)')
          RETURN
100       FORMAT(///5X,'ITALIANO (ITALIAN): 1'/
     .    5X,'(INGLESE) ENGLISH:  0')
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
          WRITE(6,1170) NCX,NCY,NTC,NTI,STP,X0,Y0 
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
          FUNCTION DFN01(Z)
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          REAL*4 SRW,SRV,ERF  
          SRW=DABS(Z)/DSQRT(2.D0)
          SRV=DBLE(ERF(SRW))  
          DFN01=.5D0+.5D0*SRV*Z/DABS(Z)  
          RETURN
          END  
C*************************    ERF
      FUNCTION ERF(X)
      IF(X.LT.0.)THEN
        ERF=-GAMMP(.5,X**2)
      ELSE
        ERF=GAMMP(.5,X**2)
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
C*************************    BESSJ0
          FUNCTION DFJ0(X) 
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)
          REAL*4 SRW,SRV,BESSJ0
	  SRW=X
	  SRV=BESSJ0(SRW)  
          DFJ0=DBLE(SRV) 
          RETURN
          END  
C
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
C*************************    BESSJ1
          FUNCTION DFJ1(X) 
          IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
          REAL*4 SRV,SRW,BESSJ1  
          SRV=X
          SRW=BESSJ1(SRV)
          DFJ1=2.D0*DBLE(SRW)
          RETURN
	  END
C
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
C*************************    LUNG
          INTEGER FUNCTION LUNG(A) 
          CHARACTER*1 A
          DIMENSION A(1)
          DO 5 I=60,1,-1
          IF(A(I).NE.' ') GOTO 10  
5         CONTINUE 
10        LUNG=I
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
   50 ILI=NOME(N+1)
      IF(ILI.EQ.1) THEN
      WRITE(6,100)
      ELSE
      WRITE(6,101)
      ENDIF
      DO 60 I=1,N,20
      I1=I
      I2=MIN0(N,I+19)
      WRITE(6,120) (J,J=I1,I2)
      WRITE(6,120) (NOME(J),J=I1,I2)
      WRITE(6,120) (LIVEL(J),J=I1,I2)
      WRITE(6,110)
   60 CONTINUE
      STOP
  100 FORMAT(///5X,'GRAFO NON CONNESSO'///)
  101 FORMAT(///5X,'NO PATH IN THE GRAPH'///)
  110 FORMAT(/)
  120 FORMAT(5X,20I5)
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
      DATA IFORM/'(5X,','2I5,','4X, ',18*'    ','4X, ','3I5)'/
      INTERO='I5, '
      ICS='5X, '
      ILI=NOME(N+1)
      IF(ILI.EQ.1) THEN
      IF(INIZIO.EQ.0) WRITE(6,100)
      IF(INIZIO.EQ.1) WRITE(6,110)
      ELSE
      IF(INIZIO.EQ.0) WRITE(6,101)
      IF(INIZIO.EQ.1) WRITE(6,111)
      ENDIF
      DO 5 L=1,18
      ISRV(L)=0
    5 CONTINUE
      DO 30 I=1,N
      DO 10 J=1,18
      IFORM(J+3)=ICS
   10 CONTINUE
      KRONA=CRONA(I+1)-CRONA(I)
      INI=CRONA(I)+1-INIZIO
      IFI=CRONA(I+1)-INIZIO
      K=MIN0(KRONA,18)
      DO 20 J=1,K
      IFORM(J+3)=INTERO
   20 CONTINUE
      IFJ=MIN0(IFI,INI+17)
      LL=IFJ-INI+1
      DO 21 L=1,LL
      ISRV(L)=GRAFO(INI+L-1)
   21 CONTINUE
      WRITE(6,IFORM) I,NOME(I),(ISRV(J),J=1,LL),KRONA,BANDI(I),BANDA(I)
      J=0
   25 J=J+1
      DO 23 L=1,18
      ISRV(L)=0
   23 CONTINUE
      INJ=INI+18*J
      IF(INJ.GT.IFI) GO TO 30
      IFJ=MIN0(IFI,INJ+17)
      LL=IFJ-INJ+1
      DO 28 L=1,LL
      ISRV(L)=GRAFO(INJ+L-1)
   28 CONTINUE
      WRITE(6,120) (ISRV(K),K=1,LL)
      GO TO 25
   30 CONTINUE
      RETURN
  100 FORMAT(/////5X,'MATRICE DISEGNO SIMBOLICA',
     .//,5X,'(NUMERO DI ORDINE E NOME DEL PUNTO - MISURE CORONA,',
     .' CON TOTALE - BANDA INDIETRO E AVANTI)',/)
  101 FORMAT(/////5X,'SYMBOLICAL DESIGN MATRIX',
     .//,5X,'(ORDERING NUMBER AND NAME OF THE POINT - RING MEASURES,',
     .' WITH TOTAL - FORWARD AND BACKWARD BANDWIDTHS)',/)
  110 FORMAT(/////5X,'MATRICE NORMALE SIMBOLICA',
     .//,5X,'(NUMERO DI ORDINE E NOME DEL PUNTO - PUNTI CORONA,',
     .' CON TOTALE - BANDA INDIETRO E AVANTI)',/)
  111 FORMAT(/////5X,'SYMBOLICAL NORMAL MATRIX',
     .//,5X,'(ORDERING NUMBER AND NAME OF THE POINT - RING POINTS,',
     .' WITH TOTAL - FORWARD AND BACKWARD BANDWIDTHS)',/)
  120 FORMAT(19X,18I5)
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
