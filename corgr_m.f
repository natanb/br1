C	**********************************************************
C	* CORGR_M:                                               *
C	* COMPUTE THE EMPIRICAL AUTOCOVARIANCE FUNCTION OF       *
C     * A REGULAR GRID                                         * 
C     * !!!!!!!!!!!!!!!!!!!!!non toglie la media!!!!!!!!!!!!!!!*
C	**********************************************************

      PARAMETER (MMAX=1050,NMAX=1050,MAXCOR=600)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Z(MMAX,NMAX),COV(MAXCOR),SRV(MAXCOR),NSRV(MAXCOR)
      CHARACTER*80 NFILE
      READ(*,'(A80)') NFILE
      READ(*,*) ZMASK,STEP
      READ(*,*) NCOR,NDX,NDY
      OPEN(UNIT=11,FILE=NFILE,STATUS='UNKNOWN')
      
      JJ=0
      ZMED=0.D0
      ZS2=0.D0
      SRV1=0.D0
      
      DO I=1,NDX
      DO J=1,NDY
      READ(11,*) X,Y,ZZ
      Z(I,J)=ZZ
      
      IF(ZZ.NE.ZMASK) THEN 
      ZMED=ZMED+ZZ
      ZS2=ZS2+ZZ*ZZ
      JJ=JJ+1
      ENDIF      
      
      ENDDO
      ENDDO
      CLOSE(11)

      ZMED=ZMED/JJ
      ZS2=(ZS2/JJ-ZMED**2)*JJ/(JJ-1)
      
      DO I=1,NDX
      DO J=1,NDY
      IF(Z(I,J).NE.ZMASK) Z(I,J)=Z(I,J)-ZMED
      ENDDO
      ENDDO
      
      PRINT *,"ZMED = ",ZMED, " SQM = ", SQRT(ZS2)
      PRINT *,"N. PUNTI NON MASCHERATI = ",JJ

      CALL GRICOV(Z,NDX,NDY,NCOR,COV,ZMASK,SRV,NSRV,MMAX,JJ)
      
      WRITE(*,'(5X,''0.000'',F15.6,5X,''1.000'')')ZS2
      
      DO I=1,NCOR
      WRITE(*,'(F10.3,F15.6,F10.3)') STEP*(I-0.5),COV(I),COV(I)/ZS2
      ENDDO
      END
***************************************************************************
*Calcola la funzione empirca di covarianza su di un grigliato a maglia    *
*quadratra.
*
*Input:
*	Z(Mmax,1) = Valori della componente ai nodi della griglia
*	NDX       = Numero di nodi sull'asse x (1...ndx)
*	NDY       = Numero di nodi sull'asse y (1...ndy)
*	SRV       = Vettore di accumulo delle componenti
*	NSRV      = Vettore di accumulo dei conteggi  
*	NCOR      = Numero dei passi della funzione empirica di covarianza
*Output:
*	COV       = Valori della funzione empirica di covarianza
*Variabili di servizio:
*	IR = Massimo spostamento a destra   rispetto al punto i,j
*	IL = Massimo spostamento a sinistra rispetto al punto i,j
*	JB = Massimo spostamento a in basso rispetto al punto i,j
*	JT = Massimo spostamento a in alto  rispetto al punto i,j
*	IBD,JBD= Coordinate del nodo in basso a sinistra
*	ISRV,JSRV = Variabili di comodo
*	NN,KCOR   = Variabili di comodo
*	I,J,K,L   = Variabili di comodo	
***************************************************************************

      SUBROUTINE GRICOV(Z,NDX,NDY,NCOR,COV,ZMASK,SRV,NSRV,MMAX,JJ)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(MMAX,*),COV(*),SRV(*),NSRV(*)

      DO I = 1,NCOR
      COV(I) =0.D0
      SRV(I) =0.D0
      NSRV(I)=0.D0
      ENDDO

      
      DO I = 1,NDX
C       print *,' processing col = ' ,I
      DO J = 1,NDY
      IF(Z(I,J).NE.ZMASK) THEN
      IL = I - MAX(I-NCOR,1)
      IR = MIN(I+NCOR,NDX) - I
      JB = J - MAX(J-NCOR,1)
      JT = MIN(J+NCOR,NDY) - J
      
      IBD = I - IL
      JBD = J - JB

      DO KK = 1,NCOR
      SRV(KK) =0.D0
      NSRV(KK)=0
      ENDDO

      DO K = 1,IL+IR+1
      ISRV = IBD + K - 1
      DO L = 1,JB+JT+1
      JSRV = JBD + L - 1
      IF (Z(ISRV, JSRV).NE.ZMASK) THEN
      DIST = SQRT( DBLE(I-ISRV)**2 + DBLE(J-JSRV)**2 )
      KCOR = INT(DIST+1.D-10) + 1
      SRV(KCOR) = SRV(KCOR) + Z(ISRV,JSRV)
          NSRV(KCOR)=NSRV(KCOR) + 1
      ENDIF
      ENDDO
      ENDDO

C      PRINT *,'NSRV(1)',NSRV(1)
      SRV(1) = SRV(1) - Z(I,J)
      NSRV(1)=NSRV(1) - 1

      IF(NSRV(1).LT.0) THEN
      PRINT *,'******'
      PRINT *,'I,J',I,J
      PRINT *,'IBD,JBD',IBD,JBD
      PRINT *,'IL,IR,JB,JT',IL,IR,JB,JT
      PRINT *,NSRV(1),NSRV(2),NSRV(3)      
      ENDIF


      DO N = 1, NCOR
      IF(NSRV(N).GT.0) COV(N)=COV(N) + Z(I,J)*SRV(N)/NSRV(N)
      ENDDO
      ENDIF
      ENDDO
      ENDDO

C	NN=NDX*NDY	
      DO N = 1, NCOR
      COV(N)=COV(N)/JJ
      ENDDO
      
C     PRINT *,"N. PUNTI NON MASCHERATI = ",JJ
      
      RETURN
      END
