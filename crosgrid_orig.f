        PARAMETER (MMAX=100,NMAX=100,MAXCOR=150)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION Z1(MMAX,NMAX),Z2(MMAX,NMAX)
        DIMENSION COV(MAXCOR),SRV(MAXCOR),SRU(MAXCOR)
        DIMENSION NSRV(MAXCOR)
        CHARACTER*80 NFILE
        ZMEDU=0.D0
        ZMEDV=0.D0
        ZUV=0.D0
        READ(*,'(A80)') NFILE
        READ(*,*) NCOR,NDX,NDY
        OPEN(UNIT=11,FILE=NFILE,STATUS='UNKNOWN')
        
        DO J=1,NDY
        DO I=1,NDX
        READ(11,*) Z1(I,J),Z2(I,J)
        ZMEU=ZMEU+Z1(I,J)
        ZMEV=ZMEV+Z2(I,J)
        ZUV =ZUV +Z1(I,J)*Z2(I,J)
        ENDDO
        ENDDO
        CLOSE(UNIT=11)

        ZMEU=ZMEU/(NDX*NDY)
        ZMEV=ZMEV/(NDX*NDY)
        PRINT *,ZMEU,ZMEV
        DO I=1,NDX
        DO J=1,NDY
        Z1(I,J)=Z1(I,J)-ZMEU
        Z2(I,J)=Z2(I,J)-ZMEV
        ENDDO
        ENDDO

        CALL GRICOV(NDX,NDY,NCOR,MMAX,NSRV,COV,SRU,SRV,Z1,Z2)
                        
        COV(1) = ZUV/(NDX*NDY) - ZMEDU*ZMEDV
        DO I=1,NCOR
        WRITE(*,*) I,COV(I)
        ENDDO
        STOP
        END
***************************************************************************
*Calcola la funzione empirca di crosscovarianza su di un grigliato a maglia*
*quadratra.
*
*Input:
*        Z1 o Z1   = Valori della I e II componente ai nodi della griglia
*        NDX       = Numero di nodi sull'asse x (1...ndx)
*        NDY       = Numero di nodi sull'asse y (1...ndy)
*        SRU       = Vettore di accumulo della I componente
*        SRV       = Vettore di accumulo della II componente
*        NSRV      = Vettore di accumulo dei conteggi  
*        NCOR      = Numero dei passi della funzione empirica di crosscovarianza
*Output:
*        COV       = Valori della funzione empirica di crosscovarianza
*Variabili di servizio:
*        IR = Massimo spostamento a destra   rispetto al punto i,j
*        IL = Massimo spostamento a sinistra rispetto al punto i,j
*        JB = Massimo spostamento a in basso rispetto al punto i,j
*        JT = Massimo spostamento a in alto  rispetto al punto i,j
*        IBD,JBD= Coordinate del nodo in basso a sinistra
*        ISRV,JSRV = Variabili di comodo
*        NN,KCOR   = Variabili di comodo
*        I,J,K,L   = Variabili di comodo        
***************************************************************************

        SUBROUTINE GRICOV(NDX,NDY,NCOR,MMAX,NSRV,COV,SRU,SRV,Z1,Z2)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Z1(MMAX,100),Z2(MMAX,100),COV(*),
     .SRU(*),SRV(*),NSRV(*)

        DO I = 1,NCOR
        COV(I) =0.D0
        SRU(I) =0.D0
        SRV(I) =0.D0
        NSRV(I)=0
        ENDDO
       

        DO I = 1,NDX
        DO J = 1,NDY
        print *,Z1(I,J),Z2(I,J)
        enddo
        enddo

        DO I = 1,NDX
        print *,' processing row = ' ,I
        IL = I - MAX(I-NCOR,1)
        IR = MIN(I+NCOR,NDX) - I
        IBD = I - IL
        DO J = 1,NDY
        JB = J - MAX(J-NCOR,1)
        JT = MIN(J+NCOR,NDY) - J
        JBD = J - JB

        DO KK = 1,NCOR
        SRU(KK) =0.D0
        SRV(KK) =0.D0
        NSRV(KK)=0
        ENDDO

        DO K = 1,IL+IR+1
        ISRV = IBD + K - 1
        DO L = 1,JB+JT+1
        JSRV = JBD + L - 1
        DIST = SQRT( DBLE(I-ISRV)**2 + DBLE(J-JSRV)**2 )
        KCOR = INT(DIST+1.D-10) + 1
        SRU(KCOR)  = SRU(KCOR) + Z1(ISRV,JSRV)
        SRV(KCOR)  = SRV(KCOR) + Z2(ISRV,JSRV)
        NSRV(KCOR) = NSRV(KCOR) + 1
        ENDDO
        ENDDO

C        PRINT *,'NSRV(1)',NSRV(1)
        SRU(1)  = SRU(1) - Z1(I,J)
        SRV(1)  = SRV(1) - Z2(I,J)
        NSRV(1) = NSRV(1) - 1

         IF(NSRV(1).LT.0) THEN
        PRINT *,'***ERRORE***'
         stop
C        PRINT *,'I,J',I,J
C        PRINT *,'IBD,JBD',IBD,JBD
C        PRINT *,'IL,IR,JB,JT',IL,IR,JB,JT
C        PRINT *,NSRV(1),NSRV(2),NSRV(3)        
         ENDIF


        DO N = 1, NCOR
        IF(NSRV(N).GT.0) THEN
        COV(N) = COV(N)+( Z1(I,J)*SRV(N)+Z2(I,J)*SRU(N) )/NSRV(N)
c        print *,I,N,SRU(N),SRV(N),Z1(i,j),z2(i,j)
        ENDIF
        ENDDO


        ENDDO
         ENDDO
        
        NN2=2*NDX*NDY
        DO N = 1, NCOR
        COV(N)=COV(N)/NN2
        ENDDO

        RETURN
        END
