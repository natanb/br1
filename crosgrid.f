        PARAMETER (MMAX=840,NMAX=1171,MAXCOR=150)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION Z1(MMAX,NMAX),Z2(MMAX,NMAX)
        DIMENSION COV(MAXCOR),SRV(MAXCOR),SRU(MAXCOR)
        DIMENSION NSRU(MAXCOR),NSRV(MAXCOR)
        CHARACTER*80 NFILE1,NFILE2
        ZMEDU=0.D0
        ZMEDV=0.D0
        ZU2=0.D0
        ZV2=0.D0
        ZUV=0.D0
        NZ1=0
        NZ2=0
        NZ12=0
        READ(*,'(A80)') NFILE1
        READ(*,'(A80)') NFILE2
        READ(*,*) ZMASK,STEP
        READ(*,*) NCOR,NDX,NDY
        OPEN(UNIT=11,FILE=NFILE1,STATUS='UNKNOWN')
        OPEN(UNIT=12,FILE=NFILE2,STATUS='UNKNOWN')
        
        DO J=1,NDY
        DO I=1,NDX
        K=0
        READ(11,*) ZZ1
        READ(12,*) ZZ2
        Z1(I,J)=ZZ1
        Z2(I,J)=ZZ2
        IF(ZZ1.NE.ZMASK.and.ZZ2.ne.zmask) THEN
c        IF(ZZ1.NE.ZMASK) THEN
        NZ1=NZ1+1
        ZMEU=ZMEU+ZZ1
        ZU2=ZU2+ZZ1**2
        K=K+1
        ENDIF
        IF(ZZ2.NE.ZMASK.and.zz1.ne.zmask) THEN
c        IF(ZZ2.NE.ZMASK) THEN
        NZ2=NZ2+1
        ZMEV=ZMEV+ZZ2
        ZV2=ZV2+ZZ2**2
        K=K+1
        ENDIF
        IF(K.EQ.2) THEN
        NZ12=NZ12+1
        ZUV =ZUV +ZZ1*ZZ2
        ENDIF
        ENDDO
        ENDDO
        CLOSE(UNIT=11)
        CLOSE(UNIT=12)

        ZMEU=ZMEU/NZ1
        ZU2=(ZU2/NZ1-ZMEU*ZMEU)*NZ1/(NZ1-1)
        ZMEV=ZMEV/NZ2
        ZV2=(ZV2/NZ2-ZMEV*ZMEV)*NZ2/(NZ2-1)
        
        DO I=1,NDX
        DO J=1,NDY
        IF(Z1(I,J).NE.ZMASK) Z1(I,J)=Z1(I,J)-ZMEU
        IF(Z2(I,J).NE.ZMASK) Z2(I,J)=Z2(I,J)-ZMEV
        ENDDO
        ENDDO


        CALL GRICOV(NDX,NDY,NCOR,MMAX,NSRU,NSRV,
     .              ZMASK,COV,SRU,SRV,Z1,Z2)
                        
        SUSV=SQRT(ZU2*ZV2)           
        COV0 = ZUV/NZ12 - ZMEU*ZMEV
        if(nz12.gt.1)        print *,nz12,cov0/susv
        if(nz12.le.1)        print *,nz12," 0"
        PRINT *,"     N = ",NZ12,     "   RHO  = ",COV0/SUSV
        PRINT *,"  M(U) = ",ZMEU,     "   M(V) = ",ZMEV
        PRINT *,"SQM(U) = ",SQRT(ZU2)," SQM(V) = ",SQRT(ZV2)
        PRINT *,"CROSCOVARIANCE FUNCTION"
        WRITE(*,'(5X,''0.000'',f15.6,f10.3)') COV0,COV0/SUSV
        DO I=1,NCOR
        WRITE(*,'(F10.3,F15.6,F10.3)') 
     .STEP*(I-0.D0),COV(I),COV(I)/SUSV
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

        SUBROUTINE GRICOV(NDX,NDY,NCOR,MMAX,NSRU,NSRV,
     .                    ZMASK,COV,SRU,SRV,Z1,Z2)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Z1(MMAX,*),Z2(MMAX,*),COV(*),
     .SRU(*),SRV(*),NSRU(*),NSRV(*)

        NZ12=0
        DO I = 1,NCOR
        COV(I) =0.D0
        SRU(I) =0.D0
        SRV(I) =0.D0
        NSRU(I)=0
        NSRV(I)=0
        ENDDO
       

        DO I = 1,NDX
        
        IL = I - MAX(I-NCOR,1)
        IR = MIN(I+NCOR,NDX) - I
        IBD = I - IL
        DO J = 1,NDY
        IF(Z1(I,J).NE.ZMASK .AND. Z2(I,J).NE.ZMASK) THEN
        NZ12=NZ12+1
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
        KCOR = INT(DIST-1.D-10) + 1
        IF(Z1(ISRV,JSRV).NE.ZMASK) THEN
        NSRU(KCOR)=NSRU(KCOR)+1
        SRU(KCOR)  = SRU(KCOR) + Z1(ISRV,JSRV)
        ENDIF
        IF(Z2(ISRV,JSRV).NE.ZMASK) THEN 
        SRV(KCOR)  = SRV(KCOR) + Z2(ISRV,JSRV)
        NSRV(KCOR) = NSRV(KCOR) + 1
        ENDIF
        ENDDO
        ENDDO

        SRU(1)  = SRU(1) - Z1(I,J)
        SRV(1)  = SRV(1) - Z2(I,J)
        NSRU(1) = NSRU(1) - 1
        NSRV(1) = NSRV(1) - 1


        DO N = 1, NCOR
        IF(NSRV(N).GT.0) COV(N) = COV(N)+Z1(I,J)*SRV(N)/NSRV(N)
        IF(NSRU(N).GT.0) COV(N) = COV(N)+Z2(I,J)*SRU(N)/NSRU(N)
        ENDDO

        ENDIF
        ENDDO
        ENDDO
        
        NN2=2*NZ12
        DO N = 1, NCOR
        COV(N)=COV(N)/NN2
        ENDDO
        RETURN
        END
