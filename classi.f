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
