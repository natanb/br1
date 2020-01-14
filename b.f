        dimension jc(15),lc(16),ipunt(15)
        data jc/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
C        data jc/1,1,2,2,3,1,3,4,5,6,5,7,6,7,8/
C       data jc/1,2,1,3,2,3,4,8*0/
        m=15
        n=8
        i=1
        do j=1,m
        if(jc(j).eq.i) then
                lc(i)=j
                i=i+1
        endif
        enddo
        n=i-1
        do i=n,1,-1
        lc(i+1)=lc(i)+1
        enddo
        
        do j=1,n+1
        print *,j,lc(j)
        enddo
        end

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
