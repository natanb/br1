      PROGRAM PREDIZIONE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAXS=3000)
      CHARACTER*80 FILEIN,FILEOUT,NFILE
      DIMENSION  COEF(NMAXS)


       PRINT *,' NOME DEL FILE DEI COEFFICIENTI '
       READ(*,'(A80)') NFILE
       OPEN(9,FILE=NFILE,STATUS='OLD')
       READ(9,*) XORI,YORI,STP,ES,NI,NIX 
       READ(9,*) (COEF(J),J=1,NI)  
       XMXS=XORI+(NIX-1)*STP
       YMXS=YORI+(NI/NIX - 1)*STP
       CLOSE(9)

	WRITE (6,*) ' SCRIVI IL FILE DI INPUT '
        READ(*,'(A80)') FILEIN
	WRITE (6,*) ' SCRIVI IL FILE DI OUTPUT '
        READ(*,'(A80)') FILEOUT   
 
      OPEN(UNIT=13, FILE=FILEIN , STATUS='OLD')
      OPEN(UNIT=14, FILE=FILEOUT , STATUS='NEW')

	
	I=0
	print *,XMXS,YMXS,XORI,YORI

20	READ(13,*,END=10) N,XS,YS
	IF(XS.LT.XMXS.AND.XS.GT.XORI.AND.YS.LT.YMXS.
     +     AND.YS.GT.YORI)THEN
	    CALL PREBIL(XS,YS,XORI,YORI,STP,COEF,NIX,ZS)
	I=I+1

	write(14,'(5X,3(F15.5))') XS,YS,ZS

	ENDIF
	GOTO 20

10	N1=I
	write(*,*) N1

	STOP
      END

C*****************************     PREBIL
          SUBROUTINE PREBIL(XS,YS,X0,Y0,STP,COEF,NI,U)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION COEF(1)
          
          U=0.D0
          X=XS-X0
          Y=YS-Y0
          ST=STP
          NIX=NI
          NX=IDINT(X/ST)+1
          NY=IDINT(Y/ST)+1
          NC=(NY-1)*(NIX-1)+NX
          SRX=(X-(NX-1)*ST)/ST
          SRY=(Y-(NY-1)*ST)/ST
          NIC=(NY-1)*NIX+NX
          U=COEF(NIC)*(1.D0-SRX-SRY+SRX*SRY)+
     .COEF(NIC+1)*(SRX-SRX*SRY)+COEF(NIC+NIX)*(SRY-SRX*SRY)+
     .COEF(NIC+NIX+1)*SRX*SRY
          RETURN
          END