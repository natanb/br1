C*****************************     PREBIL
          SUBROUTINE PREBIL(NS,XS,YS,X0,Y0,STP,COEF,NI,U)
          IMPLICIT DOUBLE PRECISION(A-H,O-Z)
          DIMENSION COEF(3,1),X0(1),Y0(1),STP(1),NI(1),U(1)
          DO 1 I=1,NS
          U(I)=0.D0
1         CONTINUE
          DO 20 I=1,NS
          X=XS-X0(I)
          Y=YS-Y0(I)
          ST=STP(I)
          NIX=NI(I)
          NX=IDINT(X/ST)+1
          NY=IDINT(Y/ST)+1
          NC=(NY-1)*(NIX-1)+NX
          SRX=(X-(NX-1)*ST)/ST
          SRY=(Y-(NY-1)*ST)/ST
          NIC=(NY-2)*NIX+NX-1
	  U(I)=COEF(I,NIC)+COEF(I,NIC+1)*SRX+COEF(I,NIC+NIX)*SRY+
     .COEF(I,NIC+NIX+1)*SRX*SRY
20        CONTINUE
          RETURN
          END
