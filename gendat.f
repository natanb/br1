	real*4 x,y,z
	k=0
	j=100000000
	do i=1,300
	x=int(ran(j)*100)+ran(j)
	y=int(ran(j)*100)+ran(j)
	z=int(ran(j)*100)+ran(j)
	u=int(ran(j)*10)+ran(j)
	v=int(ran(j)*10)+ran(j)
	w=int(ran(j)*10)+ran(j)
	k=k+1
	write(10,'(i5,6f10.3)') k,x,y,z,u,v,w
	enddo
	end
