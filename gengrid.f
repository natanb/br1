        read(*,*) x0,y0,step
        read(*,*) nx,ny
        
        k=0
        do j=1,ny
        y=y0+(j-1)*step
        do i=1,nx
        k=k+1
        x=x0+(i-1)*step
        write(*,'(i6,1x,2f8.1)') k,x,y
        enddo
        enddo
        end
