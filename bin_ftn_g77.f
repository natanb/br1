        real*4 a
        character*4 c
        equivalence (a,c)
        open(unit=21,name='bin.bin',form='unformatted')
        
        iu=21 
        do i=1,25
        a=1./float(i)
        ia=fputc(21,c(1:1))
        ia=fputc(iu,c(2:2))
        ia=fputc(iu,c(3:3))
        ia=fputc(iu,c(4:4))
        enddo

        
c        write(21) a(1)
c        write(21) a(2)
c        write(21) a(3)
c        write(21) a(4)
c        write(21) a(5)
c       write(21) a(6)
        close(21)
        end
