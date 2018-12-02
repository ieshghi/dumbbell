program runfp
use solvefp
implicit none
integer::nx,ny,i,j
real *8,allocatable::p(:,:)

nx = 100
ny = 100
allocate(p(nx+2,ny+2))
call fullrelax(p,nx,ny)

open(1,file = 'fpres/p.dat')

do i = 1,nx
    write(1,*) p(i,:)
enddo

close(1)

endprogram runfp
