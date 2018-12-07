program runfp
use solvefp
implicit none
integer::nx,ny,i,j
real *8,allocatable::p(:,:),xvals(:),yvals(:)
real *8::hx,hy

nx = 10000
ny = 100
allocate(p(nx+2,ny+2),xvals(nx),yvals(ny))
call gengrid(nx,ny,hx,hy,xvals,yvals)

do i = 1,nx
do j = 1,ny
    p(i,j) = exp(-potential(xvals(i),yvals(j))/consts(5))
enddo
enddo

p = p/sum(p)

call fullrelax(p,nx,ny)

open(1,file = 'fpres/p.dat')

do i = 1,nx
    write(1,*) p(i,:)
enddo

close(1)

endprogram runfp
