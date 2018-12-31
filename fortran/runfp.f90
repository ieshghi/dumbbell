program runfp
use fancy_timestep
!use solvefp
implicit none
integer::nx,ny,i,j
real *8,allocatable::p(:,:),xvals(:),yvals(:)
real *8::hx,hy
complex *16::testarr(10),test2(11)

nx = 200
ny = 200

nx = nx+2
ny = ny+2
allocate(p(nx,ny),xvals(nx),yvals(ny))
call gengrid(nx,ny,hx,hy,xvals,yvals)
p(:,:) = 0

do i = 5,nx-5
do j = 5,ny-5
    p(i,j) = exp(-0.1*xvals(i)**2-(yvals(j)-5)**2)
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
