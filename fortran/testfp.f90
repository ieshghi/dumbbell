program testfp
use solvefp
implicit none
integer::nx,ny,i,j,ns(6)
real *8,allocatable::x(:),y(:),calcder(:),tder(:)
real *8::hx,hy,pi,dx,err(6)

pi = 4*atan(1.0d0)
allocate(x(1),y(1),calcder(1),tder(1))
ns = (/10,100,1000,10000,100000,1000000/)
!Test 1: derivatives

do i = 1,6
    deallocate(x,y,calcder,tder)
    allocate(x(ns(i)),y(ns(i)),calcder(ns(i)),tder(ns(i)))
    x = linspace2(0.0d0,2*pi,ns(i))
    y = cos(x)*sin(x)
    dx = x(2)-x(1)
    calcder = spec2der(x(1),2*pi,ns(i),y) 
!    calcder = compder(y,ns(i))/(dx) 
    tder = -4*cos(x)*sin(x) 
    err(i) = sum((calcder(2:(ns(i)-1))-tder(2:(ns(i)-1)))**2)/(ns(i)-2)
enddo

open(1,file='err/err.dat')
    do i = 1,6
        write(1,*) err(i)
    enddo
close(1)

endprogram testfp
