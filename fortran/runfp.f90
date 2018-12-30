program runfp
use fancy_timestep
!use solvefp
implicit none
integer::nx,ny,i,j
real *8,allocatable::p(:,:),xvals(:),yvals(:)
real *8::hx,hy
complex *16::testarr(100)

nx = 200
ny = 200

do i = 1,100
    testarr(i) = exp(cos(2.0d0*i)*sin(1.0d0*i))*cmplx(1.0D0,0.0D0,kind=16)
enddo

write(*,*) abs(easy_ifft(easy_fft(testarr,100),100)-easy_fft(easy_fft(testarr,100),100))

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
