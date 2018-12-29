module fancy_timestep
use solvefp
contains






function der_oper(nx,k)
    implicit none
    complex *16::der_oper(nx+2,nx+2),im
    real *8::k,b(2,2)
    integer::nx,ny,i,j
    b = barray()/2
    
    im = complex(0.0D0,1.0D0,kind=16) !double precision imaginary unit

    der_oper(:,:) = 0
    do i = 2,nx+1
        der_oper(i,i-1) = 1/2*(b(1,1)-2*k*im*b(1,2))
        der_oper(i,i+1) = 1/2*(b(1,1)+2*k*im*b(1,2))
        der_oper(i,i) = 1-k**2*b(2,2)-b(1,1)
    enddo
    der_oper(1,1) = 1
    der_oper(nx+2,nx+2) = 1
endfunction der_oper

function gtilde(p,xvals,yvals)
    implicit none
    real *8::p(:,:),xvals(:),yvals(:),temp(2)
    real *8,allocatable,dimension(:,:)::ax,ay,step
    complex *16,allocatable::gtilde(:,:)
    integer::nx,ny,i

    nx = size(p(:,1))-2
    ny = size(p(1,:))-2
    allocate(gtilde(nx+2,ny+2),ax(nx,ny),ay(nx,ny),step(nx+2,ny+2))
        
    do i = 1,nx
        do j = 1,ny
           temp = avec(xvals(i),yvals(j))
           ax(i,j) = temp(1)
           ay(i,j) = temp(2)
        enddo
    enddo

    call update_ghosts(ax,p,hx,hy)
    step = term1(p,xvals,yvals)+term2(p,xvals,yvals)
    gtilde(:,:) = 0
    do i = 2,nx+1
        gtilde(i,2:ny+1) = easy_fft(step(i,2:ny+1),ny)
    enddo
endfunction gtilde

endmodule fancy_timestep
