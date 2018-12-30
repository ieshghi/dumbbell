module solvefp
contains

    subroutine fullrelax(p,nx,ny)
        implicit none
        integer::nx,ny
        real *8::p(nx,ny),xvals(nx),yvals(ny),hx,hy,eps,error,p2(nx,ny),pp(nx,ny)
        
        call gengrid(nx,ny,hx,hy,xvals,yvals)
        eps = 1e-12
        error = 10
        error = sqrt(sum(p**2))

        open(1,file = 'fpres/xder.dat')
        open(2,file = 'fpres/yder.dat')

        do while(error>eps)
            p2 = p
            write(2,*) p(nx/2,:)
            write(1,*) p(:,ny/2)
            
            pp = eulerstep(p,xvals,yvals)
            call relaxstep_bashforth(p,pp,xvals,yvals)

            error = sqrt(sum((p2-p)**2))
            write(*,*) 'error = ',error,'norm = ',sum(p2-p)
        enddo
        close(1)
        close(2)
        p = p2        
    endsubroutine fullrelax

    function consts(m)
       implicit none
       integer::m
       real *8::consts,lsmall,lbig,k,tau,tbar,gm,umax,kmult,vals(8)

       lsmall = 4! 
       lbig = 6  !These two set the relative width of the two parts of the potential. 
       k = 0.1 !Stiffness of spring
       tau = 0 !dimensionless temperature difference (dT/Tbar)
       tbar = 1 !average temperature
       umax = 1 !potential height
       gm = 10 !friction coefficient. Essentially acts as 1/dt for timestepping
       kmult = 5 !Since we can't sample the full spring stretching direction, we need to cut it off somewhere
       !This sets the cutoff at some multiple of the thermal length of the spring

       vals = (/lsmall,lbig,k,tau,tbar,umax,gm,kmult/)
       consts = vals(m)

    endfunction consts 
           
    subroutine relaxstep_bashforth(p,pp,xvals,yvals)
        implicit none
        integer::nx,ny
        real *8::p(:,:),pp(:,:),xvals(:),yvals(:)
        real *8,allocatable::p2(:,:),ax(:,:)
        
        nx = size(xvals)
        ny = size(yvals)
        allocate(p2(nx+2,ny+2))
        p2 = p + 3/2*step(p,xvals,yvals)-1/2*step(pp,xvals,yvals)
        p = p + p2
    endsubroutine relaxstep_bashforth
    
    function eulerstep(p,xvals,yvals)
        implicit none
        integer::nx,ny
        real *8::p(:,:),xvals(:),yvals(:),hx,temp(2),hy
        real *8,allocatable::eulerstep(:,:)
        nx = size(xvals)
        ny = size(yvals)
        allocate(eulerstep(nx+2,ny+2))
        eulerstep = p + step(p,xvals,yvals)
    endfunction eulerstep

    function step(p,xvals,yvals)
        implicit none
        real *8::p(:,:),xvals(:),yvals(:),hx,temp(2),hy
        real *8,allocatable::step(:,:),ax(:,:)
        integer::nx,ny,i,j,m

        nx = size(p(:,1))-2
        ny = size(p(1,:))-2
        hx = xvals(2)-xvals(1)
        hy = yvals(2)-yvals(1)
        allocate(step(nx+2,ny+2),ax(nx,ny))
        
        do i = 1,nx
            do j = 1,ny
                temp = avec(xvals(i),yvals(j))
                ax(i,j) = temp(1)
            enddo
        enddo

        call update_ghosts(ax,p,hx,hy)
        step = term1(p,xvals,yvals)+term2(p,xvals,yvals)+term3(p,xvals,yvals)
    endfunction step

    subroutine update_ghosts(ax,p,hx,hy)
        implicit none
        real *8::ax(:,:),p(:,:),hx,hy,b(2,2)
        real *8,allocatable:: p_ng(:,:)
        integer::nx,ny,i,j

        b = barray()

        nx = size(p(:,1))-2
        ny = size(p(1,:))-2
        
        allocate(p_ng(nx,ny))
        p_ng = p(2:nx+1,2:ny+1)

        do i = 1,ny
             !p(1,i+1) = -2*ax(1,i)*hx*p(2,i+1)/b(1,1)+p(3,i+1)+b(1,2)/b(1,1)*hx/hy*(p(2,i+2)-p(2,i))
             !p(nx+2,i+1) = 2*ax(nx,i)*hx*p(nx+1,i+1)/b(1,1)-p(nx,i+1)-b(1,2)/b(1,1)*hx/hy*(p(nx+1,i+2)-p(nx+1,i))

        enddo
        do i = 2,nx+1
            p(i,1) = p(i,ny+1)
            p(i,ny+2) = p(i,2)
        enddo
    endsubroutine update_ghosts
    
    function term3(p,xvals,yvals)
        implicit none
        real *8::p(:,:),xvals(:),yvals(:),b(2,2),hx,hy,ymax
        integer::nx,ny,i,j
        real *8,allocatable::term3(:,:),px(:,:),py(:,:),pxx(:,:),pyy(:,:),pxy(:,:),pxy2(:,:)

        ymax = consts(1)+consts(2)
        nx = size(p(:,1))-2
        ny = size(p(1,:))-2
        hx = xvals(2)-xvals(1)
        hy = yvals(2)-yvals(1)
        allocate(term3(nx+2,ny+2),px(nx+2,ny+2),pyy(nx+2,ny+2),pxx(nx+2,ny+2),pxy(nx+2,ny+2))
        px(:,:) = 0
        pxy(:,:) = 0
        pxx(:,:) = 0
        pyy(:,:) = 0
        term3(:,:) = 0
        b = barray()
        do i = 1,nx
            !pyy(i+1,:) = comp2der(p(i+1,:),ny)/(hy*hy)
            pyy(i+1,2:ny+1) = spec2der(yvals(1),ymax,ny,p(i+1,2:ny+1))
        enddo
        do i = 1,ny
            pxx(:,i+1) = comp2der(p(:,i+1),nx)/(hx*hx)
        enddo
        do i = 1,ny
            px(:,i+1) = compder(p(:,i+1),nx)/hx
        enddo
        do i = 1,nx
            !pxy(i+1,:) = compder(p(i+1,:),ny)/hy
            pxy(i+1,2:ny+1) = specder(yvals(1),ymax,ny,px(i+1,2:ny+1))
        enddo

        term3 = b(1,1)*pxx+b(2,2)*pyy+2*b(1,2)*pxy

    endfunction term3

    function term2(p,xvals,yvals)!term 2 looks like a . grad(p)
        implicit none
        real *8::p(:,:),xvals(:),yvals(:),temp(2),hx,hy,ymax
        integer::nx,ny,i,j
        real *8,allocatable::term2(:,:),ax(:,:),ay(:,:),px(:,:),py(:,:)

        ymax = consts(1)+consts(2)
        nx = size(p(:,1))-2
        ny = size(p(1,:))-2
        hx = xvals(2)-xvals(1)
        hy = yvals(2)-yvals(1)

        allocate(term2(nx+2,ny+2),ax(nx,ny),ay(nx,ny),px(nx+2,ny+2),py(nx+2,ny+2))
        px(:,:) = 0
        py(:,:) = 0
        term2(:,:) = 0

        do i=1,nx
            do j=1,ny
                temp = avec(xvals(i),yvals(j))
                ax(i,j) = temp(1)
                ay(i,j) = temp(2)
            enddo
        enddo
        
        do i = 1,nx
            py(i+1,:) = compder(p(i+1,:),ny)/hy
            py(i+1,2:ny+1) = specder(yvals(1),ymax,ny,p(i+1,2:ny+1))
        enddo
        do i = 1,ny
            px(:,i+1) = compder(p(:,i+1),nx)/hx
        enddo
        term2(2:nx+1,2:ny+1) = ax*px(2:nx+1,2:ny+1)+ay*py(2:nx+1,2:ny+1)

    endfunction term2

    function term1(p,xvals,yvals)!term 1 looks like div(a)
        implicit none
        real *8::p(:,:),xvals(:),yvals(:)
        integer::nx,ny,i,j
        real *8,allocatable::term1(:,:)

        nx = size(p(:,1))-2
        ny = size(p(1,:))-2
    
        allocate(term1(nx+2,ny+2))

        term1(:,:) = 0

        do i = 1,nx
            do j = 1,ny
                term1(i+1,j+1) = p(i+1,j+1)*diva(xvals(i),yvals(j))
            enddo
        enddo
    endfunction term1
    

    subroutine gengrid(nx,ny,hx,hy,xvals,yvals)
       implicit none
       integer::nx,ny
       real *8::xmax,ymax,hx,hy,k,tbar,xspr,xvals(nx),yvals(ny)
       
       ymax = consts(1)+consts(2)
       k = consts(3)
       if (k==0) then
           k = 1
       endif
       tbar = consts(5)

       xspr = sqrt(2*tbar/k)
       xmax = consts(8)*xspr

       xvals = linspace(-xmax,xmax,nx)
       yvals = linspace2(0.0d0,ymax,ny) !y is a periodic coordinate so we don't
       !include the endpoint

       hx = xvals(2)-xvals(1) 
       hy = yvals(2)-yvals(1)

    endsubroutine gengrid
    
    function diva(x,y)
        implicit none
        real *8::x,y,k,diva
        k = consts(3)
        diva = -2/consts(7)*(k+5.0d0/4*(upp(y+x/2)+upp(y-x/2)))

    endfunction diva

    function upp(z)
        implicit none
        real *8::z,upp
        
        upp = 0

    endfunction upp

    function uprime(z)
        implicit none
        real *8::z,lsmall,lbig,l,umax,uprime
        lsmall = consts(1)
        lbig = consts(2)
        umax = consts(6)

        l = lsmall + lbig
        if (mod(z,l)<lsmall) then
            uprime = umax/lsmall
        else
            uprime = -umax/lbig
        endif
    endfunction uprime

    function avec(x,y)
      implicit none
      real *8:: x,y,k,up,um,gm,fx,fy,avec(2)
      k = consts(3)
      gm = consts(7)

      um = uprime(y-x/2)
      up = uprime(y+x/2)
      fx = -k*x-0.5*(up-um)
      fy = -up-um
      avec(1) = 2*fx/gm
      avec(2) = fy/(2*gm)
    endfunction avec

    function barray()
      implicit none
      real *8::tau,tbar,gm,barray(2,2)

      tau = consts(4)
      tbar = consts(5)
      gm = consts(7)

      barray(1,1) = 4
      barray(1,2) = tau
      barray(2,1) = tau
      barray(2,2) = 1
      barray = (-1)*tbar/(2*gm)*barray

    endfunction barray

    function potential(x,y)
        implicit none
        real *8::potential,x,y,umax,lsmall,lbig,zp,zm,l,k
        potential = 0
        lsmall = consts(1)
        lbig = consts(2)
        l = lsmall+lbig
        umax = consts(6)
        k = consts(3)

        zp = mod(y+x/2,l)
        zm = mod(y-x/2,l)

        if (zp<lsmall) then
            potential = potential + umax*zp/lsmall
        else
            potential = potential - umax*(zp-l)/lbig
        endif
        if (zm<lsmall) then
            potential = potential + umax*zm/lsmall
        else
            potential = potential - umax*(zm-l)/lbig
        endif
        potential = potential + k*x*x/2
    endfunction potential

    function comp2der(arr,n)!second order finite difference first derivative,
      !assumes ghost cells and doesn't touch them
      implicit none
      real *8::dy,arr(n),comp2der(n)
      integer::n,i
      comp2der(1) = arr(1)
      comp2der(n) = arr(n)

      do i = 2,n-1
        comp2der(i) = (arr(i+1)-2*arr(i)+arr(i-1))
      enddo

    endfunction comp2der

    function compder(arr,n) !second order finite difference first derivative,
      !assumes ghost cells and doesn't touch them
      implicit none
      real *8::dy,arr(n),compder(n)
      integer::n,i
        
      compder(1) = arr(1) 
      compder(n) = arr(n)
      do i = 2,n-1
        compder(i) = (arr(i+1)-arr(i-1))/2
      enddo

    endfunction compder
    
    function spec2der(xmin,xmax,n,input)  !takes spectral derivatives of order n of a function evaluated at n points.
    use, intrinsic :: iso_c_binding
    implicit none
    include '/usr/local/include/fftw3.f03'
    type(c_ptr) :: plan
    integer :: n,i,ier,lensav
    integer *8::plan_forward,plan_backward
    complex(c_double_complex), dimension(n) :: cinput,output,input2,output2
    real *8,dimension(n)::x,k,der,input,spec2der
    real *8,parameter::pi =3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862
    real *8::xmin,xmax,dx,wsav(2*n+15)

    cinput = cmplx(input,0.0D0,kind=16)

    call dfftw_plan_dft_1d_(plan_forward,n, cinput,output,fftw_forward,fftw_estimate)
    call dfftw_execute_(plan_forward)
    call dfftw_destroy_plan_(plan_forward)

    k = fftfreq(n)
    do i=1,n
      input2(i) =-(2.0d0*pi/(xmax-xmin)*k(i))**2*output(i)/(n)
    end do

    call dfftw_plan_dft_1d_(plan_backward,n, input2, output2,fftw_backward,fftw_estimate)
    call dfftw_execute_(plan_backward, input2, output2)
    call dfftw_destroy_plan_(plan_backward)
    
    do i = 1,n
      spec2der(i) = real(output2(i),kind=16)
    end do
    end function spec2der

    function specder(xmin,xmax,n,input)  !takes spectral derivatives of order n of a function evaluated at n points.
    use, intrinsic :: iso_c_binding
    implicit none
    include '/usr/local/include/fftw3.f03'
    type(c_ptr) :: plan
    integer :: n,i
    integer *8::plan_forward,plan_backward
    complex(c_double_complex), dimension(n) :: cinput,output,input2,output2
    real *8,dimension(n)::x,k,der,input,specder
    real *8,parameter::pi =3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862
    real *8::xmin,xmax,dx

    cinput = cmplx(input,0.0D0,kind=16)

    call dfftw_plan_dft_1d_(plan_forward,n, cinput,output,fftw_forward,fftw_estimate)
    call dfftw_execute_(plan_forward)
    call dfftw_destroy_plan_(plan_forward)

    k = fftfreq(n)
    do i=1,n
      input2(i) =2.0d0*pi/(xmax-xmin)*k(i)*cmplx(0.0D0,1.0D0,kind=16)*output(i)/n
    end do

    call dfftw_plan_dft_1d_(plan_backward,n, input2, output2,fftw_backward,fftw_estimate)
    call dfftw_execute_(plan_backward, input2, output2)
    call dfftw_destroy_plan_(plan_backward)

    do i = 1,n
      specder(i) = real(output2(i),kind=16)
    end do

    end function specder

    function easy_fft(cinput,n)
        use,intrinsic :: iso_c_binding
        implicit none
        include '/usr/local/include/fftw3.f03'
        type(c_ptr) :: plan
        integer:: n
        integer *8::forward
        complex(c_double_complex),dimension(n)::cinput,easy_fft
            
        call dfftw_plan_dft_1d_(forward,n, cinput,easy_fft,fftw_forward,fftw_estimate)
        call dfftw_execute_(forward)
        call dfftw_destroy_plan_(forward)
    endfunction easy_fft

    function fftfreq(n)
        implicit none
        integer::n,i
        real *8::fftfreq(n)

        if (mod(n,2)==0) then
            do i = 1,n/2+1
                fftfreq(i) = i-1
            enddo
            do i = n/2+2,n
                fftfreq(i) = i-n-1
            enddo
        else
            do i = 1,(n+1)/2
                fftfreq(i) = i-1
            enddo
            do i = (n+3)/2,n
                fftfreq(i) = i-n-1
            enddo
        endif
    endfunction fftfreq

    function  linspace(a,b,n) !equivalent of python linspace, includes endpoint
        implicit none
        real *8, intent(in):: a,b !start and endpoint
        integer, intent(in):: n !number of elements
        integer:: i ! loop variable
        real *8:: dx, linspace(n)
        dx = (b-a)/(n-1) !spacing between x's
        do i = 1,n
            linspace(i)=a+(i-1)*dx ! fill the output array
        end do
    endfunction linspace

    function  linspace2(a,b,n) !equivalent of python linspace, except it doesn't include the endpoint
        implicit none
        real *8, intent(in):: a,b !start and endpoint
        integer, intent(in):: n !number of elements
        integer:: i ! loop variable
        real *8:: dx, linspace2(n)
        dx = (b-a)/(n) !spacing between x's
        do i = 1,n
            linspace2(i)=a+(i-1)*dx ! fill the output array
        end do
    endfunction linspace2

endmodule solvefp
