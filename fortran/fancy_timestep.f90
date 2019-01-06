module fancy_timestep
contains

subroutine fullrelax(p,nx,ny)
    implicit none
    integer::nx,ny,i
    real *8::p(nx,ny),xvals(nx),yvals(ny),hx,hy,eps,error,p2(nx,ny),pp(nx,ny),j(2,nx,ny),dj(nx,ny)
        
    call gengrid(nx,ny,hx,hy,xvals,yvals)
    eps = 1e-12
    error =10

    open(1,file = 'fpres/xder.dat')
    open(2,file = 'fpres/yder.dat')

    pp = timestep_eul(nx,ny,xvals,yvals,p)
    i = 0
    j = getcurrent(pp,xvals,yvals)
    dj = div(j(1,:,:),j(2,:,:),xvals,yvals)
    do i = 1,10!while(1>0)!error>eps)
        write(2,*) p(nx/2,:)
        write(1,*) p(:,ny/2)
        p2 = pp
        pp = timestep(nx,ny,xvals,yvals,pp,p)
        error = maxval(pp-p)/maxval(p)
        write(*,*) 'error = ',error,'norm = ',sum(p),i
        p = p2
       ! i = i+1
    enddo

    open(19,file = 'fpres/jx.dat')
    do i = 1,nx
        write(19,*) j(1,i,:)
    enddo
    close(19)
    open(21,file = 'fpres/jy.dat')
    do i = 1,nx
        write(21,*) j(2,i,:)
    enddo
    close(21)
    open(22,file = 'fpres/divj.dat')
    do i = 1,nx
        write(22,*) dj(i,:)
    enddo
    close(22)
    open(20,file = 'fpres/dp.dat')
    do i = 1,nx
        write(20,*) (pp(i,:)-p(i,:))
    enddo
    close(20)

    close(1)
    close(2)
    p = p2        
endsubroutine fullrelax

function div(jx,jy,xvals,yvals)
    implicit none
    real *8::jx(:,:),jy(:,:),xvals(:),yvals(:),hx,hy,ymax
    real *8,allocatable::div(:,:)
    integer::nx,ny,i
    ymax = consts(1)+consts(2)
    nx = size(xvals)
    ny = size(yvals)
    hx = xvals(2)-xvals(1)
    allocate(div(nx,ny))
    div(:,:) = 0
    do i = 1,ny
        div(:,i) = compder(jx(:,i),nx)/hx
    enddo
    do i = 1,nx
        div(i,:) = div(i,:) + specder(yvals(1),ymax,ny,jy(i,:))
    enddo
endfunction div

function getcurrent(p,xvals,yvals)
    implicit none
    real *8::p(:,:),temp(2),xvals(:),yvals(:),ymax,b(2,2),hx,hy
    real *8,allocatable::getcurrent(:,:,:),ax(:,:),ay(:,:)
    integer::nx,ny,i,j
    ymax = consts(1)+consts(2)
    hx = xvals(2)-xvals(1)
    hy = yvals(2)-yvals(1)
    b = barray()
    nx = size(p(:,1))
    ny = size(p(1,:))
    allocate(getcurrent(2,nx,ny),ax(nx,ny),ay(nx,ny))
    do i = 1,nx
    do j = 1,ny
        temp = avec(xvals(i),yvals(j))
        ax(i,j) = temp(1)
        ay(i,j) = temp(2)
    enddo
    enddo
    
    do i = 1,ny
        getcurrent(1,:,i) = ax(:,i)*p(:,i)+b(1,1)*compder(p(:,i),nx)/hx
        getcurrent(2,:,i) = ay(:,i)*p(:,i)+b(1,2)*compder(p(:,i),nx)/hx
    enddo
    do i = 1,nx
        getcurrent(1,i,:) = getcurrent(1,i,:) + b(1,2)*specder(yvals(1),ymax,ny,p(i,:))
        getcurrent(2,i,:) = getcurrent(2,i,:) + b(2,2)*specder(yvals(1),ymax,ny,p(i,:))
    enddo

endfunction getcurrent

function timestep(nx,ny,xvals,yvals,p,pold)
    implicit none
    complex *16,dimension(nx,ny)::rs,cp,cpft,tstep_cp,cpold,axp
    complex *16::ls(nx,nx),rsvec(nx)
    real *8::xvals(nx),yvals(ny),p(nx,ny),pold(nx,ny),k(ny),ymax,timestep(nx,ny),ax(nx,ny)
    real *8::temp(2),hx
    real *8, parameter:: pi = acos(-1.0d0)
    integer::nx,ny
    integer::i,info,ipiv(nx),j
    cp = cmplx(p,0.0D0,kind=16)
    cpold = cmplx(pold,0.0D0,kind=16)
    ymax = consts(1)+consts(2)
    k = 2*pi*fftfreq(ny)/(ymax-yvals(1))
    hx = xvals(2)-xvals(1)
    rs = 1/2*(3*gtilde(p,xvals,yvals)-gtilde(pold,xvals,yvals))

    do i=1,nx
        do j=1,ny
            temp = avec(xvals(i),yvals(j))
            ax(i,j) = temp(1)
        enddo
    enddo
    do i = 1,nx
        axp(i,:) = easy_fft(ax(i,:)*1/2*(3*cp(i,:)-cpold(i,:)),ny)
    enddo

    do i = 1,nx
        cpft(i,:) = easy_fft(cp(i,:),ny)
    enddo
    do i = 1,ny
        rsvec = rs(:,i) + matvec(rhs_der_oper(nx,k(i),hx),cpft(:,i),nx)
        ls = der_oper(nx,k(i),hx)
        call boundfix(ls,rsvec,k(i),cp(:,i),axp(:,i),hx)
        call zgesv(nx,1,ls,nx,ipiv,rsvec,nx,info)
        tstep_cp(:,i) = rsvec
    enddo
    do i = 1,nx
        timestep(i,:) = real(easy_ifft(tstep_cp(i,:),ny))
    enddo
endfunction timestep

function timestep_eul(nx,ny,xvals,yvals,p)
    implicit none
    complex *16,dimension(nx,ny)::rs,cp,cpft,tstep_cp,axp
    complex *16::ls(nx,nx),rsvec(nx),rsvec_old(nx)
    real *8::xvals(nx),yvals(ny),p(nx,ny),k(ny),ymax,timestep_eul(nx,ny),ax(nx,ny)
    real *8::temp(2),hx
    real *8, parameter:: pi = acos(-1.0d0)
    integer::nx,ny
    integer::i,info,ipiv(nx),j
    cp = cmplx(p,0.0D0,kind=16)
    ymax = consts(1)+consts(2)
    k = 2*pi*fftfreq(ny)/(ymax-yvals(1))
    hx = xvals(2)-xvals(1)
    
    rs = gtilde(p,xvals,yvals)

    do i=1,nx
        do j=1,ny
            temp = avec(xvals(i),yvals(j))
            ax(i,j) = temp(1)
        enddo
    enddo
    do i = 1,nx
        axp(i,:) = easy_fft(ax(i,:)*cp(i,:),ny)
    enddo

    do i = 1,nx
        cpft(i,:) = easy_fft(cp(i,:),ny)
    enddo

    do i = 1,ny
        rsvec = rs(:,i) + matvec(rhs_der_oper(nx,k(i),hx),cpft(:,i),nx)
        ls = der_oper(nx,k(i),hx)
        call boundfix(ls,rsvec,k(i),cpft(:,i),axp(:,i),hx)
        call zgesv(nx,1,ls,nx,ipiv,rsvec,nx,info)
        tstep_cp(:,i) = rsvec
    enddo
    do i = 1,nx
        timestep_eul(i,:) = real(easy_ifft(tstep_cp(i,:),ny))
    enddo
endfunction timestep_eul

subroutine boundfix(ls,rsvec,k,cp,axp,hx)
    implicit none
    complex *16::ls(:,:),rsvec(:),cp(:),ik,axp(:)
    real *8::k,b(2,2),hx
    integer::nx,i
    nx = size(cp)
    b = barray()
    ik = cmplx(0.0D0,1.0D0,kind=16)*k

    ls(1,1) = -3/(2*hx)
    ls(1,2) = 2/hx
    ls(1,3) = -1/(2*hx)

    ls(nx,nx) = 3/(2*hx)
    ls(nx,nx-1) = -2/hx
    ls(nx,nx-2) = 1/(2*hx)

    rsvec(1) = -ik*b(1,2)/b(1,1)*cp(1) + axp(1)/(b(1,1))
    rsvec(nx) = -ik*b(1,2)/b(1,1)*cp(nx) + axp(nx)/(b(1,1))
endsubroutine boundfix

function rhs_der_oper(nx,k,hx)
    implicit none
    complex *16::rhs_der_oper(nx,nx),ik
    real *8::k,b(2,2),hx
    integer::nx,ny,i,j
    b = barray()/2
    
    ik = k*cmplx(0.0D0,1.0D0,kind=16) !double precision imaginary unit, needs revision when k is
!    incorporated in the overall ODE solver

    
    rhs_der_oper(:,:) = 0
    rhs_der_oper(1,1) = 1 - 3/(2*hx)*2*ik*b(1,2) + 2/(hx**2)*b(1,1)-(ik)
    rhs_der_oper(1,2) = 4/hx*ik*b(1,2)-5/(hx**2)*b(1,1)
    rhs_der_oper(1,3) = -1/(hx)*ik*b(1,2)+4/(hx**2)*b(1,1)
    rhs_der_oper(1,4) = -1/(hx**2)*b(1,1)
    rhs_der_oper(nx,nx) = 1 + 3/(2*hx)*2*ik*b(1,2) + 2/(hx**2)*b(1,1)-(ik)
    rhs_der_oper(nx,nx-1) = - 4/hx*ik*b(1,2)-5/(hx**2)*b(1,1)
    rhs_der_oper(nx,nx-2) = 1/(hx)*ik*b(1,2)+4/(hx**2)*b(1,1)
    rhs_der_oper(nx,nx-3) = -1/(hx**2)*b(1,1)
    
    do i = 2,nx-1
        rhs_der_oper(i,i-1) = 1/2*(2*b(1,1)/(hx**2)-2*ik*b(1,2)/hx)
        rhs_der_oper(i,i+1) = 1/2*(2*b(1,1)/(hx**2)+2*ik*b(1,2)/hx)
        rhs_der_oper(i,i) = 1+(ik)**2*b(2,2)-2*b(1,1)/(hx**2)
    enddo
endfunction rhs_der_oper

function der_oper(nx,k,hx)
    implicit none
    complex *16::der_oper(nx,nx),ik
    real *8::k,b(2,2),hx
    real *8, parameter:: pi = acos(-1.0d0)
    integer::nx,ny,i,j
    b = barray()/2
    
    ik = k*cmplx(0.0D0,1.0D0,kind=16) !double precision imaginary unit

    der_oper(:,:) = 0
    do i = 2,nx-1
        der_oper(i,i-1) = -1/2*(2*b(1,1)/(hx**2)-2*ik*b(1,2)/hx)
        der_oper(i,i+1) = -1/2*(2*b(1,1)/(hx**2)+2*ik*b(1,2)/hx)
        der_oper(i,i) = 1-(ik)**2*b(2,2)+2*b(1,1)/(hx**2)
    enddo
endfunction der_oper

function gtilde(p,xvals,yvals)
    implicit none
    real *8::p(:,:),xvals(:),yvals(:)
    complex *16,allocatable::gtilde(:,:),step(:,:)
    integer::nx,ny,i,j

    nx = size(p(:,1))
    ny = size(p(1,:))
    allocate(gtilde(nx,ny),step(nx,ny))
        
    step = (term1(p,xvals,yvals)+term2(p,xvals,yvals))*cmplx(1.0D0,0.0D0,kind=16)
    gtilde(:,:) = 0
    do i = 1,nx
        gtilde(i,:) = easy_fft(step(i,:),ny)
    enddo
endfunction gtilde

   function consts(m)
       implicit none
       integer::m
       real *8::consts,lsmall,lbig,k,tau,tbar,gm,umax,kmult,vals(8)

       lsmall = 3! 
       lbig = 7  !These two set the relative width of the two parts of the potential. 
       k = 0.1 !Stiffness of spring
       tau = 0 !dimensionless temperature difference (dT/Tbar)
       tbar = 0.5 !average temperature
       umax = 1 !potential height
       gm = 0.0001 ! inverse friction coefficient. Essentially acts as dt for timestepping
       kmult = 4 !Since we can't sample the full spring stretching direction, we need to cut it off somewhere
       !This sets the cutoff at some multiple of the thermal length of the spring

       vals = (/lsmall,lbig,k,tau,tbar,umax,gm,kmult/)
       consts = vals(m)

    endfunction consts 
           
    function term3(p,xvals,yvals)
        implicit none
        real *8::p(:,:),xvals(:),yvals(:),b(2,2),hx,hy,ymax
        integer::nx,ny,i,j
        real *8,allocatable::term3(:,:),px(:,:),py(:,:),pxx(:,:),pyy(:,:),pxy(:,:),pxy2(:,:)

        ymax = consts(1)+consts(2)
        nx = size(p(:,1))
        ny = size(p(1,:))
        hx = xvals(2)-xvals(1)
        hy = yvals(2)-yvals(1)
        allocate(term3(nx,ny),px(nx,ny),pyy(nx,ny),pxx(nx,ny),pxy(nx,ny))
        px(:,:) = 0
        pxy(:,:) = 0
        pxx(:,:) = 0
        pyy(:,:) = 0
        term3(:,:) = 0
        b = barray()
        do i = 1,nx
            pyy(i,:) = comp2der(p(i,:),ny)/(hy*hy)
            !pyy(i,:) = spec2der(yvals(1),ymax,ny,p(i,:))
        enddo
        do i = 1,ny
            pxx(:,i) = comp2der(p(:,i),nx)/(hx*hx)
        enddo
        do i = 1,ny
            px(:,i) = compder(p(:,i),nx)/hx
        enddo
        do i = 1,nx
            pxy(i,:) = compder(p(i,:),ny)/hy
            !pxy(i,:) = specder(yvals(1),ymax,ny,px(i,:))
        enddo

        term3 = b(1,1)*pxx+b(2,2)*pyy+2*b(1,2)*pxy

    endfunction term3

    function term2(p,xvals,yvals)!term 2 looks like a . grad(p)
        implicit none
        real *8::p(:,:),xvals(:),yvals(:),temp(2),hx,hy,ymax
        integer::nx,ny,i,j
        real *8,allocatable::term2(:,:),ax(:,:),ay(:,:),px(:,:),py(:,:)

        ymax = consts(1)+consts(2)
        nx = size(p(:,1))
        ny = size(p(1,:))
        hx = xvals(2)-xvals(1)
        hy = yvals(2)-yvals(1)

        allocate(term2(nx,ny),ax(nx,ny),ay(nx,ny),px(nx,ny),py(nx,ny))
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
            !py(i,:) = specder(yvals(1),ymax,ny,p(i,:))
            py(i,:) = compder(p(i,:),ny)/hy
        enddo
        do i = 1,ny
            px(:,i) = compder(p(:,i),nx)/hx
        enddo
        term2 = ax*px+ay*py

    endfunction term2

    function term1(p,xvals,yvals)!term 1 looks like div(a)
        implicit none
        real *8::p(:,:),xvals(:),yvals(:)
        integer::nx,ny,i,j
        real *8,allocatable::term1(:,:)

        nx = size(p(:,1))-2
        ny = size(p(1,:))-2
    
        allocate(term1(nx,ny))

        term1(:,:) = 0

        do i = 1,nx
            do j = 1,ny
                term1(i,j) = p(i,j)*diva(xvals(i),yvals(j))
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
        diva = -2*consts(7)*(k+5.0d0/4*(upp(y+x/2)+upp(y-x/2)))
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
      avec(1) = 2*fx*gm
      avec(2) = fy*gm/(2)
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
      barray = (-1)*tbar*gm/(2)*barray

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
      implicit none
      real *8::dy,arr(n),comp2der(n)
      integer::n,i
      comp2der(1) = 2*arr(1) - 5*arr(2) + 4*arr(3) - arr(4)
      comp2der(n) = 2*arr(n) - 5*arr(n-1) + 4*arr(n-2) - arr(n-3)
      do i = 2,n-1
        comp2der(i) = (arr(i+1)-2*arr(i)+arr(i-1))
      enddo

    endfunction comp2der

    function compder(arr,n) !second order finite difference first derivative,
      implicit none
      real *8::dy,arr(n),compder(n)
      integer::n,i
      compder(1) = -3/2*arr(1) + 2*arr(2) - 1/2*arr(3)
      compder(n) = 3/2*arr(n) - 2*arr(n-1) + 1/2*arr(n-2) 
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
        complex(c_double_complex),dimension(n)::easy_fft,input_copy,cinput
            
        input_copy = cinput
        call dfftw_plan_dft_1d_(forward,n, input_copy,easy_fft,fftw_forward,fftw_estimate)
        call dfftw_execute_(forward)
        call dfftw_destroy_plan_(forward)
        call fftshift(easy_fft)
    endfunction easy_fft

    function easy_ifft(cinput,n)
        use,intrinsic :: iso_c_binding
        implicit none
        include '/usr/local/include/fftw3.f03'
        type(c_ptr) :: plan
        integer:: n
        integer *8::backward
        complex(c_double_complex),dimension(n)::easy_ifft,input_copy,cinput
            
        input_copy = cinput
        call ifftshift(input_copy)
        call dfftw_plan_dft_1d_(backward,n, input_copy,easy_ifft,fftw_backward,fftw_estimate)
        call dfftw_execute_(backward)
        call dfftw_destroy_plan_(backward)
        easy_ifft = easy_ifft/n
    endfunction easy_ifft

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
        call rfftshift(fftfreq)
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

    function matvec(mat,vec,n)
        implicit none
        complex *16::mat(n,n),vec(n),matvec(n)
        integer::i,n,j
        matvec(:) = 0
        do i = 1,n
        do j = 1,n
            matvec(i) = matvec(i) + mat(i,j)*vec(j)
        enddo
        enddo
    endfunction matvec
    
    subroutine fftshift(vec)
        implicit none
        complex *16::vec(:)
        complex *16,allocatable::dummy(:)
        integer::i,n
        n = size(vec)
        allocate(dummy(n))
        dummy = vec
        if (mod(n,2)==0) then
            do i = 1,n/2
                vec(i) = dummy(i+n/2)
            enddo
            vec(n/2+1) = dummy(1)
            do i = n/2+2,n
                vec(i) = dummy(i-n/2)
            enddo
        else
            do i = 1,(n-1)/2
                vec(i) = dummy((n+1)/2+i)
            enddo
            vec((n+1)/2) = dummy(1)
            do i = (n+3)/2,n
                vec(i) = dummy(i-(n-1)/2)
            enddo
        endif
    endsubroutine fftshift

    subroutine ifftshift(vec)
        implicit none
        complex *16::vec(:)
        complex *16,allocatable::dummy(:)
        integer::i,n
        n = size(vec)
        allocate(dummy(n))
        dummy = vec
        if (mod(n,2)==0) then
            do i = 1,n/2
                vec(i+n/2) = dummy(i)
            enddo
            vec(1) = dummy(n/2+1)
            do i = n/2+2,n
                vec(i-n/2) = dummy(i)
            enddo
        else
            do i = 1,(n-1)/2
                vec(i+(n+1)/2) = dummy(i)
            enddo
            vec(1) = dummy((n+1)/2)
            do i = (n+3)/2,n
                vec(i-(n-1)/2) = dummy(i)
            enddo
        endif   
    endsubroutine ifftshift

    subroutine rfftshift(vec)
        implicit none
        real *8::vec(:)
        real *8,allocatable::dummy(:)
        integer::i,n
        n = size(vec)
        allocate(dummy(n))
        dummy = vec
        if (mod(n,2)==0) then
            do i = 1,n/2
                vec(i) = dummy(i+n/2)
            enddo
            vec(n/2+1) = dummy(1)
            do i = n/2+2,n
                vec(i) = dummy(i-n/2)
            enddo
        else
            do i = 1,(n-1)/2
                vec(i) = dummy((n+1)/2+i)
            enddo
            vec((n+1)/2) = dummy(1)
            do i = (n+3)/2,n
                vec(i) = dummy(i-(n-1)/2)
            enddo
        endif
    endsubroutine rfftshift

    subroutine irfftshift(vec)
        implicit none
        real *8::vec(:)
        real *8,allocatable::dummy(:)
        integer::i,n
        n = size(vec)
        allocate(dummy(n))
        dummy = vec
        if (mod(n,2)==0) then
            do i = 1,n/2
                vec(i+n/2) = dummy(i)
            enddo
            vec(1) = dummy(n/2+1)
            do i = n/2+2,n
                vec(i-n/2) = dummy(i)
            enddo
        else
            do i = 1,(n-1)/2
                vec(i+(n+1)/2) = dummy(i)
            enddo
            vec(1) = dummy((n+1)/2)
            do i = (n+3)/2,n
                vec(i-(n-1)/2) = dummy(i)
            enddo
        endif   
    endsubroutine irfftshift


endmodule fancy_timestep
