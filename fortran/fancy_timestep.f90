module fancy_timestep
contains

subroutine fullrelax(p,nx,ny)
    implicit none
    integer::nx,ny,i
    real *8::p(nx,ny),xvals(nx),yvals(ny),hx,hy,eps,error,p2(nx,ny),pp(nx,ny)
        
    call gengrid(nx,ny,hx,hy,xvals,yvals)
    eps = 1e-12
    error =10

    open(1,file = 'fpres/xder.dat')
    open(2,file = 'fpres/yder.dat')

    pp = timestep_eul(nx,ny,xvals,yvals,p)
    do i = 1,5!while(error>eps)
        write(2,*) p(nx/2,:)
        write(1,*) p(:,ny/2)
            
        p2 = pp
        pp = timestep(nx,ny,xvals,yvals,pp,p)
        p = p2

        error = sqrt(sum((pp-p)**2))
        write(*,*) 'error = ',error,'norm = ',sum(p2)
    enddo

    open(1,file = 'fpres/pp.dat')
    do i = 1,nx
        write(1,*) pp(i,:)
    enddo
    close(1)


    close(1)
    close(2)
    p = p2        
endsubroutine fullrelax

function timestep(nx,ny,xvals,yvals,p,pold)
    implicit none
    complex *16,dimension(nx,ny)::rs,cp,cpft
    complex *16::ls(nx,nx),rsvec(nx)
    real *8::xvals(nx),yvals(ny),p(nx,ny),pold(nx,ny),k(ny),ymax,timestep(nx,ny)
    real *8, parameter:: pi = acos(-1.0d0)
    integer::nx,ny
    integer::i,info,ipiv(nx)
    cp = cmplx(p,0.0D0,kind=16)

    ymax = consts(1)+consts(2)
    k = 2*pi*fftfreq(ny)/(ymax-yvals(1))

    rs = 1/2*(3*gtilde(p,xvals,yvals)-gtilde(pold,xvals,yvals))
    do i = 1,nx
        cpft(i,:) = easy_fft(cp(i,:),ny)
    enddo
    do i = 1,ny
        ls = der_oper(nx,k(i))
        rsvec = rs(:,i) + matvec(rhs_der_oper(nx,k(i)),cpft(:,i),nx)
        call boundfix(ls,rsvec,k(i),cp(:,i))
        call zgesv(nx,1,ls,nx,ipiv,rsvec,nx,info)
        timestep(:,i) = abs(easy_ifft(rsvec/ny,ny))
    enddo

endfunction timestep

function timestep_eul(nx,ny,xvals,yvals,p)
    implicit none
    complex *16,dimension(nx,ny)::rs,cp,cpft,tstep_cp
    complex *16::ls(nx,nx),rsvec(nx)
    real *8::xvals(nx),yvals(ny),p(nx,ny),k(ny),ymax,timestep_eul(nx,ny)
    real *8, parameter:: pi = acos(-1.0d0)
    integer::nx,ny
    integer::i,info,ipiv(nx)
    cp = cmplx(p,0.0D0,kind=16)

    ymax = consts(1)+consts(2)
    k = 2*pi*fftfreq(ny)/(ymax-yvals(1))

    rs = gtilde(p,xvals,yvals)

    do i = 1,nx
        cpft(i,:) = easy_fft(cp(i,:),ny)
    enddo
!    open(3,file='err/rsvec.dat')
!    open(4,file='err/isvec.dat')
    do i = 1,ny
        ls = der_oper(nx,k(i))
        rsvec = rs(:,i) + matvec(rhs_der_oper(nx,k(i)),cpft(:,i),nx)
        call boundfix(ls,rsvec,k(i),cp(:,i))
        call zgesv(nx,1,ls,nx,ipiv,rsvec,nx,info)
        tstep_cp(:,i) = easy_ifft(rsvec,ny)
        timestep_eul(:,i) = abs(easy_ifft(rsvec/ny,ny))
!        write(3,*) real(tstep_cp(:,i))
!        write(4,*) imag(tstep_cp(:,i))
    enddo
!    close(3)
!    close(4)
endfunction timestep_eul

subroutine boundfix(ls,rsvec,k,cp)
    implicit none
    complex *16::ls(:,:),rsvec(:),cp(:),ik
    real *8::k,b(2,2)
    integer::nx,ny,i
    nx = size(cp)
    b = barray()
    !This is very fishy... assume p = 0 at boundaries and solve for b.c. B_xx*d_x P = - B_xy*d_y P
    ik = cmplx(0.0D0,1.0D0,kind=16)*k

    ls(1,1) = -3/2
    ls(1,2) = 2
    ls(1,3) = -1/2

    ls(nx,nx) = 3/2
    ls(nx,nx-1) = -2
    ls(nx,nx-2) = 1/2

    rsvec(1) = -ik*b(1,2)/b(1,1)*cp(1)
    rsvec(nx) = -ik*b(1,2)/b(1,1)*cp(nx)
endsubroutine boundfix

function rhs_der_oper(nx,k)
    implicit none
    complex *16::rhs_der_oper(nx,nx),ik
    real *8::k,b(2,2)
    integer::nx,ny,i,j
    b = barray()/2
    
    ik = k*cmplx(0.0D0,1.0D0,kind=16) !double precision imaginary unit, needs revision when k is
!    incorporated in the overall ODE solver

    rhs_der_oper(:,:) = 0
    do i = 2,nx-1
        rhs_der_oper(i,i-1) = 1/2*(b(1,1)-2*ik*b(1,2))
        rhs_der_oper(i,i+1) = 1/2*(b(1,1)+2*ik*b(1,2))
        rhs_der_oper(i,i) = 1-(ik)**2*b(2,2)+b(1,1)
    enddo
endfunction rhs_der_oper

function der_oper(nx,k)
    implicit none
    complex *16::der_oper(nx,nx),ik
    real *8::k,b(2,2)
    real *8, parameter:: pi = acos(-1.0d0)
    integer::nx,ny,i,j
    b = barray()/2
    
    ik = k*cmplx(0.0D0,1.0D0,kind=16) !double precision imaginary unit

    der_oper(:,:) = 0
    do i = 2,nx-1
        der_oper(i,i-1) = -1/2*(b(1,1)-2*ik*b(1,2))
        der_oper(i,i+1) = -1/2*(b(1,1)+2*ik*b(1,2))
        der_oper(i,i) = 1+(ik)**2*b(2,2)-b(1,1)
    enddo
endfunction der_oper

function gtilde(p,xvals,yvals)
    implicit none
    real *8::p(:,:),xvals(:),yvals(:),temp(2)
    real *8,allocatable,dimension(:,:)::ax,ay
    complex *16,allocatable::gtilde(:,:),step(:,:)
    integer::nx,ny,i,j

    nx = size(p(:,1))
    ny = size(p(1,:))
    allocate(gtilde(nx,ny),ax(nx,ny),ay(nx,ny),step(nx,ny))
        
    do i = 1,nx
        do j = 1,ny
           temp = avec(xvals(i),yvals(j))
           ax(i,j) = temp(1)
           ay(i,j) = temp(2)
        enddo
    enddo

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

       lsmall = 1! 
       lbig = 9  !These two set the relative width of the two parts of the potential. 
       k = 0.1 !Stiffness of spring
       tau = 0 !dimensionless temperature difference (dT/Tbar)
       tbar = 0.5 !average temperature
       umax = 1 !potential height
       gm = 100 !friction coefficient. Essentially acts as 1/dt for timestepping
       kmult = 10 !Since we can't sample the full spring stretching direction, we need to cut it off somewhere
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
            !pyy(i+1,:) = comp2der(p(i+1,:),ny)/(hy*hy)
            pyy(i,:) = spec2der(yvals(1),ymax,ny,p(i,:))
        enddo
        do i = 1,ny
            pxx(:,i) = comp2der(p(:,i),nx)/(hx*hx)
        enddo
        do i = 1,ny
            px(:,i) = compder(p(:,i),nx)/hx
        enddo
        do i = 1,nx
            !pxy(i+1,:) = compder(p(i+1,:),ny)/hy
            pxy(i,:) = specder(yvals(1),ymax,ny,px(i,:))
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
            py(i,:) = specder(yvals(1),ymax,ny,p(i,:))
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
        complex(c_double_complex),dimension(n)::cinput,easy_fft
            
        call dfftw_plan_dft_1d_(forward,n, cinput,easy_fft,fftw_forward,fftw_estimate)
        call dfftw_execute_(forward)
        call dfftw_destroy_plan_(forward)
    endfunction easy_fft

    function easy_ifft(cinput,n)
        use,intrinsic :: iso_c_binding
        implicit none
        include '/usr/local/include/fftw3.f03'
        type(c_ptr) :: plan
        integer:: n
        integer *8::backward
        complex(c_double_complex),dimension(n)::cinput,easy_ifft
            
        call dfftw_plan_dft_1d_(backward,n, cinput,easy_ifft,fftw_backward,fftw_estimate)
        call dfftw_execute_(backward)
        call dfftw_destroy_plan_(backward)

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
endmodule fancy_timestep
