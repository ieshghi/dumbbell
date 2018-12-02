module solvefp
contains

    subroutine fullrelax(p,nx,ny)
        implicit none
        real *8::p(nx,ny),xvals(nx),yvals(ny),hx,hy,eps,error,p2(nx,ny)
        
        call gengrid(nx,ny,xvals,yvals,hx,hy)
        eps = 1e-12
        error = 10

        while(error>eps)do
            p2 = p
            relaxstep(p,xvals,yvals)
            error = max(abs(p))
        enddo

        p = p2        
    endsubroutine

    subroutine relaxstep(p,xvals,yvals)
        implicit none
        real *8::p(:,:),xvals(:),yvals(:),hx,temp(2)
        real *8,allocatable::p2(:,:),ax(:,:)
        integer::nx,ny,i,j

        nx = size(p(:,1))-2
        ny = size(p(1,:))-2
        hx = xvals(2)-xvals(1)
        allocate(p2(nx+2,ny+2),ax(nx,ny))
        
        do i = 1,nx
            do j = 1,ny
                temp = avec(xvals(i),yvals(j))
                ax(i,j) = temp(1)
            enddo
        enddo

        p2 = term1(p,xvals,yvals)+term2(p,xvals,yvals)+term3(p,xvals,yvals)
        
        call update_ghosts(ax,p2,hx)
        p = p2
    endsubroutine relaxstep

    function term3(p,xvals,yvals)
        implicit none
        real *8::p(:,:),xvals(:),yvals(:),b(2,2),hx,hy
        integer::nx,ny,i,j
        real *8,allocatable::term3(:,:),px(:,:),py(:,:),pxx(:,:),pyy(:,:),pxy(:,:)

        nx = size(p(:,1))-2
        ny = size(p(1,:))-2
        hx = xvals(2)-xvals(1)
        hy = yvals(2)-yvals(1)
        allocate(term3(nx+2,ny+2),px(nx+2,ny+2),py(nx+2,ny+2),pyy(nx+2,ny+2),pxx(nx+2,ny+2),pxy(nx+2,ny+2))
        term3(:,:) = 0
        b = barray()
        do i = 1,nx
            pyy(i+1,:) = comp2der(p(i+1,:),ny)/(hy*hy)
        enddo
        do i = 1,ny
            pxx(:,i+1) = comp2der(p(:,i+1),nx)/(hx*hx)
        enddo
        do i = 1,nx
            py(i+1,:) = compder(p(i+1,:),ny)/hy
        enddo
        do i = 1,ny
            pxy(:,i+1) = compder(py(:,i+1),nx)/hx
        enddo

        term3 = b(1,1)*pxx+b(2,2)*pyy+2*b(1,2)*pxy
    endfunction term3

    function term2(p,xvals,yvals)!term 2 looks like a . grad(p)
        implicit none
        real *8::p(:,:),xvals(:),yvals(:),temp(2),hx,hy
        integer::nx,ny,i,j
        real *8,allocatable::term2(:,:),ax(:,:),ay(:,:),px(:,:),py(:,:)

        nx = size(p(:,1))-2
        ny = size(p(1,:))-2
        hx = xvals(2)-xvals(1)
        hy = yvals(2)-yvals(1)

        allocate(term2(nx+2,ny+2),ax(nx,ny),ay(nx,ny),px(nx+2,ny+2),py(nx+2,ny+2))
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
    
    subroutine update_ghosts(ax,p,hx)
        implicit none
        real *8::ax(:,:),p(:,:),hx,hy,b(2,2)
        real *8,allocatable:: p_ng(:,:),pout(:,:)
        integer::nx,ny,i,j

        b = barray()

        nx = size(p(:,1))-2
        ny = size(p(1,:))-2
        
        allocate(p_ng(nx,ny),pout(nx+2,ny+2))
        p_ng = p(2:nx+1,2:ny+1)

        do i = 1,nx
            p(i+1,ny+2) = p(i,2)
            p(i+1,1) = p(i,ny+1)
        enddo
        do i = 1,ny
            p(nx+2,i+1) = -2*hx/(b(1,2)+b(1,1))*ax(nx,i)*p_ng(nx,i)+p_ng(nx-1,i)
            p(1,i+1) = 2*hx/(b(1,2)+b(1,1))*ax(1,i)*p_ng(1,i)+p_ng(2,i)
        enddo

        p(1,1) = 0
        p(nx+2,1) = 0
        p(1,ny+2) = 0
        p(nx+2,ny+2) = 0

    endsubroutine update_ghosts

    subroutine gengrid(nx,ny,hx,hy,xvals,yvals)
       implicit none
       integer::nx,ny
       real *8::xmax,ymax,hx,hy,k,tbar,xspr,xvals(nx),yvals(ny)
       
       ymax = consts(1)+consts(2)
       k = consts(3)
       tbar = consts(5)

       xspr = sqrt(2*tbar/k)
       xmax = 2*consts(8)*xspr
       
       hy = ymax/ny
       hx = xmax/nx

       xvals = linspace(-1.0d0,1.0d0,nx)
       yvals = linspace2(0.0d0,ymax,ny) !y is a periodic coordinate so we don't
       !include the endpoint

    endsubroutine gengrid
    
    function consts(m)
       implicit none
       integer::m
       real *8::consts,lsmall,lbig,k,tau,tbar,gm,umax,kmult,vals(8)

       lsmall = 1
       lbig = 9
       k=1
       tau = 0
       tbar = 1
       umax = 1
       gm = 1
       kmult = 5
       vals = (/lsmall,lbig,k,tau,tbar,umax,gm,kmult/)
       consts = vals(m)

    endfunction consts 

    function diva(x,y)
        implicit none
        real *8::x,y,k,diva
        k = consts(3)
        diva = -k-5.0d0/4*(upp(y+x/2)+upp(y-x/2))

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
        if mod(z,l)<lsmall then
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
      barray = tbar/(2*gm)*barray

    endfunction barray

    function comp2der(arr,n)!second order finite difference first derivative,
      !assumes ghost cells and doesn't touch them
      implicit none
      real *8::dy,arr(n),comp2der(n)
      integer::n,i
      comp2der(1) = arr(1)
      comp2der(n) = arr(n)

      do i = 2,n-1
        compder(i) = (arr(i+1)-2*arr(i)+arr(i-1)
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
