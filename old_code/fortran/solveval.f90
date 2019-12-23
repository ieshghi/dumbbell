module solveval
contains

    function genlhs(n,sig,mult)
        implicit none
        integer::n,mult
        real *8::sig,genlhs(n,n),xvals(n),h
        call gengrid(n,h,xvals,sig,mult)
        
        genlhs = comp2der(n)/(h**2) - rder(n,xvals)/sig

        genlhs(1,:) = 0
        genlhs(1,3) = 1
        genlhs(1,2) = -2*h*xvals(1)/sig

        genlhs(n,:) = 0
        genlhs(n,n-2) = 1
        genlhs(n,n-1) = 2*h*xvals(n)/sig

    endfunction genlhs

    subroutine gengrid(nx,hx,xvals,sig,mult)
       implicit none
       integer::nx,ny,mult
       real *8::xmax,hx,xspr,xvals(nx),sig

       xspr = sqrt(2*sig)
       xmax = mult*xspr
       
       hx = 2*xmax/nx

       xvals = linspace(-xmax,xmax,nx)
       !include the endpoint

    endsubroutine gengrid
    
    function rder(m,arr)
        implicit none
        integer::m,i
        real *8::arr(m),rder(m,m),inter(m,m),h
        
        h = arr(2)-arr(1)
        inter = compder(m)/h

        do i = 1,m
            rder(i,i) = inter(i,i)*arr(i)
        enddo

    endfunction rder

    function comp2der(n)!second order finite difference first derivative,
      !assumes ghost cells and doesn't touch them
      implicit none
      real *8::comp2der(n,n)
      integer::n,i

      comp2der(:,:) = 0
      do i = 2,n-1
        comp2der(i,i+1) = 1
        comp2der(i,i) = -2
        comp2der(i,i-1) = 1
      enddo

    endfunction comp2der

    function compder(n) !second order finite difference first derivative,
      !assumes ghost cells and doesn't touch them
      implicit none
      real *8::compder(n,n)
      integer::n,i
        
      compder(:,:) = 0
      do i = 2,n-1
        compder(i,i+1) = 1.0d0/2
        compder(i,i-1) = -1.0d0/2
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

endmodule solveval
