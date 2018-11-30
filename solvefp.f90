module solvefp
contains


    
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

    function der_2arr(n)
        implicit none
        integer::n,i
        real *8::der_2arr(n,n)

        der_2arr(:,:) = 0
        
        der_2arr(1,1) = 2
        der_2arr(1,2) = -5
        der_2arr(1,3) = 4
        der_2arr(1,4) = 1

        do i = 2,n-1
            der_2arr(i,i-1) = 1
            der_2arr(i,i) = -2
            der_2arr(i,i+1) = 1
        enddo
        
        der_2arr(n,n) = 2
        der_2arr(n,n-1) = -5
        der_2arr(n,n-2) = 4
        der_2arr(n,n-3) = 1
    endfunction der_2arr

    function der_arr(n)
        implicit none
        integer::n,i
        real *8::der_arr(n,n)
        der_arr(:,:) = 0
        der_arr(1,1) = -3/2
        der_arr(1,2) = 2
        der_arr(1,3) = -1/2

        do i = 2,n-1
            
            der_arr(i,i+1) = 1/2
            der_arr(i,i-1) = -1/2

        enddo

        der_arr(n,n-2) = 1/2
        der_arr(n,n-1) = -2
        der_arr(n,n) = 3/2
    endfunction der_arr


    subroutine gengrid(nx,ny,grid,hx,hy)
       implicit none
       integer::nx,ny
       real *8::xmax,ymax,hx,hy,k,tbar,xspr,grid(nx+2,ny+2)
       
       ymax = consts(1)+consts(2)
       k = consts(3)
       tbar = consts(5)

       xspr = sqrt(2*tbar/k)
       xmax = 2*consts(8)*xspr
       
       hy = ymax/ny
       hx = xmax/nx
       grid(:,:) = 0

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

    function compder(arr,n,lmax)
      !Fourth order convergence tested. Error starts blowing up @ n~1000
      implicit none
      real *8::dy
      complex *16::arr(n)
      complex *16::compder(n)
      integer::n,i
        
      dy = lmax/n
        
      compder(1) = -25.0d0/12*arr(1) + 4*arr(2) - 3*arr(3) + 4.0d0/3*arr(4) - 1.0d0/4*arr(5)
      compder(2) = -25.0d0/12*arr(2) + 4*arr(3) - 3*arr(4) + 4.0d0/3*arr(5) - 1.0d0/4* arr(6)
      compder(n) = +25.0d0/12*arr(n) - 4*arr(n-1) + 3*arr(n-2) - 4.0d0/3*arr(n-3) + 1.0d0/4*arr(n-4)
      compder(n-1) = +25.0d0/12*arr(n-1) - 4*arr(n-2) + 3*arr(n-3) - 4.0d0/3*arr(n-4) + 1.0d0/4* arr(n-5)
      do i = 3,n-2
        compder(i) = 1.0d0/12*arr(i-2) - 2.0d0/3*arr(i-1) + 2.0d0/3*arr(i+1) - 1.0d0/12*arr(i+2)
      enddo

      compder = compder/dy
    endfunction compder
  
    subroutine specder(xmin,xmax,n,input,deriv)  !takes spectral derivatives of order n of a function evaluated at n points.
      use, intrinsic :: iso_c_binding
      implicit none
      include '/usr/local/include/fftw3.f03'
      type(c_ptr) :: plan
      integer :: n,i
      integer *8::plan_forward,plan_backward
      complex(c_double_complex), dimension(n) :: input,output,input2,output2
      real *8,dimension(n)::x,k,der,deriv
      real *8::xmin,xmax,dx,pi
      pi = 4.0d0*atan(1.0d0)

      call dfftw_plan_dft_1d_(plan_forward,n, input,output,fftw_forward,fftw_estimate)
      call dfftw_execute_(plan_forward)
      call dfftw_destroy_plan_(plan_forward)

      do i = 1,n/2
      k(i) = i-1.0d0
      end do

      k(n/2+1) = 0.0d0

      do i = n/2+2,n
      k(i) = (-1.0d0)*n+i-1.0d0
      end do

      do i=1,n
        input2(i) =2.0d0*pi/(xmax-xmin)*k(i)*cmplx(0.0D0,1.0D0,kind=16)*output(i)/n
      end do

      call dfftw_plan_dft_1d_(plan_backward,n, input2, output2,fftw_backward,fftw_estimate)
      call dfftw_execute_(plan_backward, input2, output2)
      call dfftw_destroy_plan_(plan_backward)

      do i = 1,n
        deriv(i) = real(output2(i))
      end do
    end subroutine specder

endmodule solvefp
