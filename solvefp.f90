module solvefp
contains
    
    

    function uprime(z)
        implicit none
        real *8::z,lsmall,lbig,l,umax,uprime
        lsmall = 1
        lbig = 9
        umax = 1

        l = lsmall + lbig
        if mod(z,l)<lsmall then
            uprime = umax/lsmall
        else
            uprime = -umax/lbig
        endif
    endfunction uprime

    function avec(x,y,k,gm)
      implicit none
      real *8:: x,y,k,up,um,gm,fx,fy,avec(2)
      um = uprime(y-x/2)
      up = uprime(y+x/2)
      fx = -k*x-0.5*(up-um)
      fy = -up-um
      avec(1) = 2*fx/gm
      avec(2) = fy/(2*gm)
    endfunction avec

    function barray(tau,tbar,gm)
      implicit none
      real *8::tau,tbar,gm,barray(2,2)
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
