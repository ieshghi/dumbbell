program runval
    use solveval
implicit none
integer::n,mult,info,i
real *8::sig
real *8,allocatable::lhs(:,:),diag(:),work(:)

n = 10000
sig = 1
mult = 3
allocate(lhs(n,n),diag(n),work(3*n-1))

lhs = genlhs(n,sig,mult)
call dsyev('N','L',n,lhs,n,diag,work,3*n-1,info)
write(*,*) 'info = ',info

open(1,file = 'fpres/evals.dat')

do i = 1,n
    write(1,*) diag(n) 
enddo

close(1)

endprogram runval
