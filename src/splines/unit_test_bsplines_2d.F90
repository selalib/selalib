program test_bsplines_2d

integer :: i,iflag,j,jj,kx,ky,lefty,mflag,nx,ny
integer, parameter :: nx=7,kx=3,ny=6,ky=4
real(8) ::  bcoef(nx,ny),taux(nx),tauy(ny),tx(nx+kx),ty(ny+ky)
real(8) ::  work1(nx,ny),work2(nx),work3(nx*ny)

!     integer i,iflag,j,jj,kp1,kx,ky,lefty,mflag,nx,ny
!     data nx,kx,ny,ky /7,3,6,4/    
!     real bcoef(7,6),taux(7),tauy(6),tx(10),ty(10)
!    *         ,work1(7,6),work2(7),work3(42)

g(x,y) = amax1(x-3.5,0.)**2 + amax1(y-3.,0.)**3

! *** set up data points and knots
!     in x, interpolate between knots by parabolic splines, using
!     not-a-knot end condition

do i=1,nx
  taux(i) = float(i)
end do
do i=1,kx
  tx(i) = taux(1)
  tx(nx+i) = taux(nx)
end do

kp1 = kx+1  
do i=kp1,nx
  tx(i) = (taux(i-kx+1) + taux(i-kx+2))/2.
end do

!     in y, interpolate at knots by cubic splines, using not-a-knot
!     end condition

do i=1,ny
  tauy(i) = float(i)
end do
do i=1,ky
  ty(i) = tauy(1)
  ty(ny+i) = tauy(ny)
end do

kp1 = ky+1  
do i=kp1,ny
  ty(i) = tauy(i-ky+2)
end do

!  *** generate and print out function values
print 620,(tauy(i),i=1,ny)
620 format(' given data'//6f12.1)
do i=1,nx
  do j=1,ny
    bcoef(i,j) = g(taux(i),tauy(j))
    print 632,taux(i),(bcoef(i,j),j=1,ny)
  end do
end do
632 format(f5.1,6e12.5)
!
!  *** construct b-coefficients of interpolant
!
call spli2d(taux,bcoef,tx,nx,kx,ny,work2,work3,work1,iflag)
call spli2d(tauy,work1,ty,ny,ky,nx,work2,work3,bcoef,iflag)

!  *** evaluate interpolation error at mesh points and print out
print 640,(tauy(j),j=1,ny)
640 format(//' interpolation error'//6f12.1)
do j=1,ny
  call interv(ty,ny+1,tauy(j),lefty,mflag)
  do i=1,nx
    do jj=1,ky
      work2(jj)=bvalue(tx,bcoef(1,lefty-ky+jj),nx,kx,taux(i),0)
    end do
    work1(i,j) = g(taux(i),tauy(j)) - &
                 bvalue(ty(lefty-ky+1),work2,ky,ky,tauy(j),0)
  end do
end do
do i=1,nx
  print 632,taux(i),(work1(i,j),j=1,ny)
end do

contains

subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )********    
!  from  * a practical guide to splines *  by c. de boor    
!alls bsplvb, banfac/slv
!  this is an extended version of  splint , for the use in tensor prod- --------    
!  uct interpolation.                                                   --------    
!     
!   spli2d  produces the b-spline coeff.s  bcoef(j,.)  of the spline of --------    
!   order  k  with knots  t (i), i=1,..., n + k , which takes on the    --------    
!   value  gtau (i,j)  at  tau (i), i=1,..., n , j=1,..., m .           --------    
!     
!******  i n p u t  ******    
!  tau.....array of length  n , containing data point abscissae.  
!    a s s u m p t i o n . . .  tau  is strictly increasing 
!  gtau(.,j)..corresponding array of length  n , containing data point  --------    
!        ordinates, j=1,...,m                                           --------    
!  t.....knot sequence, of length  n+k    
!  n.....number of data points and dimension of spline space  s(k,t)    
!  k.....order of spline
!  m.....number of data sets                                            ********    
!     
!******  w o r k   a r e a  ******                                      ********    
!  work  a vector of length  n                                          ********    
!     
!******  o u t p u t  ******  
!  q.....array of size  (2*k-1)*n , containing the triangular factoriz- 
!        ation of the coefficient matrix of the linear system for the b-
!        coefficients of the spline interpolant.
!           the b-coeffs for the interpolant of an additional data set  
!        (tau(i),htau(i)), i=1,...,n  with the same data abscissae can  
!        be obtained without going through all the calculations in this 
!        routine, simply by loading  htau  into  bcoef  and then execut-
!        ing the    call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )  
!  bcoef.....the b-coefficients of the interpolant, of length  n  
!  iflag.....an integer indicating success (= 1)  or failure (= 2)
!        the linear system to be solved is (theoretically) invertible if
!        and only if    
!              t(i) .lt. tau(i) .lt. tau(i+k),    all i.    
!        violation of this condition is certain to lead to  iflag = 2 . 
!     
!******  m e t h o d  ******  
!     the i-th equation of the linear system  a*bcoef = b  for the b-co-
!  effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
!  hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1   
!   bands (if it is invertible).    
!     the matrix  a  is generated row by row and stored, diagonal by di-
!  agonal, in the  r o w s  of the array  q , with the main diagonal go-
!  ing into row  k .  see comments in the program below.    
!     the banded system is then solved by a call to  banfac (which con- 
!  structs the triangular factorization for  a  and stores it again in  
!   q ), followed by a call to  banslv (which then obtains the solution 
!   bcoef  by substitution).  
!     banfac  does no pivoting, since the total positivity of the matrix
!  a  makes this unnecessary. 
!     
integer ::  iflag,k,m,n,   i,ilp1mx,j,jj,km1,kpkm2,left,lenq,np1      ********    
real(8) ::  bcoef(m,n),gtau(n,m),q(1),t(1),tau(n),work(n),   taui        ********    
!     dimension q(2*k-1,n), t(n+k)  
!urrent fortran standard makes it impossible to specify precisely the   
!  dimension of  q  and  t  without the introduction of otherwise super-
!  fluous additional arguments.     
      np1 = n + 1 
      km1 = k - 1 
      kpkm2 = 2*km1     
      left = k    
c                zero out all entries of q
      lenq = n*(k+km1)  
      do 5 i=1,lenq     
    5    q(i) = 0.
c     
c  ***   loop over i to construct the  n  interpolation equations 
      do 30 i=1,n 
         taui = tau(i)  
         ilp1mx = min0(i+k,np1)     
c        *** find  left  in the closed interval (i,i+k-1) such that     
c                t(left) .le. tau(i) .lt. t(left+1)   
c        matrix is singular if this is not possible   
         left = max0(left,i)  
         if (taui .lt. t(left))         go to 998     
   15       if (taui .lt. t(left+1))    go to 16
            left = left + 1   
            if (left .lt. ilp1mx)       go to 15
         left = left - 1
         if (taui .gt. t(left+1))       go to 998     
c        *** the i-th equation enforces interpolation at taui, hence    
c        a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j = 
c        left-k+1,...,left actually might be nonzero. these  k  numbers 
c        are returned, in  work  (used for temp.storage here), by the   --------    
c        following
   16    call bsplvb ( t, k, 1, taui, left, work )                      --------    
c        we therefore want  work(j) = b(left -k+j)(taui) to go into     --------    
c        a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since  
c        a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q 
c        as a two-dim. array , with  2*k-1  rows (see comments in 
c        banfac). in the present program, we treat  q  as an equivalent 
c        one-dimensional array (because of fortran restrictions on
c        dimension statements) . we therefore want  work(j) to go into  --------    
c        entry    
c            i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)   
c                   =  i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j   
c        of  q .  
         jj = i-left+1 + (left-k)*(k+km1) 
         do 30 j=1,k    
            jj = jj+kpkm2     
   30       q(jj) = work(j)                                             --------    
c     
c     ***obtain factorization of  a  , stored again in  q.  
      call banfac ( q, k+km1, n, km1, km1, iflag )    
                                        go to (40,999), iflag     
c     *** solve  a*bcoef = gtau  by backsubstitution  
   40 do 50 j=1,m                                                       ********    
         do 41 i=1,n                                                    --------    
   41       work(i) = gtau(i,j)                                         ********    
         call banslv ( q, k+km1, n, km1, km1, work )                    --------    
         do 50 i=1,n                                                    ********    
   50       bcoef(j,i) = work(i)                                        ********    
                                        return  
  998 iflag = 2   
  999 print 699   
  699 format(41h linear system in  splint  not invertible)  
                                        return  
      end   

end program test_bsplines_2d
