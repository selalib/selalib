!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: gausslobatto
!
!> @brief Gauss-Lobatto interpolation tools
!> @details Here are several of the Gauss-Lobatto tools :
!>    - Gauss-Lobatto points and weight,
!>    - Gauss-Lobatto bases functions and the integral of their product,
!>    - integral of product of Gauss-Lobatto function and their derivative.
!>
!>    To use this module you must also link to the compilation gauss.f and lob.f
!>
!>    The mass matrix (which is the integral of \f$ \phi_i \times \phi_j \f$) is simply 
!>    diag(weight), so there is no need to store it more than just the weight.
!>
!>    We also need the derivative matrix D.
!>    \f[ D_{i,j}=\int \pih_i \phi^'_j \f]
!>
!>    This module will first be limited to 1D and should extend as people will
!>    have the need for higher dimension (and so have time to write it).
!>         
!> @author Madaule Eric
!>
module sll_gausslobatto
#include "sll_working_precision.h"

  use sll_constants
  
  implicit none

  !> Gauss-Lobatto nodes and weigh on a reference element [-1;1] in 1D.
  !> This also includes the degree of polynomials and the matrix of derivatives
  !> \f[ D_{i,j}=\int(\Psi_i\Psi'_j) = 
  !>     w_i\sum_{l/=j}1/(x_j-x_l)\prod_{m/=j,m/=l}(x_i-x_m)/(x_j-x_m) \f]
  !> 
  !>  A gausslobatto1d object contains node, weight, jac and degree.
  !>  These are allocatable array. Use the construction subroutine
  !>  to build it. Only degree is a scalar. It is the degree of corresponding
  !>  polynomials. Only jac must be filled manually (but it is allocated in
  !>  the constructor)
  type gausslobatto1D
     sll_real64,dimension(:),pointer :: node,weigh
     sll_int32 :: degree
     sll_real64,dimension(:,:),pointer :: der
  end type gausslobatto1D

  interface delete
     module procedure delete_gausslobatto_1D
  end interface delete

  abstract interface
     function function_1D(x)
       use sll_working_precision ! can't pass a header file because the
                                 ! preprocessor prevents double inclusion.
                                 ! This is very rare.
       sll_real64             :: function_1D
       sll_real64, intent(in) :: x
     end function function_1D
  end interface


contains

!> @brief Construction of Gauss-Lobatto nodes and weights
!> @details Construction of Gauss-Lobatto nodes and weights.
!>          This routine fill the node and weight but you have 
!>          to fill the Jacobian with
!>          the transformation you consider
!>          This should be extended as I will complete the type.
!>          I can't deal very well with fortran 77 interface. To change between simple
!>          and double precision you have to go into the file and change some comments
!>          (very easy, it is explained in the file) and compile again (sorry)
  subroutine init_gausslobatto_1d(size,gl_obj)

    sll_int32,intent(in) :: size !< size number of Gauss-Lobatto points, 
                                 !< should be at least 2 (if 1,
                                 !< then piecewise constant => Gauss-Lobatto 
                                 !< is impossible)
    type(gausslobatto1D),intent(out) :: gl_obj !< Gauss-Lobatto object to build

    sll_int32 :: err,i
    sll_real64,dimension(0:size-1) :: alpha,beta
    sll_real64,dimension(size) :: e,a,b

    if (size<=1) then
       print*,"not enought points to build the Gauss-Lobatto interpolator"
       print*,"exiting..."
       stop
    end if

    ALLOCATE(gl_obj%node(size))
    ALLOCATE(gl_obj%weigh(size))
    ALLOCATE(gl_obj%der(size,size))

    gl_obj%node=0.0d0
    gl_obj%weigh=0.0d0
    gl_obj%degree=size-1

    do i=0,size-1
       alpha(i)=0.0d0
       beta(i)=real(i,kind(gl_obj%degree))**2/((2.0d0*i+1.0d0)*(2.0d0*i-1.0d0))
    end do
    !for Gauss-Legendre and Gauss-Lobatto, beta(0)=int(dlambda)
    !see Algorithm xxx - ORTHPOL: A package of routines for  generating orthogonal
    !polynomials and Gauss-type quadrature rules by _Walter Gautschi_
    beta(0)=2.0d0

    !for single precision, comment the first line and uncomment the second
    !for double precision, comment the second line and uncomment the first
    call dlob(size-2,alpha,beta,-1.0d0,1.0d0,gl_obj%node,gl_obj%weigh,err,e,a,b)
    !call lob(size-2,alpha,beta,-1.0,1.0,gl_obj%node,gl_obj%weigh,err,e,a,b)

    call derivative_matrix_1d(gl_obj)

  end subroutine init_gausslobatto_1d

  !> @brief delete an object of type gausslobatto1D
  !> @details delete all array in an object of type gausslobatto1D
  !> @param[INOUT] gl_obj the object to delete
  subroutine delete_gausslobatto_1d(gl_obj)

    type(gausslobatto1D),intent(inout) :: gl_obj

    DEALLOCATE(gl_obj%node)
    DEALLOCATE(gl_obj%weigh)
    DEALLOCATE(gl_obj%der)

  end subroutine delete_gausslobatto_1d

  !> Construction of the derivative matrix for Gauss-Lobatto 1D,
  !> The matrix must be already allocated of size (number of point\f$)^2\f$.
  !> \f[  der(i,j)=int(\Phi_i.\Phi'_j)_{[-1;1]^2}
  !>                  =w_i.\Phi'_j(x_i) \f]
  !> @param[INOUT] gl_obj gausslobatto1D object to build derivative
  subroutine derivative_matrix_1d(gl_obj)

    type(gausslobatto1D),intent(inout) :: gl_obj

    sll_int32 :: nb_pts,i,j,l,m
    sll_real64 :: prod

    nb_pts=gl_obj%degree+1

    gl_obj%der=0.0d0

    !loop on all element of D
    !loop on columns
    do j=1,nb_pts
       !loop on rows
       do i=1,nb_pts
          !loop on all the derivatives
          !the code is writen so there is no if
          do l=1,j-1
             prod=1.0d0
             do m=1,l-1!min(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=l+1,j-1!min(j,l)+1,max(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=j+1,nb_pts!max(j,l)+1,nb_pts
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             prod=prod/(gl_obj%node(j)-gl_obj%node(l))
             gl_obj%der(i,j)=gl_obj%der(i,j)+prod
          end do
          do l=j+1,nb_pts
             prod=1.0d0
             do m=1,j-1!min(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=j+1,l-1!min(j,l)+1,max(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=l+1,nb_pts!max(j,l)+1,nb_pts
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             prod=prod/(gl_obj%node(j)-gl_obj%node(l))
             gl_obj%der(i,j)=gl_obj%der(i,j)+prod
          end do
          gl_obj%der(i,j)=gl_obj%der(i,j)*gl_obj%weigh(i)
       end do
    end do

  end subroutine derivative_matrix_1d

  !> @brief Gauss-Lobatto Quadrature.
  !> @details To integrate the function \f$ f(x) \f$
  !> (real-valued and of a single, real-valued argument x)
  !> over the interval \f$ [a,b] \f$, we use the Gauss-Lobatto formula 
  !> \f[ \int_{-1}^1 f(x)dx \approx \sum_{k=1}^{n} w_k f(x_k) \f]
  !> where n represents the desired number of Gauss points.
  !>
  !> the function maps the interval \f$ [-1,1] \f$ into the
  !> arbitrary interval \f$ [a,b] \f$.
  !>
  !> To be considered is to split this function into degree-specific
  !> functions to avoid the select statement.
  !> @param f the function to be integrated
  !> @param[in] a left-bound of the definition interval of f  
  !> @param[in] b right-bound of the definition interval of f 
  !> @param[in] n the desired number of Gauss points
  !> @return The value of the integral
  function gauss_lobatto_integral_1D( f, a, b, n )
    sll_real64                :: gauss_lobatto_integral_1D
    procedure(function_1D)    :: f
    sll_real64, intent(in)    :: a
    sll_real64, intent(in)    :: b
    sll_int32,  intent(in)    :: n 
    sll_real64, dimension(10) :: xk, alpha
    sll_real64, dimension(10) :: wk, beta
    sll_int32                 :: k
    sll_real64                :: x
    sll_real64                :: c1
    sll_real64                :: c2
    sll_real64                :: ans

    xk(:) = 0.0_f64
    wk(:) = 0.0_f64
    c1 = 0.5_f64*(b-a)
    c2 = 0.5_f64*(b+a)

    alpha = 0.0d0
    do i = 0, 9
    beta(i+1)=i*i/((2.0d0*i+1)*(2.0d0*i-1))
    end do

    !call dlob(size-2,alpha,beta,-1.0d0,1.0d0,gl_obj%node,gl_obj%weigh,err,e,a,b)

    do k=1,n
       x = c1*xk(k) + c2
       ans = ans + f(x)*wk(k)
    end do

    gauss_lobatto_integral_1D = c1*ans


  end function gauss_lobatto_integral_1D



!> This comes from http://dl.acm.org, Algorithme 726 : ORTHPOL, appendices and supplements

!> To use those functions, READ the documentation beside and find more information 
!> about coefficients in paper *Algorithm xxx - ORTHPOL: A package of routines for 
!> generating orthogonal polynomials and Gauss-type quadrature rules* by _Walter 
!> Gautschi_ (here xxx is 726 in other references) formulas (1.1) to (1.3) page 2, 
!> and book **Numerical Mathematics** by _Alfio Quarteroni_, _Riccardo Sacco_ and 
!> _Fausto Saleri_ section 10.
!> If you just want to use Gauss-Lobatto inside Selalib, just use what is in 
!> sll_gausslobatto.F90 ans see selalib documentation
!>
!> Given  n  and a measure  dlambda, this routine generates the 
!> (n+2)-point Gauss-Lobatto quadrature formula
!> 
!>   integral over supp(dlambda) of f(x)dlambda(x)
!> 
!>      = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k))  
!> 
!>              + w(n+1)f(x(n+1)) + R(n;f).
!> 
!> The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
!> =w(k), k=0,1,...,n,n+1. The user has to supply the recursion
!> coefficients  alpha(k), beta(k), k=0,1,...,n,n+1, for the measure
!> dlambda. The nodes and weights are computed in terms of the
!> eigenvalues and first component of the normalized eigenvectors of
!> a slightly modified Jacobi matrix of order  n+2. The routine calls 
!> upon the subroutine  gauss  and the function subroutine  r1mach.
!> 
!>   Input:  
!>     - n : the number of interior points in the Gauss-Lobatto
!>                  formula; type integer
!>     - alpha,beta : arrays of dimension  n+2  to be supplied with
!>                  the recursion coefficients  alpha(k-1), beta(k-1),
!>                  k=1,2,...,n+2, of the underlying measure; the
!>                  routine does not use  alpha(n+2), beta(n+2)
!>     - aleft,right : the prescribed left and right endpoints 
!>                  x(0)  and  x(n+1)  of the Gauss-Lobatto formula
!> 
!>   Output: 
!>     - zero : an array of dimension  n+2  containing the nodes (in 
!>                  increasing order)  zero(k)=x(k), k=0,1,...,n,n+1
!>     - weight : an array of dimension  n+2  containing the weights 
!>                  weight(k)=w(k), k=0,1,...,n,n+1
!>     - ierr : an error flag inherited from the routine  gauss
!> 
!> The arrays  e,a,b  are needed for working space.
!> 
subroutine dlob(n,dalpha,dbeta,dleft,dright,dzero,dweigh,ierr,de,da,db)

sll_int32  :: n, np1, np2, k, ierr
sll_real64 :: dleft,dright,depsma,dp0l,dp0r,dp1l,dp1r,dpm1l
sll_real64 :: dpm1r,ddet,dalpha(*),dbeta(*),dzero(*),dweigh(*),de(*),da(*),db(*)

! 
! The arrays  dalpha,dbeta,dzero,dweigh,de,da,db  are assumed to have
! dimension  n+2.
! 
depsma=epsilon(1.0d0)
! 
! depsma is the machine double precision.
! 
np1=n+1
np2=n+2
do k=1,np2
   da(k)=dalpha(k)
   db(k)=dbeta(k)
end do

dp0l=0.d0
dp0r=0.d0
dp1l=1.d0
dp1r=1.d0
do k=1,np1
  dpm1l=dp0l
  dp0l=dp1l
  dpm1r=dp0r
  dp0r=dp1r
  dp1l=(dleft-da(k))*dp0l-db(k)*dpm1l
  dp1r=(dright-da(k))*dp0r-db(k)*dpm1r
end do

ddet=dp1l*dp0r-dp1r*dp0l
da(np2)=(dleft*dp1l*dp0r-dright*dp1r*dp0l)/ddet
db(np2)=(dright-dleft)*dp1l*dp1r/ddet
call dgauss(np2,da,db,depsma,dzero,dweigh,ierr,de)
return
end subroutine dlob


!!$ This comes from http://dl.acm.org, Algorithme 726 : ORTHPOL, appendices and supplements

!!$ To use those functions, READ the documentation beside and find more information 
!!$ about coefficients in paper *Algorithm xxx - ORTHPOL: A package of routines for 
!!$ generating orthogonal polynomials and Gauss-type quadrature rules* by _Walter 
!!$ Gautschi_ (here xxx is 726 in other references) formulas (1.1) to (1.3) page 2, 
!!$ and book **Numerical Mathematics** by _Alfio Quarteroni_, _Riccardo Sacco_ and 
!!$ _Fausto Saleri_ section 10.
!!$ If you just want to use Gauss-Lobatto inside Selalib, just use what is in 
!!$ sll_gausslobatto.F90 ans see selalib documentation

!!$ 
!!$ Given  n  and a measure  dlambda, this routine generates the n-point
!!$ Gaussian quadrature formula
!!$  
!!$     integral over supp(dlambda) of f(x)dlambda(x)
!!$ 
!!$        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
!!$ 
!!$ The nodes are returned as  zero(k)=x(k) and the weights as 
!!$ weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion 
!!$ coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure 
!!$ dlambda. The routine computes the nodes as eigenvalues, and the 
!!$ weights in term of the first component of the respective normalized 
!!$ eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
!!$ It uses a translation and adaptation of the algol procedure  imtql2,
!!$ Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified 
!!$ by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for 
!!$ Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
!!$ routine  imtql2.
!!$ 
!!$        Input:  n - - the number of points in the Gaussian quadrature	
!!$                      formula; type integer
!!$                alpha,beta - - arrays of dimension  n  to be filled 
!!$                      with the values of  alpha(k-1), beta(k-1), k=1,2,
!!$                      ...,n
!!$                eps - the relative accuracy desired in the nodes
!!$                      and weights
!!$ 
!!$        Output: zero- array of dimension  n  containing the Gaussian 
!!$                      nodes (in increasing order)  zero(k)=x(k), k=1,2,
!!$                      ...,n
!!$                weight - array of dimension  n  containing the 
!!$                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
!!$                ierr- an error flag equal to  0  on normal return,
!!$                      equal to  i  if the QR algorithm does not
!!$                      converge within 30 iterations on evaluating the 
!!$                      i-th eigenvalue, equal to  -1  if  n  is not in
!!$                      range, and equal to  -2  if one of the beta's is 
!!$                      negative.
!!$ 

subroutine dgauss(n,dalpha,dbeta,deps,dzero,dweigh,ierr,de)

sll_int32  :: n, k, l, i, m, ii, j, mml, ierr
sll_real64 :: deps,dp,dg,dr,ds,dc,df,db
sll_real64 :: dalpha(n),dbeta(n),dzero(n),dweigh(n),de(n)

if(n.lt.1) then
  ierr=-1
  return
end if
ierr=0
dzero(1)=dalpha(1)
if(dbeta(1).lt.0.d0) then
  ierr=-2
  return
end if
dweigh(1)=dbeta(1)
if (n.eq.1) return
dweigh(1)=1.d0
de(n)=0.d0

do k=2,n
  dzero(k)=dalpha(k)
  if(dbeta(k).lt.0.d0) then
    ierr=-2
    return
  end if
  de(k-1)=dsqrt(dbeta(k))
  dweigh(k)=0.d0
end do

do l=1,n
   j=0
   105 do m=l,n
        if(m.eq.n) exit
        if(dabs(de(m)).le.deps*(dabs(dzero(m))+dabs(dzero(m+1)))) exit
   end do
   dp=dzero(l)
   if(m.eq.l) cycle
   if(j.eq.30) goto 400
   j=j+1
   dg=(dzero(l+1)-dp)/(2.d0*de(l))
   dr=dsqrt(dg*dg+1.d0)
   dg=dzero(m)-dp+de(l)/(dg+dsign(dr,dg))
   ds=1.d0
   dc=1.d0
   dp=0.d0
   mml=m-l
   do ii=1,mml
      i=m-ii
      df=ds*de(i)
      db=dc*de(i)

      if(dabs(df) .ge. dabs(dg)) then
         dc=dg/df
         dr=dsqrt(dc*dc+1.d0)
         de(i+1)=df*dr
         ds=1.d0/dr
         dc=dc*ds
      else
         ds=df/dg
         dr=dsqrt(ds*ds+1.d0)
         de(i+1)=dg*dr
         dc=1.d0/dr
         ds=ds*dc
      end if
      
      dg=dzero(i+1)-dp
      dr=(dzero(i)-dg)*ds+2.d0*dc*db
      dp=ds*dr
      dzero(i+1)=dg+dp
      dg=dc*dr-db
      df=dweigh(i+1)
      dweigh(i+1)=ds*dweigh(i)+dc*df
      dweigh(i)=dc*dweigh(i)-ds*df

   end do

   dzero(l)=dzero(l)-dp
   de(l)=dg
   de(m)=0.d0
   goto 105

end do

do ii=2,n

   i=ii-1
   k=i
   dp=dzero(i)

   do j=ii,n
     if(dzero(j).ge.dp) cycle
     k=j
     dp=dzero(j)
   end do

   if(k.ne.i) then
      dzero(k)=dzero(i)
      dzero(i)=dp
      dp=dweigh(i)
      dweigh(i)=dweigh(k)
      dweigh(k)=dp
   end if

end do

do k=1,n
  dweigh(k)=dbeta(1)*dweigh(k)*dweigh(k)
end do

return
  400 ierr=l
return
end subroutine

end module sll_gausslobatto
