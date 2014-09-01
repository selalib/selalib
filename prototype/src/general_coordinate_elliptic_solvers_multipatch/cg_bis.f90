! -------------------------- MODULE cg.f90 ----------------------------

!************************************************************************
!*                                                                      *
!* Conjugate Gradient Method (CG Method)                                *
!* -------------------------------------                                *
!*                                                                      *
!* Programming language: ANSI C                                         *
!* Compiler:             Turbo C 2.0                                    *
!* Computer:             IBM PS/2 70 with 80387                         *
!* Sources:              [BUNS85], [SCHW], [MAES84]                     *
!* Author:               Juergen Dietel, Computer Center, RWTH Aachen   *
!* Date:                 7.31.1992                                      *
!*                                                                      *
!*             F90 version by J-P Moreau (without dynamic allocations). *
!*                                (www.jpmoreau.fr)                     *
!************************************************************************
Subroutine cg_method (    &     ! Conjugate Gradient Method
                       n, &     ! Size of the linear system
                       a, &     ! System matrix
                       y, &     ! right hand side
                       x, &     ! solution vector 
                       fehler & ! error code
                     )
! original name: cg_verfahren()
integer, parameter :: SIZE=24
real*8, parameter :: ZERO=0.d0, MACH_EPS=2.d-16
real*8  a(0:SIZE,0:SIZE),x(0:SIZE),y(0:SIZE)
integer fehler
!************************************************************************
!* CG solves the linear system                                          *
!*                         A * X = Y                                    *
!* for a symmetric, positive definite matrix A via the conjugate        *
!* gradient method.                                                     *
!*                                                                      *
!* Input parameters:                                                    *
!* =================                                                    *
!* n  Size of the linear system                                         *
!* a  [0..n-1,0..n-1] system matrix A. Only the upper triangle of A is  *
!*    used.                                                             *
!* y  [0..n-1] vector of the right hand side                            *
!*                                                                      *
!* Output parameters:                                                   *
!* ==================                                                   *
!* x  [0..n-1] vector giving the solution                               *
!*                                                                      *
!* Return value:                                                        *
!* =============                                                        *
!* = 0: all is ok                                                       *
!* = 1: n < 2 or other disallowed input parameters                      *
!* = 2: memory exceeded                                                 *
!*                                                                      *
!************************************************************************
  real*8 d(0:n-1), &   ! (0..n-1) auxiliary vectors d and g
         g(0:n-1), &
         AmalD(0:n-1)  ! (0..n-1) auxiliary vector A * d
  real*8 alpha,    &   ! coefficient
         beta,     &   ! coefficient
         dividend, &   ! numerator and denominator of a fraction
         divisor,  &   ! respectively, used to compute alpha, beta
         hilf,     &   ! auxiliary variables
         hilf2,    &
         abstand,  &   ! distance of two successive approximations
                       ! for the solution vector x (taken in the
                       ! euclidean norm)
         xnorm         ! euklidean norm of x
  integer k, i, j      ! loop variables

  if (n < 2) then      ! invalid parameter?
    fehler=1
	return
  end if

  !------------------------------------------------------------------
  ! start with x at the origin                                       
  !------------------------------------------------------------------
  do i = n - 1, 0, -1
    x(i) = ZERO
  end do

  !------------------------------------------------------------------
  ! initialize  d and g :                                            
  ! d = -g = -(a*x - y) = y (since x = 0)                            
  !------------------------------------------------------------------
  do i = n - 1, 0, -1
    hilf = y(i)
    d(i) = hilf
    g(i) = -hilf
  end do


  !------------------------------------------------------------------
  ! perform at most n steps of the CG Method                         
  !------------------------------------------------------------------
  do k = n, 0, -1

    !----------------------------------------------------------------
    ! compute new alpha:                                             
    ! alpha = -(d(transp) * g) / (d(transp) * (a * d))               
    !----------------------------------------------------------------

    dividend = ZERO
    divisor  = ZERO

    do i = n - 1, 0, -1
      dividend = dividend + d(i) * g(i)
	  hilf = ZERO
      do j = 0, i-1
        hilf = hilf + a(j,i) * d(j)
      end do
      do j = i, n-1
        hilf = hilf + a(i,j) * d(j)
      end do
      AmalD(i) = hilf
      divisor = divisor + d(i) * hilf
    end do

    if (divisor.eq.ZERO) then
	  fehler=0
      return
    end if

    alpha = -dividend / divisor

    !----------------------------------------------------------------
    ! compute the norm of x und  alpha * d  and find a new x:        
    ! x  =  x + alpha * d, then check whether x is close enough,     
    ! in order to stop the process before n complete steps           
    !----------------------------------------------------------------
    xnorm   = ZERO
    abstand = ZERO

    do i = n - 1, 0, -1
      hilf =  x(i)
      xnorm   = xnorm + hilf*hilf
      hilf2   =  alpha * d(i)
      abstand = abstand + hilf2*hilf2
      x(i)    =  hilf + hilf2
    end do

    if (abstand < MACH_EPS * xnorm) then
	  fehler=0
      return
    end if


    !----------------------------------------------------------------
    ! compute new g:   g  =  g + alpha * (a * d)                     
    !----------------------------------------------------------------
    do i = n - 1, 0, -1
      g(i) = g(i) + alpha * AmalD(i)
    end do

    !----------------------------------------------------------------
    ! compute new beta :                                             
    ! beta = (g(transp) * (a * d)) / (d(transp) * (a * d))           
    !----------------------------------------------------------------
    dividend = ZERO
    do i = n - 1, 0, -1
      dividend = dividend + g(i) * AmalD(i)
    end do

    beta = dividend / divisor

    !----------------------------------------------------------------
    ! compute new d :   d  =  - g + beta * d                         
    !----------------------------------------------------------------
    do i = n - 1, 0, -1
      d(i) = -g(i) + beta * d(i)
    end do

  end do  !k loop

  fehler=0
  return
end

! ---------------------------- END cg.cpp -------------------------- 
