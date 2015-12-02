!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @ingroup splines 
!> @brief  
!> provides capabilities for data and derivative
!> interpolation with quintic splines 
!> @details
!> inspl5 and spln5 routines come from
!> "An algorithm for the interpolation of functions using quintic splines"
!> by E.H. Mund, P. Hallet and J.P. Hennart
!> Journal of COmputational and Applied Mathematics, volume 1, number 4, 1975.
!>
!> Periodic boundary conditions are not implemented. You must set function value
!> and its derivative a the boundary
module sll_m_quintic_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none

contains

!> Calculation of the parameters of an interpolating quintic splines
subroutine inspl5(n,x,ind1,indn,cf,h)

sll_int32,  intent(in)    :: n         !< number of interpolation points
sll_real64, intent(in)    :: x(n)      !< vector of abscissae
sll_real64, intent(inout) :: cf(1:3,n) !< ordinates, first and second derivatives
sll_int32,  intent(in)    :: ind1      !< boundary conditions switches at x=1.
sll_int32,  intent(in)    :: indn      !< boundary conditions switches at x=n.
                                       !< for both switches one 
                                       !< has the correspondance
                                       !< = 1 type 1
                                       !< = 0 type 2
                                       !< =-1 type 3
sll_real64, intent(out)   :: h(6*n-3)  !< auxiliary vector

sll_real64 :: a, as, det, di1, di2, di3, di4
sll_real64 :: dp1, dp2, dp3, dp4, dpd1, dpd2
sll_real64 :: ds1, ds2, ds3, ds4
sll_real64 :: fp, fp3, fp4
sll_real64 :: p1, p2, p3, q1, q2
sll_real64 :: sf1, sf2
sll_int32  :: i, im, k


h = 0.0_f64

as  = 64.0_f64/ 3.0_f64
a   = sqrt(as)
p1  = 0.0_f64
p2  = 0.0_f64
p3  = 0.0_f64
q1  = 0.0_f64
q2  = 0.0_f64
sf1 = 0.0_f64
sf2 = 0.0_f64
im  = 0

!Forward elimination procedure for Choleskys algorithm
!applied to block-tridiagonal systems of algebraic equations.

i = 1
goto 15

 5 continue

   di1 = ds1
   di2 = ds3
   di3 = ds2
   di4 = ds4

10 continue

   sf1 = fp4
   sf2 = fp3*q1*a

15 continue

   dp1 = 64.0_f64*p3
   dp2 = 12.0_f64*p2*q1*a
   dp3 = dp2
   dp4 = 3.0_f64*p1*q2*as

   if (i >= n) then !boundary condition at xn

     if (indn < 0) then

       dp2 = 0.0_f64
       dp4 = 1.0_f64
       di2 = 0.0_f64
       di4 = 0.0_f64
       sf2 = cf(3,n)*(x(n)-x(n-1))/a

     else if (indn > 0) then

       dp1 = 1.0_f64
       dp2 = 0.0_f64
       dp3 = 0.0_f64
       dp4 = 1.0_f64
       di1 = 0.0_f64
       di2 = 0.0_f64
       di3 = 0.0_f64
       di4 = 0.0_f64
       sf1 = -cf(2,n)
       sf2 =  cf(3,n)*(x(n)-x(n-1))/a

     endif

     goto 50

   end if

   p1 = 1.0_f64/(x(i+1)-x(i))
   p2 = p1*p1
   p3 = p1*p2
   q1 = p1
   q2 = p2

   if (i < n-1) then
     q1 = 1.0_f64 / (x(i+2)-x(i+1))
     q2 = q1*q1
   end if

   fp = cf(1,i)-cf(1,i+1)
   fp3= 20.0_f64*p3*fp
   fp4=  6.0_f64*p1*fp3
   dp1= dp1+64.0_f64*p3
   dp2= dp2-12.0_f64*p3*a
   dp3= dp2
   dp4= dp4+3.0_f64*p3*as
   ds1=-56.0_f64*p3
   ds2= 8.0_f64*p3*a
   ds3= -8.0_f64*p2*q1*a
   ds4=    p2*q1*as
   sf1=sf1+fp4
   sf2=sf2-fp3*p1*a

   if (i>1) goto 50

   di1 = ds1
   di2 = ds3
   di3 = ds2
   di4 = ds4

!Boundary condition at x1

   if (ind1 < 0) then

     dp2 = 0.0_f64
     dp4 = 1.0_f64
     ds2 = 0.0_f64   
     ds4 = 0.0_f64   
     sf2 = cf(3,1)*(x(2)-x(1))/a

   else if (ind1 > 0) then

     dp1 = 1.0_f64
     dp2 = 0.0_f64
     dp3 = 0.0_f64
     dp4 = 1.0_f64
     ds1 = 0.0_f64
     ds2 = 0.0_f64
     ds3 = 0.0_f64
     ds4 = 0.0_f64
     sf1 = -cf(2,1)
     sf2 =  cf(3,1)*(x(2)-x(1))/a

   end if

   goto 55

50 continue

   dp1=dp1-di1*h(im-3)-di3*h(im-2)
   dp2=dp2-di2*h(im-3)-di4*h(im-2)
   dp3=dp3-di1*h(im-1)-di3*h(im)
   dp4=dp4-di2*h(im-1)-di4*h(im)

55 continue

   det=dp1*dp4-dp2*dp3

   if (i/=1) then

     sf1=sf1+di1*h(im-5)+di3*h(im-4)
     sf2=sf2+di2*h(im-5)+di4*h(im-4)

   end if

   h(im+1)=( dp4*sf1-dp3*sf2)/det
   h(im+2)=(-dp2*sf1+dp1*sf2)/det

   if(i>=n) goto 65


   h(im+3)=(+dp4*ds1-dp3*ds2)/det
   h(im+4)=(-dp2*ds1+dp1*ds2)/det
   h(im+5)=(+dp4*ds3-dp3*ds4)/det
   h(im+6)=(-dp2*ds3+dp1*ds4)/det
   im=im+6
   i =i +1
   
   if (i ==2) goto 10

   goto 5

   !Backward substitution and solution to the algebraic system

65 continue

   cf(2,n) = -h(im+1)
   dpd1    =  h(im+2)*a
   cf(3,n) = dpd1/(x(n)-x(n-1))
   im      = im-6

   do i=2,n
     k=n+1-i
     cf(2,k) = -h(im+1)+h(im+3)*cf(2,k+1)-h(im+5)*dpd1/a
     dpd2    =  h(im+2)*a-h(im+4)*cf(2,k+1)*a+h(im+6)*dpd1
     cf(3,k) =  dpd2/(x(k+1)-x(k))
     dpd1    =  dpd2
     im      =  im-6
   end do

end subroutine inspl5

!> Calculation of the values of an interpolating quintic pline
!> and of its first and second derivatives at any point xx
subroutine splin5(n,x,cf,xx,order,f)
sll_int32,  intent(in)  :: n          !< number of interpolation points
sll_real64, intent(in)  :: x(n)       !< vector of abscissae
sll_real64, intent(in)  :: cf(1:3,n)  !< ordinates, first and second derivatives
sll_real64, intent(in)  :: xx         !< abscissae where values of the function
                                         !< and its derivatives are computed

sll_int32,  intent(in)  :: order      !< order of derivative
sll_real64, intent(out) :: f          !< 0: value of the interpolating function
                                      !< 1: value of its first derivative
                                      !< 2: value of its second derivative

sll_real64 :: cc, h, u, w, y
sll_real64 :: xn, xr1, xr2, xr3
sll_int32  :: i


f = 0.0_f64

xn = (xx-x(1))/(x(n)-x(1))*(n-1)

i = ceiling(xn)+1

if ( xx - x(i) == 0.0_f64 ) then

  f = cf(order+1,i)

else

  cc   = cf(1,i-1)-cf(1,i)
  h    = x(i)-x(i-1)
  xn   = xx-x(i-1)
  xr1  = xn/h
  xr2  = xr1*xr1
  xr3  = xr2*xr1
  
  select case(order)

  case(0)

    y = cf(1,i-1)+cc*xr3*(-10.0_f64+xr1*(+15.0_f64-6.0_f64*xr1))
    
    w = cf(2,i-1)*xr1 &
       +cf(2,i-1)*xr3*(-6.0_f64+xr1*(+8.0_f64-3.0_f64*xr1)) &
       -cf(2,i  )*xr3*(+4.0_f64+xr1*(-7.0_f64+3.0_f64*xr1))
    
    u = cf(3,i-1)*xr2 &
       +cf(3,i-1)*xr3*(-3.0_f64+xr1*(+3.0_f64-  xr1)) &
       +cf(3,i  )*xr3*(+1.0_f64+xr1*(-2.0_f64+  xr1)) 
    
    f = y+h*(w+h*u*0.5_f64)
  
  case(1) 
  
    y =      30.0_f64*cc*xr2*(-1.0_f64+xr1*(+2.0_f64-  xr1)) 
    
    w = cf(2,i-1)     &
       +cf(2,i-1)*xr2*(-18.0_f64+xr1*(+32.0_f64-15.0_f64*xr1)) &
       -cf(2,i  )*xr2*(+12.0_f64+xr1*(-28.0_f64+15.0_f64*xr1))
    
    u = cf(3,i-1)*xr1*2.0_f64 &
       +cf(3,i-1)*xr2*(- 9.0_f64+xr1*(+12.0_f64- 5.0_f64*xr1)) &
       +cf(3,i  )*xr2*(+ 3.0_f64+xr1*(- 8.0_f64+ 5.0_f64*xr1)) 
    
    f = y/h+w+h*u*0.5_f64
  
  case(2)
  
    xn= 10.0_f64*xr1
  
    y = cc*xn*(-1.0_f64+xr1*(+3.0_f64-2.0_f64*xr1))
    
    w = +cf(2,i-1)*(-6.0_f64+xr1*(+16.0_f64-xn)) &
        -cf(2,i  )*(+4.0_f64+xr1*(-14.0_f64+xn))
    
    u =  cf(3,i-1) &
        +cf(3,i-1)*(- 9.0_f64+xr1*(+18.0_f64- xn))*xr1 &
        +cf(3,i  )*(+ 3.0_f64+xr1*(-12.0_f64+ xn))*xr1 
    
    f = (y/h+xr1*w)*6.0_f64/h+u

  end select

end if

end subroutine splin5

subroutine inspl5_periodic(n,dx,cf,h)

sll_int32,  intent(in)    :: n          !< number of interpolation points
sll_real64, intent(in)    :: dx         !< vector of abscissae
sll_real64, intent(inout) :: cf(1:3,n)  !< ordinates, first and second derivatives
sll_real64, intent(out)   :: h(6*n-3)   !< auxiliary vector

sll_real64 :: a, as, det, di1, di2, di3, di4
sll_real64 :: dp1, dp2, dp3, dp4, dpd1, dpd2
sll_real64 :: ds1, ds2, ds3, ds4
sll_real64 :: fp, fp3, fp4
sll_real64 :: p1, p2, p3, q1, q2
sll_real64 :: sf1, sf2
sll_int32  :: i, i1, i2, im, k

h = 0.0_f64

as  = 64.0_f64/ 3.0_f64
a   = sqrt(as)
p1  = 0.0_f64
p2  = 0.0_f64
p3  = 0.0_f64
q1  = 0.0_f64
q2  = 0.0_f64
sf1 = 0.0_f64
sf2 = 0.0_f64
im  = 0

!Forward elimination procedure for Choleskys algorithm
!applied to block-tridiagonal systems of algebraic equations.

do i = 1, n

   if (i>1) then

     sf1 = fp4
     sf2 = fp3*q1*a

   end if

   if (i>2) then

     di1 = ds1
     di2 = ds3
     di3 = ds2
     di4 = ds4

   end if


   dp1 = 64.0_f64*p3
   dp2 = 12.0_f64*p2*q1*a
   dp3 = dp2
   dp4 = 3.0_f64*p1*q2*as

   if (i >= n) then !boundary condition at xn

     dp2 = 0.0_f64
     dp4 = 1.0_f64
     di2 = 0.0_f64
     di4 = 0.0_f64
     sf2 = cf(3,n)*dx/a

     goto 50

   end if

   p1 = 1.0_f64/dx
   p2 = p1*p1
   p3 = p1*p2
   q1 = p1
   q2 = p2

!   if (i < n-1) then
     q1 = 1.0_f64 / dx
     q2 = q1*q1
!   end if

   i1 = modulo(i-1,n-1)+1
   i2 = modulo(i,n-1)+1
   fp = cf(1,i1)-cf(1,i2)
   fp3= 20.0_f64*p3*fp
   fp4=  6.0_f64*p1*fp3
   dp1= dp1+64.0_f64*p3
   dp2= dp2-12.0_f64*p3*a
   dp3= dp2
   dp4= dp4+3.0_f64*p3*as
   ds1=-56.0_f64*p3
   ds2= 8.0_f64*p3*a
   ds3= -8.0_f64*p2*q1*a
   ds4=    p2*q1*as
   sf1=sf1+fp4
   sf2=sf2-fp3*p1*a

   if (i>1) goto 50

   di1 = ds1
   di2 = ds3
   di3 = ds2
   di4 = ds4

   dp2 = 0.0_f64
   dp4 = 1.0_f64
   ds2 = 0.0_f64   
   ds4 = 0.0_f64   
   sf2 = cf(3,1)*dx/a

   goto 55

50 continue

   dp1=dp1-di1*h(im-3)-di3*h(im-2)
   dp2=dp2-di2*h(im-3)-di4*h(im-2)
   dp3=dp3-di1*h(im-1)-di3*h(im)
   dp4=dp4-di2*h(im-1)-di4*h(im)

55 continue

   det=dp1*dp4-dp2*dp3

   if (i>1) then

     sf1=sf1+di1*h(im-5)+di3*h(im-4)
     sf2=sf2+di2*h(im-5)+di4*h(im-4)

   end if

   h(im+1)=( dp4*sf1-dp3*sf2)/det
   h(im+2)=(-dp2*sf1+dp1*sf2)/det

   if (i<n) then

     h(im+3)=(+dp4*ds1-dp3*ds2)/det
     h(im+4)=(-dp2*ds1+dp1*ds2)/det
     h(im+5)=(+dp4*ds3-dp3*ds4)/det
     h(im+6)=(-dp2*ds3+dp1*ds4)/det
     im=im+6

   end if

 end do

   !Backward substitution and solution to the algebraic system

   cf(2,n) = -h(im+1)
   dpd1    =  h(im+2)*a
   cf(3,n) = dpd1/dx
   im      = im-6

   do i=2,n
     k=n+1-i
     cf(2,k) = -h(im+1)+h(im+3)*cf(2,k+1)-h(im+5)*dpd1/a
     dpd2    =  h(im+2)*a-h(im+4)*cf(2,k+1)*a+h(im+6)*dpd1
     cf(3,k) =  dpd2/dx
     dpd1    =  dpd2
     im      =  im-6
   end do

end subroutine inspl5_periodic

end module sll_m_quintic_splines
