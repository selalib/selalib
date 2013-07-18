!> @brief 
!> Selalib poisson solvers (1D, 2D and 3D) unit test
!> Start date: March 20, 2012
!> Last modification (decoupling the tests): October 26, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
program test_qns2d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_qn_solver_2d
use sll_boundary_condition_descriptors

implicit none

sll_int32, parameter :: NP_r = 256
sll_int32, parameter :: NP_theta = 256

sll_int32   :: BC 
sll_real64  :: dr
sll_real64  :: dtheta
sll_real64  :: rmin, rmax, Zi

sll_real64, dimension(np_r)     :: Te
sll_real64, dimension(np_r)     :: r
sll_real64, dimension(np_r)     :: c
sll_real64, dimension(np_theta) :: theta
sll_real64, dimension(np_theta) :: f
sll_real64, dimension(np_theta) :: g

sll_real64, dimension(NP_r,NP_theta)  :: rho, phi
sll_real64, dimension(NP_r,NP_theta)  :: phi_exact

sll_int32                             :: i, j
sll_real64                            :: err
sll_real64                            :: err_bound
sll_real64                            :: Mr, Mtheta

rmin = 1.d0
rmax = 10.d0
Zi = 1.d0

dtheta = 2*sll_pi/NP_theta
do j = 1, np_theta
   theta(j) = (j-1)*dtheta
end do

Te = 1.d0

do bc=1,2

   print*, ' '
   if ( BC == SLL_NEUMANN) then
      dr = (rmax-rmin)/(NP_r-1)
      do i = 1, np_r
         r(i) = rmin + (i-1)*dr
         c(i) = 2/r(i)
      end do
      print*, 'Testing with neumann'
   else
      dr = (rmax-rmin)/(NP_r+1)
      do i = 1, np_r
         r(i) = rmin + i*dr
         c(i) = (rmax+rmin-2*r(i)) / ( (rmax-r(i))*(r(i)-rmin) )
      end do
      print*, 'Testing with dirichlet'
   endif
   print*, ' '

   call test_function_1()
   call test_function_2()

enddo

print*, 'PASSED'

contains

subroutine test_function_1()

type (qn_solver_2d), pointer :: plan

plan => new(BC, rmin, rmax, NP_r, NP_theta)

if (BC==SLL_NEUMANN) then
   f = sin(rmax-rmin)*cos(theta)
else
   f = 0.d0
endif

err_bound  = 0.d0

do j=1,NP_theta

   Mr = 4*abs(cos(theta(j)))

   do i=1,NP_r
         
      phi_exact(i,j) = sin(r(i)-rmin)*sin(rmax-r(i))*cos(theta(j))

      rho(i,j) = cos(theta(j))           &
               * (2*cos(rmin+rmax-2*r(i))   &
               - c(i)*sin(rmin+rmax-2*r(i)) &
               + (1/r(i)**2+1/(Zi*Te(i)))*sin(rmax-r(i))*sin(r(i)-rmin))

      Mtheta = abs(sin(r(i)-rmin)*sin(rmax-r(i)))
      err_bound = err_bound + &
             Mr*dr**2/12 + abs(c(i))*Mr*dr**2/6 + Mtheta*dtheta**2/(r(i)**2*12)

    enddo

enddo

g = -f

call solve(plan, rho, c, Te, f, g, Zi, phi)

err = sum(abs(phi_exact-phi))/(NP_r*NP_theta)
err_bound = err_bound/(NP_r*NP_theta)

print*, 'Error =', err
print*, 'Boundary error =', err_bound
print*, ' '

if ( err > err_bound ) then
   print*, 'test_qns2d FAILED'
   stop
endif
    
call delete(plan)

end subroutine test_function_1

subroutine test_function_2()

type (qn_solver_2d), pointer :: plan

plan => new(BC, rmin, rmax, NP_r, NP_theta)

if (BC==SLL_NEUMANN) then
   f = sin(rmax-rmin) * exp(-.5*(theta-sll_pi)**2)/sqrt(2*sll_pi)
else
   f = 0.d0
endif

err_bound  = 0.d0

do j=1,NP_theta

   Mr = 4*abs(cos(theta(j)))

   do i=1,NP_r

         
      phi_exact(i,j)  = sin(r(i)-rmin)*sin(rmax-r(i)) * &
                  exp(-.5*(theta(j)-sll_pi)**2)/sqrt(2*sll_pi)

      rho(i,j) = ( 2*cos(rmax+rmin-2*r(i)) - c(i)*sin(rmax+rmin-2*r(i)) ) &
               * exp(-.5*(theta(j)-sll_pi)**2)/sqrt(2*sll_pi)       &
               + phi_exact(i,j) &
               * ( 1/(Zi*Te(i)) - ((theta(j)-sll_pi)**2-1)/r(i)**2 )

      Mtheta = abs(sin(r(i)-rmin)*sin(rmax-r(i)))
      err_bound = err_bound + &
      Mr*dr**2/12 + abs(c(i))*Mr*dr**2/6 + Mtheta*dtheta**2/(r(i)**2*12)

   enddo

enddo

g = -f

call solve(plan, rho, c, Te, f, g, Zi, phi)

err = sum(abs(phi_exact-phi))/(NP_r*NP_theta)
err_bound = err_bound/(NP_r*NP_theta)

print*, 'Error =', err
print*, 'Boundary error =', err_bound
print*, ' '

if ( err > err_bound ) then
   print*, 'test_qns2d FAILED'
   stop
endif
    
call delete(plan)

end subroutine test_function_2


end program test_qns2d
