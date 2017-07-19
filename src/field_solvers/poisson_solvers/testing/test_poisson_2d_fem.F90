program test_poisson_2d_fem
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_poisson_solvers_macros.h"

use sll_m_constants
use sll_m_poisson_2d_fem

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_int32                           :: i
sll_int32                           :: j
sll_real64                          :: xmin
sll_real64                          :: xmax
sll_real64                          :: ymin
sll_real64                          :: ymax
sll_real64                          :: dx
sll_real64                          :: dy
sll_int32                           :: nc_x
sll_int32                           :: nc_y
sll_int32                           :: error
sll_real64                          :: dpi
sll_real64, dimension(:),   pointer :: x
sll_real64, dimension(:),   pointer :: y
sll_real64, dimension(:,:), pointer :: e_x
sll_real64, dimension(:,:), pointer :: e_y
sll_real64, dimension(:,:), pointer :: rho
sll_real64, dimension(:,:), pointer :: phi

nc_x = 64
nc_y = 64

SLL_CLEAR_ALLOCATE(e_x(1:nc_x+1,1:nc_y+1),error)  
SLL_CLEAR_ALLOCATE(e_y(1:nc_x+1,1:nc_y+1),error) 
SLL_CLEAR_ALLOCATE(rho(1:nc_x+1,1:nc_y+1),error)  
SLL_CLEAR_ALLOCATE(phi(1:nc_x+1,1:nc_y+1),error)  

SLL_ALLOCATE(x(1:nc_x+1),error)  
SLL_ALLOCATE(y(1:nc_y+1),error) 

xmin = -1.0_f64
xmax =  1.0_f64
ymin = -1.0_f64
ymax =  1.0_f64

dx = (xmax-xmin) / nc_x
dy = (ymax-ymin) / nc_y

!Create an irregular mesh
do i=1,nc_x+1
  x(i) = xmin+(i-1)*dx !* 0.5 * ((i-1)*dx+1)
enddo
do j=1,nc_y+1
  y(j) = ymin+(j-1)*dy !* 0.5 * ((j-1)*dy+1)
enddo

dpi = 2*sll_p_pi
do j = 1, nc_y+1
  do i = 1, nc_x+1
   phi(i,j) = x(i)*x(i) + y(j)*y(j)
   rho(i,j) = phi(i,j)
  end do
end do

rho(2:nc_x,2:nc_y) = -4.0_f64

call test_compact()

print*, 'PASSED'

contains

subroutine test_compact()
type( sll_t_poisson_2d_fem ) :: poisson
sll_real64 :: errmax

call sll_o_create(poisson, x, y, nc_x+1, nc_y+1)
call sll_o_solve(poisson, e_x, e_y, rho)

errmax = sum(abs(rho-phi)) / real((nc_x+1)*(nc_y+1),f64)

if ( errmax > 0.01 ) then
  stop 'Compact BC : FAILED'
end if

end subroutine test_compact

function fbc( x, y)
sll_real64 :: fbc
sll_real64 :: x, y

fbc = x*x + y*y

end function fbc

end program test_poisson_2d_fem
