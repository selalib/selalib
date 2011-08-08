module advection_field
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"

  use numeric_constants
  use sll_poisson_1d_periodic
  use sll_splines
  implicit none
contains
 ! computes the solution of functional equation xi-xout = c*deltat*(b(xi)+a(xout))
  ! for initial conditions xi corresponding to all points of a uniform grid
  ! with c=1 and b=0 this is a first order method for solving backward dx/dt = a(x,t)
  ! with c=1/2 and b = a(.,t_(n+1)) this is a second order method for the same problem
  ! In practise the first order method needs to be called in order to compute b for the second order method
  subroutine compute_flow_1D_backward( a, b, c, deltat, xmin, ncx, deltax,  bt, xout ) 
    intrinsic  :: floor
    sll_real64, dimension(:)   :: a, b    ! rhs functional equation
    sll_real64 :: xmin, deltax  ! first point and cell size of uniform mesh
    sll_real64 :: c     ! real coefficient
    sll_real64 :: deltat   ! time step
    sll_int32 :: ncx   ! number of cells of uniform grid
    sll_int32  :: bt    ! boundary_type
    sll_real64, dimension(:)   :: xout   ! solution for all initial conditions
    ! local variables
    sll_int32  :: i, id, ileft, iright
    sll_real64 :: alpha, alphabar, xmax, xi

    ! compute xmax of the grid
    xmax = xmin + ncx * deltax
    do i = 1, ncx + 1
       xi = xmin + (i-1)*deltax  ! current grid point
       ! estimate the displacement alpha * deltax = xi - xout
       ! which will give the cell around the center of which a will be Taylor expanded
       alphabar = deltat / deltax * a(i)
       id = floor( -alphabar )  ! cell will be [i+id, i+id+1]
       ! handle boundary conditions
       if (bt == PERIODIC) then
          ileft = modulo(i+id-1,ncx)+1
          iright = modulo(i+id,ncx)+1
       else if (bt == COMPACT) then
          ileft = max(i+id,1)
          iright = min (i+id+1,ncx+1)
       else
          stop 'compute_flow_1D_backward : boundary_type not implemented' 
       end if
       ! compute xout using first order Taylor expansion of a to get linear equation for alpha
       alpha = c * deltat * (b(i) + a(ileft)*(1+id) - id * a(iright))/( deltax + c * deltat * (a(iright) - a(ileft)))
       !print*,i,'alpha', alpha, id, ileft, iright, a(i), (i-1)*deltat / (1+deltat)
       xout(i) = xi - alpha * deltax 
       ! handle boundary conditions
       if (bt == PERIODIC) then
          xout(i) = modulo(xout(i),xmax-xmin) 
       else if (bt == COMPACT) then
          if (xout(i) < xmin ) then
             xout(i) = xmin   ! put particles on the left of the domain on the left boundary
          elseif (xout(i) > xmax ) then
             xout(i) = xmax   ! put particles on the right of the domain on the right boundary
          end if
       else
          stop 'compute_flow_1D_backward : boundary_type not implemented' 
       end if

       SLL_ASSERT((xout(i) >= xmin ) .and. (xout(i) <= xmax)) 
    end do

  end subroutine compute_flow_1D_backward


end module advection_field
