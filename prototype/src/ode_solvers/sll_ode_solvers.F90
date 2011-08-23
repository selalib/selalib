module ode_solvers
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  
  implicit none
  enum, bind(C)
     enumerator :: PERIODIC_ODE = 0, COMPACT_ODE = 1
  end enum

contains

  ! Computes the solution of functional equation 
  ! 
  !            xi-xout = c*deltat*(b(xi)+a(xout))
  !
  ! for initial conditions xi corresponding to all points of a uniform grid
  ! with c=1 and b=0 this is a first order method for solving backward 
  !
  !           dx/dt = a(x,t)
  !
  ! with c=1/2 and b = a(.,t_(n+1)) this is a second order method for the 
  ! same problem.
  !
  ! In practise the first order method needs to be called in order to 
  ! compute b for the second order method
  subroutine implicit_ode( order,  &
                           deltat, &
                           xmin,   &
                           ncx,    &
                           deltax, &
                           bt,     &
                           xout,   &
                           a,      &
                           a_np1 ) 
    intrinsic  :: floor, present
    sll_int32  :: order
    sll_real64 :: deltat   
    sll_real64 :: xmin  
    sll_int32  :: ncx   ! number of cells of uniform grid
    sll_real64 :: deltax
    sll_int32  :: bt    ! boundary_type
    ! solution for all initial conditions:
    sll_real64, dimension(:)                     :: xout   
    sll_real64, dimension(:)                     :: a     ! rhs at t = t_n
    sll_real64, dimension(:), pointer, optional  :: a_np1 ! rhs at t = t_n+1
    ! local variables
    sll_int32  :: i, id, ileft, iright
    sll_real64 :: alpha, alphabar, xmax, xi
    sll_real64 :: c     ! real coefficient
    sll_real64, dimension(ncx+1), target     :: zeros   ! array if zeros
    sll_real64, dimension(:), pointer        :: b

    ! initialize zeros
    zeros = 0.0_f64
    ! check order. The implementation with a 'select' construct permits
    ! to extend this solver to higher orders more conveniently.
    select case (order)
    case (1)
       c = 1.0_f64
       b => zeros
    case (2)
       c = 0.5_f64
       if (present(a_np1)) then
          b => a_np1
       else
          stop 'implicit_ode: need field at time t_n+1 for higher order'
       end if
    case default
       print*, 'order = ',order, ' not implemented'
       stop
    end select

    ! compute xmax of the grid
    SLL_ASSERT(size(a)==ncx+1)
    SLL_ASSERT(size(b)==ncx+1)
    SLL_ASSERT(size(xout)==ncx+1)

    xmax = xmin + ncx * deltax
    do i = 1, ncx + 1
       xi = xmin + (i-1)*deltax  ! current grid point
       ! estimate the displacement alpha * deltax = xi - xout
       ! which will give the cell around the center of which a will be 
       ! Taylor expanded
       ! print*, 's1', i, a(i), deltax, deltat, bt
       alphabar = deltat / deltax * a(i)
       id = floor( -alphabar )  ! cell will be [i+id, i+id+1]
       !print*, 's2', alphabar, id, a(i)
       ! handle boundary conditions
       if (bt == PERIODIC_ODE) then
          ileft = modulo(i+id-1,ncx)+1
          iright = modulo(i+id,ncx)+1
       else if (bt == COMPACT_ODE) then
          ileft = min(max(i+id,1),ncx+1)
          iright = max(min(i+id+1,ncx+1),1)
       else
          stop 'compute_flow_1D_backward : boundary_type not implemented' 
       end if
       !print*, i, ileft, iright, a(iright) - a(ileft),  deltax + c * deltat * (a(iright) - a(ileft))
       SLL_ASSERT((ileft>=1).and.(ileft<= ncx+1))
       SLL_ASSERT((iright>=1).and.(iright<= ncx+1))
       SLL_ASSERT( deltax + c * deltat * (a(iright) - a(ileft)) > 0.0 )
       ! compute xout using first order Taylor expansion of a to get 
       ! linear equation for alpha
       alpha = c*deltat*(b(i) + a(ileft)*(1+id) - id*a(iright)) &
            /( deltax + c * deltat * (a(iright) - a(ileft)))
       !print*,i,'alpha', alpha, id, ileft,iright,a(i),a(iright)-a(ileft), b(i)
       !print*,i,'b ',b(i)
       xout(i) = xi - alpha * deltax 
       ! handle boundary conditions
       if (bt == PERIODIC_ODE) then
          xout(i) = modulo(xout(i)-xmin,xmax-xmin) + xmin 
       else if (bt == COMPACT_ODE) then
          if (xout(i) < xmin ) then
             ! put particles on the left of the domain on the left boundary
             xout(i) = xmin   
          elseif (xout(i) > xmax ) then
             ! put particles on the right of the domain on the right boundary
             xout(i) = xmax   
          end if
       else
          stop 'compute_flow_1D_backward : boundary_type not implemented' 
       end if
       !print*, 'compute_flow_1D_backward' , xmin, xout(i),xmax
       SLL_ASSERT((xout(i) >= xmin ) .and. (xout(i) <= xmax)) 
    end do
  end subroutine implicit_ode


end module ode_solvers
