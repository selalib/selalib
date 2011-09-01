module ode_solvers
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
  
  implicit none
  enum, bind(C)
     enumerator :: PERIODIC_ODE = 0, COMPACT_ODE = 1
  end enum

  abstract interface
     function scalar_function_1D( eta )
       use sll_working_precision
       sll_real64 :: scalar_function_1D
       sll_real64, intent(in)  :: eta
     end function scalar_function_1D
  end interface

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
          stop 'implicit_ode : boundary_type not implemented' 
       end if
       SLL_ASSERT((xout(i) >= xmin ) .and. (xout(i) <= xmax)) 
    end do
  end subroutine implicit_ode


  ! Computes the solution of functional equation obtained on a curvilinear grid
  ! 
  !            xi-x(eta_out) = c*deltat*(b(xi)+a(eta_out))
  !
  ! for initial conditions xi=x(eta_i) corresponding to the images of all points 
  ! of a uniform gridwith c=1 and b=0 this is a first order method for solving backward 
  !
  !           dx(eta)/dt = a(x(eta(t)),t)
  !
  ! with c=1/2 and b = a(.,t_(n+1)) this is a second order method for the 
  ! same problem.
  !
  ! In practise the first order method needs to be called in order to 
  ! compute b for the second order method
  subroutine implicit_ode_curv( order,       &
                                deltat,      &
                                eta_min,     &
                                nc_eta,      &
                                delta_eta,   &
                                bt,          &
                                eta_out,     &
                                xfunc,       &
                                xprimefunc,  &
                                a,           &
                                a_np1)        
    intrinsic  :: floor, present
    sll_int32  :: order
    sll_real64 :: deltat   
    sll_real64 :: eta_min  
    sll_int32  :: nc_eta   ! number of cells of uniform grid
    sll_real64 :: delta_eta
    sll_int32  :: bt    ! boundary_type
    ! solution for all initial conditions:
    sll_real64, dimension(:)                     :: eta_out   
    procedure(scalar_function_1D), pointer       :: xfunc
    procedure(scalar_function_1D), pointer       :: xprimefunc
    sll_real64, dimension(:)                     :: a     ! rhs at t = t_n
    sll_real64, dimension(:), pointer, optional  :: a_np1 ! rhs at t = t_n+1
    ! local variables
    sll_int32  :: i, id
    sll_real64 :: alpha, aint, eta_max, eta_i, eta_k, eta_kp1, xi
    sll_real64 :: c     ! real coefficient
    sll_real64, dimension(nc_eta+1), target     :: zeros   ! array if zeros
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
          stop 'implicit_ode_curv: need field at time t_n+1 for higher order'
       end if
    case default
       print*, 'order = ',order, ' not implemented'
       stop
    end select

    ! compute eta_max of the grid
    SLL_ASSERT(size(a)==nc_eta+1)
    SLL_ASSERT(size(b)==nc_eta+1)
    SLL_ASSERT(size(eta_out)==nc_eta+1)

    eta_max = eta_min + nc_eta * delta_eta
    eta_i = eta_min
    do i = 1, nc_eta + 1
       xi = xfunc( eta_i ) ! current grid point
       ! initial guess for the Newton solver: take eta_i
       eta_kp1 = eta_i
       do while (abs(eta_k-eta_kp1) > 1.0d-8)
          eta_k = eta_kp1
          ! handle boundary conditions         
          if (bt == PERIODIC_ODE) then
             eta_k = eta_min + modulo(eta_k - eta_min, eta_max - eta_min)
          else if (bt == COMPACT_ODE) then
             if (eta_k < eta_min) then
                eta_k = eta_min
                cycle
             else if (eta_k > eta_max) then
                eta_k = eta_max
                cycle
             end if
          else
             stop 'implicit_ode_curv: boundary_type not implemented' 
          end if
          ! localize cell [id, id+1] where eta_k is in 
          id = floor( (eta_k-eta_min)/delta_eta )  
          alpha = (eta_k - eta_min)/delta_eta - id 

          ! compute linear interpolation of a at eta_k
          aint = (1.0_f64 - alpha) * a(id) + alpha * a(id+1)
          ! compute next iterate of Newton's method
          eta_kp1 = eta_k - (c*deltat*(b(i)*aint) + xfunc(eta_k) - xi) / &
               ( c*deltat*(a(id+1)-a(id))/delta_eta + xprimefunc(eta_k) )
       end do

       SLL_ASSERT((eta_kp1 >= eta_min ) .and. (eta_kp1 <= eta_max)) 
       eta_out(i) = eta_kp1
       eta_i = eta_i + delta_eta
    end do
  end subroutine implicit_ode_curv

  ! Classical second order Runge-Kutta ODE solver for an ode of the form
  ! d eta/ dt = a(eta)
  ! a is known only at grid points and linear interpolation is used in between
  subroutine rk2( nsubsteps,   &
                  deltat,      &
                  eta_min,     &
                  nc_eta,      &
                  delta_eta,   &
                  bt,          &
                  eta_out,     &
                  a,           &
                  mesh) 
    intrinsic  :: floor
    sll_int32  :: nsubsteps
    sll_real64 :: deltat   
    sll_real64 :: eta_min  
    sll_int32  :: nc_eta   ! number of cells of uniform grid
    sll_real64 :: delta_eta
    sll_int32  :: bt    ! boundary_type
    ! solution for all initial conditions:
    sll_real64, dimension(:)                     :: eta_out   
    sll_real64, dimension(:)                     :: a     ! rhs of ode
    type(mesh_descriptor_2D)                     :: mesh   
    ! local variables
    sll_real64 :: eta_max
    sll_real64 :: eta_i, eta_k, eta_kp1
    sll_real64 :: deltatsub
    sll_real64 :: a_n, a_np1, alpha
    sll_int32  :: i, id, isub

    SLL_ASSERT(size(a)==nc_eta+1)
    SLL_ASSERT(size(eta_out)==nc_eta+1)
    ! compute eta_max of the grid
    eta_max = eta_min + nc_eta * delta_eta
    do i = 1, nc_eta + 1
       eta_i = eta_min + (i-1)*delta_eta  ! current grid point
       ! loop over substeps
       eta_k = eta_i
       a_n = a(i) / mesh%geom%Jacobian(eta_i,1.0_f64)
       deltatsub = deltat / nsubsteps
       do isub = 1, nsubsteps
          ! first stage
          eta_kp1 = eta_k + deltatsub * a_n
          ! handle boundary conditions         
          if (bt == PERIODIC_ODE) then
             eta_kp1 = eta_min + modulo(eta_kp1 - eta_min, eta_max - eta_min)
          else if (bt == COMPACT_ODE) then
             if (eta_kp1 < eta_min) then
                eta_kp1 = eta_min
             else if (eta_kp1 > eta_max) then
                eta_kp1 = eta_max
             end if
          else
             stop 'sll_ode_solvers, rk2: boundary_type not implemented' 
          end if

          ! localize cell [id, id+1] where eta_k is in 
          id = floor( (eta_kp1-eta_min)/delta_eta ) 
          alpha = (eta_kp1 - eta_min)/delta_eta - id 
          ! compute linear interpolation of a at eta_k
          a_np1 = (1.0_f64 - alpha) * a(id+1) + alpha * a(id+2)
          ! divide by jacobian of mesh
          a_np1 = a_np1 / mesh%geom%Jacobian(eta_kp1,1.0_f64)
          ! compute cubic Lagrange interpolation of a at eta_k
          !a_np1 = -alpha*(alpha-1)*(alpha-2)/6 * a(id) + (alpha+1)*(alpha-1)*(alpha-2)/2 * a(id+1) &
          !     - (alpha+1)*alpha*(alpha-2)/2 * a(id+2) + (alpha+1)*alpha*(alpha-1)/6 * a(id+3) 
          ! compute solution of ode for current grid point 
          eta_kp1 = eta_k + 0.5_f64 * deltatsub * (a_n + a_np1)
          ! handle boundary conditions      
          !print*, i, eta_kp1, eta_min, eta_max
          if (bt == PERIODIC_ODE) then
             eta_kp1 = eta_min + modulo(eta_kp1  - eta_min, eta_max - eta_min)
          else if (bt == COMPACT_ODE) then
             if (eta_kp1 < eta_min) then
                eta_kp1 = eta_min
             else if (eta_kp1 > eta_max) then
                eta_kp1 = eta_max
             end if
          else
             stop 'sll_ode_solvers, rk2: boundary_type not implemented' 
          end if
          !print*, i, eta_kp1
          eta_k = eta_kp1
          a_n = a_np1
       end do
       eta_out(i) = eta_k

       SLL_ASSERT((eta_out(i) >= eta_min ) .and. (eta_out(i) <= eta_max)) 
    end do
  end subroutine rk2
end module ode_solvers
