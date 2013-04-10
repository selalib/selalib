module ode_solvers
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!#include "sll_field_2d.h"
  
  implicit none
  
  integer, parameter :: PERIODIC_ODE = 0, COMPACT_ODE = 1

#ifndef STDF95
  abstract interface
     function scalar_function_1D( eta )
       use sll_working_precision
       sll_real64 :: scalar_function_1D
       sll_real64, intent(in)  :: eta
     end function scalar_function_1D
  end interface
#endif

contains

  function f_one( eta )
    use sll_working_precision
    sll_real64 :: f_one
    sll_real64, intent(in)  :: eta

    f_one = 1.0_f64
  end function f_one


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
  !
  ! based on algorithm described in Crouseilles, Mehrenberger, Sonnendrucker, JCP 2010
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
    sll_int32  :: i, id, ileft, iright, i0
    sll_real64 :: xmax, xi, alpha
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

    ! localize cell [i0, i0+1] containing origin of characteristic ending at xmin
    i = 1
    if ( a(i) + b(i) > 0 ) then
       ! search on the left
       if (bt == PERIODIC_ODE) then
          i0 = 0 
          do while ( i0 + c*deltat/deltax*( a(modulo(i0-1,ncx)+1) + b(i) ) >= i  ) 
             i0 = i0 - 1
          end do
       else if (bt == COMPACT_ODE) then
          i0 = 1
       else
          stop 'implicit_ode : boundary_type not implemented' 
       end if
    else 
       ! search on the right
       i0 = 1 
       do while ( i0 + c*deltat/deltax*( a(i0) + b(i) ) < i  ) 
          i0 = i0 + 1
       end do
       i0 = i0 - 1
    end if
    
    do i = 1, ncx + 1
       xi = xmin + (i-1)*deltax  ! current grid point
       ! find cell which contains origin of characteristic
       do while ( i0 + c*deltat/deltax*( a(modulo(i0-1,ncx)+1) + b(i) ) <= i )
          i0 = i0 + 1
       end do
       i0 = i0 - 1
       !print*,  'out ',i, i0, i0 + c*deltat/deltax*( a(modulo(i0-1,ncx)+1) + b(i) )
       id = i - i0 
       ! handle boundary conditions
       if (bt == PERIODIC_ODE) then
          ileft = modulo(i0-1,ncx) + 1
          iright = modulo(i0,ncx) + 1
       else if (bt == COMPACT_ODE) then
          ileft = min(max(i0,1),ncx+1)
          iright = max(min(i0+1,ncx+1),1)
       else
          stop 'implicit_ode : boundary_type not implemented' 
       end if
       !print*, i, ileft, iright, a(iright) - a(ileft),  deltax + c * deltat * (a(iright) - a(ileft))
       SLL_ASSERT((ileft>=1).and.(ileft<= ncx+1))
       SLL_ASSERT((iright>=1).and.(iright<= ncx+1))
       SLL_ASSERT( deltax + c * deltat * (a(iright) - a(ileft)) > 0.0 )
       ! compute xout using linear interpolation of a 
       alpha = c*deltat*(b(i) + (1-id)*a(ileft) + id*a(iright)) &
            /( deltax + c * deltat * (a(iright) - a(ileft)))
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
       !print*,'interv ', xmin + (ileft-1)*deltax , xout(i), xmin + ileft*deltax
    end do
  end subroutine implicit_ode



  subroutine implicit_ode_nonuniform( order,  &
       deltat, &
       xin,   &
       ncx,    &
       bt,     &
       xout,   &
       a,      &
       a_np1 ) 
    intrinsic  :: floor, present
    sll_int32, intent(in)  :: order
    sll_real64, intent(in) :: deltat   
    sll_int32, intent(in)  :: ncx   ! number of cells of uniform grid
    sll_int32, intent(in)  :: bt    ! boundary_type
    sll_real64, dimension(:) :: xin(:)  
    ! solution for all initial conditions
    sll_real64, dimension(:)                     :: xout   
    sll_real64, dimension(:)                     :: a     ! rhs at t = t_n
    sll_real64, dimension(:), pointer, optional  :: a_np1 ! rhs at t = t_n+1
    ! local variables
    sll_int32  :: i, ileft, iright, i0, imax
    sll_real64 :: xi, xi0, yi0, yi0p1, x1, beta
    sll_real64 :: c     ! real coefficient
    sll_real64 :: period ! period of periodic domain
    sll_real64, dimension(ncx+1), target     :: zeros   ! array of zeros
    sll_real64, dimension(:), pointer        :: b
    sll_real64, parameter                    :: eps = 1.0e-14  ! small real number
    sll_real64 :: tmp
    ! initialize zeros
    zeros(:) = 0.0_f64
    ! check order. The implementation with a 'select' construct permits
    ! to extend this solver to higher orders more conveniently.

    yi0 = 0._f64
    select case (order)
    case (1)
       c = 1.0_f64
       b => zeros(1:ncx+1)
    case (2)
       c = 0.5_f64
       if (present(a_np1)) then
          b => a_np1(1:ncx+1)
       else
          stop 'implicit_ode: need field at time t_n+1 for higher order'
       end if
    case default
       print*, 'order = ',order, ' not implemented'
       stop
    end select

    ! check array sizes
    SLL_ASSERT(size(a)==ncx+1)
    SLL_ASSERT(size(b)==ncx+1)
    SLL_ASSERT(size(xout)==ncx+1)
    SLL_ASSERT(size(xin)==ncx+1)

    ! case of periodic boundary conditions
    !-------------------------------------
    if (bt == PERIODIC_ODE) then
       !begin modif 
       tmp = 0._f64  
       do i=1,ncx
          tmp = tmp+abs(a(i))
       enddo
       tmp = tmp/real(ncx,f64)
       if(tmp<1.e-14)then
          xout = xin
          return
       endif
       !end modif
       x1 = xin(1)
       period = xin(ncx+1)-x1
       
       do i = 1, ncx
          xi = xin(i)  ! current grid point
          ! check that displacement is less than 1 period for first point
          ! localize cell [i0, i0+1] containing origin of characteristic 
          ! ending at xin(1)
          ! we consider the mesh of the same period consisting of the points 
          ! yi0 = xi0 + c*deltat*( a(i0) + b(i) )
          ! modulo n. If there was no periodicity the sequence would be 
          ! strictly increasing
          ! we first look for i0 such that y(i0+1) < y(i0) due to periodicity
          if ( a(1) + b(i) > 0 ) then
             ! y(ncx+1) > x(1) in this case so we look backward 
             i0 = ncx + 1
             yi0p1 = modulo(xin(i0) + c*deltat*( a(i0) + b(i) )-x1, period)+x1
             i0 = ncx
             yi0 = modulo(xin(i0) + c*deltat*( a(i0) + b(i) ) - x1, period) + x1
             do while ( yi0p1 > yi0 ) 
                i0 = i0 - 1
                yi0p1 = yi0
                yi0 = modulo(xin(i0) + c*deltat*( a(i0) + b(i) )-x1, period)+x1
             end do
          else 
             ! search on the right
             i0 = 1
             yi0 = modulo(xin(i0) + c*deltat*( a(i0) + b(i) ) - x1, period) + x1
             yi0p1 = modulo(xin(i0+1) + c*deltat*( a(i0+1)+b(i) )-x1, period)+x1
             do while (yi0p1 > yi0) 
                i0 = i0 + 1
                yi0 = yi0p1
                yi0p1 = modulo(xin(i0+1)+c*deltat*(a(i0+1)+b(i))-x1, period)+x1
             end do
          end if
          
          if ((i0 < 1 ) .or. (i0 > ncx)) then ! yi0 is strictly increasing
             i0 = 1
             yi0 = modulo(xin(i0) + c*deltat*( a(i0) + b(i) ) - x1, period) + x1
             yi0p1 = modulo(xin(i0+1) + c*deltat*(a(i0+1)+b(i))-x1, period)+x1
          end if
          imax = i0
          ! find cell which contains origin of characteristic
          do while ( (yi0p1 < xi + eps))
             i0 = modulo(i0,ncx) + 1
             if (i0 == imax) then
                yi0 = yi0p1
                yi0p1 = modulo(xin(i0+1)+c*deltat*(a(i0+1)+b(i))-x1, period)+x1
                exit
             else
                yi0 = yi0p1
                yi0p1 = modulo(xin(i0+1)+c*deltat*(a(i0+1)+b(i))-x1, period)+x1
             end if
          end do

          SLL_ASSERT((i0>=1).and.(i0<= ncx))
          ! compute xout using linear interpolation of a 
          if (yi0p1 > yi0) then 
             beta = (xi - yi0)/(yi0p1 - yi0)
          else if (xi >= yi0) then
             beta = (xi - yi0)/(yi0p1 - yi0 + period)
          else 
             beta = (xi - yi0 + period)/(yi0p1 - yi0 + period)
          end if
          !print*, i, i0, yi0, xi, yi0p1, beta, period, x1
          if(.not.((beta>=-eps) .and. (beta <= 1)))then
            print *,xin
            print *,a
            print *,beta,eps,i0,beta-1._f64,i
            print *,'problem'
          endif
          !SLL_ASSERT((beta>=-eps) .and. (beta < 1))
          SLL_ASSERT((beta>=-eps) .and. (beta <= 1))
          xout(i) = xin(i0) + beta * (xin(i0+1)-xin(i0))
          ! handle periodic boundary conditions
          xout(i) = modulo(xout(i)-xin(1),xin(ncx+1)-xin(1)) + xin(1) 
          SLL_ASSERT((xout(i) >= xin(1) ) .and. (xout(i) <= xin(ncx+1))) 
       end do
       ! due to periodicity, origin of last point is same as origin of first 
       ! point
       xout(ncx+1) = xout(1)
    else if (bt == COMPACT_ODE) then
       ! localize cell [i0, i0+1] containing origin of characteristic ending 
       ! at xmin
       i = 1
       if ( a(i) + b(i) > 0 ) then
          i0 = 1
       else 
          ! search on the right
          i0 = 2 
          do while ( xin(i0) + c*deltat*( a(i0) + b(i) ) < xin(i)  ) 
             i0 = i0 + 1
          end do
          i0 = i0 - 1
       end if

       do i = 1, ncx + 1
          xi = xin(i)  ! current grid point
          ! find cell which contains origin of characteristic
          xi0 = xin(i0)
          if ( yi0 <= xi ) then
             i0 = i0 + 1
             yi0 = xin(i0) + c*deltat*( a(i0) + b(i) )
             do while ( (yi0  <= xi) .and. (i0 < ncx+1) )
                i0 = i0 + 1
                yi0 = xin(i0) + c*deltat*( a(i0) + b(i) )
                !print*, 'i0 ', i, i0, yi0, modulo(xi - y1, period) + y1
             end do
             i0 = i0 - 1
          else 
             do while ( (yi0 > xi) .and. (i0 > 1))
                i0 = i0 - 1
                yi0 = xin(i0) + c*deltat*( a(i0) + b(i) )
                !print*, 'i0 ', i, i0, yi0, xin(i0), xi, modulo(xi - y1, period) + y1
             end do
          end if

          ileft = i0
          iright = i0 + 1
          SLL_ASSERT((ileft>=1).and.(ileft<= ncx))
          SLL_ASSERT((iright>=2).and.(iright<= ncx+1))
          SLL_ASSERT( xin(ileft+1)-xin(ileft) + c * deltat * (a(iright) - a(ileft)) > 0.0 )
          ! compute xout using linear interpolation of a 
          beta = (xin(i)-xin(ileft)-c*deltat*(b(i) + a(ileft))) &
               /( xin(ileft+1)-xin(ileft) + c * deltat * (a(iright) - a(ileft)))
          xout(i) = xin(ileft) + beta * (xin(ileft+1)-xin(ileft)) 

          ! Handle characteristics that have gone out of domain
          if (xout(i) < xin(1) ) then
             ! put particles on the left of the domain on the left boundary
             xout(i) = xin(1)   
          elseif (xout(i) > xin(ncx+1) ) then
             ! put particles on the right of the domain on the right boundary
             xout(i) = xin(ncx+1)   
          end if
          SLL_ASSERT((xout(i) >= xin(1) ) .and. (xout(i) <= xin(ncx+1))) 
       end do
    else
       stop 'implicit_ode : boundary_type not implemented' 
    end if

  end subroutine implicit_ode_nonuniform












  subroutine implicit_ode_nonuniform_old( order,  &
       deltat, &
       xin,   &
       ncx,    &
       bt,     &
       xout,   &
       a,      &
       a_np1 ) 
    intrinsic  :: floor, present
    sll_int32  :: order
    sll_real64 :: deltat   
    sll_int32  :: ncx   ! number of cells of uniform grid
    sll_int32  :: bt    ! boundary_type
    sll_real64, dimension(:) :: xin(:)  
    ! solution for all initial conditions
    sll_real64, dimension(:)                     :: xout   
    sll_real64, dimension(:)                     :: a     ! rhs at t = t_n
    sll_real64, dimension(:), pointer, optional  :: a_np1 ! rhs at t = t_n+1
    ! local variables
    sll_int32  :: i, ileft, iright, i0, imax
    sll_real64 :: xi, xi0, yi0, yi0p1, x1, beta
    sll_real64 :: c     ! real coefficient
    sll_real64 :: period ! period of periodic domain
    sll_real64, dimension(ncx+1), target     :: zeros   ! array of zeros
    sll_real64, dimension(:), pointer        :: b
    sll_real64, parameter                    :: eps = 1.0e-14  ! small real number
    sll_real64 :: tmp
    ! initialize zeros
    zeros(:) = 0.0_f64
    ! check order. The implementation with a 'select' construct permits
    ! to extend this solver to higher orders more conveniently.
    select case (order)
    case (1)
       c = 1.0_f64
       b => zeros(1:ncx+1)
    case (2)
       c = 0.5_f64
       if (present(a_np1)) then
          b => a_np1(1:ncx+1)
       else
          stop 'implicit_ode: need field at time t_n+1 for higher order'
       end if
    case default
       print*, 'order = ',order, ' not implemented'
       stop
    end select

    ! check array sizes
    SLL_ASSERT(size(a)==ncx+1)
    SLL_ASSERT(size(b)==ncx+1)
    SLL_ASSERT(size(xout)==ncx+1)
    SLL_ASSERT(size(xin)==ncx+1)

    ! case of periodic boundary conditions
    !-------------------------------------
    if (bt == PERIODIC_ODE) then
      !begin modif 
      tmp = 0._f64  
      do i=1,ncx
        tmp = tmp+abs(a(i))
      enddo
      tmp = tmp/real(ncx,f64)
      if(tmp<1.e-14)then
        xout = xin
        return
      endif
      !end modif
       x1 = xin(1)
       period = xin(ncx+1)-x1

       do i = 1, ncx
          xi = xin(i)  ! current grid point
          ! check that displacement is less than 1 period for first point
          ! localize cell [i0, i0+1] containing origin of characteristic ending at xin(1)
          ! we consider the mesh of the same period consisting of the points yi0 = xi0 + c*deltat*( a(i0) + b(i) )
          ! modulo n. If there was no periodicity the sequence would be strictly increasing
          ! we first look for i0 such that y(i0+1) < y(i0) due to periodicity
          if ( a(1) + b(i) > 0 ) then
             ! y(ncx+1) > x(1) in this case so we look backward 
             i0 = ncx + 1
             yi0p1 = modulo(xin(i0) + c*deltat*( a(i0) + b(i) ) - x1, period) + x1
             i0 = ncx
             yi0 = modulo(xin(i0) + c*deltat*( a(i0) + b(i) ) - x1, period) + x1
             do while ( yi0p1 > yi0 ) 
                i0 = i0 - 1
                yi0p1 = yi0
                yi0 = modulo(xin(i0) + c*deltat*( a(i0) + b(i) ) - x1, period) + x1
             end do
          else 
             ! search on the right
             i0 = 1
             yi0 = modulo(xin(i0) + c*deltat*( a(i0) + b(i) ) - x1, period) + x1
             yi0p1 = modulo(xin(i0+1) + c*deltat*( a(i0+1) + b(i) ) - x1, period) + x1
             do while (yi0p1 > yi0) 
                i0 = i0 + 1
                yi0 = yi0p1
                yi0p1 = modulo(xin(i0+1) + c*deltat*( a(i0+1) + b(i) ) - x1, period) + x1
             end do
          end if
          
          if ((i0 < 1 ) .or. (i0 > ncx)) then ! yi0 is strictly increasing
             i0 = 1
             yi0 = modulo(xin(i0) + c*deltat*( a(i0) + b(i) ) - x1, period) + x1
             yi0p1 = modulo(xin(i0+1) + c*deltat*( a(i0+1) + b(i) ) - x1, period) + x1
          end if
          imax = i0
          ! find cell which contains origin of characteristic
          do while ( (yi0p1 < xi + eps))
             i0 = modulo(i0,ncx) + 1
             if (i0 == imax) then
                yi0 = yi0p1
                yi0p1 = modulo(xin(i0+1) + c*deltat*( a(i0+1) + b(i) ) - x1, period) + x1
                exit
             else
                yi0 = yi0p1
                yi0p1 = modulo(xin(i0+1) + c*deltat*( a(i0+1) + b(i) ) - x1, period) + x1
             end if
          end do

          SLL_ASSERT((i0>=1).and.(i0<= ncx))
          ! compute xout using linear interpolation of a 
          if (yi0p1 > yi0) then 
             beta = (xi - yi0)/(yi0p1 - yi0)
          else if (xi >= yi0) then
             beta = (xi - yi0)/(yi0p1 - yi0 + period)
          else 
             beta = (xi - yi0 + period)/(yi0p1 - yi0 + period)
          end if
          !print*, i, i0, yi0, xi, yi0p1, beta, period, x1
          if(.not.((beta>=-eps) .and. (beta <= 1)))then
            print *,xin
            print *,a
            print *,beta,eps,i0,beta-1._f64,i
            print *,'problem'
          endif
          !SLL_ASSERT((beta>=-eps) .and. (beta < 1))
          SLL_ASSERT((beta>=-eps) .and. (beta <= 1))
          xout(i) = xin(i0) + beta * (xin(i0+1)-xin(i0))
          ! handle periodic boundary conditions
          xout(i) = modulo(xout(i)-xin(1),xin(ncx+1)-xin(1)) + xin(1) 
          SLL_ASSERT((xout(i) >= xin(1) ) .and. (xout(i) <= xin(ncx+1))) 
       end do
       ! due to periodicity, origin of last point is same as origin of first point
       xout(ncx+1) = xout(1)
    else if (bt == COMPACT_ODE) then
       ! localize cell [i0, i0+1] containing origin of characteristic ending at xmin
       i = 1
       if ( a(i) + b(i) > 0 ) then
          i0 = 1
       else 
          ! search on the right
          i0 = 2 
          do while ( xin(i0) + c*deltat*( a(i0) + b(i) ) < xin(i)  ) 
             i0 = i0 + 1
          end do
          i0 = i0 - 1
       end if
       do i = 1, ncx + 1
          xi = xin(i)  ! current grid point
          ! find cell which contains origin of characteristic
          xi0 = xin(i0)
          if ( yi0 <= xi ) then
             i0 = i0 + 1
             yi0 = xin(i0) + c*deltat*( a(i0) + b(i) )
             do while ( (yi0  <= xi) .and. (i0 < ncx+1) )
                i0 = i0 + 1
                yi0 = xin(i0) + c*deltat*( a(i0) + b(i) )
                !print*, 'i0 ', i, i0, yi0, modulo(xi - y1, period) + y1
             end do
             i0 = i0 - 1
          else 
             do while ( (yi0 > xi) .and. (i0 > 1))
                i0 = i0 - 1
                yi0 = xin(i0) + c*deltat*( a(i0) + b(i) )
                !print*, 'i0 ', i, i0, yi0, xin(i0), xi, modulo(xi - y1, period) + y1
             end do
          end if
          ileft = i0
          iright = i0 + 1
          SLL_ASSERT((ileft>=1).and.(ileft<= ncx))
          SLL_ASSERT((iright>=2).and.(iright<= ncx+1))
          SLL_ASSERT( xin(ileft+1)-xin(ileft) + c * deltat * (a(iright) - a(ileft)) > 0.0 )
          ! compute xout using linear interpolation of a 
          beta = (xin(i)-xin(ileft)-c*deltat*(b(i) + a(ileft))) &
               /( xin(ileft+1)-xin(ileft) + c * deltat * (a(iright) - a(ileft)))
          xout(i) = xin(ileft) + beta * (xin(ileft+1)-xin(ileft)) 

          ! Handle characteristics that have gone out of domain
          if (xout(i) < xin(1) ) then
             ! put particles on the left of the domain on the left boundary
             xout(i) = xin(1)   
          elseif (xout(i) > xin(ncx+1) ) then
             ! put particles on the right of the domain on the right boundary
             xout(i) = xin(ncx+1)   
          end if
          SLL_ASSERT((xout(i) >= xin(1) ) .and. (xout(i) <= xin(ncx+1))) 
       end do
    else
       stop 'implicit_ode : boundary_type not implemented' 
    end if

  end subroutine implicit_ode_nonuniform_old



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
#ifdef STDF95
    sll_real64 :: xfunc
    sll_real64 :: xprimefunc
#else 
    procedure(scalar_function_1D), pointer       :: xfunc
    procedure(scalar_function_1D), pointer       :: xprimefunc
#endif
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
                  jac) 
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
#ifdef STDF95
    sll_real64, optional :: jac
    !local variables
    sll_real64           :: jacobian
#else
    procedure(scalar_function_1D), pointer, optional  :: jac
    ! local variables
    procedure(scalar_function_1D), pointer  :: jacobian
#endif
    sll_real64 :: eta_max
    sll_real64 :: eta_i, eta_k, eta_kp1
    sll_real64 :: deltatsub
    sll_real64 :: a_n, a_np1, alpha
    sll_int32  :: i, id, isub

#ifndef STDF95
    if (present(jac)) then
       jacobian => jac 
    else 
       jacobian => f_one
    end if
#endif

    SLL_ASSERT(size(a)==nc_eta+1)
    SLL_ASSERT(size(eta_out)==nc_eta+1)
    ! compute eta_max of the grid
    eta_max = eta_min + nc_eta * delta_eta
    do i = 1, nc_eta + 1
       eta_i = eta_min + (i-1)*delta_eta  ! current grid point
       ! loop over substeps
       eta_k = eta_i
#ifdef STDF95
       if(present(jac)) then
          a_n = a(i)/ jac(eta_i)
       else
          a_n = a(i) / 1.0_f64
       endif
#else
       a_n = a(i) / jacobian(eta_i)
#endif
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
#ifdef STDF95
          if(present(jac)) then
             a_np1 = a_np1 / jac(eta_kp1)
          else
             a_np1 = a_np1 / 1.0_f64
          endif
#else
          a_np1 = a_np1 / jacobian(eta_kp1)
#endif
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
          ! compute a_np1 at eta_kp1 for next substeps
          ! localize cell [id, id+1] where eta_k is in 
          id = floor( (eta_kp1-eta_min)/delta_eta ) 
          alpha = (eta_kp1 - eta_min)/delta_eta - id 
          ! compute linear interpolation of a at eta_k
          a_np1 = (1.0_f64 - alpha) * a(id+1) + alpha * a(id+2)
          ! divide by jacobian of mesh
#ifdef STDF95
          if(present(jac)) then
             a_np1 = a_np1 / jac(eta_kp1)
          else
             a_np1 = a_np1 / 1.0_f64
          endif
#else
          a_np1 = a_np1 / jacobian(eta_kp1)
#endif
          ! update
          eta_k = eta_kp1
          a_n = a_np1
       end do
       eta_out(i) = eta_k

       SLL_ASSERT((eta_out(i) >= eta_min ) .and. (eta_out(i) <= eta_max)) 
    end do
  end subroutine rk2

  ! Classical fourth order Runge-Kutta ODE solver for an ode of the form
  ! d eta/ dt = a(eta)
  ! a is known only at grid points and linear interpolation is used in between
    subroutine rk4( nsubsteps, &
                  deltat,      &
                  eta_min,     &
                  nc_eta,      &
                  delta_eta,   &
                  bt,          &
                  eta_out,     &
                  a,           &
                  jac) 
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
#ifdef STDF95
    sll_real64, optional :: jac
    !local variables
    sll_real64           :: jacobian
#else
    procedure(scalar_function_1D), pointer, optional  :: jac
    ! local variables
    procedure(scalar_function_1D), pointer  :: jacobian
#endif
    sll_real64 :: eta_max
    sll_real64 :: eta_i, eta_k, eta_kp1
    sll_real64 :: deltatsub
    sll_real64 :: a_n, a_np1, alpha, k2, k3, k4
    sll_int32  :: i, id, isub

#ifndef STDF95
    if (present(jac)) then
       jacobian => jac 
    else 
       jacobian => f_one
    end if
#endif

    SLL_ASSERT(size(a)==nc_eta+1)
    SLL_ASSERT(size(eta_out)==nc_eta+1)
    ! compute eta_max of the grid
    eta_max = eta_min + nc_eta * delta_eta
    do i = 1, nc_eta + 1
       eta_i = eta_min + (i-1)*delta_eta  ! current grid point
       ! loop over substeps
       eta_k = eta_i
#ifdef STDF95
       if(present(jac)) then
          a_n = a(i)/ jac(eta_i)
       else
          a_n = a(i) / 1.0_f64
       endif
#else
       a_n = a(i) / jacobian(eta_i)
#endif
       deltatsub = deltat / nsubsteps
       do isub = 1, nsubsteps
          ! second stage
          eta_kp1 = eta_k + 0.5_f64*deltatsub * a_n
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
          k2 = (1.0_f64 - alpha) * a(id+1) + alpha * a(id+2)
          ! divide by jacobian of mesh
#ifdef STDF95
          if(present(jac))then
             k2 = k2 / jac(eta_kp1)
          else
             k2 = k2 / 1.0_f64
          end if
#else
          k2 = k2 / jacobian(eta_kp1)
#endif
 
         ! third stage
          eta_kp1 = eta_k + 0.5_f64*deltatsub * k2
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
          k3 = (1.0_f64 - alpha) * a(id+1) + alpha * a(id+2)
          ! divide by jacobian of mesh
#ifdef STDF95
          if(present(jac))then
             k3 = k3 / jac(eta_kp1)
          else
             k3 = k3 / 1.0_f64
          end if
#else
          k3 = k3 / jacobian(eta_kp1)
#endif

          ! fourth stage
          eta_kp1 = eta_k + deltatsub * k3
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
          k4 = (1.0_f64 - alpha) * a(id+1) + alpha * a(id+2)
          ! divide by jacobian of mesh
#ifdef STDF95
          if(present(jac))then
             k4 = k4 / jac(eta_kp1)
          else
             k4 = k4 / 1.0_f64
          end if
#else
          k4 = k4 / jacobian(eta_kp1)
#endif


          ! compute solution of ode for current grid point 
          eta_kp1 = eta_k + deltatsub/6.0_f64 * (a_n + 2.0_f64*(k2+k3) + k4)
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
          ! compute a_np1 at eta_kp1 for next substeps
          ! localize cell [id, id+1] where eta_k is in 
          id = floor( (eta_kp1-eta_min)/delta_eta ) 
          alpha = (eta_kp1 - eta_min)/delta_eta - id 
          ! compute linear interpolation of a at eta_k
          a_np1 = (1.0_f64 - alpha) * a(id+1) + alpha * a(id+2)
          ! divide by jacobian of mesh
#ifdef STDF95
          if(present(jac)) then
             a_np1 = a_np1 / jac(eta_kp1)
          else
             a_np1 = a_np1 / 1.0_f64
          endif
#else
          a_np1 = a_np1 / jacobian(eta_kp1)
#endif
          ! update
          eta_k = eta_kp1
          a_n = a_np1
       end do
       eta_out(i) = eta_k

       SLL_ASSERT((eta_out(i) >= eta_min ) .and. (eta_out(i) <= eta_max)) 
    end do
  end subroutine rk4

  ! Classical fourth order Runge-Kutta ODE solver for an ode of the form
  ! d eta/ dt = a(eta)
  ! a is known only at grid points and linear interpolation is used in between
    subroutine rk4_non_unif( nsubsteps, &
                  deltat,      &
                  eta_min,     &
                  nc_eta,      &
                  delta_eta,   &
                  bt,          &
                  eta_out,     &
                  a,           &
                  xin) 
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
    sll_real64, dimension(:)                     :: xin     ! rhs of ode

    sll_real64 :: eta_max
    sll_real64 :: eta_i, eta_k, eta_kp1
    sll_real64 :: deltatsub
    sll_real64 :: a_n, a_np1, alpha, k2, k3, k4
    sll_int32  :: i, id, isub


    SLL_ASSERT(size(a)==nc_eta+1)
    SLL_ASSERT(size(eta_out)==nc_eta+1)
    ! compute eta_max of the grid
    eta_max = eta_min + nc_eta * delta_eta
    do i = 1, nc_eta + 1
       eta_i = eta_min + (i-1)*delta_eta  ! current grid point
       ! loop over substeps
       eta_k = eta_i
       a_n = a(i) 
       deltatsub = deltat / nsubsteps
       do isub = 1, nsubsteps
          ! second stage
          eta_kp1 = eta_k + 0.5_f64*deltatsub * a_n
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
          k2 = (1.0_f64 - alpha) * a(id+1) + alpha * a(id+2)
          ! divide by jacobian of mesh
          
 
         ! third stage
          eta_kp1 = eta_k + 0.5_f64*deltatsub * k2
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
          k3 = (1.0_f64 - alpha) * a(id+1) + alpha * a(id+2)
          ! divide by jacobian of mesh
          

          ! fourth stage
          eta_kp1 = eta_k + deltatsub * k3
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
          k4 = (1.0_f64 - alpha) * a(id+1) + alpha * a(id+2)
          ! divide by jacobian of mesh
          



          ! compute solution of ode for current grid point 
          eta_kp1 = eta_k + deltatsub/6.0_f64 * (a_n + 2.0_f64*(k2+k3) + k4)
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
          ! compute a_np1 at eta_kp1 for next substeps
          ! localize cell [id, id+1] where eta_k is in 
          id = floor( (eta_kp1-eta_min)/delta_eta ) 
          alpha = (eta_kp1 - eta_min)/delta_eta - id 
          ! compute linear interpolation of a at eta_k
          a_np1 = (1.0_f64 - alpha) * a(id+1) + alpha * a(id+2)
          ! divide by jacobian of mesh
          
          ! update
          eta_k = eta_kp1
          a_n = a_np1
       end do
       eta_out(i) = eta_k

       SLL_ASSERT((eta_out(i) >= eta_min ) .and. (eta_out(i) <= eta_max)) 
    end do
  end subroutine rk4_non_unif


end module ode_solvers
