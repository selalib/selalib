!------------------------------------------------------------------------------
! Selalib
!------------------------------------------------------------------------------
!
! MODULE: gauss_legendre_integration
!
! DESCRIPTION:
!> @author Module Author Name and Affiliation
!> @brief This is a low-level mathematical utility that applies the 
!> Gauss-Legendre method to compute numeric integrals.
!> @details This module aims at providing a single interface to the process of 
!> integrating a function on a given interval.
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module gauss_legendre_integration
#include "sll_working_precision.h"
#include "sll_assert.h"
!  use sll_splines
  implicit none
  
  ! The following interface is supposed to represent any function of one
  ! real argument that returns a real value (double precision).
  !
  ! We are exploring the approach of having multiple versions of the
  ! integrators, each with an interface as simple as possible. Ultimately,
  ! all of these functions can be made to exist behind a single
  ! generic interface.
  ! 
  ! The ugly aspect of this approach is that it loses modularity: for 
  ! example, to specialize the integrator on the 'interpolated_function_1D'
  ! requires use of the spline module. 

  abstract interface
     function function_1D_legendre(x)
       use sll_working_precision ! can't pass a header file because the
                                 ! preprocessor prevents double inclusion.
                                 ! This is very rare.
       sll_real64             :: function_1D_legendre
       sll_real64, intent(in) :: x
     end function function_1D_legendre
  end interface

!  abstract interface
!     function interpolated_function_1D(x,spline_obj)
!       use sll_working_precision
!       use sll_splines
!       sll_real64                   :: interpolated_function_1D
!       sll_real64, intent(in)       :: x
!       type(sll_spline_1D), pointer :: spline_obj
!     end function interpolated_function_1D
!  end interface

  interface gauss_legendre_integrate_1D
     module procedure gauss_legendre_integral_1D !,gauss_legendre_integral_interpolated_1D
  end interface

contains

  ! object macro to put all the gauss points and weights in a single place.
  ! The fact that we have such case statement in a function that could be 
  ! used inside a critical loop is a flaw. Eventually we should split
  ! the integrators into individual cases. For now, we leave this as is
  ! due to the advantage of the simplified interface.

#define SELECT_CASES  \
    case(1); \
       xk(1)  = 0.0_f64 ; wk(1)  =  2.0_f64;\
    case (2); \
       xk(1)  = -1.0_f64/sqrt(3.0_f64); wk(1)  =  1.0_f64;\
       xk(2)  =  1.0_f64/sqrt(3.0_f64); wk(2)  =  1.0_f64; \
    case (3);\
       xk(1)  = -sqrt(3.0_f64/5.0_f64); wk(1)  =  5.0_f64/9.0_f64;\
       xk(2)  =  0.0_f64;               wk(2)  =  8.0_f64/9.0_f64;\
       xk(3)  =  sqrt(3.0_f64/5.0_f64); wk(3)  =  5.0_f64/9.0_f64;\
    case (4);\
       xk(1)  = -sqrt((3.0_f64 + 2.0_f64*sqrt(6.0_f64/5.0_f64))/7.0_f64);\
       xk(2)  = -sqrt((3.0_f64 - 2.0_f64*sqrt(6.0_f64/5.0_f64))/7.0_f64);\
       xk(3)  =  sqrt((3.0_f64 - 2.0_f64*sqrt(6.0_f64/5.0_f64))/7.0_f64);\
       xk(4)  =  sqrt((3.0_f64 + 2.0_f64*sqrt(6.0_f64/5.0_f64))/7.0_f64);\
       wk(1)  =  (18.0_f64 - sqrt(30.0_f64))/36.0_f64;\
       wk(2)  =  (18.0_f64 + sqrt(30.0_f64))/36.0_f64;\
       wk(3)  =  (18.0_f64 + sqrt(30.0_f64))/36.0_f64;\
       wk(4)  =  (18.0_f64 - sqrt(30.0_f64))/36.0_f64;\
    case (5);\
       xk(1)  = -0.90617984593866374_f64; wk(1)  =  0.23692688505618875_f64;\
       xk(2)  = -0.53846931010568311_f64; wk(2)  =  0.47862867049936653_f64;\
       xk(3)  =  0.0_f64;                 wk(3)  =  0.568888888888888888_f64;\
       xk(4)  =  0.53846931010568311_f64; wk(4)  =  0.47862867049936653_f64;\
       xk(5)  =  0.90617984593866374_f64; wk(5)  =  0.23692688505618875_f64;\
    case (6);\
       xk(1)  = -0.93246951420315249_f64; wk(1)  =  0.17132449237917155_f64;\
       xk(2)  = -0.66120938646626437_f64; wk(2)  =  0.36076157304813761_f64;\
       xk(3)  = -0.23861918608319668_f64; wk(3)  =  0.46791393457269115_f64;\
       xk(4)  =  0.23861918608319668_f64; wk(4)  =  0.46791393457269115_f64;\
       xk(5)  =  0.66120938646626437_f64; wk(5)  =  0.36076157304813761_f64;\
       xk(6)  =  0.93246951420315249_f64; wk(6)  =  0.17132449237917155_f64;\
    case (7);\
       xk(1)  = -0.94910791234275871_f64; wk(1)  =  0.12948496616886976_f64;\
       xk(2)  = -0.74153118559939302_f64; wk(2)  =  0.27970539148927637_f64;\
       xk(3)  = -0.4058451513773973_f64;  wk(3)  =  0.38183005050511931_f64;\
       xk(4)  =  0._f64;                  wk(4)  =  0.41795918367346824_f64;\
       xk(5)  =  0.4058451513773973_f64;  wk(5)  =  0.38183005050511931_f64;\
       xk(6)  =  0.74153118559939302_f64; wk(6)  =  0.27970539148927637_f64;\
       xk(7)  =  0.94910791234275871_f64; wk(7)  =  0.12948496616886976_f64;\
    case (8);\
       xk(1)  = -0.96028985649753706_f64; wk(1)  =  0.10122853629037747_f64;\
       xk(2)  = -0.79666647741362673_f64; wk(2)  =  0.22238103445337337_f64;\
       xk(3)  = -0.52553240991632888_f64; wk(3)  =  0.31370664587788677_f64;\
       xk(4)  = -0.18343464249565036_f64; wk(4)  =  0.36268378337836155_f64;\
       xk(5)  =  0.18343464249565036_f64; wk(5)  =  0.36268378337836155_f64;\
       xk(6)  =  0.52553240991632888_f64; wk(6)  =  0.31370664587788677_f64;\
       xk(7)  =  0.79666647741362673_f64; wk(7)  =  0.22238103445337337_f64;\
       xk(8)  =  0.96028985649753706_f64; wk(8)  =  0.10122853629037747_f64;\
    case (9);\
       xk(1)  = -0.968160239507626_f64;   wk(1)  =  0.081274388361574_f64;\
       xk(2)  = -0.836031107326636_f64;   wk(2)  =  0.180648160694857_f64;\
       xk(3)  = -0.613371432700590_f64;   wk(3)  =  0.260610696402935_f64;\
       xk(4)  = -0.324253423403809_f64;   wk(4)  =  0.312347077040003_f64;\
       xk(5)  =  0.0_f64;                 wk(5)  =  0.330239355001260_f64;\
       xk(6)  =  0.324253423403809_f64;   wk(6)  =  0.312347077040003_f64;\
       xk(7)  =  0.613371432700590_f64;   wk(7)  =  0.260610696402935_f64;\
       xk(8)  =  0.836031107326636_f64;   wk(8)  =  0.180648160694857_f64;\
       xk(9)  =  0.968160239507626_f64;   wk(9)  =  0.081274388361574_f64;\
    case (10);  \
       xk(1)  = -0.973906528517172_f64;   wk(1)  =  0.066671344308688_f64;\
       xk(2)  = -0.865063366688985_f64;   wk(2)  =  0.149451349150581_f64;\
       xk(3)  = -0.679409568299024_f64;   wk(3)  =  0.219086362515982_f64;\
       xk(4)  = -0.433395394129247_f64;   wk(4)  =  0.269266719309996_f64;\
       xk(5)  = -0.148874338981631_f64;   wk(5)  =  0.295524224714753_f64;\
       xk(6)  =  0.148874338981631_f64;   wk(6)  =  0.295524224714753_f64;\
       xk(7)  =  0.433395394129247_f64;   wk(7)  =  0.269266719309996_f64;\
       xk(8)  =  0.679409568299024_f64;   wk(8)  =  0.219086362515982_f64;\
       xk(9)  =  0.865063366688985_f64;   wk(9)  =  0.149451349150581_f64;\
       xk(10) =  0.973906528517172_f64;   wk(10) =  0.066671344308688_f64;\
    case default;\
       print *, 'gauss_legendre_integration: ';\
       print *, 'degree of integration not implemented. Exiting...';\
       stop;

  !> @brief Gauss-Legendre Quadrature.
  !> @details To integrate the function f(x) 
  !> (real-valued and of a single, real-valued argument x)
  !> over the interval [a,b], we use the Gauss-Legendre formula 
  !> \f[ \int_{-1}^1 f(x)dx \approx \sum_{k=1}^{n} w_k f(x_k) \f]
  !> where n represents the desired number of Gauss points.
  !>
  !> the function maps the interval [-1,1] into the
  !> arbitrary interval [a,b].
  !>
  !> To be considered is to split this function into degree-specific
  !> functions to avoid the select statement.
  !> @param f the function to be integrated
  !> @param[in] a left-bound of the definition interval of f  
  !> @param[in] b right-bound of the definition interval of f 
  !> @param[in] n the desired number of Gauss points
  !> @return The value of the integral
  function gauss_legendre_integral_1D( f, a, b, n )
    intrinsic                       :: sqrt
    sll_real64                      :: gauss_legendre_integral_1D
    procedure(function_1d_legendre) :: f
    sll_real64, intent(in)          :: a
    sll_real64, intent(in)          :: b
    sll_int32,  intent(in)          :: n ! needs better name
    sll_real64, dimension(1:10)     :: xk
    sll_real64, dimension(1:10)     :: wk
    sll_int32                       :: k
    sll_real64                      :: x
    sll_real64                      :: c1
    sll_real64                      :: c2
    sll_real64                      :: ans
    xk(:) = 0.0_f64
    wk(:) = 0.0_f64
    select case(n)
       SELECT_CASES
    end select
    ans = 0.0
    ! need to map the interval [-1,1] into the interval [a,b]
    c1 = 0.5_f64*(b-a)
    c2 = 0.5_f64*(b+a)
    do k=1,n
       x = c1*xk(k) + c2
       ans = ans + f(x)*wk(k)
    end do
    gauss_legendre_integral_1D = c1*ans
  end function gauss_legendre_integral_1D

  ! Consider changing this into a function that simply receives as an
  ! argument the spline_1D object and internally calls the interpolate_value()
  ! function. This would have a simpler interface. Although, there could be
  ! some advantages to have the interpolating function parametrized also, like
  ! in this case.
  !---------------------------------------------------------------------------  
  !> @author 
  !> Routine Author Name and Affiliation.
  !
  ! DESCRIPTION:
  !> @brief Integrates a function represented by a spline object.
  !> @details The function f in this case is the spline interpolation function.
  !> It looks like this interface could be simplified and we could eliminate 
  !> the first parameter and pass only the spline object.
  !> The only reason to leave the interpolation function as an argument is if
  !> we find some compelling reason to parametrize the interpolation function as well.
  !
  ! REVISION HISTORY:
  ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
  !
  !> @param f the spline interpolation function
  !> @param spline a spline
  !> @param[in] a left-bound of the definition interval of f  
  !> @param[in] b right-bound of the definition interval of f 
  !> @param[in] n the desired number of Gauss points
  !> @return The value of the integral
  !---------------------------------------------------------------------------
!!$  function gauss_legendre_integral_interpolated_1D( f, spline, a, b, n )
!!$    intrinsic                        :: sqrt
!!$    sll_real64                       :: gauss_legendre_integral_interpolated_1D
!!$    procedure(interpolated_function_1D) :: f
!!$    type(sll_spline_1D), pointer        :: spline
!!$    sll_real64, intent(in)              :: a
!!$    sll_real64, intent(in)              :: b
!!$    sll_int32,  intent(in)              :: n ! needs better name
!!$    sll_real64, dimension(1:10)         :: xk
!!$    sll_real64, dimension(1:10)         :: wk
!!$    sll_int32                           :: k
!!$    sll_real64                          :: x
!!$    sll_real64                          :: c1
!!$    sll_real64                          :: c2
!!$    sll_real64                          :: ans
!!$    xk(:) = 0.0_f64
!!$    wk(:) = 0.0_f64
!!$    select case(n)
!!$       SELECT_CASES
!!$    end select
!!$    ans = 0.0
!!$    ! need to map the interval [-1,1] into the interval [a,b]
!!$    c1 = 0.5_f64*(b-a)
!!$    c2 = 0.5_f64*(b+a)
!!$    do k=1,n
!!$       x = c1*xk(k) + c2
!!$       ans = ans + f(x,spline)*wk(k)
!!$    end do
!!$    gauss_legendre_integral_interpolated_1D = c1*ans
!!$  end function gauss_legendre_integral_interpolated_1D



  !> @brief Returns a 2d array of size (2,npoints) containing gauss-legendre 
  !> points and weights in the interval [a,b].
  !> @param[in] npoints Number of gauss points.
  !> @param[in] a OPTIONAL Minimum value of the interval.
  !> @param[in] b OPTIONAL Maximun value of the interval.


  !> gauss_points(degree) returns a real 2D array with the values of the
  !> locations of the gaussian points 'x_k' in the [-1,1] interval and
  !> their corresponding weights 'w_k'. Each column of the answer array
  !> contains the pair (x_k, w_k). Optionally, the user may provide the
  !> endpoints for the desired interval [a,b] where the gauss points should
  !> be mapped.
  function gauss_legendre_points_and_weights(npoints, a, b ) result(xw)
    sll_int32, intent(in)              :: npoints
    sll_real64, intent(in), optional   :: a
    sll_real64, intent(in), optional   :: b
    sll_real64, dimension(2,1:npoints)  :: xw
    sll_real64, dimension(1:npoints)    :: xk
    sll_real64, dimension(1:npoints)    :: wk
    sll_real64                         :: c1
    sll_real64                         :: c2
    sll_int32                          :: k

    SLL_ASSERT( npoints >= 1 )
    
    xk(:) = 0.0_f64
    wk(:) = 0.0_f64
    
    ! fill out the xk and wk arrays.
    select case(npoints)
       SELECT_CASES
    end select
    
    if(present(a) .and. present(b)) then
       ! need to map the interval [-1,1] into the interval [a,b].
       c1 = 0.5_f64*(b-a)
       c2 = 0.5_f64*(b+a)
       do k=1,npoints
          xw(1,k) = c1*xk(k) + c2
          xw(2,k) = wk(k)*c1
       end do
    else ! use default values in the [-1,1] interval
       xw(1,1:npoints) = xk(1:npoints)
       xw(2,1:npoints) = wk(1:npoints)
    end if

  end function gauss_legendre_points_and_weights

end module gauss_legendre_integration
