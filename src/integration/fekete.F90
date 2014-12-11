!> @ingroup integration
!> @author Laura S. Mendoza
!> @brief Fekete quadrature rules for a triangle
!> @details 
!> This module contains the Fekete quadrature rule
!> adapted to triangles. So far the quadrature is only
!> for the rule 1 (order 10) but can be implemented for higher orders.
!> The module should incorporate integrations but for the
!> moment only contains the quadrature points and weights.
module fekete_integration

#include "sll_working_precision.h"
#include "sll_assert.h"

implicit none



#ifndef DOXYGEN_SHOULD_SKIP_THIS
abstract interface
   !> 2d real function
   function function_2D(x, y)
      use sll_working_precision ! can't pass a header file because the
      ! preprocessor prevents double inclusion.
      ! This is very rare.
      sll_real64             :: function_2D
      sll_real64, intent(in) :: x
      sll_real64, intent(in) :: y
    end function function_2D
end interface
#endif


!> Integrate numerically with Fekete quadrature formula
interface fekete_integrate
  module procedure fekete_integral
end interface


contains
  !---------------------------------------------------------------------------
  !> @brief Fekete quadrature rule over a triangle
  !> @details To integrate the function \f$ f(x) \f$
  !> (real-valued and over a triangle) we use the Fekete formula 
  !> \f[ \int_{\Omega} f(x)dx \approx \sum_{k=1}^{N} w_k f(x_k) \f]
  !> The only quadrature rule possible for now is 1 (10 points)
  !> @param[in]  f the function to be integrated
  !> @param[in]  pxy array of dimesion (2,3) containg the coordinates 
  !>             of the edges of the triangle
  !> @return The value of the integral
  function fekete_integral( f, pxy )
    sll_real64                :: fekete_integral
    procedure(function_2D)    :: f
    sll_real64, dimension(2, 3), intent(in) :: pxy
    sll_real64, dimension(3,10)             :: xyw
    sll_int32 :: k
    sll_int32 :: N

    ! Initialiting the fekete points and weigths ......
    xyw(:,:) = 0.
    xyw = fekete_points_and_weights(pxy)

    N = 10
    fekete_integral = 0.
    do k=1,N
       fekete_integral = fekete_integral + f(xyw(1,k), xyw(2,k))*xyw(3,k)
    end do

  end function fekete_integral


  !---------------------------------------------------------------------------
  !> @brief Function to compute rule 1 fekete points and weights on triangle
  !> @details Returns a 2d array of size (2,10) containing the fekete
  !> points and weights (rule 1) of a triangle which edges' coordinates are
  !> defined on px and py.
  !> Reference:
  !> Mark Taylor, Beth Wingate, Rachel Vincent,
  !> An Algorithm for Computing Fekete Points in the Triangle,
  !> SIAM Journal on Numerical Analysis,
  !> Volume 38, Number 5, 2000, pages 1707-1720.
  !> @param[in]  pxy array of dimesion (2,3) containg the coordinates 
  !>             of the edges of the triangle
  !> @param[out] xyw array of dimesion (3,10) containg the fekete points
  !>             and weights using rule 1. wx(1,:) contains the x coordinates
  !>             of the fekete points, wx(2, :) the y coordinates and
  !>             wx(3, :) contains the associated weights.
  function fekete_points_and_weights(pxy) result(xyw)
    sll_real64, dimension(2, 3), intent(in) :: pxy
    sll_real64, dimension(3,10)             :: xyw
    ! Other variables
    sll_real64, dimension(3,3) :: ref_xyw
    sll_real64, dimension(2)   :: v1
    sll_real64, dimension(2)   :: v2
    sll_int32                  :: i
    sll_int32                  :: next
    sll_real64                 :: temp_x
    sll_real64                 :: temp_y
    logical, dimension(10)        :: cond1
    logical, dimension(10)        :: cond2
    logical, dimension(10)        :: cond3
    
    ! Coordinates of base orbits points and weights in the following reference triangle
    !    |
    !    1  3
    !    |  |\
    !    |  | \
    !    |  |  \
    !    |  |   \
    !    |  |    \
    !    0  1-----2
    !    |
    !    +--0-----1-->
!    SLL_ALLOCATE(ref_xyw(3,3),ierr)
    ref_xyw(1:3, 1) = (/ 1./3._f64,    1./3._f64,   0.45_f64 /)
    ref_xyw(1:3, 2) = (/    0._f64,       0._f64, 1./60._f64 /)  
    ref_xyw(1:3, 3) = (/    0._f64, 0.2763932023_f64, 1./12._f64 /)

    ! Initialiting array containing points and weigths
    !SLL_ALLOCATE(xyw(3,10), ierr)
    xyw(:,:) = 0.0_f64
    next = 1

    ! Defining Fekete points parting from first reference edge
    !! v1 = Vector(p1, p2)
    v1(1) = pxy(1, 2) - pxy(1, 1)
    v1(2) = pxy(2, 2) - pxy(2, 1)
    !! v2 = Vector(p1, p3)
    v2(1) = pxy(1, 3) - pxy(1, 1)
    v2(2) = pxy(2, 3) - pxy(2, 1)

    do i = 1,3
       temp_x = v1(1)*ref_xyw(1, i) + v2(1)*ref_xyw(2, i) + pxy(1, 1)
       temp_y = v1(2)*ref_xyw(1, i) + v2(2)*ref_xyw(2, i) + pxy(2, 1)

       cond1 = abs(xyw(1,:)-temp_x).le.0.1E-14
       cond2 = abs(xyw(2,:)-temp_y).le.0.1E-14
       cond3 = abs(xyw(3,:)-ref_xyw(3,i)).le.0.1E-14
       if ( .not.( ANY(cond1.and.cond2.and.cond3) ) ) then
          xyw(1, next) = temp_x
          xyw(2, next) = temp_y
          xyw(3, next) = ref_xyw(3, i)
          next = next + 1
       end if

       temp_x = v1(1)*ref_xyw(2, i) + v2(1)*ref_xyw(1, i) + pxy(1, 1)
       temp_y = v1(2)*ref_xyw(2, i) + v2(2)*ref_xyw(1, i) + pxy(2, 1)
       cond1 = abs(xyw(1,:)-temp_x).le.0.1E-14
       cond2 = abs(xyw(2,:)-temp_y).le.0.1E-14
       cond3 = abs(xyw(3,:)-ref_xyw(3,i)).le.0.1E-14
       if ( .not.( ANY(cond1.and.cond2.and.cond3) ) ) then
          xyw(1, next) = temp_x
          xyw(2, next) = temp_y
          xyw(3, next) = ref_xyw(3, i)
          next = next + 1
       end if
    end do

    ! Defining Fekete points parting from p2
    !! v1 = Vector(p2, p1)
    v1(1) = pxy(1, 1) - pxy(1, 2)
    v1(2) = pxy(2, 1) - pxy(2, 2)
    !! v2 = Vector(p2, p3)
    v2(1) = pxy(1, 3) - pxy(1, 2)
    v2(2) = pxy(2, 3) - pxy(2, 2)

    do i = 1,3 

       temp_x = v1(1)*ref_xyw(1, i) + v2(1)*ref_xyw(2, i) + pxy(1, 2)
       temp_y = v1(2)*ref_xyw(1, i) + v2(2)*ref_xyw(2, i) + pxy(2, 2)
       cond1 = abs(xyw(1,:)-temp_x).le.0.1E-14
       cond2 = abs(xyw(2,:)-temp_y).le.0.1E-14
       cond3 = abs(xyw(3,:)-ref_xyw(3,i)).le.0.1E-14
       if ( .not.( ANY(cond1.and.cond2.and.cond3) ) ) then
          xyw(1, next) = temp_x
          xyw(2, next) = temp_y
          xyw(3, next) = ref_xyw(3, i)
          next = next + 1
       end if

       temp_x = v1(1)*ref_xyw(2, i) + v2(1)*ref_xyw(1, i) + pxy(1, 2)
       temp_y = v1(2)*ref_xyw(2, i) + v2(2)*ref_xyw(1, i) + pxy(2, 2)
       cond1 = abs(xyw(1,:)-temp_x).le.0.1E-14
       cond2 = abs(xyw(2,:)-temp_y).le.0.1E-14
       cond3 = abs(xyw(3,:)-ref_xyw(3,i)).le.0.1E-14
       if ( .not.( ANY(cond1.and.cond2.and.cond3) ) ) then
          xyw(1, next) = temp_x
          xyw(2, next) = temp_y
          xyw(3, next) = ref_xyw(3, i)
          next = next + 1
       end if
    end do

    ! Defining Fekete points parting from p3
    !! v1 = Vector(p3, p2)
    v1(1) = pxy(1, 2) - pxy(1, 3)
    v1(2) = pxy(2, 2) - pxy(2, 3)
    !! v2 = Vector(p3, p1)
    v2(1) = pxy(1, 1) - pxy(1, 3)
    v2(2) = pxy(2, 1) - pxy(2, 3)

    do i = 1,3 
       temp_x = v1(1)*ref_xyw(1, i) + v2(1)*ref_xyw(2, i) + pxy(1, 3)
       temp_y = v1(2)*ref_xyw(1, i) + v2(2)*ref_xyw(2, i) + pxy(2, 3)
       cond1 = abs(xyw(1,:)-temp_x).le.0.1E-14
       cond2 = abs(xyw(2,:)-temp_y).le.0.1E-14
       cond3 = abs(xyw(3,:)-ref_xyw(3,i)).le.0.1E-14
       if ( .not.( ANY(cond1.and.cond2.and.cond3) ) ) then
          xyw(1, next) = temp_x
          xyw(2, next) = temp_y
          xyw(3, next) = ref_xyw(3, i)
          next = next + 1
       end if

       temp_x = v1(1)*ref_xyw(2, i) + v2(1)*ref_xyw(1, i) + pxy(1, 3)
       temp_y = v1(2)*ref_xyw(2, i) + v2(2)*ref_xyw(1, i) + pxy(2, 3)
       cond1 = abs(xyw(1,:)-temp_x).le.0.1E-14
       cond2 = abs(xyw(2,:)-temp_y).le.0.1E-14
       cond3 = abs(xyw(3,:)-ref_xyw(3,i)).le.0.1E-14
       if ( .not.( ANY(cond1.and.cond2.and.cond3) ) ) then
          xyw(1, next) = temp_x
          xyw(2, next) = temp_y
          xyw(3, next) = ref_xyw(3, i)
          next = next + 1
       end if
    end do

  end function fekete_points_and_weights


end module fekete_integration
