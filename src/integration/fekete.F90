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

contains

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
  !> @param[IN]  pxy array of dimesion (2,3) containg the coordinates 
  !>             of the edges of the triangle
  !> @param[OUT] wxy array of dimesion (3,10) containg the fekete points
  !>             and weights using rule 1. wx(1,:) contains the x coordinates
  !>             of the fekete points, wx(2, :) the y coordinates and
  !>             wx(3, :) contains the associated weights.
  function fekete_points_and_weights(pxy) result(wxy)
    sll_real64, dimension(2,3), intent(in) :: pxy
    sll_real64, dimension(3,10)            :: wxy
    ! Other variables
    sll_real64, dimension(3,3) :: ref_wxy
    sll_real64, dimension(2)   :: v1
    sll_real64, dimension(2)   :: v2
    sll_int32                  :: i
    sll_int32                  :: next
    sll_real64                 :: temp_x
    sll_real64                 :: temp_y
    
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
!    SLL_ALLOCATE(ref_wxy(3,3),ierr)
    ref_wxy(1:3, 1) = (/ 1./3._f64,    1./3._f64,   0.45_f64 /)
    ref_wxy(1:3, 2) = (/    0._f64,       0._f64, 1./60._f64 /)  
    ref_wxy(1:3, 3) = (/    0._f64, 0.2763932023_f64, 1./12._f64 /)

    ! Initialiting array containing points and weigths
    !SLL_ALLOCATE(wxy(3,10), ierr)
    wxy(:,:) = 0.0_f64
    next = 1

    ! Defining Fekete points parting from first reference edge
    !! v1 = Vector(p1, p2)
    v1(1) = pxy(1, 2) - pxy(1, 1)
    v1(2) = pxy(2, 2) - pxy(2, 1)
    !! v2 = Vector(p1, p3)
    v2(1) = pxy(1, 3) - pxy(1, 1)
    v2(2) = pxy(2, 3) - pxy(2, 1)

    do i = 1,3
       temp_x = v1(1)*ref_wxy(1, i) + v2(1)*ref_wxy(2, i) + pxy(1, 1)
       temp_y = v1(2)*ref_wxy(1, i) + v2(2)*ref_wxy(2, i) + pxy(2, 1)
       if (.not.( ANY( (wxy(1,:).eq.temp_x).and.(wxy(2,:).eq.temp_y) ))) then
          wxy(1, next) = temp_x
          wxy(2, next) = temp_y
          wxy(3, next) = ref_wxy(3, i)
          next = next + 1
       end if

       temp_x = v1(1)*ref_wxy(2, i) + v2(1)*ref_wxy(1, i) + pxy(1, 1)
       temp_y = v1(2)*ref_wxy(2, i) + v2(2)*ref_wxy(1, i) + pxy(2, 1)

       if (.not.( ANY( (wxy(1,:).eq.temp_x).and.(wxy(2,:).eq.temp_y) ))) then
          wxy(1, next) = temp_x
          wxy(2, next) = temp_y
          wxy(3, next) = ref_wxy(3, i)
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
       temp_x = v1(1)*ref_wxy(1, i) + v2(1)*ref_wxy(2, i) + pxy(1, 2)
       temp_y = v1(2)*ref_wxy(1, i) + v2(2)*ref_wxy(2, i) + pxy(2, 2)
       if (.not.( ANY( (wxy(1,:).eq.temp_x).and.(wxy(2,:).eq.temp_y) ))) then
          wxy(1, next) = temp_x
          wxy(2, next) = temp_y
          wxy(3, next) = ref_wxy(3, i)
          next = next + 1
       end if

       temp_x = v1(1)*ref_wxy(2, i) + v2(1)*ref_wxy(1, i) + pxy(1, 2)
       temp_y = v1(2)*ref_wxy(2, i) + v2(2)*ref_wxy(1, i) + pxy(2, 2)

       if (.not.( ANY( (wxy(1,:).eq.temp_x).and.(wxy(2,:).eq.temp_y) ))) then
          wxy(1, next) = temp_x
          wxy(2, next) = temp_y
          wxy(3, next) = ref_wxy(3, i)
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
       temp_x = v1(1)*ref_wxy(1, i) + v2(1)*ref_wxy(2, i) + pxy(1, 3)
       temp_y = v1(2)*ref_wxy(1, i) + v2(2)*ref_wxy(2, i) + pxy(2, 3)
       if (.not.( ANY( (wxy(1,:).eq.temp_x).and.(wxy(2,:).eq.temp_y) ))) then
          wxy(1, next) = temp_x
          wxy(2, next) = temp_y
          wxy(3, next) = ref_wxy(3, i)
          next = next + 1
       end if

       temp_x = v1(1)*ref_wxy(2, i) + v2(1)*ref_wxy(1, i) + pxy(1, 3)
       temp_y = v1(2)*ref_wxy(2, i) + v2(2)*ref_wxy(1, i) + pxy(2, 3)

       if (.not.( ANY( (wxy(1,:).eq.temp_x).and.(wxy(2,:).eq.temp_y) ))) then
          wxy(1, next) = temp_x
          wxy(2, next) = temp_y
          wxy(3, next) = ref_wxy(3, i)
          next = next + 1
       end if
    end do

  end function fekete_points_and_weights


! def main():
!     !Defining the vertices of the triangle
!     p1 = Point(0,0)
!     p2 = Point(1,0)
!     p3 = Point(0,1)

!     [fekPts, fekWei] = fekete3(p1, p2, p3)

!     for i in range(len(fekPts)) :
!         print '{0:2d} {1:6f} {2:6f} {3:6f}'.format(i+1, fekPts[i].X, fekPts[i].Y, fekWei[i])



end module fekete_integration
