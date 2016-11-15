!> @ingroup quadrature
!> @author Laura S. Mendoza
!> @brief Gaussian quadrature rule for a triangle
!> @details
!> This module contains the Gaussian quadrature rule adapted to triangles.
!> Some references for more information:
!>  [1] D. A. Dunavant, High degree efficient symmetrical Gaussian
!>        quadrature rules for the triangle, Int. J.
!>        Num. Meth. Engng, 21(1985):1129-1148.
!>  [2] http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
module sll_m_gauss_triangle_integration

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

  !---------------------------------------------------------------------------
  !> @brief returns the values of the points and weights on a reference triangle
  !> @details Returns the coordinates of the sll_m_gaussian points as well as their
  !> weights. They are computed in the reference tiangle :
  !>    |
  !>    1  3
  !>    |  |\
  !>    |  | \
  !>    |  |  \
  !>    |  |   \
  !>    |  |    \
  !>    0  1-----2
  !>    |
  !>    +--0-----1-->
  !> @param[IN]  order integer the order of the rule
  !> @param[OUT] points_xyw(order) real sll_m_gaussian quadrature points and weights
  subroutine gaussian_subrule(order, points_xyw)
    sll_int32,  intent(in)  :: order
    sll_real64, intent(out) :: points_xyw(3, order+1)

    if ((order .gt. 3) .or. (order .lt. 3)) then
       print *, "Error in gaussian_subrule(): Wrong order (", order, ")"
       STOP
    end if

    if (order .eq. 3) then
       points_xyw(1:2, 1) =   2._f64/15._f64
       points_xyw(  1, 2) =  11._f64/15._f64
       points_xyw(  2, 2) =   2._f64/15._f64
       points_xyw(1:2, 3) =   1._f64/ 3._f64
       points_xyw(  1, 4) =   2._f64/15._f64
       points_xyw(  2, 4) =  11._f64/15._f64
       ! Weights :
       points_xyw(3, 1:2) =  25._f64/48._f64
       points_xyw(3,   3) = -27._f64/48._f64
       points_xyw(3,   4) =  25._f64/48._f64
    end if
    
  end subroutine gaussian_subrule


    !---------------------------------------------------------------------------
  !> @brief maps T3 reference points to physical points.
  !> @details This code was first written by John Burkardt and is available
  !> online under the GNU LGPL license.
  !>    Given the vertices of an order 3 physical triangle and a point
  !>    (XSI,ETA) in the reference triangle, the routine computes the value
  !>    of the corresponding image point (X,Y) in physical space.
  !>
  !>    This routine is also appropriate for an order 4 triangle,
  !>    as long as the fourth node is the centroid of the triangle.
  !>
  !>    This routine may also be appropriate for an order 6
  !>    triangle, if the mapping between reference and physical space
  !>    is linear.  This implies, in particular, that the sides of the
  !>    image triangle are straight and that the "midside" nodes in the
  !>    physical triangle are literally halfway along the sides of
  !>    the physical triangle.
  !>
  !>  Reference Element T3:
  !>
  !>    |
  !>    1  3
  !>    |  |\
  !>    |  | \
  !>    S  |  \
  !>    |  |   \
  !>    |  |    \
  !>    0  1-----2
  !>    |
  !>    +--0--R--1-->
  !>
  !> @param[IN] node_xy(2,3) integer the coordinates of the vertices.
  !> The vertices are assumed to be the images of (0,0), (1,0) and (0,1).
  !> @param[IN] n integer the number of objects to transform
  !> @param[IN] ref(2,n) points in the reference triangle
  !> @param[OUT] phy(2,n) corresponding points in the physical triangle
  subroutine reference_to_physical_t3(node_xy, n, ref, phy)
    sll_int32,  intent(in) ::  n
    sll_real64, intent(in) :: node_xy(2,3)
    sll_real64, intent(in) :: ref(2,n)
    sll_real64, dimension(2,n), intent(out) :: phy
    sll_int32 ::  i

    do i = 1, 2
       phy(i,1:n) = node_xy(i,1) *(1.0_f64 - ref(1,1:n) - ref(2,1:n)) &
            + node_xy(i,2) * ref(1,1:n) &
            + node_xy(i,3) * ref(2,1:n)
    end do
  end subroutine reference_to_physical_t3


  !---------------------------------------------------------
  !> @brief Gives the sll_m_gaussian points coordinates and associated weights
  !> for a certain order in a given triangle
  !> @details 
  !> @param[IN]  node_xy2 array of dimesion (2,3) containg the coordinates
  !>             of the edges of the triangle
  !> @param[IN] order integer quadrature rule
  !> @return     xyw array of dimesion (3,n) containg the sll_m_gaussian points
  !>             and weights using the rule number given in parameter.
  !>             xyw(1,:) contains the x coordinates
  !>             of the sll_m_gaussian points, xyw(2, :) the y coordinates and
  !>             xyw(3, :) contains the associated weights.
  function gauss_triangle_points_and_weights(node_xy2, order) result(xyw)
    sll_real64, dimension(2,3), intent(in) :: node_xy2
    sll_int32,  intent(in) :: order
    sll_int32  :: ierr
    sll_real64, allocatable :: points_xyw(:,:)
    sll_real64, allocatable :: xy2(:,:)
    sll_real64, allocatable :: xyw(:,:)

    SLL_ALLOCATE(points_xyw(1:3,1:(order + 1)), ierr)
    SLL_ALLOCATE(xy2(1:2,1:(order + 1)), ierr)
    SLL_ALLOCATE(xyw(1:3,1:(order + 1)), ierr)

    call gaussian_subrule(order, points_xyw)
    xyw(3, :) = points_xyw(3, :)
    
    call reference_to_physical_t3 (node_xy2, order + 1, points_xyw(1:2,:), xy2)
    xyw(1:2, :) = xy2

    SLL_DEALLOCATE_ARRAY(points_xyw,  ierr)
    SLL_DEALLOCATE_ARRAY(xy2, ierr)
  end function gauss_triangle_points_and_weights
end module sll_m_gauss_triangle_integration
