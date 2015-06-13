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
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"

use sll_hex_meshes

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
    sll_real64, dimension(2, 3), intent(in) :: pxy
    sll_real64                  :: fekete_integral
    procedure(function_2D)      :: f
    sll_real64, dimension(3,10) :: xyw
    sll_real64, dimension(2)    :: v1
    sll_real64, dimension(2)    :: v2
    sll_real64, dimension(2)    :: v3
    sll_real64 :: a
    sll_real64 :: b
    sll_real64 :: c
    sll_real64 :: p
    sll_real64 :: area
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

    ! Computing the area of the triangle
    ! v1 = Vector(p1, p2)
    v1(1) = pxy(1, 2) - pxy(1, 1)
    v1(2) = pxy(2, 2) - pxy(2, 1)
    a = sqrt(v1(1)*v1(1) + v1(2)*v1(2))
    ! v2 = Vector(p1, p3)
    v2(1) = pxy(1, 3) - pxy(1, 1)
    v2(2) = pxy(2, 3) - pxy(2, 1)
    b = sqrt(v2(1)*v2(1) + v2(2)*v2(2))
    ! v3 = Vector(p2, p3)
    v3(1) = pxy(1, 3) - pxy(1, 2)
    v3(2) = pxy(2, 3) - pxy(2, 2)
    c = sqrt(v3(1)*v3(1) + v3(2)*v3(2))
    ! Computing demi-perimeter
    p = 0.5*(a+b+c)
    ! area
    area = sqrt(p*(p-a)*(p-b)*(p-c))

    fekete_integral = fekete_integral * area
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
  !> @return     xyw array of dimesion (3,10) containg the fekete points
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
    logical, dimension(10)     :: cond1
    logical, dimension(10)     :: cond2
    logical, dimension(10)     :: cond3

    ! Coordinates of base orbits points and weights
    ! in the following reference triangle
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

    ref_xyw(1:3, 1) = [ 1.0_f64/3.0_f64,  1.0_f64/3.0_f64,        0.45_f64 ]
    ref_xyw(1:3, 2) = [ 0.0_f64        ,          0.0_f64, 1.0_f64/60._f64 ]
    ref_xyw(1:3, 3) = [ 0.0_f64        , 0.2763932023_f64, 1.0_f64/12._f64 ]

    ! Initialiting array containing points and weigths
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


  !---------------------------------------------------------------------------
  !> @brief Computes fekete points coordinates on a hex-mesh
  !> @details Goes through each cell (triangle) of the hex mesh and computes
  !> the fekete points on each one. It fills in the knots table and the
  !> connectivity table LM.
  !> @param[in]  rule integer for the fekete quadrature rule
  !> @param[in]  mesh hex_mesh on which we want to evaluate quadrature points
  !> @param[out] knots real 2D array to be filled in with fekete points+weights
  !> @param[out] LM connectivity array to pass from local to global indices
  subroutine initialize_knots_hexmesh( rule, &
                                       mesh, &
                                       knots, &
                                       LM)

    sll_int32, intent(in) :: rule
    type(sll_hex_mesh_2d), intent(in), pointer :: mesh
    sll_real64, dimension(:,:), intent(out)    :: knots
    sll_int32,  dimension(:,:), intent(out)    :: LM
    logical,    dimension(:),   allocatable    :: cond1
    logical,    dimension(:),   allocatable    :: cond2
    logical,    dimension(:),   allocatable    :: cond3
    sll_int32,  dimension(:),   allocatable    :: indices
    sll_real64, dimension(3,10)                :: fekete_tri
    sll_real64, dimension(2, 3)                :: vertices_coo
    sll_int32,  dimension(1)                   :: temp
    
    sll_int32  :: i
    sll_int32  :: ierr
    sll_int32  :: ntri
    sll_int32  :: next
    sll_int32  :: num_fekete
    sll_int32  :: vertex1
    sll_int32  :: vertex2
    sll_int32  :: vertex3
    sll_real64 :: cntrx
    sll_real64 :: cntry

    if (rule .eq. 1) then

       ! When rule = 1, then there is a fekete point on each vertex of the mesh
       ! plus two fekete points on each edge and one on the middle of each cell
       num_fekete = 2*mesh%num_edges + mesh%num_pts_tot + mesh%num_triangles

       ! Allocation of helper tables:
       SLL_ALLOCATE(cond1(num_fekete), ierr)
       SLL_ALLOCATE(cond2(num_fekete), ierr)
       SLL_ALLOCATE(cond3(num_fekete), ierr)
       SLL_ALLOCATE(indices(num_fekete), ierr)

       ! Initialitazion tables and counter
       knots(1:3, 1:num_fekete) = 0._f64
       LM(1:mesh%num_triangles, 1:10) = 0
       indices(1:num_fekete) = (/ (i, i = 1,num_fekete) /)
       next = 1

       !We go through each cell (triangle) of the hexagon
       do ntri = 1,mesh%num_triangles

          !We fetch the vertices's coordinates for each triangle
          cntrx = mesh%center_cartesian_coord(1, ntri)
          cntry = mesh%center_cartesian_coord(2, ntri)
          call get_cell_vertices_index(cntrx, cntry, mesh, vertex1, vertex2, vertex3)
          vertices_coo(:, 1) = mesh%cartesian_coord(:, vertex1)
          vertices_coo(:, 2) = mesh%cartesian_coord(:, vertex2)
          vertices_coo(:, 3) = mesh%cartesian_coord(:, vertex3)

          !We compute the fekte points on that triangle :
          fekete_tri = fekete_points_and_weights(vertices_coo)

          do i = 1, 10
             !We add them if they are not already on the table:
             cond1 = abs(knots(1,:)-fekete_tri(1,i)).le.0.1E-14
             cond2 = abs(knots(2,:)-fekete_tri(2,i)).le.0.1E-14
             cond3 = abs(knots(3,:)-fekete_tri(3,i)).le.0.1E-14
             if ( .not.( ANY(cond1.and.cond2.and.cond3) ) ) then
                knots(:, next) = fekete_tri(:, i)
                LM(ntri, i) = next
                next = next + 1
             else
                temp = MAXLOC(indices, MASK = cond1.and.cond2.and.cond3)
                LM(ntri, i) = temp(1)
             end if
          end do
       end do
    end if

    ! do ntri = 1,mesh%num_triangles
    !    print *, ntri
    !    print *, LM(ntri, :)
    ! end do
  end subroutine initialize_knots_hexmesh

  !---------------------------------------------------------------------------
  !> @brief Writes fekete points coordinates of a hex-mesh reference triangle
  !> @details Takes the reference triangle of a hexmesh and computes the
  !> fekete points on it. Then it writes the results in a file following
  !> CAID nomenclature. 
  !> Output file : quadrature.txt
  !> @param[in]  rule integer for the fekete quadrature rule
  subroutine write_quadrature(rule)
    sll_int32, intent(in)       :: rule
    sll_int32                   :: out_unit
    character(len=14), parameter :: name = "quadrature.txt"
    sll_real64, dimension(2, 3) :: ref_pts
    sll_real64, dimension(3,10) :: quad_pw
    sll_int32  :: num_fek
    sll_int32  :: i
    sll_real64 :: x1, x2, x3
    sll_real64 :: y1, y2, y3
    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: w
    sll_real64 :: lambda1
    sll_real64 :: lambda2
    sll_int32  :: ierr
    ! Definition of reference triangle, such that:
    !    |
    !    1  3
    !    |  | \
    !    |  |  \2
    !    |  |   /   same cell that first cell of
    !    |  |  /    a simple hexagon of radius 1.
    !    |  | /
    !    0  1
    !    |
    !    +--0-----1-->
    ref_pts(:,1) = (/ 0._f64,          0.0_f64 /)
    ref_pts(:,2) = (/ 1./sqrt(3._f64), 0.5_f64 /)
    ref_pts(:,3) = (/ 0._f64,          1.0_f64 /)

    ! Computing fekete points on that triangle
    quad_pw = fekete_points_and_weights(ref_pts)

    call sll_new_file_id(out_unit, ierr)
    open (unit=out_unit,file=name,action="write",status="replace")
    

    if (rule .eq. 1) then
       num_fek = 10
    else
       SLL_WARNING( "write_quadrature", "Required rule is not implemented." )
       num_fek = 0
       STOP
    end if

    write(out_unit, "(i6)") num_fek
    
    do i=1,num_fek
       x = quad_pw(1,i)
       y = quad_pw(2,i)
       w = quad_pw(3,i)
       ! Transformation to barycentric coordinates
       ! for more info read the wiki on barycentric coordinates
       ! "Conversion between barycentric and Cartesian coordinates"
       x1 = ref_pts(1, 1) ; y1 = ref_pts(2, 1)
       x2 = ref_pts(1, 2) ; y2 = ref_pts(2, 2)
       x3 = ref_pts(1, 3) ; y3 = ref_pts(2, 3)
       lambda1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3)) / ((y2-y3)*(x1-x3) + (x3-x2)*(y1-y3))
       lambda2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3)) / ((y2-y3)*(x1-x3) + (x3-x2)*(y1-y3))
       write(out_unit, "(2(g25.17,a,1x),(g25.17))") lambda1, ",", lambda2, ",", w
    end do
    close(out_unit)

  end subroutine write_quadrature

end module fekete_integration
