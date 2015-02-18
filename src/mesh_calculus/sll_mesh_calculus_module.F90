module sll_mesh_calculus_2d_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"

  use sll_cartesian_meshes
  use sll_triangular_meshes
  use sll_coordinate_transformation_2d_base_module
  use gauss_legendre_integration
  implicit none

  ! --------------------------------------------------------------------------
  !
  ! The mesh calculus module aims at providing the means of computing
  ! quantities like edge-lengths, areas, volumes and normals of cell elements 
  ! in a deformed mesh.
  !
  ! --------------------------------------------------------------------------

contains

  function cell_volume( T, ic, jc, integration_degree ) result(vol)
    intrinsic :: abs
    class(sll_coordinate_transformation_2d_base), pointer :: T
    sll_int32, intent(in) :: ic
    sll_int32, intent(in) :: jc
    sll_int32, intent(in) :: integration_degree
    sll_real64            :: vol   ! volume in physical space
    sll_real64 :: jac
    sll_real64, dimension(2,integration_degree) :: pts_g1 ! gauss-legendre pts.
    sll_real64, dimension(2,integration_degree) :: pts_g2 ! gauss-legendre pts.
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: max1
    sll_real64 :: min1
    sll_real64 :: max2
    sll_real64 :: min2
    !sll_real64 :: factor1
    !sll_real64 :: factor2
    sll_int32  :: i
    sll_int32  :: j
    class(sll_cartesian_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_cartesian_mesh()
    ! verify that the indices requested are within the logical mesh.
    SLL_ASSERT(ic <= m%num_cells1)
    SLL_ASSERT(jc <= m%num_cells2)

    vol = 0.0_f64

    eta1_min = m%eta1_min
    delta1   = m%delta_eta1
    eta2_min = m%eta2_min
    delta2   = m%delta_eta2

    ! This function carries out the integral of the jacobian evaluated on
    ! gauss-legendre points within the cell.
    ! The limits of integration are the limits of the cell in eta-space
    min1    = eta1_min + (ic-1)*delta1
    max1    = eta1_min + ic*delta1
    min2    = eta2_min + (jc-1)*delta2
    max2    = eta2_min + jc*delta2
    !factor1 = 0.5_f64*(max1-min1)
    !factor2 = 0.5_f64*(max2-min2) 
    pts_g1(:,:) = &
         gauss_legendre_points_and_weights(integration_degree, min1, max1)
    !gauss_points(integration_degree, min1, max1)
    pts_g2(:,:) = &
         gauss_legendre_points_and_weights(integration_degree, min2, max2)
    ! gauss_points(integration_degree, min2, max2)

    do j=1,integration_degree
       do i=1,integration_degree
          jac = T%jacobian(pts_g1(1,i),pts_g2(1,j))
          vol = vol + abs(jac)*pts_g1(2,i)*pts_g2(2,j)
       end do
    end do
    !vol = vol*factor1*factor2

  end function cell_volume

  ! length of the 'east' edge of the cell.
  function edge_length_eta1_plus( T, ic, jc, integration_degree ) result(len)
    intrinsic :: abs
    class(sll_coordinate_transformation_2d_base), pointer :: T
    sll_int32, intent(in) :: ic
    sll_int32, intent(in) :: jc
    sll_int32, intent(in) :: integration_degree
    sll_real64            :: len   ! length of edge in physical space
    sll_real64, dimension(2,2) :: jac_mat
!    sll_real64, dimension(2,integration_degree) :: pts_g1 ! gauss-legendre pts.
    sll_real64, dimension(2,integration_degree) :: pts_g2 ! gauss-legendre pts.
    sll_real64 :: eta1
    !sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    !sll_real64 :: max1
    !sll_real64 :: min1
    sll_real64 :: max2
    sll_real64 :: min2
    !sll_real64 :: factor1
    !sll_real64 :: factor2
    sll_int32  :: j
    sll_real64 :: x1_eta2  ! derivative of x1(eta1,eta2) with respect to eta2
    sll_real64 :: x2_eta2  ! derivative of x1(eta1,eta2) with respect to eta2
    class(sll_cartesian_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_cartesian_mesh()
    ! verify that the indices requested are within the logical mesh.
    SLL_ASSERT(ic <= m%num_cells1)
    SLL_ASSERT(jc <= m%num_cells2)

    len = 0.0_f64

!    eta1_min = T%mesh%eta1_min
    delta1   = m%delta_eta1
    eta2_min = m%eta2_min
    delta2   = m%delta_eta2

    ! The limits of integration are the limits of the cell in eta-space
    eta1    = ic*delta1 ! only difference with edge_length_eta1_minus function
    min2    = eta2_min + (jc-1)*delta2
    max2    = eta2_min + jc*delta2
    pts_g2(:,:) = &
         gauss_legendre_points_and_weights(integration_degree, min2, max2)
    
    do j=1,integration_degree
       ! this can be made more efficient if we could access directly each
       ! term of the jacobian matrix independently.
       jac_mat(:,:) = T%jacobian_matrix(eta1,pts_g2(1,j))
       x1_eta2 = jac_mat(1,2)
       x2_eta2 = jac_mat(2,2)
       len = len + sqrt(x1_eta2**2 + x2_eta2**2)*pts_g2(2,j)
    end do
    !len = len*factor2
  end function edge_length_eta1_plus

  ! length of the 'west' edge of the cell.
  function edge_length_eta1_minus( T, ic, jc, integration_degree ) result(len)
    intrinsic :: abs
    class(sll_coordinate_transformation_2d_base), pointer :: T
    sll_int32, intent(in) :: ic
    sll_int32, intent(in) :: jc
    sll_int32, intent(in) :: integration_degree
    sll_real64            :: len   ! length of edge in physical space
    sll_real64, dimension(2,2) :: jac_mat
!    sll_real64, dimension(2,integration_degree) :: pts_g1 ! gauss-legendre pts.
    sll_real64, dimension(2,integration_degree) :: pts_g2 ! gauss-legendre pts.
    sll_real64 :: eta1
    !sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    !sll_real64 :: max1
    !sll_real64 :: min1
    sll_real64 :: max2
    sll_real64 :: min2
    !    sll_real64 :: factor1
    !    sll_real64 :: factor2
    sll_int32  :: j
    sll_real64 :: x1_eta2  ! derivative of x1(eta1,eta2) with respect to eta2
    sll_real64 :: x2_eta2  ! derivative of x1(eta1,eta2) with respect to eta2
    class(sll_cartesian_mesh_2d), pointer :: m
    
    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_cartesian_mesh()

    ! verify that the indices requested are within the logical mesh.
    SLL_ASSERT(ic <= m%num_cells1)
    SLL_ASSERT(jc <= m%num_cells2)

    len = 0.0_f64
    
    !    eta1_min = T%mesh%eta1_min
    delta1   = m%delta_eta1
    eta2_min = m%eta2_min
    delta2   = m%delta_eta2

    ! The limits of integration are the limits of the cell in eta-space
    eta1    = (ic-1)*delta1
    min2    = eta2_min + (jc-1)*delta2
    max2    = eta2_min + jc*delta2
    !    factor1 = 0.5_f64*(max1-min1)
    !    factor2 = 0.5_f64*(max2-min2) 
    !    pts_g1(:,:) = gauss_points(integration_degree, min1, max1)
    pts_g2(:,:) = &
         gauss_legendre_points_and_weights(integration_degree, min2, max2)
    !gauss_points(integration_degree, min2, max2)
    
    do j=1,integration_degree
       ! this can be made more efficient if we could access directly each
       ! term of the jacobian matrix independently.
       jac_mat(:,:) = T%jacobian_matrix(eta1,pts_g2(1,j))
       x1_eta2 = jac_mat(1,2)
       x2_eta2 = jac_mat(2,2)
       len = len + sqrt(x1_eta2**2 + x2_eta2**2)*pts_g2(2,j)
    end do
    !len = len*factor2
  end function edge_length_eta1_minus

  ! length of the 'north' edge of the cell.
  function edge_length_eta2_plus( T, ic, jc, integration_degree ) result(len)
    intrinsic :: abs
    class(sll_coordinate_transformation_2d_base), pointer :: T
    sll_int32, intent(in) :: ic
    sll_int32, intent(in) :: jc
    sll_int32, intent(in) :: integration_degree
    sll_real64            :: len   ! length of edge in physical space
    sll_real64, dimension(2,2) :: jac_mat
    sll_real64, dimension(2,integration_degree) :: pts_g1 ! gauss-legendre pts.
!    sll_real64, dimension(2,integration_degree) :: pts_g2 ! gauss-legendre pts.
    sll_real64 :: eta2
    sll_real64 :: eta1_min
    !sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: max1
    sll_real64 :: min1
    sll_int32  :: i
    sll_real64 :: x1_eta1  ! derivative of x1(eta1,eta2) with respect to eta1
    sll_real64 :: x2_eta1  ! derivative of x1(eta1,eta2) with respect to eta1
    class(sll_cartesian_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_cartesian_mesh()

    ! verify that the indices requested are within the logical mesh.
    SLL_ASSERT(ic <= m%num_cells1)
    SLL_ASSERT(jc <= m%num_cells2)

    len = 0.0_f64

    eta1_min = m%eta1_min
    delta1   = m%delta_eta1
    delta2   = m%delta_eta2  ! is this used?

    ! The limits of integration are the limits of the cell in eta-space
    min1    = eta1_min + (ic-1)*delta1
    max1    = eta1_min + ic*delta1
    eta2    = jc*delta2 ! only difference with edge_length_eta2_minus function
    pts_g1(:,:) = &
         gauss_legendre_points_and_weights(integration_degree, min1, max1)
    !gauss_points(integration_degree, min1, max1)
    !    pts_g2(:,:) = gauss_points(integration_degree, min2, max2)
    
    do i=1,integration_degree
       ! this can be made more efficient if we could access directly each
       ! term of the jacobian matrix independently.
       jac_mat(:,:) = T%jacobian_matrix(pts_g1(1,i),eta2)
       x1_eta1 = jac_mat(1,1)
       x2_eta1 = jac_mat(2,1)
       len = len + sqrt(x1_eta1**2 + x2_eta1**2)*pts_g1(2,i)
    end do
    ! len = len*factor1
  end function edge_length_eta2_plus
  
  ! length of the 'south' edge of the cell.
  function edge_length_eta2_minus( T, ic, jc, integration_degree ) result(len)
    intrinsic :: abs
    class(sll_coordinate_transformation_2d_base), pointer :: T
    sll_int32, intent(in) :: ic
    sll_int32, intent(in) :: jc
    sll_int32, intent(in) :: integration_degree
    sll_real64            :: len   ! length of edge in physical space
    sll_real64, dimension(2,2) :: jac_mat
    sll_real64, dimension(2,integration_degree) :: pts_g1 ! gauss-legendre pts.
    !    sll_real64, dimension(2,integration_degree) :: pts_g2 ! gauss-legendre pts.
    sll_real64 :: eta2
    sll_real64 :: eta1_min
    !sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: max1
    sll_real64 :: min1
    !sll_real64 :: max2
    !sll_real64 :: min2
    !sll_real64 :: factor1
    !    sll_real64 :: factor2
    sll_int32  :: i
    sll_real64 :: x1_eta1  ! derivative of x1(eta1,eta2) with respect to eta1
    sll_real64 :: x2_eta1  ! derivative of x1(eta1,eta2) with respect to eta1
    class(sll_cartesian_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_cartesian_mesh()

    ! verify that the indices requested are within the logical mesh.
    SLL_ASSERT(ic <= m%num_cells1)
    SLL_ASSERT(jc <= m%num_cells2)
    
    len = 0.0_f64
    
    eta1_min = m%eta1_min
    delta1   = m%delta_eta1
    delta2   = m%delta_eta2  ! is this used?

    ! The limits of integration are the limits of the cell in eta-space
    min1    = eta1_min + (ic-1)*delta1
    max1    = eta1_min + ic*delta1
    eta2    = (jc-1)*delta2 !only difference with edge_length_eta2_plus function
    !    min2    = eta2_min + (jc-1)*delta2
    !    max2    = eta2_min + jc*delta2
    !factor1 = 0.5_f64*(max1-min1)
    !    factor2 = 0.5_f64*(max2-min2) 
    pts_g1(:,:) = &
         gauss_legendre_points_and_weights(integration_degree, min1, max1)
    !gauss_points(integration_degree, min1, max1)
    !    pts_g2(:,:) = gauss_points(integration_degree, min2, max2)
    
    do i=1,integration_degree
       ! this can be made more efficient if we could access directly each
       ! term of the jacobian matrix independently.
       jac_mat(:,:) = T%jacobian_matrix(pts_g1(1,i),eta2)
       x1_eta1 = jac_mat(1,1)
       x2_eta1 = jac_mat(2,1)
       len = len + sqrt(x1_eta1**2 + x2_eta1**2)*pts_g1(2,i)
    end do
    !len = len*factor1
  end function edge_length_eta2_minus
  
  ! integral of the normal vector over the 'east' edge of the cell.
  function normal_integral_eta1_plus( T,ic,jc,integration_degree ) result(res)
    class(sll_coordinate_transformation_2d_base), pointer :: T
    sll_int32, intent(in)     :: ic
    sll_int32, intent(in)     :: jc
    sll_int32, intent(in)     :: integration_degree
    sll_real64, dimension(2)  :: res
    
    sll_real64, dimension(2,2) :: inv_jac_mat ! inverse jacobian matrix
    !    sll_real64, dimension(2,integration_degree) :: pts_g1 ! gauss-legendre pts.
    sll_real64, dimension(2,integration_degree) :: pts_g2 ! gauss-legendre pts.
    sll_real64 :: eta1
    !sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    !sll_real64 :: max1
    !sll_real64 :: min1
    sll_real64 :: max2
    sll_real64 :: min2
    !sll_real64 :: factor1
    !sll_real64 :: factor2
    sll_int32  :: j
    sll_real64 :: eta1_x1  ! derivative of eta1(x1,x2) with respect to x1
    sll_real64 :: eta1_x2  ! derivative of eta1(x1,x2) with respect to x2
    sll_real64 :: edge_length
    class(sll_cartesian_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_cartesian_mesh()

    ! verify that the indices requested are within the logical mesh.
    SLL_ASSERT(ic <= m%num_cells1)
    SLL_ASSERT(jc <= m%num_cells2)
    
    !    eta1_min = T%mesh%eta1_min
    delta1   = m%delta_eta1
    eta2_min = m%eta2_min
    delta2   = m%delta_eta2

    ! The limits of integration are the limits of the cell in eta-space
    !    min1    = eta1_min + (ic-1)*delta1
    !    max1    = eta1_min + ic*delta1
    eta1    = ic*delta1 ! <- line differs w/ normal_integral_eta1_minus()
    min2    = eta2_min + (jc-1)*delta2
    max2    = eta2_min + jc*delta2
    !    factor1 = 0.5_f64*(max1-min1)
    !factor2 = 0.5_f64*(max2-min2) 
    !    pts_g1(:,:) = gauss_points(integration_degree, min1, max1)
    pts_g2(:,:) = &
         gauss_legendre_points_and_weights(integration_degree, min2, max2)
    !gauss_points(integration_degree, min2, max2)
    
    ! For efficiency, this code should be refactored. Consider:
    ! - adding direct access functions to the jacobian matrix elements and the
    !   inverse jacobian matrix elements
    ! - use macros to eliminate the massive code redundancy in this module's
    !   functions.
    ! - same macros can be used to improve a function call like the one next.
    edge_length = edge_length_eta1_plus( T, ic, jc, integration_degree )
    res(:) = 0.0_f64
    
    do j=1,integration_degree
       ! this can be made more efficient if we could access directly each
       ! term of the jacobian matrix independently.
       inv_jac_mat(:,:) = T%inverse_jacobian_matrix(eta1,pts_g2(1,j))
       eta1_x1 = inv_jac_mat(1,1)
       eta1_x2 = inv_jac_mat(2,1)
       SLL_ASSERT(T%jacobian(eta1,pts_g2(1,j)) > 0.0_f64)
       res(1) = res(1) + eta1_x1*pts_g2(2,j)
       res(2) = res(2) + eta1_x2*pts_g2(2,j)
    end do
    res(1) = res(1)*edge_length !res(1) = res(1)*factor2*edge_length 
    res(2) = res(2)*edge_length ! res(2) = res(2)*factor2*edge_length
  end function normal_integral_eta1_plus
  
  ! integral of the normal vector over the 'west' edge of the cell.
  function normal_integral_eta1_minus( T,ic,jc,integration_degree ) result(res)
    class(sll_coordinate_transformation_2d_base), pointer :: T
    sll_int32, intent(in)    :: ic
    sll_int32, intent(in)    :: jc
    sll_int32, intent(in)    :: integration_degree
    sll_real64, dimension(2) :: res
    
    sll_real64, dimension(2,2) :: inv_jac_mat ! inverse jacobian matrix
    sll_real64, dimension(2,integration_degree) :: pts_g2 ! gauss-legendre pts.
    sll_real64 :: eta1
    !sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    !sll_real64 :: max1
    !sll_real64 :: min1
    sll_real64 :: max2
    sll_real64 :: min2
    !sll_real64 :: factor1
    !sll_real64 :: factor2
    sll_int32  :: j
    sll_real64 :: eta1_x1  ! derivative of eta1(x1,x2) with respect to x1
    sll_real64 :: eta1_x2  ! derivative of eta1(x1,x2) with respect to x2
    sll_real64 :: edge_length
    class(sll_cartesian_mesh_2d), pointer :: m
    
    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_cartesian_mesh()

    ! verify that the indices requested are within the logical mesh.
    SLL_ASSERT(ic <= m%num_cells1)
    SLL_ASSERT(jc <= m%num_cells2)

    delta1   = m%delta_eta1
    eta2_min = m%eta2_min
    delta2   = m%delta_eta2
    
    ! The limits of integration are the limits of the cell in eta-space
    !    min1    = eta1_min + (ic-1)*delta1
    !    max1    = eta1_min + ic*delta1
    eta1    = (ic-1)*delta1 ! <- line differs w/ normal_integral_eta1_plus()
    min2    = eta2_min + (jc-1)*delta2
    max2    = eta2_min + jc*delta2
    !    factor1 = 0.5_f64*(max1-min1)
    !factor2 = 0.5_f64*(max2-min2) 
    !    pts_g1(:,:) = gauss_points(integration_degree, min1, max1)
    pts_g2(:,:) = &
         gauss_legendre_points_and_weights(integration_degree, min2, max2)
    !gauss_points(integration_degree, min2, max2)
    
    ! For efficiency, this code should be refactored. Consider:
    ! - adding direct access functions to the jacobian matrix elements and the
    !   inverse jacobian matrix elements
    ! - use macros to eliminate the massive code redundancy in this module's
    !   functions.
    ! - same macros can be used to improve a function call like the one next.
    edge_length = edge_length_eta1_minus( T, ic, jc, integration_degree )
    
    do j=1,integration_degree
       ! this can be made more efficient if we could access directly each
       ! term of the jacobian matrix independently.
       inv_jac_mat(:,:) = T%inverse_jacobian_matrix(eta1,pts_g2(1,j))
       eta1_x1 = inv_jac_mat(1,1)
       eta1_x2 = inv_jac_mat(2,1)
       SLL_ASSERT(T%jacobian(eta1,pts_g2(1,j)) > 0.0_f64)
       res(1) = res(1) + eta1_x1*pts_g2(2,j)
       res(2) = res(2) + eta1_x2*pts_g2(2,j)
    end do
    ! change of sign due to the direction in which the integral is done
    res(1) = -res(1)*edge_length 
    res(2) = -res(2)*edge_length 
!!$    res(1) = -res(1)*factor2*edge_length 
!!$    res(2) = -res(2)*factor2*edge_length 
  end function normal_integral_eta1_minus
  
  ! integral of the normal vector over the 'north' edge of the cell.
  function normal_integral_eta2_plus( T,ic,jc,integration_degree ) result(res)
    class(sll_coordinate_transformation_2d_base), pointer :: T
    sll_int32, intent(in)    :: ic
    sll_int32, intent(in)    :: jc
    sll_int32, intent(in)    :: integration_degree
    sll_real64, dimension(2) :: res
    
    sll_real64, dimension(2,2) :: inv_jac_mat ! inverse jacobian matrix
    sll_real64, dimension(2,integration_degree) :: pts_g1 ! gauss-legendre pts.
    !    sll_real64, dimension(2,integration_degree) :: pts_g2 ! gauss-legendre pts.
    sll_real64 :: eta2
    sll_real64 :: eta1_min
    !sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: max1
    sll_real64 :: min1
    !sll_real64 :: max2
    !sll_real64 :: min2
    !sll_real64 :: factor1
    !sll_real64 :: factor2
    sll_int32  :: i
    sll_real64 :: eta2_x1  ! derivative of eta2(x1,x2) with respect to x1
    sll_real64 :: eta2_x2  ! derivative of eta2(x1,x2) with respect to x2
    sll_real64 :: edge_length
    class(sll_cartesian_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_cartesian_mesh()

    ! verify that the indices requested are within the logical mesh.
    SLL_ASSERT(ic <= m%num_cells1)
    SLL_ASSERT(jc <= m%num_cells2)
    
    eta1_min = m%eta1_min
    delta1   = m%delta_eta1
    delta2   = m%delta_eta2
    
    ! The limits of integration are the limits of the cell in eta-space
    min1    = eta1_min + (ic-1)*delta1
    max1    = eta1_min + ic*delta1
    eta2    = jc*delta2 ! <- line differs w/ normal_integral_eta2_minus()
    !    min2    = eta2_min + (jc-1)*delta2
    !    max2    = eta2_min + jc*delta2
   ! factor1 = 0.5_f64*(max1-min1)
    !    factor2 = 0.5_f64*(max2-min2) 
    pts_g1(:,:) = gauss_legendre_points_and_weights(integration_degree, min1, max1)
    !gauss_points(integration_degree, min1, max1)
    !    pts_g2(:,:) = gauss_points(integration_degree, min2, max2)
    
    ! For efficiency, this code should be refactored. Consider:
    ! - adding direct access functions to the jacobian matrix elements and the
    !   inverse jacobian matrix elements
    ! - use macros to eliminate the massive code redundancy in this module's
    !   functions.
    ! - same macros can be used to improve a function call like the one next.
    edge_length = edge_length_eta2_plus( T, ic, jc, integration_degree )
    res(:) = 0.0_f64
    
    do i=1,integration_degree
       ! this can be made more efficient if we could access directly each
       ! term of the jacobian matrix independently.
       inv_jac_mat(:,:) = T%inverse_jacobian_matrix(pts_g1(1,i),eta2)
       eta2_x1 = inv_jac_mat(2,1)
       eta2_x2 = inv_jac_mat(2,2)
       SLL_ASSERT(T%jacobian(pts_g1(1,i),eta2) > 0.0_f64)
       res(1) = res(1) + eta2_x1*pts_g1(2,i)
       res(2) = res(2) + eta2_x2*pts_g1(2,i)
    end do
    res(1) = res(1)*edge_length
    res(2) = res(2)*edge_length
!!$    res(1) = res(1)*factor1*edge_length
!!$    res(2) = res(2)*factor1*edge_length
  end function normal_integral_eta2_plus
  
  ! integral of the normal vector over the 'southth' edge of the cell.
  function normal_integral_eta2_minus( T,ic,jc,integration_degree ) result(res)
    class(sll_coordinate_transformation_2d_base), pointer :: T
    sll_int32, intent(in)    :: ic
    sll_int32, intent(in)    :: jc
    sll_int32, intent(in)    :: integration_degree
    sll_real64, dimension(2) :: res
    
    sll_real64, dimension(2,2) :: inv_jac_mat ! inverse jacobian matrix
    sll_real64, dimension(2,integration_degree) :: pts_g1 ! gauss-legendre pts.
    !    sll_real64, dimension(2,integration_degree) :: pts_g2 ! gauss-legendre pts.
    sll_real64 :: eta2
    sll_real64 :: eta1_min
    !sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: max1
    sll_real64 :: min1
    !sll_real64 :: max2
    !sll_real64 :: min2
    !sll_real64 :: factor1
    !sll_real64 :: factor2
    sll_int32  :: i
    sll_real64 :: eta2_x1  ! derivative of eta2(x1,x2) with respect to x1
    sll_real64 :: eta2_x2  ! derivative of eta2(x1,x2) with respect to x2
    sll_real64 :: edge_length
    class(sll_cartesian_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_cartesian_mesh()

    ! verify that the indices requested are within the logical mesh.
    SLL_ASSERT(ic <= m%num_cells1)
    SLL_ASSERT(jc <= m%num_cells2)
    
    eta1_min = m%eta1_min
    delta1   = m%delta_eta1
    delta2   = m%delta_eta2
    
    ! The limits of integration are the limits of the cell in eta-space
    min1    = eta1_min + (ic-1)*delta1
    max1    = eta1_min + ic*delta1
    eta2    = (jc-1)*delta2 ! <- line differs w/ normal_integral_eta2_pluus()
    !    min2    = eta2_min + (jc-1)*delta2
    !    max2    = eta2_min + jc*delta2
    !factor1 = 0.5_f64*(max1-min1)
    !    factor2 = 0.5_f64*(max2-min2) 
    pts_g1(:,:) = gauss_legendre_points_and_weights(integration_degree, min1, max1)
    !gauss_points(integration_degree, min1, max1)
    !    pts_g2(:,:) = gauss_points(integration_degree, min2, max2)
    
    ! For efficiency, this code should be refactored. Consider:
    ! - adding direct access functions to the jacobian matrix elements and the
    !   inverse jacobian matrix elements
    ! - use macros to eliminate the massive code redundancy in this module's
    !   functions.
    ! - same macros can be used to improve a function call like the one next.
    edge_length = edge_length_eta2_minus( T, ic, jc, integration_degree )
    res(:) = 0.0_f64

    do i=1,integration_degree
       ! this can be made more efficient if we could access directly each
       ! term of the jacobian matrix independently.
       inv_jac_mat(:,:) = T%inverse_jacobian_matrix(pts_g1(1,i),eta2)
       eta2_x1 = inv_jac_mat(2,1)
       eta2_x2 = inv_jac_mat(2,2)
       SLL_ASSERT(T%jacobian(pts_g1(1,i),eta2) > 0.0_f64)
       res(1) = res(1) + eta2_x1*pts_g1(2,i)
       res(2) = res(2) + eta2_x2*pts_g1(2,i)
    end do
    ! change of sign due to direction of integration
    res(1) = -res(1)*edge_length
    res(2) = -res(2)*edge_length
!!$    res(1) = -res(1)*factor1*edge_length
!!$    res(2) = -res(2)*factor1*edge_length
  end function normal_integral_eta2_minus
  
  

!========================================================================


!> Compute some data about a triangular mesh 
!!                                                          
!!   nbs  - nombre de noeuds                          
!!   nbt  - nombre de triangles du maillage            
!!   coor - coordonnees des noeuds           
!!   refn - numeros de references des noeuds 
!!   ntri - numeros des sommets              
!!   vois - numeros des voisins des triangles
!!   voiv - numeros des voisins des triangles
!!   aire - aires des triangles              
!!   base - integrales des fonctions de base  
!!   nusd - numeros de sous-domaine            
!!
!!   npoel1 - pointeur du tableau "npoel2"       
!!                                             
!!   nucfl - num de cotes ne satisfaisant pas    
!!            la condition CFL                           
!!                                                      
!!   petitl - petite longueur de reference         
!!   grandl - grande longueur de reference       
!!   imxref - grand nombre entier de reference  
!!                                             
!!   nbtcot - nombre total de cotes           
!!   nbcoti - nombre de cotes internes       
!!   nbcfli - nbre de cotes internes ne satisfaisant pas CFL 
!!   ncotcu - nombre de cotes internes et frontieres        
!!   ncotcu - nombre de cotes cumules :                    
!!   - 1 = nombre de cotes internes effectifs             
!!   - 2 = 1 + nombre de cotes internes fictifs          
!!   - 3 = 2 + nombre de cotes frontieres references 1  
!!   - 4 = 3 + nombre de cotes frontieres references 2 , etc... 
!!
!!   nelmatf - nombre d'elements de la matrice profil Laplacien
!!                                                                     
!<
subroutine analyze_triangular_mesh(mesh, vmsh, ntypfr)

integer, dimension(:), intent(in) :: ntypfr
type(mesh_data), intent(inout)    :: mesh
type(voronoi), intent(out)        :: vmsh

integer, dimension (:), allocatable :: indc
integer, dimension (:), allocatable :: nuctfr

integer :: i, it, ntmp, id1, nct, iel, ind, nel
integer :: i1, i2, i3, is, is1, is2, is3
integer :: j, jel1, jel2, jel3, nel1, nel2, nel3
integer :: ic, ifr
integer :: ict, ictcl, iref, jref
integer :: k, keltmp, kcttmp, kretmp, ks1tmp, ks2tmp
integer :: nbcot, ict1, ict2, ictcl1, ictcl2, ie
integer :: iel1

real(8) :: lx1, lx2, ly1, ly2, x1, y1, x2, y2
real(8) :: xtmp
real(8) :: airtot

logical :: ldefr, ltemp, lcalsu, lmainv

integer :: itab(3,3)
integer, parameter :: imxref = 9999
integer, parameter :: iout   = 6

!-----
!---- itab(isomloc1,isomloc2) est le numero local de l'arete qui a pour
!---- sommets locaux isomloc1 et isomloc2
!-----

itab = reshape(source=(/0,1,3,1,0,2,3,2,0/), shape=(/3,3/))

ldefr  = .false.
ltemp  = .true.
lcalsu = .false.
lmainv = .false.

write(iout,"(/////10x,'>>> Calcul des quantites liees au maillage <<<'/)")

!*** Calcul des longueurs de reference
#ifdef DEBUG
write(iout,*)"*** Calcul des longueurs de reference ***"
#endif /* DEBUG */

xlml = minval(mesh%coor(1,:))
xlmu = maxval(mesh%coor(1,:))
ylml = minval(mesh%coor(2,:))
ylmu = maxval(mesh%coor(2,:))

mesh%petitl = 1.e-04 * min(xlmu-xlml,ylmu-ylml)/sqrt(float(mesh%nbs))
mesh%grandl = 1.e+04 * max(xlmu-xlml,ylmu-ylmu)

!*** Correction des erreurs de precision ----------------------
!    pour les noeuds sur l'axe en axisymetrique.
      
!*** Calcul des aires des triangles
#ifdef DEBUG
write(iout,*)"*** Calcul des aires des triangles ***"
#endif /* DEBUG */

allocate(mesh%aire(mesh%nbt)); mesh%aire=0.0

airtot = 0.

do it = 1, mesh%nbt

   lx1 = mesh%coor(1,mesh%ntri(2,it))-mesh%coor(1,mesh%ntri(1,it))
   ly1 = mesh%coor(2,mesh%ntri(3,it))-mesh%coor(2,mesh%ntri(1,it))
   lx2 = mesh%coor(1,mesh%ntri(3,it))-mesh%coor(1,mesh%ntri(1,it))
   ly2 = mesh%coor(2,mesh%ntri(2,it))-mesh%coor(2,mesh%ntri(1,it))

   mesh%aire(it) = 0.5 * abs(lx1*ly1 - lx2*ly2)

   if( mesh%aire(it) <= 0. ) then
     write(iout,*) " Triangle : ", it
     write(iout,*) mesh%ntri(1,it), ":",mesh%coor(1:2,mesh%ntri(1,it))
     write(iout,*) mesh%ntri(2,it), ":",mesh%coor(1:2,mesh%ntri(2,it))
     write(iout,*) mesh%ntri(3,it), ":",mesh%coor(1:2,mesh%ntri(3,it))
     SLL_ERROR("Aire de triangle negative")
   end if

   airtot = airtot + mesh%aire(it)

end do

!Calcul des integrales des fonctions de base
write(iout,*)"*** Calcul des integrales des fonctions de base ***"

allocate(mesh%xbas(mesh%nbs)); mesh%xbas=0.0

xtmp = 2. * sll_pi / 6.

do it = 1, mesh%nbt

   is1 = mesh%ntri(1,it) 
   is2 = mesh%ntri(2,it) 
   is3 = mesh%ntri(3,it) 

   mesh%xbas(is1) = mesh%xbas(is1) + mesh%aire(it)/3.
   mesh%xbas(is2) = mesh%xbas(is2) + mesh%aire(it)/3.
   mesh%xbas(is3) = mesh%xbas(is3) + mesh%aire(it)/3.

end do

write(iout,"(/10x,'Longueurs de reference :',2E15.5/)") mesh%petitl,mesh%grandl
write(iout,"(/10x,'Limites x du domaine   :',2E15.5/    &
  &        10x,'Limites y du domaine   :',2E15.5/   &
  &        10x,'Aire des triangles     :', E15.5/)") xlml,xlmu,ylml,ylmu,airtot


! --- Gestion des triangles ayant un noeud en commun -----------
write(iout,*)"*** Gestion des triangles ayant un noeud en commun ***"
 
! ... recherche des elements ayant un sommet commun
!     creation du tableau npoel1(i+1)  contenant le nombre de 
!     triangles ayant le noeud i en commun

allocate(mesh%npoel1(mesh%nbs+1))

mesh%npoel1 = 0
do i=1,mesh%nbt
   is1 = mesh%ntri(1,i)
   is2 = mesh%ntri(2,i)
   is3 = mesh%ntri(3,i)
   mesh%npoel1(is1+1) = mesh%npoel1(is1+1)+1
   mesh%npoel1(is2+1) = mesh%npoel1(is2+1)+1
   mesh%npoel1(is3+1) = mesh%npoel1(is3+1)+1
end do

! ... le tableau npoel1 devient le tableau donnant l'adresse 
!     dans npoel2 du dernier element dans la suite des triangles
!     communs a un noeud

mesh%npoel1(1)=0
do i=3,mesh%nbs+1
   mesh%npoel1(i)=mesh%npoel1(i-1)+mesh%npoel1(i)
end do

! ... creation du tableau npoel2 contenant sequentiellement les 
!     numeros des triangles ayant un noeud en commun      
!     le premier triangle s'appuyant sur le noeud i est
!     adresse par "npoel1(i)+1" 
!     le nombre de triangles ayant le noeud i en commun est
!     "npoel1(i+1)-npoel1(i)"


allocate(mesh%npoel2(mesh%npoel1(mesh%nbs+1)))
allocate(indc(mesh%nbs))

indc   = 1  !Le tableau temporaire indc doit etre initialise a 1

do it = 1,mesh%nbt
   do k = 1,3
      is = mesh%ntri(k,it)
      mesh%npoel2(mesh%npoel1(is)+indc(is)) = it
      indc(is) = indc(is)+1
   end do
end do
 
! --- Recherche des numeros des triangles voisins d'un triangle 

write(iout,*)"*** Recherche des numeros des triangles voisins d'un triangle ***"
do i = 1, mesh%nbt
   write(iout,*) " Triangle ", i, " Voisins :", mesh%nvois(1:3,i)
end do

do iel=1,mesh%nbt

   ! ... numeros des 3 sommets du triangle

   is1=mesh%ntri(1,iel)
   is2=mesh%ntri(2,iel)
   is3=mesh%ntri(3,iel)

   ! ... boucles imbriquees sur les elements pointant vers
   !     les 2 noeuds extremites de l'arete consideree
   !     Le voisin est le triangle commun (hormis iel)

   ! ... premiere arete (entre le sommet is1 et is2)

   nel1=mesh%npoel1(is1+1)-mesh%npoel1(is1) !nb de triangles communs a is1
   nel2=mesh%npoel1(is2+1)-mesh%npoel1(is2) !nb de triangles communs a is2

   loop1:do i1=1,nel1
            jel1=mesh%npoel2(mesh%npoel1(is1)+i1) !premier triangle is1
            if(jel1.ne.iel) then
               do i2=1,nel2
                  jel2=mesh%npoel2(mesh%npoel1(is2)+i2)
                  if(jel2 == jel1) then
                     mesh%nvois(1,iel)  = jel1
                     exit loop1
              end if
           end do
            end if
     end do loop1

   ! ... deuxieme arete (entre le sommet is2 et is3)

   nel2=mesh%npoel1(is2+1)-mesh%npoel1(is2)
   nel3=mesh%npoel1(is3+1)-mesh%npoel1(is3)

   loop2:do i2=1,nel2
            jel2=mesh%npoel2(mesh%npoel1(is2)+i2)
            if(jel2 /= iel) then
               do i3=1,nel3
                  jel3=mesh%npoel2(mesh%npoel1(is3)+i3)
                  if(jel3 == jel2) then
                     mesh%nvois(2,iel)=jel2
             exit loop2
              end if
           end do
            end if
     end do loop2

   ! ... troisieme arete (entre le sommet is3 et is1)

   nel3=mesh%npoel1(is3+1)-mesh%npoel1(is3)
   nel1=mesh%npoel1(is1+1)-mesh%npoel1(is1)

   loop3:do i3=1,nel3
            jel3=mesh%npoel2(mesh%npoel1(is3)+i3)
            if(jel3 /= iel) then
               do i1=1,nel1
                  jel1=mesh%npoel2(mesh%npoel1(is1)+i1)
                  if(jel1 == jel3) then
                     mesh%nvois(3,iel)=jel3
             exit loop3
              end if
           end do
            end if
     end do loop3
end do

!Calcul de nvoif : Numero local dans le voisin de la face j commune a l'elt i

allocate(mesh%nvoif(3,mesh%nbt))
mesh%nvoif(1:3,:) = mesh%nvois(1:3,:)

do i = 1, mesh%nbt
   !write(*,"(/,4i7)") i, mesh%ntri(1:3,i)
   do j = 1,3
      jel1 = mesh%nvois(j,i)
      if (jel1 > 0) then
         do k = 1,3
            jel2 = mesh%nvois(k,jel1)
            if (jel2 == i) then 
               mesh%nvoif(j,i) = k
               exit
            end if
         end do
         !write(*,"(5i7)") jel1, mesh%ntri(1:3,jel1), mesh%nvoif(j,i)
      end if
   end do
end do

! --- Definition de nctfrt: le nombre de cotes frontieres

mesh%nctfrt=0
do i=1,mesh%nbt
   if (mesh%nvois(1,i) < 0) mesh%nctfrt=mesh%nctfrt+1
   if (mesh%nvois(2,i) < 0) mesh%nctfrt=mesh%nctfrt+1
   if (mesh%nvois(3,i) < 0) mesh%nctfrt=mesh%nctfrt+1
end do

! --- Rangement de npoel2 dans l'ordre trigonometrique ---------

do is=1,mesh%nbs

   nel =mesh%npoel1(is+1)-mesh%npoel1(is)

   if ( nel > 1 ) then

      !*** Noeuds internes (Numero de reference nul) ***

      if( mesh%refs(is) == 0) then

        ind =1
        iel1=mesh%npoel2(mesh%npoel1(is)+1)

   loop4:do iel=2,nel-1
            do j=1,3
               if(mesh%ntri(j,iel1) == is) nct=mod(j+1,3)+1
            end do

            iel1=mesh%nvois(nct,iel1)
            do id1=ind+1,nel
               if(iel1 == mesh%npoel2(mesh%npoel1(is)+id1)) then
                  ind=ind+1
                  ntmp=mesh%npoel2(mesh%npoel1(is)+ind)
                  mesh%npoel2(mesh%npoel1(is)+ind)=iel1
                  mesh%npoel2(mesh%npoel1(is)+id1)=ntmp
              cycle loop4
               end if
            end do
         end do loop4

      ! Noeuds frontieres

      else 

      ! --> Recherche du premier triangle dans l'ordre trigonometrique

   loop5:do id1=1,nel
            iel1=mesh%npoel2(mesh%npoel1(is)+id1)
            do j=1,3
               if(mesh%nvois(j,iel1).le.0 .and. mesh%ntri(j,iel1) == is) then
                  ntmp=mesh%npoel2(mesh%npoel1(is)+1)
                  mesh%npoel2(mesh%npoel1(is)+1)=iel1
                  mesh%npoel2(mesh%npoel1(is)+id1)=ntmp
              exit loop5
               end if
            end do
         end do loop5
            
      ! --> Rangement des autres triangles dans l'ordre trigonometrique
      !     (s'il y en a plus que 2) 

         if(nel  > 2) then
            ind =1
            iel1=mesh%npoel2(mesh%npoel1(is)+1)
   
      loop6:do iel=2,nel-1
               do j=1,3
              if(mesh%ntri(j,iel1)==is) then
                 nct=mod(j+1,3)+1
              end if
               end do

               iel1=mesh%nvois(nct,iel1)
   
               do id1=ind+1,nel
                  if(iel1 == mesh%npoel2(mesh%npoel1(is)+id1)) then
                     ind=ind+1
                     ntmp=mesh%npoel2(mesh%npoel1(is)+ind)
                     mesh%npoel2(mesh%npoel1(is)+ind)=iel1
                     mesh%npoel2(mesh%npoel1(is)+id1)=ntmp
                 cycle loop6
              end if
               end do

            end do loop6

         end if

     end if

  end if

end do

do iel = 1, mesh%nbt
   do j = 1, 3
   if( mesh%nvois(j,iel) == imxref ) then
      write(iout,*) " Triangle ", iel, " Voisins :", (mesh%nvois(i,iel),i=1,3)
      write(iout,*) " Coordonnees x =", mesh%coor(1,mesh%ntri(1,iel))
      write(iout,*) " Coordonnees y =", mesh%coor(2,mesh%ntri(1,iel))
      stop
   end if
   end do
end do

!Calcul des voisins pour les particules
write(iout,*)"*** Calcul des voisins pour les particules ***"

! --- Remplissage de "nvoiv" -----------------------------------
!     Identique a "nvois" sauf pour les aretes appartenant a 
!     une frontiere. Le chiffre correspond ici a un code pour 
!     le traitement des conditions aux limites sur les 
!     particules, alors que dans "nvois" ce chiffre est 
!     l'oppose du numero de reference de la frontiere concernee

allocate(mesh%nvoiv(3,mesh%nbt))

do i = 1,mesh%nbt

   ! ... Cotes internes
   do j = 1, 3
      if (mesh%nvois(j,i)>0) then
         mesh%nvoiv(j,i) = mesh%nvois(j,i)
      end if
   end do

   ! ... Cotes frontieres
   do j = 1, 3
      if (mesh%nvois(j,i)<0) then

         select case (ntypfr(-mesh%nvois(j,i)))

         case(0)
            mesh%nvoiv(j,i) = -1
         case(1)
            mesh%nvoiv(j,i) =  0
         case(2)
            mesh%nvoiv(j,i) = -1
         case(3)
            mesh%nvoiv(j,i) =  0
         case(4)
            mesh%nvoiv(j,i) =  0
         case(5)
            mesh%nvoiv(j,i) = -1
         case(6)
            mesh%nvoiv(j,i) = -1
         case(7) 
            mesh%nvoiv(j,i) =  0
         end select

      end if
   end do    

end do

!  nbtcot : nombre total de cotes                           *
!  nbcoti : nombre de cotes internes                        *
!  ncotcu : nombre de cotes internes et frontieres          *
!  ncotcu : nombre de cotes cumules :                       *
!  1 = nombre de cotes internes effectifs                   *
!  2 = 1 + nombre de cotes internes fictifs                 *
!  3 = 2 + nombre de cotes frontieres references 1          *
!  4 = 3 + nombre de cotes frontieres references 2 , etc... *

write(iout,"(/10x,a,i3)") 'Nombre maximum de frontieres referencees ', mesh%nmxfr

!*** Calcul du nb de cotes internes et total (nbcoti,nbtcot)

mesh%nbtcot = (3*mesh%nbt+mesh%nctfrt)/2

write(iout,"( 10x,a,i6)") 'Nombre total de cotes =', mesh%nbtcot

mesh%nbcoti = mesh%nbtcot - mesh%nctfrt

write(iout,"( 10x,a,i6)") 'Nombre de cotes internes =', mesh%nbcoti
write(iout,"( 10x,a,i6/)")'Nombre de cotes frontieres =', mesh%nctfrt
 
!  Calcul du nb de cotes cumules par type de traitement (ncotcu) ....

allocate(vmsh%ncotcu(mesh%nbcoti+mesh%nmxfr*mesh%nctfrt))

vmsh%ncotcu=0
vmsh%ncotcu(1)=mesh%nbcoti
vmsh%ncotcu(2)=mesh%nbcoti

do ifr=1,mesh%nmxfr
   do i=1,mesh%nbt
      do j = 1,3
         if(mesh%nvois(j,i) == -ifr) then
            vmsh%ncotcu(ifr+2) = vmsh%ncotcu(ifr+2) + 1
         end if
      end do
   end do
end do

do ifr=1,mesh%nmxfr
   vmsh%ncotcu(ifr+2)=vmsh%ncotcu(ifr+1)+vmsh%ncotcu(ifr+2)
end do

write(iout,"(10x,'Nombre de cotes cumules a partir de nbcoti :')")
do ifr=1,mesh%nmxfr
   write(iout,"(20x,'ref = ',i3,'   ncotcu(i) = ',i6)") ifr,vmsh%ncotcu(ifr+2)
end do

allocate(vtaux(mesh%nbtcot),vtauy(mesh%nbtcot))
vtaux = 0.; vtauy = 0.

!==============================================================!
!--- Calcul complementaires relatifs au lissage des champs ----!
!    dans le cas de la resolution de POISSON               !
!==============================================================!

!tableau des numeros des noeuds associes a un cote 
allocate(vmsh%nuvac(2,mesh%nbtcot)); vmsh%nuvac = 0

!tableau des longueurs des cotes 
allocate(vmsh%xlcod(mesh%nbtcot)); vmsh%xlcod = 0.0

!tableaux de lissage et des cotes tangeants 
allocate(mesh%xmal1(mesh%nbs),mesh%xmal2(mesh%nbs),mesh%xmal3(mesh%nbs))
mesh%xmal1=0.;mesh%xmal2=0.;mesh%xmal3=0.

allocate(vmsh%nbcov(mesh%nbs+1))  !pointeur des cotes pointant sur le meme noeud
allocate(vmsh%nugcv(10*mesh%nbs)) !tableau contenant les numeros de ces cotes
allocate(nuctfr(mesh%nmxfr)) !tableau temporaire                              

!Creation du maillage de Voronoi et quantites associees

call poclis(mesh, vmsh, nuctfr, vtaux, vtauy)

!On tue le tableau temporaire         I                      

deallocate(nuctfr) 


!Event: Numerotation des cotes frontieres
!
!    Numerotation des cotes frontieres dans le sens 
!    trigonometrique, par numero de reference, pour
!    effectuer tous les diagnostics relatifs aux  
!    frontieres.                                 
!                                               
!    kelfro - numero d'element auquel appartient le cote 
!    kctfro - numero local (1,2,3) de cote              
!    krefro - numero de reference du cote              
!    ksofro - numeros des 2 sommets extremite du cote 
!    vnofro - composantes du vecteur normal (vers l'interieur)
!                                                            
!    nctfrt - Nb total de cotes frontieres                  
!    nctfro - Nb de cotes frontieres par reference         
!    nctfrp - Pointeur de tableau (nb de cote par ref)    
!                                                        

! --- Initialisation (numerotation quelconque) -----------------

allocate(nctfrp(0:mesh%nmxfr))
allocate(nctfro(mesh%nmxfr))
allocate(mesh%kctfro(mesh%nctfrt))
allocate(mesh%kelfro(mesh%nctfrt))
allocate(mesh%krefro(mesh%nctfrt))
allocate(mesh%ksofro(2,mesh%nctfrt))
allocate(mesh%vnofro(2,mesh%nctfrt))

nctfro = 0

!*** Premiere numerotation des cotes frontieres ***

ifr=0
do ict=1,3
   do iel=1,mesh%nbt
      if ( mesh%nvois(ict,iel) < 0 ) then 

         iref = -mesh%nvois(ict,iel)

         ifr  = ifr+1
         mesh%kelfro(ifr) = iel  !element auquel appartient le cote
         mesh%kctfro(ifr) = ict  !numero local du cote dans cet element
         mesh%krefro(ifr) = iref !reference de ce cote
         mesh%ksofro(1,ifr) = mesh%ntri(ict,iel) !numeros des 2 sommets de ce cote
         mesh%ksofro(2,ifr) = mesh%ntri(mod(ict,3)+1,iel) 

         nctfro(iref)  = nctfro(iref)+1  !nombre de cote par reference

      end if 
   end do
end do

!... Pointeur de tableau .........................................

nctfrp(0)=0
do ifr=1,mesh%nmxfr
   nctfrp(ifr)=nctfrp(ifr-1)+nctfro(ifr)
end do

!--- Balayage sur les numeros de reference --------------------
!    Premier classement par numero de reference        

ictcl=1
do iref=1,mesh%nmxfr

   nbcot = nctfro(iref)
   if (nbcot > 0) then 

      ict1 = ictcl

      do ict=ict1,mesh%nctfrt
         jref = mesh%krefro(ict)

         if (jref == iref) then 

            keltmp = mesh%kelfro(ict)
            kcttmp = mesh%kctfro(ict)
            kretmp = mesh%krefro(ict)
            ks1tmp = mesh%ksofro(1,ict)
            ks2tmp = mesh%ksofro(2,ict)

            mesh%kelfro(ict) = mesh%kelfro(ictcl)
            mesh%kctfro(ict) = mesh%kctfro(ictcl)
            mesh%krefro(ict) = mesh%krefro(ictcl)
            mesh%ksofro(1,ict) = mesh%ksofro(1,ictcl)
            mesh%ksofro(2,ict) = mesh%ksofro(2,ictcl)

            mesh%kelfro(ictcl) = keltmp
            mesh%kctfro(ictcl) = kcttmp
            mesh%krefro(ictcl) = kretmp
            mesh%ksofro(1,ictcl) = ks1tmp
            mesh%ksofro(2,ictcl) = ks2tmp

            ictcl=ictcl+1

         end if 
      end do
   end if 
end do

! --- Rangement dans l'ordre trigonometrique -------------------

do iref=1,mesh%nmxfr

   nbcot=nctfro(iref)

   !Une seule reference pour toute la frontiere ................

   if(nbcot  ==  mesh%nctfrt) then 
      ict1=1
 35   continue
      is2=mesh%ksofro(2,ict1)
      do ict2=ict1,mesh%nctfrt
         if (ict1 /= ict2) then 
            is1=mesh%ksofro(1,ict2)
            if (is1==is2) then 
            
               keltmp = mesh%kelfro(ict1+1)
               kcttmp = mesh%kctfro(ict1+1)
               kretmp = mesh%krefro(ict1+1)
               ks1tmp = mesh%ksofro(1,ict1+1)
               ks2tmp = mesh%ksofro(2,ict1+1)

               mesh%kelfro(ict1+1) = mesh%kelfro(ict2)
               mesh%kctfro(ict1+1) = mesh%kctfro(ict2)
               mesh%krefro(ict1+1) = mesh%krefro(ict2)
               mesh%ksofro(1,ict1+1) = mesh%ksofro(1,ict2)
               mesh%ksofro(2,ict1+1) = mesh%ksofro(2,ict2)

               mesh%kelfro(ict2) = keltmp
               mesh%kctfro(ict2) = kcttmp
               mesh%krefro(ict2) = kretmp
               mesh%ksofro(1,ict2) = ks1tmp
               mesh%ksofro(2,ict2) = ks2tmp

               if (ict1<mesh%nctfrt) then 
              ict1=ict1+1
          GOTO 35
           end if

       end if
        end if
     end do

     ! ... Plusieurs references sur la frontiere ...........................

  else if (nbcot>1) then 

     ictcl1=nctfrp(iref-1)+1
     ictcl2=nctfrp(iref)

     ! ... Premier cote d'un segment ........................................

     ict1=ictcl1
 31  continue
     is1=mesh%ksofro(1,ict1)
     do ict2=ictcl1,ictcl2
        if (ict1 /= ict2) then 
           is2=mesh%ksofro(2,ict2)
           if (is1==is2) then 
              ict1=ict2
              if (ict1==ictcl1) then !Test si on tourne en rond ............
                 GOTO 37
              end if
              GOTO 31
           end if
        end if
     end do

! ... Permutation ......................................................
               
     keltmp=mesh%kelfro(ict1)
     kcttmp=mesh%kctfro(ict1)
     kretmp=mesh%krefro(ict1)
     ks1tmp=mesh%ksofro(1,ict1)
     ks2tmp=mesh%ksofro(2,ict1)

     mesh%kelfro(ict1)=mesh%kelfro(ictcl1)
     mesh%kctfro(ict1)=mesh%kctfro(ictcl1)
     mesh%krefro(ict1)=mesh%krefro(ictcl1)

     mesh%ksofro(1,ict1)=mesh%ksofro(1,ictcl1)
     mesh%ksofro(2,ict1)=mesh%ksofro(2,ictcl1)

     mesh%kelfro(ictcl1)=keltmp
     mesh%kctfro(ictcl1)=kcttmp
     mesh%krefro(ictcl1)=kretmp

     mesh%ksofro(1,ictcl1)=ks1tmp
     mesh%ksofro(2,ictcl1)=ks2tmp
        
! ... Classement des cotes suivants ....................................

 37  continue
     ict1=ictcl1
 33  continue
     is2=mesh%ksofro(2,ict1)
     do ict2=ict1,ictcl2
        if (ict1 /= ict2) then 
           is1=mesh%ksofro(1,ict2)
           if (is1==is2) then 
        
              keltmp=mesh%kelfro(ict1+1)
              kcttmp=mesh%kctfro(ict1+1)
              kretmp=mesh%krefro(ict1+1)
              ks1tmp=mesh%ksofro(1,ict1+1)
              ks2tmp=mesh%ksofro(2,ict1+1)

              mesh%kelfro(ict1+1)=mesh%kelfro(ict2)
              mesh%kctfro(ict1+1)=mesh%kctfro(ict2)
              mesh%krefro(ict1+1)=mesh%krefro(ict2)
              mesh%ksofro(1,ict1+1)=mesh%ksofro(1,ict2)
              mesh%ksofro(2,ict1+1)=mesh%ksofro(2,ict2)

              mesh%kelfro(ict2)=keltmp
              mesh%kctfro(ict2)=kcttmp
              mesh%krefro(ict2)=kretmp
              mesh%ksofro(1,ict2)=ks1tmp
              mesh%ksofro(2,ict2)=ks2tmp

              ict1=ict1+1
              GOTO 33
           end if
        end if
     end do

! ...Test s'il y a d'autres segments de la meme reference .............

     if (ict1<ictcl2) then 
        ictcl1=ict1+1
        ict1=ictcl1
        GOTO 31
     end if
   end if
end do

!  Modification de nvois ------------------------------------
! (Numero de reference change par le numero de cote)
 
do ict=1,mesh%nctfrt
   ie=mesh%kelfro(ict)
   ic=mesh%kctfro(ict)
   mesh%nvois(ic,ie)=-ict
end do

!--- Calcul des composantes du vecteur normal -----------------

do ict=1,mesh%nctfrt

   is1 = mesh%ksofro(1,ict)
   is2 = mesh%ksofro(2,ict)

   x1 = mesh%coor(1,is1)
   y1 = mesh%coor(2,is1)
   x2 = mesh%coor(1,is2)
   y2 = mesh%coor(2,is2)

   mesh%vnofro(1,ict) = -y2+y1
   mesh%vnofro(2,ict) =  x2-x1

end do

deallocate(indc)

end subroutine analyze_triangular_mesh

!**************************************************************


!Subroutine: poclis
!                                               
!   Calcul des matrices de lissage associees a chaque     
!   noeud du maillage.                                  
!   Necessaires au calcul des composantes de E1,E2     
!   a partir des composantes tangeantielles de E      
!   connues sur les cotes des triangles.             
!   Calcul des composantes des vecteurs unitaires tangeants.                                     
!                                                              
!   Variables d'entree:                                    
!  
!      nuvac  - numeros des PV associes aux cotes          
!      coor   - coordonnees des noeuds Delaunay           
!      xlcod  - longueur des cotes Delaunay              
!      npoel1 - pointeur du tableau npoel2              
!      npoel2 - numeros des triangles entourant un noeud       
!      xmal1  - somme des taux*taux entourant un noeud (/det) 
!      xmal2  - somme des tauy*tauy entourant un noeud (/det)
!      xmal3  - somme des taux*tauy entourant un noeud (/det)
!      vtaux  - composante x des vecteurs tangeants         
!      vtauy  - composante y des vecteurs tangeants        
!                                                               
!Auteur:
! A. Adolf - Version 1.0   Septembre 1994  
subroutine poclis(mesh, vmsh, nuctfr, vtaux, vtauy)

type(mesh_data) :: mesh
type(voronoi)   :: vmsh
double precision, dimension(:) :: vtaux, vtauy
integer, dimension(:) :: nuctfr
logical :: lerr
double precision :: det, s1, s2, s3, x21, y21, xa, ya, xb, yb
integer :: ic, nuctf, ind, nm1, nm2, nel1, nel2, indv1, indv2
integer :: indn1, indn2, num1, num2, n1, n2, ivois, nc, iel, nucti
integer :: nbti, is, i
integer, parameter :: iout = 6
     
!======================================================================
! --- 1.0 --- Pointeur des numeros de cotes pointant vers un noeud -----

do is=1,mesh%nbs
   nbti=mesh%npoel1(is+1)-mesh%npoel1(is)
   vmsh%nbcov(is+1)=nbti
   if(mesh%refs(is) /= 0) then
      vmsh%nbcov(is+1)=nbti+1
   end if
end do

vmsh%nbcov(1)=0
do is=1,mesh%nbs
   vmsh%nbcov(is+1)=vmsh%nbcov(is)+vmsh%nbcov(is+1)
end do

! --- 1.5 --- Tableau temporaire (cumul des cotes frontieres) ----------
nucti=0
do i=1,mesh%nmxfr
   nuctfr(i)=0
end do

! ======================================================================
! --- 2.0 --- Numerotation des cotes -----------------------------------

do iel=1,mesh%nbt
   do nc=1,3
      ivois=mesh%nvois(nc,iel)
      n1=nc
      n2=mod(nc,3)+1

      num1 =mesh%ntri(n1,iel)
      num2 =mesh%ntri(n2,iel)
      indn1=mesh%npoel1(num1)
      indn2=mesh%npoel1(num2)
      indv1=vmsh%nbcov(num1)
      indv2=vmsh%nbcov(num2)
      nel1 =mesh%npoel1(num1+1)-mesh%npoel1(num1)
      nel2 =mesh%npoel1(num2+1)-mesh%npoel1(num2)

      !Cas des cotes internes ...........................................

      if(ivois >  iel) then

         nucti=nucti+1

         !Numeros globaux de cotes pointant vers le meme noeud 

         do nm1=1,nel1
            if(mesh%npoel2(indn1+nm1) == ivois) then
               vmsh%nugcv(indv1+nm1)=nucti
            end if
         end do

         do nm2=1,nel2
            if(mesh%npoel2(indn2+nm2) == iel) then
               vmsh%nugcv(indv2+nm2)=nucti
            end if
         end do

         !Numeros des triangles ou polygones associes

         vmsh%nuvac(1,nucti)=num1
         vmsh%nuvac(2,nucti)=num2

         !Cas des cotes frontaliers ........................................

      else if(ivois < 0) then

         ind=-ivois
         nuctfr(ind)=nuctfr(ind)+1
         nuctf=vmsh%ncotcu(ind+1)+nuctfr(ind)

         vmsh%nugcv(indv1+nel1+1)=nuctf
         vmsh%nugcv(indv2+nel2  )=nuctf
         vmsh%nuvac(1,nuctf)=num1
         vmsh%nuvac(2,nuctf)=num2

      end if

   end do

end do

!======================================================================
!----------- Longueurs des cotes des triangles ------------------------

do ic=1,mesh%nbtcot
   xa=mesh%coor(1,vmsh%nuvac(1,ic))
   ya=mesh%coor(2,vmsh%nuvac(1,ic))
   xb=mesh%coor(1,vmsh%nuvac(2,ic))
   yb=mesh%coor(2,vmsh%nuvac(2,ic))
   vmsh%xlcod(ic)=sqrt((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb))
end do

!======================================================================
!--- 4.0 --- Calcul des matrices de lissage ---------------------------

do ic=1,mesh%nbtcot
   n1 = vmsh%nuvac(1,ic)
   n2 = vmsh%nuvac(2,ic)
   x21= (mesh%coor(1,n2)-mesh%coor(1,n1))/vmsh%xlcod(ic)
   y21= (mesh%coor(2,n2)-mesh%coor(2,n1))/vmsh%xlcod(ic)
   s1 = x21*x21
   s2 = y21*y21
   s3 = x21*y21
   vtaux(ic)=x21
   vtauy(ic)=y21
   mesh%xmal1(n1)=mesh%xmal1(n1)+s1
   mesh%xmal2(n1)=mesh%xmal2(n1)+s2
   mesh%xmal3(n1)=mesh%xmal3(n1)+s3
   mesh%xmal1(n2)=mesh%xmal1(n2)+s1
   mesh%xmal2(n2)=mesh%xmal2(n2)+s2
   mesh%xmal3(n2)=mesh%xmal3(n2)+s3
end do

!--- 4.5 --- Normalisation par rapport aux determinants ---------------
lerr = .false.

do is=1,mesh%nbs
   det=mesh%xmal1(is)*mesh%xmal2(is)-mesh%xmal3(is)*mesh%xmal3(is)
   if(det /= 0) then
      mesh%xmal1(is)=mesh%xmal1(is)/det
      mesh%xmal2(is)=mesh%xmal2(is)/det
      mesh%xmal3(is)=mesh%xmal3(is)/det
   else
      lerr=.TRUE.
   end if
end do
 
if(lerr) then
   write(iout,900)
   SLL_ERROR("poclis")
end if

! ======================================================================
! --- 5.0 --- Impression des tableaux ----------------------------------
 
!write(iout,902)
!do is=1,mesh%nbs
!   write(iout,903) mesh%xmal1(is),mesh%xmal2(is),mesh%xmal3(is)
!end do

!write(iout,904)
!do ic=1,mesh%nbtcot
!   write(iout,905) vtaux(ic),vtauy(ic)
!end do

if(lerr) then
   write(iout,901)
   SLL_ERROR("poclis")
end if

 900  format(//10x,'Determinant des coefficients des matrices'  &
              /10x,'de lissage nul'             &
              /10x,'Il faut modifier le maillage'//)
 901  format(//10x,'Le nombre de triangles communs a un noeud'  &
              /10x,'est superieur a 12'             &   
              /10x,'Modifiez legerement le maillage'//)

end subroutine poclis

end module sll_mesh_calculus_2d_module
