module sll_mesh_calculus_2d_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_logical_meshes
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
    type(sll_logical_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_logical_mesh()
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
    type(sll_logical_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_logical_mesh()
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
    type(sll_logical_mesh_2d), pointer :: m
    
    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_logical_mesh()

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
    type(sll_logical_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_logical_mesh()

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
    type(sll_logical_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_logical_mesh()

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
    type(sll_logical_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_logical_mesh()

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
    type(sll_logical_mesh_2d), pointer :: m
    
    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_logical_mesh()

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
    type(sll_logical_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_logical_mesh()

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
    type(sll_logical_mesh_2d), pointer :: m

    ! Verify arguments
    SLL_ASSERT(associated(T))
    m => T%get_logical_mesh()

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
  
  
end module sll_mesh_calculus_2d_module
