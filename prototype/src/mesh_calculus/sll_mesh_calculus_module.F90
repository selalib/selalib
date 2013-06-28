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
  ! quantities like edge-lengths, areas and volumes of cell elements in 
  ! a deformed mesh.
  !
  ! --------------------------------------------------------------------------

  ! Try first a module without a native type since all information is 
  ! contained in a coordinate transformation... does this imply that these
  ! should be methods of the coordinate transformation type? Probably, but
  ! not necessarily if this module behaves exclusively as a client of the
  ! coordinate transformation module.

!!$  type :: sll_mesh_calculus_2d
!!$      class(sll_coordinate_transformation_2d_base), pointer :: T
!!$   contains
!!$     procedure, pass(obj) :: cell_volume => cell_vol
!!$  end type sll_mesh_calculus_2d


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
    sll_real64 :: factor1
    sll_real64 :: factor2
    sll_int32  :: i
    sll_int32  :: j

    ! Verify arguments
    SLL_ASSERT(associated(T))
    ! verify that the indices requested are within the logical mesh.
    SLL_ASSERT(ic <= T%mesh%num_cells1)
    SLL_ASSERT(jc <= T%mesh%num_cells2)

    vol = 0.0_f64

    eta1_min = T%mesh%eta1_min
    delta1   = T%mesh%delta_eta1
    eta2_min = T%mesh%eta2_min
    delta2   = T%mesh%delta_eta2

    ! This function carries out the integral of the jacobian evaluated on
    ! gauss-legendre points within the cell.
    ! The limits of integration are the limits of the cell in eta-space
    min1    = eta1_min + (ic-1)*delta1
    max1    = eta1_min + ic*delta1
    min2    = eta2_min + (jc-1)*delta2
    max2    = eta2_min + jc*delta2
    factor1 = 0.5_f64*(max1-min1)
    factor2 = 0.5_f64*(max2-min2) 
    pts_g1(:,:) = gauss_points(integration_degree, min1, max1)
    pts_g2(:,:) = gauss_points(integration_degree, min2, max2)

    do j=1,integration_degree
       do i=1,integration_degree
          jac = T%jacobian(pts_g1(1,i),pts_g2(1,j))
          vol = vol + abs(jac)*pts_g1(2,i)*pts_g2(2,j)
       end do
    end do
    vol = vol*factor1*factor2

  end function cell_volume

end module sll_mesh_calculus_2d_module
