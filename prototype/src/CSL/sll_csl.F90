module sll_csl
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"

  use numeric_constants
  use sll_splines
  use advection_field
  use distribution_function
  implicit none

type csl_workspace
   type (sll_spline_1D), pointer :: spl_eta1
   type (sll_spline_1D), pointer :: spl_eta2
end type csl_workspace
contains
  function new_csl_workspace(dist_func_2D)
    type (csl_workspace), pointer :: new_csl_workspace
    type (sll_distribution_function_2D_t), pointer  :: dist_func_2D 
    sll_int32  :: ierr
    sll_int32  :: nc_eta1
    sll_int32  :: nc_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max

    ! allocate pointer
    SLL_ALLOCATE(new_csl_workspace,ierr)

    ! get dimensions
    nc_eta1    = get_df_nc_eta1( dist_func_2D ) 
    eta1_min   = get_df_eta1_min( dist_func_2D )
    eta1_max   = get_df_eta1_max( dist_func_2D )
    nc_eta2    = get_df_nc_eta2( dist_func_2D ) 
    eta2_min   = get_df_eta2_min( dist_func_2D )
    eta2_max   = get_df_eta2_max( dist_func_2D )

    ! initialize splines
    new_csl_workspace%spl_eta1 => new_spline_1D( nc_eta1, eta1_min, eta1_max, PERIODIC_SPLINE )
    new_csl_workspace%spl_eta2 => new_spline_1D( nc_eta2, eta2_min, eta2_max, HERMITE_SPLINE )    
  end function new_csl_workspace

! subroutine csl_first_order
! Advances the distribution function on a time step deltat using a first order conservative semi-Lagrangian scheme
  subroutine csl_first_order(csl_work, dist_func_2D, advfield, deltat)
    type (csl_workspace), pointer :: csl_work
    type (sll_distribution_function_2D_t), pointer  :: dist_func_2D  ! distribution function to be advanced
    type (field_2D_vec2), pointer  :: advfield   ! advection_field
    sll_real64  ::  deltat  ! time step

    sll_real64, dimension(:), allocatable  ::  advfield_1D_1
    sll_real64, dimension(:), allocatable  ::  advfield_1D_2
    sll_real64, dimension(:), allocatable  ::  primitive1
    sll_real64, dimension(:), allocatable  ::  primitive2
    sll_real64, dimension(:), allocatable  ::  eta1_out 
    sll_real64, dimension(:), allocatable  ::  eta2_out
    sll_real64, dimension(:), allocatable  ::  zeros1
    sll_real64, dimension(:), allocatable  ::  zeros2
    sll_int32  :: i1
    sll_int32  :: i2
    sll_int32  :: ierr
    sll_int32  :: nc_eta1
    sll_int32  :: nc_eta2
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_max
    sll_real64 :: val
    sll_real64 :: avg
    sll_real64 :: eta1

    ! parameter checking
    SLL_ASSERT(associated(csl_work))
    SLL_ASSERT(associated(dist_func_2D))
    SLL_ASSERT(associated(advfield))

    ! get dimensions
    nc_eta1    = get_df_nc_eta1( dist_func_2D ) 
    delta_eta1 = get_df_delta_eta1( dist_func_2D )
    eta1_min   = get_df_eta1_min( dist_func_2D )
    eta1_max   = get_df_eta1_max( dist_func_2D )
    nc_eta2    = get_df_nc_eta2( dist_func_2D ) 
    delta_eta2 = get_df_delta_eta2( dist_func_2D )
    eta2_min   = get_df_eta2_min( dist_func_2D )
    eta2_max   = get_df_eta2_max( dist_func_2D )
    
    ! allocation
    SLL_ALLOCATE(zeros1(nc_eta1+1),ierr)
    SLL_ALLOCATE(zeros2(nc_eta2+1),ierr)
    SLL_ALLOCATE(advfield_1D_1(nc_eta1+1),ierr)
    SLL_ALLOCATE(advfield_1D_2(nc_eta2+1),ierr)
    SLL_ALLOCATE(primitive1(nc_eta1+1),ierr)
    SLL_ALLOCATE(primitive2(nc_eta2+1),ierr)
    SLL_ALLOCATE(eta1_out(nc_eta1+1),ierr)
    SLL_ALLOCATE(eta2_out(nc_eta2+1),ierr)

    zeros1(:) = 0.0_f64
    zeros2(:) = 0.0_f64
    !print*, 'zeros2', zeros2
    ! advection along the first direction (assumed periodic)
    do i2=1, nc_eta2
       primitive1 (1) = 0.0_f64  ! set primitive to 0 on left boundary 
       advfield_1D_1 ( 1 ) = FIELD_2D_AT_I_V1( advfield, 1, i2 )
       do i1 = 2, nc_eta1+1
          ! extract subarray from advection field
          advfield_1D_1 ( i1 ) = FIELD_2D_AT_I_V1( advfield, i1, i2 )
          ! compute primitive of distribution function along this line
          primitive1 ( i1 ) = primitive1 ( i1-1 ) &
               + delta_eta1 * sll_get_df_val( dist_func_2D, i1-1, i2 )
       end do
       avg = primitive1 ( nc_eta1+1 ) / (eta1_max - eta1_min) ! average of dist func along the line
       ! modify primitive so that it becomes periodic
       eta1 = eta1_min
       do i1 = 2, nc_eta1+1
          eta1 = eta1 + delta_eta1
          primitive1 ( i1 ) = primitive1 ( i1 ) - avg * (eta1-eta1_min)
       end do
       call compute_flow_1D_backward( advfield_1D_1, zeros1, 1.0_f64, deltat, &
                                      eta1_min, nc_eta1, delta_eta1,          &
                                      PERIODIC, eta1_out ) 
       do i1 = 1, nc_eta1 
      !    print*, 'sll_csl before splines', i1,i2, primitive1 ( i1 )
       end do
       call compute_spline_1D_periodic( primitive1, nc_eta1, csl_work%spl_eta1 )
       ! interpolate primitive at origin of characteritics
       do i1 = 1, nc_eta1 
       !   print*, 'sll_c', i1,i2, eta1_out(i1), primitive1 ( i1 ), csl_work%spl_eta1%c(i1)
       end do
       call interpolate_array_values( eta1_out, primitive1, nc_eta1+1, csl_work%spl_eta1 )
       !do i1 = 1, nc_eta1 + 1         
       !   primitive1 ( i1 ) = interpolate_value( eta1_out(i1) , csl_work%spl_eta1 )
       !end do
       ! update average value of distribution function in cell using difference of primitives
       do i1 = 1, nc_eta1 
       !   print*, 'sll_csl', i1,i2, eta1_out(i1), primitive1 ( i1 ), csl_work%spl_eta1%c(i1) 
          val = (primitive1 ( i1+1 ) - primitive1 ( i1 )) / delta_eta1 + avg
          call sll_set_df_val( dist_func_2D, i1, i2, val )
          !print*, 'sll_csl', i1,i2, val
       end do
    end do
    ! advection along the second direction
    do i1=1, nc_eta1
       primitive2 (1) = 0.0_f64  ! set primitive to 0 on left boundary 
       advfield_1D_2(1) = FIELD_2D_AT_I_V1( advfield, i1, 1 )
       do i2 = 2, nc_eta2+1
          ! extract subarray from advection field
          advfield_1D_2(i2) = FIELD_2D_AT_I_V1( advfield, i1, i2 )
          ! compute primiti2e of distribution function along this line
          primitive2 (i2) = primitive2 (i2-1) &
               + delta_eta2 * sll_get_df_val( dist_func_2D, i1, i2-1 )
       end do
       call compute_flow_1D_backward( advfield_1D_2, zeros2, 1.0_f64, deltat, &
                                      eta2_min, nc_eta2, delta_eta2,          &
                                      COMPACT, eta2_out ) 
       call compute_spline_1D_hermite( primitive2, nc_eta2, csl_work%spl_eta2 )
       ! interpolate primitive at origin of characteritics
       call interpolate_array_values( eta2_out, primitive2, nc_eta2+1, csl_work%spl_eta2 )
       ! update average value of distribution function in cell using difference of primitives
       do i2 = 1, nc_eta2 
          val = ( primitive2 ( i2+1 ) - primitive2 ( i2 ) ) / delta_eta2
          call sll_set_df_val( dist_func_2D, i1, i2, val )
       end do
    end do

  end subroutine csl_first_order
end module sll_csl
