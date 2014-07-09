!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_distribution_function_4d_multipatch
!
!> @author
!> - Edwin
!
! DESCRIPTION: 
!
!> @brief
!> Encapsulates a group of distribution functions, each associated with a 
!> patch.
!>
!>@details
!>
!
! REVISION HISTORY:
! 08 jul 2014 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_distribution_function_4d_multipatch_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_file_io.h"
  use sll_coordinate_transformation_multipatch_module
  use sll_constants
  use sll_remapper
  use sll_module_interpolators_2d_base
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sll_utilities
  use sll_boundary_condition_descriptors
  use sll_gnuplot
  implicit none


  type ::  sll_distribution_function_4d_multipatch
     sll_int32 :: num_patches
     sll_int32 :: npts_x3
     sll_int32 :: npts_x4
     type(sll_coordinate_transformation_multipatch_2d), pointer :: transf
     type(data_4d_ptr), dimension(:), pointer :: distributions
   contains
     procedure, pass(df) :: allocate_memory => allocate_memory_df_4d_mp
  end type sll_distribution_function_4d_multipatch

  type :: data_4d_ptr
     sll_real64, dimension(:,:,:,:), pointer :: f
  end type data_4d_ptr

contains

  function sll_new_distribution_function_4d_multipatch( &
       transf_mp, &
       npts_x3, &
       npts_x4 ) result(df)

    type(sll_distribution_function_4d_multipatch), pointer :: df
    type(sll_coordinate_transformation_multipatch_2d), intent(in), target:: &
         transf_mp
    sll_int32 :: npts_x3
    sll_int32 :: npts_x4
    sll_int32 :: ierr

    SLL_ALLOCATE(df,ierr)
    df%transf => transf_mp
    df%npts_x3 = npts_x3
    df%npts_x4 = npts_x4
    df%num_patches = transf_mp%get_number_patches()
    call df%allocate_memory()
  end function sll_new_distribution_function_4d_multipatch

  subroutine allocate_memory_df_4d_mp( df )
    class(sll_distribution_function_4d_multipatch), intent(inout) :: df
    sll_int32 :: ierr
    sll_int32 :: i
    sll_int32 :: np_x1
    sll_int32 :: np_x2
    sll_int32 :: np_x3
    sll_int32 :: np_x4
    type(sll_logical_mesh_2d), pointer :: lm

    np_x3 = df%npts_x3
    np_x4 = df%npts_x4

    SLL_ALLOCATE( df%distributions(df%num_patches), ierr )

    do i=1,df%num_patches
       np_x1 = df%transf%get_num_cells_eta1(i)
       np_x2 = df%transf%get_num_cells_eta2(i)
       SLL_ALLOCATE(df%distributions(i)%f(np_x1,np_x2,np_x3,np_x4),ierr)
    end do
  end subroutine allocate_memory_df_4d_mp

end module sll_distribution_function_4d_multipatch_module
