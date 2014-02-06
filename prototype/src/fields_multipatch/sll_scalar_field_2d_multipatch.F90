!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_scalar_field_2d_multpatch
!
!> @author
!> - Edwin
!
! DESCRIPTION: 
!
!> @brief
!> Encapsulates a group of scalar fields, each associated with a particular
!> region of a physical domain. The ensemble is referred to as 'multipatch'.
!>
!>@details
!>
!> The use of the 2D scalar field, multipatch version, is a direct extension of
!> the use of the 2D scalar field, except that its interface routines include
!> and additional parameter, the patch number under inquiry (a number between
!> 0 and the number of patches - 1) plus additional information related with
!> the connectivity between the patches.
!
! REVISION HISTORY:
! 04 Feb 2014 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_module_scalar_field_2d_multipatch
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_file_io.h"
!  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative
  use sll_constants
  use sll_module_interpolators_2d_base
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sll_utilities
  use sll_boundary_condition_descriptors
  use sll_gnuplot
  implicit none


  type ::  sll_scalar_field_multipatch_2d
     sll_int32 :: num_patches
     type(sll_coordinate_transformation_multipatch_2d), pointer    :: transf
     type(sll_scalar_field_2d_discrete_ptr), dimension(:), pointer :: fields &
          => null()
   contains
     procedure, pass :: initialize_field => initialize_scalar_field_sfmp2d
     procedure, pass :: update_interpolation_coefficients => &
          update_interp_coeffs_sfmp2d
     procedure, pass :: get_transformation  => get_patch_transformation_sfmp2d
     procedure, pass :: get_logical_mesh    => get_patch_logical_mesh_sfmp2d
     procedure, pass :: get_jacobian_matrix => get_jacobian_matrix_sfmp2d
     procedure, pass :: value_at_point      => value_at_pt_sfmp2d
     procedure, pass :: value_at_indices    => value_at_indices_sfmp2d
     procedure, pass :: first_deriv_eta1_value_at_point => &
          first_deriv_eta1_value_at_pt_sfmp2d
     procedure, pass :: first_deriv_eta2_value_at_point => &
          first_deriv_eta2_value_at_pt_sfmp2d
     procedure, pass :: first_deriv_eta1_value_at_indices => &
          first_deriv_eta1_value_at_indices_sfmp2d
     procedure, pass :: first_deriv_eta2_value_at_indices => &
          first_deriv_eta2_value_at_indices_sfmp2d
     procedure, pass :: set_field_data => set_field_data_sfmp2d
     procedure, pass :: write_to_file => write_to_file_sfmp2d
     procedure, pass :: delete => delete_field_sfmp2d
  end type sll_scalar_field_multipatch_2d


  interface sll_delete
     module procedure delete_field_sfmp2d
  end interface sll_delete


contains   ! *****************************************************************


  function new_scalar_field_multipatch_2d( &
    field_name, &
    transformation ) result(obj)

    type(sll_scalar_field_multipatch_2d), pointer :: obj
    character(len=*), intent(in)                  :: field_name
    type(sll_coordinate_transformation_multipatch_2d), pointer :: transf
    sll_int32  :: num_patches
    sll_int32  :: ierr
    
    num_patches = transf%get_number_patches()

    SLL_ALLOCATE(obj,ierr)
    obj%num_patches = num_patches
    SLL_ALLOCATE(obj%fields(num_patches), ierr)
    ! do we have all the information to initialize the fields?
    ! - we know the connectivities, thus we know what boundary conditions to
    !   apply for each patch.
    ! - we can choose directly the interpolator type.
    ! - the transformation is available
    ! - BC's also
    ! So here we go...

  end function new_scalar_field_multipatch_2d
  
  subroutine initialize_scalar_field_sfmp2d( &
    mp, &
    patch )
    
    class(sll_scalar_field_2d_discrete_alt)              :: field
    character(len=*), intent(in)                         :: field_name
    class(sll_interpolator_2d_base), target              :: interpolator_2d
    class(sll_coordinate_transformation_2d_base), target :: transformation
    type(sll_logical_mesh_2d), pointer  :: m2d    
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: bc_bottom
    sll_int32, intent(in) :: bc_top
    sll_real64, dimension(:), optional :: point1_1d
    sll_real64, dimension(:), optional :: point2_1d
    sll_int32,optional :: sz_point1
    sll_int32,optional :: sz_point2
    !sll_int32 :: i
    sll_int32 :: ierr   

    m2d => transformation%mesh
    field%T => transformation
    field%interp_2d => interpolator_2d
    !    field%mesh%written = .false.
    field%name      = trim(field_name)
    field%bc_left   = bc_left
    field%bc_right  = bc_right
    field%bc_bottom = bc_bottom
    field%bc_top    = bc_top

    ! Allocate internal array to store locally a copy of the data.
    SLL_ALLOCATE(field%values(m2d%num_cells1+1,m2d%num_cells2+1),ierr)    
    !print*,'first line',  field%values(1,:)
    !print*, 'second line', field%values(2,:)
!!$    call field%interp_2d%compute_interpolants( &
!!$         field%values, & !array_2d, &
!!$         point1_1d, &
!!$         sz_point1, &
!!$         point2_1d, &
!!$         sz_point2 )

  end subroutine initialize_scalar_field_2d_discrete_alt
  
  ! need to do something about deallocating the field proper, when allocated
  ! in the heap...
  subroutine delete_field_2d_discrete_alt( field )
    class(sll_scalar_field_2d_discrete_alt), intent(inout) :: field
    sll_int32 :: ierr
    if(associated(field%values))    SLL_DEALLOCATE(field%values,ierr)
    if(associated(field%T))         nullify(field%T)
    if(associated(field%interp_2d)) nullify(field%interp_2d)
    if(associated(field%point1_1d)) nullify(field%point1_1d)
    if(associated(field%point2_1d)) nullify(field%point2_1d)
  end subroutine delete_field_2d_discrete_alt

  subroutine set_field_data_sfmp2d( mp, patch, values )
    type(sll_scalar_field_multipatch_2d), intent(inout) :: mp
    sll_int32, intent(in)                               :: patch
    sll_real64, dimension(:,:), intent(in)              :: values
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    call mp%fields(patch+1)%f%set_field_data(values)
  end subroutine set_field_data_sfmp2d

  subroutine update_interp_coeffs_sfmp2d( mp, patch )
    type(sll_scalar_field_multipatch_2d), intent(inout) :: mp
    sll_int32, intent(in)                               :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    call mp%fields(patch+1)%f%update_interpolation_coefficients( )
  end subroutine update_interp_coeffs_sfmp2d

  function get_patch_transformation_sfmp2d( mp, patch ) result(res)
    type(sll_coordinate_transformation_2d_nurbs), pointer :: res
    type(sll_scalar_field_multipatch_2d), intent(in)      :: mp
    sll_int32, intent(in)                                 :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    call mp%transf%get_transformation( patch )
  end function get_patch_transformation_sfmp2d

  function get_patch_logical_mesh_sfmp2d( mp, patch ) result(res)
    type(sll_logical_mesh_2d), pointer                 :: res
    type(sll_scalar_field_multipatch_2d), intent(in)   :: mp
    sll_int32, intent(in)                              :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res => mp%transf%get_logical_mesh(patch)
  end function get_patch_logical_mesh_sfmp2d

  function get_jacobian_matrix_sfmp2d( mp, eta1, eta2, patch ) result(res)
    sll_real64, dimension(2,2) :: res
    type(sll_scalar_field_multipatch_2d), intent(in)   :: mp
    sll_real64, intent(in)                             :: eta1
    sll_real64, intent(in)                             :: eta2
    sll_int32, intent(in)                              :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res(:,:) = mp%transf%jacobian_matrix(eta1,eta2,patch)
  end function get_jacobian_matrix_sfmp2d

  function value_at_pt_sfmp2d( mp, eta1, eta2, patch ) result(res)
    sll_real64                                          :: res
    type(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_real64, intent(in)                              :: eta1
    sll_real64, intent(in)                              :: eta2
    sll_int32, intent(in)                               :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%value_at_point(eta1,eta2)
  end function value_at_pt_sfmp2d

  function value_at_indices_sfmp2d( mp, i, j, patch ) result(res)
    sll_real64                                       :: res
    type(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_int32, intent(in)                            :: i
    sll_int32, intent(in)                            :: j
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%value_at_indices(i,j)
  end function value_at_indices_sfmp2d

  function first_deriv_eta1_value_at_pt_sfmp2d( mp, eta1, eta2, patch ) &
    result(res)

    sll_real64                                       :: res
    type(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_real64, intent(in)                           :: eta1
    sll_real64, intent(in)                           :: eta2
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%first_deriv_eta1_value_at_point(eta1,eta2)
  end function first_deriv_eta1_value_at_pt_sfmp2d

  function first_deriv_eta2_value_at_pt_sfmp2d( mp, eta1, eta2, patch ) &
    result(res)

    sll_real64                                       :: res
    type(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_real64, intent(in)                           :: eta1
    sll_real64, intent(in)                           :: eta2
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%first_deriv_eta2_value_at_point(eta1,eta2)
  end function first_deriv_eta2_value_at_pt_sfmp2d
  
  function first_deriv_eta1_value_at_indices_sfmp2d( mp, i, j, patch ) &
    result(res)

    sll_real64                                       :: res
    type(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_int32, intent(in)                            :: i
    sll_int32, intent(in)                            :: j
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%first_deriv_eta1_value_at_indices(i,j)
  end function first_deriv_eta1_value_at_indices_sfmp2d

  function first_deriv_eta2_value_at_indices_sfmp2d( mp, i, j, patch ) &
    result(res)

    sll_real64                                       :: res
    type(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_int32, intent(in)                            :: i
    sll_int32, intent(in)                            :: j
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%first_deriv_eta2_value_at_indices(i,j)
  end function first_deriv_eta2_value_at_indices_sfmp2d


  subroutine write_to_file_sfmp2d( mp, tag )
    type(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_int32, intent(in)                            :: tag
    sll_int32                                        :: num_patches
    sll_int32                                        :: i

    num_patches = mp%num_patches
    do i=0, num_patches-1
       call mp%fields(i)%write_to_file( tag )
    end do
  end subroutine write_to_file_sfmp2d

end module sll_module_scalar_field_2d_multipatch
