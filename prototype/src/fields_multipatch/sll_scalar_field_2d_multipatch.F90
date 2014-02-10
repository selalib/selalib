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
  use sll_module_scalar_field_2d_alternative
  use sll_coordinate_transformation_multipatch_module
  use sll_constants
  use sll_module_interpolators_2d_base
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sll_utilities
  use sll_boundary_condition_descriptors
  use sll_gnuplot
  implicit none


  type ::  sll_scalar_field_multipatch_2d
     sll_int32 :: num_patches
     character(len=128) :: field_name
     type(sll_coordinate_transformation_multipatch_2d), pointer    :: transf
     type(sll_scalar_field_2d_discrete_ptr), dimension(:), pointer :: fields &
          => null()
     type(sll_arb_deg_2d_interpolator_ptr), dimension(:), pointer :: interps &
          => null()
     ! to be used if the field is supposed to allocate its own data and so
     ! it would be responsible for deleting it too. Else, it is just a pointer
     ! that needs to be set with the proper access function.
     logical                                         :: owns_memory = .false.
     type(multipatch_data_2d), dimension(:), pointer :: patch_data => null()
   contains
     procedure, pass :: initialize => initialize_scalar_field_sfmp2d
     procedure, pass :: allocate_memory => allocate_memory_sfmp2d
     procedure, pass :: update_interpolation_coefficients => &
          update_interp_coeffs_sfmp2d
     procedure, pass :: get_transformation  => get_patch_transformation_sfmp2d
     procedure, pass :: get_logical_mesh    => get_patch_logical_mesh_sfmp2d
     procedure, pass :: get_jacobian_matrix => get_jacobian_matrix_sfmp2d
     procedure, pass :: get_number_patches  => get_number_patches_sfmp2d
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

  type multipatch_data_2d
     sll_real64, dimension(:,:), pointer :: array => null()
  end type multipatch_data_2d


  interface sll_delete
     module procedure delete_field_sfmp2d_ptr
  end interface sll_delete


contains   ! *****************************************************************


  function new_scalar_field_multipatch_2d( &
    field_name, &
    transformation ) result(obj)

    type(sll_scalar_field_multipatch_2d), pointer              :: obj
    character(len=*), intent(in)                               :: field_name
    type(sll_coordinate_transformation_multipatch_2d), pointer :: transformation
    sll_int32  :: ierr

    SLL_ALLOCATE(obj,ierr)
    call obj%initialize( field_name, transformation )
  end function new_scalar_field_multipatch_2d
  
  subroutine initialize_scalar_field_sfmp2d( &
    fmp, &
    field_name, &
    transf )
    
    class(sll_scalar_field_multipatch_2d)                      :: fmp
    character(len=*), intent(in)                               :: field_name
    type(sll_coordinate_transformation_multipatch_2d), pointer :: transf
    character(len=256)                                         :: patch_name
    type(sll_logical_mesh_2d), pointer                         :: lm
    sll_int32, dimension(1:2)                                  :: connectivity
    character(len=128) :: format_string 
    sll_int32  :: i
    sll_int32  :: num_patches
    sll_int32  :: bc_left
    sll_int32  :: bc_right
    sll_int32  :: bc_bottom
    sll_int32  :: bc_top
    sll_int32  :: ierr

    fmp%transf => transf
    
    ! to build the name of each patch.
    format_string = "(a, a, I0.3)"

    num_patches = transf%get_number_patches()

    fmp%num_patches = num_patches
    SLL_ALLOCATE(fmp%fields(num_patches), ierr)
    SLL_ALLOCATE(fmp%interps(num_patches),ierr)

    ! Build a field object for each individual patch. Interpolator type is
    ! hardwired for the moment. It would be desirable to make this an option.
    ! This could only be properly done when all the compilers we care about
    ! permit us to make arrays of a type which contains a polymorphic pointer.
    do i=0, num_patches-1
       ! create the patch-dedicated interpolator.
       lm=>fmp%transf%get_logical_mesh(i)

       !------------------------------------------------------------------
       !                      WARNING!!!!!!!!
       !------------------------------------------------------------------
       ! WARNING: Note that here the convention established in CAID is 
       ! HARDWIRED.  Specifically, the numbering of the faces in a given patch 
       ! is numbered in CAID as:
       !
       !                       (top)
       !                         2
       !                   +-----------+
       !                   |           |
       !                   |           |
       !      (left)   1   |   patch   |  3   (right)  
       !                   |           |
       !                   |           |
       !                   +-----------+
       !                         0 
       !                      (bottom)
       ! 
       ! Any issues regarding compatibility with faces should be explored here
       ! first. Not updating this commentary if the convention changes would
       ! be a troubling mistake.

       connectivity(:) = fmp%transf%get_connectivity(i,0)
       ! just being paranoid, there is no way that one of the values could be
       ! negative and not the other...
       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          bc_bottom = SLL_DIRICHLET !SLL_HERMITE
       else
          bc_bottom = SLL_DIRICHLET ! THIS IS TEMPORARY, MORE OPTIONS ARE NEEDED
       end if

       connectivity(:) = fmp%transf%get_connectivity(i,1)
       ! just being paranoid, there is no way that one of the values could be
       ! negative and not the other...
       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          bc_left = SLL_DIRICHLET !SLL_HERMITE
       else
          bc_left = SLL_DIRICHLET ! THIS IS TEMPORARY, MORE OPTIONS ARE NEEDED
       end if

       connectivity(:) = fmp%transf%get_connectivity(i,2)
       ! just being paranoid, there is no way that one of the values could be
       ! negative and not the other...
       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          bc_top = SLL_DIRICHLET !SLL_HERMITE
       else
          bc_top = SLL_DIRICHLET ! THIS IS TEMPORARY, MORE OPTIONS ARE NEEDED
       end if

       connectivity(:) = fmp%transf%get_connectivity(i,3)
       ! just being paranoid, there is no way that one of the values could be
       ! negative and not the other...
       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          bc_right = SLL_DIRICHLET !SLL_HERMITE
       else
          bc_right = SLL_DIRICHLET ! THIS IS TEMPORARY, MORE OPTIONS ARE NEEDED
       end if

       ! NOTE THAT THE SPLINE DEGREE 3 WAS HARDWIRED HERE. NEED A DECISION
       ! ABOUT WHERE/WHEN TO SPECIFY THIS.
       write(patch_name, format_string) trim(field_name), "_patch", i
       print *, 'building patch named ', patch_name

       fmp%interps(i+1)%interp => new_arbitrary_degree_spline_interp2d( &
            lm%num_cells1+1, &
            lm%num_cells2+1, &
            lm%eta1_min, &
            lm%eta1_max, &
            lm%eta2_min, &
            lm%eta2_max, &
            bc_left, &
            bc_right, &
            bc_bottom, &
            bc_top, &
            3, &
            3 )   ! <--- HARDWIRED degree of splines, not OK

       fmp%fields(i+1)%f => new_scalar_field_2d_discrete_alt( &
            patch_name, &
            fmp%interps(i+1)%interp, &
            fmp%transf%get_transformation(i), &
            bc_left, &
            bc_right, &
            bc_bottom, &
            bc_top )
    end do
  end subroutine initialize_scalar_field_sfmp2d


  function get_number_patches_sfmp2d( field ) result(res)
    sll_int32 :: res
    class(sll_scalar_field_multipatch_2d), intent(in) :: field
    res = field%num_patches
  end function get_number_patches_sfmp2d


  subroutine delete_field_sfmp2d_ptr( field )
    type(sll_scalar_field_multipatch_2d), pointer :: field
    sll_int32 :: ierr

    call field%delete()
    SLL_DEALLOCATE(field,ierr)
    field => null()
  end subroutine delete_field_sfmp2d_ptr


  subroutine delete_field_sfmp2d( field )
    class(sll_scalar_field_multipatch_2d), intent(inout) :: field
    sll_int32 :: i
    sll_int32 :: num_patches
    sll_int32 :: ierr

    num_patches = field%get_number_patches()
    do i=0,num_patches-1
       call sll_delete(field%fields(i+1)%f)
       call sll_delete(field%interps(i+1)%interp)
    end do

    SLL_DEALLOCATE( field%fields, ierr )
    SLL_DEALLOCATE( field%interps, ierr )

    if( field%owns_memory .eqv. .true. ) then
       do i=0, num_patches-1
          SLL_DEALLOCATE(field%patch_data(i+1)%array,ierr)
       end do
    end if
  end subroutine delete_field_sfmp2d

  subroutine allocate_memory_sfmp2d( field )
    class(sll_scalar_field_multipatch_2d), intent(inout) :: field
    type(sll_logical_mesh_2d), pointer :: lm
    sll_int32 :: i
    sll_int32 :: ierr
    sll_int32 :: num_patches
    sll_int32 :: numpts1
    sll_int32 :: numpts2

    if( field%owns_memory .eqv. .true. ) then
       print *, 'ERROR, allocate_memory_sfmp2d(), the memory for the ', &
            field%field_name, ' object has already been allocated.'
       stop
    end if
    
    num_patches = field%num_patches
    SLL_ALLOCATE(field%patch_data(num_patches),ierr)

    do i=0,num_patches-1
       lm => field%transf%get_logical_mesh(i)
       numpts1 = lm%num_cells1+1
       numpts2 = lm%num_cells2+1
       SLL_ALLOCATE(field%patch_data(i+1)%array(numpts1,numpts2),ierr)
       call field%set_field_data(i,field%patch_data(i+1)%array)
    end do
    field%owns_memory = .true.
  end subroutine allocate_memory_sfmp2d

  subroutine set_field_data_sfmp2d( mp, patch, values )
    class(sll_scalar_field_multipatch_2d), intent(inout) :: mp
    sll_int32, intent(in)                               :: patch
    sll_real64, dimension(:,:), intent(in)              :: values
    type(sll_logical_mesh_2d), pointer                  :: lm
    sll_int32                                           :: numpts1
    sll_int32                                           :: numpts2

    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    lm => mp%transf%get_logical_mesh(patch)
    numpts1 = lm%num_cells1+1
    numpts2 = lm%num_cells2+1

    ! this should be probably changed by something more informative which
    ! writes the name of the field, etc.
    SLL_ASSERT( size(values,1) >= numpts1 )
    SLL_ASSERT( size(values,2) >= numpts2 )
    call mp%fields(patch+1)%f%set_field_data(values)
  end subroutine set_field_data_sfmp2d

  subroutine update_interp_coeffs_sfmp2d( mp, patch )
    class(sll_scalar_field_multipatch_2d), intent(inout) :: mp
    sll_int32, intent(in)                               :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    ! This is a crucial and difficult function in the sense that its role is
    ! to convert a problem that is by nature global to something local. In
    ! other words, the spline coefficients are to be computed on a per-patch
    ! basis needing some kind of compatibility condition at the internal
    ! edges and hopefully this will result in a useful spline reconstruction
    ! of the data competitive with what would have been calculated if the 
    ! calculation had taken place globally.
    !
    ! We need to abstract a function that handles the compatibility between
    ! the patches. This is where most of the work will be for a while...
    call mp%fields(patch+1)%f%update_interpolation_coefficients( )
  end subroutine update_interp_coeffs_sfmp2d

  function get_patch_transformation_sfmp2d( mp, patch ) result(res)
    type(sll_coordinate_transformation_2d_nurbs), pointer :: res
    class(sll_scalar_field_multipatch_2d), intent(in)      :: mp
    sll_int32, intent(in)                                 :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res => mp%transf%get_transformation( patch )
  end function get_patch_transformation_sfmp2d

  function get_patch_logical_mesh_sfmp2d( mp, patch ) result(res)
    type(sll_logical_mesh_2d), pointer                 :: res
    class(sll_scalar_field_multipatch_2d), intent(in)   :: mp
    sll_int32, intent(in)                              :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res => mp%transf%get_logical_mesh(patch)
  end function get_patch_logical_mesh_sfmp2d

  function get_jacobian_matrix_sfmp2d( mp, eta1, eta2, patch ) result(res)
    sll_real64, dimension(2,2) :: res
    class(sll_scalar_field_multipatch_2d), intent(in)   :: mp
    sll_real64, intent(in)                             :: eta1
    sll_real64, intent(in)                             :: eta2
    sll_int32, intent(in)                              :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res(:,:) = mp%transf%jacobian_matrix(eta1,eta2,patch)
  end function get_jacobian_matrix_sfmp2d

  function value_at_pt_sfmp2d( mp, eta1, eta2, patch ) result(res)
    sll_real64                                          :: res
    class(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_real64, intent(in)                              :: eta1
    sll_real64, intent(in)                              :: eta2
    sll_int32, intent(in)                               :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%f%value_at_point(eta1,eta2)
  end function value_at_pt_sfmp2d

  function value_at_indices_sfmp2d( mp, i, j, patch ) result(res)
    sll_real64                                       :: res
    class(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_int32, intent(in)                            :: i
    sll_int32, intent(in)                            :: j
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%f%value_at_indices(i,j)
  end function value_at_indices_sfmp2d

  function first_deriv_eta1_value_at_pt_sfmp2d( mp, eta1, eta2, patch ) &
    result(res)

    sll_real64                                       :: res
    class(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_real64, intent(in)                           :: eta1
    sll_real64, intent(in)                           :: eta2
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%f%first_deriv_eta1_value_at_point(eta1,eta2)
  end function first_deriv_eta1_value_at_pt_sfmp2d

  function first_deriv_eta2_value_at_pt_sfmp2d( mp, eta1, eta2, patch ) &
    result(res)

    sll_real64                                       :: res
    class(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_real64, intent(in)                           :: eta1
    sll_real64, intent(in)                           :: eta2
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%f%first_deriv_eta2_value_at_point(eta1,eta2)
  end function first_deriv_eta2_value_at_pt_sfmp2d
  
  function first_deriv_eta1_value_at_indices_sfmp2d( mp, i, j, patch ) &
    result(res)

    sll_real64                                       :: res
    class(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_int32, intent(in)                            :: i
    sll_int32, intent(in)                            :: j
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%f%first_deriv_eta1_value_at_indices(i,j)
  end function first_deriv_eta1_value_at_indices_sfmp2d

  function first_deriv_eta2_value_at_indices_sfmp2d( mp, i, j, patch ) &
    result(res)

    sll_real64                                       :: res
    class(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_int32, intent(in)                            :: i
    sll_int32, intent(in)                            :: j
    sll_int32, intent(in)                            :: patch
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    res = mp%fields(patch+1)%f%first_deriv_eta2_value_at_indices(i,j)
  end function first_deriv_eta2_value_at_indices_sfmp2d


  subroutine write_to_file_sfmp2d( mp, tag )
    class(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_int32, intent(in)                            :: tag
    sll_int32                                        :: num_patches
    sll_int32                                        :: i

    num_patches = mp%num_patches
    do i=0, num_patches-1
       call mp%fields(i)%f%write_to_file( tag )
    end do
  end subroutine write_to_file_sfmp2d

end module sll_module_scalar_field_2d_multipatch
