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
     type(multipatch_data_2d_real), dimension(:), pointer :: patch_data => null()
     ! The following memory buffers are meant to mimic the type of 
     ! organization that we will need for the parallel case. The actual 
     ! size of these buffers depend on the degree of the spline used to 
     ! reconstruct the field data and of course the logical mesh size.
     type(multipatch_data_2d_real), dimension(:), pointer :: buffers0 => null()
     type(multipatch_data_2d_real), dimension(:), pointer :: buffers1 => null()
     type(multipatch_data_2d_real), dimension(:), pointer :: buffers2 => null()
     type(multipatch_data_2d_real), dimension(:), pointer :: buffers3 => null()
     type(multipatch_data_2d_real), dimension(:), pointer :: derivs0 => null()
     type(multipatch_data_2d_real), dimension(:), pointer :: derivs1 => null()
     type(multipatch_data_2d_real), dimension(:), pointer :: derivs2 => null()
     type(multipatch_data_2d_real), dimension(:), pointer :: derivs3 => null()
   contains
     procedure, pass :: initialize => initialize_scalar_field_sfmp2d
     procedure, pass :: allocate_memory => allocate_memory_sfmp2d
     procedure, pass :: update_interpolation_coefficients => &
          update_interp_coeffs_sfmp2d
     procedure, pass :: get_transformation  => get_patch_transformation_sfmp2d
     procedure, pass :: get_logical_mesh    => get_patch_logical_mesh_sfmp2d
     procedure, pass :: get_jacobian_matrix => get_jacobian_matrix_sfmp2d
     procedure, pass :: get_number_patches  => get_number_patches_sfmp2d
     ! These are the functions to aid the finite element calculation  
     !procedure, pass :: get_spline_local_index=> get_spline_local_index_sfmp2d
     !procedure, pass :: get_spline_global_index=> get_spline_global_index_sfmp2d
     !procedure, pass :: get_spline_local_to_global_index=>get_spline_local_to_global_index_sfmp2d
     procedure, pass :: value_at_point      => value_at_pt_sfmp2d
     procedure, pass :: value_at_indices    => value_at_indices_sfmp2d
     procedure, pass :: set_value_at_indices => set_value_at_indices_sfmp2d
     procedure, pass :: first_deriv_eta1_value_at_point => &
          first_deriv_eta1_value_at_pt_sfmp2d
     procedure, pass :: first_deriv_eta2_value_at_point => &
          first_deriv_eta2_value_at_pt_sfmp2d
     procedure, pass :: first_deriv_eta1_value_at_indices => &
          first_deriv_eta1_value_at_indices_sfmp2d
     procedure, pass :: first_deriv_eta2_value_at_indices => &
          first_deriv_eta2_value_at_indices_sfmp2d
     procedure, pass :: get_patch_data_pointer => get_patch_data_ptr_sfmp2d
     procedure, pass :: set_patch_data_pointer => set_patch_data_ptr_sfmp2d
     procedure, pass :: set_field_data => set_field_data_sfmp2d
     procedure, pass :: write_to_file => write_to_file_sfmp2d
     procedure, pass :: delete => delete_field_sfmp2d
  end type sll_scalar_field_multipatch_2d


  interface sll_delete
     module procedure delete_field_sfmp2d_ptr
  end interface sll_delete


contains   ! *****************************************************************


  function new_scalar_field_multipatch_2d( &
    field_name, &
    transformation ) result(obj)

    type(sll_scalar_field_multipatch_2d), pointer              :: obj
    character(len=*), intent(in)                               :: field_name
    type(sll_coordinate_transformation_multipatch_2d), target :: transformation
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
    type(sll_coordinate_transformation_multipatch_2d), target :: transf
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
    sll_int32  :: num_pts1
    sll_int32  :: num_pts2
    sll_int32  :: ierr

    fmp%field_name = field_name
    fmp%transf => transf
    
    ! to build the name of each patch.
    format_string = "(a, a, I0.3)"

    num_patches = transf%get_number_patches()

    fmp%num_patches = num_patches
    SLL_ALLOCATE(fmp%patch_data(num_patches),ierr)
    SLL_ALLOCATE(fmp%fields(num_patches), ierr)
    SLL_ALLOCATE(fmp%interps(num_patches),ierr)

    ! Build a field object for each individual patch. Interpolator type is
    ! hardwired for the moment. It would be desirable to make this an option.
    ! This could only be properly done when all the compilers we care about
    ! permit us to make arrays of a type which contains a polymorphic pointer.
    print *, 'proceeding to create the patches...'
    do i=0, num_patches-1
       ! create the patch-dedicated interpolator.
       ! The following 'fix' is just designed to support gfortran 4.6.3, once
       ! this is not an issue, this should be changed when it is decided not to
       ! support this compiler anymore.
       !       lm=>fmp%transf%get_logical_mesh(i)
       lm => fmp%transf%transfs(i+1)%t%mesh
       print *, 'extracted logical mesh from patch ', i
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
          bc_bottom = SLL_HERMITE
       else
          bc_bottom = SLL_HERMITE!SLL_DIRICHLET ! THIS IS TEMPORARY, MORE OPTIONS ARE NEEDED
       end if

       connectivity(:) = fmp%transf%get_connectivity(i,1)
       ! just being paranoid, there is no way that one of the values could be
       ! negative and not the other...
       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          bc_left = SLL_HERMITE
       else
          bc_left = SLL_HERMITE!SLL_DIRICHLET ! THIS IS TEMPORARY, MORE OPTIONS ARE NEEDED
       end if

       connectivity(:) = fmp%transf%get_connectivity(i,2)
       ! just being paranoid, there is no way that one of the values could be
       ! negative and not the other...
       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          bc_top = SLL_HERMITE
       else
          bc_top = SLL_HERMITE!SLL_DIRICHLET ! THIS IS TEMPORARY, MORE OPTIONS ARE NEEDED
       end if

       connectivity(:) = fmp%transf%get_connectivity(i,3)
       ! just being paranoid, there is no way that one of the values could be
       ! negative and not the other...
       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          bc_right = SLL_HERMITE
       else
          bc_right = SLL_HERMITE!SLL_DIRICHLET ! THIS IS TEMPORARY, MORE OPTIONS ARE NEEDED
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
       print *, 'created interpolator for patch ', i
       fmp%fields(i+1)%f => new_scalar_field_2d_discrete_alt( &
            patch_name, &
            fmp%interps(i+1)%interp, &
            fmp%transf%get_transformation(i), &
            bc_left, &
            bc_right, &
            bc_bottom, &
            bc_top )
       print *, 'created field associated to patch ', i
    end do

    ! Allocate the memory needed to work with the patch compatibility
    ! algorithms. Each buffer depends on the size of the logical mesh 
    ! associated with a given patch border.
    SLL_ALLOCATE(fmp%buffers0(num_patches),ierr)
    SLL_ALLOCATE(fmp%buffers1(num_patches),ierr)
    SLL_ALLOCATE(fmp%buffers2(num_patches),ierr)
    SLL_ALLOCATE(fmp%buffers3(num_patches),ierr)
    SLL_ALLOCATE(fmp%derivs0(num_patches),ierr)
    SLL_ALLOCATE(fmp%derivs1(num_patches),ierr)
    SLL_ALLOCATE(fmp%derivs2(num_patches),ierr)
    SLL_ALLOCATE(fmp%derivs3(num_patches),ierr)

    ! WARNING: this is temporary, then number of derivatives to be estimated
    ! should come from the degree of the spline interpolation used for the
    ! field data. For the cubic splines we only specify the first derivative,
    ! but this should be extended once we are more confident of the 
    ! soundness of this methodology.
    
#define NUM_DERIVS 1

    do i=1,num_patches
       ! Get rid of the following fix once it is concluded that gfortran 4.6
       ! should not be supported.
       !       lm => fmp%transf%get_logical_mesh(i-1)
       lm => fmp%transf%transfs(i)%t%mesh
       num_pts1 = lm%num_cells1 + 1
       num_pts2 = lm%num_cells2 + 1
       SLL_ALLOCATE(fmp%buffers0(i)%array(num_pts1,NUM_DERIVS),ierr)
       SLL_ALLOCATE(fmp%buffers1(i)%array(num_pts2,NUM_DERIVS),ierr)
       SLL_ALLOCATE(fmp%buffers2(i)%array(num_pts1,NUM_DERIVS),ierr)
       SLL_ALLOCATE(fmp%buffers3(i)%array(num_pts2,NUM_DERIVS),ierr)
       ! The calculation of the compatibility conditions between patches is
       ! being written with the future parallelization in mind, hence the
       ! redundant calculation of the slopes...
       SLL_ALLOCATE(fmp%derivs0(i)%array(num_pts1,NUM_DERIVS),ierr)
       SLL_ALLOCATE(fmp%derivs1(i)%array(num_pts2,NUM_DERIVS),ierr)
       SLL_ALLOCATE(fmp%derivs2(i)%array(num_pts1,NUM_DERIVS),ierr)
       SLL_ALLOCATE(fmp%derivs3(i)%array(num_pts2,NUM_DERIVS),ierr)
    end do
    print *, 'initialize_scalar_field_sfmp2d(): finished initialization.'
  end subroutine initialize_scalar_field_sfmp2d


  function get_number_patches_sfmp2d( field ) result(res)
    sll_int32 :: res
    class(sll_scalar_field_multipatch_2d), intent(in) :: field
    res = field%num_patches
  end function get_number_patches_sfmp2d


  subroutine delete_field_sfmp2d_ptr( field )
    class(sll_scalar_field_multipatch_2d), pointer :: field
    sll_int32 :: ierr

    call field%delete()
    SLL_DEALLOCATE(field, ierr)
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
       SLL_DEALLOCATE_ARRAY(field%buffers0(i+1)%array,ierr)
       SLL_DEALLOCATE_ARRAY(field%buffers1(i+1)%array,ierr)
       SLL_DEALLOCATE_ARRAY(field%buffers2(i+1)%array,ierr)
       SLL_DEALLOCATE_ARRAY(field%buffers3(i+1)%array,ierr)
       SLL_DEALLOCATE_ARRAY(field%derivs0(i+1)%array,ierr)
       SLL_DEALLOCATE_ARRAY(field%derivs1(i+1)%array,ierr)
       SLL_DEALLOCATE_ARRAY(field%derivs2(i+1)%array,ierr)
       SLL_DEALLOCATE_ARRAY(field%derivs3(i+1)%array,ierr)
    end do

    SLL_DEALLOCATE( field%fields, ierr )
    SLL_DEALLOCATE( field%interps, ierr )

    if( field%owns_memory .eqv. .true. ) then
       do i=0, num_patches-1
          SLL_DEALLOCATE(field%patch_data(i+1)%array,ierr)
       end do
    end if

       SLL_DEALLOCATE(field%buffers0,ierr)
       SLL_DEALLOCATE(field%buffers1,ierr)
       SLL_DEALLOCATE(field%buffers2,ierr)
       SLL_DEALLOCATE(field%buffers3,ierr)
       SLL_DEALLOCATE(field%derivs0,ierr)
       SLL_DEALLOCATE(field%derivs1,ierr)
       SLL_DEALLOCATE(field%derivs2,ierr)
       SLL_DEALLOCATE(field%derivs3,ierr)
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

    do i=0,num_patches-1
       ! get rid of the following 'fix'whenever gfortran 4.6 is not supported
       ! by Selalib
       !       lm => field%transf%get_logical_mesh(i)
       lm => field%transf%transfs(i+1)%t%mesh
       numpts1 = lm%num_cells1+1
       numpts2 = lm%num_cells2+1
       SLL_ALLOCATE(field%patch_data(i+1)%array(numpts1,numpts2),ierr)
       call field%set_field_data(i,field%patch_data(i+1)%array)
    end do
    field%owns_memory = .true.
    ! And link each patch with the newly allocated memory.
    do i=0,num_patches-1
       call field%fields(i+1)%f%set_field_data(field%patch_data(i+1)%array)
    end do

    ! Link each patch with the newly allocated memory. There is a problem here:
    ! Each upon the call to 'set_field_data()', each field COPIES the data into
    ! a local memory managed by the field. The more desirable behavior would
    ! simply be to reset the pointer of the field so that it knows instantly
    ! if the data is altered.
    do i=0,num_patches-1
       call field%fields(i+1)%f%free_internal_data_copy()
       call field%fields(i+1)%f%reset_data_pointer(field%patch_data(i+1)%array)
    end do
  end subroutine allocate_memory_sfmp2d

  subroutine set_field_data_sfmp2d( mp, patch, values )
    class(sll_scalar_field_multipatch_2d), intent(inout) :: mp
    sll_int32, intent(in)                               :: patch
    sll_real64, dimension(:,:), intent(in)              :: values
    type(sll_logical_mesh_2d), pointer                  :: lm
    sll_int32                                           :: numpts1
    sll_int32                                           :: numpts2

    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    ! Get rid of the following 'fix' whenever it is decided that gfortran 4.6
    ! is not supported.
    !    lm => mp%transf%get_logical_mesh(patch)
    lm => mp%transf%transfs(patch+1)%t%mesh
    numpts1 = lm%num_cells1+1
    numpts2 = lm%num_cells2+1

    ! this should be probably changed by something more informative which
    ! writes the name of the field, etc.
    SLL_ASSERT( size(values,1) >= numpts1 )
    SLL_ASSERT( size(values,2) >= numpts2 )
    call mp%fields(patch+1)%f%set_field_data(values)
  end subroutine set_field_data_sfmp2d

  function get_patch_data_ptr_sfmp2d( mp, patch ) result(ptr)
    sll_real64, dimension(:,:), pointer                  :: ptr
    class(sll_scalar_field_multipatch_2d), intent(inout) :: mp
    sll_int32, intent(in)                                :: patch
    ptr => mp%fields(patch+1)%f%get_data_pointer()
  end function get_patch_data_ptr_sfmp2d

  subroutine set_patch_data_ptr_sfmp2d( mp, patch, ptr )
    sll_real64, dimension(:,:), pointer                  :: ptr
    class(sll_scalar_field_multipatch_2d), intent(inout) :: mp
    sll_int32, intent(in)                                :: patch
    mp%patch_data(patch+1)%array => ptr
  end subroutine set_patch_data_ptr_sfmp2d

  subroutine update_interp_coeffs_sfmp2d( mp )
    class(sll_scalar_field_multipatch_2d), intent(inout) :: mp
    sll_int32 :: ipatch
    ! WARNING: this step would be unnecessary if the fields referred to 
    ! their data via a pointer, the following call will be copying lots of
    ! arrays...
    do ipatch=0,mp%num_patches-1
       call mp%fields(ipatch+1)%f%set_field_data(mp%patch_data(ipatch+1)%array)
    end do
    ! patches must agree on the compatibility in internal borders before
    ! individually launching the update coefficients per patch.
    call compute_compatible_derivatives_in_borders(mp)
    do ipatch=0,mp%num_patches-1
       call mp%fields(ipatch+1)%f%update_interpolation_coefficients( )
    end do
  end subroutine update_interp_coeffs_sfmp2d


  ! THIS SHOULD BE CHANGED TO INCLUDE THE CASE IN WHICH THE CELL SPACING
  ! CHANGES IN BETWEEN PATCHES!!!
  subroutine compute_compatible_derivatives_in_borders( fmp )
    class(sll_scalar_field_multipatch_2d), intent(inout) :: fmp
    type(sll_logical_mesh_2d), pointer :: m
    sll_int32 :: num_patches
    sll_int32 :: ip
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: num_pts1
    sll_int32 :: num_pts2
    sll_int32 :: other_patch
    sll_int32 :: other_face
    sll_int32 :: current_face
    sll_real64 :: rdelta1
    sll_real64 :: rdelta2
    sll_int32, dimension(2) :: connectivity
    sll_real64, dimension(:,:), pointer :: d
    sll_real64, dimension(:,:), pointer :: buf
    sll_real64, dimension(:,:), pointer :: this_buffer
    sll_real64, dimension(:,:), pointer :: other_buffer
    sll_real64, dimension(:,:), pointer :: derivs

    num_patches = fmp%num_patches

    do ip=1,num_patches
       ! Please get rid of this awful 'fix' whenever it is decided that 
       ! gfortran 4.6 should not be supported by Selalib
       !       m => fmp%transf%get_logical_mesh(ip-1)
       m => fmp%transf%transfs(ip)%t%mesh
       num_pts1 = m%num_cells1 + 1
       num_pts2 = m%num_cells2 + 1
       rdelta1  = 1.0_f64/m%delta_eta1
       rdelta2  = 1.0_f64/m%delta_eta2
       d   => fmp%patch_data(ip)%array

       ! Here we should have a select case and adjust the calculation to 
       ! the degree of the derivatives which should be calculated. Presently
       ! we only put the cubic spline compatibility condition which only 
       ! requires the calculation of the first derivative.

       ! ---------------------------------------------------------------------
       !
       ! Compute local contribution to face 0 (south face).
       !
       ! ---------------------------------------------------------------------
       connectivity(:) = fmp%transf%get_connectivity(ip-1,0)
       buf => fmp%buffers0(ip)%array

       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          ! this face is connected, thus the calculation is 'shared' between
          ! the patches. The stencil used is:
          !         
          ! f'_(0) = (1/h)(-(5/3)f_(-3)+(3/20)f_(-2)-(3/4)f_(-1) +
          !                 (3/4)f_( 1)-(3/20)f_( 2)+(5/3)f_( 3) )
          !
          ! Thus note that for the moment this is assuming that h is the same
          ! on both patches...
          do j=1,1 ! <--- this would be the number of derivatives to calculate
             do i=1,num_pts1 
                ! only compute the first derivative contribution for now
                buf(i,j) = ((3.0_f64/4.0_f64 )*d(i,2) - &
                            (3.0_f64/20.0_f64)*d(i,3) + &
                            (5.0_f64/3.0_f64 )*d(i,4))*rdelta2
             end do
          end do
       else 
          ! the face is not connected. Use different stencil:
          !
          ! f'_(0) = (1/h)(-2.45*f_(0) +  6*f_(1) -    7.5*f_(2) + (20/3)*f_(3) 
          !                -3.75*f_(4) + 1.2*f(5) - (5/30)*f_(6))
          do j=1,1 ! <--- this would be the number of derivatives to calculate
             do i=1,num_pts1 
                ! only compute the first derivative contribution for now
                buf(i,j) = (-2.45_f64*d(i,1) + 6.0_f64*d(i,2) -7.5_f64*d(i,3) &
                            + (20.0_f64/3.0_f64)*d(i,4) - 3.75*d(i,5) &
                            + 1.2*d(i,6) - (5.0_f64/30.0_f64 )*d(i,7))*rdelta2
             end do
          end do
       end if

       ! ---------------------------------------------------------------------
       !
       ! Compute local contribution to face 2 (north face).
       !
       ! ---------------------------------------------------------------------
       connectivity(:) = fmp%transf%get_connectivity(ip-1,2)
       buf => fmp%buffers2(ip)%array

       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          ! this face is connected, thus the calculation is 'shared' between
          ! the patches. The stencil used is:
          !         
          ! f'_(0) = (1/h)(-(5/3)f_(-3)+(3/20)f_(-2)-(3/4)f_(-1) +
          !                 (3/4)f_( 1)-(3/20)f_( 2)+(5/3)f_( 3) )
          !
          ! Thus note that for the moment this is assuming that h is the same
          ! on both patches...
          do j=1,1 ! <--- this would be the number of derivatives to calculate
             do i=1,num_pts1 
                ! only compute the first derivative contribution for now
                buf(i,j) = (- (5.0_f64/3.0_f64 )*d(i,num_pts2-3) &
                            + (3.0_f64/20.0_f64)*d(i,num_pts2-2) &
                            - (3.0_f64/4.0_f64 )*d(i,num_pts2-1) )*rdelta2
             end do
          end do
       else 
          ! the face is not connected. Use different stencil:
          !
          ! f'_(0) = (1/h)(+2.45*f_(0) -  6*f_(-1) + 7.5*f_(-2) - (20/3)*f_(-3) 
          !                +3.75*f_(-4) -1.2*f(-5) + (5/30)*f_(-6))
          do j=1,1 ! <--- this would be the number of derivatives to calculate
             do i=1,num_pts1 
                ! only compute the first derivative contribution for now
                buf(i,j) = (+ (5.0_f64/30.0_f64)*d(i,num_pts2-6) &
                            - 1.2               *d(i,num_pts2-5) &
                            + 3.75              *d(i,num_pts2-4) &
                            - (20.0_f64/3.0_f64)*d(i,num_pts2-3) &
                            + 7.5_f64           *d(i,num_pts2-2) &
                            - 6.0_f64           *d(i,num_pts2-1) & 
                            + 2.45_f64          *d(i,num_pts2) )*rdelta2
             end do
          end do
       end if

       ! ---------------------------------------------------------------------
       !
       ! Compute local contribution to face 1 (west face).
       !
       ! ---------------------------------------------------------------------
       connectivity(:) = fmp%transf%get_connectivity(ip-1,1)
       buf => fmp%buffers1(ip)%array

       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          ! this face is connected, thus the calculation is 'shared' between
          ! the patches. The stencil used is:
          !         
          ! f'_(0) = (1/h)(-(5/3)f_(-3)+(3/20)f_(-2)-(3/4)f_(-1) +
          !                 (3/4)f_( 1)-(3/20)f_( 2)+(5/3)f_( 3) )
          !
          ! Thus note that for the moment this is assuming that h is the same
          ! on both patches...
          do j=1,1 ! <--- this would be the number of derivatives to calculate
             do i=1,num_pts2 
                ! only compute the first derivative contribution for now
                buf(i,j) = ((3.0_f64/4.0_f64 )*d(2,i) - &
                            (3.0_f64/20.0_f64)*d(3,i) + &
                            (5.0_f64/3.0_f64 )*d(4,i))*rdelta1
             end do
          end do
       else 
          ! the face is not connected. Use different stencil:
          !
          ! f'_(0) = (1/h)(-2.45*f_(0) +  6*f_(1) -    7.5*f_(2) + (20/3)*f_(3) 
          !                -3.75*f_(4) + 1.2*f(5) - (5/30)*f_(6))
          do j=1,1 ! <--- this would be the number of derivatives to calculate
             do i=1,num_pts2
                ! only compute the first derivative contribution for now
                buf(i,j) = (-2.45_f64*d(1,i) + 6.0_f64*d(2,i) -7.5_f64*d(3,i) &
                            + (20.0_f64/3.0_f64)*d(4,i) - 3.75*d(5,i) &
                            + 1.2*d(6,i) - (5.0_f64/30.0_f64 )*d(7,i))*rdelta1
             end do
          end do
       end if

       ! ---------------------------------------------------------------------
       !
       ! Compute local contribution to face 3 (east face).
       !
       ! ---------------------------------------------------------------------
       connectivity(:) = fmp%transf%get_connectivity(ip-1,3)
       buf => fmp%buffers3(ip)%array

       if( (connectivity(1) >= 0) .and. (connectivity(2) >= 0) ) then
          ! this face is connected, thus the calculation is 'shared' between
          ! the patches. The stencil used is:
          !         
          ! f'_(0) = (1/h)(-(5/3)f_(-3)+(3/20)f_(-2)-(3/4)f_(-1) +
          !                 (3/4)f_( 1)-(3/20)f_( 2)+(5/3)f_( 3) )
          !
          ! Thus note that for the moment this is assuming that h is the same
          ! on both patches...
          do j=1,1 ! <--- this would be the number of derivatives to calculate
             do i=1,num_pts2 
                ! only compute the first derivative contribution for now
                buf(i,j) = (- (5.0_f64/3.0_f64 )*d(num_pts1-3,i) &
                            + (3.0_f64/20.0_f64)*d(num_pts1-2,i) &
                            - (3.0_f64/4.0_f64 )*d(num_pts1-1,i) )*rdelta1
             end do
          end do
       else 
          ! the face is not connected. Use different stencil:
          !
          ! f'_(0) = (1/h)(+2.45*f_(0) -  6*f_(-1) + 7.5*f_(-2) - (20/3)*f_(-3) 
          !                +3.75*f_(-4) -1.2*f(-5) + (5/30)*f_(-6))
          do j=1,1 ! <--- this would be the number of derivatives to calculate
             do i=1,num_pts2 
                ! only compute the first derivative contribution for now
                buf(i,j) = (+ (5.0_f64/30.0_f64)*d(num_pts1-6,i) &
                            - 1.2               *d(num_pts1-5,i) &
                            + 3.75              *d(num_pts1-4,i) &
                            - (20.0_f64/3.0_f64)*d(num_pts1-3,i) &
                            + 7.5_f64           *d(num_pts1-2,i) &
                            - 6.0_f64           *d(num_pts1-1,i) & 
                            + 2.45_f64          *d(num_pts1  ,i) )*rdelta1
             end do
          end do
       end if
    end do

    ! The buffers contain the contribution from each patch to the derivative
    ! calculation for the faces which are connected and the actual 
    ! derivative of the data for those faces which aren't.
    ! Proceed to assemble the derivatives information for each border.

    do ip=1,num_patches
       do current_face=0,3
          connectivity(:) = fmp%transf%get_connectivity(ip-1,current_face)
          other_patch = connectivity(1)
          other_face  = connectivity(2)
          ! select the derivs and local buffers
          select case ( current_face )
          case(0)
             derivs      => fmp%derivs0(ip)%array
             this_buffer => fmp%buffers0(ip)%array
          case(1)
             derivs      => fmp%derivs1(ip)%array
             this_buffer => fmp%buffers1(ip)%array
          case(2)
             derivs      => fmp%derivs2(ip)%array
             this_buffer => fmp%buffers2(ip)%array
          case(3)
             derivs      => fmp%derivs3(ip)%array
             this_buffer => fmp%buffers3(ip)%array
          end select

          if( (other_patch <  0) .and. (other_face <  0) ) then
             ! This is an exterior face, the slopes have been computed 
             ! locally already.
             derivs(:,:) = this_buffer(:,:)
          else
             ! Face shared between patches, need to combine two buffers.
             this_buffer => fmp%buffers0(ip)%array
             select case ( other_face )
             case(0)
                other_buffer => fmp%buffers0(other_patch+1)%array
             case(1)
                other_buffer => fmp%buffers1(other_patch+1)%array
             case(2)
                other_buffer => fmp%buffers2(other_patch+1)%array
             case(3)
                other_buffer => fmp%buffers3(other_patch+1)%array
             case default
                print *, 'ERROR in ', __FILE__, __LINE__,'face connectivity ', &
                     'found not between the expected integers 0 and 3.'
                stop
             end select
             derivs(:,:) = other_buffer(:,:) + this_buffer(:,:)
          end if
       end do
    end do
  end subroutine compute_compatible_derivatives_in_borders

  function get_patch_transformation_sfmp2d( mp, patch ) result(res)
    class(sll_coordinate_transformation_2d_nurbs), pointer :: res
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

  subroutine set_value_at_indices_sfmp2d( mp, i, j, patch, val )
    class(sll_scalar_field_multipatch_2d), intent(in) :: mp
    sll_int32, intent(in)                             :: i
    sll_int32, intent(in)                             :: j
    sll_int32, intent(in)                             :: patch
    sll_real64, intent(in)                            :: val
    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
    mp%patch_data(patch+1)%array(i,j) = val
  end subroutine set_value_at_indices_sfmp2d

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
    sll_int32                                        :: gnu_id
    sll_int32                                        :: error
    character(len=4)                                 :: ctag

    num_patches = mp%num_patches
    do i=1, num_patches
        call mp%fields(i)%f%write_to_file( tag )
    end do

    call sll_new_file_id(gnu_id, error)

    call int2string(tag,ctag)
    open(gnu_id,file=trim(mp%field_name)//".gnu")
    write(gnu_id,"(a)",advance='no') &
        "splot '"//trim(mp%fields(1)%f%name)//"_"//ctag//".dat' w l"
    if (num_patches > 1) then
       do i=2, num_patches
           write(gnu_id,"(a)",advance='no') &
        ", '"//trim(mp%fields(i)%f%name)//"_"//ctag//".dat' w l"
       end do
       write(gnu_id,*)
       close(gnu_id)
    end if

  end subroutine write_to_file_sfmp2d

  subroutine set_slope_mp(mp,patch,slope_left,&
       slope_right,&
       slope_bottom,&
       slope_top)
    
    sll_real64, dimension(:),optional :: slope_left
    sll_real64, dimension(:),optional :: slope_right
    sll_real64, dimension(:),optional :: slope_bottom
    sll_real64, dimension(:),optional :: slope_top
    class(sll_scalar_field_multipatch_2d), intent(in) :: mp
    class(arb_deg_2d_interpolator),pointer:: interpolator
    sll_int32 :: patch
    sll_int32 :: num_pts1,num_pts2

    interpolator => mp%interps(patch+1)%interp

    num_pts1 = interpolator%num_pts1
    num_pts2 = interpolator%num_pts2

    if ( (present(slope_left)) .and. &
         (present(slope_right)) .and. &
         (present(slope_bottom)) .and. &
         (present(slope_top)) ) then

       call set_slope2d(&
         interpolator,&
         slope_left,&
         slope_right,&
         slope_bottom,&
         slope_top)
    else
       
       ! ATTENTION !
       ! see convention above about the face numbering
       ! 
       call set_slope2d(&
            interpolator,&
            mp%derivs1(patch+1)%array(1:num_pts2,1),&
            mp%derivs3(patch+1)%array(1:num_pts2,1),&
            mp%derivs0(patch+1)%array(1:num_pts1,1),&
            mp%derivs2(patch+1)%array(1:num_pts1,1))
    end if

  end subroutine set_slope_mp


  !! normally in the Ahmed code local_index_array
  !! is a 2D array with a dimension ( num_splines_loc, num_cell) in a patch
  !! where 
  !! num_splines_loc = (spline_degre_eta1 +1)*(spline_degre_eta2 +1)
  !! and 
  !! num_cell is the number of cells in the patch i
  !! we stocke this 2D array such as a 1D array with a dimension
  !! num_splines_loc*num_cell
  !! i.e. k + (l-1)* num_cell for k = 1,num_cell and 
  !! l = 1, num_splines_loc
!!$  function get_spline_local_index_sfmp2d(mp,patch,splines_local,cell_i,cell_j)result(res)
!!$    class(sll_scalar_field_multipatch_2d), intent(inout) :: mp
!!$    sll_int32 :: patch
!!$    sll_int32 :: splines_local
!!$    sll_int32 :: cell_i,cell_j
!!$    sll_int32 :: num_cell
!!$    sll_int32 :: res
!!$    sll_int32 :: num_spline_loc_max
!!$    sll_int32 :: total_num_cells_in_patch
!!$    type(sll_logical_mesh_2d), pointer                         :: lm
!!$    
!!$   
!!$    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
!!$    SLL_ASSERT( (cell_i >= 1) .and. (cell_i < mp%fields(patch+1)%f%mesh%num_cells1) )
!!$    SLL_ASSERT( (cell_j >= 1) .and. (cell_j < mp%fields(patch+1)%f%mesh%num_cells2) )
!!$    SLL_ASSERT( (splines_local >= 1) .and. (splines_local < num_spline_loc_max) )
!!$
!!$    lm=>mp%transf%get_logical_mesh(patch)
!!$    num_spline_loc_max = (mp%interps(patch+1)%interp%spline_degree1 +1)*&
!!$                         (mp%interps(patch+1)%interp%spline_degree2 +1)
!!$    num_cell = cell_i + (cell_j-1)*lm%num_cells1
!!$    total_num_cells_in_patch = lm%num_cells1*lm%num_cells2
!!$    res = mp%transf%local_indices(patch+1)%array(num_cell+ (splines_local-1)*total_num_cells_in_patch)
!!$  end function get_spline_local_index_sfmp2d
!!$
!!$
!!$
!!$  !! global_indices contains the numeration of splines in all the domain
!!$  !! It is a 1D array with the following size
!!$  !!  sum_i num_splines_in_domain_patch_i 
!!$  !! where  num_splines_in_domain_patch_i=(number_cells_eta1_patch_i + deg_spline_eta1_patch_i)
!!$  !!                                   * (number_cells_eta2_patch_i + deg_spline_eta2_patch_i) 
!!$  
!!$  function get_spline_global_index_sfmp2d(mp,num_splines_global)result(res)
!!$    class(sll_scalar_field_multipatch_2d), intent(inout) :: mp
!!$    sll_int32 :: num_splines_global
!!$    sll_int32 :: res
!!$    
!!$    res = mp%transf%global_indices(num_splines_global)
!!$  end function get_spline_global_index_sfmp2d
!!$
!!$
!!$  !! normally in the Ahmed code local_to_global_index_array
!!$  !! is a 2D array with a dimension ( num_splines_loc, num_cell) in a patch
!!$  !! where 
!!$  !! num_splines_loc = (spline_degre_eta1 +1)*(spline_degre_eta2 +1)
!!$  !! and 
!!$  !! num_cell is the number of cells in the patch i
!!$  !! we stocke this 2D array such as a 1D array with a dimension
!!$  !! num_splines_loc*num_cell
!!$  !! i.e. k + (l-1)* num_cell for k = 1,num_cell and 
!!$  !! l = 1, num_splines_loc
!!$  function get_spline_local_to_global_index_sfmp2d(mp,patch,splines_local,cell_i,cell_j)&
!!$       result(res)
!!$    class(sll_scalar_field_multipatch_2d), intent(inout) :: mp
!!$    sll_int32 :: patch
!!$    sll_int32 :: splines_local
!!$    sll_int32 :: cell_i,cell_j
!!$    sll_int32 :: num_cell
!!$    sll_int32 :: res
!!$    sll_int32 :: num_spline_loc_max
!!$    sll_int32 :: index
!!$    sll_int32 :: total_num_cells_in_patch
!!$    type(sll_logical_mesh_2d), pointer                         :: lm
!!$    
!!$    
!!$    SLL_ASSERT( (patch >= 0) .and. (patch < mp%num_patches) )
!!$    SLL_ASSERT( (cell_i >= 1) .and. (cell_i < mp%fields(patch+1)%f%mesh%num_cells1) )
!!$    SLL_ASSERT( (cell_j >= 1) .and. (cell_j < mp%fields(patch+1)%f%mesh%num_cells2) )
!!$   
!!$    
!!$    lm=>mp%transf%get_logical_mesh(patch)
!!$    num_spline_loc_max = (mp%interps(patch+1)%interp%spline_degree1 +1)*&
!!$         (mp%interps(patch+1)%interp%spline_degree2 +1)
!!$
!!$    SLL_ASSERT( (splines_local >= 1) .and. (splines_local < num_spline_loc_max) )
!!$    num_cell = cell_i + (cell_j-1)*lm%num_cells1
!!$    total_num_cells_in_patch = lm%num_cells1*lm%num_cells2
!!$    index = num_cell+(splines_local-1)*total_num_cells_in_patch
!!$    res=mp%transf%local_to_global_ind(patch+1)%array(index)
!!$  end function get_spline_local_to_global_index_sfmp2d
  

end module sll_module_scalar_field_2d_multipatch
