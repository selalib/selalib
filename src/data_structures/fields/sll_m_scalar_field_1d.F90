!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_scalar_field_1d
!
!> @author
!> - Aurore
!
! DESCRIPTION:
!
!> @brief
!> Implements the geometry and mesh descriptor types
!>
!>@details
!>
!> This module depends on:
!>    - memory
!>    - precision
!>    - assert
!>    - utilities
!>    - constants
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_m_scalar_field_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_cartesian_meshes, only: &
      sll_t_cartesian_mesh_1d

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   use sll_m_scalar_field_1d_base, only: &
      sll_c_scalar_field_1d_base

   implicit none

   public :: &
      sll_f_new_scalar_field_1d_analytic, &
      sll_f_new_scalar_field_1d_discrete, &
      sll_t_scalar_field_1d_analytic, &
      sll_t_scalar_field_1d_discrete

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_scalar_field_1d_base) :: sll_t_scalar_field_1d_analytic
      procedure(one_var_parametrizable_function), pointer, nopass :: func
      procedure(one_var_parametrizable_function), pointer, nopass :: first_derivative
      sll_real64, dimension(:), pointer        :: params
      character(len=64)                        :: name
      !     sll_int32                          :: plot_counter
      sll_int32 :: bc_left
      sll_int32 :: bc_right
      type(sll_t_cartesian_mesh_1d), pointer   :: mesh
      ! allows to decide if the user put the derivative of the analiytic function: func
      logical :: present_derivative
   contains
      procedure, pass(field) :: init => initialize_scalar_field_1d_analytic
      procedure, pass(field) :: get_cartesian_mesh => &
         get_cartesian_mesh_1d_analytic
      procedure, pass(field) :: value_at_point => value_at_pt_analytic_1d
      procedure, pass(field) :: value_at_indices => value_at_index_analytic_1d
      procedure, pass(field) :: derivative_value_at_point => &
         derivative_value_at_pt_analytic_1d
      procedure, pass(field) :: derivative_value_at_indices => &
         derivative_value_at_index_analytic_1d
      procedure, pass(field) :: set_field_data => set_field_data_analytic_1d
      procedure, pass(field) :: update_interpolation_coefficients => &
         update_interp_coeffs_1d_analytic
      procedure, pass(field) :: write_to_file => write_to_file_analytic_1d
      procedure, pass(field) :: delete => delete_field_1d_analytic
   end type sll_t_scalar_field_1d_analytic

   type, extends(sll_c_scalar_field_1d_base) :: sll_t_scalar_field_1d_discrete
      sll_real64, dimension(:), pointer  :: values
      !sll_real64, dimension(:,:), pointer  :: coeff_spline
      !sll_int32                            :: sz_coeff1
      !sll_int32                            :: sz_coeff2
      character(len=64)                    :: name
      !     sll_int32                            :: plot_counter
      class(sll_c_interpolator_1d), pointer :: interp_1d !!! a implementer
      sll_real64, dimension(:), pointer :: point
      !sll_real64, dimension(:,:), pointer :: point2d
      sll_int32 :: bc_left
      sll_int32 :: bc_right
      type(sll_t_cartesian_mesh_1d), pointer :: mesh
   contains
      procedure, pass(field) :: init => initialize_scalar_field_1d_discrete
      procedure, pass(field) :: get_cartesian_mesh => &
         get_cartesian_mesh_1d_discrete
      procedure, pass(field) :: value_at_point => value_at_pt_discrete_1d
      procedure, pass(field) :: value_at_indices => value_at_index_discrete_1d
      procedure, pass(field) :: derivative_value_at_point => &
         derivative_value_at_pt_discrete_1d
      procedure, pass(field) :: derivative_value_at_indices => &
         derivative_value_at_index_discrete_1d
      procedure, pass(field) :: set_field_data => set_field_data_discrete_1d
      procedure, pass(field) :: update_interpolation_coefficients => &
         update_interp_coeffs_1d_discrete
      procedure, pass(field) :: write_to_file => write_to_file_discrete_1d
      procedure, pass(field) :: delete => delete_field_1d_discrete
   end type sll_t_scalar_field_1d_discrete

   abstract interface
      function one_var_parametrizable_function(eta, params)
         use sll_m_working_precision
         sll_real64 :: one_var_parametrizable_function
         sll_real64, intent(in) :: eta
         sll_real64, dimension(:), intent(in), optional :: params
      end function one_var_parametrizable_function
   end interface

   interface delete
      module procedure delete_field_1d_analytic, delete_field_1d_discrete
   end interface delete

contains   ! *****************************************************************

   ! **************************************************************************
   !
   !                         ANALYTIC CASE
   !
   ! **************************************************************************

   function value_at_pt_analytic_1d(field, eta)
      class(sll_t_scalar_field_1d_analytic), intent(inout) :: field
      sll_real64, intent(in) :: eta
      sll_real64             ::  value_at_pt_analytic_1d
      value_at_pt_analytic_1d = field%func(eta, field%params)
   end function value_at_pt_analytic_1d

   function value_at_index_analytic_1d(field, i)
      class(sll_t_scalar_field_1d_analytic), intent(inout) :: field
      sll_int32, intent(in) :: i
      sll_real64            :: eta
      sll_real64            :: value_at_index_analytic_1d
      eta = field%mesh%eta_min + real(i - 1, f64)*field%mesh%delta_eta
      value_at_index_analytic_1d = field%func(eta, field%params)
   end function value_at_index_analytic_1d

   function derivative_value_at_pt_analytic_1d(field, eta)
      class(sll_t_scalar_field_1d_analytic), intent(inout) :: field
      sll_real64, intent(in) :: eta
      sll_real64             :: derivative_value_at_pt_analytic_1d

      if (field%present_derivative) then
         derivative_value_at_pt_analytic_1d = &
            field%first_derivative(eta, field%params)
      else
         print *, ' first derivative is not given in the initialization'
      end if

   end function derivative_value_at_pt_analytic_1d

   function derivative_value_at_index_analytic_1d(field, i)
      class(sll_t_scalar_field_1d_analytic), intent(inout) :: field
      sll_int32, intent(in) :: i
      sll_real64            :: eta
      sll_real64            :: derivative_value_at_index_analytic_1d

      eta = field%mesh%eta_min + real(i - 1, f64)*field%mesh%delta_eta

      if (field%present_derivative) then
         derivative_value_at_index_analytic_1d = &
            field%first_derivative(eta, field%params)
      else
         print *, ' first derivative is not given in the initialization'
      end if

   end function derivative_value_at_index_analytic_1d

   function sll_f_new_scalar_field_1d_analytic( &
      func, &
      field_name, &
      bc_left, &
      bc_right, &
      mesh, &
      func_params, &
      first_derivative) result(obj)

      type(sll_t_scalar_field_1d_analytic), pointer :: obj
      procedure(one_var_parametrizable_function)      :: func
      procedure(one_var_parametrizable_function), optional :: first_derivative
      character(len=*), intent(in)                    :: field_name
      sll_real64, dimension(:), intent(in), optional, target :: func_params
      sll_int32, intent(in) :: bc_left
      sll_int32, intent(in) :: bc_right
      sll_int32  :: ierr
      type(sll_t_cartesian_mesh_1d), pointer   :: mesh

      SLL_ALLOCATE(obj, ierr)
      call obj%init( &
         func, &
         field_name, &
         bc_left, &
         bc_right, &
         mesh, &
         func_params, &
         first_derivative)
   end function sll_f_new_scalar_field_1d_analytic

   subroutine set_field_data_analytic_1d(field, values)
      class(sll_t_scalar_field_1d_analytic), intent(inout) :: field
      sll_real64, dimension(:), intent(in) :: values
      print *, 'WARNING: set_field_data_analytic_1d(): it is useless to ', &
         'call this function on an analytic scalar field.'
      print *, field%bc_left*values(1)
   end subroutine set_field_data_analytic_1d

   subroutine update_interp_coeffs_1d_analytic(field)
      class(sll_t_scalar_field_1d_analytic), intent(inout) :: field
      print *, 'WARNING: update_interpolation_coefficients_1d_analytic(): ', &
         ' it is useless to call this function on an analytic scalar field.'
      print *, field%bc_left
   end subroutine update_interp_coeffs_1d_analytic

   subroutine delete_field_1d_analytic(field)
      class(sll_t_scalar_field_1d_analytic), intent(inout) :: field
      ! nothing internal do deallocate, just nullify pointers. Can't call
      ! delete on them because the field does not 'own' these data.
      !print*, associated(field%func)
      !print*, associated(field%params)
      nullify (field%func)
      !nullify(field%params) !not associated
   end subroutine delete_field_1d_analytic

   ! For those cases in which handling pointers to field structures is not
   ! convenient, we offer the following initialization.
   subroutine initialize_scalar_field_1d_analytic( &
      field, &
      func, &
      field_name, &
      bc_left, &
      bc_right, &
      mesh, &
      func_params, &
      first_derivative)

      class(sll_t_scalar_field_1d_analytic), intent(out) :: field
      procedure(one_var_parametrizable_function)      :: func
      procedure(one_var_parametrizable_function), optional :: first_derivative
      character(len=*), intent(in)                    :: field_name
      sll_real64, dimension(:), intent(in), optional, target :: func_params
      type(sll_t_cartesian_mesh_1d), pointer   :: mesh

      sll_int32, intent(in) :: bc_left
      sll_int32, intent(in) :: bc_right

      field%func => func
      if (present(func_params)) then
         field%params => func_params
      else
         ! Allocate an empty array in order to avoid passing a non associated
         ! pointer to 'func', which has a non-pointer array dummy argument.
         ! Note that ifort segfaults without the following line!
         ! [YG - 06.10.2015]
         allocate (field%params(0))
      end if
      field%name = trim(field_name)
      field%bc_left = bc_left
      field%bc_right = bc_right
      field%mesh => mesh
      if (present(first_derivative)) then
         field%first_derivative => first_derivative
         field%present_derivative = .TRUE.
      end if
   end subroutine initialize_scalar_field_1d_analytic

   ! The following pair of subroutines are tricky. We want them as general
   ! services by the fields, hence we need this subroutine interface, yet
   ! we would also like a flexibility in how the derivatives are computed.
   ! A general interpolator interface would cover most of the cases, maybe
   ! all. It could be that a finite difference scheme would also work, if
   ! we ignore some of the interpolator services, like the ability to return
   ! values anywhere instead of at the nodes.
   ! For now, this interface would permit to have multiple implementations.
!!$  subroutine compute_eta1_derivative_on_col( field2d, ith_col, deriv_out )
!!$    type(scalar_field_2d), intent(in)    :: field2d
!!$    sll_int32, intent(in)                :: ith_col
!!$    sll_real64, dimension(:),intent(out) :: deriv_out
!!$  end subroutine compute_eta1_derivative_on_col

   function get_cartesian_mesh_1d_analytic(field) result(res)
      class(sll_t_scalar_field_1d_analytic), intent(in) :: field
      type(sll_t_cartesian_mesh_1d), pointer :: res
      res => field%mesh
   end function get_cartesian_mesh_1d_analytic

   subroutine write_to_file_analytic_1d(field, tag)
      class(sll_t_scalar_field_1d_analytic), intent(inout) :: field
      sll_int32, intent(in)                               :: tag
      sll_int32 :: nptsx
      sll_real64, dimension(:), allocatable :: xcoords
      sll_real64, dimension(:), allocatable :: values
      sll_real64                              :: eta
      sll_int32 :: i
      sll_int32 :: ierr
      ! print*, 'passed'
      ! use the logical mesh information to find out the extent of the
      ! domain and allocate the arrays for the plotter.
      !mesh   => field%get_cartesian_mesh()
      !print*, 'passed'
      nptsx = field%mesh%num_cells + 1

      !print*, 'passed',nptsx
      SLL_ALLOCATE(xcoords(nptsx), ierr)
      SLL_ALLOCATE(values(nptsx), ierr)
      ! print*, 'passed'
      ! Fill the arrays with the needed information.
      do i = 1, nptsx
         eta = field%mesh%eta_min + (i - 1)*field%mesh%delta_eta
         xcoords(i) = eta
         values(i) = field%value_at_point(eta)*tag
      end do

      print *, 'not implemented  sll_gnuplot_curv_1d'
!!$    call sll_gnuplot_curv_1d( & ! a implementer
!!$         nptsx, &
!!$         xcoords, &
!!$         values, &
!!$         trim(field%name), &
!!$         tag, &
!!$         ierr )

      SLL_DEALLOCATE_ARRAY(xcoords, ierr)
      SLL_DEALLOCATE_ARRAY(values, ierr)
   end subroutine write_to_file_analytic_1d

   ! **************************************************************************
   !
   !                         DISCRETE CASE
   !
   ! **************************************************************************

   function sll_f_new_scalar_field_1d_discrete( &
      ! array_1d, &
      field_name, &
      interpolator_1d, &
      bc_left, &
      bc_right, &
      mesh, &
      point_1d, &
      sz_point) result(obj)

      type(sll_t_scalar_field_1d_discrete), pointer :: obj
!    sll_real64, dimension(:), intent(in), target  :: array_1d
      character(len=*), intent(in)                    :: field_name
      class(sll_c_interpolator_1d), target        :: interpolator_1d ! a implementer
      type(sll_t_cartesian_mesh_1d), pointer   :: mesh
      sll_real64, dimension(:), optional :: point_1d
      sll_int32, optional :: sz_point
      ! sll_real64, dimension(:,:), optional :: point2d
      sll_int32, intent(in) :: bc_left
      sll_int32, intent(in) :: bc_right
      sll_int32  :: ierr

      SLL_ALLOCATE(obj, ierr)
      call obj%init( &
         ! array_1d, &
         field_name, &
         interpolator_1d, &
         bc_left, &
         bc_right, &
         mesh, &
         point_1d, &
         sz_point)
   end function sll_f_new_scalar_field_1d_discrete

   subroutine initialize_scalar_field_1d_discrete( &
      field, &
      ! array_1d, &
      field_name, &
      interpolator_1d, &
      bc_left, &
      bc_right, &
      mesh, &
      point_1d, &
      sz_point)

      class(sll_t_scalar_field_1d_discrete)         :: field
      ! sll_real64, dimension(:), intent(in), target  :: array_1d
      character(len=*), intent(in)                    :: field_name
      class(sll_c_interpolator_1d), target        :: interpolator_1d
      type(sll_t_cartesian_mesh_1d), pointer   :: mesh
      sll_real64, dimension(:), optional :: point_1d
      sll_int32, optional :: sz_point
      sll_int32, intent(in) :: bc_left
      sll_int32, intent(in) :: bc_right
      sll_int32 :: ierr

      ! field%values => array_1d
      field%interp_1d => interpolator_1d
      !    field%mesh%written = .false.
      field%name = trim(field_name)
      field%bc_left = bc_left
      field%bc_right = bc_right
      field%mesh => mesh

      !SLL_ALLOCATE(point(sz_point),ierr)
      SLL_ALLOCATE(field%values(field%mesh%num_cells + 1), ierr)

!!$    if ( present( point_1d) .and. present(sz_point)) then
!!$       print*, ' not implemented yet in initialize_scalar_field_1d_discrete'
!!$       stop
!!$    end if
!!$    if ( present( point_1d) .and. .not. present(sz_point)) then
!!$       print*, ' problem presence of point_1d and not the size in '
!!$       print*, 'initialize_scalar_field_1d_discrete'
!!$       stop
!!$    end if

!!$    call field%interp_1d%compute_interpolants( &  !  a implementer
!!$         array_1d, &
!!$         point_1d, &
!!$         sz_point )

      return
      SLL_ASSERT(present(point_1d))
      SLL_ASSERT(present(sz_point))

   end subroutine initialize_scalar_field_1d_discrete

   subroutine set_field_data_discrete_1d(field, values)
      class(sll_t_scalar_field_1d_discrete), intent(inout) :: field
      sll_real64, dimension(:), intent(in) :: values
      if (size(field%values, 1) .ne. size(values, 1)) then
         print *, 'WARNING, set_field_data_discrete_1d(), passed array ', &
            'is not of the size originally declared for this field.'
      end if
      field%values(:) = values(:)
   end subroutine set_field_data_discrete_1d

   ! need to do something about deallocating the field proper, when allocated
   ! in the heap...
   subroutine delete_field_1d_discrete(field)
      class(sll_t_scalar_field_1d_discrete), intent(inout) :: field
      ! just nullify pointers, nothing to deallocate that this object owns.
      !print*, associated(field%values)
      !print*, associated(field%interp_1d)
      !print*, associated(field%point)

      nullify (field%values)
      nullify (field%interp_1d)
      !nullify(field%point) !not associated
   end subroutine delete_field_1d_discrete

   subroutine update_interp_coeffs_1d_discrete(field)
      class(sll_t_scalar_field_1d_discrete), intent(inout) :: field
      call field%interp_1d%compute_interpolants(field%values)
   end subroutine update_interp_coeffs_1d_discrete

   function get_cartesian_mesh_1d_discrete(field) result(res)
      class(sll_t_scalar_field_1d_discrete), intent(in) :: field
      type(sll_t_cartesian_mesh_1d), pointer :: res
      res => field%mesh
   end function get_cartesian_mesh_1d_discrete

   function value_at_pt_discrete_1d(field, eta)
      class(sll_t_scalar_field_1d_discrete), intent(inout) :: field
      sll_real64, intent(in) :: eta
      sll_real64             :: value_at_pt_discrete_1d

      value_at_pt_discrete_1d = field%interp_1d%interpolate_from_interpolant_value(eta)
   end function value_at_pt_discrete_1d

   function value_at_index_discrete_1d(field, i)
      class(sll_t_scalar_field_1d_discrete), intent(inout) :: field
      sll_int32, intent(in) :: i
      sll_real64            :: eta
      sll_real64            :: value_at_index_discrete_1d
      eta = field%mesh%eta_min + real(i - 1, f64)*field%mesh%delta_eta
      value_at_index_discrete_1d = field%interp_1d%interpolate_from_interpolant_value(eta)
   end function value_at_index_discrete_1d

   function derivative_value_at_pt_discrete_1d(field, eta)
      class(sll_t_scalar_field_1d_discrete), intent(inout) :: field
      sll_real64, intent(in) :: eta
      sll_real64             :: derivative_value_at_pt_discrete_1d

      derivative_value_at_pt_discrete_1d = &
         field%interp_1d%interpolate_from_interpolant_derivative_eta1(eta)
   end function derivative_value_at_pt_discrete_1d

   function derivative_value_at_index_discrete_1d(field, i)
      class(sll_t_scalar_field_1d_discrete), intent(inout) :: field
      sll_int32, intent(in) :: i
      sll_real64            :: eta
      sll_real64            :: derivative_value_at_index_discrete_1d
      eta = field%mesh%eta_min + real(i - 1, f64)*field%mesh%delta_eta
      derivative_value_at_index_discrete_1d = &
         field%interp_1d%interpolate_from_interpolant_derivative_eta1(eta)
   end function derivative_value_at_index_discrete_1d

   subroutine write_to_file_discrete_1d(field, tag)
      class(sll_t_scalar_field_1d_discrete), intent(inout) :: field
      sll_int32, intent(in)                               :: tag
      sll_int32 :: nptsx
      sll_real64, dimension(:), allocatable :: xcoords
      sll_real64, dimension(:), allocatable :: values
      sll_real64                              :: eta
      type(sll_t_cartesian_mesh_1d), pointer      :: mesh
      sll_int32 :: i
      sll_int32 :: ierr

      ! use the logical mesh information to find out the extent of the
      ! domain and allocate the arrays for the plotter.
      mesh => field%get_cartesian_mesh()
      nptsx = mesh%num_cells + 1

      SLL_ALLOCATE(xcoords(nptsx), ierr)
      SLL_ALLOCATE(values(nptsx), ierr)

      ! Fill the arrays with the needed information.
      do i = 1, nptsx
         eta = mesh%eta_min + (i - 1)*mesh%delta_eta
         xcoords(i) = eta!field%x(eta)
         values(i) = field%value_at_point(eta)*tag
      end do

      print *, 'not implement the sll_gnuplot_curv_1d '
!!$    call sll_gnuplot_curv_1d( &
!!$         nptsx, &
!!$         xcoords, &
!!$         values, &
!!$         trim(field%name), &
!!$         tag, &
!!$         ierr )

      SLL_DEALLOCATE_ARRAY(xcoords, ierr)
      SLL_DEALLOCATE_ARRAY(values, ierr)
   end subroutine write_to_file_discrete_1d

end module sll_m_scalar_field_1d
