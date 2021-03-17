!> @ingroup coordinate_transformations
!> @brief
!> Analytic coordinate transformation
!> @details
!> There are two main types of coordinate transformations: analytic and
!> discrete. In the first case the transformation can be
!> specified analytically, through the two functions:
!>
!> \f[
!> \begin{matrix}
!>                     x_1 = x_1(\eta_1,\eta_2)  \\\\
!>                     x_2 = x_2(\eta_1,\eta_2)
!> \end{matrix}
!> \f]
!>
!> Where both, eta1 and eta2 should be defined on the intervals that define
!> the extent of the logical mesh (default values in logical mesh are  [0,1].
!> The same transformation
!> can be specified by the set of transformed points \f$x_1(i,j), x_2(i,j),\f$ as
!> two 2D arrays or 1D arrays that describe the transformation on each
!> direction.
!>
!> The transformation is also represented by the Jacobian matrix:
!> \f[
!> J(\eta_1,\eta_2) =
!> \begin{bmatrix}
!> \partial x_1(\eta_1,\eta_2) / \partial \eta_1 &
!> \partial x_1(\eta_1,\eta_2) / \partial \eta_2 \\\\
!> \partial x_2(\eta_1,\eta_2) / \partial \eta_1 &
!> \partial x_2(\eta_1,\eta_2) / \partial \eta_2
!> \end{bmatrix}
!> \f]
module sll_m_coordinate_transformations_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_dirichlet

   use sll_m_cartesian_meshes, only: &
      sll_f_new_cartesian_mesh_2d, &
      sll_t_cartesian_mesh_2d, &
      sll_o_delete

   use sll_m_coordinate_transformation_2d_base, only: &
      sll_c_coordinate_transformation_2d_base, &
      sll_p_io_gnuplot, &
      sll_p_io_mtv, &
      sll_p_io_xdmf, &
      sll_i_transformation_func_nopass

   use sll_m_gnuplot, only: &
      sll_o_gnuplot_2d

   use sll_m_interpolators_2d_base, only: &
      sll_c_interpolator_2d

   use sll_m_plotmtv, only: &
      sll_o_plotmtv_write

   use sll_m_utilities, only: &
      sll_s_new_file_id

   use sll_m_xdmf, only: &
      sll_s_xdmf_close, &
      sll_o_xdmf_open, &
      sll_o_xdmf_write_array

   implicit none

   public :: &
      sll_f_new_coordinate_transformation_2d_analytic, &
      sll_f_new_coordinate_transformation_2d_discrete, &
      sll_s_coordinate_transformation_2d_analytic_init, &
      sll_s_coordinate_transformation_2d_discrete_init, &
      sll_t_coordinate_transformation_2d_analytic, &
      sll_t_coordinate_transformation_2d_discrete, &
      sll_o_delete

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   abstract interface
      function j_matrix_f_nopass(eta1, eta2, params) result(val)
         use sll_m_working_precision
         sll_real64, dimension(2, 2)   :: val
         sll_real64   :: eta1
         sll_real64   :: eta2
         sll_real64, dimension(:), optional, intent(in) :: params
      end function j_matrix_f_nopass
   end interface

   !> @brief
   !> Functions array
   !> @details
   !> Here we try to represent the Jacobian matrix an actual 2D array of
   !> functions. But since fortran does not allow arrays of pointers, here
   !> we define a special type that can be used as an array element.
   type jacobian_matrix_element
      procedure(sll_i_transformation_func_nopass), pointer, nopass :: f
   end type jacobian_matrix_element

!> Analytic transformation
   type, extends(sll_c_coordinate_transformation_2d_base):: &
      sll_t_coordinate_transformation_2d_analytic
!!$     sll_real64, dimension(:,:), pointer :: x1_node   ! x1(i,j)
!!$     sll_real64, dimension(:,:), pointer :: x2_node   ! x2(i,j)
      !character(len=64) :: label
      !logical           :: written! = .false.
      !type(sll_t_cartesian_mesh_2d), pointer :: mesh => null()
      !>
      type(jacobian_matrix_element), dimension(2, 2) :: j_matrix
      !>
      procedure(sll_i_transformation_func_nopass), pointer, nopass :: x1_func  ! user
      !>
      procedure(sll_i_transformation_func_nopass), pointer, nopass :: x2_func  ! user
      !>
      procedure(two_arg_message_passing_func_analyt), pointer, pass :: jacobian_func
      !>
      procedure(j_matrix_f_nopass), pointer, nopass :: jacobian_matrix_function
      !>
      sll_real64, dimension(:), pointer :: params => null() ! transf params
   contains
      !>
      procedure, pass(transf) :: init => sll_s_coordinate_transformation_2d_analytic_init
      !>
      procedure, pass(transf) :: get_cartesian_mesh => get_cartesian_mesh_analytic
      ! Functions with integer arguments
      !>
      procedure, pass(transf) :: x1_at_node => x1_node_analytic
      !>
      procedure, pass(transf) :: x2_at_node => x2_node_analytic
      !>
      procedure, pass(transf) :: jacobian_at_node => jacobian_node_analytic
      !>
      procedure, pass(transf) :: x1_at_cell => x1_cell_analytic
      !>
      procedure, pass(transf) :: x2_at_cell => x2_cell_analytic
      !>
      procedure, pass(transf) :: jacobian_at_cell => jacobian_2d_cell_analytic
      ! Functions with real arguments
      !>
      procedure, pass(transf) :: x1 => x1_analytic
      !>
      procedure, pass(transf) :: x2 => x2_analytic
      !>
      procedure, pass(transf) :: jacobian => jacobian_2d_analytic
      !>
      procedure, pass(transf) :: jacobian_matrix => jacobian_matrix_2d_analytic
      !>
      procedure, pass(transf) :: inverse_jacobian_matrix => &
         inverse_jacobian_matrix_2d_analytic
      !>
      procedure, pass(transf) :: write_to_file => write_to_file_2d_analytic
      !>
      procedure, pass(transf) :: read_from_file => read_from_file_2d_analytic
      !>
      procedure, pass(transf) :: delete => delete_transformation_2d_analytic
   end type sll_t_coordinate_transformation_2d_analytic

   ! -----------------------------------------------------------------------
   !
   !                          Discrete case
   !
   ! -----------------------------------------------------------------------

   type, extends(sll_c_coordinate_transformation_2d_base) :: &
      sll_t_coordinate_transformation_2d_discrete
      !>
      sll_real64, dimension(:, :), pointer :: x1_node => null()   ! x1(i,j)
      !>
      sll_real64, dimension(:, :), pointer :: x2_node => null()  ! x2(i,j)
      !>
      sll_real64, dimension(:, :), pointer :: x1_cell => null()
      !>
      sll_real64, dimension(:, :), pointer :: x2_cell => null()
      !>
      sll_real64, dimension(:, :), pointer :: jacobians_n => null()
      !>
      sll_real64, dimension(:, :), pointer :: jacobians_c => null()
!     type(jacobian_matrix_element), dimension(2,2) :: j_matrix
      !>
      class(sll_c_interpolator_2d), pointer               :: x1_interp
      !>
      class(sll_c_interpolator_2d), pointer               :: x2_interp
      !type(sll_t_cartesian_mesh_2d), pointer :: mesh => null()
   contains
      !>
      procedure, pass(transf) :: init => &
         sll_s_coordinate_transformation_2d_discrete_init
      !>
      procedure, pass(transf) :: get_cartesian_mesh => get_cartesian_mesh_discrete
      !>
      procedure, pass(transf) :: x1_at_node => x1_node_discrete
      !>
      procedure, pass(transf) :: x2_at_node => x2_node_discrete
      !>
      procedure, pass(transf) :: jacobian_at_node => transf_2d_jacobian_node_discrete
      !>
      procedure, pass(transf) :: x1 => x1_discrete
      !>
      procedure, pass(transf) :: x2 => x2_discrete
      !>
      procedure, pass(transf) :: x1_at_cell => x1_cell_discrete
      !>
      procedure, pass(transf) :: x2_at_cell => x2_cell_discrete
      !>
      procedure, pass(transf) :: jacobian_at_cell => jacobian_2d_cell_discrete
      !>
      procedure, pass(transf) :: jacobian => jacobian_2d_discrete
      !>
      procedure, pass(transf) :: jacobian_matrix => jacobian_matrix_2d_discrete
      !>
      procedure, pass(transf) :: inverse_jacobian_matrix => &
         inverse_jacobian_matrix_2d_discrete
      !>
      procedure, pass(transf) :: write_to_file => write_to_file_2d_discrete
      !>
      procedure, pass(transf) :: read_from_file => read_from_file_2d_discrete
      !>
      procedure, pass(transf) :: delete => delete_transformation_2d_discrete
   end type sll_t_coordinate_transformation_2d_discrete

#ifndef DOXYGEN_SHOULD_SKIP_THIS

   abstract interface
      function two_arg_message_passing_func_discr(transf, eta1, eta2)
         use sll_m_working_precision
         import     :: sll_t_coordinate_transformation_2d_discrete
         sll_real64                      :: two_arg_message_passing_func_discr
         class(sll_t_coordinate_transformation_2d_discrete)  :: transf
         sll_real64, intent(in)          :: eta1
         sll_real64, intent(in)          :: eta2
      end function two_arg_message_passing_func_discr
   end interface

   abstract interface
      function two_arg_message_passing_func_analyt(transf, eta1, eta2)
         use sll_m_working_precision
         import     :: sll_t_coordinate_transformation_2d_analytic
         sll_real64                      :: two_arg_message_passing_func_analyt
         class(sll_t_coordinate_transformation_2d_analytic) :: transf
         sll_real64, intent(in)          :: eta1
         sll_real64, intent(in)          :: eta2
      end function two_arg_message_passing_func_analyt
   end interface

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

   !> Deallocate
   interface sll_o_delete
      module procedure &
         delete_transformation_2d_analytic, &
         delete_transformation_2d_discrete
   end interface sll_o_delete

contains

   !**************************************************************************
   !
   !       Functions for the analytic coordinate transformation
   !
   !**************************************************************************

   !> Create the analytic coordinate transformation
   function sll_f_new_coordinate_transformation_2d_analytic( &
      label, &
      mesh_2d, &
      x1_func, &
      x2_func, &
      j11_func, &
      j12_func, &
      j21_func, &
      j22_func, &
      params)

      type(sll_t_coordinate_transformation_2d_analytic), pointer :: &
         sll_f_new_coordinate_transformation_2d_analytic
      character(len=*), intent(in)                  :: label
      type(sll_t_cartesian_mesh_2d), target         :: mesh_2d
      procedure(sll_i_transformation_func_nopass)   :: x1_func
      procedure(sll_i_transformation_func_nopass)   :: x2_func
      procedure(sll_i_transformation_func_nopass)   :: j11_func
      procedure(sll_i_transformation_func_nopass)   :: j12_func
      procedure(sll_i_transformation_func_nopass)   :: j21_func
      procedure(sll_i_transformation_func_nopass)   :: j22_func
      sll_real64, dimension(:), intent(in)          :: params
      sll_int32                                     :: ierr

      SLL_ALLOCATE(sll_f_new_coordinate_transformation_2d_analytic, ierr)

      call sll_s_coordinate_transformation_2d_analytic_init( &
         sll_f_new_coordinate_transformation_2d_analytic, &
         label, &
         mesh_2d, &
         x1_func, &
         x2_func, &
         j11_func, &
         j12_func, &
         j21_func, &
         j22_func, &
         params)

   end function sll_f_new_coordinate_transformation_2d_analytic

   subroutine sll_s_coordinate_transformation_2d_analytic_init( &
      transf, &
      label, &
      mesh_2d, &
      x1_func, &
      x2_func, &
      j11_func, &
      j12_func, &
      j21_func, &
      j22_func, &
      params)

      class(sll_t_coordinate_transformation_2d_analytic), intent(inout) :: &
         transf
      character(len=*), intent(in)                   :: label
      procedure(sll_i_transformation_func_nopass)    :: x1_func
      procedure(sll_i_transformation_func_nopass)    :: x2_func
      procedure(sll_i_transformation_func_nopass)    :: j11_func
      procedure(sll_i_transformation_func_nopass)    :: j12_func
      procedure(sll_i_transformation_func_nopass)    :: j21_func
      procedure(sll_i_transformation_func_nopass)    :: j22_func
      type(sll_t_cartesian_mesh_2d), target          :: mesh_2d
      sll_real64, dimension(:), intent(in)           :: params
      sll_int32  :: ierr

      transf%label = trim(label)
      transf%mesh => mesh_2d

      ! Assign the transformation functions and parameters
      transf%x1_func => x1_func
      transf%x2_func => x2_func
      if (size(params) > 0) then
         SLL_ALLOCATE(transf%params(size(params)), ierr)
         transf%params(:) = params(:)
      else
         SLL_CLEAR_ALLOCATE(transf%params(1), ierr)
      end if

      ! Fill the jacobian matrix
      transf%j_matrix(1, 1)%f => j11_func
      transf%j_matrix(1, 2)%f => j12_func
      transf%j_matrix(2, 1)%f => j21_func
      transf%j_matrix(2, 2)%f => j22_func
      transf%jacobian_func => jacobian_2d_analytic

   end subroutine sll_s_coordinate_transformation_2d_analytic_init

   subroutine delete_transformation_2d_analytic(transf)
      class(sll_t_coordinate_transformation_2d_analytic), intent(inout) :: transf

      nullify (transf%x1_func)
      nullify (transf%x2_func)
      nullify (transf%jacobian_func)
      nullify (transf%jacobian_matrix_function)
   end subroutine delete_transformation_2d_analytic

   function get_cartesian_mesh_analytic(transf) result(res)
      class(sll_t_coordinate_transformation_2d_analytic), intent(in) :: transf
      class(sll_t_cartesian_mesh_2d), pointer :: res
      res => transf%mesh
   end function get_cartesian_mesh_analytic

   function jacobian_2d_analytic(transf, eta1, eta2) result(val)
      sll_real64                        :: val
      class(sll_t_coordinate_transformation_2d_analytic) :: transf
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      sll_real64             :: j11
      sll_real64             :: j12
      sll_real64             :: j21
      sll_real64             :: j22
      j11 = (transf%j_matrix(1, 1)%f(eta1, eta2, transf%params))
      j12 = (transf%j_matrix(1, 2)%f(eta1, eta2, transf%params))
      j21 = (transf%j_matrix(2, 1)%f(eta1, eta2, transf%params))
      j22 = (transf%j_matrix(2, 2)%f(eta1, eta2, transf%params))
      ! For debugging:
      !    print *, 'jacobian_2d_analytic: '
      !    print *, j11, j12
      !    print *, j21, j22
      val = j11*j22 - j12*j21
   end function jacobian_2d_analytic

   ! The efficiency of the following function could be improved by just
   ! passing the output array rather than returning values on the stack which
   ! need to be caught by the caller.
   function jacobian_matrix_2d_analytic(transf, eta1, eta2)
      sll_real64, dimension(1:2, 1:2)     :: jacobian_matrix_2d_analytic
      class(sll_t_coordinate_transformation_2d_analytic), intent(inout):: transf
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      sll_real64             :: j11
      sll_real64             :: j12
      sll_real64             :: j21
      sll_real64             :: j22
      j11 = (transf%j_matrix(1, 1)%f(eta1, eta2, transf%params))
      j12 = (transf%j_matrix(1, 2)%f(eta1, eta2, transf%params))
      j21 = (transf%j_matrix(2, 1)%f(eta1, eta2, transf%params))
      j22 = (transf%j_matrix(2, 2)%f(eta1, eta2, transf%params))
      ! For debugging:
      !    print *, 'jacobian_2d_analytic: '
      !    print *, j11, j12
      !    print *, j21, j22
      jacobian_matrix_2d_analytic(1, 1) = j11
      jacobian_matrix_2d_analytic(1, 2) = j12
      jacobian_matrix_2d_analytic(2, 1) = j21
      jacobian_matrix_2d_analytic(2, 2) = j22
   end function jacobian_matrix_2d_analytic

   function inverse_jacobian_matrix_2d_analytic(transf, eta1, eta2)
      sll_real64, dimension(1:2, 1:2)     :: inverse_jacobian_matrix_2d_analytic
      class(sll_t_coordinate_transformation_2d_analytic), intent(inout) :: transf
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      sll_real64             :: inv_j11
      sll_real64             :: inv_j12
      sll_real64             :: inv_j21
      sll_real64             :: inv_j22
      sll_real64             :: r_jacobian ! reciprocal of the jacobian

      r_jacobian = 1.0_f64/transf%jacobian(eta1, eta2)
      inv_j11 = (transf%j_matrix(2, 2)%f(eta1, eta2, transf%params))*r_jacobian
      inv_j12 = -(transf%j_matrix(1, 2)%f(eta1, eta2, transf%params))*r_jacobian
      inv_j21 = -(transf%j_matrix(2, 1)%f(eta1, eta2, transf%params))*r_jacobian
      inv_j22 = (transf%j_matrix(1, 1)%f(eta1, eta2, transf%params))*r_jacobian
      ! For debugging:
      !    print *, 'jacobian_2d_analytic: '
      !    print *, j11, j12
      !    print *, j21, j22
      inverse_jacobian_matrix_2d_analytic(1, 1) = inv_j11
      inverse_jacobian_matrix_2d_analytic(1, 2) = inv_j12
      inverse_jacobian_matrix_2d_analytic(2, 1) = inv_j21
      inverse_jacobian_matrix_2d_analytic(2, 2) = inv_j22
   end function inverse_jacobian_matrix_2d_analytic

   function x1_analytic(transf, eta1, eta2) result(val)
      sll_real64                         :: val
      class(sll_t_coordinate_transformation_2d_analytic) :: transf
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      val = transf%x1_func(eta1, eta2, transf%params)
   end function x1_analytic

   function x2_analytic(transf, eta1, eta2) result(val)
      sll_real64                        :: val
      class(sll_t_coordinate_transformation_2d_analytic) :: transf
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      val = transf%x2_func(eta1, eta2, transf%params)
   end function x2_analytic

   function x1_node_analytic(transf, i, j) result(val)
      class(sll_t_coordinate_transformation_2d_analytic) :: transf
      sll_real64             :: val
      sll_int32, intent(in) :: i
      sll_int32, intent(in) :: j
      sll_real64            :: eta1
      sll_real64            :: eta2

      eta1 = transf%mesh%eta1_node(i, j)
      eta2 = transf%mesh%eta2_node(i, j)
      val = transf%x1_func(eta1, eta2, transf%params)
   end function x1_node_analytic

   function x2_node_analytic(transf, i, j) result(val)
      class(sll_t_coordinate_transformation_2d_analytic) :: transf
      sll_real64             :: val
      sll_int32, intent(in) :: i
      sll_int32, intent(in) :: j
      sll_real64            :: eta1
      sll_real64            :: eta2

      eta1 = transf%mesh%eta1_node(i, j)
      eta2 = transf%mesh%eta2_node(i, j)
      val = transf%x2_func(eta1, eta2, transf%params)
   end function x2_node_analytic

   function x1_cell_analytic(transf, i, j) result(var)
      class(sll_t_coordinate_transformation_2d_analytic) :: transf
      sll_real64            :: var
      sll_int32, intent(in) :: i
      sll_int32, intent(in) :: j
      sll_real64 :: eta1
      sll_real64 :: eta2

      eta1 = transf%mesh%eta1_cell(i, j)
      eta2 = transf%mesh%eta2_cell(i, j)
      var = transf%x1_func(eta1, eta2, transf%params)
   end function x1_cell_analytic

   function x2_cell_analytic(transf, i, j) result(var)
      class(sll_t_coordinate_transformation_2d_analytic) :: transf
      sll_real64            :: var
      sll_int32, intent(in) :: i
      sll_int32, intent(in) :: j
      sll_real64 :: eta1
      sll_real64 :: eta2

      eta1 = transf%mesh%eta1_cell(i, j)
      eta2 = transf%mesh%eta2_cell(i, j)
      var = transf%x2_func(eta1, eta2, transf%params)
   end function x2_cell_analytic

   function jacobian_2d_cell_analytic(transf, i, j) result(val)
      class(sll_t_coordinate_transformation_2d_analytic) :: transf
      sll_real64            :: val
      sll_int32, intent(in) :: i
      sll_int32, intent(in) :: j
      sll_real64 :: eta1
      sll_real64 :: eta2
      sll_real64 :: j11
      sll_real64 :: j12
      sll_real64 :: j21
      sll_real64 :: j22

      eta1 = transf%mesh%eta1_cell(i, j)
      eta2 = transf%mesh%eta2_cell(i, j)
      j11 = (transf%j_matrix(1, 1)%f(eta1, eta2, transf%params))
      j12 = (transf%j_matrix(1, 2)%f(eta1, eta2, transf%params))
      j21 = (transf%j_matrix(2, 1)%f(eta1, eta2, transf%params))
      j22 = (transf%j_matrix(2, 2)%f(eta1, eta2, transf%params))
      ! For debugging:
      !    print *, 'jacobian_2d_analytic: '
      !    print *, j11, j12
      !    print *, j21, j22
      val = j11*j22 - j12*j21
   end function jacobian_2d_cell_analytic

   function jacobian_node_analytic(transf, i, j)
      class(sll_t_coordinate_transformation_2d_analytic)   :: transf
      sll_real64              :: jacobian_node_analytic
      sll_int32, intent(in)   :: i
      sll_int32, intent(in)   :: j
      sll_int32 :: num_pts_1
      sll_int32 :: num_pts_2
      sll_real64 :: eta1
      sll_real64 :: eta2
      sll_real64 :: j11
      sll_real64 :: j12
      sll_real64 :: j21
      sll_real64 :: j22

      num_pts_1 = transf%mesh%num_cells1 + 1
      num_pts_2 = transf%mesh%num_cells2 + 1
      SLL_ASSERT((i .ge. 1) .and. (i .le. num_pts_1))
      SLL_ASSERT((j .ge. 1) .and. (j .le. num_pts_2))

      eta1 = transf%mesh%eta1_cell(i, j)
      eta2 = transf%mesh%eta2_cell(i, j)
      j11 = (transf%j_matrix(1, 1)%f(eta1, eta2, transf%params))
      j12 = (transf%j_matrix(1, 2)%f(eta1, eta2, transf%params))
      j21 = (transf%j_matrix(2, 1)%f(eta1, eta2, transf%params))
      j22 = (transf%j_matrix(2, 2)%f(eta1, eta2, transf%params))
      ! For debugging:
      !    print *, 'jacobian_2d_analytic: '
      !    print *, j11, j12
      !    print *, j21, j22
      jacobian_node_analytic = j11*j22 - j12*j21
   end function jacobian_node_analytic

   subroutine write_to_file_2d_analytic(transf, output_format)
      class(sll_t_coordinate_transformation_2d_analytic), intent(inout) :: transf
      sll_int32, intent(in), optional :: output_format
      sll_int32           :: local_format
      sll_real64, dimension(:, :), allocatable :: x1mesh
      sll_real64, dimension(:, :), allocatable :: x2mesh
      sll_int32  :: i1
      sll_int32  :: i2
      sll_real64 :: eta1
      sll_real64 :: eta2
      sll_int32  :: ierr
      sll_int32  :: file_id
      sll_int32  :: nc_eta1
      sll_int32  :: nc_eta2

      nc_eta1 = transf%mesh%num_cells1
      nc_eta2 = transf%mesh%num_cells2

      if (.not. present(output_format)) then
         local_format = sll_p_io_gnuplot
      else
         local_format = output_format
      end if

      if (.not. transf%written) then
         if (local_format == sll_p_io_xdmf) then
            SLL_ALLOCATE(x1mesh(nc_eta1 + 1, nc_eta2 + 1), ierr)
            SLL_ALLOCATE(x2mesh(nc_eta1 + 1, nc_eta2 + 1), ierr)
            eta1 = transf%mesh%eta1_min
            do i1 = 1, nc_eta1 + 1
               eta2 = transf%mesh%eta2_min
               do i2 = 1, nc_eta2 + 1
                  x1mesh(i1, i2) = x1_node_analytic(transf, i1, i2)
                  x2mesh(i1, i2) = x2_node_analytic(transf, i1, i2)
                  eta2 = eta2 + transf%mesh%delta_eta2
               end do
               eta1 = eta1 + transf%mesh%delta_eta1
            end do

            call sll_o_xdmf_open(trim(transf%label)//".xmf", transf%label, &
                                 nc_eta1 + 1, nc_eta2 + 1, file_id, ierr)
            call sll_o_xdmf_write_array(transf%label, x1mesh, "x1", ierr)
            call sll_o_xdmf_write_array(transf%label, x2mesh, "x2", ierr)
            call sll_s_xdmf_close(file_id, ierr)

         else if (local_format == sll_p_io_gnuplot) then

            SLL_ALLOCATE(x1mesh(nc_eta1 + 1, nc_eta2 + 1), ierr)
            SLL_ALLOCATE(x2mesh(nc_eta1 + 1, nc_eta2 + 1), ierr)

            do i1 = 1, nc_eta1 + 1
               do i2 = 1, nc_eta2 + 1
                  x1mesh(i1, i2) = x1_node_analytic(transf, i1, i2)
                  x2mesh(i1, i2) = x2_node_analytic(transf, i1, i2)
               end do
            end do

            call sll_o_gnuplot_2d(nc_eta1 + 1, nc_eta2 + 1, x1mesh, x2mesh, &
                                  trim(transf%label), ierr)

         else if (local_format == sll_p_io_mtv) then

            SLL_ALLOCATE(x1mesh(nc_eta1 + 1, nc_eta2 + 1), ierr)
            SLL_ALLOCATE(x2mesh(nc_eta1 + 1, nc_eta2 + 1), ierr)

            do i1 = 1, nc_eta1 + 1
               do i2 = 1, nc_eta2 + 1
                  x1mesh(i1, i2) = x1_node_analytic(transf, i1, i2)
                  x2mesh(i1, i2) = x2_node_analytic(transf, i1, i2)
               end do
            end do

            call sll_o_plotmtv_write(nc_eta1 + 1, nc_eta2 + 1, &
                                     x1mesh, x2mesh, trim(transf%label), ierr)

         else
            print *, 'Not recognized format to write this mesh'
            stop
         end if
      else
         print *, ' Warning, you have already written the mesh '
      end if
      transf%written = .true.
      if (allocated(x1mesh)) deallocate (x1mesh)
      if (allocated(x2mesh)) deallocate (x2mesh)

   end subroutine write_to_file_2d_analytic

   subroutine read_from_file_2d_analytic(transf, filename)
      class(sll_t_coordinate_transformation_2d_analytic), intent(inout) :: transf
      character(len=*), intent(in) :: filename
      print *, filename
      print *, 'read_from_file_2d_analytic: not yet implemented'
      !call sll_o_display(transf%mesh)
      call transf%mesh%display()
      ! here we could put a case select to choose which analytic transformation
      ! we would like to use.
   end subroutine read_from_file_2d_analytic

   !**************************************************************************
   !
   !        Functions for the discrete general transformation
   !
   !**************************************************************************

   function get_cartesian_mesh_discrete(transf) result(res)
      class(sll_t_coordinate_transformation_2d_discrete), intent(in) :: transf
      class(sll_t_cartesian_mesh_2d), pointer :: res
      res => transf%mesh
   end function get_cartesian_mesh_discrete

   function x1_node_discrete(transf, i, j) result(val)
      class(sll_t_coordinate_transformation_2d_discrete) :: transf
      sll_real64             :: val
      sll_int32, intent(in) :: i
      sll_int32, intent(in) :: j
      val = transf%x1_node(i, j)
   end function x1_node_discrete

   function x2_node_discrete(transf, i, j) result(val)
      class(sll_t_coordinate_transformation_2d_discrete) :: transf
      sll_real64             :: val
      sll_int32, intent(in) :: i
      sll_int32, intent(in) :: j
      val = transf%x2_node(i, j)
   end function x2_node_discrete

   function x1_cell_discrete(transf, i, j) result(var)
      class(sll_t_coordinate_transformation_2d_discrete) :: transf
      sll_real64                         :: var
      sll_int32, intent(in)              :: i
      sll_int32, intent(in)              :: j
      var = transf%x1_cell(i, j)
   end function x1_cell_discrete

   function x2_cell_discrete(transf, i, j) result(var)
      class(sll_t_coordinate_transformation_2d_discrete) :: transf
      sll_real64                         :: var
      sll_int32, intent(in)              :: i
      sll_int32, intent(in)              :: j
      var = transf%x2_cell(i, j)
   end function x2_cell_discrete

   function x1_discrete(transf, eta1, eta2) result(val)
      class(sll_t_coordinate_transformation_2d_discrete) :: transf
      sll_real64             :: val
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      val = transf%x1_interp%interpolate_from_interpolant_value(eta1, eta2)
   end function x1_discrete

   function x2_discrete(transf, eta1, eta2) result(val)
      class(sll_t_coordinate_transformation_2d_discrete) :: transf
      sll_real64             :: val
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      val = transf%x2_interp%interpolate_from_interpolant_value(eta1, eta2)
   end function x2_discrete

   function jacobian_2d_discrete(transf, eta1, eta2) result(jac)
      class(sll_t_coordinate_transformation_2d_discrete) :: transf
      sll_real64             :: jac
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      sll_real64             :: j11
      sll_real64             :: j12
      sll_real64             :: j21
      sll_real64             :: j22
      j11 = transf%x1_interp%interpolate_from_interpolant_derivative_eta1(eta1, eta2)
      j12 = transf%x1_interp%interpolate_from_interpolant_derivative_eta2(eta1, eta2)
      j21 = transf%x2_interp%interpolate_from_interpolant_derivative_eta1(eta1, eta2)
      j22 = transf%x2_interp%interpolate_from_interpolant_derivative_eta2(eta1, eta2)
      ! For debugging:
      !    print *, 'jacobian_2D_discrete: '
      !    print *, j11, j12
      !    print *, j21, j22
      jac = j11*j22 - j12*j21
   end function jacobian_2d_discrete

   function jacobian_2d_cell_discrete(transf, i, j) result(var)
      class(sll_t_coordinate_transformation_2d_discrete) :: transf
      sll_real64                         :: var
      sll_int32, intent(in)              :: i
      sll_int32, intent(in)              :: j
      var = transf%jacobians_c(i, j)
   end function jacobian_2d_cell_discrete

   function jacobian_matrix_2d_discrete(transf, eta1, eta2)
      class(sll_t_coordinate_transformation_2d_discrete), intent(inout) :: transf
      sll_real64, dimension(1:2, 1:2)     :: jacobian_matrix_2d_discrete
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      sll_real64             :: j11
      sll_real64             :: j12
      sll_real64             :: j21
      sll_real64             :: j22
      j11 = transf%x1_interp%interpolate_from_interpolant_derivative_eta1(eta1, eta2)
      j12 = transf%x1_interp%interpolate_from_interpolant_derivative_eta2(eta1, eta2)
      j21 = transf%x2_interp%interpolate_from_interpolant_derivative_eta1(eta1, eta2)
      j22 = transf%x2_interp%interpolate_from_interpolant_derivative_eta2(eta1, eta2)
      ! For debugging:
      !    print *, 'jacobian_2D_discrete: '
      !    print *, j11, j12
      !    print *, j21, j22
      jacobian_matrix_2d_discrete(1, 1) = j11
      jacobian_matrix_2d_discrete(1, 2) = j12
      jacobian_matrix_2d_discrete(2, 1) = j21
      jacobian_matrix_2d_discrete(2, 2) = j22
   end function jacobian_matrix_2d_discrete

   function inverse_jacobian_matrix_2d_discrete(transf, eta1, eta2)
      class(sll_t_coordinate_transformation_2d_discrete), intent(inout) :: transf
      sll_real64, dimension(1:2, 1:2)     :: inverse_jacobian_matrix_2d_discrete
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      sll_real64             :: inv_j11
      sll_real64             :: inv_j12
      sll_real64             :: inv_j21
      sll_real64             :: inv_j22
      sll_real64             :: r_jac ! reciprocal of the jacobian
      r_jac = 1.0_f64/transf%jacobian(eta1, eta2)
      inv_j11 = transf%x1_interp%interpolate_from_interpolant_derivative_eta1(eta1, eta2)
      inv_j12 = transf%x1_interp%interpolate_from_interpolant_derivative_eta2(eta1, eta2)
      inv_j21 = transf%x2_interp%interpolate_from_interpolant_derivative_eta1(eta1, eta2)
      inv_j22 = transf%x2_interp%interpolate_from_interpolant_derivative_eta2(eta1, eta2)
      ! For debugging:
      !    print *, 'jacobian_2D_discrete: '
      !    print *, j11, j12
      !    print *, j21, j22
      inverse_jacobian_matrix_2d_discrete(1, 1) = inv_j22*r_jac
      inverse_jacobian_matrix_2d_discrete(1, 2) = -inv_j12*r_jac
      inverse_jacobian_matrix_2d_discrete(2, 1) = -inv_j21*r_jac
      inverse_jacobian_matrix_2d_discrete(2, 2) = inv_j11*r_jac
   end function inverse_jacobian_matrix_2d_discrete

   !> Create a new coordinate transformation object
   function sll_f_new_coordinate_transformation_2d_discrete( &
      mesh_2d, &
      label, &
      x1_interpolator, &
      x2_interpolator, &
      jacobians_n_interpolator, &
      x1_node, &
      x2_node, &
      jacobians_node, &
      x1_cell, &
      x2_cell, &
      jacobians_cell)

      ! INPUT VARIABLES
      type(sll_t_cartesian_mesh_2d), pointer    :: mesh_2d
      character(len=*), intent(in) :: label

      class(sll_c_interpolator_2d), target  :: x1_interpolator
      class(sll_c_interpolator_2d), target  :: x2_interpolator
      class(sll_c_interpolator_2d), target  :: jacobians_n_interpolator
      sll_real64, dimension(:, :), intent(in), optional :: x1_node
      sll_real64, dimension(:, :), intent(in), optional :: x2_node
      sll_real64, dimension(:, :), intent(in), optional :: jacobians_node
      sll_real64, dimension(:, :), intent(in), optional :: x1_cell
      sll_real64, dimension(:, :), intent(in), optional :: x2_cell
      sll_real64, dimension(:, :), intent(in), optional :: jacobians_cell

      ! LOCAL VARIABLES
      type(sll_t_coordinate_transformation_2d_discrete), pointer :: &
         sll_f_new_coordinate_transformation_2d_discrete
      sll_int32 :: ierr

      SLL_ALLOCATE(sll_f_new_coordinate_transformation_2d_discrete, ierr)
      call sll_s_coordinate_transformation_2d_discrete_init( &
         sll_f_new_coordinate_transformation_2d_discrete, &
         mesh_2d, &
         label, &
         x1_interpolator, &
         x2_interpolator, &
         jacobians_n_interpolator, &
         x1_node, &
         x2_node, &
         jacobians_node, &
         x1_cell, &
         x2_cell, &
         jacobians_cell)
   end function sll_f_new_coordinate_transformation_2d_discrete

   subroutine sll_s_coordinate_transformation_2d_discrete_init( &
      transf, &
      mesh_2d, &
      label, &
      x1_interpolator, &
      x2_interpolator, &
      jacobians_n_interpolator, &
      x1_node, &
      x2_node, &
      jacobians_node, &
      x1_cell, &
      x2_cell, &
      jacobians_cell)

      class(sll_t_coordinate_transformation_2d_discrete) :: transf
      type(sll_t_cartesian_mesh_2d), target              :: mesh_2d
      character(len=*), intent(in)                       :: label

      class(sll_c_interpolator_2d), target  :: x1_interpolator
      class(sll_c_interpolator_2d), target  :: x2_interpolator
      class(sll_c_interpolator_2d), target  :: jacobians_n_interpolator
      sll_real64, dimension(:, :), intent(in), optional :: x1_node
      sll_real64, dimension(:, :), intent(in), optional :: x2_node
      sll_real64, dimension(:, :), intent(in), optional :: jacobians_node
      sll_real64, dimension(:, :), intent(in), optional :: jacobians_cell
      sll_real64, dimension(:, :), intent(in), optional :: x1_cell
      sll_real64, dimension(:, :), intent(in), optional :: x2_cell

      sll_real64 :: eta_1
      sll_real64 :: eta_2
      sll_real64 :: eta_1_min
      sll_real64 :: eta_2_min
      sll_real64 :: delta_eta_1
      sll_real64 :: delta_eta_2
      sll_real64 :: jacobian_val
      sll_int32  :: i
      sll_int32  :: j
      sll_int32  :: ierr
      sll_int32  :: npts1
      sll_int32  :: npts2
      logical    :: x1n
      logical    :: x2n
      logical    :: x1c
      logical    :: x2c
      logical    :: jc
      logical    :: jn

      transf%mesh => mesh_2d
      transf%label = trim(label)
      x1n = present(x1_node)
      x2n = present(x2_node)
      x1c = present(x1_cell)
      x2c = present(x2_cell)
      jc = present(jacobians_cell)
      jn = present(jacobians_node)
      npts1 = mesh_2d%num_cells1 + 1
      npts2 = mesh_2d%num_cells2 + 1

      ! Check argument consistency
      ! DISCRETE_MAPs require only some of the parameters. If the mapping is
      ! defined from the nodes of the logical (eta1, eta2) mesh to the nodes
      ! of the physical mesh (x1,x2), then either:
      ! - the node arrays are required: jacobians_node, x1_node and x2_node. Or
      ! - the x1_interpolator and x2_interpolator must contain already the
      !   coefficient information that would permit the calculation of the
      !   x1 and x2 points.
      !
      ! If the transformation is done on the points at the center of the cells
      ! then these parameters are also required:
      ! jacobians_cell, x1_cell, x2_cell.
      ! node and cell values are not mutually exclusive, thus all 6 parameters
      ! can be provided in the discrete case. It is up to the caller to make
      ! sure that the data set is consistent.

      ! 1. Check that the discrete representation of the transformation is
      !    consistent with the size of the 2D array.

      if ((x1n .and. (.not. x2n)) .or. &
          ((.not. x1n) .and. x2n)) then

         print *, 'ERROR, sll_s_coordinate_transformation_2d_discrete_init():', &
            'for the moment, this function does not support specifying ', &
            'transformation only with one of the node arrays x1_node or ', &
            'x2_node. Either pass both, or none, but with the ', &
            'corresponding interpolators having their coefficients ', &
            'already set.'
         STOP
         call jacobians_n_interpolator%delete() !PN added to remove the warning
      end if

      if (x1n) then
         if (size(x1_node, 1) .lt. npts1) then
            print *, 'ERROR, sll_s_coordinate_transformation_2d_discrete_init()', &
               ' the size of the x1_node arrays is ', &
               'inconsistent with the number of points declared, ', &
               'in the logical mesh.'
            STOP
         end if
      end if

      if (x2n) then
         if (size(x1_node, 2) .lt. npts2) then
            print *, 'ERROR, sll_s_coordinate_transformation_2d_discrete_init()', &
               ' the size of the x2_node arrays is ', &
               'inconsistent with the number of points declared, ', &
               'in the logical mesh.'
            STOP
         end if
      end if

      if (jn .eqv. .true.) then
         if ( &
            (size(jacobians_node, 1) .lt. npts1 - 1) .or. &
            (size(jacobians_node, 2) .lt. npts2 - 1)) then
            print *, 'ERROR, sll_s_coordinate_transformation_2d_discrete_init()', &
               ': the size of the jacobians_node array is ', &
               'inconsistent with the number of points declared, ', &
               'npts1 or npts2.'
            STOP
         end if
      end if

      if (jc .eqv. .true.) then
         if ( &
            (size(jacobians_cell, 1) .lt. npts1 - 1) .or. &
            (size(jacobians_cell, 2) .lt. npts2 - 1)) then
            print *, 'ERROR, sll_s_coordinate_transformation_2d_discrete_init()', &
               ': the size of the jacobians_cell arrays is ', &
               'inconsistent with the number of points declared, ', &
               'npts1 or npts2.'
            STOP
         end if
      end if

      transf%x1_interp => x1_interpolator
      transf%x2_interp => x2_interpolator

      ! Allocate the arrays for precomputed jacobians.
      SLL_ALLOCATE(transf%jacobians_n(npts1, npts2), ierr)
      SLL_ALLOCATE(transf%jacobians_c(npts1 - 1, npts2 - 1), ierr)

      ! Allocation for x1 and x2 at nodes
      SLL_ALLOCATE(transf%x1_node(npts1, npts2), ierr)
      SLL_ALLOCATE(transf%x2_node(npts1, npts2), ierr)

      ! Allocation for x1 and x2 at cells
      SLL_ALLOCATE(transf%x1_cell(npts1 - 1, npts2 - 1), ierr)
      SLL_ALLOCATE(transf%x2_cell(npts1 - 1, npts2 - 1), ierr)

      ! initialize the local arrays. Note that since the map has its
      ! own copies, it owns this information locally and will destroy
      ! this information when the object is deleted. The caller is
      ! thus responsible for deallocating the arrays that were passed as
      ! arguments.

      eta_1_min = mesh_2d%eta1_min
      eta_2_min = mesh_2d%eta2_min
      delta_eta_1 = mesh_2d%delta_eta1
      delta_eta_2 = mesh_2d%delta_eta2

      if (x1n .and. x2n) then
         do j = 1, npts2
            do i = 1, npts1
               transf%x1_node(i, j) = x1_node(i, j)
               transf%x2_node(i, j) = x2_node(i, j)
            end do
         end do
      else
         if (x1_interpolator%coefficients_are_set() .eqv. .false.) then
            print *, 'ERROR, sll_s_coordinate_transformation_2d_discrete_init()', &
               ': the x1_node array was not passed and the corresponding ', &
               'interpolator has no initialized coefficients. Exiting...'
            STOP
         end if
         if (x2_interpolator%coefficients_are_set() .eqv. .false.) then
            print *, 'ERROR, sll_s_coordinate_transformation_2d_discrete_init()', &
               ': the x2_node array was not passed and the corresponding ', &
               'interpolator has no initialized coefficients. Exiting...'
            STOP
         end if
         ! now initialize the arrays starting from the interpolator information
         ! and the logical mesh information.
         do j = 0, npts2 - 1
            eta_2 = eta_2_min + real(j, f64)*delta_eta_2
            do i = 0, npts1 - 1
               eta_1 = eta_1_min + real(i, f64)*delta_eta_1
               transf%x1_node(i + 1, j + 1) = &
                  x1_interpolator%interpolate_from_interpolant_value(eta_1, eta_2)
               transf%x2_node(i + 1, j + 1) = &
                  x2_interpolator%interpolate_from_interpolant_value(eta_1, eta_2)
            end do
         end do
      end if

      ! Compute the spline coefficients
      if (x1n .and. (x1_interpolator%coefficients_are_set() .eqv. .false.)) then
         call x1_interpolator%compute_interpolants(transf%x1_node)
      end if
      if (x2n .and. (x2_interpolator%coefficients_are_set() .eqv. .false.)) then
         call x2_interpolator%compute_interpolants(transf%x2_node)
      end if

      ! The splines contain all the information to compute the
      ! jacobians everywhere; however, here we explore assigning
      ! the jacobians-at-the-nodes array with the values provided
      ! by the user if available. If there are discrepancies between
      ! the user-provided values and the predictions from the splines,
      ! then this may itself be a way to look for errors.
      !
      ! Copy the values of the jacobians at the nodes if user given:
      if (jn .eqv. .true.) then
         do j = 1, npts2
            do i = 1, npts1
               transf%jacobians_n(i, j) = jacobians_node(i, j)
            end do
         end do
      else
         ! Fill the jacobian values at the nodes calculated from the splines
         do j = 0, npts2 - 1
            eta_2 = eta_2_min + real(j, f64)*delta_eta_2
            do i = 0, npts1 - 1
               eta_1 = eta_1_min + real(i, f64)*delta_eta_1
               jacobian_val = transf%jacobian(eta_1, eta_2)
               transf%jacobians_n(i + 1, j + 1) = jacobian_val
            end do
         end do
      end if

      ! copy the cell-based transformation arrays if available
      if ((x1c .and. x2c) .eqv. .true.) then
         SLL_ALLOCATE(transf%x1_cell(npts1 - 1, npts2 - 1), ierr)
         SLL_ALLOCATE(transf%x2_cell(npts1 - 1, npts2 - 1), ierr)
         do j = 1, npts2 - 1
            do i = 1, npts1 - 1
               transf%x1_cell(i, j) = x1_cell(i, j)
               transf%x2_cell(i, j) = x2_cell(i, j)
            end do
         end do
      end if
      ! copy the cell-based jacobians if available
      if (jc .eqv. .true.) then
         do j = 1, npts2 - 1
            do i = 1, npts1 - 1
               transf%jacobians_c(i, j) = jacobians_cell(i, j)
            end do
         end do
      else ! if cell-based jacobians are not available, compute them.
         ! Fill the values at the mid-point of the cells
         do j = 0, npts2 - 2
            eta_2 = eta_2_min + delta_eta_2*(real(j, f64) + 0.5_f64)
            do i = 0, npts1 - 2
               ! it is very bad practice to invoke the mesh methods while
               ! we are not even done initializing the mesh object...
               eta_1 = eta_1_min + delta_eta_1*(real(i, f64) + 0.5_f64)
               transf%x1_cell(i + 1, j + 1) = transf%x1(eta_1, eta_2)
               transf%x2_cell(i + 1, j + 1) = transf%x2(eta_1, eta_2)
               transf%jacobians_c(i + 1, j + 1) = transf%jacobian(eta_1, eta_2)
            end do
         end do
      end if
   end subroutine sll_s_coordinate_transformation_2d_discrete_init

   function transf_2d_jacobian_node_discrete(transf, i, j)
      class(sll_t_coordinate_transformation_2d_discrete)   :: transf
      sll_real64              :: transf_2d_jacobian_node_discrete
      sll_int32, intent(in)   :: i
      sll_int32, intent(in)   :: j
      sll_int32 :: num_pts_1
      sll_int32 :: num_pts_2

      num_pts_1 = transf%mesh%num_cells1 + 1
      num_pts_2 = transf%mesh%num_cells2 + 1
      SLL_ASSERT((i .ge. 1) .and. (i .le. num_pts_1))
      SLL_ASSERT((j .ge. 1) .and. (j .le. num_pts_2))
      transf_2d_jacobian_node_discrete = transf%jacobians_n(i, j)
   end function transf_2d_jacobian_node_discrete

   subroutine write_to_file_2d_discrete(transf, output_format)
      class(sll_t_coordinate_transformation_2d_discrete), intent(inout) :: transf
      sll_int32, intent(in), optional :: output_format
      sll_int32           :: local_format
      sll_real64, dimension(:, :), pointer :: x1mesh
      sll_real64, dimension(:, :), pointer :: x2mesh
      sll_int32  :: i1
      sll_int32  :: i2
      sll_real64 :: eta1
      sll_real64 :: eta2
      sll_int32  :: ierr
      sll_int32  :: file_id
      sll_int32  :: npts_eta1
      sll_int32  :: npts_eta2
      sll_real64 :: eta1_min
      sll_real64 :: eta2_min
      sll_real64 :: delta_eta1
      sll_real64 :: delta_eta2

      npts_eta1 = transf%mesh%num_cells1 + 1
      npts_eta2 = transf%mesh%num_cells2 + 1
      eta1_min = transf%mesh%eta1_min
      eta2_min = transf%mesh%eta1_min
      delta_eta1 = transf%mesh%delta_eta1
      delta_eta2 = transf%mesh%delta_eta2

      if (.not. present(output_format)) then
         local_format = sll_p_io_gnuplot
      else
         local_format = output_format
      end if

      if (.not. transf%written) then

         if (local_format == sll_p_io_xdmf) then
            SLL_ALLOCATE(x1mesh(npts_eta1, npts_eta2), ierr)
            SLL_ALLOCATE(x2mesh(npts_eta1, npts_eta2), ierr)
            eta1 = eta1_min
            do i1 = 1, npts_eta1
               eta2 = eta2_min
               do i2 = 1, npts_eta2
                  x1mesh(i1, i2) = x1_node_discrete(transf, i1, i2)
                  x2mesh(i1, i2) = x2_node_discrete(transf, i1, i2)
                  eta2 = eta2 + delta_eta2
               end do
               eta1 = eta1 + delta_eta1
            end do

            call sll_o_xdmf_open(trim(transf%label)//".xmf", transf%label, &
                                 npts_eta1, npts_eta2, file_id, ierr)
            call sll_o_xdmf_write_array(transf%label, x1mesh, "x1", ierr)
            call sll_o_xdmf_write_array(transf%label, x2mesh, "x2", ierr)
            call sll_s_xdmf_close(file_id, ierr)

         else if (local_format == sll_p_io_gnuplot) then

            SLL_ALLOCATE(x1mesh(npts_eta1, npts_eta2), ierr)
            SLL_ALLOCATE(x2mesh(npts_eta1, npts_eta2), ierr)

            do i1 = 1, npts_eta1
               do i2 = 1, npts_eta2
                  x1mesh(i1, i2) = x1_node_discrete(transf, i1, i2)
                  x2mesh(i1, i2) = x2_node_discrete(transf, i1, i2)
               end do
            end do

            call sll_o_gnuplot_2d(npts_eta1, npts_eta2, x1mesh, x2mesh, &
                                  trim(transf%label), ierr)

         else if (local_format == sll_p_io_mtv) then

            SLL_ALLOCATE(x1mesh(npts_eta1, npts_eta2), ierr)
            SLL_ALLOCATE(x2mesh(npts_eta1, npts_eta2), ierr)

            do i1 = 1, npts_eta1
               do i2 = 1, npts_eta2
                  x1mesh(i1, i2) = x1_node_discrete(transf, i1, i2)
                  x2mesh(i1, i2) = x2_node_discrete(transf, i1, i2)
               end do
            end do

            call sll_o_plotmtv_write(npts_eta1, npts_eta2, &
                                     x1mesh, x2mesh, trim(transf%label), ierr)

         else

            print *, 'Not recognized format to write this mesh'
            stop

         end if
      else
         print *, ' Warning, you have already written the mesh '
      end if

      transf%written = .true.
   end subroutine

   subroutine delete_transformation_2d_discrete(transf)
      class(sll_t_coordinate_transformation_2d_discrete), intent(inout) :: transf
      transf%label = ""
      transf%written = .false.
      nullify (transf%x1_node)
      nullify (transf%x2_node)
      nullify (transf%x1_cell)
      nullify (transf%x2_cell)
      nullify (transf%jacobians_n)
      nullify (transf%jacobians_c)

      !call delete( transf%x1_interp)
      !call delete( transf%x2_interp)
      call sll_o_delete(transf%mesh)
      ! Fix: there is a dependency problem where these pointers are not recognized
      ! during the linking step. A similar nullification of an abstract class
      ! pointer is carried out in the fields_2d type without problems.
!    transf%x1_interp => null() this gives a different message.
!!$    nullify( transf%x1_interp )
!!$    nullify( transf%x2_interp )
   end subroutine delete_transformation_2d_discrete

   ! What do we need to initialize fully a discrete coordinate transformation?
   ! - logical mesh
   ! - label
   ! - array with x1 node positions
   ! - array with x2 node positions
   ! - array with x1 at cell-center positions
   ! - array with x2 at cell-center positions
   ! - array with jacobians at nodes
   ! - array with jacobians at cell-centers
   ! - interpolator 2d for x1
   ! - interpolator 2d for x2
   ! - the file used to initialize the transformation should allow us to
   !   initialize all this data. This routine has special rights in that it
   !   is allowed to allocate and initialize a logical mesh and the two
   !   interpolators.

   ! - Issues to decide:
   ! - Will there be a single file format? Or multiple file formats?
   !   The transformation can be specified by two 2D arrays of points, or by
   !   the spline coefficients...
   ! - The BC information is not inside the files we are currently considering,
   !   so this should be included.

   subroutine read_from_file_2d_discrete(transf, filename)
      class(sll_t_coordinate_transformation_2d_discrete), intent(inout) :: transf
      character(len=*), intent(in) :: filename
      intrinsic :: trim
      character(len=256) :: filename_local
      sll_int32 :: IO_stat
      sll_int32 :: input_file_id
      sll_int32 :: ierr
      sll_int32 :: spline_deg1
      sll_int32 :: spline_deg2
      sll_int32 :: num_pts1
      sll_int32 :: num_pts2
      sll_int32 :: is_rational
      character(len=256) :: label
      sll_real64, dimension(:), allocatable :: knots1
      sll_real64, dimension(:), allocatable :: knots2
      sll_real64, dimension(:), allocatable :: control_pts1
      sll_real64, dimension(:), allocatable :: control_pts2
      sll_real64, dimension(:), allocatable :: weights
      sll_real64, dimension(:, :), allocatable :: control_pts1_2d
      sll_real64, dimension(:, :), allocatable :: control_pts2_2d
      sll_real64, dimension(:, :), allocatable :: weights_2d
      sll_real64 :: eta1_min
      sll_real64 :: eta1_max
      sll_real64 :: eta2_min
      sll_real64 :: eta2_max
      sll_int32  :: bc_left
      sll_int32  :: bc_right
      sll_int32  :: bc_bottom
      sll_int32  :: bc_top
!    sll_real64, dimension(:,:), allocatable :: nodes1
!    sll_real64, dimension(:,:), allocatable :: nodes2
      sll_int32  :: number_cells1, number_cells2
      sll_int32 :: sz_knots1, sz_knots2
      type(sll_t_cartesian_mesh_2d), pointer      :: mesh_2d

      namelist /transf_label/ label
      namelist /degree/ spline_deg1, spline_deg2
      namelist /shape/ num_pts1, num_pts2 ! it is not the number of points but the number of coeff sdpline in each direction !!
      namelist /rational/ is_rational
      namelist /knots_1/ knots1
      namelist /knots_2/ knots2
      namelist /control_points/ control_pts1, control_pts2
      namelist /pt_weights/ weights
      namelist /logical_mesh_2d/ number_cells1, number_cells2

      if (len(filename) >= 256) then
         print *, 'ERROR, read_coefficients_from_file => ', &
            'read_from_file_discrete():', &
            'filenames longer than 256 characters are not allowed.'
         STOP
      end if
      filename_local = trim(filename)

      ! get a new identifier for the file.
      call sll_s_new_file_id(input_file_id, ierr)
      if (ierr .ne. 0) then
         print *, 'ERROR while trying to obtain an unique identifier for file ', &
            filename, '. Called from read_coeffs_ad2d().'
         stop
      end if
      open (unit=input_file_id, file=filename_local, STATUS="OLD", IOStat=IO_stat)
      if (IO_Stat .ne. 0) then
         print *, 'ERROR while opening file ', filename, &
            '. Called from read_coeffs_ad2d().'
         stop
      end if

      read (input_file_id, transf_label)
      read (input_file_id, degree)
      read (input_file_id, shape)
      read (input_file_id, rational)
      SLL_ALLOCATE(knots1(num_pts1 + spline_deg1 + 1), ierr)
      SLL_ALLOCATE(knots2(num_pts2 + spline_deg2 + 1), ierr)
      read (input_file_id, knots_1)
      read (input_file_id, knots_2)
      SLL_ALLOCATE(control_pts1(num_pts1*num_pts2), ierr)
      SLL_ALLOCATE(control_pts2(num_pts1*num_pts2), ierr)
      SLL_ALLOCATE(weights(num_pts1*num_pts2), ierr)
      SLL_ALLOCATE(control_pts1_2d(num_pts1, num_pts2), ierr)
      SLL_ALLOCATE(control_pts2_2d(num_pts1, num_pts2), ierr)
      SLL_ALLOCATE(weights_2d(num_pts1, num_pts2), ierr)

      read (input_file_id, control_points)
      control_pts1_2d = reshape(control_pts1, (/num_pts1, num_pts2/))
      control_pts2_2d = reshape(control_pts2, (/num_pts1, num_pts2/))
      read (input_file_id, pt_weights)
      weights_2d = reshape(weights, (/num_pts1, num_pts2/))
      read (input_file_id, logical_mesh_2d)
      close (input_file_id)

      eta1_min = knots1(1)
      eta2_min = knots2(1)
      eta1_max = knots1(num_pts1 + spline_deg1 + 1)
      eta2_max = knots2(num_pts2 + spline_deg2 + 1)

      ! for the moment we put the boundary condition like a dirichlet
      ! boundary condition
      ! but we must modified this part <-- this means that this info must
      ! come within the input file: ECG

      bc_left = sll_p_dirichlet
      bc_right = sll_p_dirichlet
      bc_bottom = sll_p_dirichlet
      bc_top = sll_p_dirichlet

      sz_knots1 = size(knots1)
      sz_knots2 = size(knots2)

      ! Initialization of the interpolator spline 2D in x
      ! ACHTUNG we have not delete it   <--- What???:ECG
      stop "This case in no longer supported"
   end subroutine read_from_file_2d_discrete

end module sll_m_coordinate_transformations_2d
