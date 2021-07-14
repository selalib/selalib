!> @ingroup mesh
!> @brief
!> Module interfaces for coordinate transformation
!> @details
!> xi1, xi2, xi3 are the logical coordinates which are transformed
!> to the physical coordinates x1, x2, x3
!> The Jacobian matrix is the derivate matrix of the transformation function
!> the determinant of the Jacobian is needed for the transformation theorem
!> the logical paramters are in the interval [0,1]
!> @author
!> Benedikt Perse

module sll_m_mapping_3d
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
       sll_p_pi

  use sll_m_3d_coordinate_transformations

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_linear_operator_kron, only : &
       sll_t_linear_operator_kron

  use sll_m_linear_solver_kron, only : &
       sll_t_linear_solver_kron

  use sll_m_linear_solver_mgmres, only : &
       sll_t_linear_solver_mgmres

  use sll_m_splines_pp, only :&
       sll_t_spline_pp_1d, &
       sll_t_spline_pp_3d, &
       sll_f_spline_pp_horner_1d, &
       sll_f_spline_pp_horner_3d, &
       sll_f_spline_pp_horner_3d_d1, &
       sll_f_spline_pp_horner_3d_d2, &
       sll_f_spline_pp_horner_3d_d3, &
       sll_s_spline_pp_init_3d, &
       sll_s_spline_pp_b_to_pp_3d_clamped_2full

  implicit none

  public :: &
       sll_t_mapping_3d

  private
  !> abstract interface for mapping functions
  abstract interface
     function sll_i_eval_function( xi, params ) result(res)
       use sll_m_working_precision
       sll_real64, intent(in) :: xi(3) !< logical coordinates
       sll_real64, intent(in) :: params(:) !< transformation parameters
       sll_real64             :: res
     end function sll_i_eval_function
  end interface

  type matrix_element
     procedure(sll_i_eval_function), pointer, nopass :: f
  end type matrix_element
  !> type collecting functions for analytical coordinate mapping
  type :: sll_t_mapping_3d

     type(matrix_element), dimension(:,:), pointer         :: j_matrix !< jacobi matrix
     procedure(sll_i_eval_function),       pointer, nopass :: jacob !< determinant of the jacobi matrix
     procedure(sll_i_eval_function),       pointer, nopass :: x1_func !< coordinate transformation function
     procedure(sll_i_eval_function),       pointer, nopass :: x2_func !< coordinate transformation function
     procedure(sll_i_eval_function),       pointer, nopass :: x3_func !< coordinate transformation function
     procedure(sll_i_eval_function),       pointer, nopass :: xi1_func !< coordinate transformation function
     procedure(sll_i_eval_function),       pointer, nopass :: xi2_func !< coordinate transformation function
     procedure(sll_i_eval_function),       pointer, nopass :: xi3_func !< coordinate transformation function
     sll_real64,           dimension(:),   pointer         :: params !< transformation parameters
     !sll_real64 :: cylinder_params(3) !< parameter for cylindrical mapping
     logical                                               :: flag = .false. !< flag is true if entries of the Jacobian matrix and its determinant are given analytically
     logical                                               :: flag2d = .false.  !< logical flag for mappings mixing the first two logical coordiantes
     logical                                               :: flag3d = .false.  !< logical flag for mappings mixing all three logical coordinates
     logical                                               :: singular = .false. !< logical for singular mapping
     logical                                               :: inverse = .false.  !< logical for analytical inverse of the chosen mapping

     type(sll_t_spline_pp_3d) :: spline_3d !< 3d pp-spline

     sll_real64, allocatable :: x1_func_dofs_pp(:,:) !< coefficient for spline mapping
     sll_real64, allocatable :: x2_func_dofs_pp(:,:) !< coefficient for spline mapping
     sll_real64, allocatable :: x3_func_dofs_pp(:,:) !< coefficient for spline mapping

     sll_int32 :: n_cells(3) !< number of cells (and grid points)
     sll_int32 :: deg(3)     !< spline deg
     sll_real64 :: volume = 1._f64 !< volume of physical domain
     sll_real64 :: Lx(3) = 1._f64 !< length of physical domain

   contains

     procedure :: get_x1 !< Calculate x1(xi1,xi2,xi3)

     procedure :: get_x2 !< Calculate x2(xi1,xi2,xi3)

     procedure :: get_x3 !< Calculate x3(xi1,xi2,xi3)

     procedure :: get_x  !< Calculate x(xi1,xi2,xi3)

     procedure :: get_xi !< Calculate xi(x1,x2,x3)

     procedure :: jacobian !< Calculate the determinant of the jacobian matrix

     procedure :: jacobian_matrix !< Calculate the jacobian matrix of the mapping

     procedure :: jacobian_matrix_transposed !< Calculate the transposed jacobian matrix of the mapping

     procedure :: jacobian_matrix_inverse !< Calculate the inverse of the jacobian matrix of the mapping

     procedure :: jacobian_matrix_inverse_transposed !< Calculate transposed inverse of the jacobian matrix of the mapping

     procedure :: metric !< Calculate the metric of the mapping

     procedure :: metric_inverse !< Calculate the inverse metric of the mapping

     procedure :: metric_single !< Calculate one entry of the metric of the mapping

     procedure :: metric_inverse_single !< Calculate one entry of the inverse metric of the mapping

     procedure :: metric_single_jacobian !< Calculate one entry of the metric of the mapping divided by the jacobian

     procedure :: metric_inverse_single_jacobian !< Calculate one entry of the inverse metric of the mapping multiplied by the jacobian

     procedure :: init  !< Initialize the mapping class

     procedure :: init_from_file !< Initialize the mapping class with parameters read from nml-file

     procedure :: free !< Free mapping class

  end type sll_t_mapping_3d

contains


  function get_x1( self, xi ) result(x1)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64                              :: x1   !< physical coordinate
    !local variables
    sll_int32 :: box(3)
    sll_real64 :: xbox(3)

    if(self%flag)then
       x1 = self%x1_func( xi, self%params )
    else
       call convert_x_to_xbox( self, xi, xbox, box )
       x1 = sll_f_spline_pp_horner_3d(self%deg, self%x1_func_dofs_pp, xbox, box, self%n_cells)
    end if

  end function get_x1



  function get_x2( self, xi ) result(x2)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64                              :: x2   !< physical coordinate
    !local variables
    sll_int32 :: box(3)
    sll_real64 :: xbox(3)

    if(self%flag)then
       x2 = self%x2_func( xi, self%params )
    else
       call convert_x_to_xbox( self, xi, xbox, box )
       x2 = sll_f_spline_pp_horner_3d(self%deg, self%x2_func_dofs_pp, xbox, box, self%n_cells)
    end if

  end function get_x2



  function get_x3( self, xi ) result(x3)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64                              :: x3   !< physical coordinate
    !local variables
    sll_int32 :: box(3)
    sll_real64 :: xbox(3)

    if(self%flag)then
       x3 = self%x3_func( xi, self%params )
    else
       call convert_x_to_xbox( self, xi, xbox, box )
       x3 = sll_f_spline_pp_horner_3d(self%deg, self%x3_func_dofs_pp, xbox, box, self%n_cells)
    end if

  end function get_x3



  function get_x( self, xi ) result(x)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64                              :: x(3) !< physical coordinates
    !local variables
    sll_int32 :: box(3)
    sll_real64 :: xbox(3)

    if(self%flag)then
       x(1) = self%x1_func( xi, self%params )
       x(2) = self%x2_func( xi, self%params )
       x(3) = self%x3_func( xi, self%params )
    else
       call convert_x_to_xbox( self, xi, xbox, box )
       x(1) = sll_f_spline_pp_horner_3d(self%deg, self%x1_func_dofs_pp, xbox, box, self%n_cells)
       x(2) = sll_f_spline_pp_horner_3d(self%deg, self%x2_func_dofs_pp, xbox, box, self%n_cells)
       x(3) = sll_f_spline_pp_horner_3d(self%deg, self%x3_func_dofs_pp, xbox, box, self%n_cells)
    end if

  end function get_x

  function get_xi( self, x ) result(xi)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: x(3)  !< physical coordinates
    sll_real64                              :: xi(3) !< logical coordinates

    xi(1) = self%xi1_func( x, self%params )
    xi(2) = self%xi2_func( x, self%params )
    xi(3) = self%xi3_func( x, self%params )
  end function get_xi


  function jacobian( self, xi )result(x)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64                              :: x ! jacobian
    !local variables
    sll_real64 :: y(3,3)

    if(self%flag)then
       x = self%jacob( xi, self%params )
    else
       y = self%jacobian_matrix( xi )
       if( self%flag3d) then
          x =  y(1,1) * y(2,2) * y(3,3) +&
               y(1,2) * y(2,3) * y(3,1) +&
               y(1,3) * y(2,1) * y(3,2) -&
               y(1,3) * y(2,2) * y(3,1) -&
               y(1,2) * y(2,1) * y(3,3) -&
               y(1,1) * y(2,3) * y(3,2)
       else 
          x = y(3,3) * ( y(1,1) * y(2,2) - y(1,2) * y(2,1) ) 
       end if
    end if

  end function jacobian



  function jacobian_matrix( self, xi )result(y)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64                              :: y(3,3) !< jacobian matrix
    !local variables
    sll_int32 :: box(3)
    sll_real64 :: xbox(3)

    y = 0._f64

    if(self%flag)then    
       y(1,1) = self%j_matrix(1,1)%f( xi, self%params )
       y(2,2) = self%j_matrix(2,2)%f( xi, self%params )
       y(3,3) = self%j_matrix(3,3)%f( xi, self%params )
       y(1,2) = self%j_matrix(1,2)%f( xi, self%params )
       y(2,1) = self%j_matrix(2,1)%f( xi, self%params )
       if( self%flag3d ) then
          y(1,3) = self%j_matrix(1,3)%f( xi, self%params )
          y(2,3) = self%j_matrix(2,3)%f( xi, self%params )
          y(3,1) = self%j_matrix(3,1)%f( xi, self%params )
          y(3,2) = self%j_matrix(3,2)%f( xi, self%params )
       end if
    else
       call convert_x_to_xbox( self, xi, xbox, box )
       y(1,1) = sll_f_spline_pp_horner_3d_d1( self%deg, self%x1_func_dofs_pp, xbox, box, self%n_cells)
       y(2,2) = sll_f_spline_pp_horner_3d_d2( self%deg, self%x2_func_dofs_pp, xbox, box, self%n_cells)
       y(3,3) = sll_f_spline_pp_horner_3d_d3( self%deg, self%x3_func_dofs_pp, xbox, box, self%n_cells)
       y(1,2) = sll_f_spline_pp_horner_3d_d2( self%deg, self%x1_func_dofs_pp, xbox, box, self%n_cells)
       y(2,1) = sll_f_spline_pp_horner_3d_d1( self%deg, self%x2_func_dofs_pp, xbox, box, self%n_cells)
       if( self%flag3d ) then
          y(1,3) = sll_f_spline_pp_horner_3d_d3( self%deg, self%x1_func_dofs_pp, xbox, box, self%n_cells)
          y(2,3) = sll_f_spline_pp_horner_3d_d3( self%deg, self%x2_func_dofs_pp, xbox, box, self%n_cells) 
          y(3,1) = sll_f_spline_pp_horner_3d_d1( self%deg, self%x3_func_dofs_pp, xbox, box, self%n_cells)
          y(3,2) = sll_f_spline_pp_horner_3d_d2( self%deg, self%x3_func_dofs_pp, xbox, box, self%n_cells)
       end if
    end if

  end function jacobian_matrix



  function jacobian_matrix_transposed( self, xi )result(y)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64                              :: y(3,3) !< transposed jacobian matrix
    !local variables
    sll_real64 :: w(3,3)
    w = self%jacobian_matrix(xi)

    y(1,1) = w(1,1)
    y(1,2) = w(2,1)
    y(1,3) = w(3,1)
    y(2,1) = w(1,2)
    y(2,2) = w(2,2)
    y(2,3) = w(3,2)
    y(3,1) = w(1,3)
    y(3,2) = w(2,3)
    y(3,3) = w(3,3)

  end function jacobian_matrix_transposed

  function jacobian_matrix_inverse( self, xi )result(y)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64                              :: y(3,3) !< inverse jacobian matrix
    !local variables
    sll_real64 :: w(3,3)  !< jacobian matrix 

    w = self%jacobian_matrix(xi)

    if( self%flag3d ) then
       y(1,1) = w(2,2) * w(3,3) - w(2,3) * w(3,2)
       y(1,2) = w(1,3) * w(3,2) - w(1,2) * w(3,3)
       y(1,3) = w(1,2) * w(2,3) - w(1,3) * w(2,2)
       y(2,1) = w(2,3) * w(3,1) - w(2,1) * w(3,3)
       y(2,2) = w(1,1) * w(3,3) - w(1,3) * w(3,1)
       y(2,3) = w(1,3) * w(2,1) - w(1,1) * w(2,3)
       y(3,1) = w(2,1) * w(3,2) - w(2,2) * w(3,1)
       y(3,2) = w(1,2) * w(3,1) - w(1,1) * w(3,2)
       y(3,3) = w(1,1) * w(2,2) - w(1,2) * w(2,1)

       y = y/(w(1,1) * w(2,2) * w(3,3) +&
            w(1,2) * w(2,3) * w(3,1) +&
            w(1,3) * w(2,1) * w(3,2) -&
            w(1,3) * w(2,2) * w(3,1) -&
            w(1,2) * w(2,1) * w(3,3) -&
            w(1,1) * w(2,3) * w(3,2) )


    else 
       y = 0._f64
       y(1,1) = w(2,2) /( w(1,1) * w(2,2) - w(1,2) * w(2,1) )
       y(1,2) = - w(1,2) /( w(1,1) * w(2,2) - w(1,2) * w(2,1) ) 
       y(2,1) = - w(2,1) /( w(1,1) * w(2,2) - w(1,2) * w(2,1) ) 
       y(2,2) = w(1,1) /( w(1,1) * w(2,2) - w(1,2) * w(2,1) ) 
       y(3,3) = 1._f64/w(3,3)
    end if

  end function jacobian_matrix_inverse



  function jacobian_matrix_inverse_transposed( self, xi )result(y)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64                              :: y(3,3) !< jacobian matrix
    !local variable
    sll_real64                              :: w(3,3) !< transposed inverse jacobian matrix

    w=self%jacobian_matrix_inverse(xi)

    y(1,1) = w(1,1) 
    y(1,2) = w(2,1)
    y(1,3) = w(3,1)
    y(2,1) = w(1,2) 
    y(2,2) = w(2,2) 
    y(2,3) = w(3,2) 
    y(3,1) = w(1,3) 
    y(3,2) = w(2,3)
    y(3,3) = w(3,3) 

  end function jacobian_matrix_inverse_transposed



  function metric( self, xi )result(g)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64                              :: g(3,3) !< metric
    sll_real64                              :: y(3,3) !< jacobian matrix

    y=self%jacobian_matrix( xi )
    g(1,1) = y(1,1)**2+y(2,1)**2+y(3,1)**2
    g(1,2) = y(1,1)*y(1,2)+y(2,1)*y(2,2)+y(3,1)*y(3,2)
    g(1,3) = y(1,1)*y(1,3)+y(2,1)*y(2,3)+y(3,1)*y(3,3)
    g(2,1) = g(1,2)
    g(2,2) = y(1,2)**2+y(2,2)**2+y(3,2)**2
    g(2,3) = y(1,2)*y(1,3)+y(2,2)*y(2,3)+y(3,2)*y(3,3)
    g(3,1) = g(1,3)
    g(3,2) = g(2,3)
    g(3,3) = y(1,3)**2+y(2,3)**2+y(3,3)**2

  end function metric


  function metric_single( self, xi, component1, component2 )result(g)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_int32,               intent( in )   :: component1, component2 !< components of the wished entry
    sll_real64                              :: g !< single entry of the metric
    !local variable
    sll_real64                              :: w(3,3)
    w = self%jacobian_matrix(xi)

    g = w(1,component1)*w(1,component2)+&
         w(2,component1)*w(2,component2)+&
         w(3,component1)*w(3,component2)

  end function metric_single

  function metric_inverse( self, xi )result(g)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_real64               :: g(3,3) !< inverse metric
    sll_real64               :: y(3,3) !< inverse jacobian matrix

    y=self%jacobian_matrix_inverse_transposed(xi)
    g(1,1) = y(1,1)**2+y(2,1)**2+y(3,1)**2
    g(1,2) = y(1,1)*y(1,2)+y(2,1)*y(2,2)+y(3,1)*y(3,2)
    g(1,3) = y(1,1)*y(1,3)+y(2,1)*y(2,3)+y(3,1)*y(3,3)
    g(2,1) = g(1,2)
    g(2,2) = y(1,2)**2+y(2,2)**2+y(3,2)**2
    g(2,3) = y(1,2)*y(1,3)+y(2,2)*y(2,3)+y(3,2)*y(3,3)
    g(3,1) = g(1,3)
    g(3,2) = g(2,3)
    g(3,3) = y(1,3)**2+y(2,3)**2+y(3,3)**2

  end function metric_inverse

  function metric_inverse_single( self, xi, component1, component2 )result(g)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_int32,  intent( in )   :: component1, component2 !< components of the wished entry
    sll_real64                 :: g !< single entry of the inverse metric
    sll_real64                 :: y(3,3) !< inverse jacobian matrix 

    y = self%jacobian_matrix_inverse(xi)
    g = y(component1,1)*y(component2,1)+&
         y(component1,2)*y(component2,2)+&
         y(component1,3)*y(component2,3)

  end function metric_inverse_single

  function metric_single_jacobian( self, xi, component1, component2 )result(g)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_int32,  intent( in )   :: component1, component2 !< components of the wished entry
    sll_real64                 :: g !< single entry of the metric divided by the jacobian
    !local variables
    sll_real64                 :: y(3,3) !< jacobian matrix

    y = self%jacobian_matrix(xi)

    if( self%flag3d) then
       g =  (y(1,component1)*y(1,component2)+&
            y(2,component1)*y(2,component2)+&
            y(3,component1)*y(3,component2) )/  (y(1,1) * y(2,2) * y(3,3) +&
            y(1,2) * y(2,3) * y(3,1) +&
            y(1,3) * y(2,1) * y(3,2) -&
            y(1,3) * y(2,2) * y(3,1) -&
            y(1,2) * y(2,1) * y(3,3) -&
            y(1,1) * y(2,3) * y(3,2))
    else
       g =  (y(1,component1)*y(1,component2)+&
            y(2,component1)*y(2,component2)+&
            y(3,component1)*y(3,component2) )/ (y(3,3) * ( y(1,1) * y(2,2) - y(1,2) * y(2,1) ) )
    end if



  end function metric_single_jacobian

  function metric_inverse_single_jacobian( self, xi, component1, component2 )result(g)
    class(sll_t_mapping_3d), intent( inout )   :: self !< coordinate transformation
    sll_real64,              intent( in )   :: xi(3)  !< logical coordinates
    sll_int32,  intent( in )   :: component1, component2 !< components of the wished entry
    sll_real64                 :: g !< single entry of the inverse metric multiplied by the jacobian
    !local variables
    sll_real64                 :: w(3,3), y(3,3) !< jacobian matrix


    if( self%flag3d) then
       y = self%jacobian_matrix_inverse(xi)
       w(1,1) = self%jacobian(xi)
       g =  (y(component1,1)*y(component2,1)+&
            y(component1,2)*y(component2,2)+&
            y(component1,3)*y(component2,3)) *  w(1,1)
    else
       w = self%jacobian_matrix(xi)
       y = 0._f64
       y(1,1) = w(2,2) * w(3,3) 
       y(1,2) = - w(1,2) * w(3,3)
       y(2,1) = - w(2,1) * w(3,3)
       y(2,2) = w(1,1) * w(3,3) 
       y(3,3) = w(1,1) * w(2,2) - w(1,2) * w(2,1)

       g =  (y(component1,1)*y(component2,1)+&
            y(component1,2)*y(component2,2)+&
            y(component1,3)*y(component2,3))/ (w(3,3) * ( w(1,1) * w(2,2) - w(1,2) * w(2,1) ) )
    end if

  end function metric_inverse_single_jacobian


  subroutine init(  self, params, x1_func, x2_func, x3_func, jac11, jac12, jac13, jac21, jac22, jac23, jac31, jac32, jac33, jacob, xi1_func, xi2_func, xi3_func, flag2d, flag3d, n_cells, deg, Lx, volume  )
    class(sll_t_mapping_3d),  intent(   out ) :: self !< coordinate transformation
    sll_real64, dimension(:), intent( in    ) :: params !< transformation parameters
    procedure(sll_i_eval_function)            :: x1_func !< transformation function
    procedure(sll_i_eval_function)            :: x2_func !< transformation function
    procedure(sll_i_eval_function)            :: x3_func !< transformation function
    procedure(sll_i_eval_function), optional            :: jac11 !< entry of jacobian matrix
    procedure(sll_i_eval_function), optional            :: jac12 !< entry of jacobian matrix
    procedure(sll_i_eval_function), optional            :: jac13 !< entry of jacobian matrix
    procedure(sll_i_eval_function), optional            :: jac21 !< entry of jacobian matrix
    procedure(sll_i_eval_function), optional            :: jac22 !< entry of jacobian matrix
    procedure(sll_i_eval_function), optional            :: jac23 !< entry of jacobian matrix
    procedure(sll_i_eval_function), optional            :: jac31 !< entry of jacobian matrix
    procedure(sll_i_eval_function), optional            :: jac32 !< entry of jacobian matrix
    procedure(sll_i_eval_function), optional            :: jac33 !< entry of jacobian matrix
    procedure(sll_i_eval_function), optional            :: jacob !< determinant of jacobian matrix
    procedure(sll_i_eval_function), optional            :: xi1_func !< transformation function
    procedure(sll_i_eval_function), optional            :: xi2_func !< transformation function
    procedure(sll_i_eval_function), optional            :: xi3_func !< transformation function
    logical,                        optional            :: flag2d !< logical flag for mappings mixing the first two logical coordiantes
    logical,                        optional            :: flag3d  !< logical flag for mappings mixing all three logical coordinates
    sll_int32, optional :: n_cells(3) !< number of cells (and grid points)
    sll_int32, optional :: deg(3)     !< spline deg 
    sll_real64, optional :: Lx(3)     !< length of the physical domain  
    sll_real64, optional :: volume     !< volume of the physical domain  
    !local variables
    sll_int32  :: ierr
    sll_real64, allocatable :: x1_func_dofs(:)
    sll_real64, allocatable :: x2_func_dofs(:)
    sll_real64, allocatable :: x3_func_dofs(:)
    sll_real64, allocatable :: conp(:,:)
    type(sll_t_matrix_csr) :: matrix(3)
    type(sll_t_linear_solver_mgmres) :: matrix_solver(3)
    type(sll_t_linear_solver_kron) :: kron_solver
    sll_real64, allocatable  :: rhs(:,:)
    sll_real64, allocatable :: tk(:,:), xk(:,:)
    sll_int32 :: ntotal0, n_dofs0(3), ntotal1, n_dofs1(3)
    sll_int32 :: i, j, k
    sll_real64, allocatable :: val(:)
    sll_real64 :: tau
    sll_int32 :: deg1(3)
    sll_real64, allocatable :: scratch0(:), scratch1(:), work0(:), work1(:)

    SLL_ALLOCATE(self%j_matrix(1:3,1:3), ierr)
    SLL_ALLOCATE(self%params(1:size(params)),ierr)
    self%params = params
    self%x1_func => x1_func
    self%x2_func => x2_func
    self%x3_func => x3_func

    if( present(xi1_func) .and. present(xi2_func) .and. present(xi3_func) ) then
       self%xi1_func => xi1_func
       self%xi2_func => xi2_func
       self%xi3_func => xi3_func
       !self%cylinder_params = self%params(1:3)
       self%inverse = .true.
    end if
    if( present(jac11) .and. present(jacob)) then
       self%j_matrix(1,1)%f => jac11
       self%j_matrix(1,2)%f => jac12
       self%j_matrix(1,3)%f => jac13
       self%j_matrix(2,1)%f => jac21
       self%j_matrix(2,2)%f => jac22
       self%j_matrix(2,3)%f => jac23
       self%j_matrix(3,1)%f => jac31
       self%j_matrix(3,2)%f => jac32
       self%j_matrix(3,3)%f => jac33

       self%flag= .true.
       self%jacob => jacob
    end if
    if(params(1) == 0._f64) then
       self%singular = .true.
       print*, 'Error: singular mapping not yet implemented'
    end if


    if( present(flag2d) )then
       self%flag2d = flag2d
    end if

    if( present(flag3d) )then
       self%flag2d = flag2d
       self%flag3d = flag3d
    end if

    if( present(Lx) )then
       self%Lx = Lx
    end if

    if( present(volume) )then
       self%volume = volume
    end if

    if( present(n_cells) .and. present(deg)  ) then
       self%n_cells = n_cells
       self%deg = deg
       deg1 = deg-1
       self%flag= .false.

       n_dofs0 = self%n_cells
       n_dofs0(1) = self%n_cells(1)+self%deg(1)
       n_dofs0(3) = self%n_cells(3)+self%deg(3)
       n_dofs1 = n_dofs0
       n_dofs1(1) = self%n_cells(1)+deg1(1)
       ntotal0 = product(n_dofs0)
       ntotal1 = product(n_dofs1)

       call sll_s_spline_pp_init_3d(self%spline_3d, self%deg, self%n_cells, [1,0,1])
       allocate( x1_func_dofs(1:ntotal0) )
       allocate( x2_func_dofs(1:ntotal0) )
       allocate( x3_func_dofs(1:ntotal0) )
       allocate( self%x1_func_dofs_pp(1:product(self%deg+1), product(self%n_cells)) )
       allocate( self%x2_func_dofs_pp(1:product(self%deg+1), product(self%n_cells)) )
       allocate( self%x3_func_dofs_pp(1:product(self%deg+1), product(self%n_cells)) )
       x1_func_dofs = 0._f64
       x2_func_dofs = 0._f64
       x3_func_dofs = 0._f64
       self%x1_func_dofs_pp = 0._f64
       self%x2_func_dofs_pp = 0._f64
       self%x3_func_dofs_pp = 0._f64

       allocate( tk(maxval(self%n_cells)+2*maxval(self%deg)+1, 3) )
       allocate( xk(maxval(n_dofs0), 3) )
       allocate( rhs(1:ntotal0, 3) )
       tk = 0._f64
       xk = 0._f64
       rhs = 0._f64


       ! knot sequence
       !open
       tk(1:self%deg(1),1) = 0._f64
       do i = 1, self%n_cells(1)+1
          tk(self%deg(1)+i,1) = real(i-1,f64)/real(self%n_cells(1),f64)
       end do
       tk(self%n_cells(1)+self%deg(1)+2:self%n_cells(1)+2*self%deg(1)+1,1) = 1._f64

       tk(1:self%deg(3),3) = 0._f64
       do i = 1, self%n_cells(3)+1
          tk(self%deg(3)+i,3) = real(i-1,f64)/real(self%n_cells(3),f64)
       end do
       tk(self%n_cells(3)+self%deg(3)+2:self%n_cells(3)+2*self%deg(3)+1,3) = 1._f64

       !periodic
       do j = 2, 2
          do i = 1, self%n_cells(j)
             tk(self%deg(j)+i,j) = real(i-1,f64)/real(self%n_cells(j),f64)
          end do
          do i = 1, self%deg(j)
             tk(i,j) = tk(self%n_cells(j)+i,j)-1._f64  
             tk(self%n_cells(j)+self%deg(j)+i,j) = tk(self%deg(j)+i,j) + 1._f64
          end do
       end do


       !greville points
       do i = 1, n_dofs0(1)
          xk(i,1) = sum(tk(1+i:i+self%deg(1),1))/real(self%deg(1),f64)
       end do

       do i = 1, n_dofs0(3)
          xk(i,3) = sum(tk(1+i:i+self%deg(3),3))/real(self%deg(3),f64)
       end do

       !periodic
       do j = 2, 2
          if( modulo(real(self%deg(j),f64),2._f64) == 0._f64) then
             xk(1:n_dofs0(j),j) = 0.5_f64*(tk(self%deg(j)+1:self%deg(j)+self%n_cells(j),j)+tk(self%deg(j)+2:self%deg(j)+1+self%n_cells(j),j))
          else
             xk(1:n_dofs0(j),j) = tk(self%deg(j)+1:self%deg(j)+self%n_cells(j),j)
          end if
       end do

       ! interpolation matrix
       call calculate_interpolation_matrix_1d( self%n_cells(1), self%deg(1), xk(1:n_dofs0(1),1), self%spline_3d%spline1, matrix(1) )
       call calculate_interpolation_matrix_1d_periodic( self%n_cells(2), self%deg(2), self%spline_3d%spline2, matrix(2) )
       call calculate_interpolation_matrix_1d( self%n_cells(3), self%deg(3), xk(1:n_dofs0(3),3), self%spline_3d%spline3, matrix(3) )

       do j = 1, 3
          call matrix_solver(j)%create( matrix(j) )
          matrix_solver(j)%atol = 1d-13
       end do
       !matrix_solver(2)%verbose = .true.

       call kron_solver%create( linear_solver_a=matrix_solver(1), &
            linear_solver_b=matrix_solver(2), &
            linear_solver_c=matrix_solver(3) )


       !rhs evaluation from transformation function
       do  k = 1, n_dofs0(3)
          do j  = 1, n_dofs0(2)
             do i  = 1, n_dofs0(1)
                rhs( i + (j-1)*n_dofs0(1) + (k-1)*n_dofs0(1)*n_dofs0(2), 1 ) = x1_func( [xk(i,1), xk(j,2), xk(k,3)], params )
                rhs( i + (j-1)*n_dofs0(1) + (k-1)*n_dofs0(1)*n_dofs0(2), 2 ) = x2_func( [xk(i,1), xk(j,2), xk(k,3)], params )
                rhs( i + (j-1)*n_dofs0(1) + (k-1)*n_dofs0(1)*n_dofs0(2), 3 ) = x3_func( [xk(i,1), xk(j,2), xk(k,3)], params )
             end do
          end do
       end do

       call kron_solver%solve( rhs(:,1), x1_func_dofs)
       call kron_solver%solve( rhs(:,2), x2_func_dofs)
       call kron_solver%solve( rhs(:,3), x3_func_dofs)


       call sll_s_spline_pp_b_to_pp_3d_clamped_2full( self%spline_3d, self%n_cells, x1_func_dofs, self%x1_func_dofs_pp)
       call sll_s_spline_pp_b_to_pp_3d_clamped_2full( self%spline_3d, self%n_cells, x2_func_dofs, self%x2_func_dofs_pp)
       call sll_s_spline_pp_b_to_pp_3d_clamped_2full( self%spline_3d, self%n_cells, x3_func_dofs, self%x3_func_dofs_pp)

    end if

  end subroutine init

  subroutine init_from_file(  self, filename   )
    class(sll_t_mapping_3d),  intent(   out ) :: self !< coordinate transformation
    character(len=*),         intent( in    ) :: filename
    !local variables
    sll_int32  :: input_file
    sll_int32  :: io_stat
    sll_int32  :: nparams
    sll_real64 :: params(12)
    character(len=256) :: mapping_case
    sll_int32 :: n_cells(3) = 3
    sll_int32 :: deg(3) = 1
    sll_real64 :: vol = 1._f64
    sll_real64 :: Lx(3) = 1._f64
    namelist /trafo/ nparams, params, mapping_case, n_cells, deg 

    ! Read parameters from file
    open(newunit = input_file, file=trim(filename), IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_mapping() failed to open file ', filename
       STOP
    end if

    read(input_file, trafo)
    close(input_file)

    !> mapping cases defined in sll_m_3d_coordinate_transformations
    select case(trim(mapping_case))
    case( "cylindrical_sqrt" )
       Lx(1:2) = 2._f64*(params(2)-params(1))
       Lx(3) = params(3)
       vol = sll_p_pi * (params(2)-params(1))**2 * params(3)
       call self%init(params, sll_f_cylindrical_sqrt_x1, sll_f_cylindrical_sqrt_x2, sll_f_cylindrical_sqrt_x3, &
            sll_f_cylindrical_sqrt_jac11, sll_f_cylindrical_sqrt_jac12, sll_f_cylindrical_sqrt_jac13,&
            sll_f_cylindrical_sqrt_jac21, sll_f_cylindrical_sqrt_jac22, sll_f_cylindrical_sqrt_jac23,&
            sll_f_cylindrical_sqrt_jac31, sll_f_cylindrical_sqrt_jac32, sll_f_cylindrical_sqrt_jac33,&
            sll_f_cylindrical_sqrt_jacobian, Lx=Lx, volume=vol)
    case( "cylindrical_sqrt_discrete" )
       Lx(1:2) = 2._f64*(params(2)-params(1))
       Lx(3) = params(3)
       vol = sll_p_pi * (params(2)-params(1))**2 * params(3)
       call self%init(params, sll_f_cylindrical_sqrt_x1, sll_f_cylindrical_sqrt_x2, sll_f_cylindrical_sqrt_x3, &
            n_cells=n_cells, deg=deg, Lx=Lx, volume=vol ) 
    case( "cylindrical_sqrt_inverse" )
       Lx(1:2) = 2._f64*(params(2)-params(1))
       Lx(3) = params(3)
       vol = sll_p_pi * (params(2)-params(1))**2 * params(3)
       call self%init(params, sll_f_cylindrical_sqrt_x1, sll_f_cylindrical_sqrt_x2, sll_f_cylindrical_sqrt_x3, &
            sll_f_cylindrical_sqrt_jac11, sll_f_cylindrical_sqrt_jac12, sll_f_cylindrical_sqrt_jac13,&
            sll_f_cylindrical_sqrt_jac21, sll_f_cylindrical_sqrt_jac22, sll_f_cylindrical_sqrt_jac23,&
            sll_f_cylindrical_sqrt_jac31, sll_f_cylindrical_sqrt_jac32, sll_f_cylindrical_sqrt_jac33,&
            sll_f_cylindrical_sqrt_jacobian, &
            sll_f_cylindrical_sqrt_xi1, sll_f_cylindrical_sqrt_xi2, sll_f_cylindrical_sqrt_xi3, &
            Lx=Lx, volume=vol)
    case( "cylindrical_sqrt_discrete_inverse" )
       Lx(1:2) = 2._f64*(params(2)-params(1))
       Lx(3) = params(3)
       vol = sll_p_pi * (params(2)-params(1))**2 * params(3)
       call self%init(params, sll_f_cylindrical_sqrt_x1, sll_f_cylindrical_sqrt_x2, sll_f_cylindrical_sqrt_x3, &
            xi1_func=sll_f_cylindrical_sqrt_xi1, xi2_func=sll_f_cylindrical_sqrt_xi2, xi3_func=sll_f_cylindrical_sqrt_xi3, &
            n_cells=n_cells, deg=deg, Lx=Lx, volume=vol )
    case( "cylindrical" )
       Lx(1:2) = 2._f64*(params(2)-params(1))
       Lx(3) = params(3)
       vol = sll_p_pi * (params(2)-params(1))**2 * params(3)
       call self%init(params, sll_f_cylindrical_x1, sll_f_cylindrical_x2, sll_f_cylindrical_x3, &
            sll_f_cylindrical_jac11, sll_f_cylindrical_jac12, sll_f_cylindrical_jac13,&
            sll_f_cylindrical_jac21, sll_f_cylindrical_jac22, sll_f_cylindrical_jac23,&
            sll_f_cylindrical_jac31, sll_f_cylindrical_jac32, sll_f_cylindrical_jac33,&
            sll_f_cylindrical_jacobian, Lx=Lx, volume=vol)

    case( "cylindrical_discrete" )
       Lx(1:2) = 2._f64*(params(2)-params(1))
       Lx(3) = params(3)
       vol = sll_p_pi * (params(2)-params(1))**2 * params(3)
       call self%init(params, sll_f_cylindrical_x1, sll_f_cylindrical_x2, sll_f_cylindrical_x3, &
            flag2d = .true., n_cells=n_cells, deg=deg, Lx=Lx, volume=vol )
    case( "cylindrical_inverse" )
       Lx(1:2) = 2._f64*(params(2)-params(1))
       Lx(3) = params(3)
       vol = sll_p_pi * (params(2)-params(1))**2 * params(3)
       call self%init(params, sll_f_cylindrical_x1, sll_f_cylindrical_x2, sll_f_cylindrical_x3, &
            sll_f_cylindrical_jac11, sll_f_cylindrical_jac12, sll_f_cylindrical_jac13,&
            sll_f_cylindrical_jac21, sll_f_cylindrical_jac22, sll_f_cylindrical_jac23,&
            sll_f_cylindrical_jac31, sll_f_cylindrical_jac32, sll_f_cylindrical_jac33,&
            sll_f_cylindrical_jacobian, &
            sll_f_cylindrical_xi1, sll_f_cylindrical_xi2, sll_f_cylindrical_xi3, &
            Lx=Lx, volume=vol)
    case( "cylindrical_discrete_inverse" )
       Lx(1:2) = 2._f64*(params(2)-params(1))
       Lx(3) = params(3)
       vol = sll_p_pi * (params(2)-params(1))**2 * params(3)
       call self%init(params, sll_f_cylindrical_x1, sll_f_cylindrical_x2, sll_f_cylindrical_x3, &
            xi1_func=sll_f_cylindrical_xi1, xi2_func=sll_f_cylindrical_xi2, xi3_func=sll_f_cylindrical_xi3, &
            n_cells=n_cells, deg=deg, Lx=Lx, volume=vol )
    case( "spherical" )
       vol = 4._f64/3._f64 * sll_p_pi * (params(2)-params(1))**3
       call self%init(params, sll_f_spherical_x1, sll_f_spherical_x2, sll_f_spherical_x3, &
            sll_f_spherical_jac11, sll_f_spherical_jac12, sll_f_spherical_jac13,&
            sll_f_spherical_jac21, sll_f_spherical_jac22, sll_f_spherical_jac23,&
            sll_f_spherical_jac31, sll_f_spherical_jac32, sll_f_spherical_jac33,&
            sll_f_spherical_jacobian, flag3d = .true., Lx=Lx, volume=vol)
    case( "elliptical" )
       Lx(1) = 3.1_f64*params(2)
       Lx(1) = 2.4_f64*params(2)
       Lx(3) = params(3)
       vol = sll_p_pi * 7.44_f64* params(2)**2 * params(3)
       call self%init(params, sll_f_elliptical_x1, sll_f_elliptical_x2, sll_f_elliptical_x3, &
            sll_f_elliptical_jac11, sll_f_elliptical_jac12, sll_f_elliptical_jac13,&
            sll_f_elliptical_jac21, sll_f_elliptical_jac22, sll_f_elliptical_jac23,&
            sll_f_elliptical_jac31, sll_f_elliptical_jac32, sll_f_elliptical_jac33,&
            sll_f_elliptical_jacobian, Lx=Lx, volume=vol)
    case( "parallelogram" )
       Lx = params(1:3)
       vol = product(params(1:3))
       call self%init(params, sll_f_parallelogram_x1, sll_f_parallelogram_x2, sll_f_parallelogram_x3, &
            sll_f_parallelogram_jac11, sll_f_parallelogram_jac12, sll_f_parallelogram_jac13,&
            sll_f_parallelogram_jac21, sll_f_parallelogram_jac22, sll_f_parallelogram_jac23,&
            sll_f_parallelogram_jac31, sll_f_parallelogram_jac32, sll_f_parallelogram_jac33,&
            sll_f_parallelogram_jacobian, flag2d = .true., Lx=Lx, volume=vol)
    case( "coltest" )
       Lx = params(1:3)
       vol = product(params(1:3))
       call self%init(params, sll_f_coltest_x1, sll_f_coltest_x2, sll_f_coltest_x3, &
            sll_f_coltest_jac11, sll_f_coltest_jac12, sll_f_coltest_jac13,&
            sll_f_coltest_jac21, sll_f_coltest_jac22, sll_f_coltest_jac23,&
            sll_f_coltest_jac31, sll_f_coltest_jac32, sll_f_coltest_jac33,&
            sll_f_coltest_jacobian, flag2d = .true., Lx=Lx, volume=vol)
    case( "colbound" )
       Lx = params(1:3)
       vol = product(params(1:3))
       call self%init(params, sll_f_colbound_x1, sll_f_colbound_x2, sll_f_colbound_x3, &
            sll_f_colbound_jac11, sll_f_colbound_jac12, sll_f_colbound_jac13,&
            sll_f_colbound_jac21, sll_f_colbound_jac22, sll_f_colbound_jac23,&
            sll_f_colbound_jac31, sll_f_colbound_jac32, sll_f_colbound_jac33,&
            sll_f_colbound_jacobian, flag2d = .true., Lx=Lx, volume=vol)
    case( "colella" )
       Lx = params(1:3)
       vol = product(params(1:3))
       call self%init(params, sll_f_colella_x1, sll_f_colella_x2, sll_f_colella_x3, &
            sll_f_colella_jac11, sll_f_colella_jac12, sll_f_colella_jac13,&
            sll_f_colella_jac21, sll_f_colella_jac22, sll_f_colella_jac23,&
            sll_f_colella_jac31, sll_f_colella_jac32, sll_f_colella_jac33,&
            sll_f_colella_jacobian, flag2d = .true., Lx=Lx, volume=vol)
    case( "colella_discrete" )
       Lx = params(1:3)
       vol = product(params(1:3))
       call self%init(params, sll_f_colella_x1, sll_f_colella_x2, sll_f_colella_x3, &
            flag2d = .true., n_cells=n_cells, deg=deg, Lx=Lx, volume=vol )
    case( "orthogonal" )
       Lx = params(1:3)
       vol = product(params(1:3))
       call self%init(params, sll_f_orthogonal_x1, sll_f_orthogonal_x2, sll_f_orthogonal_x3, &
            sll_f_orthogonal_jac11, sll_f_orthogonal_jac12, sll_f_orthogonal_jac13,&
            sll_f_orthogonal_jac21, sll_f_orthogonal_jac22, sll_f_orthogonal_jac23,&
            sll_f_orthogonal_jac31, sll_f_orthogonal_jac32, sll_f_orthogonal_jac33,&
            sll_f_orthogonal_jacobian, Lx=Lx, volume=vol )
    case( "orthogonal_discrete" )
       vol = product(params(1:3))
       call self%init(params, sll_f_orthogonal_x1, sll_f_orthogonal_x2, sll_f_orthogonal_x3, &
            n_cells=n_cells, deg=deg, Lx=Lx, volume=vol)
    case( "polynomial" )
       Lx = params(1:3)
       vol = product(params(1:3))
       call self%init(params, sll_f_polynomial_x1, sll_f_polynomial_x2, sll_f_polynomial_x3, &
            sll_f_polynomial_jac11, sll_f_polynomial_jac12, sll_f_polynomial_jac13,&
            sll_f_polynomial_jac21, sll_f_polynomial_jac22, sll_f_polynomial_jac23,&
            sll_f_polynomial_jac31, sll_f_polynomial_jac32, sll_f_polynomial_jac33,&
            sll_f_polynomial_jacobian, Lx=Lx, volume=vol)
    case( "polynomial_discrete" )
       Lx = params(1:3)
       vol = product(params(1:3))
       call self%init(params, sll_f_polynomial_x1, sll_f_polynomial_x2, sll_f_polynomial_x3, &
            n_cells=n_cells, deg=deg, Lx=Lx, volume=vol)
    case( "scaling")
       Lx = params(1:3)
       vol = product(params(1:3))
       call self%init(params, sll_f_scaling_x1, sll_f_scaling_x2, sll_f_scaling_x3, &
            sll_f_scaling_jac11, sll_f_scaling_jac12, sll_f_scaling_jac13,&
            sll_f_scaling_jac21, sll_f_scaling_jac22, sll_f_scaling_jac23,&
            sll_f_scaling_jac31, sll_f_scaling_jac32, sll_f_scaling_jac33,&
            sll_f_scaling_jacobian, Lx=Lx, volume=vol)
    case( "scaling_discrete")
       Lx = params(1:3)
       vol = product(params(1:3))
       call self%init(params, sll_f_scaling_x1, sll_f_scaling_x2, sll_f_scaling_x3, &
            n_cells=n_cells, deg=deg, Lx=Lx, volume=vol)
    case default
       print*, 'error init_mapping: mapping name not implemented'
    end select

  end subroutine init_from_file

  subroutine free (self )
    class(sll_t_mapping_3d), intent(inout) :: self


    self%x1_func => null()
    self%x2_func => null()
    self%x3_func => null()
    DEALLOCATE(self%j_matrix)
    DEALLOCATE(self%params)

  end subroutine free

  subroutine calculate_interpolation_matrix_1d( n_cells, deg, xk, spline, matrix )
    sll_int32,  intent( in ) :: n_cells !< number of cells (and grid points)
    sll_int32,  intent( in ) :: deg     !< spline deg
    sll_real64, intent( in ) :: xk(:)
    type( sll_t_spline_pp_1d), intent( in ) :: spline
    type(sll_t_matrix_csr), intent( out ) :: matrix
    !local variables
    sll_int32 :: row, column, ind, i
    sll_int32 :: box, n_nnz, n_dofs
    sll_real64 :: xi

    n_dofs = n_cells + deg
    n_nnz = n_dofs*(deg+1)

    call matrix%create( n_rows=n_dofs, n_cols=n_dofs, n_nnz=n_nnz )
    matrix%arr_ia(1)=1
    do row = 2, n_dofs+1
       matrix%arr_ia(row) = matrix%arr_ia(row-1) + deg+1
    end do

    ind = 1
    do row = 1, floor(0.5_f64*real(deg,f64))
       do column = 1, deg+1
          matrix%arr_ja(ind) = column
          ind = ind+1
       end do
    end do
    do row = 1, n_cells
       do column = row, row+deg
          matrix%arr_ja(ind) = column
          ind = ind+1
       end do
    end do
    do row = 1, ceiling(0.5_f64*real(deg,f64))
       do column = n_cells, n_dofs
          matrix%arr_ja(ind) = column
          ind = ind+1
       end do
    end do
    matrix%arr_a = 0.0_f64
    SLL_ASSERT( ind == matrix%n_nnz+1 )
    ind = 1

    do row = 1, n_dofs
       xi = xk(row)*real(n_cells,f64)
       box = floor( xi ) + 1
       xi = xi - real(box-1, f64)
       if( box == n_cells + 1 ) then
          xi = 1._f64
          box = n_cells 
       end if

       if(box <= deg-1)then
          do i = 1, deg+1
             matrix%arr_a(1,ind) = sll_f_spline_pp_horner_1d( deg, spline%poly_coeffs_boundary_left(:,:,box), xi, i)
             ind = ind+1
          end do
       else if(box >= n_cells-deg+2)then
          do i = 1, deg+1
             matrix%arr_a(1,ind) = sll_f_spline_pp_horner_1d( deg, spline%poly_coeffs_boundary_right(:,:,box-1-n_cells+deg), xi, i)
             ind = ind+1
          end do
       else
          do i = 1, deg+1
             matrix%arr_a(1,ind) = sll_f_spline_pp_horner_1d( deg, spline%poly_coeffs, xi, i)
             ind = ind+1
          end do
       end if
    end do
    SLL_ASSERT( ind == matrix%n_nnz+1 )

  end subroutine calculate_interpolation_matrix_1d

  subroutine calculate_interpolation_matrix_1d_periodic( n_cells, deg, spline, matrix )
    sll_int32,  intent( in ) :: n_cells !< number of cells (and grid points)
    sll_int32,  intent( in ) :: deg     !< spline deg
    type( sll_t_spline_pp_1d), intent( in ) :: spline
    type(sll_t_matrix_csr), intent( out ) :: matrix
    !local variables
    sll_int32 :: row, column, ind, i
    sll_int32 :: n_nnz
    sll_real64 :: xi, val(1:deg+1)

    n_nnz = n_cells*(deg+1)

    call matrix%create( n_rows=n_cells, n_cols=n_cells, n_nnz=n_nnz )
    matrix%arr_ia(1)=1
    do row = 2, n_cells+1
       matrix%arr_ia(row) = matrix%arr_ia(row-1) + deg+1
    end do

    ind = 1
    do row = 1, deg
       do column = 1, row
          matrix%arr_ja(ind) = column
          ind = ind+1
       end do
       do column = n_cells+row-deg, n_cells
          matrix%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    do row = deg+1, n_cells
       do column = row-deg, row
          matrix%arr_ja(ind) = column
          ind = ind+1
       end do
    end do



    matrix%arr_a = 0.0_f64
    SLL_ASSERT( ind == matrix%n_nnz+1 )

    if( modulo(real(deg,f64),2._f64) == 0._f64)then
       xi = 0.5_f64
    else
       xi = 0._f64
    end if

    do i = 1, deg+1
       val(i) = sll_f_spline_pp_horner_1d( deg, spline%poly_coeffs, xi, i)
    end do

    ind = 1
    do row = 1, deg
       do i = deg+2-row, deg+1
          matrix%arr_a(1,ind) = val(i)
          ind = ind+1
       end do
       do i = 1, deg+1-row
          matrix%arr_a(1,ind) = val(i) 
          ind = ind+1
       end do
    end do

    do row = deg+1, n_cells
       matrix%arr_a(1,ind:ind+deg) = val
       ind = ind+deg+1
    end do

    SLL_ASSERT( ind == matrix%n_nnz+1 )

  end subroutine calculate_interpolation_matrix_1d_periodic

  subroutine convert_x_to_xbox( self, position, xi, box )
    class(sll_t_mapping_3d), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(3) !< Position of the particle
    sll_real64,                               intent( out )    :: xi(3) !< Position of the particle
    sll_int32,                                intent( out )    :: box(3) !< Position of the particle

    xi = position  * real(self%n_cells,f64)
    box = floor( xi ) + 1
    xi = xi - real(box-1, f64)

    if( box(1) == self%n_cells(1) + 1 ) then
       if( xi(1) == 0._f64)then 
          xi(1) = 1._f64
          box(1) = self%n_cells(1)
       else
          print*, 'box, x, xbox', box, position(1), xi(1)
          SLL_ERROR('convert_x_to_xbox', 'box1 value to high' ) 
       end if
    else if( box(1) == 0 ) then
       print*, 'box, x, xbox', box, position(1), xi(1)
       SLL_ERROR('convert_x_to_xbox', 'box1 value to low' ) 
    end if

    if( box(2) == self%n_cells(2) + 1 ) then
       if( xi(2) == 0._f64)then 
          xi(2) = 1._f64
          box(2) = self%n_cells(2)
       else
          print*, 'box, x, xbox', box, position(2), xi(2)
          SLL_ERROR('convert_x_to_xbox', 'box2 value to high' ) 
       end if
    else if( box(2) == 0 ) then
       print*, 'box, x, xbox', box, position(2), xi(2)
       SLL_ERROR('convert_x_to_xbox', 'box2 value to low' ) 
    end if


    if( box(3) == self%n_cells(3) + 1 ) then
       if( xi(3) == 0._f64)then 
          xi(3) = 1._f64
          box(3) = self%n_cells(3)
       else
          print*, 'box, x, xbox', box, position(3), xi(3)
          SLL_ERROR('convert_x_to_xbox', 'box3 value to high' ) 
       end if
    else if( box(3) == 0 ) then
       print*, 'box, x, xbox', box, position(3), xi(3)
       SLL_ERROR('convert_x_to_xbox', 'box3 value to low' ) 
    end if

  end subroutine convert_x_to_xbox


end module sll_m_mapping_3d
