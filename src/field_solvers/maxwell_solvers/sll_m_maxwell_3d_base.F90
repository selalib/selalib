!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations in 3D
!> @details
!> Contains the abstract class to create a Maxwell solver in 3D.

module sll_m_maxwell_3d_base
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_utilities, only: &
       sll_s_int2string

  use sll_m_profile_functions, only: &
       sll_t_profile_functions

  implicit none

  public :: &
       sll_i_function_3d_real64, &
       sll_c_maxwell_3d_base
  !    sll_s_plot_two_fields_3d

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, abstract :: sll_c_maxwell_3d_base
     
     sll_real64 :: delta_x(3)       !< cell size
     sll_real64 :: volume           !< product(delta_x)
     sll_int32  :: n_cells(3)       !< number of cells (and gridpoints)
     sll_int32  :: n_dofs(3)        !< number of Degrees of Freedom 
     sll_int32  :: n_total          !< total number of gridcells
     sll_int32  :: n_total0         !< total number of Dofs for 0form
     sll_int32  :: n_total1         !< total number of Dofs for 1form
     sll_int32  :: s_deg_0(3)       !< spline degree 0-forms
     sll_int32  :: s_deg_1(3)       !< spline degree 1-forms
     sll_real64 :: Lx(3)            !< length of Periodic domain
     sll_real64 :: mass_solver_tolerance = 1d-12 !< tolerance for the mass solver
     sll_real64 :: poisson_solver_tolerance = 1d-12 !< tolerance for the Poisson solver
     sll_real64 :: solver_tolerance = 1d-12 !< tolerance for the field solver
     type(sll_t_profile_functions) :: profile !< temperature and density profiles

   contains

     procedure(compute_field1_from_field2), deferred :: &
          compute_E_from_B !< Solve E and B part of Ampere's law with B constant in time
     procedure(compute_field1_from_field2), deferred :: &
          compute_B_from_E !< Solve Faraday equation with E constant in time
     procedure(compute_field_from_field), deferred :: &
          compute_E_from_rho !< Solve E from rho using Poisson
     procedure(compute_field_from_field), deferred :: &
          compute_rho_from_E !< Compute rho from E by Gauss law
     procedure(compute_E_from_j_3d), deferred :: &
          compute_E_from_j !< Solve E from time integrated current (second part of Ampere's law)
     procedure(compute_phi_e_from_field), deferred :: &
          compute_phi_from_rho !< Compute phi from rho 
     procedure(compute_phi_e_from_field), deferred :: &
          compute_phi_from_j !< Compute phi from j 
     procedure(update_dofs_function), deferred :: &
          compute_rhs_from_function !< Compute the right-hand-side for a given function f. For Galerkin it is the inner product with the basis functions. For Collocation it is simply a function evaluation at the grid points.
     procedure(update_dofs_function), deferred :: &
          L2projection !< L2 projection
     procedure(norm_squared), deferred :: &
          L2norm_squared !< Square of the L2norm
     procedure(inner_product), deferred :: &
          inner_product
     procedure(empty), deferred :: &
          free !< destructor

     procedure(compute_field_from_field), deferred :: &
          multiply_c
     procedure(compute_field_from_field), deferred :: &
          multiply_ct
     procedure(compute_field_from_field), deferred :: &
          multiply_g
     procedure(compute_field_from_field), deferred :: &
         multiply_gt
     procedure(multiply_mass), deferred :: &
          multiply_mass
     procedure(multiply_mass_inverse), deferred :: &
          multiply_mass_inverse
     procedure :: compute_curl_part

  end type sll_c_maxwell_3d_base



  !---------------------------------------------------------------------------!
  abstract interface 
     subroutine compute_field1_from_field2( self, delta_t, field_in, field_out )
       use sll_m_working_precision
       import sll_c_maxwell_3d_base     
       class(sll_c_maxwell_3d_base) :: self
       sll_real64, intent( in    )  :: delta_t
       sll_real64, intent( in    )  :: field_in(:)
       sll_real64, intent( inout )  :: field_out(:)
       
     end subroutine compute_field1_from_field2
  end interface



  !---------------------------------------------------------------------------!
  abstract interface 
     subroutine compute_field_from_field( self, field_in, field_out )
       use sll_m_working_precision
       import sll_c_maxwell_3d_base     
       class(sll_c_maxwell_3d_base) :: self
       sll_real64, intent( in    )  :: field_in(:)
       sll_real64, intent(   out )  :: field_out(:)
       
     end subroutine compute_field_from_field
  end interface

  !---------------------------------------------------------------------------!
  abstract interface 
     subroutine compute_phi_e_from_field( self, field_in, field_out, efield_dofs )
       use sll_m_working_precision
       import sll_c_maxwell_3d_base     
       class(sll_c_maxwell_3d_base) :: self
       sll_real64, intent( in    )  :: field_in(:)
       sll_real64, intent( inout )  :: field_out(:)
       sll_real64, intent(   out )  :: efield_dofs(:)
       
     end subroutine compute_phi_e_from_field
  end interface

  !---------------------------------------------------------------------------!  
  abstract interface
     subroutine compute_E_from_j_3d( self, current, E, component )
       use sll_m_working_precision
       import sll_c_maxwell_3d_base
       class(sll_c_maxwell_3d_base)          :: self
       sll_real64, intent( in    )           :: current(:)
       sll_real64, intent( inout )           :: E(:)
       sll_int32,  intent( in    ), optional :: component
       
     end subroutine compute_E_from_j_3d
  end interface



  !---------------------------------------------------------------------------!
  abstract interface
     !> 3d real function
     function sll_i_function_3d_real64(x)
       use sll_m_working_precision ! can't pass a header file because the
       ! preprocessor prevents double inclusion.
       ! It is very rare.
       sll_real64             :: sll_i_function_3d_real64
       sll_real64, intent(in) :: x(3)
     end function sll_i_function_3d_real64
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine update_dofs_function( self, form, component, coefs_dofs, func1, func2, func3 )
       use sll_m_working_precision
       import sll_c_maxwell_3d_base
       import sll_i_function_3d_real64
       class( sll_c_maxwell_3d_base)                  :: self          !< Maxwell solver object.
       sll_int32,  intent( in    )                    :: form          !< Specify if the function is a0,1,2 or 3-form
       sll_int32,  intent( in    )                    :: component     !< Specify the component of the function
       sll_real64, intent(   out )                    :: coefs_dofs(:) !< Coefficients of the projection.
       procedure(sll_i_function_3d_real64)            :: func1         !< Function to be projected.
       procedure(sll_i_function_3d_real64), optional  :: func2         !< Function to be projected.
       procedure(sll_i_function_3d_real64), optional  :: func3         !< Function to be projected.
       
     end subroutine update_dofs_function
  end interface



  !---------------------------------------------------------------------------!
  abstract interface
     function norm_squared( self, coefs, form, component ) result( r )
       use sll_m_working_precision
       import sll_c_maxwell_3d_base
       class( sll_c_maxwell_3d_base) :: self      !< Maxwell solver object.
       sll_real64                    :: coefs(:)  !< Values of the coefficient vectors for each DoF
       sll_int32                     :: form      !< Specify 0,1,2 or 3-form
       sll_int32                     :: component !< Specify the component of the form
       sll_real64                    :: r         !< Result: squared L2 norm
     end function norm_squared
  end interface



  !---------------------------------------------------------------------------!
  abstract interface
     function inner_product( self, coefs1, coefs2, form, component ) result( r )
       use sll_m_working_precision
       import sll_c_maxwell_3d_base
       class( sll_c_maxwell_3d_base) :: self      !< Maxwell solver object.
       sll_real64                    :: coefs1(:) !< Values of the coefficient vectors for each DoF
       sll_real64                    :: coefs2(:) !< Values of the coefficient vectors for each Do
       sll_int32                     :: form           !< Specify 0,1,2 or 3-form
       sll_int32, optional           :: component      !< Specify the component of the form
       sll_real64                    :: r
       
     end function inner_product
  end interface



  !---------------------------------------------------------------------------!
  abstract interface
     subroutine empty( self ) 
       import sll_c_maxwell_3d_base
       class( sll_c_maxwell_3d_base) :: self !< Maxwell solver object.

     end subroutine empty
  end interface



  !---------------------------------------------------------------------------!
  abstract interface 
     subroutine multiply_mass(  self, deg, coefs_in, coefs_out )
       use sll_m_working_precision
       import sll_c_maxwell_3d_base     
       class(sll_c_maxwell_3d_base) :: self
       sll_int32,  intent( in    )  :: deg(:)
       sll_real64, intent( in    )  :: coefs_in(:)
       sll_real64, intent(   out )  :: coefs_out(:)
       
     end subroutine multiply_mass
  end interface
  !---------------------------------------------------------------------------!
  abstract interface 
     subroutine multiply_mass_inverse(  self, form, coefs_in, coefs_out )
       use sll_m_working_precision
       import sll_c_maxwell_3d_base     
       class(sll_c_maxwell_3d_base) :: self
       sll_int32,  intent( in    )  :: form
       sll_real64, intent( in    )  :: coefs_in(:)
       sll_real64, intent(   out )  :: coefs_out(:)
       
     end subroutine multiply_mass_inverse
  end interface
  
contains

  subroutine compute_curl_part( self, delta_t, efield, bfield, betar )
    class(sll_c_maxwell_3d_base) :: self
    sll_real64, intent(in)     :: delta_t   !< Time step
    sll_real64, intent(inout)  :: efield(:)  !< Ey
    sll_real64, intent(inout)  :: bfield(:)  !< Bz
    sll_real64, optional       :: betar      !< 1/beta


    call self%compute_B_from_E( &
         delta_t, efield, bfield )

    call self%compute_E_from_B(&
         delta_t, bfield*betar , efield )

  end subroutine compute_curl_part


end module sll_m_maxwell_3d_base
