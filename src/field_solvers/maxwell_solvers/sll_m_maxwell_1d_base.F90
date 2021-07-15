!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations in 1D
!> @details
!> Contains the abstract class to create a Maxwell solver in 1D.
!> @author
!> Katharina Kormann
module sll_m_maxwell_1d_base
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_utilities, only: &
       sll_s_int2string

  implicit none

  public :: &
       sll_i_function_1d_real64, &
       sll_c_maxwell_1d_base, &
       sll_s_plot_two_fields_1d

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, abstract :: sll_c_maxwell_1d_base
     sll_real64 :: Lx          !< length of periodic domain
     sll_int32  :: n_cells     !< number of cells (and grid points)
     sll_int32  :: n_dofs0      !< number of degrees of freedom
     sll_int32  :: n_dofs1     !< number of degrees of freedom
     sll_real64 :: delta_x     !< cell size
     sll_int32  :: s_deg_0     !< spline degree 0-forms
     sll_int32  :: s_deg_1     !< spline degree 1-forms
     sll_real64 :: mass_solver_tolerance = 1d-12 !< tolerance for the mass solver
     sll_real64 :: poisson_solver_tolerance = 1d-12 !< tolerance for the Poisson solver
     logical :: strong_ampere = .false. !< flag to switch between strong and weak Ampere formulation 

     sll_real64 :: solver_tolerance = 1d-12 !< tolerance for the field solver

   contains
     procedure(compute_field1_from_field2), deferred :: &
          compute_e_from_b !< Solve E and B part of Ampere's law with B constant in time
     procedure(compute_field1_from_field2), deferred :: &
          compute_b_from_e !< Solve Faraday equation with E constant in time
     procedure(signature_compute_field_from_field), deferred :: &
          compute_e_from_rho !< Solve E from rho using Poisson
     procedure(signature_compute_field_from_field), deferred :: &
          compute_rho_from_e !< Compute rho from efield by Gauss law
     procedure(signature_compute_E_from_j_1d), deferred :: &
          compute_e_from_j !< Solve E from time integrated current (second part of Ampere's law)
     procedure(signature_compute_phi_from_field), deferred :: &
           compute_phi_from_rho !< Solve E from rho using Poisson
       procedure(signature_compute_phi_from_field), deferred :: &
          compute_phi_from_j !< Solve E from rho using Poisson
     !procedure(signature_solve), deferred :: &
     !     solve !< Solve Amperes law and Faraday equation
     procedure(update_dofs_function), deferred :: &
          compute_rhs_from_function !< Compute the right-hand-side for a given function f. For Galerkin it is the inner product with the basis functions. For Collocation it is simply a function evaluation at the grid points.
     procedure(update_dofs_function), deferred :: &
          L2projection !< L2 projection
     procedure(norm_squared), deferred :: &
          L2norm_squared !< Square of the L2norm
     procedure(signature_inner_product), deferred :: &
          inner_product !< inner product with spline mass matrix
     procedure(empty), deferred :: &
          free !< destructor
     procedure(signature_mass), deferred :: &
          multiply_mass !< product with mass matrix
     procedure(signature_mass), deferred :: &
          invert_mass !< product with inverse mass matrix
     procedure(multiply_derivative), deferred :: &
          multiply_g !< multiplication with gradient matrix
     procedure(multiply_derivative), deferred :: &
          multiply_gt !< multiplication with divergence matrix
     procedure :: compute_curl_part !< Solve source-free Maxwell's equations

     procedure :: transform_dofs !< By default this is just a copy of input to output

  end type sll_c_maxwell_1d_base

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine empty(self) 
       import sll_c_maxwell_1d_base
       class( sll_c_maxwell_1d_base)                    :: self !< Maxwell solver object.

     end subroutine empty
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     function norm_squared(self, coefs_dofs, degree) result( r )
       use sll_m_working_precision
       import sll_c_maxwell_1d_base
       class( sll_c_maxwell_1d_base)                    :: self !< Maxwell solver object.
       sll_real64                                     :: coefs_dofs(:) !< Values of the coefficient vectors for each DoF
       sll_int32                                      :: degree !< Degree of the basis function used for whcih the DoF-coefficients are given.
       sll_real64                                     :: r
     end function norm_squared
  end interface


  !---------------------------------------------------------------------------!
  abstract interface
     subroutine signature_mass(self, in, out, degree)
       use sll_m_working_precision
       import sll_c_maxwell_1d_base
       class( sll_c_maxwell_1d_base), intent(inout)   :: self !< Maxwell solver object.
       sll_real64, intent(in)                         :: in(:) !< Values of the coefficient vectors for each DoF
       sll_real64, intent(out)                        :: out(:) !< Values of the coefficient vectors for each DoF
       sll_int32, intent(in)                          :: degree !< Degree of the basis function used for whcih the DoF-coefficients are given.
     end subroutine signature_mass
  end interface




  !---------------------------------------------------------------------------!
  abstract interface
     function signature_inner_product(self, coefs1_dofs, coefs2_dofs, degree, degree2) result( r )
       use sll_m_working_precision
       import sll_c_maxwell_1d_base
       class( sll_c_maxwell_1d_base)                    :: self !< Maxwell solver object.
       sll_real64                                     :: coefs1_dofs(:) !< Values of the coefficient vectors for each DoF
       sll_real64                                     :: coefs2_dofs(:) !< Values of the coefficient vectors for each Do
       sll_int32                                      :: degree !< Degree of the basis function used for which the DoF-coefficients are given.
       sll_int32, optional                            :: degree2 !< Degree of the basis function used for which the second DoF-coefficients are given (if not present same as degree).
       sll_real64                                     :: r
     end function signature_inner_product
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     !> 1d real function
     function sll_i_function_1d_real64(x)
       use sll_m_working_precision ! can't pass a header file because the
       ! preprocessor prevents double inclusion.
       ! It is very rare.
       sll_real64             :: sll_i_function_1d_real64
       sll_real64, intent(in) :: x
     end function sll_i_function_1d_real64
  end interface
  !---------------------------------------------------------------------------!
  abstract interface
     subroutine update_dofs_function(self, func, degree, coefs_dofs)
       use sll_m_working_precision
       import sll_c_maxwell_1d_base
       import sll_i_function_1d_real64
       class( sll_c_maxwell_1d_base)    :: self !< Maxwell solver object.
       procedure(sll_i_function_1d_real64)  :: func !< Function to be projected.
       sll_int32, intent(in)          :: degree !< Degree of the basis function that should be used for projection.
       sll_real64, intent(out)        :: coefs_dofs(:) !< Coefficients of the projection.
     end subroutine update_dofs_function
  end interface

  !---------------------------------------------------------------------------!
  abstract interface 
     subroutine compute_field1_from_field2(self, delta_t, field_in, field_out)
       use sll_m_working_precision
       import sll_c_maxwell_1d_base     
       class(sll_c_maxwell_1d_base) :: self
       sll_real64, intent(in)     :: delta_t
       sll_real64, intent(in)     :: field_in(:)
       sll_real64, intent(inout)  :: field_out(:)
     end subroutine compute_field1_from_field2
  end interface

  !---------------------------------------------------------------------------!
  abstract interface 
     subroutine signature_compute_field_from_field(self, field_in, field_out)
       use sll_m_working_precision
       import sll_c_maxwell_1d_base     
       class(sll_c_maxwell_1d_base) :: self
       sll_real64, intent(in)     :: field_in(:)
       sll_real64, intent(out)  :: field_out(:)
     end subroutine signature_compute_field_from_field
  end interface

  !---------------------------------------------------------------------------!  
  abstract interface
     subroutine signature_compute_E_from_j_1d(self, current, component, E)
       use sll_m_working_precision
       import sll_c_maxwell_1d_base
       class(sll_c_maxwell_1d_base)             :: self
       sll_real64,dimension(:),intent(in)    :: current
       sll_int32, intent(in)                 :: component
       sll_real64,dimension(:),intent(inout) :: E
     end subroutine signature_compute_E_from_j_1d
  end interface

  !---------------------------------------------------------------------------!
  abstract interface 
     subroutine signature_compute_phi_from_field(self, in, phi, efield)
     use sll_m_working_precision
     import sll_c_maxwell_1d_base     
     class(sll_c_maxwell_1d_base) :: self
     sll_real64, intent(in)   :: in(:)
     sll_real64, intent(out)  :: phi(:)
     sll_real64, intent(out)  :: efield(:)
   end subroutine signature_compute_phi_from_field
  end interface

  !---------------------------------------------------------------------------!
  abstract interface 
     subroutine multiply_derivative(self, in, out)
       use sll_m_working_precision
       import sll_c_maxwell_1d_base     
       class(sll_c_maxwell_1d_base), intent(in) :: self
       sll_real64, intent(in)   :: in(:)
       sll_real64, intent(out)  :: out(:)
     end subroutine multiply_derivative
  end interface


contains
  !> write files to visualize 1d fields with gnuplot
  subroutine sll_s_plot_two_fields_1d(fname, n1, f1, f2, iplot, time )
    character(len=*),             intent(in) :: fname !< output file name
    sll_int32,                    intent(in) :: n1    !< size of f1 and f2 
    sll_real64, dimension(n1),    intent(in) :: f1    !< first field 2d
    sll_real64, dimension(n1),    intent(in) :: f2    !< second field 2d
    sll_int32,                    intent(in) :: iplot !< plot counter
    sll_real64,                   intent(in) :: time  !< step time

    integer          :: i
    character(len=4) :: cplot

    call sll_s_int2string(iplot, cplot)

    !write domains
    open( 80, file = fname//cplot//".dat" )
    do i=1,n1
       write(80,*) i, sngl(f1(i)), sngl(f2(i))
    end do
    close(80)

!!$    open( 90, file = fname//'plots.gnu', position="append" )
!!$    if ( iplot == 1 ) then
!!$       rewind(90)
!!$       !write(90,*)"set xr[-0.1:1.1]"
!!$       !write(90,*)"set yr[-0.1:1.1]"
!!$       !write(90,*)"set zr[-1.1:1.1]"
!!$       !write(90,*)"set cbrange[-1:1]"
!!$       !write(90,*)"set pm3d"
!!$       !write(90,*)"set surf"
!!$       write(90,*)"set term x11"
!!$    end if
!!$    write(90,*)"set title 'Time = ",time,"'"
!!$    write(90,"(a)",advance='no')"plot '"//fname//cplot//".dat' w lines"
!!$    write(90,"(a)")",'"//fname//cplot//".dat' u 1:3 w lines"
!!$    write(90,*)"pause -1"
!!$    close(90)

  end subroutine sll_s_plot_two_fields_1d

  !> Transformation of the dofs (identity by default)
  subroutine transform_dofs(self, in, out, degree)
     class(sll_c_maxwell_1d_base), intent(inout) :: self
     sll_int32, intent(in) :: degree !< spline degree
     sll_real64, intent(in) :: in(:)  !< spline coefficients of projection
     sll_real64, intent(out) :: out(:)  !< spline coefficients of projection
     ! local variables

     out = in
     
   end subroutine transform_dofs

   !> Solve source-free Maxwell's equations: Default implementation by solving first B_from_E then E_from_B
   subroutine compute_curl_part( self, delta_t, efield, bfield, betar )
     class(sll_c_maxwell_1d_base) :: self
     sll_real64, intent(in)     :: delta_t   !< Time step
     sll_real64, intent(inout)  :: efield(:)  !< Ey
     sll_real64, intent(inout)  :: bfield(:)  !< Bz
     sll_real64, optional       :: betar      !< 1/beta
     !local variables
     sll_real64 :: factor
     
     if( present(betar) ) then
        factor = betar
     else
        factor = 1._f64
     end if

     
     call self%compute_B_from_E( &
          delta_t, efield, bfield )
    
     call self%compute_E_from_B(&
          delta_t, bfield*factor , efield )
     
   end subroutine compute_curl_part

end module sll_m_maxwell_1d_base
