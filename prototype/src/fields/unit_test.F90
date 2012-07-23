program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
  use numeric_constants
  use geometry_functions
  use sll_module_interpolators_1d_base
  use sll_cubic_spline_interpolator_1d
  use sll_scalar_field_2d
  use sll_module_mapped_meshes_2d_base
  use sll_module_mapped_meshes_2d
  use sll_scalar_field_initializers_base
  use sll_landau_2d_initializer
  implicit none
  
  type(sll_mapped_mesh_2d_analytic), target :: mesh
  type(scalar_field_2d)                     :: field
  class(sll_mapped_mesh_2d_base), pointer   :: m
  sll_int32 :: nc1, nc2, iplot
  procedure(polar_x1), pointer :: px1, px2, pjac11, pjac12, pjac21, pjac22
  type(init_landau_2d), target :: init_landau
  class(scalar_field_2d_initializer_base), pointer    :: pfinit
  type(cubic_spline_1d_interpolator), target  :: interp_eta1
  type(cubic_spline_1d_interpolator), target  :: interp_eta2
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr


  nc1 = 10
  nc2 = 10
  px1 => sinprod_x1
  px2 => sinprod_x2
  pjac11 => sinprod_jac11
  pjac12 => sinprod_jac12
  pjac21 => sinprod_jac21
  pjac22 => sinprod_jac22
  call mesh%initialize(&
       'sinprod',&
       nc1, &
       nc2, &
       px1, &
       px2, &
       pjac11, &
       pjac12, &
       pjac21, &
       pjac22)
  m => mesh

  call init_landau%initialize(m,NODE_CENTERED_FIELD,0.001_f64)
  pfinit => init_landau

  ! Set up the interpolators for the field
  call interp_eta1%initialize( nc1+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
  call interp_eta2%initialize( nc2+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
  interp_eta1_ptr => interp_eta1
  interp_eta2_ptr => interp_eta2

  call initialize_scalar_field_2d( &
       field, &
       "px1_field", &
       m, &
       NODE_CENTERED_FIELD, &
       pfinit, &
       interp_eta1_ptr, &
       interp_eta2_ptr )

  print*, m%x1_at_node(5,3), m%x1(.3_f64, .4_f64)

  

  do iplot = 1, 10
     field%data = exp(-(mesh%x1_node**2+mesh%x2_node**2)*iplot*0.1)
     call write_scalar_field_2d( field, multiply_by_jacobian=.true. )
     call write_scalar_field_2d( field, multiply_by_jacobian=.true., output_file_name="field" )
  end do

  call mesh%write_to_file()

  print *, 'Successful, exiting program.'
  
!  contains
!
!  function init_function( eta1, eta2) result res
!     sll_real64, intent(in) :: eta1
!     sll_real64, intent(in) :: eta2
!     sll_real64             :: res
!   
!     res = exp(-(eta1*eta1+eta2+eta2))
!  end function init_function
  
end program unit_test
