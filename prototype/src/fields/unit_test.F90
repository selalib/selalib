program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_constants
  use sll_common_coordinate_transformations
  use sll_module_interpolators_1d_base
  use sll_cubic_spline_interpolator_1d
  use sll_scalar_field_2d
 ! use sll_module_mapped_meshes_2d_base
 ! use sll_module_mapped_meshes_2d
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_scalar_field_initializers_base
  use sll_landau_2d_initializer
  implicit none
  
  type(sll_coordinate_transformation_2d_analytic), pointer :: transf
  type(sll_logical_mesh_2d), pointer :: ml
  type(scalar_field_2d)                     :: field
  class(sll_coordinate_transformation_2d_base), pointer   :: m
  sll_int32 :: nc1, nc2, iplot
  procedure(transformation_func_nopass), pointer :: px1, px2, pjac11, pjac12, &
       pjac21, pjac22
  type(init_landau_2d), target :: init_landau
  class(scalar_field_2d_initializer_base), pointer    :: pfinit
  type(cubic_spline_1d_interpolator), target  :: interp_eta1
  type(cubic_spline_1d_interpolator), target  :: interp_eta2
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr
  sll_int32 :: i,j
  sll_real64, dimension(4) :: sinprod_params 

  nc1 = 10
  nc2 = 10
  ml => new_logical_mesh_2d( nc1, nc2)
  px1 => sinprod_x1
  px2 => sinprod_x2
  pjac11 => sinprod_jac11
  pjac12 => sinprod_jac12
  pjac21 => sinprod_jac21
  pjac22 => sinprod_jac22
  sinprod_params(:) = (/0.1_f64, 0.1_f64,1.0_f64,1.0_f64/)

  transf => new_coordinate_transformation_2d_analytic(&
       'sinprod',&
       ml,  &
       px1, &
       px2, &
       pjac11, &
       pjac12, &
       pjac21, &
       pjac22, &
       sinprod_params )
  m => transf

  call init_landau%initialize(m,NODE_CENTERED_FIELD,0.001_f64)
  pfinit => init_landau

  ! Set up the interpolators for the field
  call interp_eta1%initialize( nc1+1, 0.0_f64, 1.0_f64, SLL_PERIODIC )
  call interp_eta2%initialize( nc2+1, 0.0_f64, 1.0_f64, SLL_PERIODIC )
  interp_eta1_ptr => interp_eta1
  interp_eta2_ptr => interp_eta2
 
  call initialize_scalar_field_2d( &
       field, &
       "px1_field", &
       m, &
       NODE_CENTERED_FIELD, &
       interp_eta1_ptr, &
       interp_eta2_ptr, &
       pfinit )

  print*, m%x1_at_node(5,3), m%x1(.3_f64, .4_f64)

  

  do iplot = 1, 10
     do i=1, nc1+1
        do j=1,nc2+1
           field%data(i,j) = exp(-(transf%x1_at_node(i,j)**2+transf%x2_at_node(i,j)**2)*iplot*0.1)
        end do
     end do
     call write_scalar_field_2d( field, multiply_by_jacobian=.true. )
     call write_scalar_field_2d( field, multiply_by_jacobian=.true., output_file_name="field" )
  end do

  call transf%write_to_file()

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
