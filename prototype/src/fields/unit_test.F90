program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
  use numeric_constants
  use geometry_functions
  use sll_scalar_field_2d
  !use sll_mapped_mesh_base
  use sll_mapped_meshes
  implicit none
  
  type(sll_mapped_mesh_2d_analytic), target :: mesh
  type(scalar_field_2d)                     :: field
  class(sll_mapped_mesh_2d_base), pointer   :: m
  sll_int32 :: nc1, nc2
  procedure(polar_x1), pointer :: px1, px2, pjac11, pjac12, pjac21, pjac22


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

  call initialize_scalar_field_2d( &
       field, &
       "px1_field", &
       m, &
       NODE_CENTERED_FIELD, &
       px1)

  print*, m%x1_at_node(5,3), m%x1(.3_f64, .4_f64)

  call write_scalar_field_2d( field, &
                              multiply_by_jacobian=.true. )

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
