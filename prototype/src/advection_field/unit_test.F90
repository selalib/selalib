program unit_test
#include "sll_mesh_types.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
  use numeric_constants
  use advection_field
  implicit none
  
  integer :: i,j
  sll_int32 :: ncells, ierr
  sll_real64 :: xmin, xmax, deltax, deltat, error
  sll_real64, dimension(:), pointer :: array1, array2, xout
  type(mesh_descriptor_1D), pointer :: m1D_descriptor
  type(field_1D_vec1), pointer :: field1D_scalar

  print*, 'checking compute_flow_1D'
  ncells = 10
  xmin = 0.0_f64
  xmax = 1.0_f64
  
  m1D_descriptor => new_mesh_descriptor_1D(xmin, xmax, ncells, PERIODIC)
  field1D_scalar => new_field_1D_vec1(m1D_descriptor)
  deltax = GET_FIELD_DELTA_X1( field1D_scalar )
  SLL_ALLOCATE(array1(ncells+1), ierr)
  SLL_ALLOCATE(array2(ncells+1), ierr)
  SLL_ALLOCATE(xout(ncells+1), ierr)
  print*,'testing constant advection to the right'
  do i = 1, ncells+1
     array1(i) =  1.0_f64 !xmin+ (i-1)*deltax ! 1.0_f64
     array2(i) = 0.0_f64
  end do
  deltat = 0.01_f64
  call compute_flow_1D_backward( array1, array2, 1.0_f64, deltat, xmin, ncells, deltax, PERIODIC, xout ) 
  error = 0.0_f64
  do i = 1, ncells+1
     error = max(error, abs(modulo(xmin+(i-1)*deltax - deltat, xmax-xmin) - xout(i)))
     !print*, i, modulo(xmin+(i-1)*deltax - deltat, xmax-xmin), xout(i)
  end do
  print*,'error=', error


  print *, 'Successful, exiting program.'
  
end program unit_test
