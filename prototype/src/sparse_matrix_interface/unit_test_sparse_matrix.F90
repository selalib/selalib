program unit_test_sparse_matrix
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_sparse_matrix_module
  use sll_collective
  implicit none
  
  !sll_int32 :: colptr(6)
  !sll_int32 :: row(6)
  !sll_real64 :: avals(9)
  type(sll_csr_matrix), pointer :: mat
  sll_int32 :: num_rows
  sll_int32 :: num_elements
  sll_int32, dimension(:,:), allocatable :: local_to_global_row
  sll_real64, dimension(:), allocatable :: b
  sll_real64, dimension(:), allocatable :: x
  sll_int32 :: num_local_dof_row
  sll_int32 :: ierr
  sll_int32 :: i
  sll_real64 :: val

  call sll_boot_collective()
  
  num_rows = 32
  num_elements = num_rows
  num_local_dof_row = 1
  
  SLL_ALLOCATE(local_to_global_row(num_local_dof_row,num_elements),ierr)
  SLL_ALLOCATE(b(num_rows),ierr)
  SLL_ALLOCATE(x(num_rows),ierr)

  local_to_global_row(1:num_local_dof_row,1:num_elements) = 0
  do i=1,num_elements
    local_to_global_row(1,i)=i
  enddo

  b = 1._f64

  !this is the call in the elliptic solver
!   es%sll_csr_mat => new_csr_matrix( &
!      solution_size, &
!      solution_size, &
!      num_cells_eta1*num_cells_eta2, &
!      es%local_to_global_spline_indices, &
!      es%total_num_splines_loc, &
!      es%local_to_global_spline_indices, &
!      es%total_num_splines_loc)

  
  print *,'#begin new_csr_matrix'
  mat => new_csr_matrix( &
    num_rows, &
    num_rows, &
    num_elements, &
    local_to_global_row, &
    num_local_dof_row, &
    local_to_global_row, &
    num_local_dof_row)
  print *,'#end new_csr_matrix'
  
  !do identity matrix
  print *,'#begin add_to_csr_matrix'
  val = 1._f64
  do i=1,num_rows
    call sll_add_to_csr_matrix(mat, val, i, i)
  enddo
  !print *,'#colptr=',mat%linear_solver%colptr
  print *,'#end add_to_csr_matrix'
  print *,'#begin sll_factorize_csr_matrix'
  call sll_factorize_csr_matrix(mat)
  print *,'#end sll_factorize_csr_matrix'
  print *,'#begin sll_solve_csr_matrix'
  call sll_solve_csr_matrix(mat, b, x)
  print *,'#end sll_solve_csr_matrix'
  
  
  !colptr(1:6)    = (/1,3,5,7,9,10/)
  !row(1:9)    = (/1,2,2,3,3,4,4,5,5/)
  !avals(1:9) = (/2._f64,-1._f64,2._f64,-1._f64,2._f64,-1._f64,2._f64,-1._f64,2._f64/)

  !print *,'#x=',x
  call sll_halt_collective()
  if(maxval(abs(x-1))<1.e-13)then
    print *,'#PASSED'
  else
    print *,'#system is not good resolved',x
  endif


end program
