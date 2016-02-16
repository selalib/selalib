program test_sparse_matrix_mumps

#include "sll_working_precision.h"
#include "sll_memory.h"

use sll_m_sparse_matrix
use sll_m_collective

implicit none

type(sll_t_csr_matrix),     pointer     :: mat
sll_int32                               :: num_rows
sll_int32                               :: num_elements
sll_int32,  dimension(:,:), allocatable :: local_to_global_row
sll_real64, dimension(:),   allocatable :: b
sll_real64, dimension(:),   allocatable :: x
sll_int32                               :: num_local_dof_row
sll_int32                               :: ierr
sll_int32                               :: i
sll_real64                              :: val
sll_int32                               :: prank
sll_int32                               :: psize

call sll_s_boot_collective()
prank = sll_f_get_collective_rank( sll_v_world_collective )
psize = sll_f_get_collective_size( sll_v_world_collective )
write(*,*) " #", prank, " of ", psize

num_elements      = 2000
num_rows          = num_elements
num_local_dof_row = 1

SLL_ALLOCATE(local_to_global_row(num_local_dof_row,num_elements),ierr)
SLL_ALLOCATE(b(num_rows),ierr)
SLL_ALLOCATE(x(num_rows),ierr)

local_to_global_row(1:num_local_dof_row,1:num_elements) = 0
do i=1,num_elements
  local_to_global_row(1,i)=i
enddo

b = 1._f64

print *,'#begin sll_f_new_csr_matrix with mumps solver'

mat => sll_f_new_csr_matrix(  num_rows,            &
                              num_rows,            &
                              num_elements,        &
                              local_to_global_row, &
                              num_local_dof_row,   &
                              local_to_global_row, &
                              num_local_dof_row,   &
                              sll_p_mumps )

print *,'#end sll_f_new_csr_matrix with mumps solver'

!do identity matrix
print *,'#begin add_to_csr_matrix'
val = 1._f64
do i=1,num_rows
  call sll_s_add_to_csr_matrix(mat, val, i, i)
enddo
print *,'#end add_to_csr_matrix'

call sll_s_factorize_csr_matrix(mat)

x = 0.0_f64
print *,'#begin sll_s_solve_csr_matrix with mumps'
call sll_s_solve_csr_matrix(mat, b, x)
print *,'#end sll_s_solve_csr_matrix with mumps'
  
if(maxval(abs(x-1.0_f64))<1.e-13)then
  print *,'#PASSED'
else
  print *,'#system is not good resolved',maxval(abs(x))
endif

call sll_s_free_csr_matrix(mat)

call sll_s_halt_collective()

end program test_sparse_matrix_mumps
