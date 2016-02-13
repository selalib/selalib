program unit_test_sparse_matrix

#include "sll_working_precision.h"
#include "sll_memory.h"

use sll_m_sparse_matrix
#ifdef PASTIX
use sll_m_collective
#endif /* PASTIX */

implicit none

type(sll_t_csr_matrix),     pointer     :: mat
#ifdef UMFPACK
type(sll_t_csr_matrix),     pointer     :: mat2
#endif /* UMFPACK */
sll_int32                               :: num_rows
sll_int32                               :: num_elements
sll_int32,  dimension(:,:), allocatable :: local_to_global_row
sll_real64, dimension(:),   allocatable :: b
sll_real64, dimension(:),   allocatable :: x
sll_int32                               :: num_local_dof_row
sll_int32                               :: ierr
sll_int32                               :: i
sll_real64                              :: val

#ifdef PASTIX

sll_int32                               :: prank
sll_int32                               :: psize

call sll_s_boot_collective()
prank = sll_f_get_collective_rank( sll_v_world_collective )
psize = sll_f_get_collective_size( sll_v_world_collective )
write(*,*) " #", prank, " of ", psize

#endif /* PASTIX */

num_rows = 2000
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

print *,'#begin sll_f_new_csr_matrix'

mat => sll_f_new_csr_matrix(  num_rows,            &
                              num_rows,            &
                              num_elements,        &
                              local_to_global_row, &
                              num_local_dof_row,   &
                              local_to_global_row, &
                              num_local_dof_row)

print *,'#end sll_f_new_csr_matrix'

!do identity matrix
print *,'#begin add_to_csr_matrix'
val = 1._f64
do i=1,num_rows
  call sll_s_add_to_csr_matrix(mat, val, i, i)
enddo
print *,'#end add_to_csr_matrix'

x = 0.0_f64
print *,'#begin sll_s_solve_csr_matrix'
call sll_s_solve_csr_matrix(mat, b, x)
print *,'#end sll_s_solve_csr_matrix'
  
if(maxval(abs(x-1.0_f64))<1.e-13)then
  print *,'#PASSED'
else
  print *,'#system is not good resolved',maxval(abs(x))
endif

call sll_s_delete_csr_matrix(mat)

#ifdef PASTIX
call sll_s_halt_collective()
#endif /* PASTIX */

#ifdef UMFPACK

print *,'#begin sll_f_new_csr_matrix with umfpack'

mat2 => sll_f_new_csr_matrix( num_rows,            &
                              num_rows,            &
                              num_elements,        &
                              local_to_global_row, &
                              num_local_dof_row,   &
                              local_to_global_row, &
                              num_local_dof_row,   &
                              sll_p_umfpack )

print *,'#end sll_f_new_csr_matrix with umfpack'

!do identity matrix
print *,'#begin add_to_csr_matrix with umfpack'
val = 1._f64
do i=1,num_rows
  call sll_s_add_to_csr_matrix(mat2, val, i, i)
enddo
print *,'#end add_to_csr_matrix with umfpack'

call sll_s_factorize_csr_matrix(mat2)
b = 1._f64
x = 0.0_f64
print *,'#begin sll_s_solve_csr_matrix with umfpack'
call sll_s_solve_csr_matrix(mat2, b, x)
print *,'#end sll_s_solve_csr_matrix with umfpack'
  
if(maxval(abs(x-1.0_f64))<1.e-13)then
  print *,'#PASSED'
else
  print *,'#system is not good resolved with umfpack',maxval(abs(x))
endif

call sll_s_delete_csr_matrix(mat2)

#endif /* UMFPACK */

#ifdef PASTIX
call sll_s_halt_collective()
#endif /* PASTIX */

end program
