program test_sparse_matrix_umfpack

#include "sll_working_precision.h"
#include "sll_memory.h"

   use sll_m_sparse_matrix
   use sll_m_umfpack

   implicit none

   type(sll_t_csr_matrix), pointer     :: mat
   sll_int32                               :: num_rows
   sll_int32                               :: num_elements
   sll_int32, dimension(:, :), allocatable :: local_to_global_row
   sll_real64, dimension(:), allocatable :: b
   sll_real64, dimension(:), allocatable :: x
   sll_int32                               :: num_local_dof_row
   sll_int32                               :: ierr
   sll_int32                               :: i
   sll_real64                              :: val

   num_elements = 2000
   num_rows = num_elements
   num_local_dof_row = 1

   SLL_ALLOCATE(local_to_global_row(num_local_dof_row, num_elements), ierr)
   SLL_ALLOCATE(b(num_rows), ierr)
   SLL_ALLOCATE(x(num_rows), ierr)

   local_to_global_row(1:num_local_dof_row, 1:num_elements) = 0
   do i = 1, num_elements
      local_to_global_row(1, i) = i
   end do

   b = 1._f64

   print *, '#begin sll_f_new_csr_matrix with umfpack'

   mat => sll_f_new_csr_matrix(num_rows, &
                               num_rows, &
                               num_elements, &
                               local_to_global_row, &
                               num_local_dof_row, &
                               local_to_global_row, &
                               num_local_dof_row, &
                               sll_p_umfpack)

   print *, '#end sll_f_new_csr_matrix with umfpack'

!do identity matrix
   print *, '#begin add_to_csr_matrix with umfpack'
   val = 1._f64
   do i = 1, num_rows
      call sll_s_add_to_csr_matrix(mat, val, i, i)
   end do
   print *, '#end add_to_csr_matrix with umfpack'

   call sll_s_factorize_csr_matrix(mat)
   b = 1._f64
   x = 0.0_f64
   print *, '#begin sll_s_solve_csr_matrix with umfpack'
   call sll_s_solve_csr_matrix(mat, b, x)
   print *, '#end sll_s_solve_csr_matrix with umfpack'

   if (maxval(abs(x - 1.0_f64)) < 1.e-13) then
      print *, '#PASSED'
   else
      print *, '#system is not good resolved with umfpack', maxval(abs(x))
   end if

   call sll_s_free_csr_matrix(mat)

end program test_sparse_matrix_umfpack
