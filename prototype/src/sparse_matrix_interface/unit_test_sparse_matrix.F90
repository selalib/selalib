program unit_test_sparse_matrix
  use sll_sparse_matrix_module
  implicit none
  
  type(sll_csr_matrix), pointer :: mat
  
  
  !mat => new_csr_matrix()
  
  print *,'#PASSED'

end program