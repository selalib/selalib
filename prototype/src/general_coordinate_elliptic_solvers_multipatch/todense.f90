subroutine todense_csrMatrix(csr_mat, apr_dense)
  use SparseMatrix_Module
  implicit none
  type(csr_matrix) :: csr_mat
  real(8), dimension(:,:),INTENT(INOUT) :: apr_dense		
  !local var
  integer  :: li_i,li_j,li_k	        
  
  !initialisation
  print*,'csr_mat%oi_nR=', csr_mat%oi_nR
  print*,'csr_mat%oi_nC=', csr_mat%oi_nC
  !print*, apr_dense(1,1)
  !apr_dense ( 1 : csr_mat%oi_nR , 1 : csr_mat%oi_nC ) = 0.0_8
  print*, '?'
  do li_i =1,csr_mat%oi_nR 
     
     do li_k = csr_mat%opi_ia ( li_i ) , csr_mat%opi_ia ( li_i + 1 ) - 1
        
        li_j = csr_mat%opi_ja(li_k)
        
        print*, '?'

        apr_dense ( li_i , li_j ) = csr_mat%opr_a ( li_k )		
				
     end do
			
  end do

  print*, '??'
				
end subroutine todense_csrMatrix
