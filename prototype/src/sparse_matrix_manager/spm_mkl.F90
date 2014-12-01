!
!	SparseMatrix.F90
!	SparseMatrix
!
!	Created by ahmed ratnani on 02/01/10.
!	Copyright 2010 __MyCompanyName__. All rights reserved.
!
module SparseMatrix_Module
    use sparse_matrix_def
    use SparseMatrix_Basic_Module
    use SparseMatrix_LM_Module
    use SparseMatrix_con_Module
    use SparseMatrix_Composed_Module
    use SparseMatrix_Operations_Module
    implicit none

    interface create_csr
            module procedure create_SparseMatrix    &
            , create_SparseMatrix_with_LM   &
            , create_SparseMatrix_with_connectivity &
            , create_SparseMatrixProfil_From4
    end interface

    interface copy_csr
            module procedure copy_SparseMatrix_From4
    end interface

    interface free_csr
            module procedure free_SparseMatrix
    end interface

    interface reset_csr
            module procedure reset_SparseMatrix
    end interface

    interface Mult_MV
        module procedure Mult_CSR_Matrix_Vector
    end interface

    interface Mult_MScal
        module procedure Mult_CSR_Matrix_Scalar
    end interface

    interface Add_MScal
        module procedure Add_CSR_Matrix_Scalar
    end interface

    interface add_MVal
        module procedure add_to_csr_Matrix &
        , add_to_csr_Matrix_condition
    end interface

    interface add_csr_Matrix
        module procedure add_csr_Matrix
    end interface

    interface mult_csr_Matrix
        module procedure mult_csr_Matrix
    end interface

    interface eigenval_for_circulant
        module procedure eigenval_for_circulant
    end interface eigenval_for_circulant


!    interface todense
!        module procedure todense_csrMatrix
!    end interface
!
!    interface remove_zeros
!        module procedure remove_zeros
!    end interface
!
!    interface set_direct_access_i
!        module procedure set_direct_access_i
!    end interface

end module SparseMatrix_Module
