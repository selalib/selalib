!     
! File:   sparse_matrix_operations.F90
! Author: root
!
! Created on January 18, 2012, 11:37 AM
!
module SparseMatrix_Operations_Module
    use sparse_matrix_def
    use qsort_c_module
    use connectivities_def
    use connectivity_module
    implicit none

contains
    !---------------------------------------------------------------------------------
    subroutine Mult_CSR_Matrix_Vector(self, apr_x, apr_y)
        implicit none
        type(csr_matrix) :: self
        real(wp), dimension(:) :: apr_x, apr_y
        !local var
        integer :: li_i, li_k_1, li_k_2

        do li_i = 1, self % oi_nR

            li_k_1 = self % opi_ia(li_i)
            li_k_2 = self % opi_ia(li_i + 1) - 1

            apr_y(li_i) = DOT_PRODUCT(self % opr_a(li_k_1: li_k_2), apr_x(self % opi_ja(li_k_1: li_k_2)))
            
        end do

    end subroutine Mult_CSR_Matrix_Vector
    !---------------------------------------------------------------------------------
    subroutine add_to_csr_Matrix(self, ar_value, ai_A, ai_Aprime)
        implicit none
        type(csr_matrix) :: self
        real(wp) :: ar_value
        integer :: ai_A, ai_Aprime
        !local var
        integer :: li_result, li_i, li_j, li_k

        li_result = 0

        ! THE CURRENT LINE IS self%opi_ia(ai_A)
        do li_k = self % opi_ia(ai_A), self % opi_ia(ai_A + 1) - 1
            li_j = self % opi_ja(li_k)
            if (li_j == ai_Aprime) then
                self % opr_a(li_k) = self % opr_a(li_k) + ar_value
                exit
            end if
        end do

    end subroutine add_to_csr_Matrix
    !---------------------------------------------------------------------------------
    subroutine add_csr_Matrix(self, ao_A, ao_B)
        implicit none
        type(csr_matrix) :: self
        type(csr_matrix) :: ao_A
        type(csr_matrix) :: ao_B

        self % opr_a = ao_A % opr_a + ao_B % opr_a

    end subroutine add_csr_Matrix
    !---------------------------------------------------------------------------------
    subroutine add_to_csr_Matrix_condition(self, ar_value, ai_A, ai_Aprime, al_flag)
        implicit none
        type(csr_matrix) :: self
        real(wp) :: ar_value
        integer :: ai_A, ai_Aprime
        logical :: al_flag
        !local var
        integer :: li_result, li_i, li_j, li_k

        li_result = 0

        ! THE CURRENT LINE IS self%opi_ia(ai_A)
        do li_k = self % opi_ia(ai_A), self % opi_ia(ai_A + 1) - 1
            li_j = self % opi_ja(li_k)
            if ((li_j == ai_Aprime) .AND. (al_flag) .AND. (self % opr_a(li_k) == 0.0_wp)) then
                self % opr_a(li_k) = self % opr_a(li_k) + ar_value
                exit
            end if
        end do

    end subroutine add_to_csr_Matrix_condition
    !---------------------------------------------------------------------------------
    !> this routine computes the eigenvalues of the matrices M and K
    !> M and K ARE CIRCULANT  MATRICES
    !> FIRST WE COMPUTE THE COEF m_i OF THE FIRST LINE, AND COMPUTE THE EIGENVALUES
    subroutine eigenval_for_circulant(ao_K, apr_eigenK)
        implicit none
        type(CSR_MATRIX), intent ( in) :: ao_K
        real(wp), dimension(:), intent(inout) :: apr_eigenK
        ! LOCAL VARIABLES
        real(wp), dimension ( ao_K % oi_nR) :: lpr_K
        integer :: li_i
        integer :: li_j
        integer :: li_k
        integer :: li_N
        real(wp) :: lr_ck

!        CALL printlog("eigenval_for_circulant : Start", ai_dtllevel = mi_dtllevel_base + 1)

        lpr_K = 0.0_wp

        li_i = 1
        do li_k = ao_K % opi_ia(li_i), ao_K % opi_ia(li_i + 1) - 1

            li_j = ao_K % opi_ja(li_k)

            lpr_K(li_j) = ao_K % opr_a(li_k)

        end do

        !***********************************************************************************
        ! compute the k^th eigenvalue of K
        !***********************************************************************************
        li_N = ao_K % oi_nR
        apr_eigenK = 0.0_wp

        do li_k = 1, li_N

            lr_ck = 2.0_wp * pi * (li_k - 1) / li_N

            do li_j = 1, li_N

                apr_eigenK(li_k) = apr_eigenK(li_k) + lpr_K(li_j) * cos(lr_ck * (li_j - 1))

            end do

        end do
        !***********************************************************************************

!#ifdef _DEBUG
!        call concatmsg("K-row : ")
!        call concatmsg(lpr_K)
!        call printmsg(ai_dtllevel = DEBUG_LEVEL)
!
!        call concatmsg("eigenvalues : ")
!        call concatmsg(apr_eigenK)
!        call printmsg(ai_dtllevel = DEBUG_LEVEL)
!#endif
!
!        CALL printlog("eigenval_for_circulant : End", ai_dtllevel = mi_dtllevel_base + 1)

      end subroutine eigenval_for_circulant
end module SparseMatrix_Operations_Module
