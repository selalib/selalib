!     
! File:   sparse_matrix_basic.F90
! Author: root
!
! Created on January 18, 2012, 10:36 AM
!

module SparseMatrix_Basic_Module
    use used_precision
    use sparse_matrix_def
    implicit none

contains
    !---------------------------------------------------------------------------------
    subroutine create_SparseMatrix(self, ai_nR, ai_nC, ai_nel)
        implicit none
        !> param[inout] self : CSR MATRIX STRUCTURE
        type(csr_matrix) :: self
        !> param[in] ai_nC : NUMBER OF COLUMNS
        integer :: ai_nC
        !> param[in] ai_nR : NUMBER OF ROWS
        integer :: ai_nR
        !> param[in] ai_nel : NUMBER OF NON ZERO ELEMENTS
        integer :: ai_nel
        !local var

        print*, 'tet'
        self % ol_use_mm_format = .FALSE.

        self % oi_nR = ai_nR
        self % oi_nC = ai_nC
        self % oi_nel = ai_nel

        allocate(self % opi_ia(self % oi_nR + 1))

        allocate(self % opi_ja(self % oi_nel))

        allocate(self % opr_a(self % oi_nel))

        self % opr_a(:) = 0.0_wp
    end subroutine create_SparseMatrix
    !---------------------------------------------------------------------------------
    ! self ROUTINE CREATES A NEW FULL MATRIX FROM A COMPRESSED SYMMETRIC MATRIX
    subroutine create_SparseMatrix_from_sym(self, ao_A)
        implicit none
        !> param[inout] self : CSR MATRIX STRUCTURE
        type(csr_matrix) :: self
        !> param[in] ao_A : CSR MATRIX STRUCTURE
        type(csr_matrix) :: ao_A
        !local var
        integer :: li_i, li_k,li_count, li_i_tmp, li_index
        integer, dimension(ao_A % oi_nR + 1) :: lpi_occ

        lpi_occ(:) = 0

        ! COMPUTING THE NUMBER OF NON ZERO DIAGONAL ELEMENT
        li_count = 0
        do li_i = 1, ao_A % oi_nR
            if (ao_A % opi_ja(ao_A % opi_ia(li_i)) == li_i) then
                li_count = li_count + 1
            end if
        end do

        self % ol_use_mm_format = .FALSE.

        self % oi_nel = 2 * ao_A % oi_nel - li_count
        self % oi_nR = ao_A % oi_nR
        self % oi_nC = ao_A % oi_nC

        allocate(self % opi_ia(self % oi_nR + 1))

        allocate(self % opi_ja(self % oi_nel))

        allocate(self % opr_a(self % oi_nel))

        ! COPY THE OLD MATRIX IN THE NEW MATRIX
        do li_i = 1, ao_A % oi_nR
            lpi_occ(li_i) = ao_A % opi_ia(li_i + 1) - ao_A % opi_ia(li_i)
            do li_i_tmp = 1, li_i - 1
                li_k = ao_A % opi_ia(li_i_tmp)
                do while ((li_k <= ao_A % opi_ia(li_i_tmp + 1) - 1) .AND. (ao_A % opi_ja(li_k) <= li_i))
                    if (li_i == ao_A % opi_ja(li_k)) then
                        lpi_occ(li_i) = lpi_occ(li_i) + 1
                    end if
                    li_k = li_k + 1
                end do
            end do
        end do

        self % opi_ia(1) = 1

        do li_i = 1, self % oi_nR

            self % opi_ia(li_i + 1) = self % opi_ia(1) + SUM(lpi_occ(1: li_i))

        end do

        li_index = 1
        do li_i = 1, ao_A % oi_nR
            do li_i_tmp = 1, li_i - 1
                li_k = ao_A % opi_ia(li_i_tmp)
                do while ((li_k <= ao_A % opi_ia(li_i_tmp + 1) - 1) .AND. (ao_A % opi_ja(li_k) <= li_i))
                    if (li_i == ao_A % opi_ja(li_k)) then
                        self % opi_ja(li_index) = li_i_tmp
                        li_index = li_index + 1
                    end if
                    li_k = li_k + 1
                end do
            end do

            self % opi_ja(li_index: li_index + ao_A % opi_ia(li_i + 1) - 1 - ao_A % opi_ia(li_i)) = &
            ao_A % opi_ja(ao_A % opi_ia(li_i): ao_A % opi_ia(li_i + 1) - 1)

            li_index = li_index + ao_A % opi_ia(li_i + 1) - ao_A % opi_ia(li_i)
        end do

    end subroutine create_SparseMatrix_from_sym
    !---------------------------------------------------------------------------------
    subroutine free_SparseMatrix(self)
        implicit none
        type(csr_matrix) :: self

        deallocate(self % opi_ia)
        deallocate(self % opi_ja)
        deallocate(self % opr_a)

        if (self % ol_use_mm_format) then

            DEALLOCATE ( self % opi_i)

        end if

    end subroutine free_SparseMatrix
    !---------------------------------------------------------------------------------
    !%%MatrixMarket matrix coordinate real general
    !%=================================================================================
    !%
    !% self ASCII file represents a sparse MxN matrix with L
    !% nonzeros in the following Matrix Market format:
    !%
    !% +----------------------------------------------+
    !% |%%MatrixMarket matrix coordinate real general | <--- header line
    !% |%                                             | <--+
    !% |% comments                                    |    |-- 0 or more comment lines
    !% |%                                             | <--+
    !% |    M  N  L                                   | <--- rows, columns, entries
    !% |    I1  J1  A(I1, J1)                         | <--+
    !% |    I2  J2  A(I2, J2)                         |    |
    !% |    I3  J3  A(I3, J3)                         |    |-- L lines
    !% |        . . .                                 |    |
    !% |    IL JL  A(IL, JL)                          | <--+
    !% +----------------------------------------------+
    !%
    !% Indices are 1-based, i.e. A(1,1) is the first element.
    !%
    !%=================================================================================
    subroutine print_csrMatrix_mmFormat(self, as_file)
        implicit none
        type(CSR_MATRIX) :: self
        character(len = *), intent(in) :: as_file
        !local var
        integer :: li_i, li_j, li_k
        integer :: li_file = 1
        integer :: li_ios

        open(unit = li_file, file = as_file // '.mm', iostat = li_ios)
        if (li_ios /= 0) STOP "print_csrMatrix_mmFormat : erreur d'ouverture du fichier "

        !write(li_file, *) '%%MatrixMarket matrix coordinate real general'
        !write(li_file, *) self % oi_nR, ',', self % oi_nC, ',', self % oi_nel
        do li_i = 1, self % oi_nR

            do li_k = self % opi_ia(li_i), self % opi_ia(li_i + 1) - 1

                li_j = self % opi_ja(li_k)

                write(li_file, *) self % opr_a(li_k)!li_i, ',', li_j, ',', self % opr_a(li_k)

            end do

        end do

        close(li_file)

    end subroutine print_csrMatrix_mmFormat
    !---------------------------------------------------------------------------------
    subroutine reset_SparseMatrix(self)
        implicit none
        type(csr_matrix) :: self

        self % opr_a = 0.0_wp

    end subroutine reset_SparseMatrix
end module SparseMatrix_Basic_Module

