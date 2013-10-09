!     
! File:   sparse_matrix_LM.F90
! Author: root
!
! Created on January 18, 2012, 10:33 AM
!

module SparseMatrix_LM_Module
    use sparse_matrix_def
    use qsort_c_module
    implicit none

    private
    public :: create_SparseMatrix_with_LM

contains

    !---------------------------------------------------------------------------------
    subroutine create_SparseMatrix_with_LM(this, ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, ai_COEF)
        implicit none
        !> param[inout] this : CSR MATRIX STRUCTURE
        type(csr_matrix) :: this
        !> param[in] ai_nC : NUMBER OF COLUMNS, IT IS THE DIMENSION OF THE 1st SPACE
        integer :: ai_nC
        !> param[in] ai_nR : NUMBER OF ROWS, IT IS THE DIMENSION OF THE 2nd SPACE
        integer :: ai_nR
        !> param[in] ai_nel : NUMBER OF ELEMENTS (WITH NON ZERO MEASURE) IN THE PATCH
        integer :: ai_nel
        !> param[in] api_LM_1 : LM ARRAY FOR ROWS
        integer, dimension(:,:), pointer :: api_LM_1
        !> param[in] api_LM_1 : LM ARRAY FOR COLUMNS
        integer, dimension(:,:), pointer :: api_LM_2
        !> param[in] ai_nen_1 : NUMBER OF NON VANISHING FUNCTIONS PER ELEMENT, IN THE 1st SPACE
        integer :: ai_nen_1
        !> param[in] ai_nen_2 : NUMBER OF NON VANISHING FUNCTIONS PER ELEMENT, IN THE 2nd SPACE
        integer :: ai_nen_2
        !>
        integer, optional :: ai_COEF
        !local var
        integer :: li_err, li_flag
        integer :: li_nnz
        integer, dimension(:,:), pointer :: lpi_columns
        integer, dimension(:), pointer :: lpi_occ
        integer :: li_COEF

        li_COEF = 10
        if (present(ai_COEF)) then
            li_COEF = ai_COEF
        end if

        allocate(lpi_columns(ai_nR, 0:li_COEF * ai_nen_2))

        allocate(lpi_occ(ai_nR + 1))
        lpi_columns(:,:) = 0
        lpi_occ(:) = 0

        ! COUNTING NON ZERO ELEMENTS

        li_nnz = count_non_zero_elts(ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, lpi_columns, lpi_occ)


        this % ol_use_mm_format = .FALSE.

        this % oi_nR = ai_nR
        this % oi_nC = ai_nC
        this % oi_nel = li_nnz


        allocate(this % opi_ia(this % oi_nR + 1))


        allocate(this % opi_ja(this % oi_nel))


        allocate(this % opr_a(this % oi_nel))


        call init_SparseMatrix(this, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, lpi_columns, lpi_occ)

        this % opr_a(:) = 0.0_wp

        deallocate(lpi_columns)
        deallocate(lpi_occ)
    end subroutine create_SparseMatrix_with_LM
    !---------------------------------------------------------------------------------
    subroutine init_SparseMatrix(self, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, api_columns, api_occ)
        ! _1 FOR ROWS
        ! _2 FOR COLUMNS
        implicit none
        type(csr_matrix) :: self
        integer, dimension(:,:), pointer :: api_LM_1, api_LM_2
        integer :: ai_nel, ai_nen_1, ai_nen_2
        integer, dimension(:,:), pointer :: api_columns
        integer, dimension(:), pointer :: api_occ
        !local var
        integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_index, li_i, li_size
        integer :: li_err, li_flag
        real(wp), dimension(:), pointer :: lpr_tmp

        ! INITIALIZING ia
        self % opi_ia(1) = 1

        do li_i = 1, self % oi_nR

            self % opi_ia(li_i + 1) = self % opi_ia(1) + SUM(api_occ(1: li_i))

        end do

        ! INITIALIZING ja
        do li_e = 1, ai_nel

            do li_b_1 = 1, ai_nen_1

                li_A_1 = api_LM_1(li_b_1, li_e)

                if (li_A_1 == 0) then
                    cycle
                end if

                if (api_columns(li_A_1, 0) == 0) then
                    cycle
                end if

                li_size = api_columns(li_A_1, 0)

                allocate ( lpr_tmp(li_size), stat = li_err)
                if (li_err .ne. 0) li_flag = 10

                lpr_tmp(1: li_size) = real( api_columns(li_A_1, 1: li_size))

                call QsortC(lpr_tmp)

                do li_i = 1, li_size

                    self % opi_ja(self % opi_ia(li_A_1) + li_i - 1) = int ( lpr_tmp(li_i))

                end do

                api_columns(li_A_1, 0) = 0
                deallocate ( lpr_tmp)

            end do

        end do

    end subroutine init_SparseMatrix
    !---------------------------------------------------------------------------------
    integer function count_non_zero_elts(ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, api_columns, api_occ)
        ! _1 FOR ROWS
        ! _2 FOR COLUMNS
        implicit none
        integer :: ai_nR, ai_nC
        integer, dimension(:,:), pointer :: api_LM_1, api_LM_2
        integer :: ai_nel, ai_nen_1, ai_nen_2
        integer, dimension(:,:), pointer :: api_columns
        integer, dimension(:), pointer :: api_occ
        !local var
        integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_i
        integer :: li_err, li_flag
        real(wp), dimension(:), pointer :: lpr_tmp
        integer :: li_result
        integer, dimension(2) :: lpi_size
        logical :: ll_done
        integer, dimension(:,:), pointer :: lpi_columns

        ! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED
        do li_e = 1, ai_nel

            do li_b_1 = 1, ai_nen_1

                li_A_1 = api_LM_1(li_b_1, li_e)
                if (li_A_1 == 0) then
                    cycle
                end if

                do li_b_2 = 1, ai_nen_2

                    li_A_2 = api_LM_2(li_b_2, li_e)
                    if (li_A_2 == 0) then
                        cycle
                    end if

                    ll_done = .false.
                    ! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_1, li_A_2)
                    do li_i = 1, api_columns(li_A_1, 0)

                        if (api_columns(li_A_1, li_i) /= li_A_2) then
                            cycle
                        end if

                        ll_done = .true.
                        exit

                    end do

                    if (.not.ll_done) then

                        api_occ(li_A_1) = api_occ(li_A_1) + 1

                        ! li_A_1 IS THE ROW NUM, li_A_2 THE COLUMN NUM
                        ! INITIALIZATION OF THE SPARSE MATRIX
                        api_columns(li_A_1, 0) = api_columns(li_A_1, 0) + 1
                        api_columns(li_A_1, api_columns(li_A_1, 0)) = li_A_2

                        ! resizing the array
                        lpi_size(1) = SIZE(api_columns, 1)
                        lpi_size(2) = SIZE(api_columns, 2)
                        if (lpi_size(2) < api_columns(li_A_1, 0)) then
                            ALLOCATE(lpi_columns(lpi_size(1), lpi_size(2)))
                            lpi_columns = api_columns

                            DEALLOCATE(api_columns)

                            ALLOCATE(api_columns(lpi_size(1), 2 * lpi_size(2)))
                            api_columns(1:lpi_size(1), 1:lpi_size(2)) = lpi_columns(1:lpi_size(1), 1:lpi_size(2))

                            DEALLOCATE(lpi_columns)
                        end if


                    end if

                end do

            end do

        end do

        ! COUNT NON ZERO ELEMENTS
        li_result = SUM(api_occ(1: ai_nR))

        count_non_zero_elts = li_result
    end function count_non_zero_elts

end module SparseMatrix_LM_Module

