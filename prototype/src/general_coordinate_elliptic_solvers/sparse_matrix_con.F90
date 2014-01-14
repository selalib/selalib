!     
! File:   sparse_matrix_con.F90
! Author: root
!
! Created on January 18, 2012, 10:34 AM
!

module SparseMatrix_Con_Module
    use sparse_matrix_def
    use qsort_c_module
    use connectivities_def
    use connectivity_module
    implicit none

    private
    public :: create_SparseMatrix_with_connectivity

contains
    !---------------------------------------------------------------------------------
    subroutine create_SparseMatrix_with_connectivity(self, ai_nR, ai_nC &
        , ai_npatch, api_nel    &
        , ao_con_C, ao_con_R    &
        , ai_COEF)
        implicit none
        !> param[inout] self : CSR MATRIX STRUCTURE
        type(csr_matrix) :: self
        !> param[in] ai_nC : NUMBER OF COLUMNS
        integer :: ai_nC
        !> param[in] ai_nR : NUMBER OF ROWS
        integer :: ai_nR
        !> param[in] ai_npatch : NUMBER OF PATCHS
        integer :: ai_npatch
        !> param[in] api_nel : NUMBER OF NON ZERO ELEMENTS IN EACH PATCH
        integer, dimension(:) :: api_nel
        !>
        TYPE(CONNECTIVITY) :: ao_con_C
        TYPE(CONNECTIVITY) :: ao_con_R
        !>
        integer, optional :: ai_COEF
        !local var
        !integer :: li_err, li_flag
        integer :: li_nnz
        integer, dimension(:,:), pointer :: lpi_columns
        integer, dimension(:), pointer :: lpi_occ
        integer :: li_COEF
        integer :: li_maxnen_C

        !> \todo il y a un bug ici, il faudra reduire li_COEF et voir
        li_COEF = 20
        if (present(ai_COEF)) then
            li_COEF = ai_COEF
        end if

        li_maxnen_C = MAXVAL(ao_con_C % opi_nen(:,:))

        allocate(lpi_columns(ai_nR, 0:li_COEF * li_maxnen_C))

        allocate(lpi_occ(ai_nR + 1))

        lpi_columns(:,:) = 0
        lpi_occ(:) = 0

    ! ...
    ! COUNTING NON ZERO ELEMENTS
    ! ...
        li_nnz = count_non_zero_elts(ai_nR, ai_nC   &
        , ai_npatch, api_nel    &
        , ao_con_C, ao_con_R    &
        , lpi_columns, lpi_occ)
    ! ...
        
        self % ol_use_mm_format = .FALSE.

        self % oi_nR = ai_nR
        self % oi_nC = ai_nC
        self % oi_nel = li_nnz

        allocate(self % opi_ia(self % oi_nR + 1))

        allocate(self % opi_ja(self % oi_nel))

        allocate(self % opr_a(self % oi_nel))

        call init_SparseMatrix(self, ai_nR, ai_nC   &
        , ai_npatch, api_nel    &
        , ao_con_C, ao_con_R    &
        , lpi_columns, lpi_occ)

        self % opr_a(:) = 0.0_wp

        deallocate(lpi_columns)
        deallocate(lpi_occ)

    end subroutine create_SparseMatrix_with_connectivity
    !---------------------------------------------------------------------------------
    subroutine init_SparseMatrix(self, ai_nR, ai_nC   &
        , ai_npatch, api_nel    &
        , ao_con_C, ao_con_R    &
        , api_columns, api_occ)
        ! _C FOR ROWS
        ! _R FOR COLUMNS
        implicit none
        type(csr_matrix) :: self
        integer :: ai_nC
        integer :: ai_nR
        integer :: ai_npatch
        integer, dimension(:) :: api_nel
        TYPE(CONNECTIVITY) :: ao_con_C
        TYPE(CONNECTIVITY) :: ao_con_R
        integer, dimension(:,:), pointer :: api_columns
        integer, dimension(:), pointer :: api_occ
        !local var
        integer :: li_nel
        integer :: li_id
        integer :: li_e
        integer :: li_b_C
        integer :: li_A_C
        !integer :: li_b_R
        !integer :: li_A_R
        !integer :: li_index
        integer :: li_i
        integer :: li_size
        integer :: li_nen_C
        integer :: li_err
        integer :: li_flag
        real(wp), dimension(:), pointer :: lpr_tmp

        ! INITIALIZING ia
        self % opi_ia(1) = 1

        do li_i = 1, self % oi_nR

            self % opi_ia(li_i + 1) = self % opi_ia(1) + SUM(api_occ(1: li_i))

        end do

        ! INITIALIZING ja
        DO li_id = 1, ai_npatch

            li_nel = api_nel (li_id)

            ! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED
            do li_e = 1, li_nel

                li_nen_C = ao_con_C % opi_nen(li_id, li_e)
                do li_b_C = 1, li_nen_C

                    li_A_C = ao_con_C % opi_LM(li_id, li_b_C, li_e)

                    if (li_A_C == 0) then
                        cycle
                    end if

                    if (api_columns(li_A_C, 0) == 0) then
                        cycle
                    end if

                    li_size = api_columns(li_A_C, 0)

                    allocate ( lpr_tmp(li_size), stat = li_err)
                    if (li_err .ne. 0) li_flag = 10

                    lpr_tmp(1: li_size) = real( api_columns(li_A_C, 1: li_size))

                    call QsortC(lpr_tmp)

                    do li_i = 1, li_size

                        self % opi_ja(self % opi_ia(li_A_C) + li_i - 1) = int ( lpr_tmp(li_i))

                    end do

                    api_columns(li_A_C, 0) = 0
                    deallocate ( lpr_tmp)

                end do

            end do

        end do

    end subroutine init_SparseMatrix
    !---------------------------------------------------------------------------------
    integer function count_non_zero_elts(ai_nR, ai_nC   &
        , ai_npatch, api_nel    &
        , ao_con_C, ao_con_R    &
        , api_columns, api_occ)
        ! _C FOR ROWS
        ! _R FOR COLUMNS
        implicit none
        !type(csr_matrix) :: self
        integer :: ai_nC
        integer :: ai_nR
        integer :: ai_npatch
        integer, dimension(:) :: api_nel
        TYPE(CONNECTIVITY) :: ao_con_C
        TYPE(CONNECTIVITY) :: ao_con_R
        integer, dimension(:,:), pointer :: api_columns
        integer, dimension(:), pointer :: api_occ
        !local var
        integer :: li_e
        integer :: li_nel
        integer :: li_id
        integer :: li_nen_C
        integer :: li_nen_R
        integer :: li_b_C
        integer :: li_A_C
        integer :: li_b_R
        integer :: li_A_R
        integer :: li_i
        integer :: li_result
        logical :: ll_done
        integer, dimension(2) :: lpi_size
        !real(wp), dimension(:), pointer :: lpr_tmp
        integer, dimension(:,:), pointer :: lpi_columns

        DO li_id = 1, ai_npatch

            li_nel = api_nel (li_id)
            
            ! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED
            do li_e = 1, li_nel
!                print *,"li_id, li_e=",li_id, li_e
                li_nen_C = ao_con_C % opi_nen(li_id, li_e)
                do li_b_C = 1, li_nen_C

                    li_A_C = ao_con_C % opi_LM(li_id, li_b_C, li_e)
                    if (li_A_C == 0) then
                        cycle
                    end if

                    li_nen_R = ao_con_R % opi_nen(li_id, li_e)
                    do li_b_R = 1, li_nen_R

                        li_A_R = ao_con_R % opi_LM(li_id, li_b_R, li_e)
                        if (li_A_R == 0) then
                            cycle
                        end if

                        ll_done = .false.
                        ! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_C, li_A_R)
                        do li_i = 1, api_columns(li_A_C, 0)

                            if (api_columns(li_A_C, li_i) /= li_A_R) then
                                cycle
                            end if

                            ll_done = .true.
                            exit

                        end do

                        if (.not.ll_done) then

                            api_occ(li_A_C) = api_occ(li_A_C) + 1

                            ! li_A_C IS THE ROW NUM, li_A_R THE COLUMN NUM
                            ! INITIALIZATION OF THE SPARSE MATRIX
                            api_columns(li_A_C, 0) = api_columns(li_A_C, 0) + 1
                            api_columns(li_A_C, api_columns(li_A_C, 0)) = li_A_R

                            ! resizing the array
                            lpi_size(1) = SIZE(api_columns, 1)
                            lpi_size(2) = SIZE(api_columns, 2)
                            if (lpi_size(2) < api_columns(li_A_C, 0)) then
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

        END DO

        ! COUNT NON ZERO ELEMENTS
        li_result = SUM(api_occ(1: ai_nR))

        count_non_zero_elts = li_result
    end function count_non_zero_elts

end module SparseMatrix_Con_Module
