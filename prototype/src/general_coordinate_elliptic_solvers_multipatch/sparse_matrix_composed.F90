!     
! File:   sparse_matrix_composed.F90
! Author: abdelkaderratnani
!
! Created on January 29, 2012, 9:44 AM
!


module SparseMatrix_Composed_Module
    use used_precision
    use sparse_matrix_def
    implicit none

contains
    !---------------------------------------------------------------------------------
    !> \todo a tester
    !> create the matrix :
    !>  A , B
    !>  C , D
    !> INPUTS : A, B, C, D
    !> OUTPUT : this
    !> We must have :
    !> ao_A%oi_nC = ao_C%oi_nC
    !> ao_B%oi_nC = ao_D%oi_nC
    !> ao_A%oi_nR = ao_B%oi_nR
    !> ao_C%oi_nR = ao_D%oi_nR
    subroutine create_SparseMatrixProfil_From4(this, ao_A, ao_B, ao_C, ao_D)
        implicit none
        !> param[inout] this : CSR MATRIX STRUCTURE
        type(csr_matrix) :: this
        !> param[in] ao_A, ao_B, ao_C, ao_D : INPUT MATRICES
        type(csr_matrix) :: ao_A
        type(csr_matrix) :: ao_B
        type(csr_matrix) :: ao_C
        type(csr_matrix) :: ao_D
        !local var
        integer :: li_err, li_flag
        integer :: li_i
        integer :: li_j
        integer :: li_index
        integer :: li_k
        integer :: li_iloc
        integer :: li_jloc

        if ((ao_A % oi_nC /= ao_C % oi_nC) .OR. &
            (ao_B % oi_nC /= ao_D % oi_nC) .OR. &
            (ao_A % oi_nR /= ao_B % oi_nR) .OR. &
            (ao_C % oi_nR /= ao_D % oi_nR) &
            ) then

            print*, 'Error create_Sparse: ao_A, ao_B, ao_C, ao_D are not compatibles'
            STOP

        end if

        this % ol_use_mm_format = .FALSE.

        this % oi_nR = ao_A % oi_nR + ao_C % oi_nR
        this % oi_nC = ao_A % oi_nC + ao_B % oi_nC
        this % oi_nel = ao_A % oi_nel &
        +ao_B % oi_nel &
        +ao_C % oi_nel &
        +ao_D % oi_nel

        allocate(this % opi_ia(this % oi_nR + 1), stat = li_err)
        if (li_err .ne. 0) li_flag = 10

        allocate(this % opi_ja(this % oi_nel), stat = li_err)
        if (li_err .ne. 0) li_flag = 20

        allocate(this % opr_a(this % oi_nel), stat = li_err)
        if (li_err .ne. 0) li_flag = 30

        this % opr_a(:) = 0.0_wp

        !remplissage de la matrice, a partir des matrices A,B,C et D
        !print*,'ao_A%oi_nR=',ao_A%oi_nR
        li_index = 1
        do li_i = 1, this % oi_nR

            this % opi_ia(li_i) = li_index
            ! remplissage des termes de A, puis B
            if (li_i <= ao_A % oi_nR) then

                !*********************************************************************
                ! copie des termes de A
                !*********************************************************************
                ! initialisation de la colonne
                li_k = ao_A % opi_ia(li_i)
                !print*,'-------- begin new row'
                do while (li_k < ao_A % opi_ia(li_i + 1))
                    !print*,'li_indexA=',li_index
                    !print*,'li_k=',li_k
                    ! on recupere la colonne de A
                    li_j = ao_A % opi_ja(li_k)

                    ! copie la valeur de Aij
                    this % opi_ja(li_index) = li_j
                    !                    this%opr_a  ( li_index ) = ao_A%opr_a ( li_k )

                    ! incrementation de l'indice
                    li_index = li_index + 1
                    li_k = li_k + 1

                end do
                !*********************************************************************
                !print*,'-------- done'
                !*********************************************************************
                ! copie des termes de B
                !*********************************************************************
                ! initialisation de la colonne
                li_k = ao_B % opi_ia(li_i)

                do while (li_k < ao_B % opi_ia(li_i + 1))
                    !print*,'li_indexB=',li_index
                    ! on recupere la colonne de A
                    li_j = ao_B % opi_ja(li_k)

                    li_jloc = li_j + ao_A % oi_nC

                    ! copie la valeur de Aij
                    this % opi_ja(li_index) = li_jloc
                    !                    this%opr_a  ( li_index ) = ao_B%opr_a ( li_k )

                    ! incrementation de l'indice
                    li_index = li_index + 1
                    li_k = li_k + 1

                end do
                
!                this % opi_ia(li_i) = li_jloc
                !*********************************************************************
                !print*,'done'
            end if

            ! remplissage des termes de C, puis D
            if (li_i >= ao_A % oi_nR + 1) then

                li_iloc = li_i - ao_A % oi_nR

                !*********************************************************************
                ! copie des termes de C
                !*********************************************************************
                ! initialisation de la colonne
                li_k = ao_C % opi_ia(li_iloc)
                !print*,'-------- begin new row'
                do while (li_k < ao_C % opi_ia(li_iloc + 1))
                    !print*,'li_indexC=',li_index
                    ! on recupere la colonne de A
                    li_j = ao_C % opi_ja(li_k)

                    ! copie la valeur de Aij
                    this % opi_ja(li_index) = li_j
                    !                    this%opr_a  ( li_index ) = ao_C%opr_a ( li_k )

                    ! incrementation de l'indice
                    li_index = li_index + 1
                    li_k = li_k + 1

                end do
                !*********************************************************************
                !print*,'done'
                !*********************************************************************
                ! copie des termes de D
                !*********************************************************************
                ! initialisation de la colonne
                li_k = ao_D % opi_ia(li_iloc)

                do while (li_k < ao_D % opi_ia(li_iloc + 1))
                    !print*,'li_indexD=',li_index
                    ! on recupere la colonne de A
                    li_j = ao_D % opi_ja(li_k)

                    li_jloc = li_j + ao_A % oi_nC

                    ! copie la valeur de Aij
                    this % opi_ja(li_index) = li_jloc
                    !                    this%opr_a  ( li_index ) = ao_D%opr_a ( li_k )

                    ! incrementation de l'indice
                    li_index = li_index + 1
                    li_k = li_k + 1

                end do
                
!                this % opi_ia(li_i) = li_jloc
                !*********************************************************************
                !print*,'-------- done'
            end if

        end do

        ! mise a jour du dernier terme
        this % opi_ia(this % oi_nR + 1) = li_index
!        print*,'this%opi_ia=',this%opi_ia(:)

    end subroutine create_SparseMatrixProfil_From4
    !---------------------------------------------------------------------------------
    !> \todo a tester et ameliorer la performance
    !> create the matrix :
    !>  A , B
    !>  C , D
    !> INPUTS : A, B, C, D
    !> api_id (i,j) = 1 if we want to copy the Aij matrix, 0 otherwise
    !> OUTPUT : this
    !> We must have :
    !> ao_A%oi_nC = ao_C%oi_nC
    !> ao_B%oi_nC = ao_D%oi_nC
    !> ao_A%oi_nR = ao_B%oi_nR
    !> ao_C%oi_nR = ao_D%oi_nR
    subroutine copy_SparseMatrix_From4(this, ao_A, ao_B, ao_C, ao_D, api_id)
        implicit none
        !> param[inout] this : CSR MATRIX STRUCTURE
        type(csr_matrix) :: this
        !> param[in] ao_A, ao_B, ao_C, ao_D : INPUT MATRICES
        type(csr_matrix) :: ao_A
        type(csr_matrix) :: ao_B
        type(csr_matrix) :: ao_C
        type(csr_matrix) :: ao_D
        integer, dimension(:,:) :: api_id
        !local var
        integer :: li_i
        integer :: li_j
        integer :: li_index
        integer :: li_k
        integer :: li_iloc
        integer :: li_jloc


        !remplissage de la matrice, a partir des matrices A,B,C et D
        !print*,'ao_A%oi_nR=',ao_A%oi_nR
        li_index = 1
        do li_i = 1, this % oi_nR

            ! remplissage des termes de A, puis B
            if (li_i <= ao_A % oi_nR) then

                !*********************************************************************
                ! copie des termes de A
                !*********************************************************************
                ! initialisation de la colonne
                li_k = ao_A % opi_ia(li_i)
                !print*,'-------- begin new row'
                do while (li_k < ao_A % opi_ia(li_i + 1))
                    !print*,'li_indexA=',li_index
                    !print*,'li_k=',li_k
                    ! on recupere la colonne de A
                    li_j = ao_A % opi_ja(li_k)

                    ! copie la valeur de Aij
                    if (api_id(1,1)==1) then
                        this%opr_a  ( li_index ) = ao_A%opr_a ( li_k )
                    end if

                    ! incrementation de l'indice
                    li_index = li_index + 1
                    li_k = li_k + 1

                end do
                !*********************************************************************
                !print*,'-------- done'
                !*********************************************************************
                ! copie des termes de B
                !*********************************************************************
                ! initialisation de la colonne
                li_k = ao_B % opi_ia(li_i)

                do while (li_k < ao_B % opi_ia(li_i + 1))
                    !print*,'li_indexB=',li_index
                    ! on recupere la colonne de A
                    li_j = ao_B % opi_ja(li_k)

                    li_jloc = li_j + ao_A % oi_nC

                    ! copie la valeur de Aij
                    if (api_id(1,2)==1) then
                        this%opr_a  ( li_index ) = ao_B%opr_a ( li_k )
                    end if

                    ! incrementation de l'indice
                    li_index = li_index + 1
                    li_k = li_k + 1

                end do
                !*********************************************************************
                !print*,'done'
            end if

            ! remplissage des termes de C, puis D
            if (li_i >= ao_A % oi_nR + 1) then

                li_iloc = li_i - ao_A % oi_nR

                !*********************************************************************
                ! copie des termes de C
                !*********************************************************************
                ! initialisation de la colonne
                li_k = ao_C % opi_ia(li_iloc)
                !print*,'-------- begin new row'
                do while (li_k < ao_C % opi_ia(li_iloc + 1))
                    !print*,'li_indexC=',li_index
                    ! on recupere la colonne de A
                    li_j = ao_C % opi_ja(li_k)

                    ! copie la valeur de Aij
                    if (api_id(2,1)==1) then
                        this%opr_a  ( li_index ) = ao_C%opr_a ( li_k )
                    end if

                    ! incrementation de l'indice
                    li_index = li_index + 1
                    li_k = li_k + 1

                end do
                !*********************************************************************
                !print*,'done'
                !*********************************************************************
                ! copie des termes de D
                !*********************************************************************
                ! initialisation de la colonne
                li_k = ao_D % opi_ia(li_iloc)

                do while (li_k < ao_D % opi_ia(li_iloc + 1))
                    !print*,'li_indexD=',li_index
                    ! on recupere la colonne de A
                    li_j = ao_D % opi_ja(li_k)

                    li_jloc = li_j + ao_A % oi_nC

                    ! copie la valeur de Aij
                    if (api_id(2,2)==1) then
                        this%opr_a  ( li_index ) = ao_D%opr_a ( li_k )
                    end if

                    ! incrementation de l'indice
                    li_index = li_index + 1
                    li_k = li_k + 1

                end do
                !*********************************************************************
                !print*,'-------- done'
            end if

        end do

    end subroutine copy_SparseMatrix_From4
end module SparseMatrix_Composed_Module


