!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @file sll_sparse_matrix_module.F90
!> @namespace sll_sparse_matrix
!> @brief sparse matrix


module sll_sparse_matrix_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

  !> @brief type for CSR format
  type sll_csr_matrix
    sll_int32 :: num_rows !< number of rows
    sll_int32 :: num_cols !< number of columns
    sll_int32 :: num_nz !< number of non zero elements
    sll_int32, dimension(:), pointer :: opi_ia
    sll_int32, dimension(:), pointer :: opi_ja
    sll_real64, dimension(:), pointer :: opr_a
    
        !................
    !logical :: ol_use_mm_format
    sll_int32, dimension(:), pointer :: opi_i
        !................
  end type sll_csr_matrix

  interface sll_delete
     module procedure delete_csr_matrix
  end interface sll_delete
     

contains


  subroutine delete_csr_matrix(csr_mat)
    type(sll_csr_matrix),pointer :: csr_mat

    nullify(csr_mat)
   ! SLL_DEALLOCATE_ARRAY(csr_mat%opi_ia,ierr)
   ! SLL_DEALLOCATE_ARRAY(csr_mat%opi_ja,ierr)
   ! SLL_DEALLOCATE_ARRAY(csr_mat%opr_a,ierr)
   ! SLL_DEALLOCATE_ARRAY(csr_mat%opi_i,ierr)
    
    
  end subroutine delete_csr_matrix

  !> @brief allocates the memory space for a new CSR type on the heap,
  !> initializes it with the given arguments and returns a pointer to the
  !> object.
  !> param[in] num_rows :  number of rows
  !> param[in] num_cols :  number of columns
  !> param[in] num_element :  number of elements
  !> param[in] local_to_global_row : local_to_global_row(\ell,i) gives the global 
  !> row index of the matrix, for the element i and local degree of freedom \ell
  !> param[in] num_local_dof_row : number of local degrees of freedom for the rows
  !> param[in] local_to_global_col : local_to_global_col(\ell,i) gives the global 
  !> column index of the matrix, for the element i and local degree of freedom \ell
  !> param[in] num_local_dof_col : number of local degrees of freedom for the columns
  !> return a pointer to the newly allocated object.
  function new_csr_matrix( &
    num_rows, &
    num_cols, &
    num_elements, &
    local_to_global_row, &
    num_local_dof_row, &
    local_to_global_col, &
    num_local_dof_col) &
    result(mat)
    type(sll_csr_matrix), pointer :: mat
    sll_int32, intent(in) :: num_rows
    sll_int32, intent(in) :: num_cols
    sll_int32, intent(in) :: num_elements
    sll_int32, dimension(:,:), intent(in) :: local_to_global_row
    sll_int32, intent(in) :: num_local_dof_row
    sll_int32, dimension(:,:), intent(in) :: local_to_global_col
    sll_int32, intent(in) :: num_local_dof_col
    sll_int32 :: ierr
    SLL_ALLOCATE(mat, ierr)
    call initialize_csr_matrix( &
      mat, &
      num_rows, &
      num_cols, &
      num_elements, &
      local_to_global_row, &
      num_local_dof_row, &
      local_to_global_col, &
      num_local_dof_col)
      
  end function new_csr_matrix

  !> @brief initialization of CSR matrix type
  !> thanks to the global index of each local dof of each element
  !> param[inout] mat : CSR matrix structure
  !> param[in] num_rows :  number of rows
  !> param[in] num_cols :  number of columns
  !> param[in] num_element :  number of elements
  !> param[in] local_to_global_row : local_to_global_row(\ell,i) gives the global 
  !> row index of the matrix, for the element i and local degree of freedom \ell
  !> param[in] num_local_dof_row : number of local degrees of freedom for the rows
  !> param[in] local_to_global_col : local_to_global_col(\ell,i) gives the global 
  !> column index of the matrix, for the element i and local degree of freedom \ell
  !> param[in] num_local_dof_col : number of local degrees of freedom for the columns
  subroutine initialize_csr_matrix( &
    mat, &
    num_rows, &
    num_cols, &
    num_elements, &
    local_to_global_row, &
    num_local_dof_row, &
    local_to_global_col, &
    num_local_dof_col)
    type(sll_csr_matrix), intent(inout) :: mat
    sll_int32, intent(in) :: num_rows
    sll_int32, intent(in) :: num_cols
    sll_int32, intent(in) :: num_elements
    sll_int32, dimension(:,:), intent(in) :: local_to_global_row
    sll_int32, intent(in) :: num_local_dof_row
    sll_int32, dimension(:,:), intent(in) :: local_to_global_col
    sll_int32, intent(in) :: num_local_dof_col
    !local variables
    sll_int32 :: num_nz,ierr
    sll_int32, dimension(:,:), pointer :: lpi_columns
    sll_int32, dimension(:), pointer :: lpi_occ
    sll_int32 :: li_COEF
    !print *,'#num_rows=',num_rows
    !print *,'#num_nz=',num_nz

    print *,'#initialize_csr_matrix'
    li_COEF = 10
    SLL_ALLOCATE(lpi_columns(num_rows, 0:li_COEF * num_local_dof_col),ierr)
    SLL_ALLOCATE(lpi_occ(num_rows + 1),ierr)
    lpi_columns(:,:) = 0
    lpi_occ(:) = 0
    ! COUNTING NON ZERO ELEMENTS
    num_nz = sll_count_non_zero_elts( &
      num_rows, &
      num_cols, &
      num_elements, &
      local_to_global_row, &
      num_local_dof_row, &
      local_to_global_col, &
      num_local_dof_col, &
      lpi_columns, &
      lpi_occ)      
    mat%num_rows = num_rows
    mat%num_cols = num_cols
    mat%num_nz = num_nz   
    print *,'#num_rows=',num_rows
    print *,'#num_nz=',num_nz
    SLL_ALLOCATE(mat%opi_ia(num_rows + 1),ierr)
    SLL_ALLOCATE(mat%opi_ja(num_nz),ierr)
    SLL_ALLOCATE(mat%opr_a(num_nz),ierr)
    
    
    call sll_init_SparseMatrix( &
      mat, &
      num_elements, &
      local_to_global_row, &
      num_local_dof_row, &
      local_to_global_col, &
      num_local_dof_col, &
      lpi_columns, &
      lpi_occ)
    
    mat%opr_a(:) = 0.0_f64
    SLL_DEALLOCATE_ARRAY(lpi_columns,ierr)
    SLL_DEALLOCATE_ARRAY(lpi_occ,ierr)

  end subroutine initialize_csr_matrix


  subroutine sll_factorize_csr_matrix(mat)
    type(sll_csr_matrix), intent(inout) :: mat
    
    print *,'#sll_factorize_csr_matrix does nothing here'
    
  end subroutine sll_factorize_csr_matrix

  subroutine initialize_csr_matrix_with_constraint( &
    mat, &
    mat_a)
    type(sll_csr_matrix), intent(inout) :: mat
    type(sll_csr_matrix), intent(in) :: mat_a
    sll_int32 :: ierr
    !print*,' COUNTING NON ZERO ELEMENTS'
    mat%num_nz = mat_a%num_nz + 2*mat_a%num_rows       
    print*,'num_nz mat, num_nz mat_tot', mat_a%num_nz,mat%num_nz 
    mat%num_rows = mat_a%num_rows  !+  1
    print*,'num_rows mat, num_rows mat_tot',mat_a%num_rows , mat%num_rows
    mat%num_cols = mat_a%num_cols  !+  1
    print*,'num_cols mat, num_cols mat_tot',mat_a%num_cols , mat%num_cols 

    
    SLL_ALLOCATE(mat%opi_ia(mat%num_rows),ierr)
    SLL_ALLOCATE(mat%opi_ja(mat%num_nz),ierr)
    SLL_ALLOCATE(mat%opr_a(mat%num_nz),ierr)
    mat%opr_a(:) = 0.0_f64

  end subroutine initialize_csr_matrix_with_constraint

  function new_csr_matrix_with_constraint(mat_a) result(mat)
    type(sll_csr_matrix), pointer :: mat
    type(sll_csr_matrix) :: mat_a
    sll_int32 :: ierr
    SLL_ALLOCATE(mat, ierr)
    call initialize_csr_matrix_with_constraint( &
         mat, &
         mat_a)
  end function new_csr_matrix_with_constraint
  subroutine csr_add_one_constraint( &
    ia_in, &
    ja_in, &
    a_in, &
    num_rows_in, &
    num_nz_in, &
    constraint_vec, &
    ia_out, &
    ja_out, &
    a_out)
    integer, dimension(:), intent(in) :: ia_in  
    integer, dimension(:), intent(in) :: ja_in  
    real(8), dimension(:), intent(in) :: a_in
    integer, intent(in) :: num_rows_in
    integer, intent(in) :: num_nz_in
    real(8), dimension(:), intent(in) :: constraint_vec
    integer, dimension(:), intent(out) :: ia_out
    integer, dimension(:), intent(out) :: ja_out
    real(8), dimension(:), intent(out) :: a_out
    integer :: num_rows_out
    integer :: num_nz_out
    integer :: i
    integer :: s
    integer :: k
    
    num_rows_out = num_rows_in!+1
    num_nz_out = num_nz_in+2*num_rows_in
    
    if(size(ia_in)<num_rows_in) then
      print *, '#problem of size of ia_in', size(ia_in),num_rows_in!+1
      stop
    endif
    if(size(ja_in)<num_nz_in) then
      print *, '#problem of size of ja_in', size(ja_in),num_nz_in
      stop
    endif
    if(size(a_in)<num_nz_in) then
      print *, '#problem of size of a_in', size(a_in),num_nz_in
      stop
    endif
    if(size(ia_out)<num_rows_out) then
      print *, '#problem of size of ia_out', size(ia_out),num_rows_out!+1
      stop
    endif
    if(size(ja_out)<num_nz_out) then
      print *, '#problem of size of ja_out', size(ja_out),num_nz_out
      stop
    endif
    if(size(a_out)<num_nz_out) then
      print *, '#problem of size of a_out', size(a_out),num_nz_out
      stop
    endif
    if(ia_in(num_rows_in).ne.num_nz_in)then
      print *,'#bad value of ia_in(num_rows_in+1)', ia_in(num_rows_in),num_nz_in!+1
      stop
    endif
    
    s = 1
    do i=1,num_rows_in
      ia_out(i) = s
      do k = ia_in(i), ia_in(i+1)-1
        a_out(s) = a_in(k)
        ja_out(s) = ja_in(k)
        s = s+1
      enddo
      a_out(s) = constraint_vec(i)
      ja_out(s) = num_rows_out
      s = s+1
    enddo
    ia_out(num_rows_in+1) = s
    do i=1,num_rows_in
      a_out(s) = constraint_vec(i)
      ja_out(s) = i
      s = s+1      
    enddo
    ia_out(num_rows_in+2) = s
     
    if(ia_out(num_rows_out+1).ne.num_nz_out+1)then
      print *,'#bad value of ia_out(num_rows_out+1)',ia_out(num_rows_out+1),num_nz_out+1
      stop
    endif
    
  end subroutine csr_add_one_constraint
  
  subroutine sll_mult_csr_matrix_vector(mat, input, output)
    implicit none
    type(sll_csr_matrix), intent(in) :: mat
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output
    !local var
    sll_int32 :: li_i
    sll_int32 :: li_k_1
    sll_int32 :: li_k_2

    do li_i = 1, mat % num_rows

      li_k_1 = mat % opi_ia(li_i)
      li_k_2 = mat % opi_ia(li_i + 1) - 1
      output(li_i) = &
        DOT_PRODUCT(mat % opr_a(li_k_1: li_k_2), input(mat % opi_ja(li_k_1: li_k_2)))
            
    end do

  end subroutine sll_mult_csr_matrix_vector



  subroutine sll_add_to_csr_matrix(mat, val, ai_A, ai_Aprime)
    implicit none
    type(sll_csr_matrix), intent(inout) :: mat
    sll_real64, intent(in) :: val
    sll_int32, intent(in) :: ai_A
    sll_int32, intent(in) :: ai_Aprime
    !local var
    sll_int32 :: li_j
    sll_int32 :: li_k


    ! THE CURRENT LINE IS self%opi_ia(ai_A)
    do li_k = mat % opi_ia(ai_A), mat % opi_ia(ai_A + 1) - 1
      li_j = mat % opi_ja(li_k)
      if (li_j == ai_Aprime) then
        mat % opr_a(li_k) = mat % opr_a(li_k) + val
        exit
      end if
    end do

  end subroutine sll_add_to_csr_matrix
  


  subroutine sll_solve_csr_matrix(mat, apr_B, apr_U)
    implicit none
    type(sll_csr_matrix), intent(in) :: mat
    sll_real64, dimension(:),intent(in) :: apr_B
    sll_real64, dimension(:),intent(out) :: apr_U
    !local var
    !sll_int32  :: sys
    !sll_real64, dimension(umfpack_info) :: info
    !sys = 0
    !call umf4sol(sys,apr_U,apr_B,mat%umf_numeric,mat%umf_control,info)
    !use SparseMatrix_Module
    !implicit none
    !type(csr_matrix) :: this
    !real(8), dimension(:) :: apr_U
    !real(8), dimension(:) :: apr_B
    integer  :: ai_maxIter
    real(8) :: ar_eps
    !local var
    real(8), dimension(:), pointer :: lpr_Ad
    real(8), dimension(:), pointer :: lpr_r
    real(8), dimension(:), pointer :: lpr_d
    real(8), dimension(:), pointer :: lpr_Ux
    real(8) :: lr_Norm2r1
    real(8) :: lr_Norm2r0
    real(8) :: lr_NormInfb
    real(8) :: lr_NormInfr
    real(8) :: lr_ps
    real(8) :: lr_beta
    real(8) :: lr_alpha
    logical  :: ll_continue
    integer  :: li_iter
    integer  :: li_err
    integer  :: li_flag
	
    ar_eps = 1.d-13
    ai_maxIter = 100000
	
		
    if ( mat%num_rows /= mat%num_cols ) then
            PRINT*,'#ERROR Gradient_conj: The matrix must be square'
            stop
    end if

    if ( ( dabs ( MAXVAL ( apr_B ) ) < ar_eps ) .AND. ( dabs ( MINVAL ( apr_B ) ) < ar_eps ) ) then
            apr_U = 0.0_8
            return
    end if

    allocate(lpr_Ad(mat%num_rows),stat=li_err)
    if (li_err.ne.0) li_flag=10
    allocate(lpr_r(mat%num_rows),stat=li_err)
    if (li_err.ne.0) li_flag=20
    allocate(lpr_d(mat%num_rows),stat=li_err)
    if (li_err.ne.0) li_flag=30
    allocate(lpr_Ux(mat%num_rows),stat=li_err)
    if (li_err.ne.0) li_flag=40
    !================!
    ! initialisation !
    !================!
    
    apr_U(:)  = 0.0_8
    lpr_Ux(:) = apr_U(:)
    li_iter = 0
    call sll_mult_csr_matrix_vector( mat , lpr_Ux , lpr_Ad )
    !-------------------!
    ! calcul des normes !
    !-------------------!
    lpr_r       = apr_B - lpr_Ad
    lr_Norm2r0  = DOT_PRODUCT( lpr_r , lpr_r )
    lr_NormInfb = maxval( dabs( apr_B ) )

    lpr_d = lpr_r
    !================!

 
    ll_continue=.true.
    do while(ll_continue)
            li_iter = li_iter + 1
            !--------------------------------------!
            ! calcul du ak parametre optimal local !
            !--------------------------------------!

            call sll_mult_csr_matrix_vector( mat , lpr_d , lpr_Ad )
            lr_ps = DOT_PRODUCT( lpr_Ad , lpr_d )
            lr_alpha = lr_Norm2r0 / lr_ps
            
            !==================================================!
            ! calcul de l'approximation Xk+1 et du residu Rk+1 !
            !==================================================!
            ! calcul des composantes residuelles
            !-----------------------------------
            lpr_r = lpr_r - lr_alpha * lpr_Ad
            
            !----------------------------------------!
            ! approximations ponctuelles au rang k+1 !
            !----------------------------------------!
            lpr_Ux = lpr_Ux + lr_alpha * lpr_d
            
            !-------------------------------------------------------!
            ! (a) extraction de la norme infinie du residu          !
            !     pour le test d'arret                              !
            ! (b) extraction de la norme euclidienne du residu rk+1 !
            !-------------------------------------------------------!
            lr_NormInfr = maxval(dabs( lpr_r ))
            lr_Norm2r1 = DOT_PRODUCT( lpr_r , lpr_r )
            !==================================================!
            ! calcul de la nouvelle direction de descente dk+1 !
            !==================================================!
            lr_beta = lr_Norm2r1 / lr_Norm2r0
            lr_Norm2r0 = lr_Norm2r1
            lpr_d = lpr_r + lr_beta * lpr_d
            
            !-------------------!
            ! boucle suivante ? !
            !-------------------!
            ll_continue=( ( lr_NormInfr / lr_NormInfb ) >= ar_eps ) .AND. ( li_iter < ai_maxIter )
            !print*, 'norme infr = ',lr_NormInfr, 'norme infb=',  lr_NormInfb
            
    end do
    apr_U = lpr_Ux
    if ( li_iter == ai_maxIter ) then
            print*,'Warning Gradient_conj : li_iter == ai_maxIter'
            print*,'Error after CG =',( lr_NormInfr / lr_NormInfb )
    end if

    
    deallocate(lpr_Ad)
    deallocate(lpr_d)
    deallocate(lpr_r)
    deallocate(lpr_Ux)

  end subroutine sll_solve_csr_matrix



    integer function sll_count_non_zero_elts( &
      ai_nR,&
      ai_nC,&
      ai_nel,&
      api_LM_1,&
      ai_nen_1,&
      api_LM_2,&
      ai_nen_2,&
      api_columns,&
      api_occ)
        ! _1 FOR ROWS
        ! _2 FOR COLUMNS
        implicit none
        integer :: ai_nR, ai_nC
        integer, dimension(:,:), intent(in) :: api_LM_1, api_LM_2
        integer :: ai_nel, ai_nen_1, ai_nen_2
        integer, dimension(:,:), pointer :: api_columns
        integer, dimension(:), pointer :: api_occ
        !local var
        integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_i
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

        sll_count_non_zero_elts = li_result
    end function sll_count_non_zero_elts

    subroutine sll_init_SparseMatrix(&
         self,&
         ai_nel,&
         api_LM_1,&
         ai_nen_1,&
         api_LM_2, &
         ai_nen_2,&
         api_columns,&
         api_occ)
        ! _1 FOR ROWS
        ! _2 FOR COLUMNS
        implicit none
        type(sll_csr_matrix) :: self
        integer, dimension(:,:), intent(in) :: api_LM_1
        integer, dimension(:,:), intent(in) :: api_LM_2
        integer :: ai_nel, ai_nen_1, ai_nen_2
        integer, dimension(:,:), pointer :: api_columns
        integer, dimension(:), pointer :: api_occ
        !local var
        integer :: li_e, li_b_1, li_A_1,li_i, li_size
        integer :: li_err, li_flag
        real(f64), dimension(:), pointer :: lpr_tmp

        ! INITIALIZING ia
        self % opi_ia(1) = 1

        do li_i = 1, self % num_rows

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

    end subroutine sll_init_SparseMatrix





recursive subroutine QsortC(A)
  real(f64), intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real(f64), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real(f64) :: temp
  real(f64) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition



subroutine sll_solve_csr_matrix_perper ( this, apr_B,apr_U,Masse_tot )
    implicit none
    type(sll_csr_matrix) :: this
    real(8), dimension(:) :: apr_U
    real(8), dimension(:) :: apr_B
    integer  :: ai_maxIter
    real(8) :: ar_eps
    !local var
    real(8), dimension(:), pointer :: lpr_Ad
    real(8), dimension(:), pointer :: lpr_r
    real(8), dimension(:), pointer :: lpr_d
    real(8), dimension(:), pointer :: lpr_Ux,one
    real(8), dimension(:), pointer :: Masse_tot
    real(8) :: lr_Norm2r1
    real(8) :: lr_Norm2r0
    real(8) :: lr_NormInfb
    real(8) :: lr_NormInfr
    real(8) :: lr_ps
    real(8) :: lr_beta
    real(8) :: lr_alpha
    logical  :: ll_continue
    integer  :: li_iter
    integer  :: li_err
    integer  :: li_flag

    ai_maxIter = 100000
    ar_eps = 1.d-13
		
    if ( this%num_rows /= this%num_cols ) then
            PRINT*,'ERROR Gradient_conj: The matrix must be square'
            stop
    end if

    if ( ( dabs ( MAXVAL ( apr_B ) ) < ar_eps ) .AND. ( dabs ( MINVAL ( apr_B ) ) < ar_eps ) ) then
            apr_U = 0.0_8
            return
    end if

    allocate(lpr_Ad(this%num_rows),stat=li_err)
    if (li_err.ne.0) li_flag=10
    allocate(lpr_r(this%num_rows),stat=li_err)
    if (li_err.ne.0) li_flag=20
    allocate(lpr_d(this%num_rows),stat=li_err)
    if (li_err.ne.0) li_flag=30
    allocate(lpr_Ux(this%num_rows),stat=li_err)
    if (li_err.ne.0) li_flag=40
    allocate(one(this%num_rows),stat=li_err)
    if (li_err.ne.0) li_flag=50
    !================!
    ! initialisation !
    !================!
    !apr_U(3:this%num_rows)  = 0.0_8
    !apr_U(1)  = 1.0_8
    !apr_U(2)  = -1.0_8
    apr_U(:)  = 0.0_8
    one(:) = 1.
    lpr_Ux(:) = apr_U(:)
    li_iter = 0
    call sll_mult_csr_matrix_vector( this , lpr_Ux , lpr_Ad )
    lpr_Ad = lpr_Ad - dot_product(Masse_tot, lpr_Ux)
    !print*, 'calcul des normes',dot_product(lpr_Ad, one),maxval(lpr_Ad)
    !-----------------
    !-------------------!
    ! calcul des normes !
    !-------------------!
    lpr_r       = apr_B - lpr_Ad
    lr_Norm2r0  = DOT_PRODUCT( lpr_r , lpr_r )
    lr_NormInfb = maxval( dabs( apr_B ) )

    lpr_d = lpr_r
    !================!

!    print *,'%%%%'
    ll_continue=.true.
    do while(ll_continue)
            li_iter = li_iter + 1
            !--------------------------------------!
            ! calcul du ak parametre optimal local !
            !--------------------------------------!

            call sll_mult_csr_matrix_vector( this , lpr_d , lpr_Ad )
            
            lpr_Ad = lpr_Ad - dot_product(Masse_tot, lpr_d)
            lr_ps = DOT_PRODUCT( lpr_Ad , lpr_d )
            lr_alpha = lr_Norm2r0 / lr_ps
            
            !==================================================!
            ! calcul de l'approximation Xk+1 et du residu Rk+1 !
            !==================================================!
            ! calcul des composantes residuelles
            !-----------------------------------
            lpr_r = lpr_r - lr_alpha * lpr_Ad
            
            !----------------------------------------!
            ! approximations ponctuelles au rang k+1 !
            !----------------------------------------!
            lpr_Ux = lpr_Ux + lr_alpha * lpr_d
            
            !-------------------------------------------------------!
            ! (a) extraction de la norme infinie du residu          !
            !     pour le test d'arret                              !
            ! (b) extraction de la norme euclidienne du residu rk+1 !
            !-------------------------------------------------------!
            lr_NormInfr = maxval(dabs( lpr_r ))
            lr_Norm2r1 = DOT_PRODUCT( lpr_r , lpr_r )
            
            !==================================================!
            ! calcul de la nouvelle direction de descente dk+1 !
            !==================================================!
            lr_beta = lr_Norm2r1 / lr_Norm2r0
            lr_Norm2r0 = lr_Norm2r1
            lpr_d = lpr_r + lr_beta * lpr_d
            !lpr_d(1) =  apr_B(1)
            !-------------------!
            ! boucle suivante ? !
            !-------------------!
            ll_continue=( ( lr_NormInfr / lr_NormInfb ) >= ar_eps ) .AND. ( li_iter < ai_maxIter )
            !print*, 'norme infr = ',lr_NormInfr, 'norme infb=',  lr_NormInfb
            
    end do
    apr_U = lpr_Ux
    
    if ( li_iter == ai_maxIter ) then
            print*,'Warning Gradient_conj : li_iter == ai_maxIter'
            print*,'Error after CG =',( lr_NormInfr / lr_NormInfb )
    end if
    
    deallocate(lpr_Ad)
    deallocate(lpr_d)
    deallocate(lpr_r)
    deallocate(lpr_Ux)
    deallocate(one)
  end subroutine sll_solve_csr_matrix_perper


!    call create_CSR( &
!        es%csr_mat, &
!        solution_size, &
!        solution_size, &
!        num_cells_eta1*num_cells_eta2, &
!        es%local_to_global_spline_indices, &
!        es%total_num_splines_loc, &
!        es%local_to_global_spline_indices, &
!        es%total_num_splines_loc )




end module sll_sparse_matrix_module
