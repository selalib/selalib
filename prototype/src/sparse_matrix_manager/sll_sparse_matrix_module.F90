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

module sll_sparse_matrix_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none

!> @brief type for CSR format
type sll_csr_matrix
  private
  sll_int32,  public          :: num_rows !< rows, public
  sll_int32,  public          :: num_cols !< columns
  sll_int32,  public          :: num_nz   !< non zeros
  sll_int32,  public, pointer :: opi_ia(:)
  sll_int32,  public, pointer :: opi_ja(:)
  sll_real64, public, pointer :: opr_a(:)
  sll_int32,  pointer         :: opi_i(:)
end type sll_csr_matrix

interface sll_delete
  module procedure delete_csr_matrix
end interface sll_delete
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine delete_csr_matrix(csr_mat)
  type(sll_csr_matrix),pointer :: csr_mat

  nullify(csr_mat)
    
end subroutine delete_csr_matrix

!> @brief allocates the memory space for a new CSR matrix,
!> @details
!> initializes it with the given arguments and returns a pointer to the
!> object.
!> @param[in] num_rows :  number of rows
!> @param[in] num_cols :  number of columns
!> @param[in] num_element :  number of elements
!> @param[in] local_to_global_row : local_to_global_row(\ell,i) gives the global 
!> row index of the matrix, for the element i and local degree of freedom \ell
!> @param[in] num_local_dof_row : number of local degrees of freedom for the rows
!> @param[in] local_to_global_col : local_to_global_col(\ell,i) gives the global 
!> column index of the matrix, for the element i and local degree of freedom \ell
!> @param[in] num_local_dof_col : number of local degrees of freedom for the columns
!> @returns a pointer to the newly allocated object.
function new_csr_matrix( &
  num_rows,              &
  num_cols,              &
  num_elements,          &
  local_to_global_row,   &
  num_local_dof_row,     &
  local_to_global_col,   &
  num_local_dof_col)     &
  result(mat)

  type(sll_csr_matrix), pointer :: mat
  sll_int32, intent(in) :: num_rows
  sll_int32, intent(in) :: num_cols
  sll_int32, intent(in) :: num_elements
  sll_int32, dimension(:,:), intent(in) :: local_to_global_row
  sll_int32, dimension(:,:), intent(in) :: local_to_global_col
  sll_int32, intent(in) :: num_local_dof_row
  sll_int32, intent(in) :: num_local_dof_col

  sll_int32 :: ierr

  SLL_ALLOCATE(mat, ierr)

  call initialize_csr_matrix( &
    mat,                      &
    num_rows,                 &
    num_cols,                 &
    num_elements,             &
    local_to_global_row,      &
    num_local_dof_row,        &
    local_to_global_col,      &
    num_local_dof_col)
    
end function new_csr_matrix

!> @brief initialization of CSR matrix type
!> thanks to the global index of each local dof of each element
!> @param[inout] mat : CSR matrix structure
!> @param[in] num_rows :  number of rows
!> @param[in] num_cols :  number of columns
!> @param[in] num_element :  number of elements
!> @param[in] local_to_global_row : local_to_global_row(\ell,i) gives the global 
!> row index of the matrix, for the element i and local degree of freedom \ell
!> @param[in] num_local_dof_row : number of local degrees of freedom for the rows
!> @param[in] local_to_global_col : local_to_global_col(\ell,i) gives the global 
!> column index of the matrix, for the element i and local degree of freedom \ell
!> @param[in] num_local_dof_col : number of local degrees of freedom for the columns
subroutine initialize_csr_matrix( &
  mat,                            &
  num_rows,                       &
  num_cols,                       &
  num_elements,                   &
  local_to_global_row,            &
  num_local_dof_row,              &
  local_to_global_col,            & 
  num_local_dof_col)

  type(sll_csr_matrix),      intent(inout) :: mat
  sll_int32,                 intent(in)    :: num_rows
  sll_int32,                 intent(in)    :: num_cols
  sll_int32,                 intent(in)    :: num_elements
  sll_int32, dimension(:,:), intent(in)    :: local_to_global_row
  sll_int32,                 intent(in)    :: num_local_dof_row
  sll_int32, dimension(:,:), intent(in)    :: local_to_global_col
  sll_int32,                 intent(in)    :: num_local_dof_col

  sll_int32                                :: num_nz
  sll_int32                                :: ierr
  sll_int32,  dimension(:,:), allocatable  :: lpi_columns
  sll_int32,  dimension(:),   allocatable  :: lpi_occ
  sll_int32                                :: li_COEF
  sll_int32                                :: li_err
  sll_int32                                :: li_e
  sll_int32                                :: li_b_1
  sll_int32                                :: li_A_1
  sll_int32                                :: li_b_2
  sll_int32                                :: li_A_2
  sll_int32                                :: li_i
  sll_int32                                :: li_flag
  sll_int32                                :: li_size
  sll_int32                                :: li_result
  sll_int32                                :: lpi_size(2)
  logical                                  :: ll_done

  print *,'#initialize_csr_matrix'
  li_COEF = 10

  SLL_ALLOCATE(lpi_columns(num_rows, 0:li_COEF*num_local_dof_col),ierr)
  SLL_ALLOCATE(lpi_occ(num_rows+1),ierr)

  lpi_columns(:,:) = 0
  lpi_occ(:) = 0
  
  ! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED
  do li_e = 1, num_elements

    do li_b_1 = 1, num_local_dof_row

      li_A_1 = local_to_global_row(li_b_1, li_e)

      if (li_A_1 == 0) cycle

      do li_b_2 = 1, num_local_dof_col

        li_A_2 = local_to_global_col(li_b_2, li_e)
        if (li_A_2 == 0) cycle

        ll_done = .false.
        ! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_1, li_A_2)
        do li_i = 1, lpi_columns(li_A_1, 0)
          if (lpi_columns(li_A_1, li_i) /= li_A_2) cycle
          ll_done = .true.
          exit
        end do

        if (.not.ll_done) then

          lpi_occ(li_A_1) = lpi_occ(li_A_1) + 1

          ! li_A_1 IS THE ROW NUM, li_A_2 THE COLUMN NUM
          ! INITIALIZATION OF THE SPARSE MATRIX
          lpi_columns(li_A_1, 0) = lpi_columns(li_A_1, 0) + 1
          lpi_columns(li_A_1, lpi_columns(li_A_1, 0)) = li_A_2

          ! resizing the array
          !lpi_size(1) = SIZE(lpi_columns, 1)
          !lpi_size(2) = SIZE(lpi_columns, 2)
          !if (lpi_size(2) < lpi_columns(li_A_1, 0)) then
          !  !ALLOCATE(lpi_columns(lpi_size(1), lpi_size(2)))
          !  !lpi_columns = lpi_columns
          !  !DEALLOCATE(lpi_columns)
          !  !ALLOCATE(lpi_columns(lpi_size(1), 2 * lpi_size(2)))
          !  lpi_columns(1:lpi_size(1),1:lpi_size(2)) = &
          !    lpi_columns(1:lpi_size(1), 1:lpi_size(2))
          !  DEALLOCATE(lpi_columns)
          !end if

        end if

      end do

    end do

  end do

  ! COUNT NON ZERO ELEMENTS
  num_nz = SUM(lpi_occ(1:num_rows))

  mat%num_rows = num_rows
  mat%num_cols = num_cols
  mat%num_nz   = num_nz   

  print *,'#num_rows=',num_rows
  print *,'#num_nz=',num_nz

  SLL_ALLOCATE(mat%opi_ia(num_rows + 1),ierr)
  SLL_ALLOCATE(mat%opi_ja(num_nz),ierr)
  SLL_ALLOCATE(mat%opr_a(num_nz),ierr)
  
  mat%opi_ia(1) = 1

  do li_i = 1, mat%num_rows
    mat%opi_ia(li_i + 1) = mat%opi_ia(1) + SUM(lpi_occ(1: li_i))
  end do

  do li_e = 1, num_elements

    do li_b_1 = 1, num_local_dof_row

      li_A_1 = local_to_global_row(li_b_1, li_e)

      if (li_A_1 == 0) cycle
      if (lpi_columns(li_A_1, 0) == 0) cycle

      li_size = lpi_columns(li_A_1, 0)

      if (li_err .ne. 0) li_flag = 10

      call QsortC(lpi_columns(li_A_1, 1: li_size))

      do li_i = 1, li_size
         mat%opi_ja(mat%opi_ia(li_A_1)+li_i-1) = lpi_columns(li_A_1,li_i)
      end do

      lpi_columns(li_A_1, 0) = 0

      end do

   end do

  mat%opr_a(:) = 0.0_f64
  SLL_DEALLOCATE_ARRAY(lpi_columns,ierr)
  SLL_DEALLOCATE_ARRAY(lpi_occ,ierr)

end subroutine initialize_csr_matrix

subroutine initialize_csr_matrix_with_constraint( mat, mat_a)

  type(sll_csr_matrix), intent(inout) :: mat
  type(sll_csr_matrix), intent(in) :: mat_a
  sll_int32 :: ierr

  mat%num_nz = mat_a%num_nz + 2*mat_a%num_rows       
  print*,'num_nz mat, num_nz mat_tot', mat_a%num_nz,mat%num_nz 
  mat%num_rows = mat_a%num_rows  !+  1
  print*,'num_rows mat, num_rows mat_tot',mat_a%num_rows , mat%num_rows
  mat%num_cols = mat_a%num_cols  !+  1
  print*,'num_cols mat, num_cols mat_tot',mat_a%num_cols , mat%num_cols 

  SLL_ALLOCATE(mat%opi_ia(mat%num_rows),ierr)
  SLL_ALLOCATE(mat%opi_ja(mat%num_nz),ierr)
  SLL_CLEAR_ALLOCATE(mat%opr_a(1:mat%num_nz),ierr)

end subroutine initialize_csr_matrix_with_constraint

function new_csr_matrix_with_constraint(mat_a) result(mat)

  type(sll_csr_matrix), pointer :: mat
  type(sll_csr_matrix) :: mat_a

  sll_int32 :: ierr
  SLL_ALLOCATE(mat, ierr)
  call initialize_csr_matrix_with_constraint( mat, mat_a)

end function new_csr_matrix_with_constraint

subroutine csr_add_one_constraint( &
  ia_in,                           &
  ja_in,                           &
  a_in,                            &
  num_rows_in,                     &
  num_nz_in,                       &
  constraint_vec,                  &
  ia_out,                          &
  ja_out,                          &
  a_out)

  sll_int32,  dimension(:), intent(in)  :: ia_in  
  sll_int32,  dimension(:), intent(in)  :: ja_in  
  sll_real64, dimension(:), intent(in)  :: a_in
  sll_int32,                intent(in)  :: num_rows_in
  sll_int32,                intent(in)  :: num_nz_in
  sll_real64, dimension(:), intent(in)  :: constraint_vec
  sll_int32,  dimension(:), intent(out) :: ia_out
  sll_int32,  dimension(:), intent(out) :: ja_out
  sll_real64, dimension(:), intent(out) :: a_out


  sll_int32 :: num_rows_out
  sll_int32 :: num_nz_out
  sll_int32 :: i
  sll_int32 :: s
  sll_int32 :: k
    
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
    print *,'#bad value of ia_out(num_rows_out+1)', &
      ia_out(num_rows_out+1),num_nz_out+1
    stop
  endif
  
end subroutine csr_add_one_constraint
  
subroutine sll_mult_csr_matrix_vector(mat, input, output)
    
  type(sll_csr_matrix),     intent(in)  :: mat
  sll_real64, dimension(:), intent(in)  :: input
  sll_real64, dimension(:), intent(out) :: output

  sll_int32 :: li_i
  sll_int32 :: li_k_1
  sll_int32 :: li_k_2

  do li_i = 1, mat%num_rows

    li_k_1 = mat%opi_ia(li_i)
    li_k_2 = mat%opi_ia(li_i+1)-1

    output(li_i) = & 
      dot_product(mat%opr_a(li_k_1:li_k_2),input(mat%opi_ja(li_k_1:li_k_2)))
            
  end do

end subroutine sll_mult_csr_matrix_vector

subroutine sll_add_to_csr_matrix(mat, val, ai_A, ai_Aprime)

  type(sll_csr_matrix), intent(inout) :: mat
  sll_real64, intent(in) :: val
  sll_int32, intent(in) :: ai_A
  sll_int32, intent(in) :: ai_Aprime

  sll_int32 :: li_j
  sll_int32 :: li_k

  ! THE CURRENT LINE IS self%opi_ia(ai_A)
  do li_k = mat%opi_ia(ai_A), mat%opi_ia(ai_A+1) - 1
    li_j = mat%opi_ja(li_k)
    if (li_j == ai_Aprime) then
      mat%opr_a(li_k) = mat%opr_a(li_k) + val
      exit
    end if
  end do

end subroutine sll_add_to_csr_matrix
  
subroutine sll_solve_csr_matrix(mat, apr_B, apr_U)

  type(sll_csr_matrix), intent(in) :: mat
  sll_real64, dimension(:),intent(inout) :: apr_B
  sll_real64, dimension(:),intent(out) :: apr_U

  sll_int32  :: ai_maxIter
  sll_real64 :: ar_eps

  sll_real64, dimension(:), allocatable :: lpr_Ad
  sll_real64, dimension(:), allocatable :: lpr_d

  logical    :: ll_continue
  sll_real64 :: lr_Norm2r1
  sll_real64 :: lr_Norm2r0
  sll_real64 :: lr_NormInfb
  sll_real64 :: lr_NormInfr
  sll_real64 :: lr_ps
  sll_real64 :: lr_beta
  sll_real64 :: lr_alpha
  sll_int32  :: li_iter
  sll_int32  :: li_err
  sll_int32  :: li_flag

  ar_eps = 1.d-13
  ai_maxIter = 10000

  if ( mat%num_rows /= mat%num_cols ) then
    PRINT*,'#ERROR Gradient_conj: The matrix must be square'
    stop
  end if

  if ((dabs(maxval(apr_B)) < ar_eps) .AND. (dabs(minval(apr_B)) < ar_eps)) then
    apr_U = 0.0_8
    return
  end if

  SLL_ALLOCATE(lpr_Ad(mat%num_rows),li_err)
  SLL_ALLOCATE(lpr_d(mat%num_rows),li_err)
  
  apr_U   = 0.0_8
  li_iter = 0

  lr_NormInfb = maxval(abs(apr_B))
  lr_Norm2r0  = dot_product(apr_B,apr_B)

  lpr_d = apr_B
 
  do li_iter = 1, ai_maxiter
    !--------------------------------------!
    ! calcul du ak parametre optimal local !
    !--------------------------------------!

    call sll_mult_csr_matrix_vector( mat , lpr_d , lpr_Ad )

    lr_ps = dot_product( lpr_Ad , lpr_d )
    lr_alpha = lr_Norm2r0 / lr_ps
            
    !==================================================!
    ! calcul de l'approximation Xk+1 et du residu Rk+1 !
    !==================================================!
    ! calcul des composantes residuelles
    !-----------------------------------
    apr_B = apr_B - lr_alpha * lpr_Ad
           
    !----------------------------------------!
    ! approximations ponctuelles au rang k+1 !
    !----------------------------------------!
    apr_U = apr_U + lr_alpha * lpr_d
            
    !-------------------------------------------------------!
    ! (a) extraction de la norme infinie du residu          !
    !     pour le test d'arret                              !
    ! (b) extraction de la norme euclidienne du residu rk+1 !
    !-------------------------------------------------------!
    lr_NormInfr = maxval(dabs(apr_B))
    lr_Norm2r1  = dot_product(apr_B,apr_B)

    !==================================================!
    ! calcul de la nouvelle direction de descente dk+1 !
    !==================================================!
    lr_beta    = lr_Norm2r1 / lr_Norm2r0
    lr_Norm2r0 = lr_Norm2r1
    lpr_d      = apr_B + lr_beta * lpr_d
           
    !-------------------!
    ! boucle suivante ? !
    !-------------------!
    if ( lr_NormInfr / lr_NormInfb < ar_eps ) exit
            
  end do

  if ( li_iter == ai_maxIter ) then
    print*,'Warning Gradient_conj : li_iter == ai_maxIter'
    print*,'Error after CG =',( lr_NormInfr / lr_NormInfb )
  end if
    
  deallocate(lpr_Ad)
  deallocate(lpr_d)

end subroutine sll_solve_csr_matrix


recursive subroutine QsortC(A)
  sll_int32, intent(inout), dimension(:) :: A
  sll_int32 :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  sll_int32, intent(inout), dimension(:) :: A
  sll_int32, intent(out) :: marker
  sll_int32 :: i, j
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
  type(sll_csr_matrix) :: this
  sll_real64, dimension(:) :: apr_U
  sll_real64, dimension(:) :: apr_B
  sll_int32  :: ai_maxIter
  sll_real64 :: ar_eps

  sll_real64, dimension(:), pointer :: lpr_Ad
  sll_real64, dimension(:), pointer :: lpr_r
  sll_real64, dimension(:), pointer :: lpr_d
  sll_real64, dimension(:), pointer :: lpr_Ux,one
  sll_real64, dimension(:), pointer :: Masse_tot
  sll_real64 :: lr_Norm2r1
  sll_real64 :: lr_Norm2r0
  sll_real64 :: lr_NormInfb
  sll_real64 :: lr_NormInfr
  sll_real64 :: lr_ps
  sll_real64 :: lr_beta
  sll_real64 :: lr_alpha
  logical  :: ll_continue
  sll_int32  :: li_iter
  sll_int32  :: li_err
  sll_int32  :: li_flag

  ai_maxIter = 100000
  ar_eps = 1.d-13

  if ( this%num_rows /= this%num_cols ) then
    PRINT*,'ERROR Gradient_conj: The matrix must be square'
    stop
  end if

  if ((dabs(maxval(apr_B)) < ar_eps ) .AND. (dabs(MINVAL(apr_B)) < ar_eps )) then
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

  apr_U(:)  = 0.0_8
  one(:) = 1.
  lpr_Ux(:) = apr_U(:)
  li_iter = 0
  call sll_mult_csr_matrix_vector( this , lpr_Ux , lpr_Ad )
  lpr_Ad = lpr_Ad - dot_product(Masse_tot, lpr_Ux)
  lpr_r       = apr_B - lpr_Ad
  lr_Norm2r0  = DOT_PRODUCT( lpr_r , lpr_r )
  lr_NormInfb = maxval( dabs( apr_B ) )

  lpr_d = lpr_r

  ll_continue=.true.
  do while(ll_continue)
    li_iter = li_iter + 1
    !--------------------------------------!
    ! calcul du ak parametre optimal local !
    !--------------------------------------!

    call sll_mult_csr_matrix_vector( this , lpr_d , lpr_Ad )
          
    lpr_Ad = lpr_Ad - dot_product(Masse_tot, lpr_d)
    lr_ps = dot_product( lpr_Ad , lpr_d )
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
    ll_continue=((lr_NormInfr/lr_NormInfb) >= ar_eps) .AND. (li_iter < ai_maxIter)
            
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

end module sll_sparse_matrix_module
