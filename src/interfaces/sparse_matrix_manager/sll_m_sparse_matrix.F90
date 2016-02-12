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

!> @ingroup sparse_matrix_module
!> @brief Sparse matrix linear solver utilities
!> @details This part of selalib is derived from SPM library
!> developed by Ahmed Ratnani (http://ratnani.org/spm_doc/html/)
!>
!> Michel Mehrenberger did the implementation of UMFPACK option 
!>
module sll_m_sparse_matrix
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use sll_m_qsort_partition, only : sll_s_qsortc

#ifdef UMFPACK
  use sll_m_umfpack
  use iso_fortran_env, only: output_unit
#endif /* UMFPACK */

implicit none

public :: &
  sll_f_new_csr_matrix, &
  sll_f_new_csr_matrix_with_constraint, &
  sll_s_csr_add_one_constraint, &
  sll_s_csr_todense, &
  sll_s_delete_csr_matrix, &
  !sll_s_initialize_csr_matrix, &
  sll_s_initialize_csr_matrix_with_constraint, &
  sll_s_add_to_csr_matrix, &
  sll_t_csr_matrix, &
  sll_s_factorize_csr_matrix, &
  sll_s_mult_csr_matrix_vector, &
  sll_s_solve_csr_matrix, &
  sll_s_solve_csr_matrix_perper, &
  sll_p_umfpack,                 &
  sll_o_delete

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief type for CSR format
type :: sll_t_csr_matrix
  sll_int32           :: num_rows !< rows, public
  sll_int32           :: num_cols !< columns
  sll_int32           :: num_nz   !< non zeros
  sll_int32,  pointer :: row_ptr(:)
  sll_int32,  pointer :: col_ind(:)
  sll_real64, pointer :: val(:)
  sll_int32           :: solver = 0

#ifdef UMFPACK
  sll_int32,  dimension(:), pointer :: Ai
  sll_int32,  dimension(:), pointer :: Ap

  integer(umf_void)                 :: umf_symbolic
  integer(umf_void)                 :: umf_numeric
  sll_real64, dimension(:), pointer :: umf_control
#endif /* UMFPACK */

end type sll_t_csr_matrix

sll_int32, parameter :: sll_p_umfpack = 1


interface sll_o_delete
  module procedure sll_s_delete_csr_matrix
end interface sll_o_delete
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sll_s_delete_csr_matrix(csr_mat)
type(sll_t_csr_matrix),pointer :: csr_mat
nullify(csr_mat)
end subroutine sll_s_delete_csr_matrix

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
function sll_f_new_csr_matrix( &
        num_rows,              &
        num_cols,              &
        num_elements,          &
        local_to_global_row,   &
        num_local_dof_row,     &
        local_to_global_col,   &
        num_local_dof_col,     &
        solver                 ) result(mat)

type(sll_t_csr_matrix), pointer          :: mat
sll_int32,                    intent(in) :: num_rows
sll_int32,                    intent(in) :: num_cols
sll_int32,                    intent(in) :: num_elements
sll_int32, dimension(:,:),    intent(in) :: local_to_global_row
sll_int32, dimension(:,:),    intent(in) :: local_to_global_col
sll_int32,                    intent(in) :: num_local_dof_row
sll_int32,                    intent(in) :: num_local_dof_col
sll_int32,                    optional   :: solver

sll_int32 :: ierr

SLL_ALLOCATE(mat, ierr)

if (present(solver)) then
  mat%solver = solver
else
  mat%solver = 0
endif

call sll_s_initialize_csr_matrix( &
  mat,                            &
  num_rows,                       &
  num_cols,                       &
  num_elements,                   &
  local_to_global_row,            &
  num_local_dof_row,              &
  local_to_global_col,            &
  num_local_dof_col)

    
end function sll_f_new_csr_matrix


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

subroutine sll_s_initialize_csr_matrix( mat,                 &
                                        num_rows,            &
                                        num_cols,            &
                                        num_elements,        &
                                        local_to_global_row, &
                                        num_local_dof_row,   &
                                        local_to_global_col, & 
                                        num_local_dof_col)

type(sll_t_csr_matrix),    intent(inout) :: mat
sll_int32,                 intent(in)    :: num_rows
sll_int32,                 intent(in)    :: num_cols
sll_int32,                 intent(in)    :: num_elements
sll_int32, dimension(:,:), intent(in)    :: local_to_global_row
sll_int32,                 intent(in)    :: num_local_dof_row
sll_int32, dimension(:,:), intent(in)    :: local_to_global_col
sll_int32,                 intent(in)    :: num_local_dof_col

sll_int32                                :: num_nz
sll_int32                                :: ierr
sll_int32,  dimension(:,:), allocatable  :: lpi_col
sll_int32,  dimension(:),   allocatable  :: lpi_occ
sll_int32                                :: COEF
sll_int32                                :: elt
sll_int32                                :: ii
sll_int32                                :: jj
sll_int32                                :: row
sll_int32                                :: col
sll_int32                                :: i
sll_int32                                :: sz
logical                                  :: ll_done

#ifdef DEBUG
print *,'#sll_s_initialize_csr_matrix'
#endif

COEF = 6

SLL_ALLOCATE(lpi_col(num_rows, 0:COEF*num_local_dof_col),ierr)
SLL_ALLOCATE(lpi_occ(num_rows+1),ierr)

lpi_col(:,:) = 0
lpi_occ(:) = 0

do elt = 1, num_elements  !Loop over cells
  do ii = 1, num_local_dof_row
    row = local_to_global_row(ii, elt) !Row number in matrix
    if (row /= 0) then
      do jj = 1, num_local_dof_col
        col = local_to_global_col(jj, elt) !Column number in matrix
        if (col /= 0) then
          ll_done = .false.
          ! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (row, col)
          do i = 1, lpi_col(row, 0)
            if (lpi_col(row, i) == col) then
              ll_done = .true.
              exit
            end if
          end do
          if (.not.ll_done) then
            lpi_occ(row)                  = lpi_occ(row) + 1
            lpi_col(row, 0)               = lpi_col(row, 0) + 1
            lpi_col(row, lpi_col(row, 0)) = col
          end if
        end if
      end do
    end if
  end do
end do

! COUNT NON ZERO ELEMENTS
num_nz = sum(lpi_occ(1:num_rows))

mat%num_rows = num_rows
mat%num_cols = num_cols
mat%num_nz   = num_nz   

#ifdef DEBUG
print *,'#num_rows=',num_rows
print *,'#num_cols=',num_cols
print *,'#num_nz=',num_nz
#endif

SLL_ALLOCATE(mat%row_ptr(num_rows + 1),ierr)
SLL_ALLOCATE(mat%col_ind(num_nz),ierr)
SLL_ALLOCATE(mat%val(num_nz),ierr)

mat%row_ptr(1) = 1

do i = 1, mat%num_rows
  mat%row_ptr(i+1) = mat%row_ptr(1) + sum(lpi_occ(1:i))
end do

do elt = 1, num_elements
  do ii = 1, num_local_dof_row
    row = local_to_global_row(ii, elt)
    if (row /= 0) then
      if (lpi_col(row,0) /= 0) then
        sz = lpi_col(row, 0)
        call sll_s_qsortc(lpi_col(row,1:sz))
        do i = 1, sz
          mat%col_ind(mat%row_ptr(row)+i-1) = lpi_col(row,i)
        end do
        lpi_col(row, 0) = 0
      end if
    end if
  end do
end do

#ifdef UMFPACK
if (mat%solver == sll_p_umfpack) then
  SLL_ALLOCATE(mat%umf_control(umfpack_control),ierr)
  SLL_ALLOCATE(mat%Ai(num_nz),ierr)
  SLL_ALLOCATE(mat%Ap(num_rows+1),ierr)
  
  call umf4def(mat%umf_control)  ! get the default configuration
end if
#endif /* UMFPACK */

mat%val(:) = 0.0_f64
SLL_DEALLOCATE_ARRAY(lpi_col,ierr)
SLL_DEALLOCATE_ARRAY(lpi_occ,ierr)

end subroutine sll_s_initialize_csr_matrix

  
subroutine sll_s_initialize_csr_matrix_with_constraint( mat, mat_a)

type(sll_t_csr_matrix), intent(inout) :: mat
type(sll_t_csr_matrix), intent(in) :: mat_a
sll_int32 :: ierr

mat%num_nz = mat_a%num_nz + 2*mat_a%num_rows       
print*,'num_nz mat, num_nz mat_tot', mat_a%num_nz,mat%num_nz 
mat%num_rows = mat_a%num_rows  +  1
print*,'num_rows mat, num_rows mat_tot',mat_a%num_rows , mat%num_rows
mat%num_cols = mat_a%num_cols  +  1
print*,'num_cols mat, num_cols mat_tot',mat_a%num_cols , mat%num_cols 

SLL_ALLOCATE(mat%row_ptr(mat%num_rows+1),ierr)
SLL_ALLOCATE(mat%col_ind(mat%num_nz),ierr)
SLL_CLEAR_ALLOCATE(mat%val(1:mat%num_nz),ierr)

#ifdef UMFPACK
if (mat%solver == sll_p_umfpack) then
  SLL_ALLOCATE(mat%umf_control(umfpack_control),ierr)
  SLL_ALLOCATE(mat%Ai(mat%num_nz),ierr)
  SLL_ALLOCATE(mat%Ap(mat%num_rows+1),ierr)
  call umf4def(mat%umf_control)  
end if
#endif /* UMFPACK */

end subroutine sll_s_initialize_csr_matrix_with_constraint

function sll_f_new_csr_matrix_with_constraint(mat_a) result(mat)

type(sll_t_csr_matrix), pointer :: mat
type(sll_t_csr_matrix)          :: mat_a

sll_int32 :: ierr
SLL_ALLOCATE(mat, ierr)
call sll_s_initialize_csr_matrix_with_constraint( mat, mat_a)

end function sll_f_new_csr_matrix_with_constraint

subroutine sll_s_csr_add_one_constraint( ia_in,          &
                                         ja_in,          &
                                         a_in,           &
                                         num_rows_in,    &
                                         num_nz_in,      &
                                         constraint_vec, &
                                         ia_out,         &
                                         ja_out,         &
                                         a_out)

integer, dimension(:), intent(in)  :: ia_in  
integer, dimension(:), intent(in)  :: ja_in  
real(8), dimension(:), intent(in)  :: a_in
integer,               intent(in)  :: num_rows_in
integer,               intent(in)  :: num_nz_in
real(8), dimension(:), intent(in)  :: constraint_vec
integer, dimension(:), intent(out) :: ia_out
integer, dimension(:), intent(out) :: ja_out
real(8), dimension(:), intent(out) :: a_out

integer :: num_rows_out
integer :: num_nz_out
integer :: i
integer :: s
integer :: k

num_rows_out = num_rows_in+1
num_nz_out   = num_nz_in+2*num_rows_in

SLL_ASSERT(size(ia_in)          >= num_rows_in+1)
SLL_ASSERT(size(ja_in)          >= num_nz_in)
SLL_ASSERT(size(a_in)           >= num_nz_in)
SLL_ASSERT(size(ia_out)         >= num_rows_out+1)
SLL_ASSERT(size(ja_out)         >= num_nz_out)
SLL_ASSERT(size(a_out)          >= num_nz_out)
SLL_ASSERT(ia_in(num_rows_in+1) == num_nz_in+1)

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
 
SLL_ASSERT(ia_out(num_rows_out+1) == num_nz_out+1)
  
end subroutine sll_s_csr_add_one_constraint

subroutine sll_s_factorize_csr_matrix(mat)

type(sll_t_csr_matrix), intent(inout) :: mat

#ifdef UMFPACK

sll_real64, dimension(umfpack_info) :: info
  
if (mat%solver == sll_p_umfpack) then

  mat%Ap = mat%row_ptr(:) - 1
  mat%Ai = mat%col_ind(:) - 1
  
  ! pre-order and symbolic analysis
  call umf4sym( mat%num_rows,     &
                mat%num_cols,     &
                mat%Ap,           &
                mat%Ai,           &
                mat%val,          &
                mat%umf_symbolic, &
                mat%umf_control,  &
                info)
  
  if (info(1) .lt. 0) then
    print *, '#Error occurred in umf4sym: ', info(1)
    stop
  endif
  
  ! numeric factorization
  call umf4num( mat%Ap,           &
                mat%Ai,           &
                mat%val,          &
                mat%umf_symbolic, &
                mat%umf_numeric,  &
                mat%umf_control,  &
                info)
      
  if (info(1) .lt. 0) then
    print *, '#Error occurred in umf4num: ', info(1)
    stop
  endif


end if

#else /* UMFPACK */

SLL_ASSERT(associated(mat%val))
 
#endif /* UMFPACK */


end subroutine sll_s_factorize_csr_matrix

  
subroutine sll_s_mult_csr_matrix_vector(mat, input, output)
    
type(sll_t_csr_matrix),     intent(in)  :: mat
sll_real64, dimension(:), intent(in)  :: input
sll_real64, dimension(:), intent(out) :: output

sll_int32 :: i
sll_int32 :: k_1
sll_int32 :: k_2

do i = 1, mat%num_rows

  k_1 = mat%row_ptr(i)
  k_2 = mat%row_ptr(i+1)-1

  output(i) = dot_product(mat%val(k_1:k_2),input(mat%col_ind(k_1:k_2)))
          
end do

end subroutine sll_s_mult_csr_matrix_vector

subroutine sll_s_add_to_csr_matrix(mat, val, row, col)

type(sll_t_csr_matrix), intent(inout) :: mat
sll_real64,             intent(in)    :: val
sll_int32,              intent(in)    :: row
sll_int32,              intent(in)    :: col

sll_int32 :: k

do k = mat%row_ptr(row), mat%row_ptr(row+1) - 1
  if (mat%col_ind(k) == col) then
    mat%val(k) = mat%val(k) + val
    exit
  end if
end do

end subroutine sll_s_add_to_csr_matrix
  
subroutine sll_s_solve_csr_matrix(mat, B, U)

type(sll_t_csr_matrix),   intent(in)    :: mat
sll_real64, dimension(:), intent(inout) :: B
sll_real64, dimension(:), intent(out)   :: U

#ifdef UMFPACK

sll_int32  :: sys
sll_real64, dimension(umfpack_info) :: info

#endif

sll_int32  :: maxIter
sll_real64 :: eps

sll_real64, dimension(:), allocatable :: Ad
sll_real64, dimension(:), allocatable :: d

sll_real64 :: Norm2r1
sll_real64 :: Norm2r0
sll_real64 :: NormInfb
sll_real64 :: NormInfr
sll_real64 :: ps
sll_real64 :: beta
sll_real64 :: alpha
sll_int32  :: iter
sll_int32  :: err

if (mat%solver == sll_p_umfpack) then

#ifdef UMFPACK
sys = 0
call umf4sol(sys,U,B,mat%umf_numeric,mat%umf_control,info)
#endif /* UMFPACK */

else

  eps = 1.d-13
  maxIter = 10000
  
  if ( mat%num_rows /= mat%num_cols ) then
    PRINT*,'#ERROR Gradient_conj: The matrix must be square'
    stop
  end if
  
  if ((abs(maxval(B)) < eps) .AND. (abs(minval(B)) < eps)) then
    U = 0.0_8
    return
  end if
  
  SLL_ALLOCATE(Ad(mat%num_rows),err)
  SLL_ALLOCATE(d(mat%num_rows),err)
  
  U   = 0.0_8
  iter = 0
  
  NormInfb = maxval(abs(B))
  Norm2r0  = dot_product(B,B)
  
  d = B
  
  do iter = 1, maxiter
    !--------------------------------------!
    ! calcul du ak parametre optimal local !
    !--------------------------------------!
  
    call sll_s_mult_csr_matrix_vector( mat , d , Ad )
  
    ps = dot_product( Ad , d )
    alpha = Norm2r0 / ps
            
    !==================================================!
    ! calcul de l'approximation Xk+1 et du residu Rk+1 !
    !==================================================!
    ! calcul des composantes residuelles
    !-----------------------------------
    B = B - alpha * Ad
           
    !----------------------------------------!
    ! approximations ponctuelles au rang k+1 !
    !----------------------------------------!
    U = U + alpha * d
            
    !-------------------------------------------------------!
    ! (a) extraction de la norme infinie du residu          !
    !     pour le test d'arret                              !
    ! (b) extraction de la norme euclidienne du residu rk+1 !
    !-------------------------------------------------------!
    NormInfr = maxval(abs(B))
    Norm2r1  = dot_product(B,B)
  
    !==================================================!
    ! calcul de la nouvelle direction de descente dk+1 !
    !==================================================!
    beta    = Norm2r1 / Norm2r0
    Norm2r0 = Norm2r1
    d      = B + beta * d
           
    !-------------------!
    ! boucle suivante ? !
    !-------------------!
    if ( NormInfr / NormInfb < eps ) exit
            
  end do
  
  if ( iter == maxIter ) then
    print*,'Warning Gradient_conj : iter == maxIter'
    print*,'Error after CG =',( NormInfr / NormInfb )
  end if
    
  deallocate(Ad)
  deallocate(d)

endif 

end subroutine sll_s_solve_csr_matrix

subroutine sll_s_solve_csr_matrix_perper ( mat, B,U,Masse_tot )

type(sll_t_csr_matrix)            :: mat
sll_real64, dimension(:)          :: U
sll_real64, dimension(:)          :: B
sll_real64, dimension(:), pointer :: Masse_tot

#ifdef UMFPACK

sll_real64, dimension(umfpack_info) :: info
sll_int32  :: sys

#endif /* UMFPACK */

sll_real64, dimension(:), pointer :: Ad
sll_real64, dimension(:), pointer :: r
sll_real64, dimension(:), pointer :: d
sll_real64, dimension(:), pointer :: Ux,one
sll_real64 :: Norm2r1
sll_real64 :: Norm2r0
sll_real64 :: NormInfb
sll_real64 :: NormInfr
sll_real64 :: ps
sll_real64 :: beta
sll_real64 :: alpha
sll_int32  :: iter
sll_int32  :: err
logical    :: ll_continue
sll_int32  :: maxIter
sll_real64 :: eps

if (mat%solver == sll_p_umfpack) then

#ifdef UMFPACK
  sys = 0
  call umf4sol(sys,U,B,mat%umf_numeric,mat%umf_control,info)
#endif /* UMFPACK */

else

  maxIter = 100000
  eps = 1.d-13
  
  if ( mat%num_rows /= mat%num_cols ) then
    PRINT*,'ERROR Gradient_conj: The matrix must be square'
    stop
  end if
  
  if ((abs(maxval(B)) < eps ) .AND. (abs(MINVAL(B)) < eps )) then
    U = 0.0_8
    return
  end if
  
  SLL_ALLOCATE(Ad(mat%num_rows),err)
  SLL_ALLOCATE(r(mat%num_rows),err)
  SLL_ALLOCATE(d(mat%num_rows),err)
  SLL_ALLOCATE(Ux(mat%num_rows),err)
  SLL_ALLOCATE(one(mat%num_rows),err)
  
  U(:)  = 0.0_8
  one(:) = 1.0_f64
  Ux(:) = U(:)
  iter = 0
  call sll_s_mult_csr_matrix_vector( mat , Ux , Ad )
  Ad = Ad - dot_product(Masse_tot, Ux)
  r       = B - Ad
  Norm2r0  = DOT_PRODUCT( r , r )
  NormInfb = maxval( abs( B ) )
  
  d = r
  
  ll_continue=.true.
  do while(ll_continue)
    iter = iter + 1
    !--------------------------------------!
    ! calcul du ak parametre optimal local !
    !--------------------------------------!
  
    call sll_s_mult_csr_matrix_vector( mat , d , Ad )
          
    Ad = Ad - dot_product(Masse_tot, d)
    ps = dot_product( Ad , d )
    alpha = Norm2r0 / ps
           
    !==================================================!
    ! calcul de l'approximation Xk+1 et du residu Rk+1 !
    !==================================================!
    ! calcul des composantes residuelles
    !-----------------------------------
    r = r - alpha * Ad
           
    !----------------------------------------!
    ! approximations ponctuelles au rang k+1 !
    !----------------------------------------!
    Ux = Ux + alpha * d
           
    !-------------------------------------------------------!
    ! (a) extraction de la norme infinie du residu          !
    !     pour le test d'arret                              !
    ! (b) extraction de la norme euclidienne du residu rk+1 !
    !-------------------------------------------------------!
    NormInfr = maxval(abs( r ))
    Norm2r1 = dot_product( r , r )
           
    !==================================================!
    ! calcul de la nouvelle direction de descente dk+1 !
    !==================================================!
    beta = Norm2r1 / Norm2r0
    Norm2r0 = Norm2r1
    d = r + beta * d
    !d(1) =  B(1)
    !-------------------!
    ! boucle suivante ? !
    !-------------------!
    ll_continue=((NormInfr/NormInfb) >= eps) .AND. (iter < maxIter)
            
  end do
  U = Ux
   
  if ( iter == maxIter ) then
    print*,'Warning Gradient_conj : iter == maxIter'
    print*,'Error after CG =',( NormInfr / NormInfb )
  end if
    
  deallocate(Ad)
  deallocate(d)
  deallocate(r)
  deallocate(Ux)
  deallocate(one)

endif

end subroutine sll_s_solve_csr_matrix_perper

subroutine sll_s_csr_todense( mat, dense_matrix)

  type(sll_t_csr_matrix)       :: mat
  sll_real64, dimension(:,:) :: dense_matrix
  sll_int32                  :: i, j, k, l

  l = 0
  do i = 1, mat%num_rows 
     do k = mat%row_ptr(i),mat%row_ptr(i+1)-1 
        l = l + 1
        j = mat%col_ind(l)
        dense_matrix(i,j) = mat%val(l)
     end do
  end do

end subroutine sll_s_csr_todense

!> @brief
!> Test function to initialize a CSR matrix
!> @details
!> Fill a matrix in CSR format corresponding to a constant coefficient
!> five-point stencil on a square grid. This function comes from AGMG
!> test program.
subroutine uni2d(mat,f)

type(sll_t_csr_matrix) :: mat
sll_real64           :: f(:)
sll_real64, pointer  :: a(:)
sll_int32            :: m
sll_int32, pointer   :: ia(:)
sll_int32, pointer   :: ja(:)
integer              :: k,l,i,j

real (kind(0d0)), parameter :: zero=0.0d0,cx=-1.0d0,cy=-1.0d0, cd=4.0d0

a  => mat%val
ia => mat%row_ptr
ja => mat%col_ind

m = mat%num_rows

k=0
l=0
ia(1)=1
do i=1,m
  do j=1,m
    k=k+1
    l=l+1
    a(l)=cd
    ja(l)=k
    f(k)=zero
    if(j < m) then
       l=l+1
       a(l)=cx
       ja(l)=k+1
      else
       f(k)=f(k)-cx
    end if
    if(i < m) then
       l=l+1
       a(l)=cy
       ja(l)=k+m
      else
       f(k)=f(k)-cy
    end if
    if(j > 1) then
       l=l+1
       a(l)=cx
       ja(l)=k-1
      else
       f(k)=f(k)-cx
    end if
    if(i >  1) then
       l=l+1
       a(l)=cy
       ja(l)=k-m
      else
       f(k)=f(k)-cy
    end if
    ia(k+1)=l+1
  end do
end do

mat%num_nz = l

return
end subroutine uni2D


end module sll_m_sparse_matrix
