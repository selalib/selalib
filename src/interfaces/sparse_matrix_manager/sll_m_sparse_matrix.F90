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

   use sll_m_qsort_partition, only: sll_s_qsortc

#ifdef UMFPACK
   use sll_m_umfpack
   use iso_fortran_env, only: output_unit
#endif /* UMFPACK */

#ifdef PASTIX
   use sll_m_pastix
#endif /* PASTIX */

#ifdef MUMPS
   use sll_m_mumps
#endif /* MUMPS */

   implicit none

   public :: &
      sll_f_new_csr_matrix, &
      sll_f_new_csr_matrix_with_constraint, &
      sll_s_csr_add_one_constraint, &
      sll_s_csr_todense, &
      sll_s_free_csr_matrix, &
      sll_s_csr_matrix_with_constraint_init, &
      sll_s_add_to_csr_matrix, &
      sll_t_csr_matrix, &
      sll_s_factorize_csr_matrix, &
      sll_s_mult_csr_matrix_vector, &
      sll_s_solve_csr_matrix, &
      sll_s_solve_csr_matrix_perper, &
      sll_p_umfpack, &
      sll_p_pastix, &
      sll_p_mumps

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief type for CSR format
   type :: sll_t_csr_matrix
      sll_int32           :: num_rows !< rows
      sll_int32           :: num_cols !< columns
      sll_int32           :: num_nz   !< non zeros
      sll_int32, pointer :: row_ptr(:)
      sll_int32, pointer :: col_ind(:)
      sll_real64, pointer :: val(:)
      sll_int32           :: solver = 0

#ifdef UMFPACK
      integer(umf_void)                 :: umf_symbolic
      integer(umf_void)                 :: umf_numeric
      sll_real64, dimension(:), pointer :: umf_control
#endif /* UMFPACK */

#ifdef PASTIX
      type(pastix_solver)               :: pstx
#endif
#ifdef MUMPS
      type(mumps_solver)                :: mmps
#endif

   end type sll_t_csr_matrix

   sll_int32, parameter :: sll_p_umfpack = 1
   sll_int32, parameter :: sll_p_pastix = 2
   sll_int32, parameter :: sll_p_mumps = 3

   interface sll_f_new_csr_matrix
      module procedure new_csr_matrix_with_dof
   end interface sll_f_new_csr_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine sll_s_free_csr_matrix(mat)
      type(sll_t_csr_matrix), pointer :: mat
      if (associated(mat)) nullify (mat)
   end subroutine sll_s_free_csr_matrix

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
   function new_csr_matrix_with_dof( &
      num_rows, &
      num_cols, &
      num_elts, &
      local_to_global_row, &
      num_local_dof_row, &
      local_to_global_col, &
      num_local_dof_col, &
      solver) result(mat)

      type(sll_t_csr_matrix), pointer    :: mat
      sll_int32, intent(in) :: num_rows
      sll_int32, intent(in) :: num_cols
      sll_int32, intent(in) :: num_elts
      sll_int32, dimension(:, :), intent(in) :: local_to_global_row
      sll_int32, dimension(:, :), intent(in) :: local_to_global_col
      sll_int32, intent(in) :: num_local_dof_row
      sll_int32, intent(in) :: num_local_dof_col
      sll_int32, optional   :: solver

      sll_int32 :: ierr

      SLL_ALLOCATE(mat, ierr)

      mat%solver = 0
      if (present(solver)) then
#ifdef UMFPACK
         if (solver == sll_p_umfpack) mat%solver = solver
#endif /* UMFPACK */
#ifdef PASTIX
         if (solver == sll_p_pastix) mat%solver = solver
#endif /* PASTIX */
#ifdef MUMPS
         if (solver == sll_p_mumps) mat%solver = solver
#endif /* MUMPS */
      end if

      call sll_s_csr_matrix_init( &
         mat, &
         num_rows, &
         num_cols, &
         num_elts, &
         local_to_global_row, &
         num_local_dof_row, &
         local_to_global_col, &
         num_local_dof_col)

   end function new_csr_matrix_with_dof

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

   subroutine sll_s_csr_matrix_init(mat, &
                                    num_rows, &
                                    num_cols, &
                                    num_elts, &
                                    local_to_global_row, &
                                    num_local_dof_row, &
                                    local_to_global_col, &
                                    num_local_dof_col)

      type(sll_t_csr_matrix), intent(inout) :: mat
      sll_int32, intent(in)    :: num_rows
      sll_int32, intent(in)    :: num_cols
      sll_int32, intent(in)    :: num_elts
      sll_int32, dimension(:, :), intent(in)    :: local_to_global_row
      sll_int32, intent(in)    :: num_local_dof_row
      sll_int32, dimension(:, :), intent(in)    :: local_to_global_col
      sll_int32, intent(in)    :: num_local_dof_col

      sll_int32                                :: num_nz
      sll_int32                                :: ierr
      sll_int32, dimension(:, :), allocatable  :: lpi_col
      sll_int32, dimension(:), allocatable  :: lpi_occ
      sll_int32                                :: coef
      sll_int32                                :: elt
      sll_int32                                :: ii
      sll_int32                                :: jj
      sll_int32                                :: row
      sll_int32                                :: col
      sll_int32                                :: i
      sll_int32                                :: sz
      logical                                  :: ll_done

#ifdef DEBUG
      print *, '#sll_s_csr_matrix_init'
#endif

      coef = 6

      SLL_ALLOCATE(lpi_col(num_rows, 0:coef*num_local_dof_col), ierr)
      SLL_ALLOCATE(lpi_occ(num_rows + 1), ierr)

      lpi_col(:, :) = 0
      lpi_occ(:) = 0

      do elt = 1, num_elts  !Loop over cells
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
                     if (.not. ll_done) then
                        lpi_occ(row) = lpi_occ(row) + 1
                        lpi_col(row, 0) = lpi_col(row, 0) + 1
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
      mat%num_nz = num_nz

#ifdef DEBUG
      print *, '#solver  =', mat%solver
      print *, '#num_rows=', num_rows
      print *, '#num_cols=', num_cols
      print *, '#num_nz  =', num_nz
#endif

      SLL_ALLOCATE(mat%row_ptr(num_rows + 1), ierr)
      SLL_ALLOCATE(mat%col_ind(num_nz), ierr)
      SLL_ALLOCATE(mat%val(num_nz), ierr)

      mat%row_ptr(1) = 1

      do i = 1, mat%num_rows
         mat%row_ptr(i + 1) = mat%row_ptr(1) + sum(lpi_occ(1:i))
      end do

      do elt = 1, num_elts
         do ii = 1, num_local_dof_row
            row = local_to_global_row(ii, elt)
            if (row /= 0) then
               if (lpi_col(row, 0) /= 0) then
                  sz = lpi_col(row, 0)
                  call sll_s_qsortc(lpi_col(row, 1:sz))
                  do i = 1, sz
                     mat%col_ind(mat%row_ptr(row) + i - 1) = lpi_col(row, i)
                  end do
                  lpi_col(row, 0) = 0
               end if
            end if
         end do
      end do

      mat%val(:) = 0.0_f64
      SLL_DEALLOCATE_ARRAY(lpi_col, ierr)
      SLL_DEALLOCATE_ARRAY(lpi_occ, ierr)

      select case (mat%solver)

      case (sll_p_umfpack)

#ifdef UMFPACK
         SLL_ALLOCATE(mat%umf_control(umfpack_control), ierr)
         mat%row_ptr = mat%row_ptr - 1
         mat%col_ind = mat%col_ind - 1
         call umf4def(mat%umf_control)  ! get the default configuration
#endif /* UMFPACK */

      case (sll_p_pastix)

#ifdef PASTIX
         call initialize(mat%pstx, num_rows, num_nz, mat%row_ptr, mat%col_ind, mat%val)
#endif /* PASTIX */

      case (sll_p_mumps)

#ifdef MUMPS
         call initialize(mat%mmps, num_rows, num_nz, mat%row_ptr, mat%col_ind, mat%val)
#endif /* MUMPS */

      end select

   end subroutine sll_s_csr_matrix_init

   subroutine sll_s_csr_matrix_with_constraint_init(mat, mat_a)

      type(sll_t_csr_matrix), intent(inout) :: mat
      type(sll_t_csr_matrix), intent(in) :: mat_a
      sll_int32 :: ierr

      mat%num_nz = mat_a%num_nz + 2*mat_a%num_rows
      mat%num_rows = mat_a%num_rows + 1
      mat%num_cols = mat_a%num_cols + 1
#ifdef DEBUG
      print *, '# sll_s_csr_matrix_with_constraint_init'
      print *, '# num_nz mat, num_nz mat_tot', mat_a%num_nz, mat%num_nz
      print *, '# num_rows mat, num_rows mat_tot', mat_a%num_rows, mat%num_rows
      print *, '# num_cols mat, num_cols mat_tot', mat_a%num_cols, mat%num_cols
#endif /* DEBUG */

      SLL_ALLOCATE(mat%row_ptr(mat%num_rows + 1), ierr)
      SLL_ALLOCATE(mat%col_ind(mat%num_nz), ierr)
      SLL_CLEAR_ALLOCATE(mat%val(1:mat%num_nz), ierr)

      mat%solver = mat_a%solver

      select case (mat%solver)
#ifdef UMFPACK
      case (sll_p_umfpack)
         SLL_ALLOCATE(mat%umf_control(umfpack_control), ierr)
         call umf4def(mat%umf_control)
#endif /* UMFPACK */
      case (sll_p_pastix)
#ifdef PASTIX
         call initialize(mat%pstx, mat%num_rows, mat%num_nz, mat%row_ptr, &
                         mat%col_ind, mat%val)
#endif /* PASTIX */
      case (sll_p_mumps)
#ifdef MUMPS
         call initialize(mat%mmps, mat%num_rows, mat%num_nz, mat%row_ptr, mat%col_ind, mat%val)
#endif /* MUMPS */
      end select

   end subroutine sll_s_csr_matrix_with_constraint_init

   function sll_f_new_csr_matrix_with_constraint(mat_a) result(mat)

      type(sll_t_csr_matrix), pointer :: mat
      type(sll_t_csr_matrix)          :: mat_a

      sll_int32 :: ierr
      SLL_ALLOCATE(mat, ierr)
      call sll_s_csr_matrix_with_constraint_init(mat, mat_a)

   end function sll_f_new_csr_matrix_with_constraint

   subroutine sll_s_csr_add_one_constraint(ia_in, &
                                           ja_in, &
                                           a_in, &
                                           num_rows_in, &
                                           num_nz_in, &
                                           constraint_vec, &
                                           ia_out, &
                                           ja_out, &
                                           a_out)

      integer, dimension(:), intent(in)  :: ia_in
      integer, dimension(:), intent(in)  :: ja_in
      real(8), dimension(:), intent(in)  :: a_in
      integer, intent(in)  :: num_rows_in
      integer, intent(in)  :: num_nz_in
      real(8), dimension(:), intent(in)  :: constraint_vec
      integer, dimension(:), intent(out) :: ia_out
      integer, dimension(:), intent(out) :: ja_out
      real(8), dimension(:), intent(out) :: a_out

      integer :: num_rows_out
      integer :: num_nz_out
      integer :: i
      integer :: s
      integer :: k

      num_rows_out = num_rows_in + 1
      num_nz_out = num_nz_in + 2*num_rows_in

      SLL_ASSERT(size(ia_in) >= num_rows_in + 1)
      SLL_ASSERT(size(ja_in) >= num_nz_in)
      SLL_ASSERT(size(a_in) >= num_nz_in)
      SLL_ASSERT(size(ia_out) >= num_rows_out + 1)
      SLL_ASSERT(size(ja_out) == num_nz_out)
      SLL_ASSERT(size(a_out) >= num_nz_out)

      if (ia_in(num_rows_in + 1) == num_nz_in) then ! UMFPACK

         s = 1
         do i = 1, num_rows_in
            ia_out(i) = s - 1
            do k = ia_in(i) + 1, ia_in(i + 1)
               a_out(s) = a_in(k)
               ja_out(s) = ja_in(k)
               s = s + 1
            end do
            a_out(s) = constraint_vec(i)
            ja_out(s) = num_rows_out - 1
            s = s + 1
         end do
         ia_out(num_rows_in + 1) = s - 1
         do i = 1, num_rows_in
            a_out(s) = constraint_vec(i)
            ja_out(s) = i - 1
            s = s + 1
         end do
         ia_out(num_rows_in + 2) = s - 1

         SLL_ASSERT(ia_out(num_rows_out + 1) == num_nz_out)

      else ! Default CG

         SLL_ASSERT(ia_in(num_rows_in + 1) == num_nz_in + 1)

         s = 1
         do i = 1, num_rows_in
            ia_out(i) = s
            do k = ia_in(i), ia_in(i + 1) - 1
               a_out(s) = a_in(k)
               ja_out(s) = ja_in(k)
               s = s + 1
            end do
            a_out(s) = constraint_vec(i)
            ja_out(s) = num_rows_out
            s = s + 1
         end do
         ia_out(num_rows_in + 1) = s
         do i = 1, num_rows_in
            a_out(s) = constraint_vec(i)
            ja_out(s) = i
            s = s + 1
         end do
         ia_out(num_rows_in + 2) = s

         SLL_ASSERT(ia_in(num_rows_in + 1) == num_nz_in + 1)

      end if

   end subroutine sll_s_csr_add_one_constraint

   subroutine sll_s_factorize_csr_matrix(mat)

      type(sll_t_csr_matrix), intent(inout) :: mat

#ifdef UMFPACK
      sll_real64, dimension(umfpack_info) :: info
#endif /* UMFPACK */
#ifdef MUMPS
      sll_int32 :: i, k, l
#endif /* MUMPS */

#ifdef DEBUG
      call print_solver_type(mat)
      print *, '# begin factorize_csr_matrix '
#endif /* DEBUG */

      select case (mat%solver)

      case (sll_p_umfpack)

#ifdef UMFPACK
         ! pre-order and symbolic analysis
         call umf4sym(mat%num_rows, &
                      mat%num_cols, &
                      mat%row_ptr, &
                      mat%col_ind, &
                      mat%val, &
                      mat%umf_symbolic, &
                      mat%umf_control, &
                      info)

         if (info(1) .lt. 0) then
            print *, '#Error occurred in umf4sym: ', info(1)
            stop
         end if

         ! numeric factorization
         call umf4num(mat%row_ptr, &
                      mat%col_ind, &
                      mat%val, &
                      mat%umf_symbolic, &
                      mat%umf_numeric, &
                      mat%umf_control, &
                      info)

         if (info(1) .lt. 0) then
            print *, '#Error occurred in umf4num: ', info(1)
            stop
         end if
#endif /* UMFPACK */

      case (sll_p_pastix)

#ifdef PASTIX
         call factorize(mat%pstx)
#endif /* PASTIX */

      case (sll_p_mumps)

#ifdef MUMPS
         l = 0
         do i = 1, mat%num_rows
            do k = mat%row_ptr(i), mat%row_ptr(i + 1) - 1
               l = l + 1
               mat%mmps%mumps_par%irn(l) = i
               mat%mmps%mumps_par%jcn(l) = mat%col_ind(l)
               mat%mmps%mumps_par%a(l) = mat%val(l)
            end do
         end do
         call factorize(mat%mmps)
#endif /* MUMPS */

      case default

         SLL_ASSERT(associated(mat%val)) ! Just to avoid warning in debug mode

      end select

#ifdef DEBUG
      print *, '# end factorize_csr_matrix '
#endif /* DEBUG */

   end subroutine sll_s_factorize_csr_matrix

   subroutine sll_s_mult_csr_matrix_vector(mat, input, output)

      type(sll_t_csr_matrix), intent(in)  :: mat
      sll_real64, dimension(:), intent(in)  :: input
      sll_real64, dimension(:), intent(out) :: output

      sll_int32 :: i
      sll_int32 :: k_1
      sll_int32 :: k_2

      do i = 1, mat%num_rows

         k_1 = mat%row_ptr(i)
         k_2 = mat%row_ptr(i + 1) - 1

         output(i) = dot_product(mat%val(k_1:k_2), input(mat%col_ind(k_1:k_2)))

      end do

   end subroutine sll_s_mult_csr_matrix_vector

   subroutine sll_s_add_to_csr_matrix(mat, val, row, col)

      type(sll_t_csr_matrix), intent(inout) :: mat
      sll_real64, intent(in)    :: val
      sll_int32, intent(in)    :: row
      sll_int32, intent(in)    :: col

      sll_int32 :: k

      if (mat%solver == sll_p_umfpack) then
         do k = mat%row_ptr(row) + 1, mat%row_ptr(row + 1)
            if (mat%col_ind(k) + 1 == col) then
               mat%val(k) = mat%val(k) + val
               exit
            end if
         end do
      else
         do k = mat%row_ptr(row), mat%row_ptr(row + 1) - 1
            if (mat%col_ind(k) == col) then
               mat%val(k) = mat%val(k) + val
               exit
            end if
         end do
      end if

   end subroutine sll_s_add_to_csr_matrix

   subroutine sll_s_solve_csr_matrix(mat, b, u)

      type(sll_t_csr_matrix), intent(in)    :: mat
      sll_real64, dimension(:), intent(inout) :: b
      sll_real64, dimension(:), intent(out)   :: u

#ifdef UMFPACK

      sll_int32                           :: sys
      sll_real64, dimension(umfpack_info) :: info

#endif

#ifdef DEBUG
      print *, '# begin solve_csr_matrix '
#endif /* DEBUG */

      select case (mat%solver)

      case (sll_p_umfpack)

#ifdef UMFPACK
         sys = 0
         call umf4sol(sys, u, b, mat%umf_numeric, mat%umf_control, info)
#endif /* UMFPACK */

      case (sll_p_pastix)

#ifdef PASTIX
         u = b
         call solve(mat%pstx, u)
#endif /* PASTIX */

      case (sll_p_mumps)

#ifdef MUMPS
         u = b
         call solve(mat%mmps, u)
#endif /* MUMPS */

      case default

         call solve_with_cg(mat, b, u)

      end select

#ifdef DEBUG
      print *, '# end solve_csr_matrix '
#endif /* DEBUG */

   end subroutine sll_s_solve_csr_matrix

   subroutine sll_s_solve_csr_matrix_perper(mat, b, u, masse_tot)

      type(sll_t_csr_matrix)            :: mat
      sll_real64, dimension(:)          :: u
      sll_real64, dimension(:)          :: b
      sll_real64, dimension(:), pointer :: masse_tot

#ifdef UMFPACK

      sll_real64, dimension(umfpack_info) :: info
      sll_int32  :: sys

#endif /* UMFPACK */

      select case (mat%solver)

      case (sll_p_umfpack)

#ifdef UMFPACK
         sys = 0
         call umf4sol(sys, u, b, mat%umf_numeric, mat%umf_control, info)
#endif /* UMFPACK */

      case (sll_p_pastix)

#ifdef PASTIX
         u = b
         call solve(mat%pstx, u)
#endif /* PASTIX */

      case (sll_p_mumps)

#ifdef MUMPS
         u = b
         call solve(mat%mmps, u)
#endif /* MUMPS */

      case default

         call solve_with_cg(mat, b, u)

      end select

   end subroutine sll_s_solve_csr_matrix_perper

   subroutine sll_s_csr_todense(mat, dense_matrix)

      type(sll_t_csr_matrix)     :: mat
      sll_real64, dimension(:, :) :: dense_matrix
      sll_int32                  :: i, j, k, l

      if (mat%solver == sll_p_umfpack) then
         l = 0
         do i = 1, mat%num_rows
            do k = mat%row_ptr(i) + 1, mat%row_ptr(i + 1)
               l = l + 1
               j = mat%col_ind(l) + 1
               dense_matrix(i, j) = mat%val(l)
            end do
         end do
      else
         l = 0
         do i = 1, mat%num_rows
            do k = mat%row_ptr(i), mat%row_ptr(i + 1) - 1
               l = l + 1
               j = mat%col_ind(l)
               dense_matrix(i, j) = mat%val(l)
            end do
         end do
      end if

   end subroutine sll_s_csr_todense

   subroutine solve_with_cg(mat, b, u)

      type(sll_t_csr_matrix)   :: mat
      sll_real64, dimension(:) :: u
      sll_real64, dimension(:) :: b

      sll_int32  :: maxIter
      sll_real64 :: eps

      sll_real64, dimension(:), allocatable :: ad
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

      eps = 1.d-13
      maxIter = 10000

      if (mat%num_rows /= mat%num_cols) then
         print *, '#ERROR Gradient_conj: The matrix must be square'
         stop
      end if

      if ((abs(maxval(b)) < eps) .AND. (abs(minval(b)) < eps)) then
         u = 0.0_f64
         return
      end if

      SLL_ALLOCATE(ad(mat%num_rows), err)
      SLL_ALLOCATE(d(mat%num_rows), err)

      u = 0.0_f64
      iter = 0

      NormInfb = maxval(abs(b))
      Norm2r0 = dot_product(b, b)

      d = b

      do iter = 1, maxiter
         !--------------------------------------!
         ! calcul du ak parametre optimal local !
         !--------------------------------------!

         call sll_s_mult_csr_matrix_vector(mat, d, ad)

         ps = dot_product(ad, d)
         alpha = Norm2r0/ps

         !==================================================!
         ! calcul de l'approximation Xk+1 et du residu Rk+1 !
         !==================================================!
         ! calcul des composantes residuelles
         !-----------------------------------
         b = b - alpha*ad

         !----------------------------------------!
         ! approximations ponctuelles au rang k+1 !
         !----------------------------------------!
         u = u + alpha*d

         !-------------------------------------------------------!
         ! (a) extraction de la norme infinie du residu          !
         !     pour le test d'arret                              !
         ! (b) extraction de la norme euclidienne du residu rk+1 !
         !-------------------------------------------------------!
         NormInfr = maxval(abs(b))
         Norm2r1 = dot_product(b, b)

         !==================================================!
         ! calcul de la nouvelle direction de descente dk+1 !
         !==================================================!
         beta = Norm2r1/Norm2r0
         Norm2r0 = Norm2r1
         d = b + beta*d

         !-------------------!
         ! boucle suivante ? !
         !-------------------!
         if (NormInfr/NormInfb < eps) exit

      end do

      if (iter == maxIter) then
         print *, 'Warning Gradient_conj : iter == maxIter'
         print *, 'Error after CG =', (NormInfr/NormInfb)
      end if

      deallocate (ad)
      deallocate (d)

   end subroutine solve_with_cg

   subroutine print_solver_type(mat)

      type(sll_t_csr_matrix)   :: mat

      print *, 'print_solver  =', mat%solver
      select case (mat%solver)
      case (sll_p_umfpack)
         print *, "# We use UMFPACK to solve the system "
      case (sll_p_pastix)
         print *, "# We use PASTIX to solve the system "
      case (sll_p_mumps)
         print *, "# We use MUMPS to solve the system "
      case default
         print *, "# We use CG to solve the system "
      end select

   end subroutine print_solver_type

end module sll_m_sparse_matrix
