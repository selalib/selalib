!> @ingroup splines
!> @brief   Derived type for banded matrices
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_matrix_banded

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

   use sll_m_working_precision, only: f64

   use sll_m_spline_matrix_base, only: &
      sll_c_spline_matrix

   use iso_fortran_env, only: output_unit

   implicit none

   public :: sll_t_spline_matrix_banded

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

   !-----------------------------------------------------------------------------
   !> Banded matrix type
   !-----------------------------------------------------------------------------
   type, extends(sll_c_spline_matrix) :: sll_t_spline_matrix_banded

      integer :: n
      integer :: kl
      integer :: ku
      integer, allocatable :: ipiv(:)
      real(wp), allocatable :: q(:, :)

   contains

      procedure :: init => s_spline_matrix_banded__init
      procedure :: set_element => s_spline_matrix_banded__set_element
      procedure :: factorize => s_spline_matrix_banded__factorize
      procedure :: solve_inplace => s_spline_matrix_banded__solve_inplace
      procedure :: write => s_spline_matrix_banded__write
      procedure :: free => s_spline_matrix_banded__free

   end type sll_t_spline_matrix_banded

   !-----------------------------------------------------------------------------
   ! Interfaces to LAPACK subroutines (www.netlib.org/lapack)
   !-----------------------------------------------------------------------------
   interface

      ! LU factorization of a real M-by-N band matrix A using partial pivoting
      ! with row interchanges
      subroutine dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
         integer, intent(in) :: m             ! no. of rows
         integer, intent(in) :: n             ! no. of columns
         integer, intent(in) :: kl            ! no. of subdiagonals
         integer, intent(in) :: ku            ! no. of superdiagonals
         double precision, intent(inout) :: ab(ldab, n)    ! matrix A in band storage
         integer, intent(in) :: ldab          ! leading dimension of A
         integer, intent(out) :: ipiv(min(m, n))! pivot indices
         integer, intent(out) :: info          ! 0 if successful
      end subroutine dgbtrf

      ! Solution of the linear system A*X=B or A^T*X=B with a general band
      ! matrix A using the LU factorization computed by DGBTRF
      subroutine dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
         character(len=1), intent(in) :: trans         ! form of system of eqns.
         integer, intent(in) :: n             ! order of matrix A
         integer, intent(in) :: kl            ! no. of subdiagonals
         integer, intent(in) :: ku            ! no. of superdiagonals
         integer, intent(in) :: nrhs          ! no. of right hand sides
         double precision, intent(in) :: ab(ldab, n)    ! LU factors (A=P*L*U)
         integer, intent(in) :: ldab          ! leading dimension of A
         integer, intent(in) :: ipiv(n)       ! pivot indices
         double precision, intent(inout) :: b(ldb, nrhs)   ! B on entry, X on exit
         integer, intent(in) :: ldb           ! leading dimension of B
         integer, intent(out) :: info          ! 0 if successful
      end subroutine dgbtrs

   end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !-----------------------------------------------------------------------------
   subroutine s_spline_matrix_banded__init(self, n, kl, ku)
      class(sll_t_spline_matrix_banded), intent(out) :: self
      integer, intent(in) :: n
      integer, intent(in) :: kl
      integer, intent(in) :: ku

      SLL_ASSERT(n > 0)
      SLL_ASSERT(kl >= 0)
      SLL_ASSERT(ku >= 0)
      SLL_ASSERT(kl < n)
      SLL_ASSERT(ku < n)

      self%n = n
      self%kl = kl
      self%ku = ku

      ! Given the linear system A*x=b, we assume that A is a square (n by n)
      ! with ku super-diagonals and kl sub-diagonals.
      ! All non-zero elements of A are stored in the rectangular matrix q, using
      ! the format required by DGBTRF (LAPACK): diagonals of A are rows of q.
      ! q has 2*kl rows for the subdiagonals, 1 row for the diagonal, and ku rows
      ! for the superdiagonals. (The kl additional rows are needed for pivoting.)
      ! The term A(i,j) of the full matrix is stored in q(i-j+2*kl+1,j).

      allocate (self%ipiv(n))
      allocate (self%q(2*kl + ku + 1, n))
      self%q(:, :) = 0.0_wp

   end subroutine s_spline_matrix_banded__init

   !-----------------------------------------------------------------------------
   subroutine s_spline_matrix_banded__set_element(self, i, j, a_ij)
      class(sll_t_spline_matrix_banded), intent(inout) :: self
      integer, intent(in) :: i
      integer, intent(in) :: j
      real(wp), intent(in) :: a_ij

      SLL_ASSERT(max(1, j - self%ku) <= i)
      SLL_ASSERT(i <= min(self%n, j + self%kl))

      self%q(self%kl + self%ku + 1 + i - j, j) = a_ij

   end subroutine s_spline_matrix_banded__set_element

   !-----------------------------------------------------------------------------
   subroutine s_spline_matrix_banded__factorize(self)
      class(sll_t_spline_matrix_banded), intent(inout) :: self

      integer :: info

      character(len=*), parameter :: &
         this_sub_name = "sll_t_spline_matrix_banded % factorize"
      character(len=256) :: err_msg

      ! Perform LU decomposition of matrix q with Lapack
      call dgbtrf(self%n, self%n, self%kl, self%ku, self%q, 2*self%kl + self%ku + 1, &
                  self%ipiv, info)

      if (info < 0) then
         write (err_msg, '(i0,a)') abs(info), "-th argument had an illegal value"
         SLL_ERROR(this_sub_name, err_msg)
      else if (info > 0) then
         write (err_msg, "('U(',i0,',',i0,')',a)") info, info, &
            " is exactly zero. The factorization has been completed, but the factor" &
            //" U is exactly singular, and division by zero will occur if it is used to" &
            //" solve a system of equations."
         SLL_ERROR(this_sub_name, err_msg)
      end if

   end subroutine s_spline_matrix_banded__factorize

   !-----------------------------------------------------------------------------
   subroutine s_spline_matrix_banded__solve_inplace(self, bx)
      class(sll_t_spline_matrix_banded), intent(in) :: self
      real(wp), intent(inout) :: bx(:)

      integer :: info

      character(len=*), parameter :: &
         this_sub_name = "sll_t_spline_matrix_banded % solve_inplace"
      character(len=256) :: err_msg

      SLL_ASSERT(size(bx) == self%n)

      call dgbtrs('N', self%n, self%kl, self%ku, 1, self%q, 2*self%kl + self%ku + 1, &
                  self%ipiv, bx, self%n, info)

      if (info < 0) then
         write (err_msg, '(i0,a)') abs(info), "-th argument had an illegal value"
         SLL_ERROR(this_sub_name, err_msg)
      end if

   end subroutine s_spline_matrix_banded__solve_inplace

   !-----------------------------------------------------------------------------
   subroutine s_spline_matrix_banded__write(self, unit, fmt)
      class(sll_t_spline_matrix_banded), intent(in) :: self
      integer, optional, intent(in) :: unit
      character(len=*), optional, intent(in) :: fmt

      integer :: i, j
      integer :: unit_loc
      character(len=32) :: fmt_loc

      if (present(unit)) then
         unit_loc = unit
      else
         unit_loc = output_unit
      end if

      if (present(fmt)) then
         fmt_loc = fmt
      else
         fmt_loc = 'es9.1'
      end if

      write (fmt_loc, '(a)') "("//trim(fmt_loc)//")"

      do i = 1, self%n
         do j = 1, self%n
            if (max(1, j - self%ku) <= i .and. i <= min(self%n, j + self%kl)) then
               write (unit_loc, fmt_loc, advance='no') self%q(self%kl + self%ku + 1 + i - j, j)
            else
               write (unit_loc, fmt_loc, advance='no') 0.0_wp
            end if
         end do
         write (unit_loc, *)
      end do

   end subroutine s_spline_matrix_banded__write

   !-----------------------------------------------------------------------------
   subroutine s_spline_matrix_banded__free(self)
      class(sll_t_spline_matrix_banded), intent(inout) :: self

      self%n = -1
      self%kl = -1
      self%ku = -1
      if (allocated(self%ipiv)) deallocate (self%ipiv)
      if (allocated(self%q)) deallocate (self%q)

   end subroutine s_spline_matrix_banded__free

end module sll_m_spline_matrix_banded
