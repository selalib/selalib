!> @ingroup splines
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_matrix_dense

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_spline_matrix_base, only: &
    sll_c_spline_matrix

  use iso_fortran_env, only: output_unit

  implicit none

  public :: sll_t_spline_matrix_dense

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  type, extends( sll_c_spline_matrix ) :: sll_t_spline_matrix_dense

    integer :: n
    real(wp), allocatable :: ipiv(:)
    real(wp), allocatable :: a(:,:)

  contains

    procedure :: init          => s_spline_matrix_dense__init
    procedure :: set_element   => s_spline_matrix_dense__set_element
    procedure :: factorize     => s_spline_matrix_dense__factorize
    procedure :: solve_inplace => s_spline_matrix_dense__solve_inplace
    procedure :: write         => s_spline_matrix_dense__write
    procedure :: free          => s_spline_matrix_dense__free

  end type sll_t_spline_matrix_dense

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__init( self, n )
    class(sll_t_spline_matrix_dense), intent(  out) :: self
    integer                         , intent(in   ) :: n

    SLL_ASSERT( n > 0 )

    self%n = n
    allocate( self%ipiv(n) )
    allocate( self%a(n,n) )
    self%a(:,:) = 0.0_wp

  end subroutine s_spline_matrix_dense__init

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__set_element( self, i, j, a_ij )
    class(sll_t_spline_matrix_dense), intent(inout) :: self
    integer                         , intent(in   ) :: i
    integer                         , intent(in   ) :: j
    real(wp)                        , intent(in   ) :: a_ij

    self%a(i,j) = a_ij

  end subroutine s_spline_matrix_dense__set_element

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__factorize( self )
    class(sll_t_spline_matrix_dense), intent(inout) :: self

    integer :: info

    character(len=*), parameter :: &
         this_sub_name = "sll_t_spline_matrix_dense % factorize"
    character(len=256) :: err_msg

    SLL_ASSERT( size(self%a,1) == self%n )
    SLL_ASSERT( size(self%a,2) == self%n )

    call dgetrf( self%n, self%n, self%a, self%n, self%ipiv, info )

    if ( info < 0 ) then
      write( err_msg, '(i0,a)' ) abs(info), "-th argument had an illegal value"
      SLL_ERROR(this_sub_name,err_msg)
    else if ( info > 0 ) then
      write( err_msg, "('U(',i0,',',i0,')',a)" ) info, info, &
           " is exactly zero. The factorization has been completed, but the factor" &
           //" U is exactly singular, and division by zero will occur if it is used to" &
           //" solve a system of equations."
      SLL_ERROR(this_sub_name,err_msg)
    end if

  end subroutine s_spline_matrix_dense__factorize

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__solve_inplace( self, bx )
    class(sll_t_spline_matrix_dense), intent(in   ) :: self
    real(wp)                        , intent(inout) :: bx(:)

    integer :: info

    character(len=*), parameter :: &
         this_sub_name = "sll_t_spline_matrix_dense % solve_inplace"
    character(len=256) :: err_msg

    SLL_ASSERT( size(self%a,1) == self%n )
    SLL_ASSERT( size(self%a,2) == self%n )
    SLL_ASSERT( size(bx)  == self%n )

    call dgetrs( 'N', self%n, 1, self%a, self%n, self%ipiv, bx, self%n, info )

    if ( info < 0 ) then
      write( err_msg, '(i0,a)' ) abs(info), "-th argument had an illegal value"
      SLL_ERROR(this_sub_name,err_msg)
    end if

    end subroutine s_spline_matrix_dense__solve_inplace

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__write( self, unit, fmt )
    class(sll_t_spline_matrix_dense), intent(in) :: self
    integer         , optional      , intent(in) :: unit
    character(len=*), optional      , intent(in) :: fmt

    integer :: i
    integer :: unit_loc
    character(len=32) :: fmt_loc

    if ( present( unit ) ) then
      unit_loc = unit
    else
      unit_loc = output_unit
    end if

    if ( present( fmt  ) ) then
      fmt_loc = fmt
    else
      fmt_loc = 'es9.1'
    end if

    write(fmt_loc,'(a)') "('(',i0,'" // trim(fmt_loc) // ")')"
    write(fmt_loc,fmt_loc) size(self%a,2)

    do i = 1, size(self%a,1)
      write(unit_loc,fmt_loc) self%a(i,:)
    end do

  end subroutine s_spline_matrix_dense__write

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_dense__free( self )
    class(sll_t_spline_matrix_dense), intent(inout) :: self

    self%n = -1
    if ( allocated( self%ipiv ) ) deallocate( self%ipiv )
    if ( allocated( self%a    ) ) deallocate( self%a    )

  end subroutine s_spline_matrix_dense__free

end module sll_m_spline_matrix_dense
