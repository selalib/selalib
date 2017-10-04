!> @ingroup splines
!> @brief   Derived type for periodic banded matrices
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_matrix_periodic_banded

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_spline_matrix_base, only: &
    sll_c_spline_matrix

  use schur_complement, only: &
    schur_complement_solver , &
    schur_complement_fac    , &
    schur_complement_slv    , &
    schur_complement_free

  use iso_fortran_env, only: &
    output_unit

  implicit none

  public :: sll_t_spline_matrix_periodic_banded

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  type, extends( sll_c_spline_matrix ) :: sll_t_spline_matrix_periodic_banded

    integer                       :: n
    integer                       :: kl
    integer                       :: ku
    real(wp), allocatable         :: q(:,:)
    type(schur_complement_solver) :: schur

  contains

    procedure :: init          => s_spline_matrix_periodic_banded__init
    procedure :: set_element   => s_spline_matrix_periodic_banded__set_element
    procedure :: factorize     => s_spline_matrix_periodic_banded__factorize
    procedure :: solve_inplace => s_spline_matrix_periodic_banded__solve_inplace
    procedure :: write         => s_spline_matrix_periodic_banded__write
    procedure :: free          => s_spline_matrix_periodic_banded__free

  end type sll_t_spline_matrix_periodic_banded

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_periodic_banded__init( self, n, kl, ku )
    class(sll_t_spline_matrix_periodic_banded), intent(  out) :: self
    integer                                   , intent(in   ) :: n
    integer                                   , intent(in   ) :: kl
    integer                                   , intent(in   ) :: ku

    character(len=*), parameter :: &
      this_sub_name = "sll_t_spline_matrix_periodic_banded % init"

    SLL_ASSERT( n  >  0 )
    SLL_ASSERT( kl >= 0 )
    SLL_ASSERT( ku >= 0 )
    SLL_ASSERT( ku+1+kl <= n )

    ! FIXME: current linear solver only works for kl=ku
    if ( kl /= ku ) then
      SLL_ERROR( this_sub_name, "Schur complement solver requires kl=ku.")
    end if

    self%n  = n
    self%kl = kl
    self%ku = ku

    ! Allocate matrix Q for linear system:
    !   - M is circulant banded matrix (sparse)
    !   - Q is compressed form of A (circulant diagonals of M are rows of Q)
    allocate( self%q(ku+1+kl,n) )
    self%q(:,:) = 0.0_wp

  end subroutine s_spline_matrix_periodic_banded__init

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_periodic_banded__set_element( self, i, j, a_ij )
    class(sll_t_spline_matrix_periodic_banded), intent(inout) :: self
    integer                                   , intent(in   ) :: i
    integer                                   , intent(in   ) :: j
    real(wp)                                  , intent(in   ) :: a_ij

    integer :: d

    SLL_ASSERT( 1 <= i .and. i <= self%n )
    SLL_ASSERT( 1 <= j .and. j <= self%n )

    ! Compute diagonal index:
    !  . d=0 for main diagonal,
    !  . d>0 for superdiagonals,
    !  . d<0 for subdiagonals
    d = j-i
    if (d >  self%n/2) then; d = d-self%n; else &
    if (d < -self%n/2) then; d = d+self%n; end if

    ! Verify that diagonal index is valid
    SLL_ASSERT( d >= -self%kl )
    SLL_ASSERT( d <=  self%ku )

    ! Write element in correct position
    self%q( self%ku+1-d, j ) = a_ij

  end subroutine s_spline_matrix_periodic_banded__set_element

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_periodic_banded__factorize( self )
    class(sll_t_spline_matrix_periodic_banded), intent(inout) :: self

    ! Prepare solution of linear system by performing:
    ! a) Block decomposition of M = [[A,B],[C,D]] with
    !     - A (n-k)x(n-k) banded Toeplitz matrix (square, large)
    !     - B (n-k)x(  k) low-rank Toeplitz (rectangular, empty in the center)
    !     - C (  k)x(n-k) low-rank Toeplitz (rectangular, empty in the center)
    !     - D (  k)x(  k) fully dense Toeplitz matrix (square, small)
    ! b) LU factorization of matrix A
    ! c) Calculation of Schur complement of A: H = D - C A^{-1} B
    ! d) LU factorization of Schur complement H
    call schur_complement_fac( self%schur, self%n, self%kl, self%q )

  end subroutine s_spline_matrix_periodic_banded__factorize

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_periodic_banded__solve_inplace( self, bx )
    class(sll_t_spline_matrix_periodic_banded), intent(in   ) :: self
    real(wp)                                  , intent(inout) :: bx(:)

    SLL_ASSERT( size(bx)  == self%n )

    call schur_complement_slv( self%schur, self%n, self%kl, self%q, bx )

  end subroutine s_spline_matrix_periodic_banded__solve_inplace

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_periodic_banded__write( self, unit, fmt )
    class(sll_t_spline_matrix_periodic_banded), intent(in) :: self
    integer         ,                 optional, intent(in) :: unit
    character(len=*),                 optional, intent(in) :: fmt

    ! Default file unit and print format
    integer         , parameter :: unit_default = output_unit
    character(len=*), parameter ::  fmt_default = 'es9.1'

    ! Local variables
    integer           :: i, j, d
    integer           :: unit_loc
    character(len=32) :: fmt_loc
    
    ! Overwrite defaults with optional arguments, if present
    if (present( unit )) then
      unit_loc = unit
    else
      unit_loc = unit_default
    end if

    if (present( fmt )) then
      fmt_loc = fmt
    else
      fmt_loc = fmt_default
    end if

    ! Finalize format string
    write(fmt_loc,'(a)') "(" // trim(fmt_loc) // ")"

    ! Write full matrix to selected output stream
    do i = 1, self%n
      do j = 1, self%n
        d = modulo( j-i+self%kl, self%n ) - self%kl
        if (-self%kl<=d .and. d<=self%ku) then
          write(unit_loc,fmt_loc,advance='no') self%q( self%ku+1-d, j )
        else
          write(unit_loc,fmt_loc,advance='no') 0.0_wp
        end if
      end do
      write(unit_loc,*)
    end do

  end subroutine s_spline_matrix_periodic_banded__write

  !-----------------------------------------------------------------------------
  subroutine s_spline_matrix_periodic_banded__free( self )
    class(sll_t_spline_matrix_periodic_banded), intent(inout) :: self

    self%n  = -1
    self%kl = -1
    self%ku = -1
    if (allocated( self%q        )) deallocate( self%q )
    if (allocated( self%schur%bb )) call schur_complement_free( self%schur )

  end subroutine s_spline_matrix_periodic_banded__free

end module sll_m_spline_matrix_periodic_banded
