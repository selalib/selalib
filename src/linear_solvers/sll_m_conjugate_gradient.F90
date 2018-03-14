module sll_m_conjugate_gradient
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_linear_operator_base, only: sll_c_linear_operator

  use sll_m_vector_space_base, only: sll_c_vector_space

  implicit none

  public :: sll_t_conjugate_gradient

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_conjugate_gradient_info

    integer  :: iterations
    logical  :: success
    real(wp) :: residual

  end type sll_t_conjugate_gradient_info

  type :: sll_t_conjugate_gradient

    type(sll_t_conjugate_gradient_info) :: info

  contains

    procedure :: solve => s_conjugate_gradient__solve

  end type sll_t_conjugate_gradient

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Conjugate gradient algorithm for solving linear system Ax = b
  subroutine s_conjugate_gradient__solve( self, A, b, tol, verbose, x )
    class(sll_t_conjugate_gradient), intent(inout) :: self
    class(sll_c_linear_operator)   , intent(in   ) :: A
    class(sll_c_vector_space)      , intent(in   ) :: b
    real(wp)                       , intent(in   ) :: tol
    logical                        , intent(in   ) :: verbose
    class(sll_c_vector_space)      , intent(inout) :: x ! already constructed, first guess

    integer  :: m, n(2)
    real(wp) :: am, am1, l, tol_sqr
    class(sll_c_vector_space), allocatable :: y, p, r, v, t

    ! Get shape of linear operator
    n = A % get_shape()

    ! Some checks
    SLL_ASSERT( n(1) == n(2) )

    ! Construct auxiliary vector spaces
    call x % source( y )
    call x % source( p )
    call x % source( r )
    call x % source( v )
    call x % source( t )

    ! First values
    call A % dot( x, y )     ! y =  Ax
    call y % scal( -1.0_wp ) ! y = -Ax
    call p % add( b, y )     ! p = b - Ax
    call r % copy( p )       ! r = b - Ax
    am = r % inner( r )      ! am = < r, r >

    tol_sqr = tol**2

    ! Iterate to convergence
    do m = 0, n(1)-1

      if ( am < tol_sqr ) exit

      call A % dot( p, v )          ! v = Ap
      l = v % inner( p )            ! l = < v, p >
      l = am / l                    ! l = am / < v, p >
      call t % mult( l, p )         ! t = l*p
      call x % incr( t )            ! x = x + l*p
      call t % mult( -1.0_wp*l, v ) ! t = -l*v
      call r % incr( t )            ! r = r - l*v
      am1 = r % inner( r )          ! am1 = < r, r >
      call p % scal( am1/am )       ! p = am1/am*p
      call p % incr( r )            ! p = r + am1/am*p
      am = am1

    end do

    self % info % iterations = m
    self % info % success    = .true.
    self % info % residual   = sqrt( am )

    if ( verbose ) then

      write(*,*)
      write(*,'(a,i0,a,es8.2)') " CG method converged after ", self % info % iterations, &
                                " iterations with residual " , self % info % residual

    end if

  end subroutine s_conjugate_gradient__solve

end module sll_m_conjugate_gradient
