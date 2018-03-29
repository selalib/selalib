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

  type :: sll_t_conjugate_gradient

    ! Default parameters (can be overwritten from method 'init')
    real(wp) :: tol = 1.0e-14_wp
    logical  :: verbose = .true.

    ! Output information
    integer  :: iterations
    logical  :: success
    real(wp) :: residual

    ! Temporary vector spaces needed by solver
    class(sll_c_vector_space), private, allocatable :: p
    class(sll_c_vector_space), private, allocatable :: r
    class(sll_c_vector_space), private, allocatable :: v

    ! Logical flag telling whether temporary vector spaces are (de)allocated only once
    logical :: allocate_once

  contains

    procedure :: init  => s_conjugate_gradient__init
    procedure :: solve => s_conjugate_gradient__solve
    procedure :: free  => s_conjugate_gradient__free

  end type sll_t_conjugate_gradient

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initializer
  subroutine s_conjugate_gradient__init( self, tol, verbose, template_vector )
    class(sll_t_conjugate_gradient)    , intent(inout) :: self
    real(wp)                 , optional, intent(in   ) :: tol
    logical                  , optional, intent(in   ) :: verbose
    class(sll_c_vector_space), optional, intent(in   ) :: template_vector

    if ( present( tol     ) ) self % tol     = tol
    if ( present( verbose ) ) self % verbose = verbose

    ! Construct temporary vector spaces
    if ( present( template_vector ) ) then

      call template_vector % source( self % p )
      call template_vector % source( self % r )
      call template_vector % source( self % v )

      self % allocate_once = .true.

    else

      self % allocate_once = .false.

    end if

  end subroutine s_conjugate_gradient__init

  ! Conjugate gradient algorithm for solving linear system Ax = b
  subroutine s_conjugate_gradient__solve( self, A, b, x )
    class(sll_t_conjugate_gradient), intent(inout) :: self
    class(sll_c_linear_operator)   , intent(in   ) :: A
    class(sll_c_vector_space)      , intent(in   ) :: b
    class(sll_c_vector_space)      , intent(inout) :: x ! already constructed, first guess

    integer  :: m, n(2)
    real(wp) :: am, am1, l, tol_sqr

    ! Get shape of linear operator
    n = A % get_shape()

    ! Some checks
    SLL_ASSERT( n(1) == n(2) )

    ! Construct temporary vector spaces, if not done already at initialization
    if ( .not. self % allocate_once ) then

      call x % source( self % p )
      call x % source( self % r )
      call x % source( self % v )

    end if

    ! Conjugate gradient algorithm
    associate( p => self % p, r => self % r, v => self % v )

      ! First values
      call A % dot( x, p )     ! p = Ax
      call p % scal( -1.0_wp ) ! p = -Ax
      call p % incr( b )       ! p = b - Ax
      call r % copy( p )       ! r = b - Ax
      am = r % inner( r )      ! am = < r, r >

      tol_sqr = self % tol**2

      ! Iterate to convergence
      do m = 0, n(1)-1

        if ( am < tol_sqr ) exit

        call A % dot( p, v )        ! v = Ap
        l = v % inner( p )          ! l = < v, p >
        l = am / l                  ! l = am / < v, p >
        call x % incr_mult( l, p )  ! x = x + l*p
        call r % incr_mult( -l, v ) ! r = r - l*v
        am1 = r % inner( r )        ! am1 = < r, r >
        call p % scal( am1/am )     ! p = am1/am*p
        call p % incr( r )          ! p = r + am1/am*p
        am = am1

      end do

    end associate

    ! Destroy temporary vector spaces, if not done only once by 'free' method
    if ( .not. self % allocate_once ) then

      call self % p % delete()
      call self % r % delete()
      call self % v % delete()

    end if

    ! Store info
    self % iterations = m
    self % success    = .true.
    self % residual   = sqrt( am )

    ! Verbose output
    if ( self % verbose ) write(*,'(/a,i0,a,es8.2)') " CG method converged after ", &
                          self % iterations, " iterations with residual " , self % residual

  end subroutine s_conjugate_gradient__solve

  ! Free
  subroutine s_conjugate_gradient__free( self )
    class(sll_t_conjugate_gradient), intent(inout) :: self

    if ( self % allocate_once ) then

      call self % p % delete()
      call self % r % delete()
      call self % v % delete()

    end if

  end subroutine s_conjugate_gradient__free

end module sll_m_conjugate_gradient
