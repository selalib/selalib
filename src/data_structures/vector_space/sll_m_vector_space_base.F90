!> @ingroup vector_space
!> @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!> @authors Marco Restelli - <marco.restelli@gmail.com>
!> @authors Edoardo Zoni   - <edoardo.zoni@ipp.mpg.de>
!> @brief   Abstract type implementing a generic vector space
!
module sll_m_vector_space_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  implicit none

  public :: sll_c_vector_space

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !-----------------------------------------------------------------------------
  ! Abstract type definition
  !-----------------------------------------------------------------------------

  !> @brief Abstract base class for all vector spaces
  type, abstract :: sll_c_vector_space

  contains

    !> @name Basic operations (abstract methods)
    !> Any non-abstract extended type MUST provide an implementation for these.
    procedure( i_copy ), deferred :: copy ! z = x
    procedure( i_incr ), deferred :: incr ! z+= x
    procedure( i_scal ), deferred :: scal ! z*= a

    !> @name Additional operations
    !> We provide a standard implementation that uses
    !> the operations above and creates temporary objects where needed.
    !> The user should provide a more efficient implementation.
    procedure :: add       => s_vector_space__add       ! z = x + y
    procedure :: mult      => s_vector_space__mult      ! z = a * x
    procedure :: mult_add  => s_vector_space__mult_add  ! z = a * x + y
    procedure :: incr_mult => s_vector_space__incr_mult ! z+= a * x
    procedure :: lcmb      => s_vector_space__lcmb      ! z = sum_i ( a_i*x_i )
    procedure :: incr_lcmb => s_vector_space__incr_lcmb ! z+= sum_i ( a_i*x_i )

    !> @name Optional subroutines and functions
    !> Provide a *norm* to have a normed vector space.
    !> Provide an *inner product* to have an inner product space.
    !> Provide a *show* function for debugging purposes.
    procedure :: norm  => f_vector_space__norm  !< norm(z)
    procedure :: inner => f_vector_space__inner !< inner(z,x)
    procedure :: show  => s_vector_space__show

    !> @name Constructors & destructors
    generic :: source => source_scalar, source_array
    !> @}

    !> @name Copy constructors (one vector or an array of vectors)
    procedure, private :: source_scalar => s_vector_space__source_scalar
    procedure, private :: source_array  => s_vector_space__source_array
    !> @}

  end type sll_c_vector_space

  !-----------------------------------------------------------------------------
  ! Abstract interfaces: basic operations
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> @brief        z = x
  !> @details      Copy contents of input array *x* into array *z* ('self').
  !>               Usage: \code z%copy(x) \endcode
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    x    Input vector
  !-----------------------------------------------------------------------------
  abstract interface
    subroutine i_copy( self, x )
      import sll_c_vector_space
      class( sll_c_vector_space ), intent(inout) :: self
      class( sll_c_vector_space ), intent(in   ) :: x
    end subroutine i_copy
  end interface

  !-----------------------------------------------------------------------------
  !> @brief        z += x
  !> @details      Increment the vector *z* ('self') by the vector *x*.
  !>               Usage: \code z%incr(x) \endcode
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    x    Input vector
  !-----------------------------------------------------------------------------
  abstract interface
    subroutine i_incr( self, x )
      import sll_c_vector_space
      class( sll_c_vector_space ), intent(inout) :: self
      class( sll_c_vector_space ), intent(in   ) :: x
    end subroutine i_incr
  end interface

  !-----------------------------------------------------------------------------
  !> @brief        z *= a
  !> @details      Scale the vector *z* ('self') by the real number *a*.
  !>               Usage: \code z%scal(a) \endcode
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Input scalar, double precision real
  !-----------------------------------------------------------------------------
  abstract interface
    subroutine i_scal( self, a )
      import sll_c_vector_space, wp
      class( sll_c_vector_space ), intent(inout) :: self
      real(wp)                   , intent(in   ) :: a
    end subroutine i_scal
  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  ! Default subroutines: composite operations
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> @brief        z = x + y
  !> @details      Add two vectors *x* and *y*, and write the result to the
  !>               vector *z* ('self') that invoked this subroutine.
  !>               Usage: \code call z%add(x,y) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    x    Input vector #1
  !> @param[in]    y    Input vector #2
  !-----------------------------------------------------------------------------
  subroutine s_vector_space__add( self, x, y )
    class(sll_c_vector_space), intent(inout) :: self
    class(sll_c_vector_space), intent(in   ) :: x
    class(sll_c_vector_space), intent(in   ) :: y

    call self % copy( x )
    call self % incr( y )

  end subroutine s_vector_space__add

  !-----------------------------------------------------------------------------
  !> @brief        z = a * x
  !> @details      Multiply the vector *x* by the scalar *a*, and assign the
  !>               result to the vector *z* ('self') that invoked this
  !>               subroutine.
  !>               Usage: \code call z%mult(a,x) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Input scalar
  !> @param[in]    x    Input vector
  !-----------------------------------------------------------------------------
  subroutine s_vector_space__mult( self, a, x )
    class(sll_c_vector_space), intent(inout) :: self
    real(wp)                 , intent(in   ) :: a
    class(sll_c_vector_space), intent(in   ) :: x

    call self % copy( x )
    call self % scal( a )

  end subroutine s_vector_space__mult

  !-----------------------------------------------------------------------------
  !> @brief        z = a * x + y
  !> @details      Multiply the vector *x* by the scalar *a*, sum with the
  !>               vector *y*, and assign the result to the vector *z* ('self')
  !>               that invoked this subroutine.
  !>               Usage: \code call z%mult_add(a,x,y) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Input scalar
  !> @param[in]    x    Input vector #1
  !> @param[in]    y    Input vector #2
  !-----------------------------------------------------------------------------
  subroutine s_vector_space__mult_add( self, a, x, y )
    class(sll_c_vector_space), intent(inout) :: self
    real(wp)                 , intent(in   ) :: a
    class(sll_c_vector_space), intent(in   ) :: x
    class(sll_c_vector_space), intent(in   ) :: y
    
    call self % mult( a, x )
    call self % incr( y )

  end subroutine s_vector_space__mult_add

  !-----------------------------------------------------------------------------
  !> @brief        z += a * x
  !> @details      Increment the vector *z* ('self') by the result of the
  !>               multiplication between the vector *x* and the scalar *a*.
  !>               Usage: \code call z%incr_mult(a,x) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Input scalar
  !> @param[in]    x    Input vector
  !-----------------------------------------------------------------------------
  subroutine s_vector_space__incr_mult( self, a, x )
    class(sll_c_vector_space), intent(inout) :: self
    real(wp)                 , intent(in   ) :: a
    class(sll_c_vector_space), intent(in   ) :: x
    
    class(sll_c_vector_space), allocatable :: temp

    call self % source( temp )
    call temp % mult( a, x )
    call self % incr( temp )

    deallocate( temp )

  end subroutine s_vector_space__incr_mult

  !-----------------------------------------------------------------------------
  !> @brief        z = \f$ \sum_i ( a_i x_i ) \f$
  !> @details      Compute a linear combination of N vectors \f$ x_i \f$,
  !>               according to N coefficients \f$ a_i \f$, and assign the
  !>               result to the vector *z* that invoked this subroutine.
  !>               Usage: \code call z%lcmb(a,x) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Array of N scalars
  !> @param[in]    x    Array of N vectors
  !-----------------------------------------------------------------------------
  subroutine s_vector_space__lcmb( self, a, x )
    class(sll_c_vector_space), intent(inout) :: self
    real(wp)                 , intent(in   ) :: a(:)
    class(sll_c_vector_space), intent(in   ) :: x(:)

    class(sll_c_vector_space), allocatable :: temp
    integer :: i

    call self % source( temp )
    call self % mult( a(1), x(1) )

    do i = 2, size(a)
      call temp % mult( a(i), x(i) )
      call self % incr( temp )
    end do

    deallocate( temp )
 
  end subroutine s_vector_space__lcmb

  !-----------------------------------------------------------------------------
  !> @brief        z += \f$ \sum_i ( a_i x_i ) \f$
  !> @details      Compute a linear combination of N vectors \f$ x_i \f$,
  !>               according to N coefficients \f$ a_i \f$, and use the result
  !>               to increment the vector *z* that invoked this subroutine.
  !>               Usage: \code call z%incr_lcmb(a,x) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Array of N scalars
  !> @param[in]    x    Array of N vectors
  !-----------------------------------------------------------------------------
  subroutine s_vector_space__incr_lcmb( self, a, x )
    class(sll_c_vector_space), intent(inout) :: self
    real(wp)                 , intent(in   ) :: a(:)
    class(sll_c_vector_space), intent(in   ) :: x(:)

    class(sll_c_vector_space), allocatable :: temp
    integer :: i

    call self % source( temp )

    do i = 1, size(a)
      call temp % mult( a(i), x(i) )
      call self % incr( temp )
    end do

    deallocate( temp )

  end subroutine s_vector_space__incr_lcmb

  !-----------------------------------------------------------------------------
  ! Empty functions/subroutines: optional operations and debugging
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> @brief     Norm of vector: ||z||
  !> @details   Compute the norm of the vector *z* ('self') that is invoking
  !>            the function. By default, \c norm(z)=sqrt(inner(z,z)).
  !>            Usage: \code r = z%norm() \endcode
  !>
  !> @param[in] self Vector *z*, caller
  !> @returns        Scalar value, of type *real(wp)*
  !-----------------------------------------------------------------------------
  function f_vector_space__norm( self ) result( res )
    class(sll_c_vector_space), intent(in) :: self
    real(wp) :: res

    res = sqrt( self%inner( self ) )

  end function f_vector_space__norm

  !-----------------------------------------------------------------------------
  !> @brief     Inner product: \<z,x\>
  !> @details   Compute the inner product between the vector *z* ('self') and
  !>            another vector *x* of the same type.
  !>            Usage: \code r = z%inner(x) \endcode
  !>
  !> @param[in] self Vector *z*, caller
  !> @param[in] x    Vector
  !> @returns        Scalar value, of type *real(wp)*
  !-----------------------------------------------------------------------------
  function f_vector_space__inner( self, x ) result( res )
    class(sll_c_vector_space), intent(in) :: self
    class(sll_c_vector_space), intent(in) :: x
    real(wp) :: res

    SLL_ERROR( "sll_c_vector_space % inner", "Function not implemented." )

#ifdef DEBUG
    print*, storage_size( self ), storage_size( x )
#endif

    res = 0.0_f64

  end function f_vector_space__inner

  !-----------------------------------------------------------------------------
  !> @brief     Show something, for debug
  !> @param[in] self Vector *z*, caller
  !-----------------------------------------------------------------------------
  subroutine s_vector_space__show( self )
    class( sll_c_vector_space ), intent( in ) :: self

    SLL_WARNING( "sll_c_vector_space % show", "Overload this subroutine if you need it." )

    print*, storage_size( self )

  end subroutine s_vector_space__show

  !-----------------------------------------------------------------------------
  !> @brief      Copy constructor: create one copy of vector *z*
  !> @details    Allocate vector *x* and copy contents of *z* into it. Usage:
  !>             \code
  !>               class( sll_c_vector_space ), allocatable :: x
  !>               call z%source(x) 
  !>             \endcode
  !> @param[in]  self Vector *z*, caller
  !> @param[out] x    Output copy (allocatable)
  !-----------------------------------------------------------------------------
  subroutine s_vector_space__source_scalar( self, x )
    class(sll_c_vector_space),              intent(in   ) :: self
    class(sll_c_vector_space), allocatable, intent(  out) :: x

    allocate( x, source=self )

  end subroutine s_vector_space__source_scalar

  !-----------------------------------------------------------------------------
  !> @brief      Copy constructor: create *n* copies of vector *z*
  !> @details    Allocate array *x* of *n* vectors, and copy contents of *z*
  !>             into each element *x(i)*. Usage:
  !>             \code
  !>               class( sll_c_vector_space ), allocatable :: x(:)
  !>               call z%source(x,n) 
  !>             \endcode
  !> @param[in]  self Vector *z*, caller
  !> @param[out] x    Output array of copies (allocatable, n vectors on exit)
  !> @param[in]  n    Number of identical copies required
  !-----------------------------------------------------------------------------
  subroutine s_vector_space__source_array( self, x, n )
    class(sll_c_vector_space)             , intent(in   ) :: self
    class(sll_c_vector_space), allocatable, intent(  out) :: x(:)
    integer                               , intent(in   ) :: n

    allocate( x(n), source=self )

  end subroutine s_vector_space__source_array

end module sll_m_vector_space_base
