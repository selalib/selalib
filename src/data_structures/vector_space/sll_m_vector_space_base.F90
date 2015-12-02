!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE:   sll_m_vector_space_base
!
! DESCRIPTION:
!> @ingroup vector_space
!> @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!> @authors Marco Restelli - <marco.restelli@gmail.com>
!> @brief   Abstract type implementing a generic vector space.
!> @todo    Add detailed description
!------------------------------------------------------------------------------
module sll_m_vector_space_base

#include "sll_working_precision.h"
#include "sll_errors.h"

  implicit none
!  public :: sll_vector_space_base
  private

  !============================================================================
  ! Abstract type definition
  !============================================================================

  !> @brief   Abstract base class for all vector spaces
  type, public, abstract :: sll_vector_space_base

  contains
    !> @name Basic operations (abstract methods)
    !> Any non-abstract extended type MUST provide an implementation for these.
    procedure( i_copy ), deferred :: copy       !< z = x
    procedure( i_incr ), deferred :: incr       !< z+= x
    procedure( i_scal ), deferred :: scal       !< z*= a

    !> @name Additional operations
    !> We provide a standard implementation that uses
    !> the operations above and creates temporary objects where needed.
    !> The user should provide a more efficient implementation.
    procedure          :: add       => add__base       !< z = x + y
    procedure          :: mult      => mult__base      !< z = a * x
    procedure          :: mult_add  => mult_add__base  !< z = a * x + y
    procedure          :: incr_mult => incr_mult__base !< z+= a * x
    procedure          :: lcmb      => lcmb__base      !< z = sum_i ( a_i*x_i )
    procedure          :: incr_lcmb => incr_lcmb__base !< z+= sum_i ( a_i*x_i )

    !> @name Optional subroutines and functions
    !> Provide a *norm* to have a normed vector space.
    !> Provide an *inner product* to have an inner product space.
    !> Provide a *show* function for debugging purposes.
    procedure          :: norm      => norm__base     !< norm(z)
    procedure          :: inner     => inner__base    !< inner(z,x)
    procedure          :: show      => show__base

    !> @name Constructors & destructors
    generic            :: source    => source_scalar, source_array
    procedure          :: delete    => delete__base
    !> @}

    !> @name Copy constructors (one vector or an array of vectors)
    procedure, private :: initialize_copy => initialize_copy__base
    procedure, private :: source_scalar   => source_scalar__base
    procedure, private :: source_array    => source_array__base
    !> @}

  end type sll_vector_space_base

  !============================================================================
  ! Abstract interfaces: basic operations
  !============================================================================

  !----------------------------------------------------------------------------
  !> @brief        z = x
  !> @details      Copy contents of input array *x* into array *z* ('self').
  !>               Usage: \code z%copy(x) \endcode
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    x    Input vector
  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_copy( self, x )
    import sll_vector_space_base
    class( sll_vector_space_base ), intent( inout ) :: self
    class( sll_vector_space_base ), intent( in    ) :: x
   end subroutine i_copy
  end interface

  !----------------------------------------------------------------------------
  !> @brief        z += x
  !> @details      Increment the vector *z* ('self') by the vector *x*.
  !>               Usage: \code z%incr(x) \endcode
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    x    Input vector
  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_incr( self, x )
    import sll_vector_space_base
    class( sll_vector_space_base ), intent( inout ) :: self
    class( sll_vector_space_base ), intent( in    ) :: x
   end subroutine i_incr
  end interface

  !----------------------------------------------------------------------------
  !> @brief        z *= a
  !> @details      Scale the vector *z* ('self') by the real number *a*.
  !>               Usage: \code z%scal(a) \endcode
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Input scalar, double precision real
  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_scal( self, a )
    use sll_m_working_precision
    import sll_vector_space_base
    class( sll_vector_space_base ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a
   end subroutine i_scal
  end interface

  !++++++++++++++++++++++++++++++++++
  private :: add__base  ! for Doxygen
  !++++++++++++++++++++++++++++++++++

  contains

  !============================================================================
  ! Default subroutines: composite operations
  !============================================================================

  !----------------------------------------------------------------------------
  !> @brief        z = x + y
  !> @details      Add two vectors *x* and *y*, and write the result to the
  !>               vector *z* ('self') that invoked this subroutine.
  !>               Usage: \code call z%add(x,y) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    x    Input vector #1
  !> @param[in]    y    Input vector #2
  !----------------------------------------------------------------------------
  subroutine add__base( self, x, y )
    class( sll_vector_space_base ), intent( inout ) :: self
    class( sll_vector_space_base ), intent( in    ) :: x
    class( sll_vector_space_base ), intent( in    ) :: y

    call self%copy( x )
    call self%incr( y )

  end subroutine add__base

  !----------------------------------------------------------------------------
  !> @brief        z = a * x
  !> @details      Multiply the vector *x* by the scalar *a*, and assign the
  !>               result to the vector *z* ('self') that invoked this
  !>               subroutine.
  !>               Usage: \code call z%mult(a,x) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Input scalar
  !> @param[in]    x    Input vector
  !----------------------------------------------------------------------------
  subroutine mult__base( self, a, x )
    class( sll_vector_space_base ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a
    class( sll_vector_space_base ), intent( in    ) :: x

    call self%copy( x )
    call self%scal( a )

  end subroutine mult__base

  !----------------------------------------------------------------------------
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
  !----------------------------------------------------------------------------
  subroutine mult_add__base( self, a, x, y )
    class( sll_vector_space_base ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a
    class( sll_vector_space_base ), intent( in    ) :: x
    class( sll_vector_space_base ), intent( in    ) :: y
    
    call self%mult( a, x )
    call self%incr( y )

  end subroutine mult_add__base

  !----------------------------------------------------------------------------
  !> @brief        z += a * x
  !> @details      Increment the vector *z* ('self') by the result of the
  !>               multiplication between the vector *x* and the scalar *a*.
  !>               Usage: \code call z%incr_mult(a,x) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Input scalar
  !> @param[in]    x    Input vector
  !----------------------------------------------------------------------------
  subroutine incr_mult__base( self, a, x )
    class( sll_vector_space_base ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a
    class( sll_vector_space_base ), intent( in    ) :: x
    
    class( sll_vector_space_base ), allocatable     :: temp

    call self%source( temp )
    call temp%mult( a, x )
    call self%incr( temp )
    deallocate( temp )

  end subroutine incr_mult__base

  !----------------------------------------------------------------------------
  !> @brief        z = \f$ \sum_i ( a_i x_i ) \f$
  !> @details      Compute a linear combination of N vectors \f$ x_i \f$,
  !>               according to N coefficients \f$ a_i \f$, and assign the
  !>               result to the vector *z* that invoked this subroutine.
  !>               Usage: \code call z%lcmb(a,x) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Array of N scalars
  !> @param[in]    x    Array of N vectors
  !----------------------------------------------------------------------------
  subroutine lcmb__base( self, a, x )
    class( sll_vector_space_base ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a(:)
    class( sll_vector_space_base ), intent( in    ) :: x(:)

    class( sll_vector_space_base ), allocatable     :: temp
    integer                                         :: i

    call self%source( temp )
    call self%mult( a(1), x(1) )
    do i = 2,size(a)
      call temp%mult( a(i), x(i) )
      call self%incr( temp )
    end do
    deallocate( temp )
 
  end subroutine lcmb__base

  !----------------------------------------------------------------------------
  !> @brief        z += \f$ \sum_i ( a_i x_i ) \f$
  !> @details      Compute a linear combination of N vectors \f$ x_i \f$,
  !>               according to N coefficients \f$ a_i \f$, and use the result
  !>               to increment the vector *z* that invoked this subroutine.
  !>               Usage: \code call z%incr_lcmb(a,x) \endcode
  !>
  !> @param[inout] self Result vector *z*, caller
  !> @param[in]    a    Array of N scalars
  !> @param[in]    x    Array of N vectors
  !----------------------------------------------------------------------------
  subroutine incr_lcmb__base( self, a, x )
    class( sll_vector_space_base ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a(:)
    class( sll_vector_space_base ), intent( in    ) :: x(:)

    class( sll_vector_space_base ), allocatable     :: temp
    integer                                         :: i

    call self%source( temp )
    do i = 1,size(a)
      call temp%mult( a(i), x(i) )
      call self%incr( temp )
    end do
    deallocate( temp )

  end subroutine incr_lcmb__base

  !============================================================================
  ! Empty functions/subroutines: optional operations and debugging
  !============================================================================

  !----------------------------------------------------------------------------
  !> @brief     Norm of vector: ||z||
  !> @details   Compute the norm of the vector *z* ('self') that is invoking
  !>            the function. By default, \c norm(z)=inner(z,z).
  !>            Usage: \code r = z%norm() \endcode
  !>
  !> @param[in] self Vector *z*, caller
  !> @returns        Scalar value, of type *sll_real64*
  !----------------------------------------------------------------------------
  function norm__base( self ) result( res )
    class( sll_vector_space_base ), intent( in ) :: self
    sll_real64                                   :: res

    res = self%inner( self )

  end function norm__base

  !----------------------------------------------------------------------------
  !> @brief     Inner product: \<z,x\>
  !> @details   Compute the inner product between the vector *z* ('self') and
  !>            another vector *x* of the same type.
  !>            Usage: \code r = z%inner(x) \endcode
  !>
  !> @param[in] self Vector *z*, caller
  !> @param[in] x    Vector
  !> @returns        Scalar value, of type *sll_real64*
  !----------------------------------------------------------------------------
  function inner__base( self, x ) result( res )
    class( sll_vector_space_base ), intent( in ) :: self
    class( sll_vector_space_base ), intent( in ) :: x
    sll_real64                                   :: res

    SLL_ERROR( "inner", "Function not implemented." )
#ifdef DEBUG
    print*, storage_size(self), storage_size(x)
#endif
    res = 0.0_f64

  end function inner__base

  !----------------------------------------------------------------------------
  !> @brief     Show something, for debug
  !> @param[in] self Vector *z*, caller
  !----------------------------------------------------------------------------
  subroutine show__base( self )
    class( sll_vector_space_base ), intent( in ) :: self

    SLL_WARNING( "show", "Overload this subroutine if you need it." )
    print*, storage_size(self)

  end subroutine show__base

  !----------------------------------------------------------------------------
  !> @brief        Destructor, release memory
  !> @param[inout] self Vector *z*, caller
  !----------------------------------------------------------------------------
  subroutine delete__base( self )
    class( sll_vector_space_base ), intent( inout ) :: self

    SLL_ERROR( "delete", "Subroutine not implemented." )
#ifdef DEBUG
    print*, storage_size(self)
#endif

  end subroutine delete__base

  !============================================================================
  ! Copy constructors (all private)
  !============================================================================
  ! The 'source' subroutines allocate storage for a new variable 'x' and its
  ! 'allocatable' fields, then copy the values of each field from 'self'.
  ! Moreover, pointers in 'x' will point to the same targets as the
  ! pointers in 'self'.
  !
  ! The 'initialize_copy' subroutine (empty for now) should perform additional
  ! operations on the object before making it available to the user.
  !
  ! If not all the 'allocatable' fields should be allocated, please overload
  ! both the 'source' functions.
  !
  ! If pointers should point to some separate storage (deep copy), please
  ! overload the 'initialize_copy' subroutine.

  !----------------------------------------------------------------------------
  !> @brief        Initialize a newly created vector copy *x*
  !> @details      Perform additional operations on a new copy *x* of some
  !>               vector *z*.  For example, allocate new storage for pointers
  !>               and then copy the data.
  !>               Usage: \code call x%initialize_copy() \endcode
  !>
  !> @param[inout] self Vector *x*, caller
  !----------------------------------------------------------------------------
  subroutine initialize_copy__base( self )
    class( sll_vector_space_base ), intent( inout ) :: self
    ! Do nothing for now
#ifdef DEBUG
    print*, storage_size(self)
#endif

  end subroutine initialize_copy__base

  !----------------------------------------------------------------------------
  !> @brief      Copy constructor: create one copy of vector *z*
  !> @details    Allocate vector *x* and copy contents of *z* into it. Usage:
  !>             \code
  !>               class( sll_vector_space_base ), allocatable :: x
  !>               call z%source(x) 
  !>             \endcode
  !> @param[in]  self Vector *z*, caller
  !> @param[out] x    Output copy (allocatable)
  !----------------------------------------------------------------------------
  subroutine source_scalar__base( self, x )
    class( sll_vector_space_base ),              intent( in    ) :: self
    class( sll_vector_space_base ), allocatable, intent(   out ) :: x

    allocate( x, source=self )
    call x%initialize_copy()

  end subroutine source_scalar__base

  !----------------------------------------------------------------------------
  !> @brief      Copy constructor: create *n* copies of vector *z*
  !> @details    Allocate array *x* of *n* vectors, and copy contents of *z*
  !>             into each element *x(i)*. Usage:
  !>             \code
  !>               class( sll_vector_space_base ), allocatable :: x(:)
  !>               call z%source(x,n) 
  !>             \endcode
  !> @param[in]  self Vector *z*, caller
  !> @param[out] x    Output array of copies (allocatable, n vectors on exit)
  !> @param[in]  n    Number of identical copies required
  !----------------------------------------------------------------------------
  subroutine source_array__base( self, x, n )
    class( sll_vector_space_base ),              intent( in    ) :: self
    class( sll_vector_space_base ), allocatable, intent(   out ) :: x(:)
    sll_int32                                  , intent( in    ) :: n

    sll_int32 :: i

    allocate( x(n), source=self )
    do i = 1,n
      call x(i)%initialize_copy()
    end do

  end subroutine source_array__base

  !============================================================================

end module sll_m_vector_space_base
