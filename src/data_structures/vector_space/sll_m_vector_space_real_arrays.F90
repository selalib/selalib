!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE:   sll_m_vector_space_real_arrays
!
! DESCRIPTION:
!> @ingroup vector_space
!> @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!> @authors Marco Restelli - <marco.restelli@gmail.com>
!> @brief   Vector spaces for wrapping 1D/2D/3D Fortran arrays.
!> @todo    Add detailed description
!------------------------------------------------------------------------------
module sll_m_vector_space_real_arrays

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_working_precision.h"

  use sll_m_vector_space_base, only: &
    sll_vector_space_base

  implicit none

  public :: &
    sll_vector_space_real_1d, &
    sll_vector_space_real_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !============================================================================
  ! Generic interface for real arrays - assumes contiguous array storage
  !============================================================================
  !> @brief   Generic interface for real arrays
  !> @details Abstract class providing the efficient implementation of all
  !>          basic operations on a vector space that wraps a single Fortran 
  !>          multidimensional array of double precision real numbers.
  !>          This abstract type is extended to provide specific interfaces to
  !>          1D/2D/3D arrays.
  !> @note    Internally, the data is accessed through a 1D array pointer.
  !>          In order to access multidimentional arrays, we must assume
  !>          contiguous array storage.
  !> @note    An array allocated by an *allocate* statement is always
  !>          contiguous, hence we strongly suggest to use allocatable arrays.
  !============================================================================
  type, public, abstract, &
                      extends( sll_vector_space_base ) :: sll_vector_space_real

    sll_real64, private, pointer, contiguous :: d(:)        => null()
    logical   , private                      :: owns_memory = .false.

  contains
    !> @name Basic operations (overloading the abstract methods)
    procedure      :: copy      => copy__real
    procedure      :: incr      => incr__real
    procedure      :: scal      => scal__real

    !> @name Additional operations
    procedure      :: add       => add__real
    procedure      :: mult      => mult__real
    procedure      :: mult_add  => mult_add__real
    procedure      :: incr_mult => incr_mult__real
    procedure      :: lcmb      => lcmb__real
    procedure      :: incr_lcmb => incr_lcmb__real
 
    !> @name Optional subroutines and functions
    procedure      :: norm      => norm__real
    procedure      :: inner     => inner__real
    !> @}

  end type sll_vector_space_real

  !============================================================================
  ! Error messages
  !============================================================================
  character( len=* ), parameter :: WRONG_TYPE_X = &
    "x not of class( sll_vector_space_real )"

  character( len=* ), parameter :: WRONG_TYPE_Y = &
    "y not of class( sll_vector_space_real )"

  !============================================================================
  ! 1D array
  !---------
  !> Simple vector space for 1D Fortran arrays (double precision float).
  !============================================================================
  type, extends( sll_vector_space_real ) :: sll_vector_space_real_1d

    sll_real64, pointer, contiguous :: array(:) => null()      !< 1D real array

  contains
    procedure          :: attach          => attach__real_1d !< Attach 2D array
    procedure          :: delete          => delete__real_1d !< Erase contents
    procedure, private :: initialize_copy => initialize_copy__real_1d

  end type sll_vector_space_real_1d

  !============================================================================
  ! 2D array
  !---------
  !> Simple vector space for 2D Fortran arrays (double precision float).
  !============================================================================
  type, extends( sll_vector_space_real ) :: sll_vector_space_real_2d

    sll_real64, pointer, contiguous :: array(:,:) => null()    !< 2D real array

  contains
    procedure          :: attach          => attach__real_2d !< Attach 2D array
    procedure          :: delete          => delete__real_2d !< Erase contents
    procedure, private :: initialize_copy => initialize_copy__real_2d

  end type sll_vector_space_real_2d

  !============================================================================
  ! 3D array
  !---------
  !> Simple vector space for 3D Fortran arrays (double precision float).
  !============================================================================
  type, extends( sll_vector_space_real ) :: sll_vector_space_real_3d

    sll_real64, pointer, contiguous :: array(:,:,:) => null()  !< 3D real array

  contains
    procedure          :: attach          => attach__real_3d !< Attach 3D array
    procedure          :: delete          => delete__real_3d !< Erase contents
    procedure, private :: initialize_copy => initialize_copy__real_3d

  end type sll_vector_space_real_3d

  !============================================================================
contains

  !============================================================================
  ! Generic methods
  !============================================================================

  !----------------------------------------------------------------------------
  subroutine copy__real( self, x )
    class( sll_vector_space_real ), intent( inout ) :: self
    class( sll_vector_space_base ), intent( in    ) :: x
    select type( x ); class is( sll_vector_space_real )
 
    self%d = x%d
 
    class default; SLL_ERROR( "copy__real", WRONG_TYPE_X ); end select
  end subroutine copy__real

  !----------------------------------------------------------------------------
  subroutine incr__real( self, x )
    class( sll_vector_space_real ), intent( inout ) :: self
    class( sll_vector_space_base ), intent( in    ) :: x
    select type( x ); class is( sll_vector_space_real )

    self%d = self%d + x%d

    class default; SLL_ERROR( "incr__real", WRONG_TYPE_X ); end select
  end subroutine incr__real

  !----------------------------------------------------------------------------
  subroutine scal__real( self, a )
    class( sll_vector_space_real ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a

    self%d = self%d * a

  end subroutine scal__real

  !----------------------------------------------------------------------------
  subroutine add__real( self, x, y )
    class( sll_vector_space_real ), intent( inout ) :: self
    class( sll_vector_space_base ), intent( in    ) :: x
    class( sll_vector_space_base ), intent( in    ) :: y
    select type( x ); class is( sll_vector_space_real )
    select type( y ); class is( sll_vector_space_real )

    self%d = x%d + y%d

    class default; SLL_ERROR( "add__real", WRONG_TYPE_Y); end select
    class default; SLL_ERROR( "add__real", WRONG_TYPE_X); end select
  end subroutine add__real

  !----------------------------------------------------------------------------
  subroutine mult__real( self, a, x )
    class( sll_vector_space_real ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a
    class( sll_vector_space_base ), intent( in    ) :: x
    select type( x ); class is( sll_vector_space_real )

    self%d = a * x%d

    class default; SLL_ERROR( "mult__real", WRONG_TYPE_X ); end select
  end subroutine mult__real

  !----------------------------------------------------------------------------
  subroutine mult_add__real( self, a, x, y )
    class( sll_vector_space_real ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a
    class( sll_vector_space_base ), intent( in    ) :: x
    class( sll_vector_space_base ), intent( in    ) :: y
    select type( x ); class is( sll_vector_space_real )
    select type( y ); class is( sll_vector_space_real )
    
    self%d = a * x%d + y%d

    class default; SLL_ERROR( "mult_add__real", WRONG_TYPE_Y); end select
    class default; SLL_ERROR( "mult_add__real", WRONG_TYPE_X); end select
  end subroutine mult_add__real

  !----------------------------------------------------------------------------
  subroutine incr_mult__real( self, a, x )
    class( sll_vector_space_real ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a
    class( sll_vector_space_base ), intent( in    ) :: x
    select type( x ); class is( sll_vector_space_real )
    
    self%d = self%d + a * x%d

    class default; SLL_ERROR( "incr_mult__real", WRONG_TYPE_X); end select
  end subroutine incr_mult__real

  !----------------------------------------------------------------------------
  subroutine lcmb__real( self, a, x )
    class( sll_vector_space_real ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a(:)
    class( sll_vector_space_base ), intent( in    ) :: x(:)

    integer    :: i
!    integer    :: j
!    sll_real64 :: temp

    select type( x ); class is( sll_vector_space_real )

    ! Option 1: cycle over each array position first
!    do j = 1, size( self%d )
!      temp = 0.0_f64
!      do i = 1, size( a )
!        temp = temp + a(i) * x(i)%d(j)
!      end do
!      self%d(j) = temp
!    end do

    ! Option 2: construct linear combination incrementally
    self%d = a(1) * x(1)%d
    do i = 2, size( a )
      self%d = self%d + a(i) * x(i)%d
    end do

    class default; SLL_ERROR( "lcmb__real", WRONG_TYPE_X); end select
  end subroutine lcmb__real

  !----------------------------------------------------------------------------
  subroutine incr_lcmb__real( self, a, x )
    class( sll_vector_space_real ), intent( inout ) :: self
    sll_real64                    , intent( in    ) :: a(:)
    class( sll_vector_space_base ), intent( in    ) :: x(:)

    integer  :: i

    select type( x ); class is( sll_vector_space_real )

    ! Option 2: construct linear combination incrementally
    do i = 1, size( a )
      self%d = self%d + a(i) * x(i)%d
    end do

    class default; SLL_ERROR( "incr_lcmb__real", WRONG_TYPE_X); end select
  end subroutine incr_lcmb__real

  !----------------------------------------------------------------------------
  function norm__real( self ) result( res )
    class( sll_vector_space_real ), intent( in ) :: self
    sll_real64                                   :: res

    res = dot_product( self%d, self%d )

  end function norm__real

  !----------------------------------------------------------------------------
  function inner__real( self, x ) result( res )
    class( sll_vector_space_real ), intent( in ) :: self
    class( sll_vector_space_base ), intent( in ) :: x
    sll_real64                                   :: res
    select type( x ); class is( sll_vector_space_real )

    res = dot_product( self%d, x%d )

    class default; SLL_ERROR( "inner__real", WRONG_TYPE_X); end select
  end function inner__real

  !============================================================================
  ! 1D array
  !============================================================================

  !----------------------------------------------------------------------------
  !> @brief        \c "self%array => x"
  !> @details      Attach the wrapper *z* (a 1D vector object) to an existing
  !>               1D Fortran array *x* (which must be pointer or non-pointer
  !>               target). Usage:
  !>               \code call z%attach( x ) \endcode
  !>
  !> @param[inout] self 1D vector *z*, caller
  !> @param[inout] x    1D Fortran array
  !----------------------------------------------------------------------------
  subroutine attach__real_1d( self, x )
    class( sll_vector_space_real_1d ),         intent( inout ) :: self
    sll_real64                       , target, intent( inout ) :: x(:)
    ! Safer option: Automatic pointer targetting, not yet supported by ifort:
    ! sll_real64, pointer, intent( in ) :: x(:)

    self%array => x
    self%d     => self%array
    self%owns_memory = .false.

  end subroutine attach__real_1d

  !----------------------------------------------------------------------------
  subroutine initialize_copy__real_1d( self )
    class( sll_vector_space_real_1d ), intent( inout ) :: self

    sll_real64, pointer, contiguous :: p(:)
    sll_int32                       :: sz, lb(1), ub(1)

    p  => self%array
    sz =  size  ( p )
    lb =  lbound( p )
    ub =  ubound( p )

    allocate( self%array( lb(1):ub(1) ), source=p )

    self%d(1:sz) => self%array
    self%owns_memory = .true.

  end subroutine initialize_copy__real_1d

  !----------------------------------------------------------------------------
  !> @brief        Erase vector content
  !> @details      If wrapper, detach real vector space from 1D array target.
  !>               If copy, deallocate local memory.
  !>               \code call z%delete() \endcode
  !>
  !> @param[inout] self 1D vector *z*, caller
  !----------------------------------------------------------------------------
  subroutine delete__real_1d( self )
    class( sll_vector_space_real_1d ), intent( inout ) :: self

    if( self%owns_memory ) then
      deallocate( self%d )
    endif

    self%d     => null()
    self%array => null()

  end subroutine delete__real_1d

  !============================================================================
  ! 2D array
  !============================================================================

  !----------------------------------------------------------------------------
  !> @brief        \c "self%array => x"
  !> @details      Attach the wrapper *z* (a 2D vector object) to an existing
  !>               2D Fortran array *x* (which must be pointer or non-pointer
  !>               target, and it must be contiguous). Usage:
  !>               \code call z%attach( x ) \endcode
  !>
  !> @param[inout] self 2D vector *z*, caller
  !> @param[inout] x    2D Fortran array
  !----------------------------------------------------------------------------
  subroutine attach__real_2d( self, x )
    class( sll_vector_space_real_2d ),         intent( inout ) :: self
    sll_real64,            contiguous, target, intent( inout ) :: x(:,:)
    ! Safer option: Automatic pointer targetting, not yet supported by ifort:
    ! sll_real64, contiguous, pointer, intent( in ) :: x(:,:)

    self%array        => x
    self%d(1:size(x)) => x

  end subroutine attach__real_2d

  !----------------------------------------------------------------------------
  subroutine initialize_copy__real_2d( self )
    class( sll_vector_space_real_2d ), intent( inout ) :: self

    sll_real64, pointer, contiguous :: p(:,:)
    sll_int32                       :: sz, lb(2), ub(2)

    p  => self%array
    sz =  size  ( p )
    lb =  lbound( p )
    ub =  ubound( p )

    allocate( self%array( lb(1):ub(1), lb(2):ub(2) ), source=p )

    p => self%array   !
    self%d(1:sz) => p !self%array
    self%owns_memory = .true.

  end subroutine initialize_copy__real_2d

  !----------------------------------------------------------------------------
  !> @brief        Erase vector content
  !> @details      If wrapper, detach real vector space from 2D array target.
  !>               If copy, deallocate local memory.
  !>               \code call z%delete() \endcode
  !>
  !> @param[inout] self 2D vector *z*, caller
  !----------------------------------------------------------------------------
  subroutine delete__real_2d( self )
    class( sll_vector_space_real_2d ), intent( inout ) :: self

    if( self%owns_memory ) then
      deallocate( self%d )
    endif

    self%d     => null()
    self%array => null()

  end subroutine delete__real_2d

  !============================================================================
  ! 3D array
  !============================================================================

  !----------------------------------------------------------------------------
  !> @brief        \c "self%array => x"
  !> @details      Attach the wrapper *z* (a 3D vector object) to an existing
  !>               3D Fortran array *x* (which must be pointer or non-pointer
  !>               target, and it must be contiguous). Usage:
  !>               \code call z%attach( x ) \endcode
  !>
  !> @param[inout] self 3D vector *z*, caller
  !> @param[inout] x    3D Fortran array
  !----------------------------------------------------------------------------
  subroutine attach__real_3d( self, x )
    class( sll_vector_space_real_3d ),         intent( inout ) :: self
    sll_real64,           contiguous , target, intent( inout ) :: x(:,:,:)
    ! Safer option: Automatic pointer targetting, not yet supported by ifort:
    ! sll_real64, contiguous, pointer, intent( in ) :: x(:,:,:)

    self%array        => x
    self%d(1:size(x)) => x

  end subroutine attach__real_3d

  !----------------------------------------------------------------------------
  subroutine initialize_copy__real_3d( self )
    class( sll_vector_space_real_3d ), intent( inout ) :: self

    sll_real64, pointer, contiguous :: p(:,:,:)
    sll_int32                       :: sz, lb(3), ub(3)

    p  => self%array
    sz =  size  ( p )
    lb =  lbound( p )
    ub =  ubound( p )

    allocate( self%array( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ), source=p )

    p => self%array   !
    self%d(1:sz) => p !self%array
    self%owns_memory = .true.

  end subroutine initialize_copy__real_3d

  !----------------------------------------------------------------------------
  !> @brief        Erase vector content
  !> @details      If wrapper, detach real vector space from 3D array target.
  !>               If copy, deallocate local memory.
  !>               \code call z%delete() \endcode
  !>
  !> @param[inout] self 3D vector *z*, caller
  !----------------------------------------------------------------------------
  subroutine delete__real_3d( self )
    class( sll_vector_space_real_3d ), intent( inout ) :: self

    if( self%owns_memory ) then
      deallocate( self%d )
    endif

    self%d     => null()
    self%array => null()

  end subroutine delete__real_3d

  !============================================================================

end module sll_m_vector_space_real_arrays
