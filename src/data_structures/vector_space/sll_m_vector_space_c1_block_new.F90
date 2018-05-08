!> @ingroup vector_space
!> @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!> @authors Edoardo Zoni   - <edoardo.zoni@ipp.mpg.de>
!> @brief   Vector space for wrapping 2D Fortran real arrays
!
module sll_m_vector_space_c1_block_new
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_vector_space_base, only: sll_c_vector_space

  use sll_m_vector_space_real_array_1d, only: sll_t_vector_space_real_array_1d

  use sll_m_vector_space_real_array_2d, only: sll_t_vector_space_real_array_2d

  implicit none

  public :: sll_t_vector_space_c1_block_new

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  type, extends(sll_c_vector_space) :: sll_t_vector_space_c1_block_new

    integer :: n1
    integer :: n2
    integer :: p1
    integer :: p2

    type(sll_t_vector_space_real_array_1d) :: vd ! dense component
    type(sll_t_vector_space_real_array_2d) :: vs ! stencil component

  contains

    procedure :: init => s_vector_space_c1_block_new__init

    !> @name Basic operations (overloading the abstract methods)
    procedure :: copy => s_vector_space_c1_block_new__copy
    procedure :: incr => s_vector_space_c1_block_new__incr
    procedure :: scal => s_vector_space_c1_block_new__scal

    !> @name Additional operations
    procedure :: add       => s_vector_space_c1_block_new__add
    procedure :: mult      => s_vector_space_c1_block_new__mult
    procedure :: mult_add  => s_vector_space_c1_block_new__mult_add
    procedure :: incr_mult => s_vector_space_c1_block_new__incr_mult
    procedure :: lcmb      => s_vector_space_c1_block_new__lcmb
    procedure :: incr_lcmb => s_vector_space_c1_block_new__incr_lcmb
 
    !> @name Optional subroutines and functions
    procedure :: norm  => f_vector_space_c1_block_new__norm
    procedure :: inner => f_vector_space_c1_block_new__inner
    !> @}

  end type sll_t_vector_space_c1_block_new

  ! Error messages
  character(len=*), parameter :: wrong_type_x = "x not of type 'sll_t_vector_space_c1_block_new'"
  character(len=*), parameter :: wrong_type_y = "y not of type 'sll_t_vector_space_c1_block_new'"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_vector_space_c1_block_new__init( self, n1, n2, p1, p2 )
    class(sll_t_vector_space_c1_block_new), intent(inout) :: self
    integer                               , intent(in   ) :: n1
    integer                               , intent(in   ) :: n2
    integer                               , intent(in   ) :: p1
    integer                               , intent(in   ) :: p2

    self % n1 = n1
    self % n2 = n2
    self % p1 = p1
    self % p2 = p2

  end subroutine s_vector_space_c1_block_new__init

  !-----------------------------------------------------------------------------
  subroutine s_vector_space_c1_block_new__copy( self, x )
    class(sll_t_vector_space_c1_block_new), intent(inout) :: self
    class(sll_c_vector_space)             , intent(in   ) :: x

    character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block_new % copy"

    select type ( x )

    class is ( sll_t_vector_space_c1_block_new )

      call self % vd % copy( x % vd )
      call self % vs % copy( x % vs )

    class default
      SLL_ERROR( this_sub_name, wrong_type_x )

    end select

  end subroutine s_vector_space_c1_block_new__copy

  !-----------------------------------------------------------------------------
  subroutine s_vector_space_c1_block_new__incr( self, x )
    class(sll_t_vector_space_c1_block_new), intent(inout) :: self
    class(sll_c_vector_space)             , intent(in   ) :: x

    character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block_new % incr"

    select type ( x )

    class is ( sll_t_vector_space_c1_block_new )

      call self % vd % incr( x % vd )
      call self % vs % incr( x % vs )

    class default
      SLL_ERROR( this_sub_name, wrong_type_x )

    end select

  end subroutine s_vector_space_c1_block_new__incr

  !-----------------------------------------------------------------------------
  subroutine s_vector_space_c1_block_new__scal( self, a )
    class(sll_t_vector_space_c1_block_new), intent(inout) :: self
    real(wp)                              , intent(in   ) :: a

    call self % vd % scal( a )
    call self % vs % scal( a )

  end subroutine s_vector_space_c1_block_new__scal

  !-----------------------------------------------------------------------------
  subroutine s_vector_space_c1_block_new__add( self, x, y )
    class(sll_t_vector_space_c1_block_new), intent(inout) :: self
    class(sll_c_vector_space)             , intent(in   ) :: x
    class(sll_c_vector_space)             , intent(in   ) :: y

    character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block_new % add"

    select type( x )

    class is( sll_t_vector_space_c1_block_new )

      select type( y )

      class is( sll_t_vector_space_c1_block_new )

        call self % vd % add( x % vd, y % vd )
        call self % vs % add( x % vs, y % vs )

      class default

        SLL_ERROR( this_sub_name, wrong_type_y )

      end select

    class default

      SLL_ERROR( this_sub_name, wrong_type_x )

    end select

  end subroutine s_vector_space_c1_block_new__add

  !-----------------------------------------------------------------------------
  subroutine s_vector_space_c1_block_new__mult( self, a, x )
    class(sll_t_vector_space_c1_block_new), intent(inout) :: self
    real(wp)                              , intent(in   ) :: a
    class(sll_c_vector_space)             , intent(in   ) :: x

    character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block_new % mult"

    select type( x )

    class is( sll_t_vector_space_c1_block_new )

      call self % vd % mult( a, x % vd )
      call self % vs % mult( a, x % vs )

    class default

      SLL_ERROR( this_sub_name, wrong_type_x )

    end select

  end subroutine s_vector_space_c1_block_new__mult

  !-----------------------------------------------------------------------------
  subroutine s_vector_space_c1_block_new__mult_add( self, a, x, y )
    class(sll_t_vector_space_c1_block_new), intent(inout) :: self
    real(wp)                              , intent(in   ) :: a
    class(sll_c_vector_space)             , intent(in   ) :: x
    class(sll_c_vector_space)             , intent(in   ) :: y

    character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block_new % mult_add"

    select type( x )

    class is( sll_t_vector_space_c1_block_new )

      select type( y )

      class is( sll_t_vector_space_c1_block_new )

        call self % vd % mult_add( a, x % vd, y % vd )
        call self % vs % mult_add( a, x % vs, y % vs )

      class default

        SLL_ERROR( this_sub_name, wrong_type_y )

      end select

    class default

      SLL_ERROR( this_sub_name, wrong_type_x )

    end select

  end subroutine s_vector_space_c1_block_new__mult_add

  !-----------------------------------------------------------------------------
  subroutine s_vector_space_c1_block_new__incr_mult( self, a, x )
    class(sll_t_vector_space_c1_block_new), intent(inout) :: self
    real(wp)                              , intent(in   ) :: a
    class(sll_c_vector_space)             , intent(in   ) :: x

    character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block_new % incr_mult"

    select type( x )

    class is( sll_t_vector_space_c1_block_new )

      call self % vd % incr_mult( a, x % vd )
      call self % vs % incr_mult( a, x % vs )

    class default

      SLL_ERROR( this_sub_name, wrong_type_x )

    end select

  end subroutine s_vector_space_c1_block_new__incr_mult

  !----------------------------------------------------------------------------
  subroutine s_vector_space_c1_block_new__lcmb( self, a, x )
    class(sll_t_vector_space_c1_block_new), intent(inout) :: self
    real(wp)                              , intent(in   ) :: a(:)
    class(sll_c_vector_space)             , intent(in   ) :: x(:)

    character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block_new % lcmb"

    select type( x )

    class is( sll_t_vector_space_c1_block_new )

      call self % vd % lcmb( a, x % vd )
      call self % vs % lcmb( a, x % vs )

    class default

      SLL_ERROR( this_sub_name, wrong_type_x )

    end select

  end subroutine s_vector_space_c1_block_new__lcmb

  !-----------------------------------------------------------------------------
  subroutine s_vector_space_c1_block_new__incr_lcmb( self, a, x )
    class(sll_t_vector_space_c1_block_new), intent(inout) :: self
    real(wp)                              , intent(in   ) :: a(:)
    class(sll_c_vector_space)             , intent(in   ) :: x(:)

    character(len=*), parameter :: this_sub_name = "sll_t_vector_space_c1_block_new % incr_lcmb"

    select type( x )

    class is( sll_t_vector_space_c1_block_new )

      call self % vd % incr_lcmb( a, x % vd )
      call self % vs % incr_lcmb( a, x % vs )

    class default

      SLL_ERROR( this_sub_name, wrong_type_x )

    end select

  end subroutine s_vector_space_c1_block_new__incr_lcmb

  !-----------------------------------------------------------------------------
  function f_vector_space_c1_block_new__norm( self ) result( res )
    class(sll_t_vector_space_c1_block_new), intent(in) :: self
    real(wp) :: res

    res = sqrt( self % inner( self ) )

  end function f_vector_space_c1_block_new__norm

  !-----------------------------------------------------------------------------
  function f_vector_space_c1_block_new__inner( self, x ) result( res )
    class(sll_t_vector_space_c1_block_new), intent(in) :: self
    class(sll_c_vector_space)             , intent(in) :: x
    real(wp) :: res

    character(len=*), parameter :: this_fun_name = "sll_t_vector_space_c1_block_new % inner"

    select type( x )

    class is( sll_t_vector_space_c1_block_new )

      associate( n1 => self % n1, n2 => self % n2 )
        res = self % vd % inner( x % vd ) + &
              sum( self % vs % array(1:n1,1:n2) * x % vs % array(1:n1,1:n2) )
      end associate

    class default

      SLL_ERROR( this_fun_name, wrong_type_x )

    end select

  end function f_vector_space_c1_block_new__inner

end module sll_m_vector_space_c1_block_new
