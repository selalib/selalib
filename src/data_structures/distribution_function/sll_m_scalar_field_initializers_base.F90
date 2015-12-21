! Essentially, the 'initializer' has a single objective which is to provide
! a uniform interface to a function of a given number of arguments (i.e.: 
! 4, for the 4D case). The object itself encapsulates any parameters that the 
! function might need. There is a choice to make. Either:
!
! 1. The object is limited to encapsulating the function parameters. Thus
!    the caller needs to define the loops and thus to find out the index
!    limits to carry out the initialization of the data array. Or,
!
! 2. the object also encapsulates information about the extent of the data
!    array to initialize. In the case of sequential initializers, this usually
!    means that the object must contain a reference to the underlying mesh.
!    In parallel cases, a reference to a layout object.
!
! The choice here has been the second option, to reduce the amount of work
! that the caller has to do. Hence the 'initializer' serves as the uniform
! interface to a simple function and has knowledge of the domain of 
! definition (mesh or layout).
!
! As mentioned above, whether sequential or parallel, initializers have
! the same interface. The only difference among them being the dimensionality
! of the data to be initialized. It is for this reason that we include all
! the abstract types in the same file (this one).

!>@defgroup fields  sll_m_scalar_field_initializers_base
module sll_m_scalar_field_initializers_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_p_cell_centered_field, &
    sll_p_node_centered_field, &
    sll_c_scalar_field_2d_initializer_base, &
    sll_c_scalar_field_4d_initializer_base, &
    sll_c_scalar_field_6d_initializer_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer, parameter :: sll_p_node_centered_field = 0, sll_p_cell_centered_field = 1


  ! **************************************************************************
  !
  !                              2D cases
  !
  ! **************************************************************************

  type, abstract :: sll_c_scalar_field_2d_initializer_base
     sll_int32   :: data_position
   contains
     procedure(scalar_field_2d_initializer), deferred, pass :: f_of_x1x2
  end type sll_c_scalar_field_2d_initializer_base

  abstract interface
     subroutine scalar_field_2d_initializer( init_obj, data_out )
       use sll_m_working_precision
       import sll_c_scalar_field_2d_initializer_base
       class(sll_c_scalar_field_2d_initializer_base), intent(inout) :: init_obj
       sll_real64, dimension(:,:), intent(out)                :: data_out
     end subroutine scalar_field_2d_initializer
  end interface

  ! **************************************************************************
  !
  !                              4D cases
  !
  ! **************************************************************************

  type, abstract :: sll_c_scalar_field_4d_initializer_base
     sll_int32 :: data_position
   contains
     procedure(scalar_field_4d_initializer), deferred, pass :: f_of_4args
  end type sll_c_scalar_field_4d_initializer_base


  abstract interface
     subroutine scalar_field_4d_initializer( init_obj, data_out )
       use sll_m_working_precision
       import sll_c_scalar_field_4d_initializer_base
       class(sll_c_scalar_field_4d_initializer_base), intent(inout) :: init_obj
       sll_real64, dimension(:,:,:,:), intent(out)            :: data_out
     end subroutine scalar_field_4d_initializer
  end interface

  ! **************************************************************************
  !
  !                              6D cases
  !
  ! **************************************************************************

  type, abstract :: sll_c_scalar_field_6d_initializer_base
     sll_int32 :: data_position
   contains
     procedure(scalar_field_6d_initializer), deferred, pass :: f_of_6args
  end type sll_c_scalar_field_6d_initializer_base


  abstract interface
     subroutine scalar_field_6d_initializer( init_obj, data_out )
       use sll_m_working_precision
       import sll_c_scalar_field_6d_initializer_base
       class(sll_c_scalar_field_6d_initializer_base), intent(inout) :: init_obj
       sll_real64, dimension(:,:,:,:,:,:), intent(out)        :: data_out
     end subroutine scalar_field_6d_initializer
  end interface

end module sll_m_scalar_field_initializers_base
