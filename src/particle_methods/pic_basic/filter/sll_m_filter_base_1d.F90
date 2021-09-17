!> @ingroup filter
!> @author Benedikt Perse
!> @brief filter base class
!> @details ...

module sll_m_filter_base_1d

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
       sll_c_filter_base_1d

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, abstract :: sll_c_filter_base_1d

     sll_int32 :: iterations = 0
     sll_int32 :: n_dofs = 2
     
   contains

     procedure(initialise), deferred :: init
     procedure(filter_from_field_to_field), deferred :: apply
     procedure(filter_one_field), deferred :: apply_inplace

  end type sll_c_filter_base_1d

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine initialise( self, iterations, n_dofs, mode )
       use sll_m_working_precision
       import sll_c_filter_base_1d
       class( sll_c_filter_base_1d ), intent( out ) :: self
       sll_int32, intent( in ) :: iterations
       sll_int32, intent( in ) :: n_dofs
       sll_int32, optional, intent( in ) :: mode

     end subroutine initialise
  end interface

  abstract interface
     subroutine filter_from_field_to_field( self, field_in, field_out )
       use sll_m_working_precision
       import sll_c_filter_base_1d
       class( sll_c_filter_base_1d ), intent( inout ) :: self
       sll_real64,                            intent(in)  :: field_in(:) !< array for the coefficients of the fields 
       sll_real64,                            intent(out)  :: field_out(:) !< array for the coefficients of the fields
     end subroutine filter_from_field_to_field
  end interface

  abstract interface
     subroutine filter_one_field( self, field )
       use sll_m_working_precision
       import sll_c_filter_base_1d
       class( sll_c_filter_base_1d ), intent( inout ) :: self
       sll_real64,                            intent(inout)  :: field(:) !< array for the coefficients of the fields 
     end subroutine filter_one_field
  end interface


end module sll_m_filter_base_1d
