!> @ingroup pic_weights
!> @author Katharina Kormann, IPP
!> @brief Control variate
!> @details This base class gives an abstract interface to the basic functions for accumulation of charge and current densities as well as the evaluation of a function at particle positions.

module sll_m_control_variate

#include "sll_working_precision.h"
#include "sll_memory.h"

  implicit none
  private

  public :: sll_i_control_variate

  !> Control variate object
  type, public :: sll_t_control_variate

     sll_real64, pointer :: control_variate_parameters(:) => null()

     procedure(sll_i_control_variate), pointer :: control_variate => null()

   contains

     procedure :: update_df_weight !< function defining the control variate
     procedure :: init => init_control_variate !< initialize the type
     procedure :: free => free_control_variate !< finalize

  end type sll_t_control_variate


  abstract interface
     !> 1d real function, abstract interface for function defining the control variate
     function sll_i_control_variate(self, xi, vi, time)
       use sll_m_working_precision 
       import sll_t_control_variate
       class(sll_t_control_variate)       :: self  !< Control variate object      
       sll_real64, optional, intent( in ) :: xi(:) !< particle position
       sll_real64, optional, intent( in ) :: vi(:) !< particle velocity
       sll_real64, optional, intent( in ) :: time  !< current time
       sll_real64                         :: sll_i_control_variate

     end function sll_i_control_variate
  end interface


  contains
    
    !> Update the delta f weights
    function update_df_weight(self, xi, vi, time, weight_ff, g0) result(weight_df)
      class(sll_t_control_variate) :: self !< Control variate object
      sll_real64, intent( in ) :: xi(:) !< particle position
      sll_real64, intent( in ) :: vi(:) !< particle velocity
      sll_real64, intent( in ) :: time  !< current time
      sll_real64, intent( in ) :: weight_ff !< particle weight for full f
      sll_real64, intent( in ) :: g0 !< initial sampling distribution at particle coordinates at time zero
      sll_real64 :: weight_df !< particle weight for delta f

      weight_df = weight_ff - self%control_variate(xi, vi, time)/g0

    end function update_df_weight

    !> Initialization
    subroutine init_control_variate( self, control_function, parameters ) 
      class(sll_t_control_variate), intent(out) :: self !< Control variate object
      procedure(sll_i_control_variate) :: control_function !< Function defining the control variate
      sll_real64, pointer, intent(in) :: parameters(:) !< Parameter values needed in control variate function

      self%control_variate => control_function
      
      self%control_variate_parameters => parameters

      
    end subroutine init_control_variate

    !> Destructor
    subroutine free_control_variate( self ) 
      class(sll_t_control_variate), intent(inout) :: self !< Control variate object
     
      self%control_variate_parameters => null()

    end subroutine free_control_variate



end module sll_m_control_variate
