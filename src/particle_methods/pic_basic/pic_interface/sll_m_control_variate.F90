module sll_m_control_variate

#include "sll_working_precision.h"
#include "sll_memory.h"

  implicit none
  private

  public :: sll_f_control_variate, sll_new_control_variate

  type, public :: sll_t_control_variate

     sll_real64, pointer :: control_variate_parameters(:) => null()

     procedure(sll_f_control_variate), pointer :: control_variate => null()

   contains

     procedure :: update_df_weight !< function defining the control variate

  end type sll_t_control_variate


  abstract interface
     !> 1d real function
     function sll_f_control_variate(this, xi, vi, time)
       use sll_m_working_precision 
       import sll_t_control_variate
       class(sll_t_control_variate) :: this
       sll_real64, intent( in ) :: xi(3) !< particle position
       sll_real64, intent( in ) :: vi(3) !< particle velocity
       sll_real64, intent( in ) :: time  !< current time
       sll_real64               :: sll_f_control_variate

     end function sll_f_control_variate
  end interface


  contains
    
    function update_df_weight(this, xi, vi, time, weight_ff, g0) result(weight_df)
      class(sll_t_control_variate) :: this
      sll_real64, intent( in ) :: xi(3) !< particle position
      sll_real64, intent( in ) :: vi(3) !< particle velocity
      sll_real64, intent( in ) :: time  !< current time
      sll_real64, intent( in ) :: weight_ff !< particle weight for full f
      sll_real64, intent( in ) :: g0 !< initial sampling distribution at particle coordinates at time zero
      sll_real64 :: weight_df !< particle weight for delta f

      weight_df = weight_ff - this%control_variate(xi, vi, time)/g0

    end function update_df_weight


    function sll_new_control_variate(control_function, parameters) result(control_variate)
      procedure(sll_f_control_variate) :: control_function
      class(sll_t_control_variate), pointer :: control_variate
      sll_real64, pointer, intent(in) :: parameters(:)

      control_variate%control_variate => control_function
      
      control_variate%control_variate_parameters => parameters


    end function sll_new_control_variate


end module sll_m_control_variate
