!> \file leap_frog_1st_flow_2d.F90
!> \authors                    
!> Martin CAMPOS PINTO (campos@ann.jussieu.fr) 

module sll_leap_frog_1st_flow_2d
#include "sll_working_precision.h"
#include "sll_assert.h"
  use sll_constants
  use sll_flow_base
  implicit none

  type, extends(flow_base_class) :: leap_frog_1st_flow_2d
  contains
    procedure, pass(flow) :: initialize => init_lf_1st_flow
    procedure, pass(flow) :: flow_at_xv => lf_1st_flow_at_xv
  end type leap_frog_1st_flow_2d

contains

  subroutine init_lf_1st_flow( flow, dt )
    class(leap_frog_1st_flow_2d), intent(inout)  :: flow
    sll_real64, intent(in)                       :: dt
    flow%dt = dt
  end subroutine init_lf_1st_flow

  subroutine lf_1st_flow_at_xv( flow, x,v, f_x,f_v )
    class(leap_frog_1st_flow_2d), intent(inout)       :: flow
    sll_real64, intent(in)     :: x
    sll_real64, intent(in)     :: v
    sll_real64, intent(out)    :: f_x
    sll_real64, intent(out)    :: f_v

    f_x = x + 0.5*flow%dt
    f_v = v
  end subroutine lf_1st_flow_at_xv

end module sll_leap_frog_1st_flow_2d
