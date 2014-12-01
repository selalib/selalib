!> \file flow_base.F90
!> \authors                    
!> Martin CAMPOS PINTO (campos@ann.jussieu.fr) 
!> \brief  
!> Essentially, the 'flow' has a single objective which is to provide
!> a uniform interface to flow mappings in arbitrary dimensions.

module sll_flow_base
#include "sll_working_precision.h"

  implicit none

  ! **************************************************************************
  !
  !                              Leap-frog flows
  !
  ! **************************************************************************

  type, abstract :: flow_base_class   ! MCP -- is that a class ?
     sll_real64  :: dt
   contains
     procedure(flow_evaluate), pass(flow), deferred :: flow_at_xv ! MCP -- for arbitrary dimensions, should be 'flow_at_x'
  end type flow_base_class
  

  abstract interface
     ! MCP -- here one should use arguments (x,f_x) with arbitrary dimensions (maybe later)
     subroutine flow_evaluate(flow,x,v,f_x,f_v)
       use sll_working_precision
       import flow_base_class
       class(flow_base_class)   :: flow
       sll_real64,intent(in)    :: x
       sll_real64,intent(in)    :: v
       sll_real64,intent(out)   :: f_x
       sll_real64,intent(out)   :: f_v
     end subroutine flow_evaluate
  end interface

  
end module sll_flow_base
