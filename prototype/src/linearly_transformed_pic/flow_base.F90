!> \file flow_base.F90
!> \namespace sll_flow_base
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

  type, abstract :: flow_base_class
     sll_real64   :: dt
   contains
     procedure(flow_evaluate), pass(flow), deferred :: flow_at_xv ! should be 'flow_at_x' (see below)
  end type flow_base_class

  abstract interface
     ! MCP -- todo: use arguments (x,f_x), with arbitrary dimensions
     subroutine flow_evaluate(flow,x,v,f_x,f_v)
       import flow_base_class
       class(flow_base_class)   :: flow
       sll_real64,intent(in)    :: x
       sll_real64,intent(in)    :: v
       sll_real64,intent(out)   :: f_x
       sll_real64,intent(out)   :: f_v
     end subroutine flow_evaluate
  end interface

  
end module sll_flow_base
