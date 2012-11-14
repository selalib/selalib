! Essentially, the 'flow' has a single objective which is to provide
! a uniform interface to a flow mapping in a given number of dimensions.

module sll_flow_base
#include "sll_working_precision.h"
  implicit none

  ! **************************************************************************
  !
  !                              Leap-frog flows
  !
  ! **************************************************************************

  type, abstract :: flow_base
     sll_real64   :: dt
   contains
     procedure( *** MCP - don't know what to write here *** ), deferred, pass :: flow_at_xv
     ! need to declare the initialize procedure ??
  end type flow_base

  ! do we need an interface ?
  abstract interface
  end interface

  
end module sll_flow_base
