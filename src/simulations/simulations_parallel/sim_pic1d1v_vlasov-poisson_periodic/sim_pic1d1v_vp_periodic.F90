module sll_module_simulation_pic1d1v_vp_periodic

#include "sll_working_precision.h"
#include "sll_errors.h"
  use sll_simulation_base, only: sll_simulation_base_class
  implicit none

!==============================================================================
  
  type, extends( sll_simulation_base_class ) :: &
      sll_simulation_pic1d1v_vp_periodic

  contains
    procedure :: run            => run_fake
    procedure :: init_from_file => init_fake

  end type sll_simulation_pic1d1v_vp_periodic

!==============================================================================
contains
!==============================================================================

  subroutine run_fake( sim )
    class( sll_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
    character( len=64 ), parameter :: this_sub_name = "run_fake"
    SLL_WARNING( this_sub_name, "'run' method not implemented" )
  end subroutine run_fake

  subroutine init_fake( sim, filename )
    class( sll_simulation_pic1d1v_vp_periodic ), intent( inout ) :: sim
    character( len=* )                         , intent( in    ) :: filename
    character( len=64 ), parameter :: this_sub_name = "init_fake"
    SLL_WARNING( this_sub_name, "'init_from_file' method not implemented" )
  end subroutine init_fake

end module sll_module_simulation_pic1d1v_vp_periodic
