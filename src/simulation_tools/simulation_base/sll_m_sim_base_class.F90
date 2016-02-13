!> @ingroup simulations
!> Parent class module for simulation
module sll_m_sim_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_utilities, only: &
    sll_s_int2string

  implicit none

  public :: &
    sll_c_simulation_base_class

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Basic signature for all simulations. These can be declared, initialized,
  ! executed and deleted. This abstract class provides the general interface
  ! to be implemented in all specific simulations. To consider: maybe this
  ! is general enough to not merit a distinction among the dimensionality of
  ! the simulation.

  !> Parent class for simulation
  type, abstract :: sll_c_simulation_base_class
   contains
     !> Constructor
     procedure(simulation_initializer), pass(sim), deferred :: &
          init_from_file
     !> Run method
     procedure(simulation_execute), pass(sim), deferred :: run
  end type sll_c_simulation_base_class


#ifndef DOXYGEN_SHOULD_SKIP_THIS

  abstract interface
     !> Class constructor for simulation
     !> \param[inout] simulation class
     !> \param[in]    input data file 
     subroutine simulation_initializer(sim, filename)
       import sll_c_simulation_base_class
       class(sll_c_simulation_base_class), intent(inout) :: sim
       character(len=*), intent(in)      :: filename
     end subroutine simulation_initializer
  end interface

  abstract interface
     !> Class run method for simulation
     !> \param[inout] simulation class
     subroutine simulation_execute(sim)
       import sll_c_simulation_base_class
       class(sll_c_simulation_base_class), intent(inout) :: sim
     end subroutine simulation_execute
  end interface

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

end module sll_m_sim_base
