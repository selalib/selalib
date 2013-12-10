module sll_simulation_base
  implicit none

  ! Basic signature for all simulations. These can be declared, initialized,
  ! executed and deleted. This abstract class provides the general interface
  ! to be implemented in all specific simulations. To consider: maybe this
  ! is general enough to not merit a distinction among the dimensionality of
  ! the simulation.

  type, abstract :: sll_simulation_base_class
   contains
     procedure(simulation_initializer), pass(sim), deferred :: &
          init_from_file
     procedure(simulation_execute), pass(sim), deferred :: run
  end type sll_simulation_base_class


  abstract interface
     subroutine simulation_initializer(sim, filename)
       import sll_simulation_base_class
       class(sll_simulation_base_class), intent(inout) :: sim
       character(len=*), intent(in)      :: filename
     end subroutine simulation_initializer
  end interface

  abstract interface
     subroutine simulation_execute(sim)
       import sll_simulation_base_class
       class(sll_simulation_base_class), intent(inout) :: sim
     end subroutine simulation_execute
  end interface

end module sll_simulation_base
