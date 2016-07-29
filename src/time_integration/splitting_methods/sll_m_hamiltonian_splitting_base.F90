!> @ingroup time_integration
!> @author Katharina Kormann, IPP
!> @brief Base class for Hamiltonian splittings.
!> @details  Contains deferred function for strang splitting, lie splitting and lie splitting with oposite ordering of the split-steps.
!> Moreover, composition methods based on lie and strang splitting are implemented (cf. Hairer, Lubich, Wanner, Geometric numeric integration.
module sll_m_hamiltonian_splitting_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_c_hamiltonian_splitting_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Type for Hamiltonian splittings
  type, abstract :: sll_c_hamiltonian_splitting_base

   contains 
     procedure(splitting), deferred :: lie_splitting  !< first order Lie splitting
     procedure(splitting), deferred :: lie_splitting_back !< first order Lie splitting with oposite ordering of the split steps
     procedure(splitting), deferred :: strang_splitting !< second order strang splitting
     procedure                      :: splitting_fourth !< fourth order composition method based on 3 steps of the Strang splitting
     procedure                      :: splitting_fourth_10steps !< fourth order composition method based on 5+5 steps of the two Lie splittings
     procedure                      :: splitting_second_4steps !< second order composition method based on 2+2 steps of the two Lie splittings
     procedure(empty),     deferred :: free
  end type sll_c_hamiltonian_splitting_base


  abstract interface
     subroutine splitting(self, dt, number_steps)
       use sll_m_working_precision
       import sll_c_hamiltonian_splitting_base
       class(sll_c_hamiltonian_splitting_base), intent(inout)   :: self !< time splitting object
       sll_real64,                              intent(in)      :: dt !< time step size
       sll_int32,                               intent(in)      :: number_steps !< number of time steps
     end subroutine splitting
  end interface

  abstract interface
     subroutine empty(self)
       import sll_c_hamiltonian_splitting_base
       class(sll_c_hamiltonian_splitting_base), intent(inout)   :: self !< time splitting object
     end subroutine empty
  end interface
  
contains

  subroutine splitting_fourth( self, dt, number_steps )
    class(sll_c_hamiltonian_splitting_base), intent(inout)   :: self !< time splitting object
    sll_real64,                              intent(in)      :: dt !< time step size
    sll_int32,                               intent(in)      :: number_steps !< number of time steps

    sll_real64 :: dt1, dt2
    sll_int32  :: j

    dt1 = dt / (2.0_f64 - 2.0_f64**(1.0_f64/3.0_f64))
    dt2 = - 2.0_f64**(1.0_f64/3.0_f64)* dt1

    do j=1,number_steps

       call self%strang_splitting( dt1, 1 )
       call self%strang_splitting( dt2, 1 )
       call self%strang_splitting( dt1, 1 )
       
    end do

    
  end subroutine splitting_fourth



  subroutine splitting_fourth_10steps( self, dt, number_steps )
    class(sll_c_hamiltonian_splitting_base), intent(inout)   :: self !< time splitting object
    sll_real64,                              intent(in)      :: dt !< time step size
    sll_int32,                               intent(in)      :: number_steps !< number of time steps

    sll_real64 :: a1, a2, a3, a4, a5
    sll_int32  :: j


    a1 = (146.0_f64 + 5.0_f64*sqrt(19.0_f64))/540.0_f64
    a2 = (-2.0_f64+10.0_f64*sqrt(19.0_f64))/135.0_f64
    a3 = 0.2_f64
    a4 = (-23.0_f64-20.0_f64*sqrt(19.0_f64))/270.0_f64
    a5 = (14.0_f64-sqrt(19.0_f64))/108.0_f64

    do j=1,number_steps

       call self%lie_splitting_back( a5*dt,1 )
       call self%lie_splitting( a1*dt,1 )
       call self%lie_splitting_back( a4*dt,1 )
       call self%lie_splitting( a2*dt,1 )
       call self%lie_splitting_back( a3*dt,1 )
       call self%lie_splitting( a3*dt,1 )
       call self%lie_splitting_back( a2*dt,1 )
       call self%lie_splitting( a4*dt,1 )
       call self%lie_splitting_back( a1*dt,1 )
       call self%lie_splitting( a5*dt,1 )       
    end do

    
  end subroutine splitting_fourth_10steps


  subroutine splitting_second_4steps( self, dt, number_steps )
    class(sll_c_hamiltonian_splitting_base), intent(inout)   :: self !< time splitting object
    sll_real64,                              intent(in)      :: dt !< time step size
    sll_int32,                               intent(in)      :: number_steps !< number of time steps

    sll_real64 :: a1, a2
    sll_int32  :: j


    a1 = 0.1932_f64!0.25_f64
    a2 = 0.5_f64 - a1

    do j=1,number_steps

       call self%lie_splitting_back( a1*dt,1 )
       call self%lie_splitting( a2*dt,1 )
       call self%lie_splitting_back( a2*dt,1 )
       call self%lie_splitting( a1*dt,1 )     
    end do
    
  end subroutine splitting_second_4steps
  
  
end module sll_m_hamiltonian_splitting_base
