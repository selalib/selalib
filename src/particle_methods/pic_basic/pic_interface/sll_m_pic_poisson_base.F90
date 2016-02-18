!> @ingroup pic_interface
!> @author Katharina Kormann, IPP
!> @brief Base class for Poisson solver for particle methods.
!> @details This base class gives an abstract interface to the basic functions for accumulation of charge  as well as the evaluation of a function at particle positions.
!> @todo Try to integrate also general elliptic and quasi neutral solvers.
!> @todo Allow for more than one weights
!> @todo Functions to compute derivatives
module sll_m_pic_poisson_base

#include "sll_working_precision.h"

  use sll_m_poisson_2d_base, only : &
       sll_i_fucntion_of_position

  implicit none

  public :: sll_c_pic_poisson

  private

  !> Basic type of Poisson solver for PIC simulations
  type, public, abstract :: sll_c_pic_poisson
     sll_int32              :: dim !< Dimension
     sll_int32              :: no_weights = 1 !< Number of weights used for accumulation (one per default)
     
   contains
     procedure(add_single), deferred            :: add_charge_single !< Add the contribution of one particle to the charge density
     procedure                                  :: add_charge_vector !< Add the contribution of a number of particles to the charge density
     procedure(eval_component_single), deferred :: evaluate_field_single !< Evaluate given components of the field at a given position.
     procedure                                  :: evaluate_field_vector !< Evaluate given components of the field at given positions.
     procedure(eval_single), deferred           :: evaluate_rho_single !< Evaluate charge density at given position.
     procedure                                  :: evaluate_rho_vector !< Evaluate charge density at given positions.
     procedure(eval_single), deferred           :: evaluate_phi_single !< Evaluate potential at given position.
     procedure                                  :: evaluate_phi_vector !< Evaluate potential at given positions.
     procedure(empty), deferred                 :: reset !< Reset the accumulated charge to zero
     procedure(empty), deferred                 :: solve !< Solve for the electric potential and field
     procedure(empty), deferred                 :: solve_phi !< Solve for phi
     procedure(empty), deferred                 :: solve_fields !< solve for electric field
     procedure(compute_energy), deferred        :: compute_field_energy !< Compute the L2 norm of one field component
     !procedure(update_dofs_function), deferred  :: compute_rhs_from_function
     !procedure(update_dofs_function), deferred  :: l2projection
     procedure(linear_combination), deferred    :: add_analytic_charge !< Set charge as linear combination of previously accumulated charge and previously set analytic charge.
     procedure(update_dofs_function), deferred  :: set_analytic_charge !< Set the value of the analytic charge contribution from a given function.

     procedure(empty), deferred                  :: free


     generic :: add_charge => add_charge_single, add_charge_vector
     generic :: evaluate_field => evaluate_field_single, evaluate_field_vector
     generic :: evaluate_rho => evaluate_rho_single, evaluate_rho_vector
     generic :: evaluate_phi => evaluate_phi_single, evaluate_phi_vector 

  end type sll_c_pic_poisson

  
  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_single(self, position, marker_charge) 
       use sll_m_working_precision
       import sll_c_pic_poisson
       class (sll_c_pic_poisson), intent( inout ) :: self !< Pic Poisson solver object
       sll_real64,                intent( in ) :: position(self%dim) !< Position of the particle
       sll_real64,                intent( in ) :: marker_charge !< Particle Weight times charge
     end subroutine add_single
  end interface

!!$  !---------------------------------------------------------------------------!
!!$  abstract interface
!!$     subroutine add_vector(self,n_part, position, weight) 
!!$       use sll_m_working_precision
!!$       import sll_c_pic_poisson
!!$       class (sll_c_pic_poisson), intent( in ) :: self !< Pic Poisson solver object
!!$       sll_int32,                 intent( in ) :: n_part !< Number of particles whos positions are given
!!$       sll_real64,                intent( in ) :: position(self%dim, n_part) !< Position of the particle
!!$       sll_real64,                intent( in ) :: weight(n_part) !< Weight of the particle
!!$     end subroutine add_vector
!!$  end interface


  !---------------------------------------------------------------------------!
  abstract interface
     subroutine eval_single(self, position, func_value)
       use sll_m_working_precision
       import sll_c_pic_poisson
       class (sll_c_pic_poisson), intent( inout ) :: self !< Pic Poisson solver object
       sll_real64,                intent( in ) :: position(self%dim) !< Position of the particle
       sll_real64,                intent( out) :: func_value !< Value(s) of the electric fields at given position
     end subroutine eval_single
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine eval_component_single(self, position, components, func_value)
       use sll_m_working_precision
       import sll_c_pic_poisson
       class (sll_c_pic_poisson), intent( inout ) :: self !< Pic Poisson solver object
       sll_real64,                intent( in ) :: position(self%dim) !< Position of the particle
       sll_int32,                 intent( in ) :: components(:) !< Components of the field to be computed
       sll_real64,                intent( out) :: func_value(:) !< Value(s) of the electric fields at given position
     end subroutine eval_component_single
  end interface


!!$  !---------------------------------------------------------------------------!
!!$  abstract interface
!!$     subroutine eval_vector(self, position, n_part, func_value)
!!$       use sll_m_working_precision
!!$       import sll_c_pic_poisson
!!$       class (sll_c_pic_poisson), intent( in ) :: self !< Pic Poisson solver object
!!$       sll_real64,                intent( in ) :: position(self%dim, n_part) !< Position of the particles
!!$       sll_int32,                 intent( in ) :: n_part !< Number of particles whos positions are given
!!$       sll_real64,                intent( out) :: func_value(self%dim, n_part) !< Value(s) of the electric fields at given position
!!$     end subroutine eval_vector
!!$  end interface
!!$
!!$  !---------------------------------------------------------------------------!
!!$  abstract interface
!!$     subroutine eval_component_vector(self, position, n_part, components,func_value)
!!$       use sll_m_working_precision
!!$       import sll_c_pic_poisson
!!$       class (sll_c_pic_poisson), intent( in ) :: self !< Pic Poisson solver object
!!$       sll_real64,                intent( in ) :: position(self%dim, n_part) !< Position of the particle
!!$       sll_int32,                 intent( in ) :: n_part !< Number of particles whos positions are given
!!$       sll_int32,                 intent( in ) :: components(:) !< Components of the field to be computed
!!$       sll_real64,                intent( out) :: func_value(:,:) !< Value(s) of the electric fields at given positions
!!$     end subroutine eval_component_vector
!!$  end interface


  !---------------------------------------------------------------------------!
  abstract interface
     subroutine empty(self)
       import sll_c_pic_poisson
       class (sll_c_pic_poisson), intent( inout ) :: self
     end subroutine empty
  end interface


  !---------------------------------------------------------------------------!
  abstract interface
     function compute_energy(self, component) result(energy)
       use sll_m_working_precision
       import sll_c_pic_poisson
       class (sll_c_pic_poisson), intent( in ) :: self !< PIC Poisson solver object
       sll_int32, intent( in ) :: component !< Component of the electric field for which the energy should be computed
       sll_real64 :: energy !< L2 norm squarred of 
     end function compute_energy
  end interface
  

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine update_dofs_function(self, func)
       use sll_m_working_precision
       import sll_c_pic_poisson
       import sll_i_fucntion_of_position
       class( sll_c_pic_poisson ), intent( inout )    :: self !< PIC Poisson solver object.
       procedure(sll_i_fucntion_of_position)                :: func !< Function to be projected.
     end subroutine update_dofs_function
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine linear_combination(self, factor_present, factor_analytic)
       use sll_m_working_precision
       import sll_c_pic_poisson
       class( sll_c_pic_poisson ), intent( inout ) :: self !< PIC Poisson solver object
       sll_real64, intent( in ) :: factor_present !< Factor to multiply accumulated charge with
       sll_real64, intent( in ) :: factor_analytic !< Factor to multiply added analytic charge with

     end subroutine linear_combination
  end interface

contains

    !---------------------------------------------------------------------------!
    !< Add the contribution of \a n_part particles to the charge density. Per default it is implemented as a loop over the implementation for a single particle but can be overwritten if necessary.
  subroutine add_charge_vector(self,n_part, position, marker_charge) 
    class (sll_c_pic_poisson), intent( inout ) :: self !< Pic Poisson solver object
    sll_int32,                 intent( in ) :: n_part !< Number of particles whos positions are given
    sll_real64,                intent( in ) :: position(self%dim, n_part) !< Position of the particle
    sll_real64,                intent( in ) :: marker_charge(n_part) !< Particle weights times charge

    !local variables
    sll_int32  :: i_part

    do i_part = 1, n_part
       call self%add_charge_single( position(:,i_part), marker_charge(i_part))
    end do

  end subroutine add_charge_vector


  !---------------------------------------------------------------------------!
  !< Evaluation of charge density for a vector of \a n_part particles. Per default it is implemented as a loop over the implementation for a single particle but can be overwritten if necessary.
  subroutine evaluate_rho_vector(self, position, n_part, func_value)
    class (sll_c_pic_poisson), intent( inout ) :: self !< Pic Poisson solver object
    sll_real64,                intent( in ) :: position(:, :) !< Position of the particles (size (self%dim, n_part))
    sll_int32,                 intent( in ) :: n_part !< Number of particles whos positions are given
    sll_real64,                intent( out) :: func_value(n_part) !< Value(s) of the electric fields at given position

    !local variables
    sll_int32 :: i_part

    do i_part = 1, n_part
       call self%evaluate_rho_single&
            (position(:, i_part), func_value(i_part))
    end do
    

  end subroutine evaluate_rho_vector

  !---------------------------------------------------------------------------!
  !< Evaluation of potential for a vector of \a n_part particles. Per default it is implemented as a loop over the implementation for a single particle but can be overwritten if necessary.
  subroutine evaluate_phi_vector(self, position, n_part, func_value)
    class (sll_c_pic_poisson), intent( inout ) :: self !< Pic Poisson solver object
    sll_real64,                intent( in ) :: position(:,:) !< Position of the particles(size (self%dim, n_part))
    sll_int32,                 intent( in ) :: n_part !< Number of particles whos positions are given
    sll_real64,                intent( out) :: func_value(n_part) !< Value(s) of the electric fields at given position

    !local variables
    sll_int32 :: i_part

    do i_part = 1, n_part
       call self%evaluate_phi_single&
            (position(:, i_part), func_value(i_part))
    end do
    

  end subroutine evaluate_phi_vector

  
  !---------------------------------------------------------------------------!
  subroutine evaluate_field_vector(self, position, n_part, components,func_value)
    class (sll_c_pic_poisson), intent( inout ) :: self !< Pic Poisson solver object
    sll_real64,                intent( in ) :: position(:,:) !< Position of the particle (size (self%dim, n_part))
    sll_int32,                 intent( in ) :: n_part !< Number of particles whos positions are given
    sll_int32,                 intent( in ) :: components(:) !< Components of the field to be computed
    sll_real64,                intent( out) :: func_value(:,:) !< Value(s) of the electric fields at given positions
    
    !local variables
    sll_int32 :: i_part

    do i_part = 1, n_part
       call self%evaluate_field_single&
            (position(:, i_part), components, func_value(:,i_part))
    end do
    

  end subroutine evaluate_field_vector
  

end module sll_m_pic_poisson_base
