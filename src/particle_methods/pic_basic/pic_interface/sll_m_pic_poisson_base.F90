!> @ingroup pic_interface
!> @author Katharina Kormann, IPP
!> @brief Base class for Poisson solver for particle methods.
!> @details This base class gives an abstract interface to the basic functions for accumulation of charge  as well as the evaluation of a function at particle positions.
!> @todo Try to integrate also general elliptic and quasi neutral solvers.
module sll_m_pic_poisson_base

#include "sll_working_precision.h"

  implicit none

  public :: sll_c_pic_poisson

  private

  !> Basic type of Poisson solver for PIC simulations
  type, public, abstract :: sll_c_pic_poisson
     sll_int32              :: dim !< Dimension
     
   contains
     procedure(add_single), deferred        :: add_charge_single !< Add the contribution of one particle to the charge density
     procedure(add_vector), deferred        :: add_charge_vector !< Add the contribution of a number of particles to the charge density
     procedure(eval_single), deferred       :: evaluate_field_single
     procedure(eval_vector), deferred       :: evaluate_field_vector
     procedure(empty), deferred             :: reset !< Reset the accumulated charge to zero
     procedure(empty), deferred             :: solve !< Solve for the electric field
     procedure(update_single), deferred     :: add_shape_factor_single !< Prepare for the accumulation by computing the shape factors (add just a single particle)
     procedure(update_vector), deferred     :: add_shape_factor_vector !< Prepare for the accumulation by computing the shape factors (add an array of particles)

     generic :: add_charge => add_charge_single, add_charge_vector
     generic :: evaluate_field => evaluate_field_single, evaluate_field_vector
     generic :: add_shape_factor => add_shape_factor_single, add_shape_factor_vector

  end type sll_c_pic_poisson

  
  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_single(this, weight, position) 
       use sll_m_working_precision
       import sll_c_pic_poisson
       class (sll_c_pic_poisson), intent( in ) :: this !< Pic Poisson solver object
       sll_real64,                intent( in ) :: weight !< Weight of the particle
       sll_real64, optional,      intent( in ) :: position(this%dim) !< Position of the particle
     end subroutine add_single
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_vector(this,n_part, weight, position) 
       use sll_m_working_precision
       import sll_c_pic_poisson
       class (sll_c_pic_poisson), intent( in ) :: this !< Pic Poisson solver object
       sll_int32,                 intent( in ) :: n_part !< Number of particles whos positions are given
       sll_real64,                intent( in ) :: weight(n_part) !< Weight of the particle
       sll_real64, optional,      intent( in ) :: position(this%dim, n_part) !< Position of the particle
     end subroutine add_vector
  end interface


  !---------------------------------------------------------------------------!
  abstract interface
     pure function eval_single(this, position) result(field_value)
       use sll_m_working_precision
       import sll_c_pic_poisson
       class (sll_c_pic_poisson), intent( in ) :: this !< Pic Poisson solver object
       sll_real64,                intent( in ) :: position(this%dim) !< Position of the particle
       sll_real64                              :: field_value(this%dim) !< Value(s) of the electric fields at given position
     end function eval_single
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     pure function eval_vector(this, position, n_part) result(field_value)
       use sll_m_working_precision
       import sll_c_pic_poisson
       class (sll_c_pic_poisson), intent( in ) :: this !< Pic Poisson solver object
       sll_real64,                intent( in ) :: position(this%dim, n_part) !< Position of the particles
       sll_int32,                 intent( in ) :: n_part !< Number of particles whos positions are given
       sll_real64                              :: field_value(this%dim, n_part) !< Value(s) of the electric fields at given position
     end function eval_vector
  end interface


  !---------------------------------------------------------------------------!
  abstract interface
     subroutine empty(this)
       import sll_c_pic_poisson
       class (sll_c_pic_poisson), intent( inout ) :: this
     end subroutine empty
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine update_single(this, position)
       use sll_m_working_precision
       import sll_c_pic_poisson
       class( sll_c_pic_poisson ), intent(inout) :: this !< Pic Poisson solver object.
       sll_real64, intent(in)  :: position(this%dim) !< Position of the particle
     end subroutine update_single
  end interface
  
  !---------------------------------------------------------------------------!
  abstract interface
     subroutine update_vector(this, position, n_part)
       use sll_m_working_precision
       import sll_c_pic_poisson
       class( sll_c_pic_poisson ), intent(inout) :: this !< Pic Poisson solver object.
       sll_real64, intent(in)  :: position(this%dim,n_part) !< Position of the particles
       sll_int32, intent(in)   :: n_part !< Number of particles whos positions are given
     end subroutine update_vector
  end interface
  


end module sll_m_pic_poisson_base
