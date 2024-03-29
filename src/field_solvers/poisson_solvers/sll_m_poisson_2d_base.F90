!> @ingroup poisson_solvers
!> @brief
!> Module interface to solve Poisson equation in 2D
!> @details
!> Contains the abstract class to create a Poisson solver in 2D.
module sll_m_poisson_2d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   implicit none

   public :: &
      sll_c_poisson_2d_base, &
      sll_i_function_of_position

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !>
   type, abstract :: sll_c_poisson_2d_base

   contains

      !> solves \f$ -\Delta \phi_{ij} = \rho_{ij} \f$
      procedure(signature_compute_phi_from_rho_2d), deferred, pass(poisson) :: &
         compute_phi_from_rho

      !> solves \f$ -\Delta \phi_{ij} = \rho_{ij} \f$ and \f$ E_{ij} = \nabla  \phi_{ij} \f$
      procedure(signature_compute_E_from_rho_2d), deferred, pass(poisson) :: &
         compute_E_from_rho

      !> Compute the squared L_2 for given coefficients
      procedure(signature_norm_squared), deferred :: &
         l2norm_squared
      !> Compute the right hand side from a given function
      procedure(signature_update_dofs_function), deferred :: &
         compute_rhs_from_function

      procedure(signature_empty), deferred :: &
         free

   end type sll_c_poisson_2d_base

   abstract interface
      !> nd real function
      function sll_i_function_of_position(x)
         use sll_m_working_precision
         sll_real64             :: sll_i_function_of_position
         sll_real64, intent(in) :: x(:)
      end function sll_i_function_of_position
   end interface

#ifndef DOXYGEN_SHOULD_SKIP_THIS

   abstract interface

      ! solves -\Delta phi = rho in 2d

      subroutine signature_compute_phi_from_rho_2d(poisson, phi, rho)

         use sll_m_working_precision
         import sll_c_poisson_2d_base

         class(sll_c_poisson_2d_base)           :: poisson
         sll_real64, dimension(:, :), intent(in)  :: rho
         sll_real64, dimension(:, :), intent(out) :: phi

      end subroutine signature_compute_phi_from_rho_2d

   end interface

   abstract interface
      ! solves E = -\nabla Phi with -\Delta phi = rho in 2d
      subroutine signature_compute_E_from_rho_2d(poisson, E1, E2, rho)

         use sll_m_working_precision
         import sll_c_poisson_2d_base

         class(sll_c_poisson_2d_base)            :: poisson
         sll_real64, dimension(:, :), intent(in)  :: rho
         sll_real64, dimension(:, :), intent(out) :: E1
         sll_real64, dimension(:, :), intent(out) :: E2

      end subroutine signature_compute_E_from_rho_2d
   end interface

   !---------------------------------------------------------------------------!
   abstract interface
      function signature_norm_squared(poisson, coefs_dofs) result(r)
         use sll_m_working_precision
         import sll_c_poisson_2d_base
         class(sll_c_poisson_2d_base), intent(in) :: poisson !< Poisson solver object.
         sll_real64, intent(in)                    :: coefs_dofs(:, :) !< Values of the coefficient vectors for each DoF
         sll_real64                                 :: r
      end function signature_norm_squared
   end interface

   !---------------------------------------------------------------------------!
   abstract interface
      subroutine signature_update_dofs_function(poisson, func, coefs_dofs)
         use sll_m_working_precision
         import sll_c_poisson_2d_base
         import sll_i_function_of_position
         class(sll_c_poisson_2d_base)            :: poisson !< Maxwell solver object.
         procedure(sll_i_function_of_position)    :: func !< Function to be projected.
         sll_real64, intent(out)                  :: coefs_dofs(:) !< Coefficients of the projection.
      end subroutine signature_update_dofs_function
   end interface

   abstract interface

      ! Destructor

      subroutine signature_empty(poisson)
         import sll_c_poisson_2d_base
         class(sll_c_poisson_2d_base)   :: poisson

      end subroutine signature_empty

   end interface

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

end module sll_m_poisson_2d_base
