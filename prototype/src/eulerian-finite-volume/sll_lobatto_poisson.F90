module sll_lobatto_poisson
#include "sll_working_precision.h"
   use sll_module_coordinate_transformations_2d
   use sll_common_coordinate_transformations
   use lobalap
   implicit none
   
   private
   
   type, public :: lobatto_poisson_solver
   
   end type lobatto_poisson_solver
   
   interface initialize
   module procedure initialize_lobatto_poisson
   end interface initialize
   
   interface solve
      module procedure solve_lobatto_poisson
   end interface solve
   
   interface delete
   module procedure delete_lobatto_poisson
   end interface delete
   
public :: initialize, solve, delete

contains

subroutine initialize_lobatto_poisson(this, tau)

   class(sll_coordinate_transformation_2d_base),pointer :: tau
   type(lobatto_poisson_solver) :: this
   procedure(map_function)      :: map0


  call init(30,10,2, map0)
  call assemb()
  call computeLU()

end subroutine initialize_lobatto_poisson

subroutine solve_lobatto_poisson(this, phi, rho)
  sll_real64, dimension(:,:) :: phi
  sll_real64, dimension(:,:) :: rho
  type(lobatto_poisson_solver) :: this

  call compute_phi()

end subroutine solve_lobatto_poisson

subroutine delete_lobatto_poisson(this)

   type(lobatto_poisson_solver) :: this

  call plotgmsh()
  call release()

end subroutine delete_lobatto_poisson


end module sll_lobatto_poisson
