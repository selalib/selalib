module sll_lobatto_poisson
#include "sll_working_precision.h"
   use sll_module_coordinate_transformations_2d
   use sll_common_coordinate_transformations
   use lobalap
   implicit none
   
   private
   
   class(sll_coordinate_transformation_2d_base),pointer :: tau

   type, public :: lobatto_poisson_solver
   !contains
      !procedure :: map => map_function
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

subroutine map_function( u, v, x, y) 
   real(8), intent(in) :: u, v
   real(8), intent(out) :: x, y

   x = tau%x1(u,v)
   y = tau%x2(u,v)

end subroutine map_function


subroutine initialize_lobatto_poisson(this, tau0)

   type(lobatto_poisson_solver) :: this
   class(sll_coordinate_transformation_2d_base),pointer :: tau0
   procedure(map_interface), pointer :: map0

   tau => tau0
   map0 => map_function
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
