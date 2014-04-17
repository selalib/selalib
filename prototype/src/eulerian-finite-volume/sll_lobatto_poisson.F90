module sll_lobatto_poisson

#include "sll_working_precision.h"

   use sll_module_coordinate_transformations_2d
   use sll_common_coordinate_transformations
   use map_function_module, only: set_map_function
   use lobalap
   use sll_dg_fields

   implicit none
   
   private
   
   type, public :: lobatto_poisson_solver
      class(sll_coordinate_transformation_2d_analytic),pointer :: tau
      sll_int32  :: order
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

subroutine initialize_lobatto_poisson(this, tau, order)

   type(lobatto_poisson_solver) :: this
   class(sll_coordinate_transformation_2d_analytic),pointer :: tau
   sll_int32, optional :: order
   sll_int32 :: nx0
   sll_int32 :: ny0
   sll_int32 :: order0

   this%tau => tau
   nx0 = tau%mesh%num_cells1+1
   ny0 = tau%mesh%num_cells2+1

   call set_map_function(tau)

   if (present(order)) then
      call init(nx0,ny0,order)
   else
      call init(nx0,ny0,2)
   end if
   call assemb()
   call computeLU()

end subroutine initialize_lobatto_poisson

subroutine solve_lobatto_poisson(this)!, rhs)

  type(lobatto_poisson_solver) :: this

  call compute_phi()



end subroutine solve_lobatto_poisson

subroutine delete_lobatto_poisson(this)

   type(lobatto_poisson_solver) :: this

  call plotgmsh()
  call release()

end subroutine delete_lobatto_poisson


end module sll_lobatto_poisson
