module sll_m_lobatto_poisson

#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_m_coordinate_transformation_2d_base
  use sll_m_coordinate_transformations_2d
  use sll_m_common_coordinate_transformations
  use sll_m_map_function, only: set_map_function
  use sll_m_lobalap
  use sll_m_dg_fields

  implicit none
  
  type, public :: lobatto_poisson_solver
    sll_transformation, pointer :: tau
    sll_int32                   :: order
  end type lobatto_poisson_solver
   
  interface sll_create
    module procedure initialize_lobatto_poisson
  end interface sll_create
   
  interface sll_solve
    module procedure solve_lobatto_poisson
  end interface sll_solve
   
  interface sll_delete
    module procedure delete_lobatto_poisson
  end interface sll_delete

  public :: sll_create, sll_solve, sll_delete

  private

contains

subroutine initialize_lobatto_poisson(this, tau, order)

  type(lobatto_poisson_solver) :: this
  sll_transformation, pointer :: tau
  sll_int32, optional :: order
  sll_int32 :: nx0
  sll_int32 :: ny0

  this%tau => tau
  nx0 = tau%mesh%num_cells1
  ny0 = tau%mesh%num_cells2

  call set_map_function(tau)

  if (present(order)) then
     call init(nx0,ny0,order)
  else
     call init(nx0,ny0,2)
  end if
  call assemb()
  call computeLU()

end subroutine initialize_lobatto_poisson

subroutine solve_lobatto_poisson(this, rhs, ex, ey)

  type(lobatto_poisson_solver) :: this
  type(sll_dg_field_2d)        :: rhs
  type(sll_dg_field_2d)        :: ex
  type(sll_dg_field_2d)        :: ey

  SLL_ASSERT(this%order>0)
  call assemb_rhs(rhs%array)
  call compute_phi()
  call compute_electric_field(ex%array, ey%array)

end subroutine solve_lobatto_poisson

subroutine delete_lobatto_poisson(this)

  type(lobatto_poisson_solver) :: this

  call plotgmsh()
  call release()
  SLL_ASSERT(this%order>0)

end subroutine delete_lobatto_poisson

end module sll_m_lobatto_poisson
