module sll_m_lobatto_poisson

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_dg_fields, only: &
    sll_t_dg_field_2d

  use sll_m_lobalap, only: &
    sll_s_assemb, &
    sll_s_assemb_rhs, &
    sll_s_compute_electric_field, &
    sll_s_compute_phi, &
    sll_s_computelu, &
    sll_s_init, &
    sll_s_plotgmsh, &
    sll_s_release

  use sll_m_map_function, only: &
    sll_s_set_map_function

  implicit none

  public :: &
    sll_t_lobatto_poisson_solver, &
    sll_o_create, &
    sll_o_delete, &
    sll_o_solve

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  type :: sll_t_lobatto_poisson_solver
    class(sll_c_coordinate_transformation_2d_base), pointer :: tau
    sll_int32                                             :: order
  end type sll_t_lobatto_poisson_solver
   
  interface sll_o_create
    module procedure initialize_lobatto_poisson
  end interface sll_o_create
   
  interface sll_o_solve
    module procedure solve_lobatto_poisson
  end interface sll_o_solve
   
  interface sll_o_delete
    module procedure delete_lobatto_poisson
  end interface sll_o_delete



contains

subroutine initialize_lobatto_poisson(this, tau, order)

  type(sll_t_lobatto_poisson_solver)                          :: this
  class(sll_c_coordinate_transformation_2d_base), pointer :: tau
  sll_int32, optional                                   :: order

  sll_int32 :: nx0
  sll_int32 :: ny0

  this%tau => tau
  nx0 = tau%mesh%num_cells1
  ny0 = tau%mesh%num_cells2

  call sll_s_set_map_function(tau)

  if (present(order)) then
     call sll_s_init(nx0,ny0,order)
  else
     call sll_s_init(nx0,ny0,2)
  end if
  call sll_s_assemb()
  call sll_s_computelu()

end subroutine initialize_lobatto_poisson

subroutine solve_lobatto_poisson(this, rhs, ex, ey)

  type(sll_t_lobatto_poisson_solver) :: this
  type(sll_t_dg_field_2d)        :: rhs
  type(sll_t_dg_field_2d)        :: ex
  type(sll_t_dg_field_2d)        :: ey

  SLL_ASSERT(this%order>0)
  call sll_s_assemb_rhs(rhs%array)
  call sll_s_compute_phi()
  call sll_s_compute_electric_field(ex%array, ey%array)

end subroutine solve_lobatto_poisson

subroutine delete_lobatto_poisson(this)

  type(sll_t_lobatto_poisson_solver) :: this

  call sll_s_plotgmsh()
  call sll_s_release()
  SLL_ASSERT(this%order>0)

end subroutine delete_lobatto_poisson

end module sll_m_lobatto_poisson
