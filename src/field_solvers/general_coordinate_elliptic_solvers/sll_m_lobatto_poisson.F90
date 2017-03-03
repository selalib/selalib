module sll_m_lobatto_poisson

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_cubic_spline_interpolator_2d, only: &
    sll_f_new_cubic_spline_interpolator_2d

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  use sll_m_scalar_field_2d, only: &
    sll_f_new_scalar_field_2d_discrete

  use sll_m_scalar_field_2d_base, only: &
    sll_c_scalar_field_2d_base

  use sll_m_dg_fields, only: &
    sll_t_dg_field_2d

  use sll_m_lobalap, only: &
    sll_i_2a_func, &
    sll_s_assemb, &
    ! sll_s_assemb_rhs, &
    ! sll_s_compute_electric_field, &
    sll_s_compute_phi, &
    sll_s_computelu, &
    sll_s_init, &
    sll_s_plotgmsh, &
    sll_s_release

  use sll_m_map_function, only: &
    sll_s_set_map_function

  implicit none

  public :: &
       sll_f_new_lobatto_poisson, &
       sll_t_lobatto_poisson_solver, &
       sll_o_create, &
       sll_o_delete, &
       sll_o_solve

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type :: sll_t_lobatto_poisson_solver
    class(sll_c_coordinate_transformation_2d_base), pointer :: transf
    class(sll_c_scalar_field_2d_base),              pointer :: rho_field
    class(sll_c_interpolator_2d),                   pointer :: interp_rho
    sll_int32                                               :: order
  contains
    procedure, pass(solver) :: initialize => initialize_lobatto_poisson
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


    ! fonction donnant le potentiel exact (debug) et/ou les conditions aux limites
  function potexact(x,y)
    implicit none
    sll_real64,intent(in) :: x,y
    sll_real64 :: potexact
    potexact=x*x+y*y
    ! potexact=0.0_f64+x-x+y-y
  end function potexact

  function sll_f_new_lobatto_poisson(transf, order, &
       rho_values, &
       interp_rho, &
       bc_eta1_left, &
       bc_eta1_right, &
       bc_eta2_left, &
       bc_eta2_right)

    type(sll_t_lobatto_poisson_solver),             pointer :: &
         sll_f_new_lobatto_poisson
    class(sll_c_coordinate_transformation_2d_base), target  :: transf
    class(sll_c_interpolator_2d),                   target  :: interp_rho
    sll_real64, dimension(:,:), intent(in) :: rho_values
    sll_int32, intent(in), optional  :: order
    sll_int32, intent(in)  :: bc_eta1_left
    sll_int32, intent(in)  :: bc_eta1_right
    sll_int32, intent(in)  :: bc_eta2_left
    sll_int32, intent(in)  :: bc_eta2_right
    sll_int32 :: ierr

    SLL_ALLOCATE(sll_f_new_lobatto_poisson, ierr)
    call initialize_lobatto_poisson( &
         sll_f_new_lobatto_poisson, &
         transf, order, &
         rho_values, &
         interp_rho, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
  end function sll_f_new_lobatto_poisson


subroutine initialize_lobatto_poisson(solver, transf, order, &
       rho_values, &
       interp_rho, &
       bc_eta1_left, &
       bc_eta1_right, &
       bc_eta2_left, &
       bc_eta2_right)

    class(sll_t_lobatto_poisson_solver)                     :: solver
    class(sll_c_coordinate_transformation_2d_base), target  :: transf
    class(sll_c_interpolator_2d),                   target  :: interp_rho
    sll_real64, dimension(:,:), intent(in) :: rho_values
    sll_int32, intent(in), optional  :: order
    sll_int32, intent(in)  :: bc_eta1_left
    sll_int32, intent(in)  :: bc_eta1_right
    sll_int32, intent(in)  :: bc_eta2_left
    sll_int32, intent(in)  :: bc_eta2_right

  sll_int32 :: nx0
  sll_int32 :: ny0

  ! Setting mesh coo. tranformation....................
  solver%transf => transf
  nx0 = transf%mesh%num_cells1 + 1
  ny0 = transf%mesh%num_cells2 + 1
  call sll_s_set_map_function(transf)
  ! ...................................................

  ! Init tabs and gauss points computation.............
    if (present(order)) then
     call sll_s_init(nx0, ny0, order)
     solver%order = order
  else
     call sll_s_init(nx0, ny0, 2)
     solver%order = 2
  end if
  ! ...................................................

  ! Creating rho interp and field from values..........
  solver%interp_rho => interp_rho
  solver%rho_field  => sll_f_new_scalar_field_2d_discrete( &
       "rho_field", &
       solver%interp_rho, &
       transf, &
       bc_eta1_left, &
       bc_eta1_right, &
       bc_eta2_left, &
       bc_eta2_right)
  call solver%rho_field%set_field_data(rho_values)
  call solver%rho_field%update_interpolation_coefficients( )
  !....................................................

  call sll_s_assemb(solver%rho_field, potexact)
  call sll_s_computelu()

end subroutine initialize_lobatto_poisson


subroutine solve_lobatto_poisson(this, rhs, ex, ey)

  type(sll_t_lobatto_poisson_solver) :: this
  type(sll_t_dg_field_2d)            :: rhs
  type(sll_t_dg_field_2d)            :: ex
  type(sll_t_dg_field_2d)            :: ey

  SLL_ASSERT(this%order>0)
  ! call sll_s_assemb_rhs(rhs%array)
  call sll_s_compute_phi()
  ! call sll_s_compute_electric_field(ex%array, ey%array)

end subroutine solve_lobatto_poisson

subroutine delete_lobatto_poisson(this)

  type(sll_t_lobatto_poisson_solver) :: this

  call sll_s_plotgmsh()
  call sll_s_release()
  SLL_ASSERT(this%order>0)

end subroutine delete_lobatto_poisson


end module sll_m_lobatto_poisson
