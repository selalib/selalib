!> \file leap_frog_2nd_flow_2d.F90
!> \namespace sll_leap_frog_2nd_flow_2d
!> \authors                    
!> Martin CAMPOS PINTO (campos@ann.jussieu.fr) 

module sll_leap_frog_2nd_flow_2d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use numeric_constants
  use sll_flow_base
  implicit none

  type, extends(flow_base_class) :: leap_frog_2nd_flow_2d
    sll_real64, dimension(:), pointer     :: electric_field
    sll_int32                             :: nc_poisson_mesh ! number of cells in the poisson mesh
    sll_real64                            :: xmin_poisson_mesh
    sll_real64                            :: dx_poisson_mesh

  contains
    procedure, pass(flow) :: initialize => init_lf_2nd_flow
    procedure, pass(flow) :: flow_at_xv => lf_2nd_flow_at_xv
  end type leap_frog_2nd_flow_2d

contains

  subroutine init_lf_2nd_flow( flow, dt, nc_poisson_mesh, xmin_poisson_mesh, xmax_poisson_mesh)
    class(leap_frog_2nd_flow_2d), intent(inout)  :: flow
    sll_real64, intent(in)                       :: dt
    sll_int32,  intent(in)                       :: nc_poisson_mesh
    sll_real64, intent(in)                       :: xmin_poisson_mesh
    sll_real64, intent(in)                       :: xmax_poisson_mesh
    sll_int32  :: ierr

    flow%dt = dt
    flow%xmin_poisson_mesh = xmin_poisson_mesh
    flow%nc_poisson_mesh   = nc_poisson_mesh
    flow%dx_poisson_mesh   = (xmax_poisson_mesh-xmin_poisson_mesh)/nc_poisson_mesh
    SLL_ALLOCATE( flow%electric_field(nc_poisson_mesh+1),  ierr)
  end subroutine init_lf_2nd_flow

  subroutine lf_2nd_flow_at_xv( flow, x,v, f_x,f_v )
    class(leap_frog_2nd_flow_2d), intent(inout)       :: flow
    sll_real64, intent(in)     :: x
    sll_real64, intent(in)     :: v
    sll_real64, intent(out)    :: f_x
    sll_real64, intent(out)    :: f_v
    sll_int32                  :: ix
    sll_real64                 :: dt
    sll_real64                 :: normalized_x
    sll_real64                 :: theta_x
    sll_real64                 :: elec_field_at_x
    
    ! here I evaluate the electric field at x with affine interpolation, but I think there should be a class for the E field,
    ! and we should call something like electric_field%get_value(x)
    normalized_x = (x-flow%xmin_poisson_mesh)/flow%dx_poisson_mesh
    ix = 1+int( floor(normalized_x) )
    SLL_ASSERT( ix >= 1 )
    SLL_ASSERT( ix <= flow%nc_poisson_mesh )
    theta_x = ix+1-normalized_x
    elec_field_at_x = theta_x * flow%electric_field(ix) + (1-theta_x) * flow%electric_field(ix+1)
        
    dt = flow%dt 
    f_v = v + dt*elec_field_at_x
    f_x = x + 0.5*dt*f_v

  end subroutine lf_2nd_flow_at_xv

end module sll_leap_frog_2nd_flow_2d
