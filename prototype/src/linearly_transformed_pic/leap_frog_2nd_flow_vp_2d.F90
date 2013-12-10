!> \file leap_frog_2nd_flow_2d.F90
!> \authors                    
!> Martin CAMPOS PINTO (campos@ann.jussieu.fr) 

module sll_leap_frog_2nd_flow_2d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_constants
  use sll_flow_base
  implicit none

  type, extends(flow_base_class) :: leap_frog_2nd_flow_2d
    sll_real64, dimension(:), pointer     :: electric_field
    sll_int32                             :: nc_poisson_mesh ! number of cells in the poisson mesh
    sll_real64                            :: xmin_poisson_mesh
    sll_real64                            :: dx_poisson_mesh
    sll_int32                             :: bc_type           ! periodic, open domain

  contains
    procedure, pass(flow) :: initialize => init_lf_2nd_flow
    procedure, pass(flow) :: flow_at_xv => lf_2nd_flow_at_xv
  end type leap_frog_2nd_flow_2d


#ifdef STDF95
   integer, parameter :: PERIODIC_E_FLOW = 0, OPEN_DOMAIN_E_FLOW = 1
#else
  enum, bind(C)
     enumerator :: PERIODIC_E_FLOW = 0, OPEN_DOMAIN_E_FLOW = 1
  end enum
#endif

contains

  subroutine init_lf_2nd_flow( flow, dt, nc_poisson_mesh, xmin_poisson_mesh, xmax_poisson_mesh, bc_type )
    class(leap_frog_2nd_flow_2d), intent(inout)  :: flow
    sll_real64, intent(in)                       :: dt
    sll_int32,  intent(in)                       :: nc_poisson_mesh
    sll_real64, intent(in)                       :: xmin_poisson_mesh
    sll_real64, intent(in)                       :: xmax_poisson_mesh
    sll_int32,  intent(in)                       :: bc_type
    sll_int32  :: ierr

    flow%dt                 = dt
    flow%xmin_poisson_mesh  = xmin_poisson_mesh
    flow%nc_poisson_mesh    = nc_poisson_mesh
    flow%dx_poisson_mesh    = (xmax_poisson_mesh-xmin_poisson_mesh)/nc_poisson_mesh
    flow%bc_type            = bc_type
    SLL_ALLOCATE( flow%electric_field(nc_poisson_mesh+1),  ierr)
  end subroutine init_lf_2nd_flow

  subroutine lf_2nd_flow_at_xv( flow, x,v, f_x,f_v )
    class(leap_frog_2nd_flow_2d), intent(inout)       :: flow
    sll_real64, intent(in)     :: x
    sll_real64, intent(in)     :: v
    sll_real64, intent(out)    :: f_x
    sll_real64, intent(out)    :: f_v
    sll_int32                  :: ix
    sll_int32                  :: ix_aux
    sll_real64                 :: normalized_x
    sll_real64                 :: theta_x
    sll_real64                 :: elec_field_at_x
    
    ! here I evaluate the electric field at x with affine interpolation, but I think there should be a class for the E field,
    ! and we should call something like electric_field%get_value(x)
    normalized_x = (x-flow%xmin_poisson_mesh)/flow%dx_poisson_mesh
    ix = 1+int( floor(normalized_x) )
    if( flow%bc_type .eq. PERIODIC_E_FLOW ) then    
      ix_aux = 1 + modulo( (ix-1), flow%nc_poisson_mesh )
    else
      ix_aux = ix  
    end if
    if( ix_aux >= 1 .and. ix_aux <= flow%nc_poisson_mesh ) then
      theta_x = ix-normalized_x
      elec_field_at_x = theta_x * flow%electric_field(ix_aux) + (1-theta_x) * flow%electric_field(ix_aux+1)    
    else
      elec_field_at_x = 0
    end if              
    f_v = v +     flow%dt*elec_field_at_x
    f_x = x + 0.5*flow%dt*f_v

  end subroutine lf_2nd_flow_at_xv

end module sll_leap_frog_2nd_flow_2d
