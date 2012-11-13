module sll_landau_2d_initializer
#include "sll_working_precision.h"
#include "sll_assert.h"
  use numeric_constants
#ifdef STDF95
  use sll_module_mapped_meshes_2d
#else
  use sll_module_mapped_meshes_2d_base
#endif
  use sll_scalar_field_initializers_base
  implicit none

#ifdef STDF95
  type :: init_landau_2d
    type(sll_mapped_mesh_2d_discrete), pointer :: mesh
    sll_int32   :: data_position
#else
  type, extends(scalar_field_2d_initializer_base) :: init_landau_2d
    class(sll_mapped_mesh_2d_base), pointer :: mesh
#endif
    sll_real64 :: eps
    sll_real64 :: kx
    sll_int32 :: is_delta_f
#ifndef STDF95
  contains
    procedure, pass(init_obj) :: initialize => initialize_landau_2d
    procedure, pass(init_obj) :: f_of_x1x2  => f_x1x2_landau_2d
#endif
  end type init_landau_2d

contains

#ifdef STDF95
  subroutine init_landau_2d_initialize(init_obj, mesh, data_position, eps_val, kx_val, &
       is_delta_f)
    type(init_landau_2d), intent(inout)  :: init_obj
    type(sll_mapped_mesh_2d_discrete), intent(in), target :: mesh
#else
  subroutine initialize_landau_2d( init_obj, mesh, data_position, eps_val, kx_val, &
       is_delta_f)
    class(init_landau_2d), intent(inout)  :: init_obj
    class(sll_mapped_mesh_2d_base), intent(in), target :: mesh
#endif
    sll_int32 :: data_position
    sll_real64, intent(in), optional     :: eps_val
    sll_real64, intent(in), optional     :: kx_val
    sll_int32, intent(in), optional      :: is_delta_f

    init_obj%data_position = data_position
    if( present(eps_val) ) then
       init_obj%eps = eps_val
    else
       init_obj%eps = 0.01_f64 ! just some default value
    end if
    if( present(kx_val) ) then
       init_obj%kx = kx_val
    else
       init_obj%kx = 0.5_f64 ! just some default value
    end if
    if( present(is_delta_f) ) then
       init_obj%is_delta_f = is_delta_f
    else
       init_obj%is_delta_f = 1 !  default value is false
    end if
    init_obj%mesh => mesh
    ! kx remains uninitialized because we need mesh information
  end subroutine

#ifdef STDF95
  subroutine init_landau_2d_f_of_x1x2( init_obj, data_out )
    type(init_landau_2d), intent(inout)       :: init_obj
    type(sll_mapped_mesh_2d_discrete), pointer    :: mesh
#else
  subroutine f_x1x2_landau_2d( init_obj, data_out )
    class(init_landau_2d), intent(inout)       :: init_obj
    class(sll_mapped_mesh_2d_base), pointer    :: mesh
#endif
    sll_real64, dimension(:,:), intent(out)    :: data_out

    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: num_pts1
    sll_int32  :: num_pts2
    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: x
    sll_real64 :: v

    eps = init_obj%eps
    mesh => init_obj%mesh
    if (init_obj%data_position ==  NODE_CENTERED_FIELD) then
       num_pts1 = mesh%nc_eta1+1
       num_pts2 = mesh%nc_eta2+1
    else if (init_obj%data_position ==  CELL_CENTERED_FIELD) then
       num_pts1 = mesh%nc_eta1
       num_pts2 = mesh%nc_eta2
    end if
    kx = init_obj%kx
    SLL_ASSERT( size(data_out,1) .ge. num_pts1 )
    SLL_ASSERT( size(data_out,2) .ge. num_pts2 )
    do j=1,num_pts2
       do i=1, num_pts1
          if (init_obj%data_position ==  NODE_CENTERED_FIELD) then
#ifdef STDF95
             v = x2_at_node(mesh, i,j)
             x = x1_at_node(mesh, i,j)
#else
             v = mesh%x2_at_node(i,j)
             x = mesh%x1_at_node(i,j)
#endif
          else if (init_obj%data_position ==  CELL_CENTERED_FIELD) then
#ifdef STDF95
             v = x2_cell_discrete(mesh, i,j)
             x = x1_cell_discrete(mesh, i,j)
#else
             v = mesh%x2_at_cell(i,j)
             x = mesh%x1_at_cell(i,j)
#endif
          else
             print*, 'f_x1x2_landau_2d:',  init_obj%data_position, 'not defined'
          end if
          if (init_obj%is_delta_f==0) then ! delta_f code
             data_out(i,j) = eps*cos(kx*x) /sqrt(2*sll_pi)*exp(-0.5_f64*v*v)
          else  ! full f
             data_out(i,j) = &
               ( 1.0_f64 + eps*cos(kx*x) )/sqrt(2*sll_pi)*exp(-0.5_f64*v*v)
          endif
       end do
    end do
  end subroutine

end module sll_landau_2d_initializer
