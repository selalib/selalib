module sll_landau_2d_initializer
#include "sll_working_precision.h"
#include "sll_assert.h"
  use numeric_constants
  use sll_module_mapped_meshes_2d_base
  use sll_scalar_field_initializers_base
  implicit none

  type, extends(scalar_field_2d_initializer_base) :: init_landau_2d
    sll_real64 :: eps
    class(sll_mapped_mesh_2d_base), pointer :: mesh
    sll_real64 :: kx
  contains
    procedure, pass(init_obj) :: initialize => initialize_landau_2d
    procedure, pass(init_obj) :: f_of_x1x2  => f_x1x2_landau_2d
  end type init_landau_2d

contains

  subroutine initialize_landau_2d( init_obj, mesh, eps_val )
    class(init_landau_2d), intent(inout)  :: init_obj
    class(sll_mapped_mesh_2d_base), intent(in), target :: mesh
    sll_real64, intent(in), optional     :: eps_val
    if( present(eps_val) ) then
       init_obj%eps = eps_val
    else
       init_obj%eps = 0.01_f64 ! just some default value
    end if
    init_obj%mesh => mesh
    ! kx remains uninitialized because we need mesh information
  end subroutine initialize_landau_2d

  subroutine f_x1x2_landau_2d( init_obj, data_out )
    class(init_landau_2d), intent(inout)       :: init_obj
    sll_real64, dimension(:,:), intent(out)    :: data_out

    class(sll_mapped_mesh_2d_base), pointer    :: mesh
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
    num_pts1 = mesh%nc_eta1+1
    num_pts2 = mesh%nc_eta2+1
    kx = 2.0_f64*sll_pi/(mesh%x1_at_node(num_pts1,1) - mesh%x1_at_node(1,1))
    SLL_ASSERT( size(data_out,1) .ge. num_pts1 )
    SLL_ASSERT( size(data_out,2) .ge. num_pts2 )
    do j=1,num_pts2
       do i=1, num_pts1
          v = mesh%x2_at_node(i,j)
          x = mesh%x1_at_node(i,j)
          data_out(i,j) = &
               ( 1.0_f64 + eps*cos(kx*x) )/sqrt(2*sll_pi)*exp(-0.5_f64*v*v)
       end do
    end do
  end subroutine f_x1x2_landau_2d

end module sll_landau_2d_initializer
