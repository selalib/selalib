!> Initializes a Gaussian of the form 
!> exp -((x-xc)**2/(2*sigma_x)**2 + (y-yc)**2/(2*sigma_y)**2)
!--------------------------------------------
module sll_gaussian_2d_initializer
#include "sll_working_precision.h"
#include "sll_assert.h"
  use sll_constants
  use sll_scalar_field_initializers_base
  implicit none

  type, extends(scalar_field_2d_initializer_base) :: init_gaussian_2d
     class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64 :: xc, yc
    sll_real64 :: sigma_x, sigma_y
  contains
    procedure, pass(init_obj) :: initialize => initialize_gaussian_2d
    procedure, pass(init_obj) :: f_of_x1x2  => f_x1x2_gaussian_2d
 end type init_gaussian_2d

contains

  subroutine initialize_gaussian_2d( init_obj, transf, data_position, xc, yc, sigma_x, sigma_y )
    class(init_gaussian_2d), intent(inout)  :: init_obj
    class(sll_coordinate_transformation_2d_base), target :: transf
    sll_int32 :: data_position
    sll_real64, intent(in), optional     :: xc, yc, sigma_x, sigma_y 
    init_obj%data_position = data_position

    init_obj%transf => transf
    if( present(xc) ) then
       init_obj%xc = xc
    else
       init_obj%xc = 0.0_f64 ! default center is 0
    end if
    if( present(yc) ) then
       init_obj%yc = yc
    else
       init_obj%yc = 0.0_f64 ! default center is 0
    end if
    if( present(sigma_x) ) then
       init_obj%sigma_x = sigma_x
    else
       init_obj%sigma_x = 1.0_f64 ! default radius is 1
    end if
    if( present(sigma_y) ) then
       init_obj%sigma_y = sigma_y
    else
       init_obj%sigma_y = 1.0_f64 ! default radius is 1
    end if
  end subroutine 

  subroutine f_x1x2_gaussian_2d( init_obj, data_out )
    class(init_gaussian_2d), intent(inout)       :: init_obj
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64, dimension(:,:), intent(out)    :: data_out
    type(sll_logical_mesh_2d), pointer         :: mesh
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: num_pts1
    sll_int32  :: num_pts2
    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: jac

    mesh => init_obj%transf%get_logical_mesh()

    if (init_obj%data_position == CELL_CENTERED_FIELD) then
       num_pts1 = mesh%num_cells1
       num_pts2 = mesh%num_cells2
    else if (init_obj%data_position == NODE_CENTERED_FIELD) then
       num_pts1 = mesh%num_cells1 + 1
       num_pts2 = mesh%num_cells2 + 1
    else
       print*, 'sll_gaussian_2d_initializer: Pb with data_position', init_obj%data_position
    end if
    SLL_ASSERT( size(data_out,1) .ge. num_pts1 )
    SLL_ASSERT( size(data_out,2) .ge. num_pts2 )
    if (init_obj%data_position == NODE_CENTERED_FIELD) then
       do j=1,num_pts2
          do i=1, num_pts1
             y = transf%x2_at_node(i,j)
             x = transf%x1_at_node(i,j)
             data_out(i,j) = &
                  1.0_f64/(2*sll_pi*init_obj%sigma_x*init_obj%sigma_y)*exp(-0.5_f64*( &
                  (x-init_obj%xc)**2/init_obj%sigma_x**2 + &
                  (y-init_obj%yc)**2/init_obj%sigma_y**2))
          end do
       end do
    else if (init_obj%data_position == CELL_CENTERED_FIELD) then
       ! jacobian times distribution function is stored
       do j=1,num_pts2
          do i=1, num_pts1
             y = transf%x2_at_cell(i,j)
             x = transf%x1_at_cell(i,j)
             jac = transf%jacobian_at_cell(i,j)
             data_out(i,j) = &
                  jac / (2*sll_pi*init_obj%sigma_x*init_obj%sigma_y)*exp(-0.5_f64*( &
                  (x-init_obj%xc)**2/init_obj%sigma_x**2 + &
                  (y-init_obj%yc)**2/init_obj%sigma_y**2))
          end do
       end do
    end if
  end subroutine

end module sll_gaussian_2d_initializer
