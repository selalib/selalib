!> \file poly_2d_initializer.F90
!> \authors                    
!> Martin CAMPOS PINTO (campos@ann.jussieu.fr) 
!> \brief  
!> polynomial initializer, used to test the loading of ltpic-particles
!>

module sll_poly_2d_initializer
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_module_mapped_meshes_2d_base
  use sll_scalar_field_initializers_base
  implicit none

  type, extends(scalar_field_2d_initializer_base) :: init_poly_2d
    sll_int32                               :: degree
    sll_real64, dimension(:,:), pointer     :: coefs_xv
    class(sll_mapped_mesh_2d_base), pointer :: mesh
  contains
    procedure, pass(init_obj) :: initialize => initialize_poly_2d
    procedure, pass(init_obj) :: f_of_x1x2  => f_x1x2_poly_2d
  end type init_poly_2d



contains

  subroutine initialize_poly_2d( init_obj, mesh, data_position, degree_val, coefs_xv_val )
    class(init_poly_2d), intent(inout)                  :: init_obj
    class(sll_mapped_mesh_2d_base), intent(in), target  :: mesh
    sll_int32,  intent(in)                              :: data_position
    sll_int32,  intent(in)                              :: degree_val
    sll_real64, dimension(:,:), intent(in)              :: coefs_xv_val
    
    sll_int32                                           :: ierr

    SLL_ASSERT( size(coefs_xv_val,1) .eq. degree_val+1 )
    SLL_ASSERT( size(coefs_xv_val,2) .eq. degree_val+1 )
    SLL_ALLOCATE( init_obj%coefs_xv(degree_val+1,degree_val+1),   ierr )    
    init_obj%mesh => mesh
    init_obj%data_position = data_position
    init_obj%degree        = degree_val
    init_obj%coefs_xv      = coefs_xv_val
  end subroutine initialize_poly_2d

  subroutine f_x1x2_poly_2d( init_obj, data_out )
    class(init_poly_2d),        intent(inout)  :: init_obj
    sll_real64, dimension(:,:), intent(out)    :: data_out

    class(sll_mapped_mesh_2d_base), pointer    :: mesh
    sll_real64, dimension(:,:),     pointer    :: coefs_xv
    sll_int32     :: degree   
    sll_int32     :: i
    sll_int32     :: j
    sll_int32     :: num_pts1
    sll_int32     :: num_pts2
    sll_real64    :: x
    sll_real64    :: v
    sll_int32     :: a
    sll_int32     :: b
    sll_real64    :: x_pow_a
    sll_real64    :: xv_pow_ab
        
    mesh     => init_obj%mesh
    coefs_xv => init_obj%coefs_xv
    degree   =  init_obj%degree
    
    if (init_obj%data_position ==  NODE_CENTERED_FIELD) then
       num_pts1 = mesh%nc_eta1+1
       num_pts2 = mesh%nc_eta2+1
    else if (init_obj%data_position ==  NODE_CENTERED_FIELD) then
       num_pts1 = mesh%nc_eta1
       num_pts2 = mesh%nc_eta2
    end if
    SLL_ASSERT( size(data_out,1) .ge. num_pts1 )
    SLL_ASSERT( size(data_out,2) .ge. num_pts2 )
    do j=1,num_pts2
       do i=1, num_pts1
          if (init_obj%data_position ==  NODE_CENTERED_FIELD) then
             v = mesh%x2_at_node(i,j)
             x = mesh%x1_at_node(i,j)
          else if (init_obj%data_position ==  NODE_CENTERED_FIELD) then
             v = mesh%x2_at_cell(i,j)
             x = mesh%x1_at_cell(i,j)
          else
             print*, 'f_x1x2_tsi_2d:',  init_obj%data_position, 'not defined'
          end if
          x_pow_a = 1
          do a = 0, degree
             xv_pow_ab = x_pow_a 
             do b = 0, degree
                data_out(i,j) = data_out(i,j) + coefs_xv(a+1,b+1) * xv_pow_ab
                xv_pow_ab = xv_pow_ab * v 
             end do
             x_pow_a = x_pow_a * x
          end do
       end do
    end do
  end subroutine f_x1x2_poly_2d

end module sll_poly_2d_initializer
