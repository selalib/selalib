module sll_m_linear_operator_particle_mass_3d_diag
#include "sll_working_precision.h"

  use sll_m_particle_mass_3d_base, only: &
       sll_c_particle_mass_3d_base

  implicit none

  public :: sll_t_linear_operator_particle_mass_3d_diag
  
  private

  type, extends(sll_c_particle_mass_3d_base):: sll_t_linear_operator_particle_mass_3d_diag

   contains
     procedure :: create => create_linear_operator_particle_mass_3d_diag
     procedure :: free => free_particle_mass_3d_diag
     procedure :: dot => dot_particle_mass_3d_diag
     procedure :: print_info => print_info_particle_mass_3d_diag

  end type sll_t_linear_operator_particle_mass_3d_diag


contains

  subroutine create_linear_operator_particle_mass_3d_diag( self, n_cells, degree, degree2 )
    class(sll_t_linear_operator_particle_mass_3d_diag), intent( inout ) :: self !< Particle mass
    sll_int32,  intent( in ) :: n_cells(3) !< grid cells
    sll_int32,  intent( in ) :: degree(3) !< spline degree1
    sll_int32, optional,  intent( in ) :: degree2(3) !< spline degree2
    
    self%degree = degree
    self%n_dofs = n_cells
    self%n_total = product(n_cells)
 
    allocate( self%particle_mass(  product(2*self%degree+1), self%n_total) )
    self%particle_mass = 0._f64
    self%n_rows = self%n_total
    self%n_cols = self%n_total
    
    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
    
  end subroutine create_linear_operator_particle_mass_3d_diag

  subroutine free_particle_mass_3d_diag( self )
    class(sll_t_linear_operator_particle_mass_3d_diag), intent( inout ) :: self !< Particle mass

    deallocate( self%particle_mass )

  end subroutine free_particle_mass_3d_diag
  
  
  subroutine dot_particle_mass_3d_diag ( self, x, y )
    class(sll_t_linear_operator_particle_mass_3d_diag), intent( in ) :: self !< Particle mass
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outpoutvariable
    !local variables
    sll_int32 :: ind_i, ind_j, ind
    sll_int32 :: i1, i2, i3, j1, j2, j3
    
    ind_i = 1
    do i3 = 1, self%n_dofs(3)
       do i2 = 1, self%n_dofs(2)
          do i1 = 1, self%n_dofs(1)
             y(ind_i) = 0.0_f64
             ind = 1
             do j3 = i3-self%degree(3), i3+self%degree(3)
                do j2 = i2-self%degree(2), i2+self%degree(2)
                   do j1 = i1-self%degree(1), i1+self%degree(1)
                      ind_j = index_3dto1d( self%n_dofs, j1, j2, j3)
                      y(ind_i) = y(ind_i) + &
                           self%sign * self%particle_mass( ind, ind_i) * x(ind_j)
                      ind = ind+1
                   end do
                end do
             end do
             ind_i = ind_i + 1 
          end do
       end do
    end do

    
  end subroutine dot_particle_mass_3d_diag


  function index_3dto1d( num_pts, ind1, ind2, ind3 ) result( ind1d)
    sll_int32 :: num_pts(3)
    sll_int32 :: ind1
    sll_int32 :: ind2
    sll_int32 :: ind3
    sll_int32 :: ind1d

    ind1d = 1+ modulo(ind1-1,num_pts(1)) + &
         modulo(ind2-1,num_pts(2)) * num_pts(1) + &
         modulo(ind3-1,num_pts(3))* num_pts(1)*num_pts(2)


  end function index_3dto1d


  subroutine print_info_particle_mass_3d_diag( self )
    class(sll_t_linear_operator_particle_mass_3d_diag), intent(in) :: self !< Particle mass 
  end subroutine print_info_particle_mass_3d_diag
  
end module sll_m_linear_operator_particle_mass_3d_diag
