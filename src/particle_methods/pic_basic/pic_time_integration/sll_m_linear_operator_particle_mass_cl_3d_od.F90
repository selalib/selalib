module sll_m_linear_operator_particle_mass_cl_3d_od
#include "sll_working_precision.h"
  
  use sll_m_particle_mass_3d_base, only: &
       sll_c_particle_mass_3d_base

  
  implicit none

  public :: sll_t_linear_operator_particle_mass_cl_3d_od
  
  private

  type, extends(sll_c_particle_mass_3d_base) :: sll_t_linear_operator_particle_mass_cl_3d_od

     sll_int32 :: degree2(3)  !< spline degree2
     sll_int32 :: n_dofs2(3) !< number of dofs2
     sll_int32 :: n_total2 !< product of dofs
         

   contains
     procedure :: create => create_linear_operator_particle_mass_cl_3d_od
     procedure :: free => free_particle_mass_cl_3d_od
     procedure :: dot => dot_particle_mass_cl_3d_od
     procedure :: print_info => print_info_particle_mass_cl_3d_od

  end type sll_t_linear_operator_particle_mass_cl_3d_od


contains

  subroutine create_linear_operator_particle_mass_cl_3d_od( self, n_cells, degree, degree2 )
    class(sll_t_linear_operator_particle_mass_cl_3d_od), intent( inout ) :: self !< Particle mass
    sll_int32,  intent( in ) :: n_cells(3) !< grid cells
    sll_int32,  intent( in ) :: degree(3) !< spline degree1
    sll_int32, optional,  intent( in ) :: degree2(3) !< spline degree2
  
    self%degree = degree
    self%n_dofs  = n_cells
    self%n_dofs(1)  = n_cells(1)+degree(1)
    self%n_total = product(self%n_dofs)
    self%degree2 = degree2
    self%n_dofs2  = n_cells
    self%n_dofs2(1)  = n_cells(1)+degree2(1)
    self%n_total2 = product(self%n_dofs2)
     
    allocate( self%particle_mass( product(self%degree+self%degree2+1), self%n_total))
    self%particle_mass = 0._f64

    self%n_rows = self%n_total
    self%n_cols = self%n_total2
    
    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
    
  end subroutine create_linear_operator_particle_mass_cl_3d_od

  subroutine free_particle_mass_cl_3d_od( self )
    class(sll_t_linear_operator_particle_mass_cl_3d_od), intent( inout ) :: self !< Particle mass

    deallocate( self%particle_mass )

  end subroutine free_particle_mass_cl_3d_od
  
  
  subroutine dot_particle_mass_cl_3d_od ( self, x, y )
    class(sll_t_linear_operator_particle_mass_cl_3d_od), intent( in ) :: self !< Particle mass
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outpoutvariable
    !local variables
    sll_int32 :: ind_i, ind_j, ind
    sll_int32 :: i1, i2, i3, j1, j2, j3
    
    ind_i = 1
    do i3 = 1, self%n_dofs(3)
       do i2 = 1, self%n_dofs(2)
           do i1 = 1, self%degree(1)
             y(ind_i) = 0.0_f64
             ind = 1
             do j3 = i3-self%degree2(3), i3+self%degree(3)
                do j2 = i2-self%degree2(2), i2+self%degree(2)
                   ind = ind + self%degree(1)+1-i1 
                   do j1 = 1, i1+self%degree2(1)
                      ind_j = j1 + modulo(j2-1,self%n_dofs2(2))*self%n_dofs2(1) + modulo(j3-1,self%n_dofs2(3))*self%n_dofs2(1)*self%n_dofs2(2)
                      y(ind_i) = y(ind_i) + &
                           self%sign * self%particle_mass( ind, ind_i) * x(ind_j)
                      ind = ind+1
                   end do
                end do
             end do
             ind_i = ind_i + 1 
          end do
          do i1 = self%degree(1)+1, self%n_dofs(1)-self%degree(1)
             y(ind_i) = 0.0_f64
             ind = 1
             do j3 = i3-self%degree2(3), i3+self%degree(3)
                do j2 = i2-self%degree2(2), i2+self%degree(2)
                   do j1 = i1-self%degree(1), i1+self%degree2(1)
                      ind_j = j1 + modulo(j2-1,self%n_dofs2(2))*self%n_dofs2(1) + modulo(j3-1,self%n_dofs2(3))*self%n_dofs2(1)*self%n_dofs2(2)
                      y(ind_i) = y(ind_i) + &
                           self%sign * self%particle_mass( ind, ind_i) * x(ind_j)
                      ind = ind+1
                   end do
                end do
             end do
             ind_i = ind_i + 1 
          end do
           do i1 = self%n_dofs(1)-self%degree(1)+1, self%n_dofs(1)
             y(ind_i) = 0.0_f64
             ind = 1
             do j3 = i3-self%degree2(3), i3+self%degree(3)
                do j2 = i2-self%degree2(2), i2+self%degree(2)
                   do j1 = i1-self%degree(1), self%n_dofs2(1)
                      ind_j = j1 + modulo(j2-1,self%n_dofs2(2))*self%n_dofs2(1) + modulo(j3-1,self%n_dofs2(3))*self%n_dofs2(1)*self%n_dofs2(2)
                      y(ind_i) = y(ind_i) + &
                           self%sign * self%particle_mass( ind, ind_i) * x(ind_j)
                      ind = ind+1
                   end do
                   ind = ind + self%degree(1)-self%n_dofs(1)+i1 
                end do
             end do
             ind_i = ind_i + 1 
          end do
       end do
    end do

    
  end subroutine dot_particle_mass_cl_3d_od


  subroutine print_info_particle_mass_cl_3d_od( self )
    class(sll_t_linear_operator_particle_mass_cl_3d_od), intent(in) :: self !< Particle mass
  end subroutine print_info_particle_mass_cl_3d_od
  
end module sll_m_linear_operator_particle_mass_cl_3d_od
