!> @ingroup fem_spline
!> @brief
!> Invert a circulant matrix based on diagonalization in Fourier space (3d version)
!> @details
!> 
!> @authors
!> Katharina Kormann
!>

module sll_m_linear_solver_spline_mass_fft
#include "sll_working_precision.h"
#include "sll_errors.h"

   use sll_m_linear_solver_abstract, only : &
       sll_t_linear_solver_abstract

   use sll_m_fft

  implicit none
  private

  public :: sll_t_linear_solver_spline_mass_fft
  

!> Linear solver for FFT-based inversion of 3d tensor product of circulant matrices  (extending the abstract linear solver class)
  type, extends(sll_t_linear_solver_abstract) :: sll_t_linear_solver_spline_mass_fft

     sll_int32 :: n_dofs(3)   !< no of dofs per dimension

     sll_real64 :: factor = 1.0_f64 !< factor to multiply the matrix with

     sll_real64, allocatable :: inv_eig_values_1(:)  !< eigenvalues of inverse matrix along dimension 1 
     sll_real64, allocatable :: inv_eig_values_2(:)  !< eigenvalues of inverse matrix along dimension 2
     sll_real64, allocatable :: inv_eig_values_3(:)  !< eigenvalues of inverse matrix along dimension 3

     type(sll_t_fft) :: fft1   !< data type for fft along dimension 1
     type(sll_t_fft) :: fft2   !< data type for fft along dimension 2
     type(sll_t_fft) :: fft3   !< data type for fft along dimension 3
     type(sll_t_fft) :: ifft1  !< data type for inverse fft along dimension 1
     type(sll_t_fft) :: ifft2  !< data type for inverse fft along dimension 2
     type(sll_t_fft) :: ifft3  !< data type for inverse fft along dimension 3

   contains

     procedure :: read_from_file => read_from_file_mass1
     procedure :: set_verbose => set_verbose_mass1
     procedure :: solve_real => solve_real_mass1
     procedure :: print_info => print_info_mass1
     procedure :: create => create_mass1
     procedure :: free => free_mass1
     
  end type sll_t_linear_solver_spline_mass_fft

contains

  subroutine create_mass1( self, n_dofs, inv_eig_values_1, inv_eig_values_2, inv_eig_values_3)
    class( sll_t_linear_solver_spline_mass_fft), intent( inout ) :: self !< Fourier solver 
    sll_int32, intent( in ) :: n_dofs(3) !< no of dofs per dimension
    sll_real64, intent( in ) :: inv_eig_values_1(:) !< eigenvalues of inverse matrix along dimension 1 
    sll_real64, intent( in ) :: inv_eig_values_2(:) !< eigenvalues of inverse matrix along dimension 2
    sll_real64, intent( in ) :: inv_eig_values_3(:) !< eigenvalues of inverse matrix along dimension 3
    !local variables
    sll_comp64 :: array1d_x(n_dofs(1)),  array1d_y(n_dofs(2)),  array1d_z(n_dofs(3))

    self%n_dofs = n_dofs
    self%n_global_rows = product( n_dofs)
    self%n_global_cols = self%n_global_rows
    
    allocate(self%inv_eig_values_1(n_dofs(1)))
    self%inv_eig_values_1 = inv_eig_values_1
    allocate(self%inv_eig_values_2(n_dofs(2)))
    self%inv_eig_values_2 = inv_eig_values_2
    allocate(self%inv_eig_values_3(n_dofs(3)))
    self%inv_eig_values_3 = inv_eig_values_3

    call sll_s_fft_init_c2c_1d( self%fft1, n_dofs(1), array1d_x, array1d_x, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft2, n_dofs(2), array1d_y, array1d_y, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft3, n_dofs(3), array1d_z, array1d_z, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%ifft1, n_dofs(1), array1d_x, array1d_x, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft2, n_dofs(2), array1d_y, array1d_y, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft3, n_dofs(3), array1d_z, array1d_z, &
         sll_p_fft_backward, normalized=.true.)
    
  end subroutine create_mass1

  subroutine solve_real_mass1(self, rhs, unknown)
    class( sll_t_linear_solver_spline_mass_fft), intent( inOUT ) :: self !< Fourier solver
    real(kind=f64), intent(in    ) :: rhs(:) !< given right-hand side
    real(kind=f64), intent(  out ) :: unknown(:) !< unknown-left hand side

    sll_int32 :: ind, i, j, k, n_dofs(3)
    sll_comp64 :: scratch(self%n_dofs(1), self%n_dofs(2), self%n_dofs(3))
    sll_comp64 :: array1d_x(self%n_dofs(1)),  array1d_y(self%n_dofs(2)),  array1d_z(self%n_dofs(3))
    
    n_dofs = self%n_dofs
    
    ! Fourier transform
    ind=0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             ind = ind+1
             array1d_x(i) = cmplx( rhs(ind), 0_f64, kind=f64)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft1, array1d_x, array1d_x)
          do i=1,n_dofs(1)
             scratch(i,j,k) = array1d_x(i)
          end do
       end do
    end do
    do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             array1d_y(j) = scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft2, array1d_y, array1d_y)
          do j=1,n_dofs(2)
             scratch(i,j,k) = array1d_y(j)
          end do
       end do
    end do
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          do k=1,n_dofs(3)
             array1d_z(k) = scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft3, array1d_z, array1d_z)
          do k=1,n_dofs(3)
             scratch(i,j,k) = array1d_z(k)
          end do
       end do
    end do

    ! Multiply by inverse mass
     do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)

             scratch(i,j,k) = scratch(i,j,k)* &
                  cmplx(self%inv_eig_values_1(i)* &
                  self%inv_eig_values_2(j)* &
                  self%inv_eig_values_3(k) /self%factor, 0.0_f64, f64)

          end do
       end do
    end do
    

    ! Inverse Fourier transform
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          do k=1,n_dofs(3)
             array1d_z(k) = scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft3, array1d_z, array1d_z)
          do k=1,n_dofs(3)
             scratch(i,j,k) = array1d_z(k)
          end do
       end do
    end do

    do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             array1d_y(j) = scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft2, array1d_y, array1d_y)
          do j=1,n_dofs(2)
             scratch(i,j,k) = array1d_y(j)
          end do
       end do
    end do
    
    ind=0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             array1d_x(i) = scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft1, array1d_x, array1d_x)
          
          do i=1,n_dofs(1)
             ind = ind+1
             unknown(ind) = real( array1d_x(i), kind=f64 )
          end do
       end do
    end do
    
  end subroutine solve_real_mass1
    
  subroutine read_from_file_mass1(self, filename)
    class( sll_t_linear_solver_spline_mass_fft), intent( inout ) :: self !< Fourier solver
    character(len=*), intent( in ) :: filename !< filename
    
  end subroutine read_from_file_mass1
    
  subroutine print_info_mass1(self)
    class( sll_t_linear_solver_spline_mass_fft), intent( in ) :: self !< Fourier solver
    
  end subroutine print_info_mass1

  subroutine set_verbose_mass1( self, verbose )
    class( sll_t_linear_solver_spline_mass_fft), intent( inout ) :: self !< Fourier solver
    logical, intent( in ) :: verbose !< logical for convergence information

    self%verbose = verbose
    
  end subroutine set_verbose_mass1

  
  subroutine free_mass1(self)
    class( sll_t_linear_solver_spline_mass_fft), intent( inout ) :: self !< Fourier solver
    
  end subroutine free_mass1

end module sll_m_linear_solver_spline_mass_fft
