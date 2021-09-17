!> @ingroup fem_spline
!> @brief
!> Invert a circulant matrix based on diagonalization in Fourier space (2d version)
!> @details
!> 
!> @authors
!> Katharina Kormann
!>
module sll_m_linear_solver_spline_mass_2d_fft
#include "sll_working_precision.h"
#include "sll_errors.h"

   use sll_m_linear_solver_abstract, only : &
       sll_t_linear_solver_abstract

   use sll_m_fft

  implicit none
  private

  public :: sll_t_linear_solver_spline_mass_2d_fft
  

!> Data type for a linear solver inverting a 2d tensor product of circulant matrices based on FFT
  type, extends(sll_t_linear_solver_abstract) :: sll_t_linear_solver_spline_mass_2d_fft

     sll_int32 :: n_dofs(2)  !< no of dofs per dimension

     sll_real64 :: factor = 1.0_f64  !< factor to multiply the matrix with

     sll_real64, allocatable :: eig_values_1(:) !< eigenvalues of the matrix to be inverted along dimension 1
     sll_real64, allocatable :: eig_values_2(:) !< eigenvalues of the matrix to be inverted along dimension 2

     type(sll_t_fft) :: fft1   !< data type for FFT along dimension 1
     type(sll_t_fft) :: fft2   !< data type for FFT along dimension 2
     type(sll_t_fft) :: ifft1  !< data type for inverse FFT along dimension 1
     type(sll_t_fft) :: ifft2  !< data type for inverse FFT along dimension 2

     
   contains

     procedure :: read_from_file => read_from_file_mass1
     procedure :: set_verbose => set_verbose_mass1
     procedure :: solve_real => solve_real_mass1
     procedure :: solve_complex => solve_complex_mass1
     procedure :: print_info => print_info_mass1
     procedure :: create => create_mass1
     procedure :: free => free_mass1
     
  end type sll_t_linear_solver_spline_mass_2d_fft

contains

  subroutine create_mass1( self, n_dofs, eig_values_1, eig_values_2)
    class( sll_t_linear_solver_spline_mass_2d_fft), intent( inout ) :: self !< Fourier solver 
    sll_int32, intent( in ) :: n_dofs(2) !< no of dofs per dimension
    sll_real64, intent( in ) :: eig_values_1(:) !< eigenvalues of inverse matrix along dimension 1 
    sll_real64, intent( in ) :: eig_values_2(:) !< eigenvalues of inverse matrix along dimension 2
    !local variables
    sll_comp64 :: array1d_x(n_dofs(1)),  array1d_y(n_dofs(2))

    self%n_dofs = n_dofs

    allocate(self%eig_values_1(n_dofs(1)))
    self%eig_values_1 = eig_values_1
    allocate(self%eig_values_2(n_dofs(2)))
    self%eig_values_2 = eig_values_2

    call sll_s_fft_init_c2c_1d( self%fft1, n_dofs(1), array1d_x, array1d_x, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft2, n_dofs(2), array1d_y, array1d_y, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%ifft1, n_dofs(1), array1d_x, array1d_x, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft2, n_dofs(2), array1d_y, array1d_y, &
         sll_p_fft_backward, normalized=.true.)
    
  end subroutine create_mass1

  subroutine solve_real_mass1(self, rhs, unknown)
    class( sll_t_linear_solver_spline_mass_2d_fft), intent( inout ) :: self !< Fourier solver
    real(kind=f64), intent(in   ) :: rhs(:) !< given right-hand side
    real(kind=f64), intent(  out) :: unknown(:) !< unknown-left hand side

    sll_int32 :: ind, i, j, n_dofs(2)
    sll_comp64 :: scratch(self%n_dofs(1), self%n_dofs(2))
    sll_comp64 :: array1d_x(self%n_dofs(1)),  array1d_y(self%n_dofs(2))
    
    n_dofs = self%n_dofs
    
    ! Fourier transform
    ind=0
    !do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             ind = ind+1
             array1d_x(i) = cmplx( rhs(ind), 0_f64, kind=f64)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft1, array1d_x, array1d_x)
          do i=1,n_dofs(1)
             scratch(i,j) = array1d_x(i)
          end do
       end do
    !end do
    !do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             array1d_y(j) = scratch(i,j)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft2, array1d_y, array1d_y)
          do j=1,n_dofs(2)
             scratch(i,j) = array1d_y(j)
          end do
       end do
    !end do

    ! Multiply by inverse mass
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)

             scratch(i,j) = scratch(i,j)/ &
                  cmplx(self%eig_values_1(i)* &
                  self%eig_values_2(j)* &
                  self%factor, 0.0_f64, f64)

          end do
       end do
    

    ! Inverse Fourier transform

    !do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             array1d_y(j) = scratch(i,j)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft2, array1d_y, array1d_y)
          do j=1,n_dofs(2)
             scratch(i,j) = array1d_y(j)
          end do
       end do
    !end do
    
    ind=0
    !do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             array1d_x(i) = scratch(i,j)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft1, array1d_x, array1d_x)
          
          do i=1,n_dofs(1)
             ind = ind+1
             unknown(ind) = real( array1d_x(i), kind=f64 )
          end do
       end do
    !end do
    
  end subroutine solve_real_mass1


  
  
  subroutine solve_complex_mass1(self, rhs, unknown)
    class( sll_t_linear_solver_spline_mass_2d_fft), intent( in ) :: self  !< Fourier solver
    complex(kind=f64), intent(in   ) :: rhs(:) !< given right-hand side
    complex(kind=f64), intent(  out) :: unknown(:) !< unknown-left hand side

    SLL_ERROR('solve_complex', 'Procedure not implemented.')
    
  end subroutine solve_complex_mass1
    
  subroutine read_from_file_mass1(self, filename)
    class( sll_t_linear_solver_spline_mass_2d_fft), intent( inout ) :: self  !< Fourier solver
    character(len=*), intent( in ) :: filename !< filename
    
  end subroutine read_from_file_mass1
    
  subroutine print_info_mass1(self)
    class( sll_t_linear_solver_spline_mass_2d_fft), intent( in ) :: self  !< Fourier solver
    
  end subroutine print_info_mass1

  subroutine set_verbose_mass1( self, verbose )
    class( sll_t_linear_solver_spline_mass_2d_fft), intent( inout ) :: self  !< Fourier solver
    logical, intent( in ) :: verbose !< logical for convergence information

    self%verbose = verbose
    
  end subroutine set_verbose_mass1

  
  subroutine free_mass1(self)
    class( sll_t_linear_solver_spline_mass_2d_fft), intent( inout ) :: self  !< Fourier solver
    
  end subroutine free_mass1

end module sll_m_linear_solver_spline_mass_2d_fft
