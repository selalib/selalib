!> @ingroup fem_spline
!> @brief
!> Invert a circulant matrix based on diagonalization in Fourier space (3d version)
!> @details
!> 
!> @authors
!> Benedikt Perse
!>

module sll_m_preconditioner_curl_solver_fft
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_linear_solver_abstract, only : &
       sll_t_linear_solver_abstract

  use sll_m_constants

  use sll_m_fft

  use sll_m_spline_fem_utilities

  implicit none
  private

  public :: sll_t_preconditioner_curl_solver_fft


  !> Linear solver for FFT-based inversion of 3d tensor product of circulant matrices  (extending the abstract linear solver class)
  type, extends(sll_t_linear_solver_abstract) :: sll_t_preconditioner_curl_solver_fft
     sll_int32 :: n_total  !< product of number of degrees of freedom
     sll_int32 :: n_dofs(3) !< number of degrees of freedom

     ! For Fourier variant
     sll_comp64, allocatable :: eig_values_d1(:)
     sll_comp64, allocatable :: eig_values_d1t(:)
     sll_comp64, allocatable :: eig_values_d2(:)
     sll_comp64, allocatable :: eig_values_d2t(:)
     sll_comp64, allocatable :: eig_values_d3(:)
     sll_comp64, allocatable :: eig_values_d3t(:)
     sll_real64, allocatable :: eig_values_mass_0_1(:)
     sll_real64, allocatable :: eig_values_mass_0_2(:)
     sll_real64, allocatable :: eig_values_mass_0_3(:)
     sll_real64, allocatable :: eig_values_mass_1_1(:)
     sll_real64, allocatable :: eig_values_mass_1_2(:)
     sll_real64, allocatable :: eig_values_mass_1_3(:)


     type(sll_t_fft) :: fft1   !< data type for fft along dimension 1
     type(sll_t_fft) :: fft2   !< data type for fft along dimension 2
     type(sll_t_fft) :: fft3   !< data type for fft along dimension 3
     type(sll_t_fft) :: ifft1  !< data type for inverse fft along dimension 1
     type(sll_t_fft) :: ifft2  !< data type for inverse fft along dimension 2
     type(sll_t_fft) :: ifft3  !< data type for inverse fft along dimension 3

   contains

     procedure :: read_from_file => read_from_file_preconditioner
     procedure :: set_verbose => set_verbose_preconditioner
     procedure :: solve_real => solve_real_preconditioner
     procedure :: print_info => print_info_preconditioner
     procedure :: create => create_preconditioner
     procedure :: free => free_preconditioner

  end type sll_t_preconditioner_curl_solver_fft

contains

  subroutine create_preconditioner( self, n_dofs, delta_x, degree )
    class( sll_t_preconditioner_curl_solver_fft), intent( inout ) :: self !< Fourier solver 
    sll_int32  :: n_dofs(3) !< number of degrees of freedom
    sll_real64 :: delta_x(3) !< cell size
    sll_int32  :: degree(3)
    !local variables
    sll_real64 :: angle
    sll_int32 :: j
    sll_comp64 :: array1d_x(n_dofs(1))
    sll_comp64 :: array1d_y(n_dofs(2))
    sll_comp64 :: array1d_z(n_dofs(3))
    sll_real64, allocatable :: mass_line_0(:)
    sll_real64, allocatable :: mass_line_1(:)


    allocate( self%eig_values_mass_0_1(n_dofs(1)) )
    allocate( self%eig_values_mass_0_2(n_dofs(2)) )
    allocate( self%eig_values_mass_0_3(n_dofs(3)) )
    allocate( self%eig_values_mass_1_1(n_dofs(1)) )
    allocate( self%eig_values_mass_1_2(n_dofs(2)) )
    allocate( self%eig_values_mass_1_3(n_dofs(3)) )


    self%n_dofs = n_dofs
    self%n_total = product(self%n_dofs)

    self%n_rows = self%n_total*3
    self%n_cols = self%n_total*3

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols

    ! Fourier product

    ! Eigenvalues of derivative matrices
    allocate( self%eig_values_d1(1:self%n_dofs(1)) )
    allocate( self%eig_values_d1t(1:self%n_dofs(1) ) )

    self%eig_values_d1(1) = cmplx(0.0_f64, 0.0_f64, f64)
    self%eig_values_d1t(1) = cmplx(0.0_f64, 0.0_f64, f64)
    do j=2,self%n_dofs(1)
       angle = sll_p_twopi*real(j-1,f64)/real(self%n_dofs(1), f64)

       self%eig_values_d1(j) = cmplx((1.0_f64 - cos(angle))/delta_x(1),sin(angle)/delta_x(1), f64 )
       self%eig_values_d1t(j) = cmplx((1.0_f64 - cos(angle))/delta_x(1),-sin(angle)/delta_x(1), f64 )
    end do


    allocate( self%eig_values_d2(1:self%n_dofs(2)) )
    allocate( self%eig_values_d2t(1:self%n_dofs(2) ) )

    self%eig_values_d2(1) =  cmplx(0.0_f64, 0.0_f64, f64)
    self%eig_values_d2t(1) =  cmplx(0.0_f64, 0.0_f64, f64)
    do j=2,self%n_dofs(2)
       angle = sll_p_twopi*real(j-1,f64)/real(self%n_dofs(2), f64)

       self%eig_values_d2(j) = cmplx((1.0_f64 - cos(angle))/delta_x(2),sin(angle)/delta_x(2), f64 )
       self%eig_values_d2t(j) = cmplx((1.0_f64 - cos(angle))/delta_x(2),-sin(angle)/delta_x(2), f64 )
    end do

    allocate( self%eig_values_d3(1:self%n_dofs(3)) )
    allocate( self%eig_values_d3t(1:self%n_dofs(3) ) )

    self%eig_values_d3(1) = cmplx(0.0_f64, 0.0_f64, f64)
    self%eig_values_d3t(1) = cmplx(0.0_f64, 0.0_f64, f64)
    do j=2,self%n_dofs(3)
       angle = sll_p_twopi*real(j-1,f64)/real(self%n_dofs(3), f64)

       self%eig_values_d3(j) = cmplx((1.0_f64 - cos(angle))/delta_x(3),sin(angle)/delta_x(3), f64 )
       self%eig_values_d3t(j) = cmplx((1.0_f64 - cos(angle))/delta_x(3),-sin(angle)/delta_x(3), f64 )
    end do

    allocate( mass_line_0(degree(1)+1) )
    allocate( mass_line_1(degree(1)) )
    ! Eigenvalues of mass matrices
    call sll_s_spline_fem_mass_line( degree(1), mass_line_0 )
    call sll_s_spline_fem_mass_line( degree(1)-1, mass_line_1 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs(1), degree(1),  mass_line_0*delta_x(1), self%eig_values_mass_0_1 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs(1), degree(1)-1, mass_line_1*delta_x(1), self%eig_values_mass_1_1 )
    deallocate( mass_line_0 )
    deallocate( mass_line_1 )


    allocate( mass_line_0(degree(2)+1) )
    allocate( mass_line_1(degree(2)) )
    ! Eigenvalues of mass matrices
    call sll_s_spline_fem_mass_line( degree(2), mass_line_0 )
    call sll_s_spline_fem_mass_line( degree(2)-1, mass_line_1 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs(2), degree(2),  mass_line_0*delta_x(2), self%eig_values_mass_0_2 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs(2), degree(2)-1, mass_line_1*delta_x(2), self%eig_values_mass_1_2 )
    deallocate( mass_line_0 )
    deallocate( mass_line_1 )

    allocate( mass_line_0(degree(3)+1) )
    allocate( mass_line_1(degree(3)) )
    ! Eigenvalues of mass matrices
    call sll_s_spline_fem_mass_line( degree(3), mass_line_0 )
    call sll_s_spline_fem_mass_line( degree(3)-1, mass_line_1 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs(3), degree(3),  mass_line_0*delta_x(3), self%eig_values_mass_0_3 )
    call sll_s_spline_fem_compute_mass_eig( self%n_dofs(3), degree(3)-1, mass_line_1*delta_x(3), self%eig_values_mass_1_3 )
    deallocate( mass_line_0 )
    deallocate( mass_line_1 )


    ! Initialize fft
    call sll_s_fft_init_c2c_1d( self%fft1, self%n_dofs(1), array1d_x, array1d_x, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft2, self%n_dofs(2), array1d_y, array1d_y, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft3, self%n_dofs(3), array1d_z, array1d_z, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%ifft1, self%n_dofs(1), array1d_x, array1d_x, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft2, self%n_dofs(2), array1d_y, array1d_y, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft3, self%n_dofs(3), array1d_z, array1d_z, &
         sll_p_fft_backward, normalized=.true.)


  end subroutine create_preconditioner

  subroutine solve_real_preconditioner(self, rhs, unknown)
    class( sll_t_preconditioner_curl_solver_fft), intent( inout ) :: self !< Fourier solver
    real(kind=f64), intent(in    ) :: rhs(:) !< given right-hand side
    real(kind=f64), intent(  out ) :: unknown(:) !< unknown-left hand side
    !local variables
    sll_comp64 :: scratch1(self%n_dofs(1),self%n_dofs(2), self%n_dofs(3),3)
    sll_comp64 :: scratch2(self%n_dofs(1),self%n_dofs(2), self%n_dofs(3),3)
    sll_comp64 :: mat(3,3), imat(3,3)
    sll_comp64 :: array1d_x(self%n_dofs(1))
    sll_comp64 :: array1d_y(self%n_dofs(2))
    sll_comp64 :: array1d_z(self%n_dofs(3))
    sll_int32 :: i,j,k



    ! Compute Fourier transform of input
    call fft3d( self, array1d_x, array1d_y, array1d_z, 1, rhs(1:self%n_total), scratch1 )
    call fft3d( self, array1d_x, array1d_y, array1d_z, 2, rhs(1+self%n_total:2*self%n_total), scratch1 )
    call fft3d( self, array1d_x, array1d_y, array1d_z, 3, rhs(1+2*self%n_total:3*self%n_total), scratch1 )


    ! Apply eigenvalues to Fourier coefficient
    do k=1,self%n_dofs(3)
       do j=1,self%n_dofs(2)
          do i=1,self%n_dofs(1)

             mat(1,1) =  (self%eig_values_d3t(k)*self%eig_values_d3(k)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  self%eig_values_d2t(j)*self%eig_values_d2(j)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(1,2) = - ( self%eig_values_d2t(j)*self%eig_values_d1(i)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(1,3) = - ( self%eig_values_d3t(k)*self%eig_values_d1(i)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))

             mat(2,2) = (self%eig_values_d3t(k)*self%eig_values_d3(k)*&
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  self%eig_values_d1t(i)*self%eig_values_d1(i)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(2,1) = - ( self%eig_values_d1t(i)*self%eig_values_d2(j)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(2,3) = - ( self%eig_values_d3t(k)*self%eig_values_d2(j)* &
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))

             mat(3,3) = (self%eig_values_d1t(i)*self%eig_values_d1(i)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  self%eig_values_d2t(j)*self%eig_values_d2(j)*&
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))
             mat(3,2) = - ( self%eig_values_d2t(j)*self%eig_values_d3(k)* &
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))
             mat(3,1) = - ( self%eig_values_d1t(i)*self%eig_values_d3(k)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))

             call invert3d( mat, imat )

             scratch2(i,j,k,1) = imat(1,1)*scratch1(i,j,k,1) + &
                  imat(1,2)*scratch1(i,j,k,2) + &
                  imat(1,3)*scratch1(i,j,k,3)
             scratch2(i,j,k,2) = imat(2,1)*scratch1(i,j,k,1) + &
                  imat(2,2)*scratch1(i,j,k,2) + &
                  imat(2,3)*scratch1(i,j,k,3)
             scratch2(i,j,k,3) = imat(3,1)*scratch1(i,j,k,1) + &
                  imat(3,2)*scratch1(i,j,k,2) + &
                  imat(3,3)*scratch1(i,j,k,3)

          end do
       end do
    end do
    !scratch2 = scratch1

    ! Compute inverse Fourier transform of result
    call ifft3d( self, array1d_x, array1d_y, array1d_z, 1, scratch2, unknown(1:self%n_total) )
    call ifft3d( self, array1d_x, array1d_y, array1d_z, 2, scratch2, unknown(1+self%n_total:2*self%n_total) )
    call ifft3d( self, array1d_x, array1d_y, array1d_z, 3, scratch2, unknown(1+2*self%n_total:3*self%n_total) )

  end subroutine solve_real_preconditioner


  !> Helper function
  subroutine fft3d( self, array1d_x, array1d_y, array1d_z, inde, x, scratch1 )
    class( sll_t_preconditioner_curl_solver_fft), intent( inout ) :: self !< Fourier solver
    sll_comp64, intent( inout ) :: array1d_x(:)
    sll_comp64, intent( inout ) :: array1d_y(:)
    sll_comp64, intent( inout ) :: array1d_z(:)
    sll_int32, intent( in ) :: inde
    sll_real64, intent( in ) :: x(:)
    sll_comp64, intent( out ) :: scratch1(:,:,:,:)
    !local variables
    sll_int32 :: ind, i,j,k


    ind=0
    do k=1,self%n_dofs(3)
       do j=1,self%n_dofs(2)
          do i=1,self%n_dofs(1)
             ind = ind+1
             array1d_x(i) = cmplx( x(ind), 0_f64, kind=f64)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft1, array1d_x, array1d_x)
          do i=1,self%n_dofs(1)
             scratch1(i,j,k,inde) = array1d_x(i)
          end do
       end do
    end do
    do k=1,self%n_dofs(3)
       do i=1,self%n_dofs(1)
          do j=1,self%n_dofs(2)
             array1d_y(j) = scratch1(i,j,k,inde)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft2, array1d_y, array1d_y)
          do j=1,self%n_dofs(2)
             scratch1(i,j,k,inde) = array1d_y(j)
          end do
       end do
    end do
    do j=1,self%n_dofs(2)
       do i=1,self%n_dofs(1)
          do k=1,self%n_dofs(3)
             array1d_z(k) = scratch1(i,j,k,inde)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft3, array1d_z, array1d_z)
          do k=1,self%n_dofs(3)
             scratch1(i,j,k,inde) = array1d_z(k)
          end do
       end do
    end do

  end subroutine fft3d

  !> Helper function
  subroutine ifft3d( self, array1d_x, array1d_y, array1d_z, inde,  scratch, y )
    class( sll_t_preconditioner_curl_solver_fft), intent( inout ) :: self !< Fourier solver
    sll_comp64, intent( inout ) :: array1d_x(:)
    sll_comp64, intent( inout ) :: array1d_y(:)
    sll_comp64, intent( inout ) :: array1d_z(:)
    sll_int32, intent( in ) :: inde
    sll_real64, intent( out ) :: y(:)
    sll_comp64, intent( inout ) :: scratch(:,:,:,:)
    !local variables
    sll_int32 :: ind, i,j,k


    do j=1,self%n_dofs(2)
       do i=1,self%n_dofs(1)
          do k=1,self%n_dofs(3)
             array1d_z(k) = scratch(i,j,k,inde)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft3, array1d_z, array1d_z)
          do k=1,self%n_dofs(3)
             scratch(i,j,k,inde) = array1d_z(k)
          end do
       end do
    end do

    do k=1,self%n_dofs(3)
       do i=1,self%n_dofs(1)
          do j=1,self%n_dofs(2)
             array1d_y(j) = scratch(i,j,k,inde)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft2, array1d_y, array1d_y)
          do j=1,self%n_dofs(2)
             scratch(i,j,k,inde) = array1d_y(j)
          end do
       end do
    end do

    ind=0
    do k=1,self%n_dofs(3)
       do j=1,self%n_dofs(2)
          do i=1,self%n_dofs(1)
             array1d_x(i) = scratch(i,j,k,inde)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft1, array1d_x, array1d_x)

          do i=1,self%n_dofs(1)
             ind = ind+1
             y(ind) = real( array1d_x(i), kind=f64 )
          end do
       end do
    end do

  end subroutine ifft3d

  !> Helper function to invert 3x3 matrix
  subroutine invert3d( mat, mat_inv )
    sll_comp64, intent( in ) :: mat(3,3)
    sll_comp64, intent( out ) :: mat_inv(3,3)

    sll_comp64 :: det

    det = mat(1,1)*mat(2,2)*mat(3,3) + mat(1,2)*mat(2,3)*mat(3,1) + mat(1,3)*mat(2,1)*mat(3,2) - mat(1,3)*mat(2,2)*mat(3,1) - mat(3,3)*mat(1,2)*mat(2,1) - mat(1,1)*mat(2,3)*mat(3,2)

    mat_inv(1,1) = mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)
    mat_inv(1,2) = mat(1,3)*mat(3,2) - mat(1,2)*mat(3,3)
    mat_inv(1,3) = mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)

    mat_inv(2,1) = mat(2,3)*mat(3,1) - mat(2,1)*mat(3,3)
    mat_inv(2,2) = mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)
    mat_inv(2,3) = mat(1,3)*mat(2,1) - mat(1,1)*mat(2,3)

    mat_inv(3,1) = mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)
    mat_inv(3,2) = mat(1,2)*mat(3,1) - mat(1,1)*mat(3,2)
    mat_inv(3,3) = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)

    mat_inv = mat_inv/det

  end subroutine invert3d

  subroutine read_from_file_preconditioner(self, filename)
    class( sll_t_preconditioner_curl_solver_fft), intent( inout ) :: self !< Fourier solver
    character(len=*), intent( in ) :: filename !< filename

  end subroutine read_from_file_preconditioner

  subroutine print_info_preconditioner(self)
    class( sll_t_preconditioner_curl_solver_fft), intent( in ) :: self !< Fourier solver

  end subroutine print_info_preconditioner

  subroutine set_verbose_preconditioner( self, verbose )
    class( sll_t_preconditioner_curl_solver_fft), intent( inout ) :: self !< Fourier solver
    logical, intent( in ) :: verbose !< logical for convergence information

    self%verbose = verbose

  end subroutine set_verbose_preconditioner


  subroutine free_preconditioner(self)
    class( sll_t_preconditioner_curl_solver_fft), intent( inout ) :: self !< Fourier solver

  end subroutine free_preconditioner

end module sll_m_preconditioner_curl_solver_fft
