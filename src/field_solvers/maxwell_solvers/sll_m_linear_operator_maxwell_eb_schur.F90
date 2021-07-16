!> This linear operator implements the compatible spline FEM operator for the curl part of Maxwell's equation (Schur complement operator) on uniform periodic grid 
!> The operator is implemented based on its diagonal form in Fouier space
!> It also contains a dot_inverse that applies the inverse of the matrix (by inversion in Fouier space
!> @author Katharina Kormann
module sll_m_linear_operator_maxwell_eb_schur
#include "sll_working_precision.h"
  use sll_m_linear_operator_abstract

  use sll_m_constants

  use sll_m_fft
  
  use sll_m_spline_fem_utilities
  
  use sll_m_maxwell_3d_base, only : &
       sll_c_maxwell_3d_base

  implicit none

  public :: sll_t_linear_operator_maxwell_eb_schur


  private

  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_maxwell_eb_schur
     class( sll_c_maxwell_3d_base ), pointer :: maxwell_solver => null()

     sll_real64 :: factor

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

     type(sll_t_fft) :: fft1
     type(sll_t_fft) :: fft2
     type(sll_t_fft) :: fft3
     type(sll_t_fft) :: ifft1
     type(sll_t_fft) :: ifft2
     type(sll_t_fft) :: ifft3
     
     
   contains
     procedure :: create => create_maxwell_eb
     procedure :: free => free_maxwell_eb
     procedure :: dot => dot_maxwell_eb !dot_maxwell_eb_fourier
     procedure :: print_info => print_info_maxwell_eb
     procedure :: dot_inverse => inverse_dot_maxwell_eb_fourier
     procedure :: dot_inv_mass_1_1
     procedure :: dot_inv_mass_1_2
     procedure :: dot_inv_mass_1_3

     
  end type sll_t_linear_operator_maxwell_eb_schur


contains

  subroutine create_maxwell_eb( self, maxwell_solver )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( inout ) :: self
    class( sll_c_maxwell_3d_base ), target, intent( in )  :: maxwell_solver

    sll_real64 :: angle, dx
    sll_int32 :: j, n_dofs, ndofs(3)
    sll_comp64 :: array1d_x(maxwell_solver%n_dofs(1))
    sll_comp64 :: array1d_y(maxwell_solver%n_dofs(2))
    sll_comp64 :: array1d_z(maxwell_solver%n_dofs(3))
    sll_real64, allocatable :: mass_line_0(:)
    sll_real64, allocatable :: mass_line_1(:)
    sll_int32 :: degree(3)
    self%maxwell_solver => maxwell_solver

    degree = self%maxwell_solver%s_deg_0
    self%n_rows = self%maxwell_solver%n_total*3
    self%n_cols = self%maxwell_solver%n_total*3

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols

    ! Fourier product

    ! Eigenvalues of derivative matrices
    dx = self%maxwell_solver%delta_x(1)
    n_dofs = self%maxwell_solver%n_dofs(1)
    
    allocate( self%eig_values_d1(1:n_dofs) )
    allocate( self%eig_values_d1t(1:n_dofs ) )
    
    self%eig_values_d1(1) = cmplx(0.0_f64, 0.0_f64, f64)
    self%eig_values_d1t(1) = cmplx(0.0_f64, 0.0_f64, f64)
    do j=2,n_dofs
       angle = sll_p_twopi*real(j-1,f64)/real(n_dofs, f64)
       
       self%eig_values_d1(j) = cmplx((1.0_f64 - cos(angle))/dx,sin(angle)/dx, f64 )
       self%eig_values_d1t(j) = cmplx((1.0_f64 - cos(angle))/dx,-sin(angle)/dx, f64 )
    end do
    
    dx = self%maxwell_solver%delta_x(2)
    n_dofs = self%maxwell_solver%n_dofs(2)
    
    allocate( self%eig_values_d2(1:n_dofs) )
    allocate( self%eig_values_d2t(1:n_dofs ) )
    
    self%eig_values_d2(1) =  cmplx(0.0_f64, 0.0_f64, f64)
    self%eig_values_d2t(1) =  cmplx(0.0_f64, 0.0_f64, f64)
    do j=2,n_dofs
       angle = sll_p_twopi*real(j-1,f64)/real(n_dofs, f64)
       
       self%eig_values_d2(j) = cmplx((1.0_f64 - cos(angle))/dx,sin(angle)/dx, f64 )
       self%eig_values_d2t(j) = cmplx((1.0_f64 - cos(angle))/dx,-sin(angle)/dx, f64 )
    end do

    dx = self%maxwell_solver%delta_x(3)
    n_dofs = self%maxwell_solver%n_dofs(3)
    
    allocate( self%eig_values_d3(1:n_dofs) )
    allocate( self%eig_values_d3t(1:n_dofs ) )

    self%eig_values_d3(1) = cmplx(0.0_f64, 0.0_f64, f64)
    self%eig_values_d3t(1) = cmplx(0.0_f64, 0.0_f64, f64)
    do j=2,n_dofs
       angle = sll_p_twopi*real(j-1,f64)/real(n_dofs, f64)
       
       self%eig_values_d3(j) = cmplx((1.0_f64 - cos(angle))/dx,sin(angle)/dx, f64 )
       self%eig_values_d3t(j) = cmplx((1.0_f64 - cos(angle))/dx,-sin(angle)/dx, f64 )
    end do
        
    n_dofs = self%maxwell_solver%n_dofs(1)
    dx = self%maxwell_solver%delta_x(1)
    allocate( self%eig_values_mass_0_1(n_dofs) )
    allocate( self%eig_values_mass_1_1(n_dofs) )
    allocate( mass_line_0(degree(1)+1) )
    allocate( mass_line_1(degree(1)) )
    ! Eigenvalues of mass matrices
    call sll_s_spline_fem_mass_line( degree(1), mass_line_0 )
    call sll_s_spline_fem_mass_line( degree(1)-1, mass_line_1 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs, degree(1),  mass_line_0*dx, self%eig_values_mass_0_1 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs, degree(1)-1, mass_line_1*dx, self%eig_values_mass_1_1 )
    deallocate( mass_line_0 )
    deallocate( mass_line_1 )
    
    n_dofs = self%maxwell_solver%n_dofs(2)
    dx = self%maxwell_solver%delta_x(2)
    allocate( self%eig_values_mass_0_2(n_dofs) )
    allocate( self%eig_values_mass_1_2(n_dofs) )
    allocate( mass_line_0(degree(2)+1) )
    allocate( mass_line_1(degree(2)) )
    ! Eigenvalues of mass matrices
    call sll_s_spline_fem_mass_line( degree(2), mass_line_0 )
    call sll_s_spline_fem_mass_line( degree(2)-1, mass_line_1 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs, degree(2), mass_line_0*dx, self%eig_values_mass_0_2 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs, degree(2)-1, mass_line_1*dx, self%eig_values_mass_1_2 )
    deallocate( mass_line_0 )
    deallocate( mass_line_1 )
    
    n_dofs = self%maxwell_solver%n_dofs(3)
    dx = self%maxwell_solver%delta_x(3)
    allocate( self%eig_values_mass_0_3(n_dofs) )
    allocate( self%eig_values_mass_1_3(n_dofs) )
    allocate( mass_line_0(degree(3)+1) )
    allocate( mass_line_1(degree(3)) )
    ! Eigenvalues of mass matrices
    call sll_s_spline_fem_mass_line( degree(3), mass_line_0 )
    call sll_s_spline_fem_mass_line( degree(3)-1, mass_line_1 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs, degree(3), mass_line_0*dx, self%eig_values_mass_0_3 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs, degree(3)-1, mass_line_1*dx, self%eig_values_mass_1_3 )  
    deallocate( mass_line_0 )
    deallocate( mass_line_1 )
    

    ! Initialize fft
    ndofs = self%maxwell_solver%n_dofs
    call sll_s_fft_init_c2c_1d( self%fft1, ndofs(1), array1d_x, array1d_x, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft2, ndofs(2), array1d_y, array1d_y, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft3, ndofs(3), array1d_z, array1d_z, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%ifft1, ndofs(1), array1d_x, array1d_x, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft2, ndofs(2), array1d_y, array1d_y, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft3, ndofs(3), array1d_z, array1d_z, &
         sll_p_fft_backward, normalized=.true.)
    
    
  end subroutine create_maxwell_eb

  subroutine free_maxwell_eb( self )
    class( sll_t_linear_operator_maxwell_eb_schur), intent( inout ) :: self

    self%maxwell_solver => null()
    
  end subroutine free_maxwell_eb

  
   subroutine dot_c( self, x, y )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( in ) :: self
    sll_real64, intent( in ) :: x(:)
    sll_real64, intent( inout ) :: y(:)

    call self%maxwell_solver%multiply_c( x, y )
    
  end subroutine dot_c

  subroutine dot_maxwell_eb( self, x, y )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( in ) :: self
    sll_real64, intent( in ) :: x(:)
    sll_real64, intent( out ) :: y(:)

    sll_int32 :: ndofs, comp, istart, iend, j
    sll_int32 :: deg(3)
    sll_real64 :: scratch(self%maxwell_solver%n_total*3)
    

    ndofs = self%maxwell_solver%n_total

    ! Compute C x
    call self%maxwell_solver%multiply_c( x, scratch )!y(ndofs*3+1:ndofs*6) )

    ! Compute M2 C x
    do comp=1,3
       istart = 1+(comp-1)*ndofs
       iend =  comp*ndofs
       deg = 2
       deg(comp) = 1
       call self%maxwell_solver%multiply_mass(  deg, scratch(istart:iend), y(istart:iend) )
      
    end do

    ! Compute C^T M2 C x
    call self%maxwell_solver%multiply_ct( y, scratch )

    ! Compute M1 x
    do comp=1,3
       istart = 1+(comp-1)*ndofs
       iend =  comp*ndofs
       deg = 1
       deg(comp) = 2
       call self%maxwell_solver%multiply_mass(  deg, x(istart:iend), y(istart:iend) )
      
    end do

    ! Sum up the two parts
    do j=1,ndofs*3
       y(j) = y(j) + self%factor * scratch(j)
    end do
    

    
  end subroutine dot_maxwell_eb

  !> Helper function
  subroutine fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, inde, x, scratch1 )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( in ) :: self
    sll_comp64, intent( inout ) :: array1d_x(:)
    sll_comp64, intent( inout ) :: array1d_y(:)
    sll_comp64, intent( inout ) :: array1d_z(:)
    sll_int32, intent( in ) :: n_dofs(3)
    sll_int32, intent( in ) :: inde
    sll_real64, intent( in ) :: x(:)
    sll_comp64, intent( out ) :: scratch1(:,:,:,:)

    sll_int32 :: ind, i,j,k
    
    
    ind=0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             ind = ind+1
             array1d_x(i) = cmplx( x(ind), 0_f64, kind=f64)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft1, array1d_x, array1d_x)
          do i=1,n_dofs(1)
             scratch1(i,j,k,inde) = array1d_x(i)
          end do
       end do
    end do
    do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             array1d_y(j) = scratch1(i,j,k,inde)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft2, array1d_y, array1d_y)
          do j=1,n_dofs(2)
             scratch1(i,j,k,inde) = array1d_y(j)
          end do
       end do
    end do
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          do k=1,n_dofs(3)
             array1d_z(k) = scratch1(i,j,k,inde)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft3, array1d_z, array1d_z)
          do k=1,n_dofs(3)
             scratch1(i,j,k,inde) = array1d_z(k)
          end do
       end do
    end do

  end subroutine fft3d

  !> Helper function
  subroutine ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, inde,  scratch, y )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( in ) :: self
    sll_comp64, intent( inout ) :: array1d_x(:)
    sll_comp64, intent( inout ) :: array1d_y(:)
    sll_comp64, intent( inout ) :: array1d_z(:)
    sll_int32, intent( in ) :: n_dofs(3)
    sll_int32, intent( in ) :: inde
    sll_real64, intent( out ) :: y(:)
    sll_comp64, intent( inout ) :: scratch(:,:,:,:)

    sll_int32 :: ind, i,j,k
    
    
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          do k=1,n_dofs(3)
             array1d_z(k) = scratch(i,j,k,inde)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft3, array1d_z, array1d_z)
          do k=1,n_dofs(3)
             scratch(i,j,k,inde) = array1d_z(k)
          end do
       end do
    end do

    do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             array1d_y(j) = scratch(i,j,k,inde)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft2, array1d_y, array1d_y)
          do j=1,n_dofs(2)
             scratch(i,j,k,inde) = array1d_y(j)
          end do
       end do
    end do
    
    ind=0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             array1d_x(i) = scratch(i,j,k,inde)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft1, array1d_x, array1d_x)
          
          do i=1,n_dofs(1)
             ind = ind+1
             y(ind) = real( array1d_x(i), kind=f64 )
          end do
       end do
    end do

  end subroutine ifft3d

    subroutine dot_maxwell_eb_fourier( self, x, y )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( in ) :: self
    sll_real64, intent( in ) :: x(:)
    sll_real64, intent( inout ) :: y(:)

    sll_int32 :: ndofs, n_dofs(3)
    sll_comp64 :: scratch1(self%maxwell_solver%n_dofs(1),self%maxwell_solver%n_dofs(2), self%maxwell_solver%n_dofs(3),3)
    sll_comp64 :: scratch2(self%maxwell_solver%n_dofs(1),self%maxwell_solver%n_dofs(2), self%maxwell_solver%n_dofs(3),3)
    sll_comp64 :: mat(3,3)
    sll_comp64 :: array1d_x(self%maxwell_solver%n_dofs(1))
    sll_comp64 :: array1d_y(self%maxwell_solver%n_dofs(2))
    sll_comp64 :: array1d_z(self%maxwell_solver%n_dofs(3))
    sll_int32 :: i,j,k
    

    ndofs = self%maxwell_solver%n_total
    n_dofs = self%maxwell_solver%n_dofs

    
    ! Compute Fourier transform of input
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, x(1:ndofs), scratch1 )
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 2, x(1+ndofs:2*ndofs), scratch1 )
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 3, x(1+2*ndofs:3*ndofs), scratch1 )

             
    ! Apply eigenvalues to Fourier coefficient
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)

             mat(1,1) = cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64) + &
                  cmplx(self%factor,kind=f64) * (self%eig_values_d3t(k)*self%eig_values_d3(k)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  self%eig_values_d2t(j)*self%eig_values_d2(j)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(1,2) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d2t(j)*self%eig_values_d1(i)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(1,3) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d3t(k)*self%eig_values_d1(i)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))
             
             mat(2,2) = cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64) + &
                  cmplx(self%factor,kind=f64) * (self%eig_values_d3t(k)*self%eig_values_d3(k)*&
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  self%eig_values_d1t(i)*self%eig_values_d1(i)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(2,1) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d1t(i)*self%eig_values_d2(j)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(2,3) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d3t(k)*self%eig_values_d2(j)* &
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))
             
             mat(3,3) = cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  cmplx(self%factor,kind=f64) * (self%eig_values_d1t(i)*self%eig_values_d1(i)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  self%eig_values_d2t(j)*self%eig_values_d2(j)*&
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))
             mat(3,2) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d2t(j)*self%eig_values_d3(k)* &
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))
             mat(3,1) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d1t(i)*self%eig_values_d3(k)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))

             scratch2(i,j,k,1) = mat(1,1)*scratch1(i,j,k,1) + mat(1,2)*scratch1(i,j,k,2) + &
                  mat(1,3)*scratch1(i,j,k,3)
             scratch2(i,j,k,2) = mat(2,1)*scratch1(i,j,k,1) + mat(2,2)*scratch1(i,j,k,2) + &
                  mat(2,3)*scratch1(i,j,k,3)
             scratch2(i,j,k,3) = mat(3,1)*scratch1(i,j,k,1) + mat(3,2)*scratch1(i,j,k,2) + &
                  mat(3,3)*scratch1(i,j,k,3)

          end do
       end do
    end do
    !scratch2 = scratch1
             
    ! Compute inverse Fourier transform of result
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, scratch2, y(1:ndofs) )
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 2, scratch2, y(1+ndofs:2*ndofs) )
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 3, scratch2, y(1+2*ndofs:3*ndofs) )

    
  end subroutine dot_maxwell_eb_fourier


  

  subroutine inverse_dot_maxwell_eb_fourier( self, x, y )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( in ) :: self
    sll_real64, intent( in ) :: x(:)
    sll_real64, intent( inout ) :: y(:)

    sll_int32 :: ndofs, n_dofs(3)
    sll_comp64 :: scratch1(self%maxwell_solver%n_dofs(1),self%maxwell_solver%n_dofs(2), self%maxwell_solver%n_dofs(3),3)
    sll_comp64 :: scratch2(self%maxwell_solver%n_dofs(1),self%maxwell_solver%n_dofs(2), self%maxwell_solver%n_dofs(3),3)
    sll_comp64 :: mat(3,3), imat(3,3)
    sll_comp64 :: array1d_x(self%maxwell_solver%n_dofs(1))
    sll_comp64 :: array1d_y(self%maxwell_solver%n_dofs(2))
    sll_comp64 :: array1d_z(self%maxwell_solver%n_dofs(3))
    sll_int32 :: i,j,k
    

    ndofs = self%maxwell_solver%n_total
    n_dofs = self%maxwell_solver%n_dofs

    
    ! Compute Fourier transform of input
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, x(1:ndofs), scratch1 )
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 2, x(1+ndofs:2*ndofs), scratch1 )
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 3, x(1+2*ndofs:3*ndofs), scratch1 )

             
    ! Apply eigenvalues to Fourier coefficient
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)

             mat(1,1) = cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64) + &
                  cmplx(self%factor,kind=f64) * (self%eig_values_d3t(k)*self%eig_values_d3(k)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  self%eig_values_d2t(j)*self%eig_values_d2(j)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(1,2) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d2t(j)*self%eig_values_d1(i)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(1,3) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d3t(k)*self%eig_values_d1(i)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))
             
             mat(2,2) = cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64) + &
                  cmplx(self%factor,kind=f64) * (self%eig_values_d3t(k)*self%eig_values_d3(k)*&
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  self%eig_values_d1t(i)*self%eig_values_d1(i)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(2,1) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d1t(i)*self%eig_values_d2(j)* &
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_0_3(k),kind=f64))
             mat(2,3) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d3t(k)*self%eig_values_d2(j)* &
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))
             
             mat(3,3) = cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  cmplx(self%factor,kind=f64) * (self%eig_values_d1t(i)*self%eig_values_d1(i)*&
                  cmplx(self%eig_values_mass_1_1(i)*self%eig_values_mass_0_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64) + &
                  self%eig_values_d2t(j)*self%eig_values_d2(j)*&
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))
             mat(3,2) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d2t(j)*self%eig_values_d3(k)* &
                  cmplx(self%eig_values_mass_0_1(i)*self%eig_values_mass_1_2(j)*&
                  self%eig_values_mass_1_3(k),kind=f64))
             mat(3,1) = - cmplx(self%factor,kind=f64) * ( &
                  self%eig_values_d1t(i)*self%eig_values_d3(k)* &
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
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, scratch2, y(1:ndofs) )
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 2, scratch2, y(1+ndofs:2*ndofs) )
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 3, scratch2, y(1+2*ndofs:3*ndofs) )

    
  end subroutine inverse_dot_maxwell_eb_fourier

  subroutine dot_c_fourier( self, x, y )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( in ) :: self
    sll_real64, intent( in ) :: x(:)
    sll_real64, intent( inout ) :: y(:)

    sll_int32 :: ndofs,  n_dofs(3)
    sll_comp64 :: scratch1(self%maxwell_solver%n_dofs(1),self%maxwell_solver%n_dofs(2), self%maxwell_solver%n_dofs(3),3)
    sll_comp64 :: scratch2(self%maxwell_solver%n_dofs(1),self%maxwell_solver%n_dofs(2), self%maxwell_solver%n_dofs(3),3)
    sll_comp64 :: mat(3,3)
    sll_comp64 :: array1d_x(self%maxwell_solver%n_dofs(1))
    sll_comp64 :: array1d_y(self%maxwell_solver%n_dofs(2))
    sll_comp64 :: array1d_z(self%maxwell_solver%n_dofs(3))
    sll_int32 :: i,j,k
    

    ndofs = self%maxwell_solver%n_total
    n_dofs = self%maxwell_solver%n_dofs

    
    ! Compute Fourier transform of input
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, x(1:ndofs), scratch1 )
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 2, x(1+ndofs:2*ndofs), scratch1 )
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 3, x(1+2*ndofs:3*ndofs), scratch1 )

    write(71,*) real(scratch1)
    write(72,*) aimag(scratch1)

    mat = cmplx(0.0_f64, 0.0_f64, f64)
             
    ! Apply eigenvalues to Fourier coefficient
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)

             mat(1,1) =  cmplx(0.0_f64, 0.0_f64, f64)
             mat(1,2) = -self%eig_values_d3(k)
             mat(1,3) = self%eig_values_d2(j)
             
             mat(2,2) = cmplx(0.0_f64, 0.0_f64, f64)
             mat(2,1) = self%eig_values_d3(k)
             mat(2,3) = -self%eig_values_d1(i)
             
             mat(3,3) = cmplx(0.0_f64, 0.0_f64, f64)
             mat(3,2) = self%eig_values_d1(i)
             mat(3,1) = - self%eig_values_d2(j)

             scratch2(i,j,k,1) = mat(1,1)*scratch1(i,j,k,1) + mat(1,2)*scratch1(i,j,k,2) + &
                  mat(1,3)*scratch1(i,j,k,3)
             scratch2(i,j,k,2) = mat(2,1)*scratch1(i,j,k,1) + mat(2,2)*scratch1(i,j,k,2) + &
                  mat(2,3)*scratch1(i,j,k,3)
             scratch2(i,j,k,3) = mat(3,1)*scratch1(i,j,k,1) + mat(3,2)*scratch1(i,j,k,2) + &
                  mat(3,3)*scratch1(i,j,k,3)

          end do
       end do
    end do
    !scratch2 = scratch1
             
    ! Compute inverse Fourier transform of result
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, scratch2, y(1:ndofs) )
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 2, scratch2, y(1+ndofs:2*ndofs) )
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 3, scratch2, y(1+2*ndofs:3*ndofs) )

    

  end subroutine dot_c_fourier
  
  
  subroutine print_info_maxwell_eb( self )
    class(sll_t_linear_operator_maxwell_eb_schur), intent(in) :: self 
  end subroutine print_info_maxwell_eb




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



  !> Product of inverse(M_1,1)
  subroutine dot_inv_mass_1_1( self, x, y )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( in ) :: self
    sll_real64, intent( in ) :: x(:)
    sll_real64, intent( inout ) :: y(:)

    sll_int32 :: ndofs, n_dofs(3)
    sll_comp64 :: scratch(self%maxwell_solver%n_dofs(1),self%maxwell_solver%n_dofs(2), self%maxwell_solver%n_dofs(3),1)
    sll_comp64 :: array1d_x(self%maxwell_solver%n_dofs(1))
    sll_comp64 :: array1d_y(self%maxwell_solver%n_dofs(2))
    sll_comp64 :: array1d_z(self%maxwell_solver%n_dofs(3))
    sll_int32 :: i,j,k
    

    ndofs = self%maxwell_solver%n_total
    n_dofs = self%maxwell_solver%n_dofs

    
    ! Compute Fourier transform of input
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, x(1:ndofs), scratch )
             
    ! Apply eigenvalues to Fourier coefficient
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)

             scratch(i,j,k,1) = scratch(i,j,k,1)/ &
                  cmplx(self%eig_values_mass_1_1(i)* &
                  self%eig_values_mass_0_2(j)* &
                  self%eig_values_mass_0_3(k),kind=f64 )

          end do
       end do
    end do
             
    ! Compute inverse Fourier transform of result
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, scratch, y(1:ndofs) )

  end subroutine dot_inv_mass_1_1

    !> Product of inverse(M_1,2)
  subroutine dot_inv_mass_1_2( self, x, y )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( in ) :: self
    sll_real64, intent( in ) :: x(:)
    sll_real64, intent( inout ) :: y(:)

    sll_int32 :: ndofs, n_dofs(3)
    sll_comp64 :: scratch(self%maxwell_solver%n_dofs(1),self%maxwell_solver%n_dofs(2), self%maxwell_solver%n_dofs(3),1)
    sll_comp64 :: array1d_x(self%maxwell_solver%n_dofs(1))
    sll_comp64 :: array1d_y(self%maxwell_solver%n_dofs(2))
    sll_comp64 :: array1d_z(self%maxwell_solver%n_dofs(3))
    sll_int32 :: i,j,k
    

    ndofs = self%maxwell_solver%n_total
    n_dofs = self%maxwell_solver%n_dofs

    
    ! Compute Fourier transform of input
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, x(1:ndofs), scratch )
             
    ! Apply eigenvalues to Fourier coefficient
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)

             scratch(i,j,k,1) = scratch(i,j,k,1)/ &
                  cmplx(self%eig_values_mass_0_1(i)* &
                  self%eig_values_mass_1_2(j)* &
                  self%eig_values_mass_0_3(k), kind=f64 )

          end do
       end do
    end do
             
    ! Compute inverse Fourier transform of result
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, scratch, y(1:ndofs) )

  end subroutine dot_inv_mass_1_2

    !> Product of inverse(M_1,3)
  subroutine dot_inv_mass_1_3( self, x, y )
    class(sll_t_linear_operator_maxwell_eb_schur), intent( in ) :: self
    sll_real64, intent( in ) :: x(:)
    sll_real64, intent( inout ) :: y(:)

    sll_int32 :: ndofs, n_dofs(3)
    sll_comp64 :: scratch(self%maxwell_solver%n_dofs(1),self%maxwell_solver%n_dofs(2), self%maxwell_solver%n_dofs(3),1)
    sll_comp64 :: array1d_x(self%maxwell_solver%n_dofs(1))
    sll_comp64 :: array1d_y(self%maxwell_solver%n_dofs(2))
    sll_comp64 :: array1d_z(self%maxwell_solver%n_dofs(3))
    sll_int32 :: i,j,k
    

    ndofs = self%maxwell_solver%n_total
    n_dofs = self%maxwell_solver%n_dofs

    
    ! Compute Fourier transform of input
    call fft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, x(1:ndofs), scratch )
             
    ! Apply eigenvalues to Fourier coefficient
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)

             scratch(i,j,k,1) = scratch(i,j,k,1)/ &
                  cmplx(self%eig_values_mass_0_1(i)* &
                  self%eig_values_mass_0_2(j)* &
                  self%eig_values_mass_1_3(k) , kind=f64 )

          end do
       end do
    end do
             
    ! Compute inverse Fourier transform of result
    call ifft3d( self, array1d_x, array1d_y, array1d_z, n_dofs, 1, scratch, y(1:ndofs) )

  end subroutine dot_inv_mass_1_3

  
end module sll_m_linear_operator_maxwell_eb_schur
