module sll_m_poisson_2d_fem_fft

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_low_level_bsplines, only: &
    sll_s_uniform_bsplines_eval_basis
  
  use sll_m_constants, only : &
       sll_p_twopi
  
  use sll_m_fft

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights
  
  use sll_m_spline_fem_utilities
  
  implicit none

  public :: &
    sll_t_poisson_2d_fem_fft

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  

  

  type :: sll_t_poisson_2d_fem_fft

     sll_int32 :: n_dofs(2)
     sll_int32 :: degree
     sll_real64 :: delta_x(2)

     sll_real64, allocatable :: eig_values_dtm1d_1(:)    
     sll_real64, allocatable :: eig_values_dtm1d_2(:)
     sll_comp64, allocatable :: eig_values_d1(:)    
     sll_comp64, allocatable :: eig_values_d2(:)
     sll_real64, allocatable :: eig_values_mass_0_1(:)
     sll_real64, allocatable :: eig_values_mass_0_2(:)

     type(sll_t_fft) :: fft1
     type(sll_t_fft) :: fft2
     type(sll_t_fft) :: ifft1
     type(sll_t_fft) :: ifft2

     ! Scratch data
     sll_comp64, allocatable :: array1d_x(:)
     sll_comp64, allocatable :: array1d_y(:)
     sll_comp64, allocatable :: scratch(:,:)
     sll_comp64, allocatable :: scratchx(:,:)
     sll_comp64, allocatable :: scratchy(:,:)
     
   contains
     procedure :: compute_phi_from_rho => compute_phi_from_rho_fft
     procedure :: compute_e_from_rho   => compute_e_from_rho_fft
     procedure :: free => free_fft
     procedure :: init => init_fft
     procedure :: compute_rhs_from_function => compute_rhs_from_function_fft
     
  end type sll_t_poisson_2d_fem_fft

  abstract interface
     !> 3d real function
     function sll_i_function_2d_real64(x)
       use sll_m_working_precision ! can't pass a header file because the
                                 ! preprocessor prevents double inclusion.
                                 ! It is very rare.
       sll_real64             :: sll_i_function_3d_real64
       sll_real64, intent(in) :: x(2)
     end function sll_i_function_2d_real64
  end interface

contains

  subroutine compute_e_from_rho_fft( self,  rho, efield )
    class(sll_t_poisson_2d_fem_fft), intent( inout )  :: self
    sll_real64, intent(in) :: rho(:)
    sll_real64, intent(out) :: efield(:)
    
    sll_int32 :: n_dofs(2), ntotal
    sll_int32 :: i,j
    sll_real64 :: eig_val
    
    n_dofs = self%n_dofs
    ntotal = product(n_dofs)
    
    ! Compute Fourier transform
    call fft2d( self, rho )
      
    ! Apply inverse matrix of eigenvalues on mode
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          if ( i == 1 .and. j==1  ) then
             self%scratch(i,j) = cmplx(0.0_f64, 0.0_f64, f64)
          else
             eig_val = self%eig_values_dtm1d_1(i) * self%eig_values_mass_0_2(j) + &
                  self%eig_values_mass_0_1(i) * self%eig_values_dtm1d_2(j) 
             self%scratch(i,j) = self%scratch(i,j) / cmplx(eig_val,0._f64, f64)
          end if
          
          self%scratchx(i,j) = -self%scratch(i,j)* self%eig_values_d1(i)
          self%scratchy(i,j) = -self%scratch(i,j)* self%eig_values_d2(j)
          
       end do
    end do
    
      
    ! Compute inverse Fourier transfrom
    

    call ifft2d( self, self%scratchx, efield(1:ntotal) )
    call ifft2d( self, self%scratchy, efield(1+ntotal:2*ntotal) )
   
      
  end subroutine compute_e_from_rho_fft

  subroutine compute_phi_from_rho_fft( self, rho, phi )
    class(sll_t_poisson_2d_fem_fft), intent( inout )  :: self
    sll_real64, intent(in) :: rho(:)
    sll_real64, intent(out) :: phi(:)
    
    sll_int32 :: n_dofs(2)
    sll_int32 :: i,j
    sll_real64 :: eig_val


    n_dofs = self%n_dofs
    
    ! Compute Fourier transform
    call fft2d( self, rho )

    ! Apply inverse matrix of eigenvalues on mode
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          if ( i == 1 .and. j==1 ) then
             self%scratch(i,j) = cmplx(0.0_f64, 0.0_f64, f64)
          else
             eig_val = self%eig_values_dtm1d_1(i) * self%eig_values_mass_0_2(j) + &
                  self%eig_values_mass_0_1(i) * self%eig_values_dtm1d_2(j) 
             self%scratch(i,j) = self%scratch(i,j) / cmplx(eig_val, 0.0_f64, f64)
          end if
       end do
    end do

    call ifft2d( self, self%scratch, phi )
    
      
  end subroutine compute_phi_from_rho_fft

    

  subroutine free_fft( self )
    class(sll_t_poisson_2d_fem_fft), intent( inout ) :: self
      
    call sll_s_fft_free( self%fft1 )
    
      
  end subroutine free_fft

  subroutine init_fft( self, n_dofs, degree, delta_x)
    class(sll_t_poisson_2d_fem_fft), intent( out ) :: self
    sll_int32, intent( in ) :: n_dofs(2)
    sll_int32, intent( in ) :: degree
    sll_real64, intent( in ) :: delta_x(2)

    
    sll_real64 :: mass_line_0(degree+1)
    sll_real64 :: mass_line_1(degree)
    sll_real64 :: eig_values_mass_1_1(n_dofs(1))
    sll_real64 :: eig_values_mass_1_2(n_dofs(2))
    sll_real64 :: angle
    sll_int32 :: j

    self%n_dofs = n_dofs

    allocate(self%array1d_x(n_dofs(1)))
    allocate(self%array1d_y(n_dofs(2)))
    allocate(self%scratch(n_dofs(1), n_dofs(2)))
    allocate(self%scratchx(n_dofs(1), n_dofs(2)))
    allocate(self%scratchy(n_dofs(1), n_dofs(2)))

    
    call sll_s_fft_init_c2c_1d( self%fft1, n_dofs(1), self%array1d_x, self%array1d_x, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft2, n_dofs(2), self%array1d_y, self%array1d_y, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%ifft1, n_dofs(1), self%array1d_x, self%array1d_x, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft2, n_dofs(2), self%array1d_y, self%array1d_y, &
         sll_p_fft_backward, normalized=.true.)

    ! Eigenvalues of mass matrices
    call sll_s_spline_fem_mass_line( degree, mass_line_0 )
    call sll_s_spline_fem_mass_line( degree-1, mass_line_1 )

    
    allocate( self%eig_values_mass_0_1(n_dofs(1)) )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(1), degree, mass_line_0*delta_x(1), self%eig_values_mass_0_1 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(1), degree-1, mass_line_1*delta_x(1), eig_values_mass_1_1 )  

    allocate( self%eig_values_mass_0_2(n_dofs(2)) )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(2), degree, mass_line_0*delta_x(2), self%eig_values_mass_0_2 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(2), degree-1, mass_line_1*delta_x(2), eig_values_mass_1_2 )  


    allocate( self%eig_values_d1(n_dofs(1)) )
    allocate( self%eig_values_dtm1d_1(n_dofs(1)) )
    
    self%eig_values_d1(1) = cmplx(0.0_f64, 0.0_f64, f64)
    self%eig_values_dtm1d_1(1) = 0.0_f64
    do j=2,n_dofs(1)
       angle = sll_p_twopi*real(j-1,f64)/real(n_dofs(1), f64)
       
       self%eig_values_d1(j) = cmplx((1.0_f64 - cos(angle))/delta_x(1),sin(angle)/delta_x(1), f64 )
       self%eig_values_dtm1d_1(j) = 2.0_f64/delta_x(1)**2*(1.0_f64-cos(angle))* eig_values_mass_1_1(j)
       
    end do

    allocate( self%eig_values_d2(n_dofs(2)) )
    allocate( self%eig_values_dtm1d_2(n_dofs(2)) )

    self%eig_values_d2(1) = cmplx(0.0_f64, 0.0_f64, f64)
    self%eig_values_dtm1d_2(1) = 0.0_f64
    do j=2,n_dofs(2)
       angle = sll_p_twopi*real(j-1,f64)/real(n_dofs(2), f64)
       
       self%eig_values_d2(j) = cmplx((1.0_f64 - cos(angle))/delta_x(2),sin(angle)/delta_x(2), f64 )
       self%eig_values_dtm1d_2(j) = 2.0_f64/delta_x(2)**2*(1.0_f64-cos(angle))* eig_values_mass_1_2(j)
       
    end do

    self%degree = degree
    self%delta_x = delta_x

  end subroutine init_fft


  
 !> Compute the FEM right-hand-side for a given function f and periodic splines of given degree
   !> Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$ 
   subroutine compute_rhs_from_function_fft(self, func, coefs_dofs)
     class(sll_t_poisson_2d_fem_fft)             :: self
     procedure(sll_i_function_2d_real64) :: func
     sll_real64, intent(out) :: coefs_dofs(:)  ! Finite Element right-hand-side

     ! local variables
     sll_int32 :: i1, i2, j1, j2,  k1, k2,  k, counter
     sll_real64 :: coef
     sll_real64, allocatable :: xw_gauss_d1(:,:)
     sll_real64, allocatable :: bspl_d1(:,:)
     sll_real64 :: volume

     volume = product(self%delta_x)

     ! take enough Gauss points so that projection is exact for splines of degree deg
     ! rescale on [0,1] for compatibility with B-splines
     allocate(xw_gauss_d1(2,self%degree+1))
     allocate(bspl_d1(self%degree+1, self%degree+1))
     xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights(self%degree+1, 0.0_f64, 1.0_f64)
     ! Compute bsplines at gauss_points
     do k=1,self%degree+1
        call sll_s_uniform_bsplines_eval_basis(self%degree,xw_gauss_d1(1,k), bspl_d1(k,:))
     end do


     counter = 1
     ! Compute coefs_dofs = int f(x)N_i(x) 
     do i2 = 1, self%n_dofs(2)
        do i1 = 1, self%n_dofs(1)
           coef=0.0_f64
           ! loop over support of B spline
           do j1 = 1, self%degree+1
              do j2 = 1, self%degree+1
                 ! loop over Gauss points
                 do k1=1, self%degree+1
                    do k2=1, self%degree+1
                       coef = coef + xw_gauss_d1(2,k1)* xw_gauss_d1(2,k2)* &
                                     func([self%delta_x(1)*(xw_gauss_d1(1,k1) + real(i1 + j1 - 2,f64)), self%delta_x(2)*(xw_gauss_d1(1,k2) + real(i2 + j2 - 2,f64))] ) * &
                                     bspl_d1(k1,self%degree+2-j1)*&
                                     bspl_d1(k2,self%degree+2-j2)
                    
                    end do
                 end do
              end do
           end do
           ! rescale by cell size
           coefs_dofs(counter) = coef*volume
           counter = counter+1
        end do
     end do

   end subroutine compute_rhs_from_function_fft


   !> Helper function for 3d fft

   subroutine ifft2d( self, inval, outval )
     class(sll_t_poisson_2d_fem_fft),intent( inout) :: self
     sll_comp64, intent( inout )  :: inval(:,:)
     sll_real64, intent( out ) :: outval(:)

     sll_int32 :: i,j,ind
     sll_int32 :: n_dofs(2)

     n_dofs = self%n_dofs
     
   ! Compute inverse Fourier transfrom
    !do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             self%array1d_y(j) = inval(i,j)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft2, self%array1d_y, self%array1d_y)
          do j=1,n_dofs(2)
             inval(i,j) = self%array1d_y(j)
          end do
       end do
    !end do
    
    ind=0
    !do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             self%array1d_x(i) = inval(i,j)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft1, self%array1d_x, self%array1d_x)
          
          do i=1,n_dofs(1)
             ind = ind+1
             outval(ind) = real( self%array1d_x(i), kind=f64 )
          end do
       end do
    !end do

     end subroutine ifft2d
     

   subroutine fft2d( self, rho )
     class(sll_t_poisson_2d_fem_fft),intent( inout) :: self
     sll_real64, intent( in )  :: rho(:)

     sll_int32 :: i,j,ind
     sll_int32 :: n_dofs(2)

     n_dofs = self%n_dofs
     

     ind=0
     do j=1,n_dofs(2)
        do i=1,n_dofs(1)
           ind = ind+1
           self%array1d_x(i) = cmplx( rho(ind), 0_f64, kind=f64)
        end do
        call sll_s_fft_exec_c2c_1d( self%fft1, self%array1d_x, self%array1d_x)
        do i=1,n_dofs(1)
           self%scratch(i,j) = self%array1d_x(i)
        end do
     end do
     
        
     do i=1,n_dofs(1)
        do j=1,n_dofs(2)
           self%array1d_y(j) = self%scratch(i,j)
        end do
        call sll_s_fft_exec_c2c_1d( self%fft2, self%array1d_y, self%array1d_y)
        do j=1,n_dofs(2)
           self%scratch(i,j) = self%array1d_y(j)
        end do
     end do

     
   end subroutine fft2d
  
  
end module sll_m_poisson_2d_fem_fft
