module sll_m_poisson_3d_fem_fft

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
    sll_t_poisson_3d_fem_fft

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  

  

  type :: sll_t_poisson_3d_fem_fft

     sll_int32 :: n_dofs(3)
     sll_int32 :: degree(3)
     sll_real64 :: delta_x(3)

     sll_real64, allocatable :: eig_values_dtm1d_1(:)    
     sll_real64, allocatable :: eig_values_dtm1d_2(:)
     sll_real64, allocatable :: eig_values_dtm1d_3(:)
     sll_comp64, allocatable :: eig_values_d1(:)    
     sll_comp64, allocatable :: eig_values_d2(:)
     sll_comp64, allocatable :: eig_values_d3(:)
     sll_real64, allocatable :: eig_values_mass_0_1(:)
     sll_real64, allocatable :: eig_values_mass_0_2(:)
     sll_real64, allocatable :: eig_values_mass_0_3(:)

     type(sll_t_fft) :: fft1
     type(sll_t_fft) :: fft2
     type(sll_t_fft) :: fft3
     type(sll_t_fft) :: ifft1
     type(sll_t_fft) :: ifft2
     type(sll_t_fft) :: ifft3

     ! Scratch data
     sll_comp64, allocatable :: array1d_x(:)
     sll_comp64, allocatable :: array1d_y(:)
     sll_comp64, allocatable :: array1d_z(:)
     sll_comp64, allocatable :: scratch(:,:,:)
     sll_comp64, allocatable :: scratchx(:,:,:)
     sll_comp64, allocatable :: scratchy(:,:,:)
     
   contains
     procedure :: compute_phi_from_rho => compute_phi_from_rho_fft
     procedure :: compute_e_from_rho   => compute_e_from_rho_fft
     procedure :: free => free_fft
     procedure :: init => init_fft
     procedure :: compute_rhs_from_function => compute_rhs_from_function_fft
     
  end type sll_t_poisson_3d_fem_fft

  abstract interface
     !> 3d real function
     function sll_i_function_3d_real64(x)
       use sll_m_working_precision ! can't pass a header file because the
                                 ! preprocessor prevents double inclusion.
                                 ! It is very rare.
       sll_real64             :: sll_i_function_3d_real64
       sll_real64, intent(in) :: x(3)
     end function sll_i_function_3d_real64
  end interface

contains

  subroutine compute_e_from_rho_fft( self,  rho, efield )
    class(sll_t_poisson_3d_fem_fft), intent( inout )  :: self
    sll_real64, intent(in) :: rho(:)
    sll_real64, intent(out) :: efield(:)
    
    sll_int32 :: n_dofs(3), ntotal
    sll_int32 :: i,j,k, ind
    sll_real64 :: factor1jk, factor2jk, eig_val
    
    n_dofs = self%n_dofs
    ntotal = product(n_dofs)
    
    ! Compute Fourier transform
    ind=0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             ind = ind+1
             self%array1d_x(i) = cmplx( rho(ind), 0_f64, kind=f64)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft1, self%array1d_x, self%array1d_x)
          do i=1,n_dofs(1)
             self%scratch(i,j,k) = self%array1d_x(i)
          end do
       end do
    end do
    do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             self%array1d_y(j) = self%scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft2, self%array1d_y, self%array1d_y)
          do j=1,n_dofs(2)
             self%scratch(i,j,k) = self%array1d_y(j)
          end do
       end do
    end do
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          do k=1,n_dofs(3)
             self%array1d_z(k) = self%scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft3, self%array1d_z, self%array1d_z)
          do k=1,n_dofs(3)
             self%scratch(i,j,k) = self%array1d_z(k)
          end do
       end do
    end do
      

    ! Apply inverse matrix of eigenvalues on mode
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          factor1jk = self%eig_values_dtm1d_2(j)*self%eig_values_mass_0_3(k) + &
               self%eig_values_dtm1d_3(k)*self%eig_values_mass_0_2(j)
          factor2jk = self%eig_values_mass_0_2(j) * self%eig_values_mass_0_3(k)
          do i=1,n_dofs(1)
             if ( i == 1 .and. j==1 .and. k==1 ) then
                self%scratch(i,j,k) = cmplx(0.0_f64, 0.0_f64, f64)
             else
                eig_val = factor2jk * self%eig_values_dtm1d_1(i) + &
                     factor1jk * self%eig_values_mass_0_1(i)
                self%scratch(i,j,k) = self%scratch(i,j,k) / cmplx(eig_val, 0.0_f64 , f64)
             end if

             self%scratchx(i,j,k) = -self%scratch(i,j,k)* self%eig_values_d1(i)
             self%scratchy(i,j,k) = -self%scratch(i,j,k)* self%eig_values_d2(j)
             self%scratch(i,j,k) = -self%scratch(i,j,k)* self%eig_values_d3(k)
             
          end do
       end do
    end do
    
      
    ! Compute inverse Fourier transfrom
    

    call ifft3d( self, self%scratchx, efield(1:ntotal) )
    call ifft3d( self, self%scratchy, efield(1+ntotal:2*ntotal) )
    call ifft3d( self, self%scratch, efield(1+2*ntotal:3*ntotal) )
   
      
  end subroutine compute_e_from_rho_fft

 subroutine compute_phi_from_rho_fft( self, rho, phi )
    class(sll_t_poisson_3d_fem_fft), intent( inout )  :: self
    sll_real64, intent(in) :: rho(:)
    sll_real64, intent(out) :: phi(:)
    
    sll_int32 :: n_dofs(3)
    sll_int32 :: i,j,k, ind
    sll_real64 :: factor1jk, factor2jk, eig_val


    n_dofs = self%n_dofs
    
    ! Compute Fourier transform
    ind=0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             ind = ind+1
             self%array1d_x(i) = cmplx( rho(ind), 0_f64, kind=f64)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft1, self%array1d_x, self%array1d_x)
          do i=1,n_dofs(1)
             self%scratch(i,j,k) = self%array1d_x(i)
          end do
       end do
    end do
    do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             self%array1d_y(j) = self%scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft2, self%array1d_y, self%array1d_y)
          do j=1,n_dofs(2)
             self%scratch(i,j,k) = self%array1d_y(j)
          end do
       end do
    end do
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          do k=1,n_dofs(3)
             self%array1d_z(k) = self%scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%fft3, self%array1d_z, self%array1d_z)
          do k=1,n_dofs(3)
             self%scratch(i,j,k) = self%array1d_z(k)
          end do
       end do
    end do
      

    ! Apply inverse matrix of eigenvalues on mode
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          factor1jk = self%eig_values_dtm1d_2(j)*self%eig_values_mass_0_3(k) + &
               self%eig_values_dtm1d_3(k)*self%eig_values_mass_0_2(j)
          factor2jk = self%eig_values_mass_0_2(j) * self%eig_values_mass_0_3(k)
          do i=1,n_dofs(1)
             if ( i == 1 .and. j==1 .and. k==1 ) then
                self%scratch(i,j,k) = cmplx(0.0_f64, 0.0_f64, f64)
             else
                eig_val = factor2jk * self%eig_values_dtm1d_1(i) + &
                     factor1jk * self%eig_values_mass_0_1(i)
                self%scratch(i,j,k) = self%scratch(i,j,k) / cmplx(eig_val, 0.0_f64, f64)
             end if
          end do
       end do
    end do

    call ifft3d( self, self%scratch, phi )
    
      
    end subroutine compute_phi_from_rho_fft

    

  subroutine free_fft( self )
    class(sll_t_poisson_3d_fem_fft), intent( inout ) :: self
      
    call sll_s_fft_free( self%fft1 )
    
      
  end subroutine free_fft

  subroutine init_fft( self, n_dofs, degree, delta_x)
    class(sll_t_poisson_3d_fem_fft), intent( out ) :: self
    sll_int32, intent( in ) :: n_dofs(3)
    sll_int32, intent( in ) :: degree(3)
    sll_real64, intent( in ) :: delta_x(3)

    
    sll_real64 :: mass_line0_1(degree(1)+1)
    sll_real64 :: mass_line1_1(degree(1))
    sll_real64 :: mass_line0_2(degree(2)+1)
    sll_real64 :: mass_line1_2(degree(2))
    sll_real64 :: mass_line0_3(degree(3)+1)
    sll_real64 :: mass_line1_3(degree(3))
    sll_real64 :: eig_values_mass_1_1(n_dofs(1))
    sll_real64 :: eig_values_mass_1_2(n_dofs(2))
    sll_real64 :: eig_values_mass_1_3(n_dofs(3))
    sll_real64 :: angle
    sll_int32 :: j

    self%n_dofs = n_dofs

    allocate(self%array1d_x(n_dofs(1)))
    allocate(self%array1d_y(n_dofs(2)))
    allocate(self%array1d_z(n_dofs(3)))
    allocate(self%scratch(n_dofs(1), n_dofs(2), n_dofs(3)))
    allocate(self%scratchx(n_dofs(1), n_dofs(2), n_dofs(3)))
    allocate(self%scratchy(n_dofs(1), n_dofs(2), n_dofs(3)))

    
    call sll_s_fft_init_c2c_1d( self%fft1, n_dofs(1), self%array1d_x, self%array1d_x, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft2, n_dofs(2), self%array1d_y, self%array1d_y, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft3, n_dofs(3), self%array1d_z, self%array1d_z, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%ifft1, n_dofs(1), self%array1d_x, self%array1d_x, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft2, n_dofs(2), self%array1d_y, self%array1d_y, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft3, n_dofs(3), self%array1d_z, self%array1d_z, &
         sll_p_fft_backward, normalized=.true.)

      ! Eigenvalues of mass matrices
    call sll_s_spline_fem_mass_line( degree(1), mass_line0_1 )
    call sll_s_spline_fem_mass_line( degree(2), mass_line0_2 )
    call sll_s_spline_fem_mass_line( degree(3), mass_line0_3 )
    call sll_s_spline_fem_mass_line( degree(1)-1, mass_line1_1 )
    call sll_s_spline_fem_mass_line( degree(2)-1, mass_line1_2 )
    call sll_s_spline_fem_mass_line( degree(3)-1, mass_line1_3 )

    
    allocate( self%eig_values_mass_0_1(n_dofs(1)) )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(1), degree(1), mass_line0_1*delta_x(1), self%eig_values_mass_0_1 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(1), degree(1)-1, mass_line1_1*delta_x(1), eig_values_mass_1_1 )  

    allocate( self%eig_values_mass_0_2(n_dofs(2)) )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(2), degree(2), mass_line0_2*delta_x(2), self%eig_values_mass_0_2 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(2), degree(2)-1, mass_line1_2*delta_x(2), eig_values_mass_1_2 )  

    allocate( self%eig_values_mass_0_3(n_dofs(3)) )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(3), degree(3), mass_line0_3*delta_x(3), self%eig_values_mass_0_3 )
    call sll_s_spline_fem_compute_mass_eig( n_dofs(3), degree(3)-1, mass_line1_3*delta_x(3), eig_values_mass_1_3 )


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

    allocate( self%eig_values_d3(n_dofs(3)) )
    allocate( self%eig_values_dtm1d_3(n_dofs(3)) )
    

    self%eig_values_d3(1) = cmplx(0.0_f64, 0.0_f64, f64)
    self%eig_values_dtm1d_3(1) = 0.0_f64
    do j=2,n_dofs(3)
       angle = sll_p_twopi*real(j-1,f64)/real(n_dofs(3), f64)
       
       self%eig_values_d3(j) = cmplx((1.0_f64 - cos(angle))/delta_x(3),sin(angle)/delta_x(3), f64 )
       self%eig_values_dtm1d_3(j) = 2.0_f64/delta_x(3)**2*(1.0_f64-cos(angle))* eig_values_mass_1_3(j)
       
    end do

    self%degree = degree
    self%delta_x = delta_x

  end subroutine init_fft


  
 !> Compute the FEM right-hand-side for a given function f and periodic splines of given degree
   !> Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$ 
   subroutine compute_rhs_from_function_fft(self, func, coefs_dofs)
     class(sll_t_poisson_3d_fem_fft)             :: self
     procedure(sll_i_function_3d_real64) :: func
     sll_real64, intent(out) :: coefs_dofs(:)  ! Finite Element right-hand-side

     ! local variables
     sll_int32 :: i1, i2, i3,j1, j2, j3, k1, k2, k3, k, counter
     sll_real64 :: coef
     sll_real64 :: xw_gauss_d1(2,self%degree(1)+1) 
     sll_real64 :: xw_gauss_d2(2,self%degree(2)+1)
     sll_real64 :: xw_gauss_d3(2,self%degree(3)+1)
     sll_real64 :: bspl_d1(self%degree(1)+1, self%degree(1)+1)
     sll_real64 :: bspl_d2(self%degree(2)+1, self%degree(2)+1)
     sll_real64 :: bspl_d3(self%degree(3)+1, self%degree(3)+1)
     sll_real64 :: volume

     volume = product(self%delta_x)

     ! take enough Gauss points so that projection is exact for splines of degree deg
     ! rescale on [0,1] for compatibility with B-splines
     xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights(self%degree(1)+1, 0.0_f64, 1.0_f64)
     ! Compute bsplines at gauss_points
     do k=1,self%degree(1)+1
        call sll_s_uniform_bsplines_eval_basis(self%degree(1),xw_gauss_d1(1,k), bspl_d1(k,:))
     end do

     xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights(self%degree(2)+1, 0.0_f64, 1.0_f64)
     ! Compute bsplines at gauss_points
     do k=1,self%degree(2)+1
        call sll_s_uniform_bsplines_eval_basis(self%degree(2),xw_gauss_d2(1,k), bspl_d2(k,:))
     end do

     xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights(self%degree(3)+1, 0.0_f64, 1.0_f64)
     ! Compute bsplines at gauss_points
     do k=1,self%degree(3)+1
        call sll_s_uniform_bsplines_eval_basis(self%degree(3),xw_gauss_d3(1,k), bspl_d3(k,:))
     end do



     counter = 1
     ! Compute coefs_dofs = int f(x)N_i(x) 
     do i3 = 1, self%n_dofs(3)
        do i2 = 1, self%n_dofs(2)
           do i1 = 1, self%n_dofs(1)
              coef=0.0_f64
              ! loop over support of B spline
              do j1 = 1, self%degree(1)+1
                 do j2 = 1, self%degree(2)+1
                    do j3 = 1, self%degree(3)+1
                       ! loop over Gauss points
                       do k1=1, self%degree(1)+1
                          do k2=1, self%degree(2)+1
                             do k3=1, self%degree(3)+1
                                coef = coef + xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* &
                                     xw_gauss_d3(2,k3) *&
                                     func([self%delta_x(1)*(xw_gauss_d1(1,k1) + real(i1 + j1 - 2,f64)), self%delta_x(2)*(xw_gauss_d2(1,k2) + real(i2 + j2 - 2,f64)), self%delta_x(3)*(xw_gauss_d3(1,k3) + real(i3 + j3 - 2,f64))] ) * &
                                     bspl_d1(k1,self%degree(1)+2-j1)*&
                                     bspl_d2(k2,self%degree(2)+2-j2)*&
                                     bspl_d3(k3,self%degree(3)+2-j3)
                    
                             enddo
                          enddo
                       end do
                    end do
                 end do
              end do
              ! rescale by cell size
              coefs_dofs(counter) = coef*volume
              counter = counter+1
           enddo
        end do
     end do

   end subroutine compute_rhs_from_function_fft


   !> Helper function for 3d fft

   subroutine ifft3d( self, inval, outval )
     class(sll_t_poisson_3d_fem_fft),intent( inout) :: self
     sll_comp64, intent( inout )  :: inval(:,:,:)
     sll_real64, intent( out ) :: outval(:)

     sll_int32 :: i,j,k,ind
     sll_int32 :: n_dofs(3)

     n_dofs = self%n_dofs
     
   ! Compute inverse Fourier transfrom
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          do k=1,n_dofs(3)
             self%array1d_z(k) = inval(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft3, self%array1d_z, self%array1d_z)
          do k=1,n_dofs(3)
             inval(i,j,k) = self%array1d_z(k)
          end do
       end do
    end do

    do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             self%array1d_y(j) = inval(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft2, self%array1d_y, self%array1d_y)
          do j=1,n_dofs(2)
             inval(i,j,k) = self%array1d_y(j)
          end do
       end do
    end do
    
    ind=0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             self%array1d_x(i) = inval(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( self%ifft1, self%array1d_x, self%array1d_x)
          
          do i=1,n_dofs(1)
             ind = ind+1
             outval(ind) = real( self%array1d_x(i), kind=f64 )
          end do
       end do
    end do

  end subroutine ifft3d
  
  
end module sll_m_poisson_3d_fem_fft
