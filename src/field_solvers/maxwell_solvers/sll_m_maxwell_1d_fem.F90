!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations in 1D
!> @details
!> Contains the abstract class to create a Maxwell solver in 1D.

module sll_m_maxwell_1d_fem
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_arbitrary_degree_splines, only: &
    sll_s_uniform_b_splines_at_x

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_fft, only: &
    sll_t_fft, &
    sll_s_fft_init_r2r_1d, &
    sll_s_fft_exec_r2r_1d, &
    sll_s_fft_free, &
    sll_p_fft_forward, &
    sll_p_fft_backward

  use sll_m_gauss_legendre_integration, only: &
    sll_f_gauss_legendre_points_and_weights

  use sll_m_maxwell_1d_base, only: &
    sll_i_function_1d_real64, &
    sll_c_maxwell_1d_base

  implicit none

  public :: &
    sll_t_maxwell_1d_fem

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_c_maxwell_1d_base) :: sll_t_maxwell_1d_fem

     sll_real64 :: Lx          !< length of Periodic domain
     sll_real64 :: delta_x     !< cell size
     sll_int32  :: n_dofs      !< number of cells (and grid points)
     sll_int32  :: s_deg_0     !< spline degree 0-forms
     sll_int32  :: s_deg_1     !< spline degree 1-forms
     sll_real64, allocatable :: mass_0(:)      !< coefficients of 0-form mass matrix
     sll_real64, allocatable :: mass_1(:)      !< coefficients of 1-form mass matrix
     sll_real64, allocatable :: eig_mass0(:)   !< eigenvalues of circulant 0-form mass matrix
     sll_real64, allocatable :: eig_mass1(:)   !< eigenvalues of circulant 1-form mass matrix
     sll_real64, allocatable :: eig_weak_ampere(:)  !< eigenvalues of circulant update matrix for Ampere
     sll_real64, allocatable :: eig_weak_poisson(:) !< eigenvalues of circulant update matrix for Poisson
     type(sll_t_fft) :: plan_fw !< fft plan (forward)
     type(sll_t_fft) :: plan_bw !< fft plan (backward)
     sll_real64, allocatable :: wsave(:) !< scratch data
     sll_real64, allocatable :: work(:)  !< scratch data

   contains
     procedure :: &
          compute_E_from_B => sll_s_compute_e_from_b_1d_fem!< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => sll_s_compute_b_from_e_1d_fem!< Solve Faraday equation with E constant in time
     procedure :: &
          compute_E_from_rho => sll_s_compute_e_from_rho_1d_fem!< Solve E from rho using Poisson
     procedure :: &
          compute_E_from_j => compute_E_from_j_1d_fem !< Solve E from j 
     procedure :: &
          compute_rhs_from_function => sll_s_compute_fem_rhs
     procedure :: &
          L2norm_squared => L2norm_squared_1d_fem
     procedure :: &
          inner_product => inner_product_1d_fem
     procedure :: &
          L2projection => L2projection_1d_fem
     procedure :: &
          free => free_1d_fem
     procedure :: &
          init => init_1d_fem

  end type sll_t_maxwell_1d_fem

contains
  !> compute Ey from Bz using weak Ampere formulation 
  subroutine sll_s_compute_e_from_b_1d_fem(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_1d_fem) :: self
    sll_real64, intent(in)     :: delta_t   !< Time step
    sll_real64, intent(in)     :: field_in(:)  !< Bz
    sll_real64, intent(inout)  :: field_out(:)  !< Ey
    ! local variables
    sll_real64 :: coef
    
    ! Compute potential weak curl of Bz using eigenvalue of circulant inverse matrix
    call solve_circulant(self, self%eig_weak_ampere, field_in, self%work)
    ! Update bz from self value
    coef = delta_t/self%delta_x
    field_out =  field_out + coef*self%work

  end subroutine sll_s_compute_e_from_b_1d_fem

  !> Compute Bz from Ey using strong 1D Faraday equation for spline coefficients
  !> $B_z^{new}(x_j) = B_z^{old}(x_j) - \frac{\Delta t}{\Delta x} (E_y(x_j) - E_y(x_{j-1})  $
   subroutine sll_s_compute_b_from_e_1d_fem(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_1d_fem)  :: self
    sll_real64, intent(in)     :: delta_t
    sll_real64, intent(in)     :: field_in(:)  ! ey
    sll_real64, intent(inout)  :: field_out(:) ! bz 
    ! local variables
    sll_real64 :: coef
    sll_int32 :: i

    coef = delta_t/self%delta_x
    ! relation betwen spline coefficients for strong Ampere
    do i=2,self%n_dofs
       field_out(i) = field_out(i) + coef * ( field_in(i-1) - field_in(i) )
    end do
    ! treat Periodic point
    field_out(1) = field_out(1) + coef * ( field_in(self%n_dofs) - field_in(1) )
   end subroutine sll_s_compute_b_from_e_1d_fem

   !> Compute E_i from j_i integrated over the time interval using weak Ampere formulation
   subroutine compute_E_from_j_1d_fem(self, current, component, E)
     class(sll_t_maxwell_1d_fem)             :: self !< Maxwell solver class
     sll_real64,dimension(:),intent(in)    :: current !< Component \a component of the current integrated over time interval
     sll_int32, intent(in)                 :: component !< Component of the Efield to be computed
     sll_real64,dimension(:),intent(inout) :: E !< Updated electric field
     ! local variables
     sll_int32 :: i 
     sll_real64, dimension(self%n_dofs) :: eigvals

     ! Multiply by inverse mass matrix  using the eigenvalues of the circulant inverse matrix
     eigvals=0.0_f64
     if (component == 1) then
        do i=1,self%n_dofs/2+1
           eigvals(i) = 1.0_f64 / self%eig_mass1(i)
        end do
        call solve_circulant(self, eigvals, current, self%work)
     elseif (component == 2) then
        do i=1,self%n_dofs/2+1
           eigvals(i) = 1.0_f64 / self%eig_mass0(i)
        end do
        call solve_circulant(self, eigvals, current, self%work)
     else
        print*, 'Component ', component, 'not implemented in compute_E_from_j_1d_fem.'
     end if
     

     ! Update the electric field and scale
     E = E - self%work/self%delta_x

   end subroutine compute_E_from_j_1d_fem
  
   subroutine sll_s_compute_e_from_rho_1d_fem(self, E, rho )       
     class(sll_t_maxwell_1d_fem) :: self
     sll_real64,dimension(:),intent(in) :: rho
     sll_real64,dimension(:),intent(out) :: E
     ! local variables
     sll_int32 :: i 

     ! Compute potential phi from rho, using eigenvalue of circulant inverse matrix
     call solve_circulant(self, self%eig_weak_poisson, rho, self%work)
     ! Compute spline coefficients of Ex from those of phi
     do i=2,self%n_dofs
        E(i) =  (self%work(i-1) -  self%work(i)) !* (self%delta_x)
     end do
     ! treat Periodic point
     E(1) = (self%work(self%n_dofs) - self%work(1)) !* (self%delta_x)

   end subroutine sll_s_compute_e_from_rho_1d_fem

   subroutine solve_circulant(self, eigvals, rhs, res)
     class(sll_t_maxwell_1d_fem) :: self
     sll_real64, intent(in) :: eigvals(:)    ! eigenvalues of circulant matrix
     sll_real64, intent(in) :: rhs(:)
     sll_real64, intent(out) :: res(:)
     ! local variables
     sll_int32 :: k
     sll_real64 :: re, im 

     ! Compute res from rhs, using eigenvalue of circulant  matrix
     res = rhs
     ! Forward FFT
     call sll_s_fft_exec_r2r_1d ( self%plan_fw, res, self%wsave )
     self%wsave(1) = self%wsave(1) * eigvals(1)
     do k=2, self%n_dofs/2
        re = self%wsave(k) * eigvals(k) - &
             self%wsave(self%n_dofs-k+2) * eigvals(self%n_dofs-k+2)
        im = self%wsave(k) * eigvals(self%n_dofs-k+2) + &
             self%wsave(self%n_dofs-k+2) * eigvals(k)
        self%wsave(k) = re
        self%wsave(self%n_dofs-k+2) = im
     end do
     self%wsave(self%n_dofs/2+1) = self%wsave(self%n_dofs/2+1)*eigvals(self%n_dofs/2+1)
     ! Backward FFT 
     call sll_s_fft_exec_r2r_1d( self%plan_bw, self%wsave, res )
     ! normalize
     res = res / self%n_dofs
   end subroutine solve_circulant


   !> Compute the FEM right-hand-side for a given function f and periodic splines of given degree
   !> Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$ 
   subroutine sll_s_compute_fem_rhs(self, func, degree, coefs_dofs)
     class(sll_t_maxwell_1d_fem)             :: self
     procedure(sll_i_function_1d_real64) :: func
     sll_int32, intent(in) :: degree
     sll_real64, intent(out) :: coefs_dofs(:)  ! Finite Element right-hand-side
     ! local variables
     sll_int32 :: i,j,k
     sll_real64 :: coef
     sll_real64, dimension(2,degree+1) :: xw_gauss
     sll_real64, dimension(degree+1,degree+1) :: bspl

     ! take enough Gauss points so that projection is exact for splines of degree deg
     ! rescale on [0,1] for compatibility with B-splines
     xw_gauss = sll_f_gauss_legendre_points_and_weights(degree+1, 0.0_f64, 1.0_f64)
     ! Compute bsplines at gauss_points
     do k=1,degree+1
        call sll_s_uniform_b_splines_at_x(degree,xw_gauss(1,k), bspl(k,:))
        !print*, 'bs', bspl(k,:)
     end do

     ! Compute coefs_dofs = int f(x)N_i(x) 
     do i = 1, self%n_dofs
        coef=0.0_f64
        ! loop over support of B spline
        do j = 1, degree+1
           ! loop over Gauss points
           do k=1, degree+1
              coef = coef + xw_gauss(2,k)*func(self%delta_x*(xw_gauss(1,k) + i + j - 2)) * bspl(k,degree+2-j)
              !print*, i,j,k, xw_gauss(2,k), xw_gauss(1,k),f(self%delta_x*(xw_gauss(1,k) + i + j - 2)) 
           enddo
        enddo
        ! rescale by cell size
        coefs_dofs(i) = coef*self%delta_x
     enddo

   end subroutine sll_s_compute_fem_rhs

   !> Compute the L2 projection of a given function f on periodic splines of given degree
   subroutine L2projection_1d_fem(self, func, degree, coefs_dofs)
     class(sll_t_maxwell_1d_fem) :: self
     procedure(sll_i_function_1d_real64) :: func
     sll_int32, intent(in) :: degree
     sll_real64, intent(out) :: coefs_dofs(:)  ! spline coefficients of projection
     ! local variables
     sll_int32 :: i
     !sll_real64 :: coef
     !sll_real64, dimension(2,degree+1) :: xw_gauss
     !sll_real64, dimension(degree+1,degree+1) :: bspl
     sll_real64, dimension(self%n_dofs) :: eigvals

     ! Compute right-hand-side
     call sll_s_compute_fem_rhs(self, func, degree, self%work)

     ! Multiply by inverse mass matrix (! complex numbers stored in real array with fftpack ordering)
     eigvals=0.0_f64
     if (degree == self%s_deg_0) then
        do i=1,self%n_dofs/2+1
           eigvals(i) = 1.0_f64 / self%eig_mass0(i)
        end do
     elseif  (degree == self%s_deg_0-1) then
        do i=1,self%n_dofs/2+1
           eigvals(i) = 1.0_f64 / self%eig_mass1(i)
        end do
     else
        print*, 'degree ', degree, 'not availlable in maxwell_1d_fem object' 
     endif

     call solve_circulant(self, eigvals, self%work, coefs_dofs)
     ! Account for scaling in the mass matrix by dx
     coefs_dofs = coefs_dofs/self%delta_x

   end subroutine L2projection_1d_fem

   !> Compute square of the L2norm 
   function L2norm_squared_1d_fem(self, coefs_dofs, degree) result (r)
     class(sll_t_maxwell_1d_fem) :: self !< Maxwell solver object
     sll_real64 :: coefs_dofs(:) !< Coefficient for each DoF
     sll_int32  :: degree !< Specify the degree of the basis functions
     sll_real64 :: r !< Result: squared L2 norm

     ! Multiply coefficients by mass matrix (use diagonalization FFT and mass matrix eigenvalues)
     if (degree == self%s_deg_0 ) then

        call solve_circulant(self, self%eig_mass0, coefs_dofs, self%work)

     elseif (degree == self%s_deg_1) then

        call solve_circulant(self, self%eig_mass1, coefs_dofs, self%work)

     end if
     ! Multiply by the coefficients from the left (inner product)
     r = sum(coefs_dofs*self%work)
     ! Scale by delt_x
     r = r*self%delta_x

   end function L2norm_squared_1d_fem

   function inner_product_1d_fem( self, coefs1_dofs, coefs2_dofs, degree ) result (r)
     class(sll_t_maxwell_1d_fem) :: self !< Maxwell solver object
     sll_real64 :: coefs1_dofs(:) !< Coefficient for each DoF
     sll_real64 :: coefs2_dofs(:) !< Coefficient for each DoF
     sll_int32  :: degree !< Specify the degree of the basis functions
     sll_real64 :: r !< Result: squared L2 norm

     ! Multiply coefficients by mass matrix (use diagonalization FFT and mass matrix eigenvalues)
     if (degree == self%s_deg_0 ) then

        call solve_circulant(self, self%eig_mass0, coefs2_dofs, self%work)

     elseif (degree == self%s_deg_1) then

        call solve_circulant(self, self%eig_mass1, coefs2_dofs, self%work)

     end if
     ! Multiply by the coefficients from the left (inner product)
     r = sum(coefs1_dofs*self%work)
     ! Scale by delt_x
     r = r*self%delta_x
     
   end function inner_product_1d_fem
   

   subroutine init_1d_fem( self, domain, n_dofs, s_deg_0 )
     class(sll_t_maxwell_1d_fem), intent(out) :: self !< solver object
     sll_real64, intent(in) :: domain(2)     ! xmin, xmax
     sll_int32, intent(in) :: n_dofs  ! number of degrees of freedom (here number of cells and grid points)
     !sll_real64 :: delta_x ! cell size
     sll_int32, intent(in) :: s_deg_0 ! highest spline degree

    ! local variables
     sll_int32 :: ierr
     sll_int32 :: j, k ! loop variables
     sll_real64 :: coef0, coef1, sin_mode, cos_mode 

     self%n_dofs = n_dofs
     self%Lx = domain(2) - domain(1)
     self%delta_x = self%Lx / n_dofs
     self%s_deg_0 = s_deg_0
     self%s_deg_1 = s_deg_0 - 1

     SLL_ALLOCATE(self%mass_0(s_deg_0+1), ierr)
     SLL_ALLOCATE(self%mass_1(s_deg_0), ierr)

     select case(s_deg_0)
     case(1) ! linear and constant splines
        ! Upper diagonal coeficients  of linear spline mass matrix (from Eulerian numbers)
        self%mass_0(1) = 4.0_f64/6.0_f64 
        self%mass_0(2) = 1.0_f64/6.0_f64
        ! Upper diagonal coeficients  of constant spline mass matrix
        self%mass_1(1) = 1.0_f64 
     case(2) ! quadratic and linear splines
        ! Upper diagonal coeficients  of quadratic spline mass matrix (from Eulerian numbers)
        self%mass_0(1) = 66.0_f64/120.0_f64 
        self%mass_0(2) = 26.0_f64/120.0_f64
        self%mass_0(3) = 1.0_f64/120.0_f64
        ! Upper diagonal coeficients  of linear spline mass matrix (from Eulerian numbers)
        self%mass_1(1) = 4.0_f64/6.0_f64 
        self%mass_1(2) = 1.0_f64/6.0_f64
     case(3)
        ! Upper diagonal coeficients  of cubic spline mass matrix (from Eulerian numbers)
        self%mass_0(1) = 2416.0_f64/5040.0_f64 
        self%mass_0(2) = 1191.0_f64/5040.0_f64
        self%mass_0(3) = 120.0_f64/5040.0_f64
        self%mass_0(4) = 1.0_f64/5040.0_f64
        ! Upper diagonal coeficients  of quadratic spline mass matrix (from Eulerian numbers)
        self%mass_1(1) = 66.0_f64/120.0_f64 
        self%mass_1(2) = 26.0_f64/120.0_f64
        self%mass_1(3) = 1.0_f64/120.0_f64

     case default
        print*, 'sll_t_maxwell_1d_fem init: spline degree ', s_deg_0, ' not implemented'
     end select

     SLL_ALLOCATE(self%eig_mass0(n_dofs), ierr)
     SLL_ALLOCATE(self%eig_mass1(n_dofs), ierr)
     SLL_ALLOCATE(self%eig_weak_ampere(n_dofs), ierr)
     SLL_ALLOCATE(self%eig_weak_poisson(n_dofs), ierr)
     ! Initialise FFT     
     SLL_ALLOCATE(self%work(n_dofs),ierr)
     SLL_ALLOCATE(self%wsave(n_dofs),ierr)
     call sll_s_fft_init_r2r_1d( self%plan_fw, self%n_dofs, self%work, self%wsave, sll_p_fft_forward, normalized=.false. )
     call sll_s_fft_init_r2r_1d( self%plan_bw, self%n_dofs, self%work, self%work, sll_p_fft_backward, normalized=.false. )

     ! Compute eigenvalues of circulant Ampere update matrix M_0^{-1} D^T M_1
     ! and circulant Poisson Matrix (D^T M_1 D)^{-1}
     ! zero mode vanishes due to derivative matrix D^T
     self%eig_weak_ampere(1) = 0.0_f64 
     self%eig_weak_poisson(1) = 0.0_f64  ! Matrix is not invertible: 0-mode is set to 0
     self%eig_mass0(1) = 1.0_f64  ! sum of coefficents is one
     self%eig_mass1(1) = 1.0_f64  ! sum of coefficents is one

     do k=1, n_dofs/2 - 1
        coef0 =  self%mass_0(1)
        coef1 =  self%mass_1(1)
        do j=1,s_deg_0 - 1
           cos_mode = cos(2*sll_p_pi*j*k/n_dofs)
           coef0 = coef0 + 2* self%mass_0(j+1)*cos_mode
           coef1 = coef1 + 2* self%mass_1(j+1)*cos_mode
        enddo
        ! add last term for larger matrix
        j = s_deg_0
        coef0 = coef0 + 2* self%mass_0(j+1)*cos(2*sll_p_pi*j*k/n_dofs)
        ! compute eigenvalues
        self%eig_mass0(k+1) = coef0 ! real part
        self%eig_mass0(n_dofs-k+1) = 0.0_f64 ! imaginary part
        self%eig_mass1(k+1) = coef1 ! real part
        self%eig_mass1(n_dofs-k+1) = 0.0_f64 ! imaginary part
        cos_mode = cos(2*sll_p_pi*k/n_dofs)
        sin_mode = sin(2*sll_p_pi*k/n_dofs)
        self%eig_weak_ampere(k+1) =  (coef1 / coef0) * (1-cos_mode) ! real part
        self%eig_weak_ampere(n_dofs-k+1) =  -(coef1 / coef0) * sin_mode   ! imaginary part
        self%eig_weak_poisson(k+1) = 1.0_f64 / (coef1 * ((1-cos_mode)**2 + &
             sin_mode**2))  ! real part
        self%eig_weak_poisson(n_dofs-k+1) = 0.0_f64  ! imaginary part
     enddo
     ! N/2 mode
     coef0 =  self%mass_0(1)
     coef1 =  self%mass_1(1)
     do j=1, s_deg_0 - 1
        coef0 = coef0 + 2 * self%mass_0(j+1)*cos(sll_p_pi*j)
        coef1 = coef1 + 2 * self%mass_1(j+1)*cos(sll_p_pi*j)
     enddo
     
     ! add last term for larger matrix
     j = s_deg_0
     coef0 = coef0 + 2 * self%mass_0(j+1)*cos(sll_p_pi*j)

     ! compute eigenvalues
     self%eig_mass0(n_dofs/2+1) = coef0
     self%eig_mass1(n_dofs/2+1) = coef1
     self%eig_weak_ampere(n_dofs/2+1) = 2.0_f64 * (coef1 / coef0)
     self%eig_weak_poisson(n_dofs/2+1) = 1.0_f64 / (coef1 *4.0_f64) 

   end subroutine init_1d_fem

   subroutine free_1d_fem(self)
     class(sll_t_maxwell_1d_fem) :: self

     call sll_s_fft_free( self%plan_fw )
     call sll_s_fft_free( self%plan_bw )
     deallocate(self%mass_0)
     deallocate(self%mass_1)
     deallocate(self%eig_mass0)
     deallocate(self%eig_mass1)
     deallocate(self%eig_weak_ampere)
     deallocate(self%eig_weak_poisson)
     deallocate(self%wsave)
     deallocate(self%work)

   end subroutine free_1d_fem


end module sll_m_maxwell_1d_fem
