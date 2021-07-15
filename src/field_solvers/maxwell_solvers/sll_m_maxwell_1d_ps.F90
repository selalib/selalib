!> @ingroup maxwell_solvers
!> @brief
!> Solve Maxwell's equations in 1D based on a pseudospectral solver 
!> @author
!> Katharina Kormann

module sll_m_maxwell_1d_ps
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
       sll_p_pi, &
       sll_p_twopi

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
    sll_t_maxwell_1d_ps

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Maxwell solver class with pseudospectral method
  type, extends(sll_c_maxwell_1d_base) :: sll_t_maxwell_1d_ps

     ! Fourier plans
     type(sll_t_fft) :: plan_fw !< fft plan (forward) 
     type(sll_t_fft) :: plan_bw !< fft plan (backward)

     ! Eigenvalues
     sll_real64, allocatable :: eig_d(:) !< eigen values of the derivative matrix
     sll_real64, allocatable :: eig_interp1(:) !< eigen values of the interpolation matrix
     sll_real64, allocatable :: eig_interp1_inverse(:) !< eigen values of the inverse of the interpolation matrix
     sll_real64, allocatable :: eig_interp1t_inverse(:) !< eigen values of the inverse transpose of the interpolation matrix
     sll_real64, allocatable :: eig_poisson(:)!< eigen values of the Poisson matrix
     sll_real64, allocatable :: eig_schur_curl_part(:) ! Schur complement in curl part computation

     sll_real64, allocatable :: eig_mixed(:) !< eigen values for the product of 0 and 1 forms
     
     ! Scratch arrays
     sll_real64, allocatable :: work(:) !< scratch vector
     sll_real64, allocatable :: wsave(:) !< scratch vector

     !
     sll_real64 :: domain(2) !< Domain boundaries


   contains
     procedure :: &
          compute_e_b => sll_s_compute_e_b_1d!< Solve source-free Maxwell's equations
     procedure :: &
          compute_e_from_b => compute_field_from_field_1d_ps!< Solve source-free part of Ampere's law
     procedure :: &
          compute_b_from_e => compute_field_from_field_1d_ps!< Solve Faraday's law
     procedure :: &
          compute_curl_part => compute_curl_part_1d_ps!< Solve source-free Maxwell's equations
     procedure :: &
          compute_e_from_rho => compute_e_from_rho_1d_ps!< Compute E from rho through Poisson's equation
     procedure :: &
          compute_rho_from_e => compute_rho_from_e_1d_ps!< Compute rho from E by Poisson matrix multiply
     procedure :: &
          compute_e_from_j => compute_e_from_j_1d_ps!< Source-part of Ampere's equation
     procedure :: &
          compute_rhs_from_function!< Compute the integrals over a given function tested by the basis functions
     procedure :: &
          L2norm_squared => L2norm_squared_1d_ps !< Square of the L_2 norm defined by the mass matrix
     procedure :: &
          inner_product => inner_product_1d_ps !< Mass-matrix induced inner product
     procedure :: &
          L2projection
     procedure :: transform_dofs => transform_dofs_1d_ps !< Transformation of the dofs
     procedure :: multiply_mass => multiply_mass_1d_ps !< Multiply dof vector by mass matrix
     procedure :: invert_mass => invert_mass_1d_ps !< Multiply vector by inverse mass matrix
     procedure :: multiply_g !< Multiply dof vector by gradient matrix
     procedure :: multiply_gt !< Multiply dof vector by divergence matrix
     procedure :: free => free_1d_ps !< Free the Maxwell class
     procedure :: &
          init => init_1d_ps !< Initialize the Maxwell class

     procedure :: grad_proj0 !< 
     procedure :: proj1_grad !< 

     procedure :: compute_phi_from_rho => compute_phi_from_rho_1d_ps !< Compute phi from rho by Poisson's equation
     procedure :: compute_phi_from_j => compute_phi_from_j_1d_ps !< Compute phi from j by the dynamic version of the quasineutrality equation
     

  end type sll_t_maxwell_1d_ps

contains

  subroutine sll_s_compute_e_b_1d( self, delta_t, efield_dofs, bfield_dofs )
    class(sll_t_maxwell_1d_ps) :: self
    sll_real64, intent( in    )   :: delta_t      !< Time step
    sll_real64, intent( inout )   :: efield_dofs(:)  !< E
    sll_real64, intent( inout )   :: bfield_dofs(:) !< B

  end subroutine sll_s_compute_e_b_1d

  subroutine grad_proj0( self, in, out )
    class(sll_t_maxwell_1d_ps), intent(inout) :: self
    sll_real64, intent(in) :: in(:)
    sll_real64, intent(out) :: out(:)

    self%wsave = in / sqrt(real(self%n_cells,f64))
    call sll_s_fft_exec_r2r_1d( self%plan_fw, self%wsave, self%work )

    

    call complex_product_real( self%n_cells, self%eig_d, self%work, out )
  end subroutine grad_proj0
  
  subroutine proj1_grad( self, in, out )
    class(sll_t_maxwell_1d_ps), intent(inout) :: self
    sll_real64, intent(in) :: in(:)
    sll_real64, intent(out) :: out(:)

    self%wsave = in / sqrt(real(self%n_cells,f64))
    call sll_s_fft_exec_r2r_1d( self%plan_fw, self%wsave, self%work )

    call complex_product_real( self%n_cells, self%eig_d, self%work, self%wsave )
    call complex_product_real( self%n_cells, self%eig_interp1_inverse, self%wsave, out )
    
  end subroutine proj1_grad


  subroutine init_1d_ps( self, domain, n_dofs )
    class(sll_t_maxwell_1d_ps), intent(out) :: self !< solver object
    sll_real64, intent(in) :: domain(2)     ! xmin, xmax
    sll_int32, intent(in) :: n_dofs  ! number of degrees of freedom (here number of cells and grid points)
    ! local variables
    sll_int :: ierr, k
    sll_real64 :: cos_mode, sin_mode, factor, modulus, ni

    self%s_deg_0 = -1
    
    self%n_dofs0 = n_dofs
    self%n_dofs1 = n_dofs
    self%n_cells = n_dofs
    ni = 1.0_f64/real(n_dofs,f64)
    self%domain = domain
    self%Lx = domain(2)-domain(1)
    self%delta_x = self%Lx * ni
    !print*, 'dx',self%delta_x
    ! Initialise FFT     
    SLL_ALLOCATE(self%work(n_dofs),ierr)
    SLL_ALLOCATE(self%wsave(n_dofs),ierr)
    call sll_s_fft_init_r2r_1d( self%plan_fw, self%n_cells, self%work, self%wsave, sll_p_fft_forward, normalized=.false. )
    call sll_s_fft_init_r2r_1d( self%plan_bw, self%n_cells, self%work, self%work, sll_p_fft_backward, normalized=.false. )

    ! Eigenvalues
    SLL_ALLOCATE( self%eig_d(n_dofs) , ierr )
    SLL_ALLOCATE( self%eig_interp1(n_dofs) , ierr )
    SLL_ALLOCATE( self%eig_interp1_inverse(n_dofs) , ierr )
    SLL_ALLOCATE( self%eig_interp1t_inverse(n_dofs) , ierr )
    SLL_ALLOCATE( self%eig_poisson(n_dofs) , ierr )
    SLL_ALLOCATE( self%eig_schur_curl_part(n_dofs) , ierr )
    SLL_ALLOCATE( self%eig_mixed(n_dofs), ierr )


    ! Set the eigenvalues
    self%eig_d(1) = 0.0_f64
    self%eig_interp1(1) = 1.0_f64!self%Lx/n_dofs
    self%eig_interp1_inverse(1) = 1.0_f64!n_dofs/self%Lx
    self%eig_mixed(1) = 1.0_f64
    
    do k=1, (n_dofs+1)/2 - 1
       cos_mode = cos(sll_p_twopi*ni*real(k, f64))
       sin_mode = sin(sll_p_twopi*ni*real(k, f64))
       factor = sll_p_twopi/self%Lx*real(k,f64)
       modulus = (cos_mode-1.0_f64)**2 + sin_mode**2
       ! real part first derivative
       self%eig_d(k+1) = 0.0_f64
       ! imaginary part first derivative
       self%eig_d(n_dofs-k+1) = factor

       factor = factor  * self%Lx*ni

       ! real part for interpolation 1
       self%eig_interp1(k+1) =  sin_mode/factor
       self%eig_interp1(n_dofs-k+1) =  (cos_mode-1.0_f64)/factor

       self%eig_mixed(k+1) = sin_mode/factor
       
       self%eig_interp1_inverse(k+1) =  sin_mode/modulus * factor
       self%eig_interp1_inverse(n_dofs-k+1) = (1.0_f64-cos_mode)/modulus * factor

       self%eig_mixed(n_dofs-k+1) = 0.0_f64
    end do

    if ( modulo(n_dofs,2) == 0 ) then
       self%eig_d(n_dofs/2 + 1) = 0.0_f64
       self%eig_interp1(n_dofs/2+1) =  0.0_f64

       self%eig_interp1_inverse(n_dofs/2+1) = 0.0_f64

       self%eig_interp1t_inverse(1:n_dofs/2+1) = self%eig_interp1_inverse(1:n_dofs/2+1)
       self%eig_interp1t_inverse(n_dofs/2+2:n_dofs) = -self%eig_interp1_inverse(n_dofs/2+2:n_dofs)
       self%eig_poisson(n_dofs/2 + 1) = 0.0_f64

       self%eig_mixed(n_dofs/2 + 1) = 0.0_f64
    else
       self%eig_interp1t_inverse(1:n_dofs/2+1) = self%eig_interp1_inverse(1:n_dofs/2+1)
       self%eig_interp1t_inverse(n_dofs/2+2:n_dofs) = -self%eig_interp1_inverse(n_dofs/2+2:n_dofs)
    end if

    self%eig_poisson(1) = 0.0_f64

    do k=1,(n_dofs+1)/2 - 1
       self%eig_poisson(k+1) = self%eig_interp1_inverse(n_dofs-k+1) / self%eig_d(n_dofs-k+1)
       self%eig_poisson(n_dofs-k+1) = -self%eig_interp1_inverse(k+1) / self%eig_d(n_dofs-k+1)
    end do
    
  end subroutine init_1d_ps


  subroutine set_eig_schur_curl_part( self, delta_t )
    class( sll_t_maxwell_1d_ps ), intent( inout ) :: self
    sll_real64, intent( in ) :: delta_t

    sll_int32 :: k, n_dofs
    sll_real64 :: factor

    n_dofs = self%n_cells
    factor = 0.25_f64 * delta_t**2

    self%eig_schur_curl_part = 0.0_f64
    self%eig_schur_curl_part(1) = 1.0_f64
    do k=1, (n_dofs+1)/2 - 1
       self%eig_schur_curl_part(k+1) = 1.0_f64/(1.0_f64 + factor * self%eig_d(n_dofs-k+1)**2)
       !self%eig_schur_curl_part(n_dofs-k+1) = 0.0_f64
    end do

  end subroutine set_eig_schur_curl_part

  subroutine free_1d_ps( self )
     class(sll_t_maxwell_1d_ps) :: self !< Maxwell solver object
    
   end subroutine free_1d_ps
     
  function inner_product_1d_ps( self, coefs1_dofs, coefs2_dofs, degree, degree2 ) result (r)
     class(sll_t_maxwell_1d_ps) :: self !< Maxwell solver object
     sll_real64 :: coefs1_dofs(:) !< Coefficient for each DoF
     sll_real64 :: coefs2_dofs(:) !< Coefficient for each DoF
     sll_int32  :: degree !< Specify the degree of the basis functions
     sll_int32, optional  :: degree2 !< Specify the degree of the basis functions
     sll_real64 :: r !< Result: innter product

     if (present(degree2) ) then
        if  (degree .ne. degree2 ) then
           call complex_product_real( self%n_cells, self%eig_mixed, coefs1_dofs, self%wsave )
           r = self%wsave(1) * coefs2_dofs(1) * self%delta_x
           r = r +  sum(self%wsave(2:self%n_cells) * coefs2_dofs(2:self%n_cells) ) * self%delta_x * 2.0_f64
        else
           r = coefs1_dofs(1) * coefs2_dofs(1) * self%delta_x
           r = r +  sum(coefs1_dofs(2:self%n_cells) * coefs2_dofs(2:self%n_cells) ) * self%delta_x * 2.0_f64
        end if
     else
        r = coefs1_dofs(1) * coefs2_dofs(1) * self%delta_x
        r = r +  sum(coefs1_dofs(2:self%n_cells) * coefs2_dofs(2:self%n_cells) ) * self%delta_x * 2.0_f64
     end if

   end function inner_product_1d_ps

   function L2norm_squared_1d_ps(self, coefs_dofs, degree) result (r)
     class(sll_t_maxwell_1d_ps) :: self !< Maxwell solver object
     sll_real64 :: coefs_dofs(:) !< Coefficient for each DoF
     sll_int32  :: degree !< Specify the degree of the basis functions
     sll_real64 :: r !< Result: squared L2 norm

     r = coefs_dofs(1) * coefs_dofs(1) * self%delta_x
     r = r + sum(coefs_dofs(2:self%n_cells) * coefs_dofs(2:self%n_cells) ) * self%delta_x * 2.0_f64
     
   end function L2norm_squared_1d_ps

 !> Compute E_i from j_i integrated over the time interval using weak Ampere formulation  
   subroutine compute_E_from_j_1d_ps(self, current, component, E)
     class(sll_t_maxwell_1d_ps)           :: self !< Maxwell solver class
     sll_real64,dimension(:),intent(in)    :: current !< Component \a component of the current integrated over time interval
     sll_int32, intent(in)                 :: component !< Component of the Efield to be computed
     sll_real64,dimension(:),intent(inout) :: E !< Updated electric field


     if ( component == 1 ) then
        self%wsave = current
        call sll_s_fft_exec_r2r_1d ( self%plan_fw, self%wsave, self%work )
        E = E - self%work / sqrt(real(self%n_cells,f64))
     elseif ( component == 2) then
        self%wsave = current
        call sll_s_fft_exec_r2r_1d ( self%plan_fw, self%wsave, self%work )
        call complex_product_real( self%n_cells, self%eig_interp1_inverse, self%work, self%wsave )
        E = E  - self%wsave / sqrt(real(self%n_cells,f64))
     else
        print*, 'Component ', component, 'not implemented in compute_E_from_j_1d_ps.'
     end if
     


   end subroutine compute_E_from_j_1d_ps


   subroutine compute_field_from_field_1d_ps(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_1d_ps) :: self
    sll_real64, intent(in)     :: delta_t   !< Time step
    sll_real64, intent(in)     :: field_in(:)  !< Bz
    sll_real64, intent(inout)  :: field_out(:)  !< Ey

    call complex_product_real ( self%n_cells, self%eig_d, field_in, self%work )
    field_out = field_out - delta_t * self%work

  end subroutine compute_field_from_field_1d_ps

  subroutine compute_curl_part_1d_ps( self, delta_t, efield, bfield, betar )
    class(sll_t_maxwell_1d_ps) :: self
    sll_real64, intent(in)     :: delta_t   !< Time step
    sll_real64, intent(inout)  :: efield(:)  !< Ey
    sll_real64, intent(inout)  :: bfield(:)  !< Bz
    sll_real64, optional       :: betar      !< 1/beta
    !local variables
    sll_real64 :: factor
    sll_real64 :: dth
    
    if( present(betar) ) then
       factor = betar
    else
       factor = 1._f64
    end if
     
    dth = 0.5_f64 * delta_t
    
    call set_eig_schur_curl_part( self, delta_t )

    self%wsave = bfield
    call self%compute_b_from_e ( dth, efield, self%wsave )
    call self%compute_e_from_b ( dth, bfield*factor, efield )

    call self%compute_b_from_e ( dth, efield, self%wsave )

    call complex_product_real ( self%n_cells, self%eig_schur_curl_part, self%wsave, bfield )
    call self%compute_e_from_b( dth, bfield*factor, efield )

  end subroutine compute_curl_part_1d_ps


  subroutine transform_dofs_1d_ps(self, in, out, degree)
     class(sll_t_maxwell_1d_ps), intent(inout) :: self
     sll_int32, intent(in) :: degree ! this is 0 for 0 form and 1 for 1 form
     sll_real64, intent(in) :: in(:)  ! spline coefficients of projection
     sll_real64, intent(out) :: out(:)  ! spline coefficients of projection
     ! local variables

     if ( degree == 0 ) then
        ! Backward FFT
        self%wsave = in
        call sll_s_fft_exec_r2r_1d( self%plan_bw, self%wsave, out )
        ! normalize
        out = out / sqrt(real(self%n_cells, f64)) * self%delta_x
     elseif ( degree == 1 ) then
        call complex_product_real( self%n_cells, self%eig_interp1t_inverse, in, self%wsave )
        ! Backward FFT 
        call sll_s_fft_exec_r2r_1d( self%plan_bw, self%wsave, out )
        ! normalize
        out = out / sqrt(real(self%n_cells, f64)) * self%delta_x
     elseif ( degree == 2 ) then
        !self%work = in
        call complex_product_real( self%n_cells, self%eig_interp1t_inverse, in, self%work )
        call complex_product_real( self%n_cells, self%eig_mixed, self%work, self%wsave )
        call sll_s_fft_exec_r2r_1d( self%plan_bw, self%wsave, out )
        ! normalize
        out = out / sqrt(real(self%n_cells, f64)) * self%delta_x
     elseif ( degree == 3 ) then
        call complex_product_real( self%n_cells, self%eig_interp1t_inverse, in, self%work )
        call complex_product_real( self%n_cells, self%eig_mixed, self%work, self%wsave )
        call sll_s_fft_exec_r2r_1d( self%plan_bw, self%wsave, out )
        ! normalize
        out = out / sqrt(real(self%n_cells, f64)) * self%delta_x
     else
        print*, 'Component', degree, 'not defined in maxwell_1d_ps.'
     end if
     
   end subroutine transform_dofs_1d_ps


   subroutine compute_e_from_rho_1d_ps(self, field_in, field_out )       
     class(sll_t_maxwell_1d_ps) :: self
     sll_real64, intent(in)  :: field_in(:)  !rho
     sll_real64, intent(out) :: field_out(:) !E

     ! Fourier transform rho
     self%wsave = field_in/sqrt(real(self%n_cells,f64))
     call sll_s_fft_exec_r2r_1d ( self%plan_fw, self%wsave, self%work )
     ! Invert the Interpolation
     call complex_product_real( self%n_cells, self%eig_poisson, self%work, field_out )
     
   end subroutine compute_e_from_rho_1d_ps
   
   subroutine compute_rho_from_e_1d_ps(self, field_in, field_out )       
     class(sll_t_maxwell_1d_ps) :: self
     sll_real64, intent(in)  :: field_in(:)  !rho
     sll_real64, intent(out) :: field_out(:) !E

     call complex_product_real( self%n_cells, self%eig_d, field_in, self%wsave )
     call complex_product_real( self%n_cells, self%eig_interp1, self%wsave, self%work )
     call sll_s_fft_exec_r2r_1d( self%plan_bw, self%work, field_out )
     field_out = field_out / sqrt(real(self%n_cells,f64))

   end subroutine compute_rho_from_e_1d_ps
     
   ! Helper function
   subroutine complex_product_real(n_dofs, eigvals, in, out)
     sll_int32, intent(in)  :: n_dofs
     sll_real64, intent(in) :: eigvals(:)    ! eigenvalues of circulant matrix
     sll_real64, intent(in) :: in(:)
     sll_real64, intent(out) :: out(:)
     ! local variables
     sll_int32 :: k
     sll_real64 :: re, im 


     out(1) = in(1) * eigvals(1)
     do k=2, (n_dofs+1)/2
        re = in(k) * eigvals(k) - &
             in(n_dofs-k+2) * eigvals(n_dofs-k+2)
        im = in(k) * eigvals(n_dofs-k+2) + &
             in(n_dofs-k+2) * eigvals(k)
        out(k) = re
        out(n_dofs-k+2) = im
     end do
     if ( modulo(n_dofs,2) == 0 ) then
        out(n_dofs/2+1) = in(n_dofs/2+1)*eigvals(n_dofs/2+1)
     end if
     
   end subroutine complex_product_real


   subroutine multiply_mass_1d_ps( self,  in, out, degree)
     class(sll_t_maxwell_1d_ps), intent(inout) :: self !< Maxwell solver object
     sll_real64, intent(in)  :: in(:) !< Coefficient for each DoF
     sll_real64, intent(out) :: out(:) !< Coefficient for each DoF
     sll_int32,  intent(in)  :: degree !< Specify the degree of the basis functions

     out = in * self%delta_x!self%Lx

   end subroutine multiply_mass_1d_ps

   subroutine invert_mass_1d_ps(self, in, out, degree)
     class(sll_t_maxwell_1d_ps), intent(inout) :: self
     sll_int32, intent(in) :: degree
     sll_real64, intent(in) :: in(:)  ! spline coefficients of projection
     sll_real64, intent(out) :: out(:)  ! spline coefficients of projection

     out = in / self%delta_x !self%Lx

   end subroutine invert_mass_1d_ps

    subroutine multiply_g( self,  in, out)
      class(sll_t_maxwell_1d_ps), intent(in) :: self !< Maxwell_Clamped solver object
      sll_real64, intent(in)  :: in(:) !< Coefficient for each DoF
      sll_real64, intent(out) :: out(:) !< Coefficient for each DoF

      call complex_product_real( self%n_cells, self%eig_d, in, out )

    end subroutine multiply_g

    
    subroutine multiply_gt( self,  in, out)
      class(sll_t_maxwell_1d_ps), intent(in) :: self !< Maxwell_Clamped solver object
      sll_real64, intent(in)  :: in(:) !< Coefficient for each DoF
      sll_real64, intent(out) :: out(:) !< Coefficient for each DoF

      call complex_product_real( self%n_cells, self%eig_d, in, out )
      out = - out

    end subroutine multiply_gt

    !> Compute the FEM right-hand-side for a given function f and periodic splines of given degree
   !> Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$ 
   subroutine compute_rhs_from_function(self, func, degree, coefs_dofs)
     class(sll_t_maxwell_1d_ps)             :: self
     procedure(sll_i_function_1d_real64) :: func
     sll_int32, intent(in) :: degree
     sll_real64, intent(out) :: coefs_dofs(:)  ! Finite Element right-hand-side

     ! local variables
     sll_int32 :: k
     sll_real64 :: x

     do k=1, self%n_cells
        x = self%domain(1) + self%delta_x * real(k-1, f64)
        coefs_dofs(k) = func( x )
     end do
     
     
   end subroutine compute_rhs_from_function

    !> Compute the L2 projection of a given function f on periodic splines of given degree
   subroutine L2projection(self, func, degree, coefs_dofs)
     class(sll_t_maxwell_1d_ps) :: self
     procedure(sll_i_function_1d_real64) :: func
     sll_int32, intent(in) :: degree
     sll_real64, intent(out) :: coefs_dofs(:)  ! spline coefficients of projection

     
     call self%compute_rhs_from_function( func, degree, self%wsave )
     
     call sll_s_fft_exec_r2r_1d ( self%plan_fw, self%wsave, coefs_dofs )
     coefs_dofs = coefs_dofs/sqrt(real(self%n_cells,f64))

   end subroutine L2projection

   
  !> For model with adiabatic electrons
  subroutine compute_phi_from_rho_1d_ps( self, in, phi, efield )
    class(sll_t_maxwell_1d_ps) :: self
    sll_real64, intent(in)                     :: in(:)
    sll_real64, intent(out)                    :: phi(:)
    sll_real64, intent(out)                    :: efield(:)
    
    ! Compute phi by inverting the mass matrix
    call self%invert_mass( in, phi, self%s_deg_0 )

    ! Compute the degrees of freedom of the electric field as -G phi
    call complex_product_real( self%n_cells, self%eig_d, phi, efield )
    efield = -efield
  
  end subroutine compute_phi_from_rho_1d_ps

  
  !> For model with adiabatic electrons
  subroutine compute_phi_from_j_1d_ps( self, in, phi, efield )
    class(sll_t_maxwell_1d_ps) :: self
    sll_real64, intent(in)                     :: in(:)
    sll_real64, intent(out)                    :: phi(:)
    sll_real64, intent(out)                    :: efield(:)

    ! Compute divergence of the current (G^T current) (assuming a 1v model)
    self%wsave = in
    call complex_product_real( self%n_cells, self%eig_d, self%wsave, self%work )
    
    ! Compute phi by inverting the mass matrix
    call self%invert_mass( self%work, self%wsave, self%s_deg_0 )
    phi = phi + self%wsave
    
    ! Compute the degrees of freedom of the electric field as -G phi
    call complex_product_real( self%n_cells, self%eig_d, phi, efield )
    efield = -efield
  
  end subroutine compute_phi_from_j_1d_ps
     
end module sll_m_maxwell_1d_ps
