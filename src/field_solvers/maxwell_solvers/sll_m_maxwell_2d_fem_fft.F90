!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations in 2D
!> The linear systems are solved based on FFT diagnoalization
!> @details
!> 
!> @author
!> Katharina Kormann


! TODO: Write FFT-based mass solver: There is such a solver already defined as linear_solver_mass1 in particle_methods. Reuse? Can we put the parallelization in this solver?
! Remove all parts that belong the PLF
! Add also solver for combined e and b (first step for AVF-based algorithm)


module sll_m_maxwell_2d_fem_fft
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_low_level_bsplines, only: &
    sll_s_uniform_bsplines_eval_basis

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_poisson_2d_fem_fft, only : &
       sll_t_poisson_2d_fem_fft

  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line, &
       sll_s_spline_fem_mixedmass_line, &
       sll_s_spline_fem_multiply_mass, &
       sll_s_spline_fem_multiply_massmixed, &
       sll_s_spline_fem_compute_mass_eig
  
  use sll_m_linear_solver_spline_mass_2d_fft, only : &
       sll_t_linear_solver_spline_mass_2d_fft
  
  implicit none

  public :: &
       sll_t_maxwell_2d_fem_fft

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type :: sll_t_maxwell_2d_fem_fft

     sll_real64 :: Lx(2)          !< length of Periodic domain
     sll_real64 :: delta_x(2)     !< cell size
     sll_real64 :: volume         !< product(delta_x)
     sll_int32  :: n_dofs(2)      !< number of cells (and grid points)
     sll_int32  :: n_total        !< total number of cells
     sll_int32  :: s_deg_0        !< spline degree 0-forms
     sll_int32  :: s_deg_1        !< spline degree 1-forms
     sll_real64, allocatable :: wsave(:) !< scratch data
     sll_real64, allocatable :: work(:)  !< scratch data
     sll_real64, allocatable :: work2d(:,:)  !< scratch data
     sll_real64, allocatable :: work_d1(:)
     sll_real64, allocatable :: work_d2_in(:)
     sll_real64, allocatable :: work_d2_out(:)
     sll_real64, allocatable :: work2(:)  !< scratch data
     sll_real64, allocatable :: mass_line_0(:,:)
     sll_real64, allocatable :: mass_line_1(:,:)
     sll_real64, allocatable :: mass_line_mixed(:,:)
     type(sll_t_linear_solver_spline_mass_2d_fft) :: inverse_mass_1(3)
     type(sll_t_linear_solver_spline_mass_2d_fft) :: inverse_mass_2(3)
     type(sll_t_poisson_2d_fem_fft ) :: poisson_fft

   contains
     procedure :: &
          compute_E_from_B => sll_s_compute_e_from_b_2d_fem!< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => sll_s_compute_b_from_e_2d_fem!< Solve Faraday equation with E constant in time
     procedure :: &
          compute_E_from_rho => sll_s_compute_e_from_rho_2d_fem!< Solve E from rho using Poisson
     procedure :: &
          compute_E_from_j => compute_E_from_j_2d_fem !< Solve E from j 
     procedure :: &
          compute_rhs_from_function => sll_s_compute_fem_rhs
     procedure :: &
          L2norm_squared => L2norm_squared_2d_fem
     procedure :: &
          inner_product => inner_product_2d_fem
     procedure :: &
          L2projection => L2projection_2d_fem
     procedure :: &
          free => free_2d_fem
     procedure :: &
          init => init_2d_fem
     procedure :: &
          compute_rho_from_e => compute_rho_from_e_2d_fem

     procedure :: &
          multiply_ct
     procedure :: &
          multiply_c
     procedure :: &
          multiply_mass => multiply_mass_all

  end type sll_t_maxwell_2d_fem_fft
  
!---------------------------------------------------------------------------!
  abstract interface
     !> 2d real function
     function sll_i_function_2d_real64(x)
       use sll_m_working_precision ! can't pass a header file because the
                                 ! preprocessor prevents double inclusion.
                                 ! It is very rare.
       sll_real64             :: sll_i_function_2d_real64
       sll_real64, intent(in) :: x(2)
     end function sll_i_function_2d_real64
  end interface
contains

  !> compute rho from e using weak Gauss law ( rho = G^T M_1 e ) 
  subroutine compute_rho_from_e_2d_fem(self, efield, rho)
    class(sll_t_maxwell_2d_fem_fft) :: self
    sll_real64, intent(in)     :: efield(:)  !< efield
    sll_real64, intent(inout)  :: rho(:)  !< rho
    

    call multiply_mass_1form( self, efield, self%work )
    

    call multiply_gt( self, self%work, rho )
    rho = - rho
    

  end subroutine compute_rho_from_e_2d_fem


  
  !> compute Ey from Bz using weak Ampere formulation 
  subroutine sll_s_compute_e_from_b_2d_fem(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_2d_fem_fft) :: self
    sll_real64, intent(in)     :: delta_t   !< Time step
    sll_real64, intent(in)     :: field_in(:)  !< B
    sll_real64, intent(inout)  :: field_out(:)  !< E
    ! local variables
    sll_real64 :: coef
    sll_int32  :: comp, istart, iend

    call multiply_mass_2form( self, field_in, self%work )
    
    call multiply_ct(self, self%work, self%work2)
    
    do comp=1,3
       istart = 1+(comp-1)*self%n_total
       iend =  comp*self%n_total
       call self%inverse_mass_1(comp)%solve( self%work2(istart:iend), self%work(istart:iend) ) 
    end do
    ! Update b from self value
    coef = delta_t
    field_out = field_out + coef*self%work

  end subroutine sll_s_compute_e_from_b_2d_fem

  !> Compute Bz from Ey using strong 1D Faraday equation for spline coefficients
  !> $B_z^{new}(x_j) = B_z^{old}(x_j) - \frac{\Delta t}{\Delta x} (E_y(x_j) - E_y(x_{j-1})  $
   subroutine sll_s_compute_b_from_e_2d_fem(self, delta_t, field_in, field_out)
    class(sll_t_maxwell_2d_fem_fft)  :: self
    sll_real64, intent(in)     :: delta_t
    sll_real64, intent(in)     :: field_in(:)  ! ey
    sll_real64, intent(inout)  :: field_out(:) ! bz 

    call multiply_c(self, field_in, self%work)

    field_out = field_out - delta_t * self%work
  

  end subroutine sll_s_compute_b_from_e_2d_fem
  

   !> Compute E_i from j_i integrated over the time interval using weak Ampere formulation
   subroutine compute_E_from_j_2d_fem(self, current, component, E)
     class(sll_t_maxwell_2d_fem_fft)             :: self !< Maxwell solver class
     sll_real64,dimension(:),intent(in)    :: current !< Component \a component of the current integrated over time interval
     sll_int32, intent(in)                 :: component !< Component of the Efield to be computed
     sll_real64,dimension(:),intent(inout) :: E !< Updated electric field

     call self%inverse_mass_1(component)%solve( current, self%work(1:self%n_total) )

     E = E - self%work(1:self%n_total)

     ! Account for scaling in the mass matrix by dx
     !E = E/self%volume


   end subroutine compute_E_from_j_2d_fem
  
   subroutine sll_s_compute_e_from_rho_2d_fem(self, rho, E )       
     class(sll_t_maxwell_2d_fem_fft) :: self
     sll_real64,dimension(:),intent(in) :: rho
     sll_real64,dimension(:),intent(out) :: E


     call self%poisson_fft%compute_e_from_rho( rho, E )

   end subroutine sll_s_compute_e_from_rho_2d_fem

  !> Compute the FEM right-hand-side for a given function f and periodic splines of given degree
   !> Its components are $\int f N_i dx$ where $N_i$ is the B-spline starting at $x_i$ 
   subroutine sll_s_compute_fem_rhs(self, func , component, form, coefs_dofs)
     class(sll_t_maxwell_2d_fem_fft)             :: self
     procedure(sll_i_function_2d_real64) :: func
     sll_int32  :: component !< Specify the component
     sll_int32  :: form !< Specify 0,1,2 or 3 form
     sll_real64, intent(out) :: coefs_dofs(:)  ! Finite Element right-hand-side

     ! local variables
     sll_int32 :: i1, i2, j1, j2,  k1, k2, k, counter
     sll_int32 :: degree(2)
     sll_real64 :: coef
     sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:)
     sll_real64, allocatable :: bspl_d1(:,:), bspl_d2(:,:)

     ! Define the spline degree in the 3 dimensions, depending on form and component of the form
     if ( form == 0 ) then
        degree = self%s_deg_0
     elseif (form == 1 ) then
        degree = self%s_deg_0
        if (component<3) then
           degree(component) = self%s_deg_1
        end if
     elseif( form == 2) then
        degree =  self%s_deg_1
        if (component<3) then
           degree(component) = self%s_deg_0
        end if
     elseif( form == 3) then
        degree =  self%s_deg_1
     else 
        print*, 'Wrong form.'
     end if

     ! take enough Gauss points so that projection is exact for splines of degree deg
     ! rescale on [0,1] for compatibility with B-splines
     allocate(xw_gauss_d1(2,degree(1)+1))
     allocate(bspl_d1(degree(1)+1, degree(1)+1))
     xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights(degree(1)+1, 0.0_f64, 1.0_f64)
     ! Compute bsplines at gauss_points
     do k=1,degree(1)+1
        call sll_s_uniform_bsplines_eval_basis(degree(1),xw_gauss_d1(1,k), bspl_d1(k,:))
     end do

     allocate(xw_gauss_d2(2,degree(2)+1))
     allocate(bspl_d2(degree(2)+1, degree(2)+1))
     xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights(degree(2)+1, 0.0_f64, 1.0_f64)
     ! Compute bsplines at gauss_points
     do k=1,degree(2)+1
        call sll_s_uniform_bsplines_eval_basis(degree(2),xw_gauss_d2(1,k), bspl_d2(k,:))
     end do

     counter = 1
     ! Compute coefs_dofs = int f(x)N_i(x) 
     !do i3 = 1, self%n_dofs(3)
        do i2 = 1, self%n_dofs(2)
           do i1 = 1, self%n_dofs(1)
              coef=0.0_f64
              ! loop over support of B spline
              do j1 = 1, degree(1)+1
                 do j2 = 1, degree(2)+1
                    !do j3 = 1, degree(3)+1
                       ! loop over Gauss points
                       do k1=1, degree(1)+1
                          do k2=1, degree(2)+1
                             !do k3=1, degree(3)+1
                                coef = coef + xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* &
                                     func([self%delta_x(1)*(xw_gauss_d1(1,k1) + real(i1 + j1 - 2, f64)), &
                                     self%delta_x(2)*(xw_gauss_d2(1,k2) + real(i2 + j2 - 2, f64))] ) * &
                                     bspl_d1(k1,degree(1)+2-j1)*&
                                     bspl_d2(k2,degree(2)+2-j2)
                    
                             !enddo
                          enddo
                       end do
                    !end do
                 end do
              end do
              ! rescale by cell size
              coefs_dofs(counter) = coef*self%volume
              counter = counter+1
           enddo
        end do
     !end do

   end subroutine sll_s_compute_fem_rhs


   !> Compute the L2 projection of a given function f on periodic splines of given degree
   subroutine L2projection_2d_fem(self, func, component, form, coefs_dofs)
     class(sll_t_maxwell_2d_fem_fft) :: self
     procedure(sll_i_function_2d_real64) :: func
     sll_int32  :: component !< Specify the component
     sll_int32  :: form !< Specify 0,1,2 or 3 form
     sll_real64, intent(out) :: coefs_dofs(:)  ! spline coefficients of projection
     ! local variables
     !sll_real64 :: coef
     !sll_real64, dimension(2,degree+1) :: xw_gauss
     !sll_real64, dimension(degree+1,degree+1) :: bspl

     ! Compute right-hand-side
     call sll_s_compute_fem_rhs(self, func, component, form, self%work(1:self%n_total) )

     select case( form )
     case( 1 )
        call self%inverse_mass_1(component)%solve( self%work(1:self%n_total), coefs_dofs )
     case( 2 )
        call self%inverse_mass_2(component)%solve( self%work(1:self%n_total), coefs_dofs )
     case  default
        print*, 'L2projection for', form, '-form not implemented.'
     end select

     ! Account for scaling in the mass matrix by dx
     !coefs_dofs = coefs_dofs/self%volume

   end subroutine L2projection_2d_fem

   !> Compute square of the L2norm 
   function L2norm_squared_2d_fem(self, coefs_dofs, component, form) result (r)
     class(sll_t_maxwell_2d_fem_fft) :: self !< Maxwell solver object
     sll_real64 :: coefs_dofs(:) !< Coefficient for each DoF
     sll_int32  :: component !< Specify the component
     sll_int32  :: form !< Specify 0,1,2 or 3 form
     sll_real64 :: r !< Result: squared L2 norm


     r = inner_product_2d_fem(self, coefs_dofs, coefs_dofs, component, form)


   end function L2norm_squared_2d_fem

   !tmp: OK
   function inner_product_2d_fem(self, coefs1_dofs, coefs2_dofs, component, form) result (r)
     class(sll_t_maxwell_2d_fem_fft) :: self !< Maxwell solver object
     sll_real64 :: coefs1_dofs(:) !< Coefficient for each DoF
     sll_real64 :: coefs2_dofs(:) !< Coefficient for each DoF
     sll_int32  :: component !< Specify the component
     sll_int32  :: form !< Specify 0,1,2 or 3 form
     sll_real64 :: r !< Result: squared L2 norm


     if ( form == 0 ) then
        call multiply_mass_2dkron( self, self%mass_line_0(:,1), &
             self%mass_line_0(:,2),  &
             coefs2_dofs, self%work(1:self%n_total) )
     elseif (form == 1 ) then
        select case(component)
           case (1)
              call multiply_mass_2dkron( self, self%mass_line_1(:,1), &
                   self%mass_line_0(:,2),  &
                   coefs2_dofs, self%work(1:self%n_total) )
           case(2)
              call multiply_mass_2dkron( self, self%mass_line_0(:,1), &
                   self%mass_line_1(:,2),  &
                   coefs2_dofs, self%work(1:self%n_total) )
           case(3)
              call multiply_mass_2dkron( self, self%mass_line_0(:,1), &
                   self%mass_line_0(:,2),  &
                   coefs2_dofs, self%work(1:self%n_total) )
           case default
              print*, 'wrong component.'
           end select
     elseif( form == 2) then
        select case(component)
           case (1)
              call multiply_mass_2dkron( self, self%mass_line_0(:,1), &
                   self%mass_line_1(:,2),  &
                   coefs2_dofs, self%work(1:self%n_total) )
           case(2)
              call multiply_mass_2dkron( self, self%mass_line_1(:,1), &
                   self%mass_line_0(:,2),  &
                   coefs2_dofs, self%work(1:self%n_total) )
           case(3)
              call multiply_mass_2dkron( self, self%mass_line_1(:,1), &
                   self%mass_line_1(:,2),  &
                   coefs2_dofs, self%work(1:self%n_total) )
           case default
              print*, 'wrong component.'
           end select
     elseif( form == 3) then
        call multiply_mass_2dkron( self, self%mass_line_1(:,1), &
             self%mass_line_1(:,2),  &
             coefs2_dofs, self%work(1:self%n_total) )
     else
        print*, 'Wrong form.'
     end if

     r = sum(coefs1_dofs*self%work(1:self%n_total))

     
   end function inner_product_2d_fem
   

   subroutine init_2d_fem( self, domain, n_dofs, s_deg_0 )
     class(sll_t_maxwell_2d_fem_fft), intent(out) :: self !< solver object
     sll_real64, intent(in) :: domain(2,2)     ! xmin, xmax
     sll_int32, intent(in) :: n_dofs(2)  ! number of degrees of freedom (here number of cells and grid points)
     !sll_real64 :: delta_x ! cell size
     sll_int32, intent(in) :: s_deg_0 ! highest spline degree

    ! local variables
     sll_int32 :: j
     sll_real64 :: mass_line_0(s_deg_0+1), mass_line_1(s_deg_0), mass_line_mixed(s_deg_0*2)
     sll_real64, allocatable :: eig_values_mass_0_1(:)
     sll_real64, allocatable :: eig_values_mass_0_2(:)
     sll_real64, allocatable :: eig_values_mass_1_1(:)
     sll_real64, allocatable :: eig_values_mass_1_2(:)

     self%n_dofs = n_dofs
     self%n_total = product(n_dofs)
     
     self%Lx = domain(:,2) - domain(:,1)
     self%delta_x = self%Lx /real( n_dofs, f64 )
     self%s_deg_0 = s_deg_0
     self%s_deg_1 = s_deg_0 - 1

     self%volume = product(self%delta_x)

     ! Allocate scratch data
     allocate( self%work2d(n_dofs(1), n_dofs(2)) )
     allocate( self%work(self%n_total*3) )
     allocate( self%work2(self%n_total*3) )
     allocate( self%work_d1( n_dofs(1) ) ) 
     allocate( self%work_d2_in( n_dofs(2) ) ) 
     allocate( self%work_d2_out( n_dofs(2) ) ) 
   
 
     ! Sparse matrices
     ! Assemble the mass matrices
     ! First assemble a mass line for both degrees
     allocate( self%mass_line_0(s_deg_0+1,2) )
     allocate( self%mass_line_1(s_deg_0,2) )
     allocate( self%mass_line_mixed(2*s_deg_0,2) )
     call sll_s_spline_fem_mass_line ( self%s_deg_0, mass_line_0 )
     call sll_s_spline_fem_mass_line ( self%s_deg_1, mass_line_1 )

     call sll_s_spline_fem_mixedmass_line ( self%s_deg_0, mass_line_mixed )

     ! Next put together the 1d parts of the 2d Kronecker product
     do j=1, 2
        self%mass_line_0(:,j) = self%delta_x(j) * mass_line_0
        self%mass_line_1(:,j) = self%delta_x(j) * mass_line_1
        self%mass_line_mixed(:,j) = self%delta_x(j) * mass_line_mixed
     end do

     allocate( eig_values_mass_0_1( n_dofs(1) ) )
     allocate( eig_values_mass_0_2( n_dofs(2) ) )
     allocate( eig_values_mass_1_1( n_dofs(1) ) )
     allocate( eig_values_mass_1_2( n_dofs(2) ) )
     call sll_s_spline_fem_compute_mass_eig( n_dofs(1), s_deg_0, mass_line_0*self%delta_x(1), &
          eig_values_mass_0_1 )
     call sll_s_spline_fem_compute_mass_eig( n_dofs(2), s_deg_0, mass_line_0*self%delta_x(2), &
          eig_values_mass_0_2 )
     call sll_s_spline_fem_compute_mass_eig( n_dofs(1), s_deg_0-1, mass_line_1*self%delta_x(1), &
          eig_values_mass_1_1 )
     call sll_s_spline_fem_compute_mass_eig( n_dofs(2), s_deg_0-1, mass_line_1*self%delta_x(2), &
          eig_values_mass_1_2 )

     call self%inverse_mass_1(1)%create( n_dofs, eig_values_mass_1_1, eig_values_mass_0_2 )
     call self%inverse_mass_1(2)%create( n_dofs, eig_values_mass_0_1, eig_values_mass_1_2 )
     call self%inverse_mass_1(3)%create( n_dofs, eig_values_mass_0_1, eig_values_mass_0_2 )

     call self%inverse_mass_2(1)%create( n_dofs, eig_values_mass_0_1, eig_values_mass_1_2 )
     call self%inverse_mass_2(2)%create( n_dofs, eig_values_mass_1_1, eig_values_mass_0_2 )
     call self%inverse_mass_2(3)%create( n_dofs, eig_values_mass_1_1, eig_values_mass_1_2 )


     ! Poisson solver based on fft inversion
     call self%poisson_fft%init( self%n_dofs, self%s_deg_0, self%delta_x )

     
   end subroutine init_2d_fem


   subroutine free_2d_fem(self)
     class(sll_t_maxwell_2d_fem_fft) :: self

    
     deallocate(self%work)

   end subroutine free_2d_fem

  
   !tmp:OK
   !> Multiply by discrete curl matrix
   subroutine multiply_c(self, field_in, field_out)
    class(sll_t_maxwell_2d_fem_fft)  :: self
    sll_real64, intent(in)     :: field_in(:)  
    sll_real64, intent(inout)  :: field_out(:)  
     ! local variables
    sll_real64 :: coef(2)
    sll_int32 :: stride(2), jump(2), indp(2)
    sll_int32 :: i,j, ind2d, ind2d_1, ind2d_2

   
    ! TODO: Avoid the IF for periodic boundaries
    ! First component
    coef(1) = 1.0_f64/ self%delta_x(2)
    !coef(2) = -1.0_f64/ self%delta_x(3)

    stride(1) = self%n_dofs(1)
    stride(2) = self%n_dofs(1)*self%n_dofs(2)

    !jump(1) = self%n_total
    jump(2) = 2*self%n_total

    
    ind2d = 0
    do j=1,self%n_dofs(2)
       if ( j==1 ) then
          indp(1) = stride(1)*(self%n_dofs(2)-1)
       else
          indp(1) = - stride(1)
       end if
       do i=1,self%n_dofs(1)
          ind2d = ind2d + 1
          
          ind2d_1 = ind2d +indp(1)+jump(2)
          field_out(ind2d) =  &
               coef(1) * ( field_in( ind2d+jump(2) ) -&
               field_in( ind2d_1 ))
          
       end do
    end do
    

    ! Second component
    !coef(1) = 1.0_f64/ self%delta_x(3)
    coef(2) = -1.0_f64/ self%delta_x(1)

    stride(2) = 1
    !stride(1) = self%n_dofs(1)*self%n_dofs(2)

    jump(1) = self%n_total
    !jump(2) = -self%n_total

    
    do j=1,self%n_dofs(2)
       do i=1,self%n_dofs(1)
          if ( i==1 ) then
             indp(2) = stride(2)*(self%n_dofs(1)-1)
          else
             indp(2) = - stride(2)
          end if
          ind2d = ind2d + 1

          ind2d_2 = ind2d +indp(2)+jump(1)
             
          field_out(ind2d) = &
               coef(2) * ( field_in(ind2d+jump(1) ) - &
               field_in( ind2d_2 ))
          
       end do
    end do

    ! Third component
    coef(1) = 1.0_f64/ self%delta_x(1)
    coef(2) = -1.0_f64/ self%delta_x(2)

    stride(1) = 1
    stride(2) = self%n_dofs(1)

    jump(1) = -2*self%n_total
    jump(2) = -self%n_total

    
    do j=1,self%n_dofs(2)
       if (j == 1) then
          indp(2) = stride(2)*(self%n_dofs(2)-1)
       else
          indp(2) = - stride(2)
       end if
       do i=1,self%n_dofs(1)
          if ( i==1 ) then
             indp(1) = stride(1)*(self%n_dofs(1)-1)
          else
             indp(1) = - stride(1)
          end if
          ind2d = ind2d + 1
          
          ind2d_1 = ind2d +indp(1)+jump(2)
          ind2d_2 = ind2d +indp(2)+jump(1)
          
          field_out(ind2d) =  &
               coef(1) * ( field_in( ind2d+jump(2) ) -&
               field_in( ind2d_1 ) )+ &
               coef(2) * ( field_in(ind2d+jump(1) ) - &
               field_in( ind2d_2 ) )
          
       end do
    end do
    
  end subroutine multiply_c
  
  !tmp:OK
  !> Multiply by transpose of discrete curl matrix
  subroutine multiply_ct(self, field_in, field_out)
    class(sll_t_maxwell_2d_fem_fft)  :: self !< Maxwell solver object
    sll_real64, intent(in)     :: field_in(:)  !< Matrix to be multiplied
    sll_real64, intent(inout)  :: field_out(:)  !< C*field_in
     ! local variables
    sll_real64 :: coef(2)
    sll_int32 :: stride(2), jump(2), indp(2)
    sll_int32 :: i,j, ind2d, ind2d_1, ind2d_2

   
    ! TODO: Avoid the IF for periodic boundaries
    ! First component
    coef(1) = -1.0_f64/ self%delta_x(2)
    !coef(2) = 1.0_f64/ self%delta_x(3)

    stride(1) = self%n_dofs(1)
    !stride(2) = self%n_dofs(1)*self%n_dofs(2)

    !jump(1) = self%n_total
    jump(2) = 2*self%n_total

    
    ind2d = 0
    do j=1,self%n_dofs(2)
       if ( j== self%n_dofs(2)) then
          indp(1) = -stride(1)*(self%n_dofs(2)-1)
       else
          indp(1) = stride(1)
       end if
       do i=1,self%n_dofs(1)
          ind2d = ind2d + 1

          ind2d_1 = ind2d +indp(1)+jump(2)
             
          field_out(ind2d) =  &
               coef(1) * ( field_in( ind2d+jump(2) ) -&
                  field_in( ind2d_1 ))
          
       end do
    end do


    ! Second component
    !coef(1) = -1.0_f64/ self%delta_x(3)
    coef(2) = 1.0_f64/ self%delta_x(1)

    stride(2) = 1
    !stride(1) = self%n_dofs(1)*self%n_dofs(2)

    jump(1) = self%n_total
    !jump(2) = -self%n_total

    
    do j=1,self%n_dofs(2)
       do i=1,self%n_dofs(1)
          if ( i==self%n_dofs(1) ) then
             indp(2) = -stride(2)*(self%n_dofs(1)-1)
          else
             indp(2) = stride(2)
          end if
          ind2d = ind2d + 1
          ind2d_2 = ind2d +indp(2)+jump(1)
             
          field_out(ind2d) = &
               coef(2) * ( field_in(ind2d+jump(1) ) - &
               field_in( ind2d_2 ))
             
       end do
    end do

    ! Third component
    coef(1) = -1.0_f64/ self%delta_x(1)
    coef(2) = 1.0_f64/ self%delta_x(2)

    stride(1) = 1
    stride(2) = self%n_dofs(1)

    jump(1) = -2*self%n_total
    jump(2) = -self%n_total

    
    do j=1,self%n_dofs(2)
       if (j == self%n_dofs(2)) then
          indp(2) = -stride(2)*(self%n_dofs(2)-1)
       else
          indp(2) = stride(2)
       end if
       do i=1,self%n_dofs(1)
          if ( i==self%n_dofs(1) ) then
             indp(1) = -stride(1)*(self%n_dofs(1)-1)
          else
             indp(1) = stride(1)
          end if
          ind2d = ind2d + 1

          ind2d_1 = ind2d +indp(1)+jump(2)
          ind2d_2 = ind2d +indp(2)+jump(1)
             
          field_out(ind2d) =  &
               coef(1) * ( field_in( ind2d+jump(2) ) -&
               field_in( ind2d_1 ) )+ &
               coef(2) * ( field_in(ind2d+jump(1) ) - &
               field_in( ind2d_2 ) )
             
       end do
    end do

  end subroutine multiply_ct


  !> Multiply by transpose of dicrete gradient matrix
  subroutine multiply_gt(self, field_in, field_out)
    class(sll_t_maxwell_2d_fem_fft)  :: self !< Maxwell solver object
    sll_real64, intent(in)     :: field_in(:)  !< Matrix to be multiplied
    sll_real64, intent(inout)  :: field_out(:)  !< C*field_in

    sll_real64 :: coef
    sll_int32  :: jump, jump_end
    sll_int32  :: ind2d, ind2d_1
    sll_int32  :: i,j


    coef = 1.0_f64/ self%delta_x(1)
    jump = 1
    jump_end = 1-self%n_dofs(1)

    
    ind2d = 0
    do j=1,self%n_dofs(2)
       do i=1,self%n_dofs(1)-1
          ind2d = ind2d + 1
             
          field_out(ind2d) =  &
               coef * ( field_in(ind2d) - field_in(ind2d+jump) )
             
       end do
       ind2d = ind2d + 1
       field_out(ind2d) = coef * ( field_in(ind2d) - field_in(ind2d+jump_end) )
    end do

    coef = 1.0_f64/ self%delta_x(2)
    jump = self%n_dofs(1)
    jump_end = (1-self%n_dofs(2))*self%n_dofs(1)

    ind2d_1 = 0
    do j=1,self%n_dofs(2)-1
       do i=1,self%n_dofs(1)
          ind2d = ind2d + 1
          ind2d_1 = ind2d_1 + 1
             
          field_out(ind2d_1) =  field_out(ind2d_1) + &
               coef * ( field_in(ind2d) - field_in(ind2d+jump) )
       end do
    end do
    do i=1,self%n_dofs(1)
       ind2d = ind2d + 1
       ind2d_1 = ind2d_1 + 1
          
       field_out(ind2d_1) =  field_out(ind2d_1) + &
               coef * ( field_in(ind2d) - field_in(ind2d+jump_end) )
          
    end do
    

  end subroutine multiply_gt

  !tmp:OK
  subroutine multiply_mass_1form( self, coefs_in, coefs_out )
    class(sll_t_maxwell_2d_fem_fft), intent( inout )  :: self
    sll_real64, intent(in)     :: coefs_in(:)
    sll_real64, intent(out)  :: coefs_out(:)

    sll_int32:: iend, istart

    istart = 1
    iend = self%n_total
    call multiply_mass_2dkron(  self, self%mass_line_1(:,1), &
         self%mass_line_0(:,2), &
         coefs_in(istart:iend), coefs_out(istart:iend) )
    istart = iend+1
    iend = iend + self%n_total
    call multiply_mass_2dkron(  self, self%mass_line_0(:,1), &
         self%mass_line_1(:,2),  &
         coefs_in(istart:iend), coefs_out(istart:iend) )
    istart = iend+1
    iend = iend + self%n_total
    call multiply_mass_2dkron(  self, self%mass_line_0(:,1), &
         self%mass_line_0(:,2), &
         coefs_in(istart:iend), coefs_out(istart:iend) )

    
  end subroutine multiply_mass_1form
  
  !tmp:OK
   subroutine multiply_mass_2form( self, coefs_in, coefs_out )
    class(sll_t_maxwell_2d_fem_fft), intent( inout )  :: self
    sll_real64, intent(in)     :: coefs_in(:)
    sll_real64, intent(out)  :: coefs_out(:)

    sll_int32:: iend, istart

    istart = 1
    iend = self%n_total
    call multiply_mass_2dkron(  self, self%mass_line_0(:,1), &
         self%mass_line_1(:,2),  &
         coefs_in(istart:iend), coefs_out(istart:iend) )
    istart = iend+1
    iend = iend + self%n_total
    call multiply_mass_2dkron(  self, self%mass_line_1(:,1), &
         self%mass_line_0(:,2),  &
         coefs_in(istart:iend), coefs_out(istart:iend) )
    istart = iend+1
    iend = iend + self%n_total
    call multiply_mass_2dkron(  self, self%mass_line_1(:,1), &
         self%mass_line_1(:,2),  &
         coefs_in(istart:iend), coefs_out(istart:iend) )

    
  end subroutine multiply_mass_2form
  

  !tmp:OK
  !> Multiply by the mass matrix 
  subroutine multiply_mass_2dkron(  self, mass_line_1, mass_line_2,  coefs_in, coefs_out )
    class(sll_t_maxwell_2d_fem_fft), intent( inout )  :: self
    sll_real64, intent(in)    :: mass_line_1(:)
    sll_real64, intent(in)    :: mass_line_2(:)
    sll_real64, intent(in)     :: coefs_in(:)
    sll_real64, intent(out)  :: coefs_out(:)  

    ! Local variables
    sll_int32 :: i,j,istart,iend
    sll_int32 :: deg(2)

    deg(1) = size(mass_line_1)-1
    deg(2) = size(mass_line_2)-1
    
    istart = 1
    iend = self%n_dofs(1)
    do j=1,self%n_dofs(2)
       call sll_s_spline_fem_multiply_mass ( self%n_dofs(1), deg(1), &
            mass_line_1, coefs_in(istart:iend), self%work_d1 )           
       self%work2d(:,j) = self%work_d1
       istart = iend+1
       iend = iend + self%n_dofs(1)
    end do

    istart = 1
    do i =1,self%n_dofs(1)
       self%work_d2_in = self%work2d(i,:)
       call sll_s_spline_fem_multiply_mass ( self%n_dofs(2), deg(2), &
            mass_line_2, self%work_d2_in, self%work_d2_out )
       do j=1,self%n_dofs(2)
          coefs_out(istart+(j-1)*self%n_dofs(1)) = self%work_d2_out(j)
       end do
       istart = istart+1
    end do
     
  end subroutine multiply_mass_2dkron

  !tmp:OK
  !> Multiply by the mass matrix 
  subroutine multiply_mass_all(  self, deg, coefs_in, coefs_out )
    class(sll_t_maxwell_2d_fem_fft), intent( inout )  :: self
    sll_int32, intent( in ) :: deg(2) !< \a deg(i) specifies the degree of the 1d mass matrix in dimension \a i (Note: 1 for 0-form, 2 for 1-form, 3 for 0-1-form mix)
    sll_real64, intent(in)   :: coefs_in(:)
    sll_real64, intent(out)  :: coefs_out(:)  

    ! Local variables
    sll_int32 :: i,j,istart,iend
    
    istart = 1
    iend = self%n_dofs(1)
    do j=1,self%n_dofs(2)
       select case ( deg(1) )
       case( 1 )
          call sll_s_spline_fem_multiply_mass ( self%n_dofs(1), self%s_deg_0, &
               self%mass_line_0(:,1), coefs_in(istart:iend), self%work_d1 )
       case( 2 )
          call sll_s_spline_fem_multiply_mass ( self%n_dofs(1), self%s_deg_1, &
               self%mass_line_1(:,1), coefs_in(istart:iend), self%work_d1 )
       case ( 3 )
          call sll_s_spline_fem_multiply_massmixed ( self%n_dofs(1), self%s_deg_0, &
               self%mass_line_mixed(:,1), coefs_in(istart:iend), self%work_d1 )
       case   default
          print*, 'not implemented.'
       end select
       self%work2d(:,j) = self%work_d1
       istart = iend+1
       iend = iend + self%n_dofs(1)
    end do

    istart = 1
    do i =1,self%n_dofs(1)
       self%work_d2_in = self%work2d(i,:)
       select case ( deg(2) )
       case( 1 )
          call sll_s_spline_fem_multiply_mass ( self%n_dofs(2), self%s_deg_0, &
               self%mass_line_0(:,2), self%work_d2_in, self%work_d2_out )
       case( 2 )
          call sll_s_spline_fem_multiply_mass ( self%n_dofs(2), self%s_deg_1, &
               self%mass_line_1(:,2), self%work_d2_in, self%work_d2_out )
       case ( 3 )
          call sll_s_spline_fem_multiply_massmixed ( self%n_dofs(2), self%s_deg_0, &
               self%mass_line_mixed(:,2), self%work_d2_in, self%work_d2_out )
       case   default
          print*, 'not implemented.'
       end select

           
       do j=1,self%n_dofs(2)
          coefs_out(istart+(j-1)*self%n_dofs(1)) = self%work_d2_out(j)
       end do
       istart = istart+1
    end do
     
  end subroutine multiply_mass_all

 end module sll_m_maxwell_2d_fem_fft
