!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations in 1D
!> @details
!> Contains the abstract class to create a Maxwell solver in 1D.

module sll_m_maxwell_1d_fem
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"

  use sll_constants
  use sll_m_maxwell_1d_base

  implicit none
  private
  
  public :: sll_new_maxwell_1d_fem

  type, public, extends(sll_maxwell_1d_base) :: sll_maxwell_1d_fem

     sll_real64 :: Lx          !< length of Periodic domain
     sll_real64 :: delta_x     !< cell size
     sll_int32  :: n_dofs      !< number of cells (and grid points)
     sll_int32  :: s_deg_0     !< spline degree 0-forms
     sll_int32  :: s_deg_1     !< spline degree 1-forms
     sll_real64, allocatable :: mass_0(:)      !< coefficients of 0-form mass matrix
     sll_real64, allocatable :: mass_1(:)      !< coefficients of 1-form mass matrix
     sll_real64, allocatable :: eigenvalues(:)  !< eigenvalues of circulant update matrix
     sll_real64, dimension(:), pointer :: wsave !< array used by fftpack
     sll_real64, dimension(:), pointer :: work  !< array used by fftpack

   contains
     procedure :: &
          compute_E_from_B => compute_E_from_B_1d_fem!< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => compute_B_from_E_1d_fem!< Solve Faraday equation with E constant in time
     procedure :: &
          compute_E_from_rho => compute_E_from_rho_1d_fem!< Solve E from rho using Poisson
  end type sll_maxwell_1d_fem

contains

  subroutine compute_E_from_B_1d_fem(this, delta_t, field_in, field_out)
    class(sll_maxwell_1d_fem) :: this
    sll_real64, intent(in)     :: delta_t
    sll_real64, intent(in)     :: field_in(:)
    sll_real64, intent(inout)  :: field_out(:)
    ! local variables
    sll_int32 :: i
    sll_real64 :: coef

    coef = delta_t/this%delta_x
    ! relation betwen spline coefficients for strong Ampere
    do i=2,this%n_dofs
       field_out(i) = field_out(i) + coef * ( field_in(i-1) - field_in(i) )
    end do
    ! treat Periodic point
    field_out(1) = field_out(1) + coef * ( field_in(this%n_dofs-1) - field_in(1) )
  end subroutine compute_E_from_B_1d_fem

   subroutine compute_B_from_E_1d_fem(this, delta_t, field_in, field_out)
    class(sll_maxwell_1d_fem)  :: this
    sll_real64, intent(in)     :: delta_t
    sll_real64, intent(in)     :: field_in(:)  ! ey
    sll_real64, intent(inout)  :: field_out(:) ! bz 
    ! local variables
    sll_real64 :: coef
    
    this%work = field_in
    ! Forward FFT
    call dfftf( this%n_dofs, this%work, this%wsave)
    ! multiply by eigenvalue
    this%work(:) = this%work(:) * this%eigenvalues(:)
    ! Backward FFT 
    call dfftb( this%n_dofs, this%work,  this%wsave )
    ! Update bz from this value
    coef = delta_t/this%delta_x
    field_out(:) =  field_out(:) + coef*field_in(:)
   end subroutine compute_B_from_E_1d_fem

  
   subroutine compute_E_from_rho_1d_fem(this, E, rho )       
     class(sll_maxwell_1d_fem) :: this
     sll_real64,dimension(:),intent(in) :: rho
     sll_real64,dimension(:),intent(out) :: E

     E = 0.0_f64

   end subroutine compute_E_from_rho_1d_fem

   function sll_new_maxwell_1d_fem(Lx, n_dofs, s_deg_0) result(this)
     sll_real64 :: Lx     ! length of periodic domain
     sll_int32 :: n_dofs  ! number of degrees of freedom (here number of cells and grid points)
     sll_real64 :: delta_x ! cell size
     sll_int32 :: s_deg_0 ! highest spline degree
     class(sll_maxwell_1d_fem), pointer :: this

     ! local variables
     sll_int32 :: ierr
     sll_int32 :: j, k ! loop variables
     sll_real64 :: coef0, coef1 

     SLL_ALLOCATE(this, ierr)

     this%n_dofs = n_dofs
     this%Lx = Lx
     this%delta_x = Lx / n_dofs
     this%s_deg_0 = s_deg_0
     this%s_deg_1 = s_deg_0 - 1

     SLL_ALLOCATE(this%mass_0(s_deg_0+1), ierr)
     SLL_ALLOCATE(this%mass_1(s_deg_0), ierr)

     select case(s_deg_0)
        case(3)
           ! Upper diagonal coeficients  of cubic spline mass matrix (from Eulerian numbers)
           this%mass_0(1) = 2416.0_f64/5040.0_f64 
           this%mass_0(2) = 1191.0_f64/5040.0_f64
           this%mass_0(3) = 120.0_f64/5040.0_f64
           this%mass_0(4) = 1.0_f64/5040.0_f64
           ! Upper diagonal coeficients  of quadratic spline mass matrix (from Eulerian numbers)
           this%mass_0(1) = 66.0_f64/120.0_f64 
           this%mass_0(2) = 26.0_f64/120.0_f64
           this%mass_0(3) = 1.0_f64/120.0_f64
           
        case default
           print*, 'sll_new_maxwell_1d_fem: spline degree ', s_deg_0, ' not implemented'
     end select

     SLL_ALLOCATE(this%eigenvalues(n_dofs), ierr)
     ! Initialise FFT
     SLL_ALLOCATE(this%wsave(2*this%n_dofs+15),ierr)
     call dffti(this%n_dofs,this%wsave)
     SLL_ALLOCATE(this%work(n_dofs),ierr)

     ! Compute eigenvalues of circulant update matrix M_0^{-1} D^T M_1
     ! zero mode vanishes due to derivative matrix D^T
     this%eigenvalues(1) = 0.0_f64 
     do k=1, n_dofs/2 - 1
        coef0 =  this%mass_0(1)
        coef1 =  this%mass_1(1)
    
        do j=1,s_deg_0 - 1
           coef0 = coef0 + 2* this%mass_0(j+1)*cos(2*sll_pi*j*k/n_dofs)
           coef1 = coef1 + 2* this%mass_1(j+1)*cos(2*sll_pi*j*k/n_dofs)
        enddo
        ! add last term for larger matrix
        j = s_deg_0
        coef0 = coef0 + 2* this%mass_0(j+1)*cos(2*sll_pi*j*k/n_dofs)
        ! compute eigenvalues
        this%eigenvalues(2*k) =  (coef1 / coef0) * (1+cos(2*sll_pi*k/n_dofs)) ! real part
        this%eigenvalues(2*k+1) =  (coef1 / coef0) * sin(2*sll_pi*k/n_dofs) ! imaginary part
     enddo
     ! N/2 mode
     coef0 =  this%mass_0(1)
     coef1 =  this%mass_1(1)
    
     do j=1,s_deg_0 - 1
        coef0 = coef0 + 2 * this%mass_0(j+1)*cos(sll_pi*j)
        coef1 = coef1 + 2 * this%mass_1(j+1)*cos(sll_pi*j)
     enddo
     ! add last term for larger matrix
     j = s_deg_0
     coef0 = coef0 + 2 * this%mass_0(j+1)*cos(sll_pi*j)
     ! compute eigenvalues
     this%eigenvalues(n_dofs) = 2.0_f64 * (coef1 / coef0)
   end function sll_new_maxwell_1d_fem


end module sll_m_maxwell_1d_fem
