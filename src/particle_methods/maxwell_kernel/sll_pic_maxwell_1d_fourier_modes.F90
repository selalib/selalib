!> @ingroup particle_methods
!> @brief
!> Interface for combined Maxwell solver and kernel smoother for PIC methods
!> @details
!> Derived types can be implemented directly or initialized using a Maxwell solver and a kernel smoother.

module sll_m_pic_maxwell_1d_fourier_modes
#include "sll_working_precision.h"
#include "sll_memory.h"

  use sll_m_pic_base
  use sll_m_pic_maxwell_base
  use sll_constants, only: sll_pi

  implicit none
  !private
  
  type, public, extends(sll_pic_maxwell_base) :: sll_pic_maxwell_1d_fourier_modes
     sll_int32  :: n_fmodes
     sll_real64, allocatable :: kEB(:) !(this%n_fmodes) ! (1:n_fmodes) * 2*sll_pi/L

     sll_real64, allocatable :: shape_factors(:,:)

     sll_int32 :: n_particles

   contains
     ! Maxwell solver functions
     procedure :: &
          compute_E_from_B => compute_E_from_B_1d_fm!< Solve E and B part of Amperes law with B constant in time
     procedure :: &
          compute_B_from_E => compute_E_from_B_1d_fm!< Solve Faraday equation with E constant in time
     procedure :: &
          compute_E_from_rho => compute_E_from_rho_1d_fm!< Solve E from rho using Poisson
     !procedure(signature_solve), deferred :: &
     !     solve !< Solve Amperes law and Faraday equation
     ! Kernel smoother functions
     procedure          :: &
          compute_shape_factors => compute_shape_factors_1d_fm!< Prepare for the accumulation by computing the shape factors
     procedure          :: &
          accumulate_rho_from_klimontovich => accumulate_rho_1d_fm!< Accumulate the charge density
     procedure          :: &
          accumulate_j_from_klimontovich => accumulate_j_1d_fm!< Accumulate the current density
     procedure          :: evaluate_kernel_function => evaluate_kernel_function_1d_fm!< Evaluate a function based on the given degrees of freedom
     

  end type sll_pic_maxwell_1d_fourier_modes

contains

  subroutine compute_E_from_B_1d_fm(this, delta_t, field_in, field_out)    
     class(sll_pic_maxwell_1d_fourier_modes) :: this
     sll_real64, intent(in)     :: delta_t
     sll_real64, intent(in)     :: field_in(:)
     sll_real64, intent(inout)  :: field_out(:)

     field_out(2:this%n_fmodes+1) = field_out(2:this%n_fmodes+1) - &
          delta_t * this%kEB * field_in(this%n_fmodes+2:2*this%n_fmodes+1)
     field_out(this%n_fmodes+2:2*this%n_fmodes+1) = &
          field_out(this%n_fmodes+2:2*this%n_fmodes+1) + &
          delta_t * this%kEB * field_in(2:this%n_fmodes+1)
    

   end subroutine compute_E_from_B_1d_fm


   subroutine compute_E_from_rho_1d_fm(this, E, rho )       
     class(sll_pic_maxwell_1d_fourier_modes) :: this
     sll_real64,dimension(:),intent(in) :: rho
     sll_real64,dimension(:),intent(out) :: E

     E(2:this%n_fmodes+1) = - rho(this%n_fmodes+2:2*this%n_fmodes+1)/this%kEB
     E(this%n_fmodes+2:2*this%n_fmodes+1) = rho(2:this%n_fmodes+1)/this%kEB

     E(1) = 0.0_f64
     
   end subroutine compute_E_from_rho_1d_fm
  


  ! Interfaces for kernel smoother
  !---------------------------------------------------------------------------!
   subroutine compute_shape_factors_1d_fm(this, particle_group)
     class( sll_pic_maxwell_1d_fourier_modes), intent(inout) :: this !< Kernel smoother object.
     class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.

     sll_int32 :: i_part, j
     sll_real64 :: x(3)
     
     do i_part = 1, this%n_particles
        x = particle_group%get_x(i_part)
        do j=1, this%n_fmodes
           this%shape_factors(j, i_part) = cos(this%kEB(j)*x(1))
           this%shape_factors(j+this%n_fmodes, i_part) = sin(this%kEB(j)*x(1))
        end do
     end do


   end subroutine compute_shape_factors_1d_fm
  
!---------------------------------------------------------------------------!
  subroutine accumulate_rho_1d_fm(this, particle_group, rho_dofs)       
    class( sll_pic_maxwell_1d_fourier_modes), intent(in)    :: this !< Kernel smoother object.
    class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
    sll_real64, intent(inout)                       :: rho_dofs(:) !< Degrees of freedom in kernel representation (can be point values or weights in a basis function representation).


    sll_int32 :: i_part, j
    sll_real64 :: wq
     
    do i_part = 1, this%n_particles      
        wq =  particle_group%get_charge(i_part)
        rho_dofs(1) = rho_dofs(1) + wq
        do j=1, this%n_fmodes*2
           rho_dofs(j+1) = rho_dofs(j+1) + wq * this%shape_factors(j, i_part) 
        end do
     end do
    

  end subroutine accumulate_rho_1d_fm

!---------------------------------------------------------------------------!
  subroutine accumulate_j_1d_fm(this, &
          particle_group, &
          j_dofs, &
          component)       
    class( sll_pic_maxwell_1d_fourier_modes), intent(in)    :: this !< Kernel smoother object.
    class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
    sll_real64, intent(inout)                       :: j_dofs(:)!< Degrees of freedom in kernel representation (can be point values or weights in a basis function representation).
    sll_int32, intent (in)                          :: component !< Component of the current density that should be evaluated.

    sll_int32 :: i_part, j
    sll_real64 :: wq
    sll_real64 :: v(3)

    do i_part = 1, this%n_particles      
        wq =  particle_group%get_charge(i_part)
        v = particle_group%get_v(i_part)
        j_dofs(1) = j_dofs(1) + wq
        do j=1, this%n_fmodes
           j_dofs(j+1) = j_dofs(j+1) + &
                wq * this%shape_factors(j, i_part) * v(component) / &
                (this%kEB(j)*v(1))
           j_dofs(j+1+this%n_fmodes) = j_dofs(j+1+this%n_fmodes) + &
                wq * this%shape_factors(j+this%n_fmodes, i_part) * v(component) / &
                (this%kEB(j)*v(1))
        end do
     end do


  end subroutine accumulate_j_1d_fm

!---------------------------------------------------------------------------!
  subroutine evaluate_kernel_function_1d_fm(this, particle_group, rho_dofs, particle_values)       
    class( sll_pic_maxwell_1d_fourier_modes), intent(in)    :: this !< Kernel smoother object.
    class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
    sll_real64, intent(in)                       :: rho_dofs(:) !< Degrees of freedom in kernel representation.
    sll_real64, intent(out)                      :: particle_values(:) !< Values of the function represented by \a rho_dofs at particle positions.

    sll_int32 :: i_part, j

    particle_values = rho_dofs(1)
    do i_part = 1, particle_group%n_particles
       do j=1,this%n_fmodes
          particle_values(i_part) = particle_values(i_part) + &
               rho_dofs(j+1) * this%shape_factors(j, i_part) + &
               rho_dofs(this%n_fmodes + 1 + j) * this%shape_factors(j+ this%n_fmodes, i_part) 
       end do
    end do

  end subroutine evaluate_kernel_function_1d_fm


  !< Constructor
  function sll_new_pic_maxwell_1d_fourier_modes(n_fmodes, Lx, n_particles) result(this)
    class( sll_pic_maxwell_1d_fourier_modes), pointer :: this
    sll_int32, intent(in) :: n_fmodes
    sll_int32, intent(in) :: n_particles
    sll_real64, intent(in) :: Lx

    sll_int32 :: ierr, i

    SLL_ALLOCATE(this, ierr)

    this%n_fmodes = n_fmodes
    SLL_ALLOCATE(this%kEB(this%n_fmodes), ierr)

    this%n_dofs = this%n_fmodes*2 + 1

    do i=1,this%n_fmodes
       this%kEB(i) = 2.0_f64*sll_pi/Lx * real(i, f64)
    end do

    this%n_particles = n_particles

    SLL_ALLOCATE(this%shape_factors(this%n_fmodes*2, this%n_particles), ierr)
    

  end function sll_new_pic_maxwell_1d_fourier_modes


end module sll_m_pic_maxwell_1d_fourier_modes
