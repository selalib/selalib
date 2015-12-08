!> @ingroup kernel_smoother
!> @author Katharina Kormann, IPP
!> @brief Kernel smoother for 2d with splines of arbitrary degree placed on a uniform mesh.
!> @details Spline with index i starts at point i
module sll_m_kernel_smoother_spline_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_arbitrary_degree_splines, only: &
    uniform_b_splines_at_x

  use sll_m_kernel_smoother_base, only: &
    sll_collocation, &
    sll_kernel_smoother_base

  use sll_m_particle_group_base, only: &
    sll_particle_group_base

  implicit none

  public :: &
    sll_kernel_smoother_spline_1d, &
    sll_new_smoother_spline_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !>  Spline kernel smoother in1d.
  type, extends(sll_kernel_smoother_base) :: sll_kernel_smoother_spline_1d

     ! Information about the 1d mesh
     sll_real64 :: delta_x(1)  !< Value of grid spacing along both directions.
     sll_real64 :: domain(1,2) !< Definition of the domain: domain(1,1) = x1_min  domain(1,2) = x1_max
     
     ! Information about the particles
     sll_int32 :: no_particles !< Number of particles of underlying PIC method (processor local)

     ! 
     sll_int32 :: spline_degree !< Degree of smoothing kernel spline
     sll_int32 :: n_span !< Number of intervals where spline non zero (spline_degree + 1)
     sll_int32 :: smoothing_type !< Defines whether we use Galerkin or collocation in order to decide on scaling when accumulating

     ! Sparse structure for shape factors
     sll_int32,  allocatable :: index_grid(:,:)    !< First dof index with contribution from th
     sll_real64, allocatable :: values_grid(:,:,:) !< Values of the space factors in each dimesion.
     
   contains
     procedure :: compute_shape_factors => compute_shape_factors_spline_1d !> Compute shape factors 
     procedure :: accumulate_rho_from_klimontovich => accumulate_rho_from_klimontovich_spline_1d !> accumulate density
     procedure :: accumulate_j_from_klimontovich => accumulate_j_from_klimontovich_spline_1d !> accumulate component of current density 
     procedure :: evaluate_kernel_function_particle => evaluate_kernel_function_particle_spline_1d !> Evaluate spline function with given coefficients

  end type sll_kernel_smoother_spline_1d
  
contains
  !---------------------------------------------------------------------------!
  subroutine compute_shape_factors_spline_1d(this, particle_group)
    class( sll_kernel_smoother_spline_1d), intent(inout) :: this !< kernel smoother object
    class( sll_particle_group_base), intent(in)     :: particle_group  !< particle group

    ! local variables
    sll_real64 :: xi(3)
    sll_int32 :: i_part
    sll_real64 :: spline_val(this%n_span)

    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x(i_part)
       xi(1) = (xi(1) - this%domain(1,1)) /&
            this%delta_x(1)
       ! This is the index of the last spline the particle contributes to
       this%index_grid(:,i_part) = ceiling(xi(1))
       xi(1) = xi(1) - real(this%index_grid(1,i_part) -1,f64)
       ! Now we subtract the degree of the spline to get the index of the first spline.
       this%index_grid(:,i_part) =  this%index_grid(:,i_part) - this%spline_degree
       spline_val = uniform_b_splines_at_x(this%spline_degree, xi(1))
       this%values_grid(:,1,i_part) = spline_val
    end do

  end subroutine compute_shape_factors_spline_1d

  !---------------------------------------------------------------------------!
  subroutine accumulate_rho_from_klimontovich_spline_1d(this, particle_group,&
       rho_dofs)
    class( sll_kernel_smoother_spline_1d), intent(in)    :: this !< kernel smoother object
    class( sll_particle_group_base), intent(in)     :: particle_group !< particle group
    sll_real64, intent(inout)                       :: rho_dofs(:) !< array holding the spline coefficients accumulated rho
    
    !local variables
    sll_int32 :: i_part, i1
    sll_int32 :: index1d
    sll_real64 :: wi(1)

    do i_part = 1, particle_group%n_particles
       wi = particle_group%get_charge(i_part)
       do i1 = 1, this%n_span
          index1d = modulo(this%index_grid(1,i_part)+i1-2,this%n_grid(1))+1
          rho_dofs(index1d) = rho_dofs(index1d) +&
               (wi(1) * &
               this%values_grid(i1, 1, i_part))
       end do
    end do

    if (this%smoothing_type == SLL_COLLOCATION) then
       rho_dofs = rho_dofs/this%delta_x(1)
    end if
       
  end subroutine accumulate_rho_from_klimontovich_spline_1d

  !---------------------------------------------------------------------------!
  subroutine accumulate_j_from_klimontovich_spline_1d(this, particle_group,&
       j_dofs, component)
    class( sll_kernel_smoother_spline_1d), intent(in)    :: this !< kernel smoother object
    class( sll_particle_group_base), intent(in)     :: particle_group !< particle group
    sll_real64, intent(inout)                       :: j_dofs(:) !< array holding the spline coefficients of component \a component the accumulated j
    sll_int32, intent(in)                           :: component !< component of \a j_dof to be accumulated
    
    !local variables
    sll_int32 :: i_part, i1
    sll_int32 :: index1d
    sll_real64 :: vpart(3)
    sll_real64 :: wi(1)

    do i_part = 1, particle_group%n_particles
       wi = particle_group%get_charge(i_part)
       do i1 = 1, this%n_span
          index1d = modulo(this%index_grid(1,i_part)+i1-2,this%n_grid(1))+1
          vpart = particle_group%get_v(i_part)
          j_dofs(index1d) = j_dofs(index1d) +&
               (wi(1) * &
               vpart(component) * &
               this%values_grid(i1,1,i_part))
       end do
    end do
   
    if (this%smoothing_type == SLL_COLLOCATION) then
       j_dofs = j_dofs/this%delta_x(1)
    end if

  end subroutine accumulate_j_from_klimontovich_spline_1d
  
  !---------------------------------------------------------------------------!
  subroutine evaluate_kernel_function_particle_spline_1d(this, rho_dofs, i_part, particle_value)
    class( sll_kernel_smoother_spline_1d), intent(in)    :: this !< kernel smoother object
    sll_real64, intent(in)                       :: rho_dofs(:) !< Degrees of freedom in kernel representation.
    sll_int32, intent(in)                        :: i_part !< particle number
    sll_real64, intent(out)                      :: particle_value !< Value of the function at the position of particle \a i_part
  
    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d

    particle_value = 0.0_f64
    do i1 = 1, this%n_span
       index1d = modulo(this%index_grid(1,i_part)+i1-2, this%n_grid(1))+1
       particle_value = particle_value + &
            rho_dofs(index1d) *  &
            this%values_grid(i1,1, i_part)
    end do

  end subroutine evaluate_kernel_function_particle_spline_1d


  !-------------------------------------------------------------------------------------------
  !< Constructor 
  function sll_new_smoother_spline_1d(domain, n_grid, no_particles, spline_degree, smoothing_type) result (this)
    class( sll_kernel_smoother_spline_1d), pointer   :: this !< kernel smoother object
    sll_int32, intent(in) :: n_grid(1) !< number of DoFs (spline coefficients)
    sll_real64, intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32, intent(in) :: no_particles !< number of particles
    sll_int32, intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32, intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE( this, ierr)

    ! Store grid information
    this%domain(1,:) = domain
    SLL_ALLOCATE(this%n_grid(1), ierr)
    this%n_grid = n_grid
    this%n_dofs = product(n_grid)
    this%delta_x = (this%domain(:,2)-this%domain(:,1))/real(n_grid, f64)

    ! Store basis function information
    this%no_particles = no_particles

    ! Initialize information on the spline
    this%spline_degree = spline_degree
    this%n_span = spline_degree + 1

    ! Initialize information on smoothing type
    this%smoothing_type = smoothing_type

    ! Initialize sparse structure for shape factors
    SLL_ALLOCATE(this%index_grid(1, no_particles),ierr)
    SLL_ALLOCATE(this%values_grid(this%n_span, 1, no_particles),ierr)

  end function sll_new_smoother_spline_1d

  

end module sll_m_kernel_smoother_spline_1d
