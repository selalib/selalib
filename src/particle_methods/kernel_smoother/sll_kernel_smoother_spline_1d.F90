!> @ingroup particle_methods
!> @author Katharina Kormann
!> @brief Kernel smoother for 2d with splines of arbitrary degree placed on a uniform mesh.
!> @details Spline with index i starts at point i
module sll_m_kernel_smoother_spline_1d

#include "sll_working_precision.h"
#include "sll_memory.h"

  use sll_m_kernel_smoother_base
  use sll_m_pic_base
  use sll_arbitrary_degree_splines
  
  implicit none
  private

  public :: sll_new_smoother_spline_1d
  
  !>  Spline kernel smoother in1d.
  type, public, extends(sll_kernel_smoother_base) :: sll_kernel_smoother_spline_1d

     ! Information about the 1d mesh
     sll_real64 :: delta_x(1)  !< Value of grid spacing along both directions.
     sll_real64 :: domain(1,2) !< Definition of the domain: domain(1,1) = x1_min  domain(1,2) = x1_max
     
     ! Information about the particles
     sll_int32 :: no_particles !< Number of particles of underlying PIC method (processor local)

     ! 
     sll_int32 :: spline_degree !< Degree of smoothing kernel spline
     sll_int32 :: n_span !< Number of intervals where spline non zero (spline_degree + 1)

     ! Sparse structure for shape factors
     sll_int32,  allocatable :: index_grid(:,:)    !< First dof index with contribution from th
     sll_real64, allocatable :: values_grid(:,:,:) !< Values of the space factors in each dimesion.
     
   contains
     procedure :: compute_shape_factors => compute_shape_factors_spline_1d
     procedure :: accumulate_rho_from_klimontovich => accumulate_rho_from_klimontovich_spline_1d
     procedure :: accumulate_j_from_klimontovich => accumulate_j_from_klimontovich_spline_1d
     procedure :: evaluate_kernel_function => evaluate_kernel_function_spline_1d


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
    class( sll_kernel_smoother_spline_1d), intent(in)    :: this
    class( sll_particle_group_base), intent(in)     :: particle_group
    sll_real64, intent(inout)                       :: rho_dofs(:)
    
    !local variables
    sll_int32 :: i_part, i1
    sll_int32 :: index1d

    do i_part = 1, particle_group%n_particles
       do i1 = 1, this%n_span
          index1d = modulo(this%index_grid(1,i_part)+i1-2,this%n_grid(1))+1
          rho_dofs(index1d) = rho_dofs(index1d) +&
               (particle_group%get_charge(i_part) * &
               this%values_grid(i1, 1, i_part)) *&
               real(this%n_grid(1), f64)
       end do
    end do


  end subroutine accumulate_rho_from_klimontovich_spline_1d

  !---------------------------------------------------------------------------!
  subroutine accumulate_j_from_klimontovich_spline_1d(this, particle_group,&
       j_dofs, component)
    class( sll_kernel_smoother_spline_1d), intent(in)    :: this
    class( sll_particle_group_base), intent(in)     :: particle_group
    sll_real64, intent(inout)                       :: j_dofs(:)
    sll_int32, intent(in)                           :: component
    
    !local variables
    sll_int32 :: i_part, i1
    sll_int32 :: index1d
    sll_real64 :: vpart(3)

    do i_part = 1, particle_group%n_particles
       do i1 = 1, this%n_span
          index1d = modulo(this%index_grid(1,i_part)+i1-2,this%n_grid(1))+1
          vpart = particle_group%get_v(i_part)
          j_dofs(index1d) = j_dofs(index1d) +&
               (particle_group%get_charge(i_part) * &
               vpart(component) * &
               this%values_grid(i1,1,i_part)) *&
               real(product(this%n_grid), f64)
       end do
    end do
   

  end subroutine accumulate_j_from_klimontovich_spline_1d
  
  !---------------------------------------------------------------------------!
  subroutine evaluate_kernel_function_spline_1d(this, particle_group, rho_dofs, particle_values)
    class( sll_kernel_smoother_spline_1d), intent(in)    :: this
    class( sll_particle_group_base), intent(in)     :: particle_group
    sll_real64, intent(in)                       :: rho_dofs(:)
    sll_real64, intent(out)                      :: particle_values(:)
    
    !local variables
    sll_int32 :: i_part, i1
    sll_int32 :: index1d

    do i_part = 1, particle_group%n_particles
       particle_values(i_part) = 0.0_f64
       do i1 = 1, this%n_span
          index1d = modulo(this%index_grid(1,i_part)+i1-2, this%n_grid(1))+1
          particle_values(i_part ) = particle_values(i_part) + &
               rho_dofs(index1d) *  &
               this%values_grid(i1,1, i_part)
       end do
    end do

  end subroutine evaluate_kernel_function_spline_1d

  !-------------------------------------------------------------------------------------------
  !< Constructor 
  function sll_new_smoother_spline_1d(domain, n_grid, no_particles, spline_degree) result (this)
    class( sll_kernel_smoother_spline_1d), pointer   :: this
    sll_int32, intent(in) :: n_grid(1)
    sll_real64, intent(in) :: domain(2)
    sll_int32, intent(in) :: no_particles
    sll_int32, intent(in) :: spline_degree !< Degree of smoothing kernel spline

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

    ! Initialize sparse structure for shape factors
    SLL_ALLOCATE(this%index_grid(1, no_particles),ierr)
    SLL_ALLOCATE(this%values_grid(this%n_span, 1, no_particles),ierr)

  end function sll_new_smoother_spline_1d

  

end module sll_m_kernel_smoother_spline_1d
