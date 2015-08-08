!> @ingroup particle_methods
!> @author Katharina Kormann
!> @brief Kernel smoother for 2d with splines of arbitrary degree placed on a uniform mesh.
!> @details  Spline with index i starts at point i
module sll_m_kernel_smoother_spline_2d

#include "sll_working_precision.h"
#include "sll_memory.h"

  use sll_m_kernel_smoother_base
  use sll_module_pic_base
  use sll_arbitrary_degree_splines
  
  implicit none
  
  !>  Spline kernel smoother in 2d.
  type, extends(sll_kernel_smoother_base) :: sll_kernel_smoother_spline_2d

     ! Information about the 2d mesh
     sll_real64 :: delta_x(2)  !< Value of grid spacing along both directions.
     sll_real64 :: domain(2,2) !< Definition of the domain: domain(1,1) = x1_min, domain(1,2) = x2_min,  domain(2,1) = x1_max, domain(2,2) = x2_max
     
     ! Information about the particles
     sll_int32 :: no_particles !< Number of particles of underlying PIC method (processor local)

     ! 
     sll_int32 :: spline_degree !< Degree of smoothing kernel spline
     sll_int32 :: n_span !< Number of intervals where spline non zero (spline_degree + 1)

     ! Sparse structure for shape factors
     sll_int32,  allocatable :: index_grid(:,:)    !< First dof index with contribution from th
     sll_real64, allocatable :: values_grid(:,:,:) !< Values of the space factors in each dimesion.
     
   contains
     procedure :: compute_shape_factors => compute_shape_factors_spline_2d
     procedure :: accumulate_rho_from_klimontovich => accumulate_rho_from_klimontovich_spline_2d
     procedure :: accumulate_j_from_klimontovich => accumulate_j_from_klimontovich_spline_2d
     procedure :: evaluate_kernel_function => evaluate_kernel_function_spline_2d

     !procedure :: initialize => initialize_smoother_spline_2d

  end type sll_kernel_smoother_spline_2d
  
contains
  !---------------------------------------------------------------------------!
  subroutine compute_shape_factors_spline_2d(this, particle_group)
    class( sll_kernel_smoother_spline_2d), intent(inout) :: this !< kernel smoother object
    class( sll_particle_group_base), intent(in)     :: particle_group  !< particle group

    ! local variables
    sll_real64 :: xi(3)
    sll_int32 :: i_part
    sll_real64 :: spline_val(this%n_span)

    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x(i_part)
       xi(1:2) = (xi(1:2) - this%domain(:,1)) /&
            this%delta_x
       this%index_grid(:,i_part) = ceiling(xi(1:2))
       xi(1:2) = xi(1:2) - real(this%index_grid(:,i_part) -1,f64)
       this%index_grid(:,i_part) =  this%index_grid(:,i_part) - this%spline_degree
       spline_val = uniform_b_splines_at_x(this%spline_degree, xi(1))!basis_functions(xi) ! TODO
       this%values_grid(:,1,i_part) = spline_val
       spline_val = uniform_b_splines_at_x(this%spline_degree, xi(2))
       this%values_grid(:,2,i_part) = spline_val
    end do

    !write(21,*) this%values_grid
    !write(21,*) this%index_grid

  end subroutine compute_shape_factors_spline_2d

  !---------------------------------------------------------------------------!
  subroutine accumulate_rho_from_klimontovich_spline_2d(this, particle_group,&
       rho_dofs)
    class( sll_kernel_smoother_spline_2d), intent(in)    :: this
    class( sll_particle_group_base), intent(in)     :: particle_group
    sll_real64, intent(inout)                       :: rho_dofs(:)
    
    !local variables
    sll_int32 :: i_part, i1, i2, index2d
    sll_int32 :: index1d(2)

    do i_part = 1, particle_group%n_particles
       index1d = this%index_grid(i_part,:)
       do i1 = 1, this%n_span
          index1d(1) = this%index_grid(i_part,1)+i1-2
          do i2 = 1, this%n_span
             index1d(2) = this%index_grid(i_part,2)+i2-2
             index2d = index_1dto2d_column_major(this,index1d)
             rho_dofs(index2d) = rho_dofs(index2d) +&
                  (particle_group%get_charge(i_part) * &
                  this%values_grid(i1, 1, i_part) *&
                  this%values_grid(i2, 2, i_part)) *&
                  real(product(this%n_grid), f64)
          end do
       end do
    end do


  end subroutine accumulate_rho_from_klimontovich_spline_2d

  !---------------------------------------------------------------------------!
  subroutine accumulate_j_from_klimontovich_spline_2d(this, particle_group,&
       j_dofs, component)
    class( sll_kernel_smoother_spline_2d), intent(in)    :: this
    class( sll_particle_group_base), intent(in)     :: particle_group
    sll_real64, intent(inout)                       :: j_dofs(:)
    sll_int32, intent(in)                           :: component
    
    !local variables
    sll_int32 :: i_part, i1, i2, index2d
    sll_int32 :: index1d(2)
    sll_real64 :: vpart(3)

    do i_part = 1, particle_group%n_particles
       index1d = this%index_grid(:,i_part)
       do i1 = 1, this%n_span
          index1d(1) = this%index_grid(1,i_part)+i1-2
          do i2 = 1, this%n_span
             index1d(2) = this%index_grid(2,i_part)+i2-2
             index2d = index_1dto2d_column_major(this,index1d)
             vpart = particle_group%get_v(i_part)
             j_dofs(index2d) = j_dofs(index2d) +&
                  (particle_group%get_charge(i_part) * &
                  vpart(component) * &
                  this%values_grid(i1,1,i_part) *&
                  this%values_grid(i2,2,i_part)) * &
                  real(product(this%n_grid), f64)
          end do
       end do
    end do
   

  end subroutine accumulate_j_from_klimontovich_spline_2d
  
  !---------------------------------------------------------------------------!
  subroutine evaluate_kernel_function_spline_2d(this, particle_group, rho_dofs, particle_values)
    class( sll_kernel_smoother_spline_2d), intent(in)    :: this
    class( sll_particle_group_base), intent(in)     :: particle_group
    sll_real64, intent(in)                       :: rho_dofs(:)
    sll_real64, intent(out)                      :: particle_values(:)
    
    !local variables
    sll_int32 :: i_part, i1, i2, index2d
    sll_int32 :: index1d(2)

    do i_part = 1, particle_group%n_particles
       particle_values(i_part) = 0.0_f64
       do i1 = 1, this%n_span
          index1d(1) = this%index_grid(1,i_part)+i1-2
          do i2 = 1, this%n_span
             index1d(2) = this%index_grid(2,i_part)+i2-2
             index2d = index_1dto2d_column_major(this,index1d)
             particle_values(i_part ) = particle_values(i_part) + &
                  rho_dofs(index2d) *  &
                  this%values_grid(i1,1, i_part) *&
                  this%values_grid(i2, 2, i_part)
          end do
       end do
    end do

  end subroutine evaluate_kernel_function_spline_2d

! Version as a subroutine. Not up to date with sll_new function.
!!$  subroutine initialize_smoother_spline_2d(this, domain, n_grid, no_particles, spline_degree)
!!$    class( sll_kernel_smoother_spline_2d), intent(inout)    :: this
!!$    sll_int32, intent(in) :: n_grid(2)
!!$    sll_real64, intent(in) :: domain(2,2)
!!$    sll_int32, intent(in) :: no_particles
!!$    sll_int32, intent(in) :: spline_degree !< Degree of smoothing kernel spline
!!$
!!$    !local variables
!!$    sll_int32 :: ierr
!!$
!!$    ! Store grid information
!!$    this%domain = domain
!!$    SLL_ALLOCATE(this%n_grid(2), ierr)
!!$    this%n_grid = n_grid
!!$    this%n_dofs = product(n_grid)
!!$    this%delta_x = (domain(:,2)-domain(:,1))/n_grid
!!$
!!$    ! Store basis function information
!!$    this%no_particles = no_particles
!!$
!!$    ! Initialize information on the spline
!!$    this%spline_degree = spline_degree
!!$    this%n_span = spline_degree + 1
!!$
!!$    ! Initialize sparse structure for shape factors
!!$    SLL_ALLOCATE(this%index_grid(no_particles, 2),ierr)
!!$    SLL_ALLOCATE(this%values_grid(no_particles, 2, this%n_span),ierr)
!!$
!!$  end subroutine initialize_smoother_spline_2d

  !-------------------------------------------------------------------------------------------
  !< Constructor 
  function sll_new_smoother_spline_2d(domain, n_grid, no_particles, spline_degree) result (this)
    class( sll_kernel_smoother_spline_2d), pointer   :: this
    sll_int32, intent(in) :: n_grid(2)
    sll_real64, intent(in) :: domain(2,2)
    sll_int32, intent(in) :: no_particles
    sll_int32, intent(in) :: spline_degree !< Degree of smoothing kernel spline

    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE( this, ierr)

    ! Store grid information
    this%domain = domain
    SLL_ALLOCATE(this%n_grid(2), ierr)
    this%n_grid = n_grid
    this%n_dofs = product(n_grid)
    this%delta_x = (domain(:,2)-domain(:,1))/real(n_grid, f64)

    ! Store basis function information
    this%no_particles = no_particles

    ! Initialize information on the spline
    this%spline_degree = spline_degree
    this%n_span = spline_degree + 1

    ! Initialize sparse structure for shape factors
    SLL_ALLOCATE(this%index_grid(2, no_particles),ierr)
    SLL_ALLOCATE(this%values_grid(this%n_span, 2, no_particles),ierr)

  end function sll_new_smoother_spline_2d

  
  !< This function computes the index of a 1D array that stores 2D data in column major ordering. It also takes periodic boundary conditions into account.
  function index_1dto2d_column_major(this, index1d) result(index2d)
    class( sll_kernel_smoother_spline_2d), intent(in)    :: this !< Kernel smoother object.
    sll_int32 :: index1d(2) !< 2d array with indices along each of the two directions (start counting with zero).
    sll_int32 :: index2d    !< Corresponding index in 1d array representing 2d data (start counting with one).

    index1d(1) = modulo(index1d(1), this%n_grid(1))
    index1d(2) = modulo(index1d(2), this%n_grid(2))
    index2d = index1d(1) + index1d(2)*this%n_grid(1) + 1

  end function index_1dto2d_column_major
  

end module sll_m_kernel_smoother_spline_2d
