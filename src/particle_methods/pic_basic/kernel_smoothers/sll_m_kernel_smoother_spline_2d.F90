!> @ingroup kernel_smoothers
!> @author Katharina Kormann, IPP
!> @brief Kernel smoother for 2d with splines of arbitrary degree placed on a uniform mesh.
!> @details  Spline with index i starts at point i
module sll_m_kernel_smoother_spline_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_arbitrary_degree_splines, only: &
    sll_s_uniform_b_splines_at_x

  use sll_m_kernel_smoother_base, only: &
    sll_p_collocation, &
    sll_c_kernel_smoother, &
    sll_p_galerkin

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  implicit none

  public :: &
       sll_s_new_kernel_smoother_spline_2d_ptr, & 
       sll_s_new_kernel_smoother_spline_2d, &
       sll_t_kernel_smoother_spline_2d
  
  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  !>  Spline kernel smoother in 2d.
  type, public, extends(sll_c_kernel_smoother) :: sll_t_kernel_smoother_spline_2d

     ! Information about the 2d mesh
     sll_real64 :: delta_x(2)  !< Value of grid spacing along both directions.
     sll_real64 :: domain(2,2) !< Definition of the domain: domain(1,1) = x1_min, domain(2,1) = x2_min,  domain(1,2) = x1_max, domain(2,2) = x2_max
     
     ! Information about the particles
     sll_int32 :: no_particles !< Number of particles of underlying PIC method (processor local)

     ! 
     sll_int32 :: spline_degree !< Degree of smoothing kernel spline
     sll_int32 :: n_span !< Number of intervals where spline non zero (spline_degree + 1)
!     sll_int32 :: smoothing_type !< Defines whether we use Galerkin or collocation in order to decide on scaling when accumulating
     sll_real64 :: scaling

     ! Internal work space
     sll_real64, allocatable :: spline_val(:,:)
     
   contains
     procedure :: compute_shape_factor_spline_2d !> Compute the shape factor
     procedure :: add_charge => add_charge_single_spline_2d !> Accumulate the density
     !procedure :: accumulate_j_from_klimontovich => accumulate_j_from_klimontovich_spline_2d !> Accumulate a component of the current density
     procedure :: evaluate => evaluate_field_single_spline_2d !> Evaluate the spline with given coefficients
     procedure :: evaluate_multiple => evaluate_multiple_spline_2d
     procedure :: add_current_update_v => add_current_update_v_spline_2d
     procedure :: init => init_spline_2d
     procedure :: free => free_spline_2d


  end type sll_t_kernel_smoother_spline_2d
  
contains
  !---------------------------------------------------------------------------!
  !> Helper function computing shape factor
  subroutine compute_shape_factor_spline_2d(self, position, indices)
    class( sll_t_kernel_smoother_spline_2d), intent(inout) :: self !< kernel smoother object
    sll_real64, intent( in ) :: position(2)
    sll_int32, intent( out ) :: indices(2)

    ! local variables
    sll_real64 :: xi(2)


    xi(1:2) = (position(1:2) - self%domain(:,1)) /&
         self%delta_x
    indices = ceiling(xi(1:2))
    xi(1:2) = xi(1:2) - real(indices -1,f64)
    indices =  indices - self%spline_degree
    call sll_s_uniform_b_splines_at_x(self%spline_degree, xi(1), self%spline_val(1:self%n_span,1))
    call sll_s_uniform_b_splines_at_x(self%spline_degree, xi(2), self%spline_val(1:self%n_span,2))
    !self%spline_val(1:self%n_span,1) = sll_f_uniform_b_splines_at_x(self%spline_degree, xi(1))
    !self%spline_val(1:self%n_span,2) = sll_f_uniform_b_splines_at_x(self%spline_degree, xi(2))

  end subroutine compute_shape_factor_spline_2d

  !---------------------------------------------------------------------------!
  !> Add charge of single particle
  subroutine add_charge_single_spline_2d(self, position, marker_charge, rho_dofs)
    class( sll_t_kernel_smoother_spline_2d), intent(inout)    :: self !< kernel smoother object
    sll_real64, intent( in ) :: position(self%dim)
    sll_real64, intent( in ) :: marker_charge
    sll_real64, intent(inout)                       :: rho_dofs(self%n_dofs ) !< spline coefficient of accumulated density
    
    !local variables
    sll_int32 :: i1, i2, index2d
    sll_int32 :: index1d(2)
    sll_int32  :: indices(2)
    

    call compute_shape_factor_spline_2d(self, position, indices)
    do i1 = 1, self%n_span
       index1d(1) = indices(1)+i1-2
       do i2 = 1, self%n_span
          index1d(2) = indices(2)+i2-2
          index2d = index_1dto2d_column_major(self,index1d)
          rho_dofs(index2d) = rho_dofs(index2d) +&
               ( marker_charge* self%scaling * &
               self%spline_val(i1,1) * self%spline_val(i2,2))
       end do
    end do

  end subroutine add_charge_single_spline_2d


 !---------------------------------------------------------------------------!
  !> Add current and update v for single particle
  subroutine add_current_update_v_spline_2d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_kernel_smoother_spline_2d), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in) :: position_old(self%dim)
    sll_real64, intent(in) :: position_new(self%dim)
    sll_real64, intent(in) :: marker_charge
    sll_real64, intent(in) :: qoverm
    sll_real64, intent(in) :: bfield_dofs(self%n_dofs)
    sll_real64, intent(inout) :: vi(:)
    sll_real64, intent(inout) :: j_dofs(self%n_dofs)

    print*, 'add_current_update_v_spline_2d not implemented.'

  end subroutine add_current_update_v_spline_2d

!!$  !---------------------------------------------------------------------------!
!!$  subroutine accumulate_j_from_klimontovich_spline_2d(self, particle_group,&
!!$       j_dofs, component, i_weight)
!!$    class( sll_t_kernel_smoother_spline_2d), intent(in)    :: self !< kernel smoother object
!!$    class( sll_c_particle_group_base), intent(in)     :: particle_group !< particle group
!!$    sll_real64, intent(inout)                       :: j_dofs(:) !< spline coefficients ofcomponent \a component accumulated current density
!!$    sll_int32, intent(in)                           :: component !< component of \a j_dofs to be accumulated.
!!$    sll_int32, optional             , intent( in ) :: i_weight !< No. of the weight that should be used in the accumulation.
!!$
!!$    !local variables
!!$    sll_int32 :: i_part, i1, i2, index2d
!!$    sll_int32 :: index1d(2)
!!$    sll_real64 :: vpart(3)
!!$    sll_real64 :: wi(1)
!!$    sll_int32  :: i_wi
!!$    
!!$    i_wi = 1
!!$    if(present(i_weight)) i_wi = i_weight
!!$
!!$    do i_part = 1, particle_group%n_particles
!!$       wi = particle_group%get_charge(i_part, i_wi)
!!$       do i1 = 1, self%n_span
!!$          index1d(1) = self%index_grid(1,i_part)+i1-2
!!$          do i2 = 1, self%n_span
!!$             index1d(2) = self%index_grid(2,i_part)+i2-2
!!$             index2d = index_1dto2d_column_major(self,index1d)
!!$             vpart = particle_group%get_v(i_part)
!!$             j_dofs(index2d) = j_dofs(index2d) +&
!!$                  (wi(1) * &
!!$                  vpart(component) * &
!!$                  self%values_grid(i1,1,i_part) *&
!!$                  self%values_grid(i2,2,i_part)) 
!!$          end do
!!$       end do
!!$    end do
!!$   
!!$    if (self%smoothing_type == sll_p_collocation) then
!!$       j_dofs = j_dofs /product(self%delta_x)
!!$    end if
!!$
!!$  end subroutine accumulate_j_from_klimontovich_spline_2d
  

  !---------------------------------------------------------------------------!
  !> Evaluate field with given dofs at position \a position
  subroutine evaluate_field_single_spline_2d(self, position, field_dofs, field_value)
    class( sll_t_kernel_smoother_spline_2d), intent(inout)    :: self !< kernel smoother object    
    sll_real64, intent( in ) :: position(self%dim)
    sll_real64, intent(in)                               :: field_dofs(self%n_dofs) !< Degrees of freedom in kernel representation.
    sll_real64, intent(out) :: field_value
    
    !local variables
    sll_int32 :: i1, i2, index2d
    sll_int32 :: index1d(2)
    sll_int32  :: indices(2)
    

    call compute_shape_factor_spline_2d(self, position, indices)


    field_value = 0.0_f64
    do i1 = 1, self%n_span
       index1d(1) = indices(1)+i1-2
       do i2 = 1, self%n_span
          index1d(2) = indices(2)+i2-2
          index2d = index_1dto2d_column_major(self,index1d)
          field_value = field_value + &
               field_dofs(index2d) *  &
               self%spline_val(i1,1) *&
               self%spline_val(i2,2)
       end do
    end do

  end subroutine evaluate_field_single_spline_2d

  !---------------------------------------------------------------------------!
  !> Evaluate multiple fields at position \a position
  subroutine evaluate_multiple_spline_2d(self, position, components, field_dofs, field_value)
    class( sll_t_kernel_smoother_spline_2d), intent(inout)    :: self !< kernel smoother object    
    sll_real64, intent( in ) :: position(self%dim)
    sll_int32, intent(in) :: components(:)
    sll_real64, intent(in)                               :: field_dofs(:,:) !< Degrees of freedom in kernel representation.
    sll_real64, intent(out) :: field_value(:)
    
    !local variables
    sll_int32 :: i1, i2, index2d
    sll_int32 :: index1d(2)
    sll_int32  :: indices(2)
    

    call compute_shape_factor_spline_2d(self, position, indices)


    field_value = 0.0_f64
    do i1 = 1, self%n_span
       index1d(1) = indices(1)+i1-2
       do i2 = 1, self%n_span
          index1d(2) = indices(2)+i2-2
          index2d = index_1dto2d_column_major(self,index1d)
          field_value = field_value + &
               field_dofs(index2d,components) *  &
               self%spline_val(i1,1) *&
               self%spline_val(i2,2)
       end do
    end do

  end subroutine evaluate_multiple_spline_2d

  !-------------------------------------------------------------------------------------------
  !> Destructor
    subroutine free_spline_2d(self)
    class (sll_t_kernel_smoother_spline_2d), intent( inout ) :: self !< Kernel smoother object 

    deallocate(self%spline_val)
    

  end subroutine free_spline_2d


  !-------------------------------------------------------------------------------------------
  !< Constructor 
  subroutine init_spline_2d(self, domain, n_grid, no_particles, spline_degree, smoothing_type)
    class( sll_t_kernel_smoother_spline_2d), intent(out)   :: self
    sll_int32, intent(in) :: n_grid(2) !< no. of spline coefficients
    sll_real64, intent(in) :: domain(2,2) !< lower and upper bounds of the domain
    sll_int32, intent(in) :: no_particles !< no. of particles
    sll_int32, intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32, intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr

    self%dim = 2

    ! Store grid information
    self%domain = domain
    SLL_ALLOCATE(self%n_grid(2), ierr)
    self%n_grid = n_grid
    self%n_dofs = product(n_grid)
    self%delta_x = (domain(:,2)-domain(:,1))/real(n_grid, f64)

    ! Store basis function information
    self%no_particles = no_particles

    ! Initialize information on the spline
    self%spline_degree = spline_degree
    self%n_span = spline_degree + 1

    ! Initialize information on smoothing type
    if (smoothing_type == sll_p_collocation) then
       self%scaling = 1.0_f64/product(self%delta_x)
    elseif (smoothing_type == sll_p_galerkin) then
       self%scaling = 1.0_f64
    else
       print*, 'Smoothing Type ', smoothing_type, ' not implemented for kernel_smoother_spline_2d.'
    end if

    allocate( self%spline_val(self%n_span, 2), stat = ierr)
    SLL_ASSERT( ierr == 0)

  end subroutine init_spline_2d

  !> Constructor for abstract type (pointer)
  subroutine sll_s_new_kernel_smoother_spline_2d_ptr(smoother, domain, n_grid, no_particles, spline_degree, smoothing_type) 
    class( sll_c_kernel_smoother), pointer, intent(out)   :: smoother
    sll_int32, intent(in) :: n_grid(2) !< no. of spline coefficients
    sll_real64, intent(in) :: domain(2,2) !< lower and upper bounds of the domain
    sll_int32, intent(in) :: no_particles !< no. of particles
    sll_int32, intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32, intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE( sll_t_kernel_smoother_spline_2d :: smoother, ierr)
    select type (smoother)
    type is (sll_t_kernel_smoother_spline_2d)
       call smoother%init( domain, n_grid, no_particles, spline_degree, smoothing_type)
    end select

  end subroutine sll_s_new_kernel_smoother_spline_2d_ptr


 !> Constructor for abstract type (allocatable)
  subroutine sll_s_new_kernel_smoother_spline_2d(smoother, domain, n_grid, no_particles, spline_degree, smoothing_type) 
    class( sll_c_kernel_smoother), allocatable, intent(out)   :: smoother
    sll_int32, intent(in) :: n_grid(2) !< no. of spline coefficients
    sll_real64, intent(in) :: domain(2,2) !< lower and upper bounds of the domain
    sll_int32, intent(in) :: no_particles !< no. of particles
    sll_int32, intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32, intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE( sll_t_kernel_smoother_spline_2d :: smoother, ierr)
    select type (smoother)
    type is (sll_t_kernel_smoother_spline_2d)
       call smoother%init( domain, n_grid, no_particles, spline_degree, smoothing_type)
    end select

  end subroutine sll_s_new_kernel_smoother_spline_2d


  
  !< Self function computes the index of a 1D array that stores 2D data in column major ordering. It also takes periodic boundary conditions into account.
  function index_1dto2d_column_major(self, index1d) result(index2d)
    class( sll_t_kernel_smoother_spline_2d), intent(in)    :: self !< Kernel smoother object.
    sll_int32 :: index1d(2) !< 2d array with indices along each of the two directions (start counting with zero).
    sll_int32 :: index2d    !< Corresponding index in 1d array representing 2d data (start counting with one).

    index1d(1) = modulo(index1d(1), self%n_grid(1))
    index1d(2) = modulo(index1d(2), self%n_grid(2))
    index2d = index1d(1) + index1d(2)*self%n_grid(1) + 1

  end function index_1dto2d_column_major

end module sll_m_kernel_smoother_spline_2d
