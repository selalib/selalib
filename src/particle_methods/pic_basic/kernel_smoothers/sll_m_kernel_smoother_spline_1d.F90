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

  use sll_m_gauss_legendre_integration, only : &
       gauss_legendre_points_and_weights

  use sll_m_kernel_smoother_base, only: &
    sll_collocation, &
    sll_c_kernel_smoother, &
    sll_galerkin

  use sll_m_particle_group_base, only: &
    sll_particle_group_base


  implicit none

  public :: &
    sll_t_kernel_smoother_spline_1d, &
    sll_new_smoother_spline_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !>  Spline kernel smoother in1d.
  type, extends(sll_c_kernel_smoother) :: sll_t_kernel_smoother_spline_1d

     ! Information about the 1d mesh
     sll_real64 :: delta_x(1)  !< Value of grid spacing along both directions.
     sll_real64 :: domain(1,2) !< Definition of the domain: domain(1,1) = x1_min  domain(1,2) = x1_max
     
     ! Information about the particles
     sll_int32 :: no_particles !< Number of particles of underlying PIC method (processor local)

     ! 
     sll_int32 :: spline_degree !< Degree of smoothing kernel spline
     sll_int32 :: n_span !< Number of intervals where spline non zero (spline_degree + 1)
!     sll_int32 :: smoothing_type !< Defines whether we use Galerkin or collocation in order to decide on scaling when accumulating
     sll_real64 :: scaling
     sll_int32  :: n_quad_points

     sll_real64, allocatable :: spline_val(:)
     sll_real64, allocatable :: quad_xw(:,:)

     
   contains
     procedure :: add_charge => add_charge_single_spline_1d !> 
     procedure :: add_current_update_v => add_current_update_v_spline_1d
     procedure :: evaluate => evaluate_field_single_spline_1d !> Evaluate spline function with given coefficients
     procedure :: evaluate_multiple => evaluate_multiple_spline_1d !> Evaluate multiple spline functions with given coefficients
     procedure :: update_jv !> elper function to compute the integral of j using Gauss quadrature

  end type sll_t_kernel_smoother_spline_1d
  
contains

  !---------------------------------------------------------------------------!
  subroutine add_charge_single_spline_1d(this, position, weight, rho_dofs)
    class( sll_t_kernel_smoother_spline_1d ), intent(inout) :: this !< kernel smoother object
     sll_real64, intent(in)  :: position(this%dim) !< Position of the particle
    sll_real64,                intent( in ) :: weight !< Weight of the particle
    sll_real64,                 intent( inout ) :: rho_dofs(this%n_dofs) !< Coefficient vector of the charge distribution

    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d, index
    sll_real64 :: xi(1)
   
    xi(1) = (position(1) - this%domain(1,1))/this%delta_x(1)
    index = ceiling(xi(1))
    xi(1) = xi(1) - real(index-1, f64)
    index = index - this%spline_degree
    this%spline_val = uniform_b_splines_at_x(this%spline_degree, xi(1))

    do i1 = 1, this%n_span
       index1d = modulo(index+i1-2,this%n_grid(1))+1
       rho_dofs(index1d) = rho_dofs(index1d) +&
            (weight * this%spline_val(i1)* this%scaling)
    end do

  end subroutine add_charge_single_spline_1d


  subroutine add_current_update_v_spline_1d (this, position_old, position_new, weight, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_kernel_smoother_spline_1d), intent(inout) :: this !< kernel smoother object
    sll_real64, intent(in) :: position_old(this%dim)
    sll_real64, intent(in) :: position_new(this%dim)
    sll_real64, intent(in) :: weight
    sll_real64, intent(in) :: qoverm
    sll_real64, intent(in) :: bfield_dofs(this%n_dofs)
    sll_real64, intent(inout) :: vi(:)
    sll_real64, intent(inout) :: j_dofs(this%n_dofs)

    ! local variables
    sll_real64 :: xi
    sll_int32  :: index_old, index_new, ind
    sll_real64 :: r_old, r_new

    ! Read out particle position and velocity
    ! Compute index_old, the index of the last DoF on the grid the particle contributes to, and r_old, its position (normalized to cell size one).
       xi = (position_old(1) - this%domain(1,1)) /&
            this%delta_x(1)
       index_old = floor(xi)
       r_old = xi - real(index_old,f64)

       ! Compute the new box index index_new and normalized position r_old.
       xi = (position_new(1) - this%domain(1,1)) /&
            this%delta_x(1)
       index_new = floor(xi)
       r_new = xi - real(index_new ,f64) 

       if (index_old == index_new) then
          if (r_old < r_new) then
             call this%update_jv(r_old, r_new, index_old, weight, &
                  qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          else
             call this%update_jv(r_new, r_old, index_old, weight, qoverm, &
                  -1.0_f64, vi(2), j_dofs, bfield_dofs)
          end if
       elseif (index_old < index_new) then
          call this%update_jv (r_old, 1.0_f64, index_old, weight, &
               qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          call this%update_jv (0.0_f64, r_new, index_new, weight, &
               qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          do ind = index_old+1, index_new-1
             call this%update_jv (0.0_f64, 1.0_f64, ind, weight, &
                  qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          end do
       else
          call this%update_jv (r_new, 1.0_f64, index_new, weight, qoverm, &
               -1.0_f64, vi(2), j_dofs, bfield_dofs)
          call this%update_jv (0.0_f64, r_old, index_old, weight, qoverm, &
               -1.0_f64, vi(2), j_dofs, bfield_dofs)
          do ind = index_new+1, index_old-1
             call this%update_jv (0.0_f64, 1.0_f64, ind, weight, qoverm, &
                  -1.0_f64, vi(2), j_dofs, bfield_dofs)
          end do
       end if    


     end subroutine add_current_update_v_spline_1d

 ! TODO: This is hard coded for quadratic, cubic splines. Make general.
 subroutine update_jv(this, lower, upper, index, weight, qoverm, sign, vi, j_dofs, bfield_dofs)
   class(sll_t_kernel_smoother_spline_1d), intent(inout) :: this !< time splitting object 
   sll_real64, intent(in) :: lower
   sll_real64, intent(in) :: upper
   sll_int32,  intent(in) :: index
   sll_real64, intent(in) :: weight
   sll_real64, intent(in) :: qoverm
   sll_real64, intent(in) :: sign
   sll_real64, intent(inout) :: vi
   sll_real64, intent(in) :: bfield_dofs(this%n_dofs)
   sll_real64, intent(inout) :: j_dofs(this%n_dofs)

   !Local variables
   sll_int32  :: ind, i_grid, i_mod, n_cells, j


   n_cells = this%n_grid(1)

   this%quad_xw = gauss_legendre_points_and_weights(this%n_quad_points, lower, upper)

   this%spline_val = this%quad_xw(2,1) * &
        uniform_b_splines_at_x(this%spline_degree, this%quad_xw(1,1))
   do j=2,this%n_quad_points
      this%spline_val = this%spline_val + &
           this%quad_xw(2,j) * &
           uniform_b_splines_at_x(this%spline_degree, this%quad_xw(1,j))
   end do
   this%spline_val = this%spline_val * sign*this%delta_x(1)

   ind = 1
   do i_grid = index - this%spline_degree , index
      i_mod = modulo(i_grid, n_cells ) + 1
      j_dofs(i_mod) = j_dofs(i_mod) + &
           (weight*this%spline_val(ind)* this%scaling)
      vi = vi - qoverm* this%spline_val(ind)*bfield_dofs(i_mod)
      ind = ind + 1
   end do

 end subroutine update_jv

!!$  !---------------------------------------------------------------------------!
!!$  subroutine add_shape_factor_single_spline_1d(this, position, i_part)
!!$    class( sll_kernel_smoother_spline_1d ), intent(in) :: this !< kernel smoother object  
!!$    sll_real64, intent(in)  :: position(this%dim) !< Position of the particle     
!!$    sll_int32, optional,       intent( in ) :: i_part !< no. of the particle
!!$    
!!$    ! local variables
!!$    sll_real64 :: spline_val(this%n_span)
!!$
!!$    position(1) = (position(1) - this%domain(1,1)) /&
!!$         this%delta_x(1)
!!$    ! This is the index of the last spline the particle contributes to
!!$    this%index_grid(:,i_part) = ceiling(position(1))
!!$    position(1) = position(1) - real(this%index_grid(1,i_part) -1,f64)
!!$    ! Now we subtract the degree of the spline to get the index of the first spline.
!!$    this%index_grid(:,i_part) =  this%index_grid(:,i_part) - this%spline_degree
!!$    spline_val = uniform_b_splines_at_x(this%spline_degree, position(1))
!!$    this%values_grid(:,1,i_part) = spline_val
!!$    
!!$    
!!$  end subroutine add_shape_factor_single_spline_1d



  !---------------------------------------------------------------------------!
  subroutine evaluate_field_single_spline_1d(this, position, field_dofs, field_value)
    class (sll_t_kernel_smoother_spline_1d), intent( inout ) :: this !< Kernel smoother object 
    sll_real64,                intent( in ) :: position(this%dim) !< Position of the particle
    sll_real64,                    intent( in ) :: field_dofs(this%n_dofs) !< Coefficient vector for the field DoFs
    sll_real64, intent( out)                              :: field_value !< Value(s) of the electric fields at given position
    
    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d, index
    sll_real64 :: xi(1)

    xi(1) = (position(1) - this%domain(1,1))/this%delta_x(1)
    index = ceiling(xi(1))
    xi(1) = xi(1) - real(index-1, f64)
    index = index - this%spline_degree
    this%spline_val = uniform_b_splines_at_x(this%spline_degree, xi(1))

    field_value = 0.0_f64
    do i1 = 1, this%n_span
       index1d = modulo(index+i1-2, this%n_grid(1))+1
       field_value = field_value + &
            field_dofs(index1d) *  &
            this%spline_val(i1)
    end do

  end subroutine evaluate_field_single_spline_1d



  !---------------------------------------------------------------------------!
  subroutine evaluate_multiple_spline_1d(this, position, components, field_dofs, field_value)
    class (sll_t_kernel_smoother_spline_1d), intent( inout ) :: this !< Kernel smoother object 
    sll_real64,                intent( in ) :: position(this%dim) !< Position of the particle
    sll_int32, intent(in) :: components(:)
    sll_real64,                    intent( in ) :: field_dofs(:,:) !< Coefficient vector for the field DoFs
    sll_real64, intent(out)                              :: field_value(:) !< Value(s) of the electric fields at given position
    
    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d, index
    sll_real64 :: xi(1)

    ! TODO: Add assertions on sive of field_dofs(this%n_dofs, size(field_value)

    xi(1) = (position(1) - this%domain(1,1))/this%delta_x(1)
    index = ceiling(xi(1))
    xi(1) = xi(1) - real(index-1, f64)
    index = index - this%spline_degree
    this%spline_val = uniform_b_splines_at_x(this%spline_degree, xi(1))

    field_value = 0.0_f64
    do i1 = 1, this%n_span
       index1d = modulo(index+i1-2, this%n_grid(1))+1
       field_value = field_value + &
            field_dofs(index1d,components) *  &
            this%spline_val(i1)
    end do

  end subroutine evaluate_multiple_spline_1d



  !-------------------------------------------------------------------------------------------
  !< Constructor 
  function sll_new_smoother_spline_1d(domain, n_grid, no_particles, spline_degree, smoothing_type) result (this)
    class( sll_t_kernel_smoother_spline_1d), pointer   :: this !< kernel smoother object
    sll_int32, intent(in) :: n_grid(1) !< number of DoFs (spline coefficients)
    sll_real64, intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32, intent(in) :: no_particles !< number of particles
    sll_int32, intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32, intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr


    SLL_ALLOCATE( this, ierr)

    this%dim = 1

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
    if (smoothing_type == SLL_COLLOCATION) then
       this%scaling = 1.0_f64/this%delta_x(1)
    elseif (smoothing_type == SLL_GALERKIN) then
       this%scaling = 1.0_f64
    else
       print*, 'Smoothing Type ', smoothing_type, ' not implemented for kernel_smoother_spline_1d.'
    end if
    
    this%n_quad_points = (this%spline_degree+2)/2

    ALLOCATE( this%spline_val(this%n_span))
    ALLOCATE( this%quad_xw(2,this%n_quad_points))

  end function sll_new_smoother_spline_1d

  

end module sll_m_kernel_smoother_spline_1d
