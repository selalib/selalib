!> @ingroup particle_mesh_coupling
!> @author Katharina Kormann, IPP
!> @brief Kernel smoother for 2d with splines of arbitrary degree placed on a uniform mesh.
!> @details Spline with index i starts at point i
module sll_m_particle_mesh_coupling_spline_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_arbitrary_degree_splines, only: &
    sll_s_uniform_b_splines_at_x

  use sll_m_gauss_legendre_integration, only : &
    sll_f_gauss_legendre_points_and_weights

  use sll_m_particle_mesh_coupling_base, only: &
    sll_p_collocation, &
    sll_c_particle_mesh_coupling, &
    sll_p_galerkin

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base


  implicit none

  public :: &
    sll_t_particle_mesh_coupling_spline_1d, &
    sll_s_new_particle_mesh_coupling_spline_1d_ptr, &
    sll_s_new_particle_mesh_coupling_spline_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !>  Spline kernel smoother in1d.
  type, extends(sll_c_particle_mesh_coupling) :: sll_t_particle_mesh_coupling_spline_1d

     ! Information about the 1d mesh
     sll_real64 :: delta_x(1)  !< Value of grid spacing along both directions.
     sll_real64 :: domain(1,2) !< Definition of the domain: domain(1,1) = x1_min  domain(1,2) = x1_max
     
     ! Information about the particles
     sll_int32  :: no_particles !< Number of particles of underlying PIC method (processor local)

     ! 
     sll_int32  :: spline_degree !< Degree of smoothing kernel spline
     sll_int32  :: n_span !< Number of intervals where spline non zero (spline_degree + 1)
     sll_real64 :: scaling !< Scaling factor depending on whether Galerkin or collocation
     sll_int32  :: n_quad_points !< Number of quadrature points

     sll_real64, allocatable :: spline_val(:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_val_more(:) !< more scratch data for spline evaluation
     sll_real64, allocatable :: quad_xw(:,:) !< quadrature weights and points

     
   contains
     procedure :: add_charge => add_charge_single_spline_1d !> Add charge of one particle
     procedure :: add_current_update_v => add_current_update_v_spline_1d !> Add current of one particle
     procedure :: evaluate => evaluate_field_single_spline_1d !> Evaluate spline function with given coefficients
     procedure :: evaluate_multiple => evaluate_multiple_spline_1d !> Evaluate multiple spline functions with given coefficients
     procedure :: update_jv !> elper function to compute the integral of j using Gauss quadrature
     procedure :: init => init_spline_1d !> Constructor
     procedure :: free => free_spline_1d !> Destructor

  end type sll_t_particle_mesh_coupling_spline_1d
  
contains

  !---------------------------------------------------------------------------!
  !> Add charge of one particle
  subroutine add_charge_single_spline_1d(self, position, marker_charge, rho_dofs)
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: rho_dofs(self%n_dofs) !< Coefficient vector of the charge distribution

    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d, index
    sll_real64 :: xi(1)
   
    xi(1) = (position(1) - self%domain(1,1))/self%delta_x(1)
    index = ceiling(xi(1))
    xi(1) = xi(1) - real(index-1, f64)
    index = index - self%spline_degree
    !self%spline_val = sll_f_uniform_b_splines_at_x(self%spline_degree, xi(1))
    call sll_s_uniform_b_splines_at_x(self%spline_degree, xi(1), self%spline_val)

    do i1 = 1, self%n_span
       index1d = modulo(index+i1-2,self%n_grid(1))+1
       rho_dofs(index1d) = rho_dofs(index1d) +&
            (marker_charge * self%spline_val(i1)* self%scaling)
    end do

  end subroutine add_charge_single_spline_1d


  !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting)
  subroutine add_current_update_v_spline_1d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_1d), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
    sll_real64, intent(in)    :: position_new(self%dim) !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(self%n_dofs) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(:) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(self%n_dofs) !< Coefficients of current expansion

    ! local variables
    sll_real64 :: xi
    sll_int32  :: index_old, index_new, ind
    sll_real64 :: r_old, r_new

    ! Read out particle position and velocity
    ! Compute index_old, the index of the last DoF on the grid the particle contributes to, and r_old, its position (normalized to cell size one).
       xi = (position_old(1) - self%domain(1,1)) /&
            self%delta_x(1)
       index_old = floor(xi)
       r_old = xi - real(index_old,f64)

       ! Compute the new box index index_new and normalized position r_old.
       xi = (position_new(1) - self%domain(1,1)) /&
            self%delta_x(1)
       index_new = floor(xi)
       r_new = xi - real(index_new ,f64) 

       if (index_old == index_new) then
          if (r_old < r_new) then
             call self%update_jv(r_old, r_new, index_old, marker_charge, &
                  qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          else
             call self%update_jv(r_new, r_old, index_old, marker_charge, qoverm, &
                  -1.0_f64, vi(2), j_dofs, bfield_dofs)
          end if
       elseif (index_old < index_new) then
          call self%update_jv (r_old, 1.0_f64, index_old, marker_charge, &
               qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          call self%update_jv (0.0_f64, r_new, index_new, marker_charge, &
               qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          do ind = index_old+1, index_new-1
             call self%update_jv (0.0_f64, 1.0_f64, ind, marker_charge, &
                  qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          end do
       else
          call self%update_jv (r_new, 1.0_f64, index_new, marker_charge, qoverm, &
               -1.0_f64, vi(2), j_dofs, bfield_dofs)
          call self%update_jv (0.0_f64, r_old, index_old, marker_charge, qoverm, &
               -1.0_f64, vi(2), j_dofs, bfield_dofs)
          do ind = index_new+1, index_old-1
             call self%update_jv (0.0_f64, 1.0_f64, ind, marker_charge, qoverm, &
                  -1.0_f64, vi(2), j_dofs, bfield_dofs)
          end do
       end if    


     end subroutine add_current_update_v_spline_1d

 !> Helper function for \a add_current_update_v.
 subroutine update_jv(self, lower, upper, index, marker_charge, qoverm, sign, vi, j_dofs, bfield_dofs)
   class(sll_t_particle_mesh_coupling_spline_1d), intent(inout) :: self !< time splitting object 
   sll_real64,                             intent(in)    :: lower
   sll_real64,                             intent(in)    :: upper
   sll_int32,                              intent(in)    :: index
   sll_real64,                             intent(in)    :: marker_charge
   sll_real64,                             intent(in)    :: qoverm
   sll_real64,                             intent(in)    :: sign
   sll_real64,                             intent(inout) :: vi
   sll_real64,                             intent(in)    :: bfield_dofs(self%n_dofs)
   sll_real64,                             intent(inout) :: j_dofs(self%n_dofs)

   !Local variables
   sll_int32  :: ind, i_grid, i_mod, n_cells, j
   sll_real64 :: c1, c2


   n_cells = self%n_grid(1)

   c1 =  0.5_f64*(upper-lower)
   c2 =  0.5_f64*(upper+lower)
   
   call sll_s_uniform_b_splines_at_x(self%spline_degree, c1*self%quad_xw(1,1)+c2, &
        self%spline_val)
   self%spline_val = self%spline_val * (self%quad_xw(2,1)*c1)
   do j=2,self%n_quad_points
      call sll_s_uniform_b_splines_at_x(self%spline_degree, c1*self%quad_xw(1,j)+c2, &
           self%spline_val_more)
      self%spline_val = self%spline_val + self%spline_val_more * (self%quad_xw(2,j)*c1)
   end do
   self%spline_val = self%spline_val * (sign*self%delta_x(1))

   ind = 1
   do i_grid = index - self%spline_degree , index
      i_mod = modulo(i_grid, n_cells ) + 1
      j_dofs(i_mod) = j_dofs(i_mod) + &
           (marker_charge*self%spline_val(ind)* self%scaling)
      vi = vi - qoverm* self%spline_val(ind)*bfield_dofs(i_mod)
      ind = ind + 1
   end do

 end subroutine update_jv

!!$  !---------------------------------------------------------------------------!
!!$  subroutine add_shape_factor_single_spline_1d(self, position, i_part)
!!$    class( sll_kernel_smoother_spline_1d ), intent(in) :: self !< kernel smoother object  
!!$    sll_real64, intent(in)  :: position(self%dim) !< Position of the particle     
!!$    sll_int32, optional,       intent( in ) :: i_part !< no. of the particle
!!$    
!!$    ! local variables
!!$    sll_real64 :: spline_val(self%n_span)
!!$
!!$    position(1) = (position(1) - self%domain(1,1)) /&
!!$         self%delta_x(1)
!!$    ! Self is the index of the last spline the particle contributes to
!!$    self%index_grid(:,i_part) = ceiling(position(1))
!!$    position(1) = position(1) - real(self%index_grid(1,i_part) -1,f64)
!!$    ! Now we subtract the degree of the spline to get the index of the first spline.
!!$    self%index_grid(:,i_part) =  self%index_grid(:,i_part) - self%spline_degree
!!$    spline_val = sll_f_uniform_b_splines_at_x(self%spline_degree, position(1))
!!$    self%values_grid(:,1,i_part) = spline_val
!!$    
!!$    
!!$  end subroutine add_shape_factor_single_spline_1d



  !---------------------------------------------------------------------------!
 !> Evaluate field at at position \a position
  subroutine evaluate_field_single_spline_1d(self, position, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                              intent( in )    :: field_dofs(self%n_dofs) !< Coefficient vector for the field DoFs
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the electric fields at given position
    
    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d, index
    sll_real64 :: xi(1)

    xi(1) = (position(1) - self%domain(1,1))/self%delta_x(1)
    index = ceiling(xi(1))
    xi(1) = xi(1) - real(index-1, f64)
    index = index - self%spline_degree
    !self%spline_val = sll_f_uniform_b_splines_at_x(self%spline_degree, xi(1))
    call sll_s_uniform_b_splines_at_x(self%spline_degree, xi(1), self%spline_val)

    field_value = 0.0_f64
    do i1 = 1, self%n_span
       index1d = modulo(index+i1-2, self%n_grid(1))+1
       field_value = field_value + &
            field_dofs(index1d) *  &
            self%spline_val(i1)
    end do

  end subroutine evaluate_field_single_spline_1d



  !---------------------------------------------------------------------------!
  !> Evaluate several fields at position \a position
  subroutine evaluate_multiple_spline_1d(self, position, components, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_int32,                               intent(in)      :: components(:) !< Components of field_dofs that shall be updated
    sll_real64,                              intent( in )    :: field_dofs(:,:) !< Coefficient vector for the field DoFs
    sll_real64,                              intent(out)     :: field_value(:) !< Value(s) of the electric fields at given position
    
    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d, index
    sll_real64 :: xi(1)

    SLL_ASSERT( size(field_dofs,1) == self%n_dofs )
    SLL_ASSERT( size(field_dofs,2) == size(field_value) )

    xi(1) = (position(1) - self%domain(1,1))/self%delta_x(1)
    index = ceiling(xi(1))
    xi(1) = xi(1) - real(index-1, f64)
    index = index - self%spline_degree    
    call sll_s_uniform_b_splines_at_x(self%spline_degree, xi(1), self%spline_val)
    !self%spline_val = sll_f_uniform_b_splines_at_x(self%spline_degree, xi(1))

    field_value = 0.0_f64
    do i1 = 1, self%n_span
       index1d = modulo(index+i1-2, self%n_grid(1))+1
       field_value = field_value + &
            field_dofs(index1d,components) *  &
            self%spline_val(i1)
    end do

  end subroutine evaluate_multiple_spline_1d

  !> Destructor
  subroutine free_spline_1d(self)
    class (sll_t_particle_mesh_coupling_spline_1d), intent( inout ) :: self !< Kernel smoother object 

    deallocate(self%spline_val)
    deallocate(self%spline_val_more)
    deallocate(self%quad_xw)
    

  end subroutine free_spline_1d



  !-------------------------------------------------------------------------------------------
  !< Constructor for abstract type
  subroutine sll_s_new_particle_mesh_coupling_spline_1d_ptr(smoother, domain, n_grid, no_particles, spline_degree, smoothing_type)
    class( sll_c_particle_mesh_coupling), pointer, intent(out) :: smoother !< kernel smoother object
    sll_int32,                              intent(in)  :: n_grid(1) !< number of DoFs (spline coefficients)
    sll_real64,                             intent(in)  :: domain(2) !< x_min and x_max of the domain
    sll_int32,                              intent(in)  :: no_particles !< number of particles
    sll_int32,                              intent(in)  :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                              intent(in)  :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr


    SLL_ALLOCATE( sll_t_particle_mesh_coupling_spline_1d :: smoother , ierr)
    SLL_ASSERT( ierr == 0)
    
    select type( smoother )
    type is ( sll_t_particle_mesh_coupling_spline_1d )
       call smoother%init( domain, n_grid, no_particles, spline_degree, smoothing_type )
    end select

  end subroutine sll_s_new_particle_mesh_coupling_spline_1d_ptr

  !-------------------------------------------------------------------------------------------
  !< Constructor for abstract type
  subroutine sll_s_new_particle_mesh_coupling_spline_1d(smoother, domain, n_grid, no_particles, spline_degree, smoothing_type)
    class( sll_c_particle_mesh_coupling), allocatable, intent(out):: smoother !< kernel smoother object
    sll_int32,                                  intent(in) :: n_grid(1) !< number of DoFs (spline coefficients)
    sll_real64,                                 intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32,                                  intent(in) :: no_particles !< number of particles
    sll_int32,                                  intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                                  intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr


    allocate( sll_t_particle_mesh_coupling_spline_1d :: smoother , stat=ierr)
    SLL_ASSERT( ierr == 0)
    
    select type( smoother )
    type is ( sll_t_particle_mesh_coupling_spline_1d )
       call smoother%init( domain, n_grid, no_particles, spline_degree, smoothing_type )
    end select

  end subroutine sll_s_new_particle_mesh_coupling_spline_1d

  !> Initializer
  subroutine init_spline_1d( self, domain, n_grid, no_particles, spline_degree, smoothing_type )
    class( sll_t_particle_mesh_coupling_spline_1d), intent(out)  :: self !< kernel smoother object
    sll_int32,                               intent(in) :: n_grid(1) !< number of DoFs (spline coefficients)
    sll_real64,                              intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32,                               intent(in) :: no_particles !< number of particles
    sll_int32,                               intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                               intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr


    self%dim = 1

    ! Store grid information
    self%domain(1,:) = domain
    SLL_ALLOCATE(self%n_grid(1), ierr)
    self%n_grid = n_grid
    self%n_dofs = product(n_grid)
    self%delta_x = (self%domain(:,2)-self%domain(:,1))/real(n_grid, f64)

    ! Store basis function information
    self%no_particles = no_particles

    ! Initialize information on the spline
    self%spline_degree = spline_degree
    self%n_span = spline_degree + 1

    ! Initialize information on smoothing type
    if (smoothing_type == sll_p_collocation) then
       self%scaling = 1.0_f64/self%delta_x(1)
    elseif (smoothing_type == sll_p_galerkin) then
       self%scaling = 1.0_f64
    else
       print*, 'Smoothing Type ', smoothing_type, ' not implemented for kernel_smoother_spline_1d.'
    end if
    
    self%n_quad_points = (self%spline_degree+2)/2

    ALLOCATE( self%spline_val(self%n_span), stat = ierr)
    SLL_ASSERT( ierr == 0 )
    ALLOCATE( self%spline_val_more(self%n_span), stat = ierr )
    SLL_ASSERT( ierr == 0 )
    ALLOCATE( self%quad_xw(2,self%n_quad_points), stat = ierr )
    SLL_ASSERT( ierr == 0 )

    ! normalized Gauss Legendre points and weights
    self%quad_xw = sll_f_gauss_legendre_points_and_weights(self%n_quad_points)

  end subroutine init_spline_1d

  

end module sll_m_particle_mesh_coupling_spline_1d
