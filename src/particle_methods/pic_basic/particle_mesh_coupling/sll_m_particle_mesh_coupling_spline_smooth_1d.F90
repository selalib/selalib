!> @ingroup particle_mesh_coupling
!> @author Katharina Kormann, IPP
!> @brief Kernel smoother for 2d with splines of arbitrary degree placed on a uniform mesh.
!> @details Spline with index i starts at point i
module sll_m_particle_mesh_coupling_spline_smooth_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_gauss_legendre_integration, only : &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_particle_mesh_coupling_base_1d, only: &
       sll_p_collocation, &
       sll_c_particle_mesh_coupling_1d, &
       sll_p_galerkin

  use sll_m_particle_mesh_coupling_spline_1d, only :&
       sll_t_particle_mesh_coupling_spline_1d
  
  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base

  use sll_m_splines_pp, only: &
       sll_s_spline_pp_horner_primitive_1d
  
  
  implicit none

  public :: &
    sll_t_particle_mesh_coupling_spline_smooth_1d, &
    sll_s_new_particle_mesh_coupling_spline_smooth_1d_ptr, &
    sll_s_new_particle_mesh_coupling_spline_smooth_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
  !>  Spline kernel smoother in1d.
  type, extends(sll_t_particle_mesh_coupling_spline_1d) :: sll_t_particle_mesh_coupling_spline_smooth_1d

   contains
     procedure :: add_charge => add_charge_single_smoothened_spline_1d !> Add charge of one particle
     procedure :: add_particle_mass => add_particle_mass_spline_smooth_1d !> Add the contribution of one particle to the approximate mass matrix defined by the sum over all particles using the symmetry
     procedure :: add_particle_mass_full => add_particle_mass_full_spline_smooth_1d !> Add the contribution of one particle to the approximate mass matrix defined by the sum over all particles using the symmetry
     procedure :: evaluate => evaluate_field_single_smoothened_spline_1d !> Evaluate spline function with given coefficients
     procedure :: add_current_update_v => add_current_update_v_smoothened_spline_1d !> Add contribution of pne particle to the current density and update velocity
     
     procedure :: add_current_evaluate => add_current_evaluate_smoothened_spline_1d !> Add contribution of one particle to the current density (integrated over x)

     procedure :: evaluate_int => evaluate_int_smoothened_spline_1d !> Evaluates the integral int_{poisition_old}^{position_new} field(x) d x
     
  end type sll_t_particle_mesh_coupling_spline_smooth_1d

contains
  
   !---------------------------------------------------------------------------!
!> Add charge of one particle with smoothing
  subroutine add_charge_single_smoothened_spline_1d(self, position, marker_charge, rho_dofs)
    class( sll_t_particle_mesh_coupling_spline_smooth_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: rho_dofs(self%n_dofs) !< Coefficient vector of the charge distribution

    !local variables
    sll_int32 :: i1, q
    sll_int32 :: index1d, index
    sll_real64 :: xi(1), pos, scaling, weight
    sll_int32 :: n_quad_points
    !sll_real64, allocatable :: quad_xw(:,:)
    sll_real64 :: integ(1:self%n_span)
   
    xi(1) = (position(1) - self%domain(1))/self%delta_x
    index = floor(xi(1))+1
    xi(1) = xi(1) - real(index-1, f64)
    !index = index - self%spline_degree

    ! Quadrature formula
    n_quad_points = self%n_quadp1_points!(self%spline_degree+3)/2
    !allocate( quad_xw(2, n_quad_points) )
    !quad_xw  = sll_f_gauss_legendre_points_and_weights(n_quad_points)
    !quad_xw(1,:) = (quad_xw(1,:) + 1.0_f64)*0.5_f64
    !quad_xw(2,:) = quad_xw(2,:)*0.5_f64

    scaling = 1.0_f64 - xi(1)
    weight = marker_charge *  self%scaling * scaling 
    ! First and third interval
    do q = 1, n_quad_points
       pos = xi(1)+ scaling * self%quads1_xw(1,q)
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, pos, self%spline_val)
       ! First interval
       integ = self%quads1_xw(2,q)*self%spline_val*scaling * self%quads1_xw(1,q) 
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-3,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d) + weight * integ(i1)!&
               !(marker_charge * integ(i1)* self%scaling * scaling * self%delta_x)
       end do
       ! Third interval
       integ = self%quads1_xw(2,q)*self%spline_val*(1.0_f64-scaling * self%quads1_xw(1,q)) 
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-2,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d) +weight * integ(i1) !&
               !(marker_charge * integ(i1)* self%scaling * scaling * self%delta_x)
       end do
    end do

    weight = marker_charge *  self%scaling * xi(1) 
    ! Second and fourth interval
    do q = 1, n_quad_points
       pos = xi(1) * self%quads1_xw(1,q)
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, pos, self%spline_val)
       ! Second interval
       integ = self%quads1_xw(2,q)*self%spline_val*(scaling + xi(1) * self%quads1_xw(1,q)) 
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-2,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d)  +weight * integ(i1)
       end do
       ! Fourth interval
       integ = self%quads1_xw(2,q)*self%spline_val* xi(1)*(1.0_f64-self%quads1_xw(1,q)) 
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-1,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d)  +weight * integ(i1) 
       end do
    end do
  end subroutine add_charge_single_smoothened_spline_1d

  
  !---------------------------------------------------------------------------!
 !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting)
  subroutine add_current_update_v_smoothened_spline_1d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_smooth_1d), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
    sll_real64, intent(in)    :: position_new(self%dim) !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(:) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(:) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(:) !< Coefficients of current expansion

    ! local variables
    sll_real64 :: xi
    sll_int32  :: index_old, index_new, ind
    sll_real64 :: r_old, r_new

    ! Read out particle position and velocity
    ! Compute index_old, the index of the last DoF on the grid the particle contributes to, and r_old, its position (normalized to cell size one).
    xi = (position_old(1) - self%domain(1)) /&
         self%delta_x
    index_old = floor(xi)
    r_old = xi - real(index_old,f64)
    
    ! Compute the new box index index_new and normalized position r_old.
    xi = (position_new(1) - self%domain(1)) /&
         self%delta_x
    index_new = floor(xi)
    r_new = xi - real(index_new ,f64) 
    
    ! Evaluate the integral over x with prime function of the smoothing kernel S at old and new position
    ! We only account for the support of the kernel S
    call update_smoothened_single_spline_1d(self, [r_old], index_old, marker_charge, qoverm, &
         j_dofs, bfield_dofs, vi(2))
    call update_smoothened_single_spline_1d(self, [r_new], index_new, -marker_charge, -qoverm, &
         j_dofs, bfield_dofs, vi(2))

    ! Now, we have to add the part where one of the prime functions is already 1 but the other is not yet
    ! (when both are one they just cancel out each other and we do not have any contribution)
    ! In this part, we can reuse the old update_jv routines since the prime of the kernel is just the 1
    

    ! This version does the same based on the primitive
!!$       if (index_old == index_new) then
!!$          call self%update_jv_pp(r_old, r_new, index_old+1, marker_charge, &
!!$                  qoverm,  vi(2), j_dofs, bfield_dofs)
!!$       elseif (index_old < index_new) then
!!$          call self%update_jv_pp (r_old, 1.0_f64, index_old+1, marker_charge, &
!!$               qoverm, vi(2), j_dofs, bfield_dofs)
!!$          call self%update_jv_pp (0.0_f64, r_new, index_new+1, marker_charge, &
!!$               qoverm, vi(2), j_dofs, bfield_dofs)
!!$          do ind = index_old+2, index_new
!!$             call self%update_jv_pp (0.0_f64, 1.0_f64, ind, marker_charge, &
!!$                  qoverm, vi(2), j_dofs, bfield_dofs)
!!$          end do
!!$       else
!!$          call self%update_jv_pp (1.0_f64, r_new, index_new+1, marker_charge, qoverm, &
!!$               vi(2), j_dofs, bfield_dofs)
!!$          call self%update_jv_pp (r_old, 0.0_f64, index_old+1, marker_charge, qoverm, &
!!$               vi(2), j_dofs, bfield_dofs)
!!$          do ind = index_new+2, index_old
!!$             call self%update_jv_pp (1.0_f64, 0.0_f64, ind, marker_charge, qoverm, &
!!$                  vi(2), j_dofs, bfield_dofs)
!!$          end do
!!$       end if
       ! Version without primitive
       if (index_old == index_new) then
          if (r_old < r_new) then
             call self%update_jv(r_old, r_new, index_old+1, marker_charge, &
                  qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          else
             call self%update_jv(r_new, r_old, index_old+1, marker_charge, qoverm, &
                  -1.0_f64, vi(2), j_dofs, bfield_dofs)
          end if
       elseif (index_old < index_new) then
          call self%update_jv (r_old, 1.0_f64, index_old+1, marker_charge, &
               qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          call self%update_jv (0.0_f64, r_new, index_new+1, marker_charge, &
               qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          do ind = index_old+2, index_new
             call self%update_jv (0.0_f64, 1.0_f64, ind, marker_charge, &
                  qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          end do
       else
          call self%update_jv (r_new, 1.0_f64, index_new+1, marker_charge, qoverm, &
               -1.0_f64, vi(2), j_dofs, bfield_dofs)
          call self%update_jv (0.0_f64, r_old, index_old+1, marker_charge, qoverm, &
               -1.0_f64, vi(2), j_dofs, bfield_dofs)
          do ind = index_new+2, index_old
             call self%update_jv (0.0_f64, 1.0_f64, ind, marker_charge, qoverm, &
                  -1.0_f64, vi(2), j_dofs, bfield_dofs)
          end do
       end if

     end subroutine add_current_update_v_smoothened_spline_1d

     
 !---------------------------------------------------------------------------!
     !> Integration over the prime function of the smoothing kernel S (for v update and j accumulation)
  subroutine update_smoothened_single_spline_1d(self, position_normalized, box, marker_charge, qoverm, rho_dofs, bfield_dofs, vi)
    class( sll_t_particle_mesh_coupling_spline_smooth_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_normalized(self%dim) !< Position of the particle
    sll_int32,                                intent( in )    :: box
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                             intent(in)    :: qoverm
    sll_real64,                               intent( inout ) :: rho_dofs(self%n_dofs) !< Coefficient vector of the charge distribution
   sll_real64,                             intent(inout) :: vi
   sll_real64,                             intent(in)    :: bfield_dofs(self%n_dofs)

    !local variables
    sll_int32 :: i1, q
    sll_int32 :: index1d, index
    sll_real64 :: xi(1), pos, scaling, weight, pos_hat
    sll_int32 :: n_quad_points
    !sll_real64, allocatable :: quad_xw(:,:)
    sll_real64 :: integ(1:self%n_span)
   
    xi(1) = position_normalized(1)
    index = box +1

    ! Quadrature formula
    n_quad_points = self%n_quads2_points!(self%spline_degree+4)/2
    !allocate( quad_xw(2, n_quad_points) )
    !quad_xw  = sll_f_gauss_legendre_points_and_weights(n_quad_points)
    !quad_xw(1,:) = (quad_xw(1,:) + 1.0_f64)*0.5_f64
    !quad_xw(2,:) = quad_xw(2,:)*0.5_f64

    scaling = 1.0_f64 - xi(1)
    weight = marker_charge *  self%scaling * scaling  * self%delta_x
    ! First and third interval
    do q = 1, n_quad_points
       pos = xi(1)+ scaling * self%quads2_xw(1,q)
       pos_hat = scaling * self%quads2_xw(1,q)
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, pos, self%spline_val)
       ! First interval
       integ = self%quads2_xw(2,q)*self%spline_val* pos_hat**2*0.5_f64
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-3,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d) + weight * integ(i1)
          vi = vi - qoverm * integ(i1) * bfield_dofs(index1d) * self%delta_x * scaling 
       end do
       ! Third interval
       pos_hat = self%quads2_xw(1,q)* scaling
       integ = self%quads2_xw(2,q)*self%spline_val* (pos_hat + (1.0_f64-pos_hat**2)*0.5_f64)
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-2,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d) +weight * integ(i1)
          vi = vi - qoverm * integ(i1) * bfield_dofs(index1d) * self%delta_x * scaling 
       end do
    end do

    weight = marker_charge *  self%scaling * xi(1)  * self%delta_x
    ! Second and fourth interval
    do q = 1, n_quad_points
       pos = xi(1) * self%quads2_xw(1,q)
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, pos, self%spline_val)
       ! Second interval
       pos_hat = 1.0_f64 - xi(1) * (1.0_f64 - self%quads2_xw(1,q))
       integ = self%quads2_xw(2,q)*self%spline_val* pos_hat**2*0.5_f64 
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-2,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d)  +weight * integ(i1)
          vi = vi - qoverm * integ(i1) * bfield_dofs(index1d) * self%delta_x * xi(1)
       end do
       ! Fourth interval
       !pos_hat = xi(1) * (1.0_f64 - self%quads2_xw(1,q))
       pos_hat = scaling +  xi(1) * self%quads2_xw(1,q)
       integ = self%quads2_xw(2,q)*self%spline_val*(pos_hat + (1.0_f64-pos_hat**2)*0.5_f64)
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-1,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d)  +weight * integ(i1) 
          vi = vi - qoverm * integ(i1) * bfield_dofs(index1d) * self%delta_x * xi(1)
       end do
    end do
  end subroutine update_smoothened_single_spline_1d

  
 !---------------------------------------------------------------------------!
  !> Integration over the prime function of the smoothing kernel S (for evaluation only)
  subroutine evaluate_int_prime_smoothened_single_spline_1d(self, position_normalized, box,  field_dofs, field_int )
    class( sll_t_particle_mesh_coupling_spline_smooth_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_normalized(self%dim) !< Position of the particle
    sll_int32,                                intent( in ) :: box
    sll_real64,                             intent(inout) :: field_int
    sll_real64,                             intent(in)    :: field_dofs(self%n_dofs)

    !local variables
    sll_int32 :: i1, q
    sll_int32 :: index1d, index
    sll_real64 :: xi(1), pos, scaling, weight, pos_hat
    sll_int32 :: n_quad_points
    !sll_real64, allocatable :: quad_xw(:,:)
    sll_real64 :: integ(1:self%n_span)
   
    !xi(1) = (position(1) - self%domain(1))/self%delta_x
    !index = floor(xi(1))+1
    !xi(1) = xi(1) - real(index-1, f64)
    xi(1) = position_normalized(1)
    index = box+1
    
    ! Quadrature formula
    n_quad_points = self%n_quads2_points

    scaling = 1.0_f64 - xi(1)
    weight = self%scaling * scaling  * self%delta_x
    ! First and third interval
    do q = 1, n_quad_points
       pos = xi(1)+ scaling * self%quads2_xw(1,q)
       pos_hat = scaling * self%quads2_xw(1,q)
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, pos, self%spline_val)
       ! First interval
       integ = self%quads2_xw(2,q)*self%spline_val* pos_hat**2*0.5_f64
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-3,self%n_cells)+1
          field_int = field_int + integ(i1) * field_dofs(index1d) * self%delta_x * scaling 
       end do
       ! Third interval
       pos_hat = self%quads2_xw(1,q)* scaling
       integ = self%quads2_xw(2,q)*self%spline_val* (pos_hat + (1.0_f64-pos_hat**2)*0.5_f64)
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-2,self%n_cells)+1
          field_int = field_int + integ(i1) * field_dofs(index1d) * self%delta_x * scaling 
       end do
    end do

    weight =  self%scaling * xi(1)  * self%delta_x
    ! Second and fourth interval
    do q = 1, n_quad_points
       pos = xi(1) * self%quads2_xw(1,q)
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, pos, self%spline_val)
       ! Second interval
       pos_hat = 1.0_f64 - xi(1) * (1.0_f64 - self%quads2_xw(1,q))
       integ = self%quads2_xw(2,q)*self%spline_val* pos_hat**2*0.5_f64 
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-2,self%n_cells)+1
          field_int = field_int + integ(i1) * field_dofs(index1d) * self%delta_x * xi(1)
       end do
       ! Fourth interval
       pos_hat = scaling +  xi(1) * self%quads2_xw(1,q)
       integ = self%quads2_xw(2,q)*self%spline_val*(pos_hat + (1.0_f64-pos_hat**2)*0.5_f64)
    
       do i1 = 1, self%n_span
          index1d = modulo(index- self%spline_degree+i1-1,self%n_cells)+1
          field_int = field_int +  integ(i1) * field_dofs(index1d) * self%delta_x * xi(1)
       end do
    end do
    
  end subroutine evaluate_int_prime_smoothened_single_spline_1d

  !> Helper function to evaluate the integrated efield in one interval (without smoothing kernel)
  !> (similar to update_jv).
  subroutine evaluate_int_simple_smoothened_spline_1d(self, lower, upper, index, sign, field_dofs, field_int)
   class(sll_t_particle_mesh_coupling_spline_1d), intent(inout) :: self !< time splitting object 
   sll_real64,                             intent(in)    :: lower
   sll_real64,                             intent(in)    :: upper
   sll_int32,                              intent(in)    :: index
   sll_real64,                             intent(in)    :: sign
   sll_real64,                             intent(in)    :: field_dofs(self%n_dofs)
   sll_real64,                             intent(inout) :: field_int

   !Local variables
   sll_int32  :: ind, i_grid, i_mod, n_cells, j
   sll_real64 :: c1, c2


   n_cells = self%n_cells

   c1 =  0.5_f64*(upper-lower)
   c2 =  0.5_f64*(upper+lower)

   call sll_s_uniform_bsplines_eval_basis(self%spline_degree, c1*self%quad_xw(1,1)+c2, &
        self%spline_val)
   !alternative version with horner and pp-form
   !call sll_s_spline_pp_horner_m_1d(self%spline_pp, self%spline_val, self%spline_degree, c1*self%quad_xw(1,1)+c2)
   self%spline_val = self%spline_val * (self%quad_xw(2,1)*c1)
   do j=2,self%n_quad_points
      call sll_s_uniform_bsplines_eval_basis(self%spline_degree, c1*self%quad_xw(1,j)+c2, &
           self%spline_val_more)
      !alternative version with horner and pp-form
      !call sll_s_spline_pp_horner_m_1d(self%spline_pp, self%spline_val_more, self%spline_degree,  c1*self%quad_xw(1,j)+c2)
      self%spline_val = self%spline_val + self%spline_val_more * (self%quad_xw(2,j)*c1)
   end do
   self%spline_val = self%spline_val * (sign*self%delta_x)
   
   ind = 1
   do i_grid = index - self%spline_degree , index
      i_mod = modulo(i_grid, n_cells ) + 1 
      field_int = field_int + self%spline_val(ind)*field_dofs(i_mod)
      ind = ind + 1
   end do

 end subroutine evaluate_int_simple_smoothened_spline_1d
  

  !---------------------------------------------------------------------------!
 !> Evaluate field with smoothing at at position \a position
  subroutine evaluate_field_single_smoothened_spline_1d(self, position, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_smooth_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                              intent( in )    :: field_dofs(self%n_dofs) !< Coefficient vector for the field DoFs
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the electric fields at given position
    
    !local variables
    sll_int32 :: j
    sll_int32 :: index
    sll_real64 :: dx
    sll_real64 :: field_val
    sll_real64 :: scaling
    !sll_real64 :: quad_xw(2,self%n_quadp1_points)

    !quad_xw = self%quadp1_xw
    !quad_xw(1,:) = (quad_xw(1,:) + 1.0_f64)*0.5_f64
    !quad_xw(2,:) = quad_xw(2,:)*0.5_f64


    dx = (position(1) - self%domain(1))/self%delta_x
    index = floor(dx)
    dx = dx - real(index, f64)

    field_value = 0.0_f64
    ! First interval
    scaling = 1.0_f64-dx
    do j=1,self%n_quadp1_points
       call self%evaluate_simple &
            ( [position(1)+(self%quads1_xw(1,j)*scaling-1.0_f64)*self%delta_x], field_dofs, field_val)
       field_value = field_value + self%quads1_xw(2,j)*scaling* (self%quads1_xw(1,j)*scaling) * field_val
    end do

    ! Second interval
    do j=1,self%n_quadp1_points
       call self%evaluate_simple &
            ( [position(1)+(self%quads1_xw(1,j)-1.0_f64)*dx*self%delta_x], field_dofs, field_val)
       field_value = field_value + self%quads1_xw(2,j)*dx* (self%quads1_xw(1,j)*dx+scaling) * field_val
    end do
       
    ! Third interval
    do j=1,self%n_quadp1_points
       call self%evaluate_simple &
            ( [position(1)+self%quads1_xw(1,j)*scaling*self%delta_x], field_dofs, field_val)
       field_value = field_value + self%quads1_xw(2,j)*scaling* (1.0_f64-self%quads1_xw(1,j)*scaling) * field_val
    end do
    
    ! Fourth interval
    do j=1,self%n_quadp1_points
       call self%evaluate_simple &
            ( [position(1)+(1.0_f64+dx*(self%quads1_xw(1,j)-1.0_f64))*self%delta_x], field_dofs, field_val)
       field_value = field_value + self%quads1_xw(2,j)*dx* ((1.0_f64-self%quads1_xw(1,j))*dx) * field_val
    end do
    

  end subroutine evaluate_field_single_smoothened_spline_1d

  
  !---------------------------------------------------------------------------!
 !> Evaluate field with smoothing at at position \a position (version based on primitive)
  subroutine evaluate2_field_single_smoothened_spline_1d(self, position, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_smooth_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                              intent( in )    :: field_dofs(self%n_dofs) !< Coefficient vector for the field DoFs
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the electric fields at given position
    
    !local variables
    sll_int32 :: j, i1, index1d
    sll_int32 :: index
    sll_real64 :: dx
    sll_real64 :: field_val
    sll_real64 :: scaling
    sll_int32 :: degree
    sll_int32 :: ierr
    sll_real64, allocatable :: poly_coeffs_fp(:,:)


    dx = (position(1) - self%domain(1))/self%delta_x
    index = floor(dx)
    dx = dx - real(index, f64)

    field_value = 0.0_f64

    degree = self%spline_degree
    
    SLL_ALLOCATE(poly_coeffs_fp(degree+2,degree+1), ierr)
    SLL_ASSERT( ierr == 0 )
    
    scaling = 1.0_f64-dx
    
    ! First interval
    poly_coeffs_fp(1,:) = self%spline_pp%poly_coeffs(1,:)/int(degree+2, f64)

    do j=2, degree+1
       poly_coeffs_fp(j,:) = (-dx * self%spline_pp%poly_coeffs(j-1,:) + self%spline_pp%poly_coeffs(j,:))/int(degree+3-j, f64)
    end do
    poly_coeffs_fp(degree+2,:) =  -dx * self%spline_pp%poly_coeffs(degree+1,:)

    
    call sll_s_spline_pp_horner_primitive_1d ( self%spline_val, degree+1, poly_coeffs_fp, dx )
    call sll_s_spline_pp_horner_primitive_1d ( self%spline_val_more, degree+1, poly_coeffs_fp, 1.0_f64 )
    
    
    field_value = 0.0_f64
    do i1=1, self%n_span
       index1d = modulo(index+i1-degree-2, self%n_cells)+1
       field_value = field_value + &
            field_dofs(index1d) *  &
            (self%spline_val_more(i1)-self%spline_val(i1))
    end do

    ! Second interval
    poly_coeffs_fp(1,:) = self%spline_pp%poly_coeffs(1,:)/int(degree+2, f64)

    do j=2, degree+1
       poly_coeffs_fp(j,:) = (scaling * self%spline_pp%poly_coeffs(j-1,:) + self%spline_pp%poly_coeffs(j,:))/int(degree+3-j, f64)
    end do
    poly_coeffs_fp(degree+2,:) =  scaling * self%spline_pp%poly_coeffs(degree+1,:)

    call sll_s_spline_pp_horner_primitive_1d ( self%spline_val_more, degree+1, poly_coeffs_fp, dx )
    
     do i1=1, self%n_span
       index1d = modulo(index+i1-degree-1, self%n_cells)+1
       field_value = field_value + &
            field_dofs(index1d) *  &
            self%spline_val_more(i1)
    end do

    !Third interval
    scaling = 1.0_f64+ dx
    poly_coeffs_fp(1,:) = -self%spline_pp%poly_coeffs(1,:)/int(degree+2, f64)

    do j=2, degree+1
       poly_coeffs_fp(j,:) = (scaling * self%spline_pp%poly_coeffs(j-1,:) - self%spline_pp%poly_coeffs(j,:))/int(degree+3-j, f64)
    end do
    poly_coeffs_fp(degree+2,:) =  scaling * self%spline_pp%poly_coeffs(degree+1,:)

    
    call sll_s_spline_pp_horner_primitive_1d ( self%spline_val, degree+1, poly_coeffs_fp, dx )
    call sll_s_spline_pp_horner_primitive_1d ( self%spline_val_more, degree+1, poly_coeffs_fp, 1.0_f64 )
        
    do i1=1, self%n_span
       index1d = modulo(index+i1-degree-1, self%n_cells)+1
       field_value = field_value + &
            field_dofs(index1d) *  &
            (self%spline_val_more(i1)-self%spline_val(i1))
    end do
    
    ! Fourth interval
    poly_coeffs_fp(1,:) = - self%spline_pp%poly_coeffs(1,:)/int(degree+2, f64)

    do j=2, degree+1
       poly_coeffs_fp(j,:) = (dx * self%spline_pp%poly_coeffs(j-1,:) - self%spline_pp%poly_coeffs(j,:))/int(degree+3-j, f64)
    end do
    poly_coeffs_fp(degree+2,:) =  dx * self%spline_pp%poly_coeffs(degree+1,:)

    call sll_s_spline_pp_horner_primitive_1d ( self%spline_val_more, degree+1, poly_coeffs_fp, dx )
    
     do i1=1, self%n_span
       index1d = modulo(index+i1-degree, self%n_cells)+1
       field_value = field_value + &
            field_dofs(index1d) *  &
            self%spline_val_more(i1)
    end do
    

  end subroutine evaluate2_field_single_smoothened_spline_1d



  !---------------------------------------------------------------------------!
 !> Add current for one particle and evaluate the field (according to iterative part in discrete gradient method)
  subroutine add_current_evaluate_smoothened_spline_1d (self, position_old, position_new, marker_charge, vbar, field_dofs, j_dofs, field)
    class( sll_t_particle_mesh_coupling_spline_smooth_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( in )    :: vbar !< Particle weights time charge
    sll_real64,                               intent( in    ) :: field_dofs(self%n_dofs) !< Coefficient vector of the current density
    sll_real64,                               intent( inout ) :: j_dofs(self%n_dofs) !< Coefficient vector of the current density
    sll_real64,                               intent( out   ) :: field

    sll_real64 :: fields(2)

    fields(2) = 0.0_f64
    
    call add_current_update_v_smoothened_spline_1d(self, position_old, position_new, marker_charge, &
         -1.0_f64, field_dofs,fields, j_dofs)
    field = fields(2)/vbar
    

  end subroutine add_current_evaluate_smoothened_spline_1d


  !---------------------------------------------------------------------------!
 !> Add current for one particle and evaluate the field (according to iterative part in discrete gradient method)
  subroutine evaluate_int_smoothened_spline_1d (self, position_old, position_new, vbar, field_dofs, field)
    class( sll_t_particle_mesh_coupling_spline_smooth_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: vbar !< Particle weights time charge
    sll_real64,                               intent( in    ) :: field_dofs(self%n_dofs) !< Coefficient vector of the current density
    sll_real64,                               intent( out   ) :: field

    sll_real64 :: xi
    sll_int32  :: index_old, index_new, ind
    sll_real64 :: r_old, r_new

    ! Read out particle position and velocity
    ! Compute index_old, the index of the last DoF on the grid the particle contributes to, and r_old, its position (normalized to cell size one).
       xi = (position_old(1) - self%domain(1)) /&
            self%delta_x
       index_old = floor(xi)
       r_old = xi - real(index_old,f64)

       ! Compute the new box index index_new and normalized position r_old.
       xi = (position_new(1) - self%domain(1)) /&
            self%delta_x
       index_new = floor(xi)
       r_new = xi - real(index_new ,f64)

       field = 0.0_f64
       call evaluate_int_prime_smoothened_single_spline_1d(self, [r_new], index_new,  field_dofs, field )
       field = -field;
       call evaluate_int_prime_smoothened_single_spline_1d(self, [r_old], index_old,  field_dofs, field )
       
       ! Account for the part where one of the prime function is constantly one
       if (index_old == index_new) then
          if (r_old < r_new) then
             !call self%update_jv(r_old, r_new, index_old+1, marker_charge, &
             !     qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
             call evaluate_int_simple_smoothened_spline_1d (self, r_old, r_new, index_old+1, &
                  1.0_f64, field_dofs, field)
          else
             !call self%update_jv(r_new, r_old, index_old+1, marker_charge, qoverm, &
             !     -1.0_f64, vi(2), j_dofs, bfield_dofs)
             call evaluate_int_simple_smoothened_spline_1d (self, r_new, r_old, index_old+1, &
                  -1.0_f64, field_dofs, field)
          end if
       elseif (index_old < index_new) then
          call evaluate_int_simple_smoothened_spline_1d (self, r_old, 1.0_f64, index_old+1, &
               1.0_f64, field_dofs, field)
          call evaluate_int_simple_smoothened_spline_1d (self, 0.0_f64, r_new, index_new+1, &
               1.0_f64, field_dofs, field)
          !call self%update_jv (r_old, 1.0_f64, index_old+1, marker_charge, &
          !     qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          !call self%update_jv (0.0_f64, r_new, index_new+1, marker_charge, &
          !     qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          do ind = index_old+2, index_new
             call evaluate_int_simple_smoothened_spline_1d (self, 0.0_f64, 1.0_f64, ind, &
                  1.0_f64, field_dofs, field)
             !call self%update_jv (0.0_f64, 1.0_f64, ind, marker_charge, &
             !     qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
          end do
       else
          call evaluate_int_simple_smoothened_spline_1d (self, r_new, 1.0_f64, index_new+1, &
               -1.0_f64, field_dofs, field)
          call evaluate_int_simple_smoothened_spline_1d (self, 0.0_f64, r_old, index_old+1, &
               -1.0_f64, field_dofs, field)
          !call self%update_jv (r_new, 1.0_f64, index_new+1, marker_charge, qoverm, &
          !     -1.0_f64, vi(2), j_dofs, bfield_dofs)
          !call self%update_jv (0.0_f64, r_old, index_old+1, marker_charge, qoverm, &
          !     -1.0_f64, vi(2), j_dofs, bfield_dofs)
          do ind = index_new+2, index_old
             call evaluate_int_simple_smoothened_spline_1d (self, 0.0_f64, 1.0_f64, ind, &
                  -1.0_f64, field_dofs, field)
             !call self%update_jv (0.0_f64, 1.0_f64, ind, marker_charge, qoverm, &
             !     -1.0_f64, vi(2), j_dofs, bfield_dofs)
          end do
       end if

       field = field/vbar

  end subroutine evaluate_int_smoothened_spline_1d
  
  
  !-------------------------------------------------------------------------------------------
  !< Constructor for abstract type
  subroutine sll_s_new_particle_mesh_coupling_spline_smooth_1d_ptr(smoother, domain, n_cells, no_particles, spline_degree, smoothing_type)
    class( sll_c_particle_mesh_coupling_1d), pointer, intent(out) :: smoother !< kernel smoother object
    sll_int32,                              intent(in)  :: n_cells !< number of DoFs (spline coefficients)
    sll_real64,                             intent(in)  :: domain(2) !< x_min and x_max of the domain
    sll_int32,                              intent(in)  :: no_particles !< number of particles
    sll_int32,                              intent(in)  :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                              intent(in)  :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr


    SLL_ALLOCATE( sll_t_particle_mesh_coupling_spline_smooth_1d :: smoother , ierr)
    SLL_ASSERT( ierr == 0)
    
    select type( smoother )
    type is ( sll_t_particle_mesh_coupling_spline_smooth_1d )
       call smoother%init( domain, n_cells, no_particles, spline_degree, smoothing_type )
    end select

  end subroutine sll_s_new_particle_mesh_coupling_spline_smooth_1d_ptr

  !----------------------------------------------------------------------------------
  !< Constructor for abstract type
  subroutine sll_s_new_particle_mesh_coupling_spline_smooth_1d(smoother, domain, n_cells, no_particles, spline_degree, smoothing_type)
    class( sll_c_particle_mesh_coupling_1d), allocatable, intent(out):: smoother !< kernel smoother object
    sll_int32,                                  intent(in) :: n_cells !< number of DoFs (spline coefficients)
    sll_real64,                                 intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32,                                  intent(in) :: no_particles !< number of particles
    sll_int32,                                  intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                                  intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr


    allocate( sll_t_particle_mesh_coupling_spline_smooth_1d :: smoother , stat=ierr)
    SLL_ASSERT( ierr == 0)
    
    select type( smoother )
    type is ( sll_t_particle_mesh_coupling_spline_smooth_1d )
       call smoother%init( domain, n_cells, no_particles, spline_degree, smoothing_type )
    end select

  end subroutine sll_s_new_particle_mesh_coupling_spline_smooth_1d
!---------------------------------------------------------------------------!
  !> Add particle mass for one particle (upper triangular version)
  subroutine add_particle_mass_spline_smooth_1d(self, position, marker_charge, particle_mass)
    class( sll_t_particle_mesh_coupling_spline_smooth_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution

    !local variables
    sll_int32 :: i1, column
    sll_int32 :: index1d, index
    sll_real64 :: xi(1)
    sll_real64 :: scratch(1:self%n_span+2)

    SLL_ASSERT( size(particle_mass,1) == self%n_span+2 )
    SLL_ASSERT( size(particle_mass,2) == self%n_dofs )
    
    xi(1) = (position(1) - self%domain(1))/self%delta_x
    index = ceiling(xi(1))
    xi(1) = xi(1) - real(index-1, f64)
    index = index - self%spline_degree-1

    scratch = 0.0_f64
    call add_charge_box(self, xi, scratch)


    do i1 = 1, self%n_span+2
       index1d = modulo(index+i1-2,self%n_cells)+1
       do column = i1, self%n_span+2
          particle_mass(column-i1+1, index1d) = particle_mass( column-i1+1, index1d) + &
               marker_charge * scratch(i1) * scratch(column)
         ! Note: No scaling since Galerkin function
       end do
    end do
    

  end subroutine add_particle_mass_spline_smooth_1d
!---------------------------------------------------------------------------!
  !> Add particle mass for one particle (full matrix)
  !> WARNING: THIS FUNCTION SEEMS TO HAVE BUGS
  subroutine add_particle_mass_full_spline_smooth_1d(self, position, marker_charge, particle_mass)
    class( sll_t_particle_mesh_coupling_spline_smooth_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution

    !local variables
    sll_int32 :: i1, column, ind
    sll_int32 :: index1d, index
    sll_real64 :: xi(1)
    sll_real64 :: scratch(1:self%n_span+2)

    SLL_ASSERT( size(particle_mass,1) == 2*self%n_span+3 )
    SLL_ASSERT( size(particle_mass,2) == self%n_dofs )

    xi(1) = (position(1) - self%domain(1))/self%delta_x
    index = ceiling(xi(1))
    xi(1) = xi(1) - real(index-1, f64)
    index = index - self%spline_degree-1

    scratch = 0.0_f64
    call add_charge_box(self, xi, scratch)


    do i1 = 1, self%n_span+2
       index1d = modulo(index+i1-2,self%n_cells)+1
       ind = 1+(self%n_span+2-i1)
       do column = 1, self%n_span+2
          particle_mass(ind, index1d) = particle_mass( ind, index1d) + &
               marker_charge * scratch(i1) * scratch(column)
          ! Note: No scaling since Galerkin function
          ind = ind+1
       end do
    end do
    

  end subroutine add_particle_mass_full_spline_smooth_1d

  
     !---------------------------------------------------------------------------!
!> Add charge of one particle with smoothing
  subroutine add_charge_box(self, position, rho_dofs)
    class( sll_t_particle_mesh_coupling_spline_smooth_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( inout ) :: rho_dofs(self%n_span+2) !< Coefficient vector of the charge distribution

    !local variables
    sll_int32 :: i1, q
    sll_int32 :: index1d, index
    sll_real64 :: xi(1), pos, scaling, weight
    sll_int32 :: n_quad_points
    sll_real64 :: integ(1:self%n_span)
   
    xi(1) = position(1)

    ! Quadrature formula
    n_quad_points = self%n_quadp1_points

    scaling = 1.0_f64 - xi(1)
    weight = scaling 
    ! First and third interval
    do q = 1, n_quad_points
       pos = xi(1)+ scaling * self%quads1_xw(1,q)
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, pos, self%spline_val)
       ! First interval
       integ = self%quads1_xw(2,q)*self%spline_val*scaling * self%quads1_xw(1,q) 
    
       do i1 = 1, self%n_span
          index1d = i1
          !index1d = modulo(index- self%spline_degree+i1-3,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d) + weight * integ(i1)!&
               !(marker_charge * integ(i1)* self%scaling * scaling * self%delta_x)
       end do
       ! Third interval
       integ = self%quads1_xw(2,q)*self%spline_val*(1.0_f64-scaling * self%quads1_xw(1,q)) 
    
       do i1 = 1, self%n_span
          index1d = i1+1
          !index1d = modulo(index- self%spline_degree+i1-2,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d) +weight * integ(i1) !&
               !(marker_charge * integ(i1)* self%scaling * scaling * self%delta_x)
       end do
    end do

    weight = xi(1) 
    ! Second and fourth interval
    do q = 1, n_quad_points
       pos = xi(1) * self%quads1_xw(1,q)
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, pos, self%spline_val)
       ! Second interval
       integ = self%quads1_xw(2,q)*self%spline_val*(scaling + xi(1) * self%quads1_xw(1,q)) 
    
       do i1 = 1, self%n_span
          index1d = i1+1
          !index1d = modulo(index- self%spline_degree+i1-2,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d)  +weight * integ(i1) !+&
               !(marker_charge * integ(i1)* self%scaling * xi(1) * self%delta_x)
       end do
       ! Fourth interval
       integ = self%quads1_xw(2,q)*self%spline_val* xi(1)*(1.0_f64-self%quads1_xw(1,q)) 
    
       do i1 = 1, self%n_span
          index1d = i1+2
          !index1d = modulo(index- self%spline_degree+i1-1,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d)  +weight * integ(i1) !+&
               !(marker_charge * integ(i1)* self%scaling * xi(1) * self%delta_x)
       end do
    end do

  end subroutine add_charge_box
  
end module sll_m_particle_mesh_coupling_spline_smooth_1d
