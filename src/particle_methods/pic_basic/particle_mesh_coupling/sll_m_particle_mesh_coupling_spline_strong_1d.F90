!> @ingroup particle_mesh_coupling
!> @author Katharina Kormann, IPP
!> @brief Kernel smoother for 1d with splines of arbitrary degree placed on a uniform mesh. This version is for a formulation of the Maxwell's equation with strong Ampere.
!> @details Spline with index i starts at point i
!> Reference: Campos Pinto, Kormann, Sonnendrücker: Variational Framework for Structure-Preserving Electromagnetic Particle-In-Cell Methods, arXiv 2101.09247, 2021.
module sll_m_particle_mesh_coupling_spline_strong_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_low_level_bsplines, only: &
    sll_s_uniform_bsplines_eval_basis

  use sll_m_gauss_legendre_integration, only : &
    sll_f_gauss_legendre_points_and_weights

  use sll_m_particle_mesh_coupling_base_1d, only: &
    sll_p_collocation, &
    sll_c_particle_mesh_coupling_1d, &
    sll_p_galerkin

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_splines_pp

  implicit none

  public :: &
       sll_t_particle_mesh_coupling_spline_strong_1d, &
       sll_s_new_particle_mesh_coupling_spline_strong_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Basic type of a kernel smoother used for PIC simulations
  type, extends(sll_c_particle_mesh_coupling_1d) :: sll_t_particle_mesh_coupling_spline_strong_1d
     ! Member variables form abstract type
     logical                 :: integ
     logical                 :: eval_same = .true. !< eval same is true if the solver degree is the same order as the spline degree
     sll_int32               :: index_shift !< index shift
     
     ! Information about the 1d mesh
     sll_real64 :: delta_x  !< Value of grid spacing along both directions.
     sll_real64 :: domain(1,2) !< Definition of the domain: domain(1,1) = x1_min  domain(1,2) = x1_max
     sll_int32  :: n_span !< Number of intervals where spline non zero (spline_degree + 1)
     type(sll_t_spline_pp_1d) :: spline_pp !< 1d pp-spline

     
     sll_real64, allocatable :: spline_val(:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_val_more(:) !< scratch data for spline evaluation
     
     sll_real64, allocatable :: spline_val_prim(:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_val_prim_old(:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_val_prim_new(:) !< scratch data for spline evaluation
     
     sll_int32  :: n_quad_points !< Number of quadrature points
     sll_real64, allocatable :: quad_xw(:,:) !< quadrature weights and points
     
   contains

     procedure    :: add_charge => add_charge_single_spline_strong_1d !> Add the contribution of one particle to the charge density
     procedure    :: add_particle_mass => add_particle_mass_spline_strong_1d !> Add the contribution of one particle to the approximate mass matrix defined by the sum over all particles using the symmetry
     procedure    :: add_particle_mass_full  => add_particle_mass_spline_strong_1d  !> Add the contribution of one particle to the approximate mass matrix defined by the sum over all particles
     procedure    :: evaluate => evaluate_field_single_spline_strong_1d !> Evaluate spline function with given coefficients
     procedure    :: evaluate_multiple => evaluate_multiple_spline_strong_1d !> Evaluate multiple spline functions with given coefficients
     procedure    :: add_current_update_v => add_current_update_v_spline_strong_1d_quadrature !> Add contribution of pne particle to the current density and update velocity
     !procedure    :: add_current_update_v_pp => add_current_update_v_pp_spline_strong_1d !> Add contribution of pne particle to the current density and update velocity 
     procedure    :: add_current => add_current_spline_strong_1d !> Add contribution of one particle to the current density (integrated over x)
     procedure    :: add_current_evaluate => add_current_evaluate_spline_strong_1d_quadrature

     procedure           :: init => init_spline_strong_1d !> Constructor
     procedure           :: free => free_spline_strong_1d !< Destructor


  end type sll_t_particle_mesh_coupling_spline_strong_1d

contains

  ! TODOTODOTODOTODO

  ! Check the /delta_x scaling in the add_charge function
  subroutine helper_normalized_position ( self, position, index, axi)
    class( sll_t_particle_mesh_coupling_spline_strong_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_int32, intent(out) :: index
    sll_real64, intent(out) :: axi(self%dim)

    sll_real64 :: xi(1)

    xi(1) = (position(1) - self%domain(1,1))/self%delta_x
    index = floor(xi(1))+1
    xi(1) = xi(1) - real(index-1, f64)

    ! Old version before changed to option of evaluating in other points than the grid point
!!$    if ( modulo(self%spline_degree, 2) == 0 ) then
!!$       if (xi(1) < 0.5_f64 ) then
!!$          axi(1) = 0.5_f64 - xi(1)
!!$          index = index - self%spline_degree/2 
!!$       else
!!$          axi(1) = 1.5_f64 - xi(1)
!!$          index = index - self%spline_degree/2 + 1
!!$       end if
!!$    else
!!$       axi(1) = 1.0_f64 - xi(1)
!!$       index = index - (self%n_span)/2+1
!!$    end if
    if ( self%eval_same .eqv. .true. ) then ! case odd spline with evaluation at grid points or even spline with evaluation at mid points
       axi(1) = 1.0_f64 - xi(1)
       index = index - self%index_shift!(self%spline_degree)/2
    else ! case even spline with evaluation at grid points or odd spline with evaluation at mid points
       if (xi(1) < 0.5_f64 ) then
          axi(1) = 0.5_f64 - xi(1)
          index = index - self%index_shift!(self%n_span-1)/2 
       else
          axi(1) = 1.5_f64 - xi(1)
          index = index - self%index_shift + 1!(self%n_span-1)/2 + 1
       end if
    end if

      ! print*, axi, index
  !     print*,  (self%spline_degree-1)/2,  (self%spline_degree)/2
       
  end subroutine helper_normalized_position
  
  
 !---------------------------------------------------------------------------!
!> Add charge of one particle
  subroutine add_charge_single_spline_strong_1d(self, position, marker_charge, rho_dofs)
    class( sll_t_particle_mesh_coupling_spline_strong_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: rho_dofs(self%n_dofs) !< Coefficient vector of the charge distribution

    ! NOTE THAT THIS IS THE VARIANT FOR INT=TRUE ONLY

    
    !local variables
    sll_int32 :: i1, i
    sll_int32 :: index1d, index
    sll_real64 :: xi(1), axi(1)

    call helper_normalized_position( self, position, index, axi )

    if ( self%integ .eqv. .false. ) then
       ! Evaluate the spline function
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree, axi(1), &
            self%spline_val )
       
       do i1 = 1, self%n_span
          index1d = modulo(index+i1-2, self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d) + &
               marker_charge * self%spline_val(self%n_span-i1+1)/ self%delta_x 
       end do
    else
       
       ! Old version with primitive formula
       ! Evaluate the primitive function
       do i=1, self%n_span
          self%spline_val_prim(self%n_span-i+2) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, axi(1), i) 
       end do
       
       ! fill the right-hand-side
       do i=1, self%n_span+1
          index1d = modulo(index+i-2,self%n_cells)+1
          rho_dofs(index1d) = rho_dofs(index1d) + &
               marker_charge * ( self%spline_val_prim(i+1) - self%spline_val_prim(i) )/ self%delta_x 
       end do
    end if
    
  end subroutine add_charge_single_spline_strong_1d

  !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting). Implementation based on quadrature rather than the primitive
  subroutine add_current_update_v_spline_strong_1d_quadrature (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_strong_1d), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
    sll_real64, intent(in)    :: position_new(self%dim) !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(:) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(:) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(:) !< Coefficients of current expansion

    !local variables
    sll_int32 :: ind
    sll_int32 ::  index_new, index_old
    sll_real64 :: axi_new(1), axi_old(1)

    call helper_normalized_position( self, position_old, index_old, axi_old )
    call helper_normalized_position( self, position_new, index_new, axi_new )

    if (index_old == index_new) then      
       call current_v_local( self, axi_old(1), axi_new(1), index_old, marker_charge, qoverm, vi(2), j_dofs, bfield_dofs )
    elseif(index_old < index_new ) then
       call current_v_local( self, axi_old(1), 0.0_f64, index_old, marker_charge, qoverm, vi(2), j_dofs, bfield_dofs )
       do ind = index_old+1, index_new-1
          call current_v_local( self, 1.0_f64, 0.0_f64, ind , marker_charge, qoverm, vi(2), j_dofs, bfield_dofs)
       end do
       call current_v_local( self, 1.0_f64, axi_new(1), index_new, marker_charge, qoverm, vi(2), j_dofs, bfield_dofs )
    else
       call current_v_local( self, axi_old(1), 1.0_f64, index_old, marker_charge, qoverm, vi(2), j_dofs, bfield_dofs )
       do ind = index_new+1, index_old-1
          call current_v_local( self, 0.0_f64, 1.0_f64, ind , marker_charge, qoverm, vi(2), j_dofs, bfield_dofs)
       end do
       call current_v_local( self, 0.0_f64, axi_new(1), index_new, marker_charge, qoverm, vi(2), j_dofs, bfield_dofs )  
    end if
    

  end subroutine add_current_update_v_spline_strong_1d_quadrature


  subroutine current_v_local( self, upper, lower, box, marker_charge, qoverm,  vi, j_dofs, bfield_dofs)
   class(sll_t_particle_mesh_coupling_spline_strong_1d), intent(inout) :: self !< time splitting object 
   sll_real64,                             intent(in)    :: lower
   sll_real64,                             intent(in)    :: upper
   sll_int32,                              intent(in)    :: box
   sll_real64,                             intent(in)    :: marker_charge
   sll_real64,                             intent(in)    :: qoverm
   sll_real64,                             intent(inout) :: vi
   sll_real64,                             intent(in)    :: bfield_dofs(self%n_dofs)
   sll_real64,                             intent(inout) :: j_dofs(self%n_dofs)
   
   !Local variables
   sll_int32  :: ind, i_grid, i_mod, n_cells, j, i
   sll_real64 :: c1, c2
   n_cells = self%n_cells

    c1 =  0.5_f64*(upper-lower)
    c2 =  0.5_f64*(upper+lower)

    call sll_s_uniform_bsplines_eval_basis(self%spline_degree, c1*self%quad_xw(1,1)+c2, &
         self%spline_val)
    self%spline_val = self%spline_val * (self%quad_xw(2,1)*c1)
    do j=2,self%n_quad_points
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, c1*self%quad_xw(1,j)+c2, &
            self%spline_val_more)
       self%spline_val = self%spline_val + self%spline_val_more * (self%quad_xw(2,j)*c1)
    end do
!    self%spline_val = self%spline_val * (self%delta_x)
!    do i=1, self%n_span
!       self%spline_val_prim_old(self%n_span-i+1) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, upper, i) 
!       self%spline_val_prim_new(self%n_span-i+1) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, lower, i) 
!    end do
    
    ind = self%spline_degree+1
    do i_grid = box, box+self%spline_degree
       i_mod = modulo(i_grid-1, n_cells ) +1
       j_dofs(i_mod) = j_dofs(i_mod) + marker_charge*self%spline_val(ind)
       vi = vi - qoverm * self%spline_val(ind)*bfield_dofs(i_mod)
       ind = ind-1
    end do


  end subroutine current_v_local

  
   !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting)
  subroutine add_current_update_v_spline_strong_1d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_strong_1d), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
    sll_real64, intent(in)    :: position_new(self%dim) !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(:) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(:) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(:) !< Coefficients of current expansion


    
    
    !local variables
    sll_int32 :: i, i0, in, i1, diff
    sll_int32 :: index1d, index_new, index_old
    sll_real64 :: axi_new(1), axi_old(1)

    call helper_normalized_position( self, position_old, index_old, axi_old )
    call helper_normalized_position( self, position_new, index_new, axi_new )
    
    ! Evaluate the primitive function
    do i=1, self%n_span
       self%spline_val_prim_old(self%n_span-i+1) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, axi_old(1), i) 
       self%spline_val_prim_new(self%n_span-i+1) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, axi_new(1), i) 
    end do
    
    !call sll_s_spline_pp_horner_primitive_1d(self%spline_val, self%spline_degree, self%spline_pp%poly_coeffs_fp, lower)
    !call sll_s_spline_pp_horner_primitive_1d(self%spline_val, self%spline_degree, self%spline_pp%poly_coeffs_fp, lower) 
    !print*, axi_old(1), axi_new(1)
    !print*, self%spline_val_prim_new
    !print*, self%spline_val_prim_old
    !stop
   ! do i=1, self%n_span
   !    self%spline_val_prim_old(i) =  sll_f_spline_pp_horner_1d(self%spline_degree, self%spline_pp%poly_coeffs, axi_old(1), i) / self%delta_x * (position_new(1) - position_old(1))
   !    self%spline_val_prim_new(i) =  -sll_f_spline_pp_horner_1d(self%spline_degree, self%spline_pp%poly_coeffs, axi_new(1), i) / self%delta_x * (position_new(1)-position_old(1))
   ! end do

    if (index_old<index_new) then
       diff = index_new-index_old

       do i=0,min(diff-1, self%spline_degree)
          i1 = modulo(index_old+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) + self%spline_val_prim_old(i+1)* marker_charge
          vi (2) = vi(2) - qoverm * bfield_dofs(i1) * ( self%spline_val_prim_old(i+1))
       end do
       do i = self%spline_degree+1, diff-1
          i1 = modulo(index_old+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) +  marker_charge
          vi (2) = vi(2) - qoverm * bfield_dofs(i1) 
       end do
       do i = diff, self%spline_degree
          i1 = modulo(index_old+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) + ( self%spline_val_prim_old(i+1) - self%spline_val_prim_new(i+1-diff)) &
               * marker_charge
          vi(2) = vi(2) - qoverm * bfield_dofs(i1) * ( self%spline_val_prim_old(i+1) - self%spline_val_prim_new(i+1-diff))
       end do
       do i=max(self%spline_degree+1, diff), self%spline_degree+diff
          i1 = modulo(index_old+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) + ( 1.0_f64 - self%spline_val_prim_new(i+1-diff))*marker_charge
          vi(2) = vi(2) - qoverm * bfield_dofs(i1) * ( 1.0_f64 - self%spline_val_prim_new(i+1-diff))
       end do
    else
       diff = index_old-index_new
       do i=0,min(diff-1, self%spline_degree)
          i1 = modulo(index_new+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) - self%spline_val_prim_new(i+1)*marker_charge
          vi(2) = vi(2) - qoverm * bfield_dofs(i1) * ( - self%spline_val_prim_new(i+1))
       end do
       do i=self%spline_degree+1, diff-1
          i1 = modulo(index_new+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) - marker_charge
          vi(2) = vi(2) + qoverm * bfield_dofs(i1) 
       end do
          
       do i = diff, self%spline_degree
          i1 = modulo(index_new+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) +( self%spline_val_prim_old(i+1-diff) - self%spline_val_prim_new(i+1)) &
               * marker_charge
          vi(2) = vi(2) - qoverm * bfield_dofs(i1) * ( self%spline_val_prim_old(i+1-diff) - self%spline_val_prim_new(i+1))
       end do
       do i=max(self%spline_degree+1, diff), self%spline_degree+diff
          i1 = modulo(index_new+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) + (-1.0_f64 + self%spline_val_prim_old(i+1-diff)) &
               * marker_charge
          vi(2) = vi(2) - qoverm * bfield_dofs(i1) * ( self%spline_val_prim_old(i+1-diff) - 1.0_f64 )
       end do
    end if
    
  end subroutine add_current_update_v_spline_strong_1d


    subroutine evaluate_field_single_spline_strong_1d(self, position, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_strong_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                              intent( in )    :: field_dofs(self%n_dofs) !< Coefficient vector for the field DoFs
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the electric fields at given position
    
    !local variables
    sll_int32 :: i1, i
    sll_int32 :: index1d, index
    sll_real64 :: axi(1)

    
    call helper_normalized_position( self, position, index, axi )
    
    if ( self%integ .eqv. .true. ) then
       ! Evaluate the primitive function
       do i=1, self%n_span
          self%spline_val_prim(self%n_span-i+2) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, axi(1), i) 
       end do

       field_value = 0.0_f64
       do i1 = 1, self%n_span+1
          index1d = modulo(index+i1-2, self%n_cells)+1
          field_value = field_value + &
               field_dofs(index1d) *  &
               (self%spline_val_prim(i1+1)-self%spline_val_prim(i1)) / &
               self%delta_x
       end do
    else
       ! Evaluate the spline function
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree, axi(1), &
            self%spline_val )

       field_value = 0.0_f64
       do i1 = 1, self%n_span
          index1d = modulo(index+i1-2, self%n_cells)+1
          field_value = field_value + &
               field_dofs(index1d) *  &
               self%spline_val(self%n_span-i1+1) / &
               self%delta_x
       end do
         
       
    end if
    
  end subroutine evaluate_field_single_spline_strong_1d


    !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting). Implementation based on quadrature rather than the primitive
  subroutine add_current_evaluate_spline_strong_1d_quadrature (self, position_old, position_new, marker_charge, vbar, field_dofs, j_dofs, field)
    class(sll_t_particle_mesh_coupling_spline_strong_1d), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
    sll_real64, intent(in)    :: position_new(self%dim) !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent( in )    :: vbar !< Particle weights time charge
    sll_real64, intent(in)    :: field_dofs(self%n_dofs) !< Coefficient of field expansion
    sll_real64, intent(inout) :: j_dofs(self%n_dofs) !< Coefficients of current expansion
    sll_real64, intent( out ) :: field !< Evaluated field (integrated)

    !local variables
    sll_int32 :: ind
    sll_int32 ::  index_new, index_old
    sll_real64 :: axi_new(1), axi_old(1)


    field = 0.0_f64
    call helper_normalized_position( self, position_old, index_old, axi_old )
    call helper_normalized_position( self, position_new, index_new, axi_new )

    if (index_old == index_new) then      
       call current_eval_local( self, axi_old(1), axi_new(1), index_old, marker_charge, field, j_dofs, field_dofs )
    elseif(index_old < index_new ) then
       call current_eval_local( self, axi_old(1), 0.0_f64, index_old, marker_charge, field, j_dofs, field_dofs )
       do ind = index_old+1, index_new-1
          call current_eval_local( self, 1.0_f64, 0.0_f64, ind , marker_charge, field, j_dofs, field_dofs)
       end do
       call current_eval_local( self, 1.0_f64, axi_new(1), index_new, marker_charge, field, j_dofs, field_dofs )
    else
       call current_eval_local( self, axi_old(1), 1.0_f64, index_old, marker_charge, field, j_dofs, field_dofs )
       do ind = index_new+1, index_old-1
          call current_eval_local( self, 0.0_f64, 1.0_f64, ind , marker_charge, field, j_dofs, field_dofs)
       end do
       call current_eval_local( self, 0.0_f64, axi_new(1), index_new, marker_charge, field, j_dofs, field_dofs )  
    end if

    field = field/vbar

  end subroutine add_current_evaluate_spline_strong_1d_quadrature

  

  subroutine current_eval_local( self, upper, lower, box, marker_charge, field,  j_dofs, field_dofs)
   class(sll_t_particle_mesh_coupling_spline_strong_1d), intent(inout) :: self !< time splitting object 
   sll_real64,                             intent(in)    :: lower
   sll_real64,                             intent(in)    :: upper
   sll_int32,                              intent(in)    :: box
   sll_real64,                             intent(in)    :: marker_charge
   sll_real64,                             intent(inout) :: field
   sll_real64,                             intent(in)    :: field_dofs(self%n_dofs)
   sll_real64,                             intent(inout) :: j_dofs(self%n_dofs)
   
   !Local variables
   sll_int32  :: ind, i_grid, i_mod, n_cells, j, i
   sll_real64 :: c1, c2
   n_cells = self%n_cells

    c1 =  0.5_f64*(upper-lower)
    c2 =  0.5_f64*(upper+lower)

    call sll_s_uniform_bsplines_eval_basis(self%spline_degree, c1*self%quad_xw(1,1)+c2, &
         self%spline_val)
    self%spline_val = self%spline_val * (self%quad_xw(2,1)*c1)
    do j=2,self%n_quad_points
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, c1*self%quad_xw(1,j)+c2, &
            self%spline_val_more)
       self%spline_val = self%spline_val + self%spline_val_more * (self%quad_xw(2,j)*c1)
    end do
!    self%spline_val = self%spline_val * (self%delta_x)
!    do i=1, self%n_span
!       self%spline_val_prim_old(self%n_span-i+1) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, upper, i) 
!       self%spline_val_prim_new(self%n_span-i+1) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, lower, i) 
!    end do
    
    ind = self%spline_degree+1
    do i_grid = box, box+self%spline_degree
       i_mod = modulo(i_grid-1, n_cells ) +1
       j_dofs(i_mod) = j_dofs(i_mod) + marker_charge*self%spline_val(ind)
       field = field + field_dofs(i_mod)* self%spline_val(ind)
       ind = ind-1
    end do


  end subroutine current_eval_local

  

  
  subroutine add_current_evaluate_spline_strong_1d(self, position_old, position_new, marker_charge, vbar, field_dofs, j_dofs, field)
    class( sll_t_particle_mesh_coupling_spline_strong_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( in )    :: vbar !< Particle weights time charge
    sll_real64,                               intent( in    ) :: field_dofs(self%n_dofs) !< Coefficient vector of the current density
    sll_real64,                               intent( inout ) :: j_dofs(self%n_dofs) !< Coefficient vector of the current density
    sll_real64,                               intent( out   ) :: field
 !local variables
    sll_int32 :: i, i0, in, i1, diff
    sll_int32 :: index1d, index_new, index_old
    sll_real64 :: axi_new(1), axi_old(1)

    call helper_normalized_position( self, position_old, index_old, axi_old )
    call helper_normalized_position( self, position_new, index_new, axi_new )
    
    ! Evaluate the primitive function
    do i=1, self%n_span
       self%spline_val_prim_old(self%n_span-i+1) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, axi_old(1), i) 
       self%spline_val_prim_new(self%n_span-i+1) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, axi_new(1), i) 
    end do

    field = 0.0_f64
    
    if (index_old<index_new) then
       diff = index_new-index_old

       do i=0,min(diff-1, self%spline_degree)
          i1 = modulo(index_old+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) + self%spline_val_prim_old(i+1)* marker_charge
          field = field + field_dofs(i1) * ( self%spline_val_prim_old(i+1))
       end do
       do i = self%spline_degree+1, diff-1
          i1 = modulo(index_old+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) +  marker_charge
          field = field + field_dofs(i1) 
       end do
       do i = diff, self%spline_degree
          i1 = modulo(index_old+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) + ( self%spline_val_prim_old(i+1) - self%spline_val_prim_new(i+1-diff)) &
               * marker_charge
          field = field + field_dofs(i1) * ( self%spline_val_prim_old(i+1) - self%spline_val_prim_new(i+1-diff))
       end do
       do i=max(self%spline_degree+1, diff), self%spline_degree+diff
          i1 = modulo(index_old+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) + ( 1.0_f64 - self%spline_val_prim_new(i+1-diff))*marker_charge
          field = field + field_dofs(i1) * ( 1.0_f64 - self%spline_val_prim_new(i+1-diff))
       end do
    else
       diff = index_old-index_new
       do i=0,min(diff-1, self%spline_degree)
          i1 = modulo(index_new+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) - self%spline_val_prim_new(i+1)*marker_charge
          field = field + field_dofs(i1) * ( - self%spline_val_prim_new(i+1))
       end do
       do i=self%spline_degree+1, diff-1
          i1 = modulo(index_new+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) - marker_charge
          field = field - field_dofs(i1) 
       end do
          
       do i = diff, self%spline_degree
          i1 = modulo(index_new+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) +( self%spline_val_prim_old(i+1-diff) - self%spline_val_prim_new(i+1)) &
               * marker_charge
         field = field + field_dofs(i1) * ( self%spline_val_prim_old(i+1-diff) - self%spline_val_prim_new(i+1))
       end do
       do i=max(self%spline_degree+1, diff), self%spline_degree+diff
          i1 = modulo(index_new+i-1, self%n_cells) + 1
          j_dofs(i1) = j_dofs(i1) + (-1.0_f64 + self%spline_val_prim_old(i+1-diff)) &
               * marker_charge
          field = field + field_dofs(i1) * ( self%spline_val_prim_old(i+1-diff) - 1.0_f64 )
       end do
    end if

    field = field/vbar
    
  end subroutine add_current_evaluate_spline_strong_1d
  

  !> Constructor
  subroutine init_spline_strong_1d( self, domain, n_cells, spline_degree, integ, eval_grid_points )
    class (sll_t_particle_mesh_coupling_spline_strong_1d), intent( inout ) :: self !< Particle mesh coupling object
    sll_int32,                               intent(in) :: n_cells !< number of DoFs (spline coefficients)
    sll_real64,                              intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32,                               intent(in) :: spline_degree !< Degree of smoothing kernel spline
    logical, intent( in ) :: integ
    logical, intent( in ) :: eval_grid_points
    
    !local variables
    sll_int32 :: ierr

    self%integ = .false.!integ
    self%dim = 1

    ! Store grid information
    self%domain(1,:) = domain
    self%n_cells = n_cells
    self%n_dofs = n_cells
    self%delta_x = (self%domain(1,2)-self%domain(1,1))/real(n_cells, f64)


    ! Initialize information on the spline
    self%spline_degree = spline_degree
    self%n_span = spline_degree + 1
    
!    self%n_quad_points = (self%spline_degree+2)/2

    ALLOCATE( self%spline_val(self%n_span), stat = ierr)
    SLL_ASSERT( ierr == 0 )
    ALLOCATE( self%spline_val_more(self%n_span), stat = ierr)
    SLL_ASSERT( ierr == 0 )
    ALLOCATE( self%spline_val_prim(self%n_span+2), stat = ierr )
    SLL_ASSERT ( ierr == 0 )
    self%spline_val_prim(1) = 0.0_f64
    self%spline_val_prim(self%n_span+2) = 1.0_f64
    ALLOCATE( self%spline_val_prim_old(self%n_span), stat = ierr )
    SLL_ASSERT ( ierr == 0 )
    ALLOCATE( self%spline_val_prim_new(self%n_span), stat = ierr )
    SLL_ASSERT ( ierr == 0 )

    if ( self%spline_degree < 6 ) then
       call sll_s_spline_pp_init_1d( self%spline_pp, spline_degree, n_cells)
    end if
    !> Now, we have to check the various options
    !> If we evaluate the integral, we use the following formula
    !> int_{x_i-h}^{x_i} N_{x_p}^{p-1}(x) d x = N_{x_p}^p (x_i - h/2)
    !> (Note that if we have the formula for N_j instead of N_{x_p} we evaluate at the upper bound but
    !> N_j^{p-1} and N_j^{p} start at the same point, i.e. the center is shifted by -h/2 which here
    !> appears in the argument instead).
    !> Hence, we need to evaluate at the grid point i instead of the midpoint i and
    !> at the mid point i-1 instead of the grid point i if we evaluate the integral.
    !> Now we compute the index shift that we have in each case from the index i
    !> we the particle is located. Note that in the case of eval_same = .false.
    !> the shift depends on whether the particle is located in the first or second half
    !> of the cell. The index shift is for the case that it is in the first half. Otherwise, we
    !> have to add one.
    if ( integ .eqv. .false. ) then
       if ( eval_grid_points .eqv. .true. ) then
          if ( modulo( self%spline_degree, 2 ) == 0 ) then
             self%eval_same = .false.
             self%index_shift = (self%spline_degree)/2
          else
             self%eval_same = .true.
             self%index_shift = (self%spline_degree-1)/2
          end if
       else
          if ( modulo( self%spline_degree, 2 ) == 0 ) then
             self%eval_same = .true.
             self%index_shift = (self%spline_degree)/2
          else
             self%eval_same = .false.
             self%index_shift = (self%spline_degree+1)/2
          end if
       end if
    else
       if ( eval_grid_points .eqv. .true. ) then
          if ( modulo( self%spline_degree, 2 ) == 0 ) then
             self%eval_same = .true.
             self%index_shift = (self%spline_degree-2)/2
          else
             self%eval_same = .false.
             self%index_shift = (self%spline_degree-1)/2
          end if
       else
          if ( modulo( self%spline_degree, 2 ) == 0 ) then
             self%eval_same = .false.
             self%index_shift = (self%spline_degree)/2
          else
             self%eval_same = .true.
             self%index_shift = (self%spline_degree-1)/2
          end if
       end if
    end if

    ! normalized Gauss Legendre points and weights
    self%n_quad_points = (self%spline_degree+2)/2
    self%quad_xw = sll_f_gauss_legendre_points_and_weights(self%n_quad_points)
!!$
!!$    ! For nonlinear disgradE
!!$    self%n_quadp1_points = (self%spline_degree+3)/2
!!$    allocate( self%quadp1_xw(2, self%n_quadp1_points) )
!!$    allocate( self%spline_pol_val(self%n_span) )
!!$    ! normalized Gauss Legendre points and weights
!!$    self%quadp1_xw = sll_f_gauss_legendre_points_and_weights(self%n_quadp1_points)
!!$
!!$    
!!$    ! For smoothed evaluate and add_charge
!!$    allocate( self%quads1_xw(2, self%n_quadp1_points) )
!!$    ! normalized Gauss Legendre points and weights to [0,1]
!!$    self%quads1_xw = self%quadp1_xw
!!$    self%quads1_xw(1,:) = (self%quads1_xw(1,:) + 1.0_f64)*0.5_f64
!!$    self%quads1_xw(2,:) = self%quads1_xw(2,:)*0.5_f64
!!$
!!$    
!!$    
!!$    ! For smoothed add_current
!!$    self%n_quads2_points = (self%spline_degree+4)/2
!!$    allocate( self%quads2_xw(2, self%n_quads2_points) )
!!$    ! normalized Gauss Legendre points and weights to [0,1]
!!$    self%quads2_xw = sll_f_gauss_legendre_points_and_weights(self%n_quads2_points)
!!$    self%quads2_xw(1,:) = (self%quads2_xw(1,:) + 1.0_f64)*0.5_f64
!!$    self%quads2_xw(2,:) = self%quads2_xw(2,:)*0.5_f64

  end subroutine init_spline_strong_1d
  
  !> Destructor
  subroutine free_spline_strong_1d(self)
    class (sll_t_particle_mesh_coupling_spline_strong_1d), intent( inout ) :: self !< Particle mesh coupling object 

    deallocate(self%spline_val)
    deallocate(self%spline_val_more)
    deallocate(self%quad_xw)
    deallocate(self%spline_val_prim)
    deallocate(self%spline_val_prim_old)
    deallocate(self%spline_val_prim_new)
    if (self%spline_degree < 6 ) then
       call sll_s_spline_pp_free_1d(self%spline_pp)
    end if
    
  end subroutine free_spline_strong_1d

!!!!!!!!!!!! FROM HERE NOT IMPLEMENTED !!!!!!!!!!!!!


   !---------------------------------------------------------------------------!
  !> Add charge of one particle
  subroutine add_particle_mass_spline_strong_1d(self, position, marker_charge, particle_mass)
    class( sll_t_particle_mesh_coupling_spline_strong_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution


    SLL_ERROR("add_particle_mass", "Not implemented.")

  end subroutine add_particle_mass_spline_strong_1d

    !---------------------------------------------------------------------------!
  !> Evaluate several fields at position \a position
  subroutine evaluate_multiple_spline_strong_1d(self, position, components, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_strong_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_int32,                               intent(in)      :: components(:) !< Components of field_dofs that shall be updated
    sll_real64,                              intent( in )    :: field_dofs(:,:) !< Coefficient vector for the field DoFs
    sll_real64,                              intent(out)     :: field_value(:) !< Value(s) of the electric fields at given position

        SLL_ERROR("evaluate_multiple", "Not implemented.")

  end subroutine evaluate_multiple_spline_strong_1d


    subroutine evaluate_field_single_spline_strong_pp_1d(self, position, field_dofs_pp, field_value)
    class (sll_t_particle_mesh_coupling_spline_strong_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                              intent(in)      :: field_dofs_pp(:,:) !< Degrees of freedom in kernel representation.
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the electric fields at given position
    
    SLL_ERROR("evaluate_pp", "Not implemented.")

  end subroutine evaluate_field_single_spline_strong_pp_1d


   !> Add current with integration over x
  subroutine add_current_spline_strong_1d( self, position_old, position_new, marker_charge, j_dofs )
    class( sll_t_particle_mesh_coupling_spline_strong_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: j_dofs(self%n_dofs) !< Coefficient vector of the current density


    SLL_ERROR("add_current", "Not implemented.")
    
  end subroutine add_current_spline_strong_1d

    !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting)
  subroutine add_current_update_v_pp_spline_strong_1d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_strong_1d), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
    sll_real64, intent(in)    :: position_new(self%dim) !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(self%n_dofs) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(:) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(self%n_dofs) !< Coefficients of current expansion


    SLL_ERROR("add_current_update_v_pp", "Not implemented.")

  end subroutine add_current_update_v_pp_spline_strong_1d

    !----------------------------------------------------------------------------------
  !< Constructor for abstract type
  subroutine sll_s_new_particle_mesh_coupling_spline_strong_1d( smoother, domain, n_cells, spline_degree, integ, eval_grid_points )
    class( sll_c_particle_mesh_coupling_1d), allocatable, intent(out):: smoother !< kernel smoother object
    sll_int32,                                  intent(in) :: n_cells !< number of DoFs (spline coefficients)
    sll_real64,                                 intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32,                                  intent(in) :: spline_degree !< Degree of smoothing kernel spline
    logical,   intent( in ) :: integ
    logical,   intent( in ) :: eval_grid_points

    !local variables
    sll_int32 :: ierr


    allocate( sll_t_particle_mesh_coupling_spline_strong_1d :: smoother , stat=ierr)
    SLL_ASSERT( ierr == 0)
    
    select type( smoother )
    type is ( sll_t_particle_mesh_coupling_spline_strong_1d )
       call smoother%init( domain, n_cells, spline_degree, integ, eval_grid_points )
    end select

  end subroutine sll_s_new_particle_mesh_coupling_spline_strong_1d

end module sll_m_particle_mesh_coupling_spline_strong_1d
