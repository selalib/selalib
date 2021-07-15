!> @ingroup particle_mesh_coupling
!> @author Katharina Kormann, IPP
!> @brief Particle mesh coupling for 3d with splines of arbitrary degree placed on a uniform tensor product mesh.
!> @details Spline with index i starts at point i

module sll_m_particle_mesh_coupling_spline_3d_feec

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_gauss_legendre_integration, only : &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis, &
       sll_s_uniform_bsplines_eval_deriv

  use sll_m_particle_mesh_coupling_base_3d, only: &
       sll_c_particle_mesh_coupling_3d

  use sll_m_splines_pp, only : &
       sll_t_spline_pp_1d, &
       sll_t_spline_pp_3d, &
       sll_s_spline_pp_init_1d, &
       sll_s_spline_pp_init_3d, &
       sll_s_spline_pp_free_1d, &
       sll_s_spline_pp_free_3d, &
       sll_f_spline_pp_horner_1d, &
       sll_s_spline_pp_horner_m_1d, &
       sll_f_spline_pp_horner_3d, &
       sll_s_spline_pp_horner_m_3d
  
  implicit none

  public :: &
       sll_t_particle_mesh_coupling_spline_3d_feec

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Particle mesh coupling in 3d based on (arbitrary degree) spline on a tensor product uniform mesh
  type, extends(sll_c_particle_mesh_coupling_3d) :: sll_t_particle_mesh_coupling_spline_3d_feec

     sll_real64 :: rdelta_x(3)  !< Inverse values of delta_x
     ! Information about the particles
     sll_int32  :: no_particles !< Number of particles of underlying PIC method (processor local)
     sll_int32  :: n_quad_points !< Number of quadrature points

     sll_real64, allocatable :: spline_val(:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_0(:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_0_deriv(:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline1_0(:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_1(:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_1_deriv(:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_val_more(:,:) !< more scratch data for spline evaluation
     sll_real64, allocatable :: spline_2d(:,:), spline1_2d(:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline1_3d(:,:,:), spline2_3d(:,:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: j1d(:) !< scratch data for 1d slices of j
     sll_real64, allocatable :: quad_xw(:,:) !< quadrature weights and points

     sll_int32, allocatable :: index1d(:,:) !< 1d box index in every tensor direction

     type(sll_t_spline_pp_1d) :: spline_pp1d_0(3) !< 1d pp-splines in every tensor direction
     type(sll_t_spline_pp_1d) :: spline_pp1d_1(3) !< 1d pp-splines in every tensor direction

     ! For 1d quadrature in 3d space
     sll_int32 :: n_quad_points_line !< no. of quadrature points for the line integral
     sll_real64, allocatable :: quad_xw_line(:,:) !< quadrature weights and points for the line integral

   contains

     procedure :: add_charge => add_charge_single_spline_3d_feec !> Add charge of one particle
     procedure :: add_particle_mass => add_particle_mass_spline_3d_feec !> Collect diagonal entry of particle mass matrix
     procedure :: add_particle_mass_od => add_particle_mass_od_spline_3d_feec !> Collect off-diagonal entry of particle mass matrix
     procedure :: evaluate => evaluate_field_single_spline_3d_feec !> Evaluate spline function with given coefficients
     procedure :: evaluate_multiple => evaluate_multiple_spline_3d_feec !> Evaluate multiple spline functions with given coefficients
     procedure :: add_current => add_current_3d!> Add current via line integral
     procedure :: add_current_evaluate !>  Add current via line integral and evaluate spline
     procedure :: add_current_update_v_component1 => add_current_update_v_primitive_component1_spline_3d_feec !> Add contribution of one particle to the current density in x1 direction and update velocity
     procedure :: add_current_update_v_component2 => add_current_update_v_primitive_component2_spline_3d_feec !> Add contribution of one particle to the current density in x2 direction and update velocity
     procedure :: add_current_update_v_component3 => add_current_update_v_primitive_component3_spline_3d_feec !> Add contribution of one particle to the current density in x3 direction and update velocity
     procedure :: init => init_spline_3d_feec !> Constructor
     procedure :: free => free_spline_3d_feec !> Destructor
   
     procedure :: add_current_evaluate_int !> Evaluates the integral int_{poisition_old}^{position_new} field(x) d x and the integrated current
  end type sll_t_particle_mesh_coupling_spline_3d_feec

  type :: vector
     sll_real64, allocatable :: vals(:)
  end type vector

contains


  !---------------------------------------------------------------------------!
  !> Add charge of one particle
  subroutine add_charge_single_spline_3d_feec(self, position, marker_charge, degree, rho_dofs)
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(3) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_int32,                                intent( in )    :: degree(3) !< Spline degree along each dimension
    sll_real64,                               intent( inout ) :: rho_dofs(:) !< Coefficient vector of the charge distribution
    !local variables
    sll_int32 :: box(3)
    sll_real64 :: xi(3)
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k


    call convert_x_to_xbox( self, position, xi, box )

    self%spline_0 = 0.0_f64
    do j = 1, 3
       call sll_s_uniform_bsplines_eval_basis( degree(j), xi(j), self%spline_0(1:degree(j)+1,j) )
    end do

    ! Build scaling with marker_charge into spline along first dimension
    self%spline_0(:,1) = self%spline_0(:,1)*marker_charge
    ! 2d array combining first and second dimension
    do j = 1, degree(2)+1
       do i = 1, degree(1)+1
          self%spline_2d(i,j) = self%spline_0(i,1)*self%spline_0(j,2)
       end do
    end do

    box = box-degree
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 ) 
    do k = 1, degree(3)+1       
       index3d(3) = self%index1d(k,3)!modulo(box(3)+k-2,self%n_cells(3))+1
       do j = 1, degree(2)+1
          index3d(2) = self%index1d(j,2)!modulo(box(2)+j-2,self%n_cells(2))+1
          do i = 1, degree(1)+1
             index3d(1) = self%index1d(i,1)!modulo(box(1)+i-2,self%n_cells(1))+1
             index1d = convert_index_3d_to_1d( index3d, self%n_cells )
             rho_dofs(index1d) = rho_dofs(index1d) + &
                  (self%spline_2d(i,j) * self%spline_0(k,3))
             !(marker_charge * self%spline_0(i,1) * self%spline_0(j,2) * &
             !self%spline_0(k,3) )
          end do
       end do
    end do


  end subroutine add_charge_single_spline_3d_feec

  !---------------------------------------------------------------------------!
  !> Add charge of one particle
  subroutine add_charge_single_spline_pp_3d_feec(self, position, marker_charge, degree, rho_dofs)
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(3) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_int32,                                intent( in )    :: degree(3) !< Spline degree along each dimension
    sll_real64,                               intent( inout ) :: rho_dofs(self%n_total) !< Coefficient vector of the charge distribution
    !local variables
    sll_int32 :: box(3)
    sll_real64 :: xi(3)
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k


    call convert_x_to_xbox( self, position, xi, box )

    self%spline_0 = 0.0_f64
    call sll_s_spline_pp_horner_m_3d (self%spline_pp_0, self%spline_0, degree, xi)


    ! Build scaling with marker_charge into spline along first dimension
    self%spline_0(:,1) = self%spline_0(:,1)*marker_charge
    ! 2d array combining first and second dimension
    do j = 1, degree(2)+1
       do i = 1, degree(1)+1
          self%spline_2d(i,j) = self%spline_0(i,1)*self%spline_0(j,2)
       end do
    end do


    box = box-degree
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 ) 
    do k = 1, degree(3)+1       
       index3d(3) = self%index1d(k,3)!modulo(box(3)+k-2,self%n_cells(3))+1
       do j = 1, degree(2)+1
          index3d(2) = self%index1d(j,2)!modulo(box(2)+j-2,self%n_cells(2))+1
          do i = 1, degree(1)+1
             index3d(1) = self%index1d(i,1)!modulo(box(1)+i-2,self%n_cells(1))+1
             index1d = convert_index_3d_to_1d( index3d, self%n_cells )
             rho_dofs(index1d) = rho_dofs(index1d) + &
                  (self%spline_2d(i,j) * &
                  self%spline_0(k,3) )
          end do
       end do
    end do


  end subroutine add_charge_single_spline_pp_3d_feec

  !---------------------------------------------------------------------------!
  !> Add charge of one particle
  subroutine add_particle_mass_spline_3d_feec(self, position, marker_charge, degree, particle_mass )
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(3) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_int32,                                intent( in )    :: degree(3) !< Spline degree along each dimension
    sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution
    !local variables
    sll_int32 :: box(3)
    sll_real64 :: xi(3)
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k, col1, col2, col3, ind
    sll_real64 :: splineijk, splineijkcol3


    call convert_x_to_xbox( self, position, xi, box )

    ! Old version based on arbitrary degree splines
    self%spline_0 = 0.0_f64
    do j = 1, 3
       call sll_s_uniform_bsplines_eval_basis( degree(j), xi(j), self%spline_0(1:degree(j)+1,j) )
    end do

    ! 2d array combining first and second dimension
    do j = 1, degree(2)+1
       do i = 1, degree(1)+1
          self%spline_2d(i,j) = self%spline_0(i,1)*self%spline_0(j,2)
       end do
    end do
    ! TODO: Check if also 3d array would make sense


    box = box-degree
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 ) 
    do k = 1, degree(3)+1       
       index3d(3) = self%index1d(k,3)!modulo(box(3)+k-2,self%n_cells(3))+1
       do j = 1, degree(2)+1
          index3d(2) = self%index1d(j,2)!modulo(box(2)+j-2,self%n_cells(2))+1
          do i = 1, degree(1)+1
             index3d(1) = self%index1d(i,1)!modulo(box(1)+i-2,self%n_cells(1))+1
             index1d = convert_index_3d_to_1d( index3d, self%n_cells )
             ind = 1+(degree(1)+1-i)+(degree(2)+1-j)*(2*degree(1)+1)+ &
                  (degree(3)+1-k)*(2*degree(1)+1)*(2*degree(2)+1)
             splineijk = marker_charge * self%spline_2d(i,j) *&
                  self%spline_0( k, 3)
             do col3 = 1,degree(3)+1
                splineijkcol3 = splineijk*self%spline_0(col3, 3)
                do col2 = 1,degree(2)+1
                   do col1 = 1,degree(1)+1                     
                      particle_mass(ind, index1d) = &
                           particle_mass( ind, index1d) + &
                           splineijkcol3 * &
                           self%spline_2d(col1,col2)
                      ind = ind+1
                   end do
                   ind = ind+degree(1)
                end do
                ind = ind+ degree(2) * (2*degree(1)+1)
             end do
          end do
       end do
    end do


  end subroutine add_particle_mass_spline_3d_feec


  !> Add charge of one particle
  subroutine add_particle_mass_od_spline_3d_feec(self, position, marker_charge, degree1, degree2, particle_mass )
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(3) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_int32,                                intent( in )    :: degree1(3), degree2(3) !< Spline degree along each dimension
    sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution
    !local variables
    sll_int32 :: box(3)
    sll_real64 :: xi(3)
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k, col1, col2, col3, ind
    sll_real64 :: splineijk, splineijkcol3


    call convert_x_to_xbox( self, position, xi, box )

    self%spline_0 = 0.0_f64
    self%spline1_0 = 0.0_f64
    do j = 1, 3
       call sll_s_uniform_bsplines_eval_basis( degree1(j), xi(j), self%spline_0(1:degree1(j)+1,j) )
       call sll_s_uniform_bsplines_eval_basis( degree2(j), xi(j), self%spline1_0(1:degree2(j)+1,j) )
    end do

    self%spline1_3d = 0._f64
    !3d array
    do k = 1, degree1(3)+1
       do j = 1, degree1(2)+1
          do i = 1, degree1(1)+1
             self%spline1_3d(i,j,k) = self%spline_0(i,1)*self%spline_0(j,2)*self%spline_0(k,3)
          end do
       end do
    end do

    self%spline2_3d = 0._f64
    do k = 1, degree2(3)+1
       do j = 1, degree2(2)+1
          do i = 1, degree2(1)+1
             self%spline2_3d(i,j,k) = self%spline1_0(i,1)*self%spline1_0(j,2)*self%spline1_0(k,3)
          end do
       end do
    end do


    box = box-degree1
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 ) 
    do k = 1, degree1(3)+1
       index3d(3) = self%index1d(k,3)!modulo(box(3)+k-2,self%n_cells(3))+1
       do j = 1, degree1(2)+1
          index3d(2) = self%index1d(j,2)!modulo(box(2)+j-2,self%n_cells(2))+1
          do i = 1, degree1(1)+1
             index3d(1) = self%index1d(i,1)!modulo(box(1)+i-2,self%n_cells(1))+1
             index1d = convert_index_3d_to_1d( index3d, self%n_cells )
             ind = 1+(degree1(1)+1-i)+(degree1(2)+1-j)*(degree1(1)+degree2(1)+1)+ &
                  (degree1(3)+1-k)*(degree1(1)+degree2(1)+1)*(degree1(2)+degree2(2)+1)
             do col3 = 1,degree2(3)+1
                do col2 = 1,degree2(2)+1
                   do col1 = 1,degree2(1)+1
                      particle_mass( ind, index1d)=particle_mass( ind, index1d)+&
                           self%spline1_3d(i,j,k)*self%spline2_3d(col1,col2,col3)*&
                           marker_charge
                      ind = ind+1
                   end do
                   ind = ind+degree1(1)
                end do
                ind = ind+degree1(2) * (degree1(1)+degree2(1)+1)
             end do
          end do
       end do
    end do

  end subroutine add_particle_mass_od_spline_3d_feec

  
  !---------------------------------------------------------------------------!
  !> Evaluate field at at position \a position
  subroutine evaluate_field_single_spline_3d_feec(self, position, degree, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_3d_feec), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(3) !< Position of the particle
    sll_int32 ,                              intent( in )    :: degree(3) !< Spline degree of the various components
    sll_real64,                              intent( in )    :: field_dofs(:) !< Coefficient vector for the field DoFs
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the fields at given position
    !local variables
    sll_int32 :: i,j,k
    sll_int32 :: box(3)
    sll_int32 :: index3d(3)
    sll_int32 :: index1d
    sll_real64 :: xi(3)

    ! TODO: Optimize by sum factorization

    call convert_x_to_xbox( self, position, xi, box )

    do j = 1, 3
       call sll_s_uniform_bsplines_eval_basis( degree(j), xi(j), self%spline_val(1:degree(j)+1,j) )
    end do

    field_value = 0.0_f64
    box = box-degree
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 )  
    do k = 1, degree(3)+1       
       index3d(3) = self%index1d(k,3)!modulo(box(3)+k-2,self%n_cells(3))+1
       do j = 1, degree(2)+1
          index3d(2) = self%index1d(j,2)!modulo(box(2)+j-2,self%n_cells(2))+1
          do i = 1, degree(1)+1
             index3d(1) = self%index1d(i,1)!modulo(box(1)+i-2,self%n_cells(1))+1
             index1d = convert_index_3d_to_1d( index3d, self%n_cells )
             field_value = field_value + &
                  field_dofs(index1d) * &
                  self%spline_val(i,1) * self%spline_val(j,2) * &
                  self%spline_val(k,3)
          end do
       end do
    end do




  end subroutine evaluate_field_single_spline_3d_feec

  !---------------------------------------------------------------------------!
  !> Evaluate field at at position \a position
  subroutine evaluate_field_single_spline_pp_3d_feec(self, position, degree, field_dofs_pp, field_value)
    class (sll_t_particle_mesh_coupling_spline_3d_feec), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(3) !< Position of the particle
    sll_int32 ,                              intent( in )    :: degree(3) !< Spline degree of the various components
    sll_real64,                              intent( in )    :: field_dofs_pp(:,:) !< Coefficient vector for the field DoFs
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the fields at given position

    sll_int32 :: box(3)
    sll_real64 :: xi(3)

    call convert_x_to_xbox( self, position, xi, box )

    field_value = sll_f_spline_pp_horner_3d(degree, field_dofs_pp, xi, box,self%n_cells)

  end subroutine evaluate_field_single_spline_pp_3d_feec



  !---------------------------------------------------------------------------!
  !> Evaluate several fields at position \a position
  subroutine evaluate_multiple_spline_3d_feec(self, position, components, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_3d_feec), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(3) !< Position of the particle
    sll_int32,                               intent(in)      :: components(:) !< Components of field_dofs that shall be updated
    sll_real64,                              intent( in )    :: field_dofs(:,:) !< Coefficient vector for the field DoFs
    sll_real64,                              intent(out)     :: field_value(:) !< Value(s) of the electric fields at given position



  end subroutine evaluate_multiple_spline_3d_feec

  !> Convert 3d index to 1d index (first dimension without stride)
  pure function convert_index_3d_to_1d( index3d, n_cells ) result( index1d )
    sll_int32, intent( in ) :: index3d(3)
    sll_int32, intent( in ) :: n_cells(3)
    sll_int32 :: index1d

    index1d = index3d(1) + (index3d(2)-1)*n_cells(1) + (index3d(3)-1)*n_cells(1)*n_cells(2)

  end function convert_index_3d_to_1d

  !> Identify the box in which the particle is located and its normalized position within the box
  subroutine convert_x_to_xbox( self, position, xi, box )
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(3) !< Position of the particle
    sll_real64,                               intent( out )    :: xi(3) !< Position of the particle
    sll_int32,                                intent( out )    :: box(3) !< Position of the particle

    xi = (position - self%domain(:,1)) /self%delta_x!* self%rdelta_x ! destroys perfect gauss conservation for symplectic splitting with fft solver
    box = floor( xi ) + 1
    xi = xi - real(box-1, f64)

  end subroutine convert_x_to_xbox

  !> Identify the box in which the particle is located and its normalized position within the box (only along the dimension \a component)
  subroutine convert_x_to_xbox_1d( self, component, position, xi, box )
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_int32,                                intent( in )    :: component !< direction along which the position is given
    sll_real64,                               intent( in )    :: position !< Position of the particle
    sll_real64,                               intent( out )    :: xi !< Position of the particle
    sll_int32,                                intent( out )    :: box !< Position of the particle

    xi = (position - self%domain(component,1))/self%delta_x(component)!* self%rdelta_x ! destroys perfect gauss conservation for symplectic splitting with fft solver
    box = floor( xi ) + 1
    xi = xi - real(box-1, f64)

  end subroutine convert_x_to_xbox_1d


  !> Sets the index of the splines that are involved in computations that concern splines from index \a box
  subroutine box_index( self, box, comp )
    class(sll_t_particle_mesh_coupling_spline_3d_feec), intent(inout) :: self !< kernel smoother object
    sll_int32, intent(in) :: box !> starting box number
    sll_int32, intent(in) :: comp !> direction along which to operate
    !local variables
    sll_int32 :: j

    do j = 1, self%spline_degree(comp)+1
       self%index1d(j,comp) = modulo(box+j-2,self%n_cells(comp))+1
    end do

  end subroutine box_index

  !>  Add current via line integral and evaluate spline
  subroutine add_current_evaluate ( self, position_old, position_new, xdot, efield_dofs, j_dofs, efield_val)
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(3) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(3) !< Position of the particle
    sll_real64,                               intent( in )    :: xdot(3) !< velocity
    sll_real64,                               intent( in )    :: efield_dofs(:) !< electric field dofs
    sll_real64,                               intent( inout ) :: j_dofs(:) !< current dofs
    sll_real64,                               intent(   out ) :: efield_val(3) !< electric field value
    !local variables
    sll_real64 :: xold(3), xnew(3), xnewtilde(3), xbox(3)
    sll_int32 :: boxold(3), boxnew(3), sigma_counter(3), boxdiff(3), increment(3), box(3)
    sll_real64 :: sigma_r, sigma_l, sigma, sigma_next(3), field_value(3), weight
    sll_int32 :: j, q, bl, maxind, index_is
    type(vector) :: sigmas(3)

    call convert_x_to_xbox( self, position_old, xold, boxold )
    call convert_x_to_xbox( self, position_new, xnew, boxnew )

    ! Now we need to compute the normalized 1d line along which to integrate:
    boxdiff = boxnew-boxold
    xnewtilde = xnew + real(boxdiff,f64)

    sigma_l = 0.0_f64
    box = boxold ! We start in the old box

    sigma_counter = 0

    efield_val = 0.0_f64

    do j = 1, 3
       if (boxdiff(j) > 0 ) then
          allocate ( sigmas(j)%vals(boxdiff(j)+1) )
          do bl = 1, boxdiff(j)
             sigmas(j)%vals(bl) = (real(bl,f64)  - xold(j))/(xnewtilde(j)-xold(j))
          end do
          sigmas(j)%vals(boxdiff(j)+1) = 1.0_f64
          sigma_next(j) = sigmas(j)%vals(1)
          increment(j) = 1
       elseif (boxdiff(j) < 0 ) then
          allocate ( sigmas(j)%vals(-boxdiff(j)+1) )
          do bl = boxdiff(j)+1, 0
             sigmas(j)%vals(-bl+1) = (real(bl,f64)  - xold(j))/(xnewtilde(j)-xold(j))
          end do
          sigmas(j)%vals(-boxdiff(j)+1) = 1.0_f64
          sigma_next(j) = sigmas(j)%vals(1)
          increment(j) = -1
       else
          sigma_next(j) = 1.0_f64
          increment(j) = 0
       end if
    end do

    sigma_r = 0.0_f64
    do while ( sigma_r < 1.0_f64 )
       ! Identify index of next intersection
       index_is = minloc(sigma_next, dim=1)
       sigma_r = sigma_next(index_is)

       do q = 1, self%n_quad_points_line
          sigma = sigma_l + (sigma_r-sigma_l) * self%quad_xw_line(1,q)
          xbox = xold* (1.0_f64 - sigma) + xnewtilde * sigma  - real(sigma_counter*increment, f64)
          if (maxval(xbox) > 1.0_f64 ) then
             print*, xbox, sigma, sigma_counter, increment
             SLL_ERROR( 'add_current_evaluate', 'box value too large')
          elseif (minval(xbox) < 0.0_f64 ) then
             print*, xbox, sigma, sigma_counter, increment
             print*, xold, xnewtilde, sigma_r
             SLL_ERROR( 'add_current_evaluate', 'box value too low')
          end if

          weight = self%quad_xw_line(2,q)*(sigma_r-sigma_l)

          call point_add_eval (self, box, xbox, efield_dofs, xdot*weight, j_dofs, field_value )
          efield_val = efield_val + field_value * weight
       end do
       if (sigma_r < 1.0_f64 ) then
          ! Update the
          sigma_counter(index_is) = sigma_counter(index_is)+1
          sigma_next(index_is) = sigmas(index_is)%vals(sigma_counter(index_is)+1)
          box(index_is) = box(index_is) + increment(index_is)
          sigma_l = sigma_r
       end if
    end do

  end subroutine add_current_evaluate

  !> Helper function for add_current_evaluate that takes care of the per-cell computations
  subroutine point_add_eval ( self, box_in, xbox, field_dofs, weight, j_dofs, field_value )
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_int32, intent(in) :: box_in(3)
    sll_real64, intent(in) :: xbox(3)
    sll_real64, intent(in) :: field_dofs(:)
    sll_real64, intent(in) :: weight(3)
    sll_real64, intent(inout) :: j_dofs(:)
    sll_real64, intent(out) :: field_value(3)
    !local variables
    sll_int32 :: i,j,k, index1d, n_total
    sll_int32 :: box(3), index3d(3)
    sll_real64 :: spval

    n_total = self%n_total
    field_value = 0.0_f64

    do j = 1, 3
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree(j), xbox(j), &
            self%spline_0(1:self%spline_degree(j)+1,j) )
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree(j)-1, xbox(j), &
            self%spline_1(1:self%spline_degree(j),j) )
    end do

    box = box_in-self%spline_degree

    box(1) = box(1)+1
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 )      
    do k = 1, self%spline_degree(3)+1 
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j = 1, self%spline_degree(2)+1
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i = 1, self%spline_degree(1)
             index1d = index3d(2) + self%index1d(i,1)

             spval =  self%spline_1(i,1) * &
                  self%spline_0(j,2) * self%spline_0(k,3)
             field_value(1) = field_value(1) + field_dofs(index1d) * spval
             j_dofs(index1d) = j_dofs(index1d) + &
                  weight(1) * spval
          end do
       end do
    end do

    box(1) = box(1)-1
    box(2) = box(2)+1
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    do k = 1, self%spline_degree(3)+1       
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j = 1, self%spline_degree(2)
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i = 1, self%spline_degree(1)+1
             index1d = index3d(2) + self%index1d(i,1) + n_total

             spval = self%spline_0(i,1) * &
                  self%spline_1(j,2) * self%spline_0(k,3)
             field_value(2) = field_value(2) + field_dofs(index1d ) * spval
             j_dofs(index1d) = j_dofs(index1d) + &
                  weight(2) * spval
          end do
       end do
    end do

    n_total = n_total*2

    box(2) = box(2)-1
    box(3) = box(3)+1
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 ) 
    do k = 1, self%spline_degree(3)    
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j = 1, self%spline_degree(2)+1
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i = 1, self%spline_degree(1)+1
             index1d = index3d(2) + self%index1d(i,1) + n_total
             spval = self%spline_0(i,1) * &
                  self%spline_0(j,2) * self%spline_1(k,3)
             field_value(3) = field_value(3) + field_dofs(index1d ) * spval
             j_dofs(index1d) = j_dofs(index1d) + &
                  weight(3) * spval
          end do
       end do
    end do


  end subroutine point_add_eval

  !> Add current via line integral (when x changes along all three directions)
  subroutine add_current_3d  ( self, position_old, position_new, xdot, j_dofs) 
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in    ) :: position_old(3) !< Position of the particle
    sll_real64,                               intent( in    ) :: position_new(3) !< Position of the particle
    sll_real64,                               intent( in    ) :: xdot(3) !< velocity
    sll_real64,                               intent( inout ) :: j_dofs(:) !< current dofs
    !local variables
    sll_real64 :: xold(3), xnew(3), xnewtilde(3), xbox(3)
    sll_int32 :: boxold(3), boxnew(3), sigma_counter(3), boxdiff(3), increment(3), box(3)
    sll_real64 :: sigma_r, sigma_l, sigma, sigma_next(3), field_value(3), weight
    sll_int32 :: j, q, bl, maxind, index_is
    type(vector) :: sigmas(3)


    call convert_x_to_xbox( self, position_old, xold, boxold )
    call convert_x_to_xbox( self, position_new, xnew, boxnew )

    ! Now we need to compute the normalized 1d line along which to integrate:
    boxdiff = boxnew-boxold
    xnewtilde = xnew + real(boxdiff,f64)

    sigma_l = 0.0_f64
    box = boxold ! We start in the old box

    sigma_counter = 0

    do j = 1, 3
       if (boxdiff(j) > 0 ) then
          allocate ( sigmas(j)%vals(boxdiff(j)+1) )
          do bl = 1, boxdiff(j)
             sigmas(j)%vals(bl) = (real(bl,f64)  - xold(j))/(xnewtilde(j)-xold(j))
          end do
          sigmas(j)%vals(boxdiff(j)+1) = 1.0_f64
          sigma_next(j) = sigmas(j)%vals(1)
          increment(j) = 1
       elseif (boxdiff(j) < 0 ) then
          allocate ( sigmas(j)%vals(-boxdiff(j)+1) )
          do bl = boxdiff(j)+1, 0
             sigmas(j)%vals(-bl+1) = (real(bl,f64)  - xold(j))/(xnewtilde(j)-xold(j))
          end do
          sigmas(j)%vals(-boxdiff(j)+1) = 1.0_f64
          sigma_next(j) = sigmas(j)%vals(1)
          increment(j) = -1
       else
          sigma_next(j) = 1.0_f64
          increment(j) = 0
       end if
    end do

    sigma_r = 0.0_f64
    do while ( sigma_r < 1.0_f64 )
       ! Identify index of next intersection
       index_is = minloc(sigma_next, dim=1)
       sigma_r = sigma_next(index_is)

       do q = 1, self%n_quad_points_line
          sigma = sigma_l + (sigma_r-sigma_l) * self%quad_xw_line(1,q)
          xbox = xold* (1.0_f64 - sigma) + xnewtilde * sigma  - real(sigma_counter*increment, f64)
          if (maxval(xbox) > 1.0_f64 ) then
             print*, xbox, sigma, sigma_counter, increment
             SLL_ERROR( 'add_current_3d', 'box value too large')
          elseif (minval(xbox) < 0.0_f64 ) then
             print*, xbox, sigma, sigma_counter, increment
             SLL_ERROR( 'add_current_3d', 'box value too low')
          end if

          weight = self%quad_xw_line(2,q)*(sigma_r-sigma_l)

          call integrate_spline_3d(self, box, xbox, xdot*weight, j_dofs )
       end do
       if (sigma_r < 1.0_f64 ) then
          ! Update the
          sigma_counter(index_is) = sigma_counter(index_is)+1
          sigma_next(index_is) = sigmas(index_is)%vals(sigma_counter(index_is)+1)
          box(index_is) = box(index_is) + increment(index_is)
          sigma_l = sigma_r
       end if
    end do
  end subroutine add_current_3d

  !> Helper function for add_current_3d, takes care of the per box computations
  subroutine integrate_spline_3d( self, box_in, xbox, weight, j_dofs )
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_int32,  intent( in    ) :: box_in(3)
    sll_real64, intent( in    ) :: xbox(3)
    sll_real64, intent( in    ) :: weight(3)
    sll_real64, intent( inout ) :: j_dofs(:)
    !local variables
    sll_int32 :: i,j,k, index1d
    sll_int32 :: box(3), index3d(3)

    do j = 1, 3
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree(j), xbox(j), &
            self%spline_0(1:self%spline_degree(j)+1,j) )
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree(j)-1, xbox(j), &
            self%spline_1(1:self%spline_degree(j),j) )
    end do

    box = box_in-self%spline_degree

    box(1) = box(1)+1
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 ) 
    do k = 1, self%spline_degree(3)+1 
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j = 1, self%spline_degree(2)+1
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i = 1, self%spline_degree(1)
             index1d = index3d(2) + self%index1d(i,1)

             j_dofs(index1d) = j_dofs(index1d) + &
                  weight(1) * self%spline_1(i,1) * &
                  self%spline_0(j,2) * self%spline_0(k,3)
          end do
       end do
    end do

    box(1) = box(1)-1
    box(2) = box(2)+1
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    do k = 1, self%spline_degree(3)+1       
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j = 1, self%spline_degree(2)
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i = 1, self%spline_degree(1)+1
             index1d = index3d(2) + self%index1d(i,1) + self%n_total

             j_dofs(index1d) = j_dofs(index1d) + &
                  weight(2) * self%spline_0(i,1) * &
                  self%spline_1(j,2) * self%spline_0(k,3)
          end do
       end do
    end do

    box(2) = box(2)-1
    box(3) = box(3)+1
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 ) 
    do k = 1, self%spline_degree(3)    
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j = 1, self%spline_degree(2)+1
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i = 1, self%spline_degree(1)+1
             index1d = index3d(2) + self%index1d(i,1) + 2*self%n_total

             j_dofs(index1d) = j_dofs(index1d) + &
                  weight(3) * self%spline_0(i,1) * &
                  self%spline_0(j,2) * self%spline_1(k,3)
          end do
       end do
    end do


  end subroutine integrate_spline_3d

  !> Add current for one particle 
  subroutine add_current_spline_3d_feec (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_3d_feec), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(3) !< Position at time t
    sll_real64, intent(in)    :: position_new(3) !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(self%n_total*3) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(:) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(self%n_total) !< Coefficients of current expansion


    sll_int32 :: box(3), boxnew(3), boxold(3), local_size(3)
    sll_int32  :: degree
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k
    sll_real64 :: xnew(3), xold(3)
    sll_int32  :: component

    component = 1
    degree = self%spline_degree(component)
    call convert_x_to_xbox( self, position_old, xold, boxold )
    call convert_x_to_xbox( self, position_new, xnew, boxnew )

    local_size = abs(boxnew-boxold)+degree
    local_size(2) = self%spline_degree(2)+1
    local_size(3) = self%spline_degree(3)+1

    SLL_ASSERT_ALWAYS( local_size(1) <= self%n_cells(1) )

    do i=1, degree
       self%spline_1(i,1) =  sll_f_spline_pp_horner_1d(degree, self%spline_pp1d_1(component)%poly_coeffs_fpa, xold(1), i) &
            * self%delta_x(1)
       self%spline_1(i,2) =  sll_f_spline_pp_horner_1d(degree, self%spline_pp1d_1(component)%poly_coeffs_fpa, xnew(1), i) &
            * self%delta_x(1)
    end do

    if (position_old(1) < position_new(1) ) then
       self%j1d(local_size(component)-degree+1:local_size(component)) = self%spline_1(1:degree,2)
       self%j1d(1:local_size(component)-degree) = self%delta_x(1)
       self%j1d(1:degree) = self%j1d(1:degree) - self%spline_1(1:degree,1)
    else
       self%j1d(1:local_size(component)-degree) = -self%delta_x(1)
       self%j1d(local_size(component)-degree+1:local_size(component)) = -self%spline_1(1:degree,1)
       self%j1d(1:degree) = self%j1d(1:degree) + self%spline_1(1:degree,2)
    end if

    self%spline_0 = 0.0_f64
    do j=1,3
       if (j .ne. component ) then
          call sll_s_uniform_bsplines_eval_basis( self%spline_degree(j), xold(j), self%spline_0(1:self%spline_degree(j)+1,j) )
       end if
    end do

    box = boxold-self%spline_degree
    box(component) = min(boxnew(component), boxold(component))-degree+1
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 ) 
    do k=1,local_size(3)
       index3d(3) = self%index1d(k,3)!modulo(box(3)+k-2,self%n_cells(3))+1
       do j=1,local_size(2)
          index3d(2) = self%index1d(j,2)!modulo(box(2)+j-2,self%n_cells(2))+1
          do i=1,local_size(1)
             index3d(1) = self%index1d(i,1)!modulo(box(1)+i-2,self%n_cells(1))+1
             index1d = convert_index_3d_to_1d( index3d, self%n_cells )
             !print*, index1d, self%j1d(i), self%spline_0(j,2), self%spline_0(k,3)
             j_dofs(index1d) = j_dofs(index1d) + &
                  marker_charge * self%j1d( i ) * &
                  self%spline_0(j,2) * self%spline_0(k,3)
          end do
       end do
    end do



  end subroutine add_current_spline_3d_feec


  !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting), version based on  primitive function
  subroutine add_current_update_v_primitive_component1_spline_3d_feec (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_3d_feec), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(3) !< Position at time t
    sll_real64, intent(in)    :: position_new !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(:) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(3) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(:) !< Coefficients of current expansion
    !local variables
    sll_int32 :: box(3), boxnew, boxold
    sll_real64 :: xi(3)
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k
    sll_real64 :: xnew
    sll_real64 :: vt(2), vtt2, vtt3
    sll_int32  :: start1, start2
    sll_int32  :: component
    sll_int32  :: startjk, stride
    sll_real64 :: splinejk
    sll_int32 :: local_size
    sll_int32  :: degree

    component = 1
    start1 = 2*self%n_total
    start2 = self%n_total
    stride = 1
    degree = self%spline_degree(component)

    call convert_x_to_xbox( self, position_old, xi, box )
    call convert_x_to_xbox_1d( self, component, position_new, xnew, boxnew )
    boxold = box(component)

    !-- For current along x1
    local_size = abs(boxnew-boxold)+degree

    SLL_ASSERT_ALWAYS( local_size <= self%n_cells(1)+degree )

    ! For j=component, we need the primitive
    do i=1, degree
       self%spline_1(i,1) =  sll_f_spline_pp_horner_1d(degree, &
            self%spline_pp1d_1(component)%poly_coeffs_fpa, xi(component), i) &
            * self%delta_x(component)
       self%spline_1(i,2) =  sll_f_spline_pp_horner_1d(degree, &
            self%spline_pp1d_1(component)%poly_coeffs_fpa, xnew, i) &
            * self%delta_x(component)
    end do


    if (position_old(component) .le. position_new ) then
       self%j1d(local_size-degree+1:local_size) = self%spline_1(1:degree,2)
       self%j1d(1:local_size-degree) = self%delta_x(component)
       self%j1d(1:degree) = self%j1d(1:degree) - self%spline_1(1:degree,1)
    else
       self%j1d(1:local_size-degree) = -self%delta_x(component)
       self%j1d(local_size-degree+1:local_size) = -self%spline_1(1:degree,1)
       self%j1d(1:degree) = self%j1d(1:degree) + self%spline_1(1:degree,2)
    end if
    !----


    ! Achtung wir brauchen nun splines von beidem Grad
    self%spline_0 = 0.0_f64
    self%spline_1 = 0.0_f64
    do j=1,3
       if (j .ne. component ) then
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_0(j), self%spline_0(1:self%spline_degree(j)+1,j), self%spline_degree(j), xi(j))
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_1(j), self%spline_1(1:self%spline_degree(j),j), self%spline_degree(j)-1, xi(j))
       end if
    end do

    !box(component) = box(component)-degree+1
    box(2:3) = box(2:3) - self%spline_degree(2:3)

    ! Define the range of the first component
    if (boxold<boxnew) then
       box(component) = boxold-degree+1
    else
       box(component) = boxnew-degree+1
    end if

    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 )   
    do k = 1, self%spline_degree(3)+1      
       vtt2 = 0.0_f64
       vtt3 = 0.0_f64
       index3d(3) = self%index1d(k,3)!modulo(box(3)+k-2,self%n_cells(3))+1
       do j = 1, self%spline_degree(2)+1 
          index3d(2) = self%index1d(j,2)!modulo(box(2)+j-2,self%n_cells(2))+1

          startjk  = 1 + (index3d(2)-1)*self%n_cells(1) + &
               (index3d(3)-1)*self%n_cells(1)*self%n_cells(2)
          splinejk = self%spline_0(j, 2) * self%spline_0(k,3) * marker_charge

          vt = 0.0_f64
          do i = 1, local_size
             index3d(1) = modulo(box(1)+i-2,self%n_cells(1))
             index1d = startjk+index3d(1)!convert_index_3d_to_1d( index3d, self%n_cells )
             j_dofs(index1d) = j_dofs(index1d) + self%j1d(i) * &
                  splinejk 

             vt(1) = vt(1) + bfield_dofs(start1+index1d)*self%j1d(i)
             vt(2) = vt(2) + bfield_dofs(start2+index1d)*self%j1d(i)

          end do

          if (j>1) then
             vtt2 = vtt2 + vt(1)*self%spline_1(j-1, 2)
          end if
          vtt3 = vtt3 - vt(2)*self%spline_0(j, 2)

       end do
       vi(2) = vi(2) - qoverm*vtt2*self%spline_0(k, 3)
       if ( k> 1) then
          vi(3) = vi(3) - qoverm*vtt3*self%spline_1(k-1, 3)
       end if
    end do



  end subroutine add_current_update_v_primitive_component1_spline_3d_feec

  !> Add current for one particle and update v (according to H_p2 part in Hamiltonian splitting), version based on primitive function
  subroutine add_current_update_v_primitive_component2_spline_3d_feec (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_3d_feec), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(3) !< Position at time t
    sll_real64, intent(in)    :: position_new !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(:) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(3) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(:) !< Coefficients of current expansion
    !local variables
    sll_int32 :: box(3), boxnew, boxold
    sll_real64 :: xi(3)
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k
    sll_real64 :: xnew
    sll_real64 :: vt(2), vtt2, vtt3
    sll_int32  :: start1, start2
    sll_int32  :: component, local_size
    sll_int32 :: stride, startjk
    sll_real64 :: splineik
    sll_int32  :: degree

    component = 2
    start1 = 2*self%n_total 
    start2 = 0
    stride = self%n_cells(1)
    degree = self%spline_degree(component)

    call convert_x_to_xbox( self, position_old, xi, box )
    call convert_x_to_xbox_1d( self, component, position_new, xnew, boxnew )
    boxold = box(component)

    !-- For current along x2
    local_size = abs(boxnew-boxold)+degree

    SLL_ASSERT_ALWAYS( local_size <= self%n_cells(2)+degree )

    ! For j=component, we need the primitive
    do i=1, degree
       self%spline_1(i,1) =  sll_f_spline_pp_horner_1d(degree, &
            self%spline_pp1d_1(component)%poly_coeffs_fpa, xi(component), i) &
            * self%delta_x(component)
       self%spline_1(i,2) =  sll_f_spline_pp_horner_1d(degree, &
            self%spline_pp1d_1(component)%poly_coeffs_fpa, xnew, i) &
            * self%delta_x(component)
    end do

    if (position_old(component) .le. position_new ) then
       self%j1d(local_size-degree+1:local_size) = self%spline_1(1:degree,2)
       self%j1d(1:local_size-degree) = self%delta_x(component)
       self%j1d(1:degree) = self%j1d(1:degree) - self%spline_1(1:degree,1)
    else
       self%j1d(1:local_size-degree) = -self%delta_x(component)
       self%j1d(local_size-degree+1:local_size) = -self%spline_1(1:degree,1)
       self%j1d(1:degree) = self%j1d(1:degree) + self%spline_1(1:degree,2)
    end if
    !----


    ! Achtung wir brauchen nun splines von beidem Grad
    self%spline_0 = 0.0_f64
    self%spline_1 = 0.0_f64
    do j=1,3
       if (j .ne. component ) then
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_0(j), self%spline_0(1:self%spline_degree(j)+1,j), self%spline_degree(j), xi(j))
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_1(j), self%spline_1(1:self%spline_degree(j),j), self%spline_degree(j)-1, xi(j))
       end if
    end do

    box(1) = box(1) - self%spline_degree(1)
    box(3) = box(3) - self%spline_degree(3)

    ! Define the range of the first component
    if (boxold<boxnew) then
       box(component) = boxold-degree+1
    else
       box(component) = boxnew-degree+1
    end if

    call box_index( self, box(1), 1 )
    call box_index( self, box(3), 3 )
    do k = 1, self%spline_degree(3)+1      
       vtt2 = 0.0_f64
       vtt3 = 0.0_f64
       index3d(3) = self%index1d(k,3)!modulo(box(3)+k-2,self%n_cells(3))+1
       do j = 1, self%spline_degree(1)+1 
          index3d(1) = self%index1d(j,1)!modulo(box(1)+j-2,self%n_cells(1))+1
          
          startjk  = index3d(1) + &
               (index3d(3)-1)*self%n_cells(1)*self%n_cells(2)
          splineik = self%spline_0(j, 1) * self%spline_0(k,3)  * marker_charge

          vt = 0.0_f64
          do i = 1, local_size
             index3d(2) = modulo(box(2)+i-2,self%n_cells(2))
             index1d = startjk+index3d(2)*stride
             j_dofs(index1d) = j_dofs(index1d) + self%j1d(i) * &
                  splineik

             vt(1) = vt(1) + bfield_dofs(start1+index1d)*self%j1d(i)
             vt(2) = vt(2) + bfield_dofs(start2+index1d)*self%j1d(i)
          end do

          if (j>1) then
             vtt2 = vtt2 + vt(1)*self%spline_1(j-1, 1)
          end if
          vtt3 = vtt3 + vt(2)*self%spline_0(j, 1)

       end do
       vi(1) = vi(1) + qoverm*vtt2*self%spline_0(k, 3)
       if ( k> 1) then
          vi(3) = vi(3) - qoverm*vtt3*self%spline_1(k-1, 3)
       end if
    end do


  end subroutine add_current_update_v_primitive_component2_spline_3d_feec


  
  !> Add current for one particle and update v (according to H_p3 part in Hamiltonian splitting), version based on primitive function
  subroutine add_current_update_v_primitive_component3_spline_3d_feec (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_3d_feec), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(3) !< Position at time t
    sll_real64, intent(in)    :: position_new !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(:) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(3) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(:) !< Coefficients of current expansion
    !local variables
    sll_int32 :: box(3), boxnew, boxold
    sll_real64 :: xi(3)
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k
    sll_real64 :: xnew
    sll_real64 :: vt(2), vtt2, vtt3
    sll_int32  :: start1, start2
    sll_int32  :: component
    sll_int32  :: stride, startjk
    sll_real64 :: splineij
    sll_int32  :: local_size
    sll_int32  :: degree

    component = 3
    start1 = self%n_total
    start2 = 0
    stride = self%n_cells(1)*self%n_cells(2)
    degree = self%spline_degree(component)

    call convert_x_to_xbox( self, position_old, xi, box )
    call convert_x_to_xbox_1d( self, component, position_new, xnew, boxnew )
    boxold = box(component)

    !-- For current along x3
    local_size = abs(boxnew-boxold)+degree

    SLL_ASSERT_ALWAYS( local_size <= self%n_cells(3)+degree )
    
    ! For j=component, we need the primitive
    do i=1, degree
       self%spline_1(i,1) =  sll_f_spline_pp_horner_1d(degree, &
            self%spline_pp1d_1(component)%poly_coeffs_fpa, xi(component), i) &
            * self%delta_x(component)
       self%spline_1(i,2) =  sll_f_spline_pp_horner_1d(degree, &
            self%spline_pp1d_1(component)%poly_coeffs_fpa, xnew, i) &
            * self%delta_x(component)
    end do

    if (position_old(component) .le. position_new ) then
       self%j1d(local_size-degree+1:local_size) = self%spline_1(1:degree,2)
       self%j1d(1:local_size-degree) = self%delta_x(component)
       self%j1d(1:degree) = self%j1d(1:degree) - self%spline_1(1:degree,1)
    else
       self%j1d(1:local_size-degree) = -self%delta_x(component)
       self%j1d(local_size-degree+1:local_size) = -self%spline_1(1:degree,1)
       self%j1d(1:degree) = self%j1d(1:degree) + self%spline_1(1:degree,2)
    end if
    !----

    ! Achtung wir brauchen nun splines von beidem Grad
    self%spline_0 = 0.0_f64
    self%spline_1 = 0.0_f64
    do j=1,3
       if (j .ne. component ) then
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_0(j), self%spline_0(1:self%spline_degree(j)+1,j), self%spline_degree(j), xi(j))
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_1(j), self%spline_1(1:self%spline_degree(j),j), self%spline_degree(j)-1, xi(j))
       end if
    end do

    box(1:2) = box(1:2) - self%spline_degree(1:2)

    ! Define the range of the first component
    if (boxold<boxnew) then
       box(component) = boxold-degree+1
    else
       box(component) = boxnew-degree+1
    end if


    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )
    do k = 1, self%spline_degree(2)+1     
       vtt2 = 0.0_f64
       vtt3 = 0.0_f64
       index3d(2) = self%index1d(k,2)!modulo(box(2)+k-2,self%n_cells(2))+1
       do j = 1, self%spline_degree(1)+1 
          index3d(1) = self%index1d(j,1)!modulo(box(1)+j-2,self%n_cells(1))+1

          startjk  = index3d(1) + (index3d(2)-1)*self%n_cells(1)
          splineij = self%spline_0(j, 1) * self%spline_0(k,2) * marker_charge
          
          vt = 0.0_f64
          do i = 1, local_size
             index3d(3) = modulo(box(component)+i-2,self%n_cells(component))
             index1d = startjk + index3d(3)*stride
             j_dofs(index1d) = j_dofs(index1d) + self%j1d(i) * &
                  splineij  

             vt(1) = vt(1) + bfield_dofs(start1+index1d)*self%j1d(i)
             vt(2) = vt(2) + bfield_dofs(start2+index1d)*self%j1d(i)
          end do

          if (j>1) then
             vtt2 = vtt2 + vt(1)*self%spline_1(j-1, 1)
          end if
          vtt3 = vtt3 + vt(2)*self%spline_0(j, 1)


       end do
       vi(1) = vi(1) - qoverm*vtt2*self%spline_0(k, 2)
       if ( k> 1) then
          vi(2) = vi(2) + qoverm*vtt3*self%spline_1(k-1, 2)
       end if
    end do

  end subroutine add_current_update_v_primitive_component3_spline_3d_feec


  !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting), version based on quadrature
  subroutine add_current_update_v_component1_spline_3d_feec (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_3d_feec), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(3) !< Position at time t
    sll_real64, intent(in)    :: position_new !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(:) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(3) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(:) !< Coefficients of current expansion
    !local variables
    sll_int32 :: box(3), boxnew, boxold
    sll_real64 :: xi(3)
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k
    sll_real64 :: xnew
    sll_real64 :: vt(2), vtt2, vtt3
    sll_int32  :: start1, start2
    sll_int32  :: component, local_size
    sll_int32  :: startjk, stride
    sll_real64 :: splinejk
    sll_int32  :: degree

    component = 1
    start1 = 2*self%n_total
    start2 = self%n_total
    stride = 1
    degree=self%spline_degree(component)
    call convert_x_to_xbox( self, position_old, xi, box )
    call convert_x_to_xbox_1d( self, component, position_new, xnew, boxnew )
    boxold = box(component)

    ! Achtung wir brauchen nun splines von beidem Grad
    self%spline_0 = 0.0_f64
    self%spline_1 = 0.0_f64
    do j=1,3
       if (j .ne. component ) then
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_0(j), self%spline_0(1:self%spline_degree(j)+1,j), self%spline_degree(j), xi(j))
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_1(j), self%spline_1(1:self%spline_degree(j),j), self%spline_degree(j)-1, xi(j))
       end if
    end do

    !box(component) = box(component)-degree+1
    box(2:3) = box(2:3) - self%spline_degree(2:3)
    local_size = abs(boxnew-boxold)+degree

    ! Define the range of the first component
    if (boxold<boxnew) then
       box(component) = boxold-degree+1
    else
       box(component) = boxnew-degree+1
    end if
  
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 ) 
    do k = 1, self%spline_degree(3)+1      
       vtt2 = 0.0_f64
       vtt3 = 0.0_f64
       index3d(3) = self%index1d(k,3)!modulo(box(3)+k-2,self%n_cells(3))+1
       do j = 1, self%spline_degree(2)+1
          index3d(2) = self%index1d(j,2)!modulo(box(2)+j-2,self%n_cells(2))+1

          startjk  = 1 + (index3d(2)-1)*self%n_cells(1) + &
               (index3d(3)-1)*self%n_cells(1)*self%n_cells(2)
          splinejk = self%spline_0(j, 2) * self%spline_0(k,3)

          vt = 0.0_f64
          self%j1d = 0.0_f64
          call  add_current_1d(self, component, xi(component), boxold-1, xnew, boxnew-1, marker_charge, bfield_dofs, start1+startjk, start2+startjk, stride, vt, self%j1d)
          if (j>1) then
             vtt2 = vtt2 + vt(1)*self%spline_1(j-1, 2)
          end if
          vtt3 = vtt3 + vt(2)*self%spline_0(j, 2)

          do i = 1, local_size
             index3d(1) = modulo(box(1)+i-2,self%n_cells(1))+1
            
             index1d = startjk+index3d(1)-1!convert_index_3d_to_1d( index3d, self%n_cells ) 
             j_dofs(index1d) = j_dofs(index1d) + self%j1d(index3d(1)) * &
                  splinejk 
          end do
       end do
       vi(2) = vi(2) - qoverm*vtt2*self%spline_0(k, 3)
       if ( k> 1) then
          vi(3) = vi(3) + qoverm*vtt3*self%spline_1(k-1, 3)
       end if
    end do



  end subroutine add_current_update_v_component1_spline_3d_feec

  !> Add current for one particle and update v (according to H_p2 part in Hamiltonian splitting), version based on quadrature
  subroutine add_current_update_v_component2_spline_3d_feec (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_3d_feec), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(3) !< Position at time t
    sll_real64, intent(in)    :: position_new !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(:) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(3) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(:) !< Coefficients of current expansion
    !local variables
    sll_int32 :: box(3), boxnew, boxold
    sll_real64 :: xi(3)
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k
    sll_real64 :: xnew
    sll_real64 :: vt(2), vtt2, vtt3
    sll_int32  :: start1, start2
    sll_int32  :: component, local_size
    sll_int32 :: stride, startjk
    sll_real64 :: splineik
    sll_int32  :: degree

    component = 2
    start1 = 2*self%n_total 
    start2 = 0
    stride = self%n_cells(1)
    degree=self%spline_degree(component)

    call convert_x_to_xbox( self, position_old, xi, box )
    call convert_x_to_xbox_1d( self, component, position_new, xnew, boxnew )
    boxold = box(component)


    ! Achtung wir brauchen nun splines von beidem Grad
    self%spline_0 = 0.0_f64
    self%spline_1 = 0.0_f64
    do j=1,3
       if (j .ne. component ) then
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_0(j), self%spline_0(1:self%spline_degree(j)+1,j), self%spline_degree(j), xi(j))
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_1(j), self%spline_1(1:self%spline_degree(j),j), self%spline_degree(j)-1, xi(j))
       end if
    end do

    box(1) = box(1) - self%spline_degree(1)
    box(3) = box(3) - self%spline_degree(3)
    local_size = abs(boxnew-boxold)+degree
    
    ! Define the range of the second component
    if (boxold<boxnew) then
       box(component) = boxold-degree+1
    else
       box(component) = boxnew-degree+1
    end if


    call box_index( self, box(1), 1 )
    call box_index( self, box(3), 3 ) 
    do k = 1, self%spline_degree(3)+1      
       vtt2 = 0.0_f64
       vtt3 = 0.0_f64
       index3d(3) = self%index1d(k,3)!modulo(box(3)+k-2,self%n_cells(3))+1
       do j = 1, self%spline_degree(1)+1
          index3d(1) = self%index1d(j,1)!modulo(box(1)+j-2,self%n_cells(1))+1
         
          startjk  = index3d(1) + &
               (index3d(3)-1)*self%n_cells(1)*self%n_cells(2)
          splineik = self%spline_0(j, 1) * self%spline_0(k,3) 

          vt = 0.0_f64
          self%j1d = 0.0_f64
          call  add_current_1d(self, component, xi(component), boxold-1, xnew, boxnew-1, marker_charge, bfield_dofs,  start1+startjk, start2+startjk, stride, vt, self%j1d)

          if (j>1) then
             vtt2 = vtt2 + vt(1)*self%spline_1(j-1, 1)
          end if
          vtt3 = vtt3 + vt(2)*self%spline_0(j, 1)

          do i = 1, local_size
             index3d(2) = modulo(box(2)+i-2,self%n_cells(2))+1
             index1d = startjk+(index3d(2)-1)*stride!convert_index_3d_to_1d( index3d, self%n_cells )
             j_dofs(index1d) = j_dofs(index1d) + self%j1d(index3d(2)) * &
                  splineik

          end do
       end do
       vi(1) = vi(1) + qoverm*vtt2*self%spline_0(k, 3)
       if ( k> 1) then
          vi(3) = vi(3) - qoverm*vtt3*self%spline_1(k-1, 3)
       end if
    end do


  end subroutine add_current_update_v_component2_spline_3d_feec

  
  !> Add current for one particle and update v (according to H_p3 part in Hamiltonian splitting), version based on quadrature
  subroutine add_current_update_v_component3_spline_3d_feec (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_3d_feec), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(3) !< Position at time t
    sll_real64, intent(in)    :: position_new !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(:) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(3) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(:) !< Coefficients of current expansion
    !local variables
    sll_int32 :: box(3), boxnew, boxold
    sll_real64 :: xi(3)
    sll_int32  :: index3d(3)
    sll_int32  :: index1d
    sll_int32  :: i,j,k
    sll_real64 :: xnew
    sll_real64 :: vt(2), vtt2, vtt3
    sll_int32  :: start1, start2
    sll_int32  :: component, local_size
    sll_int32  :: stride, startjk
    sll_real64 :: splineij
    sll_int32  :: degree

    component = 3
    start1 = self%n_total
    start2 = 0
    stride = self%n_cells(1)*self%n_cells(2)
    degree=self%spline_degree(component)

    call convert_x_to_xbox( self, position_old, xi, box )
    call convert_x_to_xbox_1d( self, component, position_new, xnew, boxnew )
    boxold = box(component)

    ! Achtung wir brauchen nun splines von beidem Grad
    self%spline_0 = 0.0_f64
    self%spline_1 = 0.0_f64
    do j=1,3
       if (j .ne. component ) then
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_0(j), self%spline_0(1:self%spline_degree(j)+1,j), self%spline_degree(j), xi(j))
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_1(j), self%spline_1(1:self%spline_degree(j),j), self%spline_degree(j)-1, xi(j))
       end if
    end do

    box(1:2) = box(1:2) -  self%spline_degree(1:2)
    local_size = abs(boxnew-boxold)+degree
    ! Define the range of the third component
    if (boxold<boxnew) then
       box(component) = boxold-degree+1
    else
       box(component) = boxnew-degree+1
    end if

    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )
    do k = 1, self%spline_degree(2)+1      
       vtt2 = 0.0_f64
       vtt3 = 0.0_f64
       index3d(2) = self%index1d(k,2)!modulo(box(2)+k-2,self%n_cells(2))+1
       do j = 1, self%spline_degree(1)+1
          index3d(1) = self%index1d(j,1)!modulo(box(1)+j-2,self%n_cells(1))+1

          startjk  = index3d(1) + (index3d(2)-1)*self%n_cells(1)
          splineij = self%spline_0(j, 1) * self%spline_0(k,2)

          vt = 0.0_f64
          self%j1d = 0.0_f64
          call  add_current_1d(self, component, xi(component), boxold-1, xnew, boxnew-1, marker_charge, bfield_dofs, start1+startjk, start2+startjk, stride, vt, self%j1d)
          if (j>1) then
             vtt2 = vtt2 + vt(1)*self%spline_1(j-1, 1)
          end if
          vtt3 = vtt3 + vt(2)*self%spline_0(j, 1)

          do i = 1, local_size
             index3d(3) = modulo(box(component)+i-2,self%n_cells(component))+1
             index1d = startjk + (index3d(3)-1)*stride!convert_index_3d_to_1d( index3d, self%n_cells )

             j_dofs(index1d) = j_dofs(index1d) + self%j1d(index3d(3)) * &
                  splineij  
          end do
       end do
       vi(1) = vi(1) - qoverm*vtt2*self%spline_0(k, 2)
       if ( k> 1) then
          vi(2) = vi(2) + qoverm*vtt3*self%spline_1(k-1, 2)
       end if
    end do

  end subroutine add_current_update_v_component3_spline_3d_feec


  !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting)
  subroutine add_current_1d (self, component, r_old, index_old,  r_new, index_new, marker_charge, bfield_dofs, start1, start2, stride, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_3d_feec), intent(inout) :: self !< kernel smoother object
    sll_int32,  intent(in)    :: component
    sll_real64, intent(in)    :: r_old
    sll_int32,  intent(in)    :: index_old
    sll_real64, intent(in)    :: r_new
    sll_int32,  intent(in)    :: index_new
    !sll_real64, intent(in)    :: position_old(3) !< Position at time t
    !sll_real64, intent(in)    :: position_new(3) !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: bfield_dofs(self%n_total*3) !< Coefficient of B-field expansion
    sll_int32,                              intent(in)    :: start1 !< Index where B-field component starts that updates vi(1)
    sll_int32,                              intent(in)    :: start2 !< Index where B-field component starts that updates vi(2)
    sll_int32,                              intent(in)    :: stride !< Stride for index of B-field
    sll_real64, intent(inout) :: vi(2) !< Velocity of the particles to be updated
    sll_real64, intent(inout) :: j_dofs(self%n_total) !< Coefficients of current expansion
    ! local variables
    sll_int32  ::  ind

    if (index_old == index_new) then
       if (r_old < r_new) then
          call update_jv(self, component, r_old, r_new, index_old, marker_charge, &
               1.0_f64, bfield_dofs, start1, start2, stride, vi, j_dofs)
       else
          call update_jv(self, component, r_new, r_old, index_old, marker_charge,  &
               -1.0_f64, bfield_dofs, start1, start2, stride, vi, j_dofs)
       end if
    elseif (index_old < index_new) then
       call update_jv (self, component, r_old, 1.0_f64, index_old, marker_charge, &
            1.0_f64, bfield_dofs, start1, start2, stride, vi, j_dofs)
       call update_jv (self, component, 0.0_f64, r_new, index_new, marker_charge, &
            1.0_f64, bfield_dofs, start1, start2, stride, vi, j_dofs)
       do ind = index_old+1, index_new-1
          call update_jv (self, component, 0.0_f64, 1.0_f64, ind, marker_charge, &
               1.0_f64, bfield_dofs, start1, start2, stride, vi, j_dofs)
       end do
    else
       call update_jv (self, component, r_new, 1.0_f64, index_new, marker_charge, &
            -1.0_f64, bfield_dofs, start1, start2, stride, vi, j_dofs)
       call update_jv (self, component, 0.0_f64, r_old, index_old, marker_charge, &
            -1.0_f64, bfield_dofs, start1, start2, stride, vi, j_dofs)
       do ind = index_new+1, index_old-1
          call update_jv (self, component, 0.0_f64, 1.0_f64, ind, marker_charge, &
               -1.0_f64, bfield_dofs, start1, start2, stride, vi, j_dofs)
       end do
    end if


  end subroutine add_current_1d



  !> Helper function for \a add_current_update_v.
  subroutine update_jv(self, component, lower, upper, index, marker_charge, sign, bfield_dofs, start1, start2, stride, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_3d_feec), intent(inout) :: self !< time splitting object 
    sll_int32,                              intent(in)    :: component
    sll_real64,                             intent(in)    :: lower
    sll_real64,                             intent(in)    :: upper
    sll_int32,                              intent(in)    :: index
    sll_real64,                             intent(in)    :: marker_charge
    sll_real64,                             intent(in)    :: sign
    sll_real64,                             intent(inout) :: vi(2)
    sll_real64,                             intent(in)    :: bfield_dofs(self%n_total*3)
    sll_int32,                              intent(in)    :: start1
    sll_int32,                              intent(in)    :: start2
    sll_int32,                              intent(in)    :: stride
    sll_real64,                             intent(inout) :: j_dofs(self%n_total)
    !Local variables
    sll_int32  :: ind, i_grid, i_mod, j
    sll_real64 :: c1, c2
    sll_int32 :: spline_degree
    sll_int32 :: n_cells

    self%spline_val(:,1) = 0._f64
    self%spline_val_more(:,1) = 0._f64
    
    n_cells = self%n_cells(component)
    spline_degree = self%spline_degree(component)-1

    c1 =  0.5_f64*(upper-lower)
    c2 =  0.5_f64*(upper+lower)

    call sll_s_uniform_bsplines_eval_basis(spline_degree, c1*self%quad_xw(1,1)+c2, &
         self%spline_val(1:spline_degree+1,1))
    self%spline_val(:,1) = self%spline_val(:,1) * (self%quad_xw(2,1)*c1)
    do j = 2, self%n_quad_points
       call sll_s_uniform_bsplines_eval_basis(spline_degree, c1*self%quad_xw(1,j)+c2, &
            self%spline_val_more(1:spline_degree+1,1))
       self%spline_val(:,1) = self%spline_val(:,1) + self%spline_val_more(:,1) * (self%quad_xw(2,j)*c1)
    end do
    self%spline_val(:,1) = self%spline_val(:,1) * (sign*self%delta_x(component))

    ind = 1
    do i_grid = index - spline_degree , index
       i_mod = modulo(i_grid, n_cells )
       j_dofs(i_mod+1) = j_dofs(i_mod+1) + &
            (marker_charge*self%spline_val(ind,1))
       vi(1) = vi(1) + self%spline_val(ind,1)*bfield_dofs(i_mod*stride+start1)
       vi(2) = vi(2) + self%spline_val(ind,1)*bfield_dofs(i_mod*stride+start2)
       ind = ind + 1
    end do

  end subroutine update_jv



  !> Constructor
  subroutine init_spline_3d_feec ( self, n_cells, domain, spline_degree, no_particles )
    class (sll_t_particle_mesh_coupling_spline_3d_feec), intent( out ) :: self !< Kernel smoother object 
    sll_int32,                               intent(in) :: n_cells(3) !< number of DoFs (spline coefficients)
    sll_real64,                              intent(in) :: domain(3,2) !< x_min and x_max of the domain
    sll_int32,                               intent(in) :: spline_degree(3) !< Degree of smoothing kernel spline
    sll_int32, optional,                     intent(in) :: no_particles !< number of particles
    !local variables
    sll_int32 :: maxspan, j, maxind

    ! Store grid information
    self%domain = domain
    self%n_cells = n_cells
    self%n_total = product(n_cells)
    self%delta_x = (self%domain(:,2)-self%domain(:,1))/real(n_cells, f64)
    self%rdelta_x = 1.0_f64/self%delta_x

    ! Store basis function information
    if( present(no_particles) ) then
       self%no_particles = no_particles
    else
       self%no_particles = 1
    end if

    ! Initialize information on the spline
    self%spline_degree = spline_degree
    maxspan = maxval(spline_degree + 1)

!!!Besser n_quad_points als vektor
    self%n_quad_points = (maxval(self%spline_degree)+2)/2
    allocate( self%quad_xw(2,self%n_quad_points) )
    ! normalized Gauss Legendre points and weights
    self%quad_xw = sll_f_gauss_legendre_points_and_weights(self%n_quad_points)

    ! For the line integral
    self%n_quad_points_line = sum(self%spline_degree) -1 ! Degree of the polynomial in time for the line integral
    self%n_quad_points_line = ceiling(real(self%n_quad_points_line+1, f64)*0.5_f64) ! Compute the number of Gauss points needed to get exact integrals for this polynomial

    !print*, 'Size line quadrature', self%n_quad_points_line
    !self%n_quad_points_line = 3
    allocate( self%quad_xw_line(2,self%n_quad_points_line) )
    ! normalized Gauss Legendre points and weights
    self%quad_xw_line = sll_f_gauss_legendre_points_and_weights(self%n_quad_points_line)

    self%quad_xw_line(1,:) = 0.5_f64*(self%quad_xw_line(1,:)+1.0_f64)
    self%quad_xw_line(2,:) = 0.5_f64*self%quad_xw_line(2,:)


    allocate( self%spline_val(maxspan,3) )
    allocate( self%spline_val_more(maxspan,3) )
    allocate( self%spline_0(maxspan,3) )
    allocate( self%spline_1(maxspan-1,3) )
    allocate( self%spline_0_deriv(maxspan,3) )
    allocate( self%spline_1_deriv(maxspan-1,3) )
    allocate( self%spline1_0(maxspan,3) )
    allocate( self%j1d( maxval(self%n_cells)+maxval(self%spline_degree) ))
    allocate( self%spline_2d(maxspan, maxspan) )
    allocate( self%spline1_2d(1:maxspan, 1:maxspan) )
    allocate( self%spline1_3d(1:maxspan, 1:maxspan, 1:maxspan) )
    allocate( self%spline2_3d(1:maxspan, 1:maxspan, 1:maxspan) )
    allocate( self%index1d(maxspan, 3) )

    call sll_s_spline_pp_init_3d(self%spline_pp_0, spline_degree, n_cells)
    call sll_s_spline_pp_init_3d(self%spline_pp_1, spline_degree-1, n_cells)
    do j = 1, 3
       call sll_s_spline_pp_init_1d( self%spline_pp1d_0(j), spline_degree(j), n_cells(j) )
       call sll_s_spline_pp_init_1d( self%spline_pp1d_1(j), spline_degree(j)-1, n_cells(j) )
    end do

  end subroutine init_spline_3d_feec


  !> Destructor
  subroutine free_spline_3d_feec(self)
    class (sll_t_particle_mesh_coupling_spline_3d_feec), intent( inout ) :: self !< Kernel smoother object 
    !local variable
    sll_int32 :: j
    
    deallocate( self%quad_xw)
    deallocate( self%spline_val)
    deallocate( self%spline_val_more)
    deallocate( self%spline_0)
    deallocate( self%spline_1)
    deallocate( self%j1d)
    deallocate( self%spline1_0 )
    deallocate( self%spline_2d )
    deallocate( self%spline1_2d )
    deallocate( self%spline1_3d )
    deallocate( self%spline2_3d )

    call sll_s_spline_pp_free_3d( self%spline_pp_0 )
    call sll_s_spline_pp_free_3d( self%spline_pp_1 )

    do j = 1, 3
       call sll_s_spline_pp_free_1d( self%spline_pp1d_0(j) )
       call sll_s_spline_pp_free_1d( self%spline_pp1d_1(j) )
    end do
    
  end subroutine free_spline_3d_feec


  !> Evaluates the integral int_{poisition_old}^{position_new} field(x) d x and the integrated current
  subroutine add_current_evaluate_int ( self, position_old, position_new, vbar, bfield_dofs, j_dofs, bfield_val)
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: vbar(3)
    sll_real64,                               intent( in )    :: bfield_dofs(:)
    sll_real64,                               intent( inout ) :: j_dofs(:)
    sll_real64,                               intent( out   ) :: bfield_val(3)

    sll_real64 :: xold(3), xnew(3), xnewtilde(3), xbox(3)
    sll_int32 :: boxold(3), boxnew(3), sigma_counter(3), boxdiff(3), increment(3), box(3)
    sll_real64 :: sigma_r, sigma_l, sigma, sigma_next(3), field_value(3), weight
    sll_int32 :: j, q, bl, maxind, index_is
    type(vector) :: sigmas(3)

    call convert_x_to_xbox( self, position_old, xold, boxold )
    call convert_x_to_xbox( self, position_new, xnew, boxnew )

    ! Now we need to compute the normalized 1d line along which to integrate:
    boxdiff = boxnew-boxold
    xnewtilde = xnew + real(boxdiff,f64)

    sigma_l = 0.0_f64
    box = boxold ! We start in the old box


    sigma_counter = 0

    bfield_val = 0.0_f64

    do j=1,3
       if (boxdiff(j) > 0 ) then
          allocate ( sigmas(j)%vals(boxdiff(j)+1) )
          do bl=1,boxdiff(j)
             sigmas(j)%vals(bl) = (real(bl,f64)  - xold(j))/(xnewtilde(j)-xold(j))
          end do
          sigmas(j)%vals(boxdiff(j)+1) = 1.0_f64
          sigma_next(j) = sigmas(j)%vals(1)
          increment(j) = 1
       elseif (boxdiff(j) < 0 ) then
          allocate ( sigmas(j)%vals(-boxdiff(j)+1) )
          do bl=boxdiff(j)+1,0
             sigmas(j)%vals(-bl+1) = (real(bl,f64)  - xold(j))/(xnewtilde(j)-xold(j))
          end do
          sigmas(j)%vals(-boxdiff(j)+1) = 1.0_f64
          sigma_next(j) = sigmas(j)%vals(1)
          increment(j) = -1
       else
          sigma_next(j) = 1.0_f64
          increment(j) = 0
       end if
    end do

    sigma_r = 0.0_f64
    do while ( sigma_r < 1.0_f64 )
       ! Identify index of next intersection
       index_is = minloc(sigma_next, dim=1)
       sigma_r = sigma_next(index_is)

       do q = 1, self%n_quad_points_line
          sigma = sigma_l + (sigma_r-sigma_l) * self%quad_xw_line(1,q)
          xbox = xold* (1.0_f64 - sigma) + xnewtilde * sigma  - real(sigma_counter*increment, f64)
          if (maxval(xbox)> 1.0_f64 ) then
             print*, xbox, sigma, sigma_counter, increment
             SLL_ERROR( 'add_current_evaluate', 'box value too large')
          elseif (minval(xbox)< 0.0_f64 ) then
             print*, xbox, sigma, sigma_counter, increment
             print*, xold, xnewtilde, sigma_r
             SLL_ERROR( 'add_current_evaluate', 'box value too low')
          end if

          weight = self%quad_xw_line(2,q)*(sigma_r-sigma_l)

          call point_add_eval_subcyc (self, box, xbox, bfield_dofs, vbar*weight, j_dofs, field_value )
          bfield_val = bfield_val + field_value * weight * sigma !( 1.0_f64 - sigma )
       end do
       if (sigma_r < 1.0_f64 ) then
          ! Update the
          sigma_counter(index_is) = sigma_counter(index_is)+1
          sigma_next(index_is) = sigmas(index_is)%vals(sigma_counter(index_is)+1)
          box(index_is) = box(index_is) + increment(index_is)
          sigma_l = sigma_r
       end if
    end do

  end subroutine add_current_evaluate_int


  !> Helper function for add_current_evaluate_int, takes care of per cell computations
  subroutine point_add_eval_subcyc ( self, box_in, xbox, field_dofs, weight, j_dofs, field_value )
    class( sll_t_particle_mesh_coupling_spline_3d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_int32, intent(in) :: box_in(3)
    sll_real64, intent(in) :: xbox(3)
    sll_real64, intent(in) :: field_dofs(:)
    sll_real64, intent(in) :: weight(3)
    sll_real64, intent(inout) :: j_dofs(:)
    sll_real64, intent(out) :: field_value(3)



    sll_int32 :: i,j,k, index1d, n_total
    sll_int32 :: box(3), index3d(3)
    sll_real64 :: spval

    n_total = self%n_total
    field_value = 0.0_f64

    do j=1,3
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree(j), xbox(j), &
            self%spline_0(1:self%spline_degree(j)+1,j) )
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree(j)-1, xbox(j), &
            self%spline_1(1:self%spline_degree(j),j) )
    end do

    box = box_in-self%spline_degree

    !box(1) = box(1)+1
    call box_index( self, box(1), 1 )
    call box_index( self, box(2), 2 )   
    call box_index( self, box(3), 3 )

    ! First component j_dofs
    do k=1,self%spline_degree(3)+1 
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j=1,self%spline_degree(2)+1
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i=2,self%spline_degree(1)+1
             index1d = index3d(2) + self%index1d(i,1)

             spval =  self%spline_1(i-1,1) * &
                  self%spline_0(j,2) * self%spline_0(k,3)
             j_dofs(index1d) = j_dofs(index1d) + &
                  weight(1) * spval
          end do
       end do
    end do
    ! First component field
    do k=2,self%spline_degree(3)+1 
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j=2,self%spline_degree(2)+1
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i=1,self%spline_degree(1)+1
             index1d = index3d(2) + self%index1d(i,1)
             spval =  self%spline_0(i,1) * &
                  self%spline_1(j-1,2) * self%spline_1(k-1,3)
             field_value(1) = field_value(1) + field_dofs(index1d) * spval
          end do
       end do
    end do
    ! Second component j_dofs
    do k=1,self%spline_degree(3)+1       
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j=2,self%spline_degree(2)+1
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i=1,self%spline_degree(1)+1
             index1d = index3d(2) + self%index1d(i,1) + n_total

             spval = self%spline_0(i,1) * &
                  self%spline_1(j-1,2) * self%spline_0(k,3)
             j_dofs(index1d) = j_dofs(index1d) + &
                  weight(2) * spval
             field_value(2) = field_value(2) + field_dofs(index1d ) * spval
          end do
       enddo
    end do
    ! Second component field
    do k=2,self%spline_degree(3)+1       
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j=1,self%spline_degree(2)+1
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i=2,self%spline_degree(1)+1
             spval = self%spline_1(i-1,1) * &
                  self%spline_0(j,2) * self%spline_1(k-1,3)
             field_value(2) = field_value(2) + field_dofs(index1d ) * spval
          end do
       end do
    end do

    n_total = n_total*2

    ! Third component j_dofs
    do k=2,self%spline_degree(3) +1   
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j=1,self%spline_degree(2)+1
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i=1,self%spline_degree(1)+1
             index1d = index3d(2) + self%index1d(i,1) + n_total
             spval = self%spline_0(i,1) * &
                  self%spline_0(j,2) * self%spline_1(k-1,3)
             j_dofs(index1d) = j_dofs(index1d) + &
                  weight(3) * spval
          end do
       end do
    end do
    ! Third component field
    do k=1,self%spline_degree(3)+1    
       index3d(3) = (self%index1d(k,3)-1)*self%n_cells(1)*self%n_cells(2)
       do j=2,self%spline_degree(2)+1
          index3d(2) = index3d(3) + (self%index1d(j,2)-1)*self%n_cells(1)
          do i=2,self%spline_degree(1)+1
             index1d = index3d(2) + self%index1d(i,1) + n_total
             spval = self%spline_1(i-1,1) * &
                  self%spline_1(j-1,2) * self%spline_0(k,3)
             field_value(3) = field_value(3) + field_dofs(index1d ) * spval
          end do
       end do
    end do

  end subroutine point_add_eval_subcyc

end module sll_m_particle_mesh_coupling_spline_3d_feec
