!> @ingroup particle_mesh_coupling
!> @author Katharina Kormann, IPP
!> @brief Particle mesh coupling for 3d with splines of arbitrary degree placed on a uniform tensor product mesh.
!> @details Spline with index i starts at point i

module sll_m_particle_mesh_coupling_spline_2d_feec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_low_level_bsplines, only: &
    sll_s_uniform_bsplines_eval_basis

  use sll_m_gauss_legendre_integration, only : &
    sll_f_gauss_legendre_points_and_weights

  use sll_m_splines_pp, only : &
    sll_f_spline_pp_horner_2d, &
    sll_s_spline_pp_free_3d, &
    sll_s_spline_pp_horner_m_3d, &
    sll_s_spline_pp_init_3d, &
    sll_t_spline_pp_3d, &
    sll_t_spline_pp_1d, &
    sll_s_spline_pp_init_1d, &
    sll_f_spline_pp_horner_1d, &
    sll_s_spline_pp_horner_m_1d

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base


  implicit none

  public :: &
    sll_t_particle_mesh_coupling_spline_2d_feec

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Particle mesh coupling in 3d based on (arbitrary degree) spline on a tensor product uniform mesh
  type :: sll_t_particle_mesh_coupling_spline_2d_feec
     ! Information about the 3d mesh
     sll_real64 :: delta_x(2)  !< Value of grid spacing along both directions.
     sll_real64 :: domain(2,2) !< Definition of the domain: domain(1,1) = x1_min  domain(1,2) = x1_max
     sll_real64 :: rdelta_x(2)  !< Inverse values of delta_x
     
     ! Information about the particles
     sll_int32  :: no_particles !< Number of particles of underlying PIC method (processor local)

     ! 
     sll_int32  :: spline_degree !< Degree of smoothing kernel spline
     sll_int32  :: n_span !< Number of intervals where spline non zero (spline_degree + 1)
     sll_int32  :: n_quad_points !< Number of quadrature points

     sll_real64, allocatable :: spline_val(:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_0(:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_1(:,:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_val_more(:,:) !< more scratch data for spline evaluation
     sll_real64, allocatable :: spline_2d(:,:)
     sll_real64, allocatable :: j1d(:) !< scratch data for 1d slices of j
     sll_real64, allocatable :: quad_xw(:,:) !< quadrature weights and points


     sll_int32              :: dim
     sll_int32              :: n_dofs  !< Number of degrees of freedom of the smoothing kernels.
     sll_int32              :: n_grid(2) !< Number of grid points per dimension for use on tensor product grid based smoothing kernels.
     
     type(sll_t_spline_pp_1d) :: spline_pp1d_0
     type(sll_t_spline_pp_1d) :: spline_pp1d_1

   contains

     procedure :: add_charge => add_charge_single_spline_2d_feec !> Add charge of one particle
   !  procedure :: add_charge_pp => add_charge_single_spline_pp_3d_feec !> Add charge of one particle
     procedure :: add_current => add_current_spline_2d_feec
     procedure :: add_current_update_v_component1 => add_current_update_v_primitive_component1_spline_2d_feec !> Add current of one particle
     procedure :: add_current_update_v_component2 => add_current_update_v_primitive_component2_spline_2d_feec !> Add current of one particle
     !procedure :: add_current_update_v_component3 => add_current_update_v_primitive_component3_spline_3d_feec !> Add current of one particle
     procedure :: evaluate => evaluate_field_single_spline_2d_feec !> Evaluate spline function with given coefficients
     procedure :: evaluate_pp => evaluate_field_single_spline_pp_2d_feec !> Evaluate spline function with given coefficients
     procedure :: evaluate_multiple => evaluate_multiple_spline_2d_feec !> Evaluate multiple spline functions with given coefficients
     !procedure :: update_jv !> helper function to compute the integral of j using Gauss quadrature
     procedure :: add_particle_mass => add_particle_mass_spline_2d_feec
     procedure :: init => init_spline_2d_feec !> Constructor
     procedure :: free => free_spline_2d_feec !> Destructor
  end type sll_t_particle_mesh_coupling_spline_2d_feec


contains


  !---------------------------------------------------------------------------!
  !> Add charge of one particle
  subroutine add_charge_single_spline_2d_feec(self, position, marker_charge, degree, rho_dofs)
    class( sll_t_particle_mesh_coupling_spline_2d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_int32,                                intent( in )    :: degree(2) !< Spline degree along each dimension
    sll_real64,                               intent( inout ) :: rho_dofs(self%n_dofs) !< Coefficient vector of the charge distribution

    sll_int32 :: box(2)
    sll_real64 :: xi(2)
    sll_int32  :: index2d(2)
    sll_int32  :: index1d
    sll_int32  :: i,j

    
    call convert_x_to_xbox( self, position, xi, box )
    
    self%spline_0 = 0.0_f64
    do j=1,2
          call sll_s_uniform_bsplines_eval_basis( degree(j), xi(j), self%spline_0(:,j) )
    end do
    
    ! Build scaling with marker_charge into spline along first dimension
    self%spline_0(:,1) = self%spline_0(:,1)*marker_charge
    
    box = box-degree
    do j=1,degree(2)+1
       index2d(2) = modulo(box(2)+j-2,self%n_grid(2))+1
       do i=1,degree(1)+1
          index2d(1) = modulo(box(1)+i-2,self%n_grid(1))+1
          index1d = convert_index_2d_to_1d( index2d, self%n_grid )
          rho_dofs(index1d) = rho_dofs(index1d) + &
               (self%spline_0(i,1) * self%spline_0(j,2))
       end do
    end do


  end subroutine add_charge_single_spline_2d_feec
  !---------------------------------------------------------------------------!
  !> Add charge of one particle
  subroutine add_particle_mass_spline_2d_feec(self, position, marker_charge, degree, particle_mass )
    class( sll_t_particle_mesh_coupling_spline_2d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_int32,                                intent( in )    :: degree(2) !< Spline degree along each dimension
    sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution

    sll_int32 :: box(2)
    sll_real64 :: xi(2)
    sll_int32  :: index2d(2)
    sll_int32  :: index1d
    sll_int32  :: i,j, col1, col2, col3, ind
    sll_real64 :: splineijk, splineijkcol3

    
    call convert_x_to_xbox( self, position, xi, box )
    

    ! Old version based on arbitrary degree splines
    self%spline_0 = 0.0_f64
    do j=1,3
       call sll_s_uniform_bsplines_eval_basis( degree(j), xi(j), self%spline_0(:,j) )
    end do
       
    ! 2d array combining first and second dimension
    do j=1,degree(2)+1
       do i=1,degree(1)+1
          self%spline_2d(i,j) = self%spline_0(i,1)*self%spline_0(j,2)
       end do
    end do
    ! TODO: Check if also 3d array would make sense

    
    box = box-degree
    do j=1,degree(2)+1
       index2d(2) = modulo(box(2)+j-2,self%n_grid(2))+1
       do i=1,degree(1)+1
          index2d(1) = modulo(box(1)+i-2,self%n_grid(1))+1
          index1d = convert_index_2d_to_1d( index2d, self%n_grid )
          ind = 1+(degree(1)+1-i)+(degree(2)+1-j)*(2*degree(1)+1)
          splineijk = marker_charge * self%spline_2d(i,j) 
          do col2 = 1,degree(2)+1
             do col1 = 1,degree(1)+1                     
                particle_mass(ind, index1d) = &
                     particle_mass( ind, index1d) + &
                     splineijk * &
                     self%spline_2d(col1,col2)
                ind = ind+1
             end do
             ind = ind+degree(1)
          end do
          ind = ind+( degree(2) ) * (2*degree(1)+1)
       end do
    end do
      
 


  end subroutine add_particle_mass_spline_2d_feec
  

!!$  !---------------------------------------------------------------------------!
!!$  !> Add charge of one particle
!!$  subroutine add_charge_single_spline_pp_3d_feec(self, position, marker_charge, degree, rho_dofs)
!!$    class( sll_t_particle_mesh_coupling_spline_2d_feec ), intent(inout)   :: self !< kernel smoother object
!!$    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
!!$    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
!!$    sll_int32,                                intent( in )    :: degree(3) !< Spline degree along each dimension
!!$    sll_real64,                               intent( inout ) :: rho_dofs(self%n_dofs) !< Coefficient vector of the charge distribution
!!$
!!$    sll_int32 :: box(3)
!!$    sll_real64 :: xi(3)
!!$    sll_int32  :: index3d(3)
!!$    sll_int32  :: index1d
!!$    sll_int32  :: i,j,k
!!$
!!$    
!!$    call convert_x_to_xbox( self, position, xi, box )
!!$    
!!$    self%spline_0 = 0.0_f64
!!$    call sll_s_spline_pp_horner_m_3d (self%spline_pp_0, self%spline_0, degree, xi)
!!$    ! USE THIS IF DEGREE IS NOT ALLWAYS SELF%SPLINE_DEGREE
!!$    !do j=1,3
!!$    !   if (degree(j) == self%spline_degree) then
!!$    !      call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_0, self%spline_0(:,j), self%spline_degree, xi(j))
!!$    !   else
!!$    !      call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_1, self%spline_0(:,j), self%spline_degree-1, xi(j))
!!$    !   end if
!!$    !end do
!!$    
!!$
!!$    ! Build scaling with marker_charge into spline along first dimension
!!$    self%spline_0(:,1) = self%spline_0(:,1)*marker_charge
!!$    ! 2d array combining first and second dimension
!!$    do j=1,degree(2)+1
!!$       do i=1,degree(1)+1
!!$          self%spline_2d(i,j) = self%spline_0(i,1)*self%spline_0(j,2)
!!$       end do
!!$    end do
!!$
!!$    
!!$    box = box-degree
!!$    do k=1,degree(3)+1       
!!$       index3d(3) = modulo(box(3)+k-2,self%n_grid(3))+1
!!$       do j=1,degree(2)+1
!!$          index3d(2) = modulo(box(2)+j-2,self%n_grid(2))+1
!!$          do i=1,degree(1)+1
!!$             index3d(1) = modulo(box(1)+i-2,self%n_grid(1))+1
!!$             index1d = convert_index_3d_to_1d( index3d, self%n_grid )
!!$             rho_dofs(index1d) = rho_dofs(index1d) + &
!!$                  (self%spline_2d(i,j) * &
!!$                  self%spline_0(k,3) )
!!$          end do
!!$       end do
!!$    end do
!!$
!!$
!!$  end subroutine add_charge_single_spline_pp_3d_feec


  pure function convert_index_2d_to_1d( index2d, n_grid ) result( index1d )
    sll_int32, intent( in ) :: index2d(2)
    sll_int32, intent( in ) :: n_grid(2)
    sll_int32 :: index1d
    
    index1d = index2d(1) + (index2d(2)-1)*n_grid(1)

    
  end function convert_index_2d_to_1d

  subroutine convert_x_to_xbox( self, position, xi, box )
    class( sll_t_particle_mesh_coupling_spline_2d_feec ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( out )    :: xi(self%dim) !< Position of the particle
    sll_int32,                                intent( out )    :: box(self%dim) !< Position of the particle

    xi = (position - self%domain(:,1)) * self%rdelta_x
    box = ceiling( xi )
    xi = xi - real(box-1, f64)
    
  end subroutine convert_x_to_xbox


  
  !> Add current for one particle 
  subroutine add_current_spline_2d_feec (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_2d_feec), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
    sll_real64, intent(in)    :: position_new(self%dim) !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(self%n_dofs*3) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(:) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(self%n_dofs) !< Coefficients of current expansion


    sll_int32 :: box(2), boxnew(2), boxold(2), local_size(2)
    sll_int32  :: degree
    sll_int32  :: index2d(2)
    sll_int32  :: index1d
    sll_int32  :: i,j
    sll_real64 :: xnew(2), xold(2)
    sll_int32  :: component

    degree = self%spline_degree
    component = 1
    
    call convert_x_to_xbox( self, position_old, xold, boxold )
    call convert_x_to_xbox( self, position_new, xnew, boxnew )

    local_size = abs(boxnew-boxold)+degree
    local_size(2) = degree+1
    
    do i=1, degree
       self%spline_1(i,1) =  sll_f_spline_pp_horner_1d(degree, self%spline_pp1d_1%poly_coeffs_fpa, xold(1), i) &
            * self%delta_x(1)
       self%spline_1(i,2) =  sll_f_spline_pp_horner_1d(degree, self%spline_pp1d_1%poly_coeffs_fpa, xnew(1), i) &
            * self%delta_x(1)
    end do

    if (position_old(1) < position_new(1) ) then
       self%j1d(local_size(component)-degree+1:local_size(component)) = self%spline_1(:,2)
       self%j1d(1:local_size(component)-degree) = self%delta_x(1)
       self%j1d(1:degree) = self%j1d(1:degree) - self%spline_1(:,1)
    else
       self%j1d(1:local_size(component)-degree) = -self%delta_x(1)
       self%j1d(local_size(component)-degree+1:local_size(component)) = -self%spline_1(:,1)
       self%j1d(1:degree) = self%j1d(1:degree) + self%spline_1(:,2)
    end if
    
    self%spline_0 = 0.0_f64
    do j=1,2
       if (j .ne. component ) then
          call sll_s_uniform_bsplines_eval_basis( self%spline_degree, xold(j), self%spline_0(:,j) )
       end if
    end do

    box = boxold-degree
    box(component) = min(boxnew(component), boxold(component))-degree+1
    do j=1,local_size(2)
       index2d(2) = modulo(box(2)+j-2,self%n_grid(2))+1
       do i=1,local_size(1)
          index2d(1) = modulo(box(1)+i-2,self%n_grid(1))+1
          index1d = convert_index_2d_to_1d( index2d, self%n_grid )
          j_dofs(index1d) = j_dofs(index1d) + &
               marker_charge * self%j1d( i ) * &
               self%spline_0(j,2)
       end do
    end do



  end subroutine add_current_spline_2d_feec



  !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting)
  subroutine add_current_update_v_primitive_component1_spline_2d_feec (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_2d_feec), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
    sll_real64, intent(in)    :: position_new !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(self%n_dofs*3) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(:) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(self%n_dofs) !< Coefficients of current expansion


    sll_int32 :: box(2), boxnew, boxold
    sll_real64 :: xi(2)
    sll_int32  :: index2d(2)
    sll_int32  :: index1d
    sll_int32  :: i,j
    sll_real64 :: xnew
    sll_real64 :: vt(2), vtt2, vtt3
    sll_int32  :: start1, start2
    sll_int32  :: component
    sll_int32  :: stride
    sll_int32 :: local_size
    sll_int32  :: degree

    component = 1
    start1 = 2*self%n_dofs
    start2 = self%n_dofs
    stride = 1
    degree = self%spline_degree
    
    call convert_x_to_xbox( self, position_old, xi, box )

    ! Convert position_new to xbox
    xnew = (position_new-self%domain(component,1)) * self%rdelta_x(component)
    boxnew = ceiling(xnew)
    xnew = xnew-real(boxnew-1, f64)
    boxold = box(component)

    !-- For current along x1
    local_size = abs(boxnew-boxold)+degree

    ! For j=component, we need the primitive
    do i=1, degree
       self%spline_1(i,1) =  sll_f_spline_pp_horner_1d(degree, &
            self%spline_pp1d_1%poly_coeffs_fpa, xi(component), i) &
            * self%delta_x(component)
       self%spline_1(i,2) =  sll_f_spline_pp_horner_1d(degree, &
            self%spline_pp1d_1%poly_coeffs_fpa, xnew, i) &
            * self%delta_x(component)
    end do

    !print*, boxold, boxnew
    !print*, position_old, position_new
    !print*, local_size, degree
    
    if (position_old(component) .le. position_new ) then
       self%j1d(local_size-degree+1:local_size) = self%spline_1(:,2)
       self%j1d(1:local_size-degree) = self%delta_x(component)
       self%j1d(1:degree) = self%j1d(1:degree) - self%spline_1(:,1)
    else
       self%j1d(1:local_size-degree) = -self%delta_x(component)
       self%j1d(local_size-degree+1:local_size) = -self%spline_1(:,1)
       self%j1d(1:degree) = self%j1d(1:degree) + self%spline_1(:,2)
    end if
    !----


    ! Achtung wir brauchen nun splines von beidem Grad
    self%spline_0 = 0.0_f64
    self%spline_1 = 0.0_f64
    do j=1,2
       if (j .ne. component ) then
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_0, self%spline_0(:,j), self%spline_degree, xi(j))
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_1, self%spline_1(:,j), self%spline_degree-1, xi(j))
       end if
    end do

    box(component) = box(component)-self%spline_degree+1
    box(2) = box(2) - self%spline_degree

    ! Define the range of the first component
    if (boxold<boxnew) then
       box(component) = boxold-self%spline_degree+1
    else
       box(component) = boxnew-self%spline_degree+1
    end if

      
    vtt2 = 0.0_f64
    vtt3 = 0.0_f64
    do j=1,self%spline_degree+1
       index2d(2) = modulo(box(2)+j-2,self%n_grid(2))+1
          
       vt = 0.0_f64
       
       do i=1,local_size
          index2d(1) = modulo(box(1)+i-2,self%n_grid(1))+1
          index1d = convert_index_2d_to_1d( index2d, self%n_grid )
          j_dofs(index1d) = j_dofs(index1d) + self%j1d(i) * &
               self%spline_0(j,2)*marker_charge
          
          vt(1) = vt(1) + bfield_dofs(start1+index1d)*self%j1d(i)
          vt(2) = vt(2) + bfield_dofs(start2+index1d)*self%j1d(i)
             
       end do
          
       if (j>1) then
          vtt2 = vtt2 + vt(1)*self%spline_1(j-1, 2)
       end if
       vtt3 = vtt3 - vt(2)*self%spline_0(j, 2)
          
    end do
    vi(2) = vi(2) - qoverm*vtt2
    vi(3) = vi(3) - qoverm*vtt3

  end subroutine add_current_update_v_primitive_component1_spline_2d_feec

  subroutine add_current_update_v_primitive_component2_spline_2d_feec (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_2d_feec), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
    sll_real64, intent(in)    :: position_new !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(self%n_dofs*3) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(:) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(self%n_dofs) !< Coefficients of current expansion



    sll_int32 :: box(2), boxnew, boxold
    sll_real64 :: xi(2)
    sll_int32  :: index2d(2)
    sll_int32  :: index1d
    sll_int32  :: i,j
    sll_real64 :: xnew
    sll_real64 :: vt(2), vtt2, vtt3
    sll_int32  :: start1, start2
    sll_int32  :: component, local_size
    sll_int32 :: stride
    sll_int32  :: degree
    
    component = 2
    start1 = 2*self%n_dofs !
    start2 = 0!
    stride = self%n_grid(1)
    degree = self%spline_degree
    
    call convert_x_to_xbox( self, position_old, xi, box )

    ! TODO: Umstellen von 1d function of ceiling statt floor
    ! Convert position_new to xbox
    xnew = (position_new-self%domain(component,1)) * self%rdelta_x(component)
    boxnew = ceiling(xnew)
    xnew = xnew-real(boxnew-1, f64)
    boxold = box(component)

    !-- For current along x2
    local_size = abs(boxnew-boxold)+degree

    ! For j=component, we need the primitive
    do i=1, degree
       self%spline_1(i,1) =  sll_f_spline_pp_horner_1d(degree, &
            self%spline_pp1d_1%poly_coeffs_fpa, xi(component), i) &
            * self%delta_x(component)
       self%spline_1(i,2) =  sll_f_spline_pp_horner_1d(degree, &
            self%spline_pp1d_1%poly_coeffs_fpa, xnew, i) &
            * self%delta_x(component)
    end do

    if (position_old(component) .le. position_new ) then
       self%j1d(local_size-degree+1:local_size) = self%spline_1(:,2)
       self%j1d(1:local_size-degree) = self%delta_x(component)
       self%j1d(1:degree) = self%j1d(1:degree) - self%spline_1(:,1)
    else
       self%j1d(1:local_size-degree) = -self%delta_x(component)
       self%j1d(local_size-degree+1:local_size) = -self%spline_1(:,1)
       self%j1d(1:degree) = self%j1d(1:degree) + self%spline_1(:,2)
    end if
    !----


    ! Achtung wir brauchen nun splines von beidem Grad
    self%spline_0 = 0.0_f64
    self%spline_1 = 0.0_f64
    do j=1,2
       if (j .ne. component ) then
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_0, self%spline_0(:,j), self%spline_degree, xi(j))
          call sll_s_spline_pp_horner_m_1d(self%spline_pp1d_1, self%spline_1(:,j), self%spline_degree-1, xi(j))
       end if
    end do

    box(1) = box(1) - self%spline_degree

    ! Define the range of the first component
    if (boxold<boxnew) then
       box(component) = boxold-self%spline_degree+1
    else
       box(component) = boxnew-self%spline_degree+1
    end if
      
    vtt2 = 0.0_f64
    vtt3 = 0.0_f64
    do j=1,self%spline_degree+1
       index2d(1) = modulo(box(1)+j-2,self%n_grid(1))+1
       vt = 0.0_f64
       do i=1,local_size
          index2d(2) = modulo(box(2)+i-2,self%n_grid(2))+1
          index1d = convert_index_2d_to_1d( index2d, self%n_grid )
          j_dofs(index1d) = j_dofs(index1d) + self%j1d(i) * &
               self%spline_0(j,1) * marker_charge
             
          vt(1) = vt(1) + bfield_dofs(start1+index1d)*self%j1d(i)
          vt(2) = vt(2) + bfield_dofs(start2+index1d)*self%j1d(i)
       end do

       if (j>1) then
          vtt2 = vtt2 + vt(1)*self%spline_1(j-1, 1)
       end if
       vtt3 = vtt3 + vt(2)*self%spline_0(j, 1)
          
    end do
    vi(1) = vi(1) + qoverm*vtt2
    vi(3) = vi(3) - qoverm*vtt3
    
  end subroutine add_current_update_v_primitive_component2_spline_2d_feec
  

  !---------------------------------------------------------------------------!
 !> Evaluate field at at position \a position
  subroutine evaluate_field_single_spline_2d_feec(self, position, degree, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_2d_feec), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_int32 ,                              intent( in )    :: degree(self%dim) !< Spline degree of the various components
    sll_real64,                              intent( in )    :: field_dofs(self%n_dofs) !< Coefficient vector for the field DoFs
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the electric fields at given position

    sll_int32 :: i,j
    sll_int32 :: box(2)
    sll_int32 :: index2d(2)
    sll_int32 :: index1d
    sll_real64 :: xi(2)
    
    ! TODO: Optimize by sum factorization

    call convert_x_to_xbox( self, position, xi, box )

    do j=1,2
       call sll_s_uniform_bsplines_eval_basis( degree(j), xi(j), self%spline_val(1:degree(j)+1,j) )
    end do

    field_value = 0.0_f64
    box = box-degree
    do j=1,degree(2)+1
       index2d(2) = modulo(box(2)+j-2,self%n_grid(2))+1
       do i=1,degree(1)+1
          index2d(1) = modulo(box(1)+i-2,self%n_grid(1))+1
          index1d = convert_index_2d_to_1d( index2d, self%n_grid )
          field_value = field_value + &
               field_dofs(index1d) * &
               self%spline_val(i,1) * self%spline_val(j,2)
       end do
    end do

  end subroutine evaluate_field_single_spline_2d_feec
  
 !---------------------------------------------------------------------------!
 !> Evaluate field at at position \a position
  subroutine evaluate_field_single_spline_pp_2d_feec(self, position, degree, field_dofs_pp, field_value)
    class (sll_t_particle_mesh_coupling_spline_2d_feec), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_int32 ,                              intent( in )    :: degree(self%dim) !< Spline degree of the various components
    sll_real64,                              intent( in )    :: field_dofs_pp(:,:) !< Coefficient vector for the field DoFs
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the electric fields at given position

    sll_int32 :: box(3)
    sll_real64 :: xi(3)

    call convert_x_to_xbox( self, position, xi, box )
   
    field_value = sll_f_spline_pp_horner_2d(degree, field_dofs_pp, xi, box,self%n_grid)

  end subroutine evaluate_field_single_spline_pp_2d_feec



  !---------------------------------------------------------------------------!
  !> Evaluate several fields at position \a position
  subroutine evaluate_multiple_spline_2d_feec(self, position, components, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_2d_feec), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_int32,                               intent(in)      :: components(:) !< Components of field_dofs that shall be updated
    sll_real64,                              intent( in )    :: field_dofs(:,:) !< Coefficient vector for the field DoFs
    sll_real64,                              intent(out)     :: field_value(:) !< Value(s) of the electric fields at given position
    


  end subroutine evaluate_multiple_spline_2d_feec

  !> Constructor
  subroutine init_spline_2d_feec ( self, n_grid, domain, no_particles, spline_degree )
    class (sll_t_particle_mesh_coupling_spline_2d_feec), intent( out ) :: self !< Kernel smoother object 
    sll_int32,                               intent(in) :: n_grid(2) !< number of DoFs (spline coefficients)
    sll_real64,                              intent(in) :: domain(2,2) !< x_min and x_max of the domain
    sll_int32,                               intent(in) :: no_particles !< number of particles
    sll_int32,                               intent(in) :: spline_degree !< Degree of smoothing kernel spline

    self%dim = 2

    ! Store grid information
    self%domain = domain
    self%n_grid = n_grid
    self%n_dofs = product(n_grid)
    self%delta_x = (self%domain(:,2)-self%domain(:,1))/real(n_grid, f64)
    self%rdelta_x = 1.0_f64/self%delta_x
    
    ! Store basis function information
    self%no_particles = no_particles

    ! Initialize information on the spline
    self%spline_degree = spline_degree
    self%n_span = spline_degree + 1
    
    
    self%n_quad_points = (self%spline_degree+2)/2
    allocate( self%quad_xw(2,self%n_quad_points) )
    ! normalized Gauss Legendre points and weights
    self%quad_xw = sll_f_gauss_legendre_points_and_weights(self%n_quad_points)

    
    allocate( self%spline_val(self%n_span,2) )
    allocate( self%spline_val_more(self%n_span,2) )
    allocate( self%spline_0(self%n_span,2) )
    allocate( self%spline_1(self%n_span-1,2) )
    allocate( self%j1d( maxval(self%n_grid) ))
    allocate( self%spline_2d(self%n_span, self%n_span) )


!!$    call sll_s_spline_pp_init_3d(self%spline_pp_0, [spline_degree,spline_degree,spline_degree], n_grid)
!!$    call sll_s_spline_pp_init_3d(self%spline_pp_11, [spline_degree-1,spline_degree,spline_degree], n_grid)
!!$    call sll_s_spline_pp_init_3d(self%spline_pp_12, [spline_degree,spline_degree-1,spline_degree], n_grid)
!!$    call sll_s_spline_pp_init_3d(self%spline_pp_13, [spline_degree,spline_degree,spline_degree-1], n_grid)
!!$    call sll_s_spline_pp_init_3d(self%spline_pp_21, [spline_degree,spline_degree-1,spline_degree-1], n_grid)
!!$    call sll_s_spline_pp_init_3d(self%spline_pp_22, [spline_degree-1,spline_degree,spline_degree-1], n_grid)
!!$    call sll_s_spline_pp_init_3d(self%spline_pp_23, [spline_degree-1,spline_degree-1,spline_degree], n_grid)

    call sll_s_spline_pp_init_1d( self%spline_pp1d_0, spline_degree, n_grid(1) )
    call sll_s_spline_pp_init_1d( self%spline_pp1d_1, spline_degree-1, n_grid(1) )
    
  end subroutine init_spline_2d_feec

  
  !> Destructor
  subroutine free_spline_2d_feec(self)
    class (sll_t_particle_mesh_coupling_spline_2d_feec), intent( inout ) :: self !< Kernel smoother object 

    deallocate( self%quad_xw)
    deallocate( self%spline_val)
    deallocate( self%spline_val_more)
    deallocate( self%spline_0)
    deallocate( self%spline_1)
    deallocate( self%j1d)
    
    !call sll_s_spline_pp_free_3d(self%spline_pp_0)

  end subroutine free_spline_2d_feec


end module sll_m_particle_mesh_coupling_spline_2d_feec
