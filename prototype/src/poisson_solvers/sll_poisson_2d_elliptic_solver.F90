!**************************************************************
!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
!> @author
!> Adnane Hamiaz (hamiaz@math.unistra.fr)
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
!**************************************************************



module sll_module_poisson_2d_elliptic_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!use sll_boundary_condition_descriptors
  use sll_module_poisson_2d_base
  use sll_general_coordinate_elliptic_solver_module
  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative
implicit none

  type,extends(sll_poisson_2d_base) :: poisson_2d_elliptic_solver      
    type(general_coordinate_elliptic_solver), pointer      :: elliptic_solver
    class(sll_scalar_field_2d_discrete_alt), pointer        :: phi_field
    class(sll_scalar_field_2d_base), pointer                :: rho_field
    class(sll_scalar_field_2d_base), pointer                :: a11_field
    class(sll_scalar_field_2d_base), pointer                :: a12_field
    class(sll_scalar_field_2d_base), pointer                :: a21_field
    class(sll_scalar_field_2d_base), pointer                :: a22_field
    class(sll_scalar_field_2d_base), pointer                :: b1_field
    class(sll_scalar_field_2d_base), pointer                :: b2_field
    class(sll_scalar_field_2d_base), pointer                :: c_field
    type(arb_deg_2d_interpolator)                           :: interp_rho
    type(arb_deg_2d_interpolator)                           :: interp_phi
    type(arb_deg_2d_interpolator)                           :: interp_a11
    type(arb_deg_2d_interpolator)                           :: interp_a12
    type(arb_deg_2d_interpolator)                           :: interp_a21
    type(arb_deg_2d_interpolator)                           :: interp_a22
    type(arb_deg_2d_interpolator)                           :: interp_b1
    type(arb_deg_2d_interpolator)                           :: interp_b2
    type(arb_deg_2d_interpolator)                           :: interp_c
  contains
    procedure, pass(poisson) :: initialize => &
      initialize_poisson_2d_elliptic_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_2d_elliptic_solver
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_2d_elliptic_solver
  end type poisson_2d_elliptic_solver

contains




 subroutine initialize_poisson_2d_elliptic_solver(poisson, transf,&
   spline_degree_eta1, &
   spline_degree_eta2, &
   num_cells_eta1, &
   num_cells_eta2, &
   quadrature_type1, &
   quadrature_type2, &
   bc_eta1_left, &
   bc_eta1_right, &
   bc_eta2_left, &
   bc_eta2_right, &
   eta1_min, &
   eta1_max, &
   eta2_min, &
   eta2_max, &
   a11_values, & 
   a12_values, & 
   a21_values, & 
   a22_values, & 
   b1_values,&
   b2_values,&
   c_values ) 
   
   class(poisson_2d_elliptic_solver),        target  :: poisson
   class(sll_coordinate_transformation_2d_base), pointer :: transf
   sll_int32, intent(in)  :: spline_degree_eta1
   sll_int32, intent(in)  :: spline_degree_eta2
   sll_int32, intent(in)  :: num_cells_eta1
   sll_int32, intent(in)  :: num_cells_eta2
   sll_int32, intent(in)  :: bc_eta1_left
   sll_int32, intent(in)  :: bc_eta1_right
   sll_int32, intent(in)  :: bc_eta2_left
   sll_int32, intent(in)  :: bc_eta2_right
   sll_int32, intent(in)  :: quadrature_type1
   sll_int32, intent(in)  :: quadrature_type2
   sll_real64, intent(in) :: eta1_min
   sll_real64, intent(in) :: eta1_max
   sll_real64, intent(in) :: eta2_min
   sll_real64, intent(in) :: eta2_max
   sll_real64, dimension(:,:), pointer :: phi_values
   sll_real64, dimension(:,:), pointer :: rho_values
   sll_real64, dimension(:,:)          :: a11_values
   sll_real64, dimension(:,:)          :: a12_values
   sll_real64, dimension(:,:)          :: a21_values
   sll_real64, dimension(:,:)          :: a22_values
   sll_real64, dimension(:,:)          :: b1_values
   sll_real64, dimension(:,:)          :: b2_values
   sll_real64, dimension(:,:)          :: c_values
   sll_int32 :: np_eta1
   sll_int32 :: np_eta2
   sll_int32 :: ierr


    np_eta1 = num_cells_eta1 + 1
    np_eta2 = num_cells_eta2 + 1
    
    call poisson%interp_phi%initialize( &
         np_eta1, &
         np_eta2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right,&
         spline_degree_eta1, &
         spline_degree_eta2)
   
   call poisson%interp_rho%initialize( &
         np_eta1, &
         np_eta2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right,&
         spline_degree_eta1, &
         spline_degree_eta2)
         
    call poisson%interp_a11%initialize( &
         np_eta1, &
         np_eta2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right,&
         spline_degree_eta1, &
         spline_degree_eta2)     
          
     call poisson%interp_a12%initialize( &
         np_eta1, &
         np_eta2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right,&
         spline_degree_eta1, &
         spline_degree_eta2)    
           
    call poisson%interp_a21%initialize( &
         np_eta1, &
         np_eta2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right,&
         spline_degree_eta1, &
         spline_degree_eta2) 
          
    call poisson%interp_a22%initialize( &
         np_eta1, &
         np_eta2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right,&
         spline_degree_eta1, &
         spline_degree_eta2)   

    call poisson%interp_a21%initialize( &
         np_eta1, &
         np_eta2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right,&
         spline_degree_eta1, &
         spline_degree_eta2) 
          
    call poisson%interp_b1%initialize( &
         np_eta1, &
         np_eta2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right,&
         spline_degree_eta1, &
         spline_degree_eta2)
      
    call poisson%interp_b2%initialize( &
         np_eta1, &
         np_eta2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right,&
         spline_degree_eta1, &
         spline_degree_eta2)   
                                     
    poisson%a11_field => new_scalar_field_2d_discrete_alt( &
         "a11_check", &
         poisson%interp_a11, &
         transf, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%a11_field%set_field_data( a11_values )
    call poisson%a11_field%update_interpolation_coefficients( )  
     
    poisson%a12_field => new_scalar_field_2d_discrete_alt( &
         "a12_check", &
         poisson%interp_a12, &
         transf, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%a12_field%set_field_data( a12_values )
    call poisson%a12_field%update_interpolation_coefficients( ) 
    
    poisson%a21_field => new_scalar_field_2d_discrete_alt( &
         "a21_check", &
         poisson%interp_a21, &
         transf, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%a21_field%set_field_data( a21_values )
    call poisson%a21_field%update_interpolation_coefficients( ) 
    
    poisson%a22_field => new_scalar_field_2d_discrete_alt( &
         "a22_check", &
         poisson%interp_a22, &
         transf, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%a22_field%set_field_data( a22_values )
    call poisson%a22_field%update_interpolation_coefficients( )

    poisson%b1_field => new_scalar_field_2d_discrete_alt( &
         "b1_check", &
         poisson%interp_b1, &
         transf, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%b1_field%set_field_data( b1_values )
    call poisson%b1_field%update_interpolation_coefficients( ) 
    
    poisson%b2_field => new_scalar_field_2d_discrete_alt( &
         "b2_check", &
         poisson%interp_b2, &
         transf, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%b2_field%set_field_data( b2_values )
    call poisson%b2_field%update_interpolation_coefficients( )
    
    poisson%c_field => new_scalar_field_2d_discrete_alt( &
         "c_check", &
         poisson%interp_c, &
         transf, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%c_field%set_field_data( c_values )
    call poisson%c_field%update_interpolation_coefficients( )
    
    SLL_ALLOCATE(phi_values(np_eta1,np_eta2),ierr)
    SLL_ALLOCATE(rho_values(np_eta1,np_eta2),ierr)
    phi_values(:,:) = 0.0_f64
    rho_values(:,:) = 0.0_f64
  
    poisson%phi_field => new_scalar_field_2d_discrete_alt( &
         "phi_check", &
         poisson%interp_phi, &
         transf, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%phi_field%set_field_data( phi_values )
    call poisson%phi_field%update_interpolation_coefficients( )  
    
     poisson%rho_field => new_scalar_field_2d_discrete_alt( &
         "rho_check", &
         poisson%interp_rho, &
         transf, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%rho_field%set_field_data( rho_values )
    call poisson%rho_field%update_interpolation_coefficients( )
       
    poisson%elliptic_solver => new_general_elliptic_solver( &
        spline_degree_eta1, &
        spline_degree_eta2, &
        num_cells_eta1, &
        num_cells_eta2, &
        quadrature_type1, &
        quadrature_type2, &
        bc_eta1_left, &
        bc_eta1_right, &
        bc_eta2_left, &
        bc_eta2_right, &
        eta1_min, &
        eta1_max, &
        eta2_min, &
        eta2_max )
        
    ! compute matrix the field
    print *,'Compute matrix the field'
    call factorize_mat_es(&
        poisson%elliptic_solver, &
        poisson%a11_field, &
        poisson%a12_field,&
        poisson%a21_field,&
        poisson%a22_field,&
        poisson%b1_field,&
        poisson%b2_field,&
        poisson%c_field)    
        
 end subroutine initialize_poisson_2d_elliptic_solver
 
 
 function new_poisson_2d_elliptic_solver( &
   transf,&
   spline_degree_eta1, &
   spline_degree_eta2, &
   num_cells_eta1, &
   num_cells_eta2, &
   quadrature_type1, &
   quadrature_type2, &
   bc_eta1_left, &
   bc_eta1_right, &
   bc_eta2_left, &
   bc_eta2_right, &
   eta1_min, &
   eta1_max, &
   eta2_min, &
   eta2_max, &
   a11_values, & 
   a12_values, & 
   a21_values, & 
   a22_values, & 
   b1_values,&
   b2_values,&
   c_values ) &
   result(poisson)
   
   class(poisson_2d_elliptic_solver),        pointer  :: poisson
   class(sll_coordinate_transformation_2d_base), pointer :: transf
   sll_int32, intent(in)  :: spline_degree_eta1
   sll_int32, intent(in)  :: spline_degree_eta2
   sll_int32, intent(in)  :: num_cells_eta1
   sll_int32, intent(in)  :: num_cells_eta2
   sll_int32, intent(in)  :: bc_eta1_left
   sll_int32, intent(in)  :: bc_eta1_right
   sll_int32, intent(in)  :: bc_eta2_left
   sll_int32, intent(in)  :: bc_eta2_right
   sll_int32, intent(in)  :: quadrature_type1
   sll_int32, intent(in)  :: quadrature_type2
   sll_real64, intent(in) :: eta1_min
   sll_real64, intent(in) :: eta1_max
   sll_real64, intent(in) :: eta2_min
   sll_real64, intent(in) :: eta2_max
   sll_real64, dimension(:,:), pointer :: phi_values
   sll_real64, dimension(:,:), pointer :: rho_values
   sll_real64, dimension(:,:)          :: a11_values
   sll_real64, dimension(:,:)          :: a12_values
   sll_real64, dimension(:,:)          :: a21_values
   sll_real64, dimension(:,:)          :: a22_values
   sll_real64, dimension(:,:)          :: b1_values
   sll_real64, dimension(:,:)          :: b2_values
   sll_real64, dimension(:,:)          :: c_values
   sll_int32 :: np_eta1
   sll_int32 :: np_eta2
   sll_int32 :: ierr
   
   SLL_ALLOCATE(poisson,ierr)
   call initialize_poisson_2d_elliptic_solver(poisson, transf,&
   spline_degree_eta1, &
   spline_degree_eta2, &
   num_cells_eta1, &
   num_cells_eta2, &
   quadrature_type1, &
   quadrature_type2, &
   bc_eta1_left, &
   bc_eta1_right, &
   bc_eta2_left, &
   bc_eta2_right, &
   eta1_min, &
   eta1_max, &
   eta2_min, &
   eta2_max, &
   a11_values, & 
   a12_values, & 
   a21_values, & 
   a22_values, & 
   b1_values,&
   b2_values,&
   c_values )  
  end function new_poisson_2d_elliptic_solver
  
  subroutine compute_phi_from_rho_2d_elliptic_solver(poisson,phi,rho )
    ! input variables 
    class(poisson_2d_elliptic_solver), target   :: poisson
    sll_real64, dimension(:,:), intent(in)     :: rho
    ! output variables
    sll_real64, dimension(:,:), intent(out)  :: phi
    ! local variables
    !class(general_coordinate_elliptic_solver), pointer   :: elliptic_solver
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32  :: i,j
    sll_int32  :: nc_eta1
    sll_int32  :: nc_eta2
     
    ! The supposition is that all fields use the same logical mesh
    !elliptic_solver => poisson%elliptic_solver
   
    delta1    = poisson%elliptic_solver%delta_eta1
    delta2    = poisson%elliptic_solver%delta_eta2
    eta1_min  = poisson%elliptic_solver%eta1_min
    eta2_min  = poisson%elliptic_solver%eta2_min
    nc_eta1   = poisson%elliptic_solver%num_cells1 !+ 1    
    nc_eta2   = poisson%elliptic_solver%num_cells2 !+ 1
    
    call poisson%rho_field%set_field_data(rho)
    call poisson%rho_field%update_interpolation_coefficients( )
            
    call solve_general_coordinates_elliptic_eq(&
       poisson%elliptic_solver,&
       poisson%rho_field,&
       poisson%phi_field)
   
   do j=1,nc_eta2+1
        do i=1,nc_eta1+1
           phi(i,j) = poisson%phi_field%value_at_indices(i,j)
        end do
     end do
     
  end subroutine compute_phi_from_rho_2d_elliptic_solver

    ! solves E = -\nabla Phi in 2d
!    subroutine compute_E_from_phi_2d_fft( poisson, phi, E1, E2 )
!      class(poisson_2d_fft_solver) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine compute_E_from_phi_2d_fft

    ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
    subroutine compute_E_from_rho_2d_elliptic_solver( poisson, E1, E2, rho )
      class(poisson_2d_elliptic_solver) :: poisson
      sll_real64,dimension(:,:),intent(in) :: rho
      sll_real64,dimension(:,:),intent(out) :: E1
      sll_real64,dimension(:,:),intent(out) :: E2
      
      E1 = 0.0_f64
      E2 = 0.0_f64
      print *,'#compute_E_from_rho_2d_elliptic_solver'      
      print *,'#not implemented for the moment'
      stop
      
      !call solve( poisson%poiss, E1, E2, rho)
      
    end subroutine compute_E_from_rho_2d_elliptic_solver
    
 end module sll_module_poisson_2d_elliptic_solver
  
