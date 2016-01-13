#ifndef DOXYGEN_SHOULD_SKIP_THIS
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



!> @ingroup poisson_solvers
module sll_m_poisson_2d_elliptic_solver
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use iso_fortran_env, only: &
    output_unit

  use sll_m_arbitrary_degree_spline_interpolator_2d, only: &
    sll_f_new_arbitrary_degree_spline_interp2d

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_cubic_spline_interpolator_2d, only: &
    sll_f_new_cubic_spline_interpolator_2d

  use sll_m_general_coordinate_elliptic_solver, only: &
    sll_s_factorize_mat_es_prototype, &
    sll_t_general_coordinate_elliptic_solver, &
    sll_f_new_general_elliptic_solver_prototype, &
    sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype, &
    sll_s_solve_general_coordinates_elliptic_eq_prototype

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base, &
    sll_f_function_of_position

  use sll_m_scalar_field_2d, only: &
    sll_f_new_scalar_field_2d_discrete

  use sll_m_scalar_field_2d_base, only: &
    sll_c_scalar_field_2d_base

  implicit none

  public :: &
    sll_f_new_poisson_2d_elliptic_solver

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: SLL_NO_SOLVE_ELLIPTIC_SOLVER = 0 
  sll_int32, parameter :: SLL_SOLVE_ELLIPTIC_SOLVER = 1 
  sll_int32, parameter :: SLL_DO_NOTHING_ELLIPTIC_SOLVER = 2 

  type,extends(sll_c_poisson_2d_base) :: poisson_2d_elliptic_solver      
    type(sll_t_general_coordinate_elliptic_solver), pointer      :: elliptic_solver
    !class(sll_t_scalar_field_2d_discrete), pointer        :: phi_field
    class(sll_c_scalar_field_2d_base), pointer                :: rho_field
    class(sll_c_scalar_field_2d_base), pointer                :: a11_field
    class(sll_c_scalar_field_2d_base), pointer                :: a12_field
    class(sll_c_scalar_field_2d_base), pointer                :: a21_field
    class(sll_c_scalar_field_2d_base), pointer                :: a22_field
    class(sll_c_scalar_field_2d_base), pointer                :: b1_field
    class(sll_c_scalar_field_2d_base), pointer                :: b2_field
    class(sll_c_scalar_field_2d_base), pointer                :: c_field
    class(sll_c_interpolator_2d), pointer                :: interp_rho
    class(sll_c_interpolator_2d), pointer                :: interp_a11
    class(sll_c_interpolator_2d), pointer                :: interp_a12
    class(sll_c_interpolator_2d), pointer                :: interp_a21
    class(sll_c_interpolator_2d), pointer                :: interp_a22
    class(sll_c_interpolator_2d), pointer                :: interp_b1
    class(sll_c_interpolator_2d), pointer                :: interp_b2
    class(sll_c_interpolator_2d), pointer                :: interp_c
    !type(sll_t_arbitrary_degree_spline_interpolator_2d)                           :: interp_rho
    !type(sll_t_arbitrary_degree_spline_interpolator_2d)                           :: interp_phi
!    type(sll_t_arbitrary_degree_spline_interpolator_2d)                           :: interp_a11
!    type(sll_t_arbitrary_degree_spline_interpolator_2d)                           :: interp_a12
!    type(sll_t_arbitrary_degree_spline_interpolator_2d)                           :: interp_a21
!    type(sll_t_arbitrary_degree_spline_interpolator_2d)                           :: interp_a22
!    type(sll_t_arbitrary_degree_spline_interpolator_2d)                           :: interp_b1
!    type(sll_t_arbitrary_degree_spline_interpolator_2d)                           :: interp_b2
!    type(sll_t_arbitrary_degree_spline_interpolator_2d)                           :: interp_c
    sll_int32 :: control
    logical :: precompute_rhs
  contains
    procedure, pass(poisson) :: initialize => &
      initialize_poisson_2d_elliptic_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_2d_elliptic_solver
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_2d_elliptic_solver

    procedure :: l2norm_squared => l2norm_squarred_2d_elliptic_solver
    procedure :: compute_rhs_from_function => compute_rhs_from_function_2d_elliptic_solver
    procedure :: delete => delete_poisson_2d_elliptic_solver

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
   bc_interp2d_eta1, &
   bc_interp2d_eta2, &
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
   c_values, &
   interp_rho, &
   interp_rho_case, &
   control, &
   precompute_rhs, &
   with_constraint, &
   zero_mean ) 
   
   class(poisson_2d_elliptic_solver),        target  :: poisson
   class(sll_c_coordinate_transformation_2d_base), pointer :: transf
   sll_int32, intent(in)  :: spline_degree_eta1
   sll_int32, intent(in)  :: spline_degree_eta2
   sll_int32, intent(in)  :: num_cells_eta1
   sll_int32, intent(in)  :: num_cells_eta2
   sll_int32, intent(in)  :: bc_eta1_left
   sll_int32, intent(in)  :: bc_eta1_right
   sll_int32, intent(in)  :: bc_eta2_left
   sll_int32, intent(in)  :: bc_eta2_right
   sll_int32, intent(in)  :: bc_interp2d_eta1
   sll_int32, intent(in)  :: bc_interp2d_eta2
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
   class(sll_c_interpolator_2d), pointer, optional :: interp_rho
   character(len=256), optional  :: interp_rho_case
   sll_int32, intent(in), optional :: control
   logical, intent(in), optional :: precompute_rhs
   logical, intent(in), optional :: with_constraint
   logical, intent(in), optional :: zero_mean
   sll_int32 :: np_eta1
   sll_int32 :: np_eta2
   sll_int32 :: ierr
   logical :: use_cubic_splines
   
   sll_int32 :: bc_eta1_left_interp
   
   if(present(control))then
     poisson%control = control
   else
     poisson%control =  SLL_SOLVE_ELLIPTIC_SOLVER 
   endif
   
   if(present(precompute_rhs))then
     poisson%precompute_rhs = precompute_rhs
   else
     poisson%precompute_rhs =  .false. 
   endif
   
   if(poisson%precompute_rhs)then
     use_cubic_splines = .true.
   else
     use_cubic_splines = .false.    
   endif
   
   if(poisson%control==SLL_DO_NOTHING_ELLIPTIC_SOLVER)then
     return
   endif
   bc_eta1_left_interp = bc_eta1_left
   
   if(bc_eta1_left==sll_p_neumann)then
     bc_eta1_left_interp = sll_p_dirichlet
   endif

    np_eta1 = num_cells_eta1 + 1
    np_eta2 = num_cells_eta2 + 1
    
!    call poisson%interp_phi%initialize( &
!         np_eta1, &
!         np_eta2, &
!         eta1_min, &
!         eta1_max, &
!         eta2_min, &
!         eta2_max, &
!         bc_eta1_left_interp, &
!         bc_eta1_right, &
!         bc_eta2_left, &
!         bc_eta2_right,&
!         spline_degree_eta1, &
!         spline_degree_eta2)

   if(present(interp_rho))then
     poisson%interp_rho => interp_rho
   else
     if(present(interp_rho_case))then
       select case(interp_rho_case)
         case ("SLL_CUBIC_SPLINES")
           poisson%interp_rho => sll_f_new_cubic_spline_interpolator_2d( &
             np_eta1, &
             np_eta2, &
             eta1_min, &
             eta1_max, &
             eta2_min, &
             eta2_max, &
             bc_interp2d_eta1, &
             bc_interp2d_eta2)
         case ("SLL_ARBITRARY_DEGREE_SPLINES")
           poisson%interp_rho => sll_f_new_arbitrary_degree_spline_interp2d( &
             np_eta1, &
             np_eta2, &
             eta1_min, &
             eta1_max, &
             eta2_min, &
             eta2_max, &
             bc_eta1_left_interp, &
             bc_eta1_right, &
             bc_eta2_left, &
             bc_eta2_right,&
             spline_degree_eta1, &
             spline_degree_eta2)
         case default
           print *,'interp_rho_case=',interp_rho_case
           flush( output_unit )
           SLL_ERROR('initialize_poisson_2d_elliptic_solver','bad interp_rho_case')    
       end select
     else
       poisson%interp_rho => sll_f_new_cubic_spline_interpolator_2d( &
         np_eta1, &
         np_eta2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_interp2d_eta1, &
         bc_interp2d_eta2)
     endif                         
   endif

!   if(local_precompute_rhs)then
!     poisson%interp_rho => sll_f_new_arbitrary_degree_spline_interp2d( &
!         np_eta1, &
!         np_eta2, &
!         eta1_min, &
!         eta1_max, &
!         eta2_min, &
!         eta2_max, &
!         bc_eta1_left_interp, &
!         bc_eta1_right, &
!         bc_eta2_left, &
!         bc_eta2_right,&
!         spline_degree_eta1, &
!         spline_degree_eta2)
!   else
!     poisson%interp_rho => sll_f_new_cubic_spline_interpolator_2d( &
!          np_eta1, &
!          np_eta2, &
!          eta1_min, &
!          eta1_max, &
!          eta2_min, &
!          eta2_max, &
!          bc_interp2d_eta1, &
!          bc_interp2d_eta2)                
!   endif



    
  poisson%interp_a11 => sll_f_new_cubic_spline_interpolator_2d( &
    np_eta1, &
    np_eta2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    bc_interp2d_eta1, &
    bc_interp2d_eta2)                


  poisson%interp_a12 => sll_f_new_cubic_spline_interpolator_2d( &
    np_eta1, &
    np_eta2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    bc_interp2d_eta1, &
    bc_interp2d_eta2)                

  poisson%interp_a21 => sll_f_new_cubic_spline_interpolator_2d( &
    np_eta1, &
    np_eta2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    bc_interp2d_eta1, &
    bc_interp2d_eta2)                

  poisson%interp_a22 => sll_f_new_cubic_spline_interpolator_2d( &
    np_eta1, &
    np_eta2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    bc_interp2d_eta1, &
    bc_interp2d_eta2)                


  poisson%interp_b1 => sll_f_new_cubic_spline_interpolator_2d( &
    np_eta1, &
    np_eta2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    bc_interp2d_eta1, &
    bc_interp2d_eta2)                

  poisson%interp_b2 => sll_f_new_cubic_spline_interpolator_2d( &
    np_eta1, &
    np_eta2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    bc_interp2d_eta1, &
    bc_interp2d_eta2)                

  poisson%interp_c => sll_f_new_cubic_spline_interpolator_2d( &
    np_eta1, &
    np_eta2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    bc_interp2d_eta1, &
    bc_interp2d_eta2)                





   
       
!    call poisson%interp_a11%initialize( &
!         np_eta1, &
!         np_eta2, &
!         eta1_min, &
!         eta1_max, &
!         eta2_min, &
!         eta2_max, &
!         bc_eta1_left_interp, &
!         bc_eta1_right, &
!         bc_eta2_left, &
!         bc_eta2_right,&
!         spline_degree_eta1, &
!         spline_degree_eta2)     
!          
!     call poisson%interp_a12%initialize( &
!         np_eta1, &
!         np_eta2, &
!         eta1_min, &
!         eta1_max, &
!         eta2_min, &
!         eta2_max, &
!         bc_eta1_left_interp, &
!         bc_eta1_right, &
!         bc_eta2_left, &
!         bc_eta2_right,&
!         spline_degree_eta1, &
!         spline_degree_eta2)    
!           
!    call poisson%interp_a21%initialize( &
!         np_eta1, &
!         np_eta2, &
!         eta1_min, &
!         eta1_max, &
!         eta2_min, &
!         eta2_max, &
!         bc_eta1_left_interp, &
!         bc_eta1_right, &
!         bc_eta2_left, &
!         bc_eta2_right,&
!         spline_degree_eta1, &
!         spline_degree_eta2) 
!          
!    call poisson%interp_a22%initialize( &
!         np_eta1, &
!         np_eta2, &
!         eta1_min, &
!         eta1_max, &
!         eta2_min, &
!         eta2_max, &
!         bc_eta1_left_interp, &
!         bc_eta1_right, &
!         bc_eta2_left, &
!         bc_eta2_right,&
!         spline_degree_eta1, &
!         spline_degree_eta2)   
!
!    call poisson%interp_a21%initialize( &
!         np_eta1, &
!         np_eta2, &
!         eta1_min, &
!         eta1_max, &
!         eta2_min, &
!         eta2_max, &
!         bc_eta1_left_interp, &
!         bc_eta1_right, &
!         bc_eta2_left, &
!         bc_eta2_right,&
!         spline_degree_eta1, &
!         spline_degree_eta2) 
!     
!    call poisson%interp_c%initialize( &
!         np_eta1, &
!         np_eta2, &
!         eta1_min, &
!         eta1_max, &
!         eta2_min, &
!         eta2_max, &
!         bc_eta1_left_interp, &
!         bc_eta1_right, &
!         bc_eta2_left, &
!         bc_eta2_right,&
!         spline_degree_eta1, &
!         spline_degree_eta2)     
!          
!    call poisson%interp_b1%initialize( &
!         np_eta1, &
!         np_eta2, &
!         eta1_min, &
!         eta1_max, &
!         eta2_min, &
!         eta2_max, &
!         bc_eta1_left_interp, &
!         bc_eta1_right, &
!         bc_eta2_left, &
!         bc_eta2_right,&
!         spline_degree_eta1, &
!         spline_degree_eta2)
!      
!    call poisson%interp_b2%initialize( &
!         np_eta1, &
!         np_eta2, &
!         eta1_min, &
!         eta1_max, &
!         eta2_min, &
!         eta2_max, &
!         bc_eta1_left_interp, &
!         bc_eta1_right, &
!         bc_eta2_left, &
!         bc_eta2_right,&
!         spline_degree_eta1, &
!         spline_degree_eta2)   
                                   
    poisson%a11_field => sll_f_new_scalar_field_2d_discrete( &
         "a11_check", &
         poisson%interp_a11, &
         transf, &
         bc_eta1_left_interp, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
  
    call poisson%a11_field%set_field_data( a11_values )
    call poisson%a11_field%update_interpolation_coefficients( )  
   
    poisson%a12_field => sll_f_new_scalar_field_2d_discrete( &
         "a12_check", &
         poisson%interp_a12, &
         transf, &
         bc_eta1_left_interp, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%a12_field%set_field_data( a12_values )
    call poisson%a12_field%update_interpolation_coefficients( ) 
    
    poisson%a21_field => sll_f_new_scalar_field_2d_discrete( &
         "a21_check", &
         poisson%interp_a21, &
         transf, &
         bc_eta1_left_interp, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%a21_field%set_field_data( a21_values )
    call poisson%a21_field%update_interpolation_coefficients( ) 
    
    poisson%a22_field => sll_f_new_scalar_field_2d_discrete( &
         "a22_check", &
         poisson%interp_a22, &
         transf, &
         bc_eta1_left_interp, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%a22_field%set_field_data( a22_values )
    call poisson%a22_field%update_interpolation_coefficients( )

    poisson%b1_field => sll_f_new_scalar_field_2d_discrete( &
         "b1_check", &
         poisson%interp_b1, &
         transf, &
         bc_eta1_left_interp, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%b1_field%set_field_data( b1_values )
    call poisson%b1_field%update_interpolation_coefficients( ) 
    
    poisson%b2_field => sll_f_new_scalar_field_2d_discrete( &
         "b2_check", &
         poisson%interp_b2, &
         transf, &
         bc_eta1_left_interp, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%b2_field%set_field_data( b2_values )
    call poisson%b2_field%update_interpolation_coefficients( )
  
    poisson%c_field => sll_f_new_scalar_field_2d_discrete( &
         "c_check", &
         poisson%interp_c, &
         transf, &
         bc_eta1_left_interp, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)

    call poisson%c_field%set_field_data( c_values )
    call poisson%c_field%update_interpolation_coefficients( )
 
    SLL_ALLOCATE(phi_values(np_eta1,np_eta2),ierr)
    SLL_ALLOCATE(rho_values(np_eta1,np_eta2),ierr)
    phi_values(:,:) = 0.0_f64
    rho_values(:,:) = 0.0_f64
    
    
     poisson%rho_field => sll_f_new_scalar_field_2d_discrete( &
         "rho_check", &
         poisson%interp_rho, &
         transf, &
         bc_eta1_left_interp, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right)
   
    call poisson%rho_field%set_field_data( rho_values )
    call poisson%rho_field%update_interpolation_coefficients( )
    
    ! initialize of general_elliptic_solver
    print *,'sll_f_new_general_elliptic_solver'   
    poisson%elliptic_solver => sll_f_new_general_elliptic_solver_prototype( &
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
        precompute_rhs=poisson%precompute_rhs, &
        rhs_bc1 = bc_interp2d_eta1, &
        rhs_bc2 = bc_interp2d_eta2, &
        use_cubic_splines = use_cubic_splines, &
        with_constraint = with_constraint, &
        zero_mean = zero_mean)                
 
        
    ! compute matrix the field
    print *,'#begin sll_s_factorize_mat_es'
    call sll_s_factorize_mat_es_prototype(&
!    call assemble_mat_es(&
        poisson%elliptic_solver, &
        poisson%a11_field, &
        poisson%a12_field,&
        poisson%a21_field,&
        poisson%a22_field,&
        poisson%b1_field,&
        poisson%b2_field,&
        poisson%c_field)    
    print *,'#end sll_s_factorize_mat_es'
        
 end subroutine initialize_poisson_2d_elliptic_solver
 
 
 function sll_f_new_poisson_2d_elliptic_solver( &
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
   bc_interp2d_eta1, &
   bc_interp2d_eta2, &
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
   c_values, &
   interp_rho, &
   interp_rho_case, &
   control, &
   precompute_rhs, &
   with_constraint, &
   zero_mean ) &
   result(poisson)
   
   class(poisson_2d_elliptic_solver),        pointer  :: poisson
   class(sll_c_coordinate_transformation_2d_base), pointer :: transf
   sll_int32, intent(in)  :: spline_degree_eta1
   sll_int32, intent(in)  :: spline_degree_eta2
   sll_int32, intent(in)  :: num_cells_eta1
   sll_int32, intent(in)  :: num_cells_eta2
   sll_int32, intent(in)  :: bc_eta1_left
   sll_int32, intent(in)  :: bc_eta1_right
   sll_int32, intent(in)  :: bc_eta2_left
   sll_int32, intent(in)  :: bc_eta2_right
   sll_int32, intent(in)  :: bc_interp2d_eta1
   sll_int32, intent(in)  :: bc_interp2d_eta2
   sll_int32, intent(in)  :: quadrature_type1
   sll_int32, intent(in)  :: quadrature_type2
   sll_real64, intent(in) :: eta1_min
   sll_real64, intent(in) :: eta1_max
   sll_real64, intent(in) :: eta2_min
   sll_real64, intent(in) :: eta2_max
   !sll_real64, dimension(:,:), pointer :: phi_values
   !sll_real64, dimension(:,:), pointer :: rho_values
   sll_real64, dimension(:,:)          :: a11_values
   sll_real64, dimension(:,:)          :: a12_values
   sll_real64, dimension(:,:)          :: a21_values
   sll_real64, dimension(:,:)          :: a22_values
   sll_real64, dimension(:,:)          :: b1_values
   sll_real64, dimension(:,:)          :: b2_values
   sll_real64, dimension(:,:)          :: c_values
   sll_int32, intent(in), optional     :: control
   class(sll_c_interpolator_2d), pointer, optional :: interp_rho
   character(len=256), optional  :: interp_rho_case
   logical, intent(in), optional     :: precompute_rhs
   logical, intent(in), optional     :: with_constraint
   logical, intent(in), optional     :: zero_mean
   !sll_int32 :: np_eta1
   !sll_int32 :: np_eta2
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
   bc_interp2d_eta1, &
   bc_interp2d_eta2, &
   eta1_min, &
   eta1_max, &
   eta2_min, &
   eta2_max, &
   a11_values, & 
   a12_values, & 
   a21_values, & 
   a22_values, & 
   b1_values, &
   b2_values, &
   c_values, &
   interp_rho, &
   interp_rho_case, &
   control, &
   precompute_rhs, &
   with_constraint, &
   zero_mean)  
   
  end function sll_f_new_poisson_2d_elliptic_solver
  
  ! solves \Delta phi = -rho in 2d
  subroutine compute_phi_from_rho_2d_elliptic_solver(poisson,phi,rho )
    ! input variables 
    class(poisson_2d_elliptic_solver), target   :: poisson
    sll_real64, dimension(:,:), intent(in)     :: rho
    ! output variables
    sll_real64, dimension(:,:), intent(out)  :: phi
    ! local variables
    !class(sll_t_general_coordinate_elliptic_solver), pointer   :: elliptic_solver
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    !sll_real64 :: eta1
    !sll_real64 :: eta2
    sll_int32  :: i,j
    sll_int32  :: nc_eta1
    sll_int32  :: nc_eta2
     
    ! The supposition is that all fields use the same logical mesh
    !elliptic_solver => poisson%elliptic_solver

   if(poisson%control==SLL_DO_NOTHING_ELLIPTIC_SOLVER)then
     do j=1,nc_eta2+1
       do i=1,nc_eta1+1
         phi(i,j) = 0._f64
       end do
     end do
     return
   endif
   
    delta1    = poisson%elliptic_solver%delta_eta1
    delta2    = poisson%elliptic_solver%delta_eta2
    eta1_min  = poisson%elliptic_solver%eta1_min
    eta2_min  = poisson%elliptic_solver%eta2_min
    nc_eta1   = poisson%elliptic_solver%num_cells1 !+ 1    
    nc_eta2   = poisson%elliptic_solver%num_cells2 !+ 1
    
    
    
    if(poisson%control==SLL_NO_SOLVE_ELLIPTIC_SOLVER)then
      do j=1,nc_eta2+1
        do i=1,nc_eta1+1
           phi(i,j) = 0._f64
        end do
       end do
     return
    endif
    if(poisson%precompute_rhs)then
      call sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype( &
        poisson%elliptic_solver, &
        rho_values=-rho)
      call sll_s_solve_general_coordinates_elliptic_eq_prototype(&
        poisson%elliptic_solver,&
        phi)      
    else
      call poisson%rho_field%set_field_data(-rho)
      call poisson%rho_field%update_interpolation_coefficients( )
      call sll_s_solve_general_coordinates_elliptic_eq_prototype(&
        poisson%elliptic_solver,&
        phi, &
        rho_field=poisson%rho_field)
    endif
            
       !poisson%phi_field)
   
!   do j=1,nc_eta2+1
!        do i=1,nc_eta1+1
!           phi(i,j) = poisson%phi_field%value_at_indices(i,j)
!        end do
!     end do
     
  end subroutine compute_phi_from_rho_2d_elliptic_solver

    ! solves E = -\nabla Phi in 2d
!    subroutine compute_E_from_phi_2d_fft( poisson, phi, E1, E2 )
!      class(sll_t_poisson_2d_fft_solver) :: poisson
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
      E1 = 0._f64
      E2 = 0._f64
      print *,maxval(rho)
      
      if(.not.(associated(poisson%elliptic_solver)))then
        print *,'#poisson%elliptic_solver is not associated'
      endif
      
      stop
      
      
      
      
      !call solve( poisson%poiss, E1, E2, rho)
      
    end subroutine compute_E_from_rho_2d_elliptic_solver


    function l2norm_squarred_2d_elliptic_solver(poisson, coefs_dofs) result(r)
       class( poisson_2d_elliptic_solver), intent(in) :: poisson !< Poisson solver object.
       sll_real64,intent(in)                      :: coefs_dofs(:,:) !< Values of the coefficient vectors for each DoF
       sll_real64                                   :: r
       
       print*, 'l2norm_squared not implemented for 2d elliptic solver.'

     end function l2norm_squarred_2d_elliptic_solver
    
     subroutine compute_rhs_from_function_2d_elliptic_solver(poisson, func, coefs_dofs)
       class( poisson_2d_elliptic_solver)                    :: poisson !< Maxwell solver object.
       procedure(sll_f_function_of_position)          :: func !< Function to be projected.
       sll_real64, intent(out)                        :: coefs_dofs(:) !< Coefficients of the projection.

       print*, 'compute_rhs_from_function not implemented for 2d elliptic solver.'

     end subroutine compute_rhs_from_function_2d_elliptic_solver

     subroutine delete_poisson_2d_elliptic_solver(poisson)
       class( poisson_2d_elliptic_solver)                    :: poisson !< Maxwell solver object.
     end subroutine delete_poisson_2d_elliptic_solver

 end module sll_m_poisson_2d_elliptic_solver
  
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
