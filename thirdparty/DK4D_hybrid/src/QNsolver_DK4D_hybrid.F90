!===========================================================================
!> Quasi-neutrality solver for
!>  4D drift-kinetic hybrid simulation
!> 
!> Rk1: Use the generalized elliptic poisson solver
!> 
!> Rk2: Need profiles and magnetic field in (x,y) but solve the
!>     problem in (eta1,eta2)
!>
!> \date 2014-08-20
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module QNsolver_DK4D_hybrid_module
#include "sll_working_precision.h"
#include "sll_memory.h"

  use equilibrium_DK4D_module
  use fdistribu_DK4D_module
  use field3d_DK4D_module
  use magnetconf_DK4D_module
  use mesh_DK4D_module
  use sll_m_arbitrary_degree_spline_interpolator_1d
  use sll_m_arbitrary_degree_spline_interpolator_2d
  use sll_m_general_coordinate_elliptic_solver
  use sll_m_scalar_field_2d_base
  use sll_m_scalar_field_2d
  use utils_DK4D_module

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: QNsolver_DK4D_hybrid_t

    !> Boundary conditions in (eta1,eta2)
    type(boundary_conditions_2d_t) :: bound_cond

    !> For interpolations
    type(spline_degree_2d_t) :: spline_degree

    !> Elliptic coordinate solver
    type(sll_t_general_coordinate_elliptic_solver), pointer :: QNS

    !> Number of points in (eta1,eta2) direction
    sll_int32 :: Neta1
    sll_int32 :: Neta2

    !> Different terms of the matrix for the elliptic solver
    class(sll_c_scalar_field_2d_base), pointer :: QN_A11 
    class(sll_c_scalar_field_2d_base), pointer :: QN_A12
    class(sll_c_scalar_field_2d_base), pointer :: QN_A21
    class(sll_c_scalar_field_2d_base), pointer :: QN_A22
    class(sll_c_scalar_field_2d_base), pointer :: QN_B1
    class(sll_c_scalar_field_2d_base), pointer :: QN_B2
    class(sll_c_scalar_field_2d_base), pointer :: QN_C
    
    !> Interpolators for the different terms of the matrix 
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_QN_A11
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_QN_A12
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_QN_A21
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_QN_A22
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_QN_B1
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_QN_B2
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_QN_C

    !> RHS term 
    type(field3D_DK4D_t) :: QN_RHS

    !>  Temporary 2D array in (eta1,eta2) for RHS term
    class(sll_c_scalar_field_2d_base), pointer :: QNrhs2d
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_QNrhs_eta1eta2

    !> Electrostatic potential term in (eta1,eta2)
    type(sll_t_scalar_field_2d_discrete), pointer :: Phi2d
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_Phi_eta1eta2

  end type QNsolver_DK4D_hybrid_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Quasi-Neutrality solver: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_QNsolver_DK4D_hybrid( &
      QNsolver, &
      mesh4d, &
      bound_cond, &
      spline_degree )

    type(QNsolver_DK4D_hybrid_t)  , intent(inout) :: QNsolver
    type(mesh_DK4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_3d_t), intent(in)    :: bound_cond
    type(spline_degree_3d_t)      , intent(in)    :: spline_degree

    type(sll_t_cartesian_mesh_2d), pointer :: cartesian_mesh2d
    sll_int32 :: ierr

    !*** Initialization of the boundary conditions ***
    QNsolver%bound_cond%left_eta1  = bound_cond%left_eta1
    QNsolver%bound_cond%right_eta1 = bound_cond%right_eta1
    QNsolver%bound_cond%left_eta2  = bound_cond%left_eta2
    QNsolver%bound_cond%right_eta2 = bound_cond%right_eta2

    !*** Initialization for interpolations ***
    QNsolver%spline_degree%eta1 = spline_degree%eta1
    QNsolver%spline_degree%eta2 = spline_degree%eta2

    !*** Initialization of the number of points   ***
    !***   in (eta1,eta2) direction               ***
    QNsolver%Neta1 = size(mesh4d%eta1_grid)
    QNsolver%Neta2 = size(mesh4d%eta2_grid)

    !*** Initialization of the interpolators      ***
    !***   for the different terms of the matrix  ***
    cartesian_mesh2d => mesh4d%eta1_eta2_mesh2d

    call QNsolver%interp2d_QN_A11%init( &
      cartesian_mesh2d%num_cells1 +1, &
      cartesian_mesh2d%num_cells2 +1, &
      cartesian_mesh2d%eta1_min, &
      cartesian_mesh2d%eta1_max, &
      cartesian_mesh2d%eta2_min, &
      cartesian_mesh2d%eta2_max, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2, &
      QNsolver%spline_degree%eta1, &
      QNsolver%spline_degree%eta2)    

    call QNsolver%interp2d_QN_A12%init( &
      cartesian_mesh2d%num_cells1 +1, &
      cartesian_mesh2d%num_cells2 +1, &
      cartesian_mesh2d%eta1_min, &
      cartesian_mesh2d%eta1_max, &
      cartesian_mesh2d%eta2_min, &
      cartesian_mesh2d%eta2_max, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2, &
      QNsolver%spline_degree%eta1, &
      QNsolver%spline_degree%eta2)    

    call QNsolver%interp2d_QN_A21%init( &
      cartesian_mesh2d%num_cells1 +1, &
      cartesian_mesh2d%num_cells2 +1, &
      cartesian_mesh2d%eta1_min, &
      cartesian_mesh2d%eta1_max, &
      cartesian_mesh2d%eta2_min, &
      cartesian_mesh2d%eta2_max, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2, &
      QNsolver%spline_degree%eta1, &
      QNsolver%spline_degree%eta2)    

    call QNsolver%interp2d_QN_A22%init( &
      cartesian_mesh2d%num_cells1 +1, &
      cartesian_mesh2d%num_cells2 +1, &
      cartesian_mesh2d%eta1_min, &
      cartesian_mesh2d%eta1_max, &
      cartesian_mesh2d%eta2_min, &
      cartesian_mesh2d%eta2_max, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2, &
      QNsolver%spline_degree%eta1, &
      QNsolver%spline_degree%eta2)    

    call QNsolver%interp2d_QN_B1%init( &
      cartesian_mesh2d%num_cells1 +1, &
      cartesian_mesh2d%num_cells2 +1, &
      cartesian_mesh2d%eta1_min, &
      cartesian_mesh2d%eta1_max, &
      cartesian_mesh2d%eta2_min, &
      cartesian_mesh2d%eta2_max, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2, &
      QNsolver%spline_degree%eta1, &
      QNsolver%spline_degree%eta2)
        
    call QNsolver%interp2d_QN_B2%init( &
      cartesian_mesh2d%num_cells1 +1, &
      cartesian_mesh2d%num_cells2 +1, &
      cartesian_mesh2d%eta1_min, &
      cartesian_mesh2d%eta1_max, &
      cartesian_mesh2d%eta2_min, &
      cartesian_mesh2d%eta2_max, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2, &
      QNsolver%spline_degree%eta1, &
      QNsolver%spline_degree%eta2)
    
    call QNsolver%interp2d_QN_C%init( &
      cartesian_mesh2d%num_cells1 +1, &
      cartesian_mesh2d%num_cells2 +1, &
      cartesian_mesh2d%eta1_min, &
      cartesian_mesh2d%eta1_max, &
      cartesian_mesh2d%eta2_min, &
      cartesian_mesh2d%eta2_max, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2, &
      QNsolver%spline_degree%eta1, &
      QNsolver%spline_degree%eta2)    

    !*** Allocation of the 2D fields associated to ***
    !***  A11, A12, A21, A22, B1, B2 and C         ***
    QNsolver%QN_A11 => sll_f_new_scalar_field_2d_discrete( &
        "QN_A11", &
        QNsolver%interp2d_QN_A11, &     
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2)

    QNsolver%QN_A12 => sll_f_new_scalar_field_2d_discrete( &
        "QN_A12", &
        QNsolver%interp2d_QN_A12, &     
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2)

    QNsolver%QN_A21 => sll_f_new_scalar_field_2d_discrete( &
        "QN_A21", &
        QNsolver%interp2d_QN_A21, &     
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2)

    QNsolver%QN_A22 => sll_f_new_scalar_field_2d_discrete( &
        "QN_A22", &
        QNsolver%interp2d_QN_A22, &     
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2)

    QNsolver%QN_B1 => sll_f_new_scalar_field_2d_discrete( &
        "QN_B1", &
        QNsolver%interp2d_QN_B1, &     
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2)

    QNsolver%QN_B2 => sll_f_new_scalar_field_2d_discrete( &
        "QN_B2", &
        QNsolver%interp2d_QN_B1, &     
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2)

    QNsolver%QN_C => sll_f_new_scalar_field_2d_discrete( &
        "QN_C", &
        QNsolver%interp2d_QN_C, &     
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2)

    !*** Allocation of the QNS type ***
    cartesian_mesh2d => mesh4d%transf_eta1eta2_xy%mesh

    QNsolver%QNS => sll_f_new_general_elliptic_solver( &
        QNsolver%spline_degree%eta1, & 
        QNsolver%spline_degree%eta2, & 
        cartesian_mesh2d%num_cells1, &
        cartesian_mesh2d%num_cells2, &
        sll_p_es_gauss_legendre, &  ! put in arguments
        sll_p_es_gauss_legendre, &  ! put in arguments
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2, &
        cartesian_mesh2d%eta1_min, &  
        cartesian_mesh2d%eta1_max, & 
        cartesian_mesh2d%eta2_min, & 
        cartesian_mesh2d%eta2_max ) 

    !*** Initialization of the 3D field for RHS term ***
    call new_field3D_DK4D( &
        QNsolver%QN_RHS, &
        'QN_RHS', &
        mesh4d, &
        bound_cond, &
        spline_degree )

    !*** Initialization of the interpolator and the 2D field ***
    !***  for the RHS in 2D                                  ***
    !--> Interpolator
    call QNsolver%interp2d_QNrhs_eta1eta2%init( &
      cartesian_mesh2d%num_cells1 +1, &
      cartesian_mesh2d%num_cells2 +1, &
      cartesian_mesh2d%eta1_min, &
      cartesian_mesh2d%eta1_max, &
      cartesian_mesh2d%eta2_min, &
      cartesian_mesh2d%eta2_max, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2, &
      QNsolver%spline_degree%eta1, &
      QNsolver%spline_degree%eta2)    

    !--> QNrhs2D field in (eta1,eta2)
    QNsolver%QNrhs2d => sll_f_new_scalar_field_2d_discrete( &
      "QNrhs2d_seqx1x2", &
      QNsolver%interp2d_QNrhs_eta1eta2, &     
      mesh4d%transf_eta1eta2_xy, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2)

    !*** Initialization of the interpolator and the 2D field ***
    !***  for Phi in2D                                       ***
    !--> Interpolator
    call QNsolver%interp2d_Phi_eta1eta2%init( &
      cartesian_mesh2d%num_cells1+1, &
      cartesian_mesh2d%num_cells2+1, &
      cartesian_mesh2d%eta1_min, &
      cartesian_mesh2d%eta1_max, &
      cartesian_mesh2d%eta2_min, &
      cartesian_mesh2d%eta2_max, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2, &
      QNsolver%spline_degree%eta1, &
      QNsolver%spline_degree%eta2)

    !--> Phi2D field in (eta1,eta2)
    QNsolver%phi2d => sll_f_new_scalar_field_2d_discrete( &
      "phi2d_seqx1x2", &
      QNsolver%interp2d_Phi_eta1eta2, &     
      mesh4d%transf_eta1eta2_xy, &
      QNsolver%bound_cond%left_eta1, &
      QNsolver%bound_cond%right_eta1, &
      QNsolver%bound_cond%left_eta2, &
      QNsolver%bound_cond%right_eta2)

  end subroutine new_QNsolver_DK4D_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_QNsolver_DK4D_hybrid( QNsolver )

    type(QNsolver_DK4D_hybrid_t), intent(inout) :: QNsolver

    sll_int32 :: ierr
    
    call QNsolver%interp2d_QN_A11%delete( )
    call QNsolver%interp2d_QN_A12%delete( )
    call QNsolver%interp2d_QN_A21%delete( )
    call QNsolver%interp2d_QN_A22%delete( )
    call QNsolver%interp2d_QN_B1%delete( )
    call QNsolver%interp2d_QN_B2%delete( )
    call QNsolver%interp2d_QN_C%delete( )
    call QNsolver%interp2d_QNrhs_eta1eta2%delete( )
    call QNsolver%interp2d_Phi_eta1eta2%delete( )

    call QNsolver%QN_A11%free( )
    call QNsolver%QN_A11%free( )
    call QNsolver%QN_A12%free( )
    call QNsolver%QN_A21%free( )
    call QNsolver%QN_A22%free( )
    call QNsolver%QN_B1%free( )
    call QNsolver%QN_B2%free( )
    call QNsolver%QN_C%free( )
    call QNsolver%QNrhs2d%free( )
    call QNsolver%phi2d%free( )

    call sll_o_delete( QNsolver%QNS ) 

    call delete_field3D_DK4D( QNsolver%QN_RHS )

  end subroutine delete_QNsolver_DK4D_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Initialization
  !---------------------------------------------------------------------------
  subroutine init_QNsolver_DK4D_hybrid( QNsolver, &
      mesh4d, &
      magnetconf, &
      equilibrium )

    type(QNsolver_DK4D_hybrid_t), intent(inout) :: QNsolver
    type(mesh_DK4D_t)           , intent(in)    :: mesh4d
    type(magnetconf_DK4D_t)     , intent(in)    :: magnetconf
    type(equilibrium_DK4D_t)    , intent(in)    :: equilibrium

    sll_int32 :: ierr
    sll_int32 :: Neta1, Neta2
    sll_real64, dimension(:,:), pointer :: A11
    sll_real64, dimension(:,:), pointer :: A12
    sll_real64, dimension(:,:), pointer :: A21
    sll_real64, dimension(:,:), pointer :: A22
    sll_real64, dimension(:,:), pointer :: B1
    sll_real64, dimension(:,:), pointer :: B2
    sll_real64, dimension(:,:), pointer :: C

    Neta1 = QNsolver%Neta1
    Neta2 = QNsolver%Neta2
    SLL_ALLOCATE( A11(Neta1,Neta2), ierr )
    SLL_ALLOCATE( A12(Neta1,Neta2), ierr )
    SLL_ALLOCATE( A21(Neta1,Neta2), ierr )
    SLL_ALLOCATE( A22(Neta1,Neta2), ierr )
    SLL_ALLOCATE( B1(Neta1,Neta2), ierr )
    SLL_ALLOCATE( B2(Neta1,Neta2), ierr )
    SLL_ALLOCATE( C(Neta1,Neta2), ierr )
    
    !---> Initialization of the matrices A11, A12, A21, A22, B1, B2 and C
    call initialize_matrix_A_QN_DK ( &
        QNsolver, &
        magnetconf, &
        equilibrium, &
        A11, A12, A21, A22 ) 
    call initialize_vector_B_QN_DK ( &
        QNsolver, &
        mesh4d, & 
        magnetconf, &
        equilibrium, &
        B1, B2 ) 
    call initialize_scalar_C_QN_DK ( &
        QNsolver, &
        equilibrium, &
        C )
       
    !*** Initialization of the 2D fields associated to ***
    !***  A11, A12, A21, A22, B1, B2 and C             ***
    call QNsolver%QN_A11%set_field_data( A11 )    
    call QNsolver%QN_A11%update_interpolation_coefficients( )

    call QNsolver%QN_A12%set_field_data( A12 )
    call QNsolver%QN_A12%update_interpolation_coefficients( )

    call QNsolver%QN_A21%set_field_data( A21 )
    call QNsolver%QN_A21%update_interpolation_coefficients( )

    call QNsolver%QN_A22%set_field_data( A22 )
    call QNsolver%QN_A22%update_interpolation_coefficients( )

    call QNsolver%QN_B1%set_field_data( B1 )
    call QNsolver%QN_B1%update_interpolation_coefficients( )

    call QNsolver%QN_B2%set_field_data( B1 )
    call QNsolver%QN_B2%update_interpolation_coefficients( )

    call QNsolver%QN_C%set_field_data( C )
    call QNsolver%QN_C%update_interpolation_coefficients( )

    SLL_DEALLOCATE( A11, ierr )
    SLL_DEALLOCATE( A12, ierr )
    SLL_DEALLOCATE( A21, ierr )
    SLL_DEALLOCATE( A22, ierr )
    SLL_DEALLOCATE( B1, ierr )
    SLL_DEALLOCATE( B2, ierr )
    SLL_DEALLOCATE( C, ierr )

  end subroutine init_QNsolver_DK4D_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Factorizing
  !---------------------------------------------------------------------------
  subroutine factorize_QNsolver_DK4D_hybrid( QNsolver )

    type(QNsolver_DK4D_hybrid_t), intent(inout) :: QNsolver

    call sll_s_factorize_mat_es( &
        QNsolver%QNS, & 
        QNsolver%QN_A11, &
        QNsolver%QN_A12, &
        QNsolver%QN_A21, &
        QNsolver%QN_A22, &
        QNsolver%QN_B1,  &
        QNsolver%QN_B2,  &
        QNsolver%QN_C )

  end subroutine factorize_QNsolver_DK4D_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Computation of the Right-Hand-Side 
  !>  of the equation i.e
  !>   QN_RHS(x,y,eta3) = \int delta f dvpar
  !>    with delta f = f(x,y,eta3,vpar) - feq(x,y,vpar) 
  !> 
  !>  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !>       feq3d_seqx1x2x4(x1=*,x2=*,x4=*)
  !> Out : QNsolver%QN_RHS 3D field
  !---------------------------------------------------------------------------
  subroutine compRHS_QNsolver_DK4D_hybrid( QNsolver, &
      mesh4d, &
      equilibrium, &
      fdistribu )

    type(QNsolver_DK4D_hybrid_t), intent(inout) :: QNsolver
    type(mesh_DK4D_t)           , intent(in)    :: mesh4d
    type(equilibrium_DK4D_t)    , intent(in)    :: equilibrium
    type(fdistribu_DK4D_t)      , intent(in)    :: fdistribu

    sll_int32  :: Neta1_loc, Neta2_loc
    sll_int32  :: Neta3, Nvpar
    sll_int32  :: iloc1, iloc2
    sll_int32  :: i1, i2, i3, i4
    sll_real64 :: delta_f
    sll_real64 :: delta_vpar, intf_dvpar
    sll_int32, dimension(1:4) :: glob_ind4d

    Neta1_loc  = size(fdistribu%val4d_seqx3x4,1)
    Neta2_loc  = size(fdistribu%val4d_seqx3x4,2)
    Neta3      = size(fdistribu%val4d_seqx3x4,3)
    Nvpar      = size(fdistribu%val4d_seqx3x4,4)
    delta_vpar = mesh4d%vpar_mesh1d%delta_eta

    !*** Computation of the RHS locally in (x1,x2) directions ***
    do i3 = 1,Neta3
      do iloc2 = 1,Neta2_loc
        do iloc1 = 1,Neta1_loc
          intf_dvpar = 0._f64
          do i4 = 1,Nvpar-1
            glob_ind4d(:) = sll_o_local_to_global( &
                fdistribu%layout4d_seqx3x4, &
                (/iloc1,iloc2,i3,i4/))
            i1 = glob_ind4d(1)
            i2 = glob_ind4d(2)
            delta_f = fdistribu%val4d_seqx3x4(iloc1,iloc2,i3,i4) - &
                equilibrium%feq_xyvpar(i1,i2,i4)
            intf_dvpar = intf_dvpar + delta_f*delta_vpar           
          end do
          QNsolver%QN_RHS%val3d_seqx3(iloc1,iloc2,i3) = intf_dvpar
        end do
      end do
    end do

    !--> compute rho3d_seqx1x2
    call sll_o_apply_remap_3d( &
        QNsolver%QN_RHS%seqx3_to_seqx1x2, &
        QNsolver%QN_RHS%val3d_seqx3, &
        QNsolver%QN_RHS%val3d_seqx1x2 )

  end subroutine compRHS_QNsolver_DK4D_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Solving
  !---------------------------------------------------------------------------
  subroutine solve_QNsolver_DK4D_hybrid( &
      QNsolver, &
      mesh4d, &
      equilibrium, &
      fdistribu, &
      Phi )

    type(QNsolver_DK4D_hybrid_t), intent(inout) :: QNsolver
    type(mesh_DK4D_t)           , intent(in)    :: mesh4d
    type(equilibrium_DK4D_t)    , intent(in)    :: equilibrium
    type(fdistribu_DK4D_t)      , intent(in)    :: fdistribu
    type(field3d_DK4D_t)        , intent(inout) :: Phi

    sll_int32 :: ieta1, ieta2, iloc3
    sll_int32 :: Neta1, Neta2
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3

    !*** Computation of the RHS of the quasi-neutrality equation ***
    call compRHS_QNsolver_DK4D_hybrid( QNsolver, &
      mesh4d, &
      equilibrium, &
      fdistribu )

    !*** Solving of the quasi-neutrality equation                     ***
    !*** ==> Compute Phi(x,y,eta3) assuming the parallel decomposition ***
    !***      sequential in (x,y) and parallel in eta3 direction       ***
    Neta1 = QNsolver%Neta1
    Neta2 = QNsolver%Neta2

    call sll_o_compute_local_sizes( Phi%layout3d_seqx1x2, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)

    do iloc3 = 1,loc3d_sz_x3
      call QNsolver%QNrhs2d%set_field_data( &
          QNsolver%QN_RHS%val3d_seqx1x2(:,:,iloc3) )
      call QNsolver%QNrhs2d%update_interpolation_coefficients( )
      call sll_o_solve( &
        QNsolver%QNS, &
        QNsolver%QNrhs2d, &
        QNsolver%Phi2d )
      do ieta2 = 1, Neta2
        do ieta1 = 1, Neta1 
           Phi%val3d_seqx1x2(ieta1,ieta2,iloc3) = &
               QNsolver%Phi2d%value_at_indices(ieta1,ieta2)
        end do
      end do
    end do
    
    !*** Fill Phi%val3d_seqx3 ***
    call sll_o_apply_remap_3d(Phi%seqx1x2_to_seqx3, &
        Phi%val3d_seqx1x2, Phi%val3d_seqx3)    

  end subroutine solve_QNsolver_DK4D_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Initialization of the QN coefficients of matrix A
  !----------------------------------------------------
  ! In the case of drift-Kinetic and with F the change of variables such that 
  !
  !  F( eta1,eta2) =( F_1(eta1, eta2),F_2(eta1, eta2) ) = ( x, y )
  !
  !     ( -1   0 )
  ! A = (  0  -1 ) *  n_0 ( F_1(eta1, eta2),F_2(eta1, eta2)) / 
  !               (B(F_1(eta1, eta2),F_2(eta1, eta2)))
  !---------------------------------------------------------------------------
  subroutine initialize_matrix_A_QN_DK( &
      QNsolver, &
      magnetconf, &
      equilibrium, &
      values_A11, &
      values_A12, &
      values_A21, &
      values_A22 )

    type(QNsolver_DK4D_hybrid_t), intent(inout) :: QNsolver
    type(magnetconf_DK4D_t)     , intent(in)    :: magnetconf
    type(equilibrium_DK4D_t)    , intent(in)    :: equilibrium
    sll_real64  , dimension(:,:), intent(inout) :: values_A11
    sll_real64  , dimension(:,:), intent(inout) :: values_A12
    sll_real64  , dimension(:,:), intent(inout) :: values_A21
    sll_real64  , dimension(:,:), intent(inout) :: values_A22

    sll_int32  :: ieta1, ieta2
    sll_int32  :: Neta1, Neta2
    sll_real64 :: n0_xy_tmp, B_xy_tmp

    Neta1 = QNsolver%Neta1
    Neta2 = QNsolver%Neta2

    if ( (size(values_A11,1) .ne. Neta1) .OR. &
        (size(values_A11,2) .ne. Neta2) ) then
      print*, ' Problem with the dimension of A11'
    end if
    if ( (size(values_A12,1) .ne. Neta1) .OR. &
        (size(values_A12,2) .ne. Neta2) ) then
      print*, ' Problem with the dimension of A12'
    end if
    if ( (size(values_A21,1) .ne. Neta1) .OR. &
        (size(values_A21,2) .ne. Neta2) ) then
      print*, ' Problem with the dimension of A21'
    end if
    if ( (size(values_A22,1) .ne. Neta1) .OR. &
        (size(values_A22,2) .ne. Neta2) ) then
      print*, ' Problem with the dimension of A22'
    end if

    do ieta2 = 1,Neta2
      do ieta1 = 1,Neta1
        n0_xy_tmp = equilibrium%n0_xy(ieta1,ieta2)
        B_xy_tmp  = magnetconf%B_xy(ieta1,ieta2) 
        values_A11(ieta1,ieta2) = &
            - n0_xy_tmp/B_xy_tmp
        values_A12(ieta1,ieta2) = 0.0_f64 
        values_A21(ieta1,ieta2) = 0.0_f64
        values_A22(ieta1,ieta2) = &
            - n0_xy_tmp/B_xy_tmp
      end do
    end do

  end subroutine initialize_matrix_A_QN_DK
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Initialization of the QN coefficients of vector B
  !----------------------------------------------------
  ! In the case of drift-Kinetic and with F the change of variables such that 
  !
  !       ( 0 )
  ! B =   ( 0 ) 
  !---------------------------------------------------------------------------
  subroutine initialize_vector_B_QN_DK( &
      QNsolver, &
      mesh4d, & 
      magnetconf, &
      equilibrium, &
      values_B1, &
      values_B2 )

    type(QNsolver_DK4D_hybrid_t), intent(inout) :: QNsolver
    type(mesh_DK4D_t)           , intent(in)    :: mesh4d
    type(magnetconf_DK4D_t)     , intent(in)    :: magnetconf
    type(equilibrium_DK4D_t)    , intent(in)    :: equilibrium
    sll_real64  , dimension(:,:), pointer       :: values_B1
    sll_real64  , dimension(:,:), pointer       :: values_B2

    sll_int32  :: ieta1, ieta2
    sll_int32  :: Neta1, Neta2
    sll_real64 :: x_tmp, y_tmp
    sll_real64 :: n0_xy_tmp, B_xy_tmp
    sll_real64 :: norm2_xy_tmp, val_tmp

    Neta1 = QNsolver%Neta1
    Neta2 = QNsolver%Neta2

    if ( (size(values_B1,1) .ne. Neta1) .OR. &
        (size(values_B1,2) .ne. Neta2) ) then
      print*, ' Problem with the dimension of B1'
    end if
    if ( (size(values_B2,1) .ne. Neta1) .OR. &
        (size(values_B2,2) .ne. Neta2) ) then
      print*, ' Problem with the dimension of B2'
    end if

    do ieta2 = 1,Neta2
      do ieta1 = 1,Neta1
        values_B1(ieta1,ieta2) = 0.0
        values_B2(ieta1,ieta2) = 0.0
      end do
    end do

  end subroutine initialize_vector_B_QN_DK
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Initialization of the QN coefficients of scalar C
  !----------------------------------------------------
  ! In the case of drift-Kinetic and with F the change of variables such that 
  !
  ! C = e n_0 ( F_1(eta1, eta2),F_2(eta1, eta2)) / 
  !           ( T_e (F_1(eta1, eta2),F_2(eta1, eta2)) )
  !---------------------------------------------------------------------------
  subroutine initialize_scalar_C_QN_DK( &
      QNsolver, &
      equilibrium, &
      values_C )

    type(QNsolver_DK4D_hybrid_t), intent(inout) :: QNsolver
    type(equilibrium_DK4D_t)    , intent(in)    :: equilibrium
    sll_real64  , dimension(:,:), pointer       :: values_C

    sll_int32  :: ieta1, ieta2
    sll_int32  :: Neta1, Neta2

    Neta1 = QNsolver%Neta1
    Neta2 = QNsolver%Neta2

    if ( (size(values_C,1) .ne. Neta1) .or. &
        (size(values_C,2) .ne. Neta2) ) then
      print*, ' Problem with the dimension of C'
    end if

    do ieta2 = 1,Neta2
      do ieta1 = 1,Neta1
        values_C(ieta1,ieta2) = &
            equilibrium%n0_xy(ieta1,ieta2) / &
            equilibrium%Te_xy(ieta1,ieta2)
      end do
    end do

  end subroutine initialize_scalar_C_QN_DK
  !---------------------------------------------------------------------------

end module QNsolver_DK4D_hybrid_module
!---------------------------------------------------------------------------
