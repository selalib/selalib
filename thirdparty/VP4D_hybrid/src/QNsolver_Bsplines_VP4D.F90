!===========================================================================
!> Quasi-neutrality solver for
!>  4D Vlasov-Poisson hybrid simulation
!> 
!> Rk1: Use the generalized elliptic poisson solver
!>
!> \date 2015-02-26
!> \author V. Grandgirard, A. Back
!---------------------------------------------------------------------------
module QNsolver_Bsplines_VP4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"

  use equilibrium_VP4D_module
  use fdistribu_VP4D_module
  use field2d_VP4D_module
  use mesh_VP4D_module
  use sll_m_arbitrary_degree_spline_interpolator_2d
  use sll_m_general_coordinate_elliptic_solver
  use sll_m_scalar_field_2d_base
  use sll_m_scalar_field_2d
  use utils_VP4D_module

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: QNsolver_Bsplines_VP4D_t

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
    class(sll_c_scalar_field_2d_base), pointer :: QN_A11_mat
    class(sll_c_scalar_field_2d_base), pointer :: QN_A12_mat
    class(sll_c_scalar_field_2d_base), pointer :: QN_A21_mat
    class(sll_c_scalar_field_2d_base), pointer :: QN_A22_mat
    class(sll_c_scalar_field_2d_base), pointer :: QN_B1_vect
    class(sll_c_scalar_field_2d_base), pointer :: QN_B2_vect
    class(sll_c_scalar_field_2d_base), pointer :: QN_C

    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_A11_func 
    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_A12_func
    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_A21_func
    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_A22_func
    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_B1_func
    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_deriv1_B1_func
    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_deriv2_B1_func
    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_B2_func
    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_deriv1_B2_func
    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_deriv2_B2_func
    procedure(sll_i_two_var_parametrizable_function), nopass , pointer :: QN_C_func
    sll_real64, dimension(:), pointer :: QN_A11_params
    sll_real64, dimension(:), pointer :: QN_A12_params
    sll_real64, dimension(:), pointer :: QN_A21_params
    sll_real64, dimension(:), pointer :: QN_A22_params
    sll_real64, dimension(:), pointer :: QN_B1_params
    sll_real64, dimension(:), pointer :: QN_B2_params
    sll_real64, dimension(:), pointer :: QN_C_params
    
    !> RHS term 
    type(field2D_VP4D_t) :: QN_RHS

    !>  Temporary 2D array in (eta1,eta2) for RHS term
    class(sll_c_scalar_field_2d_base), pointer :: QNrhs2d
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_QNrhs_eta1eta2

    !> Electrostatic potential term in (eta1,eta2)
    type(sll_t_scalar_field_2d_discrete), pointer :: Phi2d
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_Phi_eta1eta2

  end type QNsolver_Bsplines_VP4D_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Quasi-Neutrality solver: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_QNsolver_Bsplines_VP4D( &
      QNsolver, &
      mesh4d, &
      bound_cond, &
      spline_degree )

    type(QNsolver_Bsplines_VP4D_t), intent(inout) :: QNsolver
    type(mesh_VP4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_2d_t), intent(in)    :: bound_cond
    type(spline_degree_2d_t)      , intent(in)    :: spline_degree

    class(sll_t_cartesian_mesh_2d), pointer :: cartesian_mesh2d
    sll_int32 :: ierr
    sll_real64, dimension(1) :: func_zero_params
    sll_real64, dimension(1) :: func_minus_one_params
    sll_real64, dimension(:), pointer :: A11_params_def
    sll_real64, dimension(:), pointer :: A12_params_def
    sll_real64, dimension(:), pointer :: A21_params_def
    sll_real64, dimension(:), pointer :: A22_params_def
    sll_real64, dimension(:), pointer :: B1_params_def
    sll_real64, dimension(:), pointer :: B2_params_def
    sll_real64, dimension(:), pointer :: C_params_def

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

    !*** Allocation of the 2D fields associated to ***
    !***  A11, A12, A21, A22, B1, B2 and C         ***
    QNsolver%QN_A11_func       => func_minus_one
    QNsolver%QN_A12_func       => func_zero
    QNsolver%QN_A21_func       => func_zero
    QNsolver%QN_A22_func       => func_minus_one
    QNsolver%QN_B1_func        => func_zero
    QNsolver%QN_deriv1_B1_func => func_zero
    QNsolver%QN_deriv2_B1_func => func_zero
    QNsolver%QN_B2_func        => func_zero
    QNsolver%QN_deriv1_B2_func => func_zero
    QNsolver%QN_deriv2_B2_func => func_zero
    QNsolver%QN_C_func         => func_zero

    func_zero_params(:)      = (/0.0_f64/)
    func_minus_one_params(:) = (/-1.0_f64/)  

    SLL_ALLOCATE( A11_params_def(size(func_zero_params)), ierr )
    SLL_ALLOCATE( A12_params_def(size(func_zero_params)), ierr )
    SLL_ALLOCATE( A21_params_def(size(func_zero_params)), ierr )
    SLL_ALLOCATE( A22_params_def(size(func_zero_params)), ierr )
    SLL_ALLOCATE( B1_params_def(size(func_zero_params)), ierr )
    SLL_ALLOCATE( B2_params_def(size(func_zero_params)), ierr )
    SLL_ALLOCATE( C_params_def(size(func_zero_params)), ierr )    

    A11_params_def(:) = func_minus_one_params
    A12_params_def(:) = func_zero_params
    A21_params_def(:) = func_zero_params
    A22_params_def(:) = func_minus_one_params
    B1_params_def(:)  = func_zero_params
    B2_params_def(:)  = func_zero_params
    C_params_def(:)   = func_zero_params

    SLL_ALLOCATE( QNsolver%QN_A11_params(size(A11_params_def)), ierr )
    SLL_ALLOCATE( QNsolver%QN_A12_params(size(A12_params_def)), ierr )
    SLL_ALLOCATE( QNsolver%QN_A21_params(size(A21_params_def)), ierr )
    SLL_ALLOCATE( QNsolver%QN_A22_params(size(A22_params_def)), ierr )
    SLL_ALLOCATE( QNsolver%QN_B1_params(size(B1_params_def)), ierr )
    SLL_ALLOCATE( QNsolver%QN_B2_params(size(B2_params_def)), ierr )
    SLL_ALLOCATE( QNsolver%QN_C_params(size(C_params_def)), ierr )

    QNsolver%QN_A11_mat => sll_f_new_scalar_field_2d_analytic( &
        QNsolver%QN_A11_func, &
        "QN_A11", &
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2, &
        QNsolver%QN_A11_params )

    QNsolver%QN_A12_mat => sll_f_new_scalar_field_2d_analytic( &
        QNsolver%QN_A12_func, &
        "QN_A12", &
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2, &
        QNsolver%QN_A12_params )

    QNsolver%QN_A21_mat => sll_f_new_scalar_field_2d_analytic( &
        QNsolver%QN_A21_func, &
        "QN_A21", &
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2, &
        QNsolver%QN_A21_params )

    QNsolver%QN_A22_mat => sll_f_new_scalar_field_2d_analytic( &
        QNsolver%QN_A22_func, &
        "QN_A22", &
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2, &
        QNsolver%QN_A22_params )

    QNsolver%QN_B1_vect => sll_f_new_scalar_field_2d_analytic( &
        QNsolver%QN_B1_func, &
        "QN_B1", &
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2, &
        QNsolver%QN_B1_params, &
        QNsolver%QN_deriv1_B1_func, &
        QNsolver%QN_deriv2_B1_func )

    QNsolver%QN_B2_vect => sll_f_new_scalar_field_2d_analytic( &
        QNsolver%QN_B2_func, &
        "QN_B2", &
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2, &
        QNsolver%QN_B2_params, &
        QNsolver%QN_deriv1_B2_func, &
        QNsolver%QN_deriv2_B2_func )

    QNsolver%QN_C => sll_f_new_scalar_field_2d_analytic( &
        QNsolver%QN_C_func, &
        "QN_C", &
        mesh4d%transf_eta1eta2_xy, &
        QNsolver%bound_cond%left_eta1, &
        QNsolver%bound_cond%right_eta1, &
        QNsolver%bound_cond%left_eta2, &
        QNsolver%bound_cond%right_eta2, &
        QNsolver%QN_C_params )

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

    !*** Initialization of the 2D field for RHS term ***
    call new_field2D_VP4D( &
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

  end subroutine new_QNsolver_Bsplines_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_QNsolver_Bsplines_VP4D( QNsolver )

    type(QNsolver_Bsplines_VP4D_t), intent(inout) :: QNsolver

    call QNsolver%interp2d_QNrhs_eta1eta2%delete( )
    call QNsolver%interp2d_Phi_eta1eta2%delete( )

    call QNsolver%QN_A11_mat%free( )
    call QNsolver%QN_A11_mat%free( )
    call QNsolver%QN_A12_mat%free( )
    call QNsolver%QN_A21_mat%free( )
    call QNsolver%QN_A22_mat%free( )
    call QNsolver%QN_B1_vect%free( )
    call QNsolver%QN_B2_vect%free( )
    call QNsolver%QN_C%free( )
    call QNsolver%QNrhs2d%free( )
    call QNsolver%phi2d%free( )

    call sll_o_delete( QNsolver%QNS ) 

    call delete_field2D_VP4D( QNsolver%QN_RHS )

  end subroutine delete_QNsolver_Bsplines_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Factorizing
  !---------------------------------------------------------------------------
  subroutine factorize_QNsolver_Bsplines_VP4D( QNsolver )

    type(QNsolver_Bsplines_VP4D_t), intent(inout) :: QNsolver

    call sll_s_factorize_mat_es( &
        QNsolver%QNS, & 
        QNsolver%QN_A11_mat, &
        QNsolver%QN_A12_mat, &
        QNsolver%QN_A21_mat, &
        QNsolver%QN_A22_mat, &
        QNsolver%QN_B1_vect,  &
        QNsolver%QN_B2_vect,  &
        QNsolver%QN_C )

  end subroutine factorize_QNsolver_Bsplines_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Computation of the Right-Hand-Side 
  !>  of the equation i.e
  !>   QN_RHS(x,y,eta3) = \int delta f dvpar
  !>    with delta f = f(x,y,eta3,vpar) - feq(x,y,vpar) 
  !> 
  !>  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !>       feq2d_seqx1x2x4(x1=*,x2=*,x4=*)
  !> Out : QNsolver%QN_RHS 2D field
  !---------------------------------------------------------------------------
  subroutine compRHS_QNsolver_Bsplines_VP4D( QNsolver, &
      mesh4d, &
      equilibrium, &
      fdistribu )

    type(QNsolver_Bsplines_VP4D_t), intent(inout) :: QNsolver
    type(mesh_VP4D_t)             , intent(in)    :: mesh4d
    type(equilibrium_VP4D_t)      , intent(in)    :: equilibrium
    type(fdistribu_VP4D_t)        , intent(in)    :: fdistribu

    !--> Local variables
    !----> Used to compute rho parallelized in (x1,x2)
    sll_int32  :: Neta1, Neta2
    sll_int32  :: Neta1_loc, Neta2_loc
    sll_int32  :: Nvx, Nvy
    sll_int32  :: iloc1, iloc2
    sll_int32  :: i1, i2, i3, i4
    sll_real64 :: delta_f
    sll_real64 :: dv_space, intf_dvxdvy
    sll_int32, dimension(1:4) :: glob_ind4d
    !----> Used to compute rho sequential in (x1,x2) 
    !---->  (for allgatherv operation)
    !sll_int32 :: world_size
    !sll_int32 :: send_size   
    !sll_real64, dimension(:), allocatable :: send_buf
    !sll_real64, dimension(:), allocatable :: recv_buf
    !sll_int32 , dimension(:), allocatable :: recv_sz
    !sll_int32 , dimension(:), allocatable :: disps 
    !---> Used to compute \int rho |J| deta1 deta2/ \int |J| deta1 deta2
    sll_real64 :: jacob_tmp, val_jac
    sll_real64 :: coef_int_deta1, coef_int_deta2
    sll_real64 :: int_J_deta1deta2, int_rho_J_deta1deta2

    Neta1_loc  = size(fdistribu%val4d_seqx3x4,1)
    Neta2_loc  = size(fdistribu%val4d_seqx3x4,2)
    Nvx        = size(fdistribu%val4d_seqx3x4,3)
    Nvy        = size(fdistribu%val4d_seqx3x4,4)

    !*** Computation of the RHS locally in (x1,x2) directions ***
    do iloc2 = 1,Neta2_loc
      do iloc1 = 1,Neta1_loc
        intf_dvxdvy = 0._f64
        !\todo improve integral computation (see reduction used by Aurore)
        do i3 = 1,Nvx
          do i4 = 1,Nvy
            dv_space = mesh4d%coef_int_dvx(i3) * &
                mesh4d%coef_int_dvy(i4)
            glob_ind4d(:) = sll_o_local_to_global( &
                fdistribu%layout4d_seqx3x4, &
                (/iloc1,iloc2,i3,i4/))
            i1 = glob_ind4d(1)
            i2 = glob_ind4d(2)
!VG!            delta_f = fdistribu%val4d_seqx3x4(iloc1,iloc2,i3,i4) - &
!VG!                equilibrium%feq_vxvy(i3,i4)
            delta_f = fdistribu%val4d_seqx3x4(iloc1,iloc2,i3,i4)
            intf_dvxdvy = intf_dvxdvy + delta_f*dv_space
          end do
        end do
        QNsolver%QN_RHS%val2d_parx1x2(iloc1,iloc2) = intf_dvxdvy
      end do
    end do

    !--> Compute rho(x,y) = \int f(x,y,vx,vy) dvx dvy
    !-->   (sequential in (x1,x2)=(x,y)
    call remap_parx1x2_to_seqx1x2_field2d_VP4D(QNsolver%QN_RHS)

    !--> Compute \int rho(eta1,eta2) |J| deta1 deta2 / \int |J| deta1 deta2
    Neta1 = size(QNsolver%QN_RHS%val2d_seqx1x2,1)
    Neta2 = size(QNsolver%QN_RHS%val2d_seqx1x2,2)
    int_J_deta1deta2     = 0._f64
    int_rho_J_deta1deta2 = 0._f64
    do i2 = 1,Neta2 
      coef_int_deta2 = mesh4d%coef_int_deta2(i2)
      do i1 = 1,Neta1
        jacob_tmp      = mesh4d%transf_eta1eta2_xy%jacobian_at_node(i1,i2)
        val_jac        = abs(jacob_tmp)
        coef_int_deta1 = mesh4d%coef_int_deta1(i1)

        int_J_deta1deta2     = int_J_deta1deta2 + &
            val_jac*coef_int_deta1*coef_int_deta2
        int_rho_J_deta1deta2 = int_rho_J_deta1deta2 + &
            QNsolver%QN_RHS%val2d_seqx1x2(i1,i2) * &
            val_jac*coef_int_deta1*coef_int_deta2
      end do
    end do
    int_rho_J_deta1deta2 = int_rho_J_deta1deta2/int_J_deta1deta2

    do i2 = 1,Neta2
      do i1 = 1,Neta1
        QNsolver%QN_RHS%val2d_seqx1x2(i1,i2) = &
            QNsolver%QN_RHS%val2d_seqx1x2(i1,i2)-int_rho_J_deta1deta2
      end do
    end do

  end subroutine compRHS_QNsolver_Bsplines_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Solving
  !---------------------------------------------------------------------------
  subroutine solve_QNsolver_Bsplines_VP4D( &
      QNsolver, &
      mesh4d, &
      equilibrium, &
      fdistribu, &
      Phi )

    type(QNsolver_Bsplines_VP4D_t), intent(inout) :: QNsolver
    type(mesh_VP4D_t)             , intent(in)    :: mesh4d
    type(equilibrium_VP4D_t)      , intent(in)    :: equilibrium
    type(fdistribu_VP4D_t)        , intent(in)    :: fdistribu
    type(field2d_VP4D_t)          , intent(inout) :: Phi

    sll_int32 :: ieta1, ieta2
    sll_int32 :: Neta1, Neta2
    !sll_int32 :: loc2d_sz_x1, loc2d_sz_x2

    !*** Computation of the RHS of the quasi-neutrality equation ***
    call compRHS_QNsolver_Bsplines_VP4D( QNsolver, &
      mesh4d, &
      equilibrium, &
      fdistribu )

    !*** Solving of the quasi-neutrality equation ***
    Neta1 = QNsolver%Neta1
    Neta2 = QNsolver%Neta2

    call QNsolver%QNrhs2d%set_field_data( &
        QNsolver%QN_RHS%val2d_seqx1x2(:,:) )
    call QNsolver%QNrhs2d%update_interpolation_coefficients( )
    call sll_o_solve( &
        QNsolver%QNS, &
        QNsolver%QNrhs2d, &
        QNsolver%Phi2d )
    do ieta2 = 1, Neta2
      do ieta1 = 1, Neta1 
        Phi%val2d_seqx1x2(ieta1,ieta2) = &
            QNsolver%Phi2d%value_at_indices(ieta1,ieta2)
      end do
    end do
    
  end subroutine solve_QNsolver_Bsplines_VP4D
  !---------------------------------------------------------------------------

end module QNsolver_Bsplines_VP4D_module
!---------------------------------------------------------------------------
