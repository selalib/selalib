!===========================================================================
!> Quasi-neutrality solver for
!>  4D drift-kinetic hybrid simulation
!> 
!> Rk1: Use the polar poisson solver
!> 
!> Rk2: Need profiles and magnetic field in r 
!>
!> \date 2014-08-20
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module QNsolver_DK4D_polar_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use equilibrium_DK4D_module
  use fdistribu_DK4D_module
  use field3d_DK4D_module
  use mesh_DK4D_module
  use sll_m_poisson_2d_base
  use sll_m_poisson_2d_polar
  use sll_m_qn_solver_2d_polar
  use utils_DK4D_module

  
  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: QNsolver_DK4D_polar_t

    !> Boundary conditions in (eta1,eta2)
    type(boundary_conditions_2d_t) :: bound_cond
    !> Neumann or Dirichlet boundary conditions
    sll_int32, dimension(2) :: QNS_BC

    !> For interpolations
    type(spline_degree_2d_t) :: spline_degree

    !> 2d poisson solver
    type(sll_t_qn_solver_2d_polar), pointer :: QNS

    !> Number of points in (eta1,eta2) direction
    sll_int32 :: Neta1
    sll_int32 :: Neta2

    !> RHS term 
    type(field3D_DK4D_t) :: QN_RHS

    !> Terms used in the LHS of the solver      
    sll_real64, dimension(:), pointer :: dlog_density_r
    sll_real64, dimension(:), pointer :: inv_Te_r

  end type QNsolver_DK4D_polar_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Quasi-Neutrality solver: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_QNsolver_DK4D_polar( &
      QNsolver, &
      mesh4d, &
      bound_cond, &
      spline_degree )

    type(QNsolver_DK4D_polar_t)   , intent(inout) :: QNsolver
    type(mesh_DK4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_3d_t), intent(in)    :: bound_cond
    type(spline_degree_3d_t)      , intent(in)    :: spline_degree

    sll_int32 :: ierr

    !*** Initialization of the boundary conditions ***
    QNsolver%bound_cond%left_eta1  = bound_cond%left_eta1
    QNsolver%bound_cond%right_eta1 = bound_cond%right_eta1
    QNsolver%bound_cond%left_eta2  = bound_cond%left_eta2
    QNsolver%bound_cond%right_eta2 = bound_cond%right_eta2

!\todo Let the possibility to choose between Dirichlet and Neumann
    QNsolver%QNS_BC(1) = sll_p_dirichlet
    QNsolver%QNS_BC(2) = sll_p_neumann_mode_0
!\end todo    

    !*** Initialization for interpolations ***
    QNsolver%spline_degree%eta1 = spline_degree%eta1
    QNsolver%spline_degree%eta2 = spline_degree%eta2

    !*** Initialization of the number of points   ***
    !***   in (eta1,eta2) direction               ***
    QNsolver%Neta1 = size(mesh4d%eta1_grid)
    QNsolver%Neta2 = size(mesh4d%eta2_grid)

    !*** Allocation for the LHS of the solver ***
    SLL_ALLOCATE( QNsolver%dlog_density_r(QNsolver%Neta1), ierr )
    SLL_ALLOCATE( QNsolver%inv_Te_r(QNsolver%Neta1), ierr )

    !*** Initialization of the 3D field for RHS term ***
    call new_field3D_DK4D( &
        QNsolver%QN_RHS, &
        'QN_RHS', &
        mesh4d, &
        bound_cond, &
        spline_degree )

  end subroutine new_QNsolver_DK4D_polar
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_QNsolver_DK4D_polar( QNsolver )

    type(QNsolver_DK4D_polar_t), intent(inout) :: QNsolver

    sll_int32 :: ierr
    
    SLL_DEALLOCATE( QNsolver%dlog_density_r, ierr )
    SLL_DEALLOCATE( QNsolver%inv_Te_r, ierr )

  end subroutine delete_QNsolver_DK4D_polar
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Initialization
  !---------------------------------------------------------------------------
  subroutine init_QNsolver_DK4D_polar( QNsolver, &
      mesh4d, &
      equilibrium )


    type(QNsolver_DK4D_polar_t), intent(inout) :: QNsolver
    type(mesh_DK4D_t)          , intent(in)    :: mesh4d
    type(equilibrium_DK4D_t)   , intent(in)    :: equilibrium

    sll_int32  :: ir, Nr
    sll_int32  :: num_cells1 , num_cells2
    sll_real64 :: eta1_min, eta1_max
    sll_real64 :: r_point
    
    num_cells1 = size(mesh4d%eta1_grid) - 1
    num_cells2 = size(mesh4d%eta2_grid) - 1
    Nr       = num_cells1 + 1
    eta1_min = mesh4d%eta1_grid(1)
    eta1_max = mesh4d%eta1_grid(Nr)

    do ir = 1,Nr
      r_point = mesh4d%eta1_grid(ir)
      QNsolver%dlog_density_r(ir) = - equilibrium%inv_Ln * &
          cosh((r_point-equilibrium%r_peak)/equilibrium%deltarn)**(-2)
      QNsolver%inv_Te_r(ir) = 1._f64/equilibrium%Te_r(ir)
    end do

    !*** Allocation and initialization of the QNS type ***
    !--> Verification of the number of cells in theta direction
    !-->  which must be a power of 2
    if ( modulo(log(real(num_cells2, f64)),real(log(2._F64), f64)) .ne. 0 ) then
      print*,'Problem with the number of cells in theta direction:', &
          num_cells2, ' is not a power of 2 as required for SELALIB fft'
      stop
    end if
!!$    QNsolver%QNS =>sll_f_new_poisson_2d_polar( &
!!$        eta1_min, &
!!$        eta1_max, &
!!$        num_cells1, &
!!$        num_cells2, &
!!$        QNsolver%QNS_BC, &
!!$        dlog_density=QNsolver%dlog_density_r, &
!!$        inv_Te=QNsolver%inv_Te_r, &
!!$        poisson_case=sll_p_poisson_drift_kinetic)

    call sll_s_qn_solver_2d_polar_init(QNsolver%QNS, eta1_min, eta1_max, num_cells1, num_cells2, &
        QNsolver%QNS_BC(1), &
        QNsolver%QNS_BC(2), &
        equilibrium%n0_r, &! rho_m0
        0.0_f64*QNsolver%dlog_density_r + 1.0_f64, &! b_magn !TODO
        1/sqrt(QNsolver%inv_Te_r))! Lambda
        
         

  end subroutine init_QNsolver_DK4D_polar
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
  subroutine compRHS_QNsolver_DK4D_polar( QNsolver, &
      mesh4d, &
      equilibrium, &
      fdistribu )

    use sll_m_reduction, only : sll_s_compute_reduction_4d_to_3d_direction4

    type(QNsolver_DK4D_polar_t), intent(inout) :: QNsolver
    type(mesh_DK4D_t)          , intent(in)    :: mesh4d
    type(equilibrium_DK4D_t)   , intent(in)    :: equilibrium
    type(fdistribu_DK4D_t)     , intent(in)    :: fdistribu

    sll_int32  :: Neta1_loc, Neta2_loc
    sll_int32  :: Neta3, Nvpar
    sll_int32  :: iloc1, iloc2
    sll_int32  :: i1, i2, i3, i4
    sll_int32  :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3
    !sll_int32  :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_real64 :: delta_f
    sll_real64 :: delta_vpar, intf_dvpar
    sll_int32, dimension(1:4) :: glob_ind4d


    Neta1_loc  = size(fdistribu%val4d_seqx3x4,1)
    Neta2_loc  = size(fdistribu%val4d_seqx3x4,2)
    Neta3      = size(fdistribu%val4d_seqx3x4,3)
    Nvpar      = size(fdistribu%val4d_seqx3x4,4)
    delta_vpar = mesh4d%vpar_mesh1d%delta_eta

    call sll_o_compute_local_sizes( &
        QNsolver%QN_RHS%layout3d_seqx3, &
        loc3d_sz_x1, &
        loc3d_sz_x2, &
        loc3d_sz_x3 )        

    SLL_ASSERT( loc3d_sz_x1.eq.Neta1_loc )
    SLL_ASSERT( loc3d_sz_x2.eq.Neta2_loc )
    SLL_ASSERT( loc3d_sz_x3.eq.Neta3 )

    !*** Computation of the RHS locally sequential in x3 directions ***
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

    do iloc2 = 1,Neta2_loc
      do iloc1 = 1,Neta1_loc
        glob_ind4d(:) = sll_o_local_to_global( &
            fdistribu%layout4d_seqx3x4, &
            (/iloc1,iloc2,i3,1/))
        i1 = glob_ind4d(1)
        i2 = glob_ind4d(2)
        QNsolver%QN_RHS%val3d_seqx3(iloc1,iloc2,:) = &
            QNsolver%QN_RHS%val3d_seqx3(iloc1,iloc2,:) / &
            equilibrium%n0_xy(i1,i2)
      end do
    end do

    !--> compute rho3d_seqx1x2
    call sll_o_apply_remap_3d( &
        QNsolver%QN_RHS%seqx3_to_seqx1x2, &
        QNsolver%QN_RHS%val3d_seqx3, &
        QNsolver%QN_RHS%val3d_seqx1x2 )

  end subroutine compRHS_QNsolver_DK4D_polar
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Quasi-Neutrality solver: Solving
  !---------------------------------------------------------------------------
  subroutine solve_QNsolver_DK4D_polar( &
      QNsolver, &
      mesh4d, &
      equilibrium, &
      fdistribu, &
      Phi )

    type(QNsolver_DK4D_polar_t), intent(inout) :: QNsolver
    type(mesh_DK4D_t)          , intent(in)    :: mesh4d
    type(equilibrium_DK4D_t)   , intent(in)    :: equilibrium
    type(fdistribu_DK4D_t)     , intent(in)    :: fdistribu
    type(field3d_DK4D_t)       , intent(inout) :: Phi

    sll_int32 :: iloc3
    sll_int32 :: Neta1, Neta2
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3

    !*** Computation of the RHS of the quasi-neutrality equation ***
    call compRHS_QNsolver_DK4D_polar( QNsolver, &
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

!baremettre \todo Comprendre pourquoi on ne peut pas mettre comme input QN_RHS
!VG!    do iloc3 = 1,loc3d_sz_x3
!VG!      call QNsolver%QNS%compute_phi_from_rho( &
!VG!          QNsolver%QN_RHS%val3d_seqx1x2(:,:,iloc3), &
!VG!          Phi%val3d_seqx1x2(:,:,iloc3) )
!VG!    end do

    do iloc3 = 1,loc3d_sz_x3
      Phi%val3d_seqx1x2(:,:,iloc3) = &
           QNsolver%QN_RHS%val3d_seqx1x2(:,:,iloc3)
      call sll_s_qn_solver_2d_polar_solve( QNsolver%QNS, &
          Phi%val3d_seqx1x2(:,:,iloc3), &
          Phi%val3d_seqx1x2(:,:,iloc3) )
    end do
!earemettre
    
    !*** Fill Phi%val3d_seqx3 ***
    call sll_o_apply_remap_3d(Phi%seqx1x2_to_seqx3, &
        Phi%val3d_seqx1x2, Phi%val3d_seqx3)    

  end subroutine solve_QNsolver_DK4D_polar
  !---------------------------------------------------------------------------

end module QNsolver_DK4D_polar_module
!---------------------------------------------------------------------------
