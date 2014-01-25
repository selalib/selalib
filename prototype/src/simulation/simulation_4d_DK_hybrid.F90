! The idea of this simulation is to merge the functionalities of the qns-based
! simulation and the older, cartesian 4d simulation for the purposes of 
! debugging/understanding the behavior of the QNS one. Once this objective is
! fulfilled, this simulation can be deleted.

module sll_simulation_4d_DK_hybrid_module

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"

  use sll_collective
  use sll_remapper
  use sll_simulation_base
  use sll_logical_meshes
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_fdistribu4D_DK
  use sll_general_coordinate_elliptic_solver_module
  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative
  use sll_arbitrary_degree_spline_interpolator_1d_module
  use sll_module_scalar_field_1d_base
  use sll_module_scalar_field_1d_alternative
  use sll_timer
  use sll_module_deboor_splines_1d

  implicit none

#define PRINT_PLOTS 1
  type, extends(sll_simulation_base_class) :: sll_simulation_4d_DK_hybrid
    ! Parallel environment parameters
    sll_int32  :: world_size
    sll_int32  :: my_rank
    sll_int32  :: power2 ! 2^power2 = number of processes available
    ! Mesh parameters
    sll_int32  :: nc_x1
    sll_int32  :: nc_x2
    sll_int32  :: nc_x3
    sll_int32  :: nc_x4
    sll_real64 :: r_min
    sll_real64 :: r_max
    sll_real64 :: phi_min
    sll_real64 :: phi_max
    sll_real64 :: vpar_min
    sll_real64 :: vpar_max
    !--> Equilibrium
    sll_real64 :: tau0      !-> tau0 = Ti(rpeak)/Te(rpeak)
    sll_real64 :: rho_peak    
    sll_real64 :: kappan   
    sll_real64 :: deltarn  
    sll_real64 :: kappaTi  
    sll_real64 :: deltarTi 
    sll_real64 :: kappaTe  
    sll_real64 :: deltarTe     
    !--> Pertubation
    sll_int32  :: perturb_choice
    sll_int32  :: mmode
    sll_int32  :: nmode
    sll_real64 :: eps_perturb   
    ! Physics/numerical parameters
    sll_real64 :: dt
    sll_int32  :: nb_iter
    sll_int32  :: spline_degree_eta1
    sll_int32  :: spline_degree_eta2
    sll_int32  :: spline_degree_eta3
    sll_int32  :: spline_degree_vpar

    ! diagnostics
    sll_real64 :: diag2D_step

    !--> diagnostics for time
    sll_real64 :: iter_time
    sll_real64, dimension(:),pointer :: time_evol

    !--> diagnostics norm
    sll_real64, dimension(:),pointer :: diag_mass
    sll_real64, dimension(:),pointer :: diag_norm_L1
    sll_real64, dimension(:),pointer :: diag_norm_L2
    sll_real64, dimension(:),pointer :: diag_norm_Linf
    sll_real64, dimension(:),pointer :: diag_entropy_kin

    !--> diagnostics energy
    sll_real64, dimension(:),pointer :: diag_nrj_kin
    sll_real64, dimension(:),pointer :: diag_nrj_pot
    sll_real64, dimension(:),pointer :: diag_nrj_tot
    sll_real64, dimension(:),pointer :: diag_heat_flux

    !--> 4D logical mesh (eta1,eta2,eta3,vpar)
    sll_int32 :: Neta1, Neta2, Neta3, Nvpar
    type(sll_logical_mesh_4d), pointer :: logical_mesh4d
    sll_real64, dimension(:), pointer :: eta1_grid
    sll_real64, dimension(:), pointer :: eta2_grid
    sll_real64, dimension(:), pointer :: eta3_grid
    sll_real64, dimension(:), pointer :: vpar_grid

    !--> Coordinate transformation F
    class(sll_coordinate_transformation_2d_base), pointer :: transf_xy
    !--> Norm of norm**2 = (F_1(eta1, eta2)**2 + F_2(eta1, eta2)**2)
    sll_real64, dimension(:,:), pointer :: norm_square_xy

    !-->  1D mesh
    sll_real64, dimension(:), pointer :: r_grid
    !--> 2D physic mesh
    sll_real64, dimension(:,:), pointer :: xgrid_2d
    sll_real64, dimension(:,:), pointer :: ygrid_2d

    !--> For boundary conditions
    sll_int32 :: bc_left_eta1
    sll_int32 :: bc_right_eta1
    sll_int32 :: bc_left_eta2
    sll_int32 :: bc_right_eta2
    sll_int32 :: bc_left_eta3
    sll_int32 :: bc_right_eta3
    sll_int32 :: bc_left_vpar
    sll_int32 :: bc_right_vpar

    !--> Density and temperature profiles
    sll_real64 :: major_radius    !R0
    sll_real64, dimension(:)  , pointer :: n0_r
    sll_real64, dimension(:)  , pointer :: Ti_r
    sll_real64, dimension(:)  , pointer :: Te_r
    sll_real64, dimension(:,:), pointer :: n0_xy
    sll_real64, dimension(:,:), pointer :: Ti_xy
    sll_real64, dimension(:,:), pointer :: Te_xy

    !--> Magnetic field
    sll_real64, dimension(:,:), pointer :: B_xy

    !--> Equilibrium distribution function
    sll_real64, dimension(:,:,:), pointer :: feq_xyvpar

    !--> 4D distribution function 
    !----> sequential in (x1,x2) and parallel in (x3,x4)
    type(layout_4D), pointer :: layout4d_seqx1x2
    sll_real64, dimension(:,:,:,:), pointer :: f4d_seqx1x2 
    !----> parallel in (x1,x2) and sequential in (x3,x4) 
    type(layout_4D), pointer :: layout4d_seqx3x4
    sll_real64, dimension(:,:,:,:), pointer :: f4d_seqx3x4
    !----> for remapping
    type(remap_plan_4D_real64), pointer :: seqx1x2_to_seqx3x4
    type(remap_plan_4D_real64), pointer :: seqx3x4_to_seqx1x2
    !----> for interpolations
    type(arb_deg_2d_interpolator) :: interp2d_f_eta1eta2
    type(arb_deg_1d_interpolator) :: interp1d_f_eta3
    type(arb_deg_1d_interpolator) :: interp1d_f_vpar

    !--> 3D charge density and 3D electric potential
    !----> sequential in (x1,x2)
    type(layout_3D), pointer :: layout3d_seqx1x2
    sll_real64, dimension(:,:,:), pointer :: rho3d_seqx1x2 
    sll_real64, dimension(:,:,:), pointer :: phi3d_seqx1x2 
    !----> sequential in x3
    type(layout_3D), pointer :: layout3d_seqx3
    sll_real64, dimension(:,:,:), pointer :: rho3d_seqx3
    sll_real64, dimension(:,:,:), pointer :: phi3d_seqx3
    !----> for remapping
    type(remap_plan_3D_real64), pointer :: seqx1x2_to_seqx3
    type(remap_plan_3D_real64), pointer :: seqx3_to_seqx1x2

    !--> 3D electric field
    sll_real64, dimension(:,:,:), pointer :: E3d_eta1_seqx1x2
    sll_real64, dimension(:,:,:), pointer :: E3d_eta2_seqx1x2
    sll_real64, dimension(:,:,:), pointer :: E3d_x1_seqx1x2
    sll_real64, dimension(:,:,:), pointer :: E3d_x2_seqx1x2
    sll_real64, dimension(:,:,:), pointer :: E3d_eta3_seqx3

    !--> For general QN solver
    type(general_coordinate_elliptic_solver), pointer :: QNS
    ! interpolation any arbitrary spline
    type(arb_deg_2d_interpolator) :: interp2d_rho_eta1eta2
    type(arb_deg_2d_interpolator) :: interp2d_Phi_eta1eta2
    type(arb_deg_1d_interpolator) :: interp1d_Phi_eta3
    type(arb_deg_2d_interpolator) :: interp2d_QN_A11
    type(arb_deg_2d_interpolator) :: interp2d_QN_A12
    type(arb_deg_2d_interpolator) :: interp2d_QN_A21
    type(arb_deg_2d_interpolator) :: interp2d_QN_A22
    type(arb_deg_2d_interpolator) :: interp2d_QN_B1
    type(arb_deg_2d_interpolator) :: interp2d_QN_B2
    type(arb_deg_2d_interpolator) :: interp2d_QN_C
    class(sll_scalar_field_2d_base) , pointer :: rho2d
    type(sll_scalar_field_1d_discrete_alt), pointer :: phi1d! for derivative in eta3
    type(sll_scalar_field_2d_discrete_alt), pointer :: phi2d
    class(sll_scalar_field_2d_base), pointer :: QN_A11 
    class(sll_scalar_field_2d_base), pointer :: QN_A12
    class(sll_scalar_field_2d_base), pointer :: QN_A21
    class(sll_scalar_field_2d_base), pointer :: QN_A22
    class(sll_scalar_field_2d_base), pointer :: QN_B1
    class(sll_scalar_field_2d_base), pointer :: QN_B2
    class(sll_scalar_field_2d_base), pointer :: QN_C

    !---> For diagnostic saving
    sll_int32 :: count_save_diag

  contains
    procedure, pass(sim) :: run => run_4d_DK_hybrid
    procedure, pass(sim) :: init_from_file => init_4d_DK_hybrid
  end type sll_simulation_4d_DK_hybrid

  interface sll_delete
    module procedure delete_4d_DK_hybrid
  end interface sll_delete

  interface initialize
    module procedure initialize_4d_DK_hybrid
  end interface initialize

contains

  !************************************************************
  !  4D DRIFT-KINETIC HYBRID SIMULATION
  !************************************************************

  !----------------------------------------------------
  ! Initialization of the drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine initialize_4d_DK_hybrid( sim, &
      world_size, &
      my_rank, &
      logical_mesh4d, &
      transf_xy)

    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    sll_int32                        , intent(in)    :: world_size
    sll_int32                        , intent(in)    :: my_rank
    type(sll_logical_mesh_4d)        , pointer       :: logical_mesh4d
    class(sll_coordinate_transformation_2d_base), pointer :: transf_xy

    sll_int32 :: ierr
    sll_int32 :: ieta1, ieta2, ieta3, ivpar
    sll_int32 :: Nx, Ny
    sll_int32 :: nb_diag

    !--> Parallelization initialization
    sim%world_size = world_size
    sim%my_rank    = my_rank

    !--> Initialization of the number of points
    sim%Neta1 = sim%nc_x1+1
    sim%Neta2 = sim%nc_x2+1
    sim%Neta3 = sim%nc_x3+1
    sim%Nvpar = sim%nc_x4+1

    !--> Initialization of the boundary conditions
    sim%bc_left_eta1  = SLL_DIRICHLET
    sim%bc_right_eta1 = SLL_DIRICHLET
    sim%bc_left_eta2  = SLL_PERIODIC
    sim%bc_right_eta2 = SLL_PERIODIC
    sim%bc_left_eta3  = SLL_PERIODIC
    sim%bc_right_eta3 = SLL_PERIODIC
    sim%bc_left_vpar  = SLL_DIRICHLET
    sim%bc_right_vpar = SLL_DIRICHLET

    !--> Logical mesh initialization
    sim%logical_mesh4d => logical_mesh4d
    SLL_ALLOCATE(sim%eta1_grid(sim%Neta1),ierr)    
    SLL_ALLOCATE(sim%eta2_grid(sim%Neta2),ierr)
    SLL_ALLOCATE(sim%eta3_grid(sim%Neta3),ierr)
    SLL_ALLOCATE(sim%vpar_grid(sim%Nvpar),ierr)
    do ieta1 = 1,sim%Neta1
      sim%eta1_grid(ieta1) = sim%logical_mesh4d%eta1_min + &
          (ieta1-1)*sim%logical_mesh4d%delta_eta1
    end do
    do ieta2 = 1,sim%Neta2
      sim%eta2_grid(ieta2) = sim%logical_mesh4d%eta2_min + &
          (ieta2-1)*sim%logical_mesh4d%delta_eta2
    end do
    do ieta3 = 1,sim%Neta3
      sim%eta3_grid(ieta3) = sim%logical_mesh4d%eta3_min + &
          (ieta3-1)*sim%logical_mesh4d%delta_eta3
    end do
    do ivpar = 1,sim%Nvpar
      sim%vpar_grid(ivpar) = sim%logical_mesh4d%eta4_min + &
          (ivpar-1)*sim%logical_mesh4d%delta_eta4
    end do

    !--> Transformation initialization
    sim%transf_xy => transf_xy

    !--> Initialization of the (x,y) 2D mesh and
    !--> Initialization of the square of the norm 
    Nx = sim%Neta1
    Ny = sim%Neta2
    SLL_ALLOCATE(sim%xgrid_2d(Nx,Ny),ierr)
    SLL_ALLOCATE(sim%ygrid_2d(Nx,Ny),ierr)
    SLL_ALLOCATE(sim%norm_square_xy(Nx,Ny),ierr)
    do ieta2 = 1,sim%Neta2
      do ieta1 = 1,sim%Neta1
        sim%xgrid_2d(ieta1,ieta2) = &
            sim%transf_xy%x1_at_node(ieta1,ieta2)
        sim%ygrid_2d(ieta1,ieta2) = &
            sim%transf_xy%x2_at_node(ieta1,ieta2)
        sim%norm_square_xy(ieta1,ieta2) = sim%xgrid_2d(ieta1,ieta2)**2 &
            + sim%ygrid_2d(ieta1,ieta2)**2
      end do
    end do

    !--> Initialization diagnostics for the norm
    nb_diag  = int(sim%nb_iter*sim%dt/sim%diag2D_step) + 1

    SLL_ALLOCATE(sim%time_evol(nb_diag),ierr)

    SLL_ALLOCATE(sim%diag_mass(nb_diag),ierr)
    SLL_ALLOCATE(sim%diag_norm_L1(nb_diag),ierr)
    SLL_ALLOCATE(sim%diag_norm_L2(nb_diag),ierr)
    SLL_ALLOCATE(sim%diag_norm_Linf(nb_diag),ierr)
    SLL_ALLOCATE(sim%diag_entropy_kin(nb_diag),ierr)

    !--> Initialization diagnostics for the energy
    SLL_ALLOCATE(sim%diag_nrj_kin(nb_diag),ierr)
    SLL_ALLOCATE(sim%diag_nrj_pot(nb_diag),ierr)
    SLL_ALLOCATE(sim%diag_nrj_tot(nb_diag),ierr)
    SLL_ALLOCATE(sim%diag_heat_flux(nb_diag),ierr)

    !--> Radial profile initialisation
    call init_profiles_DK(sim)

    !*** Allocation of the distribution function ***
    call allocate_fdistribu4d_DK(sim)

    !*** Allocation of the QN solveR ***
    call allocate_QN_DK(sim)

    !*** Initialization of the QN solver ***
    call initialize_QN_DK (sim)

  end subroutine initialize_4d_DK_hybrid


  !----------------------------------------------------
  ! Run drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine run_4d_DK_hybrid( sim )

    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32 :: iter, diag_num

    !--> For timers
    type(sll_time_mark) :: t0, t1 
    double precision    :: elaps_time_advec, elaps_time_QN

    elaps_time_advec = 0.0
    elaps_time_QN    = 0.0

    do iter = 1,sim%nb_iter
      sim%iter_time = sim%iter_time + sim%dt
      if (sim%my_rank.eq.0) &
          print*,' ===> ITERATION = ',iter, '/',sim%nb_iter,'==> ', &
          sim%iter_time

      call sll_set_time_mark(t0)
      !--> Advection in vpar direction'
      call advec1D_vpar(sim,0.5_f64*sim%dt)

      ! --> Advection in eta3 direction'
      call advec1D_eta3(sim,0.5_f64*sim%dt)

      !--> Sequential for the advection in eta1eta2
      call apply_remap_4D(sim%seqx3x4_to_seqx1x2, &
          sim%f4d_seqx3x4,sim%f4d_seqx1x2 )

      !--> Advection in eta1,eta2 direction'
      call advec2D_eta1eta2(sim,sim%dt)

      !--> Sequential for the advection in eta3 and in vpar
      call apply_remap_4D(sim%seqx1x2_to_seqx3x4, &
          sim%f4d_seqx1x2,sim%f4d_seqx3x4 )

      !--> Advection in eta3 direction'
      call advec1D_eta3(sim,0.5_f64*sim%dt)

      !--> Advection in vpar direction'
      call advec1D_vpar(sim,0.5_f64*sim%dt)

      !--> Sequential to solve the quasi-neutral equation
      call apply_remap_4D(sim%seqx3x4_to_seqx1x2, &
          sim%f4d_seqx3x4,sim%f4d_seqx1x2 )

      call sll_set_time_mark(t1)
      elaps_time_advec = elaps_time_advec + &
          sll_time_elapsed_between(t0,t1)

      !--> Solve the quasi-neutral equation
      call sll_set_time_mark(t0)
      !-----> Computation of the rhs of QN 
      call compute_charge_density(sim)
      !-----> Solve QN
      call solve_QN( sim )
      call sll_set_time_mark(t1)
      elaps_time_QN = elaps_time_QN + &
          sll_time_elapsed_between(t0,t1)

      !--> Compute the new electric field
      call compute_Efield( sim )

      !--> Sequential for the advection in eta3 and in vpar
      call apply_remap_4D(sim%seqx1x2_to_seqx3x4, &
          sim%f4d_seqx1x2,sim%f4d_seqx3x4)     

      !--> Save results in HDF5 files
      if ( mod(iter,int(sim%diag2D_step/sim%dt)) == 0) then
        sim%count_save_diag = sim%count_save_diag + 1
        diag_num            = sim%count_save_diag

        !--> Compute energy kinetic, potential and total
        call compute_energy(sim,diag_num)

        !--> Compute L1 norm, L2 norm, L infini norm
        call compute_norm_L1_L2_Linf(sim,diag_num)

        !--> Save diagnostics in HDF5 format
        call writeHDF5_cross_section_diag(sim,diag_num)
      end if
    end do

    !--> Save 'conservation_laws.h5' HDF5 file
    call writeHDF5_conservation_laws( sim )

    if (sim%my_rank.eq.0) then
      print*, ' Time for advec in run_4d_DK_hybrid  = ', elaps_time_advec
      print*, ' Time for QN in run_4d_DK_hybrid     = ', elaps_time_QN
    end if

  end subroutine run_4d_DK_hybrid


  !----------------------------------------------------
  ! Initialization of the drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine delete_4d_DK_hybrid( sim )

    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32 :: ierr

    SLL_DEALLOCATE(sim%eta1_grid,ierr)    
    SLL_DEALLOCATE(sim%eta2_grid,ierr)
    SLL_DEALLOCATE(sim%eta3_grid,ierr)
    SLL_DEALLOCATE(sim%vpar_grid,ierr)
    SLL_DEALLOCATE(sim%r_grid,ierr) 
    SLL_DEALLOCATE(sim%norm_square_xy,ierr)
    SLL_DEALLOCATE(sim%xgrid_2d,ierr)
    SLL_DEALLOCATE(sim%ygrid_2d,ierr) 
    SLL_DEALLOCATE(sim%n0_r,ierr)
    SLL_DEALLOCATE(sim%Ti_r,ierr)
    SLL_DEALLOCATE(sim%Te_r,ierr)
    SLL_DEALLOCATE(sim%n0_xy,ierr)
    SLL_DEALLOCATE(sim%Ti_xy,ierr)
    SLL_DEALLOCATE(sim%Te_xy,ierr)
    SLL_DEALLOCATE(sim%B_xy,ierr)
    SLL_DEALLOCATE(sim%feq_xyvpar,ierr)
    SLL_DEALLOCATE(sim%f4d_seqx1x2,ierr)
    SLL_DEALLOCATE(sim%f4d_seqx3x4,ierr)
    call sll_delete(sim%layout4d_seqx1x2)
    call sll_delete(sim%layout4d_seqx3x4)
    SLL_DEALLOCATE(sim%rho3d_seqx1x2,ierr)
    SLL_DEALLOCATE(sim%rho3d_seqx3,ierr)
    SLL_DEALLOCATE(sim%phi3d_seqx1x2,ierr)
    SLL_DEALLOCATE(sim%phi3d_seqx3,ierr)
    SLL_DEALLOCATE(sim%E3d_eta1_seqx1x2,ierr)
    SLL_DEALLOCATE(sim%E3d_eta2_seqx1x2,ierr)
    SLL_DEALLOCATE(sim%E3d_x1_seqx1x2,ierr)
    SLL_DEALLOCATE(sim%E3d_x2_seqx1x2,ierr)
    SLL_DEALLOCATE(sim%E3d_eta3_seqx3,ierr)
    SLL_DEALLOCATE(sim%time_evol,ierr)
    SLL_DEALLOCATE(sim%diag_mass,ierr)
    SLL_DEALLOCATE(sim%diag_norm_L1,ierr)
    SLL_DEALLOCATE(sim%diag_norm_L2,ierr)
    SLL_DEALLOCATE(sim%diag_norm_Linf,ierr)
    SLL_DEALLOCATE(sim%diag_entropy_kin,ierr)
    SLL_DEALLOCATE(sim%diag_nrj_kin,ierr)
    SLL_DEALLOCATE(sim%diag_nrj_pot,ierr)
    SLL_DEALLOCATE(sim%diag_nrj_tot,ierr)
    SLL_DEALLOCATE(sim%diag_heat_flux,ierr)
    call sll_delete(sim%layout3d_seqx1x2)
    call sll_delete(sim%layout3d_seqx3)

  end subroutine delete_4d_DK_hybrid



  !************************************************************
  !  PROFILE INITIALISATION
  !************************************************************

  !----------------------------------------------------------
  !  Read the initial data file : sim4d_DK_hybrid_input.txt
  !----------------------------------------------------------
  subroutine init_4d_DK_hybrid( sim, filename)
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    character(len=*)                  , intent(in)    :: filename

    sll_int32            :: IO_stat
    sll_int32, parameter :: input_file = 99

    !--> Mesh
    sll_int32  :: num_cells_x1
    sll_int32  :: num_cells_x2
    sll_int32  :: num_cells_x3
    sll_int32  :: num_cells_x4
    sll_real64 :: r_min
    sll_real64 :: r_max
    sll_real64 :: phi_min
    sll_real64 :: phi_max
    sll_real64 :: vpar_min
    sll_real64 :: vpar_max
    !--> Equilibrium
    sll_real64 :: tau0
    sll_real64 :: rho_peak    
    sll_real64 :: kappan   
    sll_real64 :: deltarn  
    sll_real64 :: kappaTi  
    sll_real64 :: deltarTi 
    sll_real64 :: kappaTe  
    sll_real64 :: deltarTe     
    !--> Pertubation
    sll_int32  :: perturb_choice
    sll_int32  :: mmode
    sll_int32  :: nmode
    sll_real64 :: eps_perturb   
    !--> Algorithm
    sll_real64 :: dt
    sll_int32  :: number_iterations
    sll_int32  :: spline_degree
    !--> Diagnostics
    sll_real64 :: diag2D_step

    namelist /mesh/ num_cells_x1, num_cells_x2, &
      num_cells_x3, num_cells_x4, &
      r_min, r_max, phi_min, phi_max, &
      vpar_min, vpar_max
    namelist /equilibrium/ tau0, rho_peak, kappan, deltarn, &
      kappaTi, deltarTi, kappaTe, deltarTe
    namelist /perturbation/ perturb_choice, mmode, nmode, eps_perturb
    namelist /sim_params/ dt, number_iterations, spline_degree
    namelist /diagnostics/ diag2D_step

    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_4d_DK_hybrid() failed to open file ', filename
       STOP
    end if
    read(input_file,mesh)
    read(input_file,equilibrium)
    read(input_file,perturbation)
    read(input_file,sim_params)
    read(input_file,diagnostics)
    close(input_file)

    !--> Mesh
    sim%nc_x1    = num_cells_x1
    sim%nc_x2    = num_cells_x2
    sim%nc_x3    = num_cells_x3
    sim%nc_x4    = num_cells_x4
    sim%r_min    = r_min
    sim%r_max    = r_max
    sim%phi_min  = phi_min
    sim%phi_max  = phi_max
    sim%vpar_min = vpar_min
    sim%vpar_max = vpar_max
    !--> Equilibrium
    sim%tau0     = tau0
    sim%rho_peak = rho_peak 
    sim%kappan   = kappan
    sim%deltarn  = deltarn
    sim%kappaTi  = kappaTi
    sim%deltarTi = deltarTi
    sim%kappaTe  = kappaTe
    sim%deltarTe = deltarTe
    !--> Pertubation
    sim%perturb_choice = perturb_choice
    sim%mmode          = mmode
    sim%nmode          = nmode
    sim%eps_perturb    = eps_perturb
    !--> Algorithm
    sim%dt                 = dt
    sim%nb_iter            = number_iterations
    sim%spline_degree_eta1 = spline_degree
    sim%spline_degree_eta2 = spline_degree
    sim%spline_degree_eta3 = spline_degree
    sim%spline_degree_vpar = spline_degree
    !--> Diagnostics
    sim%diag2D_step   = diag2D_step

  end subroutine init_4d_DK_hybrid


  !----------------------------------------------------
  ! Initialization of the logical mesh associated
  !  to the 4D DK simulation
  !----------------------------------------------------
  subroutine init_profiles_DK( sim )
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32  :: ierr
    sll_int32  :: ir, itheta
    sll_int32  :: Nr, Ntheta 
    sll_int32  :: Nx, Ny
    sll_int32  :: ieta1,ieta2
    sll_real64 :: Lr, dr
    sll_real64 :: x,y,r_point
    sll_real64 :: theta_min, Ltheta, dtheta
    sll_real64 :: r_peak, n0_rmin
    sll_real64 :: inv_Ln, inv_LTi, inv_LTe
    sll_real64 :: Ti_rmin, Te_rmin, Ti_scal, Te_scal
    sll_real64, dimension(:), pointer   :: theta_grid_tmp
    sll_real64, dimension(:,:), pointer :: B_rtheta_tmp

    sll_real64, dimension(3) :: params_n0,params_Te,params_Ti
    
    sll_real64 :: Lz

    !--> Initialisation of R0 = major_radius
    Lz               = abs(sim%phi_max-sim%phi_min)
    sim%major_radius = Lz/(2._f64*sll_pi)

    !--> Initialization of r_grid
    Nr = sim%nc_x1+1
    SLL_ALLOCATE(sim%r_grid(Nr),ierr)
    Lr = abs(sim%r_max-sim%r_min)
    dr = Lr/float(Nr-1)
    do ir = 1,Nr
      sim%r_grid(ir) = sim%r_min + float(ir-1)*dr
    end do

   ! !--> Initialization of n0(r), Ti(r) and Te(r)
    SLL_ALLOCATE(sim%n0_r(Nr),ierr)
    SLL_ALLOCATE(sim%Ti_r(Nr),ierr)
    SLL_ALLOCATE(sim%Te_r(Nr),ierr)
    sim%n0_r(:) = 0.0_f64
    sim%Ti_r(:) = 0.0_f64
    sim%Te_r(:) = 0.0_f64

    n0_rmin = 10._f64**19
    Ti_rmin = 1.3_f64
    Ti_scal = 1._f64
    
    inv_Ln  = sim%kappan/sim%major_radius
    inv_LTi = sim%kappaTi/sim%major_radius
    inv_LTe = sim%kappaTe/sim%major_radius

    r_peak  = sim%r_min + sim%rho_peak * Lr

  !  call init_n0_r(r_peak,inv_Ln, &
  !    sim%deltarn,n0_rmin,sim%r_grid,sim%n0_r)
  !  call init_T_r(r_peak,inv_LTi, &
  !    sim%deltarTi,Ti_rmin,Ti_scal,sim%r_grid,sim%Ti_r)
  !  Te_rmin = sim%Ti_r(1)
  !  Te_scal = sim%tau0
  !  call init_T_r(r_peak,inv_LTe, &
  !    sim%deltarTe,Te_rmin,Te_scal,sim%r_grid,sim%Te_r)


    !--> Initialization of n0(x,y), Ti(x,y) and Te(x,y)
    Nx = sim%Neta1
    Ny = sim%Neta2
    SLL_ALLOCATE(sim%n0_xy(Nx,Ny),ierr)
    SLL_ALLOCATE(sim%Ti_xy(Nx,Ny),ierr)
    SLL_ALLOCATE(sim%Te_xy(Nx,Ny),ierr)
   ! call function_xy_from_r(sim%r_grid,sim%n0_r,sim%xgrid_2d, &
   !      sim%ygrid_2d,sim%n0_xy)
   ! call function_xy_from_r(sim%r_grid,sim%Ti_r,sim%xgrid_2d, &
   !      sim%ygrid_2d,sim%Ti_xy)
   ! call function_xy_from_r(sim%r_grid,sim%Te_r,sim%xgrid_2d, &
   !      sim%ygrid_2d,sim%Te_xy)

    !--> parameters for radial density profiles
    params_n0(1) = inv_Ln
    params_n0(2) = sim%deltarn*Lr
    params_n0(3) = r_peak

    !--> parameters for radial temperature electrons profiles
    params_Te(1) = inv_LTe
    params_Te(2) = sim%deltarTe*Lr
    params_Te(3) = r_peak

    !--> parameter for radial temperature ions profiles
    params_Ti(1) = inv_LTi
    params_Ti(2) = sim%deltarTi*Lr
    params_Ti(3) = r_peak

    do ir = 1,Nr
      r_point = sim%r_grid(ir)
      sim%n0_r(ir) =  init_exact_profile_r(r_point,params_n0) 
      sim%Te_r(ir) =  init_exact_profile_r(r_point,params_Te)
      sim%Ti_r(ir) =  init_exact_profile_r(r_point,params_Ti)        
    end do

    do ieta1 = 1,Nx
      do ieta2 = 1,Ny 
        x = sim%xgrid_2d(ieta1,ieta2)
        y = sim%ygrid_2d(ieta1,ieta2)

        sim%n0_xy(ieta1,ieta2) = profil_xy_exacte(x,y,params_n0)
        sim%Ti_xy(ieta1,ieta2) = profil_xy_exacte(x,y,params_Ti) 
        sim%Te_xy(ieta1,ieta2) = profil_xy_exacte(x,y,params_Te)
      end do
    end do

    !--> Initialization of B(x,y)
    SLL_ALLOCATE(sim%B_xy(Nx,Ny),ierr)

    Ntheta = sim%nc_x2+1
    SLL_ALLOCATE(theta_grid_tmp(Ntheta),ierr)
    SLL_ALLOCATE(B_rtheta_tmp(Nr,Ntheta),ierr)

    theta_min = 0._f64
    Ltheta    = 2._f64*sll_pi
    dtheta    = Ltheta/float(Ntheta-1)
    do itheta = 1,Ntheta
      theta_grid_tmp(itheta) = theta_min + &
          float(itheta-1)*dtheta
    end do

    call init_Brtheta(sim%r_grid,theta_grid_tmp,B_rtheta_tmp)
    call function_xy_from_rtheta( &
      sim%r_grid,theta_grid_tmp, &
      B_rtheta_tmp,sim%xgrid_2d,sim%ygrid_2d, &
      sim%B_xy)

    SLL_DEALLOCATE(theta_grid_tmp,ierr)
    SLL_DEALLOCATE(B_rtheta_tmp,ierr)

  end subroutine init_profiles_DK


  !************************************************************
  !  DISTRIBUTION FUNCTION INITIALISATION
  !************************************************************

  !----------------------------------------------------
  ! Allocation of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine allocate_fdistribu4d_DK( sim )
#include "sll_assert.h"
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32 :: ierr, itemp
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_int32 :: nproc_x1
    sll_int32 :: nproc_x2
    sll_int32 :: nproc_x3
    sll_int32 :: nproc_x4 
   
    ! layout for sequential operations in x3 and x4. 
    ! Make an even split for x1 and x2, or as close as 
    ! even if the power of 2 is odd. This should 
    ! be packaged in some sort of routine and set up 
    ! at initialization time.
    sim%power2 = int(log(real(sim%world_size))/log(2.0))

    !--> special case N = 1, so power2 = 0
    if(sim%power2 == 0) then
       nproc_x1 = 1
       nproc_x2 = 1
       nproc_x3 = 1
       nproc_x4 = 1
    end if
    
    if(is_even(sim%power2)) then
       nproc_x1 = 2**(sim%power2/2)
       nproc_x2 = 2**(sim%power2/2)
       nproc_x3 = 1
       nproc_x4 = 1
    else 
       nproc_x1 = 2**((sim%power2-1)/2)
       nproc_x2 = 2**((sim%power2+1)/2)
       nproc_x3 = 1
       nproc_x4 = 1
    end if

    !--> Initialization of parallel layout of f4d in (x1,x2) directions
    !-->  (x1,x2) : parallelized layout
    !-->  (x3,x4) : sequential
    sim%layout4d_seqx3x4  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%nc_x1+1, &
      sim%nc_x2+1, &
      sim%nc_x3+1, &
      sim%nc_x4+1, &
      nproc_x1, &
      nproc_x2, &
      nproc_x3, &
      nproc_x4, &
      sim%layout4d_seqx3x4 )
        
    call compute_local_sizes_4d( sim%layout4d_seqx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    
    SLL_ALLOCATE(sim%f4d_seqx3x4(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)

    !--> Initialization of parallel layout of f4d in (x3,x4) directions
    !-->  (x1,x2) : sequential
    !-->  (x3,x4) : parallelized layout
!baremettre
!VG!    ! switch x1 and x3:
!VG!    itemp        = nproc_x3
!VG!    nproc_x3 = nproc_x1
!VG!    nproc_x1 = itemp
!VG!    ! switch x2 and x4
!VG!    itemp        = nproc_x4
!VG!    nproc_x4 = nproc_x2 
!VG!    nproc_x2 = itemp
    SLL_ASSERT((sim%nc_x3+1).ge.sim%power2)
    itemp    = nproc_x1*nproc_x2
    nproc_x3 = itemp
    nproc_x1 = 1
    nproc_x2 = 1
    nproc_x4 = 1
!earemettre
    sim%layout4d_seqx1x2  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%nc_x1+1, & 
      sim%nc_x2+1, & 
      sim%nc_x3+1, &
      sim%nc_x4+1, &
      nproc_x1, &
      nproc_x2, &
      nproc_x3, &
      nproc_x4, &
      sim%layout4d_seqx1x2 )
    
    ! Allocate the array needed to store the local chunk 
    ! of the distribution function data. First compute the 
    ! local sizes. Since the remap operations
    ! are out-of-place, we will allocate two different arrays, 
    ! one for each layout.
    call compute_local_sizes_4d( sim%layout4d_seqx1x2, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    SLL_ALLOCATE(sim%f4d_seqx1x2(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)

    !---> initialization for the remappings
    sim%seqx3x4_to_seqx1x2 => &
      NEW_REMAP_PLAN(sim%layout4d_seqx3x4,sim%layout4d_seqx1x2,sim%f4d_seqx3x4)
    
    sim%seqx1x2_to_seqx3x4 => &
      NEW_REMAP_PLAN(sim%layout4d_seqx1x2,sim%layout4d_seqx3x4,sim%f4d_seqx1x2)

     !----> for interpolations
    call sim%interp2d_f_eta1eta2%initialize( &
      sim%logical_mesh4d%num_cells1+1, &
      sim%logical_mesh4d%num_cells2+1, &
      sim%logical_mesh4d%eta1_min, &
      sim%logical_mesh4d%eta1_max, &
      sim%logical_mesh4d%eta2_min, &
      sim%logical_mesh4d%eta2_max, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      sim%spline_degree_eta1, &
      sim%spline_degree_eta2)

    call sim%interp1d_f_eta3%initialize( &
       sim%Neta3, &
       sim%logical_mesh4d%eta3_min, &
       sim%logical_mesh4d%eta3_max, &
       sim%bc_left_eta3, &
       sim%bc_right_eta3, &
       sim%spline_degree_eta3)

    call sim%interp1d_f_vpar%initialize( &
       sim%Nvpar, &
       sim%logical_mesh4d%eta4_min, &
       sim%logical_mesh4d%eta4_max, &
       sim%bc_left_vpar, &
       sim%bc_right_vpar, &
       sim%spline_degree_vpar)
   
  end subroutine allocate_fdistribu4d_DK


  !----------------------------------------------------
  ! Initialization of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine initialize_fdistribu4d_DK(sim)
    use sll_common_coordinate_transformations, only : &
      polar_eta2
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32  :: ierr
    sll_int32  :: i1, i2, i3, i4
    sll_int32  :: iloc1, iloc2, iloc3, iloc4
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_int32  :: Nx, Ny, Nphi, Nvpar
    sll_real64 :: theta_j, phi_k
    sll_int32, dimension(1:4) :: glob_ind4d

    sll_real64 :: dphi, Lphi, dvpar, Lvpar
    sll_real64, dimension(:), pointer :: phi_grid_tmp
    sll_real64, dimension(:), pointer :: vpar_grid_tmp
    sll_real64 :: tmp,r_peak
    Nx    = sim%Neta1
    Ny    = sim%Neta2
    Nphi  = sim%Neta3
    Nvpar = sim%Nvpar
    
    !--> Initialization of the grid in phi direction
    SLL_ALLOCATE(phi_grid_tmp(Nphi),ierr)
    Lphi = abs(sim%phi_max-sim%phi_min)
    dphi = Lphi/float(Nphi)
    do i3 = 1,Nphi
      phi_grid_tmp(i3) = sim%phi_min + &
        float(i3-1)*dphi
    end do

    !--> Initialization of the grid in vpar direction
    SLL_ALLOCATE(vpar_grid_tmp(Nvpar),ierr)
    Lvpar = abs(sim%vpar_max-sim%vpar_min)
    dvpar = Lvpar/float(Nvpar)
    do i4 = 1,Nvpar
      vpar_grid_tmp(i4) = sim%vpar_min + &
        float(i4-1)*dvpar
    end do

    !--> Initialization of the equilibrium distribution function
    SLL_ALLOCATE(sim%feq_xyvpar(Nx,Ny,Nvpar),ierr)
    call init_fequilibrium_xy(sim%xgrid_2d,sim%ygrid_2d, &
      vpar_grid_tmp,sim%n0_xy,sim%Ti_xy,sim%feq_xyvpar)

    !--> Initialization of the distribution function f4d_seqx3x4
    call compute_local_sizes_4d( sim%layout4d_seqx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )

    r_peak  = sim%r_min + sim%rho_peak * abs(sim%r_max-sim%r_min)

    do iloc4 = 1,loc4d_sz_x4
      do iloc3 = 1,loc4d_sz_x3
        do iloc2 = 1,loc4d_sz_x2
          do iloc1 = 1,loc4d_sz_x1
            glob_ind4d(:) = local_to_global_4D(sim%layout4d_seqx3x4, &
                (/iloc1,iloc2,iloc3,iloc4/))
            i1 = glob_ind4d(1)
            i2 = glob_ind4d(2)
            i3 = glob_ind4d(3)
            i4 = glob_ind4d(4)
            theta_j = sim%logical_mesh4d%eta2_min + &
                (i2-1)*sim%logical_mesh4d%delta_eta2 
            ! theta_j = polar_eta2( &
            !      sim%xgrid_2d(i1,i2), &
            !      sim%ygrid_2d(i1,i2), &
            !      (/0.0_f64/)) ! irrelevant for polar_eta2
            phi_k   = phi_grid_tmp(i3) 
            tmp = exp(-(sim%r_grid(i1)-r_peak)**2/&
                (4._f64*sim%deltarn/sim%deltarTi)) 
            sim%f4d_seqx3x4(iloc1,iloc2,i3,i4) = &
                sim%feq_xyvpar(i1,i2,i4) * &
                (1._f64 + sim%eps_perturb*&
                cos(2*sll_pi*real(sim%mmode)*theta_j + &
                2._f64*sll_pi*real(sim%nmode)*phi_k/Lphi) * tmp)
          end do
        end do
      end do
    end do

    call apply_remap_4D(sim%seqx3x4_to_seqx1x2,sim%f4d_seqx3x4,sim%f4d_seqx1x2)

    SLL_DEALLOCATE(phi_grid_tmp,ierr)
    SLL_DEALLOCATE(vpar_grid_tmp,ierr)
     print*, 'initialize fdistribution'

  end subroutine initialize_fdistribu4d_DK



  !************************************************************
  !  QUASI-NEUTRALITY SOLVING
  !************************************************************

  !----------------------------------------------------
  ! Allocation for QN solver
  !----------------------------------------------------
  subroutine allocate_QN_DK( sim )

    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    type(sll_logical_mesh_2d), pointer :: logical_mesh2d
    sll_int32 :: ierr
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3
    sll_int32 :: nproc_x1
    sll_int32 :: nproc_x2
    sll_int32 :: nproc_x3
    type(sll_logical_mesh_1d), pointer :: logical_mesh1d

    ! layout for sequential operations in (x1,x2) 
    sim%power2 = int(log(real(sim%world_size))/log(2.0))
    !--> special case N = 1, so power2 = 0
    if(sim%power2 == 0) then
       nproc_x1 = 1
       nproc_x2 = 1
       nproc_x3 = 1
    end if
    
    if(is_even(sim%power2)) then
       nproc_x1 = 2**(sim%power2/2)
       nproc_x2 = 2**(sim%power2/2)
       nproc_x3 = 1
    else 
       nproc_x1 = 2**((sim%power2-1)/2)
       nproc_x2 = 2**((sim%power2+1)/2)
       nproc_x3 = 1
    end if

    !--> Initialization of rho3d_seqx3 and phi3d_seqx3
    !-->  (x1,x2) : parallelized layout
    !-->  x3 : sequential        
    sim%layout3d_seqx3  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( &
      sim%nc_x1+1, &
      sim%nc_x2+1, &
      sim%nc_x3+1, &
      nproc_x1, &
      nproc_x2, &
      nproc_x3, &
      sim%layout3d_seqx3 )
    call compute_local_sizes_3d( sim%layout3d_seqx3, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE(sim%rho3d_seqx3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_seqx3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%E3d_eta3_seqx3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)

    !--> Initialization of rho3d_seqx1x2 and phi3d_seqx1x2
    !-->  (x1,x2) : sequential
    !-->  x3 : parallelized layout
    ! switch x3 and x1:
    nproc_x3 = nproc_x1*nproc_x2
    nproc_x1 = 1
    nproc_x2 = 1

    sim%layout3d_seqx1x2  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( &
      sim%nc_x1+1, & 
      sim%nc_x2+1, & 
      sim%nc_x3+1, &
      nproc_x1, &
      nproc_x2, &
      nproc_x3, &
      sim%layout3d_seqx1x2 )
    call compute_local_sizes_3d( sim%layout3d_seqx1x2, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE(sim%rho3d_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%E3d_eta1_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%E3d_eta2_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%E3d_x1_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%E3d_x2_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    
    !---->
    logical_mesh2d => sim%transf_xy%mesh

    !---> For iterpolations of Phi
    call sim%interp2d_Phi_eta1eta2%initialize( &
      logical_mesh2d%num_cells1+1, &
      logical_mesh2d%num_cells2+1, &
      logical_mesh2d%eta1_min, &
      logical_mesh2d%eta1_max, &
      logical_mesh2d%eta2_min, &
      logical_mesh2d%eta2_max, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      sim%spline_degree_eta1, &
      sim%spline_degree_eta2)

    call sim%interp1d_Phi_eta3%initialize( &
       sim%Neta3, &
       sim%logical_mesh4d%eta3_min, &
       sim%logical_mesh4d%eta3_max, &
       sim%bc_left_eta3, &
       sim%bc_right_eta3, &
       sim%spline_degree_eta3)

    !---> For rho
    call sim%interp2d_rho_eta1eta2%initialize( &
      logical_mesh2d%num_cells1 +1, &
      logical_mesh2d%num_cells2 +1, &
      logical_mesh2d%eta1_min, &
      logical_mesh2d%eta1_max, &
      logical_mesh2d%eta2_min, &
      logical_mesh2d%eta2_max, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      sim%spline_degree_eta1, &
      sim%spline_degree_eta2)    

    !---> For all the matrices required for QN solver
    call sim%interp2d_QN_A11%initialize( &
      logical_mesh2d%num_cells1 +1, &
      logical_mesh2d%num_cells2 +1, &
      logical_mesh2d%eta1_min, &
      logical_mesh2d%eta1_max, &
      logical_mesh2d%eta2_min, &
      logical_mesh2d%eta2_max, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      sim%spline_degree_eta1, &
      sim%spline_degree_eta2)    

    call sim%interp2d_QN_A12%initialize( &
      logical_mesh2d%num_cells1 +1, &
      logical_mesh2d%num_cells2 +1, &
      logical_mesh2d%eta1_min, &
      logical_mesh2d%eta1_max, &
      logical_mesh2d%eta2_min, &
      logical_mesh2d%eta2_max, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      sim%spline_degree_eta1, &
      sim%spline_degree_eta2)    

    call sim%interp2d_QN_A21%initialize( &
      logical_mesh2d%num_cells1 +1, &
      logical_mesh2d%num_cells2 +1, &
      logical_mesh2d%eta1_min, &
      logical_mesh2d%eta1_max, &
      logical_mesh2d%eta2_min, &
      logical_mesh2d%eta2_max, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      sim%spline_degree_eta1, &
      sim%spline_degree_eta2)    

    call sim%interp2d_QN_A22%initialize( &
      logical_mesh2d%num_cells1 +1, &
      logical_mesh2d%num_cells2 +1, &
      logical_mesh2d%eta1_min, &
      logical_mesh2d%eta1_max, &
      logical_mesh2d%eta2_min, &
      logical_mesh2d%eta2_max, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      sim%spline_degree_eta1, &
      sim%spline_degree_eta2)    

    call sim%interp2d_QN_B1%initialize( &
      logical_mesh2d%num_cells1 +1, &
      logical_mesh2d%num_cells2 +1, &
      logical_mesh2d%eta1_min, &
      logical_mesh2d%eta1_max, &
      logical_mesh2d%eta2_min, &
      logical_mesh2d%eta2_max, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      sim%spline_degree_eta1, &
      sim%spline_degree_eta2)
        
    call sim%interp2d_QN_B2%initialize( &
      logical_mesh2d%num_cells1 +1, &
      logical_mesh2d%num_cells2 +1, &
      logical_mesh2d%eta1_min, &
      logical_mesh2d%eta1_max, &
      logical_mesh2d%eta2_min, &
      logical_mesh2d%eta2_max, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      sim%spline_degree_eta1, &
      sim%spline_degree_eta2)
    
    call sim%interp2d_QN_C%initialize( &
      logical_mesh2d%num_cells1 +1, &
      logical_mesh2d%num_cells2 +1, &
      logical_mesh2d%eta1_min, &
      logical_mesh2d%eta1_max, &
      logical_mesh2d%eta2_min, &
      logical_mesh2d%eta2_max, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      sim%spline_degree_eta1, &
      sim%spline_degree_eta2)    

    !-----> rho2D field
    sim%rho2d => new_scalar_field_2d_discrete_alt( &
      "rho2d_seqx1x2", &
      sim%interp2d_rho_eta1eta2, &     
      sim%transf_xy, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2)
    !-----> phi1D in the direction eta3
    logical_mesh1d => new_logical_mesh_1d( &
      sim%nc_x3,eta_min=sim%phi_min,eta_max=sim%phi_max)
    sim%phi1d => new_scalar_field_1d_discrete_alt( &
         "phi1d_seqx3", &
         sim%interp1d_Phi_eta3, &     
         sim%bc_left_eta1, &
         sim%bc_right_eta1,&
         logical_mesh1d)
    !-----> phi2D in the direction eta1 eta2
    sim%phi2d => new_scalar_field_2d_discrete_alt( &
      "phi2d_seqx1x2", &
      sim%interp2d_Phi_eta1eta2, &     
      sim%transf_xy, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2)

    !---> initialization for the remappings
    sim%seqx3_to_seqx1x2 => &
      NEW_REMAP_PLAN(sim%layout3d_seqx3,sim%layout3d_seqx1x2,sim%rho3d_seqx3)
    
    sim%seqx1x2_to_seqx3 => &
      NEW_REMAP_PLAN(sim%layout3d_seqx1x2,sim%layout3d_seqx3,sim%rho3d_seqx1x2)
  end subroutine allocate_QN_DK


  !----------------------------------------------------
  ! Initialization of the QN coefficients
  !----------------------------------------------------
  subroutine initialize_QN_DK ( sim )

    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32 :: ierr
    sll_int32 :: Neta1, Neta2
    sll_real64, dimension(:,:), pointer :: A11
    sll_real64, dimension(:,:), pointer :: A12
    sll_real64, dimension(:,:), pointer :: A21
    sll_real64, dimension(:,:), pointer :: A22
    sll_real64, dimension(:,:), pointer :: B1
    sll_real64, dimension(:,:), pointer :: B2
    sll_real64, dimension(:,:), pointer :: C
    type(sll_logical_mesh_2d), pointer :: logical_mesh2d

    Neta1 = sim%Neta1
    Neta2 = sim%Neta2
    SLL_ALLOCATE(A11(Neta1,Neta2),ierr)
    SLL_ALLOCATE(A12(Neta1,Neta2),ierr)
    SLL_ALLOCATE(A21(Neta1,Neta2),ierr)
    SLL_ALLOCATE(A22(Neta1,Neta2),ierr)
    SLL_ALLOCATE(B1(Neta1,Neta2),ierr)
    SLL_ALLOCATE(B2(Neta1,Neta2),ierr)
    SLL_ALLOCATE(C(Neta1,Neta2),ierr)
    
    !---> Initialization of the matrices A11, A12, A21, A22, B1, B2 and C
    call initialize_matrix_A_QN_DK (sim,A11,A12,A21,A22) 
    call initialize_vector_B_QN_DK (sim,B1,B2) 
    call initialize_scalar_C_QN_DK ( sim, C )
       
    !---> Initialization of the 2D fields associated to
    !--->  A11, A12, A21, A22, B1, B2 and C
    sim%QN_A11 => new_scalar_field_2d_discrete_alt( &
      "QN_A11", &
      sim%interp2d_QN_A11, &     
      sim%transf_xy, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2)
    
    call sim%QN_A11%set_field_data( A11 )    
    call sim%QN_A11%update_interpolation_coefficients( )

    sim%QN_A12 => new_scalar_field_2d_discrete_alt( &
      "QN_A12", &
      sim%interp2d_QN_A12, &     
      sim%transf_xy, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2)

    call sim%QN_A12%set_field_data( A12 )
    call sim%QN_A12%update_interpolation_coefficients( )

    sim%QN_A21 => new_scalar_field_2d_discrete_alt( &
      "QN_A21", &
      sim%interp2d_QN_A21, &     
      sim%transf_xy, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2)

    call sim%QN_A21%set_field_data( A21 )
    call sim%QN_A21%update_interpolation_coefficients( )

    sim%QN_A22 => new_scalar_field_2d_discrete_alt( &
      "QN_A22", &
      sim%interp2d_QN_A22, &     
      sim%transf_xy, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2)

    call sim%QN_A22%set_field_data( A22 )
    call sim%QN_A22%update_interpolation_coefficients( )

    sim%QN_B1 => new_scalar_field_2d_discrete_alt( &
         "QN_B1", &
         sim%interp2d_QN_B1, &     
         sim%transf_xy, &
         sim%bc_left_eta1, &
         sim%bc_right_eta1, &
         sim%bc_left_eta2, &
         sim%bc_right_eta2)
    
    call sim%QN_B1%set_field_data( B1 )
    call sim%QN_B1%update_interpolation_coefficients( )
    
    sim%QN_B2 => new_scalar_field_2d_discrete_alt( &
         "QN_B2", &
         sim%interp2d_QN_B1, &     
         sim%transf_xy, &
         sim%bc_left_eta1, &
         sim%bc_right_eta1, &
         sim%bc_left_eta2, &
         sim%bc_right_eta2)
    
    call sim%QN_B2%set_field_data( B1 )
    call sim%QN_B2%update_interpolation_coefficients( )
    
    sim%QN_C => new_scalar_field_2d_discrete_alt( &
      "QN_C", &
      sim%interp2d_QN_C, &     
      sim%transf_xy, &
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2)

    call sim%QN_C%set_field_data( C )
    call sim%QN_C%update_interpolation_coefficients( )

    !---> Initialization of the QNS type
    logical_mesh2d => sim%transf_xy%mesh

    sim%QNS => new_general_elliptic_solver( &
      sim%spline_degree_eta1, & 
      sim%spline_degree_eta2, & 
      logical_mesh2d%num_cells1, &
      logical_mesh2d%num_cells2, &
      ES_GAUSS_LEGENDRE, &  ! put in arguments
      ES_GAUSS_LEGENDRE, &  ! put in arguments
      sim%bc_left_eta1, &
      sim%bc_right_eta1, &
      sim%bc_left_eta2, &
      sim%bc_right_eta2, &
      logical_mesh2d%eta1_min, &  
      logical_mesh2d%eta1_max, & 
      logical_mesh2d%eta2_min, & 
      logical_mesh2d%eta2_max ) 

    SLL_DEALLOCATE(A11,ierr)
    SLL_DEALLOCATE(A12,ierr)
    SLL_DEALLOCATE(A21,ierr)
    SLL_DEALLOCATE(A22,ierr)
    SLL_DEALLOCATE(B1,ierr)
    SLL_DEALLOCATE(B2,ierr)
    SLL_DEALLOCATE(C,ierr)

  end subroutine initialize_QN_DK


  !----------------------------------------------------
  ! Initialization of the QN coefficients of matrix A
  !----------------------------------------------------
  ! In the case of drift-Kinetic and with F the change of variables such that 
  !
  !  F( eta1,eta2) =( F_1(eta1, eta2),F_2(eta1, eta2) ) = ( x, y )
  !
  !    ( -1   0 )
  ! A =(        )*  n_0 ( F_1(eta1, eta2),F_2(eta1, eta2))/(B(F_1(eta1, eta2),F_2(eta1, eta2)))
  !    ( 0   -1 ) 
  !----------------------------------------------------
  subroutine initialize_matrix_A_QN_DK (&
       sim,&
       values_A11,&
       values_A12,&
       values_A21,&
       values_A22 )
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    sll_real64, dimension(:,:), intent(inout):: values_A11
    sll_real64, dimension(:,:), intent(inout):: values_A12
    sll_real64, dimension(:,:), intent(inout):: values_A21
    sll_real64, dimension(:,:), intent(inout):: values_A22

    sll_int32  :: ieta1, ieta2
    sll_int32  :: Neta1, Neta2
    sll_real64 :: n0_xy_tmp, B_xy_tmp
    
    Neta1 = sim%Neta1
    Neta2 = sim%Neta2

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
        n0_xy_tmp = sim%n0_xy(ieta1,ieta2)
        B_xy_tmp  = sim%B_xy(ieta1,ieta2) 
        values_A11(ieta1,ieta2) = &
             - n0_xy_tmp/B_xy_tmp
        values_A12(ieta1,ieta2) = 0.0_f64 
        values_A21(ieta1,ieta2) = 0.0_f64
        values_A22(ieta1,ieta2) = &
             - n0_xy_tmp/B_xy_tmp
      end do
    end do
  end subroutine initialize_matrix_A_QN_DK


  !----------------------------------------------------
  ! Initialization of the QN coefficients of vector B
  !----------------------------------------------------
  ! In the case of drift-Kinetic and with F the change of variables such that 
  !
  !       (  F_1(eta1, eta2) / norm**2 )
  ! B =   (                            )*  n_0 ( F_1(eta1, eta2),F_2(eta1, eta2))/(B(F_1(eta1, eta2),F_2(eta1, eta2))* omega_0 )
  !       (  F_2(eta1, eta2) / norm**2 ) 
  ! with 
  ! norm**2 = (F_1(eta1, eta2)**2 + F_2(eta1, eta2)**2)
  !----------------------------------------------------
   subroutine initialize_vector_B_QN_DK (&
        sim,&
        values_B1,&
        values_B2 )

    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    sll_real64, dimension(:,:), pointer :: values_B1
    sll_real64, dimension(:,:), pointer :: values_B2

    sll_int32  :: ieta1, ieta2
    sll_int32  :: Neta1, Neta2
    sll_real64 :: x_tmp, y_tmp
    sll_real64 :: n0_xy_tmp, B_xy_tmp
    sll_real64 :: norm2_xy_tmp, val_tmp
    
    Neta1 = sim%Neta1
    Neta2 = sim%Neta2
    
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
        x_tmp        = sim%xgrid_2d(ieta1,ieta2)
        y_tmp        = sim%ygrid_2d(ieta1,ieta2)
        n0_xy_tmp    = sim%n0_xy(ieta1,ieta2)
        B_xy_tmp     = sim%B_xy(ieta1,ieta2) 
        norm2_xy_tmp = sim%norm_square_xy(ieta1,ieta2) 
        val_tmp      = n0_xy_tmp/(B_xy_tmp*norm2_xy_tmp)
        values_B1(ieta1,ieta2) = 0.0!-val_tmp*x_tmp
        values_B2(ieta1,ieta2) = 0.0!-val_tmp*y_tmp
      end do
    end do

  end subroutine initialize_vector_B_QN_DK


  !----------------------------------------------------
  ! Initialization of the QN coefficients of scalar C
  !----------------------------------------------------
  ! In the case of drift-Kinetic and with F the change of variables such that 
  !
  ! C = e n_0 ( F_1(eta1, eta2),F_2(eta1, eta2))/( T_e (F_1(eta1, eta2),F_2(eta1, eta2)) )
  !---------------------------------------------------- 
  subroutine initialize_scalar_C_QN_DK (&
       sim,&
       values_C )

    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    sll_real64, dimension(:,:), pointer :: values_C

    sll_int32  :: ieta1, ieta2
    sll_int32  :: Neta1, Neta2
    
    Neta1 = sim%Neta1
    Neta2 = sim%Neta2

    if ( (size(values_C,1) .ne. Neta1) .OR. &
         (size(values_C,2) .ne. Neta2) ) then
       print*, ' Problem with the dimension of C'
    end if

    do ieta2 = 1,Neta2
      do ieta1 = 1,Neta1
         values_C(ieta1,ieta2) = &
              sim%n0_xy(ieta1,ieta2) / &
              sim%Te_xy(ieta1,ieta2)
      end do
    end do

  end subroutine initialize_scalar_C_QN_DK


  !-----------------------------------------------------------
  ! Computation of the charge density, i.e   
  !   rho(x,y,z) = \int delta f dvpar
  !   delta f = f(x,y,z,vpar) - feq(x,y,vpar) 
  !
  ! Rk: This corresponds to the RHS of the quasi-neutrality
  !     equation
  !
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !       feq3d_seqx1x2x4(x1=*,x2=*,x4=*)
  !  Out: rho3d_seqx3(x1=distrib,x2=distrib,x3=*)
  !-----------------------------------------------------------
  subroutine compute_charge_density(sim)

    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32  :: Neta1_loc,Neta2_loc,Neta3, Nvpar
    sll_int32  :: iloc1, iloc2
    sll_int32  :: i1, i2, i3, i4
    sll_real64 :: delta_f
    sll_real64 :: delta_vpar, intf_dvpar
    sll_int32, dimension(1:4) :: glob_ind4d

    Neta1_loc  = size(sim%f4d_seqx3x4,1)
    Neta2_loc  = size(sim%f4d_seqx3x4,2)
    Neta3      = size(sim%f4d_seqx3x4,3)
    Nvpar      = size(sim%f4d_seqx3x4,4)
    delta_vpar = sim%logical_mesh4d%delta_eta4 

    !-> Computation of the charge density locally in (x1,x2) directions
    do i3 = 1,Neta3
      do iloc2 = 1,Neta2_loc
        do iloc1 = 1,Neta1_loc

          intf_dvpar = 0._f64

          do i4 = 1,Nvpar-1
            glob_ind4d(:) = local_to_global_4D(sim%layout4d_seqx3x4, &
                  (/iloc1,iloc2,i3,i4/))
            i1 = glob_ind4d(1)
            i2 = glob_ind4d(2)
            delta_f = sim%f4d_seqx3x4(iloc1,iloc2,i3,i4) - &
                sim%feq_xyvpar(i1,i2,i4)
            intf_dvpar = intf_dvpar + delta_f*delta_vpar           
          end do

          sim%rho3d_seqx3(iloc1,iloc2,i3) = intf_dvpar

        end do
      end do
    end do
    !--> compute rho3d_seqx1x2
    call apply_remap_3D(sim%seqx3_to_seqx1x2,sim%rho3d_seqx3,sim%rho3d_seqx1x2)

  end subroutine compute_charge_density


  !----------------------------------------------------
  ! Compute Phi(eta1,eta1,phi) by solving
  !  the quasi-neutrality equation by using
  !  a general elliptic equation solver
  !----------------------------------------------------
  subroutine solve_QN( sim )
    use sll_common_coordinate_transformations, only : &
      polar_eta2
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32 :: ieta1, ieta2, iloc3
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3

    call compute_local_sizes_3d( sim%layout3d_seqx1x2, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)

    !---> Compute phi3d_seqx1x2
    do iloc3 = 1,loc3d_sz_x3
      call sim%rho2d%set_field_data( sim%rho3d_seqx1x2(:,:,iloc3) )
      call sim%rho2d%update_interpolation_coefficients( )
      call solve_general_coordinates_elliptic_eq( &
        sim%QNS, &
        sim%rho2d, &
        sim%phi2d)
      do ieta2 = 1,sim%Neta2
        do ieta1 = 1,sim%Neta1 
           sim%phi3d_seqx1x2(ieta1,ieta2,iloc3) = sim%phi2d%value_at_indices(ieta1,ieta2)
        end do
      end do
     ! print*,'first', sim%phi3d_seqx1x2(1,:,iloc3)
     ! print*,'second', sim%phi3d_seqx1x2(sim%Neta1,:,iloc3)
    end do
    if (sim%my_rank.eq.0) &
      call sim%phi2d%write_to_file(iloc3)
    
    !---> Fill phi3d_seqx3
    call apply_remap_3D(sim%seqx1x2_to_seqx3, sim%phi3d_seqx1x2, sim%phi3d_seqx3)    
  end subroutine solve_QN


  !--------------------------------------------------------------
  ! Compute the electric field
  !   E = - grad Phi 
  ! where Phi is the electrostatic potential solution of
  ! quasi-neutrality equation
  !--------------------------------------------------------------
  subroutine compute_Efield( sim )
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32  :: ierr
    sll_int32  :: ieta1, ieta2, ieta3
    sll_int32  :: iloc1, iloc2, iloc3
    sll_int32  :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3
    sll_real64 :: eta1, eta2, eta3
    sll_real64, dimension(:), pointer :: phi1d_seqx3_tmp
    sll_real64, dimension(2,2)  :: matinv_jac

    !--> Compute E3d_eta1_seqx1x2 = -dPhi3d_seqx1x2/deta1 and 
    !-->  E3d_eta2_seqx1x2 = -dPhi3d_seqx1x2/deta2
    call compute_local_sizes_3d( sim%layout3d_seqx1x2, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
      
    do iloc3 = 1,loc3d_sz_x3
      call sim%phi2d%set_field_data( sim%phi3d_seqx1x2(:,:,iloc3) )
      call sim%phi2d%update_interpolation_coefficients( )
     ! print*,'first', sim%phi3d_seqx1x2(1,:,iloc3)
     ! print*,'second', sim%phi3d_seqx1x2(sim%Neta1,:,iloc3)
      do ieta2 = 1,sim%Neta2
        eta2 = sim%eta2_grid(ieta2)
        do ieta1 = 1,sim%Neta1
          eta1 = sim%eta1_grid(ieta1)
          matinv_jac(:,:) = sim%transf_xy%inverse_jacobian_matrix(eta1,eta2) 

          sim%E3d_eta1_seqx1x2(ieta1,ieta2,iloc3) = &
            - sim%phi2d%first_deriv_eta1_value_at_point(eta1,eta2)

          sim%E3d_eta2_seqx1x2(ieta1,ieta2,iloc3) = &
            - sim%phi2d%first_deriv_eta2_value_at_point(eta1,eta2)

          sim%E3d_x1_seqx1x2(ieta1,ieta2,iloc3) = &
               matinv_jac(1,1)*sim%E3d_eta1_seqx1x2(ieta1,ieta2,iloc3) &
               + matinv_jac(2,1)*sim%E3d_eta2_seqx1x2(ieta1,ieta2,iloc3)

          sim%E3d_x2_seqx1x2(ieta1,ieta2,iloc3) = &
               matinv_jac(1,2)*sim%E3d_eta1_seqx1x2(ieta1,ieta2,iloc3) &
               + matinv_jac(2,2)*sim%E3d_eta2_seqx1x2(ieta1,ieta2,iloc3)
        end do
      end do
    end do
     
    !--> Compute E3d_eta3_seqx3 = -dPhi3d_seqx3/deta3
    SLL_ALLOCATE(phi1d_seqx3_tmp(sim%Neta3),ierr)
    call compute_local_sizes_3d( sim%layout3d_seqx3, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    do iloc2 = 1,loc3d_sz_x2
      do iloc1 = 1,loc3d_sz_x1
        do ieta3 = 1,sim%Neta3
          phi1d_seqx3_tmp(ieta3) = sim%phi3d_seqx3(iloc1,iloc2,ieta3)
        end do
       
        call sim%phi1d%set_field_data(phi1d_seqx3_tmp)
        call sim%phi1d%update_interpolation_coefficients( ) 
        do ieta3 = 1,sim%Neta3 
          eta3 = sim%eta3_grid(ieta3)
          sim%E3d_eta3_seqx3(iloc1,iloc2,ieta3) = &
            -sim%phi1d%derivative_value_at_point(eta3)
        end do
      end do
    end do
    SLL_DEALLOCATE(phi1d_seqx3_tmp,ierr)
  end subroutine compute_Efield



  !************************************************************
  !  VLASOV SOLVING
  !************************************************************

  !----------------------------------------------------
  ! First step of the drift-kinetic 4D simulation :
  !  1) Initialization of the distribution function 
  !  2) Computation of the charge density 
  !      (r.h.s of the quasi-neutrality solver)
  !  3) Solving of the quasi-neutrality equation
  !----------------------------------------------------
  subroutine first_step_4d_DK_hybrid( sim )
    use sll_timer
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    
    sll_int32 :: diag_num

    !--> For timers
    type(sll_time_mark) :: t0, t1 
    double precision    :: elaps_time

    !*** Initialization of the distribution function ***
    !***  i.e f4d(t=t0)                              ***
    call sll_set_time_mark(t0)
    print*, 'initialize fdistribution'
    call initialize_fdistribu4d_DK(sim)
    call sll_set_time_mark(t1)
    elaps_time = sll_time_elapsed_between(t0,t1)
    if (sim%my_rank.eq.0) &
      print*, ' Time for initialize_fdistribu4d_DK = ', elaps_time
        
    !*** Saving of the equilibrium state in HDF5 file 'init_state.h5' ***
    call writeHDF5_equilibrium_state(sim)

    !*** Computation of the rhs of QN ***
    call sll_set_time_mark(t0)
    call compute_charge_density(sim)
    !--> compute rho3d_seqx1x2
    call apply_remap_3D(sim%seqx3_to_seqx1x2, sim%rho3d_seqx3, sim%rho3d_seqx1x2 )
    call sll_set_time_mark(t1)
    elaps_time = sll_time_elapsed_between(t0,t1)
    if (sim%my_rank.eq.0) &
      print*, ' Time for compute_charge_density = ', elaps_time

    !*** Matrix factorization for QN solver ***
    call sll_set_time_mark(t0)
    call factorize_mat_es( &
         sim%QNS, & 
         sim%QN_A11, &
         sim%QN_A12, &
         sim%QN_A21, &
         sim%QN_A22, &
         sim%QN_B1,  &
         sim%QN_B2,  &
         sim%QN_C)
    call sll_set_time_mark(t1)
    elaps_time = sll_time_elapsed_between(t0,t1)
    if (sim%my_rank.eq.0) &
      print*, ' Time for factorize_mat_es = ', elaps_time

    !*** Compute Phi(t=t0) by solving the QN equation ***
    call sll_set_time_mark(t0)
    call solve_QN(sim)
    call sll_set_time_mark(t1)
    elaps_time = sll_time_elapsed_between(t0,t1)
    if (sim%my_rank.eq.0) &
      print*, ' Time for solve_QN = ', elaps_time

    !*** Compute E = -grad Phi ***
    call sll_set_time_mark(t0)
    call compute_Efield( sim )
    call sll_set_time_mark(t1)
    elaps_time = sll_time_elapsed_between(t0,t1)
    if (sim%my_rank.eq.0) &
      print*, ' Time for compute_Efield = ', elaps_time

    !*** Diagnostics at time t=0 ***
    sim%count_save_diag = 1
    diag_num = sim%count_save_diag
    !--> Compute energy kinetic, potential and total
    call compute_energy(sim,diag_num)
    
    !--> Compute L1 norm, L2 norm, L infini norm
    call compute_norm_L1_L2_Linf(sim,diag_num)

    !--> Save initial step in HDF5 file
    call writeHDF5_cross_section_diag(sim,diag_num)    
  end subroutine first_step_4d_DK_hybrid
  

  !----------------------------------------------------
  ! 1D advection in vpar direction
  !----------------------------------------------------
  subroutine advec1D_vpar( sim , deltat_advec )
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    sll_real64                        , intent(in)    :: deltat_advec

    sll_int32 :: ierr
    sll_int32 :: iloc1, iloc2
    sll_int32 :: ieta3, ivpar
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32 :: loc4d_sz_x3, loc4d_sz_x4
    sll_real64 :: eta3, vpar, E_z, alpha4
    sll_real64, dimension(:), pointer :: f1d_vpar_tmp

    SLL_ASSERT(size(sim%f4d_seqx3x4,1).eq.size(sim%phi3d_seqx3,1))
    SLL_ASSERT(size(sim%f4d_seqx3x4,2).eq.size(sim%phi3d_seqx3,2))
    SLL_ASSERT(size(sim%f4d_seqx3x4,3).eq.size(sim%phi3d_seqx3,3))

    call compute_local_sizes_4d( sim%layout4d_seqx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    

    !---> dvpar/dt = E_z with E_z = -dPhi/dz
    SLL_ALLOCATE(f1d_vpar_tmp(sim%Nvpar),ierr)
    do iloc2 = 1,loc4d_sz_x2
      do iloc1 = 1,loc4d_sz_x1
        do ieta3 = 1,sim%Neta3
          eta3   = sim%eta3_grid(ieta3)
          E_z    = sim%E3d_eta3_seqx3(iloc1,iloc2,ieta3)
          do ivpar = 1,sim%Nvpar
            f1d_vpar_tmp(ivpar) = sim%f4d_seqx3x4(iloc1,iloc2,ieta3,ivpar)

          end do
          call sim%interp1d_f_vpar%compute_interpolants(f1d_vpar_tmp)   
          do ivpar = 1,sim%Nvpar
            ! si change of coordinates in 3D
            ! val_jac  = sim%transf_xy%jacobian(eta1,eta2,eta3) 
            alpha4 = deltat_advec*E_z 
            vpar   = sim%vpar_grid(ivpar) - alpha4
            vpar   = max(min(vpar,sim%vpar_max),sim%vpar_min)

            sim%f4d_seqx3x4(iloc1,iloc2,ieta3,ivpar) = &
              sim%interp1d_f_vpar%interpolate_value(vpar)             
          end do
        end do
      end do
    end do
    SLL_DEALLOCATE(f1d_vpar_tmp,ierr)
  end subroutine advec1D_vpar

  
  !----------------------------------------------------
  ! 1D advection in eta3 direction
  !----------------------------------------------------
  subroutine advec1D_eta3( sim, deltat_advec )

    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    sll_real64                        , intent(in)    :: deltat_advec

    sll_int32 :: ierr
    sll_int32 :: iloc1, iloc2
    sll_int32 :: ieta3, ivpar
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32 :: loc4d_sz_x3, loc4d_sz_x4
    sll_real64:: eta3, alpha3, vpar
    sll_real64, dimension(:), pointer :: f1d_eta3_tmp

    call compute_local_sizes_4d( sim%layout4d_seqx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    

    !---> deta3/dt = vpar
    SLL_ALLOCATE(f1d_eta3_tmp(sim%Neta3),ierr)
    do iloc2 = 1,loc4d_sz_x2
      do iloc1 = 1,loc4d_sz_x1
        do ivpar = 1,sim%Nvpar
          vpar   = sim%vpar_grid(ivpar)
          do ieta3 = 1,sim%Neta3
            f1d_eta3_tmp(ieta3) = sim%f4d_seqx3x4(iloc1,iloc2,ieta3,ivpar)               
          end do
          call sim%interp1d_f_eta3%compute_interpolants(f1d_eta3_tmp)  


          do ieta3 = 1,sim%Neta3
            alpha3 = deltat_advec*vpar
            eta3    = sim%eta3_grid(ieta3) - alpha3

            sim%f4d_seqx3x4(iloc1,iloc2,ieta3,ivpar) = &
              sim%interp1d_f_eta3%interpolate_value(eta3)
          end do
        end do
      end do
    end do
    SLL_DEALLOCATE(f1d_eta3_tmp,ierr)

  end subroutine advec1D_eta3


  !----------------------------------------------------
  ! 2D advection in eta1eta2 direction
  !----------------------------------------------------
  subroutine advec2D_eta1eta2( sim, deltat_advec )

    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    sll_real64                        , intent(in)    :: deltat_advec

    sll_int32 :: ierr
    sll_int32 :: iloc3, iloc4
    sll_int32 :: ieta1, ieta2
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32 :: loc4d_sz_x3, loc4d_sz_x4
    sll_real64 :: eta1,eta2, alpha1,alpha2, E_eta1, E_eta2
    sll_real64 :: val_jac,val_B
    sll_real64, dimension(:,:), pointer :: f2d_eta1eta2_tmp

    call compute_local_sizes_4d( sim%layout4d_seqx1x2, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    

    !---> ( deta1/dt )                              ( 0    -1 ) (  d phi/deta1 ) 
    !     (          )=  1 / (B* jac(eta_1,eta2) )  (         ) (              )
    !---> ( deta2/dt )                              ( 1     0 ) (  d phi/deta2 )
    SLL_ALLOCATE(f2d_eta1eta2_tmp(sim%Neta1,sim%Neta2),ierr)

    do iloc4 = 1,loc4d_sz_x4
      do iloc3 = 1,loc4d_sz_x3
        do ieta2 = 1,sim%Neta2
          do ieta1 = 1,sim%Neta1
            f2d_eta1eta2_tmp(ieta1,ieta2)=&
              sim%f4d_seqx1x2(ieta1,ieta2,iloc3,iloc4)            
          end do
        end do

        call sim%interp2d_f_eta1eta2%compute_interpolants(f2d_eta1eta2_tmp)   

        do ieta2 = 1,sim%Neta2
          eta2 = sim%eta2_grid(ieta2)
          do ieta1 = 1,sim%Neta1
            eta1    = sim%eta1_grid(ieta1)
            E_eta1  = sim%E3d_eta1_seqx1x2(ieta1,ieta2,iloc3)
            E_eta2  = sim%E3d_eta2_seqx1x2(ieta1,ieta2,iloc3)
            val_jac = sim%transf_xy%jacobian(eta1,eta2)
            !print*, val_jac
            val_B   = sim%B_xy(ieta1,ieta2)
            alpha1  =  deltat_advec*E_eta2/ (val_jac*val_B)
            alpha2  =  -deltat_advec*E_eta1/ (val_jac*val_B)
            eta1    = sim%eta1_grid(ieta1) - alpha1
            eta2    = sim%eta2_grid(ieta2) - alpha2

            sim%f4d_seqx1x2(ieta1,ieta2,iloc3,iloc4) = &
              sim%interp2d_f_eta1eta2%interpolate_value(eta1,eta2)

          end do
        end do
      end do
    end do
    SLL_DEALLOCATE(f2d_eta1eta2_tmp,ierr)

  end subroutine advec2D_eta1eta2



  !************************************************************
  !  DIAGNOSTICS
  !************************************************************

  !------------------------------------------------------------
  ! Save the initial density and temperature profiles
  !  and the equilibrium state in the file 'init_state.h5'
  !------------------------------------------------------------
  subroutine writeHDF5_equilibrium_state( sim )

    use sll_hdf5_io, only: sll_hdf5_file_create, &
        sll_hdf5_write_array_1d, sll_hdf5_file_close

    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    !--> For initial equilibrium state HDF5 saving
    integer                      :: file_err
    sll_int32                    :: file_id
    character(len=13), parameter :: filename_prof = "init_state.h5"

    if (sim%my_rank.eq.0) then
      call sll_hdf5_file_create(filename_prof,file_id,file_err)
      !--> Saving of the 1D radial profiles of 
      !-->  temperatures and density
      call sll_hdf5_write_array_1d(file_id,sim%r_grid,'r_grid',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%n0_r,'n0_r',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Ti_r,'Ti_r',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Te_r,'Te_r',file_err)

      !--> Saving of the 2D (x,y) profiles of 
      !-->  temperatures and density
      call sll_hdf5_write_array_2d(file_id,sim%xgrid_2d,'xgrid_2d',file_err)
      call sll_hdf5_write_array_2d(file_id,sim%ygrid_2d,'ygrid_2d',file_err)
      call sll_hdf5_write_array_2d(file_id,sim%n0_xy,'n0_xy',file_err)
      call sll_hdf5_write_array_2d(file_id,sim%Ti_xy,'Ti_xy',file_err)
      call sll_hdf5_write_array_2d(file_id,sim%Te_xy,'Te_xy',file_err)
      call sll_hdf5_write_array_2d(file_id,sim%B_xy,'B_xy',file_err)
      !---> Saving of the 3D equilibrium function feq(x,y,vpar)
      call sll_hdf5_write_array_3d(file_id,sim%feq_xyvpar,'feq_xyvpar',file_err)
      call sll_hdf5_file_close(file_id,file_err)
    end if

  end subroutine writeHDF5_equilibrium_state
 

  !--------------------------------------------------------------
  ! Save 2D cross-sections for:
  !  - Electrostatic potential
  !  - Electric field
  !  - Distribution function
  !
  ! Rk: One file 'DK4d_diag_d<num_diag>.h5' per diagnostic
  !--------------------------------------------------------------
  subroutine writeHDF5_cross_section_diag( sim, diag_num )

    use sll_hdf5_io, only: sll_hdf5_file_create, &
      sll_hdf5_write_array_1d, sll_hdf5_file_close
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    sll_int32                         , intent(in)    :: diag_num

    sll_int32  :: ix1_diag, ix2_diag
    sll_int32  :: ix3_diag, ivpar_diag

    sll_real64, dimension(1) :: iter_time_tmp

    !--> For initial profile HDF5 saving
    integer             :: file_err
    sll_int32           :: file_id
    character(len=80)   :: filename_HDF5
    character(20), save :: numfmt = "'_d',i5.5"
    
    ix1_diag   = int(sim%Neta1/2)
    ix2_diag   = int(sim%Neta2/3)
    ix3_diag   = int(sim%Neta3/4)
    ivpar_diag = int(sim%Nvpar/3)

    write(filename_HDF5,'(A,'//numfmt//',A)') &
      "DK4d_diag", diag_num, ".h5"

    if (sim%my_rank.eq.0) then
      print*,'--> Save HDF5 file: ',filename_HDF5
      call sll_hdf5_file_create(filename_HDF5,file_id,file_err)
      !---> Time
      iter_time_tmp = sim%iter_time
      call sll_hdf5_write_array_1d(file_id, &
        iter_time_tmp,'time_diag',file_err)
      !---> Mesh 1D
      call sll_hdf5_write_array_1d(file_id, &
        sim%eta1_grid,'eta1_grid',file_err)
      call sll_hdf5_write_array_1d(file_id, &
        sim%eta2_grid,'eta2_grid',file_err)
      call sll_hdf5_write_array_1d(file_id, &
        sim%eta3_grid,'eta3_grid',file_err)
      call sll_hdf5_write_array_1d(file_id, &
        sim%vpar_grid,'vpar_grid',file_err)
      !---> Mesh 2D
      call sll_hdf5_write_array_2d(file_id, &
        sim%xgrid_2d,'xgrid_2d',file_err)
      call sll_hdf5_write_array_2d(file_id, &
        sim%ygrid_2d,'ygrid_2d',file_err)
      !---> RHS of QN equation
      call sll_hdf5_write_array_2d(file_id, &
        sim%rho3d_seqx1x2(:,:,ix3_diag),'rho2d_xy',file_err)
      !---> Electrostatic potential
      call sll_hdf5_write_array_2d(file_id, &
        sim%Phi3d_seqx1x2(:,:,ix3_diag),'phi2d_xy',file_err)
      !---> Electric field
      call sll_hdf5_write_array_2d(file_id, &
        sim%E3d_eta1_seqx1x2(:,:,ix3_diag),'E2d_eta1_xy',file_err)
      call sll_hdf5_write_array_2d(file_id, &
        sim%E3d_eta2_seqx1x2(:,:,ix3_diag),'E2d_eta2_xy',file_err)
      call sll_hdf5_write_array_2d(file_id, &
        sim%E3d_x1_seqx1x2(:,:,ix3_diag),'E2d_x1_xy',file_err)
      call sll_hdf5_write_array_2d(file_id, &
        sim%E3d_x2_seqx1x2(:,:,ix3_diag),'E2d_x2_xy',file_err)
      !---> Distribution function
      call sll_hdf5_write_array_2d(file_id, &
        sim%f4d_seqx1x2(:,:,ix3_diag,ivpar_diag),'f2d_xy',file_err)
      call sll_hdf5_write_array_2d(file_id, &
        sim%f4d_seqx3x4(ix1_diag,ix2_diag,:,:),'f2d_zvpar',file_err)
    end if

  end subroutine writeHDF5_cross_section_diag
  

  !----------------------------------------------------
  ! 
  !----------------------------------------------------
  subroutine writeHDF5_conservation_laws( sim )

    use sll_hdf5_io, only: sll_hdf5_file_create, &
      sll_hdf5_write_array_1d, sll_hdf5_file_close

    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    !--> diagnostics norm
    sll_real64, dimension(:), pointer :: diag_mass_tmp
    sll_real64, dimension(:), pointer :: diag_norm_L1_tmp
    sll_real64, dimension(:), pointer :: diag_norm_L2_tmp
    sll_real64, dimension(:), pointer :: diag_norm_Linf_tmp
    sll_real64, dimension(:), pointer :: diag_entropy_kin_tmp
    
    !--> diagnostics energy
    sll_real64, dimension(:), pointer :: diag_nrj_kin_tmp
    sll_real64, dimension(:), pointer :: diag_nrj_pot_tmp
    sll_real64, dimension(:), pointer :: diag_nrj_tot_tmp
    sll_real64, dimension(:), pointer :: diag_heat_flux_tmp
    sll_real64, dimension(:), pointer :: diag_relative_error_nrj_tot_tmp

    sll_int32  :: ierr
    sll_int32  :: idiag, nb_diag
    sll_real64 :: dum
    
    !--> For conservation law HDF5 saving
    integer             :: file_err
    sll_int32           :: file_id
    character(len=20), parameter :: filename_CL = "conservation_laws.h5"

    nb_diag = size(sim%diag_mass)

    SLL_ALLOCATE(diag_mass_tmp(nb_diag),ierr)
    SLL_ALLOCATE(diag_norm_L1_tmp(nb_diag),ierr)
    SLL_ALLOCATE(diag_norm_L2_tmp(nb_diag),ierr)
    SLL_ALLOCATE(diag_norm_Linf_tmp(nb_diag),ierr)
    SLL_ALLOCATE(diag_entropy_kin_tmp(nb_diag),ierr)
    SLL_ALLOCATE(diag_nrj_kin_tmp(nb_diag),ierr)
    SLL_ALLOCATE(diag_nrj_pot_tmp(nb_diag),ierr)
    SLL_ALLOCATE(diag_nrj_tot_tmp(nb_diag),ierr)
    SLL_ALLOCATE(diag_heat_flux_tmp(nb_diag),ierr)
    SLL_ALLOCATE(diag_relative_error_nrj_tot_tmp(nb_diag),ierr)

    diag_mass_tmp        = 0.0_f64
    diag_norm_L1_tmp     = 0.0_f64
    diag_norm_L2_tmp     = 0.0_f64
    diag_norm_Linf_tmp   = 0.0_f64
    diag_entropy_kin_tmp = 0.0_f64
    diag_nrj_kin_tmp     = 0.0_f64
    diag_nrj_pot_tmp     = 0.0_f64
    diag_nrj_tot_tmp     = 0.0_f64
    diag_heat_flux_tmp   = 0.0_f64
    diag_relative_error_nrj_tot_tmp = 0.0_f64

    call sll_collective_reduce_real64( &
        sll_world_collective, &
        sim%diag_mass, &
        nb_diag, &
        MPI_SUM, &
        0, &
        diag_mass_tmp )

    call sll_collective_reduce_real64( &
        sll_world_collective, &
        sim%diag_norm_L1, &
        nb_diag, &
        MPI_SUM, &
        0, &
        diag_norm_L1_tmp )

    call sll_collective_reduce_real64( &
        sll_world_collective, &
        sim%diag_norm_L2, &
        nb_diag, &
        MPI_SUM, &
        0, &
        diag_norm_L2_tmp )

    call sll_collective_reduce_real64( &
        sll_world_collective, &
        sim%diag_norm_Linf, &
        nb_diag, &
        MPI_SUM, &
        0, &
        diag_norm_Linf_tmp )

    call sll_collective_reduce_real64( &
        sll_world_collective, &
        sim%diag_entropy_kin, &
        nb_diag, &
        MPI_SUM, &
        0, &
        diag_entropy_kin_tmp )

    call sll_collective_reduce_real64( &
        sll_world_collective, &
        sim%diag_nrj_kin, &
        nb_diag, &
        MPI_SUM, &
        0, &
        diag_nrj_kin_tmp )

    call sll_collective_reduce_real64( &
        sll_world_collective, &
        sim%diag_heat_flux, &
        nb_diag, &
        MPI_SUM, &
        0, &
        diag_heat_flux_tmp )

    call sll_collective_reduce_real64( &
        sll_world_collective, &
        sim%diag_nrj_pot, &
        nb_diag, &
        MPI_SUM, &
        0, &
        diag_nrj_pot_tmp )

    call sll_collective_reduce_real64( &
        sll_world_collective, &
        sim%diag_nrj_tot, &
        nb_diag, &
        MPI_SUM, &
        0, &
        diag_nrj_tot_tmp )

    do idiag = 1,nb_diag
      dum  = sqrt(.5*((diag_nrj_kin_tmp(idiag)-diag_nrj_kin_tmp(1))**2 + &
          (diag_nrj_pot_tmp(idiag)-diag_nrj_pot_tmp(1))**2)) 

      if ( dum /= 0.0_f64) &
          diag_relative_error_nrj_tot_tmp(idiag) = &
          (diag_nrj_tot_tmp(idiag)-diag_nrj_tot_tmp(1))/dum
    end do

    if (sim%my_rank.eq.0) then
      print*,'--> Save HDF5 file: ',filename_CL
      call sll_hdf5_file_create(filename_CL,file_id,file_err)
      call sll_hdf5_write_array_1d(file_id,&
          sim%time_evol(:),'time_evol',file_err)
      call sll_hdf5_write_array_1d(file_id,&
          diag_nrj_kin_tmp(:),'nrj_kin',file_err)
      call sll_hdf5_write_array_1d(file_id,&
          diag_nrj_pot_tmp(:),'nrj_pot',file_err)
      call sll_hdf5_write_array_1d(file_id,&
          diag_nrj_tot_tmp(:),'nrj_tot',file_err)
      call sll_hdf5_write_array_1d(file_id,&
          diag_relative_error_nrj_tot_tmp(:),&
          'relative_error_nrj_tot',file_err)
      call sll_hdf5_write_array_1d(file_id,&
          diag_heat_flux_tmp(:),'heat_flux',file_err)
      call sll_hdf5_write_array_1d(file_id,&
          diag_mass_tmp(:),&
          'mass',file_err)
      call sll_hdf5_write_array_1d(file_id,&
          diag_norm_L1_tmp(:),'L1_norm',file_err)
      call sll_hdf5_write_array_1d(file_id,&
          diag_norm_L2_tmp(:),'L2_norm',file_err)
      call sll_hdf5_write_array_1d(file_id,&
          diag_norm_Linf_tmp(:),'Linf_norm',file_err)
      call sll_hdf5_write_array_1d(file_id,&
          diag_entropy_kin_tmp(:),'entropy_kin',file_err)
      call sll_hdf5_file_close(file_id,file_err)
    end if

    SLL_DEALLOCATE(diag_mass_tmp,ierr)
    SLL_DEALLOCATE(diag_norm_L1_tmp,ierr)
    SLL_DEALLOCATE(diag_norm_L2_tmp,ierr)
    SLL_DEALLOCATE(diag_norm_Linf_tmp,ierr)
    SLL_DEALLOCATE(diag_entropy_kin_tmp,ierr)
    SLL_DEALLOCATE(diag_nrj_kin_tmp,ierr)
    SLL_DEALLOCATE(diag_nrj_pot_tmp,ierr)
    SLL_DEALLOCATE(diag_nrj_tot_tmp,ierr)
    SLL_DEALLOCATE(diag_heat_flux_tmp,ierr)
    SLL_DEALLOCATE(diag_relative_error_nrj_tot_tmp,ierr)

  end subroutine writeHDF5_conservation_laws


  !-----------------------------------------------------------
  ! Computation of the kinetic energy, i.e   
  !   nrj_kin = \int delta f * vpar**2 * 0.5 * jac dvpar deta1 deta2 deta3
  !   delta f = f(x,y,z,vpar) - feq(x,y,vpar) 
  !   jac     = jacobian of the change of coordinates
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !       feq3d_seqx1x2x4(x1=*,x2=*,x4=*)
  !  Out: nrj_kin
  !
  ! Computation of the potential energy, i.e   
  !   nrj_pot = \int delta f * phi * 0.5 * jac dvpar deta1 deta2 deta3
  !   delta f = f(eta1,eta2,eta3,vpar) - feq(eta1,eta2,vpar) 
  !   phi     = phi(eta1,eta2,eta3) 
  !   jac     = jacobian of the change of coordinates
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !       feq3d_seqx1x2x4(x1=*,x2=*,x4=*)
  !       phi3d_seqx3(x1=distrib,x2=distrib,x3=*)
  !  Out: nrj_pot
  !
  ! Computation of the energy total = nrj_pot + nrj_kin
  !
  ! Computation of the heat flux, i.e   
  !   nrj_pot = \int delta f * phi * 0.5 * jac dvpar deta1 deta2 deta3
  !   delta f = f(eta1,eta2,eta3,vpar) - feq(eta1,eta2,vpar) 
  !   phi     = phi(eta1,eta2,eta3) 
  !   jac     = jacobian of the change of coordinates
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !       feq3d_seqx1x2x4(x1=*,x2=*,x4=*)
  !       phi3d_seqx3(x1=distrib,x2=distrib,x3=*)
  !  Out: nrj_pot
  !
  !-----------------------------------------------------------  
  subroutine compute_energy(sim,diag_num)

    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim    
    sll_int32                         , intent(in)    :: diag_num

    ! local variables
    sll_int32  :: Neta1_loc,Neta2_loc,Neta3,Nvpar,Neta1,Neta2
    sll_real64 :: delta_eta1,delta_eta2,delta_eta3,delta_vpar
    sll_real64 :: val_jac
    sll_real64 :: vpar
    sll_real64 :: val_elec_pot
    sll_real64 :: delta_f
    sll_int32  :: iloc1, iloc2
    sll_int32  :: i1,i2,i3,i4
    sll_int32, dimension(1:4) :: glob_ind4d
    sll_real64 :: nrj_kin,nrj_pot,nrj_tot,heat_flux

    Neta1_loc  = size(sim%f4d_seqx3x4,1)
    Neta2_loc  = size(sim%f4d_seqx3x4,2)
    Neta1      = sim%nc_x1 + 1
    Neta2      = sim%nc_x2 + 1
    Neta3      = size(sim%f4d_seqx3x4,3)
    Nvpar      = size(sim%f4d_seqx3x4,4)
    delta_eta1 = sim%logical_mesh4d%delta_eta1
    delta_eta2 = sim%logical_mesh4d%delta_eta2
    delta_eta3 = sim%logical_mesh4d%delta_eta3
    delta_vpar = sim%logical_mesh4d%delta_eta4

    sim%time_evol(diag_num) = sim%iter_time

    nrj_kin   = 0.0
    nrj_pot   = 0.0
    nrj_tot   = 0.0
    heat_flux = 0.0

    !-> Computation of the energy kinetic locally in (x1,x2) directions
    do iloc2 = 1,Neta2_loc
      do iloc1 = 1,Neta1_loc
        do i4 = 1,Nvpar-1
          vpar = sim%vpar_grid(i4)
          do i3 = 1,Neta3-1
            val_elec_pot = sim%phi3d_seqx3(iloc1,iloc2,i3)

            glob_ind4d(:) = local_to_global_4D(sim%layout4d_seqx3x4, &
              (/iloc1,iloc2,i3,i4/))
            i1 = glob_ind4d(1)
            i2 = glob_ind4d(2)

            if (i1 .ne. Neta1) then
              if (i2 .ne. Neta2) then 

                val_jac = abs(sim%transf_xy%jacobian_at_node(i1,i2))

                delta_f = sim%f4d_seqx3x4(iloc1,iloc2,i3,i4) - &
                  sim%feq_xyvpar(i1,i2,i4)

                nrj_kin = nrj_kin + &
                  delta_f * vpar**2 * 0.5 * val_jac * &
                  delta_eta1*delta_eta2*delta_eta3*delta_vpar

                nrj_pot = nrj_pot + &
                  delta_f * val_elec_pot * 0.5 * val_jac * &
                  delta_eta1*delta_eta2*delta_eta3*delta_vpar

                ! definition FALSE 
                ! heat_flux = heat_flux + &
                !     delta_f * vpar**2* 0.5* val_jac * &
                !    delta_eta1*delta_eta2*delta_eta3*delta_vpar
              end if
            end if
          end do
        end do
      end do
    end do

    nrj_tot = nrj_kin + nrj_pot

    sim%diag_nrj_kin(diag_num)   = nrj_kin
    sim%diag_nrj_pot(diag_num)   = nrj_pot
    sim%diag_nrj_tot(diag_num)   = nrj_tot
    sim%diag_heat_flux(diag_num) = heat_flux

  end subroutine compute_energy


  !-----------------------------------------------------------
  ! Computation of the L1 norm , i.e   
  !   Norm_L1 = \int abs(delta f) * jac dvpar deta1 deta2 deta3
  !   delta f = f(x,y,z,vpar)
  !   jac     = jacobian of the change of coordinates
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !     
  !  Out: Norm_L1
  !
  ! Computation of the L2 norm , i.e   
  !   Norm_L2 = SQRT(\int abs(delta f)**2 * jac dvpar deta1 deta2 deta3)
  !   delta f = f(x,y,z,vpar)
  !   jac     = jacobian of the change of coordinates
  !  In : f4d_seqx3x4(x1=distrib,x2=distrib,x3=*,x4=*)
  !      
  !  Out: Norm_L2
  !
  ! Computation of the L infini = max( abs( delta f) )
  !-----------------------------------------------------------
  subroutine compute_norm_L1_L2_Linf(sim,diag_num)
    
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    sll_int32                         , intent(in)    :: diag_num
    
    ! local variables
    sll_int32  :: Neta1_loc,Neta2_loc,Neta3,Nvpar
    sll_int32  :: Neta1, Neta2
    sll_real64 :: delta_eta1,delta_eta2,delta_eta3,delta_vpar
    sll_real64 :: val_jac
    sll_real64 :: delta_f
    sll_int32  :: iloc1, iloc2
    sll_int32  :: i1,i2,i3,i4
    sll_int32, dimension(1:4) :: glob_ind4d
    sll_real64 :: mass, norm_L1,norm_L2,norm_Linf,entropy_kin

    Neta1_loc  = size(sim%f4d_seqx3x4,1)
    Neta2_loc  = size(sim%f4d_seqx3x4,2)
    Neta1      = sim%nc_x1 + 1
    Neta2      = sim%nc_x2 + 1
    Neta3      = size(sim%f4d_seqx3x4,3)
    Nvpar      = size(sim%f4d_seqx3x4,4)
    delta_eta1 = sim%logical_mesh4d%delta_eta1
    delta_eta2 = sim%logical_mesh4d%delta_eta2
    delta_eta3 = sim%logical_mesh4d%delta_eta3
    delta_vpar = sim%logical_mesh4d%delta_eta4
    
    norm_L1    = 0.0
    norm_L2    = 0.0
    norm_Linf  = 0.0
    mass       = 0.0
    entropy_kin= 0.0

    !-> Computation of the enrgy kinetic locally in (x1,x2) directions
    do iloc2 = 1,Neta2_loc
      do iloc1 = 1,Neta1_loc
        do i3 = 1,Neta3-1
          do i4 = 1,Nvpar-1

            glob_ind4d(:) = local_to_global_4D(sim%layout4d_seqx3x4, &
              (/iloc1,iloc2,i3,i4/))
            i1 = glob_ind4d(1)
            i2 = glob_ind4d(2)

            if (i1 .ne. Neta1) then
              if (i2 .ne. Neta2) then 

                val_jac = abs(sim%transf_xy%jacobian_at_node(i1,i2))

                delta_f = sim%f4d_seqx3x4(iloc1,iloc2,i3,i4)

                mass   = mass + &
                  delta_f * val_jac * &
                  delta_eta1*delta_eta2*delta_eta3*delta_vpar

                norm_L1 = norm_L1 + &
                  abs(delta_f) * val_jac * &
                  delta_eta1*delta_eta2*delta_eta3*delta_vpar

                norm_L2 = norm_L2 + &
                  abs(delta_f)**2 * val_jac * &
                  delta_eta1*delta_eta2*delta_eta3*delta_vpar

                if ( delta_f /= 0.0_f64) &
                  entropy_kin = entropy_kin - &
                  delta_f* log(abs(delta_f)) * val_jac * &
                  delta_eta1*delta_eta2*delta_eta3*delta_vpar

                norm_Linf = max(abs(delta_f),norm_Linf)
              end if
            end if
          end do
        end do
      end do
    end do

    norm_L2   = sqrt(norm_L2)

    sim%diag_mass(diag_num)        = mass
    sim%diag_norm_L1(diag_num)     = norm_L1
    sim%diag_norm_L2(diag_num)     = norm_L2
    sim%diag_norm_Linf(diag_num)   = norm_Linf
    sim%diag_entropy_kin(diag_num) = entropy_kin

  end subroutine compute_norm_L1_L2_Linf

end module sll_simulation_4d_DK_hybrid_module
