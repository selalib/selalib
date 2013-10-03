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

  implicit none

#define PRINT_PLOTS 1
  type, extends(sll_simulation_base_class) :: sll_simulation_4d_DK_hybrid
     ! Parallel environment parameters
     sll_int32  :: world_size
     sll_int32  :: my_rank
     sll_int32  :: power2 ! 2^power2 = number of processes available
     ! Processor mesh sizes
     sll_int32  :: nproc_x1
     sll_int32  :: nproc_x2
     sll_int32  :: nproc_x3
     sll_int32  :: nproc_x4 
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
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     sll_int32  :: spline_degree_eta1, spline_degree_eta2
     sll_int32  :: spline_degree_eta3, spline_degree_eta4
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

     !--> 4D logical mesh (eta1,eta2,eta3,vpar)
     sll_int32 :: Neta1, Neta2, Neta3, Neta4
     type(sll_logical_mesh_4d), pointer :: logical_mesh4d

     !--> Coordinate transformation
     class(sll_coordinate_transformation_2d_base), pointer :: transf_xy

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
     sll_real64, dimension(:)  , pointer :: n0_r
     sll_real64, dimension(:)  , pointer :: Ti_r
     sll_real64, dimension(:)  , pointer :: Te_r
     sll_real64, dimension(:,:), pointer :: n0_xy
     sll_real64, dimension(:,:), pointer :: Ti_xy
     sll_real64, dimension(:,:), pointer :: Te_xy

     !--> Equilibrium distribution function
     sll_real64, dimension(:,:,:), pointer :: feq_xyvpar

     !--> 4D distribution function 
     !----> sequential in (x1,x2) and parallel in (x3,x4)
     type(layout_4D), pointer :: layout4d_x1x2
     sll_real64, dimension(:,:,:,:), pointer :: f4d_x1x2 
     !----> parallel in (x1,x2) and sequential in (x3,x4) 
     type(layout_4D), pointer :: layout4d_x3x4
     sll_real64, dimension(:,:,:,:), pointer :: f4d_x3x4

     !--> 3D charge density and 3D electric potential
     !----> sequential in (x1,x2)
     type(layout_3D), pointer :: layout3d_x1x2
     sll_real64, dimension(:,:,:), pointer :: rho3d_x1x2 
     sll_real64, dimension(:,:,:), pointer :: phi3d_x1x2 
     !----> sequential in x3
     type(layout_3D), pointer :: layout3d_x3
     sll_real64, dimension(:,:,:), pointer :: rho3d_x3
     sll_real64, dimension(:,:,:), pointer :: phi3d_x3

     !--> For general QN solver
     ! interpolation any arbitrary spline
     type(arb_deg_2d_interpolator) :: interp_rho2d
     type(arb_deg_2d_interpolator) :: interp_phi2d
     type(arb_deg_2d_interpolator) :: interp_QN_A11
     type(arb_deg_2d_interpolator) :: interp_QN_A12
     type(arb_deg_2d_interpolator) :: interp_QN_A21
     type(arb_deg_2d_interpolator) :: interp_QN_A22
     type(arb_deg_2d_interpolator) :: interp_QN_C
     type(sll_scalar_field_2d_discrete_alt), pointer :: rho2d
     type(sll_scalar_field_2d_discrete_alt), pointer :: phi2d
     type(sll_scalar_field_2d_discrete_alt), pointer :: QN_A11 
     type(sll_scalar_field_2d_discrete_alt), pointer :: QN_A12
     type(sll_scalar_field_2d_discrete_alt), pointer :: QN_A21
     type(sll_scalar_field_2d_discrete_alt), pointer :: QN_A22
     type(sll_scalar_field_2d_discrete_alt), pointer :: QN_C

   contains
     procedure, pass(sim) :: run => run_4d_DK_hybrid
     procedure, pass(sim) :: init_from_file => init_4d_DK_hybrid
  end type sll_simulation_4d_DK_hybrid

  interface delete
     module procedure delete_4d_DK_hybrid
  end interface delete

  interface initialize
     module procedure initialize_4d_DK_hybrid
  end interface initialize

contains

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

    namelist /mesh/ num_cells_x1, num_cells_x2, &
      num_cells_x3, num_cells_x4, &
      r_min, r_max, phi_min, phi_max, &
      vpar_min, vpar_max
    namelist /equilibrium/ tau0, rho_peak, kappan, deltarn, &
      kappaTi, deltarTi, kappaTe, deltarTe
    namelist /perturbation/ perturb_choice, mmode, nmode, eps_perturb
    namelist /sim_params/ dt, number_iterations, spline_degree

    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_4d_DK_hybrid() failed to open file ', filename
       STOP
    end if
    read(input_file,mesh)
    read(input_file,equilibrium)
    read(input_file,perturbation)
    read(input_file,sim_params)
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
    sim%vpar_min = vpar_max
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
    sim%num_iterations     = number_iterations
    sim%spline_degree_eta1 = spline_degree
    sim%spline_degree_eta2 = spline_degree
    sim%spline_degree_eta3 = spline_degree
    sim%spline_degree_eta4 = spline_degree
  end subroutine init_4d_DK_hybrid


  !----------------------------------------------------
  ! Initialization of the logical mesh associated
  !  to the 4D DK simulation
  !----------------------------------------------------
  subroutine init_profiles_DK( sim )
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32  :: ierr
    sll_int32  :: ir, Nr, Nx, Ny
    sll_real64 :: Lr, dr
    sll_real64 :: r_peak, n0_rmin
    sll_real64 :: Ti_rmin, Te_rmin, Ti_scal, Te_scal
    sll_real64, dimension(:), pointer :: r_grid_tmp

    !--> Initialization of r_grid
    Nr = sim%nc_x1+1
    SLL_ALLOCATE(r_grid_tmp(Nr),ierr)
    Lr = abs(sim%r_max-sim%r_min)
    dr = Lr/float(Nr)
    do ir = 1,Nr
      r_grid_tmp(ir) = sim%r_min + float(ir-1)*dr
    end do

    !--> Initialization of n0(r), Ti(r) and Te(r)
    SLL_ALLOCATE(sim%n0_r(Nr),ierr)
    SLL_ALLOCATE(sim%Ti_r(Nr),ierr)
    SLL_ALLOCATE(sim%Te_r(Nr),ierr)
    sim%n0_r(:) = 0.0_f64
    sim%Ti_r(:) = 0.0_f64
    sim%Te_r(:) = 0.0_f64

    n0_rmin = 10._f64**19
    Ti_rmin = 1.3_f64
    Ti_scal = 1._f64
    
    r_peak  = sim%r_min + sim%rho_peak * Lr
    call init_n0_r(r_peak,sim%kappan, &
      sim%deltarn,n0_rmin,r_grid_tmp,sim%n0_r)
    call init_T_r(r_peak,sim%kappaTi, &
      sim%deltarTi,Ti_rmin,Ti_scal,r_grid_tmp,sim%Ti_r)
    Te_rmin = sim%Ti_r(1)
    Te_scal = sim%tau0
    call init_T_r(r_peak,sim%kappaTe, &
      sim%deltarTe,Te_rmin,Te_scal,r_grid_tmp,sim%Te_r)

    !--> Initialization of n0(x,y), Ti(x,y) and Te(x,y)
    Nx = sim%Neta1
    Ny = sim%Neta2
    SLL_ALLOCATE(sim%n0_xy(Nx,Ny),ierr)
    SLL_ALLOCATE(sim%Ti_xy(Nx,Ny),ierr)
    SLL_ALLOCATE(sim%Te_xy(Nx,Ny),ierr)
    call function_xy_from_r(r_grid_tmp,sim%n0_r,sim%xgrid_2d, &
      sim%ygrid_2d,sim%n0_xy)
    call function_xy_from_r(r_grid_tmp,sim%Ti_r,sim%xgrid_2d, &
      sim%ygrid_2d,sim%Ti_xy)
    call function_xy_from_r(r_grid_tmp,sim%Te_r,sim%xgrid_2d, &
      sim%ygrid_2d,sim%Te_xy)

    SLL_DEALLOCATE(r_grid_tmp,ierr)
  end subroutine init_profiles_DK


  !----------------------------------------------------
  ! Allocation of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine allocate_fdistribu4d_DK( sim )
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32 :: ierr, itemp
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4

    ! layout for sequential operations in x3 and x4. 
    ! Make an even split for x1 and x2, or as close as 
    ! even if the power of 2 is odd. This should 
    ! be packaged in some sort of routine and set up 
    ! at initialization time.
    sim%power2 = int(log(real(sim%world_size))/log(2.0))
    !--> special case N = 1, so power2 = 0
    if(sim%power2 == 0) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    end if
    
    if(is_even(sim%power2)) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 2**(sim%power2/2)
       sim%nproc_x4 = 2**(sim%power2/2)
    else 
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 2**((sim%power2-1)/2)
       sim%nproc_x4 = 2**((sim%power2+1)/2)
    end if

    !--> Initialization of parallel layout of f4d in (x3,x4) directions
    !-->  (x1,x2) : sequential
    !-->  (x3,x4) : parallelized layout
    sim%layout4d_x1x2  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%nc_x1+1, & 
      sim%nc_x2+1, & 
      sim%nc_x3+1, &
      sim%nc_x4+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_x1x2 )
    
    ! Allocate the array needed to store the local chunk 
    ! of the distribution function data. First compute the 
    ! local sizes. Since the remap operations
    ! are out-of-place, we will allocate two different arrays, 
    ! one for each layout.
    call compute_local_sizes_4d( sim%layout4d_x1x2, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    SLL_ALLOCATE(sim%f4d_x1x2(loc4d_sz_x1,loc4d_sz_x2, loc4d_sz_x3,loc4d_sz_x4),ierr)

    !--> Initialization of parallel layout of f4d in (x1,x2) directions
    !-->  (x1,x2) : parallelized layout
    !-->  (x3,x4) : sequential
    ! switch x3 and x1:
    itemp        = sim%nproc_x1
    sim%nproc_x1 = sim%nproc_x3
    sim%nproc_x3 = itemp
    ! switch x4 and x2
    itemp        = sim%nproc_x2
    sim%nproc_x2 = sim%nproc_x4 
    sim%nproc_x4 = itemp

    sim%layout4d_x3x4  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%nc_x1+1, &
      sim%nc_x2+1, &
      sim%nc_x3+1, &
      sim%nc_x4+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_x3x4 )
        
    call compute_local_sizes_4d( sim%layout4d_x3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    
    SLL_ALLOCATE(sim%f4d_x3x4(loc4d_sz_x1,loc4d_sz_x2, loc4d_sz_x3,loc4d_sz_x4),ierr)
  end subroutine allocate_fdistribu4d_DK


  !----------------------------------------------------
  ! Initialization of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine initialize_fdistribu4d_DK(sim)
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32  :: ierr
    sll_int32  :: i1, i2, i3, i4
    sll_int32  :: iloc1, iloc2, iloc3, iloc4
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_int32  :: Nx, Ny, Nvpar
    sll_real64 :: theta_j, phi_k
    sll_int32, dimension(1:4) :: glob_ind

    sll_real64 :: Lphi, dvpar, Lvpar
    sll_real64, dimension(:), pointer :: vpar_grid_tmp

    Nx    = sim%Neta1
    Ny    = sim%Neta2
    Nvpar = sim%Neta4
    
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

    !--> Initialization of the distribution function f4d_x3x4
    call compute_local_sizes_4d( sim%layout4d_x3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )

    Lphi = abs(sim%phi_max-sim%phi_min)
    do iloc4 = 1,loc4d_sz_x4
      do iloc3 = 1,loc4d_sz_x3
        do iloc2 = 1,loc4d_sz_x2
          do iloc1 = 1,loc4d_sz_x1
            glob_ind(:) = local_to_global_4D(sim%layout4d_x3x4, &
              (/iloc1,iloc2,iloc3,iloc4/))
            i1 = glob_ind(1)
            i2 = glob_ind(2)
            i3 = glob_ind(3)
            i4 = glob_ind(4)
!VG!            theta_j = sim%eta2_grid(i2)
!VG!            phi_k = sim%eta3_grid(i3)
            sim%f4d_x3x4(iloc1,iloc2,i3,i4) = &
              sim%feq_xyvpar(i1,i2,i4) * &
              (1._f64+sim%eps_perturb*cos(real(sim%mmode)*theta_j + &
              2._f64*sll_pi*real(sim%nmode) / Lphi*phi_k))
          end do
        end do
      end do
    end do
    SLL_DEALLOCATE(vpar_grid_tmp,ierr)
  end subroutine initialize_fdistribu4d_DK


  !----------------------------------------------------
  ! Allocation for QN solver
  !----------------------------------------------------
  subroutine allocate_QN_DK( sim )
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    type(sll_logical_mesh_2d), pointer :: logical_mesh2d

    sll_int32 :: ierr, itemp
    sll_int32 :: i1, i2, i3, i4
    sll_int32 :: iloc1, iloc2, iloc3, iloc4
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3
    sll_int32 :: nproc3d_x3

    ! layout for sequential operations in x3 
    sim%power2 = int(log(real(sim%world_size))/log(2.0))
    !--> special case N = 1, so power2 = 0
    if(sim%power2 == 0) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    end if
    
    if(is_even(sim%power2)) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 2**(sim%power2/2)
       sim%nproc_x4 = 2**(sim%power2/2)
    else 
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 2**((sim%power2-1)/2)
       sim%nproc_x4 = 2**((sim%power2+1)/2)
    end if

    !--> Initialization of rho3d_x1x2 and phi3d_x1x2
    !-->  (x1,x2) : sequential
    !-->  x3 : parallelized layout    
    sim%layout3d_x1x2  => new_layout_3D( sll_world_collective )
    nproc3d_x3 = sim%nproc_x3*sim%nproc_x4
    call initialize_layout_with_distributed_3D_array( &
      sim%nc_x1+1, & 
      sim%nc_x2+1, & 
      sim%nc_x3+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      nproc3d_x3, &
      sim%layout3d_x1x2 )
    call compute_local_sizes_3d( sim%layout3d_x1x2, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE(sim%rho3d_x1x2(loc3d_sz_x1,loc3d_sz_x2, loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_x1x2(loc3d_sz_x1,loc3d_sz_x2, loc3d_sz_x3),ierr)

    !--> Initialization of rho3d_x3 and phi3d_x3
    !-->  (x1,x2) : parallelized layout
    !-->  x3 : sequential
    ! switch x3 and x1:
    itemp        = sim%nproc_x1
    sim%nproc_x1 = sim%nproc_x3
    sim%nproc_x3 = itemp
    ! switch x4 and x2
    itemp        = sim%nproc_x2
    sim%nproc_x2 = sim%nproc_x4 
    sim%nproc_x4 = itemp
        
    sim%layout3d_x3  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( &
      sim%nc_x1+1, &
      sim%nc_x2+1, &
      sim%nc_x3+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%layout3d_x3 )
    call compute_local_sizes_3d( sim%layout3d_x3, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE(sim%rho3d_x3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_x3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    
    !---->
    logical_mesh2d => sim%transf_xy%mesh
    
    call sim%interp_phi2d%initialize( &
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
print*,sim%my_rank,"===> OK"

!VG!    call sim%interp_rho2d%initialize( &
!VG!      logical_mesh2d%num_cells1 +1, &
!VG!      logical_mesh2d%num_cells2 +1, &
!VG!      logical_mesh2d%eta1_min, &
!VG!      logical_mesh2d%eta1_max, &
!VG!      logical_mesh2d%eta2_min, &
!VG!      logical_mesh2d%eta2_max, &
!VG!      sim%bc_left_eta1, &
!VG!      sim%bc_right_eta1, &
!VG!      sim%bc_left_eta2, &
!VG!      sim%bc_right_eta2, &
!VG!      sim%spline_degree_eta1, &
!VG!      sim%spline_degree_eta2)    
!VG!print*,sim%my_rank,"===> OK1"
!VG!
!VG!    call sim%interp_QN_A11%initialize( &
!VG!      logical_mesh2d%num_cells1 +1, &
!VG!      logical_mesh2d%num_cells2 +1, &
!VG!      logical_mesh2d%eta1_min, &
!VG!      logical_mesh2d%eta1_max, &
!VG!      logical_mesh2d%eta2_min, &
!VG!      logical_mesh2d%eta2_max, &
!VG!      sim%bc_left_eta1, &
!VG!      sim%bc_right_eta1, &
!VG!      sim%bc_left_eta2, &
!VG!      sim%bc_right_eta2, &
!VG!      sim%spline_degree_eta1, &
!VG!      sim%spline_degree_eta2)    
!VG!print*,sim%my_rank,"===> OK2"
!VG!
!VG!    call sim%interp_QN_A12%initialize( &
!VG!      logical_mesh2d%num_cells1 +1, &
!VG!      logical_mesh2d%num_cells2 +1, &
!VG!      logical_mesh2d%eta1_min, &
!VG!      logical_mesh2d%eta1_max, &
!VG!      logical_mesh2d%eta2_min, &
!VG!      logical_mesh2d%eta2_max, &
!VG!      sim%bc_left_eta1, &
!VG!      sim%bc_right_eta1, &
!VG!      sim%bc_left_eta2, &
!VG!      sim%bc_right_eta2, &
!VG!      sim%spline_degree_eta1, &
!VG!      sim%spline_degree_eta2)    
!VG!print*,sim%my_rank,"===> OK3"
!VG!
!VG!    call sim%interp_QN_A21%initialize( &
!VG!      logical_mesh2d%num_cells1 +1, &
!VG!      logical_mesh2d%num_cells2 +1, &
!VG!      logical_mesh2d%eta1_min, &
!VG!      logical_mesh2d%eta1_max, &
!VG!      logical_mesh2d%eta2_min, &
!VG!      logical_mesh2d%eta2_max, &
!VG!      sim%bc_left_eta1, &
!VG!      sim%bc_right_eta1, &
!VG!      sim%bc_left_eta2, &
!VG!      sim%bc_right_eta2, &
!VG!      sim%spline_degree_eta1, &
!VG!      sim%spline_degree_eta2)    
!VG!print*,sim%my_rank,"===> OK4"
!VG!
!VG!    call sim%interp_QN_A22%initialize( &
!VG!      logical_mesh2d%num_cells1 +1, &
!VG!      logical_mesh2d%num_cells2 +1, &
!VG!      logical_mesh2d%eta1_min, &
!VG!      logical_mesh2d%eta1_max, &
!VG!      logical_mesh2d%eta2_min, &
!VG!      logical_mesh2d%eta2_max, &
!VG!      sim%bc_left_eta1, &
!VG!      sim%bc_right_eta1, &
!VG!      sim%bc_left_eta2, &
!VG!      sim%bc_right_eta2, &
!VG!      sim%spline_degree_eta1, &
!VG!      sim%spline_degree_eta2)    
!VG!print*,sim%my_rank,"===> OK5"
!VG!
!VG!    call sim%interp_QN_C%initialize( &
!VG!      logical_mesh2d%num_cells1 +1, &
!VG!      logical_mesh2d%num_cells2 +1, &
!VG!      logical_mesh2d%eta1_min, &
!VG!      logical_mesh2d%eta1_max, &
!VG!      logical_mesh2d%eta2_min, &
!VG!      logical_mesh2d%eta2_max, &
!VG!      sim%bc_left_eta1, &
!VG!      sim%bc_right_eta1, &
!VG!      sim%bc_left_eta2, &
!VG!      sim%bc_right_eta2, &
!VG!      sim%spline_degree_eta1, &
!VG!      sim%spline_degree_eta2)    
!VG!print*,sim%my_rank,"===> OK6"
    !----->
!!$    sim%rho2d => new_scalar_field_2d_discrete_alt( &
!!$      sim%rho_full - density_tot, &
!!$      "rho_field_check", &
!!$      sim%interp_rho, &     
!!$      sim%transf_xy, &
!!$      sim%bc_left_eta1, &
!!$      sim%bc_right_eta1, &
!!$      sim%bc_left_eta2, &
!!$      sim%bc_right_eta2)
  end subroutine allocate_QN_DK


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
    sll_int32 :: ieta1, ieta2
    sll_int32 :: Nx, Ny

    !--> Parallelization initialization
    sim%world_size = world_size
    sim%my_rank    = my_rank

    !--> Initialization of the number of points
    sim%Neta1 = sim%nc_x1+1
    sim%Neta2 = sim%nc_x2+1
    sim%Neta3 = sim%nc_x3+1
    sim%Neta4 = sim%nc_x4+1

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
    
    !--> Transformation initialization
    sim%transf_xy => transf_xy

    !--> Initialization of the (x,y) 2D mesh
    Nx = sim%Neta1
    Ny = sim%Neta2
    SLL_ALLOCATE(sim%xgrid_2d(Nx,Ny),ierr)
    SLL_ALLOCATE(sim%ygrid_2d(Nx,Ny),ierr)
    do ieta2 = 1,sim%Neta2
      do ieta1 = 1,sim%Neta1
        sim%xgrid_2d(ieta1,ieta2) = &
          sim%transf_xy%x1_at_node(ieta1,ieta2)
        sim%ygrid_2d(ieta1,ieta2) = &
          sim%transf_xy%x2_at_node(ieta1,ieta2)
      end do
    end do
      
!VG!    !--> Radial profile initialisation
!VG!    call init_profiles_DK(sim)

    !*** Allocation of the distribution function ***
    call allocate_fdistribu4d_DK(sim)
    
    !*** Allocation of the QN solver ***
    call allocate_QN_DK(sim)
  end subroutine initialize_4d_DK_hybrid


  !----------------------------------------------------
  ! First step of the drift-kinetic 4D simulation :
  !  1) Initialization of the distribution function 
  !  2) Computation of the charge density 
  !      (r.h.s of the quasi-neutrality solver)
  !  3) Solving of the quasi-neutrality equation
  !----------------------------------------------------
  subroutine first_step_4d_DK_hybrid( sim )
    use sll_hdf5_io, only: sll_hdf5_file_create, &
      sll_hdf5_write_array_1d, sll_hdf5_file_close
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    !--> For initial profile HDF5 saving
    integer                      :: file_err
    sll_int32                    :: file_id
    character(len=12), parameter :: filename_prof = "init_prof.h5"

    !*** Initialization of the distribution function ***
    call initialize_fdistribu4d_DK(sim)

    !*** Saving of the radial profiles in HDF5 file ***
    if (sim%my_rank.eq.0) then
      call sll_hdf5_file_create(filename_prof,file_id,file_err)
      call sll_hdf5_write_array_1d(file_id,sim%n0_r,'n0_r',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Ti_r,'Ti_r',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Te_r,'Te_r',file_err)
!VG!      call sll_hdf5_write_array_2d(file_id,sim%feq_2d,'feq_2d',file_err)
      call sll_hdf5_file_close(file_id,file_err)
    end if

    !*** Computation of the rhs of QN ***
    call compute_charge_density(sim%logical_mesh4d, &
      sim%f4d_x3x4,sim%rho3d_x3)
  end subroutine first_step_4d_DK_hybrid


  !----------------------------------------------------
  ! Run drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine run_4d_DK_hybrid( sim )
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

!VG!    !-->
!VG!    sll_real64, dimension(:), pointer :: send_buf
!VG!    sll_real64, dimension(:), pointer :: recv_buf
!VG!    sll_int32 , dimension(:), pointer :: recv_sz
!VG!    sll_int32 , dimension(:), pointer :: disps ! for allgatherv operation
!VG!
!VG!    !--> Initialization of the 4D parallel layout
!VG!    SLL_ALLOCATE(recv_sz(sim%world_size),ierr)
!VG!    SLL_ALLOCATE(disps(sim%world_size),ierr)
!VG!
!VG!    !--> Deallocation of temparory array
!VG!    SLL_DEALLOCATE(recv_sz,ierr)
!VG!    SLL_DEALLOCATE(disps,ierr)
  end subroutine run_4d_DK_hybrid


  !-----------------------------------------------------------
  ! Computation of the charge density, i.e
  !  rho(eta1,eta2,eta3) = \int f(eta1,eta2,eta3,vpar) dvpar
  !  In : f4d_x3x4(x1 part,x2 part,x3=*,x4=*)
  !  Out: rho3d(x1=*,x2=*,x3=*)
  !-----------------------------------------------------------
  subroutine compute_charge_density(logical_mesh4d,f4d_x3x4,rho3d_x3)
    type(sll_logical_mesh_4d)     , intent(in)    :: logical_mesh4d
    sll_real64, dimension(:,:,:,:), intent(in)    :: f4d_x3x4
    sll_real64, dimension(:,:,:)  , intent(inout) :: rho3d_x3

    sll_int32  :: Neta1_loc,Neta2_loc,Neta3, Nvpar
    sll_int32  :: iloc1, iloc2, i3, i4
    sll_real64 :: delta_vpar, intf_dvpar
    
    Neta1_loc  = size(f4d_x3x4,1)
    Neta2_loc  = size(f4d_x3x4,2)
    Neta3      = size(f4d_x3x4,3)
    Nvpar      = size(f4d_x3x4,4)
    delta_vpar = logical_mesh4d%delta_eta4 

    !-> Computation of the charge density locally in (x1,x2) directions
    do i3 = 1,Neta3
      do iloc2 = 1,Neta2_loc
        do iloc1 = 1,Neta1_loc
          intf_dvpar = 0._f64
          do i4 = 1,Nvpar
            intf_dvpar = intf_dvpar + &
              f4d_x3x4(iloc1,iloc2,i3,i4)*delta_vpar
          end do
          rho3d_x3(iloc1,iloc2,i3) = intf_dvpar
        end do
      end do
    end do
  end subroutine compute_charge_density
  

  !----------------------------------------------------
  ! Initialization of the drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine delete_4d_DK_hybrid( sim )
    class(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32 :: ierr

    SLL_DEALLOCATE(sim%n0_r,ierr)
    SLL_DEALLOCATE(sim%Ti_r,ierr)
    SLL_DEALLOCATE(sim%Te_r,ierr)
    SLL_DEALLOCATE(sim%feq_xyvpar,ierr)
    SLL_DEALLOCATE(sim%f4d_x1x2,ierr)
    SLL_DEALLOCATE(sim%f4d_x3x4,ierr)
    call delete(sim%layout4d_x1x2)
    call delete(sim%layout4d_x3x4)
  end subroutine delete_4d_DK_hybrid

end module sll_simulation_4d_DK_hybrid_module
