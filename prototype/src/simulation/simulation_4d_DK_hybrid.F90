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
!VG!  use sll_constants
!VG!  use sll_cubic_spline_interpolator_2d
!VG!  use sll_poisson_2d_periodic_cartesian_par
!VG!  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  use sll_logical_meshes
  use sll_fdistribu4D_DK
!VG!  use sll_coordinate_transformation_2d_base_module
!VG!  use sll_gnuplot_parallel
!VG!  use sll_general_coordinate_qn_solver_module
!VG!  use sll_module_scalar_field_2d_base
!VG!  use sll_module_scalar_field_2d_alternative

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
     sll_real64 :: eta1_min
     sll_real64 :: eta2_min
     sll_real64 :: eta3_min
     sll_real64 :: vpar_min
     sll_real64 :: Leta1
     sll_real64 :: Leta2
     sll_real64 :: Leta3
     sll_real64 :: Lvpar
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

     !--> 4D logical mesh (eta1,eta2,eta3,v)
     type(sll_logical_mesh_4d), pointer :: logical_mesh4d
     sll_real64, dimension(:), pointer  :: eta1_grid
     sll_real64, dimension(:), pointer  :: eta2_grid
     sll_real64, dimension(:), pointer  :: eta3_grid
     sll_real64, dimension(:), pointer  :: vpar_grid

     !--> For boundary conditions
     sll_int32 :: bc_left_eta1
     sll_int32 :: bc_right_eta1
     sll_int32 :: bc_left_eta2
     sll_int32 :: bc_right_eta2
     sll_int32 :: bc_left_eta3
     sll_int32 :: bc_right_eta3
     sll_int32 :: bc_left_vpar
     sll_int32 :: bc_right_vpar

     !--> Parallel decomposition of the mesh

     !--> Density and temperature profiles
     sll_real64, dimension(:), pointer :: n0
     sll_real64, dimension(:), pointer :: Ti
     sll_real64, dimension(:), pointer :: Te

     !--> Equilibrium distribution function
     sll_real64, dimension(:,:), pointer :: feq_2d

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
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta3_min
    sll_real64 :: vpar_min
    sll_real64 :: Leta1
    sll_real64 :: Leta2
    sll_real64 :: Leta3
    sll_real64 :: Lvpar
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
      eta1_min, eta2_min, eta3_min, vpar_min, &
      Leta1, Leta2, Leta3, Lvpar
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
    sim%eta1_min = eta1_min
    sim%eta2_min = eta2_min
    sim%eta3_min = eta3_min
    sim%vpar_min = vpar_min
    sim%Leta1    = Leta1
    sim%Leta2    = Leta2
    sim%Leta3    = Leta3
    sim%Lvpar    = Lvpar
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
  subroutine init_logicalmesh_DK( sim, logical_mesh4d)
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    type(sll_logical_mesh_4d)        , pointer :: logical_mesh4d

    sll_int32  :: ierr
    sll_int32  :: ieta1, ieta2, ieta3, ivpar
    sll_int32  :: Neta1, Neta2, Neta3, Nvpar
    sll_real64 :: eta1_min, eta2_min, eta3_min, vpar_min
    sll_real64 :: delta_eta1, delta_eta2, delta_eta3 , delta_vpar

    sim%logical_mesh4d => logical_mesh4d

    !----> Init eta1 grid
    Neta1      = sim%logical_mesh4d%num_cells1+1
    SLL_ALLOCATE(sim%eta1_grid(Neta1),ierr)
    delta_eta1 = sim%logical_mesh4d%delta_eta1
    eta1_min   = sim%logical_mesh4d%eta1_min
    do ieta1 = 1,Neta1
      sim%eta1_grid(ieta1) = eta1_min + (ieta1-1)*delta_eta1
    end do

    !----> Init eta2 grid
    Neta2      = sim%logical_mesh4d%num_cells2+1
    SLL_ALLOCATE(sim%eta2_grid(Neta2),ierr)
    delta_eta2 = sim%logical_mesh4d%delta_eta2
    eta2_min   = sim%logical_mesh4d%eta2_min
    do ieta2 = 1,Neta2
      sim%eta2_grid(ieta2) = eta2_min + (ieta2-1)*delta_eta2
    end do

    !----> Init eta3 grid
    Neta3      = sim%logical_mesh4d%num_cells3+1
    SLL_ALLOCATE(sim%eta3_grid(Neta3),ierr)
    delta_eta3 = sim%logical_mesh4d%delta_eta3
    eta3_min   = sim%logical_mesh4d%eta3_min
    do ieta3 = 1,Neta3
      sim%eta3_grid(ieta3) = eta3_min + (ieta3-1)*delta_eta3
    end do

    !----> Init vpar grid
    Nvpar      = sim%logical_mesh4d%num_cells4+1
    SLL_ALLOCATE(sim%vpar_grid(Nvpar),ierr)
    delta_vpar = sim%logical_mesh4d%delta_eta4
    vpar_min   = sim%logical_mesh4d%eta4_min
    do ivpar = 1,Nvpar
      sim%vpar_grid(ivpar) = vpar_min + (ivpar-1)*delta_vpar
    end do
  end subroutine init_logicalmesh_DK


  !----------------------------------------------------
  ! Initialization of the logical mesh associated
  !  to the 4D DK simulation
  !----------------------------------------------------
  subroutine init_profiles_DK( sim )
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32  :: ierr
    sll_real64 :: eta1_peak
    sll_real64 :: n0_eta1min, Ti_eta1min, Te_eta1min, Ti_scal, Te_scal

    SLL_ALLOCATE(sim%n0(sim%nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Ti(sim%nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Te(sim%nc_x1+1),ierr)

    sim%n0(:) = 0.0_f64
    sim%Ti(:) = 0.0_f64
    sim%Te(:) = 0.0_f64

    n0_eta1min = 10._f64**19
    Ti_eta1min = 1.3_f64
    Ti_scal    = 1._f64
    eta1_peak  = sim%logical_mesh4d%eta1_min + sim%rho_peak * &
      abs(sim%logical_mesh4d%eta1_max-sim%logical_mesh4d%eta1_min)
    call init_n0_eta1(eta1_peak,sim%kappan, &
      sim%deltarn,n0_eta1min,sim%eta1_grid,sim%n0)
    call init_T_eta1(eta1_peak,sim%kappaTi, &
      sim%deltarTi,Ti_eta1min,Ti_scal,sim%eta1_grid,sim%Ti)
    Te_eta1min = sim%Ti(1)
    Te_scal    = sim%tau0
    call init_T_eta1(eta1_peak,sim%kappaTe, &
      sim%deltarTe,Te_eta1min,Te_scal,sim%eta1_grid,sim%Te)
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
    SLL_ALLOCATE(sim%f4d_x1x2(loc4d_sz_x1,loc4d_sz_x2, &
      loc4d_sz_x3,loc4d_sz_x4),ierr)

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
    SLL_ALLOCATE(sim%f4d_x3x4(loc4d_sz_x1,loc4d_sz_x2, &
      loc4d_sz_x3,loc4d_sz_x4),ierr)
  end subroutine allocate_fdistribu4d_DK


  !----------------------------------------------------
  ! Initialization of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine initialize_fdistribu4d_DK(sim)
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

    sll_int32  :: ierr
    sll_int32  :: Neta1, Neta3, Nvpar
    sll_int32  :: i1, i2, i3, i4
    sll_int32  :: iloc1, iloc2, iloc3, iloc4
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_real64 :: eta2_j, eta3_k
    sll_int32, dimension(1:4) :: glob_ind

    !--> Initialization of the equilibrium distribution function
    Neta1 = sim%nc_x1+1
    Neta3 = sim%nc_x3+1
    Nvpar = sim%nc_x4+1
    SLL_ALLOCATE(sim%feq_2d(Neta1,Nvpar),ierr)
    call init_fequilibrium(Neta1,Nvpar,sim%eta1_grid,sim%vpar_grid, &
      sim%n0,sim%Ti,sim%feq_2d)

    !--> Initialization of the distribution function f4d_x3x4
    call compute_local_sizes_4d( sim%layout4d_x3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )

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
            eta2_j = sim%eta2_grid(i2)
            eta3_k = sim%eta3_grid(i3)
            sim%f4d_x3x4(iloc1,iloc2,i3,i4) = sim%feq_2d(i1,i4) * &
              (1._f64+sim%eps_perturb*cos(real(sim%mmode)*eta2_j + &
              2._f64*sll_pi*real(sim%nmode)/sim%logical_mesh4d%eta3_max*eta3_k))
          end do
        end do
      end do
    end do    
  end subroutine initialize_fdistribu4d_DK


  !----------------------------------------------------
  ! Allocation for QN solver
  !----------------------------------------------------
  subroutine allocate_QN_DK( sim )
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim

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
    SLL_ALLOCATE(sim%rho3d_x1x2(loc3d_sz_x1,loc3d_sz_x2, &
      loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_x1x2(loc3d_sz_x1,loc3d_sz_x2, &
      loc3d_sz_x3),ierr)

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
  end subroutine allocate_QN_DK


  !----------------------------------------------------
  ! Initialization of the drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine initialize_4d_DK_hybrid( sim, &
    world_size, &
    my_rank, &
    logical_mesh4d)
    type(sll_simulation_4d_DK_hybrid), intent(inout) :: sim
    sll_int32                        , intent(in)    :: world_size
    sll_int32                        , intent(in)    :: my_rank
    type(sll_logical_mesh_4d)        , pointer       :: logical_mesh4d

    !--> Parallelization initialization
    sim%world_size = world_size
    sim%my_rank    = my_rank

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
    call init_logicalmesh_DK(sim,logical_mesh4D)

    !--> Radial profile initialisation
    call init_profiles_DK(sim)

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
      call sll_hdf5_write_array_1d(file_id,sim%n0,'n0',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Ti,'Ti',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Te,'Te',file_err)
      call sll_hdf5_write_array_2d(file_id,sim%feq_2d,'feq_2d',file_err)
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
            intf_dvpar = intf_dvpar + f4d_x3x4(iloc1,iloc2,i3,i4)*delta_vpar
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

    SLL_DEALLOCATE(sim%eta1_grid,ierr)
    SLL_DEALLOCATE(sim%vpar_grid,ierr) 
    SLL_DEALLOCATE(sim%n0,ierr)
    SLL_DEALLOCATE(sim%Ti,ierr)
    SLL_DEALLOCATE(sim%Te,ierr)
    SLL_DEALLOCATE(sim%feq_2d,ierr)
    SLL_DEALLOCATE(sim%f4d_x1x2,ierr)
    SLL_DEALLOCATE(sim%f4d_x3x4,ierr)
    call delete(sim%layout4d_x1x2)
    call delete(sim%layout4d_x3x4)
  end subroutine delete_4d_DK_hybrid


!VG!  !----------------------------------------------------
!VG!  ! Initialization of the generalized QN solver
!VG!  !----------------------------------------------------
!VG!  subroutine initialize_QN_general( &
!VG!    sim, &
!VG!    transformation_x, &
!VG!    a11_f, &
!VG!    a12_f, &
!VG!    a21_f, &
!VG!    a22_f, &
!VG!    c_f)
!VG!
!VG!   type(sll_simulation_4d_DK_hybrid), intent(inout)      :: sim
!VG!   class(sll_coordinate_transformation_2d_base), pointer :: transformation_x
!VG!   procedure(two_var_parametrizable_function) :: a11_f
!VG!   procedure(two_var_parametrizable_function) :: a12_f
!VG!   procedure(two_var_parametrizable_function) :: a21_f
!VG!   procedure(two_var_parametrizable_function) :: a22_f
!VG!   procedure(two_var_parametrizable_function) :: c_f
!VG!
!VG!   sim%mesh2d_x  => mesh2d_x
!VG!   sim%mesh2d_v  => mesh2d_v
!VG!   sim%transfx   => transformation_x
!VG!   sim%init_func => init_func
!VG!   sim%params    => params
!VG!   sim%a11_f     => a11_f
!VG!   sim%a12_f     => a12_f
!VG!   sim%a21_f     => a21_f
!VG!   sim%a22_f     => a22_f
!VG!   sim%c_f       => c_f
!VG!   sim%spline_degree_eta1 = spline_degre1
!VG!   sim%spline_degree_eta2 = spline_degre2
!VG!
!VG!   sim%bc_left   = bc_left
!VG!   sim%bc_right  = bc_right
!VG!   sim%bc_bottom = bc_bottom
!VG!   sim%bc_top    = bc_top
!VG!
!VG!   call sim%interp_phi%initialize( &
!VG!        sim%mesh2d_x%num_cells1 +1, &
!VG!        sim%mesh2d_x%num_cells2 +1, &
!VG!        sim%mesh2d_x%eta1_min, &
!VG!        sim%mesh2d_x%eta1_max, &
!VG!        sim%mesh2d_x%eta2_min, &
!VG!        sim%mesh2d_x%eta2_max, &
!VG!        sim%bc_left, &
!VG!        sim%bc_right, &
!VG!        sim%bc_bottom, &
!VG!        sim%bc_top, &
!VG!        sim%spline_degree_eta1, &
!VG!        sim%spline_degree_eta2)
!VG!
!VG!   call sim%interp_rho%initialize( &
!VG!        sim%mesh2d_x%num_cells1 +1, &
!VG!        sim%mesh2d_x%num_cells2 +1, &
!VG!        sim%mesh2d_x%eta1_min, &
!VG!        sim%mesh2d_x%eta1_max, &
!VG!        sim%mesh2d_x%eta2_min, &
!VG!        sim%mesh2d_x%eta2_max, &
!VG!        sim%bc_left, &
!VG!        sim%bc_right, &
!VG!        sim%bc_bottom, &
!VG!        sim%bc_top, &
!VG!        sim%spline_degree_eta1, &
!VG!        sim%spline_degree_eta2)
!VG!  end subroutine initialize_QN_general

end module sll_simulation_4d_DK_hybrid_module
