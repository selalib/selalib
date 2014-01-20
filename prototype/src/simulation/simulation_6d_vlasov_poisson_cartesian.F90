module sll_simulation_6d_vlasov_poisson_cartesian
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
  use sll_collective
  use sll_remapper
  use sll_constants
  use sll_cubic_spline_interpolator_1d
  use sll_distribution_function_6d_initializer
  use sll_poisson_3d_periodic_par
  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  implicit none

  type, extends(sll_simulation_base_class) :: &
     sll_simulation_6d_vlasov_poisson_cart
     
! Parallel environment parameters
     sll_int32  :: world_size
     sll_int32  :: my_rank
     sll_int32  :: power2 ! 2^power2 = number of processes available
     ! Processor mesh sizes
     sll_int32  :: nproc_x1
     sll_int32  :: nproc_x2
     sll_int32  :: nproc_x3
     sll_int32  :: nproc_x4 
     sll_int32  :: nproc_x5
     sll_int32  :: nproc_x6
     ! Physics parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     ! Mesh parameters
     sll_int32  :: nc_x1
     sll_int32  :: nc_x2
     sll_int32  :: nc_x3
     sll_int32  :: nc_x4
     sll_int32  :: nc_x5
     sll_int32  :: nc_x6
     ! for initializers
     type(init_test_6d_par)                      :: init_6d
     type(simple_cartesian_6d_mesh), pointer     :: mesh6d
     type(poisson_3d_periodic_plan_par), pointer :: poisson_plan

     ! distribution functions. There are several because each array represents
     ! a differently shaped chunk of memory. In this example, each chunk 
     ! allows sequential operations in one given direction. f_x1x2x3 should 
     ! permit to carry out sequential operations in x1, x2 and x3 for example.
     sll_real64, dimension(:,:,:,:,:,:), pointer     :: f_x1x2x3 
     sll_real64, dimension(:,:,:,:,:,:), pointer     :: f_x4x5x6
     sll_real64, dimension(:,:,:,:,:), allocatable   :: partial_reduction
     sll_real64, dimension(:,:,:), allocatable       :: rho_x1 
     sll_real64, dimension(:,:,:), allocatable       :: rho_x2 
     sll_real64, dimension(:,:,:), allocatable       :: rho_x3 
     sll_real64, dimension(:,:,:), allocatable       :: rho_split
     sll_real64, dimension(:,:,:), allocatable       :: phi_x1
     sll_real64, dimension(:,:,:), allocatable       :: phi_x2
     sll_real64, dimension(:,:,:), allocatable       :: phi_x3
     sll_real64, dimension(:,:,:), allocatable       :: ex_x1
     sll_real64, dimension(:,:,:), allocatable       :: ex_split
     sll_real64, dimension(:,:,:), allocatable       :: ey_x2
     sll_real64, dimension(:,:,:), allocatable       :: ey_split
     sll_real64, dimension(:,:,:), allocatable       :: ez_x3
     sll_real64, dimension(:,:,:), allocatable       :: ez_split
     ! for remap
     type(layout_6D), pointer :: sequential_x1x2x3
     type(layout_6D), pointer :: sequential_x4x5x6
     type(layout_3D), pointer :: rho_seq_x1
     type(layout_3D), pointer :: rho_seq_x2
     type(layout_3D), pointer :: rho_seq_x3
     type(layout_3D), pointer :: split_rho_layout ! not sequential in any dir.
     type(remap_plan_3D_real64), pointer :: split_to_seqx1
     type(remap_plan_3D_real64), pointer :: seqx1_to_seqx2
     type(remap_plan_3D_real64), pointer :: seqx2_to_seqx3
     ! remaps for the electric field data
     type(remap_plan_3D_real64), pointer :: ex_x1_to_split
     type(remap_plan_3D_real64), pointer :: ey_x2_to_split
     type(remap_plan_3D_real64), pointer :: ez_x3_to_split
     type(remap_plan_6D_real64), pointer :: seqx1x2x3_to_seqx4x5x6
     type(remap_plan_6D_real64), pointer :: seqx4x5x6_to_seqx1x2x3
     ! interpolators and their pointers
     type(cubic_spline_1d_interpolator) :: interp_x1
     type(cubic_spline_1d_interpolator) :: interp_x2
     type(cubic_spline_1d_interpolator) :: interp_x3
     type(cubic_spline_1d_interpolator) :: interp_x4
     type(cubic_spline_1d_interpolator) :: interp_x5
     type(cubic_spline_1d_interpolator) :: interp_x6
   contains
     procedure, pass(sim) :: run => run_vp6d_cartesian
     procedure, pass(sim) :: init_from_file => init_vp6d_par_cart
  end type sll_simulation_6d_vlasov_poisson_cart

  interface sll_delete
     module procedure delete_vp6d_par_cart
  end interface sll_delete

contains

  subroutine init_vp6d_par_cart( sim, filename )
    intrinsic :: trim
    class(sll_simulation_6d_vlasov_poisson_cart), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
    sll_real64            :: dt
    sll_int32             :: number_iterations
    sll_int32             :: num_cells_x1
    sll_int32             :: num_cells_x2
    sll_int32             :: num_cells_x3
    sll_int32             :: num_cells_x4
    sll_int32             :: num_cells_x5
    sll_int32             :: num_cells_x6
    sll_int32, parameter  :: input_file = 99
    sll_int32             :: IO_stat

    namelist /sim_params/ dt, number_iterations
    namelist /grid_dims/ num_cells_x1, num_cells_x2, num_cells_x3
    namelist /grid_dims/ num_cells_x4, num_cells_x5, num_cells_x6
    ! Try to add here other parameters to initialize the mesh values like
    ! xmin, xmax and also for the distribution function initializer.
    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_vp6d_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, sim_params)
    read(input_file,grid_dims)
    close(input_file)

    sim%dt = dt
    sim%num_iterations = number_iterations
    ! In this particular simulation, since the system is periodic, the number
    ! of points is the same as the number of cells in all directions.
    sim%nc_x1 = num_cells_x1
    sim%nc_x2 = num_cells_x2
    sim%nc_x3 = num_cells_x3
    sim%nc_x4 = num_cells_x4
    sim%nc_x5 = num_cells_x5
    sim%nc_x6 = num_cells_x6
  end subroutine init_vp6d_par_cart


  subroutine run_vp6d_cartesian(sim)
    class(sll_simulation_6d_vlasov_poisson_cart), intent(inout) :: sim
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: loc_sz_x3
    sll_int32  :: loc_sz_x4
    sll_int32  :: loc_sz_x5
    sll_int32  :: loc_sz_x6
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_int32  :: m
    sll_int32  :: n
    sll_real64 :: vmin
    sll_real64 :: delta
    sll_real64 :: alpha
    sll_int32  :: itemp
    sll_int32  :: ierr
    sll_int32  :: itime
    sll_real64 :: ex
    sll_real64 :: ey
    sll_real64 :: ez
    sll_int32  :: itmp1
    sll_int32  :: itmp2

    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)

    ! allocate the layouts...
    sim%sequential_x1x2x3  => new_layout_6D( sll_world_collective )
    sim%sequential_x4x5x6  => new_layout_6D( sll_world_collective )
    sim%rho_seq_x1       => new_layout_3D( sll_world_collective )
    sim%rho_seq_x2       => new_layout_3D( sll_world_collective )
    sim%rho_seq_x3       => new_layout_3D( sll_world_collective )
    sim%split_rho_layout => new_layout_3D( sll_world_collective )

    ! layout for sequential operations in x4, x5 and x6. Make an even split for
    ! x1, x2 and x3, or as close as even if the power of 2 is not divisible by
    ! three. This should be set up at initialization.
    sim%power2 = int(log(real(sim%world_size))/log(2.0))
    call factorize_in_three_powers_of_two( sim%world_size, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%nproc_x3 )
    sim%nproc_x4 = 1
    sim%nproc_x5 = 1
    sim%nproc_x6 = 1
    
    call initialize_layout_with_distributed_6D_array( &
         sim%nc_x1, &
         sim%nc_x2, &
         sim%nc_x3, &
         sim%nc_x4, &
         sim%nc_x5, &
         sim%nc_x6, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%nproc_x3, &
         sim%nproc_x4, &
         sim%nproc_x5, &
         sim%nproc_x6, &
         sim%sequential_x4x5x6 )

    ! Use this information to initialize the layout that describes the result
    ! of computing rho (thus, a 3D layout). This layout is not useful to do 
    ! sequential operations in any of the three available directions. 
    call initialize_layout_with_distributed_3D_array( &
         sim%nc_x1, &
         sim%nc_x2, &
         sim%nc_x3, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%nproc_x3, &
         sim%split_rho_layout )
    
    ! We also initialize the other two layouts needed for both sequential 
    ! operations on x1, x2 and x3 in the 3D case.
    call factorize_in_two_powers_of_two( sim%world_size, itmp1, itmp2 )
    call initialize_layout_with_distributed_3D_array( &
         sim%nc_x1, &
         sim%nc_x2, &
         sim%nc_x3, &
         1, & 
         itmp1, &
         itmp2, &
         sim%rho_seq_x1 )
    
    call compute_local_sizes_3d(sim%rho_seq_x1, loc_sz_x1, loc_sz_x2, loc_sz_x3)
    SLL_ALLOCATE(sim%rho_x1(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi_x1(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)

    ! Since now each point has 3 electric field components, we can't use a 
    ! complex number array. So now we use three different arrays, ex, ey, 
    ! and ez and treat them separately. This is expensive and there should 
    ! be a better way to do this, but this will do for now. 
    SLL_ALLOCATE(sim%ex_x1(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)

    ! Sequential operations in x2:
    call initialize_layout_with_distributed_3D_array( &
         sim%nc_x1, &
         sim%nc_x2, &
         sim%nc_x3, &
         itmp1, &
         1, &
         itmp2, &
         sim%rho_seq_x2 )
    call compute_local_sizes_3d( sim%rho_seq_x2, loc_sz_x1,loc_sz_x2,loc_sz_x3 )
    SLL_ALLOCATE(sim%rho_x2(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi_x2(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)
    SLL_ALLOCATE(sim%ey_x2(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)
    
    ! Sequential operations in x3:
    call initialize_layout_with_distributed_3D_array( &
         sim%nc_x1, &
         sim%nc_x2, &
         sim%nc_x3, &
         itmp1, &
         itmp2, &
         1, &
         sim%rho_seq_x3 )
    call compute_local_sizes_3d( sim%rho_seq_x3, loc_sz_x1,loc_sz_x2,loc_sz_x3 )
    SLL_ALLOCATE(sim%rho_x3(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi_x3(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)
    SLL_ALLOCATE(sim%ez_x3(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)

    ! layout for sequential operations in x1, x2 and x3. This is basically 
    ! just the flipping of the values between x1 and x4, x2 and x5, and
    ! x3 and x6 on the previous layout.
    ! This is the wrong (dangerous) way to do this. We should not be changing
    ! the values of these variables, which should behave like parameters. 
    ! Choose a better name for each, initialize and don't change any further.

    ! switch x1 and x4:
    itemp = sim%nproc_x4
    sim%nproc_x4 = sim%nproc_x1
    sim%nproc_x1 = itemp
    ! switch x2 and x5
    itemp = sim%nproc_x5
    sim%nproc_x5 = sim%nproc_x2 
    sim%nproc_x2 = itemp
    ! switch x3 and x6
    itemp = sim%nproc_x6
    sim%nproc_x6 = sim%nproc_x3 
    sim%nproc_x3 = itemp
    
    call initialize_layout_with_distributed_6D_array( &
         sim%nc_x1, &
         sim%nc_x2, &
         sim%nc_x3, &
         sim%nc_x4, &
         sim%nc_x5, &
         sim%nc_x6, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%nproc_x3, &
         sim%nproc_x4, &
         sim%nproc_x5, &
         sim%nproc_x6, &
         sim%sequential_x1x2x3 )
    
    ! Allocate the array needed to store the local chunk of the distribution
    ! function data. First compute the local sizes. Since the remap operations
    ! are out-of-place, we will allocate multiple arrays, one for each layout.
    call compute_local_sizes_6d( sim%sequential_x1x2x3, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4, &
         loc_sz_x5, &
         loc_sz_x6 )
    SLL_ALLOCATE(sim%f_x1x2x3(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4,loc_sz_x5,loc_sz_x6),ierr)
    
    ! This layout is also useful to represent the charge density array. Since
    ! this is a result of a local reduction on x4, x5 and x6, the new layout is
    ! 3D but with the same dimensions of the process mesh in x1, x2 and x3.
    SLL_ALLOCATE(sim%rho_split(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)
    SLL_ALLOCATE(sim%ex_split(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)
    SLL_ALLOCATE(sim%ey_split(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)
    SLL_ALLOCATE(sim%ez_split(loc_sz_x1,loc_sz_x2,loc_sz_x3),ierr)
  
    call compute_local_sizes_6d( sim%sequential_x4x5x6, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4, &
         loc_sz_x5, &
         loc_sz_x6 )
    SLL_ALLOCATE(sim%f_x4x5x6(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4,loc_sz_x5,loc_sz_x6),ierr)
    
    ! These dimensions are also the ones needed for the array where we store
    ! the intermediate results of the charge density computation.
    SLL_ALLOCATE(sim%partial_reduction(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4,loc_sz_x5),ierr)
    
    ! Initialize the initial distribution function data. We do this with an
    ! initializer object which needs to be initialized itself! Note also that 
    ! the mesh is described in global terms. This should be fine as meshes are
    ! supposed to be read-only entities.
    ! This should be done elsewhere...
    sim%mesh6d => new_cartesian_6d_mesh( &
         sim%nc_x1, &
         sim%nc_x2, &
         sim%nc_x3, &
         sim%nc_x4, &
         sim%nc_x5, &
         sim%nc_x6, &
         0.0_f64, 1.0_f64, &
         0.0_f64, 1.0_f64, &
         0.0_f64, 1.0_f64, &
         -1.0_f64, 1.0_f64, &
         -1.0_f64, 1.0_f64, &
         -1.0_f64, 1.0_f64 )
    
    call load_test_6d_initializer( sim%init_6d, &
         NODE_CENTERED_FIELD, &
         sim%mesh6d, &
         0.1_f64, &
         sim%sequential_x4x5x6 )
    call sim%init_6d%f_of_6args(sim%f_x4x5x6)
    ! With the distribution function initialized in at least one configuration,
    ! we can proceed to carry out the computation of the electric potential.
    ! First we need to compute the charge density. Some thoughts:
    !
    ! The computation of rho is a reduction process that takes as input a 6d
    ! array and that should return a 3d array (or alternatively: a 6d array
    ! of size 1 in the reduced directions). For example, a df of dimensions
    ! np1 X np2 X np3 X np4 X np5 X np6 might effectively end up as an array 
    ! of dimensions
    ! np1 X np2 X np3 X np4 X np5 X 1  after a summation of all the values in 
    ! x6. After a
    ! second summation along x5, the dimensions would be:
    ! np1 X np2 X np3 X np4 X 1 X 1 and after a last summation:
    ! np1 X np2 X np3 X 1 X 1 X 1 
    ! One simple-minded but inefficient way to prepare the data for the triple
    ! reduction could be to have a layout in a process mesh 
    ! NP1 X NP2 X NP3 X 1 X 1 X 1
    ! where NP1xNP2xNP3 is the total number of processors available. The 
    ! problem 
    ! here is that the end result would still need a layout change to be fed
    ! into the Poisson solver...
    sim%rho_split(:,:,:) = 0.0

    call compute_charge_density_6d( &
         sim%mesh6d, &
         size(sim%f_x4x5x6,1), &
         size(sim%f_x4x5x6,2), &
         size(sim%f_x4x5x6,3), &
         sim%f_x4x5x6, &
         sim%partial_reduction, &
         sim%rho_split )
 
    ! Re-arrange rho_split in a way that permits sequential operations in x1, to
    ! feed to the Poisson solver.
    sim%split_to_seqx1 => &
         NEW_REMAP_PLAN(sim%split_rho_layout, sim%rho_seq_x1, sim%rho_split)
    call apply_remap_3D( sim%split_to_seqx1, sim%rho_split, sim%rho_x1 )
    
    ! We are in a position now to compute the electric potential.
    ! Initialize the poisson plan
    sim%poisson_plan => new_poisson_3d_periodic_plan_par( &
         sim%rho_seq_x1, &
         sim%nc_x1, &
         sim%nc_x2, &
         sim%nc_x3, &
         1.0_f64, &    ! parametrize with mesh values
         1.0_f64, &
         1.0_f64 )     ! parametrize with mesh values
    
    ! solve for the electric potential
    call solve_poisson_3d_periodic_par( &
         sim%poisson_plan, &
         sim%rho_x1, &
         sim%phi_x1)

    ! compute the values of the electric field. Since the electric field is 
    ! now stored in 3 different arrays, we need to do this in multiple steps.
    ! This is quite ugly and inefficient. Need to improve here...
    ! rho is configured for sequential operations in x1, thus we start by 
    ! computing the E_x component.
    ! The following call is inefficient and unnecessary. The local sizes for
    ! the arrays should be kept around as parameters basically and not on 
    ! variables whose content could be anything... This will have to do for now.
    call compute_local_sizes_3d( sim%rho_seq_x1,loc_sz_x1,loc_sz_x2,loc_sz_x3 )
    call compute_electric_field_x1_3d( &
         sim%phi_x1, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         sim%mesh6d%delta_x1, &
         sim%ex_x1 )
    ! note that we are 'recycling' the layouts used for the other arrays because
    ! they represent identical configurations.
    sim%ex_x1_to_split => &
         NEW_REMAP_PLAN( sim%rho_seq_x1, sim%split_rho_layout, sim%ex_x1)
    ! to reconfigure the potential to compute the electric field in X2.
    sim%seqx1_to_seqx2 => &
         NEW_REMAP_PLAN( sim%rho_seq_x1, sim%rho_seq_x2, sim%phi_x1 )
    ! it would be nice if the next call were executed by another thread...
    call apply_remap_3D( sim%ex_x1_to_split, sim%ex_x1, sim%ex_split )
    call apply_remap_3D( sim%seqx1_to_seqx2, sim%phi_x1, sim%phi_x2 )

    call compute_local_sizes_3d( sim%rho_seq_x2,loc_sz_x1,loc_sz_x2,loc_sz_x3 )
    call compute_electric_field_x2_3d( &
         sim%phi_x2, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         sim%mesh6d%delta_x2, &
         sim%ey_x2 )

    sim%ey_x2_to_split => &
         NEW_REMAP_PLAN( sim%rho_seq_x2, sim%split_rho_layout, sim%ey_x2)
    ! to reconfigure the potential to compute the electric field in X3.
    sim%seqx2_to_seqx3 => &
         NEW_REMAP_PLAN( sim%rho_seq_x2, sim%rho_seq_x3, sim%phi_x2 )
    ! it would be nice if the next call were executed by another thread...
    call apply_remap_3D( sim%ey_x2_to_split, sim%ey_x2, sim%ey_split )
    call apply_remap_3D( sim%seqx2_to_seqx3, sim%phi_x2, sim%phi_x3 )

    call compute_local_sizes_3d( sim%rho_seq_x3,loc_sz_x1,loc_sz_x2,loc_sz_x3 )
    call compute_electric_field_x3_3d( &
         sim%phi_x3, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         sim%mesh6d%delta_x3, &
         sim%ez_x3 )

    ! But now, to make the Ez field data configuration compatible with
    ! the sequential operations in x4x5x6...
    sim%ez_x3_to_split => &
         NEW_REMAP_PLAN( sim%rho_seq_x3, sim%split_rho_layout, sim%ez_x3)
    call apply_remap_3D( sim%ez_x3_to_split, sim%ez_x3, sim%ez_split )

    ! Now we proceed to reconfigure the data. This is very expensive. 
    ! There might be advantages to this approach if we avoid larger data 
    ! transfers like with an all-to-all transfer... however, we could end 
    ! up paying more if the simulation is latency-dominated.
    
    ! Proceed to carry out the advections. The following should go inside a
    ! subroutine...
    
    ! Start the interpolators... Watch out: the periodic case has equal number
    ! of cells than points. Is this properly handled by the interpolators??
    ! The interpolators need the number of points and always consider that
    ! num_cells = num_pts - 1. This is a possible source of confusion.

    call sim%interp_x1%initialize( &
         sim%nc_x1, &
         sim%mesh6d%x1_min, &
         sim%mesh6d%x1_max, &
         SLL_PERIODIC)

    call sim%interp_x2%initialize( &
         sim%nc_x2, &
         sim%mesh6d%x2_min, &
         sim%mesh6d%x2_max, &
         SLL_PERIODIC)

    call sim%interp_x3%initialize( &
         sim%nc_x3, &
         sim%mesh6d%x3_min, &
         sim%mesh6d%x3_max, &
         SLL_PERIODIC)

    call sim%interp_x4%initialize( &
         sim%nc_x4, &
         sim%mesh6d%x4_min, &
         sim%mesh6d%x4_max, &
         SLL_HERMITE)

    call sim%interp_x5%initialize( &
         sim%nc_x5, &
         sim%mesh6d%x5_min, &
         sim%mesh6d%x5_max, &
         SLL_HERMITE)

    call sim%interp_x6%initialize( &
         sim%nc_x6, &
         sim%mesh6d%x6_min, &
         sim%mesh6d%x6_max, &
         SLL_HERMITE)

    call compute_local_sizes_3d(sim%rho_seq_x1, loc_sz_x1, loc_sz_x2, loc_sz_x3)

    ! ------------------------------------------------------------------------
    !
    !                                MAIN LOOP
    !
    ! ------------------------------------------------------------------------

    sim%seqx4x5x6_to_seqx1x2x3 => &
       NEW_REMAP_PLAN(sim%sequential_x4x5x6,sim%sequential_x1x2x3,sim%f_x4x5x6)

    sim%seqx1x2x3_to_seqx4x5x6 => &
       NEW_REMAP_PLAN(sim%sequential_x1x2x3,sim%sequential_x4x5x6,sim%f_x1x2x3)


    do itime=1, sim%num_iterations
       if (sim%my_rank == 0) then
          print *, 'Iteration ', itime, ' of ', sim%num_iterations
       end if
       ! Carry out a 'dt/2' advection in the velocities.
       ! Start with vx...(x4)
       call compute_local_sizes_6d( sim%sequential_x4x5x6, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4, loc_sz_x5, loc_sz_x6 )
       do n=1, sim%mesh6d%num_cells6
          do m=1, sim%mesh6d%num_cells5
             do k=1, loc_sz_x3 
                do j=1,loc_sz_x2
                   do i=1,loc_sz_x1
                      ex    =  sim%ex_split(i,j,k)
                      alpha = -ex*0.5_f64*sim%dt
                      ! interpolate_array_disp() has an interface that must 
                      ! be changed.
                      sim%f_x4x5x6(i,j,k,:,m,n) = &
                           sim%interp_x4%interpolate_array_disp( &
                             sim%nc_x4, &
                             sim%f_x4x5x6(i,j,k,:,m,n), &
                             alpha )
                   end do
                end do
             end do
          end do
       end do

       ! Continue with vy...(x5)
       do n=1, sim%mesh6d%num_cells6
          do l=1, sim%mesh6d%num_cells4
             do k=1, loc_sz_x3
                do j=1,loc_sz_x2
                   do i=1,loc_sz_x1
                      ey    = sim%ey_split(i,j,k)
                      alpha = -ey*0.5_f64*sim%dt
                      ! interpolate_array_disp() has an interface that must 
                      ! be changed
                      sim%f_x4x5x6(i,j,k,l,:,n) = &
                           sim%interp_x5%interpolate_array_disp( &
                             sim%nc_x5, &
                             sim%f_x4x5x6(i,j,k,l,:,n), &
                             alpha )
                   end do
                end do
             end do
          end do
       end do

       ! Continue with vz...(x6)
       do m=1, sim%mesh6d%num_cells5
          do l=1, sim%mesh6d%num_cells4
             do k=1, loc_sz_x3 
                do j=1,loc_sz_x2
                   do i=1,loc_sz_x1
                      ez    =  sim%ez_split(i,j,k)
                      alpha = -ez*0.5_f64*sim%dt
                      ! interpolate_array_disp() has an interface that must 
                      ! be changed
                      sim%f_x4x5x6(i,j,k,l,m,:) = &
                           sim%interp_x6%interpolate_array_disp( &
                             sim%nc_x6, &
                             sim%f_x4x5x6(i,j,k,l,m,:), &
                             alpha )
                   end do
                end do
             end do
          end do
       end do

       ! Proceed to the advections in the spatial directions, 'x' and 'y'
       ! Reconfigure data for sequential operations in x, y, and z:
       call apply_remap_6D( sim%seqx4x5x6_to_seqx1x2x3, &
            sim%f_x4x5x6, sim%f_x1x2x3 )

       ! what are the new local limits on x4, x5 and x6? It is bothersome to 
       ! have to make these calls...
       call compute_local_sizes_6d( sim%sequential_x1x2x3, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4, loc_sz_x5, loc_sz_x6 )
       
       ! full time step advection in 'x' (x1)
       do n=1, loc_sz_x6
          do m=1, loc_sz_x5
             do l=1, loc_sz_x4 
                do k=1,sim%mesh6d%num_cells3
                   do j=1,sim%mesh6d%num_cells2
                      vmin = sim%mesh6d%x4_min
                      delta = sim%mesh6d%delta_x4
                      alpha = (vmin + (k-1)*delta)*sim%dt
                      sim%f_x1x2x3(:,j,k,l,m,n) = &
                           sim%interp_x1%interpolate_array_disp( &
                             sim%nc_x1, &
                             sim%f_x1x2x3(:,j,k,l,m,n), &
                             alpha )
                   end do
                end do
             end do
          end do
       end do

       ! full time step advection in 'y' (x2)
       do n=1, loc_sz_x6
          do m=1, loc_sz_x5
             do l=1,loc_sz_x4
                do k=1,sim%mesh6d%num_cells3
                   do i=1,sim%mesh6d%num_cells1
                      vmin = sim%mesh6d%x5_min
                      delta = sim%mesh6d%delta_x5
                      alpha = (vmin + (l-1)*delta)*sim%dt
                      sim%f_x1x2x3(i,:,k,l,m,n) = &
                           sim%interp_x2%interpolate_array_disp( &
                             sim%nc_x2, &
                             sim%f_x1x2x3(i,:,k,l,m,n), &
                             alpha )
                   end do
                end do
             end do
          end do
       end do

       ! Compute the fields:
       ! 1. Reconfigure data for sequential operations in x4, x5 and x6 in order
       !    to compute the charge density.
       ! 2. Compute charge density.
       ! 3. Reconfigure charge density to feed to Poisson solver
       call apply_remap_6D( sim%seqx1x2x3_to_seqx4x5x6, &
            sim%f_x1x2x3, sim%f_x4x5x6 )

       call compute_charge_density_6d( &
            sim%mesh6d, &
            size(sim%f_x4x5x6,1), &
            size(sim%f_x4x5x6,2), &
            size(sim%f_x4x5x6,3), &
            sim%f_x4x5x6, &
            sim%partial_reduction, &
            sim%rho_split )

       ! 3d charge density is 'fully split', no sequential operations can be
       ! fully done. Thus a remap is needed.
       call apply_remap_3D( sim%split_to_seqx1, sim%rho_split, sim%rho_x1 )

       ! Compute the electric potential.
       call solve_poisson_3d_periodic_par( &
            sim%poisson_plan, &
            sim%rho_x1, &
            sim%phi_x1)

       ! compute the values of the electric field. rho is configured for 
       ! sequential operations in x1, thus we start by computing the E_x 
       ! component.
       ! The following call is inefficient and unnecessary. The local sizes for
       ! the arrays should be kept around as parameters basically and not on 
       ! variables whose content could be anything... This will have to do for 
       ! now.
       call compute_local_sizes_3d(sim%rho_seq_x1,loc_sz_x1,loc_sz_x2,loc_sz_x3)
       call compute_electric_field_x1_3d( &
            sim%phi_x1, &
            loc_sz_x1, &
            loc_sz_x2, &
            loc_sz_x3, &
            sim%mesh6d%delta_x1, &
            sim%ex_x1 )
       call apply_remap_3D( sim%ex_x1_to_split, sim%ex_x1, sim%ex_split )

       call compute_local_sizes_3d(sim%rho_seq_x2,loc_sz_x1,loc_sz_x2,loc_sz_x3)
       call compute_electric_field_x2_3d( &
            sim%phi_x2, &
            loc_sz_x1, &
            loc_sz_x2, &
            loc_sz_x3, &
            sim%mesh6d%delta_x2, &
            sim%ey_x2 )

       ! it would be nice if the next call were executed by another thread...
       call apply_remap_3D( sim%ey_x2_to_split, sim%ey_x2, sim%ey_split )

       call compute_local_sizes_3d(sim%rho_seq_x3,loc_sz_x1,loc_sz_x2,loc_sz_x3)
       call compute_electric_field_x3_3d( &
            sim%phi_x3, &
            loc_sz_x1, &
            loc_sz_x2, &
            loc_sz_x3, &
            sim%mesh6d%delta_x3, &
            sim%ez_x3 )

       ! But now, to make the Ez field data configuration compatible with
       ! the sequential operations in x4x5x6...
       call apply_remap_3D( sim%ez_x3_to_split, sim%ez_x3, sim%ez_split )

       ! ...and another half time step advection in the velocities.
       ! Start with vx...(x4)
       call compute_local_sizes_6d( sim%sequential_x4x5x6, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4, loc_sz_x5, loc_sz_x6 )
       
       do n=1, sim%mesh6d%num_cells6
          do m=1, sim%mesh6d%num_cells5
             do k=1, loc_sz_x2 
                do j=1,loc_sz_x2
                   do i=1,loc_sz_x1
                      ex    = sim%ex_split(i,j,k)
                      alpha = -ex*0.5_f64*sim%dt
                      ! interpolate_array_disp() has an interface that must 
                      ! be changed
                      sim%f_x4x5x6(i,j,k,:,m,n) = &
                           sim%interp_x4%interpolate_array_disp( &
                             sim%nc_x4, &
                             sim%f_x4x5x6(i,j,k,:,m,n), &
                             alpha )
                   end do
                end do
             end do
          end do
       end do

       ! Continue with vy...(x5)
       do n=1, sim%mesh6d%num_cells6
          do l=1, sim%mesh6d%num_cells4
             do k=1, loc_sz_x2 
                do j=1,loc_sz_x2
                   do i=1,loc_sz_x1
                      ey    = sim%ey_split(i,j,k)
                      alpha = -ey*0.5_f64*sim%dt
                      ! interpolate_array_disp() has an interface that must 
                      ! be changed
                      sim%f_x4x5x6(i,j,k,l,:,n) = &
                           sim%interp_x5%interpolate_array_disp( &
                             sim%nc_x5, &
                             sim%f_x4x5x6(i,j,k,l,:,n), &
                             alpha )
                   end do
                end do
             end do
          end do
       end do

       ! Continue with vz...(x6)
       do m=1, sim%mesh6d%num_cells5
          do l=1, sim%mesh6d%num_cells4
             do k=1, loc_sz_x2 
                do j=1,loc_sz_x2
                   do i=1,loc_sz_x1
                      ez    =  sim%ez_split(i,j,k)
                      alpha = -ez*0.5_f64*sim%dt
                      ! interpolate_array_disp() has an interface that must 
                      ! be changed
                      sim%f_x4x5x6(i,j,k,l,m,:) = &
                           sim%interp_x6%interpolate_array_disp( &
                             sim%nc_x6, &
                             sim%f_x4x5x6(i,j,k,l,m,:), &
                             alpha )
                   end do
                end do
             end do
          end do
       end do
       ! The following function needs to be modified for the 3D case
  !     call plot_fields(itime, sim)

    end do ! main loop

    ! Test graphical output
    call test_write(sim)

  end subroutine run_vp6d_cartesian

  function divisible_by_two( num )
    logical :: divisible_by_two
    sll_int32, intent(in) :: num
    if(modulo(num,2) == 0) then
       divisible_by_two = .true.
    else
       divisible_by_two = .false.
    end if
  end function divisible_by_two

  function divisible_by_three( num )
    logical :: divisible_by_three
    sll_int32, intent(in) :: num
    if(modulo(num,3) == 0) then
       divisible_by_three = .true.
    else
       divisible_by_three = .false.
    end if
  end function divisible_by_three

  ! given a number that is a power of two, we decompose it in 3 factors 
  ! intended to be as close to each other as possible, while still keeping
  ! them factors of two as well.
  subroutine factorize_in_three_powers_of_two( num_procs, f1, f2, f3 )
    sll_int32, intent(in)  :: num_procs
    sll_int32, intent(out) :: f1
    sll_int32, intent(out) :: f2
    sll_int32, intent(out) :: f3
    sll_int32  :: exponent
    sll_int32  :: tmpi

    SLL_ASSERT( is_power_of_two(int(num_procs,f64)) )
    exponent = int(log(real(num_procs))/log(2.0))
    if( (exponent > 0) .and. divisible_by_three(exponent) ) then
       tmpi = exponent/3
       f1 = 2**tmpi
       f2 = 2**tmpi
       f3 = 2**tmpi
    else
       if( exponent == 0 ) then
          f1   = 1
          f2   = 1
          f3   = 1
       else if( exponent == 1 ) then
          f1   = 2
          f2   = 1
          f3   = 1
       else if( divisible_by_three(exponent-1) ) then
          tmpi = (exponent-1)/3
          f1   = 2**tmpi
          f2   = 2**tmpi
          f3   = 2**(tmpi+1)
       else if( divisible_by_three(exponent+1) ) then
          tmpi = (exponent+1)/3
          f1   = 2**tmpi
          f2   = 2**tmpi
          f3   = 2**(tmpi-1)
       end if
    end if
  end subroutine factorize_in_three_powers_of_two

  subroutine factorize_in_two_powers_of_two( num_procs, f1, f2 )
    sll_int32, intent(in)  :: num_procs
    sll_int32, intent(out) :: f1
    sll_int32, intent(out) :: f2
    sll_int32  :: exponent
    sll_int32  :: tmpi

    SLL_ASSERT( is_power_of_two(int(num_procs,f64)) )
    exponent = int(log(real(num_procs))/log(2.0))
    if( (exponent > 0) .and. divisible_by_two(exponent) ) then
       tmpi = exponent/2
       f1 = 2**tmpi
       f2 = 2**tmpi
    else
       if( exponent == 0 ) then
          f1   = 1
          f2   = 1
       else 
          tmpi = (exponent-1)/2
          f1   = 2**((exponent-1)/2)
          f2   = 2**((exponent+1)/2)
       end if
    end if
  end subroutine factorize_in_two_powers_of_two


  subroutine delete_vp6d_par_cart( sim )
    class(sll_simulation_6d_vlasov_poisson_cart) :: sim
    sll_int32 :: ierr
    SLL_DEALLOCATE( sim%f_x1x2x3, ierr )
    SLL_DEALLOCATE( sim%f_x4x5x6, ierr )
    SLL_DEALLOCATE_ARRAY( sim%partial_reduction, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x3, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_split, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_x3, ierr )
    SLL_DEALLOCATE_ARRAY( sim%ex_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%ex_split, ierr )
    SLL_DEALLOCATE_ARRAY( sim%ey_x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%ey_split, ierr )
    SLL_DEALLOCATE_ARRAY( sim%ez_x3, ierr )
    SLL_DEALLOCATE_ARRAY( sim%ez_split, ierr )
    call sll_delete( sim%sequential_x1x2x3 )
    call sll_delete( sim%sequential_x4x5x6 )
    call sll_delete( sim%rho_seq_x1 )
    call sll_delete( sim%rho_seq_x2 )
    call sll_delete( sim%rho_seq_x3 )
    call sll_delete( sim%split_rho_layout )
    call sll_delete( sim%split_to_seqx1 )
    call sll_delete( sim%seqx1_to_seqx2 )
    call sll_delete( sim%seqx2_to_seqx3 )
    call sll_delete( sim%ex_x1_to_split )
    call sll_delete( sim%ey_x2_to_split )
    call sll_delete( sim%ez_x3_to_split )
    call sll_delete( sim%seqx1x2x3_to_seqx4x5x6 )
    call sll_delete( sim%seqx4x5x6_to_seqx1x2x3 )
    call delete( sim%interp_x1 )
    call delete( sim%interp_x2 )
    call delete( sim%interp_x3 )
    call delete( sim%interp_x4 )
    call delete( sim%interp_x5 )
    call delete( sim%interp_x6 )
  end subroutine delete_vp6d_par_cart

  ! we put the reduction functions here for now, since we are only using
  ! simple data for the distribution function. This should go elsewhere.
  ! THIS SUBROUTINE IS JUST A PLACEHOLDER, IT IS NUMERICALLY INCORRECT.
  ! Change it later by something that uses some acceptable integrator in
  ! 1D.
  ! Design issues with this subroutine:
  ! 1. The distribution function needs to be preserved, thus this is an
  !    out-of-place operation.
  ! 2. There is probably a cleverer way to do this, but if the reduction
  !    happens in two steps a. reduction in x4 and b. reduction in x3, we
  !    need an array to store the intermediate result (after reducing in
  !    x4). This array should come as an argument.
  subroutine compute_charge_density_6d( &
    mesh, &
    numpts1, &
    numpts2, &
    numpts3, &
    f, &
    partial, &
    rho )

    type(simple_cartesian_6d_mesh), pointer        :: mesh
    sll_real64, intent(in), dimension(:,:,:,:,:,:) :: f   ! local distr. func
    sll_real64, intent(inout), dimension(:,:,:,:,:):: partial !intermediate res.
    sll_real64, intent(inout), dimension(:,:,:)    :: rho     ! local rho
    ! local sizes in the split directions have to be given by caller.
    sll_int32, intent(in)                       :: numpts1
    sll_int32, intent(in)                       :: numpts2
    sll_int32, intent(in)                       :: numpts3
    sll_real64                                  :: delta4
    sll_real64                                  :: delta5
    sll_real64                                  :: delta6
    sll_int32                                   :: numpts4
    sll_int32                                   :: numpts5
    sll_int32                                   :: numpts6
    sll_int32 :: i, j, k, l, m, n
    
    delta4   = mesh%delta_x4
    delta5   = mesh%delta_x5
    delta6   = mesh%delta_x6
    partial(:,:,:,:,:) = 0.0 
    numpts4 = mesh%num_cells4
    numpts5 = mesh%num_cells5
    numpts6 = mesh%num_cells6

    ! 'partial' must have been initialized to zero.
    ! reduction in x6:
    do m=1, numpts5
       do l=1,numpts4   
          do k=1,numpts3
             do j=1,numpts2
                do i=1,numpts1
                   ! This summation happens on a super-long stride... slow stuff
                   ! This loop should be substituted by a proper integration
                   ! function that we could use in the other directions as well.
                   do n=1,numpts4
                      partial(i,j,k,l,m) = partial(i,j,k,l,m) + &
                           f(i,j,k,l,m,n)*delta6
                   end do
                end do
             end do
          end do
       end do
    end do
    
    ! reduction on x5.
    do l=1,numpts4
       do k=1,numpts3
          do j=1,numpts2
             do i=1,numpts1
                do m=2,numpts5
                   ! This summation happens on a very-long stride... slow stuff
                   ! This loop should be substituted by a proper integration
                   ! function that we could use in the other directions as well.
                   ! See above reduction function for same problem.
                   partial(i,j,k,l,1) = partial(i,j,k,l,1) + &
                        partial(i,j,k,l,m)*delta5
                end do
             end do
          end do
       end do
    end do

    ! Final reduction on x4. Note that rho is not initialized
    ! to zero since it may already have the partial charge accumulation from
    ! other species.
    do k=1,numpts3
       do j=1,numpts2
          do i=1,numpts1
             do l=1,numpts4
                rho(i,j,k) = rho(i,j,k) + partial(i,j,k,l,1)*delta4
             end do
          end do
       end do
    end do
  end subroutine compute_charge_density_6d
  
  ! Temporary utility to compute the values of the electric field given 
  ! a pointer to an array of double precision values. It uses 
  ! forward/backward differencing schemes for the end points and a 
  ! centered one for the interior points. 
  subroutine compute_electric_field_on_line( &
       phi, &
       num_pts, &
       delta, &
       efield )
    
    sll_real64, dimension(:), intent(in) :: phi
    sll_int32                            :: num_pts
    sll_real64, intent(in)               :: delta
    sll_real64, dimension(:), intent(out):: efield
    sll_int32                            :: i
    sll_real64                           :: r_delta  ! reciprocal
    
    ! FIXME: check arrays sizes
    
    r_delta = 1.0_f64/delta
    
    ! Do first point:
    efield(1) = r_delta*(-1.5_f64*phi(1) + 2.0_f64*phi(2) - 0.5_f64*phi(3))
    
    ! Do the internal values:
    do i=2,num_pts-1
       efield(i) = r_delta*(phi(i+1) - phi(i-1))
    end do
    
    ! Do last point:
    efield(num_pts) = r_delta*( 0.5_f64*phi(num_pts-2) - &
         2.0_f64*phi(num_pts-1) + &
         1.5_f64*phi(num_pts) )
  end subroutine compute_electric_field_on_line
  
  subroutine compute_electric_field_x1_3d( &
    phi_x1, &
    num_pts_x1, &
    num_pts_x2, &
    num_pts_x3, &
    delta_x1, &
    efield_x1 )

    sll_real64, dimension(:,:,:), intent(in)  :: phi_x1
    sll_int32, intent(in)                     :: num_pts_x1
    sll_int32, intent(in)                     :: num_pts_x2
    sll_int32, intent(in)                     :: num_pts_x3
    sll_real64, intent(in)                    :: delta_x1
    sll_real64, dimension(:,:,:), intent(out) :: efield_x1
    sll_int32                                 :: i
    sll_int32                                 :: j
    sll_int32                                 :: k
    sll_real64                                :: r_delta
    sll_real64                                :: ex
    ! FIXME: arg checking
    
    r_delta = 1.0_f64/delta_x1
    
    ! Compute the electric field values on the left and right edges.
    do k=1,num_pts_x3
       do j=1,num_pts_x2
          ! i=1 plane:
          ex = r_delta*(-1.5_f64*phi_x1(1,j,k) + &
                         2.0_f64*phi_x1(2,j,k) - &
                         0.5_f64*phi_x1(3,j,k) )
          efield_x1(1,j,k) = ex  
          ! i=num_pts_x1 plane:
          ex = r_delta*(0.5_f64*phi_x1(num_pts_x1-2,j,k) - &
                        2.0_f64*phi_x1(num_pts_x1-1,j,k) + &
                        1.5_f64*phi_x1(num_pts_x1  ,j,k) )
          efield_x1(num_pts_x1,j,k) = ex
       end do
    end do
    
    ! Electric field in interior points
    do k=1,num_pts_x3
       do j=1,num_pts_x2
          do i=2, num_pts_x1-1
             ex = r_delta*0.5_f64*(phi_x1(i+1,j,k) - phi_x1(i-1,j,k))
             efield_x1(i,j,k) = ex
          end do
       end do
    end do

  end subroutine compute_electric_field_x1_3d

  
  ! This function only sets the Ey component of the electric field.
  subroutine compute_electric_field_x2_3d( &
       phi_x2, &
       num_pts_x1, &
       num_pts_x2, &
       num_pts_x3, &
       delta_x2, &
       efield_x2 )
    
    sll_real64, dimension(:,:,:), intent(in)  :: phi_x2
    sll_int32, intent(in)                     :: num_pts_x1
    sll_int32, intent(in)                     :: num_pts_x2
    sll_int32, intent(in)                     :: num_pts_x3
    sll_real64, intent(in)                    :: delta_x2
    sll_real64, dimension(:,:,:), intent(out) :: efield_x2
    sll_int32                                 :: i
    sll_int32                                 :: j
    sll_int32                                 :: k
    sll_real64                                :: r_delta
    sll_real64                                :: ey

    ! FIXME: arg checking
    
    r_delta = 1.0_f64/delta_x2
    
    ! Compute the electric field values on the bottom and top edges.
    do k=1,num_pts_x3
       do i=1,num_pts_x1
          ! bottom:
          ey = r_delta*(-1.5_f64*phi_x2(i,1,k) + 2.0_f64*phi_x2(i,2,k) - &
                         0.5_f64*phi_x2(i,3,k))
          efield_x2(i,1,k) = ey
          ! top:
          ey = r_delta*(0.5_f64*phi_x2(i,num_pts_x2-2,k) - &
                        2.0_f64*phi_x2(i,num_pts_x2-1,k) + &
                        1.5_f64*phi_x2(i,num_pts_x2  ,k))
          efield_x2(i,num_pts_x2,k) = ey 
       end do
    end do
    
    ! Electric field in interior points
    do k=1,num_pts_x3
       do j=2,num_pts_x2-1
          do i=1, num_pts_x1
             ey = r_delta*0.5_f64*(phi_x2(i,j+1,k) - phi_x2(i,j-1,k))
             efield_x2(i,j,k) = ey
          end do
       end do
    end do
  end subroutine compute_electric_field_x2_3d

  ! This function only sets the Ey component of the electric field.
  subroutine compute_electric_field_x3_3d( &
    phi_x3, &
    num_pts_x1, &
    num_pts_x2, &
    num_pts_x3, &
    delta_x3, &
    efield_x3 )
    
    sll_real64, dimension(:,:,:), intent(in)  :: phi_x3
    sll_int32, intent(in)                     :: num_pts_x1
    sll_int32, intent(in)                     :: num_pts_x2
    sll_int32, intent(in)                     :: num_pts_x3
    sll_real64, intent(in)                    :: delta_x3
    sll_real64, dimension(:,:,:), intent(out) :: efield_x3
    sll_int32                                 :: i
    sll_int32                                 :: j
    sll_int32                                 :: k
    sll_real64                                :: r_delta
    sll_real64                                :: ez

    ! FIXME: arg checking
    
    r_delta = 1.0_f64/delta_x3
    
    ! Compute the electric field values on the end faces.
    do j=1,num_pts_x2
       do i=1,num_pts_x1
          ez = r_delta*(-1.5_f64*phi_x3(i,j,1) + 2.0_f64*phi_x3(i,j,2) - &
                         0.5_f64*phi_x3(i,j,3))
          efield_x3(i,j,1) = ez
          ! top:
          ez = r_delta*(0.5_f64*phi_x3(i,j,num_pts_x3-2) - &
                        2.0_f64*phi_x3(i,j,num_pts_x3-1) + &
                        1.5_f64*phi_x3(i,j,num_pts_x3  ))
          efield_x3(i,j,num_pts_x3) = ez
       end do
    end do
    
    ! Electric field in interior points
    do k=2,num_pts_x3-1
       do j=1,num_pts_x2
          do i=1, num_pts_x1
             ez = r_delta*0.5_f64*(phi_x3(i,j,k+1) - phi_x3(i,j,k-1))
             efield_x3(i,j,k) = ez
          end do
       end do
    end do
  end subroutine compute_electric_field_x3_3d

  
  ! Tentative advection routines.
  subroutine advection_x_1d( dt, vmin, delta_v, num_pts, f_line, f_interp )
    sll_real64, intent(in)                      :: dt
    sll_real64, intent(in)                      :: vmin
    sll_real64, intent(in)                      :: delta_v
    sll_int32, intent(in)                       :: num_pts
    sll_real64, dimension(:), intent(inout)     :: f_line
    class(sll_interpolator_1d_base)             :: f_interp
    sll_int32  :: i
    sll_real64 :: displacement

    do i=1, num_pts
       displacement = (vmin + real(i-1,f64)*delta_v)*dt
       ! remember that the function interpolate_array_disp() has the wrong
       ! interface since it should be a subroutine, not a function.
       f_line = f_interp%interpolate_array_disp(num_pts, f_line, displacement)
    end do
  end subroutine advection_x_1d

  subroutine advection_v_1d(dt, efield, num_pts, f_line, f_interp)
    sll_real64, intent(in)                   :: dt
    sll_real64, dimension(:), intent(in)     :: efield
    sll_int32, intent(in)                    :: num_pts
    sll_real64, dimension(:), intent(inout)  :: f_line
    class(sll_interpolator_1d_base), pointer :: f_interp
    sll_int32                                :: i
    sll_real64                               :: displacement

    do i=1, num_pts
       ! Why is the negative sign there?
       displacement = -efield(i)*0.5_f64*dt
       f_line = f_interp%interpolate_array_disp(num_pts, f_line, displacement)
    end do
  end subroutine advection_v_1d

  subroutine test_write(sim)
    use sll_xdmf_parallel
    class(sll_simulation_6d_vlasov_poisson_cart), intent(in) :: sim
    type(layout_3D), pointer :: my_layout
    integer(HSIZE_T), dimension(3)  :: array_dims 
    integer(HSSIZE_T), dimension(3) :: offset 
    sll_int32,  dimension(3) :: global_indices
    sll_real64, dimension(:,:,:), allocatable :: x1
    sll_real64, dimension(:,:,:), allocatable :: x2
    sll_real64, dimension(:,:,:), allocatable :: x3
    sll_int32  :: file_id
    sll_int32  :: my_rank
    sll_int32  :: world_size
    sll_int32  :: error
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: delta_x3
    sll_int32  :: local_nx1
    sll_int32  :: local_nx2
    sll_int32  :: local_nx3
    sll_real64 :: x1_min
    sll_real64 :: x2_min
    sll_real64 :: x3_min
    sll_int32  :: i, j, k, gi, gj, gk

    array_dims(1) =  sim%nc_x1
    array_dims(2) =  sim%nc_x2
    array_dims(3) =  sim%nc_x3

    world_size    =  sll_get_collective_size(sll_world_collective)
    my_rank       =  sll_get_collective_rank(sll_world_collective)

    my_layout     => sim%rho_seq_x1
    offset(1)     =  get_layout_i_min( my_layout, my_rank ) - 1
    offset(2)     =  get_layout_j_min( my_layout, my_rank ) - 1
    offset(3)     =  get_layout_k_min( my_layout, my_rank ) - 1

    call compute_local_sizes_3d(my_layout, local_nx1, local_nx2, local_nx3)
    SLL_ALLOCATE(x1(local_nx1,local_nx2,local_nx3),error)
    SLL_ALLOCATE(x2(local_nx1,local_nx2,local_nx3),error)
    SLL_ALLOCATE(x3(local_nx1,local_nx2,local_nx3),error)

    x1_min   = sim%mesh6d%x1_min
    x2_min   = sim%mesh6d%x2_min
    x3_min   = sim%mesh6d%x3_min

    delta_x1 = sim%mesh6d%delta_x1
    delta_x2 = sim%mesh6d%delta_x2
    delta_x3 = sim%mesh6d%delta_x3

    do k = 1, local_nx3
       do j = 1, local_nx2
          do i = 1, local_nx1
             global_indices = local_to_global_3D( my_layout, (/i, j, k/) )
             gi = global_indices(1)
             gj = global_indices(2)
             gk = global_indices(3)
             x1(i,j,k) = x1_min + (gi-1._f64)*delta_x1
             x2(i,j,k) = x2_min + (gj-1._f64)*delta_x2
             x3(i,j,k) = x3_min + (gk-1._f64)*delta_x3
          end do
       end do
    end do
       
    call sll_xdmf_open(my_rank,"test.xmf","grid",               &
                       sim%nc_x1,sim%nc_x2,sim%nc_x3,   &
                       file_id,error)
    call sll_xdmf_write_array("grid",array_dims,offset,x1,'x1',error)
    call sll_xdmf_write_array("grid",array_dims,offset,x2,'x2',error)
    call sll_xdmf_write_array("grid",array_dims,offset,x3,'x3',error)
    call sll_xdmf_write_array("grid",array_dims,     &
                              offset,sim%rho_x1,"rho",error,file_id,"Node")
    call sll_xdmf_close(file_id,error)

    deallocate(x1)
    deallocate(x2)
    deallocate(x3)

  end subroutine test_write


  subroutine plot_fields(itime, sim)
    use sll_collective
    use hdf5
    use sll_hdf5_io_parallel
    use sll_xml_io
    sll_int32, intent(in) :: itime
    character(len=4)      :: ctime
    sll_int32             :: i_layout
    character(len=1)      :: c_layout
    class(sll_simulation_6d_vlasov_poisson_cart), intent(in) :: sim
    type(layout_3D), pointer :: my_layout
    character(len=7),  parameter :: hdf_file = "data.h5"  ! File name
    sll_real64 :: tcpu1, tcpu2
    sll_int32  :: my_rank
    sll_int32  :: world_size
    sll_int32  :: local_nx1
    sll_int32  :: local_nx2
    sll_int32  :: local_nx3
    sll_int32  :: global_nx1
    sll_int32  :: global_nx2
    sll_int32  :: global_nx3
    sll_int32  :: error
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: gi
    sll_int32  :: gj
    sll_int32  :: gk
    sll_int32,  dimension(3) :: global_indices
    sll_real64, dimension(:,:,:), allocatable :: x1
    sll_real64, dimension(:,:,:), allocatable :: x2
    sll_real64, dimension(:,:,:), allocatable :: x3
    sll_real64 :: x1_min
    sll_real64 :: x1_max
    sll_real64 :: x2_min
    sll_real64 :: x2_max
    sll_real64 :: x3_min
    sll_real64 :: x3_max
    sll_real64 :: x4_min
    sll_real64 :: x4_max
    sll_real64 :: x5_min
    sll_real64 :: x5_max
    sll_real64 :: x6_min
    sll_real64 :: x6_max
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: delta_x3
    sll_real64 :: delta_x4 
    sll_real64 :: delta_x5 
    sll_real64 :: delta_x6 

    integer(HID_T)                  :: hdf_file_id
    sll_int32                       :: xml_file_id
    integer(HSIZE_T), dimension(3)  :: array_dims 
    integer(HSSIZE_T), dimension(3) :: offset 

    array_dims(1) = sim%nc_x1
    array_dims(2) = sim%nc_x2
    array_dims(3) = sim%nc_x3
    ! The following inquiry should be on a specific layout, not on world...
    world_size    = sll_get_collective_size(sll_world_collective)
    my_rank       = sll_get_collective_rank(sll_world_collective)

    tcpu1 = MPI_WTIME()

    do i_layout = 1, 2
       if (i_layout == 1) then
          my_layout => sim%rho_seq_x1
       else
          my_layout => sim%rho_seq_x2
       end if

       call compute_local_sizes_3d(my_layout, local_nx1, local_nx2, local_nx3)
    
       offset(1) =  get_layout_i_min( my_layout, my_rank ) - 1
       offset(2) =  get_layout_j_min( my_layout, my_rank ) - 1
       offset(3) =  get_layout_k_min( my_layout, my_rank ) - 1

       if (itime == 1) then
          SLL_ALLOCATE(x1(local_nx1,local_nx2,local_nx3),error)
          SLL_ALLOCATE(x2(local_nx1,local_nx2,local_nx3),error)
          SLL_ALLOCATE(x3(local_nx1,local_nx2,local_nx3),error)
          x1_min = sim%mesh6d%x1_min
          x1_max = sim%mesh6d%x1_max
          x2_min = sim%mesh6d%x2_min
          x2_max = sim%mesh6d%x2_max
          x3_min = sim%mesh6d%x3_min
          x3_max = sim%mesh6d%x3_max
          x4_min = sim%mesh6d%x4_min
          x4_max = sim%mesh6d%x4_max
          x5_min = sim%mesh6d%x5_min
          x5_max = sim%mesh6d%x5_max
          x6_min = sim%mesh6d%x6_min
          x6_max = sim%mesh6d%x6_max
   
          delta_x1 = sim%mesh6d%delta_x1
          delta_x2 = sim%mesh6d%delta_x2
          delta_x3 = sim%mesh6d%delta_x3
          delta_x4 = sim%mesh6d%delta_x4
          delta_x5 = sim%mesh6d%delta_x5
          delta_x6 = sim%mesh6d%delta_x6
   
          do k = 1, local_nx3
             do j = 1, local_nx2
                do i = 1, local_nx1
                   global_indices = local_to_global_3D( my_layout, (/i, j, k/) )
                   gi = global_indices(1)
                   gj = global_indices(2)
                   gk = global_indices(3)
                   x1(i,j,k) = x1_min + (gi-1._f64)*delta_x1
                   x2(i,j,k) = x2_min + (gj-1._f64)*delta_x2
                   x3(i,j,k) = x3_min + (gk-1._f64)*delta_x3
                end do
             end do
          end do
       
          call sll_hdf5_file_create("mesh_x"//c_layout//"_seq.h5",hdf_file_id,error)
          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,x1,"x1",error)
          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,x2,"x2",error)
          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,x3,"x3",error)
          call sll_hdf5_file_close(hdf_file_id,error)

          deallocate(x1)
          deallocate(x2)
          deallocate(x3)
       end if

       call int2string(itime, ctime)
       c_layout = char(i_layout+48)

       call sll_hdf5_file_create("fields_x"//c_layout//"-"//ctime//".h5", &
                                 hdf_file_id,error)

       if (i_layout == 1) then
          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,sim%rho_x1, &
                                    "rho_x"//c_layout,error)
          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,sim%phi_x1, &
                                    "phi_x"//c_layout,error)
       else
          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,sim%rho_x2, &
                                    "rho_x"//c_layout,error)
          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,sim%phi_x2, &
                                    "phi_x"//c_layout,error)
       end if

       call sll_hdf5_file_close(hdf_file_id,error)
   
       if (my_rank == 0) then
          
          !Conversion int64 -> int32
          global_nx1 = transfer(array_dims(1),global_nx1)
          global_nx2 = transfer(array_dims(2),global_nx2)
          global_nx3 = transfer(array_dims(3),global_nx3)
       ! Pierre: further changes are needed below...
          call sll_xml_file_create("fields_x"//c_layout//"-"//ctime//".xmf", &
                                   xml_file_id,error)
          call sll_xml_grid_geometry(xml_file_id,          &
                                  "mesh_x"//c_layout//"_seq.h5",global_nx1, &
                                  "mesh_x"//c_layout//"_seq.h5",global_nx2, &
                                  "x1", "x2" )
          call sll_xml_field(xml_file_id,'rho_x'//c_layout,  &
                             "fields_x"//c_layout//"-"//ctime//".h5:/rho_x"//c_layout, &
                             global_nx1, global_nx2,'HDF','Node')
          call sll_xml_field(xml_file_id,'phi_x'//c_layout,  &
                             "fields_x"//c_layout//"-"//ctime//".h5:/phi_x"//c_layout, &
                          global_nx1, global_nx2,'HDF','Node')
          call sll_xml_file_close(xml_file_id,error)

       end if



    end do

    tcpu2 = MPI_WTIME()
    !if (my_rank == 0) &
    !   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*world_size
    
  end subroutine plot_fields
  
end module sll_simulation_6d_vlasov_poisson_cartesian
