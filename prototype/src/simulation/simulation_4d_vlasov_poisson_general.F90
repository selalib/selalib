module sll_simulation_4d_vlasov_poisson_general

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"

  use sll_collective
  use sll_remapper
  use sll_constants
  use sll_cubic_spline_interpolator_2d
  use sll_poisson_2d_periodic_cartesian_par
  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  use sll_logical_meshes
  use sll_parallel_array_initializer_module
  use sll_coordinate_transformation_2d_base_module
  use sll_gnuplot_parallel
  implicit none

  type, extends(sll_simulation_base_class) :: sll_simulation_4d_vp_general
  
     ! Parallel environment parameters
     sll_int32  :: world_size
     sll_int32  :: my_rank
     sll_int32  :: power2 ! 2^power2 = number of processes available
     ! Processor mesh sizes
     sll_int32  :: nproc_x1
     sll_int32  :: nproc_x2
     sll_int32  :: nproc_x3
     sll_int32  :: nproc_x4 
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     ! Mesh parameters
     sll_int32  :: nc_x1
     sll_int32  :: nc_x2
     sll_int32  :: nc_x3
     sll_int32  :: nc_x4
     ! the logical meshes are split in two one for space, one for velocity
     type(sll_logical_mesh_2d), pointer    :: mesh2d_x
     type(sll_logical_mesh_2d), pointer    :: mesh2d_v
     ! This simulation only applies a coordinate transformation to the spatial
     ! coordinates.
     class(sll_coordinate_transformation_2d_base), pointer :: transfx
     type(poisson_2d_periodic_plan_cartesian_par), pointer :: poisson_plan

     ! distribution functions. There are several because each array represents
     ! a differently shaped chunk of memory. In this example, each chunk 
     ! allows sequential operations in one given direction. f_x1x2 should 
     ! permit to carry out sequential operations in x1 and x2 for ex.
     sll_real64, dimension(:,:,:,:), pointer     :: f_x1x2 
     sll_real64, dimension(:,:,:,:), pointer     :: f_x3x4
     sll_real64, dimension(:,:,:), allocatable   :: partial_reduction
     sll_real64, dimension(:,:), allocatable     :: rho_x1 
     sll_real64, dimension(:,:), allocatable     :: rho_x2 
     sll_real64, dimension(:,:), allocatable     :: rho_split
     sll_real64, dimension(:,:), allocatable     :: phi_x1
     sll_real64, dimension(:,:), allocatable     :: phi_x2
     sll_real64, dimension(:,:), allocatable     :: phi_split
     
     ! for remap
     type(layout_4D), pointer :: sequential_x1x2
     type(layout_4D), pointer :: sequential_x3x4
     type(layout_2D), pointer :: rho_seq_x1
     type(layout_2D), pointer :: rho_seq_x2
     type(layout_2D), pointer :: split_rho_layout ! not sequential in any dir.
     type(layout_2D), pointer :: split_phi_layout ! not sequential in any dir.
     type(remap_plan_2D_real64), pointer :: split_to_seqx1
     type(remap_plan_2D_real64), pointer :: seqx1_to_seqx2
     ! remaps for the electric field data
!     type(remap_plan_2D), pointer :: efld_split_to_seqx1
     type(remap_plan_2D_comp64), pointer :: efld_seqx1_to_seqx2
     type(remap_plan_2D_comp64), pointer :: efld_seqx2_to_split
     type(remap_plan_4D_real64), pointer :: seqx1x2_to_seqx3x4
     type(remap_plan_4D_real64), pointer :: seqx3x4_to_seqx1x2
     ! interpolators and their pointers
     type(cubic_spline_2d_interpolator) :: interp_x1x2
!!$     type(cubic_spline_1d_interpolator) :: interp_x1
!!$     type(cubic_spline_1d_interpolator) :: interp_x2
     type(cubic_spline_1d_interpolator) :: interp_x3
     type(cubic_spline_1d_interpolator) :: interp_x4
     ! Field accumulator
     sll_comp64, dimension(:,:), allocatable :: efield_x1
     sll_comp64, dimension(:,:), allocatable :: efield_x2
     sll_comp64, dimension(:,:), allocatable :: efield_split
     ! for distribution function initializer:
     procedure(sll_scalar_initializer_4d), nopass, pointer :: init_func
     sll_real64, dimension(:), pointer :: params
   contains
     procedure, pass(sim) :: run => run_vp4d_cartesian_general
     procedure, pass(sim) :: init_from_file => init_vp4d_par_gen
  end type sll_simulation_4d_vp_general

  interface sll_delete
     module procedure delete_vp4d_par_gen
  end interface sll_delete

  interface initialize
     module procedure initialize_vp4d_general
  end interface initialize

contains

  ! Tentative function to initialize the simulation object 'manually'.
  subroutine initialize_vp4d_general( &
   sim, &
   mesh2d_x, &
   mesh2d_v, &
   transformation_x, &
   init_func, &
   params )

   type(sll_simulation_4d_vp_general), intent(inout)     :: sim
   type(sll_logical_mesh_2d), pointer                    :: mesh2d_x
   type(sll_logical_mesh_2d), pointer                    :: mesh2d_v
   class(sll_coordinate_transformation_2d_base), pointer :: transformation_x
   procedure(sll_scalar_initializer_4d)                  :: init_func
   sll_real64, dimension(:), target                      :: params
   sim%mesh2d_x  => mesh2d_x
   sim%mesh2d_v  => mesh2d_v
   sim%transfx   => transformation_x
   sim%init_func => init_func
   sim%params    => params
  end subroutine initialize_vp4d_general


  subroutine init_vp4d_par_gen( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_vp_general), intent(inout) :: sim
    character(len=*), intent(in)                                   :: filename
    sll_int32             :: IO_stat
    sll_real64            :: dt
    sll_int32             :: number_iterations
    sll_int32             :: num_cells_x1
    sll_int32             :: num_cells_x2
    sll_int32             :: num_cells_x3
    sll_int32             :: num_cells_x4
    sll_int32, parameter  :: input_file = 99

    namelist /sim_params/ dt, number_iterations
    namelist /grid_dims/ num_cells_x1, num_cells_x2, num_cells_x3, num_cells_x4
    ! Try to add here other parameters to initialize the mesh values like
    ! xmin, xmax and also for the distribution function initializer.
    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_vp4d_par_cart() failed to open file ', filename
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
  end subroutine init_vp4d_par_gen

  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_vp4d_cartesian_general(sim)
    class(sll_simulation_4d_vp_general), intent(inout) :: sim
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: loc_sz_x3
    sll_int32  :: loc_sz_x4
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_real64 :: vmin3
    sll_real64 :: vmax3
    sll_real64 :: vmin4
    sll_real64 :: vmax4
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_real64 :: alpha1
    sll_real64 :: alpha2
    sll_real64 :: alpha3
    sll_real64 :: alpha4
    sll_int32  :: itemp
    sll_int32  :: ierr
    sll_int32  :: itime
    sll_int32  :: nc_x1
    sll_int32  :: nc_x2
    sll_int32  :: nc_x3
    sll_int32  :: nc_x4
    sll_real64 :: ex
    sll_real64 :: ey
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: eta3
    sll_real64 :: eta4
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta3_min
    sll_real64 :: eta4_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_max
    sll_real64 :: eta3_max
    sll_real64 :: eta4_max
    sll_real64 :: eta1_new
    sll_real64 :: eta2_new
    sll_real64 :: diff
    sll_real64, dimension(1:2,1:2) :: inv_j
    sll_real64, dimension(1:2,1:2) :: jac_m
    sll_int32, dimension(1:2)      :: gi     ! for storing global indices
    sll_int32, dimension(1:4)      :: gi4d   ! for storing global indices
    sll_real64 :: efield_energy_total
    ! The following could probably be abstracted for convenience
#define BUFFER_SIZE 10
    sll_real64, dimension(BUFFER_SIZE) :: buffer
    sll_real64, dimension(BUFFER_SIZE) :: buffer_result
    sll_real64, dimension(BUFFER_SIZE) :: num_particles_local
    sll_real64, dimension(BUFFER_SIZE) :: num_particles_global
    sll_real64 :: tmp
    sll_real64, dimension(1) :: tmp1
    sll_int32 :: buffer_counter
    sll_int32 :: efield_energy_file_id
    sll_int32 :: num_particles_file_id
    sll_int32 :: global_indices(4)
    sll_int32 :: iplot
    character(len=4) :: cplot
    ! only for debugging...
!!$    sll_real64, dimension(:,:), allocatable :: ex_field
!!$    sll_real64, dimension(:,:), allocatable :: ey_field
    print *, 'executing simulation'

    buffer_counter = 1

    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)

    if( sim%my_rank == 0 ) then
       call sll_new_file_id( efield_energy_file_id, ierr )
       if( ierr == 1 ) then
          print *, 'sll_new_file_id() failed to obtain a file identifier.', &
               ' Exiting...'
          stop
       end if
    end if

    if( sim%my_rank == 0 ) then
       call sll_new_file_id( num_particles_file_id, ierr )
       if( ierr == 1 ) then
          print *, 'sll_new_file_id() failed to obtain a file identifier.', &
               ' Exiting...'
          stop
       end if
    end if


    ! allocate the layouts...
    sim%sequential_x1x2  => new_layout_4D( sll_world_collective )
    sim%sequential_x3x4  => new_layout_4D( sll_world_collective )
    sim%rho_seq_x1       => new_layout_2D( sll_world_collective )
    sim%rho_seq_x2       => new_layout_2D( sll_world_collective )
    sim%split_rho_layout => new_layout_2D( sll_world_collective )
    sim%split_phi_layout => new_layout_2D( sll_world_collective )
    print *, 'completed layout allocation'

    ! layout for sequential operations in x3 and x4. Make an even split for
    ! x1 and x2, or as close as even if the power of 2 is odd. This should 
    ! be packaged in some sort of routine and set up at initialization time.
    sim%power2 = int(log(real(sim%world_size))/log(2.0))

    ! special case N = 1, so power2 = 0
    if(sim%power2 == 0) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    end if
    
    if(is_even(sim%power2)) then
       sim%nproc_x1 = 2**(sim%power2/2)
       sim%nproc_x2 = 2**(sim%power2/2)
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    else 
       sim%nproc_x1 = 2**((sim%power2-1)/2)
       sim%nproc_x2 = 2**((sim%power2+1)/2)
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    end if
    
    nc_x1 = sim%mesh2d_x%num_cells1
    nc_x2 = sim%mesh2d_x%num_cells2
    nc_x3 = sim%mesh2d_v%num_cells1
    nc_x4 = sim%mesh2d_v%num_cells2
    delta1 = sim%mesh2d_x%delta_eta1
    delta2 = sim%mesh2d_x%delta_eta2
    delta3 = sim%mesh2d_v%delta_eta1
    delta4 = sim%mesh2d_v%delta_eta2
    eta1_min = sim%mesh2d_x%eta1_min
    eta2_min = sim%mesh2d_x%eta2_min
    eta3_min = sim%mesh2d_v%eta1_min
    eta4_min = sim%mesh2d_v%eta2_min
    eta1_max = sim%mesh2d_x%eta1_max
    eta2_max = sim%mesh2d_x%eta2_max
    eta3_max = sim%mesh2d_v%eta1_max
    eta4_max = sim%mesh2d_v%eta2_max

    call initialize_layout_with_distributed_4D_array( &
         nc_x1+1, &
         nc_x2+1, &
         nc_x3+1, &
         nc_x4+1, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%nproc_x3, &
         sim%nproc_x4, &
         sim%sequential_x3x4 )

    ! Use this information to initialize the layout that describes the result
    ! of computing rho. This layout is not useful to do sequential operations
    ! in any of the two available directions. We also initialize the other two
    ! layouts needed for both sequential operations on x1 and x2 in the 2D case.
    call initialize_layout_with_distributed_2D_array( &
         nc_x1+1, &
         nc_x2+1, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%split_rho_layout )

    call initialize_layout_with_distributed_2D_array( &
         nc_x1+1, &
         nc_x2+1, &
         1, &
         sim%world_size, &
         sim%rho_seq_x1 )
    
    call compute_local_sizes_2d( sim%rho_seq_x1, loc_sz_x1, loc_sz_x2 )
    SLL_ALLOCATE(sim%rho_x1(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%phi_x1(loc_sz_x1,loc_sz_x2),ierr)
!    SLL_ALLOCATE(ex_field(loc_sz_x1,loc_sz_x2),ierr)
    ! Experiment with a dedicated array to store the values of the electric
    ! field in each point of the grid.
    SLL_ALLOCATE(sim%efield_x1(loc_sz_x1,loc_sz_x2),ierr)

    call initialize_layout_with_distributed_2D_array( &
         nc_x1+1, &
         nc_x2+1, &
         sim%world_size, &
         1, &
         sim%rho_seq_x2 )
    call compute_local_sizes_2d( sim%rho_seq_x2, loc_sz_x1, loc_sz_x2 )
    SLL_ALLOCATE(sim%rho_x2(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%phi_x2(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x2(loc_sz_x1,loc_sz_x2),ierr)
 !   SLL_ALLOCATE(ey_field(loc_sz_x1,loc_sz_x2),ierr)
    ! layout for sequential operations in x1 and x2. This is basically just the
    ! flipping of the values between x1,x3 and x2,x4 on the previous layout.
    ! This is the wrong (dangerous) way to do this. We should not be changing
    ! the values of these variables, which should behave like parameters. 
    ! Choose a better name for each, initialize and don't change any further.
    ! switch x1 and x3:
    itemp = sim%nproc_x3
    sim%nproc_x3 = sim%nproc_x1
    sim%nproc_x1 = itemp
    ! switch x2 and x4
    itemp = sim%nproc_x4
    sim%nproc_x4 = sim%nproc_x2 
    sim%nproc_x2 = itemp
    
    call initialize_layout_with_distributed_4D_array( &
         nc_x1+1, &
         nc_x2+1, &
         nc_x3+1, &
         nc_x4+1, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%nproc_x3, &
         sim%nproc_x4, &
         sim%sequential_x1x2 )
    
    ! Allocate the array needed to store the local chunk of the distribution
    ! function data. First compute the local sizes. Since the remap operations
    ! are out-of-place, we will allocate four different arrays, one for each
    ! layout.
    call compute_local_sizes_4d( sim%sequential_x1x2, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )
    SLL_ALLOCATE(sim%f_x1x2(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)
    SLL_ALLOCATE(sim%phi_split(loc_sz_x3,loc_sz_x4),ierr)
    
    
    call compute_local_sizes_4d( sim%sequential_x3x4, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )

    call initialize_layout_with_distributed_2D_array( &
         nc_x3, &
         nc_x4, &
         sim%nproc_x3, &
         sim%nproc_x4, &
         sim%split_phi_layout )

    ! This layout is also useful to represent the charge density array. Since
    ! this is a result of a local reduction on x3 and x4, the new layout is
    ! 2D but with the same dimensions of the process mesh in x1 and x2.
    SLL_ALLOCATE(sim%rho_split(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_split(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%f_x3x4(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)
    
    ! These dimensions are also the ones needed for the array where we store
    ! the intermediate results of the charge density computation.
    SLL_ALLOCATE(sim%partial_reduction(loc_sz_x1,loc_sz_x2, loc_sz_x3),ierr)

    print *, 'completed memory allocations'
    
    call sll_4d_parallel_array_initializer( &
         sim%sequential_x3x4, &
         sim%mesh2d_x, &
         sim%mesh2d_v, &
         sim%f_x3x4, &
         sim%init_func, &
         sim%params, &
         transf_x1_x2=sim%transfx )

    print *, 'initialized the distribution function'

    delta1 = sim%mesh2d_x%delta_eta1
    delta2 = sim%mesh2d_x%delta_eta2
    delta3 = sim%mesh2d_v%delta_eta1
    delta4 = sim%mesh2d_v%delta_eta2

    call compute_charge_density( sim%mesh2d_x,           &
                                 sim%mesh2d_v,           &
                                 size(sim%f_x3x4,1),     &
                                 size(sim%f_x3x4,2),     &
                                 sim%f_x3x4,             &
                                 sim%partial_reduction,  &
                                 sim%rho_split )

    print *, 'computed charge density'

    sim%split_to_seqx1 => &
         NEW_REMAP_PLAN(sim%split_rho_layout, sim%rho_seq_x1, sim%rho_split)
    call apply_remap_2D( sim%split_to_seqx1, sim%rho_split, sim%rho_x1 )
       
    global_indices(1:2) =  &
         local_to_global_2D( sim%rho_seq_x1, (/1, 1/) )

    print *, 'proceeding to write rho_x1 to a file'
! This program is crashing at this point. Why is 
!    if(sim%my_rank == 0) then
    call compute_local_sizes_2d( sim%rho_seq_x1, loc_sz_x1, loc_sz_x2 )

       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min + &
               (global_indices(1)-1)*sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min + &
            (global_indices(2)-1)*sim%mesh2d_x%delta_eta2, &
            sim%mesh2d_x%delta_eta2, &
            loc_sz_x1, loc_sz_x2, &
            sim%rho_x1, &
            "rho_x1", &
            0, &
            ierr )
       print *, 'printed rho_x1'
 !   end if

    ! With the distribution function initialized in at least one configuration,
    ! we can proceed to carry out the computation of the electric potential.
    ! First we need to compute the charge density. Some thoughts:
    !
    ! The computation of rho is a reduction process that takes as input a 4d
    ! array and that should return a 2d array (or alternatively, a 4d array
    ! of size 1 in the reduced directions). For example, a df of dimensions
    ! np1 X np2 X np3 X np4 might effectively end up as an array of dimensions
    ! np1 X np2 X np3 X 1  after a summation of all the values in x4. After a
    ! second summation along x3, the dimensions would be np1 X np2 X 1  X 1. 
    ! One simple-minded but inefficient way to prepare the data for the double
    ! reduction could be to have a layout in a process mesh NP1 X NP2 X 1 X 1
    ! where NP1xNP2 is the total number of processors available. The problem 
    ! here is that the end result would still need a layout change to be fed
    ! into the Poisson solver...
    !
    ! So can we do better? Let try working backwards. The desired input for
    ! the Poisson solver is a 2d array where sequential operations in x1 are
    ! possible. Hence, the last reduction operation
    !
    ! Let's start with the idea that we want to maintain the same number of
    ! processors busy when we launch the Poisson solver. This means that if the
    ! original process mesh for the 4d data is NP1xNP2xNP3x1 (we start with a
    ! layout that permits a reduction in x4), then the processor mesh for the
    ! Poisson step should be NP1'xNP2'x1x1 where NP1xNP2xNP3 = NP1'xNP2'

    sim%seqx3x4_to_seqx1x2 => &
         NEW_REMAP_PLAN(sim%sequential_x3x4,sim%sequential_x1x2,sim%f_x3x4)

    sim%seqx1x2_to_seqx3x4 => &
         NEW_REMAP_PLAN(sim%sequential_x1x2,sim%sequential_x3x4,sim%f_x1x2)
    
    call apply_remap_4D( sim%seqx3x4_to_seqx1x2, sim%f_x3x4, sim%f_x1x2 )
    call compute_local_sizes_4d( sim%sequential_x1x2, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )
    
    ! Carry out the 2D advection in the eta1-eta2 plane.
    ! dt in eta1 and eta2
    vmin3  = sim%mesh2d_v%eta1_min
    vmax3  = sim%mesh2d_v%eta1_max
    vmin4  = sim%mesh2d_v%eta2_min
    vmax4  = sim%mesh2d_v%eta2_max

    ! First dt/2 advection for eta1-eta2:
    
    ! compute the spline coefficients
    ! Start the interpolators... Watch out: the periodic case has equal number
    ! of cells than points. Is this properly handled by the interpolators??
    ! The interpolators need the number of points and always consider that
    ! num_cells = num_pts - 1. This is a possible source of confusion.
    
    ! endpoints are hardcoded because the logical mesh lives always on the 
    ! unit square. But with the more general logical meshes, this is not
    ! true anymore, so this should be changed to depend on the values stored
    ! in the logical grids.
    call sim%interp_x1x2%initialize( &
         nc_x1+1, &
         nc_x2+1, &
         sim%mesh2d_x%eta1_min, &
         sim%mesh2d_x%eta1_max, &
         sim%mesh2d_x%eta2_min, &
         sim%mesh2d_x%eta2_max, &
         SLL_PERIODIC, &
         SLL_PERIODIC )

    print *, 'proceeding to carry out first advection'
    call advection_x1x2(sim,0.5*sim%dt)
    print *, 'finished first advection in x1 and x2'

    ! other test cases use periodic bc's here...        
    call sim%interp_x3%initialize( &
         nc_x3+1, &
         vmin3, &
         vmax3, &
         SLL_PERIODIC)

    call sim%interp_x4%initialize( &
         nc_x4+1, &
         vmin4, &
         vmax4, &
         SLL_PERIODIC)

    call compute_local_sizes_2d( sim%rho_seq_x1, loc_sz_x1, loc_sz_x2 )

    ! Initialize the poisson plan before going into the main loop.

    sim%poisson_plan => new_poisson_2d_periodic_plan_cartesian_par( &
            sim%rho_seq_x1, &
            nc_x1, &
            nc_x2, &
            sim%mesh2d_x%eta1_max - sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%eta2_max - sim%mesh2d_x%eta2_min )

    sim%efld_seqx1_to_seqx2 => &
          NEW_REMAP_PLAN( sim%rho_seq_x1, sim%rho_seq_x2, sim%efield_x1)

    sim%seqx1_to_seqx2 => &
         NEW_REMAP_PLAN( sim%rho_seq_x1, sim%rho_seq_x2, sim%phi_x1 )

    sim%efld_seqx2_to_split => &
         NEW_REMAP_PLAN( sim%rho_seq_x2, sim%split_rho_layout,sim%efield_x2 )
    ! ------------------------------------------------------------------------
    !
    !                                MAIN LOOP
    !
    ! ------------------------------------------------------------------------


    do itime=1, sim%num_iterations
       if(sim%my_rank == 0) then
          print *, 'Starting iteration ', itime, ' of ', sim%num_iterations
       end if
       ! The splitting scheme used here is meant to attain a dt^2 accuracy.
       ! We use:
       
       ! dt/2   in eta1 and eta2
       ! solve the poisson equation (this acts as a prediction step for E in 
       ! BSL)
       !
       ! dt in vx
       ! dt in vy
       ! dt/2   in eta1 and eta2

       ! Here we join the first and last steps in eta1 and eta2 and to one
       ! extra dt/2 step outside of the main loop:

       ! dt/2   in eta1 and eta2 (outside of loop)
       !
       ! in main loop:
       !
       ! solve the poisson equation (this acts as a prediction step for E in 
       ! BSL)
       !
       ! dt in vx
       ! dt in vy
       ! dt in eta1 and eta2

       ! Note: Since the Ex and Ey values are used separately, the proposed
       ! data structure is actually not good. These field values should be kept
       ! separate.
       call apply_remap_4D( sim%seqx1x2_to_seqx3x4, sim%f_x1x2, sim%f_x3x4 )

       call compute_local_sizes_4d( sim%sequential_x1x2, &
                                    loc_sz_x1,           &
                                    loc_sz_x2,           &
                                    loc_sz_x3,           &
                                    loc_sz_x4 )

       sim%rho_split(:,:) = 0.0_f64
       call compute_charge_density( sim%mesh2d_x,           &
                                    sim%mesh2d_v,           &
                                    size(sim%f_x3x4,1),     &
                                    size(sim%f_x3x4,2),     &
                                    sim%f_x3x4,             &
                                    sim%partial_reduction,  &
                                    sim%rho_split )



       ! Re-arrange rho_split in a way that permits sequential operations in 
       ! x1, to feed to the Poisson solver.
       call apply_remap_2D( sim%split_to_seqx1, sim%rho_split, sim%rho_x1 )

       ! It is important to remember a particular property of the periodic 
       ! poisson solver used here: For a given input charge density 
       ! configuration, there is an infinite number of solutions for the 
       ! potential, all differing by a constant. By design, the solver will
       ! return the solution with ZERO-MEAN. This implies that if given as
       ! input a uniform blob of charge, the solver will return a potential
       ! field with value equal to zero everywhere. In other words, the solver
       ! will only 'react' to nonuniform inputs.
       !
       ! The above is important when dealing with certain tests or certain
       ! model assumptions. For example in a Landau-damping test, one may
       ! want to apply a uniform neutralizing field for instance by adding or
       ! subtracting a constant value. This step becomes unnecessary since
       ! the answer returned by the solver will not be affected.
       !
       !    sim%rho_x1(:,:) = sim%rho_x1(:,:) - 1.0_f64  <--- unnecessary step

       ! solve for the electric potential
       ! but first adjust the sign of the source term (the sign should come
       ! automatically from the computation of the charge density, please fix)
!       sim%rho_x1(:,:) = - sim%rho_x1(:,:) 

       global_indices(1:2) =  &
            local_to_global_2D( sim%rho_seq_x1, (/1, 1/) )
       
       call sll_gnuplot_rect_2d_parallel( &
          sim%mesh2d_x%eta1_min+(global_indices(1)-1)*sim%mesh2d_x%delta_eta1, &
          sim%mesh2d_x%delta_eta1, &
          sim%mesh2d_x%eta2_min+(global_indices(2)-1)*sim%mesh2d_x%delta_eta2, &
          sim%mesh2d_x%delta_eta2, &
          size(sim%rho_x1,1), &
          size(sim%rho_x1,2), &
          sim%rho_x1, &
          "rho_x1", &
          itime, &
          ierr )

       call solve_poisson_2d_periodic_cartesian_par( &
            sim%poisson_plan, &
            sim%rho_x1, &
            sim%phi_x1)
 

       global_indices(1:2) =  &
            local_to_global_2D( sim%rho_seq_x1, (/1, 1/) )
       
       call sll_gnuplot_rect_2d_parallel( &
         sim%mesh2d_x%eta1_min+(global_indices(1)-1)*sim%mesh2d_x%delta_eta1, &
         sim%mesh2d_x%delta_eta1, &
         sim%mesh2d_x%eta2_min+(global_indices(2)-1)*sim%mesh2d_x%delta_eta2, &
         sim%mesh2d_x%delta_eta2, &
         size(sim%phi_x1,1), &
         size(sim%phi_x1,2), &
         sim%phi_x1, &
         "phi_x1", &
         itime, &
         ierr )


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       ! Here is where Aurore's general geometry poisson should be included...
       !
       !???????????????????????????????????????????????????????????????????????
       

       ! compute the values of the electric field. rho is configured for 
       ! sequential operations in x1, thus we start by computing the E_x 
       ! component.
       ! The following call is inefficient and unnecessary. The local sizes for
       ! the arrays should be kept around as parameters basically and not on 
       ! variables whose content could be anything... This will have to do for 
       ! now.
       call compute_local_sizes_2d( sim%rho_seq_x1, loc_sz_x1, loc_sz_x2 )

       call compute_electric_field_eta1( &
            sim%phi_x1, &
            loc_sz_x1, &
            loc_sz_x2, &
            sim%mesh2d_x%delta_eta1, &
            sim%efield_x1 )

!!$       ex_field(:,:) = real(sim%efield_x1(:,:))
!!$
!!$
!!$       global_indices(1:2) =  &
!!$            local_to_global_2D( sim%rho_seq_x1, (/1, 1/) )
!!$       
!!$       call sll_gnuplot_rect_2d_parallel( &
!!$         sim%mesh2d_x%eta1_min+(global_indices(1)-1)*sim%mesh2d_x%delta_eta1, &
!!$         sim%mesh2d_x%delta_eta1, &
!!$         sim%mesh2d_x%eta2_min+(global_indices(2)-1)*sim%mesh2d_x%delta_eta2, &
!!$         sim%mesh2d_x%delta_eta2, &
!!$         ex_field, &
!!$         "ex_x1", &
!!$         itime, &
!!$         ierr )
!!$




       ! note that we are 'recycling' the layouts used for the other arrays 
       ! because they represent an identical configuration.
       call apply_remap_2D( sim%efld_seqx1_to_seqx2, sim%efield_x1, sim%efield_x2 )
       call apply_remap_2D( sim%seqx1_to_seqx2, sim%phi_x1, sim%phi_x2 )
       call compute_local_sizes_2d( sim%rho_seq_x2, loc_sz_x1, loc_sz_x2 )

       call compute_electric_field_eta2( &
            sim%phi_x2, &
            loc_sz_x1, &
            loc_sz_x2, &
            sim%mesh2d_x%delta_eta2, &
            sim%efield_x2 )
       
!!$       ey_field(:,:) = aimag(sim%efield_x2(:,:))
!!$
!!$       global_indices(1:2) =  &
!!$            local_to_global_2D( sim%rho_seq_x2, (/1, 1/) )
!!$       
!!$       call sll_gnuplot_rect_2d_parallel( &
!!$         sim%mesh2d_x%eta1_min+(global_indices(1)-1)*sim%mesh2d_x%delta_eta1, &
!!$         sim%mesh2d_x%delta_eta1, &
!!$         sim%mesh2d_x%eta2_min+(global_indices(2)-1)*sim%mesh2d_x%delta_eta2, &
!!$         sim%mesh2d_x%delta_eta2, &
!!$         ey_field, &
!!$         "ey_x2", &
!!$         itime, &
!!$         ierr )



       ! But now, to make the electric field data configuration compatible with
       ! the sequential operations in x2x3 we need still another remap 
       ! operation.
       call apply_remap_2D( &
            sim%efld_seqx2_to_split, &
            sim%efield_x2, &
            sim%efield_split )
       ! Now we proceed to reconfigure the data. This is very expensive. 
       ! There might be advantages to this approach if we avoid larger data 
       ! transfers like with an all-to-all transfer... however, we could end 
       ! up paying more if the simulation is latency-dominated.
      
       call compute_local_sizes_4d( sim%sequential_x1x2, &
                                    loc_sz_x1,           &
                                    loc_sz_x2,           &
                                    loc_sz_x3,           &
                                    loc_sz_x4 )

       call compute_local_sizes_4d( sim%sequential_x3x4, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 

       efield_energy_total = 0.0_f64
       
       ! Start with dt in vx...(x3)
       do l=1,sim%mesh2d_v%num_cells2 + 1
          do j=1,loc_sz_x2
             do i=1,loc_sz_x1
                global_indices(1:2) = &
                     local_to_global_2D( sim%split_rho_layout, (/i,j/))
                eta1   =  real(global_indices(1)-1,f64)*delta1
                eta2   =  real(global_indices(2)-1,f64)*delta2
                inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                jac_m  =  sim%transfx%jacobian_matrix(eta1,eta2)
                ex     =  real( sim%efield_split(i,j),f64)
                ey     =  aimag(sim%efield_split(i,j))
                alpha3 = -sim%dt*(inv_j(1,1)*ex + inv_j(2,1)*ey)
                sim%f_x3x4(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                     nc_x3+1, &
                     sim%f_x3x4(i,j,:,l), &
                     alpha3 )
                ! Extra work to calculate the electric field energy. We should 
                ! consider placing this somewhere else, probably at greater
                ! expense.
!print *, 'while calculating E energy: ', sim%transfx%jacobian_at_node(global_indices(1),global_indices(2)), delta1, delta2, inv_j(1,1), inv_j(2,1), inv_j(1,2), inv_j(2,2), ex, ey
!!$                efield_energy_total = efield_energy_total + &
!!$                     sim%transfx%jacobian_at_node(global_indices(1),&
!!$                                                  global_indices(2))*&
!!$                                                  delta1*delta2*&
!!$                     sqrt(((inv_j(1,1)*ex + inv_j(2,1)*ey)**2 + &
!!$                      (inv_j(1,2)*ex + inv_j(2,2)*ey)**2))

               efield_energy_total = efield_energy_total + &
                    delta1*delta2*sqrt(((jac_m(1,1)*ex + jac_m(2,1)*ey)**2 + &
                                        (jac_m(1,2)*ex + jac_m(2,2)*ey)**2))
             end do
          end do
       end do

       call compute_local_sizes_4d( sim%sequential_x3x4, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 

       ! dt in vy...(x4)
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             do k=1,sim%mesh2d_v%num_cells1 + 1
                global_indices(1:2) = &
                     local_to_global_2D( sim%split_rho_layout, (/i,j/))
                eta1   =  real(global_indices(1)-1,f64)*delta1
                eta2   =  real(global_indices(2)-1,f64)*delta2
                inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                ex     =  real( sim%efield_split(i,j),f64)
                ey     =  aimag(sim%efield_split(i,j))
                alpha4 = -sim%dt*(inv_j(1,2)*ex + inv_j(2,2)*ey)
                sim%f_x3x4(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
                     nc_x4+1, &
                     sim%f_x3x4(i,j,k,:), &
                     alpha4 )
             end do
          end do
       end do

       call apply_remap_4D( sim%seqx3x4_to_seqx1x2, sim%f_x3x4, sim%f_x1x2 )

       call compute_local_sizes_4d( sim%sequential_x1x2, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 

!!$       do l = 1, loc_sz_x4
!!$          do k = 1, loc_sz_x3
!!$              sim%phi_split(k,l) = sum(sim%f_x1x2(:,:,k,l))
!!$          end do
!!$       end do

!!$       global_indices(1:2) =  &
!!$            local_to_global_2D( sim%split_phi_layout, (/1, 1/) )
!!$       
!!$       call sll_gnuplot_rect_2d_parallel( &
!!$         sim%mesh2d_v%eta1_min+(global_indices(1)-1)*sim%mesh2d_v%delta_eta1, &
!!$         sim%mesh2d_v%delta_eta1, &
!!$         sim%mesh2d_v%eta2_min+(global_indices(2)-1)*sim%mesh2d_v%delta_eta2, &
!!$         sim%mesh2d_v%delta_eta2, &
!!$         sim%phi_split, &
!!$         "phi_split", &
!!$         itime, &
!!$         ierr )


       ! Approximate the integral of the distribution function along all
       ! directions.
       num_particles_local(buffer_counter) = &
            sum(sim%f_x3x4)*delta1*delta2*delta3*delta4
       if( buffer_counter == BUFFER_SIZE ) then
          call sll_collective_reduce_real64( &
               sll_world_collective, &
               num_particles_local, &
               BUFFER_SIZE, &
               MPI_SUM, &
               0, &
               num_particles_global )
          if(sim%my_rank == 0) then
             open(num_particles_file_id,file="number_particles",&
                  position="append")
             if(itime == BUFFER_SIZE) then
                rewind(num_particles_file_id)
             end if
             do i=1,BUFFER_SIZE
                write(num_particles_file_id,*) num_particles_global(i)
             end do
             close(num_particles_file_id)
          end if
       end if


       buffer(buffer_counter) = efield_energy_total
       ! This should be abstracted away...
       ! Each processor keeps a local buffer, when the buffer reaches a
       ! predetermined size, we reduce ther buffer with an addition on 
       ! process 0, who appends it to a file. Then reset the buffer.
       if( buffer_counter == BUFFER_SIZE ) then
          call sll_collective_reduce_real64( &
               sll_world_collective, &
               buffer, &
               BUFFER_SIZE, &
               MPI_SUM, &
               0, &
               buffer_result )

          buffer_counter = 1
          if(sim%my_rank == 0) then
             open(efield_energy_file_id,file="electric_field_energy",&
                  position="append")
             if(itime == BUFFER_SIZE) then
                rewind(efield_energy_file_id)
             end if
             buffer_result(:) = log(buffer_result(:))
             do i=1,BUFFER_SIZE
                write(efield_energy_file_id,*) buffer_result(i)
             end do
             close(efield_energy_file_id)
          end if
       else
          buffer_counter         = buffer_counter + 1
       end if
       efield_energy_total    = 0.0_f64
       ! Proceed to the advections in the spatial directions, 'x' and 'y'
       ! Reconfigure data. 
       
       ! what are the new local limits on x3 and x4? It is bothersome to have
       ! to make these calls...

       call advection_x1x2(sim,sim%dt)
       call apply_remap_4D( sim%seqx1x2_to_seqx3x4, sim%f_x1x2, sim%f_x3x4 )
       call compute_local_sizes_4d( sim%sequential_x3x4, &
                                 loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 

!!$       do j = 1, loc_sz_x2
!!$          do i = 1, loc_sz_x1
!!$             sim%rho_split(i,j) = sum(sim%f_x3x4(i,j,:,:))
!!$          end do
!!$       end do
!!$
!!$      global_indices(1:2) =  local_to_global_2D( sim%split_rho_layout, (/1, 1/) )
!!$
!!$      call sll_gnuplot_rect_2d_parallel( &
!!$           sim%mesh2d_x%eta1_min+(global_indices(1)-1)*delta1,delta1, &
!!$           sim%mesh2d_x%eta2_min+(global_indices(2)-1)*delta2,delta2, &
!!$           sim%rho_split,"rho_split",itime,ierr )

  !     call plot_fields(itime, sim)
#undef BUFFER_SIZE

   end do ! main loop

  end subroutine run_vp4d_cartesian_general


  subroutine advection_x1x2(sim,deltat)
    class(sll_simulation_4d_vp_general) :: sim
    sll_real64, intent(in) :: deltat
    sll_int32 :: gi, gj, gk, gl
    sll_real64, dimension(1:2,1:2) :: inv_j
    sll_int32 :: loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 
    sll_int32, dimension(4) :: global_indices
    sll_real64 :: alpha1, alpha2
    sll_real64 :: eta1, eta2, eta3, eta4
    sll_int32 :: i, j, k, l

    call compute_local_sizes_4d( sim%sequential_x1x2, &
         loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 )

    do l=1,loc_sz_x4
       do k=1,loc_sz_x3
          call sim%interp_x1x2%compute_interpolants(sim%f_x1x2(:,:,k,l))
          do j=1,loc_sz_x2
             do i=1,loc_sz_x1
                global_indices = local_to_global_4D(sim%sequential_x1x2,(/i,j,k,l/))
                gi = global_indices(1)
                gj = global_indices(2)
                gk = global_indices(3)
                gl = global_indices(4)
                eta1 = sim%mesh2d_x%eta1_min + (gi-1)*sim%mesh2d_x%delta_eta1
                eta2 = sim%mesh2d_x%eta2_min + (gj-1)*sim%mesh2d_x%delta_eta2
                eta3 = sim%mesh2d_v%eta1_min + (gk-1)*sim%mesh2d_v%delta_eta1
                eta4 = sim%mesh2d_v%eta2_min + (gl-1)*sim%mesh2d_v%delta_eta2
                inv_j  = sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                alpha1 = -deltat*(inv_j(1,1)*eta3 + inv_j(1,2)*eta4)
                alpha2 = -deltat*(inv_j(2,1)*eta3 + inv_j(2,2)*eta4)

                eta1 = eta1+alpha1
                ! This is hardwiring the periodic BC, please improve this...
                if( eta1 <  sim%mesh2d_x%eta1_min ) then
                   eta1 = eta1+sim%mesh2d_x%eta1_max-sim%mesh2d_x%eta1_min
                else if( eta1 >  sim%mesh2d_x%eta1_max ) then
                   eta1 = eta1+sim%mesh2d_x%eta1_min-sim%mesh2d_x%eta1_max
                end if

                eta2 = eta2+alpha2
                if( eta2 <  sim%mesh2d_x%eta2_min ) then
                   eta2 = eta2+sim%mesh2d_x%eta2_max-sim%mesh2d_x%eta2_min
                else if( eta2 >  sim%mesh2d_x%eta2_max ) then
                   eta2 = eta2+sim%mesh2d_x%eta2_min-sim%mesh2d_x%eta2_max
                end if
                
                sim%f_x1x2(i,j,k,l) = sim%interp_x1x2%interpolate_value(eta1,eta2)

!!$             alpha1 = -(vmin3 + (k-1)*delta3)*sim%dt*0.5_f64
!!$             alpha2 = -(vmin4 + (l-1)*delta4)*sim%dt*0.5_f64
!!$             !call sim%interp_x1%compute_interpolants( sim%f_x1x2(i,:,k,l) )
!!$             ! interpolate_array_disp() has an interface that must be changed
!!$             sim%f_x1x2(:,:,k,l) = sim%interp_x1x2%interpolate_array_disp( &
!!$                  nc_x1, &
!!$                  nc_x2 , &
!!$                  sim%f_x1x2(:,:,k,l), &
!!$                  alpha1, &
!!$                  alpha2 )

             end do
          end do
       end do
    end do

  end subroutine advection_x1x2

  subroutine delete_vp4d_par_gen( sim )
    type(sll_simulation_4d_vp_general) :: sim
    sll_int32 :: ierr
    SLL_DEALLOCATE( sim%f_x1x2, ierr )
    SLL_DEALLOCATE( sim%f_x3x4, ierr )
    SLL_DEALLOCATE_ARRAY( sim%partial_reduction, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_split, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_split, ierr )
    call sll_delete( sim%sequential_x1x2 )
    call sll_delete( sim%sequential_x3x4 )
    call sll_delete( sim%rho_seq_x1 )
    call sll_delete( sim%rho_seq_x2 )
    call sll_delete( sim%split_rho_layout )
    call sll_delete( sim%split_to_seqx1 )
    call sll_delete( sim%efld_seqx1_to_seqx2 )
    call sll_delete( sim%efld_seqx2_to_split )
    call sll_delete( sim%seqx1x2_to_seqx3x4 )
    call sll_delete( sim%seqx3x4_to_seqx1x2 )
    call delete( sim%interp_x1x2 )
    call delete( sim%interp_x3 )
    call delete( sim%interp_x4 )
    SLL_DEALLOCATE_ARRAY( sim%efield_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%efield_x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%efield_split, ierr )
  end subroutine delete_vp4d_par_gen

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
  subroutine compute_charge_density( &
    mx, &
    mv, &
    numpts1, &
    numpts2, &
    f, &
    partial, &
    rho )

    type(sll_logical_mesh_2d), pointer     :: mx
    type(sll_logical_mesh_2d), pointer     :: mv
    sll_int32, intent(in)                  :: numpts1
    sll_int32, intent(in)                  :: numpts2
    sll_real64, intent(in),  dimension(:,:,:,:) :: f       ! local distr. func
    sll_real64, intent(inout),  dimension(:,:,:):: partial ! intermediate res.
    sll_real64, intent(inout), dimension(:,:)     :: rho     ! local rho
    ! local sizes in the split directions have to be given by caller.
    sll_int32                       :: numpts3
    sll_int32                       :: numpts4
    sll_real64                      :: delta3
    sll_real64                      :: delta4

    sll_int32 :: i, j, k, l
    
    numpts3 = mv%num_cells1
    numpts4 = mv%num_cells2
    delta3  = mv%delta_eta1
    delta4  = mv%delta_eta2

    partial(:,:,:) = 0.0
    
    ! This expects partial to be already initialized to zero!!!
    do k=1,numpts3
       do j=1,numpts2
          do i=1,numpts1
             ! This summation happens on a super-long stride... slow stuff
             ! This loop should be substituted by a proper integration
             ! function that we could use in the other directions as well...
             do l=1,numpts4
                partial(i,j,k) = partial(i,j,k) + f(i,j,k,l)*delta4
             end do
          end do
       end do
    end do
    
    ! Carry out the final reduction on x3. Note that rho is not initialized
    ! to zero since it may already have the partial charge accumulation from
    ! other species.
    do j=1,numpts2
       do i=1,numpts1
          do k=1,numpts3
             ! This summation happens on a very-long stride... slow stuff
             ! This loop should be substituted by a proper integration
             ! function that we could use in the other directions as well.
             ! See above reduction function for same problem.
             rho(i,j) = rho(i,j) + partial(i,j,k)*delta3
          end do
       end do
    end do
  end subroutine compute_charge_density
  
!!$  ! Super-ugly and ad-hoc but this needs to be done if we use the last point
!!$  ! to compute the electric fields.
!!$  subroutine ensure_periodicity_x(phi, how_many, nptsx)
!!$    sll_real64, dimension(:,:), intent(inout) :: phi
!!$    sll_int32, intent(in) :: how_many
!!$    sll_int32, intent(in) :: nptsx
!!$    sll_int32 :: j
!!$
!!$    do j=1,how_many
!!$       phi(nptsx,j) = phi(1,j)
!!$    end do
!!$  end subroutine ensure_periodicity_x
!!$
!!$  subroutine ensure_periodicity_y(phi, how_many, nptsy)
!!$    sll_real64, dimension(:,:), intent(inout) :: phi
!!$    sll_int32, intent(in) :: how_many
!!$    sll_int32, intent(in) :: nptsy
!!$    sll_int32 :: i
!!$
!!$    do i=1,how_many
!!$       phi(i,nptsy) = phi(i,1)
!!$    end do
!!$  end subroutine ensure_periodicity_y
!!$

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
    
    r_delta = 0.5_f64/delta

    ! Do first point:
    efield(1) = 2.0*r_delta*(-1.5_f64*phi(1) + 2.0_f64*phi(2) - 0.5_f64*phi(3))
    
    ! Do the internal values:
    do i=2,num_pts-1
       efield(i) = r_delta*(phi(i+1) - phi(i-1))
    end do
    
    ! Do last point:
    efield(num_pts) = 2.0*r_delta*( 0.5_f64*phi(num_pts-2) - &
         2.0_f64*phi(num_pts-1) + &
         1.5_f64*phi(num_pts) )
  end subroutine compute_electric_field_on_line
  
  subroutine compute_electric_field_eta1( &
       phi_x1, &
       num_pts_x1, &
       num_pts_x2, &
       delta_x1, &
       efield_x1 )

    intrinsic                               :: cmplx
    sll_real64, dimension(:,:), intent(in)  :: phi_x1
    sll_int32, intent(in)                   :: num_pts_x1
    sll_int32, intent(in)                   :: num_pts_x2
    sll_real64, intent(in)                  :: delta_x1
    sll_comp64, dimension(:,:), intent(out) :: efield_x1
    sll_int32                               :: i
    sll_int32                               :: j
    sll_real64                              :: r_delta
    sll_real64                              :: ex
    ! FIXME: arg checking
    
    r_delta = 1.0_f64/delta_x1
    
    ! Compute the electric field values on the left and right edges.
    do j=1,num_pts_x2
       ! left:
       ex = -r_delta*(-1.5_f64*phi_x1(1,j) + &
                       2.0_f64*phi_x1(2,j) - &
                       0.5_f64*phi_x1(3,j) )
       efield_x1(1,j) = cmplx(ex,f64)  
       ! right:
       ex = -r_delta*(0.5_f64*phi_x1(num_pts_x1-2,j)-&
                      2.0_f64*phi_x1(num_pts_x1-1,j)+&
                      1.5_f64*phi_x1(num_pts_x1,j) )
       efield_x1(num_pts_x1,j) = cmplx(ex,f64) 
    end do
    
    ! Electric field in interior points
    do j=1,num_pts_x2
       do i=2, num_pts_x1-1
          ex = -r_delta*0.5_f64*(phi_x1(i+1,j) - phi_x1(i-1,j))
          efield_x1(i,j) = cmplx(ex,f64)
       end do
    end do
  end subroutine compute_electric_field_eta1
  
  ! This function only sets the Ey component of the electric field. NOTE: This
  ! and the above function have a terrible flaw: they are not independent. 
  ! compute_electric_field_x2 assumes that compute_electric_field_x1 has been
  ! called first, thus the ex values have been computed. 
  ! compute_electric_field_x1() does not have a similar assumption. It would
  ! be good to change things in order to get rid of these functions.
  subroutine compute_electric_field_eta2( &
       phi_x2, &
       num_pts_x1, &
       num_pts_x2, &
       delta_x2, &
       efield_x2 )
    
    intrinsic                               :: cmplx
    sll_real64, dimension(:,:), intent(in)  :: phi_x2
    sll_int32, intent(in)                   :: num_pts_x1
    sll_int32, intent(in)                   :: num_pts_x2
    sll_real64, intent(in)                  :: delta_x2
    sll_comp64, dimension(:,:), intent(out) :: efield_x2
    sll_int32                               :: i
    sll_int32                               :: j
    sll_real64                              :: r_delta
    sll_real64                              :: ex
    sll_real64                              :: ey

    ! FIXME: arg checking
    
    r_delta = 1.0_f64/delta_x2
    
    ! Compute the electric field values on the bottom and top edges.
    do i=1,num_pts_x1
       ! bottom:
       ! ... there has to be a better way to do this in fortran :-(
       ex = real(efield_x2(i,1),f64)
       ey = -r_delta*(-1.5_f64*phi_x2(i,1) + 2.0_f64*phi_x2(i,2) - &
                      0.5_f64*phi_x2(i,3))
       efield_x2(i,1) = cmplx(ex,ey,f64)
       ! top:
       ex = real(efield_x2(i,num_pts_x1),f64)
       ey = -r_delta*(0.5_f64*phi_x2(i,num_pts_x2-2) - &
                     2.0_f64*phi_x2(i,num_pts_x2-1)+&
                     1.5_f64*phi_x2(i,num_pts_x2))
       efield_x2(i,num_pts_x2) = cmplx(ex,ey,f64) 
    end do
    
    ! Electric field in interior points
    do j=2,num_pts_x2-1
       do i=1, num_pts_x1
          ex = real(efield_x2(i,j),f64)
          ey = -r_delta*0.5_f64*(phi_x2(i,j+1) - phi_x2(i,j-1))
          efield_x2(i,j) = cmplx(ex,ey,f64) 
       end do
    end do
  end subroutine compute_electric_field_eta2
  
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
       displacement = -efield(i)*0.5_f64*dt
       f_line = f_interp%interpolate_array_disp(num_pts, f_line, displacement)
    end do
  end subroutine advection_v_1d

#if 0
  subroutine plot_fields(itime, sim)
    use sll_collective
    use hdf5
    use sll_hdf5_io_parallel
    use sll_xml_io
    sll_int32, intent(in) :: itime
    character(len=4)      :: ctime
    sll_int32             :: i_layout
    character(len=1)      :: c_layout
    type(sll_simulation_4d_vp_general), intent(in) :: sim
    type(layout_2D), pointer :: my_layout
    character(len=7),  parameter :: hdf_file = "data.h5"  ! File name
    sll_real64 :: tcpu1, tcpu2
    sll_int32  :: my_rank
    sll_int32  :: world_size
    sll_int32  :: local_nx1
    sll_int32  :: local_nx2
    sll_int32  :: global_nx1
    sll_int32  :: global_nx2
    sll_int32  :: error
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: gi
    sll_int32  :: gj
    sll_int32,  dimension(2) :: global_indices
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_real64 :: x1_min
    sll_real64 :: x1_max
    sll_real64 :: x2_min
    sll_real64 :: x2_max
    sll_real64 :: x3_min
    sll_real64 :: x3_max
    sll_real64 :: x4_min
    sll_real64 :: x4_max
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: delta_x3
    sll_real64 :: delta_x4 

    integer(HID_T)                  :: hdf_file_id
    sll_int32                       :: xml_file_id
    integer(HSIZE_T), dimension(2)  :: array_dims 
    integer(HSSIZE_T), dimension(2) :: offset 

    array_dims(1) = nc_x1
    array_dims(2) = nc_x2
    world_size    = sll_get_collective_size(sll_world_collective)
    my_rank       = sll_get_collective_rank(sll_world_collective)

    tcpu1 = MPI_WTIME()

    do i_layout = 1, 2

       if (i_layout == 1) then
          my_layout => sim%rho_seq_x1
       else
          my_layout => sim%rho_seq_x2
       end if

       call compute_local_sizes_2d( my_layout, local_nx1, local_nx2)        
    
       offset(1) =  get_layout_2D_i_min( my_layout, my_rank ) - 1
       offset(2) =  get_layout_2D_j_min( my_layout, my_rank ) - 1

       if (itime == 1) then

          SLL_ALLOCATE(x1(local_nx1,local_nx2),error)
          SLL_ALLOCATE(x2(local_nx1,local_nx2),error)
          ! These are rather values in the eta space than in x... should change
          x1_min = sim%mesh2d_x%eta1_min
          x1_max = sim%mesh2d_x%eta1_max
          x2_min = sim%mesh2d_x%eta2_min
          x2_max = sim%mesh2d_x%eta2_max
          x3_min = sim%mesh2d_v%eta1_min
          x3_max = sim%mesh2d_v%eta1_max
          x4_min = sim%mesh2d_v%eta1_min
          x4_max = sim%mesh2d_v%eta1_max
   
          delta_x1 = sim%mesh2d_x%delta_eta1
          delta_x2 = sim%mesh2d_x%delta_eta2
          delta_x3 = sim%mesh2d_v%delta_eta1
          delta_x4 = sim%mesh2d_v%delta_eta2
   
          do j = 1, local_nx2
             do i = 1, local_nx1
                global_indices =  local_to_global_2D( my_layout, (/i, j/) )
                gi = global_indices(1)
                gj = global_indices(2)
                x1(i,j) = x1_min + (gi-1._f64)*delta_x1
                x2(i,j) = x2_min + (gj-1._f64)*delta_x2
             end do
          end do
       
          call sll_hdf5_file_create("mesh_x"//c_layout//"_seq.h5",hdf_file_id,error)
          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,x1,"x1",error)
          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,x2,"x2",error)
          call sll_hdf5_file_close(hdf_file_id,error)

          deallocate(x1)
          deallocate(x2)

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
#endif
!!$  subroutine compute_electric_field_energy( efield, layout )
!!$    sll_cmplx64, dimension(:,:), intent(in) :: efield
!!$
!!$  end subroutine compute_electric_field_energy



end module sll_simulation_4d_vlasov_poisson_general
