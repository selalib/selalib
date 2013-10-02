! The idea of this simulation is to merge the functionalities of the qns-based
! simulation and the older, cartesian 4d simulation for the purposes of 
! debugging/understanding the behavior of the QNS one. Once this objective is
! fulfilled, this simulation can be deleted.

module sll_simulation_4d_qns_mixed_module

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
  use sll_general_coordinate_qn_solver_module
  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative
  implicit none

#define PRINT_PLOTS 1
  type, extends(sll_simulation_base_class) :: sll_simulation_4d_qns_mixed
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
     ! for QNS spline_degre in each direction
     sll_int32  :: spline_degree_eta1
     sll_int32  :: spline_degree_eta2
     ! for QNS boundary conditions
     sll_int32  :: bc_left
     sll_int32  :: bc_right
     sll_int32  :: bc_bottom
     sll_int32  :: bc_top
     ! the logical meshes are split in two one for space, one for velocity
     type(sll_logical_mesh_2d), pointer    :: mesh2d_x
     type(sll_logical_mesh_2d), pointer    :: mesh2d_v
     ! This simulation only applies a coordinate transformation to the spatial
     ! coordinates.
     class(sll_coordinate_transformation_2d_base), pointer :: transfx
     type(poisson_2d_periodic_plan_cartesian_par), pointer :: poisson_plan
     type(general_coordinate_qn_solver), pointer           :: qns

     ! distribution functions. There are several because each array represents
     ! a differently shaped chunk of memory. In this example, each chunk 
     ! allows sequential operations in one given direction. f_x1x2 should 
     ! permit to carry out sequential operations in x1 and x2 for ex.
     sll_real64, dimension(:,:,:,:), pointer     :: f_x1x2_diff
     sll_real64, dimension(:,:,:,:), pointer     :: f_x1x2_c ! cartesian case 
     sll_real64, dimension(:,:,:,:), pointer     :: f_x3x4_c ! cartesian case
     sll_real64, dimension(:,:,:,:), pointer     :: f_x1x2_q ! qns case
     sll_real64, dimension(:,:,:,:), pointer     :: f_x3x4_q ! qns case
     sll_real64, dimension(:,:,:), allocatable   :: partial_reduction_c
     sll_real64, dimension(:,:,:), allocatable   :: partial_reduction_q
     sll_real64, dimension(:,:), allocatable     :: rho_full 
     sll_real64, dimension(:,:), allocatable     :: rho_x2 
     sll_real64, dimension(:,:), allocatable     :: rho_diff
     ! for cartesian solver...
     sll_real64, dimension(:,:), allocatable     :: rho_x1 
     sll_real64, dimension(:,:), allocatable     :: rho_split_c
     sll_real64, dimension(:,:), allocatable     :: rho_split_q
     sll_real64, dimension(:,:), allocatable     :: phi_x1
     sll_real64, dimension(:,:), allocatable     :: phi_x2
     sll_real64, dimension(:,:), allocatable     :: phi_split
     sll_real64, dimension(:,:), allocatable     :: phi_diff
     sll_real64, dimension(:,:), allocatable     :: phi_q

     ! for remap
     type(layout_4D), pointer :: sequential_x1x2
     type(layout_4D), pointer :: sequential_x3x4
     type(layout_2D), pointer :: rho_full_layout
     type(layout_2D), pointer :: rho_seq_x1
     type(layout_2D), pointer :: rho_seq_x2
     type(layout_2D), pointer :: split_rho_layout ! not sequential in any dir.
     type(layout_2D), pointer :: split_phi_layout ! not sequential in any dir.
     type(remap_plan_2D_real64), pointer :: split_to_full

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
     ! interpolation any arbitrary spline
      type(arb_deg_2d_interpolator)     :: interp_rho
      type(arb_deg_2d_interpolator)     :: interp_phi
     ! Field accumulator
     sll_comp64, dimension(:,:), allocatable :: efield_x1
     sll_comp64, dimension(:,:), allocatable :: efield_x2
     sll_comp64, dimension(:,:), allocatable :: efield_q
     sll_real64, dimension(:,:), allocatable :: efield_exact
     sll_comp64, dimension(:,:), allocatable :: efield_split
     sll_real64, dimension(:,:), allocatable :: efield_x1_diff
     sll_real64, dimension(:,:), allocatable :: efield_x2_diff
     sll_real64, dimension(:,:), allocatable :: efield_x1_ref
     sll_real64, dimension(:,:), allocatable :: efield_x1_qns
     sll_real64, dimension(:,:), allocatable :: efield_x2_qns
     sll_real64, dimension(:,:), allocatable :: efield_x1_cart
     sll_real64, dimension(:,:), allocatable :: efield_x2_cart
     sll_real64, dimension(:,:), allocatable :: efield_x1_ref_diff_qns
     sll_real64, dimension(:,:), allocatable :: efield_x1_ref_diff_cart
     sll_real64, dimension(:,:), allocatable :: phi_ref
     ! for distribution function initializer:
     procedure(sll_scalar_initializer_4d), nopass, pointer :: init_func
     sll_real64, dimension(:), pointer :: params
     ! for general coordinate QNS, analytical fields
     procedure(two_var_parametrizable_function),nopass,pointer :: a11_f
     procedure(two_var_parametrizable_function),nopass,pointer :: a12_f
     procedure(two_var_parametrizable_function),nopass,pointer :: a21_f
     procedure(two_var_parametrizable_function),nopass,pointer :: a22_f
     procedure(two_var_parametrizable_function),nopass,pointer :: c_f
   contains
     procedure, pass(sim) :: run => run_4d_qns_mixed
     procedure, pass(sim) :: init_from_file => init_4d_qns_mixed
  end type sll_simulation_4d_qns_mixed

  interface delete
     module procedure delete_4d_qns_mixed
  end interface delete

  interface initialize
     module procedure initialize_4d_qns_mixed
  end interface initialize

contains

  ! Tentative function to initialize the simulation object 'manually'.
  subroutine initialize_4d_qns_mixed( &
   sim, &
   mesh2d_x, &
   mesh2d_v, &
   transformation_x, &
   init_func, &
   params,&
   a11_f,&
   a12_f,&
   a21_f,&
   a22_f,&
   c_f,&
   spline_degre1,&
   spline_degre2,&
   bc_left,&
   bc_right,&
   bc_bottom,&
   bc_top)

   type(sll_simulation_4d_qns_mixed), intent(inout)     :: sim
   type(sll_logical_mesh_2d), pointer                    :: mesh2d_x
   type(sll_logical_mesh_2d), pointer                    :: mesh2d_v
   class(sll_coordinate_transformation_2d_base), pointer :: transformation_x
   procedure(sll_scalar_initializer_4d)                  :: init_func
   sll_real64, dimension(:), target                      :: params
   procedure(two_var_parametrizable_function) :: a11_f
   procedure(two_var_parametrizable_function) :: a12_f
   procedure(two_var_parametrizable_function) :: a21_f
   procedure(two_var_parametrizable_function) :: a22_f
   procedure(two_var_parametrizable_function) :: c_f
   sll_int32  :: spline_degre1
   sll_int32  :: spline_degre2
   sll_int32  :: bc_left
   sll_int32  :: bc_right
   sll_int32  :: bc_bottom
   sll_int32  :: bc_top

   sim%mesh2d_x  => mesh2d_x
   sim%mesh2d_v  => mesh2d_v
   sim%transfx   => transformation_x
   sim%init_func => init_func
   sim%params    => params
   sim%a11_f     => a11_f
   sim%a12_f     => a12_f
   sim%a21_f     => a21_f
   sim%a22_f     => a22_f
   sim%c_f       => c_f
   sim%spline_degree_eta1 = spline_degre1
   sim%spline_degree_eta2 = spline_degre2

   sim%bc_left   = bc_left
   sim%bc_right  = bc_right
   sim%bc_bottom = bc_bottom
   sim%bc_top    = bc_top

   call sim%interp_phi%initialize( &
        sim%mesh2d_x%num_cells1 +1, &
        sim%mesh2d_x%num_cells2 +1, &
        sim%mesh2d_x%eta1_min, &
        sim%mesh2d_x%eta1_max, &
        sim%mesh2d_x%eta2_min, &
        sim%mesh2d_x%eta2_max, &
        sim%bc_left, &
        sim%bc_right, &
        sim%bc_bottom, &
        sim%bc_top, &
        sim%spline_degree_eta1, &
        sim%spline_degree_eta2)

   call sim%interp_rho%initialize( &
        sim%mesh2d_x%num_cells1 +1, &
        sim%mesh2d_x%num_cells2 +1, &
        sim%mesh2d_x%eta1_min, &
        sim%mesh2d_x%eta1_max, &
        sim%mesh2d_x%eta2_min, &
        sim%mesh2d_x%eta2_max, &
        sim%bc_left, &
        sim%bc_right, &
        sim%bc_bottom, &
        sim%bc_top, &
        sim%spline_degree_eta1, &
        sim%spline_degree_eta2)
 end subroutine initialize_4d_qns_mixed


 subroutine init_4d_qns_mixed( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_qns_mixed), intent(inout) :: sim
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
  end subroutine init_4d_qns_mixed

  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_4d_qns_mixed(sim)
    class(sll_simulation_4d_qns_mixed), intent(inout) :: sim
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
    sll_real64 :: efield_energy_total_q
    sll_real64 :: efield_energy_total_c
    sll_real64 :: efield_energy_total_t
    ! The following could probably be abstracted for convenience
#define BUFFER_SIZE 50
    sll_real64, dimension(BUFFER_SIZE) :: buffer_c
    sll_real64, dimension(BUFFER_SIZE) :: buffer_q
    sll_real64, dimension(BUFFER_SIZE) :: buffer_t
    sll_real64, dimension(BUFFER_SIZE) :: buffer_result_c
    sll_real64, dimension(BUFFER_SIZE) :: buffer_result_q
    sll_real64, dimension(BUFFER_SIZE) :: buffer_result_t
    sll_real64, dimension(BUFFER_SIZE) :: num_particles_local_c
    sll_real64, dimension(BUFFER_SIZE) :: num_particles_local_q
    sll_real64, dimension(BUFFER_SIZE) :: num_particles_global_c
    sll_real64, dimension(BUFFER_SIZE) :: num_particles_global_q
    sll_real64 :: tmp
    sll_real64, dimension(1) :: tmp1
    sll_int32 :: buffer_counter
    sll_int32 :: efield_energy_file_id_cart
    sll_int32 :: droite_test_pente
    sll_int32 :: efield_energy_file_id_qns
    sll_int32 :: efield_energy_file_id_teo
    sll_int32 :: num_particles_file_id_cart
    sll_int32 :: num_particles_file_id_qns
    sll_int32 :: global_indices(4)
    sll_int32 :: iplot
    character(len=4) :: cplot
    type(sll_scalar_field_2d_base_ptr)                    :: a11_field_mat
    type(sll_scalar_field_2d_base_ptr)                    :: a21_field_mat
    type(sll_scalar_field_2d_base_ptr)                    :: a12_field_mat
    type(sll_scalar_field_2d_base_ptr)                    :: a22_field_mat
    class(sll_scalar_field_2d_base), pointer              :: c_field
    class(sll_scalar_field_2d_discrete_alt), pointer      :: rho
    type(sll_scalar_field_2d_discrete_alt), pointer       :: phi
    sll_real64, dimension(:), allocatable :: send_buf
    sll_real64, dimension(:), allocatable :: recv_buf
    sll_int32, dimension(:), allocatable  :: recv_sz
    sll_real64, dimension(:,:), allocatable :: phi_values
    sll_real64 :: density_tot
    sll_int32  :: send_size   ! for allgatherv operation
    sll_int32, dimension(:), allocatable :: disps ! for allgatherv operation
    ! only for debugging...
!!$    sll_real64, dimension(:,:), allocatable :: ex_field
!!$    sll_real64, dimension(:,:), allocatable :: ey_field
    
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

    ! continue with the fields
    a11_field_mat%base => new_scalar_field_2d_analytic_alt( &
         sim%a11_f, &
         "a11", &
         sim%transfx, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top) 

    a12_field_mat%base => new_scalar_field_2d_analytic_alt( &
         sim%a12_f, &
         "a12", &
         sim%transfx, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top) 

    a21_field_mat%base => new_scalar_field_2d_analytic_alt( &
         sim%a21_f, &
         "a21", &
         sim%transfx, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top)
    
    a22_field_mat%base => new_scalar_field_2d_analytic_alt( &
         sim%a22_f, &
         "a22", &
         sim%transfx, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top) 

    print*, 'finished fields a'
    c_field => new_scalar_field_2d_analytic_alt( &
         sim%c_f, &
         "c_field", &
         sim%transfx, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top)

    print*, 'finished fields c'
    SLL_ALLOCATE(phi_values(nc_x1+1,nc_x2+1),ierr)
    phi_values(:,:) = 0.0_f64

    phi => new_scalar_field_2d_discrete_alt( &
         phi_values, &
         "phi_potential", &
         sim%interp_phi, &
         sim%transfx, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top)

    print*, 'finished fields phi'
    buffer_counter = 1

    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)

    SLL_ALLOCATE(recv_sz(sim%world_size),ierr)
    SLL_ALLOCATE(disps(sim%world_size),ierr)

    if( sim%my_rank == 0 ) then
       call sll_new_file_id( efield_energy_file_id_cart, ierr )
       if( ierr == 1 ) then
          print *, 'sll_new_file_id() failed to obtain a file identifier.', &
               ' Exiting...'
          stop
       end if
    end if

    if( sim%my_rank == 0 ) then
       call sll_new_file_id( efield_energy_file_id_qns, ierr )
       if( ierr == 1 ) then
          print *, 'sll_new_file_id() failed to obtain a file identifier.', &
               ' Exiting...'
          stop
       end if
    end if

    if( sim%my_rank == 0 ) then
       call sll_new_file_id( efield_energy_file_id_teo, ierr )
       if( ierr == 1 ) then
          print *, 'sll_new_file_id() failed to obtain a file identifier.', &
               ' Exiting...'
          stop
       end if
    end if

    if( sim%my_rank == 0 ) then
       call sll_new_file_id( num_particles_file_id_cart, ierr )
       if( ierr == 1 ) then
          print *, 'sll_new_file_id() failed to obtain a file identifier.', &
               ' Exiting...'
          stop
       end if
    end if

    if( sim%my_rank == 0 ) then
       call sll_new_file_id( num_particles_file_id_qns, ierr )
       if( ierr == 1 ) then
          print *, 'sll_new_file_id() failed to obtain a file identifier.', &
               ' Exiting...'
          stop
       end if
    end if

    ! allocate the layouts...
    sim%sequential_x1x2  => new_layout_4D( sll_world_collective )
    sim%sequential_x3x4  => new_layout_4D( sll_world_collective )
    sim%rho_full_layout  => new_layout_2D( sll_world_collective )
    sim%rho_seq_x1       => new_layout_2D( sll_world_collective )
    sim%rho_seq_x2       => new_layout_2D( sll_world_collective )
    sim%split_rho_layout => new_layout_2D( sll_world_collective )
    sim%split_phi_layout => new_layout_2D( sll_world_collective )

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

    print *, 'sequential_x3x4 mode... ', &
         'initializing layouts and allocating tables'
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
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%split_phi_layout )


    ! move this down!!! This is self-contained but we are initializing
    ! entities split in x1 and x2, so this isn't the logical place for this.
    call initialize_layout_with_distributed_2D_array( &
         nc_x1+1, &
         nc_x2+1, &
         1, &
         sim%world_size, &
         sim%rho_seq_x1 )
    
    call compute_local_sizes_2d( sim%rho_seq_x1, loc_sz_x1, loc_sz_x2 )
    SLL_ALLOCATE(sim%rho_x1(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%phi_x1(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%phi_diff(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%rho_diff(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%phi_q(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x1(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_q(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_exact(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x1_ref(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%phi_ref(loc_sz_x1,loc_sz_x2),ierr)
    ! end move this down

    call compute_local_sizes( sim%split_rho_layout, loc_sz_x1, loc_sz_x2)
    ! This layout is also useful to represent the charge density array. Since
    ! this is a result of a local reduction on x3 and x4, the new layout is
    ! 2D but with the same dimensions of the process mesh in x1 and x2.
    SLL_ALLOCATE(sim%rho_split_c(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%rho_split_q(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%phi_split(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(send_buf(loc_sz_x1*loc_sz_x2), ierr)
    SLL_ALLOCATE(sim%efield_split(loc_sz_x1,loc_sz_x2),ierr)
    ! the following only works in sequential mode! no parallel!
    SLL_ALLOCATE(sim%efield_x1_diff(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x2_diff(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x1_qns(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x2_qns(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x1_cart(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x2_cart(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x1_ref_diff_qns(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x1_ref_diff_cart(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%rho_full(nc_x1+1,nc_x2+1),ierr) ! changed from nc,nc
    SLL_ALLOCATE(recv_buf((nc_x1+1)*(nc_x2+1)),ierr) !changed from nc, nc

    ! move this down
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

    !end move down

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
    
    print *, 'nproc_x1: ', sim%nproc_x1
    print *, 'nproc_x2: ', sim%nproc_x2
    print *, 'nproc_x3: ', sim%nproc_x3
    print *, 'nproc_x4: ', sim%nproc_x4
    
    call compute_local_sizes_4d( sim%sequential_x3x4, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )
    
    SLL_ALLOCATE(sim%f_x3x4_c(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)
    SLL_ALLOCATE(sim%f_x3x4_q(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)
    ! These dimensions are also the ones needed for the array where we store
    ! the intermediate results of the charge density computation.
    SLL_ALLOCATE(sim%partial_reduction_c(loc_sz_x1,loc_sz_x2, loc_sz_x3),ierr)
    SLL_ALLOCATE(sim%partial_reduction_q(loc_sz_x1,loc_sz_x2, loc_sz_x3),ierr)
    
    
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
    
    print *, 'sequential x1x2 mode...'
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
    ! are out-of-place, we will allocate two different arrays, one for each
    ! layout.
    call compute_local_sizes_4d( sim%sequential_x1x2, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )
    SLL_ALLOCATE(sim%f_x1x2_c(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)
    SLL_ALLOCATE(sim%f_x1x2_q(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)
    SLL_ALLOCATE(sim%f_x1x2_diff(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)

    call compute_local_sizes_4d( sim%sequential_x3x4, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )

    print *, 'Initialize distributed array that represents f: the distr. func.'
    call sll_4d_parallel_array_initializer( &
         sim%sequential_x3x4, &
         sim%mesh2d_x, &
         sim%mesh2d_v, &
         sim%f_x3x4_c, &
         sim%init_func, &
         sim%params, &
         transf_x1_x2=sim%transfx )

    call sll_4d_parallel_array_initializer( &
         sim%sequential_x3x4, &
         sim%mesh2d_x, &
         sim%mesh2d_v, &
         sim%f_x3x4_q, &
         sim%init_func, &
         sim%params, &
         transf_x1_x2=sim%transfx )


    sim%rho_split_c(:,:) = 0.0_f64
    sim%rho_split_q(:,:) = 0.0_f64
    ! this only works because there is no transformation applied in the
    ! velocity space...
    call compute_charge_density( &
         sim%mesh2d_x,           &
         sim%mesh2d_v,           &
         size(sim%f_x3x4_c,1),     &
         size(sim%f_x3x4_c,2),     &
         sim%f_x3x4_c,             &
         sim%partial_reduction_c,  &
         sim%rho_split_c )

    call compute_charge_density( &
         sim%mesh2d_x,           &
         sim%mesh2d_v,           &
         size(sim%f_x3x4_q,1),     &
         size(sim%f_x3x4_q,2),     &
         sim%f_x3x4_q,             &
         sim%partial_reduction_q,  &
         sim%rho_split_q )
 
   global_indices(1:2) =  &
         local_to_global_2D( sim%split_rho_layout, (/1, 1/) )

   ! rho_split will be used by both calculations, cartesian and qns, so
   ! at this point there is no divergence in these calculations... 
   call sll_gnuplot_rect_2d_parallel( &
        sim%mesh2d_x%eta1_min, &
        sim%mesh2d_x%delta_eta1, &
        sim%mesh2d_x%eta2_min, &
        sim%mesh2d_x%delta_eta2, &
        sim%rho_split_c, &
        "rho_split_cartesian", &
        0, &
        ierr )
 
   call sll_gnuplot_rect_2d_parallel( &
        sim%mesh2d_x%eta1_min, &
        sim%mesh2d_x%delta_eta1, &
        sim%mesh2d_x%eta2_min, &
        sim%mesh2d_x%delta_eta2, &
        sim%rho_split_q, &
        "rho_split_qns", &
        0, &
        ierr )

    call sll_gnuplot_rect_2d_parallel( &
        sim%mesh2d_x%eta1_min, &
        sim%mesh2d_x%delta_eta1, &
        sim%mesh2d_x%eta2_min, &
        sim%mesh2d_x%delta_eta2, &
        sim%rho_split_q-sim%rho_split_c, &
        "rho_split_diff", &
        0, &
        ierr )


    ! for the poisson cartesian case:
    sim%split_to_seqx1 => &
         NEW_REMAP_PLAN(sim%split_rho_layout, sim%rho_seq_x1, sim%rho_split_c)
    call apply_remap_2D( sim%split_to_seqx1, sim%rho_split_c, sim%rho_x1 )

    global_indices(1:2) =  &
         local_to_global_2D( sim%rho_seq_x1, (/1, 1/) )
       
    ! rho_x1 is the input for the cartesian poisson solver
    call sll_gnuplot_rect_2d_parallel( &
         sim%mesh2d_x%eta1_min, &
         sim%mesh2d_x%delta_eta1, &
         sim%mesh2d_x%eta2_min, &
         sim%mesh2d_x%delta_eta2, &
         sim%rho_x1, &
         "rho_x1", &
         0, &
         ierr )


    ! At the same time, the sequential QNS needs that every processor
    ! contain the full information on rho, which we proceed to provide next.
    
    call load_buffer( sim%split_rho_layout, sim%rho_split_q, send_buf )

    recv_sz(:) = receive_counts_array( sim%split_rho_layout, sim%world_size )

    
    send_size = size(send_buf)

    call compute_displacements_array( &
         sim%split_rho_layout, &
         sim%world_size, &
         disps )

    call sll_collective_allgatherv_real64( &
         sll_world_collective, &
         send_buf, &
         send_size, &
         recv_sz, &
         disps, &
         recv_buf )

    call unload_buffer(sim%split_rho_layout, recv_buf, sim%rho_full)
    
    call sll_gnuplot_rect_2d_parallel( &
         sim%mesh2d_x%eta1_min, &
         sim%mesh2d_x%delta_eta1, &
         sim%mesh2d_x%eta2_min, &
         sim%mesh2d_x%delta_eta2, &
         sim%rho_full, &
         "rho_full", &
         0, &
         ierr )
    
    call compute_average_f( &
         sim,&
         sim%mesh2d_x,&
         sim%rho_full, &
         density_tot )
    

    
    rho => new_scalar_field_2d_discrete_alt( &
         sim%rho_full - density_tot, &
         "rho_field_qns", &
         sim%interp_rho, &     
         sim%transfx, &
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top)


    
    if(sim%my_rank == 0) then
       call rho%write_to_file(0)
    end if


       do j=1,nc_x2+1 ! pay attention to the +1 ...
          do i=1,nc_x1+1
             sim%rho_diff(i,j) = sim%rho_x1(i,j) -  rho%value_at_indices(i,j)-density_tot!sim%rho_full(i,j)
             !print*, sim%rho_diff(i,j)
          end do
       end do

       call sll_gnuplot_rect_2d_parallel( &
         sim%mesh2d_x%eta1_min, &
         sim%mesh2d_x%delta_eta1, &
         sim%mesh2d_x%eta2_min, &
         sim%mesh2d_x%delta_eta2, &
         sim%rho_diff, &
         "rho_diff", &
         0, &
         ierr )



    ! It should not matter that these remap plans are defined in reference
    ! to the 'cartesian' case. We should be able to use them in the qns data
    ! as well.
    sim%seqx3x4_to_seqx1x2 => &
         NEW_REMAP_PLAN(sim%sequential_x3x4,sim%sequential_x1x2,sim%f_x3x4_c)
    
    sim%seqx1x2_to_seqx3x4 => &
         NEW_REMAP_PLAN(sim%sequential_x1x2,sim%sequential_x3x4,sim%f_x1x2_c)
    
    call apply_remap_4D( sim%seqx3x4_to_seqx1x2, sim%f_x3x4_q, sim%f_x1x2_q )
    call apply_remap_4D( sim%seqx3x4_to_seqx1x2, sim%f_x3x4_c, sim%f_x1x2_c )
    
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
    
    call advection_x1x2_q(sim,0.5*sim%dt)
    call advection_x1x2_c(sim,0.5*sim%dt)
    

!!$    call sll_gnuplot_rect_2d_parallel( &
!!$         sim%mesh2d_x%eta1_min+1*sim%mesh2d_x%delta_eta1, &
!!$         sim%mesh2d_x%delta_eta1, &
!!$         sim%mesh2d_x%eta2_min+1*sim%mesh2d_x%delta_eta2, &
!!$         sim%mesh2d_x%delta_eta2, &
!!$         sim%f_x1x2_diff(:,:,1,1), &
!!$         "f_difference", &
!!$         0, &
!!$         ierr )
!!$
!!$    print *, 'sum of f_difference = ', sum(sim%f_x1x2_diff(:,:,:,:))
    ! Conclusion number 1: the space advection is not guilty!
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

    ! qns stuff
    sim%qns => new_general_qn_solver( &
         sim%spline_degree_eta1, & 
         sim%spline_degree_eta2, & 
         sim%mesh2d_x%num_cells1, &
         sim%mesh2d_x%num_cells2, &
         QNS_GAUSS_LEGENDRE, &  
         QNS_GAUSS_LEGENDRE, &  
         sim%bc_left, &
         sim%bc_right, &
         sim%bc_bottom, &
         sim%bc_top, &
         sim%mesh2d_x%eta1_min, &  
         sim%mesh2d_x%eta1_max, & 
         sim%mesh2d_x%eta2_min, & 
         sim%mesh2d_x%eta2_max ) 

    call factorize_mat_qns(&
         sim%qns, & 
         a11_field_mat, &
         a12_field_mat, &
         a21_field_mat, &
         a22_field_mat, &
         c_field, &
         rho)

    do j=1,nc_x2+1
       do i=1,nc_x1+1
          eta1 = eta1_min + (i-1)*delta1
          eta2 = eta2_min + (j-1)*delta2
          sim%efield_x1_ref(i,j) = - 0.05_f64/(0.5)*sin(0.5*eta1)
          sim%phi_ref(i,j)       = - 0.05_f64/((0.5)**2)*cos(0.5*eta1)
       end do
    end do

    call sll_gnuplot_rect_2d_parallel( &
         sim%mesh2d_x%eta1_min, &
         sim%mesh2d_x%delta_eta1, &
         sim%mesh2d_x%eta2_min, &
         sim%mesh2d_x%delta_eta2, &
         sim%efield_x1_ref, &
         "electric_field_reference", &
         0, &
         ierr )


    print*, ' ... finished initialization, entering main loop.'


    
    ! ------------------------------------------------------------------------
    !
    !                                MAIN LOOP
    !
    ! ------------------------------------------------------------------------
!#if 0    
    do itime=1,sim%num_iterations
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

       ! compute exact electric field (theoretical)

       do j=1,nc_x2
          do i=1,nc_x1
             eta1 = eta1_min + (i-1)*delta1
             eta2 = eta2_min + (j-1)*delta2
             tmp = efield_exact(sim%dt*(real(itime-1) + 0.5_f64),&
                  0.5_f64, &
                  0.05_f64, &
                  eta1, &
                  eta2 )
             sim%efield_exact(i,j)=tmp
 !            print *, 'damned energy is: ', tmp
          end do
       end do

       call apply_remap_4D( sim%seqx1x2_to_seqx3x4, sim%f_x1x2_q, sim%f_x3x4_q )
       call apply_remap_4D( sim%seqx1x2_to_seqx3x4, sim%f_x1x2_c, sim%f_x3x4_c )
       
       
       call compute_local_sizes_4d( &
            sim%sequential_x1x2, &
            loc_sz_x1,           &
            loc_sz_x2,           &
            loc_sz_x3,           &
            loc_sz_x4 )
       
       sim%rho_split_c(:,:) = 0.0_f64
       sim%rho_split_q(:,:) = 0.0_f64
       
       call compute_charge_density( &
            sim%mesh2d_x,           &
            sim%mesh2d_v,           &
            size(sim%f_x3x4_q,1),     &
            size(sim%f_x3x4_q,2),     &
            sim%f_x3x4_q,             &
            sim%partial_reduction_q,  &
            sim%rho_split_q )
       
       call compute_charge_density( &
            sim%mesh2d_x,           &
            sim%mesh2d_v,           &
            size(sim%f_x3x4_c,1),     &
            size(sim%f_x3x4_c,2),     &
            sim%f_x3x4_c,             &
            sim%partial_reduction_c,  &
            sim%rho_split_c )
       
       
       
       global_indices(1:2) =  &
            local_to_global_2D( sim%split_rho_layout, (/1, 1/) )
#if PRINT_PLOTS       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%rho_split_c, &
            "rho_split_cartesian", &
            itime, &
            ierr )
       
       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%rho_split_q-sim%rho_split_c, &
            "rho_split_qns", &
            itime, &
            ierr )
#endif
       call load_buffer( sim%split_rho_layout, sim%rho_split_q, send_buf )
       
       recv_sz(:) = receive_counts_array( sim%split_rho_layout, sim%world_size )
       
       call sll_collective_allgatherv( &
            sll_world_collective, &
            send_buf, &
            size(send_buf), &
            recv_sz, &
            disps, &
            recv_buf )
       
       call unload_buffer(sim%split_rho_layout, recv_buf, sim%rho_full)
       
       if(sim%my_rank == 0) then
          call sll_gnuplot_rect_2d_parallel( &
               sim%mesh2d_x%eta1_min, &
               sim%mesh2d_x%delta_eta1, &
               sim%mesh2d_x%eta2_min, &
               sim%mesh2d_x%delta_eta2, &
               sim%rho_full, &
               "rho_full_check", &
               itime, &
               ierr )
       end if
       
       ! for poisson cartesian:
       ! Re-arrange rho_split in a way that permits sequential operations in 
       ! x1, to feed to the Poisson solver.
       call apply_remap_2D( sim%split_to_seqx1, sim%rho_split_c, sim%rho_x1 )
       

       
       ! the rho field has a pointer to sim%rho_full so it is already 
       ! 'aware' that the data has changed. However, the interpolation
       ! coefficients are out of date.
       !
       ! It should be seriously considered to modify the interpolators
       ! in such a way that they always carry a pointer to the data to 
       ! be interpolated. In such case we would need:
       ! - a routine to (re-)set this pointer to the data
       ! - a change of the compute_coefficients interface to set the
       !   data parameter as optional...
       
       call compute_average_f( &
            sim,&
            sim%mesh2d_x,&
            sim%rho_full, &
            density_tot )
       
       !ATTENTION TO CHANGES HERE
       print*, 'density', density_tot
       call rho%update_interpolation_coefficients(sim%rho_full-density_tot)
       
       do j=1,nc_x2+1 ! pay attention to the +1 ...
          do i=1,nc_x1+1
             sim%rho_diff(i,j) = sim%rho_x1(i,j) - rho%value_at_indices(i,j)-density_tot! sim%rho_full(i,j)
             ! print*, rho%value_at_indices(i,j)
          end do
       end do
       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%rho_diff, &
            "rho_diff", &
            itime, &
            ierr )
       
       
       
       
!!$       if(sim%my_rank == 0) then
!!$          call rho%write_to_file(itime)
!!$       end if
       
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
       
!!$       global_indices(1:2) =  &
!!$            local_to_global_2D( sim%rho_seq_x1, (/1, 1/) )
!!$       
!!$       call sll_gnuplot_rect_2d_parallel( &
!!$          sim%mesh2d_x%eta1_min+(global_indices(1)-1)*sim%mesh2d_x%delta_eta1, &
!!$          sim%mesh2d_x%delta_eta1, &
!!$          sim%mesh2d_x%eta2_min+(global_indices(2)-1)*sim%mesh2d_x%delta_eta2, &
!!$          sim%mesh2d_x%delta_eta2, &
!!$          sim%rho_x1, &
!!$          "rho_x1", &
!!$          itime, &
!!$          ierr )
       !       if(sim%my_rank == 0) call rho%write_to_file(itime)
       

       call solve_qns( &
            sim%qns, &
            rho, &
            phi )
!!$       call solve_quasi_neutral_eq_general_coords( &
!!$            sim%qns, & 
!!$            a11_field_mat, &
!!$            a12_field_mat, &
!!$            a21_field_mat, &
!!$            a22_field_mat, &
!!$            c_field, &
!!$            rho, &
!!$            phi )
       
       if(sim%my_rank == 0) then
          call phi%write_to_file(itime)
       end if
       
       do j=1,nc_x2+1
          do i=1,nc_x1+1
             sim%phi_q(i,j) = phi%value_at_indices(i,j)
          end do
       end do
       ! print *, 'the sum of phi  = ', sum(sim%phi_q(1:nc_x1,1:nc_x2))
       
       
#if PRINT_PLOTS
       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%phi_q, &
            "phi_q", &
            itime, &
            ierr )
       
       
       do j=1,nc_x2+1
          do i=1,nc_x1+1
             sim%phi_diff(i,j) =  phi%value_at_indices(i,j)
             
          end do
       end do
       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%phi_diff, &
            "phi_diff_qns", &
            itime, &
            ierr )
#endif

       call solve_poisson_2d_periodic_cartesian_par( &
            sim%poisson_plan, &
            sim%rho_x1, &
            sim%phi_x1)
!!$       print *, 'What the heck is going on with phi?'
!!$       print *,'last column = '
!!$       print *, sim%phi_x1(:,1) - sim%phi_x1(:,33)
       
       global_indices(1:2) =  &
            local_to_global_2D( sim%rho_seq_x1, (/1, 1/) )
       
#if PRINT_PLOTS       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%phi_x1, &
            "phi_x1", &
            itime, &
            ierr )
       
       
       do j=1,nc_x2+1
          do i=1,nc_x1+1
             sim%phi_diff(i,j) =  sim%phi_x1(i,j)- phi%value_at_indices(i,j)
             !print*, 'value diff', sim%phi_diff(i,j)
          end do
       end do
       
       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%phi_diff, &
            "phi_diff_cart", &
            itime, &
            ierr )
#endif
       
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
       
       call compute_electric_field_eta1( &
            sim%phi_q, &
            loc_sz_x1, &
            loc_sz_x2, &
            sim%mesh2d_x%delta_eta1, &
            sim%efield_q )
       
       !print*,'tab',real(sim%efield_q(:,1),f64) 
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
       
       call compute_electric_field_eta2( &
            sim%phi_q, &
            loc_sz_x1, &
            loc_sz_x2, &
            sim%mesh2d_x%delta_eta2, &
            sim%efield_q )
       
       !       do j=1,loc_sz_x1
       !          do i=1
       
       !print *, 'sum of efield poisson = ', sum(real(sim%efield_x1(1:loc_sz_x1-1,1:loc_sz_x2-1)))
       !print *, 'sum of efield QNS     = ', sum(real(sim%efield_q(1:loc_sz_x1-1,1:loc_sz_x2-1)))
       
       !print*,'tab',real(sim%efield_q(:,1),f64) 
       ! But now, to make the electric field data configuration compatible with
       ! the sequential operations in x2x3 we need still another remap 
       ! operation.
       call apply_remap_2D( &
            sim%efld_seqx2_to_split, &
            sim%efield_x2, &
            sim%efield_split )
       
#if PRINT_PLOTS 
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             !             sim%efield_x1_diff(i,j) = sim%efield_x1_ref(i,j) - &
             !                 ( - phi%first_deriv_eta1_value_at_indices(i,j))
             sim%efield_x1_cart(i,j)= real( sim%efield_split(i,j),f64) 
             sim%efield_x2_cart(i,j)= aimag(sim%efield_split(i,j))
          end do
       end do
       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%efield_x1_cart, &
            "electric_field_x1_cart", &
            itime, &
            ierr )
       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%efield_x2_cart, &
            "electric_field_x2_cart", &
            itime, &
            ierr )
#endif

       
       call compute_local_sizes_4d( sim%sequential_x3x4, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 
       
       
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             sim%efield_x1_ref_diff_qns(i,j) = sim%efield_x1_ref(i,j) - &
                  ( - phi%first_deriv_eta1_value_at_indices(i,j))
             sim%efield_x1_qns(i,j)=-phi%first_deriv_eta1_value_at_indices(i,j)
             sim%efield_x2_qns(i,j)=-phi%first_deriv_eta2_value_at_indices(i,j)
             
          end do
       end do
       
  
#if PRINT_PLOTS
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%efield_x1_ref_diff_qns, &
            "electric_field_difference_ref_qns", &
            itime, &
            ierr )
       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%efield_x1_qns, &
            "electric_field_x1_qns", &
            itime, &
            ierr )
       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%efield_x2_qns, &
            "electric_field_x2_qns", &
            itime, &
            ierr )
       
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             sim%efield_x1_ref_diff_cart(i,j) =  real(sim%efield_x1_cart(i,j)) - sim%efield_exact(i,j)!sim%efield_x1_ref(i,j) - &
                 ! real(sim%efield_x1_cart(i,j))
             sim%efield_x1_qns(i,j)=-phi%first_deriv_eta1_value_at_indices(i,j)
             sim%efield_x2_qns(i,j)=-phi%first_deriv_eta2_value_at_indices(i,j)
          end do
       end do
       
!!$       print *, 'What the heck is going on with electric field x1?'
!!$       print *,'last column of reference = '
!!$       print *, sim%efield_x1_ref(:,33)
!!$       print *, 'last column of calculated efield x1 = '
!!$       print *, sim%efield_x1_cart(:,33)
       
       
       
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%efield_x1_ref_diff_cart(1:loc_sz_x1,1:loc_sz_x2), &
            "electric_field_difference_ref_cart", &
            itime, &
            ierr )
       
       do j=1,loc_sz_x2+1
          do i=1,loc_sz_x1
             sim%efield_x1_diff(i,j) = &!real(sim%efield_x1_cart(i,j),f64) - &
                  sim%efield_x1_qns(i,j)- sim%efield_exact(i,j)
             !             sim%efield_x2_diff(i,j) = aimag(sim%efield_split(i,j)) - &
             !                 ( - phi%first_deriv_eta2_value_at_indices(i,j))       
          end do
       end do
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2d_x%eta1_min, &
            sim%mesh2d_x%delta_eta1, &
            sim%mesh2d_x%eta2_min, &
            sim%mesh2d_x%delta_eta2, &
            sim%efield_x1_diff(1:loc_sz_x1,1:loc_sz_x2), &
            "electric_field_x1_difference", &
            itime, &
            ierr )
       
#endif
       
       efield_energy_total_c = 0.0_f64
       efield_energy_total_q = 0.0_f64
#if 1
       
       ! **********
!!! QNS CASE
       ! **********
       ! Start with dt in vx...(x3)
       do l=1,loc_sz_x4 !sim%mesh2d_v%num_cells2+1
          do j=1,loc_sz_x2
             do i=1,loc_sz_x1
                global_indices(1:2) = &
                     local_to_global_2D( sim%split_rho_layout, (/i,j/))
                eta1   =  eta1_min + real(global_indices(1)-1,f64)*delta1
                eta2   =  eta2_min + real(global_indices(2)-1,f64)*delta2
                !print*, phi%value_at_indices(i,j), 0.05/0.5**2*cos(0.5*(eta1))
                inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                jac_m  =  sim%transfx%jacobian_matrix(eta1,eta2)
!!$                ex     =  real( sim%efield_q(i,j),f64)
!!$                ey     =  aimag(sim%efield_q(i,j))
!!$                ex     =  real( sim%efield_split(i,j),f64)
!!$                ey     =  aimag(sim%efield_split(i,j))
                
!!$
                ex     =  - phi%first_deriv_eta1_value_at_point(eta1,eta2)
                ey     =  - phi%first_deriv_eta2_value_at_point(eta1,eta2)

                alpha3 = -sim%dt*(inv_j(1,1)*ex + inv_j(2,1)*ey)

                sim%f_x3x4_q(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                     nc_x3+1, &
                     sim%f_x3x4_q(i,j,:,l), &
                     alpha3 )

                
                efield_energy_total_q = efield_energy_total_q + &
                     delta1*delta2 *abs(jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1))&
                     *( (inv_j(1,1)*ex + inv_j(2,1)*ey)**2+ (inv_j(1,2)*ex + inv_j(2,2)*ey)**2)
                                  end do
          end do
       end do
       
       efield_energy_total_q = sqrt(efield_energy_total_q)
       !print*, 'nrj',efield_energy_total_q
       ! **********
!!! QNS CASE
       ! **********
       ! dt in vy...(x4)
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             do k=1,sim%mesh2d_v%num_cells1+1
                global_indices(1:2) = &
                     local_to_global_2D( sim%split_rho_layout, (/i,j/))
                eta1   =  eta1_min + real(global_indices(1)-1,f64)*delta1
                eta2   =  eta2_min + real(global_indices(2)-1,f64)*delta2
                inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
!!$                ex     =  real( sim%efield_q(i,j),f64)
!!$                ey     =  aimag(sim%efield_q(i,j))
!!$                ex     =  real( sim%efield_split(i,j),f64)
!!$                ey     =  aimag(sim%efield_split(i,j))
                
                ex     =  - phi%first_deriv_eta1_value_at_point(eta1,eta2)
                ey     =  - phi%first_deriv_eta2_value_at_point(eta1,eta2)
!!$                ex     =  - phi%first_deriv_eta1_value_at_indices(&
!!$                     global_indices(1), global_indices(2))
!!$                ey     =  - phi%first_deriv_eta2_value_at_indices(&
!!$                     global_indices(1), global_indices(2))
                alpha4 = -sim%dt*(inv_j(1,2)*ex + inv_j(2,2)*ey)

                !print*,eta1,eta2, inv_j(1,2),inv_j(2,2)
                sim%f_x3x4_q(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
                     nc_x4+1, &
                     sim%f_x3x4_q(i,j,k,:), &
                     alpha4 )
             end do
          end do
       end do
       ! **************
       ! cartesian case
       ! ************** 
       efield_energy_total_t = 0.0_f64

       do l=1,loc_sz_x4 !sim%mesh2d_v%num_cells2+1
          do j=1,loc_sz_x2
             do i=1,loc_sz_x1
                global_indices(1:2) = &
                     local_to_global_2D( sim%split_rho_layout, (/i,j/))
                eta1   =  eta1_min + real(global_indices(1)-1,f64)*delta1
                eta2   =  eta2_min + real(global_indices(2)-1,f64)*delta2

                inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                jac_m  =  sim%transfx%jacobian_matrix(eta1,eta2)
                ex     =  real( sim%efield_split(i,j),f64)
                ey     =  aimag(sim%efield_split(i,j))
                alpha3 = -sim%dt*(inv_j(1,1)*ex + inv_j(2,1)*ey)
                sim%f_x3x4_c(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                     nc_x3+1, &
                     sim%f_x3x4_c(i,j,:,l), &
                     alpha3 )
                efield_energy_total_c = efield_energy_total_c + &
                     delta1*delta2*(((jac_m(2,2)*ex - jac_m(2,1)*ey)**2 + &
                     (jac_m(1,1)*ex - jac_m(1,2)*ey)**2))&
                     /(jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1))

                efield_energy_total_t = efield_energy_total_t + &
                     delta1*delta2* sim%efield_exact(i,j)**2
             end do
          end do
       end do
       
       efield_energy_total_c = sqrt(efield_energy_total_c)
       
       efield_energy_total_t =  sqrt(efield_energy_total_t)
       
       
       ! **************
       ! cartesian case
       ! **************
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             do k=1,sim%mesh2d_v%num_cells1+1
                global_indices(1:2) = &
                     local_to_global_2D( sim%split_rho_layout, (/i,j/))
                eta1   =  eta1_min + real(global_indices(1)-1,f64)*delta1
                eta2   =  eta2_min + real(global_indices(2)-1,f64)*delta2
                inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                ex     =  real( sim%efield_split(i,j),f64)
                ey     =  aimag(sim%efield_split(i,j))
                alpha4 = -sim%dt*(inv_j(1,2)*ex + inv_j(2,2)*ey)
                sim%f_x3x4_c(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
                     nc_x4+1, &
                     sim%f_x3x4_c(i,j,k,:), &
                     alpha4 )
             end do
          end do
       end do
#endif
       call apply_remap_4D( sim%seqx3x4_to_seqx1x2, sim%f_x3x4_q, sim%f_x1x2_q )
       call apply_remap_4D( sim%seqx3x4_to_seqx1x2, sim%f_x3x4_c, sim%f_x1x2_c )
       
       
       call compute_local_sizes_4d( sim%sequential_x1x2, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 
       
       ! Approximate the integral of the distribution function along all
       ! directions.
       num_particles_local_q(buffer_counter) = &
            sum(sim%f_x3x4_q)*delta1*delta2*delta3*delta4
       if( buffer_counter == BUFFER_SIZE ) then
          call sll_collective_reduce_real64( &
               sll_world_collective, &
               num_particles_local_q, &
               BUFFER_SIZE, &
               MPI_SUM, &
               0, &
               num_particles_global_q )
          if(sim%my_rank == 0) then
             open(num_particles_file_id_qns,file="number_particles_qns",&
                  position="append")
             if(itime == BUFFER_SIZE) then
                rewind(num_particles_file_id_qns)
             end if
             do i=1,BUFFER_SIZE
                write(num_particles_file_id_qns,*) num_particles_global_q(i)
             end do
             close(num_particles_file_id_qns)
          end if
       end if
       
       num_particles_local_c(buffer_counter) = &
            sum(sim%f_x3x4_c)*delta1*delta2*delta3*delta4
       if( buffer_counter == BUFFER_SIZE ) then
          call sll_collective_reduce_real64( &
               sll_world_collective, &
               num_particles_local_c, &
               BUFFER_SIZE, &
               MPI_SUM, &
               0, &
               num_particles_global_c )
          if(sim%my_rank == 0) then
             open(num_particles_file_id_cart,file="number_particles_cart",&
                  position="append")
             if(itime == BUFFER_SIZE) then
                rewind(num_particles_file_id_cart)
             end if
             do i=1,BUFFER_SIZE
                write(num_particles_file_id_cart,*) num_particles_global_c(i)
             end do
             close(num_particles_file_id_cart)
          end if
       end if
       
       
       buffer_c(buffer_counter) = efield_energy_total_c
       buffer_q(buffer_counter) = efield_energy_total_q
       buffer_t(buffer_counter) = efield_energy_total_t
       ! This should be abstracted away...
       ! Each processor keeps a local buffer, when the buffer reaches a
       ! predetermined size, we reduce ther buffer with an addition on 
       ! process 0, who appends it to a file. Then reset the buffer.
       if( buffer_counter == BUFFER_SIZE ) then
          ! While using a sequential QNS solver, each processor contains the
          ! full information on the electric field energy, thus it is not 
          ! necessary to add the individual contributions.
          
          call sll_collective_reduce_real64( &
               sll_world_collective, &
               buffer_c, &
               BUFFER_SIZE, &
               MPI_SUM, &
               0, &
               buffer_result_c )
          
          call sll_collective_reduce_real64( &
               sll_world_collective, &
               buffer_q, &
               BUFFER_SIZE, &
               MPI_SUM, &
               0, &
               buffer_result_q )
          
          ! Use next line only if no communications are needed!
          !       buffer_result(:) = buffer(:)
          
          buffer_counter = 1
          if(sim%my_rank == 0) then
             open(efield_energy_file_id_qns,file="electric_field_energy_qns",&
                  position="append")
             if(itime == BUFFER_SIZE) then
                rewind(efield_energy_file_id_qns)
             end if
             buffer_result_q(:) = log(buffer_result_q(:))
             do i=1,BUFFER_SIZE
                write(efield_energy_file_id_qns,*) buffer_result_q(i)
             end do
             close(efield_energy_file_id_qns)
          end if
          
          if(sim%my_rank == 0) then
             open(efield_energy_file_id_cart,file="electric_field_energy_cart",&
                  position="append")
             if(itime == BUFFER_SIZE) then
                rewind(efield_energy_file_id_cart)
             end if
             buffer_result_c(:) = log(buffer_result_c(:))
             do i=1,BUFFER_SIZE
                write(efield_energy_file_id_cart,*) buffer_result_c(i)
             end do
             close(efield_energy_file_id_cart)
          end if
          
          if(sim%my_rank == 0) then
             open(efield_energy_file_id_teo,file="electric_field_energy_teo",&
                  position="append")
             if(itime == BUFFER_SIZE) then
                rewind(efield_energy_file_id_teo)
             end if
             buffer_t(:) = log(buffer_t(:))
             do i=1,BUFFER_SIZE
                write(efield_energy_file_id_teo,*) buffer_t(i)
             end do
             close(efield_energy_file_id_teo)
          end if

          
       else
          buffer_counter         = buffer_counter + 1
       end if
       efield_energy_total_c    = 0.0_f64
       efield_energy_total_q    = 0.0_f64
       efield_energy_total_t    = 0.0_f64
       ! Proceed to the advections in the spatial directions, 'x' and 'y'
       ! Reconfigure data. 
       
       ! what are the new local limits on x3 and x4? It is bothersome to have
       ! to make these calls...

       call advection_x1x2_c(sim,sim%dt)
       call advection_x1x2_q(sim,sim%dt)
       
       open(droite_test_pente,file="droite_test_pente",&
            position="append")
        write(droite_test_pente,*) -0.1533*(itime-1)*sim%dt + 2.8
        close(droite_test_pente)
    end do ! main loop
!#endif


 
#undef BUFFER_SIZE

  end subroutine run_4d_qns_mixed

  function efield_exact(t, k, eps, x, y) result(res)
    sll_real64, intent(in) :: t
    sll_real64, intent(in) :: k
    sll_real64, intent(in) :: eps
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64 :: res
    res =4.0_f64*eps*0.3677_f64*exp(-0.1533*t)*sin(k*x)*cos(-1.4156*t+0.536245)
  end function efield_exact

  subroutine advection_x1x2_c(sim,deltat)
    class(sll_simulation_4d_qns_mixed) :: sim
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
          call sim%interp_x1x2%compute_interpolants(sim%f_x1x2_c(:,:,k,l))
          do j=1,loc_sz_x2-1 ! last point excluded
             do i=1,loc_sz_x1-1 ! last point excluded
                global_indices = &
                     local_to_global_4D(sim%sequential_x1x2,(/i,j,k,l/))
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
                
                sim%f_x1x2_c(i,j,k,l) = &
                     sim%interp_x1x2%interpolate_value(eta1,eta2)

             end do
          end do
       end do
    end do

    do l=1,loc_sz_x4
       do k=1,loc_sz_x3
          do j=1,loc_sz_x2-1 ! last point excluded
             sim%f_x1x2_c(loc_sz_x1,j,k,l) = sim%f_x1x2_c(1,j,k,l)
          end do
       end do
    end do
 
   do l=1,loc_sz_x4
      do k=1,loc_sz_x3
         sim%f_x1x2_c(:,loc_sz_x2,k,l) = sim%f_x1x2_c(:,1,k,l)
      end do
   end do
 


  end subroutine advection_x1x2_c


  subroutine advection_x1x2_q(sim,deltat)
    class(sll_simulation_4d_qns_mixed) :: sim
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
          call sim%interp_x1x2%compute_interpolants(sim%f_x1x2_q(:,:,k,l))
          do j=1,loc_sz_x2-1 ! last point excluded
             do i=1,loc_sz_x1-1 ! last point excluded
                global_indices = &
                     local_to_global_4D(sim%sequential_x1x2,(/i,j,k,l/))
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
                
                sim%f_x1x2_q(i,j,k,l) = &
                     sim%interp_x1x2%interpolate_value(eta1,eta2)
                sim%f_x1x2_diff(i,j,k,l) = sim%f_x1x2_c(i,j,k,l) - &
                                           sim%f_x1x2_q(i,j,k,l)
             end do
          end do
       end do
    end do

    do l=1,loc_sz_x4
       do k=1,loc_sz_x3
          do j=1,loc_sz_x2-1 ! last point excluded
             sim%f_x1x2_q(loc_sz_x1,j,k,l) = sim%f_x1x2_q(1,j,k,l)
          end do
       end do
    end do
 
   do l=1,loc_sz_x4
      do k=1,loc_sz_x3
         sim%f_x1x2_q(:,loc_sz_x2,k,l) = sim%f_x1x2_q(:,1,k,l)
      end do
   end do



  end subroutine advection_x1x2_q

  subroutine advection_x3(sim,phi,deltat,efield_energy_total)
    class(sll_simulation_4d_qns_mixed) :: sim
    sll_real64, intent(in) :: deltat
    sll_int32 :: gi, gj, gk, gl
    sll_real64, dimension(1:2,1:2) :: inv_j
    sll_int32 :: loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 
    sll_int32, dimension(4) :: global_indices
    sll_real64 :: alpha3
    sll_real64 :: eta1, eta2, eta3, eta4
    sll_int32  :: i, j, k, l
    sll_real64, intent(out) :: efield_energy_total
    type(sll_scalar_field_2d_discrete_alt), pointer       :: phi
    sll_real64, dimension(1:2,1:2) :: jac_m
    sll_real64 :: ex
    sll_real64 :: ey
    sll_int32  :: nc_x3
    
    nc_x3 = sim%mesh2d_v%num_cells1
    call compute_local_sizes_4d( sim%sequential_x3x4, &
         loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 
    
    efield_energy_total = 0.0_f64
    
    ! Start with dt in vx...(x3)
    do l=1,loc_sz_x4 !sim%mesh2d_v%num_cells2+1
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             global_indices(1:2) = &
                  local_to_global_2D( sim%split_rho_layout, (/i,j/))
             eta1   =  sim%mesh2d_v%eta1_min + &
                  real(global_indices(1)-1,f64)*sim%mesh2d_v%delta_eta1
             eta2   =  sim%mesh2d_v%eta2_min + &
                  real(global_indices(2)-1,f64)*sim%mesh2d_v%delta_eta2
             
             inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
             jac_m  =  sim%transfx%jacobian_matrix(eta1,eta2)
             
             ex     =  - phi%first_deriv_eta1_value_at_point(eta1,eta2)
             ey     =  - phi%first_deriv_eta2_value_at_point(eta1,eta2)
             
             alpha3 = -sim%dt*(inv_j(1,1)*ex + inv_j(2,1)*ey)
             sim%f_x3x4_q(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                  nc_x3+1, &
                  sim%f_x3x4_q(i,j,:,l), &
                  alpha3 )
             ! Extra work to calculate the electric field energy. We should 
             ! consider placing this somewhere else, probably at greater
             ! expense.
             
             efield_energy_total = efield_energy_total + &
                  sim%mesh2d_v%delta_eta1*sim%mesh2d_v%delta_eta2*&
                  (((jac_m(2,2)*ex - jac_m(2,1)*ey)**2 + &
                  (jac_m(1,1)*ex - jac_m(1,2)*ey)**2))/&
                  (jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1))
             
          end do
       end do
    end do
    
    efield_energy_total = sqrt(efield_energy_total)
    
    
    
  end subroutine advection_x3
  
  
  subroutine advection_x4(sim,phi,deltat,efield_energy_total)
    class(sll_simulation_4d_qns_mixed) :: sim
    sll_real64, intent(in) :: deltat
    sll_int32 :: gi, gj, gk, gl
    sll_real64, dimension(1:2,1:2) :: inv_j
    sll_int32 :: loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 
    sll_int32, dimension(4) :: global_indices
    sll_real64 :: alpha4
    sll_real64 :: eta1, eta2, eta3, eta4
    sll_int32  :: i, j, k, l
    sll_real64, intent(out) :: efield_energy_total
    type(sll_scalar_field_2d_discrete_alt), pointer       :: phi
    sll_real64, dimension(1:2,1:2) :: jac_m
    sll_real64 :: ex
    sll_real64 :: ey
    sll_int32  :: nc_x4
    
    nc_x4 = sim%mesh2d_v%num_cells2
    call compute_local_sizes_4d( sim%sequential_x3x4, &
         loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 
    
    efield_energy_total = 0.0_f64
    ! dt in vy...(x4)
    do j=1,loc_sz_x2
       do i=1,loc_sz_x1
          do k=1,sim%mesh2d_v%num_cells1+1
             global_indices(1:2) = &
                  local_to_global_2D( sim%split_rho_layout, (/i,j/))
             eta1   =  sim%mesh2d_v%eta1_min + &
                  real(global_indices(1)-1,f64)*sim%mesh2d_v%delta_eta1
             eta2   =  sim%mesh2d_v%eta2_min + &
                  real(global_indices(2)-1,f64)*sim%mesh2d_v%delta_eta2
             inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
             
             ex     =  - phi%first_deriv_eta1_value_at_point(eta1,eta2)
             ey     =  - phi%first_deriv_eta2_value_at_point(eta1,eta2)
             alpha4 = -sim%dt*(inv_j(1,2)*ex + inv_j(2,2)*ey)
             sim%f_x3x4_q(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
                  nc_x4+1, &
                  sim%f_x3x4_q(i,j,k,:), &
                  alpha4 )
             efield_energy_total = efield_energy_total + &
                  sim%mesh2d_v%delta_eta1*sim%mesh2d_v%delta_eta2*&
                  (((jac_m(2,2)*ex - jac_m(2,1)*ey)**2 + &
                  (jac_m(1,1)*ex - jac_m(1,2)*ey)**2))/&
                  (jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1))
          end do
       end do
    end do
    
    efield_energy_total = sqrt(efield_energy_total)
    
    
  end subroutine advection_x4

  subroutine delete_4d_qns_mixed( sim )
    type(sll_simulation_4d_qns_mixed) :: sim
    sll_int32 :: ierr
    SLL_DEALLOCATE( sim%f_x1x2_q, ierr )
    SLL_DEALLOCATE( sim%f_x3x4_q, ierr )
    SLL_DEALLOCATE_ARRAY( sim%partial_reduction_q, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_full, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_split_q, ierr )
   ! SLL_DEALLOCATE_ARRAY( sim%phi_x1, ierr )
   ! SLL_DEALLOCATE_ARRAY( sim%phi_x2, ierr )
   ! SLL_DEALLOCATE_ARRAY( sim%phi_split, ierr )
    call delete( sim%sequential_x1x2 )
    call delete( sim%sequential_x3x4 )
    call delete( sim%rho_full_layout )
    call delete( sim%rho_seq_x2 )
    call delete( sim%split_rho_layout )
    call delete( sim%split_to_full )
    call delete( sim%efld_seqx1_to_seqx2 )
    call delete( sim%efld_seqx2_to_split )
    call delete( sim%seqx1x2_to_seqx3x4 )
    call delete( sim%seqx3x4_to_seqx1x2 )
    call delete( sim%interp_x1x2 )
    call delete( sim%interp_x3 )
    call delete( sim%interp_x4 )
    !SLL_DEALLOCATE_ARRAY( sim%efield_x1, ierr )
    !SLL_DEALLOCATE_ARRAY( sim%efield_x2, ierr )
    !SLL_DEALLOCATE_ARRAY( sim%efield_split, ierr )
  end subroutine delete_4d_qns_mixed

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
    sll_int32                       :: i, j, k, l
    
    numpts3 = mv%num_cells1+1
    numpts4 = mv%num_cells2+1
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
    sll_comp64, dimension(:,:), intent(inout) :: efield_x2
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
       ex = real(efield_x2(i,num_pts_x2),f64)
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
  



  subroutine compute_average_f( &
       sim,&
       mx,&
       rho, &
       density_tot )
    class(sll_simulation_4d_qns_mixed)     :: sim
    type(sll_logical_mesh_2d), pointer     :: mx
    sll_real64, intent(inout), dimension(:,:)     :: rho     ! local rho
    sll_real64, intent(out)              :: density_tot
    
    ! local sizes in the split directions have to be given by caller.
    sll_int32   :: numpts1
    sll_int32   :: numpts2
    sll_real64  :: delta1
    sll_real64  :: delta2
    sll_real64  :: lenght1
    sll_real64  :: lenght2
    sll_real64  :: length_total
    sll_int32 :: i, j
    sll_real64, dimension(1:2,1:2) :: jac_m
    sll_real64   :: eta1, eta2
    numpts1 = mx%num_cells1+1
    numpts2 = mx%num_cells2+1
    delta1  = mx%delta_eta1
    delta2  = mx%delta_eta2
    lenght1 = mx%eta1_max- mx%eta1_min
    lenght2 = mx%eta2_max- mx%eta2_min
    
    density_tot  = 0.0_f64
    length_total = 0.0_f64
    do j=1,numpts2-1
       do i=1,numpts1-1
          eta1 = sim%mesh2d_x%eta1_min + (i-1)*sim%mesh2d_x%delta_eta1
          eta2 = sim%mesh2d_x%eta2_min + (j-1)*sim%mesh2d_x%delta_eta2
          jac_m  = sim%transfx%jacobian_matrix(eta1,eta2)
          density_tot = density_tot + rho(i,j)*delta1*delta2*&
               (jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1))
          !print*, jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1)
          length_total = length_total + &
               delta1*delta2*(jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1))
       end do
    end do
    print*, length_total,(4*sll_pi)**2
    density_tot = density_tot/ (length_total)
  end subroutine compute_average_f

  ! Some ad-hoc functions to prepare data for allgather operations. Should
  ! consider putting this elsewhere.

  subroutine compute_displacements_array( layout, collective_size, disps )
    type(layout_2D), pointer                           :: layout
    sll_int32, intent(in)                              :: collective_size
    sll_int32, dimension(collective_size), intent(out) :: disps
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    sll_int32 :: i,j
    sll_int32 :: counter
    sll_int32 :: rank

    counter  = 0
    disps(1) = counter
    do rank=1,collective_size-1
       imin = get_layout_i_min( layout, rank-1 )
       imax = get_layout_i_max( layout, rank-1 )
       jmin = get_layout_j_min( layout, rank-1 )
       jmax = get_layout_j_max( layout, rank-1 )
       size_i      = imax - imin + 1
       size_j      = jmax - jmin + 1
       counter     = counter + size_i*size_j
       disps(rank+1) = counter
    end do
  end subroutine compute_displacements_array

  subroutine load_buffer( layout, data, buffer )
    type(layout_2D), pointer   :: layout
    sll_real64, dimension(:,:), intent(in) :: data
    sll_real64, dimension(:),  intent(out) :: buffer
    sll_int32 :: myrank
    sll_int32 :: data_size
    sll_int32 :: send_size
    type(sll_collective_t), pointer :: col
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    sll_int32 :: i,j
    sll_int32 :: counter

    col => get_layout_collective( layout )
    myrank = sll_get_collective_rank( col )
    data_size = size(data,1)*size(data,2)
!print *, 'data1: ', size(data,1), 'data2:', size(data,2)
    imin = get_layout_i_min( layout, myrank )
    imax = get_layout_i_max( layout, myrank )
    jmin = get_layout_j_min( layout, myrank )
    jmax = get_layout_j_max( layout, myrank )
    size_i = imax - imin + 1
    size_j = jmax - jmin + 1

!!$    if( data_size .ne. size_i*size_j ) then
!!$       print *, 'function load_buffer():'
!!$       print *, 'size(data) = ', size(data,1), size(data,2)
!!$       print *, 'warning from rank ', myrank
!!$       print *, 'there seems to be a discrepancy between the data size ', &
!!$            'passed and the size declared in the layout.'
!!$       print *, 'data size = ', data_size, 'size from layout = ', size_i*size_j
!!$    end if
!!$print *, 'size_j', size_j, 'size_i', size_i
    counter=0
    do j=1,size_j
       do i=1,size_i
          counter = counter + 1
          buffer(counter) = data(i,j)
       end do
    end do

  end subroutine load_buffer

  function receive_counts_array( layout, n ) result(rc)
    type(layout_2D), pointer :: layout
    sll_int32, intent(in)    :: n
    sll_int32, dimension(n) :: rc
    sll_int32 :: i
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j

    do i=0,n-1
       imin = get_layout_i_min( layout, i )
       imax = get_layout_i_max( layout, i )
       jmin = get_layout_j_min( layout, i )
       jmax = get_layout_j_max( layout, i )
       size_i = imax - imin + 1
       size_j = jmax - jmin + 1
       rc(i+1)  = size_i*size_j
    end do
  end function receive_counts_array

  subroutine unload_buffer( layout, buffer, data )
    type(layout_2D), pointer                :: layout
    sll_real64, dimension(:,:), intent(out) :: data
    sll_real64, dimension(:), intent(in)    :: buffer
    sll_int32 :: col_sz
    type(sll_collective_t), pointer :: col
    sll_int32 :: i, j
    sll_int32 :: box
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    sll_int32 :: pos            ! position in buffer
    col => get_layout_collective( layout )
    col_sz = sll_get_collective_size( col )

    ! loop over all the boxes in the layout and fill the data array by chunks.
    pos = 1
    do box=0,col_sz-1
       imin = get_layout_i_min( layout, box )
       imax = get_layout_i_max( layout, box )
       jmin = get_layout_j_min( layout, box )
       jmax = get_layout_j_max( layout, box )
       ! this will fill the data array in whatever order that the boxes
       ! are ordered.
       do j=jmin,jmax
          do i=imin,imax
             data(i,j) = buffer(pos)
             pos       = pos + 1
          end do
       end do
    end do
  end subroutine unload_buffer

  
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



end module sll_simulation_4d_qns_mixed_module
