module sll_simulation_4d_qns_general_multipatch_module

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"

  use sll_module_scalar_field_2d_multipatch
  use sll_general_coordinate_elliptic_solver_multipatch_module
  use sll_coordinate_transformation_multipatch_module
  use sll_distribution_function_4d_multipatch_module
!  use sll_collective
!  use sll_remapper
!  use sll_constants
  use sll_cubic_spline_interpolator_1d
!  use sll_cubic_spline_interpolator_2d
  use sll_simulation_base
!  use sll_logical_meshes
  use sll_parallel_array_initializer_module
!  use sll_coordinate_transformation_2d_base_module
  use sll_gnuplot_parallel
  use sll_general_coordinate_elliptic_solver_module
  use sll_common_array_initializers_module
!  use sll_module_scalar_field_2d_base
!  use sll_module_scalar_field_2d_alternative
!  use sll_arbitrary_degree_spline_interpolator_1d_module
!  use sll_timer
  implicit none

  type, extends(sll_simulation_base_class) :: sll_simulation_4d_qns_general_multipatch
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
     ! the degree of interpolation in vx and vy
     sll_int32  :: spline_degree_vx
     sll_int32  :: spline_degree_vy
     ! for QNS boundary conditions
     sll_int32  :: bc_eta1_0
     sll_int32  :: bc_eta1_1
     sll_int32  :: bc_eta2_0
     sll_int32  :: bc_eta2_1
     sll_int32  :: bc_vx_0
     sll_int32  :: bc_vx_1
     sll_int32  :: bc_vy_0
     sll_int32  :: bc_vy_1
     ! the logical meshes are split in two one for space, one for velocity
     type(sll_logical_mesh_2d), pointer    :: mesh2d_v
     ! This simulation only applies a coordinate transformation to the spatial
     ! coordinates.
     type(sll_coordinate_transformation_multipatch_2d), pointer :: transfx
     type(general_coordinate_elliptic_solver_mp), pointer      :: qns


     !  number dignostics for the simulation
      sll_int32  ::number_diags

      !  type quadrature for solver in direction eta1 and eta2
      sll_int32 :: quadrature_type1,quadrature_type2 

     ! distribution functions. There are several because each array represents
     ! a differently shaped chunk of memory. In this example, each chunk 
     ! allows sequential operations in one given direction. f_x1x2 should 
     ! permit to carry out sequential operations in x1 and x2 for ex.
     sll_real64, dimension(:,:,:,:), pointer     :: f_x1x2 
     sll_real64, dimension(:,:,:,:), pointer     :: f_x3x4
     sll_real64, dimension(:,:,:), allocatable   :: partial_reduction
     sll_real64, dimension(:,:), pointer         :: rho_full 
     sll_real64, dimension(:,:), allocatable     :: rho_x2 
     sll_real64, dimension(:,:), allocatable     :: rho_split


     !--> diagnostics norm
     sll_real64, dimension(:),pointer :: diag_masse
     sll_real64, dimension(:),pointer :: diag_norm_L1
     sll_real64, dimension(:),pointer :: diag_norm_L2
     sll_real64, dimension(:),pointer :: diag_norm_Linf
     sll_real64, dimension(:),pointer :: diag_entropy_kin

     !--> diagnostics energydiag_nrj_kin
     sll_real64, dimension(:),pointer :: diag_nrj_kin
     sll_real64, dimension(:),pointer :: diag_nrj_pot
     sll_real64, dimension(:),pointer :: diag_nrj_tot

     !---> For diagnostic saving
     sll_int32 :: count_save_diag

     !---> To stocke the values of Jacobian in mesh points

     ! sll_real64,dimension(:,:,:,:),pointer :: values_jacobian_mat
     ! sll_real64,dimension(:,:),    pointer :: values_jacobian
     ! sll_real64,dimension(:,:,:,:),pointer :: values_jacobian_matinv
     ! sll_real64,dimension(:,:),pointer     :: point_x,point_y
     ! sll_real64,dimension(:,:), pointer    :: values_ex,values_ey

     ! ---> point mesh logical
     sll_real64,dimension(:),pointer :: pt_eta1
     sll_real64,dimension(:),pointer :: pt_eta2
     ! for remap
     type(layout_4D), pointer :: sequential_x1x2
     type(layout_4D), pointer :: sequential_x3x4
     type(layout_2D), pointer :: rho_full_layout
     type(layout_2D), pointer :: rho_seq_x2
     type(layout_2D), pointer :: split_rho_layout ! not sequential in any dir.
     type(remap_plan_2D_real64), pointer :: split_to_full
     type(remap_plan_2D_real64), pointer :: seqx1_to_seqx2
     ! remaps for the electric field data
     !     type(remap_plan_2D), pointer :: efld_split_to_seqx1
     type(remap_plan_2D_comp64), pointer :: efld_seqx1_to_seqx2
     type(remap_plan_2D_comp64), pointer :: efld_seqx2_to_split
     type(remap_plan_4D_real64), pointer :: seqx1x2_to_seqx3x4
     type(remap_plan_4D_real64), pointer :: seqx3x4_to_seqx1x2
     ! interpolators
     type(cubic_spline_1d_interpolator) :: interp_x3
     type(cubic_spline_1d_interpolator) :: interp_x4
     ! for distribution function initializer:
     procedure(sll_scalar_initializer_4d), nopass, pointer :: init_func
     sll_real64, dimension(:), pointer :: params
     ! for general coordinate QNS, analytical fields
     procedure(two_var_parametrizable_function),nopass,pointer :: a11_f
     procedure(two_var_parametrizable_function),nopass,pointer :: a12_f
     procedure(two_var_parametrizable_function),nopass,pointer :: a21_f
     procedure(two_var_parametrizable_function),nopass,pointer :: a22_f
     procedure(two_var_parametrizable_function),nopass,pointer :: b1_f
     procedure(two_var_parametrizable_function),nopass,pointer :: b2_f
     procedure(two_var_parametrizable_function),nopass,pointer :: der1_b1_f
     procedure(two_var_parametrizable_function),nopass,pointer :: der1_b2_f
     procedure(two_var_parametrizable_function),nopass,pointer :: der2_b1_f
     procedure(two_var_parametrizable_function),nopass,pointer :: der2_b2_f
     procedure(two_var_parametrizable_function),nopass,pointer :: c_f
     procedure(two_var_parametrizable_function),nopass,pointer :: elec_field_ext_1
     procedure(two_var_parametrizable_function),nopass,pointer :: elec_field_ext_2
     sll_real64, dimension(:), pointer :: a11_f_params
     sll_real64, dimension(:), pointer :: a12_f_params
     sll_real64, dimension(:), pointer :: a21_f_params
     sll_real64, dimension(:), pointer :: a22_f_params
     sll_real64, dimension(:), pointer :: b1_f_params
     sll_real64, dimension(:), pointer :: b2_f_params
     sll_real64, dimension(:), pointer :: c_f_params
     sll_real64, dimension(:), pointer :: elec_field_ext_f_params
   contains
     procedure, pass(sim) :: run => run_4d_qns_general_mp
     procedure, pass(sim) :: init_from_file => init_4d_qns_gen_mp
     procedure, pass(sim) :: initialize => initialize_4d_qns_gen_mp
  end type sll_simulation_4d_qns_general_multipatch

  interface sll_delete
     module procedure delete_4d_qns_gen_mp
  end interface sll_delete

  ! interface initialize
  !    module procedure initialize_4d_qns_general_multipatch
  ! end interface initialize

contains

  ! Tentative function to initialize the simulation object 'manually'.
  subroutine initialize_4d_qns_gen_mp(&
       sim, &
       mesh2d_v, &
       transformation_x, &
       init_func, &
       params,&
       a11_f,&
       a11_f_params, &
       a12_f,&
       a12_f_params, &
       a21_f,&
       a21_f_params, &
       a22_f,&
       a22_f_params, &
       b1_f, &
       b1_f_params, &
       der1_b1_f,&
       der2_b1_f,&
       b2_f, &
       b2_f_params, &
       der1_b2_f,&
       der2_b2_f,&
       c_f,&
       c_f_params, &
       spline_degre_vx,&
       spline_degre_vy,&
       bc_eta1_0,&
       bc_eta1_1,&
       bc_eta2_0,&
       bc_eta2_1,&
       bc_vx_0,&
       bc_vx_1,&
       bc_vy_0,&
       bc_vy_1,&
       quadrature_type1,&
       quadrature_type2,&
       electric_field_ext_1,&
       electric_field_ext_2,&
       elec_field_ext_f_params,&
       number_diags)
    
    class(sll_simulation_4d_qns_general_multipatch),   intent(inout)      :: sim
    type(sll_coordinate_transformation_multipatch_2d), intent(in), target :: transformation_x
    type(sll_logical_mesh_2d), pointer         :: mesh2d_v
    procedure(sll_scalar_initializer_4d)       :: init_func
    sll_real64, dimension(:), target           :: params
    procedure(two_var_parametrizable_function) :: a11_f
    procedure(two_var_parametrizable_function) :: a12_f
    procedure(two_var_parametrizable_function) :: a21_f
    procedure(two_var_parametrizable_function) :: a22_f
    procedure(two_var_parametrizable_function) :: b1_f
    procedure(two_var_parametrizable_function) :: b2_f
    procedure(two_var_parametrizable_function) :: der1_b1_f
    procedure(two_var_parametrizable_function) :: der1_b2_f
    procedure(two_var_parametrizable_function) :: der2_b1_f
    procedure(two_var_parametrizable_function) :: der2_b2_f
    procedure(two_var_parametrizable_function) :: c_f
    procedure(two_var_parametrizable_function) :: electric_field_ext_1
    procedure(two_var_parametrizable_function) :: electric_field_ext_2
    sll_real64, dimension(:), intent(in) :: a11_f_params
    sll_real64, dimension(:), intent(in) :: a12_f_params
    sll_real64, dimension(:), intent(in) :: a21_f_params
    sll_real64, dimension(:), intent(in) :: a22_f_params
    sll_real64, dimension(:), intent(in) :: b1_f_params
    sll_real64, dimension(:), intent(in) :: b2_f_params
    sll_real64, dimension(:), intent(in) :: c_f_params
    sll_real64, dimension(:), intent(in) :: elec_field_ext_f_params
    sll_int32  :: spline_degre_vx
    sll_int32  :: spline_degre_vy
    sll_int32  :: number_diags
    sll_int32  :: bc_eta1_0
    sll_int32  :: bc_eta1_1
    sll_int32  :: bc_eta2_0
    sll_int32  :: bc_eta2_1
    sll_int32  :: bc_vx_0
    sll_int32  :: bc_vx_1
    sll_int32  :: bc_vy_0
    sll_int32  :: bc_vy_1
    sll_int32  :: quadrature_type1,quadrature_type2
    sll_int32 :: ierr
    
    ! sim%mesh2d_x  => transformation_x%mesh
    sim%mesh2d_v  => mesh2d_v
    sim%transfx   => transformation_x
    sim%init_func => init_func
    sim%params    => params
    sim%a11_f     => a11_f
    sim%a12_f     => a12_f
    sim%a21_f     => a21_f
    sim%a22_f     => a22_f
    sim%b1_f      => b1_f
    sim%b2_f      => b2_f
    sim%der1_b1_f  => der1_b1_f
    sim%der1_b2_f  => der1_b2_f
    sim%der2_b1_f  => der2_b1_f
    sim%der2_b2_f  => der2_b2_f
    sim%c_f       => c_f
    sim%spline_degree_vx = spline_degre_vx
    sim%spline_degree_vy = spline_degre_vy
    sim%number_diags     = number_diags
    sim%quadrature_type1 = quadrature_type1
    sim%quadrature_type2 = quadrature_type2
    sim%count_save_diag  = 0
    sim%bc_eta1_0   = bc_eta1_0
    sim%bc_eta1_1   = bc_eta1_1
    sim%bc_eta2_0   = bc_eta2_0
    sim%bc_eta2_1   = bc_eta2_1

    sim%bc_vx_0     = bc_vx_0
    sim%bc_vx_1     = bc_vx_1
    sim%bc_vy_0     = bc_vy_0
    sim%bc_vy_1     = bc_vy_1

    sim%elec_field_ext_1 => electric_field_ext_1
    sim%elec_field_ext_2 => electric_field_ext_2
    
    SLL_ALLOCATE(sim%a11_f_params(size(a11_f_params)),ierr)
    SLL_ALLOCATE(sim%a12_f_params(size(a12_f_params)),ierr)
    SLL_ALLOCATE(sim%a21_f_params(size(a21_f_params)),ierr)
    SLL_ALLOCATE(sim%a22_f_params(size(a22_f_params)),ierr)
    SLL_ALLOCATE(sim%b1_f_params(size(b1_f_params)),ierr)
    SLL_ALLOCATE(sim%b2_f_params(size(b2_f_params)),ierr)
    SLL_ALLOCATE(sim%c_f_params(size(c_f_params)),ierr)
    SLL_ALLOCATE(sim%elec_field_ext_f_params(size(elec_field_ext_f_params)),ierr)

    sim%a11_f_params(:) = a11_f_params
    sim%a12_f_params(:) = a12_f_params
    sim%a21_f_params(:) = a21_f_params
    sim%a22_f_params(:) = a22_f_params
    sim%b1_f_params(:)  = b1_f_params
    sim%b2_f_params(:)  = b2_f_params
    sim%c_f_params(:)   = c_f_params
    sim%elec_field_ext_f_params(:) = elec_field_ext_f_params
    
    ! SLL_ALLOCATE(sim%values_jacobian_mat(sim%mesh2d_x%num_cells1 +1,sim%mesh2d_x%num_cells2 +1,2,2),ierr)
    ! sim%values_jacobian_mat(:,:,:,:) = 0.0_f64
    ! SLL_ALLOCATE(sim%values_jacobian_matinv(sim%mesh2d_x%num_cells1 +1,sim%mesh2d_x%num_cells2 +1,2,2),ierr)
    ! sim%values_jacobian_matinv(:,:,:,:) = 0.0_f64
    ! SLL_ALLOCATE(sim%values_jacobian(sim%mesh2d_x%num_cells1 +1,sim%mesh2d_x%num_cells2 +1),ierr)
    ! sim%values_jacobian(:,:) = 0.0_f64
    ! SLL_ALLOCATE(sim%pt_eta1(sim%mesh2d_x%num_cells1 +1),ierr)
    ! SLL_ALLOCATE(sim%pt_eta2(sim%mesh2d_x%num_cells2 +1),ierr)
    ! SLL_ALLOCATE(sim%point_x(sim%mesh2d_x%num_cells1 +1,sim%mesh2d_x%num_cells2 +1),ierr)
    ! SLL_ALLOCATE(sim%point_y(sim%mesh2d_x%num_cells1 +1,sim%mesh2d_x%num_cells2 +1),ierr)
    ! SLL_ALLOCATE(sim%values_ex(sim%mesh2d_x%num_cells1 +1,sim%mesh2d_x%num_cells2 +1 ),ierr)
    ! sim%values_ex = 0.0_f64
    ! SLL_ALLOCATE(sim%values_ey(sim%mesh2d_x%num_cells1 +1,sim%mesh2d_x%num_cells2 +1),ierr)
    ! sim%values_ey = 0.0_f64
  end subroutine initialize_4d_qns_gen_mp


  subroutine init_4d_qns_gen_mp( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_qns_general_multipatch), intent(inout) :: sim
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
    print*, 'number iterations', number_iterations
    ! In this particular simulation, since the system is periodic, the number
    ! of points is the same as the number of cells in all directions.
    sim%nc_x1 = num_cells_x1
    sim%nc_x2 = num_cells_x2
    sim%nc_x3 = num_cells_x3
    sim%nc_x4 = num_cells_x4
  end subroutine init_4d_qns_gen_mp

  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_4d_qns_general_mp(sim)
    class(sll_simulation_4d_qns_general_multipatch), intent(inout) :: sim
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
    sll_real64 :: x
    sll_real64 :: y
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
    sll_int32,  dimension(1:2)     :: gi     ! for storing global indices
    sll_int32,  dimension(1:4)     :: gi4d   ! for storing global indices
    sll_real64 :: efield_energy_total
    ! The following could probably be abstracted for convenience
#define BUFFER_SIZE sim%number_diags
    sll_real64, dimension(BUFFER_SIZE) :: buffer_energy
    sll_real64, dimension(BUFFER_SIZE) :: buffer_energy_result
    sll_real64, dimension(BUFFER_SIZE) :: num_particles_local
    sll_real64, dimension(BUFFER_SIZE) :: num_particles_global
    sll_real64, dimension(1) :: tmp1
    sll_real64 :: tmp,numpart
    sll_int32  :: buffer_counter
    sll_int32  :: efield_energy_file_id
    sll_int32  :: num_particles_file_id
    sll_int32  :: global_indices(4)
    sll_int32  :: iplot
    character(len=4) :: cplot
    class(sll_scalar_field_multipatch_2d), pointer      :: a11_field_mat
    class(sll_scalar_field_multipatch_2d), pointer      :: a21_field_mat
    class(sll_scalar_field_multipatch_2d), pointer      :: a12_field_mat
    class(sll_scalar_field_multipatch_2d), pointer      :: a22_field_mat
    class(sll_scalar_field_multipatch_2d), pointer      :: b1_field_vect
    class(sll_scalar_field_multipatch_2d), pointer      :: b2_field_vect
    class(sll_scalar_field_multipatch_2d), pointer      :: c_field
    type(sll_scalar_field_multipatch_2d), pointer      :: elec_field_ext_1
    type(sll_scalar_field_multipatch_2d), pointer      :: elec_field_ext_2
    type(sll_scalar_field_multipatch_2d), pointer      :: rho
    type(sll_scalar_field_multipatch_2d), pointer      :: phi
    type(sll_scalar_field_multipatch_2d), pointer      :: layer_x1x2
    type(sll_logical_mesh_2d), pointer                         :: logical_m
    class(sll_coordinate_transformation_2d_nurbs), pointer     :: transf
    type(sll_distribution_function_4d_multipatch), pointer     :: f_mp

    sll_real64, dimension(:), allocatable :: send_buf
    sll_real64, dimension(:), allocatable :: recv_buf
    sll_int32,  dimension(:), allocatable :: recv_sz
    sll_int32,  dimension(:), allocatable :: disps ! for allgatherv operation
    sll_real64 :: val_a11, val_a12, val_a21, val_a22
    sll_real64 :: val_b1,  val_b2,  val_c
    sll_real64 :: val_rho, val_phi, val_phi_exacte
    sll_real64 :: density_tot
    sll_real64 :: x1
    sll_real64 :: x2
    sll_int32  :: send_size   ! for allgatherv operation
    sll_int32 :: droite_test_pente
    sll_int32 :: num_patches
    sll_int32 :: ipatch
    sll_int32 :: num_pts1
    sll_int32 :: num_pts2


    ! only for debugging...
!!$    sll_real64, dimension(:,:), allocatable :: ex_field
!!$    sll_real64, dimension(:,:), allocatable :: ey_field
    ! time variables
    type(sll_time_mark)  :: t0 
    type(sll_time_mark)  :: t1
    sll_real64 :: time 
    sll_int32 :: size_diag    
    
    ! compute Jacobian and logical mesh points
    call compute_values_jacobian_and_mesh_points(sim)

    ! Start with the fields
    a11_field_mat => new_scalar_field_multipatch_2d("a11", sim%transfx)
    a12_field_mat => new_scalar_field_multipatch_2d("a12", sim%transfx)
    a21_field_mat => new_scalar_field_multipatch_2d("a21", sim%transfx)
    a22_field_mat => new_scalar_field_multipatch_2d("a22", sim%transfx)
    b1_field_vect => new_scalar_field_multipatch_2d("b1", sim%transfx)
    b2_field_vect => new_scalar_field_multipatch_2d("b2", sim%transfx)
    c_field       => new_scalar_field_multipatch_2d("c", sim%transfx)
    layer_x1x2    => new_scalar_field_multipatch_2d("layer_x1x2", sim%transfx)
    phi => new_scalar_field_multipatch_2d("potential_field_phi", sim%transfx)
    rho => new_scalar_field_multipatch_2d("rho_field_check", sim%transfx)


    ! elec_field_ext_1 => new_scalar_field_multipatch_2d("E1_ext", sim%transfx)    
    ! elec_field_ext_2 => new_scalar_field_multipatch_2d("E2_ext", sim%transfx)    

    call a11_field_mat%allocate_memory()
    call a12_field_mat%allocate_memory()
    call a21_field_mat%allocate_memory()
    call a22_field_mat%allocate_memory()
    call b1_field_vect%allocate_memory()
    call b2_field_vect%allocate_memory()
    call c_field%allocate_memory()
    call phi%allocate_memory()

    ! call elec_field_ext_1_field_mat%allocate_memory()
    ! call elec_field_ext_2_field_mat%allocate_memory()
    
    


    num_patches = sim%transfx%get_number_patches()

    do ipatch= 0,num_patches-1
       ! Please get rid of these 'fixes' whenever it is decided that gfortran 4.6
       ! is no longer supported by Selalib.
       !     m        => sim%transfx%get_logical_mesh(ipatch)
       ! logical_m => sim%transfx%transfs(ipatch+1)%t%mesh
       !     transf   => sim%transfx%get_transformation(ipatch)
       transf => sim%transfx%transfs(ipatch+1)%t

       num_pts1 = sim%transfx%get_num_cells_eta1(ipatch) + 1! logical_m%num_cells1+1
       num_pts2 = sim%transfx%get_num_cells_eta2(ipatch) + 1!logical_m%num_cells2+1
       delta1   = sim%transfx%get_delta_eta1(ipatch)!logical_m%delta_eta1
       delta2   = sim%transfx%get_delta_eta2(ipatch)!logical_m%delta_eta2
       eta1_min = sim%transfx%get_eta1_min(ipatch)!logical_m%eta1_min
       eta2_min = sim%transfx%get_eta2_min(ipatch)!logical_m%eta2_min

       print *, "num_pts = ", num_pts1
       print *, "num_pts2 = ", num_pts2
       print *, "num_pts = ", delta1
       print *, "num_pts = ", delta2
       print *, "num_pts = ", eta1_min
       print *, "num_pts = ", eta2_min

       do j=1,num_pts1
          eta2 = eta2_min + (j-1)*delta2
          do i=1,num_pts2
             ! here it is assumed that the eta_min's are = 0. This is supposed
             ! to be the case for NURBS transformations.
             eta1 = eta1_min + (i-1)*delta1
             x1 = transf%x1(eta1,eta2)
             x2 = transf%x2(eta1,eta2)
             val_a11 = sim%a11_f(x1, x2, sim%a11_f_params)
             val_a12 = sim%a12_f(x1, x2, sim%a12_f_params)
             val_a21 = sim%a21_f(x1, x2, sim%a21_f_params)
             val_a22 = sim%a22_f(x1, x2, sim%a22_f_params)
             val_b1  = sim%b1_f(x1, x2, sim%b1_f_params)
             val_b2  = sim%b2_f(x1, x2, sim%b2_f_params)
             val_c   = sim%c_f(x1, x2, sim%c_f_params)
             val_phi = 0.0_f64
             call a11_field_mat%set_value_at_indices ( i, j, ipatch, val_a11 ) 
             call a12_field_mat%set_value_at_indices ( i, j, ipatch, val_a12 ) 
             call a21_field_mat%set_value_at_indices ( i, j, ipatch, val_a21 ) 
             call a22_field_mat%set_value_at_indices ( i, j, ipatch, val_a22 ) 
             call b1_field_vect%set_value_at_indices ( i, j, ipatch, val_b1 ) 
             call b2_field_vect%set_value_at_indices ( i, j, ipatch, val_b2 ) 
             call c_field%set_value_at_indices  ( i, j, ipatch, val_c ) 
             call phi%set_value_at_indices( i, j, ipatch, val_phi)
          end do
       end do
    end do

    call a11_field_mat%update_interpolation_coefficients()
    call a12_field_mat%update_interpolation_coefficients()
    call a21_field_mat%update_interpolation_coefficients()
    call a22_field_mat%update_interpolation_coefficients()
    call b1_field_vect%update_interpolation_coefficients()
    call b2_field_vect%update_interpolation_coefficients()
    call c_field%update_interpolation_coefficients()
    call phi%update_interpolation_coefficients()

    buffer_counter = 1
    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)
    
    ! Consider deleting the following lines :
    SLL_ALLOCATE(recv_sz(sim%world_size),ierr)
    SLL_ALLOCATE(disps(sim%world_size),ierr)
    
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

    !--> Initialization diagnostics for the norm
    size_diag  =  int(sim%num_iterations/BUFFER_SIZE) + 1
    SLL_ALLOCATE(sim%diag_masse(size_diag),ierr)
    SLL_ALLOCATE(sim%diag_norm_L1(size_diag),ierr)
    SLL_ALLOCATE(sim%diag_norm_L2(size_diag),ierr)
    SLL_ALLOCATE(sim%diag_norm_Linf(size_diag),ierr)
    SLL_ALLOCATE(sim%diag_entropy_kin(size_diag),ierr)
    
    !--> Initialization diagnostics for the energy
    SLL_ALLOCATE(sim%diag_nrj_kin(size_diag),ierr)
    SLL_ALLOCATE(sim%diag_nrj_pot(size_diag),ierr)
    SLL_ALLOCATE(sim%diag_nrj_tot(size_diag),ierr)
    
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

    f_mp => sll_new_distribution_function_4d_multipatch( sll_world_collective, &
         sim%transfx, sim%mesh2d_v, sim%nproc_x1, sim%nproc_x2 )

    call f_mp%initialize( sim%init_func, sim%params )
 

    print *, 'sequential_x3x4 mode...'

    itemp = sim%nproc_x3
    sim%nproc_x3 = sim%nproc_x1
    sim%nproc_x1 = itemp
    ! switch x2 and x4
    itemp = sim%nproc_x4
    sim%nproc_x4 = sim%nproc_x2 
    sim%nproc_x2 = itemp
    
    print *, 'sequential x1x2 mode...'
  

    call rho%allocate_memory()
    call compute_charge_density_multipatch(f_mp, rho)
    call rho%update_interpolation_coefficients( )
  
    if(sim%my_rank == 0) then
       call rho%write_to_file(0)
    end if
          
    ! call f_mp%set_to_sequential_x1x2()
    ! call f_mp%set_to_sequential_x3x4()

    print*, 'initialize qns'
    ! Initialize the poisson plan before going into the main loop.
    sim%qns => new_general_elliptic_solver_mp( &
         sim%quadrature_type1,& !ES_GAUSS_LEGENDRE, &  ! put in arguments
         sim%quadrature_type2,& !ES_GAUSS_LEGENDRE, &  ! put in arguments
         sim%transfx)

    print*, 'factorise matrice qns'
    
    call factorize_mat_es_mp(&
         sim%qns, & 
         a11_field_mat, &
         a12_field_mat, &
         a21_field_mat, &
         a22_field_mat, &
         b1_field_vect, &
         b2_field_vect, &
         c_field)

    print*, '--- end factorise matrice qns'

    call solve_general_coordinates_elliptic_eq_mp(&
       sim%qns,&
       rho,&
       phi)

    print*, '--- end solve'

    if(sim%my_rank == 0) then
       call phi%write_to_file(0)
    end if

#if 0

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
         sim%mesh2d_x%num_cells1 +1, &
         sim%mesh2d_x%num_cells2 +1, &
         sim%mesh2d_x%eta1_min, &
         sim%mesh2d_x%eta1_max, &
         sim%mesh2d_x%eta2_min, &
         sim%mesh2d_x%eta2_max, &
         sim%bc_eta1_0, &
         sim%bc_eta1_1, &
         sim%bc_eta2_0, &
         sim%bc_eta2_1)!, &
         ! sim%spline_degree_eta1, &
         ! sim%spline_degree_eta2)

    call advection_x1x2(sim,0.5*sim%dt)

   
    ! other test cases use periodic bc's here...        
!!$    call sim%interp_x3%initialize( &
!!$         nc_x3+1, &
!!$         vmin3, &
!!$         vmax3, &
!!$         sim%bc_vx_0,&
!!$         sim%bc_vx_1,&
!!$         sim%spline_degree_vx)

    call sim%interp_x3%initialize( &
         nc_x3+1, &
         vmin3, &
         vmax3, &
         SLL_HERMITE)!sim%bc_vx_0)
    
    call sim%interp_x4%initialize( &
         nc_x4+1, &
         vmin4, &
         vmax4, &
        SLL_HERMITE)! sim%bc_vy_0)
 


    print*, ' ... finished initialization, entering main loop.'
    
    !--> Compute energy kinetic, potential and total
    !print*, 'compute nrj'
    !call compute_energy_qns(sim,phi)
    
    !--> Compute L1 norm, L2 norm, L infini norm
    !print*, 'compute  L^p'
    !call compute_norm_L1_L2_Linf_qns(sim)
    
    
    !call writeHDF5_diag_qns( sim )
    ! ------------------------------------------------------------------------
    !
    !                                MAIN LOOP
    !
    ! ------------------------------------------------------------------------
    !print*, 'sim%my_rank = ',sim%my_rank, ' integral x3x4= ', sum(sim%f_x3x4(:,:,:,:))*delta1*delta2*delta3*delta4
    
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




       call apply_remap_4D( sim%seqx1x2_to_seqx3x4, sim%f_x1x2, sim%f_x3x4 )
       
       call compute_local_sizes_4d( &
            sim%sequential_x1x2, &
            loc_sz_x1,           &
            loc_sz_x2,           &
            loc_sz_x3,           &
            loc_sz_x4 )
       
       sim%rho_split(:,:) = 0.0_f64
       call compute_charge_density( &
            sim%mesh2d_x,           &
            sim%mesh2d_v,           &
            size(sim%f_x3x4,1),     &
            size(sim%f_x3x4,2),     &
            sim%f_x3x4,             &
            sim%partial_reduction,  &
            sim%rho_split )
       
       global_indices(1:2) =  &
            local_to_global_2D( sim%split_rho_layout, (/1, 1/) )
       
       call sll_gnuplot_rect_2d_parallel( &
          sim%mesh2d_x%eta1_min+(global_indices(1)-1)*sim%mesh2d_x%delta_eta1, &
          sim%mesh2d_x%delta_eta1, &
          sim%mesh2d_x%eta2_min+(global_indices(2)-1)*sim%mesh2d_x%delta_eta2, &
          sim%mesh2d_x%delta_eta2, &
          size(sim%rho_split,1), &
          size(sim%rho_split,2), &
          sim%rho_split, &
          "rho_split", &
          itime, &
          ierr )
       
       call load_buffer( sim%split_rho_layout, sim%rho_split, send_buf )
       
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
               size(sim%rho_full,1), &
               size(sim%rho_full,2), &
               sim%rho_full, &
               "rho_full_check", &
               itime, &
               ierr )
       end if
    
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
       
!!$       call compute_average_f( &
!!$            sim,&
!!$            sim%mesh2d_x,&
!!$            sim%rho_full, &
!!$            density_tot )
     
       ! print*, 'density', density_tot
       ! The subtraction of density_tot is supposed to be made inside the 
       ! elliptic solver.
       !
!       call rho%update_interpolation_coefficients(sim%rho_full-density_tot)

       call rho%set_field_data(sim%rho_full)
       call rho%update_interpolation_coefficients( )

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
       
       !call sll_set_time_mark(t0)  
       call solve_general_coordinates_elliptic_eq( &
            sim%qns, &
            rho, &
            phi )
       !time = sll_time_elapsed_since(t0)
     
       !print*, 'timer to solve QNS =', time
       if(sim%my_rank == 0) then
          call phi%write_to_file(itime)
       end if
       
       call compute_local_sizes_4d( sim%sequential_x1x2, &
            loc_sz_x1,           &
            loc_sz_x2,           &
            loc_sz_x3,           &
            loc_sz_x4 )
       
       call compute_local_sizes_4d( sim%sequential_x3x4, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 
       
       efield_energy_total = 0.0_f64
       !print*, 'advection vx'
       ! Start with dt in vx...(x3)
       do l=1,loc_sz_x4 !sim%mesh2d_v%num_cells2+1
          do j=1,loc_sz_x2
             do i=1,loc_sz_x1
                global_indices(1:2) = &
                     local_to_global_2D( sim%split_rho_layout, (/i,j/))
                !eta1   =  eta1_min + real(global_indices(1)-1,f64)*delta1
                eta1 = sim%pt_eta1(global_indices(1))
                !eta2   =  eta2_min + real(global_indices(2)-1,f64)*delta2
                eta2 = sim%pt_eta2(global_indices(2))
                x    = sim%point_x(global_indices(1),global_indices(2))
                y    = sim%point_y(global_indices(1),global_indices(2))
                inv_j  =  sim%values_jacobian_matinv(global_indices(1),global_indices(2),:,:)
                !sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                jac_m  =  sim%values_jacobian_mat(global_indices(1),global_indices(2),:,:)
                ! sim%transfx%jacobian_matrix(eta1,eta2)
                
                ex     =  - phi%first_deriv_eta1_value_at_point(eta1,eta2)
                ey     =  - phi%first_deriv_eta2_value_at_point(eta1,eta2)
                sim%values_ex (global_indices(1),global_indices(2) ) = ex
                sim%values_ey (global_indices(1),global_indices(2) ) = ey
                
                alpha3 = -sim%dt*(inv_j(1,1)*ex + inv_j(2,1)*ey)
                alpha3 = alpha3 -sim%dt*(elec_field_ext_1%value_at_point(x,y) )
                sim%f_x3x4(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                     nc_x3+1, &
                     sim%f_x3x4(i,j,:,l), &
                     alpha3 )
                efield_energy_total = efield_energy_total + &
                     delta1*delta2*&
                     abs(sim%values_jacobian(global_indices(1),global_indices(2))) &
                     *abs( ( inv_j(1,1) *inv_j(1,1) + inv_j(1,2)*inv_j(1,2))*ex**2 &
                     +2* ( inv_j(1,1) *inv_j(2,1) + inv_j(1,2)*inv_j(2,2))*abs(ex)*abs(ey) &
                     + ( inv_j(2,1)*inv_j(2,1)  + inv_j(2,2)*inv_j(2,2))*ey**2)
                
             end do
          end do
       end do
       
       !print*, 'energy total', efield_energy_total
       !efield_energy_total = sqrt(efield_energy_total)
       
       call compute_local_sizes_4d( sim%sequential_x3x4, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 
       numpart = 0.0_f64
       ! dt in vy...(x4)
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             global_indices(1:2) = local_to_global_2D( sim%split_rho_layout, (/i,j/))
             do k=1,sim%mesh2d_v%num_cells1+1
                !eta1   =  eta1_min + real(global_indices(1)-1,f64)*delta1
                !eta2   =  eta2_min + real(global_indices(2)-1,f64)*delta2
                eta1 = sim%pt_eta1(global_indices(1))
                eta2 = sim%pt_eta2(global_indices(2))
                x    = sim%point_x(global_indices(1),global_indices(2))
                y    = sim%point_y(global_indices(1),global_indices(2))
                !inv_j  =  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                inv_j  =  sim%values_jacobian_matinv(global_indices(1),global_indices(2),:,:)
                !jac_m  =  sim%transfx%jacobian_matrix(eta1,eta2)
                jac_m  =  sim%values_jacobian_mat(global_indices(1),global_indices(2),:,:)
                ex     =  sim%values_ex (global_indices(1),global_indices(2) ) 
                !- phi%first_deriv_eta1_value_at_point(eta1,eta2)
                ey     =  sim%values_ey (global_indices(1),global_indices(2) ) !
                !- phi%first_deriv_eta2_value_at_point(eta1,eta2)
                alpha4 = -sim%dt*(inv_j(1,2)*ex + inv_j(2,2)*ey)
                alpha4 = alpha4 -sim%dt*(elec_field_ext_2%value_at_point(x,y))
                sim%f_x3x4(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
                     nc_x4+1, &
                     sim%f_x3x4(i,j,k,:), &
                     alpha4 )

                
             end do
             numpart = numpart + sum(sim%f_x3x4(i,j,:,:))*&
                  abs(sim%values_jacobian(global_indices(1),global_indices(2)))
!!$             efield_energy_total = efield_energy_total + &
!!$                  delta1*delta2 *abs(jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1)) &
!!$                  *abs( ( inv_j(1,1) *inv_j(1,1) + inv_j(1,2)*inv_j(1,2))*ex**2 &
!!$                  +2* ( inv_j(1,1) *inv_j(2,1) + inv_j(1,2)*inv_j(2,2))*abs(ex)*abs(ey) &
!!$                  + ( inv_j(2,1)*inv_j(2,1)  + inv_j(2,2)*inv_j(2,2))*ey**2)
          end do
       end do

       !print*, 'sim%my_rank = ', sim%my_rank, 'numpart = ', numpart*delta1*delta2*delta3*delta4

       call apply_remap_4D( sim%seqx3x4_to_seqx1x2, sim%f_x3x4, sim%f_x1x2 )
       
       call compute_local_sizes_4d( sim%sequential_x1x2, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 
       
       ! Approximate the integral of the distribution function along all
       ! directions.
       num_particles_local(buffer_counter) = numpart*delta1*delta2*delta3*delta4
       buffer_energy(buffer_counter)       = efield_energy_total
           ! sum(sim%f_x3x4)*delta1*delta2*delta3*delta4
       if( buffer_counter == BUFFER_SIZE ) then
          ! ----------------------------------
          ! write particles buffer to disk 
          ! ----------------------------------
          num_particles_global = 0.0_f64
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
          
          ! ----------------------------------
          ! write electric field energy to disk 
          ! ----------------------------------
          buffer_energy_result = 0.0_f64
          call sll_collective_reduce_real64( &
               sll_world_collective, &
               buffer_energy, &
               BUFFER_SIZE, &
               MPI_SUM, &
               0, &
               buffer_energy_result )
          

           if(sim%my_rank == 0) then
             open(efield_energy_file_id,file="electric_field_energy_qns",&
                  position="append")
             
             if(itime == BUFFER_SIZE) then
                rewind(efield_energy_file_id)
             end if
             buffer_energy_result(:) = log(sqrt(buffer_energy_result(:)))
             do i=1,BUFFER_SIZE
                write(efield_energy_file_id,*) buffer_energy_result(i)
             end do
             close(efield_energy_file_id)
          end if
          buffer_counter = 1

       else
          buffer_counter         = buffer_counter + 1
       end if
       
       
       

!!$       if (sim%my_rank == 0) then
!!$          
!!$          call sll_new_file_id(droite_test_pente, ierr) 
!!$          open(droite_test_pente,file="droite_test_pente",&
!!$               position="append")
!!$          write(droite_test_pente,*) -0.1533*(itime-1)*sim%dt + 0.676!1.676 ! the test case 2002
!!$          close(droite_test_pente)
!!$       end if

      
       ! Proceed to the advections in the spatial directions, 'x' and 'y'
       ! Reconfigure data. 
       
       ! what are the new local limits on x3 and x4? It is bothersome to have
       ! to make these calls...

       call advection_x1x2(sim,sim%dt) 

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
       
       
      
       !--> Save results in HDF5 files
       if ( mod(itime,BUFFER_SIZE) == 0) then
          !--> Compute energy kinetic, potential and total
          !print*, 'compute nrj'
          call compute_energy_qns(sim,phi)
          
          !--> Compute L1 norm, L2 norm, L infini norm
          !print*, 'compute  L^p'
          call compute_norm_L1_L2_Linf_qns(sim)
          
          call writeHDF5_diag_qns( sim )
          
       end if
       
    end do ! main loop
#undef BUFFER_SIZE
#endif
  end subroutine run_4d_qns_general_mp
  

  subroutine advection_x1x2(sim,deltat)
    class(sll_simulation_4d_qns_general_multipatch) :: sim
    sll_real64, intent(in) :: deltat
    sll_int32 :: gi, gj, gk, gl
    sll_real64, dimension(1:2,1:2) :: inv_j,jac_m
    sll_int32 :: loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 
    sll_int32, dimension(4) :: global_indices
    sll_real64 :: alpha1, alpha2
    sll_real64 :: eta1, eta2, eta3, eta4
    sll_int32 :: i, j, k, l

    call compute_local_sizes_4d( sim%sequential_x1x2, &
         loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 )

    do l=1,loc_sz_x4
       do k=1,loc_sz_x3
          ! call sim%interp_x1x2%compute_interpolants(sim%f_x1x2(:,:,k,l))
          
          do j=1,loc_sz_x2
             do i=1,loc_sz_x1
                global_indices = &
                     local_to_global_4D(sim%sequential_x1x2,(/i,j,k,l/))
                gi = global_indices(1)
                gj = global_indices(2)
                gk = global_indices(3)
                gl = global_indices(4)
!!$                eta1 = sim%mesh2d_x%eta1_min + (gi-1)*sim%mesh2d_x%delta_eta1
!!$                eta2 = sim%mesh2d_x%eta2_min + (gj-1)*sim%mesh2d_x%delta_eta2
                eta1 = sim%pt_eta1(gi)
                eta2 = sim%pt_eta2(gj)
                eta3 = sim%mesh2d_v%eta1_min + (gk-1)*sim%mesh2d_v%delta_eta1
                eta4 = sim%mesh2d_v%eta2_min + (gl-1)*sim%mesh2d_v%delta_eta2
                ! inv_j  =  sim%values_jacobian_matinv(gi,gj,:,:) 
                !sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                !jac_m  =  sim%transfx%jacobian_matrix(eta1,eta2)
                alpha1 = -deltat*(inv_j(1,1)*eta3 + inv_j(1,2)*eta4)
                alpha2 = -deltat*(inv_j(2,1)*eta3 + inv_j(2,2)*eta4)
                
                eta1 = eta1+alpha1
                eta2 = eta2+alpha2
                ! This is hardwiring the periodic BC, please improve this..
                if (( sim%bc_eta1_0 == SLL_PERIODIC) .and.&
                     (sim%bc_eta1_1 == SLL_PERIODIC) .and.&
                     (sim%bc_eta2_0 == SLL_PERIODIC) .and.&
                     (sim%bc_eta2_1 == SLL_PERIODIC)) then
                   
                   ! PERIODIC TEST CASE.
                   ! if( eta1 <  sim%mesh2d_x%eta1_min ) then
                   !    eta1 = eta1+sim%mesh2d_x%eta1_max-sim%mesh2d_x%eta1_min
                   ! else if( eta1 >  sim%mesh2d_x%eta1_max ) then
                   !    eta1 = eta1+sim%mesh2d_x%eta1_min-sim%mesh2d_x%eta1_max
                   ! end if
                   ! if( eta2 <  sim%mesh2d_x%eta2_min ) then
                   !    eta2 = eta2+sim%mesh2d_x%eta2_max-sim%mesh2d_x%eta2_min
                   ! else if( eta2 >  sim%mesh2d_x%eta2_max ) then
                   !    eta2 = eta2+sim%mesh2d_x%eta2_min-sim%mesh2d_x%eta2_max
                   ! end if

                   ! sim%f_x1x2(i,j,k,l) = sim%interp_x1x2%interpolate_value(eta1,eta2)

                else if (( sim%bc_eta1_0 == SLL_PERIODIC) .and.&
                     (sim%bc_eta1_1 == SLL_PERIODIC) .and.&
                     (sim%bc_eta2_0 == SLL_DIRICHLET) .and.&
                     (sim%bc_eta2_1 == SLL_DIRICHLET)) then

                   ! if( eta1 <  sim%mesh2d_x%eta1_min ) then
                   !    eta1 = eta1+sim%mesh2d_x%eta1_max-sim%mesh2d_x%eta1_min
                   ! else if( eta1 >  sim%mesh2d_x%eta1_max ) then
                   !    eta1 = eta1+sim%mesh2d_x%eta1_min-sim%mesh2d_x%eta1_max
                   ! end if
                   ! if( eta2 <  sim%mesh2d_x%eta2_min ) then
                   !    !eta2 = eta2+sim%mesh2d_x%eta2_max-sim%mesh2d_x%eta2_min
                   !    sim%f_x1x2(i,j,k,l) = 0.0_f64
                   ! else if( eta2 >  sim%mesh2d_x%eta2_max ) then
                   !    !eta2 = eta2+sim%mesh2d_x%eta2_min-sim%mesh2d_x%eta2_max
                   !    sim%f_x1x2(i,j,k,l) = 0.0_f64
                   ! else
                   !    sim%f_x1x2(i,j,k,l) = sim%interp_x1x2%interpolate_value(eta1,eta2)
                   ! end if
                   
                else if (( sim%bc_eta1_0 == SLL_DIRICHLET) .and.&
                     (sim%bc_eta1_1 == SLL_DIRICHLET) .and.&
                     (sim%bc_eta2_0 == SLL_PERIODIC) .and.&
                     (sim%bc_eta2_1 == SLL_PERIODIC)) then

                   ! if( eta2 <  sim%mesh2d_x%eta2_min ) then
                   !    eta2 = eta2+sim%mesh2d_x%eta2_max-sim%mesh2d_x%eta2_min
                   ! else if( eta2 >  sim%mesh2d_x%eta2_max ) then
                   !    eta2 = eta2+sim%mesh2d_x%eta2_min-sim%mesh2d_x%eta2_max
                   ! end if
                   ! if( eta1 <  sim%mesh2d_x%eta1_min ) then
                   !    !eta1 = eta1+sim%mesh2d_x%eta1_max-sim%mesh2d_x%eta1_min
                   !    sim%f_x1x2(i,j,k,l) = 0.0_f64
                   ! else if( eta1 >  sim%mesh2d_x%eta1_max ) then
                   !    !eta1 = eta1+sim%mesh2d_x%eta1_min-sim%mesh2d_x%eta1_max
                   !    sim%f_x1x2(i,j,k,l) = 0.0_f64
                   ! else
                   !    sim%f_x1x2(i,j,k,l) = sim%interp_x1x2%interpolate_value(eta1,eta2)
                   ! end if
                else if (( sim%bc_eta1_0 == SLL_DIRICHLET) .and.&
                     (sim%bc_eta1_1 == SLL_DIRICHLET) .and.&
                     (sim%bc_eta2_0 == SLL_DIRICHLET) .and.&
                     (sim%bc_eta2_1 == SLL_DIRICHLET)) then
                   
                   ! if(  (eta2 <  sim%mesh2d_x%eta2_min ) .or.&
                   !      (eta2 >  sim%mesh2d_x%eta2_max ) .or.&
                   !      (eta1 <  sim%mesh2d_x%eta1_min ) .or.&
                   !      (eta1 >  sim%mesh2d_x%eta1_max )) then
                   !    sim%f_x1x2(i,j,k,l) = 0.0_f64
                   ! else
                   !    sim%f_x1x2(i,j,k,l) = sim%interp_x1x2%interpolate_value(eta1,eta2)
                   ! end if
                   !sim%f_x1x2(1,:,k,l) = 0.0_f64
                   !sim%f_x1x2(:,1,k,l) = 0.0_f64
                   !sim%f_x1x2(sim%mesh2d_x%num_cells1 +1,:,k,l) = 0.0_f64
                   !sim%f_x1x2(:,sim%mesh2d_x%num_cells2 +1,k,l) = 0.0_f64
                else 
                   print*, 'problem boundary conditon for particles motions in eta2'
                   stop
                end if
                !print*,sim%f_x1x2(i,j,k,l),gk,gl,gi,gj,alpha1,alpha2,eta1,eta2,eta3,eta4

                
                !print*, 'before', eta1,eta2,i,j
                !sim%f_x1x2(i,j,k,l) = sim%interp_x1x2%interpolate_value(eta1,eta2)
                
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
          !print*, sim%f_x1x2(1,1,k,l) 
       end do
    end do
   
 

  end subroutine advection_x1x2

  subroutine advection_x3(sim,phi,deltat,efield_energy_total)
    class(sll_simulation_4d_qns_general_multipatch) :: sim
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
             
             ! !FAUX
             ! inv_j  =  sim%values_jacobian_matinv(global_indices(1),global_indices(2),:,:)
             ! !sim%transfx%inverse_jacobian_matrix(eta1,eta2)
             ! jac_m  =  sim%values_jacobian_mat(global_indices(1),global_indices(2),:,:)
             ! !sim%transfx%jacobian_matrix(eta1,eta2)
             
             ex     =  - phi%first_deriv_eta1_value_at_point(eta1,eta2)
             ey     =  - phi%first_deriv_eta2_value_at_point(eta1,eta2)
             
             alpha3 = -sim%dt*(inv_j(1,1)*ex + inv_j(2,1)*ey)
             sim%f_x3x4(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                  nc_x3+1, &
                  sim%f_x3x4(i,j,:,l), &
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
    class(sll_simulation_4d_qns_general_multipatch) :: sim
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
             ! !FAUX
             ! inv_j  = sim%values_jacobian_matinv(global_indices(1),global_indices(2),:,:)  
             !sim%transfx%inverse_jacobian_matrix(eta1,eta2)
             
             ex     =  - phi%first_deriv_eta1_value_at_point(eta1,eta2)
             ey     =  - phi%first_deriv_eta2_value_at_point(eta1,eta2)
             alpha4 = -sim%dt*(inv_j(1,2)*ex + inv_j(2,2)*ey)
             sim%f_x3x4(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
                  nc_x4+1, &
                  sim%f_x3x4(i,j,k,:), &
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


  subroutine delete_4d_qns_gen_mp( sim )
    type(sll_simulation_4d_qns_general_multipatch) :: sim
    sll_int32 :: ierr
    SLL_DEALLOCATE( sim%f_x1x2, ierr )
    SLL_DEALLOCATE( sim%f_x3x4, ierr )
    SLL_DEALLOCATE_ARRAY( sim%partial_reduction, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_full, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_split, ierr )
    ! SLL_DEALLOCATE(sim%values_ex ,ierr)
    ! SLL_DEALLOCATE(sim%values_ey ,ierr)
   ! SLL_DEALLOCATE_ARRAY( sim%phi_x1, ierr )
   ! SLL_DEALLOCATE_ARRAY( sim%phi_x2, ierr )
   ! SLL_DEALLOCATE_ARRAY( sim%phi_split, ierr )
    call sll_delete( sim%sequential_x1x2 )
    call sll_delete( sim%sequential_x3x4 )
    call sll_delete( sim%rho_full_layout )
    call sll_delete( sim%rho_seq_x2 )
    call sll_delete( sim%split_rho_layout )
    call sll_delete( sim%split_to_full )
    call sll_delete( sim%efld_seqx1_to_seqx2 )
    call sll_delete( sim%efld_seqx2_to_split )
    call sll_delete( sim%seqx1x2_to_seqx3x4 )
    call sll_delete( sim%seqx3x4_to_seqx1x2 )
    ! call sll_delete( sim%interp_x1x2 )
    call delete( sim%interp_x3 )
    call delete( sim%interp_x4 )
    SLL_DEALLOCATE(sim%diag_masse,ierr)
    SLL_DEALLOCATE(sim%diag_norm_L1,ierr)
    SLL_DEALLOCATE(sim%diag_norm_L2,ierr)
    SLL_DEALLOCATE(sim%diag_norm_Linf,ierr)
    SLL_DEALLOCATE(sim%diag_entropy_kin,ierr)

    !---> DEALLOCATE mat fort jacobian

    ! SLL_DEALLOCATE(sim%values_jacobian_mat,ierr)
    ! SLL_DEALLOCATE(sim%values_jacobian_matinv,ierr)
    ! SLL_DEALLOCATE(sim%values_jacobian,ierr)
    ! SLL_DEALLOCATE(sim%point_x,ierr)
    ! SLL_DEALLOCATE(sim%point_y,ierr)
    ! ---> DEALLOCATE array 1D contains mesh points

    SLL_DEALLOCATE(sim%pt_eta1,ierr)
    SLL_DEALLOCATE(sim%pt_eta2,ierr)

    !--> DEALLOCATE diagnostics for the energy
    SLL_DEALLOCATE(sim%diag_nrj_kin,ierr)
    SLL_DEALLOCATE(sim%diag_nrj_pot,ierr)
    SLL_DEALLOCATE(sim%diag_nrj_tot,ierr)
    !SLL_DEALLOCATE_ARRAY( sim%efield_x1, ierr )
    !SLL_DEALLOCATE_ARRAY( sim%efield_x2, ierr )
    !SLL_DEALLOCATE_ARRAY( sim%efield_split, ierr )
  end subroutine delete_4d_qns_gen_mp

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


   subroutine compute_average_f( &
       sim,&
       mx,&
       rho, &
       density_tot )
    class(sll_simulation_4d_qns_general_multipatch)     :: sim
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
          !eta1 = sim%mesh2d_x%eta1_min + (i-1)*sim%mesh2d_x%delta_eta1
          eta1 = sim%pt_eta1(i)
          eta2 = sim%pt_eta2(j)
          !eta2 = sim%mesh2d_x%eta2_min + (j-1)*sim%mesh2d_x%delta_eta2
          ! jac_m  = sim%values_jacobian_mat(i,j,:,:)
          !sim%transfx%jacobian_matrix(eta1,eta2)
          
          ! density_tot = density_tot + rho(i,j)*delta1*delta2*&
          !      (sim%values_jacobian(i,j))
          !print*, jac_m(1,1)*jac_m(2,2)-jac_m(1,2)*jac_m(2,1)
          ! length_total = length_total + &
          !      delta1*delta2*(sim%values_jacobian(i,j))
       end do
    end do
    !print*, length_total
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



 !----------------------------------------------------
  subroutine writeHDF5_diag_qns( sim )
   ! use sll_collective
    use sll_hdf5_io_serial, only: sll_hdf5_file_create, &
      sll_hdf5_write_array_1d, sll_hdf5_file_close
    class(sll_simulation_4d_qns_general_multipatch), intent(inout) :: sim

    sll_int32 :: ix1_diag, ix2_diag
    sll_int32 :: iv1_diag, iv2_diag

    !--> diagnostics norm
    sll_real64, dimension(sim%count_save_diag + 1) :: diag_masse_result
    sll_real64, dimension(sim%count_save_diag + 1) :: diag_norm_L1_result
    sll_real64, dimension(sim%count_save_diag + 1) :: diag_norm_L2_result
    sll_real64, dimension(sim%count_save_diag + 1) :: diag_norm_Linf_result
    sll_real64, dimension(sim%count_save_diag + 1) :: diag_entropy_kin_result
    
    !--> diagnostics energy
    sll_real64, dimension(sim%count_save_diag + 1) :: diag_nrj_kin_result
    sll_real64, dimension(sim%count_save_diag + 1) :: diag_nrj_pot_result
    sll_real64, dimension(sim%count_save_diag + 1) :: diag_nrj_tot_result
    sll_real64, dimension(sim%count_save_diag + 1) :: diag_relative_error_nrj_tot_result
    
    !--> For initial profile HDF5 saving
    integer             :: file_err
    sll_int32           :: file_id
    character(len=80)   :: filename_HDF5
    character(20), save :: numfmt = "'_d',i5.5"
    
    ix1_diag   = int(sim%nc_x1/2)
    ix2_diag   = int(sim%nc_x2/3)
    iv1_diag   = int(sim%nc_x3/4)
    iv2_diag   = int(sim%nc_x4/3)

    diag_masse_result       = 0.0_f64
    diag_norm_L1_result     = 0.0_f64
    diag_norm_L2_result     = 0.0_f64
    diag_norm_Linf_result   = 0.0_f64
    diag_entropy_kin_result = 0.0_f64
    diag_nrj_kin_result     = 0.0_f64
    diag_nrj_pot_result     = 0.0_f64
    diag_nrj_tot_result     = 0.0_f64
    diag_relative_error_nrj_tot_result = 0.0_f64

   ! print*, 'sim%count_save_diag = ', sim%count_save_diag
 
   call sll_collective_reduce_real64( &
               sll_world_collective, &
               sim%diag_masse(1:sim%count_save_diag + 1), &
               sim%count_save_diag + 1, &
               MPI_SUM, &
               0, &
               diag_masse_result )

    call sll_collective_reduce_real64( &
               sll_world_collective, &
               sim%diag_norm_L1(1:sim%count_save_diag + 1), &
               sim%count_save_diag + 1, &
               MPI_SUM, &
               0, &
               diag_norm_L1_result )

    call sll_collective_reduce_real64( &
               sll_world_collective, &
               sim%diag_norm_L2(1:sim%count_save_diag + 1), &
               sim%count_save_diag + 1, &
               MPI_SUM, &
               0, &
               diag_norm_L2_result )

    call sll_collective_reduce_real64( &
               sll_world_collective, &
               sim%diag_norm_Linf(1:sim%count_save_diag + 1), &
               sim%count_save_diag + 1, &
               MPI_SUM, &
               0, &
               diag_norm_Linf_result )

    call sll_collective_reduce_real64( &
               sll_world_collective, &
               sim%diag_entropy_kin(1:sim%count_save_diag + 1), &
               sim%count_save_diag + 1, &
               MPI_SUM, &
               0, &
               diag_entropy_kin_result )

    call sll_collective_reduce_real64( &
               sll_world_collective, &
               sim%diag_nrj_kin(1:sim%count_save_diag + 1), &
               sim%count_save_diag + 1, &
               MPI_SUM, &
               0, &
               diag_nrj_kin_result )

    
    call sll_collective_reduce_real64( &
               sll_world_collective, &
               sim%diag_nrj_pot(1:sim%count_save_diag + 1), &
               sim%count_save_diag + 1, &
               MPI_SUM, &
               0, &
               diag_nrj_pot_result )

    call sll_collective_reduce_real64( &
               sll_world_collective, &
               sim%diag_nrj_tot(1:sim%count_save_diag + 1), &
               sim%count_save_diag + 1, &
               MPI_SUM, &
               0, &
               diag_nrj_tot_result )

    diag_relative_error_nrj_tot_result(:) = &
         (diag_nrj_tot_result(:)-diag_nrj_tot_result(1))/&
         sqrt(0.5*( (diag_nrj_kin_result(:)-diag_nrj_kin_result(1) )**2 + &
                    (diag_nrj_pot_result(:)-diag_nrj_pot_result(1) )**2 ) ) 

    write(filename_HDF5,'(A,'//numfmt//',A)') &
      "vp4D_diag", sim%count_save_diag, ".h5"

    if (sim%my_rank.eq.0) then
      print*,'--> Save HDF5 file: ',filename_HDF5
      call sll_hdf5_file_create(filename_HDF5,file_id,file_err)
      call sll_hdf5_write_array_2d(file_id, &
        sim%f_x1x2(:,:,iv1_diag,iv2_diag),'f2d_xy',file_err)
      call sll_hdf5_write_array_2d(file_id, &
        sim%f_x3x4(ix1_diag,ix2_diag,:,:),'f2d_v1v2',file_err)
      call sll_hdf5_write_array_1d(file_id,&
                          diag_nrj_kin_result(:),'nrj_kin',file_err)
      call sll_hdf5_write_array_1d(file_id,&
                          diag_nrj_pot_result(:),'nrj_pot',file_err)
      call sll_hdf5_write_array_1d(file_id,&
                          diag_nrj_tot_result(:),'nrj_tot',file_err)
      call sll_hdf5_write_array_1d(file_id,&
                          diag_relative_error_nrj_tot_result(:),&
                          'relative_error_nrj_tot',file_err)
      call sll_hdf5_write_array_1d(file_id,&
                          diag_masse_result(:),&
                          'masse',file_err)
      call sll_hdf5_write_array_1d(file_id,&
                          diag_norm_L1_result(:),'L1_norm',file_err)
      call sll_hdf5_write_array_1d(file_id,&
                          diag_norm_L2_result(:),'L2_norm',file_err)
      call sll_hdf5_write_array_1d(file_id,&
                          diag_norm_Linf_result(:),'Linf_norm',file_err)
      call sll_hdf5_write_array_1d(file_id,&
                          diag_entropy_kin_result(:),'entropy_kin',file_err)
      ! call sll_hdf5_write_array_2d(file_id,sim%point_x(:,:),'X_coord',file_err)
      ! call sll_hdf5_write_array_2d(file_id,sim%point_y(:,:),'Y_coord',file_err)
      call sll_hdf5_file_close(file_id,file_err)
      
    end if
    sim%count_save_diag = sim%count_save_diag + 1
    
  end subroutine writeHDF5_diag_qns



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
  
  subroutine compute_energy_qns(sim,phi)

    class(sll_simulation_4d_qns_general_multipatch), intent(inout) :: sim
    type(sll_scalar_field_2d_discrete_alt), pointer       :: phi

    ! local variables
    sll_int32  :: Neta1_loc,Neta2_loc,Nv2,Nv1,Neta1,Neta2
    sll_real64 :: delta_eta1,delta_eta2,delta_v1,delta_v2
    sll_real64, dimension(1:2,1:2) :: inv_j
    sll_real64 :: val_jac
    sll_real64 :: v1,v2
    sll_real64 :: ex,ey
    sll_real64 :: eta1,eta2
    sll_real64 :: delta_f
    sll_int32  :: iloc1, iloc2
    sll_int32  :: i1,i2,iv1,iv2
    sll_int32, dimension(1:4) :: glob_ind4d
    sll_real64 :: nrj_kin,nrj_pot,nrj_tot
    
    Neta1_loc  = size(sim%f_x3x4,1)
    Neta2_loc  = size(sim%f_x3x4,2)
    Neta1      = sim%nc_x1 + 1
    Neta2      = sim%nc_x2 + 1
    Nv1        = size(sim%f_x3x4,3)
    Nv2        = size(sim%f_x3x4,4)
    ! delta_eta1 = sim%mesh2d_x%delta_eta1
    ! delta_eta2 = sim%mesh2d_x%delta_eta2
    delta_v1   = sim%mesh2d_v%delta_eta1
    delta_v2   = sim%mesh2d_v%delta_eta2
    
    
    nrj_kin   = 0.0
    nrj_pot   = 0.0
    nrj_tot   = 0.0
    
    !-> Computation of the energy kinetic locally in (x1,x2) directions
    do iloc2 = 1,Neta2_loc
       do iloc1 = 1,Neta1_loc
          do iv1 = 1,Nv1-1
             v1 = sim%mesh2d_v%eta1_min + (iv1-1)*delta_eta1
                
             do iv2 = 1,Nv2-1
                v2 = sim%mesh2d_v%eta2_min + (iv2-1)*delta_eta2
                
                glob_ind4d(:) = local_to_global_4D(sim%sequential_x3x4, &
                     (/iloc1,iloc2,iv1,iv2/))
                i1 = glob_ind4d(1)
                i2 = glob_ind4d(2)

                !eta1   =  sim%mesh2d_x%eta1_min + real(i1-1,f64)*delta_eta1
                !eta2   =  sim%mesh2d_x%eta2_min + real(i2-1,f64)*delta_eta2
                eta1 = sim%pt_eta1(i1)
                eta2 = sim%pt_eta2(i2)
                ! ex     =  sim%values_ex (i1,i2 )!
                !- phi%first_deriv_eta1_value_at_point(eta1,eta2)
                ! ey     =  sim%values_ey(i1,i2)
                !- phi%first_deriv_eta2_value_at_point(eta1,eta2)
                ! inv_j  =  sim%values_jacobian_matinv(i1,i2,:,:)
                !sim%transfx%inverse_jacobian_matrix(eta1,eta2)
                
                if (i1 .ne. Neta1) then
                   if (i2 .ne. Neta2) then 
                      
                      ! val_jac = abs(sim%values_jacobian(i1,i2))
                      !abs(sim%transfx%jacobian_at_node(i1,i2))
                      
                      delta_f = sim%f_x3x4(iloc1,iloc2,iv1,iv2)
                      
                      
                      nrj_kin = nrj_kin + &
                           delta_f * (v1**2 + v2**2) * 0.5 * val_jac * &
                           delta_eta1*delta_eta2*delta_v1*delta_v2
                      
                   end if
                end if
             end do
          end do
          
          nrj_pot = nrj_pot + &
               ((inv_j(1,1)*ex + inv_j(2,1)*ey)**2+ &
               (inv_j(1,2)*ex + inv_j(2,2)*ey)**2) * 0.5 * val_jac * &
               delta_eta1*delta_eta2
       end do
    end do
    
    nrj_tot = nrj_kin + nrj_pot
    
    sim%diag_nrj_kin(sim%count_save_diag + 1)   = nrj_kin
    sim%diag_nrj_pot(sim%count_save_diag + 1)   = nrj_pot
    sim%diag_nrj_tot(sim%count_save_diag + 1)   = nrj_tot
  end subroutine compute_energy_qns
  
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
  
  subroutine compute_norm_L1_L2_Linf_qns(sim)
    
    class(sll_simulation_4d_qns_general_multipatch), intent(inout) :: sim
    
    ! local variables
    sll_int32  :: Neta1_loc,Neta2_loc,Nv1,Nv2
    sll_int32  :: Neta1, Neta2
    sll_real64 :: delta_eta1,delta_eta2,delta_v1,delta_v2
    sll_real64 :: val_jac
    sll_real64 :: delta_f
    sll_int32  :: iloc1, iloc2
    sll_int32  :: i1,i2,iv1,iv2
    sll_int32, dimension(1:4) :: glob_ind4d
    sll_real64 :: masse, norm_L1,norm_L2,norm_Linf,entropy_kin

    Neta1_loc  = size(sim%f_x3x4,1)
    Neta2_loc  = size(sim%f_x3x4,2)
    Neta1      = sim%nc_x1 + 1
    Neta2      = sim%nc_x2 + 1
    Nv1        = size(sim%f_x3x4,3)
    Nv2        = size(sim%f_x3x4,4)
    ! delta_eta1 = sim%mesh2d_x%delta_eta1
    ! delta_eta2 = sim%mesh2d_x%delta_eta2
    delta_v1   = sim%mesh2d_v%delta_eta1
    delta_v2   = sim%mesh2d_v%delta_eta2

    
    norm_L1    = 0.0
    norm_L2    = 0.0
    norm_Linf  = 0.0
    masse      = 0.0
    entropy_kin= 0.0

    !-> Computation of the enrgy kinetic locally in (x1,x2) directions
    do iloc2 = 1,Neta2_loc
       do iloc1 = 1,Neta1_loc
          do iv1 = 1,Nv1-1
             do iv2 = 1,Nv2-1
                
                glob_ind4d(:) = local_to_global_4D(sim%sequential_x3x4, &
                     (/iloc1,iloc2,iv1,iv2/))
                i1 = glob_ind4d(1)
                i2 = glob_ind4d(2)
                
                if (i1 .ne. Neta1) then
                   if (i2 .ne. Neta2) then 
                      
                      ! val_jac = abs(sim%values_jacobian(i1,i2))
                      !abs(sim%transfx%jacobian_at_node(i1,i2))
                      
                      delta_f = sim%f_x3x4(iloc1,iloc2,iv1,iv2)
                      
                      masse   = masse + &
                           delta_f * val_jac * &
                           delta_eta1*delta_eta2*delta_v1*delta_v2
                      
                      norm_L1 = norm_L1 + &
                           abs(delta_f) * val_jac * &
                           delta_eta1*delta_eta2*delta_v1*delta_v2
                      
                      norm_L2 = norm_L2 + &
                           abs(delta_f)**2 * val_jac * &
                           delta_eta1*delta_eta2*delta_v1*delta_v2
                      
                      entropy_kin = entropy_kin - &
                           delta_f* log(abs(delta_f)) * val_jac * &
                           delta_eta1*delta_eta2*delta_v1*delta_v2
                      
                      norm_Linf = max(abs(delta_f),norm_Linf)
                   end if
                end if
             end do
          end do
       end do
    end do
    !print*, 'norm_L2',norm_L2
    norm_L2   = sqrt(norm_L2)
    !print*, 'norm_L2',norm_L2
    
    sim%diag_masse(sim%count_save_diag + 1)       = masse
    sim%diag_norm_L1(sim%count_save_diag + 1)     = norm_L1
    sim%diag_norm_L2(sim%count_save_diag + 1)     = norm_L2
    sim%diag_norm_Linf(sim%count_save_diag + 1)   = norm_Linf
    sim%diag_entropy_kin(sim%count_save_diag + 1) = entropy_kin
  end subroutine compute_norm_L1_L2_Linf_qns

  subroutine compute_values_jacobian_and_mesh_points(sim)
    class(sll_simulation_4d_qns_general_multipatch), intent(inout) :: sim
    sll_real64 :: delta1,delta2
    sll_int32  :: i,j
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta1,eta2
    
    ! delta1 = sim%mesh2d_x%delta_eta1
    ! delta2 = sim%mesh2d_x%delta_eta2
    
    ! eta1_min = sim%mesh2d_x%eta1_min
    ! eta2_min = sim%mesh2d_x%eta2_min
    
    ! do j=1,sim%mesh2d_x%num_cells2 +1
    !    eta2   =  eta2_min + real(j-1,f64)*delta2
    !    sim%pt_eta2(j) = eta2
    !    do i=1,sim%mesh2d_x%num_cells1 +1
    !       eta1   =  eta1_min + real(i-1,f64)*delta1
    !       sim%pt_eta1(i) = eta1
          
    !       sim%values_jacobian_matinv(i,j,:,:)=  sim%transfx%inverse_jacobian_matrix(eta1,eta2)
    !       sim%values_jacobian_mat(i,j,:,:)   =  sim%transfx%jacobian_matrix(eta1,eta2)
    !       sim%values_jacobian(i,j)           =  sim%values_jacobian_mat(i,j,1,1)*&
    !                                             sim%values_jacobian_mat(i,j,2,2)-&
    !                                             sim%values_jacobian_mat(i,j,1,2)*&
    !                                             sim%values_jacobian_mat(i,j,2,1)
    !       sim%point_x(i,j) = sim%transfx%x1(eta1,eta2)
    !       sim%point_y(i,j) = sim%transfx%x2(eta1,eta2)
    !    end do
    ! end do
    
  end subroutine compute_values_jacobian_and_mesh_points
end module sll_simulation_4d_qns_general_multipatch_module
