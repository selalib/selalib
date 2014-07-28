module sll_simulation_4d_qns_general_multipatch_module

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"

  use sll_logical_meshes
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
  use sll_timer
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
     ! for general coordinate QNS
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
    
  end subroutine initialize_4d_qns_gen_mp


  subroutine init_4d_qns_gen_mp( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_qns_general_multipatch), intent(inout) :: sim
    character(len=*), intent(in)                                   :: filename
    sll_int32             :: IO_stat
    sll_real64            :: dt
    sll_int32             :: number_iterations
    sll_int32, parameter  :: input_file = 99

    namelist /sim_params/ dt, number_iterations
    ! Try to add here other parameters to initialize the mesh values like
    ! xmin, xmax and also for the distribution function initializer.
    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_vp4d_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, sim_params)
!    read(input_file,grid_dims)
    close(input_file)

    sim%dt = dt
    sim%num_iterations = number_iterations
    print*, 'number iterations', number_iterations
  end subroutine init_4d_qns_gen_mp


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

    print *, '******************** ENTERED THE RUN ROUTINE **************'

    ! Start with the fields
    call sll_set_time_mark(t0) 
    a11_field_mat => &
         new_scalar_field_multipatch_2d("a11", sim%transfx, owns_data=.true.)

    a12_field_mat => &
         new_scalar_field_multipatch_2d("a12", sim%transfx, owns_data=.true.)

    a21_field_mat => &
         new_scalar_field_multipatch_2d("a21", sim%transfx, owns_data=.true.)

    a22_field_mat => &
         new_scalar_field_multipatch_2d("a22", sim%transfx, owns_data=.true.)

    b1_field_vect => &
         new_scalar_field_multipatch_2d("b1", sim%transfx, owns_data=.true.)
    
    b2_field_vect => &
         new_scalar_field_multipatch_2d("b2", sim%transfx, owns_data=.true.)
    
    c_field       => &
         new_scalar_field_multipatch_2d("c", sim%transfx, owns_data=.true.)
    
    layer_x1x2    => &
         new_scalar_field_multipatch_2d("layer_x1x2", &
                                        sim%transfx, &
                                        owns_data=.false.)
    
    phi => &
         new_scalar_field_multipatch_2d("potential_field_phi", &
                                        sim%transfx, &
                                        owns_data=.true.)
    
    rho => &
         new_scalar_field_multipatch_2d("rho_field_multipatch", &
                                        sim%transfx, &
                                        owns_data=.true.)
    time = sll_time_elapsed_since(t0)
    print*, 'rank: ', sim%my_rank, 'time to create multipatch fields =', time

    ! elec_field_ext_1 => new_scalar_field_multipatch_2d("E1_ext", sim%transfx)    
    ! elec_field_ext_2 => new_scalar_field_multipatch_2d("E2_ext", sim%transfx)    


    ! call elec_field_ext_1_field_mat%allocate_memory()
    ! call elec_field_ext_2_field_mat%allocate_memory()
    
    
    num_patches = sim%transfx%get_number_patches()

    print *, 'rank:', sim%my_rank, 'initializing patches.'
    call sll_set_time_mark(t0)
    do ipatch= 0,num_patches-1
       ! Please get rid of these 'fixes' whenever it is decided that gfortran 
       ! 4.6 is no longer supported by Selalib.
       !     m        => sim%transfx%get_logical_mesh(ipatch)
       ! logical_m => sim%transfx%transfs(ipatch+1)%t%mesh
       !     transf   => sim%transfx%get_transformation(ipatch)
       transf => sim%transfx%transfs(ipatch+1)%t
       !       call sll_display(transf%mesh)
       
       num_pts1 = sim%transfx%get_num_cells_eta1(ipatch) + 1! logical_m%num_cells1+1
       num_pts2 = sim%transfx%get_num_cells_eta2(ipatch) + 1!logical_m%num_cells2+1
       delta1   = sim%transfx%get_delta_eta1(ipatch)!logical_m%delta_eta1
       delta2   = sim%transfx%get_delta_eta2(ipatch)!logical_m%delta_eta2
       eta1_min = sim%transfx%get_eta1_min(ipatch)!logical_m%eta1_min
       eta2_min = sim%transfx%get_eta2_min(ipatch)!logical_m%eta2_min

       print *, "num_patches = ", num_patches
       print *, "ipatch = ", ipatch
       print *, "num_pts1 = ", num_pts1
       print *, "num_pts2 = ", num_pts2
       print *, "delta1 = ", delta1
       print *, "delta2 = ", delta2
       print *, "eta1_min = ", eta1_min
       print *, "eta2_min = ", eta2_min

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
             val_rho = 0.0_f64
             call a11_field_mat%set_value_at_indices ( i, j, ipatch, val_a11 ) 
             call a12_field_mat%set_value_at_indices ( i, j, ipatch, val_a12 ) 
             call a21_field_mat%set_value_at_indices ( i, j, ipatch, val_a21 ) 
             call a22_field_mat%set_value_at_indices ( i, j, ipatch, val_a22 ) 
             call b1_field_vect%set_value_at_indices ( i, j, ipatch, val_b1 ) 
             call b2_field_vect%set_value_at_indices ( i, j, ipatch, val_b2 ) 
             call c_field%set_value_at_indices  ( i, j, ipatch, val_c ) 
             call phi%set_value_at_indices( i, j, ipatch, val_phi)
             call rho%set_value_at_indices( i, j, ipatch, val_rho)
          end do
       end do
    end do

    time = sll_time_elapsed_since(t0)
    print*, 'rank: ', sim%my_rank, 'time to initialize MP fields =', time
       
    print *, 'rank: ', sim%my_rank, 'updating interpolation coefficients.'

    call sll_set_time_mark(t0)

    call a11_field_mat%update_interpolation_coefficients()
    call a12_field_mat%update_interpolation_coefficients()
    call a21_field_mat%update_interpolation_coefficients()
    call a22_field_mat%update_interpolation_coefficients()
    call b1_field_vect%update_interpolation_coefficients()
    call b2_field_vect%update_interpolation_coefficients()
    call c_field%update_interpolation_coefficients()
    call phi%update_interpolation_coefficients()
    call rho%update_interpolation_coefficients()

    time = sll_time_elapsed_since(t0)
    print*, 'rank: ', sim%my_rank, 'time to update coefficients =', time

    print *, '********** INITIALIZED MULTIPATCH FIELDS & INTERPOLANTS *******'
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


    ! Creating and intializing distribution function.

    print *, 'creating and initializing distribution function MP.'
    call sll_set_time_mark(t0)
    f_mp => sll_new_distribution_function_4d_multipatch( sll_world_collective, &
         sim%transfx, sim%mesh2d_v, sim%nproc_x1, sim%nproc_x2 )
    time = sll_time_elapsed_since(t0)
    print*, 'rank: ', sim%my_rank, 'time to create MP F =', time

    call sll_set_time_mark(t0)
    call f_mp%initialize( sim%init_func, sim%params ) 
    time = sll_time_elapsed_since(t0)
    print*, 'rank: ', sim%my_rank, 'time to initialize MP F =', time

    print *, 'reconfiguring to sequential_x1x2 mode...'

    ! First dt/2 advection for eta1-eta2:
    call sll_set_time_mark(t0)
    call f_mp%set_to_sequential_x1x2()
    time = sll_time_elapsed_since(t0)
    print*, 'rank: ', sim%my_rank, 'time to reconfigure F =', time

    print *, 'First-time advection in x1x2.'
    call sll_set_time_mark(t0)
    call advection_x1x2( sim, layer_x1x2, f_mp, 0.5*sim%dt)    
    time = sll_time_elapsed_since(t0)
    print*, 'rank: ', sim%my_rank, 'time for advection x1x2 =', time

    print *, '********** COMPLETED X1X2 ADVECTION ***********'

    ! Initialize the poisson plan before going into the main loop.
    print *, 'Creating and initializing elliptic solver.'
    call sll_set_time_mark(t0)
    sim%qns => new_general_elliptic_solver_mp( &
         sim%quadrature_type1,& !ES_GAUSS_LEGENDRE, &  ! put in arguments
         sim%quadrature_type2,& !ES_GAUSS_LEGENDRE, &  ! put in arguments
         sim%transfx)
    time = sll_time_elapsed_since(t0)
    print*, 'rank: ', sim%my_rank, 'time for creating field solver', time
    print*, 'about to factorize matrix qns'
    

    print *, 'factorizing solver matrices.'
    call sll_set_time_mark(t0)
    call factorize_mat_es_mp(&
         sim%qns, & 
         a11_field_mat, &
         a12_field_mat, &
         a21_field_mat, &
         a22_field_mat, &
         b1_field_vect, &
         b2_field_vect, &
         c_field)
    time = sll_time_elapsed_since(t0)
    print*, 'rank: ', sim%my_rank, 'time for factorizing matrices', time

    print*, '--- ended factorization matrix qns'

    call sim%interp_x3%initialize( &
         sim%mesh2d_v%num_cells1+1, &
         sim%mesh2d_v%eta1_min, &
         sim%mesh2d_v%eta1_max, &
         SLL_HERMITE)!sim%bc_vx_0)
    
    call sim%interp_x4%initialize( &
         sim%mesh2d_v%num_cells2+1, &
         sim%mesh2d_v%eta2_min, &
         sim%mesh2d_v%eta2_max, &
         SLL_HERMITE)! sim%bc_vy_0)
    
    
    ! ------------------------------------------------------------------------
    !
    !                                MAIN LOOP
    !
    ! ------------------------------------------------------------------------
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

    print*, ' ... finished initialization, entering main loop.'    
    
    do itime=1,sim%num_iterations
       if(sim%my_rank == 0) then
          print *, 'Starting iteration ', itime, ' of ', sim%num_iterations
       end if

       call f_mp%set_to_sequential_x3x4()       
       
       call compute_charge_density_multipatch(f_mp, rho)

       call rho%update_interpolation_coefficients( )
       
       if(sim%my_rank == 0) then
          call rho%write_to_file(0)
       end if
       
       call sll_set_time_mark(t0)         
       call solve_general_coordinates_elliptic_eq_mp(&
            sim%qns,&
            rho,&
            phi)
       time = sll_time_elapsed_since(t0)
       print*, 'rank: ', sim%my_rank, 'time to solve QNS =', time
       
       print*, '--- end solve qns'
       
       if(sim%my_rank == 0) then
          call phi%write_to_file(0)
       end if
       
       ! advection x3
       
       ! put these comments elsewhere    
       !--> Compute energy kinetic, potential and total
       !print*, 'compute nrj'
       !call compute_energy_qns(sim,phi)
       
       !--> Compute L1 norm, L2 norm, L infini norm
       !print*, 'compute  L^p'
       !call compute_norm_L1_L2_Linf_qns(sim)
    
       !call writeHDF5_diag_qns( sim )

     
       if(sim%my_rank == 0) then
          call phi%write_to_file(itime)
       end if
       
       efield_energy_total = 0.0_f64

       print*, 'advection vx'
       call advection_x3(sim, f_mp, phi, sim%dt, efield_energy_total)

       !print*, 'energy total', efield_energy_total
       !efield_energy_total = sqrt(efield_energy_total)
       
       numpart = 0.0_f64
       ! dt in vy...(x4)
       print *, 'advection vy (x4)'
       call advection_x4(sim, f_mp, phi, sim%dt, numpart )
       print *, 'finished advection in vy (x4)'
       call f_mp%set_to_sequential_x1x2()
       print *, 'reconfigured parallelization of f'

       ! Approximate the integral of the distribution function along all
       ! directions.
       num_particles_local(buffer_counter) = f_mp%compute_moment(0)
       buffer_energy(buffer_counter)       = efield_energy_total

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

       print *, 'advection in xy'
       call advection_x1x2(sim, layer_x1x2, f_mp, sim%dt) 

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
       
       
      print *, 'writing some diagnostics...'
       !--> Save results in HDF5 files
       if ( mod(itime,BUFFER_SIZE) == 0) then
          !--> Compute energy kinetic, potential and total
          !print*, 'compute nrj'
!          call compute_energy_qns(sim,phi)
          
          !--> Compute L1 norm, L2 norm, L infini norm
          !print*, 'compute  L^p'
!          call compute_norm_L1_L2_Linf_qns(sim)
          
          call writeHDF5_diag_qns( sim )
          
       end if
       print *, 'finished iteration.'
    end do ! main loop
#undef BUFFER_SIZE

  end subroutine run_4d_qns_general_mp
  

  subroutine advection_x1x2( &
       sim, &
       field_x1x2, &
       f_mp, &
       deltat)
    
    class(sll_simulation_4d_qns_general_multipatch)        :: sim
    type(sll_scalar_field_multipatch_2d), pointer          :: field_x1x2
    type(sll_distribution_function_4d_multipatch), pointer :: f_mp
    sll_real64, intent(in)                                 :: deltat
    sll_real64, dimension(1:2,1:2) :: inv_j,jac_m
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: loc_sz_x3
    sll_int32  :: loc_sz_x4
    sll_real64, dimension(1:4) :: eta
    sll_real64 :: alpha1, alpha2
    sll_real64 :: eta1, eta2 !, eta3, eta4
    sll_real64 :: eta1_min    
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min    
    sll_real64 :: eta2_max
    sll_int32  :: i 
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_int32  :: ip
    sll_real64 :: f_interpolated
    sll_int32  :: num_patches
    type(sll_time_mark) :: t0 !delete this
    sll_real64 :: time     ! delete this
    ! Here we do something ugly, which is to call the 'get_local_sizes()'
    ! method just to obtain the local sizes in the x3 and x4 directions, which
    ! are 'cut' (parallelized), as we are carrying out a sequential x1x2
    ! advection. So we are using the information that since the multipatch
    ! approach is applied only to the space coordinates (x1, x2), the limits
    ! of the indices in x3 and x4 are shared by all the patches. 
    ! In other words, this call only helps to hide a call like
    ! compute_local_sizes() using the layout of the distribution function of
    ! any of the patches. Sorry about this... :-(
    call f_mp%get_local_data_sizes( &
         0, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )

    num_patches = f_mp%num_patches
       
    do l=1,loc_sz_x4
       do k=1,loc_sz_x3
          ! Make the 'recyclable' field point to the right data on the
          ! distribution function.
          do ip=0, num_patches-1
             call field_x1x2%set_patch_data_pointer(&
                  ip, &
                  f_mp%get_x1x2_data_slice_pointer(ip,k,l) )
          end do
          ! update_interpolation_coefficients() already implies a loop around
          ! the patches.
          call field_x1x2%update_interpolation_coefficients()

          do ip=0, num_patches-1

             call f_mp%get_local_data_sizes( &
                  ip, &
                  loc_sz_x1, &
                  loc_sz_x2, &
                  loc_sz_x3, &
                  loc_sz_x4 )

             eta1_min = f_mp%transf%get_eta1_min(ip)
             eta1_max = f_mp%transf%get_eta1_max(ip)
             eta2_min = f_mp%transf%get_eta2_min(ip)
             eta2_max = f_mp%transf%get_eta2_max(ip)

             do j=1,loc_sz_x2
                do i=1,loc_sz_x1
                   eta(:) = f_mp%get_eta_coordinates( ip, (/i,j,k,l/) )
                   call sll_set_time_mark(t0) ! delete this
                   inv_j(:,:)  = &
                     field_x1x2%transf%inverse_jacobian_matrix(eta(1),eta(2),ip)
                   time = sll_time_elapsed_since(t0) ! delete
                   print *, 'time for inverse jacobian: ', time
                   alpha1 = -deltat*(inv_j(1,1)*eta(3) + inv_j(1,2)*eta(4))
                   alpha2 = -deltat*(inv_j(2,1)*eta(3) + inv_j(2,2)*eta(4))
                   eta1   = eta1+alpha1
                   eta2   = eta2+alpha2
                   ! Apply the BCs between the patches. This involves:
                   ! - Checking if the eta coordinates are outside of the patch.
                   ! - Checking if the patch is connected in such direction.
                   ! - Treat accordingly. All this needs implementation.
                   ! - If the coordinates are inside the patch, carry out the
                   !   interpolation to update the distribution function...
                   
                   ! FIXME: NO TREATMENT OF INTERNAL OR EXTERNAL BCs!!!

                   ! ***************************
                   !
                   !    WORK NEEDED HERE!!!!!
                   !
                   ! ***************************
                   if( eta1 >= eta1_min .and. eta1 <= eta1_max .and. &
                       eta2 >= eta2_min .and. eta2 <= eta2_max ) then
                      f_interpolated = field_x1x2%value_at_point(eta1, eta2, ip)
                   ! Eatch x1x2 slice is entirely contained in the process, so
                   ! no local to global or anything else should be needed.
                      call field_x1x2%set_value_at_indices(i,j,ip,f_interpolated)
                   end if
                   
                end do
             end do
          end do
       end do
    end do
  end subroutine advection_x1x2

  subroutine advection_x3(sim, f_mp, phi, deltat, efield_energy_total)
    type(sll_distribution_function_4d_multipatch), pointer :: f_mp
    class(sll_simulation_4d_qns_general_multipatch) :: sim
    type(sll_scalar_field_multipatch_2d), pointer   :: phi
    sll_real64, intent(in) :: deltat
    sll_real64, intent(out) :: efield_energy_total
    type(sll_coordinate_transformation_multipatch_2d), pointer :: transf
    sll_int32 :: ip
    sll_int32 :: num_patches
    sll_real64, dimension(1:4) :: eta
    sll_real64 :: x
    sll_real64 :: y
    sll_real64, dimension(1:2,1:2) :: jac_mat
    sll_real64, dimension(1:2,1:2) :: inv_j
    sll_real64, dimension(:), pointer :: line
    sll_int32  :: nc_x3
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: alpha3
    sll_int32  :: loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 
    sll_int32  :: i, j, k, l
    sll_real64 :: ex
    sll_real64 :: ey

    nc_x3               = sim%mesh2d_v%num_cells1    
    efield_energy_total = 0.0_f64
    num_patches         = f_mp%num_patches
    transf              => f_mp%transf

    print *, 'inside advection x3'
    ! Start with dt in vx...(x3)
    do ip=0,num_patches-1

       call f_mp%get_local_data_sizes( &
            ip, &
            loc_sz_x1, &
            loc_sz_x2, &
            loc_sz_x3, &
            loc_sz_x4 )

       delta1 = transf%get_delta_eta1(ip)
       delta2 = transf%get_delta_eta2(ip)

       do l=1,loc_sz_x4 !sim%mesh2d_v%num_cells2+1
          do j=1,loc_sz_x2
             do i=1,loc_sz_x1

                eta(:)  = f_mp%get_eta_coordinates( ip, (/i,j,1,l/) )
                x       = transf%x1( eta(1), eta(2), ip )
                y       = transf%x2( eta(1), eta(2), ip )
                inv_j   = transf%inverse_jacobian_matrix(eta(1),eta(2),ip)
                jac_mat = transf%jacobian_matrix(eta(1),eta(2),ip)
                ex = - phi%first_deriv_eta1_value_at_point(eta(1),eta(2),ip)
                ey = - phi%first_deriv_eta2_value_at_point(eta(1),eta(2),ip)

                !sim%values_ex (global_indices(1),global_indices(2) ) = ex
                !sim%values_ey (global_indices(1),global_indices(2) ) = ey
                alpha3 = -deltat*(inv_j(1,1)*ex + inv_j(2,1)*ey)
                ! Add correction in case of external applied electric field:
                !alpha3 = alpha3 -sim%dt*(elec_field_ext_1%value_at_point(x,y) )

                ! ATTENTION, the following 'weird' construct needs to be verified
                ! in any case, we should be using some routine from the 
                ! advection module and not the deprecated ones from the
                ! interpolators!
                line => f_mp%get_x3_line_pointer(ip, i,j,l)
                line(:) = sim%interp_x3%interpolate_array_disp( &
                  nc_x3 + 1, &
                  line(:), &
                  alpha3 )
                  ! Extra work to calculate the electric field energy. Where
                  ! does this formula come from?
                  efield_energy_total = efield_energy_total + delta1*delta2* &
                       abs(jac_mat(1,1)*jac_mat(2,2) - &
                           jac_mat(1,2)*jac_mat(2,1))*&
                       abs( (inv_j(1,1)**2 + inv_j(1,2)**2)*ex**2 + &
                       2.0_f64*(inv_j(1,1)*inv_j(2,1)+inv_j(1,2)*inv_j(2,2))*&
                       abs(ex)*abs(ey) + &
                       (inv_j(2,1)**2 + inv_j(2,2)**2)*ey**2)
               end do
            end do
         end do
      end do
      efield_energy_total = sqrt(efield_energy_total)
    end subroutine advection_x3
  
  
  subroutine advection_x4( sim, f_mp, phi, deltat, integral_f )
    type(sll_distribution_function_4d_multipatch), pointer :: f_mp
    class(sll_simulation_4d_qns_general_multipatch) :: sim
    type(sll_scalar_field_multipatch_2d), pointer   :: phi
    sll_real64, intent(in) :: deltat
    sll_real64, intent(inout) :: integral_f
    type(sll_coordinate_transformation_multipatch_2d), pointer :: transf
    sll_int32 :: ip
    sll_int32 :: num_patches
    sll_real64, dimension(1:4) :: eta
    sll_real64, dimension(1:2,1:2) :: inv_j
    sll_real64, dimension(1:2,1:2) :: jac_mat
    sll_int32 :: loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 
    sll_real64 :: alpha4
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_int32  :: i, j, k, l
    sll_real64 :: ex
    sll_real64 :: ey
    sll_real64 :: x
    sll_real64 :: y
    sll_int32  :: nc_x4
    sll_real64, dimension(:), pointer :: line
    sll_real64, dimension(:,:,:,:), pointer :: f4d

    nc_x4               = sim%mesh2d_v%num_cells2
    num_patches         = f_mp%num_patches
    transf              => f_mp%transf
    ! dt in vy...(x4)
    do ip=0,num_patches-1
       
       call f_mp%get_local_data_sizes( &
            ip, &
            loc_sz_x1, &
            loc_sz_x2, &
            loc_sz_x3, &
            loc_sz_x4 )

       delta1 = transf%get_delta_eta1(ip)
       delta2 = transf%get_delta_eta2(ip)
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1

             eta(:)  = f_mp%get_eta_coordinates( ip, (/i,j,1,1/) )
             do k=1, loc_sz_x3 !sim%mesh2d_v%num_cells1+1

                x       = transf%x1( eta(1), eta(2), ip )
                y       = transf%x2( eta(1), eta(2), ip )
                inv_j   = transf%inverse_jacobian_matrix(eta(1),eta(2),ip)
                jac_mat = transf%jacobian_matrix(eta(1),eta(2),ip)
                ex = - phi%first_deriv_eta1_value_at_point(eta(1),eta(2),ip)
                ey = - phi%first_deriv_eta2_value_at_point(eta(1),eta(2),ip)
                alpha4 = -deltat*(inv_j(1,2)*ex + inv_j(2,2)*ey)
                ! Add correction in case that external field is used
                ! alpha4 = alpha4 -sim%dt*(elec_field_ext_2%value_at_point(x,y))
                line => f_mp%get_x4_line_pointer(ip,i,j,k)
                line(:) = sim%interp_x4%interpolate_array_disp( &
                  nc_x4 + 1, &
                  line(:), &
                  alpha4 )
             end do
             ! Check this calculation, why is this not multiplied by the deltas
             ! in velocity space?? Fix this, because it it wrong even...
!             f4d => f_mp%get_full_patch_data_pointer(ip)
 !            integral_f = integral_f + &
  !                sum(f4d(i,j,:,:))*abs(jac_mat(1,1)*jac_mat(2,2) - &
   !                                     jac_mat(1,2)*jac_mat(2,1))
          end do
       end do
    end do
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
    
!!$    ix1_diag   = int(sim%nc_x1/2)
!!$    ix2_diag   = int(sim%nc_x2/3)
!!$    iv1_diag   = int(sim%nc_x3/4)
!!$    iv2_diag   = int(sim%nc_x4/3)

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
  
! Save this for now as a reference for the potential energy.          
!!$          nrj_pot = nrj_pot + &
!!$               ((inv_j(1,1)*ex + inv_j(2,1)*ey)**2+ &
!!$               (inv_j(1,2)*ex + inv_j(2,2)*ey)**2) * 0.5 * val_jac * &
!!$               delta_eta1*delta_eta2
  
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
#if 0  
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
#endif
 
end module sll_simulation_4d_qns_general_multipatch_module
