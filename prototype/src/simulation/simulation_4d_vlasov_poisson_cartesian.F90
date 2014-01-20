module sll_simulation_4d_vlasov_poisson_cartesian
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
  use sll_collective
  use sll_remapper
  use sll_constants
  use sll_cubic_spline_interpolator_1d
  use sll_test_4d_initializer
  use sll_poisson_2d_periodic_cartesian_par
  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  implicit none

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_4d_vlasov_poisson_cart
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
     ! for initializers
     type(init_test_4d_par)                     :: init_4d
     type(simple_cartesian_4d_mesh), pointer    :: mesh4d
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
     
     ! for remap
     type(layout_4D), pointer :: sequential_x1x2
     type(layout_4D), pointer :: sequential_x3x4
     type(layout_2D), pointer :: rho_seq_x1
     type(layout_2D), pointer :: rho_seq_x2
     type(layout_2D), pointer :: split_rho_layout ! not sequential in any dir.
     type(remap_plan_2D_real64), pointer :: split_to_seqx1
     type(remap_plan_2D_real64), pointer :: seqx1_to_seqx2
     ! remaps for the electric field data
!     type(remap_plan_2D), pointer :: efld_split_to_seqx1
     type(remap_plan_2D_comp64), pointer :: efld_seqx1_to_seqx2
     type(remap_plan_2D_comp64), pointer :: efld_seqx2_to_split
     type(remap_plan_4D_real64), pointer :: seqx1x2_to_seqx3x4
     type(remap_plan_4D_real64), pointer :: seqx3x4_to_seqx1x2
     ! interpolators and their pointers
     type(cubic_spline_1d_interpolator) :: interp_x1
     type(cubic_spline_1d_interpolator) :: interp_x2
     type(cubic_spline_1d_interpolator) :: interp_x3
     type(cubic_spline_1d_interpolator) :: interp_x4
     ! Field accumulator
     sll_comp64, dimension(:,:), allocatable :: efield_x1
     sll_comp64, dimension(:,:), allocatable :: efield_x2
     sll_comp64, dimension(:,:), allocatable :: efield_split
   contains
     procedure, pass(sim) :: run => run_vp4d_cartesian
     procedure, pass(sim) :: init_from_file => init_vp4d_par_cart
  end type sll_simulation_4d_vlasov_poisson_cart

  interface sll_delete
     module procedure delete_vp4d_par_cart
  end interface sll_delete

contains

  subroutine init_vp4d_par_cart( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_vlasov_poisson_cart), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
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
  end subroutine init_vp4d_par_cart


  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_vp4d_cartesian(sim)
    class(sll_simulation_4d_vlasov_poisson_cart), intent(inout) :: sim
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: loc_sz_x3
    sll_int32  :: loc_sz_x4
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_real64 :: vmin
    sll_real64 :: delta
    sll_real64 :: alpha
    sll_int32  :: itemp
    sll_int32  :: ierr
    sll_int32  :: itime
    sll_real64 :: ex
    sll_real64 :: ey

    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)

    ! allocate the layouts...
    sim%sequential_x1x2  => new_layout_4D( sll_world_collective )
    sim%sequential_x3x4  => new_layout_4D( sll_world_collective )
    sim%rho_seq_x1       => new_layout_2D( sll_world_collective )
    sim%rho_seq_x2       => new_layout_2D( sll_world_collective )
    sim%split_rho_layout => new_layout_2D( sll_world_collective )

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
    
    call initialize_layout_with_distributed_4D_array( &
         sim%nc_x1+1, &
         sim%nc_x2+1, &
         sim%nc_x3+1, &
         sim%nc_x4+1, &
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
         sim%nc_x1+1, &
         sim%nc_x2+1, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%split_rho_layout )
    
    call initialize_layout_with_distributed_2D_array( &
         sim%nc_x1+1, &
         sim%nc_x2+1, &
         1, &
         sim%world_size, &
         sim%rho_seq_x1 )
    
    call compute_local_sizes_2d( sim%rho_seq_x1, loc_sz_x1, loc_sz_x2 )
    SLL_ALLOCATE(sim%rho_x1(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%phi_x1(loc_sz_x1,loc_sz_x2),ierr)
    ! Experiment with a dedicated array to store the values of the electric
    ! field in each point of the grid.
    SLL_ALLOCATE(sim%efield_x1(loc_sz_x1,loc_sz_x2),ierr)

    call initialize_layout_with_distributed_2D_array( &
         sim%nc_x1+1, &
         sim%nc_x2+1, &
         sim%world_size, &
         1, &
         sim%rho_seq_x2 )
    call compute_local_sizes_2d( sim%rho_seq_x2, loc_sz_x1, loc_sz_x2 )
    SLL_ALLOCATE(sim%rho_x2(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%phi_x2(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_x2(loc_sz_x1,loc_sz_x2),ierr)
    
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
         sim%nc_x1+1, &
         sim%nc_x2+1, &
         sim%nc_x3+1, &
         sim%nc_x4+1, &
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
    
    ! This layout is also useful to represent the charge density array. Since
    ! this is a result of a local reduction on x3 and x4, the new layout is
    ! 2D but with the same dimensions of the process mesh in x1 and x2.
    SLL_ALLOCATE(sim%rho_split(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%efield_split(loc_sz_x1,loc_sz_x2),ierr)
    
    call compute_local_sizes_4d( sim%sequential_x3x4, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3, &
         loc_sz_x4 )

    SLL_ALLOCATE(sim%f_x3x4(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)
    
    ! These dimensions are also the ones needed for the array where we store
    ! the intermediate results of the charge density computation.
    SLL_ALLOCATE(sim%partial_reduction(loc_sz_x1,loc_sz_x2, loc_sz_x3),ierr)
    
    ! Initialize the initial distribution function data. We do this with an
    ! initializer object which needs to be initialized itself! Note also that 
    ! the mesh is described in global terms. This should be fine as meshes are
    ! supposed to be read-only entities.
    ! This should be done elsewhere...
    sim%mesh4d => new_cartesian_4d_mesh( sim%nc_x1, sim%nc_x2, sim%nc_x3, &
         sim%nc_x4, &
         0.0_f64, 1.0_f64, &
         0.0_f64, 1.0_f64, &
         -1.0_f64, 1.0_f64, &
         -1.0_f64, 1.0_f64 )

    call load_test_4d_initializer( sim%init_4d, &
         NODE_CENTERED_FIELD, &
         sim%mesh4d, &
         0.1_f64, &
         sim%sequential_x3x4 )
    call sim%init_4d%f_of_4args(sim%f_x3x4)
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
    sim%rho_split(:,:) = 0.0

    call compute_charge_density( &
         sim%mesh4d, &
         size(sim%f_x3x4,1), &
         size(sim%f_x3x4,2), &
         sim%f_x3x4, &
         sim%partial_reduction, &
         sim%rho_split )

    ! Re-arrange rho_split in a way that permits sequential operations in x1, to
    ! feed to the Poisson solver.
    sim%split_to_seqx1 => &
         NEW_REMAP_PLAN(sim%split_rho_layout, sim%rho_seq_x1, sim%rho_split)
    call apply_remap_2D( sim%split_to_seqx1, sim%rho_split, sim%rho_x1 )

    ! We are in a position now to compute the electric potential.
    ! Initialize the poisson plan
    sim%poisson_plan => new_poisson_2d_periodic_plan_cartesian_par( &
         sim%rho_seq_x1, &
         sim%nc_x1, &
         sim%nc_x2, &
         1.0_f64, &    ! parametrize with mesh values
         1.0_f64 )     ! parametrize with mesh values

    ! solve for the electric potential
    call solve_poisson_2d_periodic_cartesian_par( &
         sim%poisson_plan, &
         sim%rho_x1, &
         sim%phi_x1)

    ! compute the values of the electric field. rho is configured for 
    ! sequential operations in x1, thus we start by computing the E_x component.
    ! The following call is inefficient and unnecessary. The local sizes for
    ! the arrays should be kept around as parameters basically and not on 
    ! variables whose content could be anything... This will have to do for now.
    call compute_local_sizes_2d( sim%rho_seq_x1, loc_sz_x1, loc_sz_x2 )
    call compute_electric_field_x1( &
         sim%phi_x1, &
         loc_sz_x1, &
         loc_sz_x2, &
         sim%mesh4d%delta_x1, &
         sim%efield_x1 )

    ! note that we are 'recycling' the layouts used for the other arrays because
    ! they represent an identical configuration.
    sim%efld_seqx1_to_seqx2 => &
         NEW_REMAP_PLAN( sim%rho_seq_x1, sim%rho_seq_x2, sim%efield_x1)
    sim%seqx1_to_seqx2 => &
         NEW_REMAP_PLAN( sim%rho_seq_x1, sim%rho_seq_x2, sim%phi_x1 )
    call apply_remap_2D( sim%efld_seqx1_to_seqx2, sim%efield_x1, sim%efield_x2 )
    call apply_remap_2D( sim%seqx1_to_seqx2, sim%phi_x1, sim%phi_x2 )
    call compute_local_sizes_2d( sim%rho_seq_x2, loc_sz_x1, loc_sz_x2 )
    call compute_electric_field_x2( &
         sim%phi_x2, &
         loc_sz_x1, &
         loc_sz_x2, &
         sim%mesh4d%delta_x2, &
         sim%efield_x2 )

    ! But now, to make the electric field data configuration compatible with
    ! the sequential operations in x2x3 we need still another remap operation.
    sim%efld_seqx2_to_split => &
         NEW_REMAP_PLAN( sim%rho_seq_x2, sim%split_rho_layout,sim%efield_x2 )
    call apply_remap_2D( &
         sim%efld_seqx2_to_split, &
         sim%efield_x2, &
         sim%efield_split )
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
         sim%nc_x1+1, &
         sim%mesh4d%x1_min, &
         sim%mesh4d%x1_max, &
         SLL_PERIODIC)

    call sim%interp_x2%initialize( &
         sim%nc_x2+1, &
         sim%mesh4d%x2_min, &
         sim%mesh4d%x2_max, &
         SLL_PERIODIC)

    call sim%interp_x3%initialize( &
         sim%nc_x3+1, &
         sim%mesh4d%x3_min, &
         sim%mesh4d%x3_max, &
         SLL_HERMITE)

    call sim%interp_x4%initialize( &
         sim%nc_x4+1, &
         sim%mesh4d%x4_min, &
         sim%mesh4d%x4_max, &
         SLL_HERMITE)

    call compute_local_sizes_2d( sim%rho_seq_x1, loc_sz_x1, loc_sz_x2 )


    ! ------------------------------------------------------------------------
    !
    !                                MAIN LOOP
    !
    ! ------------------------------------------------------------------------

       sim%seqx3x4_to_seqx1x2 => &
           NEW_REMAP_PLAN(sim%sequential_x3x4,sim%sequential_x1x2,sim%f_x3x4)

       sim%seqx1x2_to_seqx3x4 => &
           NEW_REMAP_PLAN(sim%sequential_x1x2,sim%sequential_x3x4,sim%f_x1x2)


    do itime=1, sim%num_iterations
       ! The splitting scheme used here is meant to attain a dt^2 accuracy.
       ! We use:
       !
       ! dt/2 in vx
       ! dt/2 in vy
       ! dt/4 in vx
       ! dt/2 in y
       ! dt   in x
       ! dt/2 in y
       ! dt/4 in vx
       ! dt/2 in vy
       !
       ! Note that two dt/4 in vx have been joined in the first step.

       ! Note: Since the Ex and Ey values are used separately, the proposed
       ! data structure is actually not good. These field values should be kept
       ! separate.
       call compute_local_sizes_4d( sim%sequential_x3x4, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 

       ! Start with dt/2 in vx...(x3)
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             do l=1,sim%mesh4d%num_cells4+1
                ex    = real(sim%efield_split(i,j),f64)
                alpha = -ex*0.5_f64*sim%dt
                ! interpolate_array_disp() has an interface that must be changed
                sim%f_x3x4(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                     sim%nc_x3+1, &
                     sim%f_x3x4(i,j,:,l), &
                     alpha )
             end do
          end do
       end do

       ! dt/2 in vy...(x4)
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             do k=1,sim%mesh4d%num_cells3+1
                ey    = aimag(sim%efield_split(i,j))
                alpha = -ey*0.5_f64*sim%dt
                ! interpolate_array_disp() has an interface that must be changed
                sim%f_x3x4(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
                     sim%nc_x4+1, &
                     sim%f_x3x4(i,j,k,:), &
                     alpha )
             end do
          end do
       end do

       ! dt/4 in vx:
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             do l=1,sim%mesh4d%num_cells4+1
                ex    = real(sim%efield_split(i,j),f64)
                alpha = -ex*0.25_f64*sim%dt
                ! interpolate_array_disp() has an interface that must be changed
                sim%f_x3x4(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                     sim%nc_x3+1, &
                     sim%f_x3x4(i,j,:,l), &
                     alpha )
             end do
          end do
       end do

       ! Proceed to the advections in the spatial directions, 'x' and 'y'
       ! Reconfigure data. 
       call apply_remap_4D( sim%seqx3x4_to_seqx1x2, sim%f_x3x4, sim%f_x1x2 )
       
       ! what are the new local limits on x3 and x4? It is bothersome to have
       ! to make these calls...
       call compute_local_sizes_4d( sim%sequential_x1x2, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 )
            
       ! dt/2 in 'y' (x2)
       do l=1,loc_sz_x4
          do k=1,loc_sz_x3
             do i=1,sim%mesh4d%num_cells1+1
                vmin = sim%mesh4d%x4_min
                delta = sim%mesh4d%delta_x4
                alpha = (vmin + (l-1)*delta)*sim%dt*0.5_f64
                !call sim%interp_x1%compute_interpolants( sim%f_x1x2(i,:,k,l) )
                ! interpolate_array_disp() has an interface that must be changed
                sim%f_x1x2(i,:,k,l) = sim%interp_x2%interpolate_array_disp( &
                     sim%nc_x2+1, &
                     sim%f_x1x2(i,:,k,l), &
                     alpha )
             end do
          end do
       end do

       ! dt in 'x' (x1)
       do l=1,loc_sz_x4
          do k=1,loc_sz_x3
             do j=1,sim%mesh4d%num_cells2+1
                vmin = sim%mesh4d%x3_min
                delta = sim%mesh4d%delta_x3
                alpha = (vmin + (k-1)*delta)*sim%dt
                !call sim%interp_x1%compute_interpolants( sim%f_x1x2(:,j,k,l) )
                ! interpolate_array_disp() has an interface that must be changed
                sim%f_x1x2(:,j,k,l) = sim%interp_x1%interpolate_array_disp( &
                     sim%nc_x1+1, &
                     sim%f_x1x2(:,j,k,l), &
                     alpha )
             end do
          end do
       end do
       
       ! dt/2 in 'y' (x2)
       do l=1,loc_sz_x4
          do k=1,loc_sz_x3
             do i=1,sim%mesh4d%num_cells1+1
                vmin = sim%mesh4d%x4_min
                delta = sim%mesh4d%delta_x4
                alpha = (vmin + (l-1)*delta)*sim%dt*0.5_f64
                !call sim%interp_x1%compute_interpolants( sim%f_x1x2(i,:,k,l) )
                ! interpolate_array_disp() has an interface that must be changed
                sim%f_x1x2(i,:,k,l) = sim%interp_x2%interpolate_array_disp( &
                     sim%nc_x2+1, &
                     sim%f_x1x2(i,:,k,l), &
                     alpha )
             end do
          end do
       end do

       ! And reconfigure the data for the last set of advections in vx and vy.
       ! This will also help for the calculation of the charge density later
       ! on.
       call apply_remap_4D( sim%seqx1x2_to_seqx3x4, sim%f_x1x2, sim%f_x3x4 )

       call compute_local_sizes_4d( sim%sequential_x3x4, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 ) 

       ! dt/4 in vx
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             do l=1,sim%mesh4d%num_cells4
                ex    = real(sim%efield_split(i,j),f64)
                alpha = -ex*0.25_f64*sim%dt
                ! interpolate_array_disp() has an interface that must be changed
                sim%f_x3x4(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                     sim%nc_x3, &
                     sim%f_x3x4(i,j,:,l), &
                     alpha )
             end do
          end do
       end do

       ! dt/2 in vy...(x4)
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             do k=1,sim%mesh4d%num_cells3
                ey    = aimag(sim%efield_split(i,j))
                alpha = -ey*0.5_f64*sim%dt
                ! interpolate_array_disp() has an interface that must be changed
                sim%f_x3x4(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
                     sim%nc_x4, &
                     sim%f_x3x4(i,j,k,:), &
                     alpha )
             end do
          end do
       end do

       ! Compute the fields:
       ! 1. Compute charge density.
       ! 2. Reconfigure charge density to feed to Poisson solver
       
       call compute_charge_density( &
            sim%mesh4d, &
            size(sim%f_x3x4,1), &
            size(sim%f_x3x4,2), &
            sim%f_x3x4, &
            sim%partial_reduction, &
            sim%rho_split )
       
       ! 2d charge density is 'fully split', no sequential operations can be
       ! fully done.
       call apply_remap_2D( sim%split_to_seqx1, sim%rho_split, sim%rho_x1 )
       
       ! Compute the electric potential.
       call solve_poisson_2d_periodic_cartesian_par( &
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
       call compute_local_sizes_2d( sim%rho_seq_x1, loc_sz_x1, loc_sz_x2 )
       call compute_electric_field_x1( &
            sim%phi_x1, &
            loc_sz_x1, &
            loc_sz_x2, &
            sim%mesh4d%delta_x1, &
            sim%efield_x1 )

       ! Prepare computation of the electric field in x2. All these steps should
       ! disappear, as the Poisson solver can compute this directly.
       call apply_remap_2D(sim%efld_seqx1_to_seqx2,sim%efield_x1,sim%efield_x2)
       call apply_remap_2D( sim%seqx1_to_seqx2, sim%phi_x1, sim%phi_x2 )
       call compute_local_sizes_2d( sim%rho_seq_x2, loc_sz_x1, loc_sz_x2 )
       call compute_electric_field_x2( &
            sim%phi_x2, &
            loc_sz_x1, &
            loc_sz_x2, &
            sim%mesh4d%delta_x2, &
            sim%efield_x2 )
       ! But now, to make the electric field data configuration compatible with
       ! the sequential operations in x2x3 we need still another remap 
       ! operation.
       call apply_remap_2D( &
            sim%efld_seqx2_to_split, &
            sim%efield_x2, &
            sim%efield_split)

       ! ...and another half time step advection in the velocities.
       ! Start with vx...(x3)
       ! Note: Since the Ex and Ey values are used separately, the proposed
       ! data structure is actually not good. These field values should be kept
       ! separate.
       call compute_local_sizes_4d( sim%sequential_x3x4, &
            loc_sz_x1, loc_sz_x2, loc_sz_x3, loc_sz_x4 )
       
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             do l=1,sim%mesh4d%num_cells4
                ex    = real(sim%efield_split(i,j),f64)
                alpha = -ex*0.5_f64*sim%dt
                ! interpolate_array_disp() has an interface that must be changed
                sim%f_x3x4(i,j,:,l) = sim%interp_x3%interpolate_array_disp( &
                     sim%nc_x3, &
                     sim%f_x3x4(i,j,:,l), &
                     alpha )
             end do
          end do
       end do

       ! Continue with vy...(x4)
       do j=1,loc_sz_x2
          do i=1,loc_sz_x1
             do k=1,sim%mesh4d%num_cells3
                ey    = real(sim%efield_split(i,j),f64)
                alpha = -ey*0.5_f64*sim%dt
                ! interpolate_array_disp() has an interface that must be changed
                sim%f_x3x4(i,j,k,:) = sim%interp_x4%interpolate_array_disp( &
                     sim%nc_x4, &
                     sim%f_x3x4(i,j,k,:), &
                     alpha )
             end do
          end do
       end do
       ! plot fields is giving some errors after the size of the data arrays was
       ! changed to include the last point (the periodic point). After this, some
       ! error messages are given by hdf5...
!       call plot_fields(itime, sim)

    end do ! main loop
  end subroutine run_vp4d_cartesian

  subroutine delete_vp4d_par_cart( sim )
    class(sll_simulation_4d_vlasov_poisson_cart) :: sim
    sll_int32 :: ierr
    SLL_DEALLOCATE( sim%f_x1x2, ierr )
    SLL_DEALLOCATE( sim%f_x3x4, ierr )
    SLL_DEALLOCATE_ARRAY( sim%partial_reduction, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_split, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_x2, ierr )
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
    call delete( sim%interp_x1 )
    call delete( sim%interp_x2 )
    call delete( sim%interp_x3 )
    call delete( sim%interp_x4 )
    SLL_DEALLOCATE_ARRAY( sim%efield_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%efield_x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%efield_split, ierr )
  end subroutine delete_vp4d_par_cart

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
  subroutine compute_charge_density( mesh, numpts1, numpts2, f, partial, rho )
    type(simple_cartesian_4d_mesh), pointer     :: mesh
    sll_real64, intent(in),  dimension(:,:,:,:) :: f       ! local distr. func
    sll_real64, intent(inout),  dimension(:,:,:):: partial ! intermediate res.
    sll_real64, intent(inout), dimension(:,:)     :: rho     ! local rho
    ! local sizes in the split directions have to be given by caller.
    sll_int32, intent(in)                       :: numpts1
    sll_int32, intent(in)                       :: numpts2
    sll_real64                                  :: delta3
    sll_real64                                  :: delta4
    sll_int32                                   :: numpts3
    sll_int32                                   :: numpts4
    sll_int32 :: i, j, k, l
    
    delta4   = mesh%delta_x4
    delta3   = mesh%delta_x3
    partial(:,:,:) = 0.0
    numpts3 = mesh%num_cells3
    numpts4 = mesh%num_cells4
    
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
  
  subroutine compute_electric_field_x1( &
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
       ex = r_delta*(-1.5_f64*phi_x1(1,j) + &
                      2.0_f64*phi_x1(2,j) - &
                      0.5_f64*phi_x1(3,j) )
       efield_x1(1,j) = cmplx(ex,f64)  
       ! right:
       ex = r_delta*(0.5_f64*phi_x1(num_pts_x1-2,j)-&
                     2.0_f64*phi_x1(num_pts_x1-1,j)+&
                     1.5_f64*phi_x1(num_pts_x1,j) )
       efield_x1(num_pts_x1,j) = cmplx(ex,f64) 
    end do
    
    ! Electric field in interior points
    do j=1,num_pts_x2
       do i=2, num_pts_x1-1
          ex = r_delta*0.5_f64*(phi_x1(i+1,j) - phi_x1(i-1,j))
          efield_x1(i,j) = cmplx(ex,f64)
       end do
    end do
  end subroutine compute_electric_field_x1
  
  ! This function only sets the Ey component of the electric field. NOTE: This
  ! and the above function have a terrible flaw: they are not independent. 
  ! compute_electric_field_x2 assumes that compute_electric_field_x1 has been
  ! called first, thus the ex values have been computed. 
  ! compute_electric_field_x1() does not have a similar assumption. It would
  ! be good to change things in order to get rid of these functions.
  subroutine compute_electric_field_x2( &
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
       ey = r_delta*(-1.5_f64*phi_x2(i,1) + 2.0_f64*phi_x2(i,2) - &
                      0.5_f64*phi_x2(i,3))
       efield_x2(i,1) = cmplx(ex,ey,f64)
       ! top:
       ex = real(efield_x2(i,num_pts_x1),f64)
       ey = r_delta*(0.5_f64*phi_x2(i,num_pts_x1-2) - &
                     2.0_f64*phi_x2(i,num_pts_x1-1)+&
                     1.5_f64*phi_x2(i,num_pts_x1))
       efield_x2(i,num_pts_x1) = cmplx(ex,ey,f64) 
    end do
    
    ! Electric field in interior points
    do j=2,num_pts_x2-1
       do i=1, num_pts_x1
          ex = real(efield_x2(i,j),f64)
          ey = r_delta*0.5_f64*(phi_x2(i,j+1) - phi_x2(i,j-1))
          efield_x2(i,j) = cmplx(ex,ey,f64) 
       end do
    end do
  end subroutine compute_electric_field_x2
  
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

  subroutine plot_fields(itime, sim)
    use sll_collective
    use hdf5
    use sll_hdf5_io_parallel
    use sll_xml_io
    sll_int32, intent(in) :: itime
    character(len=4)      :: ctime
    sll_int32             :: i_layout
    character(len=1)      :: c_layout
    class(sll_simulation_4d_vlasov_poisson_cart), intent(in) :: sim
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

    array_dims(1) = sim%nc_x1
    array_dims(2) = sim%nc_x2
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
       
          x1_min = sim%mesh4d%x1_min
          x1_max = sim%mesh4d%x1_max
          x2_min = sim%mesh4d%x2_min
          x2_max = sim%mesh4d%x2_max
          x3_min = sim%mesh4d%x3_min
          x3_max = sim%mesh4d%x3_max
          x4_min = sim%mesh4d%x4_min
          x4_max = sim%mesh4d%x4_max
   
          delta_x1 = sim%mesh4d%delta_x1
          delta_x2 = sim%mesh4d%delta_x2
          delta_x3 = sim%mesh4d%delta_x3
          delta_x4 = sim%mesh4d%delta_x4
   
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

end module sll_simulation_4d_vlasov_poisson_cartesian
