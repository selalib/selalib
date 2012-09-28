! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D cartesian: x, y, vx, vy (or x1, x2, x3, x4)
! - parallel

program vlasov_poisson_4d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_remap.h"
  use sll_collective
  use numeric_constants
  use sll_cubic_spline_interpolator_1d
  use sll_test_4d_initializer
  use sll_poisson_2d_periodic_cartesian_par
  use sll_cubic_spline_interpolator_1d
  use sll_electric_field_2d_accumulator
  implicit none

  type(cubic_spline_1d_interpolator), target :: interp_x
  type(cubic_spline_1d_interpolator), target :: interp_v
  type(cubic_spline_1d_interpolator), target :: interp_vx
  type(cubic_spline_1d_interpolator), target :: interp_vy
  class(sll_interpolator_1d_base), pointer   :: interp_x_ptr
  class(sll_interpolator_1d_base), pointer   :: interp_v_ptr
  class(sll_interpolator_1d_base), pointer   :: interp_vx_ptr
  class(sll_interpolator_1d_base), pointer   :: interp_vy_ptr
  ! for initializers
  type(init_test_4d_par)                     :: init_4d
  type(simple_cartesian_4d_mesh), pointer    :: mesh4d
  type(poisson_2d_periodic_plan_cartesian_par), pointer :: poisson_plan

  sll_int32  :: world_size
  sll_int32  :: my_rank
  sll_int32  :: nproc_x1, nproc_x2, nproc_x3, nproc_x4 
  sll_int32  :: nc_x1, nc_x2, nc_x3, nc_x4
  sll_int32  :: loc_sz_x1
  sll_int32  :: loc_sz_x2
  sll_int32  :: loc_sz_x3
  sll_int32  :: loc_sz_x4
  sll_int32  :: i
  sll_int32  :: j
  sll_int32  :: k
  sll_int32  :: l
  sll_real64 :: alpha
  sll_real64 :: dt
  ! distribution functions. There are several because each array represents
  ! a differently shaped chunk of memory. In this example, each chunk allows   
  ! sequential operations in one given direction. f1 should permit to carry
  ! out sequential operations in x1, f2 in x2, etc.
  sll_real64, dimension(:), pointer           :: line
  sll_real64, dimension(:,:,:,:), allocatable :: f_x1x2 
  sll_real64, dimension(:,:,:,:), allocatable :: f_x3x4
  sll_real64, dimension(:,:,:), allocatable   :: partial_reduction
  sll_real64, dimension(:,:), allocatable     :: rho_x1 
  sll_real64, dimension(:,:), allocatable     :: rho_x2 
  sll_real64, dimension(:,:), allocatable     :: rho_split
  sll_real64, dimension(:,:), allocatable     :: phi_x1
  sll_real64, dimension(:,:), allocatable     :: phi_x2
  sll_int32 :: ierr

  ! for remap
  type(layout_4D), pointer :: sequential_x1x2
  type(layout_4D), pointer :: sequential_x3x4
  type(layout_2D), pointer :: rho_seq_x1
  type(layout_2D), pointer :: rho_seq_x2
  type(layout_2D), pointer :: split_rho_layout ! layout that is not sequential 
  type(remap_plan_2D), pointer :: split_to_seqx1
  type(remap_plan_2D), pointer :: seqx1_to_seqx2
  ! remaps for the electric field data
  type(remap_plan_2D), pointer :: efld_split_to_seqx1
  type(remap_plan_2D), pointer :: efld_seqx1_to_seqx2
  type(remap_plan_2D), pointer :: efld_seqx2_to_split

  type(remap_plan_4D), pointer :: seqx1x2_to_seqx3x4
  type(remap_plan_4D), pointer :: seqx3x4_to_seqx1x2

  ! interpolators and their pointers
  type(cubic_spline_1d_interpolator), target :: interp_x1
  type(cubic_spline_1d_interpolator), target :: interp_x2
  type(cubic_spline_1d_interpolator), target :: interp_x3
  type(cubic_spline_1d_interpolator), target :: interp_x4
  class(sll_interpolator_1d_base), pointer   :: interp_x1_ptr
  class(sll_interpolator_1d_base), pointer   :: interp_x2_ptr
  class(sll_interpolator_1d_base), pointer   :: interp_x3_ptr
  class(sll_interpolator_1d_base), pointer   :: interp_x4_ptr

  ! Field accumulator
  type(efield_2d_point), dimension(:,:), allocatable :: efield_x1
  type(efield_2d_point), dimension(:,:), allocatable :: efield_x2
  type(efield_2d_point), dimension(:,:), allocatable :: efield_split

  sll_int32 :: power2 ! 2^power2 = number of processes available
  sll_int32 :: itemp

  dt = 0.01

  call sll_boot_collective() ! this should be inside another function
  world_size = sll_get_collective_size(sll_world_collective)
  my_rank    = sll_get_collective_rank(sll_world_collective)

  ! allocate the layouts...
  sequential_x1x2  => new_layout_4D( sll_world_collective )
  sequential_x3x4  => new_layout_4D( sll_world_collective )
  rho_seq_x1       => new_layout_2D( sll_world_collective )
  rho_seq_x2       => new_layout_2D( sll_world_collective )
  split_rho_layout => new_layout_2D( sll_world_collective )

  ! In this particular simulation, since the system is periodic, the number
  ! of points is the same as the number of cells in all directions.
  nc_x1 = 32 
  nc_x2 = 32
  nc_x3 = 32
  nc_x4 = 32

  ! layout for sequential operations in x3 and x4. Make an even split for
  ! x1 and x2, or as close as even if the power of 2 is odd. This should 
  ! be packaged in some sort of routine.
  power2 = int(log(real(world_size))/log(2.0))

  ! special case N = 1, so power2 = 0
  if(power2 == 0) then
     nproc_x1 = 1
     nproc_x2 = 1
     nproc_x3 = 1
     nproc_x4 = 1
  end if

  if(is_even(power2)) then
     nproc_x1 = 2**(power2/2)
     nproc_x2 = 2**(power2/2)
     nproc_x3 = 1
     nproc_x4 = 1
  else 
     nproc_x1 = 2**((power2-1)/2)
     nproc_x2 = 2**((power2+1)/2)
     nproc_x3 = 1
     nproc_x4 = 1
  end if

  call initialize_layout_with_distributed_4D_array( &
       nc_x1, &
       nc_x2, &
       nc_x3, &
       nc_x4, &
       nproc_x1, &
       nproc_x2, &
       nproc_x3, &
       nproc_x4, &
       sequential_x3x4 )

  ! Use this information to initialize the layout that describes the result
  ! of computing rho. This layout is not useful to do sequential operations
  ! in any of the two available directions. We also initialize the other two
  ! layouts needed for both sequential operations on x1 and x2 in the 2D case.
  call initialize_layout_with_distributed_2D_array( &
       nc_x1, &
       nc_x2, &
       nproc_x1, &
       nproc_x2, &
       split_rho_layout )

  call initialize_layout_with_distributed_2D_array( &
       nc_x1, &
       nc_x2, &
       1, &
       world_size, &
       rho_seq_x1 )

  call compute_local_sizes_2d( rho_seq_x1, loc_sz_x1, loc_sz_x2 )
  SLL_ALLOCATE(rho_x1(loc_sz_x1,loc_sz_x2),ierr)
  SLL_ALLOCATE(phi_x1(loc_sz_x1,loc_sz_x2),ierr)
  ! Experiment with a dedicated array to store the values of the electric
  ! field in each point of the grid.
  SLL_ALLOCATE(efield_x1(loc_sz_x1,loc_sz_x2),ierr)

  call initialize_layout_with_distributed_2D_array( &
       nc_x1, &
       nc_x2, &
       world_size, &
       1, &
       rho_seq_x2 )
  call compute_local_sizes_2d( rho_seq_x2, loc_sz_x1, loc_sz_x2 )
  SLL_ALLOCATE(rho_x2(loc_sz_x1,loc_sz_x2),ierr)
  SLL_ALLOCATE(phi_x2(loc_sz_x1,loc_sz_x2),ierr)
  SLL_ALLOCATE(efield_x2(loc_sz_x1,loc_sz_x2),ierr)

  ! layout for sequential operations in x1 and x2. This is basically just the
  ! flipping of the values between x1,x3 and x2,x4 on the previous layout.

  ! switch x1 and x3:
  itemp = nproc_x3
  nproc_x3 = nproc_x1
  nproc_x1 = itemp
  ! switch x2 and x4
  itemp = nproc_x4
  nproc_x4 = nproc_x2 
  nproc_x2 = itemp

  call initialize_layout_with_distributed_4D_array( &
       nc_x1, &
       nc_x2, &
       nc_x3, &
       nc_x4, &
       nproc_x1, &
       nproc_x2, &
       nproc_x3, &
       nproc_x4, &
       sequential_x1x2 )

  ! Allocate the array needed to store the local chunk of the distribution
  ! function data. First compute the local sizes. Since the remap operations
  ! are out-of-place, we will allocate four different arrays, one for each
  ! layout.
  call compute_local_sizes_4d( sequential_x1x2, &
                               loc_sz_x1, &
                               loc_sz_x2, &
                               loc_sz_x3, &
                               loc_sz_x4 )
  SLL_ALLOCATE(f_x1x2(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)

  ! This layout is also useful to represent the charge density array. Since
  ! this is a result of a local reduction on x3 and x4, the new layout is
  ! 2D but with the same dimensions of the process mesh in x1 and x2.
  SLL_ALLOCATE(rho_split(loc_sz_x1,loc_sz_x2),ierr)
  SLL_ALLOCATE(efield_split(loc_sz_x1,loc_sz_x2),ierr)

  call compute_local_sizes_4d( sequential_x3x4, &
                               loc_sz_x1, &
                               loc_sz_x2, &
                               loc_sz_x3, &
                               loc_sz_x4 )
  SLL_ALLOCATE(f_x3x4(loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4),ierr)

  ! These dimensions are also the ones needed for the array where we store
  ! the intermediate results of the charge density computation.
  SLL_ALLOCATE(partial_reduction(loc_sz_x1,loc_sz_x2, loc_sz_x3),ierr)

  ! Initialize the initial distribution function data. We do this with an
  ! initializer object which needs to be initialized itself! Note also that the
  ! mesh is described in global terms. This should be fine as meshes are
  ! supposed to be read-only entities.
  mesh4d => new_cartesian_4d_mesh( nc_x1, nc_x2, nc_x3, nc_x4, &
                                   0.0_f64, 1.0_f64, &
                                   0.0_f64, 1.0_f64, &
                                  -1.0_f64, 1.0_f64, &
                                  -1.0_f64, 1.0_f64 )

  call load_test_4d_initializer( init_4d, &
                                 NODE_CENTERED_FIELD, &
                                 mesh4d, &
                                 0.1_f64, &
                                 sequential_x3x4 )
  call init_4d%f_of_4args(f_x3x4)
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
  rho_split(:,:) = 0.0

  call compute_charge_density( &
       mesh4d, &
       size(f_x3x4,1), &
       size(f_x3x4,2), &
       f_x3x4, &
       partial_reduction, &
       rho_split )

  ! Re-arrange rho_split in a way that permits sequential operations in x1, to
  ! feed to the Poisson solver.
  split_to_seqx1 => NEW_REMAP_PLAN_2D( split_rho_layout, rho_seq_x1, rho_split )
  call apply_remap_2D( split_to_seqx1, rho_split, rho_x1 )

  ! We are in a position now to compute the electric potential.
  ! Initialize the poisson plan
  poisson_plan => new_poisson_2d_periodic_plan_cartesian_par( &
       rho_seq_x1, &
       nc_x1, &
       nc_x2, &
       1.0_f64, &    ! parametrize with mesh values
       1.0_f64 )     ! parametrize with mesh values

  ! solve for the electric potential
  call solve_poisson_2d_periodic_cartesian_par(poisson_plan, rho_x1, phi_x1)

  ! compute the values of the electric field. rho is configured for 
  ! sequential operations in x1, thus we start by computing the E_x component.
  ! The following call is inefficient and unnecessary. The local sizes for
  ! the arrays should be kept around as parameters basically and not on 
  ! variables whose content could be anything... This will have to do for now.
  call compute_local_sizes_2d( rho_seq_x1, loc_sz_x1, loc_sz_x2 )
  call compute_electric_field_x1( &
       phi_x1, &
       loc_sz_x1, &
       loc_sz_x2, &
       mesh4d%delta_x1, &
       efield_x1 )
  ! note that we are 'recycling' the layouts used for the other arrays because
  ! they represent an identical configuration.
  efld_seqx1_to_seqx2 => NEW_REMAP_PLAN_2D( rho_seq_x1, rho_seq_x2, efield_x1)
  seqx1_to_seqx2 => NEW_REMAP_PLAN_2D( rho_seq_x1, rho_seq_x2, phi_x1 )
  call apply_remap_2D( efld_seqx1_to_seqx2, efield_x1, efield_x2 )
  call apply_remap_2D( seqx1_to_seqx2, phi_x1, phi_x2 )
  call compute_local_sizes_2d( rho_seq_x2, loc_sz_x1, loc_sz_x2 )
  call compute_electric_field_x2( &
       phi_x2, &
       loc_sz_x1, &
       loc_sz_x2, &
       mesh4d%delta_x2, &
       efield_x2 )
  ! But now, to make the electric field data configuration compatible with
  ! the sequential operations in x2x3 we need still another remap operation.
  efld_seqx2_to_split => &
       NEW_REMAP_PLAN_2D( rho_seq_x2, split_rho_layout, efield_x2 )
  call apply_remap_2D( efld_seqx2_to_split, efield_x2, efield_split )
  ! Now we proceed to reconfigure the data. This is very expensive. There might
  ! be advantages to this approach if we avoid larger data transfers like with 
  ! an all-to-all transfer... however, we could end up paying more if the 
  ! simulation is latency-dominated.

  ! Proceed to carry out the advections. The following should go inside a
  ! subroutine...

  ! Start the interpolators... Watch out: the periodic case has equal number
  ! of cells than points. Is this properly handled by the interpolators??
  call interp_x1%initialize(nc_x1, mesh4d%x1_min, mesh4d%x1_max,PERIODIC_SPLINE)
  interp_x1_ptr => interp_x1
  call interp_x2%initialize(nc_x2, mesh4d%x2_min, mesh4d%x2_max,PERIODIC_SPLINE)
  interp_x2_ptr => interp_x2
  call interp_x3%initialize(nc_x3, mesh4d%x3_min, mesh4d%x3_max,HERMITE_SPLINE)
  interp_x3_ptr => interp_x3
  call interp_x4%initialize(nc_x4, mesh4d%x4_min, mesh4d%x4_max,HERMITE_SPLINE)
  interp_x4_ptr => interp_x4

  !**************   Pierre: don write between these asterisks     ***********
  call compute_local_sizes_2d( rho_seq_x1, loc_sz_x1, loc_sz_x2 )
  do j=1,loc_sz_x2
     do i=1,loc_sz_x1
        do k=1,mesh4d%num_cells3
           do l=1,mesh4d%num_cells4
              alpha = efield_split(i,j)*0.5_f64*dt
              line => f_x3x4(i,j,k,:)
              ! review with Eric's version... to 
           end do
        end do
     end do
  end do
  !***************************************************************************

  print *, 'reached end of vp4d test'
  print *, 'PASSED'
  call sll_halt_collective()


  contains

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
      delta4   = mesh%delta_x3
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

      sll_real64, dimension(:,:), intent(in)  :: phi_x1
      sll_int32, intent(in)                   :: num_pts_x1
      sll_int32, intent(in)                   :: num_pts_x2
      sll_real64, intent(in)                  :: delta_x1
      type(efield_2d_point), dimension(:,:), intent(out) :: efield_x1
      sll_int32                               :: i
      sll_int32                               :: j
      sll_real64                              :: r_delta
      ! FIXME: arg checking

      r_delta = 1.0_f64/delta_x1

      ! Compute the electric field values on the left and right edges.
      do j=1,num_pts_x2
         ! left:
         efield_x1(1,j)%ex = r_delta*( -1.5_f64*phi_x1(1,j) + &
                                        2.0_f64*phi_x1(2,j) - &
                                        0.5_f64*phi_x1(3,j) )
         ! right:
         efield_x1(num_pts_x1,j)%ex = r_delta*(0.5_f64*phi_x1(num_pts_x1-2,j)-&
                                               2.0_f64*phi_x1(num_pts_x1-1,j)+&
                                               1.5_f64*phi_x1(num_pts_x1,j) )
      end do

      ! Electric field in interior points
      do j=1,num_pts_x2
         do i=2, num_pts_x1-1
            efield_x1(i,j)%ex = r_delta*0.5_f64*(phi_x1(i+1,j) - phi_x1(i-1,j))
         end do
      end do
    end subroutine compute_electric_field_x1

    subroutine compute_electric_field_x2( &
      phi_x2, &
      num_pts_x1, &
      num_pts_x2, &
      delta_x2, &
      efield_x2 )

      sll_real64, dimension(:,:), intent(in)  :: phi_x2
      sll_int32, intent(in)                   :: num_pts_x1
      sll_int32, intent(in)                   :: num_pts_x2
      sll_real64, intent(in)                  :: delta_x2
      type(efield_2d_point), dimension(:,:), intent(out) :: efield_x2
      sll_int32                               :: i
      sll_int32                               :: j
      sll_real64                              :: r_delta

      ! FIXME: arg checking

      r_delta = 1.0_f64/delta_x2

      ! Compute the electric field values on the bottom and top edges.
      do i=1,num_pts_x1
         ! bottom:
         efield_x2(i,1)%ey = r_delta*(-1.5_f64*phi_x2(i,1) + &
                                       2.0_f64*phi_x2(i,2) - &
                                       0.5_f64*phi_x2(i,3))
         ! top:
         efield_x2(i,num_pts_x1)%ey = r_delta*(0.5_f64*phi_x2(i,num_pts_x1-2)-&
                                               2.0_f64*phi_x2(i,num_pts_x1-1)+&
                                               1.5_f64*phi_x2(i,num_pts_x1))
      end do

      ! Electric field in interior points
      do j=2,num_pts_x2-1
         do i=1, num_pts_x1
            efield_x2(i,j)%ey = r_delta*0.5_f64*(phi_x1(i,j+1) - phi_x1(i,j-1))
         end do
      end do
    end subroutine compute_electric_field_x2


end program vlasov_poisson_4d


