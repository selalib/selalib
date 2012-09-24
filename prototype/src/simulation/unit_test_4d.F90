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

  sll_int32 :: world_size
  sll_int32 :: my_rank
  sll_int32 :: nproc_x1, nproc_x2, nproc_x3, nproc_x4 
  sll_int32 :: nc_x1, nc_x2, nc_x3, nc_x4
  sll_int32 :: loc_sz_x1
  sll_int32 :: loc_sz_x2
  sll_int32 :: loc_sz_x3
  sll_int32 :: loc_sz_x4
  ! distribution functions. There are several because each array represents
  ! a differently shaped chunk of memory. In this example, each chunk allows   
  ! sequential operations in one given direction. f1 should permit to carry
  ! out sequential operations in x1, f2 in x2, etc.
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
  type(remap_plan_4D), pointer :: seqx1x2_to_seqx3x4
  type(remap_plan_4D), pointer :: seqx3x4_to_seqx1x2
  sll_int32 :: power2 ! 2^power2 = number of processes available
  sll_int32 :: itemp

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
  nc_x3 = 16
  nc_x4 = 16

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

  call initialize_layout_with_distributed_2D_array( &
       nc_x1, &
       nc_x2, &
       world_size, &
       1, &
       rho_seq_x2 )
  SLL_ALLOCATE(rho_x2(loc_sz_x1,loc_sz_x2),ierr)
  SLL_ALLOCATE(phi_x2(loc_sz_x1,loc_sz_x2),ierr)

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

  ! Proceed to carry out the advections.


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

end program vlasov_poisson_4d


