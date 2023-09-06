! ---
!  6D domain decomposition, 6d -> 3d, demo program.
! ---
!  2015
!    Klaus Reuter, khr@mpcdf.mpg.de
!    Max Planck Computing and Data Facility (MPCDF)
!    Garching, Germany
! ---

program test_decomposition
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective
  use sll_m_decomposition
  use mpi

  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! --- variable section

  sll_int32, parameter :: nd = 6
  sll_int32 :: procs_per_dimension(nd)
  logical :: keep_dim(nd), periodic(nd)!, check(2*nd+1)

  type(sll_t_cartesian_topology_6d), pointer :: topology_6d
  type(sll_t_cartesian_topology_3d), pointer :: topology_3d_velocity
  type(sll_t_cartesian_topology_3d), pointer :: topology_3d_spatial
  type(sll_t_collective_t), pointer :: collective_6d
  type(sll_t_collective_t), pointer :: collective_3d_velocity
  type(sll_t_collective_t), pointer :: collective_3d_spatial
  type(sll_t_decomposition_6d), pointer :: decomposition

  sll_int32 :: ierr
  sll_int32 :: world_size, my_rank
  sll_int32 :: global_grid_points_per_dimension(nd)
  sll_int32 :: halo_width_per_dimension(nd)
  sll_real64, dimension(:,:,:,:,:,:), pointer :: f6d
  sll_real64 :: val1, val2
  logical, parameter :: verbose = .false.

  ! --- executable section

  call sll_s_boot_collective()

  world_size = sll_f_get_collective_size(sll_v_world_collective)
  my_rank = sll_f_get_collective_rank(sll_v_world_collective)


  !> Create a distribution of the MPI processes over the dimensions.
  !> There is a convenience function provided by MPI as well, however
  !> we do this by hand for the moment to have better control.
  select case(world_size)
  case(1)
    procs_per_dimension(:) = 1
  case(8)
    procs_per_dimension(1:3) = 2
    procs_per_dimension(4:6) = 1
  case(16)
    procs_per_dimension(1:4) = 2
    procs_per_dimension(5:6) = 1
  case(32)
    procs_per_dimension(1:5) = 2
    procs_per_dimension(6) = 1
  case(64)
    procs_per_dimension(:) = 2
  case default
    write(*,*) "WARNING: No topology implemented for ", world_size, " processes.  Skipping tests."
    write(*,*) "Domain decomposition unit test: PASSED"
    go to 101
  end select



  !> (1) Create a topology based on the distribution.  In addition, a nd boolean array
  !> is passed to indicate if there is a periodic BC associated to a certain dimension.
  periodic(:) = .true.
  topology_6d => &
    sll_f_new_cartesian_topology_6d(sll_v_world_collective, procs_per_dimension, periodic)
  collective_6d => sll_f_create_collective(topology_6d%comm)

  !> Create 3D sub-topologies from the 6D topology, keeping the velocity space 3D topology.
  keep_dim(1:3) = .false.
  keep_dim(4:6) = .true.
  topology_3d_velocity => &
    sll_f_new_cartesian_topology_3d_from_6d(topology_6d, keep_dim)
  !> Derive a sll_collective
  collective_3d_velocity => sll_f_create_collective(topology_3d_velocity%comm)


  !> Create 3D sub-topologies from the 6D topology, keeping the spatial 3D topology.
  keep_dim(1:3) = .true.
  keep_dim(4:6) = .false.
  topology_3d_spatial => &
    sll_f_new_cartesian_topology_3d_from_6d(topology_6d, keep_dim)
  !> Derive a sll_collective
  collective_3d_spatial => sll_f_create_collective(topology_3d_spatial%comm)


  if (all(topology_3d_velocity%coords == 0)) then
    !write(*,*) topology_6d%rank, topology_3d_spatial%coords
  endif


  !   if (verbose) then
  !		 write(*,'(A,I0.2,A,I1,A,I1,A,I1,A,I1,A,I1,A,I1,A,Z0.3,A,I1,A,I1,A,I1,A,I1,A)') &
  !			 "6D[", topology_6d%rank, "](",&
  !			 topology_6d%coords(1), ",", topology_6d%coords(2), ",", topology_6d%coords(3), ",",&
  !			 topology_6d%coords(4), ",", topology_6d%coords(5), ",", topology_6d%coords(6),&
  !			 ")  ==>  3D[", topology_3d%info, ",", topology_3d%rank, "](",&
  !			 topology_3d%coords(1), ",", topology_3d%coords(2), ",", topology_3d%coords(3),&
  !			 ")"
  !   endif


  !go to 101


  !> (2) Create a domain decomposition, create a global array distributed over the
  !> MPI processes.  The size of the local box of the array is passed via an nd array
  !> as well as the widh of the halo to be exchanged between neighboring processes.
  !> The decomposition object contains all the information necessary to allocate and
  !> access the local array (ie it has all the indices), see below.
  global_grid_points_per_dimension(:) = 8
  halo_width_per_dimension(:) = 2
  decomposition => &
    sll_f_new_cartesian_domain_decomposition_6d(topology_6d, global_grid_points_per_dimension, halo_width_per_dimension)



  !> Allocate the local array.  For whatever reason, the SLL__ALLOCATE()
  !> macro does not support multi line arguments ...
  allocate( f6d(decomposition%local%lo(1):decomposition%local%hi(1), &
            decomposition%local%lo(2):decomposition%local%hi(2), &
            decomposition%local%lo(3):decomposition%local%hi(3), &
            decomposition%local%lo(4):decomposition%local%hi(4), &
            decomposition%local%lo(5):decomposition%local%hi(5), &
            decomposition%local%lo(6):decomposition%local%hi(6)),&
            stat=ierr )
  SLL_ASSERT( ierr == 0 )


  !> Fill the array with the global MPI rank.
  f6d = real(my_rank, kind=f64)


  ! At each {x,y,z} position in the 6d Cartesian MPI topology
  ! do a sum over all processes that cover the upper three dimensions,
  ! i.e. sum over the velocity space distribution.
  val1 = real(my_rank, kind=f64)
  call sll_o_collective_globalsum(collective_3d_velocity, val1)
  if (all(topology_3d_velocity%coords == 0)) then
    ! write(*,*) topology_6d%rank, "|", topology_3d_spatial%rank, "|", topology_3d_spatial%coords

    ! write(*,*) topology_6d%rank, "|", topology_3d_spatial%rank, "|", collective_3d_spatial%rank, &
    !                              "|", topology_3d_velocity%rank, "|", collective_3d_velocity%rank

    call sll_o_collective_globalsum(collective_3d_spatial, val1)
    if (collective_3d_spatial%rank == 0) &
      write(*,*) "sum velocity --> sum spatial is", val1
  endif

  ! For comparison, do the direct global sum.
  val2 = real(my_rank, kind=f64)
  call sll_o_collective_globalsum(collective_6d, val2)
  if (collective_6d%rank == 0) &
    write(*,*) "direct sum is", val2



  SLL_DEALLOCATE(topology_6d, ierr)
  SLL_DEALLOCATE(topology_3d_spatial, ierr)
  SLL_DEALLOCATE(topology_3d_velocity, ierr)
  SLL_DEALLOCATE(collective_3d_spatial, ierr)
  SLL_DEALLOCATE(collective_3d_velocity, ierr)
  SLL_DEALLOCATE(decomposition, ierr)
  deallocate( f6d )

  101 if (my_rank == 0) &
    write(*,*) "Domain decomposition 6d-->3d unit test stub : PASSED"
  call sll_s_halt_collective()


end program test_decomposition
