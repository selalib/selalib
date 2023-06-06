! ---
!  6D domain decomposition, 6d -> 3d, remap, Poisson solver, demo program.
! ---
!  2016
!    Katharina Kormann, katharina.kormann@ipp.mpg.de
!      Max Planck Institute for Plasma Physics (IPP)
!      Garching, Germany
!    Klaus Reuter, khr@mpcdf.mpg.de
!      Max Planck Computing and Data Facility (MPCDF)
!      Garching, Germany
! ---

program test_decomposition_6d3d_remap
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective
  use sll_m_decomposition
  use sll_m_remapper
  use sll_m_poisson_3d_periodic_par
  use sll_mpi

  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! --- variable section

  sll_int32, parameter :: nd = 6
  sll_int32 :: procs_per_dimension(nd)
  logical :: keep_dim(nd), periodic(nd), check(2*nd+1)

  type(sll_t_cartesian_topology_6d), pointer :: topology_6d
  type(sll_t_cartesian_topology_3d), pointer :: topology_3d_velocity
  type(sll_t_cartesian_topology_3d), pointer :: topology_3d_spatial
  type(sll_t_collective_t), pointer :: collective_6d
  type(sll_t_collective_t), pointer :: collective_3d_velocity
  type(sll_t_collective_t), pointer :: collective_3d_spatial
  type(sll_t_decomposition_6d), pointer :: decomposition

  ! attempt to wire the domain decomposition with the existing poisson solver
  type(sll_t_layout_3d), pointer :: layout_3d
  type(sll_t_poisson_3d_periodic_par) :: poisson

  sll_int32 :: ierr
  sll_int32 :: world_size, my_rank
  sll_int32 :: global_grid_points_per_dimension(nd)
  sll_int32 :: halo_width_per_dimension(nd)

  sll_real64, allocatable :: f6d(:,:,:,:,:,:)
  sll_real64, allocatable :: rho(:,:,:)
  sll_real64, allocatable :: phi(:,:,:)
  sll_real64, allocatable :: ex(:,:,:)
  sll_real64, allocatable :: ey(:,:,:)
  sll_real64, allocatable :: ez(:,:,:)

  sll_real64 :: val1, val2, volume_v, d_eta
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
    go to 101
  end select


  !> (1) Create a topology based on the distribution.  In addition, a nd boolean array
  !> is passed to indicate if there is a periodic BC associated to a certain dimension.
  periodic(:) = .true.
  topology_6d => &
    sll_f_new_cartesian_topology_6d(sll_v_world_collective, procs_per_dimension, periodic)
  !> derive a sll_collective from the topologies' MPI communicator
  collective_6d => sll_f_create_collective(topology_6d%comm)


  !> Create 3D sub-topologies from the 6D topology, keeping the velocity space 3D topology.
  keep_dim(1:3) = .false.
  keep_dim(4:6) = .true.
  topology_3d_velocity => &
    sll_f_new_cartesian_topology_3d_from_6d(topology_6d, keep_dim)
  !> derive a sll_collective from the topologies' MPI communicator
  collective_3d_velocity => sll_f_create_collective(topology_3d_velocity%comm)


  !> Create 3D sub-topologies from the 6D topology, keeping the spatial 3D topology.
  keep_dim(1:3) = .true.
  keep_dim(4:6) = .false.
  topology_3d_spatial => &
    sll_f_new_cartesian_topology_3d_from_6d(topology_6d, keep_dim)
  !> derive a sll_collective from the topologies' MPI communicator
  collective_3d_spatial => sll_f_create_collective(topology_3d_spatial%comm)



  !> (2) Create a domain decomposition, create a global array distributed over the
  !> MPI processes.  The size of the local box of the array is passed via an nd array
  !> as well as the widh of the halo to be exchanged between neighboring processes.
  !> The decomposition object contains all the information necessary to allocate and
  !> access the local array (ie it has all the indices), see below.
  global_grid_points_per_dimension(:) = 8
  halo_width_per_dimension(:) = 0
  decomposition => &
    sll_f_new_cartesian_domain_decomposition_6d(topology_6d, &
      global_grid_points_per_dimension, &
      halo_width_per_dimension)


  !> Allocate the local array.
  allocate( f6d(decomposition%local%lo(1):decomposition%local%hi(1), &
                decomposition%local%lo(2):decomposition%local%hi(2), &
                decomposition%local%lo(3):decomposition%local%hi(3), &
                decomposition%local%lo(4):decomposition%local%hi(4), &
                decomposition%local%lo(5):decomposition%local%hi(5), &
                decomposition%local%lo(6):decomposition%local%hi(6)),&
                stat=ierr )
  SLL_ASSERT( ierr == 0 )

  allocate( rho(decomposition%local%mn(1):decomposition%local%mx(1), &
                decomposition%local%mn(2):decomposition%local%mx(2), &
                decomposition%local%mn(3):decomposition%local%mx(3)),&
            stat=ierr )
  SLL_ASSERT( ierr == 0 )

  allocate( ex(decomposition%local%mn(1):decomposition%local%mx(1), &
               decomposition%local%mn(2):decomposition%local%mx(2), &
               decomposition%local%mn(3):decomposition%local%mx(3)),&
            stat=ierr )
  SLL_ASSERT( ierr == 0 )

  allocate( ey(decomposition%local%mn(1):decomposition%local%mx(1), &
               decomposition%local%mn(2):decomposition%local%mx(2), &
               decomposition%local%mn(3):decomposition%local%mx(3)),&
            stat=ierr )
  SLL_ASSERT( ierr == 0 )

  allocate( ez(decomposition%local%mn(1):decomposition%local%mx(1), &
               decomposition%local%mn(2):decomposition%local%mx(2), &
               decomposition%local%mn(3):decomposition%local%mx(3)),&
            stat=ierr )
  SLL_ASSERT( ierr == 0 )


  !> Fill the 6D array with the global MPI rank.
  f6d = real(my_rank, kind=f64)
  rho = 0.0_f64
  ex = 0.0_f64
  ey = 0.0_f64
  ez = 0.0_f64

  ! compute the charge density by summing over the velocity dimensions
  ! using the velocity-dimension collectives derived from "topology_3d_velocity"
  volume_v = 1.0_f64  ! dummy volume variable, to be obtained from sll_mesh, later
  call sll_s_compute_charge_density_6d(f6d, decomposition, rho, collective_3d_velocity, volume_v)

  ! TODO : remap the data to use all the available processors

  ! select a subset of the processors to do the Poisson solve step
  if (all(topology_3d_velocity%coords == 0)) then
    ! TODO : Make absolutely sure that the processor at the velocity
    !        topology coordinate (0,0,0) always has the rank 0!
    SLL_ASSERT(collective_3d_velocity%rank == 0)

    layout_3d => sll_f_new_layout_3d( collective_3d_spatial )

    call sll_o_initialize_layout_with_distributed_array( &
      global_grid_points_per_dimension(1), &
      global_grid_points_per_dimension(2), &
      global_grid_points_per_dimension(3), &
      topology_3d_spatial%procs(1), &
      topology_3d_spatial%procs(2), &
      topology_3d_spatial%procs(3), &
      layout_3d)

    d_eta = 1.0_f64  ! dumme d_eta variable, to be obtained from sll_mesh, later
    ! initialize the existing Poisson solver
    call sll_s_poisson_3d_periodic_par_init( &
      layout_3d, &
      global_grid_points_per_dimension(1), &
      global_grid_points_per_dimension(2), &
      global_grid_points_per_dimension(3), &
      d_eta, &
      d_eta, &
      d_eta, &
      poisson )
    !
    allocate( phi(1:poisson%loc_sizes(1,1), &
                  1:poisson%loc_sizes(1,2), &
                  1:poisson%loc_sizes(1,3)),&
                  stat=ierr )
    SLL_ASSERT( ierr == 0 )
    phi = 0.0_f64

    ! ! call the existing Poisson solver
    call sll_s_poisson_3d_periodic_par_solve(poisson, rho, phi)

    call sll_s_poisson_3d_periodic_par_compute_e_from_phi(poisson, phi, &
         ex, ey, ez)

    SLL_DEALLOCATE(layout_3d, ierr)
    SLL_DEALLOCATE(poisson, ierr)
    deallocate( phi )
  endif

  call sll_collective_bcast_3d_real64(collective_3d_velocity, ex,  0)
  call sll_collective_bcast_3d_real64(collective_3d_velocity, ey,  0)
  call sll_collective_bcast_3d_real64(collective_3d_velocity, ez,  0)


  SLL_DEALLOCATE(topology_6d, ierr)
  SLL_DEALLOCATE(topology_3d_spatial, ierr)
  SLL_DEALLOCATE(topology_3d_velocity, ierr)
  SLL_DEALLOCATE(collective_6d, ierr)
  SLL_DEALLOCATE(collective_3d_spatial, ierr)
  SLL_DEALLOCATE(collective_3d_velocity, ierr)
  SLL_DEALLOCATE(decomposition, ierr)
  deallocate(f6d)
  deallocate(rho)
  deallocate(ex)
  deallocate(ey)
  deallocate(ez)

  101 if (my_rank == 0) &
    write(*,*) "Domain decomposition 6d-->3d unit test stub : PASSED"
  call sll_s_halt_collective()

contains


  !> PROTOTYPE : function to derive a "collective" from a plain MPI communicator,
  !> to be propagated into sll_m_collective.F90 once the approach has proven useful.
  function sll_f_create_collective(comm)
    type(sll_t_collective_t), pointer :: sll_f_create_collective
    sll_int32, intent(in) :: comm
    sll_int32 :: ierr

    SLL_ALLOCATE(sll_f_create_collective, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)
    sll_f_create_collective%comm = comm
    sll_f_create_collective%color = 0
    sll_f_create_collective%key = 0
    ! plain MPI calls will be fine inside "sll_m_collective.F90"
    call MPI_COMM_RANK(sll_f_create_collective%comm, sll_f_create_collective%rank, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)
    call MPI_COMM_SIZE(sll_f_create_collective%comm, sll_f_create_collective%size, ierr)
    SLL_ASSERT(ierr == MPI_SUCCESS)
  end function sll_f_create_collective


  ! PROTOTYPE : Compute charge density from a domain-decomposed distribution function.
  subroutine sll_s_compute_charge_density_6d(f_6d, decomp_6d, rho_3d, coll_3d_v, volume_v)
    type(sll_t_decomposition_6d), intent(in) :: decomp_6d
    sll_real64, intent(in)  :: f_6d(decomp_6d%local%lo(1):decomp_6d%local%hi(1), &
                                    decomp_6d%local%lo(2):decomp_6d%local%hi(2), &
                                    decomp_6d%local%lo(3):decomp_6d%local%hi(3), &
                                    decomp_6d%local%lo(4):decomp_6d%local%hi(4), &
                                    decomp_6d%local%lo(5):decomp_6d%local%hi(5), &
                                    decomp_6d%local%lo(6):decomp_6d%local%hi(6))
    type(sll_t_collective_t), pointer, intent(in) :: coll_3d_v
    sll_real64, intent(in)  :: volume_v
    sll_real64, intent(out) :: rho_3d(decomp_6d%local%mn(1):decomp_6d%local%mx(1), &
                                      decomp_6d%local%mn(2):decomp_6d%local%mx(2), &
                                      decomp_6d%local%mn(3):decomp_6d%local%mx(3))
    sll_int32 :: i, j, k

    do k=decomp_6d%local%mn(3), decomp_6d%local%mx(3)
       do j=decomp_6d%local%mn(2), decomp_6d%local%mx(2)
          do i=decomp_6d%local%mn(1), decomp_6d%local%mx(1)
             rho_3d(i,j,k) = volume_v * sum(f_6d(i,j,k,decomp_6d%local%mn(4):decomp_6d%local%mx(4),&
                  decomp_6d%local%mn(5): decomp_6d%local%mx(5),decomp_6d%local%mn(6):decomp_6d%local%mx(6)))
          end do
       end do
    end do
    ! use an allreduce operation until we know where the result needs to go
    call sll_collective_allreduce_sum_3d_real64(coll_3d_v, rho_3d)
  end subroutine sll_s_compute_charge_density_6d


  ! PROTOTYPE : Sum all elements in a 3d buffer.
  ! Update: Removed MPI_IN_PLACE because of problems with OpenMPI.
  subroutine sll_collective_allreduce_sum_3d_real64( col, buffer )
    use sll_mpi
    type(sll_t_collective_t), pointer :: col
    sll_real64, dimension(:,:,:), intent(inout) :: buffer
    sll_real64, dimension(:,:,:), allocatable :: temp
    sll_int32 :: count, ierr

    count = size(buffer)
    allocate(temp(size(buffer,dim=1), size(buffer,dim=2), size(buffer,dim=3)))
    call MPI_ALLREDUCE( buffer, temp, count, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        col%comm, ierr )
    SLL_ASSERT(ierr == MPI_SUCCESS)
    buffer = temp
    deallocate(temp)
  end subroutine sll_collective_allreduce_sum_3d_real64


  ! PROTOTYPE : 3d array broadcast.
  subroutine sll_collective_bcast_3d_real64( col, buffer, root )
    type(sll_t_collective_t), pointer :: col
    sll_real64, dimension(:,:,:), intent(inout) :: buffer
    sll_int32, intent(in) :: root
    sll_int32 :: count, ierr

    count = size(buffer)
    call MPI_BCAST( buffer, count, MPI_DOUBLE_PRECISION, root, col%comm, ierr )
    SLL_ASSERT(ierr == MPI_SUCCESS)
  end subroutine sll_collective_bcast_3d_real64

end program test_decomposition_6d3d_remap
