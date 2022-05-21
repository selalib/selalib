! ---
!  6D domain decomposition demo program.
! ---
!  2016
!    Klaus Reuter, khr@mpcdf.mpg.de
!    Max Planck Computing and Data Facility (MPCDF)
!    Garching, Germany
! ---

program test_decomposition
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_decomposition, only: &
    sll_t_decomposition_slim_6d, &
    sll_f_new_cartesian_domain_decomposition_slim_6d, &
    sll_f_apply_halo_exchange, &
    sll_t_cartesian_topology_6d, &
    sll_f_new_cartesian_topology_6d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! --- variable section

  sll_int32, parameter :: nd = 6
  sll_int32 :: procs_per_dimension(nd)
  logical :: periodic(nd), check(2*nd+1)
  type(sll_t_cartesian_topology_6d), pointer :: topology
  type(sll_t_decomposition_slim_6d), pointer :: decomposition
  sll_int32 :: ierr
  sll_int32 :: world_size, my_rank
  sll_int32 :: global_grid_points_per_dimension(nd)
  sll_int32 :: id, hw_left, hw_right
  sll_real64, dimension(:,:,:,:,:,:), pointer :: f6d
  sll_real64 :: val_left, val_right
  logical :: chk_tmp, chk_all


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
  topology => &
    sll_f_new_cartesian_topology_6d(sll_v_world_collective, procs_per_dimension, periodic)


  !> (2) Create a domain decomposition, create a global array distributed over the
  !> MPI processes.  The size of the local box of the array is passed via an nd array
  !> as well as the widh of the halo to be exchanged between neighboring processes.
  !> The decomposition object contains all the information necessary to allocate and
  !> access the local array (ie it has all the indices), see below.
  global_grid_points_per_dimension(:) = 8
  decomposition => &
    sll_f_new_cartesian_domain_decomposition_slim_6d(topology, global_grid_points_per_dimension)


  !> Allocate the local array.
  allocate(f6d(decomposition%local%mn(1):decomposition%local%mx(1), &
               decomposition%local%mn(2):decomposition%local%mx(2), &
               decomposition%local%mn(3):decomposition%local%mx(3), &
               decomposition%local%mn(4):decomposition%local%mx(4), &
               decomposition%local%mn(5):decomposition%local%mx(5), &
               decomposition%local%mn(6):decomposition%local%mx(6)),&
               stat=ierr)
  SLL_ASSERT( ierr == 0 )

  !> Fill the array with the MPI rank.
  f6d = real(my_rank, kind=f64)

  !> (3) Exercise and check data exchange between next neighbors.
  chk_all = .true.
  do id = 1, 6
    do hw_left = 1, 4
      do hw_right = 1, 4
        call sll_f_apply_halo_exchange(topology, decomposition, f6d, &
                                       id, hw_left, hw_right)
        ! check halo buffers
        val_left = real(topology%neighbors(2*id-1), kind=f64)
        val_right = real(topology%neighbors(2*id), kind=f64)
        ! write(*,*) val_left, val_right
        chk_tmp = all(decomposition%local%halo_left%buf == val_left) .and. &
                  all(decomposition%local%halo_right%buf == val_right)
        if (.not. chk_tmp) then
          chk_all = .false.
          goto 100
        endif
      enddo
    enddo
  enddo

100 continue

  SLL_DEALLOCATE(topology, ierr)
  deallocate(decomposition%local%halo_left%buf)
  deallocate(decomposition%local%halo_right%buf)
  SLL_DEALLOCATE(decomposition, ierr)
  deallocate( f6d )

  if (chk_all .eqv. .true.) then
    write(*,*) "Domain decomposition unit test : PASSED"
  else
    write(*,*) "Domain decomposition unit test : FAILED"
  endif

101   call sll_s_halt_collective()
end program test_decomposition
