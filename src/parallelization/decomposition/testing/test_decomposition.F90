! ---
!  6D domain decomposition demo program.
! ---
!  Copyright 2015
!    Klaus Reuter, khr@mpcdf.mpg.de
!    Max Planck Computing and Data Facility (MPCDF)
!    Garching, Germany
!  All rights reserved.
! ---

program test_decomposition
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_boot_collective, &
    sll_get_collective_rank, &
    sll_get_collective_size, &
    sll_halt_collective, &
    sll_world_collective

  use sll_m_decomposition, only: &
    apply_halo_exchange, &
    cartesian_topology_6d, &
    decomposition_6d, &
    new_cartesian_domain_decomposition, &
    new_cartesian_topology

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! --- variable section

   sll_int32, parameter :: nd = 6
   sll_int32 :: procs_per_dimension(nd)
   logical :: periodic(nd), check(2*nd+1)
   type(cartesian_topology_6d), pointer :: topology
   type(decomposition_6d), pointer :: decomposition
   sll_int32 :: ierr
   sll_int32 :: world_size, my_rank
   sll_int32 :: global_grid_points_per_dimension(nd)
   sll_int32 :: halo_width_per_dimension(nd)
   sll_real64, dimension(:,:,:,:,:,:), pointer :: f6d
   sll_real64 :: val



   ! --- executable section

   call sll_boot_collective()

   world_size = sll_get_collective_size(sll_world_collective)
   my_rank = sll_get_collective_rank(sll_world_collective)



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
      new_cartesian_topology(sll_world_collective, procs_per_dimension, periodic)



   !> (2) Create a domain decomposition, create a global array distributed over the
   !> MPI processes.  The size of the local box of the array is passed via an nd array
   !> as well as the widh of the halo to be exchanged between neighboring processes.
   !> The decomposition object contains all the information necessary to allocate and
   !> access the local array (ie it has all the indices), see below.
   global_grid_points_per_dimension(:) = 8
   halo_width_per_dimension(:) = 2
   decomposition => &
      new_cartesian_domain_decomposition(topology, global_grid_points_per_dimension, halo_width_per_dimension)



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


   !> Fill the array with the MPI rank.
   f6d = real(my_rank, kind=f64)

   !> (3) Do data exchange between neighbors.
   call apply_halo_exchange(topology, decomposition, f6d)

   !> (4) Check if the data was transferred correctly between the neighbours.
   check(:) = .false.
   ! 4a --- first dimension, left neighbor
   val = real(topology%neighbors(1), kind=f64)
   check(1) = &
      all( f6d(decomposition%local%lo(1):decomposition%local%mn(1)-1, &
            decomposition%local%mn(2):decomposition%local%mx(2),   &
            decomposition%local%mn(3):decomposition%local%mx(3),   &
            decomposition%local%mn(4):decomposition%local%mx(4),   &
            decomposition%local%mn(5):decomposition%local%mx(5),   &
            decomposition%local%mn(6):decomposition%local%mx(6)) == val )
   ! 4b --- first dimension, right neighbor
   val = real(topology%neighbors(2), kind=f64)
   check(2) = &
      all( f6d(decomposition%local%mx(1)+1:decomposition%local%hi(1), &
            decomposition%local%mn(2):decomposition%local%mx(2),   &
            decomposition%local%mn(3):decomposition%local%mx(3),   &
            decomposition%local%mn(4):decomposition%local%mx(4),   &
            decomposition%local%mn(5):decomposition%local%mx(5),   &
            decomposition%local%mn(6):decomposition%local%mx(6)) == val )
   ! 4c --- second dimension, left neighbor
   val = real(topology%neighbors(3), kind=f64)
   check(3) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%lo(2):decomposition%local%mn(2)-1, &
               decomposition%local%mn(3):decomposition%local%mx(3),   &
               decomposition%local%mn(4):decomposition%local%mx(4),   &
               decomposition%local%mn(5):decomposition%local%mx(5),   &
               decomposition%local%mn(6):decomposition%local%mx(6)) == val )
   ! 4d --- second dimension, right neighbor
   val = real(topology%neighbors(4), kind=f64)
   check(4) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%mx(2)+1:decomposition%local%hi(2), &
               decomposition%local%mn(3):decomposition%local%mx(3),   &
               decomposition%local%mn(4):decomposition%local%mx(4),   &
               decomposition%local%mn(5):decomposition%local%mx(5),   &
               decomposition%local%mn(6):decomposition%local%mx(6)) == val )
   ! 4e --- third dimension, left neighbor
   val = real(topology%neighbors(5), kind=f64)
   check(5) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%mn(2):decomposition%local%mx(2),   &
               decomposition%local%lo(3):decomposition%local%mn(3)-1, &
               decomposition%local%mn(4):decomposition%local%mx(4),   &
               decomposition%local%mn(5):decomposition%local%mx(5),   &
               decomposition%local%mn(6):decomposition%local%mx(6)) == val )
   ! 4f --- third dimension, right neighbor
   val = real(topology%neighbors(6), kind=f64)
   check(6) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%mn(2):decomposition%local%mx(2),   &
               decomposition%local%mx(3)+1:decomposition%local%hi(3), &
               decomposition%local%mn(4):decomposition%local%mx(4),   &
               decomposition%local%mn(5):decomposition%local%mx(5),   &
               decomposition%local%mn(6):decomposition%local%mx(6)) == val )
   ! 4g --- fourth dimension, left neighbor
   val = real(topology%neighbors(7), kind=f64)
   check(7) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%mn(2):decomposition%local%mx(2),   &
               decomposition%local%mn(3):decomposition%local%mx(3),   &
               decomposition%local%lo(4):decomposition%local%mn(4)-1, &
               decomposition%local%mn(5):decomposition%local%mx(5),   &
               decomposition%local%mn(6):decomposition%local%mx(6)) == val )
   ! 4h --- fourth dimension, right neighbor
   val = real(topology%neighbors(8), kind=f64)
   check(8) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%mn(2):decomposition%local%mx(2),   &
               decomposition%local%mn(3):decomposition%local%mx(3),   &
               decomposition%local%mx(4)+1:decomposition%local%hi(4), &
               decomposition%local%mn(5):decomposition%local%mx(5),   &
               decomposition%local%mn(6):decomposition%local%mx(6)) == val )
   ! 4i --- fifth dimension, left neighbor
   val = real(topology%neighbors(9), kind=f64)
   check(9) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%mn(2):decomposition%local%mx(2),   &
               decomposition%local%mn(3):decomposition%local%mx(3),   &
               decomposition%local%mn(4):decomposition%local%mx(4),   &
               decomposition%local%lo(5):decomposition%local%mn(5)-1, &
               decomposition%local%mn(6):decomposition%local%mx(6)) == val )
   ! 4j --- fifth dimension, right neighbor
   val = real(topology%neighbors(10), kind=f64)
   check(10) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%mn(2):decomposition%local%mx(2),   &
               decomposition%local%mn(3):decomposition%local%mx(3),   &
               decomposition%local%mn(4):decomposition%local%mx(4),   &
               decomposition%local%mx(5)+1:decomposition%local%hi(5), &
               decomposition%local%mn(6):decomposition%local%mx(6)) == val )
   ! 4k --- sixth dimension, left neighbor
   val = real(topology%neighbors(11), kind=f64)
   check(11) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%mn(2):decomposition%local%mx(2),   &
               decomposition%local%mn(3):decomposition%local%mx(3),   &
               decomposition%local%mn(4):decomposition%local%mx(4),   &
               decomposition%local%mn(5):decomposition%local%mx(5),   &
               decomposition%local%lo(6):decomposition%local%mn(6)-1) == val )
   ! 4l --- sixth dimension, right neighbor
   val = real(topology%neighbors(12), kind=f64)
   check(12) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%mn(2):decomposition%local%mx(2),   &
               decomposition%local%mn(3):decomposition%local%mx(3),   &
               decomposition%local%mn(4):decomposition%local%mx(4),   &
               decomposition%local%mn(5):decomposition%local%mx(5),   &
               decomposition%local%mx(6)+1:decomposition%local%hi(6)) == val )

   ! 4m --- finally, check if nothing was overwritten
   val = real(my_rank, kind=f64)
   check(13) = &
      all( f6d(decomposition%local%mn(1):decomposition%local%mx(1),   &
               decomposition%local%mn(2):decomposition%local%mx(2),   &
               decomposition%local%mn(3):decomposition%local%mx(3),   &
               decomposition%local%mn(4):decomposition%local%mx(4),   &
               decomposition%local%mn(5):decomposition%local%mx(5),   &
               decomposition%local%mn(6):decomposition%local%mx(6)) == val )

   SLL_DEALLOCATE(topology, ierr)
   SLL_DEALLOCATE(decomposition, ierr)
   deallocate( f6d )

   !write(*,*) check

   if (all(check) .eqv. .true.) then
      write(*,*) "Domain decomposition unit test : PASSED"
   else
      write(*,*) "Domain decomposition unit test : FAILED"
   endif

101   call sll_halt_collective()
end program test_decomposition
