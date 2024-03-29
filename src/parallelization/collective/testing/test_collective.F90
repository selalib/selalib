program collective_test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_collective, only: &
      sll_s_boot_collective, &
      sll_o_collective_allgather, &
      sll_o_collective_allgatherv, &
      sll_s_collective_allreduce_logical, &
      sll_s_collective_allreduce_real32, &
      sll_s_collective_alltoall_int, &
      sll_s_collective_alltoallv_int_simple, &
      sll_s_collective_alltoallv_real, &
      sll_s_collective_barrier, &
      sll_o_collective_bcast, &
      sll_o_collective_gather, &
      sll_s_collective_gatherv_real, &
      sll_o_collective_globalsum, &
      sll_s_collective_reduce_int, &
      sll_s_collective_reduce_logical, &
      sll_s_collective_reduce_real32, &
      sll_o_collective_scatter, &
      sll_s_collective_scatterv_real, &
      sll_f_get_collective_rank, &
      sll_f_get_collective_size, &
      sll_s_halt_collective, &
      sll_v_world_collective

   use mpi, only: &
      mpi_land, &
      mpi_sum

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! NOTE: some poor choices were made when implementing this test. For example,
   ! the sendbuf, recvbuf, etc. arrays are continually allocated and
   ! deallocated throughout for each individual function test. While this is
   ! OK and probably necessary, things get complicated when doing these
   ! allocations/deallocations by hand and without documenting why some choices
   ! of sizes were made. Some compilers, like gfortran do not seem to care if
   ! some of these arrays are not allocated if they don't play an active role,
   ! like the send buffers in the receiving process of a gather operation, but
   ! an Intel compiler will complain...
   sll_int32 :: rank, size, i, ierr
   LOGICAL, DIMENSION(1) :: logic, logic2
   sll_real32, ALLOCATABLE, DIMENSION(:) :: somme

   sll_real32, allocatable, dimension(:) :: sendbuf_real, recvbuf_real
   sll_int32, allocatable, dimension(:) :: sendbuf_int!, recvbuf_int
   sll_int32, pointer, dimension(:) :: recvbuf_int
   logical, allocatable, dimension(:) :: sendbuf_log, recvbuf_log
   sll_int32, allocatable, dimension(:) :: sendcounts, recvcounts
   sll_int32, allocatable, dimension(:) :: sdispls, rdispls

   call sll_s_boot_collective()

   rank = sll_f_get_collective_rank(sll_v_world_collective)
   size = sll_f_get_collective_size(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_real(1), ierr)
   sendbuf_real(:) = real(rank, f32)
   SLL_ALLOCATE(recvbuf_real(1), ierr)

   call sll_s_collective_reduce_real32(sll_v_world_collective, sendbuf_real, &
                                       1, MPI_SUM, 0, recvbuf_real)

   if (rank == 0) then
      if (recvbuf_real(1) .eq. (size - 1)*size/2.0_f64) then
         print *, '(REDUCE REAL) PASS'
      else
         stop '(REDUCE REAL) NOT PASS'
      end if
   end if

   SLL_DEALLOCATE_ARRAY(sendbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(recvbuf_real, ierr)

   call sll_s_collective_barrier(sll_v_world_collective)
   if (rank == 0) print *, '-----------------------------'
   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_int(1), ierr)
   sendbuf_int(:) = rank + 1
   SLL_ALLOCATE(recvbuf_int(1), ierr)

   call sll_s_collective_reduce_int(sll_v_world_collective, sendbuf_int, &
                                    1, MPI_SUM, 0, recvbuf_int)

   if (rank == 0) then
      if (recvbuf_int(1) .eq. (1 + size)*size/2.0_f64) then
         print *, '(REDUCE INT) PASS'
      else
         stop '(REDUCE INT) NOT PASS'
      end if
   end if

   SLL_DEALLOCATE_ARRAY(sendbuf_int, ierr)
   SLL_DEALLOCATE_ARRAY(recvbuf_int, ierr)

   call sll_s_collective_barrier(sll_v_world_collective)
   if (rank == 0) print *, '-----------------------------'
   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_log(1), ierr)
   sendbuf_log(:) = .true.
   SLL_ALLOCATE(recvbuf_log(1), ierr)
   !recvbuf_log(:)=.false.

   call sll_s_collective_reduce_logical(sll_v_world_collective, sendbuf_log, &
                                        1, MPI_LAND, 0, recvbuf_log)

   if (rank == 0) then
      if (recvbuf_log(1)) then
         print *, '(REDUCE LOGICAL) PASS'
      else
         stop '(REDUCE LOGICAL) NOT PASS'
      end if
   end if

   SLL_DEALLOCATE_ARRAY(sendbuf_log, ierr)
   SLL_DEALLOCATE_ARRAY(recvbuf_log, ierr)

   call sll_s_collective_barrier(sll_v_world_collective)
   if (rank == 0) print *, '-----------------------------'
   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_log(1), ierr)
   sendbuf_log(:) = .true.
   SLL_ALLOCATE(recvbuf_log(1), ierr)
   recvbuf_log(1) = .false.

   call sll_s_collective_allreduce_logical(sll_v_world_collective, &
                                           sendbuf_log, 1, MPI_LAND, recvbuf_log)

   call sll_s_collective_barrier(sll_v_world_collective)

   call sll_s_collective_reduce_logical(sll_v_world_collective, recvbuf_log, &
                                        1, MPI_LAND, 0, sendbuf_log)

   if (rank == 0) then
      if (sendbuf_log(1)) then
         print *, '(ALLREDUCE LOGICAL) PASS'
      else
         stop '(ALLREDUCE LOGICAL) NOT PASS'
      end if
   end if

   SLL_DEALLOCATE_ARRAY(sendbuf_log, ierr)
   SLL_DEALLOCATE_ARRAY(recvbuf_log, ierr)

   call sll_s_collective_barrier(sll_v_world_collective)
   if (rank == 0) print *, '-----------------------------'
   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_real(1), ierr)
   sendbuf_real(:) = real(rank, f32)
   SLL_ALLOCATE(recvbuf_real(1), ierr)

   call sll_s_collective_allreduce_real32(sll_v_world_collective, &
                                          sendbuf_real, 1, MPI_SUM, recvbuf_real)

   call sll_s_collective_barrier(sll_v_world_collective)

   call sll_s_collective_reduce_real32(sll_v_world_collective, recvbuf_real, &
                                       1, MPI_SUM, 0, sendbuf_real)

   if (rank == 0) then
      if (sendbuf_real(1) - real(size*(size - 1)*size/2.0, kind=f32) < 1e-14) then
         print *, '(ALLREDUCE REAL) PASS'
      else
         stop '(ALLREDUCE REAL) NOT PASS'
      end if
   end if

   SLL_DEALLOCATE_ARRAY(sendbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(recvbuf_real, ierr)

   call sll_s_collective_barrier(sll_v_world_collective)
   if (rank == 0) print *, '-----------------------------'
   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_real(1), ierr)
   SLL_ALLOCATE(somme(1), ierr)
   if (rank == 0) then
      sendbuf_real(1) = 1.0
   end if

   call sll_o_collective_bcast(sll_v_world_collective, sendbuf_real, 1, 0)
   !PRINT *,'(BCAST) ','Me, process ',rank,', I''ve received  ',values,&
   !        ' from process 0'

   call sll_s_collective_reduce_real32(sll_v_world_collective, sendbuf_real, 1, &
                                       MPI_SUM, 0, somme)

   if (rank == 0) then
      if (somme(1) == REAL(size, f32)) then
         print *, '(BCAST REAL) PASS'
      else
         stop '(BCAST REAL) NOT PASS'
      end if
   end if

   call sll_s_collective_barrier(sll_v_world_collective)
   SLL_DEALLOCATE_ARRAY(sendbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(somme, ierr)
   if (rank == 0) print *, '-----------------------------'
   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_real(size), ierr)

   if (rank == 0) then
      ! Intel compilers complain if only process 0 allocates the send_buffer
      !   SLL_ALLOCATE(sendbuf_real(size),ierr)
      sendbuf_real(:) = (/(0.+i, i=0, size - 1)/)
      !PRINT *,'(SCATTER REAL) ', 'Me, process ',rank,'send the values : ',send_buf
   end if

   SLL_ALLOCATE(recvbuf_real(1), ierr)
   SLL_ALLOCATE(somme(1), ierr)

   call sll_o_collective_scatter(sll_v_world_collective, sendbuf_real, 1, &
                                 0, recvbuf_real)
   !PRINT *,'(SCATTER REAL) ', 'Me, process ', rank, ', I''ve received', recvbuf_real, &
   !         ' from process 0'

   call sll_s_collective_reduce_real32(sll_v_world_collective, recvbuf_real, 1, &
                                       MPI_SUM, 0, somme)

   if (rank .eq. 0) then
      if (somme(1) .eq. size*(size - 1)/2.0_f64) then
         print *, '(SCATTER REAL) PASS'
      else
         stop '(SCATTER REAL) NOT PASS'
      end if
   end if

   call sll_s_collective_barrier(sll_v_world_collective)
   SLL_DEALLOCATE_ARRAY(recvbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(somme, ierr)
   if (rank .eq. 0) then
      print *, '-----------------------------'
      SLL_DEALLOCATE_ARRAY(sendbuf_real, ierr)
   end if
   call sll_s_collective_barrier(sll_v_world_collective)

   if (rank .eq. 0) then
      SLL_ALLOCATE(sendbuf_real(size + 1), ierr)
      sendbuf_real(:) = 1.0_f64
   end if

   SLL_ALLOCATE(sendcounts(size), ierr)
   sendcounts(1:size - 1) = 1
   sendcounts(size) = 2

   !print *, 'Me process',rank,', sendcounts ',sendcounts

   SLL_ALLOCATE(sdispls(size), ierr)
   sdispls(1) = 0
   do i = 2, size
      sdispls(i) = sdispls(i - 1) + sendcounts(i - 1)
   end do

   SLL_ALLOCATE(recvcounts(1), ierr)
   recvcounts(1) = sendcounts(rank + 1)

   SLL_ALLOCATE(recvbuf_real(recvcounts(1)), ierr)

   !if(rank==0) print *, 'Me process 0, send ',sendbuf_real

   call sll_s_collective_scatterv_real(sll_v_world_collective, sendbuf_real, &
                                       sendcounts, sdispls, &
                                       recvcounts(1), 0, recvbuf_real)

   !print *, 'Me process ',rank,' I''ve receveid ', recvbuf_real

   SLL_ALLOCATE(somme(1), ierr)
   call sll_s_collective_reduce_real32(sll_v_world_collective, &
                                       (/SUM(recvbuf_real)/), &
                                       1, MPI_SUM, 0, somme)

   if (rank .eq. 0) then
      if (somme(1) .eq. size + 1.0) then
         print *, '(SCATTERV REAL) PASS'
      else
         stop '(SCATTERV REAL) NOT PASS'
      end if
   end if

   call sll_s_collective_barrier(sll_v_world_collective)
   SLL_DEALLOCATE_ARRAY(recvbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(sendcounts, ierr)
   SLL_DEALLOCATE_ARRAY(recvcounts, ierr)
   SLL_DEALLOCATE_ARRAY(somme, ierr)
   SLL_DEALLOCATE_ARRAY(sdispls, ierr)

   if (rank == 0) then
      print *, '-----------------------------'
!   SLL_DEALLOCATE_ARRAY(sendbuf_real,ierr)
   end if
   SLL_DEALLOCATE_ARRAY(sendbuf_real, ierr)
   call sll_s_collective_barrier(sll_v_world_collective)
!___________________________________________________________________

   SLL_ALLOCATE(sendcounts(1), ierr)
   if (rank == 0) then
      SLL_ALLOCATE(sendbuf_real(2), ierr)
      sendbuf_real(:) = real(rank, f32)

      sendcounts(1) = 2

      SLL_ALLOCATE(sdispls(size), ierr)
      sdispls(1) = 1
      do i = 2, size
         sdispls(i) = sdispls(i - 1) + 1
      end do
      sdispls(1) = 0

      SLL_ALLOCATE(recvbuf_real(size + 1), ierr)

      SLL_ALLOCATE(recvcounts(size), ierr)
      recvcounts(1) = 2
      recvcounts(2:size) = 1
   else
      SLL_ALLOCATE(sendbuf_real(1), ierr)
      SLL_ALLOCATE(recvbuf_real(size), ierr)
      SLL_ALLOCATE(recvcounts(size), ierr)
      SLL_ALLOCATE(sdispls(size), ierr)
      sendbuf_real(:) = real(rank, f32)
      sendcounts(1) = 1
   end if

   call sll_s_collective_barrier(sll_v_world_collective)

   call sll_s_collective_gatherv_real(sll_v_world_collective, sendbuf_real, &
                                      sendcounts(1), recvcounts, sdispls, 0, recvbuf_real)

   if (rank == 0) then
      SLL_ALLOCATE(somme(1), ierr)
      somme(1) = SUM(recvbuf_real)
      if (somme(1) == (size - 1)*size/2) then
         print *, '(GATHERV REAL) PASS'
      else
         stop '(GATHERV REAL) NOT PASS'
      end if
      SLL_DEALLOCATE_ARRAY(somme, ierr)
   end if

   call sll_s_collective_barrier(sll_v_world_collective)
   SLL_DEALLOCATE_ARRAY(sendbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(sendcounts, ierr)

   if (rank == 0) then
      print *, '-----------------------------'
   end if

   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_DEALLOCATE_ARRAY(recvbuf_real, ierr)
   SLL_ALLOCATE(sendbuf_real(2), ierr)
   SLL_ALLOCATE(recvbuf_real(size*2), ierr)
   SLL_DEALLOCATE_ARRAY(recvcounts, ierr)
   SLL_DEALLOCATE_ARRAY(sdispls, ierr)
   sendbuf_real(:) = (/rank*2., rank*2.+1.0/)
   !PRINT *,'(GATHER) ', 'Me, process ', rank, 'send the values : ', sendbuf_real, &
   !         'to the process 0'

   call sll_o_collective_gather(sll_v_world_collective, sendbuf_real, 2, 0, &
                                recvbuf_real)

   SLL_ALLOCATE(somme(1), ierr)
   call sll_s_collective_reduce_real32(sll_v_world_collective, &
                                       (/SUM(sendbuf_real)/), 1, MPI_SUM, 0, somme)

   !IF(rank==0) THEN
   ! PRINT *,'(GATHER) ', 'Me, process 0.', ' I''ve receveid the values : ', recvbuf_real
   !ENDIF

   if (rank == 0) then
      if (somme(1) == size*(2.0*size - 1)) then
         print *, '(GATHER REAL) PASS'
      else
         stop '(GATHER REAL) NOT PASS'
      end if
   end if

   call sll_s_collective_barrier(sll_v_world_collective)
   SLL_DEALLOCATE_ARRAY(sendbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(recvbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(somme, ierr)
   if (rank == 0) print *, '-----------------------------'
   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_int(2), ierr)
   SLL_ALLOCATE(recvbuf_int(size*2), ierr)
   sendbuf_int(:) = (/rank*2, rank*2 + 1/)
   !PRINT *,'(ALLGATHER) ', 'Me, process ', rank, 'send the values : ',&
   !         sendbuf_int,'to all process'

   call sll_o_collective_allgather(sll_v_world_collective, sendbuf_int, 2, &
                                   recvbuf_int, 2)
   SLL_ALLOCATE(somme(1), ierr)
   somme(1) = SUM(real(recvbuf_int, f32))
   if (somme(1) == size*(2.0*size - 1)) then
      logic(1) = .true.
   else
      logic(1) = .false.
   end if

   call sll_s_collective_reduce_logical(sll_v_world_collective, logic, &
                                        1, MPI_LAND, 0, logic2)

   !PRINT *,'(ALLGATHER) Me, process ', rank, &
   !          ' I''ve receveid the values : ', recvbuf_int

   if (rank == 0) then
      if (logic2(1)) then
         print *, '(ALLGATHER REAL) PASS'
      else
         stop '(ALLGATHER REAL) NOT PASS'
      end if
   end if

   call sll_s_collective_barrier(sll_v_world_collective)
   if (rank == 0) print *, '-----------------------------'
   SLL_DEALLOCATE_ARRAY(recvbuf_int, ierr)
   SLL_DEALLOCATE_ARRAY(sendbuf_int, ierr)
   SLL_DEALLOCATE_ARRAY(somme, ierr)
   call sll_s_collective_barrier(sll_v_world_collective)

   if (rank .eq. 0) then
      SLL_ALLOCATE(sendbuf_real(2), ierr)
      sendbuf_real(1) = 0.0_f32
      sendbuf_real(2) = 1.0_f32
   else
      SLL_ALLOCATE(sendbuf_real(1), ierr)
      sendbuf_real(1) = real(rank + 1, f32)
   end if
   SLL_ALLOCATE(sendcounts(1), ierr)
   if (rank .eq. 0) then
      sendcounts(1) = 2
   else
      sendcounts(1) = 1
   end if
   SLL_ALLOCATE(recvcounts(size), ierr)
   recvcounts(1) = 2
   recvcounts(2:size) = 1
   SLL_ALLOCATE(sdispls(size), ierr)
   sdispls(1) = 0
   if (size >= 2) then
      do i = 2, size
         sdispls(i) = sdispls(i - 1) + recvcounts(i - 1)
      end do
   end if
   SLL_ALLOCATE(recvbuf_real(SUM(recvcounts)), ierr)

   call sll_o_collective_allgatherv( &
      sll_v_world_collective, &
      sendbuf_real, &
      sendcounts(1), &
      recvcounts, &
      sdispls, &
      recvbuf_real)

   if (SUM(recvbuf_real) .eq. (size)*(size + 1)/2.0) then
      logic(1) = .true.
   else
      logic(1) = .false.
   end if

   call sll_s_collective_reduce_logical(sll_v_world_collective, logic, &
                                        1, MPI_LAND, 0, logic2)

   if (rank == 0) then
      if (logic2(1)) then
         print *, '(ALLGATHERV REAL) PASS'
      else
         stop '(ALLGATHERV REAL) NOT PASS'
      end if
   end if

   call sll_s_collective_barrier(sll_v_world_collective)
   SLL_DEALLOCATE_ARRAY(sendbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(recvbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(sendcounts, ierr)
   SLL_DEALLOCATE_ARRAY(recvcounts, ierr)
   SLL_DEALLOCATE_ARRAY(sdispls, ierr)
   if (rank == 0) print *, '-----------------------------'
   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_int(size), ierr)
   sendbuf_int(:) = (/(0 + i, i=0, size - 1)/)
   SLL_ALLOCATE(recvbuf_int(size), ierr)

   !PRINT *,'Moi, processus ',rank,'envoie mon tableau valeurs : ',&
   !                                 sendbuf_int(:)

   call sll_s_collective_alltoall_int(sendbuf_int, 1, 1, &
                                      recvbuf_int, sll_v_world_collective)

   !PRINT *,'Moi, processus ',rank,', j''ai recu ',recvbuf_int

   SLL_ALLOCATE(somme(1), ierr)
   somme(1) = REAL(SUM(recvbuf_int))
   if (somme(1) .eq. size*rank*1.0) then
      logic(1) = .true.
   else
      logic(1) = .false.
   end if

   call sll_s_collective_reduce_logical(sll_v_world_collective, logic, &
                                        1, MPI_LAND, 0, logic2)

   if (rank == 0) then
      if (logic2(1)) then
         print *, '(ALLTOALL INT) PASS'
      else
         stop '(ALLTOALL INT) NOT PASS'
      end if
   end if

   call sll_s_collective_barrier(sll_v_world_collective)
   if (rank == 0) print *, '-----------------------------'
   SLL_DEALLOCATE_ARRAY(recvbuf_int, ierr)
   SLL_DEALLOCATE_ARRAY(sendbuf_int, ierr)
   SLL_DEALLOCATE_ARRAY(somme, ierr)
   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_int(size + 1), ierr)
   sendbuf_int(:) = rank

   SLL_ALLOCATE(sendcounts(size), ierr)

   if (MOD(rank, 2) == 0) then
      sendcounts(1:size - 1) = 1
      sendcounts(size) = 2
   else
      sendcounts(2:size) = 1
      sendcounts(1) = 2
   end if

   call sll_s_collective_alltoallv_int_simple(sendbuf_int, sendcounts, &
                                              recvbuf_int, sll_v_world_collective)

   SLL_ALLOCATE(somme(1), ierr)
   call sll_s_collective_reduce_int(sll_v_world_collective, &
                                    (/SUM(recvbuf_int)/), 1, MPI_SUM, 0, sendbuf_int)
   if (rank .eq. 0) then
      if (sendbuf_int(1) .eq. (size + 1)*(size - 1)*size/2) then
         print *, '(ALLTOALLV INT) PASS'
      else
         stop '(ALLTOALLV INT) NOT PASS'
      end if
   end if

   call sll_s_collective_barrier(sll_v_world_collective)
   if (rank == 0) print *, '-----------------------------'
   SLL_DEALLOCATE_ARRAY(recvbuf_int, ierr)
   SLL_DEALLOCATE_ARRAY(sendbuf_int, ierr)
   SLL_DEALLOCATE_ARRAY(somme, ierr)
   SLL_DEALLOCATE_ARRAY(sendcounts, ierr)
   call sll_s_collective_barrier(sll_v_world_collective)

   SLL_ALLOCATE(sendbuf_real(size + 1), ierr)
   sendbuf_real(:) = real(rank, f32)

   SLL_ALLOCATE(sendcounts(size), ierr)

   if (MOD(rank, 2) == 0) then
      sendcounts(1:size - 1) = 1
      sendcounts(size) = 2
   else
      sendcounts(2:size) = 1
      sendcounts(1) = 2
   end if

   ! Define RECV_CNTS
   SLL_ALLOCATE(recvcounts(size), ierr)
   call sll_s_collective_alltoall_int(sendcounts, 1, 1, &
                                      recvcounts, sll_v_world_collective)

   ! Define RECV_BUF
   SLL_ALLOCATE(recvbuf_real(SUM(recvcounts)), ierr)

   ! Define SEND_DISPLS
   SLL_ALLOCATE(sdispls(size), ierr)
   sdispls(1) = 0
   do i = 2, size
      sdispls(i) = sdispls(i - 1) + sendcounts(i - 1)
   end do

   ! Define RECV_DISPLS
   SLL_ALLOCATE(rdispls(size), ierr)
   rdispls(1) = 0
   do i = 2, size
      rdispls(i) = rdispls(i - 1) + recvcounts(i - 1)
   end do

   call sll_s_collective_alltoallv_real(sendbuf_real, sendcounts, &
                                        sdispls, &
                                        recvbuf_real, recvcounts, &
                                        rdispls, sll_v_world_collective)

   SLL_ALLOCATE(somme(1), ierr)
   call sll_s_collective_reduce_real32(sll_v_world_collective, &
                                       (/SUM(recvbuf_real(:))/), 1, MPI_SUM, 0, somme)
   if (rank .eq. 0) then
      if (somme(1) .eq. (size + 1)*(size - 1)*size/2.0) then
         print *, '(ALLTOALLV REAL) PASS'
      else
         stop '(ALLTOALLV REAL) NOT PASS'
      end if
   end if

   SLL_DEALLOCATE_ARRAY(somme, ierr)
   SLL_DEALLOCATE_ARRAY(rdispls, ierr)
   SLL_DEALLOCATE_ARRAY(sdispls, ierr)
   SLL_DEALLOCATE_ARRAY(sendbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(recvbuf_real, ierr)
   SLL_DEALLOCATE_ARRAY(sendcounts, ierr)
   SLL_DEALLOCATE_ARRAY(recvcounts, ierr)

   call test_sll_collective_globalsum()

   call sll_s_halt_collective()

contains

!> @brief tests sll_o_collective_globalsum interface
!> If this test passes, the reduce and allreduce are also working
   subroutine test_sll_collective_globalsum
      sll_real64 :: summand_real64
      sll_real32 :: summand_real32
      sll_comp64 :: summand_comp64
      sll_comp32 :: summand_comp32
      sll_int32 ::  summand_int32

      if (rank == 0) print *, '-----------------------------'
      call sll_s_collective_barrier(sll_v_world_collective)

      summand_real64 = 1.0_f64 !or better sll_p_pi
      call sll_o_collective_globalsum(sll_v_world_collective, summand_real64, 0)

      summand_comp64 = (1.0_f64, 0.0_f64)
      call sll_o_collective_globalsum(sll_v_world_collective, summand_comp64, 0)

      summand_comp32 = (1.0_f32, 0.0_f32)
      call sll_o_collective_globalsum(sll_v_world_collective, summand_comp32, 0)

      summand_real32 = 1.0_f32
      call sll_o_collective_globalsum(sll_v_world_collective, summand_real32, 0)

      summand_int32 = 1
      call sll_o_collective_globalsum(sll_v_world_collective, summand_int32, 0)

      if (rank == 0) then

         if (summand_real64 == 1.0_f64*size .and. &
             summand_real32 == 1.0_f64*size .and. &
             summand_comp32 == 1.0_f64*size .and. &
             summand_comp64 == 1.0_f64*size .and. &
             summand_int32 == size) then
            print *, '(GLOBALSUM) PASS'
         else
            !print *, summand_real32, summand_real64, summand_comp32 , summand_comp64 , summand_int32
            stop '(GLOBALSUM) NOT PASS'
         end if

      end if

   end subroutine test_sll_collective_globalsum

end program collective_test
