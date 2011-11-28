program collective_test
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_collective
  implicit none
  intrinsic :: int

  sll_int32 :: rank, size, i
  sll_real32, ALLOCATABLE, DIMENSION(:) :: values
  sll_int32, ALLOCATABLE, DIMENSION(:) :: valuesi
  sll_real32, ALLOCATABLE, DIMENSION(:) :: recbuf
  sll_int32, ALLOCATABLE, DIMENSION(:) :: recbufi
  !sll_real32, DIMENSION(1) :: sum
  LOGICAL, DIMENSION(1) :: logic,logic2
  sll_real32, ALLOCATABLE, DIMENSION(:) :: somme
  !sll_int32, ALLOCATABLE, DIMENSION(:) :: sommei
  sll_real32, DIMENSION(2) :: summ
  sll_real32 :: siz

  CALL sll_boot_collective()
  rank = sll_get_collective_rank( sll_world_collective )
  size = sll_get_collective_size( sll_world_collective )
  ALLOCATE(values(1))
  ALLOCATE(somme(1))
  IF(rank == 0) THEN
   values(1)=1.0
  ENDIF
  CALL sll_collective_bcast( sll_world_collective, values, 1, 0)
  !PRINT *,'(BCAST) ','Me, process ',rank,', I''ve received  ',values,&
  !        ' from process 0'
  
  CALL sll_collective_reduce_real(sll_world_collective, values, 1, &
                                  MPI_SUM,0,somme)

  IF( rank == 0 ) THEN
   IF( somme(1) == REAL(size,f32)) THEN
    PRINT *,'(BCAST) ', 'PASSED'
   ELSE
    stop 'NOT PASSED'
   ENDIF
  ENDIF

  CALL sll_collective_barrier(sll_world_collective)
  DEALLOCATE(values)
  DEALLOCATE(somme)
  IF(rank==0) PRINT *,'-----------------------------'  
  CALL sll_collective_barrier(sll_world_collective)

  ALLOCATE(recbuf(1))
  ALLOCATE(somme(1))
  IF(rank==0) THEN
   ALLOCATE(values(size))
   values(:)=(/(0.+i,i=0,size-1)/)
   !PRINT *,'(SCATTER) ', 'Me, process ',rank,'send the values : ',values
  ENDIF
  
  CALL sll_collective_scatter( sll_world_collective, values, 1, &
                                 0,  recbuf)
  !PRINT *,'(SCATTER) ', 'Me, process ', rank, ', I''ve received', recbuf, &
  !         ' from process 0'

  CALL sll_collective_reduce_real(sll_world_collective, recbuf, 1, &
                                  MPI_SUM,0,somme)

  IF( rank == 0 ) THEN
   siz = REAL(size,f32)
   IF( somme(1) == siz*(siz-1)/2) THEN
    PRINT *,'(SCATTER) ', 'PASSED'
   ELSE
    stop '(SCATTER) NOT PASSED'
   ENDIF
  ENDIF

  CALL sll_collective_barrier(sll_world_collective)
  DEALLOCATE(recbuf)
  DEALLOCATE(somme)
  IF(rank==0) THEN
   PRINT *,'-----------------------------'  
   DEALLOCATE(values)
  ENDIF
  CALL sll_collective_barrier(sll_world_collective)

  ALLOCATE(values(2))
  ALLOCATE(recbuf(size*2))
  values(:)=(/ rank*2. , rank*2. + 1.0 /)
  !PRINT *,'(GATHER) ', 'Me, process ', rank, 'send the values : ', values, &
  !         'to the process 0'
  
  CALL sll_collective_gather( sll_world_collective, values, 2, 0, &
                                  recbuf )

  CALL sll_collective_reduce_real(sll_world_collective,values,2,&
                                  MPI_SUM,0,summ)

  !IF(rank==0) THEN
  ! PRINT *,'(GATHER) ', 'Me, process 0.', ' I''ve receveid the values : ', recbuf
  !ENDIF

  IF( rank == 0 ) THEN
   siz = REAL(size,f32)
   IF( summ(1)+summ(2) == siz*(2*siz-1)) THEN
    PRINT *,'(GATHER) ', 'PASSED'
   ELSE
    stop '(GATHER) NOT PASSED'
   ENDIF
  ENDIF

  CALL sll_collective_barrier(sll_world_collective)
  IF(rank==0) PRINT *,'-----------------------------'  
  DEALLOCATE(recbuf)
  CALL sll_collective_barrier(sll_world_collective)

  ALLOCATE(valuesi(2))
  ALLOCATE(recbufi(size*2))
  valuesi(:)=(/ rank*2 , rank*2 + 1 /)
  !PRINT *,'(ALLGATHER) ', 'Me, process ', rank, 'send the values : ',&
  !         values,'to all process'
  
  CALL sll_collective_allgather( sll_world_collective, valuesi(:), 2, &
                                  recbufi(:), 2 )
  ALLOCATE(somme(1))  
  somme(1) = SUM(recbufi(:))
  IF( somme(1)==size*(2*size-1) ) THEN
   logic(1)=.TRUE.
  ELSE
   logic(1)=.FALSE.
  ENDIF

  CALL sll_collective_reduce_logical(sll_world_collective,logic,&
                                      1,MPI_LAND,0,logic2)

  !PRINT *,'(ALLGATHER) Me, process ', rank, &
  !          ' I''ve receveid the values : ', recbufi

  IF( rank == 0 ) THEN
   IF( logic2(1) ) THEN
    PRINT *,'(ALLGATHER) ', 'PASSED'
   ELSE
    stop '(ALLGATHER) NOT PASSED'
   ENDIF
  ENDIF
  
  CALL sll_halt_collective()
end program collective_test
