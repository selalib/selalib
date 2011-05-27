program collective_test
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_collective
  implicit none
  intrinsic :: int

  ! test on a mesh of size NPXxNPYxNPZ
#define NPX 2
#define NPY 3
#define NPZ 4

  type(sll_collective_t), pointer :: col
  sll_int32               :: IC, JC, KC  ! node coordinates
  sll_int32               :: rank
  sll_int32               :: col_size
  sll_int32               :: ierr
  sll_int32, dimension(:), allocatable :: local_array
  sll_int32, dimension(:), allocatable :: global_array


  print *, ' '
  print *, '************* Unit test for the collective module ****************'
  print *, ' '

  call sll_boot_collective()

  ! Figure out what are this node's coordinates in the 'process mesh':
  ! IC in [0,NPX-1]
  ! JC in [0,NPY-1]
  ! KC in [0,NPZ-1]
  rank     = sll_get_collective_rank( sll_world_collective )
  col_size = sll_get_collective_size( sll_world_collective )
  KC       = int(rank/(NPX*NPY))
  JC       = int((rank - KC*NPX*NPY)/NPX)
  IC       = rank - KC*NPX*NPY - JC*NPX
  write (*,'(a i8 a i4 i4 i4 a)') 'rank: ', rank, &
       '. Has coords [', IC, JC, KC,']'
  call flush(6)

  ! Test allgather
  SLL_ALLOCATE( local_array(2), ierr )
  SLL_ALLOCATE( global_array(0:(2*col_size-1)), ierr )
  local_array(:) = (/rank, rank/)
  print *, 'initialized local array'
  print *, local_array(:)
  call flush()
  call sll_collective_allgather( sll_world_collective, local_array(:), 2, global_array(:), 2 )
  print *, global_array(:)
  call flush()


  ! Test splitting the collective into rods along 'Z', a direction in which
  ! we have 4 domains.
  col => sll_new_collective( sll_world_collective, (IC+NPX*JC), rank )

  write (*, '(a i4 a i4 a i4 a i4)') 'Process: ', rank, ': color = ', sll_get_collective_color(col), ', rank in new collective: ', sll_get_collective_rank(col), '. Size = ', sll_get_collective_size( col )
  call sll_delete_collective(col)
  write (*, '(a)') 'Deleted new collective, ready to exit...'
  call sll_halt_collective()

  print *, 'Test complete'



end program collective_test
