program unit_test_buffer_loader_utilities
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_collective
  use sll_remapper
  use sll_buffer_loader_utilities_module

implicit none

  sll_int32, dimension(:), allocatable :: collective_displs
  sll_int32, dimension(:), allocatable :: collective_recvcnts
  sll_int32 :: collective_size
  sll_real32,dimension(:,:),pointer :: f_visu 
  sll_real32,dimension(:),pointer :: f_visu_buf1d
  sll_real32,dimension(:),pointer :: f_x1_buf1d
  sll_real32,dimension(:,:),pointer :: f_x1
  sll_int32 :: np_x1
  sll_int32 :: np_x2
  sll_int32 :: ierr
  type(layout_2D), pointer :: layout_x1
  sll_int32 :: nproc_x1
  sll_int32 :: nproc_x2
  sll_int32 :: local_size_x1
  sll_int32 :: local_size_x2
  sll_int32 :: i


  
  np_x1 = 2
  np_x2 = 2

  ! Boot parallel environment
  call sll_boot_collective()


  if(sll_get_collective_rank(sll_world_collective)==0)then
    SLL_ALLOCATE(f_visu(np_x1,np_x2),ierr)
    SLL_ALLOCATE(f_visu_buf1d(np_x1*np_x2),ierr)
  else
      SLL_ALLOCATE(f_visu(1:0,1:0),ierr)          
      SLL_ALLOCATE(f_visu_buf1d(1:0),ierr)          
  endif

  collective_size = sll_get_collective_size(sll_world_collective)
  SLL_ALLOCATE(collective_displs(collective_size),ierr)
  SLL_ALLOCATE(collective_recvcnts(collective_size),ierr)

  layout_x1       => new_layout_2D( sll_world_collective )
  nproc_x1 = sll_get_collective_size( sll_world_collective )
  nproc_x2 = 1
  call initialize_layout_with_distributed_2D_array( &
    np_x1, np_x2, nproc_x2, nproc_x1, layout_x1 )


    call compute_displacements_array_2d( &
      layout_x1, &
      collective_size, &
      collective_displs )
    collective_recvcnts = receive_counts_array_2d( &
      layout_x1, &
      collective_size )

    call compute_local_sizes_2d( layout_x1, local_size_x1, local_size_x2 )
    SLL_ALLOCATE(f_x1(local_size_x1,local_size_x2),ierr)    
    SLL_ALLOCATE(f_x1_buf1d(local_size_x1*local_size_x2),ierr)    
    
    
    f_x1 = 1._f64

    print *,'#collective_recvcnts',collective_recvcnts
    
    print *,'#collective_displs',collective_displs

    call load_buffer32_2d( layout_x1, f_x1, f_x1_buf1d )
    !call sll_collective_gatherv_real64( &
    print *,'#before',sll_get_collective_rank(sll_world_collective), maxval(f_x1_buf1d),minval(f_x1_buf1d)
      do i=1,size(f_x1_buf1d,1)
       print *,'#',sll_get_collective_rank(sll_world_collective),i,f_x1_buf1d(i)
      enddo    


    !print *,'#sendcount=', 

    call sll_collective_barrier( sll_world_collective )
    
    
    
    
    call sll_collective_gatherv_real( &
      sll_world_collective, &
      f_x1_buf1d, &
      size(f_x1_buf1d,1), &
      !sll_get_collective_rank(sll_world_collective), &
      collective_recvcnts, &
      collective_displs, &
      0, &
      f_visu_buf1d )
    
    call sll_collective_barrier( sll_world_collective )
    
    if(sll_get_collective_rank(sll_world_collective)==0)then
      print *,'#after',sll_get_collective_rank(sll_world_collective), maxval(f_visu_buf1d),minval(f_visu_buf1d)
      do i=1,size(f_visu_buf1d,1)
       print *,'#',i,f_visu_buf1d(i)
      enddo    
    endif
    f_visu = reshape(f_visu_buf1d, shape(f_visu))


  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *,'#PASSED'
  endif
  call sll_halt_collective()

end program