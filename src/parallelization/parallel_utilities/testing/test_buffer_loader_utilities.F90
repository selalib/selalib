program unit_test_buffer_loader_utilities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_buffer_loader_utilities, only: &
    sll_s_compute_displacements_array_2d, &
    sll_s_load_buffer32_2d, &
    sll_f_receive_counts_array_2d

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_s_collective_barrier, &
    sll_s_collective_gatherv_real, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_remapper, only: &
    sll_o_compute_local_sizes, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_f_new_layout_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  type(sll_t_layout_2d), pointer :: layout_x1
  sll_int32 :: nproc_x1
  sll_int32 :: nproc_x2
  sll_int32 :: local_size_x1
  sll_int32 :: local_size_x2
  sll_int32 :: i


  
  np_x1 = 2
  np_x2 = 2

  ! Boot parallel environment
  call sll_s_boot_collective()


  if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
    SLL_ALLOCATE(f_visu(np_x1,np_x2),ierr)
    SLL_ALLOCATE(f_visu_buf1d(np_x1*np_x2),ierr)
  else
      SLL_ALLOCATE(f_visu(1:0,1:0),ierr)          
      SLL_ALLOCATE(f_visu_buf1d(1:0),ierr)          
  endif

  collective_size = sll_f_get_collective_size(sll_v_world_collective)
  SLL_ALLOCATE(collective_displs(collective_size),ierr)
  SLL_ALLOCATE(collective_recvcnts(collective_size),ierr)

  layout_x1       => sll_f_new_layout_2d( sll_v_world_collective )
  nproc_x1 = sll_f_get_collective_size( sll_v_world_collective )
  nproc_x2 = 1
  call sll_o_initialize_layout_with_distributed_array( &
    np_x1, np_x2, nproc_x2, nproc_x1, layout_x1 )


    call sll_s_compute_displacements_array_2d( &
      layout_x1, &
      collective_size, &
      collective_displs )
    collective_recvcnts = sll_f_receive_counts_array_2d( &
      layout_x1, &
      collective_size )

    call sll_o_compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
    SLL_ALLOCATE(f_x1(local_size_x1,local_size_x2),ierr)    
    SLL_ALLOCATE(f_x1_buf1d(local_size_x1*local_size_x2),ierr)    
    
    
    f_x1 = 1._f64

    print *,'#collective_recvcnts',collective_recvcnts
    
    print *,'#collective_displs',collective_displs

    call sll_s_load_buffer32_2d( layout_x1, f_x1, f_x1_buf1d )
    !call sll_s_collective_gatherv_real64( &
    print *,'#before',sll_f_get_collective_rank(sll_v_world_collective), maxval(f_x1_buf1d),minval(f_x1_buf1d)
      do i=1,size(f_x1_buf1d,1)
       print *,'#',sll_f_get_collective_rank(sll_v_world_collective),i,f_x1_buf1d(i)
      enddo    


    !print *,'#sendcount=', 

    call sll_s_collective_barrier( sll_v_world_collective )
    
    
    
    
    call sll_s_collective_gatherv_real( &
      sll_v_world_collective, &
      f_x1_buf1d, &
      size(f_x1_buf1d,1), &
      !sll_f_get_collective_rank(sll_v_world_collective), &
      collective_recvcnts, &
      collective_displs, &
      0, &
      f_visu_buf1d )
    
    call sll_s_collective_barrier( sll_v_world_collective )
    
    if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
      print *,'#after',sll_f_get_collective_rank(sll_v_world_collective), maxval(f_visu_buf1d),minval(f_visu_buf1d)
      do i=1,size(f_visu_buf1d,1)
       print *,'#',i,f_visu_buf1d(i)
      enddo    
    endif
    f_visu = reshape(f_visu_buf1d, shape(f_visu))


  if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
    print *,'#PASSED'
  endif
  call sll_s_halt_collective()

end program