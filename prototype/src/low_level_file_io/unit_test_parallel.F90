program test_io
  use sll_collective
#include "sll_remap.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
  implicit none

  ! Test of the 2D remapper takes a 2D array whose global size Nx*Ny,
  ! distributed among NPi*NPj processors.
  sll_real64, dimension(:,:), allocatable :: local_array
  sll_real64, dimension(:,:), allocatable ::  arrays_diff
  ! Take a 2D array of dimensions ni*nj where ni, nj are the dimensions of
  ! the full array.
  integer , parameter                       :: ni = 512
  integer , parameter                       :: nj = 128
  ! Local sizes
  integer                                   :: loc_sz_i_init
  integer                                   :: loc_sz_j_init
  integer                                   :: loc_sz_i_final
  integer                                   :: loc_sz_j_final

  ! the process mesh
  integer                                   :: npi
  integer                                   :: npj
  sll_int32                                 :: gi, gj
  integer                                   :: ierr
  integer                                   :: myrank
  sll_int64                                 :: colsz        ! collective size
  ! Remap stuff
  type(layout_2D), pointer                  :: layout
  type(remap_plan_2D), pointer              :: rmp2

  sll_real64                                :: rand_real
  integer, parameter                        :: nbtest = 25
  integer                                   :: i_test
  integer                                   :: i, j
  sll_int32, dimension(2)                   :: global_indices, g
  sll_real32   , dimension(1)               :: prod4test
  logical                                   :: test_passed
  integer                                   :: ok

  test_passed = .true.

  ! Boot parallel environment
  call sll_boot_collective()

  colsz  = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  if( myrank .eq. 0) then
     print *, ' '
     print *, '--------------- REMAP test ---------------------'
     print *, ' '
     print *, 'Running a test on ', colsz, 'processes'
     call flush()
  end if

  if (.not. is_power_of_two(colsz)) then     
     print *, 'This test needs to run in a number of processes which is ',&
          'a power of 2.'
     stop
  end if

  layout => new_layout_2D( sll_world_collective )        
  call factorize_in_random_2powers_2d(colsz, npi, npj)

  if( myrank .eq. 0 ) then
     print *, 'source configuration: ', npi, npj
  end if

  call initialize_layout_with_distributed_2D_array( &
          ni, &
          nj, &
          npi, &
          npj, &
          layout )
     
  call compute_local_sizes_2d( layout, loc_sz_i_init, loc_sz_j_init)        

  SLL_ALLOCATE(local_array(loc_sz_i_init,loc_sz_j_init),ierr)
 
  ! initialize the local data    
  do j=1,loc_sz_j_init 
     do i=1,loc_sz_i_init
        global_indices =  local_to_global_2D( layout, (/i, j/) )
        gi = global_indices(1)
        gj = global_indices(2)
        local_array(i,j) = gi + (gj-1)*ni
     enddo
  enddo
     
  if( myrank .eq. 0) then
     print *, ' '
     print *, '-------------------------------------------'
     print *, ' '
     call flush()
  end if
  call flush() 
       
  call sll_collective_barrier(sll_world_collective)
  
  call delete_layout_2D( layout )
  SLL_DEALLOCATE_ARRAY(local_array, ierr)

  call sll_halt_collective()
  
contains

  subroutine factorize_in_random_2powers_2d(n, n1, n2)
    sll_int64, intent(in) :: n
    integer, intent(out)  ::n1, n2
    integer   :: expo, expo1, expo2
    if (.not.is_power_of_two(n)) then   
       print*, 'The number of processors must be a power of 2'
       stop
    endif 
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    expo2 = expo - expo1
    n1 = 2**expo1
    n2 = 2**expo2
  end subroutine factorize_in_random_2powers_2d

  subroutine compute_local_sizes_2d( layout, loc_sz_i, loc_sz_j )
    type(layout_2D), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to new_distributed_mesh_2D'
       print *, 'Exiting...'
       STOP
    end if
    my_rank = sll_get_collective_rank(get_layout_2D_collective(layout))
    i_min = get_layout_2D_i_min( layout, my_rank )
    i_max = get_layout_2D_i_max( layout, my_rank )
    j_min = get_layout_2D_j_min( layout, my_rank )
    j_max = get_layout_2D_j_max( layout, my_rank )
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
  end subroutine compute_local_sizes_2d
      
end program test_io
