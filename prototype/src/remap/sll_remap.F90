module remapper
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "misc_utils.h"
  use sll_collective
  implicit none
  
  ! The box types contain information on the index limits contained        
  ! in a given processor.
  type box_3D
     sll_int32 :: i_min, i_max
     sll_int32 :: j_min, j_max
     sll_int32 :: k_min, k_max
  end type box_3D

  type box_4D
     sll_int32 :: i_min, i_max
     sll_int32 :: j_min, j_max
     sll_int32 :: k_min, k_max
     sll_int32 :: l_min, l_max
  end type box_4D

  type box_5D
     sll_int32 :: i_min, i_max
     sll_int32 :: j_min, j_max
     sll_int32 :: k_min, k_max
     sll_int32 :: l_min, l_max
     sll_int32 :: m_min, m_max
  end type box_5D
  
  ! The sll_layout types contain information on a collective and an
  ! array of boxes that describes the distribution of data among
  ! different nodes.
  type layout_3D_t
     type(sll_collective_t), pointer     :: collective
     type(box_3D), dimension(:), pointer :: boxes
  end type layout_3D_t

  type layout_4D_t
     type(sll_collective_t), pointer     :: collective
     type(box_4D), dimension(:), pointer :: boxes
  end type layout_4D_t

  type layout_5D_t
     type(sll_collective_t), pointer     :: collective
     type(box_5D), dimension(:), pointer :: boxes
  end type layout_5D_t

  ! Since the plan stores the information on box intersections, now
  ! we need a different type of plan for every dimension. It is also
  ! wasteful to allocate a full array of size(collective_size) to store
  ! the intersection information. Eventually this has to be made leaner
  ! by storing only pointers to boxes and possibly a linked list.
  !
  ! Coalesce with a macro. Note that regardless of the dimensionality of the
  ! data, the arrays are always linear. This implies manual packing/unpacking.
  !
  ! The 'is_uniform' slot is a logical flag that indicates whether the same
  ! amount of data is going to be sent to all the processes in a 
  ! communicator (sender included). This is important to know because we can
  ! then replace the call to alltoallv by a call to alltoall.
  type remap_plan_3D_t
     type(layout_3D_t), pointer          :: initial_layout
     type(layout_3D_t), pointer          :: final_layout
     integer, dimension(:), pointer      :: send_displs
     integer, dimension(:), pointer      :: send_counts
     integer, dimension(:), pointer      :: recv_displs
     integer, dimension(:), pointer      :: recv_counts
     type(box_3D), dimension(:), pointer :: send_boxes
     type(box_3D), dimension(:), pointer :: recv_boxes
     type(sll_collective_t), pointer     :: collective
     sll_int32, dimension(:), pointer    :: send_buffer
     sll_int32, dimension(:), pointer    :: recv_buffer
     logical                             :: is_uniform
  end type remap_plan_3D_t

 interface get_layout_i_min
     module procedure get_layout_3D_i_min
  end interface

  interface set_layout_i_min
     module procedure set_layout_3D_i_min
  end interface

  interface get_layout_i_max
     module procedure get_layout_3D_i_max
  end interface

  interface set_layout_i_max
     module procedure set_layout_3D_i_max
  end interface

  interface get_layout_j_min
     module procedure get_layout_3D_j_min
  end interface

  interface set_layout_j_min
     module procedure set_layout_3D_j_min
  end interface

  interface get_layout_j_max
     module procedure get_layout_3D_j_max
  end interface

  interface set_layout_j_max
     module procedure set_layout_3D_j_max
  end interface

 interface get_layout_k_min
    module procedure get_layout_3D_k_min
  end interface

  interface set_layout_k_min
     module procedure set_layout_3D_k_min
  end interface

  interface get_layout_k_max
     module procedure get_layout_3D_k_min
  end interface

  interface set_layout_k_max
     module procedure set_layout_3D_k_max
  end interface

  interface get_layout_num_nodes
     module procedure get_layout_3D_num_nodes
  end interface

  interface get_layout_box
     module procedure get_layout_3D_box
  end interface

  interface sll_get_num_nodes
     module procedure sll_get_num_nodes_3D
  end interface

  interface apply_remap_3D
     module procedure apply_remap_3D_int, apply_remap_3D_double
  end interface

contains  !******************************************************************



  ! EXPERIMENTAL (in the Fortran context, that is)
  ! Here we explore how to address some code redundancies with a macro.
  ! This may come in handy when we can parametrize the creation of functions
  ! and thus centralize the code generation. The macro must be essentially
  ! bug-free lest it becomes a headache. In case that this is considered
  ! an abuse of the macro facility, we can convert into the multiple
  ! explicit declarations. But the centralized nature of the macro should
  ! not be discounted easily.
  ! 
  ! One thing learned from this is that the macros called from inside a
  ! macro like this (like SLL_ALLOCATE), don't need the semicolon, as they
  ! already have one themselves... it might be a good idea to remove the
  ! semicolon from the last line of all macros, so that they don't introduce
  ! this type of inconsistencies...

#define NEW_LAYOUT_FUNCTION( func_name, layout_type )           \
  function func_name( col );                                    \
    intrinsic :: associated;                                    \
    type(layout_type), pointer          :: func_name;           \
    type(sll_collective_t),  pointer    :: col;                 \
    sll_int32                           :: n_nodes;             \
    sll_int32                           :: ierr;                \
    if( .not. associated(col) ) then;                           \
       write (*,'(a)') 'ERROR: uninitialized collective';       \
       STOP 'NEW_LAYOUT_FUNCTION';                              \
    end if;                                                     \
    SLL_ALLOCATE( func_name, ierr )                             \
    func_name%collective => col;                                \
    n_nodes              = sll_get_collective_size(col);        \
    SLL_ALLOCATE( func_name%boxes(0:(n_nodes-1)), ierr )        \
  end function func_name

NEW_LAYOUT_FUNCTION( new_layout_3D, layout_3D_t )
NEW_LAYOUT_FUNCTION( new_layout_4D, layout_4D_t )
NEW_LAYOUT_FUNCTION( new_layout_5D, layout_5D_t )

#define NEW_DELETE_LAYOUT_FUNCTION( fname, layout_type )        \
  subroutine fname( layout );                                   \
    type(layout_type), pointer :: layout;                       \
    sll_int32                  :: ierr;                         \
    nullify( layout%collective );                               \
    SLL_DEALLOCATE( layout%boxes, ierr );                       \
    SLL_DEALLOCATE( layout, ierr );                             \
  end subroutine fname

NEW_DELETE_LAYOUT_FUNCTION( delete_layout_3D, layout_3D_t )
NEW_DELETE_LAYOUT_FUNCTION( delete_layout_4D, layout_4D_t )
NEW_DELETE_LAYOUT_FUNCTION( delete_layout_5D, layout_5D_t )

  ! Access functions for the boxes. This is really an overkill... On one hand,
  ! it is nice to hide everything behind access functions so that we 
  ! preserve the freedom of changing the representation of the types if
  ! needed. Also, this is not a performance-critical process. On the other
  ! hand, the only thing there is to hide here is a pair of chained %'s...
  ! All these could be greatly reduced by a few one-line macros, but at least
  ! for now we choose the conventional approach.
  function get_layout_3D_num_nodes( layout )
    sll_int32                  :: get_layout_3D_num_nodes
    type(layout_3D_t), pointer :: layout
    get_layout_3D_num_nodes = sll_get_collective_size( layout%collective )
  end function get_layout_3D_num_nodes

  function get_layout_3D_box( layout, rank )
    type(box_3D)               :: get_layout_3D_box
    type(layout_3D_t), pointer :: layout
    sll_int32, intent(in)      :: rank
    SLL_ASSERT((rank.ge.0).and.(rank.le.(get_layout_3D_num_nodes(layout)-1)))
    get_layout_3D_box = layout%boxes(rank)
  end function get_layout_3D_box

#define MAKE_GET_LAYOUT_SLOT_FUNCTION( fname, datatype, slot )    \
  function fname( layout, rank );                                 \
    sll_int32                  :: fname;                          \
    type(datatype), pointer    :: layout;                         \
    sll_int32, intent(in)      :: rank;                           \
    fname = layout%boxes(rank)%slot;                              \
  end function fname

#define MAKE_SET_LAYOUT_SLOT_FUNCTION( fname, datatype, slot )    \
  subroutine fname( layout, rank, val );                          \
    type(datatype), pointer    :: layout;                         \
    sll_int32, intent(in)      :: rank;                           \
    sll_int32, intent(in)      :: val;                            \
    layout%boxes(rank)%slot = val;                                \
  end subroutine fname

  ! We use the macros to write the set_ get_ functions for the 3D case as 
  ! well.
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_i_min, layout_3D_t, i_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_i_max, layout_3D_t, i_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_j_min, layout_3D_t, j_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_j_max, layout_3D_t, j_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_k_min, layout_3D_t, k_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_k_max, layout_3D_t, k_max )

  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_i_min, layout_3D_t, i_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_i_max, layout_3D_t, i_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_j_min, layout_3D_t, j_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_j_max, layout_3D_t, j_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_k_min, layout_3D_t, k_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_k_max, layout_3D_t, k_max )

  ! Why should lims just give its collective nilly-willy? This is not 
  ! pretty but I have the suspicion that direct access of the collective 
  ! will be needed.
  function get_layout_3D_collective( layout )
    intrinsic                       :: associated
    type(sll_collective_t), pointer :: get_layout_3D_collective
    type(layout_3D_t), pointer      :: layout
    if( .not. associated(layout) ) then
       stop 'ERROR: uninitialized argument, get_layout_3D_collective()'
    end if
    get_layout_3D_collective => layout%collective
  end function get_layout_3D_collective

  ! get_layout_3D_size() returns the size of the collective associated
  ! with a given layout.
  function get_layout_3D_size( layout )
    intrinsic                  :: associated
    sll_int32                  :: get_layout_3D_size
    type(layout_3D_t), pointer :: layout
    if( .not. associated(layout) ) then
       STOP 'ERROR: not associated argument passed to get_layout_3D_size().'
    end if
    get_layout_3D_size = sll_get_collective_size( layout%collective )
  end function get_layout_3D_size

  ! Utility functions to help build layouts in 3D.

  function linear_index_3D(npx1, npx2, i, j, k)
    sll_int32, intent(in) :: npx1
    sll_int32, intent(in) :: npx2
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32, intent(in) :: k
    sll_int32 :: linear_index_3D
    linear_index_3D = i + npx1*(j + npx2*k)
  end function linear_index_3D



  subroutine initialize_layout_with_distributed_3D_array( &
    global_npx1, &  
    global_npx2, &
    global_npx3, &
    num_proc_x1, &
    num_proc_x2, &
    num_proc_x3, &
    layout_3D )
    
    ! layout_3D should have been allocated with new(), which means that
    ! its memory is allocated in accordance with the size of collective.
    ! This should be error-checked below for consistency.
    sll_int32, intent(in) :: global_npx1
    sll_int32, intent(in) :: global_npx2
    sll_int32, intent(in) :: global_npx3
    sll_int32, intent(in) :: num_proc_x1
    sll_int32, intent(in) :: num_proc_x2
    sll_int32, intent(in) :: num_proc_x3
    type(layout_3D_t), pointer :: layout_3D
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: k
    sll_int32 :: total_num_processors
    sll_int32 :: node
    sll_int32 :: collective_size
    sll_int32 :: err
!    sll_int32, dimension(0:1,0:num_proc_x1-1) :: intervals_x1
!    sll_int32, dimension(0:1,0:num_proc_x2-1) :: intervals_x2
!    sll_int32, dimension(0:1,0:num_proc_x3-1) :: intervals_x3
    sll_int32, dimension(:,:), allocatable :: intervals_x1
    sll_int32, dimension(:,:), allocatable :: intervals_x2
    sll_int32, dimension(:,:), allocatable :: intervals_x3

    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max

    if( &
       .not. is_power_of_two(int(num_proc_x1,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x2,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x3,i64)) ) then
       print *, 'ERROR: distribute_3D_array() needs that the integers that',&
            'describe the process mesh are powers of 2.'
       STOP
    end if

    if( &
       .not. (global_npx1 .gt. 0) .or. &
       .not. (global_npx2 .gt. 0) .or. &
       .not. (global_npx3 .gt. 0) ) then
       print *, 'ERROR: distribute_3D_array() needs that the array dimensions',&
            'be greater than zero.'
       STOP
    end if

    ! FIXME: add further error checking, like a minimum number of points
    ! needed given a processor number along a dimension. Also, num_proc_xi
    ! should be different than zero.
    SLL_ALLOCATE( intervals_x1(0:1,0:num_proc_x1-1), err )
    SLL_ALLOCATE( intervals_x2(0:1,0:num_proc_x2-1), err )
    SLL_ALLOCATE( intervals_x3(0:1,0:num_proc_x3-1), err )

    ! Allocate the layout to be returned.    
    total_num_processors = num_proc_x1*num_proc_x2*num_proc_x3
    collective_size = get_layout_3D_size(layout_3D)
    if( total_num_processors .ne. collective_size ) then
       print *, 'requested size of the processor mesh is inconsistent with ', &
            'the size of the collective.', 'number of processors = ', &
            total_num_processors, ' collective size = ', collective_size
       STOP
    end if

    ! Compute the arrays with the split index information along the different
    ! dimensions.
    intervals_x1(0:1,0:num_proc_x1-1) = &
         split_array_indices( 1, global_npx1, num_proc_x1 )

    intervals_x2(0:1,0:num_proc_x2-1) = &
         split_array_indices( 1, global_npx2, num_proc_x2 )

    intervals_x3(0:1,0:num_proc_x3-1) = &
         split_array_indices( 1, global_npx3, num_proc_x3 )

    ! Fill the layout array.
    do k=0, num_proc_x3-1
       do j=0, num_proc_x2-1
          do i=0, num_proc_x1-1
             node = linear_index_3D( num_proc_x1, num_proc_x2, i, j, k )
             i_min = intervals_x1(0,i)
             i_max = intervals_x1(1,i)
             j_min = intervals_x2(0,j)
             j_max = intervals_x2(1,j)
             k_min = intervals_x3(0,k)
             k_max = intervals_x3(1,k)
             call set_layout_i_min( layout_3D, node, i_min )
             call set_layout_i_max( layout_3D, node, i_max )
             call set_layout_j_min( layout_3D, node, j_min )
             call set_layout_j_max( layout_3D, node, j_max )
             call set_layout_k_min( layout_3D, node, k_min )
             call set_layout_k_max( layout_3D, node, k_max )
          end do
       end do
    end do
    SLL_DEALLOCATE_ARRAY( intervals_x1, err )
    SLL_DEALLOCATE_ARRAY( intervals_x2, err )
    SLL_DEALLOCATE_ARRAY( intervals_x3, err )
   end subroutine initialize_layout_with_distributed_3D_array



  function split_array_indices( min, max, num_intervals )
    sll_int32, intent(in)                       :: num_intervals
    sll_int32, dimension(0:1,0:num_intervals-1) :: split_array_indices
    sll_int32, intent(in)                       :: min
    sll_int32, intent(in)                       :: max
    call split_array_indices_aux( &
      split_array_indices, &
      0, &
      num_intervals, &
      min, &
      max, &
      num_intervals )
  end function split_array_indices

  ! split_array_indices_aux() is an auxiliary function that splits a range
  ! of indices described by 2 integers into a given number of intervals, the
  ! function tries to partition the original interval equitably. 
  recursive subroutine split_array_indices_aux( &
    intervals_array, &
    start_index, &
    interval_segment_length, &
    min, &
    max, &
    total_intervals )

    sll_int32, intent(in) :: total_intervals
    sll_int32, intent(in) :: start_index
    sll_int32, intent(in) :: interval_segment_length
    sll_int32, intent(inout), dimension(0:1,0:total_intervals-1) :: &
         intervals_array
    sll_int32, intent(in) :: min
    sll_int32, intent(in) :: max
    sll_int32 :: num_elems
    sll_int32 :: new_min1
    sll_int32 :: new_max1
    sll_int32 :: new_min2
    sll_int32 :: new_max2
    if( interval_segment_length .eq. 1 ) then
       ! terminate recursion by filling values for this interval
       intervals_array(0,start_index) = min
       intervals_array(1,start_index) = max
    else
       ! split this interval and launch new recursions
       num_elems = max - min + 1
       if( is_even(num_elems) ) then
          new_min1 = min
          new_max1 = min + (max-min+1)/2 - 1
          new_min2 = new_max1 + 1
          new_max2 = max
       else
          new_min1 = min
          new_max1 = min + int((max-min+1)/2)
          new_min2 = new_max1 + 1
          new_max2 = max
       end if
       call split_array_indices_aux( &
          intervals_array, &
          start_index, &
          interval_segment_length/2, &
          new_min1, &
          new_max1, &
          total_intervals )
       call split_array_indices_aux( &
          intervals_array, &
          start_index + interval_segment_length/2, &
          interval_segment_length/2, &
          new_min2,  &
          new_max2,  &
          total_intervals )
    end if
  end subroutine split_array_indices_aux


  ! The new_remap_plan() functions define the communication pattern in a 
  ! collective. From the perspective of an individual process, they examines
  ! the communication needs in a one-to-many and many-to-one sense. The
  ! plan produces information needed to feed a lower level function that
  ! will actually take care of the communications.
  !
  ! To achieve this, new_remap_plan examines two things:
  ! - to whom does an individual process need to send its information, and
  ! - from whom does a process need to receive its information.
  ! Thus, while all processes make this call, the resulting plan will be
  ! different among the processes, but the plan will be consistent in terms 
  ! of which process is expecting what from whom.
  !
  ! The remap plan stores the buffers where the data to be sent/received
  ! are kept. This raises the issue of type dependence. We want to make this
  ! facility as general as possible. Here we try the approach of having a
  ! single, standard format for data storage, i.e. an integer. This means
  ! that we would need to use the transfer() function to store and retrieve
  ! data from the buffers, which is inefficient. The alternative is to
  ! have type-dependent plans... 
  function new_remap_plan_3D( initial, final, int32_data_size )
    intrinsic                       :: associated
    type(remap_plan_3D_t), pointer  :: new_remap_plan_3D 
    type(layout_3D_t), pointer      :: initial
    type(layout_3D_t), pointer      :: final
    sll_int32, intent(in)           :: int32_data_size
    type(sll_collective_t), pointer :: col
    type(box_3D)                    :: ibox, fbox, inters
    sll_int32                       :: i, f          ! loop index
    sll_int32                       :: my_rank
    sll_int32                       :: col_size
    sll_int32                       :: ierr
    sll_int32                       :: disp_counter  ! displacements counter
    sll_int32                       :: send_counter
    sll_int32                       :: recv_counter
    sll_int64                       :: acc = 0       ! accumulator, for sizing

    if( (.not. associated(initial)) .or. (.not. associated(final)) ) then
       write (*,'(a)') &
            'ERROR: un-initialized arguments given to sll_new_remap_plan_3D'
       stop 'sll_new_remap_plan_3D'
    end if
    if( .not. associated(get_layout_3D_collective(initial),&
         target=get_layout_3D_collective(final)) ) then
       write (*,'(a)') &
            'ERROR: init and final configurations given to new_remap_plan do not refer to the same collective.'
       stop 'new_remap_plan_3D'
    end if
  
    col => get_layout_3D_collective(initial)
    my_rank  = sll_get_collective_rank( col )
    col_size = sll_get_collective_size( col )

    SLL_ALLOCATE( new_remap_plan_3D, ierr )
    SLL_CLEAR_ALLOCATE( new_remap_plan_3D%send_displs(0:col_size-1), ierr )
    SLL_CLEAR_ALLOCATE( new_remap_plan_3D%send_counts(0:col_size-1), ierr )
    SLL_CLEAR_ALLOCATE( new_remap_plan_3D%recv_displs(0:col_size-1), ierr )
    SLL_CLEAR_ALLOCATE( new_remap_plan_3D%recv_counts(0:col_size-1), ierr )
    ! Can't use CLEAR_ALLOCATE with types for which the '= 0' has not been
    ! defined as an operator.
    SLL_ALLOCATE( new_remap_plan_3D%send_boxes(0:col_size-1), ierr )
    SLL_ALLOCATE( new_remap_plan_3D%recv_boxes(0:col_size-1), ierr )
    new_remap_plan_3D%collective => get_layout_3D_collective(initial)
    send_counter = 0
    disp_counter = 0
    ibox = get_layout_3D_box(initial, my_rank)
    new_remap_plan_3D%initial_layout => initial
    new_remap_plan_3D%final_layout   => final
    ! Find what data to send. 
    do f = 0, col_size-1  ! loop over the final layout to look for
                          ! box intersections.
    fbox = get_layout_3D_box(final, f)
       if( intersect_boxes_3D( ibox, fbox, inters ) ) then 
          ! compute how many elements to send. 
          send_counter                     = count_elements_in_box_3D(inters)
          new_remap_plan_3D%send_counts(f)   = send_counter
          new_remap_plan_3D%send_displs(f) = disp_counter
          disp_counter                     = disp_counter + send_counter
          new_remap_plan_3D%send_boxes(f)  = inters
          acc                              = acc + send_counter
       else ! if no intersection, there is no need to send anything.
          new_remap_plan_3D%send_counts(f)     = 0
          new_remap_plan_3D%send_displs(f)   = disp_counter
          new_remap_plan_3D%send_boxes(f)    = inters
       end if
    end do
    ! Now that we know the total amount of data to send, we can allocate
    ! the 'send' buffer. Note that we are allocating an integer array, 
    ! so a size adjustment is needed.
    SLL_ALLOCATE(new_remap_plan_3D%send_buffer(0:(acc*int32_data_size-1)),ierr)
    acc          = 0
    disp_counter = 0  
    ! Find what data to receive. Now we compare it with the target layout
    ! for this node.
    fbox = get_layout_3D_box(final, my_rank)
    do i = 0, col_size-1  ! loop over the initial layout to look for box
                          ! intersections.
       ibox = get_layout_3D_box(initial,i)
       if( intersect_boxes_3D( ibox, fbox, inters ) ) then
          ! compute how many elements to receive
          recv_counter                     = count_elements_in_box_3D(inters)
          new_remap_plan_3D%recv_counts(i)   = recv_counter
          new_remap_plan_3D%recv_displs(i) = disp_counter
          disp_counter                     = disp_counter + recv_counter
          new_remap_plan_3D%recv_boxes(i)  = inters
          acc                              = acc + recv_counter
       else ! no intersection, don't expect to receive anything
          new_remap_plan_3D%recv_counts(i)   = 0
          new_remap_plan_3D%recv_displs(i) = disp_counter
          new_remap_plan_3D%recv_boxes(i)  = inters
       end if
    end do
    SLL_ALLOCATE(new_remap_plan_3D%recv_buffer(0:(acc*int32_data_size-1)),ierr)
 !   call optimize_remap_plan_3D(new_remap_plan_3D)
  end function new_remap_plan_3D

  ! The optimizer function is stand-alone. It may be used just
  ! before exiting the new_remap_plan_3D function.
  ! 
  ! The main idea behind the optimizations is the following: 
  ! An unoptimized remap plan is just a call to alltoallV() on a (possibly
  ! very large) communicator. This is concise, but inefficient. For a 
  ! given remapping operation, the re-arrangement of data will 
  ! normally only require communications between subsets of processes 
  ! within this communicator. The optimizer function thus has several 
  ! optimization opportunities:
  !
  ! 1. The first level of optimizations is to identify the subsets of 
  !    processes that communicate with one another, and then launch the
  !    alltoallV() call on this new (hopefully much smaller) communicator.
  ! 2. The second level of optimizations involves identifying those calls to
  !    alltoallV() which exhibit a regular pattern and than can be replaced
  !    by a call to alltoall(), which gives the implementors of the 
  !    communications libraries better optimization opportunities.
  ! 3. The last level of optimization is to simplify the loading of the
  !    send buffers whenever possible.
  !
  ! First optimization level:
  !
  ! We use an algorithm that works as follows:
  ! Every process is aware of the ranks of the processes with which it is
  ! supposed to directly exchange data. This permits the following steps:
  ! a. every process allocates an array of length 'collective_size'
  ! b. each process will figure out the ranks of the processes with which it
  !    communicates (send or receive), and will find the one with the lowest 
  !    rank ( hereafter denominated 'lowest_rank' ) and will set:
  !            array(my_rank) = lowest_rank 
  ! c. the array is shared between all processes with an allgather() operation.
  ! d. we save a copy of this array.
  ! e. then, each process will inspect the fields 'array(i)' of this array that 
  !    correspond to the ranks with which this process is supposed to exchange
  !    (send or receive) any data. The process will find the lowest value and
  !    will set     array(my_rank) = lowest_value.
  ! f. compare with the stored version of this array. If there was no change,
  !    we are done.
  ! g. else, repeat the process starting from d.
  ! h. Split new collectives using array(my_rank) as the color. These are
  !    the sub-collectives that we were looking for.
  !
  ! This algorithm will be able to sort out through complicated communication
  ! patterns where the exchanges are asymmetric.
  !
  ! The optimized plan also compresses the 'box' arrays that store the 
  ! information on the data that is to be sent/received, as well as the 
  ! send_counts and displacements...  This introduces a
  ! problem: The optimized plan has a notion of a reduced collective, as well
  ! as compressed box arrays, but will still need the 'global' information
  ! about the layouts, since the global_to_local function only has meaning 
  ! in the context of the global layout. I find this 'mixing' very unpleasant, 
  ! and it might invite confusing those two collectives (the parent and the
  ! reduced); for now see no clean & easy way to fix this. Fortunately, at 
  ! least, the reference to the larger collective is hidden inside the 'layout' 
  ! information and used only by 'global_to_local()'.
  subroutine optimize_remap_plan_3D( plan ) !, sub_collective )
    type(remap_plan_3D_t), pointer       :: plan
    sll_int32, dimension(:), pointer     :: send_counts
    sll_int32, dimension(:), pointer     :: send_displs
    sll_int32, dimension(:), pointer     :: recv_counts
    sll_int32, dimension(:), pointer     :: recv_displs
    type(sll_collective_t), pointer      :: col
    sll_int32                            :: col_sz
    sll_int32, dimension(:), allocatable :: lowest_color
    sll_int32, dimension(:), allocatable :: colors
    sll_int32, dimension(:), allocatable :: colors_copy
    sll_int32                            :: ierr
    sll_int32                            :: my_rank
    sll_int32                            :: i
    type(sll_collective_t), pointer      :: new_collective
    sll_int32                            :: new_col_sz
    sll_int32, dimension(:), pointer     :: new_send_counts
    sll_int32, dimension(:), pointer     :: new_send_displs
    sll_int32, dimension(:), pointer     :: new_recv_counts
    sll_int32, dimension(:), pointer     :: new_recv_displs
    type(box_3D), dimension(:), pointer  :: new_send_boxes
    type(box_3D), dimension(:), pointer  :: new_recv_boxes
    sll_int32                            :: new_i
    sll_int32                            :: my_color
    sll_int32                            :: exchange_size
    logical, dimension(1:1)              :: is_uniform_local
    logical, dimension(1:1)              :: is_uniform_collective

    col         => plan%collective
    col_sz      = sll_get_collective_size( col )
    my_rank     = sll_get_collective_rank( col )
    send_counts => plan%send_counts
    send_displs => plan%send_displs
    recv_counts => plan%recv_counts
    recv_displs => plan%recv_displs
    SLL_CLEAR_ALLOCATE( lowest_color(1), ierr ) ! awful, but I need an array.
    SLL_CLEAR_ALLOCATE( colors(0:col_sz-1), ierr )
    SLL_CLEAR_ALLOCATE( colors_copy(0:col_sz-1), ierr )

    ! FIRST LEVEL OF OPTIMIZATION: 
    ! Identify the sub-collectives in which the communication should 
    ! be divided.
    lowest_color(1) = my_rank  ! we want a starting point. This should not 
                               ! change for the lowest ranks in the 
                               ! sub-collectives.
    do
       ! Load the copy
       colors_copy(0:col_sz-1) = colors(0:col_sz-1)
       ! Find the lowest rank with which this process communicates
       do i=0,col_sz-1
          if( (send_counts(i) .ne. 0) .or. (recv_counts(i) .ne. 0) ) then
             if( i .lt. lowest_color(1) ) then
                lowest_color(1) = i
             end if  ! else we don't do anything as lowest_color is already 
                     ! my_rank
          end if
       end do
       call sll_collective_allgather( col, lowest_color(1:1), 1, &
                                      colors(0:col_sz-1), 1 )
       if(arrays_are_equal(colors, colors_copy, col_sz)) exit
    end do
    ! The results can now be used as the color for a splitting operation.
    new_collective => sll_new_collective( col, colors(my_rank), my_rank )
    new_col_sz     = sll_get_collective_size( new_collective )
    ! Allocate the new counters and displacements with the reduced 
    ! collective size.
    SLL_ALLOCATE( new_send_counts(0:new_col_sz-1), ierr )
    SLL_ALLOCATE( new_send_displs(0:new_col_sz-1), ierr )
    SLL_ALLOCATE( new_recv_counts(0:new_col_sz-1), ierr )
    SLL_ALLOCATE( new_recv_displs(0:new_col_sz-1), ierr )
    SLL_ALLOCATE( new_send_boxes( 0:new_col_sz-1), ierr )
    SLL_ALLOCATE( new_recv_boxes( 0:new_col_sz-1), ierr )
    ! Compress the 'send' and 'receive' information
    new_i = 0
    my_color = colors(my_rank)
    do i=0,col_sz-1
       if( colors(i) .eq. my_color ) then
          new_send_counts(new_i) = send_counts(i)
          new_send_displs(new_i) = send_displs(i)
          new_send_boxes(new_i)  = plan%send_boxes(i)
          new_recv_counts(new_i) = recv_counts(i)
          new_recv_displs(new_i) = recv_displs(i)
          new_recv_boxes(new_i)  = plan%recv_boxes(i)
          new_i                  = new_i + 1
       end if
    end do
    ! Change the fields of the plan to reflect the new:
    ! - collective,
    ! - send_counts,
    ! - send_displs,
    ! - recv_counts,
    ! - recv_displs,
    ! - send_boxes, and
    ! - recv_boxes.
    ! The send/receive buffers remain unchanged. Watch out for possible
    ! memory leaks...
    plan%collective => new_collective
    SLL_DEALLOCATE( plan%send_counts, ierr )
    plan%send_counts => new_send_counts
    SLL_DEALLOCATE( plan%send_displs, ierr )
    plan%send_displs => new_send_displs
    SLL_DEALLOCATE( plan%recv_counts, ierr )
    plan%recv_counts => new_recv_counts
    SLL_DEALLOCATE( plan%recv_displs, ierr )
    plan%recv_displs => new_recv_displs
    SLL_DEALLOCATE( plan%send_boxes, ierr )
    plan%send_boxes => new_send_boxes
    SLL_DEALLOCATE( plan%recv_boxes, ierr )
    plan%recv_boxes => new_recv_boxes
    SLL_DEALLOCATE_ARRAY( lowest_color, ierr )
    SLL_DEALLOCATE_ARRAY( colors, ierr )
    SLL_DEALLOCATE_ARRAY( colors_copy, ierr )

    ! SECOND LEVEL OF OPTIMIZATION:
    ! Identify whether the communication in the local collective is regular,
    ! thus permitting a call to alltoall(). Note that it is not sufficient
    ! to detect whether the exchange looks regular locally, all processes
    ! in the communicator must agree in this view.
    exchange_size = plan%send_counts(0)
    do i=0,new_col_sz-1
       if(plan%send_counts(i) .eq. exchange_size) then
          is_uniform_local(1) = .true.
       else ! plan%send_counts(i) is different than the first value
          is_uniform_local(1) = .false.
          exit  
       end if
    end do
    ! Use a reduction operation to find out if this result is shared with 
    ! the other processes in the collective. Hmmm... look at this slightly
    ! disastrous occurrence: the MPI reduction operation MPI_LAND got out of
    ! the cage... this needs to be addressed.
    call sll_collective_allreduce(plan%collective, is_uniform_local(:), 1, &
         MPI_LAND, is_uniform_collective(:) )
    plan%is_uniform = is_uniform_collective(1)
    ! This flag will be used for an optimized call in apply_remap_plan()
#if 0
    write (*,'(a,i4)') 'collective color: ', &
         sll_get_collective_color(plan%collective)
    write (*,'(a,i4)') 'my_rank in new collective: ', &
         sll_get_collective_rank(plan%collective)
    write (*,'(a,i4)') 'new collective size: ', &
         sll_get_collective_size(plan%collective)
    print *, plan%send_counts(:)
    print *, plan%send_displs(:)
    print *, plan%recv_counts(:)
    print *, plan%recv_displs(:)
    call flush()
#endif
  end subroutine optimize_remap_plan_3D

  function arrays_are_equal( a1, a2, n )
    logical :: arrays_are_equal
    sll_int32, dimension(:), intent(in) :: a1
    sll_int32, dimension(:), intent(in) :: a2
    sll_int32, intent(in)               :: n
    sll_int32                           :: i
    arrays_are_equal = .true.
    do i=1,n
       if( a1(i).ne.a2(i) ) then
          arrays_are_equal = .false. 
       end if
    end do
  end function arrays_are_equal

  function get_remap_3D_initial_layout( plan )
    type(layout_3D_t), pointer      :: get_remap_3D_initial_layout
    type(remap_plan_3D_t), pointer  :: plan
    if( .not. associated( plan ) ) then
       write (*,'(a)') 'not associated pointer argument'
       stop 'get_remap_3D_initial_layout'
    end if
    get_remap_3D_initial_layout => plan%initial_layout
  end function get_remap_3D_initial_layout

  function get_remap_3D_final_layout( plan )
    type(layout_3D_t), pointer      :: get_remap_3D_final_layout
    type(remap_plan_3D_t), pointer  :: plan
    if( .not. associated( plan ) ) then
       write (*,'(a)') 'not associated pointer argument'
       stop 'get_remap_3D_final_layout'
    end if
    get_remap_3D_final_layout => plan%final_layout
  end function get_remap_3D_final_layout

  ! In this implementation, the user provides the memory location where the
  ! result of the remap operation will end up. The remap functions are
  ! type dependent. The objective will be to keep the type dependence here
  ! only. This objective is straightforward when we want to send types that
  ! happen to be basic MPI types. A problem arises when we want to communicate
  ! something like a derived type.
  !
  ! For derived types, one would normally be required to use the MPI 
  ! derived type functions like MPI_Type_struct and all that (or whatever
  ! their Fortran equivalent is) and all of a sudden we have lost containment
  ! (read: lost modularity) of the MPI LIBRARY.
  !
  ! Here we try an approach, standard in C but apparently unusual in Fortran.
  ! The idea is to always communicate the MPI_INTEGER datatype. The only extra
  ! step is to translate the arrays for counts and displacements in terms
  ! of their integer-sizes. While we would still need a different 
  ! apply_remap_XD() function for every new type, at least this will not 
  ! affect the new_remap_plan() function.
  !
  ! For apply_remap_XD_int(), we use this approach as a test case, even though
  ! it is not necessary since 'integer' is also a native MPI type.

  subroutine convert_into_integer_sizes( sz, ai, n, bi )
    sll_int32, intent(in)                :: sz   ! in integer-size
    sll_int32, intent(in), dimension(:)  :: ai   ! array to convert
    sll_int32, intent(in)                :: n    ! size of array
    sll_int32, intent(out), dimension(:) :: bi   ! output
    sll_int32                            :: i
    do i=1,n
       bi(i) = ai(i)*sz
    end do
  end subroutine convert_into_integer_sizes

  ! **********************************************************************
  !
  !    Continue here with 3D, 4D and 5D access functions...
  !
  ! **********************************************************************
  subroutine apply_remap_3D_int( plan, data_in, data_out )
    intrinsic                                :: transfer
    type(remap_plan_3D_t), pointer           :: plan
    sll_int32, dimension(:,:,:), intent(in)  :: data_in
    sll_int32, dimension(:,:,:), intent(out) :: data_out
    sll_int32, dimension(:), pointer         :: sb       ! send buffer
    sll_int32, dimension(:), pointer         :: rb       ! receive buffer
    sll_int32, dimension(:), pointer         :: sdisp    ! send displacements
    sll_int32, dimension(:), pointer         :: rdisp    ! receive displacements
    sll_int32, dimension(:), pointer         :: scnts    ! send counts
    sll_int32, dimension(:), pointer         :: rcnts    ! receive counts
    type(sll_collective_t), pointer          :: col      ! collective
    type(layout_3D_t), pointer               :: init_layout  => NULL()
    type(layout_3D_t), pointer               :: final_layout => NULL()
    sll_int32                                :: id, jd, kd
    sll_int32                                :: i
    sll_int32                                :: col_sz
    sll_int32                                :: ierr
    sll_int32                                :: loi, loj, lok
    sll_int32                                :: hii, hij, hik
    type(box_3D)                             :: sbox
    sll_int32                                :: my_rank
    sll_int32                                :: loc
    sll_int32, dimension(1:3)                :: local_lo, local_hi
    sll_int32                                :: int32_data_size

    ! to load the MPI function and send integers, we have a separate set of
    ! arrays to store this information for now.
    sll_int32, dimension(:), allocatable     :: sdispi  ! send displacements
    sll_int32, dimension(:), allocatable     :: rdispi  ! receive displacements
    sll_int32, dimension(:), allocatable     :: scntsi  ! send counts
    sll_int32, dimension(:), allocatable     :: rcntsi  ! receive counts

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp        => plan%send_displs
    rdisp        => plan%recv_displs
    scnts        => plan%send_counts
    rcnts        => plan%recv_counts
    col          => plan%collective
    col_sz       =  sll_get_collective_size(col)
    init_layout  => get_remap_3D_initial_layout(plan)
    final_layout => get_remap_3D_final_layout(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer

    SLL_ALLOCATE(sdispi(0:col_sz-1), ierr)
    SLL_ALLOCATE(rdispi(0:col_sz-1), ierr)
    SLL_ALLOCATE(scntsi(0:col_sz-1), ierr)
    SLL_ALLOCATE(rcntsi(0:col_sz-1), ierr)

    ! Translate the amounts into integers
#if 1
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1,1)), sdisp, &
         col_sz, sdispi)
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1,1)), rdisp, &
         col_sz, rdispi)
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1,1)), scnts, &
         col_sz, scntsi)
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1,1)), rcnts, &
         col_sz, rcntsi)
#endif
    
#if 0
    write (*,'(a,i4)') 'parameters from rank ', my_rank
    print *, scntsi(:)
    print *, sdispi(:)
    print *, rcntsi(:)
    print *, rdispi(:)
#endif
    
    ! load the send buffer
    loc = 0             ! first loading is at position zero
    ! This step is obviously not needed for integers themselves. We put this
    ! here for generality.
    int32_data_size = INT32_SIZEOF( data_in(1,1,1) )
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdispi(i) ) then
             write (*,'(a,i4,a,i16)') 'apply_remap_3D_int() ERROR: discrepancy between displs(i) and the loading index for i = ', i, ' displs(i) = ', sdispi(i)
             write(*,'(a,i8)') 'col_sz = ', col_sz
             stop 'apply_remap(): loading error'
          end if
          ! get the information on the box to send, get the limits,
          ! convert to the local indices and find out where in the 
          ! buffer to start writing.
          sbox = plan%send_boxes(i)
          loi = get_box_3D_i_min(sbox)
          loj = get_box_3D_j_min(sbox)
          lok = get_box_3D_k_min(sbox)
          hii = get_box_3D_i_max(sbox)
          hij = get_box_3D_j_max(sbox)
          hik = get_box_3D_k_max(sbox)
          local_lo = global_to_local_3D( init_layout, (/loi,loj,lok/) )
          local_hi = global_to_local_3D( init_layout, (/hii,hij,hik/) )

          ! The plan to load the send buffer is to traverse the integer
          ! array with a single index (loc). When we load the buffer, each
          ! data element may occupy multiple integer 'slots', hence the
          ! loading index needs to be manually increased. As an advantage,
          ! we can do some error checking every time we send data to a 
          ! different process, as we know what is the expected value of 
          ! the index at that point.
          do kd = local_lo(3), local_hi(3)
             do jd = local_lo(2), local_hi(2)
                do id = local_lo(1), local_hi(1)
                   sb(loc:) = transfer(data_in(id,jd,kd),(/1_i32/))
                   loc      = loc + int32_data_size
                end do
             end do
          end do
       end if
    end do
    
!    write (*,'(a,i4)') 'the send buffer in rank:', my_rank
!    print *, sb(0:(size(sb)-1))
!    call flush()
 
   if( plan%is_uniform .eqv. .false. ) then 
       call sll_collective_alltoallV( sb(:),       &
                                      scntsi(0:col_sz-1), &
                                      sdispi(0:col_sz-1), &
                                      rb(:),       &
                                      rcntsi(0:col_sz-1), &
                                      rdispi(0:col_sz-1), col )
    else
       call sll_collective_alltoall ( sb(:), &
                                      scntsi(0), &
                                      rcntsi(0), &
                                      rb(:), col )
    end if
!    write (*,'(a, i4)') 'the receive buffer in rank: ', my_rank
!    print *, rb(0:size(rb)-1)
!    call flush()
    ! Unpack the plan into the outgoing buffer.
    loc = 0  ! We load first from position 0 in the receive buffer.
    do i = 0, col_sz-1
       if( rcnts(i) .ne. 0 ) then ! we expect something from rank 'i'
          if( loc .ne. rdispi(i) ) then
             write (*,'(a,i4)') &
                  'ERROR: discrepancy between rdispi(i) and index for i = ', i
             stop 'unpacking error'
          end if
          ! get the information on the box to receive, get the limits, and 
          ! convert to the local indices.
          sbox = plan%recv_boxes(i)
          loi = get_box_3D_i_min(sbox)
          loj = get_box_3D_j_min(sbox)
          lok = get_box_3D_k_min(sbox)
          hii = get_box_3D_i_max(sbox)
          hij = get_box_3D_j_max(sbox)
          hik = get_box_3D_k_max(sbox)
          local_lo = global_to_local_3D( final_layout, (/loi,loj,lok/) )
          local_hi = global_to_local_3D( final_layout, (/hii,hij,hik/) )
          do kd = local_lo(3), local_hi(3)
             do jd = local_lo(2), local_hi(2)
                do id = local_lo(1), local_hi(1)
                   data_out(id,jd,kd) = transfer(rb(loc:),data_out(1,1,1))
                   loc                = loc + int32_data_size
                end do
             end do
          end do
       end if
    end do
  end subroutine apply_remap_3D_int

  subroutine apply_remap_3D_double( plan, data_in, data_out )
    intrinsic                                 :: transfer
    type(remap_plan_3D_t), pointer            :: plan
    sll_real64, dimension(:,:,:), intent(in)  :: data_in
    sll_real64, dimension(:,:,:), intent(out) :: data_out
    sll_int32, dimension(:), pointer          :: sb     ! send buffer
    sll_int32, dimension(:), pointer          :: rb     ! receive buffer
    sll_int32, dimension(:), pointer          :: sdisp  ! send displacements
    sll_int32, dimension(:), pointer          :: rdisp  ! receive displacements 
    sll_int32, dimension(:), pointer          :: scnts  ! send counts
    sll_int32, dimension(:), pointer          :: rcnts  ! receive counts
    type(sll_collective_t), pointer           :: col    ! collective
    type(layout_3D_t), pointer                :: init_layout  => NULL()
    type(layout_3D_t), pointer                :: final_layout => NULL()
    sll_int32                                 :: id, jd, kd
    sll_int32                                 :: i
    sll_int32                                 :: col_sz
    sll_int32                                 :: ierr
    sll_int32                                 :: loi, loj, lok
    sll_int32                                 :: hii, hij, hik
    type(box_3D)                              :: sbox
    sll_int32                                 :: my_rank
    sll_int32                                 :: loc
    sll_int32, dimension(1:3)                 :: local_lo, local_hi
    sll_int32                                 :: int32_data_size

    ! to load the MPI function and send integers, we have a separate set of
    ! arrays to store this information for now.
    sll_int32, dimension(:), allocatable     :: sdispi  ! send displacements
    sll_int32, dimension(:), allocatable     :: rdispi  ! receive displacements
    sll_int32, dimension(:), allocatable     :: scntsi  ! send counts
    sll_int32, dimension(:), allocatable     :: rcntsi  ! receive counts

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp        => plan%send_displs
    rdisp        => plan%recv_displs
    scnts        => plan%send_counts
    rcnts        => plan%recv_counts
    col          => plan%collective
    col_sz       =  sll_get_collective_size(col)
    init_layout  => get_remap_3D_initial_layout(plan)
    final_layout => get_remap_3D_final_layout(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer

    SLL_ALLOCATE(sdispi(0:col_sz-1), ierr)
    SLL_ALLOCATE(rdispi(0:col_sz-1), ierr)
    SLL_ALLOCATE(scntsi(0:col_sz-1), ierr)
    SLL_ALLOCATE(rcntsi(0:col_sz-1), ierr)
    ! print *, 'from rank ', my_rank, 'loading parameters: ', sdisp, rdisp, &
    ! scnts, rcnts

    ! Translate the amounts into integers
#if 1
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1,1)), sdisp, &
         col_sz, sdispi)
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1,1)), rdisp, &
         col_sz, rdispi)
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1,1)), scnts, &
         col_sz, scntsi)
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1,1)), rcnts, &
         col_sz, rcntsi)
#endif
    
#if 1
    write (*,'(a,i4)') 'parameters from rank ', my_rank
    print *, 'scntsi', scntsi(:)
    print *, 'sdispi', sdispi(:)
    print *, 'rcntsi', rcntsi(:)
    print *, 'rdispi', rdispi(:)
    call flush()
#endif
    
    ! load the send buffer
    loc = 0             ! first loading is at position zero
    ! This step is obviously not needed for integers themselves. We put this
    ! here for generality.
    int32_data_size = INT32_SIZEOF( data_in(1,1,1) )
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdispi(i) ) then
             print *, 'ERROR DETECTED in process: ', my_rank
             print *, 'apply_remap_3D_double() ERROR: ', &
                  'discrepancy between displs(i) and the loading index for ',&
                  'i = ', i, ' displs(i) = ', sdispi(i)
             write(*,'(a,i8)') 'col_sz = ', col_sz
             call flush()
             stop 'apply_remap(): loading error'
          end if
          ! get the information on the box to send, get the limits,
          ! convert to the local indices and find out where in the 
          ! buffer to start writing.
          sbox = plan%send_boxes(i)
          loi = get_box_3D_i_min(sbox)
          loj = get_box_3D_j_min(sbox)
          lok = get_box_3D_k_min(sbox)
          hii = get_box_3D_i_max(sbox)
          hij = get_box_3D_j_max(sbox)
          hik = get_box_3D_k_max(sbox)
          local_lo = global_to_local_3D( init_layout, (/loi,loj,lok/) )
          local_hi = global_to_local_3D( init_layout, (/hii,hij,hik/) )

          ! The plan to load the send buffer is to traverse the integer
          ! array with a single index (loc). When we load the buffer, each
          ! data element may occupy multiple integer 'slots', hence the
          ! loading index needs to be manually increased. As an advantage,
          ! we can do some error checking every time we send data to a 
          ! different process, as we know what is the expected value of 
          ! the index at that point.
          do kd = local_lo(3), local_hi(3)
             do jd = local_lo(2), local_hi(2)
                do id = local_lo(1), local_hi(1)
                   sb(loc:) = transfer(data_in(id,jd,kd),(/1_i32/))
                   loc      = loc + int32_data_size
                end do
             end do
          end do
       end if
    end do
    
!    write (*,'(a,i4)') 'the send buffer in rank:', my_rank
!    print *, sb(0:(size(sb)-1))
!    call flush()
 
   if( plan%is_uniform .eqv. .false. ) then 
       call sll_collective_alltoallV( sb(:),       &
                                      scntsi(0:col_sz-1), &
                                      sdispi(0:col_sz-1), &
                                      rb(:),       &
                                      rcntsi(0:col_sz-1), &
                                      rdispi(0:col_sz-1), col )
    else
       call sll_collective_alltoall ( sb(:), &
                                      scntsi(0), &
                                      rcntsi(0), &
                                      rb(:), col )
    end if
!    write (*,'(a, i4)') 'the receive buffer in rank: ', my_rank
!    print *, rb(0:size(rb)-1)
!    call flush()
    ! Unpack the plan into the outgoing buffer.
    loc = 0  ! We load first from position 0 in the receive buffer.
    do i = 0, col_sz-1
       if( rcnts(i) .ne. 0 ) then ! we expect something from rank 'i'
          if( loc .ne. rdispi(i) ) then
             write (*,'(a,i4)') &
                  'ERROR: discrepancy between rdispi(i) and index for i = ', i
             stop 'unpacking error'
          end if
          ! get the information on the box to receive, get the limits, and 
          ! convert to the local indices.
          sbox = plan%recv_boxes(i)
          loi = get_box_3D_i_min(sbox)
          loj = get_box_3D_j_min(sbox)
          lok = get_box_3D_k_min(sbox)
          hii = get_box_3D_i_max(sbox)
          hij = get_box_3D_j_max(sbox)
          hik = get_box_3D_k_max(sbox)
          local_lo = global_to_local_3D( final_layout, (/loi,loj,lok/) )
          local_hi = global_to_local_3D( final_layout, (/hii,hij,hik/) )
          do kd = local_lo(3), local_hi(3)
             do jd = local_lo(2), local_hi(2)
                do id = local_lo(1), local_hi(1)
                   data_out(id,jd,kd) = transfer(rb(loc:),data_out(1,1,1))
                   loc                = loc + int32_data_size
                end do
             end do
          end do
       end if
    end do

    SLL_DEALLOCATE_ARRAY(sdispi, ierr)
    SLL_DEALLOCATE_ARRAY(rdispi, ierr)
    SLL_DEALLOCATE_ARRAY(scntsi, ierr)
    SLL_DEALLOCATE_ARRAY(rcntsi, ierr)
  end subroutine apply_remap_3D_double


  ! Placeholder: the load/unload subroutines for the send and receive buffers
  ! should ideally be abstracted out. This means that we need to probably
  ! define some generic interface and hide behind the types that we want to
  ! exchange. We can think of direct loading and exchanging the basic types:
  ! - single/double precision floats
  ! - integers
  ! - stay with the transfer function for other types.
#if 0
  subroutine load_send_buffer_3D( plan, data_in )
    type(remap_plan_3D_t), pointer :: plan
    sll_int32                      :: load_point
  end subroutine load_send_buffer_3D
#endif

#define MAKE_GET_BOX_SLOT_FUNCTION( fname, boxtype, slot )   \
  function fname( b );                                       \
    sll_int32                 :: fname;                      \
    type(boxtype), intent(in) :: b;                          \
    fname = b%slot;                                          \
  end function fname

  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_i_min, box_3D, i_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_i_max, box_3D, i_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_j_min, box_3D, j_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_j_max, box_3D, j_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_k_min, box_3D, k_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_k_max, box_3D, k_max )

  function sll_get_num_nodes_3D( lims )
    sll_int32                  :: sll_get_num_nodes_3D
    type(layout_3D_t), pointer :: lims
    sll_get_num_nodes_3D = sll_get_collective_size( lims%collective )
  end function sll_get_num_nodes_3D

  ! It seems that it is essential to have functions that would convert
  ! indices from a local to a global indexing and back. For instance, if we
  ! consider a 1D array as a single global unit, it has an unique indexing.
  ! If we domain-decompose the array, even if there are overlapping domains,
  ! the result of a local2global() indexing function would be unambiguous.
  !
  ! However, the global2local() indexing operation can be ambiguous if the
  ! domains are overlapping. There are thus some design choices that emerge:
  !
  ! On the one hand, we could choose a scheme in which the global2local()
  ! function only responds in the context of the calling process, thus, if
  ! a call to this function returns, say, '0' (and the array is indexed 1:N), 
  ! then we know that the global index has no presence in the layout in the 
  ! local process.
  ! 
  ! On the other hand, the global2local() function could also be made to
  ! indicate the local index and rank of the other processes that also happen
  ! to have that global index in their domain. This second part, while it
  ! seems more useful, may be way too complicated (and would return multiple
  ! things). Thus here we go for the first option.
  function local_to_global_3D( layout, triplet )
    sll_int32, dimension(1:3)             :: local_to_global_3D
    type(layout_3D_t), pointer            :: layout
    sll_int32, intent(in), dimension(1:3) :: triplet
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: my_rank
    type(box_3D)                          :: box
    ! fixme: arg checking
    col                => get_layout_3D_collective( layout )
    my_rank            =  sll_get_collective_rank( col )
    box                =  get_layout_3D_box( layout, my_rank )
    local_to_global_3D(1) = get_box_3D_i_min(box) + triplet(1) - 1
    local_to_global_3D(2) = get_box_3D_j_min(box) + triplet(2) - 1
    local_to_global_3D(3) = get_box_3D_k_min(box) + triplet(3) - 1
  end function local_to_global_3D


  ! We need to make sure that the decision of choosing '-1' as the return
  ! value when the global index is not available locally does not backfire.
  ! If one decides to use an array with an indexing that contains -1, this 
  ! would be problematic.
  function global_to_local_3D( layout, gtuple )
    intrinsic                             :: associated
    sll_int32, dimension(1:3)             :: global_to_local_3D
    type(layout_3D_t), pointer            :: layout
    sll_int32, dimension(1:3), intent(in) :: gtuple ! global indices, as array
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: my_rank
    type(box_3D)                          :: box
    if( .not. associated(get_layout_3D_collective(layout)) ) then
       write (*,'(a)') 'ERROR in global_to_local_3D(), not-associated col'
       stop 'global_to_local_3D'
    end if
    col     => get_layout_3D_collective( layout )
    my_rank =  sll_get_collective_rank( col )
    box     =  get_layout_3D_box( layout, my_rank )
    if( (gtuple(1) .ge. get_box_3D_i_min(box)) .and. &
         (gtuple(1) .le. get_box_3D_i_max(box)) .and. &
         (gtuple(2) .ge. get_box_3D_j_min(box)) .and. &
         (gtuple(2) .le. get_box_3D_j_max(box)) .and. &
         (gtuple(3) .ge. get_box_3D_k_min(box)) .and. &
         (gtuple(3) .le. get_box_3D_k_max(box)) ) then  ! the index is present
       global_to_local_3D(1) = gtuple(1) - get_box_3D_i_min(box) + 1
       global_to_local_3D(2) = gtuple(2) - get_box_3D_j_min(box) + 1
       global_to_local_3D(3) = gtuple(3) - get_box_3D_k_min(box) + 1
    else  ! the index is not present
       global_to_local_3D(1) = -1
       global_to_local_3D(2) = -1
       global_to_local_3D(3) = -1
    end if
  end function global_to_local_3D

  subroutine view_box_3D( b )
    type(box_3D), intent(in) :: b
    write(*,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') &
         '[  [', b%i_min,',', b%i_max,'], [', &
                 b%j_min,',', b%j_max,'], [', &
                 b%k_min,',', b%k_max,']  ]'
  end subroutine view_box_3D

  subroutine sll_view_lims_3D( layout )
    type(layout_3D_t), pointer :: layout
    sll_int32                  :: i
    sll_int32                  :: sz
    sz = sll_get_num_nodes( layout )
    print *, 'limits: '
    do i=0,sz-1
       call view_box_3D(get_layout_3D_box( layout, i ))
!       write(*,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') &
!            '[  [', get_layout_i_min(lims,i),',',&
!                  get_layout_i_max(lims,i),'], [', &
!                  get_layout_j_min(lims,i),',',    &
!                  get_layout_j_max(lims,i),'], [', &
!                  get_layout_k_min(lims,i),',',    &
!                  get_layout_k_max(lims,i),']  ]'
    end do
    call flush()
  end subroutine sll_view_lims_3D

  ! the return value of intersect_boxes() is 'logical', and answers
  ! the question whether the boxes intersect or not. 'ans' is a box with
  ! the actual intersection between the argument boxes. In case that there
  ! is no intersection between the boxes the value [0,0] is returned. 
  function intersect_boxes_3D( b1, b2, ans )
    intrinsic                 :: min, max
    logical                   :: intersect_boxes_3D
    type(box_3D), intent(in)  :: b1, b2
    type(box_3D), intent(out) :: ans
    sll_int32                 :: loi, hii
    sll_int32                 :: loj, hij
    sll_int32                 :: lok, hik
    sll_int32                 :: loib1, hiib1
    sll_int32                 :: lojb1, hijb1
    sll_int32                 :: lokb1, hikb1
    sll_int32                 :: loib2, hiib2
    sll_int32                 :: lojb2, hijb2
    sll_int32                 :: lokb2, hikb2
    ! FIXME: add error checking, if boxes are null, for instance.
    loib1 = get_box_3D_i_min(b1)
    hiib1 = get_box_3D_i_max(b1)
    lojb1 = get_box_3D_j_min(b1)
    hijb1 = get_box_3D_j_max(b1)
    lokb1 = get_box_3D_k_min(b1)
    hikb1 = get_box_3D_k_max(b1)

    loib2 = get_box_3D_i_min(b2)
    hiib2 = get_box_3D_i_max(b2)
    lojb2 = get_box_3D_j_min(b2)
    hijb2 = get_box_3D_j_max(b2)
    lokb2 = get_box_3D_k_min(b2)
    hikb2 = get_box_3D_k_max(b2)

    SLL_ASSERT( (loib1 .le. hiib1) .and. (loib2 .le. hiib2) )
    SLL_ASSERT( (lojb1 .le. hijb1) .and. (lojb2 .le. hijb2) )
    SLL_ASSERT( (lokb1 .le. hikb1) .and. (lokb2 .le. hikb2) )

    loi = max(loib1, loib2)
    hii = min(hiib1, hiib2)
    loj = max(lojb1, lojb2)
    hij = min(hijb1, hijb2)
    lok = max(lokb1, lokb2)
    hik = min(hikb1, hikb2)

    if( (loi .gt. hii) .or. (loj .gt. hij) .or. (lok .gt. hik) ) then 
       ans%i_min = 0
       ans%i_max = 0
       ans%j_min = 0
       ans%j_max = 0
       ans%k_min = 0
       ans%k_max = 0
       intersect_boxes_3D = .false.
    else
       ans%i_min = loi
       ans%i_max = hii
       ans%j_min = loj
       ans%j_max = hij
       ans%k_min = lok
       ans%k_max = hik
       intersect_boxes_3D = .true.
    end if
  end function intersect_boxes_3D

  function count_elements_in_box_3D( box )
    sll_int32                 :: count_elements_in_box_3D
    type(box_3D), intent (in) :: box
    sll_int32                 :: irange
    sll_int32                 :: jrange
    sll_int32                 :: krange
    irange = get_box_3D_i_max(box) - get_box_3D_i_min(box) + 1
    jrange = get_box_3D_j_max(box) - get_box_3D_j_min(box) + 1
    krange = get_box_3D_k_max(box) - get_box_3D_k_min(box) + 1
    count_elements_in_box_3D = irange*jrange*krange
  end function count_elements_in_box_3D


end module remapper
