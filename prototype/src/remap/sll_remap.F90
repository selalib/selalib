module remapper
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "misc_utils.h"
  use sll_collective
  implicit none
  
  ! The box types contain information on the index limits contained        
  ! in a given processor.

  type box_1D
     sll_int32 :: i_min, i_max
  end type box_1D
  
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

  type layout_1D_t
     type(sll_collective_t), pointer     :: collective
     type(box_1D), dimension(:), pointer :: boxes
  end type layout_1D_t

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
  type remap_plan_1D_t
     type(layout_1D_t), pointer              :: initial_layout
     type(layout_1D_t), pointer              :: final_layout
     integer, dimension(:), allocatable      :: send_displs
     integer, dimension(:), allocatable      :: send_cnts
     integer, dimension(:), allocatable      :: recv_displs
     integer, dimension(:), allocatable      :: recv_cnts
     type(box_1D), dimension(:), allocatable :: send_boxes
     type(box_1D), dimension(:), allocatable :: recv_boxes
     type(sll_collective_t), pointer         :: collective
  end type remap_plan_1D_t

  type remap_plan_3D_t
     type(layout_3D_t), pointer              :: initial_layout
     type(layout_3D_t), pointer              :: final_layout
     integer, dimension(:), allocatable      :: send_displs
     integer, dimension(:), allocatable      :: send_cnts
     integer, dimension(:), allocatable      :: recv_displs
     integer, dimension(:), allocatable      :: recv_cnts
     type(box_3D), dimension(:), allocatable :: send_boxes
     type(box_3D), dimension(:), allocatable :: recv_boxes
     type(sll_collective_t), pointer         :: collective
     sll_int32, dimension(:), allocatable    :: send_buffer
     sll_int32, dimension(:), allocatable    :: recv_buffer
  end type remap_plan_3D_t



 interface get_layout_i_min
     module procedure get_layout_1D_i_min, get_layout_3D_i_min
  end interface

  interface set_layout_i_min
     module procedure set_layout_1D_i_min, set_layout_3D_i_min
  end interface

  interface get_layout_i_max
     module procedure get_layout_1D_i_max, get_layout_3D_i_max
  end interface

  interface set_layout_i_max
     module procedure set_layout_1D_i_max, set_layout_3D_i_max
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
     module procedure get_layout_1D_num_nodes, get_layout_3D_num_nodes
  end interface

  interface get_layout_box
     module procedure get_layout_1D_box, get_layout_3D_box
  end interface

  interface sll_get_num_nodes
     module procedure sll_get_num_nodes_1D, sll_get_num_nodes_3D
  end interface


contains

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

NEW_LAYOUT_FUNCTION( new_layout_1D, layout_1D_t )
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

NEW_DELETE_LAYOUT_FUNCTION( delete_layout_1D, layout_1D_t )
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

  function get_layout_1D_num_nodes( layout )
    sll_int32                  :: get_layout_1D_num_nodes
    type(layout_1D_t), pointer :: layout
    get_layout_1D_num_nodes = sll_get_collective_size( layout%collective )
  end function get_layout_1D_num_nodes

  function get_layout_3D_num_nodes( layout )
    sll_int32                  :: get_layout_3D_num_nodes
    type(layout_3D_t), pointer :: layout
    get_layout_3D_num_nodes = sll_get_collective_size( layout%collective )
  end function get_layout_3D_num_nodes

  function get_layout_1D_box( layout, rank )
    type(box_1D)               :: get_layout_1D_box
    type(layout_1D_t), pointer :: layout
    sll_int32, intent(in)      :: rank
    SLL_ASSERT((rank.ge.0).and.(rank.le.(get_layout_1D_num_nodes(layout)-1)))
    get_layout_1D_box = layout%boxes(rank)
  end function get_layout_1D_box

  function get_layout_3D_box( layout, rank )
    type(box_3D)               :: get_layout_3D_box
    type(layout_3D_t), pointer :: layout
    sll_int32, intent(in)      :: rank
    SLL_ASSERT((rank.ge.0).and.(rank.le.(get_layout_3D_num_nodes(layout)-1)))
    get_layout_3D_box = layout%boxes(rank)
  end function get_layout_3D_box

  ! OK, so I broke the rule of only one % sign... Does this really matter
  ! if this happens inside a library function that the top-level users are 
  ! never supposed to touch (or see)? This is hidden by an access function
  ! anyway. Would it make sense to prohibit the use of nested structures
  ! in a third-party library, like say FFTW or anything else? So why do we 
  ! prohibit this to ourselves? 

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

  ! Use the macros for the 1D case, which will be deleted soon!

  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_1D_i_min, layout_1D_t, i_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_1D_i_max, layout_1D_t, i_max )

  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_1D_i_min, layout_1D_t, i_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_1D_i_max, layout_1D_t, i_max )

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
  function get_layout_1D_collective( layout )
    intrinsic                       :: associated
    type(sll_collective_t), pointer :: get_layout_1D_collective
    type(layout_1D_t), pointer      :: layout
    if( .not. associated(layout) ) then
       stop 'ERROR: uninitialized argument, get_layout_1D_collective()'
    end if
    get_layout_1D_collective => layout%collective
  end function get_layout_1D_collective

  function get_layout_3D_collective( layout )
    intrinsic                       :: associated
    type(sll_collective_t), pointer :: get_layout_3D_collective
    type(layout_3D_t), pointer      :: layout
    if( .not. associated(layout) ) then
       stop 'ERROR: uninitialized argument, get_layout_3D_collective()'
    end if
    get_layout_3D_collective => layout%collective
  end function get_layout_3D_collective

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

  function new_remap_plan_1D( initial, final )
    intrinsic                       :: associated
    type(remap_plan_1D_t), pointer  :: new_remap_plan_1D
    type(layout_1D_t), pointer      :: initial
    type(layout_1D_t), pointer      :: final
    type(sll_collective_t), pointer :: col
    type(box_1D)                    :: ibox, fbox, inters
    sll_int32                       :: i, f          ! loop index
    sll_int32                       :: my_rank
    sll_int32                       :: col_size
    sll_int32                       :: ierr
    sll_int32                       :: disp_counter  ! displacements counter
    sll_int32                       :: send_counter
    sll_int32                       :: recv_counter

    if( (.not. associated(initial)) .or. (.not. associated(final)) ) then
       write (*,'(a)') &
            'ERROR: un-initialized arguments given to sll_new_remap_plan_1D'
       stop 'sll_new_remap_plan_1D'
    end if
    if( .not. associated(get_layout_1D_collective(initial),&
         target=get_layout_1D_collective(final)) ) then
       write (*,'(a)') &
            'ERROR: init and final configurations given to new_remap_plan do not refer to the same collective.'
       stop 'new_remap_plan_1D'
    end if

    col => get_layout_1D_collective(initial)
    my_rank  = sll_get_collective_rank( col )
    col_size = sll_get_collective_size( col )

    SLL_ALLOCATE( new_remap_plan_1D, ierr )
    SLL_CLEAR_ALLOCATE( new_remap_plan_1D%send_displs(0:col_size-1), ierr )
    SLL_CLEAR_ALLOCATE( new_remap_plan_1D%send_cnts(0:col_size-1), ierr )
    SLL_CLEAR_ALLOCATE( new_remap_plan_1D%recv_displs(0:col_size-1), ierr )
    SLL_CLEAR_ALLOCATE( new_remap_plan_1D%recv_cnts(0:col_size-1), ierr )
    ! Can't use CLEAR_ALLOCATE with types for which the '= 0' has not been
    ! defined as an operator.
    SLL_ALLOCATE( new_remap_plan_1D%send_boxes(0:col_size-1), ierr )
    SLL_ALLOCATE( new_remap_plan_1D%recv_boxes(0:col_size-1), ierr )
    new_remap_plan_1D%collective => get_layout_1D_collective(initial)
    send_counter = 0
    disp_counter = 0
    ibox = get_layout_1D_box(initial, my_rank)

    new_remap_plan_1D%initial_layout => initial
    new_remap_plan_1D%final_layout   => final
    ! Find what data to send.
    do f = 0, col_size-1  ! loop over the final layout to look for
                          ! box intersections.
    fbox = get_layout_1D_box(final, f)
       if( intersect_boxes_1D( ibox, fbox, inters ) ) then 
          ! compute how many elements to send
          send_counter                       = get_box_1D_i_max(inters) - &
                                               get_box_1D_i_min(inters) + 1
          new_remap_plan_1D%send_cnts(f)     = send_counter
          new_remap_plan_1D%send_displs(f)   = disp_counter
          disp_counter                       = disp_counter + send_counter
          new_remap_plan_1D%send_boxes(f)    = inters
       else ! if no intersection, there is no need to send anything.
          new_remap_plan_1D%send_cnts(f)     = 0
          new_remap_plan_1D%send_displs(f)   = disp_counter
          new_remap_plan_1D%send_boxes(f)    = inters
       end if
       ! end if
    end do

    disp_counter = 0  
! write (*,'(a, i4)') 'Displacement counter = ', disp_counter
    ! Find what data to receive. Now we compare it with the target layout
    ! for this node.
    fbox = get_layout_1D_box(final, my_rank)
    do i = 0, col_size-1  ! loop over the initial layout to look for box
                          ! intersections.
       ibox = get_layout_1D_box(initial,i)
       if( intersect_boxes_1D( ibox, fbox, inters ) ) then
          ! compute how many elements to receive
          recv_counter =  get_box_1D_i_max(inters) - &
                          get_box_1D_i_min(inters) + 1
          new_remap_plan_1D%recv_cnts(i)   = recv_counter

!write (*,'(a,i4,i4)') 'displacement counter in rank, node: ', my_rank, i
!print *, disp_counter
!call flush()

             new_remap_plan_1D%recv_displs(i) = disp_counter
             disp_counter                     = disp_counter + recv_counter
             new_remap_plan_1D%recv_boxes(i)  = inters
          else ! no intersection, don't expect to receive anything
             new_remap_plan_1D%recv_cnts(i)   = 0
             new_remap_plan_1D%recv_displs(i) = disp_counter
             new_remap_plan_1D%recv_boxes(i)  = inters
          end if
      ! end if
    end do
!write (*,'(a,i4)') ' displacements from rank ', my_rank
!print *, new_remap_plan_1D%recv_displs(:)
!call flush
  end function new_remap_plan_1D

  ! The remap plan stores the buffers where the data to be sent/received
  ! are kept. This raises the issue of type dependence. We want to make this
  ! facility as general as possible. Here we try the approach of having a
  ! single, standard format for data storage, i.e. an integer. This means
  ! that we would need to use the transfer() function to store and retrieve
  ! data from the buffers, which is very inefficient. The alternative is to
  ! have type-dependent plans...

  function new_remap_plan_3D( initial, final, int32_data_size )
    intrinsic                       :: associated, ceiling
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
    SLL_CLEAR_ALLOCATE( new_remap_plan_3D%send_cnts(0:col_size-1), ierr )
    SLL_CLEAR_ALLOCATE( new_remap_plan_3D%recv_displs(0:col_size-1), ierr )
    SLL_CLEAR_ALLOCATE( new_remap_plan_3D%recv_cnts(0:col_size-1), ierr )
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
          new_remap_plan_3D%send_cnts(f)   = send_counter
          new_remap_plan_3D%send_displs(f) = disp_counter
          disp_counter                     = disp_counter + send_counter
          new_remap_plan_3D%send_boxes(f)  = inters
          acc                              = acc + send_counter
       else ! if no intersection, there is no need to send anything.
          new_remap_plan_3D%send_cnts(f)     = 0
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
          new_remap_plan_3D%recv_cnts(i)   = recv_counter
          new_remap_plan_3D%recv_displs(i) = disp_counter
          disp_counter                     = disp_counter + recv_counter
          new_remap_plan_3D%recv_boxes(i)  = inters
          acc                              = acc + recv_counter
       else ! no intersection, don't expect to receive anything
          new_remap_plan_3D%recv_cnts(i)   = 0
          new_remap_plan_3D%recv_displs(i) = disp_counter
          new_remap_plan_3D%recv_boxes(i)  = inters
       end if
    end do
    SLL_ALLOCATE(new_remap_plan_3D%recv_buffer(0:(acc*int32_data_size-1)),ierr)
  end function new_remap_plan_3D

  function get_remap_1D_initial_layout( plan )
    type(layout_1D_t), pointer      :: get_remap_1D_initial_layout
    type(remap_plan_1D_t), pointer  :: plan
    if( .not. associated( plan ) ) then
       write (*,'(a)') 'not associated pointer argument'
       stop 'get_remap_1D_initial_layout'
    end if
    get_remap_1D_initial_layout => plan%initial_layout
  end function get_remap_1D_initial_layout

  function get_remap_1D_final_layout( plan )
    type(layout_1D_t), pointer      :: get_remap_1D_final_layout
    type(remap_plan_1D_t), pointer  :: plan
    ! FIXME: arg checking
    get_remap_1D_final_layout => plan%final_layout
  end function get_remap_1D_final_layout

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
  ! The idea is to always communicate the MPI_BYTE datatype. The only extra
  ! step is to translate the arrays for counts and displacements in terms
  ! of their byte-sizes. While we would still need a different 
  ! apply_remap_1D() function for every new type, at least the modularization
  ! will be maintained.
  !
  ! For apply_remap_1D_int(), we use this approach as a test case, even though
  ! it is not necessary since 'integer' is also a native MPI type.

  subroutine convert_to_bytes( sz, ai, n, bi )
    sll_int32, intent(in)                :: sz   ! in bytes
    sll_int32, intent(in), dimension(:)  :: ai   ! array to convert
    sll_int32, intent(in)                :: n    ! size of array
    sll_int32, intent(out), dimension(:) :: bi   ! output
    sll_int32                            :: i
    do i=1,n
       bi(i) = ai(i)*sz
    end do
  end subroutine convert_to_bytes

  subroutine apply_remap_1D_int( plan, data_in, data_out )
    type(remap_plan_1D_t), pointer       :: plan
    sll_int32, dimension(:), intent(in)  :: data_in
    sll_int32, dimension(:), intent(out) :: data_out
    sll_int32, dimension(:), allocatable, target :: sb       ! send buffer
    sll_int32, dimension(:), pointer     :: sptr
    sll_int32, dimension(:), allocatable :: rb       ! receive buffer
    sll_int32, dimension(:), pointer     :: sdisp    ! send displacements
    sll_int32, dimension(:), pointer     :: rdisp    ! receive displacements
    sll_int32, dimension(:), pointer     :: scnts    ! send counts
    sll_int32, dimension(:), pointer     :: rcnts    ! receive counts
    type(sll_collective_t), pointer      :: col      ! collective
    type(layout_1D_t), pointer           :: init_layout => NULL()
    sll_int32                            :: i
    sll_int32                            :: col_sz
    sll_int64                            :: sacc = 0 ! send count accumulator
    sll_int64                            :: racc = 0 ! receive count accum.
    sll_int64                            :: send_sz  ! size of send buffer
    sll_int64                            :: recv_sz  ! size of recv buffer
    sll_int32                            :: ierr
    sll_int32                            :: lo
    sll_int32                            :: hi
    type(box_1D)                         :: sbox
    sll_int32                            :: j
    sll_int32                            :: my_rank
    sll_int32                            :: loc
    sll_int32                            :: num
    ! to load the MPI function and send bytes, we have a separate set of
    ! arrays to store this information for now.
    sll_int32, dimension(:), allocatable     :: sdispb  ! send displacements
    sll_int32, dimension(:), allocatable     :: rdispb  ! receive displacements
    sll_int32, dimension(:), allocatable     :: scntsb  ! send counts
    sll_int32, dimension(:), allocatable     :: rcntsb  ! receive counts

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp       => plan%send_displs
    rdisp       => plan%recv_displs
    scnts       => plan%send_cnts
    rcnts       => plan%recv_cnts
    col         => plan%collective
    init_layout => get_remap_1D_initial_layout(plan)
    my_rank     =  sll_get_collective_rank(col)

    ! estimate the size of the send and receive buffers and allocate them
    col_sz = sll_get_collective_size(col)
    do i=0,col_sz-1
       sacc = sacc + scnts(i)
       racc = racc + rcnts(i)
    end do

    ! Send/Receive buffers are zero-indexed
    SLL_ALLOCATE(sb(0:sacc-1),ierr)
    SLL_ALLOCATE(rb(0:racc-1),ierr)

    SLL_ALLOCATE(sdispb(0:col_sz-1), ierr)
    SLL_ALLOCATE(rdispb(0:col_sz-1), ierr)
    SLL_ALLOCATE(scntsb(0:col_sz-1), ierr)
    SLL_ALLOCATE(rcntsb(0:col_sz-1), ierr)

    ! Translate the amounts into bytes
#if 1
    call convert_to_bytes(BYTE_SIZEOF(data_in(1)), sdisp, col_sz, sdispb)
    call convert_to_bytes(BYTE_SIZEOF(data_in(1)), rdisp, col_sz, rdispb)
    call convert_to_bytes(BYTE_SIZEOF(data_in(1)), scnts, col_sz, scntsb)
    call convert_to_bytes(BYTE_SIZEOF(data_in(1)), rcnts, col_sz, rcntsb)
#endif
write (*,'(a,i4)') 'parameters from rank ', my_rank
print *, scntsb(:)
print *, sdispb(:)
print *, rcntsb(:)
print *, rdispb(:)



    ! load the send buffer
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          ! get the information on the box to send, get the limits and
          ! find out where in the buffer to start writing.
          sbox = plan%send_boxes(i)
          lo   = global_to_local_1D( init_layout, get_box_1D_i_min(sbox) )
          hi   = global_to_local_1D( init_layout, get_box_1D_i_max(sbox) )
          num  = scnts(i)
          loc  = sdisp(i)
          sptr => sb(loc:)
          sptr(1:num) = data_in(lo:hi)
       end if
    end do

       write (*,'(a,i4)') 'the send buffer in rank:', my_rank
       print *, sb(0:(size(sb)-1))
       call flush()

       call sll_collective_all_to_allV_int( sb(:),       &
                                            scntsb(0:col_sz-1), &
                                            sdispb(0:col_sz-1), &
                                            rb(:),       &
                                            rcntsb(0:col_sz-1), &
                                            rdispb(0:col_sz-1), col )
 write (*,'(a, i4)') 'the receive buffer in rank: ', my_rank
       print *, rb(0:size(rb)-1)
       call flush()

       ! deallocate the buffers
       SLL_DEALLOCATE_ARRAY(sb, ierr)
       SLL_DEALLOCATE_ARRAY(rb, ierr)
       
     end subroutine apply_remap_1D_int

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
    type(layout_3D_t), pointer               :: init_layout => NULL()
    sll_int32                                :: id, jd, kd
    sll_int32                                :: i
    sll_int32                                :: col_sz
    sll_int64                                :: sacc = 0 ! send accumulator
    sll_int64                                :: racc = 0 ! receive count accum.
    sll_int64                                :: send_sz  ! size of send buffer
    sll_int64                                :: recv_sz  ! size of recv buffer
    sll_int32                                :: ierr
    sll_int32                                :: loi, loj, lok
    sll_int32                                :: hii, hij, hik
    type(box_3D)                             :: sbox
    sll_int32                                :: my_rank
    sll_int32                                :: loc
    sll_int32                                :: num
    sll_int32, dimension(1:3)                :: local_lo, local_hi
    sll_int32                                :: int32_data_size

    ! to load the MPI function and send bytes, we have a separate set of
    ! arrays to store this information for now.
    sll_int32, dimension(:), allocatable     :: sdispb  ! send displacements
    sll_int32, dimension(:), allocatable     :: rdispb  ! receive displacements
    sll_int32, dimension(:), allocatable     :: scntsb  ! send counts
    sll_int32, dimension(:), allocatable     :: rcntsb  ! receive counts

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp       => plan%send_displs
    rdisp       => plan%recv_displs
    scnts       => plan%send_cnts
    rcnts       => plan%recv_cnts
    col         => plan%collective
    init_layout => get_remap_3D_initial_layout(plan)
    my_rank     =  sll_get_collective_rank(col)
    sb          => plan%send_buffer
    rb          => plan%recv_buffer
    
    SLL_ALLOCATE(sdispb(0:col_sz-1), ierr)
    SLL_ALLOCATE(rdispb(0:col_sz-1), ierr)
    SLL_ALLOCATE(scntsb(0:col_sz-1), ierr)
    SLL_ALLOCATE(rcntsb(0:col_sz-1), ierr)

    ! Translate the amounts into bytes
#if 1
    call convert_to_bytes(INT32_SIZEOF(data_in(1,1,1)), sdisp, col_sz, sdispb)
    call convert_to_bytes(INT32_SIZEOF(data_in(1,1,1)), rdisp, col_sz, rdispb)
    call convert_to_bytes(INT32_SIZEOF(data_in(1,1,1)), scnts, col_sz, scntsb)
    call convert_to_bytes(INT32_SIZEOF(data_in(1,1,1)), rcnts, col_sz, rcntsb)
#endif

#if 1
    write (*,'(a,i4)') 'parameters from rank ', my_rank
    print *, scntsb(:)
    print *, sdispb(:)
    print *, rcntsb(:)
    print *, rdispb(:)
#endif

    ! load the send buffer
    loc = 0             ! first loading is at position zero
    ! This step is obviously not needed for integers themselves. We put this
    ! here for generality.
    int32_data_size = INT32_SIZEOF( data_in(1,1,1) )
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdispb(i) ) then
             write (*,'(a,i4)') 'ERROR: discrepancy between displs(i) and the loading index for i = ', i
             stop 'loading error'
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
    
    write (*,'(a,i4)') 'the send buffer in rank:', my_rank
    print *, sb(0:(size(sb)-1))
    call flush()

    call sll_collective_all_to_allV_int( sb(:),       &
                                         scntsb(0:col_sz-1), &
                                         sdispb(0:col_sz-1), &
                                         rb(:),       &
                                         rcntsb(0:col_sz-1), &
                                         rdispb(0:col_sz-1), col )
    write (*,'(a, i4)') 'the receive buffer in rank: ', my_rank
    print *, rb(0:size(rb)-1)
    call flush()
    
  end subroutine apply_remap_3D_int



#define MAKE_GET_BOX_SLOT_FUNCTION( fname, boxtype, slot )   \
  function fname( b );                                       \
    sll_int32                 :: fname;                      \
    type(boxtype), intent(in) :: b;                          \
    fname = b%slot;                                          \
  end function fname

  MAKE_GET_BOX_SLOT_FUNCTION( get_box_1D_i_min, box_1D, i_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_1D_i_max, box_1D, i_max )

  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_i_min, box_3D, i_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_i_max, box_3D, i_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_j_min, box_3D, j_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_j_max, box_3D, j_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_k_min, box_3D, k_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_k_max, box_3D, k_max )

  function sll_get_num_nodes_1D( lims )
    sll_int32                  :: sll_get_num_nodes_1D
    type(layout_1D_t), pointer :: lims
    sll_get_num_nodes_1D = sll_get_collective_size( lims%collective )
  end function sll_get_num_nodes_1D

  function sll_get_num_nodes_3D( lims )
    sll_int32                  :: sll_get_num_nodes_3D
    type(layout_3D_t), pointer :: lims
    sll_get_num_nodes_3D = sll_get_collective_size( lims%collective )
  end function sll_get_num_nodes_3D




  ! uniform_dist_1D() gives a two-element array, the first element is the
  ! lowest index stored at the given node, and the second element is the
  ! highest. For now this is to be considered a 'helper' function and is
  ! thus not really a full member of this module. 
  function uniform_dist_1D( proc, nelem, group_sz )
    sll_int32, dimension(1:2) :: uniform_dist_1D
    sll_int32, intent(in)     :: proc     ! process number
    sll_int32, intent(in)     :: nelem    ! total number of elements
    sll_int32, intent(in)     :: group_sz ! number of processes in group    
    sll_int32                 :: nloc     ! local number of elements/process
    nloc               = nelem/group_sz
    uniform_dist_1D(1) = proc*nloc + 1
    uniform_dist_1D(2) = (proc+1)*nloc
  end function uniform_dist_1D



  ! TENTATIVE... another helper function.
  ! ugly implementation of an indexing which allows negative indices as well
  ! as indices greater than the array size, always returning an index within
  ! the array as if it were repeated indefinitely in both directions.
  function circular_array_index( i, N )
    intrinsic           :: modulo, abs
    integer             :: circular_array_index
    integer, intent(in) :: i
    integer, intent(in) :: N ! array size
    if( (i .ge. 1) .and. (i .le. N) ) then ! normal, in-bounds index
       circular_array_index = i
    else if ( i .gt. N ) then
       circular_array_index = modulo(i,N) ! only works for i<2N
    else if ( i .eq. 0 ) then
       print *, '0 index given to circular array'
       stop 'circular array: zero index'
    else
       circular_array_index = N - abs(i) + 1 ! only works for i > -2N
    end if
  end function circular_array_index

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
  function local_to_global_1D( layout, i_local )
    sll_int32                       :: local_to_global_1D
    type(layout_1D_t), pointer      :: layout
    sll_int32, intent(in)           :: i_local
    type(sll_collective_t), pointer :: col
    sll_int32                       :: my_rank
    type(box_1D)                    :: box
    ! FIXME: arg checking
print *, 'entered here?'
call flush()
    col                => get_layout_1D_collective( layout )
    my_rank            =  sll_get_collective_rank( col )
    box                =  get_layout_1D_box( layout, my_rank )
    local_to_global_1D =  get_box_1D_i_min(box) + i_local - 1
  end function local_to_global_1D

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
  function global_to_local_1D( layout, i_global )
    intrinsic                       :: associated
    sll_int32                       :: global_to_local_1D
    type(layout_1D_t), pointer      :: layout
    sll_int32, intent(in)           :: i_global
    type(sll_collective_t), pointer :: col
    sll_int32                       :: my_rank
    type(box_1D)                    :: box
    ! FIXME: some argument checking would be nice...
    if( .not. associated(get_layout_1D_collective(layout)) ) then
       write (*,'(a)') 'ERROR in global_to_local_1D(), not-associated col'
       stop 'global_to_local_1D'
    end if

    col     => get_layout_1D_collective( layout )
    my_rank =  sll_get_collective_rank( col )
    box     =  get_layout_1D_box( layout, my_rank )
    if( (i_global .ge. get_box_1D_i_min(box)) .and. &
        (i_global .le. get_box_1D_i_max(box)) ) then  ! the index is present
       global_to_local_1D = i_global - get_box_1D_i_min(box) + 1
    else  ! the index is not present
       global_to_local_1D = -1
    end if
  end function global_to_local_1D

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

  ! compute_min_max_1D returns the a 2-element array, with the lower
  ! and upper indices of the data, given a process and

  function compute_min_max_1D( rank, total_size, col_size )
    sll_int32, dimension(2) :: compute_min_max_1D
    sll_int32, intent(in)   :: rank
    sll_int32, intent(in)   :: total_size
    sll_int32, intent(in)   :: col_size
    sll_int32               :: local_size
    sll_int32               :: lo
    sll_int32               :: hi
    local_size      = total_size/col_size
    lo              = rank*local_size + 1 
    hi              = local_size*(rank + 1)
#if 1
    if( (lo .lt. 1) .or. (lo .gt. total_size) ) then 
       lo = circular_array_index( lo, total_size )
    end if
    if( (hi .lt. 1) .or. (hi .gt. total_size) ) then 
       hi = circular_array_index( hi, 16 )
    end if
#endif
    compute_min_max_1D(1) = lo
    compute_min_max_1D(2) = hi
  end function compute_min_max_1D
  
  function compute_min_max2( rank, total_size, col_size )
    sll_int32, dimension(2) :: compute_min_max2
    sll_int32, intent(in)   :: rank
    sll_int32, intent(in)   :: total_size
    sll_int32, intent(in)   :: col_size
    sll_int32               :: local_size
    sll_int32               :: lo
    sll_int32               :: hi

    local_size      = total_size/col_size
    lo              = rank*local_size + 1 
    hi              = local_size*(rank + 1)
    if( rank .eq. 0 ) then
       compute_min_max2(1) = lo
       compute_min_max2(2) = hi + 1
    else if( rank .eq. col_size - 1 ) then
       compute_min_max2(1) = lo - 1
       compute_min_max2(2) = hi 
    else
#if 0
    if( (lo .lt. 1) .or. (lo .gt. total_size) ) then 
       lo = circular_array_index( lo, total_size )
    end if
    if( (hi .lt. 1) .or. (hi .gt. total_size) ) then 
       hi = circular_array_index( hi, 16 )
    end if
#endif
    compute_min_max2(1) = lo - 1
    compute_min_max2(2) = hi + 1
 end if
end function compute_min_max2

function compute_min_max3( rank, total_size, col_size )
    sll_int32, dimension(2) :: compute_min_max3
    sll_int32, intent(in)   :: rank
    sll_int32, intent(in)   :: total_size
    sll_int32, intent(in)   :: col_size
    compute_min_max3(1) = 1
    compute_min_max3(2) = 16
  end function compute_min_max3

  subroutine initialize_layout_1D( lims_func, total_sz, local_sz, layout )

    interface

       function lims_func( rank, total_size, col_size )
         ! how to bring something like integer(kind=i32) back into scope here???
         integer, dimension(2) :: lims_func
         integer, intent(in)   :: rank
         integer, intent(in)   :: total_size
         integer, intent(in)   :: col_size
       end function lims_func
    end interface

    type(layout_1D_t), pointer :: layout
    sll_int32, intent(in)      :: total_sz
    sll_int32, intent(in)      :: local_sz
    sll_int32                  :: col_sz ! size of the collective
    sll_int32                  :: i
    sll_int32                  :: lo
    sll_int32                  :: hi
    sll_int32, dimension(1:2)  :: lh

    col_sz  = sll_get_num_nodes(layout)
    do i=0,col_sz-1
       lh = lims_func( i, total_sz, col_sz )
       lo = lh(1)
       hi = lh(2)
       call set_layout_1D_i_min( layout, i, lo )
       call set_layout_1D_i_max( layout, i, hi )
    end do
   end subroutine initialize_layout_1D



  subroutine sll_view_lims_1D( lims )
    type(layout_1D_t), pointer :: lims
    sll_int32                  :: i
    sll_int32                  :: sz
    sz = sll_get_num_nodes( lims )
    print *, 'limits:'
    do i=0,sz-1
       write (*, '(i4, i4)') get_layout_1D_i_min( lims, i ), &
                             get_layout_1D_i_max( lims, i )
    end do
    call flush()
  end subroutine sll_view_lims_1D

  subroutine sll_view_lims_3D( lims )
    type(layout_3D_t), pointer :: lims
    sll_int32                  :: i
    sll_int32                  :: sz
    sz = sll_get_num_nodes( lims )
    print *, 'limits: '
    do i=0,sz-1
       write(*,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') &
            '[[', get_layout_i_min(lims,i),',',&
                  get_layout_i_max(lims,i),'], [', &
                  get_layout_j_min(lims,i),',',    &
                  get_layout_j_max(lims,i),'], [', &
                  get_layout_k_min(lims,i),',',    &
                  get_layout_k_max(lims,i),']]'
    end do
    call flush()
  end subroutine sll_view_lims_3D

#if 0
  function box_is_ordered( b )
    logical
#endif


  ! the return value of intersect_boxes_1D() is 'logical' valued, and answers
  ! the question whether the boxes intersect or not. 'ans' is a box with
  ! the actual intersection between the argument boxes. In case that there
  ! is no intersection between the boxes the value [0,0] is returned. 
  !
  ! Here we introduce a complication, which is the possibility that a 'box'
  ! wraps around the minimum and maximum indices of an array. Thus, for
  ! an array 'a' indexed from 1 through N, we would be allowing an 
  ! interval such as [a(N-2),a(3)]. 

 function intersect_boxes_1D( b1, b2, ans )
    intrinsic                    :: min, max
    logical                      :: intersect_boxes_1D
    type(box_1D), intent(in)     :: b1, b2
    type(box_1D), intent(out)    :: ans
    sll_int32                    :: lo, hi
    sll_int32                    :: lob1, hib1, lob2, hib2

    lob1      = get_box_1D_i_min(b1) 
    hib1      = get_box_1D_i_max(b1)
    lob2      = get_box_1D_i_min(b2)
    hib2      = get_box_1D_i_max(b2)

    SLL_ASSERT( (lob1 .le. hib1) .and. (lob2 .le. hib2) )

    lo = max(lob1, lob2)
    hi = min(hib1, hib2)
    if(lo .gt. hi) then ! there is no intersection
       ans%i_min = 0
       ans%i_max = 0
       intersect_boxes_1D = .false.
    else
       ans%i_min = lo
       ans%i_max = hi
       intersect_boxes_1D = .true.
    end if
  end function intersect_boxes_1D

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


  function boxes_equal_1D( b1, b2 )
    logical :: boxes_equal_1D
    type(box_1D), intent(in) :: b1, b2
    if( (b1%i_min .eq. b2%i_min) .and. (b1%i_max .eq. b1%i_max) ) then
       boxes_equal_1D = .true.
    else
       boxes_equal_1D = .false.
    end if
  end function boxes_equal_1D

  subroutine view_box_1D( b )
    type(box_1D), intent(in) :: b
    write (*, '(a, i8, a, i8, a)') 'box limits = [', b%i_min, ', ', b%i_max, ']'
  end subroutine view_box_1D

  subroutine view_remap_plan_1D( plan )
    type(remap_plan_1D_t), pointer  :: plan
    sll_int32                       :: rank
    sll_int32                       :: i
    sll_int32                       :: col_size
    rank     = sll_get_collective_rank(plan%collective)
    col_size = sll_get_collective_size(plan%collective)

    write (*,'(a i8)') 'From rank: ', rank
    print *, plan%send_cnts(:)
    print *, plan%send_displs(:)
    print *, plan%recv_cnts(:)
    print *, plan%recv_displs(:)
    print *, 'boxes to send: '
    do i=0,col_size-1
       call view_box_1D(plan%send_boxes(i))
    end do
    print *, 'boxes to receive: '
    do i=0,col_size-1
       call view_box_1D(plan%recv_boxes(i))
    end do
  end subroutine view_remap_plan_1D


end module remapper
