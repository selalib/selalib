!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

module sll_remapper
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"
  use sll_collective
  use sll_electric_field_2d_accumulator ! terrible dependency here...
  implicit none
  
  ! The box types contain information on the index limits contained        
  ! in a given processor.
  type box_2D
     sll_int32 :: i_min, i_max
     sll_int32 :: j_min, j_max
  end type box_2D

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

  type box_6D
     sll_int32 :: i_min, i_max
     sll_int32 :: j_min, j_max
     sll_int32 :: k_min, k_max
     sll_int32 :: l_min, l_max
     sll_int32 :: m_min, m_max
     sll_int32 :: n_min, n_max
  end type box_6D

  
  ! The sll_layout types contain information on a collective and an
  ! array of boxes that describes the distribution of data among
  ! different nodes.
  type layout_2D
     type(sll_collective_t), pointer     :: collective
     type(box_2D), dimension(:), pointer :: boxes
  end type layout_2D

  type layout_3D
     type(sll_collective_t), pointer     :: collective
     type(box_3D), dimension(:), pointer :: boxes
  end type layout_3D

  type layout_4D
     type(sll_collective_t), pointer     :: collective
     type(box_4D), dimension(:), pointer :: boxes
  end type layout_4D

  type layout_5D
     type(sll_collective_t), pointer     :: collective
     type(box_5D), dimension(:), pointer :: boxes
  end type layout_5D

  type layout_6D
     type(sll_collective_t), pointer     :: collective
     type(box_6D), dimension(:), pointer :: boxes
  end type layout_6D


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
#define MAKE_REMAP_PLAN( type_name, layout_type, box_type, data_type )   \
  type type_name;                                             \
     type(layout_type), pointer            :: initial_layout; \
     type(layout_type), pointer            :: final_layout;   \
     integer, dimension(:), pointer        :: send_displs;    \
     integer, dimension(:), pointer        :: send_counts;    \
     integer, dimension(:), pointer        :: recv_displs;    \
     integer, dimension(:), pointer        :: recv_counts;    \
     type(box_type), dimension(:), pointer :: send_boxes;     \
     type(box_type), dimension(:), pointer :: recv_boxes;     \
     type(sll_collective_t), pointer       :: collective;     \
     data_type, dimension(:), pointer      :: send_buffer;    \
     data_type, dimension(:), pointer      :: recv_buffer;    \
     logical                               :: is_uniform;     \
  end type type_name

  ! 2D Remap types:
  MAKE_REMAP_PLAN(remap_plan_2D_int32, layout_2D, box_2D, sll_int32)
  MAKE_REMAP_PLAN(remap_plan_2D_real64, layout_2D, box_2D, sll_real64)
  MAKE_REMAP_PLAN(remap_plan_2D_comp64, layout_2D, box_2D, sll_comp64)
  ! 3D Remap types:
  MAKE_REMAP_PLAN(remap_plan_3D_int32, layout_3D, box_3D, sll_int32)
  MAKE_REMAP_PLAN(remap_plan_3D_real64, layout_3D, box_3D, sll_real64)
  MAKE_REMAP_PLAN(remap_plan_3D_comp64, layout_3D, box_3D, sll_comp64)
  ! 4D Remap types:
  MAKE_REMAP_PLAN(remap_plan_4D_int32, layout_4D, box_4D, sll_int32)
  MAKE_REMAP_PLAN(remap_plan_4D_real64, layout_4D, box_4D, sll_real64)
  MAKE_REMAP_PLAN(remap_plan_4D_comp64, layout_4D, box_4D, sll_comp64)
  ! 6D Remap types:
  MAKE_REMAP_PLAN(remap_plan_6D_int32, layout_6D, box_6D, sll_int32)
  MAKE_REMAP_PLAN(remap_plan_6D_real64, layout_6D, box_6D, sll_real64)
  MAKE_REMAP_PLAN(remap_plan_6D_comp64, layout_6D, box_6D, sll_comp64)

 interface get_layout_i_min
     module procedure get_layout_2D_i_min, get_layout_3D_i_min, &
          get_layout_4D_i_min, get_layout_6d_i_min
  end interface

  interface set_layout_i_min
     module procedure set_layout_2D_i_min, set_layout_3D_i_min, &
          set_layout_4D_i_min, set_layout_6D_i_min
  end interface

  interface get_layout_i_max
     module procedure get_layout_2D_i_max, get_layout_3D_i_max, &
          get_layout_4D_i_max, get_layout_6D_i_max
  end interface

  interface set_layout_i_max
     module procedure set_layout_2D_i_max, set_layout_3D_i_max, &
          set_layout_4D_i_max, set_layout_6D_i_max
  end interface

  interface get_layout_j_min
     module procedure get_layout_2D_j_min, get_layout_3D_j_min, &
          get_layout_4D_j_min, get_layout_6D_j_min
  end interface

  interface set_layout_j_min
     module procedure set_layout_2D_j_min, set_layout_3D_j_min, &
          set_layout_4D_j_min, set_layout_6D_j_min  
  end interface

  interface get_layout_j_max
     module procedure get_layout_2D_j_max, get_layout_3D_j_max, &
          get_layout_4D_j_max, get_layout_6D_j_max
  end interface

  interface set_layout_j_max
     module procedure set_layout_2D_j_max, set_layout_3D_j_max, &
          set_layout_4D_j_max, set_layout_6D_j_max
  end interface

 interface get_layout_k_min
    module procedure get_layout_3D_k_min, get_layout_4D_k_min, &
         get_layout_6D_k_min
  end interface

  interface set_layout_k_min
     module procedure set_layout_3D_k_min, set_layout_4D_k_min, &
          set_layout_6D_k_min
  end interface

  interface get_layout_k_max
     module procedure get_layout_3D_k_max, get_layout_4D_k_max, &
          get_layout_6D_k_max
  end interface

  interface set_layout_k_max
     module procedure set_layout_3D_k_max, set_layout_4D_k_max, &
          set_layout_6D_k_max
  end interface

  interface get_layout_l_min
     module procedure get_layout_4D_l_min, get_layout_6D_l_min
  end interface get_layout_l_min

  interface set_layout_l_min
     module procedure set_layout_4D_l_min, set_layout_6D_l_min
  end interface set_layout_l_min

  interface get_layout_l_max
     module procedure get_layout_4D_l_max, get_layout_6D_l_max 
  end interface get_layout_l_max

  interface set_layout_l_max
     module procedure set_layout_4D_l_max, set_layout_6D_l_max
  end interface set_layout_l_max

  interface get_layout_m_min
     module procedure get_layout_6D_m_min
  end interface get_layout_m_min

  interface get_layout_m_max
     module procedure get_layout_6D_m_max
  end interface get_layout_m_max

  interface set_layout_m_min
     module procedure set_layout_6D_m_min
  end interface set_layout_m_min

  interface set_layout_m_max
     module procedure set_layout_6D_m_max
  end interface set_layout_m_max

  interface get_layout_n_min
     module procedure get_layout_6D_n_min
  end interface get_layout_n_min

  interface get_layout_n_max
     module procedure get_layout_6D_n_max
  end interface get_layout_n_max

  interface set_layout_n_min
     module procedure set_layout_6D_n_min
  end interface set_layout_n_min

  interface set_layout_n_max
     module procedure set_layout_6D_n_max
  end interface set_layout_n_max

  interface get_layout_num_nodes
     module procedure get_layout_2D_num_nodes, get_layout_3D_num_nodes, &
          get_layout_4D_num_nodes, get_layout_6D_num_nodes
  end interface

  interface get_layout_box
     module procedure get_layout_2D_box, get_layout_3D_box, get_layout_4D_box, &
          get_layout_6D_box
  end interface

  interface get_layout_collective
     module procedure get_layout_2D_collective, get_layout_3D_collective, &
          get_layout_4D_collective, get_layout_6D_collective
  end interface get_layout_collective

  interface sll_get_num_nodes
     module procedure sll_get_num_nodes_2D, sll_get_num_nodes_3D, &
          sll_get_num_nodes_4D, sll_get_num_nodes_6D
   end interface

  interface count_elements_in_box
     module procedure count_elements_in_box_2D, count_elements_in_box_3D, &
          count_elements_in_box_4D, count_elements_in_box_6D
  end interface count_elements_in_box

  interface intersect_boxes
     module procedure intersect_boxes_2D, intersect_boxes_3D, &
          intersect_boxes_4D, intersect_boxes_6D
  end interface intersect_boxes

  interface optimize_remap_plan
     module procedure &
          optimize_remap_plan_2D_int32, &
          optimize_remap_plan_2D_real64, &
          optimize_remap_plan_2D_comp64, &
          optimize_remap_plan_3D_int32, &
          optimize_remap_plan_3D_real64, &
          optimize_remap_plan_3D_comp64, &
          optimize_remap_plan_4D_int32, &
          optimize_remap_plan_4D_real64, &
          optimize_remap_plan_4D_comp64, &
          optimize_remap_plan_6D_int32, &
          optimize_remap_plan_6D_real64, &
          optimize_remap_plan_6D_comp64
  end interface optimize_remap_plan

  interface new_remap_plan
     module procedure &
          new_remap_plan_2d_int32, &
          new_remap_plan_2d_real64, &
          new_remap_plan_2d_comp64, &
          new_remap_plan_3d_int32, &
          new_remap_plan_3d_real64, &
          new_remap_plan_3d_comp64, &
          new_remap_plan_4d_int32, &
          new_remap_plan_4d_real64, &
          new_remap_plan_4d_comp64, &
          new_remap_plan_6d_int32, &
          new_remap_plan_6d_real64, &
          new_remap_plan_6d_comp64
  end interface new_remap_plan

  interface apply_remap_2D
     module procedure apply_remap_2D_double, apply_remap_2d_complex !, &
!          apply_remap_2d_efield
  end interface apply_remap_2D

  interface apply_remap_3D
     module procedure apply_remap_3D_int, apply_remap_3D_double, &
          apply_remap_3D_complex
  end interface

  interface apply_remap_4D
     module procedure apply_remap_4D_double
  end interface apply_remap_4D

  interface apply_remap_6D
     module procedure apply_remap_6D_double, apply_remap_6D_int
  end interface apply_remap_6D

  interface delete
     module procedure delete_layout_2D, delete_layout_3D, delete_layout_4D, &
          delete_layout_5D, delete_layout_6D, &
          delete_remap_2D_int32, &
          delete_remap_2D_real64, &
          delete_remap_2D_comp64, &
          delete_remap_3D_int32, &
          delete_remap_3D_real64, &
          delete_remap_3D_comp64, &
          delete_remap_4D_int32, &
          delete_remap_4D_real64, &
          delete_remap_4D_comp64, &
          delete_remap_6D_int32, &
          delete_remap_6D_real64, &
          delete_remap_6D_comp64
  end interface delete

  interface compute_local_sizes
     module procedure compute_local_sizes_2d, compute_local_sizes_3d, &
          compute_local_sizes_4d, compute_local_sizes_6d
  end interface compute_local_sizes

  interface local_to_global
     module procedure local_to_global_2D,local_to_global_3D, &
          local_to_global_4D, local_to_global_6D
  end interface local_to_global


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

  NEW_LAYOUT_FUNCTION( new_layout_2D, layout_2D )
  NEW_LAYOUT_FUNCTION( new_layout_3D, layout_3D )
  NEW_LAYOUT_FUNCTION( new_layout_4D, layout_4D )
  NEW_LAYOUT_FUNCTION( new_layout_5D, layout_5D )
  NEW_LAYOUT_FUNCTION( new_layout_6D, layout_6D )

#define NEW_DELETE_LAYOUT_FUNCTION( fname, layout_type )        \
  subroutine fname( layout );                                   \
    type(layout_type), pointer :: layout;                       \
    sll_int32                  :: ierr;                         \
    nullify( layout%collective );                               \
    SLL_DEALLOCATE( layout%boxes, ierr );                       \
    SLL_DEALLOCATE( layout, ierr );                             \
  end subroutine fname

  NEW_DELETE_LAYOUT_FUNCTION( delete_layout_2D, layout_2D )
  NEW_DELETE_LAYOUT_FUNCTION( delete_layout_3D, layout_3D )
  NEW_DELETE_LAYOUT_FUNCTION( delete_layout_4D, layout_4D )
  NEW_DELETE_LAYOUT_FUNCTION( delete_layout_5D, layout_5D )
  NEW_DELETE_LAYOUT_FUNCTION( delete_layout_6D, layout_6D )

  ! Access functions for the boxes. This is really an overkill... On one hand,
  ! it is nice to hide everything behind access functions so that we 
  ! preserve the freedom of changing the representation of the types if
  ! needed. Also, this is not a performance-critical process. On the other
  ! hand, the only thing there is to hide here is a pair of chained %'s...
  ! All these could be greatly reduced by a few one-line macros, but at least
  ! for now we choose the conventional approach.
#define MAKE_NUM_NODES_FUNCTION( fname, layout_type ) \
  function fname( layout ); \
    sll_int32                  :: fname; \
    type(layout_type), pointer :: layout; \
    fname = sll_get_collective_size( layout%collective ); \
  end function fname

  MAKE_NUM_NODES_FUNCTION( get_layout_2D_num_nodes, layout_2D )
  MAKE_NUM_NODES_FUNCTION( get_layout_3D_num_nodes, layout_3D )
  MAKE_NUM_NODES_FUNCTION( get_layout_4D_num_nodes, layout_4D )
  MAKE_NUM_NODES_FUNCTION( get_layout_6D_num_nodes, layout_6D )

#define MAKE_GET_BOX_FUNCTION( fname, layout_type, box_type ) \
  function fname( layout, rank ); \
    type(box_type)             :: fname; \
    type(layout_type), pointer :: layout; \
    sll_int32, intent(in)      :: rank; \
    SLL_ASSERT((rank.ge.0).and.(rank.le.(get_layout_num_nodes(layout)-1))); \
    fname = layout%boxes(rank); \
  end function fname

  MAKE_GET_BOX_FUNCTION( get_layout_2D_box, layout_2D, box_2D )
  MAKE_GET_BOX_FUNCTION( get_layout_3D_box, layout_3D, box_3D )
  MAKE_GET_BOX_FUNCTION( get_layout_4D_box, layout_4D, box_4D )
  MAKE_GET_BOX_FUNCTION( get_layout_6D_box, layout_6D, box_6D )

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

  ! We use the macros to write the set_ get_ functions for the different 
  ! dimensions.

  ! 2D case:
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_2D_i_min, layout_2D, i_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_2D_i_max, layout_2D, i_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_2D_j_min, layout_2D, j_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_2D_j_max, layout_2D, j_max )

  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_2D_i_min, layout_2D, i_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_2D_i_max, layout_2D, i_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_2D_j_min, layout_2D, j_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_2D_j_max, layout_2D, j_max )

  ! 3D case:
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_i_min, layout_3D, i_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_i_max, layout_3D, i_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_j_min, layout_3D, j_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_j_max, layout_3D, j_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_k_min, layout_3D, k_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_3D_k_max, layout_3D, k_max )

  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_i_min, layout_3D, i_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_i_max, layout_3D, i_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_j_min, layout_3D, j_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_j_max, layout_3D, j_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_k_min, layout_3D, k_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_3D_k_max, layout_3D, k_max )

  ! 4D case:
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_4D_i_min, layout_4D, i_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_4D_i_max, layout_4D, i_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_4D_j_min, layout_4D, j_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_4D_j_max, layout_4D, j_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_4D_k_min, layout_4D, k_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_4D_k_max, layout_4D, k_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_4D_l_min, layout_4D, l_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_4D_l_max, layout_4D, l_max )

  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_4D_i_min, layout_4D, i_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_4D_i_max, layout_4D, i_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_4D_j_min, layout_4D, j_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_4D_j_max, layout_4D, j_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_4D_k_min, layout_4D, k_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_4D_k_max, layout_4D, k_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_4D_l_min, layout_4D, l_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_4D_l_max, layout_4D, l_max )

  ! 6D case:
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_i_min, layout_6D, i_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_i_max, layout_6D, i_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_j_min, layout_6D, j_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_j_max, layout_6D, j_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_k_min, layout_6D, k_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_k_max, layout_6D, k_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_l_min, layout_6D, l_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_l_max, layout_6D, l_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_m_min, layout_6D, m_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_m_max, layout_6D, m_max )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_n_min, layout_6D, n_min )
  MAKE_GET_LAYOUT_SLOT_FUNCTION( get_layout_6D_n_max, layout_6D, n_max )

  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_i_min, layout_6D, i_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_i_max, layout_6D, i_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_j_min, layout_6D, j_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_j_max, layout_6D, j_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_k_min, layout_6D, k_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_k_max, layout_6D, k_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_l_min, layout_6D, l_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_l_max, layout_6D, l_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_m_min, layout_6D, m_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_m_max, layout_6D, m_max )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_n_min, layout_6D, n_min )
  MAKE_SET_LAYOUT_SLOT_FUNCTION( set_layout_6D_n_max, layout_6D, n_max )


  ! Why should lims just give its collective nilly-willy? This is not 
  ! pretty but I have the suspicion that direct access of the collective 
  ! will be needed.
#define MAKE_GET_LAYOUT_COLLECTIVE_FUNCTION( fname, layout_type ) \
  function fname( layout ); \
    intrinsic                       :: associated; \
    type(sll_collective_t), pointer :: fname; \
    type(layout_type), pointer      :: layout; \
    if( .not. associated(layout) ) then; \
       stop 'ERROR: uninitialized argument, get_layout_XD_collective()'; \
    end if; \
    fname => layout%collective; \
  end function fname

  MAKE_GET_LAYOUT_COLLECTIVE_FUNCTION( get_layout_2D_collective, layout_2D )
  MAKE_GET_LAYOUT_COLLECTIVE_FUNCTION( get_layout_3D_collective, layout_3D )
  MAKE_GET_LAYOUT_COLLECTIVE_FUNCTION( get_layout_4D_collective, layout_4D )
  MAKE_GET_LAYOUT_COLLECTIVE_FUNCTION( get_layout_6D_collective, layout_6D )

  ! get_layout_XD_size() returns the size of the collective associated
  ! with a given layout.
#define MAKE_GET_LAYOUT_SIZE_FUNCTION( fname, layout_type ) \
  function fname( layout ); \
    intrinsic                  :: associated; \
    sll_int32                  :: fname; \
    type(layout_type), pointer :: layout; \
    if( .not. associated(layout) ) then; \
       STOP 'ERROR: not associated argument passed to get_layout_size().'; \
    end if; \
    fname = sll_get_collective_size( layout%collective ); \
  end function fname

  MAKE_GET_LAYOUT_SIZE_FUNCTION(get_layout_2D_size, layout_2D)
  MAKE_GET_LAYOUT_SIZE_FUNCTION(get_layout_3D_size, layout_3D)
  MAKE_GET_LAYOUT_SIZE_FUNCTION(get_layout_4D_size, layout_4D)
  MAKE_GET_LAYOUT_SIZE_FUNCTION(get_layout_6D_size, layout_6D)

  ! Utility functions to help build layouts.
  function linear_index_2D(npx1, i, j)
    sll_int32, intent(in) :: npx1
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: linear_index_2D
    linear_index_2D = i + npx1*j
  end function linear_index_2D

  function linear_index_3D(npx1, npx2, i, j, k)
    sll_int32, intent(in) :: npx1
    sll_int32, intent(in) :: npx2
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32, intent(in) :: k
    sll_int32 :: linear_index_3D
    linear_index_3D = i + npx1*(j + npx2*k)
  end function linear_index_3D

  function linear_index_4D(npx1, npx2, npx3, i, j, k, l)
    sll_int32, intent(in) :: npx1
    sll_int32, intent(in) :: npx2
    sll_int32, intent(in) :: npx3
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32, intent(in) :: k
    sll_int32, intent(in) :: l
    sll_int32 :: linear_index_4D
    linear_index_4D = i + npx1*(j + npx2*(k + npx3*l))
  end function linear_index_4D

  function linear_index_6D(npx1, npx2, npx3, npx4, npx5, i, j, k, l, m, n)
    sll_int32, intent(in) :: npx1
    sll_int32, intent(in) :: npx2
    sll_int32, intent(in) :: npx3
    sll_int32, intent(in) :: npx4
    sll_int32, intent(in) :: npx5
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32, intent(in) :: k
    sll_int32, intent(in) :: l
    sll_int32, intent(in) :: m
    sll_int32, intent(in) :: n
    sll_int32 :: linear_index_6D
    linear_index_6D = i + npx1*(j + npx2*(k + npx3*(l + npx4*(m + npx5*n))))
  end function linear_index_6D

  subroutine initialize_layout_with_distributed_2D_array( &
    global_npx1, &  
    global_npx2, &
    num_proc_x1, &
    num_proc_x2, &
    layout )
    
    ! layout_2D should have been allocated with new(), which means that
    ! its memory is allocated in accordance with the size of collective.
    ! This should be error-checked below for consistency.
    sll_int32, intent(in) :: global_npx1
    sll_int32, intent(in) :: global_npx2
    sll_int32, intent(in) :: num_proc_x1
    sll_int32, intent(in) :: num_proc_x2
    type(layout_2D), pointer :: layout
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: total_num_processors
    sll_int32 :: node
    sll_int32 :: collective_size
    sll_int32 :: err
    sll_int32, dimension(:,:), allocatable :: intervals_x1
    sll_int32, dimension(:,:), allocatable :: intervals_x2

    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max

    if( &
       .not. is_power_of_two(int(num_proc_x1,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x2,i64)) ) then
       print *, 'ERROR: distribute_2D_array() needs that the integers that',&
            'describe the process mesh are powers of 2.'
       STOP
    end if

    if( &
       .not. (global_npx1 .gt. 0) .or. &
       .not. (global_npx2 .gt. 0) ) then
       print *, 'ERROR: distribute_2D_array() needs that the array dimensions',&
            'be greater than zero.'
       STOP
    end if

    ! FIXME: add further error checking, like a minimum number of points
    ! needed given a processor number along a dimension. Also, num_proc_xi
    ! should be different than zero.
    SLL_ALLOCATE( intervals_x1(0:1,0:num_proc_x1-1), err )
    SLL_ALLOCATE( intervals_x2(0:1,0:num_proc_x2-1), err )

    ! Allocate the layout to be returned.    
    total_num_processors = num_proc_x1*num_proc_x2
    collective_size = get_layout_2D_size(layout)
    if( total_num_processors .ne. collective_size ) then
       print *, 'ERROR, initialize_layout_with_distributed_2D_array(): ',&
            'requested size of the processor mesh is inconsistent with ', &
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

    ! Fill the layout array.
    do j=0, num_proc_x2-1
       do i=0, num_proc_x1-1
          node = linear_index_2D( num_proc_x1, i, j )
          i_min = intervals_x1(0,i)
          i_max = intervals_x1(1,i)
          j_min = intervals_x2(0,j)
          j_max = intervals_x2(1,j)
          call set_layout_i_min( layout, node, i_min )
          call set_layout_i_max( layout, node, i_max )
          call set_layout_j_min( layout, node, j_min )
          call set_layout_j_max( layout, node, j_max )
       end do
    end do
    SLL_DEALLOCATE_ARRAY( intervals_x1, err )
    SLL_DEALLOCATE_ARRAY( intervals_x2, err )
  end subroutine initialize_layout_with_distributed_2D_array

  subroutine initialize_layout_with_distributed_3D_array( &
    global_npx1, &  
    global_npx2, &
    global_npx3, &
    num_proc_x1, &
    num_proc_x2, &
    num_proc_x3, &
    layout )
    
    ! layout should have been allocated with new(), which means that
    ! its memory is allocated in accordance with the size of collective.
    ! This should be error-checked below for consistency.
    sll_int32, intent(in) :: global_npx1
    sll_int32, intent(in) :: global_npx2
    sll_int32, intent(in) :: global_npx3
    sll_int32, intent(in) :: num_proc_x1
    sll_int32, intent(in) :: num_proc_x2
    sll_int32, intent(in) :: num_proc_x3
    type(layout_3D), pointer :: layout
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
    collective_size = get_layout_3D_size(layout)
    if( total_num_processors .ne. collective_size ) then
       print *, 'ERROR, initialize_layout_with_distributed_3D_array(): ', &
            'requested size of the processor mesh is inconsistent with ', &
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
             call set_layout_i_min( layout, node, i_min )
             call set_layout_i_max( layout, node, i_max )
             call set_layout_j_min( layout, node, j_min )
             call set_layout_j_max( layout, node, j_max )
             call set_layout_k_min( layout, node, k_min )
             call set_layout_k_max( layout, node, k_max )
          end do
       end do
    end do
    SLL_DEALLOCATE_ARRAY( intervals_x1, err )
    SLL_DEALLOCATE_ARRAY( intervals_x2, err )
    SLL_DEALLOCATE_ARRAY( intervals_x3, err )
   end subroutine initialize_layout_with_distributed_3D_array


  subroutine initialize_layout_with_distributed_4D_array( &
    global_npx1, &  
    global_npx2, &
    global_npx3, &
    global_npx4, &
    num_proc_x1, &
    num_proc_x2, &
    num_proc_x3, &
    num_proc_x4, &
    layout )
    
    ! layout_4D should have been allocated with new(), which means that
    ! its memory is allocated in accordance with the size of collective.
    ! This should be error-checked below for consistency.
    sll_int32, intent(in) :: global_npx1
    sll_int32, intent(in) :: global_npx2
    sll_int32, intent(in) :: global_npx3
    sll_int32, intent(in) :: global_npx4
    sll_int32, intent(in) :: num_proc_x1
    sll_int32, intent(in) :: num_proc_x2
    sll_int32, intent(in) :: num_proc_x3
    sll_int32, intent(in) :: num_proc_x4
    type(layout_4D), pointer :: layout
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: k
    sll_int32 :: l
    sll_int32 :: total_num_processors
    sll_int32 :: node
    sll_int32 :: collective_size
    sll_int32 :: err
    sll_int32, dimension(:,:), allocatable :: intervals_x1
    sll_int32, dimension(:,:), allocatable :: intervals_x2
    sll_int32, dimension(:,:), allocatable :: intervals_x3
    sll_int32, dimension(:,:), allocatable :: intervals_x4
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: l_min
    sll_int32 :: l_max

    if( &
       .not. is_power_of_two(int(num_proc_x1,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x2,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x3,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x4,i64)) ) then
       print *, 'ERROR: distribute_4D_array() needs that the integers that',&
            'describe the process mesh are powers of 2.'
       STOP
    end if

    if( &
       .not. (global_npx1 .gt. 0) .or. &
       .not. (global_npx2 .gt. 0) .or. &
       .not. (global_npx3 .gt. 0) .or. &
       .not. (global_npx4 .gt. 0) ) then
       print *, 'ERROR: distribute_4D_array() needs that the array dimensions',&
            'be greater than zero.'
       STOP
    end if

    ! FIXME: add further error checking, like a minimum number of points
    ! needed given a processor number along a dimension. Also, num_proc_xi
    ! should be different than zero.
    SLL_ALLOCATE( intervals_x1(0:1,0:num_proc_x1-1), err )
    SLL_ALLOCATE( intervals_x2(0:1,0:num_proc_x2-1), err )
    SLL_ALLOCATE( intervals_x3(0:1,0:num_proc_x3-1), err )
    SLL_ALLOCATE( intervals_x4(0:1,0:num_proc_x4-1), err )

    ! Allocate the layout to be returned.    
    total_num_processors = num_proc_x1*num_proc_x2*num_proc_x3*num_proc_x4
    collective_size = get_layout_4D_size(layout)
    if( total_num_processors .ne. collective_size ) then
       print *, 'ERROR, initialize_layout_with_distributed_4D_array():', &
            'requested size of the processor mesh is inconsistent with ', &
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

    intervals_x4(0:1,0:num_proc_x4-1) = &
         split_array_indices( 1, global_npx4, num_proc_x4 )

    ! Fill the layout array.
    do l=0, num_proc_x4-1
       do k=0, num_proc_x3-1
          do j=0, num_proc_x2-1
             do i=0, num_proc_x1-1
                node = linear_index_4D( &
                     num_proc_x1, &
                     num_proc_x2, &
                     num_proc_x3, &
                     i, &
                     j, &
                     k, &
                     l )
                i_min = intervals_x1(0,i)
                i_max = intervals_x1(1,i)

                j_min = intervals_x2(0,j)
                j_max = intervals_x2(1,j)

                k_min = intervals_x3(0,k)
                k_max = intervals_x3(1,k)

                l_min = intervals_x4(0,l)
                l_max = intervals_x4(1,l)
                call set_layout_i_min( layout, node, i_min )
                call set_layout_i_max( layout, node, i_max )
                call set_layout_j_min( layout, node, j_min )
                call set_layout_j_max( layout, node, j_max )
                call set_layout_k_min( layout, node, k_min )
                call set_layout_k_max( layout, node, k_max )
                call set_layout_l_min( layout, node, l_min )
                call set_layout_l_max( layout, node, l_max )
             end do
          end do
       end do
    end do
    SLL_DEALLOCATE_ARRAY( intervals_x1, err )
    SLL_DEALLOCATE_ARRAY( intervals_x2, err )
    SLL_DEALLOCATE_ARRAY( intervals_x3, err )
    SLL_DEALLOCATE_ARRAY( intervals_x4, err )
  end subroutine initialize_layout_with_distributed_4D_array

  subroutine initialize_layout_with_distributed_6D_array( &
    global_npx1, &  
    global_npx2, &
    global_npx3, &
    global_npx4, &
    global_npx5, &
    global_npx6, &
    num_proc_x1, &
    num_proc_x2, &
    num_proc_x3, &
    num_proc_x4, &
    num_proc_x5, &
    num_proc_x6, &
    layout )
    
    ! layout_6D should have been allocated with new(), which means that
    ! its memory is allocated in accordance with the size of collective.
    ! This should be error-checked below for consistency.
    sll_int32, intent(in) :: global_npx1
    sll_int32, intent(in) :: global_npx2
    sll_int32, intent(in) :: global_npx3
    sll_int32, intent(in) :: global_npx4
    sll_int32, intent(in) :: global_npx5
    sll_int32, intent(in) :: global_npx6
    sll_int32, intent(in) :: num_proc_x1
    sll_int32, intent(in) :: num_proc_x2
    sll_int32, intent(in) :: num_proc_x3
    sll_int32, intent(in) :: num_proc_x4
    sll_int32, intent(in) :: num_proc_x5
    sll_int32, intent(in) :: num_proc_x6

    type(layout_6D), pointer :: layout
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: k
    sll_int32 :: l
    sll_int32 :: m
    sll_int32 :: n
    sll_int32 :: total_num_processors
    sll_int32 :: node
    sll_int32 :: collective_size
    sll_int32 :: err
    sll_int32, dimension(:,:), allocatable :: intervals_x1
    sll_int32, dimension(:,:), allocatable :: intervals_x2
    sll_int32, dimension(:,:), allocatable :: intervals_x3
    sll_int32, dimension(:,:), allocatable :: intervals_x4
    sll_int32, dimension(:,:), allocatable :: intervals_x5
    sll_int32, dimension(:,:), allocatable :: intervals_x6
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: l_min
    sll_int32 :: l_max
    sll_int32 :: m_min
    sll_int32 :: m_max
    sll_int32 :: n_min
    sll_int32 :: n_max

    if( &
       .not. is_power_of_two(int(num_proc_x1,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x2,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x3,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x4,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x5,i64)) .or. &
       .not. is_power_of_two(int(num_proc_x6,i64)) ) then
       print *, 'ERROR: distribute_6D_array() needs that the integers that',&
            'describe the process mesh are powers of 2.'
       STOP
    end if

    if( &
       .not. (global_npx1 .gt. 0) .or. &
       .not. (global_npx2 .gt. 0) .or. &
       .not. (global_npx3 .gt. 0) .or. &
       .not. (global_npx4 .gt. 0) .or. &
       .not. (global_npx5 .gt. 0) .or. &
       .not. (global_npx6 .gt. 0) ) then
       print *, 'ERROR: distribute_6D_array() needs that the array ', &
            'dimensions be greater than zero. Passed:', global_npx1, &
            global_npx2, global_npx3, global_npx4, global_npx5, global_npx6
       STOP
    end if

    ! FIXME: add further error checking, like a minimum number of points
    ! needed given a processor number along a dimension. Also, num_proc_xi
    ! should be different than zero.
    SLL_ALLOCATE( intervals_x1(0:1,0:num_proc_x1-1), err )
    SLL_ALLOCATE( intervals_x2(0:1,0:num_proc_x2-1), err )
    SLL_ALLOCATE( intervals_x3(0:1,0:num_proc_x3-1), err )
    SLL_ALLOCATE( intervals_x4(0:1,0:num_proc_x4-1), err )
    SLL_ALLOCATE( intervals_x5(0:1,0:num_proc_x5-1), err )
    SLL_ALLOCATE( intervals_x6(0:1,0:num_proc_x6-1), err )

    ! Allocate the layout to be returned.    
    total_num_processors = &
         num_proc_x1*num_proc_x2*num_proc_x3*num_proc_x4*num_proc_x5*num_proc_x6
    collective_size = get_layout_6D_size(layout)
    if( total_num_processors .ne. collective_size ) then
       print *, 'ERROR, initialize_layout_with_distributed_6D_array():', &
            'requested size of the processor mesh is inconsistent with ', &
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

    intervals_x4(0:1,0:num_proc_x4-1) = &
         split_array_indices( 1, global_npx4, num_proc_x4 )

    intervals_x5(0:1,0:num_proc_x5-1) = &
         split_array_indices( 1, global_npx5, num_proc_x5 )

    intervals_x6(0:1,0:num_proc_x6-1) = &
         split_array_indices( 1, global_npx6, num_proc_x6 )

    ! Fill the layout array.
    do n=0, num_proc_x6-1
       do m=0, num_proc_x5-1
          do l=0, num_proc_x4-1
             do k=0, num_proc_x3-1
                do j=0, num_proc_x2-1
                   do i=0, num_proc_x1-1
                      node = linear_index_6D( &
                           num_proc_x1, &
                           num_proc_x2, &
                           num_proc_x3, &
                           num_proc_x4, &
                           num_proc_x5, &
                           i, &
                           j, &
                           k, &
                           l, &
                           m, &
                           n )
                      i_min = intervals_x1(0,i)
                      i_max = intervals_x1(1,i)
                      
                      j_min = intervals_x2(0,j)
                      j_max = intervals_x2(1,j)
                      
                      k_min = intervals_x3(0,k)
                      k_max = intervals_x3(1,k)
                      
                      l_min = intervals_x4(0,l)
                      l_max = intervals_x4(1,l)
                      
                      m_min = intervals_x5(0,m)
                      m_max = intervals_x5(1,m)
                      
                      n_min = intervals_x6(0,n)
                      n_max = intervals_x6(1,n)
                      
                      call set_layout_i_min( layout, node, i_min )
                      call set_layout_i_max( layout, node, i_max )
                      call set_layout_j_min( layout, node, j_min )
                      call set_layout_j_max( layout, node, j_max )
                      call set_layout_k_min( layout, node, k_min )
                      call set_layout_k_max( layout, node, k_max )
                      call set_layout_l_min( layout, node, l_min )
                      call set_layout_l_max( layout, node, l_max )
                      call set_layout_m_min( layout, node, m_min )
                      call set_layout_m_max( layout, node, m_max )
                      call set_layout_n_min( layout, node, n_min )
                      call set_layout_n_max( layout, node, n_max )
                   end do
                end do
             end do
          end do
       end do
    end do
    SLL_DEALLOCATE_ARRAY( intervals_x1, err )
    SLL_DEALLOCATE_ARRAY( intervals_x2, err )
    SLL_DEALLOCATE_ARRAY( intervals_x3, err )
    SLL_DEALLOCATE_ARRAY( intervals_x4, err )
    SLL_DEALLOCATE_ARRAY( intervals_x5, err )
    SLL_DEALLOCATE_ARRAY( intervals_x6, err )
  end subroutine initialize_layout_with_distributed_6D_array

  function split_array_indices( min, max, num_intervals )
    sll_int32, intent(in)                       :: num_intervals
    sll_int32, dimension(0:1,0:num_intervals-1) :: split_array_indices
    sll_int32, intent(in)                       :: min
    sll_int32, intent(in)                       :: max
    sll_int32                                   :: num_elements
    num_elements = max - min + 1
    if( num_elements < num_intervals ) then
       print *, 'ERROR, split_array_indices(): the array given to split ', &
            'has less elements than the number of intervals requested. ', &
            'We have not implemented how to handle this case.'
       print *, 'number of elements: ', num_elements
       print *, 'number of intervals: ', num_intervals
       STOP 
    end if
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
#define MAKE_NEW_REMAP_PLAN_FUNCTION(fname, remap_type, layout_type, box_type, array_type, array_dim) \
  function fname( initial, final, array ); \
    intrinsic                       :: associated; \
    type(remap_type), pointer       :: fname; \
    type(layout_type), pointer      :: initial; \
    type(layout_type), pointer      :: final; \
    array_type, array_dim           :: array; \
    type(sll_collective_t), pointer :: coli; \
    type(sll_collective_t), pointer :: colf; \
    type(box_type)                  :: ibox, fbox, inters; \
    sll_int32                       :: i, f; \
    sll_int32                       :: my_rank; \
    sll_int32                       :: col_size; \
    sll_int32                       :: ierr; \
    sll_int32                       :: disp_counter; \
    sll_int32                       :: send_counter; \
    sll_int32                       :: recv_counter; \
    sll_int64                       :: acc; \
    if( (.not. associated(initial)) .or. (.not. associated(final)) ) then; \
       print *, 'ERROR: un-initialized arguments given to sll_new_remap_plan'; \
       stop; \
    end if; \
    coli => get_layout_collective(initial); \
    colf => get_layout_collective(final); \
    if( .not. collectives_are_same( coli, colf ) ) then; \
       print *, 'ERROR: init and final configurations given to new_remap_plan do not refer to the same collective.'; \
       stop; \
    end if; \
    acc = 0; \
    coli => get_layout_collective(initial); \
    my_rank  = sll_get_collective_rank( coli ); \
    col_size = sll_get_collective_size( coli ); \
    SLL_ALLOCATE( fname, ierr ); \
    SLL_ALLOCATE( fname%send_displs(0:col_size-1), ierr ); \
    fname%send_displs(:) = 0; \
    SLL_ALLOCATE( fname%send_counts(0:col_size-1), ierr ); \
    fname%send_counts(:) = 0; \
    SLL_ALLOCATE( fname%recv_displs(0:col_size-1), ierr ); \
    fname%recv_displs(:) = 0; \
    SLL_ALLOCATE( fname%recv_counts(0:col_size-1), ierr ); \
    fname%recv_counts(:) = 0; \
    SLL_ALLOCATE( fname%send_boxes(0:col_size-1), ierr ); \
    SLL_ALLOCATE( fname%recv_boxes(0:col_size-1), ierr ); \
    fname%collective => get_layout_collective(initial); \
    send_counter = 0; \
    disp_counter = 0; \
    ibox = get_layout_box(initial, my_rank); \
    fname%initial_layout => initial; \
    fname%final_layout   => final; \
    do f = 0, col_size-1; \
    fbox = get_layout_box(final, f); \
       if( intersect_boxes( ibox, fbox, inters ) ) then; \
          send_counter         = count_elements_in_box(inters); \
          fname%send_counts(f) = send_counter; \
          fname%send_displs(f) = disp_counter; \
          disp_counter         = disp_counter + send_counter; \
          fname%send_boxes(f)  = inters; \
          acc                  = acc + send_counter; \
       else; \
          fname%send_counts(f) = 0; \
          fname%send_displs(f) = disp_counter; \
          fname%send_boxes(f)  = inters; \
       end if; \
    end do; \
    SLL_ALLOCATE(fname%send_buffer(0:(acc-1)),ierr); \
    acc          = 0; \
    disp_counter = 0; \
    fbox = get_layout_box(final, my_rank); \
    do i = 0, col_size-1; \
       ibox = get_layout_box(initial,i); \
       if( intersect_boxes( ibox, fbox, inters ) ) then; \
          recv_counter         = count_elements_in_box(inters); \
          fname%recv_counts(i) = recv_counter; \
          fname%recv_displs(i) = disp_counter; \
          disp_counter         = disp_counter + recv_counter; \
          fname%recv_boxes(i)  = inters; \
          acc                  = acc + recv_counter; \
       else; \
          fname%recv_counts(i)   = 0; \
          fname%recv_displs(i) = disp_counter; \
          fname%recv_boxes(i)  = inters; \
       end if; \
    end do; \
    SLL_ALLOCATE(fname%recv_buffer(0:(acc-1)),ierr); \
    fname%is_uniform = .false.; \
    call optimize_remap_plan(fname); \
  end function fname

  ! 2D remaps
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_2D_int32,remap_plan_2D_int32,layout_2D,box_2D, sll_int32, dimension(:,:) )
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_2D_real64,remap_plan_2D_real64,layout_2D,box_2D, sll_real64, dimension(:,:) )
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_2D_comp64,remap_plan_2D_comp64,layout_2D,box_2D, sll_comp64, dimension(:,:) )
  ! 3D remaps
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_3D_int32,remap_plan_3D_int32,layout_3D,box_3D, sll_int32, dimension(:,:,:) )
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_3D_real64,remap_plan_3D_real64,layout_3D, box_3D, sll_real64, dimension(:,:,:) )
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_3D_comp64,remap_plan_3D_comp64,layout_3D,box_3D, sll_comp64, dimension(:,:,:) )
  ! 4D remaps
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_4D_int32,remap_plan_4D_int32,layout_4D,box_4D, sll_int32, dimension(:,:,:,:) )
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_4D_real64,remap_plan_4D_real64,layout_4D,box_4D, sll_real64, dimension(:,:,:,:) )
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_4D_comp64,remap_plan_4D_comp64,layout_4D,box_4D, sll_comp64, dimension(:,:,:,:) )
  ! 6D remaps
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_6D_int32,remap_plan_6D_int32,layout_6D,box_6D, sll_int32, dimension(:,:,:,:,:,:) )
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_6D_real64,remap_plan_6D_real64,layout_6D,box_6D, sll_real64, dimension(:,:,:,:,:,:) )
  MAKE_NEW_REMAP_PLAN_FUNCTION(new_remap_plan_6D_comp64,remap_plan_6D_comp64,layout_6D,box_6D, sll_comp64, dimension(:,:,:,:,:,:) )


 ! Try to fix the name of the subroutine in the print statement by stringifying
 ! the name.
#define MAKE_DELETE_REMAP_SUBROUTINE( fname, plan_type ) \
 subroutine fname( plan ); \
   type(plan_type), pointer :: plan; \
   sll_int32                    :: ierr; \
   if( .not. associated(plan) ) then; \
      print *, 'ERROR, delete_remap_plan(): passed plan was not associated.'; \
      stop; \
   end if; \
   SLL_DEALLOCATE(plan%send_displs, ierr); \
   SLL_DEALLOCATE(plan%send_counts, ierr); \
   SLL_DEALLOCATE(plan%recv_displs, ierr); \
   SLL_DEALLOCATE(plan%recv_counts, ierr); \
   SLL_DEALLOCATE(plan%send_boxes, ierr); \
   SLL_DEALLOCATE(plan%recv_boxes, ierr); \
   SLL_DEALLOCATE(plan%send_buffer, ierr); \
   SLL_DEALLOCATE(plan%recv_buffer, ierr); \
   SLL_DEALLOCATE(plan, ierr); \
 end subroutine fname

 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_2D_int32, remap_plan_2D_int32 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_2D_real64, remap_plan_2D_real64 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_2D_comp64, remap_plan_2D_comp64 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_3D_int32, remap_plan_3D_int32 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_3D_real64, remap_plan_3D_real64 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_3D_comp64, remap_plan_3D_comp64 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_4D_int32, remap_plan_4D_int32 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_4D_real64, remap_plan_4D_real64 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_4D_comp64, remap_plan_4D_comp64 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_6D_int32, remap_plan_6D_int32 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_6D_real64, remap_plan_6D_real64 )
 MAKE_DELETE_REMAP_SUBROUTINE( delete_remap_6D_comp64, remap_plan_6D_comp64 )

#if 0
  ! We leave this function here for reference, as it was the original and
  ! has comments.
  function new_remap_plan_3D( initial, final, int32_data_size )
    intrinsic                       :: associated
    type(remap_plan_3D), pointer  :: new_remap_plan_3D 
    type(layout_3D), pointer      :: initial
    type(layout_3D), pointer      :: final
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
    SLL_ALLOCATE( new_remap_plan_3D%send_displs(0:col_size-1), ierr )
    new_remap_plan_3D%send_displs(:) = 0
    SLL_ALLOCATE( new_remap_plan_3D%send_counts(0:col_size-1), ierr )
    new_remap_plan_3D%send_counts(:) = 0
    SLL_ALLOCATE( new_remap_plan_3D%recv_displs(0:col_size-1), ierr )
    new_remap_plan_3D%recv_displs(:) = 0
    SLL_ALLOCATE( new_remap_plan_3D%recv_counts(0:col_size-1), ierr )
    new_remap_plan_3D%recv_counts(:) = 0
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
    call optimize_remap_plan_3D(new_remap_plan_3D)
  end function new_remap_plan_3D
#endif

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
#define MAKE_REMAP_OPTIMIZER( fname, remap_type, box_type ) \
  subroutine fname( plan ); \
    type(remap_type), pointer            :: plan; \
    sll_int32, dimension(:), pointer     :: send_counts; \
    sll_int32, dimension(:), pointer     :: send_displs; \
    sll_int32, dimension(:), pointer     :: recv_counts; \
    sll_int32, dimension(:), pointer     :: recv_displs; \
    type(sll_collective_t), pointer      :: col; \
    sll_int32                            :: col_sz; \
    sll_int32, dimension(:), allocatable :: lowest_color; \
    sll_int32, dimension(:), allocatable :: colors; \
    sll_int32, dimension(:), allocatable :: colors_copy; \
    sll_int32                            :: ierr; \
    sll_int32                            :: my_rank; \
    sll_int32                            :: i; \
    type(sll_collective_t), pointer      :: new_collective; \
    sll_int32                            :: new_col_sz; \
    sll_int32, dimension(:), pointer     :: new_send_counts; \
    sll_int32, dimension(:), pointer     :: new_send_displs; \
    sll_int32, dimension(:), pointer     :: new_recv_counts; \
    sll_int32, dimension(:), pointer     :: new_recv_displs; \
    type(box_type), dimension(:), pointer :: new_send_boxes; \
    type(box_type), dimension(:), pointer :: new_recv_boxes; \
    sll_int32                            :: new_i; \
    sll_int32                            :: my_color; \
    sll_int32                            :: exchange_size; \
    logical, dimension(1:1)              :: is_uniform_local; \
    logical, dimension(1:1)              :: is_uniform_collective; \
    sll_int32                            :: new_sdisp; \
    sll_int32                            :: new_rdisp; \
    col         => plan%collective; \
    col_sz      = sll_get_collective_size( col ); \
    my_rank     = sll_get_collective_rank( col ); \
    send_counts => plan%send_counts; \
    send_displs => plan%send_displs; \
    recv_counts => plan%recv_counts; \
    recv_displs => plan%recv_displs; \
    SLL_ALLOCATE( lowest_color(1), ierr ); \
    lowest_color(1) = 0; \
    SLL_ALLOCATE( colors(0:col_sz-1), ierr ); \
    colors(:) = 0; \
    SLL_ALLOCATE( colors_copy(0:col_sz-1), ierr ); \
    colors_copy(:) = 0; \
    lowest_color(1) = my_rank; \
    call sll_collective_allgather(col,lowest_color,1,colors(0:col_sz-1),1); \
    do; \
       colors_copy(0:col_sz-1) = colors(0:col_sz-1); \
       do i=0,col_sz-1; \
          if( (send_counts(i) .ne. 0) .or. (recv_counts(i) .ne. 0) ) then; \
             if( colors(i) .lt. lowest_color(1) ) then; \
                lowest_color(1) = colors(i); \
             end if; \
          end if; \
       end do; \
       call sll_collective_allgather(col,lowest_color,1,colors(0:col_sz-1),1); \
       if(arrays_are_equal(colors, colors_copy, col_sz)) then; \
          exit; \
       end if; \
    end do; \
    new_collective => sll_new_collective( col, colors(my_rank), my_rank ); \
    new_col_sz     = sll_get_collective_size( new_collective ); \
    SLL_ALLOCATE( new_send_counts(0:new_col_sz-1), ierr ); \
    SLL_ALLOCATE( new_send_displs(0:new_col_sz-1), ierr ); \
    SLL_ALLOCATE( new_recv_counts(0:new_col_sz-1), ierr ); \
    SLL_ALLOCATE( new_recv_displs(0:new_col_sz-1), ierr ); \
    SLL_ALLOCATE( new_send_boxes( 0:new_col_sz-1), ierr ); \
    SLL_ALLOCATE( new_recv_boxes( 0:new_col_sz-1), ierr ); \
    new_i = 0; \
    my_color = colors(my_rank); \
    new_sdisp = 0; \
    new_rdisp = 0; \
    do i=0,col_sz-1; \
       if( colors(i) .eq. my_color ) then; \
          new_send_counts(new_i) = send_counts(i); \
          new_send_displs(new_i) = new_sdisp; \
          new_send_boxes(new_i)  = plan%send_boxes(i); \
          new_sdisp              = new_sdisp + send_counts(i); \
          new_recv_counts(new_i) = recv_counts(i); \
          new_recv_displs(new_i) = new_rdisp; \
          new_recv_boxes(new_i)  = plan%recv_boxes(i); \
          new_rdisp              = new_rdisp + recv_counts(i); \
          new_i                  = new_i + 1; \
       end if; \
    end do; \
    plan%collective => new_collective; \
    SLL_DEALLOCATE( plan%send_counts, ierr ); \
    plan%send_counts => new_send_counts; \
    SLL_DEALLOCATE( plan%send_displs, ierr ); \
    plan%send_displs => new_send_displs; \
    SLL_DEALLOCATE( plan%recv_counts, ierr ); \
    plan%recv_counts => new_recv_counts; \
    SLL_DEALLOCATE( plan%recv_displs, ierr ); \
    plan%recv_displs => new_recv_displs; \
    SLL_DEALLOCATE( plan%send_boxes, ierr ); \
    plan%send_boxes => new_send_boxes; \
    SLL_DEALLOCATE( plan%recv_boxes, ierr ); \
    plan%recv_boxes => new_recv_boxes; \
    SLL_DEALLOCATE_ARRAY( lowest_color, ierr ); \
    SLL_DEALLOCATE_ARRAY( colors, ierr ); \
    SLL_DEALLOCATE_ARRAY( colors_copy, ierr ); \
    exchange_size = plan%send_counts(0); \
    do i=0,new_col_sz-1; \
       if(plan%send_counts(i) .eq. exchange_size) then; \
          is_uniform_local(1) = .true.; \
       else; \
          is_uniform_local(1) = .false.; \
          exit; \
       end if; \
    end do; \
    call sll_collective_allreduce(plan%collective,is_uniform_local(:),1,MPI_LAND, is_uniform_collective(:) ); \
    plan%is_uniform = is_uniform_collective(1); \
  end subroutine fname

  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_2D_int32, remap_plan_2D_int32,box_2D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_2D_real64, remap_plan_2D_real64,box_2D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_2D_comp64, remap_plan_2D_comp64,box_2D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_3D_int32, remap_plan_3D_int32,box_3D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_3D_real64, remap_plan_3D_real64,box_3D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_3D_comp64, remap_plan_3D_comp64,box_3D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_4D_int32, remap_plan_4D_int32,box_4D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_4D_real64, remap_plan_4D_real64,box_4D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_4D_comp64, remap_plan_4D_comp64,box_4D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_6D_int32, remap_plan_6D_int32,box_6D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_6D_real64, remap_plan_6D_real64,box_6D)
  MAKE_REMAP_OPTIMIZER(optimize_remap_plan_6D_comp64, remap_plan_6D_comp64,box_6D)

#if 0
  ! We leave this here for reference, as it was the original and has comments.
  subroutine optimize_remap_plan_3D( plan ) !, sub_collective )
    type(remap_plan_3D), pointer         :: plan
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
    sll_int32                            :: new_sdisp
    sll_int32                            :: new_rdisp
    col         => plan%collective
    col_sz      = sll_get_collective_size( col )
    my_rank     = sll_get_collective_rank( col )
    send_counts => plan%send_counts
    send_displs => plan%send_displs
    recv_counts => plan%recv_counts
    recv_displs => plan%recv_displs
    SLL_ALLOCATE( lowest_color(1), ierr ) ! awful, but I need an array.
    lowest_color(1) = 0
    SLL_ALLOCATE( colors(0:col_sz-1), ierr )
    colors(:) = 0
    SLL_ALLOCATE( colors_copy(0:col_sz-1), ierr )
    colors_copy = 0
    ! FIRST LEVEL OF OPTIMIZATION: 
    ! Identify the sub-collectives in which the communication should 
    ! be divided. The purpose is to subdivide the original communicator
    ! into multiple communicators, each one minimally sized, but meeting
    ! the condition that each process belongs to a single communicator.
    ! In some situations, this could lead to cases in which two processes
    ! can belong to a communicator even though they do not exchange data
    ! amongst themselves directly, but need to exchange data with a
    ! third process. Even this situation might still not be necessarily
    ! slower, than having the split communicators, and here we have much
    ! simpler code.

    ! we want a starting point. This should not change for the lowest ranks 
    ! in the sub-collectives.
    lowest_color(1) = my_rank  
    call sll_collective_allgather( &
       col, &
       lowest_color(1:1), &
       1, &
       colors(0:col_sz-1), &
       1 )
    do
       ! Load the copy
       colors_copy(0:col_sz-1) = colors(0:col_sz-1)
       ! Find the lowest rank with which this process communicates
       do i=0,col_sz-1
          if( (send_counts(i) .ne. 0) .or. (recv_counts(i) .ne. 0) ) then
             if( colors(i) .lt. lowest_color(1) ) then 
                lowest_color(1) = colors(i) 
             end if  ! else, nothing, as lowest_color is already my_rank
          end if
       end do
       ! Gather the information from all processes
       call sll_collective_allgather( &
            col, &
            lowest_color(1:1), &
            1, &
            colors(0:col_sz-1), &
            1 )
#if 0
       print *, my_rank, 'colors intermediate: ', colors(:) ! delete this
       ! Check which is the color of the lowest rank with which this
       ! process communicates, and reassign the color accordingly.
       lowest_color(1) = colors(lowest_color(1))

       call sll_collective_allgather( &
            col, &
            lowest_color(1:1), &
            1, &
            colors(0:col_sz-1), &
            1 )
#endif
       if(arrays_are_equal(colors, colors_copy, col_sz)) exit
    end do
#if 0
    if(my_rank .eq. 0) then
       print *, 'final colors array: ', colors(:) ! delete this
    end if
#endif
    ! The results can now be used as the color for a collective-splitting 
    ! operation.
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
    new_sdisp = 0
    new_rdisp = 0
    do i=0,col_sz-1
       if( colors(i) .eq. my_color ) then
          new_send_counts(new_i) = send_counts(i)
          new_send_displs(new_i) = new_sdisp
          new_send_boxes(new_i)  = plan%send_boxes(i)
          new_sdisp              = new_sdisp + send_counts(i)
          new_recv_counts(new_i) = recv_counts(i)
          new_recv_displs(new_i) = new_rdisp
          new_recv_boxes(new_i)  = plan%recv_boxes(i)
          new_rdisp              = new_rdisp + recv_counts(i)
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
    call flush(6)
#endif
  end subroutine optimize_remap_plan_3D
#endif

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

#define MAKE_GET_REMAP_INITIAL_LAYOUT( fname, layout_type, remap_type) \
  function fname( plan ); \
    type(layout_type), pointer :: fname; \
    type(remap_type), pointer  :: plan; \
    if( .not. associated( plan ) ) then; \
       write (*,'(a)') 'not associated pointer argument'; \
       stop 'get_remap_initial_layout'; \
    end if; \
    fname => plan%initial_layout; \
  end function fname

  MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_2D_initial_layout_int32, layout_2D, remap_plan_2D_int32)
  MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_2D_initial_layout_real64, layout_2D, remap_plan_2D_real64)
  MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_2D_initial_layout_comp64, layout_2D, remap_plan_2D_comp64)
  MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_3D_initial_layout_int32, layout_3D, remap_plan_3D_int32)
 MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_3D_initial_layout_real64, layout_3D, remap_plan_3D_real64)
 MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_3D_initial_layout_comp64, layout_3D, remap_plan_3D_comp64)
  MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_4D_initial_layout_int32, layout_4D, remap_plan_4D_int32)
  MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_4D_initial_layout_real64, layout_4D, remap_plan_4D_real64)
  MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_4D_initial_layout_comp64, layout_4D, remap_plan_4D_comp64)
  MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_6D_initial_layout_int32, layout_6D, remap_plan_6D_int32)
  MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_6D_initial_layout_real64, layout_6D, remap_plan_6D_real64)
  MAKE_GET_REMAP_INITIAL_LAYOUT(get_remap_6D_initial_layout_comp64, layout_6D, remap_plan_6D_comp64)

#define MAKE_GET_REMAP_FINAL_LAYOUT( fname, layout_type, remap_type ) \
  function fname( plan ); \
    type(layout_type), pointer :: fname; \
    type(remap_type), pointer  :: plan;  \
    if( .not. associated( plan ) ) then; \
       write (*,'(a)') 'not associated pointer argument'; \
       stop 'get_remap_final_layout'; \
    end if; \
    fname => plan%final_layout; \
  end function fname

  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_2D_final_layout_int32,layout_2D,remap_plan_2d_int32)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_2D_final_layout_real64,layout_2D,remap_plan_2d_real64)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_2D_final_layout_comp64,layout_2D,remap_plan_2d_comp64)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_3D_final_layout_int32,layout_3D,remap_plan_3d_int32)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_3D_final_layout_real64,layout_3D,remap_plan_3d_real64)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_3D_final_layout_comp64,layout_3D,remap_plan_3d_comp64)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_4D_final_layout_int32,layout_4D,remap_plan_4d_int32)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_4D_final_layout_real64,layout_4D,remap_plan_4d_real64)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_4D_final_layout_comp64,layout_4D,remap_plan_4d_comp64)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_6D_final_layout_int32,layout_6D,remap_plan_6d_int32)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_6D_final_layout_real64,layout_6D,remap_plan_6d_real64)
  MAKE_GET_REMAP_FINAL_LAYOUT(get_remap_6D_final_layout_comp64,layout_6D,remap_plan_6d_comp64)

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
    SLL_ASSERT( n <= size(ai) )
    SLL_ASSERT( n <= size(bi) )
    do i=1,n
       bi(i) = ai(i)*sz
    end do
  end subroutine convert_into_integer_sizes

  ! **********************************************************************
  !
  !    Continue here with 3D, 4D and 5D functions...
  !
  ! **********************************************************************
  subroutine apply_remap_3D_int( plan, data_in, data_out )
    type(remap_plan_3D_int32), pointer             :: plan
    sll_int32, dimension(:,:,:), intent(in)  :: data_in
    sll_int32, dimension(:,:,:), intent(out) :: data_out
    sll_int32, dimension(:), pointer         :: sb       ! send buffer
    sll_int32, dimension(:), pointer         :: rb       ! receive buffer
    sll_int32, dimension(:), pointer         :: sdisp    ! send displacements
    sll_int32, dimension(:), pointer         :: rdisp    ! receive displacements
    sll_int32, dimension(:), pointer         :: scnts    ! send counts
    sll_int32, dimension(:), pointer         :: rcnts    ! receive counts
    type(sll_collective_t), pointer          :: col      ! collective
    type(layout_3D), pointer                 :: init_layout  => NULL()
    type(layout_3D), pointer                 :: final_layout => NULL()
    sll_int32                                :: id, jd, kd
    sll_int32                                :: i
    sll_int32                                :: col_sz
    sll_int32                                :: loi, loj, lok
    sll_int32                                :: hii, hij, hik
    type(box_3D)                             :: sbox
    sll_int32                                :: my_rank
    sll_int32                                :: loc
    sll_int32, dimension(1:3)                :: local_lo, local_hi
    sll_int32, dimension(1:3)                :: tmpa

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp        => plan%send_displs
    rdisp        => plan%recv_displs
    scnts        => plan%send_counts
    rcnts        => plan%recv_counts
    col          => plan%collective
    col_sz       =  sll_get_collective_size(col)
    init_layout  => get_remap_3D_initial_layout_int32(plan)
    final_layout => get_remap_3D_final_layout_int32(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer
    ! load the send buffer
    loc = 0             ! first loading is at position zero
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdisp(i) ) then
             write (*,'(a,i4,a,i16)') 'apply_remap_3D_int() ERROR: discrepancy between displs(i) and the loading index for i = ', i, ' displs(i) = ', sdisp(i)
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
          tmpa(:)  = (/loi,loj,lok/)
          local_lo = global_to_local_3D( init_layout, tmpa)
          tmpa(:)  = (/hii,hij,hik/)
          local_hi = global_to_local_3D( init_layout, tmpa )

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
                   sb(loc) = data_in(id,jd,kd)
                   loc     = loc + 1
                end do
             end do
          end do
       end if
    end do
    
!    write (*,'(a,i4)') 'the send buffer in rank:', my_rank
!    print *, sb(0:(size(sb)-1))
!    call flush(6)
 
   if( plan%is_uniform .eqv. .false. ) then 
       call sll_collective_alltoallV( sb(:),       &
                                      scnts(0:col_sz-1), &
                                      sdisp(0:col_sz-1), &
                                      rb(:),       &
                                      rcnts(0:col_sz-1), &
                                      rdisp(0:col_sz-1), col )
    else
       call sll_collective_alltoall ( sb(:), &
                                      scnts(0), &
                                      rcnts(0), &
                                      rb(:), col )
    end if
!    write (*,'(a, i4)') 'the receive buffer in rank: ', my_rank
!    print *, rb(0:size(rb)-1)
!    call flush(6)
    ! Unpack the plan into the outgoing buffer.
    loc = 0  ! We load first from position 0 in the receive buffer.
    do i = 0, col_sz-1
       if( rcnts(i) .ne. 0 ) then ! we expect something from rank 'i'
          if( loc .ne. rdisp(i) ) then
             write (*,'(a,i4)') &
                  'ERROR: discrepancy between rdisp(i) and index for i = ', i
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
          tmpa(:)  = (/loi,loj,lok/)
          local_lo = global_to_local_3D( final_layout, tmpa )
          tmpa(:)  = (/hii,hij,hik/)
          local_hi = global_to_local_3D( final_layout, tmpa )
          do kd = local_lo(3), local_hi(3)
             do jd = local_lo(2), local_hi(2)
                do id = local_lo(1), local_hi(1)
                   data_out(id,jd,kd) = rb(loc)
                   loc                = loc + 1
                end do
             end do
          end do
       end if
    end do
  end subroutine apply_remap_3D_int

  subroutine apply_remap_2D_double( plan, data_in, data_out )
    type(remap_plan_2D_real64), pointer       :: plan
    sll_real64, dimension(:,:), intent(in)    :: data_in
    sll_real64, dimension(:,:), intent(out)   :: data_out
    sll_real64, dimension(:), pointer         :: sb     ! send buffer
    sll_real64, dimension(:), pointer         :: rb     ! receive buffer
    sll_int32, dimension(:), pointer          :: sdisp  ! send displacements
    sll_int32, dimension(:), pointer          :: rdisp  ! receive displacements 
    sll_int32, dimension(:), pointer          :: scnts  ! send counts
    sll_int32, dimension(:), pointer          :: rcnts  ! receive counts
    type(sll_collective_t), pointer           :: col    ! collective
    type(layout_2D), pointer                  :: init_layout  => NULL()
    type(layout_2D), pointer                  :: final_layout => NULL()
    sll_int32                                 :: id, jd
    sll_int32                                 :: i
    sll_int32                                 :: col_sz
    sll_int32                                 :: loi, loj
    sll_int32                                 :: hii, hij
    type(box_2D)                              :: sbox
    sll_int32                                 :: my_rank
    sll_int32                                 :: loc
    sll_int32, dimension(1:2)                 :: local_lo, local_hi
    sll_int32, dimension(1:2)                 :: tmpa

!!$    ! to load the MPI function and send integers, we have a separate set of
!!$    ! arrays to store this information for now.
!!$    sll_int32, dimension(:), allocatable     :: sdispi  ! send displacements
!!$    sll_int32, dimension(:), allocatable     :: rdispi  ! receive displacements
!!$    sll_int32, dimension(:), allocatable     :: scntsi  ! send counts
!!$    sll_int32, dimension(:), allocatable     :: rcntsi  ! receive counts

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp        => plan%send_displs
    rdisp        => plan%recv_displs
    scnts        => plan%send_counts
    rcnts        => plan%recv_counts
    col          => plan%collective
    col_sz       =  sll_get_collective_size(col)
    init_layout  => get_remap_2D_initial_layout_real64(plan)
    final_layout => get_remap_2D_final_layout_real64(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer

!!$    SLL_ALLOCATE(sdispi(0:col_sz-1), ierr)
!!$    SLL_ALLOCATE(rdispi(0:col_sz-1), ierr)
!!$    SLL_ALLOCATE(scntsi(0:col_sz-1), ierr)
!!$    SLL_ALLOCATE(rcntsi(0:col_sz-1), ierr)

    ! Translate the amounts into integers
!!$#if 1
!!$    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1)), sdisp, &
!!$         col_sz, sdispi)
!!$    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1)), rdisp, &
!!$         col_sz, rdispi)
!!$    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1)), scnts, &
!!$         col_sz, scntsi)
!!$    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1)), rcnts, &
!!$         col_sz, rcntsi)
!!$#endif
    
!!$#if 0
!!$    write (*,'(a,i4)') 'parameters from rank ', my_rank
!!$    print *, 'scntsi', scntsi(:)
!!$    print *, 'sdispi', sdispi(:)
!!$    print *, 'rcntsi', rcntsi(:)
!!$    print *, 'rdispi', rdispi(:)
!!$    call flush(6)
!!$#endif

    ! load the send buffer
    loc = 0             ! first loading is at position zero
!!$    ! This step is obviously not needed for integers themselves. We put this
!!$    ! here for generality.
!!$    int32_data_size = INT32_SIZEOF( data_in(1,1) )
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
!!$          if( loc .ne. sdispi(i) ) then
          if( loc .ne. sdisp(i) ) then
             print *, 'ERROR DETECTED in process: ', my_rank
             print *, 'apply_remap_2D_double() ERROR: ', &
                  'discrepancy between displs(i) and the loading index for ',&
!!$                  'i = ', i, ' displs(i) = ', sdispi(i)
                  'i = ', i, ' displs(i) = ', sdisp(i)
             write(*,'(a,i8)') 'col_sz = ', col_sz
             call flush(6)
             stop 'apply_remap(): loading error'
          end if
          ! get the information on the box to send, get the limits,
          ! convert to the local indices and find out where in the 
          ! buffer to start writing.
          sbox = plan%send_boxes(i)
          loi = get_box_2D_i_min(sbox)
          loj = get_box_2D_j_min(sbox)
          hii = get_box_2D_i_max(sbox)
          hij = get_box_2D_j_max(sbox)
          tmpa(:)  = (/loi,loj/)
          local_lo = global_to_local_2D( init_layout, tmpa )
          tmpa(:)  = (/hii,hij/)
          local_hi = global_to_local_2D( init_layout, tmpa )

!!$          local_lo = global_to_local_2D( init_layout, (/loi,loj/) )
!!$          local_hi = global_to_local_2D( init_layout, (/hii,hij/)  )

          ! The plan to load the send buffer is to traverse the integer
          ! array with a single index (loc). When we load the buffer, each
          ! data element may occupy multiple integer 'slots', hence the
          ! loading index needs to be manually increased. As an advantage,
          ! we can do some error checking every time we send data to a 
          ! different process, as we know what is the expected value of 
          ! the index at that point.
          do jd = local_lo(2), local_hi(2)
             do id = local_lo(1), local_hi(1)
!!$                sb(loc:) = transfer(data_in(id,jd),(/1_i32/))
                sb(loc) = data_in(id,jd)
!!$                loc      = loc + int32_data_size
                loc      = loc + 1
             end do
          end do
       end if
    end do

    ! Comment the following when not debugging    
    !   write (*,'(a,i4)') 'the send buffer in rank:', my_rank
    !  print *, sb(0:(size(sb)-1))
    ! call flush(6)
    !    print *, 'from inside remap: rank ', my_rank, 'calling communications'
    !    call flush(6)
   if( plan%is_uniform .eqv. .false. ) then 
      ! the following call can be changed from a generic to a type-specific
      ! call when right away, but especially if the apply_remap function gets
      ! specialized (i.e. gets rid of transfer() calls).
       call sll_collective_alltoallV( sb(:),       &
!!$                                      scntsi(0:col_sz-1), &
!!$                                      sdispi(0:col_sz-1), &
                                      scnts(0:col_sz-1), &
                                      sdisp(0:col_sz-1), &
                                      rb(:),       &
!!$                                      rcntsi(0:col_sz-1), &
!!$                                      rdispi(0:col_sz-1), &
                                      rcnts(0:col_sz-1), &
                                      rdisp(0:col_sz-1), &
                                      col )
    else
       call sll_collective_alltoall ( sb(:), &
!!$                                      scntsi(0), &
!!$                                      rcntsi(0), &
                                      scnts(0), &
                                      rcnts(0), &
                                      rb(:), col )
    end if
!    write (*,'(a, i4)') 'the receive buffer in rank: ', my_rank
!    print *, rb(0:size(rb)-1)
!    call flush(6)
    ! Unpack the plan into the outgoing buffer.
    loc = 0  ! We load first from position 0 in the receive buffer.
    do i = 0, col_sz-1
       if( rcnts(i) .ne. 0 ) then ! we expect something from rank 'i'
!!$          if( loc .ne. rdispi(i) ) then
          if( loc .ne. rdisp(i) ) then
             write (*,'(a,i4)') &
                  'ERROR: discrepancy between rdisp(i) and index for i = ', i
             stop 'unpacking error'
          end if
          ! get the information on the box to receive, get the limits, and 
          ! convert to the local indices.
          sbox = plan%recv_boxes(i)
          loi = get_box_2D_i_min(sbox)
          loj = get_box_2D_j_min(sbox)
          hii = get_box_2D_i_max(sbox)
          hij = get_box_2D_j_max(sbox)
          tmpa(:)  = (/loi,loj/)
          local_lo = global_to_local_2D( final_layout, tmpa )
          tmpa(:)  = (/hii,hij/)
          local_hi = global_to_local_2D( final_layout, tmpa )
          do jd = local_lo(2), local_hi(2)
             do id = local_lo(1), local_hi(1)
!!$                data_out(id,jd) = transfer(rb(loc:),data_out(1,1))
                data_out(id,jd) = rb(loc)
!!$                loc                = loc + int32_data_size
                loc                = loc + 1
             end do
          end do
       end if
    end do
!!$    ! And why weren't these arrays part of the plan anyway??
!!$    SLL_DEALLOCATE_ARRAY(sdispi, ierr)
!!$    SLL_DEALLOCATE_ARRAY(rdispi, ierr)
!!$    SLL_DEALLOCATE_ARRAY(scntsi, ierr)
!!$    SLL_DEALLOCATE_ARRAY(rcntsi, ierr)
  end subroutine apply_remap_2D_double

  subroutine apply_remap_2D_complex( plan, data_in, data_out )
    type(remap_plan_2D_comp64), pointer       :: plan
    sll_comp64, dimension(:,:), intent(in)    :: data_in
    sll_comp64, dimension(:,:), intent(out)   :: data_out
    sll_comp64, dimension(:), pointer         :: sb     ! send buffer
    sll_comp64, dimension(:), pointer         :: rb     ! receive buffer
    sll_int32, dimension(:), pointer          :: sdisp  ! send displacements
    sll_int32, dimension(:), pointer          :: rdisp  ! receive displacements 
    sll_int32, dimension(:), pointer          :: scnts  ! send counts
    sll_int32, dimension(:), pointer          :: rcnts  ! receive counts
    type(sll_collective_t), pointer           :: col    ! collective
    type(layout_2D), pointer                  :: init_layout  => NULL()
    type(layout_2D), pointer                  :: final_layout => NULL()
    sll_int32                                 :: id, jd
    sll_int32                                 :: i
    sll_int32                                 :: col_sz
    sll_int32                                 :: loi, loj
    sll_int32                                 :: hii, hij
    type(box_2D)                              :: sbox
    sll_int32                                 :: my_rank
    sll_int32                                 :: loc
    sll_int32, dimension(1:2)                 :: local_lo, local_hi
    sll_int32, dimension(1:2)                 :: tmpa

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp        => plan%send_displs
    rdisp        => plan%recv_displs
    scnts        => plan%send_counts
    rcnts        => plan%recv_counts
    col          => plan%collective
    col_sz       =  sll_get_collective_size(col)
    init_layout  => get_remap_2D_initial_layout_comp64(plan)
    final_layout => get_remap_2D_final_layout_comp64(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer

#if 0
print *, 'remap 2d complex:'
    print *, 'scnts = ', scnts(:)
    print *, 'rcnts = ', rcnts(:)
    print *, 'sdisp = ', sdisp(:)
    print *, 'rdisp = ', rdisp(:)
#endif
    ! load the send buffer
    loc = 0             ! first loading is at position zero
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdisp(i) ) then
             print *, 'ERROR DETECTED in process: ', my_rank
             print *, 'apply_remap_2D_complex() ERROR: ', &
                  'discrepancy between displs(i) and the loading index for ',&
                  'i = ', i, ' displs(i) = ', sdisp(i)
             write(*,'(a,i8)') 'col_sz = ', col_sz
             call flush(6)
             stop 'apply_remap(): loading error'
          end if
          ! get the information on the box to send, get the limits,
          ! convert to the local indices and find out where in the 
          ! buffer to start writing.
          sbox = plan%send_boxes(i)
          loi = get_box_2D_i_min(sbox)
          loj = get_box_2D_j_min(sbox)
          hii = get_box_2D_i_max(sbox)
          hij = get_box_2D_j_max(sbox)
          tmpa(:)  = (/loi,loj/)
          local_lo = global_to_local_2D( init_layout, tmpa )
          tmpa(:)  = (/hii,hij/)
          local_hi = global_to_local_2D( init_layout, tmpa )

          ! The plan to load the send buffer is to traverse the integer
          ! array with a single index (loc). When we load the buffer, each
          ! data element may occupy multiple integer 'slots', hence the
          ! loading index needs to be manually increased. As an advantage,
          ! we can do some error checking every time we send data to a 
          ! different process, as we know what is the expected value of 
          ! the index at that point.
          do jd = local_lo(2), local_hi(2)
             do id = local_lo(1), local_hi(1)
                sb(loc) = data_in(id,jd)
                loc     = loc + 1
             end do
          end do
       end if
    end do
    ! Comment the following when not debugging    
!!$    write (*,'(a,i4)') 'the send buffer in rank:', my_rank
!!$    print *, sb(0:(size(sb)-1))
!!$    call flush(6)
!!$    print *, 'from inside remap: rank ', my_rank, 'calling communications'
!!$    call flush(6)

   if( plan%is_uniform .eqv. .false. ) then 
      ! the following call can be changed from a generic to a type-specific
      ! call when right away, but especially if the apply_remap function gets
      ! specialized (i.e. gets rid of transfer() calls).
       call sll_collective_alltoallV( sb(:),       &
                                      scnts(0:col_sz-1), &
                                      sdisp(0:col_sz-1), &
                                      rb(:),       &
                                      rcnts(0:col_sz-1), &
                                      rdisp(0:col_sz-1), col )
    else
       call sll_collective_alltoall ( sb(:), &
                                      scnts(0), &
                                      rcnts(0), &
                                      rb(:), col )
    end if

    ! Comment when not debugging:
!!$    write (*,'(a, i4)') 'the receive buffer in rank: ', my_rank
!!$    print *, rb(0:size(rb)-1)
!!$    call flush(6)
    ! Unpack the plan into the outgoing buffer.
    loc = 0  ! We load first from position 0 in the receive buffer.
    do i = 0, col_sz-1
       if( rcnts(i) .ne. 0 ) then ! we expect something from rank 'i'
          if( loc .ne. rdisp(i) ) then
             write (*,'(a,i4)') &
                  'ERROR: discrepancy between rdispi(i) and index for i = ', i
             stop 'unpacking error'
          end if
          ! get the information on the box to receive, get the limits, and 
          ! convert to the local indices.
          sbox = plan%recv_boxes(i)
          loi = get_box_2D_i_min(sbox)
          loj = get_box_2D_j_min(sbox)
          hii = get_box_2D_i_max(sbox)
          hij = get_box_2D_j_max(sbox)
          tmpa(:)  = (/loi,loj/)
          local_lo = global_to_local_2D( final_layout, tmpa )
          tmpa(:)  = (/hii,hij/)
          local_hi = global_to_local_2D( final_layout, tmpa )
          do jd = local_lo(2), local_hi(2)
             do id = local_lo(1), local_hi(1)
                data_out(id,jd) = rb(loc)
                loc             = loc + 1
             end do
          end do
       end if
    end do
  end subroutine apply_remap_2D_complex

#if 0
  ! This function stopped being viable if the transfer function is taken out,
  ! unless a way is found to ship around arbitrary data types passed from a
  ! high-level call.
  subroutine apply_remap_2D_efield( plan, data_in, data_out )
    intrinsic                                 :: transfer
    type(remap_plan_2D), pointer              :: plan
    type(efield_2d_point), dimension(:,:), intent(in)    :: data_in
    type(efield_2d_point), dimension(:,:), intent(out)   :: data_out
    sll_int32, dimension(:), pointer          :: sb     ! send buffer
    sll_int32, dimension(:), pointer          :: rb     ! receive buffer
    sll_int32, dimension(:), pointer          :: sdisp  ! send displacements
    sll_int32, dimension(:), pointer          :: rdisp  ! receive displacements 
    sll_int32, dimension(:), pointer          :: scnts  ! send counts
    sll_int32, dimension(:), pointer          :: rcnts  ! receive counts
    type(sll_collective_t), pointer           :: col    ! collective
    type(layout_2D), pointer                :: init_layout  => NULL()
    type(layout_2D), pointer                :: final_layout => NULL()
    sll_int32                                 :: id, jd
    sll_int32                                 :: i
    sll_int32                                 :: col_sz
    sll_int32                                 :: ierr
    sll_int32                                 :: loi, loj
    sll_int32                                 :: hii, hij
    type(box_2D)                              :: sbox
    sll_int32                                 :: my_rank
    sll_int32                                 :: loc
    sll_int32, dimension(1:2)                 :: local_lo, local_hi
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
    init_layout  => get_remap_2D_initial_layout(plan)
    final_layout => get_remap_2D_final_layout(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer

    SLL_ALLOCATE(sdispi(0:col_sz-1), ierr)
    SLL_ALLOCATE(rdispi(0:col_sz-1), ierr)
    SLL_ALLOCATE(scntsi(0:col_sz-1), ierr)
    SLL_ALLOCATE(rcntsi(0:col_sz-1), ierr)

    ! Translate the amounts into integers
#if 1
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1)), sdisp, &
         col_sz, sdispi)
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1)), rdisp, &
         col_sz, rdispi)
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1)), scnts, &
         col_sz, scntsi)
    call convert_into_integer_sizes(INT32_SIZEOF(data_in(1,1)), rcnts, &
         col_sz, rcntsi)
#endif
    
#if 0
    write (*,'(a,i4)') 'parameters from rank ', my_rank
    print *, 'scntsi', scntsi(:)
    print *, 'sdispi', sdispi(:)
    print *, 'rcntsi', rcntsi(:)
    print *, 'rdispi', rdispi(:)
    call flush(6)
#endif
    
    ! load the send buffer
    loc = 0             ! first loading is at position zero
    ! This step is obviously not needed for integers themselves. We put this
    ! here for generality.
    int32_data_size = INT32_SIZEOF( data_in(1,1) )
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdispi(i) ) then
             print *, 'ERROR DETECTED in process: ', my_rank
             print *, 'apply_remap_2D_double() ERROR: ', &
                  'discrepancy between displs(i) and the loading index for ',&
                  'i = ', i, ' displs(i) = ', sdispi(i)
             write(*,'(a,i8)') 'col_sz = ', col_sz
             call flush(6)
             stop 'apply_remap(): loading error'
          end if
          ! get the information on the box to send, get the limits,
          ! convert to the local indices and find out where in the 
          ! buffer to start writing.
          sbox = plan%send_boxes(i)
          loi = get_box_2D_i_min(sbox)
          loj = get_box_2D_j_min(sbox)
          hii = get_box_2D_i_max(sbox)
          hij = get_box_2D_j_max(sbox)
          local_lo = global_to_local_2D( init_layout, (/loi,loj/) )
          local_hi = global_to_local_2D( init_layout, (/hii,hij/) )

          ! The plan to load the send buffer is to traverse the integer
          ! array with a single index (loc). When we load the buffer, each
          ! data element may occupy multiple integer 'slots', hence the
          ! loading index needs to be manually increased. As an advantage,
          ! we can do some error checking every time we send data to a 
          ! different process, as we know what is the expected value of 
          ! the index at that point.
          do jd = local_lo(2), local_hi(2)
             do id = local_lo(1), local_hi(1)
                sb(loc:) = transfer(data_in(id,jd),(/1_i32/))
                loc      = loc + int32_data_size
             end do
          end do
       end if
    end do
    ! Comment the following when not debugging    
    !   write (*,'(a,i4)') 'the send buffer in rank:', my_rank
    !  print *, sb(0:(size(sb)-1))
    ! call flush(6)
    !    print *, 'from inside remap: rank ', my_rank, 'calling communications'
    !    call flush(6)
   if( plan%is_uniform .eqv. .false. ) then 
      ! the following call can be changed from a generic to a type-specific
      ! call when right away, but especially if the apply_remap function gets
      ! specialized (i.e. gets rid of transfer() calls).
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
!    call flush(6)
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
          loi = get_box_2D_i_min(sbox)
          loj = get_box_2D_j_min(sbox)
          hii = get_box_2D_i_max(sbox)
          hij = get_box_2D_j_max(sbox)
          local_lo = global_to_local_2D( final_layout, (/loi,loj/) )
          local_hi = global_to_local_2D( final_layout, (/hii,hij/) )
          do jd = local_lo(2), local_hi(2)
             do id = local_lo(1), local_hi(1)
                data_out(id,jd) = transfer(rb(loc:),data_out(1,1))
                loc                = loc + int32_data_size
             end do
          end do
       end if
    end do
    ! And why weren't these arrays part of the plan anyway??
    SLL_DEALLOCATE_ARRAY(sdispi, ierr)
    SLL_DEALLOCATE_ARRAY(rdispi, ierr)
    SLL_DEALLOCATE_ARRAY(scntsi, ierr)
    SLL_DEALLOCATE_ARRAY(rcntsi, ierr)
  end subroutine apply_remap_2D_efield
#endif

  subroutine apply_remap_3D_double( plan, data_in, data_out )
    type(remap_plan_3D_real64), pointer              :: plan
    sll_real64, dimension(:,:,:), intent(in)  :: data_in
    sll_real64, dimension(:,:,:), intent(out) :: data_out
    sll_real64, dimension(:), pointer         :: sb     ! send buffer
    sll_real64, dimension(:), pointer         :: rb     ! receive buffer
    sll_int32, dimension(:), pointer          :: sdisp  ! send displacements
    sll_int32, dimension(:), pointer          :: rdisp  ! receive displacements 
    sll_int32, dimension(:), pointer          :: scnts  ! send counts
    sll_int32, dimension(:), pointer          :: rcnts  ! receive counts
    type(sll_collective_t), pointer           :: col    ! collective
    type(layout_3D), pointer                  :: init_layout  => NULL()
    type(layout_3D), pointer                  :: final_layout => NULL()
    sll_int32                                 :: id, jd, kd
    sll_int32                                 :: i
    sll_int32                                 :: col_sz
    sll_int32                                 :: loi, loj, lok
    sll_int32                                 :: hii, hij, hik
    type(box_3D)                              :: sbox
    sll_int32                                 :: my_rank
    sll_int32                                 :: loc
    sll_int32, dimension(1:3)                 :: local_lo, local_hi
    sll_int32, dimension(1:3)                 :: tmpa

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp        => plan%send_displs
    rdisp        => plan%recv_displs
    scnts        => plan%send_counts
    rcnts        => plan%recv_counts
    col          => plan%collective
    col_sz       =  sll_get_collective_size(col)
    init_layout  => get_remap_3D_initial_layout_real64(plan)
    final_layout => get_remap_3D_final_layout_real64(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer

    ! print *, 'from rank ', my_rank, 'loading parameters: ', sdisp, rdisp, &
    ! scnts, rcnts

    ! load the send buffer
    loc = 0             ! first loading is at position zero
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdisp(i) ) then
             print *, 'ERROR DETECTED in process: ', my_rank
             print *, 'apply_remap_3D_double() ERROR: ', &
                  'discrepancy between displs(i) and the loading index for ',&
                  'i = ', i, ' displs(i) = ', sdisp(i)
             write(*,'(a,i8)') 'col_sz = ', col_sz
             call flush(6)
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
          tmpa(:)  = (/loi,loj,lok/)
          local_lo = global_to_local_3D( init_layout, tmpa )
          tmpa(:)  = (/hii,hij,hik/)
          local_hi = global_to_local_3D( init_layout, tmpa )

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
                   sb(loc) = data_in(id,jd,kd)
                   loc     = loc + 1
                end do
             end do
          end do
       end if
    end do
    ! Comment the following when not debugging    
 !   write (*,'(a,i4)') 'the send buffer in rank:', my_rank
  !  print *, sb(0:(size(sb)-1))
   ! call flush(6)
!    print *, 'from inside remap: rank ', my_rank, 'calling communications'
!    call flush(6)
   if( plan%is_uniform .eqv. .false. ) then 
       call sll_collective_alltoallV( sb(:),       &
                                      scnts(0:col_sz-1), &
                                      sdisp(0:col_sz-1), &
                                      rb(:),       &
                                      rcnts(0:col_sz-1), &
                                      rdisp(0:col_sz-1), col )
    else
       call sll_collective_alltoall ( sb(:), &
                                      scnts(0), &
                                      rcnts(0), &
                                      rb(:), col )
    end if
!    write (*,'(a, i4)') 'the receive buffer in rank: ', my_rank
!    print *, rb(0:size(rb)-1)
!    call flush(6)
    ! Unpack the plan into the outgoing buffer.
    loc = 0  ! We load first from position 0 in the receive buffer.
    do i = 0, col_sz-1
       if( rcnts(i) .ne. 0 ) then ! we expect something from rank 'i'
          if( loc .ne. rdisp(i) ) then
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
          tmpa(:)  = (/loi,loj,lok/)
          local_lo = global_to_local_3D( final_layout, tmpa )
          tmpa(:)  = (/hii,hij,hik/)
          local_hi = global_to_local_3D( final_layout, tmpa )
          do kd = local_lo(3), local_hi(3)
             do jd = local_lo(2), local_hi(2)
                do id = local_lo(1), local_hi(1)
                   data_out(id,jd,kd) = rb(loc)
                   loc                = loc + 1
                end do
             end do
          end do
       end if
    end do
  end subroutine apply_remap_3D_double

  subroutine apply_remap_4D_double( plan, data_in, data_out )
    type(remap_plan_4D_real64), pointer                :: plan
    sll_real64, dimension(:,:,:,:), intent(in)  :: data_in
    sll_real64, dimension(:,:,:,:), intent(out) :: data_out
    sll_real64, dimension(:), pointer         :: sb     ! send buffer
    sll_real64, dimension(:), pointer         :: rb     ! receive buffer
    sll_int32, dimension(:), pointer          :: sdisp  ! send displacements
    sll_int32, dimension(:), pointer          :: rdisp  ! receive displacements 
    sll_int32, dimension(:), pointer          :: scnts  ! send counts
    sll_int32, dimension(:), pointer          :: rcnts  ! receive counts
    type(sll_collective_t), pointer           :: col    ! collective
    type(layout_4D), pointer                  :: init_layout  => NULL()
    type(layout_4D), pointer                  :: final_layout => NULL()
    sll_int32                                 :: id, jd, kd, ld
    sll_int32                                 :: i
    sll_int32                                 :: col_sz
    sll_int32                                 :: loi, loj, lok, lol
    sll_int32                                 :: hii, hij, hik, hil
    type(box_4D)                              :: sbox
    sll_int32                                 :: my_rank
    sll_int32                                 :: loc
    sll_int32, dimension(1:4)                 :: local_lo, local_hi
    sll_int32, dimension(1:4)                 :: tmpa

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp        => plan%send_displs
    rdisp        => plan%recv_displs
    scnts        => plan%send_counts
    rcnts        => plan%recv_counts
    col          => plan%collective
    col_sz       =  sll_get_collective_size(col)
    init_layout  => get_remap_4D_initial_layout_real64(plan)
    final_layout => get_remap_4D_final_layout_real64(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer

    ! print *, 'from rank ', my_rank, 'loading parameters: ', sdisp, rdisp, &
    ! scnts, rcnts

    ! load the send buffer
    loc = 0             ! first loading is at position zero
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdisp(i) ) then
             print *, 'ERROR DETECTED in process: ', my_rank
             print *, 'apply_remap_4D_double() ERROR: ', &
                  'discrepancy between displs(i) and the loading index for ',&
                  'i = ', i, ' displs(i) = ', sdisp(i)
             write(*,'(a,i8)') 'col_sz = ', col_sz
             call flush(6)
             stop 'apply_remap(): loading error'
          end if
          ! get the information on the box to send, get the limits,
          ! convert to the local indices and find out where in the 
          ! buffer to start writing.
          sbox = plan%send_boxes(i)
          loi = get_box_4D_i_min(sbox)
          loj = get_box_4D_j_min(sbox)
          lok = get_box_4D_k_min(sbox)
          lol = get_box_4D_l_min(sbox)
          hii = get_box_4D_i_max(sbox)
          hij = get_box_4D_j_max(sbox)
          hik = get_box_4D_k_max(sbox)
          hil = get_box_4D_l_max(sbox)
          tmpa(:)  = (/loi,loj,lok,lol/)
          local_lo = global_to_local_4D( init_layout, tmpa )
          tmpa(:)  = (/hii,hij,hik,hil/)
          local_hi = global_to_local_4D( init_layout, tmpa )

          ! The plan to load the send buffer is to traverse the integer
          ! array with a single index (loc). When we load the buffer, each
          ! data element may occupy multiple integer 'slots', hence the
          ! loading index needs to be manually increased. As an advantage,
          ! we can do some error checking every time we send data to a 
          ! different process, as we know what is the expected value of 
          ! the index at that point.
          do ld = local_lo(4), local_hi(4)
             do kd = local_lo(3), local_hi(3)
                do jd = local_lo(2), local_hi(2)
                   do id = local_lo(1), local_hi(1)
                      sb(loc) = data_in(id,jd,kd,ld)
                      loc     = loc + 1
                   end do
                end do
             end do
          end do
       end if
    end do
    ! Comment the following when not debugging    
    !   write (*,'(a,i4)') 'the send buffer in rank:', my_rank
    !  print *, sb(0:(size(sb)-1))
    ! call flush(6)
    !    print *, 'from inside remap: rank ', my_rank, 'calling communications'
    !    call flush(6)
    if( plan%is_uniform .eqv. .false. ) then 
       call sll_collective_alltoallV( sb(:),       &
                                      scnts(0:col_sz-1), &
                                      sdisp(0:col_sz-1), &
                                      rb(:),       &
                                      rcnts(0:col_sz-1), &
                                      rdisp(0:col_sz-1), col )
    else
       call sll_collective_alltoall ( sb(:), &
                                      scnts(0), &
                                      rcnts(0), &
                                      rb(:), col )
    end if
    !    write (*,'(a, i4)') 'the receive buffer in rank: ', my_rank
    !    print *, rb(0:size(rb)-1)
    !    call flush(6)
    ! Unpack the plan into the outgoing buffer.
    loc = 0  ! We load first from position 0 in the receive buffer.
    do i = 0, col_sz-1
       if( rcnts(i) .ne. 0 ) then ! we expect something from rank 'i'
          if( loc .ne. rdisp(i) ) then
             write (*,'(a,i4)') &
                  'ERROR: discrepancy between rdispi(i) and index for i = ', i
             stop 'unpacking error'
          end if
          ! get the information on the box to receive, get the limits, and 
          ! convert to the local indices.
          sbox = plan%recv_boxes(i)
          loi = get_box_4D_i_min(sbox)
          loj = get_box_4D_j_min(sbox)
          lok = get_box_4D_k_min(sbox)
          lol = get_box_4D_l_min(sbox)
          hii = get_box_4D_i_max(sbox)
          hij = get_box_4D_j_max(sbox)
          hik = get_box_4D_k_max(sbox)
          hil = get_box_4D_l_max(sbox)
          tmpa(:)  = (/loi,loj,lok,lol/)
          local_lo = global_to_local_4D( final_layout, tmpa )
          tmpa(:)  = (/hii,hij,hik,hil/)
          local_hi = global_to_local_4D( final_layout, tmpa )
          do ld = local_lo(4), local_hi(4)
             do kd = local_lo(3), local_hi(3)
                do jd = local_lo(2), local_hi(2)
                   do id = local_lo(1), local_hi(1)
                      data_out(id,jd,kd,ld) = rb(loc)
                      loc                   = loc + 1
                   end do
                end do
             end do
          end do
       end if
    end do
  end subroutine apply_remap_4D_double

  subroutine apply_remap_3D_complex( plan, data_in, data_out )
    type(remap_plan_3D_comp64), pointer              :: plan
    sll_comp64, dimension(:,:,:), intent(in)  :: data_in
    sll_comp64, dimension(:,:,:), intent(out) :: data_out
    sll_comp64, dimension(:), pointer          :: sb     ! send buffer
    sll_comp64, dimension(:), pointer          :: rb     ! receive buffer
    sll_int32, dimension(:), pointer          :: sdisp  ! send displacements
    sll_int32, dimension(:), pointer          :: rdisp  ! receive displacements 
    sll_int32, dimension(:), pointer          :: scnts  ! send counts
    sll_int32, dimension(:), pointer          :: rcnts  ! receive counts
    type(sll_collective_t), pointer           :: col    ! collective
    type(layout_3D), pointer                  :: init_layout  => NULL()
    type(layout_3D), pointer                  :: final_layout => NULL()
    sll_int32                                 :: id, jd, kd
    sll_int32                                 :: i
    sll_int32                                 :: col_sz
    sll_int32                                 :: loi, loj, lok
    sll_int32                                 :: hii, hij, hik
    type(box_3D)                              :: sbox
    sll_int32                                 :: my_rank
    sll_int32                                 :: loc
    sll_int32, dimension(1:3)                 :: local_lo, local_hi
    sll_int32, dimension(1:3)                 :: tmpa

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp        => plan%send_displs
    rdisp        => plan%recv_displs
    scnts        => plan%send_counts
    rcnts        => plan%recv_counts
    col          => plan%collective
    col_sz       =  sll_get_collective_size(col)
    init_layout  => get_remap_3D_initial_layout_comp64(plan)
    final_layout => get_remap_3D_final_layout_comp64(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer

    ! load the send buffer
    loc = 0             ! first loading is at position zero
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdisp(i) ) then
             write (*,'(a,i4,a,i16)') 'apply_remap_3D_int() ERROR: discrepancy between displs(i) and the loading index for i = ', i, ' displs(i) = ', sdisp(i)
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
          tmpa(:)  = (/loi,loj,lok/)
          local_lo = global_to_local_3D( init_layout, tmpa )
          tmpa(:)  = (/hii,hij,hik/)
          local_hi = global_to_local_3D( init_layout, tmpa )

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
                   sb(loc) = data_in(id,jd,kd)
                   loc     = loc + 1
                end do
             end do
          end do
       end if
    end do
    
!    write (*,'(a,i4)') 'the send buffer in rank:', my_rank
!    print *, sb(0:(size(sb)-1))
!    call flush(6)
 
   if( plan%is_uniform .eqv. .false. ) then 
       call sll_collective_alltoallV( sb(:),       &
                                      scnts(0:col_sz-1), &
                                      sdisp(0:col_sz-1), &
                                      rb(:),       &
                                      rcnts(0:col_sz-1), &
                                      rdisp(0:col_sz-1), col )
    else
       call sll_collective_alltoall ( sb(:), &
                                      scnts(0), &
                                      rcnts(0), &
                                      rb(:), col )
    end if
!    write (*,'(a, i4)') 'the receive buffer in rank: ', my_rank
!    print *, rb(0:size(rb)-1)
!    call flush(6)
    ! Unpack the plan into the outgoing buffer.
    loc = 0  ! We load first from position 0 in the receive buffer.
    do i = 0, col_sz-1
       if( rcnts(i) .ne. 0 ) then ! we expect something from rank 'i'
          if( loc .ne. rdisp(i) ) then
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
          tmpa(:)  = (/loi,loj,lok/)
          local_lo = global_to_local_3D( final_layout, tmpa )
          tmpa(:)  = (/hii,hij,hik/)
          local_hi = global_to_local_3D( final_layout, tmpa )
          do kd = local_lo(3), local_hi(3)
             do jd = local_lo(2), local_hi(2)
                do id = local_lo(1), local_hi(1)
                   data_out(id,jd,kd) = rb(loc)
                   loc                = loc + 1
                end do
             end do
          end do
       end if
    end do
  end subroutine apply_remap_3D_complex

  subroutine apply_remap_6D_double( plan, data_in, data_out )
    type(remap_plan_6D_real64), pointer              :: plan
    sll_real64, dimension(:,:,:,:,:,:), intent(in)  :: data_in
    sll_real64, dimension(:,:,:,:,:,:), intent(out) :: data_out
    sll_real64, dimension(:), pointer          :: sb     ! send buffer
    sll_real64, dimension(:), pointer          :: rb     ! receive buffer
    sll_int32, dimension(:), pointer          :: sdisp  ! send displacements
    sll_int32, dimension(:), pointer          :: rdisp  ! receive displacements 
    sll_int32, dimension(:), pointer          :: scnts  ! send counts
    sll_int32, dimension(:), pointer          :: rcnts  ! receive counts
    type(sll_collective_t), pointer           :: col    ! collective
    type(layout_6D), pointer                  :: init_layout  => NULL()
    type(layout_6D), pointer                  :: final_layout => NULL()
    sll_int32                                 :: id, jd, kd, ld, md, nd
    sll_int32                                 :: i
    sll_int32                                 :: col_sz
    sll_int32                                 :: loi, loj, lok, lol, lom, lon
    sll_int32                                 :: hii, hij, hik, hil, him, hin
    type(box_6D)                              :: sbox
    sll_int32                                 :: my_rank
    sll_int32                                 :: loc
    sll_int32, dimension(1:6)                 :: local_lo, local_hi
    sll_int32, dimension(1:6)                 :: tmpa

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp        => plan%send_displs
    rdisp        => plan%recv_displs
    scnts        => plan%send_counts
    rcnts        => plan%recv_counts
    col          => plan%collective
    col_sz       =  sll_get_collective_size(col)
    init_layout  => get_remap_6D_initial_layout_real64(plan)
    final_layout => get_remap_6D_final_layout_real64(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer

    ! print *, 'from rank ', my_rank, 'loading parameters: ', sdisp, rdisp, &
    ! scnts, rcnts

    ! load the send buffer
    loc = 0             ! first loading is at position zero
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdisp(i) ) then
             print *, 'ERROR DETECTED in process: ', my_rank
             print *, 'apply_remap_6D_double() ERROR: ', &
                  'discrepancy between displs(i) and the loading index for ',&
                  'i = ', i, ' displs(i) = ', sdisp(i)
             write(*,'(a,i8)') 'col_sz = ', col_sz
             call flush(6)
             stop 'apply_remap(): exchange buffer loading error'
          end if
          ! get the information on the box to send, get the limits,
          ! convert to the local indices and find out where in the 
          ! buffer to start writing.
          sbox = plan%send_boxes(i)
          loi = get_box_6D_i_min(sbox)
          loj = get_box_6D_j_min(sbox)
          lok = get_box_6D_k_min(sbox)
          lol = get_box_6D_l_min(sbox)
          lom = get_box_6D_m_min(sbox)
          lon = get_box_6D_n_min(sbox)
          hii = get_box_6D_i_max(sbox)
          hij = get_box_6D_j_max(sbox)
          hik = get_box_6D_k_max(sbox)
          hil = get_box_6D_l_max(sbox)
          him = get_box_6D_m_max(sbox)
          hin = get_box_6D_n_max(sbox)
          tmpa(:)  = (/loi,loj,lok,lol,lom,lon/)
          local_lo = global_to_local_6D( init_layout, tmpa )
          tmpa(:)  = (/hii,hij,hik,hil,him,hin/)
          local_hi = global_to_local_6D( init_layout, tmpa)

          ! The plan to load the send buffer is to traverse the send buffer
          ! array with a single index (loc). We manually increment the loading
          ! index. As an advantage,
          ! we can do some error checking every time we send data to a 
          ! different process, as we know what is the expected value of 
          ! the index at that point.
          do nd = local_lo(6), local_hi(6)
             do md = local_lo(5), local_hi(5)
                do ld = local_lo(4), local_hi(4)
                   do kd = local_lo(3), local_hi(3)
                      do jd = local_lo(2), local_hi(2)
                         do id = local_lo(1), local_hi(1)
                            sb(loc) = data_in(id,jd,kd,ld,md,nd)
                            loc     = loc + 1
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end if
    end do
    ! Comment the following when not debugging    
!!$    write (*,'(a,i4)') 'the send buffer in rank:', my_rank
!!$    print *, sb(0:(size(sb)-1))
!!$    call flush(6)
!!$    print *, 'from inside remap: rank ', my_rank, 'calling communications'
!!$    call flush(6)

    if( plan%is_uniform .eqv. .false. ) then 
       call sll_collective_alltoallV( sb(:),       &
                                      scnts(0:col_sz-1), &
                                      sdisp(0:col_sz-1), &
                                      rb(:),       &
                                      rcnts(0:col_sz-1), &
                                      rdisp(0:col_sz-1), col )
    else
       call sll_collective_alltoall ( sb(:), &
                                      scnts(0), &
                                      rcnts(0), &
                                      rb(:), col )
    end if

!!$    write (*,'(a, i4)') 'receive buffer in rank: ', my_rank
!!$    print *, rb(0:size(rb)-1)
!!$    call flush(6)
    ! Unpack the plan into the outgoing buffer.
    loc = 0  ! We load first from position 0 in the receive buffer.
    do i = 0, col_sz-1
       if( rcnts(i) .ne. 0 ) then ! we expect something from rank 'i'
          if( loc .ne. rdisp(i) ) then
             write (*,'(a,i4)') &
                  'ERROR: discrepancy between rdispi(i) and index for i = ', i
             stop 'exchange buffer unpacking error'
          end if
          ! get the information on the box to receive, get the limits, and 
          ! convert to the local indices.
          sbox = plan%recv_boxes(i)
          loi = get_box_6D_i_min(sbox)
          loj = get_box_6D_j_min(sbox)
          lok = get_box_6D_k_min(sbox)
          lol = get_box_6D_l_min(sbox)
          lom = get_box_6D_m_min(sbox)
          lon = get_box_6D_n_min(sbox)
          hii = get_box_6D_i_max(sbox)
          hij = get_box_6D_j_max(sbox)
          hik = get_box_6D_k_max(sbox)
          hil = get_box_6D_l_max(sbox)
          him = get_box_6D_m_max(sbox)
          hin = get_box_6D_n_max(sbox)
          tmpa(:)  = (/loi,loj,lok,lol,lom,lon/)
          local_lo = global_to_local_6D(final_layout, tmpa )
          tmpa(:)  = (/hii,hij,hik,hil,him,hin/)
          local_hi = global_to_local_6D(final_layout, tmpa)
          do nd = local_lo(6), local_hi(6)
             do md = local_lo(5), local_hi(5)
                do ld = local_lo(4), local_hi(4)
                   do kd = local_lo(3), local_hi(3)
                      do jd = local_lo(2), local_hi(2)
                         do id = local_lo(1), local_hi(1)
                            data_out(id,jd,kd,ld,md,nd) = rb(loc)
                            loc                         = loc + 1
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end if
    end do
  end subroutine apply_remap_6D_double

  subroutine apply_remap_6D_int( plan, data_in, data_out )
    type(remap_plan_6D_int32), pointer              :: plan
    sll_int32, dimension(:,:,:,:,:,:), intent(in)  :: data_in
    sll_int32, dimension(:,:,:,:,:,:), intent(out) :: data_out
    sll_int32, dimension(:), pointer          :: sb     ! send buffer
    sll_int32, dimension(:), pointer          :: rb     ! receive buffer
    sll_int32, dimension(:), pointer          :: sdisp  ! send displacements
    sll_int32, dimension(:), pointer          :: rdisp  ! receive displacements 
    sll_int32, dimension(:), pointer          :: scnts  ! send counts
    sll_int32, dimension(:), pointer          :: rcnts  ! receive counts
    type(sll_collective_t), pointer           :: col    ! collective
    type(layout_6D), pointer                  :: init_layout  => NULL()
    type(layout_6D), pointer                  :: final_layout => NULL()
    sll_int32                                 :: id, jd, kd, ld, md, nd
    sll_int32                                 :: i
    sll_int32                                 :: col_sz
    sll_int32                                 :: loi, loj, lok, lol, lom, lon
    sll_int32                                 :: hii, hij, hik, hil, him, hin
    type(box_6D)                              :: sbox
    sll_int32                                 :: my_rank
    sll_int32                                 :: loc
    sll_int32, dimension(1:6)                 :: local_lo, local_hi

    ! unpack the plan: There are inconsistencies here, one one hand we access
    ! directly and on the other with access functions... standardize...
    sdisp        => plan%send_displs
    rdisp        => plan%recv_displs
    scnts        => plan%send_counts
    rcnts        => plan%recv_counts
    col          => plan%collective
    col_sz       =  sll_get_collective_size(col)
    init_layout  => get_remap_6D_initial_layout_int32(plan)
    final_layout => get_remap_6D_final_layout_int32(plan)
    my_rank      =  sll_get_collective_rank(col)
    sb           => plan%send_buffer
    rb           => plan%recv_buffer

    ! load the send buffer
    loc = 0             ! first loading is at position zero
    do i = 0, col_sz-1
       if( scnts(i) .ne. 0 ) then ! send something to rank 'i'
          if( loc .ne. sdisp(i) ) then
             print *, 'ERROR DETECTED in process: ', my_rank
             print *, 'apply_remap_6D_double() ERROR: ', &
                  'discrepancy between displs(i) and the loading index for ',&
                  'i = ', i, ' displs(i) = ', sdisp(i)
             write(*,'(a,i8)') 'col_sz = ', col_sz
             call flush(6)
             stop 'apply_remap(): exchange buffer loading error'
          end if
          ! get the information on the box to send, get the limits,
          ! convert to the local indices and find out where in the 
          ! buffer to start writing.
          sbox = plan%send_boxes(i)
          loi = get_box_6D_i_min(sbox)
          loj = get_box_6D_j_min(sbox)
          lok = get_box_6D_k_min(sbox)
          lol = get_box_6D_l_min(sbox)
          lom = get_box_6D_m_min(sbox)
          lon = get_box_6D_n_min(sbox)
          hii = get_box_6D_i_max(sbox)
          hij = get_box_6D_j_max(sbox)
          hik = get_box_6D_k_max(sbox)
          hil = get_box_6D_l_max(sbox)
          him = get_box_6D_m_max(sbox)
          hin = get_box_6D_n_max(sbox)
          local_lo = &
               global_to_local_6D( init_layout,(/loi,loj,lok,lol,lom,lon/))
          local_hi = &
               global_to_local_6D( init_layout,(/hii,hij,hik,hil,him,hin/))
          ! The plan to load the send buffer is to traverse the integer
          ! array with a single index (loc). When we load the buffer, each
          ! data element may occupy multiple integer 'slots', hence the
          ! loading index needs to be manually increased. As an advantage,
          ! we can do some error checking every time we send data to a 
          ! different process, as we know what is the expected value of 
          ! the index at that point.
          do nd = local_lo(6), local_hi(6)
             do md = local_lo(5), local_hi(5)
                do ld = local_lo(4), local_hi(4)
                   do kd = local_lo(3), local_hi(3)
                      do jd = local_lo(2), local_hi(2)
                         do id = local_lo(1), local_hi(1)
                            sb(loc) = data_in(id,jd,kd,ld,md,nd)
                            loc      = loc + 1
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end if
    end do

    ! Comment the following when not debugging    
!!$    if(my_rank == 0) then
!!$       write (*,'(a,i4)') 'send buffer, rank:', my_rank, 'buffer size = ', &
!!$            size(sb)
!!$       print *, 'sb: ', sb(0:(size(sb)-1))
!!$       print *, 'scntsi: ', scntsi(:)
!!$       print *, 'sdispi: ', sdispi(:)
!!$       print *, 'rb: ', rb(:)
!!$       print *, 'rcntsi: ', rcntsi(:)
!!$       print *, 'rdispi: ', rdispi(:)
!!$       print *, 'uniformity: ', plan%is_uniform
!!$       call flush(6)
!!$    end if
!!$    print *, 'from inside remap: rank ', my_rank, 'calling communications'
!!$    call flush(6)

    if( plan%is_uniform .eqv. .false. ) then 
       call sll_collective_alltoallV( sb(:),       &
                                      scnts(0:col_sz-1), &
                                      sdisp(0:col_sz-1), &
                                      rb(:),       &
                                      rcnts(0:col_sz-1), &
                                      rdisp(0:col_sz-1), col )
    else
       call sll_collective_alltoall ( sb(:), &
                                      scnts(0), &
                                      rcnts(0), &
                                      rb(:), col )
    end if

!!$    write (*,'(a, i4)') 'receive buffer in rank: ', my_rank
!!$    print *, rb(0:size(rb)-1)
!!$    call flush(6)

    ! Unpack the plan into the outgoing buffer.
    loc = 0  ! We load first from position 0 in the receive buffer.
    do i = 0, col_sz-1
       if( rcnts(i) .ne. 0 ) then ! we expect something from rank 'i'
          if( loc .ne. rdisp(i) ) then
             write (*,'(a,i4)') &
                  'ERROR: discrepancy between rdispi(i) and index for i = ', i
             stop 'exchange buffer unpacking error'
          end if
          ! get the information on the box to receive, get the limits, and 
          ! convert to the local indices.
          sbox = plan%recv_boxes(i)
          loi = get_box_6D_i_min(sbox)
          loj = get_box_6D_j_min(sbox)
          lok = get_box_6D_k_min(sbox)
          lol = get_box_6D_l_min(sbox)
          lom = get_box_6D_m_min(sbox)
          lon = get_box_6D_n_min(sbox)
          hii = get_box_6D_i_max(sbox)
          hij = get_box_6D_j_max(sbox)
          hik = get_box_6D_k_max(sbox)
          hil = get_box_6D_l_max(sbox)
          him = get_box_6D_m_max(sbox)
          hin = get_box_6D_n_max(sbox)
          local_lo = &
               global_to_local_6D(final_layout, (/loi,loj,lok,lol,lom,lon/))
          local_hi = &
               global_to_local_6D(final_layout, (/hii,hij,hik,hil,him,hin/))
          do nd = local_lo(6), local_hi(6)
             do md = local_lo(5), local_hi(5)
                do ld = local_lo(4), local_hi(4)
                   do kd = local_lo(3), local_hi(3)
                      do jd = local_lo(2), local_hi(2)
                         do id = local_lo(1), local_hi(1)
                            data_out(id,jd,kd,ld,md,nd) = rb(loc)
                            loc                         = loc + 1
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end if
    end do
  end subroutine apply_remap_6D_int

  ! Placeholder: the load/unload subroutines for the send and receive buffers
  ! should ideally be abstracted out. This means that we need to probably
  ! define some generic interface and hide behind the types that we want to
  ! exchange. We can think of direct loading and exchanging the basic types:
  ! - single/double precision floats
  ! - integers
  ! - stay with the transfer function for other types.
#if 0
  subroutine load_send_buffer_3D( plan, data_in )
    type(remap_plan_3D), pointer :: plan
    sll_int32                      :: load_point
  end subroutine load_send_buffer_3D
#endif

#define MAKE_GET_BOX_SLOT_FUNCTION( fname, boxtype, slot )   \
  function fname( b );                                       \
    sll_int32                 :: fname;                      \
    type(boxtype), intent(in) :: b;                          \
    fname = b%slot;                                          \
  end function fname

  MAKE_GET_BOX_SLOT_FUNCTION( get_box_2D_i_min, box_2D, i_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_2D_i_max, box_2D, i_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_2D_j_min, box_2D, j_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_2D_j_max, box_2D, j_max )

  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_i_min, box_3D, i_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_i_max, box_3D, i_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_j_min, box_3D, j_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_j_max, box_3D, j_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_k_min, box_3D, k_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_3D_k_max, box_3D, k_max )

  MAKE_GET_BOX_SLOT_FUNCTION( get_box_4D_i_min, box_4D, i_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_4D_i_max, box_4D, i_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_4D_j_min, box_4D, j_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_4D_j_max, box_4D, j_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_4D_k_min, box_4D, k_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_4D_k_max, box_4D, k_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_4D_l_min, box_4D, l_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_4D_l_max, box_4D, l_max )

  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_i_min, box_6D, i_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_i_max, box_6D, i_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_j_min, box_6D, j_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_j_max, box_6D, j_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_k_min, box_6D, k_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_k_max, box_6D, k_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_l_min, box_6D, l_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_l_max, box_6D, l_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_m_min, box_6D, m_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_m_max, box_6D, m_max )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_n_min, box_6D, n_min )
  MAKE_GET_BOX_SLOT_FUNCTION( get_box_6D_n_max, box_6D, n_max )


  function sll_get_num_nodes_2D( lims )
    sll_int32                  :: sll_get_num_nodes_2D
    type(layout_2D), pointer :: lims
    sll_get_num_nodes_2D = sll_get_collective_size( lims%collective )
  end function sll_get_num_nodes_2D

  function sll_get_num_nodes_3D( lims )
    sll_int32                  :: sll_get_num_nodes_3D
    type(layout_3D), pointer :: lims
    sll_get_num_nodes_3D = sll_get_collective_size( lims%collective )
  end function sll_get_num_nodes_3D

  function sll_get_num_nodes_4D( lims )
    sll_int32                  :: sll_get_num_nodes_4D
    type(layout_4D), pointer :: lims
    sll_get_num_nodes_4D = sll_get_collective_size( lims%collective )
  end function sll_get_num_nodes_4D

  function sll_get_num_nodes_6D( lims )
    sll_int32                  :: sll_get_num_nodes_6D
    type(layout_6D), pointer :: lims
    sll_get_num_nodes_6D = sll_get_collective_size( lims%collective )
  end function sll_get_num_nodes_6D


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
  function local_to_global_2D( layout, doublet )
    sll_int32, dimension(1:2)             :: local_to_global_2D
    type(layout_2D), pointer              :: layout
    sll_int32, intent(in), dimension(1:2) :: doublet
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: my_rank
    type(box_2D)                          :: box
    ! fixme: arg checking
    col                => get_layout_2D_collective( layout )
    my_rank            =  sll_get_collective_rank( col )
    box                =  get_layout_2D_box( layout, my_rank )
    local_to_global_2D(1) = get_box_2D_i_min(box) + doublet(1) - 1
    local_to_global_2D(2) = get_box_2D_j_min(box) + doublet(2) - 1
  end function local_to_global_2D

  function local_to_global_3D( layout, triplet )
    sll_int32, dimension(1:3)             :: local_to_global_3D
    type(layout_3D), pointer            :: layout
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

  function local_to_global_4D( layout, quad )
    sll_int32, dimension(1:4)             :: local_to_global_4D
    type(layout_4D), pointer              :: layout
    sll_int32, intent(in), dimension(1:4) :: quad
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: my_rank
    type(box_4D)                          :: box
    ! fixme: arg checking
    col                => get_layout_4D_collective( layout )
    my_rank            =  sll_get_collective_rank( col )
    box                =  get_layout_4D_box( layout, my_rank )
    local_to_global_4D(1) = get_box_4D_i_min(box) + quad(1) - 1
    local_to_global_4D(2) = get_box_4D_j_min(box) + quad(2) - 1
    local_to_global_4D(3) = get_box_4D_k_min(box) + quad(3) - 1
    local_to_global_4D(4) = get_box_4D_l_min(box) + quad(4) - 1
  end function local_to_global_4D

  function local_to_global_6D( layout, hextet )
    sll_int32, dimension(1:6)             :: local_to_global_6D
    type(layout_6D), pointer              :: layout
    sll_int32, intent(in), dimension(1:6) :: hextet
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: my_rank
    type(box_6D)                          :: box
    ! fixme: arg checking
    col                => get_layout_6D_collective( layout )
    my_rank            =  sll_get_collective_rank( col )
    box                =  get_layout_6D_box( layout, my_rank )
    local_to_global_6D(1) = get_box_6D_i_min(box) + hextet(1) - 1
    local_to_global_6D(2) = get_box_6D_j_min(box) + hextet(2) - 1
    local_to_global_6D(3) = get_box_6D_k_min(box) + hextet(3) - 1
    local_to_global_6D(4) = get_box_6D_l_min(box) + hextet(4) - 1
    local_to_global_6D(5) = get_box_6D_m_min(box) + hextet(5) - 1
    local_to_global_6D(6) = get_box_6D_n_min(box) + hextet(6) - 1
  end function local_to_global_6D


  ! We need to make sure that the decision of choosing '-1' as the return
  ! value when the global index is not available locally does not backfire.
  ! If one decides to use an array with an indexing that contains -1, this 
  ! would be problematic.

  function global_to_local_2D( layout, gtuple )
    intrinsic                             :: associated
    sll_int32, dimension(1:2)             :: global_to_local_2D
    type(layout_2D), pointer            :: layout
    sll_int32, dimension(1:2), intent(in) :: gtuple ! global indices, as array
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: my_rank
    type(box_2D)                          :: box
    if( .not. associated(get_layout_2D_collective(layout)) ) then
       write (*,'(a)') 'ERROR in global_to_local_2D(), not-associated col'
       stop 'global_to_local_2D'
    end if
    col     => get_layout_2D_collective( layout )
    my_rank =  sll_get_collective_rank( col )
    box     =  get_layout_2D_box( layout, my_rank )
    if( (gtuple(1) .ge. get_box_2D_i_min(box)) .and. &
        (gtuple(1) .le. get_box_2D_i_max(box)) .and. &
        (gtuple(2) .ge. get_box_2D_j_min(box)) .and. &
        (gtuple(2) .le. get_box_2D_j_max(box)) ) then  ! the index is present
       global_to_local_2D(1) = gtuple(1) - get_box_2D_i_min(box) + 1
       global_to_local_2D(2) = gtuple(2) - get_box_2D_j_min(box) + 1
    else  ! the index is not present
       global_to_local_2D(1) = -1
       global_to_local_2D(2) = -1
    end if
  end function global_to_local_2D

  function global_to_local_3D( layout, gtuple )
    intrinsic                             :: associated
    sll_int32, dimension(1:3)             :: global_to_local_3D
    type(layout_3D), pointer            :: layout
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

  function global_to_local_4D( layout, gtuple )
    intrinsic                             :: associated
    sll_int32, dimension(1:4)             :: global_to_local_4D
    type(layout_4D), pointer            :: layout
    sll_int32, dimension(1:4), intent(in) :: gtuple ! global indices, as array
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: my_rank
    type(box_4D)                          :: box
    if( .not. associated(get_layout_4D_collective(layout)) ) then
       write (*,'(a)') 'ERROR in global_to_local_4D(), not-associated col'
       stop 'global_to_local_4D'
    end if
    col     => get_layout_4D_collective( layout )
    my_rank =  sll_get_collective_rank( col )
    box     =  get_layout_4D_box( layout, my_rank )
    if( (gtuple(1) .ge. get_box_4D_i_min(box)) .and. &
        (gtuple(1) .le. get_box_4D_i_max(box)) .and. &
        (gtuple(2) .ge. get_box_4D_j_min(box)) .and. &
        (gtuple(2) .le. get_box_4D_j_max(box)) .and. &
        (gtuple(3) .ge. get_box_4D_k_min(box)) .and. &
        (gtuple(3) .le. get_box_4D_k_max(box)) .and. &
        (gtuple(4) .ge. get_box_4D_l_min(box)) .and. &
        (gtuple(4) .le. get_box_4D_l_max(box)) ) then  ! the index is present
       global_to_local_4D(1) = gtuple(1) - get_box_4D_i_min(box) + 1
       global_to_local_4D(2) = gtuple(2) - get_box_4D_j_min(box) + 1
       global_to_local_4D(3) = gtuple(3) - get_box_4D_k_min(box) + 1
       global_to_local_4D(4) = gtuple(4) - get_box_4D_l_min(box) + 1
    else  ! the index is not present
       global_to_local_4D(1) = -1
       global_to_local_4D(2) = -1
       global_to_local_4D(3) = -1
       global_to_local_4D(4) = -1
    end if
  end function global_to_local_4D

  function global_to_local_6D( layout, gtuple )
    intrinsic                             :: associated
    sll_int32, dimension(1:6)             :: global_to_local_6D
    type(layout_6D), pointer              :: layout
    sll_int32, dimension(1:6), intent(in) :: gtuple ! global indices, as array
    type(sll_collective_t), pointer       :: col
    sll_int32                             :: my_rank
    type(box_6D)                          :: box
    if( .not. associated(get_layout_6D_collective(layout)) ) then
       write (*,'(a)') 'ERROR in global_to_local_6D(), not-associated col'
       stop 'global_to_local_6D'
    end if
    col     => get_layout_6D_collective( layout )
    my_rank =  sll_get_collective_rank( col )
    box     =  get_layout_6D_box( layout, my_rank )
    if( (gtuple(1) .ge. get_box_6D_i_min(box)) .and. &
        (gtuple(1) .le. get_box_6D_i_max(box)) .and. &
        (gtuple(2) .ge. get_box_6D_j_min(box)) .and. &
        (gtuple(2) .le. get_box_6D_j_max(box)) .and. &
        (gtuple(3) .ge. get_box_6D_k_min(box)) .and. &
        (gtuple(3) .le. get_box_6D_k_max(box)) .and. &
        (gtuple(4) .ge. get_box_6D_l_min(box)) .and. &
        (gtuple(4) .le. get_box_6D_l_max(box)) .and. &
        (gtuple(5) .ge. get_box_6D_m_min(box)) .and. &
        (gtuple(5) .le. get_box_6D_m_max(box)) .and. &
        (gtuple(6) .ge. get_box_6D_n_min(box)) .and. &
        (gtuple(6) .le. get_box_6D_n_max(box)) ) then  ! the index is present
       global_to_local_6D(1) = gtuple(1) - get_box_6D_i_min(box) + 1
       global_to_local_6D(2) = gtuple(2) - get_box_6D_j_min(box) + 1
       global_to_local_6D(3) = gtuple(3) - get_box_6D_k_min(box) + 1
       global_to_local_6D(4) = gtuple(4) - get_box_6D_l_min(box) + 1
       global_to_local_6D(5) = gtuple(5) - get_box_6D_m_min(box) + 1
       global_to_local_6D(6) = gtuple(6) - get_box_6D_n_min(box) + 1
    else  ! the index is not present
       print *, 'WARNING: from rank ', my_rank, 'index is not present.'
       call view_box_6D( box )
       print *, 'passed global indices: ', gtuple(:)
       global_to_local_6D(1) = -1
       global_to_local_6D(2) = -1
       global_to_local_6D(3) = -1
       global_to_local_6D(4) = -1
       global_to_local_6D(5) = -1
       global_to_local_6D(6) = -1
    end if
  end function global_to_local_6D


  subroutine view_box_2D( b )
    type(box_2D), intent(in) :: b
    write(*,'(a,i4,a,i4,a,i4,a,i4,a)') &
         '[  [', b%i_min,',', b%i_max,'], [', &
                 b%j_min,',', b%j_max,']  ]'
  end subroutine view_box_2D

  subroutine view_box_3D( b )
    type(box_3D), intent(in) :: b
    write(*,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') &
         '[  [', b%i_min,',', b%i_max,'], [', &
                 b%j_min,',', b%j_max,'], [', &
                 b%k_min,',', b%k_max,']  ]'
  end subroutine view_box_3D

  subroutine view_box_4D( b )
    type(box_4D), intent(in) :: b
    write(*,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') &
         '[  [', b%i_min,',', b%i_max,'], [', &
                 b%j_min,',', b%j_max,'], [', &
                 b%k_min,',', b%k_max,'], [', &
                 b%l_min,',', b%l_max,'] ]'
  end subroutine view_box_4D

  subroutine view_box_6D( b )
    type(box_6D), intent(in) :: b
    write(*,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') &
         '[  [', b%i_min,',', b%i_max,'], [', &
                 b%j_min,',', b%j_max,'], [', &
                 b%k_min,',', b%k_max,'], [', &
                 b%l_min,',', b%l_max,'], [', &
                 b%m_min,',', b%m_max,'], [', &
                 b%n_min,',', b%n_max,'] ]'
  end subroutine view_box_6D


  subroutine sll_view_lims_2D( layout )
    type(layout_2D), pointer :: layout
    sll_int32                  :: i
    sll_int32                  :: sz
    sz = sll_get_num_nodes( layout )
    print *, 'limits: '
    do i=0,sz-1
       call view_box_2D(get_layout_2D_box( layout, i ))
    end do
    call flush(6)
  end subroutine sll_view_lims_2D


  subroutine sll_view_lims_3D( layout )
    type(layout_3D), pointer :: layout
    sll_int32                  :: i
    sll_int32                  :: sz
    sz = sll_get_num_nodes( layout )
    print *, 'limits: '
    do i=0,sz-1
       call view_box_3D(get_layout_3D_box( layout, i ))
    end do
    call flush(6)
  end subroutine sll_view_lims_3D

  subroutine sll_view_lims_4D( layout )
    type(layout_4D), pointer :: layout
    sll_int32                  :: i
    sll_int32                  :: sz
    sz = sll_get_num_nodes( layout )
    print *, 'limits: '
    do i=0,sz-1
       call view_box_4D(get_layout_4D_box( layout, i ))
    end do
    call flush(6)
  end subroutine sll_view_lims_4D

  subroutine sll_view_lims_6D( layout )
    type(layout_6D), pointer :: layout
    sll_int32                  :: i
    sll_int32                  :: sz
    sz = sll_get_num_nodes( layout )
    print *, 'limits: '
    do i=0,sz-1
       call view_box_6D(get_layout_6D_box( layout, i ))
    end do
    call flush(6)
  end subroutine sll_view_lims_6D


  ! the return value of intersect_boxes() is 'logical', and answers
  ! the question whether the boxes intersect or not. 'ans' is a box with
  ! the actual intersection between the argument boxes. In case that there
  ! is no intersection between the boxes the value [0,0] is returned. 

  function intersect_boxes_2D( b1, b2, ans )
    intrinsic                 :: min, max
    logical                   :: intersect_boxes_2D
    type(box_2D), intent(in)  :: b1, b2
    type(box_2D), intent(out) :: ans
    sll_int32                 :: loi, hii
    sll_int32                 :: loj, hij
    sll_int32                 :: loib1, hiib1
    sll_int32                 :: lojb1, hijb1
    sll_int32                 :: loib2, hiib2
    sll_int32                 :: lojb2, hijb2
    ! FIXME: add error checking, if boxes are null, for instance.
    loib1 = get_box_2D_i_min(b1)
    hiib1 = get_box_2D_i_max(b1)
    lojb1 = get_box_2D_j_min(b1)
    hijb1 = get_box_2D_j_max(b1)

    loib2 = get_box_2D_i_min(b2)
    hiib2 = get_box_2D_i_max(b2)
    lojb2 = get_box_2D_j_min(b2)
    hijb2 = get_box_2D_j_max(b2)

    SLL_ASSERT( (loib1 .le. hiib1) .and. (loib2 .le. hiib2) )
    SLL_ASSERT( (lojb1 .le. hijb1) .and. (lojb2 .le. hijb2) )

    loi = max(loib1, loib2)
    hii = min(hiib1, hiib2)
    loj = max(lojb1, lojb2)
    hij = min(hijb1, hijb2)

    if( (loi .gt. hii) .or. (loj .gt. hij) ) then 
       ans%i_min = 0
       ans%i_max = 0
       ans%j_min = 0
       ans%j_max = 0
       intersect_boxes_2D = .false.
    else
       ans%i_min = loi
       ans%i_max = hii
       ans%j_min = loj
       ans%j_max = hij
       intersect_boxes_2D = .true.
    end if
  end function intersect_boxes_2D

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

  function intersect_boxes_4D( b1, b2, ans )
    intrinsic                 :: min, max
    logical                   :: intersect_boxes_4D
    type(box_4D), intent(in)  :: b1, b2
    type(box_4D), intent(out) :: ans
    sll_int32                 :: loi, hii
    sll_int32                 :: loj, hij
    sll_int32                 :: lok, hik
    sll_int32                 :: lol, hil
    sll_int32                 :: loib1, hiib1
    sll_int32                 :: lojb1, hijb1
    sll_int32                 :: lokb1, hikb1
    sll_int32                 :: lolb1, hilb1
    sll_int32                 :: loib2, hiib2
    sll_int32                 :: lojb2, hijb2
    sll_int32                 :: lokb2, hikb2
    sll_int32                 :: lolb2, hilb2

    ! FIXME: add error checking, if boxes are null, for instance.
    loib1 = get_box_4D_i_min(b1)
    hiib1 = get_box_4D_i_max(b1)
    lojb1 = get_box_4D_j_min(b1)
    hijb1 = get_box_4D_j_max(b1)
    lokb1 = get_box_4D_k_min(b1)
    hikb1 = get_box_4D_k_max(b1)
    lolb1 = get_box_4D_l_min(b1)
    hilb1 = get_box_4D_l_max(b1)

    loib2 = get_box_4D_i_min(b2)
    hiib2 = get_box_4D_i_max(b2)
    lojb2 = get_box_4D_j_min(b2)
    hijb2 = get_box_4D_j_max(b2)
    lokb2 = get_box_4D_k_min(b2)
    hikb2 = get_box_4D_k_max(b2)
    lolb2 = get_box_4D_l_min(b2)
    hilb2 = get_box_4D_l_max(b2)

    SLL_ASSERT( (loib1 .le. hiib1) .and. (loib2 .le. hiib2) )
    SLL_ASSERT( (lojb1 .le. hijb1) .and. (lojb2 .le. hijb2) )
    SLL_ASSERT( (lokb1 .le. hikb1) .and. (lokb2 .le. hikb2) )
    SLL_ASSERT( (lolb1 .le. hilb1) .and. (lolb2 .le. hilb2) )

    loi = max(loib1, loib2)
    hii = min(hiib1, hiib2)
    loj = max(lojb1, lojb2)
    hij = min(hijb1, hijb2)
    lok = max(lokb1, lokb2)
    hik = min(hikb1, hikb2)
    lol = max(lolb1, lolb2)
    hil = min(hilb1, hilb2)

    if( (loi .gt. hii) .or. (loj .gt. hij) .or. (lok .gt. hik) .or. &
        (lol .gt. hil) ) then 
       ans%i_min = 0
       ans%i_max = 0
       ans%j_min = 0
       ans%j_max = 0
       ans%k_min = 0
       ans%k_max = 0
       ans%l_min = 0
       ans%l_max = 0
       intersect_boxes_4D = .false.
    else
       ans%i_min = loi
       ans%i_max = hii
       ans%j_min = loj
       ans%j_max = hij
       ans%k_min = lok
       ans%k_max = hik
       ans%l_min = lol
       ans%l_max = hil
       intersect_boxes_4D = .true.
    end if
  end function intersect_boxes_4D

  function intersect_boxes_6D( b1, b2, ans )
    intrinsic                 :: min, max
    logical                   :: intersect_boxes_6D
    type(box_6D), intent(in)  :: b1, b2
    type(box_6D), intent(out) :: ans
    sll_int32                 :: loi, hii
    sll_int32                 :: loj, hij
    sll_int32                 :: lok, hik
    sll_int32                 :: lol, hil
    sll_int32                 :: lom, him
    sll_int32                 :: lon, hin
    sll_int32                 :: loib1, hiib1
    sll_int32                 :: lojb1, hijb1
    sll_int32                 :: lokb1, hikb1
    sll_int32                 :: lolb1, hilb1
    sll_int32                 :: lomb1, himb1
    sll_int32                 :: lonb1, hinb1

    sll_int32                 :: loib2, hiib2
    sll_int32                 :: lojb2, hijb2
    sll_int32                 :: lokb2, hikb2
    sll_int32                 :: lolb2, hilb2
    sll_int32                 :: lomb2, himb2
    sll_int32                 :: lonb2, hinb2

    ! FIXME: add error checking, if boxes are null, for instance.
    loib1 = get_box_6D_i_min(b1)
    hiib1 = get_box_6D_i_max(b1)
    lojb1 = get_box_6D_j_min(b1)
    hijb1 = get_box_6D_j_max(b1)
    lokb1 = get_box_6D_k_min(b1)
    hikb1 = get_box_6D_k_max(b1)
    lolb1 = get_box_6D_l_min(b1)
    hilb1 = get_box_6D_l_max(b1)
    lomb1 = get_box_6D_m_min(b1)
    himb1 = get_box_6D_m_max(b1)
    lonb1 = get_box_6D_n_min(b1)
    hinb1 = get_box_6D_n_max(b1)

    loib2 = get_box_6D_i_min(b2)
    hiib2 = get_box_6D_i_max(b2)
    lojb2 = get_box_6D_j_min(b2)
    hijb2 = get_box_6D_j_max(b2)
    lokb2 = get_box_6D_k_min(b2)
    hikb2 = get_box_6D_k_max(b2)
    lolb2 = get_box_6D_l_min(b2)
    hilb2 = get_box_6D_l_max(b2)
    lomb2 = get_box_6D_m_min(b2)
    himb2 = get_box_6D_m_max(b2)
    lonb2 = get_box_6D_n_min(b2)
    hinb2 = get_box_6D_n_max(b2)

    SLL_ASSERT( (loib1 .le. hiib1) .and. (loib2 .le. hiib2) )
    SLL_ASSERT( (lojb1 .le. hijb1) .and. (lojb2 .le. hijb2) )
    SLL_ASSERT( (lokb1 .le. hikb1) .and. (lokb2 .le. hikb2) )
    SLL_ASSERT( (lolb1 .le. hilb1) .and. (lolb2 .le. hilb2) )
    SLL_ASSERT( (lomb1 .le. himb1) .and. (lomb2 .le. himb2) )
    SLL_ASSERT( (lonb1 .le. hinb1) .and. (lonb2 .le. hinb2) )

    loi = max(loib1, loib2)
    hii = min(hiib1, hiib2)
    loj = max(lojb1, lojb2)
    hij = min(hijb1, hijb2)
    lok = max(lokb1, lokb2)
    hik = min(hikb1, hikb2)
    lol = max(lolb1, lolb2)
    hil = min(hilb1, hilb2)
    lom = max(lomb1, lomb2)
    him = min(himb1, himb2)
    lon = max(lonb1, lonb2)
    hin = min(hinb1, hinb2)

    if( (loi .gt. hii) .or. (loj .gt. hij) .or. (lok .gt. hik) .or. &
        (lol .gt. hil) .or. (lom .gt. him) .or. (lon .gt. hin) ) then 
       ans%i_min = 0
       ans%i_max = 0
       ans%j_min = 0
       ans%j_max = 0
       ans%k_min = 0
       ans%k_max = 0
       ans%l_min = 0
       ans%l_max = 0
       ans%m_min = 0
       ans%m_max = 0
       ans%n_min = 0
       ans%n_max = 0
       intersect_boxes_6D = .false.
    else
       ans%i_min = loi
       ans%i_max = hii
       ans%j_min = loj
       ans%j_max = hij
       ans%k_min = lok
       ans%k_max = hik
       ans%l_min = lol
       ans%l_max = hil
       ans%m_min = lom
       ans%m_max = him
       ans%n_min = lon
       ans%n_max = hin
       intersect_boxes_6D = .true.
    end if
  end function intersect_boxes_6D


  function count_elements_in_box_2D( box )
    sll_int32                 :: count_elements_in_box_2D
    type(box_2D), intent (in) :: box
    sll_int32                 :: irange
    sll_int32                 :: jrange
    irange = get_box_2D_i_max(box) - get_box_2D_i_min(box) + 1
    jrange = get_box_2D_j_max(box) - get_box_2D_j_min(box) + 1
    count_elements_in_box_2D = irange*jrange
  end function count_elements_in_box_2D

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

  function count_elements_in_box_4D( box )
    sll_int32                 :: count_elements_in_box_4D
    type(box_4D), intent (in) :: box
    sll_int32                 :: irange
    sll_int32                 :: jrange
    sll_int32                 :: krange
    sll_int32                 :: lrange
    irange = get_box_4D_i_max(box) - get_box_4D_i_min(box) + 1
    jrange = get_box_4D_j_max(box) - get_box_4D_j_min(box) + 1
    krange = get_box_4D_k_max(box) - get_box_4D_k_min(box) + 1
    lrange = get_box_4D_l_max(box) - get_box_4D_l_min(box) + 1
    count_elements_in_box_4D = irange*jrange*krange*lrange
  end function count_elements_in_box_4D

  function count_elements_in_box_6D( box )
    sll_int32                 :: count_elements_in_box_6D
    type(box_6D), intent (in) :: box
    sll_int32                 :: irange
    sll_int32                 :: jrange
    sll_int32                 :: krange
    sll_int32                 :: lrange
    sll_int32                 :: mrange
    sll_int32                 :: nrange
    irange = get_box_6D_i_max(box) - get_box_6D_i_min(box) + 1
    jrange = get_box_6D_j_max(box) - get_box_6D_j_min(box) + 1
    krange = get_box_6D_k_max(box) - get_box_6D_k_min(box) + 1
    lrange = get_box_6D_l_max(box) - get_box_6D_l_min(box) + 1
    mrange = get_box_6D_m_max(box) - get_box_6D_m_min(box) + 1
    nrange = get_box_6D_n_max(box) - get_box_6D_n_min(box) + 1
    count_elements_in_box_6D = irange*jrange*krange*lrange*mrange*nrange
  end function count_elements_in_box_6D


  subroutine compute_local_sizes_6d( &
    layout, &
    loc_sz_i, &
    loc_sz_j, &
    loc_sz_k, &
    loc_sz_l, &
    loc_sz_m, &
    loc_sz_n )

    type(layout_6D), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32, intent(out) :: loc_sz_k
    sll_int32, intent(out) :: loc_sz_l
    sll_int32, intent(out) :: loc_sz_m
    sll_int32, intent(out) :: loc_sz_n
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: l_min
    sll_int32 :: l_max
    sll_int32 :: m_min
    sll_int32 :: m_max
    sll_int32 :: n_min
    sll_int32 :: n_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to compute_local_sizes_6d'
       print *, 'Exiting...'
       STOP
    end if
    my_rank = sll_get_collective_rank(get_layout_6D_collective(layout))
    i_min = get_layout_6D_i_min( layout, my_rank )
    i_max = get_layout_6D_i_max( layout, my_rank )
    j_min = get_layout_6D_j_min( layout, my_rank )
    j_max = get_layout_6D_j_max( layout, my_rank )
    k_min = get_layout_6D_k_min( layout, my_rank )
    k_max = get_layout_6D_k_max( layout, my_rank )
    l_min = get_layout_6D_l_min( layout, my_rank )
    l_max = get_layout_6D_l_max( layout, my_rank )
    m_min = get_layout_6D_m_min( layout, my_rank )
    m_max = get_layout_6D_m_max( layout, my_rank )
    n_min = get_layout_6D_n_min( layout, my_rank )
    n_max = get_layout_6D_n_max( layout, my_rank )
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
    loc_sz_k = k_max - k_min + 1
    loc_sz_l = l_max - l_min + 1
    loc_sz_m = m_max - m_min + 1
    loc_sz_n = n_max - n_min + 1
  end subroutine compute_local_sizes_6d


  subroutine compute_local_sizes_4d( &
    layout, &
    loc_sz_i, &
    loc_sz_j, &
    loc_sz_k, &
    loc_sz_l )

    type(layout_4D), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32, intent(out) :: loc_sz_k
    sll_int32, intent(out) :: loc_sz_l
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: l_min
    sll_int32 :: l_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to compute_local_sizes_4d'
       print *, 'Exiting...'
       STOP
    end if
    my_rank = sll_get_collective_rank(get_layout_4D_collective(layout))
    i_min = get_layout_4D_i_min( layout, my_rank )
    i_max = get_layout_4D_i_max( layout, my_rank )
    j_min = get_layout_4D_j_min( layout, my_rank )
    j_max = get_layout_4D_j_max( layout, my_rank )
    k_min = get_layout_4D_k_min( layout, my_rank )
    k_max = get_layout_4D_k_max( layout, my_rank )
    l_min = get_layout_4D_l_min( layout, my_rank )
    l_max = get_layout_4D_l_max( layout, my_rank )
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
    loc_sz_k = k_max - k_min + 1
    loc_sz_l = l_max - l_min + 1
  end subroutine compute_local_sizes_4d

  subroutine compute_local_sizes_3d( &
    layout, &
    loc_sz_i, &
    loc_sz_j, &
    loc_sz_k )

    type(layout_3D), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32, intent(out) :: loc_sz_k
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to compute_local_sizes_3d'
       print *, 'Exiting...'
       STOP
    end if
    my_rank = sll_get_collective_rank(get_layout_3D_collective(layout))
    i_min = get_layout_3D_i_min( layout, my_rank )
    i_max = get_layout_3D_i_max( layout, my_rank )
    j_min = get_layout_3D_j_min( layout, my_rank )
    j_max = get_layout_3D_j_max( layout, my_rank )
    k_min = get_layout_3D_k_min( layout, my_rank )
    k_max = get_layout_3D_k_max( layout, my_rank )
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
    loc_sz_k = k_max - k_min + 1
  end subroutine compute_local_sizes_3d

  subroutine compute_local_sizes_2d( &
    layout, &
    loc_sz_i, &
    loc_sz_j )

    type(layout_2D), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to compute_local_sizes_2d'
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

  !***************************************************************************
  !
  ! Functions that operate with more than one dimension.
  !
  ! We've found that sometimes it may be necessary to interpret a layout of
  ! a given dimension, say 4D, as a layout of a different dimension, like 2D.
  ! For this operation to be valid some conditions need to be met. For
  ! example, in the 4D to 2D case, some array dimensions should be equal to 1
  ! in order to reinterpret the layout properly. Here we just proceed to add
  ! these routines as we need them.
  !
  !***************************************************************************

  ! layout_2D_from_layout_4D() takes a 4D layout that describes the distribution
  ! of a 4D array of dimensions npx1 X npx2 X 1 X 1 and returns a 2D layout,
  ! defined over the same collective, which describes the distribution of a 2D
  ! array of dimensions npx1 X npx2. Note that it assumes that it is the last
  ! two dimensions which are of size 1.
  ! 
  ! This function is special in that it allocates the new layout to be returned.
  ! So the usual interface of declaring the layout, calling new_layout() and
  ! then initializing is not followed. This irregularity is itself a bit of
  ! a problem, but may be a sign that the usual way to allocate and initialize
  ! layouts might need to be merged.
  function layout_2D_from_layout_4D( layout4d )
    type(layout_2D), pointer :: layout_2D_from_layout_4D
    type(layout_4D), pointer :: layout4d
    type(sll_collective_t), pointer :: coll
    sll_int32                :: coll_size
    sll_int32                :: process
    sll_int32                :: i_min
    sll_int32                :: i_max
    sll_int32                :: j_min
    sll_int32                :: j_max
    sll_int32                :: k_min
    sll_int32                :: k_max
    sll_int32                :: l_min
    sll_int32                :: l_max

    SLL_ASSERT( associated(layout4d) )
    coll                     => get_layout_collective( layout4d )
    coll_size                = sll_get_collective_size( coll )
    layout_2D_from_layout_4D => new_layout_2d( coll )
    ! Just copy the contents of the layout
    do process=0, coll_size-1
       i_min = get_layout_i_min( layout4d, process )
       i_max = get_layout_i_max( layout4d, process )
       j_min = get_layout_j_min( layout4d, process )
       j_max = get_layout_j_max( layout4d, process )
       call set_layout_i_min( layout_2D_from_layout_4D, process, i_min )
       call set_layout_i_max( layout_2D_from_layout_4D, process, i_max )
       call set_layout_j_min( layout_2D_from_layout_4D, process, j_min )
       call set_layout_j_max( layout_2D_from_layout_4D, process, j_max )
       ! For safety, check if there is any loss of information
       k_min = get_layout_k_min( layout4d, process )
       k_max = get_layout_k_max( layout4d, process )
       l_min = get_layout_l_min( layout4d, process )
       l_max = get_layout_l_max( layout4d, process )
       if( (k_min .ne. 1) .or. (k_max .ne. 1) .or. &
           (l_min .ne. 1) .or. (l_max .ne. 1) ) then
           print *, 'WARNING, layout_2D_from_layout_4D(): there is loss of ',&
                'information in the convertion. Printing values:'
           print *, 'k_min = ', k_min
           print *, 'k_max = ', k_max
           print *, 'l_min = ', l_min
           print *, 'l_max = ', l_max
        end if
    end do
  end function layout_2D_from_layout_4D

#if 0
  function layout_4D_from_layout_2D( layout2d )
    type(layout_4D), pointer :: layout_4D_from_layout_2D
    type(layout_2D), pointer :: layout2d
    type(sll_collective_t), pointer :: coll
    sll_int32                :: coll_size
    sll_int32                :: process
    sll_int32                :: i_min
    sll_int32                :: i_max
    sll_int32                :: j_min
    sll_int32                :: j_max
    sll_int32                :: k_min
    sll_int32                :: k_max
    sll_int32                :: l_min
    sll_int32                :: l_max

    SLL_ASSERT( associated(layout4d) )
    coll                     => get_layout_collective( layout4d )
    coll_size                = sll_get_collective_size( coll )
    layout_2D_from_layout_4D => new_layout_4d( coll )
    ! Just copy the contents of the layout
    do process=0, coll_size-1
       i_min = get_layout_i_min( layout4d, process )
       i_max = get_layout_i_max( layout4d, process )
       j_min = get_layout_j_min( layout4d, process )
       j_max = get_layout_j_max( layout4d, process )
       call set_layout_i_min( layout_2D_from_layout_4D, process, i_min )
       call set_layout_i_max( layout_2D_from_layout_4D, process, i_max )
       call set_layout_j_min( layout_2D_from_layout_4D, process, j_min )
       call set_layout_j_max( layout_2D_from_layout_4D, process, j_max )
       ! For safety, check if there is any loss of information
       k_min = get_layout_k_min( layout4d, process )
       k_max = get_layout_k_max( layout4d, process )
       l_min = get_layout_l_min( layout4d, process )
       l_max = get_layout_l_max( layout4d, process )
       if( (k_min .ne. 1) .or. (k_max .ne. 1) .or.
           (l_min .ne. 1) .or. (l_max .ne. 1) ) then
           print *, 'WARNING, layout_2D_from_layout_4D(): there is loss of ',&
                'information in the convertion. Printing values:'
           print *, 'k_min = ', k_min
           print *, 'k_max = ', k_max
           print *, 'l_min = ', l_min
           print *, 'l_max = ', l_max
        end if
    end do
  end function layout_2D_from_layout_4D
#endif

end module sll_remapper
