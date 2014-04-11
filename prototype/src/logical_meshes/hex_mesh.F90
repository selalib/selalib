!**************************************************************
!  Copyright INRIA
!  Authors : 
!     Laura Mendoza
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


module hex_mesh
#include "sll_working_precision.h"
#include "sll_memory.h"
  implicit none

  type hex_mesh_2d
     sll_int32  :: num_cells  ! number of cells in any direction parting from origin
     sll_int32  :: num_pts_tot! number of total points
     sll_real64 :: radius     ! distance between origin and external vertex
     sll_real64 :: center_x1  ! x1 cartesian coordinate of the origin
     sll_real64 :: center_x2  ! x2 cartesian coordinate of the origin
     sll_real64 :: delta      ! cell spacing
     ! generator vectors (r1, r2, r3) coordinates -- need to be scaled by delta
     sll_real64 :: r1_x1
     sll_real64 :: r1_x2
     sll_real64 :: r2_x1
     sll_real64 :: r2_x2
     sll_real64 :: r3_x1
     sll_real64 :: r3_x2
  end type hex_mesh_2d

  type hex_mesh_2d_ptr
     type(hex_mesh_2d), pointer :: hm
  end type hex_mesh_2d_ptr

  ! this should be sll_delete library-wide...
  interface delete
     module procedure delete_hex_mesh_2d
  end interface delete

  interface sll_display
      module procedure display_hex_mesh_2d
  end interface sll_display

contains

#define TEST_PRESENCE_AND_ASSIGN_VAL( obj, arg, slot, default_val ) \
  if( present(arg) ) then ; \
    obj%slot = arg; \
  else; \
    obj%slot = default_val; \
end if

  function new_hex_mesh_2d( &
    num_cells, &
    center_x1, &
    center_x2, &
    r11, &
    r12, &
    r21, &
    r22, &
    r31, &
    r32, &
    radius ) result(m)

    type(hex_mesh_2d), pointer :: m
    sll_int32, intent(in)  :: num_cells
    sll_real64, optional, intent(in) :: radius
    sll_real64, optional, intent(in) :: center_x1
    sll_real64, optional, intent(in) :: center_x2
    sll_real64, optional, intent(in) :: r11, r12
    sll_real64, optional, intent(in) :: r21, r22
    sll_real64, optional, intent(in) :: r31, r32
    sll_int32 :: ierr

    SLL_ALLOCATE(m, ierr)
    call initialize_hex_mesh_2d( &
         m, &
         num_cells, &
         radius, &
         center_x1, &
         center_x2, &
         r11, &
         r12, &
         r21, &
         r22, &
         r31, &
         r32)

  end function new_hex_mesh_2d


  subroutine initialize_hex_mesh_2d( &
    m, & 
    num_cells, &
    radius,    &
    center_x1, &
    center_x2, &
    r1_x1,     &
    r1_x2,     &
    r2_x1,     &
    r2_x2,     &
    r3_x1,     &
    r3_x2)


    type(hex_mesh_2d), pointer :: m
    sll_int32, intent(in)  :: num_cells
    sll_real64, optional, intent(in) :: radius
    sll_real64, optional, intent(in) :: center_x1
    sll_real64, optional, intent(in) :: center_x2
    sll_real64, optional, intent(in) :: r1_x1, r1_x2
    sll_real64, optional, intent(in) :: r2_x1, r2_x2
    sll_real64, optional, intent(in) :: r3_x1, r3_x2

    TEST_PRESENCE_AND_ASSIGN_VAL( m, center_x1, center_x1, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, center_x2, center_x2, 0.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, radius, radius, 1.0_f64 )
    TEST_PRESENCE_AND_ASSIGN_VAL( m, r1_x1, r1_x1, sqrt(real(3, f64))*0.5) 
    TEST_PRESENCE_AND_ASSIGN_VAL( m, r1_x2, r1_x2, 0.5)
    TEST_PRESENCE_AND_ASSIGN_VAL( m, r2_x1, r2_x1, -sqrt(real(3, f64))*0.5)
    TEST_PRESENCE_AND_ASSIGN_VAL( m, r2_x2, r2_x2, 0.5)
    TEST_PRESENCE_AND_ASSIGN_VAL( m, r3_x1, r3_x1, 0.0)
    TEST_PRESENCE_AND_ASSIGN_VAL( m, r3_x2, r3_x2, 1.0)

    m%num_cells = num_cells
    m%delta = m%radius/real(num_cells,f64)
    m%num_pts_tot = 6*m%num_cells*(m%num_cells+1)/2+1

    ! resizing :
    m%r1_x1 = m%r1_x1 * m%delta
    m%r1_x2 = m%r1_x2 * m%delta
    m%r2_x1 = m%r2_x1 * m%delta
    m%r2_x2 = m%r2_x2 * m%delta
    m%r3_x1 = m%r3_x1 * m%delta
    m%r3_x2 = m%r3_x2 * m%delta

    if ( m%radius <= 0.) then
       print*,'ERROR, initialize_hex_mesh_2d(): ', &
              'Problem to construct the mesh 2d '
       print*,'because radius <= 0.'
       STOP
    end if
    if ( m%num_cells <= 0) then
       print*,'ERROR, initialize_hex_mesh_2d(): ', &
            'Problem to construct the mesh 2d '
       print*,'because num_cells <= 0.'
       STOP
    end if
  end subroutine initialize_hex_mesh_2d


  function x1_node(mesh, k1, k2) result(val)
    type(hex_mesh_2d), pointer :: mesh
    sll_int32, intent(in)  :: k1
    sll_int32, intent(in)  :: k2
    sll_real64 :: val
    val = mesh%r1_x1*k1 + mesh%r2_x1*k2 + mesh%center_x1
  end function x1_node


  function x2_node(mesh, k1, k2) result(val)
    type(hex_mesh_2d), pointer :: mesh
    sll_int32, intent(in)  :: k1
    sll_int32, intent(in)  :: k2
    sll_real64  :: val
    val = mesh%r1_x2*k1 + mesh%r2_x2*k2 + mesh%center_x1
  end function x2_node


  function global_index(k1, k2) result(val)
    ! Takes the coordinates (k1,k2) on the (r1,r2) basis and 
    ! returns global index of that mesh point.
    sll_int32, intent(in)   :: k1
    sll_int32, intent(in)   :: k2
    sll_int32 :: val
    sll_int32 :: hex_num, first_index

    ! We compute the hexagon-ring number
    if (k1*k2 .ge. 0.) then
       hex_num = max( abs(k1), abs(k2))
    else 
       hex_num = abs(k1) + abs(k2)
    end if

    ! We get the index of the first point in this ring
    first_index = 1 + 3*(hex_num - 1)*hex_num


    if ( (k1 .eq. 0) .and. (k2 .eq. 0)) then
       ! Origin :
       val = 0
    elseif (k1 .eq. hex_num) then
       ! First edge of hexagone
       val = first_index + k2
    elseif (k2 .eq. hex_num) then
       ! Second edge of hexagone
       val = first_index + 2*hex_num - k1
    elseif (( k1 .lt. 0) .and. (k2 .gt. 0)) then
       ! Third edge of hexagone
       val = first_index + 3*hex_num - k2
    elseif (k1 .eq. -hex_num) then
       ! Forth edge of hexagone
       val = first_index + 3*hex_num + abs(k2)
    elseif (k2 .eq. -hex_num) then
       ! Fifth edge of hexagone
       val = first_index + 5*hex_num - abs(k1)
    elseif (( k1 .gt. 0) .and. (k2 .lt. 0)) then
       ! Sixth edge of hexagone
       val = first_index + 6*hex_num - abs(k2)
    else 
       print *, "ERROR : in global_index(k1,k2)"
       print *, "Not recognized combination of k1,k2 ", k1, k2
       STOP 'global_index'
    endif

  end function global_index


  function get_hex_num(index) result(hex_num)
      ! returns the number of nested hexagon on which the point is
      ! = returns the number of cells from center to point
      sll_int32 :: index
      sll_int32 :: hex_num
      sll_int32 :: flag

      
      hex_num = 0
      flag = 0

      if (index .lt. 0) then
         print *, "ERROR : in get_hex_num(index)"
         print *, "Global index cannot be negative"
         STOP 'get_hex_num'         
      else if (index .eq. 0) then
         hex_num = 0
      else
         hex_num = 1
         do while(flag .eq. 0)
             if (index .gt. 3*hex_num*(hex_num+1)) then
                hex_num = hex_num + 1
             else
                flag = 1
             endif
         enddo
      endif
  end function get_hex_num


  function from_global_index_k1(index) result(k1)
      sll_int32 :: index
      sll_int32 :: hex_num
      sll_int32 :: first_index
      sll_int32 :: last_index
      sll_int32 :: k1

      hex_num = get_hex_num(index)
      first_index = 3*(hex_num - 1)*hex_num + 1
      last_index  = 3*(hex_num + 1)*hex_num

      if (hex_num .eq. 0) then
         k1 = 0
      elseif (index .le. first_index+hex_num) then
         !index on first edge
         k1 = hex_num
      elseif (index .le. first_index+2*hex_num) then
         !index on second edge
         k1 = first_index + 2*hex_num - index
      elseif (index .le. first_index+3*hex_num-1) then
         !index on third edge
         k1 = first_index + 2*hex_num - index 
      elseif (index .le. first_index+4*hex_num) then
         !index on forth edge
         k1 = - hex_num
      elseif (index .le.  first_index+5*hex_num) then
         !index on fifth edge
         k1 = index - first_index - 5*hex_num
      elseif (index .le. last_index) then
         !index of sixth edge
         k1 = hex_num -last_index + index - 1
      else
       print *, "ERROR : in from_global_index_k1(index)"
       print *, "Not recognized index"
       STOP 'from_global_index'
    endif
  end function from_global_index_k1


  function from_global_index_k2(index) result(k2)

      sll_int32 :: index
      sll_int32 :: hex_num
      sll_int32 :: first_index
      sll_int32 :: last_index
      sll_int32 :: k2

      hex_num = get_hex_num(index)
      first_index = 3*(hex_num - 1)*hex_num + 1
      last_index  = 3*(hex_num + 1)*hex_num
      
      if (hex_num .eq. 0) then
         k2 = 0
      elseif (index .le. first_index+hex_num) then
         !index on first edge
         k2 = index - first_index
      elseif (index .le. first_index+2*hex_num) then
         !index on second edge
         k2 = hex_num
      elseif (index .le. first_index+3*hex_num-1) then
         !index on third edge
         k2 = first_index + 3*hex_num - index 
      elseif (index .le. first_index+4*hex_num) then
         !index on forth edge
         k2 = first_index + 3*hex_num - index
      elseif (index .le.  first_index+5*hex_num) then
         !index on fifth edge
         k2 = -hex_num
      elseif (index .le. last_index) then
         !index of sixth edge
         k2 = -1 + index - last_index
      else
       print *, "ERROR : in from_global_index_k2(index)"
       print *, "Not recognized index"
       STOP 'from_global_index'
    endif
  end function from_global_index_k2

  function from_global_x1(mesh, index) result(x1)
      type(hex_mesh_2d), pointer :: mesh
      sll_int32  :: index
      sll_int32  :: k1
      sll_int32  :: k2
      sll_real64 :: x1

      k1 = from_global_index_k1(index)
      k2 = from_global_index_k2(index)
      
      x1 = x1_node(mesh, k1, k2)
  end function from_global_x1

  function from_global_x2(mesh, index) result(x2)
      type(hex_mesh_2d), pointer :: mesh
      sll_int32  :: index
      sll_int32  :: k1
      sll_int32  :: k2
      sll_real64 :: x2

      k1 = from_global_index_k1(index)
      k2 = from_global_index_k2(index)
      
      x2 = x2_node(mesh, k1, k2)
  end function from_global_x2


  function from_cart_index_k1(mesh, x1, x2) result(k1)
      type(hex_mesh_2d), pointer :: mesh
      sll_real64 :: x1
      sll_real64 :: x2
      sll_int32  :: k1
      sll_real64 :: delta
      
      delta = mesh%r1_x1 * mesh%r2_x2 - mesh%r2_x1 * mesh%r1_x2
      k1 = nint((mesh%r2_x2 * x1 - mesh%r2_x1 * x2)/delta)
  end function from_cart_index_k1

  function from_cart_index_k2(mesh, x1, x2) result(k2)
      type(hex_mesh_2d), pointer :: mesh
      sll_real64 :: x1
      sll_real64 :: x2
      sll_int32  :: k2
      sll_real64 :: delta
      
      delta = mesh%r1_x1 * mesh%r2_x2 - mesh%r2_x1 * mesh%r1_x2
      k2 = nint((mesh%r1_x1 * x2 - mesh%r1_x2 * x1)/delta)
  end function from_cart_index_k2

  function local_index(ref_index,j) result(new_index)
      ! returns the index of the point Pj as if the 
      ! notation had as origin Pi
      ! ie. local_index(i,i) = 0
      sll_int32 :: ref_index, j
      sll_int32 :: k1_i, k2_i
      sll_int32 :: k1_j, k2_j
      sll_int32 :: new_index

      k1_i = from_global_index_k1(ref_index)
      k2_i = from_global_index_k2(ref_index)
      k1_j = from_global_index_k1(j)
      k2_j = from_global_index_k2(j)

      new_index = global_index(k1_i - k1_j, k2_i - k2_j)

  end function local_index


  function find_neighbour(ref_index, local_index) result(global)
      ! returns the global index of the point which has as
      ! local index local_index in the ref_index system
      ! ie. find_neighbour(0, i) = i
      sll_int32 :: ref_index, local_index
      sll_int32 :: k1_i, k2_i
      sll_int32 :: k1_j, k2_j
      sll_int32 :: global
      
      k1_i = from_global_index_k1(ref_index)
      k2_i = from_global_index_k2(ref_index)
      k1_j = from_global_index_k1(local_index)
      k2_j = from_global_index_k2(local_index)

      global = global_index(k1_i + k1_j, k2_i + k2_j)

  end function find_neighbour

      
  subroutine display_hex_mesh_2d(mesh)
    type(hex_mesh_2d), pointer :: mesh

    write(*,"(/,(a))") '2D mesh : num_cells   num_pts        center_x1       center_x2 &
     &       radius'
    write(*,"(10x,2(i4,10x),3(g13.3,1x))") mesh%num_cells,  &
                                         mesh%num_pts_tot,&
                                         mesh%center_x1,  &
                                         mesh%center_x2,  &
                                         mesh%radius
  end subroutine display_hex_mesh_2d


  subroutine write_hex_mesh_2d(mesh, name)
    type(hex_mesh_2d), pointer :: mesh
    character(len=*) :: name
    sll_int32  :: i
    sll_int32  :: hex_num
    sll_int32  :: num_pts_tot
    sll_int32  :: k1, k2
    sll_int32, parameter :: out_unit=20

    open (unit=out_unit,file=name,action="write",status="replace")

    num_pts_tot = mesh%num_pts_tot

!    write(*,"(/,(a))") 'hex mesh : num_pnt    x1     x2'

    do i=0, num_pts_tot-1
       hex_num = get_hex_num(i)
       k1 = from_global_index_k1(i)
       k2 = from_global_index_k2(i)
       write (out_unit, "(4(i2,1x),2(g13.3,1x))") i,                &
                                           hex_num,                 &
                                           k1,                      &
                                           k2,                      &
                                           from_global_x1(mesh, i), &
                                           from_global_x2(mesh, i)
    end do

    close(out_unit)
  end subroutine write_hex_mesh_2d

  subroutine write_field_hex_mesh(mesh, field, name)
    type(hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:) :: field
    character(len=*) :: name
    sll_int32  :: i
    sll_int32  :: num_pts_tot
    sll_real64 :: x1, x2
    sll_int32, parameter :: out_unit=20

    open (unit=out_unit,file=name,action="write",status="replace")

    num_pts_tot = mesh%num_pts_tot

    do i=0, num_pts_tot-1
       x1 = from_global_x1(mesh, i)
       x2 = from_global_x2(mesh, i)
      write (out_unit, "(3(g13.3,1x))") x1, &
                                        x2, &
                                        field(i+1)
    end do

    close(out_unit)
  end subroutine write_field_hex_mesh


  subroutine delete_hex_mesh_2d( mesh )
    type(hex_mesh_2d), pointer :: mesh
    sll_int32 :: ierr
    if(.not. associated(mesh))then
       print *, 'delete_hex_mesh_2d'
       print *, 'ERROR: passed argument is not associated'
       print *, '       Crash imminent...'
       STOP
    end if
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_hex_mesh_2d


  
#undef TEST_PRESENCE_AND_ASSIGN_VAL

end module hex_mesh
