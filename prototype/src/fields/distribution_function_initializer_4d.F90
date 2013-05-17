module sll_test_4d_initializer
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_constants
  use sll_scalar_field_initializers_base
  use sll_collective
  use sll_remapper
  implicit none


  ! This is a simplistic initializer aimed at a 4d cartesian distribution
  ! function, periodic in x and y, and compact-ish in vx and vy.
  !
  ! Basically:
  !
  ! f(x,y,vx,vy) = pi/(8*vt^2)*sin(pi*x)*sin(pi*y)*exp(-(vx^2+vy^2)/(2*vt^2))
  !
  ! It is meant to be used in the intervals:
  ! x:  [ 0,1]
  ! y:  [ 0,1]
  ! vx: [-1,1]
  ! vy: [-1,1]
  ! thus vthermal has to be small relative to 1 to make the support as
  ! 'compact' as possible. 
  !
  ! Issue to resolve: to initialize the data we need the 'mesh'/'mapping'
  ! information. How to deal with this in 4D? Do we pass multiple 2D
  ! meshes? As something tentative, here we just write an ad hoc 4D
  ! cartesian mesh.

#ifdef STDF95
  type :: init_test_4d_par
     sll_int32 :: data_position
#else
  type, extends(scalar_field_4d_initializer_base) :: init_test_4d_par
#endif
     sll_real64 :: v_thermal
     type(layout_4D), pointer :: data_layout
     type(simple_cartesian_4d_mesh), pointer :: mesh_4d
#ifndef STDF95
   contains
     procedure, pass(init_obj) :: initialize => load_test_4d_initializer
     procedure, pass(init_obj) :: f_of_4args => compact_4d_field
#endif
  end type init_test_4d_par
  
  ! This has to end up somewhere else, if it is to stay in the library
  type simple_cartesian_4d_mesh
     sll_int32  :: num_cells1
     sll_int32  :: num_cells2
     sll_int32  :: num_cells3
     sll_int32  :: num_cells4
     sll_real64 :: x1_min
     sll_real64 :: x1_max
     sll_real64 :: x2_min
     sll_real64 :: x2_max
     sll_real64 :: x3_min
     sll_real64 :: x3_max
     sll_real64 :: x4_min
     sll_real64 :: x4_max
     sll_real64 :: delta_x1
     sll_real64 :: delta_x2
     sll_real64 :: delta_x3
     sll_real64 :: delta_x4
  end type simple_cartesian_4d_mesh

contains

  function new_cartesian_4d_mesh( &
    num_cells1, &
    num_cells2, &
    num_cells3, &
    num_cells4, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max, &
    x3_min, &
    x3_max, &
    x4_min, &
    x4_max )

    type(simple_cartesian_4d_mesh), pointer :: new_cartesian_4d_mesh
    sll_int32, intent(in)  :: num_cells1
    sll_int32, intent(in)  :: num_cells2
    sll_int32, intent(in)  :: num_cells3
    sll_int32, intent(in)  :: num_cells4
    sll_real64, intent(in) :: x1_min
    sll_real64, intent(in) :: x1_max
    sll_real64, intent(in) :: x2_min
    sll_real64, intent(in) :: x2_max
    sll_real64, intent(in) :: x3_min
    sll_real64, intent(in) :: x3_max
    sll_real64, intent(in) :: x4_min
    sll_real64, intent(in) :: x4_max
    sll_int32 :: ierr

    SLL_ALLOCATE(new_cartesian_4d_mesh, ierr)
    new_cartesian_4d_mesh%num_cells1 = num_cells1
    new_cartesian_4d_mesh%num_cells2 = num_cells2
    new_cartesian_4d_mesh%num_cells3 = num_cells3
    new_cartesian_4d_mesh%num_cells4 = num_cells4
    new_cartesian_4d_mesh%x1_min     = x1_min
    new_cartesian_4d_mesh%x1_max     = x1_max
    new_cartesian_4d_mesh%x2_min     = x2_min
    new_cartesian_4d_mesh%x2_max     = x2_max
    new_cartesian_4d_mesh%x3_min     = x3_min
    new_cartesian_4d_mesh%x3_max     = x3_max
    new_cartesian_4d_mesh%x4_min     = x4_min
    new_cartesian_4d_mesh%x4_max     = x4_max
    new_cartesian_4d_mesh%delta_x1   = (x1_max - x1_min)/real(num_cells1,f64)
    new_cartesian_4d_mesh%delta_x2   = (x2_max - x2_min)/real(num_cells2,f64)
    new_cartesian_4d_mesh%delta_x3   = (x3_max - x3_min)/real(num_cells3,f64)
    new_cartesian_4d_mesh%delta_x4   = (x4_max - x4_min)/real(num_cells4,f64)
  end function new_cartesian_4d_mesh

  subroutine delete_cartesian_4d_mesh( mesh )
    type(simple_cartesian_4d_mesh), pointer :: mesh
    sll_int32 :: ierr
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_cartesian_4d_mesh

#ifdef STDF95
  subroutine init_test_4d_par_initialize( &
    init_obj, &
    data_position, &
    mesh_4d, &
    v_thermal, &
    layout )

    type(init_test_4d_par), intent(inout)  :: init_obj
#else
  subroutine load_test_4d_initializer( &
    init_obj, &
    data_position, &
    mesh_4d, &
    v_thermal, &
    layout )

    class(init_test_4d_par), intent(inout)  :: init_obj
#endif
    sll_int32                               :: data_position
    sll_real64, intent(in)                  :: v_thermal
    type(layout_4D), pointer                :: layout
    type(simple_cartesian_4d_mesh), pointer :: mesh_4d

    if( .not. associated(layout) ) then
       print *, 'initialize_test_4d(): ERROR, passed layout is not initialized'
       stop
    end if

    if( .not. associated(mesh_4d) ) then
       print *, 'initialize_test_4d(): ERROR, passed 4d mesh is not initialized'
       stop
    end if

    init_obj%data_position =  data_position
    init_obj%mesh_4d       => mesh_4d
    init_obj%v_thermal     =  v_thermal
    init_obj%data_layout   => layout
  end subroutine 

#ifdef STDF95
  subroutine init_test_4d_par_f_of_4args( init_obj, data_out )
    type(init_test_4d_par), intent(inout)      :: init_obj
#else
  subroutine compact_4d_field( init_obj, data_out )
    class(init_test_4d_par), intent(inout)      :: init_obj
#endif
    sll_real64, dimension(:,:,:,:), intent(out) :: data_out
    type(layout_4D), pointer                    :: layout
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_int32  :: i_min
    sll_int32  :: i_max
    sll_int32  :: j_min
    sll_int32  :: j_max
    sll_int32  :: k_min
    sll_int32  :: k_max
    sll_int32  :: l_min
    sll_int32  :: l_max
    sll_int32  :: myrank
    sll_int32  :: num_pts1
    sll_int32  :: num_pts2
    sll_int32  :: num_pts3
    sll_int32  :: num_pts4
    sll_int32, dimension(1:4) :: gi
    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: vx
    sll_real64 :: vy
    sll_real64 :: vx_min
    sll_real64 :: vx_max
    sll_real64 :: vy_min
    sll_real64 :: vy_max
    sll_real64 :: v_thermal
    type(sll_collective_t), pointer :: col
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_real64 :: factor, factor2

    ! Find the local limits of the array to initialize. This could have been
    ! done as an initialization step, but it really does not matter, since
    ! the whole process is supposed to be done only once anyway.
    layout => init_obj%data_layout
    col    => get_layout_collective(layout)
    myrank =  sll_get_collective_rank( col )
    i_min = get_layout_i_min( layout, myrank )
    i_max = get_layout_i_max( layout, myrank )
    j_min = get_layout_j_min( layout, myrank )
    j_max = get_layout_j_max( layout, myrank )
    k_min = get_layout_k_min( layout, myrank )
    k_max = get_layout_k_max( layout, myrank )
    l_min = get_layout_l_min( layout, myrank )
    l_max = get_layout_l_max( layout, myrank )
    num_pts1 = i_max - i_min + 1
    num_pts2 = j_max - j_min + 1
    num_pts3 = k_max - k_min + 1
    num_pts4 = l_max - l_min + 1
    delta1 = init_obj%mesh_4d%delta_x1
    delta2 = init_obj%mesh_4d%delta_x2
    delta3 = init_obj%mesh_4d%delta_x3
    delta4 = init_obj%mesh_4d%delta_x4

    vx_min = init_obj%mesh_4d%x3_min
    vx_max = init_obj%mesh_4d%x3_max
    vy_min = init_obj%mesh_4d%x4_min
    vy_max = init_obj%mesh_4d%x4_max
    v_thermal = init_obj%v_thermal
    factor    = sll_pi/(8.0_f64*v_thermal**2)
    factor2   = 0.5_f64/v_thermal**2

    SLL_ASSERT( size(data_out,1) .ge. num_pts1 )
    SLL_ASSERT( size(data_out,2) .ge. num_pts2 )
    SLL_ASSERT( size(data_out,3) .ge. num_pts3 )
    SLL_ASSERT( size(data_out,4) .ge. num_pts4 )

    ! At this point nothing has been done regarding the data being
    ! cell centered... add this later.
             
    do l=1, num_pts4
       do k=1, num_pts3
          do j=1, num_pts2
             do i=1, num_pts1
                ! convert to global indices
                gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
                x  = gi(1)*delta1  ! danger: implicit xmin = 0
                y  = gi(2)*delta2  ! danger: implicit ymin = 0
                vx = vx_min + gi(3)*delta3
                vy = vy_min + gi(4)*delta4
                ! The following function could be defined outside and
                ! maybe stored in the object as a function pointer.
                data_out(i,j,k,l) = factor*sin(sll_pi*x)*sin(sll_pi*y)*&
                                           exp(-(vx**2+vy**2)*factor2)
             end do
          end do
       end do
    end do

  end subroutine 


!!$#ifdef STDF95
!!$  subroutine init_test_4d_par_f_of_4args( init_obj, data_out )
!!$    type(init_test_4d_par), intent(inout)      :: init_obj
!!$#else
!!$  subroutine compact_4d_field( init_obj, data_out )
!!$    class(init_test_4d_par), intent(inout)      :: init_obj
!!$#endif
!!$    sll_real64, dimension(:,:,:,:), intent(out) :: data_out
!!$    type(layout_4D), pointer                    :: layout
!!$    sll_int32  :: i
!!$    sll_int32  :: j
!!$    sll_int32  :: k
!!$    sll_int32  :: l
!!$    sll_int32  :: i_min
!!$    sll_int32  :: i_max
!!$    sll_int32  :: j_min
!!$    sll_int32  :: j_max
!!$    sll_int32  :: k_min
!!$    sll_int32  :: k_max
!!$    sll_int32  :: l_min
!!$    sll_int32  :: l_max
!!$    sll_int32  :: myrank
!!$    sll_int32  :: num_pts1
!!$    sll_int32  :: num_pts2
!!$    sll_int32  :: num_pts3
!!$    sll_int32  :: num_pts4
!!$    sll_int32, dimension(1:4) :: gi
!!$    sll_real64 :: x
!!$    sll_real64 :: y
!!$    sll_real64 :: vx
!!$    sll_real64 :: vy
!!$    sll_real64 :: vx_min
!!$    sll_real64 :: vx_max
!!$    sll_real64 :: vy_min
!!$    sll_real64 :: vy_max
!!$    sll_real64 :: v_thermal
!!$    type(sll_collective_t), pointer :: col
!!$    sll_real64 :: delta1
!!$    sll_real64 :: delta2
!!$    sll_real64 :: delta3
!!$    sll_real64 :: delta4
!!$    sll_real64 :: factor, factor2
!!$
!!$    ! Find the local limits of the array to initialize. This could have been
!!$    ! done as an initialization step, but it really does not matter, since
!!$    ! the whole process is supposed to be done only once anyway.
!!$    layout => init_obj%data_layout
!!$    col    => get_layout_collective(layout)
!!$    myrank =  sll_get_collective_rank( col )
!!$    i_min = get_layout_i_min( layout, myrank )
!!$    i_max = get_layout_i_max( layout, myrank )
!!$    j_min = get_layout_j_min( layout, myrank )
!!$    j_max = get_layout_j_max( layout, myrank )
!!$    k_min = get_layout_k_min( layout, myrank )
!!$    k_max = get_layout_k_max( layout, myrank )
!!$    l_min = get_layout_l_min( layout, myrank )
!!$    l_max = get_layout_l_max( layout, myrank )
!!$    num_pts1 = i_max - i_min + 1
!!$    num_pts2 = j_max - j_min + 1
!!$    num_pts3 = k_max - k_min + 1
!!$    num_pts4 = l_max - l_min + 1
!!$    delta1 = init_obj%mesh_4d%delta_x1
!!$    delta2 = init_obj%mesh_4d%delta_x2
!!$    delta3 = init_obj%mesh_4d%delta_x3
!!$    delta4 = init_obj%mesh_4d%delta_x4
!!$
!!$    vx_min = init_obj%mesh_4d%x3_min
!!$    vx_max = init_obj%mesh_4d%x3_max
!!$    vy_min = init_obj%mesh_4d%x4_min
!!$    vy_max = init_obj%mesh_4d%x4_max
!!$    v_thermal = init_obj%v_thermal
!!$    factor    = sll_pi/(8.0_f64*v_thermal**2)
!!$    factor2   = 0.5_f64/v_thermal**2
!!$
!!$    SLL_ASSERT( size(data_out,1) .ge. num_pts1 )
!!$    SLL_ASSERT( size(data_out,2) .ge. num_pts2 )
!!$    SLL_ASSERT( size(data_out,3) .ge. num_pts3 )
!!$    SLL_ASSERT( size(data_out,4) .ge. num_pts4 )
!!$
!!$    ! At this point nothing has been done regarding the data being
!!$    ! cell centered... add this later.
!!$             
!!$    do l=1, num_pts4
!!$       do k=1, num_pts3
!!$          do j=1, num_pts2
!!$             do i=1, num_pts1
!!$                ! convert to global indices
!!$                gi(:) = local_to_global_4D( layout, (/i,j,k,l/) )
!!$                x  = gi(1)*delta1  ! danger: implicit xmin = 0
!!$                y  = gi(2)*delta2  ! danger: implicit ymin = 0
!!$                vx = vx_min + gi(3)*delta3
!!$                vy = vy_min + gi(4)*delta4
!!$                ! The following function could be defined outside and
!!$                ! maybe stored in the object as a function pointer.
!!$                data_out(i,j,k,l) = factor*sin(sll_pi*x)*sin(sll_pi*y)*&
!!$                                           exp(-(vx**2+vy**2)*factor2)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$  end subroutine 



end module sll_test_4d_initializer
