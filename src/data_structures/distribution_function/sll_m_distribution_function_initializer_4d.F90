module sll_m_distribution_function_initializer_4d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_collective, only: &
      sll_t_collective_t, &
      sll_f_get_collective_rank

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_remapper, only: &
      sll_o_get_layout_collective, &
      sll_o_get_layout_i_max, &
      sll_o_get_layout_i_min, &
      sll_o_get_layout_j_max, &
      sll_o_get_layout_j_min, &
      sll_o_get_layout_k_max, &
      sll_o_get_layout_k_min, &
      sll_o_get_layout_l_max, &
      sll_o_get_layout_l_min, &
      sll_t_layout_4d, &
      sll_o_local_to_global

   use sll_m_scalar_field_initializers_base, only: &
      sll_c_scalar_field_4d_initializer_base

   implicit none

   public :: &
      sll_t_init_test_4d_par, &
      sll_s_load_test_4d_initializer, &
      sll_f_new_cartesian_4d_mesh, &
      sll_t_simple_cartesian_4d_mesh

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

   type, extends(sll_c_scalar_field_4d_initializer_base) :: sll_t_init_test_4d_par
      sll_real64 :: v_thermal
      type(sll_t_layout_4d), pointer :: data_layout
      type(sll_t_simple_cartesian_4d_mesh), pointer :: mesh_4d
   contains
      procedure, pass(init_obj) :: initialize => sll_s_load_test_4d_initializer
      procedure, pass(init_obj) :: f_of_4args => compact_4d_field
   end type sll_t_init_test_4d_par

   ! This has to end up somewhere else, if it is to stay in the library
   type sll_t_simple_cartesian_4d_mesh
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
   end type sll_t_simple_cartesian_4d_mesh

contains

   function sll_f_new_cartesian_4d_mesh( &
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
      x4_max)

      type(sll_t_simple_cartesian_4d_mesh), pointer :: sll_f_new_cartesian_4d_mesh
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

      SLL_ALLOCATE(sll_f_new_cartesian_4d_mesh, ierr)
      sll_f_new_cartesian_4d_mesh%num_cells1 = num_cells1
      sll_f_new_cartesian_4d_mesh%num_cells2 = num_cells2
      sll_f_new_cartesian_4d_mesh%num_cells3 = num_cells3
      sll_f_new_cartesian_4d_mesh%num_cells4 = num_cells4
      sll_f_new_cartesian_4d_mesh%x1_min = x1_min
      sll_f_new_cartesian_4d_mesh%x1_max = x1_max
      sll_f_new_cartesian_4d_mesh%x2_min = x2_min
      sll_f_new_cartesian_4d_mesh%x2_max = x2_max
      sll_f_new_cartesian_4d_mesh%x3_min = x3_min
      sll_f_new_cartesian_4d_mesh%x3_max = x3_max
      sll_f_new_cartesian_4d_mesh%x4_min = x4_min
      sll_f_new_cartesian_4d_mesh%x4_max = x4_max
      sll_f_new_cartesian_4d_mesh%delta_x1 = (x1_max - x1_min)/real(num_cells1, f64)
      sll_f_new_cartesian_4d_mesh%delta_x2 = (x2_max - x2_min)/real(num_cells2, f64)
      sll_f_new_cartesian_4d_mesh%delta_x3 = (x3_max - x3_min)/real(num_cells3, f64)
      sll_f_new_cartesian_4d_mesh%delta_x4 = (x4_max - x4_min)/real(num_cells4, f64)
   end function sll_f_new_cartesian_4d_mesh

   subroutine delete_cartesian_4d_mesh(mesh)
      type(sll_t_simple_cartesian_4d_mesh), pointer :: mesh
      sll_int32 :: ierr
      SLL_DEALLOCATE(mesh, ierr)
   end subroutine delete_cartesian_4d_mesh

   subroutine sll_s_load_test_4d_initializer( &
      init_obj, &
      data_position, &
      mesh_4d, &
      v_thermal, &
      layout)

      class(sll_t_init_test_4d_par), intent(inout)  :: init_obj
      sll_int32                               :: data_position
      sll_real64, intent(in)                  :: v_thermal
      type(sll_t_layout_4d), pointer                :: layout
      type(sll_t_simple_cartesian_4d_mesh), pointer :: mesh_4d

      if (.not. associated(layout)) then
         print *, 'initialize_test_4d(): ERROR, passed layout is not initialized'
         stop
      end if

      if (.not. associated(mesh_4d)) then
         print *, 'initialize_test_4d(): ERROR, passed 4d mesh is not initialized'
         stop
      end if

      init_obj%data_position = data_position
      init_obj%mesh_4d => mesh_4d
      init_obj%v_thermal = v_thermal
      init_obj%data_layout => layout
   end subroutine

   subroutine compact_4d_field(init_obj, data_out)
      class(sll_t_init_test_4d_par), intent(inout)      :: init_obj
      sll_real64, dimension(:, :, :, :), intent(out) :: data_out
      type(sll_t_layout_4d), pointer                    :: layout
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
      type(sll_t_collective_t), pointer :: col
      sll_real64 :: delta1
      sll_real64 :: delta2
      sll_real64 :: delta3
      sll_real64 :: delta4
      sll_real64 :: factor, factor2

      ! Find the local limits of the array to initialize. This could have been
      ! done as an initialization step, but it really does not matter, since
      ! the whole process is supposed to be done only once anyway.
      layout => init_obj%data_layout
      col => sll_o_get_layout_collective(layout)
      myrank = sll_f_get_collective_rank(col)
      i_min = sll_o_get_layout_i_min(layout, myrank)
      i_max = sll_o_get_layout_i_max(layout, myrank)
      j_min = sll_o_get_layout_j_min(layout, myrank)
      j_max = sll_o_get_layout_j_max(layout, myrank)
      k_min = sll_o_get_layout_k_min(layout, myrank)
      k_max = sll_o_get_layout_k_max(layout, myrank)
      l_min = sll_o_get_layout_l_min(layout, myrank)
      l_max = sll_o_get_layout_l_max(layout, myrank)
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
      factor = sll_p_pi/(8.0_f64*v_thermal**2)
      factor2 = 0.5_f64/v_thermal**2

      SLL_ASSERT(size(data_out, 1) .ge. num_pts1)
      SLL_ASSERT(size(data_out, 2) .ge. num_pts2)
      SLL_ASSERT(size(data_out, 3) .ge. num_pts3)
      SLL_ASSERT(size(data_out, 4) .ge. num_pts4)

      ! At this point nothing has been done regarding the data being
      ! cell centered... add this later.

      do l = 1, num_pts4
         do k = 1, num_pts3
            do j = 1, num_pts2
               do i = 1, num_pts1
                  ! convert to global indices
                  gi(:) = sll_o_local_to_global(layout, (/i, j, k, l/))
                  x = gi(1)*delta1  ! danger: implicit xmin = 0
                  y = gi(2)*delta2  ! danger: implicit ymin = 0
                  vx = vx_min + gi(3)*delta3
                  vy = vy_min + gi(4)*delta4
                  ! The following function could be defined outside and
                  ! maybe stored in the object as a function pointer.
                  data_out(i, j, k, l) = factor*sin(sll_p_pi*x)*sin(sll_p_pi*y)* &
                                         exp(-(vx**2 + vy**2)*factor2)
               end do
            end do
         end do
      end do

   end subroutine

!!$  subroutine compact_4d_field( init_obj, data_out )
!!$    class(sll_t_init_test_4d_par), intent(inout)      :: init_obj
!!$    sll_real64, dimension(:,:,:,:), intent(out) :: data_out
!!$    type(sll_t_layout_4d), pointer                    :: layout
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
!!$    type(sll_t_collective_t), pointer :: col
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
!!$    col    => sll_o_get_layout_collective(layout)
!!$    myrank =  sll_f_get_collective_rank( col )
!!$    i_min = sll_o_get_layout_i_min( layout, myrank )
!!$    i_max = sll_o_get_layout_i_max( layout, myrank )
!!$    j_min = sll_o_get_layout_j_min( layout, myrank )
!!$    j_max = sll_o_get_layout_j_max( layout, myrank )
!!$    k_min = sll_o_get_layout_k_min( layout, myrank )
!!$    k_max = sll_o_get_layout_k_max( layout, myrank )
!!$    l_min = sll_o_get_layout_l_min( layout, myrank )
!!$    l_max = sll_o_get_layout_l_max( layout, myrank )
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
!!$    factor    = sll_p_pi/(8.0_f64*v_thermal**2)
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
!!$                gi(:) = sll_o_local_to_global( layout, (/i,j,k,l/) )
!!$                x  = gi(1)*delta1  ! danger: implicit xmin = 0
!!$                y  = gi(2)*delta2  ! danger: implicit ymin = 0
!!$                vx = vx_min + gi(3)*delta3
!!$                vy = vy_min + gi(4)*delta4
!!$                ! The following function could be defined outside and
!!$                ! maybe stored in the object as a function pointer.
!!$                data_out(i,j,k,l) = factor*sin(sll_p_pi*x)*sin(sll_p_pi*y)*&
!!$                                           exp(-(vx**2+vy**2)*factor2)
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$  end subroutine

end module sll_m_distribution_function_initializer_4d
