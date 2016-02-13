!> @ingroup parallel_utilities

module sll_m_buffer_loader_utilities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_t_collective_t, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size

  use sll_m_remapper, only: &
    sll_o_get_layout_collective, &
    sll_o_get_layout_i_max, &
    sll_o_get_layout_i_min, &
    sll_o_get_layout_j_max, &
    sll_o_get_layout_j_min, &
    sll_o_get_layout_k_max, &
    sll_o_get_layout_k_min, &
    sll_t_layout_2d, &
    sll_t_layout_3d

  implicit none

  public :: &
    sll_s_compute_displacements_array_2d, &
    sll_s_load_buffer32_2d, &
    sll_s_load_buffer_2d, &
    sll_s_load_buffer_3d, &
    sll_f_receive_counts_array_2d, &
    sll_s_unload_buffer_2d, &
    sll_s_unload_buffer_3d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

  !> the objective of this subroutine is to help to prepare data
  !> for gatherv operations for example
  subroutine sll_s_compute_displacements_array_2d( layout, collective_size, disps )
    type(sll_t_layout_2d), pointer                           :: layout
    sll_int32, intent(in)                              :: collective_size
    sll_int32, dimension(collective_size), intent(out) :: disps
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    !sll_int32 :: i,j
    sll_int32 :: counter
    sll_int32 :: rank

    counter  = 0
    disps(1) = counter
    do rank=1,collective_size-1
       imin = sll_o_get_layout_i_min( layout, rank-1 )
       imax = sll_o_get_layout_i_max( layout, rank-1 )
       jmin = sll_o_get_layout_j_min( layout, rank-1 )
       jmax = sll_o_get_layout_j_max( layout, rank-1 )
       size_i      = imax - imin + 1
       size_j      = jmax - jmin + 1
       counter     = counter + size_i*size_j
       disps(rank+1) = counter
    end do
  end subroutine sll_s_compute_displacements_array_2d

  
  !> this subroutine loads 2D array into a 1D
  !> buffer for use in collective communications
  subroutine sll_s_load_buffer_2d( layout, data, buffer )
    type(sll_t_layout_2d), pointer   :: layout
    sll_real64, dimension(:,:), intent(in) :: data
    sll_real64, dimension(:),  intent(out) :: buffer
    sll_int32 :: myrank
    sll_int32 :: data_size
    !sll_int32 :: send_size
    type(sll_t_collective_t), pointer :: col
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    sll_int32 :: i,j
    sll_int32 :: counter

    col => sll_o_get_layout_collective( layout )
    myrank = sll_f_get_collective_rank( col )
    data_size = size(data,1)*size(data,2)

    imin = sll_o_get_layout_i_min( layout, myrank )
    imax = sll_o_get_layout_i_max( layout, myrank )
    jmin = sll_o_get_layout_j_min( layout, myrank )
    jmax = sll_o_get_layout_j_max( layout, myrank )
    size_i = imax - imin + 1
    size_j = jmax - jmin + 1

    counter=0
    do j=1,size_j
       do i=1,size_i
          counter = counter + 1
          buffer(counter) = data(i,j)
       end do
    end do

  end subroutine sll_s_load_buffer_2d

  subroutine sll_s_load_buffer32_2d( layout, data, buffer )
    type(sll_t_layout_2d), pointer   :: layout
    sll_real32, dimension(:,:), intent(in) :: data
    sll_real32, dimension(:),  intent(out) :: buffer
    sll_int32 :: myrank
    sll_int32 :: data_size
    !sll_int32 :: send_size
    type(sll_t_collective_t), pointer :: col
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    sll_int32 :: i,j
    sll_int32 :: counter

    col => sll_o_get_layout_collective( layout )
    myrank = sll_f_get_collective_rank( col )
    data_size = size(data,1)*size(data,2)

    imin = sll_o_get_layout_i_min( layout, myrank )
    imax = sll_o_get_layout_i_max( layout, myrank )
    jmin = sll_o_get_layout_j_min( layout, myrank )
    jmax = sll_o_get_layout_j_max( layout, myrank )
    size_i = imax - imin + 1
    size_j = jmax - jmin + 1

    counter=0
    do j=1,size_j
       do i=1,size_i
          counter = counter + 1
          buffer(counter) = data(i,j)
       end do
    end do

  end subroutine sll_s_load_buffer32_2d



  !> this subroutine loads 3D array into a 1D
  !> buffer for use in collective communications
  subroutine sll_s_load_buffer_3d( layout, data, buffer )
    type(sll_t_layout_3d), pointer   :: layout
    sll_real64, dimension(:,:,:), intent(in) :: data
    sll_real64, dimension(:),  intent(out) :: buffer
    sll_int32 :: myrank
    sll_int32 :: data_size
    !sll_int32 :: send_size
    type(sll_t_collective_t), pointer :: col
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: kmin, kmax
    sll_int32 :: size_i, size_j, size_k
    sll_int32 :: i,j,k
    sll_int32 :: counter

    col => sll_o_get_layout_collective( layout )
    myrank = sll_f_get_collective_rank( col )
    data_size = size(data,1)*size(data,2)*size(data,3)
    imin = sll_o_get_layout_i_min( layout, myrank )
    imax = sll_o_get_layout_i_max( layout, myrank )
    jmin = sll_o_get_layout_j_min( layout, myrank )
    jmax = sll_o_get_layout_j_max( layout, myrank )
    kmin = sll_o_get_layout_k_min( layout, myrank )
    kmax = sll_o_get_layout_k_max( layout, myrank )
    size_i = imax - imin + 1
    size_j = jmax - jmin + 1
    size_k = kmax - kmin + 1
    counter=0
    do k=1,size_k
      do j=1,size_j
        do i=1,size_i
          counter = counter + 1
          buffer(counter) = data(i,j,k)
        end do
      end do
    end do
  end subroutine sll_s_load_buffer_3d


  !> this function is also a helper for collective routines for example gatherv
  !> should be change into subroutine 
  function sll_f_receive_counts_array_2d( layout, n ) result(rc)
    type(sll_t_layout_2d), pointer :: layout
    sll_int32, intent(in)    :: n
    sll_int32, dimension(n) :: rc
    sll_int32 :: i
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j

    do i=0,n-1
       imin = sll_o_get_layout_i_min( layout, i )
       imax = sll_o_get_layout_i_max( layout, i )
       jmin = sll_o_get_layout_j_min( layout, i )
       jmax = sll_o_get_layout_j_max( layout, i )
       size_i = imax - imin + 1
       size_j = jmax - jmin + 1
       rc(i+1)  = size_i*size_j
    end do
  end function sll_f_receive_counts_array_2d

  !this routine takes 1D array and stores it in a 2D array

  subroutine sll_s_unload_buffer_2d( layout, buffer, data )
    type(sll_t_layout_2d), pointer                :: layout
    sll_real64, dimension(:,:), intent(out) :: data
    sll_real64, dimension(:), intent(in)    :: buffer
    sll_int32 :: col_sz
    type(sll_t_collective_t), pointer :: col
    sll_int32 :: i, j
    sll_int32 :: box
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    !sll_int32 :: size_i, size_j
    sll_int32 :: pos            ! position in buffer
    col => sll_o_get_layout_collective( layout )
    col_sz = sll_f_get_collective_size( col )

    ! loop over all the boxes in the layout and fill the data array by chunks.
    pos = 1
    do box=0,col_sz-1
       imin = sll_o_get_layout_i_min( layout, box )
       imax = sll_o_get_layout_i_max( layout, box )
       jmin = sll_o_get_layout_j_min( layout, box )
       jmax = sll_o_get_layout_j_max( layout, box )
       ! this will fill the data array in whatever order that the boxes
       ! are ordered.
       do j=jmin,jmax
          do i=imin,imax
             data(i,j) = buffer(pos)
             pos       = pos + 1
          end do
       end do
    end do
  end subroutine sll_s_unload_buffer_2d


  !> this routine takes 1D array and stores it in a 3D array
  !> warning definition changes are done wrt sll_s_unload_buffer_2d
  subroutine sll_s_unload_buffer_3d( layout, buffer, data )
    type(sll_t_layout_3d), pointer                :: layout
    sll_real64, dimension(:,:,:), intent(out) :: data
    sll_real64, dimension(:), intent(in)    :: buffer
    sll_int32 :: col_sz
    type(sll_t_collective_t), pointer :: col
    sll_int32 :: i, j, k
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: kmin, kmax
    sll_int32 :: size_i, size_j, size_k
    sll_int32 :: pos            ! position in buffer
    sll_int32 :: myrank
    col => sll_o_get_layout_collective( layout )
    col_sz = sll_f_get_collective_size( col )
    myrank = sll_f_get_collective_rank( col )
    ! loop over all the boxes in the layout and fill the data array by chunks.
    pos = 1
    !do box=0,col_sz-1
       imin = sll_o_get_layout_i_min( layout, myrank )
       imax = sll_o_get_layout_i_max( layout, myrank )
       jmin = sll_o_get_layout_j_min( layout, myrank )
       jmax = sll_o_get_layout_j_max( layout, myrank )
       kmin = sll_o_get_layout_k_min( layout, myrank )
       kmax = sll_o_get_layout_k_max( layout, myrank )
    size_i = imax - imin + 1
    size_j = jmax - jmin + 1
    size_k = kmax - kmin + 1
       !print *,'kmin=',kmin,kmax
       ! this will fill the data array in whatever order that the boxes
       ! are ordered.
       do k=1,size_k !k=kmin,kmax
         do j=1,size_j !j=jmin,jmax
           do i=1,size_i !i=imin,imax
             data(i,j,k) = buffer(pos)
             pos       = pos + 1
           end do
         end do
       end do  
    !end do
  end subroutine sll_s_unload_buffer_3d





end module sll_m_buffer_loader_utilities
