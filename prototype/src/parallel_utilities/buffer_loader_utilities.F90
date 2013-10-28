module sll_buffer_loader_utilities_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

use sll_collective
use sll_remapper

implicit none

contains

  ! the objective of this subroutine is to help to prepare data
  ! for gatherv operations for example

  subroutine compute_displacements_array_2d( layout, collective_size, disps )
    type(layout_2D), pointer                           :: layout
    sll_int32, intent(in)                              :: collective_size
    sll_int32, dimension(collective_size), intent(out) :: disps
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    sll_int32 :: i,j
    sll_int32 :: counter
    sll_int32 :: rank

    counter  = 0
    disps(1) = counter
    do rank=1,collective_size-1
       imin = get_layout_i_min( layout, rank-1 )
       imax = get_layout_i_max( layout, rank-1 )
       jmin = get_layout_j_min( layout, rank-1 )
       jmax = get_layout_j_max( layout, rank-1 )
       size_i      = imax - imin + 1
       size_j      = jmax - jmin + 1
       counter     = counter + size_i*size_j
       disps(rank+1) = counter
    end do
  end subroutine compute_displacements_array_2d

  
  ! this subroutine loads 2D array into a 1D
  ! buffer for use in collective communications


  subroutine load_buffer_2d( layout, data, buffer )
    type(layout_2D), pointer   :: layout
    sll_real64, dimension(:,:), intent(in) :: data
    sll_real64, dimension(:),  intent(out) :: buffer
    sll_int32 :: myrank
    sll_int32 :: data_size
    sll_int32 :: send_size
    type(sll_collective_t), pointer :: col
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    sll_int32 :: i,j
    sll_int32 :: counter

    col => get_layout_collective( layout )
    myrank = sll_get_collective_rank( col )
    data_size = size(data,1)*size(data,2)
!print *, 'data1: ', size(data,1), 'data2:', size(data,2)
    imin = get_layout_i_min( layout, myrank )
    imax = get_layout_i_max( layout, myrank )
    jmin = get_layout_j_min( layout, myrank )
    jmax = get_layout_j_max( layout, myrank )
    size_i = imax - imin + 1
    size_j = jmax - jmin + 1

!!$    if( data_size .ne. size_i*size_j ) then
!!$       print *, 'function load_buffer():'
!!$       print *, 'size(data) = ', size(data,1), size(data,2)
!!$       print *, 'warning from rank ', myrank
!!$       print *, 'there seems to be a discrepancy between the data size ', &
!!$            'passed and the size declared in the layout.'
!!$       print *, 'data size = ', data_size, 'size from layout = ', size_i*size_j
!!$    end if
!!$print *, 'size_j', size_j, 'size_i', size_i
    counter=0
    do j=1,size_j
       do i=1,size_i
          counter = counter + 1
          buffer(counter) = data(i,j)
       end do
    end do

  end subroutine load_buffer_2d


  ! this function is also a helper for collective routines for example gatherv
  ! should be change into subroutine 
  function receive_counts_array_2d( layout, n ) result(rc)
    type(layout_2D), pointer :: layout
    sll_int32, intent(in)    :: n
    sll_int32, dimension(n) :: rc
    sll_int32 :: i
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j

    do i=0,n-1
       imin = get_layout_i_min( layout, i )
       imax = get_layout_i_max( layout, i )
       jmin = get_layout_j_min( layout, i )
       jmax = get_layout_j_max( layout, i )
       size_i = imax - imin + 1
       size_j = jmax - jmin + 1
       rc(i+1)  = size_i*size_j
    end do
  end function receive_counts_array_2d

  !this routine takes 1D array and stores it in a 2D array

  subroutine unload_buffer_2d( layout, buffer, data )
    type(layout_2D), pointer                :: layout
    sll_real64, dimension(:,:), intent(out) :: data
    sll_real64, dimension(:), intent(in)    :: buffer
    sll_int32 :: col_sz
    type(sll_collective_t), pointer :: col
    sll_int32 :: i, j
    sll_int32 :: box
    sll_int32 :: imin, imax
    sll_int32 :: jmin, jmax
    sll_int32 :: size_i, size_j
    sll_int32 :: pos            ! position in buffer
    col => get_layout_collective( layout )
    col_sz = sll_get_collective_size( col )

    ! loop over all the boxes in the layout and fill the data array by chunks.
    pos = 1
    do box=0,col_sz-1
       imin = get_layout_i_min( layout, box )
       imax = get_layout_i_max( layout, box )
       jmin = get_layout_j_min( layout, box )
       jmax = get_layout_j_max( layout, box )
       ! this will fill the data array in whatever order that the boxes
       ! are ordered.
       do j=jmin,jmax
          do i=imin,imax
             data(i,j) = buffer(pos)
             pos       = pos + 1
          end do
       end do
    end do
  end subroutine unload_buffer_2d




end module sll_buffer_loader_utilities_module