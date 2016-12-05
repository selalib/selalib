 sll_int32, dimension(:), pointer     :: send_counts
 sll_int32, dimension(:), pointer     :: send_displs
 sll_int32, dimension(:), pointer     :: recv_counts
 sll_int32, dimension(:), pointer     :: recv_displs
 type(sll_t_collective_t), pointer    :: col
 sll_int32                            :: col_sz
 sll_int32, dimension(:), allocatable :: lowest_color
 sll_int32, dimension(:), allocatable :: colors
 sll_int32, dimension(:), allocatable :: colors_copy
 sll_int32                            :: ierr
 sll_int32                            :: my_rank
 sll_int32                            :: i
 type(sll_t_collective_t), pointer    :: new_collective
 sll_int32                            :: new_col_sz
 sll_int32, dimension(:), pointer     :: new_send_counts
 sll_int32, dimension(:), pointer     :: new_send_displs
 sll_int32, dimension(:), pointer     :: new_recv_counts
 sll_int32, dimension(:), pointer     :: new_recv_displs
 sll_int32                            :: new_i
 sll_int32                            :: my_color
 sll_int32                            :: exchange_size_s
 sll_int32                            :: exchange_size_r
 logical, dimension(1:1)              :: is_uniform_local
 logical, dimension(1:1)              :: is_uniform_collective
 sll_int32                            :: new_sdisp
 sll_int32                            :: new_rdisp
 
 col         => plan%collective
 col_sz      = sll_f_get_collective_size( col )
 my_rank     = sll_f_get_collective_rank( col )
 send_counts => plan%send_counts
 send_displs => plan%send_displs
 recv_counts => plan%recv_counts
 recv_displs => plan%recv_displs
 SLL_ALLOCATE( lowest_color(1), ierr )
 lowest_color(1) = 0
 SLL_ALLOCATE( colors(0:col_sz-1), ierr )
 colors(:) = 0
 SLL_ALLOCATE( colors_copy(0:col_sz-1), ierr )
 colors_copy(:) = 0
 
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
 call sll_o_collective_allgather(col,lowest_color,1,colors(0:col_sz-1),1)
 do
    ! Load the copy
    colors_copy(0:col_sz-1) = colors(0:col_sz-1)
    ! Find the lowest rank with which this process communicates
    do i=0,col_sz-1
       if( (send_counts(i) .ne. 0) .or. (recv_counts(i) .ne. 0) ) then
          if( colors(i) .lt. lowest_color(1) ) then
             lowest_color(1) = colors(i)
          end if ! else, nothing, as lowest_color is already my_rank
       end if
    end do
    ! Gather the information from all processes
    call sll_o_collective_allgather(col,lowest_color,1,colors(0:col_sz-1),1)
    if(arrays_are_equal(colors, colors_copy, col_sz)) then
       exit
    end if
 end do
 ! The results can now be used as the color for a collective-splitting 
 ! operation.
 new_collective => sll_f_new_collective( col, colors(my_rank), my_rank )
 new_col_sz     = sll_f_get_collective_size( new_collective )
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

 ! First we initialize "is_uniform_local" to .true. It will become false if
 ! one count does not agree
 is_uniform_local = .true.
 exchange_size_s = plan%send_counts(0)
 do i=0,new_col_sz-1
    if(plan%send_counts(i) .eq. exchange_size_s) then
       is_uniform_local(1) = is_uniform_local(1) .and. .true.
    else ! plan%send_counts(i) is different than the first value
       is_uniform_local(1) = is_uniform_local(1) .and. .false.
       exit
    end if
 end do
 exchange_size_r = plan%recv_counts(0)
 do i=0,new_col_sz-1
    if(plan%recv_counts(i) .eq. exchange_size_r) then
       is_uniform_local(1) = is_uniform_local(1) .and. .true.
    else
       is_uniform_local(1) = is_uniform_local(1) .and. .false.
       exit
    end if
 end do
 ! Use a reduction operation to find out if this result is shared with 
 ! the other processes in the collective. Hmmm... look at this slightly
 ! disastrous occurrence: the MPI reduction operation MPI_LAND got out of
 ! the cage... this needs to be addressed.
 call sll_o_collective_allreduce(plan%collective,is_uniform_local(:),1,MPI_LAND, is_uniform_collective(:) )
 
 ! This flag will be used for an optimized call in apply_remap_plan()
 plan%is_uniform = is_uniform_collective(1)

end subroutine
