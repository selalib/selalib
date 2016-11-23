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
 lowest_color(1) = my_rank
 call sll_o_collective_allgather(col,lowest_color,1,colors(0:col_sz-1),1)
 do
    colors_copy(0:col_sz-1) = colors(0:col_sz-1)
    do i=0,col_sz-1
       if( (send_counts(i) .ne. 0) .or. (recv_counts(i) .ne. 0) ) then
          if( colors(i) .lt. lowest_color(1) ) then
             lowest_color(1) = colors(i)
          end if
       end if
    end do
    call sll_o_collective_allgather(col,lowest_color,1,colors(0:col_sz-1),1)
    if(arrays_are_equal(colors, colors_copy, col_sz)) then
       exit
    end if
 end do
 new_collective => sll_f_new_collective( col, colors(my_rank), my_rank )
 new_col_sz     = sll_f_get_collective_size( new_collective )
 SLL_ALLOCATE( new_send_counts(0:new_col_sz-1), ierr )
 SLL_ALLOCATE( new_send_displs(0:new_col_sz-1), ierr )
 SLL_ALLOCATE( new_recv_counts(0:new_col_sz-1), ierr )
 SLL_ALLOCATE( new_recv_displs(0:new_col_sz-1), ierr )
 SLL_ALLOCATE( new_send_boxes( 0:new_col_sz-1), ierr )
 SLL_ALLOCATE( new_recv_boxes( 0:new_col_sz-1), ierr )
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
 
 is_uniform_local = .true.
 
 exchange_size_s = plan%send_counts(0)
 do i=0,new_col_sz-1
    if(plan%send_counts(i) .eq. exchange_size_s) then
       is_uniform_local(1) = is_uniform_local(1) .and. .true.
    else
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
 call sll_o_collective_allreduce(plan%collective,is_uniform_local(:),1,MPI_LAND, is_uniform_collective(:) )
 plan%is_uniform = is_uniform_collective(1)

end subroutine
