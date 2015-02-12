program scatter
    use mpi
    implicit none

    integer, parameter :: gridsize = 6    ! size of array
    integer, parameter :: procgridsize = 2 ! size of process grid
    character, allocatable, dimension (:,:) :: global, local
    integer, dimension(procgridsize**2)   :: counts, displs
    integer, parameter    :: root = 0
    integer :: rank, comsize
    integer :: localsize
    integer :: i, j, row, col, ierr, p, charsize
    integer, dimension(2) :: sizes, subsizes, starts

    integer :: newtype, resizedtype
    integer, parameter :: tag = 1
    integer, dimension(MPI_STATUS_SIZE) :: rstatus
    integer(kind=MPI_ADDRESS_KIND) :: extent, begin

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    if (comsize /= procgridsize**2) then
        if (rank == root) then
            print *, 'Only works with np = ', procgridsize**2, ' for now.'
        endif
        call MPI_Finalize(ierr)
        stop
    endif

    localsize = gridsize/procgridsize
    allocate( local(localsize, localsize) )
    if (rank == root) then
        allocate( global(gridsize, gridsize) )
        forall( col=1:procgridsize, row=1:procgridsize )
            global((row-1)*localsize+1:row*localsize, &
                   (col-1)*localsize+1:col*localsize) = &
                    achar(ichar('0')+(row-1)+(col-1)*procgridsize)
        end forall

        print *, 'global array is: '
        do i=1,gridsize
            print *, global(i,:)
        enddo
    endif
    starts   = [0,0]
    sizes    = [gridsize, gridsize]
    subsizes = [localsize, localsize]

    call MPI_Type_create_subarray(2, sizes, subsizes, starts,        &
                                  MPI_ORDER_FORTRAN, MPI_CHARACTER,  &
                                  newtype, ierr)
    call MPI_Type_size(MPI_CHARACTER, charsize, ierr)
    extent = localsize*charsize
    begin  = 0
    call MPI_Type_create_resized(newtype, begin, extent, resizedtype, ierr)
    call MPI_Type_commit(resizedtype, ierr)

    counts = 1          ! we will send one of these new types to everyone
    forall( col=1:procgridsize, row=1:procgridsize )
       displs(1+(row-1)+procgridsize*(col-1)) = (row-1) + localsize*procgridsize*(col-1)
    endforall

    call MPI_Scatterv(global, counts, displs,   & ! proc i gets counts(i) types from displs(i)
            resizedtype,                        &
            local, localsize**2, MPI_CHARACTER, & ! I'm receiving localsize**2 chars
            root, MPI_COMM_WORLD, ierr)           !... from (root, MPI_COMM_WORLD)

    do p=1, comsize
        if (rank == p-1) then
            print *, 'Rank ', rank, ' received: '
            do i=1, localsize
                print *, local(i,:)
            enddo
        endif
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    enddo

    local = achar( ichar(local) + 1 )

    do p=1, comsize
        if (rank == p-1) then
            print *, 'Rank ', rank, ' sending: '
            do i=1, localsize
                print *, local(i,:)
            enddo
        endif
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    enddo

    call MPI_Gatherv( local, localsize**2, MPI_CHARACTER, & ! I'm sending localsize**2 chars
                      global, counts, displs, resizedtype,&
                      root, MPI_COMM_WORLD, ierr)

    if (rank == root) then
        print *, ' Root received: '
        do i=1,gridsize
            print *, global(i,:)
        enddo
    endif

    call MPI_Type_free(newtype,ierr)
    if (rank == root) deallocate(global)
    deallocate(local)
    call MPI_Finalize(ierr)

end program scatter
