program local_to_global

  use mpi
  implicit none

  integer  :: ncols, nrows

  integer, dimension(8)     :: neighbor
  integer, parameter        :: N =7, S =3, W =5, E =1
  integer, parameter        :: NW=6, SW=4, NE=8, SE=2
  integer, parameter        :: ndims = 2
  integer, dimension(ndims) :: dims
  integer, dimension(ndims) :: coords
  logical                   :: reorder
  logical,dimension(ndims)  :: periods

  integer                   :: statut(MPI_STATUS_SIZE)
  integer, parameter        :: tag = 1111

  integer, parameter :: gridsize_x = 30   ! size of array
  integer, parameter :: gridsize_y = 30   ! size of array
  character, allocatable, dimension (:,:) :: global
  character, allocatable, dimension (:,:) :: local
  integer, dimension(:), allocatable   :: counts, displs
  integer, parameter    :: proot = 0
  integer :: localsize_x, localsize_y
  integer :: row, col, ierr, iproc, charsize
  integer, dimension(ndims) :: sizes, subsizes, starts

  integer :: comm2d, psize, prank
  integer :: newtype, resizedtype
  integer, dimension(MPI_STATUS_SIZE) :: rstatus
  integer(kind=MPI_ADDRESS_KIND) :: extent, begin
  integer :: i
  real(8) :: tcpu1, tcpu2


  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, psize, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, prank, ierr)

  tcpu1 = MPI_WTIME()

  dims = 0
  CALL MPI_DIMS_CREATE(int(psize,4),ndims,dims,ierr)
  ncols = dims(1)
  nrows = dims(2)
  print*, 'ncols, nrows =', ncols, nrows
  periods(1) = .true.
  periods(2) = .true.
  reorder    = .true.

  CALL MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,periods,reorder,comm2d,ierr)

  neighbor(:) = MPI_PROC_NULL

  CALL MPI_CART_SHIFT(comm2d,0,1,neighbor(N),neighbor(S),ierr)
  CALL MPI_CART_SHIFT(comm2d,1,1,neighbor(W),neighbor(E),ierr)
  CALL MPI_COMM_RANK(comm2d,prank,ierr)
  CALL MPI_CART_COORDS(comm2d,prank,ndims,coords,ierr)
  CALL MPI_CART_RANK(comm2d,(/coords(1)-1,coords(2)+1/),neighbor(NW),ierr)
  CALL MPI_CART_RANK(comm2d,(/coords(1)+1,coords(2)+1/),neighbor(NE),ierr)
  CALL MPI_CART_RANK(comm2d,(/coords(1)-1,coords(2)-1/),neighbor(SW),ierr)
  CALL MPI_CART_RANK(comm2d,(/coords(1)+1,coords(2)-1/),neighbor(SE),ierr)

  call flush(6)
  print"('proc: ',i2,', coords: ',2i3,', neighbors: ',8i3)",prank,coords,neighbor
  call flush(6)
  call MPI_BARRIER(comm2d,ierr)

  if (mod(gridsize_x,nrows) == 0) then
     localsize_x = gridsize_x/nrows
  else
     print*, gridsize_x, nrows
     call MPI_FINALIZE(ierr); stop
  end if
  if (mod(gridsize_y,ncols) == 0) then
     localsize_y = gridsize_y/ncols
  else
     print*, gridsize_y, ncols
     call MPI_FINALIZE(ierr); stop
  end if

  localsize_y = gridsize_y/ncols
  allocate( local(localsize_x, localsize_y) )
  allocate( global(gridsize_x, gridsize_y) )
  global = char(48)
  local  = achar( ichar('0') + prank )

  starts   = [0,0]
  sizes    = [gridsize_x, gridsize_y]
  subsizes = [localsize_x, localsize_y]

  print*, 'create the subarray type'
  call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts,        &
                                MPI_ORDER_FORTRAN, MPI_CHARACTER,  &
                                newtype, ierr)

  print*, 'resize the subarray type'
  call MPI_TYPE_SIZE(MPI_CHARACTER, charsize, ierr)
  extent = localsize_x*charsize
  begin  = 0
  call MPI_TYPE_CREATE_RESIZED(newtype, begin, extent, resizedtype, ierr)
  call MPI_TYPE_COMMIT(resizedtype, ierr)

  allocate(counts(ncols*nrows))
  allocate(displs(ncols*nrows))

  print*, 'compute displacements'
  counts = 1          
  forall( col=1:ncols, row=1:nrows )
     displs((col-1)*ncols+row) = (row-1) + localsize_y*ncols*(col-1)
  endforall
 
  call MPI_AllGatherv( local, localsize_x*localsize_y, MPI_CHARACTER, & 
                       global, counts, displs, resizedtype,&
                       comm2d, ierr)

  do iproc=0, psize-1
     if (iproc == prank) then
     print *, ' Rank ', iproc, ' received: '
     do i=1,gridsize_x
        print *, global(i,:)
     enddo
     end if
     call MPI_Barrier(MPI_COMM_WORLD, ierr)
  enddo

  call MPI_Type_free(newtype,ierr)
  deallocate(global)
  deallocate(local)


  tcpu2 = MPI_WTIME()
  if (prank == proot) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize

  call MPI_FINALIZE(ierr)
  stop

end program local_to_global
