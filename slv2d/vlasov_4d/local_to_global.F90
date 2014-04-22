program local_to_global

  use mpi
  implicit none

  integer  :: nx, ny

  integer, dimension(8)     :: neighbor
  integer, parameter        :: N =7, S =3, W =5, E =1
  integer, parameter        :: NW=6, SW=4, NE=8, SE=2
  integer, parameter        :: ndims = 2
  integer, dimension(ndims) :: dims
  integer, dimension(ndims) :: coords
  logical                   :: reorder
  logical,dimension(ndims)  :: periods

  integer, parameter        :: tag = 1111

  integer, parameter :: gridsize_x = 12   ! size of array
  integer, parameter :: gridsize_y = 24   ! size of array
  character, allocatable, dimension (:,:) :: global
  character, allocatable, dimension (:,:) :: local
  integer, dimension(:), allocatable   :: counts, displs
  integer, parameter    :: proot = 0
  integer :: localsize_x, localsize_y
  integer :: i, j, k,  ierr, iproc, charsize
  integer, dimension(ndims) :: sizes, subsizes, starts

  integer :: comm2d, psize, prank
  integer :: newtype, resizedtype
  integer(kind=MPI_ADDRESS_KIND) :: extent, begin
  real(8) :: tcpu1, tcpu2


  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, psize, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, prank, ierr)

  tcpu1 = MPI_WTIME()

  dims = 0
  CALL MPI_DIMS_CREATE(int(psize,4),ndims,dims,ierr)
  nx = dims(1)
  ny = dims(2)
  print*, 'nx, ny =', nx, ny
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

  if (mod(gridsize_x,nx) == 0) then
     localsize_x = gridsize_x/nx
  else
     print*, 'nx=', nx, 'gridsize_x', gridsize_x
     call MPI_FINALIZE(ierr); stop
  end if
  if (mod(gridsize_y,ny) == 0) then
     localsize_y = gridsize_y/ny
  else
     print*, 'ny=', ny, 'gridsize_y', gridsize_y
     call MPI_FINALIZE(ierr); stop
  end if

  allocate( local(localsize_x, localsize_y) )
  allocate( global(gridsize_x, gridsize_y) )
  global = char(48)
  local  = achar( ichar('A') + prank )

  print*, gridsize_x, gridsize_y
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

  allocate(counts(nx*ny))
  allocate(displs(nx*ny))

  print*, 'compute displacements'
  counts = 1          
  displs = 0
  k = 0
  do j=1,ny
     do i=1,nx
        k=k+1
        displs(k) = j-1+(i-1)*localsize_y*nx
     end do
  end do

 
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
