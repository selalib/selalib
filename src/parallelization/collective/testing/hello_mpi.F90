program hello
   use mpi
   implicit none
   integer :: prank, psize, code
   integer :: nthreads
   integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
   real(8) :: tcpu1, tcpu2, time, err_l2, sum_l2, xp, yp
   real(8) :: x1, y1, r2, dtloc
   integer, dimension(MPI_STATUS_SIZE) :: statut
   integer                  :: comm2d
   integer, parameter        :: tag = 1111
   integer, dimension(8)     :: voisin
   integer, parameter        :: N = 1, S = 2, W = 3, E = 4
   integer, parameter        :: NW = 5, SW = 6, NE = 7, SE = 8
   integer, parameter        :: ndims = 2
   integer, dimension(ndims) :: dims, coords
   integer, dimension(ndims) :: coordsse, coordssw, coordsne, coordsnw
   logical                  :: reorder
   logical, dimension(ndims) :: periods
   integer                  :: nxp, nyp, mx, my
   integer                              :: comm
   integer                              :: error
   integer                              :: provided
   integer                              :: required

   comm = MPI_COMM_WORLD
   required = MPI_THREAD_MULTIPLE
   call MPI_Init_thread(required, provided, error)
   print *, 'hello, mpi', MPI_THREAD_MULTIPLE, provided
   call MPI_Comm_rank(MPI_COMM_WORLD, prank, code)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, psize, code)

   CALL MPI_BARRIER(MPI_COMM_WORLD, code)

   tcpu1 = MPI_WTIME()

   dims(:) = 0
   CALL MPI_DIMS_CREATE(psize, ndims, dims, code)
   nxp = dims(1)
   nyp = dims(2)

   periods(1) = .true.
   periods(2) = .true.
   reorder = .true.

   CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, comm2d, code)

   voisin(:) = MPI_PROC_NULL

   CALL MPI_CART_SHIFT(comm2d, 0, 1, voisin(N), voisin(S), code)
   CALL MPI_CART_SHIFT(comm2d, 1, 1, voisin(W), voisin(E), code)
   CALL MPI_COMM_RANK(comm2d, prank, code)
   CALL MPI_CART_COORDS(comm2d, prank, ndims, coords, code)

   if (voisin(N) /= MPI_PROC_NULL) then
      coordsnw(1) = coords(1) - 1; coordsnw(2) = coords(2) - 1
      coordsne(1) = coords(1) - 1; coordsne(2) = coords(2) + 1
      CALL MPI_CART_RANK(comm2d, coordsnw, voisin(NW), code)
      CALL MPI_CART_RANK(comm2d, coordsne, voisin(NE), code)
   end if

   if (voisin(S) /= MPI_PROC_NULL) then
      coordssw(1) = coords(1) + 1; coordssw(2) = coords(2) - 1
      coordsse(1) = coords(1) + 1; coordsse(2) = coords(2) + 1
      CALL MPI_CART_RANK(comm2d, coordssw, voisin(SW), code)
      CALL MPI_CART_RANK(comm2d, coordsse, voisin(SE), code)
   end if

   call MPI_BARRIER(MPI_COMM_WORLD, code)

#ifdef _OPENMP

   if (MPI_THREAD_MULTIPLE == provided) then

      !Fork a team of threads giving them their own copies of variables
!$OMP PARALLEL
      !Obtain thread number
!$OMP CRITICAL
      PRINT "('nodes:',i3,' rank:',i3,' threads:',i3,' thread_id:',i3)", &
         psize, prank, OMP_GET_NUM_THREADS(), OMP_GET_THREAD_NUM()
!$OMP END CRITICAL
!$OMP END PARALLEL

   end if

#endif // _OPENMP

   call MPI_FINALIZE(code)

end program hello
