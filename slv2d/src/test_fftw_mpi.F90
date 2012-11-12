program example
use mpi
use, intrinsic :: iso_c_binding
include 'fftw3-mpi.f03'
integer(C_INTPTR_T), parameter :: nx = 128
integer(C_INTPTR_T), parameter :: ny = 256
type(C_PTR) :: fw_plan, bw_plan, cdata
complex(C_DOUBLE_COMPLEX), pointer :: data(:,:)
integer(C_INTPTR_T) :: i, j, alloc_local, local_ny, local_j_offset
real(8) :: xmin, xmax, dx
real(8) :: ymin, ymax, dy
integer :: ierr, rank, size

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
print *, "Hello, world, I am ", rank, " of ", size

call fftw_mpi_init()

xmin = -2.; xmax = 2.
ymin = -2.; ymax = 2.

dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)

alloc_local = fftw_mpi_local_size_2d(ny, nx, MPI_COMM_WORLD, &
                                       local_ny, local_j_offset)
cdata = fftw_alloc_complex(alloc_local)
call c_f_pointer(cdata, data, [nx,local_ny])

fw_plan = fftw_mpi_plan_dft_2d(nx,ny,data,data, &
                               MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE)
bw_plan = fftw_mpi_plan_dft_2d(nx,ny,data,data, &
                               MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE)
do j = 1, local_ny
   do i = 1, nx
     data(i, j) = i*j + local_j_offset
   end do
end do



call fftw_mpi_execute_dft(fw_plan, data, data)
call fftw_mpi_execute_dft(bw_plan, data, data)

call fftw_destroy_plan(fw_plan)
call fftw_destroy_plan(bw_plan)
call fftw_free(cdata)
call fftw_mpi_cleanup()
call MPI_FINALIZE(ierr)

end program example
