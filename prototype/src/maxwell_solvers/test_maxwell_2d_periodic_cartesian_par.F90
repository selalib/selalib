#define MPI_MASTER 0
program test_maxwell_2d_periodic_cart_par

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use remapper
  use numeric_constants
  use sll_maxwell_2d_periodic_cartesian_par
  use sll_collective
  use hdf5
  use sll_hdf5_io_parallel, only: sll_hdf5_file_create, &
                                  sll_hdf5_write_array, &
                                  sll_hdf5_file_close
  implicit none

  sll_int32                                    :: ncx, ncy
  sll_int32                                    :: nx_loc, ny_loc
  sll_int32                                    :: ierr
  sll_real64                                   :: Lx, Ly
  sll_real64                                   :: dx, dy
  sll_real64                                   :: x, y
  sll_real64, dimension(:,:), allocatable      :: ex
  sll_real64, dimension(:,:), allocatable      :: ey
  sll_real64, dimension(:,:), allocatable      :: bz
  sll_real64, dimension(:,:), allocatable      :: bz_an
  sll_int32                                    :: i, j
  type (maxwell_2d_periodic_plan_cartesian_par), pointer :: plan
  sll_real64                                   :: average_err
  sll_int32, dimension(1:2)                    :: global
  sll_int32                                    :: gi, gj
  sll_int32                                    :: prank
  type(layout_2D), pointer                     :: layout_x
  sll_int64                                    :: psize ! collective size
  sll_int32                                    :: nprocx, nprocy
  sll_int32                                    :: e
  sll_real32                                   :: ok 
  sll_real32, dimension(1)                     :: prod4test
  sll_real64                                   :: dt

  ok = 1.0

  !Boot parallel environment
  call sll_boot_collective()

  ! Number of cells is equal to number of points in this case
  ncx = 512
  ncy = 512
  Lx  = 2.0*sll_pi
  Ly  = 2.0*sll_pi

  psize = sll_get_collective_size(sll_world_collective)
  prank = sll_get_collective_rank(sll_world_collective)

  dx = Lx/ncx
  dy = Ly/ncy

  psize  = sll_get_collective_size(sll_world_collective)
  e      = int(log(real(psize))/log(2.))
  print *, 'running on ', 2**e, 'processes'

  ! Layout and local sizes for FFTs in x-direction
  layout_x => new_layout_2D( sll_world_collective )
  nprocx = 1
  nprocy = 2**e
  call initialize_layout_with_distributed_2D_array( ncx, ncy, &
       nprocx, nprocy, layout_x )

  if (prank == MPI_MASTER) print *, 'Printing 2D layout: '
  call sll_view_lims_2D( layout_x )
  call flush(6)
  if (prank == MPI_MASTER) print *, '--------------------'
  call sll_halt_collective()
  stop

  plan => new_maxwell_2d_periodic_plan_cartesian_par(layout_x, ncx, ncy, Lx, Ly)

  call compute_local_sizes( layout_x, nx_loc, ny_loc )

  SLL_ALLOCATE(bz_an(nx_loc,ny_loc), ierr)
  SLL_ALLOCATE(ex(nx_loc,ny_loc), ierr)
  SLL_ALLOCATE(ey(nx_loc,ny_loc), ierr)
  SLL_ALLOCATE(bz(nx_loc,ny_loc), ierr)
  stop

  ! initialize reference array
  do j=1,ny_loc
     do i=1,nx_loc
        global = local_to_global_2D( layout_x, (/i, j/))
        gi = global(1)
        gj = global(2)
        x  = (gi-1)*dx
        y  = (gj-1)*dy
        bz_an(i,j) = cos(x)*cos(y) 
     end do
  end do

  bz = bz_an

  call parallel_hdf5_write_array_2d( 'bz_init.h5', ncx, ncy, bz, 'bz', layout_x)

  call solve_maxwell_2d_periodic_cartesian_par(plan, dt, ex, ey, bz, faraday)
  call solve_maxwell_2d_periodic_cartesian_par(plan, dt, ex, ey, bz, ampere)


  call parallel_hdf5_write_array_2d( 'bz_analytical.h5', ncx, ncy, bz_an, 'bz_an', layout_x)
  call parallel_hdf5_write_array_2d( 'bz_computed.h5', ncx, ncy, bz, 'bz', layout_x)
  average_err  = 0.d0
  do j=1,ny_loc
     do i=1,nx_loc
        average_err  = average_err + abs(bz_an(i,j) - bz(i,j))
     enddo
  enddo
     
  average_err  = average_err/(ncx*ncy)

  call flush(6); print*, ' ------------------'
  call flush(6); print*, ' prank ', prank
  call flush(6); print*, 'local average error:', average_err
  call flush(6); print*, 'dx*dy =', dx*dy
  call flush(6); print*, ' ------------------'

  if (average_err> 1.0e-15 ) then
     print*, 'Test stopped by "sll_maxwell_2d_periodic_par" failure'
     stop
  endif
! 
  SLL_DEALLOCATE_ARRAY(ex, ierr)
  SLL_DEALLOCATE_ARRAY(ey, ierr)
  SLL_DEALLOCATE_ARRAY(bz, ierr)
  SLL_DEALLOCATE_ARRAY(bz_an, ierr)

  call sll_collective_reduce_real(sll_world_collective, (/ ok /), &
       1, MPI_PROD, 0, prod4test )
  if (prank==MPI_MASTER) then

     if (prod4test(1)==1.) then
        call flush(6)
        print*, ' '
        call flush(6)
        print*, '"sll_maxwell_2d_periodic_cart_par" test: PASS'
        call flush(6)
        print*, ' '
     endif
  endif

  call delete_maxwell_2d_periodic_plan_cartesian_par(plan)


  call sll_halt_collective()

contains

  ! Experimental interface for an hdf5 array writer in parallel
  subroutine parallel_hdf5_write_array_2d( &
    filename, &
    n_pts1, &
    n_pts2, &
    array, &
    dataset_name, &
    layout )

    character(len=*), intent(in)           :: filename
    sll_int32, intent(in)                  :: n_pts1
    sll_int32, intent(in)                  :: n_pts2
    integer(HSIZE_T), dimension(1:2)       :: global_dims
    sll_real64, dimension(:,:), intent(in) :: array
    character(len=*), intent(in)           :: dataset_name
    type(layout_2D), pointer               :: layout
    sll_int32                              :: error
    sll_int32                              :: file_id
    integer(HSIZE_T), dimension(1:2)       :: offset
    type(sll_collective_t), pointer        :: col
    sll_int32                              :: prank

    SLL_ASSERT( associated(layout) )
    col => get_layout_collective( layout )
    prank = sll_get_collective_rank( col )
    global_dims(:) = (/ n_pts1,n_pts2 /)

    offset(1) = get_layout_i_min( layout, prank ) - 1
    offset(2) = get_layout_j_min( layout, prank ) - 1

    call sll_hdf5_file_create(trim(filename),file_id,error)
    call sll_hdf5_write_array(file_id,global_dims,offset, &
                              dble(array),trim(dataset_name),error)
    call sll_hdf5_file_close(file_id,error)
  end subroutine parallel_hdf5_write_array_2d


end program test_maxwell_2d_periodic_cart_par
