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

  sll_int32                                    :: ncx
  sll_int32                                    :: ncy
  sll_int32                                    :: nx_loc
  sll_int32                                    :: ny_loc
  sll_int32                                    :: error
  sll_real64                                   :: Lx, Ly
  sll_real64                                   :: dx, dy
  sll_real64                                   :: x, y
  sll_real64, dimension(:,:), allocatable      :: ex
  sll_real64, dimension(:,:), allocatable      :: ey
  sll_int32                                    :: i, j

  type (maxwell_2d_periodic_plan_cartesian_par), pointer :: plan

  type(layout_2D), pointer                     :: layout_x
  type(layout_2D), pointer                     :: layout_y

  sll_int32, dimension(1:2)                    :: global
  sll_int32                                    :: gi, gj
  sll_int32                                    :: prank
  sll_int64                                    :: psize ! collective size
  sll_int32                                    :: nprocx, nprocy
  sll_int32                                    :: e
  sll_real32                                   :: ok 
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
  layout_y => new_layout_2D( sll_world_collective )
  nprocx = 2**e
  nprocy = 2**e
  call initialize_layout_with_distributed_2D_array( ncx, &
                                                    ncy, &
                                                      1, &
                                                 nprocy, &
                                                 layout_x )

  call initialize_layout_with_distributed_2D_array( ncx, &
                                                    ncy, &
                                                 nprocx, &
                                                      1, &
                                                 layout_y )
  call flush(6)
  call sll_view_lims_2D( layout_x )
  call flush(6)
  call sll_view_lims_2D( layout_y )
  call flush(6)

  plan => new_maxwell_2d_periodic_plan_cartesian_par(layout_x, &
                                                     layout_y, &
                                                     ncx, ncy, Lx, Ly)


  call compute_local_sizes_2d(plan%layout_x,nx_loc,ny_loc)
  SLL_ALLOCATE(ey(nx_loc,ny_loc), error)
  do j=1,ny_loc
     do i=1,nx_loc
        global = local_to_global_2D( layout_x, (/i, j/))
        gi = global(1)
        gj = global(2)
        x  = (gi-1)*dx
        y  = (gj-1)*dy
        ey(i,j) = cos(x)*cos(y) 
     end do
  end do

  call parallel_hdf5_write_array_2d('ey',ncx,ncy,ey,'ey',layout_x)

  call compute_local_sizes_2d(plan%layout_y,nx_loc,ny_loc)
  SLL_ALLOCATE(ex(nx_loc,ny_loc), error)
  do j=1,ny_loc
     do i=1,nx_loc
        global = local_to_global_2D( layout_y, (/i, j/))
        gi = global(1)
        gj = global(2)
        x  = (gi-1)*dx
        y  = (gj-1)*dy
        ex(i,j) = sin(x)*sin(y) 
     end do
  end do

  call parallel_hdf5_write_array_2d('ex',ncx,ncy,ex,'ex',layout_y)

  call solve_maxwell_2d_periodic_cartesian_par(plan, dt, ex, ey, faraday)
  call solve_maxwell_2d_periodic_cartesian_par(plan, dt, ex, ey, ampere)
     
  SLL_DEALLOCATE_ARRAY(ex, error)
  SLL_DEALLOCATE_ARRAY(ey, error)

  call delete_maxwell_2d_periodic_plan_cartesian_par(plan)

  call sll_halt_collective()

contains

  !> Experimental interface for an hdf5 array writer in parallel
  subroutine parallel_hdf5_write_array_2d( &
    array_name,                            &
    n_pts1,                                &
    n_pts2,                                &
    array,                                 &
    dataset_name,                          &
    layout )

    character(len=*), intent(in)           :: array_name
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

    call sll_hdf5_file_create(array_name//'.h5',file_id,error)
    call sll_hdf5_write_array(file_id,global_dims,offset, &
                              dble(array),trim(dataset_name),error)
    call sll_hdf5_file_close(file_id,error)

  end subroutine parallel_hdf5_write_array_2d


end program test_maxwell_2d_periodic_cart_par
