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

  sll_int32   :: ncx
  sll_int32   :: ncy
  sll_int32   :: nx_loc
  sll_int32   :: ny_loc
  sll_int32   :: error
  sll_real64  :: Lx, Ly
  sll_real64  :: dx, dy
  sll_real64  :: x, y
  sll_int32   :: i, j
  sll_int32   :: gi, gj
  sll_int32   :: prank
  sll_int64   :: psize 
  sll_int32   :: nprocx
  sll_int32   :: nprocy
  sll_int32   :: e
  sll_real32  :: ok 
  sll_real64  :: dt
  sll_int32   :: istep
  sll_int32   :: nstep
  sll_real64  :: time = 0.0
  sll_real64  :: omega
  sll_real64  :: err_l2
  sll_int32   :: mode = 1
  sll_real64  :: tcpu1
  sll_real64  :: tcpu2

  type (maxwell_2d_periodic_plan_cartesian_par), pointer :: plan

  type(layout_2D), pointer                     :: layout_x
  type(layout_2D), pointer                     :: layout_y
  sll_real64, dimension(:,:), allocatable      :: ex
  sll_real64, dimension(:,:), allocatable      :: ey
  sll_int32, dimension(1:2)                    :: global
  character(len=4)                             :: cstep

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

  omega  = sqrt((mode*sll_pi/Lx)**2+(mode*sll_pi/Ly)**2)

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
                                                     ncx, ncy, &
                                                     Lx, Ly)
  plan%e_0  = 1.0_f64 
  plan%mu_0 = 1.0_f64 

  !Ex si sequential along y
  call compute_local_sizes_2d(plan%layout_y,nx_loc,ny_loc)
  SLL_CLEAR_ALLOCATE(ex(nx_loc,ny_loc), error)

  !Ey si sequential along x
  call compute_local_sizes_2d(plan%layout_x,nx_loc,ny_loc)
  SLL_CLEAR_ALLOCATE(ey(nx_loc,ny_loc), error)

  do j=1,ny_loc
     do i=1,nx_loc
        global = local_to_global_2D( layout_x, (/i, j/))
        gi = global(1); gj = global(2)
        x  = (gi-1)*dx; y  = (gj-1)*dy
        plan%bz_x(i,j) = - cos(mode*(i-1)*dx) &
                         * cos(mode*(j-1)*dy) &
                         * cos(omega*time)
     end do
  end do

  tcpu1 = MPI_WTIME()

  nstep = 100
  dt = 0.5 / sqrt (1./(dx*dx)+1./(dy*dy))
  print*, ' dt = ', dt
  call faraday_te(plan,0.5*dt,ex,ey)
  time = 0.5*dt

  do istep = 1, nstep

     call int2string(istep,cstep)
     call ampere_te(plan,dt,ex,ey)
     time = time + dt
     call faraday_te(plan,dt,ex,ey)
     call write_fields('fields-'//cstep)

     err_l2 = 0.0
     do j=1,ny_loc
     do i=1,nx_loc
        global = local_to_global_2D( layout_x, (/i, j/))
        gi = global(1); gj = global(2)
        x  = (gi-1)*dx; y  = (gj-1)*dy
        err_l2 = err_l2 + abs(plan%bz_x(i,j) + cos(mode*(i-1)*dx) &
                 * cos(mode*(j-1)*dy) * cos(omega*time))
     end do
     end do

     if ( prank == MPI_MASTER ) then
        write(*,"(10x,' istep = ',I6)",advance="no") istep
        write(*,"(' time = ',g12.3,' sec')",advance="no") time
        write(*,"(' L2 error = ',g15.5)") err_l2 / (ncx*ncy)
     end if

  end do
     
  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) then
     write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize
  end if

  SLL_DEALLOCATE_ARRAY(ex, error)
  SLL_DEALLOCATE_ARRAY(ey, error)

  call delete_maxwell_2d_periodic_plan_cartesian_par(plan)

  call sll_halt_collective()

contains

  !> Experimental interface for an hdf5 array writer in parallel
  subroutine write_fields( file_name )

    character(len=*), intent(in)     :: file_name
    integer(HSIZE_T), dimension(1:2) :: global_dims
    sll_int32                        :: error
    sll_int32                        :: file_id
    integer(HSIZE_T), dimension(1:2) :: offset

    global_dims(:) = (/ ncx,ncy /)

    call sll_hdf5_file_create(file_name//'.h5',file_id,error)
    offset(1) = get_layout_i_min(layout_y,prank) - 1
    offset(2) = get_layout_j_min(layout_y,prank) - 1
    call sll_hdf5_write_array(file_id,global_dims,offset, &
                              ex,'ex_values',error)
    call sll_hdf5_write_array(file_id,global_dims,offset, &
                              plan%bz_y,'bz_y_values',error)
    offset(1) = get_layout_i_min(layout_x,prank) - 1
    offset(2) = get_layout_j_min(layout_x,prank) - 1
    call sll_hdf5_write_array(file_id,global_dims,offset, &
                              ey,'ey_values',error)
    call sll_hdf5_write_array(file_id,global_dims,offset, &
                              plan%bz_x,'bz_x_values',error)
    call sll_hdf5_file_close(file_id,error)

  end subroutine write_fields

end program test_maxwell_2d_periodic_cart_par
