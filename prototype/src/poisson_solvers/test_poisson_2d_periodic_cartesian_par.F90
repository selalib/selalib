program test_poisson_2d_periodic_cart_par
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_remapper
  use sll_constants
  use sll_poisson_2d_periodic_cartesian_par
  use sll_collective
  use hdf5
  use sll_hdf5_io_parallel
  use sll_gnuplot_parallel

  implicit none

  type (poisson_2d_periodic_plan_cartesian_par), pointer :: plan
  type (poisson_2d_periodic_plan_cartesian_par), pointer :: plan_alt

  sll_int32                               :: ncx, ncy
  sll_int32                               :: nx_loc, ny_loc
  sll_int32                               :: error
  sll_real64                              :: Lx, Ly
  sll_real64                              :: dx, dy
  sll_real64                              :: x, y
  sll_real64, dimension(:,:), allocatable :: rho
  sll_real64, dimension(:,:), allocatable :: phi_an
  sll_real64, dimension(:,:), allocatable :: phi
  sll_int32                               :: i, j
  sll_real64                              :: average_err
  sll_int32, dimension(1:2)               :: global
  sll_int32                               :: gi, gj
  sll_int32                               :: myrank
  type(layout_2D), pointer                :: layout_x
  sll_int64                               :: colsz ! collective size
  sll_int32                               :: nprocx, nprocy
  sll_int32                               :: e
  sll_real32                              :: ok 
  sll_real32, dimension(1)                :: prod4test
  sll_int32                               :: offset(2)

  ok = 1.0

  !Boot parallel environment
  call sll_boot_collective()

  ! Number of cells is equal to number of points in this case
  ncx = 512
  ncy = 512
  Lx  = 2.0*sll_pi
  Ly  = 2.0*sll_pi

  colsz  = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  dx = Lx/ncx
  dy = Ly/ncy

  colsz  = sll_get_collective_size(sll_world_collective)
  e      = int(log(real(colsz))/log(2.))
  print *, 'running on ', 2**e, 'processes'

!########### ALTERNATIVE SOLVER ##########################################

  ! Layout and local sizes for FFTs in x-direction
  layout_x => new_layout_2D( sll_world_collective )
  nprocx = 1
  nprocy = 2**e
  call initialize_layout_with_distributed_2D_array( ncx, ncy, &
       nprocx, nprocy, layout_x )

  plan_alt => new_poisson_2d_periodic_plan_cartesian_par_alt(&
       layout_x, ncx, ncy, Lx, Ly)

  call compute_local_sizes( layout_x, nx_loc, ny_loc )
  call sll_view_lims_2D( layout_x )

  SLL_ALLOCATE(rho(nx_loc,ny_loc), error)
  SLL_ALLOCATE(phi_an(nx_loc,ny_loc), error)
  SLL_ALLOCATE(phi(nx_loc,ny_loc), error)

  ! initialize reference array
  do j=1,ny_loc
     do i=1,nx_loc
        global = local_to_global_2D( layout_x, (/i, j/))
        gi = global(1)
        gj = global(2)
        x  = (gi-1)*dx
        y  = (gj-1)*dy
        phi_an(i,j) =  cos(x)*sin(y) 
        rho(i,j)    = -2.0_f64*phi_an(i,j)
     end do
  end do

  call parallel_hdf5_write_array_2d( 'q_density.h5', &
     ncx, ncy, rho,  'rho', layout_x)

  call solve_poisson_2d_periodic_cartesian_par(plan_alt, rho, phi)

  call parallel_hdf5_write_array_2d( 'phi_analytical.h5', &
     ncx, ncy, phi_an, 'phi_an', layout_x)
  call parallel_hdf5_write_array_2d( 'phi_computed.h5',   &
     ncx, ncy, phi, 'phi', layout_x)

  average_err  = sum(abs(phi_an-phi))/(ncx*ncy)

  call flush(6); print*, ' ------------------'
  call flush(6); print*, ' myrank ', myrank
  call flush(6); print*, 'local average error:', average_err
  call flush(6); print*, 'dx*dy =', dx*dy
  call flush(6); print*, ' ------------------'

  if (average_err> 1.0e-06 ) then
     print*, 'Test stopped by "sll_poisson_2d_periodic_par" failure'
     call sll_halt_collective()
     stop
  endif
 
  SLL_DEALLOCATE_ARRAY(phi,    error)
  SLL_DEALLOCATE_ARRAY(rho,    error)
  SLL_DEALLOCATE_ARRAY(phi_an, error)

  call delete_poisson_2d_periodic_plan_cartesian_par(plan_alt)
  call delete_layout_2D( layout_x )

!#########FIRST VERSION WITH LAST PERIODIC POINT ADDED #########################

  ! Layout and local sizes for FFTs in x-direction
  layout_x => new_layout_2D( sll_world_collective )
  nprocx = 1
  nprocy = 2**e
  call initialize_layout_with_distributed_2D_array( ncx+1, ncy+1, &
       nprocx, nprocy, layout_x )

  plan => new_poisson_2d_periodic_plan_cartesian_par(&
       layout_x, ncx, ncy, Lx, Ly)

  call compute_local_sizes( layout_x, nx_loc, ny_loc )
  call sll_view_lims_2D( layout_x )

  SLL_ALLOCATE(rho(nx_loc,ny_loc), error)
  SLL_ALLOCATE(phi_an(nx_loc,ny_loc), error)
  SLL_ALLOCATE(phi(nx_loc,ny_loc), error)

  ! initialize reference array
  do j=1,ny_loc
     do i=1,nx_loc
        global = local_to_global_2D( layout_x, (/i, j/))
        gi = global(1)
        gj = global(2)
        x  = (gi-1)*dx
        y  = (gj-1)*dy
        phi_an(i,j) =  cos(x)*sin(y) 
        rho(i,j)    = -2.0_f64*phi_an(i,j)
     end do
  end do

  call parallel_hdf5_write_array_2d( 'q_density.h5', &
     ncx+1, ncy+1, rho,  'rho', layout_x)

  call solve_poisson_2d_periodic_cartesian_par(plan, rho, phi)

  offset(1) =  get_layout_2D_i_min( layout_x, myrank ) - 1
  offset(2) =  get_layout_2D_j_min( layout_x, myrank ) - 1
  call sll_gnuplot_rect_2d_parallel(dble(offset(1)), dble(1), &
                                    dble(offset(2)), dble(1), &
                                    rho, "rho", 1, error)  

  call parallel_hdf5_write_array_2d( 'phi_analytical.h5', &
     ncx+1, ncy+1, phi_an, 'phi_an', layout_x)
!  call parallel_hdf5_write_array_2d( 'phi_computed.h5',   &
!     ncx+1, ncy+1, phi, 'phi', layout_x)

  average_err  = 0.d0
  do j=1,ny_loc
     do i=1,nx_loc
        average_err  = average_err + abs(phi_an(i,j) - phi(i,j))
     enddo
  enddo
     
  average_err  = average_err/(ncx*ncy)

  call flush(6); print*, ' ------------------'
  call flush(6); print*, ' myrank ', myrank
  call flush(6); print*, 'local average error:', average_err
  call flush(6); print*, 'dx*dy =', dx*dy
  call flush(6); print*, ' ------------------'

  if (average_err> 1.0e-06 ) then
     print*, 'Test stopped by "sll_poisson_2d_periodic_par" failure'
     call sll_halt_collective()
     stop
  endif
 
  SLL_DEALLOCATE_ARRAY(phi,    error)
  SLL_DEALLOCATE_ARRAY(rho,    error)
  SLL_DEALLOCATE_ARRAY(phi_an, error)

  call sll_collective_reduce(sll_world_collective, (/ ok /), &
       1, MPI_PROD, 0, prod4test )

  if (myrank==0) then
     if (prod4test(1)==1.) then
        call flush(6)
        print*, ' '
        call flush(6)
        print*, '"sll_poisson_2d_periodic_cart_par" test: PASSED'
        call flush(6)
        print*, ' '
     endif
  endif

  call delete_poisson_2d_periodic_plan_cartesian_par(plan)

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
    integer(HID_T)                         :: file_id
    integer(HSIZE_T), dimension(1:2)       :: offset
    type(sll_collective_t), pointer        :: col
    sll_int32                              :: myrank

    SLL_ASSERT( associated(layout) )
    col => get_layout_collective( layout )
    myrank = sll_get_collective_rank( col )
    global_dims(:) = (/ n_pts1,n_pts2 /)
    
    offset(1) = get_layout_i_min( layout, myrank ) - 1
    offset(2) = get_layout_j_min( layout, myrank ) - 1

    call sll_hdf5_file_create(filename,file_id,error)
    call sll_hdf5_write_array(file_id,global_dims,offset, &
                              array,dataset_name,error)
    call sll_hdf5_file_close(file_id,error)

  end subroutine parallel_hdf5_write_array_2d

end program test_poisson_2d_periodic_cart_par
