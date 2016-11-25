program test_poisson_2d_periodic_cart_par
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use iso_fortran_env, only: &
    output_unit

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_o_collective_reduce, &
    sll_t_collective_t, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_gnuplot_parallel, only: &
    sll_s_gnuplot_rect_2d_parallel

  use sll_m_hdf5_io_parallel, only: &
    sll_t_hdf5_par_handle,      &
    sll_s_hdf5_par_file_create, &
    sll_o_hdf5_par_write_array, &
    sll_s_hdf5_par_file_close

  use sll_m_poisson_2d_periodic_cartesian_par, only: &
    sll_s_delete_poisson_2d_periodic_plan_cartesian_par, &
    sll_f_new_poisson_2d_periodic_plan_cartesian_par, &
    sll_f_new_poisson_2d_periodic_plan_cartesian_par_alt, &
    sll_t_poisson_2d_periodic_plan_cartesian_par, &
    sll_s_solve_poisson_2d_periodic_cartesian_par, &
    sll_s_solve_poisson_2d_periodic_cartesian_par_alt

  use sll_m_remapper, only: &
    sll_o_compute_local_sizes, &
    sll_o_get_layout_collective, &
    sll_o_get_layout_i_min, &
    sll_o_get_layout_j_min, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_o_view_lims

  use sll_mpi, only: &
    mpi_prod

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type (sll_t_poisson_2d_periodic_plan_cartesian_par), pointer :: plan
  type (sll_t_poisson_2d_periodic_plan_cartesian_par), pointer :: plan_alt

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
  type(sll_t_layout_2d), pointer                :: layout_x
  type(sll_t_layout_2d), pointer                :: layout_alt
  sll_int64                               :: colsz ! collective size
  sll_int32                               :: nprocx, nprocy
  sll_int32                               :: e
  sll_real32                              :: ok 
  sll_real32, dimension(1)                :: prod4test
  sll_int32                               :: offset(2)

  ok = 1.0

  !Boot parallel environment
  call sll_s_boot_collective()

  ! Number of cells is equal to number of points in this case
  ncx = 512
  ncy = 512
  Lx  = 2.0*sll_p_pi
  Ly  = 2.0*sll_p_pi

  colsz  = int(sll_f_get_collective_size(sll_v_world_collective),i64)
  myrank = sll_f_get_collective_rank(sll_v_world_collective)

  dx = Lx/ncx
  dy = Ly/ncy

  colsz  = int(sll_f_get_collective_size(sll_v_world_collective),i64)
  e      = int(log(real(colsz))/log(2.))
  print *, 'running on ', 2**e, 'processes'

!########### ALTERNATIVE SOLVER ##########################################

  ! Layout and local sizes for FFTs in x-direction
  layout_alt => sll_f_new_layout_2d( sll_v_world_collective )
  nprocx = 1
  nprocy = 2**e
  call sll_o_initialize_layout_with_distributed_array( ncx, ncy, &
       nprocx, nprocy, layout_alt )

  plan_alt => sll_f_new_poisson_2d_periodic_plan_cartesian_par_alt(&
       layout_alt, ncx, ncy, Lx, Ly)

  call sll_o_compute_local_sizes( layout_alt, nx_loc, ny_loc )
  call sll_o_view_lims( layout_alt )

  SLL_ALLOCATE(rho(nx_loc,ny_loc), error)
  SLL_ALLOCATE(phi_an(nx_loc,ny_loc), error)
  SLL_ALLOCATE(phi(nx_loc,ny_loc), error)

  ! initialize reference array
  do j=1,ny_loc
     do i=1,nx_loc
        global = sll_o_local_to_global( layout_alt, (/i, j/))
        gi = global(1)
        gj = global(2)
        x  = (gi-1)*dx
        y  = (gj-1)*dy
        phi_an(i,j) =  cos(x)*sin(y) 
        rho(i,j)    = -2.0_f64*phi_an(i,j)
     end do
  end do

  call parallel_hdf5_write_array_2d( 'q_density.h5', &
     ncx, ncy, rho,  'rho', layout_alt)

  call sll_s_solve_poisson_2d_periodic_cartesian_par_alt(plan_alt, rho, phi)

  call parallel_hdf5_write_array_2d( 'phi_analytical.h5', &
     ncx, ncy, phi_an, 'phi_an', layout_alt)
  call parallel_hdf5_write_array_2d( 'phi_computed.h5',   &
     ncx, ncy, phi, 'phi', layout_alt)

  average_err  = sum(abs(phi_an-phi))/real(ncx*ncy,f64)

  flush( output_unit ); print*, ' ------------------'
  flush( output_unit ); print*, ' myrank ', myrank
  flush( output_unit ); print*, 'local average error:', average_err
  flush( output_unit ); print*, 'dx*dy =', dx*dy
  flush( output_unit ); print*, ' ------------------'

  if (average_err> 1.0e-06 ) then
     print*, 'Test stopped by "sll_poisson_2d_periodic_par" failure'
     !call sll_s_halt_collective()
     !stop
  endif
 
  SLL_DEALLOCATE_ARRAY(phi,    error)
  SLL_DEALLOCATE_ARRAY(rho,    error)
  SLL_DEALLOCATE_ARRAY(phi_an, error)

  call sll_s_delete_poisson_2d_periodic_plan_cartesian_par(plan_alt)

!#########FIRST VERSION WITH LAST PERIODIC POINT ADDED #########################

  ! Layout and local sizes for FFTs in x-direction
  layout_x => sll_f_new_layout_2d( sll_v_world_collective )
  nprocx = 1
  nprocy = 2**e
  call sll_o_initialize_layout_with_distributed_array( ncx+1, ncy+1, &
       nprocx, nprocy, layout_x )

  plan => sll_f_new_poisson_2d_periodic_plan_cartesian_par(&
       layout_x, ncx, ncy, Lx, Ly)

  call sll_o_compute_local_sizes( layout_x, nx_loc, ny_loc )
  call sll_o_view_lims( layout_x )

  SLL_ALLOCATE(rho(nx_loc,ny_loc), error)
  SLL_ALLOCATE(phi_an(nx_loc,ny_loc), error)
  SLL_ALLOCATE(phi(nx_loc,ny_loc), error)

  ! initialize reference array
  do j=1,ny_loc
     do i=1,nx_loc
        global = sll_o_local_to_global( layout_x, (/i, j/))
        gi = global(1)
        gj = global(2)
        x  = (gi-1)*dx
        y  = (gj-1)*dy
        phi_an(i,j) =  cos(x)*sin(y) 
        rho(i,j)    = -2.0_f64*phi_an(i,j)
     end do
  end do

  call sll_s_solve_poisson_2d_periodic_cartesian_par(plan, rho, phi)

  offset(1) =  sll_o_get_layout_i_min( layout_x, myrank ) - 1
  offset(2) =  sll_o_get_layout_j_min( layout_x, myrank ) - 1
  call sll_s_gnuplot_rect_2d_parallel(dble(offset(1)), dble(1), &
                                    dble(offset(2)), dble(1), &
                                    size(rho,1), size(rho,2), &
                                    rho, "rho", 1, error)  

  average_err  = sum(abs(phi_an-phi))/real(ncx*ncy,f64)

  flush( output_unit ); print*, ' ------------------'
  flush( output_unit ); print*, ' myrank ', myrank
  flush( output_unit ); print*, 'local average error:', average_err
  flush( output_unit ); print*, 'dx*dy =', dx*dy
  flush( output_unit ); print*, ' ------------------'

  if (average_err> 1.0e-06 ) then
     print*, 'Test stopped by "sll_poisson_2d_periodic_par" failure'
     call sll_s_halt_collective()
     stop
  endif
 
  SLL_DEALLOCATE_ARRAY(phi,    error)
  SLL_DEALLOCATE_ARRAY(rho,    error)
  SLL_DEALLOCATE_ARRAY(phi_an, error)

  call sll_o_collective_reduce(sll_v_world_collective, (/ ok /), &
       1, MPI_PROD, 0, prod4test )

  if (myrank==0) then
     if (prod4test(1)==1.) then
        flush( output_unit )
        print*, ' '
        flush( output_unit )
        print*, '"sll_poisson_2d_periodic_cart_par" test: PASSED'
        flush( output_unit )
        print*, ' '
     endif
  endif

  call sll_s_delete_poisson_2d_periodic_plan_cartesian_par(plan)

  call sll_s_halt_collective()

contains

  ! Experimental interface for an hdf5 array writer in parallel
  subroutine parallel_hdf5_write_array_2d( &
    filename, &
    n_pts1, &
    n_pts2, &
    array, &
    dataset_name, &
    layout )

    character(len=*)     , intent(in) :: filename
    sll_int32            , intent(in) :: n_pts1
    sll_int32            , intent(in) :: n_pts2
    sll_real64           , intent(in) :: array(:,:)
    character(len=*)     , intent(in) :: dataset_name
    type(sll_t_layout_2d), pointer    :: layout

    type(sll_t_collective_t), pointer :: col
    type(sll_t_hdf5_par_handle)           :: handle
    integer(i64)                      :: global_dims(1:2)
    integer(i64)                      :: offset(1:2)
    sll_int32                         :: error
    sll_int32                         :: myrank
    sll_int32                         :: comm

    SLL_ASSERT( associated(layout) )
    col => sll_o_get_layout_collective( layout )
    myrank = sll_f_get_collective_rank( col )
    global_dims(:) = int( [n_pts1,n_pts2], i64 )
    
    offset(1) = int( sll_o_get_layout_i_min( layout, myrank )-1, i64 )
    offset(2) = int( sll_o_get_layout_j_min( layout, myrank )-1, i64 )

    comm = sll_v_world_collective%comm
    call sll_s_hdf5_par_file_create( filename, comm, handle, error )
    call sll_o_hdf5_par_write_array( handle, global_dims, offset, &
                              array, dataset_name, error )
    call sll_s_hdf5_par_file_close( handle, error )

  end subroutine parallel_hdf5_write_array_2d

end program test_poisson_2d_periodic_cart_par
