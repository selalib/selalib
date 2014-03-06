program poisson_2d_parallel
#include "selalib-mpi.h"

  implicit none

  type (poisson_2d_periodic_plan_cartesian_par), pointer :: plan

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

  call solve_poisson_2d_periodic_cartesian_par(plan, rho, phi)

  offset(1) =  get_layout_2D_i_min( layout_x, myrank ) - 1
  offset(2) =  get_layout_2D_j_min( layout_x, myrank ) - 1
  call sll_gnuplot_rect_2d_parallel(dble(offset(1)), dble(1), &
                                    dble(offset(2)), dble(1), &
                                    size(phi,1), size(phi,2), &
                                    phi, "phi", 1, error)  

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

end program poisson_2d_parallel
