program test_poisson_2d_periodic_cart_par
#include "sll_remap.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use numeric_constants
  use sll_poisson_2d_periodic_cartesian_par
  use sll_collective
  implicit none

  sll_int32                                    :: ncx, ncy
  sll_int32                                    :: nx_loc, ny_loc
  sll_int32                                    :: ierr
  sll_real64                                   :: Lx, Ly
  sll_real64                                   :: dx, dy
  sll_real64, dimension(:,:), allocatable      :: x, y
  sll_real64, dimension(:,:), allocatable      :: rho
  sll_real64, dimension(:,:), allocatable      :: phi_an
  sll_real64, dimension(:,:), allocatable      :: phi
  sll_int32                                    :: i, j
  type (poisson_2d_periodic_plan_cartesian_par), pointer :: plan
  sll_real64                                   :: average_err
  sll_int32, dimension(1:2)                    :: global
  sll_int32                                    :: gi, gj
  sll_int32                                    :: myrank
  type(layout_2D), pointer                     :: layout_x
  sll_int64                                    :: colsz ! collective size
  sll_int32                                    :: i_test
  sll_int32                                    :: nprocx, nprocy
  sll_int32                                    :: e
  sll_real32                                   :: ok 
  sll_real32, dimension(1)                     :: prod4test

  ok = 1.0

  !Boot parallel environment
  call sll_boot_collective()

  ! Number of cells is equal to number of points in this case
  ncx = 64
  ncy = 64
  Lx  = 2*sll_pi
  Ly  = 2*sll_pi

  colsz  = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  dx = Lx/ncx
  dy = Ly/ncy

  colsz  = sll_get_collective_size(sll_world_collective)
  e      = int(log(real(colsz))/log(2.))

  ! Layout and local sizes for FFTs in x-direction
  layout_x => new_layout_2D( sll_world_collective )
  nprocx = 1
  nprocy = 2**e
  call initialize_layout_with_distributed_2D_array( ncx, ncy, &
       nprocx, nprocy, layout_x )

  plan => new_poisson_2d_periodic_plan_cartesian_par(layout_x, ncx, ncy, Lx, Ly)

  call compute_local_sizes( layout_x, nx_loc, ny_loc )

  SLL_ALLOCATE(rho(nx_loc,ny_loc), ierr)
  SLL_ALLOCATE(x(nx_loc,ny_loc), ierr)
  SLL_ALLOCATE(y(nx_loc,ny_loc), ierr)
  SLL_ALLOCATE(phi_an(nx_loc,ny_loc), ierr)

  do j=1,ny_loc
     do i=1,nx_loc
        global = local_to_global_2D( layout_x, (/i, j/))
        gi = global(1)
        gj = global(2)
        x(i,j) = (gi-1)*dx
        y(i,j) = (gj-1)*dy
     end do
  end do
  
  do i_test=1, 2
     if (i_test==1) then
        phi_an = cos(x)*sin(y)
     else if (i_test == 2) then
        phi_an = (4/(sll_pi * sqrt(sll_pi)*Lx*Ly)) &
             * exp(-.5*(x-Lx/2)**2)                   &
             * exp(-.5*(y-Ly/2)**2)
     end if

     do j=1,ny_loc
        do i=1,nx_loc
           if (i_test == 1) then
              rho(i,j) = 3*phi_an(i,j)
           else if(i_test == 2) then
              rho(i,j) = phi_an(i,j)*(3-((x(i,j)-Lx/2)**2 + &
                                         (y(i,j)-Ly/2)**2))
           end if
        enddo
     enddo

     SLL_ALLOCATE(phi(nx_loc,ny_loc), ierr)
     call solve_poisson_2d_periodic_cartesian_par(plan, rho, phi)

     average_err  = 0.d0
     do j=1,ny_loc
        do i=1,nx_loc
           average_err  = average_err + abs(phi_an(i,j) - phi(i,j))
        enddo
     enddo
     
     average_err  = average_err/(nx_loc*ny_loc)

     call flush(); print*, ' ------------------'
     call flush(); print*, ' myrank ', myrank
     call flush(); print*, 'local average error:', average_err
     call flush(); print*, 'dx*dy =', dx*dy
     call flush(); print*, ' ------------------'

     if (average_err> dx*dy ) then
        print*, 'Test stopped by "sll_poisson_3d_periodic_par" failure'
        stop
     endif
  end do
  SLL_DEALLOCATE_ARRAY(phi, ierr)
  SLL_DEALLOCATE_ARRAY(x, ierr)
  SLL_DEALLOCATE_ARRAY(y, ierr)
  SLL_DEALLOCATE_ARRAY(rho, ierr)
  SLL_DEALLOCATE_ARRAY(phi_an, ierr)

  call sll_collective_reduce_real(sll_world_collective, (/ ok /), &
       1, MPI_PROD, 0, prod4test )
  if (myrank==0) then

     if (prod4test(1)==1.) then
        call flush()
        print*, ' '
        call flush()
        print*, '"sll_poisson_3d_periodic_par" test: PASS'
        call flush()
        print*, ' '
     endif
  endif

  call delete_poisson_2d_periodic_plan_cartesian_par(plan)

  SLL_DEALLOCATE_ARRAY(phi, ierr)
  SLL_DEALLOCATE_ARRAY(x, ierr)
  SLL_DEALLOCATE_ARRAY(y, ierr)
  SLL_DEALLOCATE_ARRAY(rho, ierr)
  SLL_DEALLOCATE_ARRAY(phi_an, ierr)

  call sll_halt_collective()

end program test_poisson_2d_periodic_cart_par
