program test_poisson_2d_dirichlet_cartesian
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_constants
  use sll_poisson_2d_dirichlet_cartesian
  use sll_gnuplot
  implicit none

  type (poisson_2d_dirichlet_cartesian), pointer :: plan

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
  sll_int64                               :: colsz ! collective size
  sll_int32                               :: nprocx, nprocy
  sll_int32                               :: e
  sll_real32                              :: ok 
  sll_real32, dimension(1)                :: prod4test
  sll_int32                               :: offset(2)


  ! Number of cells is equal to number of points in this case
  ncx = 512
  ncy = 512
  Lx  = 2.0*sll_pi
  Ly  = 2.0*sll_pi

  dx = Lx/ncx
  dy = Ly/ncy

  plan => new_poisson_2d_dirichlet_cartesian_plan(ncx, ncy, Lx, Ly)

  SLL_ALLOCATE(rho(nx_loc,ny_loc), error)
  SLL_ALLOCATE(phi_an(nx_loc,ny_loc), error)
  SLL_ALLOCATE(phi(nx_loc,ny_loc), error)

  ! initialize reference array
  do j=1,ny_loc
     do i=1,nx_loc
        x  = (i-1)*dx
        y  = (j-1)*dy
        phi_an(i,j) =  cos(x)*sin(y) 
        rho(i,j)    = -2.0_f64*phi_an(i,j)
     end do
  end do

  call solve_poisson_2d_dirichlet_cartesian(plan, rho, phi)

  average_err  = sum(abs(phi_an-phi))/(ncx*ncy)

  call flush(6); print*, ' ------------------'
  call flush(6); print*, ' myrank ', myrank
  call flush(6); print*, 'local average error:', average_err
  call flush(6); print*, 'dx*dy =', dx*dy
  call flush(6); print*, ' ------------------'

  if (average_err> 1.0e-06 ) then
     print*, 'Test stopped by "sll_poisson_2d_periodic_par" failure'
     !call sll_halt_collective()
     !stop
  endif
 
  SLL_DEALLOCATE_ARRAY(phi,    error)
  SLL_DEALLOCATE_ARRAY(rho,    error)
  SLL_DEALLOCATE_ARRAY(phi_an, error)

  call sll_delete(plan)

!!$  call sll_gnuplot_rect_2d(dble(offset(1)), dble(1), &
!!$       dble(offset(2)), dble(1), &
!!$       size(rho,1), size(rho,2), &
!!$       rho, "rho", 1, error)  

  average_err  = sum(abs(phi_an-phi))/(ncx*ncy)

  call flush(6); print*, ' ------------------'
  call flush(6); print*, ' myrank ', myrank
  call flush(6); print*, 'local average error:', average_err
  call flush(6); print*, 'dx*dy =', dx*dy
  call flush(6); print*, ' ------------------'

  print *, 'PASSED'
end program test_poisson_2d_dirichlet_cartesian
