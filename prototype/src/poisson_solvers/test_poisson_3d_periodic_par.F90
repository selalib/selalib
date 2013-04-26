!************************************************************************
!
! Selalib 2012     
! File : test_poisson_3d_par.F90
!
!> @brief 
!> Selalib poisson solvers (3D) unit test
!> Last modification: April 10, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!> Pierre NAVARO (navaro@math.unistra.fr)
!                                  
!************************************************************************

program test_poisson_3d_periodic_par
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_poisson_solvers.h"
  use sll_remapper
  use sll_constants
  use sll_poisson_3d_periodic_par
  use sll_collective

  implicit none

  sll_int32                                    :: nx, ny, nz
  sll_int32                                    :: nx_loc, ny_loc, nz_loc
  sll_int32                                    :: ierr
  sll_real64                                   :: Lx, Ly, Lz
  sll_real64                                   :: dx, dy, dz
  sll_real64, dimension(:,:,:), allocatable    :: x, y, z
  sll_real64, dimension(:,:,:), allocatable    :: rho
  sll_real64, dimension(:,:,:), allocatable    :: phi_an
  sll_real64, dimension(:,:,:), allocatable    :: phi
  sll_int32                                    :: i, j, k
  type (poisson_3d_periodic_plan_par), pointer :: plan
  sll_real64                                   :: average_err
  sll_int32, dimension(1:3)                    :: global
  sll_int32                                    :: gi, gj, gk
  sll_int32                                    :: myrank
  type(layout_3D), pointer                     :: layout
  sll_int64                                    :: colsz ! collective size
  sll_int32                                    :: i_test
  sll_int32                                    :: npx, npy, npz
  sll_int32                                    :: e
  sll_real32                                   :: ok = 1.0
  sll_real32   , dimension(1)                  :: prod4test

  !Boot parallel environment
  call sll_boot_collective()

  nx = 64
  ny = 64
  nz = 64
  Lx = 2*sll_pi
  Ly = 2*sll_pi
  Lz = 2*sll_pi

  colsz  = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)


  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz

  colsz  = sll_get_collective_size(sll_world_collective)
  e = int(log(real(colsz))/log(2.))

  ! Layout and local sizes for FFTs in x-direction
  layout => new_layout_3D( sll_world_collective )
  npx = 1
  npy = 2**(e/2)
  npz = int(colsz)/npy
  call initialize_layout_with_distributed_3D_array( nx, ny, &
                                  nz, npx, npy, npz, layout )

  plan => new_poisson_3d_periodic_plan_par(layout, nx, ny, &
                                             nz, Lx, Ly, Lz)

  call compute_local_sizes( layout, nx_loc, ny_loc, nz_loc )

  SLL_ALLOCATE(rho(nx_loc,ny_loc,nz_loc), ierr)
  SLL_ALLOCATE(x(nx_loc,ny_loc,nz_loc),ierr)
  SLL_ALLOCATE(y(nx_loc,ny_loc,nz_loc),ierr)
  SLL_ALLOCATE(z(nx_loc,ny_loc,nz_loc),ierr)
  SLL_ALLOCATE(phi_an(nx_loc,ny_loc,nz_loc),ierr)

  do k=1,nz_loc
     do j=1,ny_loc
        do i=1,nx_loc
           global = local_to_global_3D( layout, (/i, j, k/))
           gi = global(1)
           gj = global(2)
           gk = global(3)
           x(i,j,k) = (gi-1)*dx
           y(i,j,k) = (gj-1)*dy
           z(i,j,k) = (gk-1)*dz
        end do
     end do
  end do

  do i_test = 1, 2

     if (i_test==1) then
        phi_an = cos(x)*sin(y)*cos(z)
     else if (i_test == 2) then
        phi_an = (4/(sll_pi * sqrt(sll_pi)*Lx*Ly*Lz)) &
             * exp(-.5*(x-Lx/2)**2)                   &
             * exp(-.5*(y-Ly/2)**2) * sin(z)
     end if

     do k=1,nz_loc
        do j=1,ny_loc
           do i=1,nx_loc
              if (i_test == 1) then
                 rho(i,j,k) = 3*phi_an(i,j,k)
              else if(i_test == 2) then
                 rho(i,j,k) = phi_an(i,j,k)   &
                      *(3-((x(i,j,k)-Lx/2)**2 &
                      +(y(i,j,k)-Ly/2)**2))
              end if
           enddo
        enddo
     enddo

     SLL_ALLOCATE(phi(nx_loc,ny_loc,nz_loc), ierr)
     call solve_poisson_3d_periodic_par(plan, rho, phi)

     average_err  = 0.d0

     do k=1,nz_loc
        do j=1,ny_loc
           do i=1,nx_loc
              average_err  = average_err + abs( phi_an(i,j,k) &
                             - phi(i,j,k) )
           enddo
        enddo
     enddo

     average_err  = average_err  / (nx_loc*ny_loc*nz_loc)

     call flush(6); print*, ' ------------------'
     call flush(6); print*, ' myrank ', myrank
     call flush(6); print*, 'local average error:', average_err
     call flush(6); print*, 'dx*dy*dz =', dx*dy*dz
     call flush(6); print*, ' ------------------'

     if (average_err> dx*dy*dz ) then
        print*, 'Test stopped by "sll_poisson_3d_periodic_par" failure'
        stop
     endif

     SLL_DEALLOCATE_ARRAY(phi, ierr)

  end do

     call sll_collective_reduce(sll_world_collective, (/ ok /), &
          1, MPI_PROD, 0, prod4test )
     if (myrank==0) then

        if (prod4test(1)==1.) then
           call flush(6)
           print*, ' '
           call flush(6)
           print*, '"sll_poisson_3d_periodic_par" test: PASSED'
           call flush(6)
           print*, ' '
        endif
     endif           

  call delete_poisson_3d_periodic_plan_par(plan)

  SLL_DEALLOCATE_ARRAY(phi_an, ierr)
  SLL_DEALLOCATE_ARRAY(rho, ierr)
  call sll_halt_collective()

end program test_poisson_3d_periodic_par
