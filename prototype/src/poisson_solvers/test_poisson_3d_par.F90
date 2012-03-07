!***************************************************************************
!
! Selalib 2012     
! File : test_poisson_3d_par.F90
!
!> @brief 
!> Selalib poisson solvers (3D) unit test
!> Last modification: March 06, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!> Pierre NAVARO (navaro@math.unistra.fr)
!                                  
!***************************************************************************

program test_poisson_3d_par

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
#include "sll_poisson_solvers.h"
#include "sll_remap.h"

  use numeric_constants
  use sll_poisson_3d_periodic_util
  use sll_poisson_3d_periodic_par
  use sll_collective

  implicit none

  sll_int32                                 :: nx, ny, nz
  sll_int32                                 :: nx_loc, ny_loc, nz_loc
  sll_int32                                 :: ierr
  sll_real64                                :: Lx, Ly, Lz
  sll_real64                                :: dx, dy, dz
  sll_real64                                :: x, y, z
  sll_real64, dimension(:,:,:), allocatable :: rho
  sll_real64, dimension(:,:,:), allocatable :: phi_an
  sll_real64, dimension(:,:,:), allocatable :: phi
  sll_int32                                 :: i, j, k
  type (poisson_3d_periodic_plan), pointer  :: plan
  sll_real64                                :: average_err
  sll_int32, dimension(1:3)                 :: global
  sll_int32                                 :: gi, gj, gk
  sll_int32                                 :: myrank
  sll_real32                                :: ok = 1.d0
  sll_real32, dimension(1)                  :: prod4test
  sll_int32                                 :: npx, npy, npz
  type(layout_3D_t), pointer                :: layout
  sll_int64                                 :: colsz ! collective size
  sll_int32                                 :: i_test

  !Boot parallel environment
  call sll_boot_collective()

  nx = 128
  ny = 128
  nz = 128
  Lx = 2*sll_pi
  Ly = 2*sll_pi
  Lz = 2*sll_pi

  colsz  = sll_get_collective_size(sll_world_collective)

  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz

  SLL_ALLOCATE(phi_an(nx,ny,nz),ierr)
  SLL_ALLOCATE(phi(nx,ny,nz),ierr)
  SLL_ALLOCATE(rho(nx,ny,nz),ierr)

  plan => new_poisson_3d_periodic_plan(cmplx(rho, 0_f64, kind=f64), Lx, Ly, Lz)

  do i_test = 1, 2

     do k=1,nz
        z = (k-1)*dz
        do j=1,ny
           y = (j-1)*dy
           do i=1,nx
              x = (i-1)*dx
              if (i_test==1) then
                 phi_an(i,j,k) = cos(x)*sin(y)*cos(z)
                 rho(i,j,k) = 3*phi_an(i,j,k)
              else if (i_test == 2) then
                 phi_an(i,j,k) = (4/(sll_pi*sqrt(sll_pi)*Lx*Ly*Lz)) &
                                 * exp(-.5*(x-Lx/2)**2) * &
                                   exp(-.5*(y-Ly/2)**2) * sin(z)
                 rho(i,j,k) = phi_an(i,j,k) * ( 3 - ( (x-Lx/2)**2 + (y-Ly/2)**2 ) )
              end if
           enddo
        enddo
     enddo


     colsz  = sll_get_collective_size(sll_world_collective)
     myrank = sll_get_collective_rank(sll_world_collective)

     call solve_poisson_3d_periodic_par(plan, rho, phi)

     nx_loc = size(phi,1)
     ny_loc = size(phi,2)
     nz_loc = size(phi,3)
     npx = nx / nx_loc
     npy = ny / ny_loc
     npz = nz / nz_loc

     layout  => new_layout_3D( sll_world_collective ) 
     call initialize_layout_with_distributed_3D_array( nx, ny, nz, npx, npy, npz, layout )

     average_err  = 0.d0

     do k=1,nz_loc
        do j=1,ny_loc
           do i=1,nx_loc
              global = local_to_global_3D( layout, (/i, j, k/))
              gi = global(1)
              gj = global(2)
              gk = global(3)
              average_err  = average_err  + abs( phi_an (gi,gj,gk) - phi(i,j,k) )
           enddo
        enddo
     enddo

     average_err  = average_err  / (nx_loc*ny_loc*nz_loc)

     call flush()
     print*, ' myrank '
     print*, 'local average error:', average_err
     print*, 'dx*dy*dz =', dx*dy*dz

    if ( average_err > dx*dx*dy) then
       call flush()
       print*, ' '
       print*, 'Test stoppped by sll_poisson_3d_periodic_par test'
       print*, 'myrank=', myrank
       print*, ' '
       stop
    endif

    call sll_collective_reduce(sll_world_collective, (/ ok /), 1, MPI_PROD, 0, prod4test )

    if (myrank==0) then
       if (prod4test(1)==1.d0) then
          call flush()
          print*, ' '
          print*, 'sll_poisson_3d_periodic_par test: PASS'
          print*, ' '
       endif
    endif

  end do

  call delete_poisson_3d_periodic_plan(plan)


  SLL_DEALLOCATE_ARRAY(phi, ierr)
  SLL_DEALLOCATE_ARRAY(phi_an, ierr)
  SLL_DEALLOCATE_ARRAY(rho, ierr)
  call sll_halt_collective()

end program test_poisson_3d_par

