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

#include "sll_remap.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
#include "sll_poisson_solvers.h"

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
  sll_real64, dimension(:,:,:), allocatable :: x, y, z
  sll_real64, dimension(:,:,:), allocatable :: rho
  sll_real64, dimension(:,:,:), allocatable :: phi_an
  sll_real64, dimension(:,:,:), allocatable :: phi
  sll_int32                                 :: i, j, k
  type (poisson_3d_periodic_plan_par), pointer  :: plan
  sll_real64                                :: average_err
  sll_int32, dimension(1:3)                 :: global
  sll_int32                                 :: gi, gj, gk
  sll_int32                                 :: myrank
  sll_real32                                :: ok = 1.d0
  sll_real32, dimension(1)                  :: prod4test
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

  SLL_ALLOCATE(x(nx,ny,nz),ierr)
  SLL_ALLOCATE(y(nx,ny,nz),ierr)
  SLL_ALLOCATE(z(nx,ny,nz),ierr)
  SLL_ALLOCATE(phi_an(nx,ny,nz),ierr)

  plan => new_poisson_3d_periodic_plan_par(nx, ny, nz, Lx, Ly, Lz)

  layout => plan%layout_x
  call compute_local_sizes( layout, nx_loc, ny_loc, nz_loc )
  SLL_ALLOCATE(rho(nx_loc,ny_loc,nz_loc), ierr)

  layout => plan%layout_z
  call compute_local_sizes( layout, nx_loc, ny_loc, nz_loc )
  SLL_ALLOCATE(phi(nx_loc,ny_loc,nz_loc), ierr)

  do k=1,nz
     do j=1,ny
        do i=1,nx
           x(i,j,k) = (i-1)*dx
           y(i,j,k) = (j-1)*dy
           z(i,j,k) = (k-1)*dz
        end do
     end do
  end do

  do i_test = 1, 2

     if (i_test==1) then
        phi_an = cos(x)*sin(y)*cos(z)
     else if (i_test == 2) then
        phi_an = (4/(sll_pi * sqrt(sll_pi)*Lx*Ly*Lz))        &
                            * exp(-.5*(x-Lx/2)**2)           &
                            * exp(-.5*(y-Ly/2)**2) * sin(z)
     end if

     do k=1,nz_loc
        do j=1,ny_loc
           do i=1,nx_loc
              global = local_to_global_3D( layout, (/i, j, k/))
              gi = global(1)
              gj = global(2)
              gk = global(3)
              rho(i,j,k) = 3*phi_an(gi,gj,gk)
              rho(i,j,k) = phi_an(gi,gj,gk) * ( 3 - ( (x(gi,gj,gk)-Lx/2)**2 + (y(gi,gj,gk)-Ly/2)**2 ) )
           enddo
        enddo
     enddo

    call solve_poisson_3d_periodic_par(plan, rho, phi)

    average_err  = 0.d0

    layout => plan%layout_z
    call compute_local_sizes( layout, nx_loc, ny_loc, nz_loc )
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

  call delete_poisson_3d_periodic_plan_par(plan)

  SLL_DEALLOCATE_ARRAY(phi, ierr)
  SLL_DEALLOCATE_ARRAY(phi_an, ierr)
  SLL_DEALLOCATE_ARRAY(rho, ierr)
  call sll_halt_collective()

contains

subroutine compute_local_sizes( layout, loc_sz_i, loc_sz_j, loc_sz_k )
    type(layout_3D_t), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32, intent(out) :: loc_sz_k
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to new_distributed_mesh_3D'
       print *, 'Exiting...'
       STOP
    end if
    my_rank = sll_get_collective_rank(get_layout_3D_collective(layout))
    i_min = get_layout_3D_i_min( layout, my_rank )
    i_max = get_layout_3D_i_max( layout, my_rank )
    j_min = get_layout_3D_j_min( layout, my_rank )
    j_max = get_layout_3D_j_max( layout, my_rank )
    k_min = get_layout_3D_k_min( layout, my_rank )
    k_max = get_layout_3D_k_max( layout, my_rank )
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
    loc_sz_k = k_max - k_min + 1
end subroutine compute_local_sizes

end program test_poisson_3d_par
