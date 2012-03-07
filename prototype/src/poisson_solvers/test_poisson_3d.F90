!***************************************************************************
!
! Selalib 2012     
! File : test_poisson_3d.F90
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

program test_poisson_3d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
#include "sll_poisson_solvers.h"
#include "sll_remap.h"

  use numeric_constants
  use sll_poisson_3d_periodic_util
  use sll_poisson_3d_periodic_seq
  use sll_poisson_3d_periodic_par
  use sll_collective

  implicit none

  sll_int32  :: nx, ny, nz
  sll_real64 :: Lx, Ly, Lz
  sll_int64  :: colsz

  !Boot parallel environment
  call sll_boot_collective()

  nx = 128
  ny = 128
  nz = 128
  Lx = 2*sll_pi
  Ly = 2*sll_pi
  Lz = 2*sll_pi

  colsz  = sll_get_collective_size(sll_world_collective)

  call test_sll_poisson_3d_periodic(nx, ny, nz, Lx, Ly, Lz)

  call sll_halt_collective()

contains

  subroutine test_sll_poisson_3d_periodic(nx, ny, nz, Lx, Ly, Lz)

    sll_int32                                 :: nx, ny, nz
    sll_int32                                 :: nx_loc, ny_loc, nz_loc
    sll_int32                                 :: ierr
    sll_real64                                :: Lx, Ly, Lz
    sll_real64                                :: dx, dy, dz
    sll_real64                                :: x, y, z
    sll_real64, dimension(nx,ny,nz)           :: rho1, phi_an1, phi_seq1
    sll_real64, dimension(nx,ny,nz)           :: rho2, phi_an2, phi_seq2
    sll_real64, dimension(:,:,:), allocatable :: phi_par1, phi_par2
    sll_int32                                 :: i, j, k
    type (poisson_3d_periodic_plan), pointer  :: plan
    sll_real64                                :: average_err1, average_err2
    sll_real64                                :: seq_par_diff1, seq_par_diff2
    sll_int32, dimension(1:3)                 :: global
    sll_int32                                 :: gi, gj, gk
    sll_int32                                 :: myrank
    sll_real32                                :: ok = 1.d0
    sll_real32, dimension(1)                  :: prod4test
    sll_int32                                 :: npx, npy, npz
    type(layout_3D_t), pointer                :: layout
    sll_int64                                 :: colsz ! collective size

    dx = Lx/nx
    dy = Ly/ny
    dz = Lz/nz

    do k=1,nz
       z = (k-1)*dz
       do j=1,ny
          y = (j-1)*dy
          do i=1,nx
             x = (i-1)*dx
             phi_an1(i,j,k) = cos(x)*sin(y)*cos(z)
             phi_an2(i,j,k) = (4/(sll_pi*sqrt(sll_pi)*Lx*Ly*Lz)) * &
                              exp(-.5*(x-Lx/2)**2) *               &
                              exp(-.5*(y-Ly/2)**2) * sin(z)
             rho2(i,j,k) = phi_an2(i,j,k) * (3.0_f64 - ((x-Lx/2.0_f64)**2 + (y-Ly/2.0_f64)**2))
          enddo
       enddo
    enddo

    rho1 = 3*phi_an1

    colsz  = sll_get_collective_size(sll_world_collective)
    myrank = sll_get_collective_rank(sll_world_collective)

    ! Test sequential periodic 3D poisson solver
    if (myrank==0) then
       call flush()
       print*, ' '
       print*, 'Test poisson_3d in sequential'
    endif

    plan => new_poisson_3d_periodic_plan(cmplx(rho1, 0_f64, kind=f64), Lx, Ly, Lz)
    call solve_poisson_3d_periodic_seq(plan, rho1, phi_seq1)
    call solve_poisson_3d_periodic_seq(plan, rho2, phi_seq2)

    average_err1 = sum( abs(phi_an1-phi_seq1) ) / (nx*ny*nz)
    average_err2 = sum( abs(phi_an2-phi_seq2) ) / (nx*ny*nz)

    if (myrank==0) then
       call flush()
       print*, ' '
       print*, 'Average error1, Average error2:', average_err1, average_err2
       print*, 'dx*dy*dz =', dx*dy*dz
    endif

    if ( max(average_err1, average_err2) <= dx*dx*dy) then
       if (myrank==0) then
          call flush()
          print*, ' '
          print*, 'sll_poisson_3d_periodic_seq test: PASS'
       endif
    else
       call flush()
       print*, ' '
       print*, 'Test stoppped by sll_poisson_3d_periodic_seq test'
       print*, ' '
       stop
    endif

    ! Test parallel periodic 3D poisson solver

    if (myrank==0) then
       call flush()
       print*, ' '
       call flush()
       print*, 'Test poisson_3d in parallel'
    endif

    call solve_poisson_3d_periodic_par(plan, rho1, phi_par1)
    call solve_poisson_3d_periodic_par(plan, rho2, phi_par2)

    nx_loc = size(phi_par1,1)
    ny_loc = size(phi_par1,2)
    nz_loc = size(phi_par1,3)
    npx = nx / nx_loc
    npy = ny / ny_loc
    npz = nz / nz_loc

    layout  => new_layout_3D( sll_world_collective ) 
    call initialize_layout_with_distributed_3D_array( nx, ny, nz, npx, npy, npz, layout )

    average_err1  = 0.d0
    seq_par_diff1 = 0.d0
    average_err2  = 0.d0
    seq_par_diff2 = 0.d0

    do k=1,nz_loc
       do j=1,ny_loc
          do i=1,nx_loc
             global = local_to_global_3D( layout, (/i, j, k/))
             gi = global(1)
             gj = global(2)
             gk = global(3)
             average_err1  = average_err1  + abs( phi_an1 (gi,gj,gk) - phi_par1(i,j,k) )
             seq_par_diff1 = seq_par_diff1 + abs( phi_seq1(gi,gj,gk) - phi_par1(i,j,k) )
             average_err2  = average_err2  + abs( phi_an2 (gi,gj,gk) - phi_par2(i,j,k) )
             seq_par_diff2 = seq_par_diff2 + abs( phi_seq2(gi,gj,gk) - phi_par2(i,j,k) )
          enddo
       enddo
    enddo

    average_err1  = average_err1  / (nx_loc*ny_loc*nz_loc)
    seq_par_diff1 = seq_par_diff1 / (nx_loc*ny_loc*nz_loc)
    average_err2  = average_err2  / (nx_loc*ny_loc*nz_loc)
    seq_par_diff2 = seq_par_diff2 / (nx_loc*ny_loc*nz_loc)

    call flush()
    print*, ' '
    print*, 'local average error1, local average error2:', average_err1, average_err2
    print*, 'dx*dy*dz =', dx*dy*dz
    print*, 'Local average diff between seq sol 1 and par sol 2:', seq_par_diff1, seq_par_diff2

    if ( max(average_err1, average_err2) > dx*dx*dy) then
       ok = 1.d0
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

    call delete_poisson_3d_periodic_plan(plan)
    SLL_DEALLOCATE_ARRAY(phi_par1, ierr)
    SLL_DEALLOCATE_ARRAY(phi_par2, ierr)

  end subroutine test_sll_poisson_3d_periodic

end program test_poisson_3d
