
!***************************************************************************
!
! Selalib 2012     
! Module: unit_test.F90
!
!> @brief 
!> Selalib poisson solvers (1D, 2D and 3D) unit test
!> Start date: March 20, 2012
!> Last modification: March 20, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************
program test_quasi_neutral
 
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
#include "sll_remap.h"

  use sll_qns_2d_with_finite_diff_seq
  use numeric_constants
  use sll_qns_2d_with_finite_diff_util
  use sll_qns_2d_with_finite_diff_seq
 ! use sll_qns_2d_with_finite_diff_par
  use sll_collective

implicit none

  character(len=100)                      :: bc ! Boundary_conditions
  sll_int32                               :: nr, ntheta
  sll_real64                              :: rmin, rmax, Zi
  sll_real64, dimension(:),   allocatable :: Te
  sll_int32                               :: ierr
  !Boot parallel environment
  call sll_boot_collective()

  bc = 'neumann'
  nr = 64
  ntheta = 64
  rmin = 1.d0
  rmax = 10.d0
  SLL_ALLOCATE(Te(nr), ierr)
  Te = 1.d0
  Zi = 1.d0

  call test_sll_qns_2d_with_finite_diff(bc, nr, ntheta, rmin, rmax, Te, Zi)

  SLL_DEALLOCATE_ARRAY(Te, ierr)

  call sll_halt_collective()

contains


  subroutine test_sll_qns_2d_with_finite_diff(bc, nr, ntheta, rmin, rmax, Te, Zi)

    character(len=*)                                 :: bc ! Boundary_conditions
    sll_int32                                        :: nr, ntheta
    sll_real64                                       :: rmin, rmax, Zi
    sll_real64, dimension(:),   allocatable          :: c, Te, f, g    
    sll_int32                                        :: nr_loc, ntheta_loc, nz_loc
    sll_int32                                        :: ierr
    sll_real64                                       :: dr, dtheta
    sll_real64                                       :: r, theta
    sll_real64, dimension(nr,ntheta)                 :: rho_seq
    sll_real64, dimension(:,:), allocatable          :: rho_par
    sll_real64, dimension(nr,ntheta)                 :: phi_an
    sll_real64, dimension(nr,ntheta)                 :: phi_seq
    sll_real64, dimension(:,:), allocatable          :: phi_par
    sll_int32                                        :: i, j
    type (qns_2d_with_finite_diff_plan_seq), pointer :: plan_seq
    type (qns_2d_with_finite_diff_plan_par), pointer :: plan_par
    sll_real64                                       :: average_err
    sll_real64                                       :: seq_par_diff
    sll_int32, dimension(1:3)                        :: global
    sll_int32                                        :: gi, gj
    sll_int32                                        :: myrank
    sll_real32                                       :: ok = 1.d0
    sll_real32, dimension(1)                         :: prod4test
    type(layout_3D_t), pointer                       :: layout
    sll_int64                                        :: colsz ! collective size
    sll_int32                                        :: npx, npy, npz
    sll_int32                                        :: e

    dr = (rmax-rmin)/(nr-1)
    dtheta = 2*sll_pi/ntheta

    f = 0.d0
    do j=1,ntheta
       theta = (j-1)*dtheta
        if (bc=='neumann') then
           f(j) = sin(rmax-rmin)*cos(theta)
        endif
       do i=1,nr
          global = local_to_global_3D( layout, (/int(i), 1, 1/) )
          gi = global(1)
          if (bc=='neumann') then
             r = rmin + (gi-1)*dr
             c(i) = 2/r
          else ! 'dirichlet'
             r = rmin + gi*dr
             c(i) = (rmax+rmin-2*r) / ( (rmax-r)*(r-rmin) )
          endif
          phi_an(i,j)  = sin(r-rmin)*sin(rmax-r)*cos(theta)
          rho_seq(i,j) = cos(theta)*(2*cos(rmin+rmax-2*r)-c(i)*sin(rmin+rmax-2*r)&
                         + (1/r**2+1/(Zi*Te(i)))*sin(rmax-r)*sin(r-rmin))
       enddo
    enddo
    g = -f

    colsz  = sll_get_collective_size(sll_world_collective)
    myrank = sll_get_collective_rank(sll_world_collective)

    ! Test sequential periodic 3D poisson solver
    if (myrank==0) then
       call flush()
       print*, ' '
       print*, 'Test poisson_3d in sequential'
    endif

    plan_seq => new_qns_2d_with_finite_diff_plan_seq(bc, rmin, rmax, rho_seq, c, Te, f, g, Zi)
    call solve_qn_2d_with_finite_diff_seq(plan_seq, phi_seq)

    average_err = sum( abs(phi_an-phi_seq) ) / (nr*ntheta)

    if (myrank==0) then
       call flush()
       print*, ' '
       call flush()
       print*, 'Average error:', average_err
       call flush()
       print*, 'dx*dy*dz =', dr*dtheta
    endif

    if ( average_err <= dr*dtheta) then
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

    colsz  = sll_get_collective_size(sll_world_collective)
    e = int(log(real(colsz))/log(2.))

    ! Layout and local sizes for FFTs in x-direction
    layout => new_layout_3D( sll_world_collective )
    npx = 1
    npy = 2**(e/2)
    npz = int(colsz)/npy
    call initialize_layout_with_distributed_3D_array( nr, ntheta, 1, npx, npy, npz, layout )

    plan_par => new_qns_2d_with_finite_diff_plan_par(layout, bc, rmin, rmax, rho_par, c, Te, f, g, Zi)

    !call compute_local_sizes( layout, nr_loc, ntheta_loc, nz_loc )
    SLL_ALLOCATE(rho_par(nr_loc,ntheta_loc), ierr)

    do j=1,ntheta_loc
       do i=1,nr_loc
          global = local_to_global_3D( layout, (/i, j, 1/))
          gi = global(1)
          gj = global(2)
          rho_par(i,j) = rho_seq(gi,gj)
        enddo
    enddo
    
    SLL_ALLOCATE(phi_par(nr_loc,ntheta_loc), ierr)

    !call solve_qn_2d_with_finite_diff_par(plan_par, phi_par)

    average_err  = 0.d0
    seq_par_diff = 0.d0

    do j=1,ntheta_loc
       do i=1,nr_loc
          global = local_to_global_3D( layout, (/i, j, 1/))
          gi = global(1)
          gj = global(2)
          average_err  = average_err  + abs( phi_an (gi,gj) - phi_par(i,j) )
          seq_par_diff = seq_par_diff + abs( phi_seq(gi,gj) - phi_par(i,j) )
       enddo
    enddo

    average_err  = average_err  / (nr_loc*ntheta_loc)
    seq_par_diff = seq_par_diff / (nr_loc*ntheta_loc)

    call flush()
    print*, ' '
    call flush()
    print*, 'Average error in proc', myrank, ':', average_err
    call flush()
    print*, 'dx*dy*dz =', dr*dtheta
    call flush()
    print*, 'Average diff between seq sol and par sol in proc', myrank, ':', seq_par_diff


    if ( average_err > dr*dtheta) then
       ok = 1.d0
       call flush()
       print*, ' '
       call flush()
       print*, 'Test stoppped by sll_poisson_3d_periodic_par test'
       call flush()
       print*, 'myrank=', myrank
       call flush()
       print*, ' '
       stop
    endif

    call sll_collective_reduce(sll_world_collective, (/ ok /), 1, MPI_PROD, 0, prod4test )

    if (myrank==0) then
       if (prod4test(1)==1.d0) then
          call flush()
          print*, ' '
          call flush()
          print*, 'sll_poisson_3d_periodic_par test: PASS'
          call flush()
          print*, ' '
       endif
    endif

    call delete_new_qns_2d_with_finite_diff_plan_par(plan_par)
    call delete_new_qns_2d_with_finite_diff_plan_seq(plan_seq)

    SLL_DEALLOCATE_ARRAY(phi_par, ierr)
    SLL_DEALLOCATE_ARRAY(phi_par, ierr)
 
  end subroutine test_sll_qns_2d_with_finite_diff


end program test_quasi_neutral
