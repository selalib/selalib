!> @brief 
!> Selalib poisson solvers (1D, 2D and 3D) unit test
!> Start date: March 20, 2012
!> Last modification (decoupling the tests): October 26, 2012
!>   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
program test_qn_solver_2d_parallel
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use sll_remapper
use sll_constants
use sll_collective
use sll_qn_solver_2d_parallel
use sll_boundary_condition_descriptors

implicit none

  sll_int32, parameter                    :: np_r = 256
  sll_int32, parameter                    :: np_theta = 256
  sll_int32                               :: bc 
  sll_real64                              :: rmin, rmax, Zi
  sll_real64, dimension(np_r)             :: Te
  sll_real64, dimension(np_r,np_theta)    :: rho_seq
  sll_real64, dimension(np_r,np_theta)    :: phi_exact
  sll_real64, dimension(np_theta)         :: f
  sll_real64, dimension(np_theta)         :: g
  sll_int32                               :: i, j, myrank
  sll_real32, dimension(1)                :: prod4test
  sll_int32                               :: np_r_loc, np_theta_loc
  sll_int64                               :: colsz ! collective size
  type(layout_2D), pointer                :: layout
  sll_int32                               :: ierr
  sll_real64, dimension(:,:), allocatable :: rho
  sll_real64, dimension(:,:), allocatable :: phi
  sll_real64                              :: r, dr, dtheta
  sll_real64, dimension(np_theta)         :: theta

  !Boot parallel environment
  call sll_boot_collective()

  myrank = sll_get_collective_rank(sll_world_collective)
  colsz  = sll_get_collective_size(sll_world_collective)

  np_r_loc = np_r/int(colsz)
  np_theta_loc = np_theta

  SLL_ALLOCATE(rho(np_r_loc,np_theta_loc), ierr)
  SLL_ALLOCATE(phi(np_r_loc,np_theta_loc), ierr)

  layout => new_layout_2D( sll_world_collective )
  call initialize_layout_with_distributed_2D_array( np_r,       &
                                                    np_theta,   &
                                                    int(colsz), &
                                                    1,          &
                                                    layout )

  rmin = 1.d0
  rmax = 10.d0
  Zi = 1.d0

  dtheta = 2*sll_pi/np_theta
  do j = 1, np_theta
     theta(j) = (j-1)*dtheta
  end do

  do bc=1,2

     Te = 1.d0
     call test_process(Te, Zi, prod4test)

  enddo

  if (myrank==0) then
     if (prod4test(1)==1.d0) then
        call flush(6)
        print*, 'PASSED'
        call flush(6)
        print*, ' '
     endif
  endif

  call sll_halt_collective()

contains

  subroutine test_process(Te_seq, Zi, prod4test)

    type (qn_solver_2d_parallel), pointer   :: plan

    sll_real64                              :: Zi
    sll_real64, dimension(:)                :: Te_seq
    sll_real64, dimension(:),   allocatable :: c_seq
    sll_real64, dimension(:),   allocatable :: c_par
    sll_real64, dimension(:),   allocatable :: Te_par
    sll_int32                               :: i_test
    sll_real64                              :: err
    sll_real64                              :: err_bound
    sll_real64                              :: Mr
    sll_real64                              :: Mtheta
    sll_int32, dimension(2)                 :: global
    sll_int32                               :: gi
    sll_int32                               :: gj
    sll_int32                               :: myrank
    sll_real32                              :: ok = 1.d0
    sll_real32, dimension(1)                :: prod4test

    if (bc==SLL_NEUMANN) then
       dr = (rmax-rmin)/(np_r-1)
    else 
       dr = (rmax-rmin)/(np_r+1)
    endif

    SLL_ALLOCATE(c_seq(np_r), ierr)
    SLL_ALLOCATE(c_par(np_r_loc), ierr)
    SLL_ALLOCATE(Te_par(np_r_loc), ierr)

    plan => new(bc,rmin,rmax,np_r, np_theta)

    do i_test=1,2 ! 2 test functions

    f = 0.d0
    err_bound = 0.d0

    do j=1,np_theta

       Mr = 4*abs(cos(theta(j)))
       if (bc==SLL_NEUMANN) then
          if (i_test==1) then
             f(j) = sin(rmax-rmin)*cos(theta(j))
          else
             f(j)= sin(rmax-rmin)*exp(-.5*(theta(j)-sll_pi)**2)/sqrt(2*sll_pi)
          endif
       endif

       do i=1,np_r
          if (bc==SLL_NEUMANN) then
             r = rmin + (i-1)*dr
             c_seq(i) = 2/r
          else ! 'dirichlet'
             r = rmin + i*dr
             c_seq(i) = (rmax+rmin-2*r) / ( (rmax-r)*(r-rmin) )
          endif
          ! c=n_0'(r), c_seq is c in sequential
          if (i_test==1) then
             phi_exact(i,j)  = sin(r-rmin)*sin(rmax-r)*cos(theta(j))
             rho_seq(i,j) = cos(theta(j)) * ( 2*cos(rmin+rmax-2*r) - c_seq(i)* sin( &
                 rmin+rmax-2*r)+(1/r**2+1/(Zi*Te_seq(i)))*sin(rmax-r)*sin(r-rmin))
          else
             phi_exact(i,j)  = sin(r-rmin)*sin(rmax-r)*exp(-.5*(theta(j)-sll_pi)**2)/ &
                                                                   sqrt(2*sll_pi)
             rho_seq(i,j) = ( 2*cos(rmax+rmin-2*r) - c_seq(i)*sin(rmax+rmin-2*r) ) * &
                                         exp(-.5*(theta(j)-sll_pi)**2)/sqrt(2*sll_pi) + &
             phi_exact(i,j) * ( 1/(Zi*Te_seq(i)) - ((theta(j)-sll_pi)**2-1)/r**2 )
          endif
          Mtheta = abs(sin(r-rmin)*sin(rmax-r))
          err_bound = err_bound + &
          Mr*dr**2/12 + abs(c_seq(i))*Mr*dr**2/6 + Mtheta*dtheta**2/(r**2*12)
       enddo

    enddo

    g = -f

    do j=1,np_theta_loc
       do i=1,np_r_loc
          global = local_to_global_2D( layout, (/i, j/))
          gi = global(1)
          gj = global(2)
          rho(i,j) = rho_seq(gi,gj)
          c_par(i) = c_seq(gi)
          Te_par(i) = Te_seq(gi)
        enddo
    enddo
   
    call solve(plan, rho, c_par, Te_par, f, g, Zi, phi)

    err        = 0.d0
    err_bound  = 0.d0

    do j=1,np_theta_loc
       do i=1,np_r_loc
          global = local_to_global_2D(layout, (/i, j/))
          gi = global(1)
          gj = global(2)
          if (bc== SLL_NEUMANN) then
             r = rmin + (gi-1)*dr
          else 
             r = rmin + gi*dr
          endif
          err = err  + abs( phi_exact (gi,gj) - phi(i,j))
          Mr = 4*abs(cos(theta(gj)))
          Mtheta = abs(sin(r-rmin)*sin(rmax-r))
          err_bound = err_bound + Mr*dr**2/12 + &
               abs(c_par(i))*Mr*dr**2/6 + Mtheta*dtheta**2/(r**2*12)
       enddo
    enddo

    err  = err/(np_r_loc*np_theta_loc)
    err_bound = err_bound/(np_r_loc*np_theta_loc)

    call flush(6)
    print*, 'Error in proc', myrank, ':', err
    call flush(6)
    print*, 'Boundary error in proc', myrank, ':', err_bound
    call flush(6)
    print*, ' '

    if ( err > err_bound) then
       call flush(6)
       print*, 'myrank=', myrank
       call flush(6)
       print*, 'FAILED'
       stop
    endif

    call sll_collective_reduce(sll_world_collective, (/ ok /), 1,        &
                                                  MPI_PROD, 0, prod4test )
    enddo

    call delete(plan)

    SLL_DEALLOCATE_ARRAY(c_seq, ierr)
    SLL_DEALLOCATE_ARRAY(c_par, ierr)
    SLL_DEALLOCATE_ARRAY(Te_par, ierr)
 
  end subroutine test_process


end program test_qn_solver_2d_parallel
