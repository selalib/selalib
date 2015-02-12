
!***************************************************************************
!
! Selalib 2012     
! Module: unit_test.F90
!
!> @brief 
!> Selalib poisson solvers (1D, 2D and 3D) unit test
!> Start date: March 20, 2012
!> Last modification (decoupling the tests): October 26, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************
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

  sll_int32                             :: BC ! Boundary_conditions
  sll_int32                             :: NP_r, NP_theta
  ! NP_r and NP_theta are the numbers of points in directions r and 
  ! theta respectively
  sll_real64                            :: rmin, rmax, Zi
  sll_real64, dimension(:), allocatable :: Te
  sll_int32                             :: ierr, i, myrank
  sll_real32, dimension(1)              :: prod4test

  !Boot parallel environment
  call sll_boot_collective()

  myrank = sll_get_collective_rank(sll_world_collective)

  NP_r = 256
  NP_theta = 256
  rmin = 1.d0
  rmax = 10.d0
  Zi = 1.d0

  do i=1,2

     if (i==1) then
        BC = SLL_NEUMANN
     else
        BC = SLL_DIRICHLET
     endif
     if (myrank==0) then
        call flush(6)
        print*, ' '
        call flush(6)
        print*, 'Testing sll_qns2d_angular_spect_method_par with ', BC
        call flush(6)
        print*, ' '
     endif
     SLL_ALLOCATE(Te(NP_r), ierr)
     Te = 1.d0
     call test_process(BC, NP_r, NP_theta, rmin, rmax, Te, Zi, prod4test)
     SLL_DEALLOCATE_ARRAY(Te, ierr)

  enddo

  if (myrank==0) then
     if (prod4test(1)==1.d0) then
        call flush(6)
        print*, 'test_sll_qns2d_angular_spect_method_par: PASSED'
        call flush(6)
        print*, ' '
     endif
  endif

  call sll_halt_collective()

contains


  subroutine test_process(BC, NP_r, NP_theta, rmin, rmax, Te_seq, Zi, prod4test)

    sll_int32                                           :: BC ! Boundary_conditions
    sll_int32                                           :: NP_r, NP_theta
    ! NP_r and NP_theta are the numbers of points in directions r and 
    ! theta respectively
    sll_real64                                          :: rmin, rmax, Zi
    sll_real64, dimension(:)                            :: Te_seq
    sll_real64, dimension(:),   allocatable             :: c_seq, f, g, c_par, Te_par
    sll_int32                                           :: NP_r_loc, NP_theta_loc
    ! NP_r_loc and NP_theta_loc are the numbers of points locally (in the 
    ! processor) in directions r and theta respectively
    sll_int32                                           :: ierr
    sll_real64                                          :: dr, dtheta
    sll_real64                                          :: r, theta
    sll_real64, dimension(NP_r,NP_theta)                :: rho_seq
    sll_real64, dimension(NP_r,NP_theta)                :: phi_exact
    sll_real64, dimension(:,:), allocatable             :: rho_par, phi
    sll_int32                                           :: i, j, i_test
    type (qn_solver_2d_parallel), pointer :: plan
    sll_real64                                          :: average_err
    sll_real64                                          :: average_err_bound
    sll_real64                                          :: Mr, Mtheta
    sll_int32, dimension(1:3)                           :: global
    sll_int32                                           :: gi, gj
    sll_int32                                           :: myrank
    sll_real32                                          :: ok = 1.d0
    sll_real32, dimension(1)                            :: prod4test
    type(layout_3D), pointer                            :: layout
    sll_int64                                           :: colsz ! collective size

    if (BC==SLL_NEUMANN) then
       dr = (rmax-rmin)/(NP_r-1)
    else ! 'Dirichlet'
       dr = (rmax-rmin)/(NP_r+1)
    endif
    dtheta = 2*sll_pi/NP_theta

    colsz  = sll_get_collective_size(sll_world_collective)
    myrank = sll_get_collective_rank(sll_world_collective)

    NP_r_loc = NP_r/int(colsz)
    NP_theta_loc = NP_theta

    SLL_ALLOCATE(rho_par(NP_r_loc,NP_theta_loc), ierr)
    SLL_ALLOCATE(c_seq(NP_r), ierr)
    SLL_ALLOCATE(c_par(NP_r_loc), ierr)
    SLL_ALLOCATE(Te_par(NP_r_loc), ierr)
    SLL_ALLOCATE(phi(NP_r_loc,NP_theta_loc), ierr)
    SLL_ALLOCATE(f(NP_theta), ierr)
    SLL_ALLOCATE(g(NP_theta), ierr)

    plan => new_qn_solver_2d_parallel(BC,rmin,rmax,NP_r, NP_theta)

    do i_test=1,2 ! 2 test functions

    f = 0.d0
    average_err_bound = 0.d0

    do j=1,NP_theta

       theta = (j-1)*dtheta
       Mr = 4*abs(cos(theta))
       if (BC==SLL_NEUMANN) then
          if (i_test==1) then
             f(j) = sin(rmax-rmin)*cos(theta)
          else
             f(j)= sin(rmax-rmin) * exp(-.5*(theta-sll_pi)**2)/sqrt(2*sll_pi)
          endif
       endif

       do i=1,NP_r
          if (BC==SLL_NEUMANN) then
             r = rmin + (i-1)*dr
             c_seq(i) = 2/r
          else ! 'dirichlet'
             r = rmin + i*dr
             c_seq(i) = (rmax+rmin-2*r) / ( (rmax-r)*(r-rmin) )
          endif
          ! c=n_0'(r), c_seq is c in sequential
          if (i_test==1) then
             phi_exact(i,j)  = sin(r-rmin)*sin(rmax-r)*cos(theta)
             rho_seq(i,j) = cos(theta) * ( 2*cos(rmin+rmax-2*r) - c_seq(i)* sin( &
                 rmin+rmax-2*r)+(1/r**2+1/(Zi*Te_seq(i)))*sin(rmax-r)*sin(r-rmin))
          else
             phi_exact(i,j)  = sin(r-rmin)*sin(rmax-r)*exp(-.5*(theta-sll_pi)**2)/ &
                                                                   sqrt(2*sll_pi)
             rho_seq(i,j) = ( 2*cos(rmax+rmin-2*r) - c_seq(i)*sin(rmax+rmin-2*r) ) * &
                                         exp(-.5*(theta-sll_pi)**2)/sqrt(2*sll_pi) + &
                    phi_exact(i,j) * ( 1/(Zi*Te_seq(i)) - ((theta-sll_pi)**2-1)/r**2 )
          endif
          Mtheta = abs(sin(r-rmin)*sin(rmax-r))
          average_err_bound = average_err_bound + &
          Mr*dr**2/12 + abs(c_seq(i))*Mr*dr**2/6 + Mtheta*dtheta**2/(r**2*12)
       enddo

    enddo

    g = -f

    ! Test sll_qns2d_angular_spect_method_par

    layout => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( NP_r, NP_theta, 1, &
                                                int(colsz), 1, 1, layout )

    do j=1,NP_theta_loc
       do i=1,NP_r_loc
          global = local_to_global_3D( layout, (/i, j, 1/))
          gi = global(1)
          gj = global(2)
          rho_par(i,j) = rho_seq(gi,gj)
          c_par(i) = c_seq(gi)
          Te_par(i) = Te_seq(gi)
        enddo
    enddo
   
    call solve_qn_solver_2d_parallel(plan, rho_par, c_par, Te_par, f, g, Zi, phi)

    average_err        = 0.d0
    average_err_bound  = 0.d0

    do j=1,NP_theta_loc
       do i=1,NP_r_loc
          global = local_to_global_3D(layout, (/i, j, 1/))
          gi = global(1)
          gj = global(2)
          theta = (gj-1)*dtheta
          if (bc==SLL_NEUMANN) then
             r = rmin + (gi-1)*dr
          else ! 'dirichlet'
             r = rmin + gi*dr
          endif
          average_err = average_err  + abs( phi_exact (gi,gj) - phi(i,j))
          Mr = 4*abs(cos(theta))
          Mtheta = abs(sin(r-rmin)*sin(rmax-r))
          average_err_bound = average_err_bound + Mr*dr**2/12 + &
               abs(c_par(i))*Mr*dr**2/6 + Mtheta*dtheta**2/(r**2*12)
       enddo
    enddo

    average_err  = average_err/(NP_r_loc*NP_theta_loc)
    average_err_bound = average_err_bound/(NP_r_loc*NP_theta_loc)

    call flush(6)
    print*, 'Error in proc', myrank, ':', average_err
    call flush(6)
    print*, 'Boundary error in proc', myrank, ':', average_err_bound
    call flush(6)
    print*, ' '
    if ( average_err > average_err_bound) then
       call flush(6)
       print*, 'test_sll_qns2d_angular_spect_method_par: FAILED'
       call flush(6)
       print*, 'myrank=', myrank
       call flush(6)
       print*, ' '
       stop
    endif

    call sll_collective_reduce(sll_world_collective, (/ ok /), 1,        &
                                                  MPI_PROD, 0, prod4test )
    enddo

    call delete_qn_solver_2d_parallel(plan)

    SLL_DEALLOCATE_ARRAY(phi, ierr)
    SLL_DEALLOCATE_ARRAY(c_seq, ierr)
    SLL_DEALLOCATE_ARRAY(c_par, ierr)
    SLL_DEALLOCATE_ARRAY(Te_par, ierr)
    SLL_DEALLOCATE_ARRAY(f, ierr)
    SLL_DEALLOCATE_ARRAY(g, ierr)
 
  end subroutine test_process


end program test_qn_solver_2d_parallel
