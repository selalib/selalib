
!***************************************************************************
!
! Selalib 2012     
! Module: unit_test.F90
!
!> @brief 
!> Selalib poisson solvers (1D, 2D and 3D) unit test
!> Start date: March 20, 2012
!> Last modification: May 03, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************
program qns_tester
 
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
#include "sll_remap.h"

  use numeric_constants
  use sll_collective
  use sll_qns_2d_with_finite_diff_seq
  use sll_qns_2d_with_finite_diff_par
  use sll_qns2d_angular_spect_method_seq
  use sll_qns2d_angular_spect_method_par
implicit none

  character(len=100)                    :: BC ! Boundary_conditions
  sll_int32                             :: NP_r, NP_theta
  ! NP_r and NP_theta are the numbers of points in directions r and theta respectively
  sll_real64                            :: rmin, rmax, Zi
  sll_real64, dimension(:), allocatable :: Te
  sll_int32                             :: ierr, i, myrank

  !Boot parallel environment
  call sll_boot_collective()

  myrank = sll_get_collective_rank(sll_world_collective)

  NP_r = 64
  NP_theta = 64
  rmin = 1.d0
  rmax = 10.d0
  Zi = 1.d0

  do i=2,2

     if (i==1) then
        BC = 'neumann'
     else
        BC = 'dirichlet'
     endif
     SLL_ALLOCATE(Te(NP_r), ierr)
     Te = 1.d0
     if (myrank==0) then
        call flush()
        print*, ' '
        print*, 'TESTING QUASINEUTRAL SOLVERS WITH ', BC
     endif
     call test_qns_2d(BC, NP_r, NP_theta, rmin, rmax, Te, Zi)
     SLL_DEALLOCATE_ARRAY(Te, ierr)

  enddo

  call sll_halt_collective()

contains


  subroutine test_qns_2d(BC, NP_r, NP_theta, rmin, rmax, Te_seq, Zi)

    character(len=*)                                 :: BC ! Boundary_conditions
    sll_int32                                        :: NP_r, NP_theta
    ! NP_r and NP_theta are the numbers of points in directions r and theta respectively
    sll_real64                                       :: rmin, rmax, Zi
    sll_real64, dimension(:)                         :: Te_seq
    sll_real64, dimension(:),   allocatable          :: c_seq, f, g, c_par, Te_par
    sll_int32                                        :: NP_r_loc, NP_theta_loc
    ! NP_r_loc and NP_theta_loc are the numbers of points locally (in the processor) 
    ! in directions r and theta respectively
    sll_int32                                        :: ierr
    sll_real64                                       :: dr, dtheta
    sll_real64                                       :: r, theta
    sll_real64, dimension(NP_r,NP_theta)             :: rho_seq, phi_seq
    sll_real64, dimension(NP_r,NP_theta)             :: phi_spect_seq
    sll_real64, dimension(NP_r,NP_theta)             :: phi_an
    sll_real64, dimension(:,:), allocatable          :: rho_par, phi_par
    sll_real64, dimension(:,:), allocatable          :: phi_spect_par
    sll_int32                                        :: i, j
    type (qns_2d_with_finite_diff_plan_seq), pointer :: plan_seq
    type (qns_2d_with_finite_diff_plan_par), pointer :: plan_par
    type (qns2d_angular_spect_method_seq),pointer    :: plan_spect_seq
    type (qns2d_angular_spect_method_par),pointer    :: plan_spect_par
    sll_real64                                       :: average_err
    sll_real64                                       :: average_err_spect
    sll_real64                                       :: average_err_bound
    sll_real64                                       :: seq_par_diff
    sll_real64                                       :: seq_par_diff_spect
    sll_real64                                       :: Mr, Mtheta
    sll_int32, dimension(1:3)                        :: global
    sll_int32                                        :: gi, gj
    sll_int32                                        :: myrank
    sll_real32                                       :: ok = 1.d0
    sll_real32, dimension(1)                         :: prod4test
    type(layout_3D_t), pointer                       :: layout
    sll_int64                                        :: colsz ! collective size

    if (BC=='neumann') then
       dr = (rmax-rmin)/(NP_r-1)
    else ! 'Dirichlet'
       dr = (rmax-rmin)/(NP_r+1)
    endif
    dtheta = 2*sll_pi/NP_theta

    SLL_ALLOCATE(c_seq(NP_r), ierr)
    SLL_ALLOCATE(f(NP_theta), ierr)
    SLL_ALLOCATE(g(NP_theta), ierr)

    f = 0.d0
    average_err_bound  = 0.d0

    do j=1,NP_theta

       theta = (j-1)*dtheta
       Mr = 4*abs(cos(theta))
       if (BC=='neumann') then
          f(j) = sin(rmax-rmin)*cos(theta)
       endif

       do i=1,NP_r
          if (BC=='neumann') then
             r = rmin + (i-1)*dr
             c_seq(i) = 2/r
          else ! 'dirichlet'
             r = rmin + i*dr
             c_seq(i) = (rmax+rmin-2*r) / ( (rmax-r)*(r-rmin) )
          endif
          phi_an(i,j)  = sin(r-rmin)*sin(rmax-r)*cos(theta)
          rho_seq(i,j) = cos(theta) * ( 2*cos(rmin+rmax-2*r) - c_seq(i)* sin( &
              rmin+rmax-2*r)+(1/r**2+1/(Zi*Te_seq(i)))*sin(rmax-r)*sin(r-rmin))
          Mtheta = abs(sin(r-rmin)*sin(rmax-r))
          average_err_bound = average_err_bound + &
          Mr*dr**2/12 + abs(c_seq(i))*Mr*dr**2/6 + Mtheta*dtheta**2/(r**2*12)
       enddo

    enddo

    g = -f

    colsz  = sll_get_collective_size(sll_world_collective)
    myrank = sll_get_collective_rank(sll_world_collective)

    ! Test 2d qns
    if (myrank==0) then
       call flush()
       print*, ' '
       print*, '... IN SEQUENTIAL...'
    endif

    plan_seq => new_qns_2d_with_finite_diff_plan_seq(BC, rmin, rmax, NP_r, NP_theta)
    call solve_qn_2d_with_finite_diff_seq(plan_seq, rho_seq, c_seq, Te_seq, f, g, &
                                                                       Zi, phi_seq)

    plan_spect_seq => new_qns2d_angular_spect_method_seq(BC, rmin, rmax, NP_r, NP_theta)
    call solve_qns2d_angular_spect_method_seq(plan_spect_seq, rho_seq, c_seq, Te_seq,  &
                                                                f, g, Zi, phi_spect_seq)

    average_err = sum(abs(phi_an-phi_seq))/(NP_r*NP_theta)
    average_err_spect = sum(abs(phi_an-phi_spect_seq))/(NP_r*NP_theta)
    average_err_bound = average_err_bound/(NP_r*NP_theta)

    if (myrank==0) then
       call flush()
       print*, ' '
       call flush()
       print*, 'sll_qns_2d_with_finite_diff_seq average error:', average_err
       call flush()
       print*, 'sll_qns2d_angular_spect_method_seq average error:', average_err_spect
       call flush()
      print*, 'Boundary average error =', average_err_bound
    endif

    if ( average_err <= average_err_bound ) then
       if (myrank==0) then
          call flush()
          print*, ' '
          print*, 'sll_qns_2d_with_finite_diff_seq: PASSED'
       endif
    else
       call flush()
       print*, ' '
       print*, 'Test stopped by sll_qns_2d_with_finite_diff_seq failure'
       print*, ' '
       stop
    endif

    if ( average_err_spect <= average_err_bound ) then
       if (myrank==0) then
          call flush()
          print*, ' '
          print*, 'sll_qns2d_angular_spect_method_seq: PASSED'
       endif
    else
       call flush()
       print*, ' '
       print*, 'Test stopped by sll_qns2d_angular_spect_method_seq failure'
       print*, ' '
       stop
    endif

    ! Test sll_qns_2d_with_finite_diff_par

    if (myrank==0) then
       call flush()
       print*, ' '
       call flush()
       print*, '...IN PARALLEL...'
    endif

    layout => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( NP_r, NP_theta, 1, &
                                                int(colsz), 1, 1, layout )

    NP_r_loc = NP_r/int(colsz)
    NP_theta_loc = NP_theta

    SLL_ALLOCATE(rho_par(NP_r_loc,NP_theta_loc), ierr)
    SLL_ALLOCATE(c_par(NP_r_loc), ierr)
    SLL_ALLOCATE(Te_par(NP_r_loc), ierr)
    SLL_ALLOCATE(phi_par(NP_r_loc,NP_theta_loc), ierr)
    SLL_ALLOCATE(phi_spect_par(NP_r_loc,NP_theta_loc), ierr)

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
   
    plan_par => new_qns_2d_with_finite_diff_plan_par(BC,rmin,rmax,NP_r, NP_theta)
    call solve_qn_2d_with_finite_diff_par(plan_par, rho_par, c_par, Te_par, f, &
                                                                   g, Zi, phi_par)

    plan_spect_par => new_qns2d_angular_spect_method_par(BC,rmin,rmax,NP_r, NP_theta)
    call solve_qns2d_angular_spect_method_par(plan_spect_par, rho_par, c_par, &
                                                Te_par, f, g, Zi, phi_spect_par)

    average_err        = 0.d0
    average_err_spect  = 0.d0
    average_err_bound  = 0.d0
    seq_par_diff       = 0.d0
    seq_par_diff_spect = 0.d0

    do j=1,NP_theta_loc
       do i=1,NP_r_loc
          global = local_to_global_3D(layout, (/i, j, 1/))
          gi = global(1)
          gj = global(2)
          theta = (gj-1)*dtheta
          if (bc=='neumann') then
             r = rmin + (gi-1)*dr
          else ! 'dirichlet'
             r = rmin + gi*dr
          endif
          average_err = average_err  + abs( phi_an (gi,gj) - phi_par(i,j))
          average_err_spect  = average_err_spect + abs( phi_an (gi,gj) - &
                                                      phi_spect_par(i,j) )
          Mr = 4*abs(cos(theta))
          Mtheta = abs(sin(r-rmin)*sin(rmax-r))
          average_err_bound = average_err_bound + Mr*dr**2/12 + &
               abs(c_par(i))*Mr*dr**2/6 + Mtheta*dtheta**2/(r**2*12)
          seq_par_diff = seq_par_diff + abs( phi_seq(gi,gj) - phi_par(i,j) )
          seq_par_diff_spect = seq_par_diff_spect + abs(phi_spect_seq(gi,gj) & 
                                                      - phi_spect_par(i,j) )
       enddo
    enddo

    average_err  = average_err/(NP_r_loc*NP_theta_loc)
    average_err_spect  = average_err_spect/(NP_r_loc*NP_theta_loc)
    average_err_bound = average_err_bound/(NP_r_loc*NP_theta_loc)
    seq_par_diff = seq_par_diff/(NP_r_loc*NP_theta_loc)
    seq_par_diff_spect = seq_par_diff_spect/(NP_r_loc*NP_theta_loc)

    call flush()
    print*, ' '
    call flush()
    print*, 'sll_qns_2d_with_finite_diff_par average error',          &
                                    ' in proc', myrank, ':', average_err
    call flush()
    print*, 'sll_qns2d_angular_spect_method_par average error',       &
                              ' in proc', myrank, ':', average_err_spect
    call flush()
    print*, 'Boundary average error =', average_err_bound
    call flush()
    print*, 'sll_qns_2d_with_finite_diff_par average diff(seq,par)',   &
                                   ' in proc', myrank, ':', seq_par_diff
    call flush()
    print*, 'sll_qns2d_angular_spect_method_par average diff(seq,par)',&
                             ' in proc', myrank, ':', seq_par_diff_spect

    if ( average_err > average_err_bound) then
       call flush()
       print*, ' '
       call flush()
       print*, 'Test stopped by sll_qns_2d_with_finite_diff_par failure'
       call flush()
       print*, 'myrank=', myrank
       call flush()
       print*, ' '
       stop
    endif
    if ( average_err_spect > average_err_bound) then
       call flush()
       print*, ' '
       call flush()
       print*, 'Test stopped by sll_qns2d_angular_spect_method_par failure'
       call flush()
       print*, 'myrank=', myrank
       call flush()
       print*, ' '
       stop
    endif

    call sll_collective_reduce(sll_world_collective, (/ ok /), 1,        &
                                                  MPI_PROD, 0, prod4test )
    if (myrank==0) then
       if (prod4test(1)==1.d0) then
          call flush()
          print*, ' '
          call flush()
          print*, 'sll_qns_2d_with_finite_diff_par: PASSED'
          print*, 'sll_qns2d_angular_spect_method_par: PASSED'
          call flush()
          print*, ' '
       endif
    endif

    call delete_qns_2d_with_finite_diff_plan_par(plan_par)
    call delete_qns_2d_with_finite_diff_plan_seq(plan_seq)
    call delete_qns2d_angular_spect_method_seq(plan_spect_seq)
    call delete_qns2d_angular_spect_method_par(plan_spect_par)

    SLL_DEALLOCATE_ARRAY(phi_par, ierr)
    SLL_DEALLOCATE_ARRAY(c_seq, ierr)
    SLL_DEALLOCATE_ARRAY(c_par, ierr)
    SLL_DEALLOCATE_ARRAY(Te_par, ierr)
    SLL_DEALLOCATE_ARRAY(f, ierr)
    SLL_DEALLOCATE_ARRAY(g, ierr)
 
  end subroutine test_qns_2d


end program qns_tester
