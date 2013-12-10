module sll_qn_solver_2d_parallel

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_assert.h"

  use sll_constants
  use sll_fft
  use sll_collective
  use sll_remapper
  use sll_boundary_condition_descriptors
  use sll_qn_solver_2d, only: dirichlet_matrix, neumann_matrix

  implicit none


  type qn_solver_2d_parallel
     sll_int32                             :: BC 
     sll_int32                             :: NP_r 
     sll_int32                             :: NP_theta 
     sll_real64                            :: rmin
     sll_real64                            :: rmax
     type(sll_fft_plan), pointer           :: fft_plan
     type(sll_fft_plan), pointer           :: inv_fft_plan
     type(layout_3D),  pointer             :: layout_fft
     type(layout_3D),  pointer             :: layout_lin_sys
     sll_comp64, dimension(:,:,:), pointer :: array_fft
     sll_comp64, dimension(:,:,:), pointer :: array_lin_sys
     sll_comp64, dimension(:,:,:), pointer :: c_remap, Te_remap
     type(remap_plan_3D_comp64), pointer   :: rmp3_1
     type(remap_plan_3D_comp64), pointer   :: rmp3_2
  end type qn_solver_2d_parallel

contains


  function new_qn_solver_2d_parallel(BC,rmin,rmax,NP_r, NP_theta) &
                                                               result (plan)

    sll_int32                               :: BC 
    sll_real64                              :: rmin
    sll_real64                              :: rmax
    sll_comp64, dimension(:),   allocatable :: x
    sll_int32                               :: NP_r, NP_theta
    sll_int32                               :: NP_r_loc, NP_theta_loc
    sll_int32                               :: ierr
    sll_int64                               :: colsz
    type(qn_solver_2d_parallel), pointer    :: plan

    colsz  = sll_get_collective_size(sll_world_collective)
    NP_r_loc = NP_r/int(colsz)
    NP_theta_loc = NP_theta/int(colsz)

    if ( colsz > min(NP_r,NP_theta) ) then     
       print *, 'This test needs to run in a number of processes which',  &
                ' is less than or equal', min(NP_r,NP_theta)
       print *, 'Exiting...'
       stop
    end if

    if ( (.not.is_power_of_two(int(NP_theta,i64))) ) then
       print *, 'The number of points in theta-direction needs to be a power of 2'
       print *, 'Exiting...'
       stop
    end if

    if ( (.not.is_power_of_two(int(NP_r,i64))) ) then
       if (BC==SLL_NEUMANN) then
          print *, 'The number of points in r-direction needs to be a power of 2'
          print *, 'Exiting...'
          stop
       else ! dirichlet
          print *, 'The number of points in r-direction needs to be in form 2^n + 2'
          print *, 'Exiting...'
          stop
       endif 
    end if

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE( x(NP_theta), ierr )

    plan%BC       = BC
    plan%NP_r     = NP_r
    plan%NP_theta = NP_theta
    plan%rmin     = rmin
    plan%rmax     = rmax
    
    plan%fft_plan => fft_new_plan( NP_theta, x, x, FFT_FORWARD )

    plan%inv_fft_plan => fft_new_plan( NP_theta, x, x, FFT_INVERSE )

    plan%layout_fft => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( NP_r, NP_theta, 1, &
                                       int(colsz), 1, 1, plan%layout_fft )

    ! Layout for Linear systems in r-direction
    plan%layout_lin_sys => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( NP_r, NP_theta, 1, &
                                   1, int(colsz), 1, plan%layout_lin_sys )

    SLL_ALLOCATE(plan%array_fft(NP_r_loc,NP_theta,1), ierr)
    SLL_ALLOCATE(plan%array_lin_sys(NP_r,NP_theta_loc,1), ierr)
    SLL_ALLOCATE(plan%c_remap(NP_r,NP_theta_loc,1), ierr)
    SLL_ALLOCATE(plan%Te_remap(NP_r,NP_theta_loc,1), ierr)

    plan%rmp3_1 => new_remap_plan( plan%layout_fft, plan%layout_lin_sys, &
                                   plan%array_fft )
    plan%rmp3_2 => new_remap_plan( plan%layout_lin_sys,plan%layout_fft, &
                                   plan%array_lin_sys )

    SLL_DEALLOCATE_ARRAY( x, ierr )

  end function new_qn_solver_2d_parallel


 subroutine solve_qn_solver_2d_parallel(plan, rho, c, Te, f, g, Zi, phi)

    type(qn_solver_2d_parallel), pointer :: plan
    sll_real64                           :: dr, dtheta, Zi
    sll_int32                            :: NP_r, NP_theta
    sll_int32                            :: NP_r_loc, NP_theta_loc
    sll_int32                            :: i, j
    sll_real64, dimension(:,:)           :: rho, phi
    sll_real64, dimension(:)             :: c, Te, f, g 
    sll_comp64, dimension(plan%NP_theta) :: hat_f, hat_g
    sll_int32, dimension(plan%NP_r)      :: ipiv
    sll_real64, dimension(3*plan%NP_r)   :: a_resh ! 3*n
    sll_real64, dimension(7*plan%NP_r)   :: cts ! 7*n allocation
    sll_int64                            :: colsz ! collective size
    sll_int32, dimension(1:3)            :: global
    sll_int32                            :: ind

    colsz = sll_get_collective_size(sll_world_collective)

    NP_r         = plan%NP_r
    NP_theta     = plan%NP_theta
    NP_r_loc     = NP_r/int(colsz)
    NP_theta_loc = NP_theta/int(colsz)

    if (plan%BC==SLL_NEUMANN) then
       dr = (plan%rmax-plan%rmin)/(NP_r-1)
    else ! 'dirichlet'
       dr = (plan%rmax-plan%rmin)/(NP_r+1)
    endif
    dtheta = 2*sll_pi/NP_theta

    plan%array_fft(:,:,1) = cmplx(rho, 0_f64, kind=f64)
    hat_f = cmplx(f, 0_f64, kind=f64)
    hat_g = cmplx(g, 0_f64, kind=f64)

    call fft_apply_plan( plan%fft_plan, hat_f, hat_f )
    call fft_apply_plan( plan%fft_plan, hat_g, hat_g )

    do i=1,NP_r_loc

       call fft_apply_plan( plan%fft_plan, plan%array_fft(i,:,1), &
                                              plan%array_fft(i,:,1) )
       global = local_to_global_3D( plan%layout_fft, (/i, 1, 1/))
       ind = global(1)
       if (ind==1) then
          if (plan%BC==SLL_NEUMANN) then
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + &
                                             (c(i)-2/dr)*hat_f 
          else ! 'dirichlet'
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + &
                                 (1/dr**2 - c(i)/(2*dr))*hat_f
          endif
       elseif(ind==NP_r) then
          if (plan%BC==SLL_NEUMANN) then
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + &
                                             (c(i)+2/dr)*hat_g
          else ! 'dirichlet'
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + &
                                 (1/dr**2 + c(i)/(2*dr))*hat_g
          endif
       endif 
    enddo

    ! Remapping to solve linear systems
    call apply_remap_3D( plan%rmp3_1, plan%array_fft, plan%array_lin_sys ) 

    do i=1,NP_r_loc
       plan%array_fft(i,:,1) = c(i)
    enddo
    call apply_remap_3D( plan%rmp3_1, plan%array_fft, plan%c_remap )

    do i=1,NP_r_loc
       plan%array_fft(i,:,1) = Te(i)
    enddo
    call apply_remap_3D( plan%rmp3_1, plan%array_fft, plan%Te_remap )

    ! Solve linear systems (r-direction)
    do j=1,NP_theta_loc
       global = local_to_global_3D( plan%layout_lin_sys, (/1, j, 1/))
       ind = global(2)
       if (ind<=NP_theta/2) then
          ind = ind-1
       else
          ind = NP_theta-(ind-1)
       endif
       if (plan%BC==SLL_NEUMANN) then
          call neumann_matrix(plan%np_r, plan%rmin, plan%rmax, plan%np_theta, &
                              ind, real(plan%c_remap(:,1,1), f64),         &
                              real(plan%Te_remap(:,1,1), f64), Zi, a_resh)
       else ! 'dirichlet'
          call dirichlet_matrix(plan%np_r, plan%rmin, plan%rmax, plan%np_theta, &
                                ind, real(plan%c_remap(:,1,1), f64),            &
                                real(plan%Te_remap(:,1,1), f64), Zi, a_resh)
       endif 
       call setup_cyclic_tridiag( a_resh, NP_r, cts, ipiv )
       call solve_cyclic_tridiag(cts,ipiv,plan%array_lin_sys(:,j,1), &
                                       NP_r,plan%array_lin_sys(:,j,1))         
    enddo

    call apply_remap_3D( plan%rmp3_2, plan%array_lin_sys, plan%array_fft ) 

    do i=1,NP_r_loc
       !call fft_apply_plan_c2c_1d( plan%inv_fft_plan, plan%array_fft(i,:,1), &
       call fft_apply_plan( plan%inv_fft_plan, plan%array_fft(i,:,1), &
                                                  plan%array_fft(i,:,1) ) 
    enddo

    phi = real(plan%array_fft(:,:,1), f64)/NP_theta

  end subroutine solve_qn_solver_2d_parallel

  subroutine delete_qn_solver_2d_parallel(plan)

       type (qn_solver_2d_parallel), pointer :: plan
       sll_int32                             :: ierr

       SLL_ASSERT( associated(plan) )

       call fft_delete_plan(plan%fft_plan)
       call fft_delete_plan(plan%inv_fft_plan)
       call delete_layout_3D( plan%layout_fft )
       call delete_layout_3D( plan%layout_lin_sys )

       SLL_DEALLOCATE_ARRAY(plan%array_fft,ierr)
       SLL_DEALLOCATE_ARRAY(plan%array_lin_sys,ierr)
       SLL_DEALLOCATE_ARRAY(plan%c_remap,ierr)
       SLL_DEALLOCATE_ARRAY(plan%Te_remap,ierr)

       SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_qn_solver_2d_parallel


end module sll_qn_solver_2d_parallel
