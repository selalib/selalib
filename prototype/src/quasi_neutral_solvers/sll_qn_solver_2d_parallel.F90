!> @brief 
!> Selalib 2D (r, theta) quasi-neutral solver 
!>   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************

module sll_qn_solver_2d_parallel

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_assert.h"
#include "sll_constants.h"

  use sll_fft
  use sll_collective
  use sll_remapper
  use sll_boundary_condition_descriptors
  use sll_qn_solver_2d, only : neumann_matrix, dirichlet_matrix

  implicit none

  type qn_solver_2d_parallel
     sll_int32                           :: bc   !< Boundary_conditions
     sll_int32                           :: np_r !< points in r
     sll_int32                           :: np_theta !< points in theta
     sll_real64                          :: rmin
     sll_real64                          :: rmax
     type(sll_fft_plan),         pointer :: fw_fft
     type(sll_fft_plan),         pointer :: bw_fft
     type(layout_2D),            pointer :: layout_fft
     type(layout_2D),            pointer :: layout_lin_sys
     sll_comp64, dimension(:,:), pointer :: fft
     sll_comp64, dimension(:,:), pointer :: lin_sys
     sll_comp64, dimension(:,:), pointer :: c_remap, Te_remap
     type(remap_plan_2D_comp64), pointer :: rmp3_1
     type(remap_plan_2D_comp64), pointer :: rmp3_2
     sll_comp64, dimension(:),   pointer :: hat_f
     sll_comp64, dimension(:),   pointer :: hat_g
     sll_int32, dimension(:),    pointer :: ipiv
     sll_real64, dimension(:),   pointer :: a_resh ! 3*n
     sll_real64, dimension(:),   pointer :: cts ! 7*n allocation
  end type qn_solver_2d_parallel

interface new
module procedure new_qn_solver_2d_parallel
end interface new

interface solve
module procedure solve_qn_solver_2d_parallel
end interface solve

interface delete
module procedure delete_qn_solver_2d_parallel
end interface delete

contains


  !> Allocate new quasi neutral solver parallel
  function new_qn_solver_2d_parallel(bc,rmin,rmax,np_r, np_theta) &
                                                               result (plan)

    type(qn_solver_2d_parallel), pointer    :: plan

    sll_int32                               :: bc !< Boundary_conditions
    sll_real64                              :: rmin
    sll_real64                              :: rmax
    sll_int32                               :: np_r, np_theta
    sll_int32                               :: np_r_loc, np_theta_loc
    sll_int32                               :: ierr
    sll_int64                               :: colsz

    colsz  = sll_get_collective_size(sll_world_collective)
    np_r_loc = np_r/int(colsz)
    np_theta_loc = np_theta/int(colsz)

    if ( colsz > min(np_r,np_theta) ) then     
       print *, 'This test needs to run in a number of processes which',  &
                ' is less than or equal', min(np_r,np_theta)
       print *, 'Exiting...'
       stop
    end if

    if ( (.not.is_power_of_two(int(np_theta,i64))) ) then
       print *, 'The number of points in theta-direction needs to be a power of 2'
       print *, 'Exiting...'
       stop
    end if

    if ( (.not.is_power_of_two(int(np_r,i64))) ) then
       if (bc==SLL_NEUMANN) then
          print *, 'The number of points in r-direction needs to be a power of 2'
          print *, 'Exiting...'
          stop
       else 
          print *, 'The number of points in r-direction needs to be in form 2^n + 2'
          print *, 'Exiting...'
          stop
       endif 
    end if

    SLL_ALLOCATE(plan, ierr)

    plan%bc       = bc
    plan%np_r     = np_r
    plan%np_theta = np_theta
    plan%rmin     = rmin
    plan%rmax     = rmax

    SLL_ALLOCATE(plan%hat_f(np_theta), ierr)
    SLL_ALLOCATE(plan%hat_g(np_theta), ierr)
    SLL_ALLOCATE(plan%ipiv(np_r), ierr)
    SLL_ALLOCATE(plan%a_resh(3*plan%np_r), ierr)
    SLL_ALLOCATE(plan%cts(7*plan%np_r), ierr)
    
    ! For FFTs in theta-direction
    plan%fw_fft => fft_new_plan( np_theta, plan%hat_f, plan%hat_g, FFT_FORWARD )

    ! For inverse FFTs in theta-direction
    plan%bw_fft => fft_new_plan( np_theta, plan%hat_g, plan%hat_f, FFT_INVERSE )

    ! Layout for FFTs-Inv_FFT in theta-direction
    plan%layout_fft => new_layout_2D( sll_world_collective )
    call initialize_layout_with_distributed_2D_array( np_r, np_theta, &
                                       int(colsz), 1, plan%layout_fft )

    ! Layout for Linear systems in r-direction
    plan%layout_lin_sys => new_layout_2D( sll_world_collective )
    call initialize_layout_with_distributed_2D_array( np_r, np_theta, &
                                   1, int(colsz), plan%layout_lin_sys )

    SLL_ALLOCATE(plan%fft(np_r_loc,np_theta), ierr)
    SLL_ALLOCATE(plan%lin_sys(np_r,np_theta_loc), ierr)
    SLL_ALLOCATE(plan%c_remap(np_r,np_theta_loc), ierr)
    SLL_ALLOCATE(plan%Te_remap(np_r,np_theta_loc), ierr)

    plan%rmp3_1 => new_remap_plan( plan%layout_fft, plan%layout_lin_sys, &
                                   plan%fft )
    plan%rmp3_2 => new_remap_plan( plan%layout_lin_sys,plan%layout_fft, &
                                   plan%lin_sys )

  end function new_qn_solver_2d_parallel


 subroutine solve_qn_solver_2d_parallel(plan, rho, c, Te, f, g, Zi, phi)

    type(qn_solver_2d_parallel), pointer :: plan
    sll_real64                           :: dr, dtheta, Zi
    sll_int32                            :: np_r, np_theta
    sll_int32                            :: np_r_loc, np_theta_loc
    sll_int32                            :: i, j
    sll_real64, dimension(:,:)           :: rho, phi
    sll_real64, dimension(:)             :: c, Te, f, g 
    sll_int64                            :: colsz ! collective size
    sll_int32, dimension(2)              :: global
    sll_int32                            :: ind

    colsz = sll_get_collective_size(sll_world_collective)

    np_r         = plan%np_r
    np_theta     = plan%np_theta
    np_r_loc     = np_r/int(colsz)
    np_theta_loc = np_theta/int(colsz)

    if (plan%bc==SLL_NEUMANN) then
       dr = (plan%rmax-plan%rmin)/(np_r-1)
    else 
       dr = (plan%rmax-plan%rmin)/(np_r+1)
    endif
    dtheta = 2*sll_pi/np_theta

    ! FFTs (in theta-direction)

    plan%fft(:,:) = cmplx(rho, 0_f64, kind=f64)
    plan%hat_f = cmplx(f, 0_f64, kind=f64)
    plan%hat_g = cmplx(g, 0_f64, kind=f64)

    call fft_apply_plan( plan%fw_fft, plan%hat_f, plan%hat_f )
    call fft_apply_plan( plan%fw_fft, plan%hat_g, plan%hat_g )

    do i=1,np_r_loc

       call fft_apply_plan( plan%fw_fft, plan%fft(i,:), &
                                         plan%fft(i,:) )
       global = local_to_global_2D( plan%layout_fft, (/i, 1/))
       ind = global(1)
       if (ind==1) then
          if (plan%bc==SLL_NEUMANN) then
             plan%fft(i,:) = plan%fft(i,:) + &
                                             (c(i)-2/dr)*plan%hat_f 
          else 
             plan%fft(i,:) = plan%fft(i,:) + &
                                 (1/dr**2 - c(i)/(2*dr))*plan%hat_f
          endif
       elseif(ind==np_r) then
          if (plan%bc==SLL_NEUMANN) then
             plan%fft(i,:) = plan%fft(i,:) + &
                                             (c(i)+2/dr)*plan%hat_g
          else 
             plan%fft(i,:) = plan%fft(i,:) + &
                                 (1/dr**2 + c(i)/(2*dr))*plan%hat_g
          endif
       endif 
    enddo

    ! Remapping to solve linear systems
    call apply_remap_2D( plan%rmp3_1, plan%fft, plan%lin_sys ) 

    do i=1,np_r_loc
       plan%fft(i,:) = c(i)
    enddo
    call apply_remap_2D( plan%rmp3_1, plan%fft, plan%c_remap )

    do i=1,np_r_loc
       plan%fft(i,:) = Te(i)
    enddo
    call apply_remap_2D( plan%rmp3_1, plan%fft, plan%Te_remap )

    ! Solve linear systems (r-direction)
    do j=1,np_theta_loc
       global = local_to_global_2D( plan%layout_lin_sys, (/1, j/))
       ind = global(2)
       if (ind<=np_theta/2) then
          ind = ind-1
       else
          ind = np_theta-(ind-1)
       endif
       if (plan%bc==SLL_NEUMANN) then

          call neumann_matrix(np_r,                            &
                                   plan%rmin,                       &
                                   plan%rmax,                       &
                                   np_theta,                        &
                                   ind,                             &
                                   real(plan%c_remap(:,1), f64),    &
                                   real(plan%Te_remap(:,1), f64),   &
                                   Zi,                              &
                                   plan%a_resh)

       else 

          call dirichlet_matrix(np_r,                          &
                                     plan%rmin,                     &
                                     plan%rmax,np_theta,            &
                                     ind,                           &
                                     real(plan%c_remap(:,1), f64),  & 
                                     real(plan%Te_remap(:,1), f64), & 
                                     Zi,                            &
                                     plan%a_resh)
       endif 

       call setup_cyclic_tridiag( plan%a_resh, np_r, plan%cts, plan%ipiv )
       call solve_cyclic_tridiag(plan%cts,                &
                                 plan%ipiv,               &
                                 plan%lin_sys(:,j), &
                                 np_r,                    &
                                 plan%lin_sys(:,j))         
    enddo

    ! Remapping to do inverse FFTs
    call apply_remap_2D( plan%rmp3_2, plan%lin_sys, plan%fft ) 

    ! Inverse FFTs (in the theta-direction)
    do i=1,np_r_loc
       call fft_apply_plan( plan%bw_fft, plan%fft(i,:), &
                                         plan%fft(i,:) ) 
    enddo

    phi = real(plan%fft(:,:), f64)/np_theta

  end subroutine solve_qn_solver_2d_parallel

  subroutine delete_qn_solver_2d_parallel(plan)

       type (qn_solver_2d_parallel), pointer :: plan
       sll_int32                                      :: ierr

       ! Fixme: some error checking, whether the poisson pointer is 
       ! associated for instance
       SLL_ASSERT( associated(plan) )

       call fft_delete_plan(plan%fw_fft)
       call fft_delete_plan(plan%bw_fft)

       call delete_layout_2D( plan%layout_fft )
       call delete_layout_2D( plan%layout_lin_sys )

       SLL_DEALLOCATE_ARRAY(plan%fft,ierr)
       SLL_DEALLOCATE_ARRAY(plan%lin_sys,ierr)
       SLL_DEALLOCATE_ARRAY(plan%c_remap,ierr)
       SLL_DEALLOCATE_ARRAY(plan%Te_remap,ierr)

       SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_qn_solver_2d_parallel


end module sll_qn_solver_2d_parallel
