module sll_m_qn_solver_2d_parallel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_neumann

  use sll_m_collective, only: &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_fft, only: &
    sll_s_fft_apply_plan_c2c_1d, &
    sll_p_fft_backward, &
    sll_s_fft_delete_plan, &
    sll_p_fft_forward, &
    sll_f_fft_new_plan_c2c_1d, &
    sll_t_fft_plan

  use sll_m_qn_solver_2d, only: &
    sll_s_dirichlet_matrix, &
    sll_s_neumann_matrix

  use sll_m_remapper, only: &
    sll_o_apply_remap_3d, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_3d, &
    sll_o_local_to_global, &
    sll_f_new_layout_3d, &
    sll_o_new_remap_plan, &
    sll_t_remap_plan_3d_comp64, &
    sll_o_delete

  use sll_m_tridiagonal, only: &
    sll_s_setup_cyclic_tridiag, &
    sll_o_solve_cyclic_tridiag

  use sll_m_utilities, only: &
    sll_f_is_power_of_two

  implicit none

  public :: &
    sll_t_qn_solver_2d_parallel, &
    sll_f_new_qn_solver_2d_parallel, &
    sll_s_solve_qn_solver_2d_parallel, &
    sll_s_delete_qn_solver_2d_parallel

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  type sll_t_qn_solver_2d_parallel
     sll_int32                             :: BC 
     sll_int32                             :: NP_r 
     sll_int32                             :: NP_theta 
     sll_real64                            :: rmin
     sll_real64                            :: rmax
     type(sll_t_fft_plan), pointer           :: fft_plan
     type(sll_t_fft_plan), pointer           :: inv_fft_plan
     type(sll_t_layout_3d),  pointer             :: layout_fft
     type(sll_t_layout_3d),  pointer             :: layout_lin_sys
     sll_comp64, dimension(:,:,:), pointer :: array_fft
     sll_comp64, dimension(:,:,:), pointer :: array_lin_sys
     sll_comp64, dimension(:,:,:), pointer :: c_remap, Te_remap
     type(sll_t_remap_plan_3d_comp64), pointer   :: rmp3_1
     type(sll_t_remap_plan_3d_comp64), pointer   :: rmp3_2
  end type sll_t_qn_solver_2d_parallel

contains


  function sll_f_new_qn_solver_2d_parallel(BC,rmin,rmax,NP_r, NP_theta) &
                                                               result (plan)

    sll_int32                               :: BC 
    sll_real64                              :: rmin
    sll_real64                              :: rmax
    sll_comp64, dimension(:),   allocatable :: x
    sll_int32                               :: NP_r, NP_theta
    sll_int32                               :: NP_r_loc, NP_theta_loc
    sll_int32                               :: ierr
    sll_int32                               :: colsz
    type(sll_t_qn_solver_2d_parallel), pointer    :: plan

    colsz  = sll_f_get_collective_size(sll_v_world_collective)
    NP_r_loc = NP_r/int(colsz)
    NP_theta_loc = NP_theta/int(colsz)

    if ( colsz > min(NP_r,NP_theta) ) then     
       print *, 'This test needs to run in a number of processes which',  &
                ' is less than or equal', min(NP_r,NP_theta)
       print *, 'Exiting...'
       stop
    end if

    if ( (.not.sll_f_is_power_of_two(int(NP_theta,i64))) ) then
       print *, 'The number of points in theta-direction needs to be a power of 2'
       print *, 'Exiting...'
       stop
    end if

    if ( (.not.sll_f_is_power_of_two(int(NP_r,i64))) ) then
       if (BC==sll_p_neumann) then
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
    
    plan%fft_plan => sll_f_fft_new_plan_c2c_1d( NP_theta, x, x, sll_p_fft_forward )

    plan%inv_fft_plan => sll_f_fft_new_plan_c2c_1d( NP_theta, x, x, sll_p_fft_backward )

    plan%layout_fft => sll_f_new_layout_3d( sll_v_world_collective )
    call sll_o_initialize_layout_with_distributed_array( NP_r, NP_theta, 1, &
                                       int(colsz), 1, 1, plan%layout_fft )

    ! Layout for Linear systems in r-direction
    plan%layout_lin_sys => sll_f_new_layout_3d( sll_v_world_collective )
    call sll_o_initialize_layout_with_distributed_array( NP_r, NP_theta, 1, &
                                   1, int(colsz), 1, plan%layout_lin_sys )

    SLL_ALLOCATE(plan%array_fft(NP_r_loc,NP_theta,1), ierr)
    SLL_ALLOCATE(plan%array_lin_sys(NP_r,NP_theta_loc,1), ierr)
    SLL_ALLOCATE(plan%c_remap(NP_r,NP_theta_loc,1), ierr)
    SLL_ALLOCATE(plan%Te_remap(NP_r,NP_theta_loc,1), ierr)

    plan%rmp3_1 => sll_o_new_remap_plan( plan%layout_fft, plan%layout_lin_sys, &
                                   plan%array_fft )
    plan%rmp3_2 => sll_o_new_remap_plan( plan%layout_lin_sys,plan%layout_fft, &
                                   plan%array_lin_sys )

    SLL_DEALLOCATE_ARRAY( x, ierr )

  end function sll_f_new_qn_solver_2d_parallel


 subroutine sll_s_solve_qn_solver_2d_parallel(plan, rho, c, Te, f, g, Zi, phi)

    type(sll_t_qn_solver_2d_parallel), pointer :: plan
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
    sll_int32                            :: colsz ! collective size
    sll_int32, dimension(1:3)            :: global
    sll_int32                            :: ind

    colsz = sll_f_get_collective_size(sll_v_world_collective)

    NP_r         = plan%NP_r
    NP_theta     = plan%NP_theta
    NP_r_loc     = NP_r/int(colsz)
    NP_theta_loc = NP_theta/int(colsz)

    if (plan%BC==sll_p_neumann) then
       dr = (plan%rmax-plan%rmin)/(NP_r-1)
    else ! 'dirichlet'
       dr = (plan%rmax-plan%rmin)/(NP_r+1)
    endif
    dtheta = 2*sll_p_pi/NP_theta

    plan%array_fft(:,:,1) = cmplx(rho, 0_f64, kind=f64)
    hat_f = cmplx(f, 0_f64, kind=f64)
    hat_g = cmplx(g, 0_f64, kind=f64)

    call sll_s_fft_apply_plan_c2c_1d( plan%fft_plan, hat_f, hat_f )
    call sll_s_fft_apply_plan_c2c_1d( plan%fft_plan, hat_g, hat_g )

    do i=1,NP_r_loc

       call sll_s_fft_apply_plan_c2c_1d( plan%fft_plan, plan%array_fft(i,:,1), &
                                              plan%array_fft(i,:,1) )
       global = sll_o_local_to_global( plan%layout_fft, (/i, 1, 1/))
       ind = global(1)
       if (ind==1) then
          if (plan%BC==sll_p_neumann) then
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + &
                                             (c(i)-2/dr)*hat_f 
          else ! 'dirichlet'
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + &
                                 (1/dr**2 - c(i)/(2*dr))*hat_f
          endif
       elseif(ind==NP_r) then
          if (plan%BC==sll_p_neumann) then
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + &
                                             (c(i)+2/dr)*hat_g
          else ! 'dirichlet'
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + &
                                 (1/dr**2 + c(i)/(2*dr))*hat_g
          endif
       endif 
    enddo

    ! Remapping to sll_o_solve linear systems
    call sll_o_apply_remap_3d( plan%rmp3_1, plan%array_fft, plan%array_lin_sys ) 

    do i=1,NP_r_loc
       plan%array_fft(i,:,1) = cmplx(c(i),0.0_f64,f64)
    enddo
    call sll_o_apply_remap_3d( plan%rmp3_1, plan%array_fft, plan%c_remap )

    do i=1,NP_r_loc
       plan%array_fft(i,:,1) = cmplx(Te(i),0.0_f64,f64)
    enddo
    call sll_o_apply_remap_3d( plan%rmp3_1, plan%array_fft, plan%Te_remap )

    ! sll_o_solve linear systems (r-direction)
    do j=1,NP_theta_loc
       global = sll_o_local_to_global( plan%layout_lin_sys, (/1, j, 1/))
       ind = global(2)
       if (ind<=NP_theta/2) then
          ind = ind-1
       else
          ind = NP_theta-(ind-1)
       endif
       if (plan%BC==sll_p_neumann) then
          call sll_s_neumann_matrix(plan%np_r, plan%rmin, plan%rmax, plan%np_theta, &
                              ind, real(plan%c_remap(:,1,1), f64),         &
                              real(plan%Te_remap(:,1,1), f64), Zi, a_resh)
       else ! 'dirichlet'
          call sll_s_dirichlet_matrix(plan%np_r, plan%rmin, plan%rmax, plan%np_theta, &
                                ind, real(plan%c_remap(:,1,1), f64),            &
                                real(plan%Te_remap(:,1,1), f64), Zi, a_resh)
       endif 
       call sll_s_setup_cyclic_tridiag( a_resh, NP_r, cts, ipiv )
       call sll_o_solve_cyclic_tridiag(cts,ipiv,plan%array_lin_sys(:,j,1), &
                                       NP_r,plan%array_lin_sys(:,j,1))         
    enddo

    call sll_o_apply_remap_3d( plan%rmp3_2, plan%array_lin_sys, plan%array_fft ) 

    do i=1,NP_r_loc
       !call sll_s_fft_apply_plan_c2c_1d( plan%inv_fft_plan, plan%array_fft(i,:,1), &
       call sll_s_fft_apply_plan_c2c_1d( plan%inv_fft_plan, plan%array_fft(i,:,1), &
                                                  plan%array_fft(i,:,1) ) 
    enddo

    phi = real(plan%array_fft(:,:,1), f64)/real(NP_theta,f64)

  end subroutine sll_s_solve_qn_solver_2d_parallel

  subroutine sll_s_delete_qn_solver_2d_parallel(plan)

       type (sll_t_qn_solver_2d_parallel), pointer :: plan
       sll_int32                             :: ierr

       SLL_ASSERT( associated(plan) )

       call sll_s_fft_delete_plan(plan%fft_plan)
       call sll_s_fft_delete_plan(plan%inv_fft_plan)
       call sll_o_delete( plan%layout_fft )
       call sll_o_delete( plan%layout_lin_sys )

       SLL_DEALLOCATE_ARRAY(plan%array_fft,ierr)
       SLL_DEALLOCATE_ARRAY(plan%array_lin_sys,ierr)
       SLL_DEALLOCATE_ARRAY(plan%c_remap,ierr)
       SLL_DEALLOCATE_ARRAY(plan%Te_remap,ierr)

       SLL_DEALLOCATE(plan, ierr)

  end subroutine sll_s_delete_qn_solver_2d_parallel


end module sll_m_qn_solver_2d_parallel
