!***************************************************************************
!
! Selalib 2012     
! Module: sll_qns_2d_with_finite_diff.F90
!
!> @brief 
!> Selalib 2D (r, theta) quasi-neutral solver with finite differences
!> Some arrays are here in 3D for remap utilities
!> Start date: March 13, 2012
!> Last modification: May 03, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************

module sll_qns_2d_with_finite_diff_par

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
#include "sll_remap.h"

  use numeric_constants
  use sll_fft
  use sll_tridiagonal
  use sll_collective

  implicit none


  type qns_2d_with_finite_diff_plan_par
     character(len=100)                        :: BC ! Boundary_conditions
     sll_int32                                 :: NP_r ! Number of points in r-direction
     sll_int32                                 :: NP_theta!Number of points in theta-direction
     sll_real64                                :: rmin
     sll_real64                                :: rmax
     type(sll_fft_plan), pointer               :: fft_plan
     type(sll_fft_plan), pointer               :: inv_fft_plan
     type(layout_3D_t),  pointer               :: layout_fft
     type(layout_3D_t),  pointer               :: layout_lin_sys
     sll_comp64, dimension(:,:,:), allocatable :: array_fft
     sll_comp64, dimension(:,:,:), allocatable :: array_lin_sys
     sll_comp64, dimension(:,:,:), allocatable :: c_remap, Te_remap
     type(remap_plan_3D_t), pointer            :: rmp3_1!remap plan for fft to linear sytem (for rho)
     type(remap_plan_3D_t), pointer            :: rmp3_2!remap plan for linear sytem to inverse fft (for phi)
  end type qns_2d_with_finite_diff_plan_par

contains


  function new_qns_2d_with_finite_diff_plan_par(BC,rmin,rmax,NP_r, NP_theta) result (plan)

    character(len=100)                              :: BC ! boundary_conditions
    sll_real64                                      :: rmin ! minimum radius
    sll_real64                                      :: rmax ! maximum radius
    sll_comp64, dimension(:),   allocatable         :: x
    sll_int32                                       :: NP_r, NP_theta
    ! NP_r and NP_theta are the numbers of points in directions r and theta respectively
    sll_int32                                       :: NP_r_loc, NP_theta_loc
    ! NP_r_loc and NP_theta_loc are the numbers of points locally (in the processor) 
    ! in directions r and theta respectively
    sll_int32                                       :: ierr
    sll_int64                                       :: colsz
    type(qns_2d_with_finite_diff_plan_par), pointer :: plan

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
       if (bc=='neumann') then
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
    SLL_ALLOCATE( x(NP_theta_loc), ierr )

    plan%bc     = bc
    plan%NP_r   = NP_r
    plan%NP_theta = NP_theta
    plan%rmin   = rmin
    plan%rmax   = rmax

    ! For FFTs in theta-direction
    plan%fft_plan => new_plan_c2c_1d( NP_theta, x, x, FFT_FORWARD )

    ! For inverse FFTs in theta-direction
    plan%inv_fft_plan => new_plan_c2c_1d( NP_theta, x, x, FFT_INVERSE )

    ! Layout for FFTs-Inv_FFT in theta-direction
    plan%layout_fft => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( NP_r, NP_theta, &
                                 1, int(colsz), 1, 1, plan%layout_fft )

    ! Layout for Linear systems in r-direction
    plan%layout_lin_sys => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( NP_r, NP_theta, &
                             1, 1, int(colsz), 1, plan%layout_lin_sys )

    SLL_ALLOCATE(plan%array_fft(NP_r_loc,NP_theta,1), ierr)
    SLL_ALLOCATE(plan%array_lin_sys(NP_r,NP_theta_loc,1), ierr)
    SLL_ALLOCATE(plan%c_remap(NP_r,NP_theta_loc,1), ierr)
    SLL_ALLOCATE(plan%Te_remap(NP_r,NP_theta_loc,1), ierr)

    plan%rmp3_1 => NEW_REMAPPER_PLAN_3D(plan%layout_fft, plan%layout_lin_sys,plan%array_fft)
    plan%rmp3_2 => NEW_REMAPPER_PLAN_3D(plan%layout_lin_sys,plan%layout_fft,plan%array_lin_sys)

    SLL_DEALLOCATE_ARRAY( x, ierr )

  end function new_qns_2d_with_finite_diff_plan_par


 subroutine solve_qn_2d_with_finite_diff_par(plan, rho, c, Te, f, g, Zi, phi)

    type(qns_2d_with_finite_diff_plan_par), pointer :: plan
    sll_real64                                      :: dr, dtheta, Zi
    sll_int32                                       :: NP_r, NP_theta
    ! NP_r and NP_theta are the numbers of points in directions r and theta respectively
    sll_int32                                       :: NP_r_loc, NP_theta_loc
    ! NP_r_loc and NP_theta_loc are the numbers of points locally (in the processor) 
    ! in directions r and theta respectively
    sll_int32                                       :: i, j
    sll_real64, dimension(:,:)                      :: rho, phi
    sll_real64, dimension(:)                        :: c, Te, f, g 
    sll_comp64, dimension(plan%NP_theta)            :: hat_f, hat_g
    sll_int32, dimension(plan%NP_r)                 :: ipiv
    sll_real64, dimension(3*plan%NP_r)              :: a_resh ! 3*n
    sll_real64, dimension(7*plan%NP_r)              :: cts ! 7*n allocation
    sll_int64                                       :: colsz ! collective size
    sll_int32, dimension(1:3)                       :: global
    sll_int32                                       :: ind

    colsz = sll_get_collective_size(sll_world_collective)

    NP_r         = plan%NP_r
    NP_theta     = plan%NP_theta
    NP_r_loc     = NP_r/int(colsz)
    NP_theta_loc = NP_theta/int(colsz)

    if (plan%BC=='neumann') then
       dr = (plan%rmax-plan%rmin)/(NP_r-1)
    else ! 'dirichlet'
       dr = (plan%rmax-plan%rmin)/(NP_r+1)
    endif
    dtheta = 2*sll_pi/NP_theta

    ! FFTs (in theta-direction)

    plan%array_fft(:,:,1) = cmplx(rho, 0_f64, kind=f64)
    hat_f = cmplx(f, 0_f64, kind=f64)
    hat_g = cmplx(g, 0_f64, kind=f64)

    call apply_fft_c2c_1d( plan%fft_plan, hat_f, hat_f )
    call apply_fft_c2c_1d( plan%fft_plan, hat_g, hat_g )

    do i=1,NP_r_loc
       call apply_fft_c2c_1d( plan%fft_plan, plan%array_fft(i,:,1), plan%array_fft(i,:,1) )
       global = local_to_global_3D( plan%layout_fft, (/i, 1, 1/))
       ind = global(1)
       if (ind==1) then
          if (plan%bc=='neumann') then
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + (c(i)-2/dr)*hat_f 
          else ! 'dirichlet'
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + (1/dr**2 - c(i)/(2*dr))*hat_f
          endif
       elseif(ind==NP_r) then
          if (plan%bc=='neumann') then
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + (c(i)+2/dr)*hat_g
          else ! 'dirichlet'
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + (1/dr**2 + c(i)/(2*dr))*hat_g
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
       if (plan%BC=='neumann') then
          call neumann_matrix_resh_par(plan, ind-1, real(plan%c_remap(:,1,1), f64), & 
                                         real(plan%Te_remap(:,1,1), f64), Zi, a_resh)
       else ! 'dirichlet'
          call dirichlet_matrix_resh_par(plan, ind-1, real(plan%c_remap(:,1,1), f64), &
                                           real(plan%Te_remap(:,1,1), f64), Zi, a_resh)
       endif 
       call setup_cyclic_tridiag( a_resh, NP_r, cts, ipiv )
       call solve_cyclic_tridiag(cts, ipiv, plan%array_lin_sys(:,j,1), &
                                         NP_r,plan%array_lin_sys(:,j,1))         
    enddo

    ! Remapping to do inverse FFTs
    call apply_remap_3D( plan%rmp3_2, plan%array_lin_sys, plan%array_fft ) 

    ! Inverse FFTs (in the theta-direction)
    do i=1,NP_r_loc
       call apply_fft_c2c_1d( plan%inv_fft_plan, plan%array_fft(i,:,1), &
                                                  plan%array_fft(i,:,1) ) 
    enddo

    phi = real(plan%array_fft(:,:,1), f64)/NP_theta

  end subroutine solve_qn_2d_with_finite_diff_par


  subroutine delete_qns_2d_with_finite_diff_plan_par(plan)

       type (qns_2d_with_finite_diff_plan_par), pointer :: plan
       sll_int32                                        :: ierr

       ! Fixme: some error checking, whether the poisson pointer is 
       ! associated for instance
       SLL_ASSERT( associated(plan) )

       call delete_fft_plan1d(plan%fft_plan)
       call delete_fft_plan1d(plan%inv_fft_plan)

       call delete_layout_3D( plan%layout_fft )
       call delete_layout_3D( plan%layout_lin_sys )

       SLL_DEALLOCATE_ARRAY(plan%array_fft,ierr)
       SLL_DEALLOCATE_ARRAY(plan%array_lin_sys,ierr)
       SLL_DEALLOCATE_ARRAY(plan%c_remap,ierr)
       SLL_DEALLOCATE_ARRAY(plan%Te_remap,ierr)

       SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_qns_2d_with_finite_diff_plan_par


  subroutine dirichlet_matrix_resh_par(plan, j, c, Te, Zi, a_resh)

    type(qns_2d_with_finite_diff_plan_par), pointer :: plan
    sll_real64, dimension(:)                        :: c, Te
    sll_real64                                      :: dr, dtheta, Zi
    sll_real64                                      :: r, rmin, rmax
    sll_int32                                       :: i, j, NP_r
    ! NP_r is the number of points in directions r
    sll_real64, dimension(:)                        :: a_resh 

    NP_r = plan%NP_r
    rmin = plan%rmin
    rmax = plan%rmax
    dr = (rmax-rmin)/(NP_r+1)
    dtheta = 2*sll_pi / plan%NP_theta

    a_resh = 0.d0

    do i=1,NP_r
       r = rmin + i*dr
       if (i>1) then
           a_resh(3*(i-1)+1) = c(i)/(2*dr) - 1/dr**2
       endif
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2*(1-cos(j*dtheta)) &
                                                     + 1/(Zi*Te(i))
       if (i<NP_r) then
           a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
       endif
    enddo

  end subroutine dirichlet_matrix_resh_par


  subroutine neumann_matrix_resh_par(plan, j, c, Te, Zi, a_resh)

    type(qns_2d_with_finite_diff_plan_par), pointer :: plan
    sll_real64, dimension(:)                        :: c, Te
    sll_real64                                      :: dr, dtheta, Zi
    sll_real64                                      :: rmin, rmax, r
    sll_int32                                       :: i, j, NP_r
    ! NP_r is the number of points in directions r
    sll_real64, dimension(:)                        :: a_resh

    NP_r = plan%NP_r
    rmin = plan%rmin
    rmax = plan%rmax
    dr = (rmax-rmin)/(NP_r-1)
    dtheta = 2*sll_pi / plan%NP_theta

    a_resh = 0.d0

    a_resh(2) = 2/dr**2 + 2/(rmin*dtheta)**2*(1-cos(j*dtheta)) + &
                                                 1/(Zi*Te(1))
    a_resh(3) = -2/dr**2

    do i=2,NP_r-1
       r = rmin + (i-1)*dr
       a_resh(3*(i-1)+1) = c(i)/(2*dr) - 1/dr**2
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2*(1-cos(j*dtheta)) &
                                                          + 1/(Zi*Te(i))
       a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    a_resh(3*(NP_r-1)+1) = -2/dr**2
    a_resh(3*(NP_r-1)+2) = 2/dr**2 + 2/(rmax*dtheta)**2*(1-cos(j*dtheta)) &
                                                          + 1/(Zi*Te(NP_r))

  end subroutine neumann_matrix_resh_par


end module sll_qns_2d_with_finite_diff_par
