!> @brief 
!> Selalib 2D (r, theta) quasi-neutral solver with finite differences
!> Start date: March 12, 2012
!> Last modification: May 04, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
module sll_qns2d_with_finite_diff_seq

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_fft
use sll_tridiagonal
use sll_collective
use sll_boundary_condition_descriptors

implicit none

!> QNS solver object 2d
type qns2d_with_finite_diff_plan_seq
   sll_int32                   :: BC       !< boundary_conditions
   sll_int32                   :: NP_r     !< points in r-direction
   sll_int32                   :: NP_theta !< points in theta-direction
   sll_real64                  :: rmin     !< r min
   sll_real64                  :: rmax     !< r max
   type(sll_fft_plan), pointer :: fw_fft   !< forward fft plan
   type(sll_fft_plan), pointer :: bw_fft   !< backward fft plan

   sll_comp64, dimension(:,:), pointer :: hat_rho
   sll_comp64, dimension(:,:), pointer :: hat_phi
   sll_comp64, dimension(:),   pointer :: hat_f
   sll_comp64, dimension(:),   pointer :: hat_g
   sll_int32,  dimension(:),   pointer :: ipiv
   sll_real64, dimension(:),   pointer :: a_resh 
   sll_real64, dimension(:),   pointer :: cts  

end type qns2d_with_finite_diff_plan_seq

sll_int32, private :: i, j

contains

  !> Initialialize a new QNS solver in 2D space
  function new_qns2d_with_finite_diff_plan_seq(BC,     &
                                               rmin,   &
                                               rmax,   &
                                               NP_r,   &
                                               NP_theta) result (plan)

    type(qns2d_with_finite_diff_plan_seq), pointer :: plan !< qns solver

    sll_int32   :: BC       !< Boundary_conditions
    sll_real64  :: rmin     !< r min
    sll_real64  :: rmax     !< r max
    sll_int32   :: NP_r     !< points in r
    sll_int32   :: NP_theta !< points in theta
    sll_int32   :: ierr     !< error code

    SLL_ALLOCATE(plan, ierr)
    SLL_CLEAR_ALLOCATE(plan%hat_rho(1:NP_r,1:NP_theta),ierr) 
    SLL_CLEAR_ALLOCATE(plan%hat_phi(1:NP_r,1:NP_theta),ierr) 

    SLL_CLEAR_ALLOCATE(plan%hat_f(1:NP_theta),ierr)
    SLL_CLEAR_ALLOCATE(plan%hat_g(1:NP_theta),ierr)
    SLL_ALLOCATE(plan%ipiv(1:NP_r),ierr); plan%ipiv = 0
    SLL_CLEAR_ALLOCATE(plan%a_resh(1:3*NP_r),ierr) 
    SLL_CLEAR_ALLOCATE(plan%cts(1:7*NP_r),ierr)

    plan%BC       = BC
    plan%NP_r     = NP_r
    plan%NP_theta = NP_theta
    plan%rmin     = rmin
    plan%rmax     = rmax

    ! For FFTs in theta-direction
    plan%fw_fft => fft_new_plan( NP_theta, plan%hat_f, plan%hat_g, FFT_FORWARD )

    plan%bw_fft => fft_new_plan( NP_theta, plan%hat_g, plan%hat_f, FFT_INVERSE )

  end function new_qns2d_with_finite_diff_plan_seq


  subroutine solve_qns2d_with_finite_diff_seq(plan, rho, c, Te, f, g, Zi, phi)

    type(qns2d_with_finite_diff_plan_seq), pointer :: plan
    sll_real64                                     :: dr
    sll_real64                                     :: dtheta
    sll_real64                                     :: Zi 
    sll_int32                                      :: NP_r
    sll_int32                                      :: NP_theta
    sll_real64, dimension(:,:)                     :: rho, phi
    sll_real64, dimension(:)                       :: c, Te, f, g

    NP_r = plan%NP_r
    NP_theta = plan%NP_theta
    if (plan%BC==SLL_NEUMANN) then
       dr = (plan%rmax-plan%rmin)/(NP_r-1)
    else ! 'dirichlet'
       dr = (plan%rmax-plan%rmin)/(NP_r+1)
    endif
    dtheta = 2*sll_pi/NP_theta

    plan%hat_rho = cmplx(rho, 0_f64, kind=f64)
    plan%hat_f   = cmplx(f,   0_f64, kind=f64)
    plan%hat_g   = cmplx(g,   0_f64, kind=f64)

    call fft_apply_plan( plan%fw_fft, plan%hat_f, plan%hat_f )
    call fft_apply_plan( plan%fw_fft, plan%hat_g, plan%hat_g )
    
    call fft_apply_plan( plan%fw_fft, plan%hat_rho(1,:), plan%hat_rho(1,:) )
    call fft_apply_plan( plan%fw_fft, plan%hat_rho(NP_r,:), plan%hat_rho(NP_r,:) )

    if (plan%bc==SLL_NEUMANN) then
       plan%hat_rho(1,:)  = plan%hat_rho(1,:) + (c(1)-2/dr)*plan%hat_f 
       plan%hat_rho(NP_r,:) = plan%hat_rho(NP_r,:) + (c(NP_r)+2/dr)*plan%hat_g
    else ! dirichlet
       plan%hat_rho(1,:)  = plan%hat_rho(1,:) + (1/dr**2 - c(1)/(2*dr))*plan%hat_f
       plan%hat_rho(NP_r,:) = plan%hat_rho(NP_r,:) + (1/dr**2 + c(NP_r)/(2*dr))*plan%hat_g
    endif

    do i=2,NP_r-1
       call fft_apply_plan( plan%fw_fft, plan%hat_rho(i,:), plan%hat_rho(i,:) ) 
    enddo

    do j=1,NP_theta
       if (plan%BC==SLL_NEUMANN) then
          call neumann_matrix_resh_seq(plan, j-1, c, Te, Zi, plan%a_resh)
       else ! dirichlet
          call dirichlet_matrix_resh_seq(plan, j-1, c, Te, Zi,  plan%a_resh)
       endif
       call setup_cyclic_tridiag( plan%a_resh, NP_r, plan%cts, plan%ipiv )
       call solve_cyclic_tridiag( plan%cts, plan%ipiv, plan%hat_rho(:,j), &
                                  NP_r, plan%hat_phi(:,j))         
    enddo

    ! Solution phi of the Quasi-neutral equation is given by taking the inverse
    ! FFT in the k-direction of Tild_phi (storaged in phi)  
    do i=1,NP_r
       call fft_apply_plan( plan%bw_fft, plan%hat_phi(i,:), plan%hat_phi(i,:) ) 
    enddo

    phi = real(plan%hat_phi, f64)/NP_theta

  end subroutine solve_qns2d_with_finite_diff_seq


  subroutine delete_qns2d_with_finite_diff_plan_seq(plan)

     type (qns2d_with_finite_diff_plan_seq), pointer :: plan
     sll_int32                                       :: ierr

     ! Fixme: some error checking, whether the poisson pointer is associated
     ! for instance
     SLL_ASSERT( associated(plan) )

     call fft_delete_plan(plan%fw_fft)
     call fft_delete_plan(plan%bw_fft)

     SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_qns2d_with_finite_diff_plan_seq


  subroutine dirichlet_matrix_resh_seq(plan, j, c, Te, Zi, a_resh)

    type(qns2d_with_finite_diff_plan_seq), pointer :: plan
    sll_real64, dimension(:)                       :: c, Te
    sll_real64                                     :: dr, dtheta, Zi
    sll_real64                                     :: r, rmin, rmax
    sll_int32                                      :: i, j, NP_r
    sll_real64, dimension(:)                       :: a_resh

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

  end subroutine dirichlet_matrix_resh_seq


  subroutine neumann_matrix_resh_seq(plan, j, c, Te, Zi, a_resh)

    type(qns2d_with_finite_diff_plan_seq), pointer :: plan
    sll_real64, dimension(:)                       :: c, Te
    sll_real64                                     :: dr, dtheta, Zi
    sll_real64                                     :: rmin, rmax, r
    sll_int32                                      :: i, j, NP_r
    sll_real64, dimension(:)                       :: a_resh

    NP_r = plan%NP_r
    rmin = plan%rmin
    rmax = plan%rmax
    dr = (rmax-rmin)/(NP_r-1)
    dtheta = 2*sll_pi / plan%NP_theta

    a_resh = 0.d0

    a_resh(2) = 2/dr**2 + 2/(rmin*dtheta)**2*(1-cos(j*dtheta)) &
                                             + 1/(Zi*Te(1))
    a_resh(3) = -2/dr**2

    do i=2,NP_r-1
       r = rmin + (i-1)*dr
       a_resh(3*(i-1)+1) = c(i)/(2*dr) - 1/dr**2
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2* &
                       (1-cos(j*dtheta)) + 1/(Zi*Te(i))
       a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    a_resh(3*(NP_r-1)+1) = -2/dr**2
    a_resh(3*(NP_r-1)+2) = 2/dr**2 + 2/(rmax*dtheta)**2* &
                       (1-cos(j*dtheta)) + 1/(Zi*Te(NP_r))

  end subroutine neumann_matrix_resh_seq


end module sll_qns2d_with_finite_diff_seq
