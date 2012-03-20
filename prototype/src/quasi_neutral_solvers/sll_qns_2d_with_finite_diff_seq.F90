
!***************************************************************************
!
! Selalib 2012     
! Module: sll_qns_2d_with_finite_diff.F90
!
!> @brief 
!> Selalib 2D (r, theta) quasi-neutral solver with finite differences
!> Some arrays are here in 3D for remap utilities
!> Start date: March 12, 2012
!> Last modification: March 20, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************

module sll_qns_2d_with_finite_diff_seq

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
#include "sll_remap.h"

  use sll_fft
  use numeric_constants
  use sll_collective
  use sll_tridiagonal
  use sll_qns_2d_with_finite_diff_util

  implicit none

contains

  subroutine solve_qn_2d_with_finite_diff_seq(plan, phi)

    type(qns_2d_with_finite_diff_plan_seq), pointer :: plan
    sll_real64                                      :: dr, dtheta 
    sll_int32                                       :: n, ntheta, i, j, ierr
    sll_real64, dimension(:,:)                      :: phi
    sll_real64, dimension(:,:), allocatable         :: rho
    sll_comp64, dimension(:,:), allocatable         :: tild_rho, tild_phi
    sll_comp64, dimension(:),   allocatable         :: b, f, g
    sll_int32, dimension(:),    allocatable         :: ipiv
    sll_real64, dimension(:),   allocatable         :: a_resh ! 3*n
    sll_real64, dimension(:),   allocatable         :: cts  ! 7*N allocation

    dr = (plan%rmax-plan%rmin)/(plan%nr-1)
    ntheta = plan%ntheta
    dtheta = 2*sll_pi/ntheta

    SLL_ALLOCATE( f(ntheta), ierr )
    SLL_ALLOCATE( g(ntheta), ierr )

    tild_rho = cmplx(rho, 0_f64, kind=f64)
    f = cmplx(plan%f, 0_f64, kind=f64)
    g = cmplx(plan%g, 0_f64, kind=f64)

    call apply_fft_c2c_1d( plan%fft_plan, f, f )
    call apply_fft_c2c_1d( plan%fft_plan, g, g )

    if (plan%bc=='Neumann') then
       n = plan%nr
    else ! 'Dirichlet'
       n = plan%nr - 2
    endif
    
    call apply_fft_c2c_1d( plan%fft_plan, tild_rho(1,:), tild_rho(1,:) )
    call apply_fft_c2c_1d( plan%fft_plan, tild_rho(n,:), tild_rho(n,:) )

    if (plan%bc=='Neumann') then
       tild_rho(1,:) = tild_rho(1,:) + (plan%c(1)-2/dr)*f 
       tild_rho(n,:) = tild_rho(n,:) + (plan%c(n)+2/dr)*g
    else ! 'Dirichlet'
       tild_rho(1,:) = tild_rho(1,:) + (1/dr**2 - plan%c(1)/(2*dr))*f
       tild_rho(n,:) = tild_rho(n,:) + (1/dr**2 + plan%c(n)/(2*dr))*g
    endif

    do i=2,n-1
       call apply_fft_c2c_1d( plan%fft_plan, tild_rho(i,:), tild_rho(i,:) ) 
    enddo

    do j=1,ntheta
       if (plan%bc=='Neumann') then
          call neumann_matrix_resh(plan, j-1, a_resh)
       else ! 'Dirichlet'
          call dirichlet_matrix_resh(plan, j-1, a_resh)
       endif
       b = tild_rho(:,j) ! b is given by taking the FFT in the theta-direction of rho_{r,theta}      
       ! Solving the linear system: Ax = b  
       call setup_cyclic_tridiag( a_resh, n, cts, ipiv )
       call solve_cyclic_tridiag( cts, ipiv, b, n, tild_phi(:,j))         
    enddo

    ! Solution phi of the Quasi-neutral equation is given by taking the inverse FFT in the 
    ! k-direction of Tild_phi (storaged in phi)  
    do i=1,n
       call apply_fft_c2c_1d( plan%inv_fft_plan, tild_phi(i,:), tild_phi(i,:) ) 
    enddo

    phi = real(tild_phi, f64)

  end subroutine solve_qn_2d_with_finite_diff_seq

end module sll_qns_2d_with_finite_diff_seq

