
!***************************************************************************
!
! Selalib 2012     
! Module: sll_qns_2d_with_finite_diff.F90
!
!> @brief 
!> Selalib 2D (r, theta) quasi-neutral solver with finite differences
!> Some arrays are here in 3D for remap utilities
!> Start date: March 12, 2012
!> Last modification: March 23, 2012
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

  use numeric_constants
  use sll_fft
  use sll_tridiagonal
  use sll_collective

  implicit none

  type qns_2d_with_finite_diff_plan_seq
     character(len=100)                      :: bc!Boundary_conditions
     sll_int32                               :: nr!Number of points in r-direction
     sll_int32                               :: ntheta!Number of points in theta-direction
     sll_real64                              :: rmin
     sll_real64                              :: rmax
     sll_real64, dimension(:,:), allocatable :: rho
     sll_real64, dimension(:),   allocatable :: c  
     sll_real64, dimension(:),   allocatable :: Te
     sll_real64, dimension(:),   allocatable :: f
     sll_real64, dimension(:),   allocatable :: g    
     sll_real64                              :: Zi
     type(sll_fft_plan), pointer             :: fft_plan
     type(sll_fft_plan), pointer             :: inv_fft_plan
  end type qns_2d_with_finite_diff_plan_seq

contains


  function new_qns_2d_with_finite_diff_plan_seq(bc, rmin, rmax, rho, c, Te, f, g, Zi) &
                                                                          result (plan)

    character(len=100)                              :: bc ! Boundary_conditions
    sll_real64                                      :: rmin
    sll_real64                                      :: rmax
    sll_real64, dimension(:,:)                      :: rho
    sll_real64, dimension(:)                        :: c, Te, f, g    
    sll_real64                                      :: Zi
    sll_comp64, dimension(:),   allocatable         :: x
    sll_int32                                       :: nr, ntheta, ierr
    type(qns_2d_with_finite_diff_plan_seq), pointer :: plan

    nr = size(rho,1)
    ntheta = size(rho,2)

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE(plan%rho(nr,ntheta), ierr)
    SLL_ALLOCATE(plan%c(nr), ierr)
    SLL_ALLOCATE(plan%Te(nr), ierr)
    SLL_ALLOCATE(plan%f(ntheta), ierr)
    SLL_ALLOCATE(plan%g(ntheta), ierr)
    SLL_ALLOCATE( x(ntheta), ierr )

    plan%bc     = bc
    plan%nr     = nr
    plan%ntheta = ntheta
    plan%rmin   = rmin
    plan%rmax   = rmax
    plan%rho    = rho
    plan%c      = c
    plan%Te     = Te
    plan%f      = f
    plan%g      = g
    plan%Zi     = Zi 

    ! For FFTs in theta-direction
    plan%fft_plan => new_plan_c2c_1d( ntheta, x, x, FFT_FORWARD )

    ! For inverse FFTs in theta-direction
    plan%inv_fft_plan => new_plan_c2c_1d( ntheta, x, x, FFT_INVERSE )

    SLL_DEALLOCATE_ARRAY( x, ierr )

  end function new_qns_2d_with_finite_diff_plan_seq


  subroutine solve_qn_2d_with_finite_diff_seq(plan, phi)

    type(qns_2d_with_finite_diff_plan_seq), pointer :: plan
    sll_real64                                      :: dr, dtheta 
    sll_int32                                       :: nr, ntheta, i, j, ierr
    sll_real64, dimension(:,:)                      :: phi
    sll_comp64, dimension(plan%nr,plan%ntheta)      :: tild_rho, tild_phi
    sll_comp64, dimension(plan%ntheta)              :: f, g
    sll_int32, dimension(plan%nr)                   :: ipiv
    sll_real64, dimension(3*plan%nr)                :: a_resh ! 3*n
    sll_real64, dimension(7*plan%nr)                :: cts  ! 7*n allocation

    nr = plan%nr
    ntheta = plan%ntheta
    if (plan%bc=='neumann') then
       dr = (plan%rmax-plan%rmin)/(nr-1)
    else ! 'dirichlet'
       dr = (plan%rmax-plan%rmin)/(nr+1)
    endif
    dtheta = 2*sll_pi/ntheta

    tild_rho = cmplx(plan%rho, 0_f64, kind=f64)
    f = cmplx(plan%f, 0_f64, kind=f64)
    g = cmplx(plan%g, 0_f64, kind=f64)

    call apply_fft_c2c_1d( plan%fft_plan, f, f )
    call apply_fft_c2c_1d( plan%fft_plan, g, g )
    
    call apply_fft_c2c_1d( plan%fft_plan, tild_rho(1,:), tild_rho(1,:) )
    call apply_fft_c2c_1d( plan%fft_plan, tild_rho(nr,:), tild_rho(nr,:) )

    if (plan%bc=='neumann') then
       tild_rho(1,:) = tild_rho(1,:) + (plan%c(1)-2/dr)*f 
       tild_rho(nr,:) = tild_rho(nr,:) + (plan%c(nr)+2/dr)*g
    else ! 'dirichlet'
       tild_rho(1,:) = tild_rho(1,:) + (1/dr**2 - plan%c(1)/(2*dr))*f
       tild_rho(nr,:) = tild_rho(nr,:) + (1/dr**2 + plan%c(nr)/(2*dr))*g
    endif

    do i=2,nr-1
       call apply_fft_c2c_1d( plan%fft_plan, tild_rho(i,:), tild_rho(i,:) ) 
    enddo

    do j=1,ntheta
       if (plan%bc=='neumann') then
          call neumann_matrix_resh_seq(plan, j-1, a_resh)
       else ! 'dirichlet'
          call dirichlet_matrix_resh_seq(plan, j-1, a_resh)
       endif
       ! b is given by taking the FFT in the theta-direction of rho_{r,theta}      
       ! Solving the linear system: Ax = b  
       call setup_cyclic_tridiag( a_resh, nr, cts, ipiv )
       call solve_cyclic_tridiag( cts, ipiv, tild_rho(:,j), nr, tild_phi(:,j))         
    enddo

    ! Solution phi of the Quasi-neutral equation is given by taking the inverse
    ! FFT in the k-direction of Tild_phi (storaged in phi)  
    do i=1,nr
       call apply_fft_c2c_1d( plan%inv_fft_plan, tild_phi(i,:), tild_phi(i,:) ) 
    enddo

    phi = real(tild_phi, f64)/ntheta

  end subroutine solve_qn_2d_with_finite_diff_seq


  subroutine delete_new_qns_2d_with_finite_diff_plan_seq(plan)

       type (qns_2d_with_finite_diff_plan_seq), pointer :: plan
       sll_int32                                        :: ierr

       ! Fixme: some error checking, whether the poisson pointer is associated
       ! for instance
       SLL_ASSERT( associated(plan) )

       call delete_fft_plan1d(plan%fft_plan)
       call delete_fft_plan1d(plan%inv_fft_plan)

       SLL_DEALLOCATE_ARRAY(plan%rho, ierr)
       SLL_DEALLOCATE_ARRAY(plan%c, ierr)
       SLL_DEALLOCATE_ARRAY(plan%Te, ierr)
       SLL_DEALLOCATE_ARRAY(plan%f, ierr)
       SLL_DEALLOCATE_ARRAY(plan%g, ierr)
       SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_new_qns_2d_with_finite_diff_plan_seq


  subroutine dirichlet_matrix_resh_seq(plan, j, a_resh)

    type(qns_2d_with_finite_diff_plan_seq), pointer :: plan
    sll_real64                                      :: dr, dtheta, Zi
    sll_real64                                      :: r, rmin, rmax
    sll_int32                                       :: i, j, nr, ierr
    sll_real64, dimension(:)                        :: a_resh
    sll_real64, dimension(:), allocatable           :: c, Te 
    ! C & Te are the vector of the Cr & Te(i) respectively

    nr = plan%nr
    rmin = plan%rmin
    rmax = plan%rmax
    dr = (rmax-rmin)/(nr+1)
    dtheta = 2*sll_pi / plan%ntheta
    SLL_ALLOCATE( c(nr), ierr )
    c = plan%c
    SLL_ALLOCATE( Te(nr), ierr )
    Te = plan%Te
    Zi = plan%Zi        

    a_resh = 0.d0

    do i=1,nr
       r = rmin + i*dr
       if (i>1) a_resh(3*(i-1)+1) = -(1/dr**2 - c(i)/(2*dr))
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2 + 1/(Zi*Te(i)) &
                                      - 2*cos(j*dtheta)/(r*dtheta)**2
       if (i<nr) a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    SLL_DEALLOCATE_ARRAY( c, ierr )
    SLL_DEALLOCATE_ARRAY( Te, ierr )

  end subroutine dirichlet_matrix_resh_seq


  subroutine neumann_matrix_resh_seq(plan, j, a_resh)
    type(qns_2d_with_finite_diff_plan_seq), pointer :: plan
    sll_real64                                      :: dr, dtheta, Zi
    sll_real64                                      :: rmin, rmax, r
    sll_int32                                       :: i, j, nr, ierr
    sll_real64, dimension(:)                        :: a_resh
    sll_real64, dimension(:), allocatable           :: c, Te 
    ! c & Te are the vector of the cr & Te(i) respectively

    nr = plan%nr
    rmin = plan%rmin
    rmax = plan%rmax
    dr = (rmax-rmin)/(nr-1)
    dtheta = 2*sll_pi / plan%ntheta
    SLL_ALLOCATE( c(nr), ierr )
    c = plan%c
    SLL_ALLOCATE( Te(nr), ierr )
    Te = plan%Te
    Zi = plan%Zi

    a_resh = 0.d0

    a_resh(2) = 2/dr**2 + 2/(rmin*dtheta)**2 * &
                (-cos(j*dtheta)+1) + 1/(Zi*Te(1))
    a_resh(3) = -2/dr**2

    do i=2,nr-1
       r = rmin + (i-1)*dr
       a_resh(3*(i-1)+1) = -(1/dr**2 - c(i)/(2*dr))
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2 + 1/(Zi*Te(i)) - &
                                          2*cos(j*dtheta)/(r*dtheta)**2
       a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    a_resh(3*(nr-1)+1) = -2/dr**2
    a_resh(3*(nr-1)+2) = 2/dr**2 + 2/(rmax*dtheta)**2 * &
                       (-cos(j*dtheta)+1) + 1/(Zi*Te(nr))

    SLL_DEALLOCATE_ARRAY( c, ierr )
    SLL_DEALLOCATE_ARRAY( Te, ierr )

  end subroutine neumann_matrix_resh_seq


end module sll_qns_2d_with_finite_diff_seq

