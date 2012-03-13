!***************************************************************************
!
! Selalib 2012     
! Module: sll_qns_2d_with_finite_diff.F90
!
!> @brief 
!> Selalib 2D (r, theta) quasi-neutral solver with finite differences
!> Some arrays are here in 3D for remap utilities
!> Start date: March 13, 2012
!> Last modification: March 13, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************

module sll_qns_2d_with_finite_diff_util

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
#include "sll_remap.h"

  use sll_fft
  use numeric_constants
  use sll_collective
  use sll_tridiagonal

  implicit none

  type qns_2d_with_finite_diff_plan_seq
     character(len=100)                      :: bc ! Boundary_conditions
     sll_int32                               :: n_r ! Number of points in r-direction
     sll_int32                               :: n_theta ! Number of points in theta-direction
     sll_real64                              :: r_min
     sll_real64                              :: r_max
     sll_real64, dimension(:,:), allocatable :: rho
     sll_real64, dimension(:),   allocatable :: c  
     sll_real64, dimension(:),   allocatable :: Te
     sll_real64, dimension(:),   allocatable :: f
     sll_real64, dimension(:),   allocatable :: g    
     sll_real64                              :: Zi
     type(sll_fft_plan), pointer             :: fft_plan
     type(sll_fft_plan), pointer             :: inv_fft_plan
  end type qns_2d_with_finite_diff_plan_seq

  type qns_2d_with_finite_diff_plan_par
     character(len=100)                      :: bc ! Boundary_conditions
     sll_int32                               :: n_r ! Number of points in r-direction
     sll_int32                               :: n_theta ! Number of points in theta-direction
     sll_real64                              :: r_min
     sll_real64                              :: r_max
     sll_real64, dimension(:,:), allocatable :: rho
     sll_real64, dimension(:),   allocatable :: c  
     sll_real64, dimension(:),   allocatable :: Te
     sll_real64, dimension(:),   allocatable :: f
     sll_real64, dimension(:),   allocatable :: g    
     sll_real64                              :: Zi
     type(sll_fft_plan), pointer             :: fft_plan
     type(sll_fft_plan), pointer             :: inv_fft_plan
     type(layout_3D_t),  pointer             :: layout_fft
     type(layout_3D_t),  pointer             :: layout_lin_syst
     sll_int32,               dimension(2,2) :: loc_sizes 
  end type qns_2d_with_finite_diff_plan_par

contains


  function new_qns_2d_with_finite_diff_plan_seq(bc, r_min, r_max, rho, c, Te, f, g, Zi) result (plan)

    character(len=100)                              :: bc ! Boundary_conditions
    sll_real64                                      :: r_min
    sll_real64                                      :: r_max
    sll_real64, dimension(:,:), allocatable         :: rho
    sll_real64, dimension(:),   allocatable         :: c  
    sll_real64, dimension(:),   allocatable         :: Te
    sll_real64, dimension(:),   allocatable         :: f
    sll_real64, dimension(:),   allocatable         :: g    
    sll_real64                                      :: Zi
    sll_comp64, dimension(:),   allocatable         :: x
    sll_int32                                       :: n_r, n_theta, ierr
    sll_int64                                       :: colsz
    type(qns_2d_with_finite_diff_plan_seq), pointer :: plan

    n_r = size(rho,1)
    n_theta = size(rho,2)
    SLL_ALLOCATE(plan, ierr)

    plan%bc      = bc
    plan%n_r     = n_r
    plan%n_theta = n_theta
    plan%r_min   = r_min
    plan%r_max   = r_max
    plan%rho     = rho
    plan%c       = c
    plan%Te      = Te
    plan%f       = f
    plan%g       = g
    plan%Zi      = Zi 

    SLL_ALLOCATE( x(n_theta), ierr )

    ! For FFTs in theta-direction
    plan%fft_plan => new_plan_c2c_1d( n_theta, x, x, FFT_FORWARD )

    ! For inverse FFTs in theta-direction
    plan%inv_fft_plan => new_plan_c2c_1d( n_theta, x, x, FFT_INVERSE )

  end function new_qns_2d_with_finite_diff_plan_seq


  function new_qns_2d_with_finite_diff_plan_par(bc, r_min, r_max, rho, c, Te, f, g, Zi) result (plan)

    character(len=100)                              :: bc ! Boundary_conditions
    sll_real64                                      :: r_min
    sll_real64                                      :: r_max
    sll_real64, dimension(:,:), allocatable         :: rho
    sll_real64, dimension(:),   allocatable         :: c  
    sll_real64, dimension(:),   allocatable         :: Te
    sll_real64, dimension(:),   allocatable         :: f
    sll_real64, dimension(:),   allocatable         :: g    
    sll_real64                                      :: Zi
    sll_comp64, dimension(:),   allocatable         :: x
    sll_int32                                       :: n_r, n_theta, ierr
    sll_int64                                       :: colsz
    type(qns_2d_with_finite_diff_plan_par), pointer :: plan

    n_r = size(rho,1)
    n_theta = size(rho,2)
    SLL_ALLOCATE(plan, ierr)

    plan%bc      = bc
    plan%n_r     = n_r
    plan%n_theta = n_theta
    plan%r_min   = r_min
    plan%r_max   = r_max
    plan%rho     = rho
    plan%c       = c
    plan%Te      = Te
    plan%f       = f
    plan%g       = g
    plan%Zi      = Zi 

    SLL_ALLOCATE( x(n_theta), ierr )

    ! For FFTs in theta-direction
    plan%fft_plan => new_plan_c2c_1d( n_theta, x, x, FFT_FORWARD )

    ! For inverse FFTs in theta-direction
    plan%inv_fft_plan => new_plan_c2c_1d( n_theta, x, x, FFT_INVERSE )

    ! Layout and local sizes for FFTs-Inv_FFT in theta-direction
    plan%layout_fft => new_layout_3D( sll_world_collective )
    colsz  = sll_get_collective_size(sll_world_collective)
    call initialize_layout_with_distributed_3D_array( n_r, n_theta, 1, int(colsz), 1, 1, plan%layout_fft )

    plan%loc_sizes(1,1) = n_r/int(colsz)
    plan%loc_sizes(1,2) = n_theta

    ! Layout and local sizes for Linear systems in r-direction
    plan%layout_lin_syst => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( n_r, n_theta, 1, 1, int(colsz), 1, plan%layout_lin_syst )

    plan%loc_sizes(2,1) = n_r
    plan%loc_sizes(2,2) = n_theta/int(colsz)

    SLL_DEALLOCATE_ARRAY( x, ierr )

  end function new_qns_2d_with_finite_diff_plan_par


  subroutine delete_new_qns_2d_with_finite_diff_plan_seq(plan)

       type (qns_2d_with_finite_diff_plan_seq), pointer :: plan
       sll_int32                                        :: ierr

       ! Fixme: some error checking, whether the poisson pointer is associated
       ! for instance
       SLL_ASSERT( associated(plan) )

       call delete_fft_plan1d(plan%fft_plan)
       call delete_fft_plan1d(plan%inv_fft_plan)
       SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_new_qns_2d_with_finite_diff_plan_seq


  subroutine delete_new_qns_2d_with_finite_diff_plan_par(plan)

       type (qns_2d_with_finite_diff_plan_par), pointer :: plan
       sll_int32                                        :: ierr

       ! Fixme: some error checking, whether the poisson pointer is associated
       ! for instance
       SLL_ASSERT( associated(plan) )

       call delete_fft_plan1d(plan%fft_plan)
       call delete_fft_plan1d(plan%inv_fft_plan)

       call delete_layout_3D( plan%layout_fft )
       call delete_layout_3D( plan%layout_lin_syst )

       SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_new_qns_2d_with_finite_diff_plan_par


  subroutine dirichlet_matrix_resh(plan, j, a_resh)

    type(qns_2d_with_finite_diff_plan_seq), pointer :: plan ! Matrix is sequential
    sll_real64                                  :: dr, dtheta, r, r_min, Zi
    sll_int32                                   :: i, j, n, ierr
    sll_real64, dimension(:)                    :: a_resh
    sll_real64, dimension(:), allocatable       :: c, Te ! C & Te are the vector of the Cr & Te(i) respectively

    n = plan%n_r - 2
    dr = (plan%r_max - plan%r_min) / (n+1)
    dtheta = 2*sll_pi / plan%n_theta
    SLL_ALLOCATE( c(n), ierr )
    c = plan%c
    SLL_ALLOCATE( Te(n), ierr )
    Te = plan%Te
    Zi = plan%Zi       
    r_min = plan%r_min
    a_resh = 0.d0 

    do i=1,n
       r = r_min + i*dr
       if (i>1) then
          a_resh(3*(i-1)+1) = -(1/dr**2 - c(i)/(2*dr))
       endif
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2 + 1/(Zi*Te(i)) - 2*cos((j-1)*dtheta)/(r*dtheta)**2
       if (i<n) then
          a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
       endif
    enddo

    SLL_DEALLOCATE_ARRAY( c, ierr )
    SLL_DEALLOCATE_ARRAY( Te, ierr )

  end subroutine dirichlet_matrix_resh


  subroutine neumann_matrix_resh(plan, j, a_resh)
    type(qns_2d_with_finite_diff_plan_seq), pointer :: plan ! Matrix is sequential
    sll_real64                                      :: dr, dtheta, Zi, rmin, rmax, r
    sll_int32                                       :: i, j, n, ierr
    sll_real64, dimension(:)                        :: a_resh
    sll_real64, dimension(:), allocatable           :: c, Te 
    ! c & Te are the vector of the Cr & Te(i) respectively

    n = plan%n_r
    dr = (plan%r_max - plan%r_min) / (n-1)
    dtheta = 2*sll_pi / plan%n_theta
    SLL_ALLOCATE( c(n), ierr )
    c = plan%c
    SLL_ALLOCATE( Te(n), ierr )
    Te = plan%Te
    Zi = plan%Zi

    a_resh = 0.d0 
    a_resh(2) = 2/dr**2 + 2/(rmin*dtheta)**2 * (-cos((j-1)*dtheta)+1) + 1/(Zi*Te(1))
    a_resh(3) = -2/dr**2

    do i=2,n-1
       r = rmin + (i-1)*dr
       a_resh(3*(i-1)+1) = -(1/dr**2 - c(i)/(2*dr))
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2 + 1/(Zi*Te(i)) &
            - 2*cos((j-1)*dtheta)/(r*dtheta)**2
       a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    a_resh(3*(n-1)+1) = -2/dr**2
    a_resh(3*(n-1)+2) = 2/dr**2 + 2/(rmax*dtheta)**2 * (-cos((j-1)*dtheta)+1) + 1/(Zi*Te(n))

    SLL_DEALLOCATE_ARRAY( c, ierr )
    SLL_DEALLOCATE_ARRAY( Te, ierr )

  end subroutine neumann_matrix_resh


end module sll_qns_2d_with_finite_diff_util
