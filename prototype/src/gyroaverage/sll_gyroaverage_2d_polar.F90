!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************


module sll_gyroaverage_2d_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_fft
  use sll_tridiagonal
  use sll_constants
  use sll_boundary_condition_descriptors

  implicit none
  type sll_plan_gyroaverage_polar
    
     !> Domain 
     
     sll_real64          :: eta_min(2)     !< r min et theta min
     sll_real64          :: eta_max(2)     !< r max et theta max
     sll_int32           :: Nc(2)          !< number of cells in r and in theta
     
     !> Method
     
     sll_int32           :: N_points          !< number of points on the circle
     sll_int32           ::	interp_degree(2)  !< degree of interpolation in r and theta

     sll_real64, dimension(:,:), pointer    :: points
     sll_real64, dimension(:,:,:), pointer  :: deriv
     sll_int32, dimension(:), pointer       :: pre_compute_N
     sll_real64, dimension(:,:), pointer    :: pre_compute_coeff
     sll_int32, dimension(:,:), pointer     :: pre_compute_index
     sll_real64, dimension(:), pointer      :: pre_compute_coeff_spl
     sll_real64, dimension(:,:), pointer    :: lunat
     sll_real64, dimension(:), pointer      :: luper
     sll_real64, dimension(:,:,:), pointer  :: A_fft

  end type sll_plan_gyroaverage_polar
  
  !> SLL_GYROAVERAGE_PADE
  
! ... 
    
  !> SLL_GYROAVERAGE_HERMITE

! subroutine compute_gyroaverage_points_polar_hermite(gyro,f,rho)
!    type(sll_plan_gyroaverage_polar)  :: gyro
!    sll_real64,dimension(:,:),intent(inout) :: f
!    sll_real64,intent(in)::rho
!    sll_int32::i,j,k,ii(2),s
!    sll_real64::fval,sum_fval,eta_star(2),eta(2),delta_eta(2),x(2)
!    
!    fval=0._f64
!    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%N(1),f64)
!    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%N(2),f64)
!    
!    call hermite_coef_nat_per(f(1:gyro%N(1)+1,1:gyro%N(2)),gyro%deriv,gyro%N,gyro%interp_degree)
!    
!    do j=1,gyro%N(2)
!     eta(2)=gyro%eta_min(2)+real(j-1,f64)*delta_eta(2)
!      do i=1,gyro%N(1)+1
!        eta(1)=gyro%eta_min(1)+real(i-1,f64)*delta_eta(1)       
!        sum_fval = 0._f64
!        do k=1,gyro%N_points
!          x(1) = eta(1)*cos(eta(2))+rho*gyro%points(1,k)
!          x(2) = eta(1)*sin(eta(2))+rho*gyro%points(2,k)
!          call localize_polar(x,gyro%eta_min,gyro%eta_max,ii,eta_star,gyro%N)
!          call interpolate_hermite(gyro%deriv,ii,eta_star,fval,gyro%N)
!          sum_fval = sum_fval+gyro%points(3,k)*fval
!        enddo
!        f(i,j) = sum_fval
!      enddo
!    enddo
!    f(1:gyro%N(1)+1,gyro%N(2)+1)=f(1:gyro%N(1)+1,1)
!    
!  end subroutine compute_gyroaverage_points_polar_hermite

  !> SLL_GYROAVERAGE_HERMITE_C1
  
! ...    

  !> SLL_GYROAVERAGE_HERMITE_C1_PRECOMPUTE
  
! ...    
  
  !> SLL_GYROAVERAGE_HERMITE_C1_WITH_INVARIANCE
  
! ...    

  !> SLL_GYROAVERAGE_SPLINES
  
! ...    

  !> SLL_GYROAVERAGE_SPLINES_PRECOMPUTE
  
! ...    

  !> SLL_GYROAVERAGE_SPLINES_WITH_INVARIANCE
  
! ...  

  !> SLL_GYROAVERAGE_SPLINES_PRECOMPUTE_WITH_FFT
  
! ...  

end module sll_gyroaverage_2d_polar
