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
     
     sll_real64          :: eta_min(2)     !< r min et theta min
     sll_real64          :: eta_max(2)     !< r max et theta max
     sll_int32           :: Nc(2)          !< number of cells in r and in theta
     
     sll_int32           :: N_points          !< number of points on the circle
     sll_int32           :: interp_degree(2)  !< interpolation degrees in r,theta

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

contains


  function new_plan_gyroaverage_polar_hermite(eta_min,eta_max,Nc,N_points,interp_degree,deriv_size) result(this)

    implicit none

    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, intent(in)  :: N_points  
    sll_int32, intent(in)  :: interp_degree(2)
    sll_int32, intent(in)  :: deriv_size
    type(sll_plan_gyroaverage_polar), pointer :: this

    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(this%deriv(deriv_size,Nc(1)+1,Nc(2)+1),err)
    SLL_ALLOCATE(this%points(3,N_points),err)
       
    this%eta_min=eta_min
    this%eta_max=eta_max
    this%Nc=Nc
    this%N_points=N_points
    this%interp_degree=interp_degree
    
  end function new_plan_gyroaverage_polar_hermite
  
  
  
  
  
  
  
  

!subroutine compute_gyroaverage_pade_polar(gyro,f,rho)
!    type(sll_plan_gyroaverage_polar_pade)  :: gyro
!    sll_real64,dimension(:,:),allocatable :: fcomp
!    sll_real64,dimension(:,:),intent(inout) :: f
!    sll_real64,dimension(:),allocatable :: buf,diagm1,diag,diagp1
!    sll_real64,intent(in)::rho
!    sll_int32 :: i,j,k
!    sll_real64::sum_fval,eta_star(2),eta(2),x(2),dr
!    sll_int32 :: ierr
!    
!    dr=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
!    
!	SLL_ALLOCATE(buf(2*gyro%Nc(2)+15),ierr)
!	SLL_ALLOCATE(fcomp(1:gyro%Nc(1)+1,1:gyro%Nc(2)),ierr)
!	SLL_ALLOCATE(diagm1(1:gyro%Nc(1)+1),ierr)
!	SLL_ALLOCATE(diag(1:gyro%Nc(1)+1),ierr)
!	SLL_ALLOCATE(diagp1(1:gyro%Nc(1)+1),ierr)
!
!	fcomp(1:gyro%Nc(1)+1,1:gyro%Nc(2))=f(1:gyro%Nc(1)+1,1:gyro%Nc(2))
!
!    !*** Perform FFT 1D in theta direction of ***
!    !***   the system solution                ***
!	call dffti(gyro%Nc(2),buf)
!	do i=1,gyro%Nc(1)+1
!		call dfftf(gyro%Nc(2),fcomp(i,:),buf)
!	enddo
!	fcomp=fcomp/real(gyro%Nc(2),f64)
!
! 	!***POISSON
!	do k=1,gyro%Nc(2)
!	  do i=1,gyro%Nc(1)
!	    diagm1(i+1)=-(rho(1)**2/4)*(1/dr**2-1/(2*dr*(gyro%eta_min(1)+ &
!	    	(gyro%eta_max(1)-gyro%eta_min(1))*real(i,f64)/real(gyro%Nc(1),f64))))
!	    diag(i)=1-(rho(1)**2/4)*(-(2/dr**2)-((floor(k/2._f64)*1._f64)/ &
!	    	(gyro%eta_min(1)+(gyro%eta_max(1)-gyro%eta_min(1))*real(i-1,f64)/real(gyro%Nc(1),f64)))**2)
!	    diagp1(i)=-(rho(1)**2/4)*(1/dr**2+1/(2*dr*(gyro%eta_min(1)+ &
!	    	(gyro%eta_max(1)-gyro%eta_min(1))*real(i-1,f64)/real(gyro%Nc(1),f64))))
!	  enddo
!	  diagm1(1)=0._f64
!	  diagp1(gyro%Nc(1)+1)=0._f64
!	  diag(1)=1._f64
!	  diag(gyro%Nc(1)+1)=1._f64
!	  !***  Dirichlet boundary conditions ***	  
!	  diagp1(1)=0._f64
!	  diagm1(gyro%Nc(1)+1)=0._f64
!	  !***  Neumann boundary conditions ***
!!	  diagp1(1)=-1._f64
!!	  diagm1(gyro%Nc(1)+1)=-1._f64
!	  call solve_tridiag(diagm1,diag,diagp1,fcomp(1:gyro%Nc(1)+1,k),f(1:gyro%Nc(1)+1,k),gyro%Nc(1)+1)
!	enddo
!
!	!*** Perform FFT 1D inverse ***
!	do i=1,gyro%Nc(1)+1
!	  call dfftb(gyro%Nc(2),f(i+1,1:gyro%Nc(2)),buf)
!	enddo
!         
!    !*** duplicate periodic value ***
!    f(1:gyro%Nc(1)+1,gyro%Nc(2)+1)=f(1:gyro%Nc(1)+1,1)
!    
!  end subroutine compute_gyroaverage_pade_polar


end module sll_gyroaverage_2d_polar
