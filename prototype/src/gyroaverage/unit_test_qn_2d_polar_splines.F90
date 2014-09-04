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

program unit_test_qn_2d_polar_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!#include "sll_field_2d.h"
use sll_module_qn_2d_polar_splines_solver

implicit none
  
  class(sll_qn_2d_base), pointer :: qn 
  sll_real64 :: err,tmp,mu_max
  sll_real64 :: eta_min(2)
  sll_real64 :: eta_max(2)
  sll_real64 :: error(3)
  sll_int32  :: Nc(2)
  sll_int32  :: N_min(2),N_max(2)
  sll_int32  :: mode(2)
  sll_int32  :: N_points
  sll_int32  :: N_mu
  sll_int32  :: mu_case
  sll_real64,dimension(:),allocatable :: mu_points,mu_weights,lambda
  sll_real64 :: rho2d(2),val
  sll_real64,dimension(:,:),allocatable :: phi,phi_restr,phi_init,phi_qn
  sll_int32  :: ierr,i
  
  eta_min(1) = 2._f64 ! 0.1_f64
  eta_max(1) = 18._f64 ! 0.9_f64
  eta_min(2) = 0._f64
  eta_max(2) = 2._f64*sll_pi  
  
  mode(1) = 1
  mode(2) = 1
  
  Nc(1)=  64
  Nc(2)=  64
  ! Ntheta doit etre pair
  
  err = 0._f64
  
  mu_case = 2
  ! 1 : Rectangles
  ! 2 : Lin-Lee
  ! 3 : Gauss
  ! 4 : Gauss-Lobatto
  N_mu = 3
  mu_max = 3._f64
  N_points = 8
  
  val = 0._f64
  call compute_gamma0(mode,eta_min,eta_max,val)
  print *,"gamma0 = ",val
  
  SLL_ALLOCATE(phi(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_init(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_qn(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_restr(1:Nc(1)+1,1:Nc(2)),ierr)
  SLL_ALLOCATE(mu_points(1:N_mu),ierr)
  SLL_ALLOCATE(mu_weights(1:N_mu),ierr)
  SLL_ALLOCATE(lambda(1:Nc(1)+1),ierr)
  
  ! mu integration
  if (mu_case==1) then
    if (N_mu==1) then
      mu_points(1) = 1._f64
      mu_weights(1) = 1._f64
    else
      do i=1,N_mu
        mu_points(i) = real(i-1,f64)*mu_max/real(N_mu,f64)
        mu_weights(i) = mu_max/real(N_mu,f64)
      enddo
    endif
  elseif (mu_case==2) then
    if (N_mu==2) then
      mu_points(1) = 0.4167845_f64
      mu_weights(1) = 0.7194_f64
      mu_points(2) = 2.495154605_f64
      mu_weights(2) = 0.2806_f64
    elseif (N_mu==3) then
      mu_points(1) = 0.148131245_f64
      mu_weights(1) = 0.3583_f64
      mu_points(2) = 1._f64
      mu_weights(2) = 0.5004_f64
      mu_points(3) = 3.15959522_f64
      mu_weights(3) = 0.1413_f64
    endif
  endif
  
  print *,"mu_points =",mu_points
  print *,"mu_weights =",mu_weights
  
  call compute_init_f_polar(phi,mode,Nc,eta_min,eta_max)
  phi_init = phi
  phi_restr(1:Nc(1)+1,1:Nc(2)) = phi(1:Nc(1)+1,1:Nc(2))

  lambda = 1._f64

  qn => new_qn_2d_polar_splines_solver( &
    eta_min, &
    eta_max, &
    Nc, &
    N_points, &
    lambda)
  
  call qn%precompute_qn( mu_points, mu_weights , N_mu)
  call qn%solve_qn(phi_restr)

  
  phi_qn(1:Nc(1)+1,1:Nc(2)) = phi_restr(1:Nc(1)+1,1:Nc(2))
  phi_qn(1:Nc(1)+1,Nc(2)+1) = phi_restr(1:Nc(1)+1,1)
  
  call test_solve_qn_polar_new(Nc,eta_min,eta_max,mu_points,mu_weights,N_mu,mode,phi_init,phi_qn)

  if(err==0)then    
    print *, '#PASSED'
  endif

end program