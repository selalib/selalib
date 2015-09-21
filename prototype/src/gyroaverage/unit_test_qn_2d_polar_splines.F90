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
  
  Nc(1)=  32
  Nc(2)=  32
  ! Ntheta doit etre pair
  
  err = 0._f64
  
  mu_case = 1
  ! 1 : Rectangles
  ! 2 : Lin-Lee
  ! 3 : Gauss-Laguerre

  N_mu = 1
  mu_max = 0.01_f64
  N_points = 32
  
  val = 0._f64
  call compute_gamma0(mode,eta_min,eta_max,val)
  print *,"gamma0 = ",val
  
  SLL_ALLOCATE(phi(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_init(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_qn(1:Nc(1)+1,1:Nc(2)+1),ierr)
  SLL_ALLOCATE(phi_restr(1:Nc(1)+1,1:Nc(2)),ierr)
  SLL_ALLOCATE(mu_points(0:N_mu-1),ierr)
  SLL_ALLOCATE(mu_weights(0:N_mu-1),ierr)
  SLL_ALLOCATE(lambda(1:Nc(1)+1),ierr)
  
  ! mu integration
  if (mu_case==1) then
    if (N_mu==1) then
      mu_points(0) = 1._f64
      mu_weights(0) = 1._f64
    else
      do i=1,N_mu
        mu_points(i-1) = real(i-1,f64)*mu_max/real(N_mu,f64)
        mu_weights(i-1) = mu_max/real(N_mu,f64)
      enddo
    endif
  elseif (mu_case==2) then
    if (N_mu==2) then
      mu_points(0) = 0.4167845_f64
      mu_weights(0) = 0.7194_f64
      mu_points(1) = 2.495154605_f64
      mu_weights(1) = 0.2806_f64
    elseif (N_mu==3) then
      mu_points(0) = 0.148131245_f64
      mu_weights(0) = 0.3583_f64
      mu_points(1) = 1._f64
      mu_weights(1) = 0.5004_f64
      mu_points(2) = 3.15959522_f64
      mu_weights(2) = 0.1413_f64
    endif
  elseif (mu_case==3) then

    mu_points( 0)= 111.751398098 ;mu_weights( 0)= 8.22732055833e-33 ;
    mu_points( 1)= 0.0444893658333 ;mu_weights( 1)= 0.109218341952 ;
    mu_points( 2)= 0.23452610952 ;mu_weights( 2)= 0.210443107939 ;
    mu_points( 3)= 0.576884629302 ;mu_weights( 3)= 0.23521322967 ;
    mu_points( 4)= 1.07244875382 ;mu_weights( 4)= 0.195903335973 ;
    mu_points( 5)= 1.72240877644 ;mu_weights( 5)= 0.129983786286 ;
    mu_points( 6)= 2.52833670643 ;mu_weights( 6)= 0.0705786238657 ;
    mu_points( 7)= 3.49221327302 ;mu_weights( 7)= 0.0317609125092 ;
    mu_points( 8)= 4.61645676975 ;mu_weights( 8)= 0.0119182148348 ;
    mu_points( 9)= 5.90395850417 ;mu_weights( 9)= 0.00373881629461 ;
    mu_points( 10)= 7.35812673319 ;mu_weights( 10)= 0.000980803306615 ;
    mu_points( 11)= 8.98294092421 ;mu_weights( 11)= 0.000214864918801 ;
    mu_points( 12)= 10.7830186325 ;mu_weights( 12)= 3.92034196799e-05 ;
    mu_points( 13)= 12.7636979867 ;mu_weights( 13)= 5.93454161287e-06 ;
    mu_points( 14)= 14.9311397555 ;mu_weights( 14)= 7.41640457867e-07 ;
    mu_points( 15)= 17.2924543367 ;mu_weights( 15)= 7.60456787912e-08 ;
    mu_points( 16)= 19.8558609403 ;mu_weights( 16)= 6.35060222663e-09 ;
    mu_points( 17)= 22.6308890132 ;mu_weights( 17)= 4.28138297104e-10 ;
    mu_points( 18)= 25.6286360225 ;mu_weights( 18)= 2.30589949189e-11 ;
    mu_points( 19)= 28.8621018163 ;mu_weights( 19)= 9.79937928873e-13 ;
    mu_points( 20)= 32.346629154 ;mu_weights( 20)= 3.23780165773e-14 ;
    mu_points( 21)= 36.1004948058 ;mu_weights( 21)= 8.17182344342e-16 ;
    mu_points( 22)= 40.1457197715 ;mu_weights( 22)= 1.54213383339e-17 ;
    mu_points( 23)= 44.5092079958 ;mu_weights( 23)= 2.11979229016e-19 ;
    mu_points( 24)= 49.2243949873 ;mu_weights( 24)= 2.05442967379e-21 ;
    mu_points( 25)= 54.3337213334 ;mu_weights( 25)= 1.34698258664e-23 ;
    mu_points( 26)= 59.8925091621 ;mu_weights( 26)= 5.6612941304e-26 ;
    mu_points( 27)= 65.9753772879 ;mu_weights( 27)= 1.41856054546e-28 ;
    mu_points( 28)= 72.6876280907 ;mu_weights( 28)= 1.91337549445e-31 ;
    mu_points( 29)= 80.1874469779 ;mu_weights( 29)= 1.1922487601e-34 ;
    mu_points( 30)= 98.8295428683 ;mu_weights( 30)= 1.33861694201e-42 ;
    mu_points( 31)= 88.7353404179 ;mu_weights( 31)= 2.67151121924e-38 ;



    
  endif





  
  !print *,"mu_points =",mu_points
  !print *,"mu_weights =",mu_weights
  
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
  call qn%precompute_qn( mu_points(0:N_mu-1), mu_weights(0:N_mu-1) , N_mu)
  
  
  stop
  call qn%solve_qn(phi_restr)





  
  phi_qn(1:Nc(1)+1,1:Nc(2)) = phi_restr(1:Nc(1)+1,1:Nc(2))
  phi_qn(1:Nc(1)+1,Nc(2)+1) = phi_restr(1:Nc(1)+1,1)
  
  
  call test_solve_qn_polar_new(Nc,eta_min,eta_max, &
    mu_points(0:N_mu-1),mu_weights(0:N_mu-1),N_mu,mode,phi_init,phi_qn)

  if(err==0)then    
    print *, '#PASSED'
  endif

end program