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

!Here is an example of rotation test case



program rotation_2d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_module_advection_2d_BSL
use sll_module_characteristics_2d_explicit_euler
use sll_module_characteristics_2d_verlet
use sll_cubic_spline_interpolator_2d
use sll_cubic_spline_interpolator_1d

implicit none

  sll_int32 :: Nc_x1
  sll_int32 :: Nc_x2
  sll_real64 :: err
  sll_int32 :: nb_step
  sll_int32 :: step
  sll_real64,dimension(:,:),allocatable :: f
  sll_real64,dimension(:,:),allocatable :: f_exact
  sll_real64,dimension(:,:),allocatable :: f_init
  sll_real64,dimension(:,:),allocatable :: f_old
  sll_real64,dimension(:,:),allocatable :: A1
  sll_real64,dimension(:,:),allocatable :: A2
  sll_int32 ::ierr
  sll_int32 ::i1
  sll_int32 ::i2
  sll_real64 :: x1
  sll_real64 :: x2
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  sll_real64 :: x2_min
  sll_real64 :: x2_max
  sll_real64 :: delta_x1
  sll_real64 :: delta_x2
  sll_real64 :: dt
  
  class(sll_advection_2d_base), pointer :: adv
  class(sll_interpolator_2d_base), pointer :: interp
  class(sll_characteristics_2d_base), pointer :: charac
  class(sll_interpolator_2d_base), pointer   :: A1_interp_x1x2
  class(sll_interpolator_2d_base), pointer   :: A2_interp_x1x2
  class(sll_interpolator_1d_base), pointer   :: A1_interp_x1
  class(sll_interpolator_1d_base), pointer   :: A2_interp_x1
  
  
  nb_step = 1000
  dt = 0.01_f64
  Nc_x1 = 64
  Nc_x2 = 64
  
  x1_max = 10._f64
  x1_min = -x1_max
  x2_min = x1_min
  x2_max = x1_max
   
  delta_x1 = (x1_max-x1_min)/real(Nc_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(Nc_x2,f64)
  
  
  
  SLL_ALLOCATE(f(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_old(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_init(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_exact(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(A1(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(A2(Nc_x1+1,Nc_x2+1),ierr)
  
  !we initialize the distribution function
  
  do i2=1,Nc_x2
    do i1=1,Nc_x1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      f_init(i1,i2) = exp(-0.5_f64*(x1**2+x2**2))
      f_exact(i1,i2) = exp(-0.5_f64*(x1**2+x2**2))
      A1(i1,i2) = -x2
      A2(i1,i2) = x1
    enddo
  enddo
  
  f = f_init
  
  !we initialize the interpolator
  interp => new_cubic_spline_2d_interpolator( &
    Nc_x1+1, &
    Nc_x2+1, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max, &
    SLL_HERMITE, &
    SLL_HERMITE)

  



  !we initialize the characteristics
  A1_interp_x1 => new_cubic_spline_1d_interpolator( &
    Nc_x1+1, &
    x1_min, &
    x1_max, &
    SLL_HERMITE)
  A2_interp_x1 => new_cubic_spline_1d_interpolator( &
    Nc_x1+1, &
    x1_min, &
    x1_max, &
    SLL_HERMITE)
  A1_interp_x1x2 => new_cubic_spline_2d_interpolator( &
    Nc_x1+1, &
    Nc_x2+1, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max, &
    SLL_HERMITE, &
    SLL_HERMITE)
  A2_interp_x1x2 => new_cubic_spline_2d_interpolator( &
    Nc_x1+1, &
    Nc_x2+1, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max, &
    SLL_HERMITE, &
    SLL_HERMITE)
  charac => new_verlet_2d_charac(&
    Nc_x1+1, &
    Nc_x2+1, &
    A1_interp_x1x2, &
    A2_interp_x1x2, &
    A1_interp_x1, &
    A2_interp_x1, &
    bc_type_1=SLL_SET_TO_LIMIT, &
    bc_type_2=SLL_SET_TO_LIMIT)



  adv => new_BSL_2d_advector(&
    interp, &
    charac, &
    Nc_x1+1, &
    Nc_x2+1, &
    eta1_min = x1_min, &
    eta1_max = x1_max, &
    eta2_min = x2_min, &
    eta2_max = x2_max)
  






  
  
  err=0._f64
  
  do step=1,nb_step
    f_old = f
    call adv%advect_2d(A1, A2, dt, f_old, f)
  enddo
  
  err = maxval(abs(f-f_old))
  
  print *,'#err=',err

  print *,dt,Nc_x1,Nc_x2,err


end program


