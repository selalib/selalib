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

! in development; use oblic interpolation
! data are on uniform (fine) grid in x1

module sll_module_advection_2d_oblic
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
!use sll_module_advection_2d_base
use sll_module_characteristics_2d_base
use sll_module_interpolators_2d_base
use sll_fcisl_module
use lagrange_interpolation

implicit none

  type  :: oblic_2d_advector
    sll_int32 :: Nc_x1
    class(sll_advection_1d_base), pointer :: adv_x1
    sll_int32 :: Nc_x2
    sll_real64  :: x2_min
    sll_real64  :: x2_max
    sll_int32 :: stencil_r
    sll_int32 :: stencil_s
    sll_real64, dimension(:), pointer :: xx
    sll_real64, dimension(:,:), pointer :: buf
  
  end type oblic_2d_advector
   




contains
  function new_oblic_2d_advector( &
    Nc_x1, &
    adv_x1, &
    Nc_x2, &
    x2_min, &
    x2_max, &
    stencil_r, &
    stencil_s ) &
    result(adv)      
    type(oblic_2d_advector), pointer :: adv
    sll_int32, intent(in) :: Nc_x1
    class(sll_advection_1d_base), pointer :: adv_x1
    sll_int32, intent(in) :: Nc_x2
    sll_real64, intent(in) :: x2_min
    sll_real64, intent(in) :: x2_max
    sll_int32, intent(in) :: stencil_r
    sll_int32, intent(in) :: stencil_s
    !local variables
    sll_int32 :: ierr
    
    SLL_ALLOCATE(adv,ierr)
        
    call initialize_oblic_2d_advector( &
      adv, &
      Nc_x1, &
      adv_x1, &
      Nc_x2, &
      x2_min, &
      x2_max, &
      stencil_r, &
      stencil_s )
  end function  new_oblic_2d_advector


  subroutine initialize_oblic_2d_advector( &
    adv, &
    Nc_x1, &
    adv_x1, &
    Nc_x2, &
    x2_min, &
    x2_max, &
    stencil_r, &
    stencil_s )
    type(oblic_2d_advector), intent(inout) :: adv
    sll_int32, intent(in) :: Nc_x1
    class(sll_advection_1d_base), pointer :: adv_x1
    sll_int32, intent(in) :: Nc_x2
    sll_real64, intent(in) :: x2_min
    sll_real64, intent(in) :: x2_max
    sll_int32, intent(in) :: stencil_r
    sll_int32, intent(in) :: stencil_s
    !local variables
    sll_int32 :: ierr
    sll_int32 :: r
    sll_int32 :: s
    sll_int32 :: i
        
    adv%Nc_x1 = Nc_x1
    adv%adv_x1 => adv_x1
    adv%Nc_x2 = Nc_x2
    adv%x2_min = x2_min
    adv%x2_max = x2_max
    adv%stencil_r = stencil_r
    adv%stencil_s = stencil_s
    
    r = adv%stencil_r
    s = adv%stencil_s
    
    SLL_ALLOCATE(adv%xx(r:s),ierr)
    SLL_ALLOCATE(adv%buf(r:s,Nc_x1+1),ierr)
    
    do i=r,s
      adv%xx(i) = real(i,f64)
    enddo

      
  end subroutine initialize_oblic_2d_advector


!> solves \partial_t f + \nabla A \cdot f = 0, \ A = (A1,A2) for time step dt
!> interpolation in aligned along iota = A1/A2
!> iota is the number of tours that has been done at x2_max
!> the direction of iota is the line passing through 
!> (x1_min,x2_min) to (x1_max,x2_min+iota*(x2_max-x2_min))
!> here iota is real number
!> for specific cases where x2_min+iota*(x2_max-x2_min) is a mesh point
!> we refer to sll_advection_2d_integer_oblic
!> here A1 and A2 are real numbers
!> we refer to (future) oblic_advect_2d, when (A1,A2) are 2D arrays
!> that is A1 and A2 depend both on x1 and x2
!> we use here constant advection in the x1 direction
!> with an abstract 1d advector that can do the constant advection
!> in the x2 direction, we use (for the moment) Lagrange interpolation with stencil (r,s)
!> (r,s) = (0,1) : LAG1
!> (r,s) = (-1,2) : LAG3
!> (r,s) = (-2,3) : LAG5
!> periodic conditions are used in both x1 and x2 directions


  subroutine oblic_advect_2d_constant(&
    adv, &
    !iota, &
    A1, &
    A2, &
    dt, &
    input, &
    output)
    type(oblic_2d_advector) :: adv
    !sll_real64, intent(in) :: iota
    sll_real64, intent(in) :: A1
    sll_real64, intent(in) :: A2
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:,:), intent(in) :: input
    sll_real64, dimension(:,:), intent(out) :: output      
    !local variables
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i2_loc
    sll_int32 :: i0
    sll_real64 :: alpha
    sll_int32 :: d
    sll_int32 :: r
    sll_int32 :: s
    sll_real64 :: dt_loc
    sll_int32 :: ell
    sll_real64, dimension(:,:), pointer :: buf
    sll_real64, dimension(:), pointer :: xx
    class(sll_advection_1d_base), pointer :: adv_x1
    sll_real64 :: delta_x2
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    
    buf => adv%buf
    xx => adv%xx
    adv_x1 => adv%adv_x1
    Nc_x1 = adv%Nc_x1
    Nc_x2 = adv%Nc_x2
    delta_x2 = (adv%x2_max-adv%x2_min)/real(adv%Nc_x2,f64)
    r = adv%stencil_r
    s = adv%stencil_s    
    alpha = A2*dt/delta_x2
    i0 = floor(alpha)
    alpha = alpha-i0  
    d =s-r
    
    
    do i2=1,Nc_x2+1
      !choose several dt_loc so that advection in x2 is exact
      do ell=r,s
        dt_loc = real(ell+i0,f64)*delta_x2/A2         
        i2_loc = modulo(i2-ell-i0-1,Nc_x2)+1
        call adv_x1%advect_1d_constant( &
          A1, &
          dt_loc, &
          input(1:Nc_x1+1,i2_loc), &
          buf(ell,1:Nc_x1+1))
      enddo
      !interpolate between these values 
      do i1=1,Nc_x1+1
        output(i1,i2) = lagrange_interpolate(alpha, d, xx, buf(r:s,i1) )
      enddo
    enddo    



          
  end subroutine oblic_advect_2d_constant

  subroutine oblic_advect_2d(&
    adv, &
    !iota, &
    A1, &
    A2, &
    dt, &
    input, &
    output)
    type(oblic_2d_advector) :: adv
    !sll_real64, intent(in) :: iota
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:,:), intent(in) :: input
    sll_real64, dimension(:,:), intent(out) :: output      
    
    print *,'#oblic_advect_2d not implemented for the moment'
          
  end subroutine oblic_advect_2d




end module sll_module_advection_2d_oblic
