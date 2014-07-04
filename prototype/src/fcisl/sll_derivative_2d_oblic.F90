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

! in development
! compute lim_{h->0} (phi(x1+A1h,x2+A2h)-phi(x1,x2))/h
! we explore here the spaghetti to compute this derivative


module sll_module_derivative_2d_oblic
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_fcisl_module

implicit none

  type  :: oblic_2d_derivative
    sll_int32 :: Nc_x1
    class(sll_advection_1d_base), pointer :: adv_x1
    sll_int32 :: Nc_x2
    sll_real64 :: x2_min    
    sll_real64 :: x2_max
    sll_int32 :: stencil_r
    sll_int32 :: stencil_s
    sll_real64, dimension(:), pointer :: weights
    sll_real64, dimension(:,:), pointer :: buf
        
  end type oblic_2d_derivative

contains
  
  function new_oblic_2d_derivative( &
    Nc_x1, &
    adv_x1, &
    Nc_x2, &
    x2_min, &
    x2_max, &
    stencil_r, &
    stencil_s ) &
    result(deriv)      
    type(oblic_2d_derivative), pointer :: deriv
    sll_int32, intent(in) :: Nc_x1
    class(sll_advection_1d_base), pointer :: adv_x1
    sll_int32, intent(in) :: Nc_x2
    sll_real64, intent(in) :: x2_min
    sll_real64, intent(in) :: x2_max
    sll_int32, intent(in) :: stencil_r
    sll_int32, intent(in) :: stencil_s
    !local variables
    sll_int32 :: ierr
    SLL_ALLOCATE(deriv,ierr)
    
    call initialize_oblic_2d_derivative( &
      deriv, &
      Nc_x1, &
      adv_x1, & 
      Nc_x2, &
      x2_min, &
      x2_max, &
      stencil_r, &
      stencil_s )
    
  end function new_oblic_2d_derivative

  subroutine initialize_oblic_2d_derivative( &
    deriv, &
    Nc_x1, &
    adv_x1, &
    Nc_x2, &
    x2_min, &
    x2_max, &
    stencil_r, &
    stencil_s )
    type(oblic_2d_derivative), intent(inout) :: deriv
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
    
    deriv%Nc_x1 = Nc_x1
    deriv%adv_x1 => adv_x1
    deriv%Nc_x2 = Nc_x2
    deriv%x2_min = x2_min
    deriv%x2_max = x2_max
    deriv%stencil_r = stencil_r
    deriv%stencil_s = stencil_s
    r = stencil_r
    s = stencil_s
    
    SLL_ALLOCATE(deriv%weights(r:s),ierr)
    SLL_ALLOCATE(deriv%buf(r:s,Nc_x1+1),ierr)
    
    call compute_w_hermite(deriv%weights(r:s),r,s)  

    
  end subroutine initialize_oblic_2d_derivative
  
!< compute lim_{h->0} (f(x1_i+A1h,x2_j+A2h)-f(x1,x2))/h 
!< = A1*Dx1_f(i,j)+A2*Dx2_f(i,j)
!< for i=1,..,Nc_x1+1, j=1,..,Nc_x2+1 
!< note that iota= (A1/(x1_max-x1_min))/(A2/(x2_max-x2_min)) 
  subroutine compute_oblic_derivative_2d(&
    deriv, &
    A1, &
    A2, &
    input, &
    output)
    type(oblic_2d_derivative), pointer :: deriv
    sll_real64, intent(in) :: A1
    sll_real64, intent(in) :: A2
    sll_real64, dimension(:,:), intent(in) :: input
    sll_real64, dimension(:,:), intent(out) :: output
    !local variables
    sll_int32 :: i2
    sll_int32 :: ell
    sll_int32 :: i2_loc
    sll_real64, dimension(:,:), pointer :: buf
    sll_real64, dimension(:), pointer :: w
    class(sll_advection_1d_base), pointer :: adv_x1
    sll_int32 :: r
    sll_int32 :: s
    sll_real64 :: delta_x2
    sll_real64 :: dt_loc
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: length
    sll_int32 :: i1
    sll_int32 :: ind
    sll_real64 :: tmp
    
    
    Nc_x1 = deriv%Nc_x1
    Nc_x2 = deriv%Nc_x2
    length = (deriv%x2_max-deriv%x2_min)
    delta_x2 = length/real(Nc_x2,f64)
    
    buf => deriv%buf
    w => deriv%weights
    adv_x1 => deriv%adv_x1     
    r = deriv%stencil_r
    s = deriv%stencil_s
    
    
    do i2=1,Nc_x2+1
      !choose several dt_loc so that advection in x2 is exact
      do ell=r,s
        dt_loc = -real(ell,f64)*delta_x2/A2         
        i2_loc = modulo(i2+ell-1,Nc_x2)+1
        call adv_x1%advect_1d_constant( &
          A1, &
          dt_loc, &
          input(1:Nc_x1+1,i2_loc), &
          buf(ell,1:Nc_x1+1))
      enddo
      !compute the derivative using these values 
      do i1=1,Nc_x1+1
        tmp=0._f64
        do ell=r,s
          tmp=tmp+w(ell)*buf(ell,i1)
        enddo
        output(i1,i2) = tmp*A2/delta_x2
      enddo
    enddo    


    
    
  end subroutine compute_oblic_derivative_2d  
  
  
end module sll_module_derivative_2d_oblic
