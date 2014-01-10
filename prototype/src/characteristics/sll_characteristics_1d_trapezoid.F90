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

!http://en.wikipedia.org/wiki/Trapezoidal_rule_%28differential_equations%29

module sll_module_characteristics_1d_trapezoid
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_characteristics_1d_base
use sll_module_interpolators_1d_base

implicit none

  type,extends(sll_characteristics_1d_base) :: trapezoid_1d_charac_computer
    sll_int32                               :: Npts
    sll_real64                              :: eta_min   
    sll_real64                              :: eta_max  
    procedure(signature_process_outside_point), pointer, nopass    :: &
      process_outside_point
    class(sll_interpolator_1d_base), pointer               :: A_interp
    sll_int32 :: maxiter
    sll_real64 :: tol
     
  contains
    procedure, pass(charac) :: initialize => &
      initialize_trapezoid_1d_charac
    procedure, pass(charac) :: compute_characteristics => &
      compute_trapezoid_1d_charac
  end type trapezoid_1d_charac_computer

contains
  function new_trapezoid_1d_charac(&
      Npts, &
      A_interp, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point, &
      maxiter, &
      tol) &
      result(charac)
      
    type(trapezoid_1d_charac_computer),pointer :: charac
    sll_int32, intent(in) :: Npts
    sll_int32, intent(in), optional :: bc_type
    sll_real64, intent(in), optional  :: eta_min
    sll_real64, intent(in), optional  :: eta_max
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point
    class(sll_interpolator_1d_base), target :: A_interp
    sll_int32, intent(in), optional :: maxiter
    sll_real64, intent(in), optional :: tol
    sll_int32 :: ierr
      
    SLL_ALLOCATE(charac,ierr)
    call initialize_trapezoid_1d_charac(&
      charac, &
      Npts, &
      A_interp, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point, &
      maxiter, &
      tol)

    
  end function new_trapezoid_1d_charac
  subroutine initialize_trapezoid_1d_charac(&
      charac, &
      Npts, &
      A_interp, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point, &
      maxiter, &
      tol)
      
    class(trapezoid_1d_charac_computer) :: charac
    sll_int32, intent(in) :: Npts
    sll_int32, intent(in), optional :: bc_type
    sll_real64, intent(in), optional  :: eta_min
    sll_real64, intent(in), optional  :: eta_max
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point
    class(sll_interpolator_1d_base), target :: A_interp
    sll_int32, intent(in), optional :: maxiter
    sll_real64, intent(in), optional :: tol


    charac%Npts = Npts
    
    
    if(present(eta_min))then
      charac%eta_min = eta_min
    else
      charac%eta_min = 0._f64
    endif
    if(present(eta_max))then
      charac%eta_max = eta_max
    else
      charac%eta_max = 1._f64
    endif
        
    
    if(present(process_outside_point)) then
      charac%process_outside_point => process_outside_point
    else if(.not.(present(bc_type))) then
      print *,'#provide boundary condition'
      print *,'#bc_type or process_outside_point function'
      print *,'#in initialize_trapezoid_1d_charac'
      stop
    else
      select case (bc_type)
        case (SLL_PERIODIC)
          charac%process_outside_point => process_outside_point_periodic          
        case (SLL_SET_TO_LIMIT)
          charac%process_outside_point => process_outside_point_set_to_limit        
        case default
          print *,'#bad value of boundary condition'
          print *,'#in initialize_trapezoid_1d_charac'
          stop
        end select
    endif

    if((present(process_outside_point)).and.(present(bc_type)))then
      print *,'#provide either process_outside_point or bc_type'
      print *,'#and not both'
      print *,'#in initialize_trapezoid_1d_charac'
      stop
    endif
    
    charac%A_interp => A_interp
    
    if(present(maxiter))then
      charac%maxiter = maxiter
    else
      charac%maxiter = 1000  
    endif


    if(present(tol))then
      charac%tol = tol
    else
      charac%tol = 1.e-12_f64  
    endif

    
    
  end subroutine initialize_trapezoid_1d_charac

  subroutine compute_trapezoid_1d_charac( &
      charac, &
      A, &
      dt, &
      input, &
      output)
            
    class(trapezoid_1d_charac_computer) :: charac
    sll_real64, dimension(:), intent(in) :: A
    sll_real64, intent(in) :: dt
    sll_real64, dimension(:), intent(in) ::  input
    sll_real64, dimension(:), intent(out) :: output
    sll_int32 :: j
    sll_real64 :: x2
    sll_real64 :: x2_old
    sll_real64 :: x2_i !i for inside, so that interpolation is possible
    sll_int32 :: iter
    sll_int32 :: Npts
    sll_real64 :: eta_min
    sll_real64 :: eta_max
    
    Npts = charac%Npts
    eta_min = charac%eta_min
    eta_max = charac%eta_max
    
    SLL_ASSERT(size(A)>=Npts)
    SLL_ASSERT(size(input)>=Npts)
    SLL_ASSERT(size(output)>=Npts)
    
    call charac%A_interp%compute_interpolants( &
      A, &
      input, &
      Npts)
    do j=1,Npts
        !We start from Y(t_{n+1})=y_j
        !and look for Y(t_n) = Yn
        !Yn = y_j-(A(y_j)+A(Yn))*dt/2
        x2 = input(j)-dt*A(j)
        x2_old = 0._f64
        iter = 0
        do while (iter<charac%maxiter .and. abs(x2_old-x2)>charac%tol)
          x2_old = x2
          x2_i = x2
          if((x2<=eta_min).or.(x2>=eta_max))then
            x2_i= charac%process_outside_point(x2,eta_min,eta_max)
          else
            x2_i = x2  
          endif                      
          x2 = input(j)-0.5_f64*dt*(charac%A_interp%interpolate_value(x2_i)+A(j))
        end do
        if (iter==charac%maxiter .and. abs(x2_old-x2)>charac%tol) then
          print*,'#not enough iterations for compute_trapezoid_1d_charac',&
            iter,abs(x2_old-x2)
          stop
        end if
        if((x2<=eta_min).or.(x2>=eta_max))then
          x2 =  charac%process_outside_point(x2,eta_min,eta_max)
        endif                      
        output(j) = x2  
    enddo
      
  end subroutine compute_trapezoid_1d_charac



  
  
end module sll_module_characteristics_1d_trapezoid
