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

module sll_module_characteristics_1d_explicit_euler
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_characteristics_1d_base
implicit none

  type,extends(sll_characteristics_1d_base) :: explicit_euler_1d_charac_computer
    sll_int32                               :: Npts
    sll_real64                              :: eta_min   
    sll_real64                              :: eta_max  
    procedure(signature_process_outside_point_1d), pointer, nopass    :: &
      process_outside_point
     
  contains
    procedure, pass(charac) :: initialize => &
      initialize_explicit_euler_1d_charac
    procedure, pass(charac) :: compute_characteristics => &
      compute_explicit_euler_1d_charac
  end type explicit_euler_1d_charac_computer

contains
  function new_explicit_euler_1d_charac(&
      Npts, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point) &
      result(charac)
      
    type(explicit_euler_1d_charac_computer),pointer :: charac
    sll_int32, intent(in) :: Npts
    sll_int32, intent(in), optional :: bc_type
    sll_real64, intent(in), optional  :: eta_min
    sll_real64, intent(in), optional  :: eta_max
    procedure(signature_process_outside_point_1d), optional    :: &
      process_outside_point
    sll_int32 :: ierr
      
    SLL_ALLOCATE(charac,ierr)
    call initialize_explicit_euler_1d_charac(&
      charac, &
      Npts, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point)

    
  end function new_explicit_euler_1d_charac
  
  
  subroutine initialize_explicit_euler_1d_charac(&
      charac, &
      Npts, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point)
      
    class(explicit_euler_1d_charac_computer) :: charac
    sll_int32, intent(in) :: Npts
    sll_int32, intent(in), optional :: bc_type
    sll_real64, intent(in), optional  :: eta_min
    sll_real64, intent(in), optional  :: eta_max
    procedure(signature_process_outside_point_1d), optional    :: &
      process_outside_point


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
      print *,'#in initialize_explicit_euler_1d_charac'
      stop
    else
      select case (bc_type)
        case (SLL_PERIODIC)
          charac%process_outside_point => process_outside_point_periodic          
        case (SLL_SET_TO_LIMIT)
          charac%process_outside_point => process_outside_point_set_to_limit        
        case default
          print *,'#bad value of boundary condition'
          print *,'#in initialize_explicit_euler_1d_charac'
          stop
        end select
    endif
    
    if((present(process_outside_point)).and.(present(bc_type)))then
      print *,'#provide either process_outside_point or bc_type'
      print *,'#and not both'
      print *,'#in initialize_explicit_euler_2d_charac'
      stop
    endif
    


    
    
  end subroutine initialize_explicit_euler_1d_charac

  subroutine compute_explicit_euler_1d_charac( &
      charac, &
      A, &
      dt, &
      input, &
      output)
            
    class(explicit_euler_1d_charac_computer) :: charac
    sll_real64, dimension(:), intent(in) :: A
    sll_real64, intent(in) :: dt
    sll_real64, dimension(:), intent(in) ::  input
    sll_real64, dimension(:), intent(out) :: output
    sll_int32 :: i
    sll_int32 :: Npts
    sll_real64 :: eta_min
    sll_real64 :: eta_max
    
    Npts = charac%Npts
    eta_min = charac%eta_min
    eta_max = charac%eta_max
    
    SLL_ASSERT(size(A)>=charac%Npts)
    SLL_ASSERT(size(input)>=charac%Npts)
    SLL_ASSERT(size(output)>=charac%Npts)
    
    do i=1,Npts
      output(i) = input(i)-dt*A(i)
      if((output(i)<=eta_min).or.(output(i)>=eta_max))then
        output(i)=charac%process_outside_point(output(i),eta_min,eta_max)
      endif  
    enddo
      
  end subroutine compute_explicit_euler_1d_charac



  
  
end module sll_module_characteristics_1d_explicit_euler
