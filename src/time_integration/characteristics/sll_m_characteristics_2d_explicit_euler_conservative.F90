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

module sll_m_characteristics_2d_explicit_euler_conservative
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_set_to_limit, &
    sll_p_user_defined

  use sll_m_characteristics_2d_base, only: &
    sll_f_process_outside_point_periodic, &
    sll_f_process_outside_point_set_to_limit, &
    sll_i_signature_process_outside_point, &
    sll_c_characteristics_2d_base

  implicit none

  public :: &
    sll_f_new_explicit_euler_conservative_2d_charac

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type,extends(sll_c_characteristics_2d_base) :: &
    explicit_euler_conservative_2d_charac_computer
    sll_int32                               :: Npts1
    sll_int32                               :: Npts2
    sll_real64                              :: eta1_min   
    sll_real64                              :: eta1_max  
    sll_real64                              :: eta2_min   
    sll_real64                              :: eta2_max
    procedure(sll_i_signature_process_outside_point), pointer, nopass    :: &
      process_outside_point1
    procedure(sll_i_signature_process_outside_point), pointer, nopass    :: &
      process_outside_point2
    sll_real64, dimension(:,:,:), pointer :: buf_3d 
    sll_int32 :: bc_type_1
    sll_int32 :: bc_type_2

  contains
    !function, pass(charac) :: new => &
    !  sll_f_new_explicit_euler_conservative_2d_charac
    procedure, pass(charac) :: initialize => &
      initialize_explicit_euler_conservative_2d_charac
    procedure, pass(charac) :: compute_characteristics => &
      compute_explicit_euler_conservative_2d_charac
  end type explicit_euler_conservative_2d_charac_computer

contains
  function sll_f_new_explicit_euler_conservative_2d_charac(&
      Npts1, &
      Npts2, &
      bc_type_1, &
      bc_type_2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      process_outside_point1, &
      process_outside_point2) &
      result(charac)
      
    type(explicit_euler_conservative_2d_charac_computer),pointer :: charac
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_int32, intent(in), optional :: bc_type_1
    sll_int32, intent(in), optional :: bc_type_2
    sll_real64, intent(in), optional  :: eta1_min
    sll_real64, intent(in), optional  :: eta1_max
    sll_real64, intent(in), optional  :: eta2_min
    sll_real64, intent(in), optional  :: eta2_max
    procedure(sll_i_signature_process_outside_point), optional    :: &
      process_outside_point1
    procedure(sll_i_signature_process_outside_point), optional    :: &
      process_outside_point2
    sll_int32 :: ierr
      
    SLL_ALLOCATE(charac,ierr)
    call initialize_explicit_euler_conservative_2d_charac(&
      charac, &
      Npts1, &
      Npts2, &
      bc_type_1, &
      bc_type_2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      process_outside_point1, &
      process_outside_point2)

    
  end function sll_f_new_explicit_euler_conservative_2d_charac
  
  
  subroutine initialize_explicit_euler_conservative_2d_charac(&
      charac, &
      Npts1, &
      Npts2, &
      bc_type_1, &
      bc_type_2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      process_outside_point1, &
      process_outside_point2)
      
    class(explicit_euler_conservative_2d_charac_computer) :: charac
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_int32, intent(in), optional :: bc_type_1
    sll_int32, intent(in), optional :: bc_type_2
    sll_real64, intent(in), optional  :: eta1_min
    sll_real64, intent(in), optional  :: eta1_max
    sll_real64, intent(in), optional  :: eta2_min
    sll_real64, intent(in), optional  :: eta2_max
    procedure(sll_i_signature_process_outside_point), optional    :: &
      process_outside_point1
    procedure(sll_i_signature_process_outside_point), optional    :: &
      process_outside_point2
    
    sll_int32 :: ierr

    charac%Npts1 = Npts1
    charac%Npts2 = Npts2
    
    SLL_ALLOCATE(charac%buf_3d(2,0:Npts1,0:Npts2),ierr)
    
    
    
    
    if(present(eta1_min))then
      charac%eta1_min = eta1_min
    else
      charac%eta1_min = 0._f64
    endif
    if(present(eta1_max))then
      charac%eta1_max = eta1_max
    else
      charac%eta1_max = 1._f64
    endif
    if(present(eta2_min))then
      charac%eta2_min = eta2_min
    else
      charac%eta2_min = 0._f64  
    endif
    
    if(present(eta2_max))then
      charac%eta2_max = eta2_max
    else
      charac%eta2_max = 1._f64
    endif
    
    
    !charac%process_outside_point1 => process_outside_point1
    !charac%process_outside_point2 => process_outside_point2
    
    
    if(present(process_outside_point1)) then
      charac%process_outside_point1 => process_outside_point1
      charac%bc_type_1 = sll_p_user_defined
    else if(.not.(present(bc_type_1))) then
      print *,'#provide boundary condition'
      print *,'#bc_type_1 or process_outside_point1 function'
      print *,'#in initialize_explicit_euler_conservative_2d_charac'
      stop
    else
      charac%bc_type_1 = bc_type_1
      select case (bc_type_1)
        case (sll_p_periodic)
          charac%process_outside_point1 => sll_f_process_outside_point_periodic          
        case (sll_p_set_to_limit)
          charac%process_outside_point1 => sll_f_process_outside_point_set_to_limit        
        case default
          print *,'#bad value of boundary condition'
          print *,'#in initialize_explicit_euler_conservative_2d_charac_computer'
          stop
        end select
    endif
    
    if((present(process_outside_point1)).and.(present(bc_type_1)))then
      print *,'#provide either process_outside_point1 or bc_type_1'
      print *,'#and not both'
      print *,'#in initialize_explicit_euler_conservative_2d_charac_computer'
      stop
    endif
    


    if(present(process_outside_point2)) then
      charac%process_outside_point2 => process_outside_point2
      charac%bc_type_2 = sll_p_user_defined
    else if(.not.(present(bc_type_2))) then
      print *,'#provide boundary condition'
      print *,'#bc_type_2 or process_outside_point1 function'
      stop
    else
      charac%bc_type_2 = bc_type_2
    select case (bc_type_2)
        case (sll_p_periodic)
          charac%process_outside_point2 => sll_f_process_outside_point_periodic          
        case (sll_p_set_to_limit)
          charac%process_outside_point2 => sll_f_process_outside_point_set_to_limit        
        case default
          print *,'#bad value of boundary condition'
          print *,'#in initialize_explicit_euler_conservative_2d_charac_computer'
          stop
        end select
    endif

    if((present(process_outside_point2)).and.(present(bc_type_2)))then
      print *,'#provide either process_outside_point2 or bc_type_2'
      print *,'#and not both'
      print *,'#in initialize_explicit_euler_conservative_2d_charac_computer'
      stop
    endif
    
    
    
  end subroutine initialize_explicit_euler_conservative_2d_charac

  subroutine compute_explicit_euler_conservative_2d_charac( &
      charac, &
      A1, &
      A2, &
      dt, &
      input1, &
      input2, &
      output1, &
      output2)
            
    class(explicit_euler_conservative_2d_charac_computer) :: charac
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64, intent(in) :: dt
    sll_real64, dimension(:), intent(in) ::  input1
    sll_real64, dimension(:), intent(in) ::  input2
    sll_real64, dimension(:,:), intent(out) :: output1
    sll_real64, dimension(:,:), intent(out) :: output2    
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: Npts1
    sll_int32 :: Npts2
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max


    
    Npts1 = charac%Npts1
    Npts2 = charac%Npts2
    eta1_min = charac%eta1_min
    eta1_max = charac%eta1_max
    eta2_min = charac%eta2_min
    eta2_max = charac%eta2_max
    
    SLL_ASSERT(size(A1,1)>=charac%Npts1-1)
    SLL_ASSERT(size(A1,2)>=charac%Npts2-1)
    SLL_ASSERT(size(A2,1)>=charac%Npts1-1)
    SLL_ASSERT(size(A2,2)>=charac%Npts2-1)
    SLL_ASSERT(size(input1)>=charac%Npts1)
    SLL_ASSERT(size(input2)>=charac%Npts2)
    SLL_ASSERT(size(output1,1)>=charac%Npts1)
    SLL_ASSERT(size(output1,2)>=charac%Npts2)
    SLL_ASSERT(size(output2,1)>=charac%Npts1)
    SLL_ASSERT(size(output2,2)>=charac%Npts2)
    
    do j=1,Npts2-1
      do i=1,Npts1-1
        charac%buf_3d(1,i,j) = 0.5_f64*(input1(i)+input1(i+1))-dt*A1(i,j)
        !if((output1(i,j)<=eta1_min).or.(output1(i,j)>=eta1_max))then
        !  output1(i,j)=charac%process_outside_point1(output1(i,j),eta1_min,eta1_max)
        !endif  
         charac%buf_3d(2,i,j) = 0.5_f64*(input2(j)+input2(j+1))-dt*A2(i,j)  
        !if((output2(i,j)<=eta2_min).or.(output2(i,j)>=eta2_max))then
        !  output2(i,j)=charac%process_outside_point2(output2(i,j),eta2_min,eta2_max)          
        !endif
      enddo
    enddo   



    select case (charac%bc_type_1)
      case (sll_p_periodic)
        do j=1,Npts2-1 
          charac%buf_3d(1,0,j) = charac%buf_3d(1,Npts1-1,j) - (eta1_max-eta1_min)
          charac%buf_3d(2,0,j) = charac%buf_3d(2,Npts1-1,j)
          charac%buf_3d(1,Npts1,j) = charac%buf_3d(1,1,j) + (eta1_max-eta1_min)
          charac%buf_3d(2,Npts1,j) = charac%buf_3d(2,1,j)
        enddo
      case (sll_p_set_to_limit)
        do j=1,Npts2-1 
          charac%buf_3d(1,0,j) = 2._f64*eta1_min-charac%buf_3d(1,1,j)
          charac%buf_3d(2,0,j) = charac%buf_3d(2,1,j)
          charac%buf_3d(1,Npts1,j) = 2._f64*eta1_max - charac%buf_3d(1,Npts1-1,j)
          charac%buf_3d(2,Npts1,j) = charac%buf_3d(2,Npts1-1,j)
        enddo
      case default
        print *,'#bad value for charac%bc_type_1'
        stop
    end select

    select case (charac%bc_type_2)
      case (sll_p_periodic)
        do i=0,Npts1 
          charac%buf_3d(2,i,0) = charac%buf_3d(2,i,Npts2-1) - (eta2_max-eta2_min)
          charac%buf_3d(1,i,0) = charac%buf_3d(1,i,Npts2-1)
          charac%buf_3d(2,i,Npts1) = charac%buf_3d(2,i,1) + (eta2_max-eta2_min)
          charac%buf_3d(1,i,Npts1) = charac%buf_3d(1,i,1)
        enddo
      case (sll_p_set_to_limit)
        do i=0,Npts1 
          charac%buf_3d(2,i,0) = 2._f64*eta2_min-charac%buf_3d(2,i,1)
          charac%buf_3d(1,i,0) = charac%buf_3d(1,i,1)
          charac%buf_3d(2,i,Npts1) = 2._f64*eta2_max - charac%buf_3d(2,i,Npts1-1)
          charac%buf_3d(1,i,Npts1) = charac%buf_3d(1,i,Npts1-1)
        enddo
      case default
        print *,'#bad value for charac%bc_type_2'
        stop
    end select
    
    do j=1,Npts2
      do i=1,Npts1
        output1(i,j) = &
          0.25_f64*(charac%buf_3d(1,i,j) &
          +charac%buf_3d(1,i-1,j) &
          +charac%buf_3d(1,i-1,j-1) &
          +charac%buf_3d(1,i,j-1))
        output2(i,j) = &
          0.25_f64*(charac%buf_3d(2,i,j) &
          +charac%buf_3d(2,i-1,j) &
          +charac%buf_3d(2,i-1,j-1) &
          +charac%buf_3d(2,i,j-1))
      enddo
    enddo
    






      
  end subroutine compute_explicit_euler_conservative_2d_charac



  
  
end module sll_m_characteristics_2d_explicit_euler_conservative
