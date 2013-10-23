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

module sll_module_characteristics_2d_verlet
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_characteristics_2d_base
use sll_module_interpolators_2d_base
use sll_module_interpolators_1d_base

implicit none

  type,extends(sll_characteristics_2d_base) :: verlet_2d_charac_computer
    sll_int32                               :: Npts1
    sll_int32                               :: Npts2
    sll_real64                              :: eta1_min   
    sll_real64                              :: eta1_max  
    sll_real64                              :: eta2_min   
    sll_real64                              :: eta2_max
    procedure(signature_process_outside_point), pointer, nopass    :: &
      process_outside_point1
    procedure(signature_process_outside_point), pointer, nopass    :: &
      process_outside_point2
    class(sll_interpolator_2d_base), pointer               :: A1_interp_x1x2
    class(sll_interpolator_2d_base), pointer               :: A2_interp_x1x2
    class(sll_interpolator_1d_base), pointer               :: A1_interp_x1
    class(sll_interpolator_1d_base), pointer               :: A2_interp_x1
    sll_int32 :: x1_maxiter
    sll_int32 :: x2_maxiter
    sll_real64 :: x1_tol
    sll_real64 :: x2_tol
     
  contains
    procedure, pass(charac) :: initialize => &
      initialize_verlet_2d_charac
    procedure, pass(charac) :: compute_characteristics => &
      compute_verlet_2d_charac
  end type verlet_2d_charac_computer

contains
  function new_verlet_2d_charac(&
      Npts1, &
      Npts2, &
      A1_interp_x1x2, &
      A2_interp_x1x2, &
      A1_interp_x1, &
      A2_interp_x1, &
      bc_type_1, &
      bc_type_2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      process_outside_point1, &
      process_outside_point2, &
      x1_maxiter, &
      x2_maxiter, &
      x1_tol, &
      x2_tol) &
      result(charac)
      
    type(verlet_2d_charac_computer),pointer :: charac
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_int32, intent(in), optional :: bc_type_1
    sll_int32, intent(in), optional :: bc_type_2
    sll_real64, intent(in), optional  :: eta1_min
    sll_real64, intent(in), optional  :: eta1_max
    sll_real64, intent(in), optional  :: eta2_min
    sll_real64, intent(in), optional  :: eta2_max
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point1
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point2
    class(sll_interpolator_2d_base), target :: A1_interp_x1x2
    class(sll_interpolator_2d_base), target :: A2_interp_x1x2
    class(sll_interpolator_1d_base), target :: A1_interp_x1
    class(sll_interpolator_1d_base), target :: A2_interp_x1
    sll_int32, intent(in), optional :: x1_maxiter
    sll_int32, intent(in), optional :: x2_maxiter
    sll_real64, intent(in), optional :: x1_tol
    sll_real64, intent(in), optional :: x2_tol
    sll_int32 :: ierr
      
    SLL_ALLOCATE(charac,ierr)
    call initialize_verlet_2d_charac(&
      charac, &
      Npts1, &
      Npts2, &
      A1_interp_x1x2, &
      A2_interp_x1x2, &
      A1_interp_x1, &
      A2_interp_x1, &
      bc_type_1, &
      bc_type_2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      process_outside_point1, &
      process_outside_point2, &
      x1_maxiter, &
      x2_maxiter, &
      x1_tol, &
      x2_tol)

    
  end function new_verlet_2d_charac
  subroutine initialize_verlet_2d_charac(&
      charac, &
      Npts1, &
      Npts2, &
      A1_interp_x1x2, &
      A2_interp_x1x2, &
      A1_interp_x1, &
      A2_interp_x1, &
      bc_type_1, &
      bc_type_2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      process_outside_point1, &
      process_outside_point2, &
      x1_maxiter, &
      x2_maxiter, &
      x1_tol, &
      x2_tol)
      
    class(verlet_2d_charac_computer) :: charac
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_int32, intent(in), optional :: bc_type_1
    sll_int32, intent(in), optional :: bc_type_2
    sll_real64, intent(in), optional  :: eta1_min
    sll_real64, intent(in), optional  :: eta1_max
    sll_real64, intent(in), optional  :: eta2_min
    sll_real64, intent(in), optional  :: eta2_max
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point1
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point2
    class(sll_interpolator_2d_base), target :: A1_interp_x1x2
    class(sll_interpolator_2d_base), target :: A2_interp_x1x2
    class(sll_interpolator_1d_base), target :: A1_interp_x1
    class(sll_interpolator_1d_base), target :: A2_interp_x1
    sll_int32, intent(in), optional :: x1_maxiter
    sll_int32, intent(in), optional :: x2_maxiter
    sll_real64, intent(in), optional :: x1_tol
    sll_real64, intent(in), optional :: x2_tol


    charac%Npts1 = Npts1
    charac%Npts2 = Npts2
    
    
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
    else if(.not.(present(bc_type_1))) then
      print *,'#provide boundary condition'
      print *,'#bc_type_1 or process_outside_point1 function'
      print *,'#in initialize_verlet_2d_charac_computer'
      stop
    else
      select case (bc_type_1)
        case (SLL_PERIODIC)
          charac%process_outside_point1 => process_outside_point_periodic          
        case (SLL_SET_TO_LIMIT)
          charac%process_outside_point1 => process_outside_point_set_to_limit        
        case default
          print *,'#bad value of boundary condition'
          print *,'#in initialize_verlet_2d_charac_computer'
          stop
        end select
    endif

    if((present(process_outside_point1)).and.(present(bc_type_1)))then
      print *,'#provide either process_outside_point1 or bc_type_1'
      print *,'#and not both'
      print *,'#in initialize_verlet_2d_charac_computer'
      stop
    endif



    if(present(process_outside_point2)) then
      charac%process_outside_point2 => process_outside_point2
    else if(.not.(present(bc_type_2))) then
      print *,'#provide boundary condition'
      print *,'#bc_type_2 or process_outside_point1 function'
      print *,'#in initialize_verlet_2d_charac_computer'
      stop
    else
      select case (bc_type_2)
        case (SLL_PERIODIC)
          charac%process_outside_point2 => process_outside_point_periodic          
        case (SLL_SET_TO_LIMIT)
          charac%process_outside_point2 => process_outside_point_set_to_limit        
        case default
          print *,'#bad value of boundary condition'
          print *,'#in initialize_verlet_2d_charac_computer'
          stop
        end select
    endif

    if((present(process_outside_point2)).and.(present(bc_type_2)))then
      print *,'#provide either process_outside_point2 or bc_type_2'
      print *,'#and not both'
      print *,'#in initialize_verlet_2d_charac_computer'
      stop
    endif

    
    
    charac%A1_interp_x1x2 => A1_interp_x1x2
    charac%A2_interp_x1x2 => A2_interp_x1x2
    charac%A1_interp_x1 => A1_interp_x1
    charac%A2_interp_x1 => A2_interp_x1
    
    
    if(present(x1_maxiter))then
      charac%x1_maxiter = x1_maxiter
    else
      charac%x1_maxiter = 1000  
    endif
    if(present(x2_maxiter))then
      charac%x2_maxiter = x2_maxiter
    else
      charac%x2_maxiter = charac%x1_maxiter  
    endif


    if(present(x1_tol))then
      charac%x1_tol = x1_tol
    else
      charac%x1_tol = 1.e-12_f64  
    endif
    if(present(x2_tol))then
      charac%x2_tol = x2_tol
    else
      charac%x2_tol = charac%x1_tol  
    endif


    
    
  end subroutine initialize_verlet_2d_charac

  subroutine compute_verlet_2d_charac( &
      charac, &
      A1, &
      A2, &
      dt, &
      input1, &
      input2, &
      output1, &
      output2)
            
    class(verlet_2d_charac_computer) :: charac
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64, intent(in) :: dt
    sll_real64, dimension(:), intent(in) ::  input1
    sll_real64, dimension(:), intent(in) ::  input2
    sll_real64, dimension(:,:), intent(out) :: output1
    sll_real64, dimension(:,:), intent(out) :: output2    
    sll_int32 :: i
    sll_int32 :: j
    sll_real64 :: x1
    sll_real64 :: x2
    sll_real64 :: x1_old
    sll_real64 :: x2_old
    sll_real64 :: x1_i !i for inside, so that interpolation is possible
    sll_real64 :: x2_i
    sll_int32 :: iter
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
    
    SLL_ASSERT(size(A1,1)>=Npts1)
    SLL_ASSERT(size(A1,2)>=Npts2)
    SLL_ASSERT(size(A2,1)>=Npts1)
    SLL_ASSERT(size(A2,2)>=Npts2)
    SLL_ASSERT(size(input1)>=Npts1)
    SLL_ASSERT(size(input2)>=Npts2)
    SLL_ASSERT(size(output1,1)>=Npts1)
    SLL_ASSERT(size(output1,2)>=Npts2)
    SLL_ASSERT(size(output2,1)>=Npts1)
    SLL_ASSERT(size(output2,2)>=Npts2)
    
    call charac%A1_interp_x1x2%compute_interpolants( &
      A1, &
      input1, &
      Npts1, &
      input2, &
      Npts2)
    call charac%A2_interp_x1x2%compute_interpolants( &
      A2, &
      input1, &
      Npts1, &
      input2, &
      Npts2)
    do j=1,charac%Npts2

      call charac%A1_interp_x1%compute_interpolants( A1(:,j), input1, charac%Npts1 )
      call charac%A2_interp_x1%compute_interpolants( A2(:,j), input1, charac%Npts1 )
      
      do i=1,charac%Npts1
        !We start from X(t_{n+1}) = x_i, Y(t_{n+1})=y_j
        !and look for X(t_n) = Xn, Y(t_n) = Yn
        
        !X* = x_i-A1(X*,y_j)*dt/2
        x1 = input1(i)-0.5_f64*dt*A1(i,j)
        if((x1<=charac%eta1_min).or.(x1>=charac%eta1_max))then
          x1 = charac%process_outside_point1(x1,charac%eta1_min,charac%eta1_max)
        endif  
        x1_old = 0._f64
        iter = 0
        do while (iter<charac%x1_maxiter .and. abs(x1_old-x1)>charac%x1_tol)
          x1_old = x1
          if((x1<=charac%eta1_min).or.(x1>=charac%eta1_max))then
            x1_i = charac%process_outside_point1(x1,charac%eta1_min,charac%eta1_max)
          else
            x1_i = x1  
          endif            
          x1 = input1(i)-0.5_f64*dt*charac%A1_interp_x1%interpolate_value(x1_i)
        end do
        if (iter==charac%x1_maxiter .and. abs(x1_old-x1)>charac%x1_tol) then
          print*,'#not enough iterations for compute_characteristics2D_verlet x1',&
            iter,abs(x1_old-x1)
          stop
        end if
        if((x1<=charac%eta1_min).or.(x1>=charac%eta1_max))then
          x1 =  charac%process_outside_point1(x1,charac%eta1_min,charac%eta1_max)
        endif            
                
        !Yn = y_j-(A2(X*,y_j)+A2(X*,Yn))*dt/2
        x2 = input2(j)-dt*A2(i,j)
        x2_old = 0._f64
        iter = 0
        do while (iter<charac%x2_maxiter .and. abs(x2_old-x2)>charac%x2_tol)
          x2_old = x2
          x2_i = x2
          if((x2<=charac%eta2_min).or.(x2>=charac%eta2_max))then
            x2_i= charac%process_outside_point2(x2,charac%eta2_min,charac%eta2_max)
          else
            x2_i = x2  
          endif                      
          x2 = input2(j)-0.5_f64*dt*(charac%A2_interp_x1x2%interpolate_value(x1, x2_i)&
            +charac%A2_interp_x1%interpolate_value( x1))
        end do
        if (iter==charac%x2_maxiter .and. abs(x2_old-x2)>charac%x2_tol) then
          print*,'#not enough iterations for compute_characteristics2D_verlet x2',&
            iter,abs(x2_old-x2)
          stop
        end if
        if((x2<=charac%eta2_min).or.(x2>=charac%eta2_max))then
          x2 =  charac%process_outside_point2(x2,charac%eta2_min,charac%eta2_max)
        endif                      
        
        !Xn = y_j-A1(X*,Yn)*dt/2
        x1 = x1-0.5_f64*dt*charac%A1_interp_x1x2%interpolate_value( x1, x2)
        if((x1<=charac%eta1_min).or.(x1>=charac%eta1_max))then
          x1 = charac%process_outside_point1(x1,charac%eta1_min,charac%eta1_max)
        endif            
        
        output1(i,j) = x1 
        output2(i,j) = x2  
      enddo
    enddo        
      
  end subroutine compute_verlet_2d_charac



  
  
end module sll_module_characteristics_2d_verlet
