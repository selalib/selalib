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

! in development; should be at least cubic splines
! attached with computation of characteristics


module sll_module_advection_1d_CSL
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_advection_1d_base
use sll_module_characteristics_1d_base
use sll_module_interpolators_1d_base
use sll_constants
implicit none

  type,extends(sll_advection_1d_base) :: CSL_1d_advector
  
    class(sll_interpolator_1d_base), pointer  :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_real64, dimension(:), pointer :: eta_coords
    sll_real64, dimension(:), pointer :: eta_coords_unit
    sll_real64, dimension(:), pointer :: eta_coords_unit_back
    sll_real64, dimension(:), pointer :: charac_feet
    sll_real64, dimension(:), pointer :: charac_feet_i
    sll_real64, dimension(:), pointer :: buf1d
    sll_real64, dimension(:), pointer :: buf1d_out
    sll_int32 :: Npts
    sll_int32 :: bc_type  
    procedure(signature_process_outside_point_1d), pointer, nopass    :: &
      process_outside_point

  contains
    procedure, pass(adv) :: initialize => &
       initialize_CSL_1d_advector
    procedure, pass(adv) :: advect_1d => &
      CSL_advect_1d
    procedure, pass(adv) :: advect_1d_constant => &
      CSL_advect_1d_constant
    procedure, pass(adv) :: delete => delete_CSL_1d_advector
  end type CSL_1d_advector
   




contains
  function new_CSL_1d_advector( &
    interp, &
    charac, &
    Npts, &
    eta_min, &
    eta_max, &
    eta_coords, &
    bc_type, &
    process_outside_point) &
    result(adv)      
    type(CSL_1d_advector), pointer :: adv
    class(sll_interpolator_1d_base), pointer :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in), optional :: eta_min
    sll_real64, intent(in), optional :: eta_max
    sll_real64, dimension(:), pointer, optional :: eta_coords
    sll_int32, intent(in), optional :: bc_type
    procedure(signature_process_outside_point_1d), optional    :: &
      process_outside_point
    sll_int32 :: ierr
    
    SLL_ALLOCATE(adv,ierr)
        
    call initialize_CSL_1d_advector(&
      adv, &
      interp, &
      charac, &
      Npts, &
      eta_min, &
      eta_max, &
      eta_coords, &   
      bc_type, &
      process_outside_point)
    
  end function  new_CSL_1d_advector


  subroutine initialize_CSL_1d_advector(&
    adv, &
    interp, &
    charac, &
    Npts, &
    eta_min, &
    eta_max, &
    eta_coords, &   
    bc_type, &
    process_outside_point)
    class(CSL_1d_advector), intent(inout) :: adv
    class(sll_interpolator_1d_base), pointer :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in), optional :: eta_min
    sll_real64, intent(in), optional :: eta_max
    sll_real64, dimension(:), pointer, optional :: eta_coords
    sll_int32, intent(in), optional :: bc_type
    procedure(signature_process_outside_point_1d), optional    :: &
      process_outside_point
    sll_int32 :: ierr
    sll_int32 :: i
    sll_real64 :: delta_eta
    
    
    adv%Npts = Npts
    adv%interp => interp
    adv%charac => charac
    SLL_ALLOCATE(adv%eta_coords(Npts),ierr)
    SLL_ALLOCATE(adv%eta_coords_unit(Npts),ierr)
    SLL_ALLOCATE(adv%eta_coords_unit_back(Npts),ierr)
    SLL_ALLOCATE(adv%buf1d(Npts),ierr)
    SLL_ALLOCATE(adv%buf1d_out(Npts),ierr)

    SLL_ALLOCATE(adv%charac_feet(Npts),ierr)
    SLL_ALLOCATE(adv%charac_feet_i(Npts),ierr)

    if(present(eta_min).and.present(eta_max))then
      if(present(eta_coords))then
        print *,'#provide either eta_coords or eta_min and eta_max'
        print *,'#and not both in subroutine initialize_BSL_1d_advector'
        stop
      else
        delta_eta = (eta_max-eta_min)/real(Npts-1,f64)
        do i=1,Npts
          adv%eta_coords(i) = eta_min+real(i-1,f64)*delta_eta
        enddo
      endif
    else if(present(eta_coords))then
      if(size(eta_coords)<Npts)then
        print *,'#bad size for eta_coords in initialize_BSL_1d_advector'
        stop
      else
        adv%eta_coords(1:Npts) = eta_coords(1:Npts)
      endif     
    else
      print *,'#Warning, we assume eta_min = 0._f64 eta_max = 1._f64'
      delta_eta = 1._f64/real(Npts-1,f64)
      do i=1,Npts
          adv%eta_coords(i) = real(i-1,f64)*delta_eta
      enddo                      
    endif



    if(present(process_outside_point)) then
      adv%process_outside_point => process_outside_point
      adv%bc_type = SLL_USER_DEFINED
    else if(.not.(present(bc_type))) then
      print *,'#provide boundary condition'
      print *,'#bc_type or process_outside_point function'
      print *,'#in initialize_CSL_1d_advector'
      stop
    else
      adv%bc_type = bc_type
      select case (bc_type)
        case (SLL_PERIODIC)
          adv%process_outside_point => process_outside_point_periodic                    
        case (SLL_SET_TO_LIMIT)
          adv%process_outside_point => process_outside_point_set_to_limit        
        case default
          print *,'#bad value of boundary condition'
          print *,'#in initialize_CSL_1d_advector'
          stop
        end select
    endif
    
    if((present(process_outside_point)).and.(present(bc_type)))then
      print *,'#provide either process_outside_point or bc_type'
      print *,'#and not both'
      print *,'#in initialize_CSL_1d_advector'
      stop
    endif
    
    
    adv%eta_coords_unit(1:Npts) = &
      (adv%eta_coords(1:Npts)-adv%eta_coords(1)) &
      /(adv%eta_coords(Npts)-adv%eta_coords(1))
          
  end subroutine initialize_CSL_1d_advector

  subroutine CSL_advect_1d(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(CSL_1d_advector) :: adv
    sll_real64, dimension(:), intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    sll_int32 :: Npts
    sll_real64 :: mean
    sll_real64 :: x1
    sll_int32 :: i
    sll_real64 :: eta_min
    sll_real64 :: eta_max
    
    Npts = adv%Npts
    
    eta_min = adv%eta_coords(1)
    eta_max = adv%eta_coords(Npts)
    
    
    
    
    call adv%charac%compute_characteristics( &
      A, &
      dt, &
      adv%eta_coords, &
      adv%charac_feet)


    do i=1,Npts
      if((adv%charac_feet(i)<=eta_min).or.(adv%charac_feet(i)>=eta_max))then
        adv%charac_feet_i(i) = &
          adv%process_outside_point(adv%charac_feet(i),eta_min,eta_max)
      else
        adv%charac_feet_i(i) =  adv%charac_feet(i)   
      endif      
    enddo


    
    adv%buf1d(1:Npts-1) = input(1:Npts-1)

!    do i=1,Npts
!      !x1 = adv%eta_coords(1)+real(i-1,f64)*(adv%eta_coords(Npts)-adv%eta_coords(1))/real(Npts-1,f64)
!      x1 = adv%eta_coords_unit(i)
!      adv%buf1d(i) = cos(2._f64*sll_pi*x1)
!    enddo
    
    
    call function_to_primitive(adv%buf1d,adv%eta_coords_unit,Npts-1,mean)

    !adv%buf1d(1:Npts) = adv%buf1d(1:Npts)*(adv%eta_coords(Npts)-adv%eta_coords(1))
    !print *,'#mean=',mean
!
!    do i=1,Npts
!      !print *,i,input(i),adv%buf1d(i),adv%eta_coords_unit(i),adv%eta_coords(i)
!      print *,adv%eta_coords_unit(i),cos(2._f64*sll_pi*adv%eta_coords_unit(i)),adv%buf1d(i), &
!      sin(2._f64*sll_pi*adv%eta_coords_unit(i))/(2._f64*sll_pi)
!      !-sin(adv%eta_coords(i)-adv%eta_coords(1))/(2._f64*sll_pi)
!    enddo

    !do i=1,Npts
    !  print *,i,adv%eta_coords(i),adv%charac_feet(i)
    !enddo

    
!    stop
    
    
!    call adv%interp%compute_interpolants( &
!      input, &
!      adv%eta1_coords, &
!      adv%Npts1, &
!      adv%eta2_coords, &
!      adv%Npts2 )

    adv%buf1d_out = adv%interp%interpolate_array( &
      Npts, &
      adv%buf1d, &
      adv%charac_feet_i)      

    !adv%buf1d_out(1:Npts) = adv%buf1d_out(1:Npts)/(adv%eta_coords(Npts)-adv%eta_coords(1))


    adv%eta_coords_unit_back(1:Npts) = &
      (adv%charac_feet(1:Npts)-adv%charac_feet(1)) &
      /(adv%charac_feet(Npts)-adv%charac_feet(1))

    
    
    call primitive_to_function( &
      adv%buf1d_out, &
      adv%eta_coords_unit, &
      adv%eta_coords_unit_back, &
      Npts-1, &
      mean)
    
    output(1:Npts-1) = adv%buf1d_out(1:Npts-1)
    
    output(Npts) = output(1) !warning only valid in periodic case
          
  end subroutine CSL_advect_1d


  subroutine CSL_advect_1d_constant(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(CSL_1d_advector) :: adv
    sll_real64, intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    sll_real64, dimension(:), allocatable :: A1
    sll_int32 :: ierr
    
    !this version is not optimized
    
    SLL_ALLOCATE(A1(adv%Npts),ierr)
    
    A1 = A
    
    call adv%charac%compute_characteristics( &
      A1, &
      dt, &
      adv%eta_coords, &
      adv%charac_feet)

!    call adv%interp%compute_interpolants( &
!      input, &
!      adv%eta1_coords, &
!      adv%Npts1, &
!      adv%eta2_coords, &
!      adv%Npts2 )

    output = adv%interp%interpolate_array( &
      adv%Npts, &
      input, &
      adv%charac_feet)      

    SLL_DEALLOCATE_ARRAY(A1,ierr)

          
  end subroutine CSL_advect_1d_constant

  subroutine function_to_primitive(f,node_positions,N,M)
    sll_real64,dimension(:),intent(inout) :: f
    sll_real64,dimension(:),intent(in) :: node_positions
    sll_int32,intent(in):: N
    sll_real64,intent(out)::M
    sll_int32::i
    sll_real64::tmp,tmp2
        
    !from f compute the mean
    M=0._f64
    do i=1,N
      M=M+f(i)*(node_positions(i+1)-node_positions(i))
    enddo
    
    f(1)=(f(1)-M)*(node_positions(2)-node_positions(1))
    tmp=f(1)
    f(1)=0._f64
    do i=2,N!+1
      f(i)=(f(i)-M)*(node_positions(i+1)-node_positions(i))
      tmp2=f(i)
      f(i)=f(i-1)+tmp
      tmp=tmp2
    enddo    
    f(N+1)=f(N)+tmp
    
    
    !print *,f(1),f(N+1) 

  end subroutine function_to_primitive




  subroutine primitive_to_function( &
    f, &
    node_positions, &
    node_positions_back, &
    N, &
    M)
    sll_real64,dimension(:),intent(inout) :: f
    sll_real64,dimension(:),intent(in) :: node_positions
    sll_real64,dimension(:),intent(in) :: node_positions_back
    sll_int32,intent(in):: N
    sll_real64,intent(in)::M
    sll_int32::i
    sll_real64::tmp!,tmp2
    
    tmp=f(1)
    do i=1,N-1
      f(i)=f(i+1)-f(i)+M*(node_positions_back(i+1)-node_positions_back(i))
    enddo
    f(N)=tmp-f(N)+M*(node_positions_back(N+1)-node_positions_back(N))


    !from mean compute f
    do i=1,N
      f(i)=f(i)/(node_positions(i+1)-node_positions(i))
    enddo
    !f(N+1) = f(1)
  end subroutine primitive_to_function


  subroutine delete_CSL_1d_advector( adv )
    class(CSL_1d_advector), intent(inout) :: adv
    sll_int32 :: ierr
    SLL_DEALLOCATE(adv%eta_coords,ierr)
    SLL_DEALLOCATE(adv%eta_coords_unit,ierr)
    SLL_DEALLOCATE(adv%eta_coords_unit_back,ierr)
    SLL_DEALLOCATE(adv%charac_feet,ierr)
    SLL_DEALLOCATE(adv%charac_feet_i,ierr)
    SLL_DEALLOCATE(adv%buf1d,ierr)
    SLL_DEALLOCATE(adv%buf1d_out,ierr)
  end subroutine delete_CSL_1d_advector




end module sll_module_advection_1d_CSL
