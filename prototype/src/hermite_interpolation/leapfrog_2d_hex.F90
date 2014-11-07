module module_characteristics_2d_leapfrog
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
  use sll_constants
  use sll_module_characteristics_2d_base
  implicit none

  type,extends(sll_characteristics_2d_base) :: leapfrog_2d_charac_computer
  sll_int32                               :: Npts1
  sll_int32                               :: Npts2
  sll_real64                              :: eta1_min   
  sll_real64                              :: eta1_max  
  sll_real64                              :: eta2_min   
  sll_real64                              :: eta2_max
  procedure(signature_process_outside_point), pointer, nopass :: &
       process_outside_point1
  procedure(signature_process_outside_point), pointer, nopass :: &
       process_outside_point2

   contains
     
     procedure, pass(charac) :: initialize => &
          initialize_leapfrog_2d_charac
     procedure, pass(charac) :: compute_characteristics => &
          compute_leapfrog_2d_charac
     
  end type leapfrog_2d_charac_computer


contains


  function new_leapfrog_2d_charac(&
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
      
    type(leapfrog_2d_charac_computer),pointer :: charac
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_int32, intent(in), optional :: bc_type_1
    sll_int32, intent(in), optional :: bc_type_2
    sll_real64, intent(in), optional  :: eta1_min
    sll_real64, intent(in), optional  :: eta1_max
    sll_real64, intent(in), optional  :: eta2_min
    sll_real64, intent(in), optional  :: eta2_max
    procedure(signature_process_outside_point), optional  :: &
      process_outside_point1
    procedure(signature_process_outside_point), optional  :: &
      process_outside_point2
    sll_int32 :: ierr
      
    SLL_ALLOCATE(charac,ierr)

    call initialize_leapfrog_2d_charac(&
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
    
  end function new_leapfrog_2d_charac



  subroutine initialize_leapfrog_2d_charac(&
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
      
    class(leapfrog_2d_charac_computer) :: charac
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
      print *,'#in initialize_leapfrog_2d_charac'
      stop
    else
      select case (bc_type_1)
        case (SLL_PERIODIC)
          charac%process_outside_point1 => process_outside_point_periodic          
        case (SLL_SET_TO_LIMIT)
          charac%process_outside_point1 => process_outside_point_set_to_limit        
        case default
          print *,'#bad value of boundary condition'
          print *,'#in initialize_leapfrog_2d_charac_computer'
          stop
        end select
    endif
    
    if((present(process_outside_point1)).and.(present(bc_type_1)))then
      print *,'#provide either process_outside_point1 or bc_type_1'
      print *,'#and not both'
      print *,'#in initialize_leapfrog_2d_charac_computer'
      stop
    endif
    


    if(present(process_outside_point2)) then
      charac%process_outside_point2 => process_outside_point2
    else if(.not.(present(bc_type_2))) then
      print *,'#provide boundary condition'
      print *,'#bc_type_2 or process_outside_point1 function'
      stop
    else
      select case (bc_type_2)
        case (SLL_PERIODIC)
          charac%process_outside_point2 => process_outside_point_periodic          
        case (SLL_SET_TO_LIMIT)
          charac%process_outside_point2 => process_outside_point_set_to_limit        
        case default
          print *,'#bad value of boundary condition'
          print *,'#in initialize_leapfrog_2d_charac_computer'
          stop
        end select
    endif

    if((present(process_outside_point2)).and.(present(bc_type_2)))then
      print *,'#provide either process_outside_point2 or bc_type_2'
      print *,'#and not both'
      print *,'#in initialize_leapfrog_2d_charac_computer'
      stop
    endif
    
    
    
  end subroutine initialize_leapfrog_2d_charac

  subroutine compute_leapfrog_2d_charac( &
      charac, &
      A1, &
      A2, &
      dt, &
      input1, &
      input2, &
      output1, &
      output2)
            
    class(leapfrog_2d_charac_computer) :: charac
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
    sll_real64 :: dx
    sll_real64 :: dy
    
    Npts1 = charac%Npts1
    Npts2 = charac%Npts2
    eta1_min = charac%eta1_min
    eta1_max = charac%eta1_max
    eta2_min = charac%eta2_min
    eta2_max = charac%eta2_max
    dx = input1(2) - input1(1) 
    dy = input2(2) - input2(1) 

    SLL_ASSERT(size(A1,1)>=charac%Npts1)
    SLL_ASSERT(size(A1,2)>=charac%Npts2)
    SLL_ASSERT(size(A2,1)>=charac%Npts1)
    SLL_ASSERT(size(A2,2)>=charac%Npts2)
    SLL_ASSERT(size(input1)>=charac%Npts1)
    SLL_ASSERT(size(input2)>=charac%Npts2)
    SLL_ASSERT(size(output1,1)>=charac%Npts1)
    SLL_ASSERT(size(output1,2)>=charac%Npts2)
    SLL_ASSERT(size(output2,1)>=charac%Npts1)
    SLL_ASSERT(size(output2,2)>=charac%Npts2)
    
    do j=1,Npts2
      do i=1,Npts1

         call compute_characteristic_leapfrog( &
       input1(i),input2(i),A1,A2,i,j,output1(i,j),output2(i,j),dt, dx, dy )

        !output1(i,j) = input1(i)-dt*A1(i,j) ! euler

        if((output1(i,j)<=eta1_min).or.(output1(i,j)>=eta1_max))then
          output1(i,j)=charac%process_outside_point1(output1(i,j),eta1_min,eta1_max)
        endif  

        !output2(i,j) = input2(j)-dt*A2(i,j) ! euler 

        if((output2(i,j)<=eta2_min).or.(output2(i,j)>=eta2_max))then
          output2(i,j)=charac%process_outside_point2(output2(i,j),eta2_min,eta2_max)          
        endif
      enddo
    enddo   
      
  end subroutine compute_leapfrog_2d_charac



  subroutine compute_characteristic_leapfrog( &
       x1,x2,Uxn,Uyn,i,j,y1,y2,dt,dx,dy )

    sll_real64,dimension(:),intent(in):: Uxn, Uyn, Uxn_1, Uyn_1
    sll_real64, intent(in)  :: dx, dy
    sll_real64, intent(in)  :: dt
    sll_real64, intent(in)  :: x1, x2
    sll_real64, intent(out) :: y1, y2
    sll_int32, intent(in)   :: i, j
    sll_real64              :: Uxn1, Uyn1, d1x, d1y, d1j0, d1j1
    
    ! needs to be optimized

    d1x = dt * Uxn(i,j)
    d1y = dt * Uyn(i,j)

    dxUxn = (Uxn(i+1,j) - Uxn(i,j))/dx  
    dyUxn = (Uxn(i,j+1) - Uxn(i,j))/dy 
    dxUyn = (Uyn(i+1,j) - Uyn(i,j))/dx  
    dyUyn = (Uyn(i,j+1) - Uyn(i,j))/dy 

    dij0 = d1x - dt * ( d1x*dxUxn + d1y*dyUxn ) 
    dij1 = d1y - dt * ( d1y*dyUyn + d1x*dxUyn )

    y1 = x1 - 2._f64 * dij0
    y2 = x2 - 2._f64 * dij1
	  

  end subroutine compute_characteristic_leapfrog





  ! subroutine slb_compute_characteristic_adams_moulton2( &
  !      x1,x2,Uxn,Uyn,Uxn_1,Uyn_1,i,j,y1,y2,dt )

  !   sll_real64,dimension(:),intent(in):: Uxn, Uyn, Uxn_1, Uyn_1
  !   sll_real64, intent(in)  :: dt
  !   sll_real64, intent(in)  :: x1, x2
  !   sll_real64, intent(out) :: y1, y2
  !   sll_int32, intent(in)   :: i, j
  !   sll_real64              :: Uxn1, Uyn1, d1x, d1y,d1j0, d1j1
    
  !   ! needs to be optimized

  !   Uxn1 = 2._f64*Uxn(i,j) - Uxn_1(i,j)
  !   Uyn1 = 2._f64*Uyn(i,j) - Uyn_1(i,j)

  !   d1x = 0.5_f64*dt * (Uxn1 + Uxn_1(i,j))
  !   d1y = 0.5_f64*dt * (Uyn1 + Uyn_1(i,j))

  !   dij0 = d1x - 0.5_f64*dt *( d1x*dxUxn + d1y*dyUxn )
  !   dij1 = d1y - 0.5_f64*dt *( d1x*dxUyn + d1y*dyUyn )

  !   y1 = x1 - dij0
  !   y2 = x2 - dij1

  ! end subroutine slb_compute_characteristic_adams_moulton2

end module compute_characteristic
