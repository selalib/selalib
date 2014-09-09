module  compute_characteristic_euler_2d_hex
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
  use sll_constants
  use sll_module_characteristics_2d_base
  implicit none

  type,extends(sll_characteristics_2d_base) :: euler_2d_hex_charac_computer
  sll_int32                               :: Num_cells ! num_cells is the step used for a hexagonal mesh
  sll_real64                              :: radius    ! the mesh is determined by the radius ( in contrast to eta in a cartesian mesh)  

  ! since there are 6 boundary conditions
  procedure(signature_process_outside_point), pointer, nopass :: &
       process_outside_point1
  procedure(signature_process_outside_point), pointer, nopass :: &
       process_outside_point2
  procedure(signature_process_outside_point), pointer, nopass :: &
       process_outside_point3
  procedure(signature_process_outside_point), pointer, nopass :: &
       process_outside_point4
  procedure(signature_process_outside_point), pointer, nopass :: &
       process_outside_point5
  procedure(signature_process_outside_point), pointer, nopass :: &
       process_outside_point6

   contains
     
     procedure, pass(charac) :: initialize => &
          initialize_euler_2d_hex_charac
     procedure, pass(charac) :: compute_characteristics => &
          compute_euler_2d_hex_charac
     
  end type euler_2d_hex_charac_computer


contains


  function new_euler_2d_hex_charac(&
      Npts,
      bc_type_1, &
      bc_type_2, &
      bc_type_3, &
      bc_type_4, &
      bc_type_5, &
      bc_type_6, &
      radius, &
      process_outside_point1, &
      process_outside_point2, &
      process_outside_point3, &
      process_outside_point4, &
      process_outside_point5, &
      process_outside_point6) &
      result(charac)
      
    type(euler_2d_hex_charac_computer),pointer :: charac
    sll_int32, intent(in) :: Num_cells
    sll_int32, intent(in), optional :: bc_type_1
    sll_int32, intent(in), optional :: bc_type_2
    sll_int32, intent(in), optional :: bc_type_3
    sll_int32, intent(in), optional :: bc_type_4
    sll_int32, intent(in), optional :: bc_type_5
    sll_int32, intent(in), optional :: bc_type_6
    sll_real64, intent(in), optional  :: radius
    procedure(signature_process_outside_point), optional  :: &
      process_outside_point1
    procedure(signature_process_outside_point), optional  :: &
      process_outside_point2
    procedure(signature_process_outside_point), optional  :: &
      process_outside_point3
    procedure(signature_process_outside_point), optional  :: &
      process_outside_point4
    procedure(signature_process_outside_point), optional  :: &
      process_outside_point5
    procedure(signature_process_outside_point), optional  :: &
      process_outside_point6
    sll_int32 :: ierr
      
    SLL_ALLOCATE(charac,ierr)

    call initialize_euler_2d_hex_charac(&
      charac, &
      Num_cells, &
      bc_type_1, &
      bc_type_2, &
      bc_type_3, &
      bc_type_4, &
      bc_type_5, &
      bc_type_6, &
      radius, &
      process_outside_point1, &
      process_outside_point2, &
      process_outside_point3, &
      process_outside_point4, &
      process_outside_point5, &
      process_outside_point6)
    
  end function new_euler_2d_hex_charac



  subroutine initialize_euler_2d_hex_charac(&
      charac, &
      Num_cells, &
      bc_type_1, &
      bc_type_2, &
      bc_type_3, &
      bc_type_4, &
      bc_type_5, &
      bc_type_6, &
      radius, &
      process_outside_point1, &
      process_outside_point2, &
      process_outside_point3, &
      process_outside_point4, &
      process_outside_point5, &
      process_outside_point6)
      
    class(euler_2d_hex_charac_computer) :: charac
    sll_int32, intent(in) :: Num_cells
    sll_int32, intent(in), optional :: bc_type_1
    sll_int32, intent(in), optional :: bc_type_2
    sll_int32, intent(in), optional :: bc_type_3
    sll_int32, intent(in), optional :: bc_type_4
    sll_int32, intent(in), optional :: bc_type_5
    sll_int32, intent(in), optional :: bc_type_6
    sll_real64, intent(in), optional  :: radius
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point1
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point2
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point3
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point4
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point5
    procedure(signature_process_outside_point), optional    :: &
      process_outside_point6


    charac%Num_cells = Num_cells
    
    
    if(present(radius))then
      charac%eta1_min = eta1_min
    else
      charac%eta1_min = 0._f64
    endif
    
    
    if(present(process_outside_point1)) then
      charac%process_outside_point1 => process_outside_point1
    else if(.not.(present(bc_type_1))) then
      print *,'#provide boundary condition'
      print *,'#bc_type_1 or process_outside_point1 function'
      print *,'#in initialize_euler_2d_hex_charac'
      stop
    else
      select case (bc_type_1)
        case (SLL_PERIODIC)
          charac%process_outside_point1 => process_outside_point_periodic          
        case (SLL_SET_TO_LIMIT)
          charac%process_outside_point1 => process_outside_point_set_to_limit        
        case default
          print *,'#bad value of boundary condition'
          print *,'#in initialize_euler_2d_hex_charac_computer'
          stop
        end select
    endif
    
    if((present(process_outside_point1)).and.(present(bc_type_1)))then
      print *,'#provide either process_outside_point1 or bc_type_1'
      print *,'#and not both'
      print *,'#in initialize_euler_2d_hex_charac_computer'
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
          print *,'#in initialize_euler_2d_hex_charac_computer'
          stop
        end select
    endif

    if((present(process_outside_point2)).and.(present(bc_type_2)))then
      print *,'#provide either process_outside_point2 or bc_type_2'
      print *,'#and not both'
      print *,'#in initialize_euler_2d_hex_charac_computer'
      stop
    endif
    
    
    
  end subroutine initialize_euler_2d_hex_charac

  subroutine compute_euler_2d_hex_charac( &
      charac, &
      A1, &
      A2, &
      dt, &
      input1, &
      input2, &
      output1, &
      output2)
            
    class(euler_2d_hex_charac_computer) :: charac
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64, intent(in) :: dt
    sll_real64, dimension(:), intent(in) ::  input1
    sll_real64, dimension(:), intent(in) ::  input2
    sll_real64, dimension(:,:), intent(out) :: output1
    sll_real64, dimension(:,:), intent(out) :: output2    
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: Num_cells
    sll_real64 :: radius
    sll_real64 :: step
    
    Num_cells = charac%Num_cells
    radius = charac%radius
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
    
    do j = - Num_cells,Num_cells
      do i = - Num_cells,Num_cells

         if ( i*j < 0 .and. abs(i) + abs(j) == Num_cells ) then
            ! these are not mesh points
         else ! these are mesh points
         
         call compute_characteristic_euler_2d_hex( &
       input1(i),input2(i),A1,A2,i,j,output1(i,j),output2(i,j),dt,step )

        !output1(i,j) = input1(i)-dt*A1(i,j) ! euler
        !output2(i,j) = input2(j)-dt*A2(i,j) ! euler 

         ! in the case the output is outside : find through which one
         ! of the 6 sides it 'came from'   
         ! then we use the appropriate procedure


        if((output1(i,j)<=eta1_min).or.(output1(i,j)>=eta1_max))then
          output1(i,j)=charac%process_outside_point1(output1(i,j),eta1_min,eta1_max)
        endif  

        if((output2(i,j)<=eta2_min).or.(output2(i,j)>=eta2_max))then
          output2(i,j)=charac%process_outside_point2(output2(i,j),eta2_min,eta2_max)          
        endif
      enddo
    enddo   
      
  end subroutine compute_euler_2d_hex_charac



  subroutine compute_characteristic_euler_2d_hex( &
       x1,x2,Uxn,Uyn,i,j,y1,y2,dt,step )

    sll_real64,dimension(:),intent(in):: Uxn, Uyn
    sll_real64, intent(in)  :: step ! no necessary here due to no newton resolution
    sll_real64, intent(in)  :: dt
    sll_real64, intent(in)  :: x1, x2 ! point of the characteristic at tn+1 
    sll_real64, intent(out) :: y1, y2 ! point of the characteristic at tn
    sll_int32, intent(in)   :: i,j
    
    y1 = x1 - dt*Uxn(i,j) 
    y2 = x2 - dt*Uyn(i,j) 	

  end subroutine compute_characteristic_euler_2d_hex


end module compute_characteristic_euler_2d_hex
