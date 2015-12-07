module  sll_m_euler_2d_hex
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    compute_characteristic_adams2_2d_hex, &
    compute_characteristic_euler_2d_hex

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   type,extends(sll_characteristics_2d_base) :: euler_2d_hex_charac_computer
!   sll_int32                               :: Num_cells ! num_cells is the step used for a hexagonal mesh
!   sll_real64                              :: radius    ! the mesh is determined by the radius ( in contrast to eta in a cartesian mesh)  

!   ! since there are 6 boundary conditions
!   procedure(signature_process_outside_point), pointer, nopass :: &
!        process_outside_point1
!   procedure(signature_process_outside_point), pointer, nopass :: &
!        process_outside_point2
!   procedure(signature_process_outside_point), pointer, nopass :: &
!        process_outside_point3
!   procedure(signature_process_outside_point), pointer, nopass :: &
!        process_outside_point4
!   procedure(signature_process_outside_point), pointer, nopass :: &
!        process_outside_point5
!   procedure(signature_process_outside_point), pointer, nopass :: &
!        process_outside_point6

!    contains
     
!      procedure, pass(charac) :: initialize => &
!           initialize_euler_2d_hex_charac
!      procedure, pass(charac) :: compute_characteristics => &
!           compute_euler_2d_hex_charac
     
!   end type euler_2d_hex_charac_computer


 contains


!   function new_euler_2d_hex_charac(&
!       Npts,
!       bc_type_1, &
!       bc_type_2, &
!       bc_type_3, &
!       bc_type_4, &
!       bc_type_5, &
!       bc_type_6, &
!       radius, &
!       process_outside_point1, &
!       process_outside_point2, &
!       process_outside_point3, &
!       process_outside_point4, &
!       process_outside_point5, &
!       process_outside_point6) &
!       result(charac)
      
!     type(euler_2d_hex_charac_computer),pointer :: charac
!     sll_int32, intent(in) :: Num_cells
!     sll_int32, intent(in), optional :: bc_type_1
!     sll_int32, intent(in), optional :: bc_type_2
!     sll_int32, intent(in), optional :: bc_type_3
!     sll_int32, intent(in), optional :: bc_type_4
!     sll_int32, intent(in), optional :: bc_type_5
!     sll_int32, intent(in), optional :: bc_type_6
!     sll_real64, intent(in), optional  :: radius
!     procedure(signature_process_outside_point), optional  :: &
!       process_outside_point1
!     procedure(signature_process_outside_point), optional  :: &
!       process_outside_point2
!     procedure(signature_process_outside_point), optional  :: &
!       process_outside_point3
!     procedure(signature_process_outside_point), optional  :: &
!       process_outside_point4
!     procedure(signature_process_outside_point), optional  :: &
!       process_outside_point5
!     procedure(signature_process_outside_point), optional  :: &
!       process_outside_point6
!     sll_int32 :: ierr
      
!     SLL_ALLOCATE(charac,ierr)

!     call initialize_euler_2d_hex_charac(&
!       charac, &
!       Num_cells, &
!       bc_type_1, &
!       bc_type_2, &
!       bc_type_3, &
!       bc_type_4, &
!       bc_type_5, &
!       bc_type_6, &
!       radius, &
!       process_outside_point1, &
!       process_outside_point2, &
!       process_outside_point3, &
!       process_outside_point4, &
!       process_outside_point5, &
!       process_outside_point6)
    
!   end function new_euler_2d_hex_charac



!   subroutine initialize_euler_2d_hex_charac(&
!       charac, &
!       Num_cells, &
!       bc_type_1, &
!       bc_type_2, &
!       bc_type_3, &
!       bc_type_4, &
!       bc_type_5, &
!       bc_type_6, &
!       radius, &
!       process_outside_point1, &
!       process_outside_point2, &
!       process_outside_point3, &
!       process_outside_point4, &
!       process_outside_point5, &
!       process_outside_point6)
      
!     class(euler_2d_hex_charac_computer) :: charac
!     sll_int32, intent(in) :: Num_cells
!     sll_int32, intent(in), optional :: bc_type_1
!     sll_int32, intent(in), optional :: bc_type_2
!     sll_int32, intent(in), optional :: bc_type_3
!     sll_int32, intent(in), optional :: bc_type_4
!     sll_int32, intent(in), optional :: bc_type_5
!     sll_int32, intent(in), optional :: bc_type_6
!     sll_real64, intent(in), optional  :: radius
!     procedure(signature_process_outside_point), optional    :: &
!       process_outside_point1
!     procedure(signature_process_outside_point), optional    :: &
!       process_outside_point2
!     procedure(signature_process_outside_point), optional    :: &
!       process_outside_point3
!     procedure(signature_process_outside_point), optional    :: &
!       process_outside_point4
!     procedure(signature_process_outside_point), optional    :: &
!       process_outside_point5
!     procedure(signature_process_outside_point), optional    :: &
!       process_outside_point6


!     charac%Num_cells = Num_cells
    
    
!     if(present(radius))then
!       charac%eta1_min = eta1_min
!     else
!       charac%eta1_min = 0._f64
!     endif
    
    
!     if(present(process_outside_point1)) then
!       charac%process_outside_point1 => process_outside_point1
!     else if(.not.(present(bc_type_1))) then
!       print *,'#provide boundary condition'
!       print *,'#bc_type_1 or process_outside_point1 function'
!       print *,'#in initialize_euler_2d_hex_charac'
!       stop
!     else
!       select case (bc_type_1)
!         case (SLL_PERIODIC)
!           charac%process_outside_point1 => process_outside_point_periodic          
!         case (SLL_SET_TO_LIMIT)
!           charac%process_outside_point1 => process_outside_point_set_to_limit        
!         case default
!           print *,'#bad value of boundary condition'
!           print *,'#in initialize_euler_2d_hex_charac_computer'
!           stop
!         end select
!     endif
    
!     if((present(process_outside_point1)).and.(present(bc_type_1)))then
!       print *,'#provide either process_outside_point1 or bc_type_1'
!       print *,'#and not both'
!       print *,'#in initialize_euler_2d_hex_charac_computer'
!       stop
!     endif
    


!     if(present(process_outside_point2)) then
!       charac%process_outside_point2 => process_outside_point2
!     else if(.not.(present(bc_type_2))) then
!       print *,'#provide boundary condition'
!       print *,'#bc_type_2 or process_outside_point1 function'
!       stop
!     else
!       select case (bc_type_2)
!         case (SLL_PERIODIC)
!           charac%process_outside_point2 => process_outside_point_periodic          
!         case (SLL_SET_TO_LIMIT)
!           charac%process_outside_point2 => process_outside_point_set_to_limit        
!         case default
!           print *,'#bad value of boundary condition'
!           print *,'#in initialize_euler_2d_hex_charac_computer'
!           stop
!         end select
!     endif

!     if((present(process_outside_point2)).and.(present(bc_type_2)))then
!       print *,'#provide either process_outside_point2 or bc_type_2'
!       print *,'#and not both'
!       print *,'#in initialize_euler_2d_hex_charac_computer'
!       stop
!     endif
    
    
    
!   end subroutine initialize_euler_2d_hex_charac

!   subroutine compute_euler_2d_hex_charac( &
!       charac, &
!       A1, &
!       A2, &
!       dt, &
!       input1, &
!       input2, &
!       output1, &
!       output2)
            
!     class(euler_2d_hex_charac_computer) :: charac
!     sll_real64, dimension(:,:), intent(in) :: A1
!     sll_real64, dimension(:,:), intent(in) :: A2
!     sll_real64, intent(in) :: dt
!     sll_real64, dimension(:), intent(in) ::  input1
!     sll_real64, dimension(:), intent(in) ::  input2
!     sll_real64, dimension(:,:), intent(out) :: output1
!     sll_real64, dimension(:,:), intent(out) :: output2    
!     sll_int32 :: i
!     sll_int32 :: j
!     sll_int32 :: Num_cells
!     sll_real64 :: radius
!     sll_real64 :: step
    
!     Num_cells = charac%Num_cells
!     radius = charac%radius

    
!     do j = - Num_cells,Num_cells
!       do i = - Num_cells,Num_cells

!          if ( i*j < 0 .and. abs(i) + abs(j) == Num_cells ) then
!             ! these are not mesh points
!          else ! these are mesh points
         
!          call compute_characteristic_euler_2d_hex( &
!        input1(i),input2(i),A1,A2,i,j,output1(i,j),output2(i,j),dt,step )

!         !output1(i,j) = input1(i)-dt*A1(i,j) ! euler
!         !output2(i,j) = input2(j)-dt*A2(i,j) ! euler 

!          ! in the case the output is outside : find through which one
!          ! of the 6 sides it 'came from'   
!          ! then we use the appropriate procedure




!         if(  output1(i,j) output2(i,j) )then
!           output1(i,j)=charac%process_outside_point6(output2(i,j),radius)   
!           output2(i,j)=charac%process_outside_point6(output2(i,j),radius)          
!         endif

!       enddo
!     enddo   
      
!   end subroutine compute_euler_2d_hex_charac



  subroutine compute_characteristic_euler_2d_hex( x1,x2,uxn,uyn,i,y1,y2,dt)

    sll_real64,dimension(:),intent(in):: uxn, uyn
    sll_real64, intent(in)  :: dt
    sll_real64, intent(in)  :: x1, x2 ! point of the characteristic at tn+1 
    sll_real64, intent(out) :: y1, y2 ! point of the characteristic at tn
    sll_int32, intent(in)   :: i
    
    y1 = x1 - dt*uxn(i)
    y2 = x2 - dt*uyn(i)

  end subroutine compute_characteristic_euler_2d_hex

  ! subroutine compute_characteristic_verlet_2d_hex( z1,z2,uxn,uyn,dxux,dyux,dxuy,dyuy,i,zz1,zz2,dt, aire, mesh)

  !   type(sll_hex_mesh_2d), pointer    :: mesh
  !   sll_real64,dimension(:),intent(in):: uxn, uyn,dxux,dyux,dxuy,dyuy
  !   sll_real64, intent(in)  :: dt, aire
  !   sll_real64, intent(in)  :: z1, z2 ! point of the characteristic at tn+1 
  !   sll_real64, intent(out) :: zz1, zz2 ! point of the characteristic at tn
  !   sll_int32, intent(in)   :: i
  !   sll_real64              :: x1,x2,x3,y1,y2,y3,xx,ey,erreur,f
    	
  !   step = mesh%delta
    
  !   ! premier newton 

  !   erreur = 1._f64

  !   do while(erreur > 1.E-12) 

  !      f = 

  !      xx = 

  !      erreur = abs(f)

  !   enddo
    
  !   xx = z1 - uyn(n_j)
    
  !   ! deuxième newton

  !   zz2 = 


  !   ! dernière égalité

  !   call get_cell_vertices_index( xx, zz2, mesh, i1, i2, i3 )

  !   x1 = mesh%cartesian_coord(1,i1) 
  !   x2 = mesh%cartesian_coord(1,i2) 
  !   x3 = mesh%cartesian_coord(1,i3) 
  !   y1 = mesh%cartesian_coord(2,i1) 
  !   y2 = mesh%cartesian_coord(2,i2) 
  !   y3 = mesh%cartesian_coord(2,i3) 


  !   a2  = 0.5_f64/aire
  !   y3y = y3 - y
  !   y2y = y2 - y
  !   y1y = y1 - y
  !   x3x = x3 - x
  !   x2x = x2 - x
  !   x1x = x1 - x

  !   l1   = a2 * abs( x2x*y3y - x3x*y2y )    ! barycentric coordinates
  !   l2   = a2 * abs( x1x*y3y - x3x*y1y ) 
  !   l3   = 1._f64 - l1 - l2
  !   ey = l1*uxn(i1) + l2*uxn(i2) + l3*uxn(i3) 

  !   zz1 = z1 - ey * dt * 0.5_f64

  ! end subroutine compute_characteristic_verlet_2d_hex


  subroutine compute_characteristic_leapfrog_2d_hex( x1,x2,uxn,uyn,dxux,dyux,dxuy,dyuy,i,y1,y2,dt)

    sll_real64,dimension(:),intent(in):: uxn, uyn,dxux,dyux,dxuy,dyuy
    sll_real64, intent(in)  :: dt
    sll_real64, intent(in)  :: x1, x2 ! point of the characteristic at tn+1 
    sll_real64, intent(out) :: y1, y2 ! point of the characteristic at tn
    sll_int32, intent(in)   :: i
    sll_real64              :: d1x, d1y, dij0, dij1 
    	
    
    d1x = dt * uxn(i)
    d1y = dt * uyn(i)

    dij0 = d1x - dt *( d1x*dxux(i) + d1y*dyux(i) ) 
    dij1 = d1y - dt *( d1y*dyuy(i) + d1x*dxuy(i) )

    y1 = x1 - 2._f64*dij0
    y2 = x2 - 2._f64*dij1

  end subroutine compute_characteristic_leapfrog_2d_hex


  subroutine compute_characteristic_adams2_2d_hex( x1,x2,uxn,uyn,uxn_1,uyn_1,&
       dxuxn,dyuxn,dxuyn,dyuyn,i,y1,y2,dt)
    sll_real64,dimension(:),intent(in):: uxn, uyn, uxn_1, uyn_1
    sll_real64,dimension(:),intent(in):: dxuxn,dyuxn,dxuyn,dyuyn
    sll_real64, intent(in)  :: dt
    sll_real64, intent(in)  :: x1, x2 ! point of the characteristic at tn+1 
    sll_real64, intent(out) :: y1, y2 ! point of the characteristic at tn
    sll_int32, intent(in)   :: i
    sll_real64              :: d1x, d1y, dij0, dij1, uxn1, uyn1 
    
    uxn1 = 2._f64*uxn(i) - uxn_1(i)
    uyn1 = 2._f64*uyn(i) - uyn_1(i)

    d1x = 0.5_f64*dt * ( uxn1 + uxn(i) )
    d1y = 0.5_f64*dt * ( uyn1 + uyn(i) )

    dij0 = d1x - 0.5_f64*dt * ( d1x*dxuxn(i) + d1y*dyuxn(i) ) 
    dij1 = d1y - 0.5_f64*dt * ( d1x*dxuyn(i) + d1y*dyuyn(i) )

    y1 = x1 - dij0
    y2 = x2 - dij1

  end subroutine compute_characteristic_adams2_2d_hex

  subroutine compute_characteristic_adams3_2d_hex( x1,x2,uxn,uyn,uxn_1,uyn_1,&
       uxn_2,uyn_2,dxuxn,dyuxn,dxuyn,dyuyn,i,y1,y2,dt)
    sll_real64,dimension(:),intent(in):: uxn, uyn, uxn_1, uyn_1,uxn_2,uyn_2
    sll_real64,dimension(:),intent(in):: dxuxn,dyuxn,dxuyn,dyuyn
    sll_real64, intent(in)  :: dt
    sll_real64, intent(in)  :: x1, x2 ! point of the characteristic at tn+1 
    sll_real64, intent(out) :: y1, y2 ! point of the characteristic at tn
    sll_int32, intent(in)   :: i
    sll_real64              :: d1x, d1y, dij0, dij1, uxn1, uyn1, erreur 
    sll_real64              :: a, b, c, d, det, gx, gy, xn, yn, xn_1,yn_1
    sll_real64              :: uxn_loc, uyn_loc, uxn_1_loc, uyn_1_loc
    sll_real64              :: dxuxn_loc,dyuxn_loc,dxuyn_loc,dyuyn_loc
    sll_real64              :: dxuxn_1_loc,dyuxn_1_loc,dxuyn_1_loc,dyuyn_1_loc
    
#ifdef DEBUG
    sll_real64 :: dummy
    dummy = dxuxn(1)+dxuyn(1)+dyuxn(1)+dyuyn(1)
#endif

    
    uxn1 = 3._f64*uxn(i) - 3._f64*uxn_1(i) + uxn_2(i)
    uyn1 = 3._f64*uyn(i) - 3._f64*uyn_1(i) + uyn_2(i)

    ! Uxn1 = 0.5*x1*tan(tn+dt)
    ! Uyn1 = 0.5*x2*tan(tn+dt)

    d1x = dt * (5._f64*uxn1 + 8._f64*uxn_1(i) - uxn_2(i))/12._f64
    d1y = dt * (5._f64*uyn1 + 8._f64*uyn_1(i) - uyn_2(i))/12._f64
    
    erreur = 1._f64

    do while(erreur > 1.E-12) 

       xn = x1 - d1x 
       yn = x2 - d1y 

       ! interpolation de uxn(xn), uyn(yn), dxUxn,dyUxn,dxUyn,dyUyn
       ! à faire ici

       xn_1 = x1 - 2._f64*dt*(uxn1+2*uxn_loc) + 4._f64*d1x
       yn_1 = x2 - 2._f64*dt*(uyn1+2*uyn_loc) + 4._f64*d1y

       ! interpolation de uxn_1(xn_1), uyn_1(yn_1), dxUxn_1,dyUxn_1,dxUyn_1,dyUyn_1
       ! à faire ici

       a   = 1._f64 + dt*(2._f64*dxuxn_loc - dxuxn_1_loc)/3._f64
       b   =          dt*(2._f64*dyuxn_loc - dyuxn_1_loc)/3._f64
       c   =          dt*(2._f64*dxuyn_loc - dxuyn_1_loc)/3._f64
       d   = 1._f64 + dt*(2._f64*dyuyn_loc - dyuyn_1_loc)/3._f64

       det = a*d - b*c 

       a = a/det
       b = b/det
       c = c/det
       d = d/det

       gx = d1x - dt*( 8._f64*uxn_loc - uxn_1_loc + 5._f64*uxn1)/12._f64
       gy = d1y - dt*( 8._f64*uyn_loc - uyn_1_loc + 5._f64*uyn1)/12._f64

       dij0 = d1x - ( +d*gx - b*gy )
       dij1 = d1y - ( -c*gx + a*gy )

       erreur = abs(gx) + abs(gy)
       d1x = dij0
       d1y = dij1

    enddo

    y1 = x1 - d1x
    y2 = x2 - d1y

  end subroutine compute_characteristic_adams3_2d_hex

  subroutine compute_characteristic_adams4_2d_hex( x1,x2,uxn,uyn,uxn_1,uyn_1,&
       uxn_2,uyn_2,uxn_3,uyn_3,dxuxn,dyuxn,dxuyn,dyuyn,i,y1,y2,dt)
    sll_real64,dimension(:),intent(in):: uxn, uyn, uxn_1, uyn_1,uxn_2,uyn_2
    sll_real64,dimension(:),intent(in):: uxn_3,uyn_3
    sll_real64,dimension(:),intent(in):: dxuxn,dyuxn,dxuyn,dyuyn
    sll_real64, intent(in)  :: dt
    sll_real64, intent(in)  :: x1, x2 ! point of the characteristic at tn+1 
    sll_real64, intent(out) :: y1, y2 ! point of the characteristic at tn
    sll_int32, intent(in)   :: i
    sll_real64              :: d1x, d1y, dij0, dij1, uxn1, uyn1, erreur 
    sll_real64              :: a, b, c, d, det, gx, gy, xn, yn, xn_1,yn_1
    sll_real64              :: xn_2,yn_2
    sll_real64              :: uxn_loc, uyn_loc, uxn_1_loc, uyn_1_loc, uxn_2_loc, uyn_2_loc
    sll_real64              :: dxuxn_loc,dyuxn_loc,dxuyn_loc,dyuyn_loc
    sll_real64              :: dxuxn_1_loc,dyuxn_1_loc,dxuyn_1_loc,dyuyn_1_loc
    sll_real64              :: dxuxn_2_loc,dyuxn_2_loc,dxuyn_2_loc,dyuyn_2_loc
#ifdef DEBUG
    sll_real64 :: dummy
    dummy = dxuxn(1)+dxuyn(1)+dyuxn(1)+dyuyn(1)+uyn_3(1)
#endif


    uxn1 = 4._f64*uxn(i) - 6._f64*uxn_1(i) + 4._f64*uxn_2(i) - uxn_3(i)
    uyn1 = 4._f64*uyn(i) - 6._f64*uyn_1(i) + 4._f64*uxn_2(i) - uxn_3(i)

    d1x = dt * (9._f64*uxn1 + 19._f64*uxn(i) - 5._f64*uxn_1(i) + uxn_2(i))/24._f64
    d1y = dt * (9._f64*uyn1 + 19._f64*uyn(i) - 5._f64*uyn_1(i) + uyn_2(i))/24._f64

    erreur = 1._f64


    do while(erreur > 1.e-12) 

       xn = x1 - d1x 
       yn = y2 - d1y 
       ! interpolation de uxn(xn), uyn(yn), dxUxn,dyUxn,dxUyn,dyUyn
       ! à faire ici

       xn_1 = x1 + 4._f64*d1x - 2._f64*dt*(2._f64*uxn_loc+uxn1) 
       yn_1 = y2 + 4._f64*d1y - 2._f64*dt*(2._f64*uyn_loc+uyn1) 
       ! interpolation de uxn_1(xn_1), uyn_1(yn_1), dxUxn_1,dyUxn_1,dxUyn_1,dyUyn_1

       ! xn_2 = x1 + 9.*d1x - 3.*dt*(3.*uxn + uxn1)
       ! yn_2 = y2 + 9.*d1y - 3.*dt*(3.*uyn + uyn1)

       xn_2 = x1 + 27._f64*d1x - dt*(18._f64*uxn_loc + 12._f64*uxn1)
       yn_2 = x2 + 27._f64*d1y - dt*(18._f64*uyn_loc + 12._f64*uyn1)

       ! interpolation de uxn_2(xn_2), uyn_2(yn_2), dxUxn_2,dyUxn_2,dxUyn_2,dyUyn_2

       a   = 1._f64 + dt*(19._f64*dxUxn_loc + 20._f64*dxUxn_1_loc &
            - 9._f64*dxUxn_2_loc)/24._f64
       b   =          dt*(19._f64*dyUxn_loc + 20._f64*dyUxn_1_loc &
            - 9._f64*dyUxn_2_loc)/24._f64
       c   =          dt*(19._f64*dxUyn_loc + 20._f64*dxUyn_1_loc &
            - 9._f64*dxUyn_2_loc)/24._f64
       d   = 1._f64 + dt*(19._f64*dyUyn_loc + 20._f64*dyUyn_1_loc &
            - 9._f64*dyUyn_2_loc)/24._f64

       det = a*d - b*c

       a = a/det
       b = b/det
       c = c/det
       d = d/det

       gx = d1x - dt*(9.*uxn1 + 19.*uxn_loc - 5.*uxn_1_loc  + uxn_2_loc )/24._f64
       gy = d1y - dt*(9.*uyn1 + 19.*uyn_loc - 5.*uyn_1_loc  + uyn_2_loc )/24._f64

       dij0 = d1x - ( +d*gx - b*gy ) 
       dij1 = d1y - ( -c*gx + a*gy )

       erreur = abs(gx) + abs(gy)
       d1x = dij0
       d1y = dij1

    enddo

    y1 = x1 - d1x
    y2 = x2 - d1y

  end subroutine compute_characteristic_adams4_2d_hex

end module sll_m_euler_2d_hex
