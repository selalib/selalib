program test_hex_hermite_circ

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_m_constants
  use sll_m_interpolation_hex_hermite
  !use sll_m_euler_2d_hex

  implicit none

  sll_real64, dimension(:,:), allocatable :: deriv

  sll_int32    :: num_cells, n_points, n_triangle, n_edge
  sll_int32    :: i
  sll_int32    :: num_method = 15
  !character(len = 5) ::name_test = "gauss"
  sll_int32    :: nloops,ierr, EXTRA_TABLES = 1 ! put 1 for num_method = 15
  ! initial distribution
  !sll_real64   :: r_min
  sll_real64   :: gauss_x2
  sll_real64   :: gauss_x1
  sll_real64   :: gauss_sig
  sll_real64   :: gauss_amp
  sll_real64,dimension(:),allocatable :: x1
  sll_real64,dimension(:),allocatable :: x2
  sll_real64,dimension(:),allocatable :: f_init
  sll_real64,dimension(:),allocatable :: x1_char
  sll_real64,dimension(:),allocatable :: x2_char
  ! distribution at time n
  sll_real64,dimension(:),allocatable :: f_tn, center_values_tn, edge_values_tn
  ! distribution at time n + 1
  sll_real64,dimension(:),allocatable :: f_tn1,center_values_tn1,edge_values_tn1
  ! exact distribution 
  sll_real64,dimension(:),allocatable :: f_sol
  sll_real64,dimension(:),allocatable :: phi, uxn, uyn,dxuxn,dyuxn,dxuyn,dyuyn
  sll_real64   :: norm2_error, norm2_sol, norm_infinite
  sll_real64   :: norm2_error_center, norm2_error_edge, norm2_error_pt
  sll_real64   :: norm2_sol_center, norm2_sol_edge, norm2_sol_pt
 
  sll_real64   :: dt
  sll_real64   :: tmax
  sll_real64   :: t
  
  ! timer variables
  sll_real64   :: t_init, t_end
  !others
  sll_int32    :: p = 6 !-> degree of the approximation for the derivative 
  sll_real64   :: step , aire, h1, h2, f_min, x ,y
  sll_real64   :: center_mesh_x1, center_mesh_x2, radius
  sll_real64   :: xx, yy!,z, x1_temp,t1,t2,t3,t4 
  sll_real64   :: cfl
  ! character(len = 4) :: number
  logical      :: inside
  type(sll_hex_mesh_2d), pointer :: mesh
  !character(len = 50) :: filename
  !character(len = 50) :: filename2
  !character(len = 4)  :: filenum

  center_mesh_x1 = 0._f64
  center_mesh_x2 = 0._f64

  radius = 8._f64
  !r_min  = 0._f64   ! beware there are some restrictions to respect 


  call print_method(num_method)

  open(unit = 33, file="hex_errors.txt", action="write", status="replace")

  write(33,*) 

  do num_cells = 20,160,20 ! -> loop on the size of the mesh 
  
     
     ! finding the hexagone corresponding to r_min

     !k_min = int(r_min/radius*real(num_cells)+1e-6)
     !if (k_min>0) then
     !   n_min = 1+3*k_min*(k_min-1)
     !endif

     !*********************************************************
     !             allocation
     !*********************************************************
     
     n_points   = 1 + 3 * num_cells * (num_cells + 1)  
     n_triangle = 6 * num_cells * num_cells
     n_edge    = 3 * num_cells * ( 3 * num_cells + 1 )

     step = radius / real(num_cells,f64)
     aire = step**2*sqrt(3._f64)*0.25_f64

     allocate( deriv(1:6,n_points) )

     SLL_ALLOCATE(f_init( n_points),ierr)
     SLL_ALLOCATE(f_tn( n_points),ierr)
     SLL_ALLOCATE(f_tn1( n_points ),ierr)
     SLL_ALLOCATE(phi( n_points ),ierr)
     SLL_ALLOCATE(uxn( n_points ),ierr)
     SLL_ALLOCATE(uyn( n_points ),ierr)
     SLL_ALLOCATE(dxuxn( n_points ),ierr)
     SLL_ALLOCATE(dyuxn( n_points ),ierr)
     SLL_ALLOCATE(dxuyn( n_points ),ierr)
     SLL_ALLOCATE(dyuyn( n_points ),ierr)
     SLL_ALLOCATE(center_values_tn ( n_triangle),ierr)
     SLL_ALLOCATE(center_values_tn1( n_triangle),ierr)
     SLL_ALLOCATE(edge_values_tn ( n_edge),ierr)
     SLL_ALLOCATE(edge_values_tn1( n_edge),ierr)
     SLL_ALLOCATE(f_sol( n_points ),ierr)
     SLL_ALLOCATE(x1( n_points),ierr)
     SLL_ALLOCATE(x2( n_points),ierr)
     SLL_ALLOCATE(x1_char( n_points),ierr)
     SLL_ALLOCATE(x2_char( n_points),ierr)

     !*********************************************************
     !                  Mesh initialization   
     !*********************************************************
     
     mesh => new_hex_mesh_2d( num_cells, center_mesh_x1, center_mesh_x2, radius=radius, EXTRA_TABLES = EXTRA_TABLES ) 
     gauss_x1  = 2._f64
     gauss_x2  = 2._f64
     gauss_sig = 1._f64/( 2._f64 * sqrt(2._f64)) 
     gauss_amp = 1.0_f64

     do i = 1, n_points

        x1(i) = mesh%cartesian_coord(1,i)
        x2(i) = mesh%cartesian_coord(2,i)
        
        
        f_init(i) = gauss_amp*exp(-0.5_f64* &
             ((x1(i)-gauss_x1)**2 + (x2(i)-gauss_x2)**2)/ gauss_sig**2 )
        f_tn(i) = f_init(i)


     end do

     
     if ( num_method == 10 ) then

        do i = 1, n_triangle

           x = mesh%center_cartesian_coord(1,i)
           y = mesh%center_cartesian_coord(2,i)

           center_values_tn(i) = gauss_amp*exp(-0.5_f64* &
             ((x-gauss_x1)**2 + (y-gauss_x2)**2)/ gauss_sig**2 )

        enddo
     endif
     
     if ( num_method == 15 ) then

        do i = 1, n_edge

           x = mesh%edge_center_cartesian_coord(1,i)
           y = mesh%edge_center_cartesian_coord(2,i)

           edge_values_tn(i) = gauss_amp*exp(-0.5_f64* &
             ((x-gauss_x1)**2 + (y-gauss_x2)**2)/ gauss_sig**2 )

        enddo

     endif

     close(11)

     tmax  = 0.3_f64!3._f64
     dt    = 0.01_f64*20._f64 / real(num_cells,f64)  
     t     = 0._f64
     cfl   = radius * dt / ( radius / real(num_cells,f64)  )

     !*********************************************************
     !              Computing characteristics
     !*********************************************************

     x1_char(:) = x1(:)*cos(dt) - x2(:)*sin(dt)
     x2_char(:) = x1(:)*sin(dt) + x2(:)*cos(dt)


     !*********************************************************
     !                          Time loop
     !*********************************************************

     nloops = 0

     call cpu_time(t_init)

     do while (t .lt. tmax)!dt)!tmax)!

        norm2_error = 0._f64 !Error variables
        norm2_sol   = 0._f64 !Error variables
        norm_infinite = 0._f64
        f_min = 0._f64

        nloops = nloops + 1

        call  der_finite_difference( f_tn, p, step, mesh, deriv)

        t = t + dt
        !*********************************************************
        ! computation of the value at the center of the triangles
        !*********************************************************

        if ( num_method == 10 ) then 

           
           !*********************************************************
           !computation of the value at the center of the triangles
           !*********************************************************
           
          norm2_error_center = 0._f64
          norm2_sol_center = 0._f64
           
           do i=1, n_triangle 

              x = mesh%center_cartesian_coord(1,i)
              y = mesh%center_cartesian_coord(2,i)

              xx = x*cos(dt) - y*sin(dt);
              yy = x*sin(dt) + y*cos(dt);

              inside = .true.

              h1 =  xx/sqrt(3.0_f64) + yy
              h2 = -xx/sqrt(3.0_f64) + yy 

              if ( abs(h1) >  radius-1e-15 .or. abs(h2) >  radius-1e-15 ) inside = .false.
              if ( abs(xx) > (radius-1e-15)*sqrt(3._f64)*0.5_f64  ) inside = .false.

              if ( inside ) then
                 call hermite_interpolation(i, xx, yy, f_tn, center_values_tn,&
                      edge_values_tn, center_values_tn1, mesh, deriv, aire,& 
                 num_method)
              else 
                 center_values_tn1(i) = 0._f64 ! dirichlet boundary condition
              endif
              

              xx = x*cos(t) - y*sin(t);
              yy = x*sin(t) + y*cos(t);

              ! computation of the following values :
              ! norm L2 of the distribution function and its minimum
              ! norm infinit & norme L2 of the error 
                 norm2_sol_center = norm2_sol_center + &
                      abs( gauss_amp*exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64))**2

                 norm2_error_center = norm2_error_center + &
                      abs( gauss_amp*exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64) - center_values_tn1(i) )**2

                 if ( abs( gauss_amp*exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64) - center_values_tn1(i)) >  norm_infinite )&
                      norm_infinite = abs( gauss_amp*exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64) - center_values_tn1(i))

                 if ( center_values_tn1(i) < f_min ) f_min = center_values_tn1(i)

           enddo

        endif

        
        if ( num_method == 15 ) then 


           !*********************************************************
           !        computation of the value at the center
           !            of the edge of the triangles
           !*********************************************************

           norm2_error_edge = 0._f64
           norm2_sol_edge = 0._f64

           do i=1, n_edge 

              x = mesh%edge_center_cartesian_coord(1,i)
              y = mesh%edge_center_cartesian_coord(2,i)

              xx = x*cos(dt) - y*sin(dt);
              yy = x*sin(dt) + y*cos(dt);

              !             INTERPOLATION
              inside = .true.

              h1 =  xx/sqrt(3.0_f64) + yy
              h2 = -xx/sqrt(3.0_f64) + yy 
              
              if ( abs(h1) >  radius-1e-15 .or. abs(h2) >  radius-1e-15 ) inside = .false.
              if ( abs(xx) > (radius-1e-15)*sqrt(3._f64)*0.5_f64  ) inside = .false.


              if ( inside ) then
                 call hermite_interpolation(i, xx, yy, f_tn, center_values_tn,&
                      edge_values_tn, edge_values_tn1, mesh, deriv, aire,& 
                 num_method)
              else 
                 edge_values_tn1(i) = 0._f64 ! dirichlet boundary condition
              endif


              xx = x*cos(t) - y*sin(t);
              yy = x*sin(t) + y*cos(t);

              ! computation of the following values :
              ! norm L2 of the distribution function and its minimum
              ! norm infinit & norme L2 of the error 
                 norm2_sol_edge = norm2_sol_edge + &
                      abs( gauss_amp*exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64))**2

                 norm2_error_edge = norm2_error_edge + &
                      abs( gauss_amp*exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64) - edge_values_tn1(i) )**2

                 if ( abs( gauss_amp*exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64) - edge_values_tn1(i)) >  norm_infinite )&
                      norm_infinite = abs( gauss_amp*exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64) - edge_values_tn1(i))

                 if ( edge_values_tn1(i) < f_min ) f_min = edge_values_tn1(i)


           enddo

        endif


        !*********************************************************
        !       computation of the value at the mesh points
        !***********************************

        norm2_error_pt = 0._f64
        norm2_sol_pt = 0._f64

        do i=1, n_points  

           x = mesh%cartesian_coord(1,i)
           y = mesh%cartesian_coord(2,i)

           xx = x*cos(dt) - y*sin(dt);
           yy = x*sin(dt) + y*cos(dt);

        !*********************************************************
        !                INTERPOLATION
        !*********************************************************
           ! computation of the interpolation at each point
           ! F(t_(n+1),x_i,v_j) = F (t_n,X(t_n),V(t_n))

           inside = .true.

           h1 =  xx/sqrt(3._f64) + yy
           h2 = -xx/sqrt(3._f64) + yy 
           
           if ( abs(h1) >  radius-1e-14 .or. abs(h2) >  radius-1e-14 ) inside = .false.
           if ( abs(xx) >  (radius-1e-14)*sqrt(3._f64)*0.5_f64  ) inside = .false.

           if ( inside ) then

              call hermite_interpolation(i, xx, yy, f_tn, center_values_tn,&
                   edge_values_tn, f_tn1, mesh, deriv, aire,& 
                   num_method)
           else 

              f_tn1(i) = 0._f64 ! dirichlet boundary condition

           endif

           ! ******************************************************
           ! Analytical value    (-> in order to compute the error)
           ! ******************************************************

           ! Computing the characteristics' root for the exact solution
           ! in order to compute the l2 norm of the error
           
              x = x1(i)*cos(t) - x2(i)*sin(t);
              y = x1(i)*sin(t) + x2(i)*cos(t);
           
           

              f_sol(i) = gauss_amp*exp(-((x-gauss_x1)**2+(y-gauss_x2)**2)/gauss_sig**2/2._f64) 
           
           
           ! computation of the following values :
           ! norm L2 of the distribution function and its minimum
           ! norm infinit & norme L2 of the error 
           
           norm2_sol_pt = norm2_sol_pt + abs(f_sol(i))**2
           norm2_error_pt = norm2_error_pt + abs(f_sol(i) - f_tn1(i))**2

           if ( abs(f_sol(i) - f_tn1(i)) >  norm_infinite ) norm_infinite = abs(f_sol(i) - f_tn1(i))
           if ( f_tn1(i) < f_min ) f_min = f_tn1(i)

        end do
        
        
         if ( num_method == 15 ) then

           norm2_error = sqrt( ( sqrt(3._f64)*0.5_f64*norm2_error_pt + norm2_error_edge*sqrt(3._f64)*0.125_f64 )*radius**2/real(num_cells,f64)**2)
           norm2_sol = sqrt( ( sqrt(3._f64)*0.5_f64*norm2_sol_pt + norm2_sol_edge*sqrt(3._f64)*0.125_f64 )*radius**2/real(num_cells,f64)**2)


         else if ( num_method == 10 ) then

           norm2_error = sqrt( ( sqrt(3._f64)* 0.5_f64 * norm2_error_pt + &
                norm2_error_center*sqrt(3._f64)/6._f64 )*radius**2/&
                real(num_cells,f64)**2) 

           norm2_sol = sqrt( ( sqrt(3._f64)* 0.5_f64 * norm2_sol_pt + &
                norm2_sol_center*sqrt(3._f64)/6._f64 )*radius**2/&
                real(num_cells,f64)**2) 


        else
           norm2_error = sqrt(sqrt(3._f64)*0.5_f64*norm2_error_pt*radius**2/ &
                real(num_cells,f64)**2)
           norm2_sol =  sqrt(sqrt(3._f64)*0.5_f64*norm2_sol_pt*radius**2/ &
                real(num_cells,f64)**2)
        endif

        center_values_tn = center_values_tn1

        edge_values_tn  = edge_values_tn1

        
        ! call int2string(nloops,filenum)
        ! filename2 = "ana_dist"//trim(filenum)
        ! filename  = "num_dist"//trim(filenum)
        ! call write_field_hex_mesh_xmf(mesh, f_tn1, trim(filename))
        ! call write_field_hex_mesh_xmf(mesh, f_sol, trim(filename2))

        f_tn = f_tn1

     end do
     
     !call write_field_hex_mesh(mesh, f_tn, "result_hex.txt")

     SLL_DEALLOCATE_ARRAY(f_init,ierr)
     SLL_DEALLOCATE_ARRAY(f_tn,ierr)
     SLL_DEALLOCATE_ARRAY(f_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(f_sol,ierr)
     SLL_DEALLOCATE_ARRAY(phi,ierr)
     SLL_DEALLOCATE_ARRAY(uxn,ierr)
     SLL_DEALLOCATE_ARRAY(uyn,ierr)
     SLL_DEALLOCATE_ARRAY(dxuxn,ierr)
     SLL_DEALLOCATE_ARRAY(dyuxn,ierr)
     SLL_DEALLOCATE_ARRAY(dxuyn,ierr)
     SLL_DEALLOCATE_ARRAY(dyuyn,ierr)
     SLL_DEALLOCATE_ARRAY(center_values_tn,ierr)
     SLL_DEALLOCATE_ARRAY(center_values_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(edge_values_tn,ierr)
     SLL_DEALLOCATE_ARRAY(edge_values_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(x1,ierr)
     SLL_DEALLOCATE_ARRAY(x2,ierr)
     SLL_DEALLOCATE_ARRAY(x1_char,ierr)
     SLL_DEALLOCATE_ARRAY(x2_char,ierr)

     deallocate(deriv)
 
     call delete_hex_mesh_2d( mesh )

     call cpu_time(t_end)

     print*, "time used =", t_end - t_init," error_L2 = ", norm2_error


        if ( num_method == 15 ) then
           
           write(33,*) num_cells, n_points + n_edge,  dt, cfl,  norm2_error,  norm2_sol, norm_infinite, f_min, t_end - t_init,&
                tmax/dt * (3._f64*real(num_cells+1,f64)*real(num_cells,f64) + &
                9._f64*real(num_cells,f64)**2 + 3._f64*real(num_cells,f64))/(t_end - t_init)/ 1e6_f64

        elseif ( num_method == 10 ) then

           write(33,*) num_cells, n_points+n_triangle,  dt, cfl,  norm2_error,  norm2_sol, norm_infinite, f_min, t_end - t_init,&
                tmax/dt * (3._f64*real(num_cells+1,f64)*real(num_cells,f64) + &
                6._f64*real(num_cells,f64)**2)/(t_end - t_init)/ 1e6_f64

        else

           write(33,*) num_cells, n_points, dt, cfl,  norm2_error,  norm2_sol, norm_infinite, f_min, t_end - t_init,&
                tmax/dt * 3._f64*real(num_cells + 1, f64)*real(num_cells, f64)/(t_end - t_init)/ 1e6_f64

        endif


  end do

  close(33)

  
contains


  subroutine compute_hex_fields(mesh,uxn,uyn,dxux,dyux,dxuy,dyuy,phi,type)
    type(sll_hex_mesh_2d), pointer :: mesh
    sll_real64,dimension(:)        :: uxn, uyn, phi,dxux,dyux,dxuy,dyuy
    sll_int32,          intent(in) :: type
    sll_int32  :: i!,h1,h2
    !sll_real64 :: phii_2, phii_1, phii1, phii2, phij_2, phij_1, phij1, phij2
    !sll_real64 ::  uh2

    
    if (type==1) then

       do i = 1,mesh%num_pts_tot

          uxn(i) = + mesh%cartesian_coord(2,i)   ! +y
          uyn(i) = - mesh%cartesian_coord(1,i)   ! -x

          dxux(i) = + 0._f64 
          dyux(i) = + 1._f64 
          dxuy(i) = - 1._f64 
          dyuy(i) = - 0._f64 

       end do
    endif
  end subroutine compute_hex_fields


end program test_hex_hermite_circ














