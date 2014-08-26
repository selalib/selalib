program test_hex_hermite

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_constants
  use interpolation_hex_hermite

  implicit none

  sll_real64, dimension(:,:), allocatable :: deriv

  sll_int32    :: num_cells, n_points, n_triangle
  sll_int32    :: i
  sll_int32    :: nloops, num_method = 10
  sll_int32    :: ierr
  ! initial distribution
  sll_real64   :: gauss_x2
  sll_real64   :: gauss_x1
  sll_real64   :: gauss_sig
  sll_real64,dimension(:),allocatable :: x1
  sll_real64,dimension(:),allocatable :: x2
  sll_real64,dimension(:),allocatable :: f_init
  ! distribution at time n
  sll_real64,dimension(:),allocatable :: x1_char
  sll_real64,dimension(:),allocatable :: x2_char
  sll_real64,dimension(:),allocatable :: f_tn
  sll_real64,dimension(:),allocatable :: center_values_tn
  ! distribution at time n + 1
  sll_real64,dimension(:),allocatable :: f_tn1,center_values_tn1
  ! distribution at end time
  sll_real64,dimension(:),allocatable :: f_sol
  sll_real64   :: norm2_error
  ! advection
  sll_int32    :: which_advec
  sll_real64   :: advec
  sll_real64   :: dt
  sll_real64   :: tmax
  sll_real64   :: t

  sll_real64   :: t_init, t_end
  sll_int32    :: p = 6 !-> degree of the approximation for the derivative 
  sll_real64   :: step , aire, radius, center_x1, center_x2, h1, h2, f_min, x ,y
  sll_real64   :: xx, yy!,z, x1_temp,t1,t2,t3,t4 
  sll_real64   :: cfl, norm_infinite
  ! character(len = 4) :: number
  logical      :: inside
  type(hex_mesh_2d), pointer :: mesh

  
  center_x1 = 0._f64
  center_x2 = 0._f64

  radius = 8._f64

  if (num_method == 9 ) then
     print*, 
     print*, "*********************************"
     print*, " Zienkiewicz_9_degree_of_freedom "
     print*, "*********************************"
     print*, 
  else if (num_method == 10 ) then
     print*, 
     print*, "*********************************"
     print*," Zienkiewicz_10_degree_of_freedom"
     print*, "*********************************"
     print*, 
  else if (num_method == 11 ) then 
     print*, 
     print*, "*********************************"
     print*, "   Hsieh_Clough_Tocher_reduced   "
     print*, "*********************************"
     print*, 
  else if (num_method == 12 ) then 
     print*, 
     print*, "*********************************"
     print*, "   Hsieh_Clough_Tocher_complete   "
     print*, "*********************************"
     print*, 
  else if (num_method == 13 ) then 
     print*, 
     print*, "*********************************"
     print*, "  quartic element of Ganev_Dimitrov "
     print*, "*********************************"
     print*, 
  else
     print*, "specify another number correspoonding to a existing implemented method 9, 10, 11, 12 or 13"
  endif

  open(unit = 33, file="hex_errors.txt", action="write", status="replace")!,position = "append")! 

  write(33,*) 

  do num_cells = 20,340,40 ! -> loop on the size of the mesh 
  
     
     !*********************************************************
     !             allocation
     !*********************************************************
     
     n_points = 1 + 3 * num_cells * (num_cells + 1)  
     n_triangle = 6 * num_cells * num_cells

     step = radius / real(num_cells,f64)
     aire = step**2*sqrt(3._f64)*0.25_f64

     allocate( deriv(1:6,n_points) )

     SLL_ALLOCATE(f_init( n_points),ierr)
     SLL_ALLOCATE(f_tn( n_points),ierr)
     SLL_ALLOCATE(f_tn1( n_points ),ierr)
     SLL_ALLOCATE(center_values_tn( n_triangle),ierr)
     SLL_ALLOCATE(center_values_tn1( n_triangle),ierr)
     SLL_ALLOCATE(f_sol( n_points ),ierr)
     SLL_ALLOCATE(x1( n_points),ierr)
     SLL_ALLOCATE(x2( n_points),ierr)
     SLL_ALLOCATE(x1_char( n_points),ierr)
     SLL_ALLOCATE(x2_char( n_points),ierr)

     !*********************************************************
     !                  Mesh initialization   
     !*********************************************************
     
     mesh => new_hex_mesh_2d( num_cells, center_x1, center_x2, radius=radius ) 

     ! Distribution initialization

     gauss_x1  = 1._f64
     gauss_x2  = 1._f64
     gauss_sig = 1._f64/ sqrt(2._f64)!( 2._f64 * sqrt(2._f64)) 

     open(unit = 11, file="hex_hermite_init.txt", action="write", status="replace")

     do i = 1, n_points
        x1(i) = mesh%cartesian_coord(1,i)
        x2(i) = mesh%cartesian_coord(2,i)

        f_init(i) = exp( -((x1(i)-gauss_x1)**2 + (x2(i)-gauss_x2)**2)/gauss_sig**2/2._f64 )
        f_tn(i) = f_init(i)

        write(11,*) x1(i),x2(i),f_tn(i)

     end do


     if ( num_method == 10 ) then
        do i = 1, n_triangle
           x = mesh%center_cartesian_coord(1,i)
           y = mesh%center_cartesian_coord(2,i)
           center_values_tn(i) = exp( -((x-gauss_x1)**2 + (y-gauss_x2)**2)/gauss_sig**2/2._f64 )
           write(11,*) x,y,center_values_tn(i)
        enddo
     endif

     close(11)

     ! Advection initialization
     which_advec = 1  ! 0 : linear advection ; 1 : circular advection
     advec = 0.025_f64!5._f64
     tmax  = 1._f64
     dt    = 0.1_f64*20._f64 / real(num_cells,f64)  
     t     = 0._f64
     cfl   = radius * dt / ( radius / real(num_cells,f64)  )

     !*********************************************************
     !              Computing characteristics
     !*********************************************************


     if (which_advec .eq. 0) then
        ! linear advection
        x1_char(:) = x1(:) - advec*dt
        x2_char(:) = x2(:) - advec*dt
     else
        ! Circular advection
        x1_char(:) = x1(:)*cos(dt) - x2(:)*sin(dt)
        x2_char(:) = x1(:)*sin(dt) + x2(:)*cos(dt)
     end if


     !*********************************************************
     !                          Time loop
     !*********************************************************

     nloops = 0

     call cpu_time(t_init)

     do while (t .lt. tmax)!dt)!

        norm2_error = 0._f64 !Error variables
        norm_infinite = 0._f64
        f_min = 0._f64

        nloops = nloops + 1

        !*********************************************************
        !              VALIDATION DERIVATIVE
        ! let us compute the derivatives in every hexagonal direction 
        ! with p the degree of the approximation
        !*********************************************************

        call  der_finite_difference( f_tn, p, step, mesh, deriv )

        t = t + dt
        !*********************************************************
        ! computation of the value at the center of the triangles
        !*********************************************************

        if ( num_method == 10 ) then 
           
           do i=1, n_triangle  ! computation of the value at the center of the triangles

              !*********************************************************
              !  computation of the root of the characteristics
              !*********************************************************

              x = mesh%center_cartesian_coord(1,i)
              y = mesh%center_cartesian_coord(2,i)

              xx = x*cos(dt) - y*sin(dt);
              yy = x*sin(dt) + y*cos(dt);
              ! call slb_compute_characteristic_leapfrog( &
              !  x,y,E_x,E_v,xx,yy )

              !             INTERPOLATION
              inside = .true.

              h1 =  xx/sqrt(3.0_f64) + yy
              h2 = -xx/sqrt(3.0_f64) + yy 

              if ( h1 >  radius .or. h2 >  radius ) inside = .false.
              if ( h1 < -radius .or. h2 < -radius ) inside = .false.
              if ( xx  < -radius*sqrt(3._f64)*0.5_f64 .or. xx &
                   > radius*sqrt(3._f64)*0.5_f64  ) inside = .false.

              if ( inside ) then
                 call hermite_interpolation(i, xx, yy, f_tn, center_values_tn, &
                      center_values_tn1, mesh, deriv, aire,t-dt, num_method)
              else 
                 center_values_tn1(i) = 0._f64 ! dirichlet boundary condition
              endif

              if (which_advec .eq. 0) then ! linear advection
                 xx = x - advec*dt*nloops
                 yy = y - advec*dt*nloops
              else                         ! Circular advection
                 xx = x*cos(t) - y*sin(t);
                 yy = x*sin(t) + y*cos(t);
              end if

              
             ! if (center_values_tn1(i)>1.) print*, i, center_values_tn(i), center_values_tn1(i)

              !norm2_error = norm2_error + &
              !     abs( exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64) - center_values_tn1(i) )**2
              !if ( abs( exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64) - center_values_tn1(i)) >  norm_infinite )&
              !     norm_infinite = abs( exp(-((xx-gauss_x1)**2+(yy-gauss_x2)**2)/gauss_sig**2/2._f64) - center_values_tn1(i))
              !if ( center_values_tn1(i) < f_min ) f_min = center_values_tn1(i)

           enddo

        endif

        do i=1, n_points  

        !*********************************************************
        !  computation of the root of the characteristics
        !*********************************************************

        ! call slb_compute_characteristic_leapfrog( &
        !   x1(i),x2(i),E_x,E_v,x1_char(i),x2_char(i) )

        !*********************************************************
        !                INTERPOLATION
        !*********************************************************
           ! computation of the interpolation at each point
           ! F(t_(n+1),x_i,v_j) = F (t_n,X(t_n),V(t_n))

           inside = .true.

           h1 =  x1_char(i)/sqrt(3.0_f64) + x2_char(i)
           h2 = -x1_char(i)/sqrt(3.0_f64) + x2_char(i) 

           if ( h1 >  radius .or. h2 >  radius ) inside = .false.
           if ( h1 < -radius .or. h2 < -radius ) inside = .false.
           if ( x1_char(i)  < -radius*sqrt(3._f64)*0.5_f64 .or. x1_char(i) &
                > radius*sqrt(3._f64)*0.5_f64  ) inside = .false.

           if ( inside ) then
              call hermite_interpolation(i, x1_char(i), x2_char(i), f_tn, &
                   center_values_tn, f_tn1, mesh, deriv, aire,t-dt, num_method)
           else 
              
              f_tn1(i) = 0._f64 ! dirichlet boundary condition

           endif

           ! ******************************************************
           ! Analytical value    (-> in order to compute the error)
           ! ******************************************************

           ! Computing the characteristics' root for the exact solution
           ! in order to compute the l2 norm of the error
           
           if (which_advec .eq. 0) then ! linear advection
              x1(i) = mesh%cartesian_coord(1,i) - advec*dt*nloops
              x2(i) = mesh%cartesian_coord(2,i) - advec*dt*nloops
           else                         ! Circular advection
              x = x1(i)*cos(t) - x2(i)*sin(t);
              y = x1(i)*sin(t) + x2(i)*cos(t);
           end if

           f_sol(i) = exp(-((x-gauss_x1)**2+(y-gauss_x2)**2)/gauss_sig**2/2._f64) 
           norm2_error = norm2_error + abs(f_sol(i) - f_tn1(i))**2

           if ( abs(f_sol(i) - f_tn1(i)) >  norm_infinite ) norm_infinite = abs(f_sol(i) - f_tn1(i))
           if ( f_tn1(i) < f_min ) f_min = f_tn1(i)

        end do

        if ( num_method == 10 ) then
           norm2_error = sqrt(norm2_error)/(3._f64*real(num_cells+1,f64)*real(num_cells,f64) + 6*real(num_cells,f64)**2)
        else
           norm2_error = sqrt(norm2_error)/(3._f64*real(num_cells,f64)*real(num_cells+1,f64))
        endif
        !print*,"error_L2 = ", norm2_error!, "min =",f_min

        f_tn = f_tn1
        center_values_tn = center_values_tn1

     end do


     call cpu_time(t_end)

     SLL_DEALLOCATE_ARRAY(f_init,ierr)
     SLL_DEALLOCATE_ARRAY(f_tn,ierr)
     SLL_DEALLOCATE_ARRAY(f_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(f_sol,ierr)
     SLL_DEALLOCATE_ARRAY(center_values_tn,ierr)
     SLL_DEALLOCATE_ARRAY(center_values_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(x1,ierr)
     SLL_DEALLOCATE_ARRAY(x2,ierr)
     SLL_DEALLOCATE_ARRAY(x1_char,ierr)
     SLL_DEALLOCATE_ARRAY(x2_char,ierr)

     deallocate(deriv)
 
     call delete_hex_mesh_2d( mesh )

     print*, "time used =", t_end - t_init," error_L2 = ", norm2_error


        if ( num_method == 10 ) then
           
           write(33,*) num_cells, n_points+n_triangle,  dt, cfl,  norm2_error, norm_infinite, f_min, t_end - t_init,&
                tmax/dt * (3._f64*real(num_cells+1,f64)*real(num_cells,f64) + &
                6*real(num_cells,f64)**2)/(t_end - t_init)/ 1e6_f64
        else
           write(33,*) num_cells, n_points, dt, cfl,  norm2_error, norm_infinite, f_min, t_end - t_init,&
                tmax/dt * 3._f64*real(num_cells + 1, f64)*real(num_cells, f64)/(t_end - t_init)/ 1e6_f64
        endif

  end do

  close(33)

end program test_hex_hermite














