program test_hex_hermite

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_constants
  use interpolation_hex_hermite

  implicit none

  sll_real64, dimension(:,:), allocatable :: deriv, dbc, mesh_num
  sll_int32 , dimension(:,:), allocatable :: mesh_numh, mesh_coor

  sll_int32    :: num_cells, n_points
  sll_int32    :: i, j
  sll_int32    :: nloops
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
  ! distribution at time n + 1
  sll_real64,dimension(:),allocatable :: f_tn1
  ! distribution at end time
  sll_real64,dimension(:),allocatable :: f_sol
  sll_real64   :: diff_error
  sll_real64   :: norm2_error
  ! advection
  sll_int32    :: which_advec
  sll_real64   :: advec
  sll_real64   :: dt
  sll_real64   :: tmax
  sll_real64   :: t

  sll_real64   :: t_init, t_end!, t1,t2,t3,t4
  character(len = 4) :: number
  sll_int32    :: p = 6 !-> degree of the approximation for the derivative 
  sll_real64   :: step , aire, radius, x0, y0, h1, h2, f_min, x ,y, z! , x1_temp
  sll_real64   :: erl11, erl12, erl13
  logical      :: inside
  x0 = 0._f64
  y0 = 0._f64
  radius = 6._f64


  do num_cells = 80,80,10 ! -> loop on the size of the mesh 
  


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !             allocation
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     n_points = 1 + 3 * num_cells * (num_cells + 1)  

     step = radius / real( num_cells )
     aire = step**2*sqrt(3.0_f64)*0.25_f64

     allocate( mesh_num(1:2,1:n_points) )
     allocate( mesh_numh(1:2,1:n_points) )
     allocate( mesh_coor(-num_cells:num_cells,-num_cells:num_cells) )
     allocate( dbc(1:6,num_cells) )
     allocate( deriv(1:6,n_points) )

     dbc  = 0._f64

     SLL_ALLOCATE(f_init( n_points),ierr)
     SLL_ALLOCATE(f_tn( n_points),ierr)
     SLL_ALLOCATE(f_tn1( n_points),ierr)
     SLL_ALLOCATE(f_sol( n_points),ierr)
     SLL_ALLOCATE(x1( n_points),ierr)
     SLL_ALLOCATE(x2( n_points),ierr)
     SLL_ALLOCATE(x1_char( n_points),ierr)
     SLL_ALLOCATE(x2_char( n_points),ierr)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Mesh initialization   
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     

     call create_hex_mesh(mesh_num, mesh_numh, mesh_coor, num_cells, radius, x0, y0)


     ! Distribution initialization

     gauss_x1  = 1._f64
     gauss_x2  = 1._f64
     gauss_sig = 0.05_f64

     open(unit = 11, file="hex_hermite_init.txt", action="write", status="replace")

     do i = 1, n_points

        x1(i) = mesh_num(1,i)
        x2(i) = mesh_num(2,i)

        ! f_init(i) = exp(-0.5_f64*((x1(i)-gauss_x1)**2 + (x2(i)-gauss_x2)**2) / gauss_sig**2 )

        ! if (exponent(f_init(i)) .lt. -17) then
        !    f_init(i) = 0._f64
        ! end if

        f_init(i) = exp( -(x1(i)-gauss_x1)**2 - (x2(i)-gauss_x2)**2 )

        f_tn(i) = f_init(i)

        write(11,*) x1(i),x2(i),f_tn(i)

     end do

     close(11)

     ! Advection initialization
     which_advec = 1  ! 0 : linear advection ; 1 : circular advection
     advec = 0.025_f64!5_f64
     tmax  = 6.3_f64
     dt    = 0.1_f64
     t     = 0._f64

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !              Computing characteristics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     if (which_advec .eq. 0) then
        ! linear advection
        x1_char(:) = x1(:) - advec*dt
        x2_char(:) = x2(:) - advec*dt
     else
        ! Circular advection
        x1_char(:) = x1(:)*cos(dt) - x2(:)*sin(dt)
        x2_char(:) = x1(:)*sin(dt) + x2(:)*cos(dt)
     end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !                          Time loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     nloops = 0

     call cpu_time(t_init)

     do while (t .lt. tmax)
        !Error variables
        norm2_error = 0._f64
        diff_error  = 0._f64

        nloops = nloops + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! let us compute the derivatives in every hexagonal direction 
        ! with p the degree of the approximation


        !call cpu_time(t1)
        call  der_finite_difference( f_tn, p,step, mesh_numh, mesh_coor, deriv, num_cells, dbc  )
        !call cpu_time(t2)
        !print*, "temps der : ", t2-t1

        erl11 = 0._f64
        erl12 = 0._f64
        erl13 = 0._f64
        

        do i = 1, n_points
           x = mesh_num(1,i)
           y = mesh_num(2,i)  
           erl11 = erl11 + abs(deriv(1,i)+sqrt(3.0)*0.5*2.0*(x-gauss_x1)*exp(-((x-gauss_x1)**2+& 
           (y-gauss_x2)**2 ) ) +  0.5*2.0*(y-gauss_x2) * exp( -( (x-gauss_x1)**2 + &
           (y-gauss_x2)**2 ) ) )
           erl12 = erl12 + abs(deriv(2,i)-sqrt(3.0)*0.5*2.0*(x-gauss_x1)*exp(-((x-gauss_x1)**2 +&
           (y-gauss_x2)**2 ) ) +  0.5*2.0*(y-gauss_x2) * exp( -( (x-gauss_x1)**2 + &
           (y-gauss_x2)**2 ) ) )
           erl13 = erl13 + abs( deriv(3,i) + 2.0*(y-gauss_x2)*exp( -( (x-gauss_x1)**2 + (y-gauss_x2)**2 )  ) )
        enddo

        erl11 = erl11/real(n_points,f64)
        erl12 = erl12/real(n_points,f64)
        erl13 = erl13/real(n_points,f64)

        print*,"norme L1 des dérivées partielles dans les trois directions hex: "
        print*, erl11, erl12, erl13

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        t = t + dt

        f_min = 0._f64
        
        !write(number,"(4I4)")   nloops
        
        !open(unit = 11, file="hex_hermite"//number//".txt", action="write", status="replace")
        
        !call cpu_time(t3)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !              VALIDATION INTERPOLATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        norm2_error = 0
        
        do i = 1,num_cells
           do j = 1,num_cells
              x = -sqrt(3.0)*0.5*step + real(i-1,f64)*sqrt(3.0)/(num_cells + 1)*step
              y = -sqrt(3.0)*0.5*step + real(j-1,f64)*sqrt(3.0)/(num_cells + 1)*step
              z = 0
              
              if ( sqrt(x**2+y**2) >= radius*sqrt(3.0)*0.5 ) then
                 z = -1
              endif
              
              
              if ( z>-1) then
                 call hermite_interpolation(i, x1_char(i), x2_char(i), f_tn, &
                      f_tn1, mesh_num, mesh_coor, deriv, step, aire, num_cells ,radius)  
              endif
              
              norm2_error = norm2_error + abs(f_sol(i) - f_tn1(i))**2
              f_sol(i) = exp(-(x-gauss_x1)**2-(y-gauss_x2)**2) 
           enddo
        enddo

        print*, "error regular interpolation : ",   norm2_error*step**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                 INTERPOLATION for the characteristic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        do i=1, n_points   ! interpolation at each point

           ! computation of the interpolation F(t_(n+1),x_i,v_j) = F (t_n,X(t_n),V(t_n))

           inside = .true.

           h1 =  x1_char(i)/sqrt(3.0_f64) + x2_char(i)
           h2 = -x1_char(i)/sqrt(3.0_f64) + x2_char(i) 

           if ( h1 >  radius .or. h2 >  radius ) inside = .false.
           if ( h1 < -radius .or. h2 < -radius ) inside = .false.
           if ( x1_char(i)  < -radius*sqrt(3._f64)*0.5_f64 .or. x1_char(i) &
                > radius*sqrt(3._f64)*0.5_f64  ) inside = .false.

           if ( inside ) then
              call hermite_interpolation(i, x1_char(i), x2_char(i), f_tn, &
                   f_tn1, mesh_num, mesh_coor, deriv, step, aire, num_cells ,radius)  
           else 
              
              f_tn1(i) = 0._f64

           endif

           ! ******************
           ! Analytical value    (-> in order to compute the error )
           ! ******************

           ! Computing characteristics

           !if (which_advec .eq. 0) then
           ! linear advection
           !  x1(i) = mesh_num( 1, i) - advec*dt*nloops
           !  x2(i) = mesh_num( 2, i) - advec*dt*nloops
           !else
           ! Circular advection


           x = x1(i)*cos(t) - x2(i)*sin(t);
           y = x1(i)*sin(t) + x2(i)*cos(t);

           !x1_temp = sqrt(x1(i)**2 + x2(i)**2) * cos(2*sll_pi*dt + atan2(x2(i),x1(i)))
           !x2(i)   = sqrt(x1(i)**2 + x2(i)**2) * sin(2*sll_pi*dt + atan2(x2(i),x1(i)))
           !x1(i)   = x1_temp
           ! end if

           !f_sol(i) = exp(-0.5_f64*((x1(i)-gauss_x1)**2/gauss_sig**2 &
           !    + (x2(i)-gauss_x2)**2 / gauss_sig**2))

           f_sol(i) = exp(-(x-gauss_x1)**2-(y-gauss_x2)**2) 

           !if (diff_error .lt. abs(f_sol(i) - f_tn1(i)) ) then
           !   diff_error = abs(f_sol(i) - f_tn1(i))
           !end if



           norm2_error = norm2_error + abs(f_sol(i) - f_tn1(i))**2

           !write(11,*) mesh_num(1,i),mesh_num(2,i) ,f_sol(i),f_tn1(i)

           !if ( f_tn1(i) < f_min ) f_min = f_tn1(i)

        end do

        !close(11)


        !call cpu_time(t4)

        !print*, "temps interpolation + diag", t4-t3

        ! Norm2 error :

        norm2_error = sqrt(norm2_error)*radius**2/real(num_cells,f64)**2

        print*,"error_L2 = ", norm2_error!, "min =",f_min

        !write(12,*) t,  norm2_error

        f_tn = f_tn1

     end do


     call cpu_time(t_end)

     SLL_DEALLOCATE_ARRAY(f_init,ierr)
     SLL_DEALLOCATE_ARRAY(f_tn,ierr)
     SLL_DEALLOCATE_ARRAY(f_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(f_sol,ierr)
     SLL_DEALLOCATE_ARRAY(x1,ierr)
     SLL_DEALLOCATE_ARRAY(x2,ierr)
     SLL_DEALLOCATE_ARRAY(x1_char,ierr)
     SLL_DEALLOCATE_ARRAY(x2_char,ierr)


     deallocate( mesh_num, mesh_numh, mesh_coor, deriv) 

     print*, "time used =", t_end - t_init

  end do


end program test_hex_hermite














