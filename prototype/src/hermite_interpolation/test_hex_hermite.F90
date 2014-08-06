program test_hex_hermite

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use sll_constants
  use interpolation_hex_hermite

  implicit none

  sll_real64, dimension(:,:), allocatable :: mesh_num, deriv, dbc
  sll_int32 , dimension(:,:), allocatable :: mesh_numh, mesh_coor

  sll_int32    :: num_cells, n_points
  sll_int32    :: i
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
  sll_real64,dimension(:),allocatable :: f_fin
  sll_real64   :: diff_error
  sll_real64   :: norm2_error
  ! advection
  sll_int32    :: which_advec
  sll_real64   :: advec
  sll_real64   :: dt
  sll_real64   :: tmax
  sll_real64   :: t

  sll_real64   :: t_init, t_end
  !character(len = 50) :: filename
  sll_int32    :: p = 6 !-> degree of the approximation for the derivative 
  sll_real64   :: step , aire, radius, x0, y0, x1_temp
  
  x0 = 0._f64
  y0 = 0._f64
  radius = 1._f64



  do num_cells = 3,3,1 ! -> loop on the size of the mesh 
  
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Mesh initialization   
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     n_points = 1 + 3 * num_cells * (num_cells + 1)  

     allocate( mesh_num(1:2,1:n_points) )
     allocate( mesh_numh(1:2,1:n_points) )
     allocate( mesh_coor(-num_cells:num_cells,-num_cells:num_cells) )
     allocate( dbc(1:6,num_cells) )
     allocate( deriv(1:6,n_points) )

     call create_hex_mesh(mesh_num, mesh_numh, mesh_coor, num_cells, radius, x0, y0)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !             allocation
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     SLL_ALLOCATE(f_init( n_points),ierr)
     SLL_ALLOCATE(f_tn( n_points),ierr)
     SLL_ALLOCATE(f_tn1( n_points),ierr)
     SLL_ALLOCATE(f_fin( n_points),ierr)
     SLL_ALLOCATE(x1( n_points),ierr)
     SLL_ALLOCATE(x2( n_points),ierr)
     SLL_ALLOCATE(x1_char( n_points),ierr)
     SLL_ALLOCATE(x2_char( n_points),ierr)


     step = radius / real( num_cells )
     aire = step**2*sqrt(3.0_f64)*0.25_f64

     ! Distribution initialization

     gauss_x1  = -0.25_f64
     gauss_x2  = -0.25_f64
     gauss_sig = 0.05_f64

     do i = 1, n_points

        !(/x1(i), x2(i)/) = (/mesh_num(1,i),mesh_num(2,i)/)!mesh_num( 1:2 , i)

        x1(i) = mesh_num(1,i)
        x2(i) = mesh_num(2,i)

        f_init(i) = exp(-0.5_f64*((x1(i)-gauss_x1)**2 + (x2(i)-gauss_x2)**2) / gauss_sig**2 )
        if (exponent(f_init(i)) .lt. -17) then
           f_init(i) = 0._f64
        end if
        f_tn(i) = f_init(i)

     end do


     ! Advection initialization
     which_advec = 1  ! 0 : linear advection ; 1 : circular advection
     advec = 0.025_f64!5_f64
     tmax  = 2.5_f64
     dt    = 0.025_f64
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
        x1_char(2:) = sqrt(x1(2:)**2 + x2(2:)**2) * cos(2*sll_pi*dt + atan2(x2(2:),x1(2:)))
        x2_char(2:) = sqrt(x1(2:)**2 + x2(2:)**2) * sin(2*sll_pi*dt + atan2(x2(2:),x1(2:)))
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! let us compute the derivatives in every hexaedric direction with p the degree of the approximation
        ! then let us compute the derivatives in the x1 and x2 direction

        call  der_finite_difference( f_tn, p,step, mesh_numh, mesh_coor, deriv, num_cells, dbc  )  !-> CALCUL DE TOUTES LES DÉRIVÉES DANS LES DIRECTIONS H1, H2 et H3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        t = t + dt

        do i=1, n_points   ! interpolation en chaque point

           ! ******************
           ! Approximation
           ! ******************

           ! computation of the interpolation   F( t_(n+1) , x_i , v_j ) = F ( t_n , X(t_n) , V(t_n) ) 

           call hermite_interpolation(i, x1(i), x2(i), f_tn, f_tn1, mesh_num, mesh_coor, deriv, step, aire, num_cells ,radius)  

           ! if (f_tn(i) < 0.) then
           !    ! print *, "Negative value at"
           !    ! print *, "      Time  :", t
           !    ! print *, "      Loop  :", nloops
           !    ! print *, "      Point :", i
           !    ! print *, "      X1char:", x1_char(i)
           !    ! print *, "      X2char:", x2_char(i)
           !    ! STOP
           !    f_tn(i) = 0._f64
           ! end if

           ! ******************
           ! Analytical value    (-> in order to compute the error )
           ! ******************

           ! Computing characteristics

           if (which_advec .eq. 0) then
              ! linear advection
              x1(i) = mesh_num( 1, i) - advec*dt*nloops
              x2(i) = mesh_num( 2, i) - advec*dt*nloops
           else
              ! Circular advection
              x1_temp = sqrt(x1(i)**2 + x2(i)**2) * cos(2*sll_pi*dt + atan2(x2(i),x1(i)))
              x2(i)   = sqrt(x1(i)**2 + x2(i)**2) * sin(2*sll_pi*dt + atan2(x2(i),x1(i)))
              x1(i)   = x1_temp
           end if

           f_fin(i) = exp(-0.5_f64*((x1(i)-gauss_x1)**2/gauss_sig**2 &
                + (x2(i)-gauss_x2)**2 / gauss_sig**2))

           if (exponent(f_fin(i)) .lt. -17) then
              f_fin(i) = 0._f64
           end if

              ! Relative error
              if (diff_error .lt. abs(f_fin(i) - f_tn(i)) ) then
                 diff_error = abs(f_fin(i) - f_tn(i))
              end if
              ! Norm2 error :
              norm2_error = norm2_error + abs(f_fin(i) - f_tn(i))**2

           end do

        ! Norm2 error :
        norm2_error = sqrt(norm2_error)

        call cpu_time(t_end)
        print*,"  nt =", nloops, "     | error_Linf = ", diff_error
        print*,"                       | error_L2   = ", norm2_error

        f_tn = f_tn1

     end do

     print*," *    Final error  = ", diff_error, " *"


     SLL_DEALLOCATE_ARRAY(f_init,ierr)
     SLL_DEALLOCATE_ARRAY(f_tn,ierr)
     SLL_DEALLOCATE_ARRAY(f_tn1,ierr)
     SLL_DEALLOCATE_ARRAY(f_fin,ierr)
     SLL_DEALLOCATE_ARRAY(x1,ierr)
     SLL_DEALLOCATE_ARRAY(x2,ierr)
     SLL_DEALLOCATE_ARRAY(x1_char,ierr)
     SLL_DEALLOCATE_ARRAY(x2_char,ierr)
     deallocate( mesh_num, mesh_numh, mesh_coor, deriv) 


  end do


  ! contains

    
  !   subroutine der_finite_difference( f_tn, p, mesh_numh, mesh_coor, deriv, n_cell ) 
  !     !-> computation of the partial derivatives in the directions H1, H2 and H3
  !     implicit none
  !     sll_int32 , dimension(:,:)           :: mesh_coor, mesh_numh
  !     sll_real64, dimension(:,:)           :: deriv 
  !     sll_real64, dimension(:), intent(in) :: f_tn 
  !     sll_int32, intent(in)                :: n_cell, p
  !     sll_int32                            :: i, j, r, s, n_points
  !     sll_int32                            :: n1, n2, n3, n4, n5, n6
  !     sll_real64                           :: dh1, dh2, dh3, dh4, dh5, dh6
  !     sll_real64, dimension(:),allocatable :: w
  !     sll_real64, dimension(:,:),allocatable :: f
  !     sll_int32                            :: h1 ,h2, h1t ,h2t
  !     sll_int32                            :: h1n ,h2n, h1tn ,h2tn
  !     logical                              :: opsign

      
  !     n_points = size(f_tn)

  !     !*********************************************************************
  !     !computation of the coefficients in fonction of the degree p 
  !     !of the approximation required
  !     !*********************************************************************
  !     if ( (p/2)*2 == p ) then !if p is even
  !        r = -p/2       
  !        s = p/2
  !     else  !if p is odd
  !        r = -p/2 - 1 ! de-centered to the left      
  !        s = p/2 !+ 1 ! de-centered to the right
  !     endif

  !     SLL_ALLOCATE( w(r:s),ierr )
  !     allocate( f(1:6,r:s) )

  !     call compute_w_hermite(w,r,s) 

  !     !***********************************************************************
  !     ! Computation of the derivatives in both hexa directions h1 and h2,
  !     ! then for both x and y directions for every points of the mesh
  !     !***********************************************************************
      
  !     do i = 1, n_points 

  !        h1 = mesh_numh(1,i) 
  !        h2 = mesh_numh(2,i)  
  !        h1n = h1 + n_cell + 1
  !        h2n = h2 + n_cell + 1

  !        do j = r,s !find the required points in the h1, h2 and h3 directions 
            
  !           opsign = .false.

  !           h1t = h1 + j 
  !           h2t = h2 + j
  !           h1tn = h1n + j 
  !           h2tn = h2n + j

  !           if ( h2t > 0 .and. h1t < 0 .or.  h1t > 0 .and. h2t < 0 ) opsign = .true.

  !           !test if outside mesh 

  !           if (  abs(h1t) > n_cell .or. opsign .and. ( abs(h1t) + abs(h2) > n_cell) ) then
  !              f(1,j) = 0.0_f64 ! dirichlet boundary condition
  !           else
  !              n1 = mesh_coor(h1tn,h2n)
  !              f(1,j) = f_tn(n1) 
  !           endif

  !           if (  abs(h2t) > n_cell .or. opsign .and. ( abs(h1) + abs(h2t) > n_cell) ) then
  !              f(2,j) = 0.0_f64 ! dirichlet boundary condition
  !           else
  !              n2 = mesh_coor(h1n,h2tn)
  !              f(2,j) = f_tn(n2) 
  !           endif

  !           if ( abs(h1t) > n_cell .or. abs(h2t) > n_cell .or. opsign .and. ( abs(h1t) + abs(h2t) > n_cell) ) then
  !              f(3,j) = 0.0_f64 ! dirichlet boundary condition
  !           else
  !              n3 = mesh_coor(h1tn,h2tn)
  !              f(3,j) = f_tn(n1) 
  !           endif

  !           !find the required points in the -h1, -h2 and -h3 directions 
            
  !           opsign = .false.

  !           h1t = h1 - j
  !           h2t = h2 - j
  !           h1tn = h1n + j 
  !           h2tn = h2n + j

  !           if ( h2t > 0 .and. h1t < 0 .or.  h1t > 0 .and. h2t < 0 ) opsign = .true.

  !           !test if outside mesh 

  !           if (  abs(h1t) > n_cell .or. opsign .and. ( abs(h1t) + abs(h2) > n_cell) ) then
  !              f(4,j) = 0.0_f64 ! dirichlet boundary condition
  !           else
  !              n4 = mesh_coor(h1tn,h2n)
  !              f(4,j) = f_tn(n4) 
  !           endif

  !           if (  abs(h2t) > n_cell .or. opsign .and. ( abs(h1) + abs(h2t) > n_cell) ) then
  !              f(5,j) = 0.0_f64 ! dirichlet boundary condition
  !           else
  !              n5 = mesh_coor(h1n,h2tn)
  !              f(5,j) = f_tn(n5) 
  !           endif

  !           if ( abs(h1t) > n_cell .or. abs(h2t) > n_cell .or. opsign .and. ( abs(h1t) + abs(h2t) > n_cell) ) then
  !              f(6,j) = 0.0_f64 ! dirichlet boundary condition
  !           else
  !              n6 = mesh_coor(h1tn,h2tn)
  !              f(6,j) = f_tn(n6) 
  !           endif

  !        enddo


  !        dh1 = 0._f64
  !        dh2 = 0._f64
  !        dh3 = 0._f64
  !        dh4 = 0._f64
  !        dh5 = 0._f64
  !        dh6 = 0._f64

  !        do j = r,s
  !           dh1 = dh1 + w(j) * f(1,j)  
  !           dh2 = dh2 + w(j) * f(2,j)
  !           dh3 = dh3 + w(j) * f(3,j)
  !           dh4 = dh4 + w(j) * f(4,j)  
  !           dh5 = dh5 + w(j) * f(5,j)
  !           dh6 = dh6 + w(j) * f(6,j)
  !        enddo
         
  !        deriv(1,i) = dh1
  !        deriv(2,i) = dh2
  !        deriv(3,i) = dh3 
  !        deriv(4,i) = dh4
  !        deriv(5,i) = dh5
  !        deriv(6,i) = dh6 

  !     enddo


  !    deallocate(f)
  !    SLL_DEALLOCATE_ARRAY(w,ierr)

  !   end subroutine der_finite_difference


    

  !   subroutine hermite_interpolation(num, x, y, f_tn, f_tn1, mesh_num, mesh_coor, deriv, step, aire, num_cells ,radius) 

  !     implicit none
  !     sll_int32 , dimension(:,:) :: mesh_coor
  !     sll_real64, dimension(:,:) :: mesh_num, deriv 
  !     sll_real64,intent(in)      :: x, y, step, aire, radius
  !     sll_int32,intent(in)       :: num_cells, num
  !     sll_real64                 :: x1, x2, x3, y1, y2, y3, f
  !     sll_real64                 :: phi, ksi1, ksi2, ksi3, l1, l2, l3, eps = 1.d-12
  !     sll_real64                 :: ksi12, ksi13, ksi21, ksi23, ksi31, ksi32
  !     sll_real64,dimension(:), intent(in)    :: f_tn
  !     sll_real64,dimension(:), intent(out)   :: f_tn1
  !     sll_real64,dimension(1:10) :: freedom, base
  !     sll_int32                  :: i, i1, i2, i3

  !     ! find the triangle where (x,y) belongs to
  !     ! i1, i2, i3 are the indices for the S1, S2 and S3 vertexes of the triangle

  !     call search_tri( x, y, mesh_coor,step, num_cells, radius, i1, i2, i3 )

  !     ! get the ten degrees of freedom

  !     ! values at the vertexes of the triangle

  !     freedom(1) = f_tn(i1)
  !     freedom(2) = f_tn(i2)
  !     freedom(3) = f_tn(i3)

  !     x1 = mesh_num(1,i1) !(/x1,y1/) = mesh_num(1:2,i1)
  !     x2 = mesh_num(1,i2)
  !     x3 = mesh_num(1,i3)
  !     y1 = mesh_num(2,i1)
  !     y2 = mesh_num(2,i2)
  !     y3 = mesh_num(2,i3)

  !     ! value at the center of the triangle

  !     freedom(4) =  ( freedom(1) + freedom(2) + freedom(3) ) / 3.0_f64

  !     ! values of the derivatives

  !     freedom(5) = deriv(5,i1) ! derivative from S1 to S3
  !     freedom(9) = deriv(2,i3) ! derivative from S3 to S1
      
  !     !if ( (/x2-x1,y2-y1/) == h1 )then      ! if S2 S1 = h1 then we are in the first case )
  !     if ( ((x2-x1) - sqrt(3.0_f64)*0.5_f64)**2 + (y2-y1 - 0.5_f64)**2 <= eps ) then 
  !        freedom(6) = deriv(4,i1) ! derivative from S1 to S2
  !        freedom(7) = deriv(1,i2) ! derivative from S2 to S1
  !        freedom(8) = deriv(6,i2) ! derivative from S2 to S3
  !        freedom(10)= deriv(3,i3) ! derivative from S3 to S2
  !     else !if ( (/x2-x1,y2-y1/) == -h2 )! if S2 S1 = -h2 then we are in the 2nd case )
  !        freedom(6) = deriv(6,i1) ! derivative from S1 to S2 
  !        freedom(7) = deriv(3,i2) ! derivative from S2 to S1
  !        freedom(8) = deriv(4,i2) ! derivative from S2 to S3
  !        freedom(10)= deriv(1,i3) ! derivative from S3 to S2
  !     endif

  !     l1   = 0.5_f64 * abs( (x2 - x)*(y3 - y) - (x3 - x)*(y2 - y) ) / aire 
  !     l2   = 0.5_f64 * abs( (x1 - x)*(y3 - y) - (x3 - x)*(y1 - y) ) / aire 
  !     l3   = 0.5_f64 * abs( (x1 - x)*(y2 - y) - (x2 - x)*(y1 - y) ) / aire 

  !     phi = l1*l2*l3

  !     ksi1 = l1**3 - phi
  !     ksi2 = l2**3 - phi
  !     ksi3 = l3**3 - phi

  !     ksi12 = l1**2*l2 + phi*0.5_f64
  !     ksi13 = l1**2*l3 + phi*0.5_f64
  !     ksi21 = l2**2*l1 + phi*0.5_f64
  !     ksi23 = l2**2*l3 + phi*0.5_f64
  !     ksi31 = l3**2*l1 + phi*0.5_f64
  !     ksi32 = l3**2*l2 + phi*0.5_f64

  !     base(1) = 3.0_f64 * l1**2 - 2.0_f64*ksi1 - 9.0_f64*phi
  !     base(2) = 3.0_f64 * l2**2 - 2.0_f64*ksi2 - 9.0_f64*phi
  !     base(3) = 3.0_f64 * l3**2 - 2.0_f64*ksi3 - 9.0_f64*phi
  !     base(4) = 27.0_f64 * phi
  !     base(5) = step * ( ksi12 - 1.5_f64*phi )
  !     base(6) = step * ( ksi13 - 1.5_f64*phi )
  !     base(7) = step * ( ksi21 - 1.5_f64*phi )
  !     base(8) = step * ( ksi23 - 1.5_f64*phi )
  !     base(9) = step * ( ksi31 - 1.5_f64*phi )
  !     base(10)= step * ( ksi32 - 1.5_f64*phi )

  !     f = 0.0_f64

  !     do i = 1,10
  !        f = f + freedom(i)*base(i)
  !     enddo

  !     f_tn1(num) = f

  !   end subroutine hermite_interpolation


end program test_hex_hermite














