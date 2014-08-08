module interpolation_hex_hermite
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use mesh_hex_alt
  use sll_hermite_interpolation_2d_module
  implicit none

contains


  subroutine der_finite_difference( f_tn, p,step, mesh_numh, mesh_coor, deriv, n_cell, dbc ) 
    !-> computation of the partial derivatives in the directions H1, H2 and H3
    implicit none
    sll_int32 , dimension(:,:)             :: mesh_coor, mesh_numh
    sll_real64, dimension(:,:)             :: deriv 
    sll_real64, dimension(:), intent(in)   :: f_tn 
    sll_real64, dimension(:,:), intent(in) :: dbc
    sll_real64, intent(in)                 :: step
    sll_int32, intent(in)                  :: n_cell, p
    sll_int32                              :: i, j, r, s, n_points, nc1
    sll_int32                              :: n1, n2, n3, n4, n5, n6
    sll_real64                             :: dh1, dh2, dh3, dh4, dh5, dh6
    sll_int32                              :: h1 ,h2, h1t ,h2t, p2
    sll_int32                              :: h1n ,h2n, h1tn ,h2tn
    logical                                :: opsign1, opsign2, opsign3
    sll_real64, dimension(:),allocatable   :: w
    sll_real64, dimension(:,:),allocatable :: f


    n_points = size(f_tn)

    !*********************************************************************
    !computation of the coefficients in fonction of the degree p 
    !of the approximation required
    !*********************************************************************

    p2 = p/2

    if ( (p2)*2 == p ) then !if p is even
       r = -p2       
       s = p2
    else  !if p is odd
       r = -p2 - 1 ! de-centered to the left      
       s = p2 !+ 1 ! de-centered to the right
    endif

    allocate( w(r:s),f(1:6,r:s) )

    call compute_w_hermite(w,r,s) 

    !***********************************************************************
    ! Computation of the derivatives in both hexa directions h1 and h2,
    ! then for both x and y directions for every points of the mesh
    !***********************************************************************

    nc1 = n_cell + 1

    do i = 1, n_points 

       h1 = mesh_numh(1,i) 
       h2 = mesh_numh(2,i)  
       h1n = h1 + nc1
       h2n = h2 + nc1


       do j = r,s !find the required points in the h1, h2 and h3 directions 

          opsign1 = .false.
          opsign2 = .false.
          opsign3 = .false.

          h1t = h1 + j 
          h2t = h2 + j
          h1tn = h1n + j 
          h2tn = h2n + j
          
          ! if (le point est sur la frontière) on retient ses coordonnées hex
          ! 

          if ( h2  > 0 .and. h1t < 0 .or.  h1t > 0 .and. h2  < 0 ) opsign1 = .true.
          if ( h2t > 0 .and. h1  < 0 .or.  h1  > 0 .and. h2t < 0 ) opsign2 = .true.
          if ( h2t > 0 .and. h1t < 0 .or.  h1t > 0 .and. h2t < 0 ) opsign3 = .true.

          !test if outside mesh 

          if (  abs(h1t) > n_cell .or. opsign1 .and. ( abs(h1t) + abs(h2) > n_cell) ) then
             f(1,j) = dbc(1,1)!f_tn(nf) ! dirichlet boundary condition
          else
             n1 = mesh_coor(h1tn,h2n)
             f(1,j) = f_tn(n1) 
          endif


          if (  abs(h2t) > n_cell .or. opsign2 .and. ( abs(h1) + abs(h2t) > n_cell) ) then 
             f(2,j) = dbc(2,1) 
          else
             n2 = mesh_coor(h1n,h2tn)
             f(2,j) = f_tn(n2) 
          endif

          if ( abs(h1t) > n_cell .or. abs(h2t) > n_cell .or. opsign3 .and. ( abs(h1t) + abs(h2t) > n_cell) ) then
             f(3,j) = dbc(3,1) 
          else
             n3 = mesh_coor(h1tn,h2tn)
             f(3,j) = f_tn(n3) 
          endif

          !find the required points in the -h1, -h2 and -h3 directions 

          opsign1 = .false.
          opsign2 = .false.
          opsign3 = .false.

          h1t = h1 - j
          h2t = h2 - j
          h1tn = h1n - j
          h2tn = h2n - j

          if ( h2  > 0 .and. h1t < 0 .or.  h1t > 0 .and. h2  < 0 ) opsign1 = .true.
          if ( h2t > 0 .and. h1  < 0 .or.  h1  > 0 .and. h2t < 0 ) opsign2 = .true.
          if ( h2t > 0 .and. h1t < 0 .or.  h1t > 0 .and. h2t < 0 ) opsign3 = .true.

          !test if outside mesh 

          if (  abs(h1t) > n_cell .or. opsign1 .and. ( abs(h1t) + abs(h2) > n_cell) ) then
             f(4,j) = dbc(4,1) ! dirichlet boundary condition
          else
             n4 = mesh_coor(h1tn,h2n)
             f(4,j) = f_tn(n4) 
          endif

          if (  abs(h2t) > n_cell .or. opsign2 .and. ( abs(h1) + abs(h2t) > n_cell) ) then
             f(5,j) = dbc(5,1) 
          else
             n5 = mesh_coor(h1n,h2tn) 
             f(5,j) = f_tn(n5) 
          endif

          if ( abs(h1t) > n_cell .or. abs(h2t) > n_cell .or. opsign3 .and. ( abs(h1t) + abs(h2t) > n_cell) ) then
             f(6,j) = dbc(6,1) 
          else
             n6 = mesh_coor(h1tn,h2tn)
             f(6,j) = f_tn(n6) 
          endif

       enddo


       dh1 = 0._f64
       dh2 = 0._f64
       dh3 = 0._f64
       dh4 = 0._f64
       dh5 = 0._f64
       dh6 = 0._f64

       do j = r,s
          dh1 = dh1 + w(j) * f(1,j)  
          dh2 = dh2 + w(j) * f(2,j)
          dh3 = dh3 + w(j) * f(3,j)
          dh4 = dh4 + w(j) * f(4,j)  
          dh5 = dh5 + w(j) * f(5,j)
          dh6 = dh6 + w(j) * f(6,j)
       enddo

       deriv(1,i) = dh1/step
       deriv(2,i) = dh2/step
       deriv(3,i) = dh3/step
       deriv(4,i) = dh4/step
       deriv(5,i) = dh5/step
       deriv(6,i) = dh6/step

    enddo

    deallocate(f,w)

  end subroutine der_finite_difference




  subroutine hermite_interpolation(num, x, y, f_tn, f_tn1, mesh_num, mesh_coor, deriv, step, aire, num_cells ,radius) 

    implicit none
    sll_int32 , dimension(:,:) :: mesh_coor
    sll_real64, dimension(:,:) :: mesh_num, deriv 
    sll_real64,intent(in)      :: x, y, step, aire, radius
    sll_int32,intent(in)       :: num_cells, num
    sll_real64                 :: x1, x2, x3, y1, y2, y3, f
    sll_real64                 :: phi, ksi1, ksi2, ksi3, l1, l2, l3, eps = 1.d-12
    sll_real64                 :: ksi12, ksi13, ksi21, ksi23, ksi31, ksi32
    sll_real64,dimension(:), intent(in)    :: f_tn
    sll_real64,dimension(:), intent(out)   :: f_tn1
    sll_real64,dimension(1:10) :: freedom, base
    sll_int32                  :: i, i1, i2, i3
    sll_real64                :: a2, p05, p15, p9, l12, l22, l32               

    ! find the triangle where (x,y) belongs to
    ! i1, i2, i3 are the indices for the S1, S2 and S3 vertexes of the triangle
    ! with S1 the vertex with the smallest y , S2 intermediate and S3 biggest
    ! in other words , s1 bottom , s2 side, s3 top 

    call search_tri( x, y, mesh_coor,step, num_cells, radius, i1, i2, i3 )


    ! get the ten degrees of freedom

    ! values at the vertexes of the triangle

    freedom(1) = f_tn(i1)
    freedom(2) = f_tn(i2)
    freedom(3) = f_tn(i3)

    ! value at the center of the triangle

    freedom(4) =  ( freedom(1) + freedom(2) + freedom(3) ) / 3.00

    x1 = mesh_num(1,i1) 
    x2 = mesh_num(1,i2)
    x3 = mesh_num(1,i3)
    y1 = mesh_num(2,i1)
    y2 = mesh_num(2,i2)
    y3 = mesh_num(2,i3)

    ! values of the derivatives

    freedom(5) = deriv(3,i1) ! derivative from S1 to S3
    freedom(9) = deriv(6,i3) ! derivative from S3 to S1

    !if ( (/x3-x2,y3-y2/) == h1 )then      ! if S2 S3 = h1 then we are in the first case )
    if ( ((x3-x2) - sqrt(3.0_f64)*0.5_f64)**2 + (y3-y2 - 0.5_f64)**2 <= eps ) then 
       freedom(6) = deriv(2,i1) ! derivative from S1 to S2
       freedom(7) = deriv(5,i2) ! derivative from S2 to S1
       freedom(8) = deriv(1,i2) ! derivative from S2 to S3
       freedom(10)= deriv(4,i3) ! derivative from S3 to S2
    else !if ( (/x3-x2,y3-y2/) == h2 )! if S2 S3 = h2 then we are in the 2nd case )
       freedom(6) = deriv(1,i1) ! derivative from S1 to S2 
       freedom(7) = deriv(4,i2) ! derivative from S2 to S1
       freedom(8) = deriv(2,i2) ! derivative from S2 to S3
       freedom(10)= deriv(5,i3) ! derivative from S3 to S2
    endif

    a2 = 0.5_f64/aire

    l1   = a2 * abs( (x2 - x)*(y3 - y) - (x3 - x)*(y2 - y) ) 
    l2   = a2 * abs( (x1 - x)*(y3 - y) - (x3 - x)*(y1 - y) )  
    l3   = a2 * abs( (x1 - x)*(y2 - y) - (x2 - x)*(y1 - y) )  

    phi = l1*l2*l3

    ksi1 = l1**3 - phi
    ksi2 = l2**3 - phi
    ksi3 = l3**3 - phi

    ! optimization variable
    p05 = phi*0.5_f64
    p9  = 9.0_f64*phi
    p15 = 1.5_f64*phi
    l12 = l1**2
    l22 = l2**2
    l32 = l3**2

    ksi12 = l12*l2 + p05
    ksi13 = l12*l3 + p05
    ksi21 = l22*l1 + p05
    ksi23 = l22*l3 + p05
    ksi31 = l32*l1 + p05
    ksi32 = l32*l2 + p05

    base(1) = 3.0_f64 * l12 - 2.0_f64*ksi1 - p9
    base(2) = 3.0_f64 * l22 - 2.0_f64*ksi2 - p9
    base(3) = 3.0_f64 * l32 - 2.0_f64*ksi3 - p9
    base(4) = 27.0_f64 * phi
    base(5) = step * ( ksi12 - p15 )
    base(6) = step * ( ksi13 - p15 )
    base(7) = step * ( ksi21 - p15 )
    base(8) = step * ( ksi23 - p15 )
    base(9) = step * ( ksi31 - p15 )
    base(10)= step * ( ksi32 - p15 )

    f = 0.0_f64

    do i = 1,10
       f = f + freedom(i)*base(i)
    enddo

    f_tn1(num) = f

  end subroutine hermite_interpolation


end module interpolation_hex_hermite
