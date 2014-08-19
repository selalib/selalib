module interpolation_hex_hermite
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use hex_mesh
  use sll_hermite_interpolation_2d_module
  implicit none

contains


  subroutine der_finite_difference( f_tn, p,step, mesh, deriv ) 
    !-> computation of the partial derivatives in the directions H1, H2 and H3
    ! with dirichlet boundary condition
    implicit none
    type(hex_mesh_2d), pointer             :: mesh
    sll_real64, dimension(:,:), intent(out):: deriv 
    sll_real64, dimension(:), intent(in)   :: f_tn 
    sll_real64, intent(in)                 :: step
    sll_int32,  intent(in)                 :: p
    sll_int32                              :: num_cells, i, j, r, s, n_points
    sll_int32                              :: n1, n2, n3, n4, n5, n6, nc1
    sll_real64                             :: dh1, dh2, dh3, dh4, dh5, dh6
    sll_int32                              :: h1, h2, h1t, h2t, p2
    sll_int32                              :: index_boundary, index_boundary1, index_boundary2
    sll_int32                              :: index_boundary3, index_boundary4, index_boundary5
    sll_int32                              :: index_boundary6, maxrs
    logical                                :: opsign1, opsign2, opsign3
    logical                                :: boundary,boundary1,boundary2,boundary3
    logical                                :: boundary4,boundary5,boundary6
    logical                                :: secure
    sll_real64, dimension(:),allocatable   :: w
    sll_real64, dimension(:,:),allocatable :: f

    n_points = size(f_tn)
    num_cells = mesh%num_cells

    !*********************************************************************
    ! computation of the coefficients in fonction of the degree p 
    ! of the approximation required
    !*********************************************************************

    p2 = p/2

    if ( (p2)*2 == p ) then !if p is even
       r = -p2       
       s = p2
    else  !if p is odd
       r = -p2 - 1 ! de-centered to the left      
       s =  p2 !+ 1 ! de-centered to the right
    endif

    allocate( w(r:s),f(1:6,r:s) )

    call compute_w_hermite(w,r,s) 

    !***********************************************************************
    ! Computation of the derivatives in both hexa directions h1 and h2,
    ! then for both x and y directions for every points of the mesh
    !***********************************************************************
    
    maxrs = max(abs(r),s)
    nc1   = num_cells + 1

    do i = 1, n_points 

       secure = .false.
       h1 = mesh%hex_coord(1,i)
       h2 = mesh%hex_coord(2,i)
       
       ! checking if the point "i" is in the "secure zone" 
       ! where boundary conditions won't play a role 

       if ( h1 <= num_cells - maxrs .and. h2 <= num_cells - maxrs  & 
            .and. h1 >= - num_cells + maxrs .and. h2 >= - num_cells + maxrs ) then
          if ( h1 <= 0 .and. h2 >= 0 .and. h2 - h1 <= num_cells - maxrs  & 
               .or. h2 <= 0 .and. h1 >= 0 .and. h1 - h2 <= num_cells - maxrs   & 
               ) then 
             secure = .true.
          endif
       endif



       if (secure) then 

          do j = r,s !find the required points in the h1, h2 and h3 directions 

             h1t = h1 + j 
             h2t = h2 + j

             n1 = hex_to_global(mesh,h1t,h2)
             f(1,j) = f_tn(n1) 

             n2 = hex_to_global(mesh,h1,h2t)
             f(2,j) = f_tn(n2) 

             n3 = hex_to_global(mesh,h1t,h2t)
             f(3,j) = f_tn(n3) 

             !find the required points in the -h1, -h2 and -h3 directions 

             h1t = h1 - j
             h2t = h2 - j

             n4 = hex_to_global(mesh,h1t,h2)
             f(4,j) = f_tn(n4) 

             n5 = hex_to_global(mesh,h1,h2t)
             f(5,j) = f_tn(n5) 

             n6 = hex_to_global(mesh,h1t,h2t)
             f(6,j) = f_tn(n6) 

          enddo


       else


          ! checking if the point "i" is on the boundary:

          if ( h1 == num_cells .and. h2 >= 0 .and. h2 <= num_cells  & !edge 1
               .or. h1 == -num_cells .and. h2 <= 0 .and. h2 >= -num_cells  & !edge4
               .or. h2 ==  num_cells .and. h1 >= 0 .and. h1 <=  num_cells  & !edge2
               .or. h2 == -num_cells .and. h1 <= 0 .and. h1 >= -num_cells  & !edge5
               .or. h1 <= 0 .and. h2 >= 0 .and. h2 - h1 == num_cells  & !edge3
               .or. h2 <= 0 .and. h1 >= 0 .and. h1 - h2 == num_cells  & !edge6
               ) then
             boundary       = .true.
             index_boundary = i
          endif

          ! we could have kept in memory which edge it is on, but due to simplicity concern 
          ! we did not implement a test for each edge even though it would save tests
          ! it could be done for optimisation purpose in the future

          if (boundary) then 

             do j = r,s !find the required points in the h1, h2 and h3 directions 

                opsign1 = .false.
                opsign2 = .false.
                opsign3 = .false.

                h1t = h1 + j 
                h2t = h2 + j

                if ( h2*h1t < 0  ) opsign1 = .true.
                if ( h2t*h1 < 0  ) opsign2 = .true.
                if ( h2t*h1t < 0 ) opsign3 = .true.

                !test if outside mesh 

                if (  abs(h1t) > num_cells .or. opsign1 .and. ( abs(h1t) + abs(h2) > num_cells) ) then
                   f(1,j) = f_tn(index_boundary) ! dirichlet boundary condition
                else
                   n1 = hex_to_global(mesh,h1t,h2)
                   f(1,j) = f_tn(n1) 
                endif


                if (  abs(h2t) > num_cells .or. opsign2 .and. ( abs(h1) + abs(h2t) > num_cells) ) then 
                   f(2,j) = f_tn(index_boundary) 
                else
                   n2 = hex_to_global(mesh,h1,h2t)
                   f(2,j) = f_tn(n2) 
                endif

                if ( abs(h1t) > num_cells .or. abs(h2t) > num_cells .or. opsign3 .and. ( abs(h1t) + abs(h2t) > num_cells) ) then
                   f(3,j) = f_tn(index_boundary) 
                else
                   n3 = hex_to_global(mesh,h1t,h2t)
                   f(3,j) = f_tn(n3) 
                endif

                !find the required points in the -h1, -h2 and -h3 directions 

                opsign1 = .false.
                opsign2 = .false.
                opsign3 = .false.

                h1t = h1 - j
                h2t = h2 - j

                if ( h2*h1t < 0  ) opsign1 = .true.
                if ( h2t*h1 < 0  ) opsign2 = .true.
                if ( h2t*h1t < 0 ) opsign3 = .true.

                !test if outside mesh 

                if (  abs(h1t) > num_cells .or. &
                     opsign1 .and. ( abs(h1t) + abs(h2) > num_cells) ) then
                   f(4,j) = f_tn(index_boundary) ! dirichlet boundary condition
                else
                   n4 = hex_to_global(mesh,h1t,h2)
                   f(4,j) = f_tn(n4) 
                endif

                if (  abs(h2t) > num_cells .or. &
                     opsign2 .and. ( abs(h1) + abs(h2t) > num_cells) ) then
                   f(5,j) = f_tn(index_boundary) 
                else
                   n5 = hex_to_global(mesh,h1,h2t)
                   f(5,j) = f_tn(n5) 
                endif

                if ( abs(h1t) > num_cells .or. abs(h2t) > num_cells .or. &
                   opsign3 .and. ( abs(h1t) + abs(h2t) > num_cells) ) then
                   f(6,j) = f_tn(index_boundary)
                else
                   n6 = hex_to_global(mesh,h1t,h2t)
                   f(6,j) = f_tn(n6) 
                endif

             enddo

          else !case where there will be boundary points but "i" is not one of them 

             f(1:6,0) = f_tn(i)

             boundary1 = .false.
             boundary2 = .false.
             boundary3 = .false.
             boundary4 = .false.
             boundary5 = .false.
             boundary6 = .false.

             do j = 1,s

                h1t = h1 + j 
                h2t = h2 + j

                ! we need to stock the index of the boundary (ies) point(s) in any direction
                ! then test  
                
                if ( boundary1 .eqv. .false.) then
                   if ( h1t == num_cells .and. h2 >= 0 .and. h2 <= num_cells .or. &
                        h2 <= 0 .and. h1t >= 0 .and. h1t - h2 == num_cells) then
                      boundary1       = .true.
                      index_boundary1 = hex_to_global(mesh,h1t,h2)
                   endif
                endif

                if ( boundary2 .eqv. .false.) then
                   if (h2t ==  num_cells .and. h1 >= 0 .and. h1 <=  num_cells .or. &
                        h1 <= 0 .and. h2t >= 0 .and. h2t - h1 == num_cells) then
                      boundary2       = .true.
                      index_boundary2 = hex_to_global(mesh,h1,h2t)
                   endif
                endif

                if ( boundary3 .eqv. .false.) then
                   if (h1t == num_cells .and. h2t >= 0 .and. h2t <= num_cells .or. &
                        h2t ==  num_cells .and. h1t >= 0 .and. h1t <=  num_cells) then
                      boundary3       = .true.
                      index_boundary3 = hex_to_global(mesh,h1t,h2t)
                   endif
                endif

                
                if ( boundary1 ) then
                   f(1,j) = f_tn(index_boundary1) ! dirichlet boundary condition
                else
                   n1 = hex_to_global(mesh,h1t,h2)
                   f(1,j) = f_tn(n1) 
                endif

                if ( boundary2 ) then 
                   f(2,j) = f_tn(index_boundary2) 
                else
                   n2 = hex_to_global(mesh,h1,h2t)
                   f(2,j) = f_tn(n2) 
                endif

                if ( boundary3 ) then
                   f(3,j) = f_tn(index_boundary3) 
                else
                   n3 = hex_to_global(mesh,h1t,h2t)
                   f(3,j) = f_tn(n3) 
                endif

                h1t = h1 - j
                h2t = h2 - j

                if ( boundary4 .eqv. .false.) then
                   if (h1t == -num_cells .and. h2 <= 0 .and. h2 >= -num_cells .or. &
                        h1t <= 0 .and. h2 >= 0 .and. h2 - h1t == num_cells) then
                      boundary4       = .true.
                      index_boundary4 = hex_to_global(mesh,h1t,h2)
                   endif
                endif

                if ( boundary5 .eqv. .false.) then
                   if (h2t == -num_cells .and. h1 <= 0 .and. h1 >= -num_cells .or. &
                        h2t <= 0 .and. h1 >= 0 .and. h1 - h2t == num_cells) then
                      boundary5       = .true.
                      index_boundary5 = hex_to_global(mesh,h1,h2t)
                   endif
                endif

                if ( boundary6 .eqv. .false.) then
                   if (h1t == -num_cells .and. h2t <= 0 .and. h2t >= -num_cells .or. &
                        h2t == -num_cells .and. h1t <= 0 .and. h1t >= -num_cells) then
                      boundary6       = .true.
                      index_boundary6 = hex_to_global(mesh,h1t,h2t)
                   endif
                endif


                if ( boundary4 ) then
                   f(4,j) = f_tn(index_boundary4) ! dirichlet boundary condition
                else
                   n4 = hex_to_global(mesh,h1t,h2)
                   f(4,j) = f_tn(n4) 
                endif

                if ( boundary5 ) then 
                   f(5,j) = f_tn(index_boundary5) 
                else
                   n5 = hex_to_global(mesh,h1,h2t)
                   f(5,j) = f_tn(n5) 
                endif

                if ( boundary6 ) then
                   f(6,j) = f_tn(index_boundary6) 
                else
                   n6 = hex_to_global(mesh,h1t,h2t)
                   f(6,j) = f_tn(n6) 
                endif

             enddo

             boundary1 = .false.
             boundary2 = .false.
             boundary3 = .false.
             boundary4 = .false.
             boundary5 = .false.
             boundary6 = .false.

             do j = -1,r,-1

                h1t = h1 + j
                h2t = h2 + j
                
                if ( boundary1 .eqv. .false.) then
                   if (h1t == -num_cells .and. h2 <= 0 .and. h2 >= -num_cells .or. &
                        h1t <= 0 .and. h2 >= 0 .and. h2 - h1t == num_cells) then
                      boundary1       = .true.
                      index_boundary1 = hex_to_global(mesh,h1t,h2)
                   endif
                endif
                
                if ( boundary2 .eqv. .false.) then
                   if (h2t == -num_cells .and. h1 <= 0 .and. h1 >= -num_cells .or. &
                        h2t <= 0 .and. h1 >= 0 .and. h1 - h2t == num_cells) then
                      boundary2       = .true.
                      index_boundary2 = hex_to_global(mesh,h1,h2t)
                   endif
                endif

                if ( boundary3 .eqv. .false.) then
                   if (h1t == -num_cells .and. h2t <= 0 .and. h2t >= -num_cells .or. &
                        h2t == -num_cells .and. h1t <= 0 .and. h1t >= -num_cells) then
                      boundary3       = .true.
                      index_boundary3 = hex_to_global(mesh,h1t,h2t)
                   endif
                endif

                
                if ( boundary1 ) then
                   f(1,j) = f_tn(index_boundary1) ! dirichlet boundary condition
                else
                   n1 = hex_to_global(mesh,h1t,h2)
                   f(1,j) = f_tn(n1) 
                endif

                if ( boundary2 ) then 
                   f(2,j) = f_tn(index_boundary2) 
                else
                   n2 = hex_to_global(mesh,h1,h2t)
                   f(2,j) = f_tn(n2) 
                endif

                if ( boundary3 ) then
                   f(3,j) = f_tn(index_boundary3) 
                else
                   n3 = hex_to_global(mesh,h1t,h2t)
                   f(3,j) = f_tn(n3) 
                endif

                h1t = h1 - j
                h2t = h2 - j

                if ( boundary4 .eqv. .false.) then
                   if  ( h1t == num_cells .and. h2 >= 0 .and. h2 <= num_cells .or. &
                        h2 <= 0 .and. h1t >= 0 .and. h1t - h2 == num_cells) then
                      boundary4       = .true.
                      index_boundary4 = hex_to_global(mesh,h1t,h2)
                   endif
                endif

                if ( boundary5 .eqv. .false.) then
                   if (h2t ==  num_cells .and. h1 >= 0 .and. h1 <=  num_cells .or. &
                        h1 <= 0 .and. h2t >= 0 .and. h2t - h1 == num_cells) then
                      boundary5       = .true.
                      index_boundary5 = hex_to_global(mesh,h1,h2t)
                   endif
                endif

                if ( boundary6 .eqv. .false.) then
                   if (h1t == num_cells .and. h2t >= 0 .and. h2t <= num_cells .or. &
                        h2t ==  num_cells .and. h1t >= 0 .and. h1t <=  num_cells) then
                      boundary6       = .true.
                      index_boundary6 = hex_to_global(mesh,h1t,h2t)
                   endif
                endif


                if ( boundary4 ) then
                   f(4,j) = f_tn(index_boundary4) ! dirichlet boundary condition
                else
                   n4 = hex_to_global(mesh,h1t,h2)
                   f(4,j) = f_tn(n4) 
                endif

                if ( boundary5 ) then 
                   f(5,j) = f_tn(index_boundary5) 
                else
                   n5 = hex_to_global(mesh,h1,h2t)
                   f(5,j) = f_tn(n5) 
                endif

                if ( boundary6 ) then
                   f(6,j) = f_tn(index_boundary6) 
                else
                   n6 = hex_to_global(mesh,h1t,h2t)
                   f(6,j) = f_tn(n6) 
                endif

             enddo

          endif

       endif


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




  subroutine hermite_interpolation(num, x, y, f_tn, center_value, f_tn1, mesh, deriv, aire,t, num_degree_freedom) 

    implicit none
    type(hex_mesh_2d), pointer             :: mesh
    sll_real64,dimension(:), intent(in)    :: f_tn
    sll_real64,dimension(:), intent(in)    :: center_value
    sll_real64, dimension(:,:), intent(in) :: deriv 
    sll_real64,dimension(:), intent(out)   :: f_tn1 
    sll_real64,intent(in)      :: x, y, aire,t
    sll_int32,intent(in)       :: num, num_degree_freedom
    sll_real64                 :: x1, x2, x3, y1, y2, y3, f, step
    sll_real64                 :: phi, ksi1, ksi2, ksi3, l1, l2, l3, eps = 1.d-8
    sll_real64                 :: ksi12, ksi13, ksi21, ksi23, ksi31, ksi32
    sll_real64,dimension(1:10) :: freedom, base
    sll_int32                  :: i, i1, i2, i3, k11, k12, center_index
    sll_real64                 :: a2, p05, p15, p9, l12, l22, l32  
    sll_real64                 :: x1x,x2x,x3x,y1y,y2y,y3y, xx,yy
    logical                    :: test

    ! find the triangle which (x,y) belongs to
    ! i1, i2, i3 are the indices for the S1, S2 and S3 vertexes of the triangle
    ! with S1 the vertex with the smallest y , S2 intermediate and S3 biggest
    ! in other words , s1 bottom , s2 side, s3 top 

    step = mesh%delta

    call get_cell_vertices_index( x, y, mesh, i1, i2, i3 )

    x1 = mesh%cartesian_coord(1,i1) 
    x2 = mesh%cartesian_coord(1,i2) 
    x3 = mesh%cartesian_coord(1,i3) 
    y1 = mesh%cartesian_coord(2,i1) 
    y2 = mesh%cartesian_coord(2,i2) 
    y3 = mesh%cartesian_coord(2,i3) 

    k11 = mesh%hex_coord(1,i1) 
    k12 = mesh%hex_coord(2,i1) 
    
    call get_triangle_index(k11,k12,mesh,x,center_index)

    ! get the ten degrees of freedom

    ! values at the vertexes of the triangle

    freedom(1) = f_tn(i1)
    freedom(2) = f_tn(i2)
    freedom(3) = f_tn(i3)

    ! value at the center of the triangle

    ! Value given as an average of the vertices' respective value  
    !freedom(4) =  ( freedom(1) + freedom(2) + freedom(3) ) / 3._f64

    ! value computed
    if (center_index == -1) then
       print *, "problem index center l 536 interpolation hermite"
    else
       freedom(4) = center_value(center_index)
    endif

    ! exact value : 

    ! xx = x2 + ( (x3+x1)*0.5_f64 - x2 )* 2.0_f64 / 3.0_f64
    ! yy = y2 + ( (y3+y1)*0.5_f64 - y2 )* 2.0_f64 / 3.0_f64

    ! xx = (x2 + ( (x3+x1)*0.5_f64 - x2 )* 2.0_f64 / 3.0_f64) * cos(t) &
    !      - (y2 + ( (y3+y1)*0.5_f64 - y2 )* 2.0_f64 / 3.0_f64 ) * sin(t)
    ! yy = (x2 + ( (x3+x1)*0.5_f64 - x2 )* 2.0_f64 / 3.0_f64) * sin(t) &
    !      + ( y2 + ( (y3+y1)*0.5_f64 - y2 )* 2.0_f64 / 3.0_f64 ) * cos(t)

    ! freedom(4) = exp( - ( (xx-1.0_f64)**2+(yy-1.0_f64)**2) ) 


    ! values of the derivatives

    freedom(6) = deriv(3,i1) ! derivative from S1 to S3 (+h3)
    freedom(9) = deriv(6,i3) ! derivative from S3 to S1 (-h3)

    test = .false.

    ! if S2 S3 = h1 then we are in the first case 
    if ( (x3-x2 - sqrt(3._f64)*0.5_f64*step)**2 + (y3-y2 - 0.5_f64*step)**2 <= eps ) then 
       freedom(5) = deriv(2,i1) ! derivative from S1 to S2 (+h2)
       freedom(7) = deriv(5,i2) ! derivative from S2 to S1 (-h2)
       freedom(8) = deriv(1,i2) ! derivative from S2 to S3 (+h1)
       freedom(10)= deriv(4,i3) ! derivative from S3 to S2 (-h1)

       test = .true.

       ! if S2 S3 = h2 then we are in the 2nd case )
    else if ( (x3-x2 + sqrt(3._f64)*0.5_f64*step)**2 + (y3-y2 - 0.5_f64*step)**2 <= eps ) then 
       freedom(5) = deriv(1,i1) ! derivative from S1 to S2 (+h1)
       freedom(7) = deriv(4,i2) ! derivative from S2 to S1 (-h1)
       freedom(8) = deriv(2,i2) ! derivative from S2 to S3 (+h2)
       freedom(10)= deriv(5,i3) ! derivative from S3 to S2 (-h2)

       test = .true.

    endif

    if (test .eqv. .false. ) print*, "anomaly in hermite_interpolation l289"

    a2  = 0.5_f64/aire
    y3y = y3 - y
    y2y = y2 - y
    y1y = y1 - y
    x3x = x3 - x
    x2x = x2 - x
    x1x = x1 - x

    l1   = a2 * abs( x2x*y3y - x3x*y2y )    ! barycentric coordinates
    l2   = a2 * abs( x1x*y3y - x3x*y1y ) 
    l3   = 1._f64 - l1 - l2

    phi = l1*l2*l3  ! "bubble function"

    ksi1 = l1**3 - phi
    ksi2 = l2**3 - phi
    ksi3 = l3**3 - phi

    ! optimization variables

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

    if (num_degree_freedom == 10 ) then

       ! Computing the ten canonical basis functions 

       base(1) = 3._f64 * l12 - 2._f64*ksi1 - p9
       base(2) = 3._f64 * l22 - 2._f64*ksi2 - p9
       base(3) = 3._f64 * l32 - 2._f64*ksi3 - p9
       base(4) = 27._f64 * phi
       base(5) = step * ( ksi12 - p15 )
       base(6) = step * ( ksi13 - p15 )
       base(7) = step * ( ksi21 - p15 )
       base(8) = step * ( ksi23 - p15 )
       base(9) = step * ( ksi31 - p15 )
       base(10)= step * ( ksi32 - p15 )

    else if (num_degree_freedom == 9 ) then

       ! Computing the nine canonical basis functions 

       base(1) = 3._f64 * l12 - 2._f64*ksi1 
       base(2) = 3._f64 * l22 - 2._f64*ksi2 
       base(3) = 3._f64 * l32 - 2._f64*ksi3 
       base(4) = 0._f64
       base(5) = step * ksi12 
       base(6) = step * ksi13
       base(7) = step * ksi21
       base(8) = step * ksi23
       base(9) = step * ksi31
       base(10)= step * ksi32

    else

       print*, "specify another degree of freedom : 9 or 10"

    endif

    f = 0._f64

    do i = 1,10
       f = f + freedom(i)*base(i)
    enddo

    f_tn1(num) = f

  end subroutine hermite_interpolation


end module interpolation_hex_hermite
