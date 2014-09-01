module interpolation_hex_hermite
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

  use hex_mesh
  use sll_hermite_interpolation_2d_module
  implicit none

contains


  subroutine der_finite_difference( f_tn, p, step, mesh, deriv ) 
    !-> computation of the partial derivatives in the directions H1, H2 and H3
    ! with dirichlet boundary condition
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

          boundary  = .false.

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




  subroutine hermite_interpolation(num, x, y, f_tn, center_value, edge_value, output_tn1, mesh, deriv, aire, num_method)

    type(hex_mesh_2d), pointer             :: mesh
    sll_real64,dimension(:), intent(in)    :: f_tn
    sll_real64,dimension(:), intent(in)    :: edge_value
    sll_real64,dimension(:), intent(in)    :: center_value
    sll_real64, dimension(:,:), intent(in) :: deriv 
    sll_real64,dimension(:), intent(out)   :: output_tn1 
    sll_real64,intent(in)      :: x, y, aire
    sll_int32,intent(in)       :: num, num_method
    sll_real64                 :: x1, x2, x3, y1, y2, y3, f, step
    sll_real64                 :: l1, l2, l3
    sll_real64,dimension(:),allocatable :: freedom, base
    sll_real64                 :: a2
    sll_real64                 :: x1x,x2x,x3x,y1y,y2y,y3y
    sll_int32                  :: i, i1, i2, i3, k11, k12, center_index
    sll_int32                  :: num_degree
    sll_int32                  :: edge_index1,edge_index2,edge_index3
    logical                    :: test

    if ( num_method ==9 ) then
       allocate(freedom(9),base(9))
       num_degree = 9 
    elseif ( num_method == 10 ) then
       allocate(freedom(10),base(10))
       num_degree = 10 
    elseif ( num_method == 11 ) then
       allocate(freedom(9),base(9))
       num_degree = 9 
    elseif ( num_method == 12 ) then
       allocate(freedom(12),base(12))
       num_degree = 12 
    elseif ( num_method == 15 ) then
       allocate(freedom(15),base(15))
       num_degree = 15 
    endif
    
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

    call get_edge_index(k11,k12,mesh,x,edge_index1,edge_index2,edge_index3)

    ! get the first 9 degrees of freedom

    ! values at the vertices of the triangle

    freedom(1) = f_tn(i1)
    freedom(2) = f_tn(i2)
    freedom(3) = f_tn(i3)

    ! values at the middle of the edges
    
    if (num_method == 15 ) then 
       freedom(12) = edge_value(edge_index1)
       freedom(11) = edge_value(edge_index2)
       freedom(10) = edge_value(edge_index3)
    endif

    ! values of the derivatives

    freedom(5) = deriv(3,i1) ! derivative from S1 to S3 (+h3)
    freedom(8) = deriv(6,i3) ! derivative from S3 to S1 (-h3)
    
    test = .false.

    if ( x2 <= x1 ) then ! first case : triangle oriented left

       freedom(4) = deriv(2,i1) ! derivative from S1 to S2 (+h2)
       freedom(6) = deriv(5,i2) ! derivative from S2 to S1 (-h2)
       freedom(7) = deriv(1,i2) ! derivative from S2 to S3 (+h1)
       freedom(9) = deriv(4,i3) ! derivative from S3 to S2 (-h1)

       if ( num_degree == 12 ) then
          call get_normal_der(deriv,i1,i2,mesh,freedom(10:12))
       elseif ( num_degree == 15 ) then
          call get_normal_der(deriv,i1,i2,mesh,freedom(13:15))
       endif

       test = .true.

    else if ( x2 > x1 ) then ! second case triangle oriented right 
       freedom(4) = deriv(1,i1) ! derivative from S1 to S2 (+h1)
       freedom(6) = deriv(4,i2) ! derivative from S2 to S1 (-h1)
       freedom(7) = deriv(2,i2) ! derivative from S2 to S3 (+h2)
       freedom(9) = deriv(5,i3) ! derivative from S3 to S2 (-h2)
       
       if ( num_degree == 12 ) then
          call get_normal_der(deriv,i1,i2,mesh,freedom(10:12))
       elseif ( num_degree == 15 ) then
          call get_normal_der(deriv,i1,i2,mesh,freedom(13:15))
       endif

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
    
    if (num_method == 9 ) then

       ! Computing the nine canonical basis functions 
       call base_zienkiewicz_9_degree_of_freedom(base,step,l1,l2,l3)

    else if (num_method == 10 ) then

       freedom(10) = center_value(center_index)

       ! Computing the ten canonical basis functions 
       call base_zienkiewicz_10_degree_of_freedom(base,step,l1,l2,l3)

    else if (num_method == 11 ) then 

       ! Computing the basis for the cubic element of Hsieh-Clough-Tocher-reduced 
 
       call base_Hsieh_Clough_Tocher_reduced&
            (base,x1,x3,y1,y2,y3,x,y,l1,l2,l3,step)

    else if (num_method == 12 ) then 

       ! Computing the basis for the cubic element of Hsieh-Clough-Tocher-complete
 
       call base_Hsieh_Clough_Tocher_complete&
            (base,x1,x3,y1,y2,y3,x,y,l1,l2,l3,step)

    else if (num_method == 15 ) then 

       ! Computing the basis for the quartic element of Ganev_Dimitrov

       call base_ganev_dimitrov(base,l1,l2,l3,step)

    endif

    f = 0._f64

    do i = 1,num_degree
       f = f + freedom(i)*base(i)
    enddo

    output_tn1(num) = f    
    
    deallocate(freedom,base)

  end subroutine hermite_interpolation

  

  subroutine base_zienkiewicz_9_degree_of_freedom(base,step,l1,l2,l3)
  sll_real64,dimension(:),intent(out) :: base
  sll_real64,intent(in)               :: l1, l2, l3 ! barycentric coord
  sll_real64,intent(in)               :: step
  sll_real64     :: phi, p05
  sll_real64     :: l12, l22, l32, ksi1, ksi2, ksi3 
  sll_real64     :: ksi12, ksi13, ksi21, ksi23, ksi31, ksi32 
  
  
  phi = l1*l2*l3  ! "bubble function"

  ksi1 = l1**3 - phi
  ksi2 = l2**3 - phi
  ksi3 = l3**3 - phi

  ! optimization variables

  p05 = phi*0.5_f64
  l12 = l1**2
  l22 = l2**2
  l32 = l3**2

  ksi12 = l12*l2 + p05
  ksi13 = l12*l3 + p05
  ksi21 = l22*l1 + p05
  ksi23 = l22*l3 + p05
  ksi31 = l32*l1 + p05
  ksi32 = l32*l2 + p05

  base(1) = 3._f64 * l12 - 2._f64*ksi1 
  base(2) = 3._f64 * l22 - 2._f64*ksi2 
  base(3) = 3._f64 * l32 - 2._f64*ksi3 
  base(4) = step * ksi12 
  base(5) = step * ksi13
  base(6) = step * ksi21
  base(7) = step * ksi23
  base(8) = step * ksi31
  base(9) = step * ksi32

  end subroutine base_zienkiewicz_9_degree_of_freedom


  subroutine base_zienkiewicz_10_degree_of_freedom(base,step,l1,l2,l3)
  sll_real64,dimension(:),intent(out) :: base
  sll_real64,intent(in)               :: l1, l2, l3 ! barycentric coord
  sll_real64,intent(in)               :: step
  sll_real64     :: phi, p05, p9, p15
  sll_real64     :: l12, l22, l32, ksi1, ksi2, ksi3 
  sll_real64     :: ksi12, ksi13, ksi21, ksi23, ksi31, ksi32 
  
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

  base(1) = 3._f64 * l12 - 2._f64*ksi1 - p9
  base(2) = 3._f64 * l22 - 2._f64*ksi2 - p9
  base(3) = 3._f64 * l32 - 2._f64*ksi3 - p9
  base(4) = step * ( ksi12 - p15 )
  base(5) = step * ( ksi13 - p15 )
  base(6) = step * ( ksi21 - p15 )
  base(7) = step * ( ksi23 - p15 )
  base(8) = step * ( ksi31 - p15 )
  base(9)= step * ( ksi32 - p15 )
  base(10) = 27._f64 * phi

  end subroutine base_zienkiewicz_10_degree_of_freedom


  
  subroutine get_sub_triangle_hex(x1,x3,y1,y2,y3,x,y,&
       num_sub_triangle)
    sll_real64,intent(in) :: x1, y1 ! cartesian coordinates of the lowest point
    sll_real64,intent(in) :: y2     ! ordinate of the middle point
    sll_real64,intent(in) :: x3, y3 ! cartesian coordinates of the highest point
    sll_real64,intent(in) :: x , y  ! cartesian coordinates of the point to be interpolated in a hexagonal mesh
    sll_int32,intent(out) :: num_sub_triangle


    ! what is the kind of the triangle T ?

    if ( x < x1 ) then ! first kind : oriented left 

       !finding the sub triangle K_l in which the interpolated point is

       if ( y > y2 ) then! K1 or K2
          if ( y >= ( x - x3 ) * sqrt(3._f64) + y3 ) then 
             num_sub_triangle = 1
          else
             num_sub_triangle = 2
          endif
       else ! K2 or K3
          if ( y >= - ( x - x1 ) * sqrt(3._f64) + y1 ) then 
             num_sub_triangle = 2
          else
             num_sub_triangle = 3
          endif
       endif

    else ! second kind : oriented right 

       !finding the triangle K_l in which the interpolated point is

       if ( y > y2 ) then ! K1 or K2 
          if ( y >= - ( x - x3 ) * sqrt(3._f64) + y3 ) then 
             num_sub_triangle = 1
          else
             num_sub_triangle = 2
          endif
       else  ! K2 or K3
          if ( y >= ( x - x1 ) * sqrt(3._f64) + y1 ) then 
             num_sub_triangle = 2
          else
             num_sub_triangle = 3
          endif
       endif

    endif

  end subroutine get_sub_triangle_hex

  !*******************************************************************************
  !                           subroutines for hctr
  !*******************************************************************************

  subroutine product_with_sigma_hct_r(li,lj,lk,xi)
    sll_real64,dimension(:),intent(out) :: xi ! local base
    sll_real64,intent(in)               :: li, lj, lk ! barycentric coord
    ! optimisation variables 
    sll_real64                          :: li2j, li2k, lj2i
    sll_real64                          :: lj2k, lk2i, lk2j
    sll_real64                          :: li3, lj3, lk3
    sll_real64                          :: lijk

    li2k = li**2*lk 
    li2j = li**2*lj 
    lj2i = lj**2*li
    lj2k = lj**2*lk
    lk2i = lk**2*li
    lk2j = lk**2*lj
    li3  = li**3
    lj3  = lj**3
    lk3  = lk**3
    lijk = li*lj*lk

    xi(1) = 4.5_f64 * ( li2k + li2j ) 
    xi(2) = 0.5_f64*li3 + lj3 - 1.5_f64*li2k + 3.0_f64*(lj2i + lj2k + lijk )
    xi(3) = 0.5_f64*li3 + lk3 - 1.5_f64*li2j + 3.0_f64*(lk2j + lk2i + lijk )
    xi(4) = - 0.25_f64*li3 + 1.25_f64 * li2k + 0.5_f64*li2j
    xi(5) = - 0.25_f64*li3 + 0.5_f64 * li2k + 1.25_f64*li2j
    xi(6) =   0.25_f64*li3 - 0.5_f64 * li2k - 0.25_f64*li2j + lj2i + lijk 
    xi(7) = - 0.25_f64*li2k + 0.25_f64 * li2j + lj2k + 0.5_f64*lijk 
    xi(8) =   0.25_f64*li2k - 0.25_f64 * li2j + lk2j + 0.5_f64*lijk
    xi(9) =   0.25_f64*li3 - 0.25_f64*li2k - 0.5_f64*li2j + lk2i + lijk 

  end subroutine product_with_sigma_hct_r

  
  subroutine base_from_local_base_xi_htc_r(base,xi,num_sub_triangle,step)
  sll_real64,dimension(:),intent(out) :: base
  sll_real64,dimension(:),intent(in)  :: xi 
  sll_real64,             intent(in)  :: step
  sll_int32,              intent(in)  :: num_sub_triangle

  
  if ( num_sub_triangle == 1 ) then 
     base(1) = xi(1) 
     base(2) = xi(2) 
     base(3) = xi(3) 
     base(4) = xi(5)*step
     base(5) = xi(4)*step
     base(6) = xi(6)*step
     base(7) = xi(7)*step
     base(8) = xi(9)*step
     base(9) = xi(8)*step
  elseif ( num_sub_triangle == 2 ) then 
     base(1) = xi(3) 
     base(2) = xi(1) 
     base(3) = xi(2) 
     base(4) = xi(9)*step
     base(5) = xi(8)*step
     base(6) = xi(4)*step
     base(7) = xi(5)*step
     base(8) = xi(7)*step
     base(9) = xi(6)*step
  elseif ( num_sub_triangle == 3 ) then 
     base(1) = xi(2) 
     base(2) = xi(3) 
     base(3) = xi(1) 
     base(4) = xi(7)*step
     base(5) = xi(6)*step
     base(6) = xi(8)*step
     base(7) = xi(9)*step
     base(8) = xi(5)*step
     base(9) = xi(4)*step
  endif

  end subroutine base_from_local_base_xi_htc_r
  
  



  subroutine base_Hsieh_Clough_Tocher_reduced(base,x1,x3,y1,y2,y3,x,y,l1,l2,l3,step)
    sll_real64,dimension(:),intent(out) :: base
    sll_real64,intent(in)               :: x1,x3,y1,y2,y3
    sll_real64,intent(in)               :: x,y,l1,l2,l3,step
    sll_real64,dimension(:),allocatable :: xi 
    sll_real64                          :: li, lj, lk
    sll_int32                           :: num_sub_triangle

    !finding the sub triangle K_l in which the interpolated point is

    call get_sub_triangle_hex(x1,x3,y1,y2,y3,x,y,&
         num_sub_triangle)

    allocate(xi(1:9))

    if ( num_sub_triangle == 1 ) then 
       li = l1
       lj = l2
       lk = l3
    else if ( num_sub_triangle == 2 ) then 
       li = l2
       lj = l3
       lk = l1
    else if ( num_sub_triangle == 3 ) then 
       li = l3
       lj = l1
       lk = l2 
    else 
       print*, "problem with finding the sub triangle"
    endif

    call product_with_sigma_hct_r(li,lj,lk,xi)
    call base_from_local_base_xi_htc_r(base,xi,num_sub_triangle,step)

    deallocate(xi)

  end subroutine base_Hsieh_Clough_Tocher_reduced
 
  !*******************************************************************************
  !                           subroutines for hctc
  !*******************************************************************************


  subroutine product_with_sigma_hct_c(li,lj,lk,xi)
    sll_real64,dimension(:),intent(out) :: xi ! local base
    sll_real64,intent(in)               :: li, lj, lk ! barycentric coord
    ! optimisation variables 
    sll_real64                          :: li2j, li2k, lj2i, lj2k, lk2i, lk2j
    sll_real64                          :: li3, lj3, lk3
    sll_real64                          :: lijk

    li2k = li**2*lk 
    li2j = li**2*lj 
    lj2i = lj**2*li
    lj2k = lj**2*lk
    lk2i = lk**2*li
    lk2j = lk**2*lj
    li3  = li**3
    lj3  = lj**3
    lk3  = lk**3
    lijk = li*lj*lk

    xi(1) = 4.5_f64 * ( li2k + li2j ) 
    xi(2) = 0.5_f64*li3 + lj3 - 1.5_f64*li2k + 3.0_f64*(lj2i + lj2k + lijk )
    xi(3) = 0.5_f64*li3 + lk3 - 1.5_f64*li2j + 3.0_f64*(lk2j + lk2i + lijk )
    xi(4) = - 1._f64/12._f64*li3 + 1.75_f64 * li2k - 0.5_f64*li2j
    xi(5) = - 1._f64/12._f64*li3 - 0.5_f64 * li2k + 1.75_f64*li2j
    xi(6) = - 7._f64/12._f64*li3 + 0.5_f64 * li2k + 1.25_f64*li2j + lj2i - lijk 
    xi(7) =   2._f64/3._f64*li3 - 0.75_f64 * li2k - 1.25_f64*li2j + lj2k + 1.5_f64*lijk 
    xi(8) =   2._f64/3._f64*li3 - 1.25_f64 * li2k - 0.75_f64*li2j + lk2j + 1.5_f64*lijk
    xi(9) = - 7._f64/12._f64*li3 + 1.25_f64* li2k + 0.5_f64 *li2j + lk2i - lijk 
    xi(10) =   4._f64/3._f64*li3 - 2._f64*li2k - 2._f64 * li2j + 4._f64*lijk 
    xi(11) = - 2._f64/3._f64*li3 + 2._f64*li2k
    xi(12) = - 2._f64/3._f64*li3 + 2._f64*li2j

  end subroutine product_with_sigma_hct_c

  
  subroutine base_from_local_base_xi_htc_c(base,xi,num_sub_triangle,step)
  sll_real64,dimension(:),intent(out) :: base
  sll_real64,dimension(:),intent(in)  :: xi 
  sll_real64,             intent(in)  :: step
  sll_int32,              intent(in)  :: num_sub_triangle
  sll_real64                          :: step_sq       

  step_sq = step*sqrt(3._f64)*0.5_f64
  
  if ( num_sub_triangle == 1 ) then 
     base(1) = xi(1) 
     base(2) = xi(2) 
     base(3) = xi(3) 
     base(4) = xi(5)*step
     base(5) = xi(4)*step
     base(6) = xi(6)*step
     base(7) = xi(7)*step
     base(8) = xi(9)*step
     base(9) = xi(8)*step
     base(10) = xi(10)*step_sq
     base(11) = xi(11)*step_sq
     base(12) = xi(12)*step_sq
  elseif ( num_sub_triangle == 2 ) then 
     base(1) = xi(3) 
     base(2) = xi(1) 
     base(3) = xi(2) 
     base(4) = xi(9)*step
     base(5) = xi(8)*step
     base(6) = xi(4)*step
     base(7) = xi(5)*step
     base(8) = xi(7)*step
     base(9) = xi(6)*step
     base(10) = xi(12)*step_sq
     base(11) = xi(10)*step_sq
     base(12) = xi(11)*step_sq
  elseif ( num_sub_triangle == 3 ) then 
     base(1) = xi(2) 
     base(2) = xi(3) 
     base(3) = xi(1) 
     base(4) = xi(7)*step
     base(5) = xi(6)*step
     base(6) = xi(8)*step
     base(7) = xi(9)*step
     base(8) = xi(5)*step
     base(9) = xi(4)*step
     base(10) = xi(11)*step_sq
     base(11) = xi(12)*step_sq
     base(12) = xi(10)*step_sq
  endif

  end subroutine base_from_local_base_xi_htc_c

  
  subroutine base_Hsieh_Clough_Tocher_complete(base,x1,x3,y1,y2,y3,x,y,l1,l2,l3,step)
    sll_real64,dimension(:),intent(out) :: base
    sll_real64,intent(in)               :: x1,x3,y1,y2,y3
    sll_real64,intent(in)               :: x,y,l1,l2,l3,step
    sll_real64,dimension(:),allocatable :: xi 
    sll_real64                          :: li, lj, lk
    sll_int32                           :: num_sub_triangle

    !finding the sub triangle K_l in which the interpolated point is

    call get_sub_triangle_hex(x1,x3,y1,y2,y3,x,y,&
         num_sub_triangle)

    allocate(xi(1:12))

    if ( num_sub_triangle == 1 ) then 
       li = l1
       lj = l2
       lk = l3
    else if ( num_sub_triangle == 2 ) then 
       li = l2
       lj = l3
       lk = l1
    else if ( num_sub_triangle == 3 ) then 
       li = l3
       lj = l1
       lk = l2 
    else 
       print*, "problem with finding the sub triangle"
    endif

    call product_with_sigma_hct_c(li,lj,lk,xi)
    call base_from_local_base_xi_htc_c(base,xi,num_sub_triangle,step)

    deallocate(xi)

  end subroutine base_Hsieh_Clough_Tocher_complete


  !*******************************************************************************
  !                    comptuting the normal derivative 
  !*******************************************************************************

  subroutine get_normal_der(deriv,i1,i2,mesh,freedom)
    type(hex_mesh_2d), pointer             :: mesh
    sll_real64, dimension(:,:), intent(in) :: deriv 
    sll_real64, dimension(3) , intent(out) :: freedom
    sll_real64                             :: x1, x2
    sll_int32, intent(in)                  :: i1, i2
    sll_real64                             :: fi_1, fi, fi1, fi2, fi_2, fi3  
    sll_int32                              :: ni_1, ni, ni1, ni2, ni_2, ni3  
    sll_int32                              :: h1, h2, num_cells, h16
    sll_int32                              :: h11,h12,h13,h21,h22,h23,h14,h24
    sll_int32                              :: h26, h2_4, h2_5,h25,h1_5,h1_4, h15
    sll_int32                              :: h1_1,h1_3,h2_1,h2_3, h1_2, h2_2
    sll_real64,dimension(2)                :: n1_l,n2_l,n3_l,n1_r,n2_r,n3_r
    
    ! à généraliser

    n1_l = (/-1._f64/sqrt(3._f64) , -2._f64/sqrt(3._f64)/)
    n2_l = (/-1._f64/sqrt(3._f64) , +1._f64/sqrt(3._f64)/)
    n3_l = (/+2._f64/sqrt(3._f64) , +1._f64/sqrt(3._f64)/)

    n1_r = (/-2._f64/sqrt(3._f64) , -1._f64/sqrt(3._f64)/)
    n2_r = (/+1._f64/sqrt(3._f64) , -1._f64/sqrt(3._f64)/)
    n3_r = (/+1._f64/sqrt(3._f64) , +2._f64/sqrt(3._f64)/)

    num_cells = mesh%num_cells   

    x1 = mesh%cartesian_coord(1,i1) 
    x2 = mesh%cartesian_coord(1,i2) 

    h1 = mesh%hex_coord(1,i1)
    h2 = mesh%hex_coord(2,i1)

    ! optimization variables
    h1_1 = h1 - 1
    h1_2 = h1 - 2
    h1_3 = h1 - 3
    h1_4 = h1 - 4
    h1_5 = h1 - 5
    h2_1 = h2 - 1
    h2_2 = h2 - 2
    h2_3 = h2 - 3
    h2_4 = h2 - 4
    h2_5 = h2 - 5
    h11  = h1 + 1
    h12  = h1 + 2
    h13  = h1 + 3
    h14  = h1 + 4
    h15  = h1 + 5
    h16  = h1 + 6 
    h21  = h2 + 1
    h22  = h2 + 2
    h23  = h2 + 3
    h24  = h2 + 4
    h25  = h2 + 5
    h26  = h2 + 6

    ! let us determine the configuration : is the triangle oriented left or right ?   

    if ( x1 > x2 ) then ! oriented left
       
       ! the first normal derivative is oriented normal to r1 and m1 is in [S2;S3]

       if ( test_in(h13,h26,num_cells) ) then
          fi_2 = 0._f64  
       else
          ni_2 = hex_to_global(mesh,h13,h26)
          fi_2 = deriv(1,ni_2) * n1_l(1) + deriv(2,ni_2) * n1_l(2)
       endif
       
       if ( test_in(h12,h24,num_cells) ) then
          fi_1 = 0._f64 
       else
          ni_1 = hex_to_global(mesh,h12,h24)
          fi_1 = deriv(1,ni_1) * n1_l(1) + deriv(2,ni_1)* n1_l(2)
       endif

       if ( test_in(h11,h22,num_cells) ) then
          fi = 0._f64  
       else
          ni = hex_to_global(mesh,h11,h22)
          fi = deriv(1,ni) * n1_l(1) + deriv(2,ni)* n1_l(2)
       endif

       if ( test_in(h1,h2,num_cells) ) then
          fi1 = 0._f64  
       else
          ni1 = hex_to_global(mesh,h1,h2)
          fi1 = deriv(1,ni1) * n1_l(1) + deriv(2,ni1)* n1_l(2)
       endif

       if ( test_in(h1_1,h2_2,num_cells) ) then
          fi2 = 0._f64  
       else
          ni2 = hex_to_global(mesh,h1_1,h2_2)
          fi2 = deriv(1,ni2) * n1_l(1) + deriv(2,ni2)* n1_l(2)
       endif

       if ( test_in(h1_2,h2_4,num_cells) ) then
          fi3 = 0._f64  
       else
          ni3 = hex_to_global(mesh,h1_2,h2_4)
          fi3 = deriv(1,ni3) * n1_l(1) + deriv(2,ni3) * n1_l(2)
       endif

       !freedom(1) = ( -fi_1 + 9._f64*(fi + fi1) - fi2 )/16._f64
       freedom(1) = ( 3._f64*(fi_2+fi3) - 25._f64*(fi_1+fi2) + &
            150._f64*(fi+fi1)) / 256._f64 


       ! the second normal derivative is oriented normal to r3 and m2 is in [S1;S3]

       if ( test_in(h13,h2_2,num_cells) ) then
          fi_2 = 0._f64  
       else
          ni_2 = hex_to_global(mesh,h13,h2_2)
          fi_2 = deriv(1,ni_2) * n2_l(1) + deriv(2,ni_2) * n2_l(2)
       endif       

       if ( test_in(h12,h2_1,num_cells) ) then
          fi_1 = 0._f64 
       else
          ni_1 = hex_to_global(mesh,h12,h2_1)
          fi_1 = deriv(1,ni_1) * n2_l(1) + deriv(2,ni_1)* n2_l(2)
       endif

       if ( test_in(h11,h2,num_cells) ) then
          fi = 0._f64  
       else
          ni = hex_to_global(mesh,h11,h2)
          fi = deriv(1,ni) * n2_l(1) + deriv(2,ni)* n2_l(2)
       endif

       if ( test_in(h1,h2,num_cells) ) then
          fi1 = 0._f64  
       else
          ni1 = hex_to_global(mesh,h1,h21)
          fi1 = deriv(1,ni1) * n2_l(1) + deriv(2,ni1)* n2_l(2)
       endif

       if ( test_in(h1_1,h22,num_cells) ) then
          fi2 = 0._f64  
       else
          ni2 = hex_to_global(mesh,h1_1,h22)
          fi2 = deriv(1,ni2) * n2_l(1) + deriv(2,ni2)* n2_l(2)
       endif
       
       if ( test_in(h1_2,h23,num_cells) ) then
          fi3 = 0._f64  
       else
          ni3 = hex_to_global(mesh,h1_2,h23)
          fi3 = deriv(1,ni3) * n2_l(1) + deriv(2,ni3) * n2_l(2)
       endif

       !freedom(2) = ( -fi_1 + 9._f64*(fi + fi1) - fi2 )/16._f64
       freedom(2) = ( 3._f64*(fi_2+fi3) - 25._f64*(fi_1+fi2) + &
            150._f64*(fi+fi1)) / 256._f64 

       ! the third normal derivative is oriented normal to r2 and m3 is in [S1;S2]

       if ( test_in(h1_5,h2_2,num_cells) ) then
          fi_2 = 0._f64  
       else
          ni_2 = hex_to_global(mesh,h1_5,h2_2)
          fi_2 = deriv(1,ni_2) * n3_l(1) + deriv(2,ni_2) * n3_l(2)
       endif       

       if ( test_in(h1_3,h2_1,num_cells) ) then
          fi_1 = 0._f64 
       else
          ni_1 = hex_to_global(mesh,h1_3,h2_1)
          fi_1 = deriv(1,ni_1) * n3_l(1) + deriv(2,ni_1)* n3_l(2)
       endif

       if ( test_in(h1_1,h2,num_cells) ) then
          fi = 0._f64  
       else
          ni = hex_to_global(mesh,h1_1,h2)
          fi = deriv(1,ni) * n3_l(1) + deriv(2,ni)* n3_l(2)
       endif

       if ( test_in(h11,h21,num_cells) ) then
          fi1 = 0._f64  
       else
          ni1 = hex_to_global(mesh,h11,h21)
          fi1 = deriv(1,ni1) * n3_l(1) + deriv(2,ni1)* n3_l(2)
       endif

       if ( test_in(h13,h22,num_cells) ) then
          fi2 = 0._f64  
       else
          ni2 = hex_to_global(mesh,h13,h22)
          fi2 = deriv(1,ni2) * n3_l(1) + deriv(2,ni2)* n3_l(2)
       endif

       if ( test_in(h15,h23,num_cells) ) then
          fi3 = 0._f64  
       else
          ni3 = hex_to_global(mesh,h15,h23)
          fi3 = deriv(1,ni3) * n3_l(1) + deriv(2,ni3) * n3_l(2)
       endif

       !freedom(3) = ( -fi_1 + 9._f64*(fi + fi1) - fi2 )/16._f64
       freedom(3) = ( 3._f64*(fi_2+fi3) - 25._f64*(fi_1+fi2) + &
            150._f64*(fi+fi1)) / 256._f64 

    else  ! oriented right
       
       ! the first normal derivative is oriented normal to r2 and m1 is in [S2;S3]
       
       if ( test_in(h16,h23,num_cells) ) then
          fi_2 = 0._f64  
       else
          ni_2 = hex_to_global(mesh,h16,h23)
          fi_2 = deriv(1,ni_2) * n1_r(1) + deriv(2,ni_2) * n1_r(2)
       endif
       
       if ( test_in(h14,h22,num_cells) ) then
          fi_1 = 0._f64 
       else
          ni_1 = hex_to_global(mesh,h14,h22)
          fi_1 = deriv(1,ni_1) * n1_r(1) + deriv(2,ni_1)* n1_r(2)
       endif

       if ( test_in(h12,h21,num_cells) ) then
          fi = 0._f64  
       else
          ni = hex_to_global(mesh,h12,h21)
          fi = deriv(1,ni) * n1_r(1) + deriv(2,ni)* n1_r(2)
       endif

       if ( test_in(h1,h2,num_cells) ) then
          fi1 = 0._f64  
       else
          ni1 = hex_to_global(mesh,h1,h2)
          fi1 = deriv(1,ni1) * n1_r(1) + deriv(2,ni1) * n1_r(2)
       endif

       if ( test_in(h1_2,h2_1,num_cells) ) then
          fi2 = 0._f64  
       else
          ni2 = hex_to_global(mesh,h1_2,h2_1)
          fi2 = deriv(1,ni2) * n1_r(1) + deriv(2,ni2) * n1_r(2)
       endif

       
       if ( test_in(h1_4,h2_2,num_cells) ) then
          fi3 = 0._f64  
       else
          ni3 = hex_to_global(mesh,h1_4,h2_2)
          fi3 = deriv(1,ni3) * n1_r(1) + deriv(2,ni3) * n1_r(2)
       endif

       !freedom(1) = ( -fi_1 + 9._f64*(fi + fi1) - fi2 )/16._f64
       freedom(1) = ( 3._f64*(fi_2+fi3) - 25._f64*(fi_1+fi2) + &
            150._f64*(fi+fi1)) / 256._f64 

       ! the second normal derivative is oriented normal to r3 and m2 is in [S1;S3]

       
       if ( test_in(h1_2,h23,num_cells) ) then
          fi_2 = 0._f64  
       else
          ni_2 = hex_to_global(mesh,h1_2,h23)
          fi_2 = deriv(1,ni_2) * n2_r(1) + deriv(2,ni_2) * n2_r(2)
       endif
       
       if ( test_in(h1_1,h22,num_cells) ) then
          fi_1 = 0._f64 
       else
          ni_1 = hex_to_global(mesh,h1_1,h22)
          fi_1 = deriv(1,ni_1) * n2_r(1) + deriv(2,ni_1) * n2_r(2)
       endif

       if ( test_in(h1,h21,num_cells) ) then
          fi = 0._f64  
       else
          ni = hex_to_global(mesh,h1,h21)
          fi = deriv(1,ni) * n2_r(1) + deriv(2,ni) * n2_r(2)
       endif

       if ( test_in(h11,h2,num_cells) ) then
          fi1 = 0._f64  
       else
          ni1 = hex_to_global(mesh,h11,h2)
          fi1 = deriv(1,ni1) * n2_r(1) + deriv(2,ni1) * n2_r(2)
       endif

       if ( test_in(h12,h2_1,num_cells) ) then
          fi2 = 0._f64  
       else
          ni2 = hex_to_global(mesh,h12,h2_1)
          fi2 = deriv(1,ni2) * n2_r(1) + deriv(2,ni2) * n2_r(2)
       endif

       
       if ( test_in(h13,h2_2,num_cells) ) then
          fi3 = 0._f64  
       else
          ni3 = hex_to_global(mesh,h13,h2_2)
          fi3 = deriv(1,ni3) * n2_r(1) + deriv(2,ni3) * n2_r(2)
       endif

       !freedom(2) = ( -fi_1 + 9._f64*(fi + fi1) - fi2 )/16._f64! interpolation
       freedom(2) = ( 3._f64*(fi_2+fi3) - 25._f64*(fi_1+fi2) + &
            150._f64*(fi+fi1)) / 256._f64 


       ! the third normal derivative is oriented normal to r2 and m3 is in [S1;S2]

       if ( test_in(h1_2,h2_5,num_cells) ) then
          fi_2 = 0._f64  
       else
          ni_2 = hex_to_global(mesh,h1_2,h2_5)
          fi_2 = deriv(1,ni_2) * n3_r(1) + deriv(2,ni_2) * n3_r(2)
       endif
       
       if ( test_in(h1_1,h2_3,num_cells) ) then
          fi_1 = 0._f64 
       else
          ni_1 = hex_to_global(mesh,h1_1,h2_3)
          fi_1 = deriv(1,ni_1) * n3_r(1) + deriv(2,ni_1) * n3_r(2)
       endif

       if ( test_in(h1,h2_1,num_cells) ) then
          fi = 0._f64  
       else
          ni = hex_to_global(mesh,h1,h2_1)
          fi = deriv(1,ni) * n3_r(1) + deriv(2,ni) * n3_r(2)
       endif

       if ( test_in(h11,h21,num_cells) ) then
          fi1 = 0._f64  
       else
          ni1 = hex_to_global(mesh,h11,h21)
          fi1 = deriv(1,ni1) * n3_r(1) + deriv(2,ni1) * n3_r(2)
       endif

       if ( test_in(h12,h23,num_cells) ) then
          fi2 = 0._f64  
       else
          ni2 = hex_to_global(mesh,h12,h23)
          fi2 = deriv(1,ni2) * n3_r(1) + deriv(2,ni2) * n3_r(2)
       endif

       if ( test_in(h13,h25,num_cells) ) then
          fi3 = 0._f64  
       else
          ni3 = hex_to_global(mesh,h13,h25)
          fi3 = deriv(1,ni3) * n3_r(1) + deriv(2,ni3) * n3_r(2)
       endif

       !freedom(3) = ( -fi_1 + 9._f64*(fi + fi1) - fi2 )/16._f64! interpolation
       freedom(3) = ( 3._f64*(fi_2+fi3) - 25._f64*(fi_1+fi2) + &
            150._f64*(fi+fi1)) / 256._f64 

    endif


  end subroutine get_normal_der
  
  function  test_in(h1,h2,num_cells) result(bool)  
    sll_int32  :: h1, h2,num_cells
    logical    :: bool

    bool =  abs(h1) > num_cells .or. abs(h2) > num_cells .or. &
         (h1)*(h2)< 0 .and. ( abs(h1) + abs(h2) > num_cells) 

  endfunction  test_in


  subroutine product_with_sigma_ganev_dimitrov(li,lj,lk,xi)
    sll_real64,dimension(:),intent(out) :: xi ! local base
    sll_real64,intent(in)               :: li, lj, lk ! barycentric coord
    ! optimisation variables 
    sll_real64                          :: li2, lj2, lk2
    sll_real64                          :: li3, lj3, lk3
    sll_real64                          :: li4, lj4, lk4
    sll_real64                          :: li3j, li3k, lj3i, lj3k, lk3i, lk3j
    sll_real64                          :: li2j2, li2k2, lj2k2
    sll_real64                          :: li2jk,lij2k,lijk2

    !( with i = 1, j = 2 et k = 3 )

    li2   = li**2
    lj2   = lj**2
    lk2   = lk**2
    li2j2 = li2*lj2 
    li2k2 = li2*lk2 
    lj2k2 = lj2*lk2
    li3   = li**3
    lj3   = lj**3
    lk3   = lk**3
    li4   = li2**2
    lj4   = lj2**2
    lk4   = lk2**2
    li3j  = li3*lj
    li3k  = li3*lk
    lj3i  = lj3*li
    lj3k  = lj3*lk 
    lk3i  = lk3*li
    lk3j  = lk3*lj
    li2jk = li2*lj*lk
    lij2k = li*lj2*lk
    lijk2 = li*lj*lk2

    xi(1)  = li4  + 4._f64*(li3k + li3j) - 5._f64*(li2k2 + li2j2) - 4._f64*li2jk
    xi(2)  = lj4  + 4._f64*(lj3i + lj3k) - 5._f64*(lj2k2 + li2j2) - 4._f64*lij2k
    xi(3)  = lk4  + 4._f64*(lk3j + lk3i) - 5._f64*(lj2k2 + li2k2) - 4._f64*lijk2
    xi(4)  = li3k - li2k2 + 0.5_f64 * ( - li2jk - lij2k + lijk2 )
    xi(5)  = li3j - li2j2 + 0.5_f64 * ( - li2jk + lij2k - lijk2 )
    xi(6)  = lj3i - li2j2 + 0.5_f64 * ( + li2jk - lij2k - lijk2 )
    xi(7)  = lj3k - lj2k2 + 0.5_f64 * ( - li2jk - lij2k + lijk2 )
    xi(8)  = lk3j - lj2k2 + 0.5_f64 * ( - li2jk + lij2k - lijk2 )
    xi(9)  = lk3i - li2k2 + 0.5_f64 * ( + li2jk - lij2k - lijk2 )
    xi(10) = 16._f64 * ( lj2k2 - li2jk + lij2k + lijk2 )
    xi(11) = 16._f64 * ( li2k2 + li2jk - lij2k + lijk2 ) 
    xi(12) = 16._f64 * ( li2j2 + li2jk + lij2k - lijk2 )
    xi(13) = 4._f64  * (-li2jk + lij2k + lijk2 ) 
    xi(14) = 4._f64  * ( li2jk - lij2k + lijk2 )
    xi(15) = 4._f64  * ( li2jk + lij2k - lijk2 )

  end subroutine product_with_sigma_ganev_dimitrov

  
  subroutine base_ganev_dimitrov(base,l1,l2,l3,step)
    sll_real64,dimension(:),intent(out) :: base
    sll_real64,intent(in)               :: l1,l2,l3,step
    sll_real64,dimension(:),allocatable :: xi 
    sll_real64                          :: step_sq  
    
    step_sq = step*sqrt(3._f64)*0.5_f64

    allocate(xi(1:15))

    call product_with_sigma_ganev_dimitrov(l1,l2,l3,xi)

     base(1) = xi(1) 
     base(2) = xi(2) 
     base(3) = xi(3) 
     base(4) = xi(5)*step
     base(5) = xi(4)*step
     base(6) = xi(6)*step
     base(7) = xi(7)*step
     base(8) = xi(9)*step
     base(9) = xi(8)*step
     base(10) = xi(10)
     base(11) = xi(11)
     base(12) = xi(12)
     base(13) = xi(13)*step_sq
     base(14) = xi(14)*step_sq
     base(15) = xi(15)*step_sq

    deallocate(xi)

  end subroutine base_ganev_dimitrov


end module interpolation_hex_hermite
