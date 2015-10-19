program hex_hermite_compute_fekete

! program to compute the values of the different basis functions of each hermite element at fekete points on a rectangular isocele triangle for the finite element solver

! the triangle is based on the points A(0,0) -> point 1, B(1,0) -> point 2 and C(0,1) -> pt 3 

sll_real64,dimension(:),allocatable::fekete_points

open(unit = 11, file="fekete_points", action="read", status="old")

number_fekete_points = 10

allocate( fekete_points(number_fekete_points,3) )

! get list of points with their barycentric coordinates
!->fekete_points
do i = 1,number_fekete_points
   read(11,*) j,fekete_points(i,1),fekete_points(i,2),fekete_points(i,3)
enddo

close(11)

! boucle sur ces points en prenant 1,0,0,0,0,... pour les degrés de liberté 

number_basis_function = 9  -> z9
!number_basis_function = 10 -> z10 
!number_basis_function = 9  -> hctr 
!number_basis_function = 12 -> hctc
!number_basis_function = 15 -> dimitrov


open(unit = 11, file="base_function_at_fekete_points", action="write", status="replace")

allocate( values_basis_function(number_fekete_points,number_basis_function) )

do i = 1,n
   l1 = fekete_points(i,1)
   l2 = fekete_points(i,2)
   l3 = fekete_points(i,3)
   call base_zienkiewicz_9_degree_of_freedom(base,l1,l2,l3) 
   !call base_zienkiewicz_10_degree_of_freedom(base,l1,l2,l3)  
   !call base_Hsieh_Clough_Tocher_reduced(base,x1,x3,y1,y2,y3,x,y,l1,l2,l3) 
   !call base_Hsieh_Clough_Tocher_complete(base,x1,x3,y1,y2,y3,x,y,l1,l2,l3)
   !call base_ganev_dimitrov(base,l1,l2,l3)
   do j=1,number_basis_function
      values_basis_function(i,j) = base(j)
      write(11,*) values_basis_function(i,j)
   enddo
enddo




contains

  subroutine barycentric_coord(x,y,l1,l2,l3)
    sll_real64, intent(in) :: x,y
    sll_real64, intent(out):: l1, l2, l3
    sll_real64             :: x1x,x2x,x3x,y1y,y2y,y3y, aire,a2

    aire = sqrt(2._f64)

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
  end subroutine

  subroutine base_zienkiewicz_9_degree_of_freedom(base,l1,l2,l3)
  sll_real64,dimension(:),intent(out) :: base
  sll_real64,intent(in)               :: l1, l2, l3 ! barycentric coord
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
  base(4) = ksi12 
  base(5) = ksi13
  base(6) = ksi21
  base(7) = ksi23 * sqrt(2._f64) ! distance BC
  base(8) = ksi31
  base(9) = ksi32 * sqrt(2._f64) ! distance CB

  end subroutine base_zienkiewicz_9_degree_of_freedom


  subroutine base_zienkiewicz_10_degree_of_freedom(base,l1,l2,l3)
  sll_real64,dimension(:),intent(out) :: base
  sll_real64,intent(in)               :: l1, l2, l3 ! barycentric coord
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
  base(4) = ( ksi12 - p15 )
  base(5) = ( ksi13 - p15 )
  base(6) = ( ksi21 - p15 )
  base(7) = sqrt(2._f64) * ( ksi23 - p15 )
  base(8) = ( ksi31 - p15 )
  base(9)=  sqrt(2._f64) * ( ksi32 - p15 )
  base(10) = 27._f64 * phi

  end subroutine base_zienkiewicz_10_degree_of_freedom


  
  subroutine get_sub_triangle_hex(x1,x3,y1,y2,y3,x,y,&
       num_sub_triangle)
    sll_real64,intent(in) :: x1, y1 ! cartesian coordinates of the lowest point
    sll_real64,intent(in) :: y2     ! ordinate of the middle point
    sll_real64,intent(in) :: x3, y3 ! cartesian coordinates of the highest point
    sll_real64,intent(in) :: x , y  ! cartesian coordinates of the point to be interpolated
    sll_int32,intent(out) :: num_sub_triangle
    sll_real64            :: slope



       !finding the sub triangle K_l in which the interpolated point is
    
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !                    TO BE CHANGED   TO BE CHANGED   TO BE CHANGED 
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************

    if ( y > y2 ) then! K1 or K2
       if ( y >= ( x - x3 ) * slope + y3 ) then 
          num_sub_triangle = 1
       else
          num_sub_triangle = 2
       endif
    else ! K2 or K3
       if ( y >= - ( x - x1 ) * slope + y1 ) then 
          num_sub_triangle = 2
       else
          num_sub_triangle = 3
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

    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !                    TO BE CHANGED   TO BE CHANGED   TO BE CHANGED 
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
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

    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !                    TO BE CHANGED   TO BE CHANGED   TO BE CHANGED 
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
  
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
  
  



  subroutine base_Hsieh_Clough_Tocher_reduced(base,x1,x3,y1,y2,y3,x,y,l1,l2,l3)
    sll_real64,dimension(:),intent(out) :: base
    sll_real64,intent(in)               :: x1,x3,y1,y2,y3
    sll_real64,intent(in)               :: x,y,l1,l2,l3
    sll_real64,dimension(:),allocatable :: xi 
    sll_real64                          :: li, lj, lk, step
    sll_int32                           :: num_sub_triangle

    !finding the sub triangle K_l in which the interpolated point is

    call get_sub_triangle_hex(x1,x3,y1,y2,y3,x,y,&
         num_sub_triangle)

    allocate(xi(1:9))
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !                    TO BE CHANGED   TO BE CHANGED   TO BE CHANGED 
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************

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

    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !                    TO BE CHANGED   TO BE CHANGED   TO BE CHANGED 
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************

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

  
  subroutine base_from_local_base_xi_htc_c(base,xi,num_sub_triangle)
    sll_real64,dimension(:),intent(out) :: base
    sll_real64,dimension(:),intent(in)  :: xi 
    sll_int32,              intent(in)  :: num_sub_triangle

    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !                    TO BE CHANGED   TO BE CHANGED   TO BE CHANGED 
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
    !*********************************************************************************
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


  subroutine base_Hsieh_Clough_Tocher_complete(base,x1,x3,y1,y2,y3,x,y,l1,l2,l3)
    sll_real64,dimension(:),intent(out) :: base
    sll_real64,intent(in)               :: x1,x3,y1,y2,y3
    sll_real64,intent(in)               :: x,y,l1,l2,l3
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
    call base_from_local_base_xi_htc_c(base,xi,num_sub_triangle)

    deallocate(xi)

  end subroutine base_Hsieh_Clough_Tocher_complete


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


  subroutine base_ganev_dimitrov(base,l1,l2,l3)
    sll_real64,dimension(:),intent(out) :: base
    sll_real64,intent(in)               :: l1,l2,l3
    sll_real64,dimension(:),allocatable :: xi 
    sll_real64                          :: step, step_sq  


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



end program hex_hermite_compute_fekete
