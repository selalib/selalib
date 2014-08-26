program test_hctc

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
  use sll_constants

! small test to verify each part of the computation of an interpolation using 
! the hctc hermite finite element

  implicit none
  sll_int32    :: i, num_degree 
  sll_real64   :: x1, x2, x3, y1, y2, y3, f,sol,step
  sll_real64   :: x1x,x2x,x3x,y1y,y2y,y3y, aire,a2
  sll_real64   :: l1, l2, l3
  sll_real64   :: x, y
  sll_real64   :: dfni_1, dfni, dfni1, dfni2 
  sll_real64,dimension(2) :: dfi_1, dfi, dfi1, dfi2 
  sll_real64,dimension(2) :: n1_r_r,n2_r_r,n3_r_r,n1_r_l,n2_r_l,n3_r_l
  sll_real64,dimension(2) :: df, h1,h2,h3,n1_l,n2_l,n3_l,n1_r,n2_r,n3_r
  sll_real64,dimension(:),allocatable :: freedom, base
  
  sll_real64,dimension(2)   :: pt_1_i_1, pt_1_i, pt_1_i1, pt_1_i2
  sll_real64,dimension(2)   :: pt_2_i_1, pt_2_i, pt_2_i1, pt_2_i2
  sll_real64,dimension(2)   :: pt_3_i_1, pt_3_i, pt_3_i1, pt_3_i2 

  pt_1_i_1 = (/sqrt(3._f64),3._f64/)
  pt_1_i   = (/0.5_f64*sqrt(3._f64),1.5_f64/)
  pt_1_i1  = (/+0._f64,0._f64/)
  pt_1_i2  = (/-0.5_f64*sqrt(3._f64),-1.5_f64/)

  pt_2_i_1 = (/-1.5_f64*sqrt(3._f64),0.5_f64/)
  pt_2_i   = (/-0.5_f64*sqrt(3._f64),0.5_f64/)
  pt_2_i1  = (/+0.5_f64*sqrt(3._f64),0.5_f64/)
  pt_2_i2  = (/+1.5_f64*sqrt(3._f64),0.5_f64/)

  pt_3_i_1 = (/sqrt(3._f64),-2._f64/)
  pt_3_i   = (/+0.5_f64*sqrt(3._f64),-0.5_f64/)
  pt_3_i1  = (/+0._f64,1._f64/)
  pt_3_i2  = (/-0.5_f64*sqrt(3._f64),2.5_f64/)



  h1 = (/+sqrt(3._f64)*0.5_f64,+0.5_f64/)
  h3 = (/0._f64,+1._f64/)
  h2 = (/-sqrt(3._f64)*0.5_f64,+0.5_f64/)

  n1_l = (/+0.5_f64 , -sqrt(3._f64)*0.5_f64/)
  n2_l = (/-1._f64, 0._f64/)
  n3_l = (/+0.5_f64 , +sqrt(3._f64)*0.5_f64/)

  n1_r = (/-0.5_f64 , -sqrt(3._f64)*0.5_f64/)
  n2_r = (/ 1._f64, 0._f64/)
  n3_r = (/-0.5_f64 , +sqrt(3._f64)*0.5_f64/)

  
  n1_r_l = (/-1._f64/sqrt(3._f64) , -2._f64/sqrt(3._f64)/)
  n2_r_l = (/-1._f64/sqrt(3._f64) , +1._f64/sqrt(3._f64)/)
  n3_r_l = (/+2._f64/sqrt(3._f64) , +1._f64/sqrt(3._f64)/)

  n1_r_r = (/-2._f64/sqrt(3._f64) , -1._f64/sqrt(3._f64)/)
  n2_r_r = (/+1._f64/sqrt(3._f64) , -1._f64/sqrt(3._f64)/)
  n3_r_r = (/+1._f64/sqrt(3._f64) , +2._f64/sqrt(3._f64)/)


  num_degree = 12
  step = 1._f64
  aire = step**2*sqrt(3._f64)*0.25_f64
  x = 0.4_f64
  y = 0.6_f64

  allocate(freedom(num_degree),base(num_degree ))

  x1 = 0._f64
  x2 = sqrt(3._f64)*0.5_f64
  x3 = 0._f64
  y1 = 0._f64
  y2 = 0.5_f64
  y3 = 1._f64

  
  ! values at the vertexes of the triangle
  
  call f_interpolate(x1,y1,freedom(1))
  call f_interpolate(x2,y2,freedom(2))
  call f_interpolate(x3,y3,freedom(3))

  ! values of the derivatives
  call  f_grad(x1,y1,df)
  freedom(5) = df(1) * h3(1) + df(2) * h3(2) ! derivative from S1 to S3 (+h3)
  call  f_grad(x3,y3,df)
  freedom(8) = - df(1) * h3(1) - df(2) * h3(2) ! derivative from S3 to S1 (-h3)

  if ( x2 <= x1 ) then ! first case : triangle oriented left

     call  f_grad(x1,y1,df)
     freedom(4) = df(1) * h1(1) + df(2) * h1(2)! derivative from S1 to S2 (+h2)
     call  f_grad(x2,y2,df)
     freedom(6) = -df(1) * h1(1) - df(2) * h1(2)! derivative from S2 to S1 (-h2)
     call  f_grad(x2,y2,df)
     freedom(7) = df(1) * h2(1) + df(2) * h2(2)! derivative from S2 to S3 (+h1)
     call  f_grad(x3,y3,df)
     freedom(9) = -df(1) * h2(1) - df(2) * h2(2)! derivative from S3 to S2 (-h1)


     ! call  f_grad( (x2+x3)*0.5_f64,(y2+y3)*0.5_f64,df)
     ! freedom(10) = df(1)*n1_l(1) + df(2)*n1_l(2) 
     ! call  f_grad( (x1+x3)*0.5_f64,(y1+y3)*0.5_f64,df)
     ! freedom(11) = df(1)*n2_l(1) + df(2)*n2_l(2) 
     ! call  f_grad( (x2+x1)*0.5_f64,(y2+y1)*0.5_f64,df)
     ! freedom(12) = df(1)*n3_l(1) + df(2)*n3_l(2) 

  else if ( x2 > x1 ) then ! second case triangle oriented right 

     call  f_grad(x1,y1,df)
     freedom(4) = df(1) * h1(1) + df(2) * h1(2)! derivative from S1 to S2 (+h1)
     call  f_grad(x2,y2,df)
     freedom(6) = -df(1) * h1(1) - df(2) * h1(2)! derivative from S2 to S1 (-h1)
     call  f_grad(x2,y2,df)
     freedom(7) = df(1) * h2(1) + df(2) * h2(2)! derivative from S2 to S3 (+h2)
     call  f_grad(x3,y3,df)
     freedom(9) = -df(1) * h2(1) - df(2) * h2(2)! derivative from S3 to S2 (-h2)


     ! call f_grad(pt_1_i_1(1),pt_1_i_1(2),dfi_1)
     ! dfni_1 = dfi_1(1)*n1_r(1) + dfi_1(2)*n1_r(2)
     ! call f_grad(pt_1_i(1),pt_1_i(2),dfi)
     ! dfni   = dfi(1)*n1_r(1)   + dfi(2)*n1_r(2) 
     ! call f_grad(pt_1_i1(1),pt_1_i1(2),dfi1)
     ! dfni1  = dfi1(1)*n1_r(1)  + dfi1(2)*n1_r(2)
     ! call f_grad(pt_1_i2(1),pt_1_i2(2),dfi2)
     ! dfni2  = dfi2(1)*n1_r(1)  + dfi2(2)*n1_r(2)

     ! freedom(10) = ( -dfni_1 + 9._f64*dfni + 9._f64*dfni1 - dfni2 ) / 16._f64

     ! call f_grad(pt_2_i_1(1),pt_2_i_1(2),dfi_1)
     ! dfni_1 = dfi_1(1)*n2_r(1) + dfi_1(2)*n2_r(2)
     ! call f_grad(pt_2_i(1),pt_2_i(2),dfi)
     ! dfni   = dfi(1)*n2_r(1)   + dfi(2)*n2_r(2) 
     ! call f_grad(pt_2_i1(1),pt_2_i1(2),dfi1)
     ! dfni1  = dfi1(1)*n2_r(1)  + dfi1(2)*n2_r(2)
     ! call f_grad(pt_2_i2(1),pt_2_i2(2),dfi2)
     ! dfni2  = dfi2(1)*n2_r(1)  + dfi2(2)*n2_r(2)

     ! freedom(11) = ( -dfni_1 + 9._f64*dfni + 9._f64*dfni1 - dfni2 ) / 16._f64

     ! call f_grad(pt_3_i_1(1),pt_3_i_1(2),dfi_1)
     ! dfni_1 = dfi_1(1)*n3_r(1) + dfi_1(2)*n3_r(2)
     ! call f_grad(pt_3_i(1),pt_3_i(2),dfi)
     ! dfni   = dfi(1)*n3_r(1)   + dfi(2)*n3_r(2) 
     ! call f_grad(pt_3_i1(1),pt_3_i1(2),dfi1)
     ! dfni1  = dfi1(1)*n3_r(1)  + dfi1(2)*n3_r(2)
     ! call f_grad(pt_3_i2(1),pt_3_i2(2),dfi2)
     ! dfni2  = dfi2(1)*n3_r(1)  + dfi2(2)*n3_r(2)

     ! freedom(12) = ( -dfni_1 + 9._f64*dfni + 9._f64*dfni1 - dfni2 ) / 16._f64

     ! call  f_grad( (x2+x3)*0.5_f64,(y2+y3)*0.5_f64,df)
     ! print*, df(1)*n1_r(1) + df(2)*n1_r(2) , freedom(10)

     ! call  f_grad( (x1+x3)*0.5_f64,(y1+y3)*0.5_f64,df)
     ! print*, df(1)*n2_r(1) + df(2)*n2_r(2) ,  freedom(11)

     ! call  f_grad( (x2+x1)*0.5_f64,(y2+y1)*0.5_f64,df)
     ! print*, df(1)*n3_r(1) + df(2)*n3_r(2) ,  freedom(12)


     ! gradient suivant r1 et r2

     call f_grad_r1r2(pt_1_i_1(1),pt_1_i_1(2),dfi_1)
     dfni_1 = dfi_1(1)*n1_r_r(1) + dfi_1(2)*n1_r_r(2)
     call f_grad_r1r2(pt_1_i(1),pt_1_i(2),dfi)
     dfni   = dfi(1)*n1_r_r(1)   + dfi(2)*n1_r_r(2) 
     call f_grad_r1r2(pt_1_i1(1),pt_1_i1(2),dfi1)
     dfni1  = dfi1(1)*n1_r_r(1)  + dfi1(2)*n1_r_r(2)
     call f_grad_r1r2(pt_1_i2(1),pt_1_i2(2),dfi2)
     dfni2  = dfi2(1)*n1_r_r(1)  + dfi2(2)*n1_r_r(2)

     freedom(10) = ( -dfni_1 + 9._f64*dfni + 9._f64*dfni1 - dfni2 ) / 16._f64

     call f_grad_r1r2(pt_2_i_1(1),pt_2_i_1(2),dfi_1)
     dfni_1 = dfi_1(1)*n2_r_r(1) + dfi_1(2)*n2_r_r(2)
     call f_grad_r1r2(pt_2_i(1),pt_2_i(2),dfi)
     dfni   = dfi(1)*n2_r_r(1)   + dfi(2)*n2_r_r(2) 
     call f_grad_r1r2(pt_2_i1(1),pt_2_i1(2),dfi1)
     dfni1  = dfi1(1)*n2_r_r(1)  + dfi1(2)*n2_r_r(2)
     call f_grad_r1r2(pt_2_i2(1),pt_2_i2(2),dfi2)
     dfni2  = dfi2(1)*n2_r_r(1)  + dfi2(2)*n2_r_r(2)

     freedom(11) = ( -dfni_1 + 9._f64*dfni + 9._f64*dfni1 - dfni2 ) / 16._f64

     call f_grad_r1r2(pt_3_i_1(1),pt_3_i_1(2),dfi_1)
     dfni_1 = dfi_1(1)*n3_r_r(1) + dfi_1(2)*n3_r_r(2)
     call f_grad_r1r2(pt_3_i(1),pt_3_i(2),dfi)
     dfni   = dfi(1)*n3_r_r(1)   + dfi(2)*n3_r_r(2) 
     call f_grad_r1r2(pt_3_i1(1),pt_3_i1(2),dfi1)
     dfni1  = dfi1(1)*n3_r_r(1)  + dfi1(2)*n3_r_r(2)
     call f_grad_r1r2(pt_3_i2(1),pt_3_i2(2),dfi2)
     dfni2  = dfi2(1)*n3_r_r(1)  + dfi2(2)*n3_r_r(2)

     freedom(12) = ( -dfni_1 + 9._f64*dfni + 9._f64*dfni1 - dfni2 ) / 16._f64

     call  f_grad( (x2+x3)*0.5_f64,(y2+y3)*0.5_f64,df)
     print*, df(1)*n1_r(1) + df(2)*n1_r(2) , freedom(10)

     call  f_grad( (x1+x3)*0.5_f64,(y1+y3)*0.5_f64,df)
     print*, df(1)*n2_r(1) + df(2)*n2_r(2) ,  freedom(11)

     call  f_grad( (x2+x1)*0.5_f64,(y2+y1)*0.5_f64,df)
     print*, df(1)*n3_r(1) + df(2)*n3_r(2) ,  freedom(12)

     ! valeurs exactes

     ! call  f_grad( (x2+x3)*0.5_f64,(y2+y3)*0.5_f64,df)
     ! freedom(10) = df(1)*n1_r(1) + df(2)*n1_r(2) 
     ! call  f_grad( (x1+x3)*0.5_f64,(y1+y3)*0.5_f64,df)
     ! freedom(11) = df(1)*n2_r(1) + df(2)*n2_r(2) 
     ! call  f_grad( (x2+x1)*0.5_f64,(y2+y1)*0.5_f64,df)
     ! freedom(12) = df(1)*n3_r(1) + df(2)*n3_r(2) 


  endif

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


  call base_Hsieh_Clough_Tocher_complete&
       (base,x1,x2,x3,y1,y2,y3,x,y,l1,l2,l3,step)

  f = 0._f64

  do i = 1,num_degree
     f = f + freedom(i)*base(i)
  enddo

  call f_interpolate(x,y,sol)

  print*, f, "valeur exacte : " , sol , "erreur : " , abs( f -sol )

  deallocate(freedom,base)




contains


  subroutine get_sub_triangle_hex(x1,x2,x3,y1,y2,y3,x,y,&
       num_sub_triangle)
    sll_real64,intent(in) :: x1, y1 ! cartesian coordinates of the lowest point
    sll_real64,intent(in) :: x2, y2 
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



  !***************************************************************************
  !                           subroutines for hctc
  !***************************************************************************





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
    xi(2) =          0.5_f64*li3 + lj3 - 1.5_f64*li2k &
         + 3.0_f64*(lj2i + lj2k + lijk )
    xi(3) =          0.5_f64*li3 + lk3              - 1.5_f64*li2j &
         + 3.0_f64*(lk2j + lk2i + lijk )
    xi(4) = - 1._f64/12._f64*li3 + 1.75_f64 * li2k -  0.5_f64*li2j
    xi(5) = - 1._f64/12._f64*li3 -  0.5_f64 * li2k + 1.75_f64*li2j
    xi(6) = - 7._f64/12._f64*li3 +  0.5_f64 * li2k + 1.25_f64*li2j + lj2i - lijk
    xi(7) =   2._f64 /3._f64*li3 - 0.75_f64 * li2k - 1.25_f64*li2j + lj2k + 1.5_f64*lijk 
    xi(8) =   2._f64 /3._f64*li3 - 1.25_f64 * li2k - 0.75_f64*li2j + lk2j + 1.5_f64*lijk
    xi(9) = - 7._f64/12._f64*li3 + 1.25_f64* li2k  + 0.5_f64 *li2j + lk2i - lijk
    xi(10) =   4._f64/3._f64*li3    - 2._f64*li2k -  2._f64 * li2j + 4._f64*lijk
    xi(11) = - 2._f64/3._f64*li3    + 2._f64*li2k
    xi(12) = - 2._f64/3._f64*li3                  +  2._f64 * li2j

    ! reduced

    ! xi(1) = 4.5_f64 * ( li2k + li2j ) 
    ! xi(2) = 0.5_f64*li3 + lj3 - 1.5_f64*li2k + 3.0_f64*(lj2i + lj2k + lijk )
    ! xi(3) = 0.5_f64*li3 + lk3 - 1.5_f64*li2j + 3.0_f64*(lk2j + lk2i + lijk )
    ! xi(4) = - 0.25_f64*li3 + 1.25_f64 * li2k + 0.5_f64*li2j
    ! xi(5) = - 0.25_f64*li3 + 0.5_f64 * li2k + 1.25_f64*li2j
    ! xi(6) =   0.25_f64*li3 - 0.5_f64 * li2k - 0.25_f64*li2j + lj2i + lijk 
    ! xi(7) = - 0.25_f64*li2k + 0.25_f64 * li2j + lj2k + 0.5_f64*lijk 
    ! xi(8) =   0.25_f64*li2k - 0.25_f64 * li2j + lk2j + 0.5_f64*lijk
    ! xi(9) =   0.25_f64*li3 - 0.25_f64*li2k - 0.5_f64*li2j + lk2i + lijk 

    ! xi(10:12) = 0.0_f64


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


  subroutine base_Hsieh_Clough_Tocher_complete(base,x1,x2,x3,y1,y2,y3,x,y,l1,l2,l3,step)
    sll_real64,dimension(:),intent(out) :: base
    sll_real64,intent(in)               :: x1,x2,x3,y1,y2,y3
    sll_real64,intent(in)               :: x,y,l1,l2,l3,step
    sll_real64,dimension(:),allocatable :: xi 
    sll_real64                          :: li, lj, lk
    sll_int32                           :: num_sub_triangle

    !finding the sub triangle K_l in which the interpolated point is

    call get_sub_triangle_hex(x1,x2,x3,y1,y2,y3,x,y,&
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

subroutine f_interpolate(x,y,f)
  sll_real64,intent(in):: x,y
  sll_real64,intent(out)::f
  
  !f= x + 2._f64*y
  !f = x + x*y + x**2 + y**2 + y + 1
  !f = x**3 + y**3 + y**2 + x**2 + x + y + x**2*y + x*y**2 + x*y + 1._f64
  f = x**3 + 2._f64*y**3
  !f = exp( - ( (x - 1._f64)**2 + (y - 1._f64)**2 ) )  

endsubroutine f_interpolate


subroutine f_grad(x,y,df)
  sll_real64,intent(in):: x,y
  sll_real64,dimension(2),intent(out)::df
  
  ! df(1) = 1._f64
  ! df(2) = 2._f64

   !df(1) = 1._f64 + y + 2._f64 * x
   !df(2) = 1._f64 + x + 2._f64 * y

   ! df(1) = 3._f64*x**2 + 2._f64*x + 1._f64 + 2._f64*x*y + y**2 + y
   ! df(2) = 3._f64*y**2 + 2._f64*y + 1._f64 + 2._f64*x*y + x**2 + x

   df(1) =  3._f64*x**2 
   df(2) =  6._f64*y**2 

  ! df(1) = -2._f64 * (x - 1._f64) * exp( - ( (x - 1._f64)**2 + (y - 1._f64)**2 ) )
  ! df(2) = -2._f64 * (y - 1._f64) * exp( - ( (x - 1._f64)**2 + (y - 1._f64)**2 ) )  

endsubroutine f_grad



! subroutine f_interpolate_r1r2(x,y,f)
!   sll_real64,intent(in):: x,y
!   sll_real64,intent(out)::f

!  call  f_interpolate(sqrt(3._f64)*0.5_f64*(x-y),(x+y)*0.5_f64,f)

! endsubroutine f_interpolate_r1r2


subroutine f_grad_r1r2(x,y,df)
  sll_real64,intent(in):: x,y
  sll_real64,dimension(2),intent(out)::df
  sll_real64,dimension(2) :: df1

  call f_grad(x,y,df1)

  df(1) =  sqrt(3._f64)*0.5_f64*df1(1) + df1(2)*0.5_f64
  df(2) = -sqrt(3._f64)*0.5_f64*df1(1) + df1(2)*0.5_f64
  

endsubroutine f_grad_r1r2



end program
