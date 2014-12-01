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


  num_degree = 15
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

  ! values at the middle of the edges of the triangle

  call f_interpolate((x1+x2)*0.5_f64,(y1+y2)*0.5_f64,freedom(12))
  call f_interpolate((x1+x3)*0.5_f64,(y1+y3)*0.5_f64,freedom(11))
  call f_interpolate((x3+x2)*0.5_f64,(y3+y2)*0.5_f64,freedom(10))


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

     ! exact values
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


     ! gradient suivant r1 et r2

     ! call f_grad_r1r2(pt_1_i_1(1),pt_1_i_1(2),dfi_1)
     ! dfni_1 = dfi_1(1)*n1_r_r(1) + dfi_1(2)*n1_r_r(2)
     ! call f_grad_r1r2(pt_1_i(1),pt_1_i(2),dfi)
     ! dfni   = dfi(1)*n1_r_r(1)   + dfi(2)*n1_r_r(2) 
     ! call f_grad_r1r2(pt_1_i1(1),pt_1_i1(2),dfi1)
     ! dfni1  = dfi1(1)*n1_r_r(1)  + dfi1(2)*n1_r_r(2)
     ! call f_grad_r1r2(pt_1_i2(1),pt_1_i2(2),dfi2)
     ! dfni2  = dfi2(1)*n1_r_r(1)  + dfi2(2)*n1_r_r(2)

     ! freedom(13) = ( -dfni_1 + 9._f64*dfni + 9._f64*dfni1 - dfni2 ) / 16._f64

     ! call f_grad_r1r2(pt_2_i_1(1),pt_2_i_1(2),dfi_1)
     ! dfni_1 = dfi_1(1)*n2_r_r(1) + dfi_1(2)*n2_r_r(2)
     ! call f_grad_r1r2(pt_2_i(1),pt_2_i(2),dfi)
     ! dfni   = dfi(1)*n2_r_r(1)   + dfi(2)*n2_r_r(2) 
     ! call f_grad_r1r2(pt_2_i1(1),pt_2_i1(2),dfi1)
     ! dfni1  = dfi1(1)*n2_r_r(1)  + dfi1(2)*n2_r_r(2)
     ! call f_grad_r1r2(pt_2_i2(1),pt_2_i2(2),dfi2)
     ! dfni2  = dfi2(1)*n2_r_r(1)  + dfi2(2)*n2_r_r(2)

     ! freedom(14) = ( -dfni_1 + 9._f64*dfni + 9._f64*dfni1 - dfni2 ) / 16._f64

     ! call f_grad_r1r2(pt_3_i_1(1),pt_3_i_1(2),dfi_1)
     ! dfni_1 = dfi_1(1)*n3_r_r(1) + dfi_1(2)*n3_r_r(2)
     ! call f_grad_r1r2(pt_3_i(1),pt_3_i(2),dfi)
     ! dfni   = dfi(1)*n3_r_r(1)   + dfi(2)*n3_r_r(2) 
     ! call f_grad_r1r2(pt_3_i1(1),pt_3_i1(2),dfi1)
     ! dfni1  = dfi1(1)*n3_r_r(1)  + dfi1(2)*n3_r_r(2)
     ! call f_grad_r1r2(pt_3_i2(1),pt_3_i2(2),dfi2)
     ! dfni2  = dfi2(1)*n3_r_r(1)  + dfi2(2)*n3_r_r(2)

     ! freedom(15) = ( -dfni_1 + 9._f64*dfni + 9._f64*dfni1 - dfni2 ) / 16._f64

     ! call  f_grad( (x2+x3)*0.5_f64,(y2+y3)*0.5_f64,df)
     ! print*, df(1)*n1_r(1) + df(2)*n1_r(2) , freedom(13)

     ! call  f_grad( (x1+x3)*0.5_f64,(y1+y3)*0.5_f64,df)
     ! print*, df(1)*n2_r(1) + df(2)*n2_r(2) ,  freedom(14)

     ! call  f_grad( (x2+x1)*0.5_f64,(y2+y1)*0.5_f64,df)
     ! print*, df(1)*n3_r(1) + df(2)*n3_r(2) ,  freedom(15)

     ! exact values
     call  f_grad( (x2+x3)*0.5_f64,(y2+y3)*0.5_f64,df)
     freedom(13) = df(1)*n1_r(1) + df(2)*n1_r(2) 
     call  f_grad( (x1+x3)*0.5_f64,(y1+y3)*0.5_f64,df)
     freedom(14) = df(1)*n2_r(1) + df(2)*n2_r(2) 
     call  f_grad( (x2+x1)*0.5_f64,(y2+y1)*0.5_f64,df)
     freedom(15) = df(1)*n3_r(1) + df(2)*n3_r(2) 


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


  call base_ganev_dimitrov(base,l1,l2,l3,step)

  f = 0._f64

  do i = 1,num_degree
     f = f + freedom(i)*base(i)
  enddo

  call f_interpolate(x,y,sol)

  print*, f, "valeur exacte : " , sol , "erreur : " , abs( f -sol )

  deallocate(freedom,base)


contains


subroutine f_interpolate(x,y,f)
  sll_real64,intent(in):: x,y
  sll_real64,intent(out)::f
  
  !f= x + 2._f64*y
  !f = x + x*y + x**2 + y**2 + y + 1
  !f = x**3 + y**3 + y**2 + x**2 + x + y + x**2*y + x*y**2 + x*y + 1._f64
  !f = x**3 + 2._f64*y**3
  f = x**5 + 2._f64*y**5
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

   ! df(1) =  3._f64*x**2 
   ! df(2) =  6._f64*y**2 

   df(1) =  5._f64*x**4 
   df(2) =  10._f64*y**4 

  ! df(1) = -2._f64 * (x - 1._f64) * exp( - ( (x - 1._f64)**2 + (y - 1._f64)**2 ) )
  ! df(2) = -2._f64 * (y - 1._f64) * exp( - ( (x - 1._f64)**2 + (y - 1._f64)**2 ) )  

endsubroutine f_grad


subroutine f_grad_r1r2(x,y,df)
  sll_real64,intent(in):: x,y
  sll_real64,dimension(2),intent(out)::df
  sll_real64,dimension(2) :: df1

  call f_grad(x,y,df1)

  df(1) =  sqrt(3._f64)*0.5_f64*df1(1) + df1(2)*0.5_f64
  df(2) = -sqrt(3._f64)*0.5_f64*df1(1) + df1(2)*0.5_f64
  

endsubroutine f_grad_r1r2




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

end program
