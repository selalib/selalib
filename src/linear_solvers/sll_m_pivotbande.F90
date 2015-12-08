module sll_m_pivotbande
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    factolub_bande, &
    residue_bande, &
    solvlub_bande

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains

  subroutine searchband(a,n,l1,l2)
    implicit none
    sll_real64,dimension(:,:)::a
    sll_int32,intent(in )    ::n
    sll_int32,intent(out)    ::l1,l2
    sll_int32,dimension(n)   ::b1,b2
    sll_int32::i,j,ll1,ll2!,width

    l1=0
    l2=0
    ! let us find the width of the band of the matrix
    ! to the right (l2) and to the left (1)

    do j=1,n
       do i=j,n
          if (abs(a(i,j)) > 1e-16) then
             ll1 = i - j
             ll2 = i - j
          endif
       enddo
       b1(j) = ll1
       b2(j) = ll2
    enddo

    l1 = b1(1)
    l2 = b2(1)

    do i = 2,n
       if (b1(i) > l1) l1 = b1(i)
       if (b2(i) > l2) l2 = b2(i)
    enddo

    !width = l1 + l2 + 1 ! total width

  endsubroutine searchband






!***********************************************

  subroutine factolub (a,l,u,n,l1,l2)
    implicit none
    sll_real64,dimension(:,:),intent(in)::a
    sll_real64,dimension(:,:),intent(out)::l,u
    sll_int32,intent(in)::n,l1,l2
    sll_int32::i,j,k,ll1,ll2

    !init 

    u = a
    l = 0._f64
    do i = 1,n
       l(i,i) = 1._f64
    enddo

    !transformation
    
    do k = 1,n-1
       ll1 = k + l1
       if (ll1>n) ll1=n 
       do i = k + 1, ll1
          l(i,k) = u(i,k)/u(k,k)   
          ll2 = k+l2
          if (ll2>n) ll2 = n 
          do j = k,ll2
             u(i,j) = u(i,j)-l(i,k)*u(k,j)
          enddo
       enddo
    enddo

    ! do k = 1,n-1
    !    do i = k + 1, n
    !       l(i,k) = u(i,k)/u(k,k)   
    !       print*,i,k,u(k,k),u(i,k)  ,l(i,k) 
    !       do j = k,n
    !          u(i,j) = u(i,j)-l(i,k)*u(k,j)
    !       enddo
    !    enddo
    ! enddo


  end subroutine factolub
  !***********************************************

  subroutine factolub_bande(a,l,u,n,l1,l2)
    implicit none
    sll_real64,dimension(:,:),intent(in)::a
    sll_real64,dimension(:,:),intent(out)::l,u
    sll_int32,intent(in)::n,l1,l2
    sll_int32::i,j,k,ll1,band,ki,ji,ll1i

    band=l1+l2+1

    u = a
    l = 0._f64
    ll1 = l1 + 1
    do i = 1,n
       l(i,ll1) = 1._f64
    enddo


    do k = 1,n-1
       do i = 1,l1
          ki = k + i
          ll1i = ll1 - i          
          if (ki <= n) then 
             if (  abs(u(ki,ll1i)) > 1e-16 ) then
                l(ki,ll1i) = u(ki,ll1i)/u(k,ll1)
                !if (  abs(l(ki,ll1i)) > 1e-17 ) then
                do j = 1,band
                   ji  = j+i
                   if (ji<= band ) u(ki,j) = u(ki,j) - l(ki,ll1i)*u(k,ji)
                enddo
             endif
          endif
       enddo
    enddo

  end subroutine factolub_bande

  !***********************************************

  subroutine solvlub(l,u,x,b,n,l1,l2)
    implicit none
    sll_real64,dimension(:,:), intent(in)::l,u
    sll_real64,dimension(:)  , intent(in)::b
    sll_int32,                 intent(in)::l1,l2
    sll_real64,dimension(:),intent(inout)::x
    sll_int32,                 intent(in)::n
    sll_int32                            ::i,j,ll1,ll2
    sll_real64                           ::s

    x(1) = b(1)/l(1,1)

    do i = 2,n
       s = 0._f64
       ll1 = i - l1  
       if (ll1<1) ll1 = 1
       do j = ll1,i-1 
          s = s + l(i,j)*x(j)
       enddo
       x(i) = (b(i)-s)/l(i,i)
    enddo

    
    x(n) = x(n)/u(n,n)

    do i = n-1,1,-1
       s = 0._f64
       ll2 = i + l2
       if (ll2>n) ll2 = n
       do j = i+1 , ll2 
          s = s + u(i,j)*x(j)
       enddo
       x(i) = (x(i)-s)/u(i,i)
    enddo
    

  end subroutine solvlub

  !***********************************************

  subroutine solvlub_bande(l,u,x,b,n,l1,l2)
    implicit none
    sll_real64,dimension(:,:),   intent(in)::l,u
    sll_real64,dimension(:)  ,   intent(in)::b
    sll_int32,                   intent(in)::l1,l2
    sll_real64,dimension(:) , intent(inout)::x
    sll_int32,                   intent(in)::n
    sll_int32                              ::i,j,ll1,lll1,ll1j,ij
    sll_real64                             ::s

    lll1 = l1 + 1
    x(1) = b(1)/l(1,lll1)

    do i = 2,n
       s = 0._f64
       ll1 = i-lll1 
       do j = 1,l1 
          ll1j = ll1+j
          if (ll1j > 0) s = s + l(i,j)*x(ll1j)
       enddo
       x(i) = (b(i)-s)/l(i,lll1)
    enddo

    x(n) = x(n)/u(n,lll1)

    do i = n-1,1,-1
       s = 0._f64
       do j = 1,l2
          ij= i+j 
          if (ij <= n) s = s + u(i,lll1+j)*x(ij)
       enddo
       x(i) = (x(i)-s)/u(i,lll1)
    enddo



  end subroutine solvlub_bande


  !***********************************************


  subroutine gauss_seidel(a,x,b,n)
    sll_real64,dimension(:,:), intent(in) :: a
    sll_real64,dimension(:)  , intent(out):: x
    sll_real64,dimension(:)  , intent(in ):: b
    sll_int32,                 intent(in) :: n
    sll_real64,dimension(:),allocatable   :: x_k, x_k1
    sll_real64                            :: error_l2, s1, s2
    sll_int32                             :: i, j

    allocate(x_k(n),x_k1(n))

    x_k  = 1._f64
    x_k1 = x_k 

    error_l2 = 1._f64

    do while ( error_l2 > 1e-12 ) 

       do i = 1,n
          s1 = 0._f64
          s2 = 0._f64
          if (i-1>0) then
             do j = 1,i-1
                s1 = s1 + a(i,j)*x_k1(j)
             enddo
          endif
          do j = i+1,n
             s2 = s2 + a(i,j)*x_k(j)
          enddo

          x_k1(i) = ( b(i) - s1 - s2 ) / a(i,i)

       enddo

       call residue(a,x_k1,b,n,error_l2)

       x_k = x_k1
    enddo

    x = x_k

    deallocate(x_k, x_k1)

  end subroutine gauss_seidel

  subroutine residue(a,x,b,n,error_l2)
    sll_real64,dimension(:,:), intent(in)  :: a
    sll_real64,dimension(:)  , intent(in)  :: x
    sll_real64,dimension(:)  , intent(in ) :: b
    sll_int32,                 intent(in)  :: n
    sll_real64,                intent(out) :: error_l2
    sll_real64                             :: s
    sll_int32                              :: i, j

    ! computing ||Ax - b ||

    error_l2 = 0._f64

    do i = 1,n
       s = 0._f64
       do j = 1,n
          s = s + a(i,j)*x(j)
       enddo
       error_l2 = error_l2 + (s-b(i))**2
    enddo

    error_l2 = sqrt( error_l2 ) 

  end subroutine residue

  subroutine gauss_seidel_bande(a,x,b,l1,l2,n)
    sll_real64,dimension(:,:), intent(in) :: a
    sll_real64,dimension(:)  , intent(out):: x
    sll_real64,dimension(:)  , intent(in ):: b
    sll_int32,                 intent(in) :: n
    sll_real64,dimension(:),allocatable   :: x_k, x_k1
    sll_real64                            :: error_l2, s1, s2
    sll_int32,                 intent(in) :: l1, l2
    sll_int32                             :: i, j,k,l_tot

    allocate(x_k(n),x_k1(n))

    x_k  = 1._f64
    x_k1 = x_k 

    error_l2 = 1._f64
    l_tot = l1+l2+1

    do while ( error_l2 > 1e-12 ) 

       do i = 1,n

          s1 = 0._f64
          s2 = 0._f64
          k = i-l1-1

          do j = 1,l1
             if (k+j > 0 ) s1 = s1 + a(i,j)*x_k1(k+j)
          enddo

          do j = l1+2,l_tot
             if (k+j <= n) s2 = s2 + a(i,j)*x_k(k+j)
          enddo

          x_k1(i) = ( b(i) - s1 - s2 ) / a(i,1 + l1)

       enddo

       call residue_bande(a,x_k1,b,l1,l2,n,error_l2)
       x_k = x_k1

    enddo

    x = x_k

    deallocate(x_k, x_k1)

  end subroutine gauss_seidel_bande

  subroutine residue_bande(a,x,b,l1,l2,n,error_l2)
    sll_real64,dimension(:,:), intent(in)  :: a
    sll_real64,dimension(:)  , intent(in)  :: x
    sll_real64,dimension(:)  , intent(in ) :: b
    sll_int32,                 intent(in)  :: n
    sll_int32,                 intent(in)  :: l1, l2
    sll_real64,                intent(out) :: error_l2
    sll_real64                             :: s
    sll_int32                              :: i, j,k,l_tot,kj

    ! computing ||Ax - b ||

    error_l2 = 0._f64

    l_tot = l1+l2+1

    do i = 1,n
       s = 0._f64
       k = i-l1-1
       do j = 1,l_tot
          kj = k+j
          if (kj > 0 .and. kj <= n) s = s + a(i,j)*x(kj)
       enddo
       error_l2 = error_l2 + (s-b(i))**2
    enddo

    error_l2 = sqrt( error_l2 ) 

  end subroutine residue_bande

end module sll_m_pivotbande
