module pivotbande
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
  implicit none
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

    ! !init 

    ! u = a
    ! l = 0._f64
    ! do i = 1,n
    !    l(i,i) = 1._f64
    ! enddo

    ! !transformation

    ! do k = 1,n-1
    !    ll1 = k + l1
    !    if (ll1>n) ll1=n 
    !    do i = k + 1, ll1
    !       l(i,k) = u(i,k)/u(k,k)             
    !       ll2 = k+l2
    !       if (ll2>n) ll2 = n 
    !       do j = k,ll2
    !          u(i,j) = u(i,j)-l(i,k)*u(k,j)
    !       enddo
    !    enddo
    ! enddo

    ! version bande

    u = a
    l = 0._f64
    do i = 1,n
       l(i,l1+1) = 1._f64
    enddo

    ll1 = l1 + 1
    ll2 = l1 + 2

    do k = 1,n-1
       do i = 1,l1
          if (k+ll1-i <= n) l(k,i) = u(k+ll1-i,i)/u(k,ll1)
          do j = ll1,ll1+l2
             if (i+k <= n)  u(i+k,j-i) = u(i+k,j-i)-l(k,ll1-i)*u(k,j)
          enddo
       enddo
    enddo

  end subroutine factolub

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

    ! x(1) = b(1)/l(1,1)

    ! do i = 2,n
    !    s = 0._f64
    !    ll1 = i-l1  
    !    if (ll1<1) ll1 = 1
    !    do j = ll1,i-1 
    !       s = s + l(i,j)*x(j)
    !    enddo
    !    x(i) = (b(i)-s)/l(i,i)
    ! enddo

    ! x(n) = x(n)/u(n,n)
    ! do i = n-1,1,-1
    !    s = 0._f64
    !    ll2 = i+l2
    !    if (ll2>n) ll2 = n
    !    do j = (i+1),ll2 
    !       s = s + u(i,j)*x(j)
    !    enddo
    !    x(i) = (x(i)-s)/u(i,i)
    ! enddo

    !version bande
    x(1) = b(1)/l(1,1)

    do i = 2,n
       s = 0._f64
       ll1 = i-l1  
       if (ll1<1) ll1 = 1
       do j = ll1,i-1 
          s = s + l(i,j)*x(j)
       enddo
       x(i) = (b(i)-s)/l(i,i)
    enddo

    x(n) = x(n)/u(n,n)
    do i = n-1,1,-1
       s = 0._f64
       ll2 = i+l2
       if (ll2>n) ll2 = n
       do j = (i+1),ll2 
          s = s + u(i,j)*x(j)
       enddo
       x(i) = (x(i)-s)/u(i,i)
    enddo

  end subroutine solvlub

end module pivotbande
