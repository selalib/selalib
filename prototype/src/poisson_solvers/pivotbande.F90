module pivotbande
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
  implicit none
contains

  subroutine pib(a,x,b,n,l1,l2)
    implicit none
    sll_real64,dimension(:,:),intent(in)::a
    sll_real64,dimension(:),  intent(in)::b
    sll_real64,dimension(:), intent(out)::x
    sll_real64,dimension(n,n)           ::aa
    sll_real64,dimension(n)             ::bb
    sll_int32,intent(in)                ::n,l1,l2
    sll_int32                           ::i,j,k,ll1,ll2
    sll_real64                          ::p,s

    aa = a
    bb = b

    do k = 1,n-1
       ll1 = k + l1;
       if (ll1>n)  ll1 = n 
       do i = k+1,ll1  
          p = aa(i,k)/aa(k,k)             
          ll2 = k + l2
          if (ll2>n) ll2=n 
          do j = k+1,ll2 
             aa(i,j) = aa(i,j) - aa(k,j)*p 
          enddo
          bb(i) = bb(i) - bb(k) * p 
       enddo
       aa(k+1:n,k)=0._f64
    enddo

    !remontée

    x(n) = bb(n)/aa(n,n)

    do i=n-1,1,-1
       s=0.
       ll2 = i + l2
       if (ll2>n) ll2 = n
       do j = (i+1),ll2 
          s = s + aa(i,j)*x(j)
       enddo
       x(i) = (bb(i)-s)/aa(i,i)
    enddo

  endsubroutine pib


  subroutine searchband(a,n,l1,l2)
    implicit none
    sll_real64,dimension(:,:)::a
    sll_int32,intent(out)    ::l1,l2
    sll_int32,dimension(n)   ::b1,b2
    sll_int32::i,j,n,ll1,ll2!,width

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


  subroutine solvaxbameliore(a,x,b,n)
    implicit none
    sll_real64,dimension(:,:) ,intent(in)::a
    sll_real64,dimension(:)   ,intent(in)::b
    sll_real64,dimension(:),intent(inout)::x
    sll_int32,                 intent(in)::n
    sll_real64,dimension(n,n)            ::c
    sll_real64,dimension(n)              ::d,e,f
    sll_int32 ::i,j,k,l,m,p,q
    sll_real64::s,u,v, piv

    ! on résout ax=b en créant un système triangulaire équivalent
    !méthode du pivot de gauss
    c=a
    d=b

    do k=1,n-1
       q=k+1

       if (abs(c(k,k))<1e-12) then
          ! trouver le cik le + grand en val abs
          print*,"oho"
          do i=k+1,q
             p=i
             do j=i,p
                if (p==n) then
                   e=c(k,:)
                   f=c(i,:)
                   u=d(k)
                   v=d(i)
                   ! chgt de lignes
                   c(k,:)=f
                   c(i,:)=e
                   d(k)=v
                   d(i)=u
                else
                   if ((c(i,k)**2)>=(c(j,k)**2)) then
                      p=p+1
                   else 
                      q=q+1 
                   endif
                endif
             enddo
          enddo

       else

          do i=k+1,n

             !if (c(i,k)/=0) then

             piv = c(i,k)/c(k,k)
             d(i)=d(i)-d(k)*piv
             do j=n,k,-1
                if (c(k,j)/=0) then 
                   c(i,j)=c(i,j)-c(k,j)*piv
                endif
             enddo

             !endif

          enddo

       end if

    enddo


    ! remontée
    x(n)=d(n)/c(n,n)

    do i=n-1,1,-1
       s=0.
       do j=(i+1),n
          s=s+c(i,j)*x(j)
       enddo
       x(i)=(d(i)-s)/c(i,i)
    enddo

  endsubroutine solvaxbameliore



  !***********************************************


  subroutine factolub (a,l,u,n,l1,l2)
    implicit none
    sll_real64,dimension(:,:),intent(in)::a
    sll_real64,dimension(:,:),intent(out)::l,u
    sll_int32,intent(in)::n,l1,l2
    sll_int32::i,j,k,ll1,ll2

    ! init
    u = a
    l = 0._f64
    do i = 1,n
       l(i,i) = 1._f64
    enddo

    ! transformation

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
