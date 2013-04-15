!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: gausslobatto
!
! DESCRIPTION:
!> @file gauss-lobatto.F90
!> @namespace gausslobatto
!> @author Madaule Eric
!> @brief Gauss-Lobatto interpolation tools
!> @details Here are several of the Gauss-Lobatto tools :
!>            ·Gauss-Lobatto points and weight,
!>            ·Gauss-Lobatto bases functions and the integral of their product,
!>            ·integral of product of Gauss-Lobatto function and their derivative.
!>          This module will first be limited to 1D and should extend as people will have the need for
!>          higher dimension (and so have time to write it).
!>          Note that for the construction in 1D the code is based on pseudo code given by David A.
!>          Kopriva in his book <EM>Implementing Spectral Methods for Partial Differential Equations</EM> .
!------------------------------------------------------------------------------
module gausslobatto_mod
#include "sll_working_precision.h"
#include "sll_memory.h"

  use sll_constants

  implicit none

  type gausslobatto1D
     !---------------------------------------------------------------------------
     !> @brief Gauss-Lobatto
     !> @details Gauss-Lobatto nodes and weith on a reference element [-1;1] in 1D.
     !>          This also contains Jacobian of the transformation between real element
     !>          and reference element
     !> 
     !>          A gausslobatto1d object contains node, weight, jac and degree.
     !>          TheseThe firsts are allocatable array. Use the construction subroutine to build
     !>          it. Only degree is a scalar. It is the degree of corresponding polynomials.
     !>          Only jac must be filled manually (but it is allocated in the constructor)
     sll_real64,dimension(:),pointer :: node,weight,jac
     sll_int32 :: degree
  end type gausslobatto1D

contains
  ! some routines here have a number
  ! this number is the number in the book of D. A. Kopriva :
  ! Implementing Spectral Methods for Partial Differential Equations

  subroutine init_gausslobatto_1d(size,gl_obj)
    !---------------------------------------------------------------------------
    !> @brief Construction of Gauss-Lobatto nodes and weights
    !> @details Construction of Gauss-Lobatto nodes and weights.
    !>          This routine fill the node and weight but you have to fill the Jacobian with
    !>          the transformation you consider
    !>          This should be extended as I will complete the type.
    !> @param[IN] size number of Gauss-Lobatto points, should be at least 2 (if 1, then piecewise constant)
    !> @param[OUT] gl_obj Gauss-Lobatto object to build

    sll_int32,intent(in) :: size
    type(gausslobatto1D),intent(out) :: gl_obj

    sll_int32 :: err

    if (size<1) then
       print*,"not enought points to build the Gauss-Lobatto interpolator"
       print*,"exiting..."
       stop
    end if

    SLL_ALLOCATE(gl_obj%node(size),err)
    SLL_ALLOCATE(gl_obj%weight(size),err)
    SLL_ALLOCATE(gl_obj%jac(size),err)

    gl_obj%degree=size-1

    call LGL_NodesAndWeight(gl_obj%node,gl_obj%weight)

  end subroutine init_gausslobatto_1d

  !intermediary procedure for computation of Gauss-Lobatto Nodes and Weights
  !algo 24
  subroutine qAndlEvaluation(n,x,q,qp,ln)

    sll_int32, intent(in) :: n
    sll_real64, intent(in) :: x

    sll_real64, intent(out) :: q,qp,ln

    sll_int32 :: k
    sll_real64 :: lm2,lm1,lpm2,lpm1,lnp,lp1,lpp1

    k=2
    lm2=1.0d0
    lm1=x
    lpm2=0.0d0
    lpm1=1.0d0

    do k=2,n
       ln=real(2*k-1,kind(x))/real(k,kind(x))*x*lm1-real(k-1,kind(x))/real(k,kind(x))*lm2
       lnp=lpm2+real(2*k-1,kind(x))*lm1
       lm2=lm1
       lm1=ln
       lpm2=lpm1
       lpm1=lnp
    end do
    k=n+1
    lp1=real(2*k-1,kind(x))/real(k,kind(x))*x*ln-real(k-1,kind(x))/real(k,kind(x))*lm1
    lpp1=lpm1+real(2*k-1,kind(x))*ln
    q=lp1-lm1
    qp=lpp1-lpm1

  end subroutine qAndlEvaluation

  !computation of Gauss-Lobatto nodes and weights
  !algo 25
  subroutine LGL_NodesAndWeight(node,weight)

    sll_real64, dimension(0:), intent(out) :: node,weight

    sll_int32 :: n
    !for loops
    sll_int32 :: j,k,nit
    sll_real64 :: tol
    !variables for calculus
    sll_real64 :: q,qp,ln,delta

    nit=10000
    tol=4.0d0*epsilon(real(1,kind(tol)))

    n=size(node)
    if (n/=size(weight) ) then
       print*,"array dimension error in LegendreGaussLobattoNodesAndWeight"
       print*,"in and out array must have the same size"
       print*,"exiting the program"
       stop
    end if
    n=n-1 !so arrays go from 0 to n

    if (n==1) then
       node(0)=-1.0d0
       node(1)=1.0d0
       weight(0)=1.0d0
       weight(1)=weight(0)
    else
       node(0)=-1.0d0
       weight(0)=2.0d0/real(n*(n+1),kind(tol))
       node(n)=1.0d0
       weight(n)=weight(0)
       
       do j=1,(n+1)/2-1
          node(j)=-cos( (real(j,kind(tol))+0.25d0)*sll_pi/real(n,kind(tol))-3.0d0/(real(8*n,kind(tol))*sll_pi)* &
               & 1.0d0/(real(j,kind(tol))+0.25d0) )
          
          k=0
          delta=tol*2.0d0+1.0d0
          do while ( k<=nit .and. abs(delta)>=tol*abs(node(j)) )
             call qAndlEvaluation(n,node(j),q,qp,ln)
             delta=-q/qp
             node(j)=node(j)+delta
             print*,node(j)
             k=k+1
          end do

          call qAndlEvaluation(n,node(j),q,qp,ln)
          node(n-j)=-node(j)
          weight(j)=2.0d0/(real(n*(n+1),kind(tol))*ln**2)
          weight(n-j)= weight(j)
       end do

       if (modulo(n,2)==0) then
          call qAndlEvaluation(n,0.0d0,q,qp,ln)
          node(n/2)=0.0d0
          weight(n/2)=2.0d0/(real(n*n+n,kind(tol))*ln**2)
       end if

    end if

  end subroutine LGL_NodesAndWeight

  subroutine delete_gausslobatto_1d(gl_obj)
    !---------------------------------------------------------------------------
    !> @brief delete an object of type gausslobatto1D
    !> @details delete all array in an object of type gausslobatto1D
    !> @param[INOUT] gl_obj the object to delete

    type(gausslobatto1D),intent(inout) :: gl_obj

    DEALLOCATE(gl_obj%node)
    DEALLOCATE(gl_obj%weight)
    DEALLOCATE(gl_obj%jac)

  end subroutine delete_gausslobatto_1d

  subroutine gl_massmatrix_1d(gl_obj,massmat)
    !---------------------------------------------------------------------------
    !> @brief Computation of diagonal mass matrix with Gauss-Lobatto points
    !> @details Computation of diagonal mass matrix with Gauss-Lobatto points.
    !> @param[INOUT] gl_obj the object to delete

    type(gausslobatto1D),intent(in) :: gl_obj
    sll_real64,dimension(gl_obj%degree+1),intent(out) :: massmat

    sll_int32 :: i

    do i=1,size(massmat)
       massmat(i)=gl_obj%weight(i)*gl_obj%jac(i)
    end do

  end subroutine gl_massmatrix_1d

end module gausslobatto_mod
