!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: gausslobatto
!
! DESCRIPTION:
!> @file gauss-lobatto.F90
!! @author Madaule Eric
!! @brief Gauss-Lobatto interpolation tools
!! @details Here are several of the Gauss-Lobatto tools :\\
!!            ·Gauss-Lobatto points and weight,\\
!!            ·Gauss-Lobatto bases functions and the integral of their product,\\
!!            ·integral of product of Gauss-Lobatto function and their derivative.\\
!!          To use this module you must also link to the compilation gauss.f and lob.f
!!
!!          The mass matrix (which is the integral of \phi_i \times \phi_j) is simply 
!!          diag(weigh), so there is no need to store it more than just the weigh.\\
!!
!!          We also need the derivative matrix D.
!!          \f[ D_{i,j}=\int \phi_i \phi_j' \f]
!!
!!          This module will first be limited to 1D and should extend as people will
!!          have the need for higher dimension (and so have time to write it).
!!         
!------------------------------------------------------------------------------
module sll_gausslobatto
#include "sll_working_precision.h"

  use sll_constants
  
  implicit none

  type gausslobatto1D
     !---------------------------------------------------------------------------
     !< @brief Gauss-Lobatto
     !! @details Gauss-Lobatto nodes and weigh on a reference element [-1;1] in 1D.
     !!          This also includes the degree of polynomials and the matrix of derivatives
     !!          D_{i,j}=\int(\Psi_i\Psi'_j)=w_i\sum_{l/=j}1/(x_j-x_l)\prod_{m/=j,m/=l}(x_i-x_m)/(x_j-x_m)
     !! 
     !!          A gausslobatto1d object contains node, weight, jac and degree.
     !!          These are allocatable array. Use the construction subroutine
     !!          to build it. Only degree is a scalar. It is the degree of corresponding
     !!          polynomials. Only jac must be filled manually (but it is allocated in
     !!          the constructor)
     sll_real64,dimension(:),pointer :: node,weigh
     sll_int32 :: degree
     sll_real64,dimension(:,:),pointer :: der
  end type gausslobatto1D

  interface delete
     module procedure delete_gausslobatto_1D
  end interface delete

contains

  subroutine init_gausslobatto_1d(size,gl_obj)
    !---------------------------------------------------------------------------
    !< @brief Construction of Gauss-Lobatto nodes and weights
    !! @details Construction of Gauss-Lobatto nodes and weights.
    !!          This routine fill the node and weight but you have to fill the Jacobian with
    !!          the transformation you consider
    !!          This should be extended as I will complete the type.
    !!          I can't deal very well with fortran 77 interface. To change between simple
    !!          and double precision you have to go into the file and change some comments
    !!          (very easy, it is explained in the file) and compile again (sorry)
    !! @param[IN] size number of Gauss-Lobatto points, should be at least 2 (if 1,
    !!                 then piecewise constant => Gauss-Lobatto is impossible)
    !! @param[OUT] gl_obj Gauss-Lobatto object to build

    sll_int32,intent(in) :: size
    type(gausslobatto1D),intent(out) :: gl_obj

    sll_int32 :: err,i
    sll_real64,dimension(0:size-1) :: alpha,beta
    sll_real64,dimension(size) :: e,a,b

    if (size<=1) then
       print*,"not enought points to build the Gauss-Lobatto interpolator"
       print*,"exiting..."
       stop
    end if

    ALLOCATE(gl_obj%node(size))
    ALLOCATE(gl_obj%weigh(size))
    ALLOCATE(gl_obj%der(size,size))

    gl_obj%node=0.0d0
    gl_obj%weigh=0.0d0
    gl_obj%degree=size-1

    do i=0,size-1
       alpha(i)=0.0d0
       beta(i)=real(i,kind(gl_obj%degree))**2/((2.0d0*i+1.0d0)*(2.0d0*i-1.0d0))
    end do
    !for Gauss-Legendre and Gauss-Lobatto, beta(0)=int(dlambda)
    !see Algorithm xxx - ORTHPOL: A package of routines for  generating orthogonal
    !polynomials and Gauss-type quadrature rules by _Walter Gautschi_
    beta(0)=2.0d0

    !for single precision, comment the first line and uncomment the second
    !for double precision, comment the second line and uncomment the first
    call dlob(size-2,alpha,beta,-1.0d0,1.0d0,gl_obj%node,gl_obj%weigh,err,e,a,b)
    !call lob(size-2,alpha,beta,-1.0,1.0,gl_obj%node,gl_obj%weigh,err,e,a,b)

    ! this compute the Gauss points and weigh on [-1,1] for "size" points
    !call dgauss(size,alpha,beta,epsilon(1.0d0),gl_obj%node,gl_obj%weigh,err,e,a,b)

    call derivative_matrix_1d(gl_obj)

  end subroutine init_gausslobatto_1d

  subroutine delete_gausslobatto_1d(gl_obj)
    !---------------------------------------------------------------------------
    !< @brief delete an object of type gausslobatto1D
    !! @details delete all array in an object of type gausslobatto1D
    !! @param[INOUT] gl_obj the object to delete

    type(gausslobatto1D),intent(inout) :: gl_obj

    DEALLOCATE(gl_obj%node)
    DEALLOCATE(gl_obj%weigh)
    DEALLOCATE(gl_obj%der)

  end subroutine delete_gausslobatto_1d

  subroutine derivative_matrix_1d(gl_obj)
    !called by Gauss-Lobatto 1D constructor
    !---------------------------------------------------------------------------
    !< @brief construction of the derivative matrix for Gauss-Lobatto 1D
    !! @details Construction of the derivative matrix for Gauss-Lobatto 1D,
    !!          The matrix must be already allocated of size (number of point)^2.
    !!          der(i,j)=int(Phi_i.Phi'_j)_[-1;1]²
    !!                  =w_i.Phi'_j(x_i)
    !! @param[INOUT] gl_obj gausslobatto1D object to build derivative

    type(gausslobatto1D),intent(inout) :: gl_obj

    sll_int32 :: nb_pts,i,j,l,m
    sll_real64 :: prod

    nb_pts=gl_obj%degree+1

    gl_obj%der=0.0d0

    !loop on all element of D
    !loop on columns
    do j=1,nb_pts
       !loop on rows
       do i=1,nb_pts
          !loop on all the derivatives
          !the code is writen so there is no if
          do l=1,j-1
             prod=1.0d0
             do m=1,l-1!min(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=l+1,j-1!min(j,l)+1,max(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=j+1,nb_pts!max(j,l)+1,nb_pts
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             prod=prod/(gl_obj%node(j)-gl_obj%node(l))
             gl_obj%der(i,j)=gl_obj%der(i,j)+prod
          end do
          do l=j+1,nb_pts
             prod=1.0d0
             do m=1,j-1!min(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=j+1,l-1!min(j,l)+1,max(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=l+1,nb_pts!max(j,l)+1,nb_pts
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             prod=prod/(gl_obj%node(j)-gl_obj%node(l))
             gl_obj%der(i,j)=gl_obj%der(i,j)+prod
          end do
          gl_obj%der(i,j)=gl_obj%der(i,j)*gl_obj%weigh(i)
       end do
    end do

  end subroutine derivative_matrix_1d

end module sll_gausslobatto
