!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: poisson4dg
!
! DESCRIPTION:
!> @file vlasov_poisson_DG.F90
!! @namespace poisson4dg
!! @author Madaule Eric
!! @brief tools for the resolution of the Vlasov-Poisson system with Discontinuous Galerkin
!! @details Tools for the resolution of the Vlasov-Poisson system with Discontinuous Galerkin.
!!          Here is the initialization. The time step tools should arrive later.
!!
!!          This module will first be limited to 1D and should extend as people will
!!          have the need for higher dimension (and so have time to write it).
!!         
!------------------------------------------------------------------------------
module Poisson4dg
#include "sll_working_precision.h"

  use gausslobatto
  use mod_sparse
  implicit none

contains

  subroutine poisson1d_matrix(gausslob,ne,jac,c11,c12,c22,VP_mat)
    !---------------------------------------------------------------------------
    !< @brief Construction of the matrix for the resolution of Poisson equation in 1D 
    !!        for Vlasov-Poisson
    !! @details Construction of the matrix for the resolution of Poisson equations in 1D 
    !!          for Vlasov-Poisson.
    !!          
    !! @param[IN] gausslob Gauss-Lobatto initialized object
    !! @param[IN] jac vector of jacobian on each element
    !! @param[IN] c11,c12,c22 flux coefficients
    !!                        WARNING : c22 is unused <=> c22=0
    !! @param[OUT] VP_mat CSC matrix for the resolution of Poisson equations in 1D for 
    !!                    Vlasov-Poisson, type t_col

    implicit none

    type(gausslobatto1d),intent(in) :: gausslob
    sll_int32, intent(in) :: ne
    sll_real64,dimension(:),intent(in) :: jac
    sll_real64,intent(in) :: c11,c12,c22
    type(t_col),intent(out) :: VP_mat

    type(t_tri) :: m_inv,c,d
    sll_int32 :: i,j,k,ng

    ng=gausslob%degree+1

    !mass matrix
    m_inv=new_tri(ng*ne,ng*ne,ng*ne)
    do j=0,ne-1
       do i=1,gausslob%degree+1
          m_inv%ti(j*ng+i)=j*ng+i-1
          m_inv%tj(j*ng+i)=j*ng+i-1
          m_inv%tx(j*ng+i)=1.0d0/gausslob%weigh(i)*jac(j)
       end do
    end do

    !matrix C
    c=new_tri(ng*ne,ng*ne,4*ne)
    c%ti(1)=0
    c%tj(1)=0
    c%tx(1)=-c11*jac(1)
    c%ti(2)=ng-1
    c%tj(2)=ng-1
    c%tx(2)=c11*jac(1)
    c%ti(3)=0
    c%tj(3)=ne*ng-1
    c%tx(3)=c11*jac(1)
    c%ti(4)=ng-1
    c%tj(4)=ng
    c%tx(4)=-c11*jac(1)
    do i=2,ne-1
       c%ti((i-1)*4+1)=(i-1)*ng
       c%tj((i-1)*4+1)=(i-1)*ng
       c%tx((i-1)*4+1)=-c11*jac(i)
       c%ti((i-1)*4+2)=i*ng-1
       c%tj((i-1)*4+2)=i*ng-1
       c%tx((i-1)*4+2)=c11*jac(i)
       c%ti((i-1)*4+2)=(i-1)*ng
       c%tj((i-1)*4+3)=(i-1)*ng-1
       c%tx((i-1)*4+3)=c11*jac(i)
       c%ti(i*4)      =i*ng-1
       c%tj(i*4)      =i*ng
       c%tx(i*4)      =-c11*jac(i)
    end do
    c%ti(4*(ne-1))=(ne-1)*ng
    c%tj(4*(ne-1))=(ne-1)*ng
    c%tx(4*(ne-1))=-c11*jac(ne)
    c%ti(4*ne-2)  =ne*ng
    c%tj(4*ne-2)  =ne*ng
    c%tx(4*ne-2)  =c11*jac(ne)
    c%ti(4*ne-1)  =(ne-1)*ng
    c%tj(4*ne-1)  =(ne-1)*ng-1
    c%tx(4*ne-1)  =c11*jac(ne)
    c%ti(4*ne)    =ne*ng-1
    c%tj(4*ne)    =0
    c%tx(4*ne)    =-c11*jac(ne)

    !matrix D
    d=new_tri(ne*ng,ne*ng,ne*(ng**2+2))
    d%ti=2
    d%tj=2
    d%tx=0.0d0

    do i=1,ne
       !line
       do j=1,ng
          !column
          do k=1,ng
             d%ti((i-1)*ng**2+(j-1)*ng+k)=(i-1)*ng+j
             d%tj((i-1)*ng**2+(j-1)*ng+k)=(i-1)*ng+k
             d%tx((i-1)*ng**2+(j-1)*ng+k)=gausslob%der(j,k)
             if (j==1 .and. k==1) then
                d%tx((i-1)*ng**2+(j-1)*ng+k)=d%tx((i-1)*ng**2+(j-1)*ng+k)-(-0.5d0+c12)
             else if (j==ng .and. k==ng) then
                d%tx((i-1)*ng**2+(j-1)*ng+k)=d%tx((i-1)*ng**2+(j-1)*ng+k)-(0.5d0+c12)
             end if
          end do
       end do
    end do

    d%ti(ne*ng**2+1)=ng
    d%tj(ne*ng**2+1)=ng+1
    d%tx(ne*ng**2+1)=-(0.5d0-c12)

    d%ti(ne*(ng**2+2))=1
    d%tj(ne*(ng**2+2))=ne*ng
    d%tx(ne*(ng**2+2))=-(-0.5d0-c12)

    do i=2,ne-1
       d%ti(ne*ng**2+i)=i*ng
       d%tj(ne*ng**2+i)=i*ng+1
       d%tx(ne*ng**2+i)=-(0.5d0-c12)

       d%ti(ne*(ng**2+2)-i+1)=(i-1)*ng+1
       d%tj(ne*(ng**2+2)-i+1)=(i-1)*ng
       d%tx(ne*(ng**2+2)-i+1)=-(-0.5d0-c12)
    end do

    d%ti(ne*ng**2+ne)=ne*ng
    d%tj(ne*ng**2+ne)=1
    d%tx(ne*ng**2+ne)=-(0.5d0-c12)

    d%ti(ne*ng**2+ne+1)=(ne-1)*ng+1
    d%tj(ne*ng**2+ne+1)=(ne-1)*ng
    d%tx(ne*ng**2+ne+1)=-(-0.5d0-c12)

    d%ti=d%ti-1
    d%tj=d%tj-1

    !because the matrix is -d/m_inv*transpose(d)+c and we can't do A=-A
    !on object of type t_col or t_tri, it is easier to do m_inv=-m_inv
    m_inv%tx=-m_inv%tx

    !construction of vp_mat
    vp_mat=matmul(matmul(tri2col(d),tri2col(m_inv)),tri2col(transpose(d)))+c

  end subroutine poisson1d_matrix

end module Poisson4DG
