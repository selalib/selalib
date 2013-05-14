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
#include "sll_integration.h"

  use mod_sparse

  implicit none

contains

  subroutine poisson1d_matrix(gausslob,ne,jac,c11,c12,x_bound,VP_mat2,field_mat)
    !---------------------------------------------------------------------------
    !< @brief Construction of the matrix for the resolution of Poisson equation in 1D 
    !!        for Vlasov-Poisson
    !! @details Construction of the matrix for the resolution of Poisson equations in 1D 
    !!          for Vlasov-Poisson.
    !!          With this routine we compute the matrix to solve
    !!                            -Laplacian(Phi) = rho
    !!          with periodic boundary conditions
    !!          We set the constant with x_bound and bc
    !! @param[IN] gausslob Gauss-Lobatto initialized object
    !! @param[IN] jac vector of jacobian on each element
    !! @param[IN] c11,c12,c22 flux coefficients
    !!                        WARNING : c22 is unused <=> c22=0
    !!                        Curently c22 is not requiered and fixed to 0
    !! @param[IN] x_bound node number such as Phi(x_bound)=bc, where bc is a value the should
    !!                    be set in rho
    !!                    WARNING : 1 <= x_bound <=ne*ng, ng : number of GLL points/element
    !! @param[OUT] VP_mat2 CSC matrix for the resolution of Poisson equations in 1D for 
    !!                     Vlasov-Poisson, type t_col

    implicit none

    type(gausslobatto1d),intent(in) :: gausslob
    sll_int32, intent(in) :: ne
    sll_real64,dimension(:),intent(in) :: jac
    sll_real64,intent(in) :: c11,c12!,c22
    sll_int32,intent(in) :: x_bound
    type(t_col),intent(out) :: VP_mat2,field_mat

    !type(skyline_matrix) :: vp_mat3

    type(t_col) :: VP_mat
    type(t_tri) :: m_inv,c,d
    sll_int32 :: i,j,k,ng

    ng=gausslob%degree+1

    !mass matrix
    m_inv=new_tri(ng*ne,ng*ne,ng*ne)
    do j=0,ne-1
       do i=1,ng
          m_inv%ti(j*ng+i)=j*ng+i-1
          m_inv%tj(j*ng+i)=j*ng+i-1
          m_inv%tx(j*ng+i)=1.0d0/gausslob%weigh(i)*jac(j+1)
       end do
    end do

    !matrix C
    c=new_tri(ng*ne,ng*ne,4*ne)

    c%ti(1)=0
    c%tj(1)=0
    c%tx(1)=c11*jac(1)

    c%ti(2)=ng-1
    c%tj(2)=ng-1
    c%tx(2)=c11*jac(1)

    c%ti(3)=0
    c%tj(3)=ne*ng-1
    c%tx(3)=-c11*jac(1)

    c%ti(4)=ng-1
    c%tj(4)=ng
    c%tx(4)=-c11*jac(1)
    do i=2,ne-1
       c%ti((i-1)*4+1)=(i-1)*ng
       c%tj((i-1)*4+1)=(i-1)*ng
       c%tx((i-1)*4+1)=c11*jac(i)

       c%ti((i-1)*4+2)=i*ng-1
       c%tj((i-1)*4+2)=i*ng-1
       c%tx((i-1)*4+2)=c11*jac(i)

       c%ti((i-1)*4+3)=(i-1)*ng
       c%tj((i-1)*4+3)=(i-1)*ng-1
       c%tx((i-1)*4+3)=-c11*jac(i)

       c%ti(i*4)      =i*ng-1
       c%tj(i*4)      =i*ng
       c%tx(i*4)      =-c11*jac(i)
    end do
    c%ti(4*ne-3)=(ne-1)*ng
    c%tj(4*ne-3)=(ne-1)*ng
    c%tx(4*ne-3)=c11*jac(ne)

    c%ti(4*ne-2)=ne*ng-1
    c%tj(4*ne-2)=ne*ng-1
    c%tx(4*ne-2)=c11*jac(ne)

    c%ti(4*ne-1)=(ne-1)*ng
    c%tj(4*ne-1)=(ne-1)*ng-1
    c%tx(4*ne-1)=-c11*jac(ne)

    c%ti(4*ne)  =ne*ng-1
    c%tj(4*ne)  =0
    c%tx(4*ne)  =-c11*jac(ne)

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
             d%tx((i-1)*ng**2+(j-1)*ng+k)=gausslob%der(k,j)*jac(i)
             if (j==1 .and. k==1) then
                d%tx((i-1)*ng**2+(j-1)*ng+k)=d%tx((i-1)*ng**2+(j-1)*ng+k)-(-0.5d0+c12)*jac(i)
             else if (j==ng .and. k==ng) then
                d%tx((i-1)*ng**2+(j-1)*ng+k)=d%tx((i-1)*ng**2+(j-1)*ng+k)-(0.5d0+c12)*jac(i)
             end if
          end do
       end do
    end do

    d%ti(ne*ng**2+1)=ng
    d%tj(ne*ng**2+1)=ng+1
    d%tx(ne*ng**2+1)=-(0.5d0-c12)*jac(1)

    d%ti(ne*(ng**2+2))=1
    d%tj(ne*(ng**2+2))=ne*ng
    d%tx(ne*(ng**2+2))=-(-0.5d0-c12)*jac(1)

    do i=2,ne-1
       d%ti(ne*ng**2+i)=i*ng
       d%tj(ne*ng**2+i)=i*ng+1
       d%tx(ne*ng**2+i)=-(0.5d0-c12)*jac(i)

       d%ti(ne*(ng**2+2)-i+1)=(i-1)*ng+1
       d%tj(ne*(ng**2+2)-i+1)=(i-1)*ng
       d%tx(ne*(ng**2+2)-i+1)=-(-0.5d0-c12)*jac(i)
    end do

    d%ti(ne*ng**2+ne)=ne*ng
    d%tj(ne*ng**2+ne)=1
    d%tx(ne*ng**2+ne)=-(0.5d0-c12)*jac(ne)

    d%ti(ne*ng**2+ne+1)=(ne-1)*ng+1
    d%tj(ne*ng**2+ne+1)=(ne-1)*ng
    d%tx(ne*ng**2+ne+1)=-(-0.5d0-c12)*jac(ne)

    d%ti=d%ti-1
    d%tj=d%tj-1

    !because the matrix is -d*m_inv*transpose(d)+c and we can't do A=-A
    !on object of type t_col or t_tri, it is easier to do m_inv=-m_inv
    !since d=-d is not define
    !m_inv%tx=-m_inv%tx

    !construction of vp_mat
    vp_mat=matmul(matmul(tri2col(d),tri2col(m_inv)),tri2col(transpose(d)))+c
    field_mat=matmul(tri2col(m_inv),tri2col(transpose(d)))

    !we impose phi_0=0 (because the matrix is periodic) (we set the first line to (1,0...0)
    !and first value of rhs must also be 0)

!!$    vp_mat2=matmul(tri2col(new_tri( ne*ng, ne*ng, (/(i-1,i=1,ne*ng)/), (/(i-1,i=1,ne*ng)/), &
!!$         & (/( real((i-1)/max(i-1,1),8),i=1,ne*ng )/) )),vp_mat) + &
!!$         & new_tri(ne*ng,ne*ng,(/0/),(/0/),(/1.0d0/))
    field_mat=matmul(tri2col(new_tri( ne*ng,ne*ng, (/(i-1,i=1,ne*ng)/), (/(i-1,i=1,ne*ng)/), &
         & (/ ((1.0d0/jac(i),j=1,ng),i=1,ne) /) )), field_mat)

    if (x_bound>vp_mat%m .or. x_bound>vp_mat%n) then
       print*,'boudary conditions out of domain'
       print*,x_bound,vp_mat%m,vp_mat%n
       print*,'no boundary conditions imposed'
    else
       vp_mat2=matmul(tri2col(new_tri( ne*ng, ne*ng, (/(i-1,i=1,ne*ng)/), (/(i-1,i=1,ne*ng)/), &
            & (/( nequal(i,x_bound),i=1,ne*ng )/) )),vp_mat) + &
            & new_tri(ne*ng,ne*ng,(/x_bound-1/),(/x_bound-1/),(/1.0d0/))
    end if

    call clear(vp_mat)
    call clear(d)
    call clear(m_inv)
    call clear(c)

  end subroutine poisson1d_matrix

  !only used for filling vp_mat2
  function nequal(a,b) result(c)

    sll_int32,intent(in) :: a,b
    sll_real64 :: c

    if (a==b) then
       c=0.0d0
    else
       c=1.0d0
    end if

  end function nequal

end module Poisson4DG
