!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: poisson4dg
!
! DESCRIPTION:
!> @file vlasov_poisson_DG.F90
!! @namespace poisson4dg
!! @author Madaule Eric
!! @brief Tools for the resolution of the Vlasov-Poisson system with Discontinuous Galerkin.
!! @details Tools for the resolution of the Poisson equation with Discontinuous Galerkin.
!!          You must initialize your matrix(ces) with poisson1d_matrix. The resolution is done
!!          with UMFpack. You just need to call poisson_solve_4dg_1d_csc to solve the Poisson
!!          equation. If I have time I will add more flexibility to those.
!!
!!          This module will first be limited to 1D and should extend as people will
!!          have the need for higher dimension (and so have time to write it).
!!         
!------------------------------------------------------------------------------
module Poisson4dg
#include "sll_working_precision.h"

  use gausslobatto
  !those are part of FEMilaro
  use mod_sparse
  use mod_umfpack

  !use mod_octave_io_sparse

  implicit none

  type umfpack_plan
     !> plan for UMFpack
     !! This plan mainly contains array objects so they are not created
     !! every times we call UMFpack for the resolution of the Poisson problem.
     !! Note that it does not need to be initialized, it contains workink array of fixe
     !! size for UMFpack.
     sll_real64 :: control(umfpack_control), info(umfpack_info)
    integer(umf_void) :: symbolic
  end type umfpack_plan

  interface poisson_solve_4dg_1d
     module procedure poisson_solve_4dg_1d_csc
     !more could follow, it will depend on the time I have to do it
  end interface poisson_solve_4dg_1d

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
    d%ti=0
    d%tj=0
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

    !construction of vp_mat
    vp_mat=matmul(matmul(tri2col(d),tri2col(m_inv)),tri2col(transpose(d)))+c
    field_mat=matmul(tri2col(m_inv),tri2col(transpose(d)))

    field_mat=matmul(tri2col(new_tri( ne*ng,ne*ng, (/(i-1,i=1,ne*ng)/), (/(i-1,i=1,ne*ng)/), &
         & (/ ((1.0d0/jac(i),j=1,ng),i=1,ne) /) )), field_mat)

    !we impose phi_0=alpha (because the matrix is periodic) (we set the cooresponding line
    !to (0,..,0,1,0,..,0)
    !and corresponding value of rhs must be alpha)
    !alpha is set outside; the coordinate of phi_0 is set with x_bound

    if (x_bound>vp_mat%m .or. x_bound>vp_mat%n .or. x_bound<1) then
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

  subroutine poisson_solve_4dg_1d_csc(plan,a,b,x)
    !< @brief Resolution of the Poisson problem in 1D using discontinous Galerkin/
    !!        spectral element method
    !! @details This routine solve the Poisson problem AX=B for dg using UMFpack.
    !!          The matrix a is sparse, A and B must be initialized, X must be allocated before.
    !!          If this routine returns error, see UMFpack documentation at
    !!          http://www.cise.ufl.edu/research/sparse/umfpack/UMFPACK/Doc/UserGuide.pdf
    !!          In the case of DG, if you define field for v>0 and v<0 you will need to call
    !!          this rountine twice.
    !! @param[IN] plan UMFpack plan, contains array for UMFpack
    !! @param[IN] a t_col matrix (CSC matrix)
    !! @param[IN] b rhs of system AX=B, 1D array
    !! @param[OUT] x unknown of system AX=B, 1D array

    implicit none

    !umfpack variables
    type(umfpack_plan),intent(inout) :: plan
    type(t_col),intent(in) :: a
    sll_real64,dimension(:),intent(in) :: b
    sll_real64,dimension(:),intent(out) :: x

    integer :: sys,sx
    sll_int64 :: numeric

    sx=size(x)
    if (sx/=a%m .or. size(b)/=a%n .or. a%m/=a%n) then
       print*,'error in array or matrix size for the resolution of Poisson problem'
       print*,'expected square problem, passed',a%m,'×',a%n,'problem'
       print*,'expected size(x)=',a%m,' ; expected size(b)=',a%n
       print*,'passed',sx,'for x and',size(b),'for b'
       print*,'exiting...'
       stop
    end if

    sys=UMFPACK_A

    call umf4def(plan%control)
    plan%control(1)=2
    call umf4sym (int(sx,umf_int),int(sx,umf_int),a%ap,a%ai,a%ax, &
         & plan%symbolic,plan%control,plan%info)
    if (plan%info(1) .lt. 0) then
       print *, 'Error occurred in umf4sym: ', plan%info (1)
       stop
    end if
    !call umf4pinf(plan%control,plan%info);stop ! this is only for debuging
    call umf4num (a%ap,a%ai,a%ax,plan%symbolic,numeric,plan%control,plan%info)
    call umf4solr(sys,a%ap,a%ai,a%ax,x,b,numeric,plan%control,plan%info)

  end subroutine poisson_solve_4dg_1d_csc

end module Poisson4DG
