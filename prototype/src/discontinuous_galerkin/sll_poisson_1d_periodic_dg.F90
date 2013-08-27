!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: poisson4dg
!
! DESCRIPTION:
!> @file vlasov_poisson_DG.F90
!! @author Madaule Eric
!! @brief Resolution of the Poisson's problem using discontinuous Galerkin.
!! @details Resolution of the Poisson's problem using discontinuous Galerkin. The 
!!          resolution is done using UMFpack. This module solve 
!!          \f[ -\Delta(\Phi)=f, \f]
!!          \f[ E=\nabla(\Phi). \f]
!!          with the corresponding matrix 
!!          \f[ \begin{pmatrix} M &,& D^T
!!                           \\-D &,& C
!!          \end{pmatrix}
!!          \begin{pmatrix} E \\ \Phi \end{pmatrix} =
!!          \begin{pmatrix} 0 \\ M.\rho \end{pmatrix} \f]
!!          according to article _Discontinuous Galerkin methods for the one-dimensional 
!!          Vlasov-Poisson system_, Blanca Ayuso, 2011.\n
!!          But you can of course just use the UMFpack part to solve your Poisson's problem.
!!
!!          This module will first be limited to 1D and should extend as people will
!!          have the need for higher dimension (and so have time to write it).
!!         
!------------------------------------------------------------------------------
module sll_poisson4dg_1d_periodic_dg
#include "sll_working_precision.h"

  use sll_gausslobatto
  !those are part of FEMilaro
  use mod_sparse
  use mod_umfpack

  !use mod_octave_io_sparse

  implicit none

  type,public :: poisson_1d_periodic_dg
     !> plan for the Poisson's problem using discontinuous Galerkin
     !! This plan mainly contains array objects so they are not created
     !! every times we call UMFpack for the resolution of the Poisson proble.
     !! Build it with new (see corresponding section) and delete it with
     !! delete().
     sll_real64 :: control(umfpack_control), info(umfpack_info)
     integer(umf_void) :: symbolic
     sll_int64 :: numeric
     type(t_col) :: mat_poisson,mat_field
  end type poisson_1d_periodic_dg

  interface new
     module procedure new_poisson_1d_periodic_dg
  end interface new

  interface delete
     module procedure delete_poisson_1d_periodic_dg
  end interface delete

  interface solve
     module procedure poisson_solve_4dg_1d_csc
     !more could follow, it will depend on the time I have to do it
  end interface solve

contains

  subroutine poisson1d_matrix(gausslob,ne,jac,c11,c12,x_bound,VP_mat0,field_mat)
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
    !! @param[IN] c11,c12 flux coefficients
    !! @param[IN] x_bound node number such as Phi(x_bound)=bc, where bc is a value the should
    !!                    be set in rho
    !!                    WARNING : 1 <= x_bound <=ne*ng, ng : number of GLL points/element
    !! @param[OUT] VP_mat0 CSC matrix for the resolution of Poisson equations in 1D for 
    !!                     Vlasov-Poisson, type t_col
    !! @param[OUT,OPTIONAL] field_mat matrix to compute 
    !!                                \f[ E=-\nabla\Phi \f]
    !!                                in Vlasov-Poisson's problem

    implicit none

    type(gausslobatto1d),intent(in) :: gausslob
    sll_int32, intent(in) :: ne
    sll_real64,dimension(:),intent(in) :: jac
    sll_real64,intent(in) :: c11,c12
    sll_int32,intent(in) :: x_bound
    type(t_col),intent(out) :: VP_mat0
    type(t_col),intent(out),optional :: field_mat

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
          m_inv%tx(j*ng+i)=jac(j+1)/gausslob%weigh(i)
       end do
    end do

    !matrix C
    c=new_tri(ng*ne,ng*ne,4*ne)

    c%ti(1)=0
    c%tj(1)=0
    c%tx(1)=c11

    c%ti(2)=ng-1
    c%tj(2)=ng-1
    c%tx(2)=c11

    c%ti(3)=0
    c%tj(3)=ne*ng-1
    c%tx(3)=-c11

    c%ti(4)=ng-1
    c%tj(4)=ng
    c%tx(4)=-c11
    do i=2,ne-1
       c%ti((i-1)*4+1)=(i-1)*ng
       c%tj((i-1)*4+1)=(i-1)*ng
       c%tx((i-1)*4+1)=c11

       c%ti((i-1)*4+2)=i*ng-1
       c%tj((i-1)*4+2)=i*ng-1
       c%tx((i-1)*4+2)=c11

       c%ti((i-1)*4+3)=(i-1)*ng
       c%tj((i-1)*4+3)=(i-1)*ng-1
       c%tx((i-1)*4+3)=-c11

       c%ti(i*4)      =i*ng-1
       c%tj(i*4)      =i*ng
       c%tx(i*4)      =-c11
    end do
    c%ti(4*ne-3)=(ne-1)*ng
    c%tj(4*ne-3)=(ne-1)*ng
    c%tx(4*ne-3)=c11

    c%ti(4*ne-2)=ne*ng-1
    c%tj(4*ne-2)=ne*ng-1
    c%tx(4*ne-2)=c11

    c%ti(4*ne-1)=(ne-1)*ng
    c%tj(4*ne-1)=(ne-1)*ng-1
    c%tx(4*ne-1)=-c11

    c%ti(4*ne)  =ne*ng-1
    c%tj(4*ne)  =0
    c%tx(4*ne)  =-c11

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
             d%tx((i-1)*ng**2+(j-1)*ng+k)=gausslob%der(k,j)
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

    !construction of vp_mat
    vp_mat=(matmul(matmul(tri2col(d),tri2col(m_inv)),tri2col(transpose(d)))+c)
    field_mat=matmul(tri2col(m_inv),tri2col(transpose(d)))
    !field_mat%ax=-field_mat%ax

    !we impose phi_0=alpha (because the matrix is periodic) (we set the cooresponding line
    !to (0,..,0,1,0,..,0)
    !and corresponding value of rhs must be alpha)
    !alpha is set outside; the coordinate of phi_0 is set with x_bound

    if (x_bound>vp_mat%m .or. x_bound>vp_mat%n .or. x_bound<0) then
       print*,'boudary conditions out of domain'
       print*,x_bound,vp_mat%m,vp_mat%n
       print*,'no boundary conditions imposed'
    else
       vp_mat0=matmul(tri2col(new_tri( ne*ng, ne*ng, (/(i-1,i=1,ne*ng)/), (/(i-1,i=1,ne*ng)/), &
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
    !basically, this function is the same as 1-delta(i,j) with delta the 
    !Kronecker symbol

    sll_int32,intent(in) :: a,b
    sll_real64 :: c

    if (a==b) then
       c=0.0d0
    else
       c=1.0d0
    end if

  end function nequal

  subroutine poisson_solve_4dg_1d_csc(plan,b,x)
    !< @brief Resolution of the Poisson problem in 1D using discontinous Galerkin/
    !!        spectral element method
    !! @details This routine solve the Poisson problem AX=B for dg using UMFpack.
    !!          All the data are incuded in plan
    !!          If this routine returns error, see UMFpack documentation at
    !!          http://www.cise.ufl.edu/research/sparse/umfpack/UMFPACK/Doc/UserGuide.pdf
    !!          In the case of DG, if you define field for v>0 and v<0 you will need to call
    !!          this rountine twice.
    !! @param[IN] plan poisson_1d_periodic_dg object, contains arrays and variables
    !! @param[IN] b rhs of system AX=B, 1D array
    !! @param[OUT] x unknown of system AX=B, 1D array

    implicit none

    !umfpack variables
    type(poisson_1d_periodic_dg),intent(inout) :: plan
    sll_real64,dimension(:),intent(in) :: b
    sll_real64,dimension(:),intent(out) :: x

    integer :: sys,sx

    sys=0

    sx=size(x)
    if (sx/=plan%mat_poisson%m .or. size(b)/=plan%mat_poisson%n .or. &
         & plan%mat_poisson%m/=plan%mat_poisson%n) then
       print*,'error in array or matrix size for the resolution of Poisson problem'
       print*,'expected square problem, passed',plan%mat_poisson%m,'×',plan%mat_poisson%n, &
            & 'problem'
       print*,'expected size(x)=',plan%mat_poisson%m,' ; expected size(b)=',plan%mat_poisson%n
       print*,'passed',sx,'for x and',size(b),'for b'
       print*,'exiting...'
       stop
    end if

    call umf4solr(sys,plan%mat_poisson%ap,plan%mat_poisson%ai,plan%mat_poisson%ax,x,b, &
         & plan%numeric,plan%control,plan%info)
    !call umf4pinf(plan%control,plan%info);stop ! this is only for debuging

  end subroutine poisson_solve_4dg_1d_csc

  subroutine new_poisson_1d_periodic_dg(plan,gll,ne,jac,c11,c12,x_bound,flag)
    !< @brief Initialization of poisson_1d_periodic_dg object.
    !! @details Initialization of poisson_1d_periodic_dg object.
    !! @param[OUT] plan poisson_1d_periodic_dg object to initialize
    !! @param[IN] gll gausslobatto1d object
    !! @param[IN] ne number of element of your mesh
    !! @param[IN] jac real 1D array, contains the jacobian of each element
    !! @param[IN] c11,c12 flux coefficients
    !! @param[IN] x_bound integer, point to set \f[ \Phi(x0)=\alpha \f] in Poisson's problem
    !! @param[IN] flag logical, set it to \code .true. \encode to get the matrix 
    !!                 \f[ M^{-1}D^T \f] for the Vlasov-Poisson's equations

    implicit none

    type(poisson_1d_periodic_dg),intent(out) :: plan
    type(gausslobatto1d),intent(in) :: gll
    sll_int32,intent(in) :: ne,x_bound
    sll_real64,dimension(:) :: jac
    sll_real64,intent(in) :: c11,c12
    logical,intent(in) :: flag

    if (flag) then
       call poisson1d_matrix(gll,ne,jac,c11,c12,x_bound,plan%mat_poisson,plan%mat_field)
    else
       call poisson1d_matrix(gll,ne,jac,c11,c12,x_bound,plan%mat_poisson)
    end if

    call umf4def(plan%control)
    plan%control(1)=2
    call umf4sym (int(plan%mat_poisson%m,umf_int),int(plan%mat_poisson%n,umf_int), &
         & plan%mat_poisson%ap,plan%mat_poisson%ai,plan%mat_poisson%ax, &
         & plan%symbolic,plan%control,plan%info)
    if (plan%info(1) .lt. 0) then
       print *,'Error occurred in umf4sym: ', plan%info (1)
       print*,'See documentation for more information'
       stop
    end if
    call umf4num (plan%mat_poisson%ap,plan%mat_poisson%ai,plan%mat_poisson%ax,plan%symbolic, &
         & plan%numeric,plan%control,plan%info)

  end subroutine new_poisson_1d_periodic_dg

  subroutine delete_poisson_1d_periodic_dg(plan)
    !< @brief Deletion of poisson_1d_periodic_dg object.
    !! @details Deletion of poisson_1d_periodic_dg object. You can just write \n
    !!          call delete(plan)
    !! @param[INOUT] plan poisson_1d_periodic_dg  object to delete

    implicit none

    type(poisson_1d_periodic_dg),intent(inout) :: plan

    call umf4fsym(plan%symbolic)
    call umf4fnum(plan%numeric)
    call clear(plan%mat_poisson)
    if (plan%mat_field%check()) then
       call clear(plan%mat_field)
    end if

  end subroutine delete_poisson_1d_periodic_dg

end module sll_poisson4dg_1d_periodic_dg
