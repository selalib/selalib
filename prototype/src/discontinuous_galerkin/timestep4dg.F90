!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: poisson4dg
!
! DESCRIPTION:
!> @file timestep4dg.F90
!! @namespace timestep4dg
!! @author Madaule Eric
!! @brief Time step tools for discontinous Galerkin
!! @details This module contains the time steping for the Vlasov equation using 
!!          discontinuous Galerkin. The initialization is done with poisson4dg (should 
!!          be transfert to Poisson library, currently in file vlasov_poisson_DG.F90).
!!
!!          This module will first be limited to 1D and should extend as people will
!!          have the need for higher dimension (and time to write it).
!!         
!------------------------------------------------------------------------------
module timestep4dg
#include "sll_working_precision.h"

  use gausslobatto
  use poisson4dg
  use sll_nu_cart_mesh
  use sll_constants

  implicit none

  integer,parameter :: SLL_RK3=0,SLL_RK4=1

  type dg_time_steping
     !< @brief Plan for time stepping with discontinuous Galerkin.
     !! @details Plan for time stepping with discontinuous Galerkin.
     !!          It contains a working array, the alectrostatic and electric fields, the charge
     !!          distribution, the right hand side of Vlasov equation, the time step, the node
     !!          to apply \Phi(x=0)=alpha and array for UMFpack so they are not created every
     !!          time. Must be build with init_timesteping_4dg (see corresponding section) and
     !!          deleted with clear
     sll_real64,dimension(:,:,:),allocatable :: k
     sll_real64,dimension(:),allocatable :: phi,field,rho
     sll_real64,dimension(:,:),allocatable :: rhs
     type(t_col) :: matvp,matvm ! to compute phi
     type(t_col) :: fieldvp,fieldvm ! to compute E
     sll_real64 :: dt,t
     sll_int32 :: x0 ! place to apply phi(x=0)=0
     !umfpack variables
     type(umfpack_plan) :: umfpack_data
  end type dg_time_steping

  interface clear
     module procedure delete_dg_step
  end interface clear

contains

  subroutine init_timesteping_4dg(plan,method,gll,jac,dt,xbound,nx,nv,ng,c11,c12)
    !< @brief This routine is the interface to use to build interpolators for DG
    !! @details This routine only is the interface to use to build interpolators for DG.
    !!          It will extend as time stepping method will be added
    !! @param[OUT] plan dg_time_steping object to build
    !! @param[IN] method SLL_RK3 or SLL_RK4 to choose your interpolator
    !! @param[IN] gll gausslobatto1d object
    !! @param[IN] jac real array with Jacobian for each cell to convert from reference element
    !! @param[IN] dt time step, real double precision
    !! @param[IN] xbound coordinate to set \Phi(x=0)=alpha, integer, see Poisson4DG
    !! @param[IN] c11 coefficient c11 for flux
    !! @param[IN] c12 coefficient c12(v>0) for flux with assumption that c12(v>0)=-c12(v<0)

    implicit none

    type(dg_time_steping),intent(out) :: plan
    type(gausslobatto1d),intent(in) :: gll
    sll_real64,dimension(:),intent(in) :: jac
    sll_int32,intent(in) :: method,xbound
    sll_real64,intent(in) :: dt,c11,c12
    sll_int32,intent(in) :: nx,nv,ng

    if (method==0) then
       !RK3
       print*,'RK3 not done'
       stop
    else if(method==1) then
       !RK4
       call init_rk4_4dg(plan,gll,jac,dt,xbound,nx,nv,ng,c11,c12)
    end if

  end subroutine init_timesteping_4dg

  subroutine init_rk4_4dg(plan,gll,jac,dt,xbound,nx,nv,ng,c11,c12)
    !< @brief Do not use it, call the interface routine init_timesteping_4dg

    implicit none

    type(dg_time_steping),intent(out) :: plan
    type(gausslobatto1d),intent(in) :: gll
    sll_real64,dimension(:),intent(in) :: jac
    sll_int32,intent(in) :: xbound
    sll_real64,intent(in) :: dt,c11,c12
    sll_int32,intent(in) :: nx,nv,ng

    allocate(plan%k(2,nx*ng,nv*ng),plan%phi(nx*ng),plan%rho(nx*ng),plan%field(nx*ng), &
         & plan%rhs(nx*ng,nv*ng))
    plan%dt=dt
    plan%x0=xbound

    !build the matrixes for the Poisson-problem
    !((D+F_E).M-¹.(D-F_\Phi)^T - C).\Phi = M.(\rho-1)
    !this is for v>0
    call poisson1d_matrix(gll,nx,jac,c11,c12,xbound,plan%matvp,plan%fieldvp)
    !this is for v<0
    call poisson1d_matrix(gll,nx,jac,c11,-c12,xbound,plan%matvm,plan%fieldvm)

  end subroutine init_rk4_4dg

  subroutine delete_dg_step(plan)

    implicit none

    type(dg_time_steping),intent(inout) :: plan

    deallocate(plan%k,plan%field,plan%phi,plan%rho,plan%rhs)
    call clear(plan%matvm)
    call clear(plan%matvp)
    call clear(plan%fieldvm)
    call clear(plan%fieldvp)
    
  end subroutine delete_dg_step

  subroutine rhs4dg_1d(mesh,gll,field_e,dist,rhs,t)
    !< @brief Computation of rhs for Vlasov equation in 2D phase space (1D physical space)
    !! @details Computation of rhs for Vlasov equation in 2D phase space (1D physical space).
    !!          All variable must already be initialized
    !! @param[IN] mesh non_unif_cart_mesh object, gives the Jacobian on each cell and the 
    !!                 number of step
    !! @param[IN] gll gausslobatto1D object, gives most of the needed information for computation
    !! @param[IN] field_e electric field, 1D array of real
    !! @param[IN] dist distribution function, 2D array of real
    !! @param[OUt] rhs right hand side of Vlasov equation, 2D array of real

    implicit none

    type(non_unif_cart_mesh),intent(in) :: mesh
    type(gausslobatto1D),intent(in) :: gll
    sll_real64,dimension(:),intent(in) :: field_e
    sll_real64,dimension(:,:),intent(in) :: dist
    sll_real64,dimension(:,:),intent(out) :: rhs
    
    sll_real64,intent(in) :: t

    sll_int32 :: nx,nv,ng
    sll_int32 :: x1,x2,v1,v2,i
    sll_real64 :: som1,som2,x,v

    ng=gll%degree+1
    nx=mesh%n_etat1
    nv=mesh%n_etat2

    if (size(field_e)/=nx*ng .or. size(dist,1)/=nx*ng .or. size(dist,2)/=nv*ng) then
       print*,'size error in computation of rhs'
       print*,'please, check your size to finc the problem'
       print*,'exiting...'
    end if

    !construction or rhs
    !rhs=-v.d_x(f)+d_x(Phi).d_v(f), with d=\partial and d_x(Phi)=E
    rhs=0.0d0
    do v1=1,nv !loop on elements in direction v
       do v2=1,ng !loop on GLL points in direction v
          v=mesh%etat2(v1)+(1.0d0+gll%node(v2))/mesh%jac(nx+1,v1)
          do x1=1,nx !loop on elements in direction x
             do x2=1,ng !loop on GLL points in direction x
                x=mesh%etat1(x1)+(1.0d0+gll%node(x2))/mesh%jac(x1,nv+1)

                !interior part
                som1=0.0d0
                som2=0.0d0

                do i=1,ng
                   som1=som1+dist((x1-1)*ng+i,(v1-1)*ng+v2)*gll%der(i,x2)
                   som2=som2+dist((x1-1)*ng+x2,(v1-1)*ng+i)*gll%der(i,v2)
                end do
                som1=som1*mesh%jac(x1,nv+1)
                som2=som2*mesh%jac(nx+1,v1)
                rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=(som1*gll%weigh(v2)*v- &
                     & som2*gll%weigh(x2)*field_e((x1-1)*ng+x2))/mesh%jac(x1,v1)

                !boudaries part
                if (x2==ng) then 
                   if (mesh%etat2(v1)+(1.0d0+gll%node(v2))/mesh%jac(nx+1,v1) >= 0.0d0) then
                      rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                           & v*gll%weigh(v2)*dist(x1*ng,(v1-1)*ng+v2)/mesh%jac(nx+1,v1)
                   else
                      rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)-v* &
                           & gll%weigh(v2)*dist(modulo(x1*ng+1-1,nx*ng)+1,(v1-1)*ng+v2)/ &
                           & mesh%jac(nx+1,v1)
                   end if
                else if (x2==1) then
                   if (mesh%etat2(v1)+(1.0d0+gll%node(v2))/mesh%jac(nx+1,v1) >= 0.0d0) then
                      rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+v* &
                           & gll%weigh(v2)*dist(modulo((x1-1)*ng-1,nx*ng)+1,(v1-1)*ng+v2)/ &
                           & mesh%jac(nx+1,v1) 
                   else
                      rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                           & v*gll%weigh(v2)*dist((x1-1)*ng+1,(v1-1)*ng+v2)/mesh%jac(nx+1,v1) 
                   end if
                end if

                if (v2==ng) then
                   if (field_e((x1-1)*ng+x2) >= 0.0d0) then
                      if (v1<nv) then
                         rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                              & gll%weigh(x2)*field_e((x1-1)*ng+x2)* &
                              & dist((x1-1)*ng+x2,v1*ng+1)/mesh%jac(x1,nv+1) 
                      end if
                   else
                      rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                           & gll%weigh(x2)*field_e((x1-1)*ng+x2)* &
                           & dist((x1-1)*ng+x2,v1*ng)/mesh%jac(x1,nv+1) 
                   end if
                else if (v2==1) then
                   if (field_e((x1-1)*ng+x2) >= 0.0d0) then
                      rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                           & gll%weigh(x2)*field_e((x1-1)*ng+x2)* &
                           & dist((x1-1)*ng+x2,(v1-1)*ng+1)/mesh%jac(x1,nv+1)
                   else 
                      if (v1>=2) then
                         rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                              & gll%weigh(x2)*field_e((x1-1)*ng+x2)* &
                              & dist((x1-1)*ng+x2,(v1-1)*ng)/mesh%jac(x1,nv+1)
                      end if
                   end if
                end if

                !!!>>only to test Blanca's case
                rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                     & exp(-0.25d0*(4.0d0*v-1.0d0)**2)*(((4.0d0*sqrt(sll_pi)+2.0d0)*v- &
                     & (2.0d0*sll_pi+sqrt(sll_pi)))*sin(2.0d0*x-2.0d0*sll_pi*t)+ &
                     & sqrt(sll_pi)*(0.25d0-v)*sin(4.0d0*x-4.0d0*sll_pi*t))* &
                     & gll%node(x2)*gll%node(v2)/mesh%jac(x1,v1)
                !!!<<

                rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)*mesh%jac(x1,v1) &
                     & /(gll%weigh(x2)*gll%weigh(v2))

             end do
          end do
       end do
    end do

  end subroutine rhs4dg_1d

  subroutine rk4_4dg_1d(plan,gll,mesh,dist,distp1)
    !< @brief Computation of RK4 steps for Vlasov-Poisson with DG, returns the distribution
    !!        at time n+1
    !! @details Computation of RK4 steps for Vlasov-Poisson with DG, returns the distribution
    !!        at time n+1. This is the routine to call at each time step. Also see 
    !!        initialization of dg_time_steping object for more details \n
    !!        Be aware that this code does not check the size of objects. It might gives you
    !!        segmentation faults if there are error in object size.
    !! @param[INOUT] plan dg_time_steping object, contains work array and data for RK4
    !! @param[IN] gll gausslobatto1D object, necessary for the computation
    !! @param[IN] mesh non_unif_cart_mesh object, necessary for the computation
    !! @param[IN] dist distribution function at time n, 2D real array
    !! @param[OUT] distp1 distribution function at time n+1, 2D real array, same size as dist

    implicit none
    
    type(dg_time_steping),intent(inout) :: plan
    type(gausslobatto1D),intent(in) :: gll
    type(non_unif_cart_mesh),intent(in) :: mesh
    sll_real64,dimension(:,:),intent(in) :: dist
    sll_real64,dimension(:,:),intent(out) :: distp1

    !first step of RK4
    plan%k(2,:,:)=dist
    call rk4dg_step(plan,gll,mesh)
    distp1=dist+plan%k(1,:,:)/6.0d0*plan%dt

    !second step of RK4
    plan%k(2,:,:)=dist+plan%k(1,:,:)/2.0d0*plan%dt
    call rk4dg_step(plan,gll,mesh)
    distp1=distp1+plan%k(1,:,:)/3.0d0*plan%dt

    !third step of RK4
    plan%k(2,:,:)=dist+plan%k(1,:,:)/2.0d0*plan%dt
    call rk4dg_step(plan,gll,mesh)
    distp1=distp1+plan%k(1,:,:)/3.0d0

    !fourth step of RK4
    plan%k(2,:,:)=dist+plan%k(1,:,:)*plan%dt
    call rk4dg_step(plan,gll,mesh)
    distp1=distp1+plan%k(1,:,:)/6.0d0

  end subroutine rk4_4dg_1d

  subroutine rk4dg_step(plan,gll,mesh)

    implicit none
    
    type(dg_time_steping),intent(inout) :: plan
    type(gausslobatto1D),intent(in) :: gll
    type(non_unif_cart_mesh),intent(in) :: mesh

    sll_int32 :: x1,x2,v1,v2,nx,nv,ng,zero

    !!!>>only to test Blanca's case 
    sll_real64 :: x,t

    t=plan%t
    !!!<<

    nx=mesh%n_etat1
    nv=mesh%n_etat2
    ng=gll%degree+1

    zero=1

    plan%rho=0.0d0
    do x1=1,nx
       do x2=1,ng
          do v1=1,nv
             do v2=1,ng
                plan%rho((x1-1)*ng+x2)=plan%rho((x1-1)*ng+x2)+ &
                     & plan%k(2,(x1-1)*ng+x2,(v1-1)*ng+v2)*gll%weigh(v2)/mesh%jac(nx+1,v1)
             end do
          end do
          plan%rho((x1-1)*ng+x2)=plan%rho((x1-1)*ng+x2)*gll%weigh(x2)*mesh%jac(x1,nv+1)
       end do
    end do
!!$    do i=1,nx*ng
!!$       do j=1,nv
!!$          do k=1,ng
!!$             plan%rho(i)=plan%rho(i)+plan%k(2,i,(j-1)*ng+k)*gll%weigh(k)/mesh%jac(nx+1,j)
!!$          end do
!!$       end do
!!$    end do
!!$    !rho=M*rho
!!$    do i=1,nx
!!$       do j=1,ng
!!$          plan%rho((i-1)*ng+j)=plan%rho((i-1)*ng+j)*gll%weigh(j)*mesh%jac(i,nv+1)
!!$       end do
!!$    end do
    plan%rho(plan%x0)=0.0d0 ! boudary condition = \Phi(0,t) = 0

!!$    !we first consider v<0
!!$    if (mesh%etat2(1)<0.0d0) then
!!$       !v<0 for some nodes
!!$       !solve the Poisson problem
!!$       call poisson_solve_4dg_1d(plan%umfpack_data,plan%matvm,plan%rho,plan%phi)
!!$       plan%field=matmul(plan%fieldvm,plan%phi)
!!$       !construction or rhs
!!$       !rhs=-v.d_x(f)+d_x(Phi).d_v(f), with d=\partial
!!$       call rhs4dg_1d(mesh,gll,plan%field,plan%k(2,:,:),plan%rhs)
!!$
!!$       if (mesh%etat2(nv)>=0) then
!!$          do v1=1,nv
!!$             if (abs(mesh%etat2(v1))<=epsilon( max(abs(mesh%etat2(1)), & 
!!$                  & abs(mesh%etat2(mesh%n_etat2+1))) )) then
!!$                zero=v1
!!$             end if
!!$             plan%k(1,:,1:(zero-1)*ng)=plan%rhs(:,1:(zero-1)*ng)
!!$          end do
!!$       else
!!$          plan%k(1,:,:)=plan%rhs
!!$       end if
!!$    end if
!!$    
!!$    !then we consider v>=0
!!$    if (mesh%etat2(nv)>=0.0d0) then
!!$       !v>=0 for some nodes
!!$       !solve the Poisson problem
!!$       call poisson_solve_4dg_1d(plan%umfpack_data,plan%matvp,plan%rho,plan%phi)
!!$       plan%field=matmul(plan%fieldvp,plan%phi)
!!$       !construction or rhs
!!$       !rhs=-v.d_x(f)+d_x(Phi).d_v(f), with d=\partial
!!$       call rhs4dg_1d(mesh,gll,plan%field,plan%k(2,:,:),plan%rhs)
!!$
!!$        plan%k(1,:,(zero-1)*ng+1:nv*ng)=plan%rhs(:,(zero-1)*ng+1:nv*ng)
!!$    end if

    !!!>>>only to test Blanca's case, we set the field manually
    do x1=1,nx
       do x2=1,ng
          x=mesh%etat1(x1)+(1.0d0+gll%node(x2))/mesh%jac(x1,nv+1)
          plan%field((x1-1)*ng+x2)=sqrt(sll_pi)/4.0d0*sin(2.0d0*x-2.0d0*t)
       end do
    end do
    call rhs4dg_1d(mesh,gll,plan%field,plan%k(2,:,:),plan%rhs,t)
    plan%k(1,:,:)=plan%rhs
    !!!<<<

  end subroutine rk4dg_step

end module timestep4dg
