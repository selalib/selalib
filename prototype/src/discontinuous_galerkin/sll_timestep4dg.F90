!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: poisson4dg
!
! DESCRIPTION:
!> @file timestep4dg.F90
!! @author Madaule Eric
!! @brief Time step tools for discontinous Galerkin
!! @details This module contains the time steping for the Vlasov equation using 
!!          discontinuous Galerkin. The initialization is done with poisson4dg (should 
!!          be transfert to Poisson library, currently in file vlasov_poisson_DG.F90).\n
!!          We chose to have a good readability of the code rather than high performance.\n
!!
!!          In this module, I refere to Blanca Ayuso and her article Discontinuous
!!          Galerkin methods for the one-dimensional Vlasov-Poisson system.\n
!!
!!          This module will first be limited to 1D and should extend as people will
!!          have the need for higher dimension (and time to write it).
!!         
!------------------------------------------------------------------------------
module sll_timestep4dg
#include "sll_working_precision.h"
#include "sll_boundary_condition_descriptors.h"

  use sll_gausslobatto
  use sll_poisson4dg_1d_periodic_dg
  use sll_nu_cart_mesh
  use sll_constants
  !those are part of FEMilaro
  use mod_sparse
  use mod_umfpack

  !use mod_octave_io_sparse

  implicit none

  logical,parameter :: x_upwind=.true.,v_upwind=.true.
  integer,parameter :: bc_v=sll_periodic ! sll_dirichlet
  integer,parameter :: SLL_RK3=0,SLL_RK4=1,SLL_EE=2,SLL_TVDRK2=3,sll_SSP_RK=4
  !<parameters to compute the time integration

  type dg_time_steping
     !< @brief Plan for time stepping with discontinuous Galerkin.
     !! @details Plan for time stepping with discontinuous Galerkin.
     !!          It contains a working array, the alectrostatic and electric fields, the charge
     !!          distribution, the right hand side of Vlasov equation, the time step, the node
     !!          to apply \Phi(x=0)=alpha and array for UMFpack so they are not created every
     !!          time. Must be build with init_timesteping_4dg (see corresponding section) and
     !!          deleted with clear.
     sll_real64,dimension(:,:,:),allocatable :: k ! working array
     sll_real64,dimension(:),allocatable :: phi,field,rho
     sll_real64,dimension(:,:),allocatable :: rhs ! right hand side of Vlasov
     sll_real64 :: dt,t,bound,norma ! time step and time, value for Vlasov-Poisson problem
     sll_int32 :: x0 ! place to apply phi(x0)=alpha
     sll_int32 :: method,zero ! computation method, node such as j<zero => v<0
     ! Poisson's problem
     type(poisson_1d_periodic_dg) :: poisson_vp,poisson_vm
     !for SSP RK + maximumu/minimum limiters
     sll_real64 :: max0,min0
  end type dg_time_steping

  interface new
     module procedure init_timesteping_4dg
  end interface new

  interface delete
     module procedure delete_dg_step
  end interface delete

contains

  subroutine init_timesteping_4dg(plan,method,gll,mesh,dt,xbound,c11,c12,norma,alpha)
    !< @brief This routine is the interface to use to build interpolators for DG
    !! @details This routine only is the interface to use to build interpolators for DG.
    !!          It will extend as time stepping method will be added\n
    !!          Do not forget to considere the Courant number when choosing dt\n
    !!          Courant number : |v|dt/h<C/k, |v|=max(|E|)+max(|v|) for Vlasov
    !! @param[OUT] plan dg_time_steping object to build
    !! @param[IN] method SLL_RK3 or SLL_RK4 to choose your interpolator
    !! @param[IN] gll gausslobatto1d object
    !! @param[IN] mesh non_unif_cart_mesh object, needed for various calculation
    !! @param[IN] dt time step, real double precision
    !! @param[IN] xbound coordinate to set \Phi(x=0)=alpha, integer, see Poisson4DG
    !! @param[IN] c11 coefficient c11 for flux
    !! @param[IN] c12 coefficient c12(v>0) for flux with assumption that c12(v>0)=-c12(v<0)
    !! @param[IN,OPTIONAL] norma coefficient such that that the Poisson's problem is
    !!                           -Lapalacian(Phi)=rho-norma
    !! @param[IN,OPTIONAL] alpha coefficient to set such that Phi(x0)=alpha, default value is 0

    implicit none

    type(dg_time_steping),intent(out) :: plan
    type(gausslobatto1d),intent(in) :: gll
    type(non_unif_cart_mesh),intent(in) :: mesh
    sll_int32,intent(in) :: method,xbound
    sll_real64,intent(in) :: dt,c11,c12
    sll_real64,intent(in),optional :: norma
    sll_real64,intent(in),optional :: alpha

    sll_real64 :: norma0,alpha0

    plan%method=method

    if (present(alpha)) then
       alpha0=alpha
    else
       alpha0=0.0d0
    end if
    if (present(norma)) then
       norma0=norma
    else
       norma0=0.0d0
    end if

    !due to evolution in the code, rk4, ee and tvd_rk2 have the same plan

    if (plan%method==SLL_RK3) then
       !RK3
       print*,'time integration using RK3'
       call init_k3_4dg(plan,gll,mesh,dt,xbound,c11,c12,norma=norma0,alpha=alpha0)
    else if (plan%method==SLL_RK4) then
       !RK4
       print*,'time integration using RK4'
       call init_k2_4dg(plan,gll,mesh,dt,xbound,c11,c12,norma=norma0,alpha=alpha0)
    else if (plan%method==SLL_EE) then
       !Euler explicit
       print*,'time integration using Euler explicit'
       call init_k2_4dg(plan,gll,mesh,dt,xbound,c11,c12,norma=norma0,alpha=alpha0)
    else if (plan%method==SLL_TVDRK2) then
       !TVD RK2
       print*,'time integration using TVD RK2'
       call init_k2_4dg(plan,gll,mesh,dt,xbound,c11,c12,norma=norma0,alpha=alpha0)
    else if (plan%method==sll_ssp_rk) then
       !SSP RK (5,4)
       print*,"time integration using SSP RK (5,4)"
       call init_k2_4dg(plan,gll,mesh,dt,xbound,c11,c12,norma=norma0,alpha=alpha0)
    end if

  end subroutine init_timesteping_4dg

  subroutine init_k2_4dg(plan,gll,mesh,dt,xbound,c11,c12,norma,alpha)
    !< @brief Do not use it, call the interface routine init_timesteping_4dg

    implicit none

    type(dg_time_steping),intent(inout) :: plan
    type(gausslobatto1d),intent(in) :: gll
    type(non_unif_cart_mesh),intent(in) :: mesh
    sll_int32,intent(in) :: xbound
    sll_real64,intent(in) :: dt,c11,c12
    sll_real64,intent(in),optional :: norma
    sll_real64,intent(in),optional :: alpha

    sll_int32 :: nx,nv,ng,v1,v2
    sll_real64 :: vv1,vv2

    nx=mesh%n_etat1
    nv=mesh%n_etat2
    ng=gll%degree+1

    allocate(plan%k(nx*ng,nv*ng,2),plan%phi(nx*ng),plan%rho(nx*ng),plan%field(nx*ng), &
         & plan%rhs(nx*ng,nv*ng))
    plan%dt=dt
    plan%x0=xbound
    if (present(alpha)) then
       plan%bound=alpha
    else
       plan%bound=0.0d0
    end if
    if (present(norma)) then
       plan%norma=norma
    else
       plan%norma=0.0d0
    end if

    plan%zero=0
    if ( mesh%etat2(nv+1)<=epsilon(1.0d0)*abs(mesh%etat2(1))) then
       plan%zero=nv*ng+1
    else if(mesh%etat2(1)<epsilon(1.0d0)*abs(mesh%etat2(nv+1)) .and. & 
         & mesh%etat2(nv+1)>epsilon(1.0d0)*abs(mesh%etat2(1))) then
       vv1=mesh%etat2(1)
       do v1=1,nv
          do v2=1,ng
             vv2=vv1
             vv1=mesh%etat2(v1)+(1.0d0+gll%node(v2))/mesh%jac(nx+1,v1)
             if (vv2<0.0d0 .and. vv1>=0.0d0 .and. plan%zero==0) then
                plan%zero=(v1-1)*ng+v2+1
             end if
          end do
       end do
    else if (mesh%etat2(1)>=epsilon(1.0d0)*abs(mesh%etat2(nv+1))) then
       plan%zero=1
    else
       print*,'error in the source code'
       print*,'can not tell whether 0 is in velocity domain'
       plan%zero=-huge(nv)
    end if

    ! initialization of Poisson's
    call new(plan%poisson_vp,gll,nx,mesh%jac(1:nx,nv+1),c11,c12,xbound,.true.)
    call new(plan%poisson_vm,gll,nx,mesh%jac(1:nx,nv+1),c11,c12,xbound,.true.)

    plan%max0=-2.0d0
    plan%min0=2.0d0

  end subroutine init_k2_4dg

  subroutine init_k3_4dg(plan,gll,mesh,dt,xbound,c11,c12,norma,alpha)
    !< @brief Do not use it, call the interface routine init_timesteping_4dg

    implicit none

    type(dg_time_steping),intent(inout) :: plan
    type(gausslobatto1d),intent(in) :: gll
    type(non_unif_cart_mesh),intent(in) :: mesh
    sll_int32,intent(in) :: xbound
    sll_real64,intent(in) :: dt,c11,c12
    sll_real64,intent(in),optional :: norma
    sll_real64,intent(in),optional :: alpha

    sll_int32 :: nx,nv,ng,v1,v2
    sll_real64 :: vv1,vv2

    nx=mesh%n_etat1
    nv=mesh%n_etat2
    ng=gll%degree+1

    allocate(plan%k(nx*ng,nv*ng,3),plan%phi(nx*ng),plan%rho(nx*ng),plan%field(nx*ng), &
         & plan%rhs(nx*ng,nv*ng))
    plan%dt=dt
    plan%x0=xbound
    if (present(alpha)) then
       plan%bound=alpha
    else
       plan%bound=0.0d0
    end if
    if (present(norma)) then
       plan%norma=norma
    else
       plan%norma=0.0d0
    end if

    plan%zero=0
    if ( mesh%etat2(nv+1)<=epsilon(1.0d0)*abs(mesh%etat2(1))) then
       plan%zero=nv*ng+1
    else if(mesh%etat2(1)<epsilon(1.0d0)*abs(mesh%etat2(nv+1)) .and. & 
         & mesh%etat2(nv+1)>epsilon(1.0d0)*abs(mesh%etat2(1))) then
       vv1=mesh%etat2(1)
       do v1=1,nv
          do v2=1,ng
             vv2=vv1
             vv1=mesh%etat2(v1)+(1.0d0+gll%node(v2))/mesh%jac(nx+1,v1)
             if (vv2<0.0d0 .and. vv1>=0.0d0 .and. plan%zero==0) then
                plan%zero=(v1-1)*ng+v2+1
             end if
          end do
       end do
    else if (mesh%etat2(1)>=epsilon(1.0d0)*abs(mesh%etat2(nv+1))) then
       plan%zero=1
    else
       print*,'error in the source code'
       print*,'can not tell whether 0 is in velocity domain'
       plan%zero=-huge(nv)
    end if

    ! initialization of Poisson's
    call new(plan%poisson_vp,gll,nx,mesh%jac(1:nx,nv+1),c11,c12,xbound,.true.)
    call new(plan%poisson_vm,gll,nx,mesh%jac(1:nx,nv+1),c11,c12,xbound,.true.)

    plan%max0=-2.0d0
    plan%min0=2.0d0

  end subroutine init_k3_4dg

  subroutine delete_dg_step(plan)
    !< @brief deletion of dg_time_steping object
    !! @details Deletion  of dg_time_steping object. It deallocate and clear all working 
    !!         array inside.
    !! @param[INOUT] plan dg_time_steping object to clear

    implicit none

    type(dg_time_steping),intent(inout) :: plan

    deallocate(plan%k,plan%field,plan%phi,plan%rho,plan%rhs)
    call delete(plan%poisson_vp)
    call delete(plan%poisson_vm)
    
  end subroutine delete_dg_step

  subroutine time_integration_4dg(plan,gll,mesh,dist,distp1)
    !< @brief Interface for the time integration for DG.
    !! @details Interface for the time integration for DG. It will call the routine that 
    !!          correspond to the time itegrator you choosed at initialization of plan. Also see 
    !!          initialization of dg_time_steping object for more details \n
    !!          Be aware that this code does not check the size of objects. It might gives you
    !!          segmentation faults if there are error in object size.
    !! @param[INOUT] plan dg_time_steping object, contains work array and data
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

    if (plan%method==SLL_RK3) then
       !RK3
       call rk3_4dg_1d(plan,gll,mesh,dist,distp1)
    else if (plan%method==SLL_RK4) then
       !RK4
       call rk4_4dg_1d(plan,gll,mesh,dist,distp1)
    else if (plan%method==SLL_EE) then
       !ee
       call ee4dg_1d(plan,gll,mesh,dist,distp1)
    else if (plan%method==SLL_TVDRK2) then
       !TVD RK2
       call tvdrk2_4dg_1d(plan,gll,mesh,dist,distp1)
    else if (plan%method==sll_ssp_rk) then
       !SSP RK (5,4)
       call ssp_rk_4dg_1d(plan,gll,mesh,dist,distp1)
    else
       print*,'no time integration chosen'
       print*,'exiting...'; stop
    end if

  end subroutine time_integration_4dg

  subroutine rk4_4dg_1d(plan,gll,mesh,dist,distp1)
    !< @brief Computation of RK4 steps for Vlasov-Poisson with DG, returns the distribution
    !!        at time n+1. \n
    !!        You shouldn't call it but call time_integration_4dg which is the interface
    !! @details Computation of RK4 steps for Vlasov-Poisson with DG, returns the distribution
    !!          at time n+1. This routine is called at each time step. Also see 
    !!          initialization of dg_time_steping object and time_integration_4dg for more 
    !!          details. \n
    !!          Be aware that this code does not check the size of objects. It might gives you
    !!          segmentation faults if there are error in object size.
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
    plan%k(:,:,2)=dist
    call dg_step_time_integration(plan,gll,mesh)
    distp1=dist+plan%k(:,:,1)/6.0d0*plan%dt

    !second step of RK4
    plan%k(:,:,2)=dist+plan%k(:,:,1)/2.0d0*plan%dt
    plan%t=plan%t+plan%dt/2.0d0
    call dg_step_time_integration(plan,gll,mesh)
    distp1=distp1+plan%k(:,:,1)/3.0d0*plan%dt

    !third step of RK4
    plan%k(:,:,2)=dist+plan%k(:,:,1)/2.0d0*plan%dt
    call dg_step_time_integration(plan,gll,mesh)
    distp1=distp1+plan%k(:,:,1)/3.0d0*plan%dt

    !fourth step of RK4
    plan%k(:,:,2)=dist+plan%k(:,:,1)*plan%dt
    plan%t=plan%t+plan%dt/2.0d0
    call dg_step_time_integration(plan,gll,mesh)
    distp1=distp1+plan%k(:,:,1)/6.0d0*plan%dt

  end subroutine rk4_4dg_1d

  subroutine ee4dg_1d(plan,gll,mesh,dist,distp1)
    !< @brief Computation of explicit Euler steps for Vlasov-Poisson with DG, returns the
    !!        distribution at time n+1. \n
    !!        You shouldn't call it but call time_integration_4dg which is the interface
    !! @details Computation of explicit Eule steps for Vlasov-Poisson with DG, returns the 
    !!          distribution at time n+1. This routine is called at each time step. Also see 
    !!          initialization of dg_time_steping object and time_integration_4dg for more 
    !!          details. \n
    !!          Be aware that this code does not check the size of objects. It might gives you
    !!          segmentation faults if there are error in object size.
    !! @param[INOUT] plan dg_time_steping object, contains work array and data for the scheme
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

    plan%k(:,:,2)=dist
    call dg_step_time_integration(plan,gll,mesh)
    distp1=dist+plan%k(:,:,1)*plan%dt

  end subroutine ee4dg_1d

  subroutine tvdrk2_4dg_1d(plan,gll,mesh,dist,distp1)
    !< @brief Computation of TVD RK2 steps for Vlasov-Poisson with DG, returns the
    !!        distribution at time n+1. \n
    !!        You shouldn't call it but call time_integration_4dg which is the interface
    !! @details Computation of TVD RK2 steps for Vlasov-Poisson with DG, returns the 
    !!          distribution at time n+1. This routine is called at each time step. Also see 
    !!          initialization of dg_time_steping object and time_integration_4dg for more 
    !!          details. \n
    !!          Be aware that this code does not check the size of objects. It might gives you
    !!          segmentation faults if there are error in object size.
    !! @param[INOUT] plan dg_time_steping object, contains work array and data for the scheme
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

    !first step of TVD RK2
    plan%k(:,:,2)=dist
    call dg_step_time_integration(plan,gll,mesh)
    plan%k(:,:,2)=plan%dt*plan%k(:,:,1)+dist

    !second step of TVD RK2
    plan%t=plan%t+plan%dt
    call dg_step_time_integration(plan,gll,mesh)
    distp1=0.5d0*dist+0.5d0*plan%k(:,:,2)+0.5d0*plan%dt*plan%k(:,:,1)

  end subroutine tvdrk2_4dg_1d

  subroutine rk3_4dg_1d(plan,gll,mesh,dist,distp1)
    !< @brief Computation of RK3 steps for Vlasov-Poisson with DG, returns the
    !!        distribution at time n+1. \n
    !!        You shouldn't call it but call time_integration_4dg which is the interface
    !! @details Computation of TVD RK3 steps for Vlasov-Poisson with DG, returns the 
    !!          distribution at time n+1. This routine is called at each time step. Also see 
    !!          initialization of dg_time_steping object and time_integration_4dg for more 
    !!          details. \n
    !!          Be aware that this code does not check the size of objects. It might gives you
    !!          segmentation faults if there are error in object size.
    !! @param[INOUT] plan dg_time_steping object, contains work array and data for the scheme
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

    !first step of RK3
    plan%k(:,:,2)=dist
    call dg_step_time_integration(plan,gll,mesh)
    plan%k(:,:,3)=plan%k(:,:,1)
    plan%k(:,:,2)=dist+plan%k(:,:,1)*plan%dt/2.0d0
    distp1=dist+plan%k(:,:,1)*plan%dt/6.0d0

    !second step of RK3
    plan%t=plan%t+plan%dt/2.0d0
    call dg_step_time_integration(plan,gll,mesh)
    distp1=distp1+plan%k(:,:,1)*plan%dt*2.0d0/3.0d0
    plan%k(:,:,2)=dist+plan%dt*(-plan%k(:,:,3)+2.0d0*plan%k(:,:,1))

    !third step of RK3
    plan%t=plan%t+plan%dt/2.0d0
    call dg_step_time_integration(plan,gll,mesh)
    distp1=distp1+plan%k(:,:,1)*plan%dt/6.0d0

  end subroutine rk3_4dg_1d

  subroutine ssp_rk_4dg_1d(plan,gll,mesh,dist,distp1)
    !< @brief Computation of SSP RK(5,4) steps for Vlasov-Poisson with DG, returns the
    !!        distribution at time n+1. \n
    !!        You shouldn't call it but call time_integration_4dg which is the interface
    !! @details Computation of SSP RK(5,4) steps for Vlasov-Poisson with DG, returns the 
    !!          distribution at time n+1. This routine is called at each time step. Also see 
    !!          initialization of dg_time_steping object and time_integration_4dg for more 
    !!          details. The coefficit for SSP are taken in "On High Order Strong Stability
    !!          Preserving Runge–Kutta and Multi Step Time Discretizations" by Sigal Gottlieb \n
    !!          Be aware that this code does not check the size of objects. It might gives you
    !!          segmentation faults if there are error in object size.
    !! @param[INOUT] plan dg_time_steping object, contains work array and data for the scheme
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

    !first step
    plan%k(:,:,2)=dist
    call dg_step_time_integration(plan,gll,mesh)
    plan%k(:,:,2)=dist+0.391752226571890d0*plan%dt*plan%k(:,:,1)
    call maximum_minimum_limiters(gll,mesh,plan%min0,plan%max0,plan%k(:,:,2))

    !2nd step
    call dg_step_time_integration(plan,gll,mesh)
    plan%k(:,:,2)=0.444370493651235d0*dist+0.555629506348765d0*plan%k(:,:,2)+ &
         & 0.368410593050371d0*plan%dt*plan%k(:,:,1)
    call maximum_minimum_limiters(gll,mesh,plan%min0,plan%max0,plan%k(:,:,2))
    distp1=0.517231671970585d0*plan%k(:,:,2)

    !3rd step
    call dg_step_time_integration(plan,gll,mesh)
    plan%k(:,:,2)=0.620101851488403d0*dist+0.379898148511597d0*plan%k(:,:,2)+ &
         & 0.251891774271694d0*plan%dt*plan%k(:,:,1)
    call maximum_minimum_limiters(gll,mesh,plan%min0,plan%max0,plan%k(:,:,2))
    distp1=distp1+0.096059710526147d0*plan%k(:,:,2)

    !4th step
    call dg_step_time_integration(plan,gll,mesh)
    plan%k(:,:,2)=0.178079954393132d0*dist+0.821920045606868d0*plan%k(:,:,2)+ &
         & 0.544974750228521d0*plan%dt*plan%k(:,:,1)
    call maximum_minimum_limiters(gll,mesh,plan%min0,plan%max0,plan%k(:,:,2))
    distp1=distp1+0.063692468666290d0*plan%dt*plan%k(:,:,1)+0.386708617503269d0*plan%k(:,:,2)

    !5th step
    !most of the 5th step has already been done at the end of step 2 to 4
    !when updating distp1
    call dg_step_time_integration(plan,gll,mesh)
    distp1=distp1+0.226007483236906d0*plan%dt*plan%k(:,:,1)
    call maximum_minimum_limiters(gll,mesh,plan%min0,plan%max0,distp1)
    
  end subroutine ssp_rk_4dg_1d

  subroutine maximum_minimum_limiters(gll,mesh,min0,max0,dist)
    !< @brief Averaging of  a 2D function to reduce ocsillations according to the maximum-minimum
    !!        conservation.
    !! @details Averaging of  a 2D function to reduce ocsillations according to the maximum-
    !!          minimum conservation. This routine is called in the SSP-RK at each stage. The
    !!          idea comes from "On maximum-principle-satisfying high order schemes for 
    !!          scalar conservation laws" by Xiangxiong Zhang and Chi-Wang Shu.
    !! @param[IN] gll gausslobatto1D object, necessary for the computation
    !! @param[IN] mesh non_unif_cart_mesh object, necessary for the computation
    !! @param[IN] min real, minimum value allowed to the distribution
    !! @param[IN] max real, maximum value allowed to the distribution, if max <= min 
    !!                nothing will be done
    !! @param[INOUT] dist 2d real array to rescale

    !we compute p°(x)=theta(p(x)-p_mean)+p_mean, 
    !theta=min{|(max-p_mean)/(maxloc-p_mean)|,|(min-p_mean)/(minloc-p_mean)|,1}

    implicit none

    type(gausslobatto1D),intent(in) :: gll
    type(non_unif_cart_mesh),intent(in) :: mesh
    sll_real64,intent(in) :: min0,max0
    sll_real64,dimension(:,:),intent(inout) :: dist

    sll_int32 :: v1,v2,x1,x2,ng
    sll_real64 :: mean,max_loc,min_loc,theta

    ng=gll%degree+1

    if (min0<max0) then
    !if (.false.) then
       do v1=1,mesh%n_etat2
          do x1=1,mesh%n_etat1
             max_loc=min0
             min_loc=max0
             mean=0.0d0
             do v2=1,ng
                do x2=1,ng
                   min_loc=min(min_loc,dist((x1-1)*ng+x2,(v1-1)*ng+v2))
                   max_loc=max(max_loc,dist((x1-1)*ng+x2,(v1-1)*ng+v2))
                   mean=mean+gll%weigh(x2)*gll%weigh(v2)*dist((x1-1)*ng+x2,(v1-1)*ng+v2)
                end do
             end do

             mean=mean/4.0d0

             if (abs(max_loc-mean)>epsilon(1.0d0)*max_loc .and. &
                  & abs(min_loc-mean)>epsilon(1.0d0)*min_loc) then
                theta=min(abs((max0-mean)/(max_loc-mean)),abs((min0-mean)/(min_loc-mean)),1.0d0)
             else if(abs(max_loc-mean)>epsilon(1.0d0)*max_loc .and. &
                  & abs(min_loc-mean)<=epsilon(1.0d0)*min_loc) then
                theta=min(abs((max0-mean)/(max_loc-mean)),1.0d0)
                theta=min(abs((max0-mean)/(max_loc-mean)),abs((min0-mean)/(min_loc-mean)),1.0d0)
             else if(abs(max_loc-mean)<=epsilon(1.0d0)*max_loc .and. &
                  & abs(min_loc-mean)>epsilon(1.0d0)*min_loc) then
                theta=min(abs((min0-mean)/(min_loc-mean)),1.0d0)
             else
                theta=1.0d0
             end if

             do v2=1,gll%degree+1
                do x2=1,gll%degree+1
                   dist((x1-1)*ng+x2,(v1-1)*ng+v2)= &
                        & theta*(dist((x1-1)*ng+x2,(v1-1)*ng+v2)-mean) + mean
                end do
             end do

          end do
       end do
    end if

  end subroutine maximum_minimum_limiters

  subroutine dg_step_time_integration(plan,gll,mesh)
    !< @brief Computation for time integration with DG. You do not need to call it, it is used
    !!        as it should in time integration routines.
    !! @details Computation for time integration with DG. You do not need to call it unless 
    !!          you want to use a time integration method which is not already implmented. 
    !!          It is used as it should in time integration routines. In the case you would like 
    !!          to use it directly, please chek the source code of this routine and the source
    !!          code of routine ee4dg_1d (Euler explicit, which is probably easier to understand
    !!          than RK4).

    implicit none
    
    type(dg_time_steping),intent(inout) :: plan
    type(gausslobatto1D),intent(in) :: gll
    type(non_unif_cart_mesh),intent(in) :: mesh

    sll_int32 :: x1,x2,v1,v2,nx,nv,ng,zero

!!$    !!!>>>only to test Blanca's case 
!!$    sll_real64 :: x
!!$    !!!<<<

    nx=mesh%n_etat1
    nv=mesh%n_etat2
    ng=gll%degree+1
    zero=plan%zero

    plan%rho=0.0d0
    plan%norma=0.0d0
    plan%k(:,:,1)=0.0d0
    do x1=1,nx
       do x2=1,ng
          do v1=1,nv
             do v2=1,ng
                plan%rho((x1-1)*ng+x2)=plan%rho((x1-1)*ng+x2)+ &
                     & plan%k((x1-1)*ng+x2,(v1-1)*ng+v2,2)*gll%weigh(v2)/mesh%jac(nx+1,v1)
             end do
          end do
          plan%norma=plan%norma+plan%rho((x1-1)*ng+x2)*gll%weigh(x2)/mesh%jac(x1,nv+1)
       end do
    end do
    plan%norma=plan%norma/(mesh%etat1(nx+1)-mesh%etat1(1))
    do x1=1,nx
       do x2=1,ng
          plan%rho((x1-1)*ng+x2)=(plan%rho((x1-1)*ng+x2)-plan%norma)* &
               & gll%weigh(x2)/mesh%jac(x1,nv+1)
       end do
    end do
    plan%rho(plan%x0)=plan%bound ! boudary condition = \Phi(x0,t) = alpha !we first consider v<0

    !we first consider v<0
    if (mesh%etat2(1)<0.0d0) then
       !v<0 for some nodes
       !solve the Poisson problem
       call solve(plan%poisson_vm,plan%rho,plan%phi)
       plan%field=matmul(plan%poisson_vm%mat_field,plan%phi)
       !construction of rhs
       !rhs=-v.d_x(f)+d_x(Phi).d_v(f), with d=\partial
       call rhs4dg_1d(mesh,gll,plan%field,plan%k(:,:,2),plan%rhs)

       plan%k(:,1:zero-1,1)=plan%rhs(:,1:zero-1)
    end if
    
    !then we consider v>=0
    if (mesh%etat2(nv)>=0.0d0) then
       !v>=0 for some nodes
       !solve the Poisson problem
       call solve(plan%poisson_vp,plan%rho,plan%phi)
       plan%field=matmul(plan%poisson_vp%mat_field,plan%phi)
       !construction or rhs
       !rhs=-v.d_x(f)+d_x(Phi).d_v(f), with d=\partial
       call rhs4dg_1d(mesh,gll,plan%field,plan%k(:,:,2),plan%rhs)

        plan%k(:,zero:nv*ng,1)=plan%rhs(:,zero:nv*ng)
    end if

    !!!>>>only to test Blanca's case, we set the field manually
!!$    do x1=1,nx
!!$       do x2=1,ng
!!$          x=mesh%etat1(x1)+(1.0d0+gll%node(x2))/mesh%jac(x1,nv+1)
!!$          plan%field((x1-1)*ng+x2)=sqrt(sll_pi)/4.0d0*sin(2.0d0*x-2.0d0*sll_pi*plan%t)
!!$       end do
!!$    end do
!!$   do x1=1,nx
!!$       do x2=1,ng
!!$          x=mesh%etat1(x1)+(1.0d0+gll%node(x2))/mesh%jac(x1,nv+1)
!!$          plan%field((x1-1)*ng+x2)=sqrt(sll_pi)/4.0d0*sin(2.0d0*x)*10.0d0
!!$       end do
!!$    end do
!!$    call rhs4dg_1d(mesh,gll,plan%field,plan%k(:,:,2),plan%rhs,plan%t)
!!$    plan%k(:,:,1)=plan%rhs
    !!!<<<

  end subroutine dg_step_time_integration

  subroutine rhs4dg_1d(mesh,gll,field_e,dist,rhs,t)
    !< @brief Computation of rhs for Vlasov equation in 2D phase space (1D physical space)
    !! @details Computation of rhs for Vlasov equation in 2D phase space (1D physical space).
    !!          All variable must already be initialized.\n
    !!          If you want to add a forcing term in Vlasov equation you should add it there.
    !! @param[IN] mesh non_unif_cart_mesh object, gives the Jacobian on each cell and the 
    !!                 number of step
    !! @param[IN] gll gausslobatto1D object, gives most of the needed information for computation
    !! @param[IN] field_e electric field, 1D array of real
    !! @param[IN] dist distribution function, 2D array of real
    !! @param[OUT] rhs right hand side of Vlasov equation, 2D array of real
    !! @param[IN,OPTIONAL] t time, may be needed, especially for Blanca's case

    implicit none

    type(non_unif_cart_mesh),intent(in) :: mesh
    type(gausslobatto1D),intent(in) :: gll
    sll_real64,dimension(:),intent(in) :: field_e
    sll_real64,dimension(:,:),intent(in) :: dist
    sll_real64,dimension(:,:),intent(out) :: rhs
    
    sll_real64,intent(in),optional :: t

    sll_int32 :: nx,nv,ng
    sll_int32 :: x1,x2,v1,v2,i
    sll_real64 :: som1,som2,x,v

    ng=gll%degree+1
    nx=mesh%n_etat1
    nv=mesh%n_etat2

    if (size(field_e)/=nx*ng .or. size(dist,1)/=nx*ng .or. size(dist,2)/=nv*ng) then
       print*,'size error in computation of rhs'
       print*,'please, check your size to find the problem'
       print*,'exiting...'
       stop
    end if

    !construction or rhs
    !rhs=-v.d_x(f)+d_x(Phi).d_v(f) (+forcing term), with d=\partial and d_x(Phi)=E
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
                rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=(som1*v*gll%weigh(v2)- &
                     & som2*field_e((x1-1)*ng+x2)*gll%weigh(x2))/mesh%jac(x1,v1)

                !boudaries part
                !upwind
                if (x_upwind) then
                   if (x2==ng) then 
                      if (mesh%etat2(v1)+(1.0d0+gll%node(v2))/mesh%jac(nx+1,v1) >= 0.0d0) then
                         rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                              & v*dist(x1*ng,(v1-1)*ng+v2)*gll%weigh(v2)/mesh%jac(nx+1,v1)
                      else
                         rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                              & v*dist(modulo(x1*ng+1-1,nx*ng)+1,(v1-1)*ng+v2)* &
                              & gll%weigh(v2)/mesh%jac(nx+1,v1)
                      end if
                   else if (x2==1) then
                      if (mesh%etat2(v1)+(1.0d0+gll%node(v2))/mesh%jac(nx+1,v1) >= 0.0d0) then
                         rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                              & v*dist(modulo((x1-1)*ng-1,nx*ng)+1,(v1-1)*ng+v2)* &
                              & gll%weigh(v2)/mesh%jac(nx+1,v1) 
                      else
                         rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                              & v*dist((x1-1)*ng+1,(v1-1)*ng+v2)*gll%weigh(v2)/mesh%jac(nx+1,v1) 
                      end if
                   end if
                end if

                if (v_upwind) then
                   if (v2==ng) then
                      if (field_e((x1-1)*ng+x2) >= 0.0d0) then
                         if (v1<nv) then
                            rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                                 & field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,v1*ng+1)* &
                                 & gll%weigh(x2)/mesh%jac(x1,nv+1) 
                         else if (bc_v==sll_periodic) then
                            ! activate only for periodicity in v
                            rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                                 & field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,1)* &
                                 & gll%weigh(x2)/mesh%jac(x1,nv+1)
                         end if
                      else
                         rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                              & field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,v1*ng)* &
                              & gll%weigh(x2)/mesh%jac(x1,nv+1) 
                      end if
                   else if (v2==1) then
                      if (field_e((x1-1)*ng+x2) >= 0.0d0) then
                         rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                              & field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,(v1-1)*ng+1)* &
                              & gll%weigh(x2)/mesh%jac(x1,nv+1)
                      else 
                         if (v1>=2) then
                            rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                                 & field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,(v1-1)*ng)* &
                                 & gll%weigh(x2)/mesh%jac(x1,nv+1)
                         else if (bc_v==sll_periodic) then
                            ! activate only for periodicity in v
                            rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                                 & field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,nv*ng)* &
                                 & gll%weigh(x2)/mesh%jac(x1,nv+1)
                         end if
                      end if
                   end if
                end if

                !centered
                if (.not. x_upwind) then
                   if (x2==ng) then 
                      rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                           & 0.5d0*(v*dist(x1*ng,(v1-1)*ng+v2)+ &
                           & v*dist(modulo(x1*ng+1-1,nx*ng)+1,(v1-1)*ng+v2))* &
                           & gll%weigh(v2)/mesh%jac(nx+1,v1)
                   else if (x2==1) then
                      rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                           & 0.5d0*(v*dist(modulo((x1-1)*ng-1,nx*ng)+1,(v1-1)*ng+v2)+ &
                           & v*dist((x1-1)*ng+1,(v1-1)*ng+v2))*gll%weigh(v2)/mesh%jac(nx+1,v1)
                   end if
                end if

                if (.not. v_upwind) then
                   if (v2==ng) then
                      if (v1<nv) then
                         rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                              & 0.5d0*(field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,v1*ng+1)+ &
                              & field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,v1*ng))* &
                              & gll%weigh(x2)/mesh%jac(x1,nv+1)
                      else if (v1==nv) then
                         if (bc_v==sll_periodic) then
                            ! activate only for periodicity in v
                            rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                                 & 0.5d0*(field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,1)+ &
                                 & field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,nv*ng))* &
                                 & gll%weigh(x2)/mesh%jac(x1,nv+1)
                         else
                            rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                                 & 0.5d0*field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,v1*ng)* &
                                 & gll%weigh(x2)/mesh%jac(x1,nv+1)
                         end if
                      end if
                   else if (v2==1) then
                      if (v1>1) then
                         rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                              & 0.5d0*(field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,(v1-1)*ng+1)+ &
                              & field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,(v1-1)*ng))* &
                              & gll%weigh(x2)/mesh%jac(x1,nv+1)
                      else if (v1==1) then
                         if (bc_v==sll_periodic) then
                            ! activate only for periodicity in v
                            rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                                 & 0.5d0*(field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,1)+ &
                                 & field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,nv*ng))* &
                                 & gll%weigh(x2)/mesh%jac(x1,nv+1)
                         else
                            rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
                                 & 0.5d0*field_e((x1-1)*ng+x2)*dist((x1-1)*ng+x2,(v1-1)*ng+1)* &
                                 & gll%weigh(x2)/mesh%jac(x1,nv+1)
                         end if
                      end if
                   end if
                end if

                if (present(t)) then
                   !!!>>only to test Blanca's case
                   rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
                        & exp(-0.25d0*(4.0d0*v-1.0d0)**2)*( ((4.0d0*sqrt(sll_pi)+2.0d0)*v- &
                        & (2.0d0*sll_pi+sqrt(sll_pi)))*sin(2.0d0*x-2.0d0*sll_pi*t)+ &
                        & sqrt(sll_pi)*(0.25d0-v)*sin(4.0d0*x-4.0d0*sll_pi*t) )* &
                        & gll%weigh(x2)*gll%weigh(v2)/mesh%jac(x1,v1)
                   !!!<<
                end if

                rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                     & mesh%jac(x1,v1)/(gll%weigh(x2)*gll%weigh(v2))

             end do
          end do
       end do
    end do
  end subroutine rhs4dg_1d

end module sll_timestep4dg
