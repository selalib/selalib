!**************************************************************
!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************


!solves \sum_{i,j=1}^2 A_{i,j}\partial_{i,j} phi
!       +\sum_{i=1}^2B_i\partial_i phi
!       +C \phi = rho
!in polar coordinates
!this leads when A_{1,2}=A_{2,1}=0 and B_2 = 0
! A_11\partial_{1,1}\hat{phi}+B_1\partial_{1}\hat{phi}+(C+A_{2,2}k^2)\hat{phi} = \hat{rho}


module sll_module_poisson_2d_mudpack_curvilinear_solver_old
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!use sll_boundary_condition_descriptors
use sll_module_poisson_2d_base
use sll_mudpack_base
use sll_cubic_spline_interpolator_1d
use sll_cubic_spline_interpolator_2d
use sll_coordinate_transformation_2d_base_module
use sll_module_coordinate_transformations_2d
!use sll_poisson_2d_polar
implicit none


  

  type,extends(sll_poisson_2d_base) :: poisson_2d_mudpack_curvilinear_solver     
  
  !type(sll_plan_poisson_polar), pointer                   :: poiss
  sll_real64, dimension(:,:), pointer :: cxx_2d
  sll_real64, dimension(:,:), pointer :: cxy_2d
  sll_real64, dimension(:,:), pointer :: cyy_2d
  sll_real64, dimension(:,:), pointer :: cx_2d
  sll_real64, dimension(:,:), pointer :: cy_2d
  sll_real64, dimension(:,:), pointer :: ce_2d
  class(sll_interpolator_2d_base), pointer   :: cxx_2d_interp
  class(sll_interpolator_2d_base), pointer   :: cyy_2d_interp
  class(sll_interpolator_2d_base), pointer   :: cxy_2d_interp
  class(sll_interpolator_2d_base), pointer   :: cx_2d_interp
  class(sll_interpolator_2d_base), pointer   :: cy_2d_interp
  class(sll_interpolator_2d_base), pointer   :: ce_2d_interp
  class(sll_interpolator_2d_base), pointer   :: a11_interp
  class(sll_interpolator_2d_base), pointer   :: a22_interp
  class(sll_interpolator_2d_base), pointer   :: a12_interp
  class(sll_interpolator_2d_base), pointer   :: a21_interp
  sll_int32  :: mudpack_curvilinear_case
  sll_real64, dimension(:), pointer :: work !< array for tmp data
  sll_int32  :: mgopt(4) !< Option to control multigrid
  sll_int32  :: iprm(16) !< Indices to control grid sizes
  sll_real64 :: fprm(6)  !< Real to set boundary conditions
  sll_int32  :: iguess   !< Initial solution or loop over time


  
  contains
    procedure, pass(poisson) :: initialize => &
      initialize_poisson_2d_mudpack_curvilinear_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_2d_mudpack_curvilinear
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_2d_mudpack_curvilinear
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar
      
  end type poisson_2d_mudpack_curvilinear_solver

  class(poisson_2d_mudpack_curvilinear_solver), pointer   :: mudpack_curvilinear_wrapper => null()


contains
  function new_poisson_2d_mudpack_curvilinear_solver( &
    transf, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2, &    
    bc_eta1_left, &
    bc_eta1_right, &
    bc_eta2_left, &
    bc_eta2_right, &
    b11, &
    b12, &
    b21, &
    b22, &
    c) &
    result(poisson)
      
    type(poisson_2d_mudpack_curvilinear_solver),pointer :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: nc_eta2
    sll_int32, intent(in) :: bc_eta1_left
    sll_int32, intent(in) :: bc_eta1_right
    sll_int32, intent(in) :: bc_eta2_left
    sll_int32, intent(in) :: bc_eta2_right
    !sll_int32, intent(in) :: mudpack_curvilinear_case
    sll_real64, dimension(:,:), intent(in) :: b11
    sll_real64, dimension(:,:), intent(in) :: b12
    sll_real64, dimension(:,:), intent(in) :: b21
    sll_real64, dimension(:,:), intent(in) :: b22
    sll_real64, dimension(:,:), intent(in) :: c
    class(sll_coordinate_transformation_2d_base), pointer, intent(in) :: transf

    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
      
                    
    call initialize_poisson_2d_mudpack_curvilinear_solver( &
      poisson, &
      transf, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2, &    
      bc_eta1_left, &
      bc_eta1_right, &
      bc_eta2_left, &
      bc_eta2_right, &
      b11, &
      b12, &
      b21, &
      b22, &
      c)
    
  end function new_poisson_2d_mudpack_curvilinear_solver
  
  
  subroutine initialize_poisson_2d_mudpack_curvilinear_solver( &
    poisson, &
    transf, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2, &    
    bc_eta1_left, &
    bc_eta1_right, &
    bc_eta2_left, &
    bc_eta2_right, &
    b11, &
    b12, &
    b21, &
    b22, &
    c)
    class(poisson_2d_mudpack_curvilinear_solver), target :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: nc_eta2
    sll_int32, intent(in) :: bc_eta1_left
    sll_int32, intent(in) :: bc_eta1_right
    sll_int32, intent(in) :: bc_eta2_left
    sll_int32, intent(in) :: bc_eta2_right
    sll_real64, dimension(:,:), intent(in) :: b11
    sll_real64, dimension(:,:), intent(in) :: b12
    sll_real64, dimension(:,:), intent(in) :: b21
    sll_real64, dimension(:,:), intent(in) :: b22
    sll_real64, dimension(:,:), intent(in) :: c
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64,dimension(:,:),allocatable :: a12_array
    sll_real64,dimension(:,:),allocatable :: a21_array
    sll_int32 :: ierr
    !!!! begin variables for mudpack_curvilinear 
    sll_int32,  parameter   :: iixp = 2 , jjyq = 2
    sll_int32               :: icall, iiex, jjey, llwork
    sll_real64, pointer :: phi(:) !< electric potential
    sll_real64, pointer :: rhs(:) !< charge density
    !put integer and floating point argument names in contiguous
    !storeage for labelling in vectors iprm,fprm
    sll_int32  :: iprm(16)
    sll_real64 :: fprm(6)
    sll_int32  :: i,error
    sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
    sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
    common/itmud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
    sll_real64 :: xa,xb,yc,yd,tolmax,relmax
    common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax
    equivalence(intl,iprm)
    equivalence(xa,fprm)
    !declare coefficient and boundary condition input subroutines external
    !external mudpack_curvilinear_cofx,mudpack_curvilinear_cofy,mudpack_curvilinear_bndsp
    !external mudpack_curvilinear_cof
    external mudpack_curvilinear_cofcr
    external mudpack_curvilinear_bndcr
    !!!! end variables for mudpack_curvilinear 
    sll_real64 :: delta1,delta2


    nx = nc_eta1+1
    ny = nc_eta2+1

    delta1   = (eta1_max - eta1_min)/real(nc_eta1,f64)
    delta2   = (eta2_max - eta2_min)/real(nc_eta2,f64)
    ! set minimum required work space
    llwork=(7*(nx+2)*(ny+2)+44*nx*ny)/3

    allocate(poisson%work(llwork))
    icall = 0
    iiex = ceiling(log((nx-1.)/iixp)/log(2.))+1
    jjey = ceiling(log((ny-1.)/jjyq)/log(2.))+1

    !set input integer arguments
    intl = 0

    !set boundary condition flags
    nxa  = bc_eta1_left
    nxb  = bc_eta1_right
    nyc  = bc_eta2_left
    nyd  = bc_eta2_right

    !set grid sizes from parameter statements
    ixp  = iixp
    jyq  = jjyq
    iex  = iiex
    jey  = jjey

    nx = ixp*(2**(iex-1))+1
    ny = jyq*(2**(jey-1))+1

    if (nx /= nc_eta1+1 .or. ny /= nc_eta2+1) then
      print*, "nx,nc_eta1+1=", nx, nc_eta1+1
      print*, "ny,nc_eta2+1=", ny, nc_eta2+1
      stop ' nx or ny different in sll_mudpack_curvilinear_cartesian '
    end if

    !set multigrid arguments (w(2,1) cycling with fully weighted
    !residual restriction and cubic prolongation)
    poisson%mgopt(1) = 2
    poisson%mgopt(2) = 2
    poisson%mgopt(3) = 1
    poisson%mgopt(4) = 3

    !set for three cycles to ensure second-order approximation is computed
    maxcy = 3

    !set no initial guess forcing full multigrid cycling
    poisson%iguess = 0
    iguess = poisson%iguess

    !set work space length approximation from parameter statement
    nwork = llwork

    !set point relaxation
    method = 0

    !set end points of solution rectangle in (x,y) space
    xa = eta1_min
    xb = eta1_max
    yc = eta2_min
    yd = eta2_max

    !set for no error control flag
    tolmax = 0.0

!    write(*,101) (iprm(i),i=1,15)
!    write(*,102) (poisson%mgopt(i),i=1,4)
!    write(*,103) xa,xb,yc,yd,tolmax
!    write(*,104) intl

!call mud2sp(iprm,fprm,this%work,cofx,cofy,bndsp,rhs,phi,this%mgopt,error)       

    poisson%cxx_2d_interp => null()
    poisson%cyy_2d_interp => null()
    poisson%cx_2d_interp => null()
    poisson%cy_2d_interp => null()
    poisson%ce_2d_interp => null()
    poisson%a12_interp => null()
    poisson%a21_interp => null()
    
    poisson%mudpack_curvilinear_case =  SLL_NON_SEPARABLE_WITH_CROSS_TERMS                  
     !******SLL_NON_SEPARABLE_WITH_CROSS_TERMS)
    select case (poisson%mudpack_curvilinear_case)  
      case (SLL_NON_SEPARABLE_WITH_CROSS_TERMS)
        
       SLL_ALLOCATE(poisson%cxx_2d(nc_eta1+1,nc_eta2+1),ierr)       
       SLL_ALLOCATE(poisson%cxy_2d(nc_eta1+1,nc_eta2+1),ierr)
       SLL_ALLOCATE(poisson%cyy_2d(nc_eta1+1,nc_eta2+1),ierr)
       SLL_ALLOCATE(poisson%cx_2d(nc_eta1+1,nc_eta2+1),ierr) 
       SLL_ALLOCATE(poisson%cy_2d(nc_eta1+1,nc_eta2+1),ierr)
       SLL_ALLOCATE(poisson%ce_2d(nc_eta1+1,nc_eta2+1),ierr)      
        
    
       poisson%a12_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)        
       poisson%a21_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC) 
        SLL_ALLOCATE(a12_array(nx,ny),ierr)
        SLL_ALLOCATE(a21_array(nx,ny),ierr)  
        call a12_a21_array(b11,b12,b21,b22,transf,eta1_min, &
          eta2_min,delta1,delta2,nx,ny,a12_array,a21_array)
        call poisson%a12_interp%compute_interpolants( a12_array ) 
        call poisson%a21_interp%compute_interpolants( a21_array ) 
          
        poisson%cxx_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)              
        poisson%cyy_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)                            
        call coefxxyy_array(b11,b12,b21,b22,transf,eta1_min,eta2_min, & 
           delta1,delta2,nx,ny,poisson%cxx_2d ,poisson%cyy_2d)
        
        poisson%cxx_2d = -poisson%cxx_2d   
                      
        call poisson%cxx_2d_interp%compute_interpolants( poisson%cxx_2d )

        poisson%cyy_2d = -poisson%cyy_2d   
            
        call poisson%cyy_2d_interp%compute_interpolants( poisson%cyy_2d )       
           
         poisson%cxy_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
        call coefxy_array(b11,b12,b21,b22,transf,eta1_min,eta2_min, &
          delta1,delta2,nx,ny,poisson%cxy_2d)  

        poisson%cxy_2d = -poisson%cxy_2d   

        call poisson%cxy_2d_interp%compute_interpolants( poisson%cxy_2d )  
        
        poisson%cx_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
        call coefx_array(eta1_min,eta2_min,delta1,delta2,nx,ny, &
          poisson%cxx_2d_interp,poisson%a21_interp,poisson%cx_2d)

        poisson%cx_2d = -poisson%cx_2d   
            
        call poisson%cx_2d_interp%compute_interpolants( poisson%cx_2d )          



        poisson%cy_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
        call coefy_array(eta1_min,eta2_min,delta1,delta2,nx,ny, &
          poisson%cyy_2d_interp,poisson%a12_interp,poisson%cy_2d)  

        poisson%cy_2d = -poisson%cy_2d   

        call poisson%cy_2d_interp%compute_interpolants( poisson%cy_2d )          

        poisson%ce_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)       
        poisson%ce_2d = -c 
        call poisson%ce_2d_interp%compute_interpolants( poisson%ce_2d )          

                 
        if(associated(mudpack_curvilinear_wrapper))then
          print *,'#Problem mudpack_curvilinear_wrapper is not null()'
          stop
        endif
        mudpack_curvilinear_wrapper => poisson
        call mud2cr(iprm,fprm,poisson%work, &
          mudpack_curvilinear_cofcr, &
          mudpack_curvilinear_bndcr, &
          rhs, &
          phi, &
          poisson%mgopt, &
          error)
        mudpack_curvilinear_wrapper => null() 
      case default
        print *,'#bad mudpack_curvilinear_case',poisson%mudpack_curvilinear_case
        print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear_solver'
        stop 
     end select
        
  end subroutine initialize_poisson_2d_mudpack_curvilinear_solver
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_2d_mudpack_curvilinear( poisson, phi, rho )
    class(poisson_2d_mudpack_curvilinear_solver), target :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    !sll_real64        :: phi(:,:)  !< Electric potential
    !sll_real64        :: rhs(:,:)  !< Charge density
    !put integer and floating point argument names in contiguous
    !storeage for labelling in vectors iprm,fprm
    sll_int32  :: iprm(16)
    sll_real64 :: fprm(6)
    sll_int32  :: error
    sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
    sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
    common/itmud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
    sll_real64 :: xa,xb,yc,yd,tolmax,relmax
    common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

    equivalence(intl,iprm)
    equivalence(xa,fprm)

    !declare coefficient and boundary condition input subroutines external
    !external mudpack_curvilinear_cofx,mudpack_curvilinear_cofy,mudpack_curvilinear_bndsp
    !external mudpack_curvilinear_cof
    external mudpack_curvilinear_cofcr
    external mudpack_curvilinear_bndcr
    !set initial guess because solve should be called every time step in a
    !time dependent problem and the elliptic operator does not depend on time.
    iguess = poisson%iguess

    !attempt solution
    intl = 1
    !write(*,106) intl,method,iguess


    select case (poisson%mudpack_curvilinear_case)
          
      case (SLL_NON_SEPARABLE_WITH_CROSS_TERMS)
        if(associated(mudpack_curvilinear_wrapper))then
          print *,'#Problem mudpack_curvilinear_wrapper is not null()'
          stop
        endif
        mudpack_curvilinear_wrapper => poisson
        call mud2cr(iprm, &
          fprm, &
          poisson%work, &
          mudpack_curvilinear_cofcr, &
          mudpack_curvilinear_bndcr, &
          rho, &
          phi, &
          poisson%mgopt, &
          error)
        !write(*,107) error
        if (error > 0) call exit(0)
        ! attempt to improve approximation to fourth order
        ! seems not to work for the moment
        !call mud24cr(poisson%work,phi,error)
        !write (*,108) error
        if (error > 0) call exit(0)
        
         mudpack_curvilinear_wrapper => null()
         
      case default
        print *,'#bad mudpack_curvilinear_case',poisson%mudpack_curvilinear_case
        print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear_solver'
        stop 
    end select



    
    !call solve( poisson%poiss, rho, phi)
    
  end subroutine compute_phi_from_rho_2d_mudpack_curvilinear

    ! solves E = -\nabla Phi in 2d
!    subroutine compute_E_from_phi_2d_fft( poisson, phi, E1, E2 )
!      class(poisson_2d_fft_solver) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine compute_E_from_phi_2d_fft

    ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
    subroutine compute_E_from_rho_2d_mudpack_curvilinear( poisson, rho, E1, E2 )
      class(poisson_2d_mudpack_curvilinear_solver) :: poisson
      sll_real64,dimension(:,:),intent(in) :: rho
      sll_real64,dimension(:,:),intent(out) :: E1
      sll_real64,dimension(:,:),intent(out) :: E2
      
      print *,'#compute_E_from_rho_2d_mudpack_curvilinear'      
      print *,'#not implemented for the moment'
      stop
      
      !call solve( poisson%poiss, E1, E2, rho)
      
    end subroutine compute_E_from_rho_2d_mudpack_curvilinear
  
subroutine coefxxyy_array(b11,b12,b21,b22,transf,eta1_min,eta2_min, &
                         delta1,delta2,nx,ny,cxx_array,cyy_array)
  implicit none                     
    sll_real64                :: eta1,eta1_min,eta2_min
    sll_real64                :: eta2,delta1,delta2
    sll_int32                 :: i,j,nx,ny
    sll_real64, dimension(:,:):: cxx_array,cyy_array
    sll_real64, dimension(1:2,1:2) :: jac_m
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64, dimension(:,:) :: b11
    sll_real64, dimension(:,:) :: b12
    sll_real64, dimension(:,:) :: b21 
    sll_real64, dimension(:,:) :: b22
    
do j=1,ny
 eta2 = eta2_min + real(j-1,f64)*delta2
 do i=1,nx
   eta1 = eta1_min + real(i-1,f64)*delta1
   jac_m  =  transf%jacobian_matrix(eta1,eta2) 
   cxx_array(i,j)= (b11(i,j)*(jac_m(1,2)*jac_m(1,2)+jac_m(2,2)*jac_m(2,2))- &
                   & b12(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2))) &
                   /transf%jacobian(eta1,eta2)
   cyy_array(i,j)= (b22(i,j)*(jac_m(2,1)*jac_m(2,1)+jac_m(1,1)*jac_m(1,1))- &
                   & b21(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2))) &
                   /transf%jacobian(eta1,eta2)                
 enddo
enddo 
end subroutine coefxxyy_array

subroutine coefxy_array(b11,b12,b21,b22,transf,eta1_min,eta2_min, &
                         delta1,delta2,nx,ny,cxy_array)
  implicit none                     
    sll_real64                :: eta1,eta1_min,eta2_min
    sll_real64                :: eta2,delta1,delta2
    sll_real64                :: a12,a21
    sll_int32                 :: i,j,nx,ny
    sll_real64, dimension(:,:):: cxy_array
    sll_real64, dimension(1:2,1:2) :: jac_m
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64, dimension(:,:) :: b11
    sll_real64, dimension(:,:) :: b12
    sll_real64, dimension(:,:) :: b21 
    sll_real64, dimension(:,:) :: b22
    
do j=1,ny
 eta2 = eta2_min + real(j-1,f64)*delta2
 do i=1,nx
   eta1 = eta1_min + real(i-1,f64)*delta1
   jac_m  =  transf%jacobian_matrix(eta1,eta2) 
   a12= b12(i,j)*(jac_m(2,1)*jac_m(2,1)+jac_m(1,1)*jac_m(1,1))- &
                   & b11(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2))
                   
   a21= b21(i,j)*(jac_m(1,2)*jac_m(1,2)+jac_m(2,2)*jac_m(2,2))- &
                   & b22(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2))  
   cxy_array(i,j)= (a12+a21)/transf%jacobian(eta1,eta2)                              
 enddo
enddo 
end subroutine coefxy_array

subroutine a12_a21_array(b11,b12,b21,b22,transf,eta1_min,eta2_min,delta1,delta2,nx,ny,a12_array,a21_array)
  implicit none                     
    sll_real64                :: eta1,eta1_min,eta2_min
    sll_real64                :: eta2,delta1,delta2
    sll_real64                :: a12,a21
    sll_int32                 :: i,j,nx,ny
    sll_real64, dimension(:,:):: a12_array
    sll_real64, dimension(:,:):: a21_array
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64, dimension(:,:) :: b11
    sll_real64, dimension(:,:) :: b12
    sll_real64, dimension(:,:) :: b21 
    sll_real64, dimension(:,:) :: b22
    sll_real64, dimension(1:2,1:2) :: jac_m
    
do j=1,ny
 eta2 = eta2_min + real(j-1,f64)*delta2
 do i=1,nx
   eta1 = eta1_min + real(i-1,f64)*delta1
   jac_m  =  transf%jacobian_matrix(eta1,eta2) 
   a12= b12(i,j)*(jac_m(2,1)*jac_m(2,1)+jac_m(1,1)*jac_m(1,1))- &
                   & b11(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2))
   a12_array(i,j) = a12/transf%jacobian(eta1,eta2)                
   a21= b21(i,j)*(jac_m(1,2)*jac_m(1,2)+jac_m(2,2)*jac_m(2,2))- &
                   & b22(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2)) 
   a21_array(i,j) = a21/transf%jacobian(eta1,eta2)                                             
 enddo
enddo 
end subroutine a12_a21_array

subroutine coefx_array(eta1_min,eta2_min, &
                         delta1,delta2,nx,ny,cxx_2d_interp,a21_interp,cx_array)
  implicit none                     
    sll_real64                :: eta1,eta1_min,eta2_min
    sll_real64                :: eta2,delta1,delta2
    sll_int32                 :: i,j,nx,ny
    sll_real64, dimension(:,:):: cx_array
    class(sll_interpolator_2d_base), pointer   :: cxx_2d_interp
    class(sll_interpolator_2d_base), pointer   :: a21_interp
    
    
do j=1,ny
 eta2 = eta2_min + real(j-1,f64)*delta2
 do i=1,nx
   eta1 = eta1_min + real(i-1,f64)*delta1   
   cx_array(i,j)= cxx_2d_interp%interpolate_derivative_eta1(eta1,eta2)+ &
                  a21_interp%interpolate_derivative_eta2(eta1,eta2)                         
 enddo
enddo 
end subroutine coefx_array
subroutine coefy_array(eta1_min,eta2_min, &
                         delta1,delta2,nx,ny,cyy_2d_interp,a12_interp,cy_array)
  implicit none                     
    sll_real64                :: eta1,eta1_min,eta2_min
    sll_real64                :: eta2,delta1,delta2
    sll_int32                 :: i,j,nx,ny
    sll_real64, dimension(:,:):: cy_array
    class(sll_interpolator_2d_base), pointer   :: cyy_2d_interp
    class(sll_interpolator_2d_base), pointer   :: a12_interp
    
do j=1,ny
 eta2 = eta2_min + real(j-1,f64)*delta2
 do i=1,nx
   eta1 = eta1_min + real(i-1,f64)*delta1    
   cy_array(i,j)= cyy_2d_interp%interpolate_derivative_eta2(eta1,eta2)+ &
                  a12_interp%interpolate_derivative_eta1(eta1,eta2)                         
 enddo
enddo 
end subroutine coefy_array  
  
  
end module sll_module_poisson_2d_mudpack_curvilinear_solver_old


subroutine mudpack_curvilinear_cof(x,y,cxx,cyy,cx,cy,ce)
use sll_module_poisson_2d_mudpack_curvilinear_solver_old
implicit none
real(8)  :: x,cxx,cx
real(8)  :: y,cyy,cy,ce
cxx = mudpack_curvilinear_wrapper%cxx_2d_interp%interpolate_value(x,y)
cyy = mudpack_curvilinear_wrapper%cyy_2d_interp%interpolate_value(x,y)
cx  = mudpack_curvilinear_wrapper%cx_2d_interp%interpolate_value(x,y)
cy  = mudpack_curvilinear_wrapper%cy_2d_interp%interpolate_value(x,y)
ce  = mudpack_curvilinear_wrapper%ce_2d_interp%interpolate_value(x,y)
return
end

subroutine mudpack_curvilinear_cofcr(x,y,cxx,cxy,cyy,cx,cy,ce)
use sll_module_poisson_2d_mudpack_curvilinear_solver_old
implicit none
real(8)  :: x,cxx,cx,cxy
real(8)  :: y,cyy,cy,ce
cxx = mudpack_curvilinear_wrapper%cxx_2d_interp%interpolate_value(x,y)
cxy = mudpack_curvilinear_wrapper%cxy_2d_interp%interpolate_value(x,y)
cyy = mudpack_curvilinear_wrapper%cyy_2d_interp%interpolate_value(x,y)
cx  = mudpack_curvilinear_wrapper%cx_2d_interp%interpolate_value(x,y)
cy  = mudpack_curvilinear_wrapper%cy_2d_interp%interpolate_value(x,y)
ce  = mudpack_curvilinear_wrapper%ce_2d_interp%interpolate_value(x,y)

return
end
!> input mixed derivative b.c. to mud2sp
subroutine mudpack_curvilinear_bndcr(kbdy,xory,alfa,gbdy)
use sll_module_poisson_2d_mudpack_curvilinear_solver_old
implicit none
integer  :: kbdy
real(8)  :: xory,alfa,gbdy,x,y,pe,px,py
real(8)  :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

!subroutine not used in periodic case
if (kbdy == 1) then  ! x=xa boundary
   y = xory
   x = xa
   alfa = -1.0
   gbdy = px + alfa*pe
   return
end if

if (kbdy == 4) then  ! y=yd boundary
   y = yd
   x = xory
   alfa = 1.0
   gbdy = py + alfa*pe
   return
end if
end















