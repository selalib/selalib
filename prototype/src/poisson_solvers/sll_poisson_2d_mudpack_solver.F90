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
!> @author
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
!> Adnane Hamiaz (hamiaz@math.unistra.fr)
!**************************************************************

!solves \sum_{i,j=1}^2 A_{i,j}\partial_{i,j} phi
!       +\sum_{i=1}^2B_i\partial_i phi
!       +C \phi = rho
!in polar coordinates
!this leads when A_{1,2}=A_{2,1}=0 and B_2 = 0
! A_11\partial_{1,1}\hat{phi}+B_1\partial_{1}\hat{phi}+(C+A_{2,2}k^2)\hat{phi} = \hat{rho}


module sll_module_poisson_2d_mudpack_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!use sll_boundary_condition_descriptors
use sll_module_poisson_2d_base
use sll_mudpack_base
use sll_cubic_spline_interpolator_1d
use sll_cubic_spline_interpolator_2d

!use sll_poisson_2d_polar
implicit none

  !integer, parameter :: SLL_SEPARABLE  = 1    !< type of equation
  !integer, parameter :: SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS = 2    !< type of equation
  !integer, parameter :: SLL_NON_SEPARABLE_WITH_CROSS_TERMS = 3    !< type of equation

  

  type,extends(sll_poisson_2d_base) :: poisson_2d_mudpack_solver     
  
  !type(sll_plan_poisson_polar), pointer                   :: poiss
  sll_real64, dimension(:,:), pointer :: cxx_2d
  sll_real64, dimension(:,:), pointer :: cxy_2d
  sll_real64, dimension(:,:), pointer :: cyy_2d
  sll_real64, dimension(:,:), pointer :: cx_2d
  sll_real64, dimension(:,:), pointer :: cy_2d
  sll_real64, dimension(:,:), pointer :: ce_2d
  sll_real64, dimension(:), pointer :: cxx_1d
  sll_real64, dimension(:), pointer :: cyy_1d
  sll_real64, dimension(:), pointer :: cx_1d
  sll_real64, dimension(:), pointer :: cy_1d
  sll_real64, dimension(:), pointer :: cex_1d
  sll_real64, dimension(:), pointer :: cey_1d
  sll_real64 :: cxx
  sll_real64 :: cyy
  sll_real64 :: cx
  sll_real64 :: cy
  sll_real64 :: ce
  sll_int32  :: mudpack_case
  class(sll_interpolator_2d_base), pointer   :: cxx_2d_interp
  class(sll_interpolator_2d_base), pointer   :: cxy_2d_interp
  class(sll_interpolator_2d_base), pointer   :: cyy_2d_interp
  class(sll_interpolator_2d_base), pointer   :: cx_2d_interp
  class(sll_interpolator_2d_base), pointer   :: cy_2d_interp
  class(sll_interpolator_2d_base), pointer   :: ce_2d_interp
  class(sll_interpolator_1d_base), pointer   :: cxx_1d_interp
  class(sll_interpolator_1d_base), pointer   :: cyy_1d_interp
  class(sll_interpolator_1d_base), pointer   :: cx_1d_interp
  class(sll_interpolator_1d_base), pointer   :: cy_1d_interp
  class(sll_interpolator_1d_base), pointer   :: cex_1d_interp
  class(sll_interpolator_1d_base), pointer   :: cey_1d_interp


  sll_real64, dimension(:), pointer :: work !< array for tmp data
  sll_int32  :: mgopt(4) !< Option to control multigrid
  sll_int32  :: iprm(16) !< Indices to control grid sizes
  sll_real64 :: fprm(6)  !< Real to set boundary conditions
  sll_int32  :: iguess   !< Initial solution or loop over time


  
  
  contains
    procedure, pass(poisson) :: initialize => &
      initialize_poisson_2d_mudpack_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_2d_mudpack
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_2d_mudpack
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar
      
  end type poisson_2d_mudpack_solver

  class(poisson_2d_mudpack_solver), pointer   :: mudpack_wrapper => null()


contains
  function new_poisson_2d_mudpack_solver( &
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
    mudpack_case, &
    cxx_2d, &
    cxy_2d, &
    cyy_2d, &
    cx_2d, &
    cy_2d, &
    ce_2d, &
    cxx_1d, &
    cyy_1d, &
    cx_1d, &
    cy_1d, &
    cex_1d, &
    cey_1d, &
    cxx, &
    cxy, &
    cyy, &
    cx, &
    cy, &
    ce) &
    result(poisson)
      
    type(poisson_2d_mudpack_solver),pointer :: poisson
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
    sll_int32, intent(in) :: mudpack_case
    sll_real64, dimension(:,:), intent(in), optional :: cxx_2d
    sll_real64, dimension(:,:), intent(in), optional :: cxy_2d
    sll_real64, dimension(:,:), intent(in), optional :: cyy_2d
    sll_real64, dimension(:,:), intent(in), optional :: cx_2d
    sll_real64, dimension(:,:), intent(in), optional :: cy_2d
    sll_real64, dimension(:,:), intent(in), optional :: ce_2d
    sll_real64, dimension(:), intent(in), optional :: cxx_1d
    sll_real64, dimension(:), intent(in), optional :: cyy_1d
    sll_real64, dimension(:), intent(in), optional :: cx_1d
    sll_real64, dimension(:), intent(in), optional :: cy_1d
    sll_real64, dimension(:), intent(in), optional :: cex_1d
    sll_real64, dimension(:), intent(in), optional :: cey_1d
    sll_real64, intent(in), optional  :: cxx
    sll_real64, intent(in), optional  :: cxy
    sll_real64, intent(in), optional  :: cyy
    sll_real64, intent(in), optional  :: cx
    sll_real64, intent(in), optional  :: cy
    sll_real64, intent(in), optional  :: ce

    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_2d_mudpack_solver( &
      poisson, &
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
      mudpack_case, &
      cxx_2d, &
      cxy_2d, &
      cyy_2d, &
      cx_2d, &
      cy_2d, &
      ce_2d, &
      cxx_1d, &
      cyy_1d, &
      cx_1d, &
      cy_1d, &
      cex_1d, &
      cey_1d, &
      cxx, &
      cxy, &
      cyy, &
      cx, &
      cy, &
      ce)
    
  end function new_poisson_2d_mudpack_solver
  
  
  subroutine initialize_poisson_2d_mudpack_solver( &
    poisson, &
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
    mudpack_case, &
    cxx_2d, &
    cxy_2d, &
    cyy_2d, &
    cx_2d, &
    cy_2d, &
    ce_2d, &
    cxx_1d, &
    cyy_1d, &
    cx_1d, &
    cy_1d, &
    cex_1d, &
    cey_1d, &
    cxx, &
    cxy, &
    cyy, &
    cx, &
    cy, &
    ce)
    class(poisson_2d_mudpack_solver), target :: poisson
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
    sll_int32, intent(in) :: mudpack_case
    sll_real64, dimension(:,:), intent(in), optional :: cxx_2d
    sll_real64, dimension(:,:), intent(in), optional :: cxy_2d
    sll_real64, dimension(:,:), intent(in), optional :: cyy_2d
    sll_real64, dimension(:,:), intent(in), optional :: cx_2d
    sll_real64, dimension(:,:), intent(in), optional :: cy_2d
    sll_real64, dimension(:,:), intent(in), optional :: ce_2d
    sll_real64, dimension(:), intent(in), optional :: cxx_1d
    sll_real64, dimension(:), intent(in), optional :: cyy_1d
    sll_real64, dimension(:), intent(in), optional :: cx_1d
    sll_real64, dimension(:), intent(in), optional :: cy_1d
    sll_real64, dimension(:), intent(in), optional :: cex_1d
    sll_real64, dimension(:), intent(in), optional :: cey_1d
    sll_real64, intent(in), optional  :: cxx
    sll_real64, intent(in), optional  :: cxy
    sll_real64, intent(in), optional  :: cyy
    sll_real64, intent(in), optional  :: cx
    sll_real64, intent(in), optional  :: cy
    sll_real64, intent(in), optional  :: ce
    sll_int32 :: ierr
    !!!! begin variables for mudpack 
    sll_int32,  parameter   :: iixp = 2 , jjyq = 2
    sll_int32               :: icall, iiex, jjey, llwork
    sll_real64, pointer :: phi(:) !< electric potential
    sll_real64, pointer :: rhs(:) !< charge density
    !put integer and floating point argument names in contiguous
    !storeage for labelling in vectors iprm,fprm
    sll_int32  :: iprm(16)
    sll_real64 :: fprm(6)
    !sll_int32  :: i
    sll_int32 :: error
    sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
    sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
    common/itmud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
    sll_real64 :: xa,xb,yc,yd,tolmax,relmax
    common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax
    equivalence(intl,iprm)
    equivalence(xa,fprm)
    !declare coefficient and boundary condition input subroutines external
    external mudpack_cofx,mudpack_cofy,mudpack_bndsp
    external mudpack_cof,mudpack_cofcr
    !!!! end variables for mudpack 



    nx = nc_eta1+1
    ny = nc_eta2+1
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
      stop ' nx or ny different in sll_mudpack_cartesian '
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







        
    poisson%mudpack_case = mudpack_case 

    poisson%cxx_2d_interp => null()
    poisson%cxy_2d_interp => null()
    poisson%cyy_2d_interp => null()
    poisson%cx_2d_interp => null()
    poisson%cy_2d_interp => null()
    poisson%ce_2d_interp => null()

    poisson%cxx_1d_interp => null()
    poisson%cyy_1d_interp => null()
    poisson%cx_1d_interp => null()
    poisson%cy_1d_interp => null()
    poisson%cex_1d_interp => null()
    poisson%cey_1d_interp => null()


    select case (mudpack_case)
      case (SLL_SEPARABLE)
        if(present(cxx_2d).or.present(cxy_2d).or.present(cyy_2d)&
          .or.present(cx_2d).or.present(cy_2d).or.present(ce_2d)) then
          print *,'#2d arrays should not be here'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop
        endif
        
        if((.not.(present(cxx_1d))).and.(.not.(present(cxx)))) then
          print *,'#1d/0d array should be here for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cxx_1d).and.present(cxx))then
          print *,'#please choose between 1d or 0d array for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cxx_1d))then
          if(size(cxx_1d)<nc_eta1+1)then
            print *,'#Bad size for cxx_1d',size(cxx_1d),nc_eta1+1
            stop
          endif
          SLL_ALLOCATE(poisson%cxx_1d(nc_eta1+1),ierr)
          poisson%cxx_1d(1:nc_eta1+1)=cxx_1d(1:nc_eta1+1)
        endif
        if(present(cxx))then
          SLL_ALLOCATE(poisson%cxx_1d(nc_eta1+1),ierr)
          poisson%cxx_1d(1:nc_eta1+1)=cxx
        endif

        if((.not.(present(cyy_1d))).and.(.not.(present(cyy)))) then
          print *,'#1d/0d array should be here for cyy !'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cyy_1d).and.present(cyy))then
          print *,'#please choose between 1d or 0d array for cyy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cyy_1d))then
          if(size(cyy_1d)<nc_eta2+1)then
            print *,'#Bad size for cyy_1d',size(cyy_1d),nc_eta2+1
            stop
          endif
          SLL_ALLOCATE(poisson%cyy_1d(nc_eta2+1),ierr)
          poisson%cyy_1d(1:nc_eta2+1)=cyy_1d(1:nc_eta2+1)
        endif
        if(present(cyy))then
          SLL_ALLOCATE(poisson%cyy_1d(nc_eta2+1),ierr)
          poisson%cyy_1d(1:nc_eta2+1)=cyy
        endif

        if((.not.(present(cx_1d))).and.(.not.(present(cx)))) then
          SLL_ALLOCATE(poisson%cx_1d(nc_eta1+1),ierr)
          poisson%cx_1d(1:nc_eta1+1)=0._f64
        endif
        if(present(cx_1d).and.present(cx))then
          print *,'#please choose between 1d or 0d array for cx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cx_1d))then
          if(size(cx_1d)<nc_eta1+1)then
            print *,'#Bad size for cx_1d',size(cx_1d),nc_eta1+1
            stop
          endif
          SLL_ALLOCATE(poisson%cx_1d(nc_eta1+1),ierr)
          poisson%cx_1d(1:nc_eta1+1)=cx_1d(1:nc_eta1+1)
        endif
        if(present(cx))then
          SLL_ALLOCATE(poisson%cx_1d(nc_eta1+1),ierr)
          poisson%cx_1d(1:nc_eta1+1)=cx
        endif


        if((.not.(present(cy_1d))).and.(.not.(present(cy)))) then
          SLL_ALLOCATE(poisson%cy_1d(nc_eta2+1),ierr)
          poisson%cy_1d(1:nc_eta2+1)=0._f64
        endif
        if(present(cy_1d).and.present(cy))then
          print *,'#please choose between 1d or 0d array for cy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cy_1d))then
          if(size(cy_1d)<nc_eta2+1)then
            print *,'#Bad size for cy_1d',size(cy_1d),nc_eta2+1
            stop
          endif
          SLL_ALLOCATE(poisson%cy_1d(nc_eta2+1),ierr)
          poisson%cy_1d(1:nc_eta2+1)=cy_1d(1:nc_eta2+1)
        endif
        if(present(cy))then
          SLL_ALLOCATE(poisson%cy_1d(nc_eta2+1),ierr)
          poisson%cy_1d(1:nc_eta2+1)=cy
        endif

        if((.not.(present(cex_1d))).and.(.not.(present(ce)))) then
          SLL_ALLOCATE(poisson%cex_1d(nc_eta1+1),ierr)
          poisson%cex_1d = 0._f64          
        endif
        if(present(cex_1d).and.present(ce))then
          print *,'#please choose between 1d or 0d array for cex'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cex_1d))then
          if(size(cex_1d)<nc_eta1+1)then
            print *,'#Bad size for cex_1d',size(cex_1d),nc_eta1+1
            stop
          endif
          SLL_ALLOCATE(poisson%cex_1d(nc_eta1+1),ierr)
          poisson%cex_1d(1:nc_eta1+1)=cex_1d(1:nc_eta1+1)
        endif
        if(present(ce))then
          SLL_ALLOCATE(poisson%cex_1d(nc_eta1+1),ierr)
          poisson%cex_1d(1:nc_eta1+1)=0.5_f64*ce
        endif

        if((.not.(present(cey_1d))).and.(.not.(present(ce)))) then
          SLL_ALLOCATE(poisson%cey_1d(nc_eta2+1),ierr)
          poisson%cey_1d(1:nc_eta2+1)=0._f64
        endif
        if(present(cey_1d).and.present(ce))then
          print *,'#please choose between 1d or 0d array for cey'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cey_1d))then
          if(size(cey_1d)<nc_eta2+1)then
            print *,'#Bad size for cey_1d',size(cey_1d),nc_eta2+1
            stop
          endif
          SLL_ALLOCATE(poisson%cey_1d(nc_eta2+1),ierr)
          poisson%cey_1d(1:nc_eta2+1)=cey_1d(1:nc_eta2+1)
        endif
        if(present(ce))then
          SLL_ALLOCATE(poisson%cey_1d(nc_eta2+1),ierr)
          poisson%cey_1d(1:nc_eta2+1)=0.5_f64*ce
        endif

        poisson%cxx_1d_interp => new_cubic_spline_1d_interpolator( &
          nx, &
          eta1_min, &
          eta1_max, &
          SLL_PERIODIC)          
        call poisson%cxx_1d_interp%compute_interpolants( poisson%cxx_1d )          

        poisson%cyy_1d_interp => new_cubic_spline_1d_interpolator( &
          ny, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC)          
        call poisson%cyy_1d_interp%compute_interpolants( poisson%cyy_1d )          

        poisson%cx_1d_interp => new_cubic_spline_1d_interpolator( &
          nx, &
          eta1_min, &
          eta1_max, &
          SLL_PERIODIC)          
        call poisson%cx_1d_interp%compute_interpolants( poisson%cx_1d )          

        poisson%cy_1d_interp => new_cubic_spline_1d_interpolator( &
          ny, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC)          
        call poisson%cy_1d_interp%compute_interpolants( poisson%cy_1d )          

        poisson%cex_1d_interp => new_cubic_spline_1d_interpolator( &
          nx, &
          eta1_min, &
          eta1_max, &
          SLL_PERIODIC)          
        call poisson%cex_1d_interp%compute_interpolants( poisson%cex_1d )          

        poisson%cey_1d_interp => new_cubic_spline_1d_interpolator( &
          ny, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC)          
        call poisson%cey_1d_interp%compute_interpolants( poisson%cey_1d )          



        
        if(associated(mudpack_wrapper))then
          print *,'#Problem mudpack_wrapper is not null()'
          stop
        endif
        mudpack_wrapper => poisson
        call mud2sp(iprm,fprm,poisson%work, &
          mudpack_cofx, &
          mudpack_cofy, &
          mudpack_bndsp, &
          rhs, &
          phi, &
          poisson%mgopt, &
          error)
        mudpack_wrapper => null() 

                
        
      case (SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS)
      
        if(present(cxx_1d).or.present(cyy_1d).or.present(cx_1d)& 
          .or.present(cy_1d).or.present(cex_1d).or.present(cey_1d)) then
          print *,'#1d arrays should not be here'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop
        endif
        
        if((.not.(present(cxx_2d))).and.(.not.(present(cxx)))) then
          print *,'#2d/0d array should be here for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cxx_2d).and.present(cxx))then
          print *,'#please choose between 2d or 0d array for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cxx_2d))then
          if(size(cxx_2d)< (nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for cxx_2d',size(cxx_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%cxx_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cxx_2d(1:nc_eta1+1,1:nc_eta2+1)=cxx_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(cxx))then
          SLL_ALLOCATE(poisson%cxx_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cxx_2d(1:nc_eta1+1,1:nc_eta2+1)=cxx
        endif

        if((.not.(present(cyy_2d))).and.(.not.(present(cyy)))) then
          print *,'#2d/0d array should be here for cyy !'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cyy_2d).and.present(cyy))then
          print *,'#please choose between 2d or 0d array for cyy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cyy_2d))then
          if(size(cyy_2d)<(nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for cyy_2d',size(cyy_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%cyy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cyy_2d(1:nc_eta1+1,1:nc_eta2+1)=cyy_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(cyy))then
          SLL_ALLOCATE(poisson%cyy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cyy_2d(1:nc_eta1+1,1:nc_eta2+1)=cyy
        endif

        if((.not.(present(cx_2d))).and.(.not.(present(cx)))) then
          SLL_ALLOCATE(poisson%cx_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cx_2d(1:nc_eta1+1,1:nc_eta2+1)=0._f64
        endif
        if(present(cx_2d).and.present(cx))then
          print *,'#please choose between 2d or 0d array for cx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cx_2d))then
          if(size(cx_2d)<(nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for cx_2d',size(cx_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%cx_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cx_2d(1:nc_eta1+1,1:nc_eta2+1)=cx_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(cx))then
          SLL_ALLOCATE(poisson%cx_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cx_2d(1:nc_eta1+1,1:nc_eta2+1)=cx
        endif

        if((.not.(present(cy_2d))).and.(.not.(present(cy)))) then
          SLL_ALLOCATE(poisson%cy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cy_2d(1:nc_eta1+1,1:nc_eta2+1)=0._f64
        endif
        if(present(cy_2d).and.present(cy))then
          print *,'#please choose between 2d or 0d array for cy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cy_2d))then
          if(size(cy_2d)<(nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for cy_2d',size(cy_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%cy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cy_2d(1:nc_eta1+1,1:nc_eta2+1)=cy_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(cy))then
          SLL_ALLOCATE(poisson%cy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cy_2d(1:nc_eta1+1,1:nc_eta2+1)=cy
        endif
        
        if((.not.(present(ce_2d))).and.(.not.(present(ce)))) then
          SLL_ALLOCATE(poisson%ce_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%ce_2d = 0._f64          
        endif
        if(present(ce_2d).and.present(ce))then
          print *,'#please choose between 2d or 0d array for ce'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(ce_2d))then
          if(size(ce_2d)<(nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for ce_2d',size(ce_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%ce_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%ce_2d(1:nc_eta1+1,1:nc_eta2+1)=ce_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(ce))then
          SLL_ALLOCATE(poisson%ce_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%ce_2d(1:nc_eta1+1,1:nc_eta2+1)=ce
        endif


        poisson%cxx_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
        call poisson%cxx_2d_interp%compute_interpolants( poisson%cxx_2d )          

        poisson%cyy_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
        call poisson%cyy_2d_interp%compute_interpolants( poisson%cyy_2d )          

        poisson%cx_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
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
        call poisson%ce_2d_interp%compute_interpolants( poisson%ce_2d )          

                 
        if(associated(mudpack_wrapper))then
          print *,'#Problem mudpack_wrapper is not null()'
          stop
        endif
        mudpack_wrapper => poisson
        call mud2(iprm,fprm,poisson%work, &
          mudpack_cof, &
          mudpack_bndsp, &
          rhs, &
          phi, &
          poisson%mgopt, &
          error)
        mudpack_wrapper => null() 

      case (SLL_NON_SEPARABLE_WITH_CROSS_TERMS)
      
        if(present(cxx_1d).or.present(cyy_1d).or.present(cx_1d)& 
          .or.present(cy_1d).or.present(cex_1d).or.present(cey_1d)) then
          print *,'#1d arrays should not be here'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop
        endif
        
        if((.not.(present(cxx_2d))).and.(.not.(present(cxx)))) then
          print *,'#2d/0d array should be here for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cxx_2d).and.present(cxx))then
          print *,'#please choose between 2d or 0d array for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cxx_2d))then
          if(size(cxx_2d)< (nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for cxx_2d',size(cxx_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%cxx_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cxx_2d(1:nc_eta1+1,1:nc_eta2+1)=cxx_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(cxx))then
          SLL_ALLOCATE(poisson%cxx_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cxx_2d(1:nc_eta1+1,1:nc_eta2+1)=cxx
        endif
        
        if((.not.(present(cxy_2d))).and.(.not.(present(cxy)))) then
          print *,'#2d array should be here for cxy'
          print *,'# cxy=0. use another method'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cxy_2d).and.present(cxy))then
          print *,'#please choose between 2d or 0d array for cxy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif        
        if(present(cxy_2d))then
          if(size(cxy_2d)< (nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for cxy_2d',size(cxy_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%cxy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cxy_2d(1:nc_eta1+1,1:nc_eta2+1)=cxy_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(cxy))then
          SLL_ALLOCATE(poisson%cxy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cxy_2d(1:nc_eta1+1,1:nc_eta2+1)=cxy
        endif
        if((.not.(present(cyy_2d))).and.(.not.(present(cyy)))) then
          print *,'#2d/0d array should be here for cyy !'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cyy_2d).and.present(cyy))then
          print *,'#please choose between 2d or 0d array for cyy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cyy_2d))then
          if(size(cyy_2d)<(nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for cyy_2d',size(cyy_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%cyy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cyy_2d(1:nc_eta1+1,1:nc_eta2+1)=cyy_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(cyy))then
          SLL_ALLOCATE(poisson%cyy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cyy_2d(1:nc_eta1+1,1:nc_eta2+1)=cyy
        endif

        if((.not.(present(cx_2d))).and.(.not.(present(cx)))) then
          SLL_ALLOCATE(poisson%cx_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cx_2d(1:nc_eta1+1,1:nc_eta2+1)=0._f64
        endif
        if(present(cx_2d).and.present(cx))then
          print *,'#please choose between 2d or 0d array for cx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cx_2d))then
          if(size(cx_2d)<(nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for cx_2d',size(cx_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%cx_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cx_2d(1:nc_eta1+1,1:nc_eta2+1)=cx_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(cx))then
          SLL_ALLOCATE(poisson%cx_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cx_2d(1:nc_eta1+1,1:nc_eta2+1)=cx
        endif

        if((.not.(present(cy_2d))).and.(.not.(present(cy)))) then
          SLL_ALLOCATE(poisson%cy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cy_2d(1:nc_eta1+1,1:nc_eta2+1)=0._f64
        endif
        if(present(cy_2d).and.present(cy))then
          print *,'#please choose between 2d or 0d array for cy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cy_2d))then
          if(size(cy_2d)<(nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for cy_2d',size(cy_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%cy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cy_2d(1:nc_eta1+1,1:nc_eta2+1)=cy_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(cy))then
          SLL_ALLOCATE(poisson%cy_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%cy_2d(1:nc_eta1+1,1:nc_eta2+1)=cy
        endif
        
        if((.not.(present(ce_2d))).and.(.not.(present(ce)))) then
          SLL_ALLOCATE(poisson%ce_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%ce_2d = 0._f64          
        endif
        if(present(ce_2d).and.present(ce))then
          print *,'#please choose between 2d or 0d array for ce'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(ce_2d))then
          if(size(ce_2d)<(nc_eta1+1)*(nc_eta2+1))then
            print *,'#Bad size for ce_2d',size(ce_2d),(nc_eta1+1)*(nc_eta2+1)
            stop
          endif
          SLL_ALLOCATE(poisson%ce_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%ce_2d(1:nc_eta1+1,1:nc_eta2+1)=ce_2d(1:nc_eta1+1,1:nc_eta2+1)
        endif
        if(present(ce))then
          SLL_ALLOCATE(poisson%ce_2d(nc_eta1+1,nc_eta2+1),ierr)
          poisson%ce_2d(1:nc_eta1+1,1:nc_eta2+1)=ce
        endif


        poisson%cxx_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
        call poisson%cxx_2d_interp%compute_interpolants( poisson%cxx_2d )   
               
        poisson%cxy_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
        call poisson%cxy_2d_interp%compute_interpolants( poisson%cxy_2d ) 
        
        poisson%cyy_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
        call poisson%cyy_2d_interp%compute_interpolants( poisson%cyy_2d )          

        poisson%cx_2d_interp => new_cubic_spline_2d_interpolator( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
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
        call poisson%ce_2d_interp%compute_interpolants( poisson%ce_2d )          

                 
        if(associated(mudpack_wrapper))then
          print *,'#Problem mudpack_wrapper is not null()'
          stop
        endif
        mudpack_wrapper => poisson
        call mud2cr(iprm,fprm,poisson%work, &
          mudpack_cofcr, &
          mudpack_bndsp, &
          rhs, &
          phi, &
          poisson%mgopt, &
          error)
        mudpack_wrapper => null() 

      case default
        print *,'#bad mudpack_case',mudpack_case
        print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
        stop 
    end select

        
  end subroutine initialize_poisson_2d_mudpack_solver
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_2d_mudpack( poisson, phi, rho )
    class(poisson_2d_mudpack_solver), target :: poisson
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
    external mudpack_cofx,mudpack_cofy,mudpack_bndsp
    external mudpack_cof,mudpack_cofcr
    !set initial guess because solve should be called every time step in a
    !time dependent problem and the elliptic operator does not depend on time.
    iguess = poisson%iguess

    !attempt solution
    intl = 1
    !write(*,106) intl,method,iguess


    select case (poisson%mudpack_case)
      case (SLL_SEPARABLE)
        if(associated(mudpack_wrapper))then
          print *,'#Problem mudpack_wrapper is not null()'
          stop
        endif
        mudpack_wrapper => poisson


        call mud2sp(iprm, &
          fprm, &
          poisson%work, &
          mudpack_cofx, &
          mudpack_cofy, &
          mudpack_bndsp, &
          rho, &
          phi, &
          poisson%mgopt, &
          error)
        !write(*,107) error
        if (error > 0) call exit(0)
        ! attempt to improve approximation to fourth order
        call mud24sp(poisson%work,phi,error)
        !write (*,108) error
        if (error > 0) call exit(0)
        
         mudpack_wrapper => null()
        
      case (SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS)
      if(associated(mudpack_wrapper))then
          print *,'#Problem mudpack_wrapper is not null()'
          stop
        endif
        mudpack_wrapper => poisson

        call mud2(iprm, &
          fprm, &
          poisson%work, &
          mudpack_cof, &
          mudpack_bndsp, &
          rho, &
          phi, &
          poisson%mgopt, &
          error)
        !write(*,107) error
        if (error > 0) call exit(0)
        ! attempt to improve approximation to fourth order
        call mud24(poisson%work,phi,error)
        !write (*,108) error
        if (error > 0) call exit(0)
        
         mudpack_wrapper => null()
        
      case (SLL_NON_SEPARABLE_WITH_CROSS_TERMS)
        if(associated(mudpack_wrapper))then
          print *,'#Problem mudpack_wrapper is not null()'
          stop
        endif
        mudpack_wrapper => poisson

        call mud2cr(iprm, &
          fprm, &
          poisson%work, &
          mudpack_cofcr, &
          mudpack_bndsp, &
          rho, &
          phi, &
          poisson%mgopt, &
          error)
        !write(*,107) error
        if (error > 0) call exit(0)
        ! attempt to improve approximation to fourth order
        call mud24cr(poisson%work, &
          mudpack_cofcr, &
          mudpack_bndsp, &
          phi, &
          error)
        !write (*,108) error
        if (error > 0) call exit(0)
        
         mudpack_wrapper => null()
         
      case default
        print *,'#bad mudpack_case',poisson%mudpack_case
        print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
        stop 
    end select



    
    !call solve( poisson%poiss, rho, phi)
    
  end subroutine compute_phi_from_rho_2d_mudpack

    ! solves E = -\nabla Phi in 2d
!    subroutine compute_E_from_phi_2d_fft( poisson, phi, E1, E2 )
!      class(poisson_2d_fft_solver) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine compute_E_from_phi_2d_fft

    ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
    subroutine compute_E_from_rho_2d_mudpack( poisson, E1, E2, rho )
      class(poisson_2d_mudpack_solver) :: poisson
      sll_real64,dimension(:,:),intent(in) :: rho
      sll_real64,dimension(:,:),intent(out) :: E1
      sll_real64,dimension(:,:),intent(out) :: E2
      
      print *,'#compute_E_from_rho_2d_mudpack'      
      print *,'#not implemented for the moment'

      E1 = 0._f64
      E2 = 0._f64
      print *,maxval(rho)
      
      if(.not.(associated(poisson%cxx_2d)))then
        print *,'#poisson%cxx_2d is not associated'
      endif

      stop
      
      !call solve( poisson%poiss, E1, E2, rho)
      
    end subroutine compute_E_from_rho_2d_mudpack
  
  
  
  
end module sll_module_poisson_2d_mudpack_solver

!> input x dependent coefficients
subroutine mudpack_cofx(x,cxx,cx,cex)
use sll_module_poisson_2d_mudpack_solver
implicit none
real(8)  :: x,cxx,cx,cex
cxx = mudpack_wrapper%cxx_1d_interp%interpolate_value(x)
cx  = mudpack_wrapper%cx_1d_interp%interpolate_value(x)
cex = mudpack_wrapper%cex_1d_interp%interpolate_value(x)
return
end

!> input y dependent coefficients
subroutine mudpack_cofy(y,cyy,cy,cey)
use sll_module_poisson_2d_mudpack_solver
implicit none
real(8)  :: y,cyy,cy,cey
cyy = mudpack_wrapper%cyy_1d_interp%interpolate_value(y)
cy  = mudpack_wrapper%cy_1d_interp%interpolate_value(y)
cey = mudpack_wrapper%cey_1d_interp%interpolate_value(y)
return
end

subroutine mudpack_cof(x,y,cxx,cyy,cx,cy,ce)
use sll_module_poisson_2d_mudpack_solver
implicit none
real(8)  :: x,cxx,cx
real(8)  :: y,cyy,cy,ce
cxx = mudpack_wrapper%cxx_2d_interp%interpolate_value(x,y)
cyy = mudpack_wrapper%cyy_2d_interp%interpolate_value(x,y)
cx  = mudpack_wrapper%cx_2d_interp%interpolate_value(x,y)
cy  = mudpack_wrapper%cy_2d_interp%interpolate_value(x,y)
ce  = mudpack_wrapper%ce_2d_interp%interpolate_value(x,y)
return
end

subroutine mudpack_cofcr(x,y,cxx,cxy,cyy,cx,cy,ce)
use sll_module_poisson_2d_mudpack_solver
implicit none
real(8)  :: x,cxx,cx,cxy
real(8)  :: y,cyy,cy,ce
cxx = mudpack_wrapper%cxx_2d_interp%interpolate_value(x,y)
cxy = mudpack_wrapper%cxy_2d_interp%interpolate_value(x,y)
cyy = mudpack_wrapper%cyy_2d_interp%interpolate_value(x,y)
cx  = mudpack_wrapper%cx_2d_interp%interpolate_value(x,y)
cy  = mudpack_wrapper%cy_2d_interp%interpolate_value(x,y)
ce  = mudpack_wrapper%ce_2d_interp%interpolate_value(x,y)
return
end
!> input mixed derivative b.c. to mud2sp
subroutine mudpack_bndsp(kbdy,xory,alfa,gbdy)
use sll_module_poisson_2d_mudpack_solver
implicit none
integer  :: kbdy
real(8)  :: xory,alfa,gbdy,x,y,pe,px,py
real(8)  :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

pe = 0.0
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















