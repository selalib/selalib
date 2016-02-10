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

!> @ingroup poisson_solvers
!> @brief
!> Solves Poisson equation on 2d curvilinear mesh
!> @details
!> We use mudpack library for multigrid method.
!>
!> <b> Equations </b>
!>
!> \f[ \sum_{i,j=1}^2 A_{i,j}\partial_{i,j} \phi(x,y)
!>        +  \sum_{i=1}^2B_i\partial_i \phi(x,y)
!>        +C \phi(x,y) = \rho(x,y)
!> \f]
!> in polar coordinates. This leads when 
!> \f[ 
!>      A_{1,2}=A_{2,1}=0 \\
!> \f]
!> \f[ 
!>      B_2 = 0 \\
!> \f]
!> \f[ 
!>      A_{1,1}\partial_{1,1}\hat{\phi}+B_1\partial_{1}\hat{\phi}+(C+A_{2,2}k^2)\hat{\phi} = \hat{\rho}
!> \f]
module sll_m_poisson_2d_mudpack_curvilinear
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

! use F77_mudpack, only: &
!   mud2, &
!   mud24, &
!   mud24cr, &
!   mud24sp, &
!   mud2cr, &
!   mud2sp

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_periodic

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_f_new_cubic_spline_interpolator_1d

  use sll_m_cubic_spline_interpolator_2d, only: &
    sll_f_new_cubic_spline_interpolator_2d

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  use sll_m_mudpack_curvilinear, only: &
    sll_p_non_separable_with_cross_terms, &
    sll_p_non_separable_without_cross_terms, &
    sll_p_separable

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base

  implicit none

  public :: &
    sll_f_new_poisson_2d_mudpack_curvilinear

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Derived type to solve Poisson equation on 2d curvilinear mesh
  type, extends(sll_c_poisson_2d_base) :: poisson_2d_mudpack_curvilinear
  
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:,:), pointer :: cxx_2d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:,:), pointer :: cxy_2d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:,:), pointer :: cyy_2d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:,:), pointer :: cx_2d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:,:), pointer :: cy_2d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:,:), pointer :: ce_2d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:), pointer :: cxx_1d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:), pointer :: cyy_1d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:), pointer :: cx_1d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:), pointer :: cy_1d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:), pointer :: cex_1d
    !> PLEASE ADD DOCUMENTATION
    sll_real64, dimension(:), pointer :: cey_1d
    !> PLEASE ADD DOCUMENTATION
    sll_real64 :: cxx
    !> PLEASE ADD DOCUMENTATION
    sll_real64 :: cyy
    !> PLEASE ADD DOCUMENTATION
    sll_real64 :: cx
    !> PLEASE ADD DOCUMENTATION
    sll_real64 :: cy
    !> PLEASE ADD DOCUMENTATION
    sll_real64 :: ce
    !> PLEASE ADD DOCUMENTATION
    sll_int32  :: mudpack_case
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_2d), pointer   :: cxx_2d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_2d), pointer   :: cxy_2d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_2d), pointer   :: cyy_2d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_2d), pointer   :: cx_2d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_2d), pointer   :: cy_2d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_2d), pointer   :: ce_2d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_1d), pointer   :: cxx_1d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_1d), pointer   :: cyy_1d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_1d), pointer   :: cx_1d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_1d), pointer   :: cy_1d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_1d), pointer   :: cex_1d_interp
    !> PLEASE ADD DOCUMENTATION
    class(sll_c_interpolator_1d), pointer   :: cey_1d_interp

    sll_real64, dimension(:), pointer :: work !< array for tmp data
    sll_int32  :: mgopt(4) !< Option to control multigrid
    sll_int32  :: iprm(16) !< Indices to control grid sizes
    sll_real64 :: fprm(6)  !< Real to set boundary conditions
    sll_int32  :: iguess   !< Initial solution or loop over time

  contains

    !> PLEASE ADD DOCUMENTATION
    procedure, pass(poisson) :: initialize => initialize_poisson_2d_mudpack_curvilinear
    !> PLEASE ADD DOCUMENTATION
    procedure, pass(poisson) :: compute_phi_from_rho => compute_phi_from_rho_2d_mudpack
    !> PLEASE ADD DOCUMENTATION
    procedure, pass(poisson) :: compute_E_from_rho => compute_E_from_rho_2d_mudpack
      
  end type poisson_2d_mudpack_curvilinear

  !> PLEASE ADD DOCUMENTATION
  class(poisson_2d_mudpack_curvilinear), pointer :: mudpack_wrapper => null()


contains

  !> PLEASE ADD DOCUMENTATION
  function sll_f_new_poisson_2d_mudpack_curvilinear( &
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
      
    type(poisson_2d_mudpack_curvilinear),pointer :: poisson
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
    call initialize_poisson_2d_mudpack_curvilinear( &
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
    
  end function sll_f_new_poisson_2d_mudpack_curvilinear
  
  
  subroutine initialize_poisson_2d_mudpack_curvilinear( &
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
    class(poisson_2d_mudpack_curvilinear), target :: poisson
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
    !!!! end variables for mudpack 

    nx = nc_eta1+1
    ny = nc_eta2+1
    allocate(phi(nx*ny))
    allocate(rhs(nx*ny))
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
    tolmax = 0.0_f64

!    write(*,101) (iprm(i),i=1,15)
!    write(*,102) (poisson%mgopt(i),i=1,4)
!    write(*,103) xa,xb,yc,yd,tolmax
!    write(*,104) intl

!call mud2sp(iprm,fprm,self%work,cofx,cofy,bndsp,rhs,phi,self%mgopt,error)

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
      case (sll_p_separable)
        if(present(cxx_2d).or.present(cxy_2d).or.present(cyy_2d)&
          .or.present(cx_2d).or.present(cy_2d).or.present(ce_2d)) then
          print *,'#2d arrays should not be here'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
          stop
        endif
        
        if((.not.(present(cxx_1d))).and.(.not.(present(cxx)))) then
          print *,'#1d/0d array should be here for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
          stop        
        endif
        if(present(cxx_1d).and.present(cxx))then
          print *,'#please choose between 1d or 0d array for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
          stop        
        endif
        if(present(cyy_1d).and.present(cyy))then
          print *,'#please choose between 1d or 0d array for cyy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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

        poisson%cxx_1d_interp => sll_f_new_cubic_spline_interpolator_1d( &
          nx, &
          eta1_min, &
          eta1_max, &
          sll_p_periodic)          
        call poisson%cxx_1d_interp%compute_interpolants( poisson%cxx_1d )          

        poisson%cyy_1d_interp => sll_f_new_cubic_spline_interpolator_1d( &
          ny, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic)          
        call poisson%cyy_1d_interp%compute_interpolants( poisson%cyy_1d )          

        poisson%cx_1d_interp => sll_f_new_cubic_spline_interpolator_1d( &
          nx, &
          eta1_min, &
          eta1_max, &
          sll_p_periodic)          
        call poisson%cx_1d_interp%compute_interpolants( poisson%cx_1d )          

        poisson%cy_1d_interp => sll_f_new_cubic_spline_interpolator_1d( &
          ny, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic)          
        call poisson%cy_1d_interp%compute_interpolants( poisson%cy_1d )          

        poisson%cex_1d_interp => sll_f_new_cubic_spline_interpolator_1d( &
          nx, &
          eta1_min, &
          eta1_max, &
          sll_p_periodic)          
        call poisson%cex_1d_interp%compute_interpolants( poisson%cex_1d )          

        poisson%cey_1d_interp => sll_f_new_cubic_spline_interpolator_1d( &
          ny, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic)          
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

                
        
      case (sll_p_non_separable_without_cross_terms)
      
        if(present(cxx_1d).or.present(cyy_1d).or.present(cx_1d)& 
          .or.present(cy_1d).or.present(cex_1d).or.present(cey_1d)) then
          print *,'#1d arrays should not be here'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
          stop
        endif
        
        if((.not.(present(cxx_2d))).and.(.not.(present(cxx)))) then
          print *,'#2d/0d array should be here for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
          stop        
        endif
        if(present(cxx_2d).and.present(cxx))then
          print *,'#please choose between 2d or 0d array for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
          stop        
        endif
        if(present(cyy_2d).and.present(cyy))then
          print *,'#please choose between 2d or 0d array for cyy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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


        poisson%cxx_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
        call poisson%cxx_2d_interp%compute_interpolants( poisson%cxx_2d )          

        poisson%cyy_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
        call poisson%cyy_2d_interp%compute_interpolants( poisson%cyy_2d )          

        poisson%cx_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
        call poisson%cx_2d_interp%compute_interpolants( poisson%cx_2d )          

        poisson%cy_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
        call poisson%cy_2d_interp%compute_interpolants( poisson%cy_2d )          

        poisson%ce_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
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

      case (sll_p_non_separable_with_cross_terms)
      
        if(present(cxx_1d).or.present(cyy_1d).or.present(cx_1d)& 
          .or.present(cy_1d).or.present(cex_1d).or.present(cey_1d)) then
          print *,'#1d arrays should not be here'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
          stop
        endif
        
        if((.not.(present(cxx_2d))).and.(.not.(present(cxx)))) then
          print *,'#2d/0d array should be here for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
          stop        
        endif
        if(present(cxx_2d).and.present(cxx))then
          print *,'#please choose between 2d or 0d array for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
          stop        
        endif
        if(present(cxy_2d).and.present(cxy))then
          print *,'#please choose between 2d or 0d array for cxy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
          stop        
        endif
        if(present(cyy_2d).and.present(cyy))then
          print *,'#please choose between 2d or 0d array for cyy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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
          print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
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


        poisson%cxx_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
        call poisson%cxx_2d_interp%compute_interpolants( poisson%cxx_2d )   
               
        poisson%cxy_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
        call poisson%cxy_2d_interp%compute_interpolants( poisson%cxy_2d ) 
        
        poisson%cyy_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
        call poisson%cyy_2d_interp%compute_interpolants( poisson%cyy_2d )          

        poisson%cx_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
        call poisson%cx_2d_interp%compute_interpolants( poisson%cx_2d )          

        poisson%cy_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
        call poisson%cy_2d_interp%compute_interpolants( poisson%cy_2d )          

        poisson%ce_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          sll_p_periodic, &
          sll_p_periodic)    
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
        print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
        stop 
    end select

  end subroutine initialize_poisson_2d_mudpack_curvilinear
  
  ! solves \Delta phi = -rho in 2d
  subroutine compute_phi_from_rho_2d_mudpack( poisson, phi, rho )
    class(poisson_2d_mudpack_curvilinear), target :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    !sll_real64        :: phi(:,:)  !< Electric potential
    !sll_real64        :: rhs(:,:)  !< Charge density
    !put integer and floating point argument names in contiguous
    !storeage for labelling in vectors iprm,fprm
    sll_int32  :: iprm(16)
    sll_real64 :: fprm(6)
    sll_int32  :: error
    sll_int32  :: i1
    sll_int32  :: i2
    sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
    sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
    common/itmud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
    sll_real64 :: xa,xb,yc,yd,tolmax,relmax
    common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

    equivalence(intl,iprm)
    equivalence(xa,fprm)

    !set initial guess because solve should be called every time step in a
    !time dependent problem and the elliptic operator does not depend on time.
    iguess = poisson%iguess

    !attempt solution
    intl = 1
    !write(*,106) intl,method,iguess
        
    if(nxa == sll_p_dirichlet) then
       do i2=1,ny
          phi(1,i2) = 0._f64
       end do
    endif
    if(nxb == sll_p_dirichlet) then
       do i2=1,ny
          phi(nx,i2) = 0._f64
       end do
    endif
    if(nyc == sll_p_dirichlet) then
       do i1=1,nx
          phi(i1,1) = 0._f64
       end do
    endif
    if(nyd == sll_p_dirichlet) then
       do i1=1,nx
          phi(i1,ny) = 0._f64
       end do
    endif 
    select case (poisson%mudpack_case)
      case (sll_p_separable)
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
          -rho, &
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
        
      case (sll_p_non_separable_without_cross_terms)
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
          -rho, &
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
        
      case (sll_p_non_separable_with_cross_terms)
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
          -rho, &
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
        print *,'#in subroutine initialize_poisson_2d_mudpack_curvilinear'
        stop 
    end select

  end subroutine compute_phi_from_rho_2d_mudpack

  ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
  subroutine compute_E_from_rho_2d_mudpack( poisson, E1, E2, rho )
    class(poisson_2d_mudpack_curvilinear) :: poisson
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
  
  !> input x dependent coefficients
  subroutine mudpack_cofx(x,cxx,cx,cex)
    real(8)  :: x,cxx,cx,cex
    cxx = mudpack_wrapper%cxx_1d_interp%interpolate_from_interpolant_value(x)
    cx  = mudpack_wrapper%cx_1d_interp%interpolate_from_interpolant_value(x)
    cex = mudpack_wrapper%cex_1d_interp%interpolate_from_interpolant_value(x)
  end subroutine mudpack_cofx

  !> input y dependent coefficients
  subroutine mudpack_cofy(y,cyy,cy,cey)
    real(8)  :: y,cyy,cy,cey
    cyy = mudpack_wrapper%cyy_1d_interp%interpolate_from_interpolant_value(y)
    cy  = mudpack_wrapper%cy_1d_interp%interpolate_from_interpolant_value(y)
    cey = mudpack_wrapper%cey_1d_interp%interpolate_from_interpolant_value(y)
  end subroutine mudpack_cofy

  subroutine mudpack_cof(x,y,cxx,cyy,cx,cy,ce)
    real(8)  :: x,cxx,cx
    real(8)  :: y,cyy,cy,ce
    cxx = mudpack_wrapper%cxx_2d_interp%interpolate_from_interpolant_value(x,y)
    cyy = mudpack_wrapper%cyy_2d_interp%interpolate_from_interpolant_value(x,y)
    cx  = mudpack_wrapper%cx_2d_interp%interpolate_from_interpolant_value(x,y)
    cy  = mudpack_wrapper%cy_2d_interp%interpolate_from_interpolant_value(x,y)
    ce  = mudpack_wrapper%ce_2d_interp%interpolate_from_interpolant_value(x,y)
  end subroutine mudpack_cof

  subroutine mudpack_cofcr(x,y,cxx,cxy,cyy,cx,cy,ce)
    real(8)  :: x,cxx,cx,cxy
    real(8)  :: y,cyy,cy,ce
    cxx = mudpack_wrapper%cxx_2d_interp%interpolate_from_interpolant_value(x,y)
    cxy = mudpack_wrapper%cxy_2d_interp%interpolate_from_interpolant_value(x,y)
    cyy = mudpack_wrapper%cyy_2d_interp%interpolate_from_interpolant_value(x,y)
    cx  = mudpack_wrapper%cx_2d_interp%interpolate_from_interpolant_value(x,y)
    cy  = mudpack_wrapper%cy_2d_interp%interpolate_from_interpolant_value(x,y)
    ce  = mudpack_wrapper%ce_2d_interp%interpolate_from_interpolant_value(x,y)
  end subroutine mudpack_cofcr

  !> input mixed derivative b.c. to mud2sp
  subroutine mudpack_bndsp(kbdy,xory,alfa,gbdy)
    integer  :: kbdy
    real(8)  :: xory,alfa,gbdy,x,y,pe,px,py
    real(8)  :: xa,xb,yc,yd,tolmax,relmax
    common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

    pe = 0.0_8
    !subroutine not used in periodic case
    if (kbdy == 1) then  ! x=xa boundary
       y = xory
       x = xa
       alfa = -1.0_8
       gbdy = px + alfa*pe
       return
    end if

    if (kbdy == 4) then  ! y=yd boundary
       y = yd
       x = xory
       alfa = 1.0_8
       gbdy = py + alfa*pe
       return
    end if
  end subroutine mudpack_bndsp

end module sll_m_poisson_2d_mudpack_curvilinear
