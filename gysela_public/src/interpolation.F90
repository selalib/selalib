!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
      
!-------------------------------------------------------
! file : interpolation.f90
! date : 02/11/2000
! - linear interpolation in 3D 
! - interpolation in 4D by using spline
!   coefficients 
!-------------------------------------------------------
module interpolation_module
  use prec_const
  use geometry_class
  use MPIutils_module
  use spline1d_class
  use utils_module, only : locate
  implicit none
      
  !******************************
  contains
  !******************************
      
  !--------------------------------------------- 
  !  function interpolation by using the 
  !  2D-spline coefficients calculated 
  !  in r and theta directions
  !   - coef2d contains the 2D cubic splines 
  !  in r and theta
  ! (f must be divided according to (r,theta) direction)
  !---------------------------------------------
  subroutine interpol2d_rtheta(geom,coef2d,x1,x2, &
    xstar1,xstar2,finterpol)
    type(geometry)                     , intent(in)    :: geom
    real(RKIND)   , dimension(0:)      , intent(in)    :: x1, x2
    real(RKIND)   , dimension(-1:,-1:) , intent(in)    :: coef2d
    real(RKIND)                        , intent(in)    :: xstar1
    real(RKIND)                        , intent(in)    :: xstar2
    real(RKIND)                        , intent(inout) :: finterpol
    
    integer                      :: n1, n2
    integer                      :: ipos1, ipos2
    integer                      :: i, j
    real(RKIND)                  :: h1, h2
    real(RKIND), dimension(-1:2) :: sbase1, sbase2
    
    !*** local variable initialisation ***
    n1 = geom%Nr
    n2 = geom%Ntheta
    h1 = geom%dr
    h2 = geom%dtheta
    
    !*** array position location ***
    call locate(xstar1,x1,n1,h1, &
      " interpol2d_rtheta 1 "//char(0),ipos1)
    call locate(xstar2,x2,n2,h2, &
      " interpol2d_rtheta 2 "//char(0),ipos2)
      
    !*** calculation of cubic spline basis ***
    call spline_basis(x1(ipos1),xstar1,x1(ipos1+1),h1,sbase1)
    call spline_basis(x2(ipos2),xstar2,x2(ipos2+1),h2,sbase2)
      
    !*** computation of f(x1*,x2*,k,l) ***
    finterpol = 0._RKIND
    finterpol = finterpol + &
      coef2d(ipos1-1,ipos2-1)*sbase1(-1)*sbase2(-1) + &
      coef2d(ipos1,ipos2-1)*sbase1(0)*sbase2(-1) + &
      coef2d(ipos1+1,ipos2-1)*sbase1(1)*sbase2(-1) + &
      coef2d(ipos1+2,ipos2-1)*sbase1(2)*sbase2(-1)
    finterpol = finterpol + &
      coef2d(ipos1-1,ipos2)*sbase1(-1)*sbase2(0) + &
      coef2d(ipos1,ipos2)*sbase1(0)*sbase2(0) + &
      coef2d(ipos1+1,ipos2)*sbase1(1)*sbase2(0) + &
      coef2d(ipos1+2,ipos2)*sbase1(2)*sbase2(0)
    finterpol = finterpol + &
      coef2d(ipos1-1,ipos2+1)*sbase1(-1)*sbase2(1) + &
      coef2d(ipos1,ipos2+1)*sbase1(0)*sbase2(1) + &
      coef2d(ipos1+1,ipos2+1)*sbase1(1)*sbase2(1) + &
      coef2d(ipos1+2,ipos2+1)*sbase1(2)*sbase2(1)
    finterpol = finterpol + &
      coef2d(ipos1-1,ipos2+2)*sbase1(-1)*sbase2(2) + &
      coef2d(ipos1,ipos2+2)*sbase1(0)*sbase2(2) + &
      coef2d(ipos1+1,ipos2+2)*sbase1(1)*sbase2(2) + &
      coef2d(ipos1+2,ipos2+2)*sbase1(2)*sbase2(2)
  end subroutine interpol2d_rtheta
      
  !--------------------------------------------- 
  !  function interpolation by using the 
  !  1D-spline coefficients calculated 
  !  in phi direction
  !   - coef1d contains the 1D cubic splines 
  !  in phi direction
  ! (f must be divided according to (r,theta) direction)
  !---------------------------------------------
  subroutine interpol1d_phi(geom,coef1d,x3,xstar3,finterpol)
    type(geometry)                , intent(in)  :: geom
    real(RKIND)   , dimension(:)  , pointer     :: coef1d
    real(RKIND)   , dimension(:)  , pointer     :: x3
    real(RKIND)                   , intent(in)  :: xstar3
    real(RKIND)                   , intent(out) :: finterpol
    
    integer                      :: n3
    integer                      :: ipos3
    integer                      :: k
    real(RKIND)                  :: h3
    real(RKIND), dimension(-1:2) :: sbase3
    
    !*** local variable initialisation ***
    n3 = geom%Nphi
    h3 = geom%dphi
    
    !*** array position location ***
    call locate(xstar3,x3,n3,h3, &
      "interpol1d_phi"//char(0),ipos3)
      
    !*** calculation of cubic spline basis ***
    call spline_basis(x3(ipos3),xstar3,x3(ipos3+1),h3,sbase3)
      
    !*** computation of f(i,j,x3*,l) ***
    finterpol = 0._RKIND
    finterpol = finterpol + coef1d(ipos3-1)*sbase3(-1) + &
      coef1d(ipos3)*sbase3(0) + coef1d(ipos3+1)*sbase3(1) + &
      coef1d(ipos3+2)*sbase3(2)
  end subroutine interpol1d_phi
      
  !--------------------------------------------- 
  !  function interpolation by using the 
  !  1D-spline coefficients calculated 
  !  in vpar direction
  !   - coef1d contains the 1D cubic splines 
  !  in vpar direction
  ! (f must be divided according to (r,theta) direction)
  !---------------------------------------------
  subroutine interpol1d_vpar(geom,coef1d,x4,xstar4,finterpol)
    type(geometry)                , intent(in)  :: geom
    real(RKIND)   , dimension(:)  , pointer     :: coef1d
    real(RKIND)   , dimension(:)  , pointer     :: x4
    real(RKIND)                   , intent(in)  :: xstar4
    real(RKIND)                   , intent(out) :: finterpol
    
    integer                      :: n4
    integer                      :: ipos4
    integer                      :: l
    real(RKIND)                  :: h4
    real(RKIND), dimension(-1:2) :: sbase4
    
    !*** local variable initialisation ***
    n4 = geom%Nvpar
    h4 = geom%dvpar
    
    !*** array position location ***
    call locate(xstar4,x4,n4,h4, &
      "interpol1d_vpar"//char(0),ipos4)
      
    !*** calculation of cubic spline basis ***
    call spline_basis(x4(ipos4),xstar4,x4(ipos4+1),h4,sbase4)
      
    !*** computation of f(i,j,k,x4*) ***
    finterpol = 0._RKIND
    finterpol = finterpol + coef1d(ipos4-1)*sbase4(-1) + &
      coef1d(ipos4)*sbase4(0) + coef1d(ipos4+1)*sbase4(1) + &
      coef1d(ipos4+2)*sbase4(2)
  end subroutine interpol1d_vpar
end module interpolation_module
