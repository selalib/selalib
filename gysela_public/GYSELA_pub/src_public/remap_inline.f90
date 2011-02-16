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
!------------------------------------------------------
! file : remap_inline.f90
! date : 27/04/2010
!------------------------------------------------------
      
!-----------------------------------------------------
! Transform cylindrical coordinates (r,theta) 
! in cartesian coordinates (x,y) for one point
!-----------------------------------------------------
subroutine cyl2cart(r,theta,x,y)
  real(RKIND), intent(in)    :: r
  real(RKIND), intent(in)    :: theta
  real(RKIND), intent(inout) :: x
  real(RKIND), intent(inout) :: y
      
  x = r*cos(theta)
  y = r*sin(theta)
end subroutine cyl2cart
      
!-----------------------------------------------------
! Transform cartesian coordinates (x,y)
! in cylindrical coordinates (r,theta) :
!   r     = sqrt(x**2+y**2)
!   theta = arctan(y/x)
!-----------------------------------------------------
subroutine cart2cyl(x,y,r,theta)
  real(RKIND), intent(in)  :: x
  real(RKIND), intent(in)  :: y
  real(RKIND), intent(out) :: r
  real(RKIND), intent(out) :: theta
      
  r     = sqrt(x*x+y*y)
  theta = modulo(atan2(y,x),TWOPI)
end subroutine cart2cyl
      
!----------------------------------------------------------
! Verify r coordinate, i.e.
!  . if r<rmin set r=rmin 
!  . if r>rmax set r=rmax 
!----------------------------------------------------------
subroutine r_verif(rmin,rmax,r)
  real(RKIND), intent(in)    :: rmin, rmax
  real(RKIND), intent(inout) :: r
      
  !*** r verification ***
  if (r.lt.rmin) then
    r = rmin
  else 
    if (r.gt.rmax) then
      r = rmax
    end if
  end if
end subroutine r_verif
