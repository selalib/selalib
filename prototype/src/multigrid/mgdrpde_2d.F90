!------------------------------------------------------------------------
! Discretize the pde: set the coefficients of the cof matrix. Works
! for periodic, Neumann, and Dirichlet boundary conditions.
!
! cof array:
!
!         cof(4)
!           |
!           |
! cof(1)--cof(5)--cof(2)
!           |
!           |
!         cof(3)
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : --
!------------------------------------------------------------------------
subroutine mgdrpde(sxm,exm,sym,eym,nxm,nym,cof,xl,yl,bd)
#include "sll_working_precision.h"
#include "mgd2.h"
sll_int32  :: sxm,exm,sym,eym,nxm,nym,bd(8)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6),xl,yl
sll_real64 :: dlx,odlxx,dly,odlyy
sll_int32  :: i,j

! calculate off-diagonal terms
dlx=xl/float(nxm-1)
odlxx=1.0d0/(dlx*dlx)
dly=yl/float(nym-1)
odlyy=1.0d0/(dly*dly)

#ifdef POLAR

do j=sym,eym
  do i=sxm,exm
    cof(i,j,1)=odlxx
    cof(i,j,2)=odlxx
    cof(i,j,3)=odlyy
    cof(i,j,4)=odlyy
  end do
end do

#else

do j=sym,eym
  do i=sxm,exm
    cof(i,j,1)=odlxx
    cof(i,j,2)=odlxx
    cof(i,j,3)=odlyy
    cof(i,j,4)=odlyy
  end do
end do

#endif
!
! enforce Neumann BCs
!
if (bd(1).eq.1) then
  do j=sym,eym
    cof(exm,j,2)=0.0d0
  end do
end if
if (bd(5).eq.1) then
  do j=sym,eym
    cof(sxm,j,1)=0.0d0
  end do
end if
if (bd(3).eq.1) then
  do i=sxm,exm
    cof(i,sym,3)=0.0d0
  end do
end if
if (bd(7).eq.1) then
  do i=sxm,exm
    cof(i,eym,4)=0.0d0
  end do
end if
!
! calculate diagonal term
!
do j=sym,eym
  do i=sxm,exm
    cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
  end do
end do

end subroutine
