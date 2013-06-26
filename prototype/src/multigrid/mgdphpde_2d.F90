!> For the new version of the multigrid code, determine the coefficients
!> for the pressure equation at all grid levels except the finest one.
!> The coefficients are determined directly from the density array r
!> through some manipulation of indices and are values at (i+1/2,j+1/2)
!> points. Works for periodic, Neumann, and Dirichlet boundary
!> conditions.
!>
!> cof array:
!>
!>         cof(4)
!>           |
!>           |
!> cof(1)--cof(5)--cof(2)
!>           |
!>           |
!>         cof(3)
!>
!> Code      : mgd2, 2-D parallel multigrid solver
!> Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
subroutine mgdphpde(sxm,exm,sym,eym,nxm,nym,cof,         &
                    sx,ex,sy,ey,nxf,nyf,r,bd,xl,yl)

#include "sll_working_precision.h"
#include "mgd2.h"

sll_int32  :: sxm,exm,sym,eym,nxm,nym,sx,ex,sy,ey,nxf,nyf,bd(8)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6)
sll_real64 :: r(sx-1:ex+1,sy-1:ey+1),xl,yl
sll_real64 :: dlx,fodlxx,dly,fodlyy
sll_int32  :: i,j,im,jm,is,js,istep,jstep

! calculate off-diagonal terms

dlx=xl/float(nxm-1)
dly=yl/float(nym-1)


istep=(nxf-1)/(nxm-1)
jstep=(nyf-1)/(nym-1)

#ifdef POLAR

fodlxx=1.0d0/(dlx*dlx)
fodlyy=1.0d0/(dly*dly)

do j=sym,eym
  jm=2*jstep*j-3*(jstep-1)
  do i=sxm,exm
    im=2*istep*i-3*(istep-1)
    is=(im-istep)/2
    js=jm/2
    cof(i,j,1)=fodlxx+1./((r(is,js)+r(is+1,js))*dlx)
    is=(im+istep)/2
    cof(i,j,2)=fodlxx-1./((r(is,js)+r(is+1,js))*dlx)
    is=im/2
    js=(jm-jstep)/2
    cof(i,j,3)=fodlyy/(r(is,js)*r(is,js))
    js=(jm+jstep)/2
    cof(i,j,4)=fodlyy/(r(is,js)*r(is,js))
  end do
end do

#else

fodlxx=4.0d0/(dlx*dlx)
fodlyy=4.0d0/(dly*dly)

do j=sym,eym
  jm=2*jstep*j-3*(jstep-1)
  do i=sxm,exm
    im=2*istep*i-3*(istep-1)
    is=(im-istep)/2
    js=jm/2
    cof(i,j,1)=fodlxx/(r(is,js)+r(is+1,js)+r(is,js+1)+r(is+1,js+1))
    is=(im+istep)/2
    cof(i,j,2)=fodlxx/(r(is,js)+r(is+1,js)+r(is,js+1)+r(is+1,js+1))
    is=im/2
    js=(jm-jstep)/2
    cof(i,j,3)=fodlyy/(r(is,js)+r(is+1,js)+r(is,js+1)+r(is+1,js+1))
    js=(jm+jstep)/2
    cof(i,j,4)=fodlyy/(r(is,js)+r(is+1,js)+r(is,js+1)+r(is+1,js+1))
  end do
end do

#endif

! enforce wall BCs

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

! calculate diagonal term

do j=sym,eym
  do i=sxm,exm
    cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
  end do
end do

end subroutine
