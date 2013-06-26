!> For the old version of the multigrid code, determine coefficients 
!> for the pressure equation at all grid levels but the finest one.
!> The coefficients are determined from the values of the density
!> at integer nodes (i,j). Works only for periodic boundary conditions.
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
!> Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
subroutine mgdppde(sxm,exm,sym,eym,nxm,nym,cof,     &
                   sxf,exf,syf,eyf,rf,xl,yl,bd)

#include "sll_working_precision.h"

sll_int32  :: sxm,exm,sym,eym,nxm,nym,sxf,exf,syf,eyf,bd(8)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6)
sll_real64 :: rf(sxf-1:exf+1,syf-1:eyf+1),xl,yl
sll_real64 :: dlx,odlxx,dly,odlyy
sll_int32  :: i,j,is,js
sll_real64 :: c1, c2, c3, c4

! calculate off-diagonal terms
dlx=xl/float(nxm-1)
odlxx=1.0d0/(dlx*dlx)
dly=yl/float(nym-1)
odlyy=1.0d0/(dly*dly)
do j=sym,eym
  js=2*j-1
  do i=sxm,exm
    is=2*i-1
    C1=odlxx/rf(is-1,js)
    C2=odlxx/rf(is+1,js)
    C3=odlyy/rf(is,js-1)
    C4=odlyy/rf(is,js+1)
    cof(i,j,1)=C1
    cof(i,j,2)=C2
    cof(i,j,3)=C3
    cof(i,j,4)=C4
    cof(i,j,5)=-(C1+C2+C3+C4)
  end do
end do

end subroutine
