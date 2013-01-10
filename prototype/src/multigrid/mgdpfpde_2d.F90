!------------------------------------------------------------------------
!> Determine coefficients for the pressure equation at the finest grid
!> level. These coefficients involve densities half-way between the
!> pressure and density nodes. Works for periodic, Neumann, and
!> Dirichlet boundary conditions.
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
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
!------------------------------------------------------------------------
subroutine mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,cof,r,xl,yl,bd)
#include "sll_working_precision.h"
#include "mgd2.h"

sll_int32  :: sxf,exf,syf,eyf,nxf,nyf,bd(8)
sll_real64 :: cof(sxf-1:exf+1,syf-1:eyf+1,6)
sll_real64 :: r(sxf-1:exf+1,syf-1:eyf+1),xl,yl
sll_real64 :: dlx,todlxx,dly,todlyy,rij
sll_int32  :: i,j
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif

dlx=xl/float(nxf-1)
todlxx=2.0d0/(dlx*dlx)
dly=yl/float(nyf-1)
todlyy=2.0d0/(dly*dly)
do j=syf,eyf
  do i=sxf,exf
    rij=r(i,j)
    cof(i,j,1)=todlxx/(rij+r(i-1,j))
    cof(i,j,2)=todlxx/(rij+r(i+1,j))
    cof(i,j,3)=todlyy/(rij+r(i,j-1))
    cof(i,j,4)=todlyy/(rij+r(i,j+1))
  end do
end do

! enforce wall BCs 

if (bd(1).eq.1) then
  do j=syf,eyf
    cof(exf,j,2)=0.0d0
  end do
end if
if (bd(5).eq.1) then
  do j=syf,eyf
    cof(sxf,j,1)=0.0d0
  end do
end if
if (bd(3).eq.1) then
  do i=sxf,exf
    cof(i,syf,3)=0.0d0
  end do
end if
if (bd(7).eq.1) then
  do i=sxf,exf
    cof(i,eyf,4)=0.0d0
  end do
end if

! calculate diagonal term

do j=syf,eyf
  do i=sxf,exf
    cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
  end do
end do

# if cdebug
timing(84)=timing(84)+MPI_WTIME()-tinitial
# endif

end subroutine
