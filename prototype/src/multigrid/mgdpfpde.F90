module mgdpfpde
#include "sll_working_precision.h"
use mpi
implicit none 

contains

subroutine mgdpfpde_3d(sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,cof,r,xl,yl,zl)
integer :: sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf
real(8) :: cof(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1,8)
real(8) :: r(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1),xl,yl,zl
!------------------------------------------------------------------------
! Determine coefficients for the pressure equation at the finest grid
! level. These coefficients involve densities half-way between the
! pressure and density nodes.
!
! cof array: 1 -> i-1
!            2 -> i+1
!            3 -> j-1
!            4 -> j+1
!            5 -> k-1
!            6 -> k+1
!            7 -> central
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : --
!------------------------------------------------------------------------
real(8) :: dlx,todlxx,dly,todlyy,dlz,todlzz,rijk,c1,c2,c3,c4,c5,c6
integer :: i,j,k

dlx=xl/float(nxf-1)
todlxx=2.0d0/(dlx*dlx)
dly=yl/float(nyf-1)
todlyy=2.0d0/(dly*dly)
dlz=zl/float(nzf-1)
todlzz=2.0d0/(dlz*dlz)
do k=szf,ezf
  do j=syf,eyf
    do i=sxf,exf
      rijk=r(i,j,k)
      c1=todlxx/(rijk+r(i-1,j,k))
      c2=todlxx/(rijk+r(i+1,j,k))
      c3=todlyy/(rijk+r(i,j-1,k))
      c4=todlyy/(rijk+r(i,j+1,k))
      c5=todlzz/(rijk+r(i,j,k-1))
      c6=todlzz/(rijk+r(i,j,k+1))
      cof(i,j,k,1)=c1
      cof(i,j,k,2)=c2
      cof(i,j,k,3)=c3
      cof(i,j,k,4)=c4
      cof(i,j,k,5)=c5
      cof(i,j,k,6)=c6
      cof(i,j,k,7)=-(c1+c2+c3+c4+c5+c6)
    end do
  end do
end do

end subroutine


subroutine mgdpfpde_2d(sxf,exf,syf,eyf,nxf,nyf,cof,r,xl,yl,bd)
# include "mgd2.h"

integer sxf,exf,syf,eyf,nxf,nyf,bd(8)
REALN cof(sxf-1:exf+1,syf-1:eyf+1,6)
REALN r(sxf-1:exf+1,syf-1:eyf+1),xl,yl
!------------------------------------------------------------------------
! Determine coefficients for the pressure equation at the finest grid
! level. These coefficients involve densities half-way between the
! pressure and density nodes. Works for periodic, Neumann, and
! Dirichlet boundary conditions.
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
REALN dlx,todlxx,dly,todlyy,rij
integer i,j
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

end module mgdpfpde
