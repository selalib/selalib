module mgdppde
#include "sll_working_precision.h"
use mpi
implicit none 

contains

subroutine mgdppde_3d(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,cof, &
                   sxf,exf,syf,eyf,szf,ezf,rf,xl,yl,zl)


integer :: sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm 
integer :: sxf,exf,syf,eyf,szf,ezf
real(8) :: cof(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1,8)
real(8) :: rf(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1),xl,yl,zl
!------------------------------------------------------------------------
! Determine coefficients for the pressure equation at all grid levels
! but the finest one (this is done in mgdpfpde).
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
real(8) :: dlx,odlxx,dly,odlyy,dlz,odlzz,c1,c2,c3,c4,c5,c6
integer :: i,j,k,is,js,ks

dlx=xl/float(nxm-1)
odlxx=1.0d0/(dlx*dlx)
dly=yl/float(nym-1)
odlyy=1.0d0/(dly*dly)
dlz=zl/float(nzm-1)
odlzz=1.0d0/(dlz*dlz)
do k=szm,ezm
  ks=2*k-1
  do j=sym,eym
    js=2*j-1
    do i=sxm,exm
      is=2*i-1
      c1=odlxx/rf(is-1,js,ks)
      c2=odlxx/rf(is+1,js,ks)
      c3=odlyy/rf(is,js-1,ks)
      c4=odlyy/rf(is,js+1,ks)
      c5=odlzz/rf(is,js,ks-1)
      c6=odlzz/rf(is,js,ks+1)
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

subroutine mgdppde_2d(sxm,exm,sym,eym,nxm,nym,cof,     &
                   sxf,exf,syf,eyf,rf,xl,yl,bd)
#include "mgd2.h"

integer sxm,exm,sym,eym,nxm,nym,sxf,exf,syf,eyf,bd(8)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6)
sll_real64 :: rf(sxf-1:exf+1,syf-1:eyf+1),xl,yl
!------------------------------------------------------------------------
! For the old version of the multigrid code, determine coefficients 
! for the pressure equation at all grid levels but the finest one.
! The coefficients are determined from the values of the density
! at integer nodes (i,j). Works only for periodic boundary conditions.
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
sll_real64 :: dlx,odlxx,dly,odlyy
sll_int32  :: i,j,is,js
sll_real64 :: c1, c2, c3, c4
# if cdebug
sll_real64 ::  tinitial
tinitial=MPI_WTIME()
# endif
!
! calculate off-diagonal terms
!
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

# if cdebug
timing(87)=timing(87)+MPI_WTIME()-tinitial
# endif

end subroutine

end module
