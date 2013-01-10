subroutine mgdpfpde(sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,cof,r,xl,yl,zl,IOUT)


use mpi
implicit none 
integer :: sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,IOUT
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

return
end
