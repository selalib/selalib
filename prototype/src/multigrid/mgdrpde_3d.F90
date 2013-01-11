subroutine mgdrpde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,cof,xl,yl,zl,IOUT)

use mpi
implicit none 
integer sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,IOUT
real(8) :: cof(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1,8),xl,yl,zl
!------------------------------------------------------------------------
! Discretize the pde: set the coefficients of the cof matrix
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
integer i,j,k

dlx=xl/float(nxm-1)
odlxx=1.0d0/(dlx*dlx)
dly=yl/float(nym-1)
odlyy=1.0d0/(dly*dly)
dlz=zl/float(nzm-1)
odlzz=1.0d0/(dlz*dlz)
do k=szm,ezm
  do j=sym,eym
    do i=sxm,exm
      c1=odlxx
      c2=odlxx
      c3=odlyy
      c4=odlyy
      c5=odlzz
      c6=odlzz
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
