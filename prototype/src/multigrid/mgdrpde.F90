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
      subroutine mgdrpde(sxm,exm,sym,eym,nxm,nym,cof,xl,yl,bd,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,nxm,nym,bd(8),IOUT
      REALN cof(sxm-1:exm+1,sym-1:eym+1,6),xl,yl
c------------------------------------------------------------------------
c Discretize the pde: set the coefficients of the cof matrix. Works
c for periodic, Neumann, and Dirichlet boundary conditions.
c
c cof array:
c
c         cof(4)
c           |
c           |
c cof(1)--cof(5)--cof(2)
c           |
c           |
c         cof(3)
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      REALN dlx,odlxx,dly,odlyy
      integer i,j
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c calculate off-diagonal terms
c
      dlx=xl/float(nxm-1)
      odlxx=1.0d0/(dlx*dlx)
      dly=yl/float(nym-1)
      odlyy=1.0d0/(dly*dly)
      do j=sym,eym
        do i=sxm,exm
          cof(i,j,1)=odlxx
          cof(i,j,2)=odlxx
          cof(i,j,3)=odlyy
          cof(i,j,4)=odlyy
        end do
      end do
c
c enforce Neumann BCs
c
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
c
c calculate diagonal term
c
      do j=sym,eym
        do i=sxm,exm
          cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
        end do
      end do
c
# if cdebug
      timing(82)=timing(82)+MPI_WTIME()-tinitial
# endif
      return
      end
