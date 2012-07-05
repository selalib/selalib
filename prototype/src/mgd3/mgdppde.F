      subroutine mgdppde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,cof,
     1                   sxf,exf,syf,eyf,szf,ezf,rf,xl,yl,zl,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm 
      integer sxf,exf,syf,eyf,szf,ezf,IOUT
      REALN cof(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1,8)
      REALN rf(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1),xl,yl,zl
c------------------------------------------------------------------------
c Determine coefficients for the pressure equation at all grid levels
c but the finest one (this is done in mgdpfpde).
c
c cof array: 1 -> i-1
c            2 -> i+1
c            3 -> j-1
c            4 -> j+1
c            5 -> k-1
c            6 -> k+1
c            7 -> central
c
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      REALN dlx,odlxx,dly,odlyy,dlz,odlzz,c1,c2,c3,c4,c5,c6
      integer i,j,k,is,js,ks
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
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
c
# if cdebug
      timing(87)=timing(87)+MPI_WTIME()-tinitial
# endif
      return
      end
