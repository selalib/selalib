      subroutine mgdphpde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,cof,
     1                    sx,ex,sy,ey,sz,ez,nxf,nyf,nzf,r,bd,
     2                    xl,yl,zl,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm
      integer sx,ex,sy,ey,sz,ez,nxf,nyf,nzf,bd(26),IOUT
      REALN cof(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1,8)
      REALN r(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1),xl,yl,zl
c------------------------------------------------------------------------
c For the new version of the multigrid code, determine the coefficients
c for the pressure equation at all grid levels except the finest one.
c The coefficients are determined directly from the density array r
c through some manipulation of indices and are values at (i+1/2,j+1/2)
c points. Works for periodic, Neumann, and Dirichlet boundary
c conditions. Should work even when there is no coarsifying in one
c direction.
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
      REALN dlx,eidlxx,dly,eidlyy,dlz,eidlzz
      integer i,j,k,im,jm,km,is,js,ks,istep,jstep,kstep
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c calculate off-diagonal terms
c
      dlx=xl/float(nxm-1)
      eidlxx=8.0d0/(dlx*dlx)
      dly=yl/float(nym-1)
      eidlyy=8.0d0/(dly*dly)
      dlz=zl/float(nzm-1)
      eidlzz=8.0d0/(dlz*dlz)
      istep=(nxf-1)/(nxm-1)
      jstep=(nyf-1)/(nym-1)
      kstep=(nzf-1)/(nzm-1)
      do k=szm,ezm
        km=2*kstep*k-3*(kstep-1)
        do j=sym,eym
          jm=2*jstep*j-3*(jstep-1)
          do i=sxm,exm
            im=2*istep*i-3*(istep-1)
            is=(im-istep)/2
            js=jm/2
            ks=km/2
            cof(i,j,k,1)=eidlxx/(r(is,js,ks)+r(is,js+1,ks)
     1                          +r(is,js,ks+1)+r(is,js+1,ks+1)
     2                          +r(is+1,js,ks)+r(is+1,js+1,ks)
     3                          +r(is+1,js,ks+1)+r(is+1,js+1,ks+1))
            is=(im+istep)/2
            cof(i,j,k,2)=eidlxx/(r(is,js,ks)+r(is,js+1,ks)
     1                          +r(is,js,ks+1)+r(is,js+1,ks+1)
     2                          +r(is+1,js,ks)+r(is+1,js+1,ks)
     3                          +r(is+1,js,ks+1)+r(is+1,js+1,ks+1))
            is=im/2
            js=(jm-jstep)/2
            cof(i,j,k,3)=eidlyy/(r(is,js,ks)+r(is,js+1,ks)
     1                          +r(is,js,ks+1)+r(is,js+1,ks+1)
     2                          +r(is+1,js,ks)+r(is+1,js+1,ks)
     3                          +r(is+1,js,ks+1)+r(is+1,js+1,ks+1))
            js=(jm-jstep)/2
            cof(i,j,k,4)=eidlyy/(r(is,js,ks)+r(is,js+1,ks)
     1                          +r(is,js,ks+1)+r(is,js+1,ks+1)
     2                          +r(is+1,js,ks)+r(is+1,js+1,ks)
     3                          +r(is+1,js,ks+1)+r(is+1,js+1,ks+1))
            js=jm/2
            ks=(km-kstep)/2
            cof(i,j,k,5)=eidlzz/(r(is,js,ks)+r(is,js+1,ks)
     1                          +r(is,js,ks+1)+r(is,js+1,ks+1)
     2                          +r(is+1,js,ks)+r(is+1,js+1,ks)
     3                          +r(is+1,js,ks+1)+r(is+1,js+1,ks+1))
            ks=(km+kstep)/2
            cof(i,j,k,6)=eidlzz/(r(is,js,ks)+r(is,js+1,ks)
     1                          +r(is,js,ks+1)+r(is,js+1,ks+1)
     2                          +r(is+1,js,ks)+r(is+1,js+1,ks)
     3                          +r(is+1,js,ks+1)+r(is+1,js+1,ks+1))
          end do
        end do
      end do
c
c enforce wall BCs
c
      if (bd(1).eq.1) then
        do k=szm,ezm
          do i=sxm,exm
            cof(i,eym,k,4)=0.0d0
          end do
        end do
      end if
      if (bd(3).eq.1) then
        do j=sym,eym
          do i=sxm,exm
            cof(i,j,szm,5)=0.0d0
          end do
        end do
      end if
      if (bd(5).eq.1) then
        do k=szm,ezm
          do i=sxm,exm
            cof(i,sym,k,3)=0.0d0
          end do
        end do
      end if
      if (bd(7).eq.1) then
        do j=sym,eym
          do i=sxm,exm
            cof(i,j,ezm,6)=0.0d0
          end do
        end do
      end if
      if (bd(9).eq.1) then
        do k=szm,ezm
          do j=sym,eym
            cof(sxm,j,k,1)=0.0d0
          end do
        end do
      end if
      if (bd(18).eq.1) then
        do k=szm,ezm
          do j=sym,eym
            cof(exm,j,k,2)=0.0d0
          end do
        end do
      end if
c
c calculate diagonal term
c
      do k=szm,ezm
        do j=sym,eym
          do i=sxm,exm
            cof(i,j,k,7)=-(cof(i,j,k,1)+cof(i,j,k,2)+cof(i,j,k,3)
     1                    +cof(i,j,k,4)+cof(i,j,k,5)+cof(i,j,k,6))
          end do
        end do
      end do
c
# if cdebug
      timing(83)=timing(83)+MPI_WTIME()-tinitial
# endif
      return
      end
