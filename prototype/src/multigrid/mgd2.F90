!> Multigrid solver in 2 dimensions on a cartesian grid
!> \authors Pierre Navaro
!> \remark this a refactored version of Bernard Bunner code 
!> (http://cs-www.cs.yale.edu/homes/douglas-craig/mgnet-codes-bunner.html)
module mgd2
#include "sll_working_precision.h"
use mpi
use gxch1_2d

implicit none

sll_int32, private :: sx, ex, sy, ey, i, j, k

logical ::  WMGD 
sll_int32 :: nxk(20),nyk(20),sxk(20),exk(20),syk(20),eyk(20)
sll_int32 :: kpbgn(20),kcbgn(20),ikdatatype(20),jkdatatype(20)
sll_int32 :: ijkdatatype(20),sxi(20),exi(20),syi(20),eyi(20)
sll_int32 :: nxr(20),nyr(20),sxr(20),exr(20),syr(20),eyr(20)
sll_int32 :: irdatatype(20),jrdatatype(20),ijrdatatype(20)

type, public :: block
   sll_int32  :: id
   sll_int32  :: sx, ex
   sll_int32  :: sy, ey
   sll_int32  :: ixp, jyq
   sll_int32  :: iex, jey
   sll_int32  :: ngrid
   sll_int32  :: neighbor(8),bd(8)
end type block

contains

!> Define the 3 derived datatypes needed to communicate the boundary
!> data of (sx-1:ex+1,sy-1:ey+1) arrays between 'myid' and its 8
!> neighbors
subroutine grid1_type(itype,jtype,ijtype,sx,ex,sy,ey)
sll_int32  :: itype,jtype,ijtype,realtype,sx,ex,sy,ey
sll_int32  :: ierr

realtype = MPI_REAL8
!
! datatype for one row
!
call MPI_TYPE_CONTIGUOUS(ex-sx+1,realtype,itype,ierr)
call MPI_TYPE_COMMIT(itype,ierr)
!
! datatype for one column
!
call MPI_TYPE_VECTOR(ey-sy+1,1,ex-sx+3,realtype,jtype,ierr)
call MPI_TYPE_COMMIT(jtype,ierr)
!
! datatype for one 1*1 corner
!
call MPI_TYPE_CONTIGUOUS(1,realtype,ijtype,ierr)
call MPI_TYPE_COMMIT(ijtype,ierr)

return
end subroutine


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
subroutine mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,cof,r,xl,yl,bd)

sll_int32  :: sxf,exf,syf,eyf,nxf,nyf,bd(8)
sll_real64 :: cof(sxf-1:exf+1,syf-1:eyf+1,6)
sll_real64 :: r(sxf-1:exf+1,syf-1:eyf+1),xl,yl
sll_real64 :: dlx,todlxx,dly,todlyy,rij
sll_int32  :: i,j

dlx=xl/float(nxf-1)
dly=yl/float(nyf-1)

   todlxx=2.0d0/(dlx*dlx)
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

end subroutine


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
subroutine mgdphpde(sxm,exm,sym,eym,nxm,nym,cof,sx,ex,sy,ey,nxf,nyf,bd,xl,yl)


sll_int32  :: sxm,exm,sym,eym,nxm,nym,sx,ex,sy,ey,nxf,nyf,bd(8)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6)
sll_real64 :: xl,yl
sll_real64 :: dlx,fodlxx,dly,fodlyy
sll_int32  :: i,j,im,jm,is,js,istep,jstep

! calculate off-diagonal terms

dlx=xl/float(nxm-1)
dly=yl/float(nym-1)

istep=(nxf-1)/(nxm-1)
jstep=(nyf-1)/(nym-1)

fodlxx=1.0d0/(dlx*dlx)
fodlyy=1.0d0/(dly*dly)

do j=sym,eym
  jm=2*jstep*j-3*(jstep-1)
  do i=sxm,exm
    im=2*istep*i-3*(istep-1)
    is=(im-istep)/2
    js=jm/2
    cof(i,j,1)=fodlxx
    is=(im+istep)/2
    cof(i,j,2)=fodlxx
    is=im/2
    js=(jm-jstep)/2
    cof(i,j,3)=fodlyy
    js=(jm+jstep)/2
    cof(i,j,4)=fodlyy
  end do
end do

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
subroutine mgdppde(sxm,exm,sym,eym,nxm,nym,cof,     &
                   sxf,exf,syf,eyf,rf,xl,yl,bd)

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

!> Discretize the pde: set the coefficients of the cof matrix. Works
!> for periodic, Neumann, and Dirichlet boundary conditions.
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
subroutine mgdrpde(sxm,exm,sym,eym,nxm,nym,cof,xl,yl,bd)
sll_int32  :: sxm,exm,sym,eym,nxm,nym,bd(8)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6),xl,yl
sll_real64 :: dlx,odlxx,dly,odlyy
sll_int32  :: i,j

! calculate off-diagonal terms
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

!
! enforce Neumann BCs
!
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
!
! calculate diagonal term
!
do j=sym,eym
  do i=sxm,exm
    cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
  end do
end do

end subroutine

!> For the old version of the multigrid code, set the fine grid values 
!> of the density in the work vector
subroutine mgdrsetf(sxf,exf,syf,eyf,rf,r)

sll_int32  :: sxf,exf,syf,eyf
sll_real64 :: rf(sxf-1:exf+1,syf-1:eyf+1),r(sxf-1:exf+1,syf-1:eyf+1)
sll_int32  :: i,j

do j=syf-1,eyf+1
  do i=sxf-1,exf+1
    rf(i,j)=r(i,j)
  end do
end do

end subroutine

!> For the old version of the multigrid code, transfer values of the 
!> density from a finer to a coarser grid level. It is necessary to
!> exchange the boundary density data because the grid "shifts" to
!> the right as it becomes coarser. (In the new version of the
!> multigrid code, there is no such shift, hence no communication is 
!> needed).
subroutine mgdrtrsf(sxc,exc,syc,eyc,nxc,nyc,rc,            &
                       sxf,exf,syf,eyf,nxf,nyf,rf,            &
                       comm2d,myid,neighbor,bd,itype,jtype)
sll_int32  :: sxc,exc,syc,eyc,nxc,nyc,sxf,exf,syf,eyf,nxf,nyf
sll_int32  :: comm2d,myid,neighbor(8),bd(8),itype,jtype

sll_real64 :: rc(sxc-1:exc+1,syc-1:eyc+1)
sll_real64 :: rf(sxf-1:exf+1,syf-1:eyf+1)

sll_int32  :: i,j,ic,jc,i1,i2,j1,j2

if (nxc.lt.nxf) then
  i1=1
  i2=0
else 
  i1=0
  i2=1
end if
if (nyc.lt.nyf) then
  j1=1
  j2=0
else
  j1=0
  j2=1
end if
do jc=syc,eyc
  j=j1*(2*jc-1)+j2*jc
  do ic=sxc,exc
    i=i1*(2*ic-1)+i2*ic
    rc(ic,jc)=rf(i,j)
  end do
end do

! exchange the boundary values (need only lines, not corner)
call gxch1lin(rc,comm2d,sxc,exc,syc,eyc,neighbor,bd,itype,jtype)

end subroutine

!> Set the fine grid values in the work vector
subroutine mgdsetf(sxf,exf,syf,eyf,phi,rhs,phif,rhsf)
sll_int32  :: sxf,exf,syf,eyf
sll_real64 :: phi(sxf-1:exf+1,syf-1:eyf+1),rhs(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: phif(sxf-1:exf+1,syf-1:eyf+1),rhsf(sxf-1:exf+1,syf-1:eyf+1)
sll_int32  :: i,j

do j=syf-1,eyf+1
  do i=sxf-1,exf+1
    phi(i,j)=phif(i,j)
    rhs(i,j)=rhsf(i,j)
  end do
end do

end subroutine

!> Free the MPI datatypes associated witht the multigrid code
subroutine mgdend(ngrid)

sll_int32 :: k,ierr
sll_int32 :: ngrid

do k=1,ngrid-1
   call MPI_TYPE_FREE(ikdatatype(k),ierr)
   call MPI_TYPE_FREE(jkdatatype(k),ierr)
   call MPI_TYPE_FREE(ijkdatatype(k),ierr)
 if (.not. WMGD) then
   call MPI_TYPE_FREE(irdatatype(k),ierr)
   call MPI_TYPE_FREE(jrdatatype(k),ierr)
   call MPI_TYPE_FREE(ijrdatatype(k),ierr)
 endif
end do

!call MPI_FINALIZE(ierr)


end subroutine

!> Add correction from coarse grid level to fine grid level. Uses
!> bilinear interpolation for the old version of the multigrid code,
!> and area weighting for its new version.
!>
!> Tested for the case where coarsifying takes place in all directions
subroutine mgdcor(sxf,exf,syf,eyf,nxf,nyf,phif,  &
                  sxc,exc,syc,eyc,nxc,nyc,phic,  &
                  sx1,ex1,sy1,ey1,bd,phibc)

sll_int32  :: sxf,exf,syf,eyf,nxf,nyf
sll_int32  :: sxc,exc,syc,eyc,nxc,nyc,sx1,ex1,sy1,ey1,bd(8)
sll_real64 :: phif(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: phic(sxc-1:exc+1,syc-1:eyc+1),phibc(4)

sll_int32 :: i,j,ic,jc,i1,i2,j1,j2

if (WMGD) then
!------------------------------------------------------------------------
! new version: the correction is the weighted average of either two 
! or four points at the coarser grid level depending on whether 
! coarsifying takes place in all directions or not
!
if (nxf.eq.nxc) then
  do jc=syc-1,eyc
    j=2*jc-1
    do i=sxf-1,exf+1
      ic=i
      phif(i,j)=phif(i,j)+(3.0d0*phic(ic,jc)+phic(ic,jc+1))/4.0d0
      phif(i,j+1)=phif(i,j+1)+(phic(ic,jc)+3.0d0*phic(ic,jc+1))/4.0d0
    end do
  end do
else if (nyf.eq.nyc) then
  do j=syf-1,eyf+1
    jc=j
    do ic=sxc-1,exc
      i=2*ic-1
      phif(i,j)=phif(i,j)+(3.0d0*phic(ic,jc)+phic(ic+1,jc))/4.0d0
      phif(i+1,j)=phif(i+1,j)+(phic(ic,jc)+3.0d0*phic(ic+1,jc))/4.0d0
    end do
  end do
else
  do jc=syc-1,eyc
    j=2*jc-1
    do ic=sxc-1,exc
      i=2*ic-1
      phif(i,j)=phif(i,j)+  &
        (9.0d0*phic(ic,jc)+3.0d0*phic(ic+1,jc)  &
        +3.0d0*phic(ic,jc+1)+phic(ic+1,jc+1))/16.0d0
      phif(i+1,j)=phif(i+1,j)+  &
        (3.0d0*phic(ic,jc)+9.0d0*phic(ic+1,jc)  &
        +phic(ic,jc+1)+3.0d0*phic(ic+1,jc+1))/16.0d0
      phif(i,j+1)=phif(i,j+1)+  &
        (3.0d0*phic(ic,jc)+phic(ic+1,jc)  &
        +9.0d0*phic(ic,jc+1)+3.0d0*phic(ic+1,jc+1))/16.0d0
      phif(i+1,j+1)=phif(i+1,j+1)+  &
        (phic(ic,jc)+3.0d0*phic(ic+1,jc)  &
        +3.0d0*phic(ic,jc+1)+9.0d0*phic(ic+1,jc+1))/16.0d0
    end do
  end do
end if
!
! impose Neumann and Dirichlet boundary conditions
! TEMP: periodicity is not enforced to save one call to gxch1lin;
! check whether it has an impact or not...
!
      call mgdbdry(sxf,exf,syf,eyf,phif,bd,phibc)
else  
!------------------------------------------------------------------------
! old version
!
if (nxc.lt.nxf) then
  i1=1
  i2=0
else
  i1=0
  i2=1
end if
if (nyc.lt.nyf) then
  j1=1
  j2=0
else
  j1=0
  j2=1
end if
!
! identity at the points of the fine grid which have odd indices 
! in i and j
!
do jc=sy1,ey1
  j=j1*(2*jc-1)+j2*jc
  do ic=sx1,ex1
    i=i1*(2*ic-1)+i2*ic
    phif(i,j)=phif(i,j)+phic(ic,jc)
  end do
end do
!
! interpolation of the two neighboring values for the points of the 
! fine grid with even index for i and odd index for j
!
if (nxc.lt.nxf) then
  do jc=sy1,ey1
    j=j1*(2*jc-1)+j2*jc
    do ic=sxc-1,exc
      i=2*ic
      phif(i,j)=phif(i,j)+0.5d0*(phic(ic,jc)+phic(ic+1,jc))
    end do
  end do
end if
!
! interpolation of the two neighboring values for the points of the 
! fine grid with odd index for i and even index for j
!
if (nyc.lt.nyf) then
  do jc=syc-1,eyc
    j=2*jc
    do ic=sx1,ex1
      i=i1*(2*ic-1)+i2*ic
      phif(i,j)=phif(i,j)+0.5d0*(phic(ic,jc)+phic(ic,jc+1))
    end do
  end do
end if
!
! interpolation of the four neighboring values for the points of the
! fine grid with even indices for i and j
!
if (nxc.lt.nxf.and.nyc.lt.nyf) then
  do jc=syc-1,eyc
    j=2*jc
    do ic=sxc-1,exc
      i=2*ic
      phif(i,j)=phif(i,j)+0.25d0*(phic(ic,jc)+phic(ic+1,jc) &
                                 +phic(ic,jc+1)+phic(ic+1,jc+1))
    end do
  end do
end if
endif

end subroutine


!> Gauss-Seidel point relaxation with Red & Black ordering. Works for
!> periodic, Neumann, and Dirichlet boundary conditions.
subroutine mgdrelax(sxm,exm,sym,eym,phi,cof,iters,comm2d,myid, &
                    neighbor,bd,phibc,itype,jtype)


sll_int32 :: sxm,exm,sym,eym,iters
sll_int32 :: comm2d,myid,neighbor(8),bd(8),itype,jtype
sll_real64 :: phi(sxm-1:exm+1,sym-1:eym+1)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6),phibc(4)
sll_int32 :: rb,it,ipass,i,j
!
! do iters sweeps in the subdomain
!
do it=1,iters
  rb=mod(sxm,2)
  do ipass=1,2
    do j=sym,eym
      do i=sxm+rb,exm,2
        phi(i,j)=(cof(i,j,6)-(cof(i,j,1)*phi(i-1,j)               &
                             +cof(i,j,2)*phi(i+1,j)               &
                             +cof(i,j,3)*phi(i,j-1)               &
                             +cof(i,j,4)*phi(i,j+1)))/cof(i,j,5)
      end do
      rb=1-rb
    end do
    rb=1-mod(sxm,2)
if (WMGD) then
!
! new version: impose Neumann and Dirichlet boundary conditions
!
    call mgdbdry(sxm,exm,sym,eym,phi,bd,phibc)
endif
  end do
end do
!
! Exchange boundary data only once at the end. Since the number
! of relaxation sweeps at each level is characteristically small
! (1 or 2 are common values), this does not damage the convergence
! rate too badly. Overall, I have found a significant reduction
! in execution time. This also imposes the periodic BCs.
!
call gxch1lin(phi,comm2d,sxm,exm,sym,eym,neighbor,bd,itype,jtype)
if (WMGD) then
!
! new version: impose Neumann and Dirichlet boundary conditions
!
call mgdbdry(sxm,exm,sym,eym,phi,bd,phibc)
endif 

end subroutine

!> Calculate the residual and restrict it to the coarser level. In
!> the new version, the restriction involves 4 points. In the old
!> version, it involves 9 points (5 if half-weighting).
subroutine mgdrestr(sxc,exc,syc,eyc,nxc,nyc,phic,rhsc,   &
                    sxf,exf,syf,eyf,nxf,nyf,phif,cof,    &
                    resf,iresw,comm2d,myid,neighbor,bd,  &
                    itype,jtype,ijtype)

sll_int32  :: sxc,exc,syc,eyc,nxc,nyc,iresw
sll_int32  :: sxf,exf,syf,eyf,nxf,nyf
sll_int32  :: comm2d,myid,neighbor(8),bd(8),itype,jtype,ijtype
sll_real64 :: phic(sxc-1:exc+1,syc-1:eyc+1)
sll_real64 :: rhsc(sxc-1:exc+1,syc-1:eyc+1)
sll_real64 :: phif(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: resf(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: cof(sxf-1:exf+1,syf-1:eyf+1,6)
sll_int32  :: i,j,isrt,jsrt,iinc,jinc,ic,jc
!------------------------------------------------------------------------
do jc=syc-1,eyc+1
  do ic=sxc-1,exc+1
    phic(ic,jc)=0.0d0
    rhsc(ic,jc)=0.0d0
  end do
end do
!
! calculate residual
!
do j=syf,eyf
  do i=sxf,exf
    resf(i,j)=cof(i,j,6)-(cof(i,j,1)*phif(i-1,j) &
                         +cof(i,j,2)*phif(i+1,j) &
                         +cof(i,j,3)*phif(i,j-1) &
                         +cof(i,j,4)*phif(i,j+1) &
                         +cof(i,j,5)*phif(i,j))
  end do
end do
if (WMGD) then
!------------------------------------------------------------------------
! new version: calculate the right-hand side at the coarser grid
! level from the averages of the values at the 4 surrounding points;
! if there is no coarsifying in one direction, only 2 points are
! used; no exchange of boundary data is necessary
!
if (nxc.eq.nxf) then
  do jc=syc,eyc
    j=2*jc-2
    do ic=sxc,exc
      i=ic
      rhsc(ic,jc)=0.5d0*(resf(i,j)+resf(i,j+1))
    end do
  end do
else if (nyc.eq.nyf) then
  do jc=syc,eyc
    j=jc
    do ic=sxc,exc
      i=2*ic-2
      rhsc(ic,jc)=0.5d0*(resf(i,j)+resf(i+1,j))
    end do
  end do
else
  do jc=syc,eyc
    j=2*jc-2
    do ic=sxc,exc
      i=2*ic-2
      rhsc(ic,jc)=0.25d0*(resf(i,j)+resf(i+1,j) &
                   +resf(i,j+1)+resf(i+1,j+1))
    end do
  end do
end if
else  
!------------------------------------------------------------------------
! old version: have to exchange boundary data; if full-weighting, 
! need to exchange also corner points
!
call gxch1lin(resf,comm2d,sxf,exf,syf,eyf,neighbor,bd,itype,jtype)
if (iresw.eq.1) then
  call gxch1cor(resf,comm2d,sxf,exf,syf,eyf,neighbor,bd,ijtype)
end if
!
! restrict it to coarser level
!
if (nxc.lt.nxf) then
  isrt=2*sxc-1
  iinc=2
else
  isrt=sxc
  iinc=1
end if
if (nyc.lt.nyf) then
  jsrt=2*syc-1
  jinc=2
else
  jsrt=syc
  jinc=1
end if

if (iresw.eq.1) then
!
! use full weighting
!
  j=jsrt
  do jc=syc,eyc
    i=isrt
    do ic=sxc,exc
      rhsc(ic,jc)=0.25d0*resf(i,j)                        &
                 +0.125d0*(resf(i+1,j)+resf(i-1,j)        &
                          +resf(i,j+1)+resf(i,j-1))       &
                 +0.0625d0*(resf(i+1,j+1)+resf(i+1,j-1)   &
                           +resf(i-1,j-1)+resf(i-1,j+1))
      i=i+iinc
    end do
    j=j+jinc
  end do
else if (iresw.eq.2) then
!
! use half-weighting
!
  j=jsrt
  do jc=syc,eyc
    i=isrt
    do ic=sxc,exc
      rhsc(ic,jc)=0.5d0*resf(i,j)                     &
                 +0.125d0*(resf(i+1,j)+resf(i-1,j)    &
                          +resf(i,j+1)+resf(i,j-1))
      i=i+iinc
    end do
    j=j+jinc
  end do
end if
endif

end subroutine


!> Calculate the error between the new and old iterates of phi and 
!> save the new iterate into the phio array.
subroutine mgderr(relmax,sxm,exm,sym,eym,phio,phin,comm2d)

sll_int32  :: sxm,exm,sym,eym,comm2d
sll_real64 :: relmax
sll_real64 :: phio(sxm-1:exm+1,sym-1:eym+1),phin(sxm-1:exm+1,sym-1:eym+1)
sll_real64 :: phloc,reloc
sll_int32  :: i,j,ierr
!
! calculate local error
!
phloc=0.0d0
reloc=0.0d0
do j=sym,eym
  do i=sxm,exm
    phloc=max(phloc,abs(phin(i,j)))
    reloc=max(reloc,abs(phin(i,j)-phio(i,j)))
  end do
end do
if (phloc.gt.0.0) then
  reloc=reloc/phloc
else
  reloc=0.0d0
end if

! global reduce across all processes
call MPI_ALLREDUCE(reloc,relmax,1,MPI_REAL8,MPI_MAX,comm2d,ierr)

! save new values into ouput array
do j=sym-1,eym+1
  do i=sxm-1,exm+1
    phio(i,j)=phin(i,j)
  end do
end do

end subroutine

!> Enforce the Neumann and Dirichlet boundary conditions
subroutine mgdbdry(sxm,exm,sym,eym,phi,bd,phibc)


sll_int32  :: sxm,exm,sym,eym,bd(8)
sll_real64 :: phi(sxm-1:exm+1,sym-1:eym+1),phibc(4)
sll_int32 :: i,j

if (bd(1).eq.1) then
   do j=sym-1,eym+1
      phi(exm+1,j)=phi(exm,j)
   end do
else if (bd(1).eq.2) then
   do j=sym-1,eym+1
      phi(exm+1,j)=2.0d0*phibc(1)-phi(exm,j)
   end do
end if
if (bd(5).eq.1) then
   do j=sym-1,eym+1
      phi(sxm-1,j)=phi(sxm,j)
   end do
else if (bd(5).eq.2) then
   do j=sym-1,eym+1
      phi(sxm-1,j)=2.0d0*phibc(3)-phi(sxm,j)
   end do
end if
if (bd(3).eq.1) then
   do i=sxm-1,exm+1
      phi(i,sym-1)=phi(i,sym)
   end do
else if (bd(3).eq.2) then
   do i=sxm-1,exm+1
      phi(i,sym-1)=2.0d0*phibc(2)-phi(i,sym)
   end do
end if
if (bd(7).eq.1) then
   do i=sxm-1,exm+1
      phi(i,eym+1)=phi(i,eym)
   end do
else if (bd(7).eq.2) then
   do i=sxm-1,exm+1
      phi(i,eym+1)=2.0d0*phibc(4)-phi(i,eym)
   end do
end if

end subroutine

!> Do one multigrid K-cycle
!> K=1 -> V-cycle
!> K=2 -> W-cycle
subroutine mgdkcyc(work,res,kcur,kcycle,iprer,ipost,iresw, &
                   comm2d,myid,neighbor,bd,phibc)

sll_int32  :: kcur,kcycle,iprer,ipost,iresw
sll_int32  :: comm2d,myid,neighbor(8),bd(8)
sll_real64 :: work(*),res(*),phibc(4,*)

sll_int32 :: sxf,exf,syf,eyf,nxf,nyf,ipf,icf
sll_int32 :: sxc,exc,syc,eyc,nxc,nyc,ipc,irc
sll_int32 :: klevel,itype,jtype,ijtype,kount(20),l,nrel
sll_int32 :: sx1,ex1,sy1,ey1

klevel=kcur
!
! pre-relax at current finest grid level
!
sxf=sxk(klevel)
exf=exk(klevel)
syf=syk(klevel)
eyf=eyk(klevel)
nxf=nxk(klevel)
nyf=nyk(klevel)
ipf=kpbgn(klevel)
icf=kcbgn(klevel)
itype=ikdatatype(klevel)
jtype=jkdatatype(klevel)
call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),iprer, &
              comm2d,myid,neighbor,bd,phibc(1,klevel),   &
              itype,jtype)
!
! if at coarsest level, post-relax
!
if (kcur.eq.1) goto 5
!
! calculate residual and restrict it to kcur-1 
! (note: the residuals are stored in the memory space used by the
! the rhs array rhsf, which is therefore erased)
!
sxc=sxk(klevel-1)
exc=exk(klevel-1)
syc=syk(klevel-1)
eyc=eyk(klevel-1)
nxc=nxk(klevel-1)
nyc=nyk(klevel-1)
ipc=kpbgn(klevel-1)
irc=kcbgn(klevel-1)+5*(exc-sxc+3)*(eyc-syc+3)
itype=ikdatatype(klevel)
jtype=jkdatatype(klevel)
ijtype=ijkdatatype(klevel)
call mgdrestr(sxc,exc,syc,eyc,nxc,nyc,work(ipc),work(irc), &
              sxf,exf,syf,eyf,nxf,nyf,work(ipf),work(icf), &
              res,iresw,comm2d,myid,neighbor,bd,           &
              itype,jtype,ijtype)
!
! set counter for grid levels to zero
!
do l=1,kcur
  kount(l)=0
end do
!
! set new level and continue K-cycling
!
klevel=kcur-1
nrel=iprer
!
!------------------------------------------------------------------------
! K-cycle control point
!
10    continue
!
! post-relax when kcur revisited
!
if (klevel.eq.kcur) goto 5
!
! count "hit" at current level
!
kount(klevel)=kount(klevel)+1
!
! relax at current level
!
sxf=sxk(klevel)
exf=exk(klevel)
syf=syk(klevel)
eyf=eyk(klevel)
nxf=nxk(klevel)
nyf=nyk(klevel)
ipf=kpbgn(klevel)
icf=kcbgn(klevel)
itype=ikdatatype(klevel)
jtype=jkdatatype(klevel)
call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),nrel, &
              comm2d,myid,neighbor,bd,phibc(1,klevel),itype,jtype)

if (kount(klevel).eq.(kcycle+1)) then
!
! K-cycle(iprer,ipost) complete at klevel
! inject correction to finer grid
!      
  sxf=sxk(klevel+1)
  exf=exk(klevel+1)
  syf=syk(klevel+1)
  eyf=eyk(klevel+1)
  nxf=nxk(klevel+1)
  nyf=nyk(klevel+1)
  ipf=kpbgn(klevel+1)
  sxc=sxk(klevel)
  exc=exk(klevel)
  syc=syk(klevel)
  eyc=eyk(klevel)
  nxc=nxk(klevel)
  nyc=nyk(klevel)
  ipc=kpbgn(klevel)
  sx1=sxi(klevel)
  ex1=exi(klevel)
  sy1=syi(klevel)
  ey1=eyi(klevel)
  call mgdcor(sxf,exf,syf,eyf,nxf,nyf,work(ipf),         &
              sxc,exc,syc,eyc,nxc,nyc,work(ipc),         &
              sx1,ex1,sy1,ey1,bd,phibc(1,klevel+1))
!
! reset counter to zero at klevel
!
  kount(klevel)=0
!
! ascend to next higher level and set to post-relax there
!
  klevel=klevel+1
  nrel=ipost
  goto 10
  
else

! K-cycle not complete so descend unless at coarsest

  if (klevel.gt.1) then
    sxc=sxk(klevel-1)
    exc=exk(klevel-1)
    syc=syk(klevel-1)
    eyc=eyk(klevel-1)
    nxc=nxk(klevel-1)
    nyc=nyk(klevel-1)
    ipc=kpbgn(klevel-1)
    irc=kcbgn(klevel-1)+5*(exc-sxc+3)*(eyc-syc+3)
    itype=ikdatatype(klevel)
    jtype=jkdatatype(klevel)
    ijtype=ijkdatatype(klevel)
    call mgdrestr(sxc,exc,syc,eyc,nxc,nyc,work(ipc),work(irc),  &
                  sxf,exf,syf,eyf,nxf,nyf,work(ipf),work(icf),  &
                  res,iresw,comm2d,myid,neighbor,bd,            &
                  itype,jtype,ijtype)
!
! pre-relax at next coarser level
!
    klevel=klevel-1
    nrel=iprer
    goto 10

  else
!
! post-relax at coarsest level
!
    sxf=sxk(klevel)
    exf=exk(klevel)
    syf=syk(klevel)
    eyf=eyk(klevel)
    ipf=kpbgn(klevel)
    icf=kcbgn(klevel)
    itype=ikdatatype(klevel)
    jtype=jkdatatype(klevel)
    call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),ipost, &
                  comm2d,myid,neighbor,bd,phibc(1,klevel),   &
                  itype,jtype)
!
! inject correction to grid level 2
!
    sxf=sxk(2)
    exf=exk(2)
    syf=syk(2)
    eyf=eyk(2)
    nxf=nxk(2)
    nyf=nyk(2)
    ipf=kpbgn(2)
    sxc=sxk(1)
    exc=exk(1)
    syc=syk(1)
    eyc=eyk(1)
    nxc=nxk(1)
    nyc=nyk(1)
    ipc=kpbgn(1)
    sx1=sxi(1)
    ex1=exi(1)
    sy1=syi(1)
    ey1=eyi(1)
    call mgdcor(sxf,exf,syf,eyf,nxf,nyf,work(ipf),  &
                sxc,exc,syc,eyc,nxc,nyc,work(ipc),  &
                sx1,ex1,sy1,ey1,bd,phibc(1,2))
!
! set to post-relax at level 2
!
    nrel=ipost
    klevel=2
    goto 10

  end if

end if

5     continue
!
!------------------------------------------------------------------------
! post-relax at kcur level
!
sxf=sxk(kcur)
exf=exk(kcur)
syf=syk(kcur)
eyf=eyk(kcur)
ipf=kpbgn(kcur)
icf=kcbgn(kcur)
itype=ikdatatype(kcur)
jtype=jkdatatype(kcur)

call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),ipost, &
              comm2d,myid,neighbor,bd,phibc(1,kcur),     &
              itype,jtype)

end subroutine


!> Rescale the field a so that its average inside the domain
!> remains constant and equal to avo. For the density,avo should
!> be rro, this ensures conservation of mass. For the pressure,
!> avo should be 0 so that the average pressure does not drift
!> away from 0, which is the initial value.
subroutine gscale(sx,ex,sy,ey,a,avo,acorr,comm2d,nx,ny)

sll_int32  :: sx,ex,sy,ey,nx,ny
sll_real64 :: a(sx-1:ex+1,sy-1:ey+1),avo,acorr
sll_int32  :: comm2d
sll_real64 :: avloc,av
sll_int32  :: i,j,ierr

! determine average value
avloc=0.0d0
do j=sy,ey
   do i=sx,ex
      avloc=avloc+a(i,j)
   end do
end do
!
! global reduce across all process
!
call MPI_ALLREDUCE(avloc,av,1,MPI_REAL8,MPI_SUM,comm2d,ierr)

av=av/float(nx*ny)
!
! do correction
!
acorr=avo-av
do j=sy,ey
  do i=sx,ex
    a(i,j)=a(i,j)+acorr
  end do
end do

end subroutine 


!> Determine coefficients at the finest grid
!> level. These coefficients involve radius half-way between the
!> nodes. Works for periodic, Neumann, and Dirichlet boundary conditions.
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
subroutine polar_pfpde(sxf,exf,syf,eyf,nxf,nyf,cof,r,xl,yl,bd)

sll_int32  :: sxf,exf,syf,eyf,nxf,nyf,bd(8)
sll_real64 :: cof(sxf-1:exf+1,syf-1:eyf+1,6)
sll_real64 :: r(sxf-1:exf+1),xl,yl
sll_real64 :: dlx,todlxx,dly,todlyy,rij
sll_int32  :: i,j

dlx=xl/float(nxf-1)
dly=yl/float(nyf-1)

todlxx=1.0d0/(dlx*dlx)
todlyy=1.0d0/(dly*dly)
do j=syf,eyf
   do i=sxf,exf
     rij=r(i)
     cof(i,j,1)=todlxx + 1./(dlx*(rij+r(i-1)))
     cof(i,j,2)=todlxx - 1./(dlx*(rij+r(i+1)))
     cof(i,j,3)= 1./(rij*rij*dly)
     cof(i,j,4)= 1./(rij*rij*dly)
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

end subroutine


!> Determine the coefficients for the equation at all grid levels 
!> except the finest one.
!> The coefficients are determined directly from the array r
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
subroutine polar_phpde(sxm,exm,sym,eym,nxm,nym,cof,         &
                       sx,ex,sy,ey,nxf,nyf,r,bd,xl,yl)

sll_int32  :: sxm,exm,sym,eym,nxm,nym,sx,ex,sy,ey,nxf,nyf,bd(8)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6)
sll_real64 :: r(sx-1:ex+1),xl,yl
sll_real64 :: dlx,fodlxx,dly,fodlyy
sll_int32  :: i,j,im,jm,is,js,istep,jstep

! calculate off-diagonal terms

dlx=xl/float(nxm-1)
dly=yl/float(nym-1)

istep=(nxf-1)/(nxm-1)
jstep=(nyf-1)/(nym-1)

fodlxx=1.0d0/(dlx*dlx)
fodlyy=1.0d0/(dly*dly)

do j=sym,eym
  jm=2*jstep*j-3*(jstep-1)
  do i=sxm,exm
    im=2*istep*i-3*(istep-1)
    is=(im-istep)/2
    js=jm/2
    cof(i,j,1)=fodlxx+1./((r(is)+r(is+1))*dlx)
    is=(im+istep)/2
    cof(i,j,2)=fodlxx-1./((r(is)+r(is+1))*dlx)
    is=im/2
    js=(jm-jstep)/2
    cof(i,j,3)=fodlyy/(r(is)*r(is))
    js=(jm+jstep)/2
    cof(i,j,4)=fodlyy/(r(is)*r(is))
  end do
end do

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
subroutine polar_mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,cof,r,xl,yl,bd)

sll_int32  :: sxf,exf,syf,eyf,nxf,nyf,bd(8)
sll_real64 :: cof(sxf-1:exf+1,syf-1:eyf+1,6)
sll_real64 :: r(sxf-1:exf+1,syf-1:eyf+1),xl,yl
sll_real64 :: dlx,todlxx,dly,todlyy,rij
sll_int32  :: i,j

dlx=xl/float(nxf-1)
dly=yl/float(nyf-1)

todlxx=1.0d0/(dlx*dlx)
todlyy=1.0d0/(dly*dly)
do j=syf,eyf
  do i=sxf,exf
    rij=r(i,j)
    cof(i,j,1)=1./(dlx*(rij+r(i-1,j)))+todlxx
    cof(i,j,2)=1./(dlx*(rij+r(i+1,j)))+todlxx
    cof(i,j,3)=4./((rij+r(i,j-1))*(rij+r(i,j-1))*dly)
    cof(i,j,4)=4./((rij+r(i,j+1))*(rij+r(i,j-1))*dly)
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

end subroutine

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
subroutine polar_mgdphpde(sxm,exm,sym,eym,nxm,nym,cof,         &
                          sx,ex,sy,ey,nxf,nyf,r,bd,xl,yl)

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


end module mgd2
