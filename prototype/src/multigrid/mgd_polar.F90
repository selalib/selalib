!> Multigrid solver in 2 dimensions on a cartesian grid
!> \authors Pierre Navaro
!> \remark this a refactored version of Bernard Bunner code 
!> (http://cs-www.cs.yale.edu/homes/douglas-craig/mgnet-codes-bunner.html)
module mgd_polar
#include "sll_working_precision.h"
use mpi
use gxch1_2d
use mgd2

implicit none

sll_int32, private :: sx, ex, sy, ey

contains

subroutine initialize_mgd_polar(my_block, my_mg, nerror)

sll_int32       :: nerror
type(block)     :: my_block
type(mg_solver) :: my_mg

sll_int32  :: i,j,k,nxf,nyf,nxm,nym,kps,sxm,exm,sym,eym

!------------------------------------------------------------------------
! set /mgd/ variables to zero
!
do k=1,20
  nxk(k)         = 0
  nyk(k)         = 0
  sxk(k)         = 0
  exk(k)         = 0
  syk(k)         = 0
  eyk(k)         = 0
  kpbgn(k)       = 0
  kcbgn(k)       = 0
  ikdatatype(k)  = MPI_DATATYPE_NULL
  jkdatatype(k)  = MPI_DATATYPE_NULL
  ijkdatatype(k) = MPI_DATATYPE_NULL
  irdatatype(k)  = MPI_DATATYPE_NULL
  jrdatatype(k)  = MPI_DATATYPE_NULL
  ijrdatatype(k) = MPI_DATATYPE_NULL
  sxi(k)         = 0
  exi(k)         = 0
  syi(k)         = 0
  eyi(k)         = 0
  nxr(k)         = 0
  nyr(k)         = 0
  sxr(k)         = 0
  exr(k)         = 0
  syr(k)         = 0
  eyr(k)         = 0
end do
!
if (mod(my_mg%nx,my_mg%nxprocs).ne.0) then
  write(6,100) my_mg%nx,my_mg%nxprocs
  nerror=1
  return
end if
if (mod(my_mg%ny,my_mg%nyprocs).ne.0) then
  write(6,110) my_mg%ny,my_mg%nyprocs
  nerror=1
  return
end if

!
! check that the dimensions are correct
!
i=my_block%ixp*2**(my_block%iex-1)+1
if ((my_mg%nx+1).ne.i) then
  write(6,130) my_mg%nx+1,i
  nerror=1
  return
end if
j=my_block%jyq*2**(my_block%jey-1)+1
if ((my_mg%ny+1).ne.j) then
  write(6,140) my_mg%ny+1,j
  nerror=1
  return
end if
!
! check that the number of points at the coarser level is not smaller
! than the number of processes in either direction
!
if (my_block%ixp.lt.my_mg%nxprocs) then
  write(6,150) my_block%ixp,my_mg%nxprocs
  nerror=1
  return
end if
if (my_block%jyq.lt.my_mg%nyprocs) then
  write(6,160) my_block%jyq,my_mg%nyprocs
  nerror=1
  return
end if
!
! check that coarsifying takes place in all directions at the finest
! grid level
!
if (my_block%ngrid.gt.1) then
  if (my_block%iex.eq.1) then
    write(6,170) my_block%ngrid,my_block%iex
    nerror=1
    return
  end if
  if (my_block%jey.eq.1) then
    write(6,180) my_block%ngrid,my_block%jey
    nerror=1
    return
  end if
end if
!------------------------------------------------------------------------
! define all grid levels
! I have adopted the same notations as in Mudpack as far as possible.
! When a confusion was possible, I added a suffix 'm' to the name
! of the variables. For example, nxm is nx+1 for the multigrid
! code whereas nx means nxp2-2 in the rest of the code.
!
do k=1,my_block%ngrid
  nxk(k)=my_block%ixp*2**(max(k+my_block%iex-my_block%ngrid,1)-1)+1
  nyk(k)=my_block%jyq*2**(max(k+my_block%jey-my_block%ngrid,1)-1)+1
end do

! for all grid levels, set the indices of the subdomain the process
! 'myid' work on, as well as the datatypes needed for the exchange
! of boundary data

nxf=nxk(my_block%ngrid)
nyf=nyk(my_block%ngrid)
sxk(my_block%ngrid)=my_block%sx
exk(my_block%ngrid)=my_block%ex
syk(my_block%ngrid)=my_block%sy
eyk(my_block%ngrid)=my_block%ey
call grid1_type(my_block,ikdatatype(my_block%ngrid), &
                         jkdatatype(my_block%ngrid), &
                         ijkdatatype(my_block%ngrid),&
                         my_block%sx, &
                         my_block%ex, &
                         my_block%sy, &
                         my_block%ey)
do k=my_block%ngrid-1,1,-1
  nxm=nxk(k)
  nym=nyk(k)
  if (nxm.lt.nxf) then
    sxk(k)=sxk(k+1)/2+1
    exk(k)=(exk(k+1)-1)/2+1
  else
    sxk(k)=sxk(k+1)
    exk(k)=exk(k+1)
  end if
  if (nym.lt.nyf) then
    syk(k)=syk(k+1)/2+1
    eyk(k)=(eyk(k+1)-1)/2+1
  else
    syk(k)=syk(k+1)
    eyk(k)=eyk(k+1)
  end if
  nxf=nxm
  nyf=nym
  call grid1_type(my_block,ikdatatype(k),jkdatatype(k),ijkdatatype(k), &
                  sxk(k),exk(k),syk(k),eyk(k))
end do
# if xdebug1
!
! print out the indices and determine the size of the MPI messages
! as a rough check
! 
write(6,*) 'size of the multigrid phi-messages'
do k=my_block%ngrid,1,-1
  call MPI_TYPE_SIZE(ikdatatype(k),nsiz1,ierr)
  call MPI_TYPE_SIZE(jkdatatype(k),nsiz2,ierr)
  call MPI_TYPE_SIZE(ijkdatatype(k),nsiz3,ierr)
  write(6,*) 'myid: ',my_block%id,' k=',k,' sxk=',sxk(k),' exk=', &
                exk(k),' syk=',syk(k),' eyk=',eyk(k), &
                ' size of datatypes: ',nsiz1,nsiz2,nsiz3
end do
# endif
!
! set work space indices for phi, cof at each grid level
! check that there is sufficient work space
!
kps=1
do k=my_block%ngrid,1,-1
  sxm=sxk(k)
  exm=exk(k)
  sym=syk(k)
  eym=eyk(k)
  kpbgn(k)=kps
  kcbgn(k)=kpbgn(k)+(exm-sxm+3)*(eym-sym+3)
  kps=kcbgn(k)+6*(exm-sxm+3)*(eym-sym+3)
end do

!------------------------------------------------------------------------
! Set the boundary values 
! to be used for the Dirichlet boundaries. It is possible to assign
! 4 different constant values to the 4 different sides. The values are
! assigned at the finest grid level, zero is assigned at all levels
! below
! 
! vbc, phibc:
!
!       -----vbc(4)------ 
!       |               |
!       |               |
!     vbc(3)----------vbc(1)
!       |               |
!       |               |
!       ------vbc(2)-----
!
do j=1,4
  my_mg%phibc(j,my_block%ngrid)=my_mg%vbc(j)
end do
do k=my_block%ngrid-1,1,-1
  do j=1,4
    my_mg%phibc(j,k)=0.0d0
  end do
end do

100 format(/,'ERROR in mgdinit: nx=',i3,' is not a multiple of ', &
           'nxprocs=',i3,/,'cannot use the new version of the ',  &
           'multigrid code',/)
110 format(/,'ERROR in mgdinit: ny=',i3,' is not a multiple of ', &
           'nyprocs=',i3,/,'cannot use the new version of the ',  &
           'multigrid code',/)
130  format(/,'ERROR in mgdinit: nxp1=',i3,' <> ixp*2**(iex-1)+1=', &
            i3,/,'-> adjust the multigrid parameters ixp and iex', &
            ' in main',/)
140  format(/,'ERROR in mgdinit: nyp1=',i3,' <> jyq*2**(jey-1)+1=', &
            i3,/,'-> adjust the multigrid parameters jyq and jey', &
            ' in main',/)
170  format(/,'ERROR in mgdinit: ngrid=',i3,' iex=',i3, &
            /,'no coarsifying at the finest grid level in x-direction', &
            /,'this is not allowed by the mutligrid code',/)
180  format(/,'ERROR in mgdinit: ngrid=',i3,' jey=',i3, &
            /,'no coarsifying at the finest grid level in y-direction', &
            /,'this is not allowed by the mutligrid code',/)
150  format(/,'ERROR in mgdinit: ixp=',i3,' < nxprocs=',i3,/, &
            ' there must be at least one grid point at the ', &
            'coarsest grid level',/, &
            '-> increase ixp and decrease iex correspondingly', &
            ' in main',/)
160  format(/,'ERROR in mgdinit: jyq=',i3,' < nyprocs=',i3,/, &
            ' there must be at least one grid point at the ', &
            'coarsest grid level',/, &
            '-> increase jyq and decrease jey correspondingly', &
            ' in main',/)
return
end subroutine

!> Parallel multigrid solver in polar coordinates for the
!> Poisson equation
subroutine mgd_polar_solver(my_block,my_mg,phif,rhsf,r,work, &
                            rro,iter,nprscr,nerror)
sll_real64      :: phif(:,:),rhsf(:,:)
sll_real64      :: r(:)
sll_real64      :: work(*),rro
sll_int32       :: iter,nerror
type(block)     :: my_block
type(mg_solver) :: my_mg

logical    :: nprscr

sll_int32  :: k, icf
sll_real64 :: avo,acorr, relmax
sll_int32  :: sxf,exf,syf,eyf,nxf,nyf
sll_int32  :: ipf,irf,sxm,exm,sym,eym,nxm,nym,ip,ic,kcur

! discretize pde at all levels

sxf=sxk(my_block%ngrid)
exf=exk(my_block%ngrid)
syf=syk(my_block%ngrid)
eyf=eyk(my_block%ngrid)
nxf=nxk(my_block%ngrid)
nyf=nyk(my_block%ngrid)
icf=kcbgn(my_block%ngrid)

call mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,work(icf),r,my_mg%xl,my_mg%yl,my_block%bd)

! determine coefficients at coarser grid levels from
! four neighboring density points; no communication of boundary data
! is involved because of the grid setup, under the condition that
! mod(nx,nxprocs)=mod(ny,nyprocs)=0
!
do k=my_block%ngrid-1,1,-1
   sxm=sxk(k)
   exm=exk(k)
   sym=syk(k)
   eym=eyk(k)
   nxm=nxk(k)
   nym=nyk(k)
   ic=kcbgn(k)
   call polar_phpde(sxm,exm,sym,eym,nxm,nym,work(ic),  &
                    sx,ex,sy,ey,nxf,nyf,r,my_block%bd, &
                    my_mg%xl,my_mg%yl)
end do

!------------------------------------------------------------------------
! set phi,rhsf in work at the finest grid level
!
sxf = sxk(my_block%ngrid)
exf = exk(my_block%ngrid)
syf = syk(my_block%ngrid)
eyf = eyk(my_block%ngrid)
ipf = kpbgn(my_block%ngrid)
irf = kcbgn(my_block%ngrid)+5*(exf-sxf+3)*(eyf-syf+3)
call mgdsetf(sxf,exf,syf,eyf,work(ipf),work(irf),phif,rhsf)
!------------------------------------------------------------------------
! cycling at kcur=ngrid level
!
kcur=my_block%ngrid
do iter=1,my_mg%maxcy
  call mgdkcyc(work,rhsf,kcur,my_mg%kcycle,my_mg%iprer,     &
               my_mg%ipost,my_mg%iresw, my_mg%comm2d,       &
               my_block%id,my_block%neighbor,my_block%bd,   &
               my_mg%phibc)
  sxm=sxk(my_block%ngrid)
  exm=exk(my_block%ngrid)
  sym=syk(my_block%ngrid)
  eym=eyk(my_block%ngrid)
  ip=kpbgn(my_block%ngrid)
  call mgderr(relmax,sxm,exm,sym,eym,phif,work(ip),my_mg%comm2d)
  if (relmax.le.my_mg%tolmax) goto 1000
end do
!------------------------------------------------------------------------
! if not converged in maxcy cycles, issue an error message and quit
!
if (my_block%id.eq.0) write(6,100) my_mg%maxcy,relmax
nerror=1
return
!------------------------------------------------------------------------
! converged
!
1000  continue
!
! rescale phif
!
avo=0.0d0
call gscale(my_block%sx,my_block%ex,my_block%sy,my_block%ey, &
            phif,avo,acorr,my_mg%comm2d,my_mg%nx,my_mg%ny)
!
! exchange boundary data and impose periodic BCs
!
call gxch1lin(phif,my_mg%comm2d,my_block%sx,my_block%ex,     &
              my_block%sy,my_block%ey,my_block%neighbor,     &
              my_block%bd, ikdatatype(my_block%ngrid),       &
              jkdatatype(my_block%ngrid))

call gxch1cor(phif,my_mg%comm2d,sxm,exm,sym,eym,             &
              my_block%neighbor,my_block%bd,                 &
              ijkdatatype(my_block%ngrid))
!
! impose wall and Dirichlet BCs
!
call mgdbdry(sx,ex,sy,ey,phif,my_block%bd,my_mg%vbc)

if (nprscr.and.my_block%id.eq.0) write(6,120) relmax,iter,acorr

return
100   format('WARNING: failed to achieve convergence in ',i5, &
       ' cycles  error=',e12.5)
110     format('  R MGD     err=',e8.3,' iters=',i5,' rcorr=',e9.3)
120     format('  P MGD     err=',e8.3,' iters=',i5,' pcorr=',e9.3)
end subroutine


!------------------------------------------------------------------------
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
#include "sll_working_precision.h"

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

#include "sll_working_precision.h"

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
end module mgd_polar

