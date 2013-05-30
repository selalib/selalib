module sll_mudpack_cartesian
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_file_io.h"

implicit none

sll_real64, dimension(:), allocatable :: work
sll_int32 :: mgopt(4)

enum, bind(C)
   enumerator :: PERIODIC=0, DIRICHLET=1
   enumerator :: NEUMANN_RIGHT=2, NEUMANN_LEFT=3, NEUMANN=4
end enum

contains

subroutine initialize_mudpack_cartesian(phi, rhs,                    &
                                        eta1_min, eta1_max, nc_eta1, &
                                        eta2_min, eta2_max, nc_eta2, &
                                        bc_eta1_left, bc_eta1_right, &
                                        bc_eta2_left, bc_eta2_right )
implicit none

! set grid size params
sll_real64, intent(in)  :: eta1_min, eta1_max
sll_real64, intent(in)  :: eta2_min, eta2_max
sll_int32,  intent(in)  :: nc_eta1, nc_eta2
sll_int32,  intent(in)  :: bc_eta1_left
sll_int32,  intent(in)  :: bc_eta1_right
sll_int32,  intent(in)  :: bc_eta2_left
sll_int32,  intent(in)  :: bc_eta2_right
sll_int32,  parameter   :: iixp = 2 , jjyq = 2
sll_int32 :: icall, iiex, jjey, nnx, nny, llwork

sll_real64, intent(inout) :: phi(nc_eta1+1,nc_eta2+1)
sll_real64, intent(inout) :: rhs(nc_eta1+1,nc_eta2+1)

!put integer and floating point argument names in contiguous
!storeage for labelling in vectors iprm,fprm
sll_int32  :: iprm(16)
sll_real64 :: fprm(6)
sll_int32  :: i,j,error
sll_real64 :: dlx,dly,x,y,cxx,cyy,cx,cy,ce,pxx,pyy,px,py,pe,errmax
sll_real64 :: cex,cey

sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
common/itmud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

equivalence(intl,iprm)
equivalence(xa,fprm)

!declare coefficient and boundary condition input subroutines external
external cofx,cofy,bndsp

nnx = nc_eta1+1
nny = nc_eta2+1

! set minimum required work space
llwork=(7*(nnx+2)*(nny+2)+44*nnx*nny)/3
      
allocate(work(llwork))
icall = 0
iiex = ceiling(log((nnx-1.)/iixp)/log(2.))+1
jjey = ceiling(log((nny-1.)/jjyq)/log(2.))+1

! set multigrid arguments (w(2,1) cycling with fully weighted
! residual restriction and cubic prolongation)
mgopt(1) = 2
mgopt(2) = 2
mgopt(3) = 1

!set input integer arguments
intl = 0

!set boundary condition flags
nxa  = bc_eta1_left
nxb  = bc_eta1_right
nyc  = bc_eta2_left
nyd  = bc_eta2_right

!set grid sizes from parameter statements
ixp  = iixp
jyq  = jjyq
iex  = iiex
jey  = jjey

nx = ixp*(2**(iex-1))+1
ny = jyq*(2**(jey-1))+1

if (nx /= nc_eta1+1 .or. ny /= nc_eta2+1) then
   print*, "nx,nc_eta1+1=", nx, nc_eta1+1
   print*, "ny,nc_eta2+1=", ny, nc_eta2+1
   stop ' nx or ny different in sll_mudpack_cartesian '
end if

!set multigrid arguments (w(2,1) cycling with fully weighted
!residual restriction and cubic prolongation)
mgopt(1) = 2
mgopt(2) = 2
mgopt(3) = 1
mgopt(4) = 3

!set for three cycles to ensure second-order approximation is computed
maxcy = 3

!set no initial guess forcing full multigrid cycling
iguess = 0

!set work space length approximation from parameter statement
nwork = llwork

!set point relaxation
method = 0

!set end points of solution rectangle in (x,y) space
xa = eta1_min
xb = eta1_max
yc = eta2_min
yd = eta2_max

!set for no error control flag
tolmax = 0.0

write(*,101) (iprm(i),i=1,15)
write(*,102) (mgopt(i),i=1,4)
write(*,103) xa,xb,yc,yd,tolmax
write(*,104) intl

call mud2sp(iprm,fprm,work,cofx,cofy,bndsp,rhs,phi,mgopt,error)

!print error parameter and minimum work space requirement
write (*,200) error,iprm(16)
if (error > 0) call exit(0)

101 format(/' integer input arguments ', &
    &/'intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2, &
    &/' ixp = ',i2,' jyq = ',i2,' iex = ',i2,' jey = ',i2 &
    &/' nx = ',i3,' ny = ',i3,' iguess = ',i2,' maxcy = ',i2, &
    &/' method = ',i2, ' work space estimate = ',i7)
102 format(/' multigrid option arguments ', &
    &/' kcycle = ',i2, &
    &/' iprer = ',i2, &
    &/' ipost = ',i2, &
    &/' intpol = ',i2)
103 format(/' floating point input parameters ', &
    &/' xa = ',f6.3,' xb = ',f6.3,' yc = ',f6.3,' yd = ',f6.3, &
    &/' tolerance (error control) =   ',e10.3)
104 format(/' discretization call to mud2sp', ' intl = ', i2)
200 format(' error = ',i2, ' minimum work space = ',i7)
     
end subroutine initialize_mudpack_cartesian


subroutine solve_mudpack_cartesian(phi, rhs)

sll_real64 :: phi(:,:),rhs(:,:)

!put integer and floating point argument names in contiguous
!storeage for labelling in vectors iprm,fprm
sll_int32  :: iprm(16)
sll_real64 :: fprm(6)
sll_int32  :: error

sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
common/itmud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

equivalence(intl,iprm)
equivalence(xa,fprm)

!declare coefficient and boundary condition input subroutines external
external cofx,cofy,bndsp

!set no initial guess forcing full multigrid cycling
iguess = 0

!attempt solution
intl = 1
write(*,106) intl,method,iguess

call mud2sp(iprm,fprm,work,cofx,cofy,bndsp,rhs,phi,mgopt,error)
write(*,107) error
if (error > 0) call exit(0)

! attempt to improve approximation to fourth order
call mud24sp(work,phi,error)

write (*,108) error
if (error > 0) call exit(0)

106 format(/' approximation call to mud2sp', &
    &/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
107 format(' error = ',i2)
108 format(/' mud24sp test ', ' error = ',i2)
     
end subroutine solve_mudpack_cartesian

end module sll_mudpack_cartesian




!> input x dependent coefficients
subroutine cofx(x,cxx,cx,cex)
implicit none
real(8)  :: x,cxx,cx,cex
cxx = 1.0
cx  = 0.0
cex = 0.0
return
end

!> input y dependent coefficients
subroutine cofy(y,cyy,cy,cey)
implicit none
real(8)  :: y,cyy,cy,cey
cyy = 1.0
cy  = 0.0
cey = 0.0
return
end

!> input mixed derivative b.c. to mud2sp
subroutine bndsp(kbdy,xory,alfa,gbdy)
implicit none
integer  :: kbdy
real(8)  :: xory,alfa,gbdy,x,y,pe,px,py,pxx,pyy
real(8)  :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

if (kbdy == 1) then  ! x=xa boundary
   y = xory
   x = xa
   call exact(x,y,pxx,pyy,px,py,pe)
   alfa = -1.0
   gbdy = px + alfa*pe
   return
end if

if (kbdy == 4) then  ! y=yd boundary
   y = yd
   x = xory
   call exact(x,y,pxx,pyy,px,py,pe)
   alfa = 1.0
   gbdy = py + alfa*pe
   return
end if
end

