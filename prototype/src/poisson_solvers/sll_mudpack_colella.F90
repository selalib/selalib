! *********************************************************
!
! Poisson solver in colella coordinates
!
! **********************************************************
!
! red/black gauss-seidel point relaxation is used along with the
! the default multigrid options.  first mud2cr is called to generate
! a second-order approximation.  then mud24cr is called to improve
! the estimate to fourth-order.
!
! **********************************************************

module sll_mudpack_colella
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_mudpack_base
implicit none

contains

subroutine initialize_poisson_colella_mudpack(this, phi, rhs, &
                                              eta1_min, eta1_max,  &
                                              eta2_min, eta2_max, &
                                              nc_eta1, nc_eta2)
implicit none

! set grid size params

type(mudpack_2d) :: this
sll_real64, intent(in) :: eta1_min, eta1_max
sll_real64, intent(in) :: eta2_min, eta2_max
sll_int32, intent(in)  :: nc_eta1, nc_eta2
sll_int32 :: icall, error
sll_int32 :: iiex,jjey,nnx,nny,llwork
sll_int32, parameter :: iixp = 2 , jjyq = 2

sll_real64, intent(inout) ::  phi(nc_eta1+1,nc_eta2+1)
sll_real64, intent(in) ::  rhs(nc_eta1+1,nc_eta2+1)

! put sll_int32 and floating point argument names in contiguous
! storeage for labelling in vectors iprm,fprm

sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax
sll_int32  :: i,j,ierror
sll_real64 :: dlx,dly
sll_int32  :: iprm(16)
sll_real64 :: fprm(6)

equivalence(intl,iprm)
equivalence(xa,fprm)

! declare coefficient and boundary condition input subroutines external
external coef,bnd

nnx = nc_eta1+1
nny = nc_eta2+1

! set minimum required work space
llwork=(7*(nnx+2)*(nny+2)+44*nnx*nny)/3
      
SLL_ALLOCATE(this%work(llwork),error)
icall = 0
iiex = ceiling(log((nnx-1.)/iixp)/log(2.))+1
jjey = ceiling(log((nny-1.)/jjyq)/log(2.))+1

! set input sll_int32 arguments
intl = 0

! set boundary condition flags
nxa = 1
nxb = 1
nyc = 1
nyd = 1

! set grid sizes from parameter statements
ixp = iixp
jyq = jjyq
iex = iiex
jey = jjey

nx = ixp*(2**(iex-1))+1
ny = jyq*(2**(jey-1))+1

if (nx /= nc_eta1+1) then
   print*, "nx,nc_eta1=", nx, nc_eta1+1
   stop ' nx different de nc_eta1+1 '
end if
if (ny /= nc_eta2+1) then
   print*, "ny,nc_eta2+1=", ny, nc_eta2+1
   stop ' ny different de nc_eta2+1 '
end if

! set multigrid arguments (w(2,1) cycling with fully weighted
! residual restriction and cubic prolongation)
this%mgopt(1) = 2
this%mgopt(2) = 2
this%mgopt(3) = 1
this%mgopt(4) = 3

! set three cycles to ensure second-order approx
maxcy = 3

! set no initial guess forcing full multigrid cycling
iguess = 0

! set work space length approximation from parameter statement
nwork = llwork

! set point relaxation
method = 0

! set mesh increments
xa = eta1_min
xb = eta1_max
yc = eta2_min
yd = eta2_max
dlx = (xb-xa)/float(nc_eta1)
dly = (yd-yc)/float(nc_eta2)

! set for no error control flag
tolmax = 0.0

write(*,100)
write(*,101) (iprm(i),i=1,15)
write(*,102) (this%mgopt(i),i=1,4)
write(*,103) xa,xb,yc,yd,tolmax
write(*,104) intl
call mud2cr(iprm,fprm,this%work,coef,bnd,rhs,phi,this%mgopt,ierror)
write (*,200) ierror,iprm(16)
if (ierror > 0) call exit(0)

100 format(//' multigrid poisson solver in colella mesh ')
101 format(/' integer input arguments ', &
    /'intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2, &
    /' ixp = ',i2,' jyq = ',i2,' iex = ',i2,' jey = ',i2 &
    /' nx = ',i3,' ny = ',i3,' iguess = ',i2,' maxcy = ',i2, &
    /' method = ',i2, ' work space estimate = ',i7)
102 format(/' multigrid option arguments ',&
    /' kcycle = ',i2,  &
    /' iprer = ',i2,   &
    /' ipost = ',i2    &
    /' intpol = ',i2)
103 format(/' floating point input parameters ', &
    /' xa = ',f6.3,' xb = ',f6.3,' yc = ',f6.3,' yd = ',f6.3, &
    /' tolerance (error control) =   ',e10.3)
104 format(/' discretization call to mud2cr', ' intl = ', i2)
200 format(' ierror = ',i2, ' minimum work space = ',i7)

return
end subroutine initialize_poisson_colella_mudpack


subroutine solve_poisson_colella_mudpack(this, phi, rhs)
implicit none

! set grid size params

type(mudpack_2d) :: this
sll_int32        :: icall

sll_real64, intent(inout) ::  phi(:,:)
sll_real64, intent(inout) ::  rhs(:,:)

! put sll_int32 and floating point argument names in contiguous
! storeage for labelling in vectors iprm,fprm

sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
sll_int32  :: ierror
sll_int32  :: iprm(16)
sll_real64 :: fprm(6)

common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
                iguess,maxcy,method,nwork,lwrkqd,itero
common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax

equivalence(intl,iprm)
equivalence(xa,fprm)

! declare coefficient and boundary condition input subroutines external
external coef,bnd

iguess = 1
icall  = 1
intl   = 1
write(*,106) intl,method,iguess
! attempt solution
call mud2cr(iprm,fprm,this%work,coef,bnd,rhs,phi,this%mgopt,ierror)
SLL_ASSERT(ierror == 0)
! attempt fourth order approximation
call mud24cr(this%work,coef,bnd,phi,ierror)
SLL_ASSERT(ierror == 0)

106 format(/' approximation call to mud2cr', &
    /' intl = ',i2, ' method = ',i2,' iguess = ',i2)

return
end subroutine solve_poisson_colella_mudpack

end module sll_mudpack_colella


