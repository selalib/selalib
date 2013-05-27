! *********************************************************
!
! Poisson solver in polar coordinates
!
! **********************************************************
!
! red/black gauss-seidel point relaxation is used along with the
! the default multigrid options.  first mud2cr is called to generate
! a second-order approximation.  then mud24cr is called to improve
! the estimate to fourth-order.
!
! **********************************************************

module sll_mudpack_polar
#include "sll_working_precision.h"
sll_real64 ,allocatable :: work(:)

contains

subroutine initialize_poisson_polar_mudpack(phi, rhs, &
                                            r_min, r_max,  &
                                            theta_min, theta_max, &
                                            nr, nth)
implicit none

! set grid size params

sll_real64, intent(in) :: r_min, r_max
sll_real64, intent(in) :: theta_min, theta_max
sll_int32, intent(in)  :: nr, nth
sll_int32 :: icall
sll_int32 :: iiex,jjey,nnx,nny,llwork
sll_int32, parameter :: iixp = 2 , jjyq = 2

sll_real64, intent(inout) ::  phi(nr,nth)
sll_real64, intent(in) ::  rhs(nr,nth)

! put sll_int32 and floating point argument names in contiguous
! storeage for labelling in vectors iprm,fprm

sll_int32 :: iprm(16),mgopt(4)
sll_real64 :: fprm(6)
sll_int32 intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,&
              iguess,maxcy,method,nwork,lwrkqd,itero
common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax
sll_int32 :: i,j,ierror
sll_real64 :: dlx,dly

equivalence(intl,iprm)
equivalence(xa,fprm)

! declare coefficient and boundary condition input subroutines external
external cofcr,bndcr


nnx = nr
nny = nth

! set minimum required work space
llwork=(7*(nnx+2)*(nny+2)+44*nnx*nny)/3
      
if (.not. allocated(work)) then
   allocate(work(llwork))
   icall = 0
else
   icall = 1
end if
iiex = ceiling(log((nnx-1.)/iixp)/log(2.))+1
jjey = ceiling(log((nny-1.)/jjyq)/log(2.))+1

! set input sll_int32 arguments
intl = 0

! set boundary condition flags
nxa = 1
nxb = 1
nyc = 0
nyd = 0

! set grid sizes from parameter statements
ixp = iixp
jyq = jjyq
iex = iiex
jey = jjey

nx = ixp*(2**(iex-1))+1
ny = jyq*(2**(jey-1))+1

if (nx /= nr) then
   print*, "nx,nr=", nx, nr
   stop ' nx different de nr dans mg_polar_poisson'
end if
if (ny /= nth) then
   print*, "ny,nth=", ny, nth
   stop ' ny different de nth dans mg_polar_poisson'
end if

! set multigrid arguments (w(2,1) cycling with fully weighted
! residual restriction and cubic prolongation)
mgopt(1) = 2
mgopt(2) = 2
mgopt(3) = 1
mgopt(4) = 3

! set three cycles to ensure second-order approx
maxcy = 3

! set no initial guess forcing full multigrid cycling
iguess = 0

! set work space length approximation from parameter statement
nwork = llwork

! set point relaxation
method = 0

! set mesh increments
xa = r_min
xb = r_max
yc = theta_min
yd = theta_max
dlx = (xb-xa)/float(nx-1)
dly = (yd-yc)/float(ny-1)

! set for no error control flag
tolmax = 0.0

! set specified boundaries in phi at x=xa 
do j=1,ny
   phi(1,j) = 0.0
   phi(nx,j) = 0.0
end do
write(*,100)

write(*,102) (mgopt(i),i=1,4)
write(*,103) xa,xb,yc,yd,tolmax
write(*,104) intl
call mud2cr(iprm,fprm,work,cofcr,bndcr,rhs,phi,mgopt,ierror)
write (*,200) ierror,iprm(16)
if (ierror.gt.0) call exit(0)

100 format(//' multigrid poisson solver in polar coordinates ')
    write (*,101) (iprm(i),i=1,15)
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
end subroutine initialize_poisson_polar_mudpack


subroutine solve_poisson_polar_mudpack(phi, rhs)
implicit none

! set grid size params

sll_int32 :: icall
sll_int32 :: iiex,jjey,nnx,nny,llwork
sll_int32, parameter :: iixp = 2 , jjyq = 2

sll_real64, intent(inout) ::  phi(nr,nth)
sll_real64, intent(in) ::  rhs(nr,nth)

! put sll_int32 and floating point argument names in contiguous
! storeage for labelling in vectors iprm,fprm

sll_int32 :: iprm(16),mgopt(4)
sll_real64 :: fprm(6)
sll_int32 intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,&
              iguess,maxcy,method,nwork,lwrkqd,itero
common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax
sll_int32 :: i,j,ierror
sll_real64 :: dlx,dly

equivalence(intl,iprm)
equivalence(xa,fprm)

! declare coefficient and boundary condition input subroutines external
external cofcr,bndcr

icall = 1
! attempt solution
intl = 1
write(*,106) intl,method,iguess
call mud2cr(iprm,fprm,work,cofcr,bndcr,rhs,phi,mgopt,ierror)
write (*,107) ierror
if (ierror.gt.0) call exit(0)

! attempt fourth order approximation
call mud24cr(work,cofcr,bndcr,phi,ierror)
write (*,108) ierror
if (ierror.gt.0) call exit(0)

106 format(/' approximation call to mud2cr', &
    /' intl = ',i2, ' method = ',i2,' iguess = ',i2)
107 format(' ierror = ',i2)
108 format(/' mud24cr test', '  ierror = ',i2)

return
end subroutine solve_poisson_polar_mudpack

end module sll_mudpack_polar


!> input pde coefficients at any grid point (x,y) in the solution region
!> (xa.le.x.le.xb,yc.le.y.le.yd) to mud2cr
subroutine cofcr(x,y,cxx,cxy,cyy,cx,cy,ce)
implicit none
real(8) :: x,y,cxx,cxy,cyy,cx,cy,ce
cxx = 1.0 
cxy = 0.0 
cyy = 1.0 / (x*x) 
cx  = 1.0 / x 
cy  = 0.0 
ce  = 0.0 
return
end subroutine

!> input mixed "oblique" derivative b.c. to mud2cr
!> at upper y boundary
subroutine bndcr(kbdy,xory,alfa,beta,gama,gbdy)
implicit none
integer  :: kbdy
real(8)  :: xory,alfa,beta,gama,gbdy

if (kbdy.eq.2) then

   ! x=xb boundary.
   ! b.c. has the form alfyd(x)*px+betyd(x)*py+gamyd(x)*pe = gbdyd(x)
   ! where x = yorx.   alfa,beta,gama,gbdy corresponding to alfyd(x),
   ! betyd(x),gamyd(x),gbdyd(y) must be output.

   alfa = 1.0
   beta = 0.0
   gama = 0.0
   gbdy = 0.0

end if

return
end subroutine
