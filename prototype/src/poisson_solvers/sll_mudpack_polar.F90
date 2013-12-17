!>
!> Poisson solver in polar coordinates using mudpack library
!>
!> red/black gauss-seidel point relaxation is used along with the
!> the default multigrid options.  first mud2cr is called to generate
!> a second-order approximation.  then mud24cr is called to improve
!> the estimate to fourth-order.
!>
module sll_mudpack_polar
#include "sll_working_precision.h"
#include "sll_assert.h"

use sll_mudpack_base

implicit none

contains

!> Initialize the Poisson solver in polar coordinates using MUDPACK
!> library
subroutine initialize_poisson_polar_mudpack(this,                      &
                                            r_min, r_max, nr,          &
                                            theta_min, theta_max, nth, &
                                            bc_r_min, bc_r_max,        &
                                            bc_theta_min, bc_theta_max )
implicit none

type(mudpack_2d) :: this            !< Solver object
sll_real64, intent(in) :: r_min     !< radius min
sll_real64, intent(in) :: r_max     !< radius min
sll_real64, intent(in) :: theta_min !< theta min
sll_real64, intent(in) :: theta_max !< theta max
sll_int32, intent(in)  :: nr        !< radius number of points
sll_int32, intent(in)  :: nth       !< theta number of points
sll_int32 :: icall
sll_int32 :: iiex,jjey,llwork
sll_int32 :: bc_r_min               !< left boundary condition r
sll_int32 :: bc_r_max               !< right boundary condition r
sll_int32 :: bc_theta_min           !< left boundary condition theta
sll_int32 :: bc_theta_max           !< right boundary condition theta

sll_real64 ::  phi(nr,nth)          !< electric potential
sll_real64 ::  rhs(nr,nth)          !< charge density

! put sll_int32 and floating point argument names in contiguous
! storeage for labelling in vectors iprm,fprm

sll_int32 :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,&
              iguess,maxcy,method,nwork,lwrkqd,itero
common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax
sll_int32  :: i,j,ierror
sll_int32  :: iprm(16)
sll_real64 :: fprm(6)

equivalence(intl,iprm)
equivalence(xa,fprm)

! declare coefficient and boundary condition input subroutines external
external coef_polar,bndcr

nx = nr
ny = nth

! set minimum required work space
llwork=(7*(nx+2)*(ny+2)+44*nx*ny)/3
      
allocate(this%work(llwork))
icall = 0

! set input sll_int32 arguments
intl = 0

! set boundary condition flags
nxa = bc_r_min
nxb = bc_r_max
nyc = bc_theta_min
nyd = bc_theta_max

! set grid sizes from parameter statements
ixp = 2
jyq = 2
iex = ceiling(log((nx-1.)/ixp)/log(2.))+1
jey = ceiling(log((ny-1.)/jyq)/log(2.))+1

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
xa = r_min
xb = r_max
yc = theta_min
yd = theta_max

! set for no error control flag
tolmax = 0.0

write(*,100)
write(*,101) (iprm(i),i=1,15)
write(*,102) (this%mgopt(i),i=1,4)
write(*,103) xa,xb,yc,yd,tolmax
write(*,104) intl

call mud2cr(iprm,fprm,this%work,coef_polar,bndcr,rhs,phi,this%mgopt,ierror)

write (*,200) ierror,iprm(16)
if (ierror > 0) call exit(0)

100 format(//' multigrid poisson solver in polar coordinates ')
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


!> Solve the Poisson equation and get the potential
subroutine solve_poisson_polar_mudpack(this, phi, rhs)
implicit none

! set grid size params
type(mudpack_2d) :: this  !< solver data object
sll_int32 :: icall
sll_int32, parameter :: iixp = 2 , jjyq = 2

sll_real64, intent(inout) ::  phi(:,:) !< electric potential
sll_real64, intent(inout) ::  rhs(:,:) !< charge density

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
external coef_polar,bndcr

icall = 1
intl  = 1
write(*,106) intl,method,iguess
! attempt solution
call mud2cr(iprm,fprm,this%work,coef_polar,bndcr,rhs,phi,this%mgopt,ierror)
SLL_ASSERT(ierror == 0)
! attempt fourth order approximation
call mud24cr(this%work,coef_polar,bndcr,phi,ierror)
SLL_ASSERT(ierror == 0)

106 format(/' approximation call to mud2cr', &
    /' intl = ',i2, ' method = ',i2,' iguess = ',i2)

return
end subroutine solve_poisson_polar_mudpack

end module sll_mudpack_polar


!> input pde coefficients at any grid point (x,y) in the solution region
!> (xa.le.x.le.xb,yc.le.y.le.yd) to mud2cr
subroutine coef_polar(x,y,cxx,cxy,cyy,cx,cy,ce)
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
