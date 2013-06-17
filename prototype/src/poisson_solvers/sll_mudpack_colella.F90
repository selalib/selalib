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
implicit none

sll_real64 ,allocatable :: work(:)

contains

subroutine initialize_poisson_colella_mudpack(phi, rhs, &
                                              eta1_min, eta1_max,  &
                                              eta2_min, eta2_max, &
                                              nc_eta1, nc_eta2)
implicit none

! set grid size params

sll_real64, intent(in) :: eta1_min, eta1_max
sll_real64, intent(in) :: eta2_min, eta2_max
sll_int32, intent(in)  :: nc_eta1, nc_eta2
sll_int32 :: icall
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
sll_int32  :: mgopt(4)

equivalence(intl,iprm)
equivalence(xa,fprm)

! declare coefficient and boundary condition input subroutines external
external coef,bnd

nnx = nc_eta1+1
nny = nc_eta2+1

! set minimum required work space
llwork=(7*(nnx+2)*(nny+2)+44*nnx*nny)/3
      
allocate(work(llwork))
icall = 0
iiex = ceiling(log((nnx-1.)/iixp)/log(2.))+1
jjey = ceiling(log((nny-1.)/jjyq)/log(2.))+1

! set input sll_int32 arguments
intl = 0

! set boundary condition flags
nxa = 0
nxb = 0
nyc = 0
nyd = 0

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
xa = eta1_min
xb = eta1_max
yc = eta2_min
yd = eta2_max
dlx = (xb-xa)/float(nc_eta1)
dly = (yd-yc)/float(nc_eta2)

! set for no error control flag
tolmax = 0.0

write(*,100)
write(*,102) (mgopt(i),i=1,4)
write(*,103) xa,xb,yc,yd,tolmax
write(*,104) intl
call mud2cr(iprm,fprm,work,coef,bnd,rhs,phi,mgopt,ierror)
write (*,200) ierror,iprm(16)
if (ierror > 0) call exit(0)

100 format(//' multigrid poisson solver in colella mesh ')
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
end subroutine initialize_poisson_colella_mudpack


subroutine solve_poisson_colella_mudpack(phi, rhs)
implicit none

! set grid size params

sll_int32 :: icall

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
sll_int32  :: mgopt(4)

common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
                iguess,maxcy,method,nwork,lwrkqd,itero
common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax

equivalence(intl,iprm)
equivalence(xa,fprm)

! declare coefficient and boundary condition input subroutines external
external coef,bnd

iguess = 0
icall  = 1
intl   = 1
write(*,106) intl,method,iguess
! attempt solution
call mud2cr(iprm,fprm,work,coef,bnd,rhs,phi,mgopt,ierror)
SLL_ASSERT(ierror == 0)
! attempt fourth order approximation
call mud24cr(work,coef,bnd,phi,ierror)
SLL_ASSERT(ierror == 0)

106 format(/' approximation call to mud2cr', &
          /' intl = ',i2, ' method = ',i2,' iguess = ',i2)

return
end subroutine solve_poisson_colella_mudpack

end module sll_mudpack_colella


!> input pde coefficients at any grid point (x,y) in the solution region
!> (xa.le.x.le.xb,yc.le.y.le.yd) to mud2cr
!subroutine coef(x,y,cxx,cxy,cyy,cx,cy,ce)
!
!implicit none
!
!real(8):: x,y,cxx,cxy,cyy,cx,cy,ce
!real(8) :: alpha,det,a11,a12,a21,a22,a11_x,a11_y,a12_x
!real(8) :: a21_y,a22_x,a22_y,det_x,det_y
!real(8) :: sll_pi
!
!!! CHANGE HERE FOR COLELLA MESH
!sll_pi = 3.14159265359
!alpha=0.0_8
!
!xx = x + alpha * sin(2.0_8*sll_pi*x)*sin(2.0_8*sll_pi*y)
!yy = y + alpha * sin(2.0_8*sll_pi*x)*sin(2.0_8*sll_pi*y)
!
!xx_x = 1 + alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*x)*sin(2.0_8*sll_pi*y)
!xx_y = alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*x)*cos(2.0_8*sll_pi*y)
!yy_x = alpha*2.0_8*sll_pi*cos(2.0_8*sll_pi*x)*sin(2.0_8*sll_pi*y)
!yy_y = 1+alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*x)*cos(2.0_8*sll_pi*y)
!
!
!
!
!
!a11 = 
!    
!a11_x= 
!
!a11_y= 
!    
!a12 = 
!a12_x= 
!a21 = a12
!
!a21_y=
!
!a22 = 
!    
!a22_x =
!a22_y = 
!
!det = a11*a22-a12*a21
!
!det_x = a11_x * a22 + a11 * a22_x - 2.0_8 * a12_x * a12
!
!det_y = a11_y * a22 + a11 * a22_y - 2.0_8 * a21_y * a21
!
!a11 = a11/det
!a22 = a22/det
!a12 = -a12/det
!a21 = -a21/det
!
!a11_x = (a11_x*det - a11 * det_x)/det**2
!a12_x =-(a12_x*det - a12 * det_x)/det**2
!a21_y =-(a21_y*det - a21 * det_y)/det**2
!a22_y = (a22_y*det - a22 * det_y)/det**2
!
!cxx = -a11
!cyy = -a22
!cxy= -a12-a21
!cx  = -(a11_x+a21_y) 
!cy  = -(a12_x+a22_y)
!ce  = 0.0_8
!
!
!return
!end subroutine
!subroutine coef(x,y,cxx,cxy,cyy,cx,cy,ce)

implicit none

real(8):: x,y,cxx,cxy,cyy,cx,cy,ce
real(8) :: alpha,det,a11,a12,a21,a22,a11_x,a11_y,a12_x
real(8) :: a21_y,a22_x,a22_y,det_x,det_y
real(8) :: x_min,x_max,y_min,y_max,xx,yy,sll_pi

!! CHANGE HERE FOR COLELLA MESH
sll_pi = 3.14159265359
alpha=0.0_8
x_min=0.0_8
x_max=2.0_8 * sll_pi
y_min=0.0_8
y_max=2.0_8 * sll_pi

xx = (x-x_min)/(x_max-x_min)
yy = (y-y_min)/(y_max-y_min)

a11 = (alpha*2*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2*sll_pi*yy)/(y_max-y_min))**2 & 
    & + (1+alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min))**2
    
a11_x= 2.0_8*(alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min))* &
      &(alpha*cos(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy))*(2.0_8*sll_pi)**2/((x_max-x_min)*(y_max-y_min)) + &
      & 2.0_8*(1+alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min))* &
      & (alpha*cos(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy))*(2.0_8*sll_pi)**2/((x_max-x_min)*(y_max-y_min))

a11_y= 2.0_8*(alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min))* &
      &(alpha*cos(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy))*(2.0_8*sll_pi)**2/((x_max-x_min)*(y_max-y_min)) + &
      & 2.0_8*(1+alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min))* &
      & (alpha*cos(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy))*(2.0_8*sll_pi)**2/((x_max-x_min)*(y_max-y_min))
    
a12 = (1+alpha*2.0_8*sll_pi*cos(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy)/(x_max-x_min))* &
    & (alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min))+ &
    & (alpha*2.0_8*sll_pi*cos(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy)/(x_max-x_min))* &
    & (1+alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min))
    
a12_x= -( alpha*sin(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy)*(2.0_8*sll_pi/(x_max-x_min))**2)* &
    & (alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min)) + &
    & (1+alpha*2.0_8*sll_pi*cos(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy)/(x_max-x_min)) * &
    & (alpha*cos(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy))*(2*sll_pi)**2/((x_max-x_min)*(y_max-y_min))- &
    & (alpha*sin(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy)*(2.0_8*sll_pi/(x_max-x_min))**2)* &
    & (1+alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min)) + &
    & (alpha*2.0_8*sll_pi*cos(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy)/(x_max-x_min)) * &
    & ( alpha*cos(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy))*(2.0_8*sll_pi)**2.0_8/((x_max-x_min)*(y_max-y_min))  

a21 = a12

a21_y=-(alpha*cos(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy))*(2.0_8*sll_pi)**2/((x_max-x_min)*(y_max-y_min))* &
    &  (alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min))+ &
    & (1+alpha*2.0_8*sll_pi*cos(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy)/(x_max-x_min)) * &
    & (alpha*sin(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy))*(2.0_8*sll_pi/(y_max-y_min))**2- &
    & (alpha*cos(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy))*(2.0_8*sll_pi)**2/((x_max-x_min)*(y_max-y_min))* &
    & (1+alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*xx)*cos(2.0_8*sll_pi*yy)/(y_max-y_min)) + &
    & (alpha*2.0_8*sll_pi*cos(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy)/(x_max-x_min)) * &
    & (alpha*sin(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy))*(2.0_8*sll_pi/(y_max-y_min))**2

a22 = (1+alpha*2.0_8*sll_pi*cos(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy)/(x_max-x_min))**2 & 
    & + (alpha*2.0_8*sll_pi*cos(2.0_8*sll_pi*xx)*sin(2.0_8*sll_pi*yy)/(x_max-x_min))**2
    
a22_x = 2.0_8*(alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*yy)*cos(2.0_8*sll_pi*xx)/(x_max-x_min))* &
      &(alpha*cos(2.0_8*sll_pi*yy)*cos(2.0_8*sll_pi*xx))*(2.0_8*sll_pi)**2/((x_max-x_min)*(y_max-y_min)) + &
      & 2.0_8*(1+alpha*2.0_8*sll_pi*sin(2*sll_pi*yy)*cos(2*sll_pi*xx)/(x_max-x_min))* &
      & (alpha*cos(2*sll_pi*yy)*cos(2.0_8*sll_pi*xx))*(2.0_8*sll_pi)**2/((x_max-x_min)*(y_max-y_min))
    
a22_y = 2.0_8*(alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*yy)*cos(2.0_8*sll_pi*xx)/(x_max-x_min))* &
      &(alpha*cos(2.0_8*sll_pi*yy)*cos(2.0_8*sll_pi*xx))*(2.0_8*sll_pi)**2/((x_max-x_min)*(y_max-y_min)) + &
      & 2.0_8*(1+alpha*2.0_8*sll_pi*sin(2.0_8*sll_pi*yy)*cos(2.0_8*sll_pi*xx)/(x_max-x_min))* &
      & (alpha*cos(2.0_8*sll_pi*yy)*cos(2.0_8*sll_pi*xx))*(2.0_8*sll_pi)**2/((x_max-x_min)*(y_max-y_min))   
    
det = a11*a22-a12*a21

det_x = a11_x * a22 + a11 * a22_x - 2.0_8 * a12_x * a12

det_y = a11_y * a22 + a11 * a22_y - 2.0_8 * a21_y * a21

a11 = a11/det
a22 = a22/det
a12 = -a12/det
a21 = -a21/det

a11_x = (a11_x*det - a11 * det_x)/det**2
a12_x =-(a12_x*det - a12 * det_x)/det**2
a21_y =-(a21_y*det - a21 * det_y)/det**2
a22_y = (a22_y*det - a22 * det_y)/det**2

cxx = -a11
cyy = -a22
cxy= -a12-a21
cx  = -(a11_x+a21_y) 
cy  = -(a12_x+a22_y)
ce  = 0.0_8


return
end subroutine

!> at upper y boundary
subroutine bnd(kbdy,xory,alfa,beta,gama,gbdy)
implicit none
integer  :: kbdy
real(8)  :: xory,alfa,beta,gama,gbdy

!! Set bounday condition value
if (kbdy.eq.2) then

   ! x=xb boundary.
   ! b.c. has the form alfyd(x)*px+betyd(x)*py+gamyd(x)*pe = gbdyd(x)
   ! where xory = x.   alfa,beta,gama,gbdy corresponding to alfyd(x),
   ! betyd(x),gamyd(x),gbdyd(y) must be output.

   alfa = 1.0
   beta = 0.0
   gama = 0.0
   gbdy = 0.0

end if

return
end subroutine
