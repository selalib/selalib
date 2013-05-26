!
!     a sample program/test driver for mudpack is below. it can be
!     executed as an initial test.  output is listed for the case
!     described.
!
!     test the driver below by solving the elliptic pde
!
!          (1.+y**2)*pxx + exp(-(x+y))*(pyy-py) - (x+y)*pe = r(x,y)
!
!     on the unit square with specified boundary conditions at
!     xb = 1.0, yc = 0.0 and mixed boundary conditions
!
!          dp/dx - y*p(xa,y) = g(y)  (at x=0)
!
!     and
!
!          dp/dy + x*p(x,yd) = h(x)  (at y=1).
!
!     use line relaxation in the y direction and choose a grid as close
!     to 50 by 100 as the grid size parameters allow. use the exact
!     solution
!
!          pe(x,y) = x**5 + y**5 + 1.0
!
!     for testing.
!
!     the default multigrid options with no initial guess and two
!     cycles are used when first calling mud2 to yield a second-
!     order estimate.  then mud24 is called to produce a fourth-order
!     approximation.
!
!
! ************************************************************
!
!
program test_mudpack_2d
implicit none
integer  :: iixp,jjyq,iiex,jjey,nnx,nny,llwork
!
!     set grid sizes with parameter statements
!
parameter (iixp = 3 , jjyq = 3 , iiex = 5, jjey = 6)
parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
!
!     set minimal required work space length (see tmud2.f)
!
parameter (llwork=70048)
real(8) :: phi(nnx,nny),rhs(nnx,nny),work(llwork)
!
!     put integer and floating point argument names in contiguous
!     storeage for labelling in vectors iprm,fprm
!
integer :: iprm(16),mgopt(4)
real(8) :: fprm(6)
integer :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
common/itmud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
real(8) :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2/xa,xb,yc,yd,tolmax,relmax
equivalence(intl,iprm)
equivalence(xa,fprm)
integer :: i,j,ierror
real(8) :: dlx,dly,x,y,cxx,cyy,cx,cy,ce,pxx,pyy,px,py,pe,errmax

!declare coefficient and boundary condition input subroutines external
external cof,bndc

!set input integer arguments
intl = 0

!set boundary condition flags
nxa = 2
nxb = 1
nyc = 1
nyd = 2

!set grid sizes from parameter statements
ixp = iixp
jyq = jjyq
iex = iiex
jey = jjey
nx = nnx
ny = nny

!set multigrid arguments (w(2,1) cycling with fully weighted
!residual restriction and cubic prolongation)
mgopt(1) = 2
mgopt(2) = 2
mgopt(3) = 1
mgopt(4) = 3

!set two mg cycles for second-order approximation
maxcy = 2

!set no initial guess forcing full multigrid cycling
!on initial call
iguess = 0

!set work space length approximation from parameter statement
nwork = llwork

!set line-y relaxation
method = 2

!set end points of solution rectangle in (x,y) space
xa = 0.0
xb = 1.0
yc = 0.0
yd = 1.0

!set mesh increments
dlx = (xb-xa)/float(nx-1)
dly = (yd-yc)/float(ny-1)

!set for no error control flag
tolmax = 0.0

!set right hand side in rhs
!initialize phi to zero
do i=1,nx
   x = xa+float(i-1)*dlx
   do j=1,ny
      y = yc+float(j-1)*dly
      call cof(x,y,cxx,cyy,cx,cy,ce)
      call exact(x,y,pxx,pyy,px,py,pe)
      rhs(i,j) = cxx*pxx+cyy*pyy+cx*px+cy*py+ce*pe
      phi(i,j) = 0.0
   end do
end do

!set specified boundaries in phi
x = xb
do j=1,ny
   y = yc+float(j-1)*dly
   call exact(x,y,pxx,pyy,px,py,pe)
   phi(nx,j) = pe
end do
y = yc
do i=1,nx
   x = xa+float(i-1)*dlx
   call exact(x,y,pxx,pyy,px,py,pe)
   phi(i,1) = pe
end do
write(*,100)
write(*,101) (iprm(i),i=1,15)
write(*,102) (mgopt(i),i=1,4)
write(*,103) xa,xb,yc,yd,tolmax
write(*,104) intl
write(*,200) ierror,iprm(16)
write(*,106) intl,method,iguess

call mud2(iprm,fprm,work,cof,bndc,rhs,phi,mgopt,ierror)

if (ierror.gt.0) call exit(0)
intl = 1
call mud2(iprm,fprm,work,cof,bndc,rhs,phi,mgopt,ierror)
write (*,107) ierror
if (ierror.gt.0) call exit(0)

!attempt to improve approximation to fourth order
call mud24(work,phi,ierror)
write (*,108) ierror
if (ierror.gt.0) call exit(0)

!compute and print maximum norm of error
errmax = 0.0
do j=1,ny
   y = yc+(j-1)*dly
   do i=1,nx
      x = xa+(i-1)*dlx
      call exact(x,y,pxx,pyy,px,py,pe)
      errmax = dmax1(errmax,abs((phi(i,j)-pe)))
   end do
end do
write(*,201) errmax
print*,"PASSED"

100 format(//' mud2 test ')
101 format(/' integer input arguments ', &
    &/' intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2, &
    &' nyd = ',i2, &
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
104 format(/' discretization call to mud2', ' intl = ', i2)
106 format(/' approximation call to mud2', &
    &/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
107 format(' ierror = ',i2)
108 format(/' mud24 test ', ' ierror = ',i2)

200 format(' ierror = ',i2, ' minimum work space = ',i7)
201 format(' maximum error  =  ',e10.3)

end program test_mudpack_2d

!> input pde coefficients at any grid point (x,y) in the solution region
subroutine cof(x,y,cxx,cyy,cx,cy,ce)
implicit none
real(8) :: x,y,cxx,cyy,cx,cy,ce
cxx = 1.+y*y
cyy = exp(-(x+y))
cx = 0.
cy = -cyy
ce = -(x+y)
return
end

!> input mixed derivative b.c. to mud2
subroutine bndc(kbdy,xory,alfa,gbdy)
implicit none
integer :: kbdy
real(8) :: xory,alfa,gbdy,x,y,pe,px,py,pxx,pyy
real(8) :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2/xa,xb,yc,yd,tolmax,relmax

if (kbdy.eq.1) then  ! x=xa boundary
   y = xory
   x = xa
   call exact(x,y,pxx,pyy,px,py,pe)
   alfa = -y
   gbdy = px + alfa*pe
   return
end if

if (kbdy.eq.4) then  ! y=yd boundary
   y = yd
   x = xory
   call exact(x,y,pxx,pyy,px,py,pe)
   alfa = x
   gbdy = py + alfa*pe
   return
end if

end

!> this subroutine is used to set an exact solution for testing mud2
subroutine exact(x,y,pxx,pyy,px,py,pe)
implicit none
real(8) :: x,y,pxx,pyy,px,py,pe

pe  = x**5+y**5+1.
px  = 5.*x**4
py  = 5.*y**4
pxx = 20.*x**3
pyy = 20.*y**3
return

end
