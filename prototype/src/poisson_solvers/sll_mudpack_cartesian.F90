module sll_mudpack_cartesian
#include "sll_working_precision.h"
#include "sll_utilities.h"
implicit none


contains

subroutine solve_mudpack_cartesian()
!set grid sizes with parameter statements
sll_int32, parameter :: iixp = 2 , jjyq = 3 , iiex = 6, jjey = 5 
sll_int32, parameter :: nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1 

!set minimal required work space (see tmud2sp.f)
sll_int32, parameter :: llwork = 13264 
sll_real64 :: phi(nnx,nny),rhs(nnx,nny),work(llwork)

!put integer and floating point argument names in contiguous
!storeage for labelling in vectors iprm,fprm
sll_int32 :: iprm(16),mgopt(4)
sll_real64 :: fprm(6)
sll_int32 :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
           iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
sll_int32 :: i,j,ierror
sll_real64 :: dlx,dly,x,y,cxx,cyy,cx,cy,ce,pxx,pyy,px,py,pe,errmax
sll_real64 :: cex,cey

common/itmud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
equivalence(intl,iprm)
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax
equivalence(xa,fprm)

!declare coefficient and boundary condition input subroutines external
external cofx,cofy,bndsp

!set input integer arguments
intl = 0

!set boundary condition flags
nxa  = 2
nxb  = 1
nyc  = 1
nyd  = 2

!set grid sizes from parameter statements
ixp  = iixp
jyq  = jjyq
iex  = iiex
jey  = jjey
nx   = nnx
ny   = nny

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
xa = 0.0
xb = 1.0
yc = 0.0
yd = 1.0

!set mesh increments
dlx = (xb-xa)/float(nx-1)
dly = (yd-yc)/float(ny-1)

!set for no error control flag
tolmax = 0.0

!set right hand side in rhs and initialize phi to zero
do i=1,nx
   x = xa+float(i-1)*dlx
   call cofx(x,cxx,cx,cex)
   do j=1,ny
      y = yc+float(j-1)*dly
      call cofy(y,cyy,cy,cey)
      ce = cex+cey
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

!intiialization call
write(*,104) intl
call mud2sp(iprm,fprm,work,cofx,cofy,bndsp,rhs,phi,mgopt,ierror)

!print error parameter and minimum work space requirement
write (*,200) ierror,iprm(16)
if (ierror.gt.0) call exit(0)

!attempt solution
intl = 1
write(*,106) intl,method,iguess

call mud2sp(iprm,fprm,work,cofx,cofy,bndsp,rhs,phi,mgopt,ierror)
write(*,107) ierror
if (ierror.gt.0) call exit(0)

if (ierror .le. 0) then

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

end if

! attempt to improve approximation to fourth order

call mud24sp(work,phi,ierror)

call plot_field(rhs, "rhs", 0)
call plot_field(phi, "phi", 0)

write (*,108) ierror
if (ierror.gt.0) call exit(0)

if (ierror .le. 0) then

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

end if

print*,"PASSED"

100 format(//' test_mudpack_2d ')
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
    write(*,103) xa,xb,yc,yd,tolmax
103 format(/' floating point input parameters ', &
    &/' xa = ',f6.3,' xb = ',f6.3,' yc = ',f6.3,' yd = ',f6.3, &
    &/' tolerance (error control) =   ',e10.3)
104 format(/' discretization call to mud2sp', ' intl = ', i2)
106 format(/' approximation call to mud2sp', &
    &/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
107 format(' ierror = ',i2)
108 format(/' mud24sp test ', ' ierror = ',i2)
200 format(' ierror = ',i2, ' minimum work space = ',i7)
201 format(' maximum error  =  ',e10.3)
     
end subroutine solve_mudpack_cartesian


subroutine plot_field(f, fname, iplot)

   sll_int32 :: iplot, i, j
   sll_real64 :: x, y
   sll_real64, dimension(:,:) :: f
   character(len=*) :: fname
   character(len=4) :: cplot
 
   call int2string(iplot,cplot)

   open(11, file=fname//cplot//".dat")
   do i = 1, size(f,1)
      do j = 1, size(f,2)
         x = xa + (i-1)*dlx
         y = yc + (j-1)*dly
         write(11,*) x,y,f(i,j)
      end do
      write(11,*)
   end do
   close(11)
   
   open( 90, file = fname//'.gnu', position="append" )
   if ( iplot == 1 ) then
      rewind(90)
      !write(90,*)"set cbrange[-1:1]"
      !write(90,*)"set pm3d"
      write(90,*)"set surf"
      write(90,*)"set term x11"
   end if

   write(90,*)"set title 'step = ",iplot,"'"
   write(90,"(a)")"splot '"//fname//cplot//".dat' u 1:2:3 w lines"
   close(90)

end subroutine plot_field


subroutine cofx(x,cxx,cx,cex)
!
!     input x dependent coefficients
!
implicit none
real(8)  :: x,cxx,cx,cex
cxx = 1.0+x*x
cx = 0.0
cex = -x
return
end

subroutine cofy(y,cyy,cy,cey)
!
!     input y dependent coefficients
!
implicit none
real(8)  :: y,cyy,cy,cey
cyy = exp(1.0-y)
cy = -cyy
cey = -y
return
end

subroutine bndsp(kbdy,xory,alfa,gbdy)
!
!     input mixed derivative b.c. to mud2sp
!
implicit none
integer  :: kbdy
real(8)  :: xory,alfa,gbdy,x,y,pe,px,py,pxx,pyy
real(8)  :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax
if (kbdy.eq.1) then  ! x=xa boundary
y = xory
x = xa
call exact(x,y,pxx,pyy,px,py,pe)
alfa = -1.0
gbdy = px + alfa*pe
return
end if
if (kbdy.eq.4) then  ! y=yd boundary
y = yd
x = xory
call exact(x,y,pxx,pyy,px,py,pe)
alfa = 1.0
gbdy = py + alfa*pe
return
end if
end

subroutine exact(x,y,pxx,pyy,px,py,pe)
!
!     set an exact solution for testing mud2sp
!
implicit none
real(8)  :: x,y,pxx,pyy,px,py,pe
pe = (x**3+y**3+1.0)/3.0
px = x*x
py = y*y
pxx = x+x
pyy = y+y
return
end

end module sll_mudpack_cartesian

