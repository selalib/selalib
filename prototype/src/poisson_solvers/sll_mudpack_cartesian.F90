!> @brief
!> Module to use MUDPACK library to solve Poisson equation 
!> @details We consider a regular cartesian mesh
!>
!> MUDPACK library is available on this website
!> http://www2.cisl.ucar.edu/resources/legacy/mudpack
module sll_mudpack_cartesian
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_file_io.h"

use sll_mudpack_base
!use sll_cubic_spline_interpolator_1d

implicit none

interface initialize
   module procedure  initialize_mudpack_cartesian
end interface initialize

interface solve
   module procedure  solve_mudpack_cartesian
end interface solve

!class(sll_interpolator_1d_base), pointer   :: cxx_interp
contains

!> Initialize the Poisson solver using mudpack library
subroutine initialize_mudpack_cartesian( this,                        &
                                        eta1_min, eta1_max, nc_eta1, &
                                        eta2_min, eta2_max, nc_eta2, &
                                        bc_eta1_left, bc_eta1_right, &
                                        bc_eta2_left, bc_eta2_right )
implicit none

! set grid size params
type(mudpack_2d)        :: this          !< Data structure for solver
sll_real64, intent(in)  :: eta1_min      !< left corner x direction
sll_real64, intent(in)  :: eta1_max      !< right corner x direction
sll_real64, intent(in)  :: eta2_min      !< left corner x direction
sll_real64, intent(in)  :: eta2_max      !< right corner x direction
sll_int32,  intent(in)  :: nc_eta1       !< number of cells x direction
sll_int32,  intent(in)  :: nc_eta2       !< number of cells y direction
sll_int32,  optional    :: bc_eta1_left  !< boundary condtion
sll_int32,  optional    :: bc_eta1_right !< boundary condtion
sll_int32,  optional    :: bc_eta2_left  !< boundary condtion
sll_int32,  optional    :: bc_eta2_right !< boundary condtion
sll_int32,  parameter   :: iixp = 2 , jjyq = 2
sll_int32               :: iiex, jjey, llwork

sll_real64 :: phi(nc_eta1+1,nc_eta2+1) !< electric potential
sll_real64 :: rhs(nc_eta1+1,nc_eta2+1) !< charge density

!put integer and floating point argument names in contiguous
!storeage for labelling in vectors iprm,fprm
sll_int32  :: iprm(16)
sll_real64 :: fprm(6)
sll_int32  :: i,error
sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
common/itmud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax
!sll_real64,dimension(:),allocatable :: cxx_array

equivalence(intl,iprm)
equivalence(xa,fprm)

!declare coefficient and boundary condition input subroutines external
external cofx,cofy,bndsp

nx = nc_eta1+1
ny = nc_eta2+1

! set minimum required work space
llwork=(7*(nx+2)*(ny+2)+44*nx*ny)/3

!cxx_interp => new_cubic_spline_1d_interpolator( &
!          nx, &
!          eta1_min, &
!          eta1_max, &
!          SLL_PERIODIC)
!allocate(cxx_array(nx))          
!
!cxx_array = 1._f64          
!          
!call cxx_interp%compute_interpolants( cxx_array )          

      
allocate(this%work(llwork))
iiex = ceiling(log((nx-1.)/iixp)/log(2.))+1
jjey = ceiling(log((ny-1.)/jjyq)/log(2.))+1

!set input integer arguments
intl = 0

!set boundary condition flags
nxa  = merge(bc_eta1_left,  0, present(bc_eta1_left))
nxb  = merge(bc_eta1_right, 0, present(bc_eta1_right))
nyc  = merge(bc_eta2_left,  0, present(bc_eta2_left))
nyd  = merge(bc_eta2_right, 0, present(bc_eta2_right))

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
this%mgopt(1) = 2
this%mgopt(2) = 2
this%mgopt(3) = 1
this%mgopt(4) = 3

!set for three cycles to ensure second-order approximation is computed
maxcy = 5

!set no initial guess forcing full multigrid cycling
this%iguess = 0
iguess = this%iguess

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

#ifdef DEBUG
write(*,101) (iprm(i),i=1,15)
write(*,102) (this%mgopt(i),i=1,4)
write(*,103) xa,xb,yc,yd,tolmax
write(*,104) intl
#endif

call mud2sp(iprm,fprm,this%work,cofx,cofy,bndsp,rhs,phi,this%mgopt,error)

!print error parameter and minimum work space requirement
write (*,200) error,iprm(16)
if (error > 0) call exit(0)

101 format(/'# integer input arguments ', &
    &/'#intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2, &
    &/'# ixp = ',i2,' jyq = ',i2,' iex = ',i2,' jey = ',i2 &
    &/'# nx = ',i3,' ny = ',i3,' iguess = ',i2,' maxcy = ',i2, &
    &/'# method = ',i2, ' work space estimate = ',i7)
102 format(/'# multigrid option arguments ', &
    &/'# kcycle = ',i2, &
    &/'# iprer = ',i2, &
    &/'# ipost = ',i2, &
    &/'# intpol = ',i2)
103 format(/'# floating point input parameters ', &
    &/'# xa = ',f6.3,' xb = ',f6.3,' yc = ',f6.3,' yd = ',f6.3, &
    &/'# tolerance (error control) =   ',e10.3)
104 format(/'# discretization call to mud2sp', ' intl = ', i2)
200 format('# error = ',i2, ' minimum work space = ',i7)
     
end subroutine initialize_mudpack_cartesian


!> Compute the potential using mudpack library on cartesian mesh
subroutine solve_mudpack_cartesian(this, phi, rhs, ex, ey, nrj)

type(mudpack_2d)     :: this      !< Data structure for Poisson solver
sll_real64           :: phi(:,:)  !< Electric potential
sll_real64           :: rhs(:,:)  !< Charge density
sll_real64, optional :: ex(:,:)   !< Electric field
sll_real64, optional :: ey(:,:)   !< Electric field
sll_real64, optional :: nrj

!put integer and floating point argument names in contiguous
!storeage for labelling in vectors iprm,fprm
sll_int32  :: iprm(16)
sll_real64 :: fprm(6)
sll_int32  :: error
sll_int32  :: i, j
sll_real64 :: dx, dy, avg

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

!set initial guess because solve should be called every time step in a
!time dependent problem and the elliptic operator does not depend on time.
iguess = this%iguess

!attempt solution
intl = 1
#ifdef DEBUG
write(*,106) intl,method,iguess
#endif

rhs = rhs - sum(rhs) / (nx*ny)
call mud2sp(iprm,fprm,this%work,cofx,cofy,bndsp,rhs,phi,this%mgopt,error)

#ifdef DEBUG
write(*,107) error
if (error > 0) call exit(0)
#endif

iguess = 1
! attempt to improve approximation to fourth order
call mud24sp(this%work,phi,error)

#ifdef DEBUG
write (*,108) error
if (error > 0) call exit(0)
#endif

106 format(/'#approximation call to mud2sp', &
    &/'# intl = ',i2, ' method = ',i2,' iguess = ',i2)
107 format('#error = ',i2)
108 format(/'# mud24sp test ', ' error = ',i2)

if (present(ex) .and. present(ey)) then
   dx = (xb-xa)/(nx-1)
   dy = (yc-yd)/(ny-1)
   do j = 1, nx-1
      do i = 1, ny-1
         ex(i,j) = (phi(i+1,j)-phi(i,j)) / dx
         ey(i,j) = (phi(i,j+1)-phi(i,j)) / dy
      end do
   end do

   if (present(nrj)) then 
      nrj=sum(ex(1:nx-1,1:ny-1)*ex(1:nx-1,1:ny-1)+ &
              ey(1:nx-1,1:ny-1)*ey(1:nx-1,1:ny-1))*dx*dy
   end if

end if

     
end subroutine solve_mudpack_cartesian

end module sll_mudpack_cartesian

!the form of the pde solved is:
!
!
!          cxx(x)*pxx + cx(x)*px + cex(x)*p(x,y) +
!
!          cyy(y)*pyy + cy(y)*py + cey(y)*p(x,y) = r(x,y)
!
!     pxx,pyy,px,py are second and first partial derivatives of the
!     unknown real solution function p(x,y) with respect to the
!     independent variables x,y.  cxx,cx,cex,cyy,cy,cey are the known
!     real coefficients of the elliptic pde and r(x,y) is the known
!     real right hand side of the equation.  cxx and cyy should be
!     positive for all x,y in the solution region.  If some of the
!     coefficients depend on both x and y then the PDE is nonseparable.

!> input x dependent coefficients
subroutine cofx(x,cxx,cx,cex)
use sll_mudpack_cartesian
implicit none
real(8)  :: x,cxx,cx,cex
cxx = 1.0 +0.0*x !cxx_interp%interpolate_value(x)
cx  = 0.0
cex = 0.0
return
end

!> input y dependent coefficients
subroutine cofy(y,cyy,cy,cey)
implicit none
real(8)  :: y,cyy,cy,cey
cyy = 1.0 +0.0*y
cy  = 0.0
cey = 0.0
return
end

!> input mixed derivative b.c. to mud2sp
subroutine bndsp(kbdy,xory,alfa,gbdy)
implicit none
integer  :: kbdy
real(8)  :: xory,alfa,gbdy,x,y,pe,px,py
real(8)  :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

!subroutine not used in periodic case
if (kbdy == 1) then  ! x=xa boundary
   y = xory
   x = xa
   alfa = -1.0
   gbdy = px + alfa*pe
   return
end if

if (kbdy == 4) then  ! y=yd boundary
   y = yd
   x = xory
   alfa = 1.0
   gbdy = py + alfa*pe
   return
end if
end

