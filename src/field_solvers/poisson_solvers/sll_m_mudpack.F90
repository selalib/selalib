#ifndef DOXYGEN_SHOULD_SKIP_THIS

!> @ingroup poisson_solvers
!> @brief
!> Base module to provide interface to mudpack library
!> @details
!> This library contains multigrid solvers for PDE equations
!> You have to download and install mudpack 
!> http://www2.cisl.ucar.edu/resources/legacy/mudpack
module sll_m_mudpack
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

! use F77_mudpack, only: &
!   mud24cr, &
!   mud24sp, &
!   mud2cr, &
!   mud2sp

  implicit none

  public :: &
    sll_s_delete_mudpack_cartesian, &
    sll_s_initialize_mudpack_cartesian, &
    sll_s_initialize_mudpack_polar, &
    sll_o_create, &
    sll_o_delete, &
    sll_t_mudpack_solver, &
    sll_o_solve, &
    sll_s_solve_mudpack_cartesian, &
    sll_s_solve_mudpack_polar

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Mudpack solver cartesian 2d
type :: sll_t_mudpack_solver

   sll_real64, dimension(:), allocatable :: work !< array for tmp data
   sll_int32  :: mgopt(4)           !< Option to control multigrid
   sll_int32  :: iprm(16)           !< Indices to control grid sizes
   sll_real64 :: fprm(6)            !< Real to set boundary conditions
   sll_int32  :: iguess             !< Initial solution or loop over time
   sll_int32, pointer :: iwork(:,:) !< Internal work array for mudpack library

end type sll_t_mudpack_solver

integer, parameter :: CARTESIAN_2D = 2    !< geometry parameter
integer, parameter :: CARTESIAN_3D = 3    !< geometry parameter
integer, parameter :: POLAR        = 11   !< geometry parameter
integer, parameter :: CYLINDRICAL  = 12   !< geometry parameter

interface sll_o_create
  module procedure sll_s_initialize_mudpack_cartesian
end interface sll_o_create

interface sll_o_solve
  module procedure sll_s_solve_mudpack_cartesian
end interface sll_o_solve

interface sll_o_delete
  module procedure sll_s_delete_mudpack_cartesian
end interface sll_o_delete

contains

!> Initialize the Poisson solver using mudpack library
subroutine sll_s_initialize_mudpack_cartesian( self,                        &
                                        eta1_min, eta1_max, nc_eta1, &
                                        eta2_min, eta2_max, nc_eta2, &
                                        bc_eta1_left, bc_eta1_right, &
                                        bc_eta2_left, bc_eta2_right )
implicit none

! set grid size params
type(sll_t_mudpack_solver):: self          !< Data structure for solver
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
sll_int32,  parameter   :: iixp = 4 , jjyq = 4
sll_int32               :: iiex, jjey, llwork

sll_real64 :: phi(nc_eta1+1,nc_eta2+1) !< electric potential
sll_real64 :: rhs(nc_eta1+1,nc_eta2+1) !< charge density

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
!sll_real64,dimension(:),allocatable :: cxx_array

#ifdef DEBUG
sll_int32 :: i
#endif

equivalence(intl,iprm)
equivalence(xa,fprm)

nx = nc_eta1+1
ny = nc_eta2+1

! set minimum required work space
llwork=(7*(nx+2)*(ny+2)+44*nx*ny)/3

allocate(self%work(llwork))
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
self%mgopt(1) = 2
self%mgopt(2) = 2
self%mgopt(3) = 1
self%mgopt(4) = 3

!set for three cycles to ensure second-order approximation is computed
maxcy = 5

!set no initial guess forcing full multigrid cycling
self%iguess = 0
iguess = self%iguess

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
tolmax = 0.0_f64

#ifdef DEBUG
write(*,101) (iprm(i),i=1,15)
write(*,102) (self%mgopt(i),i=1,4)
write(*,103) xa,xb,yc,yd,tolmax
write(*,104) intl
#endif

call mud2sp(iprm,fprm,self%work,cofx,cofy,bndsp,rhs,phi,self%mgopt,error)

#ifdef DEBUG
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
     
#endif
end subroutine sll_s_initialize_mudpack_cartesian


!> Compute the potential using mudpack library on cartesian mesh
subroutine sll_s_solve_mudpack_cartesian(self, phi, rhs, ex, ey, nrj)

type(sll_t_mudpack_solver)     :: self      !< Data structure for Poisson solver
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
sll_real64 :: dx, dy

sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
common/itmud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

equivalence(intl,iprm)
equivalence(xa,fprm)

!set initial guess because solve should be called every time step in a
!time dependent problem and the elliptic operator does not depend on time.
if (self%iguess == 0) then
   iguess = 0
else
   self%iguess = 1
    iguess = 1
endif

!attempt solution
intl = 1
#ifdef DEBUG
write(*,106) intl,method,iguess
#endif

if ( nxa == 0 .and. nyc == 0 ) &
   rhs = rhs - sum(rhs) / real(nx*ny,f64)

call mud2sp(iprm,fprm,self%work,cofx,cofy,bndsp,rhs,phi,self%mgopt,error)

#ifdef DEBUG
write(*,107) error
if (error > 0) call exit(0)
#endif

if ( nxa == 0 .and. nyc == 0 ) &
   phi = phi - sum(phi) / real(nx*ny,f64)

iguess = 1
! attempt to improve approximation to fourth order
call mud24sp(self%work,phi,error)

#ifdef DEBUG
write (*,108) error
if (error > 0) call exit(0)

106 format(/'#approximation call to mud2sp', &
    &/'# intl = ',i2, ' method = ',i2,' iguess = ',i2)
107 format('#error = ',i2)
108 format(/'# mud24sp test ', ' error = ',i2)

#endif

if (present(ex) .and. present(ey)) then
   dx = (xb-xa)/(nx-1)
   dy = (yd-yc)/(ny-1)
   do i = 2, nx-1
      ex(i,:) = (phi(i+1,:)-phi(i-1,:)) / (2*dx)
   end do
   do j = 2, ny-1
      ey(:,j) = (phi(:,j+1)-phi(:,j-1)) / (2*dy)
   end do

   if (nxa == 0 ) then
      ex(nx,:) = (phi(2,:)-phi(nx-1,:))/(2*dx)
      ex(1,:)  = ex(nx,:)
   end if
   if (nyc == 0 ) then
      ey(:,ny) = (phi(:,2)-phi(:,ny-1))/(2*dy)
      ey(:,1) = ey(:,ny)
   end if

   if (present(nrj)) then 
      nrj=sum(ex*ex+ey*ey)*dx*dy
   end if

end if

end subroutine sll_s_solve_mudpack_cartesian

!> deallocate the mudpack solver
subroutine sll_s_delete_mudpack_cartesian(self)
type(sll_t_mudpack_solver)        :: self          !< Data structure for solver

   deallocate(self%work)

end subroutine sll_s_delete_mudpack_cartesian

!> Initialize the Poisson solver in polar coordinates using MUDPACK
!> library
subroutine sll_s_initialize_mudpack_polar(self,                      &
                                    r_min, r_max, nr,          &
                                    theta_min, theta_max, nth, &
                                    bc_r_min, bc_r_max,        &
                                    bc_theta_min, bc_theta_max )
implicit none

type(sll_t_mudpack_solver)       :: self      !< Solver object
sll_real64, intent(in) :: r_min     !< radius min
sll_real64, intent(in) :: r_max     !< radius min
sll_real64, intent(in) :: theta_min !< theta min
sll_real64, intent(in) :: theta_max !< theta max
sll_int32, intent(in)  :: nr        !< radius number of points
sll_int32, intent(in)  :: nth       !< theta number of points
sll_int32 :: icall
!sll_int32 :: iiex,jjey
sll_int32 :: llwork
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
sll_int32  :: i
!sll_int32 :: j
sll_int32 :: ierror
sll_int32  :: iprm(16)
sll_real64 :: fprm(6)

equivalence(intl,iprm)
equivalence(xa,fprm)

nx = nr
ny = nth

! set minimum required work space
llwork=(7*(nx+2)*(ny+2)+44*nx*ny)/3
      
allocate(self%work(llwork))
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
self%mgopt(1) = 2
self%mgopt(2) = 2
self%mgopt(3) = 1
self%mgopt(4) = 3

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
tolmax = 0.0_f64

write(*,100)
write(*,101) (iprm(i),i=1,15)
write(*,102) (self%mgopt(i),i=1,4)
write(*,103) xa,xb,yc,yd,tolmax
write(*,104) intl

call mud2cr(iprm,fprm,self%work,coef_polar,bndcr,rhs,phi,self%mgopt,ierror)

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
end subroutine sll_s_initialize_mudpack_polar


!> Solve the Poisson equation and get the potential
subroutine sll_s_solve_mudpack_polar(self, phi, rhs)
implicit none

! set grid size params
type(sll_t_mudpack_solver) :: self  !< solver data object
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

icall = 1
intl  = 1
write(*,106) intl,method,iguess
! attempt solution
call mud2cr(iprm,fprm,self%work,coef_polar,bndcr,rhs,phi,self%mgopt,ierror)
SLL_ASSERT(ierror == 0)
! attempt fourth order approximation
call mud24cr(self%work,coef_polar,bndcr,phi,ierror)
SLL_ASSERT(ierror == 0)

106 format(/' approximation call to mud2cr', &
    /' intl = ',i2, ' method = ',i2,' iguess = ',i2)

end subroutine sll_s_solve_mudpack_polar

!> input pde coefficients at any grid point (x,y) in the solution region
!> (xa.le.x.le.xb,yc.le.y.le.yd) to mud2cr
subroutine coef_polar(x,y,cxx,cxy,cyy,cx,cy,ce)
real(8) :: x,y,cxx,cxy,cyy,cx,cy,ce
cxx = 1.0_8 +0.0_8*y
cxy = 0.0_8 
cyy = 1.0_8 / (x*x) 
cx  = 1.0_8 / x 
cy  = 0.0_8 
ce  = 0.0_8 
end subroutine coef_polar

!> input mixed "oblique" derivative b.c. to mud2cr
!> at upper y boundary
subroutine bndcr(kbdy,xory,alfa,beta,gama,gbdy)
integer  :: kbdy
real(8)  :: xory,alfa,beta,gama,gbdy

if (kbdy.eq.2) then

   ! x=xb boundary.
   ! b.c. has the form alfyd(x)*px+betyd(x)*py+gamyd(x)*pe = gbdyd(x)
   ! where x = yorx.   alfa,beta,gama,gbdy corresponding to alfyd(x),
   ! betyd(x),gamyd(x),gbdyd(y) must be output.

   alfa = 1.0_8+0.0_8*xory
   beta = 0.0_8
   gama = 0.0_8
   gbdy = 0.0_8

end if

end subroutine bndcr

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
real(8)  :: x,cxx,cx,cex
cxx = 1.0_8 +0.0_8*x 
cx  = 0.0_8
cex = 0.0_8
end subroutine cofx

!> input y dependent coefficients
subroutine cofy(y,cyy,cy,cey)
real(8)  :: y,cyy,cy,cey
cyy = 1.0_8 +0.0_8*y
cy  = 0.0_8
cey = 0.0_8
end subroutine cofy

!> input mixed derivative b.c. to mud2sp
subroutine bndsp(kbdy,xory,alfa,gbdy)
integer  :: kbdy
real(8)  :: xory,alfa,gbdy,x,y,pe,px,py
real(8)  :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax

!subroutine not used in periodic case
if (kbdy == 1) then  ! x=xa boundary
   y = xory
   x = xa
   alfa = -1.0_8
   gbdy = px + alfa*pe
   return
end if

if (kbdy == 4) then  ! y=yd boundary
   y = yd
   x = xory
   alfa = 1.0_8
   gbdy = py + alfa*pe
   return
end if
end subroutine bndsp

end module sll_m_mudpack

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
