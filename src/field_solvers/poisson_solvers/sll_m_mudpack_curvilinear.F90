#ifndef DOXYGEN_SHOULD_SKIP_THIS
!> @ingroup poisson_solvers
!> @brief
!> Poisson solver in general coordinates using mudpack library
!> @details
!> red/black gauss-seidel point relaxation is used along with the
!> the default multigrid options.  
!>  - first mud2cr is called to generate a second-order approximation.  
!>  - then mud24cr is called to improve the estimate to fourth-order.
module sll_m_mudpack_curvilinear
#include "sll_working_precision.h"
#include "sll_assert.h"
use sll_m_common_coordinate_transformations
use sll_m_coordinate_transformation_2d_base
use sll_m_interpolators_2d_base
use sll_m_cubic_spline_interpolator_2d

implicit none
private

!> Mudpack solver cartesian 2d
type, public :: mudpack_2d

   sll_real64, dimension(:), allocatable :: work !< array for tmp data
   sll_int32  :: mgopt(4)           !< Option to control multigrid
   sll_int32  :: iprm(16)           !< Indices to control grid sizes
   sll_real64 :: fprm(6)            !< Real to set boundary conditions
   sll_int32  :: iguess             !< Initial solution or loop over time
   sll_int32, pointer :: iwork(:,:) !< Internal work array for mudpack library

end type mudpack_2d

integer, parameter, public :: SLL_SEPARABLE  = 1                        !< type of equation
integer, parameter, public :: SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS = 2 !< type of equation
integer, parameter, public :: SLL_NON_SEPARABLE_WITH_CROSS_TERMS = 3    !< type of equation

!> Interpolator to compute derivative xx
class(sll_c_interpolator_2d), pointer :: cxx_interp
!> Interpolator to compute derivative yy
class(sll_c_interpolator_2d), pointer :: cyy_interp
!> Interpolator to compute derivative xy
class(sll_c_interpolator_2d), pointer :: cxy_interp
!> Interpolator to compute derivative x
class(sll_c_interpolator_2d), pointer :: cx_interp
!> Interpolator to compute derivative y
class(sll_c_interpolator_2d), pointer :: cy_interp
!> Interpolator to compute rhs coefficient
class(sll_c_interpolator_2d), pointer :: ce_interp
!> PLEASE ADD DOCUMENTATION
class(sll_c_interpolator_2d), pointer :: a12_interp
!> PLEASE ADD DOCUMENTATION
class(sll_c_interpolator_2d), pointer :: a21_interp

!> Coordinate transformation of the mesh
class(sll_coordinate_transformation_2d_base), pointer :: transformation

interface sll_create
  module procedure initialize_poisson_curvilinear_mudpack
end interface sll_create

public :: sll_create
public :: cxx_interp
public :: cyy_interp
public :: cxy_interp
public :: cx_interp
public :: cy_interp
public :: ce_interp
public :: a12_interp
public :: a21_interp
public :: solve_poisson_curvilinear_mudpack

contains

!> Initialize the Poisson solver in curvilinear coordinates using MUDPACK
!> library
subroutine initialize_poisson_curvilinear_mudpack( &
   this,          &
   transf,        &
   b11,           &
   b12,           &
   b21,           &
   b22,           &
   c,             &
   eta1_min,      &
   eta1_max,      &
   nc_eta1,       &
   eta2_min,      &
   eta2_max,      &
   nc_eta2,       &
   bc_eta1_left,  &
   bc_eta1_right, &
   bc_eta2_left,  &
   bc_eta2_right)

type(mudpack_2d) :: this              !< Solver object
sll_real64, intent(in) :: eta1_min    !< eta1 min
sll_real64, intent(in) :: eta1_max    !< eta1 min
sll_real64, intent(in) :: eta2_min    !< eta2 min
sll_real64, intent(in) :: eta2_max    !< eta2 max
sll_int32, intent(in)  :: nc_eta1     !<  number of cells
sll_int32, intent(in)  :: nc_eta2     !<  number of cells
sll_int32 :: icall
!sll_int32 :: iiex,jjey
sll_int32 :: llwork
sll_int32 :: bc_eta1_left             !< left boundary condition r
sll_int32 :: bc_eta1_right            !< right boundary condition r
sll_int32 :: bc_eta2_left             !< left boundary condition theta
sll_int32 :: bc_eta2_right            !< right boundary condition theta

sll_real64, pointer ::  phi(:,:)      !< electric potential
sll_real64, pointer ::  rhs(:,:)      !< charge density

sll_real64, dimension(:,:), pointer :: b11 !< for general coordinate solver
sll_real64, dimension(:,:), pointer :: b12 !< for general coordinate solver
sll_real64, dimension(:,:), pointer :: b21 !< for general coordinate solver
sll_real64, dimension(:,:), pointer :: b22 !< for general coordinate solver
sll_real64, dimension(:,:), pointer :: c   !< for general coordinate solver

class(sll_coordinate_transformation_2d_base), pointer :: transf !< coordinate transformation

! put integer and floating point argument names in contiguous
! storeage for labelling in vectors iprm,fprm
sll_int32 :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,&
              iguess,maxcy,method,nwork,lwrkqd,itero
common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
              iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax
sll_int32  :: i,ierror
sll_int32  :: iprm(16)
sll_real64 :: fprm(6)
sll_real64,dimension(:,:),allocatable :: cxx_array
sll_real64,dimension(:,:),allocatable :: cyy_array
sll_real64,dimension(:,:),allocatable :: cxy_array
sll_real64,dimension(:,:),allocatable :: cx_array
sll_real64,dimension(:,:),allocatable :: cy_array
sll_real64,dimension(:,:),allocatable :: ce_array
sll_real64,dimension(:,:),allocatable :: a12_array
sll_real64,dimension(:,:),allocatable :: a21_array
sll_real64 :: delta1,delta2
sll_int32,  parameter   :: iixp = 2 , jjyq = 2

equivalence(intl,iprm)
equivalence(xa,fprm)

nx = nc_eta1+1
ny = nc_eta2+1

delta1   = (eta1_max - eta1_min)/real(nc_eta1,f64)
delta2   = (eta2_max - eta2_min)/real(nc_eta2,f64)
! set minimum required work space
llwork=(7*(nx+2)*(ny+2)+44*nx*ny)/3

allocate(cxx_array(nx,ny)) 
allocate(cyy_array(nx,ny)) 
allocate(cxy_array(nx,ny)) 
allocate(cx_array(nx,ny)) 
allocate(cy_array(nx,ny)) 
allocate(ce_array(nx,ny)) 
allocate(a12_array(nx,ny))
allocate(a21_array(nx,ny))
allocate(phi(nx,ny))

transformation => transf
cxx_interp => new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)
          
cyy_interp => new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC) 
          
 cxy_interp => new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)  
          
 cx_interp => new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC) 
 cy_interp => new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)    
                                         
ce_interp => new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)   
a12_interp => new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC) 
a21_interp => new_cubic_spline_interpolator_2d( &
          nx, &
          ny, &
          eta1_min, &
          eta1_max, &
          eta2_min, &
          eta2_max, &
          SLL_PERIODIC, &
          SLL_PERIODIC)                             
!cxx_array = 1._f64          
call coefxxyy_array(b11,b12,b21,b22,transf,eta1_min,eta2_min,delta1,delta2,nx,ny,cxx_array,cyy_array)          
call cxx_interp%compute_interpolants( cxx_array )  
call cyy_interp%compute_interpolants( cyy_array ) 

call coefxy_array(b11,b12,b21,b22,transf,eta1_min,eta2_min,delta1,delta2,nx,ny,cxy_array)
call cxy_interp%compute_interpolants( cxy_array ) 

call a12_a21_array(b11,b12,b21,b22,transf,eta1_min,eta2_min,delta1,delta2,nx,ny,a12_array,a21_array)
call a12_interp%compute_interpolants( a12_array ) 
call a21_interp%compute_interpolants( a21_array ) 

call coefx_array(eta1_min,eta2_min,delta1,delta2,nx,ny,cx_array)
call cx_interp%compute_interpolants( cx_array ) 

call coefy_array(eta1_min,eta2_min,delta1,delta2,nx,ny,cy_array)
call cy_interp%compute_interpolants( cy_array ) 
ce_array = -c
call ce_interp%compute_interpolants( ce_array ) 
     
allocate(this%work(llwork))
icall = 0

! set input sll_int32 arguments
intl = 0

! set boundary condition flags
nxa = bc_eta1_left  
nxb = bc_eta1_right 
nyc = bc_eta2_left  
nyd = bc_eta2_right 
print*,nxa,nxb,nyc,nyd
! set grid sizes from parameter statements
ixp = iixp 
jyq = jjyq 
iex = ceiling(log((nx-1.)/ixp)/log(2.))+1
jey = ceiling(log((ny-1.)/jyq)/log(2.))+1

nx = ixp*(2**(iex-1))+1
ny = jyq*(2**(jey-1))+1
allocate(this%iwork(ixp+1,jyq+1))
if (nx /= nc_eta1+1 .or. ny /= nc_eta2+1) then
   print*, "nx,nc_eta1+1=", nx, nc_eta1+1
   print*, "ny,nc_eta2+1=", ny, nc_eta2+1
   stop ' nx or ny different in sll_m_mudpack_curvilinear '
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
iguess =  0 !1

! set work space length approximation from parameter statement
nwork = llwork

! set point relaxation
method = 0

! set mesh increments
xa = eta1_min
xb = eta1_max
yc = eta2_min
yd = eta2_max

! set for no error control flag
tolmax = 0.0_8

write(*,100)
write(*,101) (iprm(i),i=1,15)
write(*,102) (this%mgopt(i),i=1,4)
write(*,103) xa,xb,yc,yd,tolmax
write(*,104) intl

!call mud2cr(iprm,fprm,this%work,coefcr,bndcr,rhs,phi,this%mgopt,ierror)
 call muh2cr(iprm,fprm,this%work,this%iwork,coefcr,bndcr,rhs,phi,this%mgopt,ierror)
!call mud2sp(iprm,fprm,this%work,cofx,cofy,bndcr,rhs,phi,this%mgopt,ierror)
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

end subroutine initialize_poisson_curvilinear_mudpack


!> Solve the Poisson equation and get the potential
subroutine solve_poisson_curvilinear_mudpack(this, phi, rho)
! set grid size params
type(mudpack_2d) :: this  !< solver data object
sll_int32 :: icall
sll_int32, parameter :: iixp = 2 , jjyq = 2

sll_real64, intent(inout) ::  phi(:,:) !< electric potential
sll_real64, intent(inout) ::  rho(:,:) !< charge density
sll_real64, pointer       ::  rhs(:,:) !< charge density

! put sll_int32 and floating point argument names in contiguous
! storeage for labelling in vectors iprm,fprm

sll_int32  :: intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
sll_int32  :: iguess,maxcy,method,nwork,lwrkqd,itero
sll_real64 :: xa,xb,yc,yd,tolmax,relmax
sll_int32  :: ierror
sll_int32  :: iprm(16)
sll_real64 :: fprm(6),eta1,eta2
sll_int32 :: i1,i2

common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
                iguess,maxcy,method,nwork,lwrkqd,itero
common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax

equivalence(intl,iprm)
equivalence(xa,fprm)

allocate(rhs(nx,ny))
rhs=0._f64
    do i2=1,ny
      eta2=yc+real(i2-1,f64)*(yd-yc)/real(ny-1,8)
      do i1=1,nx
        eta1=xa+real(i1-1,f64)*(xb-xa)/real(nx-1,8)
        rhs(i1,i2)=-rho(i1,i2)*transformation%jacobian(eta1,eta2)
      end do
    end do
  if(nxa == SLL_DIRICHLET) then
       do i2=1,ny
          phi(1,i2) = 0._f64
       end do
    endif
    if(nxb == SLL_DIRICHLET) then
       do i2=1,ny
          phi(nx,i2) = 0._f64
       end do
    endif
    if(nyc == SLL_DIRICHLET) then
       do i1=1,nx
          phi(i1,1) = 0._f64
       end do
    endif
    if(nyd == SLL_DIRICHLET) then
       do i1=1,nx
          phi(i1,ny) = 0._f64
       end do
    endif 

icall = 1
intl  = 1
!YG write(*,106) intl,method,iguess

! attempt solution
!call mud2cr(iprm,fprm,this%work,coefcr,bndcr,rhs,phi,this%mgopt,ierror)
call muh2cr(iprm,fprm,this%work,this%iwork,coefcr,bndcr,rhs,phi,this%mgopt,ierror)
!call mud2sp(iprm,fprm,this%work,cofx,cofy,bndcr,rhs,phi,this%mgopt,ierror)
!SLL_ASSERT(ierror == 0)
!YG write(*,107) ierror
if (ierror > 0) call exit(0)

! attempt fourth order approximation
!call mud24cr(this%work,coefcr,bndcr,phi,ierror)
call muh24cr(this%work,this%iwork,coefcr,bndcr,phi,ierror)
!call mud24sp(this%work,phi,ierror)
!SLL_ASSERT(ierror == 0)
!YG write (*,108) ierror
if (ierror > 0) call exit(0)

!YG 106 format(/' approximation call to mud2sp', &
!YG     &/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
!YG 107 format(' error = ',i2)
!YG 108 format(/' mud24cr test ', ' error = ',i2)

deallocate(rhs)
return
end subroutine solve_poisson_curvilinear_mudpack

subroutine coefxxyy_array(b11,b12,b21,b22,transf,eta1_min,eta2_min, &
                         delta1,delta2,nx,ny,cxx_array,cyy_array)
  implicit none                     
    sll_real64                :: eta1,eta1_min,eta2_min
    sll_real64                :: eta2,delta1,delta2
    sll_int32                 :: i,j,nx,ny
    sll_real64, dimension(:,:):: cxx_array,cyy_array
    sll_real64, dimension(1:2,1:2) :: jac_m
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64, dimension(:,:) :: b11
    sll_real64, dimension(:,:) :: b12
    sll_real64, dimension(:,:) :: b21 
    sll_real64, dimension(:,:) :: b22
    
do j=1,ny
 eta2 = eta2_min + real(j-1,f64)*delta2
 do i=1,nx
   eta1 = eta1_min + real(i-1,f64)*delta1
   jac_m  =  transf%jacobian_matrix(eta1,eta2) 
   cxx_array(i,j)= (b11(i,j)*(jac_m(1,2)*jac_m(1,2)+jac_m(2,2)*jac_m(2,2))- &
                   & b12(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2))) &
                   /transf%jacobian(eta1,eta2)
   cyy_array(i,j)= (b22(i,j)*(jac_m(2,1)*jac_m(2,1)+jac_m(1,1)*jac_m(1,1))- &
                   & b21(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2))) &
                   /transf%jacobian(eta1,eta2)                
 enddo
enddo 
end subroutine coefxxyy_array

subroutine coefxy_array(b11,b12,b21,b22,transf,eta1_min,eta2_min, &
                         delta1,delta2,nx,ny,cxy_array)
  implicit none                     
    sll_real64                :: eta1,eta1_min,eta2_min
    sll_real64                :: eta2,delta1,delta2
    sll_real64                :: a12,a21
    sll_int32                 :: i,j,nx,ny
    sll_real64, dimension(:,:):: cxy_array
    sll_real64, dimension(1:2,1:2) :: jac_m
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64, dimension(:,:) :: b11
    sll_real64, dimension(:,:) :: b12
    sll_real64, dimension(:,:) :: b21 
    sll_real64, dimension(:,:) :: b22
    
do j=1,ny
 eta2 = eta2_min + real(j-1,f64)*delta2
 do i=1,nx
   eta1 = eta1_min + real(i-1,f64)*delta1
   jac_m  =  transf%jacobian_matrix(eta1,eta2) 
   a12= b12(i,j)*(jac_m(2,1)*jac_m(2,1)+jac_m(1,1)*jac_m(1,1))- &
                   & b11(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2))
                   
   a21= b21(i,j)*(jac_m(1,2)*jac_m(1,2)+jac_m(2,2)*jac_m(2,2))- &
                   & b22(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2))  
   cxy_array(i,j)= (a12+a21)/transf%jacobian(eta1,eta2)                              
 enddo
enddo 
end subroutine coefxy_array

subroutine a12_a21_array(b11,b12,b21,b22,transf,eta1_min,eta2_min,delta1,delta2,nx,ny,a12_array,a21_array)
    sll_real64                :: eta1,eta1_min,eta2_min
    sll_real64                :: eta2,delta1,delta2
    sll_real64                :: a12,a21
    sll_int32                 :: i,j,nx,ny
    sll_real64, dimension(:,:):: a12_array
    sll_real64, dimension(:,:):: a21_array
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    sll_real64, dimension(:,:) :: b11
    sll_real64, dimension(:,:) :: b12
    sll_real64, dimension(:,:) :: b21 
    sll_real64, dimension(:,:) :: b22
    sll_real64, dimension(1:2,1:2) :: jac_m
    
do j=1,ny
 eta2 = eta2_min + real(j-1,f64)*delta2
 do i=1,nx
   eta1 = eta1_min + real(i-1,f64)*delta1
   jac_m  =  transf%jacobian_matrix(eta1,eta2) 
   a12= b12(i,j)*(jac_m(2,1)*jac_m(2,1)+jac_m(1,1)*jac_m(1,1))- &
                   & b11(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2))
   a12_array(i,j) = a12/transf%jacobian(eta1,eta2)                
   a21= b21(i,j)*(jac_m(1,2)*jac_m(1,2)+jac_m(2,2)*jac_m(2,2))- &
                   & b22(i,j)*(jac_m(2,1)*jac_m(2,2)+jac_m(1,1)*jac_m(1,2)) 
   a21_array(i,j) = a21/transf%jacobian(eta1,eta2)                                             
 enddo
enddo 
end subroutine a12_a21_array

subroutine coefx_array(eta1_min,eta2_min, &
                         delta1,delta2,nx,ny,cx_array)
    sll_real64                :: eta1,eta1_min,eta2_min
    sll_real64                :: eta2,delta1,delta2
    sll_int32                 :: i,j,nx,ny
    sll_real64, dimension(:,:):: cx_array
    
    
do j=1,ny
 eta2 = eta2_min + real(j-1,f64)*delta2
 do i=1,nx
   eta1 = eta1_min + real(i-1,f64)*delta1   
   cx_array(i,j)= cxx_interp%interpolate_from_interpolant_derivative_eta1(eta1,eta2)+ &
                  a21_interp%interpolate_from_interpolant_derivative_eta2(eta1,eta2)                         
 enddo
enddo 
end subroutine coefx_array

subroutine coefy_array(eta1_min,eta2_min, &
                         delta1,delta2,nx,ny,cy_array)
    sll_real64                :: eta1,eta1_min,eta2_min
    sll_real64                :: eta2,delta1,delta2
    sll_int32                 :: i,j,nx,ny
    sll_real64, dimension(:,:):: cy_array
    
do j=1,ny
 eta2 = eta2_min + real(j-1,f64)*delta2
 do i=1,nx
   eta1 = eta1_min + real(i-1,f64)*delta1    
   cy_array(i,j)= cyy_interp%interpolate_from_interpolant_derivative_eta2(eta1,eta2)+ &
                  a12_interp%interpolate_from_interpolant_derivative_eta1(eta1,eta2)                         
 enddo
enddo 
end subroutine coefy_array

!> input pde coefficients at any grid point (x,y) in the solution region
!> (xa.le.x.le.xb,yc.le.y.le.yd) to mud2cr
subroutine coefcr(x,y,cxx,cxy,cyy,cx,cy,ce)
real(8)  :: x,cxx,cx,cxy
real(8)  :: y,cyy,cy,ce
cxx = cxx_interp%interpolate_from_interpolant_value(x,y)
cxy = cxy_interp%interpolate_from_interpolant_value(x,y) 
cyy = cyy_interp%interpolate_from_interpolant_value(x,y) 
cx  = cx_interp%interpolate_from_interpolant_value(x,y)
cy  = cy_interp%interpolate_from_interpolant_value(x,y) 
ce  = ce_interp%interpolate_from_interpolant_value(x,y)
end subroutine coefcr

!!> input x dependent coefficients
!subroutine cofx(x,cxx,cx,cex)
!implicit none
!real(8)  :: x,cxx,cx,cex
!cxx = 1.0_8  !cxx_interp%interpolate_from_interpolant_value(x)
!cx  = 0.0_8 + x - x
!cex = 0.0_8
!end subroutine cofx
!
!!> input y dependent coefficients
!subroutine cofy(y,cyy,cy,cey)
!real(8)  :: y,cyy,cy,cey
!cyy = 1.0_8
!cy  = 0.0_8 + y - y
!cey = 0.0_8
!end subroutine cofy
!
!> input mixed "oblique" derivative b.c. to mud2cr
!> at upper y boundary
subroutine bndcr(kbdy,xory,alfa,beta,gama,gbdy)
integer  :: kbdy
real(8)  :: xory,alfa,beta,gama,gbdy

if (kbdy.eq.2) then

   ! x=xb boundary.
   ! b.c. has the form alfxb(y)*px+betxb(y)*py+gamxb(y)*pe = gbdxb(y)
   ! where xory= y.   alfa,beta,gama,gbdxb corresponding to alfxb(y),
   ! betxb(y),gamxb(y),gbdxb(y) must be output.

   alfa = 0.0_8+0_8*xory
   beta = 0.0_8
   gama = 1.0_8
   gbdy = 0.0_8

end if

if (kbdy.eq.1) then

   ! x=xa boundary.
   ! b.c. has the form alfxa(y)*px+betxa(y)*py+gamxa(y)*pe = gbdxa(y)
   ! where xory= y.   alfa,beta,gama,gbdxb corresponding to alfxa(y),
   ! betxa(y),gamxa(y),gbdxa(y) must be output.

   alfa = 0.0_8
   beta = 0.0_8
   gama = 1.0_8
   gbdy = 0.0_8

end if

end subroutine bndcr

end module sll_m_mudpack_curvilinear

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
