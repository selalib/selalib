!> Solve Maxwell equations on cartesian domain with disconituous Galerkine method:
!> * Gauss Lobatto
!> * Periodic boundary conditions.
module sll_maxwell_2d_diga
#include "sll_maxwell_solvers_macros.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"
#include "sll_file_io.h"
#include "sll_integration.h"

use sll_logical_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_mesh_calculus_2d_module

implicit none

type :: w_vector
   sll_real64, dimension(:), pointer :: ex
   sll_real64, dimension(:), pointer :: ey
   sll_real64, dimension(:), pointer :: hz
end type w_vector

type :: A_matrix
   sll_real64, dimension(3,3) :: v
end type A_matrix

!> Data object to solve Maxwell equations with discontinuous Galerkine
!> method in two dimensions with general coordinates
type :: maxwell_2d_diga
 type(sll_logical_mesh_2d), pointer :: mesh !< Logical mesh
 class(sll_coordinate_transformation_2d_analytic), pointer :: tau  !< Geometric transformation
 sll_int32           :: polarization !< TE or TM
 sll_int32           :: d       !< d of gauss integration
 type(w_vector), dimension(:), pointer :: W
 type(A_matrix), dimension(2)          :: A
end type maxwell_2d_diga


interface initialize
   module procedure initialize_maxwell_2d_diga
end interface initialize
interface solve
   module procedure solve_maxwell_2d_diga
end interface solve

sll_int32, private :: nx, ny
sll_int32, private :: error

contains

!> Initialize Poisson solver object using finite elements method.
!> Indices are shifted from [1:n+1] to [0:n] only inside this 
!> subroutine
!>
!> Construction of the derivative matrix for Gauss-Lobatto 2D
!> \f[ mdiag(i,j)  =  int(Phi_{i,j}.Phi_{k,l})_{[-1;1]}^2  \f]
!> \f[             =  w_i w_j det(\tau'(\hat{\eta}_{i,j} \delta_{i,k}.\delta_{j,l}\f]
subroutine initialize_maxwell_2d_diga( this, tau, d, polarization)
type( maxwell_2d_diga ) :: this !< solver data object
class(sll_coordinate_transformation_2d_analytic), pointer :: tau
sll_int32  :: polarization
sll_int32  :: d
sll_int32  :: nddl
sll_int32  :: ncells
sll_real64 :: xgalo(d+1)
sll_real64 :: wgalo(d+1)
sll_real64 :: dlag(d+1,d+1)
sll_real64 :: mdiag((d+1)*(d+1))
sll_real64 :: Dx((d+1)*(d+1),(d+1)*(d+1))
sll_real64 :: Dy((d+1)*(d+1),(d+1)*(d+1))
sll_real64 :: eta1_p, eta2_p
sll_int32  :: i, j, k, l, ii, jj, kk, ll
sll_real64 :: integral, det
sll_real64 :: jac_mat(2,2)
sll_real64 :: east(2), west(2), north(2), south(2)


this%tau  => tau
this%d = d
this%mesh => tau%mesh
this%polarization = polarization

call tau%write_to_file(SLL_IO_MTV)

nddl   = (d+1)*(d+1)
ncells = this%mesh%num_cells1*this%mesh%num_cells2
xgalo  = gauss_lobatto_points(d+1,-1._f64,1._f64)
write(*,*) " GL points ", xgalo
wgalo  = gauss_lobatto_weights(d+1)
write(*,*) " GL weights ", wgalo
dlag   = gauss_lobatto_derivative_matrix(d+1, -1._f64, 1._f64) 
write(*,*) "Derivative matrix"
do k = 1,d+1
   write(*,*)(dlag(k,l),l=1,d+1)
end do

print*, " test integration "

integral = 0.0_f64
do j = 1, this%mesh%num_cells2
do i = 1, this%mesh%num_cells1
   do ii = 1, d+1
   do jj = 1, d+1
      eta1_p  = (i-0.5+0.5*xgalo(ii))*tau%mesh%delta_eta1+tau%mesh%eta1_min
      eta2_p  = (j-0.5+0.5*xgalo(jj))*tau%mesh%delta_eta2+tau%mesh%eta2_min
      integral = integral + wgalo(ii)*wgalo(jj)*f(eta1_p,eta2_p)
   end do
   end do
end do
end do

write(*,*) "pi*erf(1)**2", sll_pi*erf(1._f64)**2, &
           integral/(this%mesh%num_cells1*this%mesh%num_cells2)


do j = 1, this%mesh%num_cells2+1
   do i = 1, this%mesh%num_cells1+1
   write(11,*) tau%mesh%eta1_min+(i-1)*tau%mesh%delta_eta1, &
               tau%mesh%eta2_min+(j-1)*tau%mesh%delta_eta2, 1.
   end do
   write(11,*) 
end do

mdiag = 0.0

integral = 0.0_f64

do j = 1, this%mesh%num_cells2
do i = 1, this%mesh%num_cells1

   do ii = 1, d+1
   do jj = 1, d+1
      eta1_p  = (i-0.5+0.5*xgalo(ii))*tau%mesh%delta_eta1+tau%mesh%eta1_min
      eta2_p  = (j-0.5+0.5*xgalo(jj))*tau%mesh%delta_eta2+tau%mesh%eta2_min
      jac_mat = tau%jacobian_matrix(eta1_p,eta2_p)
      det     = (jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1))
      mdiag((ii-1)*(d+1)+jj) = wgalo(ii)*wgalo(jj)*det

      
      do ll = 1, d+1
      do kk = 1, d+1
         Dx((ii-1)*(d+1)+jj,(kk-1)*(d+1)+ll) = 0
         Dy((ii-1)*(d+1)+jj,(kk-1)*(d+1)+ll) = 0
      end do
      end do
   end do
   end do

   

end do
end do


this%A(1)%v = reshape((/ 0., 0., 0., 0., 0., 1., 0., 1., 0./), (/3,3/))
this%A(2)%v = reshape((/ 0., 0.,-1., 0., 0., 0.,-1., 0., 0./), (/3,3/))


!k = 0
!do j = 1, this%mesh%num_cells2
!do i = 1, this%mesh%num_cells1
!   k = k + 1
!   do ii = 1, d+1
!   do jj = 1, d+1
!      this%x1(ii,jj,k) = tau%x1(this%eta1(ii,k),this%eta2(jj,k))
!      this%x2(ii,jj,k) = tau%x2(this%eta1(ii,k),this%eta2(jj,k))
!   end do
!   end do
!end do
!end do

SLL_ALLOCATE(this%W(ncells),error)
do k = 1, ncells
   SLL_ALLOCATE(this%W(k)%ex(nddl),error)
   SLL_ALLOCATE(this%W(k)%ey(nddl),error)
   SLL_ALLOCATE(this%W(k)%hz(nddl),error)
end do
 stop

end subroutine initialize_maxwell_2d_diga

!> Solve the maxwell equation
subroutine solve_maxwell_2d_diga( this, ex, ey, bz, dt, jx, jy, rho )

type( maxwell_2d_diga )    :: this !< Maxwell solver object

sll_real64, dimension(:,:) :: ex   !< x electric field
sll_real64, dimension(:,:) :: ey   !< y electric field
sll_real64, dimension(:,:) :: bz   !< z magnetic field
sll_real64, intent(in)     :: dt

sll_real64, dimension(:,:), optional :: jx   !< x current field
sll_real64, dimension(:,:), optional :: jy   !< y current field
sll_real64, dimension(:,:), optional :: rho  !< charge density

end subroutine solve_maxwell_2d_diga


sll_real64 function f( x , y)
sll_real64 :: x, y

f = exp(-(x*x+y*y))

end function f

end module sll_maxwell_2d_diga
