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
#include "sll_utilities.h"

use sll_logical_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_mesh_calculus_2d_module

implicit none

type :: cell_type
   sll_real64 :: eta1_min, eta1_max
   sll_real64 :: eta2_min, eta2_max
   sll_int32  :: i            !             ^
   sll_int32  :: j            !             | n2+
   sll_real64 :: n1_plus(2)   !         ---------
   sll_real64 :: n1_minus(2)  !         |       |
   sll_real64 :: n2_plus(2)   !  n1- <--|       |--> n1+
   sll_real64 :: n2_minus(2)  !         |_______|
   sll_real64 :: l_edge(4)    !             |
end type cell_type            !             V n2-

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
 type(sll_logical_mesh_2d), pointer       :: mesh !< Logical mesh
 class(sll_coordinate_transformation_2d_analytic), pointer :: tau  !< Geometric transformation
 sll_int32                                :: polarization !< TE or TM
 sll_int32                                :: d       !< d of gauss integration
 type(w_vector), dimension(:), pointer    :: W
 type(A_matrix), dimension(2)             :: A
 sll_real64, dimension(:,:),   pointer    :: MassMatrix !< Mass Matrix
 sll_real64, dimension(:,:,:), pointer    :: DxMatrix   !< X Stiffness Matrix
 sll_real64, dimension(:,:,:), pointer    :: DyMatrix   !< Y Stiffness Matrix
 type(cell_type), dimension(:,:), pointer :: cell
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
sll_int32  :: i, j, k, l, ii, jj, kk, ll
sll_real64 :: det
sll_real64 :: jac_mat(2,2)
sll_real64 :: inv_jac_mat(2,2)
sll_real64 :: east(2), west(2), north(2), south(2)
sll_real64 :: eta1_p, eta2_p
sll_real64 :: mdiag

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

call display_matrix(dlag, "f8.3")

do j = 1, this%mesh%num_cells2+1
   do i = 1, this%mesh%num_cells1+1
   write(11,*) tau%mesh%eta1_min+(i-1)*tau%mesh%delta_eta1, &
               tau%mesh%eta2_min+(j-1)*tau%mesh%delta_eta2, 1.
   end do
   write(11,*) 
end do

SLL_CLEAR_ALLOCATE(this%MassMatrix(1:(d+1)*(d+1),1:ncells), error)
SLL_CLEAR_ALLOCATE(this%DxMatrix(1:(d+1)*(d+1),1:(d+1)*(d+1),1:ncells), error)
SLL_CLEAR_ALLOCATE(this%DyMatrix(1:(d+1)*(d+1),1:(d+1)*(d+1),1:ncells), error)
SLL_ALLOCATE(this%cell(this%mesh%num_cells1,this%mesh%num_cells2), error)

k = 0
do j = 1, this%mesh%num_cells2
do i = 1, this%mesh%num_cells1

   k = k+1

   do ii = 1, d+1
   do jj = 1, d+1
      eta1_p  = (i-0.5+0.5*xgalo(ii))*tau%mesh%delta_eta1+tau%mesh%eta1_min
      eta2_p  = (j-0.5+0.5*xgalo(jj))*tau%mesh%delta_eta2+tau%mesh%eta2_min
      jac_mat = tau%jacobian_matrix(eta1_p,eta2_p)
      det     = (jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1))

      mdiag   = wgalo(ii)*wgalo(jj)*det
      this%MassMatrix((ii-1)*(d+1)+jj,k) = mdiag

      do ll = 1, d+1
      do kk = 1, d+1
         if ( jj == ll) &
            this%DxMatrix((ii-1)*(d+1)+jj,(kk-1)*(d+1)+ll,k) = mdiag*dlag(ii,kk)
         if ( ii == kk) &
            this%DyMatrix((ii-1)*(d+1)+jj,(kk-1)*(d+1)+ll,k) = mdiag*dlag(jj,ll)
      end do
      end do
   end do
   end do

   call compute_normals(tau, i, j, d+1, xgalo, wgalo, this%cell(i,j) )

end do
end do
stop

call display_matrix(this%MassMatrix(:,:),"f7.2")
call display_matrix(this%DxMatrix(:,:,1),"f7.2")
call display_matrix(this%DyMatrix(:,:,1),"f7.2")

this%A(1)%v = reshape((/ 0., 0., 0., 0., 0., 1., 0., 1., 0./), (/3,3/))
this%A(2)%v = reshape((/ 0., 0.,-1., 0., 0., 0.,-1., 0., 0./), (/3,3/))

call display_matrix(this%A(1)%v,"f7.2")
call display_matrix(this%A(2)%v,"f7.2")

SLL_ALLOCATE(this%W(ncells),error)
do k = 1, ncells
   SLL_ALLOCATE(this%W(k)%ex(nddl),error)
   SLL_ALLOCATE(this%W(k)%ey(nddl),error)
   SLL_ALLOCATE(this%W(k)%hz(nddl),error)
end do

end subroutine initialize_maxwell_2d_diga

!> Solve the maxwell equation
subroutine solve_maxwell_2d_diga( this, ex, ey, bz, dt, jx, jy, rho )

type( maxwell_2d_diga )    :: this           !< Maxwell solver object

sll_real64, dimension(:,:) :: ex             !< x electric field
sll_real64, dimension(:,:) :: ey             !< y electric field
sll_real64, dimension(:,:) :: bz             !< z magnetic field
sll_real64, intent(in)     :: dt             !< time step

sll_real64, dimension(:,:), optional :: jx   !< x current field
sll_real64, dimension(:,:), optional :: jy   !< y current field
sll_real64, dimension(:,:), optional :: rho  !< charge density

end  subroutine solve_maxwell_2d_diga

subroutine compute_normals(tau, i, j, n, x, w, cell )
class(sll_coordinate_transformation_2d_analytic), pointer :: tau
sll_int32       :: i, j, k, n
sll_real64      :: x(n), w(n)
sll_real64      :: edge_length
sll_real64      :: a, b, c1, c2
sll_real64      :: xk, det
sll_real64      :: eta1, eta2
sll_real64      :: jac_mat(2,2), inv_jac_mat(2,2)
type(cell_type) :: cell

cell%i = i
cell%j = j
!cell corners
cell%eta1_min   = tau%mesh%eta1_min + (i-1)*tau%mesh%delta_eta1
cell%eta1_max   = tau%mesh%eta1_min +  i   *tau%mesh%delta_eta1
cell%eta2_min   = tau%mesh%eta2_min + (j-1)*tau%mesh%delta_eta2
cell%eta2_max   = tau%mesh%eta2_min +  j   *tau%mesh%delta_eta2

cell%l_edge     = 0._f64
cell%n1_minus   = 0._f64
cell%n1_plus    = 0._f64
cell%n2_minus   = 0._f64
cell%n2_plus    = 0._f64

!South
a  = cell%eta1_min
b  = cell%eta1_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, n
   xk = c1*x(k) + c2
   jac_mat          = tau%jacobian_matrix(xk, cell%eta2_min)
   inv_jac_mat      = tau%inverse_jacobian_matrix(xk, cell%eta2_min)
   det              = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   cell%n2_minus    = det*matmul(inv_jac_mat,(/0._f64,-1._f64/))
   cell%l_edge(1)   = cell%l_edge(1) + sqrt(jac_mat(1,1)**2+jac_mat(2,1)**2)*w(k)
end do

!East
a  = cell%eta2_min 
b  = cell%eta2_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, n
   xk = c1*x(k) + c2
   jac_mat          = tau%jacobian_matrix(cell%eta1_max, xk)
   det              = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   inv_jac_mat      = tau%inverse_jacobian_matrix(cell%eta1_max, xk)
   cell%n1_plus     = det*matmul(inv_jac_mat,(/1._f64,0._f64/))
   cell%l_edge(2)   = cell%l_edge(2) + sqrt(jac_mat(1,2)**2+jac_mat(2,2)**2)*w(k)
end do

!North
a  = cell%eta1_min 
b  = cell%eta1_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, n
   xk = c1*x(k) + c2
   jac_mat          = tau%jacobian_matrix(xk, cell%eta2_max)
   det              = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   inv_jac_mat      = tau%inverse_jacobian_matrix(xk, cell%eta2_max)
   cell%n2_plus     = det*matmul(inv_jac_mat,(/0._f64,1._f64/))
   cell%l_edge(3)   = cell%l_edge(3) + sqrt(jac_mat(1,1)**2+jac_mat(2,1)**2)*w(k)
end do

!West
a  = cell%eta2_min 
b  = cell%eta2_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, n
   xk = c1*x(k) + c2
   jac_mat          = tau%jacobian_matrix(cell%eta1_min, xk)
   det              = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   inv_jac_mat      = tau%inverse_jacobian_matrix(cell%eta2_min, xk)
   cell%n1_minus    = det*matmul(inv_jac_mat,(/-1._f64,0._f64/))
   cell%l_edge(4)   = cell%l_edge(4) + sqrt(jac_mat(1,2)**2+jac_mat(2,2)**2)*w(k)
end do

print*, cell%l_edge
print"(2f8.3,8X,a)", cell%n2_minus, " 0 -1"
print"(2f8.3,8X,a)", cell%n1_plus,  " 1  0"
print"(2f8.3,8X,a)", cell%n2_plus,  " 0  1"
print"(2f8.3,8X,a)", cell%n1_minus, "-1  0"


end subroutine compute_normals


end module sll_maxwell_2d_diga


