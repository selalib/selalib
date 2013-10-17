!> Solve Maxwell equations on cartesian domain with Disconituous Galerkine method:
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

implicit none
private

type :: W_vector
   sll_real64, dimension(:), pointer :: ex
   sll_real64, dimension(:), pointer :: ey
   sll_real64, dimension(:), pointer :: hz
end type W_vector

interface operator(/)
  module procedure W_divide_by_DiagMatrix
end interface operator(/)

interface operator(*)
  module procedure W_multiply_by_Matrix
end interface operator(*)

!> This derived type contains information about a mesh cell
type :: cell_type
   sll_int32                             :: i,j         !< indices 
   sll_real64                            :: eta1_min    !< left side
   sll_real64                            :: eta1_max    !< right side
   sll_real64                            :: eta2_min    !< bottom side
   sll_real64                            :: eta2_max    !< top side
   sll_real64                            :: n(4,2)      !< normal vectors  
   sll_real64                            :: l_edge(4)   !< edge length 
   sll_real64, dimension(:), pointer     :: MassMatrix  !< Mass Matrix
   sll_real64, dimension(:,:), pointer   :: DxMatrix    !< X Stiffness Matrix
   sll_real64, dimension(:,:), pointer   :: DyMatrix    !< Y Stiffness Matrix
   type(W_vector)                        :: W           !< solution fields
end type cell_type


!> Data object to solve Maxwell equations with discontinuous Galerkine
!> method in two dimensions with general coordinates
type, public :: maxwell_2d_diga
 type(sll_logical_mesh_2d), pointer                        :: mesh !< Logical mesh
 class(sll_coordinate_transformation_2d_analytic), pointer :: tau  !< Geometric transformation
 sll_int32                                                 :: polarization !< TE or TM
 sll_int32                                                 :: d    !< d of gauss integration
 type(cell_type), dimension(:,:), pointer                  :: cell !< mesh cells
end type maxwell_2d_diga

interface initialize
   module procedure initialize_maxwell_2d_diga
end interface initialize

interface solve
   module procedure solve_maxwell_2d_diga
end interface solve

sll_int32                  :: nx, ny
sll_int32                  :: error
sll_int32                  :: i, j, k, l, ii, jj, kk, ll
type(w_vector)             :: V
sll_real64, dimension(3,3) :: A1
sll_real64, dimension(3,3) :: A2

public :: initialize, solve

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
sll_real64 :: det
sll_real64 :: jac_mat(2,2)
sll_real64 :: inv_jac_mat(2,2)
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

call sll_display(dlag, "f8.3")

do j = 1, this%mesh%num_cells2+1
   do i = 1, this%mesh%num_cells1+1
   write(11,*) tau%mesh%eta1_min+(i-1)*tau%mesh%delta_eta1, &
               tau%mesh%eta2_min+(j-1)*tau%mesh%delta_eta2, 1.
   end do
   write(11,*) 
end do

SLL_ALLOCATE(this%cell(this%mesh%num_cells1,this%mesh%num_cells2), error)


do i = 1, this%mesh%num_cells1
do j = 1, this%mesh%num_cells2

   SLL_CLEAR_ALLOCATE(this%cell(i,j)%MassMatrix(1:(d+1)*(d+1))            , error)
   SLL_CLEAR_ALLOCATE(this%cell(i,j)%DxMatrix(1:(d+1)*(d+1),1:(d+1)*(d+1)), error)
   SLL_CLEAR_ALLOCATE(this%cell(i,j)%DyMatrix(1:(d+1)*(d+1),1:(d+1)*(d+1)), error)

   do ii = 1, d+1
   do jj = 1, d+1
      eta1_p  = (i-0.5+0.5*xgalo(ii))*tau%mesh%delta_eta1+tau%mesh%eta1_min
      eta2_p  = (j-0.5+0.5*xgalo(jj))*tau%mesh%delta_eta2+tau%mesh%eta2_min
      jac_mat = tau%jacobian_matrix(eta1_p,eta2_p)
      det     = (jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1))

      mdiag   = wgalo(ii)*wgalo(jj)*det
      this%cell(i,j)%MassMatrix((ii-1)*(d+1)+jj) = mdiag

      do ll = 1, d+1
      do kk = 1, d+1
         if ( jj == ll) &
            this%cell(i,j)%DxMatrix((ii-1)*(d+1)+jj,(kk-1)*(d+1)+ll) = mdiag*dlag(ii,kk)
         if ( ii == kk) &
            this%cell(i,j)%DyMatrix((ii-1)*(d+1)+jj,(kk-1)*(d+1)+ll) = mdiag*dlag(jj,ll)
      end do
      end do
   end do
   end do

   call compute_normals(tau, i, j, d+1, xgalo, wgalo, this%cell(i,j) )

   !this%A1 = reshape((/ 0., 0., 0., 0., 0., this%cell%n1_plus, 0., this%cell%n1_plus, 0./), (/3,3/))
   !this%A2 = reshape((/ 0., 0.,-this%n2_plus, 0., 0., 0.,-this%n2_plus, 0., 0./), (/3,3/))

   SLL_CLEAR_ALLOCATE(this%cell(i,j)%W%ex(1:(d+1)*(d+1)),error)
   SLL_CLEAR_ALLOCATE(this%cell(i,j)%W%ey(1:(d+1)*(d+1)),error)
   SLL_CLEAR_ALLOCATE(this%cell(i,j)%W%hz(1:(d+1)*(d+1)),error)

end do
end do

call sll_display(this%cell(1,1)%MassMatrix(:),"f7.2")
call sll_display(this%cell(1,1)%DxMatrix(:,:),"f7.2")
call sll_display(this%cell(1,1)%DyMatrix(:,:),"f7.2")

SLL_CLEAR_ALLOCATE(V%ex(1:d+1),error)
SLL_CLEAR_ALLOCATE(V%ey(1:d+1),error)
SLL_CLEAR_ALLOCATE(V%hz(1:d+1),error)

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
sll_int32 :: edge, R, L

do i = 1, this%mesh%num_cells1
do j = 1, this%mesh%num_cells2

   !this%W = this%W + dt * (  matmul(this%A1,this%DxMatrix(:,:,k))*this%W &
   !                        + matmul(this%A2*this%DyMatrix(:,:,k))*this%W )

   do edge = 1, 4 ! Loop over edges
 
      do k = 1, this%d+1
         L = dof_local(edge, k, this%d)
         R = dof_neighbor(edge, k, this%d)
         V%ex(k) = 0.5*(this%cell(i,j)%W%ex(L)+this%cell(i,j)%W%ex(R))
      end do

   end do

   this%cell(i,j)%W = this%cell(i,j)%W / this%cell(i,j)%MassMatrix

end do
end do
stop


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
cell%eta1_min = tau%mesh%eta1_min + (i-1)*tau%mesh%delta_eta1
cell%eta2_min = tau%mesh%eta2_min + (j-1)*tau%mesh%delta_eta2

cell%eta1_max = cell%eta1_min + tau%mesh%delta_eta1
cell%eta2_max = cell%eta2_min + tau%mesh%delta_eta2

cell%l_edge   = 0._f64
cell%n        = 0._f64

a  = cell%eta1_min
b  = cell%eta1_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, n
   xk = c1*x(k) + c2
   jac_mat            = tau%jacobian_matrix(xk, cell%eta2_min)
   inv_jac_mat        = tau%inverse_jacobian_matrix(xk, cell%eta2_min)
   det                = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   cell%n(SOUTH,:)    = det*matmul(inv_jac_mat,(/0._f64,-1._f64/))
   cell%l_edge(SOUTH) = cell%l_edge(SOUTH) + sqrt(jac_mat(1,1)**2+jac_mat(2,1)**2)*w(k)
end do
cell%l_edge(SOUTH) = cell%l_edge(SOUTH) * c1

a  = cell%eta2_min 
b  = cell%eta2_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, n
   xk = c1*x(k) + c2
   jac_mat           = tau%jacobian_matrix(cell%eta1_max, xk)
   det               = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   inv_jac_mat       = tau%inverse_jacobian_matrix(cell%eta1_max, xk)
   cell%n(EAST,:)    = det*matmul(inv_jac_mat,(/1._f64, 0._f64/))
   cell%l_edge(EAST) = cell%l_edge(EAST) + sqrt(jac_mat(1,2)**2+jac_mat(2,2)**2)*w(k)
end do
cell%l_edge(EAST) = cell%l_edge(EAST) * c1

a  = cell%eta1_min 
b  = cell%eta1_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, n
   xk = c1*x(k) + c2
   jac_mat            = tau%jacobian_matrix(xk, cell%eta2_max)
   det                = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   inv_jac_mat        = tau%inverse_jacobian_matrix(xk, cell%eta2_max)
   cell%n(NORTH,:)    = det*matmul(inv_jac_mat,(/0._f64, 1._f64/))
   cell%l_edge(NORTH) = cell%l_edge(NORTH) + sqrt(jac_mat(1,1)**2+jac_mat(2,1)**2)*w(k)
end do
cell%l_edge(NORTH) = cell%l_edge(NORTH) * c1

a  = cell%eta2_min 
b  = cell%eta2_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, n
   xk = c1*x(k) + c2
   jac_mat           = tau%jacobian_matrix(cell%eta1_min, xk)
   det               = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   inv_jac_mat       = tau%inverse_jacobian_matrix(cell%eta2_min, xk)
   cell%n(WEST,:)    = det*matmul(inv_jac_mat,(/-1._f64, 0._f64/))
   cell%l_edge(WEST) = cell%l_edge(WEST) + sqrt(jac_mat(1,2)**2+jac_mat(2,2)**2)*w(k)
end do
cell%l_edge(WEST) = cell%l_edge(WEST) * c1

print"(/,4f8.3)",      cell%eta1_min, cell%eta1_max, cell%eta2_min, cell%eta2_max
print"(4f8.3)",        cell%l_edge
print"(2f8.3,8X,a)",   cell%n(:,SOUTH),  " 0 -1"
print"(2f8.3,8X,a)",   cell%n(:,EAST) ,  " 1  0"
print"(2f8.3,8X,a)",   cell%n(:,NORTH),  " 0  1"
print"(2f8.3,8X,a,/)", cell%n(:,WEST) ,  "-1  0"


end subroutine compute_normals

function W_divide_by_DiagMatrix( W1, D) result(W2)

  type(W_vector), intent(in)           :: W1
  sll_real64, dimension(:), intent(in) :: D
  type(W_vector)                       :: W2
  
  W2%Ex = W1%Ex / D
  W2%Ey = W1%Ey / D
  W2%Hz = W1%Hz / D

end function W_divide_by_DiagMatrix

function W_multiply_by_Matrix( Mat, W1) result(W2)

  type(W_vector), intent(in)             :: W1
  sll_real64, dimension(:,:), intent(in) :: Mat
  type(W_vector)                         :: W2
 
  W2%Ex = matmul(Mat,W1%Ex)
  W2%Ey = matmul(Mat,W1%Ey)
  W2%Hz = matmul(Mat,W1%Hz)

end function W_multiply_by_Matrix

function dof_local(edge,dof,degree)
sll_int32 :: dof_local
sll_int32 :: edge, dof, degree

select case(edge)
case(SOUTH)
   dof_local = dof
case(EAST)
   dof_local = dof*(degree+1)
case(NORTH)
   dof_local = (degree+1)*(degree+1)-dof+1 
case(WEST)
   dof_local = degree*(degree+1)+1-(dof-1)*(degree+1)
end select

end function dof_local

function dof_neighbor(edge,dof,degree)
sll_int32 :: dof_neighbor
sll_int32 :: edge, dof, degree

select case(edge)
case(SOUTH)
   dof_neighbor = degree*(degree+1)+dof
case(EAST)
   dof_neighbor = (dof-1)*(degree+1)+1
case(NORTH)
   dof_neighbor = degree+1-dof+1
case(WEST)
   dof_neighbor = (degree+1)*(degree+1)-(dof-1)*(degree+1)
end select

end function dof_neighbor

end module sll_maxwell_2d_diga

