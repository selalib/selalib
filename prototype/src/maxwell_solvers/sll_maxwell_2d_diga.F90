!> Solve Maxwell equations on cartesian domain with Disconituous Galerkine method:
!> * Gauss Lobatto for integration formula
!> * Periodic boundary conditions.
module sll_maxwell_2d_diga
#include "sll_maxwell_solvers_macros.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"
#include "sll_file_io.h"
#include "sll_integration.h"
#include "sll_utilities.h"
#include "sll_assert.h"

use sll_logical_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations

implicit none
private

type :: W_vector
   sll_real64, dimension(:,:), pointer :: f
end type W_vector

interface operator(/)
  module procedure W_divide_by_DiagMatrix
end interface operator(/)

interface operator(*)
  module procedure W_multiply_by_Matrix
end interface operator(*)

interface operator(+)
  module procedure W_vector_addition
end interface operator(+)

type :: edge_type
   sll_real64                            :: length 
   sll_real64, dimension(:,:), pointer   :: n
end type edge_type

!> This derived type contains information about a mesh cell
type :: cell_type
   sll_int32                           :: i,j         !< indices 
   sll_real64                          :: eta1_min    !< left side
   sll_real64                          :: eta1_max    !< right side
   sll_real64                          :: eta2_min    !< bottom side
   sll_real64                          :: eta2_max    !< top side
   type(edge_type)                     :: edge(4)     !< normal vectors  
   type(W_vector)                      :: W           !< fields values
   sll_real64, dimension(:), pointer   :: MassMatrix  !< Mass Matrix
   sll_real64, dimension(:,:), pointer :: DxMatrix    !< X Stiffness Matrix
   sll_real64, dimension(:,:), pointer :: DyMatrix    !< Y Stiffness Matrix
end type cell_type

!> Data object to solve Maxwell equations with discontinuous Galerkine
!> method in two dimensions with general coordinates
type, public :: maxwell_2d_diga

   type(sll_logical_mesh_2d), pointer              :: mesh !< Logical mesh
   class(sll_coordinate_transformation_2d_analytic), &
      pointer                                      :: tau  !< transformation
   sll_int32                                       :: polarization !< TE or TM
   sll_int32                                       :: d    !< d of gauss integration
   type(cell_type), dimension(:,:), pointer        :: cell !< mesh cells
   sll_int32                                       :: nc_eta1
   sll_int32                                       :: nc_eta2
   sll_real64                                      :: eta1_min
   sll_real64                                      :: eta1_max
   sll_real64                                      :: delta_eta1
   sll_real64                                      :: eta2_min
   sll_real64                                      :: eta2_max
   sll_real64                                      :: delta_eta2
   sll_real64, dimension(:,:,:,:), pointer         :: bz
   sll_real64, dimension(:), allocatable           :: xgalo
   sll_real64, dimension(:), allocatable           :: wgalo

end type maxwell_2d_diga

!> Create a Maxwell solver object using Discontinuous Galerkine 
interface initialize
   module procedure initialize_maxwell_2d_diga
end interface initialize

!> Solve Maxell system
interface solve
   module procedure solve_maxwell_2d_diga
end interface solve

sll_int32                  :: error
type(w_vector)             :: V
type(w_vector)             :: F
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
subroutine initialize_maxwell_2d_diga( this, tau, d, init_function, polarization)
type( maxwell_2d_diga ) :: this !< solver data object
class(sll_coordinate_transformation_2d_analytic), pointer :: tau
sll_int32  :: polarization
sll_int32  :: d
sll_int32  :: nddl
sll_int32  :: ncells
sll_real64 :: dlag(d+1,d+1)
sll_real64 :: det
sll_real64 :: jac_mat(2,2)
sll_real64 :: inv_jac_mat(2,2)
sll_real64 :: eta1_p
sll_real64 :: eta2_p
sll_real64 :: mdiag
sll_int32  :: i, j, k, l, ii, jj, kk, ll
sll_real64, external ::  init_function

this%tau   => tau
this%mesh  => tau%mesh

this%nc_eta1    = tau%mesh%num_cells1
this%nc_eta2    = tau%mesh%num_cells2
this%eta1_min   = tau%mesh%eta1_min
this%eta2_min   = tau%mesh%eta2_min
this%eta1_max   = tau%mesh%eta1_max
this%eta2_max   = tau%mesh%eta2_max
this%delta_eta1 = tau%mesh%delta_eta1
this%delta_eta2 = tau%mesh%delta_eta2

this%d            =  d
this%polarization = polarization

call tau%write_to_file(SLL_IO_MTV)

nddl   = (d+1)*(d+1)
ncells = this%nc_eta1*this%nc_eta2

SLL_ALLOCATE(this%xgalo(d+1),error)
SLL_ALLOCATE(this%wgalo(d+1),error)

this%xgalo  = gauss_lobatto_points(d+1,-1._f64,1._f64)
this%wgalo  = gauss_lobatto_weights(d+1)
dlag   = gauss_lobatto_derivative_matrix(d+1, -1._f64, 1._f64) 

this%bz => new_diga_field_2d(this, init_function)

call diga_plot_2d(this, this%bz)
stop
call sll_display(dlag, "f8.3")

SLL_ALLOCATE(this%cell(this%nc_eta1,this%nc_eta2), error)

do i = 1, this%nc_eta1
do j = 1, this%nc_eta2

   SLL_CLEAR_ALLOCATE(this%cell(i,j)%MassMatrix(1:(d+1)*(d+1))            , error)
   SLL_CLEAR_ALLOCATE(this%cell(i,j)%DxMatrix(1:(d+1)*(d+1),1:(d+1)*(d+1)), error)
   SLL_CLEAR_ALLOCATE(this%cell(i,j)%DyMatrix(1:(d+1)*(d+1),1:(d+1)*(d+1)), error)

   do ii = 1, d+1
   do jj = 1, d+1
      eta1_p  = (i-0.5+0.5*this%xgalo(ii))*this%delta_eta1+this%eta1_min
      eta2_p  = (j-0.5+0.5*this%xgalo(jj))*this%delta_eta2+this%eta2_min
      jac_mat = tau%jacobian_matrix(eta1_p,eta2_p)
      det     = (jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1))

      mdiag   = this%wgalo(ii)*this%wgalo(jj)*det
      this%cell(i,j)%MassMatrix((ii-1)*(d+1)+jj) = mdiag

      do ll = 1, d+1
      do kk = 1, d+1
         if ( jj == ll) &
            this%cell(i,j)%DxMatrix((ii-1)*(d+1)+jj,(kk-1)*(d+1)+ll) &
               = mdiag*dlag(ii,kk)
         if ( ii == kk) &
            this%cell(i,j)%DyMatrix((ii-1)*(d+1)+jj,(kk-1)*(d+1)+ll) &
               = mdiag*dlag(jj,ll)
      end do
      end do
   end do
   end do

   call compute_normals(tau, i, j, d, this%xgalo, this%wgalo, this%cell(i,j) )

   SLL_CLEAR_ALLOCATE(this%cell(i,j)%W%f(1:(d+1)*(d+1),3),error)


end do
end do

!call sll_display(this%cell(1,1)%MassMatrix(:),"f7.2")
!call sll_display(this%cell(1,1)%DxMatrix(:,:),"f7.2")
!call sll_display(this%cell(1,1)%DyMatrix(:,:),"f7.2")
SLL_CLEAR_ALLOCATE(F%f(1:(d+1)*(d+1),3),error)
SLL_CLEAR_ALLOCATE(V%f(1:d+1,3),error)

A1 = reshape((/ 0._f64, 0._f64, 0._f64,  &
                0._f64, 0._f64, 1._f64,  &
                0._f64, 1._f64, 0._f64/), (/3,3/))

A2 = reshape((/ 0._f64, 0._f64, -1._f64, &
                0._f64, 0._f64,  0._f64, &
               -1._f64, 0._f64,  0._f64/), (/3,3/))

end subroutine initialize_maxwell_2d_diga

!> Compute cell normals
subroutine compute_normals(tau, i, j, d, x, w, cell )

class(sll_coordinate_transformation_2d_analytic), pointer :: tau
sll_int32       :: i, j, k, d
sll_real64      :: x(d+1), w(d+1)
sll_real64      :: edge_length
sll_real64      :: a, b, c1, c2
sll_real64      :: xk, det
sll_real64      :: eta1, eta2
sll_real64      :: jac_mat(2,2), inv_jac_mat(2,2)
type(cell_type) :: cell
sll_real64      :: length
sll_int32       :: side

cell%i = i
cell%j = j
cell%eta1_min = tau%mesh%eta1_min + (i-1)*tau%mesh%delta_eta1
cell%eta2_min = tau%mesh%eta2_min + (j-1)*tau%mesh%delta_eta2

cell%eta1_max = cell%eta1_min + tau%mesh%delta_eta1
cell%eta2_max = cell%eta2_min + tau%mesh%delta_eta2

print"(/,4f8.3)", cell%eta1_min, cell%eta1_max, cell%eta2_min, cell%eta2_max

do side = 1, 4
   SLL_CLEAR_ALLOCATE(cell%edge(side)%n(1:d+1,1:2),error)
end do

length = 0._f64
a  = cell%eta1_min
b  = cell%eta1_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, d+1
   xk = c1*x(k) + c2
   jac_mat            = tau%jacobian_matrix(xk, cell%eta2_min)
   inv_jac_mat        = tau%inverse_jacobian_matrix(xk, cell%eta2_min)
   det                = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   cell%edge(SOUTH)%n(k,:) = det*matmul(inv_jac_mat,(/0._f64,-1._f64/))
   length             = length + sqrt(jac_mat(1,1)**2+jac_mat(2,1)**2)*w(k)
end do
cell%edge(SOUTH)%length = length * c1

length = 0._f64
a  = cell%eta2_min 
b  = cell%eta2_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, d+1
   xk = c1*x(k) + c2
   jac_mat                = tau%jacobian_matrix(cell%eta1_max, xk)
   det                    = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   inv_jac_mat            = tau%inverse_jacobian_matrix(cell%eta1_max, xk)
   cell%edge(EAST)%n(k,:) = det*matmul(inv_jac_mat,(/1._f64, 0._f64/))
   length                 = length + sqrt(jac_mat(1,2)**2+jac_mat(2,2)**2)*w(k)
end do
cell%edge(EAST)%length = length * c1

length = 0._f64
a  = cell%eta1_min 
b  = cell%eta1_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, d+1
   xk = c1*x(k) + c2
   jac_mat                 = tau%jacobian_matrix(xk, cell%eta2_max)
   det                     = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   inv_jac_mat             = tau%inverse_jacobian_matrix(xk, cell%eta2_max)
   cell%edge(NORTH)%n(k,:) = det*matmul(inv_jac_mat,(/0._f64, 1._f64/))
   length                  = length + sqrt(jac_mat(1,1)**2+jac_mat(2,1)**2)*w(k)
end do
cell%edge(NORTH)%length = length * c1

length = 0._f64
a  = cell%eta2_min 
b  = cell%eta2_max 
c1 = 0.5_f64 * (b-a)
c2 = 0.5_f64 * (b+a)
do k = 1, d+1
   xk = c1*x(k) + c2
   jac_mat                = tau%jacobian_matrix(cell%eta1_min, xk)
   det                    = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1) 
   inv_jac_mat            = tau%inverse_jacobian_matrix(cell%eta2_min, xk)
   cell%edge(WEST)%n(k,:) = det*matmul(inv_jac_mat,(/-1._f64, 0._f64/))
   length                 = length + sqrt(jac_mat(1,2)**2+jac_mat(2,2)**2)*w(k)
end do
cell%edge(WEST)%length = length * c1

end subroutine compute_normals

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
sll_int32 :: edge, left, right, node, side
sll_int32 :: i, j, k, l

!Loop over cells
do i = 1, this%nc_eta1
do j = 1, this%nc_eta2

   F%f(:,1) = matmul(this%cell(i,j)%DxMatrix,this%cell(i,j)%W%f(:,1)) &
            + matmul(this%cell(i,j)%DyMatrix,this%cell(i,j)%W%f(:,1))
   F%f(:,2) = matmul(this%cell(i,j)%DxMatrix,this%cell(i,j)%W%f(:,2)) &
            + matmul(this%cell(i,j)%DyMatrix,this%cell(i,j)%W%f(:,2))
   F%f(:,3) = matmul(this%cell(i,j)%DxMatrix,this%cell(i,j)%W%f(:,3)) &
            + matmul(this%cell(i,j)%DyMatrix,this%cell(i,j)%W%f(:,3))

   do side = 1, 4 ! Loop over edges
 
      !boundary conditions are periodic
      select case(side)
      case(SOUTH)
         k = i
         l = 1+modulo(j-2,this%nc_eta2) 
      case(EAST)
         k = 1+modulo(i  ,this%nc_eta1)
         l = j
      case(NORTH)
         k = i
         l = 1+modulo(j  ,this%nc_eta2)
      case(WEST)
         k = 1+modulo(i-2,this%nc_eta1)
         l = j
      end select

      !Compute the fluxes on edge points
      do node = 1, this%d+1

         left  = dof_local(side, node, this%d)
         right = dof_neighbor(side, node, this%d)

         !V%f(k,:) = 0.5*(this%cell(i,j)%W%f(left,:)+this%cell(k,l)%W%f(right,:))

         V%f(node,:) = this%cell(i,j)%edge(side)%n(node,1)*matmul(A1,V%f(node,:)) &
                     + this%cell(i,j)%edge(side)%n(node,2)*matmul(A2,V%f(node,:))

      end do

   end do

   !this%cell(i,j)%W = this%cell(i,j)%W / this%cell(i,j)%MassMatrix

end do
end do

end  subroutine solve_maxwell_2d_diga

function W_divide_by_DiagMatrix( W1, D) result(W2)

  type(W_vector), intent(in)           :: W1
  sll_real64, dimension(:), intent(in) :: D
  type(W_vector)                       :: W2
  
  W2%f(:,1) = W1%f(:,1) / D(:)
  W2%f(:,2) = W1%f(:,2) / D(:)
  W2%f(:,3) = W1%f(:,3) / D(:)

end function W_divide_by_DiagMatrix

function W_multiply_by_Matrix( AMat, W1) result(W2)

  type(W_vector), intent(in)             :: W1
  sll_real64, dimension(3,3), intent(in) :: AMat
  type(W_vector)                         :: W2
  sll_int32                              :: i, n

  n = size(W2%f,1)
 
  do i = 1, n
     W2%f(i,:) = matmul(AMat,W1%f(i,:))
  end do

end function W_multiply_by_Matrix

function W_vector_addition( W1, W2) result(W3)

  type(W_vector), intent(in) :: W1
  type(W_vector), intent(in) :: W2
  type(W_vector)             :: W3

  W3%f = W1%f + W2%f

end function W_vector_addition

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


subroutine diga_plot_2d( this, field )

   type(maxwell_2d_diga)  :: this
   sll_real64, intent(in) :: field(:,:,:,:)
   sll_int32              :: file_id
   sll_int32              :: gnu_id
   sll_real64             :: eta1, eta2
   sll_real64             :: offset(2)
   sll_int32              :: i, j, ii, jj
   sll_int32              :: icell
   character(len=4)       :: ccell

   call sll_ascii_file_create("test.gnu", gnu_id, error)

   icell = 0
   do i = 1, this%nc_eta1
   do j = 1, this%nc_eta2
 
      icell = icell+1
      call int2string(icell, ccell)
      if (icell == 1) then
         write(gnu_id,"(a)",advance='no') "splot 'test"//ccell//".dat' w l"
      else
         write(gnu_id,"(a)",advance='no') ",'test"//ccell//".dat' w l "
      end if

      call sll_ascii_file_create("test"//ccell//".dat", file_id, error)

      offset(1) = this%eta1_min + (i-1)*this%delta_eta1
      offset(2) = this%eta2_min + (j-1)*this%delta_eta2
      do ii = 1, this%d+1
      do jj = 1, this%d+1
         eta1 = offset(1) + 0.5 * (this%xgalo(ii) + 1.0) * this%delta_eta1
         eta2 = offset(2) + 0.5 * (this%xgalo(jj) + 1.0) * this%delta_eta2
         write(file_id,*) this%tau%x1(eta1,eta2), &
                          this%tau%x2(eta1,eta2), &
                          sngl(field(i,j,ii,jj))
      end do
      write(file_id,*)
      end do
      close(file_id)

   end do
   end do

   write(gnu_id,*)
   close(gnu_id)
   
end subroutine diga_plot_2d

function new_diga_field_2d( this, init_function) result(field)

   type(maxwell_2d_diga)          :: this
   sll_real64, external, optional :: init_function
   sll_real64, pointer            :: field(:,:,:,:)
   
   SLL_CLEAR_ALLOCATE(field(1:this%nc_eta1,1:this%nc_eta2,1:this%d+1,1:this%d+1),error)

   if (present(init_function)) then
      call initialize_diga_field_2d( this, field, init_function, 0.0_f64) 
   end if

end function new_diga_field_2d
   
subroutine initialize_diga_field_2d( this, field, init_function, time) 

   type(maxwell_2d_diga)          :: this
   sll_real64, pointer            :: field(:,:,:,:)
   sll_real64, external           :: init_function
   sll_real64                     :: time
   sll_real64                     :: offset(2)
   sll_real64                     :: eta1
   sll_real64                     :: eta2
   sll_int32                      :: i, j, ii, jj
   
   SLL_ASSERT(associated(field))

   do i = 1, this%nc_eta1
   do j = 1, this%nc_eta1
      offset(1) = this%eta1_min + (i-1)*this%delta_eta1
      offset(2) = this%eta2_min + (j-1)*this%delta_eta2
      do ii = 1, this%d+1
      do jj = 1, this%d+1
         eta1 = offset(1) + 0.5 * (this%xgalo(ii) + 1.0) * this%delta_eta1
         eta2 = offset(2) + 0.5 * (this%xgalo(jj) + 1.0) * this%delta_eta2
         field(i,j,ii,jj) = init_function(this%tau%x1(eta1,eta2), &
                                          this%tau%x2(eta1,eta2), &
                                          time)
      end do
      end do
   end do
   end do

end subroutine initialize_diga_field_2d



end module sll_maxwell_2d_diga

