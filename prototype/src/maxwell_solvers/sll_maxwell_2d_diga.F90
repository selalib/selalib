
#define sll_transformation class(sll_coordinate_transformation_2d_analytic)

!> Solve Maxwell equations on cartesian domain with Disconituous Galerkine method:
!> * Gauss Lobatto for integration formula
!> * Periodic boundary conditions.
module sll_maxwell_2d_diga

#include "sll_maxwell_solvers_macros.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"
#include "sll_integration.h"
#include "sll_utilities.h"
#include "sll_assert.h"

use sll_logical_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_dg_fields

implicit none
private

type :: edge_type

   sll_real64                            :: length 
   sll_real64, dimension(:,:), pointer   :: vec_norm

end type edge_type

!> This derived type contains information about a mesh cell
type :: cell_type

   sll_int32                           :: i,j         !< indices 
   sll_real64                          :: eta1_min    !< left side
   sll_real64                          :: eta1_max    !< right side
   sll_real64                          :: eta2_min    !< bottom side
   sll_real64                          :: eta2_max    !< top side
   type(edge_type)                     :: edge(4)     !< normal vectors  
   sll_real64, dimension(:),   pointer :: MassMatrix  !< Mass Matrix
   sll_real64, dimension(:,:), pointer :: DxMatrix    !< X Stiffness Matrix
   sll_real64, dimension(:,:), pointer :: DyMatrix    !< Y Stiffness Matrix

end type cell_type


!> method in two dimensions with general coordinates
type, public :: maxwell_2d_diga

   sll_transformation, pointer              :: tau  !< transformation
   type(sll_logical_mesh_2d), pointer       :: mesh !< Logical mesh
   sll_int32                                :: polarization !< TE or TM
   sll_int32                                :: degree !< degree of gauss integration
   type(cell_type), dimension(:,:), pointer :: cell !< mesh cells
   sll_int32                                :: nc_eta1
   sll_int32                                :: nc_eta2
   sll_real64                               :: eta1_min
   sll_real64                               :: eta1_max
   sll_real64                               :: delta_eta1
   sll_real64                               :: eta2_min
   sll_real64                               :: eta2_max
   sll_real64                               :: delta_eta2
   sll_real64, dimension(:,:), pointer      :: f
   sll_real64, dimension(:,:), pointer      :: w
   sll_real64, dimension(:,:), pointer      :: r

end type maxwell_2d_diga

!> Create a Maxwell solver object using Discontinuous Galerkine 
interface initialize
   module procedure initialize_maxwell_2d_diga
end interface initialize

!> Solve Maxell system
interface solve
   module procedure solve_maxwell_2d_diga
end interface solve

sll_int32                                :: error
type(dg_field), pointer                  :: po
sll_real64, parameter                    :: xi = 0.0_f64
sll_real64, dimension(4,4)               :: A1
sll_real64, dimension(4,4)               :: A2

public :: initialize, solve

contains

!> Initialize Maxwell solver object using DG method.
subroutine initialize_maxwell_2d_diga( this, tau, degree, polarization)

   type( maxwell_2d_diga )     :: this !< solver data object
   sll_transformation, pointer :: tau
   sll_int32                   :: polarization
   sll_int32                   :: degree
   sll_int32                   :: nddl
   sll_int32                   :: ncells
   sll_real64                  :: x(degree+1)
   sll_real64                  :: y(degree+1)
   sll_real64                  :: wx(degree+1)
   sll_real64                  :: wy(degree+1)
   sll_real64                  :: dlagx(degree+1,degree+1)
   sll_real64                  :: dlagy(degree+1,degree+1)
   sll_real64                  :: det
   sll_real64                  :: jac_mat(2,2)
   sll_real64                  :: inv_jac_mat(2,2)
   sll_real64                  :: mdiag
   sll_int32                   :: i, j, ii, jj, kk, ll
   sll_real64                  :: xa, xb, ya, yb

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

   this%degree       =  degree
   this%polarization = polarization

   nddl   = (degree+1)*(degree+1)
   ncells = this%nc_eta1*this%nc_eta2

   SLL_ALLOCATE(this%cell(this%nc_eta1,this%nc_eta2), error)


   do i = 1, this%nc_eta1
   do j = 1, this%nc_eta2

      call compute_normals(tau,i,j,degree,this%cell(i,j))

      SLL_CLEAR_ALLOCATE(this%cell(i,j)%MassMatrix(1:nddl)     , error)
      SLL_CLEAR_ALLOCATE(this%cell(i,j)%DxMatrix(1:nddl,1:nddl), error)
      SLL_CLEAR_ALLOCATE(this%cell(i,j)%DyMatrix(1:nddl,1:nddl), error)

      xa = this%cell(i,j)%eta1_min ; xb  = this%cell(i,j)%eta1_max 
      ya = this%cell(i,j)%eta2_min ; yb  = this%cell(i,j)%eta2_max 
      x     = gauss_lobatto_points(degree+1,xa,xb)
      y     = gauss_lobatto_points(degree+1,ya,yb)
      wx    = gauss_lobatto_weights(degree+1,xa,xb)
      wy    = gauss_lobatto_weights(degree+1,ya,yb)
      dlagx = gauss_lobatto_derivative_matrix(degree+1,x)
      dlagy = gauss_lobatto_derivative_matrix(degree+1,y)

      call sll_display(dlagx,"f9.4")

      do ii = 1, degree+1
      do jj = 1, degree+1

         jac_mat     = tau%jacobian_matrix(x(ii),y(jj))
         inv_jac_mat = tau%inverse_jacobian_matrix(x(ii),y(jj))
         det         = (jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1))
         mdiag       = wx(ii)*wy(jj)*det

         this%cell(i,j)%MassMatrix((ii-1)*(degree+1)+jj) = mdiag

         do ll = 1, degree+1
         do kk = 1, degree+1

            if (jj == ll) &
               this%cell(i,j)%DxMatrix((ii-1)*(degree+1)+jj,(kk-1)*(degree+1)+ll) &
                  = mdiag*dlagx(ii,kk)

            if (ii == kk) &
               this%cell(i,j)%DyMatrix((ii-1)*(degree+1)+jj,(kk-1)*(degree+1)+ll) &
                  = mdiag*dlagy(jj,ll)
         end do
         end do

      end do
      end do


   end do
   end do

   call sll_display(this%cell(2,2)%MassMatrix,"f9.4")
   call sll_display(this%cell(2,2)%DxMatrix,"f9.4")
   call sll_display(this%cell(2,2)%DyMatrix,"f9.4")

   SLL_CLEAR_ALLOCATE(this%w((degree+1)*(degree+1),4),error)
   SLL_CLEAR_ALLOCATE(this%r((degree+1)*(degree+1),4),error)
   SLL_CLEAR_ALLOCATE(this%f((degree+1)*(degree+1),4),error)

   A1 = reshape((/ 0.0_f64, 0.0_f64, 0.0_f64, xi,       &
                   0.0_f64, 0.0_f64, 1.0_f64, 0.0_f64,  &
                   0.0_f64, 1.0_f64, 0.0_f64, 0.0_f64,  &
                             xi, 1.0_f64, 0.0_f64, 0.0_f64   /), (/4,4/))
   
   A2 = reshape((/ 0.0_f64, 0.0_f64, -1.0_f64, 0.0_f64, &
                   0.0_f64, 0.0_f64,  0.0_f64,      xi, &
                  -1.0_f64, 0.0_f64,  0.0_f64, 0.0_f64, &
                   0.0_f64,      xi,  0.0_f64, 0.0_f64  /), (/4,4/))

   po => new_dg_field( degree, tau) 

end subroutine initialize_maxwell_2d_diga

!> Solve the maxwell equation
subroutine solve_maxwell_2d_diga( this, ex, ey, bz, dex, dey, dbz )

   type( maxwell_2d_diga )  :: this !< Maxwell solver object

   type(dg_field)  :: ex   !< x electric field
   type(dg_field)  :: ey   !< y electric field
   type(dg_field)  :: bz   !< z magnetic field

   type(dg_field)  :: dex  
   type(dg_field)  :: dey  
   type(dg_field)  :: dbz  

   sll_int32  :: left, right, node, side
   sll_int32  :: i, j, k, l, ii, jj, kk
   sll_real64 :: vec_n1
   sll_real64 :: vec_n2
   sll_real64 :: flux(4), surf
   sll_real64 :: w(this%degree+1)

   w = gauss_lobatto_weights(this%degree+1, 0.0_f64, 1.0_f64)


   do i = 1, this%nc_eta1
   do j = 1, this%nc_eta2

      surf = this%mesh%delta_eta1 * this%mesh%delta_eta2

      do jj = 1, this%degree+1
      do ii = 1, this%degree+1
         kk = (ii-1)*(this%degree+1)+jj
         this%w(kk,1) = ex%array(ii,jj,i,j)
         this%w(kk,2) = ey%array(ii,jj,i,j)
         this%w(kk,3) = bz%array(ii,jj,i,j)
         this%w(kk,4) = po%array(ii,jj,i,j)
      end do
      end do
         
      this%f(:,1) = - matmul(this%cell(i,j)%DyMatrix,this%w(:,3)) &
               - xi * matmul(this%cell(i,j)%DxMatrix,this%w(:,4))

      this%f(:,2) =   matmul(this%cell(i,j)%DxMatrix,this%w(:,3)) &
               - xi * matmul(this%cell(i,j)%DyMatrix,this%w(:,4))

      this%f(:,3) = - matmul(this%cell(i,j)%DyMatrix,this%w(:,1)) &
                    + matmul(this%cell(i,j)%DxMatrix,this%w(:,2))

      this%f(:,4) = - xi * matmul(this%cell(i,j)%DxMatrix,this%w(:,1)) &
                    + xi * matmul(this%cell(i,j)%DyMatrix,this%w(:,2))

      do side = 1, 4 ! Loop over each side of the cell
 
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

         do jj = 1, this%degree+1
         do ii = 1, this%degree+1
            kk = (ii-1)*(this%degree+1)+jj
            this%r(kk,1) = ex%array(ii,jj,k,l)
            this%r(kk,2) = ey%array(ii,jj,k,l)
            this%r(kk,3) = bz%array(ii,jj,k,l)
            this%r(kk,4) = po%array(ii,jj,k,l)
         end do
         end do
   
         !Compute the fluxes on edge points
         do node = 1, this%degree+1
   
            left   = dof_local(side, node, this%degree)
            right  = dof_neighbor(side, node, this%degree)

            vec_n1 = this%cell(i,j)%edge(side)%vec_norm(node,1)
            vec_n2 = this%cell(i,j)%edge(side)%vec_norm(node,2)
   
            flux = (this%r(right,:)-this%w(left,:))
            flux = flux * 0.5 * (this%r(right,:)+this%w(left,:))

            flux = flux * w(node) * this%cell(i,j)%edge(side)%length
  
            this%f(left,1) = this%f(left,1)+vec_n2*flux(3)!+xi*vec_n1*flux(4)
            this%f(left,2) = this%f(left,2)-vec_n1*flux(3)!+xi*vec_n2*flux(4)
            this%f(left,3) = this%f(left,3)+vec_n2*flux(1)-vec_n1*flux(2)
            !this%f(left,4) = this%f(left,4)+vec_n1*flux(1)+vec_n2*flux(2)
            
         end do
   
      end do
   
      this%f(:,1) = this%f(:,1) / this%cell(i,j)%MassMatrix(:)
      this%f(:,2) = this%f(:,2) / this%cell(i,j)%MassMatrix(:)
      this%f(:,3) = this%f(:,3) / this%cell(i,j)%MassMatrix(:)
      this%f(:,4) = this%f(:,4) / this%cell(i,j)%MassMatrix(:)

      do jj = 1, this%degree+1
      do ii = 1, this%degree+1
         kk = (ii-1)*(this%degree+1)+jj
         dex%array(ii,jj,i,j) = - this%f(kk,1)
         dey%array(ii,jj,i,j) = - this%f(kk,2)
         dbz%array(ii,jj,i,j) = - this%f(kk,3)
      end do
      end do
   
   end do
   end do
   
end subroutine solve_maxwell_2d_diga

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
      !global_ddl = degree*(degree+1)+local_ddl
   case(WEST)
      dof_local = degree*(degree+1)+1-(dof-1)*(degree+1)
      !dof_local = (local_ddl-1)*(degree+1) + 1
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

!> Compute cell normals
subroutine compute_normals(tau, i, j, d, cell )

   class(sll_coordinate_transformation_2d_analytic), pointer :: tau
   sll_int32       :: i, j, k, d
   sll_real64      :: x(d+1), w(d+1)
   sll_real64      :: a, b, c1, c2
   sll_real64      :: xk, det
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
   
   do side = 1, 4
      SLL_CLEAR_ALLOCATE(cell%edge(side)%vec_norm(1:d+1,1:2),error)
   end do
   
   x = gauss_lobatto_points(d+1)
   w = gauss_lobatto_weights(d+1)
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
      cell%edge(SOUTH)%vec_norm(k,:) = det*matmul(inv_jac_mat,(/0._f64,-1._f64/))
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
      cell%edge(EAST)%vec_norm(k,:) = det*matmul(inv_jac_mat,(/1._f64, 0._f64/))
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
      cell%edge(NORTH)%vec_norm(k,:) = det*matmul(inv_jac_mat,(/0._f64, 1._f64/))
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
      cell%edge(WEST)%vec_norm(k,:) = det*matmul(inv_jac_mat,(/-1._f64, 0._f64/))
      length                 = length + sqrt(jac_mat(1,2)**2+jac_mat(2,2)**2)*w(k)
   end do
   cell%edge(WEST)%length = length * c1

end subroutine compute_normals

end module sll_maxwell_2d_diga

