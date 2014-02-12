
#define sll_transformation class(sll_coordinate_transformation_2d_analytic)

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
use sll_dg_fields

implicit none
private

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
   sll_real64, dimension(:), pointer   :: MassMatrix  !< Mass Matrix
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
   sll_real64, dimension(:,:), pointer      :: w_vector              
   sll_real64, dimension(:,:), pointer      :: f_vector              
   sll_real64, dimension(:,:), pointer      :: r_vector              
   sll_real64, dimension(3,3)               :: A1
   sll_real64, dimension(3,3)               :: A2

end type maxwell_2d_diga

!> Create a Maxwell solver object using Discontinuous Galerkine 
interface initialize
   module procedure initialize_maxwell_2d_diga
end interface initialize

!> Solve Maxell system
interface solve
   module procedure solve_maxwell_2d_diga
end interface solve

sll_int32, private :: error

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
   sll_real64                  :: dlag(degree+1,degree+1)
   sll_real64                  :: det
   sll_real64                  :: jac_mat(2,2)
   sll_real64                  :: eta1_p
   sll_real64                  :: eta2_p
   sll_real64                  :: mdiag
   sll_int32                   :: i, j, ii, jj, kk, ll
   sll_real64                  :: xgalo(degree+1)
   sll_real64                  :: wgalo(degree+1)

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

   call tau%write_to_file(SLL_IO_MTV)

   nddl   = (degree+1)*(degree+1)
   ncells = this%nc_eta1*this%nc_eta2

   dlag   = gauss_lobatto_derivative_matrix(degree+1, -1._f64, 1._f64) 

   !call sll_display(dlag, "f8.3")

   SLL_ALLOCATE(this%cell(this%nc_eta1,this%nc_eta2), error)

   xgalo  = gauss_lobatto_points(degree+1,-1._f64,1._f64)
   wgalo  = gauss_lobatto_weights(degree+1)

   do i = 1, this%nc_eta1
   do j = 1, this%nc_eta2

      SLL_CLEAR_ALLOCATE(this%cell(i,j)%MassMatrix(1:nddl)     , error)
      SLL_CLEAR_ALLOCATE(this%cell(i,j)%DxMatrix(1:nddl,1:nddl), error)
      SLL_CLEAR_ALLOCATE(this%cell(i,j)%DyMatrix(1:nddl,1:nddl), error)

      do ii = 1, degree+1
      do jj = 1, degree+1

         eta1_p  = (i-0.5+0.5*xgalo(ii))*this%delta_eta1+this%eta1_min
         eta2_p  = (j-0.5+0.5*xgalo(jj))*this%delta_eta2+this%eta2_min
         jac_mat = tau%jacobian_matrix(eta1_p,eta2_p)
         det     = (jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1))

         mdiag   = wgalo(ii)*wgalo(jj)*det
         this%cell(i,j)%MassMatrix((ii-1)*(degree+1)+jj) = mdiag

         do ll = 1, degree+1
         do kk = 1, degree+1

            if ( jj == ll) &
               this%cell(i,j)%DxMatrix((ii-1)*(degree+1)+jj,(kk-1)*(degree+1)+ll) &
                  = mdiag*dlag(ii,kk)

            if ( ii == kk) &
               this%cell(i,j)%DyMatrix((ii-1)*(degree+1)+jj,(kk-1)*(degree+1)+ll) &
                  = mdiag*dlag(jj,ll)
         end do
         end do

      end do
      end do

      call compute_normals(tau, i, j, degree, xgalo,wgalo, this%cell(i,j) )

   end do
   end do

   !call sll_display(this%cell(1,1)%MassMatrix(:),"f7.2")
   !call sll_display(this%cell(1,1)%DxMatrix(:,:),"f7.2")
   !call sll_display(this%cell(1,1)%DyMatrix(:,:),"f7.2")

   SLL_CLEAR_ALLOCATE(this%w_vector((degree+1)*(degree+1),3),error)
   SLL_CLEAR_ALLOCATE(this%f_vector((degree+1)*(degree+1),3),error)
   SLL_CLEAR_ALLOCATE(this%r_vector((degree+1)*(degree+1),3),error)

   this%A1 = reshape((/ 0._f64, 0._f64,  0._f64,  &
                        0._f64, 0._f64,  1._f64,  &
                        0._f64, 1._f64,  0._f64/), (/3,3/))
   
   this%A2 = reshape((/ 0._f64, 0._f64, -1._f64,  &
                        0._f64, 0._f64,  0._f64,  &
                       -1._f64, 0._f64,  0._f64/), (/3,3/))

end subroutine initialize_maxwell_2d_diga


!> Solve the maxwell equation
subroutine solve_maxwell_2d_diga( this, ex, ey, bz, dt, jx, jy, rho )

   type( maxwell_2d_diga )  :: this !< Maxwell solver object

   type(dg_field), target   :: ex   !< x electric field
   type(dg_field), target   :: ey   !< y electric field
   type(dg_field), target   :: bz   !< z magnetic field
   sll_real64, intent(in)   :: dt   !< time step

   type(dg_field), optional :: jx   !< x current field
   type(dg_field), optional :: jy   !< y current field
   type(dg_field), optional :: rho  !< charge density

   sll_int32                :: left, right, node, side
   sll_int32                :: i, j, k, l, ii, jj
   sll_real64               :: flux(3)


   print*, dt
   !Loop over cells
   do i = 1, this%nc_eta1
   do j = 1, this%nc_eta2

      do jj = 1, this%degree+1
      do ii = 1, this%degree+1
         k = (ii-1)*(this%degree+1)+jj
         this%w_vector(k,1) = ex%array(ii,jj,i,j)
         this%w_vector(k,2) = ey%array(ii,jj,i,j)
         this%w_vector(k,3) = bz%array(ii,jj,i,j)
      end do
      end do
         
      if(present(jx) .and. present(jy) .and. present(rho)) then
         do jj = 1, this%degree+1
         do ii = 1, this%degree+1
            this%r_vector(k,1) = jx%array(ii,jj,i,j)
            this%r_vector(k,2) = jy%array(ii,jj,i,j)
            this%r_vector(k,3) = rho%array(ii,jj,i,j)
         end do
         end do
      end if

      this%f_vector(:,1) =  &
             matmul(this%cell(i,j)%DxMatrix,this%w_vector(:,1)) &
           + matmul(this%cell(i,j)%DyMatrix,this%w_vector(:,1))
      this%f_vector(:,2) =  &
             matmul(this%cell(i,j)%DxMatrix,this%w_vector(:,2)) &
           + matmul(this%cell(i,j)%DyMatrix,this%w_vector(:,2))
      this%f_vector(:,3) =  &
             matmul(this%cell(i,j)%DxMatrix,this%w_vector(:,3)) &
           + matmul(this%cell(i,j)%DyMatrix,this%w_vector(:,3))

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
         do node = 1, this%degree+1
   
            left  = dof_local(side, node, this%degree)
            right = dof_neighbor(side, node, this%degree)
   
            flux (:) = 0.5*(this%w_vector(left, :) &
                     +      this%w_vector(right,:))
   
            !V%field(node,:) = &
            !   this%cell(i,j)%edge(side)%n(node,1)*matmul(A1,V%field(node,:)) &
            ! + this%cell(i,j)%edge(side)%n(node,2)*matmul(A2,V%field(node,:))
   
         end do
   
      end do
   
      !this%cell(i,j)%W = this%cell(i,j)%W / this%cell(i,j)%MassMatrix

   
   end do
   end do


end  subroutine solve_maxwell_2d_diga

!function W_divide_by_DiagMatrix( W1, D) result(W2)
!
!  type(dg_field), intent(in)           :: W1
!  sll_real64, dimension(:), intent(in) :: D
!  type(dg_field)                       :: W2
!  
!  W2%field(:,1) = W1%field(:,1) / D(:)
!  W2%field(:,2) = W1%field(:,2) / D(:)
!  W2%field(:,3) = W1%field(:,3) / D(:)
!
!end function W_divide_by_DiagMatrix
!
!function W_multiply_by_Matrix( AMat, W1) result(W2)
!
!  type(dg_field), intent(in)             :: W1
!  sll_real64, dimension(3,3), intent(in) :: AMat
!  type(dg_field)                         :: W2
!  sll_int32                              :: i, n
!
!  do i = 1, W2%nddl
!     W2%field(i,:) = matmul(AMat,W1%field(i,:))
!  end do
!
!end function W_multiply_by_Matrix

!   sll_int32                               :: degree
!   sll_transformation, pointer             :: tau  
!   sll_real64, dimension(:,:,:,:), pointer :: array
!   sll_real64, dimension(:), pointer       :: xgalo
!   sll_real64, dimension(:), pointer       :: wgalo


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

!> Compute cell normals
subroutine compute_normals(tau, i, j, d, x, w, cell )

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

sll_int32 function global_ddl( degree, iface, local_ddl)

   sll_int32, intent(in)  :: degree
   sll_int32, intent(in)  :: iface
   sll_int32, intent(in)  :: local_ddl
   
   select case(iface)
   case(SOUTH)
   global_ddl = local_ddl
   case(EAST)
   global_ddl = local_ddl*(degree+1)
   case(NORTH)
   global_ddl = degree*(degree+1)+local_ddl
   case(WEST)
   global_ddl = (local_ddl-1)*(degree+1) + 1
   end select 

end function global_ddl

end module sll_maxwell_2d_diga

