
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
use sll_boundary_condition_descriptors

implicit none
private

type :: edge_type

   sll_real64                            :: length 
   sll_real64, dimension(:,:), pointer   :: vec_norm
   sll_int32                             :: bc_type

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
   sll_int32                                :: bc_south
   sll_int32                                :: bc_east
   sll_int32                                :: bc_north
   sll_int32                                :: bc_west
   sll_int32                                :: flux_type

end type maxwell_2d_diga

!> Create a Maxwell solver object using Discontinuous Galerkine 
interface initialize
   module procedure initialize_maxwell_2d_diga
end interface initialize

!> Solve Maxwell system
interface solve
   module procedure solve_maxwell_2d_diga
end interface solve

sll_int32               :: error
type(dg_field), pointer :: po
sll_real64, parameter   :: xi = 0.0_f64

public :: initialize, solve

contains

!> Initialize Maxwell solver object using DG method.
subroutine initialize_maxwell_2d_diga( this,         &
                                       tau,          &
                                       degree,       &
                                       polarization, &
                                       bc_south,     &
                                       bc_east,      &
                                       bc_north,     &
                                       bc_west,      &
                                       flux_type)

   type( maxwell_2d_diga )     :: this !< solver data object
   sll_transformation, pointer :: tau
   sll_int32                   :: polarization
   sll_int32                   :: degree
   sll_int32                   :: nddl
   sll_int32                   :: ncells
   sll_real64                  :: x(degree+1)
   sll_real64                  :: w(degree+1)
   sll_real64                  :: dlag(degree+1,degree+1)
   sll_real64                  :: det, dfx, dfy
   sll_real64                  :: j_mat(2,2)
   sll_real64                  :: inv_j(2,2)
   sll_real64                  :: dtau_ij_mat(2,2)
   sll_int32                   :: i, j, k, l, ii, jj, kk, ll
   sll_real64                  :: xa, xb, ya, yb
   sll_int32, intent(in)       :: bc_east
   sll_int32, intent(in)       :: bc_west
   sll_int32, intent(in)       :: bc_north
   sll_int32, intent(in)       :: bc_south
   sll_int32, optional         :: flux_type

   this%tau        => tau
   this%mesh       => tau%mesh
   this%bc_south   =  bc_south
   this%bc_east    =  bc_east
   this%bc_north   =  bc_north
   this%bc_west    =  bc_west

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

   if (present(flux_type)) then
      this%flux_type = flux_type
   else
      this%flux_type = SLL_CENTERED
   end if

   nddl   = (degree+1)*(degree+1)
   ncells = this%nc_eta1*this%nc_eta2

   SLL_ALLOCATE(this%cell(this%nc_eta1,this%nc_eta2), error)

   x    = gauss_lobatto_points(degree+1,0.0_f64,1.0_f64)
   w    = gauss_lobatto_weights(degree+1,0.0_f64,1.0_f64)
   dlag = gauss_lobatto_derivative_matrix(degree+1,x)

   dtau_ij_mat(1,1) = tau%mesh%delta_eta1
   dtau_ij_mat(1,2) = 0.0_f64
   dtau_ij_mat(2,1) = 0.0_f64
   dtau_ij_mat(2,2) = tau%mesh%delta_eta2

   do j = 1, this%nc_eta2   !Loop over cells
   do i = 1, this%nc_eta1

      call compute_normals(tau,bc_south,bc_east,bc_north,bc_west, &
                           i,j,degree,this%cell(i,j))

      SLL_CLEAR_ALLOCATE(this%cell(i,j)%MassMatrix(1:nddl)     , error)
      SLL_CLEAR_ALLOCATE(this%cell(i,j)%DxMatrix(1:nddl,1:nddl), error)
      SLL_CLEAR_ALLOCATE(this%cell(i,j)%DyMatrix(1:nddl,1:nddl), error)

      xa = this%cell(i,j)%eta1_min ; xb  = this%cell(i,j)%eta1_max 
      ya = this%cell(i,j)%eta2_min ; yb  = this%cell(i,j)%eta2_max 

      do jj = 1, degree+1
      do ii = 1, degree+1

         j_mat = tau%jacobian_matrix(xa+x(ii)*this%delta_eta1,&
                                     ya+x(jj)*this%delta_eta2)
         j_mat(:,1) = j_mat(:,1) * this%delta_eta1
         j_mat(:,2) = j_mat(:,2) * this%delta_eta2

         det   = (j_mat(1,1)*j_mat(2,2)-j_mat(1,2)*j_mat(2,1))

         this%cell(i,j)%MassMatrix((ii-1)*(degree+1)+jj) = w(ii)*w(jj)*det

      end do
      end do

      do jj = 1, degree+1
      do ii = 1, degree+1

         do ll = 1, degree+1
         do kk = 1, degree+1


            j_mat = tau%jacobian_matrix(xa+x(kk)*this%delta_eta1,&
                                        ya+x(ll)*this%delta_eta2)
            j_mat(:,1) = j_mat(:,1) * this%delta_eta1
            j_mat(:,2) = j_mat(:,2) * this%delta_eta2

            det   = (j_mat(1,1)*j_mat(2,2)-j_mat(1,2)*j_mat(2,1))

            inv_j(1,1) =   j_mat(2,2) / det
            inv_j(1,2) = - j_mat(1,2) / det
            inv_j(2,1) = - j_mat(2,1) / det
            inv_j(2,2) =   j_mat(1,1) / det

            k = (ii-1)*(degree+1)+jj  
            l = (kk-1)*(degree+1)+ll
             
            dfx = 0.0_f64
            dfy = 0.0_f64

            if (jj == ll) dfx = dlag(kk,ii) !Here we use transpose
            if (ii == kk) dfy = dlag(ll,jj) !dlag (don't know why)

            this%cell(i,j)%DxMatrix(k,l) = & 
               det*w(kk)*w(ll)*(inv_j(1,1)*dfx + inv_j(2,1)*dfy)

            this%cell(i,j)%DyMatrix(k,l) = &
               det*w(kk)*w(ll)*(inv_j(1,2)*dfx + inv_j(2,2)*dfy)

         end do
         end do

      end do
      end do

   end do
   end do

   !call sll_display(this%cell(1,1)%MassMatrix,"f9.4")
   !call sll_display(this%cell(1,1)%DxMatrix,"f9.4")
   !call sll_display(this%cell(1,1)%DyMatrix,"f9.4")

   SLL_CLEAR_ALLOCATE(this%w((degree+1)*(degree+1),4),error)
   SLL_CLEAR_ALLOCATE(this%r((degree+1)*(degree+1),4),error)
   SLL_CLEAR_ALLOCATE(this%f((degree+1)*(degree+1),4),error)

   po => new_dg_field( degree, tau) 

end subroutine initialize_maxwell_2d_diga

!> Solve the maxwell equation
subroutine solve_maxwell_2d_diga( this, fx, fy, fz, dx, dy, dz )

   type( maxwell_2d_diga )  :: this !< Maxwell solver object

   type(dg_field)  :: fx   !< x electric field
   type(dg_field)  :: fy   !< y electric field
   type(dg_field)  :: fz   !< z magnetic field

   type(dg_field)  :: dx  
   type(dg_field)  :: dy  
   type(dg_field)  :: dz  

   sll_int32  :: left, right, node, side, bc_type, flux_type
   sll_int32  :: i, j, k, l, ii, jj, kk
   sll_real64 :: n1
   sll_real64 :: n2
   sll_real64 :: r
   sll_real64 :: flux(4), A_plus(4,4), A_minus(4,4)

   do i = 1, this%nc_eta1
   do j = 1, this%nc_eta2

      do jj = 1, this%degree+1
      do ii = 1, this%degree+1
         kk = (ii-1)*(this%degree+1)+jj
         this%w(kk,1) = fx%array(ii,jj,i,j)
         this%w(kk,2) = fy%array(ii,jj,i,j)
         this%w(kk,3) = fz%array(ii,jj,i,j)
         this%w(kk,4) = po%array(ii,jj,i,j)
      end do
      end do
         
      this%f(:,1) = - matmul(this%cell(i,j)%DyMatrix,this%w(:,3)) &
                 + xi*matmul(this%cell(i,j)%DxMatrix,this%w(:,4))
      this%f(:,2) =   matmul(this%cell(i,j)%DxMatrix,this%w(:,3)) &
                 + xi*matmul(this%cell(i,j)%DyMatrix,this%w(:,4))
      this%f(:,3) = - matmul(this%cell(i,j)%DyMatrix,this%w(:,1)) &
                    + matmul(this%cell(i,j)%DxMatrix,this%w(:,2))
      this%f(:,4) = + matmul(this%cell(i,j)%DxMatrix,this%w(:,1)) &
                    + matmul(this%cell(i,j)%DyMatrix,this%w(:,2))


#ifdef VERBOSE

      print*,"####################################################"
      print*,' (i,j) ', i, j
      print"('Ex=',9f7.3)", this%w(:,1)
      print"('Ey=',9f7.3)", this%w(:,2)
      print"('Bz=',9f7.3)", this%w(:,3)

      print"('dEx=',9f7.3)", this%f(:,1)
      print"('dEy=',9f7.3)", this%f(:,2)
      print"('dBz=',9f7.3)", this%f(:,3)

      print"(a)", 'side'//'  bc '//'left'//'    w(left) ' &
                 &//' f(left) '//'  n1 '//'       n2 '// '       flux '
#endif

      do side = 1, 4 ! Loop over each side of the cell
 
         bc_type = this%cell(i,j)%edge(side)%bc_type
         flux_type = this%flux_type
         !periodic boundary conditions
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
            this%r(kk,1) = fx%array(ii,jj,k,l)
            this%r(kk,2) = fy%array(ii,jj,k,l)
            this%r(kk,3) = fz%array(ii,jj,k,l)
            this%r(kk,4) = po%array(ii,jj,k,l)
         end do
         end do
   
#ifdef VERBOSE
         print*,'--'
#endif
         do node = 1, this%degree+1
   
            left   = dof_local(side, node, this%degree)
            right  = dof_neighbor(side, node, this%degree)

            n1 = this%cell(i,j)%edge(side)%vec_norm(node,1)
            n2 = this%cell(i,j)%edge(side)%vec_norm(node,2)

            bc_type = this%cell(i,j)%edge(side)%bc_type
            
            if (this%flux_type == SLL_CENTERED) then

               flux = 0.5 * (this%w(left,:)+this%r(right,:)) 
               this%f(left,1) = this%f(left,1)+n2*flux(3)-xi*n1*flux(4)
               this%f(left,2) = this%f(left,2)-n1*flux(3)-xi*n2*flux(4)
               this%f(left,3) = this%f(left,3)+n2*flux(1)-n1*flux(2)
               this%f(left,4) = this%f(left,4)-n1*flux(1)-n2*flux(2)

            else
               
               r = sqrt(n1*n1+n2*n2)

               A_plus(1,1) = (n2*n2+xi*n1*n1)/r
               A_plus(1,2) = n2*n1*(xi-1.)/r
               A_plus(1,3) = -n2 
               A_plus(1,4) = xi*n1
               A_plus(2,1) = n2*n1*(xi-1.)/r
               A_plus(2,2) = (n1*n1+xi*n2*n2)/r
               A_plus(2,3) = n1
               A_plus(2,4) = xi*n2
               A_plus(3,1) = -n2
               A_plus(3,2) = n1
               A_plus(3,3) = r 
               A_plus(3,4) = 0.0_f64
               A_plus(4,1) = n1*xi
               A_plus(4,2) = n2*xi
               A_plus(4,3) = 0.0_f64
               A_plus(4,4) = xi*r

               A_minus(1,1) = -(n2*n2+xi*n1*n1)/r
               A_minus(1,2) = -n2*n1*(xi-1.)/r
               A_minus(1,3) = -n2
               A_minus(1,4) = xi*n1
               A_minus(2,1) = -n2*n1*(xi-1.)/r
               A_minus(2,2) = -(n1*n1+xi*n2*n2)/r
               A_minus(2,3) = n1 
               A_minus(2,4) = xi*n2
               A_minus(3,1) = -n2
               A_minus(3,2) =  n1
               A_minus(3,3) = -r 
               A_minus(3,4) = 0.0_f64
               A_minus(4,1) = n1*xi
               A_minus(4,2) = n2*xi
               A_minus(4,3) = 0.0_f64
               A_minus(4,4) = -xi*r

               select case (bc_type)
               case(SLL_CONDUCTOR)
                  flux = this%w(left,:)
                  this%f(left,1) = this%f(left,1)+n2*flux(3)+xi*n1*flux(4)
                  this%f(left,2) = this%f(left,2)-n1*flux(3)+xi*n2*flux(4)
                  this%f(left,3) = this%f(left,3)-n2*flux(1)+n1*flux(2)
                  this%f(left,4) = this%f(left,4)-n1*flux(1)-n2*flux(2)
               case(SLL_SILVER_MULLER)
                  A_plus(1,1) = 0.0_f64
                  A_plus(1,2) = 0.0_f64
                  A_plus(1,3) = -n2 
                  A_plus(1,4) = xi*n1
                  A_plus(2,1) = 0.0_f64
                  A_plus(2,2) = 0.0_f64
                  A_plus(2,3) = n1
                  A_plus(2,4) = xi*n2
                  A_plus(3,1) = -n2
                  A_plus(3,2) = n1
                  A_plus(3,3) = 0.0_f64
                  A_plus(3,4) = 0.0_f64
                  A_plus(4,1) = n1*xi
                  A_plus(4,2) = n2*xi
                  A_plus(4,3) = 0.0_f64
                  A_plus(4,4) = 0.0_f64
                  this%f(left,:) = this%f(left,:)-matmul(A_plus,this%w(left,:))
               case default
                  A_plus(1,1) = 0.0_f64
                  A_plus(1,2) = 0.0_f64
                  A_plus(1,3) = -n2 
                  A_plus(1,4) = xi*n1
                  A_plus(2,1) = 0.0_f64
                  A_plus(2,2) = 0.0_f64
                  A_plus(2,3) = n1
                  A_plus(2,4) = xi*n2
                  A_plus(3,1) = -n2
                  A_plus(3,2) = n1
                  A_plus(3,3) = 0.0_f64
                  A_plus(3,4) = 0.0_f64
                  A_plus(4,1) = n1*xi
                  A_plus(4,2) = n2*xi
                  A_plus(4,3) = 0.0_f64
                  A_plus(4,4) = 0.0_f64
                  flux = 0.5 * (this%w(left,:)+this%r(right,:)) 
                  this%f(left,:) = this%f(left,:)-matmul(A_plus,flux)
               end select

            end if

#ifdef VERBOSE
            print"(3i4,5f10.4)",side, bc_type, left, &
                               this%w(left,3), this%r(right,3), &
                               n1, n2, -n1*this%w(left,3)
#endif
         end do

      end do

#ifdef VERBOSE
      print"('fEx=',9f7.3)", this%f(:,1)
      print"('fEy=',9f7.3)", this%f(:,2)
      print"('fBz=',9f7.3)", this%f(:,3)
#endif

      kk = 0
      do jj = 1, this%degree+1
      do ii = 1, this%degree+1
         kk = (ii-1)*(this%degree+1)+jj
         dx%array(ii,jj,i,j) = this%f(kk,1)/this%cell(i,j)%MassMatrix(kk)
         dy%array(ii,jj,i,j) = this%f(kk,2)/this%cell(i,j)%MassMatrix(kk)
         dz%array(ii,jj,i,j) = this%f(kk,3)/this%cell(i,j)%MassMatrix(kk)
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
      dof_local = (dof-1)*(degree+1)+1
   case(EAST)
      dof_local = degree*(degree+1)+dof
   case(NORTH)
      dof_local = dof*(degree+1)
   case(WEST)
      dof_local = dof
   end select

end function dof_local

function dof_neighbor(edge,dof,degree)

   sll_int32 :: dof_neighbor
   sll_int32 :: edge, dof, degree
   
   select case(edge)
   case(SOUTH)
      dof_neighbor = dof*(degree+1)
   case(EAST)
      dof_neighbor = dof 
   case(NORTH)
      dof_neighbor = (dof-1)*(degree+1)+1
   case(WEST)
      dof_neighbor = degree*(degree+1)+dof
      !dof_neighbor = (degree+1)*(degree+1)-(dof-1)*(degree+1)
   end select

end function dof_neighbor

!> Compute cell normals
subroutine compute_normals(tau, bc_south, bc_east, bc_north, bc_west, &
                           i, j, d, cell )

   class(sll_coordinate_transformation_2d_analytic), pointer :: tau
   type(cell_type) :: cell
   sll_int32       :: i, j, d
   sll_real64      :: x(d+1), w(d+1), vec_norm(d+1,2)
   sll_real64      :: a, b, c1, c2
   sll_real64      :: xk, wk
   sll_real64      :: jac_mat(2,2)
   sll_real64      :: co_jac_mat(2,2)
   sll_real64      :: dtau_ij_mat(2,2)
   sll_real64      :: jac_mat_sll(2,2)
   sll_real64      :: length
   sll_int32       :: side
   sll_int32       :: bc_south
   sll_int32       :: bc_east
   sll_int32       :: bc_north
   sll_int32       :: bc_west
   sll_int32       :: k
   
   cell%i = i
   cell%j = j
   cell%eta1_min = tau%mesh%eta1_min + (i-1)*tau%mesh%delta_eta1
   cell%eta2_min = tau%mesh%eta2_min + (j-1)*tau%mesh%delta_eta2
   
   cell%eta1_max = cell%eta1_min + tau%mesh%delta_eta1
   cell%eta2_max = cell%eta2_min + tau%mesh%delta_eta2

   dtau_ij_mat(1,1) = tau%mesh%delta_eta1
   dtau_ij_mat(1,2) = 0.0_f64
   dtau_ij_mat(2,1) = 0.0_f64
   dtau_ij_mat(2,2) = tau%mesh%delta_eta2

   do side = 1, 4
      SLL_CLEAR_ALLOCATE(cell%edge(side)%vec_norm(1:d+1,1:2),error)
      cell%edge(side)%bc_type = SLL_INTERIOR
   end do

   if (j ==                   1) cell%edge(SOUTH)%bc_type = bc_south
   if (i == tau%mesh%num_cells1) cell%edge(EAST)%bc_type  = bc_east
   if (j == tau%mesh%num_cells2) cell%edge(NORTH)%bc_type = bc_north
   if (i ==                   1) cell%edge(WEST)%bc_type  = bc_west
   
   x = gauss_lobatto_points(d+1)
   w = gauss_lobatto_weights(d+1)

   length = 0._f64
   a  = cell%eta1_min
   b  = cell%eta1_max 
   c1 = 0.5_f64 * (b-a)
   c2 = 0.5_f64 * (b+a)
   do k = 1, d+1
      xk = c1*x(k) + c2
      wk = c1*w(k)
      jac_mat_sll     = tau%jacobian_matrix(xk, cell%eta2_min)
      jac_mat         = matmul(jac_mat_sll, dtau_ij_mat)
      co_jac_mat(1,1) =  jac_mat(2,2)
      co_jac_mat(1,2) = -jac_mat(2,1)
      co_jac_mat(2,1) = -jac_mat(1,2)
      co_jac_mat(2,2) =  jac_mat(1,1)
      length          = length + sqrt(jac_mat_sll(1,1)**2+jac_mat_sll(2,1)**2)*wk
      vec_norm(k,:)   = matmul(co_jac_mat,(/0._f64,-1._f64/))*wk
   end do

   cell%edge(SOUTH)%length = length
   cell%edge(SOUTH)%vec_norm = vec_norm/tau%mesh%delta_eta1
   
   length = 0._f64
   a  = cell%eta2_min 
   b  = cell%eta2_max 
   c1 = 0.5_f64 * (b-a)
   c2 = 0.5_f64 * (b+a)
   do k = 1, d+1
      xk = c1*x(k) + c2
      wk = c1*w(k)
      jac_mat_sll     = tau%jacobian_matrix(cell%eta1_max, xk)
      jac_mat         = matmul(jac_mat_sll, dtau_ij_mat)
      co_jac_mat(1,1) =  jac_mat(2,2)
      co_jac_mat(1,2) = -jac_mat(2,1)
      co_jac_mat(2,1) = -jac_mat(1,2)
      co_jac_mat(2,2) =  jac_mat(1,1)
      length          = length + sqrt(jac_mat_sll(1,2)**2+jac_mat_sll(2,2)**2)*wk
      vec_norm(k,:)   = matmul(co_jac_mat,(/1._f64, 0._f64/))*wk
   end do

   cell%edge(EAST)%length = length
   cell%edge(EAST)%vec_norm = vec_norm/tau%mesh%delta_eta2
   
   length = 0._f64
   a  = cell%eta1_min 
   b  = cell%eta1_max 
   c1 = 0.5_f64 * (b-a)
   c2 = 0.5_f64 * (b+a)
   do k = 1, d+1
      xk = c1*x(k) + c2
      wk = c1*w(k)
      jac_mat_sll     = tau%jacobian_matrix(xk, cell%eta2_max)
      jac_mat         = matmul(jac_mat_sll, dtau_ij_mat)
      co_jac_mat(1,1) =  jac_mat(2,2)
      co_jac_mat(1,2) = -jac_mat(2,1)
      co_jac_mat(2,1) = -jac_mat(1,2)
      co_jac_mat(2,2) =  jac_mat(1,1)
      length          = length + sqrt(jac_mat_sll(1,1)**2+jac_mat_sll(2,1)**2)*wk
      vec_norm(k,:)   = matmul(co_jac_mat,(/0._f64, 1._f64/))*wk
   end do

   cell%edge(NORTH)%length = length
   cell%edge(NORTH)%vec_norm = vec_norm/tau%mesh%delta_eta1
   
   length = 0._f64
   a  = cell%eta2_min 
   b  = cell%eta2_max 
   c1 = 0.5_f64 * (b-a)
   c2 = 0.5_f64 * (b+a)
   do k = 1, d+1
      xk = c1*x(k) + c2
      wk = c1*w(k)
      jac_mat_sll     = tau%jacobian_matrix(cell%eta1_min, xk)
      jac_mat         = matmul(jac_mat_sll, dtau_ij_mat)
      co_jac_mat(1,1) =  jac_mat(2,2)
      co_jac_mat(1,2) = -jac_mat(2,1)
      co_jac_mat(2,1) = -jac_mat(1,2)
      co_jac_mat(2,2) =  jac_mat(1,1)
      length          = length + sqrt(jac_mat_sll(1,2)**2+jac_mat_sll(2,2)**2)*wk
      vec_norm(k,:)   = matmul(co_jac_mat,(/-1._f64, 0._f64/))*wk
   end do

   cell%edge(WEST)%length = length
   cell%edge(WEST)%vec_norm = vec_norm/tau%mesh%delta_eta2

end subroutine compute_normals

end module sll_maxwell_2d_diga
