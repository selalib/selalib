
#define sll_transformation class(sll_coordinate_transformation_2d_base)

!> Solve Maxwell equations on cartesian domain with Disconituous Galerkine method:
!> * Gauss Lobatto for integration formula
!> * Periodic boundary conditions.
module sll_maxwell_2d_diga_parallel

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
use sll_maxwell_2d_diga_parallel

implicit none

!> method in two dimensions with general coordinates
type, public :: maxwell_2d_diga_parallel

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
   type(dg_field), pointer                  :: po
   sll_real64                               :: xi 

end type maxwell_2d_diga_parallel

interface new
   module procedure new_maxwell_2d_diga_parallel
end interface new

!> Create a Maxwell solver object using Discontinuous Galerkine 
interface initialize
   module procedure initialize_maxwell_2d_diga_parallel
end interface initialize


!> Solve Maxwell system
interface solve
   module procedure solve_maxwell_2d_diga_parallel
end interface solve

contains

function new_maxwell_2d_diga_parallel(       &
                               tau,          &
                               degree,       &
                               polarization, &
                               bc_south,     &
                               bc_east,      &
                               bc_north,     &
                               bc_west,      &
                               flux_type) result(this)

   type( maxwell_2d_diga ), pointer :: this !< solver data object
   sll_transformation, pointer      :: tau
   sll_int32                        :: polarization
   sll_int32                        :: degree
   sll_int32, intent(in)            :: bc_east
   sll_int32, intent(in)            :: bc_west
   sll_int32, intent(in)            :: bc_north
   sll_int32, intent(in)            :: bc_south
   sll_int32, optional              :: flux_type

   SLL_ALLOCATE(this,error)

   call initialize_maxwell_2d_diga( this,         &
                                    tau,          &
                                    degree,       &
                                    polarization, &
                                    bc_south,     &
                                    bc_east,      &
                                    bc_north,     &
                                    bc_west,      &
                                    flux_type)

 end function new_maxwell_2d_diga_parallel

!> Initialize Maxwell solver object using DG method.
subroutine initialize_maxwell_2d_diga_parallel(      &
                                       this,         &
                                       tau,          &
                                       degree,       &
                                       polarization, &
                                       bc_south,     &
                                       bc_east,      &
                                       bc_north,     &
                                       bc_west,      &
                                       flux_type)

   type(maxwell_2d_diga)       :: this !< solver data object
   sll_transformation, pointer :: tau
   sll_int32                   :: polarization
   sll_int32                   :: degree
   sll_int32, intent(in)       :: bc_east
   sll_int32, intent(in)       :: bc_west
   sll_int32, intent(in)       :: bc_north
   sll_int32, intent(in)       :: bc_south
   sll_int32, optional         :: flux_type

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

   this%tau        => tau
   ! Please undo this 'fix' whenever it is decided that gfortran 4.6 is no
   ! longer supported.
   !   this%mesh       => tau%get_logical_mesh()
   this%mesh => tau%mesh
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

   this%xi           = 0.0_f64
   this%degree       = degree
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

   this%po => new_dg_field( degree, tau) 

end subroutine initialize_maxwell_2d_diga_parallel


!> Solve the maxwell equation
subroutine solve_maxwell_2d_diga_parallel( this, fx, fy, fz, dx, dy, dz )

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
   sll_real64 :: flux(4)
   sll_real64 :: A(4,4)
   sll_real64 :: A_p(4,4)
   sll_real64 :: A_m(4,4)
   sll_real64 :: xi

   xi = this%xi

   do i = 1, this%nc_eta1
   do j = 1, this%nc_eta2

      do jj = 1, this%degree+1
      do ii = 1, this%degree+1
         kk = (ii-1)*(this%degree+1)+jj
         this%w(kk,1) = fx%array(ii,jj,i,j)
         this%w(kk,2) = fy%array(ii,jj,i,j)
         this%w(kk,3) = fz%array(ii,jj,i,j)
         this%w(kk,4) = this%po%array(ii,jj,i,j)
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
            this%r(kk,4) = this%po%array(ii,jj,k,l)
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
            r  = sqrt(n1*n1+n2*n2)

            bc_type = this%cell(i,j)%edge(side)%bc_type

            A(1,:) = [0.0_f64, 0.0_f64,     -n2,   xi*n1]
            A(2,:) = [0.0_f64, 0.0_f64,      n1,   xi*n2]
            A(3,:) = [    -n2,      n1, 0.0_f64, 0.0_f64]
            A(4,:) = [  xi*n1,   xi*n2, 0.0_f64, 0.0_f64]

            A_p(1,:)=[ (n2*n2+xi*n1*n1)/r,    n2*n1*(xi-1.)/r,    -n2, xi*n1]
            A_p(2,:)=[    n2*n1*(xi-1.)/r, (n1*n1+xi*n2*n2)/r,     n1, xi*n2]
            A_p(3,:)=[                -n2,                 n1,      r,0._f64]
            A_p(4,:)=[              n1*xi,              n2*xi,0.0_f64,  xi*r]
    
            A_m(1,:)=[-(n2*n2+xi*n1*n1)/r,   -n2*n1*(xi-1.)/r,    -n2, xi*n1]
            A_m(2,:)=[   -n2*n1*(xi-1.)/r,-(n1*n1+xi*n2*n2)/r,     n1, xi*n2]
            A_m(3,:)=[                -n2,                 n1,     -r,0._f64]
            A_m(4,:)=[              n1*xi,              n2*xi,0.0_f64, -xi*r]

            select case (bc_type)

            case(SLL_CONDUCTOR)

               A(1,:) = [-1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64]
               A(2,:) = [ 0.0_f64,-1.0_f64, 0.0_f64, 0.0_f64]
               A(3,:) = [ 0.0_f64, 0.0_f64, 1.0_f64, 0.0_f64]
               A(4,:) = [-2*xi*n1,-2*xi*n2, 0.0_f64,-1.0_f64]

               flux = matmul(A,this%w(left,:))
               flux = matmul(A_m, flux) + matmul(A_p,this%w(left,:))

               this%f(left,:) = this%f(left,:)-flux

            case(SLL_SILVER_MULLER)

               A(1,:) = [(n2*n2+xi*n1*n1)/r,    n2*n1*(xi-1)/r, 0.0_f64, 0.0_f64]
               A(2,:) = [    n2*n1*(xi-1)/r,(n1*n1+xi*n2*n2)/r, 0.0_f64, 0.0_f64]
               A(3,:) = [           0.0_f64,           0.0_f64,       r, 0.0_f64]
               A(4,:) = [           0.0_f64,           0.0_f64, 0.0_f64,    xi*r]

               this%f(left,:) = this%f(left,:)-0.5*matmul(A,this%w(left,:)) 

            case default

               if (this%flux_type == SLL_UNCENTERED) then

                  this%f(left,:) = this%f(left,:) &
                                   -0.5*matmul(A_p,this%w(left,:)) &
                                   -0.5*matmul(A_m,this%r(right,:)) 
               else

                  flux = 0.5*(this%w(left,:)+this%r(right,:))

                  this%f(left,:) = this%f(left,:)-matmul(A,flux)

               end if

            end select

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
   
end subroutine solve_maxwell_2d_diga_parallel

end module sll_maxwell_2d_diga_parallel
