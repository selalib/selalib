#define sll_transformation class(sll_c_coordinate_transformation_2d_base)

!> @ingroup maxwell_solvers
!> @brief DG for Maxwell
!> @details
!> Solve Maxwell equations on cartesian domain with Discontinuous Galerkine method:
!> * Gauss Lobatto for integration formula
!> * Periodic boundary conditions.
module sll_m_maxwell_2d_diga

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_conductor, &
    sll_p_interior, &
    sll_p_silver_muller

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_dg_fields, only: &
    sll_t_dg_field_2d, &
    sll_o_new

  use sll_m_gauss_lobatto_integration, only: &
    sll_f_gauss_lobatto_derivative_matrix, &
    sll_f_gauss_lobatto_points, &
    sll_f_gauss_lobatto_weights

  implicit none

  public :: &
    sll_o_create, &
    sll_t_maxwell_2d_diga, &
    sll_o_new, &
    sll_o_solve, &
    sll_p_uncentered

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Local type with edge properties
type :: edge_type

   sll_real64                            :: length 
   sll_real64, dimension(:,:), pointer   :: vec_norm
   sll_int32                             :: bc_type

end type edge_type

!> Information about a mesh cell
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

!> DG method in 2D with general coordinates
type :: sll_t_maxwell_2d_diga
   private
   sll_int32                           :: nc_eta1      !< x cells number
   sll_int32                           :: nc_eta2      !< y cells number
   sll_int32                           :: polarization !< TE or TM
   sll_real64                          :: e_0          !< electric conductivity
   sll_real64                          :: mu_0         !< magnetic permeability
   sll_real64                          :: c            !< speed of light
   sll_real64                          :: eta1_min     !< left side 
   sll_real64                          :: eta1_max     !< right side
   sll_real64                          :: delta_eta1   !< step size
   sll_real64                          :: eta2_min     !< bottom side
   sll_real64                          :: eta2_max     !< top side
   sll_real64                          :: delta_eta2   !< step size
   sll_transformation, pointer         :: tau          !< transformation
   type(sll_t_cartesian_mesh_2d), pointer  :: mesh         !< Logical mesh
   sll_int32                           :: degree       !< degree of gauss integration
   type(cell_type), pointer            :: cell(:,:)    !< mesh cells
   sll_real64, pointer                 :: f(:,:)       !< cell flux
   sll_real64, pointer                 :: w(:,:)       !< edge flux
   sll_real64, pointer                 :: r(:,:)       !< source flux
   sll_int32                           :: bc_south
   sll_int32                           :: bc_east
   sll_int32                           :: bc_north
   sll_int32                           :: bc_west
   sll_int32                           :: flux_type
   type(sll_t_dg_field_2d), pointer      :: po           !< Potential
   sll_real64                          :: xi 

end type sll_t_maxwell_2d_diga

!> Create a Maxwell solver using DG method in 2D
interface sll_o_new
   module procedure new_maxwell_2d_digal
end interface sll_o_new

!> Create a Maxwell solver object using Discontinuous Galerkine 
interface sll_o_create
   module procedure initialize_maxwell_2d_diga
end interface sll_o_create

!> Solve Maxwell system
interface sll_o_solve
   module procedure solve_maxwell_2d_diga
end interface sll_o_solve

!> Flux parameter
sll_int32, parameter :: SLL_CENTERED       = 20
!> Flux parameter
sll_int32, parameter :: sll_p_uncentered     = 21



contains

function new_maxwell_2d_digal( tau,          &
                               degree,       &
                               polarization, &
                               bc_south,     &
                               bc_east,      &
                               bc_north,     &
                               bc_west,      &
                               flux_type) result(self)

   type( sll_t_maxwell_2d_diga ), pointer :: self !< solver data object
   sll_transformation, pointer      :: tau
   sll_int32                        :: polarization
   sll_int32                        :: degree
   sll_int32, intent(in)            :: bc_east
   sll_int32, intent(in)            :: bc_west
   sll_int32, intent(in)            :: bc_north
   sll_int32, intent(in)            :: bc_south
   sll_int32, optional              :: flux_type

   sll_int32                        :: error

   SLL_ALLOCATE(self,error)

   call initialize_maxwell_2d_diga( self,         &
                                    tau,          &
                                    degree,       &
                                    polarization, &
                                    bc_south,     &
                                    bc_east,      &
                                    bc_north,     &
                                    bc_west,      &
                                    flux_type)

 end function new_maxwell_2d_digal

!> Initialize Maxwell solver object using DG method.
subroutine initialize_maxwell_2d_diga( self,         &
                                       tau,          &
                                       degree,       &
                                       polarization, &
                                       bc_south,     &
                                       bc_east,      &
                                       bc_north,     &
                                       bc_west,      &
                                       flux_type)

   type(sll_t_maxwell_2d_diga)       :: self !< solver data object
   sll_transformation, pointer :: tau  !< transformation
   sll_int32                   :: polarization !< TE or TM
   sll_int32                   :: degree !< degree of DG method
   sll_int32, intent(in)       :: bc_east !< Boundary condition
   sll_int32, intent(in)       :: bc_west !< Boundary condition
   sll_int32, intent(in)       :: bc_north !< Boundary condition
   sll_int32, intent(in)       :: bc_south !< Boundary condition
   sll_int32, optional         :: flux_type !< centered or not

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

   sll_int32  :: error

   self%tau        => tau
   ! Please undo self 'fix' whenever it is decided that gfortran 4.6 is no
   ! longer supported.
   !   self%mesh       => tau%get_cartesian_mesh()
   self%mesh       => tau%mesh
   self%bc_south   =  bc_south
   self%bc_east    =  bc_east
   self%bc_north   =  bc_north
   self%bc_west    =  bc_west

   self%nc_eta1    = self%mesh%num_cells1
   self%nc_eta2    = self%mesh%num_cells2
   self%eta1_min   = self%mesh%eta1_min
   self%eta2_min   = self%mesh%eta2_min
   self%eta1_max   = self%mesh%eta1_max
   self%eta2_max   = self%mesh%eta2_max
   self%delta_eta1 = self%mesh%delta_eta1
   self%delta_eta2 = self%mesh%delta_eta2

   self%xi           = 0.0_f64
   self%degree       = degree
   self%polarization = polarization

   if (present(flux_type)) then
      self%flux_type = flux_type
   else
      self%flux_type = SLL_CENTERED
   end if

   nddl   = (degree+1)*(degree+1)
   ncells = self%nc_eta1*self%nc_eta2

   SLL_ALLOCATE(self%cell(self%nc_eta1,self%nc_eta2), error)

   x    = sll_f_gauss_lobatto_points(degree+1,0.0_f64,1.0_f64)
   w    = sll_f_gauss_lobatto_weights(degree+1,0.0_f64,1.0_f64)
   dlag = sll_f_gauss_lobatto_derivative_matrix(degree+1,x)

   dtau_ij_mat(1,1) = self%mesh%delta_eta1
   dtau_ij_mat(1,2) = 0.0_f64
   dtau_ij_mat(2,1) = 0.0_f64
   dtau_ij_mat(2,2) = self%mesh%delta_eta2

   do j = 1, self%nc_eta2   !Loop over cells
   do i = 1, self%nc_eta1

      call compute_normals(tau,bc_south,bc_east,bc_north,bc_west, &
                           i,j,degree,self%cell(i,j))

      allocate(self%cell(i,j)%MassMatrix(nddl))
      self%cell(i,j)%MassMatrix(nddl) = 0.0_f64
      allocate(self%cell(i,j)%DxMatrix(nddl, nddl))
      self%cell(i,j)%DxMatrix = 0.0_f64
      allocate(self%cell(i,j)%DyMatrix(nddl, nddl))
      self%cell(i,j)%DyMatrix = 0.0_f64

      xa = self%cell(i,j)%eta1_min ; xb  = self%cell(i,j)%eta1_max 
      ya = self%cell(i,j)%eta2_min ; yb  = self%cell(i,j)%eta2_max 

      do jj = 1, degree+1
      do ii = 1, degree+1

         j_mat = tau%jacobian_matrix(xa+x(ii)*self%delta_eta1,&
                                     ya+x(jj)*self%delta_eta2)
         j_mat(:,1) = j_mat(:,1) * self%delta_eta1
         j_mat(:,2) = j_mat(:,2) * self%delta_eta2

         det   = (j_mat(1,1)*j_mat(2,2)-j_mat(1,2)*j_mat(2,1))

         self%cell(i,j)%MassMatrix((ii-1)*(degree+1)+jj) = w(ii)*w(jj)*det

      end do
      end do

      do jj = 1, degree+1
      do ii = 1, degree+1

         do ll = 1, degree+1
         do kk = 1, degree+1


            j_mat = tau%jacobian_matrix(xa+x(kk)*self%delta_eta1,&
                                        ya+x(ll)*self%delta_eta2)
            j_mat(:,1) = j_mat(:,1) * self%delta_eta1
            j_mat(:,2) = j_mat(:,2) * self%delta_eta2

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

            self%cell(i,j)%DxMatrix(k,l) = & 
               det*w(kk)*w(ll)*(inv_j(1,1)*dfx + inv_j(2,1)*dfy)

            self%cell(i,j)%DyMatrix(k,l) = &
               det*w(kk)*w(ll)*(inv_j(1,2)*dfx + inv_j(2,2)*dfy)

         end do
         end do

      end do
      end do

   end do
   end do

   !call sll_o_display(self%cell(1,1)%MassMatrix,"f9.4")
   !call sll_o_display(self%cell(1,1)%DxMatrix,"f9.4")
   !call sll_o_display(self%cell(1,1)%DyMatrix,"f9.4")

   SLL_CLEAR_ALLOCATE(self%w((degree+1)*(degree+1),4),error)
   SLL_CLEAR_ALLOCATE(self%r((degree+1)*(degree+1),4),error)
   SLL_CLEAR_ALLOCATE(self%f((degree+1)*(degree+1),4),error)

   self%po => sll_o_new( degree, tau) 

end subroutine initialize_maxwell_2d_diga


!> Solve the maxwell equation
subroutine solve_maxwell_2d_diga( self, fx, fy, fz, dx, dy, dz )

   type( sll_t_maxwell_2d_diga )  :: self !< Maxwell solver object

   type(sll_t_dg_field_2d)  :: fx   !< x electric field
   type(sll_t_dg_field_2d)  :: fy   !< y electric field
   type(sll_t_dg_field_2d)  :: fz   !< z magnetic field

   type(sll_t_dg_field_2d)  :: dx   !< x size step
   type(sll_t_dg_field_2d)  :: dy   !< y size step
   type(sll_t_dg_field_2d)  :: dz   !< z size step

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

   xi = self%xi

   do i = 1, self%nc_eta1
   do j = 1, self%nc_eta2

      do jj = 1, self%degree+1
      do ii = 1, self%degree+1
         kk = (ii-1)*(self%degree+1)+jj
         self%w(kk,1) = fx%array(ii,jj,i,j)
         self%w(kk,2) = fy%array(ii,jj,i,j)
         self%w(kk,3) = fz%array(ii,jj,i,j)
         self%w(kk,4) = self%po%array(ii,jj,i,j)
      end do
      end do
         
      self%f(:,1) = - matmul(self%cell(i,j)%DyMatrix,self%w(:,3)) &
                 + xi*matmul(self%cell(i,j)%DxMatrix,self%w(:,4))
      self%f(:,2) =   matmul(self%cell(i,j)%DxMatrix,self%w(:,3)) &
                 + xi*matmul(self%cell(i,j)%DyMatrix,self%w(:,4))
      self%f(:,3) = - matmul(self%cell(i,j)%DyMatrix,self%w(:,1)) &
                    + matmul(self%cell(i,j)%DxMatrix,self%w(:,2))
      self%f(:,4) = + matmul(self%cell(i,j)%DxMatrix,self%w(:,1)) &
                    + matmul(self%cell(i,j)%DyMatrix,self%w(:,2))


#ifdef VERBOSE

      print*,"####################################################"
      print*,' (i,j) ', i, j
      print"('Ex=',9f7.3)", self%w(:,1)
      print"('Ey=',9f7.3)", self%w(:,2)
      print"('Bz=',9f7.3)", self%w(:,3)

      print"('dEx=',9f7.3)", self%f(:,1)
      print"('dEy=',9f7.3)", self%f(:,2)
      print"('dBz=',9f7.3)", self%f(:,3)

      print"(a)", 'side'//'  bc '//'left'//'    w(left) ' &
                 &//' f(left) '//'  n1 '//'       n2 '// '       flux '
#endif

      do side = 1, 4 ! Loop over each side of the cell
 
         bc_type = self%cell(i,j)%edge(side)%bc_type
         flux_type = self%flux_type
         !periodic boundary conditions
         select case(side)
         case(SOUTH)
            k = i
            l = 1+modulo(j-2,self%nc_eta2) 
         case(EAST)
            k = 1+modulo(i  ,self%nc_eta1)
            l = j
         case(NORTH)
            k = i
            l = 1+modulo(j  ,self%nc_eta2)
         case(WEST)
            k = 1+modulo(i-2,self%nc_eta1)
            l = j
         end select

         do jj = 1, self%degree+1
         do ii = 1, self%degree+1
            kk = (ii-1)*(self%degree+1)+jj
            self%r(kk,1) = fx%array(ii,jj,k,l)
            self%r(kk,2) = fy%array(ii,jj,k,l)
            self%r(kk,3) = fz%array(ii,jj,k,l)
            self%r(kk,4) = self%po%array(ii,jj,k,l)
         end do
         end do
   
#ifdef VERBOSE
         print*,'--'
#endif
         do node = 1, self%degree+1
   
            left   = dof_local(side, node, self%degree)
            right  = dof_neighbor(side, node, self%degree)

            n1 = self%cell(i,j)%edge(side)%vec_norm(node,1)
            n2 = self%cell(i,j)%edge(side)%vec_norm(node,2)
            r  = sqrt(n1*n1+n2*n2)

            bc_type = self%cell(i,j)%edge(side)%bc_type

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

            case(sll_p_conductor)

               A(1,:) = [-1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64]
               A(2,:) = [ 0.0_f64,-1.0_f64, 0.0_f64, 0.0_f64]
               A(3,:) = [ 0.0_f64, 0.0_f64, 1.0_f64, 0.0_f64]
               A(4,:) = [-2*xi*n1,-2*xi*n2, 0.0_f64,-1.0_f64]

               flux = matmul(A,self%w(left,:))
               flux = matmul(A_m, flux) + matmul(A_p,self%w(left,:))

               self%f(left,:) = self%f(left,:)-flux

            case(sll_p_silver_muller)

               A(1,:) = [(n2*n2+xi*n1*n1)/r,    n2*n1*(xi-1)/r, 0.0_f64, 0.0_f64]
               A(2,:) = [    n2*n1*(xi-1)/r,(n1*n1+xi*n2*n2)/r, 0.0_f64, 0.0_f64]
               A(3,:) = [           0.0_f64,           0.0_f64,       r, 0.0_f64]
               A(4,:) = [           0.0_f64,           0.0_f64, 0.0_f64,    xi*r]

               self%f(left,:) = self%f(left,:)-0.5*matmul(A,self%w(left,:)) 

            case default

               if (self%flux_type == sll_p_uncentered) then

                  self%f(left,:) = self%f(left,:) &
                                   -0.5*matmul(A_p,self%w(left,:)) &
                                   -0.5*matmul(A_m,self%r(right,:)) 
               else

                  flux = 0.5*(self%w(left,:)+self%r(right,:))

                  self%f(left,:) = self%f(left,:)-matmul(A,flux)

               end if

            end select

#ifdef VERBOSE
            print"(3i4,5f10.4)",side, bc_type, left, &
                               self%w(left,3), self%r(right,3), &
                               n1, n2, -n1*self%w(left,3)
#endif
         end do

      end do

#ifdef VERBOSE
      print"('fEx=',9f7.3)", self%f(:,1)
      print"('fEy=',9f7.3)", self%f(:,2)
      print"('fBz=',9f7.3)", self%f(:,3)
#endif

      kk = 0
      do jj = 1, self%degree+1
      do ii = 1, self%degree+1
         kk = (ii-1)*(self%degree+1)+jj
         dx%array(ii,jj,i,j) = self%f(kk,1)/self%cell(i,j)%MassMatrix(kk)
         dy%array(ii,jj,i,j) = self%f(kk,2)/self%cell(i,j)%MassMatrix(kk)
         dz%array(ii,jj,i,j) = self%f(kk,3)/self%cell(i,j)%MassMatrix(kk)
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

   sll_transformation, pointer :: tau
   type(cell_type)             :: cell
   sll_int32                   :: i, j, d
   sll_real64                  :: x(d+1)
   sll_real64                  :: w(d+1)
   sll_real64                  :: vec_norm(d+1,2)
   sll_real64                  :: a, b, c1, c2
   sll_real64                  :: xk, wk
   sll_real64                  :: jac_mat(2,2)
   sll_real64                  :: co_jac_mat(2,2)
   sll_real64                  :: dtau_ij_mat(2,2)
   sll_real64                  :: jac_mat_sll(2,2)
   sll_real64                  :: length
   sll_int32                   :: side
   sll_int32                   :: bc_south
   sll_int32                   :: bc_east
   sll_int32                   :: bc_north
   sll_int32                   :: bc_west
   sll_int32                   :: k
   sll_int32                   :: error
   class(sll_t_cartesian_mesh_2d), pointer :: lm

   lm => tau%get_cartesian_mesh()
   
   cell%i = i
   cell%j = j
   cell%eta1_min = lm%eta1_min + (i-1)*lm%delta_eta1
   cell%eta2_min = lm%eta2_min + (j-1)*lm%delta_eta2
   
   cell%eta1_max = cell%eta1_min + lm%delta_eta1
   cell%eta2_max = cell%eta2_min + lm%delta_eta2

   dtau_ij_mat(1,1) = lm%delta_eta1
   dtau_ij_mat(1,2) = 0.0_f64
   dtau_ij_mat(2,1) = 0.0_f64
   dtau_ij_mat(2,2) = lm%delta_eta2

   do side = 1, 4
      SLL_CLEAR_ALLOCATE(cell%edge(side)%vec_norm(1:d+1,1:2),error)
      cell%edge(side)%bc_type = sll_p_interior
   end do

   if (j ==                   1) cell%edge(SOUTH)%bc_type = bc_south
   if (i == lm%num_cells1) cell%edge(EAST)%bc_type  = bc_east
   if (j == lm%num_cells2) cell%edge(NORTH)%bc_type = bc_north
   if (i ==                   1) cell%edge(WEST)%bc_type  = bc_west
   
   x = sll_f_gauss_lobatto_points(d+1)
   w = sll_f_gauss_lobatto_weights(d+1)

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
   cell%edge(SOUTH)%vec_norm = vec_norm/lm%delta_eta1
   
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
   cell%edge(EAST)%vec_norm = vec_norm/lm%delta_eta2
   
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
   cell%edge(NORTH)%vec_norm = vec_norm/lm%delta_eta1
   
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
   cell%edge(WEST)%vec_norm = vec_norm/lm%delta_eta2

end subroutine compute_normals

end module sll_m_maxwell_2d_diga
