!>Solve Poisson equation on cartesian domain with finit elements.
!> * Compact boundary conditions.
!> * Linear system solve with lapack (Choleski)
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

implicit none

type :: w_vector
   sll_real64, dimension(:), pointer :: ex
   sll_real64, dimension(:), pointer :: ey
   sll_real64, dimension(:), pointer :: hz
end type w_vector

!> Data object to solve Maxwell equations with discontinuous Galerkine
!> method in two dimensions with general coordinates
type :: maxwell_2d_diga
 type(sll_logical_mesh_2d), pointer :: mesh !< Logical mesh
 class(sll_coordinate_transformation_2d_analytic), pointer :: tau  !< Geometric transformation
 sll_int32           :: polarization !< TE or TM
 sll_real64, pointer :: eta1(:,:)
 sll_real64, pointer :: eta2(:,:)
 sll_real64, pointer :: x1(:,:,:)    !< position of gauss points on physical mesh
 sll_real64, pointer :: x2(:,:,:)    !< position of gauss points on physical mesh
 sll_int32           :: degree       !< degree of gauss integration
 type(w_vector), dimension(:), pointer :: W
end type maxwell_2d_diga


interface initialize
   module procedure initialize_maxwell_2d_diga
end interface initialize
interface solve
   module procedure solve_maxwell_2d_diga
end interface solve

sll_int32, private :: nx, ny
sll_int32, private :: i, j, k, l, ii, jj, kk, ll
sll_int32, private :: error

contains

!> Initialize Poisson solver object using finite elements method.
!> Indices are shifted from [1:n+1] to [0:n] only inside this 
!> subroutine
subroutine initialize_maxwell_2d_diga( this, tau, degree, polarization)
type( maxwell_2d_diga ) :: this !< solver data object
class(sll_coordinate_transformation_2d_analytic), pointer :: tau
sll_int32  :: polarization
sll_int32  :: degree
sll_int32  :: nddl
sll_int32  :: nquads
sll_real64 :: xgalo(degree+1)
sll_real64 :: wgalo(degree+1)
sll_real64 :: dlag(degree+1,degree+1)
sll_real64 :: xk(degree+1)

this%tau  => tau
this%degree = degree
this%mesh => tau%mesh
this%polarization = polarization

call tau%write_to_file(SLL_IO_MTV)

nddl   = (degree+1)*(degree+1)
nquads = this%mesh%num_cells1*this%mesh%num_cells2
xgalo  = gauss_lobatto_points(degree+1,-1._f64,1._f64)
wgalo  = gauss_lobatto_weights(degree+1)
dlag   = gauss_lobatto_derivative_matrix(degree+1, -1._f64, 1._f64) 

SLL_CLEAR_ALLOCATE(this%eta1(1:degree+1,1:nquads),  error)
SLL_CLEAR_ALLOCATE(this%eta2(1:degree+1,1:nquads),  error)
SLL_CLEAR_ALLOCATE(this%x1(1:nddl,1:nddl,1:nquads), error)
SLL_CLEAR_ALLOCATE(this%x2(1:nddl,1:nddl,1:nquads), error)

k = 0
do j = 1, this%mesh%num_cells2
   do i = 1, this%mesh%num_cells1
      k = k+1
      this%eta1(:,k) = this%mesh%eta1_min + (i-0.5+0.5*xgalo(k))*this%mesh%delta_eta1
      this%eta2(:,k) = this%mesh%eta2_min + (j-0.5+0.5*xgalo(k))*this%mesh%delta_eta2
   end do
end do

k = 0
do j = 1, this%mesh%num_cells2
do i = 1, this%mesh%num_cells1
   k = k + 1
   do ii = 1, degree+1
   do jj = 1, degree+1
      this%x1(ii,jj,k) = tau%x1(this%eta1(ii,k),this%eta2(jj,k))
      this%x2(ii,jj,k) = tau%x2(this%eta1(ii,k),this%eta2(jj,k))
   end do
   end do
end do
end do

SLL_ALLOCATE(this%W(nquads),error)
do k = 1, nquads
   SLL_ALLOCATE(this%W(k)%ex(nddl),error)
   SLL_ALLOCATE(this%W(k)%ey(nddl),error)
   SLL_ALLOCATE(this%W(k)%hz(nddl),error)
end do

end subroutine initialize_maxwell_2d_diga

!need to be replace by a macro
integer function som(i, j, k)

   integer :: i, j, k

   if (k == 1) then
      som = i-1+(j-2)*(nx-1)
   else if (k == 2) then
      som = i-1+(j-2)*(nx-1)+1
   else if (k == 3) then
      som = i-1+(j-2)*(nx-1)+nx
   else if (k == 4) then
      som = i-1+(j-2)*(nx-1)+nx-1
   end if 

end function som

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

!> Construction of the derivative matrix for Gauss-Lobatto 2D
!> 
!>          der(i,j)=int(Phi_{i,j}.Phi_{k,l})_[-1;1]²
!>                  =w_i.Phi'_j(x_i)

subroutine mass_matrix(degree, x, w)

 sll_int32  :: degree
 sll_real64 :: x(degree+1)
 sll_real64 :: w(degree+1)
 sll_real64 :: prod

 nb_pts=gl_obj%degree+1

 gl_obj%der=0.0d0

    !loop on all element of D
    !loop on columns
    do j=1,nb_pts
       !loop on rows
       do i=1,nb_pts
          !loop on all the derivatives
          !the code is writen so there is no if
          do l=1,j-1
             prod=1.0d0
             do m=1,l-1!min(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=l+1,j-1!min(j,l)+1,max(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=j+1,nb_pts!max(j,l)+1,nb_pts
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             prod=prod/(gl_obj%node(j)-gl_obj%node(l))
             gl_obj%der(i,j)=gl_obj%der(i,j)+prod
          end do
          do l=j+1,nb_pts
             prod=1.0d0
             do m=1,j-1!min(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=j+1,l-1!min(j,l)+1,max(j,l)-1
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             do m=l+1,nb_pts!max(j,l)+1,nb_pts
                prod=prod*(gl_obj%node(i)-gl_obj%node(m))/(gl_obj%node(j)-gl_obj%node(m))
             end do
             prod=prod/(gl_obj%node(j)-gl_obj%node(l))
             gl_obj%der(i,j)=gl_obj%der(i,j)+prod
          end do
          gl_obj%der(i,j)=gl_obj%der(i,j)*gl_obj%weigh(i)
       end do
    end do

  end subroutine derivative_matrix_1d

end module sll_gausslobatto

end module sll_maxwell_2d_diga
