#define sll_transformation class(sll_coordinate_transformation_2d_analytic)

!> Solve Maxwell equations on cartesian domain with Disconituous Galerkine method:
!> * Gauss Lobatto for integration formula
!> * Periodic boundary conditions.
module sll_dg_fields

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_file_io.h"
#include "sll_integration.h"
#include "sll_utilities.h"
#include "sll_assert.h"

use sll_logical_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations

implicit none
private

type, public :: dg_field

   sll_int32                               :: degree
   sll_transformation, pointer             :: tau  
   sll_real64, dimension(:,:,:,:), pointer :: array
   sll_real64, dimension(:), pointer       :: xgalo
   sll_real64, dimension(:), pointer       :: wgalo

end type dg_field

!interface operator(+)
!  module procedure dg_field_add
!end interface operator(+)

!interface operator(-)
!  module procedure dg_field_sub
!end interface operator(-)

public :: new_dg_field, plot_dg_field, operator(-)

sll_int32, private :: error

contains

function new_dg_field( degree, tau, init_function ) result (this) 

   sll_transformation, pointer    :: tau           !< transformation 
   sll_real64, external, optional :: init_function !< function
   sll_int32, intent(in)          :: degree        !< degree integration
   sll_int32                      :: nc_eta1
   sll_int32                      :: nc_eta2
   type(dg_field), pointer        :: this
   sll_int32                      :: error

   SLL_ALLOCATE(this, error)
   this%tau    => tau
   this%degree =  degree

   SLL_ALLOCATE(this%xgalo(degree+1),error)
   SLL_ALLOCATE(this%wgalo(degree+1),error)

   this%xgalo  = gauss_lobatto_points(degree+1,-1._f64,1._f64)
   this%wgalo  = gauss_lobatto_weights(degree+1)
   nc_eta1 = tau%mesh%num_cells1
   nc_eta2 = tau%mesh%num_cells2
   SLL_CLEAR_ALLOCATE(this%array(1:nc_eta1,1:nc_eta2,1:degree+1,1:degree+1),error)

   if (present(init_function)) then
      call initialize_dg_field( this, init_function, 0.0_f64) 
   end if

end function new_dg_field

subroutine initialize_dg_field( this, init_function, time) 

   type(dg_field)          :: this
   sll_real64, external    :: init_function
   sll_real64              :: time
   sll_real64              :: offset(2)
   sll_real64              :: eta1
   sll_real64              :: eta2
   sll_int32               :: i, j, ii, jj
   
   SLL_ASSERT(associated(this%array))

   do i = 1, this%tau%mesh%num_cells1
   do j = 1, this%tau%mesh%num_cells2
      offset(1) = this%tau%mesh%eta1_min + (i-1)*this%tau%mesh%delta_eta1
      offset(2) = this%tau%mesh%eta2_min + (j-1)*this%tau%mesh%delta_eta2
      do ii = 1, this%degree+1
      do jj = 1, this%degree+1
         eta1 = offset(1) + 0.5 * (this%xgalo(ii) + 1.0) * this%tau%mesh%delta_eta1
         eta2 = offset(2) + 0.5 * (this%xgalo(jj) + 1.0) * this%tau%mesh%delta_eta2
         this%array(i,j,ii,jj) = init_function(this%tau%x1(eta1,eta2), &
                                               this%tau%x2(eta1,eta2), &
                                               time)
      end do
      end do
   end do
   end do

end subroutine initialize_dg_field

subroutine plot_dg_field( this, field_name )

   type(dg_field)         :: this
   character(len=*)       :: field_name
   sll_int32              :: file_id
   sll_int32              :: gnu_id
   sll_real64             :: eta1, eta2
   sll_real64             :: offset(2)
   sll_int32              :: i, j, ii, jj
   sll_int32              :: icell
   character(len=4)       :: ccell

   call sll_ascii_file_create(field_name//".gnu", gnu_id, error)

   icell = 0
   do i = 1, this%tau%mesh%num_cells1
   do j = 1, this%tau%mesh%num_cells2
 
      icell = icell+1

      call int2string(icell, ccell)

      if (icell == 1) then
         write(gnu_id,"(a)",advance='no') "splot '"//field_name//ccell//".dat' w l"
      else
         write(gnu_id,"(a)",advance='no') ",'"//field_name//ccell//".dat' w l "
      end if

      call sll_ascii_file_create(field_name//ccell//".dat", file_id, error)

      offset(1) = this%tau%mesh%eta1_min + (i-1)*this%tau%mesh%delta_eta1
      offset(2) = this%tau%mesh%eta2_min + (j-1)*this%tau%mesh%delta_eta2
      do ii = 1, this%degree+1
      do jj = 1, this%degree+1
         eta1 = offset(1) + 0.5 * (this%xgalo(ii) + 1.0) * this%tau%mesh%delta_eta1
         eta2 = offset(2) + 0.5 * (this%xgalo(jj) + 1.0) * this%tau%mesh%delta_eta2
         write(file_id,*) this%tau%x1(eta1,eta2), &
                          this%tau%x2(eta1,eta2), &
                          sngl(this%array(i,j,ii,jj))
      end do
      write(file_id,*)
      end do
      close(file_id)

   end do
   end do

   write(gnu_id,*)
   close(gnu_id)
   
end subroutine plot_dg_field


!function dg_field_add( W1, W2) result(W3)
!
!  type(dg_field), intent(in) :: W1
!  type(dg_field), intent(in) :: W2
!  type(dg_field)             :: W3
!
!  SLL_ASSERT(W1%degree == W2%degree)
!  SLL_ASSERT(associated(W1%array))
!  SLL_ASSERT(associated(W2%array))
!
!  W3%array  = W1%array + W2%array
!
!end function dg_field_add

!function dg_field_sub( W1, W2) result(W3)
!
!  type(dg_field), intent(in) :: W1
!  type(dg_field), intent(in) :: W2
!  type(dg_field)             :: W3
!
!  SLL_ASSERT(W1%degree == W2%degree)
!  SLL_ASSERT(associated(W1%array))
!  SLL_ASSERT(associated(W2%array))
!
!  W3%array  = W1%array - W2%array
!
!end function dg_field_sub

end module sll_dg_fields

