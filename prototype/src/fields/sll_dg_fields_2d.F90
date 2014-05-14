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
   sll_int32                               :: tag
   sll_int32                               :: file_id

contains

   procedure, pass :: write_to_file => write_dg_field_2d_to_file
   procedure, pass :: set_value => initialize_dg_field 

end type dg_field


!interface operator(+)
!  module procedure dg_field_add
!end interface operator(+)

!interface operator(-)
!  module procedure dg_field_sub
!end interface operator(-)

public :: new_dg_field

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
   SLL_CLEAR_ALLOCATE(this%array(1:degree+1,1:degree+1,1:nc_eta1,1:nc_eta2),error)

   if (present(init_function)) then
      call initialize_dg_field( this, init_function, 0.0_f64) 
   end if

   this%tag = 0
   this%file_id = 0

end function new_dg_field

subroutine initialize_dg_field( this, init_function, time) 

   class(dg_field)      :: this
   sll_real64, external :: init_function
   sll_real64           :: time
   sll_real64           :: offset(2)
   sll_real64           :: eta1
   sll_real64           :: eta2
   sll_int32            :: i, j, ii, jj
   
   SLL_ASSERT(associated(this%array))

   do j = 1, this%tau%mesh%num_cells2
   do i = 1, this%tau%mesh%num_cells1
      offset(1) = this%tau%mesh%eta1_min + (i-1)*this%tau%mesh%delta_eta1
      offset(2) = this%tau%mesh%eta2_min + (j-1)*this%tau%mesh%delta_eta2
      do jj = 1, this%degree+1
      do ii = 1, this%degree+1
         eta1 = offset(1) + 0.5 * (this%xgalo(ii) + 1.0) * this%tau%mesh%delta_eta1
         eta2 = offset(2) + 0.5 * (this%xgalo(jj) + 1.0) * this%tau%mesh%delta_eta2
         this%array(ii,jj,i,j) = init_function(this%tau%x1(eta1,eta2), &
                                               this%tau%x2(eta1,eta2), &
                                               time)
      end do
      end do
   end do
   end do

end subroutine initialize_dg_field

subroutine write_dg_field_2d_to_file( this, field_name, file_format )

   class(dg_field)        :: this
   character(len=*)       :: field_name
   sll_int32, optional    :: file_format

   if (present(file_format)) then
      select case(file_format)
      case(SLL_IO_GMSH)
      call plot_dg_field_2d_with_gmsh(this, field_name)
      case(SLL_IO_MTV)
      call plot_dg_field_2d_with_plotmtv(this, field_name)
      case(SLL_IO_XDMF)
      call plot_dg_field_2d_with_xdmf(this, field_name)
      end select
   else
      call plot_dg_field_2d_with_gnuplot( this, field_name )
   endif

end subroutine write_dg_field_2d_to_file

subroutine plot_dg_field_2d_with_gnuplot( this, field_name )

   class(dg_field)        :: this
   character(len=*)       :: field_name
   sll_int32              :: file_id
   sll_int32              :: gnu_id
   sll_real64             :: eta1, eta2
   sll_real64             :: offset(2)
   sll_int32              :: i, j, ii, jj
   sll_int32              :: icell
   character(len=4)       :: ccell
   character(len=4)       :: ctag
   character(len=9)       :: label

   call int2string(this%tag, ctag)

   gnu_id = this%file_id

   if (gnu_id == 0) then
      call sll_ascii_file_create(field_name//".gnu", gnu_id, error)
      rewind(gnu_id)
   else
      open(unit=gnu_id, file=field_name//".gnu", position="append")
   end if
   write(gnu_id,"(a)") "unset key"
   write(gnu_id,"(a)") "set title '"//field_name//" step "//ctag//"'"


   icell = 0
   do j = 1, this%tau%mesh%num_cells2
   do i = 1, this%tau%mesh%num_cells1
 
      icell = icell+1

      call int2string(icell, ccell)
      label = ccell//"-"//ctag

      if (icell == 1) then
         write(gnu_id,"(a)",advance='no') "splot '"//field_name//label//".dat' w l"
      else
         write(gnu_id,"(a)",advance='no') ",'"//field_name//label//".dat' w l "
      end if

      call sll_ascii_file_create(field_name//label//".dat", file_id, error)

      offset(1) = this%tau%mesh%eta1_min + (i-1)*this%tau%mesh%delta_eta1
      offset(2) = this%tau%mesh%eta2_min + (j-1)*this%tau%mesh%delta_eta2
      do jj = 1, this%degree+1
      do ii = 1, this%degree+1
         eta1 = offset(1) + 0.5 * (this%xgalo(ii) + 1.0) * this%tau%mesh%delta_eta1
         eta2 = offset(2) + 0.5 * (this%xgalo(jj) + 1.0) * this%tau%mesh%delta_eta2
         write(file_id,*) this%tau%x1(eta1,eta2), &
                          this%tau%x2(eta1,eta2), &
                          sngl(this%array(ii,jj,i,j))
      end do
      write(file_id,*)
      end do
      close(file_id)

   end do
   end do

   write(gnu_id,*)
   write(gnu_id,'(a)') 'pause -1'
   close(gnu_id)

   this%tag = this%tag+1
   this%file_id = gnu_id
   
end subroutine plot_dg_field_2d_with_gnuplot

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

subroutine plot_dg_field_2d_with_gmsh(this, field_name)

   class(dg_field)        :: this
   character(len=*)       :: field_name
   sll_int32              :: file_id

   sll_int32, parameter :: nlocmax=16
   sll_int32 :: typelem
   sll_int32, dimension(4)  :: invpermut1=(/1,2,4,3/)
   sll_int32, dimension(9)  :: invpermut2=(/1,5,2,8,9,6,4,7,3/)
   sll_int32, dimension(16) :: invpermut3=(/1,5,6,2,12,13,14,7,11,16,15,8,4,10,9,3/)
   sll_int32 :: invpermut(nlocmax),permut(nlocmax)
   sll_int32 :: nloc, nel, neq, iel, ino, inoloc
   sll_int32 :: i, j, k, l, ii, jj, ni, nj
   sll_int32,  allocatable, dimension(:,:) :: connec
   sll_real64, allocatable, dimension(:,:) :: coords
   sll_real64, allocatable, dimension(:)   :: values
   sll_real64 :: offset(2), eta1, eta2
   character(len=32) :: my_fmt

   if (this%degree > 3) then
      write(*,*) 'ordre non prévu'
      stop
   end if

   ni   = this%tau%mesh%num_cells1
   nj   = this%tau%mesh%num_cells2
   nloc = (this%degree+1)*(this%degree+1)
   nel  = ni * nj
   neq  = (ni*this%degree+1)*(nj*this%degree+1)

   SLL_ALLOCATE(connec(1:nloc,1:neq),error)
   SLL_ALLOCATE(coords(1:2,1:neq),error)
   SLL_ALLOCATE(values(1:neq),error)

   do j=0, nj-1
   do i=0, ni-1
      iel=1+i+j*ni
      offset(1) = this%tau%mesh%eta1_min+i*this%tau%mesh%delta_eta1
      offset(2) = this%tau%mesh%eta2_min+j*this%tau%mesh%delta_eta2
      do jj=0,this%degree
      do ii=0,this%degree
         inoloc=1+ii+(this%degree+1)*jj
         ino=1+ii+this%degree*i+(this%degree*ni+1)*(jj+this%degree*j)
         connec(inoloc,iel)=ino
         eta1 = offset(1) + .5*(this%xgalo(ii+1)+1.)*this%tau%mesh%delta_eta1
         eta2 = offset(2) + .5*(this%xgalo(jj+1)+1.)*this%tau%mesh%delta_eta2
         coords(1,ino) = this%tau%x1(eta1,eta2)
         coords(2,ino) = this%tau%x2(eta1,eta2)
         values(ino) = this%array(ii+1,jj+1,i+1,j+1)
      end do
      end do
   end do
   end do

   select case(this%degree)
   case(1)
      invpermut(1:nloc)=invpermut1
      typelem=3
   case(2)
      invpermut(1:nloc)=invpermut2
      typelem=10
   case(3)
      invpermut(1:nloc)=invpermut3
      typelem=36
   case default
      write(*,*) 'ordre non prévu'
      stop
   end select

   do i=1,nloc
      permut(invpermut(i))=i
   end do
      
   call sll_ascii_file_create(field_name//".msh", file_id, error)
   write(file_id,'(A)') '$MeshFormat'
   write(file_id,'(A)') '2 0 8'
   write(file_id,'(A)') '$EndMeshFormat'
   write(file_id,'(A)') '$Nodes'
   write(file_id,*) neq
   do k = 1, neq
      write(file_id,'(i6,2x,3f10.5)') k, coords(1:2,k), 0.0
   end do
   write(file_id,'(A)') '$EndNodes'
   write(file_id,'(A)') '$Elements'
   write(file_id,*) nel

   write(my_fmt, '(a, i0, a)') '(3i6,', nloc, 'i5)'

   k = 0
   do i = 1, this%tau%mesh%num_cells1
   do j = 1, this%tau%mesh%num_cells2
      k = k+1
      write(file_id,trim(my_fmt)) k,typelem,0,(connec(l,k),l=1,nloc)
   end do
   end do
   write(file_id,'(A)') '$EndElements'
   write(file_id,'(A)') '$NodeData'
   write(file_id,*) 1
   write(file_id,*) 'values'
   write(file_id,*) 1
   write(file_id,*) 0
   write(file_id,*) 3
   write(file_id,*) 0
   write(file_id,*) 1
   write(file_id,*) neq
   do k=1,neq
      write(file_id,'(i6,2x,f10.5)') k, values(k)
   end do
   write(file_id,'(A)') '$EndNodeData'

   close(file_id)

end subroutine plot_dg_field_2d_with_gmsh

subroutine plot_dg_field_2d_with_plotmtv(this, field_name)

   class(dg_field)        :: this
   character(len=*)       :: field_name
   sll_int32              :: file_id
   sll_int32              :: ni, nj, ino
   sll_int32              :: i, j, ii, jj
   sll_real64             :: offset(2)
   sll_real64             :: eta1, eta2

   call sll_ascii_file_create(field_name//".mtv", file_id, error)

   write(file_id,*)"$DATA=CONTCURVE"
   write(file_id,*)"%equalscale=T"
   write(file_id,*)"%contfill"
   write(file_id,*)"%toplabel='"//field_name//"'"

   do j=1,this%tau%mesh%num_cells2
      offset(2) = this%tau%mesh%eta2_min+(j-1)*this%tau%mesh%delta_eta2
      do jj=1,merge(this%degree,this%degree+1,j<this%tau%mesh%num_cells2)
         do i=1,this%tau%mesh%num_cells1
            offset(1) = this%tau%mesh%eta1_min+(i-1)*this%tau%mesh%delta_eta1
            do ii=1,merge(this%degree,this%degree+1,i<this%tau%mesh%num_cells1)
               eta1 = offset(1) + .5*(this%xgalo(ii)+1.)*this%tau%mesh%delta_eta1
               eta2 = offset(2) + .5*(this%xgalo(jj)+1.)*this%tau%mesh%delta_eta2
               write(file_id,*) sngl(this%tau%x1(eta1,eta2)), &
                                sngl(this%tau%x2(eta1,eta2)), &
                                sngl(this%array(ii,jj,i,j))
            end do
         end do
      end do
   end do

   write(file_id,*)
   
   write(file_id,*)"$DATA=CURVE3D"
   write(file_id,*)"%equalscale=T"
   write(file_id,*)"%meshplot"
   write(file_id,*)"%toplabel='"//field_name//"'"

   do j=1,this%tau%mesh%num_cells2
   do i=1,this%tau%mesh%num_cells1
      offset(1) = this%tau%mesh%eta1_min+(i-1)*this%tau%mesh%delta_eta1
      offset(2) = this%tau%mesh%eta2_min+(j-1)*this%tau%mesh%delta_eta2
      do jj=1,this%degree
      do ii=1,this%degree
         eta1 = offset(1) + .5*(this%xgalo(ii)+1.)*this%tau%mesh%delta_eta1
         eta2 = offset(2) + .5*(this%xgalo(jj)+1.)*this%tau%mesh%delta_eta2
         write(file_id,*) sngl(this%tau%x1(eta1,eta2)), &
                          sngl(this%tau%x2(eta1,eta2)), &
                          sngl(this%array(ii,jj,i,j))
         eta1 = offset(1) + .5*(this%xgalo(ii+1)+1.)*this%tau%mesh%delta_eta1
         eta2 = offset(2) + .5*(this%xgalo(jj)+1.)*this%tau%mesh%delta_eta2
         write(file_id,*) sngl(this%tau%x1(eta1,eta2)), &
                          sngl(this%tau%x2(eta1,eta2)), &
                          sngl(this%array(ii+1,jj,i,j))
         eta1 = offset(1) + .5*(this%xgalo(ii+1)+1.)*this%tau%mesh%delta_eta1
         eta2 = offset(2) + .5*(this%xgalo(jj+1)+1.)*this%tau%mesh%delta_eta2
         write(file_id,*) sngl(this%tau%x1(eta1,eta2)), &
                          sngl(this%tau%x2(eta1,eta2)), &
                          sngl(this%array(ii+1,jj+1,i,j))
         eta1 = offset(1) + .5*(this%xgalo(ii)+1.)*this%tau%mesh%delta_eta1
         eta2 = offset(2) + .5*(this%xgalo(jj+1)+1.)*this%tau%mesh%delta_eta2
         write(file_id,*) sngl(this%tau%x1(eta1,eta2)), &
                          sngl(this%tau%x2(eta1,eta2)), &
                          sngl(this%array(ii,jj+1,i,j))
         eta1 = offset(1) + .5*(this%xgalo(ii)+1.)*this%tau%mesh%delta_eta1
         eta2 = offset(2) + .5*(this%xgalo(jj)+1.)*this%tau%mesh%delta_eta2
         write(file_id,*) sngl(this%tau%x1(eta1,eta2)), &
                          sngl(this%tau%x2(eta1,eta2)), &
                          sngl(this%array(ii,jj,i,j))
         write(file_id,*)
      end do
      end do
   end do
   end do

   ni   = this%tau%mesh%num_cells1
   nj   = this%tau%mesh%num_cells2

   do j=0, nj-1
   do i=0, ni-1
      offset(1) = this%tau%mesh%eta1_min+i*this%tau%mesh%delta_eta1
      offset(2) = this%tau%mesh%eta2_min+j*this%tau%mesh%delta_eta2
      do jj=0,this%degree
      do ii=0,this%degree
         ino=1+ii+this%degree*i+(this%degree*ni+1)*(jj+this%degree*j)
         eta1 = offset(1)+.5*(this%xgalo(ii+1)+1.)*this%tau%mesh%delta_eta1
         eta2 = offset(2)+.5*(this%xgalo(jj+1)+1.)*this%tau%mesh%delta_eta2
         write(file_id,"(a)"    ,advance="no") "@text x1="
         write(file_id,"(g15.3)",advance="no") this%tau%x1(eta1,eta2)
         write(file_id,"(a)"    ,advance="no")" y1="
         write(file_id,"(g15.3)",advance="no") this%tau%x2(eta1,eta2)
         write(file_id,"(a)"    ,advance="no")" z1=0. lc=5 ll='"
         write(file_id,"(i4)"   ,advance="no") ino
         write(file_id,"(a)")"'"
      end do
      end do
   end do
   end do

   write(file_id,*)"$end"
   close(file_id)
   
end subroutine plot_dg_field_2d_with_plotmtv

subroutine plot_dg_field_2d_with_xdmf(this, field_name)

   class(dg_field)   :: this
   character(len=*)  :: field_name
   sll_int32         :: file_id
   sll_int32         :: i, j, k, ii, jj
   sll_real64        :: offset(2)
   sll_real64        :: eta1, eta2

   SLL_ASSERT(this%degree < 10)

   call sll_ascii_file_create(field_name//".xmf", file_id, error)

   write(file_id,"(a)") "<?xml version='1.0' ?>"
   write(file_id,"(a)") "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>"
   write(file_id,"(a)") "<Xdmf Version='2.'>"
   write(file_id,"(a)") "<Domain>"
   write(file_id,"(a)") "<Grid Name='Mesh' GridType='Collection'>"
   k = 0
   do i=1,this%tau%mesh%num_cells1
   do j=1,this%tau%mesh%num_cells2
      offset(1) = this%tau%mesh%eta1_min+(i-1)*this%tau%mesh%delta_eta1
      offset(2) = this%tau%mesh%eta2_min+(j-1)*this%tau%mesh%delta_eta2
      k = k+1
      write(file_id,"(a,i6,a)") "<Grid Name='Mesh",k,"' GridType='Uniform'>"
      write(file_id,"(a,2i5,a)") &
         "<Topology TopologyType='2DSMesh' NumberOfElements='", & 
         this%degree+1,this%degree+1,"'/>"
      write(file_id,"(a)")"<Geometry GeometryType='X_Y'>"
      write(file_id,"(a,2i3,a)")"<DataItem Dimensions='",this%degree+1, &
                         this%degree+1, "' NumberType='Float' Format='XML'>"
      do jj=1,this%degree+1
      do ii=1,this%degree+1
         eta1 = offset(1)+.5*(this%xgalo(ii)+1.)*this%tau%mesh%delta_eta1
         eta2 = offset(2)+.5*(this%xgalo(jj)+1.)*this%tau%mesh%delta_eta2
         write(file_id,"(f7.3)",advance='no') sngl(this%tau%x1(eta1,eta2))
      end do
      write(file_id,*)
      end do
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a,2i3,a)")"<DataItem Dimensions='",this%degree+1, &
                         this%degree+1, "' NumberType='Float' Format='XML'>"
      do jj=1,this%degree+1
      do ii=1,this%degree+1
         eta1 = offset(1)+.5*(this%xgalo(ii)+1.)*this%tau%mesh%delta_eta1
         eta2 = offset(2)+.5*(this%xgalo(jj)+1.)*this%tau%mesh%delta_eta2
         write(file_id,"(f7.3)", advance='no') sngl(this%tau%x2(eta1,eta2))
      end do
      write(file_id,*)
      end do
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"</Geometry>"
      write(file_id,"(a)")"<Attribute Name='"//field_name// &
           & "' AttributeType='Scalar' Center='Node'>"
      write(file_id,"(a,2i3,a)")"<DataItem Dimensions='", &
           this%degree+1,this%degree+1, &
           "' NumberType='Float' Precision='4' Format='XML'>"
      do jj=1,this%degree+1
      do ii=1,this%degree+1
         eta1 = offset(1)+.5*(this%xgalo(ii)+1.)*this%tau%mesh%delta_eta1
         eta2 = offset(2)+.5*(this%xgalo(jj)+1.)*this%tau%mesh%delta_eta2
         write(file_id,"(f7.3)", advance='no') sngl(this%array(ii,jj,i,j))
      end do
      write(file_id,*)
      end do
      write(file_id,"(a)")"</DataItem>"
      write(file_id,"(a)")"</Attribute>"
      write(file_id,"(a)") "</Grid>"
   end do
   end do
   write(file_id,"(a)") "</Grid>"
   write(file_id,"(a)") "</Domain>"
   write(file_id,"(a)") "</Xdmf>"
   
end subroutine plot_dg_field_2d_with_xdmf

end module sll_dg_fields

