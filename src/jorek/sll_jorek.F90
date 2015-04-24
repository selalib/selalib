!> @ingroup <DIRECTORY_NAME>
!> @author <MODULE_AUTHOR_NAME_AND_AFFILIATION>
!> @brief <BRIEF_DESCRIPTION>
!> @details <DETAILED_DESCRIPTION>
module sll_jorek

#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_file_io.h"
#include "sll_memory.h"

use typedef
use jorek_model
use jorek_model
use jorek_interface
use linear_solver
use field
use mesh
use spm_def
use matrix
use space
use fem

implicit none

type, public :: sll_jorek_solver

  type(def_fem_quad_2d), pointer      :: fem_model 
  sll_int32                           :: n_quads
  sll_int32                           :: n_nodes
  sll_real64, dimension(:,:), pointer :: coord
  sll_int32,  dimension(:,:), pointer :: nodes
  sll_real64, dimension(:),   pointer :: phi
  sll_real64, dimension(:),   pointer :: rho
  sll_real64, dimension(:),   pointer :: e_x
  sll_real64, dimension(:),   pointer :: e_y
  logical                             :: created = .false.
  logical                             :: deleted = .true.

end type sll_jorek_solver

interface sll_create
module procedure initialize_jorek
end interface sll_create

interface sll_solve
module procedure solve_jorek
end interface sll_solve

interface sll_delete
module procedure delete_jorek
end interface sll_delete

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------
!> @brief Create your jorek model to solve Poisson in 2d
!> @details 
!> When you call the subroutine create_jorek_model, the program
!> command line arguments you added.
!> --solver solver.txt
!> --parameters parameter.txt
!> --geometry square_periodic/n16x16
!> Different geometries are available in jorek-prefix/src/jorek/geometries 
!> directory located in your cmake build directory.
!> @param[OUT] jorek the solver derived that contains data useful for selalib
subroutine initialize_jorek(jorek)

type(sll_jorek_solver)     :: jorek
type(def_element), pointer :: elmt => null()
sll_int32                  :: iq
sll_int32                  :: ierr

call create_jorek_model()

jorek%fem_model => fem_model
jorek%n_nodes = fem_mesh%oi_n_nodes
jorek%n_quads = fem_mesh%oi_n_elmts

SLL_ALLOCATE(jorek%coord(2,jorek%n_nodes),ierr)
SLL_ALLOCATE(jorek%nodes(4,jorek%n_quads),ierr)

SLL_ALLOCATE(jorek%rho(jorek%n_nodes),ierr)
SLL_ALLOCATE(jorek%phi(jorek%n_nodes),ierr)
SLL_ALLOCATE(jorek%e_x(jorek%n_nodes),ierr)
SLL_ALLOCATE(jorek%e_y(jorek%n_nodes),ierr)

jorek%coord = fem_mesh%thenodes%coor2d(1:2, one, :)
do iq = 1, jorek%n_quads
  elmt => fem_mesh%opo_elements(iq)
  jorek%nodes(:,iq) = elmt%opi_vertices(1:elmt%oi_n_vtex)
end do

jorek%created = .true.  
jorek%deleted = .false.  

rho => jorek%rho
e_x => jorek%e_x
e_y => jorek%e_y

end subroutine initialize_jorek

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------
!> @brief Run the created jorek model
!> @details This subroutine run the jorek model and solve the Poisson
!> equation in 2D
!> @param[INOUT] jorek the jorek data type
subroutine solve_jorek(jorek)

type(sll_jorek_solver) :: jorek

call run_jorek_model()

  jorek%phi = fem_model%opr_global_var(one,:)


end subroutine solve_jorek

!---------------------------------------------------------------------------
!> @brief Delete the jorek object and the jorek model
!> @details free the memory allocated by jorek.
!> @param[INOUT>] jorek the jorek derived type linked to the jorek model.
subroutine delete_jorek(jorek)

type(sll_jorek_solver) :: jorek

call delete_jorek_model()
jorek%created = .false.  
jorek%deleted = .true.  

end subroutine delete_jorek


!---------------------------------------------------------------------------
!> @brief Write potential computed by Jorek in xdmf format file.
!> @details This file is readable by VisIt (LLNL)
!> @param[IN] jorek the solver object that conatins all data
!> @param[IN] filename prefix of the xdmf file
!> @param[IN] label field name (string)
subroutine plot_field(jorek, filename, label)

type(sll_jorek_solver)       :: jorek
sll_int32                    :: is
sll_int32                    :: iq
sll_int32, parameter         :: n_dim=2
character(len=*), intent(in) :: filename
character(len=*), intent(in) :: label
integer                      :: nbs
integer                      :: nbq
integer                      :: xmf
integer                      :: j
sll_int32                    :: ierr

call sll_new_file_id(xmf, ierr)

nbs = jorek%n_nodes
nbq = jorek%n_quads

open(xmf, file=filename//".xmf")

write(xmf,"(a)") "<?xml version='1.0'?>"
write(xmf,"(a)") "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []> "
write(xmf,"(a)") "<Xdmf Version='2.0'>"
write(xmf,"(a)") "<Domain>"
write(xmf,"(a)") "<Grid Name='Mesh' GridType='Uniform' >"
write(xmf,"(a,i6,a)") "<Topology Type='Quadrilateral' NumberOfElements='",nbq,"'>"
write(xmf,"(a,i6,a)") "<DataItem Name='Connections' Format='XML' DataType='Int' Dimensions='", &
                      nbq, " 4'>" 
do iq=1, nbq
  write(xmf,"(4i6)") (jorek%nodes(j,iq)-1,j=1,4)
end do
write(xmf,"(a)") "</DataItem>"
write(xmf,"(a)") "</Topology>"
write(xmf,"(a)") "<Geometry GeometryType='XY'>"
write(xmf,"(a,i6,a)") "<DataItem Format='XML' Dimensions='", nbs, " 2'>"
do is=1, nbs
  write(xmf,"(2f12.6)") jorek%coord(1:2,is)
end do
write(xmf,"(a)") "</DataItem>"
write(xmf,"(a)") "</Geometry>"
write(xmf,"(a)") "<Attribute Name='"//label//"' Center='Node'>"
write(xmf,"(a,i6,a)") "<DataItem Format='XML' Datatype='Float' Dimensions='", nbs, "'>"
do is=1, nbs
  write(xmf,"(f12.6)") jorek%phi(is)
end do
write(xmf,"(a)") "</DataItem>"
write(xmf,"(a)") "</Attribute>"
write(xmf,"(a)") "</Grid>"
write(xmf,"(a)") "</Domain>"
write(xmf,"(a)") "</Xdmf>" 

close(xmf)

end subroutine plot_field

!---------------------------------------------------------------------------
!> @brief Write a plotmtv file
!> @details this file is readable plotmtv, it plots the mesh with nodes number
!> @param[IN] jorek the solver object that conatins all data
!> @param[IN] filename prefix of the mtv file
subroutine plot_jorek_field_2d_with_plotmtv(jorek, field_name)

type(sll_jorek_solver) :: jorek
character(len=*)       :: field_name
sll_int32              :: file_id
sll_real64             :: x1, y1
sll_int32              :: is1, is2, is3, is4
sll_int32              :: is, iq, nbs, nbq
sll_int32              :: error

nbs = jorek%n_nodes
nbq = jorek%n_quads

call sll_ascii_file_create(field_name//".mtv", file_id, error)

write(file_id,*)"$DATA=CURVE3D"
write(file_id,*)"%equalscale=T"
write(file_id,*)"%contfill"
write(file_id,*)"%toplabel='"//field_name//"'"

do iq = 1, nbq

  is1 = jorek%nodes(1,iq)
  is2 = jorek%nodes(2,iq)
  is3 = jorek%nodes(3,iq)
  is4 = jorek%nodes(4,iq)

  write(file_id,*) sngl(jorek%coord(1:2,is1)), sngl(jorek%phi(is1))
  write(file_id,*) sngl(jorek%coord(1:2,is2)), sngl(jorek%phi(is2))
  write(file_id,*) sngl(jorek%coord(1:2,is3)), sngl(jorek%phi(is3))
  write(file_id,*) sngl(jorek%coord(1:2,is4)), sngl(jorek%phi(is4))
  write(file_id,*) sngl(jorek%coord(1:2,is1)), sngl(jorek%phi(is1))
  write(file_id,*) 

end do

write(file_id,*)
do is = 1, nbs
  x1 = jorek%coord(1,is)
  y1 = jorek%coord(2,is)
  write(file_id,"(a)"   ,  advance="no")"@text x1="
  write(file_id,"(g15.3)", advance="no") x1
  write(file_id,"(a)"   ,  advance="no")" y1="
  write(file_id,"(g15.3)", advance="no") y1
  write(file_id,"(a)"   ,  advance="no")" z1=0. lc=5 ll='"
  write(file_id,"(i4)"  ,  advance="no") is
  write(file_id,"(a)")"'"
end do

write(file_id,*)"$end"
close(file_id)
   
end subroutine plot_jorek_field_2d_with_plotmtv


end module sll_jorek
