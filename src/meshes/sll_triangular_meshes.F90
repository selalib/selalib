!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: sll_triangular_meshes
!
! DESCRIPTION:
!> @file sll_triangular_meshes.F90
!>
!> @author
!> - Aurore Back
!> - Pierre Navaro
!>
!> @brief
!>  This module defines a triangular mesh.
!>
!> @details
!>  Triangular mesh 
!------------------------------------------------------------------------------
module sll_triangular_meshes
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_memory.h"
#include "sll_utilities.h"

use sll_meshes_base
use sll_tri_mesh_xmf

implicit none

integer :: imxref=99999999 

!> @brief 2d hexagonal mesh
!  vtaux  - composante x des vecteurs tangeants         
!  vtauy  - composante y des vecteurs tangeants        
type :: sll_triangular_mesh_2d

  sll_int32           :: num_nodes  
  sll_int32           :: num_cells 
  sll_int32           :: num_edges     
  sll_int32           :: num_bound     

  sll_real64, pointer :: coord(:,:)
  sll_int32,  pointer :: nodes(:,:)

  sll_real64          :: eta1_min
  sll_real64          :: eta1_max
  sll_real64          :: eta2_min
  sll_real64          :: eta2_max

  sll_int32           :: nbcoti
  sll_int32           :: nbtcot
  sll_int32           :: nmxfr
  sll_int32           :: nelfr
  sll_int32           :: nmxsd
  sll_int32           :: nctfrt
  sll_real64          :: petitl
  sll_real64          :: grandl

  sll_real64, dimension(:),   pointer :: aire
  sll_int32,  dimension(:),   pointer :: refs
  sll_int32,  dimension(:),   pointer :: reft
  sll_int32,  dimension(:,:), pointer :: nvois
  sll_int32,  dimension(:),   pointer :: nusd
  sll_int32,  dimension(:),   pointer :: npoel1
  sll_int32,  dimension(:),   pointer :: npoel2
  sll_int32,  dimension(:),   pointer :: krefro
  sll_int32,  dimension(:),   pointer :: kctfro
  sll_int32,  dimension(:),   pointer :: kelfro
  sll_int32,  dimension(:,:), pointer :: ksofro
  sll_real64, dimension(:,:), pointer :: vnofro
  sll_real64, dimension(:),   pointer :: xmal1
  sll_real64, dimension(:),   pointer :: xmal2
  sll_real64, dimension(:),   pointer :: xmal3
  sll_int32,  dimension(:,:), pointer :: nuvac
  sll_int32,  dimension(:),   pointer :: nugcv
  sll_int32,  dimension(:),   pointer :: nbcov
  sll_real64, dimension(:),   pointer :: xlcod

  sll_real64, dimension(:),   allocatable :: vtaux
  sll_real64, dimension(:),   allocatable :: vtauy
  sll_int32,  dimension (:),  allocatable :: nctfro
  sll_int32,  dimension (:),  allocatable :: nctfrp

!   contains

!     procedure, pass(mesh) :: eta1_node => eta1_node_hex
!     procedure, pass(mesh) :: eta2_node => eta2_node_hex
!     procedure, pass(mesh) :: eta1_cell_one_arg => eta1_cell_hex
!     procedure, pass(mesh) :: eta1_cell_two_arg => eta1_cell_triangular_two_arg
!     procedure, pass(mesh) :: eta2_cell_one_arg => eta2_cell_hex
!     procedure, pass(mesh) :: eta2_cell_two_arg => eta2_cell_triangular_two_arg
!     procedure, pass(mesh) :: display => display_triangular_mesh_2d
!     procedure, pass(mesh) :: delete => delete_triangular_mesh_2d

end type sll_triangular_mesh_2d

interface sll_create
   module procedure initialize_triangular_mesh_2d
end interface sll_create

interface sll_delete
   module procedure delete_triangular_mesh_2d
end interface sll_delete

interface sll_display
   module procedure display_triangular_mesh_2d
end interface sll_display

! interface eta1_cell
!    module procedure eta1_cell_one_arg, eta1_cell_two_arg
! end interface eta1_cell

! interface eta2_cell
!    module procedure eta2_cell_one_arg, eta2_cell_two_arg
! end interface eta2_cell


!Type: sll_triangular_mesh_2d
!Structure du maillage
!
!nbs    - nombre de sommets
!nbt    - nombre de triangles
!nbtcot - nombre total de cotes
!nbcoti - nombre de cotes internes
!num_bound  - nombre total de frontieres referencees
!nelfr  - nombre de triangles sur une frontiere
!coor   - coordonnees des sommets
!aire   - aires de elements
!xbas   - integrales des fonctions de base
!refs   - references des sommets
!reft   - references des elements
!ntri   - table de connectivite
!nvois  - numeros des voisins (solveur)
!nvoiv  - numeros des voisins (particules)
!nvoif  - numero local dans le voisin de la face commune
!nusd   - references du sous-domaine
!petitl - petite longueur de reference   
!grandl - grande longueur de reference    
!ncfrt  - nombre de cotes situes sur une frontiere
!nbcfli - nombre de cotes internes ne satisfaisant pas la CFL
!nndfnt - noeuds Dirichlet sur les frontieres internes       
!noefnt - noeuds Dirichlet sur les frontieres internes 
!irffnt - numeros de reference de ces noeuds Dirichlet
!nmxsd  - nombre de sous domaines references

contains

!> @brief allocates the memory space for a new 1D cartesian mesh on the heap,
!> initializes it with the given arguments and returns a pointer to the
!> object.
!> @param num_cells integer denoting the number of cells.
!> @param eta_min optional double precision value which represents the 
!> minimum value of the eta1 parameter in the cartesian mesh.
!> @param eta_max optional double precision value which represents the 
!> maximum value of the eta1 parameter in the cartesian mesh.
!> @return a pointer to the newly allocated object.
function new_triangular_mesh_2d( maafil ) result(m)

  type(sll_triangular_mesh_2d), pointer :: m
  character(len=*), intent(in)          :: maafil

  sll_int32 :: ierr
  SLL_ALLOCATE(m, ierr)
  call read_from_file(m, maafil)

end function new_triangular_mesh_2d

subroutine initialize_triangular_mesh_2d( mesh,     &
                                          nc_eta1,  &
                                          eta1_min, &
                                          eta1_max, &
                                          nc_eta2,  &
                                          eta2_min, &
                                          eta2_max)

type(sll_triangular_mesh_2d) :: mesh
sll_int32,  intent(in)       :: nc_eta1
sll_real64, intent(in)       :: eta1_min
sll_real64, intent(in)       :: eta1_max
sll_int32,  intent(in)       :: nc_eta2
sll_real64, intent(in)       :: eta2_min
sll_real64, intent(in)       :: eta2_max

!*-----------------------------------------------------------------------
!*  Generation du maillage pour une plaque rectangulaire dans le plan xOy.
!*  Maillage regulier de type TRIANGLES P1 - Q4T.
!*-----------------------------------------------------------------------

sll_int32 :: i
sll_int32 :: nd1, nd2, nd3
sll_int32 :: nelt
sll_int32 :: l, l1, ll1, ll2
sll_int32 :: in, iin, n, noeud, nsom, nxm, nym, neltot, ndd
sll_int32 :: nbox, nboy

sll_real64 :: xx1, xx2, pasx0, pasx1, pasy0, pasy1
sll_real64 :: alx, aly

mesh%eta1_min = eta1_min
mesh%eta1_max = eta1_max
mesh%eta2_min = eta2_min
mesh%eta2_max = eta2_max

nbox = nc_eta1+1
nboy = nc_eta2+1
ndd = max0(nbox,nboy)
alx = eta1_max-eta1_min
aly = eta2_max-eta2_min

!*---- Determination de quelques utilitaires

nxm    = nbox - 1      !nombre de quadrangles sur l'axe Ox du maillage
nym    = nboy - 1      !nombre de quadrangles sur l'axe Oy du maillage
neltot = 2 * nxm * nym !nombre total de triangles du maillage 
nsom   = nbox * nboy   !nombre de "sommets" du maillage quadrangulaire
noeud  = nsom          !nombre total de noeuds

mesh%num_cells = neltot
mesh%num_nodes = noeud
allocate(mesh%coord(1:2,noeud))
allocate(mesh%nodes(1:3,1:neltot))

!*---- Ecriture des coordonnees des sommets ----*!
                  
pasx0 = alx / nxm                   !pas le long de OX
pasx1 = pasx0 * 0.5                 !demi - pas       
pasy0 = aly / nym                   !pas le long de OY
pasy1 = pasy0 * 0.5                 !demi - pas       

do n = 0 , nym       !Coordonnees des noeuds
   in = nbox * n
   do i = 1 , nbox
      iin = in + i                  !numero du noeud
      xx1 = eta1_min + alx - (i - 1) * pasx0   !coordonnees du noeud
      xx2 = eta2_min + n * pasy0 
      mesh%coord(1,iin) = xx1
      mesh%coord(2,iin) = xx2
   end do
end do

!*---- Table de connectivite ----*!

nelt = 0                        !initialisation
do l = 1 , nym                  !boucle sur les lignes

   l1  = l - 1
   ll1 = l1 * nbox
   ll2 = l  * nbox 

   do i = 1 , nxm               !boucle sur les colonnes

      nd1 = ll1 + i             
      nd2 = ll2 + i             
      nd3 = ll2 + i + 1         

      nelt = nelt + 1           !premier e.f. du "quadrangle"

      mesh%nodes(1,nelt) = nd1       !premier  noeud
      mesh%nodes(2,nelt) = nd2       !deuxieme noeud
      mesh%nodes(3,nelt) = nd1 + 1   !troisieme noeud

      nelt = nelt + 1           !second e.f. du "quadrangle"
       
      mesh%nodes(1,nelt) = nd1 + 1   !premier noeud        
      mesh%nodes(2,nelt) = nd2       !deuxieme noeud
      mesh%nodes(3,nelt) = nd3       !troisieme noeud

   end do

end do

end subroutine initialize_triangular_mesh_2d

!  function eta1_node_triangular(mesh, i, j) result(res)
!    ! The coordinates (i, j) correspond to the (r1, r2) basis
!    ! This function returns the 1st coordinate on the cartesian system
!    class(sll_triangular_mesh_2d), intent(in) :: mesh
!    sll_int32, intent(in)  :: i
!    sll_int32, intent(in)  :: j
!    sll_real64 :: res
!
!    res = mesh%r1_x1*i + mesh%r2_x1*j + mesh%center_x1
!  end function eta1_node_triangular
!
!  !> @brief Computes the second coordinate of a given point
!  !> @details Computes the second coordinate on the cartesian system 
!  !> of a point which has for hexagonal coordinates (i,j)
!  !> @param i sll_int32 denoting the first hexagonal coordinate of a point
!  !> @param j sll_int32 denoting the second hexagonal coordinate of a point
!  !> returns res real containing the coordinate "eta2"
!  function eta2_node_triangular(mesh, i, j) result(res)
!    ! The coordinates (k1, k2) correspond to the (r1, r2) basis
!    ! This function the 2nd coordinate on the cartesian system
!    class(sll_triangular_mesh_2d), intent(in)     :: mesh
!    sll_int32, intent(in)  :: i
!    sll_int32, intent(in)  :: j
!    sll_real64  :: res
!
!    res = mesh%r1_x2*i + mesh%r2_x2*j + mesh%center_x2
!  end function eta2_node_triangular
!
!
!  function eta1_cell_triangular(mesh, cell_num) result(res)
!    ! The index num_ele corresponds to the index of triangle
!    ! This function returns the 1st coordinate on the cartesian system
!    ! of the center of the triangle at num_ele
!    class(sll_triangular_mesh_2d),intent(in)     :: mesh
!    sll_int32, intent(in)      :: cell_num
!    sll_real64 :: res
!
!    res = mesh%center_cartesian_coord(1, cell_num)
!  end function eta1_cell_triangular
!
!  function eta2_cell_triangular(mesh, cell_num) result(res)
!    ! The index num_ele corresponds to the index of triangle
!    ! This function returns the 2nd coordinate on the cartesian system
!    ! of the center of the triangle at num_ele
!    class(sll_triangular_mesh_2d),intent(in)     :: mesh
!    sll_int32, intent(in)      :: cell_num
!    sll_real64 :: res
!
!    res = mesh%center_cartesian_coord(2, cell_num)
!  end function eta2_cell_triangular
!
!
!  function eta1_cell_triangular_two_arg(mesh, i, j) result(res)
!    class(sll_triangular_mesh_2d),intent(in)     :: mesh
!    sll_int32, intent(in)      :: i, j
!    sll_real64 :: res
!
!    res = 0.0_f64
!    print *, "Error : eta1_cell for a triangular mesh only works with ONE parameter (num_cell)"
!    STOP
!  end function eta1_cell_triangular_two_arg
!
!  function eta2_cell_triangular_two_arg(mesh, i, j) result(res)
!    class(sll_triangular_mesh_2d),intent(in)     :: mesh
!    sll_int32, intent(in)      :: i, j
!    sll_real64 :: res
!
!    res = 0.0_f64
!    print *, "Error : eta2_cell for a triangular mesh only works with ONE parameter (num_cell)"
!    STOP
!  end function eta2_cell_triangular_two_arg
!
!
!  function cells_to_origin(k1, k2) result(val)
!    ! Takes the coordinates (k1,k2) on the (r1,r2) basis and 
!    ! returns the number of cells between that point and
!    ! the origin. If (k1, k2) = 0, val = 0
!    sll_int32, intent(in)   :: k1
!    sll_int32, intent(in)   :: k2
!    sll_int32               :: val
!
!    ! We compute the number of cells from point to center 
!    if (k1*k2 .gt. 0) then
!       val = max(abs(k1),abs(k2))
!    else
!       val = abs(k1) + abs(k2)
!    end if
!
!  end function cells_to_origin
!
!
!  function hex_to_global(mesh, k1, k2) result(val)
!    ! Takes the coordinates (k1,k2) on the (r1,r2) basis and 
!    ! returns global index of that mesh point.
!    ! By default the index of the center of the mesh is 0
!    ! Then following the r1 direction and a counter-clockwise motion
!    ! we assing an index to every point of the mesh.
!    class(sll_triangular_mesh_2d)      :: mesh
!    sll_int32, intent(in)   :: k1
!    sll_int32, intent(in)   :: k2
!    sll_int32               :: distance
!    sll_int32               :: index_tab
!    sll_int32               :: val
!    
!    distance = cells_to_origin(k1,k2)
!
!    ! Test if we are in domain
!    if (distance .le. mesh%num_cells) then
!
!       call index_triangular_to_global(mesh, k1, k2,index_tab)
!       val = mesh%global_indices(index_tab)
!
!    else
!       val = -1
!    end if
!
!  end function hex_to_global


!> Displays mesh information on the terminal
subroutine display_triangular_mesh_2d(mesh)
class(sll_triangular_mesh_2d), intent(in) :: mesh

write(*,"(/,(a))") '2D mesh : num_cells   num_nodes   num_edges '
write(*,"(10x,3(i6,9x),4(g13.3,1x))") &
     mesh%num_cells,  &
     mesh%num_nodes,  &
     mesh%num_edges,  &
     mesh%eta1_min,   &
     mesh%eta1_max,   &
     mesh%eta2_min,   &
     mesh%eta2_max

end subroutine display_triangular_mesh_2d

!  subroutine write_triangular_mesh_2d(mesh, name)
!    ! Writes the mesh information in a file named "name"
!    type(sll_triangular_mesh_2d), pointer :: mesh
!    character(len=*) :: name
!    sll_int32  :: i
!    sll_int32  :: num_nodes
!    sll_int32  :: k1, k2
!    sll_int32, parameter :: out_unit=20
!
!    open (unit=out_unit,file=name,action="write",status="replace")
!
!    num_nodes = mesh%num_nodes
!
!    ! Optional writing every mesh point and its cartesian coordinates :
!    !    write(*,"(/,(a))") 'hex mesh : num_pnt    x1     x2'
!
!    do i=1, num_nodes
!       k1 = mesh%global_to_hex1(i)
!       k2 = mesh%global_to_hex2(i)
!       write (out_unit, "(3(i6,1x),2(g13.3,1x))") i,                &
!            k1,                      &
!            k2,                      &
!            mesh%global_to_x1(i), &
!            mesh%global_to_x2(i)
!    end do
!
!    close(out_unit)
!  end subroutine write_triangular_mesh_2d
!
!  subroutine write_field_triangular_mesh(mesh, field, name)
!    ! Writes the points cartesian coordinates and
!    ! field(vector) values in a file named "name"
!    type(sll_triangular_mesh_2d), pointer :: mesh
!    sll_real64,dimension(:) :: field
!    character(len=*) :: name
!    sll_int32  :: i
!    sll_int32  :: num_nodes
!    sll_real64 :: x1, x2
!    sll_int32, parameter :: out_unit=20
!
!    open (unit=out_unit,file=name,action="write",status="replace")
!
!    num_nodes = mesh%num_nodes
!    do i=1, num_nodes
!       x1 = mesh%global_to_x1(i)
!       x2 = mesh%global_to_x2(i)
!       write (out_unit, "(3(g13.3,1x))") x1, &
!            x2, &
!            field(i)
!    end do
!
!    close(out_unit)
!  end subroutine write_field_triangular_mesh
!
!  subroutine write_field_triangular_mesh_xmf(mesh, field, name)
!    ! Writes the points cartesian coordinates and
!    ! field(vector) values in a file named "name"
!    type(sll_triangular_mesh_2d), pointer :: mesh
!    sll_real64,dimension(:) :: field
!    character(len=*) :: name
!    sll_int32  :: i
!    sll_int32  :: num_cells
!    sll_int32  :: num_nodes
!    sll_int32  :: out_unit
!    sll_real64, allocatable :: coord(:,:)
!    sll_int32,  allocatable :: nodes(:,:)
!    sll_int32  :: error
!    sll_real64 :: x1, x2
!
!    call sll_new_file_id(out_unit, error)
!
!    num_nodes = mesh%num_nodes
!    num_cells = mesh%num_cells
!    SLL_ALLOCATE(coord(2,num_nodes),error)
!    SLL_ALLOCATE(nodes(3,num_cells),error)
!
!    do i=1, num_nodes
!       coord(1,i) = mesh%global_to_x1(i)
!       coord(2,i) = mesh%global_to_x2(i)
!    end do
!
!    do i=1, num_cells
!       x1      = mesh%center_cartesian_coord(1, i)
!       x2      = mesh%center_cartesian_coord(2, i)
!       call get_cell_vertices_index( x1, x2, mesh, nodes(1,i), nodes(2,i), nodes(3,i))
!    end do
!
!    call write_tri_mesh_xmf(name, coor, nodes, num_pts_tot, num_cells, field, 'values')
!
!    close(out_unit)
!
!  end subroutine write_field_triangular_mesh_xmf
!

subroutine write_triangular_mesh_mtv(mesh, mtv_file)

type(sll_triangular_mesh_2d) :: mesh
sll_real64                   :: x1
sll_real64                   :: y1
character(len=*)             :: mtv_file
sll_int32                    :: out_unit
sll_int32                    :: error
sll_int32                    :: i

call sll_new_file_id(out_unit, error)

open( out_unit, file=mtv_file)

!--- Trace du maillage ---

write(out_unit,"(a)")"$DATA=CURVE3D"
write(out_unit,"(a)")"%equalscale=T"
write(out_unit,"(a)")"%toplabel='Maillage' "

do i = 1, mesh%num_cells

   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(2,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(3,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,*)

end do

!--- Numeros des noeuds et des triangles

write(out_unit,"(a)")"$DATA=CURVE3D"
write(out_unit,"(a)")"%equalscale=T"
write(out_unit,"(a)")"%toplabel='Numeros des noeuds et des triangles' "

do i = 1, mesh%num_cells
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(2,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(3,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,*)
end do

do i = 1, mesh%num_cells
   x1 = (  mesh%coord(1,mesh%nodes(1,i))  &
         + mesh%coord(1,mesh%nodes(2,i))  &
     + mesh%coord(1,mesh%nodes(3,i))    )/3.
   y1 = (  mesh%coord(2,mesh%nodes(1,i))  &
         + mesh%coord(2,mesh%nodes(2,i))  &
     + mesh%coord(2,mesh%nodes(3,i))    )/3.
   write(out_unit,"(a)"   , advance="no")"@text x1="
   write(out_unit,"(f8.5)", advance="no") x1
   write(out_unit,"(a)"   , advance="no")" y1="
   write(out_unit,"(f8.5)", advance="no") y1
   write(out_unit,"(a)"   , advance="no")" z1=0. lc=4 ll='"
   write(out_unit,"(i4)"  , advance="no") i
   write(out_unit,"(a)")"'"
end do

do i = 1, mesh%num_nodes
   x1 = mesh%coord(1,i)
   y1 = mesh%coord(2,i)
   write(out_unit,"(a)"   ,  advance="no")"@text x1="
   write(out_unit,"(g15.3)", advance="no") x1
   write(out_unit,"(a)"   ,  advance="no")" y1="
   write(out_unit,"(g15.3)", advance="no") y1
   write(out_unit,"(a)"   ,  advance="no")" z1=0. lc=5 ll='"
   write(out_unit,"(i4)"  ,  advance="no") i
   write(out_unit,"(a)")"'"
end do

!--- Numeros des noeuds 

write(out_unit,*)"$DATA=CURVE3D"
write(out_unit,*)"%equalscale=T"
write(out_unit,*)"%toplabel='Numeros des noeuds' "

do i = 1, mesh%num_cells
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(2,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(3,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,*)
end do

do i = 1, mesh%num_nodes
   x1 = mesh%coord(1,i)
   y1 = mesh%coord(2,i)
   write(out_unit,"(a)"   ,  advance="no")"@text x1="
   write(out_unit,"(g15.3)", advance="no") x1
   write(out_unit,"(a)"   ,  advance="no")" y1="
   write(out_unit,"(g15.3)", advance="no") y1
   write(out_unit,"(a)"   ,  advance="no")" z1=0. lc=5 ll='"
   write(out_unit,"(i4)"  ,  advance="no") i
   write(out_unit,"(a)")"'"
end do

!--- Numeros des triangles

write(out_unit,*)"$DATA=CURVE3D"
write(out_unit,*)"%equalscale=T"
write(out_unit,*)"%toplabel='Numeros des triangles' "

do i = 1, mesh%num_cells
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(2,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(3,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,*)
end do

do i = 1, mesh%num_cells
   x1 = (  mesh%coord(1,mesh%nodes(1,i))  &
         + mesh%coord(1,mesh%nodes(2,i))  &
     + mesh%coord(1,mesh%nodes(3,i))    )/3.
   y1 = (  mesh%coord(2,mesh%nodes(1,i))  &
         + mesh%coord(2,mesh%nodes(2,i))  &
     + mesh%coord(2,mesh%nodes(3,i))    )/3.
   write(out_unit,"(a)"   , advance="no")"@text x1="
   write(out_unit,"(g15.3)", advance="no") x1
   write(out_unit,"(a)"   , advance="no")" y1="
   write(out_unit,"(g15.3)", advance="no") y1
   write(out_unit,"(a)"   , advance="no")" z1=0. lc=4 ll='"
   write(out_unit,"(i4)"  , advance="no") i
   write(out_unit,"(a)")"'"
end do

write(out_unit,*)"$END"
close(out_unit)
 
end subroutine write_triangular_mesh_mtv

subroutine delete_triangular_mesh_2d( mesh )

class(sll_triangular_mesh_2d), intent(inout) :: mesh

print*, 'delete mesh'

end subroutine delete_triangular_mesh_2d

subroutine read_from_file(mesh, maafil)

type(sll_triangular_mesh_2d) :: mesh
character(len=*), intent(in)  :: maafil

sll_int32        :: nfmaa = 12
sll_int32        :: sll_err
sll_int32        :: i
sll_int32        :: j
integer          :: nelin
integer          :: nefro
integer          :: nmaill
integer          :: iout = 6

write(iout,"(/////10x,'>>> Read mesh from file <<<'/)")
open(nfmaa,file=maafil,status='OLD',err=80)
write(*,1050,advance='no') trim(maafil)

write(iout,"(10x,'Open the file'                      &
&        /10x,'Unit number',i3,'    fichier  ',a14)") &
&        nfmaa,maafil

read(nfmaa,*) 
read(nfmaa,*) nmaill,imxref
read(nfmaa,*) mesh%num_nodes, &
              mesh%num_cells, &
              mesh%nmxfr,     &
              mesh%nmxsd,     &
              nefro,          &
              nelin,          &
              mesh%nelfr

SLL_ALLOCATE(mesh%coord(1:2,mesh%num_nodes),   sll_err)
SLL_ALLOCATE(mesh%nodes(1:3,1:mesh%num_cells), sll_err)
SLL_ALLOCATE(mesh%refs(mesh%num_nodes),        sll_err)
SLL_ALLOCATE(mesh%reft(mesh%num_cells),        sll_err)
SLL_ALLOCATE(mesh%nusd(mesh%num_cells),        sll_err)
SLL_ALLOCATE(mesh%nvois(3,mesh%num_cells),     sll_err)

read(nfmaa,*) ((mesh%coord(i,j),i=1,2),j=1,mesh%num_nodes)
read(nfmaa,*)  (mesh%refs(i)   ,i=1,mesh%num_nodes) 
read(nfmaa,*) ((mesh%nodes(i,j),i=1,3),j=1,mesh%num_cells)
read(nfmaa,*) ((mesh%nvois(i,j),i=1,3),j=1,mesh%num_cells)
read(nfmaa,*)  (mesh%nusd(i)   ,i=1,mesh%num_cells) 

close(nfmaa)

write(iout,"(//,10x,'Nb de noeuds                : ',i10/       &
&              ,10x,'Nb de triangles             : ',i10/       &
&              ,10x,'Nb max de front referencees : ',i10/       &
&              ,10x,'Nb max de SD references     : ',i10/       &
&              ,10x,/'Nb de triangles ayant au moins 1 sommet'  &
&              ,' sur une frontiere : ',i10/)") mesh%num_nodes, &
                   mesh%num_cells,mesh%nmxfr,mesh%nmxsd,nefro

write(iout,"(//10x,'Nb d''elements internes     : ',i10/    &
          &   ,10x,'Nb d''elements frontieres   : ',i10/)") nelin,mesh%nelfr

1050 format(/' Read mesh from file  ', A, ' ?  Y')

return
80 continue
SLL_ERROR(' Input file  '//maafil//'  not found')

end subroutine read_from_file

end module sll_triangular_meshes
