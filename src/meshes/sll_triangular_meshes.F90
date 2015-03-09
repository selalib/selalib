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
#include "sll_assert.h"
#include "sll_boundary_condition_descriptors.h"

use sll_meshes_base
use sll_tri_mesh_xmf

implicit none

integer :: imxref=99999999 

!> @brief 2d hexagonal mesh
!  vtaux  - composante x des vecteurs tangeants         
!  vtauy  - composante y des vecteurs tangeants        
type :: sll_triangular_mesh_2d

  sll_int32           :: num_nodes  
  sll_int32           :: num_triangles 
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

  logical :: analyzed = .false.

contains

   procedure, pass(mesh) :: global_to_x1
   procedure, pass(mesh) :: global_to_x2
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

interface new_triangular_mesh_2d
  module procedure new_triangular_mesh_2d_from_file
  module procedure new_triangular_mesh_2d_from_hex_mesh
  module procedure new_triangular_mesh_2d_from_square
end interface new_triangular_mesh_2d

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

!> @brief allocates the memory space for a new 2D triangular mesh on the heap,
!> initializes it with the given arguments and returns a pointer to the
!> object.
!> @param maafil file name with data
!> @return a pointer to the newly allocated object.
function new_triangular_mesh_2d_from_file( maafil ) result(m)

  type(sll_triangular_mesh_2d), pointer :: m
  character(len=*), intent(in)          :: maafil

  sll_int32 :: ierr
  SLL_ALLOCATE(m, ierr)
  call read_from_file(m, maafil)

end function new_triangular_mesh_2d_from_file

!> @brief allocates the memory space for a new 2D triangular mesh on the heap,
!> initializes it with the given hexagonal mesh and returns a pointer to the
!> object.
!> @param hex_mesh hexagonal mesh
!> @return a pointer to the newly allocated object.
function new_triangular_mesh_2d_from_hex_mesh( hex_mesh ) result(tri_mesh)

  use sll_hex_meshes
  type(sll_hex_mesh_2d), intent(in), pointer :: hex_mesh
  type(sll_triangular_mesh_2d),      pointer :: tri_mesh

  sll_int32  :: ierr
  sll_int32  :: is1, iv1 
  sll_int32  :: is2, iv2
  sll_int32  :: is3, iv3
  sll_int32  :: i
  sll_int32  :: nctfr
  sll_real64 :: x1
  sll_real64 :: y1
  sll_real64 :: xa, xb, xc
  sll_real64 :: ya, yb, yc
  sll_real64 :: det


  SLL_ALLOCATE(tri_mesh, ierr)

    
  tri_mesh%num_nodes     = hex_mesh%num_pts_tot
  tri_mesh%num_triangles = hex_mesh%num_triangles
  tri_mesh%nmxfr         = 1
  tri_mesh%nmxsd         = 1
  SLL_ALLOCATE(tri_mesh%coord(1:2,tri_mesh%num_nodes),       ierr)
  SLL_ALLOCATE(tri_mesh%nodes(1:3,1:tri_mesh%num_triangles), ierr)
  SLL_ALLOCATE(tri_mesh%refs(tri_mesh%num_nodes),            ierr)
  SLL_ALLOCATE(tri_mesh%nvois(1:3,1:tri_mesh%num_triangles), ierr)
  tri_mesh%refs      =  0
  tri_mesh%nvois     =  0

  do i = 1, hex_mesh%num_pts_tot
    tri_mesh%coord(1,i) = hex_mesh%global_to_x1(i)
    tri_mesh%coord(2,i) = hex_mesh%global_to_x2(i)
  end do

  nctfr = 0
  do i = 1, hex_mesh%num_triangles
    
    x1 = hex_mesh%center_cartesian_coord(1, i)
    y1 = hex_mesh%center_cartesian_coord(2, i)
    
    call get_cell_vertices_index( x1, y1, hex_mesh, is1, is2, is3)
    call get_neighbours(hex_mesh, i, iv1, iv2, iv3)

    xa = tri_mesh%coord(1,is1)
    ya = tri_mesh%coord(2,is1)
    xb = tri_mesh%coord(1,is2)
    yb = tri_mesh%coord(2,is2)
    xc = tri_mesh%coord(1,is3)
    yc = tri_mesh%coord(2,is3)
   
    det = 2.*((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
    
    if ( det > 0) then
      tri_mesh%nodes(:,i) = [is1,is2,is3]
      tri_mesh%nvois(:,i) = [iv1,iv2,iv3]
    else
      tri_mesh%nodes(:,i) = [is1,is3,is2]
      tri_mesh%nvois(:,i) = [iv3,iv2,iv1]
    end if

    if (tri_mesh%nvois(1,i) < 0) then
      tri_mesh%refs(tri_mesh%nodes(1,i)) = 1
      tri_mesh%refs(tri_mesh%nodes(2,i)) = 1
      nctfr = nctfr+1
    end if
    if (tri_mesh%nvois(2,i) < 0) then
      tri_mesh%refs(tri_mesh%nodes(2,i)) = 1
      tri_mesh%refs(tri_mesh%nodes(3,i)) = 1
      nctfr = nctfr+1
    end if
    if (tri_mesh%nvois(3,i) < 0) then
      tri_mesh%refs(tri_mesh%nodes(3,i)) = 1
      tri_mesh%refs(tri_mesh%nodes(1,i)) = 1
      nctfr = nctfr+1
    end if

  end do

end function new_triangular_mesh_2d_from_hex_mesh


!> @brief
!> Allocates the memory space for a new 2D triangular mesh on the heap,
!> @details
!> initializes it with the given arguments and returns a pointer to the
!> object.
!> @param[in] nc_eta1 number of cells along first direction
!> @param[in] eta1_min left edge first direction
!> @param[in] eta1_max right edge first direction
!> @param[in] nc_eta2 number of cells along second direction
!> @param[in] eta2_min left edge second direction
!> @param[in] eta2_max right edge second direction
!> @return a pointer to the newly allocated object.
function new_triangular_mesh_2d_from_square( nc_eta1,  &
                                             eta1_min, & 
                                             eta1_max, & 
                                             nc_eta2,  & 
                                             eta2_min, & 
                                             eta2_max, &
                                             bc ) result(m)

  type(sll_triangular_mesh_2d), pointer :: m
  sll_int32,  intent(in)                :: nc_eta1
  sll_real64, intent(in)                :: eta1_min
  sll_real64, intent(in)                :: eta1_max
  sll_int32,  intent(in)                :: nc_eta2
  sll_real64, intent(in)                :: eta2_min
  sll_real64, intent(in)                :: eta2_max
  sll_int32,  optional                  :: bc

  sll_int32 :: ierr
  SLL_ALLOCATE(m, ierr)

  if (present(bc)) then
    SLL_ASSERT(bc == SLL_PERIODIC)
  end if
 
  call initialize_triangular_mesh_2d( m,        &
                                      nc_eta1,  &
                                      eta1_min, &
                                      eta1_max, &
                                      nc_eta2,  &
                                      eta2_min, &
                                      eta2_max)


end function new_triangular_mesh_2d_from_square


subroutine initialize_triangular_mesh_2d( mesh,     &
                                          nc_eta1,  &
                                          eta1_min, &
                                          eta1_max, &
                                          nc_eta2,  &
                                          eta2_min, &
                                          eta2_max, &
                                          bc )

type(sll_triangular_mesh_2d) :: mesh
sll_int32,  intent(in)       :: nc_eta1
sll_real64, intent(in)       :: eta1_min
sll_real64, intent(in)       :: eta1_max
sll_int32,  intent(in)       :: nc_eta2
sll_real64, intent(in)       :: eta2_min
sll_real64, intent(in)       :: eta2_max
sll_int32,  optional         :: bc

!*-----------------------------------------------------------------------
!*  Generation du maillage pour une plaque rectangulaire dans le plan xOy.
!*  Maillage regulier de type TRIANGLES P1 - Q4T.
!*-----------------------------------------------------------------------

sll_int32 :: i, j
sll_int32 :: nd1, nd2, nd3
sll_int32 :: nelt
sll_int32 :: l, l1, ll1, ll2
sll_int32 :: in, iin, n, noeud, nsom, nxm, nym, neltot, ndd
sll_int32 :: nbox, nboy
sll_int32 :: iel, iev
sll_int32 :: nelin, nefro, nelfr, nhp, nlp
sll_int32 :: ierr

sll_real64 :: xx1, xx2, pasx0, pasx1, pasy0, pasy1
sll_real64 :: alx, aly
sll_int32, allocatable :: nar(:)

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

mesh%num_triangles = neltot
mesh%num_nodes     = noeud
SLL_CLEAR_ALLOCATE(mesh%coord(1:2,1:noeud),ierr)
SLL_ALLOCATE(mesh%nodes(1:3,1:neltot),ierr); mesh%nodes = 0
SLL_ALLOCATE(mesh%nvois(3,1:neltot),ierr); mesh%nvois = 0
SLL_ALLOCATE(mesh%refs(1:noeud),ierr); mesh%refs = 0

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

!*** Recherche des elements sur la frontiere
nefro = 4*(nbox-2) + 4*(nboy-2)
nelfr = 2*(nbox-1) + 2*(nboy-1) - 2
nelin = neltot - nelfr 

allocate(nar(2*(nbox+nboy)))

!---- Description des frontieres ----

nhp = nbox + 1
do i = 1 , nbox
   nar(i) = nhp - i 
   mesh%refs(nar(i)) = 1
end do        

nhp = nbox * nym
do i = 1 , nbox
   nar(i) = nhp + i 
   mesh%refs(nar(i)) = 3
end do        

nlp = nboy + 1
do i = 1 , nboy
   nar(i) = nbox * (nlp - i)
   mesh%refs(nar(i)) = 4
end do        

do i = 1 , nboy
   nar(i) = i * nbox - nxm
   mesh%refs(nar(i)) = 2
end do        

mesh%nmxfr = 0
!*** Initialisation des voisins ***

if (present(bc)) then

   do i = 1, nym
      iel = (i-1) * 2 * nym + 1
      iev =  iel + 2 * nym - 1 
      mesh%nvois(1,iel) = iev
      mesh%nvois(3,iev) = iel
   end do
   
   do i = 1, nxm
      iel = (i-1) * 2 + 1
      iev =  iel + 2 * nxm * (nym-1) + 1 
      mesh%nvois(3,iel) = iev
      mesh%nvois(2,iev) = iel
   end do

   nelfr = 0
   nefro = 0

else

   !*** Voisins en bas et en haut
   
   do i = 1, 2*(nbox-1), 2
      mesh%nvois(3,i      ) = -1
      mesh%nvois(2,neltot-i+1) = -3
   end do
   
   !*** Voisins de chaque cote lateral
   
   mesh%nvois(1,1) = -2
   nhp = 2*(nbox-1)
   do i = 1, nboy - 2
      j = i * nhp
      mesh%nvois(3,j  ) = -4
      mesh%nvois(1,j+1) = -2 
   end do
   mesh%nvois(3,neltot) = -4

   do i = 1, noeud
      if (mesh%refs(i) > mesh%nmxfr) mesh%nmxfr = mesh%nmxfr+1
   end do

end if

deallocate(nar)

end subroutine initialize_triangular_mesh_2d

!> Displays mesh information on the terminal
subroutine display_triangular_mesh_2d(mesh)
class(sll_triangular_mesh_2d), intent(in) :: mesh

sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
sll_int32  :: i, nctfrt, nbtcot

eta1_min = minval(mesh%coord(1,:))
eta1_max = maxval(mesh%coord(1,:))
eta2_min = minval(mesh%coord(2,:))
eta2_max = maxval(mesh%coord(2,:))

nctfrt=0
do i=1,mesh%num_triangles
   if (mesh%nvois(1,i)<0) nctfrt=nctfrt+1
   if (mesh%nvois(2,i)<0) nctfrt=nctfrt+1
   if (mesh%nvois(3,i)<0) nctfrt=nctfrt+1
end do

nbtcot = (3*mesh%num_triangles+nctfrt)/2

write(*,"(/,(a))") '2D mesh : num_triangles   num_nodes   num_edges '
write(*,"(10x,3(i6,9x),/,'Frame',/,4(g13.3,1x))") &
     mesh%num_triangles,  &
     mesh%num_nodes,      &
     nbtcot,              &
     eta1_min,            &
     eta1_max,            &
     eta2_min,            &
     eta2_max

end subroutine display_triangular_mesh_2d

function global_to_x1(mesh, i) result(x1)
! Takes the global index of the point (see hex_to_global(...) for conventions)
! returns the first coordinate (x1) on the cartesian basis
class(sll_triangular_mesh_2d) :: mesh
sll_int32  :: i
sll_real64 :: x1

x1 = mesh%coord(1, i)

end function global_to_x1

function global_to_x2(mesh, i) result(x2)
! Takes the global index of the point (see hex_to_global(...) for conventions)
! returns the second coordinate (x2) on the cartesian basis
class(sll_triangular_mesh_2d) :: mesh
sll_int32  :: i
sll_real64 :: x2

x2 = mesh%coord(2, i)

end function global_to_x2

subroutine write_triangular_mesh_mtv(mesh, mtv_file)

type(sll_triangular_mesh_2d) :: mesh
sll_real64                   :: x1, x2
sll_real64                   :: y1, y2
character(len=*)             :: mtv_file
sll_int32                    :: out_unit
sll_int32                    :: error
sll_int32                    :: i, j, k, l

call sll_new_file_id(out_unit, error)

open( out_unit, file=mtv_file)

!--- Trace du maillage ---

write(out_unit,"(a)")"$DATA=CURVE3D"
write(out_unit,"(a)")"%equalscale=T"
write(out_unit,"(a)")"%toplabel='Maillage' "

do i = 1, mesh%num_triangles

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

do i = 1, mesh%num_triangles
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(2,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(3,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,*)
end do

do i = 1, mesh%num_triangles
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

do i = 1, mesh%num_triangles
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

do i = 1, mesh%num_triangles
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(2,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(3,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nodes(1,i)),0.
   write(out_unit,*)
end do

do i = 1, mesh%num_triangles
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

write(out_unit,*)

write(out_unit,*)"$DATA=CURVE3D"
write(out_unit,*)"%equalscale=T"
write(out_unit,*)"%toplabel='References des frontieres' "

do i = 1, mesh%num_triangles
   write(out_unit,*) mesh%coord(1:2,mesh%nodes(1,i)), 0.0
   write(out_unit,*) mesh%coord(1:2,mesh%nodes(2,i)), 0.0
   write(out_unit,*) mesh%coord(1:2,mesh%nodes(3,i)), 0.0
   write(out_unit,*) mesh%coord(1:2,mesh%nodes(1,i)), 0.0
   write(out_unit,*)
end do

do i = 1, mesh%num_triangles
  do j = 1, 3
    if (mesh%nvois(j,i) < 0) then
      k = mod(j-1,3)+1   !Premier sommet
      l = mod(j  ,3)+1   !Deuxieme sommet
      x1 = mesh%coord(1,mesh%nodes(k,i))
      y1 = mesh%coord(2,mesh%nodes(k,i))
      x2 = mesh%coord(1,mesh%nodes(l,i))
      y2 = mesh%coord(2,mesh%nodes(l,i))
      write(out_unit,"(a)"    , advance="no")"@text x1="
      write(out_unit,"(g15.3)", advance="no") 0.5*(x1+x2)
      write(out_unit,"(a)"    , advance="no")" y1="
      write(out_unit,"(g15.3)", advance="no") 0.5*(y1+y2)
      write(out_unit,"(a)"    , advance="no")" z1=0. lc=5 ll='"
      write(out_unit,"(i1)"   , advance="no") -mesh%nvois(j,i)
      write(out_unit,"(a)")"'"
    end if
  end do
end do

write(out_unit,*)

write(out_unit,*)"$DATA=CURVE3D"
write(out_unit,*)"%equalscale=T"
write(out_unit,*)"%toplabel='References des noeuds' "

do i = 1, mesh%num_triangles
   write(out_unit,*) mesh%coord(1:2,mesh%nodes(1,i)), 0.0
   write(out_unit,*) mesh%coord(1:2,mesh%nodes(2,i)), 0.0
   write(out_unit,*) mesh%coord(1:2,mesh%nodes(3,i)), 0.0
   write(out_unit,*) mesh%coord(1:2,mesh%nodes(1,i)), 0.0
   write(out_unit,*)
end do

do i = 1, mesh%num_nodes
   x1 = mesh%coord(1,i)
   y1 = mesh%coord(2,i)
   write(out_unit,"(a)"    , advance="no")"@text x1="
   write(out_unit,"(g15.3)", advance="no") x1
   write(out_unit,"(a)"    , advance="no")" y1="
   write(out_unit,"(g15.3)", advance="no") y1
   write(out_unit,"(a)"    , advance="no")" z1=0. lc=5 ll='"
   write(out_unit,"(i1)"   , advance="no") mesh%refs(i)
   write(out_unit,"(a)")"'"
end do

if (mesh%analyzed) then

write(out_unit,*)

write(out_unit,*)"$DATA=CURVE3D"
write(out_unit,*)"%equalscale=T"
write(out_unit,*)"%toplabel='Polygones de Voronoi'"

do i = 1, mesh%num_triangles
  write(out_unit,*)"%linetype   = 1 # Solid Linetype (default=1)"
  write(out_unit,*)"%linewidth  = 1 # Linewidth      (default=1)"
  write(out_unit,*)"%linecolor  = 1 # Line Color     (default=1)"
  do j = 1, 3
    if( mesh%nvois(j,i) > 0 ) then
      call get_cell_center(mesh, i, x1, y1)
      write(out_unit,*) x1,y1,0.
      call get_cell_center(mesh, mesh%nvois(j,i), x1, y1)
      write(out_unit,*) x1,y1,0.
      write(out_unit,*)
    end if
  end do
end do
   
do i = 1, mesh%num_triangles

  write(out_unit,*) "%linetype  = 1 # Solid Linetype (default=1)"
  write(out_unit,*) "%linewidth = 1 # Linewidth      (default=1)"
  write(out_unit,*) "%linecolor = 2 # Line Color     (default=1)"
  write(out_unit,*) mesh%coord(1:2,mesh%nodes(1,i)),0.
  write(out_unit,*) mesh%coord(1:2,mesh%nodes(2,i)),0.
  write(out_unit,*) mesh%coord(1:2,mesh%nodes(3,i)),0.
  write(out_unit,*) mesh%coord(1:2,mesh%nodes(1,i)),0.
  write(out_unit,*)

end do


write(out_unit,*)
write(out_unit,*)"$DATA=CURVE3D"
write(out_unit,*)"%equalscale=T"
write(out_unit,*)"%toplabel='Numeros des noeuds et des cotes' "

do i = 1, mesh%nbtcot
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nuvac(1,i)),0.
   write(out_unit,"(3f10.5)")mesh%coord(:,mesh%nuvac(2,i)),0.
   write(out_unit,*)
   x1 = 0.5*(mesh%coord(1,mesh%nuvac(1,i))+mesh%coord(1,mesh%nuvac(2,i)))
   y1 = 0.5*(mesh%coord(2,mesh%nuvac(1,i))+mesh%coord(2,mesh%nuvac(2,i)))
   write(out_unit,"(a)"   ,  advance="no")"@text x1="
   write(out_unit,"(g15.3)", advance="no") x1
   write(out_unit,"(a)"   ,  advance="no")" y1="
   write(out_unit,"(g15.3)", advance="no") y1
   write(out_unit,"(a)"   ,  advance="no")" z1=0. lc=5 ll='"
   write(out_unit,"(i4)"  ,  advance="no") i
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

end if

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
read(nfmaa,*) mesh%num_nodes,     &
              mesh%num_triangles, &
              mesh%nmxfr,         &
              mesh%nmxsd,         &
              nefro,              &
              nelin,              &
              mesh%nelfr

SLL_ALLOCATE(mesh%coord(1:2,mesh%num_nodes),       sll_err)
SLL_ALLOCATE(mesh%nodes(1:3,1:mesh%num_triangles), sll_err)
SLL_ALLOCATE(mesh%refs(mesh%num_nodes),            sll_err)
SLL_ALLOCATE(mesh%reft(mesh%num_triangles),        sll_err)
SLL_ALLOCATE(mesh%nusd(mesh%num_triangles),        sll_err)
SLL_ALLOCATE(mesh%nvois(3,mesh%num_triangles),     sll_err)

read(nfmaa,*) ((mesh%coord(i,j),i=1,2),j=1,mesh%num_nodes)
read(nfmaa,*)  (mesh%refs(i)   ,i=1,mesh%num_nodes) 
read(nfmaa,*) ((mesh%nodes(i,j),i=1,3),j=1,mesh%num_triangles)
read(nfmaa,*) ((mesh%nvois(i,j),i=1,3),j=1,mesh%num_triangles)
read(nfmaa,*)  (mesh%nusd(i)   ,i=1,mesh%num_triangles) 

close(nfmaa)

write(iout,"(//,10x,'Nb de noeuds                : ',i10/       &
&              ,10x,'Nb de triangles             : ',i10/       &
&              ,10x,'Nb max de front referencees : ',i10/       &
&              ,10x,'Nb max de SD references     : ',i10/       &
&              ,10x,/'Nb de triangles ayant au moins 1 sommet'  &
&              ,' sur une frontiere : ',i10/)") mesh%num_nodes, &
                   mesh%num_triangles,mesh%nmxfr,mesh%nmxsd,nefro

write(iout,"(//10x,'Nb d''elements internes     : ',i10/    &
          &   ,10x,'Nb d''elements frontieres   : ',i10/)") nelin,mesh%nelfr

1050 format(/' Read mesh from file  ', A, ' ?  Y')

return
80 continue
SLL_ERROR(' Input file  '//maafil//'  not found')

end subroutine read_from_file

subroutine get_cell_center( mesh, iel, x1, x2)

class(sll_triangular_mesh_2d) :: mesh
sll_int32,  intent(in)        :: iel
sll_real64, intent(out)       :: x1
sll_real64, intent(out)       :: x2

sll_real64 :: xa, ya, xb, yb, xc, yc
sll_real64 :: det, syca, syba

xa = mesh%coord(1,mesh%nodes(1,iel))
ya = mesh%coord(2,mesh%nodes(1,iel))
xb = mesh%coord(1,mesh%nodes(2,iel))
yb = mesh%coord(2,mesh%nodes(2,iel))
xc = mesh%coord(1,mesh%nodes(3,iel))
yc = mesh%coord(2,mesh%nodes(3,iel))

det  = 2.*((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya))
syca = (yc-ya)*(xb*xb-xa*xa+yb*yb-ya*ya)
syba = (xb-xa)*(xc*xc-xa*xa+yc*yc-ya*ya)

x1 = (syca-(yb-ya)*(xc*xc-xa*xa+yc*yc-ya*ya))/det
x2 = (syba-(xc-xa)*(xb*xb-xa*xa+yb*yb-ya*ya))/det

end subroutine get_cell_center

!> Map an hexagonal mesh on circle
!> param[inout] mesh the triangular mesh built fron an hexagonal mesh
!> param[in]    num_cells is the num_cells parameter of the hexagonal mesh
subroutine map_to_circle( mesh, num_cells )

class(sll_triangular_mesh_2d), intent(inout) :: mesh
sll_int32,                     intent(in)    :: num_cells

sll_int32  :: i, j, cell
sll_real64 :: r, alpha

i = 2
do cell = 1, num_cells
  r     = cell * 1.0_f64 / num_cells
  alpha = sll_pi / 6.0_f64
  do j = 1, cell*6
     mesh%coord(1,i+j-1) = r * cos(alpha)
     mesh%coord(2,i+j-1) = r * sin(alpha)
     alpha =  alpha + sll_pi / (3.0_f64 * cell)
  end do
  i = i + cell*6
end do

end subroutine map_to_circle

end module sll_triangular_meshes
