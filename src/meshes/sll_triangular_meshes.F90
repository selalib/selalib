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
  sll_int32,  dimension(:,:), pointer :: nvoiv
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
!nodes   - table de connectivite
!nvois  - numeros des voisins (solveur)
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

  call compute_aires( m )

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

  call compute_aires( tri_mesh )

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

  call compute_aires( m )

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
integer          :: imxref

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
SLL_ERROR( 'Read from file', 'Input file  '//maafil//'  not found.' )

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





!> @brief 
!> Compute characterisitic origin in triangular mesh
!> @details
!! Localisation de l'origine de la caracteristique le maillage       
!! et application des conditions aux limites.         
!!              Les numeros des particules a retirer sont          
!!              contenus dans le tableau numpt                     
!!              (Numero d'element et de cote dans nelet et ncott)  
!!                                                                      
!!            pos1   - position x des particules                       
!!            pos2   - position y des particules                       
!!            vit1   - vitesse vx des particules                       
!!            vit2   - vitesse vy des particules                       
!!            xlm1   - 1ere coordonnee barycentrique                   
!!            xlm2   - 2eme coordonnee barycentrique                   
!!            xlm3   - 3eme coordonnee barycentrique                   
!!            nlpa   - numero des elements contenant les particules    
!!                                                                      
!!            coor   - coordonnees des noeuds                          
!!            nodes   - numero des sommets des triangles                
!!            nvois  - numero des voisins des elements                 
!!            aire   - aire de chaque element                          
!!                                                                      
!!            coef1  - tableau temporaire des determinants             
!!            coef2  - tableau temporaire des determinants             
!!            coef3  - tableau temporaire des determinants             
!!            coef4  - tableau temporaire                              
!!                                                                      
!!            p1loc  - tableau temporaire des nouvelles positions      
!!            p2loc  - tableau temporaire des nouvelles positions     
!!            v1loc  - tableau temporaire des nouvelles vitesses      
!!            v2loc  - tableau temporaire des nouvelles vitesses     
!!                                                                      
!!            numpt  - tableau auxiliaire contenant les numeros des    
!!                      particules a eliminer                           
!!            nelet  - tableau auxiliaire contenant les numeros des   
!!                      elements qui contenaient ces particules         
!!            ncott  - tableau auxiliaire contenant les numeros des    
!!                      cotes qui ont absorbe ces particules            
!!            v1pert, v2pert  - composantes des vitesses des particules  absorbees
!!                                                                      
!!            nlmloc - tableau auxiliaire contenant les numeros des    
!!                      elements ou l'on cherche les particules         
!!            numres - tableau auxiliaire contenant les numeros des    
!!                      particules non encore localisees                
!!                                                                      
!!            itest  - tableau auxiliaire pour preciser le            
!!                      comportement d'une particule :                  
!!            - si  itest=0 la particules reste dans son triangle       
!!            - si  itest=1,2,3 la particule traverse le cote 1,2ou3   
!!                               mais ne traverse pas de frontiere      
!!            - si  itest=11,12,13 la particule est absorbee par le     
!!                                 cote frontiere 1,2,3                
!!            - si  itest=21,22,23 la particule est reflechie par       
!!                                 le cote frontiere 1,2,3             
!!          Dans tous les cas le chiffre des unites de itest designe    
!!          le numero du cote traverse par la particule.                
!!          Le type de frontiere est defini dans le tableau nvoiv       
!!          qui a ete cree dans "RECVOI" :                              
!!            - si nvoiv(i,n) > 0  le cote i n'est pas une frontiere    
!!            - si nvoiv(i,n) = 0  le cote i absorbe les particules     
!!            - si nvoiv(i,n) =-1  le cote i reflechit les particules  
!!                                                                      
!!            nbpert - nombre de particules a eliminer                 
!!            nbp    - nombre de particules de l'espece consideree     
!!                                                                      
!!                                                                      
!!Auteurs:
!!    731-POSPAR      
!!
!! J. Segre - Version 1.0  Juillet 1988  
!! L. Arnaud/A. Adolf - Version 2.0  Octobre 1991 
!!                                 (ordre trigo + optimisation)         
!!----------------------------------------------------------------------
subroutine positions(mesh, f, ex, ey, dt )

type(sll_triangular_mesh_2d), intent(in)    :: mesh !< mesh
real(8), dimension(:),        intent(inout) :: f    !< distribution function on nodes
real(8), dimension(:),        intent(in)    :: ex   !< electric field in x1
real(8), dimension(:),        intent(in)    :: ey   !< electric field in x2
real(8),                      intent(in)    :: dt   !< time step

real(8) :: eps
integer :: nbpert
integer :: nbpres
integer :: nbp
integer :: num
integer :: nrest
integer :: nfin
integer :: nbpr

integer, dimension(:), allocatable :: numres

integer :: ip
integer :: jp

integer, dimension(:),   allocatable :: nbpama
integer, dimension(:),   allocatable :: iad1
integer, dimension(:),   allocatable :: indice
integer, dimension(:),   allocatable :: itabor
integer, dimension(:),   allocatable :: itest
integer, dimension(:),   allocatable :: nlmloc
integer, dimension(:),   allocatable :: nlpa
real(8), dimension(:,:), allocatable :: xlm
real(8), dimension(:,:), allocatable :: coef
real(8), dimension(:),   allocatable :: xp
real(8), dimension(:),   allocatable :: yp

real(8) :: pa1x, pa1y, pa2x, pa2y, pa3x, pa3y

eps   = -mesh%petitl**2
allocate(numres(mesh%num_nodes)); numres = 0

!Compute number of nodes inside the domain
nbp   = 0
do ip = 1, mesh%num_nodes
  if (mesh%refs(ip) == 0) then
    nbp        = nbp+1
    numres(ip) = ip
  end if
end do

allocate(nlpa(nbp));   nlpa   = 0
allocate(itest(nbp));  itest  = 0
allocate(coef(4,nbp))
allocate(xlm(3,nbp))
allocate(nlmloc(nbp))
allocate(xp(nbp))
allocate(yp(nbp))

!Set arbitrary the position of the characteristic origin
do ip = 1, nbp
  jp         = numres(ip)
  nlpa(ip)   = mesh%npoel2(mesh%npoel1(jp)+1)
end do

!Set new positions
do ip = 1, nbp
  jp         = numres(ip)
  xp(ip)     = mesh%coord(1,jp) + ex(jp) * dt
  yp(ip)     = mesh%coord(2,jp) + ey(jp) * dt
  nlmloc(ip) = nlpa(ip)
  numres(ip) = ip
end do
     
nbpres = nbp
nbpert = 0
num    = 0

do while( nbpres > 0 )

   !*** Premiere boucle dans le meme element

   nfin  = 0
   nrest = 0
   num   = num + 1

   if ( num == 100 ) then

      write(*,"(//2x,'Arret dans POSITIONS:')",advance='no')
      write(*,"('on n''arrive pas a positionner ')",advance='no')
      write(*,"(i4)",advance='no')nbpres
      write(*,"('  particules')",advance='no')
      write(*,"(/2x,'Particules non positionnees :')",advance='no')
      write(*,"(5x,'numero - position - vitesse - element'/)") 
      do ip=1,min(50,nbpres)
         write(*,"(i10,2x,4e12.4,i10)") numres(ip), ip
      end do
      write(*,"(/5x,a)") & 
      'Il y a certainement un probleme avec le solveur de champ'
      write(*,"(/5x,'(en general la CFL n''est pas verifiee)'/)")

      stop

   end if

   do ip = 1, nbpres    

      jp = numres(ip)

      pa1x = mesh%coord(1,mesh%nodes(1,nlmloc(jp))) - xp(jp)
      pa1y = mesh%coord(2,mesh%nodes(1,nlmloc(jp))) - yp(jp)
      pa2x = mesh%coord(1,mesh%nodes(2,nlmloc(jp))) - xp(jp)
      pa2y = mesh%coord(2,mesh%nodes(2,nlmloc(jp))) - yp(jp)
      pa3x = mesh%coord(1,mesh%nodes(3,nlmloc(jp))) - xp(jp)
      pa3y = mesh%coord(2,mesh%nodes(3,nlmloc(jp))) - yp(jp)

      coef(1,ip) = pa1x*pa2y - pa1y*pa2x
      coef(2,ip) = pa2x*pa3y - pa2y*pa3x
      coef(3,ip) = pa3x*pa1y - pa3y*pa1x

      if(      coef(1,jp) >= eps    &
         .and. coef(2,jp) >= eps    &
         .and. coef(3,jp) >= eps ) then

         nfin = nfin + 1

         xlm(1,jp) = 0.5 * coef(1,ip) / mesh%aire(nlmloc(jp))
         xlm(2,jp) = 0.5 * coef(2,ip) / mesh%aire(nlmloc(jp))
         xlm(3,jp) = 0.5 * coef(3,ip) / mesh%aire(nlmloc(jp))

         nlpa(jp)  = nlmloc(jp)
         itest(ip) = 0
         
         write(*,*) ip, ' found in ', nlmloc(ip)

      end if

   end do

   !*** Deuxieme boucle pour celles qui sont sorties

   nbpr = nbpres - nfin

   if( nbpr .ne. 0 ) then

      do ip = 1, nbpres

         jp = numres(ip)

         if (       coef(1,ip) <  eps   &
              .and. coef(2,ip) >= eps   &
              .and. coef(3,ip) >= eps   ) then

            !La particule a traverse le cote 1 = (A1-A2)
            itest(ip) = 1 + 10*(1-min(1,mesh%nvoiv(1,nlmloc(jp))))

         end if

         if (       coef(1,ip) >= eps   &
              .and. coef(2,ip) <  eps   &
              .and. coef(3,ip) >= eps   ) then
   
            !La particule a traverse le cote 2 = (A2-A3)
            itest(ip) = 2 + 10*(1-min(1,mesh%nvoiv(2,nlmloc(jp))))
 
         end if
   
         if (       coef(1,ip) >= eps   &
              .and. coef(2,ip) >= eps   &
              .and. coef(3,ip) <  eps   ) then
   
            !La particule a traverse le cote 3 = (A3-A1)
            itest(ip) = 3 + 10*(1-min(1,mesh%nvoiv(3,nlmloc(jp))))

         end if
      
         if (   coef(1,ip) < eps    &
          .and. coef(2,ip) < eps )  then

            !La particule a traverse le cote 1 ou 2 

            pa2x = mesh%coord(1,mesh%nodes(2,nlmloc(jp)))-xp(jp)
            pa2y = mesh%coord(2,mesh%nodes(2,nlmloc(jp)))-yp(jp)
       
            coef(4,ip) = pa2x*ey(jp) - pa2y*ex(jp)

            itest(ip) = 1 + max(0,nint(sign(1d0,coef(4,ip))))
            itest(ip) = itest(ip)  &
             + 10*(1-min(1,mesh%nvoiv(itest(ip),nlmloc(jp))))

         end if

         if (   coef(2,ip) < eps     &
          .and. coef(3,ip) < eps )  then

            !La particule a traverse le cote 2 ou 3 

            pa3x = mesh%coord(1,mesh%nodes(3,nlmloc(jp)))-xp(jp) 
            pa3y = mesh%coord(2,mesh%nodes(3,nlmloc(jp)))-yp(jp)
   
            coef(4,ip) = pa3x*ey(jp) - pa3y*ex(jp)

            itest(ip) = 2 + max(0,nint(sign(1d0,coef(4,ip))))
            itest(ip) = itest(ip)  &
                + 10*(1-min(1,mesh%nvoiv(itest(ip),nlmloc(jp))))
         end if

         if (    coef(3,ip) < eps    &
           .and. coef(1,ip) < eps )  then

            !La particule a traverse le cote 3 ou 1 

            pa1x = mesh%coord(1,mesh%nodes(1,nlmloc(jp)))-xp(jp) 
            pa1y = mesh%coord(2,mesh%nodes(1,nlmloc(jp)))-yp(jp)

            coef(4,ip) = pa1x*ey(jp) - pa1y*ex(jp)

            itest(ip) = 1 +mod(2+max(0,nint(sign(1d0,coef(4,ip)))),3)
            itest(ip) = itest(ip)   &
                + 10*(1-min(1,mesh%nvoiv(itest(ip),nlmloc(jp))))

         end if

      end do

      !*** Particules absorbees par une frontiere -------------------
     
      do ip=1,nbpres
   
        if( itest(ip) > 10 .and. itest(ip) < 14 )  then
           nbpert = nbpert + 1
        end if
   
      end do
   
      !*** Reorganisation des tableaux temporaires ------------------
      !    Particules traversant un cote interne ou reflechie
   
      do ip=1,nbpres
   
        if( (itest(ip) >  0 .and. itest(ip) <  4) .or.     &
            (itest(ip) > 20 .and. itest(ip) < 24)) then
          nrest = nrest + 1
          numres(nrest) = numres(ip)
          if (itest(ip)> 20) then
            nlmloc(nrest) = mesh%nvois(itest(ip)-20,nlmloc(jp))  &
                            * max(0,sign(1,10-itest(ip)))    &
                            + nlmloc(jp)* max(0,sign(1,itest(ip)-20))    
          else
            nlmloc(nrest) = mesh%nvois(itest(ip),nlmloc(jp)) &
                            * max(0,sign(1,10-itest(ip)))    &
                            + nlmloc(jp)* max(0,sign(1,itest(ip)-20))    
          end if
        end if
   
      end do

   end if

   nbpres = nrest

end do

      
!!             xlm1   - 1ere coordonnee barycentrique des particules 
!!             xlm2   - 2eme coordonnee barycentrique des particules
!!             xlm3   - 3eme coordonnee barycentrique des particules 
!!             nlpa   - numeros des triangles contenant les particules
!!             coor   - coordonnees des noeuds du maillage        
!!             nodes  - numero des sommets des triangles            
!!             xbas   - integrales des fonctions de base           
!!             rho    - densite de charge                        
!!             ad1    - tableau temporaire (adresse de la 1ere    
!!                      particule de chaque maille dans le tableau 
!!                      ordonne des particules)                     
!!             indice - tableau temporaire (incrementation du nombre 
!!                      de particules dja reperees)               
!!             itabor - tableau temporaire (numeros des particules 
!!                      ordonnes suivant les numeros des mailles)
!!             nbpama - tableau temporaire (nombre de particules  
!!                      par maille)                                 
!!             xaux1,xaux4 - tableaux auxiliaires utilises pour assurer la
!!                          continuite pour les frontieres internes      
!!                          "doubles" transparentes             
!!                                                            
!!             nflst  - etiquette logique du fichier "listing" 
!!             petitl - petite longueur de reference             
!!             nbt    - nombre de triangle du maillage             
!!             nbs    - nombre de noeuds   du maillage              
!!                                                                  
!!             nbpa   - nombre de particules par espece              
!!             pcharg - charge d'une particule elementaire de l'espece
!!             alprjt - Coefficient de compensation de charge d'espace 
!!                      par particule a appliquer a la projection       
!!             ldbprj - Activation du debugger en ligne             
!!             idbprj - Iteration pour le debut des impressions      
!!             jdbprj - Frequence des impressions                     
!!             kdbprj - Increment sur les noeuds du maillage           
!!                                                                      
!!Auteur:
!!   751-DENCOU  
!!
!
!real(8) :: phi1, phi2, phi3
!real(8) :: xm11, xm21, xm31
!real(8) :: xm12, xm22, xm32
!real(8) :: xm13, xm23, xm33
!real(8) :: charge
!
!integer :: nprest
!integer :: mpa, inum, ip1, ip2, ks, ind
!
!!Recherche du nombre de particules de chaque maille -------
allocate(nbpama(mesh%num_triangles))

nbpama = 0
do ip = 1 , nbp
   mpa         = nlpa(ip)
   nbpama(mpa) = nbpama(mpa) + 1
end do

!--- Determination de l'adresse de la premiere particule 
!    de chaque maille dans le tableau ordonne

allocate(iad1(mesh%num_triangles))

iad1 = 0
ks = 1
do it = 1 , mesh%num_triangles
   if ( nbpama(it) .ne. 0 )  then
     iad1(it) = ks
     ks       = ks + nbpama(it)
   end if
end do

!--- Construction du tableau ordonne des particules -----------

allocate(indice(mesh%num_triangles))
allocate(itabor(nbp))

indice = 0; itabor = 0

do ip = 1, nbp

   mpa         = nlpa(ip)
   ind         = iad1(mpa) + indice(mpa)
   itabor(ind) = ip
   indice(mpa) = indice(mpa) + 1

end do

nprest = nbp

f_out  = 0d0

do it = 1 , mesh%nbt
      
   xm11 = 0.
   xm12 = 0.
   xm13 = 0.

   !nbpama(it)  !Nombre de particules dans la maille numero it
            
   if (nbpama(it) .ne. 0) then

      ip1 = iad1(it)
      ip2 = ip1 + nbpama(it) - 1

      !Boucle sur les particules situees dans la maille courante

      do ip = ip1 , ip2 
            
         inum = itabor(ip)
       
         phi1 = xlm(2,inum) * f(inum) 
         phi2 = xlm(3,inum) * f(inum)
         phi3 = xlm(1,inum) * f(inum)
               
         xm11 = xm11 + phi1
         xm12 = xm12 + phi2
         xm13 = xm13 + phi3

      end do

      f_out(1,it) = f_out(1,it) + xm11 
      f_out(2,it) = f_out(2,it) + xm12 
      f_out(3,it) = f_out(3,it) + xm13 
   
      if (nbpama(it) == nprest ) exit

      nprest = nprest - nbpama(it)

   end if

end do

end subroutine positions

!=======================================================================


subroutine compute_aires( mesh )

type(sll_triangular_mesh_2d), intent(inout) :: mesh !< mesh

integer, dimension(:), allocatable :: indc

sll_int32 :: i, j, k
sll_int32 :: it, is, is1, is2, is3
sll_int32 :: ntmp, id1, nct, iel, ind, iel1, nel
sll_int32 :: i1, i2, i3
sll_int32 :: jel1, jel2, jel3, nel1, nel2, nel3

real(8)   :: airtot
real(8)   :: xlml, xlmu, ylml, ylmu
real(8)   :: lx1, lx2, ly1, ly2

!*** Calcul des longueurs de reference
#ifdef DEBUG
write(6,*)"*** Calcul des longueurs de reference ***"
#endif /* DEBUG */

xlml = minval(mesh%coord(1,:))
xlmu = maxval(mesh%coord(1,:))
ylml = minval(mesh%coord(2,:))
ylmu = maxval(mesh%coord(2,:))

mesh%petitl = 1.e-04 * min(xlmu-xlml,ylmu-ylml)/sqrt(float(mesh%num_nodes))
mesh%grandl = 1.e+04 * max(xlmu-xlml,ylmu-ylmu)

!*** Calcul des aires des triangles
#ifdef DEBUG
write(6,*)"*** Calcul des aires des triangles ***"
#endif /* DEBUG */

allocate(mesh%aire(mesh%num_triangles)); mesh%aire=0.0

airtot = 0.

do it = 1, mesh%num_triangles

   lx1 = mesh%coord(1,mesh%nodes(2,it))-mesh%coord(1,mesh%nodes(1,it))
   ly1 = mesh%coord(2,mesh%nodes(3,it))-mesh%coord(2,mesh%nodes(1,it))
   lx2 = mesh%coord(1,mesh%nodes(3,it))-mesh%coord(1,mesh%nodes(1,it))
   ly2 = mesh%coord(2,mesh%nodes(2,it))-mesh%coord(2,mesh%nodes(1,it))

   mesh%aire(it) = 0.5 * abs(lx1*ly1 - lx2*ly2)

   if( mesh%aire(it) <= 0. ) then
     write(6,*) " Triangle : ", it
     write(6,*) mesh%nodes(1,it), ":",mesh%coord(1:2,mesh%nodes(1,it))
     write(6,*) mesh%nodes(2,it), ":",mesh%coord(1:2,mesh%nodes(2,it))
     write(6,*) mesh%nodes(3,it), ":",mesh%coord(1:2,mesh%nodes(3,it))
     stop "Aire de triangle negative"
   end if

   airtot = airtot + mesh%aire(it)

end do

#ifdef DEBUG
write(6,"(/10x,'Longueurs de reference :',2E15.5/)") mesh%petitl,mesh%grandl
write(6,"(/10x,'Limites x du domaine   :',2E15.5/   &
&          10x,'Limites y du domaine   :',2E15.5/   &
&          10x,'Aire des triangles     :',E15.5 /)") xlml,xlmu,ylml,ylmu,airtot

write(6,*)"*** Calcul des voisins pour les particules ***"
#endif /* DEBUG */

! --- Gestion des triangles ayant un noeud en commun -----------
!if (ldebug) &
!write(iout,*)"*** Gestion des triangles ayant un noeud en commun ***"
 
! ... recherche des elements ayant un sommet commun
!     creation du tableau npoel1(i+1)  contenant le nombre de 
!     triangles ayant le noeud i en commun

allocate(mesh%npoel1(mesh%num_nodes+1))

mesh%npoel1 = 0
do i=1,mesh%num_triangles
   is1 = mesh%nodes(1,i)
   is2 = mesh%nodes(2,i)
   is3 = mesh%nodes(3,i)
   mesh%npoel1(is1+1) = mesh%npoel1(is1+1)+1
   mesh%npoel1(is2+1) = mesh%npoel1(is2+1)+1
   mesh%npoel1(is3+1) = mesh%npoel1(is3+1)+1
end do

! ... le tableau npoel1 devient le tableau donnant l'adresse 
!     dans npoel2 du dernier element dans la suite des triangles
!     communs a un noeud

mesh%npoel1(1)=0
do i=3,mesh%num_nodes+1
   mesh%npoel1(i)=mesh%npoel1(i-1)+mesh%npoel1(i)
end do

! ... creation du tableau npoel2 contenant sequentiellement les 
!     numeros des triangles ayant un noeud en commun      
!     le premier triangle s'appuyant sur le noeud i est
!     adresse par "npoel1(i)+1" 
!     le nombre de triangles ayant le noeud i en commun est
!     "npoel1(i+1)-npoel1(i)"


allocate(mesh%npoel2(mesh%npoel1(mesh%num_nodes+1)))
allocate(indc(mesh%num_nodes))

indc   = 1  !Le tableau temporaire indc doit etre initialise a 1

do it = 1,mesh%num_triangles
   do k = 1,3
      is = mesh%nodes(k,it)
      mesh%npoel2(mesh%npoel1(is)+indc(is)) = it
      indc(is) = indc(is)+1
   end do
end do

deallocate(indc)

! --- Recherche des numeros des triangles voisins d'un triangle 

do iel=1,mesh%num_triangles

  ! ... numeros des 3 sommets du triangle

  is1=mesh%nodes(1,iel)
  is2=mesh%nodes(2,iel)
  is3=mesh%nodes(3,iel)
  
  ! ... boucles imbriquees sur les elements pointant vers
  !     les 2 noeuds extremites de l'arete consideree
  !     Le voisin est le triangle commun (hormis iel)

  ! ... premiere arete (entre le sommet is1 et is2)

  nel1=mesh%npoel1(is1+1)-mesh%npoel1(is1) !nb de triangles communs a is1
  nel2=mesh%npoel1(is2+1)-mesh%npoel1(is2) !nb de triangles communs a is2

  loop1:do i1=1,nel1
    jel1=mesh%npoel2(mesh%npoel1(is1)+i1) !premier triangle is1
    if(jel1.ne.iel) then
      do i2=1,nel2
        jel2=mesh%npoel2(mesh%npoel1(is2)+i2)
        if(jel2 == jel1) then
          mesh%nvois(1,iel)  = jel1
          exit loop1
        end if
      end do
    end if
  end do loop1

  ! ... deuxieme arete (entre le sommet is2 et is3)

  nel2=mesh%npoel1(is2+1)-mesh%npoel1(is2)
  nel3=mesh%npoel1(is3+1)-mesh%npoel1(is3)

  loop2:do i2=1,nel2
    jel2=mesh%npoel2(mesh%npoel1(is2)+i2)
    if(jel2 /= iel) then
      do i3=1,nel3
        jel3=mesh%npoel2(mesh%npoel1(is3)+i3)
        if(jel3 == jel2) then
          mesh%nvois(2,iel)=jel2
          exit loop2
        end if
      end do
    end if
  end do loop2

  ! ... troisieme arete (entre le sommet is3 et is1)

  nel3=mesh%npoel1(is3+1)-mesh%npoel1(is3)
  nel1=mesh%npoel1(is1+1)-mesh%npoel1(is1)

  loop3:do i3=1,nel3
    jel3=mesh%npoel2(mesh%npoel1(is3)+i3)
    if(jel3 /= iel) then
      do i1=1,nel1
        jel1=mesh%npoel2(mesh%npoel1(is1)+i1)
        if(jel1 == jel3) then
          mesh%nvois(3,iel)=jel3
          exit loop3
        end if
      end do
    end if
  end do loop3

end do


! --- Rangement de npoel2 dans l'ordre trigonometrique ---------

do is=1,mesh%num_nodes

  nel = mesh%npoel1(is+1)-mesh%npoel1(is)

  if ( nel > 1 ) then

    !*** Noeuds internes (Numero de reference nul) ***

    if( mesh%refs(is) == 0) then

      ind =1
      iel1=mesh%npoel2(mesh%npoel1(is)+1)

      loop4:do iel=2,nel-1
        do j=1,3
          if(mesh%nodes(j,iel1) == is) nct=mod(j+1,3)+1
        end do

        iel1=mesh%nvois(nct,iel1)
        do id1=ind+1,nel
          if(iel1 == mesh%npoel2(mesh%npoel1(is)+id1)) then
            ind=ind+1
            ntmp=mesh%npoel2(mesh%npoel1(is)+ind)
            mesh%npoel2(mesh%npoel1(is)+ind)=iel1
            mesh%npoel2(mesh%npoel1(is)+id1)=ntmp
            cycle loop4
          end if
        end do
      end do loop4

     ! Noeuds frontieres

     else 

       ! --> Recherche du premier triangle dans l'ordre trigonometrique
       loop5:do id1=1,nel
         iel1=mesh%npoel2(mesh%npoel1(is)+id1)
         do j=1,3
           if(mesh%nvois(j,iel1).le.0 .and. mesh%nodes(j,iel1) == is) then
             ntmp=mesh%npoel2(mesh%npoel1(is)+1)
             mesh%npoel2(mesh%npoel1(is)+1)=iel1
             mesh%npoel2(mesh%npoel1(is)+id1)=ntmp
             exit loop5
           end if
         end do
       end do loop5
           
       ! --> Rangement des autres triangles dans l'ordre trigonometrique
       !     (s'il y en a plus que 2) 
       if(nel  > 2) then
         ind =1
         iel1=mesh%npoel2(mesh%npoel1(is)+1)
  
         loop6:do iel=2,nel-1
           do j=1,3
             if(mesh%nodes(j,iel1)==is) then
               nct=mod(j+1,3)+1
             end if
           end do

           iel1=mesh%nvois(nct,iel1)
  
           do id1=ind+1,nel
             if(iel1 == mesh%npoel2(mesh%npoel1(is)+id1)) then
               ind=ind+1
               ntmp=mesh%npoel2(mesh%npoel1(is)+ind)
               mesh%npoel2(mesh%npoel1(is)+ind)=iel1
               mesh%npoel2(mesh%npoel1(is)+id1)=ntmp
               cycle loop6
             end if
           end do

         end do loop6

      end if

    end if

  end if

end do

! --- Remplissage de "nvoiv" -----------------------------------
!     Identique a "nvois" sauf pour les aretes appartenant a 
!     une frontiere. Le chiffre correspond ici a un code pour 
!     le traitement des conditions aux limites sur les 
!     particules, alors que dans "nvois" ce chiffre est 
!     l'oppose du numero de reference de la frontiere concernee

allocate(mesh%nvoiv(3,mesh%num_triangles))
do i = 1,mesh%num_triangles
  do j = 1, 3
    if (mesh%nvois(j,i)>0) then
      mesh%nvoiv(j,i) = mesh%nvois(j,i)
    else
      mesh%nvoiv(j,i) =  0
    end if
  end do    
end do

end subroutine compute_aires

end module sll_triangular_meshes
