!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: sll_m_triangular_meshes
!
! DESCRIPTION:
!> @file sll_m_triangular_meshes.F90
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
module sll_m_triangular_meshes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_hexagonal_meshes, only: &
      sll_s_get_cell_vertices_index, &
      sll_t_hex_mesh_2d

   implicit none

   public :: &
      sll_s_analyze_triangular_mesh, &
      sll_s_map_to_circle, &
      sll_o_new_triangular_mesh_2d, &
      sll_s_triangular_mesh_2d_init_from_file, &
      sll_s_triangular_mesh_2d_init_from_square, &
      sll_s_triangular_mesh_2d_free, &
      sll_s_read_from_file, &
      sll_o_create, &
      sll_o_delete, &
      sll_o_display, &
      sll_t_triangular_mesh_2d, &
      sll_s_write_triangular_mesh_mtv

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @brief 2d hexagonal mesh
!  vtaux  - coordinate x vectors tan
!  vtauy  - coordinate y vectors tan
   type :: sll_t_triangular_mesh_2d

      sll_int32           :: num_nodes
      sll_int32           :: num_triangles
      sll_int32           :: num_edges
      sll_int32           :: num_bound

      sll_real64, pointer :: coord(:, :)
      sll_int32, pointer :: nodes(:, :)

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

      sll_real64, dimension(:), pointer :: area
      sll_int32, dimension(:), pointer :: refs
      sll_int32, dimension(:), pointer :: reft
      sll_int32, dimension(:, :), pointer :: nvois
      sll_int32, dimension(:), pointer :: nusd
      sll_int32, dimension(:), pointer :: npoel1
      sll_int32, dimension(:), pointer :: npoel2
      sll_int32, dimension(:), pointer :: krefro
      sll_int32, dimension(:), pointer :: kctfro
      sll_int32, dimension(:), pointer :: kelfro
      sll_int32, dimension(:, :), pointer :: ksofro
      sll_real64, dimension(:, :), pointer :: vnofro
      sll_real64, dimension(:), pointer :: xmal1
      sll_real64, dimension(:), pointer :: xmal2
      sll_real64, dimension(:), pointer :: xmal3
      sll_int32, dimension(:, :), pointer :: nuvac
      sll_int32, dimension(:), pointer :: nugcv
      sll_int32, dimension(:), pointer :: nbcov
      sll_real64, dimension(:), pointer :: xlcod

      sll_real64, dimension(:), pointer :: vtaux
      sll_real64, dimension(:), pointer :: vtauy
      sll_int32, dimension(:), pointer :: nctfro
      sll_int32, dimension(:), pointer :: nctfrp

      logical :: analyzed = .false.

   contains

      procedure, pass(mesh) :: global_to_x1
      procedure, pass(mesh) :: global_to_x2

   end type sll_t_triangular_mesh_2d

   interface sll_o_create
      module procedure initialize_triangular_mesh_2d
   end interface sll_o_create

   interface sll_o_delete
      module procedure sll_s_triangular_mesh_2d_free
   end interface sll_o_delete

   interface sll_o_display
      module procedure display_triangular_mesh_2d
   end interface sll_o_display

   interface sll_o_new_triangular_mesh_2d
      module procedure new_triangular_mesh_2d_from_file
      module procedure new_triangular_mesh_2d_from_hex_mesh
      module procedure new_triangular_mesh_2d_from_square
   end interface sll_o_new_triangular_mesh_2d

!Type: sll_t_triangular_mesh_2d
!
!nbs    - nodes
!nbt    - triangles
!nbtcot - edges
!nbcoti - inside edges
!num_bound  - total references for BC
!nelfr  - triangles on boundary
!coor   - coords nodes
!area   - areas elements
!refs   - references nodes
!reft   - references elements
!nodes  - connectivity
!nvois  - neighbors
!nvoif  - number of edge in neighbor triangle
!nusd   - references of subdomain
!petitl - smaller length
!grandl - bigger length
!ncfrt  - edges on boundary
!nbcfli - edges inside with CFL problem
!nndfnt - nodes Dirichlet inside boubdary
!noefnt - nodes Dirichlet boundary
!irffnt - reference nodes Dirichlet
!nmxsd  - sub domains number

contains

!> @brief allocates the memory space for a new 2D triangular mesh on the heap,
!> initializes it with the given arguments and returns a pointer to the
!> object.
!> @param maafil file name with data
!> @return a pointer to the newly allocated object.
   function new_triangular_mesh_2d_from_file(maafil) result(m)

      type(sll_t_triangular_mesh_2d), pointer :: m
      character(len=*), intent(in)          :: maafil

      sll_int32 :: ierr
      SLL_ALLOCATE(m, ierr)
      call sll_s_read_from_file(m, maafil)

      call compute_areas(m)

   end function new_triangular_mesh_2d_from_file

!> @brief Initialize a new 2D triangular mesh
!> @param maafil file name with data
!> @return a pointer to the newly allocated object.
   subroutine sll_s_triangular_mesh_2d_init_from_file(m, maafil)

      type(sll_t_triangular_mesh_2d) :: m
      character(len=*), intent(in)   :: maafil

      call sll_s_read_from_file(m, maafil)
      call compute_areas(m)

   end subroutine sll_s_triangular_mesh_2d_init_from_file

!> @brief allocates the memory space for a new 2D triangular mesh on the heap,
!> initializes it with the given hexagonal mesh and returns a pointer to the
!> object.
!> @param hex_mesh hexagonal mesh
!> @return a pointer to the newly allocated object.
   function new_triangular_mesh_2d_from_hex_mesh(hex_mesh) result(tri_mesh)

      type(sll_t_hex_mesh_2d), intent(in), pointer :: hex_mesh
      type(sll_t_triangular_mesh_2d), pointer :: tri_mesh

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

      tri_mesh%num_nodes = hex_mesh%num_pts_tot
      tri_mesh%num_triangles = hex_mesh%num_triangles
      tri_mesh%nmxfr = 1
      tri_mesh%nmxsd = 1
      SLL_ALLOCATE(tri_mesh%coord(1:2, tri_mesh%num_nodes), ierr)
      SLL_ALLOCATE(tri_mesh%nodes(1:3, 1:tri_mesh%num_triangles), ierr)
      SLL_ALLOCATE(tri_mesh%refs(tri_mesh%num_nodes), ierr)
      SLL_ALLOCATE(tri_mesh%nvois(1:3, 1:tri_mesh%num_triangles), ierr)
      tri_mesh%refs = 0
      tri_mesh%nvois = 0

      do i = 1, hex_mesh%num_pts_tot
         tri_mesh%coord(1, i) = hex_mesh%global_to_x1(i)
         tri_mesh%coord(2, i) = hex_mesh%global_to_x2(i)
      end do

      nctfr = 0
      do i = 1, hex_mesh%num_triangles

         x1 = hex_mesh%center_cartesian_coord(1, i)
         y1 = hex_mesh%center_cartesian_coord(2, i)

         call sll_s_get_cell_vertices_index(x1, y1, hex_mesh, is1, is2, is3)
         call hex_mesh%get_neighbours(i, iv1, iv2, iv3)

         xa = tri_mesh%coord(1, is1)
         ya = tri_mesh%coord(2, is1)
         xb = tri_mesh%coord(1, is2)
         yb = tri_mesh%coord(2, is2)
         xc = tri_mesh%coord(1, is3)
         yc = tri_mesh%coord(2, is3)

         det = 2.0_f64*((xb - xa)*(yc - ya) - (xc - xa)*(yb - ya))

         if (det > 0) then
            tri_mesh%nodes(:, i) = [is1, is2, is3]
            tri_mesh%nvois(:, i) = [iv1, iv2, iv3]
         else
            tri_mesh%nodes(:, i) = [is1, is3, is2]
            tri_mesh%nvois(:, i) = [iv3, iv2, iv1]
         end if

         if (tri_mesh%nvois(1, i) < 0) then
            tri_mesh%refs(tri_mesh%nodes(1, i)) = 1
            tri_mesh%refs(tri_mesh%nodes(2, i)) = 1
            nctfr = nctfr + 1
         end if
         if (tri_mesh%nvois(2, i) < 0) then
            tri_mesh%refs(tri_mesh%nodes(2, i)) = 1
            tri_mesh%refs(tri_mesh%nodes(3, i)) = 1
            nctfr = nctfr + 1
         end if
         if (tri_mesh%nvois(3, i) < 0) then
            tri_mesh%refs(tri_mesh%nodes(3, i)) = 1
            tri_mesh%refs(tri_mesh%nodes(1, i)) = 1
            nctfr = nctfr + 1
         end if

      end do

      call compute_areas(tri_mesh)

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
   function new_triangular_mesh_2d_from_square(nc_eta1, &
                                               eta1_min, &
                                               eta1_max, &
                                               nc_eta2, &
                                               eta2_min, &
                                               eta2_max, &
                                               bc) result(m)

      type(sll_t_triangular_mesh_2d), pointer :: m
      sll_int32, intent(in)                :: nc_eta1
      sll_real64, intent(in)                :: eta1_min
      sll_real64, intent(in)                :: eta1_max
      sll_int32, intent(in)                :: nc_eta2
      sll_real64, intent(in)                :: eta2_min
      sll_real64, intent(in)                :: eta2_max
      sll_int32, optional                  :: bc

      sll_int32 :: ierr
      SLL_ALLOCATE(m, ierr)

      if (present(bc)) then
         SLL_ASSERT(bc == sll_p_periodic)
      end if

      call initialize_triangular_mesh_2d(m, &
                                         nc_eta1, &
                                         eta1_min, &
                                         eta1_max, &
                                         nc_eta2, &
                                         eta2_min, &
                                         eta2_max)

      call compute_areas(m)

   end function new_triangular_mesh_2d_from_square

!> @brief
!> Initialize a new 2D triangular mesh
!> @param[in] nc_eta1 number of cells along first direction
!> @param[in] eta1_min left edge first direction
!> @param[in] eta1_max right edge first direction
!> @param[in] nc_eta2 number of cells along second direction
!> @param[in] eta2_min left edge second direction
!> @param[in] eta2_max right edge second direction
!> @return a pointer to the newly allocated object.
   subroutine sll_s_triangular_mesh_2d_init_from_square(m, &
                                                        nc_eta1, eta1_min, eta1_max, nc_eta2, &
                                                        eta2_min, eta2_max, bc)

      type(sll_t_triangular_mesh_2d) :: m
      sll_int32, intent(in)         :: nc_eta1
      sll_real64, intent(in)         :: eta1_min
      sll_real64, intent(in)         :: eta1_max
      sll_int32, intent(in)         :: nc_eta2
      sll_real64, intent(in)         :: eta2_min
      sll_real64, intent(in)         :: eta2_max
      sll_int32, optional           :: bc

      if (present(bc)) then
         SLL_ASSERT(bc == sll_p_periodic)
      end if

      call initialize_triangular_mesh_2d(m, &
                                         nc_eta1, &
                                         eta1_min, &
                                         eta1_max, &
                                         nc_eta2, &
                                         eta2_min, &
                                         eta2_max)

      call compute_areas(m)

   end subroutine sll_s_triangular_mesh_2d_init_from_square

   subroutine initialize_triangular_mesh_2d(mesh, &
                                            nc_eta1, &
                                            eta1_min, &
                                            eta1_max, &
                                            nc_eta2, &
                                            eta2_min, &
                                            eta2_max, &
                                            bc)

      type(sll_t_triangular_mesh_2d) :: mesh
      sll_int32, intent(in)         :: nc_eta1
      sll_real64, intent(in)         :: eta1_min
      sll_real64, intent(in)         :: eta1_max
      sll_int32, intent(in)         :: nc_eta2
      sll_real64, intent(in)         :: eta2_min
      sll_real64, intent(in)         :: eta2_max
      sll_int32, optional           :: bc

!*-----------------------------------------------------------------------
!*  Mesh generation of rectangular domain
!*  TRIANGLES P1 - Q4T.
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

      nbox = nc_eta1 + 1
      nboy = nc_eta2 + 1
      ndd = max(nbox, nboy)
      alx = eta1_max - eta1_min
      aly = eta2_max - eta2_min

      nxm = nbox - 1      !number quadrangles x
      nym = nboy - 1      !number quadrangles y
      neltot = 2*nxm*nym !number total triangles
      nsom = nbox*nboy   !number nodes
      noeud = nsom          !number total nodes

      mesh%num_triangles = neltot
      mesh%num_nodes = noeud
      SLL_CLEAR_ALLOCATE(mesh%coord(1:2, 1:noeud), ierr)
      SLL_ALLOCATE(mesh%nodes(1:3, 1:neltot), ierr); mesh%nodes = 0
      SLL_ALLOCATE(mesh%nvois(3, 1:neltot), ierr); mesh%nvois = 0
      SLL_ALLOCATE(mesh%refs(1:noeud), ierr); mesh%refs = 0

      pasx0 = alx/real(nxm, f64)
      pasx1 = pasx0*0.5_f64
      pasy0 = aly/real(nym, f64)
      pasy1 = pasy0*0.5_f64

      do n = 0, nym
         in = nbox*n
         do i = 1, nbox
            iin = in + i
            xx1 = eta1_min + alx - real(i - 1, f64)*pasx0
            xx2 = eta2_min + real(n, f64)*pasy0
            mesh%coord(1, iin) = xx1
            mesh%coord(2, iin) = xx2
         end do
      end do

      nelt = 0
      do l = 1, nym

         l1 = l - 1
         ll1 = l1*nbox
         ll2 = l*nbox

         do i = 1, nxm

            nd1 = ll1 + i
            nd2 = ll2 + i
            nd3 = ll2 + i + 1

            nelt = nelt + 1

            mesh%nodes(1, nelt) = nd1
            mesh%nodes(2, nelt) = nd2
            mesh%nodes(3, nelt) = nd1 + 1

            nelt = nelt + 1

            mesh%nodes(1, nelt) = nd1 + 1
            mesh%nodes(2, nelt) = nd2
            mesh%nodes(3, nelt) = nd3

         end do

      end do

      nefro = 4*(nbox - 2) + 4*(nboy - 2)
      nelfr = 2*(nbox - 1) + 2*(nboy - 1) - 2
      nelin = neltot - nelfr

      allocate (nar(2*(nbox + nboy)))

      nhp = nbox + 1
      do i = 1, nbox
         nar(i) = nhp - i
         mesh%refs(nar(i)) = 1
      end do

      nhp = nbox*nym
      do i = 1, nbox
         nar(i) = nhp + i
         mesh%refs(nar(i)) = 3
      end do

      nlp = nboy + 1
      do i = 1, nboy
         nar(i) = nbox*(nlp - i)
         mesh%refs(nar(i)) = 4
      end do

      do i = 1, nboy
         nar(i) = i*nbox - nxm
         mesh%refs(nar(i)) = 2
      end do

      mesh%nmxfr = 0

      if (present(bc)) then

         do i = 1, nym
            iel = (i - 1)*2*nym + 1
            iev = iel + 2*nym - 1
            mesh%nvois(1, iel) = iev
            mesh%nvois(3, iev) = iel
         end do

         do i = 1, nxm
            iel = (i - 1)*2 + 1
            iev = iel + 2*nxm*(nym - 1) + 1
            mesh%nvois(3, iel) = iev
            mesh%nvois(2, iev) = iel
         end do

         nelfr = 0
         nefro = 0

      else

         do i = 1, 2*(nbox - 1), 2
            mesh%nvois(3, i) = -1
            mesh%nvois(2, neltot - i + 1) = -3
         end do

         mesh%nvois(1, 1) = -2
         nhp = 2*(nbox - 1)
         do i = 1, nboy - 2
            j = i*nhp
            mesh%nvois(3, j) = -4
            mesh%nvois(1, j + 1) = -2
         end do
         mesh%nvois(3, neltot) = -4

         do i = 1, noeud
            if (mesh%refs(i) > mesh%nmxfr) mesh%nmxfr = mesh%nmxfr + 1
         end do

      end if

      deallocate (nar)

   end subroutine initialize_triangular_mesh_2d

!> Displays mesh information on the terminal
   subroutine display_triangular_mesh_2d(mesh)
      class(sll_t_triangular_mesh_2d), intent(in) :: mesh

      sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
      sll_int32  :: i, nctfrt, nbtcot

      eta1_min = minval(mesh%coord(1, :))
      eta1_max = maxval(mesh%coord(1, :))
      eta2_min = minval(mesh%coord(2, :))
      eta2_max = maxval(mesh%coord(2, :))

      nctfrt = 0
      do i = 1, mesh%num_triangles
         if (mesh%nvois(1, i) < 0) nctfrt = nctfrt + 1
         if (mesh%nvois(2, i) < 0) nctfrt = nctfrt + 1
         if (mesh%nvois(3, i) < 0) nctfrt = nctfrt + 1
      end do

      nbtcot = (3*mesh%num_triangles + nctfrt)/2

      write (*, "(/,(a))") '2D mesh : num_triangles   num_nodes   num_edges '
      write (*, "(10x,3(i6,9x),/,'Frame',/,4(g13.3,1x))") &
         mesh%num_triangles, &
         mesh%num_nodes, &
         nbtcot, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max

   end subroutine display_triangular_mesh_2d

   function global_to_x1(mesh, i) result(x1)
! Takes the global index of the point (see hex_to_global(...) for conventions)
! returns the first coordinate (x1) on the cartesian basis
      class(sll_t_triangular_mesh_2d) :: mesh
      sll_int32  :: i
      sll_real64 :: x1

      x1 = mesh%coord(1, i)

   end function global_to_x1

   function global_to_x2(mesh, i) result(x2)
! Takes the global index of the point (see hex_to_global(...) for conventions)
! returns the second coordinate (x2) on the cartesian basis
      class(sll_t_triangular_mesh_2d) :: mesh
      sll_int32  :: i
      sll_real64 :: x2

      x2 = mesh%coord(2, i)

   end function global_to_x2

   subroutine sll_s_write_triangular_mesh_mtv(mesh, mtv_file)

      type(sll_t_triangular_mesh_2d), intent(in) :: mesh
      character(len=*), intent(in) :: mtv_file

      sll_real64 :: x1, x2
      sll_real64 :: y1, y2
      sll_int32  :: out_unit
      sll_int32  :: i, j, k, l

      open (file=mtv_file, status='replace', form='formatted', newunit=out_unit)

!--- Trace du maillage ---

      write (out_unit, "(a)") "$DATA=CURVE3D"
      write (out_unit, "(a)") "%equalscale=T"
      write (out_unit, "(a)") "%toplabel='Maillage' "

      do i = 1, mesh%num_triangles

         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(1, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(2, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(3, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(1, i)), 0.
         write (out_unit, *)

      end do

!--- Numeros des noeuds et des triangles

      write (out_unit, "(a)") "$DATA=CURVE3D"
      write (out_unit, "(a)") "%equalscale=T"
      write (out_unit, "(a)") "%toplabel='Numeros des noeuds et des triangles' "

      do i = 1, mesh%num_triangles
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(1, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(2, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(3, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(1, i)), 0.
         write (out_unit, *)
      end do

      do i = 1, mesh%num_triangles
         x1 = (mesh%coord(1, mesh%nodes(1, i)) &
               + mesh%coord(1, mesh%nodes(2, i)) &
               + mesh%coord(1, mesh%nodes(3, i)))/3.0_f64
         y1 = (mesh%coord(2, mesh%nodes(1, i)) &
               + mesh%coord(2, mesh%nodes(2, i)) &
               + mesh%coord(2, mesh%nodes(3, i)))/3.0_f64
         write (out_unit, "(a)", advance="no") "@text x1="
         write (out_unit, "(f8.5)", advance="no") x1
         write (out_unit, "(a)", advance="no") " y1="
         write (out_unit, "(f8.5)", advance="no") y1
         write (out_unit, "(a)", advance="no") " z1=0. lc=4 ll='"
         write (out_unit, "(i4)", advance="no") i
         write (out_unit, "(a)") "'"
      end do

      do i = 1, mesh%num_nodes
         x1 = mesh%coord(1, i)
         y1 = mesh%coord(2, i)
         write (out_unit, "(a)", advance="no") "@text x1="
         write (out_unit, "(g15.3)", advance="no") x1
         write (out_unit, "(a)", advance="no") " y1="
         write (out_unit, "(g15.3)", advance="no") y1
         write (out_unit, "(a)", advance="no") " z1=0. lc=5 ll='"
         write (out_unit, "(i4)", advance="no") i
         write (out_unit, "(a)") "'"
      end do

!--- Numeros des noeuds

      write (out_unit, *) "$DATA=CURVE3D"
      write (out_unit, *) "%equalscale=T"
      write (out_unit, *) "%toplabel='Numeros des noeuds' "

      do i = 1, mesh%num_triangles
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(1, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(2, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(3, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(1, i)), 0.
         write (out_unit, *)
      end do

      do i = 1, mesh%num_nodes
         x1 = mesh%coord(1, i)
         y1 = mesh%coord(2, i)
         write (out_unit, "(a)", advance="no") "@text x1="
         write (out_unit, "(g15.3)", advance="no") x1
         write (out_unit, "(a)", advance="no") " y1="
         write (out_unit, "(g15.3)", advance="no") y1
         write (out_unit, "(a)", advance="no") " z1=0. lc=5 ll='"
         write (out_unit, "(i4)", advance="no") i
         write (out_unit, "(a)") "'"
      end do

!--- Numeros des triangles

      write (out_unit, *) "$DATA=CURVE3D"
      write (out_unit, *) "%equalscale=T"
      write (out_unit, *) "%toplabel='Numeros des triangles' "

      do i = 1, mesh%num_triangles
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(1, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(2, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(3, i)), 0.
         write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nodes(1, i)), 0.
         write (out_unit, *)
      end do

      do i = 1, mesh%num_triangles
         x1 = (mesh%coord(1, mesh%nodes(1, i)) &
               + mesh%coord(1, mesh%nodes(2, i)) &
               + mesh%coord(1, mesh%nodes(3, i)))/3.0_f64
         y1 = (mesh%coord(2, mesh%nodes(1, i)) &
               + mesh%coord(2, mesh%nodes(2, i)) &
               + mesh%coord(2, mesh%nodes(3, i)))/3.0_f64
         write (out_unit, "(a)", advance="no") "@text x1="
         write (out_unit, "(g15.3)", advance="no") x1
         write (out_unit, "(a)", advance="no") " y1="
         write (out_unit, "(g15.3)", advance="no") y1
         write (out_unit, "(a)", advance="no") " z1=0. lc=4 ll='"
         write (out_unit, "(i4)", advance="no") i
         write (out_unit, "(a)") "'"
      end do

      write (out_unit, *)

      write (out_unit, *) "$DATA=CURVE3D"
      write (out_unit, *) "%equalscale=T"
      write (out_unit, *) "%toplabel='References des frontieres' "

      do i = 1, mesh%num_triangles
         write (out_unit, *) mesh%coord(1:2, mesh%nodes(1, i)), 0.0
         write (out_unit, *) mesh%coord(1:2, mesh%nodes(2, i)), 0.0
         write (out_unit, *) mesh%coord(1:2, mesh%nodes(3, i)), 0.0
         write (out_unit, *) mesh%coord(1:2, mesh%nodes(1, i)), 0.0
         write (out_unit, *)
      end do

      do i = 1, mesh%num_triangles
         do j = 1, 3
            if (mesh%nvois(j, i) < 0) then
               k = mod(j - 1, 3) + 1   !Premier sommet
               l = mod(j, 3) + 1   !Deuxieme sommet
               x1 = mesh%coord(1, mesh%nodes(k, i))
               y1 = mesh%coord(2, mesh%nodes(k, i))
               x2 = mesh%coord(1, mesh%nodes(l, i))
               y2 = mesh%coord(2, mesh%nodes(l, i))
               write (out_unit, "(a)", advance="no") "@text x1="
               write (out_unit, "(g15.3)", advance="no") 0.5_f64*(x1 + x2)
               write (out_unit, "(a)", advance="no") " y1="
               write (out_unit, "(g15.3)", advance="no") 0.5_f64*(y1 + y2)
               write (out_unit, "(a)", advance="no") " z1=0. lc=5 ll='"
               write (out_unit, "(i1)", advance="no") - mesh%nvois(j, i)
               write (out_unit, "(a)") "'"
            end if
         end do
      end do

      write (out_unit, *)

      write (out_unit, *) "$DATA=CURVE3D"
      write (out_unit, *) "%equalscale=T"
      write (out_unit, *) "%toplabel='References des noeuds' "

      do i = 1, mesh%num_triangles
         write (out_unit, *) mesh%coord(1:2, mesh%nodes(1, i)), 0.0
         write (out_unit, *) mesh%coord(1:2, mesh%nodes(2, i)), 0.0
         write (out_unit, *) mesh%coord(1:2, mesh%nodes(3, i)), 0.0
         write (out_unit, *) mesh%coord(1:2, mesh%nodes(1, i)), 0.0
         write (out_unit, *)
      end do

      do i = 1, mesh%num_nodes
         x1 = mesh%coord(1, i)
         y1 = mesh%coord(2, i)
         write (out_unit, "(a)", advance="no") "@text x1="
         write (out_unit, "(g15.3)", advance="no") x1
         write (out_unit, "(a)", advance="no") " y1="
         write (out_unit, "(g15.3)", advance="no") y1
         write (out_unit, "(a)", advance="no") " z1=0. lc=5 ll='"
         write (out_unit, "(i1)", advance="no") mesh%refs(i)
         write (out_unit, "(a)") "'"
      end do

      if (mesh%analyzed) then

         write (out_unit, *)

         write (out_unit, *) "$DATA=CURVE3D"
         write (out_unit, *) "%equalscale=T"
         write (out_unit, *) "%toplabel='Polygones de Voronoi'"

         do i = 1, mesh%num_triangles
            write (out_unit, *) "%linetype   = 1 # Solid Linetype (default=1)"
            write (out_unit, *) "%linewidth  = 1 # Linewidth      (default=1)"
            write (out_unit, *) "%linecolor  = 1 # Line Color     (default=1)"
            do j = 1, 3
               if (mesh%nvois(j, i) > 0) then
                  call get_cell_center(mesh, i, x1, y1)
                  write (out_unit, *) x1, y1, 0.
                  call get_cell_center(mesh, mesh%nvois(j, i), x1, y1)
                  write (out_unit, *) x1, y1, 0.
                  write (out_unit, *)
               end if
            end do
         end do

         do i = 1, mesh%num_triangles

            write (out_unit, *) "%linetype  = 1 # Solid Linetype (default=1)"
            write (out_unit, *) "%linewidth = 1 # Linewidth      (default=1)"
            write (out_unit, *) "%linecolor = 2 # Line Color     (default=1)"
            write (out_unit, *) mesh%coord(1:2, mesh%nodes(1, i)), 0.
            write (out_unit, *) mesh%coord(1:2, mesh%nodes(2, i)), 0.
            write (out_unit, *) mesh%coord(1:2, mesh%nodes(3, i)), 0.
            write (out_unit, *) mesh%coord(1:2, mesh%nodes(1, i)), 0.
            write (out_unit, *)

         end do

         write (out_unit, *)
         write (out_unit, *) "$DATA=CURVE3D"
         write (out_unit, *) "%equalscale=T"
         write (out_unit, *) "%toplabel='Numeros des noeuds et des cotes' "

         do i = 1, mesh%nbtcot
            write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nuvac(1, i)), 0.
            write (out_unit, "(3f10.5)") mesh%coord(:, mesh%nuvac(2, i)), 0.
            write (out_unit, *)
            x1 = 0.5_f64*(mesh%coord(1, mesh%nuvac(1, i)) + mesh%coord(1, mesh%nuvac(2, i)))
            y1 = 0.5_f64*(mesh%coord(2, mesh%nuvac(1, i)) + mesh%coord(2, mesh%nuvac(2, i)))
            write (out_unit, "(a)", advance="no") "@text x1="
            write (out_unit, "(g15.3)", advance="no") x1
            write (out_unit, "(a)", advance="no") " y1="
            write (out_unit, "(g15.3)", advance="no") y1
            write (out_unit, "(a)", advance="no") " z1=0. lc=5 ll='"
            write (out_unit, "(i4)", advance="no") i
            write (out_unit, "(a)") "'"
         end do

         do i = 1, mesh%num_nodes
            x1 = mesh%coord(1, i)
            y1 = mesh%coord(2, i)
            write (out_unit, "(a)", advance="no") "@text x1="
            write (out_unit, "(g15.3)", advance="no") x1
            write (out_unit, "(a)", advance="no") " y1="
            write (out_unit, "(g15.3)", advance="no") y1
            write (out_unit, "(a)", advance="no") " z1=0. lc=5 ll='"
            write (out_unit, "(i4)", advance="no") i
            write (out_unit, "(a)") "'"
         end do

      end if

      write (out_unit, *) "$END"
      close (out_unit)

   end subroutine sll_s_write_triangular_mesh_mtv

   subroutine sll_s_triangular_mesh_2d_free(mesh)

      class(sll_t_triangular_mesh_2d), intent(inout) :: mesh

      SLL_ASSERT(mesh%num_nodes > 0)

   end subroutine sll_s_triangular_mesh_2d_free

   subroutine sll_s_read_from_file(mesh, maafil)

      type(sll_t_triangular_mesh_2d) :: mesh
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

      write (iout, "(/////10x,'>>> Read mesh from file <<<'/)")
      print *, " will open : ", maafil
      open (nfmaa, file=maafil, status='OLD', err=80)
      print *, " opened    : ", maafil
      write (*, 1050, advance='no') trim(maafil)

      write (iout, "(10x,'Open the file'                      &
      &        /10x,'Unit number',i3,'    fichier  ',a14)") &
      &        nfmaa, maafil

      read (nfmaa, *)
      read (nfmaa, *) nmaill, imxref
      read (nfmaa, *) mesh%num_nodes, &
         mesh%num_triangles, &
         mesh%nmxfr, &
         mesh%nmxsd, &
         nefro, &
         nelin, &
         mesh%nelfr

      SLL_ALLOCATE(mesh%coord(1:2, mesh%num_nodes), sll_err)
      SLL_ALLOCATE(mesh%nodes(1:3, 1:mesh%num_triangles), sll_err)
      SLL_ALLOCATE(mesh%refs(mesh%num_nodes), sll_err)
      SLL_ALLOCATE(mesh%reft(mesh%num_triangles), sll_err)
      SLL_ALLOCATE(mesh%nusd(mesh%num_triangles), sll_err)
      SLL_ALLOCATE(mesh%nvois(3, mesh%num_triangles), sll_err)

      read (nfmaa, *) ((mesh%coord(i, j), i=1, 2), j=1, mesh%num_nodes)
      read (nfmaa, *) (mesh%refs(i), i=1, mesh%num_nodes)
      read (nfmaa, *) ((mesh%nodes(i, j), i=1, 3), j=1, mesh%num_triangles)
      read (nfmaa, *) ((mesh%nvois(i, j), i=1, 3), j=1, mesh%num_triangles)
      read (nfmaa, *) (mesh%nusd(i), i=1, mesh%num_triangles)

      close (nfmaa)

      write (iout, "(//,10x,'Nb of nodes                 : ',i10/       &
      &              ,10x,'Nb of triangles             : ',i10/       &
      &              ,10x,'Nb max of BC                : ',i10/       &
      &              ,10x,'Nb max of subdomain         : ',i10/       &
      &              ,10x,/'Nb triangles with 1 node'  &
      &              ,' on boundary : ',i10/)") mesh%num_nodes, &
                         mesh%num_triangles, mesh%nmxfr, mesh%nmxsd, nefro

      write (iout, "(//10x,'Nb elements inside     : ',i10/    &
                &   ,10x,'Nb elements boundary   : ',i10/)") nelin, mesh%nelfr

1050  format(/' Read mesh from file  ', A, ' ?  Y')

      return
80    continue
      SLL_ERROR('Read from file', 'Input file  '//maafil//'  not found.')

   end subroutine sll_s_read_from_file

   subroutine get_cell_center(mesh, iel, x1, x2)

      class(sll_t_triangular_mesh_2d) :: mesh
      sll_int32, intent(in)        :: iel
      sll_real64, intent(out)       :: x1
      sll_real64, intent(out)       :: x2

      sll_real64 :: xa, ya, xb, yb, xc, yc
      sll_real64 :: det, syca, syba

      xa = mesh%coord(1, mesh%nodes(1, iel))
      ya = mesh%coord(2, mesh%nodes(1, iel))
      xb = mesh%coord(1, mesh%nodes(2, iel))
      yb = mesh%coord(2, mesh%nodes(2, iel))
      xc = mesh%coord(1, mesh%nodes(3, iel))
      yc = mesh%coord(2, mesh%nodes(3, iel))

      det = 2.0_f64*((xb - xa)*(yc - ya) - (xc - xa)*(yb - ya))
      syca = (yc - ya)*(xb*xb - xa*xa + yb*yb - ya*ya)
      syba = (xb - xa)*(xc*xc - xa*xa + yc*yc - ya*ya)

      x1 = (syca - (yb - ya)*(xc*xc - xa*xa + yc*yc - ya*ya))/det
      x2 = (syba - (xc - xa)*(xb*xb - xa*xa + yb*yb - ya*ya))/det

   end subroutine get_cell_center

!> Map an hexagonal mesh on circle
!> param[inout] mesh the triangular mesh built fron an hexagonal mesh
!> param[in]    num_cells is the num_cells parameter of the hexagonal mesh
!> param[in]    order, optional if order=1 we move only points on the
!> boundary
   subroutine sll_s_map_to_circle(mesh, num_cells, order)

      class(sll_t_triangular_mesh_2d), intent(inout) :: mesh
      sll_int32, intent(in)    :: num_cells
      sll_int32, optional                          :: order

      sll_int32  :: i, j, cell, layer
      sll_real64 :: r, alpha

      if (present(order)) then
         layer = order
      else
         layer = num_cells
      end if

      i = 2
      do cell = 1, num_cells
         if (cell > num_cells - layer) then
            r = real(cell, f64)*1.0_f64/real(num_cells, f64)
            alpha = sll_p_pi/6.0_f64
            do j = 1, cell*6
               mesh%coord(1, i + j - 1) = r*cos(alpha)
               mesh%coord(2, i + j - 1) = r*sin(alpha)
               alpha = alpha + sll_p_pi/real(3*cell, f64)
            end do
         end if
         i = i + cell*6
      end do

   end subroutine sll_s_map_to_circle

!=======================================================================

   subroutine compute_areas(mesh)

      type(sll_t_triangular_mesh_2d), intent(inout) :: mesh !< mesh

      integer, dimension(:), allocatable :: indc

      sll_int32 :: i, j, k
      sll_int32 :: it, is, is1, is2, is3
      sll_int32 :: ntmp, id1, nct, iel, ind, iel1, nel
      sll_int32 :: i1, i2, i3
      sll_int32 :: jel1, jel2, jel3, nel1, nel2, nel3

      real(8)   :: airtot
      real(8)   :: xlml, xlmu, ylml, ylmu
      real(8)   :: lx1, lx2, ly1, ly2

      xlml = minval(mesh%coord(1, :))
      xlmu = maxval(mesh%coord(1, :))
      ylml = minval(mesh%coord(2, :))
      ylmu = maxval(mesh%coord(2, :))

      mesh%petitl = 1.d-04*min(xlmu - xlml, ylmu - ylml)/sqrt(real(mesh%num_nodes, 8))
      mesh%grandl = 1.d+04*max(xlmu - xlml, ylmu - ylmu)

      allocate (mesh%area(mesh%num_triangles)); mesh%area = 0.0_f64

      airtot = 0._f64

      do it = 1, mesh%num_triangles

         lx1 = mesh%coord(1, mesh%nodes(2, it)) - mesh%coord(1, mesh%nodes(1, it))
         ly1 = mesh%coord(2, mesh%nodes(3, it)) - mesh%coord(2, mesh%nodes(1, it))
         lx2 = mesh%coord(1, mesh%nodes(3, it)) - mesh%coord(1, mesh%nodes(1, it))
         ly2 = mesh%coord(2, mesh%nodes(2, it)) - mesh%coord(2, mesh%nodes(1, it))

         mesh%area(it) = 0.5*abs(lx1*ly1 - lx2*ly2)

         if (mesh%area(it) <= 0.) then
            write (6, *) " Triangle : ", it
            write (6, *) mesh%nodes(1, it), ":", mesh%coord(1:2, mesh%nodes(1, it))
            write (6, *) mesh%nodes(2, it), ":", mesh%coord(1:2, mesh%nodes(2, it))
            write (6, *) mesh%nodes(3, it), ":", mesh%coord(1:2, mesh%nodes(3, it))
            stop "triangle negative area"
         end if

         airtot = airtot + mesh%area(it)

      end do

! --- triangles with one node in common  -----------

!search for elements with a common summit
!creation of table npoel1(i+1) containing the number of
!triangles having node i in common

      allocate (mesh%npoel1(mesh%num_nodes + 1))

      mesh%npoel1 = 0
      do i = 1, mesh%num_triangles
         is1 = mesh%nodes(1, i)
         is2 = mesh%nodes(2, i)
         is3 = mesh%nodes(3, i)
         mesh%npoel1(is1 + 1) = mesh%npoel1(is1 + 1) + 1
         mesh%npoel1(is2 + 1) = mesh%npoel1(is2 + 1) + 1
         mesh%npoel1(is3 + 1) = mesh%npoel1(is3 + 1) + 1
      end do

! table npoel1 becomes the table giving the address
!  in npoel2 of the last element in the sequence of triangles
!  common to a node

      mesh%npoel1(1) = 0
      do i = 3, mesh%num_nodes + 1
         mesh%npoel1(i) = mesh%npoel1(i - 1) + mesh%npoel1(i)
      end do

! creation of the npoel2 table containing sequentially the
! triangles numbers with a common node
! the first triangle leaning on node i is
! address by "npoel1(i)+1"
! the number of triangles having node i in common is
! "npoel1(i+1)-npoel1(i)"

      allocate (mesh%npoel2(mesh%npoel1(mesh%num_nodes + 1)))
      allocate (indc(mesh%num_nodes))

      indc = 1

      do it = 1, mesh%num_triangles
         do k = 1, 3
            is = mesh%nodes(k, it)
            mesh%npoel2(mesh%npoel1(is) + indc(is)) = it
            indc(is) = indc(is) + 1
         end do
      end do

      deallocate (indc)

! --- search neighbors

      do iel = 1, mesh%num_triangles

         is1 = mesh%nodes(1, iel)
         is2 = mesh%nodes(2, iel)
         is3 = mesh%nodes(3, iel)

         !nested loops on elements pointing towards
         !the 2 end nodes of the considered ridge
         !The neighbour is the common triangle (except iel)
         !first ridge (between vertex is1 and is2)

         nel1 = mesh%npoel1(is1 + 1) - mesh%npoel1(is1) !nb de triangles communs a is1
         nel2 = mesh%npoel1(is2 + 1) - mesh%npoel1(is2) !nb de triangles communs a is2

         loop1: do i1 = 1, nel1
            jel1 = mesh%npoel2(mesh%npoel1(is1) + i1) !premier triangle is1
            if (jel1 .ne. iel) then
               do i2 = 1, nel2
                  jel2 = mesh%npoel2(mesh%npoel1(is2) + i2)
                  if (jel2 == jel1) then
                     mesh%nvois(1, iel) = jel1
                     exit loop1
                  end if
               end do
            end if
         end do loop1

         ! ... second edge (is2-is3)

         nel2 = mesh%npoel1(is2 + 1) - mesh%npoel1(is2)
         nel3 = mesh%npoel1(is3 + 1) - mesh%npoel1(is3)

         loop2: do i2 = 1, nel2
            jel2 = mesh%npoel2(mesh%npoel1(is2) + i2)
            if (jel2 /= iel) then
               do i3 = 1, nel3
                  jel3 = mesh%npoel2(mesh%npoel1(is3) + i3)
                  if (jel3 == jel2) then
                     mesh%nvois(2, iel) = jel2
                     exit loop2
                  end if
               end do
            end if
         end do loop2

         ! ... third edge (is3-is1)

         nel3 = mesh%npoel1(is3 + 1) - mesh%npoel1(is3)
         nel1 = mesh%npoel1(is1 + 1) - mesh%npoel1(is1)

         loop3: do i3 = 1, nel3
            jel3 = mesh%npoel2(mesh%npoel1(is3) + i3)
            if (jel3 /= iel) then
               do i1 = 1, nel1
                  jel1 = mesh%npoel2(mesh%npoel1(is1) + i1)
                  if (jel1 == jel3) then
                     mesh%nvois(3, iel) = jel3
                     exit loop3
                  end if
               end do
            end if
         end do loop3

      end do

! --- npoel2 ---------

      do is = 1, mesh%num_nodes

         nel = mesh%npoel1(is + 1) - mesh%npoel1(is)

         if (nel > 1) then

            !*** internal nodes (ref=0) ***

            if (mesh%refs(is) == 0) then

               ind = 1
               iel1 = mesh%npoel2(mesh%npoel1(is) + 1)

               loop4: do iel = 2, nel - 1
                  do j = 1, 3
                     if (mesh%nodes(j, iel1) == is) nct = mod(j + 1, 3) + 1
                  end do

                  iel1 = mesh%nvois(nct, iel1)
                  do id1 = ind + 1, nel
                     if (iel1 == mesh%npoel2(mesh%npoel1(is) + id1)) then
                        ind = ind + 1
                        ntmp = mesh%npoel2(mesh%npoel1(is) + ind)
                        mesh%npoel2(mesh%npoel1(is) + ind) = iel1
                        mesh%npoel2(mesh%npoel1(is) + id1) = ntmp
                        cycle loop4
                     end if
                  end do
               end do loop4

               ! boundary nodes

            else

               ! --> first triangle
               loop5: do id1 = 1, nel
                  iel1 = mesh%npoel2(mesh%npoel1(is) + id1)
                  do j = 1, 3
                     if (mesh%nvois(j, iel1) .le. 0 .and. mesh%nodes(j, iel1) == is) then
                        ntmp = mesh%npoel2(mesh%npoel1(is) + 1)
                        mesh%npoel2(mesh%npoel1(is) + 1) = iel1
                        mesh%npoel2(mesh%npoel1(is) + id1) = ntmp
                        exit loop5
                     end if
                  end do
               end do loop5

               if (nel > 2) then
                  ind = 1
                  iel1 = mesh%npoel2(mesh%npoel1(is) + 1)

                  loop6: do iel = 2, nel - 1
                     do j = 1, 3
                        if (mesh%nodes(j, iel1) == is) then
                           nct = mod(j + 1, 3) + 1
                        end if
                     end do

                     iel1 = mesh%nvois(nct, iel1)

                     do id1 = ind + 1, nel
                        if (iel1 == mesh%npoel2(mesh%npoel1(is) + id1)) then
                           ind = ind + 1
                           ntmp = mesh%npoel2(mesh%npoel1(is) + ind)
                           mesh%npoel2(mesh%npoel1(is) + ind) = iel1
                           mesh%npoel2(mesh%npoel1(is) + id1) = ntmp
                           cycle loop6
                        end if
                     end do

                  end do loop6

               end if

            end if

         end if

      end do

   end subroutine compute_areas

!**************************************************************

!Subroutine: poclis
!
!  Calculation of the smoothing matrices associated with each
! mesh node.
!  Required to calculate the components of E1,E2
! from the tangeantial components of E
! known on the dimensions of the triangles.
!  Calculation of the components of the tangling unit vectors.
!
!  Input variables:
!
! nuvac - PV numbers associated with the odds
! coor - coordinates of Delaunay nodes
! xlcod - Delaunay dimension length
! npoel1 - pointer of the npoel2 array
! npoel2 - number of triangles surrounding a node
! xmal1 - sum of the rates*rates surrounding a node (/det)
! xmal2 - sum of tauy*tauy surrounding a node (/det)
! xmal3 - sum of rates*tauy surrounding a node (/det)

   subroutine poclis(mesh, ncotcu, nuctfr)

      type(sll_t_triangular_mesh_2d) :: mesh
      integer, dimension(:)        :: ncotcu
      integer, dimension(:)        :: nuctfr
      logical                      :: lerr
      double precision             :: det, s1, s2, s3, x21, y21
      double precision             :: xa, ya, xb, yb
      integer                      :: nel1, nel2, indv1, indv2
      integer                      :: ic, nuctf, ind, nm1, nm2
      integer                      :: indn1, indn2, num1, num2
      integer                      :: n1, n2, ivois, nc, iel, nucti
      integer                      :: nbti, is
      integer, parameter           :: iout = 6
      integer                      :: iac, nbc

!======================================================================
! --- 1.0 --- Pointer of edge nodes -----

      do is = 1, mesh%num_nodes
         nbti = mesh%npoel1(is + 1) - mesh%npoel1(is)
         mesh%nbcov(is + 1) = nbti
         if (mesh%refs(is) /= 0) then
            mesh%nbcov(is + 1) = nbti + 1
         end if
      end do

      mesh%nbcov(1) = 0
      do is = 1, mesh%num_nodes
         mesh%nbcov(is + 1) = mesh%nbcov(is) + mesh%nbcov(is + 1)
      end do

      do is = 1, mesh%num_nodes
         nbc = mesh%nbcov(is + 1) - mesh%nbcov(is)
         iac = mesh%nbcov(is)
      end do

      nucti = 0
      nuctfr = 0

! ======================================================================
! --- 2.0 --- set edge numbers -----------------------------------

      do iel = 1, mesh%num_triangles

         do nc = 1, 3

            ivois = mesh%nvois(nc, iel)
            n1 = nc
            n2 = mod(nc, 3) + 1

            num1 = mesh%nodes(n1, iel)
            num2 = mesh%nodes(n2, iel)
            indn1 = mesh%npoel1(num1)
            indn2 = mesh%npoel1(num2)
            indv1 = mesh%nbcov(num1)
            indv2 = mesh%nbcov(num2)
            nel1 = mesh%npoel1(num1 + 1) - mesh%npoel1(num1)
            nel2 = mesh%npoel1(num2 + 1) - mesh%npoel1(num2)

            !internal edges

            if (ivois > iel) then

               nucti = nucti + 1

               !Numbers of edges with the same node

               do nm1 = 1, nel1
                  if (mesh%npoel2(indn1 + nm1) == ivois) then
                     mesh%nugcv(indv1 + nm1) = nucti
                  end if
               end do

               do nm2 = 1, nel2
                  if (mesh%npoel2(indn2 + nm2) == iel) then
                     mesh%nugcv(indv2 + nm2) = nucti
                  end if
               end do

               !number of triangles

               mesh%nuvac(1, nucti) = num1
               mesh%nuvac(2, nucti) = num2

            else if (ivois < 0) then !boundaries

               ind = -ivois
               nuctfr(ind) = nuctfr(ind) + 1
               nuctf = ncotcu(ind + 1) + nuctfr(ind)
               mesh%nugcv(indv1 + nel1 + 1) = nuctf
               mesh%nugcv(indv2 + nel2) = nuctf
               mesh%nuvac(1, nuctf) = num1
               mesh%nuvac(2, nuctf) = num2

            end if

         end do

      end do

!======================================================================
!----------- length triangles edges ------------------------

      do ic = 1, mesh%nbtcot
         xa = mesh%coord(1, mesh%nuvac(1, ic))
         ya = mesh%coord(2, mesh%nuvac(1, ic))
         xb = mesh%coord(1, mesh%nuvac(2, ic))
         yb = mesh%coord(2, mesh%nuvac(2, ic))
         mesh%xlcod(ic) = sqrt((xa - xb)*(xa - xb) + (ya - yb)*(ya - yb))
      end do

!======================================================================
!--- 4.0 --- matrices to compute electric field from potential --------

      do ic = 1, mesh%nbtcot
         n1 = mesh%nuvac(1, ic)
         n2 = mesh%nuvac(2, ic)
         x21 = (mesh%coord(1, n2) - mesh%coord(1, n1))/mesh%xlcod(ic)
         y21 = (mesh%coord(2, n2) - mesh%coord(2, n1))/mesh%xlcod(ic)
         s1 = x21*x21
         s2 = y21*y21
         s3 = x21*y21
         mesh%vtaux(ic) = x21
         mesh%vtauy(ic) = y21
         mesh%xmal1(n1) = mesh%xmal1(n1) + s1
         mesh%xmal2(n1) = mesh%xmal2(n1) + s2
         mesh%xmal3(n1) = mesh%xmal3(n1) + s3
         mesh%xmal1(n2) = mesh%xmal1(n2) + s1
         mesh%xmal2(n2) = mesh%xmal2(n2) + s2
         mesh%xmal3(n2) = mesh%xmal3(n2) + s3
      end do

!--- 4.5 --- Normalization ---------------
      lerr = .false.

      do is = 1, mesh%num_nodes
         det = mesh%xmal1(is)*mesh%xmal2(is) - mesh%xmal3(is)*mesh%xmal3(is)
         if (det /= 0) then
            mesh%xmal1(is) = mesh%xmal1(is)/det
            mesh%xmal2(is) = mesh%xmal2(is)/det
            mesh%xmal3(is) = mesh%xmal3(is)/det
         else
            lerr = .TRUE.
         end if
      end do

      if (lerr) then
         write (iout, 900)
         stop "poclis"
      end if

      if (lerr) then
         write (iout, 901)
         stop "poclis"
      end if

900   format(//10x, 'Determinant des coefficients des matrices' &
              /10x, 'de lissage nul' &
              /10x, 'Il faut modifier le maillage'//)
901   format(//10x, 'Le nombre de triangles communs a un noeud' &
              /10x, 'est superieur a 12' &
              /10x, 'Modifiez legerement le maillage'//)

   end subroutine poclis

!========================================================================
!> Compute unstructured mesh quantities
!>
!>  num_nodes - nodes
!>  num_cells - triangles
!>  coord - coord nodes
!>  refs - references
!>  nodes - number nodes
!>  nvois - number neighbors
!>  area - areas triangles
!>  base - integral bsis function
!>  nusd - num subdomains
!>
!>  npoel1 - pointer to "npoel2" array
!>  npoel2 - triangles indices with same node
!>
!>  petitl - small reference length
!>  grandl - big reference length
!>
!>  nbtcot - number of edges
!>  nbcoti - number of edges not on the boundary
!>
!> From subroutine written by A. Adolf / L. Arnaud - Octobre 1991
   subroutine sll_s_analyze_triangular_mesh(mesh)

      integer, parameter :: iout = 6
      integer            :: i, j, ifr
      type(sll_t_triangular_mesh_2d), intent(inout) :: mesh

      integer, dimension(:), allocatable :: nuctfr
      integer, dimension(:), allocatable :: ncotcu

      integer :: iel, is1, is2
      integer :: ict, ictcl, iref, jref
      integer :: keltmp, kcttmp, kretmp, ks1tmp, ks2tmp
      integer :: nbcot, ict1, ict2, ictcl1, ictcl2
      real(8) :: x1, y1, x2, y2

!-----

      write (iout, "(//10x,'>>> Compute come quantities from mesh <<<'/)")

! --- Definition de nctfrt: le nombre de cotes frontieres

      mesh%nctfrt = 0
      do i = 1, mesh%num_triangles
         if (mesh%nvois(1, i) < 0) mesh%nctfrt = mesh%nctfrt + 1
         if (mesh%nvois(2, i) < 0) mesh%nctfrt = mesh%nctfrt + 1
         if (mesh%nvois(3, i) < 0) mesh%nctfrt = mesh%nctfrt + 1
      end do

!*** total edges and internal (nbcoti,nbtcot)

      mesh%nbtcot = (3*mesh%num_triangles + mesh%nctfrt)/2

      mesh%nbcoti = mesh%nbtcot - mesh%nctfrt

!  cumsum of edges ....

      allocate (ncotcu(mesh%nbcoti + mesh%nmxfr*mesh%nctfrt))

      ncotcu = 0
      ncotcu(1) = mesh%nbcoti
      ncotcu(2) = mesh%nbcoti

      do ifr = 1, mesh%nmxfr
         do i = 1, mesh%num_triangles
            do j = 1, 3
               if (mesh%nvois(j, i) == -ifr) then
                  ncotcu(ifr + 2) = ncotcu(ifr + 2) + 1
               end if
            end do
         end do
      end do

      do ifr = 1, mesh%nmxfr
         ncotcu(ifr + 2) = ncotcu(ifr + 1) + ncotcu(ifr + 2)
      end do

      allocate (mesh%vtaux(mesh%nbtcot), mesh%vtauy(mesh%nbtcot))
      mesh%vtaux = 0._f64; mesh%vtauy = 0._f64

      allocate (mesh%nuvac(2, mesh%nbtcot)); mesh%nuvac = 0

      allocate (mesh%xlcod(mesh%nbtcot)); mesh%xlcod = 0.0_f64

      allocate (mesh%xmal1(mesh%num_nodes)); mesh%xmal1 = 0._f64
      allocate (mesh%xmal2(mesh%num_nodes)); mesh%xmal2 = 0._f64
      allocate (mesh%xmal3(mesh%num_nodes)); mesh%xmal3 = 0._f64

      allocate (mesh%nbcov(mesh%num_nodes + 1))
      allocate (mesh%nugcv(10*mesh%num_nodes))
      allocate (nuctfr(mesh%nmxfr))

      call poclis(mesh, ncotcu, nuctfr)

      deallocate (nuctfr)

!quantities for boundary edges
!
!    Numbering of border dimensions in the direction
! trigonometric, by reference number, for
! carry out all diagnoses relating to
! borders.
!
! kelfro - element number to which the dimension belongs
! kctfro - local number (1,2,3) of rating
! krefro - reference number of the dimension
! ksofro - numbers of the 2 vertices at the end of the side
! vnofro - normal vector components (inward)
!
! nctfrt - Total number of border dimensions
! nctfro - Number of border dimensions by reference
! nctfrp - Array pointer (nb of dimension by ref)

! --- Initialisation (numerotation quelconque) -----------------

      allocate (mesh%nctfrp(0:mesh%nmxfr))
      allocate (mesh%nctfro(mesh%nmxfr))
      allocate (mesh%kctfro(mesh%nctfrt))
      allocate (mesh%kelfro(mesh%nctfrt))
      allocate (mesh%krefro(mesh%nctfrt))
      allocate (mesh%ksofro(2, mesh%nctfrt))
      allocate (mesh%vnofro(2, mesh%nctfrt))

      mesh%nctfro = 0

      ifr = 0
      do ict = 1, 3
         do iel = 1, mesh%num_triangles
            if (mesh%nvois(ict, iel) < 0) then

               iref = -mesh%nvois(ict, iel)

               ifr = ifr + 1
               mesh%kelfro(ifr) = iel  !element
               mesh%kctfro(ifr) = ict  !local edge
               mesh%krefro(ifr) = iref !reference edge
               mesh%ksofro(1, ifr) = mesh%nodes(ict, iel) ! 2 edge nodes
               mesh%ksofro(2, ifr) = mesh%nodes(mod(ict, 3) + 1, iel)
               mesh%nctfro(iref) = mesh%nctfro(iref) + 1  ! edge per reference

            end if
         end do
      end do

      mesh%nctfrp(0) = 0
      do ifr = 1, mesh%nmxfr
         mesh%nctfrp(ifr) = mesh%nctfrp(ifr - 1) + mesh%nctfro(ifr)
      end do

!--- Loop over reference --------------------

      ictcl = 1
      do iref = 1, mesh%nmxfr

         nbcot = mesh%nctfro(iref)
         if (nbcot > 0) then

            ict1 = ictcl

            do ict = ict1, mesh%nctfrt
               jref = mesh%krefro(ict)

               if (jref == iref) then

                  keltmp = mesh%kelfro(ict)
                  kcttmp = mesh%kctfro(ict)
                  kretmp = mesh%krefro(ict)
                  ks1tmp = mesh%ksofro(1, ict)
                  ks2tmp = mesh%ksofro(2, ict)

                  mesh%kelfro(ict) = mesh%kelfro(ictcl)
                  mesh%kctfro(ict) = mesh%kctfro(ictcl)
                  mesh%krefro(ict) = mesh%krefro(ictcl)
                  mesh%ksofro(1, ict) = mesh%ksofro(1, ictcl)
                  mesh%ksofro(2, ict) = mesh%ksofro(2, ictcl)

                  mesh%kelfro(ictcl) = keltmp
                  mesh%kctfro(ictcl) = kcttmp
                  mesh%krefro(ictcl) = kretmp
                  mesh%ksofro(1, ictcl) = ks1tmp
                  mesh%ksofro(2, ictcl) = ks2tmp

                  ictcl = ictcl + 1

               end if
            end do
         end if
      end do

      do iref = 1, mesh%nmxfr

         nbcot = mesh%nctfro(iref)

         !One reference for the whole boundary

         if (nbcot == mesh%nctfrt) then
            ict1 = 1
35          continue
            is2 = mesh%ksofro(2, ict1)
            do ict2 = ict1, mesh%nctfrt
               if (ict1 /= ict2) then
                  is1 = mesh%ksofro(1, ict2)
                  if (is1 == is2) then

                     keltmp = mesh%kelfro(ict1 + 1)
                     kcttmp = mesh%kctfro(ict1 + 1)
                     kretmp = mesh%krefro(ict1 + 1)
                     ks1tmp = mesh%ksofro(1, ict1 + 1)
                     ks2tmp = mesh%ksofro(2, ict1 + 1)

                     mesh%kelfro(ict1 + 1) = mesh%kelfro(ict2)
                     mesh%kctfro(ict1 + 1) = mesh%kctfro(ict2)
                     mesh%krefro(ict1 + 1) = mesh%krefro(ict2)
                     mesh%ksofro(1, ict1 + 1) = mesh%ksofro(1, ict2)
                     mesh%ksofro(2, ict1 + 1) = mesh%ksofro(2, ict2)

                     mesh%kelfro(ict2) = keltmp
                     mesh%kctfro(ict2) = kcttmp
                     mesh%krefro(ict2) = kretmp
                     mesh%ksofro(1, ict2) = ks1tmp
                     mesh%ksofro(2, ict2) = ks2tmp

                     if (ict1 < mesh%nctfrt) then
                        ict1 = ict1 + 1
                        goto 35
                     end if

                  end if
               end if
            end do

            ! ... several references for one boundary

         else if (nbcot > 1) then

            ictcl1 = mesh%nctfrp(iref - 1) + 1
            ictcl2 = mesh%nctfrp(iref)

            ict1 = ictcl1
31          continue
            is1 = mesh%ksofro(1, ict1)
            do ict2 = ictcl1, ictcl2
               if (ict1 /= ict2) then
                  is2 = mesh%ksofro(2, ict2)
                  if (is1 == is2) then
                     ict1 = ict2
                     if (ict1 == ictcl1) then !Test si on tourne en rond ............
                        goto 37
                     end if
                     goto 31
                  end if
               end if
            end do

            keltmp = mesh%kelfro(ict1)
            kcttmp = mesh%kctfro(ict1)
            kretmp = mesh%krefro(ict1)
            ks1tmp = mesh%ksofro(1, ict1)
            ks2tmp = mesh%ksofro(2, ict1)

            mesh%kelfro(ict1) = mesh%kelfro(ictcl1)
            mesh%kctfro(ict1) = mesh%kctfro(ictcl1)
            mesh%krefro(ict1) = mesh%krefro(ictcl1)

            mesh%ksofro(1, ict1) = mesh%ksofro(1, ictcl1)
            mesh%ksofro(2, ict1) = mesh%ksofro(2, ictcl1)

            mesh%kelfro(ictcl1) = keltmp
            mesh%kctfro(ictcl1) = kcttmp
            mesh%krefro(ictcl1) = kretmp

            mesh%ksofro(1, ictcl1) = ks1tmp
            mesh%ksofro(2, ictcl1) = ks2tmp

37          continue
            ict1 = ictcl1
33          continue
            is2 = mesh%ksofro(2, ict1)
            do ict2 = ict1, ictcl2
               if (ict1 /= ict2) then
                  is1 = mesh%ksofro(1, ict2)
                  if (is1 == is2) then

                     keltmp = mesh%kelfro(ict1 + 1)
                     kcttmp = mesh%kctfro(ict1 + 1)
                     kretmp = mesh%krefro(ict1 + 1)
                     ks1tmp = mesh%ksofro(1, ict1 + 1)
                     ks2tmp = mesh%ksofro(2, ict1 + 1)

                     mesh%kelfro(ict1 + 1) = mesh%kelfro(ict2)
                     mesh%kctfro(ict1 + 1) = mesh%kctfro(ict2)
                     mesh%krefro(ict1 + 1) = mesh%krefro(ict2)
                     mesh%ksofro(1, ict1 + 1) = mesh%ksofro(1, ict2)
                     mesh%ksofro(2, ict1 + 1) = mesh%ksofro(2, ict2)

                     mesh%kelfro(ict2) = keltmp
                     mesh%kctfro(ict2) = kcttmp
                     mesh%krefro(ict2) = kretmp
                     mesh%ksofro(1, ict2) = ks1tmp
                     mesh%ksofro(2, ict2) = ks2tmp

                     ict1 = ict1 + 1
                     goto 33
                  end if
               end if
            end do

            if (ict1 < ictcl2) then
               ictcl1 = ict1 + 1
               ict1 = ictcl1
               goto 31
            end if
         end if
      end do

!PN commented out because it is used only for particles simulation
!!  Modification de nvois ------------------------------------
!! (Numero de reference change par le numero de cote)
!
!do ict=1,mesh%nctfrt
!   ie=mesh%kelfro(ict)
!   ic=mesh%kctfro(ict)
!   mesh%nvois(ic,ie)=-ict
!end do

!--- Calcul des composantes du vecteur normal -----------------

      do ict = 1, mesh%nctfrt

         is1 = mesh%ksofro(1, ict)
         is2 = mesh%ksofro(2, ict)

         x1 = mesh%coord(1, is1)
         y1 = mesh%coord(2, is1)
         x2 = mesh%coord(1, is2)
         y2 = mesh%coord(2, is2)

         mesh%vnofro(1, ict) = -y2 + y1
         mesh%vnofro(2, ict) = x2 - x1

      end do

      mesh%analyzed = .true.

   end subroutine sll_s_analyze_triangular_mesh

end module sll_m_triangular_meshes
