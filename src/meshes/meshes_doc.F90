! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory
! in file Doxyfile.in (line 691) if it is excluded.
! Type 'make doc' in build directory.
! To check the results, open :
! selalib/doc/build/html/doxygen/html/
! The following lines will be read by doxygen to generate documentation:


!> @defgroup meshes sll_meshes
!> @brief
!> Meshes types: cartisian, triangular and hexagonal meshes.
!> @authors Edwin Chacon-Golcher, Pierre Navaro, Aurore Back and Laura Mendoza.
!> @details
!> Cartesian meshes: logical orthogonal regular meshes.
!> Hexagonal meshes: three directional meshes discretized by equilateral
!> triangles.
!> Triangular meshes: triangular meshes.
!>
!> <b> How to use it </b>
!> - Link with   <code>-lsll_meshes</code>
!> - Add <code> use sll_meshes </code>
!>
!> <b> Examples </b>
!> @snippet See meshes/testing/unit_test.F90 example
