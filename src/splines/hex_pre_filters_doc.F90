! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory
! in file Doxyfile.in (line 691) if it is excluded.
! Type 'make doc' in build directory.
! To check the results, open :
! selalib/doc/build/html/doxygen/html/
! The following lines will be read by doxygen to generate documentation:


!> @defgroup meshes hex_pre_filters
!> @brief
!> Box splines library: splines that can defined on a hexagonal mesh
!> (see hex_meshes).
!> @author Laura Mendoza
!> @details
!>
!> <b> How to use it </b>
!> - Link with   <code>-lhex_pre_filters</code>
!> - Add <code> use hex_pre_filters </code>
!>
!> <b> Examples </b>
!> @snippet meshes/test_box_splines_deriv.F90 example
