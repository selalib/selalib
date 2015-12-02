! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory
! in file Doxyfile.in (line 691) if it is excluded.
! Type 'make doc' in build directory.
! To check the results, open :
! selalib/doc/build/html/doxygen/html/
! The following lines will be read by doxygen to generate documentation:


!> @defgroup splines sll_splines
!> @brief
!> Library to use splines.
!> @authors Edwin Chacon-Golcher, Pierre Navaro and Laura S. Mendoza.
!> @details
!> Low level modules for sll_interpolators.
!> Library to use splines, contains:
!>    - cubic_splines            : B-splines of degree 3
!>    - sll_m_cubic_non_uniform_splines: non-uniform B-splines of degree 3
!>    - quintic_splines          : B-splines of degree 5
!>    - bsplines                 : B-splines 1D and 2D
!>    - arbitrary_degree_splines : De Boor splines of arbitrary degree
!>    - box_splines              : Box-splines for hexagonal mesh
!>    - hex_pre_filters          : Hexagonal prefilters associated to boxsplines
!>
!> <b> How to use it </b>
!> - Link with   <code>-lsll_splines</code>
!> - Add <code> use sll_<spline_moudle_name> </code>
!>
!> <b> Examples </b>
!> @snippet see splines/testing directory
!>
!> Some fortran modules with deboor prefix are designed to
!> compute Bsplines of arbitrary degree using deboor algorithm.
!> @author Aurore Back, Pierre Navaro.
!> @details
!> Original F77 files are available
!> on Carl de Boor webpage http://pages.cs.wisc.edu/~deboor/
!>
