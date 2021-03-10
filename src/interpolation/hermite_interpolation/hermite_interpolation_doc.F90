! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory
! in file Doxyfile.in (line 691) if it is excluded.
! Type 'make doc' in build directory.
! To check the results, open :
! selalib/doc/build/html/doxygen/html/namespaces.html
! The following lines will be read by doxygen to generate documentation:

!> @namespace sll_hermite_interpolation
!> @brief
!> Description of hermite_interpolation library (72 characters)
!> @author Michel Mehrenberger Charles Prouveur, Laura Mendoza
!> Christophe Steiner
!> You can add a contact, do not put your email to prevent spam.
!> @details
!> Long description of  hermite_interpolation, you can add some references or math equations.
!> sll_hermite_interpolation_1d.F90
!> -> new suggested name: sll_m_hermite_interpolation_1d.F90
!> -> Hermite interpolation in 1d (as indicated!)
!> sll_hermite_interpolation_2d.F90
!> -> new suggested name: sll_m_hermite_interpolation_2d.F90
!> -> Hermite interpolation in 2d (as indicated!)
!> -> reference: A. Hamiaz, M. Mehrenberger, Sellama H., and E. Sonnendrucker. The semi-lagrangian method on curvilinear grids,
!> -> 2015,submitted.
!> -> C. Steiner
!> -> Thesis 2014, Chapter 7
!>
!> sll_hermite_aligned_interpolation_2d.F90
!> not finished, Michel
!> sll_hex_interpolation_hermite.F90
!> interpolation on hexagonal mesh
!> in progress
!> stuff developed in simulations
!> redundandly should be added here
!> (Charles, Michel, Laura)
!> sll_hex_interpolation_hermite_with_hole.F90
!> interpolation on hexagonal mesh with a hole
!> (Charles) not sure that it will remain
!> other files
!> programs that test the method
!> for rotation and guiding center
!> reference:
!> M. Mehrenberger, L. Mendoza, C. Prouveur,
!> E. Sonnendrücker,
!> Solving the guiding-center model on a regular
!> hexagonal mesh, ESAIM Proceedings, CEMRACS 2014, submitted

!> <b> Headers file available </b>
!>  - sll_hermite_interpolation.h
!>
!> <b> Modules available </b>
!>  List fortran module available
!>  - sll_hermite_interpolation
!>
!> <b> How to use it </b>
!> - Header file : \code #include 'sll_hermite_interpolation.h' \endcode
!> - Link with   <code>-lsll_%s</code>
!> - Add <code> use sll_hermite_interpolation </code>
!>
!> <b> Examples </b>
!> -Add some fortran lines to explain how ti use the library
!> \code
!> call initialize(my_type, arg_1, arg_2)
!> call solve(my_type, your_result)
!> \endcode
!>
