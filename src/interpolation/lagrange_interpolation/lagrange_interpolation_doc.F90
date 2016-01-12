! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/doc/build/html/doxygen/html/defgroup lagrange_interpolations.html 
! The following lines will be read by doxygen to generate documentation:


!> @defgroup lagrange_interpolation sll_lagrange_interpolation 
!> @brief 
!> Description of sll_m_lagrange_interpolation library (72 characters)
!> @author Raphael Blanchard, Klaus Reuter, Katharina Kormann;
!> Contact: Katharina Kormann
!> @details
!> The library sll_lagrange_interpolation currently only offers two implementations of the interpolator routine \a interpolate_array_disp on uniform 1D grids. 
!> The module sll_m_lagrange_fast is a fast implementation that also works in combination with domain decomposition. It implements a Lagrange interpolation that is always centered around the point that is displaced (odd number of points).
!> The module sll_m_lagrange_interpolation_1d, on the other hand, uses an interpolation that is centered around interval to which the point is displaced (even number of points).
!> @todo
!> - Improve efficiency of sll_m_lagrange_interpolation_1d.F90
!> - Implement other interpolation routines.

