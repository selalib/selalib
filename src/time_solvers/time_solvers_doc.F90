! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/doc/build/html/doxygen/html/defgroup time_solvers.html 
! The following lines will be read by doxygen to generate documentation:


!> @defgroup operator_splitting sll_operator_splitting 
!> @brief 
!> The time_solvers library implements operator splitting methods
!> @author Selalib team 
!> Contact: Eric Sonnendrucker
!> @details
!> Operator splitting is used to solve equations of the form \f$\frac{dU}{dt} = (T+V)U \f$, 
!> where \f$T \f$ and  \f$V \f$ are two differential operators by solving successively the simpler
!> equations  
!>
!> <b> References </b> 
!>
!> [HLW] E. Hairer, C. Lubich, G. Wanner, Geometrical numerical integration, Springer 2006
!> Composition of more that two operators can be reduce to the composition ot two
!> operators by taking there first order composition and its adjoint (See (HLW) section II.5)
!>
!> <b> Headers file available </b>
!>  - no header files
!>
!> <b> Modules available </b>
!>  List fortran module available
!>  - sll_time_splitting
!>  - sll_operator_splitting
!>
!> <b> How to use it </b>
!> - Header file : \code #include 'sll_time_solvers.h' \endcode
!> - Link with   <code>-lsll_%s</code>
!> - Add <code> use sll_operator_splitting </code>
!>

