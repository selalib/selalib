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
!> The sll_t_operator_splitting library implements operator splitting methods
!> @author Selalib team <br>
!> Contact: Eric Sonnendr&uuml;cker
!> @details
!> Operator splitting is used to solve equations of the form \f$\frac{dU}{dt} = (T+V)U \f$, 
!> where \f$T \f$ and  \f$V \f$ are two differential operators by solving successively the simpler
!> equations  
!>
!> The base class in sll_operator_splitting_base.F90 implements the composition form
!> different kinds of composition methods defined by their coefficients.
!>
!> The application of an operator splitting method to a concrete problem is done
!> by extending the sll_t_operator_splitting splitting base class by a new type
!> containing on the one hand the data on which the operators act
!> and a specific implementation of the two operators
!>
!> <b>Examples of applications are provided for </b>
!>  - The linear pendulum problem: the type definition is in sll_m_linear_pendulum_operators.F90
!>    and the main program performing the unit tests in unit_test_linear_pendulum.F90
!>  - The Vlasov equation with constant coefficients advection field:
!> type definition in sll_const_coef_adv_2d.F90 and main program in unit_test_const_coef_adv_2d.F90
!>  - The non linear Vlasov-Poisson equations in cartesian coordinates   
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
! !> <b> Modules available </b>
! !>  List fortran module available
! !>  - sll_m_time_splitting
! !>  - sll_m_operator_splitting
!>
!> <b> How to use it </b>
! !> - Header file : \code #include 'sll_time_solvers.h' \endcode
!> - Link with   <code>-lsll_%s</code>
!> - Add <code> use sll_m_operator_splitting </code>
!>

