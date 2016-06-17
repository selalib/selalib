!> @defgroup operator_splitting sll_operator_splitting 
!> @brief 
!> Operator splitting methods
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
!> <b> Modules available </b>
!>  - sll_m_time_splitting
!>  - sll_m_operator_splitting

