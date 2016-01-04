!> @defgroup ode_integrators sll_ode_integrators
!!
!! @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!! @authors Marco Restelli - <marco.restelli@gmail.com>
!!
!! @brief
!! ODE time integrators, with a common interface.  
!!
!! @details
!! The user should provide
!! a state vector of \c class(sll_vector_space_base), and an ODE system of
!! \c class(sll_ode_base).
!! This library defines an abstract type \c sll_ode_base for a generic ODE
!! system \f$ dy/dt = f(t,y) \f$, and an abstract type
!! \c sll_ode_integrator_base for a generic ODE time integrator.
!! Both types rely on the assumption that the state vector variable
!! \f$ y \f$ is described by an object of \c class(sll_vector_space_base).
!!
!! This library also defines a wide selection of time integrators, in the form
!! of concrete types extending \c sll_ode_integrator_base (WIP).
!!
!! @note
!! This library is derived from the *mod_time_integrators_base* module in
!! FEMilaro (Finite Element Method toolkit), by Marco Restelli.
!!
!! <h4> Header files available </h4>
!!   + No header files.
!! <!--
!!   + *sll_ode_integrators.h*
!! -->
!!
!! <h4> How to use it </h4>
!! <!--
!!  + Include header file:
!!    \code #include "sll_ode_integrators.h" \endcode
!! -->
!!  + Import Fortran modules:
!!    \code
!!      use sll_m_ode_integrators_base
!!    \endcode
!!  + Add dependency to *CMakeLists.txt*:
!!    \code
!!      target_link_libraries( <my_lib/exec> sll_ode_integrators ... )
!!    \endcode
!!
!! <h4> Examples </h4>
!! Add some fortran lines to explain how ti use the library
!! \code
!! call initialize(my_type, arg_1, arg_2)
!! call solve(my_type, your_result)
!! \endcode
!!


