!> @defgroup vector_space sll_vector_space
!!
!! @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!! @authors Marco Restelli - <marco.restelli@gmail.com>
!!
!! @brief
!! Wrappers for complex objects with addition/multiplication/copy.
!! Real arrays already provided.
!!
!! @details
!! This library defines an abstract type \c sll_vector_space_base which
!! can be used as a building block for the implementation of linear and
!! non-linear solvers, time integrators and so on.
!! The main purpose of this type is decoupling the details of the "vector"
!! object representation, such as internal fields and arrays with an
!! arbitrary number of dimensions, from the implementation of some general
!! purpose algorithm.
!!
!! The operators defined for \c sll_vector_space_base are essentially those
!! required by a <a href="http://en.wikipedia.org/wiki/Vector_space"
!! target="_blank"> vector space</a>, optionally including an inner product
!! and norm.
!!
!! In a general situation, the user should extend the abstract type
!! \c sll_vector_space_base and provide the necessary operators.
!! Nevertheless, in many cases of interests the working vector space simply
!! consists of a single multidimensional array.
!! For this situation we provide the derived types \c sll_vector_space_real_1d,
!! \c sll_vector_space_real_2d and \c sll_vector_space_real_3d, which easily
!! wrap an existing 1D/2D/3D array through pointer association.
!!
!! @note
!! This library is derived from the *mod_state_vars* module in FEMilaro
!! (Finite Element Method toolkit), by Marco Restelli.
!!
!! <h4> Header files available </h4>
!!   + No header files.
!! <!--
!!   + *sll_vector_space.h*
!! -->
!!
!! <h4> How to use it </h4>
!! <!--
!!  + Include header file:
!!    \code #include "sll_vector_space.h" \endcode
!! -->
!!  + Import Fortran modules:
!!    \code
!!      use sll_m_vector_space_base
!!      use sll_m_vector_space_real_arrays
!!    \endcode
!!  + Add dependency to *CMakeLists.txt*:
!!    \code
!!      target_link_libraries( <my_lib/exec> sll_vector_space ... )
!!    \endcode
!!
!! <h4> Examples </h4>
!! Add some fortran lines to explain how ti use the library
!! \code
!! call initialize(my_type, arg_1, arg_2)
!! call solve(my_type, your_result)
!! \endcode
!!

