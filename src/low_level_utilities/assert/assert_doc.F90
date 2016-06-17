!> @defgroup assert sll_assert
!!
!! @authors
!! Edwin Chacon-Golcher <br>
!! Contact: Yaman Güçlü - <yaman.guclu@gmail.com>
!!
!! @brief
!! Provides macro for assertion
!! @details
!! - `SLL_ASSERT(x)`
!!
!! <h4> Header files available </h4>
!!  + *sll_m_assert.h*
!!
!! <h4> How to use it </h4>
!!  + Include header file:
!!    \code #include "sll_m_assert.h" \endcode
!!  + Add dependency to *CMakeLists.txt*:
!!    \code
!!      target_link_libraries( <my_lib/exec> sll_m_assert ... )
!!    \endcode
!!
!! <h4> Examples </h4>
!! \code
!! SLL_ASSERT( b /= 0 )
!! c = a / b
!! \endcode
