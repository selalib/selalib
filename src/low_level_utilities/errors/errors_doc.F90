!> @defgroup errors sll_errors
!!
!! @author Yaman Güçlü - <yaman.guclu@gmail.com>
!!
!! @brief
!! Provides macros for warnings and errors.
!!
!! @details
!! - `SLL_ERROR(fun,msg)` 
!! - `SLL_WARNING(fun,msg)`
!!
!! <h4> Header files available </h4>
!!  + *sll_m_errors.h*
!!
!! <h4> How to use it </h4>
!!  + Include header file:
!!    \code #include "sll_m_errors.h" \endcode
!!  + Add dependency to *CMakeLists.txt*:
!!    \code
!!      target_link_libraries( <my_lib/exec> sll_m_errors ... )
!!    \endcode
!!
!! <h4> Examples </h4>
!! \include test_crash.F90
