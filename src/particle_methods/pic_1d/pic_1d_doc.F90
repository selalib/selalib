!> @defgroup pic_1d sll_pic_1d
!!
!! @authors Jakob Ameres   - <jakob.ameres@tum.de>
!! @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!!
!! @brief
!! Object-oriented particle-in-cell (PIC) code for the solution of the 1D-1V
!! Vlasov-Poisson equation in a periodic domain (electrons + uniform ions).
!!
!! @details
!! This is an attempt at restructuring the original code by Jakob.
!! The code is split into library side (this directory) + simulation side.
!!
!! Here (library side) we have:
!!   * a wrapper for different field solvers
!!   * particle description
!!   * functions for particle loading
!!   * reconstruction of the distribution function
!!   * other post-processing utilities
!!
!! On the simulation side we have:
!!   * A simulation type is defined
!!   * Input parameters are read from .nml files rather than hard-coded
!!   * It is easy to switch between different algorithms (FEM/FD, various
!!     time integrators, order of splines, etc.)
!! 
!! @todo
!! * Time integrators from the sll_ode_integrators library should be used.
!! * The PIC abstract class should be extended
!!
!! <h4> How to use it </h4>
!!  + Import Fortran modules:
!!    \code
!!      use sll_pic_1d
!!    \endcode
!!  + Add dependency to *CMakeLists.txt*:
!!    \code
!!      target_link_libraries( <my_lib/exec> sll_pic_1d ... )
!!    \endcode
!!
!! <h4> Examples </h4>
!! @todo
!! Add some fortran lines to explain how to use the library

