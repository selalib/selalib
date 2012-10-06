!I know this file is in the wrong directory, but I couldn't find the way to have it included in the documentation.
!I just pushed it, so if others find a way to include it in the documentation, there will be a better documentation for the Poisson polary library I wrote

!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: poisson_polar
!
! DESCRIPTION:
!> @file tridiagonal.F90
!> @namespace sll_tridiagonal
!> @author Eric MADAULE
!> @brief Poisson equation solver in polar coordinate
!> @details Solver for the Poisson equation -\Delta \phi = f in polar coordinate
!!using a fft in direction \theta and final differencies in direction r.
!!This way we solve a tridiagonal system with the solver from SELALIB.
!>
!>\section how How to use the Poisson polar solver?
!>
!>You must add \code sll_poisson_2d_polar \endcode to the list of linked libraries.
!>The Poisson solver uses the FFT, so you also need to link the FFT
!>
!>1. Declare a Poisson polar plan
!>\code type(sll_plan_poisson_polar), pointer :: plan \endcode
!>2. Initialize the plan
!>\code plan => new_plan_poisson_polar(dr,rmin,nr,ntheta,boundary_conditions) \endcode
!>nr and ntheta are the number of step in direction r and theta
!> 3. Execute the plan
!> \code call poisson_solve_polar(plan,in,out) \endcode
!> 4. Delete the plan
!> \code call delete_plan_poisson_polar(plan) \endcode
!>
!>\section bc Boundary conditions :
!>
!>The boundary conditions define the potential behaviour in r_{min} (BOT_) and r_{max} (TOP_)
!>They can take the value DIRICHLET or NEUMANN
!>
!>Summary
!>The boundary conditions define the potential behaviour in r_{min} (BOT_) and r_{max} (TOP_)
!>They can take the value DIRICHLET or NEUMANN
!>
!>Summary :
!>
!>The different boundary conditions are :
!>BOT_DIRICHLET
!>BOT_NEUMANN
!>TOP_DIRICHLET
!>TOP_NEUMANN
!>
!>You must combine BOT_ and TOP_ conditions with '+'.
!>
!>\section examples Example :
!>
!>Full code :
!>\code
!>integer :: nr, ntheta
!>real :: dr, rmin
!>integer :: bc
!>type(sll_plan_poisson_polar), pointer :: plan
!>real, dimension(:,:), allocatable :: in, out
!>
!>!define all parameters
!>nr     = 100
!>ntheta = 64 !ntheta must be a power of 2
!>rmin   = 1.0
!>dr     = 0.05
!>bc     = BOT_NEUMANN+TOP_DIRICHLET
!>
!>allocate(in(nr+1,ntheta+1))
!>allocate(out(nr+1,ntheta+1))
!>
!>in= !definition of in
!>
!>!initialization of plan
!>plan => new_plan_poisson_polar(dr,rmin,nr,ntheta,bc)
!>
!>!computation of Poisson
!>call poisson_solve_polar(plan,in,out)
!>
!>!deletion of plan
!>call delete_plan_poisson_polar(plan)
!>\endcode
