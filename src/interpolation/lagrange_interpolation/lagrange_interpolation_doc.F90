!> @defgroup lagrange_interpolation sll_lagrange_interpolation
!> @brief
!> Lagrange interpolation
!> @author Raphael Blanchard, Klaus Reuter, Katharina Kormann;
!> Contact: Katharina Kormann
!> @details
!> The library sll_lagrange_interpolation currently only offers two
!> implementations of the interpolator routine \a interpolate_array_disp
!> on uniform 1D grids.  The module sll_m_lagrange_fast is a fast
!> implementation that also works in combination with domain decomposition.
!> In particular, for the MPI domain decomposition use case, it
!> implements a Lagrange interpolation that is always centered around
!> the point that is displaced (odd number of points).  
!> The module
!> sll_m_lagrange_interpolation_1d, on the other hand, uses an
!> interpolation that is centered around interval to which the point
!> is displaced (even number of points).  For the boundary closure,
!> three possiblities are implemented: Either halo cells are provided
!> or periodic boundary conditions can be used. A third alternative
!> is to use a one-sided stencil at the boundary.  
!> @todo 
!> - Improve efficiency of sll_m_lagrange_interpolation_1d.F90 
!> - Implement other interpolation routines.

