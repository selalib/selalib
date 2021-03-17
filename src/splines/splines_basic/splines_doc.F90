!> @defgroup splines sll_splines
!> @brief
!> Splines computation and interpolation.
!>
!> @authors Edwin Chacon-Golcher
!> @authors Yaman Güçlü
!> @authors Katharina Kormann
!> @authors Laura S. Mendoza
!> @authors Pierre Navaro
!> @authors Benedikt Perse
!> @authors Eric Sonnendrücker
!> @authors Edoardo Zoni
!>
!> <b> How to use it </b>
!> - Link with <code> -lsll_splines </code>
!> - Add <code> use sll_m_<module_name> </code>
!>
!> @details
!>
!> #### List of Components ####
!>
!> 1. General object-oriented library for 1D and 2D tensor-product splines
!> 2. Box splines on 2D hexagonal mesh with quasi-interpolation property
!> 3. Optimized cubic splines (1D and 2D tensor-product) on uniform grid
!> 4. 1D cubic splines on non-uniform grid
!> 5. 1D quintic splines on uniform grid
!> 6. Low-level evaluation of 1D B-splines on uniform and non-uniform grids
!> 7. De Boor's implementation of 1D and 2D tensor-product splines
!>
!> #### 1. General object-oriented library for 1D and 2D tensor-product splines  ####
!>
!> - sll_m_bsplines
!> - sll_m_bsplines_base
!> - sll_m_bsplines_uniform
!> - sll_m_bsplines_non_uniform
!>
!> - sll_m_spline_1d
!> - sll_m_spline_2d
!>
!> - sll_m_spline_interpolator_1d
!> - sll_m_spline_interpolator_2d
!>
!> - sll_m_spline_matrix
!> - sll_m_spline_matrix_base
!> - sll_m_spline_matrix_dense
!> - sll_m_spline_matrix_banded
!> - sll_m_spline_matrix_periodic_banded
!>
!> ##### Usage example #####
!>
!> \code
!> use sll_m_bsplines, only: &
!>   sll_c_bsplines        , &
!>   sll_s_bsplines_new
!>
!> use sll_m_spline_2d             , only: sll_t_spline_2d
!> use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d
!>
!> ...
!>
!> class(sll_c_bsplines), allocatable :: basis_x1
!> class(sll_c_bsplines), allocatable :: basis_x2
!> type(sll_t_spline_2d)              :: spline_2d
!> type(sll_t_spline_interpolator_2d) :: interp_2d
!>
!> real(wp), allocatable :: breaks_x2(:)
!> real(wp), allocatable :: tau_x1(:)
!> real(wp), allocatable :: tau_x2(:)
!> real(wp), allocatable :: gtau(:,:)
!>
!> ...
!>
!> ! Create spline bases along x1 and x2 (x1 grid is uniform, x2 grid is non-uniform)
!> call sll_s_bsplines_new( basis_x1, degree=7, periodic=.false., xmin=0.0_wp, xmax=1.0_wp, ncells=100 )
!> call sll_s_bsplines_new( basis_x2, degree=4, periodic=.true. , xmin=0.0_wp, xmax=1.0_wp, ncells= 40, breaks_x2 )
!>
!> ! Initialize 2D spline as object of 2D tensor-product space
!> call spline_2d % init( basis_x1, basis_x2 )
!>
!> ! Initialize 2D spline interpolator with bases and boundary conditions
!> call interp_2d % init( &
!>   basis_x1, &
!>   basis_x2, &
!>   bc_xmin = [sll_p_greville, sll_p_periodic], &
!>   bc_xmax = [sll_p_greville, sll_p_periodic] )
!>
!> ! Get x1 and x2 coordinates of interpolation points from interpolator (Greville averaging)
!> call interp_2d % get_interp_points( tau_x1, tau_x2 )
!> allocate( gtau( size(tau_x1), size(tau_x2) ) )
!>
!> ...
!>
!> ! Compute 2D spline that interpolates solution gtau(:,:) at interpolation points
!> call interp_2d % compute_interpolant( spline_2d, gtau )
!>
!> ! Test: evaluate spline at some location (x1,x2)
!> write(*,*) spline_2d % eval( x1=0.5_wp, x2=0.5_wp )
!>
!> ! Test: print spline coefficient c_{1,1}
!> write(*,*) spline_2d % bcoef(1,1)
!> \endcode
!>
!> #### 2. Box splines on 2D hexagonal mesh with quasi-interpolation property ####
!>
!> - sll_m_box_splines
!> - sll_m_hex_pre_filters
!>
!> #### 3. Optimized cubic splines (1D and 2D tensor-product) on uniform grid ####
!>
!> - sll_m_cubic_splines
!>
!> #### 4. 1D cubic splines on non-uniform grid ####
!>
!> - sll_m_cubic_non_uniform_splines
!>
!> #### 5. 1D quintic splines on uniform grid ####
!>
!> - sll_m_quintic_splines
!>
!> #### 6. Low-level evaluation of 1D B-splines on uniform and non-uniform grids ####
!>
!> - sll_m_low_level_bsplines
!>
