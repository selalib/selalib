!> @defgroup advection sll_advection 
!> @brief 
!> Solve advection equation
!> @authors Michel Mehrenberger, Pierre Navaro
!
!> @details
!> 
!>  sll_advection_1d_base.F90
!>  -> new suggested name: sll_c_advection_1d.F90
!>  -> solves advection equation in 1d
!>    sll_advection_1d_BSL.F90 
!>    -> new suggested name: sll_m_advection_1d_BSL.F90
!>    -> standard BSL method
!>    -> which needs an interpolator
!>    -> and a characteristics solver
!>    sll_advection_1d_CSL.F90
!>    -> CSL method (not finished, Michel) 
!>    sll_advection_1d_CSL_periodic.F90
!>    -> (not finished, Michel)
!>    sll_advection_1d_PSM.F90
!>    -> (not finished, Michel)
!>    sll_advection_1d_ampere.F90
!>    -> developped by Pierre
!>    sll_advection_1d_non_uniform_cubic_splines.F90
!>    -> new suggested name: sll_m_advection_1d_non_uniform_cubic_splines.F90
!>    -> only constant advection
!>    -> with non uniform cubic splines
!>    -> used for example in KEEN wave simulation
!>    -> when two grid method is used
!>
!>    sll_advection_1d_periodic.F90
!>    -> new suggested name: sll_m_advection_1d_periodic.F90
!>    -> only constant advection
!>    -> periodic interpolation
!>    -> different choices are available:
!>    -> lagrange of odd order arbitrary
!>    -> order = 4 corresponds to Lagrange of degree 3
!>    -> order = 18 corresponds to Lagrange of degree 17
!>    -> splines of arbitrary order
!>   -> order = 4 corresponds to cubic splines
!>    
!>    sll_advection_1d_spectral.F90
!>    -> developed by Pierre
!>  sll_advection_2d_base.F90
!>  -> new suggested name: sll_c_advection_2d.F90  
!>    sll_advection_2d_BSL.F90
!>    -> new suggested name: sll_m_advection_2d_bsl.F90
!>    -> 2d extension of advection_1d_BSL    
!>    sll_advection_2d_tensor_product.F90
!>    -> new suggested name: sll_m_advection_2d_tensor_product.F90
!>    -> 2D advection using succession
!>    -> of 1D advection and STRANG splitting
!>    -> that is 1D advection in first direction over dt/2
!>    -> then 1D advection in second direction over dt
!>    -> finally 1D advection in first direction over dt/2
!>    -> initialized with 1D advections
!>     
!>    sll_advection_2d_CSL.F90
!>    -> not finished (Michel)
!>    sll_advection_2d_integer_oblic.F90
!>    -> not finished (Michel)    
!>    sll_advection_2d_oblic.F90
!>    -> not finished (Michel)
!>  
!>  other modules:
!>  sll_m_advection_2d_tri_mesh.F90
!>  -> new suggested name: sll_m_advection_2d_tri_mesh.F90
!>  -> developed by Pierre
!>  
!>  sll_m_lagrange_interpolation.F90
!>  -> see redundancy with other such files.
!>  
!>  programs:
!>  aligned_translation_2d.F90
!>  aligned_translation_2d_spaghetti.F90
!>  rotation_2d.F90
!>  -> may be put in simulations
!>  -> or can still live here
!>  -> as basic examples    
!>
