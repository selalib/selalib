!> @defgroup reduction sll_reduction 
!> @brief 
!> Compute averages in some directions
!> @author Michel Mehrenberger 
!> @details
!> Generic function which can be used for computing charge density.
!> By default we suppose a uniform mesh node based.
!> Trapezoid quadrature formula is used.
!> We can also integrate by providing a function whose signature is 
!> sll_integration_discrete_1d.
