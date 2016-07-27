!> @defgroup splines sll_splines
!> @brief
!> Splines computation and interpolation.
!> @authors Edwin Chacon-Golcher, Pierre Navaro and Laura S. Mendoza.
!> @details
!> Low level modules for sll_interpolators.
!> Library to use splines, contains:
!>    - sll_m_cubic_splines            : B-splines of degree 3
!>    - sll_m_cubic_non_uniform_splines: non-uniform B-splines of degree 3
!>    - sll_m_quintic_splines          : B-splines of degree 5
!>    - sll_m_bsplines                 : B-splines 1D and 2D
!>    - sll_m_arbitrary_degree_splines : De Boor splines of arbitrary degree
!>    - sll_m_box_splines              : Box-splines for hexagonal mesh
!>    - sll_m_hex_pre_filters          : Hexagonal prefilters associated to boxsplines
!>
!> <b> How to use it </b>
!> - Link with   <code>-lsll_splines</code>
!> - Add <code> use sll_m_<module_name> </code>
!>
!> <b> Examples </b>
!> @snippet see splines/testing directory
!>
