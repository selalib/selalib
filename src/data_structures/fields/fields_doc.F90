! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/doc/build/html/doxygen/html/defgroup fieldss.html 
! The following lines will be read by doxygen to generate documentation:


!> @defgroup fields sll_fields 
!> @brief 
!> Field objects.
!> @author Edwin Chacon-Golcher and Pierre Navaro.
!> @details
!>
!> <b> How to use it </b>
!> - Link with   <code>-lsll_fields</code>
!> - Add <code> use sll_fields </code>
!>
!> <b> Examples </b>
!> -Add some fortran lines to explain how ti use the library
!> \code
!> type(sll_cartesian_mesh_2d),                       pointer :: mesh_2d
!> class(sll_coordinate_transformation_2d_base),      pointer :: tau
!> class(sll_scalar_field_2d_base),                   pointer :: rho
!> 
!> 
!> mesh_2d => new_cartesian_mesh_2d( NUM_CELLS1, &
!>                                   NUM_CELLS2, &
!>                                   ETA1MIN,    &
!>                                   ETA1MAX,    &
!>                                   ETA2MIN,    &
!>                                   ETA2MAX )
!> 
!> tau => new_coordinate_transformation_2d_analytic( &
!> &      "analytic",                                &
!> &      mesh_2d,                                   &
!> &      identity_x1,                               &
!> &      identity_x2,                               &
!> &      identity_jac11,                            &
!> &      identity_jac12,                            &
!> &      identity_jac21,                            &
!> &      identity_jac22,                            &
!> &      [0.0_f64]       
!> 
!> 
!> rho => new_scalar_field_2d_analytic( &
!> &    rhs,                            &
!> &    "rho",                          &     
!> &    tau,                            &
!> &    SLL_DIRICHLET,                  &
!> &    SLL_DIRICHLET,                  &
!> &    SLL_DIRICHLET,                  &
!> &    SLL_DIRICHLET,                  &
!> &    [0.0_f64]                       )
!> 
!> 
!> do j=1,npts2
!> do i=1,npts1
!>    eta1  = (i-1)*mesh_2d%delta_eta1 + mesh2d%eta1_min
!>     eta2  = (j-1)*mesh_2d%delta_eta2 + mesh_2d%eta2_min
!> 
!> 
!>   point_val        = rho%value_at_point(eta1,eta2)
!>   node_val        = rho%value_at_indices(i,j)
!>   grad1_node_val  = rho%first_deriv_eta1_value_at_point(eta1, eta2)
!>   grad2_node_val  =rho%first_deriv_eta2_value_at_point(eta1, eta2)
!> 
!> end do
!> end do
!> 
!> CONTAINS
!> 
!> function rhs( eta1, eta2, params ) result(res)
!> real(8), intent(in) :: eta1
!> real(8), intent(in) :: eta2
!> real(8), dimension(:), intent(in) :: params
!> real(8) :: res
!> real(8) :: pi
!> pi = 4.0_f64*atan(1.0_f64)
!> res = -8*pi*pi*sin(2*pi*eta1)*sin(2*pi*eta2)
!> 
!> end function rhs
!> 
!> end program
!> \endcode
!>
