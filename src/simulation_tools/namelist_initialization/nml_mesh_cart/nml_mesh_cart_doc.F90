!> @defgroup nml_mesh_cart sll_nml_mesh_cart
!> @brief
!> Library to deal with initialization of cartesian mesh
!> from namelist
!> <br>
!> @details
!> We propose a uniform way of initializing cartesian meshes
!> from namelists
!> <br>
!> <br>
!> FIRST EXAMPLE
!> <br>
!> in namelist test.nml
!> <br>
!> <code>
!> &mesh_1d_unif_cart
!> <br>
!>  num_cells = 32
!> <br>
!>  eta_min = 0.
!> <br>
!>  eta_max = 2.
!> <br>
!> /
!> <br>
!> </code>
!> \code
!>!ADD_EXECUTABLE(my_prog my_prog.f90)
!>!TARGET_LINK_LIBRARIES(my_prog sll_nml_mesh_cart)
!>program my_prog
!>use sll_m_cartesian_meshes, only : &
!>    sll_cartesian_mesh_1d
!>use sll_m_nml_mesh_1d_unif_cart, only : &
!>    sll_s_nml_mesh_1d_unif_cart
!>  implicit none
!>
!>  type(sll_cartesian_mesh_1d), pointer :: mesh
!>  call sll_s_nml_mesh_1d_unif_cart( "test", mesh )
!>  print *,'#mesh%num_cells=',mesh%num_cells
!>end program
!> \endcode
!> Information from mesh_1d_unif_cart in namelist file test.nml is stored 
!> in <code>mesh</code> 
!> <br>
!> <br>
!> SECOND EXAMPLE
!> <br>
!> in namelist test.nml, we can  also have
!> another cartesian mesh corresponding for example to the
!> second dimension
!> <br>
!> <code>
!> &mesh_1d_unif_cart_2
!> <br>
!>  num_cells_2 = 32
!> <br>
!>  eta_min_2 = 0.
!> <br>
!>  eta_max_2 = 2.
!> <br>
!> /
!> <br>
!> </code>
!> \code
!>!ADD_EXECUTABLE(my_prog my_prog.f90)
!>!TARGET_LINK_LIBRARIES(my_prog sll_nml_mesh_cart)
!>program my_prog
!>use sll_m_cartesian_meshes, only : &
!>    sll_cartesian_mesh_1d
!>use sll_m_nml_mesh_1d_unif_cart, only : &
!>    sll_s_nml_mesh_1d_unif_cart
!>  implicit none
!>
!>  type(sll_cartesian_mesh_1d), pointer :: mesh
!>  call sll_s_nml_mesh_1d_unif_cart( "test", mesh, clone="_2" )
!>  print *,'#mesh%num_cells=',mesh%num_cells
!>end program
!> \endcode
!> Information from mesh_1d_unif_cart_2 in namelist file test.nml is stored in 
!> <code> mesh </code> 
!> <br>
!> <br>
!> THIRD EXAMPLE
!> <br>
!> We can also choose the mesh we want to initialize
!> when we want to initialize an array
!> (not a <code>sll_cartesian_mesh_1d</code> which is uniform) 
!> <br>
!> Suppose that we have the following namelist file test.nml
!> <br>
!> <br>
!> \code
!> &mesh_1d_cart_1
!>  choice_1 = "landau" 
!> /
!> &mesh_1d_landau_cart_1
!>  num_cells_1 = 32
!>  eta_min_1 = 0.
!>  nbox_1 = 1
!>  kmode_1 = 0.5
!> /
!> 
!> &mesh_1d_cart_2
!>  choice_2 = "unif" 
!> /
!> &mesh_1d_unif_cart_2
!>  num_cells_2 = 32
!>  eta_min_2 = -6.
!>  eta_max_2 = 6.
!> /
!>  
!> \endcode
!> <br>
!> and the following code
!> \code
!>!ADD_EXECUTABLE(my_prog my_prog.f90)
!>!TARGET_LINK_LIBRARIES(my_prog sll_nml_mesh_cart)
!>program my_prog
!>use sll_m_nml_mesh_1d_cart, only : &
!>    sll_s_nml_mesh_1d_cart
!>  implicit none
!>  sll_real64, pointer :: x(:)
!>  sll_real64, pointer :: v(:)
!>  call sll_s_nml_mesh_1d_cart( "test", x, clone="_1" )
!>  call sll_s_nml_mesh_1d_cart( "test", v, clone="_2" )
!>  print *,'#num_cells=',size(x)-1,size(v)-1
!>end program
!> \endcode
!> Information from namelist file test.nml is stored in 
!> <code> x </code> and <code> v </code>.
!> @author Michel Mehrenberger
!> @todo add initialization in 2D,3D,4D if needed
