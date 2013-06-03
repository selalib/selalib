program test_general_qns
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_general_coordinate_qn_solver_module
  use sll_module_scalar_field_2d_alternative
  use sll_constants
  use sll_arbitrary_degree_spline_interpolator_2d_module
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none

#define SPLINE_DEG1 1
#define SPLINE_DEG2 1
#define NUM_CELLS1  64
#define NUM_CELLS2  64
#define ETA1MIN  0.0_f64
#define ETA1MAX  1.0_f64
#define ETA2MIN  0.0_f64
#define ETA2MAX  1.0_f64

  type(sll_logical_mesh_2d), pointer          :: mesh_2d
  class(sll_coordinate_transformation_2d_base), pointer :: T
  type(general_coordinate_qn_solver), pointer :: qns
  type(arb_deg_2d_interpolator), target :: interp_2d
  class(sll_interpolator_2d_base), pointer :: interp_2d_ptr
!  class(sll_scalar_field_2d_analytic_alt), dimension(2,2) :: a_field_mat
  type(sll_scalar_field_2d_base_ptr), dimension(2,2) :: a_field_mat
  class(sll_scalar_field_2d_base), pointer    :: c_field
  class(sll_scalar_field_2d_base), pointer    :: rho
  type(sll_scalar_field_2d_discrete_alt), pointer    :: phi
  real(8), external :: func_zero, func_one, source_term_perper
  sll_real64, dimension(:,:), allocatable :: values
  sll_real64 :: acc
  sll_real64, dimension(:,:), allocatable    :: calculated
  sll_real64, dimension(:,:), allocatable    :: difference
  sll_int32 :: ierr
  sll_int32  :: i, j
  sll_real64 :: h1,h2,eta1,eta2,node_val,ref
  sll_int32 :: npts1,npts2
  real(8), external :: sol_exacte_perper
  
  



  !--------------------------------------------------------------------
  
  !     first test case without chane of coordinates 
  !      periodic-periodic boundary conditions

  !--------------------------------------------------------------------
  
  print*, "---------------------"
  print*, "first test case witout change of coordinates"
  print*, "---------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = 1.0_f64/real(NPTS1-1,f64)
  h2 = 1.0_f64/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2

  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  
  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       0.0_f64, 1.0_f64, 0.0_f64,1.0_f64 )

  ! Second, initialize the coordinate transformation associated with this 
  ! problem.
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22 )
  print *, 'initialized coordinate transformation'

  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a_field_mat(1,1)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 

  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 

  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 

  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 


  c_field => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_perper, &
       "rho", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SPLINE_DEG1, &
       SPLINE_DEG2 )

  interp_2d_ptr => interp_2d

  phi => new_scalar_field_2d_discrete_alt( &
       values, &
       "phi", &
       interp_2d_ptr, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  print *, 'initialized fields...'

  qns => new_general_qn_solver( &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)

  print *, 'Initialized QNS object'

  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )

  print *, 'Completed solution'

  call  set_coefficients_ad2d( &
       interp_2d, &
       qns%phi_vec)

  print*, 'reorganizaton of splines coefficients of solution'
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc = 0.0_f64

  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1
        eta2       = real(j,f64)*h2
        node_val   = interp_2d%interpolate_value(eta1,eta2)
        print*, 'coucou'
        ref        = sol_exacte_perper(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
             'theoretical = ', ref
        acc        = acc + abs(node_val-ref)
     end do
  end do
  

  print *,'Average error in nodes (per-per) without change of coordinates='&
       ,acc/(npts1*npts2)


  ! delete things...
  call delete(qns)
  print *, 'PASSED'
end program test_general_qns

function func_one( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in), optional :: params
  real(8) :: res
  res = 1.0_8
end function func_one

function func_zero( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in), optional :: params
  real(8) :: res
  res = 0.0_8
end function func_zero


!----------------------------------------------------------
!  Solution for a identity change of coordinates 
!   and periodic-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------
function source_term_perper( eta1, eta2, params ) result(res)
  use sll_constants
  intrinsic :: cos

  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in), optional :: params
  real(8) :: res
  res = (2.0*sll_pi)**2*cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
end function source_term_perper

real(8) function sol_exacte_perper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper = cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
end function sol_exacte_perper


!----------------------------------------------------------
!  Solution for a identity change of coordinates 
!   and periodic-dirichlet conditions
!   and also dirivhlet-dirichlet conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_perdir(eta1,eta2,params) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8), dimension(:), intent(in), optional :: params
  source_term_perdir = &
      (16.0*sll_pi**2*eta2**4 &
      - 16.0*sll_pi**2*eta2**2 &
      - 12.0*eta2**2 + 2.0)*sin(2*sll_pi*eta1)*cos(2*sll_pi*eta1)
end function source_term_perdir


real(8) function sol_exacte_perdir(eta1,eta2,params)
  use sll_constants
  real(8) :: eta1,eta2
  real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir = eta2 ** 2 * (eta2**2-1)&
       * cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta1)
  

end function sol_exacte_perdir



!----------------------------------------------------------
!  Solution for a r theta change of coordinates 
!   and periodic-dirichlet conditions
!   and also dirivhlet-dirichlet conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------


real(8) function rho_rtheta(eta1,eta2,params) ! in the path
  use sll_constants
  intrinsic :: cos
  intrinsic :: sin
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  real(8), dimension(:), intent(in), optional :: params

  x = eta2*cos(2*sll_pi*eta1)
  y = eta2*sin(2*sll_pi*eta1)
  
  rho_rtheta = x*y*(-32.0*x**2 - 32.0*y**2 + 15.0)  
  
    
end function rho_rtheta


real(8) function sol_exacte_rtheta(eta1,eta2,params) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  intrinsic :: cos
  intrinsic :: sin
  real(8), dimension(:), intent(in), optional :: params
  
  
  sol_exacte_rtheta = ( eta2**2-1)*(eta2**2-0.5**2)*eta2**2&
       *cos(2*sll_pi*eta1)*sin(2*sll_pi*eta1)
  
    
end function sol_exacte_rtheta


!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and periodic-periodic conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------




real(8) function source_term_chgt_perper(eta1,eta2,params) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  intrinsic :: cos
  intrinsic :: sin
  real(8), dimension(:), intent(in), optional :: params
  
  x =   eta1 + 0.1_8*sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1_8*sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  source_term_chgt_perper = 8.0*pi**2*cos(2*sll_pi*x)*cos(2*sll_pi*y) 
  
end function source_term_chgt_perper




real(8) function sol_exacte_chgt_perper(eta1,eta2,params)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  real(8), dimension(:), intent(in), optional :: params
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perper = cos(2*sll_pi*x)*cos(2*sll_pi*y)
  
  
end function sol_exacte_chgt_perper




!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and periodic-dirichlet conditons
!   and dircihlet-diichlet conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------




real(8) function source_term_chgt_perdir(eta1,eta2,params) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  intrinsic :: cos
  intrinsic :: sin
  real(8), dimension(:), intent(in), optional :: params
    
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  source_term_chgt_perdir=&
       (16.0*sll_pi**2*y**4 - 16.0*sll_pi**2*y**2 &
       - 12.0*y**2 + 2.0)*sin(2*sll_pi*x)*cos(2*sll_pi*x)
  
  
end function source_term_chgt_perdir



real(8) function sol_exacte_chgt_perdir(eta1,eta2,params)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  real(8), dimension(:), intent(in), optional :: params
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perdir = &
       y** 2 * (y**2-1)* cos(2*sll_pi*x)*sin(2*sll_pi*x)
  
  
end function sol_exacte_chgt_perdir




!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and dirchlet-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------



real(8) function source_term_chgt_dirper(eta1,eta2) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  intrinsic :: cos
  intrinsic :: sin
  ! -------------------------------------------------
  ! In the case without change of coordinates
  ! -------------------------------------------------
  x =   eta1 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  source_term_chgt_dirper = &
        (16.0*sll_pi**2*x**4 - 16.0*sll_pi**2*x**2 &
       - 12.0*x**2 + 2.0)*sin(2*sll_pi*y)*cos(2*sll_pi*y)
  
end function source_term_chgt_dirper




real(8) function sol_exacte_chgt_dirper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin

  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirper = &
       x ** 2 * (x**2-1)* cos(2* sll_pi*y)*sin(2* sll_pi*y)
  
end function sol_exacte_chgt_dirper

