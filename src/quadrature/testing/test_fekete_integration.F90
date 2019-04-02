program test_fekete_integration

#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_hexagonal_meshes, only: &
    sll_t_hex_mesh_2d

  use sll_m_fekete_integration, only: &
       sll_f_fekete_integral, &
       sll_s_write_all_django_files, &
       sll_f_fekete_points_and_weights, &
       sll_s_fekete_order_num

  use test_function_module, only: &
      one_2D

  implicit none

  sll_int32 :: i
  sll_int32 :: j
  sll_int32 :: n
  sll_int32 :: ierr
  sll_int32 :: degree

  character(len=18) :: string

  ! For the fekete quadrature:
  sll_real64, dimension(2, 3) :: pxy1
  sll_real64, dimension(2, 3) :: pxy2
  sll_real64, dimension(:,:), allocatable :: xyw
  sll_real64 :: app_res
  sll_int32  :: rule
  sll_int32  :: num_cells
  type(sll_t_hex_mesh_2d), pointer :: mesh
  sll_real64, dimension(:,:), allocatable     :: knots
  sll_int32,  dimension(:,:), allocatable     :: LM
  character(len=10) :: transf

  write(*,"(/,a)") "*********************************** "
  write(*,"(a)") "       FEKETE QUAD TEST       "
  write(*,"(a)") "*********************************** "

  !Definition of first triangle
  pxy1(:,1) = (/ 0._f64, 0._f64 /)
  pxy1(:,2) = (/ 1._f64, 0._f64 /)
  pxy1(:,3) = (/ 0._f64, 1._f64 /)

  !Definition of first triangle
  pxy2(:,1) = (/ 1._f64, 0._f64 /)
  pxy2(:,2) = (/ 1._f64, 1._f64 /)
  pxy2(:,3) = (/ 0._f64, 1._f64 /)

  rule = 2
  call sll_s_fekete_order_num ( rule, n )
  SLL_ALLOCATE(xyw(1:3, 1:n), ierr)

  write(*,"(a)") " Computing Fekete points and weights on reference triangle "
  write(*,"(/,a)") "           x                   y                    w"
  xyw = sll_f_fekete_points_and_weights(pxy1, rule)

  write (string, '( "(",I2,"f20.15)" )' )  n
  do j = 1, n
     write(*, string) (xyw(i,j), i = 1, 3)
  end do

  print *, "sum weights = ", SUM(xyw(3,:))
  write(*,"(/,a)") " --Test for a constant real function (=1) "
  write(*,"(a)") "    on the squared domain [0,1]^2 divided on 2 triangles "

  app_res = sll_f_fekete_integral(one_2D, pxy1) + sll_f_fekete_integral(one_2D, pxy2)

  write (*,"(a, f20.12, a ,/)") " aprox = ", app_res, " (expected = 1.)"


  ! Writing all django files....................................................
  num_cells = 10
  degree = 1
  rule = 1
  transf="TOKAMAK"
  call sll_s_write_all_django_files(num_cells, degree, rule, trim(transf))
  print *, ""
  print *, "*********** wrote all django files ***********"
  print *, "   - number of cells    : ", num_cells
  print *, "   - degree of splines  : ", degree
  print *, "   - rule of quadrature : ", rule
  print *, "   - coo transformation : ", trim(transf)
  print *, ""

  print*, 'PASSED'

end program test_fekete_integration
