program unit_test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_arbitrary_degree_spline_interpolator_2d, only: &
    sll_t_arbitrary_degree_spline_interpolator_2d, &
    sll_o_delete

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_hermite, &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define NPTS1 64
#define NPTS2 64 
#define SPL_DEG 3
#define X1MIN (-2.0_f64*sll_p_pi)
#define X1MAX (+2.0_f64*sll_p_pi)
#define X2MIN (-2.0_f64*sll_p_pi)
#define X2MAX (+2.0_f64*sll_p_pi)
#define TOLERANCE_NODE 1.0E-7_f64
#define TOLERANCE_DER  3.0e-5_f64

sll_real64, dimension(NPTS1,NPTS2) :: f
sll_real64, dimension(NPTS1,NPTS2) :: df_eta1
sll_real64, dimension(NPTS1,NPTS2) :: df_eta2

sll_real64, dimension(NPTS1) :: eta1
sll_real64, dimension(NPTS2) :: eta2
sll_real64, dimension(NPTS2) :: slope_left
sll_real64, dimension(NPTS2) :: slope_right
sll_real64, dimension(NPTS1) :: slope_top
sll_real64, dimension(NPTS1) :: slope_bottom
sll_real64, dimension(NPTS2) :: value_left
sll_real64, dimension(NPTS2) :: value_right
sll_real64, dimension(NPTS1) :: value_top
sll_real64, dimension(NPTS1) :: value_bottom

sll_int32  :: i, j, k=0
sll_real64 :: h1, h2
sll_real64 :: acc(9)
sll_real64 :: acc_der1(9)
sll_real64 :: acc_der2(9)
sll_real64 :: normL2(9)
sll_real64 :: normH1(9)

h1 = (X1MAX-X1MIN)/(NPTS1-1)
h2 = (X2MAX-X2MIN)/(NPTS2-1)
print *, 'h1 = ', h1
print *, 'h2 = ', h2
  
do i=1,NPTS1
  eta1(i) = X1MIN + (i-1)*h1
end do
do j=1,NPTS2
  eta2(j) = X2MIN + (j-1)*h2
end do

print *, '***********************************************************'
print *, '              periodic-periodic case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       =  cos(eta1(i))*cos(eta2(j))
    df_eta1(i,j) = -sin(eta1(i))*cos(eta2(j))
    df_eta2(i,j) = -sin(eta2(j))*cos(eta1(i))
  end do
end do

call check_interpolation( sll_p_periodic, &
                          sll_p_periodic, &
                          sll_p_periodic, &
                          sll_p_periodic  )

print *, '***********************************************************'
print *, '              periodic-dirichlet case'
print *, '***********************************************************'

do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       =  cos(eta1(i))*sin(eta2(j))
    df_eta1(i,j) = -sin(eta1(i))*sin(eta2(j))
    df_eta2(i,j) =  cos(eta1(i))*cos(eta2(j))
  end do
end do

call check_interpolation( sll_p_periodic,  &
                          sll_p_periodic,  &
                          sll_p_dirichlet, &
                          sll_p_dirichlet  )

print *, '***********************************************************'
print *, '              dirichlet-periodic case'
print *, '***********************************************************'

do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       =  sin(eta1(i))*cos(eta2(j))
    df_eta1(i,j) =  cos(eta1(i))*cos(eta2(j))
    df_eta2(i,j) = -sin(eta1(i))*sin(eta2(j))
  end do
end do

call check_interpolation( sll_p_dirichlet,  &
                          sll_p_dirichlet,  &
                          sll_p_periodic, &
                          sll_p_periodic  )

print *, '***********************************************************'
print *, '              dirichlet-dirichlet case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(eta1(i))*sin(eta2(j)) 
    df_eta1(i,j) = cos(eta1(i))*sin(eta2(j))
    df_eta2(i,j) = sin(eta1(i))*cos(eta2(j))
  end do
end do

call check_interpolation( sll_p_dirichlet, &
                          sll_p_dirichlet, &
                          sll_p_dirichlet, &
                          sll_p_dirichlet  )

print *, '***********************************************************'
print *, '              Hermite-dirichlet-dirichlet-hermite case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(eta1(i))*sin(eta2(j)) 
    df_eta1(i,j) = cos(eta1(i))*sin(eta2(j))
    df_eta2(i,j) = sin(eta1(i))*cos(eta2(j))
  end do
end do

call check_interpolation( sll_p_hermite,   &
                          sll_p_dirichlet, &
                          sll_p_dirichlet, &
                          sll_p_hermite    )

print *, '***********************************************************'
print *, '              Hermite-dirichlet-hermite-dirichlet case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(eta1(i))*sin(eta2(j)) 
    df_eta1(i,j) = cos(eta1(i))*sin(eta2(j))
    df_eta2(i,j) = sin(eta1(i))*cos(eta2(j))
  end do
end do

call check_interpolation( sll_p_hermite,   &
                          sll_p_dirichlet, &
                          sll_p_hermite,   &
                          sll_p_dirichlet  )

print *, '***********************************************************'
print *, '              dirichlet-Hermite-hermite-dirichlet case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(eta1(i))*sin(eta2(j)) 
    df_eta1(i,j) = cos(eta1(i))*sin(eta2(j))
    df_eta2(i,j) = sin(eta1(i))*cos(eta2(j))
  end do
end do

call check_interpolation( sll_p_dirichlet, &
                          sll_p_hermite,   &
                          sll_p_hermite,   &
                          sll_p_dirichlet  )

print *, '***********************************************************'
print *, '              dirichlet-Hermite-dirichlet-hermite case'
print *, '***********************************************************'
 
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(eta1(i))*sin(eta2(j)) 
    df_eta1(i,j) = cos(eta1(i))*sin(eta2(j))
    df_eta2(i,j) = sin(eta1(i))*cos(eta2(j))
  end do
end do

call check_interpolation( sll_p_dirichlet, &
                          sll_p_hermite,   &
                          sll_p_dirichlet, &
                          sll_p_hermite  )

print *, '***********************************************************'
print *, '              Hermite-Hermite-Hermite-hermite case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       =   cos(eta1(i))*cos(eta2(j)) 
    df_eta1(i,j) = - sin(eta1(i))*cos(eta2(j))
    df_eta2(i,j) = - cos(eta1(i))*sin(eta2(j))
  end do
end do

call check_interpolation( sll_p_hermite, &
                          sll_p_hermite, &
                          sll_p_hermite, &
                          sll_p_hermite  )

print*, '--------------------------------------------'
print*, ' Average error in nodes'
print*, '--------------------------------------------'
print*, 'Average error (periodic-periodic)                   = ', acc(1)
print*, 'Average error (periodic-dirichlet)                  = ', acc(2)
print*, 'Average error (dirichlet-periodic)                  = ', acc(3)
print*, 'Average error (dirichlet-dirichlet)                 = ', acc(4)
print*, 'Average error (hermite-dirichlet-dirichlet-hermite) = ', acc(5)
print*, 'Average error (hermite-dirichlet-hermite-dirichlet) = ', acc(6)
print*, 'Average error (dirichlet-hermite-hermite-dirichlet) = ', acc(7)
print*, 'Average error (dirichlet-hermite-dirichlet-hermite) = ', acc(8)
print*, 'Average error (hermite-hermite-hermite-hermite)     = ', acc(9)

print*, '--------------------------------------------'
print*, ' Average error in nodes first derivative eta1'
print*, '--------------------------------------------'
print*, 'Average error (periodic-periodic)                   = ',acc_der1(1)
print*, 'Average error (periodic-dirichlet)                  = ',acc_der1(2)
print*, 'Average error (dirichlet-periodic)                  = ',acc_der1(3)
print*, 'Average error (dirichlet-dirichlet)                 = ',acc_der1(4)
print*, 'Average error (hermite-dirichlet-dirichlet-hermite) = ',acc_der1(5)
print*, 'Average error (hermite-dirichlet-hermite-dirichlet) = ',acc_der1(6)
print*, 'Average error (dirichlet-hermite-hermite-dirichlet) = ',acc_der1(7)
print*, 'Average error (dirichlet-hermite-dirichlet-hermite) = ',acc_der1(8)
print*, 'Average error (hermite-hermite-hermite-hermite)     = ',acc_der1(9)

print*, '--------------------------------------------'
print*, ' Average error in nodes first derivative eta2'
print*, '--------------------------------------------'
print*, 'Average error (periodic-periodic)                   = ',acc_der2(1)
print*, 'Average error (periodic-dirichlet)                  = ',acc_der2(2)
print*, 'Average error (dirichlet-periodic)                  = ',acc_der2(3)
print*, 'Average error (dirichlet-dirichlet)                 = ',acc_der2(4)
print*, 'Average error (hermite-dirichlet-dirichlet-hermite) = ',acc_der2(5)
print*, 'Average error (hermite-dirichlet-hermite-dirichlet) = ',acc_der2(6)
print*, 'Average error (dirichlet-hermite-hermite-dirichlet) = ',acc_der2(7)
print*, 'Average error (dirichlet-hermite-dirichlet-hermite) = ',acc_der2(8)
print*, 'Average error (hermite-hermite-hermite-hermite)     = ',acc_der2(9)

print*, '--------------------------------------------'
print*, ' Error norm L2'
print*, '--------------------------------------------'
print*, 'Error norm L2 (periodic-periodic)                   = ',sqrt(normL2(1)),h1**(SPL_DEG)
print*, 'Error norm L2 (periodic-dirichlet)                  = ',sqrt(normL2(2)),h1**(SPL_DEG)
print*, 'Error norm L2 (dirichlet-periodic)                  = ',sqrt(normL2(3)),h1**(SPL_DEG)
print*, 'Error norm L2 (dirichlet-dirichlet)                 = ',sqrt(normL2(4)),h1**(SPL_DEG)
print*, 'Error norm L2 (hermite-dirichlet-dirichlet-hermite) = ',sqrt(normL2(5)),h1**(SPL_DEG)
print*, 'Error norm L2 (hermite-dirichlet-hermite-dirichlet) = ',sqrt(normL2(6)),h1**(SPL_DEG)
print*, 'Error norm L2 (dirichlet-hermite-hermite-dirichlet) = ',sqrt(normL2(7)),h1**(SPL_DEG)
print*, 'Error norm L2 (dirichlet-hermite-dirichlet-hermite) = ',sqrt(normL2(8)),h1**(SPL_DEG)
print*, 'Error norm L2 (hermite-hermite-hermite-hermite)     = ',sqrt(normL2(9)),h1**(SPL_DEG)
print*, '--------------------------------------------'
print*, ' Error norm H1'
print*, '--------------------------------------------'
print*, 'Error norm H1 (periodic-periodic)                   = ',sqrt(normH1(1)),h1**(SPL_DEG-3)
print*, 'Error norm H1 (periodic-dirichlet)                  = ',sqrt(normH1(2)),h1**(SPL_DEG-3)
print*, 'Error norm H1 (dirichlet-periodic)                  = ',sqrt(normH1(3)),h1**(SPL_DEG-3)
print*, 'Error norm H1 (dirichlet-dirichlet)                 = ',sqrt(normH1(4)),h1**(SPL_DEG-3)
print*, 'Error norm H1 (hermite-dirichlet-dirichlet-hermite) = ',sqrt(normH1(5)),h1**(SPL_DEG-3)
print*, 'Error norm H1 (hermite-dirichlet-hermite-dirichlet) = ',sqrt(normH1(6)),h1**(SPL_DEG-3)
print*, 'Error norm H1 (dirichlet-hermite-hermite-dirichlet) = ',sqrt(normH1(7)),h1**(SPL_DEG-3)
print*, 'Error norm H1 (dirichlet-hermite-dirichlet-hermite) = ',sqrt(normH1(8)),h1**(SPL_DEG-3)
print*, 'Error norm H1 (hermite-hermite-hermite-hermite)     = ',sqrt(normH1(9)),h1**(SPL_DEG-3)


if (( sqrt(normL2(1)) <= h1*h1*h1     )   .and. &
    ( sqrt(normL2(2)) <= h1*h1*h1     )   .and. &
    ( sqrt(normL2(3)) <= h1*h1*h1     )   .and. &
    ( sqrt(normL2(4)) <= h1*h1*h1     )   .and. &
    ( sqrt(normL2(5)) <= h1*h1*h1     )   .and. &
    ( sqrt(normL2(6)) <= h1*h1*h1     )   .and. &
    ( sqrt(normL2(7)) <= h1*h1*h1     )   .and. &
    ( sqrt(normL2(8)) <= h1*h1*h1     )   .and. &
    ( sqrt(normL2(9)) <= h1*h1*h1     )   .and. &
    ( sqrt(normH1(1)) <= h1**(SPL_DEG-3)) .and. &
    ( sqrt(normH1(2)) <= h1**(SPL_DEG-3)) .and. &
    ( sqrt(normH1(3)) <= h1**(SPL_DEG-3)) .and. &
    ( sqrt(normH1(4)) <= h1**(SPL_DEG-3)) .and. &
    ( sqrt(normH1(5)) <= h1**(SPL_DEG-3)) .and. &
    ( sqrt(normH1(6)) <= h1**(SPL_DEG-3)) .and. &
    ( sqrt(normH1(7)) <= h1**(SPL_DEG-3)) .and. &
    ( sqrt(normH1(8)) <= h1**(SPL_DEG-3)) .and. &
    ( sqrt(normH1(9)) <= h1**(SPL_DEG-3))) then
    
  print *, 'PASSED'

else

  print *, 'FAILED'

end if

contains

subroutine check_interpolation( bc_eta1_min, &
                                bc_eta1_max, &
                                bc_eta2_min, &
                                bc_eta2_max  )

sll_int32,  intent(in)  :: bc_eta1_min
sll_int32,  intent(in)  :: bc_eta1_max
sll_int32,  intent(in)  :: bc_eta2_min
sll_int32,  intent(in)  :: bc_eta2_max

type(sll_t_arbitrary_degree_spline_interpolator_2d) :: ad2d

sll_real64 :: deriv1_val 
sll_real64 :: deriv2_val
sll_real64 :: node_val 

k = k+1

call ad2d%initialize( NPTS1,         &
                      NPTS2,         &
                      X1MIN,         &
                      X1MAX,         &
                      X2MIN,         &
                      X2MAX,         &
                      bc_eta1_min,   &
                      bc_eta1_max,   &
                      bc_eta2_min,   &
                      bc_eta2_max,   &
                      SPL_DEG,       &
                      SPL_DEG )
  
slope_left   = df_eta1(1,:)
slope_right  = df_eta1(NPTS1,:)
slope_top    = df_eta2(:,1)
slope_bottom = df_eta2(:,NPTS2)
value_left   = f(1,:)
value_right  = f(NPTS1,:)
value_top    = f(:,NPTS2)
value_bottom = f(:,1)

call ad2d%set_values_at_boundary(value_left,value_right,value_bottom,value_top)
call ad2d%set_slopes_at_boundary(slope_left,slope_right,slope_bottom,slope_top)
call ad2d%compute_interpolants( f, eta1, NPTS1, eta2, NPTS2)
  
acc(k)      = 0.0_f64
acc_der1(k) = 0.0_f64
acc_der2(k) = 0.0_f64
normL2(k)   = 0.0_F64
normH1(k)   = 0.0_f64
  
do j=1,NPTS2
  do i=1,NPTS1-1
    node_val        = ad2d%interpolate_from_interpolant_value(eta1(i),eta2(j))
    normL2(k)       = normL2(k) + (node_val-f(i,j))**2 *h1*h2
    acc(k)          = acc(k) + abs(node_val-f(i,j))
    deriv1_val      = ad2d%interpolate_from_interpolant_derivative_eta1(eta1(i),eta2(j))   
    acc_der1(k)     = acc_der1(k) + abs(deriv1_val-df_eta1(i,j))
    normH1(k)       = normH1(k) + (deriv1_val-df_eta1(i,j))**2 *h1*h2
    deriv2_val      = ad2d%interpolate_from_interpolant_derivative_eta2(eta1(i),eta2(j))
    acc_der2(k)     = acc_der2(k) + abs(deriv2_val-df_eta2(i,j))
    normH1(k)       = normH1(k) + (deriv2_val-df_eta2(i,j))**2 *h1*h2
  end do
end do

acc(k)      = acc(k)/(NPTS1*NPTS2)
acc_der1(k) = acc_der1(k)/(NPTS1*NPTS2)
acc_der2(k) = acc_der2(k)/(NPTS1*NPTS2)

call sll_o_delete(ad2d)

end subroutine check_interpolation

end program unit_test

 
