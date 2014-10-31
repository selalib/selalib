program unit_test
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_interpolators.h"
#include "sll_file_io.h"

implicit none

#define NPTS1 64
#define NPTS2 64 
#define SPL_DEG 3
#define X1MIN -2.0_f64*sll_pi
#define X1MAX +2.0_f64*sll_pi
#define X2MIN -2.0_f64*sll_pi
#define X2MAX +2.0_f64*sll_pi
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

sll_int32  :: ierr
sll_int32  :: i, j, k=0
sll_real64 :: h1, h2
sll_real64 :: acc(9)
sll_real64 :: acc_der1(9)
sll_real64 :: acc_der2(9)
sll_real64 :: normL2(9)
sll_real64 :: normH1(9)
logical :: result

result = .true.

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
    df_eta2(i,j) = -sin(eta2(j))*cos(eta1(j))
  end do
end do

call check_interpolation( SLL_PERIODIC, &
                          SLL_PERIODIC, &
                          SLL_PERIODIC, &
                          SLL_PERIODIC  )

print *, '***********************************************************'
print *, '              periodic-dirichlet case'
print *, '***********************************************************'

do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       =  sin(eta1(i))*cos(eta2(j))
    df_eta1(i,j) =  cos(eta1(i))*cos(eta2(j))
    df_eta2(i,j) = -sin(eta1(i))*sin(eta2(j))
  end do
end do

call check_interpolation( SLL_PERIODIC,  &
                          SLL_PERIODIC,  &
                          SLL_DIRICHLET, &
                          SLL_DIRICHLET  )

print *, '***********************************************************'
print *, '              dirichlet-periodic case'
print *, '***********************************************************'

do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       =  sin(2*sll_pi*eta1(i))*cos(2*sll_pi*eta2(j))
    df_eta1(i,j) =  2*sll_pi*cos(2*sll_pi*eta1(i))*cos(2*sll_pi*eta2(j))
    df_eta2(i,j) = -2*sll_pi*sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))
  end do
end do

call check_interpolation( SLL_PERIODIC,  &
                          SLL_PERIODIC,  &
                          SLL_DIRICHLET, &
                          SLL_DIRICHLET  )

print *, '***********************************************************'
print *, '              dirichlet-dirichlet case'
print *, '***********************************************************'

print *, '***********************************************************'
print *, '              dirichlet-dirichlet case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j)) 
    df_eta1(i,j) = 2*sll_pi*cos(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))
    df_eta2(i,j) = 2*sll_pi*sin(2*sll_pi*eta1(i))*cos(2*sll_pi*eta2(j))
  end do
end do

call check_interpolation( SLL_DIRICHLET, &
                          SLL_DIRICHLET, &
                          SLL_DIRICHLET, &
                          SLL_DIRICHLET  )

print *, '***********************************************************'
print *, '              Hermite-dirichlet-dirichlet-hermite case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j)) 
    df_eta1(i,j) = 2*sll_pi*cos(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))
    df_eta2(i,j) = 2*sll_pi*sin(2*sll_pi*eta1(i))*cos(2*sll_pi*eta2(j))
  end do
end do

call check_interpolation( SLL_HERMITE,   &
                          SLL_DIRICHLET, &
                          SLL_DIRICHLET, &
                          SLL_HERMITE    )

print *, '***********************************************************'
print *, '              Hermite-dirichlet-hermite-dirichlet case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j)) 
    df_eta1(i,j) = 2*sll_pi*cos(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))
    df_eta2(i,j) = 2*sll_pi*sin(2*sll_pi*eta1(i))*cos(2*sll_pi*eta2(j))
  end do
end do

call check_interpolation( SLL_HERMITE,   &
                          SLL_DIRICHLET, &
                          SLL_HERMITE,   &
                          SLL_DIRICHLET  )

print *, '***********************************************************'
print *, '              dirichlet-Hermite-hermite-dirichlet case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j)) 
    df_eta1(i,j) = 2*sll_pi*cos(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))
    df_eta2(i,j) = 2*sll_pi*sin(2*sll_pi*eta1(i))*cos(2*sll_pi*eta2(j))
  end do
end do

call check_interpolation( SLL_DIRICHLET, &
                          SLL_HERMITE,   &
                          SLL_HERMITE,   &
                          SLL_DIRICHLET  )

print *, '***********************************************************'
print *, '              dirichlet-Hermite-dirichlet-hermite case'
print *, '***********************************************************'
 
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j)) 
    df_eta1(i,j) = 2*sll_pi*cos(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))
    df_eta2(i,j) = 2*sll_pi*sin(2*sll_pi*eta1(i))*cos(2*sll_pi*eta2(j))
  end do
end do

call check_interpolation( SLL_DIRICHLET, &
                          SLL_HERMITE,   &
                          SLL_DIRICHLET, &
                          SLL_HERMITE  )

print *, '***********************************************************'
print *, '              Hermite-Hermite-Hermite-hermite case'
print *, '***********************************************************'
  
do j=1,NPTS2
  do i=1,NPTS1
    f(i,j)       = sin(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j)) 
    df_eta1(i,j) = 2*sll_pi*cos(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))
    df_eta2(i,j) = 2*sll_pi*sin(2*sll_pi*eta1(i))*cos(2*sll_pi*eta2(j))
  end do
end do

call check_interpolation( SLL_HERMITE, &
                          SLL_HERMITE, &
                          SLL_HERMITE, &
                          SLL_HERMITE  )

print*, '--------------------------------------------'
print*, ' Average error in nodes'
print*, '--------------------------------------------'
print *, 'Average error in nodes (periodic-periodic) = ', acc/(NPTS1*NPTS2)
call test_value_for_acceptable_error(acc/(NPTS1*NPTS2), TOLERANCE_NODE,result)
print *, 'Average error in nodes (dirichlet-dirichlet) = ', acc3/(NPTS1*NPTS2)
call test_value_for_acceptable_error(acc3/(NPTS1*NPTS2), TOLERANCE_NODE,result)
print *, 'Average error in nodes (dirichlet-periodic) = ', acc2/(NPTS1*NPTS2)
call test_value_for_acceptable_error(acc2/(NPTS1*NPTS2), TOLERANCE_NODE,result)
print *, 'Average error in nodes (periodic-dirichlet) = ', acc1/(NPTS1*NPTS2)
call test_value_for_acceptable_error(acc1/(NPTS1*NPTS2), TOLERANCE_NODE,result)
print *, 'Average error in nodes (hermite-dirichlet-dirichlet-hermite) = ', acc4/(NPTS1*NPTS2)
call test_value_for_acceptable_error(acc4/(NPTS1*NPTS2), TOLERANCE_NODE,result)
print *, 'Average error in nodes (hermite-dirichlet-hermite-dirichlet) = ', acc5/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc5/(NPTS1*NPTS2), TOLERANCE_NODE,result)
print *, 'Average error in nodes (dirichlet-hermite-hermite-dirichlet) = ', acc6/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc6/(NPTS1*NPTS2), TOLERANCE_NODE,result)
print *, 'Average error in nodes (dirichlet-hermite-dirichlet-hermite) = ', acc7/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc7/(NPTS1*NPTS2), TOLERANCE_NODE,result)
print *, 'Average error in nodes (hermite-hermite-hermite-hermite) = ', acc8/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc8/(NPTS1*NPTS2), TOLERANCE_NODE,result)

print*, '--------------------------------------------'
print*, ' Average error in nodes first derivative eta1'
print*, '--------------------------------------------'
print *,'Average error first derivative eta1(dirichlet-dirichlet)=',acc3_der1/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc3_der1/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta1(dirichlet-periodic)=',acc2_der1/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc2_der1/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta1(periodic-dirichlet)=',acc1_der1/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc1_der1/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta1(periodic-periodic)=',acc_der1/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc_der1/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta1(hermite-dirichlet-dirichlet-hermite)=',acc4_der1/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc4_der1/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta1(hermite-dirichlet-hermite-dirichlet)=',acc5_der1/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc5_der1/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta1(dirichlet-hermite-hermite-dirichlet)=',acc6_der1/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc6_der1/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta1(dirichlet-hermite-dirichlet-hermite)=',acc7_der1/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc7_der1/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta1(hermite-hermite-hermite-hermite)=',acc8_der1/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc8_der1/(NPTS1*NPTS2),TOLERANCE_DER,result)

print*, '--------------------------------------------'
print*, ' Average error in nodes first derivative eta2'
print*, '--------------------------------------------'
print *,'Average error first derivative eta2(dirichlet-dirichlet)=',acc3_der2/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc3_der2/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta2(dirichlet-periodic)=',acc2_der2/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc2_der2/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta2(periodic-dirichlet)=',acc1_der2/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc1_der2/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta2(periodic-periodic)=',acc_der2/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc_der2/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta2(hermite-dirichlet-dirichlet-hermite)=',acc4_der2/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc4_der2/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta2(hermite-dirichlet-hermite-dirichlet)=',acc5_der2/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc5_der2/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta2(dirichlet-hermite-hermite-dirichlet)=',acc6_der2/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc6_der2/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error in nodes first derivative eta2(dirichlet-hermite-dirichlet-hermite)=',acc7_der2/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc7_der2/(NPTS1*NPTS2),TOLERANCE_DER,result)
print *,'Average error first derivative eta2(hermite-hermite-hermite-hermite)=',acc8_der2/(NPTS1*NPTS2)
call test_value_for_acceptable_error( acc8_der2/(NPTS1*NPTS2),TOLERANCE_DER,result)


print*, '--------------------------------------------'
print*, ' Error norm L2'
print*, '--------------------------------------------'
print *,'Error norm L2 (dirichlet-dirichlet)=',                sqrt(normL2_3),h1**(SPL_DEG)*(2*sll_pi)
print *,'Error norm L2 (dirichlet-periodic)=',                 sqrt(normL2_2),h1**(SPL_DEG)*(2*sll_pi)
print *,'Error norm L2 (periodic-dirichlet)=',                 sqrt(normL2_1),h1**(SPL_DEG)*(2*sll_pi)
print *,'Error norm L2 (periodic-periodic)=',                  sqrt(normL2_0),h1**(SPL_DEG)*(2*sll_pi)
print *,'Error norm L2 (hermite-dirichlet-dirichlet-hermite)=',sqrt(normL2_4),h1**(SPL_DEG)*(2*sll_pi)
print *,'Error norm L2 (hermite-dirichlet-hermite-dirichlet)=',sqrt(normL2_5),h1**(SPL_DEG)*(2*sll_pi)
print *,'Error norm L2 (dirichlet-hermite-hermite-dirichlet)=',sqrt(normL2_6),h1**(SPL_DEG)*(2*sll_pi)
print *,'Error norm L2 (dirichlet-hermite-dirichlet-hermite)=',sqrt(normL2_7),h1**(SPL_DEG)*(2*sll_pi)
print *,'Error norm L2 (hermite-hermite-hermite-hermite)=',    sqrt(normL2_8),h1**(SPL_DEG)*(2*sll_pi)
print*, '--------------------------------------------'
print*, ' Error norm H1'
print*, '--------------------------------------------'
print *,'Error norm H1 (dirichlet-dirichlet)=',                sqrt(normH1_3),h1**(SPL_DEG-3)*(2*sll_pi)**2
print *,'Error norm H1 (dirichlet-periodic)=',                 sqrt(normH1_2),h1**(SPL_DEG-3)*(2*sll_pi)**2
print *,'Error norm H1 (periodic-dirichlet)=',                 sqrt(normH1_1),h1**(SPL_DEG-3)*(2*sll_pi)**2
print *,'Error norm H1 (periodic-periodic)=',                  sqrt(normH1_0),h1**(SPL_DEG-3)*(2*sll_pi)**2
print *,'Error norm H1 (hermite-dirichlet-dirichlet-hermite)=',sqrt(normH1_4),h1**(SPL_DEG-3)*(2*sll_pi)**2
print *,'Error norm H1 (hermite-dirichlet-hermite-dirichlet)=',sqrt(normH1_5),h1**(SPL_DEG-3)*(2*sll_pi)**2
print *,'Error norm H1 (dirichlet-hermite-hermite-dirichlet)=',sqrt(normH1_6),h1**(SPL_DEG-3)*(2*sll_pi)**2
print *,'Error norm H1 (dirichlet-hermite-dirichlet-hermite)=',sqrt(normH1_7),h1**(SPL_DEG-3)*(2*sll_pi)**2
print *,'Error norm H1 (hermite-hermite-hermite-hermite)=',    sqrt(normH1_8),h1**(SPL_DEG-3)*(2*sll_pi)**2

if (  ( sqrt(normL2_0) <= h1**(SPL_DEG)*(2*sll_pi))   .AND. &
      ( sqrt(normL2_1) <= h1**(SPL_DEG)*(2*sll_pi))   .AND. &
      ( sqrt(normL2_2) <= h1**(SPL_DEG)*(2*sll_pi))   .AND. &
      ( sqrt(normL2_3) <= h1**(SPL_DEG)*(2*sll_pi))   .AND. &
      ( sqrt(normL2_4) <= h1**(SPL_DEG)*(2*sll_pi))   .AND. &
      ( sqrt(normL2_5) <= h1**(SPL_DEG)*(2*sll_pi))   .AND. &
      ( sqrt(normL2_6) <= h1**(SPL_DEG)*(2*sll_pi))   .AND. &
      ( sqrt(normL2_7) <= h1**(SPL_DEG)*(2*sll_pi))   .AND. &
      ( sqrt(normH1_0) <= h1**(SPL_DEG-3)*(2*sll_pi)**2) .AND. &
      ( sqrt(normH1_1) <= h1**(SPL_DEG-3)*(2*sll_pi)**2) .AND. &
      ( sqrt(normH1_2) <= h1**(SPL_DEG-3)*(2*sll_pi)**2) .AND. &
      ( sqrt(normH1_4) <= h1**(SPL_DEG-3)*(2*sll_pi)**2) .AND. &
      ( sqrt(normH1_5) <= h1**(SPL_DEG-3)*(2*sll_pi)**2) .AND. &
      ( sqrt(normH1_6) <= h1**(SPL_DEG-3)*(2*sll_pi)**2) .AND. &
      ( sqrt(normH1_7) <= h1**(SPL_DEG-3)*(2*sll_pi)**2) .AND. &
      ( sqrt(normH1_8) <= h1**(SPL_DEG-3)*(2*sll_pi)**2) .AND. &
      ( sqrt(normH1_3) <= h1**(SPL_DEG-3)*(2*sll_pi)**2)) then
    
     print *, 'PASSED'
end if

contains

subroutine test_value_for_acceptable_error( value, max_error, boolean )
sll_real64, intent(in) :: value
sll_real64, intent(in) :: max_error
logical, intent(inout) :: boolean

if( value <= max_error ) then
   boolean = boolean .and. .true.
else
   boolean = boolean .and. .false.
end if

end subroutine test_value_for_acceptable_error

subroutine check_interpolation( bc_eta1_min, &
                                bc_eta1_max, &
                                bc_eta2_min, &
                                bc_eta2_max  )

sll_int32,  intent(in)  :: bc_eta1_min
sll_int32,  intent(in)  :: bc_eta1_max
sll_int32,  intent(in)  :: bc_eta2_min
sll_int32,  intent(in)  :: bc_eta2_max

type(sll_arbitrary_degree_spline_interpolator_2d) :: ad2d

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
    node_val        = ad2d%interpolate_value(eta1(i),eta2(j))
    normL2(k)       = normL2(k) + (node_val-f(i,j))**2 *h1*h2
    acc(k)          = acc(k) + abs(node_val-f(i,j))
    deriv1_val      = ad2d%interpolate_derivative_eta1(eta1(i),eta2(j))   
    acc_der1(k)     = acc_der1(k) + abs(deriv1_val-df_eta1(i,j))
    normH1(k)       = normH1(k) + (deriv1_val-df_eta1(i,j))**2 *h1*h2
    deriv2_val      = ad2d%interpolate_derivative_eta2(eta1(i),eta2(j))
    acc_der2(k)     = acc_der2(k) + abs(deriv2_val-df_eta2(i,j))
    normH1(k)       = normH1(k) + (deriv2_val-df_eta2(i,j))**2 *h1*h2
  end do
end do

call sll_delete(ad2d)

end subroutine check_interpolation

end program unit_test

 
