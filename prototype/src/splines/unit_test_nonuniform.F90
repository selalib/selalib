program nonuniform_spline_tester
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
!contact: mehrenbe@math.unistra.fr for this unit_test program

  use cubic_nonuniform_splines
  use numeric_constants
  implicit none



  sll_int32 :: err
  sll_int32 :: N,i,n_array
  sll_real64,dimension(:), pointer :: node_positions,f,a_in,a_out
  type(cubic_nonunif_spline_1D), pointer :: spl_per  
  sll_real64 :: x,val
  N=32
  n_array = 10
  
  !definition of the mesh
  SLL_ALLOCATE(node_positions(N+1), err)
  do i=1,N+1
    node_positions(i)=real(i-1)/real(N)!sin(2.0_f64*3.1415_f64*real(i-1)/real(N))
  enddo
  
  !definition of the function
  SLL_ALLOCATE(f(N+1), err)  
  f = 1.0_f64
  f(1) = 1.0_f64
  do i=1,N+1
    f(i)=sin(2.0_f64*sll_pi*real(i-1)/real(N))
  enddo
  

  print *, '#proceed to allocate the spline...'  
  spl_per =>  new_cubic_nonunif_spline_1D( N, PERIODIC_SPLINE)
  call compute_spline_nonunif( f, spl_per, node_positions)
  !call compute_spline_nonunif( f, spl_per)
  !print *,spl_per%coeffs
  
  x = spl_per%xmin!1.0_f64!+1.0e-12_f64
  val = interpolate_value_nonunif( x, spl_per )
  print *,val,x
  
  !print *,node_positions
  !stop
  !allocation for a_in and a_out
  SLL_ALLOCATE(a_in(n_array), err)  
  SLL_ALLOCATE(a_out(n_array), err)  
  
  do i=1, n_array
    a_in(i) = spl_per%xmin+real(i-1)/real(n_array-1)*(spl_per%xmax-spl_per%xmin) 
  enddo
  a_in(1)=0.001_f64
  a_in(2)=1.0_f64
  a_in(3)=0.9_f64
  a_in(4)=0.01_f64
  a_in(5)=0.0_f64
  a_in(6)=0.1_f64
  a_in(7)=0.05_f64
  a_in(8)=0.025_f64
  a_in(9)=0.3_f64
  a_in(10)=0.2_f64
  call interpolate_array_value_nonunif( a_in, a_out,n_array, spl_per)
  do i=1, n_array
    print *, a_in(i),a_out(i)
  enddo
  print *, '#proceed to delete the spline...'  
  call delete_cubic_nonunif_spline_1D( spl_per, err)
  
  

end program nonuniform_spline_tester
