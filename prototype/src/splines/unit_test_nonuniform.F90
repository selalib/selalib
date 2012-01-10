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
  sll_real64,dimension(:), pointer :: node_positions,f_per,f_hrmt,a_in,a_out_per,a_out_hrmt,sol_per,sol_hrmt
  type(cubic_nonunif_spline_1D), pointer :: spl_per, spl_hrmt  
  sll_real64 :: x,val,sl,sr
  N=16
  n_array = 10000
  
  !definition of the mesh
  SLL_ALLOCATE(node_positions(N+1), err)
  do i=1,N+1
    !node_positions(i)=real(i-1)/real(N)+1.0_f64/real(N)*sin(2.0_f64*sll_pi*real(i-1)/real(N))
    node_positions(i)=2.0_f64*real(i-1)/real(N)!sin(2.0_f64*3.1415_f64*real(i-1)/real(N))
  enddo
  node_positions(1)=5.1_f64
  do i=2,N+1
    node_positions(i)=node_positions(i-1)+1.0_f64/real(N)*(1.0_f64+1.0_f64/real(N*N)+sin(2.0_f64*sll_pi*real(i-1)/real(N)))
  enddo
  
  !print *,node_positions
  !stop
  !definition of the function
  SLL_ALLOCATE(f_hrmt(N+1), err)  
  f_hrmt = 0.0_f64
  !f_hrmt(1) = 0.1_f64
  !f_hrmt(N+1) = 0.2_f64
  do i=1,N+1
    x=node_positions(i)
    !f_hrmt(i)=sin(2.0_f64*sll_pi*x/(node_positions(N+1)-node_positions(1)))
    f_hrmt(i)=exp(x)
  enddo
  sl=2.0_f64*sll_pi
  sr=2.0_f64*sll_pi
  sl=exp(node_positions(1))
  sr=exp(node_positions(N+1))
  !sl=0.0_f64
  !sr=0.0_f64
  
  SLL_ALLOCATE(f_per(N+1), err)  
  f_per = 1.0_f64
  f_per(1) = 1.0_f64
  do i=1,N+1
    x=node_positions(i)
    !f_per(i)=sin(2.0_f64*sll_pi*x)
    f_per(i)=sin(2.0_f64*sll_pi*x/(node_positions(N+1)-node_positions(1)))
  enddo

  !f_per = f_hrmt

  print *, '#proceed to allocate the spline...'  
  spl_per =>  new_cubic_nonunif_spline_1D( N, PERIODIC_SPLINE)
  spl_hrmt =>  new_cubic_nonunif_spline_1D( N, HERMITE_SPLINE)
  call compute_spline_nonunif( f_per, spl_per, node_positions)
  call compute_spline_nonunif( f_hrmt, spl_hrmt, node_positions,sl,sr)

  
  !interpolating points
  !allocation for a_in and a_out
  SLL_ALLOCATE(a_in(n_array), err)  
  SLL_ALLOCATE(a_out_per(n_array), err)  
  SLL_ALLOCATE(a_out_hrmt(n_array), err)  
  SLL_ALLOCATE(sol_per(n_array), err)  
  SLL_ALLOCATE(sol_hrmt(n_array), err)  


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
  do i=1,10
    a_in(i) = spl_per%xmin+a_in(i)*(spl_per%xmax-spl_per%xmin)
  enddo

  do i=1,n_array
    x=a_in(i)
    sol_per(i)=sin(2.0_f64*sll_pi*x/(node_positions(N+1)-node_positions(1)))
  enddo

  do i=1,n_array
    x=a_in(i)
    sol_hrmt(i)=exp(x)
  enddo

  !sol_per=sol_hrmt
  !call compute_spline_nonunif( f, spl_per)
  
  x = spl_per%xmin!1.0_f64!+1.0e-12_f64
  x = x+0.1234589_f64*(spl_per%xmax-spl_per%xmin)
  !i=N+1
  !x = node_positions(i)
  val = interpolate_value_nonunif( x, spl_per )
  print *,'#',val,x,sin(2.0_f64*sll_pi*x)-val
  val = interpolate_value_nonunif( x, spl_hrmt )
  print *,'#',val,x,sin(2.0_f64*sll_pi*x)-val
  

  call interpolate_array_value_nonunif( a_in, a_out_per,n_array, spl_per)
  call interpolate_array_value_nonunif( a_in, a_out_hrmt,n_array, spl_hrmt)
  do i=1, n_array
    print *, a_in(i),a_out_per(i),a_out_hrmt(i),sol_per(i),sol_hrmt(i)
  enddo
  print *, '#proceed to delete the spline...'  
  call delete_cubic_nonunif_spline_1D( spl_per, err)
  call delete_cubic_nonunif_spline_1D( spl_hrmt, err)
  
  

end program nonuniform_spline_tester
