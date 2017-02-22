program test_splines_pp_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"


  use sll_m_arbitrary_degree_splines, only: &
       sll_s_uniform_b_splines_at_x

  use sll_m_arbitrary_degree_spline_interpolator_1d, only: &
       sll_s_initialize_ad1d_interpolator, &
       sll_t_arbitrary_degree_spline_interpolator_1d

  use sll_m_boundary_condition_descriptors, only: &
       sll_p_periodic

  use sll_m_constants, only: &
       sll_p_twopi
  
  use sll_m_splines_pp 


  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  sll_real64, parameter :: inv_2  = 1._f64/2._f64
  sll_real64, parameter :: inv_3  = 1._f64/3._f64
  sll_real64, parameter :: inv_4  = 1._f64/4._f64
  sll_real64, parameter :: inv_6  = 1._f64/6._f64
  sll_real64, parameter :: inv_8  = 1._f64/8._f64
  sll_real64, parameter :: inv_10  = 1._f64/10._f64
  sll_real64, parameter :: inv_24 = 1._f64/24._f64
  
  type(sll_t_spline_pp_2d) :: spline_pp
  sll_int32 :: degree(2)
  sll_int32 :: n_cells(2)
  logical   :: fail
 
  n_cells=[50,50]
  degree=[3,3]
  
  call spline_test(spline_pp,degree,n_cells,fail)
     
  
  if(fail .eqv. .false.) then
     write(*,*) 'PASSED'
  else  
     write(*,*)'FAILED'
     stop
  end if

contains

  subroutine spline_test(spline_pp,degree,n_cells,fail)
    type(sll_t_spline_pp_2d), intent(inout) :: spline_pp !arbitrary degree spline
    sll_int32, intent(in) :: degree(2) !spline degree
    sll_int32, intent(in) :: n_cells(2) !grid cells
    logical, intent(out)    :: fail
    
    sll_real64 :: b_coeffs(n_cells(1)*n_cells(2))
    sll_real64 :: pp_coeffs((degree(1)+1)*(degree(2)+1),n_cells(1)*n_cells(2))
    sll_real64 :: xp(2) !particle position
    sll_real64 :: res
    sll_real64 :: res2
    sll_real64 :: domain(2,2) !left and right bound of the 1d grid
    sll_real64 :: delta_x(2) !size of gridcells
    sll_real64, allocatable :: val1(:)
    sll_real64, allocatable :: val2(:)
    sll_int32 :: index2d
    sll_int32 :: index1d(2)
    sll_int32 :: indices(2)
    sll_real64 :: xi(2)
    sll_int32 :: i,j
       
    fail=.false.
    allocate(val1(degree(1)+1))
    allocate(val2(degree(2)+1))
    domain(1,:)=0._f64
    domain(2,:)=sll_p_twopi 
    delta_x=(domain(2,:)-domain(1,:))/real(n_cells,f64)

    call random_seed()
    call random_number(b_coeffs)
    
    call sll_s_spline_pp_init_2d( spline_pp, degree, n_cells)
    call sll_s_spline_pp_b_to_pp_2d(spline_pp, n_cells, b_coeffs, pp_coeffs)

    call random_seed()
    call random_number(xp)
    !xp=[0.8_f64,0.2_f64]
    xp=xp*(domain(2,:)-domain(1,:))
    
    ! calculate index of gridcell
    xi = (xp - domain(1,:))/delta_x
    indices = floor(xi)+1
    xi = xi - real(indices-1, f64)
    res= sll_f_spline_pp_horner_2d(degree, pp_coeffs, xi, indices, n_cells)
    
    indices = indices - degree
    
    call sll_s_uniform_b_splines_at_x(degree(1), xi(1), val1(:))
    call sll_s_uniform_b_splines_at_x(degree(2), xi(2), val2(:))
    
    
    res2 = 0.0_f64
    do i = 1, degree(1)+1
       index1d(1) = modulo(indices(1)+i-2, n_cells(1)) 
       do j = 1, degree(2)+1
          index1d(2) = modulo( indices(2)+j-2, n_cells(2))
          index2d = index1d(1) + index1d(2)*n_cells(1) + 1
          res2 = res2 + b_coeffs(index2d) * val1(i) * val2(j)
       end do
    end do
    !print*, 'res-res2=', res-res2
    !write(*,*) 'Fehler horner vs normal:', abs(res-res2)
    if(abs(res-res2)>1E-15) then
       fail=.true.
       print*,'error in evaluate'
    end if
    
    !test horner for arbitrary polynomials
    call random_seed()
    call random_number(xp)
    res=sll_f_spline_pp_horner_2d(degree, pp_coeffs, xp,[1,1],[1,1])
    res2=0._f64
    do i=1, degree(1)+1
       do j=1, degree(2)+1
          res2=res2+pp_coeffs(i+(j-1)*(degree(1)+1),1)*xp(1)**((degree(1)+1)-i)*xp(2)**((degree(2)+1)-j)
       end do
    end do
    if(abs(res-res2)>1E-12) then
       fail=.true.
       print*, xp
       print*,'error in horner'
    end if

    !print*, 'res-res2=', res-res2
  
   
   
    deallocate(val1)
    deallocate(val2)
    call sll_s_spline_pp_free_2d(spline_pp)   
  end subroutine spline_test
  
end program test_splines_pp_2d
