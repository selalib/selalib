program test_splines_pp_3d
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

  type(sll_t_spline_pp_3d) :: spline_pp
  sll_int32 :: degree(3)
  sll_int32 :: n_cells(3)
  logical   :: fail
 
  n_cells=[10,10,10]
  degree=[3,3,3]
  
  call spline_test(spline_pp,degree,n_cells,fail)
     
  
  if(fail .eqv. .false.) then
     write(*,*) 'PASSED'
  else  
     write(*,*)'FAILED'
     stop
  end if

contains

  subroutine spline_test(spline_pp,degree,n_cells,fail)

    type(sll_t_spline_pp_3d), intent(inout) :: spline_pp !arbitrary degree spline
    sll_int32, intent(in) :: degree(3) !spline degree
    sll_int32, intent(in) :: n_cells(3) !grid cells
    logical, intent(out)    :: fail
    
    sll_real64 :: b_coeffs(n_cells(1)*n_cells(2)*n_cells(3))
    sll_real64 :: pp_coeffs((degree(1)+1)*(degree(2)+1)*(degree(3)+1),n_cells(1)*n_cells(2)*n_cells(3))
    sll_real64 :: xp(3) !particle position
    sll_real64 :: res
    sll_real64 :: res2
    sll_real64 :: domain(2,3) !left and right bound of the 1d grid
    sll_real64 :: delta_x(3) !size of gridcells
    sll_real64, allocatable :: val1(:)
    sll_real64, allocatable :: val2(:)
    sll_real64, allocatable :: val3(:)
    sll_int32 :: index3d(3)
    sll_int32 :: index1d
    sll_int32 :: indices(3)
    sll_real64 :: xi(3)
    sll_int32 :: i,j,k
    
    fail=.false.
    
    allocate(val1(degree(1)+1))
    allocate(val2(degree(2)+1))
    allocate(val3(degree(3)+1))
    domain(1,:)=0._f64
    domain(2,:)=sll_p_twopi 
    delta_x=(domain(2,:)-domain(1,:))/real(n_cells,f64)

    call random_seed()
    call random_number(b_coeffs)

    call sll_s_spline_pp_init_3d( spline_pp, degree, n_cells)
    call sll_s_spline_pp_b_to_pp_3d(spline_pp, n_cells, b_coeffs, pp_coeffs)
 

    call random_seed()
    call random_number(xp)
    !xp=[0._f64,1._f64,1._f64]
   
    xp=xp*sll_p_twopi
    ! calculate index of gridcell
    xi = (xp - domain(1,:))/delta_x
    indices = floor(xi)+1
    xi = xi - real(indices-1, f64)
    res= sll_f_spline_pp_horner_3d(degree, pp_coeffs, xi, indices, n_cells)
    indices = indices - degree
    
    call sll_s_uniform_b_splines_at_x(degree(1), xi(1), val1(:))
    call sll_s_uniform_b_splines_at_x(degree(2), xi(2), val2(:))
    call sll_s_uniform_b_splines_at_x(degree(3), xi(3), val3(:))

        
    res2 = 0.0_f64
    do i = 1, degree(1)+1
       index3d(1) = modulo(indices(1)+i-2, n_cells(1)) 
       do j = 1, degree(2)+1
          index3d(2) = modulo(indices(2)+j-2, n_cells(2))
          do k = 1, degree(3)+1
             index3d(3) = modulo(indices(3)+k-2, n_cells(3))
             index1d = index3d(1) + 1 + index3d(2)*n_cells(1) + index3d(3)*n_cells(1)*n_cells(2)
             res2 = res2 + b_coeffs(index1d) * val1(i) * val2(j) * val3(k)
          end do
       end do
    end do
    !print*, 'res-res2=', res-res2
    !write(*,*) 'Fehler horner vs normal:', abs(res-res2)
    if(abs(res-res2)>1E-13) then
       fail=.true.
       print*,'error in evaluate'
    end if

 
 
    !test horner for arbitrary polynomials
    call random_seed()
    call random_number(xp)
    res=sll_f_spline_pp_horner_3d(degree, pp_coeffs, xp,[1,1,1],[1,1,1])
    res2=0._f64
    do i=1, degree(1)+1
       do j=1, degree(2)+1
          do k=1, degree(3)+1
             res2=res2+pp_coeffs(i+(j-1)*(degree(1)+1)+(k-1)*(degree(1)+1)*(degree(2)+1),1)*xp(1)**((degree(1)+1)-i)*xp(2)**((degree(2)+1)-j)*xp(3)**((degree(3)+1)-k)
          end do
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
    deallocate(val3)
    call sll_s_spline_pp_free_3d(spline_pp)
  end subroutine spline_test
  
end program test_splines_pp_3d
