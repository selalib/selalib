!-------------------------------------------------------------------
!  test 4D sparse grid
!-------------------------------------------------------------------
program test_interpolation_4d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_sparse_grid_2d, only: &
    sll_t_sparse_grid_interpolator_2d

  use sll_m_sparse_grid_4d, only: &
    sll_t_sparse_grid_interpolator_4d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, dimension(:), allocatable :: order
  sll_int32 :: levels
  sll_int32 :: ierr, j,k, counter, iter, its
  sll_real64 :: disp, error, dispt
  sll_real64, dimension(:), allocatable :: errorvec

  sll_real64, dimension(4)   :: eta_min,eta_max,dx
  sll_real64, dimension(:), allocatable :: f, fref, finterp

  sll_int32, dimension(4) :: dorder, levelsini
  sll_real64, dimension(:), allocatable :: displacement
  sll_real64, dimension(:), allocatable :: displacement2d

  type(sll_t_sparse_grid_interpolator_4d), target   :: interp
  type(sll_t_sparse_grid_interpolator_2d), target   :: interp2d

  sll_real64 :: tol(8)
  logical :: fail

  ! 
  its = 2;
  ALLOCATE(order(its));
  ALLOCATE(errorvec(4*its));

  ! Set parameters
  eta_min(1) = 0.0_f64; eta_max(1) = 4.0_f64*sll_p_pi;
  eta_min(2) = 0.0_f64; eta_max(2) = 4.0_f64*sll_p_pi;
  eta_min(3) = 0.0_f64; eta_max(3) = 4.0_f64*sll_p_pi;
  eta_min(4) = 0.0_f64; eta_max(4) = 4.0_f64*sll_p_pi;

  levels =  8; 
  order(1) = 1;
  order(2) = 3;
  print*, 'Constant displacement'
  do iter = 1,its
     print*, 'Levels:', levels, ' Order: ' , order(iter)

     levelsini(1) = levels; levelsini(2) = levels;
     levelsini(3) = levels; levelsini(4) = levels;
     ! Initialize sparse grid
     call interp%initialize(levelsini,order(iter), order(iter)+1,0, eta_min, eta_max);
 
     ! Allocate 
     SLL_ALLOCATE(f(interp%size_basis),ierr);
     SLL_ALLOCATE(fref(interp%size_basis), ierr);
     SLL_ALLOCATE(finterp(interp%size_basis), ierr);
     
     ! Set displacement
     disp = 0.01_f64;
     
     print*, 'Displacement in eta1'
     ! Set initial value and reference
     do j=1,interp%size_basis
        f(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2)) *&
             sin(0.5_f64*interp%hierarchy(j)%coordinate(3)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(4));
        fref(j) = sin(0.5_f64*(interp%hierarchy(j)%coordinate(1))) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2))*&
             sin(0.5_f64*interp%hierarchy(j)%coordinate(3)) *&
             cos(0.5_f64*(interp%hierarchy(j)%coordinate(4)+disp));
     end do
     
     ! Standard sparse grid interpolation
     error = 0.0_f64;
     call interp%compute_hierarchical_surplus(f)
     do j=1,interp%size_basis
        dx(4) = eta_min(4) + modulo(interp%hierarchy(j)%coordinate(4)&
             -eta_min(4)+disp,&
             eta_max(4)-eta_min(4));
        dx(2) = interp%hierarchy(j)%coordinate(2);
        dx(3) = interp%hierarchy(j)%coordinate(3);
        dx(1) = interp%hierarchy(j)%coordinate(1);
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error standard interpolation:', error
     errorvec(4*(iter-1)+1) = error;
     
     
     ! Sparse grid interpolation with interpolate_disp (constant interpolation along one direction)
     dorder(1) = 4; dorder(2) = 1; dorder(3) = 2; dorder(4) = 3;
     call interp%interpolate_const_disp(dorder,disp,f, finterp,.FALSE.);
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error constant displacement along one direction:', error
     errorvec(4*(iter-1)+2) = error;
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     print*, 'Displacment in eta 2'
     ! Set initial value and reference
     do j=1,interp%size_basis
        f(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2))*&
             sin(0.5_f64*interp%hierarchy(j)%coordinate(3)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(4));
        fref(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*(interp%hierarchy(j)%coordinate(2)+disp))*&
             sin(0.5_f64*interp%hierarchy(j)%coordinate(3)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(4));
     end do
     
     ! Standard sparse grid interpolation
     error = 0.0_f64;
     call interp%compute_hierarchical_surplus(f)
     do j=1,interp%size_basis
        dx(2) = eta_min(2) + modulo(interp%hierarchy(j)%coordinate(2)&
             -eta_min(2)+disp,&
             eta_max(2)-eta_min(2));
        dx(1) = interp%hierarchy(j)%coordinate(1);
        dx(3) = interp%hierarchy(j)%coordinate(3);
        dx(4) = interp%hierarchy(j)%coordinate(4);
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error standard interpolation:', error
     errorvec(4*(iter-1)+3) = error;
     
     
     ! Sparse grid interpolation with interpolate_disp (constant interpolation along one direction)
     dorder(1) = 2; dorder(2) = 1; dorder(3) = 3; dorder(4) = 4;
     call interp%interpolate_const_disp(dorder,disp,f, finterp,.FALSE.);
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error constant displacement along one direction:', error
     errorvec(4*(iter-1)+4) = error;
     
     ! Deallocate
     DEALLOCATE(f);
     DEALLOCATE(fref);
     DEALLOCATE(finterp);
  end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*, 'Displacement non-constant in 1d'

 do iter = 1,its
     print*, 'Levels:', levels, ' Order: ' , order(iter)

     levelsini(1) = levels; levelsini(2) = levels;
     levelsini(3) = levels; levelsini(4) = levels;
     ! Initialize sparse grid
     call interp%initialize(levelsini,order(iter), order(iter)+1,0, eta_min, eta_max);
 
     ! Allocate 
     SLL_ALLOCATE(f(interp%size_basis),ierr);
     SLL_ALLOCATE(fref(interp%size_basis), ierr);
     SLL_ALLOCATE(finterp(interp%size_basis), ierr);
     
     ! Set displacement
     disp = 0.01_f64;
     
     print*, 'Displacement in eta1'
     ! Set initial value and reference
     do j=1,interp%size_basis
        f(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2))*&
             sin(0.5_f64*interp%hierarchy(j)%coordinate(3)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(4));
        fref(j) = sin(0.5_f64*(interp%hierarchy(j)%coordinate(1)+&
             disp*interp%hierarchy(j)%coordinate(2))) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2))*&
             sin(0.5_f64*interp%hierarchy(j)%coordinate(3)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(4));
     end do
     
     ! Standard sparse grid interpolation
     error = 0.0_f64;
     call interp%compute_hierarchical_surplus(f)
     do j=1,interp%size_basis
        dx(1) = eta_min(1) + modulo(interp%hierarchy(j)%coordinate(1)&
             -eta_min(1)+disp*interp%hierarchy(j)%coordinate(2),&
             eta_max(1)-eta_min(1));
        dx(2) = interp%hierarchy(j)%coordinate(2);
        dx(3) = interp%hierarchy(j)%coordinate(3);
        dx(4) = interp%hierarchy(j)%coordinate(4);
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error standard interpolation:', error
     errorvec(4*(iter-1)+1) = error;
     
     
     ! Sparse grid interpolation with interpolate_disp (constant interpolation along one direction)
     ALLOCATE(displacement(2**levels));
     counter = 2;
     displacement(1) = eta_min(2);
     do j=1,levels
        do k=1,2**(j-1)
           displacement(counter) = ((real(2*k,f64)-1.0_f64)/real(2**j,f64)*&
                interp%length(2)-eta_min(2))*disp;
           counter = counter+1;
        end do
     end do
     dorder(1) = 1; dorder(2) = 2; dorder(3) = 3; dorder(4) = 4;
     call interp%interpolate_disp_nconst_in_1d(displacement,dorder,&
          f,finterp);
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error constant displacement varying along eta 2:', error
     errorvec(4*(iter-1)+2) = error;
     
!!$     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     print*, 'Displacment in eta 2'
     ! Set initial value and reference
     do j=1,interp%size_basis
        f(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(4)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2));
        fref(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(4)) *&
             cos(0.5_f64*(interp%hierarchy(j)%coordinate(2)+disp*&
             interp%hierarchy(j)%coordinate(4)));
     end do
     
     ! Standard sparse grid interpolation
     call interp%compute_hierarchical_surplus(f)
     do j=1,interp%size_basis
        dx(2) = eta_min(2) + modulo(interp%hierarchy(j)%coordinate(2)&
             -eta_min(2)+disp*interp%hierarchy(j)%coordinate(4),&
             eta_max(2)-eta_min(2));
        dx(1) = interp%hierarchy(j)%coordinate(1);
        dx(3) = interp%hierarchy(j)%coordinate(3);
        dx(4) = interp%hierarchy(j)%coordinate(4);
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
     end do
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error standard interpolation:', error
     errorvec(4*(iter-1)+3) = error;
     
     ! Sparse grid interpolation with interpolate_disp (constant interpolation along one direction)
     counter = 2;
     displacement(1) = eta_min(1);
     do j=1,levels
        do k=1,2**(j-1)
           displacement(counter) = ((real(2*k,f64)-1.0_f64)/real(2**j,f64)*&
                interp%length(1)-eta_min(1))*disp;
           counter = counter+1;
        end do
     end do
     dorder(1) = 2; dorder(2) = 4; dorder(3) = 1; dorder(4) = 3;
     call interp%interpolate_disp_nconst_in_1d(displacement,dorder,&
          f,finterp);
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error constant displacement varying along eta 4:', error
     errorvec(4*(iter-1)+4) = error;     


     ! Deallocate
     DEALLOCATE(f);
     DEALLOCATE(fref);
     DEALLOCATE(finterp);
     DEALLOCATE(displacement);
  end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*, 'Displacement non-constant in 2d'

 do iter = 1,its
     print*, 'Levels:', levels, ' Order: ' , order(iter)

     levelsini(1) = levels; levelsini(2) = levels;
     levelsini(3) = levels; levelsini(4) = levels;
     ! Initialize sparse grid
     call interp%initialize(levelsini,order(iter), order(iter)+1,0,&
          eta_min, eta_max);
     call interp2d%initialize(levelsini(1:2),order(iter), order(iter)+1,0, &
          eta_min(1:2), eta_max(1:2),0,0);
 
     ! Allocate 
     SLL_ALLOCATE(f(interp%size_basis),ierr);
     SLL_ALLOCATE(fref(interp%size_basis), ierr);
     SLL_ALLOCATE(finterp(interp%size_basis), ierr);
     SLL_ALLOCATE(displacement2d(interp2d%size_basis),ierr);
     
     ! Set displacement
     disp = 0.01_f64;
     do j=1,interp2d%size_basis
        displacement2d(j) = (interp2d%hierarchy(j)%coordinate(1) +&
             interp2d%hierarchy(j)%coordinate(2))*disp;
     end do
     
     print*, 'Displacement in eta1'
     ! Set initial value and reference
     do j=1,interp%size_basis
        f(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2)) *&
             sin(0.5_f64*interp%hierarchy(j)%coordinate(3)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(4));
        dispt = (interp%hierarchy(j)%coordinate(1) +&
             interp%hierarchy(j)%coordinate(2))*disp;
        fref(j) = sin(0.5_f64*(interp%hierarchy(j)%coordinate(1))) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2))*&
             sin(0.5_f64*(interp%hierarchy(j)%coordinate(3)+dispt)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(4));
     end do
     
     call interp%compute_hierarchical_surplus(f)
     do j=1,interp%size_basis
        dispt = (interp%hierarchy(j)%coordinate(1) +&
             interp%hierarchy(j)%coordinate(2))*disp;
        dx(3) = eta_min(3) + modulo(interp%hierarchy(j)%coordinate(3)&
             -eta_min(3)+dispt,&
             eta_max(3)-eta_min(3));
        dx(2) = interp%hierarchy(j)%coordinate(2);
        dx(1) = interp%hierarchy(j)%coordinate(1);
        dx(4) = interp%hierarchy(j)%coordinate(4);
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
     end do
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error standard interpolation:', error
     errorvec(4*(iter-1)+1) = error;
     
     ! Sparse grid interpolation with displacement depending on two variables.
     dorder(1) = 3; dorder(2) = 1; dorder(3) = 2; dorder(4) = 4;
     call interp%interpolate_disp_nconst_in_2d(displacement2d,dorder,&
          f, finterp);
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error constant displacement varying along eta 2:', error
     errorvec(4*(iter-1)+2) = error;
         


     ! Deallocate
     DEALLOCATE(f);
     DEALLOCATE(fref);
     DEALLOCATE(finterp);
     DEALLOCATE(displacement2d);
  end do

  tol = [4D-2, 4D-2, 2.5D-3, 2.5D-3, 3D-2, 2D-2, 2.5D-5, 6D-7]
  fail = .FALSE.
  do j=1,4*its
     if( errorvec(j) > tol(j)) then
        fail = .TRUE.
     end if
  end do

  if (fail .EQV. .FALSE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
  end if


end program test_interpolation_4d
