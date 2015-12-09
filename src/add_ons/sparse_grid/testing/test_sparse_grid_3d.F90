!-------------------------------------------------------------------
!  test 3D sparse grid
!-------------------------------------------------------------------

program test_interpolation_3d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_pi

  use sll_m_sparse_grid_3d, only: &
    sparse_grid_interpolator_3d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, dimension(:), allocatable :: order
  sll_int32:: levels
  sll_int32 :: ierr, j, iter, its
  sll_real64 :: disp, error
  sll_real64, dimension(:), allocatable :: errorvec, tol

  sll_real64, dimension(3)   :: eta_min,eta_max,dx
  sll_real64, dimension(:), allocatable :: f, fref, finterp
  sll_comp64, dimension(:), allocatable :: fftcoeffs

  sll_int32, dimension(3) :: dorder, levelsini

  type(sparse_grid_interpolator_3d)   :: interp

  logical :: fail

  ! 
  its = 2;
  ALLOCATE(order(its));
  ALLOCATE(errorvec(4*its));
  ALLOCATE(tol(4*its));

  ! Set parameters
  eta_min(1) = 0.0_f64; eta_max(1) = 4.0_f64*sll_pi;
  eta_min(2) = 0.0_f64; eta_max(2) = 4.0_f64*sll_pi;
  eta_min(3) = 0.0_f64; eta_max(3) = 4.0_f64*sll_pi;

  levels =  8; order(1) = 1; order(2) = 3;

  print*, 'Constant displacement'
  do iter = 1,its
     print*, 'Levels:', levels, ' Order: ' , order(iter)

     levelsini(1) = levels; levelsini(2) = levels;
     levelsini(3) = levels; 
     ! Initialize sparse grid
     call interp%initialize(levelsini,order(iter), order(iter)+1,0, eta_min, eta_max,0,0);
 
     ! Allocate 
     SLL_ALLOCATE(f(interp%size_basis),ierr);
     SLL_ALLOCATE(fref(interp%size_basis), ierr);
     SLL_ALLOCATE(finterp(interp%size_basis), ierr);
     SLL_ALLOCATE(fftcoeffs(interp%size_basis), ierr);
     
     ! Set displacement
     disp = 0.01_f64;
     
     print*, 'Displacement in eta3'
     ! Set initial value and reference
     do j=1,interp%size_basis
        f(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2)) *&
             sin(0.5_f64*interp%hierarchy(j)%coordinate(3));
        fref(j) = sin(0.5_f64*(interp%hierarchy(j)%coordinate(1))) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2))*&
             sin(0.5_f64*(interp%hierarchy(j)%coordinate(3)+disp));
     end do
     
!!$    ! Trigonometic sparse grid interpolation
!!$     call interp%SPFFT(f,fftcoeffs)
!!$     call interp%interpolate_array_disp_sgfft(3,disp,fftcoeffs,finterp)
!!$     error = 0.0_f64;
!!$     do j=1,interp%size_basis
!!$        error = max(error,finterp(j)-fref(j));
!!$     end do
!!$     print*, 'Error constant displacement along one direction (trigonometric):', error

     ! Standard sparse grid interpolation
     error = 0.0_f64;
     call interp%compute_hierarchical_surplus(f)
     do j=1,interp%size_basis
        dx(3) = eta_min(3) + modulo(interp%hierarchy(j)%coordinate(3)&
             -eta_min(3)+disp,&
             eta_max(3)-eta_min(3));
        dx(2) = interp%hierarchy(j)%coordinate(2);
        dx(1) = interp%hierarchy(j)%coordinate(1);
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error standard interpolation:', error
     errorvec(4*(iter-1)+1) = error;
     
     
     ! Sparse grid interpolation with interpolate_disp (constant interpolation along one direction)
     dorder(1) = 3; dorder(2) = 1; dorder(3) = 2;
     !call interp%interpolate_disp(3,disp,f, finterp,.FALSE.);
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
             sin(0.5_f64*interp%hierarchy(j)%coordinate(3));
        fref(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*(interp%hierarchy(j)%coordinate(2)+disp))*&
             sin(0.5_f64*interp%hierarchy(j)%coordinate(3));
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
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error standard interpolation:', error
     errorvec(4*(iter-1)+3) = error;
     
     
     ! Sparse grid interpolation with interpolate_disp (constant interpolation along one direction)
     dorder(1) = 2; dorder(2) = 1; dorder(3) = 3;
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
     DEALLOCATE(fftcoeffs);
  end do

 tol = [2D-3, 2D-3, 2D-3, 2D-3, 2D-4, 2D-5, 1D-4, 6D-6]
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




end program test_interpolation_3d
