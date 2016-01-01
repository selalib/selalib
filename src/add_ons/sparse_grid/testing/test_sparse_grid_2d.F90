!-------------------------------------------------------------------
!  test 2D sparse grid
!-------------------------------------------------------------------
program test_interpolation_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_sparse_grid_2d, only: &
    sll_t_sparse_grid_interpolator_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, dimension(:), allocatable :: order
  sll_int32 :: levels
  sll_int32 :: ierr, j,k, counter, iter, its
  sll_real64 :: disp, error
  sll_real64, dimension(:), allocatable :: errorvec,errorvecFFT, tol

  sll_real64, dimension(2)   :: eta_min,eta_max,dx
  sll_real64, dimension(:), allocatable :: f, fref, finterp
  sll_comp64, dimension(:), allocatable :: fftcoeffs

  sll_int32, dimension(2) :: dorder,levelsini
  sll_real64, dimension(:), allocatable :: displacement

  type(sll_t_sparse_grid_interpolator_2d), target   :: interp

  logical :: fail

  ! 
  its = 2;
  ALLOCATE(order(its));
  ALLOCATE(errorvec(4*its));
  ALLOCATE(tol(4*its));
  ALLOCATE(errorvecFFT(its));

  ! Set parameters
  eta_min(1) = 0.0_f64; eta_max(1) = 4.0_f64*sll_p_pi;
  eta_min(2) = 0.0_f64; eta_max(2) = 4.0_f64*sll_p_pi;

  levels = 8;
  order(1) = 1;
  order(2) = 3;

  print*, 'Constant displacement'
  do iter = 1,2
     print*, 'Levels:', levels, ' Order: ' , order(iter)
     levelsini(1) = levels; levelsini(2) = levels;
     ! Initialize sparse grid
     call interp%initialize(levelsini,order(iter), order(iter)+1,0, eta_min, eta_max,0,1);
 
     ! Allocate 
     SLL_ALLOCATE(f(interp%size_basis),ierr);
     SLL_ALLOCATE(fref(interp%size_basis), ierr);
     SLL_ALLOCATE(finterp(interp%size_basis), ierr);
     SLL_ALLOCATE(fftcoeffs(interp%size_basis), ierr);
     
     ! Set displacement
     disp = 0.01_f64;
     
     print*, 'Displacement in eta1'
     ! Set initial value and reference
     do j=1,interp%size_basis
        f(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2));
        fref(j) = sin(0.5_f64*(interp%hierarchy(j)%coordinate(1)+disp)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2));
     end do

    ! Trigonometic sparse grid interpolation
     call interp%SPFFT(f,fftcoeffs)
     call interp%interpolate_array_disp_sgfft(1,disp,fftcoeffs,finterp)
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error constant displacement along one direction (trigonometric):', error
     errorvecFFT(iter) = error;

     
     ! Standard sparse grid interpolation
     
     call interp%compute_hierarchical_surplus(f)
     do j=1,interp%size_basis
        dx(1) = eta_min(1) + modulo(interp%hierarchy(j)%coordinate(1)&
             -eta_min(1)+disp,&
             eta_max(1)-eta_min(1));
        dx(2) = interp%hierarchy(j)%coordinate(2);
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
     end do
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error standard interpolation:', error
     errorvec(4*(iter-1)+1) = error;
     
     
     ! Sparse grid interpolation with interpolate_disp (constant interpolation along one direction)
     call interp%interpolate_disp(1,disp,f, finterp,.FALSE.);
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
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2));
        fref(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*(interp%hierarchy(j)%coordinate(2)+disp));
     end do
     
     ! Standard sparse grid interpolation
     error = 0.0_f64;
     call interp%compute_hierarchical_surplus(f)
     do j=1,interp%size_basis
        dx(2) = eta_min(2) + modulo(interp%hierarchy(j)%coordinate(2)&
             -eta_min(2)+disp,&
             eta_max(2)-eta_min(2));
        dx(1) = interp%hierarchy(j)%coordinate(1);
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error standard interpolation:', error
     errorvec(4*(iter-1)+3) = error;
     
     
     ! Sparse grid interpolation with interpolate_disp (constant interpolation along one direction)
     call interp%interpolate_disp(2,disp,f, finterp,.FALSE.);
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*, 'Non-constant displacement'

 do iter = 1,2
     print*, 'Levels:', levels, ' Order: ' , order(iter)
     levelsini(1) = levels; levelsini(2) = levels;
     ! Initialize sparse grid
     call interp%initialize(levelsini,order(iter), order(iter)+1,0, eta_min, eta_max,0,0);
 
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
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2));
        fref(j) = sin(0.5_f64*(interp%hierarchy(j)%coordinate(1)+&
             disp*interp%hierarchy(j)%coordinate(2))) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2));
     end do
     
     ! Standard sparse grid interpolation
     call interp%compute_hierarchical_surplus(f)
     do j=1,interp%size_basis
        dx(1) = eta_min(1) + modulo(interp%hierarchy(j)%coordinate(1)&
             -eta_min(1)+disp*interp%hierarchy(j)%coordinate(2),&
             eta_max(1)-eta_min(1));
        dx(2) = interp%hierarchy(j)%coordinate(2);
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
     end do
     error = 0.0_f64;
     do j=1,interp%size_basis
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
     dorder(1) = 1; dorder(2) = 2;
     call interp%interpolate_disp_nconst_in_1d(displacement,dorder,&
          f,finterp);
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error constant displacement varying along eta 2:', error
     errorvec(4*(iter-1)+2) = error;
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     print*, 'Displacment in eta 2'
     ! Set initial value and reference
     do j=1,interp%size_basis
        f(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*interp%hierarchy(j)%coordinate(2));
        fref(j) = sin(0.5_f64*interp%hierarchy(j)%coordinate(1)) *&
             cos(0.5_f64*(interp%hierarchy(j)%coordinate(2)+disp*&
             interp%hierarchy(j)%coordinate(1)));
     end do
     
     ! Standard sparse grid interpolation
     error = 0.0_f64;
     call interp%compute_hierarchical_surplus(f)
     do j=1,interp%size_basis
        dx(2) = eta_min(2) + modulo(interp%hierarchy(j)%coordinate(2)&
             -eta_min(2)+disp*interp%hierarchy(j)%coordinate(1),&
             eta_max(2)-eta_min(2));
        dx(1) = interp%hierarchy(j)%coordinate(1);
        finterp(j) = interp%interpolate_from_interpolant_value(f,dx);
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
     dorder(1) = 2; dorder(2) = 1;
     call interp%interpolate_disp_nconst_in_1d(displacement,dorder,&
          f,finterp);
     error = 0.0_f64;
     do j=1,interp%size_basis
        error = max(error,finterp(j)-fref(j));
     end do
     print*, 'Error constant displacement varying along eta 1:', error
     errorvec(4*(iter-1)+2) = error;     
     
     ! Deallocate
     DEALLOCATE(f);
     DEALLOCATE(fref);
     DEALLOCATE(finterp);
     DEALLOCATE(displacement);
  end do

  ! Check accuracy:
  tol(1:3) = 3D-3;
  tol(4:8) = [2D-4, 2D-5, 6D-7, 2.5D-5, 3.2D-9];
  fail = .FALSE.
  do j=1,4*its
     if( errorvec(j) > tol(j)) then
        fail = .TRUE.
     end if
  end do
  do j=1,2
     if (errorvecFFT(j) > 1D-14) then
        fail = .TRUE.
     end if
  end do

  if (fail .EQV. .FALSE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
  end if


end program test_interpolation_2d
