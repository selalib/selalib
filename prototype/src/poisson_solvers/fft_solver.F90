program fft_solver
implicit none
#include "fftw3.f"

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                   Variable declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer(8)               :: forward, backward
PetscInt                 :: m = 512, n = 512
real(8), allocatable     :: farray(:,:), uarray(:,:)

allocate(farray(0:m+1,n+1), uarray(0:m+1,n))
uarray(:,:) = farray(:,1:n) / n
!
!Forward FFT on the right hand side term
!
call dfftw_plan_r2r_1d(forward,n,uarray(1,1:n),farray(1,1:n),FFTW_R2HC,FFTW_ESTIMATE )
do i = 1, m
   call dfftw_execute_r2r(forward,uarray(i,1:n),farray(i,1:n));
end do
call dfftw_destroy_plan(forward)
!
!  Backward FFT on the solution
!
call dfftw_plan_r2r_1d(backward,n,uarray(1,1:n),farray(1,1:n),FFTW_HC2R,FFTW_ESTIMATE)
do i = 1, m
   call dfftw_execute_r2r(backward,uarray(i,1:n),farray(i,1:n))
end do
call dfftw_destroy_plan(backward)

end program fft_solver
