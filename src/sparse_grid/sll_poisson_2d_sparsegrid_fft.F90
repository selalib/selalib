module sll_poisson_2d_sparsegrid_fft

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use sparse_grid_2d
!use  sll_sparsegrid_fft_interpolator_2d
use, intrinsic :: iso_c_binding
  implicit none
 !include 'fftw3.f03'



!private
!public ::new,solve


type,public:: sll_fft_derivative
   sll_real64,dimension(:),pointer :: kx,ky
   sll_real64,dimension(:),pointer :: kpot,kex,key
   sll_comp64, dimension(:),pointer :: fcoeffs,fcoeffs2

contains
   procedure,pass(this) :: initialize=>new_poisson_2d_sparsegrid_fft
   !procedure :: delete=>free_poisson
   procedure,pass(this) :: solve=>solve_for_electric_field!test!poisson_2d_sparsegrid
   !generic   :: solve => solve_for_electric_field, solve_potential
end type sll_fft_derivative
!!$
!!$interface solve
!!$   module procedure solve_potential
!!$   module procedure solve_for_electric_field
!!$end interface
!!$interface delete
!!$   module procedure free_poisson
!!$end interface

contains

subroutine new_poisson_2d_sparsegrid_fft(this,interpolator)
  class(sll_fft_derivative),intent(inout) ::this
  type(sparse_grid_interpolator_2d),intent(in) ::interpolator
  sll_int32 :: ierr,i,j
  sll_real64, dimension(:), allocatable :: data1d
  sll_real64 :: size_factor

  SLL_ALLOCATE(this%kx(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%ky(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%kpot(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%kex(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%key(interpolator%size_basis),ierr)


  
  SLL_ALLOCATE(this%fcoeffs(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%fcoeffs2(interpolator%size_basis),ierr)

  SLL_ALLOCATE(data1d(2**(interpolator%max_level)),ierr)

  size_factor = 2.0_f64*sll_pi/interpolator%length(1)
  do i = 0,interpolator%levels(2)
     do j = interpolator%index(0,i),&
          interpolator%index(0,i) + max(2**(i-1),1)-1
        call derivative_coeffs_1d(interpolator,1,&
             min(interpolator%max_level-i,interpolator%levels(1)),&
             j,size_factor,data1d,this%kx)
     end do
  end do

  size_factor = 2.0_f64*sll_pi/interpolator%length(2)
  do i = 0,interpolator%levels(1)
     do j = interpolator%index(i,0),&
          interpolator%index(i,0) + max(2**(i-1),1)-1
        call derivative_coeffs_1d(interpolator,2,min(interpolator%levels(2),&
             interpolator%max_level-i),&
             j,size_factor,data1d,this%ky)
     end do
  end do


  this%kpot(1) = 0.0;
  this%kex(1) = 0.0;
  this%key(1) = 0.0;
  do i=2,interpolator%size_basis
     this%kpot(i) = 1.0_f64/(this%kx(i)**2+this%ky(i)**2);
     this%kex(i) = this%kpot(i)*this%kx(i);
     this%key(i) = this%kpot(i)*this%ky(i);
  end do

!!$  do i=max(interpolator%level_mapping(interpolator%max_level-2),1),interpolator%level_mapping(interpolator%max_level+1)-1
!!$     this%kex(i) = 0.0_f64;
!!$     this%key(i) = 0.0_f64;
!!$  end do

! TODO: Initialize the data fields used by unidirectional here. 
end subroutine new_poisson_2d_sparsegrid_fft


subroutine solve_potential(this,interpolator,rho,phi)
  class(sll_fft_derivative),intent(inout) ::this
  type(sparse_grid_interpolator_2d), intent(inout)   :: interpolator
  sll_real64,dimension(:),intent(inout) ::phi
  sll_real64,dimension(:),intent(inout) ::rho
  sll_int32 :: i

call SPFFT(interpolator,rho,this%fcoeffs)

do i=1,interpolator%size_basis
   this%fcoeffs(i) =  this%fcoeffs(i)*this%kpot(i)
   !this%fcoeffs(i)*this%kx(i)
end do

call ISPFFT(interpolator,this%fcoeffs,phi)


end subroutine solve_potential


subroutine solve_for_electric_field(this,interpolator,rho,ex,ey)
  class(sll_fft_derivative),intent(inout) ::this
  type(sparse_grid_interpolator_2d), intent(inout)   :: interpolator
  sll_real64,dimension(:),intent(inout) ::ex,ey
  sll_real64,dimension(:),intent(inout) ::rho
  sll_int32 :: i

call SPFFT(interpolator,rho,this%fcoeffs)

do i=1,interpolator%size_basis
   this%fcoeffs2(i) = cmplx(- aimag(this%fcoeffs(i))*this%key(i), &
          real(this%fcoeffs(i))*this%key(i), kind=f64)
   this%fcoeffs(i) = cmplx(- aimag(this%fcoeffs(i))*this%kex(i), &
          real(this%fcoeffs(i))*this%kex(i), kind=f64)
end do

call ISPFFT(interpolator,this%fcoeffs,ex)
call ISPFFT(interpolator,this%fcoeffs2,ey)

end subroutine solve_for_electric_field

!subroutine free_poisson(self)
!type(sll_fft_derivative) :: self
!
!call fft_finalize(self%fft_object,1)
!end subroutine free_poisson

subroutine derivative_coeffs_1d(interpolator,dim,max_level,index,size_factor,data1d,data)
  class(sparse_grid_interpolator_2d), intent(in) :: interpolator
  sll_real64, dimension(:), intent(inout) :: data, data1d
  sll_int32, intent(in) :: dim,max_level,index
  sll_real64 , intent(in) :: size_factor
  sll_int32 :: k, size

  size = 2**max_level

  ! Derivative
  data1d(1) = 0
  data1d(size) = size/2*size_factor
  do k=1,size/2-1
     data1d(2*k) = k*size_factor;
     data1d(2*k+1) = -data1d(2*k);
  end do
  call insert_fourier_real(interpolator,dim,max_level,&
       index,data1d,data)
  !print*, data_out
end subroutine derivative_coeffs_1d


subroutine insert_fourier_real(sparsegrid,dim,max_level,index,data_in,data_out)
  type(sparse_grid_interpolator_2d), intent(in) :: sparsegrid
  sll_int32, intent(in) :: dim,max_level,index
  sll_int32             :: n_points,index_running
  sll_real64,dimension(:),intent(in) :: data_in
  sll_real64,dimension(:),intent(out) :: data_out

  
  n_points = 2**(max_level)
  data_out(index) = data_in(1)
  if (max_level>0) then
     index_running = sparsegrid%hierarchy(index)%children(dim*2)   
     data_out(index_running) = data_in(2)
     if (max_level>1) then
        call insert_recursive_fourier_real(sparsegrid,index_running,0,&
             2,max_level,dim,data_in,data_out)
     end if
  end if
end subroutine insert_fourier_real

recursive subroutine insert_recursive_fourier_real(sparsegrid,index_sg,ind,level,max_level,dim,data_in,data_out)
  type(sparse_grid_interpolator_2d), intent(in) :: sparsegrid
  sll_int32, intent(in) :: level,max_level,index_sg,dim,ind
  sll_real64,dimension(:),intent(in) :: data_in
  sll_real64,dimension(:),intent(inout) :: data_out

  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2-1)) = &
       data_in(2**(level-1)+1+2*ind)
  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2)) = &
       data_in(2**(level-1)+1+2*ind+1)
  if (level<max_level) then
     call insert_recursive_fourier_real(sparsegrid,&
          sparsegrid%hierarchy(index_sg)%children(dim*2-1),2*ind,&
          level+1,max_level,dim,data_in,data_out)
     call insert_recursive_fourier_real(sparsegrid,&
          sparsegrid%hierarchy(index_sg)%children(dim*2),2*ind+1,&
          level+1,max_level,dim,data_in,data_out)
  end if
end subroutine insert_recursive_fourier_real

end module sll_poisson_2d_sparsegrid_fft
