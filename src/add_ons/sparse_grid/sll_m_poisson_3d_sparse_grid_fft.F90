!> @ingroup sparse_grid
!> @author Katharina Kormann, IPP 
!> @brief Implementation of a 3D pseudospectral Poisson solver on sparse grid.
!> @details <DETAILED_DESCRIPTION>

module sll_m_poisson_3d_sparse_grid_fft

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_pi

  use sll_m_sparse_grid_3d, only: &
    sparse_grid_interpolator_3d

  implicit none

  public :: &
    sll_fft3d_derivative

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> \a sll_fft3d_derivative is the Poisson solver object to solve Poisson's problem in 2d with pseudospectral on a sparse grid
type:: sll_fft3d_derivative
   sll_real64,dimension(:),pointer :: kx  !< Fourier coefficients for first derivative along x
   sll_real64,dimension(:),pointer :: ky !< Fourier coefficients for first derivative along y   
   sll_real64,dimension(:),pointer :: kz !< Fourier coefficients for first derivative along z
   sll_real64,dimension(:),pointer :: kpot !< Fourier coefficient to compute potential
   sll_real64,dimension(:),pointer :: kex !< Fourier coefficients to compute E_x
   sll_real64,dimension(:),pointer :: key !< Fourier coefficients to compute E_y
   sll_real64,dimension(:),pointer :: kez !< Fourier coefficients to compute E_z
   sll_comp64, dimension(:),pointer :: fcoeffs,fcoeffs2 !< Internal array for SGFFT computations

contains
   procedure :: initialize=>new_poisson_3d_sparse_grid_fft
   !procedure :: delete=>free_poisson
   procedure :: solve=>solve_for_electric_field
   procedure :: solve_potential
end type sll_fft3d_derivative

contains

!> Create Poisson solver object with Fourier spectral method on 3d sparse grid
subroutine new_poisson_3d_sparse_grid_fft(this,interpolator)
  class(sll_fft3d_derivative),intent(inout) ::this !< Poisson solver object
  type(sparse_grid_interpolator_3d),intent(in) ::interpolator !< Underlying sparse grid

  sll_int32 :: ierr,i1,i2,i3,i,j
  sll_real64, dimension(:), allocatable :: data1d
  sll_real64 :: size_factor

  SLL_ALLOCATE(this%kx(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%ky(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%kz(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%kpot(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%kex(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%key(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%kez(interpolator%size_basis),ierr)

  
  SLL_ALLOCATE(this%fcoeffs(interpolator%size_basis),ierr)
  SLL_ALLOCATE(this%fcoeffs2(interpolator%size_basis),ierr)

  SLL_ALLOCATE(data1d(2**(interpolator%max_level)),ierr)


  size_factor = 2.0_f64*sll_pi/interpolator%length(1)
  do i2 = 0,interpolator%levels(2)
     do i3 =0,min(interpolator%max_level - i2,interpolator%levels(3))
        do j = interpolator%index(0,i2,i3),&
             interpolator%index(0,i2,i3)+max(2**(i2-1),1)*max(2**(i3-1),1)-1
           call derivative_coeffs_1d(interpolator,1,&
                min(interpolator%levels(1),interpolator%max_level-i2-i3),&
                j,size_factor,data1d,this%kx)
        end do
     end do
  end do

  size_factor = 2.0_f64*sll_pi/interpolator%length(2)
  do i1 = 0,interpolator%levels(1)
     do i3=0,min(interpolator%max_level - i1,interpolator%levels(3))
        do j = interpolator%index(i1,0,i3),&
             interpolator%index(i1,0,i3)+max(2**(i1-1),1)*max(2**(i3-1),1)-1
           call derivative_coeffs_1d(interpolator,2,&
                min(interpolator%levels(2),interpolator%max_level-i1-i3),j,&
                size_factor,data1d,this%ky)
        end do
     end do
  end do

  size_factor = 2.0_f64*sll_pi/interpolator%length(3)
  do i1 = 0,interpolator%levels(1)
     do i2 = 0, min(interpolator%max_level - i1,interpolator%levels(2))
        do j = interpolator%index(i1,i2,0),&
             interpolator%index(i1,i2,0)+max(2**(i1-1),1)*max(2**(i2-1),1)-1
           call derivative_coeffs_1d(interpolator,3,&
                min(interpolator%levels(3),interpolator%max_level-i1-i2),j,&
                size_factor,data1d,this%kz)
        end do
     end do
  end do

  this%kpot(1) = 0.0_f64;
  this%kex(1) = 0.0_f64;
  this%key(1) = 0.0_f64;
  this%kez(1) = 0.0_f64;
  do i=2,interpolator%size_basis
     this%kpot(i) = 1.0_f64/(this%kx(i)**2+this%ky(i)**2+this%kz(i)**2);
     this%kex(i) = this%kpot(i)*this%kx(i);
     this%key(i) = this%kpot(i)*this%ky(i);
     this%kez(i) = this%kpot(i)*this%kz(i);
  end do

end subroutine new_poisson_3d_sparse_grid_fft

!> Solve for potential
subroutine solve_potential(this,interpolator,rho,phi)
  class(sll_fft3d_derivative),intent(inout) ::this !< Poisson solver object
  type(sparse_grid_interpolator_3d), intent(inout)   :: interpolator !< Underlying sparse grid
  sll_real64,dimension(:),intent(inout) ::phi !< Solution of Poisson's equation
  sll_real64,dimension(:),intent(inout) ::rho !< Right-hand-side for Poisson's equation
  sll_int32 :: i

call interpolator%SPFFT(rho,this%fcoeffs)

do i=1,interpolator%size_basis
   this%fcoeffs(i) =  this%fcoeffs(i)*this%kpot(i)
end do

call interpolator%ISPFFT(this%fcoeffs,phi)


end subroutine solve_potential

!> Compute the electric fields from rho
subroutine solve_for_electric_field(this,interpolator,rho,ex,ey,ez)
  class(sll_fft3d_derivative),intent(inout) ::this !< Poisson solver object
  type(sparse_grid_interpolator_3d), intent(inout)   :: interpolator !< underlying sparse grid
  sll_real64, dimension(:), intent(inout) :: ex !< x component of electric field
  sll_real64, dimension(:), intent(inout) :: ey !< y component of electric field
  sll_real64, dimension(:), intent(inout) :: ez !< z component of electric field
  sll_real64, dimension(:), intent(inout) :: rho !< given density as rhs

  sll_int32 :: i

call interpolator%SPFFT(rho,this%fcoeffs)

do i=1,interpolator%size_basis
   this%fcoeffs2(i) = cmplx(- aimag(this%fcoeffs(i))*this%kez(i), &
          real(this%fcoeffs(i))*this%kez(i), kind=f64)
end do

call interpolator%ISPFFT(this%fcoeffs2,ez)

do i=1,interpolator%size_basis
   this%fcoeffs2(i) = cmplx(- aimag(this%fcoeffs(i))*this%key(i), &
          real(this%fcoeffs(i))*this%key(i), kind=f64)
   this%fcoeffs(i) = cmplx(- aimag(this%fcoeffs(i))*this%kex(i), &
          real(this%fcoeffs(i))*this%kex(i), kind=f64)
end do

call interpolator%ISPFFT(this%fcoeffs,ex)
call interpolator%ISPFFT(this%fcoeffs2,ey)

end subroutine solve_for_electric_field

!> Helper function to compute the Fourier coefficients for the derivative operation
subroutine derivative_coeffs_1d(interpolator,dim,max_level,index,size_factor,data1d,data)
  class(sparse_grid_interpolator_3d), intent(in) :: interpolator
  sll_real64, dimension(:), intent(inout) :: data, data1d
  sll_int32, intent(in) :: dim,max_level,index
  sll_real64 , intent(in) :: size_factor
  sll_int32 :: k, size

  size = 2**max_level

  ! Derivative
  data1d(1) = 0.0_f64
  data1d(size) = size/2*size_factor
  do k=1,size/2-1
     data1d(2*k) = k*size_factor;
     data1d(2*k+1) = -data1d(2*k);
  end do
  call insert_fourier_real(interpolator,dim,max_level,&
       index,data1d,data)
end subroutine derivative_coeffs_1d

!> Helper function to insert the real Fourier coefficient into the sparse grid data structure
subroutine insert_fourier_real(sparsegrid,dim,max_level,index,data_in,data_out)
  type(sparse_grid_interpolator_3d), intent(in) :: sparsegrid
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

!> Helper function to insert the real Fourier coefficient into the sparse grid data structure (recursive part)
recursive subroutine insert_recursive_fourier_real(sparsegrid,index_sg,ind,level,max_level,dim,data_in,data_out)
  type(sparse_grid_interpolator_3d), intent(in) :: sparsegrid
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

end module sll_m_poisson_3d_sparse_grid_fft
