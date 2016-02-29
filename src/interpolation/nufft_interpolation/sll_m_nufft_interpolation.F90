module sll_m_nufft_interpolation
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_m_fft
use sll_m_constants

implicit none

!> Nufft object for 2d interpolation.
!> It contains fft plan and 1d array to pass data to
!> nufft2d subroutine from nufft package
type :: sll_t_nufft_2d
  private
  type(sll_t_fft)         :: fft
  sll_comp64, allocatable :: f1d(:)
  sll_comp64, allocatable :: f1(:)
  sll_comp64, allocatable :: f2(:)
  sll_comp64, allocatable :: fcmplx(:,:)
  sll_int32,  allocatable :: i(:)
  sll_int32,  allocatable :: j(:)
  sll_real64, allocatable :: x(:)
  sll_real64, allocatable :: y(:)
  sll_real64              :: epsnufft=1.0d-9
  sll_int32               :: nc_eta1
  sll_int32               :: nc_eta2
  sll_real64              :: eta1_min, eta1_max
  sll_real64              :: eta2_min, eta2_max

end type sll_t_nufft_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Allocate and initialize data and prepare the fft plan
subroutine sll_s_nufft_2d_init( self,                        &
                                nc_eta1, eta1_min, eta1_max, &
                                nc_eta2, eta2_min, eta2_max)

type(sll_t_nufft_2d)   :: self      !< nufft 2d object
sll_int32,  intent(in) :: nc_eta1   !< number of cells on x1
sll_int32,  intent(in) :: nc_eta2   !< number of cells on x2
sll_real64, intent(in) :: eta1_min  !< left 
sll_real64, intent(in) :: eta1_max  !< right 
sll_real64, intent(in) :: eta2_min  !< bottom
sll_real64, intent(in) :: eta2_max  !< top
sll_int32              :: err

self%nc_eta1  = nc_eta1      
self%nc_eta2  = nc_eta2      

self%eta1_min = eta1_min
self%eta2_min = eta2_min
self%eta1_max = eta1_max
self%eta2_max = eta2_max

SLL_ALLOCATE(self%x       (1:nc_eta1*nc_eta2),   err)
SLL_ALLOCATE(self%y       (1:nc_eta1*nc_eta2),   err)
SLL_ALLOCATE(self%i       (1:nc_eta1*nc_eta2),   err)
SLL_ALLOCATE(self%j       (1:nc_eta1*nc_eta2),   err)
SLL_ALLOCATE(self%fcmplx  (1:nc_eta1,1:nc_eta2), err)
SLL_ALLOCATE(self%f1      (1:nc_eta1),           err)
SLL_ALLOCATE(self%f2      (1:nc_eta2),           err)
SLL_ALLOCATE(self%f1d     (1:nc_eta1*nc_eta2),   err)

!$OMP CRITICAL
call sll_s_fft_init_c2c_2d(self%fft,nc_eta1,nc_eta2, &
  self%fcmplx,self%fcmplx,sll_p_fft_forward)
!$OMP END CRITICAL

end subroutine sll_s_nufft_2d_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Delete the nufft object.
subroutine sll_s_nufft_2d_free( self)

type(sll_t_nufft_2d) :: self

  DEALLOCATE(self%f1d )
  DEALLOCATE(self%x )
  DEALLOCATE(self%y )
  DEALLOCATE(self%j  )
  DEALLOCATE(self%i  )
  DEALLOCATE(self%fcmplx  )
  DEALLOCATE(self%f1      )
  DEALLOCATE(self%f2      )

end subroutine sll_s_nufft_2d_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Compute the fft and prepare data for nufft call.
subroutine sll_s_nufft_2d_compute_fft( self, f_in )

type(sll_t_nufft_2d)                 :: self
sll_real64, intent(in)               :: f_in(:,:)

sll_int32                            :: i, j, m, n1, n2, p

n1 = self%nc_eta1
n2 = self%nc_eta2

self%fcmplx = cmplx(f_in(1:n1,1:n2),0.,f64)
call sll_s_fft_exec_c2c_2d(self%fft, self%fcmplx, self%fcmplx)
self%fcmplx = self%fcmplx / cmplx(n1*n2,0.,f64)

m = n2/2
do i=1,n1
  self%f2=self%fcmplx(i,:)
  self%fcmplx(i,:)=self%f2([[(p,p=m+1,n2)],[(p,p=1,m)]] )
enddo
m = n1/2
do j=1,n2
  self%f1=self%fcmplx(:,j)
  self%fcmplx(:,j)=self%f1([[(p,p=m+1,n1)], [(p,p=1,m)]])
enddo

end subroutine sll_s_nufft_2d_compute_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Interpolate single value when the fft is already compute.
function sll_f_nufft_2d_interpolate_value_from_fft( self, x, y ) result(f_out)

type(sll_t_nufft_2d)   :: self
sll_real64, intent(in) :: x
sll_real64, intent(in) :: y
sll_real64             :: f_out
sll_int32              :: error
sll_real64             :: xij, yij

xij = (x - self%eta1_min) / (self%eta1_max-self%eta1_min)
yij = (y - self%eta2_min) / (self%eta2_max-self%eta2_min)

if ( 0.0_f64 < xij .and. xij < 1.0_f64 .and. &
     0.0_f64 < yij .and. yij < 1.0_f64 ) then

  self%i(1) = 1
  self%j(1) = 1
  self%x(1) = xij * 2.0 * sll_p_pi
  self%y(1) = yij * 2.0 * sll_p_pi

  call nufft2d2f90(1,               &
                   self%x(1), &
                   self%y(1), &
                   self%f1d(1),     &
                   1,               &
                   self%epsnufft,   &
                   1,               &
                   1,               &
                   self%fcmplx,     &
                   error)

  f_out = real(self%f1d(1))

else

  f_out = 0.0_f64

end if

end function sll_f_nufft_2d_interpolate_value_from_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Compute the fft and interpolate array values. 
subroutine sll_s_nufft_2d_interpolate_array_values( self, f_in, x, y, f_out )

type(sll_t_nufft_2d)           :: self
sll_real64, intent(in)         :: f_in(:,:)
sll_real64, intent(in)         :: x(:,:)
sll_real64, intent(in)         :: y(:,:)
sll_real64, intent(out)        :: f_out(:,:)

call sll_s_nufft_2d_compute_fft( self, f_in )
f_out = f_in
call sll_s_nufft_2d_interpolate_array_from_fft( self, x, y, f_out ) 

end subroutine sll_s_nufft_2d_interpolate_array_values

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Compute the fft and interpolate array values when x and y could be
!> outside the domaine.  Useful for rotation and axisymetric geometry.
subroutine sll_s_nufft_2d_interpolate_array_values_axi( self, f, x, y ) 

type(sll_t_nufft_2d)           :: self
sll_real64, intent(inout)      :: f(:,:)
sll_real64, intent(in)         :: x(:,:)
sll_real64, intent(in)         :: y(:,:)
sll_int32                      :: i, j, k, error
sll_int32                      :: n1, n2,  p
sll_real64                     :: xij, yij

n1 = self%nc_eta1
n2 = self%nc_eta2

call sll_s_nufft_2d_compute_fft( self, f )

p  = 0
do j=1,n2
  do i=1,n1
    xij = (x(i,j) - self%eta1_min) / (self%eta1_max-self%eta1_min)
    yij = (y(i,j) - self%eta2_min) / (self%eta2_max-self%eta2_min)
    if ( 0.0_f64 < xij .and. xij < 1.0_f64 .and. &
         0.0_f64 < yij .and. yij < 1.0_f64 ) then
      p = p + 1
      self%i(p)  = i
      self%j(p)  = j
      self%x(p)  = xij * 2.0 * sll_p_pi
      self%y(p)  = yij * 2.0 * sll_p_pi
    else
      f(i,j) = 0.0_f64
    end if
  enddo
enddo

call nufft2d2f90(p,                 &
                 self%x(1:p),       &
                 self%y(1:p),       &
                 self%f1d(1:p),     &
                 1,                 &
                 self%epsnufft,     &
                 n1,                &
                 n2,                &
                 self%fcmplx,       &
                 error)

do k=1,p
  f(self%i(k),self%j(k))=real(self%f1d(k))
end do

end subroutine sll_s_nufft_2d_interpolate_array_values_axi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Compute the fft and interpolate array values when x and y are 
!> surely inside the domain. Optimized version for advector  
subroutine sll_s_nufft_2d_interpolate_array_from_fft( self, x, y, f ) 

type(sll_t_nufft_2d)    :: self
sll_real64, intent(in)  :: x(:,:)
sll_real64, intent(in)  :: y(:,:)
sll_real64, intent(out) :: f(:,:)
sll_int32               :: n1, n2
sll_int32               :: i, j
sll_int32               :: error

n1 = self%nc_eta1
n2 = self%nc_eta2

SLL_ASSERT(size(f,1) == n1 .and. size(f,2) == n2)
SLL_ASSERT(size(x,1) == n1 .and. size(x,2) == n2)
SLL_ASSERT(size(y,1) == n1 .and. size(y,2) == n2)

do j = 1, n2
  do i = 1, n1
    self%x((j-1)*n1+i) = (-self%eta1_min + x(i,j)) * 2.0 * sll_p_pi / (self%eta1_max-self%eta1_min)
    self%y((j-1)*n1+i) = (-self%eta2_min + y(i,j)) * 2.0 * sll_p_pi / (self%eta2_max-self%eta2_min)
  end do
end do

call nufft2d2f90(n1*n2,             &
                 self%x,            &
                 self%y,            &
                 self%f1d,          &
                 1,                 &
                 self%epsnufft,     &
                 n1,                &
                 n2,                &
                 self%fcmplx,       &
                 error)
do j = 1, n2
  do i = 1, n1
    f(i,j) = real(self%f1d((j-1)*n1+i))
  end do
end do


end subroutine sll_s_nufft_2d_interpolate_array_from_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module sll_m_nufft_interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
