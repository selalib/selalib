module sll_m_nufft_interpolation
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_m_fft
use sll_m_constants

implicit none

type :: sll_t_nufft_2d
  private
  type(sll_t_fft)         :: fft
  sll_comp64, allocatable :: f1d(:)
  sll_comp64, allocatable :: f1(:)
  sll_comp64, allocatable :: f2(:)
  sll_comp64, allocatable :: fcmplx(:,:)
  sll_int32,  allocatable :: ntrace(:)
  sll_int32,  allocatable :: mtrace(:)
  sll_real64, allocatable :: x_array(:)
  sll_real64, allocatable :: y_array(:)
  sll_real64              :: epsnufft=1.0d-9
  sll_int32               :: nc_eta1
  sll_int32               :: nc_eta2
  sll_real64              :: eta1_min, eta1_max
  sll_real64              :: eta2_min, eta2_max

end type sll_t_nufft_2d

contains

subroutine sll_s_nufft_2d_init( self,                        &
                                nc_eta1, eta1_min, eta1_max, &
                                nc_eta2, eta2_min, eta2_max)

type(sll_t_nufft_2d)   :: self
sll_int32,  intent(in) :: nc_eta1
sll_int32,  intent(in) :: nc_eta2
sll_real64, intent(in) :: eta1_min, eta1_max
sll_real64, intent(in) :: eta2_min, eta2_max
sll_int32              :: err

self%nc_eta1  = nc_eta1      
self%nc_eta2  = nc_eta2      

self%eta1_min = eta1_min
self%eta2_min = eta2_min
self%eta1_max = eta1_max
self%eta2_max = eta2_max

SLL_ALLOCATE(self%x_array (1:nc_eta1*nc_eta2),   err)
SLL_ALLOCATE(self%y_array (1:nc_eta1*nc_eta2),   err)
SLL_ALLOCATE(self%mtrace  (1:nc_eta1*nc_eta2),   err)
SLL_ALLOCATE(self%ntrace  (1:nc_eta1*nc_eta2),   err)
SLL_ALLOCATE(self%fcmplx  (1:nc_eta1,1:nc_eta2), err)
SLL_ALLOCATE(self%f1      (1:nc_eta1),           err)
SLL_ALLOCATE(self%f2      (1:nc_eta2),           err)
SLL_ALLOCATE(self%f1d     (1:nc_eta1*nc_eta2),   err)

!$OMP CRITICAL
call sll_s_fft_init_c2c_2d(self%fft,nc_eta1,nc_eta2, &
  self%fcmplx,self%fcmplx,sll_p_fft_forward)
!$OMP END CRITICAL

end subroutine sll_s_nufft_2d_init

subroutine sll_s_nufft_2d_free( self)

type(sll_t_nufft_2d) :: self

  DEALLOCATE(self%f1d )
  DEALLOCATE(self%x_array )
  DEALLOCATE(self%y_array )
  DEALLOCATE(self%mtrace  )
  DEALLOCATE(self%ntrace  )
  DEALLOCATE(self%fcmplx  )
  DEALLOCATE(self%f1      )
  DEALLOCATE(self%f2      )

end subroutine sll_s_nufft_2d_free

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

function sll_f_nufft_2d_interpolate_value( self, x, y ) result(f_out)

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

  self%ntrace(1)   = 1
  self%mtrace(1)   = 1
  self%x_array(1)  = xij * 2.0 * sll_p_pi
  self%y_array(1)  = yij * 2.0 * sll_p_pi

  call nufft2d2f90(1,               &
                   self%x_array(1), &
                   self%y_array(1), &
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

end function sll_f_nufft_2d_interpolate_value

subroutine sll_s_nufft_2d_interpolate_from_fft( self, f_in, x, y, f_out )

type(sll_t_nufft_2d)              :: self
sll_real64, intent(in)            :: f_in(:,:)
sll_real64, intent(in)            :: x(:,:)
sll_real64, intent(in)            :: y(:,:)
sll_real64, intent(out)           :: f_out(:,:)

sll_int32                         :: i, j, error
sll_int32                         :: n1, n2,  p
sll_real64                        :: xij, yij

n1 = self%nc_eta1
n2 = self%nc_eta2

p  = 0
do j=1,n2
  do i=1,n1
    xij = (x(i,j) - self%eta1_min) / (self%eta1_max-self%eta1_min)
    yij = (y(i,j) - self%eta2_min) / (self%eta2_max-self%eta2_min)
    SLL_ASSERT( 0.0_f64 <= xij .and. xij <= 1.0_f64)
    SLL_ASSERT( 0.0_f64 <= yij .and. yij <= 1.0_f64)
    p=p+1
    self%ntrace(p)   = i
    self%mtrace(p)   = j
    self%x_array(p)  = xij * 2.0 * sll_p_pi
    self%y_array(p)  = yij * 2.0 * sll_p_pi
  enddo
enddo

call nufft2d2f90(p,                 &
                 self%x_array(1:p), &
                 self%y_array(1:p), &
                 self%f1d(1:p), &
                 1,                 &
                 self%epsnufft,     &
                 n1,                &
                 n2,                &
                 self%fcmplx,       &
                 error)

do i=1,p
  f_out(self%ntrace(i),self%mtrace(i))=real(self%f1d(i))
enddo

end subroutine sll_s_nufft_2d_interpolate_from_fft

subroutine sll_s_nufft_2d_interpolate_array_values( self, f_in, x, y, f_out )
type(sll_t_nufft_2d)           :: self
sll_real64, intent(in)         :: f_in(:,:)
sll_real64, intent(in)         :: x(:,:)
sll_real64, intent(in)         :: y(:,:)
sll_real64, intent(out)        :: f_out(:,:)

f_out = f_in
call sll_s_nufft_2d_interpolate_array_values_inplace( self, f_out, x, y ) 

end subroutine sll_s_nufft_2d_interpolate_array_values

subroutine sll_s_nufft_2d_interpolate_array_values_inplace( self, f, x, y ) 

type(sll_t_nufft_2d)           :: self
sll_real64, intent(inout)      :: f(:,:)
sll_real64, intent(in)         :: x(:,:)
sll_real64, intent(in)         :: y(:,:)
sll_int32                      :: i, j, error
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
      self%ntrace(p)   = i
      self%mtrace(p)   = j
      self%x_array(p)  = xij * 2.0 * sll_p_pi
      self%y_array(p)  = yij * 2.0 * sll_p_pi
    else
      f(i,j) = 0.0_f64
    end if
  enddo
enddo

call nufft2d2f90(p,                 &
                 self%x_array(1:p), &
                 self%y_array(1:p), &
                 self%f1d(1:p),     &
                 1,                 &
                 self%epsnufft,     &
                 n1,                &
                 n2,                &
                 self%fcmplx,       &
                 error)

do i=1,p
  f(self%ntrace(i),self%mtrace(i))=real(self%f1d(i))
end do

end subroutine sll_s_nufft_2d_interpolate_array_values_inplace

end module sll_m_nufft_interpolation
