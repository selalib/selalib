module sll_m_nufft
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_m_fft
use sll_m_constants

type :: sll_t_nufft_2d

  type(sll_t_fft)         :: fft
  sll_comp64, allocatable :: ftmp_1d(:)
  sll_comp64, allocatable :: f1(:)
  sll_comp64, allocatable :: f2(:)
  sll_comp64, allocatable :: fcmplx(:,:)
  sll_int32,  allocatable :: ntrace(:)
  sll_int32,  allocatable :: mtrace(:)
  sll_real64, allocatable :: r1(:)
  sll_real64, allocatable :: r2(:)
  sll_real64, allocatable :: x_array(:)
  sll_real64, allocatable :: y_array(:)
  sll_real64              :: epsnufft=1.0d-9
  sll_int32               :: nc_eta1
  sll_int32               :: nc_eta2
  sll_real64              :: eta1_min, eta1_max
  sll_real64              :: eta2_min, eta2_max

end type sll_t_nufft_2d

contains

subroutine sll_s_init_nufft_2d( self,                        &
                                nc_eta1, eta1_min, eta1_max, &
                                nc_eta2, eta2_min, eta2_max)

type(sll_t_nufft_2d)   :: self
sll_int32,  intent(in) :: nc_eta1
sll_int32,  intent(in) :: nc_eta2
sll_real64, intent(in) :: eta1_min, eta1_max
sll_real64, intent(in) :: eta2_min, eta2_max
sll_int32              :: err
sll_int32              :: i, j

  self%nc_eta1  = nc_eta1      
  self%nc_eta2  = nc_eta1      

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
  SLL_ALLOCATE(self%ftmp_1d (1:nc_eta1*nc_eta2),   err)
  SLL_ALLOCATE(self%r1      (1:nc_eta1),           err)
  SLL_ALLOCATE(self%r2      (1:nc_eta2),           err)

  call sll_s_fft_init_c2c_2d(self%fft,nc_eta1,nc_eta2, &
    self%fcmplx,self%fcmplx,sll_p_fft_forward)

  do i=1,nc_eta1
    self%r1(i) = (i-1)*2.0*sll_p_pi/ real(nc_eta1,f64) - sll_p_pi
  end do
  do j=1,nc_eta2
    self%r2(j) = (j-1)*2.0*sll_p_pi/ real(nc_eta2,f64) - sll_p_pi
  end do

end subroutine sll_s_init_nufft_2d

subroutine sll_s_free_nufft_2d( self)

type(sll_t_nufft_2d) :: self

  DEALLOCATE(self%ftmp_1d )
  DEALLOCATE(self%x_array )
  DEALLOCATE(self%y_array )
  DEALLOCATE(self%mtrace  )
  DEALLOCATE(self%ntrace  )
  DEALLOCATE(self%fcmplx  )
  DEALLOCATE(self%f1      )
  DEALLOCATE(self%f2      )

end subroutine sll_s_free_nufft_2d

subroutine sll_s_exec_nufft_2d( self, fh_fsl, t )

type(sll_t_nufft_2d)                 :: self
sll_real64, intent(inout)            :: fh_fsl(:,:)
sll_real64, intent(in)               :: t

sll_real64                           :: ct, st
sll_real64                           :: x, y
sll_real64                           :: xj, yj
sll_int32                            :: i, j
                                                         !
!#########################################################

sll_int32  :: m, n1, n2,  p

n1 = self%nc_eta1
n2 = self%nc_eta1

self%fcmplx = cmplx(fh_fsl(1:n1,1:n2),0.,f64)
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

ct = cos(t)
st = sin(t)
p  = 0
do j=1,n2
  xj = st*self%r2(j)
  yj = ct*self%r2(j)
  do i=1,n1
    x = ct*self%r1(i)-xj + sll_p_pi
    y = st*self%r1(i)+yj + sll_p_pi
    if ( 0.0_f64 < x .and. x < 2.0*sll_p_pi .and. &
         0.0_f64 < y .and. y < 2.0*sll_p_pi ) then
      p=p+1
      self%ntrace(p)   = i
      self%mtrace(p)   = j
      self%x_array(p)  = x
      self%y_array(p)  = y
    else
      fh_fsl(i,j) = 0.0_f64
    endif
  enddo
enddo

call nufft2d2f90(p,                 &
                 self%x_array(1:p), &
                 self%y_array(1:p), &
                 self%ftmp_1d(1:p), &
                 1,                 &
                 self%epsnufft,     &
                 n1,                &
                 n2,                &
                 self%fcmplx,       &
                 error)
do i=1,p
  fh_fsl(self%ntrace(i),self%mtrace(i))=real(self%ftmp_1d(i))
enddo

end subroutine sll_s_exec_nufft_2d

end module sll_m_nufft
