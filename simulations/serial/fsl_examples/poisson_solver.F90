module poisson_solver
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors
use sll_m_cubic_splines
use sll_m_constants
use sll_m_fft
use deposit_cubic_splines

implicit none

type :: poisson

  sll_comp64, allocatable              :: lx(:)
  type(sll_t_fft)                      :: fw
  type(sll_t_fft)                      :: bw
  sll_real64                           :: L
  sll_real64, allocatable              :: r(:)
  sll_int32                            :: num_cells
  type(sll_t_cubic_spline_2d), pointer :: spl_2d
  sll_real64, allocatable              :: f0(:,:)

contains

  procedure :: init   => init_poisson_solver 
  procedure :: free   => free_poisson_solver 
  procedure :: solve1 => poisson_solver_1
  procedure :: interp => poisson_interp

end type poisson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_poisson_solver( self, xmin, xmax, num_cells)

class(poisson)   :: self
sll_int32 , intent(in)  :: num_cells
sll_real64, intent(in)  :: xmin
sll_real64, intent(in)  :: xmax

sll_int32               :: i, j, m
sll_comp64, allocatable :: tmp(:)

self%num_cells = num_cells
m = num_cells/2
allocate(self%lx(num_cells))
self%lx = [ (cmplx(j,0.,f64), j=0,m-1), (cmplx(j,0.,f64), j=-m,-1 ) ]
self%lx = self%lx * 2.0d0*sll_p_pi/(xmax-xmin) * sll_p_i1 * num_cells

self%L = 4.0_f64

allocate(tmp(num_cells))
call sll_s_fft_init_c2c_1d(self%fw, num_cells, tmp, tmp, sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(self%bw, num_cells, tmp, tmp, sll_p_fft_backward)
deallocate(tmp)

self%spl_2d => sll_f_new_cubic_spline_2D(num_cells+1, &
                                         num_cells+1, &
                                                xmin, &
                                                xmax, &
                                                xmin, &
                                                xmax, &
                                      SLL_P_PERIODIC, &
                                      SLL_P_PERIODIC)

allocate(self%f0(num_cells+1,num_cells+1))
self%f0 = 0.0_f64

allocate(self%r(num_cells))
do i = 1, num_cells
  self%r(i) = xmin + (i-1) * (xmax-xmin)/real(num_cells,f64)
end do

end subroutine init_poisson_solver

subroutine free_poisson_solver( self )

class(poisson)   :: self

deallocate(self%lx)
call sll_s_fft_free(self%fw)
call sll_s_fft_free(self%bw)
call sll_o_delete(self%spl_2d)

end subroutine free_poisson_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine poisson_solver_1(self,fh_fsl,tau,num_cells,ntau,En,Enr,Ent)

class(poisson)            :: self
sll_int32,  intent(in)    :: num_cells,ntau
sll_real64, intent(in)    :: fh_fsl(num_cells,num_cells)
sll_real64, intent(in)    :: tau(0:ntau-1)
sll_real64, intent(inout) :: En(0:ntau-1,1:num_cells)
sll_real64, intent(inout) :: Enr(0:ntau-1,1:num_cells)
sll_real64, intent(inout) :: Ent(0:ntau-1,1:num_cells)

sll_real64 :: x(1:num_cells)
sll_real64 :: fvr(num_cells,num_cells)
sll_real64 :: ftv(num_cells,num_cells)
sll_real64 :: ftr(num_cells,num_cells)
sll_real64 :: ftmp1(num_cells,num_cells)
sll_real64 :: ftmp2(num_cells,num_cells)
sll_real64 :: xi1(0:ntau-1,num_cells+1,num_cells+1)
sll_real64 :: xi2(0:ntau-1,num_cells+1,num_cells+1)
sll_real64 :: v(num_cells+1)
sll_comp64 :: tmp(num_cells)
sll_comp64 :: sum0(num_cells)
sll_comp64 :: vctmp(num_cells)
sll_comp64 :: uctmp(num_cells)
sll_real64 :: gn(0:ntau-1,num_cells+1,num_cells+1)
sll_real64 :: f1, f2, f3
sll_int32  :: m, i, j

x=self%r
x(num_cells/2+1)=1.0d0
do i=0,ntau-1

  call self%interp(fh_fsl,tau(i),num_cells,fvr)
  do j=1,num_cells
    tmp = cmplx(fvr(j,:),0.,f64)
    call sll_s_fft_exec_c2c_1d(self%fw, tmp, tmp)
    sum0(j)=tmp(1)*cmplx(2.0*self%L*self%r(j)/num_cells,0.0,f64) 
  enddo
  call sll_s_fft_exec_c2c_1d(self%fw, sum0, tmp)
  
  tmp(2:num_cells)=tmp(2:num_cells)/self%lx(2:num_cells)
  tmp(1)=cmplx(0.0d0,0.0d0,kind=f64)
  call sll_s_fft_exec_c2c_1d(self%bw, tmp,tmp)

  En (i,:)=real(tmp-tmp(num_cells/2+1))/x
  Enr(i,:)=real(sum0-En(i,:))/x

enddo

do j=1,num_cells

  vctmp = cmplx(fh_fsl(j,:),0.,f64)
  uctmp = cmplx(fh_fsl(:,j),0.,f64)

  call sll_s_fft_exec_c2c_1d(self%fw, vctmp, vctmp)
  call sll_s_fft_exec_c2c_1d(self%fw, uctmp, uctmp)
  
  uctmp = uctmp/cmplx(num_cells**2,0.,f64)*self%lx
  vctmp = vctmp/cmplx(num_cells**2,0.,f64)*self%lx
 
  call sll_s_fft_exec_c2c_1d(self%bw, vctmp, tmp)
  ftv(j,:)=real(tmp)  !\partial_\x1 f_tilde(\xi1,\xi2)
  call sll_s_fft_exec_c2c_1d(self%bw, uctmp, tmp)
  ftr(:,j)=real(tmp)  !\partial_\x2 f_tilde(\xi1,\xi2)

enddo

v(1:num_cells)=self%r
v(num_cells+1)=self%L
do i=0,ntau-1
  do j=1,num_cells+1
    do m=1,num_cells+1
      xi1(i,j,m)=v(j)*cos(tau(i))-v(m)*sin(tau(i))
      xi2(i,j,m)=v(j)*sin(tau(i))+v(m)*cos(tau(i))
    enddo
  enddo
enddo

call ge2(num_cells,ntau,tau,xi1,xi2,En,gn)

do i=0,ntau-1

  call self%interp(ftv,tau(i),num_cells,ftmp1)
  call self%interp(ftr,tau(i),num_cells,ftmp2)

  do j=1,num_cells

    do m=1,num_cells
      f1 =  xi1(i,j,m)*cos(tau(i))+xi2(i,j,m)*sin(tau(i))
      f2 = -sin(tau(i))*ftmp1(j,m)+cos(tau(i))*ftmp2(j,m)
      f3 = (cos(2.0_f64*tau(i))**2 * f1 + gn(i,j,m))*f2
      vctmp(m) = cmplx(f3,0.,f64)
    enddo

    call sll_s_fft_exec_c2c_1d(self%fw, vctmp, tmp)
    sum0(j)=tmp(1)*cmplx(2.0d0*self%L*self%r(j)/num_cells,0.,f64)

  enddo

  call sll_s_fft_exec_c2c_1d(self%fw, sum0, tmp)

  tmp(2:num_cells)=tmp(2:num_cells)/self%lx(2:num_cells)
  tmp(1)=cmplx(0.0d0,0.0d0,kind=f64)

  call sll_s_fft_exec_c2c_1d(self%bw, tmp, tmp)

  Ent(i,:)=real(tmp-tmp(num_cells/2+1))/x

enddo

end subroutine poisson_solver_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine poisson_interp(self, fh_fsl,t,num_cells,fvr)

class(poisson)              :: self
sll_int32,  intent(in)      :: num_cells
sll_real64, intent(in)      :: fh_fsl(num_cells,num_cells)
sll_real64, intent(in)      :: t
sll_real64, intent(inout)   :: fvr(num_cells,num_cells)

sll_real64                  :: x, y
sll_int32                   :: i, j

self%f0(1:num_cells,1:num_cells) = fh_fsl
self%f0(num_cells+1,:)   = 0.0d0
self%f0(:,num_cells+1)   = 0.0d0

call sll_s_compute_cubic_spline_2d(self%f0, self%spl_2d)

do j=1,num_cells
  do i=1,num_cells
    x=cos(t)*self%r(i)-sin(t)*self%r(j)
    y=sin(t)*self%r(i)+cos(t)*self%r(j)
    if (abs(x)<self%L .and. abs(y)<self%L) then
      fvr(i,j)=sll_f_interpolate_value_2d(x,y,self%spl_2d)
    else
      fvr(i,j)=0.0d0
    endif
  enddo
enddo

end subroutine poisson_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module poisson_solver
