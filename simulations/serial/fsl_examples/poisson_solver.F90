module poisson_solver
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors
use sll_m_cubic_splines
use sll_m_constants
use sll_m_fft
use deposit_cubic_splines

implicit none

type :: poisson

  sll_real64, allocatable              :: lx(:)
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
  procedure :: solve2 => poisson_solver_2
  procedure :: interp => poisson_interp

end type poisson

contains

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
self%lx = [ (real(j,f64), j=0,m-1), (real(j,f64), j=-m,-1 ) ]
self%lx = self%lx * 2.0d0*sll_p_pi/(xmax-xmin)

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

!> PoissonSolver
subroutine poisson_solver_1(self,fh_fsl,tau,Nn,Ntau,En,Enr,Ent)

class(poisson)              :: self
sll_int32,  intent(in)      :: Nn,Ntau
sll_real64, intent(in)      :: fh_fsl(Nn,Nn),tau(0:Ntau-1)
sll_real64, intent(inout)   :: En(0:Ntau-1,1:Nn),Enr(0:Ntau-1,1:Nn),Ent(0:Ntau-1,1:Nn)
sll_real64    :: x(1:Nn),fvr(Nn,Nn),ftv(Nn,Nn),ftr(Nn,Nn)
sll_real64    :: ftemp1(Nn,Nn),ftemp2(Nn,Nn),xi1(0:Ntau-1,Nn+1,Nn+1),xi2(0:Ntau-1,Nn+1,Nn+1),v(Nn+1)
sll_int32 :: n,m,i
sll_comp64 :: fvptilde(1:Nn),fvptilde0(Nn),temp(Nn),sum0(Nn)
sll_comp64 :: vctmp(Nn),uctmp(Nn)
sll_real64 :: gn(0:Ntau-1,Nn+1,Nn+1)

x=self%r
x(Nn/2+1)=1.0d0
do i=0,Ntau-1
    call self%interp(fh_fsl,tau(i),Nn,fvr)
    do n=1,Nn
        vctmp = fvr(n,:)
        call sll_s_fft_exec_c2c_1d(self%fw, vctmp, fvptilde)
        sum0(n)=fvptilde(1)/real(Nn)*(2.0d0*self%L)*self%r(n) 
    enddo
    call sll_s_fft_exec_c2c_1d(self%fw,sum0, fvptilde)
    do n=2,Nn
        fvptilde(n)=fvptilde(n)/sll_p_i1/self%lx(n)/real(Nn)
    enddo
    fvptilde(1)=cmplx(0.0d0,0.0d0,kind=f64)
    call sll_s_fft_exec_c2c_1d(self%bw, fvptilde,temp)
    do n=1,Nn
        En(i,n)=real(temp(n)-temp(Nn/2+1))/x(n) 
        Enr(i,n)=real(sum0(n)-En(i,n))/x(n)
    enddo
enddo

do n=1,Nn
    vctmp = fh_fsl(n,:)
    uctmp = fh_fsl(:,n)
    call sll_s_fft_exec_c2c_1d(self%fw, vctmp, fvptilde)
    call sll_s_fft_exec_c2c_1d(self%fw, uctmp, fvptilde0)
    do m=1,Nn
    fvptilde0(m)=fvptilde0(m)/real(Nn)*sll_p_i1*self%lx(m)
    fvptilde(m)=fvptilde(m)/real(Nn)*sll_p_i1*self%lx(m)
    enddo
    call sll_s_fft_exec_c2c_1d(self%bw, fvptilde, temp)
    ftv(n,:)=real(temp)  !\partial_\xi1 f_filde(\xi1,\xi2)
    call sll_s_fft_exec_c2c_1d(self%bw, fvptilde0, temp)
    ftr(:,n)=real(temp)  !\partial_\xi2 f_filde(\xi1,\xi2)
enddo

v(1:Nn)=self%r
v(Nn+1)=self%L
do i=0,Ntau-1
    do n=1,Nn+1
    do m=1,Nn+1
        xi1(i,n,m)=v(n)*cos(tau(i))-v(m)*sin(tau(i))
        xi2(i,n,m)=v(n)*sin(tau(i))+v(m)*cos(tau(i))
    enddo
    enddo
enddo

call ge2(Nn,Ntau,tau,xi1,xi2,En,gn)

do i=0,Ntau-1
    call self%interp(ftv,tau(i),Nn,ftemp1)
    call self%interp(ftr,tau(i),Nn,ftemp2)
    do n=1,Nn
        do m=1,Nn
        vctmp(m)=(cos(2.0d0*tau(i))**2*(xi1(i,n,m)*cos(tau(i))+xi2(i,n,m)*sin(tau(i)))+gn(i,n,m))*(-sin(tau(i))*ftemp1(n,m)+cos(tau(i))*ftemp2(n,m))!partial_t f_tilde(xi1,xi2)
        enddo
        call sll_s_fft_exec_c2c_1d(self%fw, vctmp, fvptilde)
        sum0(n)=fvptilde(1)/real(Nn)*(2.0d0*self%L)*self%r(n)
    enddo
    call sll_s_fft_exec_c2c_1d(self%fw,sum0, fvptilde)
    do n=2,Nn
        fvptilde(n)=fvptilde(n)/sll_p_i1/self%lx(n)/real(Nn)
    enddo
    fvptilde(1)=cmplx(0.0d0,0.0d0,kind=f64)
    call sll_s_fft_exec_c2c_1d(self%bw, fvptilde,temp)
    do n=1,Nn
        Ent(i,n)=real(temp(n)-temp(Nn/2+1))/x(n)!E_tilde(tau,r)
    enddo
enddo
end subroutine poisson_solver_1

! ---PoissonSolver-------
subroutine poisson_solver_2(self,fh_fsl,tau,Nn,Ntau,En)

class(poisson)   :: self
sll_int32,       intent(in)    :: Nn,Ntau
sll_real64,      intent(in)    :: fh_fsl(Nn,Nn)
sll_real64,      intent(in)    :: tau(0:Ntau-1)
sll_real64,      intent(inout) :: En(0:Ntau-1,Nn)
sll_real64                     :: x(Nn)
sll_real64                     :: fvr(Nn,Nn)
sll_int32                      :: n,m,i
sll_comp64                     :: fvptilde(Nn)
sll_comp64                     :: temp(Nn)
sll_comp64                     :: sum0(Nn)
sll_comp64                     :: vctmp(Nn)

x=self%r
x(Nn/2+1)=1.0d0
do i=0,Ntau-1
  call self%interp(fh_fsl,tau(i),Nn,fvr)
  do n=1,Nn
    do m=1,Nn
        vctmp(m) = fvr(n,m)
    enddo
    call sll_s_fft_exec_c2c_1d(self%fw, vctmp, fvptilde)
    sum0(n)=fvptilde(1)/real(Nn)*(2.0d0*self%L)*self%r(n) !r*int_R fdv
  enddo
  call sll_s_fft_exec_c2c_1d(self%fw,sum0, fvptilde)
  do n=2,Nn
    fvptilde(n)=fvptilde(n)/sll_p_i1/self%lx(n)/real(Nn)
  enddo
  fvptilde(1)=cmplx(0.0d0,0.0d0,kind=f64)
  call sll_s_fft_exec_c2c_1d(self%bw, fvptilde,temp)
  do n=1,Nn
    En(i,n)=real(temp(n)-temp(Nn/2+1))/x(n) 
  enddo
enddo
end subroutine poisson_solver_2

subroutine poisson_interp(self, fh_fsl,t,n,fvr)
class(poisson)              :: self
sll_int32,  intent(in)      :: n
sll_real64, intent(in)      :: fh_fsl(n,n)
sll_real64, intent(in)      :: t
sll_real64, intent(inout)   :: fvr(n,n)

sll_real64                  :: x, y
sll_int32                   :: i, j

self%f0(1:n,1:n) = fh_fsl
self%f0(n+1,:)   = 0.0d0
self%f0(:,n+1)   = 0.0d0

call sll_s_compute_cubic_spline_2d(self%f0, self%spl_2d)

do j=1,n
  do i=1,n
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

end module poisson_solver
