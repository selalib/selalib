module deposit_cubic_splines
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors
use sll_m_cubic_splines
use sll_m_constants
use sll_m_fft

implicit none

contains

!> Interpolation by cubic spline
subroutine fvrinterp(fh_fsl,t,r,Nn,fvr)
type(sll_t_cubic_spline_2D), pointer :: spl_fsl1
sll_int32, intent(in)   :: Nn
sll_real64, intent(in)      :: fh_fsl(Nn,Nn),t,r(Nn)!,lx(1:Nn)
sll_real64, intent(inout)   :: fvr(1:Nn,1:Nn)
sll_real64    :: x,y,f0(Nn+1,Nn+1),L
sll_int32 :: n,m
L=4.0d0
spl_fsl1 => sll_f_new_cubic_spline_2D(Nn+1, &
                                      Nn+1, &
                                        -L, &
                                         L, &
                                        -L, &
                                         L, &
                            SLL_P_PERIODIC, &
                            SLL_P_PERIODIC)

f0(1:Nn,1:Nn)=fh_fsl
do n=1,Nn
  f0(n,Nn+1) = 0.0d0
  f0(Nn+1,n) = 0.0d0
enddo

call sll_s_compute_cubic_spline_2D(f0,spl_fsl1)

do n=1,Nn
  do m=1,Nn
    x=dcos(t)*r(n)-dsin(t)*r(m)
    y=dsin(t)*r(n)+dcos(t)*r(m)
    if (dabs(x)<L .and. dabs(y)<L) then
        fvr(n,m)=sll_f_interpolate_value_2D(x,y,spl_fsl1)
    else
    fvr(n,m)=0.0d0
    endif
  enddo
enddo
call sll_o_delete(spl_fsl1)
end subroutine fvrinterp

!> PoissonSolver
subroutine poissonsolver(fh_fsl,r,tau,lx,Nn,Ntau,En,Enr,Ent,PlnF,PlnB)

type(sll_t_fft), intent(in) :: PlnF,PlnB
sll_int32, intent(in)   :: Nn,Ntau
sll_real64, intent(in)      :: fh_fsl(Nn,Nn),r(Nn),lx(1:Nn),tau(0:Ntau-1)
sll_real64, intent(inout)   :: En(0:Ntau-1,1:Nn),Enr(0:Ntau-1,1:Nn),Ent(0:Ntau-1,1:Nn)
sll_real64    :: x(1:Nn),L,fvr(Nn,Nn),ftv(Nn,Nn),ftr(Nn,Nn)
sll_real64    :: ftemp1(Nn,Nn),ftemp2(Nn,Nn),xi1(0:Ntau-1,Nn+1,Nn+1),xi2(0:Ntau-1,Nn+1,Nn+1),v(Nn+1)
sll_int32 :: n,m,i
sll_comp64 :: fvptilde(1:Nn),fvptilde0(Nn),temp(Nn),sum0(Nn)
sll_comp64 :: vctmp(Nn),uctmp(Nn)
sll_real64 :: gn(0:Ntau-1,Nn+1,Nn+1)
L=4.0d0
x=r
x(Nn/2+1)=1.0d0
do i=0,Ntau-1
    call fvrinterp(fh_fsl,tau(i),r,Nn,fvr)
    do n=1,Nn
        vctmp = fvr(n,:)
        call sll_s_fft_exec_c2c_1d(PlnF, vctmp, fvptilde)
        sum0(n)=fvptilde(1)/dble(Nn)*(2.0d0*L)*r(n) !r*int_R fdv
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF,sum0, fvptilde)
    do n=2,Nn
        fvptilde(n)=fvptilde(n)/sll_p_i1/lx(n)/dble(Nn)
    enddo
    fvptilde(1)=cmplx(0.0d0,0.0d0,kind=f64)
    call sll_s_fft_exec_c2c_1d(PlnB, fvptilde,temp)
    do n=1,Nn
        En(i,n)=dreal(temp(n)-temp(Nn/2+1))/x(n) !g(tau,r)
        Enr(i,n)=dreal(sum0(n)-En(i,n))/x(n)
    enddo
enddo

do n=1,Nn
    vctmp = fh_fsl(n,:)
    uctmp = fh_fsl(:,n)
    call sll_s_fft_exec_c2c_1d(PlnF, vctmp, fvptilde)
    call sll_s_fft_exec_c2c_1d(PlnF, uctmp, fvptilde0)
    do m=1,Nn
    fvptilde0(m)=fvptilde0(m)/dble(Nn)*sll_p_i1*lx(m)
    fvptilde(m)=fvptilde(m)/dble(Nn)*sll_p_i1*lx(m)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnB, fvptilde, temp)
    ftv(n,:)=dreal(temp)  !\partial_\xi1 f_filde(\xi1,\xi2)
    call sll_s_fft_exec_c2c_1d(PlnB, fvptilde0, temp)
    ftr(:,n)=dreal(temp)  !\partial_\xi2 f_filde(\xi1,\xi2)
enddo

v(1:Nn)=r
v(Nn+1)=L
do i=0,Ntau-1
    do n=1,Nn+1
    do m=1,Nn+1
        xi1(i,n,m)=v(n)*dcos(tau(i))-v(m)*dsin(tau(i))
        xi2(i,n,m)=v(n)*dsin(tau(i))+v(m)*dcos(tau(i))
    enddo
    enddo
enddo

call ge2(Nn,Ntau,tau,xi1,xi2,En,gn)

do i=0,Ntau-1
    call fvrinterp(ftv,tau(i),r,Nn,ftemp1)
    call fvrinterp(ftr,tau(i),r,Nn,ftemp2)
    do n=1,Nn
        do m=1,Nn
        vctmp(m)=(dcos(2.0d0*tau(i))**2*(xi1(i,n,m)*dcos(tau(i))+xi2(i,n,m)*dsin(tau(i)))+gn(i,n,m))*(-dsin(tau(i))*ftemp1(n,m)+dcos(tau(i))*ftemp2(n,m))!partial_t f_tilde(xi1,xi2)
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, vctmp, fvptilde)
        sum0(n)=fvptilde(1)/dble(Nn)*(2.0d0*L)*r(n)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF,sum0, fvptilde)
    do n=2,Nn
        fvptilde(n)=fvptilde(n)/sll_p_i1/lx(n)/dble(Nn)
    enddo
    fvptilde(1)=cmplx(0.0d0,0.0d0,kind=f64)
    call sll_s_fft_exec_c2c_1d(PlnB, fvptilde,temp)
    do n=1,Nn
        Ent(i,n)=dreal(temp(n)-temp(Nn/2+1))/x(n)!E_tilde(tau,r)
    enddo
enddo
!call fftw_destroy_plan(PlnF)
!call fftw_destroy_plan(PlnB)
end subroutine poissonsolver

! ---PoissonSolver-------
subroutine poissonsolver2(fh_fsl,r,tau,lx,Nn,Ntau,En,PlnF,PlnB)

type(sll_t_fft), intent(in)    :: PlnF,PlnB
sll_int32,       intent(in)    :: Nn,Ntau
sll_real64,      intent(in)    :: fh_fsl(Nn,Nn)
sll_real64,      intent(in)    :: r(Nn)
sll_real64,      intent(in)    :: lx(Nn)
sll_real64,      intent(in)    :: tau(0:Ntau-1)
sll_real64,      intent(inout) :: En(0:Ntau-1,Nn)
sll_real64                     :: x(Nn)
sll_real64                     :: L
sll_real64                     :: fvr(Nn,Nn)
sll_int32                      :: n,m,i
sll_comp64                     :: fvptilde(Nn)
sll_comp64                     :: temp(Nn)
sll_comp64                     :: sum0(Nn)
sll_comp64                     :: vctmp(Nn)

L=4.0d0
x=r
x(Nn/2+1)=1.0d0
do i=0,Ntau-1
    call fvrinterp(fh_fsl,tau(i),r,Nn,fvr)
    do n=1,Nn
    do m=1,Nn
        vctmp(m) = fvr(n,m)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, vctmp, fvptilde)
    sum0(n)=fvptilde(1)/dble(Nn)*(2.0d0*L)*r(n) !r*int_R fdv
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF,sum0, fvptilde)
    do n=2,Nn
    fvptilde(n)=fvptilde(n)/sll_p_i1/lx(n)/dble(Nn)
    enddo
    fvptilde(1)=cmplx(0.0d0,0.0d0,kind=f64)
    call sll_s_fft_exec_c2c_1d(PlnB, fvptilde,temp)
    do n=1,Nn
    En(i,n)=dreal(temp(n)-temp(Nn/2+1))/x(n) !g(tau,r)
    enddo
enddo
!call fftw_destroy_plan(PlnF)
!call fftw_destroy_plan(PlnB)
end subroutine poissonsolver2

!----- gglobal(tau,xi1,xi2,E,lx)---
subroutine ge0(Nn,Ntau,tau,w0,En,Ent,Enr,gn,gnt,gnr)

type(sll_t_cubic_spline_1D), pointer :: spl_fsl0
type(sll_t_cubic_spline_1D), pointer :: spl_fsl1
type(sll_t_cubic_spline_1D), pointer :: spl_fsl2
sll_int32,             intent(in)    :: Nn,Ntau
sll_real64,            intent(in)    :: w0(Nn+1)
sll_real64,            intent(in)    :: tau(0:Ntau-1)
sll_real64,            intent(in)    :: En( 0:Ntau-1,Nn)
sll_real64,            intent(in)    :: Enr(0:Ntau-1,Nn)
sll_real64,            intent(in)    :: Ent(0:Ntau-1,Nn)
sll_real64,            intent(inout) :: gn( 0:Ntau-1,Nn+1,Nn+1)
sll_real64,            intent(inout) :: gnt(0:Ntau-1,Nn+1,Nn+1)
sll_real64,            intent(inout) :: gnr(0:Ntau-1,Nn+1,Nn+1)

sll_real64 :: x,L
sll_int32  :: n,i,j
sll_real64 :: E0(Nn+1)
sll_real64 :: E1(Nn+1)
sll_real64 :: E2(Nn+1)

L=4.0d0
spl_fsl1 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)
spl_fsl0 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)
spl_fsl2 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)

do n=0,Ntau-1
  do j=1,Nn
  E0(j)=En(n,j)
  E1(j)=Enr(n,j)
  E2(j)=Ent(n,j)
  enddo
  E0(Nn+1)=0.0d0
  E1(Nn+1)=0.0d0
  E2(Nn+1)=0.0d0
  call sll_s_compute_cubic_spline_1D(E0,spl_fsl0)
  call sll_s_compute_cubic_spline_1D(E1,spl_fsl1)
  call sll_s_compute_cubic_spline_1D(E2,spl_fsl2)
  do j=1,Nn+1
  do i=1,Nn+1
      x=dcos(tau(n))*w0(i)+dsin(tau(n))*w0(j)
      if (dabs(x)<L) then
      gn(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl0)
      gnr(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl1)
      gnt(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl2)
      else
      gn(n,i,j)=0.0d0
      gnr(n,i,j)=0.0d0
      gnt(n,i,j)=0.0d0
      endif
  enddo
  enddo
enddo
call sll_o_delete(spl_fsl1)
call sll_o_delete(spl_fsl0)
call sll_o_delete(spl_fsl2)
end subroutine ge0

subroutine ge1(Nn,Ntau,tau,w1,w2,En,Ent,Enr,gn,gnt,gnr)

type(sll_t_cubic_spline_1D), pointer :: spl_fsl0
type(sll_t_cubic_spline_1D), pointer :: spl_fsl1
type(sll_t_cubic_spline_1D), pointer :: spl_fsl2

sll_int32,  intent(in)    :: Nn,Ntau
sll_real64, intent(in)    :: tau(0:Ntau-1)
sll_real64, intent(in)    :: En(0:Ntau-1,Nn)
sll_real64, intent(in)    :: Enr(0:Ntau-1,Nn)
sll_real64, intent(in)    :: Ent(0:Ntau-1,Nn)
sll_real64, intent(in)    :: w1(0:Ntau-1,Nn+1,Nn+1)
sll_real64, intent(in)    :: w2(0:Ntau-1,Nn+1,Nn+1)
sll_real64, intent(inout) :: gn(0:Ntau-1,Nn+1,Nn+1)
sll_real64, intent(inout) :: gnt(0:Ntau-1,Nn+1,Nn+1)
sll_real64, intent(inout) :: gnr(0:Ntau-1,Nn+1,Nn+1)
sll_real64                :: x
sll_real64                :: L
sll_int32                 :: n
sll_int32                 :: i
sll_int32                 :: j
sll_real64                :: E0(Nn+1)
sll_real64                :: E1(Nn+1)
sll_real64                :: E2(Nn+1)

L=4.0d0
spl_fsl1 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)
spl_fsl0 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)
spl_fsl2 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)

do n=0,Ntau-1

  do j=1,Nn
    E0(j)=En(n,j)
    E1(j)=Enr(n,j)
    E2(j)=Ent(n,j)
  enddo

  E0(Nn+1)=0.0d0
  E1(Nn+1)=0.0d0
  E2(Nn+1)=0.0d0

  call sll_s_compute_cubic_spline_1D(E0,spl_fsl0)
  call sll_s_compute_cubic_spline_1D(E1,spl_fsl1)
  call sll_s_compute_cubic_spline_1D(E2,spl_fsl2)

  do j=1,Nn+1
  do i=1,Nn+1
    x=dcos(tau(n))*w1(n,i,j)+dsin(tau(n))*w2(n,i,j)
    if (dabs(x)<L) then
      gn( n,i,j) = sll_f_interpolate_from_interpolant_value(x,spl_fsl0)
      gnr(n,i,j) = sll_f_interpolate_from_interpolant_value(x,spl_fsl1)
      gnt(n,i,j) = sll_f_interpolate_from_interpolant_value(x,spl_fsl2)
    else
      gn( n,i,j) = 0.0d0
      gnr(n,i,j) = 0.0d0
      gnt(n,i,j) = 0.0d0
    endif
  enddo
  enddo

enddo

call sll_o_delete(spl_fsl1)
call sll_o_delete(spl_fsl0)
call sll_o_delete(spl_fsl2)

end subroutine ge1

subroutine ge2(Nn,Ntau,tau,w1,w2,En,gn)

type(sll_t_cubic_spline_1D), pointer :: spl_fsl0

sll_int32,  intent(in)    :: Nn
sll_int32,  intent(in)    :: Ntau
sll_real64, intent(in)    :: w1(0:Ntau-1,Nn+1,Nn+1)
sll_real64, intent(in)    :: w2(0:Ntau-1,Nn+1,Nn+1)
sll_real64, intent(in)    :: tau(0:Ntau-1)
sll_real64, intent(in)    :: En(0:Ntau-1,Nn)
sll_real64, intent(inout) :: gn(0:Ntau-1,Nn+1,Nn+1)
sll_real64                :: x
sll_real64                :: L
sll_real64                :: E0(Nn+1)
sll_int32                 :: n
sll_int32                 :: i
sll_int32                 :: j

L = 4.0d0
spl_fsl0=>sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)
do n=0,Ntau-1
  do j=1,Nn
    E0(j)=En(n,j)
  enddo
  E0(Nn+1)=0.0d0
  call sll_s_compute_cubic_spline_1D(E0,spl_fsl0)
  do j=1,Nn+1
  do i=1,Nn+1
    x=dcos(tau(n))*w1(n,i,j)+dsin(tau(n))*w2(n,i,j)
    if (dabs(x)<L) then
      gn(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl0)
    else
      gn(n,i,j)=0.0d0
    endif
  enddo
  enddo
enddo
call sll_o_delete(spl_fsl0)
end subroutine ge2

end module deposit_cubic_splines
