module deposit_cubic_splines
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors
use sll_m_cubic_splines
use sll_m_constants
use sll_m_fft

implicit none

contains

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

L=4.0_f64
spl_fsl1 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)
spl_fsl0 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)
spl_fsl2 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)

do n=0,Ntau-1

  E0(1:Nn)=En (n,1:Nn)
  E1(1:Nn)=Enr(n,1:Nn)
  E2(1:Nn)=Ent(n,1:Nn)

  E0(Nn+1)=0.0_f64
  E1(Nn+1)=0.0_f64
  E2(Nn+1)=0.0_f64

  call sll_s_compute_cubic_spline_1D(E0,spl_fsl0)
  call sll_s_compute_cubic_spline_1D(E1,spl_fsl1)
  call sll_s_compute_cubic_spline_1D(E2,spl_fsl2)

  do j=1,Nn+1
    do i=1,Nn+1
      x=cos(tau(n))*w0(i)+sin(tau(n))*w0(j)
      if (abs(x)<L) then
        gn(n,i,j)  = sll_f_interpolate_from_interpolant_value(x,spl_fsl0)
        gnr(n,i,j) = sll_f_interpolate_from_interpolant_value(x,spl_fsl1)
        gnt(n,i,j) = sll_f_interpolate_from_interpolant_value(x,spl_fsl2)
      else
        gn(n,i,j)  = 0.0_f64
        gnr(n,i,j) = 0.0_f64
        gnt(n,i,j) = 0.0_f64
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

L=4.0_f64
spl_fsl1 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)
spl_fsl0 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)
spl_fsl2 => sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)

do n=0,Ntau-1

  E0(1:Nn)=En (n,1:Nn)
  E1(1:Nn)=Enr(n,1:Nn)
  E2(1:Nn)=Ent(n,1:Nn)

  E0(Nn+1)=0.0_f64
  E1(Nn+1)=0.0_f64
  E2(Nn+1)=0.0_f64

  call sll_s_compute_cubic_spline_1D(E0,spl_fsl0)
  call sll_s_compute_cubic_spline_1D(E1,spl_fsl1)
  call sll_s_compute_cubic_spline_1D(E2,spl_fsl2)

  do j=1,Nn+1
  do i=1,Nn+1
    x=cos(tau(n))*w1(n,i,j)+sin(tau(n))*w2(n,i,j)
    if (abs(x)<L) then
      gn( n,i,j) = sll_f_interpolate_from_interpolant_value(x,spl_fsl0)
      gnr(n,i,j) = sll_f_interpolate_from_interpolant_value(x,spl_fsl1)
      gnt(n,i,j) = sll_f_interpolate_from_interpolant_value(x,spl_fsl2)
    else
      gn( n,i,j) = 0.0_f64
      gnr(n,i,j) = 0.0_f64
      gnt(n,i,j) = 0.0_f64
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

L = 4.0_f64
spl_fsl0=>sll_f_new_cubic_spline_1D(Nn+1, -L, L, SLL_P_PERIODIC)
do n=0,Ntau-1
  E0(1:Nn)=En(n,1:Nn)
  E0(Nn+1)=0.0_f64
  call sll_s_compute_cubic_spline_1D(E0,spl_fsl0)
  do j=1,Nn+1
  do i=1,Nn+1
    x=cos(tau(n))*w1(n,i,j)+sin(tau(n))*w2(n,i,j)
    if (abs(x)<L) then
      gn(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl0)
    else
      gn(n,i,j)=0.0_f64
    endif
  enddo
  enddo
enddo
call sll_o_delete(spl_fsl0)

end subroutine ge2

end module deposit_cubic_splines
