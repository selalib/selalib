! it is a test of 1d nufft
program test_nufft_1d
#include "sll_working_precision.h"
use sll_m_constants
use sll_m_fft

implicit none

sll_int32  ,parameter :: Nn=64
sll_real64            :: fh_fsl(1:Nn)
sll_real64            :: t
sll_real64            :: r(1:Nn)
sll_real64            :: v(1:Nn)
sll_real64            :: lx(1:Nn)
sll_real64            :: fh
sll_real64            :: fvr(1:Nn**2)
sll_real64            :: fnufft(1:Nn**2)
complex(8)            :: ftemp(1:Nn)
complex(8)            :: ftilde(1:Nn)
complex(8)            :: sum0
complex(8)            :: cj(1:Nn**2)
type(sll_t_fft)       :: PlnF
sll_real64            :: x
sll_real64            :: y(1:Nn**2)
sll_real64            :: L
sll_real64            :: h
sll_int32             :: n,m,p,q,i,j,ier

L=4.0_f64
h=2.0_f64*L/Nn
t=0.4_f64/0.5_f64
m=Nn/2
lx=(/ (real(n,f64), n=0,m-1), (real(n,f64), n=-m,-1 )/)*sll_p_pi/L
y=0.0_f64
do i=1,Nn
  r(i)      = -L+(i-1)*h
  fh_fsl(i) = dexp(-1.5_f64*r(i)**2)
enddo
call sll_s_fft_init_c2c_1d(PlnF,Nn,ftemp,ftilde,SLL_P_FFT_FORWARD)
ftemp = cmplx(fh_fsl,0.,f64)
call sll_s_fft_exec_c2c_1d(PlnF, ftemp, ftilde)
ftilde=ftilde/cmplx(Nn,0.,f64)
j=0
do n=1,Nn
  do m=1,Nn
    x=dcos(t)*r(n)-dsin(t)*r(m)
    q=m+(n-1)*Nn
    sum0= (0.0_f64, 0.0_f64)
    if (dabs(x)<L) then
      j=j+1
      y(j)=x
      do p=1,Nn
        sum0=sum0+cdexp(sll_p_i1*lx(p)*(x+L))*ftilde(p)
      enddo
    endif
    fvr(q)=dreal(sum0)
  enddo
enddo
call sll_s_fft_free(PlnF)
ftilde=ftilde((/ (/(i,i=Nn/2+1,Nn)/), (/(i,i=1,Nn/2)/)/) )
call nufft1d2f90(j,(y+L)*sll_p_pi/L,cj, 1,1.0d-9, Nn,ftilde,ier)
j=0
do n=1,Nn
  do m=1,Nn
    x=dcos(t)*r(n)-dsin(t)*r(m)
    q=m+(n-1)*Nn
    if (dabs(x)<L) then
      j=j+1
      fnufft(q)=dreal(cj(j))
    else
      fnufft(q)=0.0_f64
    endif
  enddo
enddo

fh = maxval(abs(fvr-fnufft))

print *,'# max error', fh

end program test_nufft_1d
