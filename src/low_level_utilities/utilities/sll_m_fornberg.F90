!Copyright (c) 2012-2013, Björn Ingvar Dahlgren
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions
!are met:
!
!Redistributions of source code must retain the above copyright
!notice, this list of conditions and the following disclaimer.
!
!Redistributions in binary form must reproduce the above copyright
!notice, this list of conditions and the following disclaimer in the
!documentation and/or other materials provided with the distribution.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!POSSIBILITY OF SUCH DAMAGE.  
!
!AUTHORS
!Björn Dahlgren 
!Ondřej Čertík 

module sll_m_fornberg
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_s_apply_fd, &
    sll_s_populate_weights

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

  !> Apply finite difference formula to compute derivative
  !> @param[in] nin number of points where the function is evaluated
  !> @param[in] maxorder maximum order of the derivative
  !> @param[in] xdata abscissae vector
  !> @param[in] ydata ordinates vector
  !> @param[in] xtgt position where the derivatives will be evaluated
  !> @param[out] out(0:n) values of n th derivatives at xtgt
  subroutine sll_s_apply_fd(nin, maxorder, xdata, ydata, xtgt, out)
    sll_int32,  intent(in)  :: nin, maxorder
    sll_real64, intent(in)  :: xdata(0:nin-1), ydata(0:nin-1), xtgt
    sll_real64, intent(out) :: out(0:maxorder)

    sll_int32 :: j
    sll_real64 :: c(0:nin-1, 0:maxorder)

    call sll_s_populate_weights(xtgt, xdata, nin-1, maxorder, c)
    forall(j=0:maxorder) out(j) = sum(c(0:,j)*ydata)
    
  end subroutine

  !> @param[in] z location where approximations are to be accurate,
  !> @param[in] x(0:nd) grid point locations, found in x(0:n)
  !> @param[in] nd dimension of x- and c-arrays in calling
  !>  program x(0:nd) and c(0:nd,0:m), respectively,
  !> @param[in]  m highest derivative for which weights are sought,
  !> @param[out] c(0:nd,0:m) weights at grid locations x(0:n) for
  !>  derivatives of order 0:m, found in c(0:nd,0:m)
  !>
  !> @details
  !>  See:
  !>      Generation of Finite Difference Formulas on Arbitrarily
  !>          Spaced Grids, Bengt Fornberg,
  !>          Mathematics of compuation, 51, 184, 1988, 699-706
  subroutine sll_s_populate_weights (z, x, nd, m, c)

    sll_real64, intent(in)  :: z
    sll_int32,  intent(in)  :: nd, m
    sll_real64, intent(in)  :: x(0:nd)
    sll_real64, intent(out) :: c(0:nd, 0:m)

    sll_int32  :: i, j, k, mn
    sll_real64 :: c1, c2, c3, c4, c5

    c1     = 1.0_f64
    c4     = x(0)-z
    c      = 0.0_f64
    c(0,0) = 1.0_f64
    do i=1,nd
      mn = min(i, m)
      c2 = 1.0_f64
      c5 = c4
      c4 = x(i)-z
      do j=0,i-1
        c3 = x(i)-x(j)
        c2 = c2*c3
        if (j == i-1) then
          do k = mn, 1, -1
            c(i, k) = c1*(k*c(i-1, k-1) - c5*c(i-1, k))/c2
          end do
          c(i, 0) = -c1*c5*c(i-1, 0)/c2
        endif
        do k=mn,1,-1
          c(j, k) = (c4*c(j, k) - k*c(j, k-1))/c3
        end do
        c(j, 0) = c4*c(j, 0)/c3
      end do
      c1 = c2
    end do
  end subroutine sll_s_populate_weights

end module
