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


module test_fornberg
#include "sll_working_precision.h"

  use sll_m_fornberg, only: populate_weights

  implicit none

  private
  public test_weights, test_forward_5pts, test_backward_5pts

contains

  subroutine test_weights()
    sll_real64 :: z
    sll_int32 :: m, nd
    sll_real64, allocatable :: c(:,:)
    sll_real64, parameter :: x(0:2) = [-1.0_f64, 0.0_f64, 1.0_f64]
    nd = size(x)-1
    m = 2
    z = 0.0_f64
    allocate(c(0:nd, 0:m))
    call populate_weights(z, x, nd, m, c)
    print *, c
  end subroutine

  subroutine test_backward_5pts()
    sll_real64 :: z
    sll_int32 :: m, nd
    sll_real64, allocatable :: c(:,:)
    sll_real64, parameter :: x(0:4) = [-4.0_f64, -3.0_f64, -2.0_f64, -1.0_f64, 0.0_f64]
    nd = size(x)-1
    m = 2
    z = 0.0_f64
    allocate(c(0:nd, 0:m))
    call populate_weights(z, x, nd, m, c)
    print *, c(0,:)
    print *, c(:,1)
    print *, c(:,2)
  end subroutine

  subroutine test_forward_5pts()
    sll_real64 :: z
    sll_int32 :: m, nd
    sll_real64, allocatable :: c(:,:)
    sll_real64, parameter :: x(0:4) = [0.0_f64, 1.0_f64, 2.0_f64, 3.0_f64, 4.0_f64]
    nd = size(x)-1
    m = 2
    z = 0.0_f64
    allocate(c(0:nd, 0:m))
    call populate_weights(z, x, nd, m, c)
    write(*,*) " "
    print *, c(:,1)
    write(*,*) " "
    write(*,*) [-(25.0_f64/12.0_f64),4.0_f64,-3.0_f64,4.0_f64/3.0_f64,-0.25_f64]
  end subroutine

end module test_fornberg

program main
#include "sll_working_precision.h"
#include "sll_constants.h"
use test_fornberg
use sll_m_fornberg

sll_int32, parameter :: n_pts = 65
sll_real64 :: eta_min = - sll_pi
sll_real64 :: eta_max = + sll_pi
sll_real64 :: delta_eta 
sll_real64 :: xdata(1:n_pts)
sll_real64 :: ydata(1:n_pts)
sll_int32, parameter  :: maxorder = 1
sll_real64 :: zdata(0:maxorder,n_pts)
sll_int32  :: nin = 5
sll_int32  :: i
sll_real64 :: error

!call test_weights()
!call test_backward_5pts()
!call test_forward_5pts()

delta_eta = (eta_max-eta_min) / n_pts
do i = 1, n_pts
  xdata(i) = eta_min + i * delta_eta
  ydata(i) = sin(xdata(i))
end do

error = 0.0_f64
do i = 1, n_pts

  if (i <= 3) then 
    call apply_fd(nin,maxorder, &
                  xdata(1:5),ydata(1:5), &
                  xdata(i),zdata(0:1,i))
  else if (i >= n_pts-2) then
    call apply_fd(nin,maxorder, &
                  xdata(n_pts-4:n_pts),ydata(n_pts-4:n_pts), &
                  xdata(i),zdata(0:1,i))
  else
    call apply_fd(nin,maxorder, &
                  xdata(i-2:i+2),ydata(i-2:i+2), &
                  xdata(i),zdata(0:1,i))
  end if

  error = error + (zdata(1,i)-cos(xdata(i)))**2*delta_eta
  write(20,*) xdata(i), cos(xdata(i)), zdata(1,i)

end do

print*, "L2 error = ", sqrt(error), delta_eta*delta_eta

if (sqrt(error) < 1e-5) print*, 'PASSED'

end program
