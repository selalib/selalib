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

  use fornberg, only: populate_weights

  implicit none

  private
  public test_weights

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

end module test_fornberg

program main
use test_fornberg, only: test_weights
call test_weights()
end program
