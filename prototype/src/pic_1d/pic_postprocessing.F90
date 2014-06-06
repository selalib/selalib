
module pic_postprocessing
#include "sll_working_precision.h"
#include "sll_assert.h"
    implicit none
contains


  subroutine det_landau_damping( time, energy)
    sll_real64, dimension(:), intent(in) :: energy
    sll_real64, dimension(:), intent(in)  :: time
    sll_real64, dimension(:,:), allocatable :: mintab, maxtab
    sll_int32 :: N
    sll_real64 :: deltaminmax
    sll_real64, dimension(3):: coefficients
    sll_real64 :: offset_energy
    N=size(energy)
    SLL_ASSERT(size(time)==N)

    deltaminmax=0.01_f64 !* log(maxval(energy))

    call peakdet(maxtab,mintab,N,energy,deltaminmax,time)

    !Determine damping factor via linear regression
    coefficients=sll_linear_regression(maxtab(:,1),log(maxtab(:,2)))

    print *, "Damping Factor: ", coefficients(2)/2 , coefficients(3)/2

  endsubroutine

 !< Linear regression with least squares
 !< y = coeffs(2)*x + coeffs(1) + epsilon
 function sll_linear_regression(x,y) result(results)
    sll_real64, dimension(:), intent(in) ::x
    sll_real64, dimension(:), intent(in) ::y
    sll_real64, dimension(3) :: results
    sll_int32 :: N
    sll_real64 :: meanx, meany
    sll_real64 :: slope_x,slope_y, off_set, slope_var,  off_set_var
    sll_real64 :: cov_xy, cov_xx,cov_yy
    N=size(x)
    SLL_ASSERT(size(y)==N)

    meanx=sum(x)/N
    meany=sum(x)/N

!         print *,"Time:" ,x
!         print *,"value:", y

    cov_xy=dot_product(x,y) - sum(x)*sum(y)/N
    cov_xx= (sum(x**2) - sum(x)**2/N)
    cov_yy= (sum(y**2) - sum(y)**2/N)

    slope_x= cov_xy/cov_xx
    slope_y= cov_xy/cov_yy
    !slope= (N*dot_product(x,y) - sum(x)*sum(y))/(N*sum(x**2) - sum(x)**2)
    off_set=meany -  slope_x*meanx

    !slope_var=
    results(3)=1.0_f64/slope_y
    results(2)=slope_x
    results(1)=off_set

 endfunction

  !PEAKDET Detect peaks in a vector
  !
  !        call PEAKDET(MAXTAB, MINTAB, N, V, DELTA) finds the local
  !        maxima and minima ("peaks") in the vector V of size N.
  !        MAXTAB and MINTAB consists of two columns. Column 1
  !        contains indices in V, and column 2 the found values.
  !
  !        call PEAKDET(MAXTAB, MINTAB, N, V, DELTA, X) replaces the
  !        indices in MAXTAB and MINTAB with the corresponding X-values.
  !
  !        A point is considered a maximum peak if it has the maximal
  !        value, and was preceded (to the left) by a value lower by
  !        DELTA relative to the preceding maximum.
  !
  ! Eli Billauer, 3.4.05 (http://billauer.co.il)
  ! Translated into Fortran by Brian McNoldy (http://andrew.rsmas.miami.edu/bmcnoldy)
  ! This function is released to the public domain; Any use is allowed.
  ! Modified for GNU Fortran by Jakob Ameres
  subroutine peakdet(maxtab,mintab,n,v,delta,x)

!    use, intrinsic :: ieee_arithmetic
    implicit none

    integer, intent(in)            :: n
    sll_real64, intent(in)               :: v(n), delta
    sll_real64, intent(in), optional     :: x(n)
    sll_real64, intent(out), allocatable :: maxtab(:,:), mintab(:,:)
    integer                        :: lookformax, i, j, c, d
    sll_real64                           :: a, NaN, Pinf, Minf, &
                                      mn, mx, mnpos, mxpos, this, &
                                      x2(n), maxtab_tmp(n,2), mintab_tmp(n,2)

    if (present(x) ) then
      x2=x
      if (size(v) /= size(x)) then
        print*,'Input vectors v and x must have same length'
      end if
    else
      forall(j=1:n) x2(j)=j*1.0_f64

    end if


    if (size((/ delta /)) > 1) then
      print*,'Input argument DELTA must be a scalar'
    end if

    if (delta <= 0) then
      print*,'Input argument DELTA must be positive'
    end if

    !NaN=ieee_value(a, ieee_quiet_nan)
    !Pinf=ieee_value(a, ieee_positive_inf)
    !Minf=ieee_value(a, ieee_negative_inf)
    NaN=-1
    Pinf=huge(a)
    Minf=tiny(a)

    mn=Pinf
    mx=Minf
    mnpos=NaN
    mxpos=NaN
    lookformax=0

    c=0 !maxima
    d=0 !minima
    maxtab_tmp(:,:)=NaN
    mintab_tmp(:,:)=NaN
    do i=1,n
      this = v(i)
      if (this > mx) then
        mx = this
        mxpos = x2(i)
      end if
      if (this < mn) then
        mn = this
        mnpos = x2(i)
      end if
      if (lookformax==1) then
        if (this < mx*(1.0_f64-delta)) then
          c=c+1
          maxtab_tmp(c,:)=(/ mxpos, mx /)
          mn = this
          mnpos = x2(i)
          lookformax = 0
        end if
      else
        if (this > mn*(1.0_f64+delta)) then
          d=d+1
          mintab_tmp(d,:)=(/ mnpos, mn /)
          mx = this
          mxpos = x2(i)
          lookformax = 1
        end if
      end if
    end do

    allocate(maxtab(c,2))
    allocate(mintab(d,2))
    !where (.not.isnan(maxtab_tmp))
    j=1
    do i=1,n

       if(maxtab_tmp(i,1) /=NaN) then
           maxtab(j,:)=maxtab_tmp(i,:)
           j=j+1
        endif
    enddo


   j=1
    do i=1,n

       if(mintab_tmp(i,1) /=NaN) then
           mintab(j,:)=mintab_tmp(i,:)
           j=j+1
        endif
    enddo


  end subroutine peakdet

end module pic_postprocessing
