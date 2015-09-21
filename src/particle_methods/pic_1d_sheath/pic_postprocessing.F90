
module pic_postprocessing
#include "sll_working_precision.h"
#include "sll_assert.h"
    implicit none
contains


    subroutine det_landau_damping( time, energy)
        sll_real64, dimension(:), intent(in) :: energy
        sll_real64, dimension(:), intent(in)  :: time
        sll_real64, dimension(:,:), allocatable :: mintab, maxtab
        sll_int32 :: N, nmaxima
        sll_real64 :: deltaminmax
        sll_real64 :: confidence=0.95_f64    !Can be changed

        sll_real64, dimension(3):: coefficients
        N=size(energy)
        SLL_ASSERT(size(time)==N)

        deltaminmax=0.01_f64 !* log(maxval(energy))

        call peakdet(maxtab,mintab,N,energy,deltaminmax,time)

        !Determine number of maxima
        nmaxima=size(maxtab,1)

        !Determine damping factor via linear regression
        !Here we loose the last maximum for cosmetic reasons
        coefficients=sll_linear_regression(maxtab(1:nmaxima-1,1),log(maxtab(1:nmaxima-1,2)),confidence)

!        print *, "Damping Factor: ", coefficients(2)/2 , "in [", &
!                        (coefficients(2)-coefficients(3))/2,",",(coefficients(2)+coefficients(3))/2,"] in 95%"
        write (*, "(A,F8.5,A,F8.5,A,F8.5,A,F4.2,A)") "Damping Factor (fit): ", coefficients(2)/2 , " in [", &
                        (coefficients(2)-coefficients(3))/2,",",(coefficients(2)+coefficients(3))/2,&
                        "] by ", confidence, "%"



    endsubroutine

    !< Linear regression with least squares
    !< See Hans-Otto Georgii - Stochastik, pp. 317
    !< y = coeffs(2)*x + coeffs(1) + epsilon
    function sll_linear_regression(x,y,confidence) result(results)
        sll_real64, dimension(:), intent(in) ::x
        sll_real64, dimension(:), intent(in) ::y
        sll_real64, dimension(3) :: results
        sll_real64, intent(in) :: confidence
        sll_real64 :: slope_x_min, slope_x_max, slope_x_error
        sll_int32 :: N
        sll_real64 :: meanx, meany
        sll_real64 :: slope_x,slope_y, off_set, slope_var,  off_set_var, var

        sll_real64 :: cov_xy, cov_xx,cov_yy
        N=size(x)
        SLL_ASSERT(size(y)==N)
        SLL_ASSERT(confidence>0)
        SLL_ASSERT(confidence< 1.0_f64)
        if (N>1) then

            meanx=sum(x)/N    !Sample mean x
            meany=sum(y)/N    !Sample mean y

            !cov_xy=dot_product(x,y) - sum(x)*sum(y)/N
            !cov_xx= (sum(x**2) - sum(x)**2/N)
            !cov_yy= (sum(y**2) - sum(y)**2/N)
            cov_xx= sum((x-meanx)**2)/(N-1)  !Sample Variance x
            cov_yy= sum((y-meany)**2)/(N-1)  !Sample Variance y
            cov_xy= sum((y-meany)*(x-meanx))/(N-1)  !Sample Covariance x,y


            slope_x= cov_xy/cov_xx
            slope_y= cov_xy/cov_yy

            off_set=meany - slope_x*meanx

            results(2)=slope_x
            results(1)=off_set

        else
            results(2)=0
            results(1)=0
        end if

        if (N>2) then

            !Variance on the residuals
            var=(cov_yy - cov_xy**2/cov_xx)/(N-2)

            slope_var= var/cov_xx
            off_set_var= var*( 1.0_f64/N + meanx**2/cov_xx )

            !Estimate confidence with Tchebyshev, since we have
            !no t-distribution available
            slope_x_error=sqrt( slope_var/(1.0_f64-confidence))
            results(3)=slope_x_error

            slope_x_min=slope_x-slope_x_error
            slope_x_max=slope_x-slope_x_error
        else
            !Could not estimate variance
            results(3)=0
        endif

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
