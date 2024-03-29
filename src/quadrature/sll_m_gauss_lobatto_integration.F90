!> @ingroup quadrature
!> @brief
!> Gauss-Lobatto integration
!> @details
!> Low-level mathematical utility
!> that applies the
!> Gauss-Lobatto method to compute numeric integrals.
module sll_m_gauss_lobatto_integration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   implicit none

   public :: &
      sll_f_gauss_lobatto_derivative_matrix, &
      sll_o_gauss_lobatto_integrate_1d, &
      sll_f_gauss_lobatto_points, &
      sll_f_gauss_lobatto_points_and_weights, &
      sll_f_gauss_lobatto_weights

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   abstract interface
      !> 1d real function
      function function_1d(x)
         use sll_m_working_precision ! can't pass a header file because the
         ! preprocessor prevents double inclusion.
         ! This is very rare.
         sll_real64             :: function_1d
         sll_real64, intent(in) :: x
      end function function_1d
   end interface
#endif

!> Integrate numerically with Gauss-Lobatto formula
   interface sll_o_gauss_lobatto_integrate_1d
      module procedure gauss_lobatto_integral_1d
   end interface

contains

   !> @brief Gauss-Lobatto Quadrature.
   !> @details To integrate the function \f$ f(x) \f$
   !> (real-valued and of a single, real-valued argument x)
   !> over the interval \f$ [a,b] \f$, we use the Gauss-Lobatto formula
   !> \f[ \int_{-1}^1 f(x)dx \approx \sum_{k=1}^{n} w_k f(x_k) \f]
   !> where n represents the desired number of Gauss points.
   !>
   !> the function maps the interval \f$ [-1,1] \f$ into the
   !> arbitrary interval \f$ [a,b] \f$.
   !>
   !> To be considered is to split this function into degree-specific
   !> functions to avoid the select statement.
   !> @param f the function to be integrated
   !> @param[in] a left-bound of the definition interval of f
   !> @param[in] b right-bound of the definition interval of f
   !> @param[in] n the desired number of Gauss points
   !> @return The value of the integral
   function gauss_lobatto_integral_1d(f, a, b, n)
      sll_real64                :: gauss_lobatto_integral_1d
      procedure(function_1d)    :: f
      sll_real64, intent(in)    :: a
      sll_real64, intent(in)    :: b
      sll_int32, intent(in)    :: n
      sll_real64, dimension(n)  :: xk
      sll_real64, dimension(n)  :: wk
      sll_int32                 :: k
      sll_int32                 :: err
      sll_real64                :: alpha(0:n - 1), beta(0:n - 1)
      sll_real64                :: de(n), da(n), db(n)
      sll_real64                :: ans
      !sll_real64                :: x
      sll_real64                :: c1
      sll_real64                :: c2

      xk(:) = 0.0_f64
      wk(:) = 0.0_f64

      alpha = 0.0_f64
      do k = 0, n - 1
         beta(k) = real(k**2, f64)/real((2*k + 1)*(2*k - 1), f64)
      end do

      !for Gauss-Legendre and Gauss-Lobatto, beta(0)=int(dlambda)
      !see Algorithm xxx - ORTHPOL: A package of routines for  generating orthogonal
      !polynomials and Gauss-type quadrature rules by _Walter Gautschi_
      beta(0) = 2.0_f64

      call dlob(n - 2, alpha, beta, -1.0_f64, 1.0_f64, xk, wk, err, de, da, db)

      xk(1) = -1.0_f64
      xk(n) = 1.0_f64

      c1 = 0.5_f64*(b - a)
      c2 = 0.5_f64*(b + a)
      xk = c1*xk + c2
      wk = c1*wk

      ans = 0.0_f64
      do k = 1, n
         ans = ans + f(xk(k))*wk(k)
      end do
      gauss_lobatto_integral_1D = ans

   end function gauss_lobatto_integral_1d

   !> Returns a 2d array of size (2,n) containing gauss-lobatto
   !> points and weights in the interval [a,b].
   !> @param[in] n Number of gauss points.
   !> @param[in] a OPTIONAL Minimum value of the interval.
   !> @param[in] b OPTIONAL Maximun value of the interval.
   !> @return array containing points (1,:) and weights (2,:)
   function sll_f_gauss_lobatto_points_and_weights(n, a, b) result(wx)
      sll_int32, intent(in)     :: n
      sll_real64, intent(in), optional    :: a
      sll_real64, intent(in), optional    :: b
      sll_real64, dimension(2, n) :: wx !< wx points and weights

      wx(1, 1:n) = sll_f_gauss_lobatto_points(n, a, b)
      wx(2, 1:n) = sll_f_gauss_lobatto_weights(n, a, b)

   end function sll_f_gauss_lobatto_points_and_weights

   !> Returns a 1d array of size (n) containing gauss-lobatto
   !> points in the interval [a,b].
   !> @param[in] n Number of gauss points.
   !> @param[in] a OPTIONAL Minimum value of the interval.
   !> @param[in] b OPTIONAL Maximun value of the interval.
   !> @return xk array containing points
   function sll_f_gauss_lobatto_points(n, a, b) result(xk)
      sll_int32, intent(in)          :: n
      sll_real64, intent(in), optional :: a
      sll_real64, intent(in), optional :: b
      sll_real64, dimension(n)        :: xk
      sll_real64, dimension(n)        :: wk
      sll_real64                      :: c1, c2
      sll_int32                       :: k
      sll_int32                       :: err
      sll_real64                      :: alpha(0:n - 1), beta(0:n - 1)
      sll_real64                      :: de(n), da(n), db(n)

      xk(:) = 0.0_f64
      wk(:) = 0.0_f64

      alpha = 0.0_f64
      do k = 0, n - 1
         beta(k) = real(k**2, f64)/real((2*k + 1)*(2*k - 1), f64)
      end do

      !for Gauss-Legendre and Gauss-Lobatto, beta(0)=int(dlambda)
      !see Algorithm xxx - ORTHPOL: A package of routines for  generating orthogonal
      !polynomials and Gauss-type quadrature rules by _Walter Gautschi_
      beta(0) = 2.0_f64

      call dlob(n - 2, alpha, beta, -1._f64, 1._f64, xk, wk, err, de, da, db)
      ! The results of this call can yield values that are beyond
      ! the [-1;1] interval. Therefore we try to correct it
      ! by forcing the boundary values.
      ! FIXME : IT NEEDS TO BE CORRECTED. (TODO)
      xk(1) = -1.0_f64
      xk(n) = 1.0_f64

      if (present(a) .and. present(b)) then
         c1 = 0.5_f64*(b - a)
         c2 = 0.5_f64*(b + a)
         xk = c1*xk + c2
      end if

   end function sll_f_gauss_lobatto_points

   !> Returns a 1d array of size (n) containing gauss-lobatto
   !> weights in the interval [a,b].
   !> @param[in] n Number of gauss points.
   !> @param[in] a OPTIONAL Minimum value of the interval.
   !> @param[in] b OPTIONAL Maximun value of the interval.
   !> @return wk array containing points
   function sll_f_gauss_lobatto_weights(n, a, b) result(wk)
      sll_int32, intent(in)           :: n
      sll_real64, intent(in), optional :: a
      sll_real64, intent(in), optional :: b
      sll_real64, dimension(n)         :: xk
      sll_real64, dimension(n)         :: wk
      sll_int32                        :: k
      sll_int32                        :: err
      sll_real64                       :: alpha(0:n - 1), beta(0:n - 1)
      sll_real64                       :: de(n), da(n), db(n)
      sll_real64                       :: c1

      xk(:) = 0.0_f64
      wk(:) = 0.0_f64

      alpha = 0.0_f64
      do k = 0, n - 1
         beta(k) = real(k**2, f64)/real((2*k + 1)*(2*k - 1), f64)
      end do

      !for Gauss-Legendre and Gauss-Lobatto, beta(0)=int(dlambda)
      !see Algorithm xxx - ORTHPOL: A package of routines for  generating orthogonal
      !polynomials and Gauss-type quadrature rules by _Walter Gautschi_
      beta(0) = 2.0_f64

      call dlob(n - 2, alpha, beta, -1._f64, 1._f64, xk, wk, err, de, da, db)

      if (present(a) .and. present(b)) then
         c1 = 0.5_f64*(b - a)
         wk = c1*wk
      end if
   end function sll_f_gauss_lobatto_weights

   !> Construction of the derivative matrix for Gauss-Lobatto,
   !> The matrix must be already allocated of size \f$ n^2 \f$.
   !> \f[
   !>  der(i,j)=(Phi_i.Phi'_j)
   !> \f]
   function sll_f_gauss_lobatto_derivative_matrix(n, y) result(d)
      sll_int32, intent(in) :: n
      sll_real64, optional :: y(n)
      sll_real64 :: x(n)
      sll_int32  :: i, j, l, m
      sll_real64 :: prod
      sll_real64 :: d(n, n)

      if (present(y)) then
         x = y
      else
         x = sll_f_gauss_lobatto_points(n)
      end if

      d = 0.0_f64

      do j = 1, n
         do i = 1, n
            do l = 1, j - 1
               prod = 1.0_f64
               do m = 1, l - 1!min(j,l)-1
                  prod = prod*(x(i) - x(m))/(x(j) - x(m))
               end do
               do m = l + 1, j - 1!min(j,l)+1,max(j,l)-1
                  prod = prod*(x(i) - x(m))/(x(j) - x(m))
               end do
               do m = j + 1, n!max(j,l)+1,n
                  prod = prod*(x(i) - x(m))/(x(j) - x(m))
               end do
               prod = prod/(x(j) - x(l))
               d(i, j) = d(i, j) + prod
            end do
            do l = j + 1, n
               prod = 1.0_f64
               do m = 1, j - 1!min(j,l)-1
                  prod = prod*(x(i) - x(m))/(x(j) - x(m))
               end do
               do m = j + 1, l - 1!min(j,l)+1,max(j,l)-1
                  prod = prod*(x(i) - x(m))/(x(j) - x(m))
               end do
               do m = l + 1, n!max(j,l)+1,n
                  prod = prod*(x(i) - x(m))/(x(j) - x(m))
               end do
               prod = prod/(x(j) - x(l))
               d(i, j) = d(i, j) + prod
            end do
         end do
      end do

   end function sll_f_gauss_lobatto_derivative_matrix

!> This comes from http://dl.acm.org, Algorithme 726 : ORTHPOL, appendices and supplements
!> To use those functions, READ the documentation beside and find more information
!> about coefficients in paper *Algorithm 726 - ORTHPOL: A package of routines for
!> generating orthogonal polynomials and Gauss-type quadrature rules* by _Walter
!> Gautschi_ (here xxx is 726 in other references) formulas (1.1) to (1.3) page 2,
!> and book **Numerical Mathematics** by _Alfio Quarteroni_, _Riccardo Sacco_ and
!> _Fausto Saleri_ section 10.
!>
!> Given  n  and a measure  dlambda, this routine generates the
!> (n+2)-point Gauss-Lobatto quadrature formula
!>
!>   integral over supp(dlambda) of f(x)dlambda(x)
!>
!>      = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k))
!>
!>              + w(n+1)f(x(n+1)) + R(n;f).
!>
!> The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
!> =w(k), k=0,1,...,n,n+1. The user has to supply the recursion
!> coefficients  alpha(k), beta(k), k=0,1,...,n,n+1, for the measure
!> dlambda. The nodes and weights are computed in terms of the
!> eigenvalues and first component of the normalized eigenvectors of
!> a slightly modified Jacobi matrix of order  n+2. The routine calls
!> upon the subroutine  gauss  and the function subroutine  r1mach.
!>
!>   Input:
!>           - n - -  the number of interior points in the Gauss-Lobatto
!>                  formula; type integer
!>           - alpha,beta - arrays of dimension  n+2  to be supplied with
!>                  the recursion coefficients  alpha(k-1), beta(k-1),
!>                  k=1,2,...,n+2, of the underlying measure; the
!>                  routine does not use  alpha(n+2), beta(n+2)
!>           - aleft,right - the prescribed left and right endpoints
!>                  x(0)  and  x(n+1)  of the Gauss-Lobatto formula
!>
!>   Output:
!>           - zero - an array of dimension  n+2  containing the nodes (in
!>                  increasing order)  zero(k)=x(k), k=0,1,...,n,n+1
!>           - weight-an array of dimension  n+2  containing the weights
!>                  weight(k)=w(k), k=0,1,...,n,n+1
!>           - ierr - an error flag inherited from the routine  gauss
!>
!> The arrays  e,a,b  are needed for working space.
!>
   subroutine dlob(n, dalpha, dbeta, dleft, dright, dzero, dweigh, ierr, de, da, db)

      sll_int32 :: n, ierr, k, np1, np2
      sll_real64 :: dleft, dright, depsma, dp0l, dp0r, dp1l, dp1r, dpm1l
      sll_real64 :: dpm1r, ddet, dalpha(*), dbeta(*), dzero(*), dweigh(*), de(*), da(*), db(*)

      depsma = epsilon(1.0_f64)
      np1 = n + 1
      np2 = n + 2
      do 10 k = 1, np2
         da(k) = dalpha(k)
         db(k) = dbeta(k)
10       continue
         dp0l = 0.0_f64
         dp0r = 0.0_f64
         dp1l = 1.0_f64
         dp1r = 1.0_f64
         do 20 k = 1, np1
            dpm1l = dp0l
            dp0l = dp1l
            dpm1r = dp0r
            dp0r = dp1r
            dp1l = (dleft - da(k))*dp0l - db(k)*dpm1l
            dp1r = (dright - da(k))*dp0r - db(k)*dpm1r
20          continue
            ddet = dp1l*dp0r - dp1r*dp0l
            da(np2) = (dleft*dp1l*dp0r - dright*dp1r*dp0l)/ddet
            db(np2) = (dright - dleft)*dp1l*dp1r/ddet
            call dgauss(np2, da, db, depsma, dzero, dweigh, ierr, de)
            return
            end subroutine dlob

!>
!> Given  n  and a measure  dlambda, this routine generates the n-point
!> Gaussian quadrature formula
!>
!>     integral over supp(dlambda) of f(x)dlambda(x)
!>
!>        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
!>
!> The nodes are returned as  zero(k)=x(k) and the weights as
!> weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion
!> coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure
!> dlambda. The routine computes the nodes as eigenvalues, and the
!> weights in term of the first component of the respective normalized
!> eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
!> It uses a translation and adaptation of the algol procedure  imtql2,
!> Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
!> by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
!> Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
!> routine  imtql2.
!>
!>        Input:
!>              - n - - the number of points in the Gaussian quadrature
!>                      formula; type integer
!>              - alpha,beta - - arrays of dimension  n  to be filled
!>                      with the values of  alpha(k-1), beta(k-1), k=1,2,
!>                      ...,n
!>              - eps - the relative accuracy desired in the nodes
!>                      and weights
!>
!>        Output:
!>              - zero- array of dimension  n  containing the Gaussian
!>                      nodes (in increasing order)  zero(k)=x(k), k=1,2,
!>                      ...,n
!>              - weight - array of dimension  n  containing the
!>                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
!>              - ierr- an error flag equal to  0  on normal return,
!>                      equal to  i  if the QR algorithm does not
!>                      converge within 30 iterations on evaluating the
!>                      i-th eigenvalue, equal to  -1  if  n  is not in
!>                      range, and equal to  -2  if one of the beta's is
!>                      negative.
!>
!> The array  e  is needed for working space.
!>
            subroutine dgauss(n, dalpha, dbeta, deps, dzero, dweigh, ierr, de)
               sll_int32 :: n, ierr, i, ii, j, k, l, m, mml
               sll_real64 :: deps
               sll_real64 :: dp, dg, dr, ds, dc, df, db
               sll_real64 :: dalpha(n), dbeta(n), dzero(n), dweigh(n), de(n)

               if (n .lt. 1) then
                  ierr = -1
                  return
               end if
               ierr = 0
               dzero(1) = dalpha(1)
               if (dbeta(1) .lt. 0.0_f64) then
                  ierr = -2
                  return
               end if
               dweigh(1) = dbeta(1)
               if (n .eq. 1) return
               dweigh(1) = 1.0_f64
               de(n) = 0.0_f64
               do 100 k = 2, n
                  dzero(k) = dalpha(k)
                  if (dbeta(k) .lt. 0.0_f64) then
                     ierr = -2
                     return
                  end if
                  de(k - 1) = sqrt(dbeta(k))
                  dweigh(k) = 0.0_f64
100               continue
                  do 240 l = 1, n
                     j = 0
105                  do 110 m = l, n
                        if (m .eq. n) goto 120
                        if (abs(de(m)) .le. deps*(abs(dzero(m)) + abs(dzero(m + 1)))) goto 120
110                     continue
120                     dp = dzero(l)
                        if (m .eq. l) goto 240
                        if (j .eq. 30) goto 400
                        j = j + 1
                        dg = (dzero(l + 1) - dp)/(2.0_f64*de(l))
                        dr = sqrt(dg*dg + 1.0_f64)
                        dg = dzero(m) - dp + de(l)/(dg + sign(dr, dg))
                        ds = 1.0_f64
                        dc = 1.0_f64
                        dp = 0.0_f64
                        mml = m - l
                        do 200 ii = 1, mml
                           i = m - ii
                           df = ds*de(i)
                           db = dc*de(i)
                           if (abs(df) .lt. abs(dg)) goto 150
                           dc = dg/df
                           dr = sqrt(dc*dc + 1.0_f64)
                           de(i + 1) = df*dr
                           ds = 1.0_f64/dr
                           dc = dc*ds
                           goto 160
150                        ds = df/dg
                           dr = sqrt(ds*ds + 1.0_f64)
                           de(i + 1) = dg*dr
                           dc = 1.0_f64/dr
                           ds = ds*dc
160                        dg = dzero(i + 1) - dp
                           dr = (dzero(i) - dg)*ds + 2.0_f64*dc*db
                           dp = ds*dr
                           dzero(i + 1) = dg + dp
                           dg = dc*dr - db
                           df = dweigh(i + 1)
                           dweigh(i + 1) = ds*dweigh(i) + dc*df
                           dweigh(i) = dc*dweigh(i) - ds*df
200                        continue
                           dzero(l) = dzero(l) - dp
                           de(l) = dg
                           de(m) = 0.0_f64
                           goto 105
240                        continue
                           do 300 ii = 2, n
                              i = ii - 1
                              k = i
                              dp = dzero(i)
                              do 260 j = ii, n
                                 if (dzero(j) .ge. dp) goto 260
                                 k = j
                                 dp = dzero(j)
260                              continue
                                 if (k .eq. i) goto 300
                                 dzero(k) = dzero(i)
                                 dzero(i) = dp
                                 dp = dweigh(i)
                                 dweigh(i) = dweigh(k)
                                 dweigh(k) = dp
300                              continue
                                 do 310 k = 1, n
                                    dweigh(k) = dbeta(1)*dweigh(k)*dweigh(k)
310                                 continue
                                    return
400                                 ierr = l
                                    return
                                    end subroutine dgauss

                                    end module sll_m_gauss_lobatto_integration
