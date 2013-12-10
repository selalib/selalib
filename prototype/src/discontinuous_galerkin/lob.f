!!$ This comes from http://dl.acm.org, Algorithme 726 : ORTHPOL, appendices and supplements

!!$ To use those functions, READ the documentation beside and find more information 
!!$ about coefficients in paper *Algorithm xxx - ORTHPOL: A package of routines for 
!!$ generating orthogonal polynomials and Gauss-type quadrature rules* by _Walter 
!!$ Gautschi_ (here xxx is 726 in other references) formulas (1.1) to (1.3) page 2, 
!!$ and book **Numerical Mathematics** by _Alfio Quarteroni_, _Riccardo Sacco_ and 
!!$ _Fausto Saleri_ section 10.
!!$ If you just want to use Gauss-Lobatto inside Selalib, just use what is in 
!!$ gauss-lobatto.F90 ans see selalib documentation

      subroutine lob(n,alpha,beta,aleft,right,zero,weight,ierr,e,a,b)
!!$
!!$Given  n  and a measure  dlambda, this routine generates the 
!!$(n+2)-point Gauss-Lobatto quadrature formula
!!$
!!$  integral over supp(dlambda) of f(x)dlambda(x)
!!$
!!$     = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k))  
!!$
!!$             + w(n+1)f(x(n+1)) + R(n;f).
!!$
!!$The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
!!$=w(k), k=0,1,...,n,n+1. The user has to supply the recursion
!!$coefficients  alpha(k), beta(k), k=0,1,...,n,n+1, for the measure
!!$dlambda. The nodes and weights are computed in terms of the
!!$eigenvalues and first component of the normalized eigenvectors of
!!$a slightly modified Jacobi matrix of order  n+2. The routine calls 
!!$upon the subroutine  gauss  and the function subroutine  r1mach.
!!$
!!$  Input:  n - -  the number of interior points in the Gauss-Lobatto
!!$                 formula; type integer
!!$          alpha,beta - arrays of dimension  n+2  to be supplied with
!!$                 the recursion coefficients  alpha(k-1), beta(k-1),
!!$                 k=1,2,...,n+2, of the underlying measure; the
!!$                 routine does not use  alpha(n+2), beta(n+2)
!!$          aleft,right - the prescribed left and right endpoints 
!!$                 x(0)  and  x(n+1)  of the Gauss-Lobatto formula
!!$
!!$  Output: zero - an array of dimension  n+2  containing the nodes (in 
!!$                 increasing order)  zero(k)=x(k), k=0,1,...,n,n+1
!!$          weight-an array of dimension  n+2  containing the weights 
!!$                 weight(k)=w(k), k=0,1,...,n,n+1
!!$          ierr - an error flag inherited from the routine  gauss
!!$
!!$The arrays  e,a,b  are needed for working space.
!!$
      dimension alpha(*),beta(*),zero(*),weight(*),e(*),a(*),b(*)
!!$
!!$The arrays  alpha,beta,zero,weight,e,a,b  are assumed to have
!!$dimension  n+2.
!!$
      epsma=epsilon(1.0)
!!$
!!$epsma is the machine single precision.
!!$
      np1=n+1
      np2=n+2
      do 10 k=1,np2
        a(k)=alpha(k)
        b(k)=beta(k)
   10 continue
      p0l=0.
      p0r=0.
      p1l=1.
      p1r=1.
      do 20 k=1,np1
        pm1l=p0l
        p0l=p1l
        pm1r=p0r
        p0r=p1r
        p1l=(aleft-a(k))*p0l-b(k)*pm1l
        p1r=(right-a(k))*p0r-b(k)*pm1r
   20 continue
      det=p1l*p0r-p1r*p0l
      a(np2)=(aleft*p1l*p0r-right*p1r*p0l)/det
      b(np2)=(right-aleft)*p1l*p1r/det
      call gauss(np2,a,b,epsma,zero,weight,ierr,e)
      return
      end



      subroutine dlob(n,dalpha,dbeta,dleft,dright,dzero,dweigh,
     *ierr,de,da,db)
!!$
!!$This is a double-precision version of the routine  lob.
!!$
      double precision dleft,dright,depsma,dp0l,dp0r,dp1l,dp1r,dpm1l,
     *dpm1r,ddet,dalpha(*),dbeta(*),dzero(*),dweigh(*),de(*),da(*),
     *db(*)
!!$
!!$The arrays  dalpha,dbeta,dzero,dweigh,de,da,db  are assumed to have
!!$dimension  n+2.
!!$
      depsma=epsilon(1.0d0)
!!$
!!$depsma is the machine double precision.
!!$
      np1=n+1
      np2=n+2
      do 10 k=1,np2
        da(k)=dalpha(k)
        db(k)=dbeta(k)
   10 continue
      dp0l=0.d0
      dp0r=0.d0
      dp1l=1.d0
      dp1r=1.d0
      do 20 k=1,np1
        dpm1l=dp0l
        dp0l=dp1l
        dpm1r=dp0r
        dp0r=dp1r
        dp1l=(dleft-da(k))*dp0l-db(k)*dpm1l
        dp1r=(dright-da(k))*dp0r-db(k)*dpm1r
   20 continue
      ddet=dp1l*dp0r-dp1r*dp0l
      da(np2)=(dleft*dp1l*dp0r-dright*dp1r*dp0l)/ddet
      db(np2)=(dright-dleft)*dp1l*dp1r/ddet
      call dgauss(np2,da,db,depsma,dzero,dweigh,ierr,de)
      return
      end
