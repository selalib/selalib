c
c     file resm2sp.f
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1999 by UCAR                 .
c  .                                                             .
c  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                      MUDPACK version 5.0                    .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c     subroutine resm2sp(nx,ny,nxa,nxb,nyc,nyd,work,res)
c
c
c ... purpose
c
c
c     subroutine resm2sp computes the fine grid residual in the nx by ny array
c     res after calling mud2sp.  if
c
c          l * p = f
c
c     is the n by n (n = nx*ny) block tri-diagonal linear system resulting
c     from the pde discretization (done internally in mud2sp) and phi is the
c     approximation to p obtained by calling mud2sp, then resm2sp computes the
c     nx by ny residual array
c
c          res = f - l * phi.
c
c     one of the vector norms of res,
c
c          || res ||
c
c     can be computed as a "measure" of how well phi satisfies the
c     discretization equations.  for example, the following statements
c     will compute the location and size of the maximum residual in res
c     on cray computers:
c
c          ij = isamax(nx*ny,res,1)
c
c          jmax = (ij-1)/nx + 1
c
c          imax = ij - (jmax-1)*nx
c
c          resmax = abs(res(imax,jmax))
c
c
c *** please note:
c
c          let pe be the exact continuous solution to the elliptic pde
c          evaluated on the nx by ny discretization grid
c
c          let p be the exact solution to the linear discretization
c
c          let phi be the approximation to p generated by the mudpack solver
c
c     then discretization level error is defined by the condition
c
c          || phi - p || < || p - pe ||.
c                        =
c
c     a common measure of multigrid efficieny is that discretization level
c     error is reached in one full multigrid cycle (see references [2,9] in
c     the mudpack file "readme").  this can happen before the residual is
c     reduced to the level of roundoff error.  consequently, || res || is
c     a conservative measure of accuracy which can be wasteful if multi-
c     grid cycles are executed until it reaches the level of roundoff error.
c
c     || res || can be used to estimate the convergence rate of multigrid
c     iteration.  let r(n) be the residual and e(n) be the error after
c     executing n cycles.  they are related by the residual equation
c
c          l * e(n) = r(n).
c
c     it follows that the ratio
c
c          || r(n+1) || / || r(n) ||
c
c     estimates
c
c          || e(n+1) || / || e(n) ||
c
c     which in turn estimates the convergence rate
c
c          c = max || e(k+1) || / || e(k) ||.
c               k
c
c     notice
c                         n
c          || e(n) || <  c  || e(0) ||.
c
c
c ... assumptions (see mud2sp.d)
c
c     (1) nx,ny have the same values as iparm(10),iparm(11) (used
c         to set the fine grid resolution when calling mud2sp)
c
c     (2) nxa,nxb,nyc,nyd have the same values as iparm(2),iparm(3),
c         iparm(4),iparm(5) (boundary condition flags) used to call
c         mud2sp
c
c     (3) work is the same work space argument used in calling mud2sp.
c
c     (4) work has not changed since the last call to mud2sp.
c
c     If (1)-(4) are not true then resm2sp cannot compute the residual
c     in res.  (3),(4) assure a copy of the last computed phi is in work.
c
      subroutine resm2sp(nx,ny,nxa,nxb,nyc,nyd,work,res)
      implicit none
      integer nx,ny,nxa,nxb,nyc,nyd,irh,icx,icy
      real work(*),res(nx,ny)
c
c     set pointer for fine grid coefficients in work
c
      irh = 1 + (nx+2)*(ny+2)
      icx = irh + nx*ny
      icy = icx + 3*nx
      call rem2sp(nx,ny,nxa,nxb,nyc,nyd,work,work(irh),
     +            work(icx),work(icy),res)
      return
      end

      subroutine rem2sp(nx,ny,nxa,nxb,nyc,nyd,phi,rhs,cofx,cofy,res)
      implicit none
      integer nx,ny,nxa,nxb,nyc,nyd,i,j,ist,ifn,jst,jfn
      real phi(0:nx+1,0:ny+1),rhs(nx,ny)
      real cofx(nx,3),cofy(ny,3),res(nx,ny)
c
c     intialize residual to zero and set limits
c
      do j=1,ny
	do i=1,nx
	  res(i,j) = 0.0
	end do
      end do
c
c     set limits
c
      ist = 1
      if (nxa.eq.1) ist = 2
      ifn = nx
      if (nxb.eq.1) ifn = nx-1
      jst = 1
      if (nyc.eq.1) jst = 2
      jfn = ny
      if (nyd.eq.1) jfn = ny-1
c
c     compute residual on nonspecified grid points
c
      do j=jst,jfn
	do i=ist,ifn
	  res(i,j) =  rhs(i,j)-(
     +                cofx(i,1)*phi(i-1,j)+
     +                cofx(i,2)*phi(i+1,j)+
     +                cofy(j,1)*phi(i,j-1)+
     +                cofy(j,2)*phi(i,j+1)+
     +                (cofx(i,3)+cofy(j,3))*phi(i,j))
	end do
      end do
      return
      end
