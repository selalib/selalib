c
c     file tmud24.f
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
c ... author and specialist
c
c          John C. Adams (National Center for Atmospheric Research)
c          email: johnad@ucar.edu, phone: 303-497-1213

c ... For MUDPACK 5.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c
c ... purpose
c
c     test program for the mudpack solver mud24
c
c ... required files
c
c     mud24.f, mud2.f, mudcom.f
c
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for mud24
c
c **********************************************************
c **********************************************************
c
c
c     a sample program/test driver for mud24 is below. it can be
c     executed as an initial test.  output is listed for the case
c     described.
c
c     test the driver below by solving the elliptic pde
c
c          (1.+y**2)*pxx + exp(-(x+y))*(pyy-py) - (x+y)*pe = r(x,y)
c
c     on the unit square with specified boundary conditions at
c     xb = 1.0, yc = 0.0 and mixed boundary conditions
c
c          dp/dx - y*p(xa,y) = g(y)  (at x=0)
c
c     and
c
c          dp/dy + x*p(x,yd) = h(x)  (at y=1).
c
c     use line relaxation in the y direction and choose a grid as close
c     to 50 by 100 as the grid size parameters allow. use the exact
c     solution
c
c          pe(x,y) = x**5 + y**5 + 1.0
c
c     for testing.
c
c     the default multigrid options with no initial guess and two
c     cycles are used when first calling mud2 to yield a second-
c     order estimate.  then mud24 is called to produce a fourth-order
c     approximation.
c
c
c ***************************************************************
c     output (64 bit floating point arithmetic)
c ****************************************************************
c
c     mud2 test
c
c     integer input arguments
c     intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
c     ixp =  3 jyq =  3 iex =  5 jey =  6
c     nx =  49 ny =  97 iguess =  0 maxcy =  2
c     method =  2 work space estimate =   70048

c
c     multigrid option arguments
c     kcycle =  2
c     iprer =  2
c     ipost =  1
c     intpol =  3
c
c     floating point input parameters
c     xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
c     tolerance (error control) =    0.000E+00
c
c     discretization call to mud2 intl =  0
c     ierror =  0 minimum work space =   70048
c
c     approximation call to mud2
c     intl =  1 method =  2 iguess =  0
c     ierror =  0
c     maximum error  =   0.342E-03
c
c     mud24 test  ierror =  0
c     maximum error  =   0.337E-05
c
c ************************************************************
c     end of output
c ************************************************************
c
c
      program tmud24
      implicit none
      integer iixp,jjyq,iiex,jjey,nnx,nny,llwork
c
c     set grid sizes with parameter statements
c
      parameter (iixp = 3 , jjyq = 3 , iiex = 5, jjey = 6)
      parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
c
c     set minimal required work space length (see tmud2.f)
c
      parameter (llwork=70048)
      real phi(nnx,nny),rhs(nnx,nny),work(llwork)
c
c     put integer and floating point argument names in contiguous
c     storeage for labelling in vectors iprm,fprm
c
      integer iprm(16),mgopt(4)
      real fprm(6)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      common/itmud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,tolmax,relmax
      common/ftmud2/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)
      integer i,j,ierror
      real dlx,dly,x,y,cxx,cyy,cx,cy,ce,pxx,pyy,px,py,pe,errmax
c
c     declare coefficient and boundary condition input subroutines external
c
      external cof,bndc
c
c
c     set input integer arguments
c
      intl = 0
c
c     set boundary condition flags
c
      nxa = 2
      nxb = 1
      nyc = 1
      nyd = 2
c
c     set grid sizes from parameter statements
c
      ixp = iixp
      jyq = jjyq
      iex = iiex
      jey = jjey
      nx = nnx
      ny = nny
c
c     set multigrid arguments (w(2,1) cycling with fully weighted
c     residual restriction and cubic prolongation)
c
      mgopt(1) = 2
      mgopt(2) = 2
      mgopt(3) = 1
      mgopt(4) = 3
c
c     set two mg cycles for second-order approximation
c
      maxcy = 2
c
c     set no initial guess forcing full multigrid cycling
c     on initial call
c
      iguess = 0
c
c     set work space length approximation from parameter statement
c
      nwork = llwork
c
c     set line-y relaxation
c
      method = 2
c
c     set end points of solution rectangle in (x,y) space
c
      xa = 0.0
      xb = 1.0
      yc = 0.0
      yd = 1.0
c
c     set mesh increments
c
      dlx = (xb-xa)/float(nx-1)
      dly = (yd-yc)/float(ny-1)
c
c     set for no error control flag
c
      tolmax = 0.0
c
c     set right hand side in rhs
c     initialize phi to zero
c
      do i=1,nx
	x = xa+float(i-1)*dlx
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call cof(x,y,cxx,cyy,cx,cy,ce)
	  call exact(x,y,pxx,pyy,px,py,pe)
	  rhs(i,j) = cxx*pxx+cyy*pyy+cx*px+cy*py+ce*pe
	  phi(i,j) = 0.0
	end do
      end do
c
c     set specified boundaries in phi
c
      x = xb
      do j=1,ny
	y = yc+float(j-1)*dly
	call exact(x,y,pxx,pyy,px,py,pe)
	phi(nx,j) = pe
      end do
      y = yc
      do i=1,nx
	x = xa+float(i-1)*dlx
	call exact(x,y,pxx,pyy,px,py,pe)
	phi(i,1) = pe
      end do
      write(*,100)
  100 format(//' mud2 test ')
      write (*,101) (iprm(i),i=1,15)
  101 format(/' integer input arguments ',
     +/' intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,
     +' nyd = ',i2,
     +/' ixp = ',i2,' jyq = ',i2,' iex = ',i2,' jey = ',i2
     +/' nx = ',i3,' ny = ',i3,' iguess = ',i2,' maxcy = ',i2,
     +/' method = ',i2, ' work space estimate = ',i7)
      write (*,102) (mgopt(i),i=1,4)
  102 format(/' multigrid option arguments ',
     +/' kcycle = ',i2,
     +/' iprer = ',i2,
     +/' ipost = ',i2,
     +/' intpol = ',i2)
      write(*,103) xa,xb,yc,yd,tolmax
  103 format(/' floating point input parameters ',
     +/' xa = ',f6.3,' xb = ',f6.3,' yc = ',f6.3,' yd = ',f6.3,
     +/' tolerance (error control) =   ',e10.3)
c
c     intiialization call
c
      write(*,104) intl
  104 format(/' discretization call to mud2', ' intl = ', i2)
      call mud2(iprm,fprm,work,cof,bndc,rhs,phi,mgopt,ierror)
c
c     print error parameter and minimum work space requirement
c
      write (*,200) ierror,iprm(16)
  200 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     attempt solution
c
      intl = 1
      write(*,106) intl,method,iguess
  106 format(/' approximation call to mud2',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
      call mud2(iprm,fprm,work,cof,bndc,rhs,phi,mgopt,ierror)
      write (*,107) ierror
  107 format(' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
c
c     compute and print maximum norm of error
c
      errmax = 0.0
      do j=1,ny
	y = yc+(j-1)*dly
	do i=1,nx
	  x = xa+(i-1)*dlx
	  call exact(x,y,pxx,pyy,px,py,pe)
	  errmax = amax1(errmax,abs((phi(i,j)-pe)))
	end do
      end do
      write(*,201) errmax
  201 format(' maximum error  =  ',e10.3)
c
c      attempt to improve approximation to fourth order
c
      call mud24(work,phi,ierror)
      write (*,108) ierror
  108 format(/' mud24 test ', ' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
c
c     compute and print maximum norm of error
c
      errmax = 0.0
      do j=1,ny
	y = yc+(j-1)*dly
	do i=1,nx
	  x = xa+(i-1)*dlx
	  call exact(x,y,pxx,pyy,px,py,pe)
	  errmax = amax1(errmax,abs((phi(i,j)-pe)))
	end do
      end do
      write(*,201) errmax
      print*,"PASSED"
      end

      subroutine cof(x,y,cxx,cyy,cx,cy,ce)
c
c     input pde coefficients at any grid point (x,y) in the solution region
c
      implicit none
      real x,y,cxx,cyy,cx,cy,ce
      cxx = 1.+y*y
      cyy = exp(-(x+y))
      cx = 0.
      cy = -cyy
      ce = -(x+y)
      return
      end

      subroutine bndc(kbdy,xory,alfa,gbdy)
c
c     input mixed derivative b.c. to mud2
c
      implicit none
      integer kbdy
      real xory,alfa,gbdy,x,y,pe,px,py,pxx,pyy
      real xa,xb,yc,yd,tolmax,relmax
      common/ftmud2/xa,xb,yc,yd,tolmax,relmax
      if (kbdy.eq.1) then  ! x=xa boundary
	y = xory
	x = xa
	call exact(x,y,pxx,pyy,px,py,pe)
	alfa = -y
	gbdy = px + alfa*pe
	return
      end if
      if (kbdy.eq.4) then  ! y=yd boundary
	y = yd
	x = xory
	call exact(x,y,pxx,pyy,px,py,pe)
	alfa = x
	gbdy = py + alfa*pe
	return
      end if
      end

      subroutine exact(x,y,pxx,pyy,px,py,pe)
c
c     this subroutine is used to set an exact solution for testing mud2
c
      implicit none
      real x,y,pxx,pyy,px,py,pe
      pe = x**5+y**5+1.
      px = 5.*x**4
      py = 5.*y**4
      pxx = 20.*x**3
      pyy = 20.*y**3
      return
      end

