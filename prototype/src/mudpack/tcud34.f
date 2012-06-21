c
c     file tcud34.f
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
c     test program for the MUDPACK solver cud34
c
c ... required MUDPACK files
c
c     cud34.f, cud3.f, cudcom.f, cud3ln.f, cud3pn.f
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for cud34
c
c **********************************************************
c **********************************************************
c
c
c     a sample program/test driver for cud34 is below. it can be
c     executed as an initial test.  the output is listed for the
c     test case described.
c
c     test the driver below by solving the complex nonseparable
c     linear elliptic pde
c
c          cmplx(y,z)*d(dp/dx)/dx + cmplx(x,z)*d(dp/dy)/dy +
c
c          cmplx(x,y)*d(dp/dz)/dz - cmplx(y+z,x+y)*pe(x,y,z)
c
c          = r(x,y,z)
c
c     on the box [0.5,1.0] x [0.5,1.0] x [0.0,1.0] with Dirchlet
c     boundary conditions at the x and y boundaries and periodic
c     boundary conditions in the z direction.  use point relaxation
c     and choose a grid 33 by 33 by 97 grid.  use the exact solution
c
c          pe(x,y,z) = exp(x*y)*cmplx(cos(tpi*z),sin(tpi*z))
c
c     for testing (tpi = 8.0*atan(1.0)).  First cud3 is called with
c     iguess=0, maxcy=1 to force one full multigrid cycle as in tcud3.f.
c     Then, to ensure second order discretization level error is reached
c     (a requirement for fourth-order mudpack solvers), 3 additional
c     cycles are executed from the finest grid level (iguess=1) with
c     cud3.  Finally, cud34 is called for a fourth-order estimate.
c
c ******************************************************
c     output (64 bit floating point arithmetic)
c *******************************************************
c
c     cud3 test
c
c     input arguments
c     intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  1
c     nze =  0 nzf =  0
c     ixp =  2 jyq =  2 kzr =  3
c     iex =  5 jey =  5 kez =  6
c     nx =  33 ny =  33 nz =  97 iguess =  0 maxcy =  1
c     method =  3 work space length input = 1854498
c     xa =  0.50 xb =  1.00
c     yc =  0.50 yd =  1.00
c     ze =  0.00 zf =  1.00
c     tolmax =  0.000E+00
c
c     multigrid options
c     kcycle =  2
c     iprer =  2
c     ipost =  1
c     intpol =  3
c
c     discretization call to cud3 intl =  0
c     ierror =  0 minimum work space = 1854498
c
c     approximation call to cud3
c     intl =  1 method =  3 iguess =  0 maxcy =  1
c     ierror =  0
c     maximum error  =   0.380E-03
c
c     approximation call to cud3
c     intl =  1 method =  3 iguess =  1 maxcy =  3
c     ierror =  0
c     maximum error  =   0.381E-03
c
c     cud24 test  ierror =  0
c     maximum error  =   0.632E-06
c
c ************************************************************
c     end of output
c ************************************************************
c
      program tcud34
      implicit none
c
c     set grid sizes with parameter statements
c
      integer iixp,jjyq,kkzr,iiex,jjey,kkez,llwork,nnx,nny,nnz
      parameter(iixp=2,jjyq=2,kkzr=3)
      parameter(iiex=5,jjey=5,kkez=6)
      parameter (nnx = iixp*2**(iiex-1)+1)
      parameter (nny = jjyq*2**(jjey-1)+1)
      parameter (nnz = kkzr*2**(kkez-1)+1)
c
c     set minimal required work space (see tcud3.f) for this grid
c
      parameter (llwork = 1854498 )
c
c
c     dimension solution,right hand side, and work arrays
c
      complex phi(nnx,nny,nnz),rhs(nnx,nny,nnz),work(llwork)
      integer iprm(23),mgopt(4)
      real fprm(8)
c
c     put integer and floating point arguments names in contiguous
c     storeage labelling
c
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +              kez,nx,ny,nz,iguess,maxcy,method,meth2,nwork,
     +              lwrkqd,itero
      common/itcud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +              kez,nx,ny,nz,iguess,maxcy,method,meth2,nwork,
     +              lwrkqd,itero

      real xa,xb,yc,yd,ze,zf,tolmax,relmax,tpi
      common/ftcud3/xa,xb,yc,yd,ze,zf,tolmax,relmax,tpi
      real dlx,dly,dlz,x,y,z,errm
      complex cxx,cyy,czz,cx,cy,cz,ce, pxx,pyy,pzz,px,py,pz,pe
      integer i,j,k,ierror
      equivalence(intl,iprm)
      equivalence(xa,fprm)
c
c     declare coefficient and boundary condition input subroutines external
c
      external cof,bndc
      tpi = 8.0*atan(1.0)
c
c     set for initial call
c
      intl = 0
c
c     set boundary condition flags
c
      nxa = 1
      nxb = 1
      nyc = 1
      nyd = 1
      nze = 0
      nzf = 0
c
c     set grid sizes from parameter statements
c
      ixp = iixp
      jyq = jjyq
      kzr = kkzr
      iex = iiex
      jey = jjey
      kez = kkez
      nx = nnx
      ny = nny
      nz = nnz
c
c     set for one multigrid cycles
c
      maxcy = 1
c
c     set work space length approximation from parameter statement
c
      nwork = llwork
c
c     set method of relaxation--line in the z direction
c
      method = 3
c
c     meth2 only used in planar relaxation--but set
c
      meth2 = 0
c
c     set full multigrid cycling by flagging no initial guess at the finest
c     grid level
c
      iguess = 0
c
c     set end points of solution region in (x,y,z) space
c
      xa = 0.5
      xb = 1.0
      yc = 0.5
      yd = 1.0
      ze = 0.0
      zf = 1.0
c
c     set default multigrid options
c
      mgopt(1) = 2
      mgopt(2) = 2
      mgopt(3) = 1
      mgopt(4) = 3
c
c     set mesh increments
c
      dlx = (xb-xa)/float(nx-1)
      dly = (yd-yc)/float(ny-1)
      dlz = (zf-ze)/float(nz-1)
c
c     set for no error control
c
      tolmax = 0.0
c
c     set right hand side in rhs and phi to zero
c
      do k=1,nz
	z = ze+(k-1)*dlz
	do j=1,ny
	  y = yc+float(j-1)*dly
	  do i=1,nx
	    x = xa+float(i-1)*dlx
	    call cof(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
	    call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	    rhs(i,j,k) = cxx*pxx+cyy*pyy+czz*pzz+cx*px+cy*py+cz*pz+ce*pe
	    phi(i,j,k) = (0.,0.)
	  end do
	end do
      end do
c
c     set specified values at x and y boundaries in phi
c
      x = xa
      do k=1,nz
	z = ze+(k-1)*dlz
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	  phi(1,j,k) = pe
	end do
      end do
      x = xb
      do k=1,nz
	z = ze+(k-1)*dlz
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	  phi(nx,j,k) = pe
	end do
      end do
      y = yc
      do k=1,nz
	z = ze+(k-1)*dlz
	do i=1,nx
	  x = xa+float(i-1)*dlx
	  call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	  phi(i,1,k) = pe
	end do
      end do
      y = yd
      do k=1,nz
	z = ze+(k-1)*dlz
	do i=1,nx
	  x = xa+float(i-1)*dlx
	  call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	  phi(i,ny,k) = pe
	end do
      end do
      write(*,50)
   50 format(//' cud3 test ')
c
c     print input arguments
c
      write(*,100)intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +            nx,ny,nz,iguess,maxcy,method,nwork,xa,xb,yc,yd,ze,zf,
     +            tolmax,mgopt(1),mgopt(2),mgopt(3),mgopt(4)
  100 format(/' input arguments ',
     +/'intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2,
     +/' nze = ',i2, ' nzf = ',i2,
     +/' ixp = ',i2,' jyq = ',i2,' kzr = ',i2,
     +/' iex = ',i2, ' jey = ',i2, ' kez = ',i2,
     +/' nx = ',i3,' ny = ',i3,' nz = ',i3, ' iguess = ',i2,
     +/' maxcy = ',i2,
     +/' method = ',i2, ' work space length input = ',i7,
     +/' xa = ',f5.2,' xb = ',f5.2,
     +/' yc = ',f5.2,' yd = ',f5.2,
     +/' ze = ',f5.2,' zf = ',f5.2,
     +/' tolmax = ' ,e10.3
     +//' multigrid options '
     +/' kcycle = ',i2
     +/' iprer = ',i2
     +/' ipost = ',i2
     +/' intpol = ',i2 )
c
c     discretize pde
c
      write(*,104) intl
  104 format(/' discretization call to cud3', ' intl = ', i2)
      call cud3(iprm,fprm,work,cof,bndc,rhs,phi,mgopt,ierror)
      write (*,105) ierror,iprm(22)
  105 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     approximate pde
c
      intl = 1
      write(*,106) intl,method,iguess,maxcy
  106 format(/' approximation call to cud3 ',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2, ' maxcy = ',i2)
      call cud3(iprm,fprm,work,cof,bndc,rhs,phi,mgopt,ierror)
      write (*,107) ierror
  107 format(' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
c
c     compute and print maximum error
c
      call error(nx,ny,nz,phi,errm)
      write(*,108) errm
  108 format(' maximum error  =  ',e10.3)
c
c     execut three more cycles from finest grid level to ensure
c     second-order approximation is reached
c
      iguess = 1
      maxcy = 3
      write(*,106) intl,method,iguess,maxcy
      call cud3(iprm,fprm,work,cof,bndc,rhs,phi,mgopt,ierror)
      write (*,107) ierror
      if (ierror.gt.0) call exit(0)
      call error(nx,ny,nz,phi,errm)
      write(*,108) errm
c
c     attempt to compute fourth-order estimate
c
      call cud34(work,phi,ierror)
      write(*,109) ierror
  109 format (/ ' cud24 test ', ' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
      call error(nx,ny,nz,phi,errm)
      write(*,108) errm
      end

      subroutine error(nx,ny,nz,phi,errm)
c
c     compute the error in the estimate in phi
c
      implicit none
      integer nx,ny,nz
      complex phi(nx,ny,nz)
      real xa,xb,yc,yd,ze,zf,tolmax,relmax,errm,tpi
      common/ftcud3/xa,xb,yc,yd,ze,zf,tolmax,relmax,tpi
      real dlx,dly,dlz,x,y,z
      complex pxx,pyy,pzz,px,py,pz,pe
      integer i,j,k
      dlx = (xb-xa)/(nx-1)
      dly = (yd-yc)/(ny-1)
      dlz = (zf-ze)/(nz-1)
      errm = 0.0
      do k=1,nz
	z = ze+(k-1)*dlz
	do j=1,ny
	  y = yc+float(j-1)*dly
	  do i=1,nx
	    x = xa+float(i-1)*dlx
	    call exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
	    errm = amax1(errm,abs(phi(i,j,k)-pe))
	  end do
	end do
      end do
      return
      end

      subroutine cof(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
c
c     input pde coefficients at grid point (x,y,z) in the solution region
c     to subroutine cud3
c
      implicit none
      complex cxx,cyy,czz,cx,cy,cz,ce
      real x,y,z
      cxx = cmplx(y,z)
      cyy = cmplx(x,z)
      czz = cmplx(x,y)
      cx = (0.,0.)
      cy = cx
      cz = cz
      ce = -cmplx(y+z,x+z)
      return
      end

      subroutine bndc(kbdy,xory,yorz,alfa,gbdy)
c
c     dummy subroutine since no mixed derivative b.c.
c
      implicit none
      integer kbdy
      real xory,yorz
      complex alfa,gbdy
      return
      end

      subroutine exact(x,y,z,pxx,pyy,pzz,px,py,pz,pe)
c
c     this subroutine is used to set an exact solution for testing cud3
c
      implicit none
      real x,y,z
      complex pxx,pyy,pzz,px,py,pz,pe
      real xa,xb,yc,yd,ze,zf,tolmax,relmax,tpi
      common/ftcud3/xa,xb,yc,yd,ze,zf,tolmax,relmax,tpi
      pe = exp(x*y)*cmplx(cos(tpi*z),sin(tpi*z))
      px = y*pe
      py = x*pe
      pz = tpi*exp(x*y)*cmplx(-sin(tpi*z),cos(tpi*z))
      pxx = y*y*pe
      pyy = x*x*pe
      pzz = -tpi*tpi*pe
      return
      end

