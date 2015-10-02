c
c     file cud3ln.f
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
c ... author and specialist
c
c          John C. Adams (National Center for Atmospheric Research)
c          email: johnad@ucar.edu, phone: 303-497-1213

c ... For MUDPACK 5.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c ... purpose
c
c     cud3ln.f contains subroutines for line relaxation in the x and y
c     and z direction.  This file must be loaded with cud3.f and cuh3.f
c
      subroutine slxcd3(nx,ny,nz,phi,cof,tx,sum,nxa,nyc,nze)
c
c     x line relaxation thru red and then black points in the
c     (y,z) plane for periodic or nonperiodic x b.c.
c
      implicit none
      integer nx,ny,nz,i,ib,j,k
      integer nxa,nyc,nze,nper
      complex phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8),tx(nx,ny,nz,*)
      complex sum(ny,nz)
c
c     set periodic indicator
c
      nper =  nxa*nyc*nze
c
c     set periodic virtual boundaries as necessary
c
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
      if (nxa.ne.0) then
c
c     x direction not periodic
c     first solve for x lines thru red points in (y,z) plane
c
C$OMP PARALLEL DO PRIVATE(i,ib,j,k), SHARED(phi,cof,tx,nx,ny,nz)
	do k=1,nz,2
	  do j=1,ny,2
	    do i=1,nx
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     forward sweep
c
	  do i=2,nx
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
c
c     backward sweep
c
	  do j=1,ny,2
	    phi(nx,j,k) = phi(nx,j,k)/tx(nx,j,k,2)
	  end do
	  do ib=2,nx
	    i = nx-ib+1
	    do j=1,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k))
     +                     /tx(i,j,k,2)
	    end do
	  end do
c
c     end of k odd loop
c
	end do
C$OMP PARALLEL DO PRIVATE(i,ib,j,k), SHARED(phi,cof,tx,nx,ny,nz)
	do k=2,nz,2
	  do j=2,ny,2
	    do i=1,nx
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do j=2,ny,2
	    phi(nx,j,k) = phi(nx,j,k)/tx(nx,j,k,2)
	  end do
	  do ib=2,nx
	    i = nx-ib+1
	    do j=2,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k))
     +                     /tx(i,j,k,2)
	    end do
	  end do
c
c     end of k even loop
c
	end do

	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     solve for x lines thru black points in (y,z) plane
c
C$OMP PARALLEL DO PRIVATE(i,ib,j,k), SHARED(phi,cof,tx,nx,ny,nz)
	do k=1,nz,2
	  do j=2,ny,2
	    do i=1,nx
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do j=2,ny,2
	    phi(nx,j,k) = phi(nx,j,k)/tx(nx,j,k,2)
	  end do
	  do ib=2,nx
	    i = nx-ib+1
	    do j=2,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k))
     +                     /tx(i,j,k,2)
	    end do
	  end do
	end do
C$OMP PARALLEL DO PRIVATE(i,ib,j,k), SHARED(phi,cof,tx,nx,ny,nz)
	do k=2,nz,2
	  do j=1,ny,2
	    do i=1,nx
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i-1,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do j=1,ny,2
	    phi(nx,j,k) = phi(nx,j,k)/tx(nx,j,k,2)
	  end do
	  do ib=2,nx
	    i = nx-ib+1
	    do j=1,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k))
     +                     /tx(i,j,k,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      else
c
c     x direction periodic
c
	do k=1,nz
	  do j=1,ny
	    sum(j,k) = (0.0,0.0)
	  end do
	end do
c
c      sweep x lines thru red (y,z) forward and back
c
C$OMP PARALLEL DO PRIVATE(i,ib,j,k), SHARED(sum,phi,cof,tx,nx,ny,nz)
	do k=1,nz,2
	  do j=1,ny,2
	    do i=1,nx-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     forward sweep
c
	  do i=2,nx-2
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
	    do j=1,ny,2
	      sum(j,k) = sum(j,k)+tx(i,j,k,5)*phi(i,j,k)
	    end do
	  end do
	  do j=1,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)-sum(j,k)
	  end do
c
c     backward sweep
c
	  do j=1,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)/tx(nx-1,j,k,2)
	    phi(nx-2,j,k) = (phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))
     +                      /tx(nx-2,j,k,2)
	  end do
	  do ib=4,nx
	    i = nx-ib+1
	    do j=1,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k)-
     +                      tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-even
c
C$OMP PARALLEL DO PRIVATE(i,ib,j,k), SHARED(sum,phi,cof,tx,nx,ny,nz)
	do k=2,nz,2
	  do j=2,ny,2
	    do i=1,nx-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx-2
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
	    do j=2,ny,2
	      sum(j,k) = sum(j,k)+tx(i,j,k,5)*phi(i,j,k)
	    end do
	  end do
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)-sum(j,k)
	  end do
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)/tx(nx-1,j,k,2)
	    phi(nx-2,j,k) = (phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))/
     +                     tx(nx-2,j,k,2)
	  end do
	  do ib=4,nx
	    i = nx-ib+1
	    do j=2,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k)-
     +                      tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     now solve x lines thru black points in (y,z) plane
c
C$OMP PARALLEL DO PRIVATE(i,ib,j,k), SHARED(sum,phi,cof,tx,nx,ny,nz)
	do k=1,nz,2
	  do j=2,ny,2
	    do i=1,nx-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx-2
	    do j=2,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
	    do j=2,ny,2
	      sum(j,k) = sum(j,k)+tx(i,j,k,5)*phi(i,j,k)
	    end do
	  end do
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)-sum(j,k)
	  end do
	  do j=2,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)/tx(nx-1,j,k,2)
	    phi(nx-2,j,k) = (phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))/
     +                     tx(nx-2,j,k,2)
	  end do
	  do ib=4,nx
	    i = nx-ib+1
	    do j=2,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k)-
     +                      tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-odd
c
C$OMP PARALLEL DO PRIVATE(i,ib,j,k), SHARED(sum,phi,cof,tx,nx,ny,nz)
	do k=2,nz,2
	  do j=1,ny,2
	    do i=1,nx-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do i=2,nx-2
	    do j=1,ny,2
	      phi(i,j,k) = phi(i,j,k)-tx(i,j,k,1)*phi(i-1,j,k)
	    end do
	  end do
	  do i=1,nx-2
	    do j=1,ny,2
	      sum(j,k) = sum(j,k)+tx(i,j,k,5)*phi(i,j,k)
	    end do
	  end do
	  do j=1,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)-sum(j,k)
	  end do
	  do j=1,ny,2
	    phi(nx-1,j,k) = phi(nx-1,j,k)/tx(nx-1,j,k,2)
	    phi(nx-2,j,k) = (phi(nx-2,j,k)-tx(nx-2,j,k,4)*phi(nx-1,j,k))/
     +                     tx(nx-2,j,k,2)
	  end do
	  do ib=4,nx
	    i = nx-ib+1
	    do j=1,ny,2
	      phi(i,j,k) = (phi(i,j,k)-tx(i,j,k,3)*phi(i+1,j,k)-
     +                      tx(i,j,k,4)*phi(nx-1,j,k))/tx(i,j,k,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      end if
      end

      subroutine slycd3(nx,ny,nz,phi,cof,ty,sum,nxa,nyc,nze)
c
c     y line relaxation thru red and then black points in the
c     (x,z) plane for periodic or nonperiodic y b.c.
c
      implicit none
      integer nx,ny,nz,i,j,jb,k
      integer nxa,nyc,nze,nper
      complex phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8),ty(ny,nx,nz,*)
      complex sum(nx,nz)
c
c     set periodic indicator
c
      nper =  nxa*nyc*nze
c
c     set periodic virtual boundaries as necessary
c
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
      if (nyc.ne.0) then
c
c     y direction not periodic
c     first solve for y lines thru red points in (x,z) plane
c
C$OMP PARALLEL DO PRIVATE(i,jb,j,k), SHARED(phi,cof,ty,nx,ny,nz)
	do k=1,nz,2
	  do i=1,nx,2
	    do j=1,ny
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     forward sweep
c
	  do j=2,ny
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
c
c     backward sweep
c
	  do i=1,nx,2
	    phi(i,ny,k) = phi(i,ny,k)/ty(ny,i,k,2)
	  end do
	  do jb=2,ny
	    j = ny-jb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k))
     +                     /ty(j,i,k,2)
	    end do
	  end do
	end do
C$OMP PARALLEL DO PRIVATE(i,jb,j,k), SHARED(phi,cof,ty,nx,ny,nz)
	do k=2,nz,2
	  do i=2,nx,2
	    do j=1,ny
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,ny,k) = phi(i,ny,k)/ty(ny,i,k,2)
	  end do
	  do jb=2,ny
	    j = ny-jb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k))
     +                     /ty(j,i,k,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     solve for x lines thru black points in (y,z) plane
c
C$OMP PARALLEL DO PRIVATE(i,jb,j,k), SHARED(phi,cof,ty,nx,ny,nz)
	do k=1,nz,2
	  do i=2,nx,2
	    do j=1,ny
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,ny,k) = phi(i,ny,k)/ty(ny,i,k,2)
	  end do
	  do jb=2,ny
	    j = ny-jb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k))
     +                     /ty(j,i,k,2)
	    end do
	  end do
	end do
C$OMP PARALLEL DO PRIVATE(i,jb,j,k), SHARED(phi,cof,ty,nx,ny,nz)
	do k=2,nz,2
	  do i=1,nx,2
	    do j=1,ny
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j-1,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do i=1,nx,2
	    phi(i,ny,k) = phi(i,ny,k)/ty(ny,i,k,2)
	  end do
	  do jb=2,ny
	    j = ny-jb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k))
     +                     /ty(j,i,k,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      else
c
c     y direction periodic
c
	do k=1,nz
	  do i=1,nx
	    sum(i,k) = (0.0,0.0)
	  end do
	end do
c
c      sweep y lines thru red (x,z) forward and back
c
C$OMP PARALLEL DO PRIVATE(i,jb,j,k), SHARED(sum,phi,cof,ty,nx,ny,nz)
	do k=1,nz,2
	  do i=1,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
c
c     forward sweep
c
	  do j=2,ny-2
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
	    do i=1,nx,2
	      sum(i,k) = sum(i,k)+ty(j,i,k,5)*phi(i,j,k)
	    end do
	  end do
	  do i=1,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)-sum(i,k)
	  end do
c
c     backward sweep
c
	  do i=1,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)/ty(ny-1,i,k,2)
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))
     +                      /ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
	    j = ny-jb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-even
c
C$OMP PARALLEL DO PRIVATE(i,jb,j,k), SHARED(sum,phi,cof,ty,nx,ny,nz)
	do k=2,nz,2
	  do i=2,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny-2
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
	    do i=2,nx,2
	      sum(i,k) = sum(i,k)+ty(j,i,k,5)*phi(i,j,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)-sum(i,k)
	  end do
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)/ty(ny-1,i,k,2)
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))
     +                      /ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
	    j = ny-jb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     now solve x lines thru black points in (y,z) plane
c
C$OMP PARALLEL DO PRIVATE(i,jb,j,k), SHARED(sum,phi,cof,ty,nx,ny,nz)
	do k=1,nz,2
	  do i=2,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny-2
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
	    do i=2,nx,2
	      sum(i,k) = sum(i,k)+ty(j,i,k,5)*phi(i,j,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)-sum(i,k)
	  end do
	  do i=2,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)/ty(ny-1,i,k,2)
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))
     +                      /ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
	    j = ny-jb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     forward even-even
c
C$OMP PARALLEL DO PRIVATE(i,jb,j,k), SHARED(sum,phi,cof,ty,nx,ny,nz)
	do k=2,nz,2
	  do i=1,nx,2
	    do j=1,ny-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,5)*phi(i,j,k-1)+cof(i,j,k,6)*phi(i,j,k+1))
	    end do
	  end do
	  do j=2,ny-2
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-ty(j,i,k,1)*phi(i,j-1,k)
	    end do
	  end do
	  do j=1,ny-2
	    do i=1,nx,2
	      sum(i,k) = sum(i,k)+ty(j,i,k,5)*phi(i,j,k)
	    end do
	  end do
	  do i=1,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)-sum(i,k)
	  end do
	  do i=1,nx,2
	    phi(i,ny-1,k) = phi(i,ny-1,k)/ty(ny-1,i,k,2)
	    phi(i,ny-2,k) = (phi(i,ny-2,k)-ty(ny-2,i,k,4)*phi(i,ny-1,k))
     +                      /ty(ny-2,i,k,2)
	  end do
	  do jb=4,ny
	    j = ny-jb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-ty(j,i,k,3)*phi(i,j+1,k)-
     +                      ty(j,i,k,4)*phi(i,ny-1,k))/ty(j,i,k,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      end if
      return
      end

      subroutine slzcd3(nx,ny,nz,phi,cof,tz,sum,nxa,nyc,nze)
c
c     z line relaxation thru red and then black points in the
c     (x,y) plane for periodic or nonperiodic z b.c.
c
      implicit none
      integer nx,ny,nz,i,j,k,kb
      integer nxa,nyc,nze,nper
      complex phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8),tz(nz,nx,ny,*)
      complex sum(nx,ny)
c
c     set periodic indicator
c
      nper =  nxa*nyc*nze
c
c     set periodic virtual boundaries as necessary
c
      if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
      if (nze.ne.0) then
c
c     z direction not periodic
c     first solve for z lines thru red points in (x,y) plane
c
C$OMP PARALLEL DO PRIVATE(i,j,k,kb), SHARED(phi,cof,tz,nx,ny,nz)
	do j=1,ny,2
	  do i=1,nx,2
	    do k=1,nz
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     fo  rward sweep
c
	  do k=2,nz
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
c
c     backward sweep
c
	  do i=1,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
	    k = nz-kb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
        end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
C$OMP PARALLEL DO PRIVATE(i,j,k,kb), SHARED(phi,cof,tz,nx,ny,nz)
	do j=2,ny,2
	  do i=2,nx,2
	    do k=1,nz
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     forward
c
	  do k=2,nz
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
c
c      backward
c
	  do i=2,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
	    k = nz-kb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     solve for x lines thru black points in (y,z) plane
c
C$OMP PARALLEL DO PRIVATE(i,j,k,kb), SHARED(phi,cof,tz,nx,ny,nz)
	do j=1,ny,2
	  do i=2,nx,2
	    do k=1,nz
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     forward
c
	  do k=2,nz
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
c
c     backward sweep
c
	  do i=2,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
	    k = nz-kb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
C$OMP PARALLEL DO PRIVATE(i,j,k,kb), SHARED(phi,cof,tz,nx,ny,nz)
	do j=2,ny,2
	  do i=1,nx,2
	    do k=1,nz
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
	  do k=2,nz
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k-1,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do i=1,nx,2
	    phi(i,j,nz) = phi(i,j,nz)/tz(nz,i,j,2)
	  end do
	  do kb=2,nz
	    k = nz-kb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1))
     +                     /tz(k,i,j,2)
	    end do
	  end do
	end do
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      else
c
c     z direction periodic
c
	do j=1,ny
	  do i=1,nx
	    sum(i,j) = (0.0,0.0)
	  end do
	end do
c
c      sweep z lines thru red (x,y) forward and back
c
C$OMP PARALLEL DO PRIVATE(i,j,k,kb), SHARED(sum,phi,cof,tz,nx,ny,nz)
	do j=1,ny,2
	  do i=1,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
c
c     forward sweep
c
	  do k=2,nz-2
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
	    do i=1,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
	  end do
	  do i=1,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)-sum(i,j)
	  end do
c
c     backward sweep
c
	  do i=1,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)/tz(nz-1,i,j,2)
	    phi(i,j,nz-2) = (phi(i,j,nz-2)-tz(nz-2,i,j,4)*phi(i,j,nz-1))
     +                      /tz(nz-2,i,j,2)
	  end do
	  do kb=4,nz
	    k = nz-kb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     even-even
c
C$OMP PARALLEL DO PRIVATE(i,j,k,kb), SHARED(sum,phi,cof,tz,nx,ny,nz)
	do j=2,ny,2
	  do i=2,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
	  do k=2,nz-2
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
	    do i=2,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)-sum(i,j)
	  end do
	  do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)/tz(nz-1,i,j,2)
	    phi(i,j,nz-2) = (phi(i,j,nz-2)-tz(nz-2,i,j,4)*phi(i,j,nz-1))
     +                      /tz(nz-2,i,j,2)
	  end do
	  do kb=4,nz
	    k = nz-kb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     now solve x lines thru black points in (y,z) plane
c
C$OMP PARALLEL DO PRIVATE(i,j,k,kb), SHARED(sum,phi,cof,tz,nx,ny,nz)
	do j=1,ny,2
	  do i=2,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
	  do k=2,nz-2
	    do i=2,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
	    do i=2,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
	  end do
	  do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)-sum(i,j)
	  end do
	  do i=2,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)/tz(nz-1,i,j,2)
	    phi(i,j,nz-2) = (phi(i,j,nz-2)-tz(nz-2,i,j,4)*phi(i,j,nz-1))
     +                      /tz(nz-2,i,j,2)
	  end do
	  do kb=4,nz
	    k = nz-kb+1
	    do i=2,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nper.eq.0) call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
C$OMP PARALLEL DO PRIVATE(i,j,k,kb), SHARED(sum,phi,cof,tz,nx,ny,nz)
	do j=2,ny,2
	  do i=1,nx,2
	    do k=1,nz-1
	      phi(i,j,k) = cof(i,j,k,8) - (
     +        cof(i,j,k,1)*phi(i-1,j,k)+cof(i,j,k,2)*phi(i+1,j,k) +
     +        cof(i,j,k,3)*phi(i,j-1,k)+cof(i,j,k,4)*phi(i,j+1,k))
	    end do
	  end do
	  do k=2,nz-2
	    do i=1,nx,2
	      phi(i,j,k) = phi(i,j,k)-tz(k,i,j,1)*phi(i,j,k-1)
	    end do
	  end do
	  do k=1,nz-2
	    do i=1,nx,2
	      sum(i,j) = sum(i,j)+tz(k,i,j,5)*phi(i,j,k)
	    end do
	  end do
	  do i=1,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)-sum(i,j)
	  end do
	  do i=1,nx,2
	    phi(i,j,nz-1) = phi(i,j,nz-1)/tz(nz-1,i,j,2)
	    phi(i,j,nz-2) = (phi(i,j,nz-2)-tz(nz-2,i,j,4)*phi(i,j,nz-1))
     +                      /tz(nz-2,i,j,2)
	  end do
	  do kb=4,nz
	    k = nz-kb+1
	    do i=1,nx,2
	      phi(i,j,k) = (phi(i,j,k)-tz(k,i,j,3)*phi(i,j,k+1)-
     +                      tz(k,i,j,4)*phi(i,j,nz-1))/tz(k,i,j,2)
	    end do
	  end do
	end do
	call cper3vb(nx,ny,nz,phi,nxa,nyc,nze)
	return
      end if
      return
      end
