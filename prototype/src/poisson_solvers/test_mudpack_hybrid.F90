!     test the hybrid multigrid/direct method solver by approximating
!     the solution to the nonseparable elliptic pde in divergence form on
!     the 2.5 degree grid on the full surface of a sphere of radius one.
!
!          div(sigma(t,p)*grad(u(t,p))) - lambda(t,p)*u(t,p) = f(t,p)
!
!     t and p are colatitude and longitude in radians.  expanding
!     the pde and multiplying thru by sin(t) puts it in the following
!     form suitable for muh2:
!
!          ctt * utt + cpp*upp + ct*ut + cp *up + ce * u = r(t,p)
!
!     the coefficients are given by
!
!          ctt(t,p) = sin(t)*sigma(t,p)
!
!          cpp(t,p) = (1/sin(t))*sigma(t,p)
!
!          ct(t,p)  = sin(t)*d(sigma)/dt + cos(t)*sigma(t,p)
!
!          cp(t,p)  = (1/sin(t))*d(sigma)/dp
!
!          ce(t,p)  = - sin(t)*lambda(t,p)
!
!          r(t,p)   = sin(t)* f(t,p)
!
!     for testing use the coefficients and exact solution:
!
!          sigma(t,p) = 1.5 + (sin(t)*cos(p))**2
!
!          lambda(t,p) = - sigma(t,p)
!
!          u(t,p) =  [sin(t)*(sin(t)*cos(t)*sin(p)*cos(p))]**2
!
!     (the exact solution is the restriction of the solution u(x,y,z) =
!     (x*y*z)**2 in cartesian coordinates to the surface of the sphere
!     in spherical coordinates).  assume the solution u(t,p) is specified
!     at the poles and is periodic in longitude.  choosing the grid size
!     parameters:
!
!          itp = iparm(6) = 9, jpq = iparm(7) = 9
!
!          iet = iparm(8) = 4,  jep = iparm(9) = 5
!
!     fits the required five degree 73 by 145 grid exactly.  the 10 x 10
!     coarsest grid is too large for effective error reduction with multgrid
!     iteration using relaxation only.  the subgrid sizes in the coarsening
!     muh2 and muh24 uses are:
!
!          73 X 145 > 37 X 73 > 19 X 37 > 10 X 19 > 10 X 10
!
!     Guassian elimination is used whenever the coarsest 10 X 10 subgrid
!     is encountered within multigrid iteration and line relaxation in
!     the phi direction is used at the higher resolution grids.  the
!     default multigrid options are used. one full multigrid cycle with
!     no initial guess is executed with muh2.  then, to ensure a second-
!     order approximation is reached, two more cycles are executed calling
!     muh2 with iguess = 1.  Due to fortuitous error cancellations, the
!     additional cycles actually increase the error measure with the
!     continuous solution but do give a better approximation to the
!     solution of the second-order finite difference equations.  Finally
!     muh24 is called to produce a fourth-order estimate.
!
! ************************************************************
!     output (64 bit floating point arithmetic)
! ***********************************************************
!
!     muh2 test
!
!     integer input parameters
!     intl =  0 nta =  1 ntb =  1 npc =  0 npd =  0
!     itp =  9 jpq =  9 iet =  4 jep =  5
!     nt =  73 np = 145 iguess =  0 maxcy =  1
!     method =  2 work space estimate =  214396
!
!     multigrid option parameters
!     kcycle =  2
!     iprer =  2
!     ipost =  1
!     intpol =  3
!
!     floating point input parameters
!     ta =  0.000 tb =  3.142 pc =  0.000 pd =  6.283
!     tolerance (error control) =   0.000E+00
!
!     discretization call to muh2  intl =  0
!     ierror =  0 minimum work space =  186661
!
!     approximation call to muh2
!     intl =  1 method =  2 iguess =  0 maxcy =  1
!     tolmax =  0.00
!     ierror =  0
!     maximum error  =  0.795E-04
!
!     approximation call to muh2
!     intl =  1 method =  2 iguess =  1 maxcy =  2
!     tolmax =  0.00
!     ierror =  0
!     maximum error  =  0.839E-04
!
!     muh24 test  ierror =  0
!     maximum error  =  0.209E-06
!
! ***********************************************************
!

program tmuh24
#include "sll_working_precision.h"

implicit none
sll_int32 :: iitp,jjpq,iiet,jjep,nnt,nnp,llwrk,lldir,llwork
sll_int32 :: iitp1,jjpq1

!set grid size with parameter statements
parameter (iitp = 9, jjpq =  9, iiet = 4, jjep = 5)
parameter (iitp1 = iitp+1, jjpq1 = jjpq+1)
parameter ( nnt = iitp*2**(iiet-1) + 1)
parameter ( nnp = jjpq*2**(jjep-1) + 1)

!set work space estimate (see muh2.d)
parameter (llwrk=4*(15*nnt*nnp+8*(nnt+nnp+2))/3 )
parameter(lldir=(2*(iitp+1)*(2*jjpq-1)+jjpq+1))
parameter (llwork = llwrk+lldir)

!dimension solution,right hand side, and work arrays
sll_real64 :: u(nnt,nnp),r(nnt,nnp),w(llwork)
sll_int32  :: iw(iitp1,jjpq1)

!dimension input argument vectors and set up continguous storage
!for labelling entries
sll_int32  :: iparm(17),mgopt(4)
sll_real64 :: fparm(6)
sll_int32  :: intl,nta,ntb,npc,npd,itp,jpq,iet,jep,nt,np, &
                iguess,maxcy,method,nwork,lwork,iter
sll_real64 :: ta,tb,pc,pd,tolmax,sinat,cosat,sinap,cosap,dlt,dlp
sll_int32  :: i,j,ierror
sll_real64 :: pi,sint,cost,sinp,cosp
sll_real64 :: ctt,cpp,ct,cp,ce,tmp,dt,dp,dtt,dpp,ue,ut,up,utt,upp
sll_real64 :: errm,p,t

common / iprm / intl,nta,ntb,npc,npd,itp,jpq,iet,jep,nt,np, &
                iguess,maxcy,method,nwork,lwork,iter
common / fprm / ta,tb,pc,pd,tolmax,sinat(73),cosat(73), &
                sinap(145),cosap(145),dlt,dlp
!equivlance iparm,fparm with labelled commons iprm,fprm
equivalence (intl,iparm)
equivalence (ta,fparm)

!declare coefficient and boundary condition input subroutines external
external cof,bndc

!set input integer parameters
!initialization
intl = 0

!set boundary condition flags: poles specified, longitude periodic
nta = 1
ntb = 1
npc = 0
npd = 0

!set grid sizes from parameter statements
itp = iitp
jpq = jjpq
iet = iiet
jep = jjep
nt = nnt
np = nnp

!full multigrid cycling (no initial guess at finest grid)
iguess = 0

!set one multigrid cycle
maxcy = 1

!set line relaxation in the longitude direction
method = 2

!set work space estimate
nwork = llwork

!set multigrid parameters (w(2,1) cycling with fully weighted
!residual restriction and cubic prolongation).  this can also
!be set by inputting mgopt(1) = 0 to muh2
mgopt(1) = 2
mgopt(2) = 2
mgopt(3) = 1
mgopt(4) = 3

!set floating point input parameters
pi = 4.0*atan(1.0)

!interval end points (in radians)
ta = 0.0
tb = pi
pc = 0.
pd = pi+pi

!no error control
tolmax = 0.0

!set mesh increments
dlt = (tb-ta)/float(nt-1)
dlp = (pd-pc)/float(np-1)

!preset sin,cos vectors to save computation on grid points
do i=1,nt
   t = ta+(i-1)*dlt
   sinat(i) = sin(t)
   cosat(i) = cos(t)
end do
do j=1,np
   p = pc+(j-1)*dlp
   sinap(j) = sin(p)
   cosap(j) = cos(p)
end do

!initialize right hand side and solution array except at poles
do i=2,nt-1
   t = ta+(i-1)*dlt
   sint = sinat(i)
   cost = cosat(i)
   do j=1,np
      p = pc+(j-1)*dlp
      call cof(t,p,ctt,cpp,ct,cp,ce)

      !set intermediate variables for exact solution
      sinp = sinap(j)
      cosp = cosap(j)
      tmp = (sint*cosp*sint*sinp*cost)
      dt = (2.*sint*cost*cost-sint**3)*(cosp*sinp)
      dp = (cosp**2-sinp**2)*(sint**2*cost)
      dtt = (2.*cost**3-4.*cost*sint**2-3.*sint**2*cost)*(cosp*sinp)
      dpp = (-4.*cosp*sinp)*(sint**2*cost)

      !set continuous solution and partial derivatives
      ue = tmp*tmp
      ut = 2.*tmp*dt
      up = 2.*tmp*dp
      utt = 2.*(dt*dt+tmp*dtt)
      upp = 2.*(dp*dp+tmp*dpp)

      !set right hand side of continuous pde on grid
      r(i,j) = ctt*utt+cpp*upp+ct*ut+cp*up+ce*ue

      !initialize solution array to zero
      u(i,j) = 0.0
   end do
end do

!set u, r(unused) at poles
do j=1,np
   u(1,j) = 0.0
   r(1,j) = 0.0
   u(nt,j) = 0.0
   r(nt,j) = 0.0
end do

!print input parameters
write(*,100)
write(*,101) (iparm(i),i=1,15)
write(*,102) (mgopt(i),i=1,4)
write(*,103) (fparm(i),i=1,5)
write(*,104) intl

!discretization call to muh2
call muh2(iparm,fparm,w,iw,cof,bndc,r,u,mgopt,ierror)
write (*,105) ierror,iparm(16)
if (ierror.gt.0) call exit(0)

!aprroximation call to muh2
intl = 1
write(*,106) intl, method, iguess, maxcy, tolmax
call muh2(iparm,fparm,w,iw,cof,bndc,r,u,mgopt,ierror)
write (*,107) ierror
if (ierror.gt.0) call exit(0)
if (ierror .le. 0) then
   !compute and print exact maximum error
   errm = 0.0
   do j=1,np
      sinp = sinap(j)
      cosp = cosap(j)
      do i=1,nt
         sint = sinat(i)
         cost = cosat(i)
         ue = (sint*cosp*sint*sinp*cost)**2
         errm = dmax1(errm,abs((u(i,j)-ue)))
      end do
   end do
   write(*,108) errm
end if

!execute two more cycles with iguess=1 to ensure second order
maxcy = 2
iguess = 1
write(*,106) intl, method, iguess, maxcy, tolmax
call muh2(iparm,fparm,w,iw,cof,bndc,r,u,mgopt,ierror)
write (*,107) ierror
if (ierror.gt.0) call exit(0)
if (ierror .le. 0) then
   !compute and print exact maximum error
   errm = 0.0
   do j=1,np
      sinp = sinap(j)
      cosp = cosap(j)
      do i=1,nt
         sint = sinat(i)
         cost = cosat(i)
         ue = (sint*cosp*sint*sinp*cost)**2
         errm = dmax1(errm,abs((u(i,j)-ue)))
      end do
   end do
   write(*,108) errm
end if

!attempt to improve approximation to fourth order
call muh24(w,iw,u,ierror)
write (*,109) ierror

if (ierror.gt.0) call exit(0)

if (ierror .le. 0) then
   !compute and print exact maximum error
   errm = 0.0
   do j=1,np
      sinp = sinap(j)
      cosp = cosap(j)
      do i=1,nt
         sint = sinat(i)
         cost = cosat(i)
         ue = (sint*cosp*sint*sinp*cost)**2
         errm = dmax1(errm,abs((u(i,j)-ue)))
      end do
   end do
   write(*,108) errm
end if

100 format(//' muh2 test' )
101 format(/' integer input parameters ', &
    &/'intl = ',i2,' nta = ',i2,' ntb = ',i2,' npc = ',i2,' npd = ',i2, &
    &/' itp = ',i2,' jpq = ',i2,' iet = ',i2,' jep = ',i2 &
    &/' nt = ',i3,' np = ',i3,' iguess = ',i2,' maxcy = ',i2, &
    &/' method = ',i2, ' work space estimate = ',i7)

102 format(/' multigrid option parameters ', &
    &/' kcycle = ',i2, &
    &/' iprer = ',i2, &
    &/' ipost = ',i2 &
    &/' intpol = ',i2)

103 format(/' floating point input parameters ', &
    &/' ta = ',f6.3,' tb = ',f6.3,' pc = ',f6.3,' pd = ',f6.3, &
    &/' tolerance (error control) =  ' ,e10.3)
104 format(/' discretization call to muh2 ', ' intl = ',i2)
105 format(' ierror = ',i2, ' minimum work space = ',i7)
106 format(/' approximation call to muh2 ', &
    &/' intl = ',i2, ' method = ',i2 , ' iguess = ',i2, ' maxcy = ',i2 &
    &/' tolmax = ',f5.2)
107 format(' ierror = ',i2 )
108 format(' maximum error  = ',e10.3 )
109 format(/' muh24 test ', ' ierror = ',i2)

end

subroutine cof(t,p,ctt,cpp,ct,cp,ce)
!coefficient subroutine
implicit none
real(8) :: t,p,ctt,cpp,ct,cp,ce
integer :: intl,nta,ntb,npc,npd,itp,jpq,iet,jep,nt,np, &
                iguess,maxcy,method,nwork,lwork,iter
real(8) :: ta,tb,pc,pd,tolmax,sinat,cosat,sinap,cosap,dlt,dlp
common / iprm / intl,nta,ntb,npc,npd,itp,jpq,iet,jep,nt,np, &
                iguess,maxcy,method,nwork,lwork,iter
common / fprm / ta,tb,pc,pd,tolmax,sinat(73),cosat(73), &
                sinap(145),cosap(145),dlt,dlp
integer :: i,j
real(8) :: sinp,cosp,sint,cost,sigma,dsigdt,dsigdp

!set subscripts for current grid point (t,p)
i = int((t-ta)/dlt+0.5)+1
j = int((p-pc)/dlp+0.5)+1

!avoid poles where solution is specified and coefficients are not used

if (i.gt. 1 .and. i.lt. nt) then

   !set sin,cos at (t,p) from precomputed vectors
   sinp = sinap(j)
   cosp = cosap(j)
   sint = sinat(i)
   cost = cosat(i)

   !set sigma and its t,p derivatives
   sigma = 1.5 + (sint*cosp)**2
   dsigdt = 2.0*cosp*cosp*sint*cost
   dsigdp = -2.0*sint*sint*cosp*sinp

   !set coefficients
   ctt = sint*sigma
   cpp = sigma/sint
   ct = sint*dsigdt + cost*sigma
   cp = dsigdp/sint
   ce = -sint*sigma
   return

else

   !set unused coefs at poles arbitrarily
   ctt = 1.0
   cpp = 1.0
   ct = 0.0
   cp = 0.0
   ce = 0.0
   return

end if

end

subroutine bndc(kbdy,torp,alfa,gbdy)

!this subroutine must be provided as a dummy argument even though
!there are no mixed derivative b.c.

return
end

subroutine exact(t,p,utt,upp,ut,up,ue)
!the exact solution used is the restriction of u(x,y,z) = (x*y*z)**2
!in cartesian coordinates to the surface of the sphere of radius one
!using the standard spherical coordinate transforms

common / fprm / ta,tb,pc,pd,tolmax,sinat(73),cosat(73), &
                sinap(145),cosap(145),dlt,dlp

!set subscripts for current grid point (t,p)
i = int((t-ta)/dlt+0.5)+1
j = int((p-pc)/dlp+0.5)+1

!set sin,cos from precomputed vectors
sinp = sinap(j)
cosp = cosap(j)
sint = sinat(i)
cost = cosat(i)

!set intermediate variables

tmp = (sint*cosp*sint*sinp*cost)
dt  = (2.*sint*cost*cost-sint**3)*(cosp*sinp)
dp  = (cosp**2-sinp**2)*(sint**2*cost)
dtt = (2.*cost**3-4.*cost*sint**2-3.*sint**2*cost)*(cosp*sinp)
dpp = (-4.*cosp*sinp)*(sint**2*cost)

!set solution and partial derivatives
ue  = tmp*tmp
ut  = 2.*tmp*dt
up  = 2.*tmp*dp
utt = 2.*(dt*dt+tmp*dtt)
upp = 2.*(dp*dp+tmp*dpp)

return
end
