module mgd3
 
   integer, dimension(20)   :: nxk,nyk,nzk,sxk,exk,syk,eyk,szk,ezk
   integer, dimension(20)   :: kpbgn,kcbgn
   integer, dimension(7,20) :: kdatatype
   integer, dimension(20)   :: sxi,exi,syi,eyi,szi,ezi
   integer, dimension(20)   :: nxr,nyr,nzr,sxr,exr,syr,eyr,szr,ezr
   integer, dimension(7,20) :: rdatatype

contains


subroutine initialize()
!>
!> initialize mgd3
!>
!call mgdinit(vbc,phibc,ixp,jyq,kzr,iex,jey,kez,ngrid,nxp2,	&
!             nyp2,nzp2,sx,ex,sy,ey,sz,ez,realtype,nxprocs,	&
!             nyprocs,nzprocs,nwork,ibdry,jbdry,kbdry,myid,	&
!             iout,nerror)
!if (nerror.eq.1) goto 1000
!>
!> initialize problem
!> xl,yl,zl are the dimensions of the domain
!> wk is the wavenumber (must be an integer value)
!> rro is the average density
!> 1/hxi,1/hyi,1/hzi are the spatial resolutions
!
!xl=1.0d0
!yl=1.0d0
!zl=1.0d0
!wk=5.0d0
!rro=1.0d0
!hxi=float(nxp2-2)/xl
!hyi=float(nyp2-2)/yl
!hzi=float(nzp2-2)/zl
!write(iout,*) 'hxi=',hxi,' hyi=',hyi,' hzi=',hzi
!call ginit(sx,ex,sy,ey,sz,ez,p,r,f,wk,hxi,hyi,hzi,pi,iout)
end subroutine initialize

subroutine solve()
!>
!> solve using mgd3
!>
!call mgdsolver(2,sx,ex,sy,ey,sz,ez,p,f,r,ngrid,work,	&
!               maxcy,tolmax,kcycle,iprer,ipost,iresw,	&
!               xl,yl,zl,rro,nx,ny,nz,comm3d,comm3dp,	&
!               comm3dl,comm3dc,myid,neighbor,bd,phibc,	&
!               iter,.true.,iout,nerror)
!if (nerror.eq.1) goto 1000
end subroutine solve

end module mgd3
