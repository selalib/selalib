module silo_module
!=========================================
!    File:          silo_module.F90
!    Project:       vlasov
!    Author(s):     Pierre Navaro
!    Creation:      08.01.2010
!=========================================

use used_precision  
use geometry_module
!use Module_MPI, my_num=>numero_processeur, num_threads=>nombre_processeurs

#ifdef _SILO

implicit none
include "silo.inc"
integer, private :: kk0, kk1, kk2, kk3, kk4
character(len=4), private :: fin
character(len=1), private :: aa,bb,cc,dd
character(len=2) :: outdir = './'

contains

subroutine write_domains(imesh,geomx,geomv,f,rho,ex,ey,bz,jx,jy, &
			 iplot,jstartx,jendx,jstartv,jendv)
type (geometry) :: geomx, geomv
integer :: err, ierr, dom, ndims, imesh, iplot
integer :: dims(2), mx, my, dbfile
real(8), dimension(:,:,:,jstartv:), intent(in) :: f
real(8), dimension(:,:), intent(in) :: ex, ey, bz, jx, jy
real(8), dimension(:,:), intent(in) :: rho
integer :: i, j, k, iv, jv, kv, n
real(4), allocatable :: x(:), y(:), vx(:), vy(:)
real(4), allocatable :: z(:,:)
real(4) :: sumloc
character(len=40) :: filename 
integer :: lfilename, jstartx, jendx, jstartv, jendv
integer :: mpierror

kk0 = iplot
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd

! Poke a number into the filename.
filename = trim(outdir)//char(48+imesh)//"/"//fin
lfilename = len_trim(filename)
! Create a new silo file.
ierr = dbcreate(filename, lfilename, DB_CLOBBER, DB_LOCAL,	&
                "multimesh data", 14, DB_HDF5, dbfile)
if (dbfile.eq.-1) then
   write (6,*) 'Could not create Silo file!\n'
   return
end if
! Displace the coordinates
n = max(geomx%nx,geomx%ny,geomv%nx,geomv%ny)
allocate(x(n),y(n),z(n,n),vx(n),vy(n))

k = 0
do i=jstartx, jendx
   k = k + 1
   x(k) = sngl(geomx%x0+(i-1)*geomx%dx)
end do
do j=1,geomx%ny
   y(j) = sngl(geomx%y0+(j-1)*geomx%dy)
end do
do iv=1,geomv%nx
   vx(iv) = sngl(geomv%x0+(iv-1)*geomv%dx)
end do
kv = 0
do jv=jstartv,jendv
   kv = kv + 1
   vy(kv) = sngl(geomv%y0+(jv-1)*geomv%dy)
end do

! Write the multimeshes
ndims = 2
do iv=1,geomv%nx
   k = 0
   do i=jstartx,jendx
      k = k + 1
      z(k,iv) = 0.
      do jv=jstartv,jendv
         do j=1,geomx%ny
            z(k,iv) = z(k,iv) + sngl(f(i,j,iv,jv))
         end do
      end do
  end do
end do

dims(1) = jendx - jstartx + 1
dims(2) = geomv%nx
err = dbputqm (dbfile, "geomxvx", 7, "x", 1, "vx", 2, "z", 1, 	&
	       x, vx, DB_F77NULL, dims, ndims,			&	
     	       DB_FLOAT, DB_COLLINEAR, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "fxvx", 4, "geomxvx", 7, z(:,:), dims, 	&
               ndims, DB_F77NULL, 0, DB_FLOAT, DB_NODECENT, 	&
	       DB_F77NULL, ierr)

do j=1,geomx%ny
   kv = 0
   do jv=jstartv,jendv
      kv = kv + 1
      z(j,kv) = 0.
      do iv=1,geomv%nx
         do i=jstartx,jendx
           z(j,kv) = z(j,kv) + sngl(f(i,j,iv,jv))
         end do
      end do
   end do
end do

dims(1) = geomx%ny
dims(2) = jendv - jstartv + 1
err = dbputqm (dbfile, "geomyvy", 7, "y", 1, "vy", 2, "z", 1, 	&
	       y, vy, DB_F77NULL, dims, ndims,			&	
     	       DB_FLOAT, DB_COLLINEAR, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "fyvy", 4, "geomyvy", 7, z(:,:), dims, 	&
               ndims, DB_F77NULL, 0, DB_FLOAT, DB_NODECENT, 	&
	       DB_F77NULL, ierr)

do iv=1,geomv%nx
   kv = 0
   do jv=jstartv,jendv
      kv = kv+1
      z(iv,kv) = 0.  
      do i=1,geomx%nx
         do j=1,geomx%ny
            z(iv,kv) = z(iv,kv) + sngl(f(i,j,iv,jv))
         end do
      end do
   end do
end do

dims(1) = geomv%nx
dims(2) = jendv - jstartv + 1
err = dbputqm (dbfile, "geomvxvy", 8, "vx", 2,			&
     	       "vy", 2, "z", 1, vx, vy, DB_F77NULL, dims, ndims,&	
     	       DB_FLOAT, DB_COLLINEAR, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "fvxvy", 5, "geomvxvy", 8, z(:,:), dims, &
               ndims, DB_F77NULL, 0, DB_FLOAT, DB_NODECENT, 	&
	       DB_F77NULL, ierr)

do j=1,geomx%ny
   k = 0
   do i=jstartx,jendx
      k = k+1
      z(k,j) = 0.  
      do jv=jstartv,jendv
         do iv=1,geomv%nx
            z(k,j) = z(k,j) + sngl(f(i,j,iv,jv))
         end do
      end do
   end do
end do

dims(1) = jendx - jstartx + 1
dims(2) = geomx%ny
err = dbputqm (dbfile, "geomxy", 6, "x", 1,			&
      "y", 1, "z", 1, x, y, DB_F77NULL, dims, ndims,	&	
      DB_FLOAT, DB_COLLINEAR, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "fxy", 3, "geomxy", 6, z(:,:), dims, 	&
      ndims, DB_F77NULL, 0, DB_FLOAT, DB_NODECENT, DB_F77NULL, ierr)

!Ecritures des champs electro-magnetiques
err = dbputqv1(dbfile, "Ex", 2, "geomxy", 6, ex(jstartx:jendx,:), dims, &
      ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "Ey", 2, "geomxy", 6, ey(jstartx:jendx,:), dims, &
      ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "Bz", 2, "geomxy", 6, bz(jstartx:jendx,:), dims, &
      ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "Jx", 2, "geomxy", 6, jx(jstartx:jendx,:), dims, &
      ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "Jy", 2, "geomxy", 6, jy(jstartx:jendx,:), dims, &
      ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
err = dbputqv1(dbfile, "Ro", 2, "geomxy", 6, rho(jstartx:jendx,:), dims, &
      ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)

! Close the Silo file
ierr = dbclose(dbfile)

end subroutine write_domains

subroutine write_master(nmesh, iplot)

integer :: err, ierr, nmesh, oldlen, nvar, dbfile, iplot, i
character(len=40) :: filename
integer :: lfilename
character(len=20), dimension(nmesh) :: meshnames
character(len=20), dimension(nmesh) :: varnames1, varnames2, varnames3
character(len=20), dimension(nmesh) :: varnames4, varnames5, varnames6
integer, dimension(nmesh) :: meshtypes, lmeshnames
integer, dimension(nmesh) :: lvarnames, vartypes

kk0 = iplot
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd

! Create a new silo file
filename = trim(outdir)//fin//".root"
lfilename = len_trim(filename)
ierr = dbcreate(filename,lfilename, DB_CLOBBER, DB_LOCAL, &
                "multimesh root", 14, DB_HDF5, dbfile)

if(dbfile.eq.-1) then
   write (6,*) "Could not create Silo file!\n"
   return
endif

! Set the maximum string length to 20 since that's how long our
! strings are
oldlen = dbget2dstrlen()
err = dbset2dstrlen(20)

do i = 1, nmesh
   meshnames(i)  = char(i+47)//"/"//fin//":geomxy"
   meshtypes(i)  = DB_QUAD_RECT
   lmeshnames(i) = len_trim(meshnames(i))
   varnames1(i)  = char(i+47)//"/"//fin//":Ex"
   varnames2(i)  = char(i+47)//"/"//fin//":Ey"
   varnames3(i)  = char(i+47)//"/"//fin//":Bz"
   varnames4(i)  = char(i+47)//"/"//fin//":Jx"
   varnames5(i)  = char(i+47)//"/"//fin//":Jy"
   varnames6(i)  = char(i+47)//"/"//fin//":Ro"
   lvarnames(i)  = 9
   vartypes(i)   = DB_QUADVAR
end do

! Write the multimesh object.
err = dbputmmesh(dbfile, "geomxy", 6, nmesh, meshnames,	&
                 lmeshnames, meshtypes, DB_F77NULL, ierr)
nvar = nmesh
err = dbputmvar(dbfile, "Ex", 2, nvar, varnames1, lvarnames,	&
                vartypes, DB_F77NULL, ierr)
err = dbputmvar(dbfile, "Ey", 2, nvar, varnames2, lvarnames,	&
                vartypes, DB_F77NULL, ierr)
err = dbputmvar(dbfile, "Bz", 2, nvar, varnames3, lvarnames,	&
                vartypes, DB_F77NULL, ierr)
err = dbputmvar(dbfile, "Jx", 2, nvar, varnames4, lvarnames,	&
                vartypes, DB_F77NULL, ierr)
err = dbputmvar(dbfile, "Jy", 2, nvar, varnames5, lvarnames,	&
                vartypes, DB_F77NULL, ierr)
err = dbputmvar(dbfile, "Ro", 2, nvar, varnames6, lvarnames,	&
                vartypes, DB_F77NULL, ierr)

do i = 1, nmesh
   meshnames(i)  = char(i+47)//"/"//fin//":geomxvx"
   meshtypes(i)  = DB_QUAD_RECT
   lmeshnames(i) = len_trim(meshnames(i))
   varnames1(i)  = char(i+47)//"/"//fin//":fxvx"
   lvarnames(i)  = len_trim(varnames1(i))
   vartypes(i)   = DB_QUADVAR
end do

err = dbputmmesh(dbfile, "geomxvx", 7, nmesh, meshnames,	&
                 lmeshnames, meshtypes, DB_F77NULL, ierr)
nvar = nmesh
err = dbputmvar(dbfile, "fxvx", 4, nvar, varnames1, lvarnames,	&
                vartypes, DB_F77NULL, ierr)

do i = 1, nmesh
   meshnames(i)  = char(i+47)//"/"//fin//":geomyvy"
   meshtypes(i)  = DB_QUAD_RECT
   lmeshnames(i) = len_trim(meshnames(i))
   varnames1(i)  = char(i+47)//"/"//fin//":fyvy"
   lvarnames(i)  = len_trim(varnames1(i))
   vartypes(i)   = DB_QUADVAR
end do

err = dbputmmesh(dbfile, "geomyvy", 7, nmesh, meshnames,	&
                 lmeshnames, meshtypes, DB_F77NULL, ierr)
nvar = nmesh
err = dbputmvar(dbfile, "fyvy", 4, nvar, varnames1, lvarnames,	&
                vartypes, DB_F77NULL, ierr)

do i = 1, nmesh
   meshnames(i)  = char(i+47)//"/"//fin//":geomvxvy"
   meshtypes(i)  = DB_QUAD_RECT
   lmeshnames(i) = len_trim(meshnames(i))
   varnames1(i)  = char(i+47)//"/"//fin//":fvxvy"
   lvarnames(i)  = len_trim(varnames1(i))
   vartypes(i)   = DB_QUADVAR
end do

err = dbputmmesh(dbfile, "geomvxvy", 8, nmesh, meshnames,	&
                 lmeshnames, meshtypes, DB_F77NULL, ierr)
nvar = nmesh
err = dbputmvar(dbfile, "fvxvy", 5, nvar, varnames1, lvarnames,	&
                vartypes, DB_F77NULL, ierr)

! Restore the previous value for maximum string length
err = dbset2dstrlen(oldlen)
! Close the Silo file
ierr = dbclose(dbfile)

end subroutine write_master

#endif

end module silo_module
