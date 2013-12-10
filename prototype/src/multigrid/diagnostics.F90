module diagnostics
implicit none

integer, private :: i, j, k

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gp_plot2d(rang, nproc, f, nx, ny, hx, hy, xp, yp )

integer, intent(in) :: rang, nx, ny, nproc
real, intent(in), dimension(:,:) :: f
real :: xp, yp, hx, hy

!write domains
open( 80, file = "p"//char(rang+48)//".dat" )
   do i= 1, nx
      do j= 1, ny
         write(80,"(3e15.5)") xp+(i-1)/hx, yp+(j-1)/hy, f(i,j)
      end do
      write(80,*) 
   end do
close(80)
   
!write master file
if (rang == 0) then
   open( 90, file = 'p.gnu', position="append" )
   write(90,"(a)",advance='no')"splot 'p0.dat' w lines"
   do j = 1, nproc - 1
      write(90,"(a)",advance='no') ",'p"//char(j+48)//".dat' w lines" 
   end do
   write(90,*)
   close(90)
end if

end subroutine gp_plot2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!subroutine silo_plot2d(rank, var, hx, hy, nx, ny, xp, yp)
!include "silo.inc"
!integer, intent(in) :: nx, ny
!integer :: dbfile, err, ierr, i,j, dims(2), ndims, rank
!real, intent(in), dimension(nx,ny) :: var
!real :: xc(nx), yc(ny)
!real, intent(in) :: hx, hy, xp, yp
!character(len=10) :: filename='multivar.X'
!
!ndims = 2
!dims(1) = nx
!dims(2) = ny
!! Poke a number into the filename.
!filename(10:) = char(49 + rank)
!! Create a new silo file.
!ierr = dbcreate(filename, 10, DB_CLOBBER, DB_LOCAL, &
!                "multivar data", 13, DB_HDF5, dbfile)
!if(dbfile.eq.-1) then
!   write (6,*) 'Could not create Silo file!\n'
!   return
!endif
!! Displace the coordinates
!
!do i=1,nx
!   xc(i) = xp + (i-1)/hx
!end do
!do j=1,ny
!   yc(j) = yp + (j-1)/hy
!end do
!write(*,*) "rang = ", rank, " nx, ny = ", nx, ny
!write(*,*) "rang = ", rank, kind(xp)
!! Write the quadmesh and the quadvar
!if (kind(xp) == 8) then
!   err = dbputqm (dbfile, "quadmesh", 8, "xc", 2, "yc", 2, "zc", 2, &
!         xc, yc, DB_F77NULL, dims, ndims, DB_DOUBLE, DB_COLLINEAR,  &
!         DB_F77NULL, ierr)
!   err = dbputqv1(dbfile, "var", 3, "quadmesh", 8, var, dims, ndims,&
!         DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
!else
!   err = dbputqm (dbfile, "quadmesh", 8, "xc", 2, "yc", 2, "zc", 2, &
!         xc, yc, DB_F77NULL, dims, ndims, DB_FLOAT, DB_COLLINEAR,   &
!         DB_F77NULL, ierr)
!   err = dbputqv1(dbfile, "var", 3, "quadmesh", 8, var, dims, ndims,&
!         DB_F77NULL, 0, DB_FLOAT, DB_NODECENT, DB_F77NULL, ierr)
!end if
!if (rank == 0) call write_master()
!! Close the Silo file
!ierr = dbclose(dbfile)
!end subroutine silo_plot2d


!subroutine write_multimesh(dbfile)
!implicit none
!include "silo.inc"
!integer, intent(in) :: dbfile
!integer :: err, ierr, nmesh
!character*20, dimension(4) :: meshnames
!integer, dimension(4) :: lmeshnames 
!integer, dimension(4) :: meshtypes 
!nmesh = 4
!meshnames =(/'multivar.1:quadmesh ', &
!             'multivar.2:quadmesh ', &
!             'multivar.3:quadmesh ', &
!             'multivar.4:quadmesh '/)
!lmeshnames = 19
!meshtypes  = DB_QUAD_RECT
!err = dbputmmesh(dbfile, "quadmesh", 8, nmesh, meshnames, &
!      lmeshnames, meshtypes, DB_F77NULL, ierr)
!end subroutine write_multimesh
!
!subroutine write_multivar(dbfile)
!implicit none
!include "silo.inc"
!integer, intent(in) :: dbfile
!integer err, ierr, nvar
!character*20 varnames(4) 
!integer :: lvarnames(4)
!integer :: vartypes(4)
!nvar = 4
!varnames =(/'multivar.1:var     ', &
!            'multivar.2:var     ', &
!            'multivar.3:var     ', &
!            'multivar.4:var     '/)
!lvarnames = 14
!vartypes = DB_QUADVAR
!err = dbputmvar(dbfile, "var", 3, nvar, varnames, lvarnames,vartypes, DB_F77NULL, ierr)
!end subroutine write_multivar
!
!subroutine write_master()
!implicit none
!include "silo.inc"
!integer err, ierr, dbfile, oldlen
!!Create a new silo file
!ierr = dbcreate("multivar.root", 13, DB_CLOBBER, DB_LOCAL, &
! "multimesh root", 14, DB_HDF5, dbfile)
!if(dbfile.eq.-1) then
!   write (6,*) 'Could not create Silo file!\n'
!   return
!endif
!!Set the maximum string length to 20 since that's how long our strings are
!oldlen = dbget2dstrlen()
!err = dbset2dstrlen(20)
!
!! Write the multimesh and multivar objects
!call write_multimesh(dbfile)
!call write_multivar(dbfile)
!
!! Restore the previous value for maximum string length
!err = dbset2dstrlen(oldlen)
!! Close the Silo file
!ierr = dbclose(dbfile)
!end subroutine write_master
!
end module diagnostics
