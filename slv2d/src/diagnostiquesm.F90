module diagnostiquesm_module
!=========================================
!    File:          diagnostiquem.f90
!    Project:       vlasov
!    Author(s):     Pierre Navaro
!    Creation:      08.01.2010
!=========================================
#include "selalib.h"
use used_precision  
use geometry_module
use Module_MPI, my_num=>numero_processeur, num_threads=>nombre_processeurs

implicit none

integer, parameter :: idata = 10
logical :: dir_e
integer, private :: kk0, kk1, kk2, kk3, kk4
character(len=4), private :: fin
character(len=1), private :: aa,bb,cc,dd
character(len=2) :: outdir = './'

contains

subroutine diagnostiquesm(f,rho,ex,ey,bz,jx,jy,geomx,geomv, &
	                  jstartx,jendx,jstartv,jendv,numdiag)

real(wp), dimension(:,:,:,jstartv:) :: f ! fonc de distribution
!    real(wp), dimension(:,jstartx:)     :: rho  ! densite de charge
!    real(wp), dimension(:,jstartx:)     :: ex,ey ! champ electrique
real(wp), dimension(:,:)     :: ex,ey ! champ electrique
real(wp), dimension(:,:)     :: jx,jy ! courants
real(wp), dimension(:,:)     :: bz    ! champ magnetique
real(wp), dimension(:,:)     :: rho   
type(geometry) :: geomx, geomv
integer :: jstartx,jendx,jstartv,jendv
integer :: numdiag

! variables locales
integer :: iex,iey,irho,ifxvx,ifyvy,ifvxvy,mpierror
integer :: i,j,iv,jv
real(wp) :: sum,sumloc
character(2) :: icnum,icpro

integer :: ifile, k, file_id, error
character(len=2), dimension(5) :: which 

#ifdef _SILO
!call write_domains(my_num,geomx,geomv,f,rho,ex,ey,bz,jx,jy,	&
!		   numdiag,jstartx,jendx,jstartv,jendv)
!
!if (my_num == 0) call write_master(num_threads,numdiag)
#endif
 
which(1) = 'Ex'
which(2) = 'Ey'
which(3) = 'Bz'
which(4) = 'Jx'
which(5) = 'Jy'

if (numdiag == 0) then
   inquire(file=trim(outdir)//char(my_num+48)//"/"".", exist=dir_e)
   if (dir_e) then
     write(*,*) "directory "//trim(outdir)//char(my_num+48)//" exists!"
   else
     !call system("mkdir -p "//trim(outdir)//char(my_num+48))
   end if
end if

kk0 = numdiag
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd

call sll_new_file_id(file_id, error)

open(file_id, file=trim(outdir)//char(my_num+48)//"/"//fin//".dat")

do i=jstartx,jendx
   do j=1,geomx%ny
      write(file_id,*) sngl(geomx%x0+(i-1)*geomx%dx), &
                       sngl(geomx%y0+(j-1)*geomx%dy), &
                       sngl(ex(i,j)), sngl(ey(i,j)),  &
                       sngl(bz(i,j)),		    &
                       sngl(jx(i,j)), sngl(jy(i,j))
   end do
   write(file_id,*) 
end do
close(file_id)
   
if (my_num == MPI_MASTER) then
   do k = 1, 5
      ifile = 90+k
      open( ifile, file = which(k)//'.gnu', position="append" )
      if ( numdiag == 1 ) then
         rewind(ifile)
         !write(ifile,*)"set xr[-0.1:1.1]"
         !write(ifile,*)"set yr[-0.1:1.1]"
         !write(ifile,*)"set zr[-1.1:1.1]"
         !write(ifile,*)"set cbrange[-1:1]"
         !write(ifile,*)"set pm3d"
         !write(ifile,*)"set surf"
         !write(ifile,*)"set term x11"
      end if
      write(ifile,*)"set title 'Time = ",numdiag,"'"
      write(ifile,"(a)",advance='no')	&
      "splot '"//trim(outdir)//"0/"//fin//".dat' u 1:2:"//char(50+k)//" w l "
      
      do j = 1, num_threads - 1
         write(ifile,"(a)",advance='no')	&
         ", '"//trim(outdir)//char(j+48)//"/"//fin//".dat' u 1:2:"//char(50+k)//" w l "
      end do
      write(ifile,*)
      !write(ifile,*)"set term gif"
      !write(ifile,*)"set output 'image"//fin//".gif'"
      !write(ifile,*)"replot"
      close(ifile)
   end do

end if

ifvxvy = my_num*100+7
if (my_num == 0) then
   open(ifvxvy, file = 'fvxvy.gnu', position="append" )
   if ( numdiag == 1 ) then
      rewind(ifvxvy)
   end if
   write(ifvxvy,*)"set title 'Time = ",numdiag,"'"
   write(ifvxvy,"(a)",advance='no') "splot '"//trim(outdir)//"0/fvxvy-"//fin//".dat' w l"
   do j = 1, num_threads - 1
      write(ifvxvy,"(a)",advance='no')	&
      ", '"//trim(outdir)//char(j+48)//"/fvxvy-"//fin//".dat' w l "
   end do
   write(ifvxvy,*)
   close(ifvxvy)
end if
open(ifvxvy,file=trim(outdir)//char(my_num+48)//"/fvxvy-"//fin//'.dat')
do jv=jstartv,jendv
   do iv=1,geomv%nx
      sum=0.  ! initialisation de la somme
      do i=1,geomx%nx
         do j=1,geomx%ny
            sum=sum + f(i,j,iv,jv)
         end do
      end do
      write(ifvxvy,*) iv, jv, sum !/(geomx%nx*geomv%ny)
   end do
   write(ifvxvy,*) 
end do
close(ifvxvy)

end subroutine diagnostiquesm

end module diagnostiquesm_module
