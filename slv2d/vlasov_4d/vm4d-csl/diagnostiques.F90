!=========================================
!   
!    File:          diagnostiquep.F90
!    Project:       vlasov
!    Author(s):     Eric Sonnendrucker
!    Creation:      09.03.1999
!    Last modified: 10.03.1999
!   
!=========================================
!=========================================
module diagnostiques_module

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_file_io.h"
use sll_collective
use used_precision  
use geometry_module

implicit none
integer, parameter :: ithf=50
integer, parameter :: idata = 10
logical :: dir_e
integer, private :: kk0, kk1, kk2, kk3, kk4
character(len=4), private :: fin
character(len=1), private :: aa,bb,cc,dd
character(len=2) :: outdir = './'
sll_int32, private :: i, j, k, l


contains

  subroutine diagnostiques(f,rho,ex,ey,geomx,geomv, &
                           jstartx,jendx,jstartv,jendv,numdiag)

    real(wp), dimension(:,:,:,jstartv:) :: f ! fonc de distribution
!    real(wp), dimension(:,jstartx:)     :: rho  ! densite de charge
!    real(wp), dimension(:,jstartx:)     :: ex,ey ! champ electrique
    real(wp), dimension(:,:)     :: rho  ! densite de charge
    real(wp), dimension(:,:)     :: ex,ey ! champ electrique
    type(geometry) :: geomx, geomv
    integer :: jstartx,jendx,jstartv,jendv
    integer :: numdiag

    ! variables locales
    integer :: iex,iey,irho,ifxvx,ifyvy,ifvxvy
    !integer :: i,j,iv!,jv
    !real(wp) :: sum,sumloc
    !character(2) :: icnum,icpro
    sll_int32 :: comm, my_num, num_threads

    num_threads  = sll_get_collective_size(sll_world_collective)
    my_num = sll_get_collective_rank(sll_world_collective)
    comm = sll_world_collective%comm

    iex=my_num*1000+1
    iey=my_num*1000+2
    irho=my_num*1000+3
    ifxvx=my_num*1000+4
    ifyvy=my_num*1000+5
    ifvxvy=my_num*1000+6

!    WRITE(ICNUM,'(I2.2)') numdiag
!    write(icpro,'(I2.2)') my_num
!    open(iex,file='Runs/ex.'//icnum//'.'//icpro//'.dat')
!    open(iey,file='Runs/ey.'//icnum//'.'//icpro//'.dat')
!    open(irho,file='Runs/rho.'//icnum//'.'//icpro//'.dat')
!    open(ifyvy,file='Runs/fyvy.'//icnum//'.'//icpro//'.dat')
!    open(ifvxvy,file='Runs/fvxvy.'//icnum//'.'//icpro//'.dat')
!    if (my_num.eq.0) then
!       open(ifxvx,file='Runs/fxvx.'//icnum//'.dat')
!    end if
!       
!    ! ecriture du champ electrique
!    write(iex,'(1x,g13.6e3)') ex(1:geomx%nx,jstartx:jendx)
!    write(iey,'(1x,g13.6e3)') ey(1:geomx%nx,jstartx:jendx)
!
!    ! ecriture de rho
!    write(irho,'(1x,g13.6e3)') rho(1:geomx%nx,jstartx:jendx)
!
!    ! ecriture de f(x,vx) - moyenne sur les autres dimensions
!    do iv=1,geomv%nx
!       do i=1,geomx%nx
!          sum=0.  ! initialisation de la somme
!#ifdef _MPI
!          sumloc=0.
!          do jv=jstartv,jendv
!             do j=1,geomx%ny
!                sumloc=sumloc + f(i,j,iv,jv)
!             end do
!          end do
!          call  mpi_reduce(sumloc,sum,1,MPI_REAL8,MPI_SUM,0, &
!               MPI_COMM_WORLD,mpierror)
!          if (my_num.eq.0) then
!             write(ifxvx,'(1x,g13.6e3)') sum !/(geomv%ny*geomx%ny)
!          end if
!#else
!          if(my_num.eq.0) then
!             do jv=1,geomv%ny
!                do j=1,geomx%ny
!                   sum=sum + f(i,j,iv,jv)
!                end do
!             end do
!             write(ifxvx,'(1x,g13.6e3)') sum !/(geomv%ny*geomx%ny)
!          end if
!#endif
!      end do
!    end do
!
!    ! ecriture de f(y,vy) - moyenne sur les autres dimensions
!    do jv=jstartv,jendv
!       do j=1,geomx%ny
!          sum=0.  ! initialisation de la somme
!          do iv=1,geomv%nx
!             do i=1,geomx%nx
!               sum=sum + f(i,j,iv,jv)  
!             end do
!          end do
!          write(ifyvy,'(1x,g13.6e3)') sum !/(geomv%ny*geomx%ny)
!       end do
!    end do
!    ! ecriture de f(vx,vy) - moyenne sur les autres dimensions
!    do jv=jstartv,jendv
!       do iv=1,geomv%nx
!          sum=0.  ! initialisation de la somme
!          do i=1,geomx%nx
!             do j=1,geomx%ny
!                sum=sum + f(i,j,iv,jv)
!             end do
!          end do
!          write(ifvxvy,'(1x,g13.6e3)') sum !/(geomx%nx*geomv%ny)
!       end do
!    end do
!
!    close(iex)
!    close(iey)
!    close(irho)
!    close(ifxvx)
!    close(ifyvy)
!    close(ifvxvy)

  end subroutine diagnostiques

  subroutine fichinit
    character(len=72) :: filename
    integer :: IO_stat ! indicateur d'erreur
    
    call getarg( 1, filename)

    open(idata,file=trim(filename),IOStat=IO_stat)
    if (IO_stat/=0) STOP "Miss argument file.nml"
    open(ithf,file="thf.dat",IOStat=IO_stat)
    if (IO_stat/=0) STOP "erreur d'ouverture du fichier thf.dat"
  end subroutine fichinit
  
  subroutine time_history(desc,format,array)
    character(3) :: desc
    character(14) :: format
    real(wp), dimension(:) :: array
    
    if (desc(1:3)=="thf") then
       !print *,'array', array
       write(ithf,format) array
    else
       write(*,*) desc," not recognized"
    endif
    
  end subroutine time_history



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
!integer :: iex,iey,irho,ifxvx,ifyvy,ifvxvy,mpierror
integer :: ifvxvy
integer :: i,j,iv,jv
real(wp) :: sum !,sumloc
!character(2) :: icnum,icpro

integer :: ifile, k, file_id, error
character(len=2), dimension(5) :: which 
sll_int32 :: comm, my_num, num_threads

num_threads  = sll_get_collective_size(sll_world_collective)
my_num = sll_get_collective_rank(sll_world_collective)
comm = sll_world_collective%comm

#ifdef _SILO
!call write_domains(my_num,geomx,geomv,f,rho,ex,ey,bz,jx,jy,    &
!          numdiag,jstartx,jendx,jstartv,jendv)
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
     call system("mkdir -p "//trim(outdir)//char(my_num+48))
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
                       sngl(bz(i,j)),           &
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
      write(ifile,"(a)",advance='no')   &
      "splot '"//trim(outdir)//"0/"//fin//".dat' u 1:2:"//char(50+k)//" w l "
      
      do j = 1, num_threads - 1
         write(ifile,"(a)",advance='no')    &
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
      write(ifvxvy,"(a)",advance='no')  &
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

end module diagnostiques_module
