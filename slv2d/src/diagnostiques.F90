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

#include "selalib.h"
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
sll_real64, dimension(:,:),pointer :: fxy
sll_real64, dimension(:,:),pointer :: fxvx
sll_real64, dimension(:,:),pointer :: fyvy
sll_real64, dimension(:,:),pointer :: fvxvy
sll_int32, private :: i, j, k, l

enum, bind(C)
enumerator :: XY_VIEW = 0, XVX_VIEW  = 1, YVY_VIEW = 2, VXVY_VIEW = 3
end enum

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

 subroutine plot_mesh4d(geomx, geomv, jstartx,jendx,jstartv,jendv)

  use sll_hdf5_io
  type(geometry) :: geomx, geomv
  sll_int32    :: file_id
  sll_int32, intent(in) :: jstartx, jendx, jstartv, jendv
  sll_int32 :: nx, nvx, ny, nvy
  sll_int32 :: error

  nx  = geomx%nx; ny  = geomx%ny
  nvx = geomv%nx; nvy = geomv%ny
  SLL_ALLOCATE(fxy(nx,ny),error)
  SLL_ALLOCATE(fxvx(nx,nvx),error)
  SLL_ALLOCATE(fyvy(ny,jstartv:jendv),error)
  SLL_ALLOCATE(fvxvy(nvx,jstartv:jendv),error)

  call sll_hdf5_file_create("mesh4d.h5",file_id,error)
  call sll_hdf5_write_array(file_id,geomx%xgrid,"/x",error)
  call sll_hdf5_write_array(file_id,geomx%ygrid,"/y",error)
  call sll_hdf5_write_array(file_id,geomv%xgrid,"/vx",error)
  call sll_hdf5_write_array(file_id,geomv%ygrid,"/vy",error)
  call sll_hdf5_file_close(file_id, error)


 end subroutine plot_mesh4d

subroutine plot_df(f4d,iplot,geomx,geomv,jstartx,jendx,jstartv,jendv,choice)
  
  use sll_hdf5_io
  type(geometry), intent(in) :: geomx, geomv
  sll_real64, intent(in) :: f4d(:,:,:,jstartv:)
  sll_int32, intent(in)  :: jstartx, jendx, jstartv, jendv
  
  sll_int32, intent(in) :: iplot
  character(len=4)      :: cplot
  sll_int32         :: file_id
  sll_real64        :: sumloc
  sll_int32         :: nx , nvx, ny, nvy
  sll_int32         :: my_num
  sll_int32         :: error
  sll_int32         :: comm
  sll_int32         :: choice
  
  my_num = sll_get_collective_rank(sll_world_collective)
  comm = sll_world_collective%comm
  
  nx  = geomx%nx; ny  = geomx%ny
  nvx = geomv%nx; nvy = geomv%ny
  
  call int2string(iplot,cplot)
  
  
  select case (choice)
  case(XY_VIEW)
     SLL_ASSERT(size(f4d,1) == size(fxy,1))
     SLL_ASSERT(size(f4d,2) == size(fxy,2))
     do j=1,ny
        do i=1,nx
           sumloc= sum(f4d(i,j,:,jstartv:jendv))
           call mpi_reduce(sumloc,fxy(i,j),1,MPI_REAL8,MPI_SUM,0,comm,error)
        end do
     end do
     if (my_num == 0) then
        call sll_hdf5_file_create('fxy'//cplot//".h5",file_id,error)
        call sll_hdf5_write_array(file_id,fxy,"/fxy",error)
        call write_xdmf("fxy"//cplot//".xmf", &
                        'fxy'//cplot//".h5","x","y","fxy",nx,ny)
        call sll_hdf5_file_close(file_id, error)
     end if
  case(XVX_VIEW)
     SLL_ASSERT(size(f4d,1) == size(fxvx,1))
     SLL_ASSERT(size(f4d,3) == size(fxvx,2))
     do k=1,nvx
        do i=1,nx
           sumloc= sum(f4d(i,:,k,jstartv:jendv))
           call mpi_reduce(sumloc,fxvx(i,k),1,MPI_REAL8,MPI_SUM,0,comm,error)
        end do
     end do
     if (my_num == 0) then
        call sll_hdf5_file_create('fxvx'//cplot//".h5",file_id,error)
        call sll_hdf5_write_array(file_id,fxvx,"/fxvx",error)
        call write_xdmf("fxvx"//cplot//".xmf", &
                        'fxvx'//cplot//".h5","x","vx","fxvx",nx,nvx)
        call sll_hdf5_file_close(file_id, error)
     end if
  case(YVY_VIEW)
     SLL_ASSERT(size(f4d,2) == size(fyvy,1))
     SLL_ASSERT(size(f4d,4) == size(fyvy,2))
     do l=jstartv,jendv
     do j=jstartx,jendx
        fyvy(j,l)= sum(f4d(:,j,:,l))
     end do
     end do
     call write_fyvy(ny,nvy,cplot,jstartv,jendv)
     if (my_num == 0) &
     call write_xdmf("fyvy"//cplot//".xmf", &
                     'fyvy'//cplot//".h5","y","vy","fyvy",ny,nvy)
  case(VXVY_VIEW)
     SLL_ASSERT(size(f4d,3) == size(fvxvy,1))
     SLL_ASSERT(size(f4d,4) == size(fvxvy,2))
     do l=jstartv,jendv
        do k=1,nvx
           fvxvy(k,l)=sum(f4d(:,:,k,l))
        end do
     end do
     call write_fvxvy(nvx,nvy,cplot,jstartv,jendv)
     if (my_num == 0) &
     call write_xdmf("fvxvy"//cplot//".xmf",'fvxvy'//cplot//".h5", &
                     "vx","vy","fvxvy",nvx,nvy)
  end select
  
end subroutine plot_df

 subroutine write_fyvy(ny,nvy,cplot,jstartv, jendv)

 use hdf5
 use sll_hdf5_io_parallel
 character(len=4) :: cplot
 integer(HID_T)    :: pfile_id
 integer(HSSIZE_T) :: offset(2)
 integer(HSIZE_T)  :: global_dims(2)
 sll_int32, intent(in)  :: ny,nvy,jstartv, jendv
 sll_int32 :: error

 global_dims = (/ny,nvy/)
 offset      = (/0, jstartv-1/)
 call sll_hdf5_file_create('fyvy'//cplot//".h5",pfile_id,error)
 call sll_hdf5_write_array_2d(pfile_id,global_dims,offset,fyvy,"/fyvy",error)
 call sll_hdf5_file_close(pfile_id, error)

 end subroutine write_fyvy

 subroutine write_fvxvy(nvx,nvy,cplot,jstartv,jendv)

 use hdf5
 use sll_hdf5_io_parallel
 character(len=4) :: cplot
 integer(HID_T)    :: pfile_id
 integer(HSSIZE_T) :: offset(2)
 integer(HSIZE_T)  :: global_dims(2)
 sll_int32, intent(in)  :: nvx, nvy, jstartv, jendv
 sll_int32 :: error

 global_dims = (/nvx,nvy/)
 offset      = (/0, jstartv-1/)
 call sll_hdf5_file_create('fvxvy'//cplot//".h5",pfile_id,error)
 call sll_hdf5_write_array_2d(pfile_id,global_dims,offset,fvxvy,"/fvxvy",error)
 call sll_hdf5_file_close(pfile_id, error)

 end subroutine write_fvxvy

 subroutine write_xdmf(xdmffilename, datafilename, xname, yname, fname, nx, ny)

  character(len=*), intent(in) :: xdmffilename
  character(len=*), intent(in) :: datafilename
  character(len=*), intent(in) :: xname
  character(len=*), intent(in) :: yname
  character(len=*), intent(in) :: fname

  sll_int32, intent(in) :: nx
  sll_int32, intent(in) :: ny
  sll_int32 :: file_id
  sll_int32 :: error

  call sll_xml_file_create(xdmffilename,file_id,error)
  write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
  write(file_id, &
   "(a,2i5,a)")"<Topology TopologyType='2DRectMesh' NumberOfElements='",ny,nx,"'/>"
  write(file_id,"(a)")"<Geometry GeometryType='VXVY'>"
  write(file_id,"(a,i5,a)")"<DataItem Dimensions='",nx, &
                           "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")"mesh4d.h5:/"//xname
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a,i5,a)")"<DataItem Dimensions='",ny, &
                           "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")"mesh4d.h5:/"//yname
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a)")"</Geometry>"
  write(file_id,"(a)")"<Attribute Name='"//fname//"' AttributeType='Scalar' Center='Node'>"
  write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ny,nx, &
                            "' NumberType='Float' Precision='8' Format='HDF'>"
  write(file_id,"(a)")datafilename//":/"//fname
  write(file_id,"(a)")"</DataItem>"
  write(file_id,"(a)")"</Attribute>"
  call sll_xml_file_close(file_id,error)

 end subroutine write_xdmf

end module diagnostiques_module
