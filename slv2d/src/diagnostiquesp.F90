module diagnostiques_module
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

  use used_precision  
  use geometry_module
#ifdef _MPI
  use Module_MPI, my_num=>numero_processeur, num_threads=>nombre_processeurs
#endif
  implicit none
  integer,parameter :: ithf=50, idata=10
contains
  subroutine diagnostiques(f,rho,ex,ey,geomx,geomv, &
	jstartx,jendx,jstartv,jendv,numdiag)
#ifdef _MPI
    real(wp), dimension(:,:,:,jstartv:) :: f ! fonc de distribution
!    real(wp), dimension(:,jstartx:)     :: rho  ! densite de charge
!    real(wp), dimension(:,jstartx:)     :: ex,ey ! champ electrique
    real(wp), dimension(:,:)     :: rho  ! densite de charge
    real(wp), dimension(:,:)     :: ex,ey ! champ electrique
#else
    real(wp), dimension(:,:,:,:) :: f ! fonc de distribution
    real(wp), dimension(:,:)     :: rho  ! densite de charge
    real(wp), dimension(:,:)     :: ex,ey ! champ electrique
#endif
    type(geometry) :: geomx, geomv
    integer :: jstartx,jendx,jstartv,jendv
    integer :: numdiag

    ! variables locales
    integer :: iex,iey,irho,ifxvx,ifyvy,ifvxvy,mpierror
    integer :: i,j,iv,jv
    real(wp) :: sum,sumloc
    character(2) :: icnum,icpro
#ifdef _OPENMP
    integer my_num

    my_num = omp_get_thread_num()
#endif

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
!          call  mpi_reduce(sumloc,sum,1,MPI_realtype,MPI_SUM,0, &
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
    integer :: IO_stat ! indicateur d'erreur
    
    open(idata,file="slv2d.dat",IOStat=IO_stat)
    if (IO_stat/=0) STOP "erreur d'ouverture du fichier mabo2d.dat"
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

end module diagnostiques_module
