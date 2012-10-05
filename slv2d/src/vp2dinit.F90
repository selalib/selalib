module vp2dinit
  use mpi
  use geometry_module
  use vlasov2d_module
  use splinepx_class
  use splinepy_class
  use poisson2dpp_module
  use diagnostiques_module

  implicit none
contains
  subroutine initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)
    !-------------------------------------------------------
    ! sert a l'initialisation globale du programme VP2D
    !-------------------------------------------------------
    type(geometry)  :: geomx, geomv  ! geometrie globale du probleme
    real(wp)     :: dt       ! pas de temps
    integer         :: nbiter   ! nombre d'iterations en temps
    integer         :: fdiag    ! frequences des diagnostiques
    integer         :: fthdiag    ! frequences des historiques en temps
    integer         :: nx, ny   ! dimensions de l'espace physique
    integer         :: nvx, nvy ! dimensions de l'espace des vitesses
    real(wp)     :: x0, y0   ! coordonnees debut du maillage espace physique
    real(wp)     :: vx0, vy0 ! coordonnees debut du maillage espace vitesses
    real(wp)     :: x1, y1   ! coordonnees fin du maillage espace physique
    real(wp)     :: vx1, vy1 ! coordonnees fin du maillage espace vitesses
    integer :: iflag,ierr  ! indicateur d'erreur
    ! definition of namelists
    namelist /time/ dt, nbiter
    namelist /diag/ fdiag, fthdiag! freq. of diags and time hist diags in steps
    namelist /phys_space/ x0,x1,y0,y1,nx,ny
    namelist /vel_space/ vx0,vx1,vy0,vy1,nvx,nvy
   
#ifdef _MPI
    if (my_num.eq.0) then
#endif
       ! open files
       call fichinit
       ! read input data
       read(idata,NML=time)
       read(idata,NML=diag)
       read(idata,NML=phys_space)
       read(idata,NML=vel_space)
#ifdef _MPI
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)

    call mpi_bcast(dt,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nbiter,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(fdiag,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(fthdiag,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(x0,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(y0,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(x1,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(y1,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nx,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ny,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vx0,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vy0,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vx1,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vy1,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nvx,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nvy,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
#endif

! initialisation de la geometrie de l'espace physique
    call new_geometry2(geomx,x0,y0,x1,y1,nx,ny,iflag,"perxy")
    if (iflag.ne.0) stop 'erreur dans l initialisation de geomx'

    ! initialisation de la geometrie de l'espace des vitesses
    call new(geomv,vx0,vy0,vx1,vy1,nvx,nvy,iflag,"natxy")
    if (iflag.ne.0) stop 'erreur dans l initialisation de geomv'

  end subroutine initglobal

  subroutine initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
       f,rho,ex,ey,vlas2d,poiss2dpp,splx,sply)
    !-------------------------------------------------------
    ! sert a l'initialisation en parallele du programme VP2D
    !-------------------------------------------------------
    type(geometry) :: geomx, geomv  ! geometrie globale du probleme
    integer :: jstartv,jendv,jstartx,jendx ! definition de la bande du proc
    real(wp), dimension(:,:,:,:),pointer :: f
    real(wp), dimension(:,:),pointer :: rho,ex,ey
    integer :: proc, mpierror
    type(vlasov2d) :: vlas2d
    type(splinepx) :: splx
    type(splinepy) :: sply
    type(poisson2dpp) :: poiss2dpp
    !variables locales
#if defined _OPENMP || defined _MPI
    integer :: ipiece_size_v, ipiece_size_x
#endif
    real(wp) :: xi,vx,vy,v2,x,y,eps,kx,ky
    integer :: i,j,iv,jv,iflag

    ! cas sequentiel
    jstartv=1
    jendv=geomv%nx
    jstartx=1
    jendx=geomx%ny

#if defined _OPENMP || defined _MPI
    ! initialisation of size of parallel zones 
    ! the total size of the vx zone is nvx
    ! the total size of the y zone is ny
    ! ipiece_size = n/num_threads rounded up
    ipiece_size_v = (geomv%ny + num_threads - 1) / num_threads
    ipiece_size_x = (geomx%ny + (num_threads-1)) / num_threads
    ! zone a traiter en fonction du numero de process
    jstartv = my_num * ipiece_size_v + 1
    jendv = min(jstartv - 1 + ipiece_size_v, geomv%ny)
    jstartx = my_num * ipiece_size_x + 1
    jendx = min(jstartx - 1 + ipiece_size_x, geomx%ny)
    
    if (jstartv.gt.jendv) stop "error in vp2ddinit: negative size zone"
    if (jstartx.gt.jendx) stop "error in vp2dinit: negative size zone"
      print*,'init zone ',my_num,jstartx,jendx,ipiece_size_x, jstartv,jendv,ipiece_size_v
#endif
#ifdef _MPI
! allocation dynamique des tableaux
  allocate(f(geomx%nx,geomx%ny,geomv%nx,jstartv:jendv),stat=iflag)
  if (iflag.ne.0) stop 'erreur dans l allocation de f'
#endif
! Poisson n'est pas parallele pour l'instant
  allocate(rho(geomx%nx,geomx%ny),stat=iflag)
  if (iflag.ne.0) stop 'erreur dans l allocation de rho'
  allocate(ex(geomx%nx,geomx%ny),stat=iflag)
  if (iflag.ne.0) stop 'erreur dans l allocation de ex'
  allocate(ey(geomx%nx,geomx%ny),stat=iflag)
  if (iflag.ne.0) stop 'erreur dans l allocation de ey'

    ! initialisation parallele des tableaux globaux, 
    ! ce qui permet  de les distribuer sur les processeurs
    ! initialisation de la fonction de distribution 
    xi=0.90
    eps=0.05
    kx=2*pi/((geomx%nx)*geomx%dx)
    ky=2*pi/((geomx%ny)*geomx%dy)
    do jv=jstartv,jendv
       vy = geomv%y0+(jv-1)*geomv%dy
       do iv=1,geomv%nx
          vx = geomv%x0+(iv-1)*geomv%dx
          v2 = vx*vx+vy*vy
          do j=1,geomx%ny
             y=geomx%y0+(j-1)*geomx%dy
             do i=1,geomx%nx
                x=geomx%x0+(i-1)*geomx%dx
                f(i,j,iv,jv)=(1+eps*cos(kx*x))*1/(2*pi)*exp(-.5*v2)
             end do
          end do
       end do
    end do

    ! initialisation de ex,ey, rho
    ex = 0.0
    ey = 0.0
    rho = 0.0

    ! initialisation du calcul du champ electrique
    call new(poiss2dpp,rho,geomx,iflag, jstartx, jendx)
    ! initialisation du module vlasov
    call new(vlas2d,geomx,geomv,iflag,jstartx, jendx,jstartv,jendv)
    ! initialisation du module vlasov1D	       
    call new(splx,geomx,geomv,iflag,jstartx,jendx,jstartv,jendv)
    call new(sply,geomx,geomv,iflag,jstartx,jendx,jstartv,jendv)
  end subroutine initlocal

end module vp2dinit
