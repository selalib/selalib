program VP2D
  !-------------------------------------------------------------------
  !  programme de simulation numerique d'un plasma electrostatique 2D
  !  modelise par les equations de Vlasov-Poisson
  !-------------------------------------------------------------------
  use used_precision  
  use geometry_module
  use diagnostiques_module
  use poisson2dpp_module
  use vlasov2d_module
  use splinepx_class
  use splinepy_class
  use vp2dinit
  use Clock
!#ifdef _MPI
!  use Module_MPI, my_num=>numero_processeur, num_threads=>nombre_processeurs
!#endif
  implicit none

  type (geometry)    :: geomx ! geometrie dans l'espace physique
  type (geometry)    :: geomv ! geometrie dans l'espace des vitesses
  type (poisson2dpp) :: poiss2dpp !champ electrique
  type (vlasov2d)    :: vlas2d ! vlasov
  type (splinepx)    :: splx       ! vlasov1d
  type (splinepy)    :: sply       ! vlasov1d


  real(wp), dimension(:,:,:,:), pointer :: f ! fonc de distribution
  real(wp), dimension(:,:),     pointer :: rho  ! densite de charge
  real(wp), dimension(:,:),     pointer :: ex,ey ! champ electrique

  ! donnees du probleme
  integer      :: nbiter   ! nombre d'iterations en temps
  real(wp)     :: dt       ! pas de temps
  integer      :: fdiag, fthdiag    ! frequences des diagnostiques

  integer      :: iflag    ! indicateur d'erreur
  integer      :: iter,i,j,iv,jv ! variables de boucles       
  character(2) :: ichar   
#ifdef _OPENMP
  integer omp_get_thread_num,omp_num_threads, my_num
#endif
  integer  :: jstartx, jendx, jstartv, jendv
  real(wp) :: nrj,durat,vtime(0:12),cumultime(0:12)
  integer  :: g_a,g_b,l_a,l_b

  ! initialisation global
#ifdef _OPENMP
  my_num = 0
  print*,'OPENMP version of  slv2d running on ',omp_num_threads, ' processors'
#endif
#ifdef _MPI
  call initialise_moduleMPI
  if (my_num.eq.0) then
     print*,'MPI Version of slv2d running on ',num_threads, ' processors'
  end if
#endif

  !call clck_temps(l_a)
  call initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)
  !call clck_temps(l_b)
  !call clck_diff(l_a,l_b,durat)
  !vtime(1)=durat
  
  if (my_num.eq.0) then
     ! write some run data
     write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
     write(*,"(2(i3,1x),6(g13.3,1x))") geomx%nx, geomx%ny, geomx%x0, &
          geomx%x0+(geomx%nx)*geomx%dx, &
          geomx%y0, geomx%y0+(geomx%ny)*geomx%dy, geomx%dx, geomx%dy   
     write(*,*) 'velocity space: nvx, nvy, vx0, vx1, vy0, vy1, dvx, dvy'
     write(*,"(2(i3,1x),6(g13.3,1x))") geomv%nx, geomv%ny, geomv%x0, &
          geomv%x0+(geomv%nx-1)*geomv%dx, &
          geomv%y0, geomv%y0+(geomv%ny-1)*geomv%dy, geomv%dx, geomv%dy
     write(*,*) 'dt,nbiter,fdiag,fthdiag'
     write(*,"(g13.3,1x,3i3)") dt,nbiter,fdiag,fthdiag
  endif

  !call clck_temps(l_a)
#ifndef _MPI
! allocation dynamique des variables globales
  allocate(f(geomx%nx,geomx%ny,geomv%nx,geomv%ny),stat=iflag)
  if (iflag.ne.0) stop 'erreur dans l allocation de f'
! on alloue une dimension nx+2 a rho,ex, ey pour la FFT
  allocate(rho(geomx%nx+2,geomx%ny),stat=iflag)
  if (iflag.ne.0) stop 'erreur dans l allocation de rho'
  allocate(ex(geomx%nx+2,geomx%ny),stat=iflag)
  if (iflag.ne.0) stop 'erreur dans l allocation de ex'
  allocate(ey(geomx%nx+2,geomx%ny),stat=iflag)
  if (iflag.ne.0) stop 'erreur dans l allocation de ey'

!$OMP PARALLEL                      &
!$OMP+ PRIVATE(jstartx,jendx,jstartv,jendv,iter,ichar,poiss2dpp,vlas2d)    &
!$OMP+ SHARED(f,rho,ex,ey,geomx,geomv,dt,nbiter,fdiag)
  my_num = omp_get_thread_num()
#endif
  !call clck_temps(l_b)
  !call clck_diff(l_a,l_b,durat)
  !vtime(2)=durat

  !initialisation en en parallele 
  !call clck_temps(l_a)
  call initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
       f,rho,ex,ey,vlas2d,poiss2dpp,splx,sply)
  ! ecriture des resultats par le processeur 0 a l'instant initial
#ifdef _MPI
  call mpi_barrier(MPI_COMM_WORLD,iflag)
#else
  !$OMP BARRIER 
#endif
  !call clck_temps(l_b)
  !call clck_diff(l_a,l_b,durat)
  !vtime(3)=durat

  !call clck_temps(l_a)
  iter = 0
  call diagnostiques(f,rho,ex,ey,geomx,geomv,jstartx,jendx,jstartv,jendv,iter)
 
#ifdef _MPI
  call mpi_barrier(MPI_COMM_WORLD,iflag)
#else
  !$OMP BARRIER 
#endif
  !call clck_temps(l_b)
  !call clck_diff(l_a,l_b,durat)
  !vtime(4)=durat

!!$     write(*,'(A,I3,I5,10(A,F8.3))') '[', my_num,&
!!$          iter, "] INIT(sec) t1:",vtime(1)," t2:",vtime(2)," t3:",vtime(3),&
!!$          " t4:",vtime(4)
		

  ! advection d'un demi-pas de temps en espace
  call advection_x(vlas2d,f,.5*dt)

  ! boucle en temps
  do iter=1,nbiter

     !call clck_temps(g_a)
     !vtime=0

     !call clck_temps(l_a)
     ! transposition de la fonction de distribution
     call transposexv(vlas2d,f)
     !call clck_temps(l_b)
     !call clck_diff(l_a,l_b,durat)
     !vtime(1)=durat

     !call clck_temps(l_a)
     ! calcul de la densite de charge
     call densite_charge(vlas2d,rho)

     ! calcul du champ electrique
     call solve(poiss2dpp,ex,ey,rho,nrj)
     !call clck_temps(l_b)
     !call clck_diff(l_a,l_b,durat)
     !vtime(2)=durat

     !call clck_temps(l_a)
     ! advection d'un pas de temps en vitesse
     call advection_v(vlas2d,ex,ey,dt)
     !call clck_temps(l_b)
     !call clck_diff(l_a,l_b,durat)
     !vtime(3)=durat

     !call clck_temps(l_a)
     ! transposition de la fonction de distribution
     call transposevx(vlas2d,f)
     !call clck_temps(l_b)
     !call clck_diff(l_a,l_b,durat)
     !vtime(4)=durat

     if (mod(iter,fdiag).eq.0) then !2 demi advections si on fait un diag.
        ! advection d'un demi-pas de temps en espace
        !call clck_temps(l_a)
        call advection_x(vlas2d,f,.5*dt)
        ! ecriture des resultats par le processeur 0
        !if (my_num.eq.0) then
           !print*, 'iteration ', iter
        !endif
        call diagnostiques(f,rho,ex,ey,geomx,geomv, &
             jstartx,jendx,jstartv,jendv,iter/fdiag)
        !!call clck_temps(l_b)
        !call clck_diff(l_a,l_b,durat)
        !vtime(6)=durat
        if (mod(iter,fthdiag).eq.0) then
           !call clck_temps(l_a)
           call thdiag(vlas2d,f,nrj,iter*dt)    
           !call clck_temps(l_b)
           !call clck_diff(l_a,l_b,durat)
           !vtime(5)=durat
        end if
        ! advection d'un demi-pas de temps en espace
        !call clck_temps(l_a)
        call advection_x(vlas2d,f,.5*dt)
        !call clck_temps(l_b)
        !call clck_diff(l_a,l_b,durat)
        !vtime(6)=vtime(6)+durat
     else   ! dans le cas ou l'on ne fait pas de diagnostique
        ! advection d'un pas de temps en espace
        !call clck_temps(l_a)
        call advection_x(vlas2d,f,dt)
        !call clck_temps(l_b)
        !call clck_diff(l_a,l_b,durat)
        !vtime(6)=durat
     end if

     !call clck_temps(g_b)
     !call clck_diff(g_a,g_b,durat)
     !vtime(0)=durat


!     write(*,'(A,I3,I5,10(A,F8.3))') '[', my_num,&
!          iter, "] TEMPS(sec) champ:",vtime(2)," transpositions:",vtime(1)+vtime(4)," velo:",vtime(3),&
!          " spac:",vtime(6)," diag:", vtime(5)," tot:",vtime(0)," sum:",SUM(vtime(1:6))
     
  end do
!$OMP END PARALLEL
  call termine_moduleMPI
end program VP2D
