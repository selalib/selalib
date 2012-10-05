program VM2D
!-------------------------------------------------------------------
!  programme de simulation numerique d'un plasma electromagnetique 2D
!  modelise par les equations de Vlasov-Maxwell
!-------------------------------------------------------------------
use used_precision  
use geometry_module
use maxwell2dfdtd_module
use poisson2dpp_module
use diagnostiques_module
use diagnostiquesm_module
use vlasov2d_module
use splinepx_class
use splinepy_class
use vlasov1d_module
use vm2dinit

implicit none

type (geometry)      :: geomx      ! geometrie dans l'espace physique
type (geometry)      :: geomv      ! geometrie dans l'espace des vitesses
type (maxwell2dfdtd) :: maxw2dfdtd ! champ electromagnetique
type (poisson2dpp)   :: poiss2dpp  ! potentiel pour la correction
type (vlasov2d)      :: vlas2d     ! vlasov
type (splinepx)      :: splx       ! vlasov1d
type (splinepy)      :: sply       ! vlasov1d


real(wp), dimension(:,:,:,:), pointer :: f,f1     ! fonc de distribution
real(wp), dimension(:,:),     pointer :: ex,ey ! champ electrique
real(wp), dimension(:,:),     pointer :: ex1,ey1 ! champ electrique
real(wp), dimension(:,:),     pointer :: jx,jy ! courant
real(wp), dimension(:,:),     pointer :: bz,bz1    ! champ magnetique
real(wp), dimension(:,:),     pointer :: rho   ! charge
real(wp), dimension(:,:),     pointer :: Jx1,Jx2,Jy1,Jy2   ! courant partiel 
real(wp), dimension(:,:),     pointer :: div   ! divergence E

! donnees du probleme
integer      :: nbiter         ! nombre d'iterations en temps
real(wp)     :: dt             ! pas de temps
integer      :: fdiag, fthdiag ! frequences des diagnostiques
integer      :: iflag          ! indicateur d'erreur
integer      :: iter,i,j,iv,jv ! variables de boucles       
character(2) :: ichar   

integer  :: jstartx, jendx, jstartv, jendv
real(wp) :: nrj, omega, nd, md

! initialisation global

call initialise_moduleMPI
if (my_num.eq.0) then
   print*,'MPI Version of slv2d running on ',num_threads, ' processors'
end if

call initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)
  
if (my_num.eq.0) then
   ! write some run data
   write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
   write(*,"(2(i3,1x),6(g13.3,1x))") 			&
	geomx%nx, geomx%ny, geomx%x0, 			&
        geomx%x0+(geomx%nx-1)*geomx%dx, 		&
        geomx%y0, geomx%y0+(geomx%ny-1)*geomx%dy, 	&
	geomx%dx, geomx%dy   
   write(*,*) 'velocity space: nvx, nvy, vx0, vx1, vy0, vy1, dvx, dvy'
   write(*,"(2(i3,1x),6(g13.3,1x))") 			&
	geomv%nx, geomv%ny, geomv%x0, 			&
        geomv%x0+(geomv%nx-1)*geomv%dx, 		&
        geomv%y0, geomv%y0+(geomv%ny-1)*geomv%dy, 	&
	geomv%dx, geomv%dy
   write(*,*) 'dt,nbiter,fdiag,fthdiag'
   write(*,"(g13.3,1x,3i3)") dt,nbiter,fdiag,fthdiag
endif

call initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
       		f,f1,rho,ex,ey,ex1,ey1,bz,bz1,jx,jy,    &
                vlas2d,maxw2dfdtd,poiss2dpp,splx,sply)

! ecriture des resultats par le processeur 0 a l'instant initial
call mpi_barrier(MPI_COMM_WORLD,iflag)

iter = 0
call diagnostiquesm(f,rho,ex,ey,bz,jx,jy,geomx,geomv,	&
		    jstartx,jendx,jstartv,jendv,iter)

call mpi_barrier(MPI_COMM_WORLD,iflag)

call transposevx(vlas2d,f)

allocate(div(geomx%nx,geomx%ny))
allocate(Jx1(geomx%nx,geomx%ny),Jx2(geomx%nx,geomx%ny))
allocate(Jy1(geomx%nx,geomx%ny),Jy2(geomx%nx,geomx%ny))



do iter=1,nbiter	!Loop over time

   ! advection d'un demi-pas de temps en espace
   !call advection_x(vlas2d,f,.5*dt)
   call advection1d_x(splx,f,0.5_wp*dt,Jx1)
   call advection1d_y(sply,f,0.5_wp*dt,Jy1)

!   print *,'fin advecx'

   !Recopie de f^(n,2) dans f1
   f1=f;ex1=ex;ey1=ey;bz1=bz

   ! transposition de la fonction de distribution
   call transposexv(vlas2d,f)

   !###########################
   !Phase prediction a la VALIS
   !########################### 

   ! calcul de la densite de courant et de charge
   ! Les deux operations doivent etre faites dans cet ordre
   ! Chercher la raison ...
   call densite_charge(vlas2d,rho)
   call densite_courant(vlas2d,jx,jy)
!   print *,'calcul jx normal',sum(jx),sum(jy)
     
   ! calcul du champ magnetique Bz(k+1/2)
   if (iter == 1) then
      call solve_faraday(maxw2dfdtd,ex,ey,bz,0.5_wp*dt)
   else
      call solve_faraday(maxw2dfdtd,ex,ey,bz,dt)
   end if

   call cl_periodiques(maxw2dfdtd,ex,ey,bz,jx,jy,dt)
   !call silver_muller(maxw2dfdtd,ex,ey,bz,jx,jy,dt)

   call solve_ampere(maxw2dfdtd,ex,ey,bz,jx,jy,nrj,dt)

   ! advection d'un pas de temps en vitesse
   call advection_v(vlas2d,ex,ey,dt,bz)
!   call advection_v(vlas2d,ex,ey,dt)

   ! transposition de la fonction de distribution
   call transposevx(vlas2d,f)

   ! advection d'un demi-pas de temps en espace     
!      call advection_x(vlas2d,f,0.5_wp*dt)
   call advection1d_x(splx,f,0.5_wp*dt,Jx2)
   call advection1d_y(sply,f,0.5_wp*dt,Jy2)

   !Calcul du J qui conserve la charge
!   jx=geomx%dx*(Jx1+Jx2)/dt;jy=geomx%dy*(Jy1+Jy2)/dt

!   print *,'fin prediction',sum(jx),sum(jy)
!   print *,'j partiel',sum(Jx1),sum(Jx2),sum(Jy1),sum(Jy2)
!stop
   f=f1;ex=ex1;ey=ey1;

   !################
   !Phase correction
   !################
   call cl_periodiques(maxw2dfdtd,ex,ey,bz,jx,jy,dt)
   !call silver_muller(maxw2dfdtd,ex,ey,bz,jx,jy,dt)
   call solve_ampere(maxw2dfdtd,ex,ey,bz,jx,jy,nrj,dt)
!   print *,'ap ampere',sum(ex),sum(ey)

   ! transposition de la fonction de distribution
   call transposexv(vlas2d,f)

   ! advection d'un pas de temps en vitesse
   call advection_v(vlas2d,ex,ey,dt,bz)

   ! transposition de la fonction de distribution
   call transposevx(vlas2d,f)

   ! advection d'un demi-pas de temps en espace     
!      call advection_x(vlas2d,f,0.5_wp*dt)
   call advection1d_x(splx,f,0.5_wp*dt,Jx1)
   call advection1d_y(sply,f,0.5_wp*dt,Jy1)

   !print *,'fin iter'

   !#############
   !Fin iteration
   !#############

   !Diagnostic
   if (mod(iter,fdiag).eq.0) then 
      ! ecriture des resultats par le processeur 0
      call diagnostiquesm(f,rho,ex,ey,bz,jx,jy,geomx,geomv, &
                          jstartx,jendx,jstartv,jendv,iter/fdiag)
   endif

   !write(*,*,advance='no') iter
   if (mod(iter,fthdiag).eq.0) then 
      call thdiag(vlas2d,f,nrj,iter*dt)
   endif

end do

call termine_moduleMPI

end program VM2D
