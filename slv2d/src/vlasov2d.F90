module Vlasov2d_module

  use used_precision
  use splinenn_class
  use splinepp_class
  use geometry_module
  use diagnostiques_module
  use Clock

!#ifdef _MPI
!  use MODULE_MPI, my_num=>numero_processeur, num_threads=>nombre_processeurs
!#endif

  implicit none
  private
  public :: new, dealloc,advection_x, advection_v,&
            densite_charge,transposexv,transposevx,thdiag,densite_courant
  type, public :: vlasov2d
     real(wp), dimension(:,:,:,:), pointer :: ft
     type (splinepp) :: splinex ! spline periodique pour X
     type (splinenn) :: splinev ! spline naturel pour V
     type (geometry) :: geomx, geomv
     logical :: transpose       ! permet de definir si f ou ft est derniere 
                                ! fonction de distribution mise a jour
     integer :: jstartx, jendx
     integer :: jstartv, jendv  ! definition de la bande de calcul
  end type vlasov2d
  ! variables globales 
  real(wp),dimension(0:9) :: diag
  real(wp),dimension(13) :: aux
  real(wp) :: vbeam  ! beam velocity
  real(wp), dimension(:,:),allocatable :: P_x, P_y
  interface new
     module procedure new_vlasov2d
  end interface
  interface dealloc
     module procedure dealloc_vlasov2d
  end interface

contains

  subroutine new_vlasov2d(this,geomx,geomv,iflag, jstartx, jendx, &
       jstartv, jendv,vz)
    type(vlasov2d),intent(out) :: this
    type(geometry),intent(in)  :: geomx, geomv
    real(wp), optional  :: vz
    integer, intent(out) :: iflag
    integer, intent(in), optional ::  jstartx, jendx, jstartv, jendv

    integer :: ierr ! indicateur d'erreur

    if (.not.(present(vz))) then
       vbeam = 1.0
    else
    ! beam velocity
       vbeam = vz
    end if
    ! indicateur d'erreur
    iflag = 0

    ! on commence par utiliser la fonction f(x,y,vx,vy)
    this%transpose=.false.
    ! definition des bandes de calcul (en n'oubliant pas le cas sequentiel)
    if (.not.(present(jstartx))) then
       this%jstartx = 1
    else
       this%jstartx = jstartx
    end if
    if (.not.(present(jendx))) then
       this%jendx = geomx%ny
    else
       this%jendx = jendx
    end if
    if (.not.(present(jstartv))) then 
       this%jstartv = 1
    else
       this%jstartv = jstartv
    end if
    if (.not.(present(jendv))) then
       this%jendv = geomv%nx
    else
       this%jendv = jendv
    end if
   
    ! initialisation de la geometrie
    this%geomx=geomx
    this%geomv=geomv
    ! initialisation des splines de l'espace des vitesses
    call new(this%splinev,geomv,iflag)
    ! initialisation des splines de l'espace physique
    call new(this%splinex,geomx,iflag)  

    ! allocation memoire
    allocate(this%ft(geomv%nx,geomv%ny,geomx%nx, &
         this%jstartx:this%jendx),stat=ierr)
    if (ierr.ne.0) iflag=20
  end subroutine new_vlasov2d

  subroutine dealloc_vlasov2d(this)
    type(vlasov2d),intent(out) :: this
    deallocate(this%ft)
  end subroutine dealloc_vlasov2d


  subroutine advection_x(this,f,dt)
    !-----------------------------------------------
    ! fait une advection en x sur un pas de temps dt
    !-----------------------------------------------
    type(vlasov2d),intent(inout) :: this
#ifdef _MPI
    real(wp), dimension(:,:,:,this%jstartv:) :: f
#else
    real(wp), dimension(:,:,:,:) :: f
#endif
    real(wp), intent(in) :: dt
    ! variables locales
    real(wp) :: depx, depy   ! deplacement par rapport au maillage
    real(wp) :: vx, vy       ! vitesse du point courant
    integer :: iv, jv ! indices de boucle

    ! verifier que la transposition est a jours
    if (this%transpose) stop 'advection_x: on travaille sur f et pas ft'
    do jv=this%jstartv,this%jendv
       vy = this%geomv%y0+(jv-1)*this%geomv%dy
       depy = vy*dt
       do iv=1,this%geomv%nx
          vx = this%geomv%x0+(iv-1)*this%geomv%dx
          depx = vx*dt

!          call interpole(this%splinex,f(:,:,iv,jv),depx,depy,(jv .eq. 3) .and. (iv .eq. 3))
          call interpole(this%splinex,f(:,:,iv,jv),depx,depy,jv==0)
       end do
    end do

!print*,'maxvaladx ',maxval(f),minval(f)
  end subroutine advection_x

  subroutine advection_v(this,fx,fy,dt,bz)
    !-----------------------------------------------
    ! fait une advection en v sur un pas de temps dt
    !-----------------------------------------------
    type(vlasov2d),intent(inout) :: this
#ifdef _MPI
!    real(wp), dimension(:,this%jstartx:),intent(in)  :: fx, fy
    real(wp), dimension(:,:), intent(in) :: fx, fy
#else
    real(wp), dimension(:,:),intent(in)  :: fx,fy
#endif
    real(wp), dimension(:,:), optional, intent(in) :: bz
    real(wp), intent(in) :: dt
    ! variables locales
    real(wp) :: depvx, depvy   ! deplacement par rapport au maillage
    integer :: i, j ,iv, jv, im1! indices de boucle
    real(wp) :: ctheta, stheta, px, py

    ! verifier que la transposition est a jour
    if (.not.(this%transpose)) &
         stop 'advection_v: on travaille sur ft et pas f'
    !print*,'maxvaladvavant ',maxval(this%ft),minval(this%ft)
    !print*,'maxvalchamp ',maxval(fx),minval(fy)

    if (present(bz)) then

       if (.not. allocated(P_x)) allocate(P_x(this%geomv%nx,this%geomv%ny))
       if (.not. allocated(P_y)) allocate(P_y(this%geomv%nx,this%geomv%ny))

       do j=this%jstartx,this%jendx
          do i=1,this%geomx%nx
             ctheta = cos(bz(i,j)*dt)
	     stheta = sin(bz(i,j)*dt)
             depvx  = -0.5*dt*fx(i,j)
             depvy  = -0.5*dt*fy(i,j)
             do jv=1,this%geomv%ny-1
                py = this%geomv%y0+(jv-1)*this%geomv%dy
                do iv=1,this%geomv%nx-1 
                   px = this%geomv%x0+(iv-1)*this%geomv%dx
                   P_x(iv,jv) = depvx+(px+depvx)*ctheta-(py+depvy)*stheta
                   P_y(iv,jv) = depvy+(px+depvx)*stheta+(py+depvy)*ctheta
                end do
             end do
             call interpole(this%splinev,this%ft(:,:,i,j),this%ft(:,:,i,j), &
                  P_x,P_y)
          end do
       end do

    else

       do j=this%jstartx,this%jendx
          do i=1,this%geomx%nx
	    im1=mod(i-1+this%geomx%nx,this%geomx%nx)
!             depvx = fx(i,j)*dt
!             depvy = fy(i,j)*dt
             depvx =  fx(i,j)*dt;depvy=0._wp
             call interpole(this%splinev,this%ft(:,:,i,j),depvx,depvy,j==0)
             !call interpole(this%splinev,this%ft(:,:,i,j),depvx,depvy,(j .eq. 3) .and. (i .eq. 3))
          end do
       end do

     end if
       
  end subroutine advection_v

  subroutine densite_charge(this, rho)
    !------------------------------------------------
    ! calcule la densite de charge rho a partir de ft
    ! en fait le moment d'ordre 0 de f. Les constantes
    ! ne sont pas introduites
    !------------------------------------------------
    ! Poisson n'est pas parallele on transmet donc rho
    ! a tous les processeurs

    type(vlasov2d),intent(inout) :: this
#ifdef _MPI
    !    real(wp), dimension(:,this%jstartx:), intent(out)  :: rho
    integer :: mpierror
    real(wp), dimension(:,:), intent(out)  :: rho
#else
    real(wp), dimension(:,:), intent(out)  :: rho
#endif
    real(wp), dimension(this%geomx%nx,this%geomx%ny) :: locrho
    ! variables locales
    integer :: i,j,iv,jv,c   ! indices de boucles

    ! verifier que la transposition est a jour
    if (.not.(this%transpose)) &
	stop 'densite_charge: on travaille sur ft et pas f'
    !    rho(:,this%jstartx:this%jendx)=0.
    rho(:,:) = 0.
    locrho(:,:) = 0.
    do j=this%jstartx,this%jendx
       do i=1,this%geomx%nx
          do jv=1,this%geomv%ny-1
             do iv=1,this%geomv%nx-1 
                locrho(i,j) = locrho(i,j) + this%geomv%dx*this%geomv%dy* &
                     this%ft(iv,jv,i,j) 
!                print*,'ft ',iv,jv,i,j,rho(i,j),this%ft(iv,jv,i,j)
             end do
          end do
!          print*,'rho ',i,j,rho(i,j),this%jstartx,this%jendx
       end do
    end do
!    print *,'mpireduce',this%jstartx,this%jendx,this%geomx%nx*this%geomx%ny,sum(locrho)*this%geomx%dx*this%geomx%dy
    !print*,'maxvaldc ',maxval(this%ft),minval(this%ft)
    !print*,'maxvalrho ',maxval(rho),minval(rho)
!!$    print*, 'jstart, jend ', this%jstartx,this%jendx
!!$    rho(:,this%jstartx:this%jendx)= 10.0+my_num
!!$    do i=1,this%geomx%nx
!!$       do j=1,this%geomx%ny
!!$          print*, my_num, 'avant dens ', i,j,rho(i,j)
!!$       end do
!!$    end do
    call mpi_barrier(MPI_COMM_WORLD,i)
    c=this%geomx%nx*this%geomx%ny
    call mpi_allreduce(locrho,rho,c, &
         MPI_realtype,MPI_SUM,MPI_COMM_WORLD,mpierror)
!    print *,'verif mpireduce',sum(rho)*this%geomx%dx*this%geomx%dy
!!$    do i=1,this%geomx%nx
!!$       do j=1,this%geomx%ny
!!$          print*, my_num, 'dens ', i,j,rho(i,j)
!!$       end do
!!$    end do
  end subroutine densite_charge

  subroutine densite_courant(this, jx, jy)
    !------------------------------------------------
    ! calcule la densite de courant jx et jy a partir de ft
    ! en fait le moment d'ordre 0 de f. Les constantes
    ! ne sont pas introduites
    !------------------------------------------------

    type(vlasov2d),intent(inout) :: this
#ifdef _MPI
    !    real(wp), dimension(:,this%jstartx:), intent(out)  :: rho
    integer :: mpierror
    real(wp), dimension(:,:), intent(out)  :: jx, jy
    real(wp) :: vx, vy       ! vitesse du point courant
#else
    real(wp), dimension(:,:), intent(out)  :: jx, jy
#endif
    real(wp), dimension(this%geomx%nx,this%geomx%ny) :: locjx
    real(wp), dimension(this%geomx%nx,this%geomx%ny) :: locjy
    ! variables locales
    integer :: i,j,iv,jv,c   ! indices de boucles

    ! verifier que la transposition est a jour
    if (.not.(this%transpose)) &
	stop 'densite_courant: on travaille sur ft et pas f'
    !    rho(:,this%jstartx:this%jendx)=0.
    jx(:,:) = 0.; jy(:,:) = 0.
    locjx(:,:) = 0.; locjy(:,:) = 0.
    do j=this%jstartx,this%jendx
       do i=1,this%geomx%nx
          do jv=1,this%geomv%ny-1
             vy = this%geomv%y0+(jv-1)*this%geomv%dy
             do iv=1,this%geomv%nx-1 
                vx = this%geomv%x0+(iv-1)*this%geomv%dx
                locjx(i,j) = locjx(i,j) + this%geomv%dx*this%geomv%dy* &
                     this%ft(iv,jv,i,j) * vx
                locjy(i,j) = locjy(i,j) + this%geomv%dx*this%geomv%dy* &
                     this%ft(iv,jv,i,j) * vy
!                print*,'ft ',iv,jv,i,j,rho(i,j),this%ft(iv,jv,i,j)
             end do
          end do
!          print*,'rho ',i,j,rho(i,j),this%jstartx,this%jendx
       end do
    end do
!    print *,'mpireduce',this%jstartx,this%jendx,this%geomx%nx*this%geomx%ny,sum(locrho)*this%geomx%dx*this%geomx%dy
    !print*,'maxvaldc ',maxval(this%ft),minval(this%ft)
    !print*,'maxvalrho ',maxval(rho),minval(rho)
!!$    print*, 'jstart, jend ', this%jstartx,this%jendx
!!$    rho(:,this%jstartx:this%jendx)= 10.0+my_num
!!$    do i=1,this%geomx%nx
!!$       do j=1,this%geomx%ny
!!$          print*, my_num, 'avant dens ', i,j,jx(i,j),jy(i,j)
!!$       end do
!!$    end do
    call mpi_barrier(MPI_COMM_WORLD,i)
    c=this%geomx%nx*this%geomx%ny
    call mpi_allreduce(locjx,jx,c, &
         MPI_realtype,MPI_SUM,MPI_COMM_WORLD,mpierror)
    call mpi_allreduce(locjy,jy,c, &
         MPI_realtype,MPI_SUM,MPI_COMM_WORLD,mpierror)
!!$    print *,'verif mpireduce',sum(jx)*this%geomx%dx*this%geomx%dy
!!$    print *,'verif mpireduce',sum(jy)*this%geomx%dx*this%geomx%dy
!!$    do i=1,this%geomx%nx
!!$       do j=1,this%geomx%ny
!!$          print*, my_num, 'jx,jy ', i,j,jx(i,j),jy(i,j)
!!$       end do
!!$    end do
  end subroutine densite_courant

  subroutine transposexv(this,f)
    !---------------------------------------------
    ! transpose la fonction de distribution
    ! ATTENTION: cet fonction fait intervenir des 
    ! communications entre les processeurs.
    !---------------------------------------------
    type(vlasov2d),intent(inout) :: this
    real(wp), dimension(:,:,:,:),intent(in) :: f
#ifndef _MPI
    ! variables locales
    integer :: i,j,iv,jv   ! indices de boucles
#endif
    
#ifdef _MPI
    integer :: sizexy, sizevxvy

    sizexy = this%geomx%nx * this%geomx%ny
    sizevxvy = this%geomv%nx * this%geomv%ny
    call transpose(f, this%ft, sizexy, sizevxvy, num_threads)
#else
    ! on attend que tous les processeurs aient fait leur calcul
!$OMP BARRIER  
    do jv=1,this%geomv%ny
       do iv=1,this%geomv%nx
          do j=this%jstartx,this%jendx
             do i=1,this%geomx%nx  
                this%ft(iv,jv,i,j) = f(i,j,iv,jv)
             end do
          end do
       end do
    end do
   ! on attend que tous les processeurs aient fait leur transposition
!$OMP BARRIER 
#endif
! a partir de maintenant on travaille sur ft
! 
    this%transpose=.true.
!
  end subroutine transposexv

  subroutine transposevx(this,f)
! transpose la fonction de distribution
    type(vlasov2d),intent(inout) :: this
    real(wp), dimension(:,:,:,:),intent(out) :: f
#ifndef _MPI
    ! variables locales
    integer :: i,j,iv,jv   ! indices de boucles
#endif
#ifdef _MPI
    integer :: sizexy, sizevxvy

    sizexy = this%geomx%nx * this%geomx%ny
    sizevxvy = this%geomv%nx * this%geomv%ny
    call transpose(this%ft,f, sizevxvy, sizexy, num_threads)
#else
    ! on attend que tous les processeurs aient fait leur calcul
!$OMP BARRIER
    do j=this%jstartx,this%jendx
       do i=1,this%geomx%nx
          do jv=1,this%geomv%ny
             do iv=1,this%geomv%nx
                f(i,j,iv,jv) = this%ft(iv,jv,i,j)
             end do
          end do
       end do
    end do
    ! on attend que tous les processeurs aient fait leur transposition
!$OMP BARRIER 
#endif
!
! a partir de maintenant on travaille sur f
    this%transpose=.false.
!
  end subroutine transposevx

  subroutine thdiag(this,f,nrj,t)
    ! time history diagnostics
    type(vlasov2d),intent(inout) :: this
#ifdef _MPI
    real(wp), dimension(:,:,:,this%jstartv:),intent(in) :: f
    integer :: mpierror
#else
    real(wp), dimension(:,:,:,:),intent(in) :: f
    integer :: my_num,mp_get_thread_num
#endif
    real(wp), intent(in) :: t,nrj   ! current time
    ! variables locales
    integer :: i,iv, j,jv   ! indices de boucles
    real :: x, vx, y, vy, dxdv
    real(wp),dimension(7) :: diagloc
    real(wp),dimension(11) :: auxloc

! initialisation des variables globales
#ifdef _OPENMP
    my_num = omp_get_thread_num()
#endif
    if (my_num.eq.0) then
       diag=0.
       aux=0.
    end if

    diagloc = 0._wp
    auxloc  = 0._wp
    do i = 1,this%geomx%nx
       x = this%geomx%x0+(i-1)*this%geomx%dx
       do j = 1,this%geomx%ny
          y= this%geomx%y0+(j-1)*this%geomx%dy
          do iv=1,this%geomv%nx
             vx = this%geomv%x0+(iv-1)*this%geomv%dx
             do jv=this%jstartv,this%jendv
                vy = this%geomv%y0+(jv-1)*this%geomv%dy
                diagloc(2) = diagloc(2) + f(i,j,iv,jv)*f(i,j,iv,jv)

                auxloc(1) = auxloc(1) + f(i,j,iv,jv)         ! avg(f)
                auxloc(2) = auxloc(2) + x*f(i,j,iv,jv)       ! avg(x)
                auxloc(3) = auxloc(3) + vx*f(i,j,iv,jv)      ! avg(vx)
                auxloc(4) = auxloc(4) + x*x*f(i,j,iv,jv)     ! avg(x^2)
                auxloc(5) = auxloc(5) + vx*vx*f(i,j,iv,jv)   ! avg(vx^2)
                auxloc(6) = auxloc(6) + x*vx*f(i,j,iv,jv)    ! avg(x*vx)
                auxloc(7) = auxloc(7) + y*f(i,j,iv,jv)       ! avg(y)
                auxloc(8) = auxloc(8) + vy*f(i,j,iv,jv)      ! avg(vy)
                auxloc(9) = auxloc(9) + y*y*f(i,j,iv,jv)     ! avg(y^2)
                auxloc(10) = auxloc(10) + vy*vy*f(i,j,iv,jv) ! avg(vy^2)
                auxloc(11) = auxloc(11) + y*vy*f(i,j,iv,jv)  ! avg(y*vy)
             end do
          end do
       end do
    end do
    auxloc=auxloc!*this%geomx%dx*this%geomx%dy*this%geomv%dx*this%geomv%dy
#ifdef _MPI
    call mpi_reduce(auxloc,aux,11,MPI_realtype,MPI_SUM,0,  &
         MPI_COMM_WORLD, mpierror)
    call mpi_reduce(diagloc(2),diag(2),1,MPI_realtype,MPI_SUM,0, &
         MPI_COMM_WORLD, mpierror)
#else
!$OMP CRITICAL
diag(2)=diag(2)+diagloc(2)
aux=aux+auxloc
!$OMP END CRITICAL

!$OMP BARRIER
#endif

!if (my_num.eq.0) then
!    ! multiply by dxdv
!    dxdv=this%geomx%dx*this%geomx%dy*this%geomv%dx*this%geomv%dy
!    diag(0)= t*vbeam
!    diag(1)=aux(1)*dxdv
!    diag(2)=diag(2)*dxdv
    ! compute rms emittance in x-xp
!    diag(3)=4*sqrt((aux(4)/aux(1)-(aux(2)/aux(1))**2)*(aux(5)/aux(1) &
!         -(aux(3)/aux(1))**2)-(aux(6)/aux(1))**2)/vbeam
!    diag(4)=sqrt(aux(4)/aux(1))   ! xrms
!    diag(5)=sqrt(aux(5)/aux(1))/vbeam ! x'rms
    ! compute rms emittance in y-yp
!    diag(6)=4*sqrt((aux(9)/aux(1)-(aux(7)/aux(1))**2)*(aux(10)/aux(1) &
!         -(aux(8)/aux(1))**2)-(aux(11)/aux(1))**2)/vbeam
!    diag(7)=sqrt(aux(9)/aux(1))     ! yrms
!    diag(8)=sqrt(aux(10)/aux(1))/vbeam  ! y'rms

!    call time_history("thf","(9(1x,e15.6))",diag(0:8))  
!endif

if (my_num==0) then
   aux(13)=t
   aux(12)=nrj
   write(*,"('time ', g8.3,' test nrj',f10.5)") t, nrj
   call time_history("thf","(13(1x,e15.6))",aux(1:13))
end if

  end subroutine thdiag
end module Vlasov2d_module


