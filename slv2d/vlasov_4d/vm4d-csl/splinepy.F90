module splinepy_class
  use used_precision
  use geometry_module!1d_module
  implicit none
!  private
  public :: initialize, interpole
  type, public :: splinepy
     type (geometry) :: geomx,geomv
     real(wp) :: a1y, a2y, a3y, a4y ! coef de la matrice 2x2 per.
     real(wp), dimension(:), pointer :: ayd ! termes diagonaux de ax
     real(wp), dimension(:), pointer :: ayod ! termes sous-diag. de ax
     real(wp), dimension(:,:), pointer :: aym1gamma
     real(wp), dimension(:), pointer :: dper,lper,mper !LU coefficient of the circulant matrix
     real(wp), dimension(:), pointer :: coef,tmpcoef! coefficients des splines
     logical :: transpose       ! permet de definir si f ou ft est derniere 
                                ! fonction de distribution mise a jour
     integer :: jstartx, jendx
     integer :: jstartv, jendv  ! definition de la bande de calcul
  end type splinepy
  interface initialize
     module procedure new_splinepy
  end interface
  interface interpole
     module procedure interpole_splinepdep
  end interface
contains
  subroutine new_splinepy(this,geom,geomv,iflag,jstartx,jendx,jstartv,jendv,vz)
    type(splinepy), intent(out) :: this
    type(geometry), intent(in) :: geom,geomv  ! geometry of problem
    real(wp), optional  :: vz
    integer, intent(out)       :: iflag ! error flag
    integer, intent(in), optional ::  jstartx,jendx,jstartv,jendv

    ! local variables
    integer :: err ! error flag
    integer :: ny! dimensions
    integer :: nyp1, nyp2,nyp3,Nbdr
    real(wp) :: aa1y, aa2y
    integer i,j, info

    ! initialisation des variables locales
    iflag = 0
    ny = geom%ny  !ES  +1
    nyp1=ny+1 ; nyp2=ny+2 ; nyp3=ny+3
    Nbdr=10

   ! on commence par utiliser la fonction f(x,y,vx,vy)
    this%transpose=.false.
    ! definition des bandes de calcul (en n'oubliant pas le cas sequentiel)
    if (.not.(present(jstartx))) then
       this%jstartx = 1
    else
       this%jstartx = jstartx
    end if
    if (.not.(present(jendx))) then
       this%jendx = geom%ny
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

    ! initialisation de geom
    this%geomx = geom
    this%geomv = geomv

    ! memory allocation
    allocate(this%aym1gamma(nyp1,2), stat=err)
    if (err.ne.0) then
       iflag = 10
       return
    end if
    allocate(this%ayd(nyp1), stat=err)
    if (err.ne.0) then
       iflag = 20
       return
    end if
    allocate(this%ayod(ny), stat=err)
    if (err.ne.0) then
       iflag = 30
       return
    end if
    allocate(this%dper(ny), stat=err)
    if (err.ne.0) then
       iflag = 40
       return
    end if
    allocate(this%lper(ny-1), stat=err)
    if (err.ne.0) then
       iflag = 41
       return
    end if
    allocate(this%mper(ny), stat=err)
    if (err.ne.0) then
       iflag = 42
       return
    end if
    allocate(this%coef(nyp3), stat=err)
    if (err.ne.0) then
       iflag = 50
       return
    end if
    allocate(this%tmpcoef(ny+2*Nbdr), stat=err)
    if (err.ne.0) then
       iflag = 51
       return
    end if


    ! factorisation de la  matrice  Ay 
    this%ayd = 4.0_wp
    this%ayod = 1.0_wp
    call dpttrf(nyp1,this%ayd,this%ayod,err)
    if (err.ne.0) then
       iflag = 70
       return
    end if
    ! compute Ay-1.gamma
    this%aym1gamma = 0.0_wp
    this%aym1gamma(1,2) = 1.0_wp
    this%aym1gamma(nyp1,1) = 1.0_wp
    call dpttrs(nyp1,2,this%ayd,this%ayod,this%aym1gamma,nyp1,err)
    if (err.ne.0) then
       iflag = 90
       return
    end if

    aa1y=3.0_wp/geom%dy
    aa2y=6.0_wp/(geom%dy*geom%dy)
    ! assemblage de la matrice 2x2 pour la spline periodique
       this%a1y = -aa1y*(1._wp + this%aym1gamma(2,1)+this%aym1gamma(ny,1))
       this%a2y = -aa1y*(1._wp + this%aym1gamma(2,2)+this%aym1gamma(ny,2))
       this%a3y = aa2y*(-1._wp + 2._wp*this%aym1gamma(1,1) - this%aym1gamma(2,1) &
            + this%aym1gamma(ny,1) - 2._wp*this%aym1gamma(nyp1,1))
       this%a4y = aa2y*(1._wp + 2._wp*this%aym1gamma(1,2) - this%aym1gamma(2,2) &
            + this%aym1gamma(ny,2) - 2._wp*this%aym1gamma(nyp1,2))

    this%dper(1)=4._wp;this%mper(1)=0.25_wp
    do i=1,ny-1
       this%lper(i)=1._wp/this%dper(i)
       this%dper(i+1)=4._wp-this%lper(i)
       this%mper(i+1)=-this%mper(i)/this%dper(i+1)
    enddo
    this%dper(ny)=this%dper(ny)-this%lper(ny-1)+2._wp*this%mper(ny-1)
    do i=1,ny
       this%dper(i)=1._wp/this%dper(i)
    enddo


  end subroutine new_splinepy


  subroutine interpole_splinepdep(this,f,depy,G,meth) 
    !----------------------------------------------------------------
    ! interpolation par spline periodique.
    ! Les points d'interpolation sont definis grace a depx 
    ! qui definissent le deplacement par rapport au maillage.
    !----------------------------------------------------------------
    type(splinepy), intent(inout) :: this
    ! f contient les valeurs de la fonction de distribution
    real(wp), dimension(:), intent(inout) :: f
    real(wp), dimension(:), intent(inout) :: G
    ! dans depx on trouve les deplacements par rapport au maillage
    ! des points dans les quels on veut evaluer la spline.
    real(wp), intent(in) :: depy
    integer, intent(in) :: meth
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    real(wp), dimension(1:this%geomx%ny) :: tmpf
    integer :: ierr,i,ny

    if (meth==0) then 
       !BSL formulation
       call per_y(this,f,ierr)
       call evaldep(this,depy,f)     
    endif

    if (meth>=1) then 
       !methodes conservatives 

       !construction de G, primitive de f
       G=0._wp
       ny=this%geomx%ny
       
       G(1)=f(1)
       do i=2,ny
          G(i)=G(i-1)+f(i)
       enddo

       if (meth==1) then 
          G(1) = G(1)+(1._wp/6._wp)*G(ny)
          G(ny)= (5._wp/6._wp)*G(ny)
          
          call per_primy(this,G,ierr)
          
          if (ierr.ne.0) then
             iflag = 10
             return
          end if

          G(ny)=(6._wp/5._wp)*G(ny)
          G(1)=G(1)-(1._wp/6._wp)*G(ny)

          tmpf(1:ny)=f(1:ny)

          call evaldepcons(this,depy,f,G)     
          !VF formulation
          G(1:ny)=G(1:ny)*this%geomx%dy/depy
          do i=2,ny
             f(i)=tmpf(i)-depy*(G(i)-G(i-1))/this%geomx%dy
          enddo
          f(1)=tmpf(1)-depy*(G(1)-G(ny))/this%geomx%dy
          
       endif
       
       if (meth>=2) then 
          tmpf(1:ny)=f(1:ny)
          call evaldepconshermite(this,depy,f,G,meth)     
          
          !VF formulation
          G(1:ny)=G(1:ny)*this%geomx%dy/depy
          do i=2,ny
             f(i)=tmpf(i)-depy*(G(i)-G(i-1))/this%geomx%dy
          enddo
          f(1)=tmpf(1)-depy*(G(1)-G(ny))/this%geomx%dy
       endif
    endif
          
  end subroutine interpole_splinepdep

  !
  ! calcul des splines periodiques
  !
  subroutine per_y(this,gtau,iflag)
    type(splinepy), intent(inout) :: this     ! objet de type splinep
    real(wp), dimension(:), intent(in) :: gtau ! valeur de la fonction 
    ! aux points du maillage
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i,j ! indices de boucle
    integer ny, nyp2, nyp3  !  ny+2, ny+3
    integer :: err ! error flag

    real(wp) :: aym1f(this%geomx%ny+1)  !ES +1-> +2
    real(wp) :: det, gamma1, gamma2

    ! initialisations
    iflag =0
    ny=this%geomx%ny  !ES +1
    nyp2=ny+2
    nyp3=ny+3
    det=this%a1y*this%a4y - this%a2y*this%a3y

    ! Calcul de Ay^-1 f
    ! assemblage du membre de droite pour le calcul de Ay^-1 f
    do i=1,ny
       aym1f(i) = 6._wp*gtau(i)
    end do
    aym1f(ny+1) = 6._wp*gtau(1) ! Calcul par periodicite   
    ! on calcule Ay^-1f
    call dpttrs(ny+1,1,this%ayd,this%ayod,aym1f,ny+1,err)
    if (err.ne.0) then
       iflag = 10
       return
    end if
    ! assemblage du second membre du systeme 2x2 
 
    gamma1 =  - (3.0_wp/this%geomx%dy)*(aym1f(2) +aym1f(ny))
    gamma2 =  (6.0_wp/(this%geomx%dy)**2._wp)*(2._wp*aym1f(1) - aym1f(2) &
         + aym1f(ny) - 2._wp*aym1f(ny+1))
    
    this%coef(nyp3)= (gamma1*this%a4y - gamma2*this%a2y)/det
    this%coef(1)= (gamma2*this%a1y - gamma1*this%a3y)/det
    do  i=2,nyp2
       this%coef(i)= aym1f(i-1) &
            - this%aym1gamma(i-1,1)*this%coef(nyp3) &
            - this%aym1gamma(i-1,2)*this%coef(1)
    end do
  end subroutine per_y
  

  !
  ! calcul des splines periodiques par inversion de la matrice circulante
  !
  subroutine per_primy(this,gtau,iflag)
    type(splinepy), intent(inout) :: this     ! objet de type splinep
    real(wp), dimension(:), intent(in) :: gtau ! valeur de la fonction 
    ! aux points du maillage
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i! indices de boucle
    integer ny, nyp2, nyp3  !  ny+2, ny+3

    real(wp) :: aym1f(this%geomx%ny+1)  !ES +1-> +2
    real(wp) :: M

    ! initialisations
    iflag =0
    ny=this%geomx%ny  !ES +1

    ! Calcul de Ax^-1 f
    ! assemblage du membre de droite pour le calcul de Ax^-1 f
    M=(6._wp/5._wp)*gtau(ny)

    aym1f=0._wp
    do i=1,ny
       aym1f(i) = 6._wp*gtau(i)
    end do
    do i=2,ny
       aym1f(i)=aym1f(i)-aym1f(i-1)*this%lper(i-1)
    enddo
    do i=1,ny-1
       aym1f(ny)=aym1f(ny)-this%mper(i)*aym1f(i)
    enddo
    aym1f(ny)=aym1f(ny)*this%dper(ny)
    aym1f(ny-1)=this%dper(ny-1)*(aym1f(ny-1)-(1._wp-this%mper(ny-2))*aym1f(ny))
    do i=ny-2,2,-1
       aym1f(i)=this%dper(i)*(aym1f(i)-aym1f(i+1)+this%mper(i-1)*aym1f(ny))
    enddo
    aym1f(1)=this%dper(1)*(aym1f(1)-aym1f(2)-aym1f(ny))

    do  i=1,ny
       this%coef(i)= aym1f(i) 
    end do
    !calcul des coef aux bords par "periodicite"
    !coef(ny+1)
    this%coef(ny+1)=this%coef(1)+M !gtau(ny)
    !coef(0)
    this%coef(ny+2)=this%coef(ny)-M !gtau(ny)
    !masse totale
    this%coef(ny+3)=M !gtau(ny)

  end subroutine per_primy


  subroutine evaldep(this,alphay,fout)
    type(splinepy) :: this
    ! deplacement par rapport aux maillages des points ou la spline est evaluee
    real(wp), intent(in) :: alphay  !valeur du deplacement
    ! fout(i) contient en sortie la valeur de la spline au 
    ! point (xi-alphax), les (xi) etant les points du maillage
    real(wp), dimension(:), intent(out) :: fout

    ! variables locales
    real(wp) :: bvaly1,bvaly2,bvaly3,bvaly4,dyy,dyyy,dyyy6,yd1,ydp1
    integer :: intaysdy
    integer i1,i

    ! debut du code
    dyy=this%geomx%dy*this%geomx%dy
    dyyy=dyy*this%geomx%dy
    dyyy6=1._wp/(6._wp*dyyy)

    if (alphay.gt.0) then
       intaysdy=int(-alphay/this%geomx%dy+epsilon)-1  
! intaysdy=int(-alphay/this%geomx%dy)-1 
!       intaysdy = -1
    else
       intaysdy=int(-alphay/this%geomx%dy)
    end if
    
    yd1=-alphay-intaysdy*this%geomx%dy
    ydp1=this%geomx%dy-yd1

    bvaly1=ydp1*ydp1*ydp1
    bvaly2=dyyy+3._wp*dyy*ydp1+3._wp*this%geomx%dy*ydp1*ydp1-3._wp*ydp1*ydp1*ydp1
    bvaly3=dyyy+3._wp*dyy*yd1+3._wp*this%geomx%dy*yd1*yd1-3._wp*yd1*yd1*yd1
    bvaly4=yd1*yd1*yd1

    do i=1,this%geomx%ny
       i1=mod(this%geomx%ny+i-1+intaysdy,this%geomx%ny)
       fout(i) = dyyy6*(bvaly1*this%coef(i1+1)+bvaly2*this%coef(i1+2) &
               + bvaly3*this%coef(i1+3) + bvaly4*this%coef(i1+4))
    end do
    ! mise a jour de la derniere ligne par periodicite
 !   fout(this%geomx%ny)=fout(1)
  end subroutine evaldep

  subroutine evaldepcons(this,alphay,fout,flux)
    type(splinepy) :: this
    ! deplacement par rapport aux maillages des points ou la spline est evaluee
    real(wp), intent(in) :: alphay  !valeur du deplacement
    ! fout(i) contient en sortie la valeur de la spline au 
    ! point (xi-alphax), les (xi) etant les points du maillage
    real(wp), dimension(:), intent(out) :: fout
    real(wp), dimension(:), intent(inout) :: flux

    ! variables locales
    integer  :: iy,i,ny,Nbdr
    real(wp) :: w(0:3)
    real(wp) :: alf,dy,y,M,dtmp

    Nbdr=10
    ny=this%geomx%ny

    M=this%coef(ny+3)
    
    do i=1,ny
       this%tmpcoef(i+Nbdr) = this%coef(i)
    enddo
    do i=1,Nbdr
       this%tmpcoef(i)=this%tmpcoef(i+ny)-M
    enddo
    do i=1,Nbdr
       this%tmpcoef(i+ny+Nbdr)=this%tmpcoef(i+Nbdr)+M
    enddo
   !dep normalise
    alf=alphay/(this%geomx%y1-this%geomx%y0);dy=1._wp/real(ny,wp)

    do i=1,this%geomx%ny
       y=i*dy-alf;iy=floor(y*ny);y=y*ny-iy
       w(0)=(1._wp/6._wp)*(1._wp-y)*(1._wp-y)*(1._wp-y);
       w(1)=(1._wp/6._wp)+0.5_wp*(1._wp-y)*(-(1._wp-y)*(1._wp-y)+(1._wp-y)+1._wp)
       w(2)=(1._wp/6._wp)+0.5_wp*y*(-y*y+y+1._wp)
       w(3)=(1._wp/6._wp)*y*y*y
       
       fout(i) = w(0)*this%tmpcoef(iy-1+Nbdr) + w(1)*this%tmpcoef(iy+Nbdr) &
            + w(2)*this%tmpcoef(iy+1+Nbdr) + w(3)*this%tmpcoef(iy+2+Nbdr) 
    end do

    do i=1,this%geomx%ny
       flux(i)=flux(i)-fout(i)
    enddo

    dtmp=fout(ny)
    do i=ny,2,-1
       fout(i)=fout(i)-fout(i-1)
    enddo
    fout(1)=fout(1)+M-dtmp
    

    ! mise a jour de la derniere ligne par periodicite
 !   fout(this%geomx%ny)=fout(1)
  end subroutine evaldepcons



  subroutine evaldepconshermite(this,alphay,fout,flux,meth)
    type(splinepy) :: this
    ! deplacement par rapport aux maillages des points ou la spline est evaluee
    real(wp), intent(in) :: alphay  !valeur du deplacement
    ! fout(i) contient en sortie la valeur de la spline au 
    ! point (xi-alphax), les (xi) etant les points du maillage
    real(wp), dimension(:), intent(out) :: fout
    real(wp), dimension(:), intent(inout) :: flux
    integer, intent(in) :: meth

    ! variables locales
    real(wp), dimension(0:this%geomx%ny) :: ftmp
    real(wp), dimension(-1:this%geomx%ny) :: fluxtmp
    real(wp), dimension(0:this%geomx%ny) :: buf
    integer  :: ix,i0,iy,i,ny,Nbdr,i1,i2,i3,im1,im2,im3
    real(wp) :: w(0:3),fbar,df0,df1
    real(wp) :: alf,dy,x,y,M,dtmp,flux0


    Nbdr=10
    ny=this%geomx%ny

    fluxtmp(0:this%geomx%ny-1)=flux(1:this%geomx%ny)
    fluxtmp(this%geomx%ny)=fluxtmp(0)+fluxtmp(this%geomx%ny-1)
    fluxtmp(-1)=0._8

    ftmp(0:this%geomx%ny-1)=fout(1:this%geomx%ny)
    ftmp(this%geomx%ny)=ftmp(0)

    buf(0:this%geomx%ny-1)=fout(1:this%geomx%ny)
    buf(this%geomx%ny)=buf(0)

   !dep normalise
    alf=alphay/(this%geomx%y1-this%geomx%y0);dy=1._wp/real(ny,wp)

    x=-alf
    x=x-real(floor(x),8)
    x=x*real(ny,8)
    ix=floor(x)
    if(ix==ny)then
      x=0._8;ix=0
    endif
    x=x-real(ix,8)
    
    w(0)=x*(1._8-x)*(1._8-x)
    w(1)=x*x*(x-1._8)
    w(2)=x*x*(3._8-2._8*x)
    i0=ix
    i1=mod(ix+1,ny)
    i2=mod(ix+2,ny)
    i3=mod(ix+3,ny)
    im1=mod(ix+ny-1,ny)
    im2=mod(ix+ny-2,ny)
    im3=mod(ix+ny-3,ny)    

    do i=0,this%geomx%ny-1
       fbar=buf(i0)
       if (meth==2) then 
          !lag3
          df0=(5._8/6._8)*buf(i0)-(1._8/6._8)*buf(i1)+(2._8/6._8)*buf(im1)
          df1=(5._8/6._8)*buf(i0)+(2._8/6._8)*buf(i1)-(1._8/6._8)*buf(im1)
       elseif (meth==3) then 
          !ppm1
          df0=(7._8/12._8)*(buf(i0)+buf(im1))-(1._8/12._8)*(buf(i1)+buf(im2))
          df1=(7._8/12._8)*(buf(i1)+buf(i0)) -(1._8/12._8)*(buf(i2)+buf(im1))
       elseif (meth==4) then 
          !ppm2
          df0=(37._8/60._8)*(buf(i0)+buf(im1))-(8._8/60._8)*(buf(i1)+buf(im2))+(1._8/60._8)*(buf(i2)+buf(im3))
          df1=(37._8/60._8)*(buf(i1)+buf(i0)) -(8._8/60._8)*(buf(i2)+buf(im1))+(1._8/60._8)*(buf(i3)+buf(im2))
       endif

       ftmp(i)=w(0)*df0+w(1)*df1+w(2)*fbar;
       im3=im2;im2=im1;im1=i0;i0=i1;i1=i1+1;if(i1>=ny)i1=i1-ny;i2=i2+1;if(i2>=ny)i2=i2-ny;i3=i3+1;if(i3>=ny)i3=i3-ny;
    enddo

    i0=ix;im1=mod(i0-1+ny,ny)
    fbar=ftmp(0);

    do i=0,ny-2
       !schema de base
!      ftmp(i)=buf(i0)+ftmp(i+1)-ftmp(i)  !ok for i=0,ny-2

       !schema volumes finis
      ftmp(i)=buf(i)-(fluxtmp(i)-fluxtmp(i0)-ftmp(i+1)-(fluxtmp(i-1)-fluxtmp(i0-1)-ftmp(i)))
      flux(i+1)=fluxtmp(i)-fluxtmp(i0)-ftmp(i+1)

      im1=i0;i0=i0+1;if(i0>=ny)i0=i0-ny;

    enddo

    !schema de base (ny-1)
!    ftmp(ny-1)=buf(i0)+fbar-ftmp(ny-1)
    !schema volumes finis (ny-1)
    ftmp(ny-1)=buf(ny-1)-(fluxtmp(ny-1)-fbar-fluxtmp(i0)-(fluxtmp(ny-2)-ftmp(ny-1)-fluxtmp(i0-1)))

    !boundary fluxes
    flux(ny)=fluxtmp(ny-1)-fbar-fluxtmp(i0)
    flux(1)=buf(0)-ftmp(0)+flux(ny)


    fout(1:ny)=ftmp(0:ny-1)

    ! mise a jour de la derniere ligne par periodicite
 !   fout(this%geomx%ny)=fout(1)
  end subroutine evaldepconshermite

end module splinepy_class
