module splinepx_class
  use used_precision
  use geometry_module!1d_module
  implicit none
  private
  public :: new, interpole
  type, public :: splinepx
     type (geometry) :: geomx,geomv
     real(wp) :: a1x, a2x, a3x, a4x ! coef de la matrice 2x2 per.
     real(wp), dimension(:), pointer :: axd ! termes diagonaux de ax
     real(wp), dimension(:), pointer :: axod ! termes sous-diag. de ax
     real(wp), dimension(:,:), pointer :: axm1gamma
     real(wp), dimension(:), pointer :: dper,lper,mper !LU coefficient of the circulant matrix
     real(wp), dimension(:), pointer :: coef,tmpcoef! coefficients des splines
     logical :: transpose       ! permet de definir si f ou ft est derniere 
                                ! fonction de distribution mise a jour
     integer :: jstartx, jendx
     integer :: jstartv, jendv  ! definition de la bande de calcul
  end type splinepx
  interface new
     module procedure new_splinepx
  end interface
  interface interpole
     module procedure interpole_splinepdep
  end interface
contains
  subroutine new_splinepx(this,geom,geomv,iflag,jstartx,jendx,jstartv,jendv,vz)
    type(splinepx), intent(out) :: this
    type(geometry), intent(in) :: geom,geomv  ! geometry of problem
    real(wp), optional  :: vz
    integer, intent(out)       :: iflag ! error flag
    integer, intent(in), optional ::  jstartx,jendx,jstartv,jendv
    ! local variables
    integer :: err ! error flag
    integer :: nx! dimensions
    integer :: nxp1, nxp2,nxp3,Nbdr
    real(wp) :: aa1x, aa2x
    integer i

    ! initialisation des variables locales
    iflag = 0
    nx = geom%nx  !ES  +1
    nxp1=nx+1 ; nxp2=nx+2 ; nxp3=nx+3
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
    allocate(this%axm1gamma(nxp1,2), stat=err)
    if (err.ne.0) then
       iflag = 10
       return
    end if
    allocate(this%axd(nxp1), stat=err)
    if (err.ne.0) then
       iflag = 20
       return
    end if
    allocate(this%axod(nx), stat=err)
    if (err.ne.0) then
       iflag = 30
       return
    end if
    allocate(this%dper(nx), stat=err)
    if (err.ne.0) then
       iflag = 40
       return
    end if
    allocate(this%lper(nx-1), stat=err)
    if (err.ne.0) then
       iflag = 41
       return
    end if
    allocate(this%mper(nx), stat=err)
    if (err.ne.0) then
       iflag = 42
       return
    end if
    allocate(this%coef(nxp3), stat=err)
    if (err.ne.0) then
       iflag = 50
       return
    end if
    allocate(this%tmpcoef(nx+2*Nbdr), stat=err)
    if (err.ne.0) then
       iflag = 51
       return
    end if

    ! factorisation de la  matrice  Ax 
    this%axd = 4.0_wp
    this%axod = 1.0_wp
    call dpttrf(nxp1,this%axd,this%axod,err)
    if (err.ne.0) then
       iflag = 70
       return
    end if
    ! compute Ax-1.gamma
    this%axm1gamma = 0.0_wp
    this%axm1gamma(1,2) = 1.0_wp
    this%axm1gamma(nxp1,1) = 1.0_wp
    call dpttrs(nxp1,2,this%axd,this%axod,this%axm1gamma,nxp1,err)
    if (err.ne.0) then
       iflag = 90
       return
    end if

    aa1x=3.0_wp/geom%dx
    aa2x=6.0_wp/(geom%dx*geom%dx)
    ! assemblage de la matrice 2x2 pour la spline periodique
    this%a1x = -aa1x*(1._wp + this%axm1gamma(2,1)+this%axm1gamma(nx,1))
    this%a2x = -aa1x*(1._wp + this%axm1gamma(2,2)+this%axm1gamma(nx,2))
    this%a3x = aa2x*(-1._wp + 2._wp*this%axm1gamma(1,1) - this%axm1gamma(2,1) &
         + this%axm1gamma(nx,1) - 2._wp*this%axm1gamma(nxp1,1))
    this%a4x = aa2x*(1._wp + 2._wp*this%axm1gamma(1,2) - this%axm1gamma(2,2) &
         + this%axm1gamma(nx,2) - 2._wp*this%axm1gamma(nxp1,2))

    this%dper(1)=4._wp;this%mper(1)=0.25_wp
    do i=1,nx-1
       this%lper(i)=1._wp/this%dper(i)
       this%dper(i+1)=4._wp-this%lper(i)
       this%mper(i+1)=-this%mper(i)/this%dper(i+1)
    enddo
    this%dper(nx)=this%dper(nx)-this%lper(nx-1)+2._wp*this%mper(nx-1)
    do i=1,nx
       this%dper(i)=1._wp/this%dper(i)
    enddo


 end subroutine new_splinepx


subroutine interpole_splinepdep(this,f,depx,G,meth) 
    !----------------------------------------------------------------
    ! interpolation par spline periodique.
    ! Les points d'interpolation sont definis grace a depx 
    ! qui definissent le deplacement par rapport au maillage.
    !----------------------------------------------------------------
    type(splinepx), intent(inout) :: this
    ! f contient les valeurs de la fonction de distribution
    real(wp), dimension(:), intent(inout) :: f
    real(wp), dimension(:), intent(inout) :: G
    ! dans depx on trouve les deplacements par rapport au maillage
    ! des points dans les quels on veut evaluer la spline.
    real(wp), intent(in) :: depx
    integer, intent(in) :: meth
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    real(wp), dimension(1:this%geomx%nx) :: tmpf
    integer :: ierr,i,nx
    
    if (meth==0) then 
       !BSL formulation cubic splines 
       call per_x(this,f,ierr)
       
       call evaldep(this,depx,f)    
    end if

    if (meth>=1) then 
       !methodes conservatives 


       !construction de G, primitive de f
       nx=this%geomx%nx
       G(1)=f(1)
       do i=2,nx
          G(i)=G(i-1)+f(i)
       enddo

       if (meth==1) then 
          !splines
          G(1)=G(1)+(1._wp/6._wp)*G(nx)
          G(nx)=(5._wp/6._wp)*G(nx)
          
          call per_primx(this,G,ierr)
          
          if (ierr.ne.0) then
             iflag = 10
             return
          end if
          
          G(nx)=(6._wp/5._wp)*G(nx)
          G(1)=G(1)-(1._wp/6._wp)*G(nx)
          
          tmpf(1:nx)=f(1:nx)
          call evaldepcons(this,depx,f,G)    
          !VF formulation
          G(1:nx)=G(1:nx)*this%geomx%dx/depx
          do i=2,nx
             f(i)=tmpf(i)-depx*(G(i)-G(i-1))/this%geomx%dx
          enddo
          f(1)=tmpf(1)-depx*(G(1)-G(nx))/this%geomx%dx
       endif

       if (meth>=2) then 
          !ppm methods

          tmpf(1:nx)=f(1:nx)
          call evaldepconshermite(this,depx,f,G,meth)     
          
          !VF formulation
          G(1:nx)=G(1:nx)*this%geomx%dx/depx
          do i=2,nx
             f(i)=tmpf(i)-depx*(G(i)-G(i-1))/this%geomx%dx
          enddo
          f(1)=tmpf(1)-depx*(G(1)-G(nx))/this%geomx%dx

       endif

    endif

  end subroutine interpole_splinepdep

  !
  ! calcul des splines periodiques
  !
  subroutine per_x(this,gtau,iflag)
    type(splinepx), intent(inout) :: this     ! objet de type splinep
    real(wp), dimension(:), intent(in) :: gtau ! valeur de la fonction 
    ! aux points du maillage
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i ! indices de boucle
    integer nx, nxp2, nxp3  !  nx+2, nx+3
    integer :: err ! error flag

    real(wp) :: axm1f(this%geomx%nx+1)  !ES +1-> +2
    real(wp) :: det, gamma1, gamma2

    ! initialisations
    iflag =0
    nx=this%geomx%nx  !ES +1
    nxp2=nx+2
    nxp3=nx+3
    det=this%a1x*this%a4x - this%a2x*this%a3x

    ! Calcul de Ax^-1 f
    ! assemblage du membre de droite pour le calcul de Ax^-1 f
    do i=1,nx
       axm1f(i) = 6._wp*gtau(i)
    end do
    axm1f(nx+1) = 6._wp*gtau(1) ! Calcul par periodicite   
    ! on calcule Ax^-1f
    call dpttrs(nx+1,1,this%axd,this%axod,axm1f,nx+1,err)
    if (err.ne.0) then
       iflag = 10
       return
    end if
    ! assemblage du second membre du systeme 2x2 
 
    gamma1 =  - (3._wp/this%geomx%dx)*(axm1f(2) +axm1f(nx))
    gamma2 =  (6._wp/(this%geomx%dx)**2._wp)*(2._wp*axm1f(1) - axm1f(2) &
         + axm1f(nx) - 2._wp*axm1f(nx+1))
    
    this%coef(nxp3)= (gamma1*this%a4x - gamma2*this%a2x)/det
    this%coef(1)= (gamma2*this%a1x - gamma1*this%a3x)/det
    do  i=2,nxp2
       this%coef(i)= axm1f(i-1) &
            - this%axm1gamma(i-1,1)*this%coef(nxp3) &
            - this%axm1gamma(i-1,2)*this%coef(1)
    end do
  end subroutine per_x

  !
  ! calcul des splines periodiques par inversion de la matrice circulante
  !
  subroutine per_primx(this,gtau,iflag)
    type(splinepx), intent(inout) :: this     ! objet de type splinep
    real(wp), dimension(:), intent(in) :: gtau ! valeur de la fonction 
    ! aux points du maillage
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i! indices de boucle
    integer nx   !  nx+2, nx+3

    real(wp) :: axm1f(this%geomx%nx+1)  !ES +1-> +2
    real(wp) :: det,M

    ! initialisations
    iflag =0
    nx=this%geomx%nx  !ES +1

    ! Calcul de Ax^-1 f
    ! assemblage du membre de droite pour le calcul de Ax^-1 f

    M=(6._wp/5._wp)*gtau(nx)

    axm1f=0._wp
    do i=1,nx
       axm1f(i) = 6._wp*gtau(i)
    end do
    do i=2,nx
       axm1f(i)=axm1f(i)-axm1f(i-1)*this%lper(i-1)
    enddo
    do i=1,nx-1
       axm1f(nx)=axm1f(nx)-this%mper(i)*axm1f(i)
    enddo
    axm1f(nx)=axm1f(nx)*this%dper(nx)
    axm1f(nx-1)=this%dper(nx-1)*(axm1f(nx-1)-(1._wp-this%mper(nx-2))*axm1f(nx))
    do i=nx-2,2,-1
       axm1f(i)=this%dper(i)*(axm1f(i)-axm1f(i+1)+this%mper(i-1)*axm1f(nx))
    enddo
    axm1f(1)=this%dper(1)*(axm1f(1)-axm1f(2)-axm1f(nx))

    do  i=1,nx
       this%coef(i)= axm1f(i) 
    end do
    !calcul des coef aux bords par "periodicite"
    !coef(nx+1)
    this%coef(nx+1)=this%coef(1)+M 
    !coef(0)
    this%coef(nx+2)=this%coef(nx)-M 
    !masse totale
    this%coef(nx+3)=M !gtau(nx)

    do  i=2,nx-1
       if (abs(this%coef(i-1)+4._wp*this%coef(i)+this%coef(i+1)-6._wp*gtau(i))>1.e-13) then 
          print *,'verif lu per_primx',abs(this%coef(i-1)+4._wp*this%coef(i)+this%coef(i+1)-gtau(i))
          print *,'coef',this%coef(i-1),this%coef(i),this%coef(i+1),i,axm1f(i)
          stop
       endif
    enddo

  end subroutine per_primx
  
  subroutine evaldep(this,alphax,fout)
    type(splinepx) :: this
    ! deplacement par rapport aux maillages des points ou la spline est evaluee
    real(wp), intent(in) :: alphax  !valeur du deplacement
    ! fout(i) contient en sortie la valeur de la spline au 
    ! point (xi-alphax), les (xi) etant les points du maillage
    real(wp), dimension(:), intent(out) :: fout

    ! variables locales
    real(wp) :: sval   ! valeur de la fonction au point d'evaluation
    real(wp) :: bvalx1,bvalx2,bvalx3,bvalx4,dxx,dxxx,dxxx6,xd1,xdp1
    real(wp) :: sval1, sval2, sval3, sval4
    integer :: intaxsdx
    integer i1,i

    ! debut du code
    dxx=this%geomx%dx*this%geomx%dx
    dxxx=dxx*this%geomx%dx
    dxxx6=1._wp/(6._wp*dxxx)

    if (alphax.gt.0) then
       intaxsdx=int(-alphax/this%geomx%dx+epsilon)-1  
! intaxsdx=int(-alphax/this%geomx%dx)-1 
!       intaxsdx = -1
    else
       intaxsdx=int(-alphax/this%geomx%dx)
    end if
    
    xd1=-alphax-intaxsdx*this%geomx%dx
    xdp1=this%geomx%dx-xd1

    bvalx1=xdp1*xdp1*xdp1
    bvalx2=dxxx+3._wp*dxx*xdp1+3._wp*this%geomx%dx*xdp1*xdp1-3._wp*xdp1*xdp1*xdp1
    bvalx3=dxxx+3._wp*dxx*xd1+3._wp*this%geomx%dx*xd1*xd1-3._wp*xd1*xd1*xd1
    bvalx4=xd1*xd1*xd1

    do i=1,this%geomx%nx
       i1=mod(this%geomx%nx+i-1+intaxsdx,this%geomx%nx)
       fout(i) = dxxx6*(bvalx1*this%coef(i1+1)+bvalx2*this%coef(i1+2) &
               + bvalx3*this%coef(i1+3) + bvalx4*this%coef(i1+4))
    end do
    ! mise a jour de la derniere ligne par periodicite
 !   fout(this%geomx%nx)=fout(1)
  end subroutine evaldep


  subroutine evaldepcons(this,alphax,fout,flux)
    type(splinepx) :: this
    ! deplacement par rapport aux maillages des points ou la spline est evaluee
    real(wp), intent(in) :: alphax  !valeur du deplacement
    ! fout(i) contient en sortie la valeur de la spline au 
    ! point (xi-alphax), les (xi) etant les points du maillage
    real(wp), dimension(:), intent(out) :: fout
    real(wp), dimension(:), intent(inout) :: flux

    ! variables locales
    integer :: ix,i,nx,Nbdr
    real(wp) :: w(0:3)
    real(wp) :: alf,dx,x,M,dtmp

    Nbdr=10
    nx=this%geomx%nx

    M=this%coef(nx+3)

    do i=1,nx
       this%tmpcoef(i+Nbdr) = this%coef(i)
    enddo
    do i=1,Nbdr
       this%tmpcoef(i)=this%tmpcoef(i+nx)-M
    enddo
    do i=1,Nbdr
       this%tmpcoef(i+nx+Nbdr)=this%tmpcoef(i+Nbdr)+M
    enddo
    !dep normalise
    alf=alphax/(this%geomx%x1-this%geomx%x0);dx=1._wp/real(nx,wp)
    do i=1,this%geomx%nx
       x=i*dx-alf;ix=floor(x*nx);x=x*nx-ix
       w(0)=(1._wp/6._wp)*(1._wp-x)*(1._wp-x)*(1._wp-x);
       w(1)=(1._wp/6._wp)+0.5_wp*(1._wp-x)*(-(1._wp-x)*(1._wp-x)+(1._wp-x)+1._wp)
       w(2)=(1._wp/6._wp)+0.5_wp*x*(-x*x+x+1._wp)
       w(3)=(1._wp/6._wp)*x*x*x
       
       fout(i) = w(0)*this%tmpcoef(ix-1+Nbdr) + w(1)*this%tmpcoef(ix+Nbdr) &
            + w(2)*this%tmpcoef(ix+1+Nbdr) + w(3)*this%tmpcoef(ix+2+Nbdr) 
    end do
    !fout(i)=F(x_{i+1/2}^*)

    !flux(i)=F(x_{i+1/2})-F(x_{i+1/2}^*)
    do i=1,this%geomx%nx
       flux(i)=flux(i)-fout(i)
    enddo

    !update de f^{n+1}_i=F(x_{i+1/2}^*)-F(x_{i-1/2}^*)
    dtmp=fout(nx)
    do i=nx,2,-1
       fout(i)=fout(i)-fout(i-1)
    enddo
    fout(1)=fout(1)+M-dtmp

 
  end subroutine evaldepcons


  subroutine evaldepconshermite(this,alphax,fout,flux,meth)
    type(splinepx) :: this
    ! deplacement par rapport aux maillages des points ou la spline est evaluee
    real(wp), intent(in) :: alphax  !valeur du deplacement
    ! fout(i) contient en sortie la valeur de la spline au 
    ! point (xi-alphax), les (xi) etant les points du maillage
    real(wp), dimension(:), intent(inout) :: fout
    real(wp), dimension(:), intent(inout) :: flux
    integer, intent(in) :: meth

    ! variables locales
    real(wp), dimension(0:this%geomx%nx) :: ftmp
    real(wp), dimension(0:this%geomx%nx) :: buf
    real(wp), dimension(-1:this%geomx%nx) :: fluxtmp
    integer :: ix,i,nx,Nbdr,im3,im2,im1,i0,i1,i2,i3
    real(wp) :: w(0:3),df0,df1,fbar
    real(wp) :: alf,dx,x,M,dtmp,mass,flux0

    Nbdr=10
    nx=this%geomx%nx
    
    fluxtmp(0:this%geomx%nx-1)=flux(1:this%geomx%nx)
    fluxtmp(this%geomx%nx)=fluxtmp(0)+fluxtmp(this%geomx%nx-1)
    fluxtmp(-1)=0._8

    ftmp(0:this%geomx%nx-1)=fout(1:this%geomx%nx)
    ftmp(this%geomx%nx)=ftmp(0)

    buf(0:nx-1)=fout(1:nx)
    buf(nx)=buf(0)

    !dep normalise
    alf=alphax/(this%geomx%x1-this%geomx%x0);dx=1._wp/real(nx,wp)

    x=-alf
    x=x-real(floor(x),8)
    x=x*real(nx,8)
    ix=floor(x)
    if(ix==nx)then
      x=0._8;ix=0
    endif
    x=x-real(ix,8)
    
    w(0)=x*(1._8-x)*(1._8-x)
    w(1)=x*x*(x-1._8)
    w(2)=x*x*(3._8-2._8*x)
    i0=ix
    i1=mod(ix+1,nx)
    i2=mod(ix+2,nx)
    i3=mod(ix+3,nx)
    im1=mod(ix+nx-1,nx)
    im2=mod(ix+nx-2,nx)
    im3=mod(ix+nx-3,nx)    


    do i=0,this%geomx%nx-1
       
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
       im3=im2;im2=im1;im1=i0;i0=i1;i1=i1+1;if(i1>=nx)i1=i1-nx;i2=i2+1;if(i2>=nx)i2=i2-nx;i3=i3+1;if(i3>=nx)i3=i3-nx;

    end do

    i0=ix;im1=mod(i0-1+nx,nx)
    fbar=ftmp(0);


    do i=0,nx-2
       !schema de base
!      ftmp(i)=buf(i0)+ftmp(i+1)-ftmp(i)

       !schema volumes finis
       ftmp(i)=buf(i)-(fluxtmp(i)-fluxtmp(i0)-ftmp(i+1)-(fluxtmp(i-1)-fluxtmp(i0-1)-ftmp(i)))

       flux(i+1)=fluxtmp(i)-fluxtmp(i0)-ftmp(i+1)

       im1=i0;i0=i0+1;if(i0>=nx)i0=i0-nx
    enddo

    !schema de base (nx-1)
!    ftmp(nx-1)=buf(i0)+fbar-ftmp(nx-1)
    !schema volumes finis (nx-1)
    ftmp(nx-1)=buf(nx-1)-(fluxtmp(nx-1)-fbar-fluxtmp(i0)-(fluxtmp(nx-2)-ftmp(nx-1)-fluxtmp(i0-1)))

    !boundary fluxes
    flux(nx)=fluxtmp(nx-1)-fbar-fluxtmp(i0)
    flux(1)=buf(0)-ftmp(0)+flux(nx)

    fout(1:nx)=ftmp(0:nx-1)


!!verif
!!$    mass=abs(ftmp(1)-(buf(1)-(flux(2)-flux(1))))
!!$    mass=abs(ftmp(0)-(buf(0)-(flux(1)-flux(nx))))
!!$    if (mass>1.e-12) then 
!!$       print *,'err i',i,mass,nx-1,i0,sum(buf(0:nx-1)),ftmp(i),flux(i)!,fluxtmp(i0),fbar
!!$       stop
!!$    end if
!!$
!!$    do i=1,nx-1
!!$       mass=abs(ftmp(i)-(buf(i)-(flux(i+1)-flux(i))))
!!$       if (mass>1.e-12) then 
!!$          print *,'err i',i,mass,nx-1,i0,sum(buf(0:nx-1)),ftmp(i),flux(i)!,fluxtmp(i0),fbar
!!$          stop
!!$       end if
!!$    enddo



  end subroutine evaldepconshermite

end module splinepx_class
