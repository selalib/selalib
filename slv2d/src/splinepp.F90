module splinepp_class
  use used_precision
  use geometry_module
  !use clock
  implicit none
  private
  public :: new, interpole
  type, public :: splinepp
     type (geometry) :: geom
     real(wp) :: a1x, a2x, a3x, a4x ! coef de la matrice 2x2 per.
     real(wp) :: a1y, a2y, a3y, a4y ! coef de la matrice 2x2 nat.
     real(wp), dimension(:), pointer :: axd, ayd ! termes diagonaux de ax
     real(wp), dimension(:), pointer :: axod, ayod ! termes sous-diag. de ax
     real(wp), dimension(:,:), pointer :: axm1gamma, aym1gamma 
     real(wp), dimension(:,:), pointer :: coef, bcoef ! coefficients des splines
  end type splinepp
  interface new
     module procedure new_splinepp
  end interface
  interface interpole
     module procedure interpole_splinepp,interpole_splineppdep
  end interface
contains
  subroutine new_splinepp(this,geom,iflag)
    type(splinepp), intent(out) :: this
    type(geometry), intent(in) :: geom  ! geometry of problem
    integer, intent(out)       :: iflag ! error flag
    ! local variables
    integer :: err ! error flag
    integer :: nx,ny ! dimensions
    integer :: nxp1, nxp2, nxp3, nyp1, nyp2, nyp3
    real(wp) :: aa1x, aa1y, aa2x, aa2y
    integer i,j, info
    ! initialisation des variables locales
    iflag = 0
    nx = geom%nx
    ny = geom%ny
    nxp1=nx+1
    nxp2=nx+2
    nxp3=nx+3
    nyp1=ny+1
    nyp2=ny+2
    nyp3=ny+3

    ! initialisation de geom
    this%geom = geom

    ! memory allocation
    allocate(this%axm1gamma(nxp1,2), stat=err)
    if (err.ne.0) then
       iflag = 10
       return
    end if
    allocate(this%aym1gamma(nyp1,2), stat=err)
    if (err.ne.0) then
       iflag = 15
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
    allocate(this%ayd(nyp1), stat=err)
    if (err.ne.0) then
       iflag = 40
       return
    end if
    allocate(this%ayod(ny), stat=err)
    if (err.ne.0) then
       iflag = 45
       return
    end if
    allocate(this%coef(nxp3,nyp3), stat=err)
    if (err.ne.0) then
       iflag = 50
       return
    end if
    allocate(this%bcoef(nxp3,ny), stat=err)
    if (err.ne.0) then
       iflag = 60
       return
    end if

    ! factorize matrices  Ax and Ay
    this%axd = 4_wp
    this%axod = 1_wp
#ifdef _CRAY
    call spttrf(nxp1,this%axd,this%axod,err)
#else
    call dpttrf(nxp1,this%axd,this%axod,err)
#endif
    if (err.ne.0) then
       iflag = 70
       return
    end if
    this%ayd = 4_wp
    this%ayod = 1_wp
#ifdef _CRAY
    call spttrf(nyp1,this%ayd,this%ayod,err)
#else
    call dpttrf(nyp1,this%ayd,this%ayod,err)
#endif
    if (err.ne.0) then
       iflag = 80
       return
    end if
    ! compute Ax-1.gamma
    this%axm1gamma = 0_wp
    this%axm1gamma(1,2) = 1_wp
    this%axm1gamma(nxp1,1) = 1_wp
#ifdef _CRAY
    call spttrs(nxp1,2,this%axd,this%axod,this%axm1gamma,nxp1,err)
#else
    call dpttrs(nxp1,2,this%axd,this%axod,this%axm1gamma,nxp1,err)
#endif
    if (err.ne.0) then
       iflag = 90
       return
    end if
    ! compute Ay-1.gamma
    this%aym1gamma = 0_wp
    this%aym1gamma(1,2) = 1_wp
    this%aym1gamma(nyp1,1) = 1_wp
#ifdef _CRAY
    call spttrs(nyp1,2,this%ayd,this%ayod,this%aym1gamma,nyp1,err)
#else
    call dpttrs(nyp1,2,this%ayd,this%ayod,this%aym1gamma,nyp1,err)
#endif
    if (err.ne.0) then
       iflag = 100
       return
    end if

    aa1x=3_wp/geom%dx
    aa1y=3_wp/geom%dy
    aa2x=6_wp/(geom%dx*geom%dx)
    aa2y=6_wp/(geom%dy*geom%dy)
    ! assemblage de la matrice 2x2 pour la spline dans la direction Ox
       this%a1x = -aa1x*(1. + this%axm1gamma(2,1)+this%axm1gamma(nx,1))
       this%a2x = -aa1x*(1. + this%axm1gamma(2,2)+this%axm1gamma(nx,2))
       this%a3x = aa2x*(-1. + 2*this%axm1gamma(1,1) - this%axm1gamma(2,1) &
            + this%axm1gamma(nx,1) - 2*this%axm1gamma(nxp1,1))
       this%a4x = aa2x*(1. + 2*this%axm1gamma(1,2) - this%axm1gamma(2,2) &
            + this%axm1gamma(nx,2) - 2*this%axm1gamma(nxp1,2))
    ! assemblage de la matrice 2x2 pour spline naturels (direction Oy)
       this%a1y = -aa1y*(1. + this%aym1gamma(2,1)+this%aym1gamma(ny,1))
       this%a2y = -aa1y*(1. + this%aym1gamma(2,2)+this%aym1gamma(ny,2))
       this%a3y = aa2y*(-1. + 2*this%aym1gamma(1,1) - this%aym1gamma(2,1) &
            + this%aym1gamma(ny,1) - 2*this%aym1gamma(nyp1,1))
       this%a4y = aa2y*(1. + 2*this%aym1gamma(1,2) - this%aym1gamma(2,2) &
            + this%aym1gamma(ny,2) - 2*this%aym1gamma(nyp1,2))
 end subroutine new_splinepp

  subroutine interpole_splinepp(this,fin,fout,x,y) 
    type(splinepp), intent(inout) :: this
    ! fin contient les valeurs de la fonction dans la grille precedente
    real(wp), dimension(:,:), intent(in) :: fin
    ! fout est destine a contenir la nouvelle valeur de f
    real(wp), dimension(:,:), intent(out):: fout
    ! dans x et y on trouve les points dans les quels on veut 
    ! evaluer la spline.
    real(wp), dimension(:,:), intent(in) :: x, y 
    ! dans fout on trouve en sortie les valeurs de f(i,j) 
    ! dans les points x(i),y(i).
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    integer ierr

    call per_x(this,fin,ierr)

    if (ierr.ne.0) then
       iflag = 10
       return
    end if

    call per_y(this,ierr)

    if (ierr.ne.0) then
       iflag = 20
       return
    end if

    call evaltab(this,x,y,fout)     

  end subroutine interpole_splinepp

  subroutine interpole_splineppdep(this,f,depx,depy,aff) 
    !----------------------------------------------------------------
    ! interpolation par spline periodique dans les deux directions.
    ! Les points d'interpolation sont definis grace a depx et depy
    ! qui definissent le deplacement par rapport au maillage.
    !----------------------------------------------------------------
    type(splinepp), intent(inout) :: this
    ! f contient les valeurs de la fonction de distribution
    real(wp), dimension(:,:), intent(inout) :: f
    ! dans depx et depy on trouve les deplacements par rapport au maillage
    ! des points dans les quels on veut evaluer la spline.
    real(wp), intent(in) :: depx, depy 
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    logical  :: aff
    integer ierr, l_a, l_b
    real(wp) vtime(1:4)

    !if (aff) then 
    !   call clck_temps(l_a)
    !end if

    call per_x(this,f,ierr)
    if (ierr.ne.0) then
       iflag = 10
       return
    end if
    !if (aff) then 
    !   call clck_temps(l_b)
    !   call clck_diff(l_a,l_b,vtime(1))
    !   call clck_temps(l_a)
    !end if

    call per_y(this,ierr)
    if (ierr.ne.0) then
       iflag = 20
       return
    end if
    !if (aff) then 
    !   call clck_temps(l_b)
    !   call clck_diff(l_a,l_b,vtime(2))
    !   call clck_temps(l_a)
    !end if

    call evaldep(this,depx,depy,f)     

    !if (aff) then 
    !   call clck_temps(l_b)
    !   call clck_diff(l_a,l_b,vtime(3))
    !   write(*,'(A,3(1x,3E14.5))') "splinepp ",vtime(1:3)
    !end if

  end subroutine interpole_splineppdep

  !
  ! calcul des splines periodiques
  !
  subroutine per_x(this,gtau,iflag)
    type(splinepp), intent(inout) :: this     ! objet de type splinepp
    real(wp), dimension(:,:), intent(in) :: gtau ! valeur de la fonction 
    ! aux points du maillage
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i,j ! indices de boucle
    integer nx, ny, nxp2, nxp3  !  nx+2, nx+3
    integer :: err ! error flag

    real(wp) :: axm1f(this%geom%nx+1,this%geom%ny)
    real(wp) :: det, gamma1, gamma2

    ! initialisations
    iflag =0
    nx=this%geom%nx
    ny=this%geom%ny
    nxp2=nx+2
    nxp3=nx+3
    det=this%a1x*this%a4x - this%a2x*this%a3x

    ! Calcul de Ax^-1 f
    ! assemblage du membre de droite pour le calcul de Ax^-1 f
    do j=1,ny
       do i=1,nx
          axm1f(i,j) = 6*gtau(i,j)
       end do
       axm1f(nx+1,j) = 6*gtau(1,j)  ! calcul par periodicite
    end do

#ifdef _CRAY
    call spttrs(nx+1,ny,this%axd,this%axod,axm1f,nx+1,err)
#else 
    call dpttrs(nx+1,ny,this%axd,this%axod,axm1f,nx+1,err)
#endif
    if (err.ne.0) then
       iflag = 10
       return
    end if
    do  j=1,ny
       ! assemblage du second membre du systeme 2x2 
       gamma1 =  - (3.0/this%geom%dx)*(axm1f(2,j) +axm1f(this%geom%nx,j))
       gamma2 =  (6.0/(this%geom%dx)**2)*(2*axm1f(1,j) - axm1f(2,j) &
            + axm1f(this%geom%nx,j) - 2*axm1f(this%geom%nx+1,j))

       this%bcoef(nxp3,j)= (gamma1*this%a4x - gamma2*this%a2x)/det
       this%bcoef(1,j)= (gamma2*this%a1x - gamma1*this%a3x)/det
       do  i=2,nxp2
          this%bcoef(i,j)= axm1f(i-1,j) &
               - this%axm1gamma(i-1,1)*this%bcoef(nxp3,j) &
               - this%axm1gamma(i-1,2)*this%bcoef(1,j)
       end do
    end do
  end subroutine per_x
  !
  subroutine per_y(this,iflag)
    type(splinepp), intent(inout) :: this     ! objet de type spline
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i,j ! indices de boucle
    integer ny, nxp3, nyp3  !  nx+3, ny+3
    integer :: err ! error flag

    real(wp) :: aym1f(this%geom%ny+1,this%geom%nx+3)
    real(wp) :: det, gamma1, gamma2

    ! initialisations
    iflag =0
    ny = this%geom%ny
    nxp3=this%geom%nx+3
    nyp3=this%geom%ny+3
    det=this%a1y*this%a4y - this%a2y*this%a3y

    ! calcul de coef par resolution de nxp2 systemes lineaires.

    ! Calcul de Ay^-1 f
    ! assemblage du membre de droite pour le calcul de Ay^-1 f
    do i=1,nxp3
       do j=1,ny
          aym1f(j,i) = 6.*this%bcoef(i,j)
       end do
       aym1f(ny+1,i) = 6.*this%bcoef(i,1)
    end do
#ifdef _CRAY
    call spttrs(ny+1,nxp3,this%ayd,this%ayod,aym1f,ny+1,err)
#else
    call dpttrs(ny+1,nxp3,this%ayd,this%ayod,aym1f,ny+1,err)
#endif
    if (err.ne.0) then
       iflag = 10
       return
    end if
    do i=1,nxp3
       ! assemblage du second membre du systeme 2x2 
       gamma1 =  - (3.0/this%geom%dy)*(aym1f(2,i) +aym1f(this%geom%ny,i))
       gamma2 =  (6.0/(this%geom%dy)**2)*(2*aym1f(1,i) - aym1f(2,i) &
            + aym1f(this%geom%ny,i) - 2*aym1f(this%geom%ny+1,i))
       ! resolution du syteme lineaire 2x2
       this%coef(i,nyp3) = (gamma1*this%a4y - gamma2*this%a2y)/det
       this%coef(i,1) = (gamma2*this%a1y - gamma1*this%a3y)/det
       do  j=2,this%geom%ny + 2
          this%coef(i,j)= aym1f(j-1,i)                   &
               - this%aym1gamma(j-1,1)*this%coef(i,nyp3) &
               - this%aym1gamma(j-1,2)*this%coef(i,1)
       end do
    end do
  end subroutine per_y

  subroutine evaltab(this,xd,yd,fout)
    ! ------------------------------------------------------
    ! Evalue la spline en tous les points (xd(i,j), yd(i,j))
    !-------------------------------------------------------
    type(splinepp) :: this
    ! coordonnees du point ou les valeurs sont calculees
    real(wp), dimension(:,:) :: xd, yd 
    ! fout(i,j)contient la valeur de la spline au point (xd(i,j), yd(i,j))
    real(wp), dimension(:,:) :: fout

    ! variables locales
    real(wp) :: sval   ! valeur de la fonction au point (xd,yd)
    real(wp) :: bvalx1,bvalx2,bvalx3,bvalx4,bvaly1,bvaly2,bvaly3,bvaly4
    real(wp) :: a1,dxx,dxxx,dxxx6,dyy,dyyy,dyyy6,xd1,xdp1,yd1,ydp1
    real(wp) :: sval1, sval2, sval3, sval4, idx, idy, lx, ly
    integer i1,j1,i,j
    !
    dxx=this%geom%dx*this%geom%dx
    dxxx=dxx*this%geom%dx
    dxxx6=1./(6.*dxxx)
    !
    !
    dyy=this%geom%dy*this%geom%dy
    dyyy=dyy*this%geom%dy
    dyyy6=1./(6.*dyyy)

    idx = 1/this%geom%dx
    idy = 1/this%geom%dy

    lx = (this%geom%nx-1)*this%geom%dx
    ly = (this%geom%ny-1)*this%geom%dy
    !print*,'l ',lx,ly

    ! evaluation des splines pour le calcul des nx-1 lignes 
    ! et ny-1 colonnes de fout
    do j=1,this%geom%ny
       do i=1,this%geom%nx
!!$          do while (xd(i,j)<this%geom%x0)
!!$             !print*,i,j,xd(i,j),this%geom%x0
!!$             xd(i,j)=xd(i,j)+lx
!!$          end do
!!$          do while (yd(i,j)<this%geom%y0)
!!$             !print*,i,j,yd(i,j),this%geom%y0
!!$             yd(i,j)=yd(i,j)+ly
!!$          end do
!!$          do while (xd(i,j)>this%geom%x0+lx)
!!$             !print*,i,j,xd(i,j),this%geom%x0
!!$             xd(i,j)=xd(i,j)-lx
!!$          end do
!!$          do while (yd(i,j)>this%geom%y0+ly)
!!$             !print*,i,j,yd(i,j),this%geom%y0
!!$             yd(i,j)=yd(i,j)-ly
!!$          end do

          i1=(xd(i,j)-this%geom%x0)*idx
          j1=(yd(i,j)-this%geom%y0)*idy

          if ((i1<0).or.(i1>this%geom%nx-1)) print*,'i1',i1
          if ((j1<0).or.(j1>this%geom%ny-1)) print*,'j1',j1
          
          !i1=mod(i1+this%geom%nx,this%geom%nx)
          !j1=mod(j1+this%geom%ny,this%geom%ny)
          !print*,'rrr ',i,j,xd(i,j),yd(i,j),i1,j1,this%geom%nx,this%geom%ny

          xdp1=this%geom%xgrid(i1+2)-xd(i,j)
          bvalx1=xdp1*xdp1*xdp1
          bvalx2=dxxx+3.*dxx*xdp1+3.*this%geom%dx* &
               &xdp1*xdp1-3.*xdp1*xdp1*xdp1
          xd1=xd(i,j)-this%geom%xgrid(i1+1)
          bvalx3=dxxx+3.*dxx*xd1+3.*this%geom%dx* &
               &xd1*xd1-3.*xd1*xd1*xd1
          bvalx4=xd1*xd1*xd1
          ydp1=this%geom%ygrid(j1+2)-yd(i,j)
          bvaly1=ydp1*ydp1*ydp1
          bvaly2=dyyy+3.*ydp1*(dyy+ydp1*(this%geom%dy-ydp1))          
          yd1=yd(i,j)-this%geom%ygrid(j1+1)
          bvaly3=dyyy+3.*yd1*(dyy+yd1*(this%geom%dy-yd1))
          bvaly4=yd1*yd1*yd1

          sval=0.
          sval1=this%coef(i1+1,j1+1)*bvaly1
          sval1=sval1+this%coef(i1+1,j1+2)*bvaly2
          sval1=sval1+this%coef(i1+1,j1+3)*bvaly3
          sval1=sval1+this%coef(i1+1,j1+4)*bvaly4        
          sval= sval+sval1*bvalx1
          sval2=this%coef(i1+2,j1+1)*bvaly1
          sval2=sval2+this%coef(i1+2,j1+2)*bvaly2
          sval2=sval2+this%coef(i1+2,j1+3)*bvaly3
          sval2=sval2+this%coef(i1+2,j1+4)*bvaly4
          sval=sval+sval2*bvalx2
          sval3=this%coef(i1+3,j1+1)*bvaly1
          sval3=sval3+this%coef(i1+3,j1+2)*bvaly2
          sval3=sval3+this%coef(i1+3,j1+3)*bvaly3
          sval3=sval3+this%coef(i1+3,j1+4)*bvaly4
          sval=sval+sval3*bvalx3
          sval4=this%coef(i1+4,j1+1)*bvaly1
          sval4=sval4+this%coef(i1+4,j1+2)*bvaly2
          sval4=sval4+this%coef(i1+4,j1+3)*bvaly3
          sval4=sval4+this%coef(i1+4,j1+4)*bvaly4
          sval=sval+sval4*bvalx4
         
          fout(i,j) = dxxx6*dyyy6*sval
       end do
    end do
 
  end subroutine evaltab

  subroutine evaldep(this,alphax,alphay,fout)
    type(splinepp) :: this
    ! deplacement par rapport aux maillages des points ou la spline est evaluee
    real(wp), intent(in) :: alphax,alphay 
    ! fout(i,j) contient en sortie la valeur de la spline au 
    ! point (xi-alphax,yj-alphay), les (xi,yj) etant les points du maillage
    real(wp), dimension(:,:), intent(out) :: fout

    ! variables locales
    real(wp) :: sval   ! valeur de la fonction au point d'evaluation
    real(wp) bvalx1,bvalx2,bvalx3,bvalx4,bvaly1,bvaly2,bvaly3, &
         &bvaly4,dxx,dxxx,dxxx6,dyy,dyyy,dyyy6,xd1,xdp1,yd1,ydp1
    real(wp) :: sval1, sval2, sval3, sval4
    integer :: intaxsdx, intaysdy 
    integer i1,j1,i,j

    ! debut du code
    dxx=this%geom%dx*this%geom%dx
    dxxx=dxx*this%geom%dx
    dxxx6=1./(6.*dxxx)

    dyy=this%geom%dy*this%geom%dy
    dyyy=dyy*this%geom%dy
    dyyy6=1./(6.*dyyy)

    if (alphax.gt.0) then
       intaxsdx=int(-alphax/this%geom%dx+epsilon)-1  
! intaxsdx=int(-alphax/this%geom%dx)-1 
!       intaxsdx = -1
    else
       intaxsdx=int(-alphax/this%geom%dx)
    end if
    
    if (alphay.gt.0) then
       intaysdy=int(-alphay/this%geom%dy+epsilon)-1
!intaysdy=int(-alphax/this%geom%dx)-1    
!intaysdy=-1
    else
       intaysdy=int(-alphay/this%geom%dy)
    end if
    xd1=-alphax-intaxsdx*this%geom%dx
    xdp1=this%geom%dx-xd1
    yd1=-alphay-intaysdy*this%geom%dy
    ydp1=this%geom%dy-yd1

    bvalx1=xdp1*xdp1*xdp1
    bvalx2=dxxx+3.*dxx*xdp1+3.*this%geom%dx* &
         &xdp1*xdp1-3.*xdp1*xdp1*xdp1
    bvalx3=dxxx+3.*dxx*xd1+3.*this%geom%dx* &
         &xd1*xd1-3.*xd1*xd1*xd1
    bvalx4=xd1*xd1*xd1
    bvaly1=ydp1*ydp1*ydp1
    bvaly2=dyyy+3.*ydp1*(dyy+ydp1*(this%geom%dy-ydp1))          
    bvaly3=dyyy+3.*yd1*(dyy+yd1*(this%geom%dy-yd1))
    bvaly4=yd1*yd1*yd1

    do j=1,this%geom%ny
       j1=mod(this%geom%ny+j-1+intaysdy,this%geom%ny)

       do i=1,this%geom%nx
          i1=mod(this%geom%nx+i-1+intaxsdx,this%geom%nx)

          fout(i,j) = dxxx6*dyyy6* ( &
               bvalx1*( this%coef(i1+1,j1+1)*bvaly1 &
               +this%coef(i1+1,j1+2)*bvaly2 &
               +this%coef(i1+1,j1+3)*bvaly3 &
               +this%coef(i1+1,j1+4)*bvaly4) &
               + bvalx2* (this%coef(i1+2,j1+1)*bvaly1 &
               +this%coef(i1+2,j1+2)*bvaly2 &
               +this%coef(i1+2,j1+3)*bvaly3 &
               +this%coef(i1+2,j1+4)*bvaly4) &
               + bvalx3* (this%coef(i1+3,j1+1)*bvaly1 &
               +this%coef(i1+3,j1+2)*bvaly2 &
               +this%coef(i1+3,j1+3)*bvaly3 &
               +this%coef(i1+3,j1+4)*bvaly4) &
               + bvalx4* (this%coef(i1+4,j1+1)*bvaly1 &
               +this%coef(i1+4,j1+2)*bvaly2 &
               +this%coef(i1+4,j1+3)*bvaly3 &
               +this%coef(i1+4,j1+4)*bvaly4))
       end do
    end do

  end subroutine evaldep
end module splinepp_class
