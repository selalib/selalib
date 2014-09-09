module splinenn_class
  use used_precision
  use geometry_module
  implicit none

  type, public :: splinenn
     type (geometry) :: geom
     real(wp) :: a1x, a2x, a3x, a4x ! coef de la matrice 2x2 per.
     real(wp) :: a1y, a2y, a3y, a4y ! coef de la matrice 2x2 nat.
     real(wp), dimension(:), pointer :: axd, ayd ! termes diagonaux de ax
     real(wp), dimension(:), pointer :: axod, ayod ! termes sous-diag. de ax
     real(wp), dimension(:), pointer :: aym1gamma1, aym1gamma2
     real(wp), dimension(:), pointer :: axm1gamma1, axm1gamma2
     real(wp), dimension(:,:), pointer :: coef, bcoef ! coefficients des splines
  end type splinenn
  interface initialize
     module procedure new_splinenn
  end interface
  interface interpole
     module procedure interpole_splinenn, interpole_splinenndep
  end interface
  public :: initialize, interpole
contains
  subroutine new_splinenn(this,geom,iflag)
    type(splinenn), intent(out) :: this
    type(geometry), intent(in) :: geom  ! geometry of problem
    integer, intent(out)       :: iflag ! error flag
    ! local variables
    integer :: err ! error flag
    integer :: nx,ny ! dimensions
    integer :: nxp1, nxp2, nyp1, nyp2
    real(wp) :: aa1x, aa1y, aa2x, aa2y
    integer i,j, info
    ! initialisation des variables locales
    iflag = 0
    nx = geom%nx
    ny = geom%ny
    nxp1=nx+1
    nxp2=nx+2
    nyp1=ny+1
    nyp2=ny+2

    ! initialisation de geom
    this%geom = geom

    ! memory allocation
    allocate(this%axm1gamma1(nx), stat=err)
    if (err.ne.0) then
       iflag = 10
       return
    end if
    allocate(this%axm1gamma2(nx), stat=err)
    if (err.ne.0) then
       iflag = 11
       return
    end if
    allocate(this%aym1gamma1(ny), stat=err)
    if (err.ne.0) then
       iflag = 15
       return
    end if
    allocate(this%aym1gamma2(ny), stat=err)
    if (err.ne.0) then
       iflag = 16
       return
    end if
    allocate(this%axd(nx), stat=err)
    if (err.ne.0) then
       iflag = 20
       return
    end if
    allocate(this%axod(nx-1), stat=err)
    if (err.ne.0) then
       iflag = 30
       return
    end if
    allocate(this%ayd(ny), stat=err)
    if (err.ne.0) then
       iflag = 40
       return
    end if
    allocate(this%ayod(ny-1), stat=err)
    if (err.ne.0) then
       iflag = 45
       return
    end if
    allocate(this%coef(nxp2,nyp2), stat=err)
    if (err.ne.0) then
       iflag = 50
       return
    end if
    allocate(this%bcoef(nyp2,nxp2), stat=err)
    if (err.ne.0) then
       iflag = 60
       return
    end if

    ! factorize matrices  Ax and Ay
    this%axd = 4_wp
    this%axod = 1_wp
#ifdef _CRAY
    call spttrf(nx,this%axd,this%axod,err)
#else
    call dpttrf(nx,this%axd,this%axod,err)
#endif
    if (err.ne.0) then
       iflag = 70
       return
    end if
    this%ayd = 4_wp
    this%ayod = 1_wp
#ifdef _CRAY
    call spttrf(ny,this%ayd,this%ayod,err)
#else
    call dpttrf(ny,this%ayd,this%ayod,err)
#endif
    if (err.ne.0) then
       iflag = 80
       return
    end if
    ! compute Ax-1.gamma
    this%axm1gamma1 = 0_wp
    this%axm1gamma2 = 0_wp
    this%axm1gamma2(1) = 1_wp
    this%axm1gamma1(nx) = 1_wp
#ifdef _CRAY
    call spttrs(nx,1,this%axd,this%axod,this%axm1gamma1,nx,err)
    call spttrs(nx,1,this%axd,this%axod,this%axm1gamma2,nx,err)
#else
    call dpttrs(nx,1,this%axd,this%axod,this%axm1gamma1,nx,err)
    call dpttrs(nx,1,this%axd,this%axod,this%axm1gamma2,nx,err)
#endif

    if (err.ne.0) then
       iflag = 90
       return
    end if
    ! compute Ay-1.gamma
    this%aym1gamma1 = 0_wp
    this%aym1gamma2 = 0_wp
    this%aym1gamma2(1) = 1_wp
    this%aym1gamma1(ny) = 1_wp
#ifdef _CRAY
    call spttrs(ny,1,this%ayd,this%ayod,this%aym1gamma1,ny,err)
    call spttrs(ny,1,this%ayd,this%ayod,this%aym1gamma2,ny,err)
#else
    call dpttrs(ny,1,this%ayd,this%ayod,this%aym1gamma1,ny,err)
    call dpttrs(ny,1,this%ayd,this%ayod,this%aym1gamma2,ny,err)
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
    this%a1x = aa2x*(1. - this%axm1gamma1(nx-1) + 2*this%axm1gamma1(nx))
    this%a2x = aa2x*(-this%axm1gamma2(nx-1) + 2*this%axm1gamma2(nx))
    this%a3x = aa2x*(2*this%axm1gamma1(1)-this%axm1gamma1(2))
    this%a4x = aa2x*(1. + 2*this%axm1gamma2(1) - this%axm1gamma2(2))
    ! assemblage de la matrice 2x2 pour spline naturels (direction Oy)
    this%a1y = aa2y*(1. - this%aym1gamma1(ny-1) + 2*this%aym1gamma1(ny))
    this%a2y = aa2y*(-this%aym1gamma2(ny-1) + 2*this%aym1gamma2(ny))
    this%a3y = aa2y*(2*this%aym1gamma1(1)-this%aym1gamma1(2))
    this%a4y = aa2y*(1. + 2*this%aym1gamma2(1) - this%aym1gamma2(2))
  end subroutine new_splinenn

  subroutine interpole_splinenn(this,fin,fout,x,y) 
    type(splinenn), intent(inout) :: this
    ! fin contient les valeurs de la fonction dans la grille precedente
    real(wp), dimension(:,:), intent(in) :: fin
    ! fout est destine a contenir la nouvelle valeur de f
    real(wp), dimension(:,:), intent(out):: fout
    ! dans x et y on trouve les points dans les quels on veut 
    ! evaluer la spline.
    real(wp), dimension(:,:), intent(in) :: x, y 
    ! dans fout on trouve en sortie les valeurs de f(i,j) 
    ! dans les points x(i),y(i).
    integer :: iflag ! error flag
    ! variables locales
    integer i,j, ierr

    ! initialisation des variables locales
    iflag = 0


    call nat_x(this,fin,ierr)
    if (ierr.ne.0) then
       iflag = 10
       return
    end if

    call nat_y(this,ierr)
    if (ierr.ne.0) then
       iflag = 20
       return
    end if

    call evaltab(this,x,y,fout)     

    return

  end subroutine interpole_splinenn

  subroutine interpole_splinenndep(this,f,depx,depy,aff) 
    !----------------------------------------------------------------
    ! interpolation par spline periodique dans les deux directions.
    ! Les points d'interpolation sont definis grace a depx et depy
    ! qui definissent le deplacement par rapport au maillage.
    !----------------------------------------------------------------
    type(splinenn), intent(inout) :: this
    ! f contient les valeurs de la fonction de distribution
    real(wp), dimension(:,:), intent(inout) :: f
    ! dans depx et depy on trouve les deplacements par rapport au maillage
    ! des points dans les quels on veut evaluer la spline.
    real(wp), intent(in) :: depx, depy 
    logical :: aff
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    integer ierr, l_a, l_b
    real(wp) durat, vtime(1:4)

    !if (aff) then 
    !   call clck_temps(l_a)
    !end if
    call nat_x(this,f,ierr)
    if (ierr.ne.0) then
       iflag = 10
       return
    end if
    !if (aff) then 
    !   call clck_temps(l_b)
    !   call clck_diff(l_a,l_b,vtime(1))
    !   call clck_temps(l_a)
    !end if

    call nat_y(this,ierr)
    
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
    !   write(*,'(A,3(1x,3E14.5))') "splinenn ",vtime(1:3)
    !end if

  end subroutine interpole_splinenndep

  !
  ! calcul des "natural splines"
  ! 
  subroutine nat_x(this,gtau,iflag)
    type(splinenn), intent(inout) :: this     ! objet de type spline
    real(wp), dimension(:,:), intent(in) :: gtau ! valeur de la fonction 
    ! aux points du maillage
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i,j ! indices de boucle
    integer nx, ny, nxp1, nxp2  !  nx+1, nx+2
    integer :: err ! error flag

    real(wp) :: axm1f(this%geom%nx,this%geom%ny)
    real(wp) :: det, gamma1, gamma2, coef1, coefnp2

    ! initialisations
    iflag =0
    nx=this%geom%nx
    ny=this%geom%ny
    nxp2=nx+2
    nxp1=nx+1
    det=this%a1x*this%a4x - this%a2x*this%a3x

    ! Calcul de Ax^-1 f
    ! assemblage du membre de droite pour le calcul de Ax^-1 f
    do j=1,ny
       do i=1,nx
          axm1f(i,j) = 6*gtau(i,j)
       end do
    end do
#ifdef _CRAY
    call spttrs(nx,ny,this%axd,this%axod,axm1f,nx,err)
#else
    call dpttrs(nx,ny,this%axd,this%axod,axm1f,nx,err)
#endif
    if (err.ne.0) then
       iflag = 10
       return
    end if
    !print*,'axmf1',axm1f
    do  j=1,ny
       ! assemblage du second membre du systeme 2x2 
       gamma1 =  (6.0/(this%geom%dx)**2)*(-axm1f(this%geom%nx-1,j) &
            + 2*axm1f(this%geom%nx,j))
       gamma2 = (6.0/(this%geom%dx)**2)*(2*axm1f(1,j) - axm1f(2,j))

       coefnp2 = (gamma1*this%a4x - gamma2*this%a2x)/det
       coef1 = (gamma2*this%a1x - gamma1*this%a3x)/det
       this%bcoef(j,nxp2)=coefnp2
       this%bcoef(j,1)=coef1

       do  i=2,nxp1
          this%bcoef(j,i)= axm1f(i-1,j) &
               - this%axm1gamma1(i-1)*coefnp2 &
               - this%axm1gamma2(i-1)*coef1
       end do
    end do
  end subroutine nat_x

  subroutine nat_y(this,iflag)
    type(splinenn), intent(inout) :: this     ! objet de type spline
    integer, intent(out) :: iflag    ! indicateur d erreur
    ! variables locales
    integer i,j          ! indices de boucle
    integer :: ny, nxp2,nyp2 ! nx+2, ny+2
    integer :: err       ! indicateur d erreur

    real(wp) :: aym1f(this%geom%ny,this%geom%nx+2)
    real(wp) :: det, gamma1, gamma2, coef1, coefnp2

    ! initialisations
    iflag =0
    ny = this%geom%ny
    nxp2=this%geom%nx+2
    nyp2=this%geom%ny+2
    det=this%a1y*this%a4y - this%a2y*this%a3y

    ! calcul de coef par resolution de nxp2 systemes lineaires.

    ! Calcul de Ay^-1 f
    ! assemblage du membre de droite pour le calcul de Ay^-1 f
    do  i=1,nxp2
       do j=1,this%geom%ny
          aym1f(j,i) = 6.*this%bcoef(j,i)
       end do
    end do
#ifdef _CRAY
    call spttrs(ny,nxp2,this%ayd,this%ayod,aym1f,ny,err)
#else
    call dpttrs(ny,nxp2,this%ayd,this%ayod,aym1f,ny,err)
#endif
    if (err.ne.0) then
       iflag = 10
       return
    end if
    ! resolution du syteme lineaire 2x2
    do i=1,nxp2
       gamma1 =  (6.0/(this%geom%dy)**2)*(-aym1f(this%geom%ny-1,i) &
            + 2*aym1f(this%geom%ny,i))
       gamma2 = (6.0/(this%geom%dy)**2)*(2*aym1f(1,i) - aym1f(2,i))

       coefnp2 = (gamma1*this%a4y - gamma2*this%a2y)/det
       coef1 = (gamma2*this%a1y - gamma1*this%a3y)/det
       this%coef(i,nyp2) = coefnp2
       this%coef(i,1) = coef1
       do  j=2,this%geom%ny + 1
          this%coef(i,j)= aym1f(j-1,i)                   &
               - this%aym1gamma1(j-1)*coefnp2 &
               - this%aym1gamma2(j-1)*coef1
       end do
    end do
  end subroutine nat_y
  !

  !
  ! tsvaleur
  !
  subroutine evaltab(this,xd,yd,fout)
    type(splinenn) :: this
    real(wp), dimension(:,:) :: xd, yd ! coordonnees du point ou les valeurs sont calculees

    real(wp) :: sval,idx,idy   ! valeur de la fonction au point (xd,yd)
    real(wp), dimension(:,:) :: fout

    real(wp) bvalx1,bvalx2,bvalx3,bvalx4,bvaly1,bvaly2,bvaly3, &
         &bvaly4,a1,dxx,dxxx,dxxx6,dyy,dyyy,dyyy6,xd1,xdp1,yd1,ydp1
    real(wp) :: sval1, sval2, sval3, sval4
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
    do j=2,this%geom%ny-1
       do i=2,this%geom%nx-1

          i1=(xd(i,j)-this%geom%x0)*idx
          j1=(yd(i,j)-this%geom%y0)*idy

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

    fout(1,:)=0
    fout(this%geom%nx,:)=0 

    fout(:,1)=0
    fout(:,this%geom%ny)=0 
  end subroutine evaltab

  subroutine evaldep(this,alphax,alphay,fout)
    type(splinenn) :: this
    real(wp) :: alphax,alphay ! deplacement par rapport aux maillages des points ou la spline est evaluee
    real(wp) :: sval   ! valeur de la fonction au point d'evaluation
    real(wp), dimension(:,:) :: fout

    real(wp) bvalx1,bvalx2,bvalx3,bvalx4,bvaly1,bvaly2,bvaly3, &
         &bvaly4,dxx,dxxx,dxxx6,dyy,dyyy,dyyy6,xd1,xdp1,yd1,ydp1
    real(wp) :: sval1, sval2, sval3, sval4
    integer :: intaxsdx, intaysdy
    integer i1,j1,i,j,ideb,ifin,jdeb,jfin
    !
    dxx=this%geom%dx*this%geom%dx
    dxxx=dxx*this%geom%dx
    dxxx6=1./(6.*dxxx)
    !
    !
    dyy=this%geom%dy*this%geom%dy
    dyyy=dyy*this%geom%dy
    dyyy6=1./(6.*dyyy)

!!$    if (abs(alphax).gt.this%geom%dx) then
!!$       print*, 'deplacement en x trop grand, alphax=',alphax
!!$       print*, 'dx=',this%geom%dx
!!$       stop
!!$    end if
!!$    if (abs(alphay).gt.this%geom%dy) then
!!$       print*,'deplacement en y trop grand, alphay=',alphay
!!$       print*, 'dy=',this%geom%dy
!!$       stop
!!$    end if

    if (alphax.gt.0) then
       intaxsdx=int(-alphax/this%geom%dx+epsilon)-1 
!       i1=this%geom%nx-2
!       intaxsdx=-1
       ideb = int(alphax/this%geom%dx)+2
       ifin = this%geom%nx - 1
    else
!       intaxsdx=0
       intaxsdx=int(-alphax/this%geom%dx)
!       i1=this%geom%nx-1
       ideb = 2
       ifin = -int(-alphax/this%geom%dx)+this%geom%nx-1
    end if

    if (alphay.gt.0) then
!       intaysdy=-1
       intaysdy=int(-alphay/this%geom%dy+epsilon)-1
!       j1=-1
       jdeb = int(alphay/this%geom%dy)+2
       jfin = this%geom%ny - 1
    else
!       intaysdy=0
       intaysdy=int(-alphay/this%geom%dy)
!       j1=0
        jdeb = 2
        jfin = -int(-alphay/this%geom%dy)+this%geom%ny-1
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
!print*,'eval ',xd1,yd1,xdp1,ydp1,intaxsdx,intaysdy
!print*,'eval ',this%geom%dx,this%geom%dy
!print*,'eval ',bvalx1,bvalx2,bvalx3,bvalx4
!print*,'eval ',bvaly1,bvaly2,bvaly3,bvaly4
!    do j=2,this%geom%ny-1
    do j=jdeb,jfin
!       j1=j1+1
!       i1=i1-(this%geom%nx-2) ! remise de i1 a sa valeur de depart
       j1=j-1+intaysdy
!j1=mod(this%geom%ny+j-2+intaysdy,this%geom%ny-1)
!       do i=2,this%geom%nx-1
       do i=ideb,ifin
!          i1=i1+1
          i1=i-1+intaxsdx
! i1=mod(this%geom%nx+i-2+intaxsdx,this%geom%nx-1)

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

    fout(1:ideb-1,:)=0
    fout(ifin+1:this%geom%nx,:)=0

    fout(:,1:jdeb-1)=0
    fout(:,jfin+1:this%geom%ny)=0

  end subroutine evaldep
end module splinenn_class
