module csl2dpp_class

  use used_precision
  use geometry_module
  !use clock
  implicit none
  private
  public :: new, interpole
  type, public :: csl2dpp
     type (geometry) :: geom
  end type csl2dpp
  interface new
     module procedure new_csl2dpp
  end interface
  interface interpole
     module procedure interpole_csl2dpp,interpole_csl2dppdep
  end interface
contains

  subroutine new_csl2dpp(this,geom,iflag)
    type(csl2dpp), intent(out) :: this
    type(geometry), intent(in) :: geom  ! geometry of problem
    integer, intent(out)       :: iflag ! error flag

    ! initialisation de geom
    this%geom = geom


  end subroutine new_csl2dpp

  subroutine interpole_csl2dpp(this,fin,fout,x,y) 
    type(csl2dpp), intent(inout) :: this
    ! fin contient les valeurs de la fonction dans la grille précedente
    real(wp), dimension(:,:), intent(in) :: fin
    ! fout est destiné à contenir la nouvelle valeur de f
    real(wp), dimension(:,:), intent(out):: fout
    ! dans x et y on trouve les points auxquels on veut 
    ! évaluer la spline.
    real(wp), dimension(:,:), intent(in) :: x, y 
    ! dans fout, on trouve en sortie les valeurs de f(i,j) 
    ! aux points x(i),y(i).
    ! indicateur d'erreur
    integer :: iflag ! error flag
    ! variables locales
    integer ierr



  end subroutine interpole_csl2dpp

  subroutine interpole_csl2dppdep(this,f,depx,depy,aff) 
    !----------------------------------------------------------------
    ! interpolation par spline périodique dans les deux directions.
    ! Les points d'interpolation sont définis gr^ace à depx et depy
    ! qui définissent le déplacement par rapport au maillage.
    !----------------------------------------------------------------
    type(csl2dpp), intent(inout) :: this
    ! f contient les valeurs de la fonction de distribution
    real(wp), dimension(:,:), intent(inout) :: f
    ! dans depx et depy on trouve les deplacements par rapport au maillage
    ! des points auxquels on veut évaluer la spline.
    real(wp), intent(in) :: depx, depy!, dt 
    ! indicateur d'erreur
    integer::i,j
    logical :: aff
    integer::timecase(0:1),interp_case,ppm_order
    ! variables locales
    real(wp),dimension(:,:,:),allocatable::carac
    real(wp),dimension(:,:),allocatable::buf2d,E0,E1
    real(wp)::dom(0:1,0:1)

    timecase(0)=3 !computation of the characteristics: 1 for Euler, 3 for symplectic Verlet
    timecase(1)=2 !number of steps fixed point algo (symplectic Verlet case)
    interp_case=2 !1:Lauritzen 2:LAG3 3:PPM CD 4:LAG(2d+1)
    ppm_order=4   !if interp_case=3: PPM0, PPM1 or PPM2; if interp_case=4: ppm_order=d for LAG(2d+1)


    allocate(E0(1:this%geom%nx,1:this%geom%ny),E1(1:this%geom%nx,1:this%geom%ny))
    allocate(carac(2,-1:this%geom%nx,-1:this%geom%ny))
    allocate(buf2d(0:this%geom%nx-1,0:this%geom%ny-1))

    dom(0,0)=this%geom%x0
    dom(0,1)=this%geom%y0
    dom(1,0)=this%geom%nx*this%geom%dx
    dom(1,1)=this%geom%ny*this%geom%dy

    do j=1,this%geom%ny
      do i=1,this%geom%nx
        E0(i,j)=depx
        E1(i,j)=depy
      enddo
    enddo

    call compute_carac_per_per(dom,E0,E1,this%geom%nx,this%geom%ny,timecase,carac)

    call advect2d_CSL(dom,f,this%geom%nx,this%geom%ny,buf2d,interp_case,carac,ppm_order)

    deallocate(E0,E1,carac,buf2d)

  end subroutine interpole_csl2dppdep

  subroutine advect2d_CSL(dom,f,N0,N1,buf2d,interp_case,carac,ppm_order) !conservative 2d remapping algorithm
    real(wp),dimension(0:1,0:1),intent(in)::dom
    integer,intent(in)::N0,N1,interp_case,ppm_order
    real(wp),dimension(0:N0-1,0:N1-1)::buf2d
    real(wp),dimension(1:N0,1:N1)::f
    real(wp),dimension(2,-1:N0,-1:N1)::carac
    real(wp)::xx(4),yy(4),xA,yA,xB,yB,res,xx0,yy0,x,xxn,y,yyn,dx,dy
    real(wp)::xA_loc,yA_loc,xB_loc,yB_loc
    integer::i,j,ii(4),jj(4),im1,jm1,i0,j0,i1,j1,s,k,sx,sy,minfl,maxfl,ix,iii,iy,iiii,jjjj
    real(wp),dimension(:),allocatable::intx,inty
    integer,dimension(:,:),allocatable::tnbr
    integer,dimension(:,:,:),allocatable::cell
    real(wp),dimension(:,:),allocatable::tt,tcell,dir,aretesh,aretesv,sommets,aretesvg,aretesvd,areteshb,areteshh,&
    sommetsbg,sommetsbd,sommetshg,sommetshd
    real(wp),dimension(:,:,:,:),allocatable::tpts
    integer::nbx(4),nby(4),nbmax,dirx,diry,ell,ell1,imin,jmin,imax,jmax,i0_loc,j0_loc
    real(wp)::xx1,xx2,yy1,yy2,w00,w10,w01,w20,w02,w11,w21,&
		w12,w22,c00,c10,c01,c20,c02,c11,c12,c21,c22,&
		fij,fim1jm1,fim1j,fim1jp1,fijm1,fijp1,fip1jm1,fip1j,fip1jp1,xxx,yyy
    integer::im2,ib,ip1,ip2,jm2,jb,jp1,jp2   

    nbmax=30
    !print*,f(1,33)
    !return
    allocate(tt(nbmax,2),cell(2,nbmax,4),tcell(nbmax,4),intx(0:nbmax),inty(0:nbmax),dir(nbmax,2))
    do j=0,N1-1
      do i=0,N0-1
        buf2d(i,j)=f(i+1,j+1)
      enddo
    enddo
    dx=dom(1,0)/real(N0,wp)
    dy=dom(1,1)/real(N1,wp)

    if ((interp_case==3) .or. (interp_case==6) .or. (interp_case==7)) then
      allocate(aretesh(0:N0-1,0:N1-1),aretesv(0:N0-1,0:N1-1),sommets(0:N0-1,0:N1-1))
      call aux(N0,N1,buf2d,aretesh,aretesv,sommets,ppm_order,carac,dom)
    endif
    if (interp_case==4) then
      allocate(aretesvg(0:N0-1,0:N1-1),aretesvd(0:N0-1,0:N1-1),areteshb(0:N0-1,0:N1-1),areteshh(0:N0-1,0:N1-1),&
      sommetsbg(0:N0-1,0:N1-1),sommetsbd(0:N0-1,0:N1-1),sommetshg(0:N0-1,0:N1-1),sommetshd(0:N0-1,0:N1-1))
      call aux2(N0,N1,buf2d,areteshb,areteshh,aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,ppm_order,carac,dom)
    endif

    do j=0,N1-1
      do i=0,N0-1

	f(i+1,j+1)=0._wp

        im1=modulo(i-1,N0)
	i1=modulo(i+1,N0)
	jm1=modulo(j-1,N1)
	j1=modulo(j+1,N1)
        
	!computation of the feet of the characteristics
        xx(1)=0.25_wp*(carac(1,i-1,j-1)+carac(1,i,j-1)+carac(1,i,j)+carac(1,i-1,j))
        xx(2)=0.25_wp*(carac(1,i,j-1)+carac(1,i+1,j-1)+carac(1,i+1,j)+carac(1,i,j))
        xx(3)=0.25_wp*(carac(1,i,j)+carac(1,i+1,j)+carac(1,i+1,j+1)+carac(1,i,j+1))
        xx(4)=0.25_wp*(carac(1,i-1,j)+carac(1,i,j)+carac(1,i,j+1)+carac(1,i-1,j+1))

        yy(1)=0.25_wp*(carac(2,i-1,j-1)+carac(2,i,j-1)+carac(2,i,j)+carac(2,i-1,j))
        yy(2)=0.25_wp*(carac(2,i,j-1)+carac(2,i+1,j-1)+carac(2,i+1,j)+carac(2,i,j))
        yy(3)=0.25_wp*(carac(2,i,j)+carac(2,i+1,j)+carac(2,i+1,j+1)+carac(2,i,j+1))
        yy(4)=0.25_wp*(carac(2,i-1,j)+carac(2,i,j)+carac(2,i,j+1)+carac(2,i-1,j+1))      


        !normalization
	xx(1)=(xx(1)-dom(0,0))/dom(1,0)*real(N0,wp)
	xx(2)=(xx(2)-dom(0,0))/dom(1,0)*real(N0,wp)
	xx(3)=(xx(3)-dom(0,0))/dom(1,0)*real(N0,wp)
	xx(4)=(xx(4)-dom(0,0))/dom(1,0)*real(N0,wp)
	
	yy(1)=(yy(1)-dom(0,1))/dom(1,1)*real(N1,wp)
	yy(2)=(yy(2)-dom(0,1))/dom(1,1)*real(N1,wp)
	yy(3)=(yy(3)-dom(0,1))/dom(1,1)*real(N1,wp)
	yy(4)=(yy(4)-dom(0,1))/dom(1,1)*real(N1,wp)

        xx=xx+0.5_wp
        yy=yy+0.5_wp 

        
        ii(1)=floor(xx(1))
        ii(2)=floor(xx(2))
        ii(3)=floor(xx(3))
        ii(4)=floor(xx(4))
	
        jj(1)=floor(yy(1))
        jj(2)=floor(yy(2))
        jj(3)=floor(yy(3))
        jj(4)=floor(yy(4))
	
	imin=min(ii(1),ii(2),ii(3),ii(4))
	jmin=min(jj(1),jj(2),jj(3),jj(4))
	imax=max(ii(1),ii(2),ii(3),ii(4))
	jmax=max(jj(1),jj(2),jj(3),jj(4))

!        allocate(tnbr(imin:imax,jmin:jmax))
!        allocate(tpts(imin:imax,jmin:jmax,100,2))

!        tnbr=0    
	!Computation of external edges

    do ell=1,4

	ell1=ell+1;if(ell1==5)ell1=1  
	xA=xx(ell);yA=yy(ell);xB=xx(ell1);yB=yy(ell1)
	i0=ii(ell);j0=jj(ell);i1=ii(ell1);j1=jj(ell1)

	s=1;
	if(i0<i1)then
	  do k=i0+1,i1
	    tt(s,1)=(real(k,wp)-xA)/(xB-xA)
	    s=s+1
	  enddo	  
	  dirx=1
	endif
	if(i0>i1)then
	  do k=i0,i1+1,-1
	    tt(s,1)=(real(k,wp)-xA)/(xB-xA)
	    s=s+1
	  enddo
	  dirx=-1
	endif
        nbx(ell)=s-1;
	s=1;
	if(j0<j1)then
	  do k=j0+1,j1
	    tt(s,2)=(real(k,wp)-yA)/(yB-yA)
	    s=s+1
	  enddo
	  diry=1
	endif
	if(j0>j1)then
	  do k=j0,j1+1,-1
	    tt(s,2)=(real(k,wp)-yA)/(yB-yA)
	    s=s+1
	  enddo
	  diry=-1
	endif
	nby(ell)=s-1
	
	cell(1,1,ell)=i0
	cell(2,1,ell)=j0
	tcell(1,ell)=0._wp
	sx=1;sy=1
	s=1

	do while((sx<=nbx(ell)).and.(sy<=nby(ell)))
	  if(tt(sx,1)<tt(sy,2))then
	    s=s+1
	    cell(1,s,ell)=cell(1,s-1,ell)+dirx
	    cell(2,s,ell)=cell(2,s-1,ell)
	    tcell(s,ell)=tt(sx,1)
	    intx(2*cell(1,s-1,ell)+(dirx-1)/2-2*imin)=yA+tt(sx,1)*(yB-yA)
	    dir(s,1)=1
	    dir(s,2)=dirx  
	    sx=sx+1	
	  else
	    s=s+1
	    cell(1,s,ell)=cell(1,s-1,ell)
	    cell(2,s,ell)=cell(2,s-1,ell)+diry
	    tcell(s,ell)=tt(sy,2)
	    inty(2*cell(2,s-1,ell)+(diry-1)/2-2*jmin)=xA+tt(sy,2)*(xB-xA)
	    dir(s,1)=2
	    dir(s,2)=diry	  
	    sy=sy+1	    
	  endif
	enddo
	do while(sx<=nbx(ell))
	  s=s+1
	  cell(1,s,ell)=cell(1,s-1,ell)+dirx
	  cell(2,s,ell)=cell(2,s-1,ell)
	  tcell(s,ell)=tt(sx,1)
	  intx(2*cell(1,s-1,ell)+(dirx-1)/2-2*imin)=yA+tt(sx,1)*(yB-yA)
	  dir(s,1)=1
	  dir(s,2)=dirx
	  sx=sx+1
	enddo  
	do while(sy<=nby(ell))
	  s=s+1
	  cell(1,s,ell)=cell(1,s-1,ell)
	  cell(2,s,ell)=cell(2,s-1,ell)+diry
	  tcell(s,ell)=tt(sy,2)
	  inty(2*cell(2,s-1,ell)+(diry-1)/2-2*jmin)=xA+tt(sy,2)*(xB-xA)
	  dir(s,1)=2
	  dir(s,2)=diry	  
	  sy=sy+1
	enddo

!Computation of extern edges

if ((ell==1) .or. (ell==4) .or. ((ell==2) .and. (i==N0-1)) .or. ((ell==3) .and. (j==N1-1))) then

        xB_loc=xx(ell)
        yB_loc=yy(ell)
        do k=2,nbx(ell)+nby(ell)+2
          xA_loc=xB_loc
          yA_loc=yB_loc
          i0_loc=cell(1,k-1,ell)
          j0_loc=cell(2,k-1,ell)
          if(k==nbx(ell)+nby(ell)+2)then
            xB_loc=xx(ell1)
            yB_loc=yy(ell1)
          else
            if(dir(k,1)==1)then
              xB_loc=real(cell(1,k,ell)+(1-dir(k,2))/2,wp)
              yB_loc=yy(ell)+tcell(k,ell)*(yy(ell1)-yy(ell))        	  
            else
              xB_loc=xx(ell)+tcell(k,ell)*(xx(ell1)-xx(ell))					
              yB_loc=real(cell(2,k,ell)+(1-dir(k,2))/2,wp)       	  
            endif
          endif

          call calcule_coeff(N0,N1,buf2d,i0_loc,j0_loc,xA_loc,yA_loc,xB_loc,yB_loc,res,aretesh,aretesv,sommets,areteshb,areteshh,&
          aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,interp_case)

	  f(i+1,j+1)=f(i+1,j+1)+res

          if ((i-1>=0) .and. (ell==4)) then
            f(i,j+1)=f(i,j+1)-res
          endif
          if ((j-1>=0) .and. (ell==1)) then
            f(i+1,j)=f(i+1,j)-res
          endif
        enddo
endif
        
      end do  	

!Computation of vertical intern edges
	
	do ell=0,imax-imin-1
	minfl=min(floor(intx(2*ell)),floor(intx(2*ell+1)))
	maxfl=max(floor(intx(2*ell)),floor(intx(2*ell+1)))

        i0_loc=imin+ell
        yB_loc=min(intx(2*ell),intx(2*ell+1))
        do k=0,maxfl-minfl
          yA_loc=yB_loc
          j0_loc=minfl+k  
          if(k==maxfl-minfl)then
            yB_loc=max(intx(2*ell),intx(2*ell+1))
          else
            yB_loc=real(minfl+k+1,wp)
          endif
          call calcule_coeffv(N0,N1,buf2d,i0_loc,j0_loc,yA_loc,yB_loc,res,aretesh,aretesv,sommets,areteshb,areteshh,&
	  aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,interp_case)
          f(i+1,j+1)=f(i+1,j+1)+res

	  enddo
	enddo

!Computation of horizontal intern edges
	
	do ell=0,jmax-jmin-1
	minfl=min(floor(inty(2*ell)),floor(inty(2*ell+1)))
	maxfl=max(floor(inty(2*ell)),floor(inty(2*ell+1)))

        j0_loc=jmin+ell
        xA_loc=min(inty(2*ell),inty(2*ell+1))
        do k=0,maxfl-minfl
          xB_loc=xA_loc
          i0_loc=minfl+k  
          if(k==maxfl-minfl)then
            xA_loc=max(inty(2*ell),inty(2*ell+1))
          else
            xA_loc=real(minfl+k+1,wp)
          endif
          
          call calcule_coeffh(N0,N1,buf2d,i0_loc,j0_loc,xA_loc,xB_loc,res,aretesh,aretesv,sommets,areteshb,areteshh,&
	  aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,interp_case)
	  f(i+1,j+1)=f(i+1,j+1)+res

	  enddo
	enddo

!        deallocate(tnbr)
!        deallocate(tpts)

      enddo

    enddo
   
!#ifdef DIAG_TIME
!      f=buf2d
!#endif
 
    deallocate(tt,cell,tcell,intx,inty,dir)
    if ((interp_case==3) .or. (interp_case==6) .or. (interp_case==7)) then
      deallocate(aretesh,aretesv,sommets)
    endif
    if (interp_case==4) then
      deallocate(areteshb,areteshh,aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd)
    endif

!print*,f(1,:),'suite',buf2d(0,:)
!stop
    
  end subroutine advect2d_CSL


	subroutine aux(N0,N1,f,aretesh,aretesv,sommets,ordre,carac,dom) !used in PPM case

		integer,intent(in)::N0,N1,ordre
                real(wp),dimension(0:1,0:1),intent(in)::dom
		real(wp),dimension(0:N0-1,0:N1-1),intent(in)::f
                real(wp),dimension(2,-1:N0+1,-1:N1+1),intent(in)::carac
		integer::i,j,im3,im2,im1,ib,ip1,ip2,jm3,jm2,jm1,jb,jp1,jp2
		real(wp),dimension(0:N0-1,0:N1-1),intent(inout)::aretesh,aretesv,sommets
!		print*,N0,N1,ordre
!stop

		do i=0,N0-1
			do j=0,N1-1
				im3=modulo(i-3,N0)
				im2=modulo(i-2,N0)
				im1=modulo(i-1,N0)
				ib=modulo(i,N0)
				ip1=modulo(i+1,N0)
				ip2=modulo(i+2,N0)
				jm3=modulo(j-3,N1)
				jm2=modulo(j-2,N1)
				jm1=modulo(j-1,N1)
				jb=modulo(j,N1)
				jp1=modulo(j+1,N1)
				jp2=modulo(j+2,N1)

				if (ordre==1) then !PPM1
					aretesv(i,j)=7._wp/12._wp*(f(im1,jb)+f(ib,jb)) &
					-1._wp/12._wp*(f(im2,jb)+f(ip1,jb))
					aretesh(i,j)=7._wp/12._wp*(f(ib,jm1)+f(ib,jb)) &	
					-1._wp/12._wp*(f(ib,jm2)+f(ib,jp1))
				else if (ordre==2) then !PPM2
					aretesv(i,j)=1._wp/60._wp*(f(ip2,jb)+f(im3,jb)) &
					-8._wp/60._wp*(f(ip1,jb)+f(im2,jb)) &
					+37._wp/60._wp*(f(ib,jb)+f(im1,jb))
					aretesh(i,j)=1._wp/60._wp*(f(ib,jp2)+f(ib,jm3)) &
					-8._wp/60._wp*(f(ib,jp1)+f(ib,jm2)) &
					+37._wp/60._wp*(f(ib,jb)+f(ib,jm1))
				else if (ordre==0) then !PPM0
					aretesv(i,j)=1._wp/2._wp*(f(ib,jb)+f(im1,jb))
					aretesh(i,j)=1._wp/2._wp*(f(ib,jb)+f(ib,jm1))
				end if
			end do
		end do	

		do i=0,N0-1
			do j=0,N1-1
				im3=modulo(i-3,N0)
				im2=modulo(i-2,N0)
				im1=modulo(i-1,N0)
				ib=modulo(i,N0)
				ip1=modulo(i+1,N0)
				ip2=modulo(i+2,N0)
				jm3=modulo(j-3,N1)
				jm2=modulo(j-2,N1)
				jm1=modulo(j-1,N1)
				jb=modulo(j,N1)
				jp1=modulo(j+1,N1)
				jp2=modulo(j+2,N1)

				if (ordre==1) then !PPM1
					sommets(i,j)=7._wp/12._wp*(aretesv(ib,jm1)+aretesv(ib,jb)) &
					-1._wp/12._wp*(aretesv(ib,jm2)+aretesv(ib,jp1))
				else if (ordre==2) then !PPM2
					sommets(i,j)=1._wp/60._wp*(aretesv(ib,jp2)+aretesv(ib,jm3)) &
					-8._wp/60._wp*(aretesv(ib,jp1)+aretesv(ib,jm2)) &
					+37._wp/60._wp*(aretesv(ib,jb)+aretesv(ib,jm1))
				else if (ordre==0) then !PPM0
					sommets(i,j)=1._wp/2._wp*(aretesv(ib,jb)+aretesv(ib,jm1))
				end if
			end do
		end do

	end subroutine aux


	subroutine aux2(N0,N1,f,areteshb,areteshh,aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,ordre,carac,dom) !used in PPM case

		integer,intent(in)::N0,N1,ordre
                real(wp),dimension(0:1,0:1),intent(in)::dom
		real(wp),dimension(0:N0-1,0:N1-1),intent(in)::f
                real(wp),dimension(2,-1:N0+1,-1:N1+1),intent(in)::carac
		integer::i,j,im3,im2,im1,ib,ip1,ip2,jm3,jm2,jm1,jb,jp1,jp2
		real(wp),dimension(0:N0-1,0:N1-1),intent(inout)::areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd
    real(wp) ::w(-ordre:ordre+1),tmp,ww(-ordre:ordre)
    
    !f2py intent(in)::buf,f
    integer::r,s,ii,d   
!		print*,N0,N1,ordre
!stop

    d=ordre
    r=-d
    s=d+1
    
    !maple code for generation of w
    !for k from r to -1 do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..k-1)*product((-j),j=k+1..-1)*product((-j),j=1..s):
    !od:
    !for k from 1 to s do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..-1)*product((-j),j=1..k-1)*product((-j),j=k+1..s):
    !od:
    !C[0]:=-add(C[k],k=r..-1)-add(C[k],k=1..s):
    
    do i=r,-1
      tmp=1._wp
      do j=r,i-1
        tmp=tmp*real(i-j,wp)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,wp)
      enddo
      tmp=1._wp/tmp
      do j=r,i-1
        tmp=tmp*real(-j,wp)
      enddo
      do j=i+1,-1
        tmp=tmp*real(-j,wp)
      enddo
      do j=1,s
        tmp=tmp*real(-j,wp)
      enddo
      w(i)=tmp      
    enddo

!    do i=r,-1
!      tmp=1._wp
!      !do j=r,i-1
!      !  tmp=tmp*real(i-j,wp)
!      !enddo
!      !do j=i+1,s
!      !  tmp=tmp*real(i-j,wp)
!      !enddo
!      !tmp=1._wp/tmp
!      do j=r,i-1 !-j/(i-j)=j/(j-i)=1/(1-i/j)
!        tmp=tmp*(1._wp-real(i,wp)/real(j,wp))
!      enddo
!      do j=i+1,-1
!        tmp=tmp*(1._wp-real(i,wp)/real(j,wp))
!      enddo
!      do j=1,s
!        tmp=tmp*(1._wp-real(i,wp)/real(j,wp))
!      enddo
!      tmp=tmp*real(i,wp)
!      w(i)=1._wp/tmp      
!    enddo



    do i=1,s
      tmp=1._wp
      do j=r,i-1
        tmp=tmp*real(i-j,wp)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,wp)
      enddo
      tmp=1._wp/tmp
      do j=r,-1
        tmp=tmp*real(-j,wp)
      enddo
      do j=1,i-1
        tmp=tmp*real(-j,wp)
      enddo
      do j=i+1,s
        tmp=tmp*real(-j,wp)
      enddo
      w(i)=tmp      
    enddo

    tmp=0._wp
    do i=r,-1
      tmp=tmp+w(i)
    enddo
    do i=1,s
      tmp=tmp+w(i)
    enddo
    w(0)=-tmp
    
    
    
    !print *,'w',w
    !do ii=r,s
    !  print *,ii,w(r+s-ii)
    !enddo
    
    !compute now ww
    !maple code
    !#for conservative formulation
    !tmp:=0:
    !for k from r to -1 do
    !tmp:=tmp+C[k]:
    !CC[k]:=-tmp:
    !od:
    !tmp:=0:
    !for k from s to 1 by -1 do
    !  tmp:=tmp+C[k]:
    !  CC[k-1]:=tmp:
    !od:
    !seq(CC[k],k=r..s-1);
    !evalf(%);

    tmp=0._wp
    do i=r,-1
      tmp=tmp+w(i)
      ww(i)=-tmp
    enddo
    tmp=0._wp
    do i=s,1,-1
      tmp=tmp+w(i)
      ww(i-1)=tmp
    enddo

    !print *,'ww',ww
    !do ii=r,s-1
    !  print *,ii,ww(r+s-1-ii)
    !enddo
    !stop
    do j=0,N1-1
      do i=0,N0-1
        tmp=0._wp
        do ii=r,s-1
          tmp=tmp+ww(r+s-1-ii)*f(modulo(i+ii-1,N0),j)
        enddo
        aretesvg(i,j)=tmp
        tmp=0._wp
        do ii=r,s-1
          tmp=tmp+ww(ii)*f(modulo(i+ii,N0),j)
        enddo
        aretesvd(i,j)=tmp
        tmp=0._wp
        do ii=r,s-1
          tmp=tmp+ww(ii)*f(i,modulo(j+ii,N1))
        enddo
        areteshh(i,j)=tmp
        tmp=0._wp
        do ii=r,s-1
          tmp=tmp+ww(r+s-1-ii)*f(i,modulo(j+ii-1,N1))
        enddo
        areteshb(i,j)=tmp
      enddo
    enddo

    do j=0,N1-1
      do i=0,N0-1
        tmp=0._wp
        do ii=r,s-1
          tmp=tmp+ww(ii)*areteshh(modulo(i+ii,N0),j)
        enddo
        sommetshd(i,j)=tmp
        tmp=0._wp
        do ii=r,s-1
          tmp=tmp+ww(r+s-1-ii)*areteshh(modulo(i+ii-1,N0),j)
        enddo
        sommetshg(i,j)=tmp
        tmp=0._wp
        do ii=r,s-1
          tmp=tmp+ww(ii)*areteshb(modulo(i+ii,N0),j)
        enddo
        sommetsbd(i,j)=tmp
        tmp=0._wp
        do ii=r,s-1
          tmp=tmp+ww(r+s-1-ii)*areteshb(modulo(i+ii-1,N0),j)
        enddo
        sommetsbg(i,j)=tmp
      enddo
    enddo



    if(0==1)then


		do i=0,N0-1
		  do j=0,N1-1
	  	    im3=modulo(i-3,N0)
		    im2=modulo(i-2,N0)
		    im1=modulo(i-1,N0)
		    ib=modulo(i,N0)
		    ip1=modulo(i+1,N0)
		    ip2=modulo(i+2,N0)
		    jm3=modulo(j-3,N1)
		    jm2=modulo(j-2,N1)
		    jm1=modulo(j-1,N1)
		    jb=modulo(j,N1)
		    jp1=modulo(j+1,N1)
		    jp2=modulo(j+2,N1)

		    if (ordre==1) then
                      aretesvg(i,j)=-1._wp/6._wp*f(im2,jb)+5._wp/6._wp*f(im1,jb) &
		      +1._wp/3._wp*f(ib,jb)
		      aretesvd(i,j)=1._wp/3._wp*f(im1,jb)+5._wp/6._wp*f(ib,jb) &
		      -1._wp/6._wp*f(ip1,jb)
		      areteshb(i,j)=-1._wp/6._wp*f(ib,jm2)+5._wp/6._wp*f(ib,jm1) &
		      +1._wp/3._wp*f(ib,jb)
		      areteshh(i,j)=1._wp/3._wp*f(ib,jm1)+5._wp/6._wp*f(ib,jb) &
		      -1._wp/6._wp*f(ib,jp1)
		    else if (ordre==2) then
                      aretesvg(i,j)=1._wp/30._wp*f(im3,jb)-13._wp/60._wp*f(im2,jb) &
		      +47._wp/60._wp*f(im1,jb)+9._wp/20._wp*f(ib,jb) &
		      -1._wp/20._wp*f(ip1,jb)
		      aretesvd(i,j)=-1._wp/20._wp*f(im2,jb)+9._wp/20._wp*f(im1,jb) &
		      +47._wp/60._wp*f(ib,jb)-13._wp/60._wp*f(ip1,jb) &
		      +1._wp/30._wp*f(ip2,jb)
		      areteshb(i,j)=1._wp/30._wp*f(ib,jm3)-13._wp/60._wp*f(ib,jm2) &
		      +47._wp/60._wp*f(ib,jm1)+9._wp/20._wp*f(ib,jb) &
		      -1._wp/20._wp*f(ib,jp1)
		      areteshh(i,j)=-1._wp/20._wp*f(ib,jm2)+9._wp/20._wp*f(ib,jm1) &
		      +47._wp/60._wp*f(ib,jb)-13._wp/60._wp*f(ib,jp1) &
		      +1._wp/30._wp*f(ib,jp2)
		    end if
		  end do
		end do	

		do i=0,N0-1
		  do j=0,N1-1
		    im3=modulo(i-3,N0)
		    im2=modulo(i-2,N0)
		    im1=modulo(i-1,N0)
		    ib=modulo(i,N0)
		    ip1=modulo(i+1,N0)
		    ip2=modulo(i+2,N0)
		    jm3=modulo(j-3,N1)
		    jm2=modulo(j-2,N1)
		    jm1=modulo(j-1,N1)
		    jb=modulo(j,N1)
		    jp1=modulo(j+1,N1)
		    jp2=modulo(j+2,N1)

		    if (ordre==1) then
		      sommetsbg(i,j)=-1._wp/6._wp*aretesvg(ib,jm2)+5._wp/6._wp*aretesvg(ib,jm1) &
		      +1._wp/3._wp*aretesvg(ib,jb)
		      sommetshg(i,j)=1._wp/3._wp*aretesvg(ib,jm1)+5._wp/6._wp*aretesvg(ib,jb) &
		      -1._wp/6._wp*aretesvg(ib,jp1)
		      sommetsbd(i,j)=-1._wp/6._wp*aretesvd(ib,jm2)+5._wp/6._wp*aretesvd(ib,jm1) &
		      +1._wp/3._wp*aretesvd(ib,jb)
		      sommetshd(i,j)=1._wp/3._wp*aretesvd(ib,jm1)+5._wp/6._wp*aretesvd(ib,jb) &
		      -1._wp/6._wp*aretesvd(ib,jp1)
		    else if (ordre==2) then
		      sommetsbg(i,j)=1._wp/30._wp*aretesvg(ib,jm3)-13._wp/60._wp*aretesvg(ib,jm2) &
		      +47._wp/60._wp*aretesvg(ib,jm1)+9._wp/20._wp*aretesvg(ib,jb) &
		      -1._wp/20._wp*aretesvg(ib,jp1)
		      sommetshg(i,j)=-1._wp/20._wp*aretesvg(ib,jm2)+9._wp/20._wp*aretesvg(ib,jm1) &
		      +47._wp/60._wp*aretesvg(ib,jb)-13._wp/60._wp*aretesvg(ib,jp1) &
		      +1._wp/30._wp*aretesvg(ib,jp2)
		      sommetsbd(i,j)=1._wp/30._wp*aretesvd(ib,jm3)-13._wp/60._wp*aretesvd(ib,jm2) &
		      +47._wp/60._wp*aretesvd(ib,jm1)+9._wp/20._wp*aretesvd(ib,jb) &
		      -1._wp/20._wp*aretesvd(ib,jp1)
		      sommetshd(i,j)=-1._wp/20._wp*aretesvd(ib,jm2)+9._wp/20._wp*aretesvd(ib,jm1) &
		      +47._wp/60._wp*aretesvd(ib,jb)-13._wp/60._wp*aretesvd(ib,jp1) &
		      +1._wp/30._wp*aretesvd(ib,jp2)
		     end if
	          end do
		end do
      endif
	end subroutine aux2
	
	
	subroutine calcule_coeff(N0,N1,a_moyenne,i,j,x1,y1,x2,y2,res,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd,cas)

	        integer,intent(in)::N0,N1,i,j,cas
		real(wp),intent(in)::x1,y1,x2,y2
		real(wp),dimension(0:N0-1,0:N1-1),intent(in)::a_moyenne,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd
		real(wp)::xx1,xx2,yy1,yy2,aux,dax,day,bx,by,c,w00,w10,w01,w20,w02,w11,w21,w12,w22,&
                c00,c10,c01,c20,c02,c11,c12,c21,c22,sij,sip1j,sijp1,sip1jp1,avij,avip1j,ahij,ahijp1,&
		fij,fim1jm1,fim1j,fim1jp1,fijm1,fijp1,fip1jm1,fip1j,fip1jp1
		integer::im2,im1,ib,ip1,ip2,jm2,jm1,jb,jp1,jp2
		real(wp),intent(out)::res
		
		if(cas==0)then
		  return
		endif

		im2=modulo(i-2,N0)
		im1=modulo(i-1,N0)
		ib=modulo(i,N0)
		ip1=modulo(i+1,N0)
		ip2=modulo(i+2,N0)
		jm2=modulo(j-2,N1)
		jm1=modulo(j-1,N1)
		jb=modulo(j,N1)
		jp1=modulo(j+1,N1)
		jp2=modulo(j+2,N1)

!1: Lauritzen, 2: Lag3, 3: PPM

		xx1=x1-real(i,wp)
		xx2=x2-real(i,wp)
		yy1=y1-real(j,wp)
		yy2=y2-real(j,wp)

		w00=1._wp/2._wp*(xx2+xx1)*(yy2-yy1)
		w10=1._wp/6._wp*(xx2**2+xx2*xx1+xx1**2)*(yy2-yy1)
		w01=-1._wp/6._wp*(yy2**2+yy2*yy1+yy1**2)*(xx2-xx1)
		w20=1._wp/12._wp*(xx2+xx1)*(xx2**2+xx1**2)*(yy2-yy1)
		w02=-1._wp/12._wp*(yy2+yy1)*(yy2**2+yy1**2)*(xx2-xx1)
		w11=1._wp/24._wp*(yy2*(3._wp*xx2**2+2._wp*xx2*xx1+xx1**2)+yy1*(xx2**2+2._wp*xx2*xx1+3._wp*xx1**2))*(yy2-yy1)
		w21=-1._wp/60._wp*(6._wp*xx1**2*yy1**2+6._wp*xx2**2*yy2**2+xx1**2*yy2**2+xx2**2*yy1**2+3._wp*xx1**2*yy1*yy2+ &
		3._wp*xx2**2*yy1*yy2+3._wp*xx1*xx2*yy1**2+3._wp*xx1*xx2*yy2**2+4._wp*xx1*xx2*yy1*yy2)*(xx2-xx1)
		w12=1._wp/60._wp*(6._wp*xx1**2*yy1**2+6._wp*xx2**2*yy2**2+xx1**2*yy2**2+xx2**2*yy1**2+3._wp*xx1**2*yy1*yy2+ &
		3._wp*xx2**2*yy1*yy2+3._wp*xx1*xx2*yy1**2+3._wp*xx1*xx2*yy2**2+4._wp*xx1*xx2*yy1*yy2)*(yy2-yy1)
		w22=-1._wp/180._wp*(10._wp*xx1**2*yy1**3+10._wp*xx2**2*yy2**3+xx2**2*yy1**3+xx1**2*yy2**3+4._wp*xx1*xx2*yy1**3 &
		+4._wp*xx1*xx2*yy2**3+6._wp*xx1**2*yy1**2*yy2+3._wp*xx1**2*yy1*yy2**2+3._wp*xx2**2*yy1**2*yy2+6._wp*xx2**2*yy1*yy2**2 &
		+6._wp*xx1*xx2*yy1**2*yy2+6._wp*xx1*xx2*yy1*yy2**2)*(xx2-xx1)

	if (cas==1) then !Lauritzen

		dax=1._wp/12._wp*(-a_moyenne(ip2,jb)+8._wp*a_moyenne(ip1,jb)-8._wp*a_moyenne(im1,jb)+a_moyenne(im2,jb))
		day=1._wp/12._wp*(-a_moyenne(ib,jp2)+8._wp*a_moyenne(ib,jp1)-8._wp*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		bx=1._wp/4._wp*(a_moyenne(ip2,jb)-6._wp*a_moyenne(ip1,jb)+10._wp*a_moyenne(ib,jb)-6._wp*a_moyenne(im1,jb)+a_moyenne(im2,jb))
		by=1._wp/4._wp*(a_moyenne(ib,jp2)-6._wp*a_moyenne(ib,jp1)+10._wp*a_moyenne(ib,jb)-6._wp*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		c=1._wp/4._wp*(a_moyenne(ip1,jp1)-a_moyenne(im1,jp1)-a_moyenne(ip1,jm1)+a_moyenne(im1,jm1))

		c00=a_moyenne(ib,jb)-0.5_wp*dax-0.5_wp*day-1._wp/6._wp*bx-1._wp/6._wp*by+1._wp/4._wp*c
		c10=dax+bx-0.5_wp*c
		c01=day+by-0.5_wp*c
		c20=-bx
		c02=-by
		c11=c

		res=w00*c00+w10*c10+w01*c01+w20*c20+w02*c02+w11*c11

	else if (cas==2) then !Lagrange 3

		fim1jm1=a_moyenne(im1,jm1)
		fim1j=a_moyenne(im1,jb)
		fim1jp1=a_moyenne(im1,jp1)
		fijm1=a_moyenne(ib,jm1)
		fij=a_moyenne(ib,jb)
		fijp1=a_moyenne(ib,jp1)
		fip1jm1=a_moyenne(ip1,jm1)
		fip1j=a_moyenne(ip1,jb)
		fip1jp1=a_moyenne(ip1,jp1)

		c00=1._wp/9._wp*fim1jm1+5._wp/18._wp*fim1j-1._wp/18._wp*fim1jp1+5._wp/18._wp*fijm1+25._wp/36._wp*fij-5._wp/36._wp*fijp1 &
		-1._wp/18._wp*fip1jm1-5._wp/36._wp*fip1j+1._wp/36._wp*fip1jp1
		c10=-1._wp/3._wp*fim1jm1-5._wp/6._wp*fim1j+1._wp/6._wp*fim1jp1+1._wp/3._wp*fijm1+5._wp/6._wp*fij-1._wp/6._wp*fijp1
		c20=1._wp/6._wp*fim1jm1+5._wp/12._wp*fim1j-1._wp/12._wp*fim1jp1-1._wp/3._wp*fijm1-5._wp/6._wp*fij+1._wp/6._wp*fijp1 &
		+1._wp/6._wp*fip1jm1+5._wp/12._wp*fip1j-1._wp/12._wp*fip1jp1
		c01=-1._wp/3._wp*fim1jm1+1._wp/3._wp*fim1j-5._wp/6._wp*fijm1+5._wp/6._wp*fij+1._wp/6._wp*fip1jm1-1._wp/6._wp*fip1j
		c11=fim1jm1-fim1j-fijm1+fij
		c21=-1._wp/2._wp*fim1jm1+1._wp/2._wp*fim1j+fijm1-fij-1._wp/2._wp*fip1jm1+1._wp/2._wp*fip1j
		c02=1._wp/6._wp*fim1jm1-1._wp/3._wp*fim1j+1._wp/6._wp*fim1jp1+5._wp/12._wp*fijm1-5._wp/6._wp*fij+5._wp/12._wp*fijp1 &
		-1._wp/12._wp*fip1jm1+1._wp/6._wp*fip1j-1._wp/12._wp*fip1jp1
		c12=-1._wp/2._wp*fim1jm1+fim1j-1._wp/2._wp*fim1jp1+1._wp/2._wp*fijm1-fij+1._wp/2._wp*fijp1
		c22=1._wp/4._wp*fim1jm1-1._wp/2._wp*fim1j+1._wp/4._wp*fim1jp1-1._wp/2._wp*fijm1+fij-1._wp/2._wp*fijp1+1._wp/4._wp*fip1jm1 &
		-1._wp/2._wp*fip1j+1._wp/4._wp*fip1jp1

		res=w00*c00+w10*c10+w01*c01+w20*c20+w02*c02+w11*c11+w21*c21+w12*c12+w22*c22

	else if (cas==3) then !PPM

		sij=sommets(ib,jb)
		sip1j=sommets(ip1,jb)
		sijp1=sommets(ib,jp1)
		sip1jp1=sommets(ip1,jp1)
		avij=aretesv(ib,jb)
		avip1j=aretesv(ip1,jb)
		ahij=aretesh(ib,jb)
		ahijp1=aretesh(ib,jp1)
		fij=a_moyenne(ib,jb)

		!sij=49._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jb)+a_moyenne(ib,jb))&
		!-7._wp/144._wp*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jb)+a_moyenne(ip1,jb))&
		!-7._wp/144._wp*(a_moyenne(im1,jm2)+a_moyenne(ib,jm2)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!+1._wp/144._wp*(a_moyenne(im2,jm2)+a_moyenne(ip1,jm2)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))
		!sip1j=49._wp/144._wp*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jb)+a_moyenne(ip1,jb))&
		!-7._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jb)+a_moyenne(ip2,jb))&
		!-7._wp/144._wp*(a_moyenne(ib,jm2)+a_moyenne(ip1,jm2)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!+1._wp/144._wp*(a_moyenne(im1,jm2)+a_moyenne(ip2,jm2)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))
		!sijp1=49._wp/144._wp*(a_moyenne(im1,jb)+a_moyenne(ib,jb)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!-7._wp/144._wp*(a_moyenne(im2,jb)+a_moyenne(ip1,jb)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))&
		!-7._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jp2)+a_moyenne(ib,jp2))&
		!+1._wp/144._wp*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jp2)+a_moyenne(ip1,jp2))
		!sip1jp1=49._wp/144._wp*(a_moyenne(ib,jb)+a_moyenne(ip1,jb)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!-7._wp/144._wp*(a_moyenne(im1,jb)+a_moyenne(ip2,jb)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))&
		!-7._wp/144._wp*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jp2)+a_moyenne(ip1,jp2))&
		!+1._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jp2)+a_moyenne(ip2,jp2))
		!avij=7._wp/12._wp*(a_moyenne(im1,jb)+a_moyenne(ib,jb))-1._wp/12._wp*(a_moyenne(im2,jb)+a_moyenne(ip1,jb))
		!avip1j=7._wp/12._wp*(a_moyenne(ib,jb)+a_moyenne(ip1,jb))-1._wp/12._wp*(a_moyenne(im1,jb)+a_moyenne(ip2,jb))
		!ahij=7._wp/12._wp*(a_moyenne(ib,jm1)+a_moyenne(ib,jb))-1._wp/12._wp*(a_moyenne(ib,jm2)+a_moyenne(ib,jp1))
		!ahijp1=7._wp/12._wp*(a_moyenne(ib,jb)+a_moyenne(ib,jp1))-1._wp/12._wp*(a_moyenne(ib,jm1)+a_moyenne(ib,jp2))


		c00=sij
		c10=-4._wp*sij-2._wp*sip1j+6._wp*ahij
		c20=3._wp*sij+3._wp*sip1j-6._wp*ahij
		c01=-4._wp*sij-2._wp*sijp1+6._wp*avij
		c11=16._wp*sij+8._wp*sip1j-24._wp*ahij+8._wp*sijp1+4._wp*sip1jp1-12._wp*ahijp1-24._wp*avij-12._wp*avip1j+36._wp*fij
		c21=-12._wp*sij-12._wp*sip1j+24._wp*ahij-6._wp*sijp1-6*sip1jp1+12._wp*ahijp1+18._wp*avij+18._wp*avip1j-36._wp*fij
		c02=3._wp*sij+3._wp*sijp1-6._wp*avij
		c12=-12._wp*sij-6._wp*sip1j+18._wp*ahij-12._wp*sijp1-6._wp*sip1jp1+18._wp*ahijp1+24._wp*avij+12._wp*avip1j-36._wp*fij
		c22=9._wp*sij+9._wp*sip1j-18._wp*ahij+9._wp*sijp1+9._wp*sip1jp1-18._wp*ahijp1-18._wp*avij-18*avip1j+36._wp*fij

		res=w00*c00+w10*c10+w01*c01+w20*c20+w02*c02+w11*c11+w21*c21+w12*c12+w22*c22

	else if (cas==4) then

		sij=sommetshd(ib,jb)
		sip1j=sommetshg(ip1,jb)
		sijp1=sommetsbd(ib,jp1)
		sip1jp1=sommetsbg(ip1,jp1)
		avij=aretesvd(ib,jb)
		avip1j=aretesvg(ip1,jb)
		ahij=areteshh(ib,jb)
		ahijp1=areteshb(ib,jp1)
		fij=a_moyenne(ib,jb)

		c00=sij
		c10=-4._wp*sij-2._wp*sip1j+6._wp*ahij
		c20=3._wp*sij+3._wp*sip1j-6._wp*ahij
		c01=-4._wp*sij-2._wp*sijp1+6._wp*avij
		c11=16._wp*sij+8._wp*sip1j-24._wp*ahij+8._wp*sijp1+4._wp*sip1jp1-12._wp*ahijp1-24._wp*avij-12._wp*avip1j+36._wp*fij
		c21=-12._wp*sij-12._wp*sip1j+24._wp*ahij-6._wp*sijp1-6*sip1jp1+12._wp*ahijp1+18._wp*avij+18._wp*avip1j-36._wp*fij
		c02=3._wp*sij+3._wp*sijp1-6._wp*avij
		c12=-12._wp*sij-6._wp*sip1j+18._wp*ahij-12._wp*sijp1-6._wp*sip1jp1+18._wp*ahijp1+24._wp*avij+12._wp*avip1j-36._wp*fij
		c22=9._wp*sij+9._wp*sip1j-18._wp*ahij+9._wp*sijp1+9._wp*sip1jp1-18._wp*ahijp1-18._wp*avij-18*avip1j+36._wp*fij

		res=w00*c00+w10*c10+w01*c01+w20*c20+w02*c02+w11*c11+w21*c21+w12*c12+w22*c22

	  endif

	end subroutine calcule_coeff



	subroutine calcule_coeffh(N0,N1,a_moyenne,i,j,x1,x2,res,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd,cas)

	        integer,intent(in)::N0,N1,i,j,cas
		real(wp),intent(in)::x1,x2
		real(wp),dimension(0:N0-1,0:N1-1),intent(in)::a_moyenne,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd
		real(wp)::xx1,xx2,yy,aux,day,by,w01,w02,w21,w22,w00,w10,w20,w11,w12,c01,c02,c21,c22,c00,c10,c20,c11,c12,&
		sij,sip1j,sijp1,sip1jp1,avij,avip1j,ahij,ahijp1,&
		fij,fim1jm1,fim1j,fim1jp1,fijm1,fijp1,fip1jm1,fip1j,fip1jp1,c,yy1,yy2
		integer::im2,im1,ib,ip1,ip2,jm2,jm1,jb,jp1,jp2
		real(wp),intent(out)::res
		
		if(cas==0)then
		  return
		endif

		im2=modulo(i-2,N0)
		im1=modulo(i-1,N0)
		ib=modulo(i,N0)
		ip1=modulo(i+1,N0)
		ip2=modulo(i+2,N0)
		jm2=modulo(j-2,N1)
		jm1=modulo(j-1,N1)
		jb=modulo(j,N1)
		jp1=modulo(j+1,N1)
		jp2=modulo(j+2,N1)

!1: Lauritzen, 2: Lag3, 3: PPM

		xx1=x1-real(i,wp)
		xx2=x2-real(i,wp)

		w01=-1._wp/2._wp*(xx2-xx1)
		w02=-1._wp/3._wp*(xx2-xx1)
		w21=-1._wp/6._wp*(xx1**2+xx2**2+xx1*xx2)*(xx2-xx1)
		w22=-1._wp/9._wp*(xx1**2+xx2**2+xx1*xx2)*(xx2-xx1)

	if (cas==1) then !Lauritzen

		day=1._wp/12._wp*(-a_moyenne(ib,jp2)+8._wp*a_moyenne(ib,jp1)-8._wp*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		by=1._wp/4._wp*(a_moyenne(ib,jp2)-6._wp*a_moyenne(ib,jp1)+10._wp*a_moyenne(ib,jb)-6._wp*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		c=1._wp/4._wp*(a_moyenne(ip1,jp1)-a_moyenne(im1,jp1)-a_moyenne(ip1,jm1)+a_moyenne(im1,jm1))

		c01=day+by-0.5_wp*c
		c02=-by

		res=w01*c01+w02*c02

	else if (cas==2) then !Lagrange 3

		fim1jm1=a_moyenne(im1,jm1)
		fim1j=a_moyenne(im1,jb)
		fim1jp1=a_moyenne(im1,jp1)
		fijm1=a_moyenne(ib,jm1)
		fij=a_moyenne(ib,jb)
		fijp1=a_moyenne(ib,jp1)
		fip1jm1=a_moyenne(ip1,jm1)
		fip1j=a_moyenne(ip1,jb)
		fip1jp1=a_moyenne(ip1,jp1)

		c01=-1._wp/3._wp*fim1jm1+1._wp/3._wp*fim1j-5._wp/6._wp*fijm1+5._wp/6._wp*fij+1._wp/6._wp*fip1jm1-1._wp/6._wp*fip1j
		c21=-1._wp/2._wp*fim1jm1+1._wp/2._wp*fim1j+fijm1-fij-1._wp/2._wp*fip1jm1+1._wp/2._wp*fip1j
		c02=1._wp/6._wp*fim1jm1-1._wp/3._wp*fim1j+1._wp/6._wp*fim1jp1+5._wp/12._wp*fijm1-5._wp/6._wp*fij+5._wp/12._wp*fijp1 &
		-1._wp/12._wp*fip1jm1+1._wp/6._wp*fip1j-1._wp/12._wp*fip1jp1
		c22=1._wp/4._wp*fim1jm1-1._wp/2._wp*fim1j+1._wp/4._wp*fim1jp1-1._wp/2._wp*fijm1+fij-1._wp/2._wp*fijp1+1._wp/4._wp*fip1jm1 &
		-1._wp/2._wp*fip1j+1._wp/4._wp*fip1jp1

		res=w01*c01+w02*c02+w21*c21+w22*c22

	else if (cas==3) then !PPM

		sij=sommets(ib,jb)
		sip1j=sommets(ip1,jb)
		sijp1=sommets(ib,jp1)
		sip1jp1=sommets(ip1,jp1)
		avij=aretesv(ib,jb)
		avip1j=aretesv(ip1,jb)
		ahij=aretesh(ib,jb)
		ahijp1=aretesh(ib,jp1)
		fij=a_moyenne(ib,jb)

		!sij=49._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jb)+a_moyenne(ib,jb))&
		!-7._wp/144._wp*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jb)+a_moyenne(ip1,jb))&
		!-7._wp/144._wp*(a_moyenne(im1,jm2)+a_moyenne(ib,jm2)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!+1._wp/144._wp*(a_moyenne(im2,jm2)+a_moyenne(ip1,jm2)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))
		!sip1j=49._wp/144._wp*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jb)+a_moyenne(ip1,jb))&
		!-7._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jb)+a_moyenne(ip2,jb))&
		!-7._wp/144._wp*(a_moyenne(ib,jm2)+a_moyenne(ip1,jm2)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!+1._wp/144._wp*(a_moyenne(im1,jm2)+a_moyenne(ip2,jm2)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))
		!sijp1=49._wp/144._wp*(a_moyenne(im1,jb)+a_moyenne(ib,jb)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!-7._wp/144._wp*(a_moyenne(im2,jb)+a_moyenne(ip1,jb)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))&
		!-7._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jp2)+a_moyenne(ib,jp2))&
		!+1._wp/144._wp*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jp2)+a_moyenne(ip1,jp2))
		!sip1jp1=49._wp/144._wp*(a_moyenne(ib,jb)+a_moyenne(ip1,jb)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!-7._wp/144._wp*(a_moyenne(im1,jb)+a_moyenne(ip2,jb)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))&
		!-7._wp/144._wp*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jp2)+a_moyenne(ip1,jp2))&
		!+1._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jp2)+a_moyenne(ip2,jp2))
		!avij=7._wp/12._wp*(a_moyenne(im1,jb)+a_moyenne(ib,jb))-1._wp/12._wp*(a_moyenne(im2,jb)+a_moyenne(ip1,jb))
		!avip1j=7._wp/12._wp*(a_moyenne(ib,jb)+a_moyenne(ip1,jb))-1._wp/12._wp*(a_moyenne(im1,jb)+a_moyenne(ip2,jb))
		!ahij=7._wp/12._wp*(a_moyenne(ib,jm1)+a_moyenne(ib,jb))-1._wp/12._wp*(a_moyenne(ib,jm2)+a_moyenne(ib,jp1))
		!ahijp1=7._wp/12._wp*(a_moyenne(ib,jb)+a_moyenne(ib,jp1))-1._wp/12._wp*(a_moyenne(ib,jm1)+a_moyenne(ib,jp2))

		c01=-4._wp*sij-2._wp*sijp1+6._wp*avij
		c21=-12._wp*sij-12._wp*sip1j+24._wp*ahij-6._wp*sijp1-6*sip1jp1+12._wp*ahijp1+18._wp*avij+18._wp*avip1j-36._wp*fij
		c02=3._wp*sij+3._wp*sijp1-6._wp*avij
		c22=9._wp*sij+9._wp*sip1j-18._wp*ahij+9._wp*sijp1+9._wp*sip1jp1-18._wp*ahijp1-18._wp*avij-18*avip1j+36._wp*fij

		res=w01*c01+w02*c02+w21*c21+w22*c22

	else if (cas==4) then

		sij=sommetshd(ib,jb)
		sip1j=sommetshg(ip1,jb)
		sijp1=sommetsbd(ib,jp1)
		sip1jp1=sommetsbg(ip1,jp1)
		avij=aretesvd(ib,jb)
		avip1j=aretesvg(ip1,jb)
		ahij=areteshh(ib,jb)
		ahijp1=areteshb(ib,jp1)
		fij=a_moyenne(ib,jb)

		c01=-4._wp*sij-2._wp*sijp1+6._wp*avij
		c21=-12._wp*sij-12._wp*sip1j+24._wp*ahij-6._wp*sijp1-6*sip1jp1+12._wp*ahijp1+18._wp*avij+18._wp*avip1j-36._wp*fij
		c02=3._wp*sij+3._wp*sijp1-6._wp*avij
		c22=9._wp*sij+9._wp*sip1j-18._wp*ahij+9._wp*sijp1+9._wp*sip1jp1-18._wp*ahijp1-18._wp*avij-18*avip1j+36._wp*fij

		res=w01*c01+w02*c02+w21*c21+w22*c22

	  endif

	end subroutine calcule_coeffh



	subroutine calcule_coeffv(N0,N1,a_moyenne,i,j,y1,y2,res,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,aretesvd,&
                sommetsbg,sommetsbd,sommetshg,sommetshd,cas)

	        integer,intent(in)::N0,N1,i,j,cas
		real(wp),intent(in)::y1,y2
		real(wp),dimension(0:N0-1,0:N1-1),intent(in)::a_moyenne,aretesh,aretesv,sommets,areteshb,areteshh,aretesvg,&
                aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd
		real(wp)::xx,yy1,yy2,aux,dax,bx,by,c,w00,w10,w20,w11,w12,w01,w02,w21,w22,&
		c00,c10,c20,c11,c12,c01,c02,c21,c22,sij,sip1j,sijp1,sip1jp1,avij,avip1j,ahij,ahijp1,&
		fij,fim1jm1,fim1j,fim1jp1,fijm1,fijp1,fip1jm1,fip1j,fip1jp1,day,xx1,xx2
		integer::im2,im1,ib,ip1,ip2,jm2,jm1,jb,jp1,jp2
		real(wp),intent(out)::res
		
		if(cas==0)then
		  return
		endif

		im2=modulo(i-2,N0)
		im1=modulo(i-1,N0)
		ib=modulo(i,N0)
		ip1=modulo(i+1,N0)
		ip2=modulo(i+2,N0)
		jm2=modulo(j-2,N1)
		jm1=modulo(j-1,N1)
		jb=modulo(j,N1)
		jp1=modulo(j+1,N1)
		jp2=modulo(j+2,N1)

!1: Lauritzen, 2: Lag3, 3: PPM
		yy1=y1-real(j,wp)
		yy2=y2-real(j,wp)

		w00=yy2-yy1
		w10=1._wp/2._wp*(yy2-yy1)
		w20=1._wp/3._wp*(yy2-yy1)
		w11=1._wp/4._wp*(yy1+yy2)*(yy2-yy1)
		w12=1._wp/6._wp*(yy1**2+yy2**2+yy1*yy2)*(yy2-yy1)

	if (cas==1) then !Lauritzen

		dax=1._wp/12._wp*(-a_moyenne(ip2,jb)+8._wp*a_moyenne(ip1,jb)-8._wp*a_moyenne(im1,jb)+a_moyenne(im2,jb))
		day=1._wp/12._wp*(-a_moyenne(ib,jp2)+8._wp*a_moyenne(ib,jp1)-8._wp*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		bx=1._wp/4._wp*(a_moyenne(ip2,jb)-6._wp*a_moyenne(ip1,jb)+10._wp*a_moyenne(ib,jb)-6._wp*a_moyenne(im1,jb)+a_moyenne(im2,jb))
		by=1._wp/4._wp*(a_moyenne(ib,jp2)-6._wp*a_moyenne(ib,jp1)+10._wp*a_moyenne(ib,jb)-6._wp*a_moyenne(ib,jm1)+a_moyenne(ib,jm2))
		c=1._wp/4._wp*(a_moyenne(ip1,jp1)-a_moyenne(im1,jp1)-a_moyenne(ip1,jm1)+a_moyenne(im1,jm1))

		c00=a_moyenne(ib,jb)-0.5_wp*dax-0.5_wp*day-1._wp/6._wp*bx-1._wp/6._wp*by+1._wp/4._wp*c
		c10=dax+bx-0.5_wp*c
		c20=-bx
		c11=c

		res=w00*c00+w10*c10+w20*c20+w11*c11

	else if (cas==2) then !Lagrange 3

		fim1jm1=a_moyenne(im1,jm1)
		fim1j=a_moyenne(im1,jb)
		fim1jp1=a_moyenne(im1,jp1)
		fijm1=a_moyenne(ib,jm1)
		fij=a_moyenne(ib,jb)
		fijp1=a_moyenne(ib,jp1)
		fip1jm1=a_moyenne(ip1,jm1)
		fip1j=a_moyenne(ip1,jb)
		fip1jp1=a_moyenne(ip1,jp1)

		c00=1._wp/9._wp*fim1jm1+5._wp/18._wp*fim1j-1._wp/18._wp*fim1jp1+5._wp/18._wp*fijm1+25._wp/36._wp*fij-5._wp/36._wp*fijp1 &
		-1._wp/18._wp*fip1jm1-5._wp/36._wp*fip1j+1._wp/36._wp*fip1jp1
		c10=-1._wp/3._wp*fim1jm1-5._wp/6._wp*fim1j+1._wp/6._wp*fim1jp1+1._wp/3._wp*fijm1+5._wp/6._wp*fij-1._wp/6._wp*fijp1
		c20=1._wp/6._wp*fim1jm1+5._wp/12._wp*fim1j-1._wp/12._wp*fim1jp1-1._wp/3._wp*fijm1-5._wp/6._wp*fij+1._wp/6._wp*fijp1 &
		+1._wp/6._wp*fip1jm1+5._wp/12._wp*fip1j-1._wp/12._wp*fip1jp1
		c11=fim1jm1-fim1j-fijm1+fij
		c12=-1._wp/2._wp*fim1jm1+fim1j-1._wp/2._wp*fim1jp1+1._wp/2._wp*fijm1-fij+1._wp/2._wp*fijp1

		res=w00*c00+w10*c10+w20*c20+w11*c11+w12*c12

	else if (cas==3) then !PPM

		sij=sommets(ib,jb)
		sip1j=sommets(ip1,jb)
		sijp1=sommets(ib,jp1)
		sip1jp1=sommets(ip1,jp1)
		avij=aretesv(ib,jb)
		avip1j=aretesv(ip1,jb)
		ahij=aretesh(ib,jb)
		ahijp1=aretesh(ib,jp1)
		fij=a_moyenne(ib,jb)

		!sij=49._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jb)+a_moyenne(ib,jb))&
		!-7._wp/144._wp*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jb)+a_moyenne(ip1,jb))&
		!-7._wp/144._wp*(a_moyenne(im1,jm2)+a_moyenne(ib,jm2)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!+1._wp/144._wp*(a_moyenne(im2,jm2)+a_moyenne(ip1,jm2)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))
		!sip1j=49._wp/144._wp*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jb)+a_moyenne(ip1,jb))&
		!-7._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jb)+a_moyenne(ip2,jb))&
		!-7._wp/144._wp*(a_moyenne(ib,jm2)+a_moyenne(ip1,jm2)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!+1._wp/144._wp*(a_moyenne(im1,jm2)+a_moyenne(ip2,jm2)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))
		!sijp1=49._wp/144._wp*(a_moyenne(im1,jb)+a_moyenne(ib,jb)+a_moyenne(im1,jp1)+a_moyenne(ib,jp1))&
		!-7._wp/144._wp*(a_moyenne(im2,jb)+a_moyenne(ip1,jb)+a_moyenne(im2,jp1)+a_moyenne(ip1,jp1))&
		!-7._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ib,jm1)+a_moyenne(im1,jp2)+a_moyenne(ib,jp2))&
		!+1._wp/144._wp*(a_moyenne(im2,jm1)+a_moyenne(ip1,jm1)+a_moyenne(im2,jp2)+a_moyenne(ip1,jp2))
		!sip1jp1=49._wp/144._wp*(a_moyenne(ib,jb)+a_moyenne(ip1,jb)+a_moyenne(ib,jp1)+a_moyenne(ip1,jp1))&
		!-7._wp/144._wp*(a_moyenne(im1,jb)+a_moyenne(ip2,jb)+a_moyenne(im1,jp1)+a_moyenne(ip2,jp1))&
		!-7._wp/144._wp*(a_moyenne(ib,jm1)+a_moyenne(ip1,jm1)+a_moyenne(ib,jp2)+a_moyenne(ip1,jp2))&
		!+1._wp/144._wp*(a_moyenne(im1,jm1)+a_moyenne(ip2,jm1)+a_moyenne(im1,jp2)+a_moyenne(ip2,jp2))
		!avij=7._wp/12._wp*(a_moyenne(im1,jb)+a_moyenne(ib,jb))-1._wp/12._wp*(a_moyenne(im2,jb)+a_moyenne(ip1,jb))
		!avip1j=7._wp/12._wp*(a_moyenne(ib,jb)+a_moyenne(ip1,jb))-1._wp/12._wp*(a_moyenne(im1,jb)+a_moyenne(ip2,jb))
		!ahij=7._wp/12._wp*(a_moyenne(ib,jm1)+a_moyenne(ib,jb))-1._wp/12._wp*(a_moyenne(ib,jm2)+a_moyenne(ib,jp1))
		!ahijp1=7._wp/12._wp*(a_moyenne(ib,jb)+a_moyenne(ib,jp1))-1._wp/12._wp*(a_moyenne(ib,jm1)+a_moyenne(ib,jp2))

		c00=sij
		c10=-4._wp*sij-2._wp*sip1j+6._wp*ahij
		c20=3._wp*sij+3._wp*sip1j-6._wp*ahij
		c11=16._wp*sij+8._wp*sip1j-24._wp*ahij+8._wp*sijp1+4._wp*sip1jp1-12._wp*ahijp1-24._wp*avij-12._wp*avip1j+36._wp*fij
		c12=-12._wp*sij-6._wp*sip1j+18._wp*ahij-12._wp*sijp1-6._wp*sip1jp1+18._wp*ahijp1+24._wp*avij+12._wp*avip1j-36._wp*fij

		res=w00*c00+w10*c10+w20*c20+w11*c11+w12*c12

	else if (cas==4) then

		sij=sommetshd(ib,jb)
		sip1j=sommetshg(ip1,jb)
		sijp1=sommetsbd(ib,jp1)
		sip1jp1=sommetsbg(ip1,jp1)
		avij=aretesvd(ib,jb)
		avip1j=aretesvg(ip1,jb)
		ahij=areteshh(ib,jb)
		ahijp1=areteshb(ib,jp1)
		fij=a_moyenne(ib,jb)
		c00=sij

		c10=-4._wp*sij-2._wp*sip1j+6._wp*ahij
		c20=3._wp*sij+3._wp*sip1j-6._wp*ahij
		c11=16._wp*sij+8._wp*sip1j-24._wp*ahij+8._wp*sijp1+4._wp*sip1jp1-12._wp*ahijp1-24._wp*avij-12._wp*avip1j+36._wp*fij
		c12=-12._wp*sij-6._wp*sip1j+18._wp*ahij-12._wp*sijp1-6._wp*sip1jp1+18._wp*ahijp1+24._wp*avij+12._wp*avip1j-36._wp*fij

		res=w00*c00+w10*c10+w20*c20+w11*c11+w12*c12

	  endif

	end subroutine calcule_coeffv

  subroutine compute_carac_dir_dir(dom,E0,E1,N0,N1,dt,timecase,eps,carac,marge)
    integer,intent(in)::N0,N1,timecase(0:1),marge
    real(wp),dimension(0:1,0:1),intent(in)::dom
    real(wp),intent(in)::dt,eps
    real(wp),dimension(-marge:N0+marge,-marge:N1+marge)::E0,E1
    real(wp),dimension(2,-1:N0,-1:N1)::carac  
    integer::i,j,ii,nstep=10,ix,iy,jj
    real(wp)::fval,xx,yy,xx0,yy0,xxn,yyn,x,y,ax,ay
   
    select case(timecase(0))
      case(1) !Euler
        do j=0,N1-1
          do i=0,N0-1
            xx=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)-E0(i,j)*dt
            yy=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)-E1(i,j)*dt
	    !print *,i,j,E0(i,j)!xx/dom(1,0)        
            carac(1,i,j)=xx
            carac(2,i,j)=yy
          enddo
        enddo
      case(2) !symplectic Euler
        do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)
            yy0=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)        
            x=0._wp
	    ix=i
	    y=0._wp
	    iy=j
	    xx=xx0-E0(i,j)*dt
	    do ii=1,timecase(1)
	      x=(xx-dom(0,0))/dom(1,0)
              !x=x-real(floor(x),wp)
              x=x*real(N0,wp)
              ix=floor(x)
              x=x-real(ix,wp)
              !if(ix==N0)then
              !  ix=0
              !  x=0._wp
              !endif
	      xxn=xx
              xx=xx0-dt*((1._wp-x)*E0(ix,j)+x*E0(ix+1,j))
            enddo  
	    !if(abs(xxn-xx)>eps)then
	    !  print *,ii,xx,xxn,xx-xxn,ix,ix+1,x,E0(ix,j),E0(ix+1,j)
	    !  stop
	    !endif
            yy=yy0-dt*((1._wp-x)*E1(ix,j)+x*E1(ix+1,j))  
            carac(1,i,j)=xx
            carac(2,i,j)=yy            
          enddo
        enddo	
       case(3) !symplectic Verlet
	do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)
            yy0=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)
            x=0._wp
	    ix=i
	    y=0._wp
	    iy=j
	    xx=xx0-E0(i,j)*0.5_wp*dt
	    yy=yy0-E1(i,j)*dt
	    do ii=1,timecase(1)
	      x=(xx-dom(0,0))/dom(1,0)
              !x=x-real(floor(x),wp)
              x=x*real(N0,wp)
              ix=floor(x)
              x=x-real(ix,wp)
              !if(ix==N0)then
              !  ix=0
              !  x=0._wp
              !endif
	      xxn=xx
              xx=xx0-0.5_wp*dt*((1._wp-x)*E0(ix,j)+x*E0(ix+1,j))
            enddo  
	    !if(abs(xxn-xx)>eps)then
	    !  print *,'no convergence x',ii,xx,xxn,xx-xxn,ix,ix+1,x,E0(ix,j),E0(ix+1,j)
	    !  stop
	    !endif
	    do ii=1,timecase(1)
	      y=(yy-dom(0,1))/dom(1,1)
              !y=y-real(floor(y),wp)
              y=y*real(N1,wp)
              iy=floor(y)
              y=y-real(iy,wp)
              !if(iy==N1)then
              !  iy=0
              !  y=0._wp
              !endif
	      yyn=yy
              yy=yy0-0.5_wp*dt*((1._wp-y)*((1._wp-x)*E1(ix,iy)+x*E1(ix+1,iy))&
	      +y*((1._wp-x)*E1(ix,iy+1)+x*E1(ix+1,iy+1)))
	      yy=yy-0.5_wp*dt*((1._wp-x)*E1(ix,j)+x*E1(ix+1,j))
            enddo
	    !if(abs(yyn-yy)>eps)then
	    !  print *,'no convergence y',ii,yy,yyn,yy-yyn,iy,iy+1,y
	    !  stop
	    !endif
            xx=xx-0.5_wp*dt*((1._wp-y)*((1._wp-x)*E0(ix,iy)+x*E0(ix+1,iy))&
	      +y*((1._wp-x)*E0(ix,iy+1)+x*E0(ix+1,iy+1)))  
            carac(1,i,j)=xx
            carac(2,i,j)=yy
          enddo
        enddo
      case(4)
        do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)
            yy0=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)        
            x=0._wp
	    ix=i
	    y=0._wp
	    iy=j
	    ax=0._wp
	    ay=0._wp
	    do ii=1,timecase(1)
	      xx=xx0-ax
	      yy=yy0-ay
	      x=(xx-dom(0,0))/dom(1,0)
              !x=x-real(floor(x),wp)
              x=x*real(N0,wp)
              ix=floor(x)
              x=x-real(ix,wp)
              !if(ix==N0)then
              !  ix=0
              !  x=0._wp
              !endif
	      y=(yy-dom(0,1))/dom(1,1)
              !y=y-real(floor(y),wp)
              y=y*real(N1,wp)
              iy=floor(y)
              y=y-real(iy,wp)
              !if(iy==N1)then
              !  iy=0
              !  y=0._wp
              !endif
	      ax=0.5_wp*dt*((1._wp-y)*((1._wp-x)*E0(ix,iy)+x*E0(ix+1,iy))&
	      +y*((1._wp-x)*E0(ix,iy+1)+x*E0(ix+1,iy+1)))
	      ay=0.5_wp*dt*((1._wp-y)*((1._wp-x)*E1(ix,iy)+x*E1(ix+1,iy))&
	      +y*((1._wp-x)*E1(ix,iy+1)+x*E1(ix+1,iy+1)))
	      xxn=xx
	      yyn=yy
	      xx=xx0-ax
	      yy=yy0-ay
            enddo  
	    !if(abs(xxn-xx)+abs(yyn-yy)>eps)then
	    !  print *,ii,xx,xxn,xx-xxn,yy,yyn,yy-yyn
	    !  stop
	    !endif
	    xx=xx0-2._wp*ax
	    yy=yy0-2._wp*ay	    
            carac(1,i,j)=xx
            carac(2,i,j)=yy
          enddo
        enddo
        		
    end select
    
    !use of Neumann conditions to complete carac array:    
    do i=0,N0-1
      carac(1,i,N1) = carac(1,i,N1-1)
      carac(2,i,N1) = carac(2,i,N1-1)+dom(1,1)/real(N1,wp)
      !carac(1,i,N1+1) = carac(1,i,N1)
      !carac(2,i,N1+1) = carac(2,i,N1)+dom(1,1)/real(N1,wp)
      carac(1,i,-1) = carac(1,i,0)
      carac(2,i,-1) = carac(2,i,0)-dom(1,1)/real(N1,wp)
    enddo

    do j=-1,N1
      carac(1,N0,j) = carac(1,N0-1,j)+dom(1,0)/real(N0,wp)
      carac(2,N0,j) = carac(2,N0-1,j)
      !carac(1,N0+1,j) = carac(1,N0,j)+dom(1,0)/real(N0,wp)
      !carac(2,N0+1,j) = carac(2,N0,j)      
      carac(1,-1,j) = carac(1,0,j)-dom(1,0)/real(N0,wp)
      carac(2,-1,j) = carac(2,0,j)    
    enddo 

  end subroutine compute_carac_dir_dir


  subroutine compute_carac_per_dir(dom,E0,E1,N0,N1,dt,timecase,eps,carac,marge)
    integer,intent(in)::N0,N1,timecase(0:1),marge
    real(wp),dimension(0:1,0:1),intent(in)::dom
    real(wp),intent(in)::dt,eps
    real(wp),dimension(-marge:N0+marge,-marge:N1+marge)::E0,E1
    real(wp),dimension(2,-1:N0,-1:N1)::carac  
    integer::i,j,ii,nstep=10,ix,iy,jj
    real(wp)::fval,xx,yy,xx0,yy0,xxn,yyn,x,y,ax,ay
   

!print*,'entrée'
    select case(timecase(0))
      case(1) !Euler
        do j=0,N1-1
          do i=0,N0-1
            xx=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)-E0(i,j)*dt
            yy=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)-E1(i,j)*dt       
            carac(1,i,j)=xx
            carac(2,i,j)=yy
          enddo
        enddo
      case(2) !symplectic Euler
        do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)
            yy0=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)       
            x=0._wp
	    ix=i
	    xx=xx0-E0(i,j)*dt
	    do ii=1,timecase(1)
	      x=(xx-dom(0,0))/dom(1,0)
              !x=x-real(floor(x),wp)
              x=x*real(N0,wp)
              ix=floor(x)
              x=x-real(ix,wp)
              !if(ix==N0)then
              !  ix=0
              !  x=0._wp
              !endif
	      xxn=xx
              xx=xx0-dt*((1._wp-x)*E0(ix,j)+x*E0(ix+1,j))
            enddo  
	    !if(abs(xxn-xx)>eps)then
	    !  print *,ii,xx,xxn,xx-xxn,ix,ix+1,x,E0(ix,j),E0(ix+1,j)
	    !  stop
	    !endif
            yy=yy0-dt*((1._wp-x)*E1(ix,j)+x*E1(ix+1,j))  
            carac(1,i,j)=xx
            carac(2,i,j)=yy            
          enddo
        enddo	
       case(3) !symplectic Verlet
	do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)
            yy0=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)
            x=0._wp
	    ix=i
	    y=0._wp
	    iy=j
	    xx=xx0-E0(i,j)*0.5*dt
	    yy=yy0-E1(i,j)*dt
	    do ii=1,timecase(1)
	      x=(xx-dom(0,0))/dom(1,0)
              !x=x-real(floor(x),wp)
              x=x*real(N0,wp)
              ix=floor(x)
              x=x-real(ix,wp)
              !if(ix==N0)then
              !  ix=0
              !  x=0._wp
              !endif
	      xxn=xx
              xx=xx0-0.5_wp*dt*((1._wp-x)*E0(ix,j)+x*E0(ix+1,j))
            enddo  
	    !if(abs(xxn-xx)>eps)then
	    !  print *,'no convergence x',ii,xx,xxn,xx-xxn,ix,ix+1,x,E0(ix,j),E0(ix+1,j)
	    !  stop
	    !endif
	    do ii=1,timecase(1)
	      y=(yy-dom(0,1))/dom(1,1)
              !y=y-real(floor(y),wp)
              y=y*real(N1,wp)
              iy=floor(y)
              y=y-real(iy,wp)
              !if(iy==N1)then
              !  iy=0
              !  y=0._wp
              !endif
	      yyn=yy
              yy=yy0-0.5_wp*dt*((1._wp-y)*((1._wp-x)*E1(ix,iy)+x*E1(ix+1,iy))&
	      +y*((1._wp-x)*E1(ix,iy+1)+x*E1(ix+1,iy+1)))
	      yy=yy-0.5_wp*dt*((1._wp-x)*E1(ix,j)+x*E1(ix+1,j))
            enddo
	    !if(abs(yyn-yy)>eps)then
	    !  print *,'no convergence y',ii,yy,yyn,yy-yyn,iy,iy+1,y
	    !  stop
	    !endif
            xx=xx-0.5_wp*dt*((1._wp-y)*((1._wp-x)*E0(ix,iy)+x*E0(ix+1,iy))&
	      +y*((1._wp-x)*E0(ix,iy+1)+x*E0(ix+1,iy+1)))  
            carac(1,i,j)=xx
            carac(2,i,j)=yy
          enddo
        enddo
      case(4)
        do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)
            yy0=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)       
            x=0._wp
	    ix=i
	    ax=0._wp
	    ay=0._wp
	    do ii=1,timecase(1)
	      xx=xx0-ax
	      yy=yy0-ay
	      x=(xx-dom(0,0))/dom(1,0)
              !x=x-real(floor(x),wp)
              x=x*real(N0,wp)
              ix=floor(x)
              x=x-real(ix,wp)
              !if(ix==N0)then
              !  ix=0
              !  x=0._wp
              !endif
	      y=(yy-dom(0,1))/dom(1,1)
              !y=y-real(floor(y),wp)
              y=y*real(N1,wp)
              iy=floor(y)
              y=y-real(iy,wp)
              !if(iy==N1)then
              !  iy=0
              !  y=0._wp
              !endif
	      ax=0.5_wp*dt*((1._wp-y)*((1._wp-x)*E0(ix,iy)+x*E0(ix+1,iy))&
	      +y*((1._wp-x)*E0(ix,iy+1)+x*E0(ix+1,iy+1)))
	      ay=0.5_wp*dt*((1._wp-y)*((1._wp-x)*E1(ix,iy)+x*E1(ix+1,iy))&
	      +y*((1._wp-x)*E1(ix,iy+1)+x*E1(ix+1,iy+1)))
	      xxn=xx
	      yyn=yy
	      xx=xx0-ax
	      yy=yy0-ay
            enddo  
	    !if(abs(xxn-xx)+abs(yyn-yy)>eps)then
	    !  print *,ii,xx,xxn,xx-xxn,yy,yyn,yy-yyn
	    !  stop
	    !endif
	    xx=xx0-2._wp*ax
	    yy=yy0-2._wp*ay	    
            carac(1,i,j)=xx
            carac(2,i,j)=yy
          enddo
        enddo
        		
    end select

    !use of periodic/Nauman conditions to complete carac array:   
    do i=0,N0-1
      carac(1,i,N1) = carac(1,i,N1-1)
      carac(2,i,N1) = carac(2,i,N1-1)+dom(1,1)/real(N1,wp)
      carac(1,i,-1) = carac(1,i,0)
      carac(2,i,-1) = carac(2,i,0)-dom(1,1)/real(N1,wp)
    enddo
   
    do j=-1,N1
      carac(1,N0,j) = carac(1,0,j)+dom(1,0)
      carac(2,N0,j) = carac(2,0,j)          
      carac(1,-1,j) = carac(1,N0-1,j)-dom(1,0)
      carac(2,-1,j) = carac(2,N0-1,j)      
    enddo
      	   
  end subroutine compute_carac_per_dir

  subroutine compute_carac_per_per(dom,E0,E1,N0,N1,timecase,carac)
    integer,intent(in)::N0,N1,timecase(0:1)
    real(wp),dimension(0:1,0:1),intent(in)::dom
    !real(wp),intent(in)::dt
    real(wp),dimension(1:N0,1:N1)::E0,E1
    real(wp),dimension(2,-1:N0,-1:N1)::carac  
    integer::i,j,ii,nstep=10,ix,iy,jj
    real(wp)::fval,xx,yy,xx0,yy0,xxn,yyn,x,y,ax,ay

   
    select case(timecase(0))
      case(1) !Euler
        do j=0,N1-1
          do i=0,N0-1
            xx=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)-E0(i+1,j+1)!*dt
            yy=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)-E1(i+1,j+1)!*dt       
            carac(1,i,j)=xx
            carac(2,i,j)=yy
          enddo
        enddo
      case(2) !symplectic Euler
        do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)
            yy0=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)       
            x=0._wp
	    ix=i
	    xx=xx0-E0(i+1,j+1)!*dt
	    do ii=1,timecase(1)
	      x=(xx-dom(0,0))/dom(1,0)
              x=x-real(floor(x),wp)
              x=x*real(N0,wp)
              ix=floor(x)
              x=x-real(ix,wp)
              if(ix==N0)then
                ix=0
                x=0._wp
              endif
	      xxn=xx
              xx=xx0-((1._wp-x)*E0(ix+1,j+1)+x*E0(modulo(ix+1,N0)+1,j+1))!*dt
            enddo  
	    !if(abs(xxn-xx)>eps)then
	    !  print *,ii,xx,xxn,xx-xxn,ix,ix+1,x,E0(ix+1,j+1),E0(ix+2,j+1)
	    !  stop
	    !endif
            yy=yy0-((1._wp-x)*E1(ix+1,j+1)+x*E1(modulo(ix+1,N0)+1,j+1))!*dt
            carac(1,i,j)=xx
            carac(2,i,j)=yy            
          enddo
        enddo	
      case(3) !symplectic Verlet
	do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,wp)*dom(1,0)/real(N0,wp)
            yy0=dom(0,1)+real(j,wp)*dom(1,1)/real(N1,wp)
            x=0._wp
	    ix=i
	    y=0._wp
	    iy=j
	    xx=xx0-E0(i+1,j+1)*0.5!*dt
	    yy=yy0-E1(i+1,j+1)!*dt
	    do ii=1,timecase(1)
	      x=(xx-dom(0,0))/dom(1,0)
              x=x-real(floor(x),wp)
              x=x*real(N0,wp)
              ix=floor(x)
              x=x-real(ix,wp)
              if(ix==N0)then
                ix=0
                x=0._wp
              endif
	      xxn=xx
              xx=xx0-0.5_wp*((1._wp-x)*E0(ix+1,j+1)+x*E0(modulo(ix+1,N0)+1,j+1))!*dt
            enddo  
	    !if(abs(xxn-xx)>eps)then
	    !  print *,'no convergence x',ii,xx,xxn,xx-xxn,ix,ix+1,x,E0(ix+1,j+1),E0(ix+2,j+1)
	    !  stop
	    !endif
	    do ii=1,timecase(1)
	      y=(yy-dom(0,1))/dom(1,1)
              y=y-real(floor(y),wp)
              y=y*real(N1,wp)
              iy=floor(y)
              y=y-real(iy,wp)
              if(iy==N1)then
                iy=0
                y=0._wp
              endif
	      yyn=yy
              yy=yy0-0.5_wp*((1._wp-y)*((1._wp-x)*E1(ix+1,iy+1)+x*E1(modulo(ix+1,N0)+1,iy+1))&
	      +y*((1._wp-x)*E1(ix+1,modulo(iy+1,N1)+1)+x*E1(modulo(ix+1,N0)+1,modulo(iy+1,N1)+1)))!*dt
	      yy=yy-0.5_wp*((1._wp-x)*E1(ix+1,j+1)+x*E1(modulo(ix+1,N0)+1,j+1))!*dt
            enddo
	    !if(abs(yyn-yy)>eps)then
	    !  print *,'no convergence y',ii,yy,yyn,yy-yyn,iy,iy+1,y
	    !  stop
	    !endif
            xx=xx-0.5_wp*((1._wp-y)*((1._wp-x)*E0(ix+1,iy+1)+x*E0(modulo(ix+1,N0)+1,iy+1))&
	      +y*((1._wp-x)*E0(ix+1,modulo(iy+1,N1)+1)+x*E0(modulo(ix+1,N0)+1,modulo(iy+1,N1)+1)))!*dt  
            carac(1,i,j)=xx
            carac(2,i,j)=yy
          enddo
        enddo
        		
    end select

    !use of periodic conditions to complete carac array:   
    do i=0,N0-1
      carac(1,i,N1) = carac(1,i,0)
      carac(2,i,N1) = carac(2,i,0)+dom(1,1)
      carac(1,i,-1) = carac(1,i,N1-1)
      carac(2,i,-1) = carac(2,i,N1-1)-dom(1,1)
    enddo
   
    do j=-1,N1
      carac(1,N0,j) = carac(1,0,j)+dom(1,0)
      carac(2,N0,j) = carac(2,0,j)          
      carac(1,-1,j) = carac(1,N0-1,j)-dom(1,0)
      carac(2,-1,j) = carac(2,N0-1,j)      
    enddo 


 
  end subroutine compute_carac_per_per



end module csl2dpp_class
