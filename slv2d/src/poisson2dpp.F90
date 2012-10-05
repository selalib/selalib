module poisson2dpp_module
  use used_precision
  use geometry_module
  use fft_module

  implicit none
  private
  public :: new, dealloc, solve, transposexy, transposeyx
  type, public:: poisson2dpp
     complex(wp), dimension(:,:), pointer :: rhot, ext, eyt
     type(fftclass) :: fftx, ffty
     type(geometry) :: geomx
     logical :: transpose
     integer :: jstartx,jendx,istartk,iendk
  end type poisson2dpp
  interface new
     module procedure new_poisson2dpp
  end interface
  interface dealloc
     module procedure dealloc_poisson2dpp
  end interface
  interface solve
     module procedure solve_poisson2dpp,solve_poisson1d
  end interface
contains
  subroutine new_poisson2dpp(this,rho,geomx,iflag, jstartx, jendx)
    type(poisson2dpp),intent(out) :: this
    real(wp), dimension(:,:), intent(in) :: rho
    type(geometry),intent(in)  :: geomx
    integer, intent(out) :: iflag
    integer, intent(in), optional ::  jstartx, jendx

    integer :: ierr ! indicateur d'erreur
    integer :: nxh1, i, j
#ifdef _OPENMP
    integer :: ipiece_size                 ! size of parallel region
    integer :: my_num, omp_get_thread_num    ! process number
    integer :: num_threads, omp_get_num_threads  ! number of threads
#endif
    ! indicateur d'erreur
    iflag = 0
    ! on commence par utiliser la fonction rho(x,y)
    this%transpose=.false.
    ! definition des bandes de calcul (en n'oubliant pas le cas sequentiel)
    ! le decoupage des tableaux exterieurs n'est pas gere par le module
    ! jstart et jend sont donnees en parametre d'entree dans le cas parallele
!!$    if (.not.(present(jstartx))) then
!!$       this%jstartx = 1
!!$    else
!!$       this%jstartx = jstartx
!!$    end if
!!$    if (.not.(present(jendx))) then
!!$       this%jendx = geomx%ny
!!$    else
!!$       this%jendx = jendx
!!$    end if
    ! on ne resout pas Poisson en parallele pour l'instant
    this%jstartx = 1
    this%jendx = geomx%ny
    ! la taille totale de la zone en kx est nxh1
    ! ATTENTION : les tableaux concernees par ce decoupage sont rhot,
    ! ext et eyt. Ce sont des complexes. Pour cette raison leur taille est
    ! la moitie de celle des tableaux reels correspondants.
    nxh1 = geomx%nx/2
    this%istartk = 1  ! cas sequentiel
    this%iendk = nxh1+1 ! cas sequentiel
#ifdef _OPENMP
    num_threads = omp_get_num_threads()
    ! ipiece_size = n/num_threads rounded up
    ipiece_size = (nxh1 + num_threads - 1) / num_threads
    ! zone a traiter en fonction du numero de process
    my_num = omp_get_thread_num()
    this%istartk = my_num * ipiece_size + 1
    this%iendk = min(this%istartk - 1 + ipiece_size, nxh1+1)
!print*,'zones k ',my_num,this%istartk,this%iendk, ipiece_size,nxh1
!print*,'zones x ',my_num,this%jstartx,this%jendx
#endif
!print*,'zones k ',this%istartk,this%iendk
!print*,'zones x ',this%jstartx,this%jendx

    ! initialisation de la geometrie
    this%geomx=geomx
    ! allocation memoire
    allocate(this%rhot(this%geomx%ny,this%istartk:this%iendk),stat=ierr)
    if (ierr.ne.0) iflag=20
    this%rhot=0.
    allocate(this%ext(this%geomx%ny,this%istartk:this%iendk),stat=ierr)
    if (ierr.ne.0) iflag=30
    this%ext=0.
    allocate(this%eyt(this%geomx%ny,this%istartk:this%iendk),stat=ierr)
    if (ierr.ne.0) iflag=40
    this%eyt=0.
    ! initialisation des fft (le nombre de points pris en compte est n)
    call initfft(this%fftx, rho, this%geomx%nx)
    call initfft(this%ffty, this%rhot, this%geomx%ny)
  end subroutine new_poisson2dpp

  subroutine dealloc_poisson2dpp(this)
    type(poisson2dpp),intent(out) :: this
    deallocate(this%rhot,this%ext,this%eyt)
  end subroutine dealloc_poisson2dpp

  subroutine solve_poisson2dpp(this,ex,ey,rho,nrj)
    type(poisson2dpp) :: this
    real(wp), dimension(:,:) :: ex,ey,rho
    real(wp), intent(out) :: nrj
    integer :: ik,jk,i,j ! indices de boucle
    integer :: istart
    real(wp) :: kx0,ky0, kx,ky, kx2, k2
    real(wp), dimension(size(rho,1),size(rho,2)) :: phi

    
    ! faire une FFT dans la direction x de rho
    call fft(this%fftx,rho(:,this%jstartx:this%jendx))

    ! transposition de rho stockee dans rhot
    call transposexy(this,rho)

    ! faire une FFT dans la direction y de rhot
    call fft(this%ffty,this%rhot(:,this%istartk:this%iendk))

    ! calcul de la transformee de Fourier de E a partir de celle de rho
    kx0=2*pi/(this%geomx%nx*this%geomx%dx)
    ky0=2*pi/(this%geomx%ny*this%geomx%dy)

    ! pour jk=1
    jk = 1    
    if (this%istartk.eq.1) then
       this%ext(1,1)=0.
       this%eyt(1,1)=0.
       istart = 2
    else
       istart = this%istartk
    end if
    do ik=istart,this%iendk
       kx= (ik-1)*kx0
       k2 = kx*kx
       this%ext(jk,ik)=-dcmplx(zero,kx/k2)*this%rhot(jk,ik)
       this%eyt(jk,ik)= 0.
       this%rhot(jk,ik)= this%rhot(jk,ik)/k2
    end do

    ! pour jk > 1
    do ik=this%istartk,this%iendk
       kx= (ik-1)*kx0
       kx2 = kx*kx
       do jk =  2, this%geomx%ny/2+1
          ky= (jk-1)*ky0
          k2= kx2 +ky*ky

          this%ext(jk,ik)=-dcmplx(zero,kx/k2)*this%rhot(jk,ik)
          this%eyt(jk,ik)=-dcmplx(zero,ky/k2)*this%rhot(jk,ik)
          this%rhot(jk,ik)= this%rhot(jk,ik)/k2
       end do
       do jk = this%geomx%ny/2+2 , this%geomx%ny     
          ky= (jk-1-this%geomx%ny)*ky0
          k2= kx2 +ky*ky

          this%ext(jk,ik)=-dcmplx(zero,kx/k2)*this%rhot(jk,ik)
          this%eyt(jk,ik)=-dcmplx(zero,ky/k2)*this%rhot(jk,ik)
          this%rhot(jk,ik)= this%rhot(jk,ik)/k2
       end do
    end do

    ! faire une FFT inverse  dans la direction y de E
!print*,'ext ', shape(this%ext(:,this%istartk:this%iendk))
    call fftinv(this%ffty,this%ext(:,this%istartk:this%iendk))
    call fftinv(this%ffty,this%eyt(:,this%istartk:this%iendk))
    call fftinv(this%ffty,this%rhot(:,this%istartk:this%iendk))

    ! transposition de E
    !call transposeyx(this,ex,ey)
    call transposeyx(this,ex,ey,phi) 
!!$do i= 1,this%geomx%nx
!!$   print*, 'rho fin solve ', i, (ex(i,j), j=this%jstartx,this%jendx)
!!$enddo
!print*, 'decoupage ', this%jstartx, this%jendx
    ! faire une FFT inverse  dans la direction x de E
    call fftinv(this%fftx,ex(:,this%jstartx:this%jendx))
    call fftinv(this%fftx,ey(:,this%jstartx:this%jendx))
    call fftinv(this%fftx,phi(:,this%jstartx:this%jendx))

    nrj=0._wp
    do i=1,this%geomx%nx
       do j=1,this%geomx%ny
          nrj=nrj+ex(i,j)*ex(i,j)+ey(i,j)*ey(i,j)          
       enddo
    enddo
    !print *,'size',this%geomx%nx, this%geomx%ny
    nrj=nrj*this%geomx%dx*this%geomx%dy
    if (nrj>1.e-30) then 
       nrj=0.5_wp*log(nrj)
    else
       nrj=-10**9
    endif
    print *,'test fft',nrj,0.5_wp*log(0.08_wp*pi*pi),abs(nrj-0.5_wp*log(0.08_wp*pi*pi))

    open(unit=10,file='exinitfft')
    do j=1,this%geomx%ny
    do i=1,this%geomx%nx
       write(10,*) (i-1)*this%geomx%dx,(j-1)*this%geomx%dy,ex(i,j), ey(i,j)
    enddo
    write(10,*) ' '
    enddo
    close(10)

    
  end subroutine solve_poisson2dpp

  subroutine solve_poisson1d(this,ex,rho,nrj)
    type(poisson2dpp) :: this
    real(wp), dimension(:,:) :: ex,rho
    real(wp), intent(out) :: nrj
    real(wp), dimension(1:this%geomx%nx,1:this%geomx%ny) :: tmp
    real(wp) :: avg,mass
    integer  :: i,j,im1

    ex(1,2) = 0._8
    avg=0._8
    do i=2,this%geomx%nx
       ex(i,2)=ex(i-1,2)+this%geomx%dx*(rho(i,2)-1._wp)                 !ex en x_{i+1/2} (cas 1D)
       avg=avg+ex(i,2)
    enddo
    avg=-avg/(1._8*this%geomx%nx)

    do i=1,this%geomx%nx
       ex(i,2)=avg+ex(i,2)
    enddo

    tmp(1:this%geomx%nx,1:this%geomx%ny)=ex(1:this%geomx%nx,1:this%geomx%ny)	
    do i=this%geomx%nx,2,-1
       im1=mod(i-1+this%geomx%nx,this%geomx%nx)
       ex(i,2)=0.5_wp*(tmp(im1,2)+tmp(i,2))
    enddo	
    ex(1,2)=0._wp 

    do i=1,this%geomx%nx
       do j=1,this%geomx%ny
          ex(i,j)=ex(i,2)
       enddo
    enddo
        
    nrj=0._wp
    do i=1,this%geomx%nx
       do j=1,this%geomx%ny
          nrj=nrj+ex(i,j)*ex(i,j)
       enddo
    enddo

    nrj=nrj*this%geomx%dx*this%geomx%dy!/this%geomx%y1
    if (nrj>1.e-30) then 
       nrj=0.5_wp*log(nrj)
    else
       nrj=-10**9
    endif
    print *,'nrj INIT',nrj,0.5_wp*log(0.08_wp*pi*pi),abs(nrj-0.5_wp*log(0.08_wp*pi*pi))

    open(unit=10,file='exinit')
    do j=1,this%geomx%ny
    do i=1,this%geomx%nx
       write(10,*) (i-1)*this%geomx%dx,(j-1)*this%geomx%dy,ex(i,j), rho(i,j)
    enddo
    write(10,*) ' '
    enddo
    close(10)
    
  end subroutine solve_poisson1d

  subroutine transposexy(this,rho)
    type(poisson2dpp) :: this
    real(wp), dimension(:,:) :: rho
    integer :: i,j ! indices de boucle
    integer :: nxh1
    ! on attend que tous les processeurs aient fait leur calcul
!$OMP BARRIER
    ! transposition
    nxh1 = this%geomx%nx/2
    do i=max(this%istartk,2),this%iendk-1
       do j=1,this%geomx%ny
          this%rhot(j,i)=dcmplx(rho(2*i-2,j),rho(2*i-1,j))
!          print*, i,j,this%rhot(j,i),rho(2*i-2,j),rho(2*i-1,j)
       end do
    end do
    if (this%istartk.eq.1) then
       do j=1,this%geomx%ny
          this%rhot(j,1) = dcmplx(rho(1,j),0.0)
       end do
    end if
    if (this%iendk.eq.nxh1+1) then
       do j=1,this%geomx%ny
          this%rhot(j,nxh1+1) = dcmplx(rho(2*nxh1,j),0.0)
       end do
    end if
    
   ! on attend que tous les processeurs aient fait leur transposition
!$OMP BARRIER
  end subroutine transposexy
  subroutine transposeyx(this,ex,ey,rho)
    type(poisson2dpp) :: this
    real(wp), dimension(:,:) :: ex,ey,rho
    integer :: i,j ! indices de boucle
    integer :: nxh1
    ! on attend que tous les processeurs aient fait leur calcul
!$OMP BARRIER
    ! transposition
    nxh1 = this%geomx%nx/2
    do j=1,this%geomx%ny
       if (this%istartk.eq.1) then
          ex(1,j)= real (this%ext(j,1))
          ey(1,j)= real (this%eyt(j,1))
          rho(1,j)= real (this%rhot(j,1))
       end if
       if (this%iendk.eq.nxh1+1) then
          ex(this%geomx%nx,j)= real (this%ext(j,nxh1+1))
          ey(this%geomx%nx,j)= real (this%eyt(j,nxh1+1))
          rho(this%geomx%nx,j)= real (this%rhot(j,nxh1+1))
       end if
       do i=max(this%istartk,2),this%iendk-1
          ex(2*i-2,j)=real(this%ext(j,i))
          ex(2*i-1,j)=imag(this%ext(j,i))
          ey(2*i-2,j)=real(this%eyt(j,i))
          ey(2*i-1,j)=imag(this%eyt(j,i))
          rho(2*i-2,j)=real(this%rhot(j,i))
          rho(2*i-1,j)=imag(this%rhot(j,i))
       end do
    end do

    ! on attend que tous les processeurs aient fait leur transposition
 !$OMP BARRIER
 end subroutine transposeyx

end module poisson2dpp_module
