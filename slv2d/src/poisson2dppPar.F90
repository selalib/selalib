module poisson2dpp_module
  use used_precision
  use geometry_module
  use fft_module
#ifdef _MPI
  use Module_MPI, my_num=>numero_processeur, num_threads=>nombre_processeurs
#endif
  implicit none
  private
  public :: new, dealloc, solve
#ifndef _MPI
  public :: transposexy, transposeyx
#endif
  type, public:: poisson2dpp
     real(wp), dimension(:,:), pointer :: rho, work1,work2
     complex(wp), dimension(:,:), pointer :: ext,eyt,rhot
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
     module procedure solve_poisson2dpp
  end interface
contains
  subroutine new_poisson2dpp(this,rho,geomx,iflag, jstartx, jendx)
    type(poisson2dpp),intent(out) :: this
#ifdef _MPI
    real(wp), dimension(:,jstartx:), intent(in) :: rho
#else
    real(wp), dimension(:,:), intent(in) :: rho
#endif
    type(geometry),intent(in)  :: geomx
    integer, intent(out) :: iflag
    integer, intent(in) ::  jstartx, jendx

    integer :: ierr ! indicateur d'erreur
    integer :: nxh1
#ifdef _OPENMP
    integer :: ipiece_size                 ! size of parallel region
    integer :: my_num, omp_get_thread_num    ! process number
    integer :: num_threads, omp_get_num_threads  ! number of threads
#endif
#ifdef _MPI
    integer :: ipiece_size             
#endif

    ! indicateur d'erreur
    iflag = 0
    ! on commence par utiliser la fonction rho(x,y)
    this%transpose=.false.
    ! definition des bandes de calcul (en n'oubliant pas le cas sequentiel)
    ! le decoupage des tableaux exterieurs n'est pas gere par le module
    ! jstart et jend sont donnees en parametre d'entree dans le cas parallele
    this%jstartx = jstartx
    this%jendx = jendx
    ! la taille totale de la zone en kx est nxh1
    ! ATTENTION : les tableaux concernees par ce decoupage sont rhot,
    ! ext et eyt. Ce sont des complexes. Pour cette raison leur taille est
    ! la moitie de celle des tableaux reels correspondants.
    nxh1 = geomx%nx/2
    this%istartk = 1  ! cas sequentiel
    this%iendk = nxh1+1 ! cas sequentiel
#ifdef _OPENMP
    num_threads = omp_get_num_threads()
    my_num = omp_get_thread_num()
#endif
#if defined _OPENMP || defined _MPI
    ! ipiece_size = n/num_threads rounded up
    ipiece_size = (nxh1 + num_threads - 1) / num_threads
    ! zone a traiter en fonction du numero de process
    if (my_num == 0 ) then
       this%istartk = 1
    else
       this%istartk = my_num * ipiece_size + 2
    end if
    this%iendk = min(this%istartk  + ipiece_size, nxh1+1)
print*,'zones ',my_num,this%istartk,this%iendk, ipiece_size,nxh1
#endif
    ! initialisation de la geometrie
    this%geomx=geomx
    ! allocation memoire
  
    allocate(this%rho(this%geomx%nx,this%jstartx:this%jendx),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, RHO, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,ierr)
#else
       stop
#endif 
    end if
    this%rho=0.
    allocate(this%rhot(this%geomx%ny,this%istartk:this%iendk),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, RHOT, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,ierr)
#else
       stop
#endif 
    end if  
    this%rhot=0.
!    allocate(this%work1(this%jstartx:this%jendx,this%geomx%nx),stat=ierr)
    allocate(this%work1(this%geomx%ny,this%geomx%nx),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, WORK, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,ierr)
#else
       stop
#endif 
    end if
    this%work1 = 0
    if (my_num == 0) then 
       allocate(this%work2(this%geomx%ny,1:2*this%iendk),stat=ierr)
    else
       allocate(this%work2(this%geomx%ny,2*this%istartk-1:2*this%iendk),stat=ierr)
    endif
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, WORK, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,ierr)
#else
       stop
#endif 
    end if
    this%work2 = 0
    allocate(this%ext(this%geomx%ny,this%istartk:this%iendk),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, ext, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,ierr)
#else
       stop
#endif
    end if
    this%ext=0.
    allocate(this%eyt(this%geomx%ny,this%istartk:this%iendk),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, eyt, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,ierr)
#else
       stop
#endif
    end if
    this%eyt=0.
    ! initialisation des fft (le nombre de points pris en compte est n)
print*, 'av initfft', rho(1,jstartx)
    call initfft(this%fftx, rho, this%geomx%nx)
    call initfft(this%ffty, this%rhot, this%geomx%ny)
print*, 'ap initfft', rho(1,jstartx)
  end subroutine new_poisson2dpp
!--------------------------------------------------
  subroutine dealloc_poisson2dpp(this)
    type(poisson2dpp),intent(out) :: this
    deallocate(this%rho,this%ext,this%eyt,this%rhot)
  end subroutine dealloc_poisson2dpp
!__________________________________________________
  subroutine solve_poisson2dpp(this,ex,ey,rho)
    type(poisson2dpp) :: this
#ifdef _MPI
    real(wp), dimension(:,this%jstartx:) :: ex,ey,rho
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: tag, mpierror
    real(wp),dimension(this%geomx%nx) :: phib
#else
    real(wp), dimension(:,:) :: ex,ey,rho
#endif
#ifdef _OPENMP
    integer :: my_num, num_threads, omp_get_num_threads, omp_get_thread_num

    num_threads = omp_get_num_threads()
    my_num = omp_get_thread_num()
#endif
    integer :: ik,jk,i,j ! indices de boucle
    integer :: istart
    real(wp) :: kx0,ky0, kx,ky, kx2, k2

    ! sauvegarde de rho pour les diagnostiques
    this%rho(:,this%jstartx:this%jendx) = rho(:,this%jstartx:this%jendx) 

    ! faire une FFT dans la direction x de rho
    call fft(this%fftx,rho(:,this%jstartx:this%jendx))

    ! transposition de rho stockee dans rhot
do i=1,this%geomx%nx
   write(*,"(i2,i2,' rho  ', 100(f6.1,1x))") my_num, i, (rho(i,j) , j=this%jstartx,this%jendx)
end do
    call transposexy(this,rho)
do j=1,this%geomx%ny
 !  write(*,"(i2,i2,' rhot  ', 100(f6.1,1x))") my_num, j, (this%rhot(j,i) , i=this%istartk,this%iendk)
   do i=this%istartk,this%iendk
      write(*,*) 'rhot ', my_num, j,i, this%rhot(j,i)
   end do
end do
    ! faire une FFT dans la direction y de rhot
    call fft(this%ffty,this%rhot(:,this%istartk:this%iendk))

do j=1,this%geomx%ny
   do i=this%istartk,this%iendk
      write(*,*) 'rhotfft ', my_num, j,i, this%rhot(j,i)
   end do
end do
    ! fa

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
print*,'apres ----------------'
do j=1,this%geomx%ny
   write(*,"(i2,i2,' rhot  ', 100(f6.1,1x))") my_num, j, (this%rhot(j,i) , i=this%istartk,this%iendk)
end do

    ! faire une FFT inverse  dans la direction y de Ex, Ey et rho 
    call fftinv(this%ffty,this%ext(:,this%istartk:this%iendk))
    call fftinv(this%ffty,this%eyt(:,this%istartk:this%iendk))
    call fftinv(this%ffty,this%rhot(:,this%istartk:this%iendk))

do j=1,this%geomx%ny
   write(*,"(i2,i2,' rhot  ', 100(f6.1,1x))") my_num, j, (this%rhot(j,i) , i=this%istartk,this%iendk)
end do

    ! transposition de Ex, Ey, rho
    call transposeyx(this,ex,ey,rho) 

do i=1,this%geomx%nx
   write(*,"(i2,i2,' rho  ', 100(f6.1,1x))") my_num, i, (rho(i,j) , j=this%jstartx,this%jendx)
end do

print*,' ----------------'
    ! faire une FFT inverse  dans la direction x de E
    call fftinv(this%fftx,ex(:,this%jstartx:this%jendx))
    call fftinv(this%fftx,ey(:,this%jstartx:this%jendx))
    call fftinv(this%fftx,rho(:,this%jstartx:this%jendx))

     ! on recopie les valeurs de rho dans rho pour les diagnostiques
!!$     do i=1,this%geomx%nx
!!$        do j=this%jstartx,this%jendx
!!$           rho(i,j) = this%rho(i,j)
!!$        end do
!!$     end do
#ifdef _MPI
    call mpi_barrier(MPI_COMM_WORLD,mpierror)
#elif defined _OPENMP
!$OMP BARRIER  
#endif
  end subroutine solve_poisson2dpp
  subroutine transposexy(this,rho)
    type(poisson2dpp) :: this
    real(wp), dimension(:,:) :: rho
    integer :: i,j ! indices de boucle
    integer :: nxh1,iflag
#ifdef _MPI
    call mpi_barrier(MPI_COMM_WORLD,iflag)
    ! transpose rho dans this%work
    call transpose(rho,this%work1, this%geomx%nx, this%geomx%ny, num_threads)
    ! convertit this%work en complexe
    nxh1 = this%geomx%nx/2
    do i=max(2,this%istartk),min(this%iendk,nxh1)
       do j=1,this%geomx%ny
          this%rhot(j,i)=dcmplx(this%work1(j,2*i-2),this%work1(j,2*i-1))
!print*, 'transp ', j,i,this%rhot(j,i),this%work1(j,2*i-2),this%work1(j,2*i-1)
       end do
    end do
    if (this%istartk == 1) then
       do j=1,this%geomx%ny
          this%rhot(j,1) = dcmplx(this%work1(j,1),0.0)
       end do
    end if
    if (this%iendk.eq.nxh1+1) then
       do j=1,this%geomx%ny
          this%rhot(j,nxh1+1) = dcmplx(this%work1(j,2*nxh1),0.0)
       end do
    end if
    call mpi_barrier(MPI_COMM_WORLD,iflag)
#else
    ! on attend que tous les processeurs aient fait leur calcul
!$OMP BARRIER
    ! transposition
    do i=this%istartk,this%iendk
       do j=1,this%geomx%ny
          this%rhot(j,i)=dcmplx(rho(2*i-1,j),rho(2*i,j))
       end do
    end do
    ! on attend que tous les processeurs aient fait leur transposition
!$OMP BARRIER
#endif
  end subroutine transposexy
  subroutine transposeyx(this,ex,ey,rho)
    type(poisson2dpp) :: this
    real(wp), dimension(:,:) :: ex,ey,rho
    integer :: i,j ! indices de boucle
    integer :: nxh1,iflag
#ifdef _MPI
    call mpi_barrier(MPI_COMM_WORLD,iflag)
    nxh1 = this%geomx%nx/2
    ! convertit ext en reel 
    do j=1,this%geomx%ny
       if (this%istartk.eq.1) then
          this%work2(j,1)= real (this%ext(j,1))
       end if
       if (this%iendk.eq.nxh1+1) then
          this%work2(j,this%geomx%nx)= real (this%ext(j,nxh1+1))
       end if
       do i=max(this%istartk,2),min(this%iendk,nxh1)
          this%work2(j,2*i-2)=real(this%ext(j,i))
          this%work2(j,2*i-1)=imag(this%ext(j,i))
       end do
    end do

    ! transpose this%work2 dans ex
    call transpose(this%work2,ex, this%geomx%ny, this%geomx%nx, num_threads)
    call mpi_barrier(MPI_COMM_WORLD,iflag)
    ! convertit eyt en reel
     do j=1,this%geomx%ny
       if (this%istartk.eq.1) then
          this%work2(j,1)= real (this%eyt(j,1))
       end if
       if (this%iendk.eq.nxh1+1) then
          this%work2(j,this%geomx%nx)= real (this%eyt(j,nxh1+1))
       end if
       do i=max(this%istartk,2),this%iendk-1
          this%work2(j,2*i-2)=real(this%eyt(j,i))
          this%work2(j,2*i-1)=imag(this%eyt(j,i))
       end do
    end do
    ! transpose this%work2 dans ey
    call transpose(this%work2,ey, this%geomx%ny, this%geomx%nx, num_threads)
    call mpi_barrier(MPI_COMM_WORLD,iflag)
    ! convertit rhot en reel
    do j=1,this%geomx%ny
       if (this%istartk.eq.1) then
          this%work2(j,1)= real (this%rhot(j,1))
       end if
       if (this%iendk.eq.nxh1+1) then
          this%work2(j,this%geomx%nx)= real (this%rhot(j,nxh1+1))
       end if
       do i=max(this%istartk,2),min(this%iendk,nxh1)
          this%work2(j,2*i-2)=real(this%rhot(j,i))
          this%work2(j,2*i-1)=imag(this%rhot(j,i))
       end do
    end do
do j= 1,this%geomx%ny
   do i= 1,this%geomx%nx
      print*, my_num, 'work2 ',j,i, this%work2(j,i)
   end do
end do
    ! transpose this%work2 dans rho
    call transpose(this%work2,rho, this%geomx%ny, this%geomx%nx, num_threads)
    call mpi_barrier(MPI_COMM_WORLD,iflag)
#else
    ! on attend que tous les processeurs aient fait leur calcul
!$OMP BARRIER
    ! transposition
    do i=this%istartk,this%iendk
       do j=1,this%geomx%ny
          ex(2*i-1,j)=real(this%ext(j,i))
          ex(2*i,j)=imag(this%ext(j,i))
          ey(2*i-1,j)=real(this%eyt(j,i))
          ey(2*i,j)=imag(this%eyt(j,i))
          rho(2*i-1,j)=real(this%rhot(j,i))
          rho(2*i,j)=imag(this%rhot(j,i))
       end do
    end do
    ! on attend que tous les processeurs aient fait leur transposition
 !$OMP BARRIER
#endif
 end subroutine transposeyx
end module poisson2dpp_module
