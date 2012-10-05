module poisson2ddir_module
  use used_precision
  use geometry_module
#ifdef _MPI
  use mpi
  use Module_MPI, my_num=>numero_processeur, num_threads=>nombre_processeurs
#endif
  implicit none
  private
  public :: new, dealloc, solve
#ifndef _MPI
  public :: transposexy, transposeyx
#else
  integer :: code
#endif
  type, public:: poisson2ddir
     real(wp), dimension(:,:), pointer :: rho, rhot1, rhot2,  work
     real(wp), dimension(:), pointer :: wsavex, wsavey
     type(geometry) :: geomx
     integer :: jstartx,jendx,istartk,iendk
  end type poisson2ddir
  interface new
     module procedure new_poisson2ddir
  end interface
  interface dealloc
     module procedure dealloc_poisson2ddir
  end interface
  interface solve
     module procedure solve_poisson2ddir
  end interface
contains
  subroutine new_poisson2ddir(this,rho,geomx,iflag, jstartx, jendx)
    type(poisson2ddir),intent(out) :: this
    real(wp), dimension(:,:), intent(in) :: rho
    type(geometry),intent(in)  :: geomx
    integer, intent(out) :: iflag
    integer, intent(in), optional ::  jstartx, jendx

    integer :: ierr ! indicateur d'erreur
    integer :: nxh1 ! taille du domaine en x
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

    ! definition des bandes de calcul (en n'oubliant pas le cas sequentiel)
    ! le decoupage des tableaux exterieurs n'est pas gere par le module
    ! jstart et jend sont donnees en parametre d'entree dans le cas parallele

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
!#ifndef _MPI
!    this%jstartx = 1
!    this%jendx = geomx%ny
!#endif
    ! la taille totale de la zone en kx est nxh1
    nxh1 = geomx%nx
    this%istartk = 1  ! cas sequentiel
    this%iendk = nxh1 ! cas sequentiel

#ifdef _OPENMP
    num_threads = omp_get_num_threads()
    my_num = omp_get_thread_num()
#endif
#if defined _OPENMP || defined _MPI
    ! ipiece_size = n/num_threads rounded up
    ipiece_size = (nxh1 + num_threads - 1) / num_threads
    ! zone a traiter en fonction du numero de process
    this%istartk = my_num * ipiece_size + 1
    this%iendk = min(this%istartk - 1 + ipiece_size, nxh1)
    if (this%istartk.gt.this%iendk) stop "error in new_poisson: zero size zone"
!print*,'zones ',my_num,this%istartk,this%iendk, ipiece_size,nxh1,num_threads
#endif
    ! initialisation de la geometrie
    this%geomx=geomx
    ! allocation memoire
    allocate(this%wsavex(3*(this%geomx%nx-2)+15),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, WSAVEX, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,code,ierr)
#else
       stop
#endif 
    end if
    allocate(this%wsavey(3*(this%geomx%ny-2)+15),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, WSAVEY, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,code,ierr)
#else
       stop
#endif 
    end if
    allocate(this%work(this%geomx%ny,this%geomx%nx),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, WORK, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,code,ierr)
#else
       stop
#endif 
    end if
    allocate(this%rho(this%geomx%nx,this%jstartx:this%jendx),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, RHO, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,code,ierr)
#else
       stop
#endif 
    end if
    allocate(this%rhot1(this%jstartx:this%jendx,this%geomx%nx),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, RHOT1, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,code,ierr)
#else
       stop
#endif 
    end if
    allocate(this%rhot2(this%istartk:this%iendk,this%geomx%ny),stat=ierr)
    if (ierr.ne.0) then
       write(*,*) 'error in memory allocation, RHOT2, new_poisson'
#ifdef _MPI
       call mpi_abort(MPI_COMM_WORLD,code,ierr)
#else
       stop
#endif 
    end if
    ! initialisation of sine transforms with n-2 (B points not included)
    call vsinti(this%geomx%nx-2,this%wsavex)
    call vsinti(this%geomx%ny-2,this%wsavey)

  end subroutine new_poisson2ddir
  subroutine dealloc_poisson2ddir(this)
    type(poisson2ddir),intent(out) :: this
    deallocate(this%rhot1,this%rhot2,this%wsavex,this%wsavey,this%work)
  end subroutine dealloc_poisson2ddir
  subroutine solve_poisson2ddir(this,ex,ey,rho)
    type(poisson2ddir) :: this
#ifdef _MPI
    real(wp), dimension(:,this%jstartx:) :: ex,ey,rho
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: tag, mpierror
    real(wp),dimension(this%geomx%nx) :: phib
#else
    real(wp), dimension(:,:) :: ex,ey,rho
#endif
    integer :: ik,jk,i,j,j1,j2 ! indices de boucle
    real(wp) :: kx0,ky0, kx,ky, kx2, k2
#ifdef _OPENMP
    integer :: my_num, num_threads, omp_get_num_threads, omp_get_thread_num

    num_threads = omp_get_num_threads()
    my_num = omp_get_thread_num()
#endif

    ! sauvegarde de rho pour les diagnostiques
    this%rho(:,this%jstartx:this%jendx) = rho(:,this%jstartx:this%jendx) 

    ! faire une FFST dans la direction x de rho
    ! on transpose d'abord localement rho pour vsint
    do i=2,this%geomx%nx-1
       do j=this%jstartx,this%jendx
          this%rhot1(j,i-1) = rho(i,j)
       end do
    end do

    call vsint (this%jendx-this%jstartx+1, this%geomx%nx-2, this%rhot1, &
         this%work, this%jendx-this%jstartx+1, this%wsavex) 

    ! on recopie localement le resultat dans rho pour transferer les donnees
    do i=2,this%geomx%nx-1
       do j=this%jstartx,this%jendx
          rho(i,j) = this%rhot1(j,i-1)
       end do
    end do
 
#ifdef _MPI
    ! transposition de rho dans work
    call transpose(rho, this%work(:,this%istartk:this%iendk), &
         this%geomx%nx, this%geomx%ny, num_threads)
    ! transposition locale de work dans rhot2
    do i=this%istartk,this%iendk
       do j=1,this%geomx%ny
          this%rhot2(i,j)=this%work(j,i)
       end do
    end do
#else
   
    ! transposition de rho stockee dans rhot2
    call transposexy(this,rho)
#endif
    ! faire une FSFT dans la direction y de rhot2
    call vsint (this%iendk-this%istartk+1, this%geomx%ny-2, this%rhot2(:,2:), &
         this%work, this%iendk-this%istartk+1, this%wsavey)

    ! calcul de la transformee de Fourier de phi a partir de celle de rho
    kx0=pi/((this%geomx%nx-1)*this%geomx%dx)
    ky0=pi/((this%geomx%ny-1)*this%geomx%dy)

    do ik=max(this%istartk,2),this%iendk
       kx= (ik-1)*kx0
       kx2 = kx*kx
       do jk =  2, this%geomx%ny
          ky= (jk-1)*ky0
          k2= kx2 +ky*ky
          this%rhot2(ik,jk)= this%rhot2(ik,jk)/k2
       end do       
    end do

    ! faire une FSFT inverse  dans la direction y de phi
    call vsint (this%iendk-this%istartk+1, this%geomx%ny-2, this%rhot2(:,2:), &
         this%work, this%iendk-this%istartk+1, this%wsavey)

    ! transposition de phi 
#ifdef _MPI
    ! transposition locale de rhot2 dans work
    do i=this%istartk,this%iendk
       do j=1,this%geomx%ny
          this%work(j,i)=this%rhot2(i,j)
       end do
    end do
    ! transposition MPI de work dans rho
    call transpose(this%work(:,this%istartk:this%iendk), rho, this%geomx%ny, this%geomx%nx, num_threads)

    ! on transpose localement rho pour vsint
    do i=2,this%geomx%nx-1
       do j=this%jstartx,this%jendx
          this%rhot1(j,i-1) = rho(i,j)
       end do
    end do
#else
    call transposeyx(this,rho) 
    ! faire une FSFT inverse  dans la direction x de phi
    ! on transpose d'abord localement rho pour vsint

    do i=2,this%geomx%nx-1
       do j=this%jstartx,this%jendx
          this%rhot1(j,i-1) = rho(i,j)
       end do
    end do
#endif
    call vsint (this%jendx-this%jstartx+1, this%geomx%nx-2, this%rhot1, &
         this%work, this%jendx-this%jstartx+1, this%wsavex)

    ! on recopie localement le resultat dans rho
    do i=2,this%geomx%nx
       do j=this%jstartx,this%jendx
          rho(i,j) = this%rhot1(j,i-1)
       end do
    end do

    ! les valeurs de phi sur les bords sont nulles
!!!$OMP BARRIER  
    rho(1,:) = 0.
    rho(this%geomx%nx,:) = 0.
#if defined _OPENMP || defined _MPI
    if (my_num.eq.0) then
       rho(:,1) = 0.
    end if
    if (my_num.eq.num_threads-1) then
       rho(:,this%geomx%ny) = 0.
    end if
#else
    rho(:,1) = 0.
 rho(:,this%geomx%ny) = 0.
#endif

    ! Calcul de E a partir de phi en prenant l'oppose du gradient
#ifdef _MPI
    call mpi_barrier(MPI_COMM_WORLD,mpierror)
#elif defined _OPENMP
!$OMP BARRIER  
#endif
!print*,mynum,'maxval',maxval(rho),rho(33,32),rho(34,32)

     do i = 2,this%geomx%nx-1
        ex(i,this%jstartx:this%jendx) = (rho(i-1,this%jstartx:this%jendx)- &
             rho(i+1,this%jstartx:this%jendx))/(2*this%geomx%dx)
     end do
     ex(1,this%jstartx:this%jendx)= (rho(1,this%jstartx:this%jendx) &
          -rho(2,this%jstartx:this%jendx))/this%geomx%dx
     ex(this%geomx%nx,this%jstartx:this%jendx)=         &
          (rho(this%geomx%nx-1,this%jstartx:this%jendx) &
          -rho(this%geomx%nx,this%jstartx:this%jendx))  &
          /this%geomx%dx

#ifdef _MPI
     ! compute Ey for available data
     do j=this%jstartx+1,this%jendx-1
        do i=1,this%geomx%nx
           ey(i,j)=(rho(i,j-1)-rho(i,j+1))/(2*this%geomx%dy)
        end do
     end do
     ! transfer boundary data and compute remaining values of Ey
     if (my_num.eq.0) then    ! first proc
       ey(:,1)= (rho(:,1) - rho(:,2))/this%geomx%dy
       tag=1000*my_num+my_num+1
       call mpi_send(rho(:,this%jendx),this%geomx%nx,MPI_realtype,1, &
            tag,mpi_comm_world,mpierror)
       tag=1000*(my_num+1)+my_num
       call mpi_recv(phib,this%geomx%nx,MPI_realtype,1, &
            tag,mpi_comm_world,status,mpierror)
       ey(:,this%jendx)=(rho(:,this%jendx-1)-phib)/(2*this%geomx%dy)
!print*,'phib',my_num,phib,ey(:,this%jendx)
    else if (my_num.eq.num_threads-1)  then  ! last proc
       ey(:,this%geomx%ny)= (rho(:,this%geomx%ny-1) - &
            rho(:,this%geomx%ny))/this%geomx%dy
       tag=1000*(my_num-1)+my_num
       call mpi_recv(phib,this%geomx%nx,MPI_realtype,my_num-1, &
            tag,mpi_comm_world,status,mpierror)
       tag=1000*my_num+my_num-1
       call mpi_send(rho(:,this%jstartx),this%geomx%nx,MPI_realtype,my_num-1, &
            tag,mpi_comm_world,mpierror)
       ey(:,this%jstartx)=(phib-rho(:,this%jstartx+1))/(2*this%geomx%dy)
!print*,'phib',my_num,phib
!print*,'phisent',rho(:,this%jstartx)
    else
       tag=1000*(my_num-1)+my_num
       call mpi_recv(phib,this%geomx%nx,MPI_realtype,my_num-1, &
            tag,mpi_comm_world,status,mpierror)
       ey(:,this%jstartx)=(phib-rho(:,this%jstartx+1))/(2*this%geomx%dy)
       tag=1000*my_num+my_num+1
       call mpi_send(rho(:,this%jendx),this%geomx%nx,MPI_realtype,my_num+1, &
            tag,mpi_comm_world,mpierror)
       tag=1000*my_num+my_num-1
       call mpi_send(rho(:,this%jstartx),this%geomx%nx,MPI_realtype,my_num-1, &
            tag,mpi_comm_world,mpierror)
       tag=1000*(my_num+1)+my_num
       call mpi_recv(phib,this%geomx%nx,MPI_realtype,my_num+1, &
            tag,mpi_comm_world,status,mpierror)
       ey(:,this%jendx)=(rho(:,this%jendx-1)-phib)/(2*this%geomx%dy)
    end if
#else
     ! copie rho dans rhot2 pour avoir localement tous les j pour le calcul de Ey
     do i=this%istartk,this%iendk
        do j=1,this%geomx%ny
           this%rhot2(i,j)=rho(i,j)
        end do
     end do
!     print*,mynum,'maxvalt',maxval(this%rhot2),rho(33,32),rho(34,32)
     
     do j=2,this%geomx%ny-1
        do i=this%istartk,this%iendk
           ey(i,j)=(this%rhot2(i,j-1)-this%rhot2(i,j+1))/(2*this%geomx%dy)
        end do
     end do
        ey(this%istartk:this%iendk,1)= (this%rhot2(this%istartk:this%iendk,1) &
             -this%rhot2(this%istartk:this%iendk,2))/this%geomx%dy
        ey(this%istartk:this%iendk,this%geomx%ny)= &
             (this%rhot2(this%istartk:this%iendk,this%geomx%ny-1) &
             -this%rhot2(this%istartk:this%iendk,this%geomx%ny)) &
             /this%geomx%dy
#endif

#ifdef _MPI
    call mpi_barrier(MPI_COMM_WORLD,mpierror)
#elif defined _OPENMP
!$OMP BARRIER  
#endif
     ! on recopie les valeurs de rho dans rho pour les diagnostiques
     do i=1,this%geomx%nx
        do j=this%jstartx,this%jendx
 !          rho(i,j) = this%rho(i,j)
        end do
     end do
#ifdef _MPI
    call mpi_barrier(MPI_COMM_WORLD,mpierror)
#elif defined _OPENMP
!$OMP BARRIER  
#endif
  end subroutine solve_poisson2ddir
#ifndef _MPI
  subroutine transposexy(this,rho)
    type(poisson2ddir) :: this
    real(wp), dimension(:,:) :: rho
    integer :: i,j ! indices de boucle

    ! on attend que tous les processeurs aient fait leur calcul
#ifdef _OPENMP
!$OMP BARRIER  
#endif
    ! transposition
    do i=this%istartk,this%iendk
       do j=1,this%geomx%ny
	  this%rhot2(i,j)=rho(i,j)
       end do
    end do
    ! on attend que tous les processeurs aient fait leur transposition
#ifdef _OPENMP
!$OMP BARRIER 
#endif
  end subroutine transposexy
  subroutine transposeyx(this,rho)
    type(poisson2ddir) :: this
    real(wp), dimension(:,:) :: rho
    integer :: i,j ! indices de boucle

    ! on attend que tous les processeurs aient fait leur calcul
#ifdef _OPENMP
!$OMP BARRIER
#endif
    ! transposition
    do i=this%istartk,this%iendk
       do j=1,this%geomx%ny
          rho(i,j)=this%rhot2(i,j)
       end do
    end do
    ! on attend que tous les processeurs aient fait leur transposition
#ifdef _OPENMP
!$OMP BARRIER
#endif
 end subroutine transposeyx
#endif
end module poisson2ddir_module
