module Vlasov1d_module

  use used_precision
  use splinepx_class
  use splinepy_class
  use geometry_module
  use diagnostiques_module

  implicit none
  
contains

  subroutine advection1d_x(this,f,dt,Jx1)
    !-----------------------------------------------
    ! fait une advection en x sur un pas de temps dt
    !-----------------------------------------------
    type(splinepx),intent(inout) :: this
    real(wp), dimension(:,:,:,this%jstartv:) :: f
    integer :: mpierror
    real(wp), dimension(:,:), intent(out)  :: Jx1
    real(wp), intent(in) :: dt
    ! variables locales
    real(wp), dimension(1:this%geomx%nx) :: f1dx,flux
    real(wp), dimension(this%geomx%nx,this%geomx%ny) :: locJx1
    real(wp) :: depx      ! deplacement par rapport au maillage
    real(wp) :: vx       ! vitesse du point courant
    integer :: i, jx, iv, jv, c ! indices de boucle

    ! verifier que la transposition est a jours
    if (this%transpose) stop 'advection_x: on travaille sur f et pas ft'

    locJx1=0._wp;Jx1(:,:)=0._wp
    do jv=this%jstartv,this%jendv
       do iv=1,this%geomv%nx
          vx = this%geomv%x0+(iv-1)*this%geomv%dx
          depx = vx*dt
          do jx = 1,this%geomx%ny
             f1dx(:)=f(:,jx,iv,jv)
             call interpole(this,f1dx,depx,flux)
             f(:,jx,iv,jv)=f1dx(:)
             locJx1(:,jx)=locJx1(:,jx)+flux(:)*this%geomv%dx*this%geomv%dy!*0.5_wp*vx
          enddo    
       end do
    end do
!parallel
    call mpi_barrier(MPI_COMM_WORLD,i)
    c=this%geomx%nx*this%geomx%ny
!    call mpi_allreduce(locJx1,Jx1,c, &
!         MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierror)
    call mpi_allreduce(locJx1,Jx1,c, &	
         MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierror)
!sequentiel
!    Jx1(:,:)=locJx1(:,:)

    !print*,'max val 1dx ',maxval(f),minval(f)
  end subroutine advection1d_x

 subroutine advection1d_y(this,f,dt,Jy1)
    !-----------------------------------------------
    ! fait une advection en y sur un pas de temps dt
    !-----------------------------------------------
    type(splinepy),intent(inout) :: this
    real(wp), dimension(:,:,:,this%jstartv:) :: f
    integer :: mpierror
    real(wp), dimension(:,:), intent(out)  :: Jy1
    real(wp), intent(in) :: dt
    ! variables locales
    real(wp), dimension(1:this%geomx%ny) :: f1dy,flux
    real(wp), dimension(this%geomx%nx,this%geomx%ny) :: locJy1
    real(wp) :: depy      ! deplacement par rapport au maillage
    real(wp) :: vy       ! vitesse du point courant
    integer :: i, ix, iv, jv, c ! indices de boucle

    ! verifier que la transposition est a jours
    if (this%transpose) stop 'advection_x: on travaille sur f et pas ft'

    Jy1(:,:)=0._wp;locJy1=0._wp
    do jv=this%jstartv,this%jendv
       vy = this%geomv%y0+(jv-1)*this%geomv%dy
       depy = vy*dt
       do iv=1,this%geomv%nx
          do ix=1,this%geomx%nx
             f1dy(:)=f(ix,:,iv,jv)
             call interpole(this,f1dy,depy,flux)
             f(ix,:,iv,jv)=f1dy(:)
             locJy1(ix,:)=locJy1(ix,:)+flux(:)*this%geomv%dx*this%geomv%dy!*vy*0.5_wp!
          enddo
       end do
    end do

!parallel    
    call mpi_barrier(MPI_COMM_WORLD,i)
    c=this%geomx%nx*this%geomx%ny
!    call mpi_allreduce(locJy1,Jy1,c, &
!         MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierror)
    call mpi_allreduce(locJy1,Jy1,c, &
         MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierror)
!seq
!    Jy1(:,:)=locJy1(:,:)
    
    !print*,'max val 1dy ',maxval(f),minval(f)
  end subroutine advection1d_y

end module Vlasov1d_module
