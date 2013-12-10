module Vlasov1d_module

#include "sll_working_precision.h"
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
    sll_real64, dimension(:,:,:,this%jstartv:) :: f
    sll_int32 :: mpierror
    sll_real64, dimension(:,:), intent(out)  :: Jx1
    sll_real64, intent(in) :: dt
    ! variables locales
    sll_real64, dimension(1:this%geomx%nx) :: f1dx,flux
    sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locJx1
    sll_real64 :: depx      ! deplacement par rapport au maillage
    sll_real64 :: vx       ! vitesse du point courant
    sll_int32 :: i, jx, iv, jv, c ! indices de boucle

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
    sll_real64, dimension(:,:,:,this%jstartv:) :: f
    sll_int32 :: mpierror
    sll_real64, dimension(:,:), intent(out)  :: Jy1
    sll_real64, intent(in) :: dt
    ! variables locales
    sll_real64, dimension(1:this%geomx%ny) :: f1dy,flux
    sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locJy1
    sll_real64 :: depy      ! deplacement par rapport au maillage
    sll_real64 :: vy       ! vitesse du point courant
    sll_int32 :: i, ix, iv, jv, c ! indices de boucle

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
