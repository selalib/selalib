module vlasov2d_csl2d_module

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

 use used_precision
 use splinenn_class
 use splinepp_class
 use csl2dpp_class
 use geometry_module
 use diagnostiques_module
 !use clock

 implicit none
 private
 public :: new, dealloc,advection_x, advection_v,&
           densite_charge,transposexv,transposevx,thdiag,densite_courant

 type, public :: vlasov2d
   sll_real64, dimension(:,:,:,:), pointer :: ft
   !uncomment type (csl2dpp) for new csl2d method
   ! or uncomment type(splinepp) for classical spline method
   
   !type (splinepp) :: interpx ! spline periodique pour X
   type (csl2dpp) ::  interpx ! new csl2d method
   type (csl2dpp) ::  interpv ! new csl2d method
   !type (splinenn) :: interpv ! spline naturel pour V
   type (geometry) :: geomx, geomv
   logical :: transposed       ! permet de definir si f ou ft est derniere 
                               ! fonction de distribution mise a jour
   sll_int32 :: jstartx, jendx
   sll_int32 :: jstartv, jendv  ! definition de la bande de calcul
 end type vlasov2d

 ! variables globales 
 sll_real64, dimension(:,:),allocatable :: P_x, P_y

 interface new
   module procedure new_vlasov2d
 end interface

 interface dealloc
   module procedure dealloc_vlasov2d
 end interface

contains

 subroutine new_vlasov2d(this,geomx,geomv,error, jstartx, jendx, &
                        jstartv, jendv)

  type(vlasov2d),intent(out)      :: this
  type(geometry),intent(in)       :: geomx
  type(geometry),intent(in)       :: geomv
  sll_int32, intent(out)          :: error
  sll_int32, intent(in), optional ::  jstartx
  sll_int32, intent(in), optional ::  jendx
  sll_int32, intent(in), optional ::  jstartv
  sll_int32, intent(in), optional ::  jendv

  error = 0

  ! on commence par utiliser la fonction f(x,y,vx,vy)
  this%transposed=.false.
  ! definition des bandes de calcul (en n'oubliant pas le cas sequentiel)
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
  
  ! initialisation de la geometrie
  this%geomx=geomx
  this%geomv=geomv
  ! initialisation des splines de l'espace des vitesses
  call new(this%interpv,geomv,error)
  ! initialisation de l'interpolation de l'espace physique
  call new(this%interpx,geomx,error)  

  ! allocation memoire
  SLL_ALLOCATE(this%ft(geomv%nx,geomv%ny,geomx%nx,this%jstartx:this%jendx),error)

  SLL_ALLOCATE(P_x(this%geomv%nx,this%geomv%ny),error)
  SLL_ALLOCATE(P_y(this%geomv%nx,this%geomv%ny),error)

 end subroutine new_vlasov2d

 subroutine dealloc_vlasov2d(this)
  type(vlasov2d),intent(out) :: this
  sll_int32 :: error
  SLL_DEALLOCATE(this%ft,error)
 end subroutine dealloc_vlasov2d

 !>
 !> fait une advection en x sur un pas de temps dt
 !>
 subroutine advection_x(this,f,dt)
  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:,:,this%jstartv:) :: f
  sll_real64, intent(in) :: dt

  sll_real64 :: depx, depy   ! deplacement par rapport au maillage
  sll_real64 :: vx, vy       ! vitesse du point courant
  sll_int32 :: iv, jv ! indices de boucle
!  sll_int32::timecase(0:1),interp_case,ppm_order


!  timecase(0)=3 !computation of the characteristics: 1 for Euler, 2 or 3 for symplectic Verlet
!  timecase(1)=2 !number of steps fixed point algo (symplectic Verlet case)
!  interp_case=1 !1:Lauritzen 2:LAG3 3:PPM CD 4:PPM up
!  ppm_order=2 !if interp_case=3: PPM0, PPM1 or PPM2 / if interp_case=4: up3 (1) or up5 (2)
  ! verifier que la transposition est a jour
  if (this%transposed) stop 'advection_x: on travaille sur f et pas ft'
  do jv=this%jstartv,this%jendv
     vy = this%geomv%y0+(jv-1)*this%geomv%dy
     depy = vy*dt
     do iv=1,this%geomv%nx
        vx = this%geomv%x0+(iv-1)*this%geomv%dx
        depx = vx*dt
        call interpole(this%interpx,f(:,:,iv,jv),depx,depy,jv==0)
     end do
  end do

 end subroutine advection_x

 !>
 !> fait une advection en v sur un pas de temps dt
 !>
 subroutine advection_v(this,fx,fy,dt,bz)

  type(vlasov2d),intent(inout) :: this

  sll_real64, dimension(:,:), intent(in) :: fx, fy
  sll_real64, dimension(:,:), optional, intent(in) :: bz
  sll_real64, intent(in) :: dt

  sll_real64 :: depvx, depvy   ! deplacement par rapport au maillage
  sll_int32 :: i, j ,iv, jv, im1! indices de boucle
  sll_real64 :: ctheta, stheta, px, py

  ! verifier que la transposition est a jour
  if (.not.(this%transposed)) stop 'advection_v: on travaille sur ft et pas f'

  if (present(bz)) then

   do j=this%jstartx,this%jendx
    do i=1,this%geomx%nx
     ctheta = cos(bz(i,j)*dt)
     stheta = sin(bz(i,j)*dt)
     depvx  = -0.5*dt*fx(i,j)
     depvy  = -0.5*dt*fy(i,j)
     do jv=1,this%geomv%ny-1
      py = this%geomv%y0+(jv-1)*this%geomv%dy
      do iv=1,this%geomv%nx-1 
       px = this%geomv%x0+(iv-1)*this%geomv%dx
       P_x(iv,jv) = depvx+(px+depvx)*ctheta-(py+depvy)*stheta
       P_y(iv,jv) = depvy+(px+depvx)*stheta+(py+depvy)*ctheta
      end do
     end do

     call interpole(this%interpv,this%ft(:,:,i,j),this%ft(:,:,i,j),P_x,P_y)

    end do
   end do

  else

   do j=this%jstartx,this%jendx
    do i=1,this%geomx%nx
     im1=mod(i-1+this%geomx%nx,this%geomx%nx)
     depvx = fx(i,j)*dt
     depvy = fy(i,j)*dt
!print*,i,j,depvx,depvy
     !depvx =  fx(i,j)*dt;depvy=0._wp
     call interpole(this%interpv,this%ft(:,:,i,j),depvx,depvy,j==0)
     !call interpole(this%interpv,this%ft(:,:,i,j),depvx,depvy,(j .eq. 3) .and. (i .eq. 3))
    end do
   end do
!stop
  end if
       
 end subroutine advection_v

!>------------------------------------------------
!> calcule la densite de charge rho a partir de ft
!> en fait le moment d'ordre 0 de f. Les constantes
!> ne sont pas introduites
!>------------------------------------------------
!> Poisson n'est pas parallele on transmet donc rho
!> a tous les processeurs
subroutine densite_charge(this, rho)

   type(vlasov2d),intent(inout) :: this
   sll_int32 :: error
   sll_real64, dimension(:,:), intent(out)  :: rho
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locrho
   sll_int32 :: i,j,iv,jv,c
   sll_int32 :: comm

   comm   = sll_world_collective%comm

   SLL_ASSERT(this%transposed)
   
   rho(:,:) = 0.
   locrho(:,:) = 0.
   do j=this%jstartx,this%jendx
      do i=1,this%geomx%nx
         do jv=1,this%geomv%ny-1
            do iv=1,this%geomv%nx-1 
               locrho(i,j) = locrho(i,j) + this%geomv%dx*this%geomv%dy* &
                    this%ft(iv,jv,i,j) 
            end do
         end do
      end do
   end do
   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locrho,rho,c,MPI_REAL8,MPI_SUM,comm,error)

end subroutine densite_charge

!>------------------------------------------------
!> calcule la densite de courant jx et jy a partir de ft
!> en fait le moment d'ordre 0 de f. Les constantes
!> ne sont pas introduites
!>------------------------------------------------
subroutine densite_courant(this, jx, jy)

   type(vlasov2d),intent(inout) :: this
   sll_int32 :: error
   sll_real64, dimension(:,:), intent(out)  :: jx, jy
   sll_real64 :: vx, vy       ! vitesse du point courant
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjx
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjy
   sll_int32 :: i,j,iv,jv,c
   sll_int32 :: comm

   comm   = sll_world_collective%comm

   ! verifier que la transposition est a jour
   if (.not.(this%transposed)) &
      stop 'densite_courant: on travaille sur ft et pas f'
   !    rho(:,this%jstartx:this%jendx)=0.
   jx(:,:) = 0.; jy(:,:) = 0.
   locjx(:,:) = 0.; locjy(:,:) = 0.
   do j=this%jstartx,this%jendx
      do i=1,this%geomx%nx
         do jv=1,this%geomv%ny-1
            vy = this%geomv%y0+(jv-1)*this%geomv%dy
            do iv=1,this%geomv%nx-1 
               vx = this%geomv%x0+(iv-1)*this%geomv%dx
               locjx(i,j) = locjx(i,j) + this%geomv%dx*this%geomv%dy* &
                    this%ft(iv,jv,i,j) * vx
               locjy(i,j) = locjy(i,j) + this%geomv%dx*this%geomv%dy* &
                    this%ft(iv,jv,i,j) * vy
            end do
         end do
      end do
   end do

   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locjx,jx,c, MPI_REAL8,MPI_SUM,comm,error)
   call mpi_allreduce(locjy,jy,c, MPI_REAL8,MPI_SUM,comm,error)

end subroutine densite_courant

!>---------------------------------------------
!> transpose la fonction de distribution
!> ATTENTION: cet fonction fait intervenir des 
!> communications entre les processeurs.
!>---------------------------------------------
subroutine transposexv(this,f)
   type(vlasov2d),intent(inout) :: this
   sll_real64, dimension(:,:,:,:),intent(in) :: f
   sll_int32 :: sizexy, sizevxvy
   sll_int32 :: my_num, num_threads

   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   sizexy = this%geomx%nx * this%geomx%ny
   sizevxvy = this%geomv%nx * this%geomv%ny
   call transpose(f, this%ft, sizexy, sizevxvy, num_threads)

this%transposed=.true.

end subroutine transposexv

subroutine transposevx(this,f)

   type(vlasov2d),intent(inout) :: this
   sll_real64, dimension(:,:,:,:),intent(out) :: f

   sll_int32 :: sizexy, sizevxvy
   sll_int32 :: my_num, num_threads

   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   sizexy = this%geomx%nx * this%geomx%ny
   sizevxvy = this%geomv%nx * this%geomv%ny
   call transpose(this%ft,f, sizevxvy, sizexy, num_threads)

   this%transposed=.false.

end subroutine transposevx

subroutine thdiag(this,f,nrj,t,jstartv)

   type(vlasov2d),intent(inout) :: this
   sll_int32, intent(in) :: jstartv
   sll_real64, dimension(:,:,:,jstartv:),intent(in) :: f
   !sll_int32 :: error
   sll_real64, intent(in) :: t,nrj   ! current time
   ! variables locales
   !sll_int32 :: i,iv, j,jv
   !sll_real64 :: x, vx, y, vy
   !sll_real64,dimension(7) :: diagloc
   !sll_real64,dimension(11) :: auxloc
   sll_real64,dimension(13) :: aux
   sll_real64,dimension(0:9) :: diag
   sll_int32 :: my_num, num_threads
   sll_int32 :: comm

   comm   = sll_world_collective%comm
   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   if (my_num == MPI_MASTER) then
      diag=0.
      aux=0.
   end if

!   diagloc = 0._wp
!   auxloc  = 0._wp
!   do i = 1,this%geomx%nx
!      x = this%geomx%x0+(i-1)*this%geomx%dx
!      do j = 1,this%geomx%ny
!         y= this%geomx%y0+(j-1)*this%geomx%dy
!         do iv=1,this%geomv%nx
!            vx = this%geomv%x0+(iv-1)*this%geomv%dx
!            do jv=this%jstartv,this%jendv
!               vy = this%geomv%y0+(jv-1)*this%geomv%dy
!               diagloc(2) = diagloc(2) + f(i,j,iv,jv)*f(i,j,iv,jv)
!
!               auxloc(1) = auxloc(1) + f(i,j,iv,jv)         ! avg(f)
!               auxloc(2) = auxloc(2) + x*f(i,j,iv,jv)       ! avg(x)
!               auxloc(3) = auxloc(3) + vx*f(i,j,iv,jv)      ! avg(vx)
!               auxloc(4) = auxloc(4) + x*x*f(i,j,iv,jv)     ! avg(x^2)
!               auxloc(5) = auxloc(5) + vx*vx*f(i,j,iv,jv)   ! avg(vx^2)
!               auxloc(6) = auxloc(6) + x*vx*f(i,j,iv,jv)    ! avg(x*vx)
!               auxloc(7) = auxloc(7) + y*f(i,j,iv,jv)       ! avg(y)
!               auxloc(8) = auxloc(8) + vy*f(i,j,iv,jv)      ! avg(vy)
!               auxloc(9) = auxloc(9) + y*y*f(i,j,iv,jv)     ! avg(y^2)
!               auxloc(10) = auxloc(10) + vy*vy*f(i,j,iv,jv) ! avg(vy^2)
!               auxloc(11) = auxloc(11) + y*vy*f(i,j,iv,jv)  ! avg(y*vy)
!
!            end do
!         end do
!      end do
!   end do
!   auxloc=auxloc!*this%geomx%dx*this%geomx%dy*this%geomv%dx*this%geomv%dy
!
!   call mpi_reduce(auxloc,aux,11,MPI_REAL8,MPI_SUM,0,  &
!                   comm, error)
!   call mpi_reduce(diagloc(2),diag(2),1,MPI_REAL8,MPI_SUM,0, &
!                   comm, error)

if (my_num==MPI_MASTER) then
   aux(13)=t
   aux(12)=nrj
   write(*,"('time ', g8.3,' test nrj',f10.5)") t, nrj
   call time_history("thf","(13(1x,e15.6))",aux(1:13))
end if

end subroutine thdiag

end module vlasov2d_csl2d_module


